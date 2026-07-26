// Contracted J/K construction using the DAG Obara-Saika primitive ERI kernels,
// libint2-style contraction (loop primitive quartets, accumulate), and BLAS
// (DGEMM/DGEMV) digestion with the density matrix. Validates against a PySCF
// reference and times the J and K build on the alkane series.
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <chrono>
#include <cblas.h>

#include "recursum_eri_scalars.h"
#include "jk_dispatch.h"

// ---- libint2 path (canonical permute -> build_eri -> unpermute) -----------
#include <libint2.h>
#include <libint2/boys.h>
libint2::FmEval_Chebyshev7<double> fmeval_chebyshev(28);
libint2::FmEval_Taylor<double,6> fmeval_taylor(1e-15);
#include "prep_libint2.h"
static Libint_t* g_eval=nullptr;
static RandomShellSet<4u>* g_rs=nullptr;   // preallocated (veclen=1, contrdepth=1), mutated per quartet

using std::vector;
static double now(){using namespace std::chrono;return duration<double>(steady_clock::now().time_since_epoch()).count();}

struct Shell { int l; double c[3]; vector<double> e, d, N; int off, ncart; };

static double radial_norm(int l, double a){
    // (2a/pi)^{3/4} (4a)^{l/2} / sqrt((2l-1)!!)
    long df=1; for(int i=1;i<=l;i++) df*=2*i-1;
    return std::pow(2*a/M_PI,0.75)*std::pow(4*a,l/2.0)/std::sqrt((double)df);
}
// Boys F_m(T), m=0..mmax: series for the top order then downward recursion.
static void boys(int mmax,double T,double*F){
    if(T<1e-13){ for(int m=0;m<=mmax;m++)F[m]=1.0/(2*m+1); return; }
    double eT=std::exp(-T);
    if(T>30.0){ // large T: F0 from erf, upward recursion is stable
        F[0]=0.5*std::sqrt(M_PI/T)*std::erf(std::sqrt(T));
        for(int m=0;m<mmax;m++) F[m+1]=((2*m+1)*F[m]-eT)/(2*T);
        return;
    }
    // moderate T: top order by Kummer series, then downward recursion
    double sum=0,term=1.0/(2*mmax+1); int k=0;
    do{ sum+=term; k++; term*=T/(mmax+0.5+k);}while(term>1e-18*sum && k<500);
    F[mmax]=eT*sum;
    for(int m=mmax;m>0;m--) F[m-1]=(2*T*F[m]+eT)/(2*m-1);
}
static int ncart(int l){return (l+1)*(l+2)/2;}

static vector<Shell> read_mol(const std::string&path,int&nao,vector<double>&D,
                              vector<double>&Jr,vector<double>&Kr){
    std::ifstream f(path); std::string line; std::getline(f,line);
    int nsh=std::stoi(line); vector<Shell> sh(nsh); int off=0;
    for(int i=0;i<nsh;i++){ std::getline(f,line); std::istringstream s(line);
        Shell&S=sh[i]; s>>S.l>>S.c[0]>>S.c[1]>>S.c[2]; int np; s>>np;
        S.e.resize(np);S.d.resize(np);S.N.resize(np);
        for(int p=0;p<np;p++)s>>S.e[p];
        for(int p=0;p<np;p++){s>>S.d[p];S.N[p]=radial_norm(S.l,S.e[p]);}
        S.ncart=ncart(S.l); S.off=off; off+=S.ncart;
    }
    std::getline(f,line); nao=std::stoi(line);
    auto readM=[&](vector<double>&M){ M.resize((size_t)nao*nao); std::getline(f,line);
        std::istringstream s(line); for(auto&x:M)s>>x; };
    readM(D); readM(Jr); readM(Kr);
    return sh;
}

// per-primitive scalars + prescaled base integrals (mirror oracle.compute_scalars)
static void prim_scalars(double a,const double*A,double b,const double*B,
                         double g,const double*C,double d,const double*Dp,
                         ScalarPack&s,double*kf,int maxm){
    double zp=a+b,zq=g+d,zpq=zp+zq,rho=zp*zq/zpq,pAB=a*b/zp,pCD=g*d/zq;
    double P[3],Q[3],W[3];
    for(int i=0;i<3;i++){P[i]=(a*A[i]+b*B[i])/zp;Q[i]=(g*C[i]+d*Dp[i])/zq;W[i]=(zp*P[i]+zq*Q[i])/zpq;}
    s.PAx=P[0]-A[0];s.PAy=P[1]-A[1];s.PAz=P[2]-A[2];
    s.QCx=Q[0]-C[0];s.QCy=Q[1]-C[1];s.QCz=Q[2]-C[2];
    s.WPx=W[0]-P[0];s.WPy=W[1]-P[1];s.WPz=W[2]-P[2];
    s.WQx=W[0]-Q[0];s.WQy=W[1]-Q[1];s.WQz=W[2]-Q[2];
    s.ABx=B[0]-A[0];s.ABy=B[1]-A[1];s.ABz=B[2]-A[2];
    s.CDx=Dp[0]-C[0];s.CDy=Dp[1]-C[1];s.CDz=Dp[2]-C[2];
    s.inv_2zp=0.5/zp;s.inv_2zq=0.5/zq;s.inv_2zpq=0.5/zpq;
    s.frac_q_over_pq=zq/zpq;s.frac_p_over_pq=zp/zpq;
    double r2ab=0,r2cd=0,r2pq=0;
    for(int i=0;i<3;i++){r2ab+=(A[i]-B[i])*(A[i]-B[i]);r2cd+=(C[i]-Dp[i])*(C[i]-Dp[i]);r2pq+=(P[i]-Q[i])*(P[i]-Q[i]);}
    double K=34.986836655249725693/(zp*zq*std::sqrt(zpq))*std::exp(-pAB*r2ab)*std::exp(-pCD*r2cd);
    double F[32]; fmeval_chebyshev.eval(F, rho*r2pq, maxm);   // same fast Boys as libint2 (fair)
    for(int m=0;m<=maxm;m++)kf[m]=K*F[m];
}

// libint2 primitive ERI block for (SA,SB,SC,SD) primitives, written into out[a,b,c,d]
// (same layout as the DAG kernels) via canonical ordering + unpermute.
static void libint2_prim_block(const Shell&SA,const Shell&SB,const Shell&SC,const Shell&SD,
                               int pa,int pb,int pc,int pd,double*out){
    const Shell* sh[4]={&SA,&SB,&SC,&SD}; int prim[4]={pa,pb,pc,pd};
    int L[4]={SA.l,SB.l,SC.l,SD.l};
    int ax[4]={0,1,2,3};                    // canonical: la>=lb, lc>=ld, bra<=ket
    if(L[ax[0]]<L[ax[1]]) std::swap(ax[0],ax[1]);
    if(L[ax[2]]<L[ax[3]]) std::swap(ax[2],ax[3]);
    if(L[ax[0]]+L[ax[1]] > L[ax[2]]+L[ax[3]]){ std::swap(ax[0],ax[2]); std::swap(ax[1],ax[3]); }
    int cl[4]={L[ax[0]],L[ax[1]],L[ax[2]],L[ax[3]]};
    RandomShellSet<4u>&rs=*g_rs;            // reuse preallocated buffers (no heap alloc per call)
    for(int p=0;p<4;p++){ const Shell*S=sh[ax[p]];
        rs.l[p]=cl[p];
        for(int k=0;k<3;k++) rs.R[p][k]=S->c[k];
        rs.exp[p][0][0]=S->e[prim[ax[p]]]; rs.coef[p][0][0]=1.0; }
    prep_libint2(g_eval,rs,0,0); g_eval->contrdepth=1;
    LIBINT2_PREFIXED_NAME(libint2_build_eri)[cl[0]][cl[1]][cl[2]][cl[3]](g_eval);
    const double* cb=g_eval->targets[0];
    int nc[4],no[4]; for(int p=0;p<4;p++){nc[p]=ncart(cl[p]); no[p]=ncart(L[p]);}
    for(int c0=0;c0<no[0];c0++)for(int c1=0;c1<no[1];c1++)
    for(int c2=0;c2<no[2];c2++)for(int c3=0;c3<no[3];c3++){
        int comp[4]={c0,c1,c2,c3};
        int i0=comp[ax[0]],i1=comp[ax[1]],i2=comp[ax[2]],i3=comp[ax[3]];
        out[((c0*no[1]+c1)*no[2]+c2)*no[3]+c3]=cb[((i0*nc[1]+i1)*nc[2]+i2)*nc[3]+i3];
    }
}

int main(int argc,char**argv){
    if(argc<3){fprintf(stderr,"usage: %s <recursum|libint2> alkane.txt [reps]\n",argv[0]);return 2;}
    std::string impl=argv[1]; bool use_l2=(impl=="libint2");
    int reps = argc>3?atoi(argv[3]):5;
    if(use_l2){ LIBINT2_PREFIXED_NAME(libint2_static_init)();
        g_eval=new Libint_t; LIBINT2_PREFIXED_NAME(libint2_init_eri)(g_eval,LIBINT2_MAX_AM_eri,0);
        static uint am0[4]={1,1,1,1}; g_rs=new RandomShellSet<4u>(am0,1,1); }
    int nao; vector<double> D,Jr,Kr;
    auto sh=read_mol(argv[2],nao,D,Jr,Kr);

    vector<double> J((size_t)nao*nao), K((size_t)nao*nao);
    // block scratch
    vector<double> blk(81), tmp(81), out(81);
    double tJ=0,tK=0,tlast=0;

    for(int rep=0; rep<reps; ++rep){
        std::fill(J.begin(),J.end(),0.0);
        std::fill(K.begin(),K.end(),0.0);
        double t0=now();
        for(auto&SA:sh)for(auto&SB:sh)for(auto&SC:sh)for(auto&SD:sh){
            int la=SA.l,lb=SB.l,lc=SC.l,ld=SD.l;
            int na=SA.ncart,nb=SB.ncart,ncc=SC.ncart,nd=SD.ncart;
            int nbra=na*nb, nket=ncc*nd, nout=nbra*nket;
            rfn fn=jk_dispatch(la,lb,lc,ld);
            // contracted block (na,nb,ncc,nd), row-major a,b,c,d
            std::fill(blk.begin(),blk.begin()+nout,0.0);
            int maxm=la+lb+lc+ld;
            for(size_t pa=0;pa<SA.e.size();++pa)for(size_t pb=0;pb<SB.e.size();++pb)
            for(size_t pc=0;pc<SC.e.size();++pc)for(size_t pd=0;pd<SD.e.size();++pd){
                if(use_l2){ libint2_prim_block(SA,SB,SC,SD,pa,pb,pc,pd,out.data()); }
                else { ScalarPack s; double kf[32];
                  prim_scalars(SA.e[pa],SA.c,SB.e[pb],SB.c,SC.e[pc],SC.c,SD.e[pd],SD.c,s,kf,maxm);
                  fn(s,kf,out.data()); }
                double w=SA.d[pa]*SA.N[pa]*SB.d[pb]*SB.N[pb]*SC.d[pc]*SC.N[pc]*SD.d[pd]*SD.N[pd];
                for(int i=0;i<nout;i++) blk[i]+=w*out[i];
            }
            // Digest. blk index = ((a*nb+b)*ncc+c)*nd+d
            // J[A,B] += sum_{c,d} blk[a,b,c,d] * D[C+c, D+d]   (GEMV: (nbra x nket) * Dket)
            for(int c=0;c<ncc;c++)for(int dd=0;dd<nd;dd++) tmp[c*nd+dd]=D[(size_t)(SC.off+c)*nao+(SD.off+dd)];
            // J_AB (nbra) += blk(nbra,nket) * tmp(nket)
            static vector<double> jab(81);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,nbra,nket,1.0,blk.data(),nket,tmp.data(),1,0.0,jab.data(),1);
            for(int a=0;a<na;a++)for(int b=0;b<nb;b++) J[(size_t)(SA.off+a)*nao+(SB.off+b)]+=jab[a*nb+b];
        }
        tJ+=now()-t0;
        // K in a separate timed pass (index pattern (mu lam|nu sig))
        t0=now();
        for(auto&SA:sh)for(auto&SB:sh)for(auto&SC:sh)for(auto&SD:sh){
            int la=SA.l,lb=SB.l,lc=SC.l,ld=SD.l;
            int na=SA.ncart,nb=SB.ncart,ncc=SC.ncart,nd=SD.ncart;
            int nout=na*nb*ncc*nd, maxm=la+lb+lc+ld;
            rfn fn=jk_dispatch(la,lb,lc,ld);
            std::fill(blk.begin(),blk.begin()+nout,0.0);
            for(size_t pa=0;pa<SA.e.size();++pa)for(size_t pb=0;pb<SB.e.size();++pb)
            for(size_t pc=0;pc<SC.e.size();++pc)for(size_t pd=0;pd<SD.e.size();++pd){
                if(use_l2){ libint2_prim_block(SA,SB,SC,SD,pa,pb,pc,pd,out.data()); }
                else { ScalarPack s; double kf[32];
                  prim_scalars(SA.e[pa],SA.c,SB.e[pb],SB.c,SC.e[pc],SC.c,SD.e[pd],SD.c,s,kf,maxm);
                  fn(s,kf,out.data()); }
                double w=SA.d[pa]*SA.N[pa]*SB.d[pb]*SB.N[pb]*SC.d[pc]*SC.N[pc]*SD.d[pd]*SD.N[pd];
                for(int i=0;i<nout;i++) blk[i]+=w*out[i];
            }
            // K[A,C] += sum_{b,d} blk[a,b,c,d] * D[B+b, D+d]
            for(int a=0;a<na;a++)for(int c=0;c<ncc;c++)for(int b=0;b<nb;b++)for(int dd=0;dd<nd;dd++){
                double v=blk[((a*nb+b)*ncc+c)*nd+dd];
                K[(size_t)(SA.off+a)*nao+(SC.off+c)]+=v*D[(size_t)(SB.off+b)*nao+(SD.off+dd)];
            }
        }
        tK+=now()-t0;
    }
    // validate against reference (last J/K)
    double je=0,ke=0,jm=0,km=0;
    for(size_t i=0;i<J.size();i++){ jm=std::max(jm,std::fabs(Jr[i])); km=std::max(km,std::fabs(Kr[i]));
        je=std::max(je,std::fabs(J[i]-Jr[i])); ke=std::max(ke,std::fabs(K[i]-Kr[i])); }
    double Jms=tJ/reps*1e3, Kms=tK/reps*1e3;
    printf("%-9s %-14s nao=%d  J=%.3f ms  K=%.3f ms   Jrelerr=%.1e Krelerr=%.1e %s\n",
           impl.c_str(), argv[2], nao, Jms, Kms, je/jm, ke/km,
           (je/jm<1e-9&&ke/km<1e-9)?"OK":"MISMATCH");
    return 0;
}
