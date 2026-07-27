// Contracted one-electron S/T/V matrix construction from the RECURSUM-emitted
// Obara-Saika kernels, using the SAME primitive-loop contraction approach as the
// ERI/JK drivers (loop primitive pairs, weight by contraction coeff x per-shell
// norm, accumulate; canonical la>=lb kernels with swap+transpose for la<lb).
// Validates against a PySCF int1e_ovlp/kin/nuc_cart reference.
//   g++ -O3 -march=native -std=c++17 onee_driver.cpp onee_kernels.cpp -o onee_driver
//   ./onee_driver mol.txt
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <array>
#include <algorithm>
#include "recursum_onee_scalars.h"
#include "onee_dispatch.h"

using std::vector;
static int ncart(int l){ return (l+1)*(l+2)/2; }

struct Shell { int l; double c[3]; vector<double> e, d, N; int off, ncart; };

// libcint per-shell Cartesian normalization: norm_l (L<=1) else gto_norm.
static double dfact(int n){ long r=1; for(int i=1;i<=n;i++) r*=2*i-1; return n>0?(double)r:1.0; }
static double norm_l(int l,double a){ return std::pow(2*a/M_PI,0.75)*std::pow(4*a,l/2.0)/std::sqrt(dfact(l)); }
static double gto_norm(int l,double a){
    double num=std::pow(2.0,2*l+3); for(int k=2;k<=l+1;k++) num*=k;      // (l+1)!
    double den=1; for(int k=2;k<=2*l+2;k++) den*=k;                       // (2l+2)!
    return std::sqrt(num*std::pow(2*a,l+1.5)/(den*std::sqrt(M_PI)));
}
static double cart_norm(int l,double a){ return l<=1?norm_l(l,a):gto_norm(l,a); }

// Boys F_0..F_mmax(T): series for the top order + downward recursion (large-T
// erf/upward branch). Same robust routine as the J/K driver.
static void boys(int mmax,double T,double*F){
    if(T<1e-13){ for(int m=0;m<=mmax;m++)F[m]=1.0/(2*m+1); return; }
    double eT=std::exp(-T);
    if(T>30.0){ F[0]=0.5*std::sqrt(M_PI/T)*std::erf(std::sqrt(T));
        for(int m=0;m<mmax;m++) F[m+1]=((2*m+1)*F[m]-eT)/(2*T); return; }
    double sum=0,term=1.0/(2*mmax+1); int k=0;
    do{ sum+=term; k++; term*=T/(mmax+0.5+k); }while(term>1e-18*sum && k<500);
    F[mmax]=eT*sum;
    for(int m=mmax;m>0;m--) F[m-1]=(2*T*F[m]+eT)/(2*m-1);
}

static vector<Shell> read_mol(const std::string&path,int&nao,
        vector<double>&Z,vector<std::array<double,3>>&Rc,
        vector<double>&S,vector<double>&Tk,vector<double>&Vn){
    std::ifstream f(path); std::string line; std::getline(f,line);
    int nsh=std::stoi(line); vector<Shell> sh(nsh); int off=0;
    for(int i=0;i<nsh;i++){ std::getline(f,line); std::istringstream s(line);
        Shell&S=sh[i]; s>>S.l>>S.c[0]>>S.c[1]>>S.c[2]; int np; s>>np;
        S.e.resize(np);S.d.resize(np);S.N.resize(np);
        for(int p=0;p<np;p++)s>>S.e[p];
        for(int p=0;p<np;p++){s>>S.d[p];S.N[p]=cart_norm(S.l,S.e[p]);}
        S.ncart=ncart(S.l);S.off=off;off+=S.ncart; }
    std::getline(f,line); nao=std::stoi(line);
    std::getline(f,line); int natm=std::stoi(line);
    Z.resize(natm); Rc.resize(natm);
    for(int k=0;k<natm;k++){ std::getline(f,line); std::istringstream s(line);
        s>>Z[k]>>Rc[k][0]>>Rc[k][1]>>Rc[k][2]; }
    auto readM=[&](vector<double>&M){ M.resize((size_t)nao*nao); std::getline(f,line);
        std::istringstream s(line); for(auto&x:M)s>>x; };
    readM(S); readM(Tk); readM(Vn);
    return sh;
}

// Contracted block (ncA x ncB) for op in {S,T,V}; kernels need la>=lb so for
// lA<lB we compute the (B,A) block and transpose. `blk[a*ncB+b]`.
enum Op { OVLP, KIN, NUC };
static void block(const Shell&SA,const Shell&SB,Op op,
                  const vector<double>&Z,const vector<std::array<double,3>>&Rc,
                  double*blk){
    bool swap = SA.l < SB.l;
    const Shell&F = swap?SB:SA;   // "a" shell (la>=lb)
    const Shell&G = swap?SA:SB;   // "b" shell
    int la=F.l, lb=G.l, nF=F.ncart, nG=G.ncart;
    int nout=nF*nG, maxm=la+lb;
    vector<double> raw(nout,0.0), tmp(nout);
    ovlp_fn ov=onee_ovlp_dispatch(la,lb);
    nuc_fn  nu=onee_nuc_dispatch(la,lb);
    kin_fn  ki=onee_kin_dispatch(la,lb);
    for(size_t ip=0;ip<F.e.size();++ip){ double al=F.e[ip], wF=F.d[ip]*F.N[ip];
      for(size_t jp=0;jp<G.e.size();++jp){ double be=G.e[jp], wG=G.d[jp]*G.N[jp];
        double zeta=al+be, w=wF*wG;
        double P[3],AB[3];
        for(int t=0;t<3;t++){ P[t]=(al*F.c[t]+be*G.c[t])/zeta; AB[t]=F.c[t]-G.c[t]; }
        double r2ab=AB[0]*AB[0]+AB[1]*AB[1]+AB[2]*AB[2];
        double K=std::exp(-al*be/zeta*r2ab);
        OneEScalars s{};
        s.PAx=P[0]-F.c[0]; s.PAy=P[1]-F.c[1]; s.PAz=P[2]-F.c[2];
        s.ABx=AB[0]; s.ABy=AB[1]; s.ABz=AB[2];
        s.inv_2zeta=0.5/zeta; s.beta=be;
        if(op==OVLP){ double S00=K*std::pow(M_PI/zeta,1.5);
            ov(s,S00,tmp.data()); for(int i=0;i<nout;i++) raw[i]+=w*tmp[i]; }
        else if(op==KIN){ double S00=K*std::pow(M_PI/zeta,1.5);
            ki(s,S00,tmp.data()); for(int i=0;i<nout;i++) raw[i]+=w*tmp[i]; }
        else { // nuclear: sum over nuclei
            for(size_t k=0;k<Z.size();++k){
                double PC[3]; for(int t=0;t<3;t++) PC[t]=P[t]-Rc[k][t];
                double Tn=zeta*(PC[0]*PC[0]+PC[1]*PC[1]+PC[2]*PC[2]);
                double F0[32]; boys(maxm,Tn,F0);
                double pref=-2.0*M_PI/zeta*K*Z[k];
                double kf[32]; for(int m=0;m<=maxm;m++) kf[m]=pref*F0[m];
                s.PCx=PC[0]; s.PCy=PC[1]; s.PCz=PC[2];
                nu(s,kf,tmp.data()); for(int i=0;i<nout;i++) raw[i]+=w*tmp[i];
            }
        }
      }
    }
    // place into blk[a(SA)*ncB(SB)+b(SB)]
    int ncB=SB.ncart;
    if(!swap) for(int a=0;a<nF;a++)for(int b=0;b<nG;b++) blk[a*ncB+b]=raw[a*nG+b];
    else      for(int a=0;a<SA.ncart;a++)for(int b=0;b<SB.ncart;b++) blk[a*ncB+b]=raw[b*nG+a];
}

int main(int argc,char**argv){
    if(argc<2){ fprintf(stderr,"usage: %s mol.txt\n",argv[0]); return 2; }
    int nao; vector<double> Z; vector<std::array<double,3>> Rc;
    vector<double> Sr,Tr,Vr;
    auto sh=read_mol(argv[1],nao,Z,Rc,Sr,Tr,Vr);
    int maxL=0; for(auto&s:sh) maxL=std::max(maxL,s.l);
    const char* names[3]={"S (overlap)","T (kinetic)","V (nuclear)"};
    Op ops[3]={OVLP,KIN,NUC}; const vector<double>* refs[3]={&Sr,&Tr,&Vr};
    printf("RECURSUM one-electron S/T/V vs PySCF int1e_*_cart  (nao=%d, maxL=%d)\n",nao,maxL);
    bool allok=true;
    for(int o=0;o<3;o++){
        vector<double> M((size_t)nao*nao,0.0), blk(256);
        for(auto&SA:sh)for(auto&SB:sh){
            if(SB.off>SA.off) continue;                 // lower triangle by offset
            block(SA,SB,ops[o],Z,Rc,blk.data());
            for(int a=0;a<SA.ncart;a++)for(int b=0;b<SB.ncart;b++){
                double v=blk[a*SB.ncart+b];
                M[(size_t)(SA.off+a)*nao+(SB.off+b)]=v;
                M[(size_t)(SB.off+b)*nao+(SA.off+a)]=v;   // symmetric
            }
        }
        const vector<double>&R=*refs[o];
        double e=0,m=0;
        for(size_t i=0;i<M.size();i++){ m=std::max(m,std::fabs(R[i])); e=std::max(e,std::fabs(M[i]-R[i])); }
        double rel=e/(m>1e-300?m:1e-300); bool ok=rel<1e-10; allok&=ok;
        printf("  %-12s max abs err=%.2e  rel err=%.2e  %s\n",names[o],e,rel,ok?"PASS":"FAIL");
    }
    printf(allok?"ALL PASS\n":"FAILURES\n");
    return allok?0:1;
}
