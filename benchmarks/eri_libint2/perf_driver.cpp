// perf-instrumented driver: run ONE kernel variant over a STREAM of many
// distinct primitive quartets, so per-quartet working-set / cache effects are
// real (not a single hot-in-L1 quartet). Wrap with `perf stat` to attribute
// cache misses / cycles to exactly this kernel.
//
//   ./perf_driver <impl> <class> <nquartets> <reps>
//     impl  = recursum | recursum_naive | libint2
//
// Prereqs for all quartets are precomputed into large arrays (sized to exceed
// cache) BEFORE timing, so the measured region is only the recurrence kernel
// streaming over them -- identical protocol for all three impls.
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "recursum_eri_scalars.h"
#include "recursum_eri_naive_dispatch.h"
#include "eri_classes.h"     // generated Cls table + dispatchers

#include <libint2.h>
#include <libint2/boys.h>
libint2::FmEval_Chebyshev7<double> fmeval_chebyshev(28);
libint2::FmEval_Taylor<double, 6> fmeval_taylor(1e-15);
#include "prep_libint2.h"

static const Cls* find_cls(const char* n){
    for (auto& c : ERI_CLASSES) if(!strcmp(c.name,n)) return &c; return nullptr;
}
static recursum_fn_t recursum_fn(const char* c, bool naive){
    return naive ? recursum_dispatch_naive(c) : recursum_dispatch(c);
}

// deterministic pseudo-random quartet i (vary geometry + exponents)
static void gen_quartet(int i, double* A,double* B,double* C,double* D,
                        double& a,double& b,double& g,double& d){
    auto f=[&](int k){ double x=std::sin(0.1*i+k*1.7)+0.5*std::cos(0.03*i*k+2.0);
                       return x; };
    for(int k=0;k<3;k++){A[k]=0.3*f(k);B[k]=0.5+0.3*f(k+3);C[k]=0.1+0.3*f(k+6);D[k]=-0.2+0.3*f(k+9);}
    a=0.6+0.5*(1+std::sin(0.07*i)); b=0.5+0.4*(1+std::cos(0.05*i));
    g=0.7+0.5*(1+std::sin(0.09*i+1)); d=0.4+0.4*(1+std::cos(0.06*i+2));
}

static void recursum_scalars(double a,double b,double g,double d,
        const double*A,const double*B,const double*C,const double*D,
        ScalarPack& s, double* kf, int max_m){
    double zp=a+b,zq=g+d,zpq=zp+zq,rho=zp*zq/zpq,pAB=a*b/zp,pCD=g*d/zq;
    double P[3],Q[3],W[3];
    for(int i=0;i<3;i++){P[i]=(a*A[i]+b*B[i])/zp;Q[i]=(g*C[i]+d*D[i])/zq;W[i]=(zp*P[i]+zq*Q[i])/zpq;}
    s.PAx=P[0]-A[0];s.PAy=P[1]-A[1];s.PAz=P[2]-A[2];
    s.QCx=Q[0]-C[0];s.QCy=Q[1]-C[1];s.QCz=Q[2]-C[2];
    s.WPx=W[0]-P[0];s.WPy=W[1]-P[1];s.WPz=W[2]-P[2];
    s.WQx=W[0]-Q[0];s.WQy=W[1]-Q[1];s.WQz=W[2]-Q[2];
    s.ABx=B[0]-A[0];s.ABy=B[1]-A[1];s.ABz=B[2]-A[2];
    s.CDx=D[0]-C[0];s.CDy=D[1]-C[1];s.CDz=D[2]-C[2];
    s.inv_2zp=0.5/zp;s.inv_2zq=0.5/zq;s.inv_2zpq=0.5/zpq;
    s.frac_q_over_pq=zq/zpq;s.frac_p_over_pq=zp/zpq;
    double R2AB=0,R2CD=0,R2PQ=0;
    for(int i=0;i<3;i++){R2AB+=(A[i]-B[i])*(A[i]-B[i]);R2CD+=(C[i]-D[i])*(C[i]-D[i]);R2PQ+=(P[i]-Q[i])*(P[i]-Q[i]);}
    double K=34.986836655249725693/(zp*zq*std::sqrt(zpq))*std::exp(-pAB*R2AB)*std::exp(-pCD*R2CD);
    double F[32]; libint2::FmEval_Reference<double>::eval(F,rho*R2PQ,max_m);
    for(int m=0;m<=max_m;m++) kf[m]=K*F[m];
}

int main(int argc,char**argv){
    if(argc<5){fprintf(stderr,"usage: %s <impl> <class> <nquartets> <reps>\n",argv[0]);return 2;}
    std::string impl=argv[1]; const Cls* cl=find_cls(argv[2]);
    int NQ=atoi(argv[3]), REPS=atoi(argv[4]);
    if(!cl){fprintf(stderr,"bad class\n");return 2;}
    int max_m=cl->l[0]+cl->l[1]+cl->l[2]+cl->l[3];

    LIBINT2_PREFIXED_NAME(libint2_static_init)();

    // ---- precompute all prereqs for NQ quartets (outside timed region) ----
    std::vector<ScalarPack> S(NQ);
    std::vector<std::vector<double>> KF(NQ);
    std::vector<double> out(cl->nout);
    for(int i=0;i<NQ;i++){
        double A[3],B[3],C[3],D[3],a,b,g,d; gen_quartet(i,A,B,C,D,a,b,g,d);
        KF[i].resize(max_m+1);
        recursum_scalars(a,b,g,d,A,B,C,D,S[i],KF[i].data(),max_m);
    }

    volatile double sink=0.0;
    if(impl=="recursum" || impl=="recursum_naive"){
        recursum_fn_t fn=recursum_fn(cl->name, impl=="recursum_naive");
        if(!fn){ fprintf(stderr,"no %s variant for %s (skipped)\n",impl.c_str(),cl->name); return 3; }
        for(int r=0;r<REPS;r++) for(int i=0;i<NQ;i++){
            fn(S[i],KF[i].data(),out.data()); sink+=out[0];
        }
    } else if(impl=="libint2"){
        // build NQ Libint_t evals (prereqs prepped outside timed loop)
        uint am[4]={(uint)cl->l[0],(uint)cl->l[1],(uint)cl->l[2],(uint)cl->l[3]};
        int lmax=LIBINT2_MAX_AM_eri;
        std::vector<Libint_t*> ev(NQ);
        for(int i=0;i<NQ;i++){
            double A[3],B[3],C[3],D[3],a,b,g,d; gen_quartet(i,A,B,C,D,a,b,g,d);
            RandomShellSet<4u> rs(am,1,1);
            for(int k=0;k<3;k++){rs.R[0][k]=A[k];rs.R[1][k]=B[k];rs.R[2][k]=C[k];rs.R[3][k]=D[k];}
            rs.exp[0][0][0]=a;rs.exp[1][0][0]=b;rs.exp[2][0][0]=g;rs.exp[3][0][0]=d;
            rs.coef[0][0][0]=1;rs.coef[1][0][0]=1;rs.coef[2][0][0]=1;rs.coef[3][0][0]=1;
            ev[i]=new Libint_t; LIBINT2_PREFIXED_NAME(libint2_init_eri)(ev[i],lmax,0);
            prep_libint2(ev[i],rs,0,0); ev[i]->contrdepth=1;
        }
        auto fn=LIBINT2_PREFIXED_NAME(libint2_build_eri)[am[0]][am[1]][am[2]][am[3]];
        for(int r=0;r<REPS;r++) for(int i=0;i<NQ;i++){
            fn(ev[i]); sink+=ev[i]->targets[0][0];
        }
    } else { fprintf(stderr,"bad impl\n"); return 2; }
    if(sink==12345.6789) printf("x");   // prevent DCE
    return 0;
}
