// Per-Hamiltonian-term timing: RECURSUM one-electron kernels (+ contraction)
// vs libint2, for overlap S, kinetic T, nuclear V. Same molecule/basis; both
// build the full contracted-Cartesian matrix. libint2 side uses its high-level
// Engine (which requires Eigen — a libint2-internal, header-only dependency);
// RECURSUM uses no linear-algebra library at all.
//   g++ -O3 -march=native -ffast-math -std=c++17 onee_bench.cpp onee_kernels.cpp \
//       -I<libint2>/include -I<libint2>/build/include -I/usr/include/eigen3 \
//       <libint2>/build/libint2.a -o onee_bench
//   taskset -c 8 ./onee_bench mol.txt [reps]
#include "onee_common.h"
#include <chrono>
#include <libint2.hpp>
#include <libint2/boys.h>

static double now(){ using namespace std::chrono;
    return duration<double>(steady_clock::now().time_since_epoch()).count(); }

// Shared fast Boys (libint2's Chebyshev) for RECURSUM's nuclear path, so the V
// comparison times recurrence-vs-recurrence, not two different Boys routines
// (exactly the fairness protocol of the ERI benchmark).
static libint2::FmEval_Chebyshev7<double> g_fm(7);
static void boys_cheb(int m,double T,double*F){ g_fm.eval(F,T,m); }

int main(int argc,char**argv){
    if(argc<2){ fprintf(stderr,"usage: %s mol.txt [reps]\n",argv[0]); return 2; }
    int reps = argc>2?atoi(argv[2]):50;
    int nao; vector<double> Z; vector<std::array<double,3>> Rc; vector<double> Sr,Tr,Vr;
    auto sh=read_mol(argv[1],nao,Z,Rc,Sr,Tr,Vr);
    int maxL=0; for(auto&s:sh) maxL=std::max(maxL,s.l);

    // ---- libint2 shells (same exponents/coeffs/centers, cartesian) ----------
    libint2::initialize();
    g_boys = &boys_cheb;   // fair: both sides use libint2's Chebyshev Boys
    std::vector<libint2::Shell> lsh;
    for(auto&s:sh){
        libint2::svector<double> e(s.e.begin(),s.e.end());
        libint2::svector<double> d(s.d.begin(),s.d.end());
        lsh.push_back(libint2::Shell{ e, {{s.l,false,d}}, {{s.c[0],s.c[1],s.c[2]}} });
    }
    size_t maxnp=1; for(auto&s:lsh) maxnp=std::max(maxnp,s.nprim());
    std::vector<std::pair<double,std::array<double,3>>> charges;
    for(size_t k=0;k<Z.size();++k) charges.push_back({Z[k],{Rc[k][0],Rc[k][1],Rc[k][2]}});

    Op ops[3]={OVLP,KIN,NUC};
    libint2::Operator lops[3]={libint2::Operator::overlap,libint2::Operator::kinetic,
                               libint2::Operator::nuclear};
    const char* names[3]={"S (overlap)","T (kinetic)","V (nuclear)"};

    printf("Per-term timing  nao=%d maxL=%d  reps=%d  (ms per full matrix build)\n",
           nao,maxL,reps);
    printf("  %-12s %10s %10s %9s\n","term","RECURSUM","libint2","speedup");
    for(int o=0;o<3;o++){
        vector<double> M((size_t)nao*nao), blk(256);
        // RECURSUM
        double t0=now();
        for(int r=0;r<reps;r++){
            for(auto&SA:sh)for(auto&SB:sh){ if(SB.off>SA.off) continue;
                block(SA,SB,ops[o],Z,Rc,blk.data());
                for(int a=0;a<SA.ncart;a++)for(int b=0;b<SB.ncart;b++){
                    double v=blk[a*SB.ncart+b];
                    M[(size_t)(SA.off+a)*nao+(SB.off+b)]=v;
                    M[(size_t)(SB.off+b)*nao+(SA.off+a)]=v; } }
        }
        double tR=(now()-t0)/reps*1e3;
        // libint2
        libint2::Engine eng(lops[o], maxnp, maxL, 0);
        if(o==2) eng.set_params(charges);
        const auto& buf = eng.results();
        double chk=0, t1=now();
        for(int r=0;r<reps;r++){
            for(size_t i=0;i<lsh.size();++i)for(size_t j=0;j<=i;++j){
                eng.compute(lsh[i],lsh[j]);
                if(buf[0]) chk+=buf[0][0];
            }
        }
        double tL=(now()-t1)/reps*1e3;
        printf("  %-12s %10.4f %10.4f %8.2fx\n",names[o],tR,tL,tL/tR);
        if(chk==12345.6) printf("");   // keep chk live
    }
    libint2::finalize();
    return 0;
}
