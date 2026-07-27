#include "onee_common.h"
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
