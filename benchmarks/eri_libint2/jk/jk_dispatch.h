#pragma once
#include "recursum_eri_scalars.h"
#include "jk_kernels.h"
typedef void(*rfn)(const ScalarPack&,const double*,double*);
static inline rfn jk_dispatch(int la,int lb,int lc,int ld){
  if(la==0&&lb==0&&lc==0&&ld==0) return recursum_eri_ssss;
  if(la==0&&lb==0&&lc==0&&ld==1) return recursum_eri_sssp;
  if(la==0&&lb==0&&lc==1&&ld==0) return recursum_eri_ssps;
  if(la==0&&lb==0&&lc==1&&ld==1) return recursum_eri_sspp;
  if(la==0&&lb==1&&lc==0&&ld==0) return recursum_eri_spss;
  if(la==0&&lb==1&&lc==0&&ld==1) return recursum_eri_spsp;
  if(la==0&&lb==1&&lc==1&&ld==0) return recursum_eri_spps;
  if(la==0&&lb==1&&lc==1&&ld==1) return recursum_eri_sppp;
  if(la==1&&lb==0&&lc==0&&ld==0) return recursum_eri_psss;
  if(la==1&&lb==0&&lc==0&&ld==1) return recursum_eri_pssp;
  if(la==1&&lb==0&&lc==1&&ld==0) return recursum_eri_psps;
  if(la==1&&lb==0&&lc==1&&ld==1) return recursum_eri_pspp;
  if(la==1&&lb==1&&lc==0&&ld==0) return recursum_eri_ppss;
  if(la==1&&lb==1&&lc==0&&ld==1) return recursum_eri_ppsp;
  if(la==1&&lb==1&&lc==1&&ld==0) return recursum_eri_ppps;
  if(la==1&&lb==1&&lc==1&&ld==1) return recursum_eri_pppp;
  return nullptr; }