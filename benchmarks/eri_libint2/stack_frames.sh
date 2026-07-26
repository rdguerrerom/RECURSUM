#!/usr/bin/env bash
# Direct measurement of each kernel's STACK FRAME from the compiled object --
# the concrete realization of the peak-liveness slot count. Reads the
# `sub $N, %rsp` prologue of each recursum_eri_<class>[_naive] function.
cd "$(dirname "$0")"
printf "%-14s %12s %12s %10s\n" kernel frame_ON frame_naive reduce
for cls in pppp dddd ssff fsfs ffff; do
  fon=$(objdump -d k_${cls}.o 2>/dev/null | awk '/<_?Z?.*recursum_eri_'"${cls}"'>:/{f=1} f&&/sub .*0x[0-9a-f]+,%rsp/{gsub(/[$,]/,"",$0); for(i=1;i<=NF;i++) if($i ~ /0x/ && $(i+1) ~ /rsp/){print strtonum($i); exit}}')
  fnv=$(objdump -d k_${cls}_naive.o 2>/dev/null | awk '/recursum_eri_'"${cls}"'_naive>:/{f=1} f&&/sub .*0x[0-9a-f]+,%rsp/{gsub(/[$,]/,"",$0); for(i=1;i<=NF;i++) if($i ~ /0x/ && $(i+1) ~ /rsp/){print strtonum($i); exit}}')
  red=$(awk -v a="$fon" -v b="$fnv" 'BEGIN{if(b>0)printf"%.0f%%",100*(1-a/b);else print"-"}')
  printf "%-14s %12s %12s %10s\n" "$cls" "${fon:-?}" "${fnv:-?}" "$red"
done
