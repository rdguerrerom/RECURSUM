#!/usr/bin/env bash
# Direct perf instrumentation backing the liveness/cache-friendliness claim.
# For each class and impl (recursum / recursum_naive / libint2), stream NQ
# distinct quartets REPS times through the kernel under `perf stat`, pinned to
# one P-core. Reports cycles, instructions, IPC, L1-dcache load misses, and
# LLC misses -- the direct evidence for whether peak-liveness slot reuse
# actually reduces memory traffic.
set -uo pipefail
cd "$(dirname "$0")"

CORE=${CORE:-8}                 # a P-core (5400 MHz) per lscpu
NQ=${NQ:-2000}                  # distinct quartets streamed (working set > L1)
REPS=${REPS:-2000}
EVENTS=cycles,instructions,L1-dcache-loads,L1-dcache-load-misses,LLC-loads,LLC-load-misses
CLASSES=${CLASSES:-"pppp dddd ssff fsfs ffff"}
IMPLS="recursum recursum_naive libint2"

echo "# perf: pinned core $CORE, NQ=$NQ quartets x REPS=$REPS, events=$EVENTS"
printf "%-8s %-15s %14s %14s %6s %16s %10s %14s %10s\n" \
  class impl cycles instructions IPC L1dloads L1miss% LLCloads LLCmiss%
for cls in $CLASSES; do
  for impl in $IMPLS; do
    OUT=$(taskset -c $CORE perf stat -x, -e $EVENTS -- ./perf_driver $impl $cls $NQ $REPS 2>&1)
    if echo "$OUT" | grep -q "skipped"; then
      printf "%-8s %-15s %14s\n" "$cls" "$impl" "(no naive variant — ON-only class)"; continue
    fi
    # hybrid CPU: pinned to a P-core, read the cpu_core/<event>/ rows
    get(){ echo "$OUT" | awk -F, -v e="cpu_core/$1/" '$3==e{print $1}'; }
    cyc=$(get cycles); ins=$(get instructions)
    l1l=$(get L1-dcache-loads); l1m=$(get L1-dcache-load-misses)
    lll=$(get LLC-loads); llm=$(get LLC-load-misses)
    ipc=$(awk -v i="$ins" -v c="$cyc" 'BEGIN{if(c>0)printf"%.2f",i/c;else print"-"}')
    l1p=$(awk -v m="$l1m" -v l="$l1l" 'BEGIN{if(l>0)printf"%.2f",100*m/l;else print"-"}')
    llp=$(awk -v m="$llm" -v l="$lll" 'BEGIN{if(l>0)printf"%.2f",100*m/l;else print"-"}')
    printf "%-8s %-15s %14s %14s %6s %16s %6s%% %14s %8s%%\n" \
      "$cls" "$impl" "$cyc" "$ins" "$ipc" "$l1l" "$l1p" "$lll" "$llp"
  done
  echo
done
