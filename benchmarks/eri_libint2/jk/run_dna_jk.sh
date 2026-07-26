#!/usr/bin/env bash
# Reliable rerun: min-of-reps, per-rep timing to stderr log. Pinned P-core 8.
set -e
cd /home/ruben/Research/Science/Projects/RECURSUM_Research/eri_libint_comparison
OUT=jk_dna_results2.txt; LOG=jk_dna_perrep.log
: > $OUT; : > $LOG
echo "# DNA J/K, 6-31G cart, 8-fold sym, screening off, reps=2 (min-of-rep), taskset -c 8, $(date -u +%FT%TZ)" >> $OUT
for mol in dna_nucleoside dna_2mer; do
  for impl in recursum libint2; do
    echo "### $impl $mol" >> $LOG
    echo ">>> $impl $mol" >&2
    taskset -c 8 ./jk_driver $impl $mol.txt 2 sym >> $OUT 2>> $LOG
  done
done
echo "ALL DONE" >> $OUT