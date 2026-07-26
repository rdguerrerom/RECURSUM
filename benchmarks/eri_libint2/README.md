# RECURSUM vs libint2 — primitive Cartesian ERI benchmark

Apples-to-apples comparison of RECURSUM-generated two-electron repulsion integral
(ERI) kernels against **libint2** on the Head-Gordon–Pople / Obara–Saika (HGP-OS)
recurrence. This is the code and data behind the ERI comparison in the manuscript.

**Result:** on 14 canonical Cartesian shell classes (s–f), under an identical
compiler and flags, RECURSUM's DAG-topological kernels match or exceed libint2 on
12 of 14 classes (medians of 20 reps in `comparison_table.txt`), agreeing to a
relative error of ≤3×10⁻¹¹, with the mechanism confirmed by hardware counters
(`perf_results.txt`). See `REPORT.md` for the full methodology and findings.

## What's here

| File | Role |
|---|---|
| `recursum_dag_emit.py` | **The DAG-topological emitter** — expands a recurrence into its integral DAG, emits each node once (global CSE), with peak-liveness slot reuse and a bounded-chunk splitter for very large kernels. Generates the `k_*.cpp` kernels. |
| `eri_dag/` | HGP-OS DAG: `integral.py`, `recurrence.py` (VRR/HRR rules), `dag.py`, `chunker.py`. |
| `oracle.py`, `pyscf_reference.py` | Correctness oracle (evaluates the DAG) + PySCF/libcint cross-check. |
| `bench_eri.cpp` | Google Benchmark harness: RECURSUM kernels vs `libint2_build_eri`. |
| `perf_driver.cpp`, `run_perf.sh`, `stack_frames.py` | `perf` instrumentation + stack-frame measurement. |
| `validate.py`, `validate_kernels.cpp` | Kernel-vs-oracle validation. |
| `Makefile`, `build_and_run.sh` | Build (kernels compiled in parallel) + run. |
| `*.txt`, `REPORT.md` | Result artifacts and the written report. |

The generated `k_*.cpp` kernels and the libint2 build tree are **not** committed
(regenerable / large); see below.

## Reproduce

Prerequisites: `python3` (numpy, scipy, pyscf), `g++` (C++17), Google Benchmark,
and **libint2 built from source with matched flags** (do *not* use a package
binary — its flags are unknown). We used the pre-generated libint2 2.8.1 export
and rebuilt it with `-O3 -march=native -ffast-math`.

```bash
python3 recursum_dag_emit.py     # generate k_*.cpp + eri_classes.h + decls
python3 validate.py              # kernels == oracle (exact)
make -j                          # compile kernels (parallel) + binaries
taskset -c <P-core> ./bench_eri --benchmark_repetitions=20 \
    --benchmark_report_aggregates_only=false \
    --benchmark_out=results_eri.json --benchmark_out_format=json
python3 analyze.py results_eri.json     # comparison table
CORE=<P-core> bash run_perf.sh          # perf instrumentation
python3 stack_frames.py                 # liveness ON vs naive stack frames
```

Machine used: Intel i7-14700 (Raptor Lake, **no AVX-512** → both libraries run
scalar, libint2 `VECLEN=1`), timing pinned to one P-core.

## Third-party

`prep_libint2.h` is copied from the libint2 test suite (Valeev group,
https://github.com/evaleev/libint) and is used unmodified to populate the
`Libint_t` prerequisites so both libraries consume identical inputs. It is
distributed under libint2's LGPL-3.0 license; see the libint2 project for terms.
All other files here are part of RECURSUM (MPL-2.0).
