# Contracted J/K matrix construction — RECURSUM vs libint2

End-to-end **Coulomb (J)** and **exchange (K)** matrix build using the
DAG-topological Obara–Saika kernels, compared head-to-head against **libint2**
under *identical* contraction and BLAS density digestion. This tests whether the
per-quartet kernel advantage of the primitive benchmark (one directory up)
survives in a full contracted algorithm.

Only the primitive integral kernel differs between the two builds (RECURSUM's
DAG kernels vs libint2's `build_eri`); both are driven through the same primitive
loop, the same Chebyshev Boys evaluation (`FmEval_Chebyshev7`, shared for
fairness), the same `dgemv`/`dgemm` (BLAS) digestion, and — for the larger
systems — the same 8-fold permutational-symmetry driver. Both builds reproduce a
PySCF reference J and K to a relative error of ≤2×10⁻¹².

## Result

**RECURSUM is ~2× faster than libint2 for both J and K**, and the advantage is
not an artifact of the (near-1D) alkane test geometry — it holds on compact 3D
DNA fragments, which simultaneously recover the full theoretical O(N⁴) scaling
(the flat alkanes give a sub-quartic O(N^3.8) fit).

| System | shells | J speedup | K speedup |
|---|---|---|---|
| alkanes CH₄–C₄H₁₀ (6-31G) | 13–40 | ~2.0× | ~2.0× |
| DNA nucleoside (6-31G)     | 116   | 2.07× | 2.20× |
| DNA strand, 2-mer (6-31G)  | 232   | 1.90× | 1.99× |

Nucleoside→strand local scaling exponent: b_J = b_K ≈ 4.13 (theoretical O(N⁴)).
Per-system numbers in `jk_results.txt` (alkanes) and `jk_dna_results.txt` (DNA).

## What's here

| File | Role |
|---|---|
| `jk_driver.cpp` | The timed C++ driver. Both impls; `<recursum\|libint2> mol.txt [reps] [sym]`. `sym` enables the 8-fold permutational-symmetry build (each unique shell quartet evaluated once, all symmetry images digested). Reports J/K time + validation vs the file's PySCF reference. |
| `jk_kernels.h`, `jk_dispatch.h` | The DAG s,p ERI kernels (all (la,lb,lc,ld)∈{0,1}⁴) + dispatch. |
| `jk_contracted.py` | Python correctness prototype (DAG primitives → contraction → BLAS J/K), validated vs PySCF `int2e_cart` / `get_jk`. |
| `export_alkanes.py` | Builds CₙH₂ₙ₊₂ (6-31G, cart), writes shells + density + reference J/K. |
| `export_mol.py` | Generic exporter: any `.xyz` → shells + density + reference J/K (used for the DNA fragments). |
| `run_dna_jk.sh` | Runs the DNA head-to-head (min-of-reps, pinned P-core, symmetric path). |

## Reproduce

Prerequisites: `python3` (numpy, pyscf), `g++` (C++17), OpenBLAS, and **libint2
built from source with matched flags** (see the parent README).

```bash
# 1. export test molecules (writes *.txt: shells + density + reference J/K)
python3 export_alkanes.py                      # alkane_C{1..4}.txt
python3 export_mol.py dna_nucleoside.xyz dna_nucleoside.txt 6-31g
python3 export_mol.py dna_strand_2mer.xyz   dna_2mer.txt       6-31g

# 2. build the driver (adjust libint2 paths)
g++ -O3 -march=native -ffast-math -funroll-loops -std=c++17 -no-pie jk_driver.cpp \
    -I <libint2>/include -I <libint2-build>/include <libint2-build>/libint2.a \
    -lopenblas -o jk_driver

# 3. run (alkanes: full loop; DNA: 8-fold symmetry, min of reps)
taskset -c <P-core> ./jk_driver recursum alkane_C4.txt 20
taskset -c <P-core> ./jk_driver libint2  alkane_C4.txt 20
bash run_dna_jk.sh
```

Screening is disabled throughout, to expose the O(N⁴) cost of unscreened
four-centre evaluation. Machine: Intel i7-14700 (no AVX-512 → both libraries run
scalar), timing pinned to one P-core.
