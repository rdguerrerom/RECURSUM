# One-electron Hartree–Fock integral stack (Obara–Saika DAGs)

Extends RECURSUM's DAG-topological code generation from the four-centre ERI
(`../recursum_dag_emit.py`, `../eri_dag/`) to the **one-electron integrals a
Hartree–Fock build also needs**: overlap **S**, kinetic energy **T**, and
nuclear attraction **V**. Together with the ERI/`J`/`K` machinery one directory
up, this gives the full HF integral set from a single Obara–Saika DAG framework.

All three are the Head-Gordon–Pople / Obara–Saika (HGP-OS) scheme, sharing the
same horizontal transfer (HRR, A→B) as the ERI:

| integral | recurrence | base case |
|---|---|---|
| **S** overlap | VRR on A `(a+1_i,0)=PA_i(a,0)+a_i/2ζ (a-1_i,0)`, then HRR A→B | `(0,0)=S₀₀=K(π/ζ)^{3/2}` |
| **V** nuclear | VRR on A with the P−C term and Boys index m, then HRR A→B | `(0,0)^{(m)}=−2π/ζ·K·Z·F_m(ζ\|P−C\|²)` |
| **T** kinetic | a fixed β-weighted combination of overlap nodes (no own recurrence) | — (built from S) |

## What's here

| File | Role |
|---|---|
| `onee_dag/` | The one-electron OS DAG — `node.py` (`OneE` node), `recurrence.py` (VRR-A overlap/nuclear, shared HRR, kinetic combination, multipole VRR), `dag.py` (backward DAG + numeric evaluator). Ported from the `osx-gpu` `onee` codegen; the sibling of `../eri_dag/`. |
| `recursum_onee_emit.py` | **The emitter.** Generates C++ `recursum_ovlp_<ab>`, `recursum_nuc_<ab>`, `recursum_kin_<ab>` kernels (ss…ff) with global-CSE topological emission + the same peak-liveness slot allocator as the ERI emitter. Writes `onee_kernels.cpp`, `onee_decls.h`, `recursum_onee_scalars.h`. |
| `onee_oracle.py` | Builds full contracted-Cartesian S/T/V matrices **from the DAG** and validates them against PySCF `int1e_ovlp/kin/nuc_cart`. |
| `onee_export_expected.py` + `onee_validate.cpp` | Numeric check of the **emitted C++ kernels** against the PySCF-validated DAG blocks, per class. |

## Validation

```bash
python3 onee_oracle.py            # DAG vs PySCF int1e_*_cart (up to f-shells)
#   H2O/6-31G   S=8.9e-16 T=8.9e-16 V=9.2e-16  PASS
#   NH/cc-pVDZ  S=2.7e-16 T=2.8e-15 V=2.3e-15  PASS
#   CO/cc-pVTZ  S=2.4e-15 T=4.8e-15 V=6.5e-15  PASS   (s,p,d,f)

python3 recursum_onee_emit.py     # emit kernels (ss..ff) + print DAG sizes
python3 onee_export_expected.py   # export PySCF-validated expected blocks
g++ -O2 -std=c++17 onee_validate.cpp onee_kernels.cpp -o onee_validate && ./onee_validate
#   emitted kernels reproduce the DAG bit-exactly (relerr 0) for S, T, V, ss..ff
```

Kernels consume precomputed per-primitive-pair scalars (a `OneEScalars` struct:
PA, PC, AB, 1/2ζ, β) plus the base value (`S₀₀` for S/T, `kf[m]=K·F_m` for V, one
call per nucleus), exactly as the ERI kernels consume their scalars + `kf[m]`.
Normalization follows libcint/PySCF `_cart` (per-shell `norm_l` for L≤1,
`gto_norm` for L≥2). Generated `.cpp`/`.h` are git-ignored (regenerable).
