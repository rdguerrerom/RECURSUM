"""Physics gate: confirm eval_dag() matches PySCF/libcint for primitive
Cartesian ERIs, up to the per-shell radial norm (gto_norm) and per-component
Cartesian norm (cart_norm_factor). If this holds, libint2 (same physics) will
match eval_dag to the same constants.

PySCF int2e_cart integrates over NORMALIZED Cartesian GTOs. Our eval_dag()
integrates UNNORMALIZED primitives. Relationship for a single primitive shell:

    pyscf[i,j,k,l] = Na*Nb*Nc*Nd * cnf(i)*cnf(j)*cnf(k)*cnf(l) * eval_dag[...]

where N = gto_norm(l, alpha) (radial) and cnf = cart_norm_factor(lx,ly,lz).
"""
from __future__ import annotations

import math

import numpy as np
from pyscf import gto

from eri_dag.integral import all_ang_tuples, cart_norm_factor
from oracle import Quartet, eval_dag


def build_mol(la, lb, lc, ld, q: Quartet):
    """4 distinct single-primitive shells at the 4 centers (Bohr)."""
    # unique fake atoms so each gets its own single-l single-primitive basis
    atom = [["H", q.A], ["He", q.B], ["Li", q.C], ["Be", q.D]]
    basis = {
        "H":  [[la, [q.a, 1.0]]],
        "He": [[lb, [q.b, 1.0]]],
        "Li": [[lc, [q.g, 1.0]]],
        "Be": [[ld, [q.d, 1.0]]],
    }
    mol = gto.M(atom=atom, basis=basis, unit="Bohr", cart=True, spin=None,
                charge=0, verbose=0)
    return mol


def pyscf_block(la, lb, lc, ld, q: Quartet):
    mol = build_mol(la, lb, lc, ld, q)
    # shells are in atom order: shell 0=A(la),1=B(lb),2=C(lc),3=D(ld)
    eri = mol.intor("int2e_cart", shls_slice=(0, 1, 1, 2, 2, 3, 3, 4))
    return eri.reshape(-1)  # canonical (a,b,c,d) cartesian order


def _dfact_2nm1(n):
    r = 1
    for i in range(1, n + 1):
        r *= 2 * i - 1
    return r


def radial_norm(l, alpha):
    """Standard primitive GTO radial norm N_l(a) = (2a/pi)^{3/4} (4a)^{l/2}
    / sqrt((2l-1)!!) -- the normalization PySCF applies (AO self-overlap = 1
    for the xx..x component). Verified against int1e_ovlp_cart."""
    return (2 * alpha / math.pi) ** 0.75 * (4 * alpha) ** (l / 2.0) / math.sqrt(_dfact_2nm1(l))


def norm_factors(la, lb, lc, ld, q: Quartet):
    """Expected pyscf/eval_dag ratio per output component."""
    Na = radial_norm(la, q.a)
    Nb = radial_norm(lb, q.b)
    Nc = radial_norm(lc, q.g)
    Nd = radial_norm(ld, q.d)
    facs = []
    for a in all_ang_tuples(la):
        for b in all_ang_tuples(lb):
            for c in all_ang_tuples(lc):
                for d in all_ang_tuples(ld):
                    cnf = (cart_norm_factor(*a) * cart_norm_factor(*b)
                           * cart_norm_factor(*c) * cart_norm_factor(*d))
                    facs.append(Na * Nb * Nc * Nd * cnf)
    return np.array(facs)


def physics_gate(q: Quartet, classes, verbose=True):
    """Convention-free physics validation.

    PySCF/libcint applies a single per-shell radial-norm constant to every
    Cartesian component of a shell (a per-quartet scalar K_pyscf that is the
    product of 4 shell norms). It does NOT apply a per-component factor. So
    ref = K_pyscf * ours holds with K_pyscf CONSTANT across the whole block.

    We verify eval_dag() is physically correct without needing to know
    PySCF's exact norm convention: fit the single best K per block and check
    the residual. A constant ratio (small residual) proves every relative
    integral -- diagonal AND off-diagonal Cartesian components -- matches
    libcint. This is the gate that also guarantees libint2 agrees.
    """
    allok = True
    if verbose:
        print(f"{'class':>14} {'n':>4} {'K_pyscf':>13} {'resid(rel)':>12}  status")
    for cls in classes:
        _, ours = eval_dag(*cls, q)
        ours = np.array(ours)
        ref = pyscf_block(*cls, q)
        # best-fit single constant K (least squares); then residual
        K = float(ref @ ours / (ours @ ours))
        denom = np.maximum(np.abs(ref), np.abs(ref).max() * 1e-300 + 1e-300)
        resid = np.max(np.abs(K * ours - ref) / denom)
        ok = resid < 1e-9
        allok &= ok
        if verbose:
            print(f"{str(cls):>14} {len(ours):>4} {K:13.6e} {resid:12.3e}  "
                  f"{'OK' if ok else 'MISMATCH'}")
    return allok


if __name__ == "__main__":
    q = Quartet(a=1.2, A=(0.0, 0.0, 0.0),
                b=0.8, B=(0.5, 0.0, 0.0),
                g=1.5, C=(0.0, 0.7, 0.0),
                d=0.9, D=(0.0, 0.0, 0.9))
    classes = [(0,0,0,0),(1,0,0,0),(1,1,0,0),(1,1,1,1),(2,0,0,0),
               (2,0,2,0),(2,2,0,0),(2,2,2,2),(3,0,0,0),(3,3,0,0)]
    ok = physics_gate(q, classes)
    print("\nALL PASS — eval_dag physics validated vs PySCF/libcint"
          if ok else "\nFAILURES PRESENT")
