"""Contracted-shell J/K via the DAG Obara-Saika primitive ERIs + BLAS digestion.

Correctness prototype:
  primitive ERIs   -> oracle.eval_dag (the same recurrence the C kernels compute)
  contraction      -> libint2-style: loop primitive quartets, weight by
                      contraction coeff x radial norm, accumulate the block
  J/K digestion    -> BLAS (numpy .dot / einsum, i.e. DGEMM) with the density
Validated against PySCF int2e_cart and get_jk on the same density.
"""
from __future__ import annotations
import numpy as np
from pyscf import gto, scf

from oracle import Quartet, eval_dag
from pyscf_reference import radial_norm
from eri_dag.integral import all_ang_tuples


def get_shells(mol):
    """List of shells: (l, center_bohr(3,), exps(np), coeffs(np), ao_off, ncart)."""
    sh = []
    off = 0
    for i in range(mol.nbas):
        l = mol.bas_angular(i)
        c = mol.atom_coord(mol.bas_atom(i))            # bohr
        exps = mol.bas_exp(i)
        coeffs = mol.bas_ctr_coeff(i).ravel()
        ncart = (l + 1) * (l + 2) // 2
        sh.append((l, np.asarray(c), np.asarray(exps), np.asarray(coeffs), off, ncart))
        off += ncart
    return sh, off


def block(shA, shB, shC, shD):
    """Contracted Cartesian ERI block (ncA,ncB,ncC,ncD) for a shell quartet."""
    lA, A, eA, cA, _, ncA = shA
    lB, B, eB, cB, _, ncB = shB
    lC, C, eC, cC, _, ncC = shC
    lD, D, eD, cD, _, ncD = shD
    out = np.zeros((ncA, ncB, ncC, ncD))
    # per-primitive weight = contraction coeff x radial norm
    wA = cA * np.array([radial_norm(lA, a) for a in eA])
    wB = cB * np.array([radial_norm(lB, a) for a in eB])
    wC = cC * np.array([radial_norm(lC, a) for a in eC])
    wD = cD * np.array([radial_norm(lD, a) for a in eD])
    for pa, a in enumerate(eA):
        for pb, b in enumerate(eB):
            for pc, g in enumerate(eC):
                for pd, d in enumerate(eD):
                    q = Quartet(a=a, A=tuple(A), b=b, B=tuple(B),
                                g=g, C=tuple(C), d=d, D=tuple(D))
                    _, vals = eval_dag(lA, lB, lC, lD, q)   # canonical (a,b,c,d) order
                    w = wA[pa] * wB[pb] * wC[pc] * wD[pd]
                    out += w * np.asarray(vals).reshape(ncA, ncB, ncC, ncD)
    return out


def build_eri(mol):
    """Full (nao,nao,nao,nao) Cartesian ERI tensor from the DAG primitives."""
    sh, nao = get_shells(mol)
    eri = np.zeros((nao, nao, nao, nao))
    for shA in sh:
        oa, na = shA[4], shA[5]
        for shB in sh:
            ob, nb = shB[4], shB[5]
            for shC in sh:
                oc, nc = shC[4], shC[5]
                for shD in sh:
                    od, nd = shD[4], shD[5]
                    blk = block(shA, shB, shC, shD)
                    eri[oa:oa+na, ob:ob+nb, oc:oc+nc, od:od+nd] = blk
    return eri, nao


def jk_from_eri(eri, D):
    """J/K by BLAS contraction (einsum -> BLAS GEMM). Chemists' notation
    (μν|λσ): J_μν = Σ (μν|λσ) D_σλ ; K_μν = Σ (μλ|νσ) D_σλ."""
    J = np.einsum('pqrs,rs->pq', eri, D, optimize=True)
    K = np.einsum('prqs,rs->pq', eri, D, optimize=True)
    return J, K


if __name__ == "__main__":
    mol = gto.M(atom="C 0 0 0; H 0.629 0.629 0.629; H -0.629 -0.629 0.629; "
                     "H 0.629 -0.629 -0.629; H -0.629 0.629 -0.629",
                basis="sto-3g", unit="Angstrom", cart=True, verbose=0)
    eri_ref = mol.intor('int2e_cart')
    eri_me, nao = build_eri(mol)
    rel = np.max(np.abs(eri_me - eri_ref) / (np.abs(eri_ref).max()))
    print(f"CH4/STO-3G: nao={nao}  max rel err (my ERI vs PySCF int2e_cart) = {rel:.2e}")

    mf = scf.RHF(mol); D = mf.get_init_guess()
    Jr, Kr = scf.hf.get_jk(mol, D)
    Jm, Km = jk_from_eri(eri_me, D)
    jrel = np.max(np.abs(Jm - Jr) / np.abs(Jr).max())
    krel = np.max(np.abs(Km - Kr) / np.abs(Kr).max())
    print(f"J max rel err vs PySCF get_jk = {jrel:.2e}")
    print(f"K max rel err vs PySCF get_jk = {krel:.2e}")
    print("PASS" if max(rel, jrel, krel) < 1e-10 else "FAIL")
