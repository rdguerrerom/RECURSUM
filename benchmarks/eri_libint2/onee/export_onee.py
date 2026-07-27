"""Export a molecule (shells + nuclei + reference int1e S/T/V) for onee_driver.cpp.
Usage: python3 export_onee.py "<atom>" <basis> <out.txt>"""
import sys
import numpy as np
from pyscf import gto


def export(atom, basis, out):
    mol = gto.M(atom=atom, basis=basis, unit="Bohr", cart=True, verbose=0)
    nao = mol.nao_nr()
    S = mol.intor("int1e_ovlp_cart")
    T = mol.intor("int1e_kin_cart")
    V = mol.intor("int1e_nuc_cart")
    # one emitted shell per (bas, contraction column), libcint AO order
    # (ctr outer, Cartesian inner) so the matrix matches int1e_*_cart.
    rows = []
    for i in range(mol.nbas):
        l = mol.bas_angular(i); c = mol.bas_coord(i)
        e = mol.bas_exp(i); ctr = mol.bas_ctr_coeff(i)   # (nprim, nctr)
        for col in range(ctr.shape[1]):
            co = ctr[:, col]
            rows.append(f"{l} {c[0]:.15g} {c[1]:.15g} {c[2]:.15g} {len(e)} "
                        + " ".join(f"{x:.15g}" for x in e) + " "
                        + " ".join(f"{x:.15g}" for x in co))
    with open(out, "w") as f:
        f.write(f"{len(rows)}\n")
        for r in rows:
            f.write(r + "\n")
        f.write(f"{nao}\n")
        f.write(f"{mol.natm}\n")
        for k in range(mol.natm):
            Z = mol.atom_charge(k); r = mol.atom_coord(k)
            f.write(f"{Z:.15g} {r[0]:.15g} {r[1]:.15g} {r[2]:.15g}\n")
        for M in (S, T, V):
            f.write(" ".join(f"{x:.15g}" for x in M.ravel()) + "\n")
    maxL = max(mol.bas_angular(i) for i in range(mol.nbas))
    print(f"{out}: {mol.natm} atoms, {mol.nbas} shells, {nao} AOs, maxL={maxL}")


if __name__ == "__main__":
    export(sys.argv[1], sys.argv[2], sys.argv[3])
