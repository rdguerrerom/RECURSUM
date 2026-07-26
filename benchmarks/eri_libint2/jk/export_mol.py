"""Export any molecule (xyz) + 6-31G basis + density + reference J/K for jk_driver.
Usage: python3 export_mol.py <in.xyz> <out.txt> [basis]"""
import sys
import numpy as np
from pyscf import gto, scf


def export(xyz, out, basis="6-31g"):
    mol = gto.M(atom=xyz, basis=basis, cart=True, verbose=0)
    nao = mol.nao_nr()
    mf = scf.RHF(mol); D = mf.get_init_guess()
    J, K = scf.hf.get_jk(mol, D)
    with open(out, "w") as f:
        f.write(f"{mol.nbas}\n")
        for i in range(mol.nbas):
            l = mol.bas_angular(i); c = mol.atom_coord(mol.bas_atom(i))
            e = mol.bas_exp(i); co = mol.bas_ctr_coeff(i).ravel()
            f.write(f"{l} {c[0]:.12g} {c[1]:.12g} {c[2]:.12g} {len(e)} "
                    + " ".join(f"{x:.12g}" for x in e) + " "
                    + " ".join(f"{x:.12g}" for x in co) + "\n")
        f.write(f"{nao}\n")
        for M in (D, J, K):
            f.write(" ".join(f"{x:.15g}" for x in M.ravel()) + "\n")
    maxL = max(mol.bas_angular(i) for i in range(mol.nbas))
    print(f"{out}: {mol.natm} atoms, {mol.nbas} shells, {nao} AOs, maxL={maxL}")


if __name__ == "__main__":
    export(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "6-31g")
