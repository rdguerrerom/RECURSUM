"""Export alkane basis + density + reference J/K for the C++ contracted-J/K driver.

Geometry: carbons on a line (1.54 A), hydrogens placed tetrahedrally per carbon
(3 on terminal C, 2 on internal). Any non-degenerate geometry is fine for a
timing/scaling benchmark; integrals are well defined for any set of Gaussians.
Basis: 6-31G (s,p only for C,H -> the existing s,p DAG kernels suffice).
File format (text): per shell (l cx cy cz nprim exps... coeffs...); then nao;
then D (row-major); then reference J and K (row-major).
"""
import numpy as np
from pyscf import gto, scf

BOHR = 1.8897259886


def alkane(n_c):
    """CnH2n+2 cartesian (Angstrom)."""
    atoms = []
    ccx = 1.54
    ch = 1.09
    # tetrahedral-ish H directions
    dirs = np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]], float)
    dirs /= np.linalg.norm(dirs[0])
    for i in range(n_c):
        cx = i * ccx
        atoms.append(("C", (cx, 0.0, 0.0)))
        # hydrogens: CH4 gets 4; else terminal C get 3, internal get 2 (CnH2n+2)
        if n_c == 1:
            hd = dirs
        elif i == 0:
            hd = dirs[[0,1,2]]
        elif i == n_c-1:
            hd = dirs[[1,2,3]]
        else:
            hd = dirs[[2,3]] if i % 2 else dirs[[0,1]]
        for d in hd:
            atoms.append(("H", (cx + ch*d[0], ch*d[1], ch*d[2])))
    return atoms


def export(n_c, path):
    atoms = alkane(n_c)
    mol = gto.M(atom=[[a, tuple(c)] for a, c in atoms], basis="6-31g",
                unit="Angstrom", cart=True, verbose=0)
    nao = mol.nao_nr()
    mf = scf.RHF(mol); D = mf.get_init_guess()
    J, K = scf.hf.get_jk(mol, D)
    with open(path, "w") as f:
        f.write(f"{mol.nbas}\n")
        for i in range(mol.nbas):
            l = mol.bas_angular(i); c = mol.atom_coord(mol.bas_atom(i))  # bohr
            e = mol.bas_exp(i); co = mol.bas_ctr_coeff(i).ravel()
            f.write(f"{l} {c[0]:.12g} {c[1]:.12g} {c[2]:.12g} {len(e)} "
                    + " ".join(f"{x:.12g}" for x in e) + " "
                    + " ".join(f"{x:.12g}" for x in co) + "\n")
        f.write(f"{nao}\n")
        for M in (D, J, K):
            f.write(" ".join(f"{x:.15g}" for x in M.ravel()) + "\n")
    maxL = max(mol.bas_angular(i) for i in range(mol.nbas))
    return mol.nbas, nao, maxL


if __name__ == "__main__":
    for n in (1, 2, 3, 4):
        nbas, nao, maxL = export(n, f"alkane_C{n}.txt")
        print(f"C{n}H{2*n+2}: {nbas} shells, {nao} AOs, maxL={maxL} -> alkane_C{n}.txt")
