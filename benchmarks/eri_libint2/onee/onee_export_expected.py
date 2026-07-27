"""Export a fixed test primitive-pair's OneEScalars + expected S/T/V blocks
(from the PySCF-validated onee_dag) as a C header, so onee_validate.cpp can
check the *emitted C++ kernels* numerically."""
import math
import numpy as np
from onee_dag import OneE, build_dag, evaluate, kinetic_terms
from onee_dag.dag import DAG
from onee_dag import all_ang_tuples
from onee_oracle import boys_array
from recursum_onee_emit import LADDER, _pair_name

# test primitive pair + nucleus
al, A = 1.2, np.array([0.0, 0.0, 0.0])
be, B = 0.8, np.array([0.5, 0.1, -0.2])
Z, C = 3.0, np.array([0.3, -0.4, 0.6])

zeta = al + be
P = (al * A + be * B) / zeta
AB = A - B
K = math.exp(-al * be / zeta * float(np.dot(AB, AB)))
PA = P - A
PC = P - C
inv2z = 1.0 / (2.0 * zeta)
S00 = K * (math.pi / zeta) ** 1.5
Targ = zeta * float(np.dot(PC, PC))
NMAX = max(la + lb for la, lb in LADDER)
F = boys_array(NMAX, Targ)
kf = -2.0 * math.pi / zeta * K * Z * F

sc = {"PAx": PA[0], "PAy": PA[1], "PAz": PA[2],
      "PCx": PC[0], "PCy": PC[1], "PCz": PC[2],
      "ABx": AB[0], "ABy": AB[1], "ABz": AB[2],
      "inv_2zeta": inv2z, "one": 1.0}


def ovlp_block(la, lb):
    outs = [OneE(a=a, b=b, m=0) for a in all_ang_tuples(la) for b in all_ang_tuples(lb)]
    v = evaluate(build_dag(outs, op="overlap"), sc, base_eval=lambda n: S00)
    return [v[o] for o in outs]


def nuc_block(la, lb):
    outs = [OneE(a=a, b=b, m=0) for a in all_ang_tuples(la) for b in all_ang_tuples(lb)]
    v = evaluate(build_dag(outs, op="nuclear"), sc, base_eval=lambda n: kf[n.m])
    return [v[o] for o in outs]


def kin_block(la, lb):
    ca, cb = all_ang_tuples(la), all_ang_tuples(lb)
    need = set()
    for a in ca:
        for b in cb:
            for kt in kinetic_terms(a, b):
                need.add(kt.source)
    S = evaluate(build_dag(list(need), op="overlap"), sc, base_eval=lambda n: S00)
    out = []
    for a in ca:
        for b in cb:
            t = 0.0
            for kt in kinetic_terms(a, b):
                if kt.kind == "beta_2b1":
                    coef = be * (2 * kt.bi + 1)
                elif kt.kind == "beta2":
                    coef = 2.0 * be * be
                else:
                    coef = 0.5 * kt.bi * (kt.bi - 1)
                t += kt.sign * coef * S[kt.source]
            out.append(t)
    return out


def carr(name, vals):
    return f"static const double {name}[] = {{{', '.join(f'{x:.17g}' for x in vals)}}};"


L = ["#pragma once", "// AUTO-GENERATED expected one-electron blocks (PySCF-validated DAG).",
     f"static const OneEScalars TEST_S = {{{PA[0]:.17g},{PA[1]:.17g},{PA[2]:.17g},"
     f"{PC[0]:.17g},{PC[1]:.17g},{PC[2]:.17g},{AB[0]:.17g},{AB[1]:.17g},{AB[2]:.17g},"
     f"{inv2z:.17g},{be:.17g}}};",
     f"static const double TEST_S00 = {S00:.17g};",
     carr("TEST_KF", kf)]
for la, lb in LADDER:
    nm = _pair_name(la, lb)
    L.append(carr(f"EXP_ovlp_{nm}", ovlp_block(la, lb)))
    L.append(carr(f"EXP_nuc_{nm}", nuc_block(la, lb)))
    L.append(carr(f"EXP_kin_{nm}", kin_block(la, lb)))
open("onee_expected.h", "w").write("\n".join(L) + "\n")
print("wrote onee_expected.h for", len(LADDER), "classes; NMAX =", NMAX)
