"""RECURSUM DAG-topological emitter for the Hartree-Fock ONE-ELECTRON integral
stack: overlap S, kinetic T, nuclear attraction V. Sibling of
recursum_dag_emit.py (four-centre ERIs); reuses the same peak-liveness slot
allocator and global-CSE topological emission, on the one-electron OS DAG
(onee_dag).

Kernels receive precomputed per-primitive-pair scalars (a `OneEScalars` struct)
plus the base value(s):
  * overlap  recursum_ovlp_<ab>(s, S00, out)      base (0,0) = S00
  * nuclear  recursum_nuc_<ab>(s, kf, out)        base (0,0)^m = kf[m]  (per nucleus)
  * kinetic  recursum_kin_<ab>(s, S00, out)       built from overlap nodes + beta

Both overlap and nuclear share the HRR A->B; nuclear adds the P-C term and Boys
index m; kinetic is a fixed beta-weighted combination of overlap nodes (no own
recurrence). Validated against PySCF int1e_*_cart via onee_oracle.py.
"""
from __future__ import annotations
from typing import Dict, List, Tuple

from onee_dag import OneE, build_dag, kinetic_terms
from onee_dag.dag import DAG
from onee_dag import all_ang_tuples

_CH = "spdfghik"


def _pair_name(la, lb) -> str:
    return f"{_CH[la]}{_CH[lb]}"


def _outputs(la, lb) -> List[OneE]:
    return [OneE(a=a, b=b, m=0) for a in all_ang_tuples(la) for b in all_ang_tuples(lb)]


# ---- coefficient string -> C expression -----------------------------------
def coef_to_c(coef: str) -> str:
    if coef == "one":
        return "1.0"
    return f"s.{coef}"          # PAx.. PCx.. ABx.. inv_2zeta


def node_rhs(n: OneE, dag: DAG, reffn, base_c) -> str:
    terms = dag.nodes[n]
    if not terms:                       # base node
        return base_c(n)
    pieces = []
    for t in terms:
        c = coef_to_c(t.coef)
        if t.imul != 1:
            c = f"{t.imul}.0" if c == "1.0" else f"{t.imul}.0 * {c}"
        r = reffn(t.source)
        sign = "-" if t.sign < 0 else "+"
        pieces.append(f"{sign} {r}" if c == "1.0" else f"{sign} {c} * {r}")
    expr = " ".join(pieces)
    return expr[2:] if expr.startswith("+ ") else expr


# ---- peak-liveness slot allocator (port of recursum_dag_emit.liveness_slots) --
def liveness_slots(topo, dag, outputs, reuse=True):
    pos = {n: i for i, n in enumerate(topo)}
    last_use: Dict[OneE, int] = {}
    for n in topo:
        for t in dag.nodes[n]:
            last_use[t.source] = max(last_use.get(t.source, -1), pos[n])
    die_at: Dict[int, List[OneE]] = {}
    for v, p in last_use.items():
        die_at.setdefault(p, []).append(v)
    slot, free, nxt = {}, [], 0
    for i, n in enumerate(topo):
        if reuse:
            for v in die_at.get(i - 1, []):
                if v in slot:
                    free.append(slot[v])
        if n in outputs or not dag.nodes[n] or n not in last_use:
            continue
        if reuse and free:
            slot[n] = free.pop()
        else:
            slot[n] = nxt; nxt += 1
    return slot, nxt


# ---- emit an overlap OR nuclear kernel (they differ only in base + scalars) --
def _emit_vrr_kernel(la, lb, op, fn_prefix, base_c, sig_base, comment) -> Tuple[str, int, int, int]:
    outs = _outputs(la, lb)
    outset = set(outs)
    oidx = {o: i for i, o in enumerate(outs)}
    dag = build_dag(outs, op=op)
    topo = dag.topological_order()
    slot, peak = liveness_slots(topo, dag, outset)
    name = _pair_name(la, lb)

    def ref(src):
        if src in oidx:
            return f"out[{oidx[src]}]"
        if not dag.nodes[src]:
            return base_c(src)
        return f"sc[{slot[src]}]"

    L = [f"// {comment} ({name}): {len(outs)} Cartesian integrals, "
         f"{len(dag.nodes)} DAG nodes, peak liveness {peak} slots",
         f"__attribute__((noinline)) void {fn_prefix}_{name}(",
         f"    const OneEScalars& __restrict__ s, {sig_base},",
         "    double* __restrict__ out) {"]
    if peak > 0:
        L.append(f"    double sc[{peak}];")
    for n in topo:
        if not dag.nodes[n] and n not in outset:
            continue
        rhs = node_rhs(n, dag, ref, base_c)
        if n in outset:
            L.append(f"    out[{oidx[n]}] = {rhs};")
        elif n in slot:
            L.append(f"    sc[{slot[n]}] = {rhs};")
    L.append("}")
    return "\n".join(L), len(outs), len(dag.nodes), peak


def emit_overlap(la, lb):
    return _emit_vrr_kernel(la, lb, "overlap", "recursum_ovlp",
                            base_c=lambda n: "S00",
                            sig_base="double S00",
                            comment="Overlap S")


def emit_nuclear(la, lb):
    return _emit_vrr_kernel(la, lb, "nuclear", "recursum_nuc",
                            base_c=lambda n: f"kf[{n.m}]",
                            sig_base="const double* __restrict__ kf",
                            comment="Nuclear attraction V (one nucleus)")


# ---- kinetic: overlap DAG (expanded b-set) as intermediates + T combination --
def emit_kinetic(la, lb) -> Tuple[str, int, int, int]:
    ca, cb = all_ang_tuples(la), all_ang_tuples(lb)
    outs = [OneE(a=a, b=b, m=0) for a in ca for b in cb]
    oidx = {o: i for i, o in enumerate(outs)}
    # overlap nodes the T-combination consumes
    need = set()
    for a in ca:
        for b in cb:
            for kt in kinetic_terms(a, b):
                need.add(kt.source)
    dag = build_dag(list(need), op="overlap")
    topo = dag.topological_order()
    # every overlap node is an intermediate here (outputs are T, not S), so give
    # each a slot via liveness against the T consumers (last_use over kinetic_terms)
    pos = {n: i for i, n in enumerate(topo)}
    # S-node -> its own slot (kept until the T-combination reads it at the end)
    snode = {n: i for i, n in enumerate(topo)}
    name = _pair_name(la, lb)

    def ref(src):
        if not dag.nodes[src]:
            return "S00"
        return f"S[{snode[src]}]"

    L = [f"// Kinetic T ({name}): {len(outs)} integrals, {len(dag.nodes)} overlap DAG nodes",
         f"__attribute__((noinline)) void recursum_kin_{name}(",
         "    const OneEScalars& __restrict__ s, double S00,",
         "    double* __restrict__ out) {",
         f"    const double beta = s.beta;",
         f"    double S[{len(topo)}];"]
    for n in topo:
        if not dag.nodes[n]:
            L.append(f"    S[{snode[n]}] = S00;")
        else:
            L.append(f"    S[{snode[n]}] = {node_rhs(n, dag, ref, lambda x: 'S00')};")
    # T combination
    for a in ca:
        for b in cb:
            terms = []
            for kt in kinetic_terms(a, b):
                s_ref = "S00" if not dag.nodes[kt.source] else f"S[{snode[kt.source]}]"
                if kt.kind == "beta_2b1":
                    coef = f"beta * {2*kt.bi+1}.0"
                elif kt.kind == "beta2":
                    coef = "2.0 * beta * beta"
                else:  # half_bb1
                    coef = f"{0.5*kt.bi*(kt.bi-1)}"
                sign = "-" if kt.sign < 0 else "+"
                terms.append(f"{sign} ({coef}) * {s_ref}")
            expr = " ".join(terms)
            expr = expr[2:] if expr.startswith("+ ") else expr
            L.append(f"    out[{oidx[OneE(a=a,b=b,m=0)]}] = {expr};")
    L.append("}")
    return "\n".join(L), len(outs), len(dag.nodes), len(topo)


# ---- scalars header --------------------------------------------------------
SCALARS_H = """#pragma once
// AUTO-GENERATED by recursum_onee_emit.py.  One-electron OS scalars.
struct OneEScalars {
    double PAx, PAy, PAz;    // P - A
    double PCx, PCy, PCz;    // P - C   (nuclear only)
    double ABx, ABy, ABz;    // A - B   (HRR)
    double inv_2zeta;        // 1 / (2 zeta)
    double beta;             // exponent on B (kinetic only)
};
"""

# canonical (la >= lb) pairs up to f
LADDER = [(0,0),(1,0),(1,1),(2,0),(2,1),(2,2),(3,0),(3,1),(3,2),(3,3)]


def generate_project(pairs=LADDER, outdir="."):
    import os
    with open(os.path.join(outdir, "recursum_onee_scalars.h"), "w") as f:
        f.write(SCALARS_H)
    parts = ['#include "recursum_onee_scalars.h"']
    decls = ["#pragma once", '#include "recursum_onee_scalars.h"']
    for la, lb in pairs:
        for emit, pref, sig in ((emit_overlap, "recursum_ovlp", "double"),
                                (emit_nuclear, "recursum_nuc", "const double*"),
                                (emit_kinetic, "recursum_kin", "double")):
            src = emit(la, lb)[0]
            parts.append(src)
            nm = f"{pref}_{_pair_name(la,lb)}"
            decls.append(f"void {nm}(const OneEScalars&, {sig}, double*);")
    with open(os.path.join(outdir, "onee_kernels.cpp"), "w") as f:
        f.write("\n\n".join(parts) + "\n")
    with open(os.path.join(outdir, "onee_decls.h"), "w") as f:
        f.write("\n".join(decls) + "\n")
    with open(os.path.join(outdir, "onee_dispatch.h"), "w") as f:
        f.write(emit_dispatch(pairs))
    return len(pairs)


def emit_dispatch(pairs=LADDER) -> str:
    """(la,lb) -> kernel function pointer, for the contracted driver. Only
    canonical la>=lb classes exist; the driver swaps+transposes for la<lb."""
    def cases(pref):
        return "\n".join(f"    case {la*8+lb}: return {pref}_{_pair_name(la,lb)};"
                         for la, lb in pairs)
    return f"""#pragma once
// AUTO-GENERATED (la,lb) kernel dispatch. DO NOT EDIT.
#include "recursum_onee_scalars.h"
#include "onee_decls.h"
typedef void (*ovlp_fn)(const OneEScalars&, double, double*);
typedef void (*nuc_fn)(const OneEScalars&, const double*, double*);
typedef void (*kin_fn)(const OneEScalars&, double, double*);
static inline ovlp_fn onee_ovlp_dispatch(int la, int lb) {{
  switch (la*8+lb) {{
{cases("recursum_ovlp")}
    default: return nullptr;
  }}
}}
static inline nuc_fn onee_nuc_dispatch(int la, int lb) {{
  switch (la*8+lb) {{
{cases("recursum_nuc")}
    default: return nullptr;
  }}
}}
static inline kin_fn onee_kin_dispatch(int la, int lb) {{
  switch (la*8+lb) {{
{cases("recursum_kin")}
    default: return nullptr;
  }}
}}
"""


if __name__ == "__main__":
    n = generate_project()
    print(f"generated overlap/kinetic/nuclear kernels for {n} shell-pair classes "
          f"-> onee_kernels.cpp + onee_decls.h + recursum_onee_scalars.h")
    print(f"{'pair':>5} {'nout':>5} {'S_nodes':>8} {'S_peak':>7} {'V_nodes':>8} {'V_peak':>7}")
    for la, lb in LADDER:
        _, no, sn, sp = emit_overlap(la, lb)
        _, _, vn, vp = emit_nuclear(la, lb)
        print(f"{_pair_name(la,lb):>5} {no:>5} {sn:>8} {sp:>7} {vn:>8} {vp:>7}")
