"""RECURSUM DAG-topological emitter for CARTESIAN MULTIPOLE-MOMENT integrals
(overlap S, dipole D, quadrupole Q) and their geometric gradients.

Sibling of recursum_onee_emit.py (S/T/V). Two things are new here.

1. A FUSED two-layer DAG. The moment VRR (HJO Eq. 9.5.43) climbs the moment
   index e and bottoms out at e == 0 on the OVERLAP value S(a,b), which is
   itself the root of the Obara-Saika overlap DAG. Rather than evaluating the
   two layers separately, the bridge M(a,b,0) -> S(a,b) is made an explicit DAG
   edge, so one memoised graph spans both layers and global CSE applies ACROSS
   them: the overlap sub-graph is built once and shared by every dipole and
   quadrupole output, instead of once per moment component. For a (d,d) class
   that is the difference between one overlap DAG and ten.

2. C99 emission. The generated kernels are straight-line scalar double code with
   no templates, no SIMD types and no allocation, so they are emitted as C99 --
   pointer-to-const-struct with `restrict` rather than a C++ reference -- for
   consumers whose core must remain C with extern "C" linkage. The C++ form
   remains available for parity with the existing one-electron kernels.

Everything else follows the established pattern: memoised backward DAG, global
common-subexpression elimination, topological emission with every node written
exactly once, and the peak-liveness slot allocator sizing a single scratch
array.

Gradients use the same augmented-kernel scheme as §8.3/§8.6: ONE kernel per
(la,lb) emits the value block AND the raised(+1)/lowered(-1) blocks on both
centres, for every moment component, over the SINGLE shared DAG. The consumer
forms

    d/dA_i M(a,b,e) = 2*alpha_A * M(a+1_i, b, e) - a_i * M(a-1_i, b, e)

and analogously on B, with the exponent factors applied per primitive at
contraction. The shift identity is independent of e, so moments cost no new
recurrence -- which is the point.
"""
from __future__ import annotations

from typing import Dict, List, Sequence, Tuple

from onee_dag import MomE, OneE, all_ang_tuples, expand, expand_moment
from onee_dag.recurrence import Term

_CH = "spdfghik"
_XYZ = "xyz"


def pair_name(la: int, lb: int) -> str:
    return f"{_CH[la]}{_CH[lb]}"


# Moment component sets, in the order the consumer expects.
#   dipole    : x, y, z
#   quadrupole: xx, yy, zz, xy, xz, yz   (the Cartesian d ordering)
E_DIPOLE: Tuple[Tuple[int, int, int], ...] = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
E_QUADRUPOLE: Tuple[Tuple[int, int, int], ...] = (
    (2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1))
E_OVERLAP: Tuple[Tuple[int, int, int], ...] = ((0, 0, 0),)


# Cartesian component orderings. RECURSUM's own enumeration is lexicographic
# descending; the "reference" ordering below is the one used by the GFN
# parameterisations and by most integral libraries, where the pure powers come
# first and the mixed ones follow:
#
#   RECURSUM  d: (200) (110) (101) (020) (011) (002)
#   reference d: (200) (020) (002) (110) (101) (011)
#
# The two coincide for s and p, which is precisely why a mismatch survives
# every cheap test and only shows up once d shells appear. Emitting directly in
# the consumer's ordering is better than permuting at the call site: the
# permutation would otherwise sit in the hot path and in every caller.
_REFERENCE_CART = {
    0: [(0, 0, 0)],
    1: [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
    2: [(2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (1, 0, 1), (0, 1, 1)],
    3: [(3, 0, 0), (0, 3, 0), (0, 0, 3), (2, 1, 0), (2, 0, 1), (1, 2, 0),
        (0, 2, 1), (1, 0, 2), (0, 1, 2), (1, 1, 1)],
}


def cart_order(l: int, convention: str = "reference"):
    """Cartesian components of a shell, in the requested ordering."""
    if convention == "recursum":
        return list(all_ang_tuples(l))
    if convention == "reference":
        order = _REFERENCE_CART.get(l)
        if order is None:
            raise ValueError(f"no reference ordering for l={l}")
        got, want = sorted(order), sorted(all_ang_tuples(l))
        assert got == want, f"reference ordering for l={l} is not a permutation"
        return order
    raise ValueError(f"unknown convention {convention!r}")


def moment_set(max_e: int) -> Tuple[Tuple[int, int, int], ...]:
    """All components with |e| <= max_e, overlap first then dipole then
    quadrupole. One kernel serves all three because they share a DAG."""
    out = E_OVERLAP
    if max_e >= 1:
        out = out + E_DIPOLE
    if max_e >= 2:
        out = out + E_QUADRUPOLE
    if max_e > 2:
        raise ValueError("moments beyond quadrupole are not wired up")
    return out


# --------------------------------------------------------------------------
# Fused DAG: moment layer over overlap layer, one memoised graph.
# --------------------------------------------------------------------------

def _expand_any(n) -> List[Term]:
    """Expand either node type. A moment at e == 0 is bridged to the overlap
    value at the same (a,b) with a unit coefficient, which is what fuses the
    two layers into one graph."""
    if isinstance(n, MomE):
        if n.is_base():
            return [Term("one", OneE(a=n.a, b=n.b, m=0), sign=+1)]
        return expand_moment(n)
    return expand(n, "overlap")


def _sort_key(n):
    """Total, deterministic order across both node types. MomE and OneE are not
    mutually comparable, so rank by type first; within a type fall back to the
    node's own ordering. Determinism matters: the emitted source must be
    reproducible byte for byte."""
    if isinstance(n, MomE):
        return (1, n.n_mom, -n.L_total, n.a, n.b, n.e)
    return (0, 0, -(n.L_a + n.L_b), n.a, n.b, n.m)


class FusedDAG:
    def __init__(self) -> None:
        self.nodes: Dict[object, List[Term]] = {}

    def add(self, n) -> None:
        if n in self.nodes:
            return
        terms = _expand_any(n)
        self.nodes[n] = terms
        for t in terms:
            self.add(t.source)

    def topological_order(self) -> List:
        in_deg = {n: 0 for n in self.nodes}
        consumers: Dict[object, List] = {n: [] for n in self.nodes}
        for node, terms in self.nodes.items():
            for t in terms:
                in_deg[node] += 1
                consumers[t.source].append(node)
        ready = sorted([n for n, d in in_deg.items() if d == 0], key=_sort_key)
        out: List = []
        while ready:
            n = ready.pop(0)
            out.append(n)
            for c in consumers[n]:
                in_deg[c] -= 1
                if in_deg[c] == 0:
                    ready.append(c)
                    ready.sort(key=_sort_key)
        assert len(out) == len(self.nodes), "topological sort dropped nodes"
        return out


def build_fused_dag(outputs: Sequence) -> FusedDAG:
    dag = FusedDAG()
    for o in outputs:
        dag.add(o)
    return dag


# --------------------------------------------------------------------------
# Peak-liveness slot allocator (same scheme as the S/T/V emitter).
# --------------------------------------------------------------------------

def liveness_slots(topo, dag: FusedDAG, outputs) -> Tuple[Dict[object, int], int]:
    pos = {n: i for i, n in enumerate(topo)}
    last_use: Dict[object, int] = {}
    for n in topo:
        for t in dag.nodes[n]:
            last_use[t.source] = max(last_use.get(t.source, -1), pos[n])
    die_at: Dict[int, List] = {}
    for v, p in last_use.items():
        die_at.setdefault(p, []).append(v)
    slot: Dict[object, int] = {}
    free: List[int] = []
    nxt = 0
    for i, n in enumerate(topo):
        for v in die_at.get(i - 1, []):
            if v in slot:
                free.append(slot[v])
        if n in outputs or not dag.nodes[n] or n not in last_use:
            continue
        if free:
            slot[n] = free.pop()
        else:
            slot[n] = nxt
            nxt += 1
    return slot, nxt


# --------------------------------------------------------------------------
# Emission
# --------------------------------------------------------------------------

def _coef_to_c(coef: str, lang: str) -> str:
    if coef == "one":
        return "1.0"
    return f"s->{coef}" if lang == "c99" else f"s.{coef}"


def _node_rhs(n, dag: FusedDAG, ref, base_expr: str, lang: str) -> str:
    terms = dag.nodes[n]
    if not terms:
        return base_expr
    pieces = []
    for t in terms:
        c = _coef_to_c(t.coef, lang)
        if t.imul != 1:
            c = f"{t.imul}.0" if c == "1.0" else f"{t.imul}.0 * {c}"
        r = ref(t.source)
        sign = "-" if t.sign < 0 else "+"
        pieces.append(f"{sign} {r}" if c == "1.0" else f"{sign} {c} * {r}")
    expr = " ".join(pieces)
    return expr[2:] if expr.startswith("+ ") else expr


def _outputs_for(la: int, lb: int, es, convention: str = "reference"
                 ) -> Tuple[List[MomE], Dict]:
    """Output nodes and their section layout: moment component major, then the
    Cartesian (a,b) block, so a consumer can take a contiguous block per
    component."""
    outs: List[MomE] = []
    sections: Dict = {}
    for e in es:
        start = len(outs)
        for a in cart_order(la, convention):
            for b in cart_order(lb, convention):
                outs.append(MomE(a=a, b=b, e=e))
        sections[e] = (start, len(outs) - start)
    return outs, sections


def emit_value_kernel(la: int, lb: int, max_e: int, prefix: str = "scdt",
                      lang: str = "c99") -> Tuple[str, dict]:
    """One kernel emitting S, D and Q blocks for a (la,lb) class from one DAG."""
    es = moment_set(max_e)
    outs, sections = _outputs_for(la, lb, es)
    outset = set(outs)
    oidx = {o: i for i, o in enumerate(outs)}
    dag = build_fused_dag(outs)
    topo = dag.topological_order()
    slot, peak = liveness_slots(topo, dag, outset)
    name = pair_name(la, lb)

    def ref(src):
        if src in oidx:
            return f"out[{oidx[src]}]"
        if not dag.nodes[src]:
            return "S00"
        return f"sc[{slot[src]}]"

    scal = "const scdt_onee_scalars *restrict s" if lang == "c99" \
        else "const OneEScalars& __restrict__ s"
    outp = "double *restrict out" if lang == "c99" else "double* __restrict__ out"
    layout = ", ".join(f"e={''.join(map(str, e))}:[{st}:{st+ct}]"
                       for e, (st, ct) in sections.items())

    L = [f"/* Multipole ({name}), |e| <= {max_e}: {len(outs)} outputs,"
         f" {len(dag.nodes)} fused DAG nodes, peak liveness {peak}.",
         f" * Sections: {layout}",
         " * Overlap sub-graph is shared by every moment component (global CSE).",
         " */",
         f"void {prefix}_moment_{name}({scal}, double S00, {outp})",
         "{"]
    if peak > 0:
        L.append(f"    double sc[{peak}];")
    for n in topo:
        if not dag.nodes[n] and n not in outset:
            continue
        rhs = _node_rhs(n, dag, ref, "S00", lang)
        if n in outset:
            L.append(f"    out[{oidx[n]}] = {rhs};")
        elif n in slot:
            L.append(f"    sc[{slot[n]}] = {rhs};")
    L.append("}")
    stats = {"outputs": len(outs), "nodes": len(dag.nodes), "peak": peak,
             "sections": sections}
    return "\n".join(L), stats


def grad_layout(la: int, lb: int, es, convention: str = "reference"):
    """value ++ upA ++ dnA ++ upB ++ dnB, each over all moment components.

    The shift identity is independent of the moment index, so the same five
    sections serve S, D and Q at once."""
    outs: List[MomE] = []
    sec: Dict[str, Tuple[int, int, Tuple[int, int]]] = {}

    def add(tag, a, b):
        start = len(outs)
        for e in es:
            for aa in cart_order(a, convention):
                for bb in cart_order(b, convention):
                    outs.append(MomE(a=aa, b=bb, e=e))
        sec[tag] = (start, len(outs) - start, (a, b))

    add("val", la, lb)
    add("upA", la + 1, lb)
    if la >= 1:
        add("dnA", la - 1, lb)
    add("upB", la, lb + 1)
    if lb >= 1:
        add("dnB", la, lb - 1)
    return outs, sec


def emit_grad_kernel(la: int, lb: int, max_e: int, prefix: str = "scdt",
                     lang: str = "c99") -> Tuple[str, dict]:
    """Augmented value+gradient kernel: one shared DAG, all shifts, all moments.

    The value sub-DAG is CSE-shared with every shift output rather than
    recomputed per shift, which is the whole reason this is one kernel."""
    es = moment_set(max_e)
    outs, sec = grad_layout(la, lb, es)
    outset = set(outs)
    oidx = {o: i for i, o in enumerate(outs)}
    dag = build_fused_dag(outs)
    topo = dag.topological_order()
    slot, peak = liveness_slots(topo, dag, outset)
    name = pair_name(la, lb)

    def ref(src):
        if src in oidx:
            return f"out[{oidx[src]}]"
        if not dag.nodes[src]:
            return "S00"
        return f"sc[{slot[src]}]"

    scal = "const scdt_onee_scalars *restrict s" if lang == "c99" \
        else "const OneEScalars& __restrict__ s"
    outp = "double *restrict out" if lang == "c99" else "double* __restrict__ out"
    hdr = ", ".join(f"{k}:[{st}:{st+ct}] ({cl[0]},{cl[1]})"
                    for k, (st, ct, cl) in sec.items())

    L = [f"/* Multipole GRADIENT ({name}), |e| <= {max_e}: {len(outs)} outputs,"
         f" {len(dag.nodes)} fused DAG nodes, peak liveness {peak}.",
         f" * Sections: {hdr}",
         " * Consumer forms  d/dA_i = 2*alpha_A * upA[+1_i] - a_i * dnA[-1_i],",
         " * and analogously on B. One shared DAG; no per-shift recompute.",
         " */",
         f"void {prefix}_moment_grad_{name}({scal}, double S00, {outp})",
         "{"]
    if peak > 0:
        L.append(f"    double sc[{peak}];")
    for n in topo:
        if not dag.nodes[n] and n not in outset:
            continue
        rhs = _node_rhs(n, dag, ref, "S00", lang)
        if n in outset:
            L.append(f"    out[{oidx[n]}] = {rhs};")
        elif n in slot:
            L.append(f"    sc[{slot[n]}] = {rhs};")
    L.append("}")
    stats = {"outputs": len(outs), "nodes": len(dag.nodes), "peak": peak,
             "sections": sec}
    return "\n".join(L), stats


SCALARS_C99 = """/* AUTO-GENERATED by recursum_moment_emit.py -- do not edit.
 * Per-primitive-pair scalars for the one-electron multipole kernels. */
#ifndef SCDT_ONEE_SCALARS_H
#define SCDT_ONEE_SCALARS_H

typedef struct {
    double PAx, PAy, PAz;   /* P - A, used by the overlap VRR   */
    double PCx, PCy, PCz;   /* P - C, C the moment origin       */
    double ABx, ABy, ABz;   /* A - B, used by the HRR A->B      */
    double inv_2zeta;       /* 1 / (2 zeta)                     */
} scdt_onee_scalars;
/* P - B is deliberately absent: the overlap VRR builds on A and the HRR
 * transfers A->B, so no rule in this DAG references it. */

#endif
"""


# --------------------------------------------------------------------------
# Project emission
# --------------------------------------------------------------------------

def generate_project(pairs, max_e: int = 2, prefix: str = "scdt",
                     lang: str = "c99", outdir: str = ".",
                     basename: str = "moment_kernels"):
    """Emit the kernel source, its declarations and the scalars header.

    `pairs` are canonical (la >= lb) classes; the consumer handles the
    transposed case by swapping arguments and transposing the block, which
    halves the number of kernels."""
    import os

    ext = "c" if lang == "c99" else "cpp"
    src = [f'/* AUTO-GENERATED by recursum_moment_emit.py -- do not edit.',
           ' *',
           ' * Cartesian multipole integrals (overlap, dipole, quadrupole) and their',
           ' * geometric gradients, emitted from one fused Obara-Saika DAG per class',
           ' * with global common-subexpression elimination across the moment layers.',
           ' */',
           f'#include "{basename}.h"',
           '']
    decl = ['/* AUTO-GENERATED by recursum_moment_emit.py -- do not edit. */',
            f'#ifndef {basename.upper()}_H',
            f'#define {basename.upper()}_H',
            '',
            '#include "scdt_onee_scalars.h"',
            '']
    stats = {}
    for la, lb in pairs:
        v, sv = emit_value_kernel(la, lb, max_e, prefix, lang)
        g, sg = emit_grad_kernel(la, lb, max_e, prefix, lang)
        src.append(v)
        src.append("")
        src.append(g)
        src.append("")
        nm = pair_name(la, lb)
        scal = "const scdt_onee_scalars *restrict s" if lang == "c99" \
            else "const OneEScalars& __restrict__ s"
        outp = "double *restrict out" if lang == "c99" else "double* __restrict__ out"
        decl.append(f"void {prefix}_moment_{nm}({scal}, double S00, {outp});")
        decl.append(f"void {prefix}_moment_grad_{nm}({scal}, double S00, {outp});")
        stats[nm] = {"value": sv, "grad": sg}

    decl.append("")
    decl.append("/* Output counts per class, so a caller can size its buffers. */")
    for la, lb in pairs:
        nm = pair_name(la, lb)
        decl.append(f"#define {prefix.upper()}_MOMENT_NOUT_{nm.upper()} "
                    f"{stats[nm]['value']['outputs']}")
        decl.append(f"#define {prefix.upper()}_MOMENT_GRAD_NOUT_{nm.upper()} "
                    f"{stats[nm]['grad']['outputs']}")
    decl.append("")
    decl.append("#endif")

    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, f"{basename}.{ext}"), "w") as f:
        f.write("\n".join(src))
    with open(os.path.join(outdir, f"{basename}.h"), "w") as f:
        f.write("\n".join(decl) + "\n")
    with open(os.path.join(outdir, "scdt_onee_scalars.h"), "w") as f:
        f.write(SCALARS_C99 if lang == "c99" else SCALARS_C99)
    return stats


GFN2_PAIRS = [(0, 0), (1, 0), (1, 1), (2, 0), (2, 1), (2, 2)]


if __name__ == "__main__":
    import sys
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    st = generate_project(GFN2_PAIRS, max_e=2, outdir=outdir)
    tot_v = sum(s["value"]["nodes"] for s in st.values())
    tot_g = sum(s["grad"]["nodes"] for s in st.values())
    print(f"emitted {len(st)} classes into {outdir}/")
    for nm, s in st.items():
        print("  %-3s value %4d out / %5d nodes / peak %3d   "
              "grad %5d out / %6d nodes / peak %4d"
              % (nm, s["value"]["outputs"], s["value"]["nodes"], s["value"]["peak"],
                 s["grad"]["outputs"], s["grad"]["nodes"], s["grad"]["peak"]))
    print(f"total DAG nodes: value {tot_v}, gradient {tot_g}")
