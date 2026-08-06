"""RECURSUM DAG-topological layered emitter for HGP-OS ERIs, with the osx
peak-liveness slot-reuse optimization ported in.

This is the new RECURSUM generator (prototype of recursum/codegen/
dag_layered_generator.py). Unlike the tower-based LayeredCppGenerator, it:

  1. Expands the DSL recurrence into a memoized integral DAG (4 coupled
     angular-momentum groups + Boys aux index) -- the osx dag.py strategy.
  2. Emits every DAG node EXACTLY ONCE in topological order (true CSE:
     no exponential recomputation of shared intermediates).
  3. LIVENESS OPTIMIZATION (ported from osx emit.py:1138-1165): a linear-scan
     slot allocator assigns each intermediate a scratch slot and REUSES the
     slot once the value is dead. The kernel's stack frame is bounded by the
     DAG's PEAK LIVENESS, not its node count -- the key to matching libint2's
     register pressure at high angular momentum.

Fairness with libint2's low-level API: the kernel receives the same
precomputed prerequisites libint2 does -- the per-quartet recurrence scalars
(PA, WP, inv_2z*, ...) and the prescaled base integrals KF[m] = K * F_m(T).
It runs only the VRR/HRR contraction. So we time recurrence-vs-recurrence.
"""
from __future__ import annotations

import re
from typing import Dict, List, Tuple

from eri_dag.integral import Integral, output_set, all_ang_tuples
from eri_dag.dag import build_dag


# ---- coefficient string -> C expression (port of osx emit.py _scalar_for) --
# Base scalars live in a struct `s`; composite names expand to products.
def coef_to_c(coef: str) -> str:
    sfx = None
    if coef.endswith("_xq"):
        sfx = "s.frac_q_over_pq"; coef = coef[:-3]
    elif coef.endswith("_xp"):
        sfx = "s.frac_p_over_pq"; coef = coef[:-3]

    if coef.startswith("ai") and "_inv_2zp" in coef:
        n = int(coef[2:coef.index("_inv_2zp")]); base = f"{n}.0 * s.inv_2zp"
    elif coef.startswith("ci") and "_inv_2zq" in coef:
        n = int(coef[2:coef.index("_inv_2zq")]); base = f"{n}.0 * s.inv_2zq"
    elif coef.startswith("ai") and "_inv_2zpq" in coef:
        n = int(coef[2:coef.index("_inv_2zpq")]); base = f"{n}.0 * s.inv_2zpq"
    elif coef.startswith("ci") and "_inv_2zpq" in coef:
        n = int(coef[2:coef.index("_inv_2zpq")]); base = f"{n}.0 * s.inv_2zpq"
    elif coef == "one":
        base = "1.0"
    else:
        base = f"s.{coef}"  # PAx, WPx, QCx, WQx, ABx, CDx

    if sfx is not None:
        base = f"({base}) * {sfx}"
    return base


# All base scalar field names the struct must carry.
SCALAR_FIELDS = [
    "PAx", "PAy", "PAz", "QCx", "QCy", "QCz",
    "WPx", "WPy", "WPz", "WQx", "WQy", "WQz",
    "ABx", "ABy", "ABz", "CDx", "CDy", "CDz",
    "inv_2zp", "inv_2zq", "inv_2zpq", "frac_q_over_pq", "frac_p_over_pq",
]


def class_name(la, lb, lc, ld) -> str:
    ch = "spdfghik"
    return f"{ch[la]}{ch[lb]}{ch[lc]}{ch[ld]}"


# ---------------------------------------------------------------------------
def liveness_slots(topo: List[Integral], dag,
                   outputs: set, reuse: bool = True) -> Tuple[Dict[Integral, int], int]:
    """Linear-scan slot allocation (port of osx emit.py:1138-1165).

    Only INTERMEDIATE nodes get a scratch slot; outputs go straight to out[]
    and base integrals are always read from kf[]. With reuse=True a slot is
    freed the moment its value dies and recycled (peak-liveness minimization).
    With reuse=False every intermediate keeps its own slot (naive baseline) --
    used to isolate the liveness optimization's effect under perf.
    Returns (slot_of_node, peak_slot_count)."""
    pos = {n: i for i, n in enumerate(topo)}

    # last_use[v] = last topo position at which v is consumed as a source.
    last_use: Dict[Integral, int] = {}
    for n in topo:
        for t in dag.nodes[n]:
            last_use[t.source] = max(last_use.get(t.source, -1), pos[n])

    # die_at[p] = list of values whose last consumer is at position p.
    die_at: Dict[int, List[Integral]] = {}
    for v, p in last_use.items():
        die_at.setdefault(p, []).append(v)

    slot: Dict[Integral, int] = {}
    free: List[int] = []
    next_slot = 0
    for i, n in enumerate(topo):
        if reuse:
            # free slots of values that died at the PREVIOUS position
            for v in die_at.get(i - 1, []):
                if v in slot:
                    free.append(slot[v])
        if n in outputs:
            continue  # outputs go to out[], not a scratch slot
        if not dag.nodes[n]:
            continue  # base integral: always referenced via kf[], no slot
        # a node that is never consumed and isn't an output = dead code; skip slot
        if n not in last_use:
            continue
        if reuse and free:
            slot[n] = free.pop()
        else:
            slot[n] = next_slot; next_slot += 1
    return slot, next_slot


# ---- SymPy algebraic-simplification pass (spec §5.3; open-source, NOT Maple) --
# Off by default. When enabled, each node's RHS is rebuilt as a SymPy expression,
# simplified (shared-multiplier factoring + constant folding via factor_terms),
# and PER-NODE NUMERICALLY VERIFIED equivalent before being emitted; any node
# that fails verification falls back to the plain expression. This is the
# fidelity gate the spec requires ("random numeric check per node").
SYMPY_PASS = False
_SYMPY_STATS = {"nodes": 0, "simplified": 0, "rejected": 0}


def coef_to_sympy(coef, sym):
    """SymPy analog of coef_to_c: build the coefficient as a SymPy expression
    over scalar-field symbols (from the `sym` cache). Mirrors coef_to_c exactly."""
    import sympy
    def S(name):
        if name not in sym:
            sym[name] = sympy.Symbol(name)
        return sym[name]
    sfx = None
    if coef.endswith("_xq"):
        sfx = S("frac_q_over_pq"); coef = coef[:-3]
    elif coef.endswith("_xp"):
        sfx = S("frac_p_over_pq"); coef = coef[:-3]
    if coef.startswith("ai") and "_inv_2zp" in coef:
        n = int(coef[2:coef.index("_inv_2zp")]); base = sympy.Integer(n) * S("inv_2zp")
    elif coef.startswith("ci") and "_inv_2zq" in coef:
        n = int(coef[2:coef.index("_inv_2zq")]); base = sympy.Integer(n) * S("inv_2zq")
    elif coef.startswith("ai") and "_inv_2zpq" in coef:
        n = int(coef[2:coef.index("_inv_2zpq")]); base = sympy.Integer(n) * S("inv_2zpq")
    elif coef.startswith("ci") and "_inv_2zpq" in coef:
        n = int(coef[2:coef.index("_inv_2zpq")]); base = sympy.Integer(n) * S("inv_2zpq")
    elif coef == "one":
        base = sympy.Integer(1)
    else:
        base = S(coef)
    return base * sfx if sfx is not None else base


def _print_c(expr, srcmap):
    """Emit a restricted SymPy expr (Add/Mul/Symbol/Number/Pow) as C. Scalar
    symbols -> s.<name>; source symbols -> their caller ref (srcmap)."""
    import sympy
    if expr.is_Add:
        return "(" + " + ".join(_print_c(a, srcmap) for a in expr.as_ordered_terms()) + ")"
    if expr.is_Mul:
        return "(" + " * ".join(_print_c(a, srcmap) for a in expr.as_ordered_factors()) + ")"
    if expr.is_Pow:
        b, e = expr.as_base_exp()
        assert e.is_Integer and int(e) >= 1 and int(e) <= 4
        return "(" + " * ".join([_print_c(b, srcmap)] * int(e)) + ")"
    if expr.is_Symbol:
        nm = expr.name
        return srcmap[nm] if nm in srcmap else f"s.{nm}"
    if expr.is_Integer:
        return f"{int(expr)}.0"
    if expr.is_Rational:
        return f"({int(expr.p)}.0 / {int(expr.q)}.0)"
    if expr.is_Float or expr.is_Number:
        return repr(float(expr))
    raise ValueError(f"unprintable node in SymPy pass: {expr!r}")


def _sympy_node_rhs(terms, reffn):
    """Simplified C RHS for a node's term list, or None if the pass rejects it."""
    import sympy
    sym = {}
    srcmap = {}       # source-symbol name -> caller ref C string
    orig = sympy.Integer(0)
    for i, t in enumerate(terms):
        sname = f"__S{i}"
        srcmap[sname] = reffn(t.source)
        ssym = sympy.Symbol(sname)
        sym[sname] = ssym
        orig += sympy.Integer(t.sign) * coef_to_sympy(t.coef, sym) * ssym
    simplified = sympy.factor_terms(orig)
    # per-node numeric fidelity gate: random-substitution equivalence check.
    syms = sorted(sym.values(), key=lambda x: x.name)
    import random
    rng = random.Random(0xC0FFEE)
    for _ in range(6):
        subs = {s: sympy.Float(rng.uniform(-1.7, 1.9), 17) for s in syms}
        a = float(orig.xreplace(subs)); b = float(simplified.xreplace(subs))
        if abs(a - b) > 1e-12 * (abs(a) + 1e-300):
            return None  # not value-preserving under our printer's semantics -> reject
    try:
        return _print_c(simplified, srcmap)
    except (ValueError, AssertionError):
        return None


def node_rhs(n: Integral, dag, reffn, kf="kf") -> str:
    """C expression for node n's value; reffn(src) names each source."""
    terms = dag.nodes[n]
    if not terms:  # base integral (0000)^(m)
        return f"{kf}[{n.m}]"
    if SYMPY_PASS:
        _SYMPY_STATS["nodes"] += 1
        c = _sympy_node_rhs(terms, reffn)
        if c is not None:
            _SYMPY_STATS["simplified"] += 1
            return c
        _SYMPY_STATS["rejected"] += 1
    pieces = []
    for t in terms:
        c = coef_to_c(t.coef)
        r = reffn(t.source)
        sign = "-" if t.sign < 0 else "+"
        if c == "1.0":
            pieces.append(f"{sign} {r}")
        else:
            pieces.append(f"{sign} {c} * {r}")
    expr = " ".join(pieces)
    if expr.startswith("+ "):
        expr = expr[2:]
    return expr


def _frontier_slots(chunks, chunk_of, dag, topo):
    """Assign frontier-scratch slots to boundary values (interior to one chunk,
    consumed by a later chunk), reusing a slot once its last-consuming chunk has
    run. Port of osx emit.py:1138-1165 at chunk granularity. Returns
    (slot_of_boundary_value, peak_slot_count)."""
    pos = {n: i for i, n in enumerate(topo)}
    consumers = {n: [] for n in dag.nodes}
    for n, terms in dag.nodes.items():
        for t in terms:
            consumers[t.source].append(n)
    # boundary value -> (producing chunk, last consuming chunk)
    prod_chunk, last_use_chunk = {}, {}
    for n in dag.nodes:
        pc = chunk_of[n]
        for c in consumers[n]:
            cc = chunk_of[c]
            if cc > pc:  # crosses a chunk boundary => n is a boundary value
                prod_chunk[n] = pc
                last_use_chunk[n] = max(last_use_chunk.get(n, pc), cc)
    die_at = {}
    for v, cc in last_use_chunk.items():
        die_at.setdefault(cc, []).append(v)
    slot, free, nxt = {}, [], 0
    n_chunks = len(chunks)
    for ci in range(n_chunks):
        for v in die_at.get(ci - 1, []):      # free values whose last consumer was ci-1
            if v in slot:
                free.append(slot[v])
        for v in [b for b in prod_chunk if prod_chunk[b] == ci]:
            if v in slot:
                continue
            slot[v] = free.pop() if free else nxt
            if not free or slot[v] == nxt:
                nxt = max(nxt, slot[v] + 1)
    return slot, nxt


def emit_kernel_chunked(la, lb, lc, ld, suffix="", interior_budget=1200):
    """Emit a class as many small chunk functions + a driver, so a huge kernel
    (e.g. ffff, 77k nodes) compiles fast (bounded per-function size) instead of
    as one pathological function. Ported from the osx DAG chunker.

    Each chunk is a static noinline function: interior nodes are `const double`
    locals (compiler does register allocation), cross-chunk boundary values pass
    through a frontier scratch array fr[] (liveness-reused), final targets go to
    out[]. Returns (source, n_out, n_nodes, fr_peak, n_chunks)."""
    from eri_dag.chunker import partition
    outputs_list = output_set(la, lb, lc, ld)
    outputs = set(outputs_list)
    out_index = {o: i for i, o in enumerate(outputs_list)}
    dag = build_dag(outputs_list)
    topo = dag.topological_order()
    chunks = partition(dag, topo, interior_budget=interior_budget,
                       frontier_budget=10**9, import_budget=10**9,
                       working_set_budget=10**9)
    chunk_of = {}
    for ci, ch in enumerate(chunks):
        for n in ch.interior:
            chunk_of[n] = ci
    fr_slot, fr_peak = _frontier_slots(chunks, chunk_of, dag, topo)

    name = class_name(la, lb, lc, ld) + suffix
    pos = {n: i for i, n in enumerate(topo)}
    L = [f"// ERI class ({name}) [chunked]: {len(outputs_list)} integrals, "
         f"{len(dag.nodes)} nodes, {len(chunks)} chunks, frontier {fr_peak} slots",
         '#include "recursum_eri_scalars.h"']

    def reffn(src, local):  # local = set of nodes emitted as locals in this chunk
        if src in out_index:  return f"out[{out_index[src]}]"   # (never happens)
        if not dag.nodes[src]: return f"kf[{src.m}]"
        if src in local:      return f"v{pos[src]}"
        return f"fr[{fr_slot[src]}]"   # boundary value from an earlier chunk

    # chunk functions
    for ci, ch in enumerate(chunks):
        L.append(f"static __attribute__((noinline)) void {name}_c{ci}("
                 f"const ScalarPack& __restrict__ s, const double* __restrict__ kf, "
                 f"double* __restrict__ fr, double* __restrict__ out) {{")
        local = set()
        for n in ch.interior:
            if not dag.nodes[n] and n not in outputs:
                continue  # base -> kf[m], never emitted (unless it IS an output)
            rhs = node_rhs(n, dag, lambda src: reffn(src, local))
            if n in outputs:
                L.append(f"    out[{out_index[n]}] = {rhs};")
            else:
                L.append(f"    const double v{pos[n]} = {rhs};")
                local.add(n)
                if n in fr_slot:                      # boundary value -> persist
                    L.append(f"    fr[{fr_slot[n]}] = v{pos[n]};")
        L.append("}")

    # driver — Uniform scratch (§10.6): the frontier buffer is the caller-owned
    # `sc` parameter (arena Frame), sized by the emitted `_nscratch` == fr_peak.
    L.append(f"__attribute__((noinline)) void recursum_eri_{name}("
             "const ScalarPack& __restrict__ s, const double* __restrict__ kf, "
             "double* __restrict__ out, double* __restrict__ sc) {")
    for ci in range(len(chunks)):
        L.append(f"    {name}_c{ci}(s, kf, sc, out);")
    L.append("}")
    return "\n".join(L), len(outputs_list), len(dag.nodes), fr_peak, len(chunks)


def emit_kernel(la, lb, lc, ld, reuse=True, suffix="", style="array"):
    """Emit the C kernel for class (la lb|lc ld).

    style="array": intermediates in a scratch array sc[] with peak-liveness
        slot reuse (reuse=True) or one-slot-per-node (reuse=False, naive).
    style="ssa": intermediates as named `const double` locals, one per DAG
        node (global CSE, each node once), letting the C compiler do liveness/
        register allocation -- the LayeredCodegen-style codegen. This is
        RECURSUM's fastest CPU variant (gcc's stack coloring is provably
        optimal here, so explicit slot reuse buys nothing).

    Returns (source, n_out, n_nodes, peak_slots)."""
    outputs_list = output_set(la, lb, lc, ld)
    outputs = set(outputs_list)
    out_index = {o: i for i, o in enumerate(outputs_list)}
    dag = build_dag(outputs_list)
    topo = dag.topological_order()
    slot, peak = liveness_slots(topo, dag, outputs, reuse=reuse)

    name = class_name(la, lb, lc, ld) + suffix
    # Uniform scratch handoff (SHRIKE spec §10.6): scratch is a caller-owned
    # __restrict pointer `sc` (from the arena Frame), never a stack-local array.
    L = [f"// ERI class ({name}) [{style}]: {len(outputs_list)} Cartesian integrals, "
         f"{len(dag.nodes)} DAG nodes, peak liveness {peak} slots",
         f"__attribute__((noinline)) void recursum_eri_{name}(",
         "    const ScalarPack& __restrict__ s,",
         "    const double* __restrict__ kf,   // kf[m] = K * F_m(T)",
         "    double* __restrict__ out,",
         "    double* __restrict__ sc) {   // caller-owned scratch, >= _nscratch"]

    if style == "ssa":
        # name every intermediate node by topo index; compiler handles liveness
        L.append("    (void)sc;  // ssa style uses named locals, not caller scratch")
        vid = {}
        def refssa(src):
            if src in out_index: return f"out[{out_index[src]}]"
            if not dag.nodes[src]: return f"kf[{src.m}]"
            return vid[src]
        for i, n in enumerate(topo):
            if not dag.nodes[n] and n not in outputs:
                continue  # base -> referenced as kf[m] (unless it IS an output)
            if n in outputs:
                L.append(f"    out[{out_index[n]}] = {node_rhs(n, dag, refssa)};")
            elif n in slot:
                vid[n] = f"v{i}"
                L.append(f"    const double v{i} = {node_rhs(n, dag, refssa)};")
        L.append("}")
        return "\n".join(L), len(outputs_list), len(dag.nodes), peak

    # style == "array"
    def refarr(src):
        if src in out_index: return f"out[{out_index[src]}]"
        if not dag.nodes[src]: return f"kf[{src.m}]"
        return f"sc[{slot[src]}]"
    # `sc` is now the caller-owned parameter (Uniform, §10.6) — no local array.
    if peak == 0:
        L.append("    (void)sc;")
    for n in topo:
        if not dag.nodes[n] and n not in outputs:
            continue  # base -> kf[m] (unless it IS an output, e.g. ssss)
        rhs = node_rhs(n, dag, refarr)
        if n in outputs:
            L.append(f"    out[{out_index[n]}] = {rhs};")
        elif n in slot:
            L.append(f"    sc[{slot[n]}] = {rhs};")
    L.append("}")
    return "\n".join(L), len(outputs_list), len(dag.nodes), peak


def grad_output_layout(la, lb, lc, ld):
    """AUGMENTED value+gradient output layout (spec §8.3). Returns
    (outputs_flat, sections) where outputs_flat is the value block followed by the
    raised (+1) and lowered (-1) blocks for the three explicit centres A,B,C
    (the 4th centre D is recovered by translational invariance §8.2). sections
    maps name -> (start, count, class) so the caller/digestion knows each block's
    offset. Centres with l==0 have no lowered block. Building ONE DAG over the
    whole list lets global CSE share the value sub-DAG across every derivative
    output -- the codegen win a per-derivative kernel throws away."""
    Ls = [la, lb, lc, ld]
    outs = []
    sections = {}
    def add(name, cls):
        o = output_set(*cls)
        sections[name] = (len(outs), len(o), cls)
        outs.extend(o)
    add("val", (la, lb, lc, ld))
    for ci, cn in enumerate(["A", "B", "C"]):
        up = list(Ls); up[ci] += 1
        add("up" + cn, tuple(up))
        if Ls[ci] >= 1:
            dn = list(Ls); dn[ci] -= 1
            add("dn" + cn, tuple(dn))
    return outs, sections


def emit_grad_kernel(la, lb, lc, ld, suffix=""):
    """Emit recursum_grad_<cls>: ONE augmented value+gradient kernel (spec §8.3).

    Produces the value block AND the raised/lowered shift blocks for centres
    A,B,C in a single straight-line body over ONE DAG -- global CSE shares the
    value sub-DAG with every derivative output (only the top angular-momentum
    layer of the VRR is new; marginal cost << the naive per-derivative 9x). The
    C++ side forms d/dX_c = 2*alpha * up(a+1_c) - a_c * dn(a-1_c) from the
    sections (the 2*alpha is applied at contraction, spec §8.4). ssa style
    (named locals; gcc does liveness). Returns (source, sections, nout, nnodes)."""
    outputs_list, sections = grad_output_layout(la, lb, lc, ld)
    outputs = set(outputs_list)
    out_index = {o: i for i, o in enumerate(outputs_list)}
    dag = build_dag(outputs_list)
    topo = dag.topological_order()
    slot, _peak = liveness_slots(topo, dag, outputs, reuse=True)
    name = class_name(la, lb, lc, ld) + suffix
    hdr = ", ".join(f"{k}:out[{s}:{s+c}]" for k, (s, c, _cls) in sections.items())
    L = [f"// AUGMENTED value+gradient kernel ({name}) [§8.3]: {len(outputs_list)} outputs, "
         f"{len(dag.nodes)} DAG nodes",
         f"//   sections: {hdr}",
         f"__attribute__((noinline)) void recursum_grad_{name}(",
         "    const ScalarPack& __restrict__ s,",
         "    const double* __restrict__ kf,   // kf[m] = K * F_m(T)",
         "    double* __restrict__ out) {",
         "    (void)s;"]
    vid = {}
    def refssa(src):
        if src in out_index: return f"out[{out_index[src]}]"
        if not dag.nodes[src]: return f"kf[{src.m}]"
        return vid[src]
    for i, n in enumerate(topo):
        if not dag.nodes[n] and n not in outputs:
            continue
        if n in outputs:
            L.append(f"    out[{out_index[n]}] = {node_rhs(n, dag, refssa)};")
        elif n in slot:
            vid[n] = f"v{i}"
            L.append(f"    const double v{i} = {node_rhs(n, dag, refssa)};")
    L.append("}")
    return "\n".join(L), sections, len(outputs_list), len(dag.nodes)

def emit_gradctr_kernel(la, lb, lc, ld, suffix=""):
    """Emit recursum_gradctr_<cls>: the CONTRACTED augmented value+gradient kernel
    (spec §8.3 + §8.4). FUSED form: the in-kernel primitive loop computes the
    augmented DAG outputs (value + shift sections) ONCE per primitive from the
    shared straight-line body (VRR sub-DAG CSE preserved), then accumulates each
    section with its contraction weight — the RAISED sections get an extra factor
    `2*alpha_centre` (the exponent-dependent raise weight, §8.4), the value and
    LOWERED sections get weight 1 (the contraction coeff is already in kf). The
    2*alpha per differentiated centre per primitive is passed in `w2` (cd x 3:
    A,B,C). deriv_X_c = up(a+1_c) - a_c*dn(a-1_c) is then formed by the C++ caller
    from the sections. Verified bit-exact vs finite-difference of the contracted
    value. Returns (source, sections, nout, nnodes)."""
    outputs_list, sections = grad_output_layout(la, lb, lc, ld)
    outputs = set(outputs_list)
    out_index = {o: i for i, o in enumerate(outputs_list)}
    nout = len(outputs_list)
    dag = build_dag(outputs_list)
    topo = dag.topological_order()
    slot, peak = liveness_slots(topo, dag, outputs, reuse=True)
    name = class_name(la, lb, lc, ld) + suffix
    # section -> weight source: raised sections use w2[.][ci]; others weight 1.
    up_w = {"upA": 0, "upB": 1, "upC": 2}
    hdr = ", ".join(f"{k}:out[{s}:{s+c}]" for k, (s, c, _cls) in sections.items())
    L = [f"// CONTRACTED augmented value+gradient kernel ({name}) [§8.3/§8.4]: "
         f"{nout} outputs, {len(dag.nodes)} DAG nodes",
         f"//   sections: {hdr}   (raised sections weighted by 2*alpha, w2[c*3+{{0:A,1:B,2:C}}])",
         f"__attribute__((noinline)) void recursum_gradctr_{name}(",
         "    const ScalarPack* __restrict__ sp,",
         "    const double* __restrict__ kfa,",
         "    const double* __restrict__ w2,   // cd x 3: 2*alpha for centres A,B,C",
         "    int cd, int mstride,",
         "    double* __restrict__ out,",
         "    double* __restrict__ sc) {",
         f"    for(int k=0;k<{nout};++k) out[k]=0.0;",
         "    for(int c=0;c<cd;++c) {",
         "        const ScalarPack& s = sp[c];",
         "        const double* kf = kfa + (long)c*mstride;",
         f"        double o[{nout}];"]
    def fref(src):
        if src in out_index:   return f"o[{out_index[src]}]"
        if not dag.nodes[src]: return f"kf[{src.m}]"
        return f"sc[{slot[src]}]"
    for n in topo:
        if not dag.nodes[n] and n not in out_index:
            continue
        rhs = node_rhs(n, dag, fref)
        if n in out_index:  L.append(f"        o[{out_index[n]}] = {rhs};")
        elif n in slot:     L.append(f"        sc[{slot[n]}] = {rhs};")
    # weighted accumulation per section
    L.append("        const double wA=w2[(long)c*3+0], wB=w2[(long)c*3+1], wC=w2[(long)c*3+2];")
    for k, (st, cnt, _cls) in sections.items():
        if k in up_w:
            w = {"upA": "wA", "upB": "wB", "upC": "wC"}[k]
            L.append(f"        for(int k={st};k<{st+cnt};++k) out[k]+={w}*o[k];")
        else:
            L.append(f"        for(int k={st};k<{st+cnt};++k) out[k]+=o[k];")
    L += ["    }", "}"]
    return "\n".join(L), sections, nout, len(dag.nodes)

def grad_has_split(la, lb, lc, ld) -> bool:
    """The augmented gradient DAG benefits from the split iff some section carries
    an HRR stage. Since upB raises l_b (and upC raises l_d), an HRR stage exists
    for essentially every non-(ss|ss) class; gate on real HRR node count instead.
    SHRIKE_GRAD_NOSPLIT=1 forces the fused form for A/B split-vs-fused benchmarking."""
    import os as _os
    if _os.environ.get("SHRIKE_GRAD_NOSPLIT") == "1":
        return False
    outs, _ = grad_output_layout(la, lb, lc, ld)
    dag = build_dag(outs)
    _v, hrr, _b, _i = split_sets(dag, set(outs))
    return len(hrr) >= 8   # tiny-HRR classes stay fused (call/setup not worth it)


def emit_gradctr_split_kernel(la, lb, lc, ld, suffix=""):
    """SPLIT augmented value+gradient kernel (spec §8.3/§8.4 + the §5.1 split applied
    to the gradient DAG — task #13). Instead of recomputing the WHOLE augmented DAG
    per primitive (the fused recursum_gradctr), it runs the VRR (a0|c0) tower per
    primitive, contracts the boundary into FOUR weight-class accumulators, then does
    HRR ONCE per weight class:
      weight 1  : value + all lowered sections (dnA/dnB/dnC)         → e0f0c_1
      weight 2αA: upA section                                        → e0f0c_A
      weight 2αB: upB section                                        → e0f0c_B
      weight 2αC: upC section                                        → e0f0c_C
    HRR is linear and geometry-only, so Σ_p w_p·HRR(bnd_p)=HRR(Σ_p w_p·bnd_p); the
    per-primitive raise weight 2α_centre is therefore applied at the boundary
    contraction, exactly matching the fused kernel's per-primitive weighting.
    Same signature as recursum_gradctr. Returns (source, sections, nout, nnodes)."""
    outputs_list, sections = grad_output_layout(la, lb, lc, ld)
    out_index = {o: i for i, o in enumerate(outputs_list)}
    nout = len(outputs_list)
    dag = build_dag(outputs_list)
    topo = dag.topological_order()
    vrr_nodes, hrr_nodes, boundary_list, bidx = split_sets(dag, set(outputs_list))
    nb = len(boundary_list); boundary = set(boundary_list)
    name = class_name(la, lb, lc, ld) + suffix

    # weight classes → the section names they own, and the raise-weight expression.
    wclasses = [("1", "1.0", [k for k in sections if k == "val" or k.startswith("dn")]),
                ("A", "wA", ["upA"] if "upA" in sections else []),
                ("B", "wB", ["upB"] if "upB" in sections else []),
                ("C", "wC", ["upC"] if "upC" in sections else [])]
    wclasses = [wc for wc in wclasses if wc[2]]

    # VRR tower liveness (boundary treated as outputs — no slot reuse).
    topo_vrr = [n for n in topo if n in vrr_nodes]
    vslot, vpeak = liveness_slots(topo_vrr, dag, boundary, reuse=True)

    L = [f"// SPLIT augmented value+gradient kernel ({name}) [§8.3/§8.4 + §5.1 split]: "
         f"{nout} outputs, {len(dag.nodes)} nodes, {nb} boundary, {len(wclasses)} weight-classes",
         f"__attribute__((noinline)) void recursum_gradctr_{name}(",
         "    const ScalarPack* __restrict__ sp,",
         "    const double* __restrict__ kfa,",
         "    const double* __restrict__ w2,   // cd x 3: 2*alpha for centres A,B,C",
         "    int cd, int mstride,",
         "    double* __restrict__ out,",
         "    double* __restrict__ sc) {",
         f"    for(int k=0;k<{nout};++k) out[k]=0.0;"]
    for wid, _w, _s in wclasses:
        L += [f"    double e0f0c_{wid}[{nb}];", f"    for(int k=0;k<{nb};++k) e0f0c_{wid}[k]=0.0;"]
    L += ["    for(int c=0;c<cd;++c) {",
          "        const ScalarPack& s = sp[c];",
          "        const double* kf = kfa + (long)c*mstride;",
          f"        double e0f0[{nb}];",
          "        const double wA=w2[(long)c*3+0], wB=w2[(long)c*3+1], wC=w2[(long)c*3+2];"]
    def vref(src):
        if src in bidx:        return f"e0f0[{bidx[src]}]"
        if not dag.nodes[src]: return f"kf[{src.m}]"
        return f"sc[{vslot[src]}]"
    for n in topo_vrr:
        if not dag.nodes[n] and n not in boundary:
            continue
        rhs = node_rhs(n, dag, vref)
        if n in bidx:    L.append(f"        e0f0[{bidx[n]}] = {rhs};")
        elif n in vslot: L.append(f"        sc[{vslot[n]}] = {rhs};")
    for wid, w, _s in wclasses:
        if w == "1.0": L.append(f"        for(int k=0;k<{nb};++k) e0f0c_{wid}[k]+=e0f0[k];")
        else:          L.append(f"        for(int k=0;k<{nb};++k) e0f0c_{wid}[k]+={w}*e0f0[k];")
    L += ["    }", "    { const ScalarPack& s = sp[0]; (void)s;"]

    # HRR once per weight class: reachable HRR nodes from that class's outputs.
    hpeak_max = 0
    for wid, _w, secnames in wclasses:
        cls_outs = [o for k in secnames for o in
                    outputs_list[sections[k][0]:sections[k][0]+sections[k][1]]]
        # backward reachability through hrr_nodes (stop at boundary/vrr).
        need = set(); stack = list(cls_outs)
        while stack:
            x = stack.pop()
            if x in need or x not in hrr_nodes: continue
            need.add(x)
            for t in dag.nodes[x]: stack.append(t.source)
        topo_h = [n for n in topo if n in need]
        hslot, hpk = liveness_slots(topo_h, dag, set(cls_outs), reuse=True)
        hpeak_max = max(hpeak_max, hpk)
        oi_cls = {o: out_index[o] for o in cls_outs}
        def href(src, wid=wid, hslot=hslot, oi_cls=oi_cls):
            if src in bidx:      return f"e0f0c_{wid}[{bidx[src]}]"
            if src in oi_cls:    return f"out[{oi_cls[src]}]"
            if not dag.nodes[src]: return f"e0f0c_{wid}[{bidx[src]}]"
            return f"sc[{hslot[src]}]"
        # boundary-only outputs (no HRR): out = weighted contracted boundary directly.
        for o in cls_outs:
            if o in boundary and o not in need:
                L.append(f"      out[{out_index[o]}] = e0f0c_{wid}[{bidx[o]}];")
        for n in topo_h:
            rhs = node_rhs(n, dag, href)
            if n in oi_cls:  L.append(f"      out[{oi_cls[n]}] = {rhs};")
            elif n in hslot: L.append(f"      sc[{hslot[n]}] = {rhs};")
    L += ["    }", "}"]
    return "\n".join(L), sections, nout, len(dag.nodes)


# ---------------------------------------------------------------------------
# Tier-1 FLAGSHIP: VRR/HRR contraction-boundary split (SHRIKE spec §5.1).
# The class DAG is partitioned so the primitive-dependent VRR tower runs per
# primitive while the geometry-only HRR runs ONCE per contracted quartet.
# Validated equivalent to the fused eval at Python level (plan.md L7).
# ---------------------------------------------------------------------------
def split_sets(dag, outputs):
    """Partition the DAG for the contraction-boundary split.

    Returns (vrr_nodes, hrr_nodes, boundary_list, bidx) where:
      vrr_nodes  : nodes with L_b==0 and L_d==0 (the (a,0|c,0)^(m) tower;
                   primitive-dependent — recomputed per primitive).
      hrr_nodes  : nodes with L_b>0 or L_d>0 (transfer A->B / C->D; geometry
                   scalars AB/CD only; all at m==0 — contraction-invariant).
      boundary   : the (a,0|c,0)^(0) intermediates consumed by an HRR node (or
                   the outputs themselves when lb==ld==0). Contracted across
                   primitives into e0f0_c[], then read once by the HRR kernel.
      bidx       : canonical {boundary node -> index} shared by both kernels."""
    out_set = set(outputs)
    vrr_nodes = {n for n in dag.nodes if n.L_b == 0 and n.L_d == 0}
    hrr_nodes = {n for n in dag.nodes if n.L_b > 0 or n.L_d > 0}
    boundary = set()
    for n in hrr_nodes:
        for t in dag.nodes[n]:
            s = t.source
            if s.L_b == 0 and s.L_d == 0 and s.m == 0:
                boundary.add(s)
    for o in outputs:
        if o.L_b == 0 and o.L_d == 0 and o.m == 0:
            boundary.add(o)
    boundary_list = sorted(boundary)  # Integral is order=True -> deterministic
    bidx = {n: i for i, n in enumerate(boundary_list)}
    return vrr_nodes, hrr_nodes, boundary_list, bidx


def emit_vrr_kernel(la, lb, lc, ld, suffix=""):
    """Emit recursum_vrr_<cls>(s, kf, e0f0, sc): builds the (a,0|c,0) tower for
    ONE primitive and writes the boundary values to e0f0[bidx]. Non-boundary
    tower nodes use caller scratch sc[] (Uniform, §10.6). Returns
    (source, n_boundary, vrr_peak)."""
    outputs = output_set(la, lb, lc, ld)
    dag = build_dag(outputs)
    topo = dag.topological_order()
    vrr_nodes, _hrr, boundary_list, bidx = split_sets(dag, outputs)
    boundary = set(boundary_list)
    topo_vrr = [n for n in topo if n in vrr_nodes]
    # Boundary nodes go to e0f0[]; treat as "outputs" for liveness (no slot reuse).
    slot, peak = liveness_slots(topo_vrr, dag, boundary, reuse=True)

    name = class_name(la, lb, lc, ld) + suffix
    def ref(src):
        if src in bidx:        return f"e0f0[{bidx[src]}]"
        if not dag.nodes[src]: return f"kf[{src.m}]"
        return f"sc[{slot[src]}]"
    L = [f"// VRR tower ({name}): {len(boundary_list)} boundary (a0|c0) ints, "
         f"peak {peak} slots  [per-primitive; SHRIKE split §5.1]",
         f"__attribute__((noinline)) void recursum_vrr_{name}(",
         "    const ScalarPack& __restrict__ s,",
         "    const double* __restrict__ kf,",
         "    double* __restrict__ e0f0,   // out: contracted-boundary layout",
         "    double* __restrict__ sc) {"]
    if peak == 0:
        L.append("    (void)sc;")
    for n in topo_vrr:
        if not dag.nodes[n] and n not in boundary:
            continue  # base -> kf[m] (never a boundary unless consumed by HRR)
        rhs = node_rhs(n, dag, ref)
        if n in bidx:
            L.append(f"    e0f0[{bidx[n]}] = {rhs};")
        elif n in slot:
            L.append(f"    sc[{slot[n]}] = {rhs};")
    L.append("}")
    return "\n".join(L), len(boundary_list), peak


def emit_hrr_kernel(la, lb, lc, ld, suffix=""):
    """Emit recursum_hrr_<cls>(s, e0f0, out, sc): transfers A->B / C->D ONCE on
    the contracted boundary e0f0[] (read like base integrals) to produce the
    class outputs. Geometry scalars (AB/CD) only. Returns (source, hrr_peak)."""
    outputs = output_set(la, lb, lc, ld)
    out_index = {o: i for i, o in enumerate(outputs)}
    dag = build_dag(outputs)
    topo = dag.topological_order()
    _vrr, hrr_nodes, boundary_list, bidx = split_sets(dag, outputs)
    topo_hrr = [n for n in topo if n in hrr_nodes]
    slot, peak = liveness_slots(topo_hrr, dag, set(outputs), reuse=True)

    name = class_name(la, lb, lc, ld) + suffix
    def ref(src):
        if src in bidx:          return f"e0f0[{bidx[src]}]"
        if src in out_index:     return f"out[{out_index[src]}]"
        return f"sc[{slot[src]}]"
    L = [f"// HRR transfer ({name}): {len(hrr_nodes)} nodes, peak {peak} slots  "
         f"[once/contracted-quartet; SHRIKE split §5.1]",
         f"__attribute__((noinline)) void recursum_hrr_{name}(",
         "    const ScalarPack& __restrict__ s,",
         "    const double* __restrict__ e0f0,   // in: contracted boundary",
         "    double* __restrict__ out,",
         "    double* __restrict__ sc) {"]
    if peak == 0:
        L.append("    (void)sc;")
    for n in topo_hrr:
        rhs = node_rhs(n, dag, ref)
        if n in out_index:
            L.append(f"    out[{out_index[n]}] = {rhs};")
        elif n in slot:
            L.append(f"    sc[{slot[n]}] = {rhs};")
    L.append("}")
    return "\n".join(L), peak


def emit_contracted_kernel(la, lb, lc, ld, suffix=""):
    """Emit recursum_ctr_<cls>(sp, kfa, cd, mstride, out, sc): the FULL contracted
    ERI with the primitive-contraction loop INSIDE one generated function — the
    libint2 `for(c<contrdepth)` structure, but with SHRIKE's DAG-CSE body kept
    register-resident and HRR done once on the contracted boundary. Replaces the
    C++ per-primitive orchestration in contract_eri (SHRIKE Engine Lever A), so
    the per-kernel advantage propagates instead of being spent on call overhead
    and the e0f0->e0f0c memory round-trip. `sp[c]`/`kfa+c*mstride` are the
    per-primitive scalars + prescaled Boys; `sc` is caller scratch. Returns
    (source, nout, nbnd, peak)."""
    outputs = output_set(la, lb, lc, ld)
    nout = len(outputs)
    out_index = {o: i for i, o in enumerate(outputs)}
    dag = build_dag(outputs)
    topo = dag.topological_order()
    name = class_name(la, lb, lc, ld) + suffix
    split = has_split(la, lb, lc, ld)
    L = [f"// Contracted kernel ({name}): in-kernel primitive loop (contrdepth) + "
         f"{'VRR-contract + HRR once' if split else 'fused-accumulate'} "
         f"[SHRIKE Engine Lever A]",
         f"__attribute__((noinline)) void recursum_ctr_{name}(",
         "    const ScalarPack* __restrict__ sp,",
         "    const double* __restrict__ kfa,",
         "    int cd, int mstride,",
         "    double* __restrict__ out,",
         "    double* __restrict__ sc) {"]
    if split:
        vrr_nodes, hrr_nodes, boundary_list, bidx = split_sets(dag, outputs)
        nbnd = len(boundary_list)
        boundary = set(boundary_list)
        topo_vrr = [n for n in topo if n in vrr_nodes]
        topo_hrr = [n for n in topo if n in hrr_nodes]
        vslot, vpeak = liveness_slots(topo_vrr, dag, boundary, reuse=True)
        hslot, hpeak = liveness_slots(topo_hrr, dag, set(outputs), reuse=True)
        L += [f"    double e0f0c[{nbnd}];",
              f"    for(int k=0;k<{nbnd};++k) e0f0c[k]=0.0;",
              "    for(int c=0;c<cd;++c) {",
              "        const ScalarPack& s = sp[c];",
              "        const double* kf = kfa + (long)c*mstride;",
              f"        double e0f0[{nbnd}];"]
        def vref(src):
            if src in bidx:        return f"e0f0[{bidx[src]}]"
            if not dag.nodes[src]: return f"kf[{src.m}]"
            return f"sc[{vslot[src]}]"
        for n in topo_vrr:
            if not dag.nodes[n] and n not in boundary:
                continue
            rhs = node_rhs(n, dag, vref)
            if n in bidx:       L.append(f"        e0f0[{bidx[n]}] = {rhs};")
            elif n in vslot:    L.append(f"        sc[{vslot[n]}] = {rhs};")
        L += [f"        for(int k=0;k<{nbnd};++k) e0f0c[k]+=e0f0[k];",
              "    }",
              "    { const ScalarPack& s = sp[0]; (void)s;"]
        def href(src):
            if src in bidx:      return f"e0f0c[{bidx[src]}]"
            if src in out_index: return f"out[{out_index[src]}]"
            if not dag.nodes[src]: return f"kf[{src.m}]"
            return f"sc[{hslot[src]}]"
        for n in topo_hrr:
            rhs = node_rhs(n, dag, href)
            if n in out_index:  L.append(f"      out[{out_index[n]}] = {rhs};")
            elif n in hslot:    L.append(f"      sc[{hslot[n]}] = {rhs};")
        L.append("    }")
        peak, nbnd_ret = max(vpeak, hpeak), nbnd
    else:
        slot, peak = liveness_slots(topo, dag, set(outputs), reuse=True)
        L += [f"    for(int k=0;k<{nout};++k) out[k]=0.0;",
              "    for(int c=0;c<cd;++c) {",
              "        const ScalarPack& s = sp[c];",
              "        const double* kf = kfa + (long)c*mstride;",
              f"        double o[{nout}];"]
        def fref(src):
            if src in out_index:   return f"o[{out_index[src]}]"
            if not dag.nodes[src]: return f"kf[{src.m}]"
            return f"sc[{slot[src]}]"
        for n in topo:
            if not dag.nodes[n] and n not in out_index:
                continue
            rhs = node_rhs(n, dag, fref)
            if n in out_index:  L.append(f"        o[{out_index[n]}] = {rhs};")
            elif n in slot:     L.append(f"        sc[{slot[n]}] = {rhs};")
        L += [f"        for(int k=0;k<{nout};++k) out[k]+=o[k];",
              "    }"]
        nbnd_ret = 0
    L.append("}")
    return "\n".join(L), nout, nbnd_ret, peak


def emit_contracted_kernel_simd(la, lb, lc, ld, V=4, suffix=""):
    """Emit recursum_ctrv_<cls>: the contracted kernel VECTORISED across primitive
    quartets (SHRIKE Engine Lever B). Lanes = V primitives; every ScalarPack field
    / Boys value / intermediate becomes a `recursum_vdbl` (GCC vector ext), so the
    SAME DAG-CSE body runs V-wide. libint2 is scalar (VECLEN=1) — this is the jump
    past it. The VRR contraction accumulates a vector boundary; lanes are
    horizontally reduced ONCE before the (scalar) HRR. `spv`/`kfa` are AoSoA blocks
    of V primitives; `sgeom` gives the (primitive-invariant) HRR geometry. Returns
    (source, nout)."""
    outputs = output_set(la, lb, lc, ld)
    nout = len(outputs)
    out_index = {o: i for i, o in enumerate(outputs)}
    dag = build_dag(outputs)
    topo = dag.topological_order()
    name = class_name(la, lb, lc, ld) + suffix
    split = has_split(la, lb, lc, ld)
    L = [f"// Contracted SIMD kernel ({name}): {V}-wide across primitives (GCC "
         f"vector ext); libint2 is scalar VECLEN=1 [SHRIKE Engine Lever B]",
         f"__attribute__((noinline)) void recursum_ctrv_{name}(",
         "    const ScalarPackV* __restrict__ spv,",
         "    const recursum_vdbl* __restrict__ kfa,",
         "    int nblk, int mstride,",
         "    const ScalarPack& __restrict__ sgeom,",
         "    double* __restrict__ out,",
         "    double* __restrict__ sc) {",
         "    const recursum_vdbl VZERO = {0.0,0.0,0.0,0.0};"]
    if split:
        vrr_nodes, hrr_nodes, boundary_list, bidx = split_sets(dag, outputs)
        nbnd = len(boundary_list)
        boundary = set(boundary_list)
        topo_vrr = [n for n in topo if n in vrr_nodes]
        topo_hrr = [n for n in topo if n in hrr_nodes]
        vslot, vpeak = liveness_slots(topo_vrr, dag, boundary, reuse=True)
        hslot, hpeak = liveness_slots(topo_hrr, dag, set(outputs), reuse=True)
        L += [f"    recursum_vdbl e0f0c[{nbnd}];",
              f"    for(int k=0;k<{nbnd};++k) e0f0c[k]=VZERO;",
              "    for(int b=0;b<nblk;++b) {",
              "        const ScalarPackV& s = spv[b];",
              "        const recursum_vdbl* kf = kfa + (long)b*mstride;",
              f"        recursum_vdbl e0f0[{nbnd}];"]
        if vpeak > 0:
            L.append(f"        recursum_vdbl scv[{vpeak}];")
        def vref(src):
            if src in bidx:        return f"e0f0[{bidx[src]}]"
            if not dag.nodes[src]: return f"kf[{src.m}]"
            return f"scv[{vslot[src]}]"
        for n in topo_vrr:
            if not dag.nodes[n] and n not in boundary:
                continue
            rhs = node_rhs(n, dag, vref)
            if n in bidx:    L.append(f"        e0f0[{bidx[n]}] = {rhs};")
            elif n in vslot: L.append(f"        scv[{vslot[n]}] = {rhs};")
        L += [f"        for(int k=0;k<{nbnd};++k) e0f0c[k]+=e0f0[k];",
              "    }",
              f"    double e0f0cs[{nbnd}];",
              f"    for(int k=0;k<{nbnd};++k){{ recursum_vdbl v=e0f0c[k]; double t=0;"
              f" for(int jL=0;jL<{V};++jL) t+=v[jL]; e0f0cs[k]=t; }}",
              "    { const ScalarPack& s = sgeom; (void)s;"]
        def href(src):
            if src in bidx:      return f"e0f0cs[{bidx[src]}]"
            if src in out_index: return f"out[{out_index[src]}]"
            return f"sc[{hslot[src]}]"
        for n in topo_hrr:
            rhs = node_rhs(n, dag, href)
            if n in out_index: L.append(f"      out[{out_index[n]}] = {rhs};")
            elif n in hslot:   L.append(f"      sc[{hslot[n]}] = {rhs};")
        L.append("    }")
    else:
        slot, peak = liveness_slots(topo, dag, set(outputs), reuse=True)
        L += ["    (void)sgeom; (void)sc;",
              f"    recursum_vdbl accv[{nout}];",
              f"    for(int k=0;k<{nout};++k) accv[k]=VZERO;",
              "    for(int b=0;b<nblk;++b) {",
              "        const ScalarPackV& s = spv[b];",
              "        const recursum_vdbl* kf = kfa + (long)b*mstride;",
              f"        recursum_vdbl o[{nout}];"]
        if peak > 0:
            L.append(f"        recursum_vdbl scv[{peak}];")
        def fref(src):
            if src in out_index:   return f"o[{out_index[src]}]"
            if not dag.nodes[src]: return f"kf[{src.m}]"
            return f"scv[{slot[src]}]"
        for n in topo:
            if not dag.nodes[n] and n not in out_index:
                continue
            rhs = node_rhs(n, dag, fref)
            if n in out_index: L.append(f"        o[{out_index[n]}] = {rhs};")
            elif n in slot:    L.append(f"        scv[{slot[n]}] = {rhs};")
        L += [f"        for(int k=0;k<{nout};++k) accv[k]+=o[k];",
              "    }",
              f"    for(int k=0;k<{nout};++k){{ recursum_vdbl v=accv[k]; double t=0;"
              f" for(int jL=0;jL<{V};++jL) t+=v[jL]; out[k]=t; }}"]
    L.append("}")
    return "\n".join(L), nout


def emit_simd_header() -> str:
    """recursum_eri_simd.h: the V-wide vector type + AoSoA ScalarPackV (fields
    mirror ScalarPack exactly, so the emitted vector body is textually identical
    to the scalar one)."""
    fields = ("PAx PAy PAz QCx QCy QCz WPx WPy WPz WQx WQy WQz "
              "ABx ABy ABz CDx CDy CDz inv_2zp inv_2zq inv_2zpq "
              "frac_q_over_pq frac_p_over_pq").split()
    decl = "    recursum_vdbl " + ", ".join(fields) + ";"
    return ("#pragma once\n"
            "// AUTO-GENERATED SIMD types for the contracted vector kernels "
            "(Lever B). DO NOT EDIT.\n"
            '#include "recursum_eri_scalars.h"\n'
            "// 4-wide double (AVX2 under -march=native; splits to SSE otherwise). "
            "Both are correct.\n"
            "// aligned(8): use UNALIGNED loads (vmovupd) so the AoSoA scratch need "
            "not be 32-byte aligned\n"
            "// (std::vector storage isn't guaranteed over-aligned). Negligible cost "
            "on AVX2+.\n"
            "typedef double recursum_vdbl __attribute__((vector_size(32), aligned(8)));\n"
            "#define RECURSUM_VLEN 4\n"
            "struct ScalarPackV {\n" + decl + "\n};\n")


def has_split(la, lb, lc, ld) -> bool:
    """The split only helps when there is an HRR stage (lb>0 or ld>0). For
    (a0|c0) classes it degenerates (boundary==outputs) — use the fused kernel."""
    return lb > 0 or ld > 0


# Canonical class ladder (a>=b, c>=d, a+b<=c+d) up to f-shells. Single source
# of truth: kernels, C class-table, and all dispatchers are generated from it.
CANON_LADDER = [
    (0,0,0,0), (0,0,1,0), (1,0,1,0), (0,0,1,1), (1,0,1,1), (1,1,1,1),
    (0,0,2,0), (2,0,2,0), (0,0,2,2), (2,2,2,2),
    (0,0,3,0), (3,0,3,0), (0,0,3,3), (3,3,3,3),
]

# Classes with more DAG nodes than this use chunked emission (many small
# noinline functions) instead of one monolithic function. Only ffff (77,756
# nodes) exceeds it; everything through dddd (7,623) stays single-function.
CHUNK_THRESHOLD = 10000


def emit_dispatch_header(classes, split_info=None, ctr_info=None, gradctr_info=None) -> str:
    """Generate eri_classes.h: the Cls[] table + RECURSUM/naive/name dispatch,
    shared by bench_eri.cpp, perf_driver.cpp and validate_kernels.cpp so the
    ladder is defined in exactly one place (CANON_LADDER).

    split_info: {class name -> (nbnd, vrr_ns, hrr_ns)} for classes with a Tier-1
    VRR/HRR split kernel (§5.1); drives the split dispatch + contraction driver."""
    split_info = split_info or {}
    ctr_info = ctr_info or {}
    gradctr_info = gradctr_info or {}
    rows, disp, vdisp = [], [], []
    for cls in classes:
        n = class_name(*cls)
        no = emit_kernel(*cls)[1]
        has_sp = 1 if n in split_info else 0
        nbnd = split_info[n][0] if n in split_info else 0
        rows.append(f'    {{"{n}",{{{cls[0]},{cls[1]},{cls[2]},{cls[3]}}},{no},'
                    f'recursum_eri_{n}_nscratch,{has_sp},{nbnd}}},')
        disp.append(f'    if(!strcmp(c,"{n}")) return recursum_eri_{n};')
        vdisp.append(f'    if(std::string(c)=="{n}") {{ recursum_eri_{n}(s,kf,out,sc); return; }}')
    # split dispatch: name -> vrr/hrr function pointers
    vrr_disp = "\n".join(
        f'    if(!strcmp(c,"{n}")) return recursum_vrr_{n};' for n in sorted(split_info))
    hrr_disp = "\n".join(
        f'    if(!strcmp(c,"{n}")) return recursum_hrr_{n};' for n in sorted(split_info))
    ctr_disp = "\n".join(
        f'    if(!strcmp(c,"{n}")) return recursum_ctr_{n};' for n in sorted(ctr_info))
    ctrv_disp = "\n".join(
        f'    if(!strcmp(c,"{n}")) return recursum_ctrv_{n};' for n in sorted(ctr_info))
    gradctr_disp = "\n".join(
        f'    if(!strcmp(c,"{n}")) return recursum_gradctr_{n};' for n in sorted(gradctr_info))
    grad_nout_disp = "\n".join(
        f'    if(!strcmp(c,"{n}")) return recursum_grad_{n}_nout;' for n in sorted(gradctr_info))
    grad_mmax_disp = "\n".join(
        f'    if(!strcmp(c,"{n}")) return recursum_grad_{n}_mmax;' for n in sorted(gradctr_info))
    return f"""#pragma once
// AUTO-GENERATED class table + ON dispatch (from CANON_LADDER). DO NOT EDIT.
#include <cstring>
#include <string>
#include "recursum_eri_scalars.h"
#include "recursum_eri_decls.h"    // extern kernel declarations (bodies in per-class TUs)
struct Cls {{ const char* name; int l[4]; int nout; int nscratch;
             int has_split; int nbnd; }};
static const Cls ERI_CLASSES[] = {{
{chr(10).join(rows)}
}};
static const int N_ERI_CLASSES = {len(classes)};
// Uniform scratch (§10.6): kernels take a caller-owned sc buffer sized by
// _nscratch; the max over all classes is RECURSUM_ERI_MAX_NSCRATCH (decls.h).
typedef void (*recursum_fn_t)(const ScalarPack&, const double*, double*, double*);
static inline recursum_fn_t recursum_dispatch(const char* c){{
{chr(10).join(disp)}
    return nullptr;
}}
// Convenience wrapper: owns a worst-case sc buffer so callers that don't manage
// an arena (test/bench harnesses) keep a simple 3-arg call. SHRIKE's Engine
// instead passes arena scratch directly to the 4-arg kernel symbol.
static inline void recursum_call(const char* c, const ScalarPack& s,
                                 const double* kf, double* out){{
    double sc[RECURSUM_ERI_MAX_NSCRATCH > 0 ? RECURSUM_ERI_MAX_NSCRATCH : 1];
{chr(10).join(vdisp)}
}}

// ---- Tier-1 VRR/HRR contraction-boundary split (§5.1) --------------------
// vrr: build (a0|c0) boundary for ONE primitive quartet -> e0f0[nbnd].
// hrr: transfer A->B / C->D ONCE on the contracted boundary -> out[nout].
typedef void (*recursum_vrr_t)(const ScalarPack&, const double*, double*, double*);
typedef void (*recursum_hrr_t)(const ScalarPack&, const double*, double*, double*);
static inline recursum_vrr_t recursum_vrr_dispatch(const char* c){{
{vrr_disp}
    return nullptr;  // null => no split for this class (use fused)
}}
static inline recursum_hrr_t recursum_hrr_dispatch(const char* c){{
{hrr_disp}
    return nullptr;
}}
// ---- Contracted kernels (SHRIKE Engine Lever A) --------------------------
// One call per contracted quartet: in-kernel primitive loop (contrdepth) with
// VRR-contract + HRR once (split classes) or fused-accumulate (a0|c0 classes).
// sp/kfa are arrays of `cd` per-primitive scalars/Boys; kfa is row-major with
// `mstride` doubles per primitive. sc = caller scratch (>= nscratch doubles).
typedef void (*recursum_ctr_t)(const ScalarPack*, const double*, int, int, double*, double*);
static inline recursum_ctr_t recursum_ctr_dispatch(const char* c){{
{ctr_disp}
    return nullptr;
}}
// SIMD contracted kernels (Lever B): V-wide across primitives (AoSoA ScalarPackV).
typedef void (*recursum_ctrv_t)(const ScalarPackV*, const recursum_vdbl*, int, int,
                                const ScalarPack&, double*, double*);
static inline recursum_ctrv_t recursum_ctrv_dispatch(const char* c){{
{ctrv_disp}
    return nullptr;
}}
// ---- Augmented value+gradient contracted kernels (§8.3/§8.4) --------------
// One call per contracted quartet -> value + shift sections in out[]; raised
// sections α-weighted via w2[cd*3] (2*alpha for centres A,B,C). Section offsets
// are computed caller-side from (la,lb,lc,ld). nout/mmax dispatched below.
typedef void (*recursum_gradctr_t)(const ScalarPack*, const double*, const double*,
                                   int, int, double*, double*);
static inline recursum_gradctr_t recursum_gradctr_dispatch(const char* c){{
{gradctr_disp}
    return nullptr;
}}
static inline int recursum_grad_nout(const char* c){{
{grad_nout_disp}
    return 0;
}}
static inline int recursum_grad_mmax(const char* c){{
{grad_mmax_disp}
    return 0;
}}
static inline const Cls* recursum_find_cls(const char* c){{
    for(int i=0;i<N_ERI_CLASSES;i++) if(!strcmp(ERI_CLASSES[i].name,c)) return &ERI_CLASSES[i];
    return nullptr;
}}
// Contracted driver: the SHRIKE §5.1 dispatch. For np primitive quartets,
// D==1 (segmented) uses the fused kernel BYTE-IDENTICALLY; np>1 with a split
// available runs VRR np times + HRR ONCE. `s`/`kf` are arrays of length np
// (per-primitive scalars + prescaled Boys); `sgeom` supplies AB/CD (identical
// across primitives). sc must hold >= max(nscratch, nbnd) doubles; e0f0_c holds
// >= nbnd doubles (caller-provided arena scratch).
static inline void recursum_contracted(const char* c, int np,
        const ScalarPack* s, const double* const* kf, const ScalarPack& sgeom,
        double* out, double* sc, double* e0f0_c){{
    const Cls* cl = recursum_find_cls(c);
    recursum_vrr_t vrr = recursum_vrr_dispatch(c);
    if(np == 1 || !vrr || !cl || !cl->has_split){{
        // segmented / no-split: sum the fused kernel over primitive quartets.
        for(int i=0;i<cl->nout;i++) out[i]=0.0;
        double tmp[RECURSUM_ERI_MAX_NSCRATCH>0?RECURSUM_ERI_MAX_NSCRATCH:1];
        recursum_fn_t fn = recursum_dispatch(c);
        double acc[16384];
        for(int p=0;p<np;p++){{
            fn(s[p], kf[p], acc, tmp);
            for(int i=0;i<cl->nout;i++) out[i]+=acc[i];
        }}
        return;
    }}
    recursum_hrr_t hrr = recursum_hrr_dispatch(c);
    for(int k=0;k<cl->nbnd;k++) e0f0_c[k]=0.0;
    double e0f0[16384];
    for(int p=0;p<np;p++){{
        vrr(s[p], kf[p], e0f0, sc);
        for(int k=0;k<cl->nbnd;k++) e0f0_c[k]+=e0f0[k];
    }}
    hrr(sgeom, e0f0_c, out, sc);
}}
"""


def emit_naive_dispatch(classes) -> str:
    skip = globals().get("_SKIP_NAIVE", set())
    dispn = []
    for cls in classes:
        n = class_name(*cls)
        if n in skip:
            continue  # no naive variant for this class (too large; ON-only)
        dispn.append(f'    if(!strcmp(c,"{n}")) return recursum_eri_{n}_naive;')
    return ("static inline recursum_fn_t recursum_dispatch_naive(const char* c){\n"
            + "\n".join(dispn) + "\n    return nullptr;  // null => no naive variant\n}\n")


def emit_header(classes: List[Tuple[int, int, int, int]]) -> str:
    parts = [f"""#pragma once
// AUTO-GENERATED by RECURSUM DAG-topological emitter (HGP-OS ERI).
// Peak-liveness slot reuse ported from osx emit.py.  DO NOT EDIT.
#include <recursum_eri_scalars.h>

#ifndef RECURSUM_FORCEINLINE
  #if defined(__GNUC__) || defined(__clang__)
    #define RECURSUM_FORCEINLINE inline __attribute__((always_inline))
  #else
    #define RECURSUM_FORCEINLINE inline
  #endif
#endif
"""]
    for cls in classes:
        parts.append(emit_kernel(*cls, reuse=True, suffix="")[0])           # liveness-ON
        parts.append(emit_kernel(*cls, reuse=False, suffix="_naive")[0])    # liveness-OFF
    return "\n\n".join(parts)


SCALAR_H = '#include "recursum_eri_scalars.h"\n'


def generate_project(classes, outdir=".", naive_cap=20000):
    """Emit the full project: one .cpp per (class,variant), a decls header, the
    dispatch header, and a manifest of the TU file names. Returns the list of
    generated .cpp basenames (for the build to compile in parallel).

    naive_cap: skip the liveness-OFF (_naive) variant for classes whose naive
    slot count exceeds this (a single giant function is pathological to compile
    and the ON-vs-naive comparison is degenerate there anyway -- e.g. ffff).
    The ON variant is always emitted."""
    import os
    tus = []
    decls = ["#pragma once", "// AUTO-GENERATED extern kernel declarations. DO NOT EDIT.",
             '#include "recursum_eri_scalars.h"', '#include "recursum_eri_simd.h"']
    skip_naive = set()
    nscratch = {}   # kernel name -> peak scratch slots (Uniform sc[], §10.6)
    for cls in classes:
        for reuse, suf in ((True, ""), (False, "_naive")):
            n_nodes = emit_kernel(*cls)[2]
            # Large classes (ffff-scale): chunked emission -- many small
            # noinline chunk functions instead of one pathological function
            # (fast compile AND faster runtime; see REPORT). ON variant only.
            if reuse and n_nodes > CHUNK_THRESHOLD:
                src, no, nn, frp, nch = emit_kernel_chunked(*cls, suffix=suf)
                nm = class_name(*cls) + suf
                with open(os.path.join(outdir, f"k_{nm}.cpp"), "w") as f:
                    f.write(src + "\n")   # chunked source already #includes scalars
                tus.append(f"k_{nm}.cpp")
                decls.append(f"void recursum_eri_{nm}(const ScalarPack&, const double*, double*, double*);")
                nscratch[nm] = frp
                continue
            if not reuse and n_nodes > CHUNK_THRESHOLD:
                skip_naive.add(class_name(*cls))
                continue  # no naive variant for chunked classes
            src, no, nn, pk = emit_kernel(*cls, reuse=reuse, suffix=suf)
            if not reuse and pk > naive_cap:
                skip_naive.add(class_name(*cls))
                continue  # ffff-scale: ON-only (see docstring)
            nm = class_name(*cls) + suf
            fn = f"k_{nm}.cpp"
            with open(os.path.join(outdir, fn), "w") as f:
                f.write(SCALAR_H + "\n" + src + "\n")
            tus.append(fn)
            decls.append(f"void recursum_eri_{nm}(const ScalarPack&, const double*, double*, double*);")
            nscratch[nm] = pk
    globals()["_SKIP_NAIVE"] = skip_naive

    # --- Tier-1 split kernels (§5.1): VRR tower + HRR transfer, emitted for
    # classes with a real HRR stage (lb>0 or ld>0). Oversized HRR (ffff-scale,
    # > CHUNK_THRESHOLD nodes) is deferred to fused-only for now — chunked split
    # emission is future work; those classes still work via the fused path.
    split_info = {}   # class name -> (nbnd, vrr_nscratch, hrr_nscratch)
    for cls in classes:
        if not has_split(*cls):
            continue
        n_nodes = emit_kernel(*cls)[2]
        if n_nodes > CHUNK_THRESHOLD:
            continue  # defer split for ffff-scale (single-function HRR too big)
        nm = class_name(*cls)
        vsrc, nbnd, vperk = emit_vrr_kernel(*cls)
        hsrc, hperk = emit_hrr_kernel(*cls)
        with open(os.path.join(outdir, f"k_{nm}_split.cpp"), "w") as f:
            f.write(SCALAR_H + "\n" + vsrc + "\n\n" + hsrc + "\n")
        tus.append(f"k_{nm}_split.cpp")
        decls.append(f"void recursum_vrr_{nm}(const ScalarPack&, const double*, double*, double*);")
        decls.append(f"void recursum_hrr_{nm}(const ScalarPack&, const double*, double*, double*);")
        nscratch[f"vrr_{nm}"] = vperk
        nscratch[f"hrr_{nm}"] = hperk
        split_info[nm] = (nbnd, vperk, hperk)
    globals()["_SPLIT_INFO"] = split_info

    # --- Contracted kernels (SHRIKE Engine Lever A): in-kernel primitive loop
    # (libint2 contrdepth structure) for EVERY class up to CHUNK_THRESHOLD. The
    # Engine calls ONE ctr kernel per contracted quartet instead of orchestrating
    # per-primitive VRR calls + accumulation in C++.
    ctr_info = {}   # class name -> ctr scratch peak
    for cls in classes:
        n_nodes = emit_kernel(*cls)[2]
        if n_nodes > CHUNK_THRESHOLD:
            continue  # ffff-scale: fused/C++-contracted path only
        nm = class_name(*cls)
        csrc, cno, cnbnd, cperk = emit_contracted_kernel(*cls)
        with open(os.path.join(outdir, f"k_{nm}_ctr.cpp"), "w") as f:
            f.write(SCALAR_H + "\n" + csrc + "\n")
        tus.append(f"k_{nm}_ctr.cpp")
        decls.append(f"void recursum_ctr_{nm}(const ScalarPack*, const double*, "
                     f"int, int, double*, double*);")
        nscratch[f"ctr_{nm}"] = cperk
        ctr_info[nm] = cperk
        # SIMD variant (Lever B): V-wide across primitives.
        vsrc, _vno = emit_contracted_kernel_simd(*cls)
        with open(os.path.join(outdir, f"k_{nm}_ctrv.cpp"), "w") as f:
            f.write('#include "recursum_eri_simd.h"\n' + vsrc + "\n")
        tus.append(f"k_{nm}_ctrv.cpp")
        decls.append(f"void recursum_ctrv_{nm}(const ScalarPackV*, const recursum_vdbl*, "
                     f"int, int, const ScalarPack&, double*, double*);")
    globals()["_CTR_INFO"] = ctr_info

    # --- Augmented value+gradient contracted kernels (§8.3/§8.4): ONE kernel per
    # class producing value + shift sections, raised sections α-weighted at
    # contraction. Guarded like ctr (skip ffff-scale augmented kernels; they need
    # the liveness chunker, future work). The C++ grad-Engine calls one per quartet.
    from eri_dag.dag import build_dag as _bdag
    gradctr_info = {}   # class name -> (nout, mmax, nscratch)
    for cls in classes:
        if emit_kernel(*cls)[2] > CHUNK_THRESHOLD:
            continue
        nm = class_name(*cls)
        # Task #13: SPLIT augmented gradient kernel (VRR/primitive + HRR once/class)
        # where the HRR stage is non-trivial; else the fused form. Same name+signature.
        if grad_has_split(*cls):
            gsrc, gsec, gnout, gnn = emit_gradctr_split_kernel(*cls)
        else:
            gsrc, gsec, gnout, gnn = emit_gradctr_kernel(*cls)
        _aug_outs = grad_output_layout(*cls)[0]
        _mmax = _bdag(_aug_outs).max_m()
        # sc sizing = FULL augmented-DAG peak (safe upper bound; the split's
        # max(vrr_peak, hrr_peak) is <= this, so the caller's scratch is sufficient).
        _slot, _gperk = liveness_slots(_bdag(_aug_outs).topological_order(),
                                       _bdag(_aug_outs), set(_aug_outs), reuse=True)
        with open(os.path.join(outdir, f"k_{nm}_gradctr.cpp"), "w") as f:
            f.write(SCALAR_H + "\n" + gsrc + "\n")
        tus.append(f"k_{nm}_gradctr.cpp")
        decls.append(f"void recursum_gradctr_{nm}(const ScalarPack*, const double*, "
                     f"const double*, int, int, double*, double*);")
        nscratch[f"gradctr_{nm}"] = max(_gperk, 1)
        gradctr_info[nm] = (gnout, _mmax, max(_gperk, 1))
    globals()["_GRADCTR_INFO"] = gradctr_info
    with open(os.path.join(outdir, "recursum_eri_simd.h"), "w") as f:
        f.write(emit_simd_header())

    # Per-kernel scratch sizes + the global MAX (spec §10.3 MAX_PEAK_LIVENESS:
    # the codegen-time constant a SHRIKE arena uses to size its scratch window).
    for nm in sorted(nscratch):
        decls.append(f"constexpr int recursum_eri_{nm}_nscratch = {nscratch[nm]};")
    for nm in sorted(split_info):
        decls.append(f"constexpr int recursum_{nm}_nbnd = {split_info[nm][0]};")
    for nm in sorted(gradctr_info):
        gn, gm, gs = gradctr_info[nm]
        decls.append(f"constexpr int recursum_grad_{nm}_nout = {gn};")
        decls.append(f"constexpr int recursum_grad_{nm}_mmax = {gm};")
    _max_ns = max(nscratch.values()) if nscratch else 0
    _max_nbnd = max((v[0] for v in split_info.values()), default=0)
    decls.append(f"constexpr int RECURSUM_ERI_MAX_NSCRATCH = {_max_ns};")
    decls.append(f"constexpr int RECURSUM_ERI_MAX_NBND = {_max_nbnd};")
    with open(os.path.join(outdir, "recursum_eri_decls.h"), "w") as f:
        f.write("\n".join(decls) + "\n")
    with open(os.path.join(outdir, "eri_classes.h"), "w") as f:
        f.write(emit_dispatch_header(classes, split_info, ctr_info, gradctr_info))
    # naive dispatch as its own tiny header (decls already cover the bodies)
    with open(os.path.join(outdir, "recursum_eri_naive_dispatch.h"), "w") as f:
        f.write('#pragma once\n#include "eri_classes.h"\n#include "recursum_eri_decls.h"\n'
                + emit_naive_dispatch(classes))
    with open(os.path.join(outdir, "kernel_tus.txt"), "w") as f:
        f.write("\n".join(tus) + "\n")
    return tus


if __name__ == "__main__":
    tus = generate_project(CANON_LADDER)
    print(f"generated {len(tus)} kernel TUs + decls + dispatch for "
          f"{len(CANON_LADDER)} classes")
    print(f"{'class':>6} {'nout':>6} {'nodes':>8} {'peak':>6} {'naive':>6} {'reduce':>7}")
    for cls in CANON_LADDER:
        _, no, nn, pk = emit_kernel(*cls, reuse=True)
        nv = emit_kernel(*cls, reuse=False)[3]
        print(f"{class_name(*cls):>6} {no:>6} {nn:>8} {pk:>6} {nv:>6} "
              f"{100*(1-pk/max(nv,1)):>6.0f}%")
