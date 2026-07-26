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


def node_rhs(n: Integral, dag, reffn, kf="kf") -> str:
    """C expression for node n's value; reffn(src) names each source."""
    terms = dag.nodes[n]
    if not terms:  # base integral (0000)^(m)
        return f"{kf}[{n.m}]"
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

    # driver
    L.append(f"__attribute__((noinline)) void recursum_eri_{name}("
             "const ScalarPack& __restrict__ s, const double* __restrict__ kf, "
             "double* __restrict__ out) {")
    L.append(f"    double fr[{max(fr_peak,1)}];")
    for ci in range(len(chunks)):
        L.append(f"    {name}_c{ci}(s, kf, fr, out);")
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
    L = [f"// ERI class ({name}) [{style}]: {len(outputs_list)} Cartesian integrals, "
         f"{len(dag.nodes)} DAG nodes, peak liveness {peak} slots",
         f"__attribute__((noinline)) void recursum_eri_{name}(",
         "    const ScalarPack& __restrict__ s,",
         "    const double* __restrict__ kf,   // kf[m] = K * F_m(T)",
         "    double* __restrict__ out) {"]

    if style == "ssa":
        # name every intermediate node by topo index; compiler handles liveness
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
    if peak > 0:
        L.append(f"    double sc[{peak}];")
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


def emit_dispatch_header(classes) -> str:
    """Generate eri_classes.h: the Cls[] table + RECURSUM/naive/name dispatch,
    shared by bench_eri.cpp, perf_driver.cpp and validate_kernels.cpp so the
    ladder is defined in exactly one place (CANON_LADDER)."""
    rows, disp, vdisp = [], [], []
    for cls in classes:
        n = class_name(*cls)
        no = emit_kernel(*cls)[1]
        rows.append(f'    {{"{n}",{{{cls[0]},{cls[1]},{cls[2]},{cls[3]}}},{no}}},')
        disp.append(f'    if(!strcmp(c,"{n}")) return recursum_eri_{n};')
        vdisp.append(f'    if(std::string(c)=="{n}") {{ recursum_eri_{n}(s,kf,out); return; }}')
    return f"""#pragma once
// AUTO-GENERATED class table + ON dispatch (from CANON_LADDER). DO NOT EDIT.
#include <cstring>
#include <string>
#include "recursum_eri_scalars.h"
#include "recursum_eri_decls.h"    // extern kernel declarations (bodies in per-class TUs)
struct Cls {{ const char* name; int l[4]; int nout; }};
static const Cls ERI_CLASSES[] = {{
{chr(10).join(rows)}
}};
static const int N_ERI_CLASSES = {len(classes)};
typedef void (*recursum_fn_t)(const ScalarPack&, const double*, double*);
static inline recursum_fn_t recursum_dispatch(const char* c){{
{chr(10).join(disp)}
    return nullptr;
}}
static inline void recursum_call(const char* c, const ScalarPack& s,
                                 const double* kf, double* out){{
{chr(10).join(vdisp)}
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
             '#include "recursum_eri_scalars.h"']
    skip_naive = set()
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
                decls.append(f"void recursum_eri_{nm}(const ScalarPack&, const double*, double*);")
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
            decls.append(f"void recursum_eri_{nm}(const ScalarPack&, const double*, double*);")
    globals()["_SKIP_NAIVE"] = skip_naive
    with open(os.path.join(outdir, "recursum_eri_decls.h"), "w") as f:
        f.write("\n".join(decls) + "\n")
    with open(os.path.join(outdir, "eri_classes.h"), "w") as f:
        f.write(emit_dispatch_header(classes))
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
