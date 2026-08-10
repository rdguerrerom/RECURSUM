"""Generic value-DAG front-end for the RECURSUM straight-line C emitter.

The ERI emitter (`recursum_dag_emit.emit_kernel`) is hardwired to the
McMurchie-Davidson / HGP-OS integral recurrence: it builds a FIXED value-DAG
of `Integral` nodes via `build_dag(output_set(...))` and emits CSE'd C. This
module adds a GENERAL value-DAG IR so arbitrary scalar value-DAGs — any straight
-line arithmetic over loads / constants / +,-,*,neg / stores — can be emitted as
CSE-optimized C through the SAME backend.

WHAT IS SHARED WITH THE ERI PATH.  The CSE is structural: a value-DAG is
hash-consed at build time (identical subexpressions intern to one node), then
`recursum_dag_emit.emit_cse_body` — the exact straight-line-emission driver the
ERI array/ssa kernels use — walks the topological order emitting each node once.
The array style additionally reuses `recursum_dag_emit.liveness_slots` (the osx
peak-liveness scratch-slot allocator) verbatim. This module reimplements NONE of
the CSE or C-emission logic; it only supplies a generic IR + `rhs`/`ref` closures.

Public API:
    g = ValueDAG()
    a = g.input("A", k)            # a load  A[k]
    c = g.const(0.5)              # a literal
    t = g.mul(a, c); s = g.add(t, ...)
    g.output("R", idx, s)         # a store R[idx] = <expr>
    src = emit_value_dag(g, inputs=["A"], outputs=["R"], fname="kern")
"""
from __future__ import annotations

from collections import namedtuple
from typing import Dict, List

from recursum_dag_emit import emit_cse_body, liveness_slots

_Term = namedtuple("_Term", ["source"])   # quacks like recurrence.Term for liveness_slots


class _Adapter:
    """Minimal DAG shim exposing `.nodes` so `liveness_slots` can be reused
    unchanged: nodes[id] -> list of `_Term(source=child_id)` ([] for a leaf)."""
    def __init__(self, nodes: Dict[int, List[_Term]]):
        self.nodes = nodes


class ValueDAG:
    """Hash-consed value-DAG builder. Every constructor returns an integer node
    id; structurally identical subexpressions return the SAME id (this IS the
    CSE). Light 0/1 constant folding keeps the emitted code — and the temp count
    — tight. Commutative ops (add, mul) canonicalize operand order so a+b and
    b+a share a node."""

    def __init__(self):
        self._intern: Dict[tuple, int] = {}
        self._specs: List[tuple] = []          # id -> spec tuple
        self._outputs: List[tuple] = []        # (out_node_id, name, index) in emit order

    # -- interning ---------------------------------------------------------
    def _get(self, spec: tuple) -> int:
        i = self._intern.get(spec)
        if i is None:
            i = len(self._specs)
            self._specs.append(spec)
            self._intern[spec] = i
        return i

    # -- leaves ------------------------------------------------------------
    def input(self, name: str, index: int) -> int:
        return self._get(("input", name, int(index)))

    def const(self, value: float) -> int:
        return self._get(("const", float(value)))

    # -- helpers for folding ----------------------------------------------
    def _const_val(self, i: int):
        s = self._specs[i]
        return s[1] if s[0] == "const" else None

    # -- arithmetic --------------------------------------------------------
    def add(self, a: int, b: int) -> int:
        av, bv = self._const_val(a), self._const_val(b)
        if av == 0.0: return b
        if bv == 0.0: return a
        if av is not None and bv is not None: return self.const(av + bv)
        a, b = sorted((a, b))                     # commutative canonicalization
        return self._get(("add", a, b))

    def sub(self, a: int, b: int) -> int:
        av, bv = self._const_val(a), self._const_val(b)
        if bv == 0.0: return a
        if av is not None and bv is not None: return self.const(av - bv)
        if a == b: return self.const(0.0)
        return self._get(("sub", a, b))

    def mul(self, a: int, b: int) -> int:
        av, bv = self._const_val(a), self._const_val(b)
        if av == 1.0: return b
        if bv == 1.0: return a
        if av == 0.0 or bv == 0.0: return self.const(0.0)
        if av is not None and bv is not None: return self.const(av * bv)
        a, b = sorted((a, b))
        return self._get(("mul", a, b))

    def neg(self, a: int) -> int:
        av = self._const_val(a)
        if av is not None: return self.const(-av)
        return self._get(("neg", a))

    # -- reductions --------------------------------------------------------
    def sum(self, xs: List[int]) -> int:
        xs = list(xs)
        if not xs: return self.const(0.0)
        acc = xs[0]
        for x in xs[1:]:
            acc = self.add(acc, x)
        return acc

    def dot(self, xs: List[int], ys: List[int]) -> int:
        assert len(xs) == len(ys)
        return self.sum([self.mul(x, y) for x, y in zip(xs, ys)])

    # -- stores ------------------------------------------------------------
    def output(self, name: str, index: int, val: int) -> int:
        oid = len(self._specs)
        self._specs.append(("output", name, int(index), val))
        self._outputs.append((oid, name, int(index)))
        return oid


# ---------------------------------------------------------------------------
def _reachable(g: ValueDAG) -> set:
    """Nodes reachable from any output (dead-code elimination)."""
    seen, stack = set(), [oid for oid, _, _ in g._outputs]
    while stack:
        n = stack.pop()
        if n in seen:
            continue
        seen.add(n)
        for c in _children(g._specs[n]):
            stack.append(c)
    return seen


def _children(spec) -> List[int]:
    op = spec[0]
    if op in ("input", "const"):
        return []
    if op == "neg":
        return [spec[1]]
    if op == "output":
        return [spec[3]]
    return [spec[1], spec[2]]              # add / sub / mul


def _topo_order(g: ValueDAG, keep: set) -> List[int]:
    """Sources-before-consumers, deterministic (child ids ascending)."""
    order, seen = [], set()
    def visit(n):
        if n in seen:
            return
        seen.add(n)
        for c in _children(g._specs[n]):
            if c in keep:
                visit(c)
        order.append(n)
    for n in sorted(keep):
        visit(n)
    return order


def cse_stats(g: ValueDAG) -> dict:
    """CSE bookkeeping: interned (deduped) node count, arithmetic op count with
    CSE, and the naive fully-inlined arithmetic op count (tree size) that a
    non-CSE expansion would emit — the ratio is the CSE win."""
    keep = _reachable(g)
    interior = [n for n in keep if g._specs[n][0] in ("add", "sub", "mul", "neg")]
    ops_cse = len(interior)
    # naive: expand each output into a tree, counting interior ops with sharing OFF
    from functools import lru_cache
    @lru_cache(maxsize=None)
    def tree_ops(n):
        sp = g._specs[n]
        op = sp[0]
        if op in ("input", "const"):
            return 0
        if op == "output":
            return tree_ops(sp[3])
        if op == "neg":
            return 1 + tree_ops(sp[1])
        return 1 + tree_ops(sp[1]) + tree_ops(sp[2])
    ops_naive = sum(tree_ops(oid) for oid, _, _ in g._outputs)
    return {"nodes": len(keep), "temps": ops_cse,
            "ops_cse": ops_cse, "ops_naive": ops_naive,
            "cse_factor": (ops_naive / ops_cse) if ops_cse else 1.0,
            "n_outputs": len(g._outputs)}


def _fmt_const(v: float) -> str:
    return repr(float(v))


def emit_value_dag(g: ValueDAG, inputs: List[str], outputs: List[str],
                   fname: str, style: str = "ssa", reuse: bool = True,
                   restrict: bool = True) -> str:
    """Emit CSE'd C for the value-DAG `g`.

    inputs/outputs: array names for the function signature (`const double*`
    inputs, `double*` outputs). style="ssa" -> `const double vN =` temporaries
    (compiler does register allocation; RECURSUM's fastest CPU variant). style=
    "array" -> a scratch `double t[peak]` with peak-liveness slot reuse (reuses
    `liveness_slots`). Returns the C source string."""
    keep = _reachable(g)
    topo = _topo_order(g, keep)
    out_nodes = {oid: (name, idx) for oid, name, idx in g._outputs}
    outputs_set = set(out_nodes)

    # adapter.nodes: leaves -> [], interior/output -> child _Terms (for liveness)
    nodes: Dict[int, List[_Term]] = {}
    for n in topo:
        sp = g._specs[n]
        if sp[0] in ("input", "const"):
            nodes[n] = []
        else:
            nodes[n] = [_Term(c) for c in _children(sp)]
    adapter = _Adapter(nodes)

    slot, peak = liveness_slots(topo, adapter, outputs_set, reuse=reuse)
    pos = {n: i for i, n in enumerate(topo)}

    def ref(x: int) -> str:
        sp = g._specs[x]
        op = sp[0]
        if op == "input":
            return f"{sp[1]}[{sp[2]}]"
        if op == "const":
            return _fmt_const(sp[1])
        if style == "array":
            return f"t[{slot[x]}]"
        return f"v{pos[x]}"

    def rhs_of(n: int) -> str:
        sp = g._specs[n]
        op = sp[0]
        if op == "output":
            return ref(sp[3])
        if op == "neg":
            return f"-{ref(sp[1])}"
        a, b = ref(sp[1]), ref(sp[2])
        return {"add": f"{a} + {b}", "sub": f"{a} - {b}", "mul": f"{a} * {b}"}[op]

    def out_fmt(n, r):
        name, idx = out_nodes[n]
        return f"    {name}[{idx}] = {r};"

    if style == "array":
        temp_fmt = lambda n, r: f"    t[{slot[n]}] = {r};"
    else:
        temp_fmt = lambda n, r: f"    const double v{pos[n]} = {r};"

    R = "__restrict__ " if restrict else ""
    params = [f"const double* {R}{nm}" for nm in inputs] + \
             [f"double* {R}{nm}" for nm in outputs]
    st = cse_stats(g)
    L = [f"// value-DAG kernel `{fname}` [{style}]: {st['n_outputs']} outputs, "
         f"{st['nodes']} DAG nodes, {st['ops_cse']} ops (CSE) vs {st['ops_naive']} "
         f"naive (x{st['cse_factor']:.2f}), peak {peak} slots",
         f"void {fname}({', '.join(params)}) {{"]
    if style == "array" and peak > 0:
        L.append(f"    double t[{peak}];")
    L += emit_cse_body(topo, nodes, outputs_set,
                       rhs_of=rhs_of, stored=lambda n: n in slot,
                       out_fmt=out_fmt, temp_fmt=temp_fmt)
    L.append("}")
    return "\n".join(L)
