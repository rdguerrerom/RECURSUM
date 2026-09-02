"""RECURSUM DAG-topological emitter for FIXED-RADIX FFT CODELETS.

Sibling of recursum_moment_emit.py. The recurrence here is Cooley-Tukey:

    DFT_R(x)[k]       = DFT_{R/2}(x_even)[k] + W_R^k DFT_{R/2}(x_odd)[k]
    DFT_R(x)[k + R/2] = DFT_{R/2}(x_even)[k] - W_R^k DFT_{R/2}(x_odd)[k]

which bottoms out at DFT_1(x) = x. For a FIXED R the recursion unrolls to a
straight-line graph in which every twiddle is a compile-time constant, and for
R <= 8 the only irrational one is 1/sqrt(2). That is what makes a codelet worth
emitting rather than looping: the generic transform must load a twiddle per
butterfly, while the codelet has them in the instruction stream.

WHY THIS IS A C12 TARGET AND THE OUTER PASSES ARE NOT. C12 asks for a
straight-line DAG over extents known at generation time. A transform of runtime
length n is not that -- it is a quadrature over a grid and stays hand-written.
A radix-8 butterfly is exactly that, and gets one specialised kernel per radix
with no branch on the radix inside the body.

C13, ONE DAG BOTH DIRECTIONS. The inverse DFT is the conjugate of the forward
DFT of the conjugate, so forward and inverse codelets are the SAME graph with
the sign of every imaginary twiddle component flipped. They are emitted from
one build of the graph rather than from two specifications, so they cannot
drift apart.

Everything else follows the established pattern: hash-consed expression DAG so
common subexpressions are shared by construction, topological emission with
every node written exactly once, and a peak-liveness slot allocator sizing a
single scratch array.
"""
from __future__ import annotations

import math
from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------
# A hash-consed real-valued expression DAG.
#
# Hash-consing IS the common-subexpression elimination: two structurally equal
# expressions are the same object, so a shared butterfly is built once no
# matter how many outputs reach it.
# ---------------------------------------------------------------------------

_POOL: Dict[tuple, "Node"] = {}


class Node:
    __slots__ = ("op", "args", "const", "uid")

    def __init__(self, op, args, const, uid):
        self.op = op          # 'in' | 'add' | 'sub' | 'mul' | 'neg' | 'zero'
        self.args = args
        self.const = const
        self.uid = uid

    def __repr__(self):
        return f"<{self.op}#{self.uid}>"


def _mk(op, args=(), const=None) -> Node:
    key = (op, tuple(a.uid for a in args), const)
    n = _POOL.get(key)
    if n is None:
        n = Node(op, tuple(args), const, len(_POOL))
        _POOL[key] = n
    return n


def reset_pool():
    _POOL.clear()


ZERO = None


def zero() -> Node:
    return _mk("zero", (), 0.0)


def inp(name: str) -> Node:
    return _mk("in", (), name)


def add(a: Node, b: Node) -> Node:
    if a.op == "zero":
        return b
    if b.op == "zero":
        return a
    return _mk("add", (a, b))


def sub(a: Node, b: Node) -> Node:
    if b.op == "zero":
        return a
    if a.op == "zero":
        return neg(b)
    return _mk("sub", (a, b))


def neg(a: Node) -> Node:
    if a.op == "zero":
        return a
    if a.op == "neg":
        return a.args[0]
    return _mk("neg", (a,))


def mul(a: Node, c: float) -> Node:
    if a.op == "zero" or c == 0.0:
        return zero()
    if c == 1.0:
        return a
    if c == -1.0:
        return neg(a)
    return _mk("mul", (a,), c)


# ---------------------------------------------------------------------------
# Complex arithmetic over pairs of real nodes.
# ---------------------------------------------------------------------------

def cadd(x, y):
    return (add(x[0], y[0]), add(x[1], y[1]))


def csub(x, y):
    return (sub(x[0], y[0]), sub(x[1], y[1]))


def cmul_const(x, c: complex):
    """(a + i b)(c + i s). Exact special cases first: for R <= 8 most twiddles
    are 0, +-1 or +-i, and emitting a multiply by one of those would be pure
    waste in the instruction stream the codelet exists to control."""
    cr, ci = c.real, c.imag
    # snap values that are exactly representable to kill floating-point fuzz
    for v in (0.0, 1.0, -1.0):
        if abs(cr - v) < 1e-15:
            cr = v
        if abs(ci - v) < 1e-15:
            ci = v
    a, b = x
    if ci == 0.0:
        return (mul(a, cr), mul(b, cr))
    if cr == 0.0:
        # (a + ib) * i s = -b s + i a s
        return (mul(neg(b), ci), mul(a, ci))
    return (sub(mul(a, cr), mul(b, ci)), add(mul(a, ci), mul(b, cr)))


def dft(xs: List[Tuple[Node, Node]], sign: int) -> List[Tuple[Node, Node]]:
    """Cooley-Tukey on a fixed, power-of-two length. `sign` is -1 for the
    forward transform and +1 for the inverse; it enters ONLY through the
    twiddle, which is why one graph serves both directions."""
    R = len(xs)
    if R == 1:
        return list(xs)
    ev = dft(xs[0::2], sign)
    od = dft(xs[1::2], sign)
    out = [None] * R
    for k in range(R // 2):
        th = sign * 2.0 * math.pi * k / R
        t = cmul_const(od[k], complex(math.cos(th), math.sin(th)))
        out[k] = cadd(ev[k], t)
        out[k + R // 2] = csub(ev[k], t)
    return out


# ---------------------------------------------------------------------------
# Topological order and peak-liveness slot allocation.
# ---------------------------------------------------------------------------

def topo_order(roots: List[Node]) -> List[Node]:
    seen, order = set(), []

    def visit(n):
        if n.uid in seen:
            return
        seen.add(n.uid)
        for a in n.args:
            visit(a)
        order.append(n)

    for r in roots:
        visit(r)
    return order


def liveness_slots(topo: List[Node], roots: List[Node]):
    """Identical in spirit to the moment emitter's allocator: a value's slot is
    returned to the free list one step after its last consumer, so the scratch
    array is sized by PEAK liveness rather than by node count."""
    pos = {n.uid: i for i, n in enumerate(topo)}
    last_use: Dict[int, int] = {}
    for n in topo:
        for a in n.args:
            last_use[a.uid] = max(last_use.get(a.uid, -1), pos[n.uid])
    die_at: Dict[int, List[Node]] = {}
    for n in topo:
        p = last_use.get(n.uid)
        if p is not None:
            die_at.setdefault(p, []).append(n)
    rootset = {r.uid for r in roots}
    slot: Dict[int, int] = {}
    free: List[int] = []
    nxt = 0
    for i, n in enumerate(topo):
        for v in die_at.get(i - 1, []):
            if v.uid in slot:
                free.append(slot[v.uid])
        if n.op in ("in", "zero") or n.uid in rootset:
            continue
        if free:
            slot[n.uid] = free.pop()
        else:
            slot[n.uid] = nxt
            nxt += 1
    return slot, nxt


# ---------------------------------------------------------------------------
# Emission
# ---------------------------------------------------------------------------

def _fmt(c: float) -> str:
    if abs(c - math.sqrt(0.5)) < 1e-15:
        return "SCDT_FFT_RSQ2"
    if abs(c + math.sqrt(0.5)) < 1e-15:
        return "-SCDT_FFT_RSQ2"
    return repr(c)


def emit_codelet(R: int, sign: int, name: str) -> Tuple[str, int, int]:
    reset_pool()
    xs = [(inp(f"xr[{i}]"), inp(f"xi[{i}]")) for i in range(R)]
    outs = dft(xs, sign)
    roots = []
    for (a, b) in outs:
        roots.append(a)
        roots.append(b)
    topo = topo_order(roots)
    slot, peak = liveness_slots(topo, roots)

    ref: Dict[int, str] = {}
    lines: List[str] = []
    nops = 0

    def r(n: Node) -> str:
        return ref[n.uid]

    for n in topo:
        if n.op == "zero":
            ref[n.uid] = "0.0"
            continue
        if n.op == "in":
            ref[n.uid] = n.const
            continue
        if n.op == "add":
            rhs = f"{r(n.args[0])} + {r(n.args[1])}"
        elif n.op == "sub":
            rhs = f"{r(n.args[0])} - {r(n.args[1])}"
        elif n.op == "neg":
            rhs = f"-{r(n.args[0])}"
        elif n.op == "mul":
            rhs = f"{_fmt(n.const)} * {r(n.args[0])}"
        else:
            raise AssertionError(n.op)
        nops += 1
        if n.uid in slot:
            tgt = f"t[{slot[n.uid]}]"
            ref[n.uid] = tgt
            lines.append(f"    {tgt} = {rhs};")
        else:
            ref[n.uid] = f"__r{n.uid}"
            lines.append(f"    const double __r{n.uid} = {rhs};")

    # Outputs are written last, after every read of the inputs, so the codelet
    # is safe in place.
    tail = []
    for i, (a, b) in enumerate(outs):
        tail.append(f"    __o{2*i} = {r(a)};")
        tail.append(f"    __o{2*i+1} = {r(b)};")
    decl = "    double " + ", ".join(f"__o{i}" for i in range(2 * R)) + ";"
    store = []
    for i in range(R):
        store.append(f"    xr[{i}] = __o{2*i}; xi[{i}] = __o{2*i+1};")

    body = []
    if peak:
        body.append(f"    double t[{peak}];")
    body.append(decl)
    body += lines
    body += tail
    body += store

    src = (f"/* radix-{R} {'forward' if sign < 0 else 'inverse'} codelet: "
           f"{nops} operations, peak liveness {peak}.\n"
           f" *\n"
           f" * STATIC INLINE, IN THE HEADER, DELIBERATELY. Emitted as an\n"
           f" * ordinary function the kernel is a call per butterfly, and its\n"
           f" * two argument arrays have to live in memory across it: measured\n"
           f" * that way the codelet transform was SLOWER than the radix-2\n"
           f" * loop it replaced, 0.95 ns per butterfly against 0.75. Inlined,\n"
           f" * the arrays never escape and stay in registers, which is the\n"
           f" * whole reason a straight-line codelet is worth generating. */\n"
           f"static inline void {name}(double *restrict xr, "
           f"double *restrict xi)\n{{\n"
           + "\n".join(body) + "\n}\n")
    return src, nops, peak


def generate(radices=(2, 4, 8), prefix="scdt_fft_dft"):
    hdr = ["/* AUTO-GENERATED by recursum_fft_emit.py -- do not edit.",
           " *",
           " * Fixed-radix DFT codelets: straight-line, twiddles as literal",
           " * constants, emitted from the Cooley-Tukey recurrence with",
           " * hash-consed common subexpressions and peak-liveness slots.",
           " *",
           " * Forward and inverse are the SAME graph with the imaginary part",
           " * of every twiddle negated (C13), so they cannot drift apart.",
           " */",
           "#ifndef SCDT_FFT_CODELETS_H",
           "#define SCDT_FFT_CODELETS_H", "",
           "#ifdef __cplusplus", 'extern "C" {', "#endif", "",
           "#define SCDT_FFT_RSQ2 0.70710678118654752440084436210485", ""]
    src = ["/* AUTO-GENERATED by recursum_fft_emit.py -- do not edit.",
           " *",
           " * The codelets are static inline in the header so they can be",
           " * inlined into the pass loops; this file exists so the generated",
           " * unit still appears in the build and in PROVENANCE. */",
           '#include "fft_codelets.h"', "",
           "int scdt_fft_codelets_present(void) { return 1; }", ""]
    stats = []
    for R in radices:
        for sign, tag in ((-1, "fwd"), (+1, "inv")):
            name = f"{prefix}{R}_{tag}"
            code, nops, peak = emit_codelet(R, sign, name)
            hdr.append(code)
            _ = code
            stats.append((R, tag, nops, peak))
    hdr += ["", "int scdt_fft_codelets_present(void);", "",
            "#ifdef __cplusplus", "}", "#endif", "#endif", ""]
    return "\n".join(hdr), "\n".join(src), stats


if __name__ == "__main__":
    import sys
    outdir = sys.argv[1] if len(sys.argv) > 1 else "."
    h, c, stats = generate()
    open(f"{outdir}/fft_codelets.h", "w").write(h)
    open(f"{outdir}/fft_codelets.c", "w").write(c)
    for R, tag, nops, peak in stats:
        print(f"  radix {R:2d} {tag}: {nops:3d} operations, peak liveness {peak}")
