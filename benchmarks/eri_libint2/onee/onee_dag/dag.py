"""Backward DAG over the one-electron OS recurrences + a numerical
evaluator used to validate the rules against the reference (Step A2).

Mirrors gpu/codegen/eri/dag.py. The evaluator is NOT the device path; it
exists so the symbolic recurrence can be checked bit-for-bit against
reference.py before any CUDA is emitted.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Dict, List

from .node import OneE
from .recurrence import Term, expand


@dataclass
class DAG:
    op: str
    nodes: Dict[OneE, List[Term]] = field(default_factory=dict)

    def add(self, n: OneE) -> None:
        if n in self.nodes:
            return
        terms = expand(n, self.op)
        self.nodes[n] = terms
        for t in terms:
            self.add(t.source)

    def max_m(self) -> int:
        return max((n.m for n in self.nodes), default=0)

    def topological_order(self) -> List[OneE]:
        in_deg: Dict[OneE, int] = {n: 0 for n in self.nodes}
        consumers: Dict[OneE, List[OneE]] = {n: [] for n in self.nodes}
        for node, terms in self.nodes.items():
            for t in terms:
                in_deg[node] += 1
                consumers[t.source].append(node)
        ready = sorted([n for n, d in in_deg.items() if d == 0])
        out: List[OneE] = []
        while ready:
            n = ready.pop(0)
            out.append(n)
            for c in consumers[n]:
                in_deg[c] -= 1
                if in_deg[c] == 0:
                    _insort(ready, c)
        assert len(out) == len(self.nodes), "topological sort dropped nodes"
        return out


def _insort(lst: List[OneE], item: OneE) -> None:
    for idx, x in enumerate(lst):
        if x > item:
            lst.insert(idx, item)
            return
    lst.append(item)


def build_dag(outputs: List[OneE], op: str) -> DAG:
    dag = DAG(op=op)
    for o in outputs:
        dag.add(o)
    return dag


def evaluate(dag: DAG, scalars: Dict[str, float],
             base_eval: Callable[[OneE], float]) -> Dict[OneE, float]:
    """Evaluate every node bottom-up. `scalars` maps coef names
    (PAx, PCx, ABx, inv_2zeta, one) to numeric values; `base_eval`
    returns the value of a base node (a==b==0)."""
    vals: Dict[OneE, float] = {}
    for n in dag.topological_order():
        terms = dag.nodes[n]
        if not terms:
            vals[n] = base_eval(n)
            continue
        acc = 0.0
        for t in terms:
            acc += t.sign * t.imul * scalars[t.coef] * vals[t.source]
        vals[n] = acc
    return vals
