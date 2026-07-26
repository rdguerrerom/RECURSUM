"""Backward DAG construction over the HGP-OS recurrence rules.

Starting from a list of output integrals (the L_a, L_b, L_c, L_d
bucket at m=0), expand each via `recurrence.expand` until all leaves
are base integrals `(0 0 | 0 0)^(m)`. Memoize: each unique
Integral is expanded at most once.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set

from .integral import Integral
from .recurrence import Term, expand


@dataclass
class DAG:
    """A directed acyclic graph of integrals.

    nodes: integral → list of Terms summing to it. Empty list ⇒ base integral.
    """
    nodes: Dict[Integral, List[Term]] = field(default_factory=dict)

    def add(self, integral: Integral) -> None:
        if integral in self.nodes:
            return
        terms = expand(integral)
        self.nodes[integral] = terms
        for t in terms:
            self.add(t.source)

    def base_integrals(self) -> List[Integral]:
        """All base (0000)^(m) integrals, sorted by m ascending."""
        bases = [n for n in self.nodes if n.is_base()]
        return sorted(bases, key=lambda i: i.m)

    def max_m(self) -> int:
        """Highest auxiliary index — drives how many F_n values we need."""
        return max(i.m for i in self.nodes)

    def topological_order(self) -> List[Integral]:
        """Topological sort with deterministic tie-break.

        Order: sources before consumers; within ties, sort by
        (L_total ascending, m descending, a, b, c, d). This emits
        base integrals first, then climbs angular momentum and falls
        in m as we go.
        """
        # Kahn's algorithm with deterministic tie-break.
        in_deg: Dict[Integral, int] = {n: 0 for n in self.nodes}
        consumers: Dict[Integral, List[Integral]] = {n: [] for n in self.nodes}
        for node, terms in self.nodes.items():
            for t in terms:
                in_deg[node] += 1
                consumers[t.source].append(node)

        # Ready set: nodes with no remaining sources (base integrals at start).
        ready = sorted(
            [n for n, d in in_deg.items() if d == 0],
            key=_dag_sort_key,
        )
        out: List[Integral] = []
        while ready:
            n = ready.pop(0)
            out.append(n)
            for cons in consumers[n]:
                in_deg[cons] -= 1
                if in_deg[cons] == 0:
                    # Insert in sorted position to preserve determinism.
                    _insort(ready, cons)
        assert len(out) == len(self.nodes), (
            f"topological sort dropped nodes: {len(out)} vs {len(self.nodes)}"
        )
        return out


def _dag_sort_key(i: Integral):
    """Tie-break for topological sort. Bases (L_total=0) first, then
    climb. Within the same L_total, lower m first (since lower m is
    consumed later in our recurrences)."""
    return (i.L_total, i.m, i.a, i.b, i.c, i.d)


def _insort(lst: List[Integral], item: Integral):
    """Insert into a sorted list. Linear; OK for synthetic-stage DAGs."""
    k = _dag_sort_key(item)
    for idx, x in enumerate(lst):
        if _dag_sort_key(x) > k:
            lst.insert(idx, item)
            return
    lst.append(item)


def build_dag(outputs: List[Integral]) -> DAG:
    """Convenience: build the full DAG from a list of output integrals."""
    dag = DAG()
    for o in outputs:
        dag.add(o)
    return dag
