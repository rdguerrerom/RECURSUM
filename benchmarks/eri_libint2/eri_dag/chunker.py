"""DAG chunker for register-budget-bounded subkernel emission.

The integral DAG for (ff|ff) is ~78k nodes — far beyond what fits in
~250 doubles of register working set per thread (the sm_89 budget
that keeps occupancy ≥ 25%). The chunker partitions the topological
order into sequential "chunks" such that each chunk's live-value set
(everything computed inside it but consumed later) stays under a
caller-specified budget.

Strategy (greedy topological partitioning):

  1. Walk topo order.
  2. Add each node to the current chunk's "interior" set.
  3. Track the chunk's "frontier" — interior nodes whose consumers
     are NOT all in the current chunk's interior (they're future
     consumers and must persist across chunk boundary).
  4. If frontier size > FRONTIER_BUDGET, close the current chunk
     and start a new one.

Frontier values are passed between chunks via a per-thread scratch
array in shared memory (or HBM for high-occupancy modes). Within a
chunk, all interior values are `const double` locals — register-
resident as long as the chunk size respects the budget.

Determinism: the topo order from dag.py is deterministic; the
chunker is a deterministic walk over it. Same DAG ⇒ same chunks ⇒
same scratch layout ⇒ same FP arithmetic (modulo the FMA-fusion
budget already documented).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Set

from .dag import DAG
from .integral import Integral


# Default budgets. FRONTIER_BUDGET caps the number of values that
# cross a chunk boundary; CHUNK_BUDGET caps the chunk's interior
# size. The product roughly determines the live-set inside the chunk:
# (frontier_in + interior + frontier_out) ≤ ~256 doubles for
# safety on sm_89 with 256 threads/block.

DEFAULT_INTERIOR_BUDGET     = 96     # nodes computed per chunk
DEFAULT_FRONTIER_BUDGET     = 64     # max boundary values OUT per chunk
DEFAULT_IMPORT_BUDGET       = 96     # max boundary values IN per chunk
DEFAULT_WORKING_SET_BUDGET  = 192    # inputs + interior (rough live-set bound)


@dataclass
class Chunk:
    """One subkernel's scope.

    chunk_id  : sequential index (0..n_chunks-1)
    interior  : ordered list of Integrals computed inside this chunk
    inputs    : Integrals consumed by this chunk that were computed in
                an earlier chunk (must be read from scratch)
    outputs   : interior Integrals that are consumed by any later chunk
                (must be written to scratch)
    """
    chunk_id: int
    interior: List[Integral] = field(default_factory=list)
    inputs:   List[Integral] = field(default_factory=list)
    outputs:  List[Integral] = field(default_factory=list)

    def __repr__(self):
        return (f"Chunk(id={self.chunk_id} "
                f"interior={len(self.interior)} "
                f"in={len(self.inputs)} out={len(self.outputs)})")


def partition(dag: DAG, topo: List[Integral],
              interior_budget: int = DEFAULT_INTERIOR_BUDGET,
              frontier_budget: int = DEFAULT_FRONTIER_BUDGET,
              import_budget:   int = DEFAULT_IMPORT_BUDGET,
              working_set_budget: int = DEFAULT_WORKING_SET_BUDGET,
              ) -> List[Chunk]:
    """Partition topo order into chunks under the given budgets.

    Live-set proxy (per thread, for nvcc register allocation):
       L  ≈  |chunk.inputs|  +  |chunk.interior|
    We close a chunk whenever ANY of the following exceeds budget:
      - interior_size > interior_budget           (interior nodes)
      - inputs_size   > import_budget             (cross-chunk imports)
      - outgoing frontier > frontier_budget       (cross-chunk exports)
      - inputs + interior > working_set_budget    (combined live set)

    The working_set_budget is the tightest constraint for sm_89's
    256-register-per-thread sweet spot: 256 doubles ≈ 256 registers
    in the worst case (no lifetime reuse).
    """
    # Build consumer lookup: integral → list of consumer integrals.
    consumers: Dict[Integral, List[Integral]] = {n: [] for n in dag.nodes}
    for node, terms in dag.nodes.items():
        for t in terms:
            consumers[t.source].append(node)

    # Position in topo for fast lookup.
    pos: Dict[Integral, int] = {n: i for i, n in enumerate(topo)}

    # For each node, find max consumer position (or -1 if no consumer).
    last_consumer_pos: Dict[Integral, int] = {}
    for n in dag.nodes:
        cs = consumers[n]
        if not cs:
            last_consumer_pos[n] = -1
        else:
            last_consumer_pos[n] = max(pos[c] for c in cs)

    chunks: List[Chunk] = []
    current = Chunk(chunk_id=0)
    # Set of integrals already computed in CLOSED chunks (must be
    # imported as inputs if used in this or later chunks).
    closed: Set[Integral] = set()
    # Set of interior integrals in current chunk.
    in_current: Set[Integral] = set()

    def project_frontier(chunk_end_pos: int) -> int:
        """Projected outgoing frontier size if we close after pos."""
        return sum(1 for n in in_current
                   if last_consumer_pos[n] > chunk_end_pos)

    for i, n in enumerate(topo):
        # Inputs needed from closed chunks.
        terms = dag.nodes[n]
        for t in terms:
            if t.source in closed and t.source not in current.inputs:
                current.inputs.append(t.source)
        # Add to interior.
        current.interior.append(n)
        in_current.add(n)
        # Should we close the chunk?
        close = False
        if len(current.interior) >= interior_budget:
            close = True
        elif len(current.inputs) > import_budget:
            close = True
        elif (len(current.inputs) + len(current.interior)
              > working_set_budget):
            close = True
        elif project_frontier(i) > frontier_budget:
            close = True
        # Never close on the last node — let it tail.
        if i == len(topo) - 1:
            close = True
        if close:
            # Compute outputs: interior values whose last consumer is
            # beyond the current chunk's end.
            current.outputs = [
                m for m in current.interior
                if last_consumer_pos[m] > i
            ]
            chunks.append(current)
            # Promote interior → closed.
            closed.update(in_current)
            # Start a new chunk.
            current = Chunk(chunk_id=len(chunks))
            in_current = set()

    return chunks


def chunk_stats(chunks: List[Chunk]) -> dict:
    """Aggregate metrics for inspection / chunker tuning."""
    return {
        "n_chunks": len(chunks),
        "max_interior": max((len(c.interior) for c in chunks), default=0),
        "max_in":       max((len(c.inputs)   for c in chunks), default=0),
        "max_out":      max((len(c.outputs)  for c in chunks), default=0),
        "total_interior": sum(len(c.interior) for c in chunks),
        "total_boundary": sum(len(c.outputs)  for c in chunks),
    }
