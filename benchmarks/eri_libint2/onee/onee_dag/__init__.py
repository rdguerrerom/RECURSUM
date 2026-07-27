"""One-electron Obara-Saika DAG (overlap S, kinetic T, nuclear V, multipole M).
Ported from osx-gpu gpu/codegen/onee for the RECURSUM HF integral stack —
the sibling of eri_dag (four-centre ERIs). See recursum_onee_emit.py."""
from .node import OneE, MomE, add_at, first_nonzero, all_ang_tuples
from .recurrence import Term, KinTerm, expand, kinetic_terms, expand_moment
from .dag import DAG, build_dag, evaluate
