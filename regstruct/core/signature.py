"""The ``Signature`` — the typed vocabulary the whole library is parametric over.

It encodes every generality axis as data: dimension, scaling, the integration
operator(s) with their Schauder order, the noise types with their regularities,
and the structural rule (which child edges each node type admits, with caps).
Trees, generation and homogeneity all consult a ``Signature``; the algorithms
never branch on scalar-vs-system or one-vs-many-noises — those are just data here.

Phase 1 populates a single second-order operator and a scalar unknown.  The
fields already accommodate systems / multiple operators (Phase 2).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

from .homogeneity import Homogeneity, MultiIndex, Scaling

# An edge rule: (operator label, edge decoration p, max multiplicity or None).
EdgeRule = tuple[str, MultiIndex, "int | None"]


@dataclass
class Signature:
    dim: int
    scaling: Scaling
    op_order: int
    op_label: str
    noise_homog: dict[str, Homogeneity]
    node_types: tuple[str, ...]
    allowed: dict[str, tuple[EdgeRule, ...]]

    @property
    def width(self) -> int:
        """Length of a spacetime multi-index (1 time + dim space)."""
        return 1 + self.dim

    def node_homogeneity(self, node_type: str) -> Homogeneity:
        if node_type in self.noise_homog:
            return self.noise_homog[node_type]
        return Homogeneity(Fraction(0))  # ● and red nodes have homogeneity 0

    def edge_gain(self, p: MultiIndex) -> Homogeneity:
        # |I_p τ| = |τ| + (m − |p|_𝔰); m is the Schauder order (=2 for the heat kernel).
        return Homogeneity(Fraction(self.op_order - self.scaling.scaled(p)))

    def scaled(self, k: MultiIndex) -> int:
        return self.scaling.scaled(k)
