"""The ``Signature`` — the typed vocabulary the whole library is parametric over.

It encodes every generality axis as data: dimension, scaling, the number of
components (unknowns), each component's integration-operator order, the noise
types with their regularities, and the structural rule (which child edges each
node type admits, with caps).  Trees, generation and homogeneity all consult a
``Signature``; scalar-vs-system and one-vs-many-noises are *data* here, not code
branches.

Edge / component identity (𝔗_e) rides on the **edge** — each child edge carries
the component index of the kernel that plants it (tourist_guide.tex 3827, 3965).
Node types stay ``{● , ∘_j}``.  The scaling is global; per-component operator
*order* is supported (mixed *scaling* is out of scope).
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .homogeneity import Homogeneity, MultiIndex, Scaling

# An edge rule: (component, edge decoration p, max multiplicity or None).
EdgeRule = tuple[int, MultiIndex, "int | None"]


@dataclass
class Signature:
    dim: int
    scaling: Scaling
    n_components: int
    comp_order: tuple[int, ...]               # integration-operator order per component
    noise_homog: dict[str, Homogeneity]
    node_types: tuple[str, ...]               # ('bullet',) + noise names
    allowed: dict[str, tuple[EdgeRule, ...]]  # node_type -> allowed (component, p, cap)

    @property
    def width(self) -> int:
        return 1 + self.dim

    def node_homogeneity(self, node_type: str) -> Homogeneity:
        if node_type in self.noise_homog:
            return self.noise_homog[node_type]
        return Homogeneity(Fraction(0))       # ● (and red) nodes have homogeneity 0

    def edge_gain(self, comp: int, p: MultiIndex) -> Homogeneity:
        # |I^{(c)}_p τ| = |τ| + (m_c − |p|_𝔰); m_c is component c's Schauder order.
        return Homogeneity(Fraction(self.comp_order[comp] - self.scaling.scaled(p)))

    def scaled(self, k: MultiIndex) -> int:
        return self.scaling.scaled(k)
