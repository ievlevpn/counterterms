"""Decorated rooted trees — the one concrete ``Symbol`` of the library.

``τ = b^n ⋆ ⨉ᵢ I^{(cᵢ)}_{p_i}(τ_i)`` (tourist_guide.tex 4506): a root of type
``b`` with node decoration ``n``, and a multiset of child edges, each an
integration operator ``I^{(c)}_p`` (component ``c``, derivative ``p``) over a
subtree.  Children are stored **canonically sorted**, so a frozen dataclass gives
canonical equality/hashing for free — the load-bearing correctness invariant
(equality, dedup, symmetry factor all depend on it).
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from math import factorial

from ..core.homogeneity import Homogeneity, MultiIndex

# An edge: (component, edge decoration p, subtree).
Edge = tuple[int, MultiIndex, "DecoratedTree"]


@dataclass(frozen=True)
class DecoratedTree:
    node_type: str
    node_dec: MultiIndex
    children: tuple[Edge, ...]  # canonically sorted

    def _sortkey(self):
        return (self.node_type, self.node_dec,
                tuple((c, p, sub._sortkey()) for (c, p, sub) in self.children))

    def homogeneity(self, sig) -> Homogeneity:
        h = sig.node_homogeneity(self.node_type) + Homogeneity(sig.scaled(self.node_dec))
        for (c, p, sub) in self.children:
            h = h + sig.edge_gain(c, p) + sub.homogeneity(sig)
        return h

    def symmetry_factor(self) -> int:
        # S(τ) = n! · Πⱼ S(σⱼ)^{mⱼ} · mⱼ!   (tourist_guide.tex 3982)
        s = 1
        for x in self.node_dec:
            s *= factorial(x)
        for (c, p, sub), m in Counter(self.children).items():
            s *= sub.symmetry_factor() ** m * factorial(m)
        return s

    def nodes(self) -> int:
        return 1 + sum(sub.nodes() for (_, _, sub) in self.children)


def tree(node_type: str, node_dec: MultiIndex, children=()) -> DecoratedTree:
    """Build a tree, canonicalising the child multiset by sorting."""
    ch = tuple(sorted(children, key=lambda e: (e[0], e[1], e[2]._sortkey())))
    return DecoratedTree(node_type, tuple(node_dec), ch)
