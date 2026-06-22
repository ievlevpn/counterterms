"""Decorated rooted trees — the one concrete ``Symbol`` of the library.

``τ = b^n ⋆ ⨉ᵢ I^{(cᵢ)}_{p_i}(τ_i)`` (tourist_guide.tex 4506): a root of type
``b`` with node decoration ``n``, and a multiset of child edges, each an
integration operator ``I^{(c)}_p`` (component ``c``, derivative ``p``) over a
subtree.  Children are stored **canonically sorted**, so a frozen dataclass gives
canonical equality/hashing for free — the load-bearing correctness invariant
(equality, dedup, symmetry factor all depend on it).

Phase 3 adds the *extended decoration*: a node ``color`` (``black`` for ordinary
nodes, ``red`` for contracted nodes produced by ``Δ⁻``, ``blue`` for ``T⁺`` roots)
and, on red nodes, an ``o`` value (``∈ ℤ[β₀]``, reusing ``Homogeneity``) equal to
the naive homogeneity of the extracted component (tourist_guide.tex 5495–5502).
Phase 1–2 trees are all black with ``o = 0``, so their behaviour is unchanged.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from math import factorial

from ..core.homogeneity import Homogeneity, MultiIndex

# An edge: (component, edge decoration p, subtree).
Edge = tuple[int, MultiIndex, "DecoratedTree"]

_ZERO = Homogeneity(0)


@dataclass(frozen=True)
class DecoratedTree:
    node_type: str
    node_dec: MultiIndex
    children: tuple[Edge, ...]            # canonically sorted
    color: str = "black"                  # 'black' | 'red' | 'blue'
    o: Homogeneity = field(default_factory=lambda: _ZERO)   # extended decoration (red nodes)

    def _sortkey(self):
        return (self.node_type, self.color, (self.o.std, self.o.kap), self.node_dec,
                tuple((c, p, sub._sortkey()) for (c, p, sub) in self.children))

    def homogeneity(self, sig) -> Homogeneity:
        """Naive homogeneity ``|τ|'`` (ignores the o-decoration; red nodes count as 0)."""
        base = _ZERO if self.color == "red" else sig.node_homogeneity(self.node_type)
        h = base + Homogeneity(sig.scaled(self.node_dec))
        for (c, p, sub) in self.children:
            h = h + sig.edge_gain(c, p) + sub.homogeneity(sig)
        return h

    def _o_sum(self) -> Homogeneity:
        s = self.o if self.color == "red" else _ZERO
        for (_c, _p, sub) in self.children:
            s = s + sub._o_sum()
        return s

    def extended_homogeneity(self, sig) -> Homogeneity:
        """Extended homogeneity ``|τ|`` = naive + Σ over red nodes of o (tex 5502)."""
        return self.homogeneity(sig) + self._o_sum()

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


def tree(node_type: str, node_dec: MultiIndex, children=(),
         color: str = "black", o: Homogeneity = _ZERO) -> DecoratedTree:
    """Build a tree, canonicalising the child multiset by sorting."""
    ch = tuple(sorted(children, key=lambda e: (e[0], e[1], e[2]._sortkey())))
    return DecoratedTree(node_type, tuple(node_dec), ch, color, o)


def red_node(o: Homogeneity, node_dec: MultiIndex = None, width: int = None) -> DecoratedTree:
    """A contracted (red) node with extended decoration ``o`` and no children."""
    if node_dec is None:
        node_dec = (0,) * (width if width is not None else 1)
    return tree("bullet", node_dec, (), color="red", o=o)
