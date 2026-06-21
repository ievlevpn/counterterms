"""Enumerate the decorated trees that carry counterterms: ``{τ : |τ| < 0}``.

Bounded fixpoint generation.  We grow a pool of valid trees with homogeneity
standard-part ``< 1`` (subtrees of any negative tree have standard-part ``< 0``,
so this bound is complete), respecting the structural rule (which child edges a
node admits — `(component, p)` with a cap on derivative-edge multiplicity = the
nonlinearity's degree in that slot — what makes "g quadratic in ∂u" terminate).

Children are selected by a **budget-pruned DFS**: edge atoms are sorted by their
homogeneity contribution, so once one overshoots the budget all later ones do too
(``break``).  Termination: an uncapped edge (a field slot, ``p = 0``) contributes
``> 0`` to the standard part when ``β₀ > −2``, so it is bounded by the budget;
every edge that can contribute ``≤ 0`` is a derivative slot with a finite cap.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction

from ..core.homogeneity import Homogeneity
from ..trees.tree import tree

_BOUND = Fraction(1)          # keep trees with homogeneity.std < 1
_MAX_NODE_SCALED = 2          # safety cap on |n|_𝔰 of a node decoration


def _multiindices(length: int, scaling, max_scaled: int):
    def rec(i, rem, acc):
        if i == length:
            yield tuple(acc)
            return
        w = scaling.weights[i]
        ki = 0
        while w * ki <= rem:
            yield from rec(i + 1, rem - w * ki, acc + [ki])
            ki += 1
    yield from rec(0, max_scaled, [])


def generate_counterterms(sig):
    decs = list(_multiindices(sig.width, sig.scaling, _MAX_NODE_SCALED))
    pool: dict = {}

    def add(t) -> bool:
        if t in pool or t.homogeneity(sig).std >= _BOUND:
            return False
        pool[t] = t
        return True

    for b in sig.node_types:           # seeds: bare decorated nodes
        for n in decs:
            add(tree(b, n, ()))

    guard = 0
    changed = True
    while changed:
        changed = False
        guard += 1
        assert guard < 50, "tree generation did not terminate (subcriticality?)"
        current = list(pool.values())

        for b in sig.node_types:
            rules = sig.allowed[b]
            caps = {(comp, p): cap for (comp, p, cap) in rules}
            atoms = [(comp, p, sub, (sig.edge_gain(comp, p) + sub.homogeneity(sig)))
                     for (comp, p, _cap) in rules for sub in current]
            atoms.sort(key=lambda a: (a[3].std, a[3].kap, a[0], a[1], a[2]._sortkey()))

            for n in decs:
                base_h = sig.node_homogeneity(b) + Homogeneity(sig.scaled(n))
                if base_h.std >= _BOUND:
                    continue
                changed |= _emit(b, n, base_h, atoms, caps, add)

    return [t for t in pool.values() if t.homogeneity(sig).is_negative()]


def _emit(b, n, base_h, atoms, caps, add) -> bool:
    """DFS over child multisets with homogeneity budget; returns whether pool grew."""
    grew = False
    counts: Counter = Counter()

    def dfs(start, chosen, h):
        nonlocal grew
        if add(tree(b, n, chosen)):
            grew = True
        for idx in range(start, len(atoms)):
            comp, p, sub, contrib = atoms[idx]
            if caps[(comp, p)] is not None and counts[(comp, p)] >= caps[(comp, p)]:
                continue
            nh = h + contrib
            if nh.std >= _BOUND:          # atoms sorted ascending ⇒ rest overshoot too
                break
            counts[(comp, p)] += 1
            dfs(idx, chosen + [(comp, p, sub)], nh)
            counts[(comp, p)] -= 1

    dfs(0, [], base_h)
    return grew
