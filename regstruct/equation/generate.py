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
from typing import TYPE_CHECKING

from ..core.homogeneity import Homogeneity
from ..trees.tree import tree

if TYPE_CHECKING:
    from collections.abc import Callable, Iterator

    from ..core.homogeneity import MultiIndex, Scaling
    from ..core.signature import Signature
    from ..trees.tree import DecoratedTree

_BOUND = Fraction(1)          # keep trees with homogeneity.std < 1
_MAX_NODE_SCALED = 2          # safety cap on |n|_𝔰 of a node decoration
_MAX_POOL = 5000              # runaway-growth backstop (largest legitimate case ≈ 191 trees)


def _multiindices(length: int, scaling: Scaling, max_scaled: int) -> Iterator[MultiIndex]:
    def rec(i: int, rem: int, acc: list[int]) -> Iterator[MultiIndex]:
        if i == length:
            yield tuple(acc)
            return
        w = scaling.weights[i]
        ki = 0
        while w * ki <= rem:
            yield from rec(i + 1, rem - w * ki, acc + [ki])
            ki += 1
    yield from rec(0, max_scaled, [])


def generate_trees(sig: Signature, bound: Fraction = _BOUND) -> list[DecoratedTree]:
    """All rule-conforming decorated trees ``τ`` with ``homogeneity.std < bound`` — the
    γ-truncated basis of the model space ``T``.  ``generate_counterterms`` is the
    ``|τ|<0`` subset (the counterterm carriers); the rest are the positive sector."""
    decs = list(_multiindices(sig.width, sig.scaling, _MAX_NODE_SCALED))
    pool: dict[DecoratedTree, DecoratedTree] = {}

    def add(t: DecoratedTree) -> bool:
        if t in pool or t.homogeneity(sig).std >= bound:
            return False
        pool[t] = t
        if len(pool) > _MAX_POOL:      # fail fast instead of hanging (see note below)
            raise RuntimeError(
                f"tree generation exceeded {_MAX_POOL} trees — the negative-homogeneity "
                f"set is intractably large for this equation. The rule is subcritical (so "
                f"the set is finite in principle), but a small Schauder gain keeps too many "
                f"trees singular: typically high/fractional operator order combined with a "
                f"rich nonlinearity (e.g. quadratic in ∂u). See notes/swapping_the_operator.md.")
        return True

    for b in sig.node_types:           # seeds: bare decorated nodes
        for n in decs:
            add(tree(b, n, ()))

    changed = True
    while changed:                     # pool grows monotonically and is capped ⇒ terminates
        changed = False
        current = list(pool.values())

        for b in sig.node_types:
            rules = sig.allowed[b]
            caps = {(comp, p): cap for (comp, p, cap) in rules}
            atoms = [(comp, p, sub, (sig.edge_gain(comp, p) + sub.homogeneity(sig)))
                     for (comp, p, _cap) in rules for sub in current]
            atoms.sort(key=lambda a: (a[3].std, a[3].kap, a[0], a[1], a[2]._sortkey()))

            for n in decs:
                base_h = sig.node_homogeneity(b) + Homogeneity(sig.scaled(n))
                # NB: do NOT skip when base_h.std ≥ bound — capped negative-contribution
                # children (a gradient edge I_{(0,1)} over a noise contributes m−|p|+β₀ < 0
                # once β₀ < −(m−|p|), e.g. β₀=−3/2) can pull the full tree below the bound.
                # The DFS rejects the bare node and the budget `break` still terminates it
                # (negative-contribution children are derivative slots with a finite cap).
                changed |= _emit(b, n, base_h, atoms, caps, add, bound)

    return list(pool.values())


def generate_counterterms(sig: Signature) -> list[DecoratedTree]:
    """The divergent trees ``|τ| < 0`` — the counterterm-carrying subset of the basis."""
    return [t for t in generate_trees(sig) if t.homogeneity(sig).is_negative()]


def _emit(
    b: str,
    n: MultiIndex,
    base_h: Homogeneity,
    atoms: list[tuple[int, MultiIndex, DecoratedTree, Homogeneity]],
    caps: dict[tuple[int, MultiIndex], int | None],
    add: Callable[[DecoratedTree], bool],
    bound: Fraction = _BOUND,
) -> bool:
    """DFS over child multisets with homogeneity budget; returns whether pool grew."""
    grew = False
    counts: Counter[tuple[int, MultiIndex]] = Counter()

    def dfs(
        start: int,
        chosen: list[tuple[int, MultiIndex, DecoratedTree]],
        h: Homogeneity,
    ) -> None:
        nonlocal grew
        if add(tree(b, n, chosen)):
            grew = True
        for idx in range(start, len(atoms)):
            comp, p, sub, contrib = atoms[idx]
            if caps[(comp, p)] is not None and counts[(comp, p)] >= caps[(comp, p)]:
                continue
            nh = h + contrib
            if nh.std >= bound:           # atoms sorted ascending ⇒ rest overshoot too
                break
            counts[(comp, p)] += 1
            dfs(idx, chosen + [(comp, p, sub)], nh)
            counts[(comp, p)] -= 1

    dfs(0, [], base_h)
    return grew
