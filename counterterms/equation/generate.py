"""Enumerate the decorated trees that carry counterterms: ``{τ : |τ| < 0}``.

Bounded fixpoint generation.  We grow a pool of valid trees with homogeneity
standard-part ``< 1`` (subtrees of any negative tree have standard-part ``< 0``,
so this bound is complete), respecting the structural rule (which child edges a
node admits — `(component, p)` with a cap on derivative-edge multiplicity = the
nonlinearity's degree in that slot — what makes "g quadratic in ∂u" terminate).

Children are selected by a **budget-pruned DFS**: edge atoms are sorted by their
homogeneity contribution, so once a *positive-contribution* atom overshoots the
budget all later ones do too (``break``).  Non-positive contributions are never
pruned on intermediate sums — a partial sum above the bound can still be pulled
back under by further capped negative children.  Termination: an uncapped edge
(a field slot, ``p = 0``) contributes ``> 0`` to the standard part by
subcriticality (``order + β₀ > 0``), so it is bounded by the budget; every edge
that can contribute ``≤ 0`` is a derivative slot with a finite cap.
"""
from __future__ import annotations

import math
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
_MAX_POOL = 5000              # runaway-growth backstop (largest legitimate case ≈ 191 trees)


def _max_node_scaled(sig: Signature, bound: Fraction) -> int:
    """Largest node-decoration weight ``|n|_𝔰`` that can appear in the pool.

    Every rule-conforming tree has ``std ≥ β₀_min``: by subcriticality each edge
    contributes ``edge_gain + child_std ≥ edge_gain + β₀_min > 0`` (induction), so a
    root decoration ``n`` can only enter the pool when ``|n|_𝔰 + β₀_min < bound``.
    At order 2 (``β₀ > −2``, ``bound = 1``) this reproduces the old fixed cap of 2;
    for higher operator order it grows with ``−β₀`` (e.g. cap 3 at order 4, β₀ = −3,
    admitting the counterterm ``Ξ^{(0,3)}`` the fixed cap silently dropped)."""
    beta_min = min(sig.node_homogeneity(b).std for b in sig.node_types)
    return max(0, math.ceil(bound - beta_min) - 1)


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
    decs = list(_multiindices(sig.width, sig.scaling, _max_node_scaled(sig, bound)))
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
                changed |= _emit(b, n, base_h, atoms, caps, add, bound,
                                 sig.grad_budget.get(b))

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
    grad_budget: int | None = None,
) -> bool:
    """DFS over child multisets with homogeneity budget; returns whether pool grew.

    ``grad_budget`` caps the TOTAL number of gradient (``p≠0``) children of the node by the
    nonlinearity's total degree in ∂u (Assumption D2 ⇒ ≤2) — the paper's rule admits at most
    ``I_{e_i}, I_{e_j}`` gradient edges (tex 5337-5340).  ``None`` ⇒ unbounded."""
    grew = False
    counts: Counter[tuple[int, MultiIndex]] = Counter()

    def dfs(
        start: int,
        chosen: list[tuple[int, MultiIndex, DecoratedTree]],
        h: Homogeneity,
        grad_used: int,
    ) -> None:
        nonlocal grew
        if add(tree(b, n, chosen)):
            grew = True
        for idx in range(start, len(atoms)):
            comp, p, sub, contrib = atoms[idx]
            if caps[(comp, p)] is not None and counts[(comp, p)] >= caps[(comp, p)]:
                continue
            is_grad = any(p)
            if is_grad and grad_budget is not None and grad_used >= grad_budget:
                continue                  # D2: ≤ grad_budget gradient edges per node
            nh = h + contrib
            if nh.std >= bound and contrib.std > 0:
                break                     # sorted ascending ⇒ the rest overshoot too
            # contrib.std ≤ 0: recurse even above the bound — further copies of a capped
            # negative slot can pull the sum back under (e.g. ●^{(0,2)}I'(Ξ)² at β₀=−7/4);
            # `add` re-checks the bound, and negative slots are finite-capped, so this
            # terminates.
            counts[(comp, p)] += 1
            dfs(idx, chosen + [(comp, p, sub)], nh, grad_used + int(is_grad))
            counts[(comp, p)] -= 1

    dfs(0, [], base_h, 0)
    return grew
