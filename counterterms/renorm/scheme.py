"""Phase 4 / Track B (B1) — the canonical (BPHZ) character, symbolic half.

A centered Gaussian noise law turns the elementary expectation ``h(σ)=𝔼[Π^ζσ](0)``
into a **Wick-pairing sum** (Isserlis):

* ``0`` when σ has an *odd* number of noise vertices (mean zero) — so the corresponding
  canonical constant vanishes;
* otherwise ``Σ`` over perfect matchings of the noise vertices of an *unevaluated*
  integral — kernel factors ``K^{(p)}`` on the tree edges times covariance factors ``C``
  on the matched pairs, over the internal vertices (root at 0).

The integrals are the genuinely divergent objects — evaluating them is Track B2; here
they stay symbolic.  Composed with the twisted antipode ``S'₋`` (which we have), the
parity rule alone already determines which canonical renormalisation constants are zero.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Protocol

import sympy

from ..trees.coproducts import _explode

if TYPE_CHECKING:
    from ..core.signature import Signature
    from ..structures import RenormalizationStructure
    from ..trees.tree import DecoratedTree


@dataclass(frozen=True)
class NoiseLaw:
    """A centered Gaussian noise law — the covariance label for the symbolic kernel
    (``𝔼[ξ(x)ξ(y)] = C(x−y)``; white noise ⇒ ``C = δ``, collapsed only at evaluation)."""

    covariance: str = "C"
    name: str = "white noise"


WHITE_NOISE = NoiseLaw()


def wick_pairings(items: Iterable[int]) -> list[list[tuple[int, int]]]:
    """All perfect matchings of ``items`` (Isserlis' theorem); ``[]`` if the count is odd."""
    items = list(items)
    if len(items) % 2 == 1:
        return []
    if not items:
        return [[]]
    first, rest = items[0], items[1:]
    out = []
    for i, other in enumerate(rest):
        for sub in wick_pairings(rest[:i] + rest[i + 1:]):
            out.append([(first, other)] + sub)
    return out


def has_odd_noise(tree: DecoratedTree, sig: Signature) -> bool:
    """Does *any* noise type appear an odd number of times? (⇒ canonical expectation 0.)

    For independent centered noises the expectation factors over noise types, so a single
    type with odd count already kills it — e.g. ``Ξ_ξ·I[Ξ_η]`` (one ξ, one η) vanishes even
    though the total count (2) is even."""
    from collections import Counter
    nodes, _ = _explode(tree)
    counts = Counter(n.node_type for n in nodes if n.node_type in sig.noise_homog)
    return any(c % 2 == 1 for c in counts.values())


@dataclass
class Expectation:
    """``𝔼[Π^ζτ](0)`` as a sum of unevaluated Wick-pairing integrals (empty list ⇒ 0).
    Each term is ``(integrand, integration_variables)``."""

    terms: list[tuple[sympy.Expr, tuple[sympy.Symbol, ...]]]

    @property
    def is_zero(self) -> bool:
        return not self.terms

    def __str__(self) -> str:
        if self.is_zero:
            return "0"
        return " + ".join(f"∫ {ig} d{vs}" for ig, vs in self.terms)


def is_extended(tree: DecoratedTree) -> bool:
    """Does the tree carry a red (contraction) node — an extended ``o``-decoration?"""
    return tree.color == "red" or any(is_extended(s) for (_c, _p, s) in tree.children)


def is_bare(tree: DecoratedTree) -> bool:
    """Is the tree *bare* — **every node decoration is zero**?  This (not the presence of red
    contraction nodes) is exactly the domain where the naive ``Π^ζ`` reduces to kernels⊗noises:
    ``Π^ζ(●^{n,α})(x)=x^n`` is *independent of the extended decoration* ``α`` (tex 4003–4004), so
    a red node with ``n=0`` is just the ``1``-vertex `expectation` already builds.  Only a
    non-zero ``n`` (the ``x^n`` factor) breaks the integral — see `expectation`."""
    return (not any(tree.node_dec)) \
        and all(is_bare(s) for (_c, _p, s) in tree.children)


def expectation_key(tree: DecoratedTree) -> tuple:
    """Identity of ``h(τ)=𝔼[Π^ζτ](0)`` **ignoring the o/colour decoration**.  ``Π^ζ(●^{n,α})``
    is independent of the extended decoration ``α`` (tex 4003–4004), so trees that differ only
    in their red nodes' o-decorations have the *same* elementary expectation.  Used only to
    *flag* duplicate ``h``-symbols in the report (display-only); the symbolic character that
    `structures.canonical_character` returns is unchanged."""
    return (tree.node_type, tree.node_dec,
            tuple(sorted((c, p, expectation_key(s)) for (c, p, s) in tree.children)))


def zero_note(tree: DecoratedTree, sig: Signature) -> str | None:
    """A short reason iff ``h(τ)=𝔼[Π^ζτ](0)`` is **provably 0** (for display-only marking),
    else ``None``.  Three rigorous cases, needing no kernel evaluation:

    * a non-zero **root** decoration ``X^n`` ⇒ ``Π(X^nτ)(0)=0^n=0`` (tex 1809, 5083);
    * an **odd** count of some noise type ⇒ 𝔼 of an odd number of centered Gaussians is 0;
    * a **no-noise** tree with a *derivative-edge leaf* (``p≠0``) ⇒ that leaf integrates a total
      derivative ``∫∂^pK=0`` (``p≠0``, decaying/compactly-supported ``K``), so the whole
      product vanishes.  (A no-noise tree with only ``p=0`` leaves needs the kernel's order-0
      moment condition, so it is *not* flagged here — conservative.)"""
    if any(tree.node_dec):
        return "root X^n vanishes at the base point"
    if has_odd_noise(tree, sig):
        return "odd noise parity"
    nodes, edges = _explode(tree)
    if not any(n.node_type in sig.noise_homog for n in nodes):
        parents = {a for (a, _b, _c, _p) in edges}
        for (_a, b, _c, p) in edges:
            if b not in parents and any(p):          # derivative edge into a leaf
                return "pure-kernel total derivative"
    return None


def expectation(tree: DecoratedTree, sig: Signature, law: NoiseLaw = WHITE_NOISE) -> Expectation:
    """The Wick expansion of ``h(τ)=𝔼[Π^ζτ](0)`` (symbolic; integrals unevaluated).

    **Bare trees only** (every node decoration zero).  The canonical model ``Π^ζ`` (tex 4000–4008)
    is multiplicative with ``Π^ζ(∘^n)(x)=x^n ζ(x)``, ``Π^ζ(●^n)(x)=Π^ζ(●^{n,α})(x)=x^n`` and
    ``Π^ζ(I_pτ)=∂^pK(Π^ζτ)``.  For a *bare* tree (all ``n=0``) this collapses to the kernel⊗noise
    integral built below, and ``𝔼[·](0)`` is exactly the Wick-pairing sum.  **Red contraction
    nodes are fine**: ``Π^ζ(●^{n,α})`` is independent of the extended decoration ``α`` (tex 4003),
    so a red node with ``n=0`` is the ``1``-vertex this builds — *not* something needing the full
    model.  The **only** thing that breaks the integral is a **non-zero node decoration** ``n``:
    ``Π^ζ`` then carries ``x^n`` — at the root (``x=0``) that is ``0^n`` (so ``h=0``), at an
    internal vertex a polynomial factor ``z^n``; neither is captured here, so we refuse.

    Independent centered noises ⇒ the expectation **factors over noise types**: vertices are
    paired only with same-type vertices (no cross-noise covariances, since ``𝔼[ξη]=0``), and
    any type with an odd count makes the whole expectation 0.

    Order: the parity rule (𝔼 of an odd number of centered Gaussians is 0, whatever the
    deterministic kernel/``x^n`` factors, and red nodes carry no noise) is always valid, so it
    returns 0 first; only a non-zero, *bare* tree reaches the integral."""
    from collections import defaultdict
    from itertools import product

    nodes, edges = _explode(tree)
    by_type: dict[str, list[int]] = defaultdict(list)
    for n in nodes:
        if n.node_type in sig.noise_homog:
            by_type[n.node_type].append(n.id)
    if any(len(ids) % 2 == 1 for ids in by_type.values()):
        return Expectation([])                      # odd count ⇒ E = 0 (always valid)

    if not is_bare(tree):                            # non-zero node decoration ⇒ value ≠ this integral
        raise NotImplementedError(
            "expectation() is the naive Wick integral, valid only for bare trees (all node "
            "decorations 0); this σ has a non-zero X^n node-decoration. By Π(X^n)(y)=y^n "
            "(tex 1809, 4003) a root X^n gives a 0^n factor (h=0) and internal X^n give "
            "polynomial factors — the full canonical model (Track B). (Red contraction nodes "
            "are fine: Π^ζ(●^{n,α}) is α-independent.)")

    z = {n.id: (sympy.Integer(0) if n.id == 0 else sympy.Symbol(f"z{n.id}")) for n in nodes}
    intvars = tuple(z[n.id] for n in nodes if n.id != 0)
    kernels = sympy.Integer(1)
    for (a, b, _comp, p) in edges:                  # K^{(p)}(z_parent − z_child) per edge
        kernels *= sympy.Function("K_" + "_".join(map(str, p)))(z[a] - z[b])
    C = sympy.Function(law.covariance)

    # a global pairing = one within-type matching per noise type (Cartesian product)
    terms = []
    for combo in product(*(wick_pairings(ids) for ids in by_type.values())):
        cov = sympy.Integer(1)
        for matching in combo:
            for (i, j) in matching:
                cov *= C(z[i] - z[j])
        terms.append((sympy.expand(kernels * cov), intvars))
    return Expectation(terms)


# --------------------------------------------------------------------------- #
# the renormalization-scheme seam (architecture §3.10)
#
# The renormalized *family* is, by definition, the orbit of the equation under a choice
# of renormalization character.  A `RenormalizationScheme` makes that choice: it turns a
# `RenormalizationStructure` into a `Character` (the constants c_τ).  Two schemes:
#   FreeConstants — the default: each c_τ a free symbol (what `renormalize` emits today).
#   BPHZ          — the canonical k^ζ = h^ζ∘S'₋: symbolic (parity-reduced) now; the
#                   *numeric* values need an `IntegralEvaluator` (Track B2 — not built).
# --------------------------------------------------------------------------- #

@dataclass
class Character:
    """A choice of renormalization constants — ``τ ↦ value`` (a free symbol, a symbolic
    expression in the elementary expectations, or — once Track B2 exists — a number)."""

    values: dict[DecoratedTree, object]

    def __call__(self, tree: DecoratedTree) -> object:
        return self.values[tree]


class RenormalizationScheme(Protocol):
    """Turn a built `RenormalizationStructure` into a `Character`."""

    def character(self, structure: RenormalizationStructure) -> Character:
        ...


class FreeConstants:
    """The default scheme: each ``c_τ`` is a free symbol (the safe, complete answer — a
    solution *is* the family indexed by these)."""

    def character(self, structure: RenormalizationStructure) -> Character:
        return Character({t: sympy.Symbol(f"c_{i}")
                          for i, t in enumerate(structure.divergent)})


class BPHZ:
    """The canonical scheme ``k^ζ = h^ζ ∘ S'₋``.

    `character` gives the **symbolic, parity-reduced** constant per tree (via
    `RenormalizationStructure.canonical_character`) — odd-noise trees vanish.  The
    surviving constants are polynomials in the elementary expectations ``h(σ)``; turning
    those into *numbers* needs an `IntegralEvaluator` (Track B2).  [SOCKET: the numeric
    path is wired but unbuilt — `numeric_character` raises until B2 lands.]
    """

    def __init__(self, noise: NoiseLaw = WHITE_NOISE, evaluator: IntegralEvaluator | None = None) -> None:
        self.noise = noise
        self.evaluator = evaluator

    def character(self, structure: RenormalizationStructure) -> Character:
        return Character({t: structure.canonical_character(t, self.noise)
                          for t in structure.divergent})

    def numeric_character(self, structure: RenormalizationStructure) -> Character:
        """The canonical constants as numbers — requires an `IntegralEvaluator` (B2)."""
        ev = self.evaluator or UnbuiltEvaluator()
        return Character({t: ev.evaluate(expectation(t, structure.sig, self.noise),
                                         noise=self.noise)
                          for t in structure.divergent})


class IntegralEvaluator(Protocol):
    """[SOCKET — Track B2] Evaluate a Wick-pairing `Expectation` to a value, given the
    concrete kernel and covariance.  The singular/divergent integrals (and their
    regularised renormalisation) live behind this seam — the analysis wall."""

    def evaluate(self, expectation: Expectation, *, noise: NoiseLaw) -> object:
        ...


class UnbuiltEvaluator:
    """Placeholder `IntegralEvaluator` — Track B2 is not implemented (the analysis wall)."""

    def evaluate(self, expectation: Expectation, *, noise: NoiseLaw = WHITE_NOISE) -> object:
        raise NotImplementedError(
            "evaluating the divergent Wick-pairing integrals is Phase 4 / Track B2 — it "
            "needs singular-integral machinery (closed forms / regularised numerics); "
            "see notes/phase4_plan.md B2.")

