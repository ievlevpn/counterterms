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
from typing import Protocol

import sympy

from ..trees.coproducts import _explode


@dataclass(frozen=True)
class NoiseLaw:
    """A centered Gaussian noise law — the covariance label for the symbolic kernel
    (``𝔼[ξ(x)ξ(y)] = C(x−y)``; white noise ⇒ ``C = δ``, collapsed only at evaluation)."""

    covariance: str = "C"
    name: str = "white noise"


WHITE_NOISE = NoiseLaw()


def wick_pairings(items):
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


def _noise_vertices(tree, sig):
    nodes, _ = _explode(tree)
    return [n.id for n in nodes if n.node_type in sig.noise_homog]


def has_odd_noise(tree, sig) -> bool:
    """Does the tree carry an odd number of noise vertices? (⇒ canonical expectation 0.)"""
    return len(_noise_vertices(tree, sig)) % 2 == 1


@dataclass
class Expectation:
    """``𝔼[Π^ζτ](0)`` as a sum of unevaluated Wick-pairing integrals (empty list ⇒ 0).
    Each term is ``(integrand, integration_variables)``."""

    terms: list

    @property
    def is_zero(self) -> bool:
        return not self.terms

    def __str__(self) -> str:
        if self.is_zero:
            return "0"
        return " + ".join(f"∫ {ig} d{vs}" for ig, vs in self.terms)


def expectation(tree, sig, law: NoiseLaw = WHITE_NOISE) -> Expectation:
    """The Wick expansion of ``h(τ)=𝔼[Π^ζτ](0)`` (symbolic; integrals unevaluated)."""
    nodes, edges = _explode(tree)
    noise_ids = [n.id for n in nodes if n.node_type in sig.noise_homog]
    if len(noise_ids) % 2 == 1:
        return Expectation([])                      # mean-zero ⇒ E = 0

    z = {n.id: (sympy.Integer(0) if n.id == 0 else sympy.Symbol(f"z{n.id}")) for n in nodes}
    intvars = tuple(z[n.id] for n in nodes if n.id != 0)
    kernels = sympy.Integer(1)
    for (a, b, _comp, p) in edges:                  # K^{(p)}(z_parent − z_child) per edge
        kernels *= sympy.Function("K_" + "_".join(map(str, p)))(z[a] - z[b])
    C = sympy.Function(law.covariance)

    terms = []
    for pairing in wick_pairings(noise_ids):
        cov = sympy.Integer(1)
        for (i, j) in pairing:
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

    values: dict

    def __call__(self, tree):
        return self.values[tree]


class RenormalizationScheme(Protocol):
    """Turn a built `RenormalizationStructure` into a `Character`."""

    def character(self, structure) -> Character:
        ...


class FreeConstants:
    """The default scheme: each ``c_τ`` is a free symbol (the safe, complete answer — a
    solution *is* the family indexed by these)."""

    def character(self, structure) -> Character:
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

    def __init__(self, noise: NoiseLaw = WHITE_NOISE, evaluator=None):
        self.noise = noise
        self.evaluator = evaluator

    def character(self, structure) -> Character:
        return Character({t: structure.canonical_character(t, self.noise)
                          for t in structure.divergent})

    def numeric_character(self, structure) -> Character:
        """The canonical constants as numbers — requires an `IntegralEvaluator` (B2)."""
        ev = self.evaluator or UnbuiltEvaluator()
        return Character({t: ev.evaluate(expectation(t, structure.sig, self.noise),
                                         noise=self.noise)
                          for t in structure.divergent})


class IntegralEvaluator(Protocol):
    """[SOCKET — Track B2] Evaluate a Wick-pairing `Expectation` to a value, given the
    concrete kernel and covariance.  The singular/divergent integrals (and their
    regularised renormalisation) live behind this seam — the analysis wall."""

    def evaluate(self, expectation: Expectation, *, noise: NoiseLaw):
        ...


class UnbuiltEvaluator:
    """Placeholder `IntegralEvaluator` — Track B2 is not implemented (the analysis wall)."""

    def evaluate(self, expectation: Expectation, *, noise: NoiseLaw = WHITE_NOISE):
        raise NotImplementedError(
            "evaluating the divergent Wick-pairing integrals is Phase 4 / Track B2 — it "
            "needs singular-integral machinery (closed forms / regularised numerics); "
            "see notes/phase4_plan.md B2.")

