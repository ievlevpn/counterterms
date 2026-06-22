"""Assembled renormalization structure for an SPDE (Phase 3).

Bundles the validated coproducts (`δ`, `δ⁻`), the negative twisted antipode `S'₋`,
and the **symbolic BHZ character** `k = h∘S'₋` for a given equation.  The analytic
input `h(τ) = 𝔼[Π^ζ τ](0)` is left as a free symbol per tree — its numeric value
needs Gaussian/Wick machinery (Phase 4); here `k(τ)` is returned as the exact
symbolic combination of those `h`-values that the twisted antipode prescribes.

Note (tex 5710/5717): `δ⁻` and `Δ` are individually coassociative and the
cointeraction `(Id⊗Δ)δ = M¹³(δ⊗δ⁺)Δ` holds, including for the singular β₀=−3/2
trees (the e'-recentering runs over the full ∂(A,F); see `tests/test_coproducts.py`
and `notes/cointeraction_singular_noise.md`).
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from fractions import Fraction

import sympy

from .equation.dsl import build_context
from .equation.generate import generate_counterterms, generate_trees
from .trees.coproducts import (
    delta_minus, delta_minus_group, delta_plus, twisted_antipode)


@dataclass
class RenormalizationStructure:
    sig: object
    divergent: tuple                       # the SC⁻ generators (|τ|' < 0)
    _h: dict = field(default_factory=dict)  # tree -> its h-symbol

    def coaction(self, t):                 # δ : U → U⁻ ⊗ U
        return delta_minus(t, self.sig)

    def coproduct(self, t):                # δ⁻ : U⁻ → U⁻ ⊗ U⁻
        return delta_minus_group(t, self.sig)

    def twisted_antipode(self, t):         # S'₋ : U⁻ → ℝ[U]
        return twisted_antipode(t, self.sig)

    def h_symbol(self, tree) -> sympy.Symbol:
        if tree not in self._h:
            self._h[tree] = sympy.Symbol(f"h{len(self._h)}")
        return self._h[tree]

    def bhz_character(self, t):
        """k(τ) = h(S'₋(τ)), a symbolic expression in the (multiplicative) h-values."""
        expr = sympy.Integer(0)
        for forest, coeff in self.twisted_antipode(t).items():
            term = sympy.Rational(coeff.numerator, coeff.denominator)
            for tr in forest:
                term *= self.h_symbol(tr)
            expr += term
        return sympy.expand(expr)


def build_renormalization(spde) -> RenormalizationStructure:
    sig, _base, _unknowns = build_context(spde)
    return RenormalizationStructure(sig, tuple(generate_counterterms(sig)))


@dataclass
class RegularityStructure:
    """A γ-truncated regularity structure ``(T, T⁺)`` for an SPDE (the positive side).

    ``model_basis`` is the basis of the model space ``T``: every rule-conforming tree
    with ``|τ|_std < γ``, graded by homogeneity (``A`` = the homogeneity set, locally
    finite and bounded below).  The recentering ``Δ : T → T ⊗ T⁺`` and the
    structure-group coproduct ``Δ⁺ : T⁺ → T⁺ ⊗ T⁺`` are the validated `delta_plus`;
    ``T⁺`` is spanned by the positive planted trees that ``Δ`` produces.  The negative
    (renormalisation) side is `RenormalizationStructure`; the two cointeract.

    ponytail: the abstract polynomial sector ``{X^k}`` is carried as node decorations,
    not as standalone basis vectors (there is no separate ``X`` node type).
    """

    sig: object
    gamma: Fraction
    model_basis: tuple

    def grades(self) -> dict:
        """``T = ⊕_α T_α`` — basis trees grouped by homogeneity."""
        g = defaultdict(list)
        for t in self.model_basis:
            g[t.homogeneity(self.sig)].append(t)
        return dict(g)

    def homogeneities(self) -> list:
        """The homogeneity set ``A`` (sorted, ascending)."""
        return sorted({t.homogeneity(self.sig) for t in self.model_basis},
                      key=lambda h: h._key())

    @property
    def divergent(self) -> tuple:
        """The negative-homogeneity subspace (the counterterm carriers)."""
        return tuple(t for t in self.model_basis if t.homogeneity(self.sig).is_negative())

    def recentering(self, t):
        """``Δ τ ∈ T ⊗ T⁺`` (tourist_guide.tex 5613)."""
        return delta_plus(t, self.sig)

    def structure_coproduct(self, b):
        """``Δ⁺ b ∈ T⁺ ⊗ T⁺`` (the structure-group Hopf algebra, tex 5709)."""
        return delta_plus(b, self.sig, project_left=True)

    def structure_antipode(self):
        """The structure-group antipode ``S : T⁺ → T⁺`` (the convolution inverse of the
        identity), via the generic connected-graded recursion in `core.hopf`."""
        from .core.hopf import antipode
        from .trees.tree import tree
        unit = tree("bullet", (0,) * self.sig.width, (), color="blue")

        def mul(a, b):                                   # X^j·X^k · (merge branches)
            nd = tuple(x + y for x, y in zip(a.node_dec, b.node_dec))
            return tree("bullet", nd, a.children + b.children, color="blue")

        return antipode(self.structure_coproduct, mul, unit)

    def positive_basis(self) -> set:
        """The ``T⁺`` elements that appear as right factors of ``Δ`` over the basis."""
        return {right for t in self.model_basis for (_l, right) in self.recentering(t)}


def build_regularity_structure(spde, gamma=Fraction(1)) -> RegularityStructure:
    """Build the γ-truncated ``(T, T⁺)`` for `spde` (``gamma`` = the std homogeneity cutoff)."""
    sig, _base, _unknowns = build_context(spde)
    basis = tuple(generate_trees(sig, gamma))
    return RegularityStructure(sig, Fraction(gamma), basis)
