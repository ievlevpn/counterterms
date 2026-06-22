"""Assembled renormalization structure for an SPDE (Phase 3).

Bundles the validated coproducts (`δ`, `δ⁻`), the negative twisted antipode `S'₋`,
and the **symbolic BHZ character** `k = h∘S'₋` for a given equation.  The analytic
input `h(τ) = 𝔼[Π^ζ τ](0)` is left as a free symbol per tree — its numeric value
needs Gaussian/Wick machinery (Phase 4); here `k(τ)` is returned as the exact
symbolic combination of those `h`-values that the twisted antipode prescribes.

Note (tex 5710/5717): `δ⁻` and `Δ` are individually coassociative and the
cointeraction holds for the gKPZ class; the cointeraction has a known residual
for more singular noise (β₀ ≤ −3/2) — see `tests/test_coproducts.py`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import sympy

from .equation.dsl import build_context
from .equation.generate import generate_counterterms
from .trees.coproducts import delta_minus, delta_minus_group, twisted_antipode


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
