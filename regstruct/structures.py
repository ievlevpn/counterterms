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
from typing import TYPE_CHECKING, Callable, Hashable

import sympy

from .core.hopf import antipode, convolve, counit
from .equation.dsl import build_context
from .equation.generate import generate_counterterms, generate_trees
from .trees.coproducts import (
    _delta_group_forest, delta_minus, delta_minus_group, delta_plus, twisted_antipode)

if TYPE_CHECKING:
    from .core.signature import Signature
    from .equation.dsl import SPDE
    from .renorm.scheme import NoiseLaw
    from .trees.tree import DecoratedTree


@dataclass
class RenormalizationStructure:
    sig: Signature
    divergent: tuple                       # the SC⁻ generators (|τ|' < 0)
    _h: dict = field(default_factory=dict)  # tree -> its h-symbol

    def coaction(self, t: DecoratedTree) -> dict:   # δ : U → U⁻ ⊗ U
        return delta_minus(t, self.sig)

    def coproduct(self, t: DecoratedTree) -> dict:  # δ⁻ : U⁻ → U⁻ ⊗ U⁻
        return delta_minus_group(t, self.sig)

    def twisted_antipode(self, t: DecoratedTree) -> dict:  # S'₋ : U⁻ → ℝ[U]
        return twisted_antipode(t, self.sig)

    def h_symbol(self, tree: DecoratedTree) -> sympy.Symbol:
        if tree not in self._h:
            self._h[tree] = sympy.Symbol(f"h{len(self._h)}")
        return self._h[tree]

    def bhz_character(self, t: DecoratedTree) -> sympy.Expr:
        """k(τ) = h(S'₋(τ)), a symbolic expression in the (multiplicative) h-values."""
        expr = sympy.Integer(0)
        for forest, coeff in self.twisted_antipode(t).items():
            term = sympy.Rational(coeff.numerator, coeff.denominator)
            for tr in forest:
                term *= self.h_symbol(tr)
            expr += term
        return sympy.expand(expr)

    def canonical_character(self, t: DecoratedTree, law: NoiseLaw | None = None) -> sympy.Expr:
        """The canonical (Gaussian) BHZ constant ``k(τ) = h(S'₋τ)`` with the parity rule
        of a centered Gaussian noise applied: ``h(σ)=0`` whenever σ has an odd number of
        noise vertices (mean zero).  Returns the parity-reduced symbolic combination in
        the surviving (even-noise) ``h``-symbols — many trees collapse to ``0`` (their
        canonical renormalisation constant vanishes).  The surviving symbols' explicit
        Wick-pairing integrals are `renorm.scheme.expectation`."""
        from .renorm.scheme import has_odd_noise
        expr = sympy.Integer(0)
        for forest, coeff in self.twisted_antipode(t).items():
            term = sympy.Rational(coeff.numerator, coeff.denominator)
            for tr in forest:
                if has_odd_noise(tr, self.sig):
                    term = sympy.Integer(0)
                    break
                term *= self.h_symbol(tr)
            expr += term
        return sympy.expand(expr)


def build_renormalization(spde: SPDE) -> RenormalizationStructure:
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

    sig: Signature
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

    def recentering(self, t: DecoratedTree) -> dict:
        """``Δ τ ∈ T ⊗ T⁺`` (tourist_guide.tex 5613)."""
        return delta_plus(t, self.sig)

    def structure_coproduct(self, b: DecoratedTree) -> dict:
        """``Δ⁺ b ∈ T⁺ ⊗ T⁺`` (the structure-group Hopf algebra, tex 5709)."""
        return delta_plus(b, self.sig, project_left=True)

    def structure_antipode(self) -> Callable[[DecoratedTree], dict]:
        """The structure-group antipode ``S : T⁺ → T⁺`` (the convolution inverse of the
        identity), via the generic connected-graded recursion in `core.hopf`."""
        from .core.hopf import antipode
        from .trees.tree import tree
        unit = tree("bullet", (0,) * self.sig.width, (), color="blue")

        def mul(a: DecoratedTree, b: DecoratedTree) -> DecoratedTree:  # X^j·X^k · (merge branches)
            nd = tuple(x + y for x, y in zip(a.node_dec, b.node_dec))
            return tree("bullet", nd, a.children + b.children, color="blue")

        return antipode(self.structure_coproduct, mul, unit)

    def positive_basis(self) -> set:
        """The ``T⁺`` elements that appear as right factors of ``Δ`` over the basis."""
        return {right for t in self.model_basis for (_l, right) in self.recentering(t)}


def build_regularity_structure(spde: SPDE, gamma: Fraction = Fraction(1)) -> RegularityStructure:
    """Build the γ-truncated ``(T, T⁺)`` for `spde` (``gamma`` = the std homogeneity cutoff)."""
    sig, _base, _unknowns = build_context(spde)
    basis = tuple(generate_trees(sig, gamma))
    return RegularityStructure(sig, Fraction(gamma), basis)


@dataclass
class RenormalizationGroup:
    """The renormalization group ``G⁻`` of an SPDE — the character group of the negative
    Hopf algebra ``U⁻``.

    An element is a **character** ``k : U⁻ → ℝ`` (an algebra morphism), determined by its
    free values ``c_τ`` on the negative-tree `generators` and extended multiplicatively
    over forests.  The group law is convolution ``(f⋆g) = (f⊗g)δ⁻`` (`core.hopf.convolve`),
    the unit is the counit ``ε⁻``, and the inverse is ``f ↦ f∘S⁻`` with ``S⁻`` the antipode
    of ``U⁻``.  The renormalized **family** of the equation is the orbit under this group;
    its canonical (BPHZ) element is `RenormalizationStructure.bhz_character`.

    The admissible subgroup ``G⁻_ad ⊂ G⁻`` is a K-admissibility (analytic model) notion —
    it depends on the kernel's vanishing moments and the Π-map — so it is deferred to the
    model layer (Phase 4); a symbolic-only reduction would be unsound.
    """

    sig: Signature
    generators: tuple                       # the negative-tree free parameters c_τ

    def _coproduct(self, forest: tuple[DecoratedTree, ...]) -> dict:
        return _delta_group_forest(forest, self.sig)

    @staticmethod
    def _mul(a: tuple[DecoratedTree, ...], b: tuple[DecoratedTree, ...]) -> tuple[DecoratedTree, ...]:
        # the algebra product on U⁻ (forest concat)
        return tuple(sorted(a + b, key=lambda t: t._sortkey()))

    def identity(self) -> Callable[[Hashable], object]:
        """The neutral element ``ε⁻`` (no renormalisation)."""
        return counit(())

    def character(self, values: dict) -> Callable[[tuple[DecoratedTree, ...]], object]:
        """The character with free constants `values` (``tree → scalar``), extended as an
        algebra morphism ``U⁻ → ℝ`` (a product over the forest, ``1`` on the unit)."""
        def chi(forest: tuple[DecoratedTree, ...]) -> object:
            out = 1
            for t in forest:
                out = out * values.get(t, 0)
            return out
        return chi

    def convolve(self, f: Callable, g: Callable) -> Callable[[Hashable], object]:
        """The group product ``(f⋆g)(τ) = (f⊗g)δ⁻τ``."""
        return convolve(f, g, self._coproduct)

    def inverse(self, f: Callable) -> Callable[[tuple[DecoratedTree, ...]], object]:
        """The group inverse ``f⋆⁻¹ = f∘S⁻`` (valid for characters)."""
        S = antipode(self._coproduct, self._mul, ())
        return lambda forest: sum((c * f(g) for g, c in S(forest).items()), 0)

    def admissible(self) -> object:
        """[SOCKET — Track B3] The admissible subgroup ``G⁻_ad ⊂ G⁻``.

        K-admissibility (``Π(I_nτ)=∂^nK∗Πτ``) constrains which renormalisation characters
        a model admits; it is an *analytic model* notion (kernel vanishing moments + the
        Π-map), not symbolic. Deferred to the model layer — see notes/phase4_plan.md B3.
        """
        raise NotImplementedError(
            "G⁻_ad reduction is Phase 4 / Track B3 — admissibility is a K-admissibility "
            "(model) notion (kernel vanishing moments + the Π-map), not a symbolic filter; "
            "see notes/phase4_plan.md B3.")


def build_renormalization_group(spde: SPDE) -> RenormalizationGroup:
    """Build the renormalization group ``G⁻`` for `spde` (free constant per negative tree)."""
    sig, _base, _unknowns = build_context(spde)
    return RenormalizationGroup(sig, tuple(generate_counterterms(sig)))
