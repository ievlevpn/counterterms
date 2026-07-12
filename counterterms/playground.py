"""Standalone tree/algebra playground — the Hopf-algebra layer without an SPDE.

``playground(...)`` builds a minimal `Signature` (noise regularities + parabolic
scaling only; the structural rule is empty, so tree *generation* is not available —
that needs a real `SPDE`).  Hand-build trees with ``.Xi`` / ``.I(...)`` /
``.tree(...)``, and homogeneity, S(τ), δ⁻, δ⁺, the twisted antipode, the BHZ
character and the Wick expectation become one-liners:

    pg = playground("-3/2 - kappa")           # KPZ noise, d = 1
    t = pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"))    # ●·𝓘ₓ[Ξ]²
    t.homogeneity(pg.sig), t.symmetry_factor(), pg.character(t), pg.expectation(t)
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from .core.homogeneity import Homogeneity, MultiIndex, Scaling
from .core.signature import Signature
from .equation.dsl import _split_kappa
from .render.tree import coord_names
from .trees.tree import DecoratedTree, Edge, tree

if TYPE_CHECKING:
    import sympy

    from .renorm.scheme import Expectation
    from .structures import RenormalizationStructure
    from .trees.coproducts import ForestSum, PlusSum, TensorSum


@dataclass
class Playground:
    """An algebra sandbox around a hand-made `Signature` — see the module docstring."""

    sig: Signature
    _rs: RenormalizationStructure | None = field(default=None, repr=False)

    # --- building trees -------------------------------------------------------- #

    def _mi(self, spec: object) -> MultiIndex:
        """A multi-index from ``"x"``/``"xx"``/``"tx"``/``""`` (coordinate letters — t, x,
        y, z) or a ready tuple (padded to width)."""
        if isinstance(spec, tuple):
            return spec + (0,) * (self.sig.width - len(spec))
        coords = coord_names(self.sig.width)
        out = [0] * self.sig.width
        for ch in str(spec):
            out[coords.index(ch)] += 1
        return tuple(out)

    def noise(self, name: str | None = None) -> DecoratedTree:
        """The noise leaf ``Ξ`` (``name`` may be omitted when there is a single noise)."""
        names = list(self.sig.noise_homog)
        if name is None:
            if len(names) != 1:
                raise ValueError(f"several noises {names}; pass a name")
            name = names[0]
        if name not in self.sig.noise_homog:
            raise ValueError(f"unknown noise {name!r} (have {names})")
        return tree(name, (0,) * self.sig.width)

    @property
    def Xi(self) -> DecoratedTree:
        """Sugar for the single-noise leaf, ``pg.noise()``."""
        return self.noise()

    def I(self, sub: DecoratedTree, p: object = "", comp: int = 0) -> Edge:  # noqa: E743 — the paper's 𝓘
        """An integration edge ``𝓘^{(comp)}_p[sub]`` — feed to `tree`."""
        return (comp, self._mi(p), sub)

    def tree(self, *edges: Edge, X: object = "") -> DecoratedTree:
        """A ``●`` node with child ``edges`` (from `I`) and node decoration ``X``
        (e.g. ``X="x"`` for ``Xₓ·``); canonicalized."""
        return tree("bullet", self._mi(X), tuple(edges))

    # --- the algebra, one-liners ------------------------------------------------ #

    def structure(self) -> RenormalizationStructure:
        """The `RenormalizationStructure` over this signature (no generated trees —
        `divergent` is empty; the coproducts/antipode/character work on any tree)."""
        if self._rs is None:
            from .structures import RenormalizationStructure
            self._rs = RenormalizationStructure(self.sig, ())
        return self._rs

    def homogeneity(self, t: DecoratedTree) -> Homogeneity:
        return t.homogeneity(self.sig)

    def delta_minus(self, t: DecoratedTree) -> TensorSum:
        """The coaction δ⁻ : U → U⁻ ⊗ U."""
        from .trees.coproducts import delta_minus
        return delta_minus(t, self.sig)

    def delta_plus(self, t: DecoratedTree) -> PlusSum:
        """The T⁺ coproduct Δ⁺."""
        from .trees.coproducts import delta_plus
        return delta_plus(t, self.sig)

    def antipode(self, t: DecoratedTree) -> ForestSum:
        """The negative twisted antipode S'₋ : U⁻ → ℝ[U]."""
        from .trees.coproducts import twisted_antipode
        return twisted_antipode(t, self.sig)

    def character(self, t: DecoratedTree) -> sympy.Expr:
        """The canonical (BPHZ, centered-Gaussian) constant ``k(τ) = h(S'₋τ)``, symbolic
        in the elementary expectations ``h``; the h-symbols are stable across calls."""
        return self.structure().canonical_character(t)

    def expectation(self, t: DecoratedTree) -> Expectation:
        """The Wick-pairing integral ``h_ε(τ) = 𝔼[Π^ε τ](0)`` (bare trees; unevaluated)."""
        from .renorm.scheme import WHITE_NOISE, expectation
        return expectation(t, self.sig, WHITE_NOISE)


def playground(regularity: object = "-3/2 - kappa", dim: int = 1, order: int = 2,
               noises: dict[str, object] | None = None) -> Playground:
    """One-call algebra sandbox: a single noise ``xi`` of the given ``regularity`` (or
    several via ``noises={name: regularity}``), parabolic ``(order, 1, …, 1)`` scaling in
    ``dim`` space dimensions.  β₀-conventions as everywhere: |Ξ| = the regularity itself."""
    spec = noises if noises is not None else {"xi": regularity}
    homog = {n: Homogeneity(*_split_kappa(r)) for n, r in spec.items()}
    node_types = ("bullet", *homog)
    sig = Signature(
        dim=dim, scaling=Scaling(tuple([order] + [1] * dim)), n_components=1,
        comp_order=(order,), noise_homog=homog, node_types=node_types,
        allowed={b: () for b in node_types})     # empty rule: hand-built trees only
    return Playground(sig)
