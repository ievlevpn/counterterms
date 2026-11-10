"""The output object: the family of renormalized equations.

For a system, each component ``a`` gets its own renormalized equation; the trees
(and their constants ``k_τ``) are **shared** across components — the same `k_τ`
weights `F_a(τ*)` in every equation where it is nonzero (tourist_guide.tex 4915).
A scalar problem is the one-component case.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import sympy

from ..core.jets import is_jet, jet, jet_parts

if TYPE_CHECKING:
    from pathlib import Path

    from ..core.homogeneity import Homogeneity
    from ..core.signature import Signature
    from ..equation.dsl import SPDE, Unknown
    from ..render.display import Picture
    from ..renorm.scheme import Expectation, NoiseLaw
    from ..structures import RenormalizationStructure
    from ..trees.tree import DecoratedTree


@dataclass(frozen=True)
class Counterterm:
    tree: DecoratedTree
    homogeneity: Homogeneity
    symmetry_factor: int
    elem_diff: sympy.Expr     # F_a(τ*) in jet variables
    constant: sympy.Symbol   # the free renormalization constant k_τ (shared across components)

    @property
    def coefficient(self) -> sympy.Expr:
        return self.constant / self.symmetry_factor      # k_τ / S(τ)

    @property
    def term(self) -> sympy.Expr:
        return self.coefficient * self.elem_diff

    @property
    def divergence(self) -> sympy.Expr:
        """Leading ε→0 rate of the canonical value of ``k_τ`` (`scheme.divergence_rate`):
        ``ε^{|τ|₀}``, or ``log(1/ε)`` for a κ-critical tree.  Up to log corrections."""
        from .scheme import divergence_rate
        return divergence_rate(self.homogeneity)

    def __str__(self) -> str:
        from .scheme import divergence_str
        return (f"{sympy.sstr(sympy.expand(self.term))}   "
                f"[τ = {self.tree}, |τ| = {self.homogeneity}, S = {self.symmetry_factor}, "
                f"~ {divergence_str(self.homogeneity)}]")

    def __repr__(self) -> str:
        return self.__str__()

    def _repr_latex_(self) -> str:
        from ..render.tree import shorthand_latex
        hom = str(self.homogeneity).replace("κ", r"\kappa")
        return (f"${sympy.latex(sympy.expand(self.term))} \\qquad "
                f"\\tau = {shorthand_latex(self.tree)},\\ |\\tau| = {hom},"
                f"\\ S = {self.symmetry_factor}$")


def _collected_rhs(eq: RenormalizedEquation, comp: int,
                   canon: dict[DecoratedTree, sympy.Expr]) -> sympy.Expr:
    """The renormalized RHS of component ``comp`` with each free constant replaced by its
    value from ``canon``: original equation + Σ (value/S(τ))·F(τ*); zero values drop out."""
    expr = eq.original_rhs(comp)
    for ct in eq.per_component[comp]:
        v = canon.get(ct.tree, sympy.Integer(0))
        if v != 0:
            expr += sympy.Rational(1, ct.symmetry_factor) * v * eq._pretty(ct.elem_diff)
    return sympy.expand(expr)


@dataclass
class CanonicalFamily:
    """The canonical (BPHZ) specialization ``k_τ = h(S'₋τ)`` of the free-constant family,
    as **computable sympy objects** (the report renders *from* this).

    ``constants[τ]`` is the canonical value of the free constant ``k_τ`` — a polynomial in
    the elementary expectations ``hᵢ = h_ε(σᵢ) = 𝔼[Π^ε σᵢ](0)`` (``0`` when it provably
    vanishes); ``legend[hᵢ]`` is the tree ``σᵢ``.  With ``reduced`` the exact identities of
    `renorm.scheme.reduction_subs` are folded in (`symmetric` additionally applies the
    spatial-reflection identity — valid only for a spatially-symmetric noise)."""

    eq: RenormalizedEquation
    structure: RenormalizationStructure               # the underlying Phase-3 algebra
    reduced: bool
    symmetric: bool
    constants: dict[DecoratedTree, sympy.Expr]        # τ → value of k_τ (homogeneity order)
    free_constants: dict[DecoratedTree, sympy.Symbol]  # τ → the free symbol k_τ it replaces
    legend: dict[sympy.Symbol, DecoratedTree]         # hᵢ → σᵢ

    @property
    def substitution(self) -> dict[sympy.Symbol, sympy.Expr]:
        """``{k_τ: canonical value}`` — apply with ``.xreplace`` to any expression built
        from the free family (e.g. ``eq.counterterm_rhs(0)``)."""
        return {self.free_constants[t]: v for t, v in self.constants.items()}

    def rhs(self, comp: int = 0) -> sympy.Expr:
        """The canonically renormalized RHS of component ``comp`` (collected; vanishing
        counterterms dropped)."""
        return _collected_rhs(self.eq, comp, self.constants)

    def expectations(self, law: NoiseLaw | None = None) -> dict[sympy.Symbol, Expectation]:
        """``hᵢ → `` its Wick-pairing integral (`renorm.scheme.expectation`), for the bare
        σᵢ (the non-bare ones need the full canonical model — Track B — and are skipped)."""
        from .scheme import WHITE_NOISE, expectation, is_bare
        return {sym: expectation(tr, self.eq.sig, law or WHITE_NOISE)
                for sym, tr in self.legend.items() if is_bare(tr)}

    def __str__(self) -> str:
        from ..render.tree import shorthand
        from .scheme import divergence_str
        sig = self.eq.sig                             # name the noises as the family does
        tag = ("canonical (BPHZ) constants"
               + (", reduced" if self.reduced else "")
               + (" [symmetric noise]" if self.reduced and self.symmetric else ""))
        lines = [tag + ":"]
        for t, v in self.constants.items():
            k = self.free_constants[t]
            rate = f", ~ {divergence_str(t.homogeneity(self.eq.sig))}" if v != 0 else ""
            lines.append(f"  {k} = {sympy.sstr(v)}   [τ = {shorthand(t, sig)}{rate}]")
        for sym, tr in self.legend.items():
            lines.append(f"  {sym} = h_ε({shorthand(tr, sig)})")
        return "\n".join(lines)

    def _repr_html_(self) -> str:
        from ..render.display import canonical_html
        return canonical_html(self)

    def __repr__(self) -> str:
        nz = sum(1 for v in self.constants.values() if v == 0)
        return (f"<CanonicalFamily: {len(self.constants)} constants "
                f"({nz} vanish), {len(self.legend)} elementary expectations, "
                f"reduced={self.reduced}, symmetric={self.symmetric}>")


@dataclass
class RenormalizedEquation:
    spde: SPDE
    sig: Signature
    base: dict[int, dict[str, sympy.Expr]]
    unknowns: list[Unknown]
    per_component: dict[int, list[Counterterm]]      # component index -> its counterterms
    all_trees: list[DecoratedTree] = field(default_factory=list)  # full 𝓑_<0 (incl. F(τ*)=0 trees)

    @property
    def n_components(self) -> int:
        return self.sig.n_components

    @property
    def counterterms(self) -> list[Counterterm]:
        """Counterterms of component 0 (the scalar case)."""
        return self.per_component[0]

    def _display_subs(self, comp: int) -> dict[sympy.Symbol, sympy.Expr]:
        u = self.unknowns[comp]
        width = self.sig.width
        subs = {jet(comp, (0,) * width): u.field}
        for i in range(width):
            e_i = tuple(1 if j == i else 0 for j in range(width))
            subs[jet(comp, e_i)] = sympy.Derivative(u.field, u.coords[i])
        return subs

    def _all_display_subs(self) -> dict[sympy.Symbol, sympy.Expr]:
        subs: dict[sympy.Symbol, sympy.Expr] = {}
        for c in range(self.n_components):
            subs.update(self._display_subs(c))
        return subs

    def counterterm_rhs(self, comp: int = 0) -> sympy.Expr:
        """Σ_τ (k_τ/S(τ)) F_comp(τ*), with jets rendered as u and its derivatives."""
        subs = self._all_display_subs()
        out = sympy.Integer(0)
        for ct in self.per_component[comp]:
            out += ct.coefficient * ct.elem_diff.xreplace(subs)
        return sympy.expand(out)

    # --- pretty substitution: jet u^c_k → compact symbol u, u_x, u_xx, … ------ #

    def _pretty(self, expr: sympy.Expr) -> sympy.Expr:
        """Render jet variables as compact subscripted field symbols for output."""
        from ..render.tree import coord_names
        coords = coord_names(self.sig.width)
        subs = {}
        for s in expr.free_symbols:
            if is_jet(s):
                comp, k = jet_parts(s)
                name = self.unknowns[comp].name
                if not any(k):
                    subs[s] = sympy.Symbol(name)
                else:
                    sub = "".join(coords[i] for i, e in enumerate(k) for _ in range(e))
                    subs[s] = sympy.Symbol(f"{name}_{sub}")
        return expr.xreplace(subs)

    def original_rhs(self, comp: int = 0) -> sympy.Expr:
        """The un-renormalized RHS  f(u)ζ + g(u,∂u)  of component ``comp``, prettified."""
        expr = self.base[comp]["bullet"]
        for nz in self.spde.noises:
            expr = expr + self.base[comp][nz.name] * nz.symbol
        return self._pretty(expr)

    # --- the canonical (BPHZ) specialization, as computable objects ------------ #

    def canonical(self, reduced: bool = False, symmetric: bool = True) -> CanonicalFamily:
        """The canonical (BPHZ) values ``k_τ = h(S'₋τ)`` of the free constants, for a
        centered Gaussian noise (parity-reduced: odd-noise trees vanish).  ``reduced``
        folds in the exact identities of `renorm.scheme.reduction_subs`; ``symmetric``
        (default) additionally applies spatial reflection — set ``False`` for an
        anisotropic noise, where that identity is invalid.

        Raises ``NotImplementedError`` for a da Prato–Debussche-lifted equation: its Wick-power
        noises are correlated and non-Gaussian (`Signature.require_gaussian`)."""
        from ..structures import build_renormalization
        rs = build_renormalization(self.spde)
        cmap = {ct.tree: ct.constant
                for cts in self.per_component.values() for ct in cts}
        trees = [t for t in sorted(self.all_trees, key=lambda t: t.homogeneity(self.sig)._key())
                 if t in cmap]
        values = [rs.canonical_character(t) for t in trees]
        if reduced:
            from .scheme import reduction_subs
            subs = reduction_subs(rs._h, self.sig, symmetric)
            values = [sympy.expand(v.xreplace(subs)) for v in values]
        # keep only the h-symbols that survive, renumbered h0, h1, … contiguously
        used = sorted(set().union(*(v.free_symbols for v in values)) if values else set(),
                      key=lambda s: int(str(s)[1:]))
        relabel = {s: sympy.Symbol(f"h{i}") for i, s in enumerate(used)}
        values = [v.xreplace(relabel) for v in values]
        legend = dict(sorted(((relabel[sym], tr) for tr, sym in rs._h.items() if sym in relabel),
                             key=lambda x: int(str(x[0])[1:])))
        return CanonicalFamily(
            eq=self, structure=rs, reduced=reduced, symmetric=symmetric,
            constants=dict(zip(trees, values)),
            free_constants={t: cmap[t] for t in trees}, legend=legend)

    def bphz(self, symmetric: bool = True) -> CanonicalFamily:
        """Sugar: the fully reduced canonical family, ``canonical(reduced=True, …)``."""
        return self.canonical(reduced=True, symmetric=symmetric)

    # --- rendering (notes/output.md) ------------------------------------------ #

    def render(self, fmt: str = "text", canonical: bool = False, reduced: bool = False,
               symmetric: bool = True) -> str:
        from ..render.report import render as _render
        return _render(self, fmt, canonical, reduced, symmetric)

    def report(self, canonical: bool = False, reduced: bool = False, symmetric: bool = True) -> str:
        return self.render("text", canonical, reduced, symmetric)

    def latex_document(self, canonical: bool = False, reduced: bool = False,
                       symmetric: bool = True) -> str:
        from ..render.latex import latex_document
        return latex_document(self, canonical, reduced, symmetric)

    def to_json(self, canonical: bool = False, reduced: bool = False, symmetric: bool = True) -> str:
        return self.render("json", canonical, reduced, symmetric)

    def save(self, stem: str = "equation", outdir: str = "output",
             formats: tuple[str, ...] = ("text", "markdown", "json", "latex"), pdf: bool = True,
             canonical: bool = False, reduced: bool = False, symmetric: bool = True) -> list[Path]:
        """Write the exports to ``outdir/`` (default ``output/``); compile the PDF
        if ``pdflatex`` is on PATH.  Returns the paths written."""
        import shutil
        import subprocess
        from pathlib import Path

        ext = {"text": "txt", "markdown": "md", "json": "json", "latex": "tex", "html": "html"}
        d = Path(outdir)
        d.mkdir(parents=True, exist_ok=True)
        written = []
        for fmt in formats:
            p = d / f"{stem}.{ext[fmt]}"
            p.write_text(self.render(fmt, canonical, reduced, symmetric), encoding="utf-8")
            written.append(p)
        if pdf and "latex" in formats and shutil.which("pdflatex"):
            # twice: longtable writes its column widths to .aux on the first pass and
            # only lays them out correctly on the second (else wide tree cells overflow).
            for _ in range(2):
                subprocess.run(
                    ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"{stem}.tex"],
                    cwd=d, capture_output=True)
            for ext in ("aux", "log", "out", "toc", "lof", "lot", "fls",     # tidy up:
                        "fdb_latexmk", "synctex.gz", "nav", "snm", "vrb",     # drop every
                        "bbl", "blg", "idx", "ilg", "ind", "spl"):            # LaTeX aux file
                (d / f"{stem}.{ext}").unlink(missing_ok=True)
            if (d / f"{stem}.pdf").exists():
                written.append(d / f"{stem}.pdf")
        return written

    def summary(self) -> str:
        """Compact per-component listing (kept for error messages / quick demos)."""
        from ..render.tree import shorthand
        from .scheme import divergence_str
        lines = []
        for comp in range(self.n_components):
            cts = self.per_component[comp]
            lines.append(f"[{self.unknowns[comp].name}]  {len(cts)} counterterm(s)")
            for ct in sorted(cts, key=lambda c: c.homogeneity._key()):
                term = sympy.expand(ct.coefficient * self._pretty(ct.elem_diff))
                lines.append(f"   |τ|={str(ct.homogeneity):>8}  S={ct.symmetry_factor}  "
                             f"~{divergence_str(ct.homogeneity):<9}  "
                             f"+ {sympy.sstr(term):<12}  [τ = {shorthand(ct.tree, self.sig)}]")
        return "\n".join(lines)

    # --- self-display: summary at the prompt, markdown report in Jupyter ------ #

    def __str__(self) -> str:
        return self.summary()

    def __repr__(self) -> str:
        total = sum(len(v) for v in self.per_component.values())
        return (f"<RenormalizedEquation: {self.n_components} component(s), "
                f"{total} counterterm(s), {len(self.all_trees)} trees in 𝓑_<0>")

    def _repr_markdown_(self) -> str:
        return self.render("markdown")

    def _repr_html_(self) -> str:
        # what a notebook shows (HTML outranks Markdown): the family and every tree, drawn
        from ..render.display import summary_html
        return summary_html(self)

    def show(self, t: DecoratedTree) -> Picture:
        """``t`` drawn with this equation's signature — a lone tree's own display cannot
        know the other noise types or that it belongs to a system."""
        from ..render.display import tree_picture
        return tree_picture(t, self.sig)

    def explore(self, canonical: bool = False, reduced: bool = False, symmetric: bool = True,
                height: int = 640) -> Picture:
        """The interactive explorer (`render("html")`) embedded in a notebook cell."""
        from ..render.display import explore_html
        return explore_html(self, canonical, reduced, symmetric, height)
