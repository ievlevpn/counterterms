"""The output object: the family of renormalized equations.

For a system, each component ``a`` gets its own renormalized equation; the trees
(and their constants ``k_τ``) are **shared** across components — the same `k_τ`
weights `F_a(τ*)` in every equation where it is nonzero (tourist_guide.tex 4915).
A scalar problem is the one-component case.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import sympy

from ..core.jets import is_jet, jet, jet_parts


@dataclass(frozen=True)
class Counterterm:
    tree: object
    homogeneity: object
    symmetry_factor: int
    elem_diff: object        # F_a(τ*) in jet variables
    constant: sympy.Symbol   # the free renormalization constant k_τ (shared across components)

    @property
    def coefficient(self):
        return self.constant / self.symmetry_factor      # k_τ / S(τ)

    @property
    def term(self):
        return self.coefficient * self.elem_diff


@dataclass
class RenormalizedEquation:
    spde: object
    sig: object
    base: dict
    unknowns: list
    per_component: dict          # component index -> list[Counterterm]
    all_trees: list = field(default_factory=list)   # full 𝓑_<0 (incl. F(τ*)=0 trees)

    @property
    def n_components(self) -> int:
        return self.sig.n_components

    @property
    def counterterms(self):
        """Counterterms of component 0 (the scalar case)."""
        return self.per_component[0]

    def _display_subs(self, comp: int):
        u = self.unknowns[comp]
        width = self.sig.width
        subs = {jet(comp, (0,) * width): u.field}
        for i in range(width):
            e_i = tuple(1 if j == i else 0 for j in range(width))
            subs[jet(comp, e_i)] = sympy.Derivative(u.field, u.coords[i])
        return subs

    def _all_display_subs(self):
        subs = {}
        for c in range(self.n_components):
            subs.update(self._display_subs(c))
        return subs

    def counterterm_rhs(self, comp: int = 0):
        """Σ_τ (k_τ/S(τ)) F_comp(τ*), with jets rendered as u and its derivatives."""
        subs = self._all_display_subs()
        out = sympy.Integer(0)
        for ct in self.per_component[comp]:
            out += ct.coefficient * ct.elem_diff.xreplace(subs)
        return sympy.expand(out)

    # --- pretty substitution: jet u^c_k → compact symbol u, u_x, u_xx, … ------ #

    def _pretty(self, expr):
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

    def original_rhs(self, comp: int = 0):
        """The un-renormalized RHS  f(u)ζ + g(u,∂u)  of component ``comp``, prettified."""
        expr = self.base[comp]["bullet"]
        for nz in self.spde.noises:
            expr = expr + self.base[comp][nz.name] * nz.symbol
        return self._pretty(expr)

    # --- rendering (notes/output.md) ------------------------------------------ #

    def render(self, fmt: str = "text", sections="all") -> str:
        from ..render.report import render as _render
        return _render(self, fmt, sections)

    def report(self) -> str:
        return self.render("text")

    def latex_document(self) -> str:
        from ..render.latex import latex_document
        return latex_document(self)

    def to_json(self) -> str:
        return self.render("json")

    def save(self, stem: str = "equation", outdir: str = "output",
             formats=("text", "markdown", "json", "latex"), pdf: bool = True) -> list:
        """Write the exports to ``outdir/`` (default ``output/``); compile the PDF
        if ``pdflatex`` is on PATH.  Returns the paths written."""
        import shutil
        import subprocess
        from pathlib import Path

        ext = {"text": "txt", "markdown": "md", "json": "json", "latex": "tex"}
        d = Path(outdir)
        d.mkdir(parents=True, exist_ok=True)
        written = []
        for fmt in formats:
            p = d / f"{stem}.{ext[fmt]}"
            p.write_text(self.render(fmt), encoding="utf-8")
            written.append(p)
        if pdf and "latex" in formats and shutil.which("pdflatex"):
            subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", f"{stem}.tex"],
                cwd=d, capture_output=True)
            for junk in (f"{stem}.aux", f"{stem}.log"):    # keep output/ tidy
                (d / junk).unlink(missing_ok=True)
            if (d / f"{stem}.pdf").exists():
                written.append(d / f"{stem}.pdf")
        return written

    def summary(self) -> str:
        """Compact per-component listing (kept for error messages / quick demos)."""
        lines = []
        for comp in range(self.n_components):
            cts = self.per_component[comp]
            lines.append(f"[{self.unknowns[comp].name}]  {len(cts)} counterterm(s)")
            for ct in sorted(cts, key=lambda c: c.homogeneity._key()):
                term = sympy.expand(ct.coefficient * self._pretty(ct.elem_diff))
                lines.append(f"   |τ|={str(ct.homogeneity):>8}  S={ct.symmetry_factor}  "
                             f"+ {sympy.sstr(term)}")
        return "\n".join(lines)
