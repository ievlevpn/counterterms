"""The output object: the family of renormalized equations.

For a system, each component ``a`` gets its own renormalized equation; the trees
(and their constants ``k_τ``) are **shared** across components — the same `k_τ`
weights `F_a(τ*)` in every equation where it is nonzero (tourist_guide.tex 4915).
A scalar problem is the one-component case.
"""
from __future__ import annotations

from dataclasses import dataclass

import sympy

from ..core.jets import jet


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

    def summary(self) -> str:
        subs = self._all_display_subs()
        lines = []
        for comp in range(self.n_components):
            cts = self.per_component[comp]
            name = self.unknowns[comp].name
            lines.append(f"[{name}]  {len(cts)} counterterm(s)")
            for ct in sorted(cts, key=lambda c: c.homogeneity._key()):
                term = ct.coefficient * ct.elem_diff.xreplace(subs)
                lines.append(f"   |τ|={str(ct.homogeneity):>8}  S={ct.symmetry_factor}  "
                             f"+ {sympy.simplify(term)}")
        return "\n".join(lines)
