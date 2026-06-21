"""The output object: the family of renormalized equations.

``RenormalizedEquation`` holds the base PDE plus one ``Counterterm`` per
negative-homogeneity tree.  The family is parametrised by the free constants
``k_τ`` (a chart of G⁻_ad): the τ-th counterterm contributes
``(k_τ / S(τ)) · F(τ*)`` to the right-hand side (tourist_guide.tex 4915).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import sympy

from ..core.jets import jet


@dataclass(frozen=True)
class Counterterm:
    tree: object
    homogeneity: object
    symmetry_factor: int
    elem_diff: object        # F(τ*) in jet variables
    constant: sympy.Symbol   # the free renormalization constant k_τ

    @property
    def coefficient(self):
        """The PDE coefficient k_τ / S(τ)."""
        return self.constant / self.symmetry_factor

    @property
    def term(self):
        """The counterterm contribution (k_τ / S(τ)) · F(τ*), in jet variables."""
        return self.coefficient * self.elem_diff


@dataclass
class RenormalizedEquation:
    spde: object
    sig: object
    base: dict
    width: int
    counterterms: list

    def _display_subs(self):
        u = self.spde.unknown
        subs = {jet((0,) * self.width): u.field}
        for i in range(self.width):
            e_i = tuple(1 if j == i else 0 for j in range(self.width))
            subs[jet(e_i)] = sympy.Derivative(u.field, u.coords[i])
        return subs

    def counterterm_rhs(self):
        """Σ_τ (k_τ/S(τ)) F(τ*), with jets rendered as u and its derivatives."""
        subs = self._display_subs()
        out = sympy.Integer(0)
        for ct in self.counterterms:
            out += ct.coefficient * ct.elem_diff.xreplace(subs)
        return sympy.expand(out)

    def summary(self) -> str:
        lines = [f"Renormalized equation: {len(self.counterterms)} counterterm(s)"]
        subs = self._display_subs()
        for ct in sorted(self.counterterms, key=lambda c: (c.homogeneity._key())):
            term = (ct.coefficient * ct.elem_diff.xreplace(subs))
            lines.append(f"  |τ|={str(ct.homogeneity):>8}  S={ct.symmetry_factor}  "
                         f"+ {sympy.simplify(term)}")
        return "\n".join(lines)
