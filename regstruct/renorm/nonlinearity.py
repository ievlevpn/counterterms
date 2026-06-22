"""The elementary-differential (Υ) map: ``τ ↦ F_a(τ*)`` — the counterterm engine.

Base cases (user data, per output equation ``a``): ``F_a(∘_j*) = f_{a,j}(u)``,
``F_a(●*) = g_a(u, ∂u)``, ``F_a(red*) = 0``.  Recursion (tourist_guide.tex 4524 /
4915), with each child edge carrying the component ``c`` of the equation it plants:

    F_a(τ*) = ( Πᵢ F_{cᵢ}(τ_i*) ) · ( D^n  Πᵢ ∂_{(cᵢ, p_i)} ) F_a(b*)

``∂_{(c,p)} = ∂/∂u^c_p`` (component ``c``'s jet), and the total derivative
``D_l = Σ_{c,k} u^c_{k+e_l} ∂_{u^c_k}`` runs over all components.  The ``∂_{(cᵢ,p_i)}``
are applied *before* ``D^n`` (they do not commute) and hit **all** slots of ``g``.
The child equation index ``cᵢ`` comes from the edge — this is how systems couple.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import sympy

from ..core.jets import is_jet, jet, jet_parts

if TYPE_CHECKING:
    from ..core.homogeneity import MultiIndex
    from ..core.signature import Signature
    from ..trees.tree import DecoratedTree


def _D(expr: sympy.Expr, ell: int, width: int) -> sympy.Expr:
    """Total derivative ``D_ℓ = Σ_{c,k} u^c_{k+e_ℓ} ∂_{u^c_k}`` (over all components)."""
    e_ell = tuple(1 if j == ell else 0 for j in range(width))
    res = sympy.Integer(0)
    for s in list(expr.free_symbols):
        if is_jet(s):
            c, k = jet_parts(s)
            shifted = tuple(a + b for a, b in zip(k, e_ell))
            res += jet(c, shifted) * sympy.diff(expr, s)
    return res


def _Dn(expr: sympy.Expr, n: MultiIndex, width: int) -> sympy.Expr:
    for ell in range(width):
        for _ in range(n[ell]):
            expr = _D(expr, ell, width)
    return expr


def elem_diff(t: DecoratedTree, comp: int, base: dict, sig: Signature) -> sympy.Expr:
    """``F_comp(t*)`` — the elementary differential of ``t`` for output equation ``comp``."""
    expr = base[comp][t.node_type]
    for (c, p, _sub) in t.children:                  # Πᵢ ∂_{(cᵢ, p_i)}
        expr = sympy.diff(expr, jet(c, p))
    expr = _Dn(expr, t.node_dec, sig.width)          # D^n
    for (c, _p, sub) in t.children:                  # × Πᵢ F_{cᵢ}(τ_i*)
        expr = expr * elem_diff(sub, c, base, sig)
    return sympy.expand(expr)
