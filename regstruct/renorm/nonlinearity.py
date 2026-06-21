"""The elementary-differential (Υ) map: ``τ ↦ F(τ*)`` — the counterterm engine.

Base cases (user data): ``F(∘_j*) = f_j(u₀)``, ``F(●*) = g(u₀, ∂u)``,
``F(red*) = 0``.  Recursion (tourist_guide.tex 4524 / 4915):

    F(τ*) = ( Πᵢ F(τ_i*) ) · ( D^n  Πᵢ ∂_{p_i} ) F(b*)

with ``∂_p = ∂/∂u_p`` and the total derivative ``D_i = Σ_k u_{k+e_i} ∂_k``.
The ``∂_{p_i}`` are applied *before* ``D^n`` (they do not commute), and the
``∂_{p_i}`` hit **all** slots of ``g`` (function and derivative arguments).
SymPy does the differentiation; this is the only place it enters the math.
"""
from __future__ import annotations

import sympy

from ..core.jets import is_jet, jet, jet_index


def _D(expr, i: int, width: int):
    """Total derivative ``D_i = Σ_k u_{k+e_i} ∂_{u_k}``."""
    e_i = tuple(1 if j == i else 0 for j in range(width))
    res = sympy.Integer(0)
    for s in list(expr.free_symbols):
        if is_jet(s):
            k = jet_index(s)
            shifted = tuple(a + b for a, b in zip(k, e_i))
            res += jet(shifted) * sympy.diff(expr, s)
    return res


def _Dn(expr, n, width: int):
    for i in range(width):
        for _ in range(n[i]):
            expr = _D(expr, i, width)
    return expr


def elem_diff(t, base_F, width: int):
    expr = base_F[t.node_type]
    for (_op, p, _sub) in t.children:        # Πᵢ ∂_{p_i}
        expr = sympy.diff(expr, jet(p))
    expr = _Dn(expr, t.node_dec, width)      # D^n
    for (_op, _p, sub) in t.children:        # × Πᵢ F(τ_i*)
        expr = expr * elem_diff(sub, base_F, width)
    return sympy.expand(expr)
