"""Jet variables ``u_k`` (k a spacetime multi-index) as SymPy symbols.

``u_(0,...,0)`` is the unknown itself; ``u_{e_i}`` is its ``∂_{x_i}`` derivative.
These are the variables in which the base nonlinearities and the elementary
differentials ``F(τ*)`` are expressed (tourist_guide.tex 4524, 4915).
"""
from __future__ import annotations

import sympy

from .homogeneity import MultiIndex


def jet(k: MultiIndex) -> sympy.Symbol:
    return sympy.Symbol("u_" + "_".join(str(i) for i in k))


def is_jet(s) -> bool:
    return isinstance(s, sympy.Symbol) and s.name.startswith("u_")


def jet_index(s: sympy.Symbol) -> MultiIndex:
    return tuple(int(x) for x in s.name[2:].split("_"))
