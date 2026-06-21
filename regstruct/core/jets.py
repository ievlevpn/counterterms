"""Jet variables ``u^c_k`` (component c, spacetime multi-index k) as SymPy symbols.

``u^c_(0,...,0)`` is component c of the unknown; ``u^c_{e_i}`` is its ``∂_{x_i}``
derivative.  These are the variables in which the per-equation nonlinearities and
the elementary differentials ``F_a(τ*)`` are expressed (tourist_guide.tex 4524,
4915).  For a scalar equation there is a single component ``c = 0``.
"""
from __future__ import annotations

import sympy

from .homogeneity import MultiIndex


def jet(comp: int, k: MultiIndex) -> sympy.Symbol:
    return sympy.Symbol(f"u{comp}_" + "_".join(str(i) for i in k))


def is_jet(s) -> bool:
    return (isinstance(s, sympy.Symbol) and len(s.name) > 1
            and s.name[0] == "u" and s.name[1].isdigit())


def jet_parts(s: sympy.Symbol) -> tuple[int, MultiIndex]:
    parts = s.name[1:].split("_")
    return int(parts[0]), tuple(int(x) for x in parts[1:])
