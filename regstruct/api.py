"""Top-level pipeline: SPDE → family of renormalized equations (free constants)."""
from __future__ import annotations

import sympy

from .equation.dsl import build_context
from .equation.generate import generate_counterterms
from .renorm.equation import Counterterm, RenormalizedEquation
from .renorm.nonlinearity import elem_diff


def renormalize(spde) -> RenormalizedEquation:
    sig, base, width = build_context(spde)
    trees = generate_counterterms(sig)

    counterterms = []
    for t in trees:
        ed = elem_diff(t, base, width)
        if ed == 0:                      # Assumption-D2 safety: drop F(τ*) = 0 trees
            continue
        k = sympy.Symbol(f"k_{len(counterterms)}")
        counterterms.append(Counterterm(
            tree=t,
            homogeneity=t.homogeneity(sig),
            symmetry_factor=t.symmetry_factor(),
            elem_diff=ed,
            constant=k,
        ))
    return RenormalizedEquation(spde=spde, sig=sig, base=base, width=width,
                                counterterms=counterterms)
