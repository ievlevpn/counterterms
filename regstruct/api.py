"""Top-level pipeline: SPDE → family of renormalized equations (free constants)."""
from __future__ import annotations

from typing import TYPE_CHECKING

import sympy

from .equation.dsl import build_context
from .equation.generate import generate_counterterms
from .renorm.equation import Counterterm, RenormalizedEquation
from .renorm.nonlinearity import elem_diff

if TYPE_CHECKING:
    from .equation.dsl import SPDE


def renormalize(spde: SPDE) -> RenormalizedEquation:
    sig, base, unknowns = build_context(spde)
    trees = generate_counterterms(sig)

    per_component = {a: [] for a in range(sig.n_components)}
    idx = 0
    for t in trees:
        # one tree may contribute to several components; the constant k_τ is shared.
        contribs = {}
        for a in range(sig.n_components):
            ed = elem_diff(t, a, base, sig)
            if ed != 0:                       # Assumption-D2 safety: drop F_a(τ*) = 0
                contribs[a] = ed
        if not contribs:
            continue
        k = sympy.Symbol(f"k_{idx}")
        idx += 1
        S = t.symmetry_factor()
        hom = t.homogeneity(sig)
        for a, ed in contribs.items():
            per_component[a].append(Counterterm(
                tree=t, homogeneity=hom, symmetry_factor=S, elem_diff=ed, constant=k))

    return RenormalizedEquation(spde=spde, sig=sig, base=base, unknowns=unknowns,
                                per_component=per_component, all_trees=trees)
