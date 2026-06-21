"""regstruct — symbolic renormalized equations for singular SPDEs.

Phase 1: SPDE → family of renormalized equations with free constants.
See notes/architecture.md and notes/initial_plan.md.
"""
from .equation.dsl import SPDE, Noise, Parabolic, Unknown, kappa
from .core.jets import jet

__all__ = ["SPDE", "Noise", "Parabolic", "Unknown", "kappa", "jet", "renormalize"]


def renormalize(spde):
    from .api import renormalize as _r
    return _r(spde)
