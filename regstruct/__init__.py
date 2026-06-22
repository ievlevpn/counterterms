"""regstruct — symbolic renormalized equations for singular SPDEs.

Phase 1: SPDE → family of renormalized equations with free constants.
See notes/architecture.md and notes/initial_plan.md.
"""
from .equation.dsl import SPDE, Noise, Parabolic, Unknown, kappa
from .core.jets import jet

__all__ = ["SPDE", "Noise", "Parabolic", "Unknown", "kappa", "jet", "renormalize",
           "build_renormalization", "build_regularity_structure",
           "build_renormalization_group", "daprato_lift"]


def renormalize(spde):
    from .api import renormalize as _r
    return _r(spde)


def build_renormalization(spde):
    """The negative (renormalisation) structure: δ, δ⁻, S'₋, and the BHZ character."""
    from .structures import build_renormalization as _b
    return _b(spde)


def build_regularity_structure(spde, gamma=1):
    """The γ-truncated regularity structure (T, T⁺): graded model basis + Δ, Δ⁺."""
    from .structures import build_regularity_structure as _b
    return _b(spde, gamma)


def build_renormalization_group(spde):
    """The renormalization group G⁻: characters on U⁻ with the convolution group law."""
    from .structures import build_renormalization_group as _b
    return _b(spde)


def daprato_lift(spde):
    """da Prato–Debussche lift: turn a supercritical additive-noise polynomial SPDE
    into the subcritical equation for the remainder v = u − X (then `.renormalize()`)."""
    from .equation.daprato import daprato_lift as _l
    return _l(spde)
