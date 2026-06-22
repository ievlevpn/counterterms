"""Subcriticality of the derived rule (equation/rule.py) — the principled replacement
for the old hardcoded ``β₀∈(−2,0)`` bound.  Subcritical ⟺ ``β₀ > −(operator order)``,
so higher-order operators legitimately relax the threshold (the hardcoded ``−2`` was a
latent over-restriction).
"""
import warnings

import pytest
from sympy import Function, Rational

from regstruct import Noise, Parabolic, SPDE, Unknown, kappa
from regstruct.equation.dsl import build_context

f = Function("f")


def _build(order, beta0_std):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")           # silence the order≠2 advisory
        u = Unknown("u", 1)
        xi = Noise("xi", regularity=Rational(beta0_std) - kappa)
        return build_context(SPDE(noises=[xi], operator=Parabolic(dim=1, order=order),
                                  unknown=u, rhs=f(u.field) * xi.symbol))


def test_second_order_subcritical_threshold():
    _build(2, -1)                                  # β₀=−1 > −2: subcritical, fine
    with pytest.raises(ValueError, match="subcritical"):
        _build(2, -2)                              # β₀=−2: supercritical


def test_higher_operator_order_relaxes_the_threshold():
    # order 4 ⇒ subcritical down to β₀ > −4 (the hardcoded −2 wrongly rejected this).
    _build(4, -3)                                  # now accepted
    with pytest.raises(ValueError, match="subcritical"):
        _build(4, -5)                              # β₀=−5 ≤ −4: supercritical


def test_phi43_reported_as_supercritical():
    # β₀=−5/2, 2nd order (Φ⁴₃ flavour): supercritical, points at da Prato–Debussche.
    u = Unknown("u", 3)
    xi = Noise("xi", regularity=Rational(-5, 2) - kappa)
    with pytest.raises(ValueError, match="da Prato"):
        build_context(SPDE(noises=[xi], operator=Parabolic(dim=3), unknown=u,
                           rhs=-u.field ** 3 + xi.symbol))


def test_nonsingular_noise_rejected():
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=Rational(1))       # β₀ ≥ 0
    with pytest.raises(ValueError, match="not singular"):
        build_context(SPDE(noises=[xi], operator=Parabolic(dim=1), unknown=u,
                           rhs=f(u.field) * xi.symbol))
