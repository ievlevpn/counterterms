"""Subcriticality of the derived rule (equation/rule.py) — the principled replacement
for the old hardcoded ``β₀∈(−2,0)`` bound.  Subcritical ⟺ ``β₀ > −(operator order)``,
so higher-order operators legitimately relax the threshold (the hardcoded ``−2`` was a
latent over-restriction).
"""
import warnings

import pytest
from sympy import Derivative, Function, Rational

from counterterms import Noise, Parabolic, SPDE, Unknown, kappa
from counterterms.equation.dsl import build_context
from counterterms.equation.generate import generate_counterterms

f, g = Function("f"), Function("g")


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


def test_d2_total_gradient_bound_in_dim2():
    """Assumption D2 bounds a node by TOTAL degree 2 in ∂u (tex 5337-5340: at most two
    gradient edges I_{e_i}, I_{e_j}).  In d≥2 a direction-mixing nonlinearity g(u)(∂₁u+∂₂u)²
    has per-direction caps that sum to 4; only the total `grad_budget` enforces the real
    bound.  Without it, generation over-produced trees with 3–4 gradient edges (all Υ-zero,
    so the renormalised equation was unaffected, but the raw basis was wrong)."""
    u = Unknown("u", 2)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    rhs = g(u.field) * (Derivative(u.field, u.x[0]) + Derivative(u.field, u.x[1])) ** 2
    sig, _b, _u = build_context(
        SPDE(noises=[xi], operator=Parabolic(dim=2), unknown=u, rhs=rhs))
    assert sig.grad_budget["bullet"] == 2

    def max_grad_edges(t):
        here = sum(1 for (_c, p, _s) in t.children if any(p))
        return max([here, *(max_grad_edges(s) for (_c, _p, s) in t.children)])

    trees = generate_counterterms(sig)
    assert all(max_grad_edges(t) <= 2 for t in trees)   # no node exceeds 2 gradient edges
    assert len(trees) == 8                              # was 11 (3 spurious) before the fix


def test_nonsingular_noise_rejected():
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=Rational(1))       # β₀ ≥ 0
    with pytest.raises(ValueError, match="not singular"):
        build_context(SPDE(noises=[xi], operator=Parabolic(dim=1), unknown=u,
                           rhs=f(u.field) * xi.symbol))
