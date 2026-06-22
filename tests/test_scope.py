"""Scope enforcement — every rejection the DSL must raise, and the one warning.

CLAUDE.md fixes the MVP scope (scalar→systems, β₀∈(−2,0), affine-in-noise, g ≤
quadratic in ∂u, |p|_𝔰≤1, 2nd-order parabolic L).  Each out-of-scope input must
fail with a clear `ValueError` (not a wrong answer); each path is pinned here.
"""
import pytest
from sympy import Derivative, Function, Rational, sin

from counterterms import FractionalHeat, Noise, Parabolic, SPDE, Unknown, kappa

f = Function("f")
g = Function("g")


def _u_xi(reg=Rational(-1) - kappa):
    return Unknown("u", 1), Noise("xi", regularity=reg)


def test_regularity_must_be_affine_in_kappa():
    with pytest.raises(ValueError, match="affine in κ"):
        Noise("xi", regularity=kappa ** 2)


def test_noise_must_be_affine():
    u, xi = _u_xi()
    with pytest.raises(ValueError, match="affine in noise"):
        SPDE(noises=[xi], operator=Parabolic(dim=1), unknown=u,
             rhs=f(u.field) * xi.symbol ** 2).renormalize()


def test_noise_must_appear_only_affinely():
    # a non-polynomial noise dependence (here sin ξ) leaves noise in the g-part.
    u, xi = _u_xi()
    with pytest.raises(ValueError, match="noise-free part still has"):
        SPDE(noises=[xi], operator=Parabolic(dim=1), unknown=u,
             rhs=sin(xi.symbol)).renormalize()


def test_runaway_tree_set_fails_fast():
    # fractional order + quadratic-gradient nonlinearity makes the negative-tree set
    # intractably large; the generator must raise (not hang) via the pool-size backstop.
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    rhs = f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2
    with pytest.raises(RuntimeError, match="intractably large"):
        SPDE(noises=[xi], operator=FractionalHeat(dim=1, sigma=Rational(3, 4)),
             unknown=u, rhs=rhs).renormalize()


def test_beta0_at_most_minus_two_rejected():
    # β₀≤−2 (Φ⁴₃, sine-Gordon, …) is supercritical — needs a da Prato–Debussche pre-pass.
    u = Unknown("u", 3)
    xi = Noise("xi", regularity=Rational(-5, 2) - kappa)
    with pytest.raises(ValueError, match="da Prato"):
        SPDE(noises=[xi], operator=Parabolic(dim=3), unknown=u,
             rhs=-u.field ** 3 + xi.symbol).renormalize()


def test_noise_coefficient_must_depend_on_u_only():
    # a derivative-of-u noise coefficient (∂ₓu · ξ) is out of scope.
    u, xi = _u_xi()
    with pytest.raises(ValueError, match="must depend on"):
        SPDE(noises=[xi], operator=Parabolic(dim=1), unknown=u,
             rhs=Derivative(u.field, u.x[0]) * xi.symbol).renormalize()


def test_singular_derivative_factor_rejected():
    # a |p|_𝔰>1 factor in g (here ∂ₓ²u) is too singular for the MVP.
    u, xi = _u_xi()
    with pytest.raises(ValueError, match=r"\|p\|_𝔰>1"):
        SPDE(noises=[xi], operator=Parabolic(dim=1), unknown=u,
             rhs=xi.symbol + Derivative(u.field, u.x[0], 2)).renormalize()


def test_cubic_in_gradient_rejected():
    # g more than quadratic in ∂u violates Assumption D2.
    u, xi = _u_xi()
    with pytest.raises(ValueError, match="quadratic"):
        SPDE(noises=[xi], operator=Parabolic(dim=1), unknown=u,
             rhs=xi.symbol + Derivative(u.field, u.x[0]) ** 3).renormalize()


def test_nonparabolic_order_warns():
    # order≠2 computes homogeneities but warns: the RS theory is unverified there.
    with pytest.warns(UserWarning, match="2nd-order"):
        Parabolic(dim=1, order=4)
