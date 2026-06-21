"""Breadth checks: KPZ homogeneity rows, gPAM, and scope rejections."""
import pytest
import sympy
from sympy import Derivative, Function, Rational

from regstruct import Noise, Parabolic, SPDE, Unknown, kappa
from regstruct.core.homogeneity import Homogeneity


def test_kpz_homogeneity_rows():
    # KPZ: (∂_t − Δ)u = (∂_x u)^2 + ξ, space-time WN d=1, β0 = −3/2 − κ.
    # The counterterm trees must populate exactly the negative-homogeneity rows
    # of the paper's KPZ table (tourist_guide.tex 6028-6063).
    u = Unknown("u", dim=1)
    xi = Noise("xi", regularity=Rational(-3, 2) - kappa)
    res = SPDE(operator=Parabolic(dim=1), unknown=u, noises=[xi],
               rhs=1 * xi.symbol + Derivative(u.field, u.x[0]) ** 2).renormalize()

    rows = {
        Homogeneity(Rational(-3, 2), -1),   # β0
        Homogeneity(-1, -2),                 # 2β0 + 2
        Homogeneity(Rational(-1, 2), -3),    # 3β0 + 4
        Homogeneity(Rational(-1, 2), -1),    # β0 + 1
        Homogeneity(0, -2),                  # 2β0 + 3
        Homogeneity(0, -4),                  # 4β0 + 6
    }
    assert {c.homogeneity for c in res.counterterms} == rows
    # the bare noise ∘ is the classical KPZ constant counterterm
    assert any(c.tree.nodes() == 1 and not c.tree.children and not any(c.tree.node_dec)
               for c in res.counterterms)


def test_gpam_d2():
    # gPAM in d=2: (∂_t − Δ)u = f(u) ξ, space WN, β0 = −1 − κ; g ≡ 0.
    f = Function("f")
    u = Unknown("u", dim=2)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    res = SPDE(operator=Parabolic(dim=2), unknown=u, noises=[xi],
               rhs=f(u.field) * xi.symbol).renormalize()
    assert len(res.counterterms) == 4
    assert {c.homogeneity for c in res.counterterms} == {
        Homogeneity(-1, -1), Homogeneity(0, -1), Homogeneity(0, -2)}


def test_phi43_rejected():
    # β0 = −5/2 ≤ −2 is out of the direct rule reader (needs da Prato–Debussche).
    u = Unknown("u", dim=3)
    xi = Noise("xi", regularity=Rational(-5, 2) - kappa)
    with pytest.raises(ValueError, match="da Prato"):
        SPDE(operator=Parabolic(dim=3), unknown=u, noises=[xi],
             rhs=-u.field ** 3 + xi.symbol).renormalize()


def test_nonaffine_noise_rejected():
    u = Unknown("u", dim=1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    f = Function("f")
    with pytest.raises(ValueError, match="affine in noise"):
        SPDE(operator=Parabolic(dim=1), unknown=u, noises=[xi],
             rhs=f(u.field) * xi.symbol ** 2).renormalize()


def test_cubic_in_gradient_rejected():
    # g more than quadratic in ∂u violates Assumption D2.
    u = Unknown("u", dim=1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    with pytest.raises(ValueError, match="quadratic"):
        SPDE(operator=Parabolic(dim=1), unknown=u, noises=[xi],
             rhs=xi.symbol + Derivative(u.field, u.x[0]) ** 3).renormalize()
