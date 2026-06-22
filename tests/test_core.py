"""Convention #1 — the homogeneity vocabulary on the `Signature`:

* ``|Ξ| = β₀`` *directly* (the noise's Hölder regularity, **not** ``α−2``),
* ``|I^{(c)}_p τ| = |τ| + (m_c − |p|_𝔰)`` — the Schauder gain is operator order minus
  the scaled derivative,
* ``width = 1 + dim``;

plus the jet symbols (``u^c_k``) round-tripping through `jet`/`jet_parts`/`is_jet`.
Getting ``|Ξ|`` wrong by the off-by-2 ``α−2`` convention is silent and fatal, so it
is pinned here directly on the built `Signature`.
"""
from fractions import Fraction

from sympy import Rational, Symbol

from counterterms import jet, kappa
from counterterms.core.homogeneity import Homogeneity
from counterterms.core.jets import is_jet, jet_parts
from counterterms.equation.dsl import build_context

from tests.conftest import gkpz, gpam2, kpz


def test_noise_homogeneity_is_beta0_directly():
    sig, _b, _u = build_context(gkpz(reg=Rational(-1) - kappa))
    assert sig.node_homogeneity("xi") == Homogeneity(-1, -1)      # β₀, not β₀−2
    assert sig.node_homogeneity("bullet") == Homogeneity(0)       # ● (kernel) node is 0
    sig2, _b2, _u2 = build_context(kpz())
    assert sig2.node_homogeneity("xi") == Homogeneity(Rational(-3, 2), -1)


def test_edge_gain_is_order_minus_scaled_derivative():
    # 2nd-order parabolic, d=1: 𝔰=(2,1), m=2.
    sig, _b, _u = build_context(gkpz())
    assert sig.edge_gain(0, (0, 0)) == Homogeneity(2)     # plain kernel: +2
    assert sig.edge_gain(0, (0, 1)) == Homogeneity(1)     # ∂ₓ kernel: 2−1
    assert sig.edge_gain(0, (0, 2)) == Homogeneity(0)     # ∂ₓ² kernel: 2−2
    assert sig.edge_gain(0, (1, 0)) == Homogeneity(0)     # ∂_t kernel: 2−2 (t weighs 2)


def test_width_is_one_plus_dim():
    assert build_context(gkpz())[0].width == 2            # d=1 ⇒ (t,x)
    assert build_context(gpam2())[0].width == 3           # d=2 ⇒ (t,x,y)


def test_jet_round_trip():
    s = jet(1, (0, 2))
    assert is_jet(s)
    assert jet_parts(s) == (1, (0, 2))
    # component 0, undifferentiated and ∂ₓ are distinct symbols
    assert jet(0, (0, 0)) != jet(0, (0, 1))


def test_is_jet_rejects_ordinary_symbols():
    assert not is_jet(Symbol("x"))
    assert not is_jet(Symbol("u"))            # bare name, no component digit
    assert not is_jet(Fraction(1))
