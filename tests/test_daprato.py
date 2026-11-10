"""da Prato–Debussche pre-pass (`equation/daprato.py`): lift a supercritical
additive-noise polynomial SPDE into the subcritical equation for the remainder
`v = u − X`, which the engine then renormalises. This unlocks the Φ⁴ family
(β₀ ≤ −order), previously rejected outright.
"""
from fractions import Fraction

import pytest
from sympy import Function, Rational, sin

from counterterms import Noise, Parabolic, SPDE, Unknown, daprato_lift, kappa


def _phi4(d, beta0_std):
    u = Unknown("u", d)
    xi = Noise("xi", regularity=Rational(beta0_std) - kappa)
    return SPDE(noises=[xi], operator=Parabolic(dim=d), unknown=u, rhs=xi.symbol - u.field ** 3)


def test_phi43_is_rejected_directly_but_lifts_and_renormalises():
    phi43 = _phi4(3, Rational(-5, 2))                 # β₀ = −5/2 − κ, supercritical
    with pytest.raises(ValueError, match="subcritical"):
        phi43.renormalize()
    res = daprato_lift(phi43).renormalize()
    assert res.counterterms                            # nonempty renormalised family


def test_lift_wick_power_regularities():
    # :X^k: has regularity k·(β₀+order); Φ⁴₃: α_X = −5/2+2 = −1/2.
    lifted = daprato_lift(_phi4(3, Rational(-5, 2)))
    regs = {n.name: (n.std, n.kap) for n in lifted.noises}
    assert regs == {
        "X1": (Fraction(-1, 2), Fraction(-1)),         # X      = K∗ξ
        "X2": (Fraction(-1), Fraction(-2)),            # :X²:
        "X3": (Fraction(-3, 2), Fraction(-3)),         # :X³:
    }


def test_lift_is_the_taylor_expansion_of_the_nonlinearity():
    # −(v+X)³ = −v³ − 3v²X − 3v:X²: − :X³:
    lifted = daprato_lift(_phi4(3, Rational(-5, 2)))
    (v, _op, rhs) = lifted.equations[0]
    X1, X2, X3 = (n.symbol for n in lifted.noises)
    assert rhs.coeff(X1, 1) == -3 * v.field ** 2       # P'(v)
    assert rhs.coeff(X2, 1) == -3 * v.field            # P''(v)/2
    assert rhs.coeff(X3, 1) == -1                       # P'''(v)/6
    deterministic = rhs.subs({X1: 0, X2: 0, X3: 0})
    assert deterministic == -v.field ** 3              # P(v)


def test_phi42_lifts_to_a_subcritical_equation():
    # Φ⁴₂ (β₀=−2−κ): α_X = −κ, so the Wick noises sit just below 0 and v is subcritical.
    res = daprato_lift(_phi4(2, Rational(-2))).renormalize()
    assert res.counterterms
    assert all(c.homogeneity.is_negative() for c in res.counterterms)


def test_lift_rejects_multiplicative_noise():
    u = Unknown("u", 3)
    xi = Noise("xi", regularity=Rational(-5, 2) - kappa)
    f = Function("f")
    with pytest.raises(ValueError, match="additive"):
        daprato_lift(SPDE(noises=[xi], operator=Parabolic(dim=3), unknown=u,
                          rhs=f(u.field) * xi.symbol - u.field ** 3))


def test_lift_rejects_non_polynomial_nonlinearity():
    u = Unknown("u", 3)
    xi = Noise("xi", regularity=Rational(-5, 2) - kappa)
    with pytest.raises(ValueError, match="polynomial"):
        daprato_lift(SPDE(noises=[xi], operator=Parabolic(dim=3), unknown=u,
                          rhs=xi.symbol + sin(u.field)))


def test_lifted_noises_are_tagged_and_gaussian_only_paths_refuse():
    # The Wick powers :X^k: are correlated (one underlying ξ) and non-Gaussian for k≥2, so the
    # independent-centered-Gaussian parity rule must not be applied to them silently.
    lifted = daprato_lift(_phi4(3, Rational(-5, 2)))
    assert [n.chaos for n in lifted.noises] == [1, 2, 3]
    res = lifted.renormalize()
    assert res.sig.non_gaussian == {"X2", "X3"}
    for call in (res.canonical, res.bphz, lambda: res.report(canonical=True)):
        with pytest.raises(NotImplementedError, match="Wick powers"):
            call()
    # the law-agnostic pieces still work: the free family and the raw h∘S'₋ character
    from counterterms import build_renormalization
    rs = build_renormalization(lifted)
    assert len(res.counterterms) == 9 and "canonical=True" in res.report()
    assert all(rs.bhz_character(t).free_symbols for t in rs.divergent)
    with pytest.raises(NotImplementedError, match="Wick powers"):
        rs.canonical_character(rs.divergent[0])


def test_canonical_view_names_the_noises():
    # CanonicalFamily.__str__ labels trees with the family's signature (Ξ_xi vs Ξ_eta), not a
    # per-tree guess that prints every single-noise tree as a bare Ξ.
    from sympy import Function
    u = Unknown("u", 1)
    xi, eta = Noise("xi", regularity=-1 - kappa), Noise("eta", regularity=-1 - kappa)
    f, g = Function("f"), Function("g")
    eq = SPDE(noises=[xi, eta], operator=Parabolic(dim=1), unknown=u,
              rhs=f(u.field) * xi.symbol + g(u.field) * eta.symbol).renormalize()
    txt = str(eq.canonical())
    assert "Ξ_xi" in txt and "Ξ_eta" in txt and "τ = Ξ]" not in txt and "τ = Ξ·" not in txt
