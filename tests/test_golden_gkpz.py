"""Golden test: the gKPZ example, tourist_guide.tex 5996-6012.

    (∂_t − Δ + 1) u = f(u) ζ + g(u) (∂_x u)²,   d = 1,   ζ ∈ C^{−1−κ}

must produce exactly five counterterms (Theorem ThmRenormPDEs):

    k(∘)            f(u)
    k(∘1)           f'(u) ∂_x u
    k(τ)/1 · 2      f(u) g(u) ∂_x u
    k(τ)            f(u) f'(u)
    k(τ)/2 · 2      f²(u) g(u)

We check the (homogeneity, S(τ), F(τ*)) triple of each, in jet variables.
The factor-2 in the third term (S=1) and the S=2 cancellation in the fifth are
the load-bearing combinatorial checks.
"""
from fractions import Fraction

import sympy
from sympy import Derivative, Function, Rational

from regstruct import Noise, Parabolic, SPDE, Unknown, jet, kappa
from regstruct.core.homogeneity import Homogeneity

U0 = jet((0, 0))
U1 = jet((0, 1))
f = Function("f")
g = Function("g")


def _build():
    u = Unknown("u", dim=1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    x = u.x[0]
    rhs = f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, x) ** 2
    spde = SPDE(operator=Parabolic(dim=1, mass=1), unknown=u, noises=[xi], rhs=rhs)
    return spde.renormalize()


def test_gkpz_counterterms():
    res = _build()

    expected = [
        (Homogeneity(-1, -1), 1, f(U0)),                          # ∘
        (Homogeneity(0, -1), 1, U1 * Derivative(f(U0), U0)),      # ∘1
        (Homogeneity(0, -1), 1, 2 * f(U0) * g(U0) * U1),          # ●—I_{(0,1)}—∘
        (Homogeneity(0, -2), 1, f(U0) * Derivative(f(U0), U0)),   # ∘—I—∘
        (Homogeneity(0, -2), 2, 2 * f(U0) ** 2 * g(U0)),          # ●—two I_{(0,1)}—∘∘
    ]

    assert len(res.counterterms) == len(expected), (
        f"expected {len(expected)} counterterms, got {len(res.counterterms)}:\n{res.summary()}"
    )

    produced = list(res.counterterms)
    for hom, S, ed in expected:
        match = None
        for ct in produced:
            if ct.homogeneity == hom and ct.symmetry_factor == S \
                    and sympy.simplify(ct.elem_diff - ed) == 0:
                match = ct
                break
        assert match is not None, (
            f"no counterterm matched (|τ|={hom}, S={S}, F={ed}).\nProduced:\n{res.summary()}"
        )
        produced.remove(match)

    assert not produced, f"unexpected extra counterterms:\n{produced}"


def test_homogeneities_are_negative():
    res = _build()
    for ct in res.counterterms:
        assert ct.homogeneity.is_negative()


if __name__ == "__main__":
    res = _build()
    print(res.summary())
    test_gkpz_counterterms()
    test_homogeneities_are_negative()
    print("OK")
