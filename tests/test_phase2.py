"""Phase 2: systems, multiple noises, general operator order."""
import warnings

import sympy
from sympy import Derivative, Function, Rational

from regstruct import Noise, Parabolic, SPDE, Unknown, jet, kappa
from regstruct.core.homogeneity import Homogeneity

f, g, h, a, b = (Function(s) for s in "fghab")


def _gkpz_expected(comp):
    """The five gKPZ counterterms (homog, S, F) in component `comp`'s jets."""
    U0, U1 = jet(comp, (0, 0)), jet(comp, (0, 1))
    return [
        (Homogeneity(-1, -1), 1, f(U0)),
        (Homogeneity(0, -1), 1, U1 * Derivative(f(U0), U0)),
        (Homogeneity(0, -1), 1, 2 * f(U0) * g(U0) * U1),
        (Homogeneity(0, -2), 1, f(U0) * Derivative(f(U0), U0)),
        (Homogeneity(0, -2), 2, 2 * f(U0) ** 2 * g(U0)),
    ]


def _match(counterterms, expected):
    produced = list(counterterms)
    assert len(produced) == len(expected)
    for hom, S, ed in expected:
        m = next((c for c in produced if c.homogeneity == hom and c.symmetry_factor == S
                  and sympy.simplify(c.elem_diff - ed) == 0), None)
        assert m is not None, f"missing (|τ|={hom}, S={S}, F={ed})"
        produced.remove(m)
    assert not produced


def test_decoupled_system_reduces_to_scalar():
    # u driven by ξ, v driven by η, no coupling ⇒ each component = its own gKPZ.
    u, v = Unknown("u", 1), Unknown("v", 1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    eta = Noise("eta", regularity=Rational(-1) - kappa)
    op = Parabolic(dim=1, mass=1)
    res = SPDE(noises=[xi, eta], equations=[
        (u, op, f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2),
        (v, op, f(v.field) * eta.symbol + g(v.field) * Derivative(v.field, v.x[0]) ** 2),
    ]).renormalize()
    assert res.n_components == 2
    _match(res.per_component[0], _gkpz_expected(0))
    _match(res.per_component[1], _gkpz_expected(1))
    # decoupled ⇒ the two components share no renormalization constant
    cs0 = {c.constant for c in res.per_component[0]}
    cs1 = {c.constant for c in res.per_component[1]}
    assert cs0.isdisjoint(cs1)


def test_coupled_shares_constants():
    # u̇ = a(v) ξ, v̇ = b(u) ξ (shared noise) — the bare ∘ tree is one object,
    # so its constant is shared, weighting a(v) in eq u and b(u) in eq v.
    u, v = Unknown("u", 1), Unknown("v", 1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    op = Parabolic(dim=1, mass=1)
    res = SPDE(noises=[xi], equations=[
        (u, op, a(v.field) * xi.symbol),
        (v, op, b(u.field) * xi.symbol),
    ]).renormalize()

    def bare(cts):
        return next(c for c in cts if c.tree.nodes() == 1 and not c.tree.children
                    and not any(c.tree.node_dec))

    bu, bv = bare(res.per_component[0]), bare(res.per_component[1])
    assert bu.constant == bv.constant                          # shared k_τ
    assert sympy.simplify(bu.elem_diff - a(jet(1, (0, 0)))) == 0   # a(v) in eq u
    assert sympy.simplify(bv.elem_diff - b(jet(0, (0, 0)))) == 0   # b(u) in eq v


def test_multiple_noises():
    # (∂_t−Δ+1)u = f(u)ξ + h(u)η ⇒ 8 counterterms; both bare primitives present.
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    eta = Noise("eta", regularity=Rational(-1) - kappa)
    res = SPDE(noises=[xi, eta], operator=Parabolic(dim=1, mass=1), unknown=u,
               rhs=f(u.field) * xi.symbol + h(u.field) * eta.symbol).renormalize()
    assert len(res.counterterms) == 8
    U0 = jet(0, (0, 0))
    bare = {sympy.simplify(c.elem_diff) for c in res.counterterms if c.tree.nodes() == 1}
    assert f(U0) in bare and h(U0) in bare


def test_operator_order_changes_homogeneities():
    # gPAM d=1, f(u)ξ. The ∘—I₀—∘ tree sits at 2β₀+m: negative for m=2, positive for m=4.
    def homs(order):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = Unknown("u", 1)
            xi = Noise("xi", regularity=Rational(-1) - kappa)
            res = SPDE(noises=[xi], operator=Parabolic(dim=1, order=order), unknown=u,
                       rhs=f(u.field) * xi.symbol).renormalize()
        return {str(c.homogeneity) for c in res.counterterms}

    assert "-2κ" in homs(2)
    assert "-2κ" not in homs(4)
    assert homs(2) != homs(4)
