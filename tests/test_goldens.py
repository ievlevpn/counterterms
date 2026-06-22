"""End-to-end goldens: ``SPDE.renormalize()`` must reproduce the paper's exact
counterterm tables.  These are the oracle (the paper is the only reference) and the
backbone regression for the whole pipeline.

* gKPZ — the exact five counterterms with their ``(|τ|, S(τ), F(τ*))`` triples
  (tourist_guide.tex 5996–6012), incl. the factor-2 (S=1) and S=2 cases (conv #3).
* KPZ — the negative-homogeneity rows of the table (tex 6028–6063).
* gPAM (d=2), systems (decoupled & coupled, conv #7), multiple noises, operator order.
"""
import warnings

from sympy import Derivative, Rational, simplify

from counterterms import Noise, Parabolic, SPDE, Unknown, jet, kappa
from counterterms.core.homogeneity import Homogeneity

from tests.conftest import a_coef, b_coef, f, g, gkpz, gpam2, h, multinoise

U0, U1 = jet(0, (0, 0)), jet(0, (0, 1))


def _match(counterterms, expected):
    """Assert the produced counterterms equal `expected` as a multiset of
    (homogeneity, symmetry factor, elementary differential) triples."""
    produced = list(counterterms)
    assert len(produced) == len(expected), \
        f"expected {len(expected)} counterterms, got {len(produced)}"
    for hom, S, ed in expected:
        m = next((c for c in produced if c.homogeneity == hom and c.symmetry_factor == S
                  and simplify(c.elem_diff - ed) == 0), None)
        assert m is not None, f"no counterterm matched (|τ|={hom}, S={S}, F={ed})"
        produced.remove(m)
    assert not produced, f"unexpected extra counterterms: {produced}"


def test_gkpz_five_counterterms():
    res = gkpz().renormalize()
    _match(res.counterterms, [
        (Homogeneity(-1, -1), 1, f(U0)),                          # ∘
        (Homogeneity(0, -1), 1, U1 * Derivative(f(U0), U0)),      # ∘1
        (Homogeneity(0, -1), 1, 2 * f(U0) * g(U0) * U1),          # ●—I_{(0,1)}—∘  (factor 2)
        (Homogeneity(0, -2), 1, f(U0) * Derivative(f(U0), U0)),   # ∘—I—∘
        (Homogeneity(0, -2), 2, 2 * f(U0) ** 2 * g(U0)),          # ●—two I_{(0,1)}—∘∘  (S=2)
    ])
    assert all(c.homogeneity.is_negative() for c in res.counterterms)


def test_kpz_negative_rows():
    # KPZ: (∂_t−Δ)u = (∂ₓu)² + ξ, space-time WN, β₀=−3/2−κ (tourist_guide.tex 6028–6063).
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=Rational(-3, 2) - kappa)
    res = SPDE(operator=Parabolic(dim=1), unknown=u, noises=[xi],
               rhs=1 * xi.symbol + Derivative(u.field, u.x[0]) ** 2).renormalize()
    assert {c.homogeneity for c in res.counterterms} == {
        Homogeneity(Rational(-3, 2), -1),   # β₀
        Homogeneity(-1, -2),                # 2β₀+2
        Homogeneity(Rational(-1, 2), -3),   # 3β₀+4
        Homogeneity(Rational(-1, 2), -1),   # β₀+1
        Homogeneity(0, -2),                 # 2β₀+3
        Homogeneity(0, -4),                 # 4β₀+6
    }


def test_gpam_d2_four_counterterms():
    res = gpam2().renormalize()
    assert len(res.counterterms) == 4
    assert {c.homogeneity for c in res.counterterms} == {
        Homogeneity(-1, -1), Homogeneity(0, -1), Homogeneity(0, -2)}


def test_decoupled_system_reduces_to_scalar():
    # u driven by ξ, v by η, no coupling ⇒ each component is its own gKPZ, with
    # disjoint renormalization constants.
    u, v = Unknown("u", 1), Unknown("v", 1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    eta = Noise("eta", regularity=Rational(-1) - kappa)
    op = Parabolic(dim=1, mass=1)
    res = SPDE(noises=[xi, eta], equations=[
        (u, op, f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2),
        (v, op, f(v.field) * eta.symbol + g(v.field) * Derivative(v.field, v.x[0]) ** 2),
    ]).renormalize()
    assert res.n_components == 2

    def gkpz_triples(comp):
        c0, c1 = jet(comp, (0, 0)), jet(comp, (0, 1))
        return [
            (Homogeneity(-1, -1), 1, f(c0)),
            (Homogeneity(0, -1), 1, c1 * Derivative(f(c0), c0)),
            (Homogeneity(0, -1), 1, 2 * f(c0) * g(c0) * c1),
            (Homogeneity(0, -2), 1, f(c0) * Derivative(f(c0), c0)),
            (Homogeneity(0, -2), 2, 2 * f(c0) ** 2 * g(c0)),
        ]
    _match(res.per_component[0], gkpz_triples(0))
    _match(res.per_component[1], gkpz_triples(1))
    assert {c.constant for c in res.per_component[0]}.isdisjoint(
        {c.constant for c in res.per_component[1]})


def test_coupled_system_shares_constants():
    # u̇=a(v)ξ, v̇=b(u)ξ (shared noise): the bare ∘ tree is one object, so its
    # constant is shared — weighting a(v) in eq u and b(u) in eq v (conv #7).
    u, v = Unknown("u", 1), Unknown("v", 1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    op = Parabolic(dim=1, mass=1)
    res = SPDE(noises=[xi], equations=[
        (u, op, a_coef(v.field) * xi.symbol),
        (v, op, b_coef(u.field) * xi.symbol),
    ]).renormalize()

    def bare(cts):
        return next(c for c in cts if c.tree.nodes() == 1 and not c.tree.children
                    and not any(c.tree.node_dec))
    bu, bv = bare(res.per_component[0]), bare(res.per_component[1])
    assert bu.constant == bv.constant                                  # shared k_τ
    assert simplify(bu.elem_diff - a_coef(jet(1, (0, 0)))) == 0        # a(v) in eq u
    assert simplify(bv.elem_diff - b_coef(jet(0, (0, 0)))) == 0        # b(u) in eq v


def test_multiple_noises_eight_counterterms():
    res = multinoise().renormalize()
    assert len(res.counterterms) == 8
    bare = {simplify(c.elem_diff) for c in res.counterterms if c.tree.nodes() == 1}
    assert f(U0) in bare and h(U0) in bare           # both bare primitives present


def test_operator_order_changes_homogeneities():
    # gPAM d=1, f(u)ξ: the ∘—I₀—∘ tree sits at 2β₀+m — negative for m=2, positive for m=4.
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
