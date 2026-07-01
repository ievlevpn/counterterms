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
from counterterms.equation.generate import generate_counterterms, generate_trees

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


def test_order4_node_decoration_counterterm_not_dropped():
    """The node-decoration cap must scale with −β₀, not sit at the order-2 value 2.
    At order 4, β₀=−3−κ, the bare-noise counterterm Ξ^{(0,3)} has homogeneity −κ < 0
    and belongs to 𝓑_{<0}; the old fixed cap |n|_𝔰 ≤ 2 silently dropped it."""
    sig, _b, _u = _build(4, -3)
    trees = generate_counterterms(sig)
    decs = {t.node_dec for t in trees if t.node_type == "xi" and not t.children}
    assert (0, 3) in decs                               # was missing before the fix
    assert all(t.homogeneity(sig).is_negative() for t in trees)


def test_pool_keeps_trees_reached_through_overshooting_partial_sums():
    """The DFS must not prune on intermediate homogeneity sums: a decorated node above
    the bound can be pulled back under by several capped negative-contribution children.
    KPZ at β₀=−7/4−κ: ●^{(0,2)}I'(Ξ)² has std 2 + 2(1+β₀) = 1/2 < 1 and is
    rule-conforming, but the old `break` pruned it after the first child (2 − 3/4 ≥ 1)."""
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=Rational(-7, 4) - kappa)
    rhs = xi.symbol + Derivative(u.field, u.x[0]) ** 2
    sig, _b, _u = build_context(
        SPDE(noises=[xi], operator=Parabolic(dim=1), unknown=u, rhs=rhs))
    pool = generate_trees(sig)
    hits = [t for t in pool
            if t.node_type == "bullet" and t.node_dec == (0, 2) and len(t.children) == 2
            and all(p == (0, 1) and s.node_type == "xi" and not s.children
                    for (_c, p, s) in t.children)]
    assert len(hits) == 1                               # was 0 before the fix
    assert hits[0].homogeneity(sig).std == Rational(1, 2)


def test_nonsingular_noise_rejected():
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=Rational(1))       # β₀ ≥ 0
    with pytest.raises(ValueError, match="not singular"):
        build_context(SPDE(noises=[xi], operator=Parabolic(dim=1), unknown=u,
                           rhs=f(u.field) * xi.symbol))
