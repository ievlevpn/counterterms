"""Phase 3: the extraction-contraction coproduct δ (negative/renormalization side)."""
from sympy import Function, Derivative, Rational

from regstruct import Noise, Parabolic, SPDE, Unknown, kappa
from regstruct.core.homogeneity import Homogeneity
from regstruct.equation.dsl import build_context
from regstruct.equation.generate import generate_counterterms
from regstruct.trees.coproducts import (
    coassoc_lhs, coassoc_rhs, delta_minus, delta_plus, twisted_antipode)
from regstruct.trees.tree import red_node, tree


def _gkpz_ctx():
    f, g = Function("f"), Function("g")
    u = Unknown("u", dim=1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    spde = SPDE(noises=[xi], operator=Parabolic(dim=1, mass=1), unknown=u,
                rhs=f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2)
    return build_context(spde)


def test_delta_circ_golden():
    # tourist_guide.tex 6170:  D⁻∘ = 𝟙₋ ⊗ ∘ + ∘ ⊗ ●(β₀)
    sig, _base, _u = _gkpz_ctx()
    circ = tree("xi", (0, 0))
    beta0 = Homogeneity(-1, -1)
    expected = {
        ((), circ): 1,                                   # 𝟙₋ ⊗ ∘
        ((circ,), red_node(beta0, width=2)): 1,          # ∘ ⊗ ●(β₀)
    }
    assert delta_minus(circ, sig) == expected


def test_delta_stability_invariants():
    # The paper's Lemma: every right factor has extended homogeneity = |τ|, every
    # extracted component is divergent, and the counit term 𝟙₋ ⊗ τ appears once.
    sig, _base, _u = _gkpz_ctx()
    for t in generate_counterterms(sig):
        res = delta_minus(t, sig)
        ext = t.extended_homogeneity(sig)
        assert res.get(((), t)) == 1
        assert all(r.extended_homogeneity(sig) == ext for (_f, r) in res)
        assert all(c.extended_homogeneity(sig).is_negative()
                   for (forest, _r) in res for c in forest)


def test_delta_minus_coassociative():
    # (δ⁻ ⊗ id) δ⁻ = (id ⊗ δ⁻) δ⁻  on U⁻ — the load-bearing algebraic invariant
    # (tourist_guide.tex 5710), needs no probabilistic input.
    sig, _base, _u = _gkpz_ctx()
    for t in generate_counterterms(sig):
        assert coassoc_lhs(t, sig) == coassoc_rhs(t, sig)


def test_twisted_antipode_circ():
    # S'₋(∘) = −∘  (∘ is primitive).
    sig, _base, _u = _gkpz_ctx()
    circ = tree("xi", (0, 0))
    assert twisted_antipode(circ, sig) == {(circ,): -1}


def test_twisted_antipode_terminates_on_all():
    sig, _base, _u = _gkpz_ctx()
    for t in generate_counterterms(sig):
        result = twisted_antipode(t, sig)
        assert result          # nonempty forest-sum


def _blue_unit():
    return tree("bullet", (0, 0), (), color="blue")   # 𝟙₊


def test_delta_plus_circ():
    # Δ∘ = ∘ ⊗ 𝟙₊  (the noise does not recenter).
    sig, _base, _u = _gkpz_ctx()
    circ = tree("xi", (0, 0))
    assert delta_plus(circ, sig) == {(circ, _blue_unit()): 1}


def test_delta_plus_comodule_coassociative():
    # (Δ ⊗ id) Δ = (id ⊗ Δ⁺) Δ   (tourist_guide.tex 5708): T is a right T⁺-comodule.
    from collections import defaultdict
    from fractions import Fraction
    sig, _base, _u = _gkpz_ctx()

    def lhs(t):
        out = defaultdict(Fraction)
        for (A, B), c in delta_plus(t, sig).items():
            for (A1, A2), c2 in delta_plus(A, sig).items():
                out[(A1, A2, B)] += c * c2
        return {k: v for k, v in out.items() if v}

    def rhs(t):
        out = defaultdict(Fraction)
        for (A, B), c in delta_plus(t, sig).items():
            for (B1, B2), c2 in delta_plus(B, sig, project_left=True).items():
                out[(A, B1, B2)] += c * c2
        return {k: v for k, v in out.items() if v}

    for t in generate_counterterms(sig):
        assert lhs(t) == rhs(t)
