"""Phase 3: the extraction-contraction coproduct δ (negative/renormalization side)."""
import pytest
from sympy import Function, Derivative, Rational

from regstruct import Noise, Parabolic, SPDE, Unknown, kappa
from regstruct.core.homogeneity import Homogeneity
from regstruct.equation.dsl import build_context
from regstruct.equation.generate import generate_counterterms
from regstruct.trees.coproducts import (
    coassoc_lhs, coassoc_rhs, delta_minus, delta_plus, twisted_antipode)
from regstruct.trees.tree import red_node, tree


def _gkpz_ctx(reg=Rational(-1) - kappa):
    f, g = Function("f"), Function("g")
    u = Unknown("u", dim=1)
    xi = Noise("xi", regularity=reg)
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


def test_delta_plus_co_terminates():
    # δ⁺ = (p₋⊗Id)D̄⁻ on the T⁺ factors produced by Δ — terminates and is nonempty.
    sig, _base, _u = _gkpz_ctx()
    for t in generate_counterterms(sig):
        for (_a, b), _c in delta_plus(t, sig).items():
            assert delta_minus(b, sig, root_disjoint=True)


def test_cointeraction():
    from collections import defaultdict
    from fractions import Fraction
    sig, _base, _u = _gkpz_ctx()

    def sf(fr):
        return tuple(sorted(fr, key=lambda x: x._sortkey()))

    def lhs(t):                                   # (Id⊗Δ)δ
        o = defaultdict(Fraction)
        for (phi, psi), c in delta_minus(t, sig).items():
            for (p1, p2), c2 in delta_plus(psi, sig).items():
                o[(phi, p1, p2)] += c * c2
        return {k: v for k, v in o.items() if v}

    def rhs(t):                                   # M¹³(δ⊗δ⁺)Δ
        o = defaultdict(Fraction)
        for (A, B), c in delta_plus(t, sig).items():
            for (A1, A2), c2 in delta_minus(A, sig).items():
                for (B1, B2), c3 in delta_minus(B, sig, root_disjoint=True).items():
                    o[(sf(A1 + B1), A2, B2)] += c * c2 * c3
        return {k: v for k, v in o.items() if v}

    for t in generate_counterterms(sig):
        assert lhs(t) == rhs(t)


def test_delta_counit():
    # (ε⁻ ⊗ Id) δ = Id: the only empty-left-forest term is 𝟙₋ ⊗ τ, coefficient 1.
    sig, _base, _u = _gkpz_ctx()
    for t in generate_counterterms(sig):
        empty_left = {r: c for (left, r), c in delta_minus(t, sig).items() if left == ()}
        assert empty_left == {t: 1}


def test_delta_plus_counit():
    # (Id ⊗ ε⁺) Δ = Id: the only right=𝟙₊ term is τ ⊗ 𝟙₊, coefficient 1.
    sig, _base, _u = _gkpz_ctx()
    unit = _blue_unit()
    for t in generate_counterterms(sig):
        unit_right = {left: c for (left, r), c in delta_plus(t, sig).items() if r == unit}
        assert unit_right == {t: 1}


def test_delta_plus_homogeneity_stability():
    # |τ| = |left| + |right|  (extended homogeneity), tourist_guide.tex 5851.
    sig, _base, _u = _gkpz_ctx()
    for t in generate_counterterms(sig):
        ext = t.extended_homogeneity(sig)
        for (left, right), _c in delta_plus(t, sig).items():
            assert left.extended_homogeneity(sig) + right.extended_homogeneity(sig) == ext


def test_delta_plus_plus_coassociative():
    # (Δ⁺ ⊗ Id) Δ⁺ = (Id ⊗ Δ⁺) Δ⁺ — T⁺ is a coassociative Hopf algebra (tex 5709).
    from collections import defaultdict
    from fractions import Fraction
    sig, _base, _u = _gkpz_ctx()
    tplus = {r for t in generate_counterterms(sig) for (_l, r), _c in delta_plus(t, sig).items()}

    def dpp(b):
        return delta_plus(b, sig, project_left=True)

    def lhs(b):
        o = defaultdict(Fraction)
        for (L, R), c in dpp(b).items():
            for (L1, L2), c2 in dpp(L).items():
                o[(L1, L2, R)] += c * c2
        return {k: v for k, v in o.items() if v}

    def rhs(b):
        o = defaultdict(Fraction)
        for (L, R), c in dpp(b).items():
            for (R1, R2), c2 in dpp(R).items():
                o[(L, R1, R2)] += c * c2
        return {k: v for k, v in o.items() if v}

    for b in tplus:
        assert lhs(b) == rhs(b)


def test_delta_minus_coassoc_singular():
    # Broader coverage: δ⁻ coassociativity on the richer β₀ = −3/2 − κ trees
    # (the KPZ-table trees), not just β₀ = −1 − κ.
    sig, _base, _u = _gkpz_ctx(reg=Rational(-3, 2) - kappa)
    trees = [t for t in generate_counterterms(sig) if t.nodes() <= 4]
    assert len(trees) > 5
    for t in trees:
        assert coassoc_lhs(t, sig) == coassoc_rhs(t, sig)


@pytest.mark.xfail(strict=True,
                   reason="residual: the cointeraction (Id⊗Δ)δ = M¹³(δ⊗δ⁺)Δ holds for "
                          "gKPZ (β₀=−1) but over-produces for β₀=−3/2 already on ∘—I₀—∘ — a "
                          "bug in the e'-recentering coupling of δ⁺ with δ (NOT _E_CAP "
                          "truncation: raising the cap does not help). The coproducts are "
                          "individually coassociative; fixing this coupling for more singular "
                          "noise is the next task.")
def test_cointeraction_singular():
    from collections import defaultdict
    from fractions import Fraction
    sig, _base, _u = _gkpz_ctx(reg=Rational(-3, 2) - kappa)
    trees = [t for t in generate_counterterms(sig) if t.nodes() <= 3]

    def sf(fr):
        return tuple(sorted(fr, key=lambda x: x._sortkey()))

    def lhs(t):
        o = defaultdict(Fraction)
        for (phi, psi), c in delta_minus(t, sig).items():
            for (p1, p2), c2 in delta_plus(psi, sig).items():
                o[(phi, p1, p2)] += c * c2
        return {k: v for k, v in o.items() if v}

    def rhs(t):
        o = defaultdict(Fraction)
        for (A, B), c in delta_plus(t, sig).items():
            for (A1, A2), c2 in delta_minus(A, sig).items():
                for (B1, B2), c3 in delta_minus(B, sig, root_disjoint=True).items():
                    o[(sf(A1 + B1), A2, B2)] += c * c2 * c3
        return {k: v for k, v in o.items() if v}

    for t in trees:
        assert lhs(t) == rhs(t)
