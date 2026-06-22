"""The renormalisation coproducts and their algebraic laws.

Two layers:

* a few **golden values** pinned on the gKPZ example (`δ∘`, `Δ∘`, `S'₋∘`);
* the **structural laws**, parametrized over the whole `ctx` corpus so they are
  checked across scalar/system, one/many noises, every dimension and β₀ — not on a
  single equation:
    - counits      `(ε⁻⊗Id)δ = Id`,  `(Id⊗ε⁺)Δ = Id`
    - stability    every right factor has extended homogeneity ``|τ|``; extracted
                   components are divergent (tourist_guide.tex 5851)
    - coassoc      `δ⁻`, `Δ⁺` associative; `Δ` a `T⁺`-comodule (tex 5708–5710)
    - cointeraction `(Id⊗Δ)δ = M¹³(δ⊗δ⁺)Δ` (tex 5717) — incl. the singular β₀=−3/2 (`kpz`)

Tree size is capped per law to keep the suite fast; the bound exceeds every
divergent tree the corpus actually produces at the stated β₀.
"""
from collections import defaultdict
from fractions import Fraction

from counterterms.core.homogeneity import Homogeneity
from counterterms.equation.dsl import build_context
from counterterms.equation.generate import generate_trees
from counterterms.trees.coproducts import (
    _delta_group_forest, coassoc_lhs, coassoc_rhs, delta_minus, delta_plus, twisted_antipode)
from counterterms.trees.tree import red_node, tree

from tests.conftest import gkpz

SIG = build_context(gkpz())[0]
CIRC = tree("xi", (0, 0))


def _sf(fr):
    return tuple(sorted(fr, key=lambda x: x._sortkey()))


# --------------------------------------------------------------------------- #
# golden values (gKPZ)
# --------------------------------------------------------------------------- #

def test_delta_minus_circ_golden():
    # tourist_guide.tex 6170:  δ∘ = 𝟙₋ ⊗ ∘ + ∘ ⊗ ●(β₀)
    assert delta_minus(CIRC, SIG) == {
        ((), CIRC): 1,
        ((CIRC,), red_node(Homogeneity(-1, -1), width=2)): 1,
    }


def test_delta_plus_circ_golden():
    # Δ∘ = ∘ ⊗ 𝟙₊  (the noise does not recenter)
    blue_unit = tree("bullet", (0, 0), (), color="blue")
    assert delta_plus(CIRC, SIG) == {(CIRC, blue_unit): 1}


def test_twisted_antipode_circ():
    # S'₋(∘) = −∘  (∘ is primitive)
    assert twisted_antipode(CIRC, SIG) == {(CIRC,): -1}


# --------------------------------------------------------------------------- #
# counits
# --------------------------------------------------------------------------- #

def test_delta_minus_counit(ctx):
    """(ε⁻⊗Id)δ = Id: the only empty-left-forest term is 𝟙₋⊗τ, coefficient 1."""
    for t in ctx.trees(max_nodes=4):
        empty = {r: c for (left, r), c in delta_minus(t, ctx.sig).items() if left == ()}
        assert empty == {t: 1}


def test_delta_plus_counit(ctx):
    """(Id⊗ε⁺)Δ = Id: the only right=𝟙₊ term is τ⊗𝟙₊, coefficient 1."""
    unit = tree("bullet", (0,) * ctx.sig.width, (), color="blue")
    for t in ctx.trees(max_nodes=4):
        right_unit = {left: c for (left, r), c in delta_plus(t, ctx.sig).items() if r == unit}
        assert right_unit == {t: 1}


# --------------------------------------------------------------------------- #
# stability of homogeneity (the grading the whole theory rests on)
# --------------------------------------------------------------------------- #

def test_delta_minus_stability(ctx):
    """Every δ right factor has extended homogeneity |τ|; every extracted component
    is divergent; the counit term 𝟙₋⊗τ appears exactly once."""
    for t in ctx.trees(max_nodes=4):
        res = delta_minus(t, ctx.sig)
        ext = t.extended_homogeneity(ctx.sig)
        assert res.get(((), t)) == 1
        assert all(r.extended_homogeneity(ctx.sig) == ext for (_f, r) in res)
        assert all(c.extended_homogeneity(ctx.sig).is_negative()
                   for (forest, _r) in res for c in forest)


def test_delta_plus_stability(ctx):
    """|τ| = |left| + |right| for every Δ term (tourist_guide.tex 5851)."""
    for t in ctx.trees(max_nodes=4):
        ext = t.extended_homogeneity(ctx.sig)
        for (left, right), _c in delta_plus(t, ctx.sig).items():
            assert (left.extended_homogeneity(ctx.sig)
                    + right.extended_homogeneity(ctx.sig)) == ext


# --------------------------------------------------------------------------- #
# coassociativities
# --------------------------------------------------------------------------- #

def test_delta_minus_coassociative(ctx):
    """(δ⁻⊗Id)δ⁻ = (Id⊗δ⁻)δ⁻ on U⁻ (tourist_guide.tex 5710) — on every divergent tree."""
    for t in ctx.trees():
        assert coassoc_lhs(t, ctx.sig) == coassoc_rhs(t, ctx.sig)


def test_delta_comodule_coassociative(ctx):
    """(Δ⊗Id)Δ = (Id⊗Δ⁺)Δ: T is a right T⁺-comodule (tourist_guide.tex 5708)."""
    def lhs(t):
        o = defaultdict(Fraction)
        for (A, B), c in delta_plus(t, ctx.sig).items():
            for (A1, A2), c2 in delta_plus(A, ctx.sig).items():
                o[(A1, A2, B)] += c * c2
        return {k: v for k, v in o.items() if v}

    def rhs(t):
        o = defaultdict(Fraction)
        for (A, B), c in delta_plus(t, ctx.sig).items():
            for (B1, B2), c2 in delta_plus(B, ctx.sig, project_left=True).items():
                o[(A, B1, B2)] += c * c2
        return {k: v for k, v in o.items() if v}

    for t in ctx.trees(max_nodes=4):
        assert lhs(t) == rhs(t)


def test_delta_plus_plus_coassociative(ctx):
    """(Δ⁺⊗Id)Δ⁺ = (Id⊗Δ⁺)Δ⁺: T⁺ is a coassociative Hopf algebra (tex 5709)."""
    sig = ctx.sig
    tplus = {r for t in ctx.trees(max_nodes=4) for (_l, r), _c in delta_plus(t, sig).items()}

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


def test_delta_plus_comodule(ctx):
    """(δ⁻⊗Id)δ⁺ = (Id⊗δ⁺)δ⁺: T⁺ is a left U⁻-comodule (compatibility condition (b),
    tourist_guide.tex 3452 / 5711) — over the T⁺ basis, incl. the singular β₀=−3/2 (`kpz`).
    Distinct from the cointeraction (c) below; guards the `coproducts.py` between-edge
    convention (the dropped "force internal under red" rule used to break this law)."""
    sig = ctx.sig

    def dplus(b):                      # δ⁺ = (p₋⊗Id)D̄⁻ : T⁺ → U⁻ ⊗ T⁺
        return delta_minus(b, sig, root_disjoint=True)

    def lhs(b):                        # (δ⁻ ⊗ Id) δ⁺
        o = defaultdict(Fraction)
        for (A, B), c in dplus(b).items():
            for (A1, A2), c2 in _delta_group_forest(A, sig).items():
                o[(A1, A2, B)] += c * c2
        return {k: v for k, v in o.items() if v}

    def rhs(b):                        # (Id ⊗ δ⁺) δ⁺
        o = defaultdict(Fraction)
        for (A, B), c in dplus(b).items():
            for (B1, B2), c2 in dplus(B).items():
                o[(A, B1, B2)] += c * c2
        return {k: v for k, v in o.items() if v}

    # T⁺ basis from the full (positive-sector) model basis, not just the divergent trees —
    # the failing elements are *planted* trees I(σ), which have positive homogeneity.
    model = [t for t in generate_trees(sig, Fraction(2)) if t.nodes() <= 4]
    tplus = {r for t in model for (_l, r), _c in delta_plus(t, sig).items()}
    for b in tplus:
        assert lhs(b) == rhs(b)


# --------------------------------------------------------------------------- #
# the cointeraction (the deepest law; includes the singular β₀=−3/2 corpus member)
# --------------------------------------------------------------------------- #

def test_cointeraction(ctx):
    """(Id⊗Δ)δ = M¹³(δ⊗δ⁺)Δ (tourist_guide.tex 5717)."""
    sig = ctx.sig

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
                    o[(_sf(A1 + B1), A2, B2)] += c * c2 * c3
        return {k: v for k, v in o.items() if v}

    for t in ctx.trees(max_nodes=4):
        assert lhs(t) == rhs(t)


def test_twisted_antipode_terminates(ctx):
    """S'₋ is total and returns a nonempty forest-sum on every divergent tree."""
    for t in ctx.trees(max_nodes=4):
        assert twisted_antipode(t, ctx.sig)
