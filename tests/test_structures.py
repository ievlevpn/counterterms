"""The assembled `RenormalizationStructure` and the symbolic BHZ character.

`build_renormalization(spde)` bundles the validated coproducts with the negative
twisted antipode `S'₋` and returns the BHZ character ``k(τ) = h(S'₋ τ)`` as an
exact symbolic combination of the (free) analytic values ``h(σ) = 𝔼[Π^ζ σ](0)``.
The numeric `h` are out of scope (need Wick/Gaussian input, Phase 4); here we check
the *symbolic* combination the twisted antipode prescribes.
"""
import pytest

from regstruct.structures import build_regularity_structure, build_renormalization
from regstruct.trees.tree import tree

from tests.conftest import CORPUS, gkpz


def test_structure_assembles_the_divergent_generators():
    rs = build_renormalization(gkpz())
    assert len(rs.divergent) == 5
    assert all(t.homogeneity(rs.sig).is_negative() for t in rs.divergent)


def test_bhz_character_on_the_noise():
    # S'₋(∘) = −∘, so k(∘) = h(S'₋∘) = h(−∘) = −h(∘).
    rs = build_renormalization(gkpz())
    circ = tree("xi", (0, 0))
    assert rs.bhz_character(circ) == -rs.h_symbol(circ)


@pytest.mark.parametrize("name", list(CORPUS))
def test_bhz_character_is_a_polynomial_in_h(name):
    # Every character value is total and is a ℚ-polynomial in the h-symbols only.
    rs = build_renormalization(CORPUS[name]())
    assert rs.divergent
    for t in rs.divergent:
        expr = rs.bhz_character(t)
        assert expr.free_symbols <= set(rs._h.values())


# --------------------------------------------------------------------------- #
# RegularityStructure (T, T⁺) — the positive side
# --------------------------------------------------------------------------- #

def _branch_positive(blue_tree, sig):
    """Every planted branch I_p(σ) of a T⁺ (blue-rooted) tree has homogeneity > 0."""
    for (comp, p, sub) in blue_tree.children:
        h = sig.edge_gain(comp, p) + sub.extended_homogeneity(sig)
        if not (h.std > 0 or (h.std == 0 and h.kap > 0)):
            return False
    return True


def test_regularity_structure_is_graded_and_bounded_below():
    rs = build_regularity_structure(gkpz())
    assert rs.model_basis
    assert all(t.homogeneity(rs.sig).std < rs.gamma for t in rs.model_basis)
    # grades partition the basis; the homogeneity set A is sorted & bounded below
    assert sum(len(v) for v in rs.grades().values()) == len(rs.model_basis)
    A = rs.homogeneities()
    assert A and A == sorted(A, key=lambda h: h._key())


def test_divergent_subspace_matches_counterterms():
    from regstruct.equation.generate import generate_counterterms
    rs = build_regularity_structure(gkpz())
    assert set(rs.divergent) == set(generate_counterterms(rs.sig))
    assert len(rs.divergent) == 5            # the five gKPZ counterterms live inside T


@pytest.mark.parametrize("name", list(CORPUS))
def test_recentering_is_triangular_into_Tplus(name):
    # Δ : T → T ⊗ T⁺ — co-graded (|left|+|right|=|τ|), triangular (|left| ≤ |τ|),
    # and every right factor is a positive blue tree (lands in T⁺).
    rs = build_regularity_structure(CORPUS[name]())
    for t in rs.model_basis:
        ext = t.extended_homogeneity(rs.sig)
        for (left, right), _c in rs.recentering(t).items():
            lh = left.extended_homogeneity(rs.sig)
            assert lh + right.extended_homogeneity(rs.sig) == ext
            assert lh._key() <= ext._key()
            assert right.color == "blue" and _branch_positive(right, rs.sig)


def test_Tplus_generators_are_positive():
    rs = build_regularity_structure(gkpz())
    pos = rs.positive_basis()
    assert pos
    assert all(_branch_positive(b, rs.sig) for b in pos)


def test_structure_antipode_via_generic_hopf():
    from fractions import Fraction
    rs = build_regularity_structure(gkpz())
    S = rs.structure_antipode()
    unit = tree("bullet", (0,) * rs.sig.width, (), color="blue")
    assert S(unit) == {unit: Fraction(1)}          # S(𝟙₊) = 𝟙₊
    assert all(S(b) for b in rs.positive_basis())  # total & nonempty on T⁺


# --------------------------------------------------------------------------- #
# RenormalizationGroup G⁻ — the character group (group law via core/hopf)
# --------------------------------------------------------------------------- #

def test_renormalization_group_axioms():
    from fractions import Fraction
    from regstruct.structures import build_renormalization_group
    G = build_renormalization_group(gkpz())
    gens1 = [(t,) for t in G.generators]
    eps = G.identity()
    f = G.character({t: Fraction(i + 1, 2) for i, t in enumerate(G.generators)})
    g = G.character({t: Fraction(i + 2, 3) for i, t in enumerate(G.generators)})
    h = G.character({t: Fraction(i - 1, 5) for i, t in enumerate(G.generators)})

    ef, fe = G.convolve(eps, f), G.convolve(f, eps)
    assert all(ef(x) == f(x) == fe(x) for x in gens1)               # unit

    ff = G.convolve(f, G.inverse(f))
    assert all(ff(x) == eps(x) for x in gens1) and ff(()) == 1      # inverse

    lhs = G.convolve(G.convolve(f, g), h)
    rhs = G.convolve(f, G.convolve(g, h))
    assert all(lhs(x) == rhs(x) for x in gens1)                     # associativity


def test_character_is_multiplicative():
    from fractions import Fraction
    from regstruct.structures import build_renormalization_group
    G = build_renormalization_group(gkpz())
    chi = G.character({t: Fraction(i + 1) for i, t in enumerate(G.generators)})
    t0, t1 = G.generators[0], G.generators[1]
    assert chi((t0, t1)) == chi((t0,)) * chi((t1,))                 # χ(ab)=χ(a)χ(b)
    assert chi(()) == 1                                             # χ(𝟙)=1
