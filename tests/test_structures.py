"""The assembled `RenormalizationStructure` and the symbolic BHZ character.

`build_renormalization(spde)` bundles the validated coproducts with the negative
twisted antipode `S'₋` and returns the BHZ character ``k(τ) = h(S'₋ τ)`` as an
exact symbolic combination of the (free) analytic values ``h(σ) = 𝔼[Π^ζ σ](0)``.
The numeric `h` are out of scope (need Wick/Gaussian input, Phase 4); here we check
the *symbolic* combination the twisted antipode prescribes.
"""
import pytest

from regstruct.structures import build_renormalization
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
