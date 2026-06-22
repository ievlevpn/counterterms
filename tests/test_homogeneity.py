"""Convention #2 — homogeneities live in the ordered ring ``ℚ ⊕ ℚ·κ`` (κ a positive
infinitesimal), compared lexicographically, **never as floats**.  Critical trees sit
exactly at homogeneity ``−kκ`` (tourist_guide.tex 6066), so the sign decision at the
``std == 0`` boundary is the load-bearing, silently-fatal-if-wrong operation.
"""
from fractions import Fraction

from regstruct.core.homogeneity import Homogeneity, Scaling


def test_is_negative_at_the_kappa_boundary():
    # The whole point of the κ-graded ring: distinguish 0, −κ and +κ.
    assert not Homogeneity(0, 0).is_negative()          # 0 is not negative
    assert Homogeneity(0, -1).is_negative()             # −κ  IS negative (a critical tree)
    assert not Homogeneity(0, 1).is_negative()          # +κ is positive
    assert Homogeneity(-1, +5).is_negative()            # std dominates: −1+5κ < 0
    assert not Homogeneity(1, -5).is_negative()         # std dominates: 1−5κ > 0


def test_ordering_is_lexicographic_std_then_kappa():
    # std part first, then the κ-coefficient.
    assert Homogeneity(-1, 0) < Homogeneity(0, 0)
    assert Homogeneity(0, -1) < Homogeneity(0, 0) < Homogeneity(0, 1)
    assert sorted([Homogeneity(0, 1), Homogeneity(-1, 5), Homogeneity(0, -2)]) == \
        [Homogeneity(-1, 5), Homogeneity(0, -2), Homogeneity(0, 1)]


def test_addition_and_negation_are_componentwise():
    assert Homogeneity(1, 2) + Homogeneity(Fraction(1, 2), -3) == Homogeneity(Fraction(3, 2), -1)
    assert -Homogeneity(2, -1) == Homogeneity(-2, 1)
    assert Homogeneity(2, -1) + (-Homogeneity(2, -1)) == Homogeneity(0, 0)


def test_value_is_normalised_to_fractions_for_exact_equality():
    # ints / floats-as-rationals collapse to one canonical value (frozen, hashable).
    assert Homogeneity(1) == Homogeneity(Fraction(1), Fraction(0))
    assert hash(Homogeneity(1, 2)) == hash(Homogeneity(1, 2))
    assert len({Homogeneity(0, -1), Homogeneity(0, -1)}) == 1


def test_str_round_trips_the_three_shapes():
    assert str(Homogeneity(0, 0)) == "0"
    assert str(Homogeneity(Fraction(-3, 2), 0)) == "-3/2"
    assert str(Homogeneity(0, -1)) == "-κ"
    assert str(Homogeneity(-1, -1)) == "-1 - κ"
    assert str(Homogeneity(0, -2)) == "-2κ"


def test_scaling_scaled_degree():
    # parabolic scaling 𝔰=(2,1): |∂_t|=2, |∂_x|=1, |∂_x²|=2.
    s = Scaling((2, 1))
    assert s.scaled((0, 0)) == 0
    assert s.scaled((0, 1)) == 1
    assert s.scaled((1, 0)) == 2
    assert s.scaled((0, 2)) == 2
