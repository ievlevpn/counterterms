"""Phase 4 / B1 — the canonical (BPHZ) character, symbolic half (`renorm/scheme.py`).

Tests the Wick-pairing combinatorics, the mean-zero parity rule, and that — composed
with the twisted antipode — it correctly zeroes the canonical constants of odd-noise
trees while keeping the genuinely divergent (even-noise) ones.
"""
from math import factorial

from regstruct.equation.dsl import build_context
from regstruct.renorm.scheme import WHITE_NOISE, expectation, has_odd_noise, wick_pairings
from regstruct.structures import build_renormalization
from regstruct.trees.tree import tree

from tests.conftest import gkpz

CIRC = tree("xi", (0, 0))
SIG = build_context(gkpz())[0]


def _double_factorial_odd(m):           # (2m−1)!! = number of perfect matchings of 2m items
    return factorial(2 * m) // (2 ** m * factorial(m))


def test_wick_pairings_count():
    assert wick_pairings([]) == [[]]
    assert wick_pairings([1]) == []                       # odd ⇒ no matching
    assert wick_pairings([1, 2]) == [[(1, 2)]]
    for m in (1, 2, 3):
        items = list(range(2 * m))
        assert len(wick_pairings(items)) == _double_factorial_odd(m)


def test_expectation_parity():
    # ∘ (1 noise) and ●—I—∘ (1 noise) vanish; ∘—I—∘ and ●—2I—∘∘ (2 noises) do not.
    one = tree("bullet", (0, 0), [(0, (0, 1), CIRC)])                       # ●—I—∘, 1 noise
    two = tree("xi", (0, 0), [(0, (0, 0), CIRC)])                           # ∘—I—∘, 2 noises
    twob = tree("bullet", (0, 0), [(0, (0, 1), CIRC), (0, (0, 1), CIRC)])   # ●—2I—∘∘, 2 noises
    assert expectation(CIRC, SIG).is_zero
    assert expectation(one, SIG).is_zero
    assert not expectation(two, SIG).is_zero and len(expectation(two, SIG).terms) == 1
    assert not expectation(twob, SIG).is_zero and len(expectation(twob, SIG).terms) == 1


def test_expectation_integrand_structure():
    # E[Π(∘—I₀—∘)](0) = ∫ K_0_0(−z) C(−z) dz : one pairing, kernel × covariance.
    two = tree("xi", (0, 0), [(0, (0, 0), CIRC)])
    (integrand, variables), = expectation(two, SIG).terms
    assert len(variables) == 1
    assert "K_0_0" in str(integrand) and "C(" in str(integrand)   # kernel × covariance


def test_canonical_character_parity_on_gkpz():
    # Of the five gKPZ counterterms, the three odd-noise trees have canonical constant 0;
    # the two even-noise trees (∘—I—∘ and ●—2I—∘∘) survive.
    rs = build_renormalization(gkpz())
    zero = [t for t in rs.divergent if rs.canonical_character(t) == 0]
    nonzero = [t for t in rs.divergent if rs.canonical_character(t) != 0]
    assert len(zero) == 3 and len(nonzero) == 2
    assert all(not has_odd_noise(t, rs.sig) for t in nonzero)
    assert all(has_odd_noise(t, rs.sig) for t in zero)
