"""Validation against the paper + independent cross-checks.

The library has no reference implementation — *the paper is the oracle*. The headline
literature golden is the **full gKPZ strongly-conforming-tree table at β₀=−3/2−κ**
(tourist_guide.tex 6024–6163): matching it exactly (43 trees, six homogeneity rows)
caught a tree-generation undercount that only manifests at β₀<−1 (a node whose bare
homogeneity ≥ the bound, pulled below it by capped negative-contribution gradient
children) — this test guards against its return.

Also: an independent symmetry-factor check (brute permutation count, not the product
formula), and benchmark regression counts across the in-scope equation classes.
"""
from collections import defaultdict
from itertools import permutations
from math import factorial

import pytest
from sympy import Derivative, Function, Rational

from counterterms import Noise, Parabolic, SPDE, Unknown, daprato_lift, kappa
from counterterms.equation.dsl import build_context
from counterterms.equation.generate import generate_counterterms

f, g = Function("f"), Function("g")


def _counts(spde) -> dict:
    sig = build_context(spde)[0]
    by: dict = defaultdict(int)
    for t in generate_counterterms(sig):
        by[str(t.homogeneity(sig))] += 1
    return dict(by)


def _gkpz(beta0):
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=beta0 - kappa)
    return SPDE(noises=[xi], operator=Parabolic(dim=1, mass=1), unknown=u,
                rhs=f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2)


def test_gkpz_sc_tree_table_matches_paper_at_beta0_minus_3_2():
    # tourist_guide.tex 6024–6163 — the SC trees (no red nodes) of gKPZ at β₀=−3/2−κ.
    assert _counts(_gkpz(Rational(-3, 2))) == {
        "-3/2 - κ": 1,          # β₀
        "-1 - 2κ": 2,           # 2β₀+2
        "-1/2 - 3κ": 6,         # 3β₀+4
        "-1/2 - κ": 2,          # β₀+1
        "-2κ": 9,               # 2β₀+3   (the row whose 9th tree the bug dropped)
        "-4κ": 23,              # 4β₀+6
    }                           # total 43


def test_gkpz_renormalized_equation_at_beta0_minus_1():
    # tourist_guide.tex 6004–6012 — the five gKPZ counterterms (β₀=−1−κ).
    assert sum(_counts(_gkpz(Rational(-1))).values()) == 5


# --- independent symmetry factor: brute child-permutation count (≠ the m! formula) --- #
def _S_independent(t) -> int:
    s = 1
    for x in t.node_dec:
        s *= factorial(x)
    ch = list(t.children)
    s *= sum(1 for perm in permutations(range(len(ch)))      # |stabiliser| = Π m_j!, brute
             if all(ch[i] == ch[perm[i]] for i in range(len(ch))))
    for (_c, _p, sub) in ch:
        s *= _S_independent(sub)
    return s


@pytest.mark.parametrize("name", ["gkpz", "kpz", "gpam1", "gpam2", "system", "multinoise"])
def test_symmetry_factor_matches_independent_automorphism_count(name, request):
    from tests.conftest import CORPUS
    sig = build_context(CORPUS[name]())[0]
    for t in generate_counterterms(sig):
        assert t.symmetry_factor() == _S_independent(t)


def test_benchmark_counterterm_counts():
    # Regression counts across the in-scope classes (post bug-fix). gKPZ/KPZ rows are
    # paper-validated; PAM/gPAM/Φ⁴ are sanity-checked (subcritical, S=Aut) predictions.
    u1, u2, u3 = Unknown("u", 1), Unknown("u", 2), Unknown("u", 3)
    x1 = Noise("xi", regularity=Rational(-1) - kappa)
    x32 = Noise("xi", regularity=Rational(-3, 2) - kappa)
    cases = {
        "KPZ": SPDE(noises=[x32], operator=Parabolic(dim=1), unknown=u1,
                    rhs=x32.symbol + Derivative(u1.field, u1.x[0]) ** 2),
        "gPAM_d2": SPDE(noises=[x1], operator=Parabolic(dim=2), unknown=u2, rhs=f(u2.field) * x1.symbol),
        "PAM_d2": SPDE(noises=[x1], operator=Parabolic(dim=2), unknown=u2, rhs=u2.field * x1.symbol),
        "PAM_d3": SPDE(noises=[x32], operator=Parabolic(dim=3), unknown=u3, rhs=u3.field * x32.symbol),
    }
    expected = {"KPZ": 11, "gPAM_d2": 4, "PAM_d2": 4, "PAM_d3": 17}
    for name, spde in cases.items():
        assert sum(_counts(spde).values()) == expected[name], name


def test_phi4_via_daprato_counts():
    u2, u3 = Unknown("u", 2), Unknown("u", 3)
    x2 = Noise("xi", regularity=Rational(-2) - kappa)
    x52 = Noise("xi", regularity=Rational(-5, 2) - kappa)
    phi42 = daprato_lift(SPDE(noises=[x2], operator=Parabolic(dim=2), unknown=u2, rhs=x2.symbol - u2.field ** 3))
    phi43 = daprato_lift(SPDE(noises=[x52], operator=Parabolic(dim=3), unknown=u3, rhs=x52.symbol - u3.field ** 3))
    assert sum(_counts(phi42).values()) == 3
    assert sum(_counts(phi43).values()) == 13
