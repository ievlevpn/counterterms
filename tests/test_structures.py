"""Phase 3: the assembled renormalization structure + symbolic BHZ character."""
from sympy import Function, Derivative, Rational

from regstruct import Noise, Parabolic, SPDE, Unknown, kappa
from regstruct.structures import build_renormalization
from regstruct.trees.tree import tree


def _gkpz():
    f, g = Function("f"), Function("g")
    u = Unknown("u", dim=1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    return SPDE(noises=[xi], operator=Parabolic(dim=1, mass=1), unknown=u,
                rhs=f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2)


def test_renormalization_structure_assembles():
    rs = build_renormalization(_gkpz())
    # the gKPZ renormalization structure has the five negative-homogeneity generators
    assert len(rs.divergent) == 5
    assert all(t.homogeneity(rs.sig).is_negative() for t in rs.divergent)


def test_bhz_character_on_noise():
    # S'₋(∘) = −∘, so k(∘) = h(−∘) = −h(∘) = −h0.
    rs = build_renormalization(_gkpz())
    circ = tree("xi", (0, 0))
    assert rs.bhz_character(circ) == -rs.h_symbol(circ)


def test_bhz_character_is_linear_in_h():
    # Every BHZ character value is a polynomial in the h-symbols with rational coeffs.
    rs = build_renormalization(_gkpz())
    for t in rs.divergent:
        expr = rs.bhz_character(t)
        assert expr.free_symbols <= set(rs._h.values())
