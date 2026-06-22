"""Phase 4 infrastructure — the seams ('sockets') are in place: the real pieces wired
through them, the unbuilt meat (B2 integrals, B3 G⁻_ad, A3 second basis) stubbed with
clear errors. No analysis/numerics implemented — just the contracts.
"""
import pytest

from regstruct.core.symbol import Symbol
from regstruct.renorm.scheme import (
    BPHZ, FreeConstants, UnbuiltEvaluator, WHITE_NOISE, expectation)
from regstruct.structures import build_renormalization, build_renormalization_group
from regstruct.trees.tree import tree

from tests.conftest import gkpz

CIRC = tree("xi", (0, 0))


def test_decorated_tree_satisfies_the_symbol_protocol():
    # A2 seam: a second basis only needs `homogeneity` + `canonical_key`.
    assert isinstance(CIRC, Symbol)
    assert CIRC.canonical_key() == CIRC._sortkey()


def test_free_constants_scheme_is_real():
    rs = build_renormalization(gkpz())
    chi = FreeConstants().character(rs)
    assert len(chi.values) == len(rs.divergent)
    assert all(v.is_Symbol for v in chi.values.values())          # free symbols


def test_bphz_scheme_symbolic_is_real_with_parity():
    # BPHZ.character is the parity-reduced canonical character (real, symbolic):
    # 3 of the 5 gKPZ constants vanish (odd noise).
    rs = build_renormalization(gkpz())
    chi = BPHZ().character(rs)
    assert sum(1 for v in chi.values.values() if v == 0) == 3
    assert sum(1 for v in chi.values.values() if v != 0) == 2


def test_bphz_numeric_path_is_an_unbuilt_socket():
    rs = build_renormalization(gkpz())
    with pytest.raises(NotImplementedError, match="B2"):
        BPHZ().numeric_character(rs)


def test_integral_evaluator_socket_raises():
    rs = build_renormalization(gkpz())
    even = next(t for t in rs.divergent if not expectation(t, rs.sig).is_zero)
    with pytest.raises(NotImplementedError, match="B2"):
        UnbuiltEvaluator().evaluate(expectation(even, rs.sig), noise=WHITE_NOISE)


def test_g_minus_ad_socket_raises():
    G = build_renormalization_group(gkpz())
    with pytest.raises(NotImplementedError, match="B3"):
        G.admissible()
