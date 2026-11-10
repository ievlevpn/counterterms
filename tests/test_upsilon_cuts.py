"""The two explanatory views: the Υ-map factor by factor (`nonlinearity.node_factors`) and the
cuts behind each coproduct term (`coproducts.cuts`), plus their drawing in the explorer."""
import xml.etree.ElementTree as ET
from fractions import Fraction
from html.parser import HTMLParser

import sympy

from counterterms import playground
from counterterms.equation.dsl import build_context
from counterterms.render import shorthand, svg
from counterterms.renorm.nonlinearity import elem_diff, node_factors
from counterterms.trees.coproducts import cuts, delta_minus
from tests.conftest import gkpz, kpz, multinoise, system


class _P(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tags = []
    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))


def _count(doc, tag, cls):
    p = _P()
    p.feed(doc)
    return sum(1 for (t, a) in p.tags if t == tag and cls in (a.get("class") or "").split())


def _els(s, cls):
    return [e for e in ET.fromstring(s).iter() if cls in (e.get("class") or "").split()]


# ---- Υ-map, node by node ---------------------------------------------------- #

def test_node_factors_multiply_to_the_elementary_differential():
    for spde in (gkpz(), kpz(), system(), multinoise()):
        eq = spde.renormalize()
        sig, base, _ = build_context(spde)
        for t in eq.all_trees:
            for a in range(sig.n_components):
                fs = node_factors(t, a, base, sig)
                assert len(fs) == t.nodes()
                assert sympy.expand(sympy.Mul(*fs)) == elem_diff(t, a, base, sig)


def test_node_factors_are_in_drawing_order_and_show_the_vanishing_node():
    spde = kpz()
    eq = spde.renormalize()
    sig, base, _ = build_context(spde)
    by = {shorthand(t, sig): t for t in eq.all_trees}
    # ●·𝓘ₓ[Ξ]²: root g=(∂ₓu)² differentiated twice in uₓ → 2; each noise leaf → its f = 1
    assert node_factors(by["●·𝓘ₓ[Ξ]²"], 0, base, sig) == [2, 1, 1]
    # Xₓ·Ξ: D_x of the constant noise coefficient is 0 — the node that kills the counterterm
    assert node_factors(by["Xₓ·Ξ"], 0, base, sig) == [0]


# ---- the cuts behind a coproduct term ---------------------------------------- #

def test_cuts_account_for_every_coefficient():
    eq = kpz().renormalize()
    for t in eq.all_trees:
        d, cs = delta_minus(t, eq.sig), cuts(t, eq.sig)
        assert set(d) <= set(cs)
        for key, coeff in d.items():
            assert sum(c for c, _blocks in cs[key]) == coeff
            for _c, blocks in cs[key]:
                flat = [i for b in blocks for i in b]
                assert len(flat) == len(set(flat)) and set(flat) <= set(range(t.nodes()))
        assert cs[((), t)] == [(Fraction(1), ())]                 # the unit term: nothing cut


def test_a_coefficient_two_is_two_symmetric_cuts():
    pg = playground("-3/2 - kappa")
    t = pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"))              # nodes: 0 root, 1 and 2 noises
    d, cs = delta_minus(t, pg.sig), cuts(t, pg.sig)
    key = next(k for k in d if k[0] == (pg.Xi,))                  # extract one noise leaf
    assert d[key] == 2
    assert sorted(blocks for _c, blocks in cs[key]) == [((1,),), ((2,),)]


def test_tracing_does_not_change_the_coproduct():
    eq = kpz().renormalize()
    t = max(eq.all_trees, key=lambda t: t.nodes())
    assert delta_minus(t, eq.sig, trace={}) == delta_minus(t, eq.sig)


# ---- drawing ------------------------------------------------------------------ #

def test_svg_marks_cut_blocks_and_numbers_nodes():
    pg = playground("-3/2 - kappa")
    t = pg.tree(pg.I(pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x")), "x"), pg.I(pg.Xi, "x"))
    plain = svg(t, pg.sig)
    s = svg(t, pg.sig, mark={1: 0, 2: 0, 4: 1}, number=True)
    assert len(_els(s, "n")) == len(_els(plain, "n")) and len(_els(s, "e")) == len(_els(plain, "e"))
    assert len(_els(s, "cut-node")) == 3 and len(_els(s, "cut-edge")) == 1   # edge 1–2 is internal
    assert {e.get("class").split()[-1] for e in _els(s, "cut-node")} == {"cut-0", "cut-1"}
    assert sorted(e.text for e in _els(s, "num")) == [str(i) for i in range(t.nodes())]
    assert not _els(plain, "cut-node") and not _els(plain, "num")


def test_explorer_shows_the_upsilon_factors_per_node():
    eq = kpz().renormalize()
    doc = eq.render("html")
    assert _count(doc, "div", "upsilon") == len(eq.all_trees)
    assert _count(doc, "tr", "factor") == sum(t.nodes() for t in eq.all_trees)
    assert doc.count("this factor vanishes") == 3                 # KPZ's three F(τ*) = 0 trees
    sys_eq = system().renormalize()
    assert _count(sys_eq.render("html"), "tr", "factor") == 2 * sum(t.nodes() for t in sys_eq.all_trees)


def test_explorer_draws_the_cuts_behind_each_term():
    eq = gkpz().renormalize()
    doc = eq.render("html")
    want = 0
    for t in eq.all_trees:
        d, cs = delta_minus(t, eq.sig), cuts(t, eq.sig)
        for key in d:
            configs = {blocks for _c, blocks in cs[key] if blocks}
            want += min(4, len(configs))
    assert want > 0 and _count(doc, "span", "cutview") == want
