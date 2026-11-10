"""The SVG tree drawer (`render/svg.py`) — notes/visual_layer_design.md Part 1."""
import xml.etree.ElementTree as ET
from fractions import Fraction

from counterterms import Noise, Parabolic, SPDE, Unknown, kappa, playground
from counterterms.core.homogeneity import Homogeneity
from counterterms.render import shorthand, svg
from counterterms.trees.tree import red_node
from tests.conftest import gkpz, kpz, multinoise, system

NS = "{http://www.w3.org/2000/svg}"


def _els(s, cls):
    root = ET.fromstring(s)
    return [e for e in root.iter() if cls in (e.get("class") or "").split()]


def _by_shorthand(spde):
    eq = spde.renormalize()
    return eq, {shorthand(t, eq.sig): t for t in eq.all_trees}


def test_svg_counts_match_the_tree():
    eq, _ = _by_shorthand(kpz())
    for t in eq.all_trees:
        s = svg(t, eq.sig)
        assert len(_els(s, "n")) == t.nodes()
        assert len(_els(s, "e")) == t.nodes() - 1
        assert sorted(int(e.get("data-id")) for e in _els(s, "n")) == list(range(t.nodes()))


def test_svg_gkpz_glyphs():
    eq, trees = _by_shorthand(gkpz())
    s = svg(trees["●·𝓘ₓ[Ξ]²"], eq.sig)
    assert len(_els(s, "n-bullet")) == 1 and len(_els(s, "n-noise")) == 2
    assert len(_els(s, "e-d")) == 2                      # two derivative kernels, dotted
    s = svg(trees["Ξ·𝓘[Ξ]"], eq.sig)
    assert len(_els(s, "n-noise")) == 2 and len(_els(s, "e-d")) == 0
    s = svg(trees["Xₓ·Ξ"], eq.sig)
    assert any("X" in "".join(e.itertext()) for e in _els(s, "lbl"))   # polynomial decoration


def test_svg_root_is_at_the_bottom():
    eq, trees = _by_shorthand(gkpz())
    s = svg(trees["●·𝓘ₓ[Ξ]²"], eq.sig)
    y = {int(e.get("data-id")): float(e.get("cy") or e.get("y")) for e in _els(s, "n")}
    assert y[0] > y[1] and y[0] > y[2]                   # SVG y grows downward


def test_svg_red_node_and_o_label():
    pg = playground("-3/2 - kappa")
    t = pg.tree(pg.I(red_node(Homogeneity(Fraction(-3, 2), Fraction(-1)), width=2), "x"))
    s = svg(t, pg.sig)
    assert len(_els(s, "n-red")) == 1
    assert any("o=" in "".join(e.itertext()) for e in _els(s, "lbl"))


def test_svg_distinguishes_noise_types():
    eq = multinoise().renormalize()
    classes = set()
    for t in eq.all_trees:
        for e in _els(svg(t, eq.sig), "n-noise"):
            classes |= {c for c in e.get("class").split() if c.startswith("n-noise-")}
    assert classes == {"n-noise-0", "n-noise-1"}


def test_svg_system_edges_carry_component_labels():
    eq = system().renormalize()
    t = next(t for t in eq.all_trees if t.children)
    assert any("(" in "".join(e.itertext()) for e in _els(svg(t, eq.sig), "lbl"))


def test_svg_pairing_arcs():
    pg = playground("-3/2 - kappa")
    t = pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"))
    assert len(_els(svg(t, pg.sig, pairing=((1, 2),)), "pair")) == 1
    assert len(_els(svg(t, pg.sig), "pair")) == 0


def test_svg_everything_inside_viewbox():
    eq, _ = _by_shorthand(kpz())
    pg = playground("-3/2 - kappa")
    wide = pg.tree(pg.I(pg.Xi), pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"), X="x")
    for t, sig in [(t, eq.sig) for t in eq.all_trees] + [(wide, pg.sig)]:
        root = ET.fromstring(svg(t, sig))
        x0, y0, w, h = map(float, root.get("viewBox").split())
        for e in root.iter():
            for ax, lo, hi in (("cx", x0, x0 + w), ("x", x0, x0 + w), ("x1", x0, x0 + w),
                               ("x2", x0, x0 + w), ("cy", y0, y0 + h), ("y", y0, y0 + h),
                               ("y1", y0, y0 + h), ("y2", y0, y0 + h)):
                if e.get(ax) is not None:
                    assert lo <= float(e.get(ax)) <= hi, (shorthand(t, sig), ax, e.attrib)


def test_svg_escapes_names():
    u = Unknown("u", 1)
    a, b, c, d = (Noise(n, regularity=-1 - kappa) for n in ("a<b", "c&d", "e", "f"))
    eq = SPDE(noises=[a, b, c, d], operator=Parabolic(dim=1), unknown=u,
              rhs=a.symbol + b.symbol + c.symbol + d.symbol).renormalize()
    for t in eq.all_trees:
        ET.fromstring(svg(t, eq.sig))                     # parses ⇒ escaped


def test_tree_repr_html_is_the_svg():
    eq, trees = _by_shorthand(gkpz())
    t = trees["●·𝓘ₓ[Ξ]²"]
    assert t._repr_html_() == svg(t) and t._repr_latex_().startswith("$")
