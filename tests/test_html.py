"""The HTML explorer (`render/html.py`) — notes/visual_layer_design.md Part 2."""
import re
from html.parser import HTMLParser

import pytest
from sympy import Rational

from counterterms import Noise, Parabolic, SPDE, Unknown, daprato_lift, kappa
from tests.conftest import gkpz, kpz, system


class _P(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tags = []
    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))


def _tags(doc, tag, cls):
    p = _P()
    p.feed(doc)
    return [a for (t, a) in p.tags if t == tag and cls in (a.get("class") or "").split()]


def test_html_gkpz_cards_and_links():
    doc = gkpz().renormalize().render("html")
    assert doc.lstrip().lower().startswith("<!doctype html>")
    cards, cts = _tags(doc, "article", "card"), _tags(doc, "span", "ct")
    assert len(cards) == 5 and len(cts) == 5
    assert {c["data-tree"] for c in cards} == {c["data-tree"] for c in cts}
    assert len(_tags(doc, "section", "detail")) == 5
    assert doc.count("<svg") >= 10                       # card + detail per tree


def test_html_is_offline():
    doc = gkpz().renormalize().render("html")
    urls = set(re.findall(r"https?://[^\s\"'<>]+", doc))
    assert urls <= {"http://www.w3.org/2000/svg", "http://www.w3.org/1998/Math/MathML"}


def test_html_kpz_marks_non_contributors():
    doc = kpz().renormalize().render("html")
    cards = _tags(doc, "article", "card")
    inert = [c for c in cards if "inert" in c["class"].split()]
    assert len(cards) == 11 and len(inert) == 3
    linked = {c["data-tree"] for c in _tags(doc, "span", "ct")}
    assert len(linked) == 8 and not linked & {c["data-tree"] for c in inert}
    assert doc.count("F(τ*) = 0 — no counterterm") >= 3
    assert len(_tags(doc, "section", "detail")) == 11


def test_html_groups_by_homogeneity_most_singular_first():
    eq = kpz().renormalize()
    homs = [g["data-hom"] for g in _tags(eq.render("html"), "div", "grade")]
    want = []
    for t in sorted(eq.all_trees, key=lambda t: t.homogeneity(eq.sig)._key()):
        h = str(t.homogeneity(eq.sig))
        if not want or want[-1] != h:
            want.append(h)
    assert homs == want


def test_html_system_links_each_component():
    eq = system().renormalize()
    cts = _tags(eq.render("html"), "span", "ct")
    assert len(cts) == sum(len(v) for v in eq.per_component.values())
    assert {c["data-comp"] for c in cts} == {"0", "1"}


def test_html_escapes_names():
    u = Unknown("u", 1)
    a = Noise("a<b", regularity=-1 - kappa)
    doc = SPDE(noises=[a], operator=Parabolic(dim=1), unknown=u, rhs=a.symbol).renormalize().render("html")
    assert "a<b" not in doc and "a&lt;b" in doc


def test_html_lifted_equation_renders_without_canonical():
    u = Unknown("u", 3)
    xi = Noise("xi", regularity=Rational(-5, 2) - kappa)
    eq = daprato_lift(SPDE(noises=[xi], operator=Parabolic(dim=3), unknown=u,
                           rhs=xi.symbol - u.field ** 3)).renormalize()
    doc = eq.render("html")
    assert len(_tags(doc, "article", "card")) == len(eq.all_trees)
    assert not _tags(doc, "div", "canonical")
    with pytest.raises(NotImplementedError, match="Wick powers"):
        eq.render("html", canonical=True)


def test_save_writes_html_only_when_asked(tmp_path):
    eq = gkpz().renormalize()
    default = eq.save(outdir=str(tmp_path / "a"), pdf=False)
    assert not any(p.suffix == ".html" for p in default)
    (p,) = eq.save(stem="g", outdir=str(tmp_path / "b"), formats=("html",), pdf=False)
    assert p.name == "g.html" and p.read_text(encoding="utf-8").startswith("<!doctype html>")


def test_coproduct_terms_match_delta_minus():
    from counterterms.trees.coproducts import delta_minus
    eq = gkpz().renormalize()
    doc = eq.render("html")
    want = sum(len(delta_minus(t, eq.sig)) for t in eq.all_trees)
    assert len(_tags(doc, "span", "term")) == want
    assert len(_tags(doc, "div", "coproduct")) == len(eq.all_trees)


def test_coproduct_overflow_is_kept():
    from counterterms.trees.coproducts import delta_minus
    eq = kpz().renormalize()
    sizes = [len(delta_minus(t, eq.sig)) for t in eq.all_trees]
    assert max(sizes) > 12                               # the fixture really exercises overflow
    doc = eq.render("html")
    assert len(_tags(doc, "span", "term")) == sum(sizes)
    assert doc.count("<details") == sum(1 for n in sizes if n > 12)


def test_canonical_section_and_pairing_diagrams():
    from counterterms.renorm.scheme import is_bare, pairings
    eq = gkpz().renormalize()
    assert not _tags(eq.render("html"), "div", "canonical")
    doc = eq.render("html", canonical=True)
    fam = eq.canonical()
    assert len(_tags(doc, "div", "canonical")) == len(fam.constants)
    arcs = sum(len(p) for sym, tr in fam.legend.items() if is_bare(tr)
               for v in fam.constants.values() if sym in v.free_symbols
               for p in pairings(tr, eq.sig))
    assert arcs > 0 and len(_tags(doc, "path", "pair")) == arcs
    assert "canonically renormalized" in doc


def test_canonical_reduced_kpz_renders():
    doc = kpz().renormalize().render("html", reduced=True)
    assert len(_tags(doc, "div", "canonical")) == 8
