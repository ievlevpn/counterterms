"""Notebook display (`render/display.py`): algebra results, the equation and the canonical
family draw themselves; a lone tree can be drawn with its equation's signature."""
from html.parser import HTMLParser

from counterterms import (Noise, Parabolic, SPDE, Unknown, build_regularity_structure,
                          build_renormalization, kappa, playground)
from counterterms.render import shorthand
from counterterms.trees.coproducts import delta_minus, delta_minus_group, delta_plus, twisted_antipode
from tests.conftest import gkpz, kpz, multinoise


class _P(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tags = []
    def handle_starttag(self, tag, attrs):
        self.tags.append((tag, dict(attrs)))


def _count(doc, tag, cls=None):
    p = _P()
    p.feed(doc)
    return sum(1 for (t, a) in p.tags
               if t == tag and (cls is None or cls in (a.get("class") or "").split()))


def _two():
    pg = playground("-3/2 - kappa")
    return pg, pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"))


def test_playground_sums_are_still_the_same_dicts():
    pg, t = _two()
    assert pg.delta_minus(t) == delta_minus(t, pg.sig) and isinstance(pg.delta_minus(t), dict)
    assert pg.delta_plus(t) == delta_plus(t, pg.sig)
    assert pg.antipode(t) == twisted_antipode(t, pg.sig)


def test_sums_draw_one_term_per_entry():
    pg, t = _two()
    for s in (pg.delta_minus(t), pg.delta_plus(t), pg.antipode(t)):
        html = s._repr_html_()
        assert _count(html, "span", "rs-term") == len(s) > 0
        assert "<script" not in html and "<style" not in html      # notebook-safe: inline only
    assert "⊗" in pg.delta_minus(t)._repr_html_()
    assert "⊗" not in pg.antipode(t)._repr_html_()                 # S'₋τ is a sum of forests


def test_sums_print_as_readable_text():
    pg, t = _two()
    txt = str(pg.delta_minus(t))
    assert "⊗" in txt and "●·𝓘ₓ[Ξ]²" in txt and "{" not in txt and repr(pg.delta_minus(t)) == txt
    assert str(pg.antipode(pg.Xi)) == "−Ξ"                         # S'₋(∘) = −∘


def test_structure_methods_return_drawable_sums():
    spde = gkpz()
    rn, rs = build_renormalization(spde), build_regularity_structure(spde)
    t = max(rn.divergent, key=lambda t: t.nodes())
    for got, want in ((rn.coaction(t), delta_minus(t, rn.sig)),
                      (rn.coproduct(t), delta_minus_group(t, rn.sig)),
                      (rn.twisted_antipode(t), twisted_antipode(t, rn.sig)),
                      (rs.recentering(t), delta_plus(t, rs.sig)),
                      (rs.structure_coproduct(t), delta_plus(t, rs.sig, project_left=True))):
        assert got == want and _count(got._repr_html_(), "span", "rs-term") == len(want)


def test_long_sums_are_capped_and_say_so():
    rn = build_renormalization(kpz())
    s = max((rn.twisted_antipode(t) for t in rn.divergent), key=len)
    assert len(s) > 40
    html = s._repr_html_()
    assert _count(html, "span", "rs-term") == 40 and f"{len(s) - 40} more terms" in html


def test_show_draws_a_tree_with_the_equation_signature():
    eq = multinoise().renormalize()
    t = next(t for t in eq.all_trees if shorthand(t, eq.sig) == "Ξ_eta")
    assert "n-noise-1" in eq.show(t)._repr_html_()        # eta is the second noise of the equation
    assert "n-noise-0" in t._repr_html_()                 # alone, the tree cannot know that
    pg = playground(noises={"xi": "-1 - kappa", "eta": "-1 - kappa"})
    assert "n-noise-1" in pg.show(pg.noise("eta"))._repr_html_()


def test_equation_displays_as_a_drawn_summary():
    eq = kpz().renormalize()
    html = eq._repr_html_()
    assert _count(html, "svg") == len(eq.all_trees) == 11
    assert _count(html, "math") >= 1 + 8                  # the family + one constant per card
    assert "<script" not in html and "<style" not in html
    assert eq._repr_markdown_().startswith("# Renormalized equation")   # still available


def test_explore_embeds_the_explorer_in_an_iframe():
    eq = gkpz().renormalize()
    html = eq.explore()._repr_html_()
    assert html.startswith("<iframe") and "srcdoc=" in html and "&lt;!doctype html&gt;" in html
    assert _count(html, "iframe") == 1 and _count(html, "article") == 0   # the page is escaped
    assert "canonically renormalized" in eq.explore(reduced=True)._repr_html_()


def test_canonical_family_displays_its_trees():
    fam = gkpz().renormalize().canonical()
    html = fam._repr_html_()
    assert _count(html, "tr", "rs-const") == len(fam.constants) == 5
    assert _count(html, "tr", "rs-h") == len(fam.legend)
    assert _count(html, "svg") == len(fam.constants) + len(fam.legend)


def test_display_escapes_names():
    u = Unknown("u", 1)
    a, b = Noise("a<b", regularity=-1 - kappa), Noise("c&d", regularity=-1 - kappa)
    eq = SPDE(noises=[a, b], operator=Parabolic(dim=1), unknown=u,
              rhs=a.symbol + b.symbol).renormalize()
    rn = build_renormalization(eq.spde)
    docs = [eq._repr_html_(), eq.canonical()._repr_html_()]
    docs += [rn.coaction(t)._repr_html_() for t in rn.divergent]
    assert all("a<b" not in d and "c&d" not in d for d in docs)
