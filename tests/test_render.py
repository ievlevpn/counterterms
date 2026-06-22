"""Rendering backbone (notes/output.md): the gKPZ report must show the five
divergent trees, in homogeneity order, across every output format."""
import json

from regstruct.render import ascii_art, forest, shorthand
from tests.conftest import gkpz


def _build():
    """The gKPZ renormalized equation — the rendering fixture."""
    return gkpz().renormalize()

# The five gKPZ trees in canonical shorthand (notes/output.md §3a).
GKPZ_SHORTHANDS = {"Ξ", "Xₓ·Ξ", "●·𝓘ₓ[Ξ]", "Ξ·𝓘[Ξ]", "●·𝓘ₓ[Ξ]²"}


def test_text_report_has_all_trees():
    eq = _build()
    txt = eq.report()
    for sh in GKPZ_SHORTHANDS:
        assert f"τ = {sh} " in txt, f"missing tree {sh}\n{txt}"
    # exact-5 sanity + the placeholder frontier is honest
    assert "5 trees, 5 counterterms, 0 dropped" in txt
    assert "CANONICAL" in txt and "Phase 4" in txt
    assert "canonical=True" in txt           # default report points to the opt-in section


def test_canonical_bhz_section():
    eq = _build()
    txt = eq.report(canonical=True)
    # the bare-noise counterterm renormalizes to k = -h(∘)  (S'₋∘ = -∘)
    assert "k_0 = -h0" in txt
    # the h-legend distinguishes a contracted (red) node from its black twin
    assert "[contracted node, o=" in txt
    import json
    data = json.loads(eq.render("json", canonical=True))
    assert data["bhz"][0]["value"] == "-h0"
    assert any(h["contracted"] for h in data["h_legend"])


def test_trees_sorted_by_homogeneity():
    eq = _build()
    keys = [t.homogeneity(eq.sig)._key() for t in
            sorted(eq.all_trees, key=lambda t: t.homogeneity(eq.sig)._key())]
    assert keys == sorted(keys)
    # most-divergent first: the bare noise Ξ (std −1) leads
    first = sorted(eq.all_trees, key=lambda t: t.homogeneity(eq.sig)._key())[0]
    assert shorthand(first, eq.sig) == "Ξ"


def test_single_tree_drawers():
    eq = _build()
    sig = eq.sig
    shapes = {shorthand(t, sig) for t in eq.all_trees}
    assert shapes == GKPZ_SHORTHANDS
    # the S=2 tree draws two child edges
    t2 = next(t for t in eq.all_trees if shorthand(t, sig) == "●·𝓘ₓ[Ξ]²")
    art = ascii_art(t2, sig)
    assert art.count("▸Ξ") == 2 and art.splitlines()[0] == "●"
    assert "\\begin{forest}" in forest(t2, sig)


def test_json_and_other_formats():
    eq = _build()
    data = json.loads(eq.render("json"))
    assert len(data["trees"]) == 5
    assert {t["shorthand"] for t in data["trees"]} == GKPZ_SHORTHANDS
    assert all(t["contributes"] for t in data["trees"])
    assert eq.render("markdown").startswith("# Renormalized equation")
    assert eq.latex_document().rstrip().endswith(r"\end{document}")


def test_save_writes_exports(tmp_path):
    eq = _build()
    written = eq.save(stem="gkpz", outdir=str(tmp_path), pdf=False)
    names = {p.name for p in written}
    assert names == {"gkpz.txt", "gkpz.md", "gkpz.json", "gkpz.tex"}
    assert all(p.exists() and p.stat().st_size > 0 for p in written)


if __name__ == "__main__":
    print(_build().report())
