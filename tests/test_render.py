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


def test_canonical_bphz_section():
    eq = _build()
    txt = eq.report(canonical=True)
    # centered-Gaussian parity: the bare noise ∘ (one noise vertex) vanishes
    assert "vanishes" in txt
    import json
    data = json.loads(eq.render("json", canonical=True))
    rows = {r["tree"]: r for r in data["canonical_bphz"]}
    assert rows["Ξ"]["vanishes"] and rows["Ξ"]["value"] == "0"   # odd parity
    # gKPZ: exactly two canonical constants survive (3 of 5 vanish)
    assert sum(not r["vanishes"] for r in data["canonical_bphz"]) == 2
    # the legend lists only h-symbols that survive in the reduced constants
    assert {h["symbol"] for h in data["h_legend"]} == {"h0", "h1"}
    # the canonically renormalized equation is printed, with the surviving counterterms
    # substituted (h-coefficients, no free k's) and the vanishing ones dropped
    assert "Canonically renormalized equation" in txt
    fam = data["canonical_family_latex"]["u"]
    assert "h_{0}" in fam and "h_{1}" in fam and "k_" not in fam


def test_canonical_shows_epsilon_regularized_integrals():
    # the surviving h_ε(σ) must be spelled out as the explicit (divergent) Wick integral,
    # with the ε-regularization made explicit (ξ_ε, C_ε), in text and LaTeX.
    eq = _build()
    txt = eq.report(canonical=True)
    assert "ξ_ε" in txt and "DIVERGE as ε→0" in txt
    # gKPZ h0 = ∫ C_ε(z1−z2)·∂^(0,1)K(−z1)·∂^(0,1)K(−z2) dz1 dz2
    assert "∫ C_ε(z1 - z2)·∂^(0,1)K(-z1)·∂^(0,1)K(-z2) dz1 dz2" in txt
    assert "h_ε(" in txt
    tex = eq.latex_document(canonical=True)
    assert r"\xi_\varepsilon" in tex and r"C_\varepsilon" in tex
    assert r"\int_{(\mathbb R^{2})^{2}}" in tex and r"\partial^{(0,1)}K" in tex


def test_forest_draws_red_contraction_node():
    # the drawer renders Phase-3 red (contraction) nodes distinctly + their o-decoration
    from fractions import Fraction
    from regstruct.core.homogeneity import Homogeneity
    from regstruct.trees.tree import red_node
    from tests.conftest import gkpz
    sig = gkpz().renormalize().sig
    tex = forest(red_node(Homogeneity(Fraction(-1)), width=sig.width), sig)
    assert "redvertex" in tex and "o{=}" in tex


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


def test_markdown_table_well_formed():
    import re
    from regstruct.render.report import _md_cell
    assert _md_cell("a|b") == "a\\|b"          # a literal pipe is escaped for GFM
    md = _build().render("markdown")
    # in the divergent-trees table every row must have the same number of *unescaped*
    # pipes (column delimiters) — a stray | in a cell would make it ragged
    rows = [ln for ln in md.splitlines() if ln.startswith("|")]
    counts = {len(re.findall(r"(?<!\\)\|", ln)) for ln in rows}
    assert len(counts) == 1, f"ragged markdown table: {counts}"


def test_save_writes_exports(tmp_path):
    eq = _build()
    written = eq.save(stem="gkpz", outdir=str(tmp_path), pdf=False)
    names = {p.name for p in written}
    assert names == {"gkpz.txt", "gkpz.md", "gkpz.json", "gkpz.tex"}
    assert all(p.exists() and p.stat().st_size > 0 for p in written)


if __name__ == "__main__":
    print(_build().report())
