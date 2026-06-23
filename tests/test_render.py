"""Rendering backbone (notes/output.md): the gKPZ report must show the five
divergent trees, in homogeneity order, across every output format."""
import json

from counterterms.render import ascii_art, forest, shorthand
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
    assert "CANONICAL" in txt and "left symbolic" in txt
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
    assert "ξ_ε" in txt and "diverge as ε→0" in txt
    # gKPZ h0 = ∫ C_ε(z1−z2)·∂^(0,1)K(−z1)·∂^(0,1)K(−z2) dz1 dz2
    assert "∫ C_ε(z1 - z2)·∂^(0,1)K(-z1)·∂^(0,1)K(-z2) dz1 dz2" in txt
    assert "h_ε(" in txt
    tex = eq.latex_document(canonical=True)
    assert r"\xi_\varepsilon" in tex and r"C_\varepsilon" in tex
    assert r"\int_{(\mathbb R^{2})^{2}}" in tex and r"\partial^{(0,1)}K" in tex


def test_forest_draws_red_contraction_node():
    # the drawer renders Phase-3 red (contraction) nodes distinctly + their o-decoration
    from fractions import Fraction
    from counterterms.core.homogeneity import Homogeneity
    from counterterms.trees.tree import red_node
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
    from counterterms.render.report import _md_cell
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


def test_canonical_legend_marks_zeros_and_duplicates_display_only():
    """The KPZ canonical legend *marks* provably-zero and duplicate h(σ) — but the marking is
    display-only: the symbolic constants are NOT reduced (a zero-valued symbol like `h7` still
    appears in a constant), so no algebra changed."""
    from tests.conftest import kpz
    md = kpz().renormalize().render("markdown", canonical=True)
    assert "pure-kernel total derivative" in md          # provable-zero mark on the legend
    assert "same value" in md                            # duplicate mark on the legend
    # display-only: a provably-zero symbol still appears in the (unreduced) constants
    consts = [ln for ln in md.splitlines() if ln.startswith("- `k_")]
    assert any("h7" in ln or "h2" in ln for ln in consts), \
        "constants must stay unreduced (display-only) — a zero h-symbol should still appear"


def test_reduced_drops_kpz_drift_terms():
    """`reduced=True` applies the spatial-reflection reduction (symmetric noise): KPZ's two ∂h drift
    counterterms — nonzero in the general (raw) canonical character — vanish, collapsing KPZ to a
    single diverging constant (Hairer). Strictly more constants vanish than in the raw view."""
    import json
    from tests.conftest import kpz
    eq = kpz().renormalize()
    raw = {r["tree"]: r["vanishes"] for r in json.loads(eq.render("json", canonical=True))["canonical_bphz"]}
    red = {r["tree"]: r["vanishes"] for r in json.loads(eq.render("json", reduced=True))["canonical_bphz"]}
    for drift in ("●·𝓘ₓ[●·𝓘ₓ[Ξ]]·𝓘ₓ[Ξ]", "●·𝓘ₓ[●·𝓘ₓ[Ξ]²]"):
        assert raw[drift] is False and red[drift] is True
    assert sum(red.values()) > sum(raw.values())


def test_reduced_pam_is_mass_only():
    """Reduced 2D PAM matches the literature `u(ξ − C)`: the u-mass counterterm survives, the
    gradient (drift) counterterms vanish (root Xⁿ)."""
    import json
    from sympy import Rational
    from counterterms import Noise, Parabolic, SPDE, Unknown, kappa
    u = Unknown("u", 2)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    eq = SPDE(noises=[xi], operator=Parabolic(dim=2), unknown=u, rhs=u.field * xi.symbol).renormalize()
    red = {r["tree"]: r["vanishes"] for r in json.loads(eq.render("json", reduced=True))["canonical_bphz"]}
    assert red["Ξ·𝓘[Ξ]"] is False                          # u·C mass survives
    assert red["Xₓ·Ξ"] is True and red["X_y·Ξ"] is True     # gradient drift vanishes


def test_reduced_does_not_claim_reflection_for_asymmetric_noise():
    """The spatial-reflection reduction is gated on `symmetric`: with `symmetric=False` (anisotropic
    noise) it is NOT applied — KPZ's drift counterterms stay non-zero — while the noise-independent
    identities (zeros/duplicates) still fold. The JSON records the assumption honestly."""
    import json
    from tests.conftest import kpz
    eq = kpz().renormalize()
    sym = json.loads(eq.render("json", reduced=True, symmetric=True))
    asy = json.loads(eq.render("json", reduced=True, symmetric=False))
    vsym = {r["tree"]: r["vanishes"] for r in sym["canonical_bphz"]}
    vasy = {r["tree"]: r["vanishes"] for r in asy["canonical_bphz"]}
    drift = "●·𝓘ₓ[●·𝓘ₓ[Ξ]]·𝓘ₓ[Ξ]"
    assert vsym[drift] is True and vasy[drift] is False   # claimed 0 only when symmetric
    assert sum(vasy.values()) < sum(vsym.values())        # fewer vanish without reflection
    assert sym["reduction_assumes_symmetric_noise"] is True
    assert asy["reduction_assumes_symmetric_noise"] is False


def test_general_noise_reduction_is_noise_independent_only():
    """Safety invariant: with `symmetric=False`, the `reduced` view must apply ONLY
    noise-independent identities. It may zero h(σ) exactly when `zero_note` does (root Xⁿ /
    parity / pure-kernel total derivative — valid for any centered Gaussian), and every merge
    must be a genuine o-duplicate (`expectation_key`). The spatial-reflection identity must NOT
    leak into the general-noise reduction (no accidental reduction for an anisotropic noise)."""
    from counterterms.structures import build_renormalization
    from counterterms.renorm.scheme import expectation_key, zero_note
    from counterterms.render.report import _reduction_subs
    from tests.conftest import CORPUS
    for name, mk in CORPUS.items():
        rs = build_renormalization(mk())
        sig = rs.sig
        for t in rs.divergent:
            rs.canonical_character(t)
        keyof = {sym: expectation_key(tr) for tr, sym in rs._h.items()}
        zn = {sym for tr, sym in rs._h.items() if zero_note(tr, sig)}
        asym = _reduction_subs(rs._h, sig, symmetric=False)
        zeros = {s for s, v in asym.items() if v == 0}
        merges = {s: v for s, v in asym.items() if v != 0}
        assert zeros == zn, f"{name}: general-noise reduction zeroed a symmetry-dependent h"
        assert all(keyof[s] == keyof[d] for s, d in merges.items()), \
            f"{name}: general-noise merge is not a genuine o-duplicate"
