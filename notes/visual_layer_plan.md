# Visual Layer Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Trees draw themselves as SVG, and `eq.render("html")` produces a single offline HTML explorer (equation, tree gallery, per-tree detail with coproduct pictures and Wick-pairing diagrams).

**Architecture:** Two new read-only renderers under `counterterms/render/`: `svg.py` (one tree → one `<svg>` string) and `html.py` (one `RenormalizedEquation` → one HTML document, all panels rendered in Python, a few lines of vanilla JS only toggling visibility). One small helper, `renorm.scheme.pairings`, exposes the Wick matchings `expectation` already enumerates.

**Tech Stack:** Python ≥ 3.11, SymPy (presentation MathML printer), hand-written SVG/HTML/CSS/JS strings, pytest, ruff. Run everything with `uv run`.

**Spec:** `notes/visual_layer_design.md` — read it first; this plan argues from it.

## Global Constraints

- No new runtime dependency; no JS framework; SymPy stays the only dependency.
- The HTML file must open offline: no `http://` or `https://` reference except the SVG/MathML namespace URIs.
- No edit to the mathematics in `trees/coproducts.py`, the generator, or the existing functions of `renorm/scheme.py` (only: add `pairings`, make `expectation` call it).
- No change to the text / markdown / JSON / LaTeX reports or to `save`'s default `formats`.
- Never show canonical constants for a lifted equation (`eq.sig.non_gaussian` non-empty).
- Every function in `counterterms/` fully type-annotated; `from __future__ import annotations` at the top of each module; domain types imported under `if TYPE_CHECKING:`. Tests are exempt.
- `uv run ruff check .` clean and `uv run pytest -q` green (195 tests before this work) after every task.
- All text placed into SVG or HTML goes through `html.escape` (noise and unknown names are user strings).
- Work on a branch `visual-layer`; one commit per task. End each commit message with
  `Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>`.

## Review Focus

1. **A noise or unknown named with markup characters** (`Noise("a<b", …)`): the SVG and HTML must still parse. Test in Task 1 (`test_svg_escapes_names`) and Task 2 (`test_html_escapes_names`).
2. **Labels or nodes outside the drawing**: a deep chain (KPZ's 4-node spine), a node with three children, and a decorated root must all lie inside the `viewBox`. Test in Task 1 (`test_svg_everything_inside_viewbox`).
3. **A two-component system**: each tree's card must link to its summand in *every* component it contributes to, and edges must carry the `(c)` label. Test in Task 2 (`test_html_system_links_each_component`).
4. **An equation with trees but few or no counterterms, and cards without a constant** (KPZ's three `F(τ*) = 0` trees): the page must render and those cards must be inert, not broken links. Test in Task 2 (`test_html_kpz_marks_non_contributors`).
5. **A tree whose coproduct has more than 12 terms**: the overflow must sit inside a `<details>` and no term may be dropped. Test in Task 3 (`test_coproduct_overflow_is_kept`).

---

## File structure

| File | Responsibility |
|---|---|
| `counterterms/render/svg.py` (new) | `svg(t, sig, pairing, scale)`: layout + glyphs for one tree |
| `counterterms/render/html.py` (new) | `html_document(eq, …)`: page assembly, MathML printing, CSS/JS constants |
| `counterterms/render/tree.py` | unchanged; `svg.py` imports `_noises`, `_display_sig`, `coord_names`, `_poly_text` |
| `counterterms/render/__init__.py` | export `svg`, `html_document` |
| `counterterms/render/report.py` | `render()` gains the `"html"` branch only |
| `counterterms/renorm/equation.py` | `save()` learns the `html` extension |
| `counterterms/renorm/scheme.py` | new `pairings`; `expectation` calls it |
| `counterterms/trees/tree.py` | `DecoratedTree._repr_html_` |
| `tests/test_svg.py`, `tests/test_html.py` (new) | the checks below |

Test fixtures already exist in `tests/conftest.py`: `gkpz()`, `kpz()`, `system()`, `multinoise()` (each returns an `SPDE`).

---

### Task 1: The SVG tree drawer

**Files:**
- Create: `counterterms/render/svg.py`, `tests/test_svg.py`
- Modify: `counterterms/render/__init__.py`, `counterterms/trees/tree.py` (next to `_repr_latex_`, line 91)

**Interfaces:**
- Produces:
  `svg(t: DecoratedTree, sig: Signature | None = None, *, pairing: tuple[tuple[int, int], ...] = (), scale: float = 1.0) -> str`
  — a standalone `<svg xmlns="http://www.w3.org/2000/svg" viewBox="…" width="…" height="…">…</svg>`.
  `DecoratedTree._repr_html_(self) -> str` returning `svg(self)`.
- Element contract later tasks and tests rely on:
  - every node: one element with `class="n n-bullet"` | `"n n-noise"` | `"n n-red"` | `"n n-blue"` and `data-id="{k}"`; a noise additionally gets `n-noise-{j}` where `j` is its index in `_noises(sig)`;
  - every edge: one `<line>` with `class="e"` (kernel) or `class="e e-d"` (derivative kernel, `any(p)`);
  - every pairing arc: one `<path class="pair">`;
  - labels are `<text class="lbl">`.
- Node ids: pre-order over `t.children` in stored (canonical) order, root = 0 — the numbering of `trees.coproducts._explode`.

- [ ] **Step 1: Write the failing tests** in `tests/test_svg.py`

```python
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
```

- [ ] **Step 2: Run to verify it fails**

Run: `uv run pytest tests/test_svg.py -q` — Expected: collection error, `cannot import name 'svg'`.

- [ ] **Step 3: Implement `svg` in `counterterms/render/svg.py`**

Decisions the tests do not fix:
- Layout: walk in pre-order assigning ids; leaves take consecutive x-slots, an inner node sits at the mean x of its children, `y = (max_depth − depth) · STEP` so the root is lowest. `STEP = 34`, slot width `30`, margin `14` plus room for labels; multiply all by `scale`.
- Glyphs: `●` → `<circle r="3.2" fill="currentColor">`; noise → `<circle r="4.2" fill="none" stroke="currentColor">` for type 0, `fill="#9aa0a6"` for type 1, `fill="#4a7bd0"` for type 2, and beyond three a `Ξ_name` `<text class="lbl">` below the node (mirrors `_NOISE_STYLE` in `render/tree.py:205`); red → `<rect>` side 7, `fill="#c0392b"`, with `x`/`y` attributes; blue → `<circle fill="#2c6fbb">`.
- Edges: `<line stroke="currentColor">`, derivative kernels add `stroke-dasharray="1.5 2.5"`. Mid-edge label when `sig.n_components > 1` (`(c)`) and/or when `any(p) and sig.dim > 1` (direction letters) — the same conditions as `_edge_opts` (`render/tree.py:258`).
- Labels: polynomial decoration (`_poly_text(node_dec, coords)`) upper right of the node; `o=…` above a red node when `o.std != 0 or o.kap != 0`. `font-size="8"`.
- Pairing arc: quadratic Bézier between the two nodes bulging upward by 14, `stroke-dasharray="3 2"`, `fill="none"`, `stroke="#c0392b"`.
- The `viewBox` is computed from the extremes of everything drawn, labels included (estimate a label's width as `0.6 · font_size · len(text)`).

- [ ] **Step 4: Add `_repr_html_` to `DecoratedTree` and export `svg` from `counterterms/render/__init__.py`** (add to the import line and `__all__`).

- [ ] **Step 5: Run the tests**

Run: `uv run pytest tests/test_svg.py -q && uv run pytest -q && uv run ruff check .` — Expected: all pass, `All checks passed!`.

- [ ] **Step 6: Look at it.** Write the five gKPZ trees and the eleven KPZ trees into one scratch HTML file under `output/` and open it; fix anything that overlaps or is unreadable. (`output/` is gitignored.)

- [ ] **Step 7: Commit** — `feat(render): SVG tree drawer; trees display as pictures in Jupyter`

---

### Task 2: The explorer shell

**Files:**
- Create: `counterterms/render/html.py`, `tests/test_html.py`
- Modify: `counterterms/render/report.py` (`render`, lines 70–86), `counterterms/render/__init__.py`, `counterterms/renorm/equation.py` (`save`, the `ext` dict)

**Interfaces:**
- Consumes: `svg` (Task 1); from `render/report.py`: `_sorted_trees(eq)`, `const_map(eq)`, `operator_str(op)`, `_prime`; `eq._pretty(expr)`, `eq.original_rhs(comp)`, `eq.per_component`, `eq.unknowns`, `eq.spde.equations`, `Counterterm.divergence`.
- Produces:
  `html_document(eq: RenormalizedEquation, canonical: bool = False, reduced: bool = False, symmetric: bool = True) -> str`
  `mathml(expr: sympy.Expr) -> str` — presentation MathML wrapped in `<math>…</math>`.
  `eq.render("html", …)`; `eq.save(formats=("html",))` writes `{stem}.html`.
- DOM contract (Tasks 3–4 and the tests rely on it):
  - `<article class="card" data-tree="{i}">` per tree, `i` = index in `_sorted_trees(eq)`; non-contributing trees add class `inert` and contain the text `F(τ*) = 0 — no counterterm`;
  - `<span class="ct" data-tree="{i}" data-comp="{a}">` per (tree, component) summand in the equation block;
  - `<section class="detail" id="tree-{i}">` per tree (all trees, inert ones included), containing `<div class="facts">`; Tasks 3 and 4 append `<div class="coproduct">` and `<div class="canonical">` inside it;
  - gallery groups: `<div class="grade" data-hom="{str(homogeneity)}">`, most singular first.
- A private hook for the next tasks: `_detail_extras(eq, t, fam) -> str` returning `""` in this task (`fam` is the `CanonicalFamily` or `None`).

- [ ] **Step 1: Write the failing tests** in `tests/test_html.py`

```python
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
```

- [ ] **Step 2: Run to verify it fails** — `uv run pytest tests/test_html.py -q` → `unknown render format 'html'`.

- [ ] **Step 3: Implement `mathml`** in `html.py`: a `MathMLPresentationPrinter` subclass whose `_print_Derivative` uses `report._prime` to emit `f′(u)` (as `<mrow><mi>f</mi><mo>′</mo>…`) and otherwise defers to the parent. Look at the output for `eq._pretty(ct.term)` on gKPZ and KPZ in a browser; if it is unacceptable, switch `mathml` to `<code>` + `html.escape(fstr(expr))` and note the reason in the module docstring. Either way every caller goes through `mathml`.

- [ ] **Step 4: Implement `html_document`.** One function assembling: `<!doctype html>` head with inline `_CSS`; header (title, domain/noise line: each noise name and regularity); equation block per component — `operator_str(op) u = ` + `mathml(original_rhs)` + one `<span class="ct">` per counterterm of that component holding `mathml(eq._pretty(ct.term))`; gallery (`div.grade` per homogeneity, `article.card` holding `svg(t, eq.sig)` and `k_τ` or the inert label); one `section.detail` per tree with a larger `svg(t, eq.sig, scale=2)` and `div.facts` (`|τ|`, `S(τ)`, and for contributing trees `k_τ`, the ε-rate from `Counterterm.divergence`, `F(τ*)` per component) followed by `_detail_extras(eq, t, fam)`; inline `_JS`.
  - `fam = eq.canonical(reduced=reduced, symmetric=symmetric) if (canonical or reduced) else None` — so a lifted equation raises exactly when the canonical view was requested.
  - `_JS` (vanilla): on `mouseenter`/`mouseleave` of any `[data-tree]`, toggle class `hl` on every element with the same `data-tree`; on `click`, add class `sel` to that tree's card and show only `#tree-{i}`; `ArrowLeft`/`ArrowRight` move the selection through the non-inert and inert cards alike; on load select the first card. Before the script runs, all `section.detail` are visible (no-script fallback): the script adds a class to `<body>` that enables the hiding rule.

- [ ] **Step 5: Wire it.** `render/report.py::render`: add `if fmt == "html": from .html import html_document; return html_document(eq, canonical, reduced, symmetric)` and extend the error message's format list. `renorm/equation.py::save`: add `"html": "html"` to `ext`. Export `html_document` from `render/__init__.py`.

- [ ] **Step 6: Run** — `uv run pytest tests/test_html.py -q && uv run pytest -q && uv run ruff check .` → all pass.

- [ ] **Step 7: Look at it.** `uv run python -u -c "from tests.conftest import kpz; kpz().renormalize().save(stem='kpz', formats=('html',), pdf=False)"`, open `output/kpz.html`, check hover linking, selection, arrow keys, a narrow window, and dark mode (`prefers-color-scheme`).

- [ ] **Step 8: Commit** — `feat(render): single-file HTML explorer — equation, tree gallery, linked detail panels`

---

### Task 3: The coproduct in the detail panel

**Files:**
- Modify: `counterterms/render/html.py`, `tests/test_html.py`

**Interfaces:**
- Consumes: `trees.coproducts.delta_minus(t: DecoratedTree, sig: Signature) -> TensorSum`, i.e. `dict[tuple[tuple[DecoratedTree, ...], DecoratedTree], Fraction]` (extracted forest, contracted tree) → coefficient.
- Produces: `_coproduct_html(t: DecoratedTree, sig: Signature) -> str`, called from `_detail_extras`; emits `<div class="coproduct">` containing one `<span class="term">` per term.

- [ ] **Step 1: Write the failing tests** (append to `tests/test_html.py`)

```python
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
```

(Checked when this plan was written: `kpz()` has 11 trees and its largest `δ` has 76 terms; gKPZ's largest has 7.)

- [ ] **Step 2: Run to verify it fails** — `uv run pytest tests/test_html.py -q -k coproduct`.

- [ ] **Step 3: Implement `_coproduct_html`.** Heading `δτ =`. Terms sorted by `(len(forest), sum of nodes in forest)`. Each term: the coefficient (omitted when 1, as a fraction otherwise), the forest's trees drawn side by side with `svg(·, sig)` — the empty forest as `𝟙` — then `⊗`, then the contracted tree. The first 12 terms inline; the remainder inside `<details><summary>{n} more terms</summary>…</details>`.

- [ ] **Step 4: Run** — `uv run pytest -q && uv run ruff check .` → all pass.

- [ ] **Step 5: Look at it** — regenerate `output/kpz.html`; check that a red contraction node with its `o=` label is legible at card size.

- [ ] **Step 6: Commit** — `feat(render): draw the extraction coproduct δτ in the explorer`

---

### Task 4: Wick pairings and the canonical section

**Files:**
- Modify: `counterterms/renorm/scheme.py` (`expectation`, lines ~205–262), `counterterms/render/html.py`, `tests/test_scheme.py`, `tests/test_html.py`

**Interfaces:**
- Produces:
  `pairings(tree: DecoratedTree, sig: Signature) -> list[tuple[tuple[int, int], ...]]` in `renorm/scheme.py` — each element one global Wick matching as `(i, j)` node-id pairs (ids as in `_explode`). Order of checks, identical to `expectation`: `sig.require_gaussian("pairings()")`; any noise type with an odd count → `[]`; non-bare tree → `NotImplementedError` (same message as today); otherwise the Cartesian product over noise types of `wick_pairings(ids)`, each combination flattened to one tuple.
  `_canonical_html(eq: RenormalizedEquation, t: DecoratedTree, fam: CanonicalFamily) -> str` → `<div class="canonical">`, called from `_detail_extras` when `fam is not None` and `t in fam.constants`.
- Consumes: `CanonicalFamily.constants`, `.legend`, `.rhs(comp)`; `scheme.is_bare`; `svg(..., pairing=...)`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_scheme.py`:

```python
def test_pairings_agree_with_expectation():
    import pytest
    from counterterms import playground
    from counterterms.renorm.scheme import expectation, pairings
    pg = playground("-3/2 - kappa")
    two = pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"))
    four = pg.tree(pg.I(two, "x"), pg.I(two, "x"))
    for t in (two, four):
        ps = pairings(t, pg.sig)
        assert len(ps) == len(expectation(t, pg.sig).terms)
        assert all(len({i for ij in p for i in ij}) == 2 * len(p) for p in ps)   # a perfect matching
    assert len(pairings(four, pg.sig)) == 3 and pairings(two, pg.sig) == [((1, 2),)]
    assert pairings(pg.tree(pg.I(pg.Xi, "x")), pg.sig) == []            # odd ⇒ []
    odd_nonbare = pg.tree(pg.I(pg.Xi, "x"), X="x")
    assert pairings(odd_nonbare, pg.sig) == []                           # parity before bareness
    with pytest.raises(NotImplementedError):
        pairings(pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"), X="x"), pg.sig)
```

Append to `tests/test_html.py`:

```python
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
```

- [ ] **Step 2: Run to verify they fail** — `uv run pytest tests/test_scheme.py tests/test_html.py -q -k "pairings or canonical"`.

- [ ] **Step 3: Implement `pairings` and make `expectation` call it.** `expectation` keeps its signature, docstring and output exactly: it builds `z`, `intvars` and `kernels` as now and loops `for matching in pairings(tree, sig)`. The existing `tests/test_scheme.py` and `tests/test_render.py` must pass unchanged — they are the regression check that the refactor is behaviour-preserving.

- [ ] **Step 4: Implement `_canonical_html`.** Shows `k_τ = ` `mathml(value)`; if the value is 0, the words `vanishes canonically`. For each `h_i` in `value.free_symbols` (ordered by index): `h_i = h_ε(` `svg(σ_i)` `)`, then for a bare `σ_i` one `svg(σ_i, sig, pairing=m)` per `m in pairings(σ_i, sig)`; for a non-bare `σ_i` the line `left symbolic — needs the full canonical model`. Also, in `html_document`, when `fam is not None` add under the equation block one unlinked line per component: `canonically renormalized: ` + `mathml(fam.rhs(a))`.

- [ ] **Step 5: Run** — `uv run pytest -q && uv run ruff check .` → all pass.

- [ ] **Step 6: Look at it** — `kpz().renormalize().save(stem='kpz_reduced', formats=('html',), pdf=False, reduced=True)`; check arcs do not hide nodes on the 4- and 6-noise trees.

- [ ] **Step 7: Commit** — `feat(render): Wick-pairing diagrams and the canonical section in the explorer; scheme.pairings`

---

### Task 5: Documentation and the demo

**Files:**
- Create: `docs/guide/seeing-trees.md`, `docs/kpz_explorer.html`
- Modify: `mkdocs.yml` (nav), `README.md`, `docs/api.md`, `CHANGELOG.md`, `ROADMAP.md`, `notes/output.md`, `CLAUDE.md`, `ENTRYPOINTS.md`, `tests/test_render.py` (`__main__` demo), `notes/visual_layer_design.md` (status line)

- [ ] **Step 1:** Generate `docs/kpz_explorer.html` with `kpz().renormalize().render("html", reduced=True)` and link it from the README's "Example output" section next to the two KPZ PDFs.
- [ ] **Step 2:** Write `docs/guide/seeing-trees.md`: the glyph table from the spec, `svg(t, sig)`, the notebook behaviour (`_repr_html_`), `eq.render("html")` / `eq.save(formats=("html",))`, what the explorer shows, and its limits (all trees listed; no canonical section for lifted equations; coproduct terms are `forest ⊗ tree` pairs with no provenance on the original tree). Add it to `mkdocs.yml` under "User guide" after "Reports & output". Every snippet in the page must be run before it is written in.
- [ ] **Step 3:** `docs/api.md`: add `svg`, `html_document`, `"html"` in the `render` format list, `pairings`. README capability tables: one row each for the SVG drawer and the explorer. `ENTRYPOINTS.md`: one line each for `render/svg.py` and `render/html.py`. `CLAUDE.md`: add `svg`, `html` to the `render/{…}` list and the new test count.
- [ ] **Step 4:** `CHANGELOG.md` `[Unreleased]` → `### Added` with the three items (SVG drawer + Jupyter display, HTML explorer, `scheme.pairings`). `ROADMAP.md` "Output & rendering": add `- [x] **O5:**` describing both renderers. `notes/output.md`: add the two renderers to its renderer list. `notes/visual_layer_design.md`: change the status line to "built" with the date.
- [ ] **Step 5:** Update every stated test count (`README.md`, `docs/validation.md`, `docs/getting-started.md`, `CLAUDE.md`, `ROADMAP.md`) to the number `uv run pytest -q` now reports.
- [ ] **Step 6:** `tests/test_render.py` `__main__`: also write `output/gkpz.html`.
- [ ] **Step 7: Verify** — `uv run pytest -q && uv run ruff check . && uv run --group docs mkdocs build --strict -d /tmp/counterterms_site` and run all six `examples/*.py`. Expected: all green, docs build with no warnings.
- [ ] **Step 8: Commit** — `docs: the visual layer — SVG trees and the HTML explorer`
