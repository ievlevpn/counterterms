# Visual layer — design

Status: **built 2026-11-10** (plan: `notes/visual_layer_plan.md`). Written 2026-11-10. Decisions marked *(call)* were made by the
designer, not requested by the user; they are open to veto.

## Intent

The package computes objects people think about as pictures (decorated trees, cuts, pairings)
and shows them as text. A tree prints as `●·𝓘ₓ[●·𝓘ₓ[Ξ]²]·𝓘ₓ[Ξ]`; a real drawing needs `pdflatex`;
the coproducts are deliberately *not* rendered at all because their text form explodes
(`ROADMAP.md`, O3). The visual layer makes the same objects visible without changing any
mathematics.

Two deliverables, in order:

1. **Trees draw themselves as SVG** — in Jupyter, and as a function any renderer can call.
2. **A single-file HTML explorer** for a renormalized equation — the equation, a gallery of its
   divergent trees, and a detail view per tree including its extraction coproduct drawn as
   pictures and its Wick pairings drawn as diagrams.

## Goals

- A tree is recognisable at a glance, in the paper's convention (tex 329–337), matching the
  existing LaTeX `forest` drawer glyph for glyph.
- `eq.save(formats=("html",))` produces **one file** that opens offline in a browser.
- In the explorer, the link between a counterterm in the equation and its tree is visible
  (hover one, the other lights up).
- The extraction coproduct `δ` of a tree is readable as a list of pictures.

## Non-goals (this round)

- A `β₀` slider / homogeneity ladder that regenerates trees (needs Python behind the page).
- A Pyodide playground on the docs site; a click-to-build tree editor.
- *(Built later, 2026-11-10, via an opt-in trace — `coproducts.cuts`; a term now shows all its cuts.)*
  Highlighting the extracted subforest *on the original tree*. `delta_minus` returns
  `{(forest, contracted tree): coeff}` and sums coefficients over different node subsets that
  give the same pair, so a term has no unique provenance. Each term is drawn as
  `forest ⊗ contracted tree` instead.
- Rendering `Δ`/`Δ⁺` (recentering) or the full twisted antipode `S'₋` expansion.
- Changing the text / markdown / JSON / LaTeX reports.

## Anti-goals (a failure even if it works)

- A new **runtime** dependency, or any JS framework. SymPy stays the only dependency.
- A file that needs the network to open (no CDN MathJax, no web fonts).
- Any edit to the mathematics in `trees/coproducts.py`, `renorm/scheme.py`'s existing
  functions, or the generator. This layer only *reads*.
- A picture that asserts something the engine refuses: the explorer must not show canonical
  constants for a da Prato–Debussche-lifted equation (`Signature.require_gaussian`).

## Constraints (from the project)

- Full type annotations, `from __future__ import annotations`, domain types (`CLAUDE.md`).
- `ruff check .` clean; every non-trivial piece leaves a runnable check.
- Lazy/minimal: hand-written SVG and HTML strings, fewest files.

---

## Part 1 — `render/svg.py`: the SVG tree drawer

**API** *(call)*

```python
def svg(t: DecoratedTree, sig: Signature | None = None, *,
        pairing: tuple[tuple[int, int], ...] = (), scale: float = 1.0) -> str
```

Returns a standalone `<svg>…</svg>` string. `sig=None` falls back to `_display_sig(t)` exactly
as `shorthand` does. Exported from `counterterms.render`. `DecoratedTree` gains
**`_repr_html_`** returning the inline SVG *(call)*; the existing `_repr_latex_` stays as the
fallback. Not `_repr_svg_`: JupyterLab ranks `text/html` (50) ahead of `text/latex` (70) ahead
of `image/svg+xml` (80), so an SVG repr alone would lose to the LaTeX shorthand and nothing
would change. Inline SVG in HTML also lets `currentColor` follow the notebook theme.

**Glyphs** — identical semantics to `forest` (`render/tree.py:233–270`), reuse its helpers
(`_noises`, `coord_names`, the label logic of `_node_opts` / `_edge_opts`):

| Object | Drawing |
|---|---|
| `●` integration node | small filled dot |
| noise `∘_j` | open circle; second/third noise types grey / blue fill, as `_NOISE_STYLE`; a text label `Ξ_name` beyond three |
| red contraction node | red filled square, `o=…` label above when `o ≠ 0` |
| blue `T⁺` root | blue filled dot |
| kernel edge `𝓘` | solid line |
| derivative kernel `𝓘_p`, `p ≠ 0` | dotted line; direction letters as a mid-edge label when `dim > 1` |
| component (systems) | `(c)` mid-edge label when `n_components > 1` |
| polynomial decoration `X^n` | small label to the upper right of the node |

**Layout** *(call)*: root at the **bottom**, children above (the paper's upward trees). Leaves
get consecutive x-slots in canonical child order; an inner node sits at the mean x of its
children; y is depth × a fixed step. No collision-avoidance beyond that — the trees are small.
Labels are SVG `<text>`; unicode for `κ`, sub/superscripts as `<tspan>` shifts.

**Styling** *(call)*: strokes and fills use `currentColor` for black elements so the picture
follows the surrounding text colour (light and dark notebooks); red/blue/grey are fixed. Every
node element carries `class="n n-{bullet|noise|red|blue}"` and `data-id="{k}"`, edges
`class="e"` / `class="e e-d"`; the explorer styles and targets these.

**Node ids.** `data-id` is the pre-order index over canonical children — the same numbering
`trees.coproducts._explode` produces (it is what `scheme.expectation` uses for `z1, z2, …`).
`pairing=((i, j), …)` draws a dashed arc between noise nodes `i` and `j`.

**Check.** `tests/test_svg.py`: output parses as XML; node and edge element counts equal
`t.nodes()` and `t.nodes() − 1`; the five gKPZ goldens have the expected glyph classes
(e.g. `●·𝓘ₓ[Ξ]²` → one `n-bullet`, two `n-noise`, two `e-d`); a red-node tree carries its `o`
label; a two-noise tree distinguishes the two noise styles.

## Part 2 — `render/html.py`: the explorer

**API** *(call)*

```python
def html_document(eq: RenormalizedEquation, canonical: bool = False,
                  reduced: bool = False, symmetric: bool = True) -> str
```

wired as `eq.render("html", …)` and into `eq.save` (`ext["html"] = "html"`). `"html"` is **not**
added to `save`'s default `formats` *(call)* — existing callers keep their four files.

**Page layout** (one column, max-width, system fonts):

```
┌──────────────────────────────────────────────────────────────┐
│ (∂t − Δ)u = f(u)ξ + g(u)(∂ₓu)²  + [k₀·…] + [k₁/2·…] + …      │  ← equation; each
│ domain · noise regularities · N trees                        │    counterterm a span
├──────────────────────────────────────────────────────────────┤
│ |τ| = −1 − 2κ      |τ| = −1/2 − 3κ        |τ| = −κ …         │  ← gallery, grouped
│ ┌─────┐ ┌─────┐    ┌─────┐                ┌─────┐ ┌─────┐    │    by homogeneity,
│ │ svg │ │ svg │    │ svg │                │ svg │ │ svg │    │    most singular first
│ │ k₀  │ │ k₁  │    │ k₂  │                │ k₃  │ │ k₄  │    │
│ └─────┘ └─────┘    └─────┘                └─────┘ └─────┘    │
├──────────────────────────────────────────────────────────────┤
│ detail (the selected tree)                                   │
│  large svg │ |τ|, S(τ), k_τ, ε-rate, F(τ*)                   │
│            │ canonical value (if canonical/reduced)          │
│  δτ =  1·(∅ ⊗ [svg])  +  1·([svg] ⊗ [svg])  + …              │
│  pairings of h(σ): [svg+arcs] [svg+arcs]   = ∫ …             │
└──────────────────────────────────────────────────────────────┘
```

- **Which trees.** The gallery shows **all of `eq.all_trees`** *(call)*, as the existing
  reports do (KPZ: 11 trees, 8 counterterms). A tree with `F(τ*) = 0` in every component gets
  a dimmed card labelled "F(τ*) = 0 — no counterterm": it has no `k_τ`, no linked summand and
  no canonical value, but still has a detail panel with `|τ|`, `S(τ)` and its coproduct.
- **Equation block.** Always the **free-constant** family, one summand per (tree, component)
  *(call)* — this is what keeps the one-span-per-tree linking. With `canonical`/`reduced` a
  second, unlinked line below it shows the canonically renormalized right-hand side
  (`CanonicalFamily.rhs`, which expands and merges terms); per-tree canonical values live in
  the detail panels.
- **Equation ↔ tree linking.** Each counterterm summand is a `<span data-tree="k">`; each
  gallery card has the same `data-tree`. Hovering either adds a highlight class to both;
  clicking a card (or a summand) selects it and fills the detail panel. For systems the
  equation block repeats per component, and a tree's card links to its summand in each.
- **Detail panel.** All panels are rendered server-side (by Python) into hidden `<section>`s;
  selection only toggles visibility. No client-side computation.
- **Coproduct.** The terms of `delta_minus(τ, sig)` (the coaction `δ : U → U⁻ ⊗ U`), each as
  `coeff · (forest drawn side by side) ⊗ (contracted tree)`. The empty forest is drawn as `𝟙`.
  Terms sorted by forest size; the first 12 shown, the rest behind a `<details>` *(call)*.
- **Canonical section** (only with `canonical`/`reduced`): the tree's `k_τ` value from
  `eq.canonical(...)`, and for each `h_i` in it the tree `σ_i` drawn, with one small diagram
  per Wick pairing (bare `σ_i` only; non-bare ones show the tree and a "left symbolic — needs
  the full canonical model" note). That wording exists only inline in the text and markdown
  reporters (`report.py:374`, `:499`); `html.py` writes its own line rather than refactoring
  them.
- **Lifted equations.** If `eq.sig.non_gaussian` is non-empty, `html_document` renders
  everything except the canonical section and prints one line saying why, instead of
  propagating the `NotImplementedError` — *but only when `canonical`/`reduced` was not
  requested*. If it was requested, the error propagates, as in every other format.

**Formulas** *(call)*: SymPy presentation MathML (`sympy.mathml(expr, printer="presentation")`)
of the already-prettified expressions (`eq._pretty`, as the text report uses). Native MathML
renders in current Chrome, Firefox and Safari with no script. Expressions must go through
`_pretty` first (raw jets print as `u₀ ₀ ₀`), and the MathML printer needs the same derivative
handling `fstr` has (SymPy prints `f'` as a `d/du` fraction by default) — a small printer
subclass. If the result is still ugly, fall back to the report's `fstr` strings in `<code>`;
the builder decides after looking at gKPZ and KPZ.

**Script and style**: one inline `<style>` and one inline `<script>` (vanilla, ≈60 lines:
hover-highlight, select, keyboard left/right through the gallery). The page is fully readable
with scripts disabled — all detail sections simply show, stacked.

**One new read-only helper** in `renorm/scheme.py`:

```python
def pairings(tree: DecoratedTree, sig: Signature) -> list[tuple[tuple[int, int], ...]]
```

the global Wick matchings that `expectation` already enumerates internally (product over
noise types of `wick_pairings`), as node-id pairs. `expectation` is refactored to call it, so
the diagrams and the integrals cannot disagree. It keeps `expectation`'s exact order of
checks (`scheme.py:234–243`): `sig.require_gaussian` first, then **odd parity returns `[]`**,
then a non-bare tree **raises** `NotImplementedError`. The explorer calls it only for bare
`σ_i`.

**Check.** `tests/test_html.py`: the gKPZ document has five cards and five linked summands;
KPZ has 11 cards of which 8 are linked and 3 are marked "no counterterm"; contains no `http://`/`https://` reference other than the SVG/MathML namespace
URIs; `canonical=True` adds the canonical section and pairing diagrams; a lifted Φ⁴₃ renders
without it and raises with `canonical=True`; `pairings` agrees term-for-term in count with
`expectation(...).terms`. The demo in `tests/test_render.py`'s `__main__` also writes
`output/gkpz.html`.

## Build order

Each step is shippable on its own.

1. `render/svg.py` + `_repr_html_` + tests.
2. `render/html.py`: equation, gallery, detail with the scalar data, linking; `render`/`save`
   wiring.
3. Coproduct terms in the detail panel.
4. `scheme.pairings` + pairing diagrams in the canonical section.
5. Docs: a "Seeing the trees" page in `docs/guide/`, one explorer file linked from the README
   next to the two KPZ PDFs, CHANGELOG, `notes/output.md` (new renderers), `ROADMAP.md` O-list.

## Left to the builder

Exact pixel sizes, colours and typography; the CSS; how the detail panel is laid out on a
narrow window; label placement tweaks; whether `Counterterm` also gets an SVG `_repr_html_`;
the wording of notes; how many trees per gallery row.
