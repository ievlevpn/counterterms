# `regstruct` — output & rendering design

A plan for **what we can emit and how to make it pretty**. Companion to
`architecture.md` §3.11/§5 (which already reserves a `render/` package and a
`Renderer` seam). This document is the spec for that package: the *information
inventory*, the *tree notation*, the *output channels*, and the *full-report
layout*. Nothing here is implemented yet — Phase-3/4 items are explicit
placeholders so the report degrades gracefully to "not yet computed".

Guiding principle (ponytail): **reuse the printers we already have.** SymPy
ships a 2D-unicode pretty-printer (`sympy.pretty`) and a LaTeX printer
(`sympy.latex`); we own only (a) the *tree* drawing, which SymPy can't do, and
(b) the *report assembly*. No new runtime dependency — box-drawing is unicode
literals, LaTeX trees use the `forest` package (a `.tex` consumer concern, not a
Python one). `rich`/ANSI colour is a later rung, added only if asked.

---

## 0. What "all useful information" means here

Three audiences, one pipeline:

1. **The researcher at the terminal** — wants to *read* the renormalized family
   and *see* the divergent trees, fast, in a scrollback-friendly form.
2. **The paper / notebook** — wants compilable LaTeX (trees that look like
   Bailleul–Hoshino's) and Markdown.
3. **A downstream tool** — wants structured JSON (the clean handoff at the
   symbolic boundary, `architecture.md` §3.11).

Everything below serves all three from one information inventory.

---

## 1. Information inventory — what we can emit

Status legend: **now** = computable from Phase-1 objects today · **easy** =
trivially derivable, just not surfaced · **P3/P4** = needs a later phase
(rendered as a labelled placeholder until then).

| # | Item | Source object | Status |
|---|------|---------------|--------|
| 1 | The SPDE as parsed (echo: `L u = f(u)ζ + g`) | `spde.rhs`, `spde.unknown`, `spde.operator` | now |
| 2 | Operator `L`, scaling 𝔰, spatial dim `d`, mass | `Signature.scaling`, `Parabolic` | now |
| 3 | Noise(s) and regularity `β_j = |Ξ_j|` | `Noise.homogeneity` | now |
| 4 | Subcriticality verdict (+ why it terminated) | generation terminating + bound | easy |
| 5 | Derived structural **rule** `R` (per node: allowed child edges + caps) | `Signature.allowed` | easy |
| 6 | `𝓑_{<0}` — **every** negative-homogeneity tree | `generate_counterterms` | now |
| 7 | per tree: drawing + shorthand | `DecoratedTree` (§3) | easy |
| 8 | per tree: `|τ|` (homogeneity, exact in `ℚ+ℚκ`) | `tree.homogeneity(sig)` | now |
| 9 | per tree: `S(τ)` symmetry factor | `tree.symmetry_factor()` | now |
| 10 | per tree: node count, depth | `tree.nodes()` (+ a `depth()` one-liner) | easy |
| 11 | per tree: `F(τ*)` (the Υ-map vector field, in jet vars) | `Counterterm.elem_diff` | now |
| 12 | per tree: free constant `k_τ` and coefficient `k_τ/S(τ)` | `Counterterm.constant/.coefficient` | now |
| 13 | per tree: contributes a counterterm? (or `F(τ*)=0`, D2-dropped) | api filter (§6 note) | easy |
| 14 | The assembled **renormalized family** RHS | `counterterm_rhs()` | now |
| 15 | Counts: #trees, #counterterms, #dropped | derived | easy |
| 16 | Coproducts `Δτ`, `Δ⁻τ` per tree | `trees/coproducts` | **P3** |
| 17 | Twisted antipode `S'₋ τ`, BHZ character (symbolic) | `core/hopf`, `structures/` | **P3** |
| 18 | The regularity structure `(T, T⁺)` homogeneity spectrum | `structures/regularity` | **P3** |
| 19 | Canonical constant *values* `c_τ` under a `NoiseLaw` | `BPHZ` scheme | **P4** |

Items 1–15 are a *complete, useful report today*. 16–19 are reserved slots.

---

## 2. Output channels

One dispatcher, four backends. Default is `text`.

| fmt | use | engine |
|-----|-----|--------|
| `text` | terminal report (default) | unicode box-drawing + `sympy.pretty` |
| `latex` | paper / standalone compilable `.tex` | `sympy.latex` + `forest` trees |
| `markdown` | README / Jupyter / GitHub | `$…$` math + fenced ASCII trees |
| `json` | machine handoff | `dataclasses.asdict` + a small tree encoder |

`text` may emit ANSI colour **only** when `sys.stdout.isatty()` (dim the
homogeneity column, bold the tree glyphs). Colour is cosmetic and stripped in
pipes — no `rich` dependency unless the user asks for tables/panels.

---

## 3. Tree notation — the crux

A `DecoratedTree` is `b^n ⋆ ⨉ᵢ I_{p_i}(τ_i)`: a root of `node_type` (noise `∘_j`
/ bullet `●` / red), a polynomial node decoration `n` (a multi-index over the
unknown's coords), and a multiset of edges `(operator, p, subtree)`. We render it
**three ways**, all from the same datatype.

### Glyph table

| object | glyph (text) | LaTeX | notes |
|--------|-------------|-------|-------|
| noise `∘_j` | `Ξ` (subscript name if >1 noise) | `\Xi` / `\Xi_{j}` | base map `F=f_j` |
| bullet `●` | `●` | `\bullet` | base map `F=g` |
| red | `▲` | `\textcolor{red}{\bullet}` | `F=0`; never in a counterterm |
| node decoration `n` | `X^{n}` factor, named per coord | `X^{n}` | e.g. `d=1`: `(a,b) → X_t^a X_x^b` |
| edge `I_p` | `𝓘` (subscript = `p` named) | `\mathcal{I}` | the integration/planting kernel |

Coordinate naming comes from `spde.unknown.coords` (`t, x, y, …`), so `p=(0,1)`
prints as the *x-derivative* edge `𝓘ₓ`, not `𝓘_{(0,1)}`. This is what makes it
"nice" instead of index soup.

### 3a. Linear shorthand (one tree → one string)

For inline use, tables, JSON keys, log lines. Grammar:

```
tree   := factors                         # product of node-dec · edge subtrees
node   := Ξ | Ξ_j | ●                      # node_type glyph
dec    := X_c^k …                          # node decoration, omitted if 0
edge   := 𝓘[ tree ] | 𝓘_c[ tree ]          # 𝓘_c = derivative edge on coord c
power  := edge^m                           # repeated identical edge → exponent
```

The five gKPZ counterterms (golden test) in shorthand:

| tree | shorthand | `|τ|` | `S` |
|------|-----------|------|-----|
| bare noise | `Ξ` | `−1−κ` | 1 |
| noise + poly | `X_x·Ξ` | `−κ` | 1 |
| bullet→∂ₓ noise | `●·𝓘ₓ[Ξ]` | `−κ` | 1 |
| noise→noise | `Ξ·𝓘[Ξ]` | `−2κ` | 1 |
| bullet→two ∂ₓ noises | `●·𝓘ₓ[Ξ]²` | `−2κ` | 2 |

### 3b. 2D terminal drawing (root at top, `tree(1)`-style)

The readable form for scrollback. Edge label rides the connector; node
decoration trails the glyph.

```
●                         ●                    Ξ
└─𝓘ₓ─ Ξ                   ├─𝓘ₓ─ Ξ              └─𝓘─ Ξ
                          └─𝓘ₓ─ Ξ
  ●·𝓘ₓ[Ξ]                  ●·𝓘ₓ[Ξ]²            Ξ·𝓘[Ξ]
```

Pure recursion over `children`; `├─/└─/│` chosen by last-child. ~25 lines. The
existing `RenormalizedEquation.summary()` (renorm/equation.py:60) is the seed —
it already sorts by homogeneity and substitutes jets; we swap its one-line term
for `shorthand + drawing + pretty(F(τ*))`.

### 3c. LaTeX (the paper form) — `forest`

Bailleul–Hoshino draw rooted trees growing upward with glyph nodes. Emit a
`forest` snippet per tree plus a shared style preamble:

```latex
% preamble (emitted once)
\usepackage{forest}
\forestset{rstree/.style={for tree={grow'=north, parent anchor=north,
  child anchor=south, s sep=6pt, l sep=10pt, inner sep=1pt}}}

% one tree: ●·𝓘ₓ[Ξ]²
\begin{forest} rstree
  [{$\bullet$}
    [{$\Xi$}, edge label={\scriptsize $\mathcal I_x$}]
    [{$\Xi$}, edge label={\scriptsize $\mathcal I_x$}]
  ]
\end{forest}
```

`forest` is a documented, ubiquitous package; no Python side cost. (A `tikz`
fallback is possible but `forest` is the lazy correct choice — it does the layout.)

---

## 4. Rendering the renormalized family

The headline object. Reuse SymPy printers for every formula; we only assemble.

**text** (`sympy.pretty` for the RHS, jets shown as `u, ∂ₓu, …`):

```
(∂_t − Δ + 1) u =  f(u)·ξ  +  g(u, ∂ₓu)
                +  c₀·f(u)                        [τ=Ξ,        |τ|=−1−κ]
                +  c₁·f'(u)·∂ₓu                    [τ=X_x·Ξ,    |τ|=−κ ]
                +  2·c₂·f(u)·g(u)·∂ₓu             [τ=●·𝓘ₓ[Ξ],  |τ|=−κ ]
                +  c₃·f(u)·f'(u)                   [τ=Ξ·𝓘[Ξ],  |τ|=−2κ]
                +  c₄·f(u)²·g(u)                   [τ=●·𝓘ₓ[Ξ]², |τ|=−2κ]
```

where `c_τ = k_τ/S(τ)` (so `c₂` already carries the ÷1·2 and `c₄` the ÷2·2 of
the symmetry factor — the load-bearing combinatorics, shown not hidden).

**latex**: `sympy.latex` on `counterterm_rhs()`, with each counterterm's `τ`
typeset as the inline `forest` glyph (§3c) in an annotation column. Wrap as an
`align` environment.

---

## 5. The divergent-tree table

The paper's homogeneity table (tex 6028–6063) is the canonical artifact. Columns:

```
 τ (drawing)      shorthand     |τ|      S(τ)   k_τ      F(τ*)
 ────────────────────────────────────────────────────────────────────
 Ξ                Ξ             −1−κ     1      k₀       f(u)
 X_x─Ξ            X_x·Ξ         −κ       1      k₁       f'(u)·∂ₓu
 ●─𝓘ₓ─Ξ           ●·𝓘ₓ[Ξ]       −κ       1      k₂       g(u)·∂ₓu·f(u)
 …
```

Sorted by `Homogeneity._key()` (exact `ℚ+ℚκ` order — never floats). In `text`
the `τ` column is the multi-line drawing (§3b) and rows are separated by blank
lines; in `latex`/`markdown` it collapses to the shorthand or the `forest` glyph.

**Show the dropped trees too.** `𝓑_{<0}` contains trees with `F(τ*)=0`
(Assumption-D2 reddened / red-node trees) that `api.renormalize` filters out
(api.py:19). For completeness the table lists them in a dim "no counterterm
(F(τ*)=0)" sub-section — "all divergent trees" should mean *all*, with the reason
each non-contributor drops. (Needs: a flag to keep zero-F trees for display;
one line in the pipeline.)

---

## 6. The full report — one command

`render(eq, fmt="text", sections=ALL)` assembles, in order:

```
┌─ HEADER ────────────────────────────────────────────────
│  Equation     (∂_t − Δ + 1) u = f(u) ξ + g(u, ∂ₓu)
│  Domain       d = 1   scaling 𝔰 = (2,1)   spacetime 1+1
│  Noise        ξ ∈ C^{−1−κ}   (β₀ = −1−κ)
│  Subcritical  yes — generation terminated at |τ|<1 bound
├─ RULE (derived) ────────────────────────────────────────
│  ●  →  𝓘[·] (any) , 𝓘ₓ[·] (≤2)        # caps = degree in each ∂u slot
│  Ξ  →  𝓘[·] (any)
├─ DIVERGENT TREES  (𝓑_{<0}: 5 trees, 5 counterterms, 0 dropped)
│  …the table of §5, each row with its §3b drawing…
├─ RENORMALIZED FAMILY ───────────────────────────────────
│  …the §4 block…
├─ ALGEBRA (Phase 3 — not yet computed) ──────────────────
│  Δτ, Δ⁻τ, S'₋τ, BHZ character          [placeholder]
└─ CANONICAL VALUES (Phase 4 — needs NoiseLaw) ───────────
   c_τ numeric                            [placeholder]
```

Placeholder sections print a single greyed line naming what would go there and
which phase delivers it — so the report is *honest about its frontier* and the
layout is stable as phases land.

`latex_document(eq)` wraps the same sections in a standalone, `pdflatex`-ready
preamble (forest + amsmath). `markdown(eq)` for notebooks/README.

---

## 7. Placeholders for later phases (reserved slots, no stubs)

- **§ALGEBRA (P3):** per tree, pretty-print `Δτ ∈ T⊗T⁺` and `Δ⁻τ` as sums of
  `τ' ⊗ τ''` using the §3 tree renderer on both tensor legs; the twisted antipode
  `S'₋τ` likewise. The regularity structure shows its homogeneity spectrum
  (sorted multiset of `|τ|` over the basis of `T`). These reuse §3 entirely — a
  coproduct renders as "tree ⊗ tree" pairs.
- **§CANONICAL VALUES (P4):** once a `NoiseLaw` exists, replace each free `k_τ`
  column with its BPHZ value `k^ζ = h^ζ ∘ S'₋`. Same table, the `k_τ` column
  gains a "value" twin.

Until then each renders the one-line placeholder of §6.

---

## 8. Module / API plan

Smallest surface that covers §1–§6 (`architecture.md` §5 already lists the dir):

```
regstruct/render/
  tree.py     # shorthand(tree,sig)  · ascii(tree,sig)  · forest(tree,sig)
  report.py   # render(eq, fmt, sections) dispatcher + the §6 assembler
              #   text/markdown inline; latex via latex.py
  latex.py    # latex_document(eq); forest preamble; align block of §4
```

Public API (thin wrappers on `RenormalizedEquation`, mirroring `summary()`):

```python
eq.render(fmt="text", sections="all")   # → str   (text/markdown/latex/json)
eq.report()                             # → str   alias for render("text")
eq.latex_document()                     # → str   compilable standalone
render_tree(tree, sig, fmt="text")      # → str   single tree, any format
```

`summary()` (renorm/equation.py:60) becomes a thin call into `render(...,
fmt="text")` so there's one code path, not two.

JSON shape (handoff): `{equation, domain, noises, rule, trees:[{shorthand,
homogeneity:{std,kap}, S, nodes, F_latex, k, contributes}], family_latex}`.

---

## 9. Reuse & dependencies (ladder)

- **Formulas** → `sympy.pretty` (2D unicode) and `sympy.latex`. Already a hard
  dep. We write *zero* math-printing code.
- **Trees** → our own ~60 lines (shorthand + ascii + forest). SymPy can't draw
  rooted decorated trees; this is the only genuinely novel renderer.
- **LaTeX trees** → `forest` (consumer's TeX distro, not a Python dep).
- **Tables** → plain f-string columns. *Skipped:* `rich`/`tabulate` — add only
  if the user wants boxed/coloured tables interactively. `# ponytail: f-string
  columns; rich if interactive tables are wanted.`
- **Colour** → bare ANSI guarded by `isatty()`, off in pipes. No dep.

---

## 10. Phasing of the output work

- **O1 (with Phase 1, now):** `render/tree.py` (all three tree forms) +
  `report.py` text/markdown covering inventory items 1–15; fold `summary()` in.
  Golden check: the gKPZ report reproduces the five trees of §3a and the
  homogeneity order of §5.
- **O2:** `latex.py` standalone document + JSON export.
- **O3 (with Phase 3):** wire the §ALGEBRA section (coproducts/antipode) into the
  existing renderer — pure reuse of §3.
- **O4 (with Phase 4):** the canonical-values column.

One runnable check per renderer (assert-based `demo()` on gKPZ), per project
convention — the report's tree set and homogeneity ordering are the backbone
assertions.
```
