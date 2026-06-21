# Roadmap

Where `regstruct` is going. Phases mirror [`notes/architecture.md`](notes/architecture.md) §7;
the mathematics and conventions are fixed in [`notes/initial_plan.md`](notes/initial_plan.md).
Legend: ✅ done · 🔨 in progress · ⬜ planned.

---

## Phase 1 — equation pipeline ✅

**Goal:** `SPDE → family of renormalized equations` with free symbolic constants. No Hopf
coproducts (the family-with-free-constants needs only rule → trees → `S(τ)` → `F(τ*)`).

- [x] `core/` — ordered-ring homogeneities (`ℚ⊕ℚκ`), jet variables, the `Signature`
- [x] `trees/` — decorated trees, canonical form, `S(τ)`, homogeneity
- [x] `equation/` — SymPy DSL, rule derivation from the nonlinearity, tree enumeration
- [x] `renorm/` — the Υ map `F(τ*)`, the `RenormalizedEquation` output
- [x] golden test (gKPZ, exact 5 counterterms) + KPZ/gPAM + scope rejections

**Acceptance (met):** the gKPZ family reproduces tourist_guide.tex 6004–6012 exactly,
including the factor-2 and `S=2` cases; KPZ lands on the table's 6 negative-homogeneity rows.

---

## Phase 2 — generality ✅

**Goal:** widen the input class with (mostly) `Signature` enrichment + validation; little new
algorithmic structure, by design.

- [x] **Systems / vector unknowns** — multiple components sharing spacetime coords; each child
      **edge carries the component** `𝔗_e` of the kernel that plants it (node types stay
      `{●, ∘_j}`). The elementary differential threads the child equation index from the edges
      (`F_a(τ*) = (Πᵢ F_{cᵢ}(τ_i*))·(D^n Πᵢ ∂_{(cᵢ,p_i)}) F_a(b*)`), and **constants are shared
      across components** (one `k_τ` per tree, weighting each equation's vector field).
- [x] **Multiple noises** — one `∘_j` node type per noise; base map `{●:g, ∘_j:f_{a,j}}`.
- [x] **General operator order** — per-component Schauder order `m` carried into the homogeneity
      recursion (`|I^{(c)}_p τ| = |τ| + (m_c − |p|_𝔰)`), with a warning outside the proven
      second-order parabolic regime. *(Global scaling; genuinely mixed scalings stay out of scope.)*
- [x] Tests: decoupled system reduces to per-component gKPZ; coupled system shares constants;
      multi-noise (8 counterterms); operator order changes the homogeneities.

**Acceptance (met):** a 2-component system and a 2-noise equation produce the correct coupled
counterterms; the decoupled system reproduces the scalar gKPZ result per component; scalar
single-noise results are unchanged (regression green).

> Moved to Phase 3: the **`G⁻` vs `G⁻_ad` toggle**. Admissibility is genuinely a model
> (K-admissibility) notion — the paper's own gKPZ example is stated for `k∈G⁻_ad` yet keeps the
> decorated `∘1` primitive, so a naive "drop `X^n⋆σ`" filter would be wrong. We emit the full free
> family (the safe superset) until the structures exist to characterise `G⁻_ad` precisely.

---

## Phase 3 — the algebraic-renormalization core ⬜

**Goal:** build the genuine regularity-structures machinery. This is the famous-hard part.
New modules: `core/hopf.py`, `trees/coproducts.py`, `structures/`.

- [ ] **Extended-decoration trees** — red nodes with `𝔬 : N^red → ℤ[β₀]`; the contraction
      `τ /^red φ`.
- [ ] **Coproducts** — `Δ`, `Δ⁺` (recentering) and `Δ⁻`, `δ` (extraction/contraction), built as
      `Coproduct` objects over trees; the cointeraction.
- [ ] **Generic Hopf machinery** (`core/hopf.py`, written against the `Symbol` protocol) —
      character convolution, the connected-graded antipode, the comodule action `k̃`, and the
      **negative twisted antipode `S'₋`** recursion (BPHZ).
- [ ] **Structures** — `RegularityStructure (T, T⁺, Δ, Δ⁺)` and
      `RenormalizationStructure (U, U⁻, δ, δ⁻)`; the (symbolic) **BHZ character** `k^ζ = h^ζ∘S'₋`.
- [ ] **Rule completion** (BHZ Prop 5.21) — close the structural rule under the contractions
      `Δ⁻` performs (replacing the current generate-then-filter shortcut).
- [ ] **`G⁻` vs `G⁻_ad` toggle** (moved from Phase 2) — impose admissibility (`k(I_p τ)=0`,
      `k(X^n⋆τ)=0` in the precise model sense) to reduce the free family to `G⁻_ad`.
- [ ] Test invariants: coassociativity and cointeraction of the coproducts; homogeneity tables
      reproduced (e.g. the full KPZ basis count, not just the negative rows).

**Acceptance:** the regularity structure for KPZ/Φ⁴₃-after-DPD is built and its coproducts satisfy
coassociativity/cointeraction; `S'₋` reproduces known forest formulas on small trees.

---

## Phase 4 — deferred seams (optional) ⬜

**Goal:** the boundaries we intentionally left without implementations.

- [ ] **`NoiseLaw` + canonical BPHZ values** — Gaussian/Wick expectations to evaluate
      `h^ζ(τ) = 𝔼[Π^ζ τ](0)` and hence the canonical `k_τ` (needs covariance kernels / integrals).
- [ ] **da Prato–Debussche pre-pass** — handle `β₀ ≤ −2` (Φ⁴₂, Φ⁴₃, sine-Gordon) by introducing
      one noise symbol per power to lift into a `β₀ > −2` structure.
- [ ] **Multi-index symbol basis** — a second `Symbol` implementation (Linares–Otto–Tempelmayr)
      plugging into the same generic algebra, validating the basis seam.
- [ ] **Analytic / numerical export** — emit a built structure (trees, homogeneities, coproducts)
      for an external analytic or numerical consumer.

---

## Output & rendering 🔨

**Goal:** emit *all useful information* — the renormalized family, and **every divergent tree**
drawn nicely — across terminal / LaTeX / Markdown / JSON. Full spec, information inventory, tree
notation, and report layout in [`notes/output.md`](notes/output.md). Principle: reuse SymPy's
pretty/LaTeX printers; we own only the **tree** drawer and the report assembler. No new runtime
dep (`forest` is a `.tex` concern; ANSI colour is `isatty`-gated; `rich` only if asked).

- [x] **O1 (with Phase 1, now):** `render/tree.py` — three faithful renderings of a
      `DecoratedTree`: linear shorthand (`●·𝓘ₓ[Ξ]²`), 2D terminal drawing, LaTeX `forest`. A
      `render/report.py` text/markdown report covering inventory items 1–15 (echoed SPDE, domain,
      noise regularities, derived rule, the `𝓑_{<0}` table sorted by `|τ|` in the exact `ℚ+ℚκ`
      order, the assembled family). `summary()` now shares the prettifier. The full `𝓑_<0`
      (incl. `F(τ*)=0` non-contributors) is on `eq.all_trees`. *(D2 produces no zero-F trees in
      current scope, so "0 dropped" — the slot is wired and labelled for when red nodes arrive.)*
- [x] **O2:** `render/latex.py` — standalone `pdflatex`-ready document; JSON export (`to_json`,
      the clean handoff at the symbolic boundary).
- [ ] **O3 (with Phase 3):** wire the algebra section — pretty-print `Δτ`, `Δ⁻τ`, `S'₋τ` as
      tree⊗tree sums (pure reuse of the O1 tree drawer) and the structure's homogeneity spectrum.
- [ ] **O4 (with Phase 4):** the canonical-values column — replace each free `k_τ` with its BPHZ
      value once a `NoiseLaw` exists.

**Acceptance:** the gKPZ report reproduces the five trees and their homogeneity ordering; one
assert-based `demo()` per renderer is the backbone check.

---

## Cross-cutting (ongoing) ⬜

- [ ] Performance: memoize `F(τ*)` and `S(τ)` by canonical key once tree sets grow.
- [ ] Packaging/CI: lockfile, lint, a CI run of `uv run pytest`.
