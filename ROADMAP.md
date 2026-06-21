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

## Phase 2 — generality ⬜

**Goal:** widen the input class with (mostly) `Signature` enrichment + validation; little new
algorithmic structure, by design.

- [ ] **Systems / vector unknowns** — multiple components, each a sector with its own planting
      operator `I^{(a)}`. Equation/component identity rides the **edge type** `𝔗_e`, not the node
      type. `F` becomes a tuple `(F_a)` with cross-component partials `∂_{u_j}F_i`.
- [ ] **Multiple noises** — one `∘_j` node type per noise; base map `{●:g, ∘_j:f_j}`.
- [ ] **General operator order & scaling** — carry `m` and `𝔰` as first-class inputs in the
      homogeneity recursion (`|I_p τ| = |τ| + (m − |p|_𝔰)`), with scope warnings outside the
      proven second-order parabolic regime; mixed-order systems.
- [ ] **`G⁻` vs `G⁻_ad` toggle** — optionally impose the admissibility constraints
      (`k(I_p τ)=0`, `k(X^n⋆τ)=0`) that drop integral/polynomial-multiple counterterms.
- [ ] Golden tests for a coupled system and a multi-noise example.

**Acceptance:** a 2-component system and a 2-noise equation each produce the correct coupled
counterterms; scalar/single-noise results are unchanged.

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

## Cross-cutting (ongoing) ⬜

- [ ] Pretty/LaTeX tree rendering grouped by homogeneity (mirror the paper's tables); promote the
      `summary()` helper into a `render/` package.
- [ ] Performance: memoize `F(τ*)` and `S(τ)` by canonical key once tree sets grow.
- [ ] Packaging/CI: lockfile, lint, a CI run of `uv run pytest`.
