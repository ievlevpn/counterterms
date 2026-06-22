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

## Phase 3 — the algebraic-renormalization core 🔨

**Goal:** build the genuine regularity-structures machinery. This is the famous-hard part.
Modules: `trees/coproducts.py`, `structures.py`, `core/hopf.py` (the basis-agnostic layer).

- [x] **Extended-decoration trees** — red nodes with `𝔬 ∈ ℤ[β₀]` (reusing `Homogeneity`); the
      contraction `τ /^red φ`, `extended_homogeneity`. (`trees/tree.py`)
- [x] **Extraction-contraction `δ⁻`** (`Δ⁻`, tex 5636) + the group coproduct, **coassociative**
      (`(δ⁻⊗id)δ⁻ = (id⊗δ⁻)δ⁻`). Golden `δ∘ = 𝟙₋⊗∘ + ∘⊗●(β₀)`; stability invariants.
- [x] **Recentering `Δ`/`Δ⁺`** (tex 5613) + **comodule coassociativity**
      `(Δ⊗id)Δ = (id⊗Δ⁺)Δ`, `Δ⁺` coassociativity, counits, homogeneity stability.
- [x] **Negative twisted antipode `S'₋`** (tex 5034); `S'₋(∘) = −∘`.
- [x] **`δ⁺` (`D̄⁻`)** + the **cointeraction** `(id⊗Δ)δ = M¹³(δ⊗δ⁺)Δ` — **holds for the gKPZ
      class (β₀=−1) AND the singular β₀=−3/2 trees**. Fix: the e'-Taylor recentering runs over
      the full extraction boundary `∂(A,F)` (between-component edges too, not only φ→outside),
      and edges below a pre-existing red node stay internal. See
      [`notes/cointeraction_singular_noise.md`](notes/cointeraction_singular_noise.md) §8.
- [x] **Structures** — `RenormalizationStructure` (`δ`, `δ⁻`, `S'₋`) + the **symbolic BHZ
      character** `k = h∘S'₋` (`h` left symbolic). (`structures.py`)
- [x] **Generic Hopf layer** (`core/hopf.py`) — basis-agnostic `convolve` (character group
      law), connected-graded `antipode`, and `comodule_action` `k̃`, over plain
      coproduct/`mul`/`unit` maps. Reused verbatim for **both** `U⁻` (forests, `δ⁻`) and `T⁺`
      (blue trees, `Δ⁺`); the antipode axiom `S⋆id = η∘ε` is the test on each.
      `RegularityStructure.structure_antipode()` is the first consumer.
- [x] **`RegularityStructure (T, T⁺)`** — γ-bounded model basis (`generate_trees`, positive
      sector included), graded by homogeneity; `Δ : T → T⁺`, `Δ⁺` on `T⁺`; tested graded +
      triangular into `T⁺`, with the divergent subspace = the counterterms. (`structures.py`)
- [x] **Subcriticality check** (BHZ §5.5, `equation/rule.py`) — the rule-based criterion (every
      uncapped field edge strictly raises homogeneity over the most singular subtree, i.e.
      `β₀ > −order`) replaces the hardcoded `β₀∈(−2,0)` and the termination guard, and
      generalises to operator order. With it the saturating generator is a sound, terminating
      enumeration of `𝓑_{<0}`.
- [ ] **Formal rule completion** (BHZ Prop 5.21) — the smallest complete rule closed under `Δ⁻`;
      needed only for hand-supplied/adversarial rules (the derived `(f,g)`-rules saturate correctly).
- [x] **Renormalization group `G⁻`** — the character group of `U⁻` as a first-class object
      (`structures.py`): convolution group law (via `core/hopf`), counit unit, antipode inverse;
      a character carries the free constants `c_τ`, the renormalized family is its orbit, and
      `bhz_character` is the canonical element. Group axioms tested. `build_renormalization_group`.
- [ ] **`G⁻_ad` reduction** *(→ Phase 4)* — the admissible subgroup is a K-admissibility (analytic
      model) notion (kernel vanishing moments + the Π-map), not symbolic; a symbolic-only filter
      would be unsound (CLAUDE.md warning), so it is deferred to the model layer.

**Acceptance:** coassociativities + cointeraction (incl. singular β₀=−3/2) green as property
tests over the SPDE corpus; `S'₋` validated; Hopf group axioms tested. **132 tests pass (~4s).**

---

## Phase 4 — deferred seams (optional) ⬜

**Goal:** the boundaries we intentionally left without implementations.

- [~] **`NoiseLaw` + canonical BPHZ values** — *B1 done* (`renorm/scheme.py`): the Wick-pairing
      *combinatorics* of `h^ζ(τ)=𝔼[Π^ζτ](0)` (Isserlis), the mean-zero **parity** rule, and
      `RenormalizationStructure.canonical_character` (parity-reduced `h∘S'₋` — odd-noise trees'
      constants vanish; gKPZ: 3 of 5 are canonically 0). **B2 (evaluating the divergent
      integrals) remains** — the analysis wall.
- [ ] **`G⁻_ad` reduction** *(moved from Phase 3)* — the admissible subgroup `G⁻_ad ⊂ G⁻`; needs
      K-admissibility (kernel vanishing moments + the Π-map), so it belongs with the model layer.
- [x] **da Prato–Debussche pre-pass** (`equation/daprato.py`, `daprato_lift`) — lift a
      supercritical additive-noise **polynomial** SPDE `u=X+v` and Taylor-expand `P(X+v)`,
      replacing `X^k` by a Wick-power noise `:X^k:` of regularity `k·(β₀+order)`; the
      remainder equation is subcritical and `.renormalize()` handles it. **Unlocks Φ⁴₂ (3
      counterterms) and Φ⁴₃ (9)**, previously rejected. sine-Gordon (non-polynomial → Wick
      exponentials) remains out of scope.
- [ ] **Multi-index symbol basis** — a second `Symbol` implementation (Linares–Otto–Tempelmayr)
      plugging into the same generic algebra, validating the basis seam.
- [x] **Full-structure export** (`render/export.py`) — reconstructible, versioned JSON of the whole
      structure: canonical round-trippable trees (`tree_to_dict`/`tree_from_dict`), the graded `T`
      basis, the divergent trees with `S(τ)`/`F(τ*)`, the coproducts `Δ`/`δ⁻` as tensor sums, the
      signature/rule, and the BHZ character — for an external analytic/numerical consumer.

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
- [x] **O3 (with Phase 3 + Phase-4 B1):** the **canonical (BPHZ) renormalization** section
      (`canonical=True`, opt-in): each free `k_τ` at its canonical value `k_τ = h(S'₋τ)` for a
      centered Gaussian noise — **parity-reduced** (`canonical_character`), so odd-noise trees
      vanish (gKPZ: 3/5 → 0) and survivors collapse to short polynomials in the surviving `h(σ)`
      (KPZ: ~90 chars / 466 ms, vs 144 `h`-terms / 2.3 s for the raw `bhz_character`). An
      `h`-legend draws each surviving `σ`; the tree drawer renders **red contraction nodes** (with
      `o`-decoration) and blue `T⁺` roots. *Deliberately not rendered:* the raw `Δτ`/`Δ⁻τ`/`S'₋τ`
      tree⊗tree expansions — the character distills them, and they explode; they stay programmatic
      (`structures.py`, and the machine-readable `render/export.py`). The model-space spectrum is
      likewise available
      via `RegularityStructure.grades()` but not auto-rendered.
- [x] **O3.5 (B1 integrals):** the `h`-legend now spells out each surviving expectation as its
      **explicit (divergent) Wick integral** `h_ε(σ) = 𝔼[Π^ε σ](0) = Σ ∫ ∏∂^pK · ∏C_ε dz` (built
      by `renorm.scheme.expectation`, typeset by `render.report.expectation_latex`/`_str`) — the
      precise objects Track B2 must evaluate, with the **ε-regularization explicit** (`ξ_ε`, `C_ε`,
      constants diverging as `ε→0`). Text / markdown / LaTeX (`pdflatex`-verified). `K` is the
      abstract **singular kernel** of `L⁻¹` (Hairer's `K` in `K̄ = K + R`: the diagonal-singular,
      compactly-supported part of the Green's function `K̄`; it explodes on the diagonal, which is
      why the integrals diverge — tex 2105, 5683), *not* the Green's function itself.
      **Validity domain (rigorous, not ad hoc):** `expectation` is the *complete* `𝔼[Π^ζσ](0)`
      only for **bare** σ — no polynomial `X^n` node-decoration and no red contraction node —
      because `Π(X^n)(y)=y^n` (tex 1809), so a root `X^n` forces `h=0`, internal `X^n` add
      polynomial factors, and a red node is a contracted subtree carrying the extended
      `o`-decoration; all three need the full canonical model `Π^ζ` (Track B). `expectation`
      refuses non-bare σ (after the always-valid parity rule: odd noise count ⇒ 0), and the
      render leaves them symbolic. So `S'₋(τ)∈ℝ[U]` (tex 5028) and `k^ζ=h^ζ∘S'₋` (tex 5060) stay
      correct symbolically; only the *bare* h's get an explicit integral.
- [ ] **O4 (with Phase 4):** the numeric column — substitute the `h(σ)` values from a `NoiseLaw`
      (Wick) so `k_τ` becomes a number, not a symbol. **Needs Track B2** (evaluating the O3.5
      integrals).

**Acceptance:** the gKPZ report reproduces the five trees and their homogeneity ordering; one
assert-based `demo()` per renderer is the backbone check.

---

## Cross-cutting (ongoing) 🔨

- [x] **Packaging/CI** — `.github/workflows/ci.yml` runs `ruff check regstruct` + `uv run
      pytest` on push/PR; `ruff` pinned in the dev group (`E741` ignored — `l`/`r` are the
      left/right tensor legs). Lockfile present.
- [x] **Generator backstop** — `generate_trees` now **fails fast** with a clear `RuntimeError`
      past `_MAX_POOL=5000` trees (largest legitimate case ≈191), replacing the silent
      `assert guard<50` (stripped under `python -O`). Trips on the fractional/high-order +
      quadratic-gradient blow-up (see below) instead of hanging.
- [ ] **Performance: memoize `homogeneity` / `F(τ*)` / `S(τ)`** by canonical key. Would speed
      large-but-tractable sets; will *not* tame the fractional-order explosion (≈25×/iteration —
      genuinely astronomical, not a perf bug). Do when a real case needs the headroom.

### Known limitation — fractional/high operator order

A subcritical rule can still have an **intractably large** negative-tree set: at fractional or
high operator order the Schauder gain `m−|p|_𝔰` shrinks, so many trees stay singular. gKPZ at
order 3/2 (the `(∂u)²` term) explodes past the backstop in a couple of iterations. The engine now
refuses it with a pointed error rather than hanging; mild nonlinearities (gPAM at order 3/2 → 5
counterterms) are fine. Taming it for real (tighter `std<0` bound for counterterm-only runs,
memoization) is deferred — see `notes/swapping_the_operator.md`.
