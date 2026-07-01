# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to adhere to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Phase 3 — the algebraic-renormalization core** (`trees/coproducts.py`, `core/hopf.py`,
  `structures.py`): extraction–contraction `δ`/`δ⁻`, recentering `Δ`/`Δ⁺`, the negative twisted
  antipode `S'₋`, the symbolic BHZ character `k = h∘S'₋`, `RegularityStructure (T, T⁺)`, the
  rule-based subcriticality check (`β₀ > −order`), and the renormalization group `G⁻`. The
  cointeraction holds over the corpus **including β₀ = −3/2**.
- **da Prato–Debussche lift** (`equation/daprato.py`, `daprato_lift`) — supercritical polynomial
  additive-noise equations (Φ⁴₂: 3 counterterms, Φ⁴₃: 9) via the subcritical remainder equation.
- **Canonical BPHZ scheme, symbolic half** (`renorm/scheme.py`): Wick-pairing combinatorics with
  within-noise-type Isserlis matching, the centered-Gaussian parity rule, and each surviving
  `h(σ)` spelled out as its explicit ε-regularized divergent integral in the report legend.
- **`reduced=True` report view**: the exact identities (parity, root `Xⁿ ⇒ 0`, pure-kernel total
  derivative, o-duplicate merging) substituted into the canonical constants, plus the
  **spatial-reflection identity** gated on a `symmetric` flag (default white noise) — collapses
  canonical KPZ to Hairer's single constant; never claimed for an anisotropic noise.
- **Full-structure JSON export** (`render/export.py`, `structure_json`) — round-trippable trees,
  graded basis, coproducts, antipode, characters.
- `FractionalHeat(dim, sigma)` operator; generator fail-fast backstop (`RuntimeError` past 5000
  pool trees instead of hanging).
- MIT license; CI (ruff + pytest); MkDocs documentation site (`mkdocs.yml`, `docs/`).
- Package renamed `regstruct` → `counterterms`; the paper untracked (copyright) with an arXiv
  pointer in `references/`.

### Fixed
- **Tree generation: node-decoration cap now scales with `−β₀`** instead of the fixed
  order-2-only `|n|_𝔰 ≤ 2` — at operator order > 2 with `β₀ ≤ −2` genuine counterterms (e.g.
  `Ξ^{(0,3)}` at order 4, β₀ = −3−κ) were silently missing from the family. Pools at order 2 are
  provably unchanged.
- **Tree generation: the DFS no longer prunes on intermediate partial sums** — a decorated node
  above the homogeneity bound can be pulled back under by several capped negative-contribution
  children (KPZ at β₀ = −7/4 was missing `●^{(0,2)}I'(Ξ)²` from the positive sector). Counterterm
  sets at order 2 unaffected.
- **Assumption-D2 total gradient bound** (`grad_budget`): in `d ≥ 2` a direction-mixing gradient
  nonlinearity could exceed total degree 2 in `∂u` per node in the *raw* tree basis (the
  renormalized equation was unaffected — spurious trees were Υ-zero).
- **δ⁺ comodule law**: removed the redundant force-internal-under-red rule that broke
  compatibility condition (b); the cointeraction is carried by the between-edge `e'` recentering
  alone. BHZ character byte-identical. Regression: `test_delta_plus_comodule`.
- **`p₋` projection**: only the *red* `●^{0,α}` maps to `𝟙₋` (a bare black `●` has homogeneity 0
  and maps to 0), per tex 5760.
- **`expectation` validity domain**: independent noises factorize (no cross-type Wick pairing);
  red contraction nodes are fine, only `Xⁿ` node decorations are refused.
- **LaTeX divergent-trees table overflow.** Two causes: (1) the `forest` trees carried a
  coordinate `baseline`, making TikZ report a wrong hbox width — now each tree is wrapped in
  `$\vcenter{\hbox{…}}$` (correct width, vertically centered on the math axis); (2) `longtable`
  computes its column widths from `.aux` and needs **two** compiler passes, but `save()` ran
  `pdflatex` once — so the trees overflowed their column on the single-pass PDF. `save()` now runs
  `pdflatex` twice. The `canonical=True` equation is also broken one term per line so its long
  combined RHS no longer runs past the right margin.

### Changed
- **Report's canonical section now shows the parity-reduced canonical (BPHZ) constants**
  (`canonical_character`) instead of the raw `bhz_character`. For a centered Gaussian noise the
  mean-zero parity rule makes odd-noise trees vanish — gKPZ collapses from 7-term polynomials in
  14 `h`-symbols to `k_3 = −h₀`, `k_4 = −h₁` with 3 of 5 constants `= 0`, and KPZ becomes
  tractable (~90-char expressions / 0.5 s instead of 144 `h`-terms / 2.3 s). Vanishing constants
  are labelled "(odd noise parity)"; the `h`-legend lists only the surviving expectations,
  renumbered contiguously. JSON key `bhz` → `canonical_bphz` (with a `vanishes` flag).
  `canonical=True` now also prints the **canonically renormalized equation** itself — the
  surviving counterterms substituted in (`(k_τ/S)·F(τ*)` with `k_τ` at its canonical `h(σ)`
  value), vanishing ones dropped, no free constants. (JSON: `canonical_family_latex`.)

### Added
- **Canonical (BHZ) renormalization in the report** (`render(..., canonical=True)`,
  `eq.report/latex_document/to_json/save(canonical=True)`). Wires Phase 3 to the renderer: each
  free constant printed as its exact canonical value `k_τ = h(S'₋τ)` — a symbolic polynomial in
  the elementary expectations `h(σ) = 𝔼[Πσ](0)` — with an `h`-legend drawing every `σ`. Opt-in
  because the twisted antipode explodes for deep trees (KPZ: 144 `h`-terms). The forest drawer
  now renders **red contraction nodes** (square, with their `o`-decoration) and blue `T⁺` roots.
  The raw coproducts/antipode stay programmatic (`structures.py`) — the BHZ character distils them.
- **Output & rendering** (`counterterms/render/`, see `notes/output.md`). A full report of the
  renormalized family in four formats via `eq.render(fmt)` / `eq.report()` /
  `eq.latex_document()` / `eq.to_json()`: echoed equation(s), domain/noise, the derived rule,
  **every** divergent tree (`𝓑_<0`, including `F(τ*)=0` non-contributors) with homogeneity,
  `S(τ)`, constant and `F(τ*)`, and the assembled per-component family. Phase-3/4 sections
  (coproducts, canonical values) render as honest placeholders.
- **Tree drawings** — three faithful renderings of a `DecoratedTree`: linear shorthand
  (`●·𝓘ₓ[Ξ]²`), a 2D terminal drawing, and a LaTeX `forest` snippet (`shorthand` / `ascii_art`
  / `forest`). Reuses SymPy's printers for all formulas; only the tree drawer is new.
- `RenormalizedEquation.all_trees` (the full `𝓑_<0`) and `.original_rhs(comp)`.

### Planned
- Numeric BPHZ constant values (Track B2 — evaluating the ε-regularized integrals).
- The `G⁻_ad` admissible-subgroup reduction (needs the analytic model layer).
- A multi-index symbol basis (Linares–Otto–Tempelmayr) on the generic-algebra seam.
- Formal rule completion (BHZ Prop 5.21) for hand-supplied rules.

## [0.2.0] — 2026-06-22

**Phase 2 — generality.** Systems, multiple noises, and general operator order.

### Added
- **Systems / vector unknowns.** Components share spacetime coordinates; each child **edge
  carries the component** of the kernel that plants it (node types stay `{●, ∘_j}`). The
  elementary differential threads the child equation index from the edges, and renormalization
  **constants are shared across components** (one `k_τ` per tree). A scalar problem is the
  one-component case.
- **Multiple noises** — one `∘_j` node type per noise.
- **General operator order** — per-component Schauder order in the homogeneity recursion, with a
  warning outside the proven 2nd-order parabolic regime (global scaling).
- `SPDE(equations=[(unknown, operator, rhs), …], noises=[…])` for systems; the scalar keyword
  form is unchanged.
- Tests (`tests/test_phase2.py`): decoupled system ⇒ per-component gKPZ; coupled system shares
  constants; multi-noise (8 counterterms); operator order changes homogeneities.

### Changed
- Jet variables are now component-indexed: `jet(comp, k)` (was `jet(k)`).
- `RenormalizedEquation` exposes `.per_component`; `.counterterms` returns component 0.

### Note
- The `G⁻`/`G⁻_ad` toggle was deferred to Phase 3 (admissibility is a model notion; the full
  free family is emitted as the safe superset).

## [0.1.0] — 2026-06-21

First implementation: **Phase 1 — `SPDE → family of renormalized equations`** with free
symbolic constants (no Hopf coproducts required).

### Added
- `core/` — ordered-ring homogeneities (`ℚ ⊕ ℚ·κ`, κ infinitesimal; never floats), jet
  variables, and the `Signature` (the typed vocabulary the library is parametric over).
- `trees/` — decorated rooted trees with canonical form, homogeneity, and symmetry factor
  `S(τ) = n!·Πⱼ S(σⱼ)^{mⱼ}·mⱼ!`.
- `equation/` — a SymPy-based DSL (`Unknown`, `Noise`, `Parabolic`, `SPDE`), rule derivation
  from the nonlinearity, and budget-pruned strongly-conforming tree enumeration.
- `renorm/` — the Υ elementary-differential map `F(τ*)` and the `RenormalizedEquation` output
  assembler (counterterm `= (k_τ / S(τ))·F(τ*)`).
- `api.renormalize` — the end-to-end pipeline.
- Tests: the gKPZ golden example (exact five counterterms, tex 6004–6012, including the
  factor-2 and `S = 2` cases), KPZ homogeneity rows, gPAM in `d = 2`, and scope rejections
  (`β₀ ≤ −2`, non-affine noise, cubic-in-`∂u`).
- Documentation: `notes/initial_plan.md` (mathematics), `notes/architecture.md` (modular
  design), `CLAUDE.md`, and the source paper under `references/`.

### Scope
Scalar unknown, single noise, second-order parabolic operator, `β₀ ∈ (−2,−1)`, `g` at most
quadratic in `∂u`, `|p|_𝔰 ≤ 1`. Covers gKPZ, KPZ, gPAM, PAM.
