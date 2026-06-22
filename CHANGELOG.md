# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to adhere to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Canonical (BHZ) renormalization in the report** (`render(..., canonical=True)`,
  `eq.report/latex_document/to_json/save(canonical=True)`). Wires Phase 3 to the renderer: each
  free constant printed as its exact canonical value `k_τ = h(S'₋τ)` — a symbolic polynomial in
  the elementary expectations `h(σ) = 𝔼[Πσ](0)` — with an `h`-legend drawing every `σ`. Opt-in
  because the twisted antipode explodes for deep trees (KPZ: 144 `h`-terms). The forest drawer
  now renders **red contraction nodes** (square, with their `o`-decoration) and blue `T⁺` roots.
  The raw coproducts/antipode stay programmatic (`structures.py`) — the BHZ character distils them.
- **Output & rendering** (`regstruct/render/`, see `notes/output.md`). A full report of the
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
- **Phase 3** — the coproducts `Δ, Δ⁺, Δ⁻`, the regularity/renormalization structures
  (`T, T⁺, U, U⁻`), the twisted antipode `S'₋`, the symbolic BHZ character, rule completion,
  and the `G⁻`/`G⁻_ad` toggle.
- **Phase 4 (deferred seams)** — `NoiseLaw` + canonical BPHZ constant values; a multi-index
  symbol basis; da Prato–Debussche pre-pass for `β₀ ≤ −2`.

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
