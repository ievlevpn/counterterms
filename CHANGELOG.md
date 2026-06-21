# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to adhere to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- **Phase 2** — systems / vector unknowns, multiple noises, and general operator order &
  scaling (mostly `Signature` enrichment + validation).
- **Phase 3** — the coproducts `Δ, Δ⁺, Δ⁻`, the regularity/renormalization structures
  (`T, T⁺, U, U⁻`), the twisted antipode `S'₋`, and the symbolic BHZ character.
- **Phase 4 (deferred seams)** — `NoiseLaw` + canonical BPHZ constant values; a multi-index
  symbol basis; da Prato–Debussche pre-pass for `β₀ ≤ −2`.

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
