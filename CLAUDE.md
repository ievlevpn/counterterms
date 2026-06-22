# CLAUDE.md — regularity_structures_symbolic

## What this project is

A symbolic package implementing Bailleul & Hoshino, *"A tourist's guide to regularity
structures and singular stochastic PDEs"* (arXiv:2006.03524 v3).

**Goal:** input a subcritical singular SPDE → output the *family of renormalized equations*
(the BCCH / `ThmRenormPDEs` formula): the original PDE plus one tree-indexed counterterm per
negative-homogeneity decorated tree, with free renormalization constants.

```
(∂_t − Δ + 1)u⁽ᵏ⁾ = f(u⁽ᵏ⁾)ζ + g(u⁽ᵏ⁾,∂u⁽ᵏ⁾) + Σ_{τ∈𝓑, |τ|<0} (k(τ)/S(τ)) F(τ*)(u⁽ᵏ⁾,∂u⁽ᵏ⁾)
```

**Status:** **Phases 1–2 implemented** — `SPDE → family of renormalized equations` with free
constants (rule → trees → S(τ) → Υ-map → assembly), now for **systems** (component on the edge
type; shared constants across components), **multiple noises**, and **general operator order**.
Golden-tested: gKPZ (exact 5 counterterms, tex 6004–6012), KPZ/gPAM, decoupled/coupled systems,
multi-noise, operator order, scope rejections. A `render/` package emits the full report (trees
drawn as shorthand / ascii / LaTeX-`forest`) in text/markdown/json/latex — see `notes/output.md`.
`uv run pytest` (163 tests, ~5s). Phase 3 complete & green: coproducts (cointeraction holds
**including β₀=−3/2**), `RegularityStructure (T,T⁺)`, the generic `core/hopf` layer,
subcriticality check, twisted antipode + BHZ character, renormalization group `G⁻`. Phase 4
started: `daprato_lift` (da Prato–Debussche) unlocks polynomial Φ⁴₂/Φ⁴₃ — see
`notes/phase4_plan.md`. Still deferred: canonical BPHZ *values* (Wick/integrals), `G⁻_ad`
reduction, formal rule completion. See `notes/architecture.md` §7 / `ROADMAP.md`.

## Layout

- `references/tourist_guide.{pdf,tex}` — the source paper (full LaTeX, 7852 lines). Cite by tex
  line number. The `.tex` is the fastest way to look things up (grep it).
- `notes/initial_plan.md` — authoritative for the **mathematics** (pipeline, conventions, scope,
  golden tests, reuse).
- `notes/architecture.md` — authoritative for the **module structure** (layered stack, the
  `Signature`-parametric design, interfaces, extension cookbook, phasing).
- `counterterms/` — package (Phase 1 modules populated): `core/{homogeneity,jets,signature}`,
  `trees/tree`, `equation/{dsl,generate}`, `renorm/{nonlinearity,equation}`, `render/{tree,report,latex}`, `api`. Phase-3
  modules (`core/{module,hopf,symbol}`, `trees/coproducts`, `structures/`) not yet created.
  `tests/` — layered by concern: `conftest.py` (the SPDE corpus + `ctx` fixture), `test_{homogeneity,
  core,trees}` (the non-negotiable conventions as unit tests), `test_pipeline` (generate + Υ-map),
  `test_coproducts` (the algebraic laws, parametrized over the corpus), `test_{structures,goldens,
  scope,render}`. Full target layout: `architecture.md` §5.
- Run: `uv run pytest`. Quick demo: `uv run python -u tests/test_render.py` (prints the
  renormalized gKPZ report).

## Environment

- Package manager: **`uv`** (not conda, not pip directly). Run Python as `uv run python -u`.
- Symbolic backend: **SymPy** (planned hard dependency). Keep tree combinatorics in a
  dedicated dict-of-trees / free-module representation; convert to SymPy only at the leaves and
  for final output (SymPy auto-flattening makes large tree sums slow/awkward).
- SageMath / Julia (`BSeries.jl`, `kauri`) are **oracles and design references only**, never
  runtime dependencies (Sage isn't `uv`-friendly).

## The pipeline (what to build)

`SPDE (DSL) → rule → strongly-conforming trees with |τ|<0 → S(τ) + F(τ*) → assembled family`.
**The MVP needs no Hopf coproducts or twisted antipode** — those only compute the *canonical*
constant values, which require probabilistic (Wick) input that is out of scope. Emit free
symbolic constants `c_τ = k(τ)/S(τ)`.

## Non-negotiable math conventions (getting these wrong is silent and fatal)

1. **`|Ξ| = β₀`** = the noise's Hölder regularity *directly* (user-supplied, negative). Do **not**
   use a "regularity = α−2" convention — it's an off-by-2. `|I_p τ| = |τ| + 2 − |p|_𝔰`,
   `|X^n| = |n|_𝔰`. Verified against the paper's example tables (tex 6004–6063).
2. **Compare homogeneities in the ordered ring** `ℚ + ℚ·κ` (κ>0 infinitesimal), never floats —
   critical trees sit at homogeneity `−kκ`.
3. Coefficient is **`k(τ)/S(τ)`**, never `k(τ)` alone. `S(τ) = n!·Πⱼ S(σⱼ)^{mⱼ}·mⱼ!`.
4. **Canonical tree isomorphism** is load-bearing (equality, dict keys, dedup, `S`). Key:
   `(node_type, n, sorted (edge_type, p, canonical(child)))`.
5. The Υ-map `F(τ*) = (Πᵢ F(τ_i*))·(D^n Πᵢ ∂_{p_i}) F(b*)` differentiates **all** slots of `g`
   (function and derivative args), `∂_{p_i}` before `D^n` (they don't commute). Base:
   `F(∘_j*)=f_j`, `F(●*)=g`, `F(red*)=0`. `D_i = Σ_k u_{k+e_i} ∂_k`.
6. `𝓑_{<0}` **includes bare primitives** `∘^n, ●^n` (the KPZ `k(∘)`, `k(∘1)` terms).
7. For **systems**, equation/component identity lives on the **edge type** `𝔗_e`, not the node
   type (node types stay `{●, ∘_j, red}`).

## Scope (enforce with explicit errors)

In scope: scalar **or systems**, single **or multiple** noises, 2nd-order parabolic `L`
(general operator order with a warning), `β₀ > −order` (rule-based subcriticality), `g` ≤
quadratic in `∂u` (Assumption D2), `|n|_𝔰 ≤ 1`. Covers gKPZ, KPZ, gPAM, PAM, coupled systems,
multi-noise. **Supercritical `β₀ ≤ −order` polynomial additive-noise equations (Φ⁴₂, Φ⁴₃) are
handled via `daprato_lift` (the da Prato–Debussche change of variables) → subcritical remainder.**

Reject (with clear messages): supercritical rules not liftable here (sine-Gordon /
non-polynomial — needs Wick exponentials); noise nonlinearities not affine in the noise; `g`
more than quadratic in `∂u`; singular derivative factors with `|n|_𝔰 > 1`; quasilinear /
non-parabolic `L`.

## Validation anchors (the paper is the oracle — no reference code exists)

Golden test = the gKPZ example's exact 5 counterterms (tex 6004–6012) and the KPZ homogeneity
table (tex 6028–6063). See plan §9. Build these as tests before trusting the engine.

## Working style

- Lazy/minimal (ponytail): reuse SymPy, build only the novel core, fewest files (plan §7),
  shortest working diff. Mark deliberate simplifications with `# ponytail:` comments.
- Non-trivial logic leaves one runnable check behind (assert-based `demo()`/`__main__` or a small
  `test_*.py`). The golden tests above are the backbone.
- This is research math: precision over speed. When unsure about a formula, **grep the `.tex`**
  and cite the line — don't guess.
- **Type annotations are required.** Every function/method in `counterterms/` carries full type
  hints — all parameters (except `self`/`cls`) and the return type (`-> None` for `__init__`);
  annotate dataclass fields too. Put `from __future__ import annotations` at the top of each
  module (annotations stay lazy strings) and import domain types under `if TYPE_CHECKING:` to
  avoid cycles. Use the project vocabulary: `Signature`, `DecoratedTree`, `Homogeneity`,
  `MultiIndex` (`=tuple[int,...]`), `SPDE`, `sympy.Expr`/`sympy.Symbol`, `Fraction`,
  `tuple[DecoratedTree, ...]` for forests, `TensorSum` for coproduct results. Prefer a precise
  domain type; use `object` over a wrong guess. **Tests (`tests/`) are exempt.** New code must
  keep this invariant — `uv run python -m ast`-clean and annotated.
