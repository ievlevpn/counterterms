# Phase 4 — plan

Phases 1–3 are the **symbolic** renormalization engine (renormalized equation, the
divergent-tree structure, the full Hopf-algebraic machinery). Phase 4 is where the
library either **broadens its symbolic reach** or **crosses into the analytic /
probabilistic half** — and these are genuinely different commitments. This plan
separates them honestly and orders the work by *value ÷ tractability*.

The central fact to keep in view: **the renormalized family with free constants is
already complete.** Phase 4 adds (a) more equations it applies to, (b) the *numeric*
values of the constants, (c) interop. Only (b) crosses the analysis wall.

```
                         Phase 4
            ┌───────────────────────┴───────────────────────┐
        TRACK A (symbolic)                        TRACK B (analytic / probabilistic)
   stays in the comfort zone                  crosses the wall — numbers & integrals
   ─ A1 da Prato–Debussche lift              ─ B1 NoiseLaw + canonical model + Wick comb.
   ─ A2 full-structure export                ─ B2 canonical BPHZ integrals  ⚠ research-grade
   ─ A3 multi-index basis (optional)         ─ B3 G⁻_ad reduction (needs the model)
```

---

## Track A — symbolic extensions (tractable, high value)

### A1. da Prato–Debussche pre-pass  — *unlock β₀ ≤ −2 (Φ⁴₂, Φ⁴₃, …)*  **[recommended first]**

**Goal.** Make the headline use case work for the *famous* equations currently rejected
as supercritical (Φ⁴₂, Φ⁴₃, …) — by a purely symbolic transformation, then reuse the
existing pipeline. No integrals.

**Math.** For `(∂_t − Δ)u = −u^m + ξ` with `β₀ ≤ −2`, split `u = X + v` where `X = K∗ξ`
is the stochastic convolution (regularity `α_X = β₀ + order`, e.g. Φ⁴₃: `β₀=−5/2 ⇒
α_X=−1/2`). Expand the nonlinearity in `X+v`; replace each pure power `X^k` by a
**Wick power** `:X^k:`, a *new noise symbol* of regularity `k·α_X` (e.g. Φ⁴₃:
`:X²:∈C^{−1}`, `:X³:∈C^{−3/2}`). The resulting equation for `v`,
```
(∂_t − Δ)v = −v³ − 3v²·X − 3v·:X²: − :X³:        (Φ⁴₃)
```
has roughest noise `:X³:∈C^{−3/2−} ⇒ β₀_eff = −3/2 > −2`: **subcritical**, so the
existing engine renormalises it. The coefficients (`−3v²`, `−3v`, `−1`) are polynomials
in `v` — already handled (a polynomial is a fine `f`).

**Design.** `equation/daprato.py`: `lift(spde) -> SPDE` — detect `β₀ ≤ −order`, do the
`u=X+v` substitution symbolically (SymPy), introduce `Noise` symbols for the Wick powers
with computed regularities, and return the `v`-equation. `renormalize` gains an opt-in
`daprato=True` (or auto-detect on the supercritical rejection). Output: the renormalised
`v`-equation, free constants — exactly the library's current contract.

**Feasibility.** High for **polynomial** nonlinearities (Φ⁴ family). **Caveats marked:**
the Wick-power noises are *correlated* (all from one `ξ`); for the counterterm *structure*
(trees, `F(τ*)`) correlation is irrelevant and the existing multi-noise machinery is
exactly right — only the *canonical values* (B-track) see the correlation. **sine-Gordon**
(`sin(βu)`) needs Wick *exponentials* `:e^{±iβX}:` with charge-dependent regularity
(`β²` thresholds) — a harder, separate step; scope polynomial first.

**Boundary.** Delivers the symbolic renormalised equation for Φ⁴₃ etc.; does **not**
give the Wick-power renormalisation *values* (those are B-track).

**Value: HIGH. Effort: MEDIUM. Risk: LOW (stays symbolic).**

---

### A2. Full-structure export  — *interop*

**Goal.** Serialize the whole built object (trees, homogeneities, `S(τ)`, `F(τ*)`, the
coproducts `Δ/δ⁻/Δ⁺`, the `T`/`T⁺` bases, the symbolic BHZ character) to JSON/HDF5 for an
external analytic-estimate or numerical-scheme consumer. Extends the existing
`render(\"json\")` (which already emits the renormalised family) to the algebra.

**Design.** `render/export.py`: a stable schema (versioned). Trees → canonical-key +
shorthand; coproducts → lists of `(coeff, left, right)`; homogeneities → `(std, kap)`.

**Feasibility/Value/Effort.** Easy; useful glue for anyone wanting to *use* the structure
downstream. **Value: MEDIUM. Effort: LOW. Risk: LOW.**

---

### A3. Multi-index symbol basis (Linares–Otto–Tempelmayr)  — *architectural, optional*

**Goal.** A *second* `Symbol` implementation — the multi-index/jet basis — plugging into
the generic `core/hopf` layer, validating the basis seam (the payoff promised by the
generic Hopf work).

**Design.** `multiindex/` package: a `MultiIndex` symbol (`homogeneity`, `canonical_key`)
+ its coproducts; the `core/hopf` operations are reused verbatim.

**Feasibility.** Substantial (a whole alternative combinatorics with its own coproducts).
Payoff is mostly *architectural validation* + research interest, not new user output.
**Value: LOW–MEDIUM (niche). Effort: HIGH. Risk: MEDIUM.** Do last, if at all.

---

## Track B — the analytic / probabilistic crossing

This is the conceptual completion (free constants → numbers) but it is a **different
kind of work**: Gaussian/Wick probability and singular-integral analysis, not tree
algebra. Split B1 (symbolic-tractable) from B2 (research-grade).

### B1. NoiseLaw + canonical model + Wick combinatorics  — *the symbolic half of the constants*

**Goal.** Express the canonical constant `k_τ = h(S'₋(τ))` with `h(σ)=𝔼[Π^ζσ](0)` as far
as is **symbolic**: build the canonical model `Π^ζ` and reduce `𝔼[Π^ζσ](0)` to an explicit
**sum over Wick pairings**, each a product of covariances × kernel integrals — *without yet
evaluating the integrals*.

**Math/design.**
- `NoiseLaw` (architecture §3.10, the named-but-unbuilt seam): the covariance `𝔼[ξ(x)ξ(y)]
  = C(x−y)` (e.g. `WhiteNoise → δ`, or a mollified `C_ε`).
- `CanonicalModel`: the recursive map `Π^ζ(Ξ)=ξ`, `Π^ζ(I_p τ)=∂^p K∗Π^ζ(τ)`,
  `Π^ζ(τσ)=Π^ζτ·Π^ζσ`, `Π^ζ(X^k)=(·)^k` (tex 2300 admissibility).
- **Wick expansion**: `𝔼[Π^ζτ](0) = Σ_{pairings of the noise leaves} ∏C · ∏(kernel
  factors)` → for each pairing a labelled graph whose value is an integral over the
  internal vertices. This step is pure combinatorics (Isserlis/Wick) — tractable.
- Compose with `S'₋` (we have it) to get `k_τ` as a (symbolic) combination of these
  pairing-integrals — the BPHZ-renormalised value, with the subdivergence subtractions
  supplied by the twisted antipode.

**Boundary.** Produces `k_τ = (explicit renormalised integral expression)` and its
**divergence degree** (already known from the homogeneity). Does not evaluate the
integrals.

**Feasibility.** Wick combinatorics: tractable. **Value: MEDIUM–HIGH (bridges to numbers).
Effort: MEDIUM–HIGH. Risk: MEDIUM.** Depends on nothing but the existing `S'₋`.

### B2. Canonical BPHZ integrals  — ⚠ **research-grade**

**Goal.** Actually *evaluate* the renormalised integrals from B1 to numbers.

**The wall.** These integrals are **singular and divergent** (UV at coinciding points) —
that is the entire reason renormalization exists. Finiteness comes only after the `S'₋`
subtractions; evaluating them means either closed forms (rare, simplest trees only) or
careful numerical integration of regularised kernels with the counterterms (`ε→0`), which
is delicate and itself a topic of active analytic/numerical research.

**Scope realistically.** (i) closed forms for the few simplest trees / specific kernels;
(ii) a numerical backend for *mollified* covariances `C_ε` returning the
regularisation-dependent constant + its divergence rate; (iii) the divergence *degree*
exactly (free, from homogeneity). The *general* exact constant is **out of reach** without
a dedicated singular-integral engine.

**Feasibility.** Partial only. **Value: HIGH (the dividend). Effort: VERY HIGH. Risk:
HIGH.** Treat as a research sub-project, not a sprint. Depends on B1.

### B3. G⁻_ad reduction  — *the admissible subgroup*

**Goal.** Reduce the free family to the admissible subgroup `G⁻_ad ⊂ G⁻` — drop/relate the
trees a K-admissible model renormalises trivially.

**Math.** K-admissibility (`Π(I_nτ)=∂^nK∗Πτ`, tex 2300) forces relations among the
renormalisation of `I_nτ` and `τ`; `G⁻_ad` is the characters respecting them. Once the
canonical model (B1) exists, these relations are computable; a symbolic-only filter is
unsound (CLAUDE.md), so this **depends on B1**.

**Feasibility/Value.** Medium once B1 lands; refines the output (fewer/related constants).
**Effort: MEDIUM (post-B1). Risk: MEDIUM.**

---

## Recommended ordering & milestones

| # | Item | Track | Value | Effort | Depends on |
|---|------|-------|-------|--------|-----------|
| 1 | **A1 da Prato–Debussche** (polynomial) | A | High | Med | — |
| 2 | **A2 full-structure export** | A | Med | Low | — |
| 3 | **B1 NoiseLaw + canonical model + Wick** | B | Med–High | Med–High | `S'₋` (done) |
| 4 | **B3 G⁻_ad** | B | Med | Med | B1 |
| 5 | **B2 canonical integrals** ⚠ | B | High | Very high | B1 |
| 6 | **A3 multi-index basis** | A | Low–Med | High | core/hopf (done) |

**Milestone M-A (symbolic reach):** A1 + A2 → the library handles Φ⁴₃ symbolically and
exports the full structure. Self-contained, high value, low risk. *Do this first.*

**Milestone M-B (the numbers):** B1 → B3 → (B2). B1 is the honest, tractable step toward
canonical values (the constants as explicit renormalised integrals); B2 is where it
becomes a research project. Stop at B1+B3 unless a real consumer needs the numbers.

---

## What stays out, even after Phase 4
- Analytic estimates / convergence / well-posedness proofs (the theory's analysis half).
- Construction of the *model* as a probabilistic object and the reconstruction theorem.
- Solving / simulating the SPDE.

These are not library goals; the library's job ends at *what to renormalise and the
structure of the renormalisation*. Phase 4's honest ceiling is: **every subcritical
(and DPD-liftable) equation's renormalised form symbolically (A1), exported (A2), with the
constants reduced to explicit renormalised integrals (B1) and, where evaluable, to numbers
(B2).**

## Suggested next action
**A1 (da Prato–Debussche, polynomial nonlinearities).** It is the one Phase-4 item that
delivers new *user-facing* capability (Φ⁴₂/Φ⁴₃ — the equations everyone actually asks
about), stays entirely symbolic, reuses the whole existing engine, and carries low risk.
