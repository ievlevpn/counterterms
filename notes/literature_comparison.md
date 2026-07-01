# Comparison of `counterterms` output against known cases

**Date:** 2026-06-23
**Purpose:** check the engine's output against (Tier 1) the source paper Bailleul & Hoshino
(arXiv:2006.03524, cited by `.tex` line) and (Tier 2) the published renormalized equations for the
classical models — KPZ, PAM, gPAM, Φ⁴₂, Φ⁴₃ — using only results that actually exist in the
literature.
**Method:** ran the engine on each equation (`spde.renormalize()`, `daprato_lift(...).renormalize()`),
listed the counterterms `(τ, |τ|, S, F(τ*), k/S)`, then compared to the paper line-by-line and to the
literature structurally. Not committed.

> **Key framing.** The engine emits the **free-constant family** — *every* tree that could carry a
> counterterm, for *any* renormalization character `k ∈ G⁻`. A published renormalized equation is the
> **canonical (BPHZ, symmetric-noise)** specialization. They agree once the appropriate reductions are
> applied. The engine performs three (noise parity; root `Xⁿ` ⇒ 0, tex 5083; pure-kernel total
> derivative ⇒ 0); it does **not** yet apply **spatial-reflection** symmetry — see Caveat 1.

---

## Tier 1 — the paper (primary oracle): exact matches

| Case | Paper | Engine | Verdict |
|---|---|---|---|
| **gKPZ**, ζ∈C^{−1−κ} | 5 counterterms (tex 6004–6012): `k(∘)f`, `k(∘1)f'∂ₓu`, `2k·f g ∂ₓu`, `k·f f'`, `k·f²g` | the same 5, with `F(τ*)` and the factor‑2/`S(τ)` bookkeeping identical | **exact** |
| **KPZ tree table**, β₀=−3/2−κ | 43 strongly-conforming trees, rows of size (1,2,6,2,23,9) at homogeneities (−3/2−κ, −1−2κ, −1/2−3κ, −1/2−κ, −4κ, −2κ) (tex 6028–6063) | 43 trees, identical per-row counts and homogeneities | **exact** |
| **Φ⁴₃ da Prato lift** | `(∂t−Δ)v = −v³ − 3v²X − 3vX² − X³`, `Xⁿ ∈ C^{(−n/2)⁻}`, reduced β₀=(−3/2)⁻ (tex 2031–2034) | lifted remainder `−v³ − 3X1·v² − 3X2·v − X3`, with X1,X2,X3 at std −1/2, −1, −3/2 | **exact** |

The gKPZ identity is the strongest single check: the factor 2 *survives* on `●·𝓘ₓ[Ξ]` (S=1) and is
*cancelled* on `●·𝓘ₓ[Ξ]²` (S=2), exactly as the paper writes `2k(τ₃)` vs plain `k(τ₅)`.

---

## Tier 2 — literature: the classical renormalized equations

For each: the published equation, the engine's free family, and what survives after reduction.

### KPZ — `∂t h = Δh + (∂x h)² + ξ`, ξ∈C^{−3/2−κ}
**Published** (Hairer, *Solving the KPZ equation*): a **single diverging constant**,
`∂t hε = ∂x²hε + λ(∂x hε)² − λCε + ξε`, `Cε ≈ ε⁻¹∫φ²`. **No gradient counterterm** (symmetric mollifier).

**Engine** (8 trees):

    Ξ                                  |Ξ|=−3/2−κ   S=1  F=1        k₀     -> 0  (odd noise parity)
    ●·𝓘ₓ[Ξ]                            −1/2−κ       S=1  F=2 h_x     k₁     -> 0  (parity)
    ●·𝓘ₓ[Ξ]²                           −1−2κ        S=2  F=2        k₂     -> const   (= −h0)
    ●·𝓘ₓ[●·𝓘ₓ[Ξ]²]·𝓘ₓ[Ξ]              −1/2−3κ      S=2  F=4        k₃     -> 0  (parity)
    ●·𝓘ₓ[●·𝓘ₓ[Ξ]]·𝓘ₓ[Ξ]               −2κ          S=1  F=4 h_x     k₄     -> −h3
    ●·𝓘ₓ[●·𝓘ₓ[Ξ]²]                     −2κ          S=2  F=4 h_x     k₅     -> −h1
    ●·𝓘ₓ[●·𝓘ₓ[Ξ]²]²                    −4κ          S=8  F=8        k₆     -> const
    ●·𝓘ₓ[●·𝓘ₓ[●·𝓘ₓ[Ξ]²]·𝓘ₓ[Ξ]]·𝓘ₓ[Ξ] −4κ          S=2  F=8        k₇     -> const

The two **`∂h` (drift) counterterms** `k₄=−h3`, `k₅=−h1` vanish for a symmetric mollifier: `h1`, `h3`
each carry **three** `∂ₓK` factors, so under spatial reflection x→−x the integrand is odd ⇒ `h1=h3=0`.
The constants `k₂,k₆,k₇` combine into the single `−Cε`. **Result: `∂t h = ∂x²h + (∂x h)² − C + ξ`** —
matches Hairer **after the spatial-reflection reduction** (which the engine does not yet apply; Caveat 1).

### 2D PAM — `∂t u = Δu + u·ξ`, ξ space white noise C^{−1−κ}
**Published**: `∂t u = Δu + u(ξ − Cε)` — a constant times `u`.

**Engine** (4 trees): `Ξ`(F=u, k₀→0 parity), `Xₓ·Ξ`/`X_y·Ξ`(F=u_x/u_y, k₁,k₂→0 **root Xⁿ**),
`Ξ·𝓘[Ξ]`(F=u, k₃→const). **Result: `u·const`.** ✅ The gradient terms vanish by root‑`Xⁿ`=0 (a
reduction the engine *does* apply); bare `Ξ` by parity.

### 2D gPAM — `∂t u = Δu + f(u)ξ`
**Published** (BCCH; Gerencsér–Hsu confirm the `f,f'` form): `∂t u = Δu + f(u)ξ − Cε·f'(u)f(u)`.

**Engine** (4 trees): `Ξ`(F=f, →0 parity), `Xₓ·Ξ`/`X_y·Ξ`(F=f'·u_x/u_y, →0 root‑`Xⁿ`),
**`Ξ·𝓘[Ξ]`(F=f·f', k₃)** survives. **Result: `−C·f(u)f'(u)`.** ✅ exact structure (and no gradient
counterterm, correctly — gPAM has no gradient nonlinearity).

### Φ⁴₂ — `(∂t−Δ)φ = −φ³ + ξ`, d=2 spacetime WN, via da Prato lift
**Published**: `∂t φ = Δφ − φ³ + 3Cε φ + ξ` (Wick / mass renormalization).

**Engine** (lifted, 3 trees): `Ξ_X1`(F=−3v², k₀), `Ξ_X2`(F=−3v, k₁), `Ξ_X3`(F=−1, k₂). The **mass term
∝ v** (`k₁`) is the renormalization; `k₀` (∝`:X:`, mean zero) and `k₂` (∝`:X³:`) drop. ✅ structure.

### Φ⁴₃ — `(∂t−Δ)φ = −φ³ + ξ`, d=3 spacetime WN, via da Prato lift
**Published** (Jagannath–Perkowski): `(∂t−Δ)φ = −φ³ + ∞·φ + ξ` — **two** diverging constants, both ∝ φ (mass).

**Engine** (lifted, 9 counterterms from 13 divergent trees): gradient `X_·Ξ_X2`(F=−3v_·, k₂,k₃,k₄ → 0 **root Xⁿ**); **mass ∝ v**
`Ξ_X2`(−3v), `Ξ_X1·𝓘[Ξ_X3]`(6v), `Ξ_X2·𝓘[Ξ_X2]`(9v) (k₁,k₆,k₈ — the two-constant mass structure);
`v²`/const `Ξ_X1`,`Ξ_X3`,`Ξ_X2·𝓘[Ξ_X3]` (k₀,k₅,k₇, parity/mean-zero). **Result: mass ∝ φ**, gradient
terms vanish. ✅ structure (two diverging mass constants).

The general renormalized-equation formula the engine implements (`ThmRenormPDEs`) is BCCH's.

---

## Caveats (where "match" carries an asterisk)

1. **Spatial-reflection reduction — now available as `reduced=True`** (2026-06-23). The default
   `canonical=True` view applies only the noise-independent reductions (parity, root‑`Xⁿ`=0,
   pure-kernel total-derivative=0), so it keeps `k₄=−h3`, `k₅=−h1` (the `∂h` drift terms) — correct
   for a *general* noise, since a non-symmetric mollifier genuinely produces a drift. The new
   `reduced=True` option additionally folds the **spatial-reflection** identity (`h(σ)=0` when the
   total spatial-derivative order `D=Σ|p_e|_space` is odd), collapsing canonical KPZ to Hairer's
   single constant. This identity is valid **only for a spatially-symmetric noise**, so it is gated
   on a `symmetric` flag (default `True` = white noise / symmetric mollifier); with `symmetric=False`
   it is **not** applied and the report says so (`reduction_assumes_symmetric_noise: false`). So the
   engine never claims the reflection reduction for an anisotropic noise.

2. **Φ⁴ canonical values are structure-only.** The da Prato lift turns the Wick powers `:Xᵏ:` into
   symbols `X1,X2,X3` that the engine treats as **independent** noises. The **free-constant family is
   correct** and matches the literature structure above, but the *numeric/canonical* constants for Φ⁴
   would need the genuine correlations among the `:Xᵏ:` (a single underlying Gaussian) — exactly the
   deferred Track B2. So Φ⁴ was compared by **structure**, not canonical values.

---

## Verdict

Every case that exists in the literature **agrees** with the engine: exactly against the paper (gKPZ
5 counterterms, KPZ 43-tree table, Φ⁴₃ da Prato reduction), and structurally against
KPZ/PAM/gPAM/Φ⁴₂/Φ⁴₃ once the appropriate symmetry/parity reduction is applied. The one genuinely
actionable gap is the **missing spatial-reflection reduction** (the KPZ drift terms), which sits
alongside the `Xⁿ`/total-derivative reductions.

## Sources

- M. Hairer, *Solving the KPZ equation* — https://www.hairer.org/papers/KPZ.pdf (renormalized KPZ:
  `∂t hε = ∂x²hε + λ(∂x hε)² − λCε + ξε`, single constant).
- Bruned, Chandra, Chevyrev, Hairer, *Renormalising SPDEs in regularity structures*, JEMS 2021 —
  https://arxiv.org/abs/1711.10239 (the general renormalized-equation formula; gPAM/KPZ examples).
- Jagannath, Perkowski, *A simple construction of the dynamical Φ⁴₃ model* —
  https://arxiv.org/abs/2108.13335 (`(∂t−Δ)φ = −φ³ + ∞·φ + ξ`, two diverging mass constants).
- Gerencsér, Hsu, *Variance renormalisation in regularity structures — 2d gPAM* —
  https://arxiv.org/abs/2602.17369 (the `f, f'` counterterm for gPAM).
- Hairer, Labbé, *A simple construction of the continuum parabolic Anderson model on R²* —
  https://arxiv.org/abs/1501.00692 (PAM `u(ξ−Cε)`, logarithmic renormalization).
- Source paper of this project: Bailleul & Hoshino, *A tourist's guide…* — arXiv:2006.03524
  (gKPZ tex 6004–6012; KPZ table tex 6028–6063; Φ⁴₃ da Prato tex 2026–2034).
