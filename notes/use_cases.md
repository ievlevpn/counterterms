# What `counterterms` can actually do — real use cases

An honest map of what this library computes and the genuine questions it answers.
Every example below was verified by running it; outputs are real.

## The precise boundary

It computes the **algebra and combinatorics** of BPHZ renormalization for subcritical
singular SPDEs, and it stops exactly at the analysis/probability wall. It tells you
**what** to renormalize and the **structure** of the renormalization — never the
**numeric values** of the constants, and never whether anything converges.

---

## What it genuinely solves

### 1. Derive the renormalized equation for an arbitrary in-scope SPDE
The core use. The BCCH formula by hand is error-prone — the symmetry factors `S(τ)`,
the Υ-map with its non-commuting `D^n` and `∂_p`, the tree enumeration. The library
returns the exact counterterm family. Example (gKPZ in **d=2**, no worked example in
the literature): 6 counterterms, each with its explicit elementary differential —
`f(u)` (S=1), `2·f(u)²g(u)` (S=2), `f(u)f'(u)`, `u_{∂x}·f'(u)`, `u_{∂y}·f'(u)`,
`2·u_{∂x}·f(u)g(u)`. The real value is for **new/variant equations, systems, and
multiple noises** — exactly where there is no paper to copy and hand derivation is
hardest.
*Who:* anyone studying a specific subcritical SPDE who needs its counterterms (for
analysis, or to set up a numerical scheme).

### 2. Subcriticality triage and the criticality boundary
Decide renormalizability-in-framework *before* investing in analysis, and map the
boundary as you vary β₀, dimension, operator order. Verified (2nd-order, single noise):

    β₀ = −1     subcritical → renormalizable here
    β₀ = −3/2   subcritical → renormalizable here
    β₀ = −2     REJECTED (supercritical → da Prato–Debussche)
    β₀ = −5/2   REJECTED (supercritical → da Prato–Debussche)

*Who:* anyone classifying a family of models, or checking "is this even in the
renormalizable regime?"

### 3. A computational oracle for the algebra itself
Compute `Δ`, `δ⁻`, `Δ⁺`, the twisted antipode `S'₋`, the BHZ character, the
cointeraction — **and verify the identities**. A real research testbed: this is how the
β₀=−3/2 cointeraction bug was found and the fix verified — a computation infeasible by
hand (see `cointeraction_singular_noise.md`). Concretely, the **BHZ constant decomposes
symbolically** into the elementary expectations `hᵢ`; for gKPZ:

    |τ|=−1−κ (1 node):  k(τ) = −h0
    |τ|=−κ   (1 node):  k(τ) = h0·h2 − h1
    |τ|=−κ   (2 nodes): k(τ) = h0·h4 − h3
    |τ|=−2κ  (2 nodes): k(τ) = h0²·h10·h2 − h0²·h9 − h0·h1·h10 − h0·h2·h7 + h0·h6 + h0·h8 + h1·h7 − h5
    |τ|=−2κ  (3 nodes): k(τ) = 2·h0²·h10·h4 − h0²·h13 − 2·h0·h10·h3 + 2·h0·h12 − 2·h0·h4·h7 − h11 + 2·h3·h7

That polynomial *is* the negative-renormalization recursion made explicit. Supply the
`hᵢ` from a noise law (computed elsewhere) and you get the actual constant.
*Who:* theorists developing or checking regularity-structures algebra.

---

## Real but bounded

### 4. Input to renormalized numerical schemes
The symbolic counterterm structure (which terms, their `F(τ*)`, `S(τ)`) is precisely
what a counterterm-subtracting discretization needs. You supply the scheme-dependent
constants from your own discrete computation; the library does the symbolic bookkeeping
half. *Boundary:* it does not compute those constants.

### 5. Complexity / scaling studies
Quantify how the renormalization grows with roughness — gPAM `f(u)ξ`, d = 1:

    β₀ = −1     3 counterterms   homogeneities {−1−κ, −κ, −2κ}
    β₀ = −3/2   11 counterterms  homogeneities {−3/2−κ, −1−2κ, −1/2−κ, −1/2−3κ, −2κ, −4κ}

and near criticality (β₀→−2) the tree set explodes (generation fails fast with a
`RuntimeError` past the 5000-tree backstop at β₀=−1.9 — a faithful reflection of
the math, not a bug).
*Who:* anyone gauging "how hard is this equation."

### 6. Combinatorial skeleton for an analytic proof / reproducing the literature
Before the hard estimates and model construction for a new equation, you need the exact
tree set, homogeneities, and counterterms — the library produces that skeleton exactly,
golden-tested against the paper's tables (gKPZ, KPZ homogeneity table, gPAM).

The companion object `build_regularity_structure(spde)` exposes the graded model space
`T = ⊕_α T_α` and `T⁺`; gKPZ spectrum: `{−1−κ:1, −κ:2, −2κ:2, 0:1}`.

---

## What it cannot do (the wall, explicitly)

- **No numeric renormalization constants** — no Gaussian/Wick expectations, no
  covariance integrals. The `hᵢ` / `k(τ)` *numbers* are not computed.
- **No analysis** — no model construction, no estimates, no convergence/well-posedness.
- **No solving or simulation.**
- **Out of scope:** β₀ ≤ −order for *non-polynomial* noise couplings (sine-Gordon — needs
  Wick exponentials); `g` beyond quadratic in `∂u` or `|p|_𝔰>1`; quasilinear /
  non-parabolic operators. *(Polynomial additive-noise supercritical equations — Φ⁴₂, Φ⁴₃ —
  **are** handled, via the built `daprato_lift` pre-pass; see `equation/daprato.py`.)*
- `G⁻_ad` reduction and formal rule completion are deferred (model-dependent / niche).

---

## The one sentence

It is the **symbolic front-end of singular-SPDE renormalization**: for any in-scope
equation it gives — exactly and automatically — the renormalized equation, the
divergent-tree structure, and the Hopf-algebraic machinery (coproducts, antipode, BHZ
character, the renormalization group), but it deliberately stops where the probability
and analysis begin. Its sweet spot is a researcher or numericist who needs the
*structure* of the renormalization for a new/non-textbook equation and wants it derived
correctly in seconds rather than by hand over days.

**On "free constants":** not a gap. A solution to a singular SPDE *is* the family
indexed by the renormalization group, so emitting free `c_τ` is the complete symbolic
answer; only the *canonical* (BPHZ) numeric values need the Phase-4 model layer.
