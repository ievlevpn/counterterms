# The cointeraction residual for singular noise — research + diagnosis

**Status: diagnosed precisely; not yet fixed.** This documents the bug, the
authoritative resolution mechanism from the literature, a principled fix attempt
that *failed* (and exactly why), and the concrete path to a correct fix.

Test: `tests/test_coproducts.py::test_cointeraction_singular` (xfail, strict).
Everything else is green (`δ⁻` coassociativity, `Δ` comodule coassociativity, the
**gKPZ** cointeraction, twisted antipode, structures): 33 passed, 1 xfailed.

---

## 1. The bug

The cointeraction `(Id⊗Δ)δ = M¹³(δ⊗δ⁺)Δ` (tourist_guide.tex 5717, appendix proof
§C.3.3) holds for the gKPZ class (β₀=−1) but **fails at β₀=−3/2** on the smallest
*decorated* tree:

    τ = ∘ —I₀— ∘^{(0,1)}        (root ∘, an I₀ edge, child = noise with 𝔫-decoration (0,1))

`τ` is only divergent (`|τ|′ = 2β₀+3 = −2κ < 0`) at β₀ ≤ −3/2, which is why
β₀=−1 never exposes it.

The mismatch is a **decoration-placement** difference between the two composition
orders, for the same left forest `∘·∘^{(0,1)}`:

| path | mid (trunk red node) | branch leaf red node |
|---|---|---|
| LHS `(Id⊗Δ)δ` (extract→recenter) | `●^{(0,1)}(o=−3/2−κ)` | `●(o=−1/2−κ)` |
| RHS `M¹³(δ⊗δ⁺)Δ` (recenter→extract) | `●(o=−1/2−κ)` | `●^{(0,1)}(o=−3/2−κ)` |

i.e. the leftover `(0,1)` decoration + the `o`-value sit on **different red nodes**
on the two sides. The two red trees `●^{(0,1)}(o=−3/2−κ)` and `●(o=−1/2−κ)` have the
**same extended homogeneity** (`−1/2−κ`) but different **naive** homogeneity
(`1` vs `0`) and different `(𝔫, o)`, so they are distinct basis elements and the
term-sums do not cancel. It is **not** a truncation artifact (raising `_E_CAP`
does not help).

How each side gets there (`δ`/`δ⁺` use the binomial node-split `n_φ`; `Δ` uses the
boundary-edge Taylor recentering `e'` with `πe'` pushed onto the trunk):
- **RHS:** `Δ` recenters the I₀ edge (`e'=(0,1)` → edge `I_{(0,1)}`, `(0,1)` pushed
  onto the trunk so `A=∘^{(0,1)}`), then `δ⁺` extracts the child with split
  `n_φ=(0,0)`, leaving `[n−n_φ]=(0,1)` on the **branch** red node → `●^{(0,1)}(o=−3/2−κ)`.
- **LHS:** `δ` extracts root and child as **separate** subtrees joined by a
  *between-edge* (which gets **no** `e'` by the formula — `∂φ` is only φ→outside),
  with split `n_φ=(0,1)` so the child is extracted as `∘^{(0,1)}` and the **trunk**
  red node later receives `(0,1)` from `Δ`'s `πe'`.

Both `δ` and `Δ` individually match the tourist-guide formulas (tex 5613/5636) and
pass their **own** coassociativity even at β₀=−3/2. The bug is purely in their
**coupling**.

---

## 2. The authoritative resolution (BHZ 1610.08468 + Chevyrev survey 2206.14557)

This is the *documented* failure mode the extended decoration `o` exists to cure
(BHZ Remark 5.38; Chevyrev survey **Example 3.17** is essentially this tree).

- The reconciliation is **not combinatorial** (no extra term, no redistribution).
  It is a **degree/grading device**: `Δ⁺` is truncated by the degree
  `|·|₊ = |·|_𝔰 + Σ_x o(x)` (BHZ Def `homog`; survey eq. 3.16 calls this the
  *crucial* `|·|₊`-degree-preservation `|C(τ;A,nₐ,eₐ)|₊ = |τ|₊`). The contracted
  node's `o` is exactly `|extracted piece|_𝔰` (survey Remark 3.18 / eq. 3.14:
  `(o(A)+[n_A+πe_A]_A)(x₁) = |E(τ;A,n_A,e_A)|_𝔰`).
- The whole identity reduces to BHZ eq. `(e:crucoin)`:
  `Δ⁻_ex ∘ proj⁺_ex = (id ⊗ proj⁺_ex) Δ⁻_ex`, which holds *because* `Δ⁻_ex` preserves
  `|·|₊` on the right factor. The single coproduct `Δ_i` on coloured forests
  projects to both `Δ⁺_ex = (id⊗proj⁺_ex)Δ₂` and `Δ⁻_ex = (proj⁻_ex⊗id)Δ₁`
  (BHZ Cor `cor:domains`); the master identity is
  `M^{(13)(2)(4)}(Δ_i⊗Δ_i)Δ_j = (id⊗Δ_j)Δ_i` for `i<j` (BHZ Prop `prop:doublecoass`).
- **Footnote-3 landmine** (survey, p.32): BHZ's printed `T⁻_ex` (their 5.23) has a
  typo — the basis must **not** allow polynomial decorations at the root, single
  noise symbols `Ξ_ℓ`, or `1^p I_m τ` with `p≠0`. Stray polynomial decorations on
  contracted/extracted roots are a known source of exactly this kind of bug.

Sources (PDFs were fetched to /tmp during research; cite by `\label`):
BHZ arXiv:1610.08468 — `def:Deltabar`, `∂(A,F)` (part0), `CKop`/`e:defhatn`,
`homog`, Thm 5.37 (`e:propWanted1`, `e:def1324`), `e:crucoin`, `e:degreeDelta`,
Remarks `explanation`/`rem:explanation`/`rem:fails`. Chevyrev arXiv:2206.14557 —
Examples 3.17/3.24/3.28, eqs 3.13–3.17, Theorem 3.31, footnote 3.

---

## 3. The fix attempt that FAILED (and why)

Following the survey's prescription literally — *a contracted node is `●^{0,α}`,
carrying no polynomial decoration; absorb the leftover `[n−n_φ]` into `o`* (so
`●^{(0,1)}(o=−3/2−κ) ↦ ●(o=−1/2−κ)`, both becoming `●^{0,−1/2−κ}`) — does make the
two cointeraction sides agree on the failing tree. **But it breaks `δ⁻`
coassociativity**: identifying those red nodes makes two distinct coassociativity
terms collapse, so a coefficient that should be `1` becomes `2`
(`{((∘,),(),()):2} ≠ {…:1}`).

So the tension is real and is the crux:
- **`δ⁻` coassociativity** needs `●^{(0,1)}(o=−3/2−κ)` and `●(o=−1/2−κ)` **distinct**.
- **The cointeraction** needs them **reconciled**.

A blanket "absorb decoration into `o`" is therefore wrong. The correct structure
keeps them distinct *as basis elements* but makes the cointeraction sums match via
the `|·|₊`-truncation of `Δ⁺` — which my implementation already uses
(`_blue_positive` tests `extended_homogeneity = naive + Σo = |·|₊`), yet the
spurious branch `I_{(0,1)}(●^{(0,1)}(o=−3/2−κ))` has `|·|₊ = 1/2−κ > 0` and so is
*kept* by `p₊`. That means the discrepancy is **not** removed by the truncation as
I currently apply it — pointing at a deeper mismatch between my direct enumeration
and BHZ's `Δ_i`-via-projections construction.

---

## 4. Path to a correct fix

The robust route is to stop computing `δ`/`Δ`/`δ⁺` by separate ad-hoc enumerations
and instead follow the **proven-correct** construction:

1. **Single master coproduct.** Implement BHZ's `Δ_i` on coloured decorated forests
   (`def:Deltabar`), or equivalently the tourist-guide appendix's factorisation
   `D(𝔽τ) = (𝔻𝔽)(*Dτ)` (tex 7020–7072): an undecorated graph part `*D` plus the
   uniform decoration-operator coproducts `𝔻X` (eq. 7026) and `𝔻I` (eq. 7027). The
   appendix proves coassociativity *and* the cointeraction on this factorisation
   (the "either a or c is 0" argument, 7247, and the `o` change-of-variables
   `ζ↔φψ, η↔σ/φ`, 7277–7289). Re-deriving my decoration distribution to match
   `𝔻X`/`𝔻I` exactly should remove the discrepancy without the coassociativity
   regression.
2. **Corrected `T⁻_ex`/`C⁺` basis** (footnote-3): forbid the disallowed decorated
   roots, so the projections `p₊/p₋` quotient correctly.
3. Re-validate against the *current* green invariants (they must stay green) plus
   the singular cointeraction; then flip the xfail.

Estimated effort: a focused re-implementation of `coproducts.py`'s decoration
handling around the `(𝔻𝔽)(*D)` factorisation — moderate, but it touches the
load-bearing core, so it should be done test-first (keep all current invariants).

## 4b. Update — second attempt + the decisive-check result

Two further things were established:

- **The "absorb into o" lever is confirmed wrong** by experiment: making red nodes
  `●^{0,α}` (absorb the leftover `[n−n_φ]` into `o`) *does* fix the singular
  cointeraction but **breaks `δ⁻` coassociativity** (a coefficient goes 1→2). The
  appendix's proven construction keeps `●^{[n−n_φ],α}` **distinct**, so absorption
  is provably incorrect. A targeted variant (absorb only on `Δ`'s push onto red
  nodes) fixes the *trunk* discrepancy but not the *branch* one. Neither is right.

- **The reduced (`*D`, no-`o`, no-projection) cointeraction test is inconclusive.**
  The reduced `°D` maps are genuinely *infinite* sums (the boundary-edge `e'`
  Taylor sum is unbounded; they are triangular/bigraded — tex 5656 Remark). A naive
  degree-truncation makes both sides finite but introduces boundary mismatches that
  *grow* with the cutoff, so it cannot localise the bug. The only **finite** tests
  are the *projected* ones (`δ`,`δ⁻`,`Δ` with `p±`), which is why `δ⁻`
  coassociativity and `Δ` comodule pass cleanly but the projected cointeraction
  (the one three-way coupling) fails.

**Conclusion on effort.** A correct fix is **not** a localized patch. It requires
implementing the construction the way BHZ actually proves it: the single master
coproduct `Δ_i` on coloured decorated forests over **bigraded spaces** (so the
infinite `e'` sums are well-defined triangular maps), with the corrected `T⁻_ex`
basis (footnote-3), then deriving `Δ⁺/Δ⁻/δ⁺` as projections. That is a research-
grade, multi-session rebuild of the algebraic core — its payoff is a *structural-
consistency* property that feeds **no** library output (see §5). Recommended
either as a dedicated effort or to stand as a documented limitation.

## 5. Why this is not on the critical path (scope note)

The cointeraction is the **structural compatibility** check between `𝒯` and `𝒰`. It
does **not** feed any output the library computes: the renormalized-equation family
(Phases 1–2) uses no coproducts, and the symbolic BHZ character `k=h∘S'₋` uses only
`S'₋`. So this residual is a gap in *structural self-consistency validation* for
decorated-node-under-edge trees at β₀≤−3/2 — worth closing for a faithful RS
library, but not a blocker for the package's computations.
