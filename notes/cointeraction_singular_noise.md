# The cointeraction residual for singular noise — research + diagnosis + FIX

**Status: FIXED.** The cointeraction `(Id⊗Δ)δ = M¹³(δ⊗δ⁺)Δ` now holds at β₀=−3/2
(0 failures, ≤4-node trees) with no regression (34 passed). The root cause and the
two-line fix are in **§8** below; §1–§7 are the (long) diagnostic journey that led
there, kept for the record.

**The fix, in one sentence:** the `e'` Taylor recentering must be applied on the
*full* extraction boundary `∂(A,F)` = {edges with parent in the extracted forest,
not internal to a component} — i.e. **between-component edges too**, not only
`φ→outside` edges. (Plus: edges below a pre-existing red node stay internal.)

Test: `tests/test_coproducts.py::test_cointeraction_singular` (now a real, passing
test). All algebraic invariants green.

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

## 6. Rebuild-branch findings (localization) — `cointeraction-bigraded-rebuild`

Systematic bisection of the failure (all on the gKPZ equation, β₀=−3/2):

- **Multi-index `o` hypothesis: REJECTED.** Research (BHZ vs Chevyrev survey) shows
  the scalar `o ∈ ℤ[β₀]` is *lossless for the cointeraction*: Chevyrev proves the
  identical cointeraction with the scalar `o` = homogeneity of the extracted piece,
  and our contracted-node `o` already equals that. The multi-index `o` only matters
  for type-level reconstruction (Π^ζ), which the MVP never does.
- **Pure-graph cointeraction (no `n`/`e`/`o`, finite): HOLDS.** So the subforest /
  subtree / contraction combinatorics are correct.
- **`δ⁻` coassociativity, `Δ` comodule (pure and red-containing): HOLD.**
- **`δ⁺` (`D̄⁻`, root-disjoint) comodule coassociativity `(δ⁻⊗Id)δ⁺ =
  (Id⊗δ⁺)δ⁺`: FAILS at β₀=−3/2** (1 of 10 T⁺ factors). This is the bug. It is the
  only map exercised solely inside the cointeraction. The failing case is the
  re-application of `δ⁺` to a right factor carrying **nested red nodes under the
  blue root** (e.g. `δ⁺(blue—I₀—∘—I₀—∘)` re-applied): a `(∘^{(0,1)}, ∘, …)` term
  arises on the `(Id⊗δ⁺)δ⁺` side with no partner on `(δ⁻⊗Id)δ⁺`.

**So the remaining fix is narrow:** correct `δ⁺`'s root-disjoint re-extraction so
its comodule coassociativity holds (likely the interaction of `root_disjoint`
with the `red_ids ⊆ φ` requirement and the `p₋` red→𝟙₋ rule on re-extraction of
nested-red/blue trees). Once `δ⁺` comodule is green, the cointeraction should
follow (it is the last un-validated piece; graph + `δ⁻` + `Δ` are all verified).

## 7. Rebuild progress (branch) and the irreducible remaining crux

**Landed (commit 697612f, no regression — 33 passed):** force edges below a
pre-existing red node to be *internal*. A red node is a contracted placeholder;
on re-extraction its extracted children must merge into it, not split off as
sibling red nodes via a between-edge. This removed a whole class of spurious `δ⁺`
terms (the δ⁺-comodule-coassociativity failure on nested-red trees) and cut the
β₀=−3/2 cointeraction failures to **2 trees**, both the original
`∘—I₀—∘^{(0,1)}` recentering case.

**The irreducible crux (still open).** For `∘—I₀—∘^{(0,1)}` the `e'=(0,1)`
recentering produces, in the middle (T) leg:

| path | middle red node |
|---|---|
| LHS `(Id⊗Δ)δ` (extract→recenter) | `●^{(0,1)}(o=−3/2−κ)` — `Δ` pushes `πe''` onto the red trunk as a **node decoration** |
| RHS `M¹³(δ⊗δ⁺)Δ` (recenter→extract) | `●(o=−1/2−κ)` — `Δ` recenters the root first, so `δ` extracts it with `(0,1)` folded into **`o`** |

Both have `|·|₊ = −1/2−κ` but are **distinct basis elements** `(𝔫=(0,1),o=−3/2−κ)`
vs `(𝔫=0,o=−1/2−κ)`, sitting in the **same tensor leg** — so the term-sums do not
match. Identifying them (absorb 𝔫 into o on red nodes) fixes the cointeraction but
breaks `δ⁻` coassociativity (§3, §4b). Confirmed (research §5/§6): the scalar `o`
is *lossless*, so the resolution is **not** the multi-index `o`; it is that these
two genuinely-distinct elements must be **paired with different leg-1/leg-3
partners that balance in the full sum** — i.e. there is a missing/extra term in
one of the two composition orders, in how the recentering `X^k` of `Δ` is shared
with the extraction `e'` of `δ` on the **same physical edge** (root→child, which
both maps touch in the composition). Pinning that exact term-level balance is the
remaining work.

## 8. THE FIX (root cause + resolution)

Two combinatorial corrections to the extraction coproduct `δ` (`delta_minus`),
both about how a *subforest* `A` (the "between-edge" subtlety) interacts with the
`e'` Taylor recentering. The decisive references are BHZ arXiv:1610.08468:
`def:Deltabar` (the master `Δ_i`) and the boundary set `∂(A,F)`.

**(a) The boundary `∂(A,F)` includes between-edges.** BHZ define
`∂(A,F) = { e ∈ E_F∖E_A : e₊ ∈ N_A }` — every edge whose *parent* lies in the
extracted forest `A` and which is *not internal to a component* of `A`. That set
is **both** the `φ→outside` edges **and** the between-component edges (parent and
child both in `A` but in different subtrees). The `ε_A^F` (our `e'`) recentering
runs over all of `∂(A,F)`: `πe'` onto the parent, `e+e'` on the contracted edge.

My implementation had recentered only the `φ→outside` edges (the notes in §1
literally asserted "between-edges get no e'" — that was the bug). With between-
edges also recentered, `δ` can recenter a node *whose child is also extracted*,
which is exactly the term the recenter-then-extract order (`M¹³(δ⊗δ⁺)Δ`) produces
via `Δ` and the extract-then-recenter order (`(Id⊗Δ)δ`) was missing.

**(b) Edges below a pre-existing red node stay internal.** A red node is a
contracted placeholder; on *re-extraction* (the second coproduct in the
cointeraction / the δ⁺ comodule) its extracted children must merge into it, never
split off as a sibling red node via a between-edge. Without this, `δ⁺` over-
produced nested-red terms and failed its own comodule coassociativity.

Both live in `delta_minus`'s subforest enumeration; `_assemble` already applied
`πe'`/`e+e'` to *any* cross-component edge, so (a) needed only to add the between-
edges to the recentering list. Net diff: a few lines. Commits on this branch:
`697612f` (b), and the between-edge recentering (a).

**Why the earlier "absorb into o" attempts were wrong** (§3, §4b): they tried to
*reconcile* the two genuinely-distinct basis elements `●^{(0,1)}(o=−3/2−κ)` and
`●(o=−1/2−κ)` by identifying them. But the cointeraction does not identify them —
it produces the *same* term on both sides once the recentering combinatorics are
correct. The discrepancy was never in the `o`-decoration (the scalar `o` is
lossless, confirmed by research §5/§6); it was a missing set of `δ` terms.