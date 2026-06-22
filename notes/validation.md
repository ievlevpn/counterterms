# Validation report

The library has **no reference implementation** — *the paper is the oracle*. This
records what is checked against the literature, what is cross-checked independently,
and the bug the pass found and fixed. Tests live in `tests/test_validation.py`
(+ `test_goldens.py`).

## 1. Paper-validated (exact)

| Equation | What | Source | Result |
|---|---|---|---|
| gKPZ, β₀=−1−κ | the **5 counterterms** (|τ|, S(τ), F(τ*)) | tex 6004–6012 | exact ✓ |
| gKPZ, β₀=−3/2−κ | the **full SC-tree table**, 43 trees over 6 homogeneity rows | tex 6024–6163 | **exact ✓** |

The gKPZ@−3/2 row counts — `−3/2−κ:1, −1−2κ:2, −1/2−3κ:6, −1/2−κ:2, −2κ:9, −4κ:23` —
match the paper tree-for-tree (count per row).

## 2. The bug this pass found (and fixed)

Matching the gKPZ@−3/2 table exposed an **undercount**: the engine produced 42 trees,
the paper 43 — one missing at `−2κ`, the tree `●^{(0,1)}·𝓘ₓ[Ξ]²` (a `●` with an
`Xₓ` node-decoration and two gradient-noise children).

**Cause.** Tree generation skipped any root whose *bare* homogeneity already reached the
budget (`if base_h.std >= bound: continue`). But a node like `●^{(0,1)}` (std=1=bound)
is pulled **below** the bound by capped negative-contribution children: a gradient edge
`I_{(0,1)}` over a noise contributes `(m−|p|_𝔰)+β₀ = 1+β₀`, which is `< 0` once `β₀<−1`.
The skip wrongly assumed children only *raise* homogeneity. It bit only at `β₀<−1`
(at β₀=−1 those children contribute exactly 0), so the β₀=−1 golden never caught it.

**Fix** (`equation/generate.py`): drop the premature skip — the DFS already rejects the
bare node and terminates via the budget `break` (negative-contribution children are
derivative slots with a finite cap). Corrected counts: gKPZ@−3/2 42→**43**, KPZ 10→**11**,
Φ⁴₃ (via DPD) 9→**13**; PAM (g=0) and all β₀=−1 cases unchanged.

## 3. Independent cross-checks

- **Symmetry factor.** `S(τ)` recomputed by a brute child-permutation stabiliser count
  (`Π m_j!` enumerated, not the closed formula) × node-decoration factorials, across the
  whole corpus — matches `symmetry_factor()` everywhere.
- **Subcriticality / homogeneity sanity.** Every benchmark equation is subcritical and
  its counterterm homogeneities are the expected `kβ₀ + j` ladder.

## 4. Benchmark (sanity-checked predictions)

Post-fix counterterm counts (paper-validated where noted; others are subcritical,
`S=Aut`-consistent predictions — the paper gives no explicit table for them):

| Equation | β₀ | counterterms | note |
|---|---|---|---|
| gKPZ d1 | −1−κ | 5 | paper-exact (renormalised eqn) |
| gKPZ d1 | −3/2−κ | 43 | paper-exact (SC table) |
| KPZ d1 | −3/2−κ | 11 | homogeneity rows ⊂ paper table |
| gPAM d2 | −1−κ | 4 | |
| PAM d2 | −1−κ | 4 | |
| PAM d3 | −3/2−κ | 17 | |
| Φ⁴₂ (DPD) | −2 → lift | 3 | lifted remainder eqn |
| Φ⁴₃ (DPD) | −5/2 → lift | 13 | lifted remainder eqn |

## 5. Honest scope of the validation

- Exact literature validation is anchored on gKPZ (the paper's worked example); other
  equations are sanity-checked predictions, not cross-checked against published tables.
- This validates the **symbolic structure** (which trees, homogeneities, S(τ), F(τ*)) —
  not the **numeric** renormalisation constants (Track B2, unbuilt) nor convergence.
