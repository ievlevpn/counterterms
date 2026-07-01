# Mathematical Correctness Audit — `counterterms/`

**Date:** 2026-06-22
**Scope:** *Mathematical correctness only* — every formula, convention, sign, factor, and
combinatorial rule in `counterterms/`, checked against the source paper
(Bailleul & Hoshino, *A tourist's guide to regularity structures…*, arXiv:2006.03524 v3,
`references/tourist_guide.tex`, the sole oracle). This is the math counterpart to the
code-quality `audit.md` (which explicitly did *not* check math).
**Method:** 14 independent verification agents (one per math concern), each reading the
relevant code *and* the relevant `.tex` sections and **running the engine to compare its
output to hand-derivations from the paper**; every flagged finding was independently
re-derived by adversarial skeptics; a completeness critic swept for blind spots. The
auditor (me) additionally hand-checked the three most load-bearing pieces (homogeneity
ring, Υ-map, assembly) against tex 6004–6012, and ran a direct before/after experiment on
the one finding whose second skeptic died mid-run. Read-only except for one temporary,
reverted experiment (§F2). `uv run pytest` → **169 passed**.

> Working notes — a point-in-time snapshot (now tracked in-repo; the F1/F2 fixes landed in
> commit d3adcaf, and the two generation gaps found later were fixed with regression tests in
> tests/test_rule.py).

---

## Verdict

**The mathematics of the delivered pipeline is correct and faithful to the paper.** Every
load-bearing convention and formula was positively verified, in most cases by executing the
engine and matching the paper's own worked numbers:

- The **gKPZ renormalized equation** (tex 6004–6012) reproduces **term-for-term**, including
  the asymmetric factor-2 (see the stringent cross-check below).
- The **KPZ homogeneity table** (tex 6028–6063) reproduces **entry-for-entry** — all 43
  trees, per-row counts (1, 2, 6, 2, 23, 9), at homogeneities β₀, 2β₀+2, 3β₀+4, β₀+1,
  4β₀+6, 2β₀+3.
- The **cointeraction identity holds including the singular β₀ = −3/2 case** (verified
  uncapped over all KPZ divergent trees, not a `max_nodes` artifact).

**No error affects any counterterm the engine outputs.** Two genuine defects exist; **both
live in layers that do not feed the delivered renormalized equation**:

| # | Severity | What | Affects output? | Status |
|---|----------|------|-----------------|--------|
| **F1** | minor (real) | Rule generation over-generates trees for `d ≥ 2` multi-direction gradient nonlinearities | No — masked by the Υ-map; no corpus test | **Fixed 2026-06-22** |
| **F2** | minor (latent) | δ⁺ comodule law (compatibility condition **(b)**) fails in the deferred Phase-3 coproduct layer | No — that map is used only in a test | **Fixed 2026-06-22** |

The rest are info-level latent/robustness notes (§Info), plus correctly-deferred analytic
pieces (§Deferred).

> **Both fixed 2026-06-22** (suite now 176 passed; the canonical/BHZ character is byte-identical
> to before the fixes). **F1**: `dsl.py`/`signature.py`/`generate.py` now cap the *total* gradient
> edges per node by the D2 total-degree (the paper's rule, tex 5337-5340) — `grad_budget` on
> `Signature`, enforced in `_emit`; regression test `test_rule.py::test_d2_total_gradient_bound_in_dim2`.
> **F2**: dropped the redundant force-internal-under-red rule (`coproducts.py:198`); regression test
> `test_coproducts.py::test_delta_plus_comodule` (condition (b) over the T⁺ basis incl. β₀=−3/2).
> Each test was confirmed to fail on the unfixed code and pass on the fixed code. Not committed.

---

## The stringent cross-check (why "it runs" ≠ "it's right")

The gKPZ family at ζ ∈ C^{−1−κ} (tex 6004–6012) is a *non-accidental* test of the S(τ) and
Υ-map machinery simultaneously, because the factor 2 behaves differently in two of its terms:

- **τ₃ = ●—∂ₓ—∘** : `S = 1`, and `Υ(τ₃*) = ∂_{∂ₓu}[g(u)(∂ₓu)²]·f = 2·f·g·∂ₓu`.
  Coefficient `k/S · Υ = 2k·f·g·∂ₓu` → the paper's explicit **`2k(τ₃)`**.
- **τ₅ = ● with two ∂ₓ—∘ branches** : `S = 2` (two identical branches ⇒ `m! = 2`), and
  `Υ(τ₅*) = ∂²_{∂ₓu}[g(u)(∂ₓu)²]·f² = 2·f²·g`. Coefficient `k/2 · 2f²g = k·f²·g` → the
  paper's plain **`k(τ₅)`**.

The factor 2 *survives* in τ₃ and is *cancelled* in τ₅. Both `S(τ)` (`trees/tree.py:66-73`)
and the Υ-map (`renorm/nonlinearity.py:47-55`) must be individually correct for this to come
out as the paper writes it. The code reproduces both exactly. I confirmed this by hand and by
running `tests/test_render.py`.

---

## Verified correct (the bulk of the math)

Each item below was checked against the paper *and* exercised in code. Citations are tex line
numbers.

**Homogeneity arithmetic** (`core/{homogeneity,signature}.py`, `trees/tree.py`)
- `|Ξ| = β₀` directly — **no off-by-2**; the "regularity = α−2" phrasing at tex 705 is an
  informal intro, the operative definition `|∘| = β₀` is tex 5276–5283. The single most
  dangerous convention, and it is correct.
- `|I_p τ| = |τ| + order − |p|_𝔰` (tex 2213/2224), with the operator order read per-component
  (`signature.py:45-47`), correctly generalizing the heat operator's "2".
- `|X^n| = |n|_𝔰` with parabolic scaling 𝔰 = (2,1,…) (tex 1153); time counts double.
- Ordered ring ℚ ⊕ ℚ·κ, lexicographic (std then κ); `is_negative` correctly **includes** the
  critical −kκ trees (std 0, kap < 0) and **excludes** exact 0 (tex 6066). Floats would
  mis-sign these.

**Decorated trees & canonical isomorphism** (`trees/tree.py`)
- Data model (node = 𝔱×ℕ^{d+1}, edge = 𝔗_e×ℕ^{d+1}) matches tex 3963–3970; **component identity
  rides on the edge type**, not the node (tex 3826-3827). Canonical key is a genuine recursive
  total order (verified by a 120-permutation collapse test); load-bearing for `S`, dedup, dict
  keys.

**Symmetry factor** `S(τ) = n!·Π_j S(σ_j)^{m_j}·m_j!` (`trees/tree.py:66-73`) — exactly tex 3982,
including the multi-index `n!` and grouping identical child-edges by `(component, p, σ)`.

**Υ-map** `F(τ*) = (Π_i F(τ_i*))·(D^n Π_i ∂_{p_i}) F(b*)` (`renorm/nonlinearity.py:47-55`) — tex
4337. The two subtle points both hold: `∂_{p_i}` is applied **before** `D^n` (they don't
commute — verified the orders genuinely differ), and `D_i = Σ_k u_{k+e_i} ∂_k` (tex 4316) is
correctly iterated (Faà di Bruno). Base cases `F(∘_j*) = f_j`, `F(●*) = g`. The recursion
property `Υ(↑_i τ*) = D_i Υ(τ*)` (tex 4597) holds on all test trees.

**Assembly / ThmRenormPDEs** (`renorm/equation.py:36-40`, `api.py`) — coefficient is `k(τ)/S(τ)`,
never `k(τ)` alone (tex 4915); one free constant per tree, **shared across components**;
sum over exactly 𝓑_{<0} **including bare primitives ∘ⁿ** (the `k(∘)`, `k(∘1)` terms), **excluding**
bare ● (|●| = 0). Trees with `F(τ*) = 0` are dropped per-component (a genuine no-op — verified
the dropped KPZ trees have vanishing Υ).

**Coproducts Δ⁻, Δ⁺** (`trees/coproducts.py`, `core/hopf.py`) — extraction/contraction,
subforest combinatorics, the 𝔬-decoration of the contracted node (`[𝔬]_φ + |component|'`,
matching the worked example tex 6170-6176), the 𝔢′ Taylor recentering, the p₋/p₊ projections
(naive homogeneity for U/U⁻, **extended** for T/T⁺ — a non-trivial split, tex 5750), all
counits, stability, and **all three coassociativities** (δ⁻, Δ⁺, and Δ-as-T⁺-comodule) pass over
the whole corpus including β₀ = −3/2. The connected-graded antipode satisfies S⋆id = η∘ε.

**Cointeraction (compatibility condition (c))** `(Id⊗Δ)δ = M¹³(δ⊗δ⁺)Δ` (tex 5717) — the test is
a faithful, leg-by-leg-correct transcription with the right maps; **holds including β₀ = −3/2**,
verified uncapped on all 11 KPZ divergent trees and on the historically-failing synthetic tree
∘—I₀—∘^{(0,1)}. The CLAUDE.md / `structures.py` "cointeraction holds incl. β₀=−3/2" claim is
substantiated.

**Twisted antipode / BHZ character / G⁻** (`coproducts.py:423`, `structures.py`) — S′₋ implements
the Dyson–Salam characterizing relation (tex 5034) exactly (sign, the `τ⊗●^{0,|τ|}` subtraction
using the **extended** homogeneity, the M₋(S′₋⊗Id) structure); verified via the identity
`M₋(S′₋⊗Id)δτ = S′₋(τ)·(●^{0,|τ|}−1)` on every KPZ/gKPZ tree, and via the BPHZ mean-zero property
𝔼[^{kζ}Π^ζτ](0) = 0. S′₋ is correctly kept **distinct from a Hopf antipode** (tex 5047). The
renormalization group G⁻ (convolution, unit ε⁻, inverse f∘S) satisfies all group axioms,
including on multi-component forests.

**Subcriticality** (`structures.py`, `rule.py`) — `β₀ > −order` faithfully implements tex
5485/2024; the heat window β₀ ∈ (−2,0) reproduced; Φ⁴₂/Φ⁴₃ correctly rejected and routed to the
lift; `min` over multiple noises binds (tex 1160).

**da Prato–Debussche lift** (`equation/daprato.py`) — matches the paper's own Φ⁴₃ example (tex
2026-2034) line-for-line: `u = X + v`, remainder `−v³−3v²X−3vX²−X³`, Wick powers X^k ∈ C^{(−k/2)⁻}
via the Schauder gain α_X = β₀ + order, reduced β₀ = (−3/2)⁻ > −2. Substitution algebra verified
symbolically for cubic and mixed polynomials. Subcriticality of the reduced problem is enforced
downstream (an insufficient lift is correctly rejected).

**BPHZ scheme / Wick** (`renorm/scheme.py`) — canonical model Π^ζ (tex 5042-5047), k^ζ = h^ζ∘S′₋,
multiplicative factorization over forests, centered-Gaussian parity (odd ⇒ 0), Isserlis
perfect-matching combinatorics (each matching once, no spurious prefactor), heat-kernel
integrand (one ∂^p K per edge), within-noise-type-only pairing, X^n-decoration refusal — all
correct. Hand-derived integrands for ∘—I₀—∘ and a 4-noise tree reproduce.

**Golden recomputes** — gKPZ (5 counterterms) and KPZ/PAM/gPAM/Φ⁴ homogeneity tables all match
the paper *and* the golden test assertions; the goldens encode the **correct** values (a golden
encoding a wrong value would be the most dangerous bug — none found).

---

## Findings that warrant action

### F1 — Rule generation over-generates for `d ≥ 2` multi-direction gradient nonlinearities  *(minor, real — both skeptics confirmed)*

**Location:** `counterterms/equation/dsl.py:243-245`.
**Paper:** tex 5337-5340 (a node carries at most two gradient edges I_{e_i}, I_{e_j}); tex 4555
(a node with ≥ 3 leaving edges of |p|_𝔰 = 1 is *not* strongly conforming).

The rule construction caps each derivative-jet edge **independently** (`caps[(comp,p)] =
degree(g, that single jet)`). The *scope check* (`dsl.py:227`) correctly uses
`Poly(...).total_degree() ≤ 2`, but the rule caps do not preserve that total bound. For
`g(u)(∂₁u + ∂₂u)²` in `d = 2`, each of the two gradient jets gets degree 2, so a ● node admits
up to `2 + 2 = 4` gradient children — degree-4-in-∂u, which violates Assumption D2 and is absent
from the paper's rule.

**Executed evidence (both skeptics, independently):** this equation generates 3 spurious trees
(e.g. `●[I₀₀₁, I₀₀₁, I₀₁₀]` at −3κ and `●[I₀₀₁, I₀₀₁, I₀₁₀, I₀₁₀]` at −4κ) — the correct count is
8, not 11. **Invisible in `d = 1`** (one gradient direction ⇒ per-jet cap = total cap), so all
gKPZ/KPZ goldens are unaffected. **Masked in the final PDE**: the Υ-map differentiates a degree-2
`g`, so a 3rd/4th derivative factor gives 0, and `api.py:28` drops `F(τ*) = 0` trees — the
renormalized equation still has the correct counterterms, 0 spurious.

**Impact:** the *raw negative-tree basis* (`generate_trees` / `generate_counterterms`, exposed as
`RenormalizedEquation.all_trees` and consumed by the Phase-3 coproduct/regularity-structure
layer, `structures.py:88,156,228`) is wrong for any `d ≥ 2` direction-mixing gradient
nonlinearity. No corpus test exercises such an equation, so the suite does not catch it. **The
canonical / BHZ character** is computed over this same `generate_counterterms` divergent set, so for
such a `d ≥ 2` equation it would gain spurious constant entries — one per Υ-zero spurious tree
(the character map does not see Υ). **No in-scope/corpus equation is affected** (all are `d = 1`
gradient or `d ≥ 2` without a gradient nonlinearity), and every *genuine* counterterm's character
is still correct; only the constant *list* for that out-of-scope class would be polluted.

**Fix:** cap the **total** number of gradient (`p ≠ 0`) children of a node by
`Poly(g, *gjets).total_degree()` (≤ 2 under D2), enforced in `generate._emit`'s DFS, rather than
per-edge. Still terminates. Add a `d = 2` gradient-nonlinearity test asserting no surviving tree
has > 2 gradient edges.

### F2 — δ⁺ comodule law (compatibility condition (b)) fails; a verified one-line fix exists  *(minor, latent — off the output path)*

**Location:** `counterterms/trees/coproducts.py:198-199` (the "force-internal-under-red" rule).
**Paper:** tex 3451-3454 condition **(b)**: `(Id⊗δ⁺)δ⁺ = (δ⁻⊗Id)δ⁺`; the unprojected source is
tex 5711 `(D⁻⊗Id)D̄⁻ = (Id⊗D̄⁻)D̄⁻`. This is the comodule axiom a *compatible* renormalization
structure must satisfy; it is **distinct from the cointeraction (c)** and is **not exercised by
any test**.

Here δ⁺ = `delta_minus(root_disjoint=True)` — the map used (only) on the RHS of the cointeraction
test. Two independent agents found condition (b) **fails** (66/195 mismatches at β₀ = −3/2; 99/119
on the gKPZ T⁺ basis), localized to the `forced = [edges below a red node]` rule at line 198,
which `notes/cointeraction_singular_noise.md` §7-8 says was added to make the cointeraction pass.

Because one of those findings lost its second skeptic to an API error, **I verified it directly**
(temporary patch, reverted):

```
forced = []            # was: [e for e in within if nodes[e[0]].color == "red"]
choice = list(within)  # was: [e for e in within if nodes[e[0]].color != "red"]
```

| state | condition (b), gKPZ T⁺ | full test suite |
|-------|------------------------|-----------------|
| **original** (rule present) | **4 / 17 fail** | 169 passed |
| **patched** (rule removed)  | **0 / 17 fail** | **169 passed** (incl. `test_cointeraction`) |

So removing the rule **fixes condition (b) at no cost to the tested cointeraction (c) or any
other test.** The note's claim that the rule is needed for the cointeraction is **stale** — it
predates the between-edge-recentering fix (the `recenter = boundary + between` at line 217, see
§Info), which alone now carries the cointeraction. The agent that reconstructed (b) reached the
same conclusion (its patch drove β₀=−3/2 failures 66→0 while keeping cointeraction green on all 43
divergent trees).

**Impact: none on any delivered quantity today.** `delta_minus(root_disjoint=True)` is used only
in `tests/test_coproducts.py:179`; the renormalized-equation family uses no coproducts at all
(free symbolic constants), and the BHZ character uses S′₋ over `delta_minus` (not the root-disjoint
variant). The defect would only bite a future use of the T⁺ comodule structure.

**The canonical / BHZ character is *not* affected** (verified directly). Although the
`forced`-under-red rule lives inside `delta_minus`, which `twisted_antipode` (= S′₋) calls,
`canonical_character` / `bhz_character` / S′₋ are **byte-identical with and without the fix** across
all six corpus equations — identical SHA-256 fingerprint of every S′₋ forest-sum, including KPZ's
1435-term S′₋ at β₀ = −3/2, and the same `canonical-char == 0` counts (gkpz 3/5, kpz 4/11, …). The
reason: S′₋ recursively applies δ only to the **black** extracted subforests, where the rule is
inert (`forced = []` regardless, since `red_ids` is empty there); the rule fires only when
extracting from *below an already-red node*, which the twisted-antipode recursion never does. So
F2 is fully contained to the root-disjoint δ⁺ comodule path used only by the cointeraction test.

**Recommendation:** adopt the one-line change **and** add a permanent condition-(b) regression test
(mirroring `test_delta_comodule_coassociative` but for δ⁺) before relying on it — the convergent
evidence (my gKPZ experiment + the agent's β₀=−3/2 reconstruction + full suite green) is strong but
a pinned test is what should guard it. Alternatively, if the Phase-3 T⁺ comodule structure is not
going to be used soon, mark `root_disjoint=True` as an unverified socket so it isn't mistaken for a
validated δ⁺. Update the note headline so "cointeraction holds" is not read as "all of conditions
(a)–(d) hold".

---

## Info-level / latent (no wrong output; noted for completeness)

These are *not* bugs in any computed result. Listed so nothing is silently carried forward.

- **`F(red*) = 0` not enforced in `elem_diff`** (`nonlinearity.py:49`): a red node's `node_type` is
  `"bullet"`, so a red node would return `g` not `0`. **Unreachable** — the Phase-1/2 generator
  emits only black trees. Add `if t.color == "red": return 0` if `elem_diff` is ever reused on the
  F^{(k)} path. Per tex 4329.
- **`_p_minus` maps a bare *black* ● to 𝟙₋** (`coproducts.py:363`): tex 5760 says only the *red*
  ●^{0,α} → 𝟙₋; a black ● (homogeneity 0) → 0. **Unreachable** — no δ produces a bare-black-bullet
  right factor (verified over the whole corpus); dropping the `or node_type=="bullet"` clause keeps
  all 169 tests green. Likely a conflation of the capital-P₋ (tex 4988) and lowercase-p₋ (tex 5760)
  maps. Tighten to `color == "red"` only.
- **Subcriticality compares standard parts only** (`rule.py:29,34-35`): `β₀ = −order + κ` (e.g.
  −2+κ) is `> −order` in the ordered ring (subcritical), but the std-only test rejects it.
  **Unreachable by physical noise** — space/spacetime white noise always has a *negative* κ-part,
  so std is the binding constraint and the check agrees with the paper there; Φ⁴₂ at −2−κ is
  genuinely supercritical and correctly rejected. One-sided (never wrongly *accepts*). For full
  fidelity to convention #2, compare `(edge_gain + β₀_min)` in the ring.
- **Multi-noise covariance label** (`scheme.py:146,154`): independent noises of *different* laws
  both render with the single symbol `C`, conflating C_ξ ≠ C_η in the (deferred) Wick integrand.
  The pairing *combinatorics* are correct (within-type only). Correct for the paper's only
  multinoise example (identically-distributed noises), and the API cannot currently express distinct
  laws, so no representable input is wrong today. Tag `C` by the matched pair's noise type for
  forward-compat.
- **`canonical_character` keeps provably-zero/one h-symbols**: only the parity reduction (odd-noise
  ⇒ 0) is applied; the further `h(●^{0,α}) = 1`, α-independence (tex 4003), and root-X^n ⇒ 0 (tex
  5083) collapses are deferred (Track B2). Symbolically correct; a consumer must not treat the
  surviving symbols as the final reduced constants. Documented in `elementary_expectations.md`.
- **Between-edge 𝔢′ recentering vs the paper's literal ∂φ** (`coproducts.py:217`): the code recenters
  between-component edges, which the tourist guide's *literal* ∂φ (tex 5566) excludes — but this is
  **correct and deliberate**, not a bug. **Both skeptics refuted the "major" framing** and the
  original auditor recommended keeping it: the paper's literal ∂φ *breaks its own cointeraction* at
  β₀ ≤ −3/2, whereas the code follows BHZ (arXiv:1610.08468)'s broader ∂(A,F), which makes the
  cointeraction hold. Documented in `cointeraction_singular_noise.md` §8. Worth a code comment + a
  regression test pinning δ on the affected trees so the convention can't silently drift back.
- **`kpz()` fixture is *classical* KPZ** (`ξ + (∂u)²`), a strictly smaller rule than the paper's
  *generalized* KPZ table at tex 6028 (`f(u)ξ + g(u)(∂u)²`). The golden checks only the homogeneity
  *set* (correct), but the comment cites the gKPZ table — worth clarifying. The true gKPZ-at-β₀=−3/2
  tree-by-tree match (43 trees) was verified separately by the engine.
- **da Prato κ-deficit `k·κ` vs the paper's single `⁻`** (`daprato.py:64-66`): immaterial — the
  counterterm tree set is identical under either convention (verified); only printed κ-coefficients
  differ. Consistent with how the paper's own gKPZ table accumulates κ per noise leaf.

---

## Correctly-deferred (analytic/probabilistic, out of symbolic scope — not bugs)

- **`G⁻_ad` admissible subgroup** (`structures.py:212`, `NotImplementedError`): k^ζ ∈ G⁻_ad needs
  kernel vanishing moments + the Π-map (tex 5083-5091) — analytic, not symbolic.
- **Canonical BPHZ *values***: evaluating the ε-regularized divergent integrals is Track B2; the
  scheme correctly leaves the kernel symbolic and unevaluated.

---

## Audit caveat

One skeptic (the second adversarial check on the cointeraction concern, lens "code execution")
**died on an API error** mid-run, leaving the δ⁺ condition-(b) finding with a single adversarial
vote. I treated that as a reason to verify the finding's core claim **myself** (§F2,
experiment), so the conclusion does not rest on a single unverified agent. Every other finding
received its full adversarial pass.

**Bottom line: the delivered SPDE → renormalized-counterterm-family math is correct. The two real
defects (F1, F2) are confined to the raw-basis / deferred-Phase-3 layers and change no output the
engine produces; F2 has a verified one-line fix.**
