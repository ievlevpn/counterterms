# Does `regstruct` do justice to the *Tourist's Guide*? — coverage report

**Date:** 2026-06-22
**Question asked:** does the codebase cover everything in the paper (arXiv:2006.03524 v3),
or what's missing?
**Method:** mapped every `\section`/`\subsection` of `references/tourist_guide.tex` against
the modules in `regstruct/`. Read-only; nothing implemented. Cross-checked against
`notes/{initial_plan,architecture}.md`, `ROADMAP.md`, `notes/validation.md`.

This is a **paper-coverage** map, complementary to `notes/audit.md` (which is code-quality only).

---

## One-paragraph verdict

The paper has two halves: an **algebraic/combinatorial** half (rules → decorated trees →
homogeneities → coproducts → twisted antipode → BHZ character → the BCCH renormalized-equation
formula) and an **analytic/probabilistic** half (models, modelled distributions, the
reconstruction theorem, Schauder estimates, the fixed-point solver, and the *numeric* values of
the BPHZ constants via Wick integrals). **`regstruct` implements the algebraic half essentially
in full, and the analytic half essentially not at all — by deliberate, documented design.**
The project's stated goal — *SPDE in → family of renormalized equations with free constants out*
(`ThmRenormPDEs`, the §5 / §6 material) — is met and golden-tested. Everything missing is either
(a) explicitly out of scope because it needs probabilistic/analytic input a symbolic package
can't supply, or (b) a named, socketed Phase-4 seam. There are **no silent gaps** in the part the
project claims; the omissions are the analysis, on purpose.

---

## Section-by-section map

Legend: ✅ implemented & tested · 🟡 partial / combinatorics-only · 🔌 socket (named stub,
raises) · ⬜ not built (out of scope) · n/a expository.

| Paper § (tex) | Topic | Status | Where / why |
|---|---|---|---|
| §1 Introduction (666) | motivation | n/a | — |
| §2.1 Algebra of local expansions (914) | pre-Lie / Hopf intuition | n/a | expository; the real algebra lands in §5/§7 below |
| §2.2 Regularity structures (1030) | the graded `(T, G)` object | 🟡 | `structures.py RegularityStructure (T,T⁺)` builds the **graded model basis** `T` and `T⁺`; the **group `G` action on a model** is analytic (⬜) |
| §2.3 **Models & modelled distributions** (1241) | `(Π, Γ)`, `D^γ` | ⬜ | **not built** — no `Π`, `Γ`, `‖·‖_γ`. The analytic core; out of scope (needs realizations of the noise) |
| §2.4 Products & derivatives (1879) | `⋆` on `D^γ`, `∂` | ⬜ | analytic; the *tree-level* derivative/product lives in the Υ-map & coproducts instead |
| §3.1 Operators on ℝ×ℝᵈ (2043) | kernels `K`, scaling | 🟡 | scaling/Schauder **exponent** bookkeeping is in `core/{signature,homogeneity}` (`|I_pτ|=|τ|+(m−|p|_𝔰)`); the **kernel `K` itself** ⬜ |
| §3.2 Abstract integration `I` (2179) | `I_p`, abstract integration | ✅ | edges with `(edge_type, p)` decoration in `trees/tree.py`; homogeneity recursion exact (golden-tested) |
| §3.3 **Admissible models** (2291) | `Π I = K ∗ Π` | 🔌 | `RenormalizationGroup.admissible()` raises `NotImplementedError` (Track B3) — K-admissibility is a *model* notion |
| §3.4 Lifting **K** as `D^γ→D^γ` (2452) | Schauder on modelled dist. | ⬜ | analytic; only the homogeneity gain is modeled symbolically |
| §3.5 Building smooth models (2663) | canonical lift of smooth ξ | ⬜ | needs an actual noise realization |
| §4.1 Periodic models (2860) | — | ⬜ | analytic |
| §4.2 Singularity at `x₀=0` (2898) | weighted `D^γ` | ⬜ | analytic |
| §4.3 Non-anticipative ops (3019) | — | ⬜ | analytic |
| §4.4 **Fixed-point solution** (3139) | the abstract solution map | ⬜ | analytic; the project outputs the *renormalized equation*, not the *solution* |
| §5.1 Renorm. structure: definition (3328) | `(T, Δ⁻, ...)` | ✅ | `structures.py RenormalizationStructure` |
| §5.2 Compatible renorm/reg structures (3420) | cointeraction | ✅ | `trees/coproducts.py` δ/δ⁻/δ⁺/Δ/Δ⁺; cointeraction holds **incl. singular β₀=−3/2** (`notes/cointeraction_singular_noise.md`) |
| §6.1 **Free E-multi-pre-Lie algebra** (3765) | the Υ / `F(τ*)` map | ✅ | `renorm/nonlinearity.py` — the elementary differential, all-slots, `∂_p` before `Dⁿ`. The heart of the deliverable |
| §6.2 Modelled-dist. solutions of SPDEs (4302) | — | ⬜ | analytic |
| §6.3 Renorm. over a multi-pre-Lie algebra (4631) | **the renormalized equation** `ThmRenormPDEs` | ✅ | `renorm/equation.py` — assembles `Σ (k(τ)/S(τ)) F(τ*)`. **Golden test: gKPZ 5 counterterms exact** |
| §7 **The BHZ character** (4972) | `k = h∘S'₋` | ✅🟡 | `structures.py bhz_character` builds it **symbolically** (`h` left free); twisted antipode `S'₋` ✅. The **numeric** `h(σ)=𝔼[Πσ](0)` ⬜ (only Wick *parity* done, `renorm/scheme.py`) |
| §8 The manifold of solutions (5157) | — | ⬜ | analytic |
| §9.1 **Rules & extended decoration** (5302) | rules, red nodes, subcriticality | ✅🟡 | `equation/rule.py` derives the `(f,g)`-rule + rule-based subcriticality; red/extended decorations in `trees/tree.py`. **Formal rule completion (Prop 5.21) ⬜** (only needed for adversarial hand-rules) |
| §9.2 Coproducts (5539) | `Δ⁻`, `Δ⁺`, `δ` | ✅ | `trees/coproducts.py` — coassociativity + cointeraction property-tested over the corpus |
| §9.3 Reg/renorm structures (5729) | `G⁻`, `G⁻_ad` | ✅🔌 | `G⁻` is first-class (`build_renormalization_group`, group axioms tested); **`G⁻_ad` 🔌** (admissible subgroup, deferred to model layer) |
| §9.4 Examples (5953) | gKPZ, KPZ tables | ✅ | reproduced **exactly**: gKPZ@−1 (5 cts, tex 6004–6012), gKPZ@−3/2 (43-tree SC table) |
| App. A–D (6217+) | notations / algebra / proofs | n/a | — |

---

## What's there (the algebraic engine — done)

The full pipeline the project set out to build, end to end:

1. **DSL → rule** (`equation/{dsl,rule}.py`): a SymPy SPDE, scope-checked, → an `(f,g)`-rule with
   rule-based subcriticality (BHZ §5.5), generalized to systems / multiple noises / operator order.
2. **Rule → trees** (`equation/generate.py`): saturating enumeration of strongly-conforming
   decorated trees with `|τ|<0`, terminating via the subcriticality budget.
3. **Trees → invariants** (`trees/tree.py`, `core/homogeneity.py`): canonical isomorphism, `S(τ)`,
   homogeneity in the ordered ring `ℚ⊕ℚκ`.
4. **Υ-map** (`renorm/nonlinearity.py`): `F(τ*)`, the elementary differential.
5. **Assembly** (`renorm/equation.py`): the BCCH family `Σ (k(τ)/S(τ)) F(τ*)` with **free**
   constants — the headline output.
6. **The Hopf core** (`core/hopf.py`, `trees/coproducts.py`, `structures.py`): `Δ⁻/Δ⁺/δ/δ⁺`,
   coassociativity + cointeraction (incl. the singular β₀=−3/2 trees), twisted antipode `S'₋`,
   symbolic BHZ character, `RegularityStructure (T,T⁺)`, renormalization group `G⁻`.
7. **da Prato–Debussche lift** (`equation/daprato.py`): supercritical polynomial additive-noise
   equations (Φ⁴₂, Φ⁴₃) → subcritical remainder, unlocking otherwise-rejected inputs.
8. **Rendering / export** (`render/`): trees as shorthand/ascii/`forest`; report in
   txt/md/json/latex; reconstructible structure JSON for an external analytic consumer.

Validated against the paper as oracle (`notes/validation.md`): gKPZ exact at both β₀=−1 and
β₀=−3/2 (a real undercount bug was caught and fixed there).

## What's missing (the analysis — deliberately out of scope)

Everything below needs a *realization of the noise* or *kernel analysis* that a symbolic package
can't supply. None is a defect; each is named in `ROADMAP.md` / CLAUDE.md scope.

| Missing piece | Paper § | Status in repo |
|---|---|---|
| **Models `(Π, Γ)` and modelled distributions `D^γ`** | §2.3–2.4 | ⬜ never started — the single biggest "absent half" of the paper |
| **Reconstruction theorem** | §2.3 (proof App. D, 6418) | ⬜ |
| **Schauder estimates / lifting `K` on `D^γ`** | §3.4 | ⬜ (only the homogeneity *gain* `m−|p|_𝔰` is symbolic) |
| **Admissible smooth models / canonical lift** | §3.3, §3.5 | 🔌 `admissible()` raises |
| **Fixed-point solver** (the actual `u`) | §4.4 | ⬜ — out of scope by design; the project outputs equations, not solutions |
| **Numeric BPHZ constants** `h(σ)=𝔼[Πσ](0)` (Wick + divergent integrals) | §7 | 🟡 combinatorics only: Isserlis/Wick **parity** done (`renorm/scheme.py`), so odd-noise constants vanish; **evaluating the integrals (`BPHZ.numeric_character`) raises** — "the analysis wall" |
| **`G⁻_ad` admissible subgroup** | §9.3 | 🔌 deferred to model layer (a symbolic-only filter would be *unsound* — CLAUDE.md warning) |
| **Formal rule completion** (Prop 5.21) | §9.1 | ⬜ — derived `(f,g)`-rules already saturate; only needed for hand-supplied adversarial rules |
| **Multi-index (Otto–Linares) symbol basis** | — (companion line of work) | ⬜ — the `core/symbol.py` Protocol is the seam, single-impl today |

## Bottom line

> **Does it do justice to the paper?** To the half the paper is *named for* in this project — the
> renormalization combinatorics culminating in the BCCH renormalized equation — **yes, faithfully
> and with paper-exact validation.** To the analytic half — models, reconstruction, Schauder,
> fixed points, and the numeric renormalization constants — **no, and intentionally not**: those
> require probabilistic/analytic input outside a symbolic engine's reach, and each is left as an
> honest, named socket rather than a silent omission.

The honest one-line scope statement: *`regstruct` computes **which** counterterms appear and
**what** their vector fields are (symbolically, with free constants); it does not compute the
**numeric** constants, nor does it solve or analyze the resulting SPDE.*
