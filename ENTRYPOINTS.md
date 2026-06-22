# ENTRYPOINTS — how to read this code

A guided reading order for understanding `regstruct` in detail. Each file has a single job;
this is the order that makes them click. (For "what can I *do*, and which module?", see the
capability map in [`README.md`](README.md); this doc is for reading the *code*.)

If you only read two things: **`regstruct/api.py`** (the whole Phase-1 pipeline in ~25 lines) and
**`tests/test_goldens.py`** (the input/output contract on the paper's worked examples).

For the mathematics behind any file, see [`notes/initial_plan.md`](notes/initial_plan.md); for
why the modules are split the way they are, see [`notes/architecture.md`](notes/architecture.md).
Formulas are cited by line number in `references/tourist_guide.tex`.

---

## The 30-second mental model

```
SPDE (SymPy) ──parse──▶ Signature + base nonlinearities
             ──rule + enumerate──▶ {decorated trees τ : |τ| < 0}
             ──S(τ), F(τ*)──▶ counterterms
             ──assemble──▶ (∂_t−Δ+1)u = … + Σ_τ (k_τ/S(τ)) F(τ*)
```

Everything below fills in one arrow.

---

## Read in this order

### 1. `tests/test_golden_gkpz.py` — what the library promises
The gKPZ example (tex 5996–6012): an SPDE in, five specific counterterms out. Read this first to
see the **input format** and the **exact expected output**. The assertions encode the three
things every counterterm carries: homogeneity `|τ|`, symmetry factor `S(τ)`, and the elementary
differential `F(τ*)`.

### 2. `regstruct/api.py` — the pipeline spine
`renormalize(spde)` is the whole thing: parse → generate trees → for each, compute `F(τ*)` and
`S(τ)` → package. Twenty lines; everything else is a subroutine called here. Note the one
non-obvious line: trees with `F(τ*) = 0` are dropped (Assumption D2 safety).

### 3. `regstruct/equation/dsl.py` — input and parsing
The user-facing classes (`Unknown`, `Noise`, `Parabolic`, `SPDE`) and `build_context`, which
turns the SymPy right-hand side into (a) the **base nonlinearities** `F(∘_j*) = f_j`,
`F(●*) = g` in jet variables, and (b) the structural **rule** (which child edges each node type
admits, with caps). This is also where **scope is enforced** — read the `raise ValueError(...)`
lines to see exactly what is in and out of scope (affine-in-noise, `g` ≤ quadratic in `∂u`,
`β₀ ∈ (−2,0)`, `|p|_𝔰 ≤ 1`).

### 4. `regstruct/core/signature.py` — the vocabulary
The `Signature` is the single object every algorithm is parametric over: dimension, scaling,
the integration operator's Schauder order, the noise regularities, and the rule. Crucially,
`node_homogeneity`, `edge_gain` (`= m − |p|_𝔰`) and `scaled` live here — the homogeneity bookkeeping.
Read this to see how "scalar vs system" and "one vs many noises" become *data*, not code branches.

### 5. `regstruct/core/homogeneity.py` — the ordered ring
`Homogeneity` is `std + kap·κ` with `κ` a positive infinitesimal; `is_negative()` and the
ordering are the load-bearing detail (critical trees sit at homogeneity `−kκ`, so comparisons
must not use floats). `Scaling` gives the scaled degree `|k|_𝔰`. Small file, but understanding
`is_negative` is essential — it decides which trees become counterterms.

### 6. `regstruct/trees/tree.py` — the central datatype
`DecoratedTree` (`τ = b^n ⋆ ⨉ᵢ I_{p_i}(τ_i)`): node type + decoration, and a canonically-sorted
multiset of child edges. The three methods to study:
- `_sortkey` / the `tree()` builder — **canonicalisation**, on which equality, hashing and
  dedup all depend (the #1 correctness invariant);
- `homogeneity(sig)` — the recursion `|τ| = |b| + |n|_𝔰 + Σ (edge_gain + |subtree|)`;
- `symmetry_factor()` — `S(τ) = n!·Πⱼ S(σⱼ)^{mⱼ}·mⱼ!` (tex 3982).

### 7. `regstruct/core/jets.py` — jet variables (tiny)
`jet(k)` maps a spacetime multi-index to a SymPy symbol `u_k`; `u_(0,…,0)` is the unknown,
`u_{e_i}` its `∂_{x_i}` derivative. These are the variables the nonlinearities and `F(τ*)` are
written in. Read it right before §8–9.

### 8. `regstruct/equation/generate.py` — which trees exist
`generate_counterterms(sig)` enumerates the strongly-conforming trees with `|τ| < 0` by a
**budget-pruned fixpoint DFS**. The comments explain termination: every uncapped edge (`I_0`)
adds positive homogeneity when `β₀ > −2`, so it is bounded by the budget; every edge that can
add `≤ 0` is a derivative edge with a finite cap (this is the "`g` quadratic in `∂u`" cutoff).
The `break` in `_emit` relies on the child atoms being sorted by homogeneity contribution.

### 9. `regstruct/renorm/nonlinearity.py` — the counterterm engine (Υ map)
`elem_diff(τ, base_F, width)` computes `F(τ*) = (Πᵢ F(τ_i*))·(D^n Πᵢ ∂_{p_i}) F(b*)`
(tex 4524 / 4915). The two operators to understand: `∂_p` (a plain slot derivative) and the
total derivative `D_i = Σ_k u_{k+e_i} ∂_k` (`_D`). Order matters: the `∂_{p_i}` are applied
*before* `D^n`, and they hit **all** slots of `g`. This is the only place SymPy does real work.

### 10. `regstruct/renorm/equation.py` — the output
`Counterterm` (tree, `|τ|`, `S(τ)`, `F(τ*)`, free constant `k_τ`; `.coefficient = k_τ/S(τ)`)
and `RenormalizedEquation` (`.counterterms`, `.summary()`, `.counterterm_rhs()`). `summary()`
renders jet variables back to `u` and its derivatives for display.

---

## Two cross-cutting things worth grepping

- **The `Signature` is threaded everywhere** — `tree.homogeneity(sig)`, `generate_counterterms(sig)`.
  It is the seam that keeps the algorithms equation-agnostic.
- **Jet symbols are the lingua franca** between the parser (`dsl._to_jet`), the engine
  (`nonlinearity`) and the renderer (`equation._display_subs`). Grep `jet(` to follow them.

## Beyond the Phase-1 spine (now built)

The reading order above is the Phase-1 pipeline (`SPDE → counterterm family`). The algebraic
layer (Phase 3) sits on top of the same `Signature`/`DecoratedTree` core:

- `trees/coproducts.py` — extraction `δ`/`δ⁻`, recentering `Δ`/`Δ⁺`, twisted antipode `S'₋`
  (the load-bearing, most intricate file; the cointeraction subtlety is in `notes/cointeraction_singular_noise.md`).
- `core/hopf.py` — generic `convolve`/`antipode`/`comodule_action` over plain coproduct/mul maps.
- `structures.py` — bundles the above into `RegularityStructure (T,T⁺)`, `RenormalizationStructure`
  (BHZ character), and the group `G⁻`.
- `renorm/scheme.py` — the canonical-character Wick/parity machinery (symbolic; integrals unbuilt).
- `equation/daprato.py` — the da Prato–Debussche lift; `render/export.py` — JSON export.

For *what each does and how to call it*, see the capability map in [`README.md`](README.md).
What is **not** built: numeric BPHZ constants (Gaussian integrals) and the `G⁻_ad` reduction —
the analytic/probability wall (`notes/phase4_plan.md`).
