# `counterterms` — symbolic renormalized equations for singular SPDEs

**High-level plan.** A Python package that takes a *singular stochastic PDE* as input and
returns the *family of renormalized equations* as output, following Bailleul & Hoshino,
*"A tourist's guide to regularity structures and singular stochastic PDEs"*
(arXiv:2006.03524, v3; `references/tourist_guide.{pdf,tex}`).

This is a planning document. **No code is written until the plan is approved.**

---

## 1. Goal, in one line

> Input: a subcritical singular SPDE `(∂_t − Δ + 1)u = f(u)ζ + g(u,∂u)`.
> Output: the equation
> `(∂_t − Δ + 1)u⁽ᵏ⁾ = f(u⁽ᵏ⁾)ζ + g(u⁽ᵏ⁾,∂u⁽ᵏ⁾) + Σ_{τ∈𝓑, |τ|<0} (c_τ) · F(τ*)(u⁽ᵏ⁾,∂u⁽ᵏ⁾)`
> i.e. the original PDE plus one counterterm per negative-homogeneity decorated tree,
> with free renormalization constants `c_τ = k(τ)/S(τ)`. This is **Theorem `ThmRenormPDEs`**
> of the paper (tex line 4914), the BCCH renormalized-equation formula.

The "family" is the set of these equations as `k` ranges over the renormalization group
`G⁻_ad` — equivalently, as the finitely many constants `{c_τ}` range freely over ℝ. The
specific *canonical* (BPHZ/BHZ) values of `c_τ` require probabilistic input (Gaussian/Wick
expectations of the noise) and are **out of scope** for the symbolic core (see §8 v3).

---

## 2. The key realization that keeps this lazy

**To emit the family of renormalized equations with free constants, you do *not* need the
Hopf-algebra coproducts, the twisted antipode, or the BHZ character.** Those machines exist
only to compute the *canonical numerical values* of the constants `c_τ` — which need the
analytic/probabilistic theory we are deliberately not implementing (paper tex lines 99/118
of our research digest; §5 "Renormalization structures" and §7 "BHZ character").

The symbolic family is fully determined by three combinatorial ingredients:

1. the finite set `𝓑_{<0}` of decorated trees with negative homogeneity (from the *rule*);
2. the symmetry factor `S(τ)` of each tree;
3. the elementary-differential map `τ ↦ F(τ*)` (the SPDE analogue of a B-series elementary
   differential).

So the MVP (§8) skips the hardest combinatorics (the two cointeracting coproducts `Δ, Δ⁻`)
entirely and still delivers the user's stated goal. The coproducts/antipode become an
optional v3 plugin only needed for the canonical constants.

---

## 3. The mathematics (precise)

### 3.1 Spacetime, scaling, homogeneity

- Spacetime `ℝ × ℝᵈ`, coordinate `x = (x₀, x₁,…,x_d)` = (time, space). Multi-indices live in
  `ℕ^{1+d}`; `e_i` is the i-th unit multi-index (`e₀` = time).
- Parabolic scaling `𝔰 = (2,1,…,1)`. Scaled length of a multi-index: `|k|_𝔰 = 2k₀ + k₁+…+k_d`.
- **Homogeneity convention (paper-authoritative, tex line 2827; verified against the example
  tables).** The noise `ζ ∈ C^{β₀}`; the noise symbol has `|Ξ| = β₀` *directly* (β₀ is the
  Hölder exponent, a negative number).
  - `|Ξ| = β₀`           (one per noise; for several noises, `|Ξ_j| = β_j`)
  - `|X^n| = |n|_𝔰`      (polynomial / node decoration)
  - `|I_p(τ)| = |τ| + 2 − |p|_𝔰`   (heat-kernel integration; the gain `+2` is the parabolic
    Schauder gain, **hardcoded in the paper's theory**; `p` is the derivative on the kernel)
  - naive homogeneity of a tree: `|τ|′ = Σ_{edges} (2 − |p_e|_𝔰) + Σ_{nodes} |n_v|_𝔰 + (#noise leaves)·β₀`

  Reference check (KPZ, d=1, `β₀ = −1−κ`, tex lines 6004–6012):
  `|∘| = β₀`, `|∘ with n=(0,1)| = β₀+1`, `|I_{(0,1)}(∘)| = β₀+2−1 = β₀+1`. ✔

> ⚠️ **The single most error-prone input.** Do *not* introduce a separate "α = regularity+2"
> convention; set `|Ξ| = (regularity supplied by user)`. An off-by-2 here silently produces
> wrong subcriticality verdicts. (One research agent used the α convention and was wrong.)

- **Ordered ring for comparisons.** Homogeneities are of the form `q + m·β₀` with `q ∈ ℚ`,
  `m ∈ ℤ≥0` (= number of noise leaves), and `β₀ = r − κ` with `r ∈ ℚ` and `κ > 0` an
  *infinitesimal*. So every homogeneity reduces to a pair `(rational, κ-coefficient)`.
  Comparison to 0 is lexicographic: compare the rational part first, break ties by the sign
  of the κ-coefficient (κ>0). **Never use floats** — critical trees sit exactly at homogeneity
  `−4κ, −2κ, …` (tex lines 6066, 6140) and a float β₀ would mis-classify them. (Adversarial
  review flagged this as essential for correctness.)

### 3.2 Decorated trees — the central data structure

A basis element of `T` (and of the multi-pre-Lie space `V`) is a **rooted, non-planar,
decorated tree** `τ = b^n ⋆ ⨉ᵢ I_{p_i}(τ_i)` (tex line 4506 / digest):

- **node type** `b ∈ {●, ∘₁,…,∘_m}` — `●` = "polynomial / nonlinearity" node (carries `g`),
  `∘_j` = "noise j" node (carries `f_j`). (Plus `red ●^{0,α}` extraction nodes used only in
  the v3 BHZ plugin.)
- **node decoration** `n : N_τ → ℕ^{1+d}` — the polynomial factor `X^n` at a node.
- **edge type + edge decoration** `(t_e, p) ∈ 𝔗_e × ℕ^{1+d}` — each child edge is an
  integration operator `I^{(t_e)}_p` (kernel `∂^p K`). For a single scalar equation `𝔗_e`
  is a singleton; **for systems, `𝔗_e` carries the equation/component identity** (which
  operator `I^{(a)}` planted the subtree). This is the correct place for component identity —
  *not* the node type (adversarial-review correction).
- restriction in this paper: `|p|_𝔰 ≤ 1`, i.e. only `I` and `∂_i I` exist (tex line 2222).

**Canonicalization up to isomorphism is mandatory and is the #1 correctness concern.**
Canonical key: `(node_type, n, sorted multiset of (edge_type, p, canonical(child)))`. Equality,
dict keys, dedup during enumeration, and `S(τ)` all depend on it.

### 3.3 The rule (derived from the SPDE)

A **rule** assigns, to each node type `b`, the set of allowed child-edge multisets (tex
lines 5306–5384). It is read off the right-hand side of the equation written in fixed-point
form `v = I(f(v)Ξ + g(v,∂v)) + (T_X)` (tex line 5308):

- each monomial of the nonlinearity that multiplies a noise → an allowed configuration under a
  `∘` node; each noise-free monomial → a configuration under a `●` node;
- a factor `∂_{x_i} u` contributes a child edge `I_{e_i}` (double line); a factor `u`
  contributes `I_{(0)}`; powers give repeated children.

Definitions we need: **conforming** (rule holds at every node except maybe the root),
**strongly conforming** `SC` (holds everywhere), **normal**, **subcritical** (for every γ,
only finitely many `SC` trees have `|τ|′ < γ` — this is what makes the algorithm terminate),
**complete** (closed under the contractions used by `Δ⁻`; existence guaranteed by BHZ Prop 5.21).

For the MVP we derive a rule from `(f,g)` and **check subcriticality**; if it fails we error.
We also let the user supply an explicit rule as an escape hatch when auto-derivation cannot be
certified normal+complete (adversarial-review correction).

### 3.4 Tree enumeration

Generate `SC` trees by saturating from the noise/polynomial generators under the rule, pruning
by the homogeneity cutoff (`|τ|′ < 0` for counterterms; `< γ` for the full structure). Terminates
iff subcritical. Output `𝓑_{<0} = { τ ∈ SC : |τ|′ < 0 }`. **`𝓑_{<0}` includes bare decorated
primitives** `∘^n` and `●^n` (when their homogeneity is `< 0`) — these are the KPZ `k(∘)` and
`k(∘1)` counterterms (tex line 6004); do not describe `𝓑_{<0}` as only "`I_p`-built trees".

### 3.5 Symmetry factor

`S(τ) = n! · Πⱼ S(σⱼ)^{mⱼ} · mⱼ!` where the distinct subtrees `σⱼ` appear with multiplicity
`mⱼ` and `n!` is the multi-index factorial of the root decoration (tex line 3982 / line 75 of
digest). Recursive; memoized on the canonical key. Validate against automorphism counts on the
gKPZ table, including the factor-2 cases (`S = 2`).

### 3.6 The elementary differential `F(τ*)` (the Υ map) — the counterterm engine

Jet variables `u_k`, `k ∈ ℕ^{1+d}` (with `u₀ = u`, `u_{e_i} = ∂_{x_i}u`). Two operators on
smooth functions of the jet:

- slot derivative `∂_p F = ∂F/∂u_p`;
- total derivative `D_i F = Σ_k u_{k+e_i} ∂_k F`,  `D^n = Π_i D_i^{n_i}`.

**Base cases (user data):**

```
F(∘_j*) = f_j(u₀)                                         # noise-j coefficient
F(●*)   = Σ_{i,j} g₂^{ij}(u₀) u_{e_i} u_{e_j}             # at most QUADRATIC in ∂u
        + Σ_i g₁^i(u₀) u_{e_i} + g₀(u₀)                   #   (Assumption D2)
F(red ●^{0,α}*) = 0
```

**Recursion (tex lines 4524 / 4915; the morphism property):** for `τ = b^n ⋆ ⨉ᵢ I_{p_i}(τ_i)`,

```
F(τ*) = ( Πᵢ F(τ_i*) ) · ( D^n Πᵢ ∂_{p_i} ) F(b*)
```

i.e. differentiate the root nonlinearity once per child edge (`∂_{p_i}`) and once per node
decoration (`D^n`), then multiply by the children's differentials. **Order matters**: apply the
`∂_{p_i}` first, then `D^n` (`∂_p` and `D_i` do not commute; tex line 112). SymPy performs the
differentiation. Final substitution `u₀ → u(x)`, `u_{e_i} → ∂_{x_i}u(x)` produces the PDE term.

> The `∂_{p_i}` must hit **both** the function-argument slot and the derivative-argument slots
> of `g` (e.g. `g(u, ∂_x u)` is differentiated in both). A scalar `sympy.diff(f, u)` is *not*
> enough (adversarial-review correction).

**Correctness invariant** (good for tests): `F(τ* ⌢_e σ*) = F(τ*) ▷_e F(σ*)` where
`G ▷_e H = G·∂_{u_e}H`, and `F(↑_i τ*) = D_i F(τ*)` (tex line 107).

### 3.7 The renormalized equation (output)

```
(∂_t − Δ + 1) u⁽ᵏ⁾ = f(u⁽ᵏ⁾)ζ + g(u⁽ᵏ⁾,∂u⁽ᵏ⁾) + Σ_{τ∈𝓑, |τ|<0} (k(τ)/S(τ)) F(τ*)(u⁽ᵏ⁾,∂u⁽ᵏ⁾)
```

We expose `c_τ := k(τ)/S(τ)` as free symbolic constants. `#counterterms = #𝓑_{<0}` =
`dim G⁻_ad` (the renormalization group is a finite-dim Lie group; manifold-of-solutions
picture, tex line 5184). Restricting to `G⁻_ad` (vs the bigger `G⁻`) imposes the *admissibility*
constraints `k(I_p τ)=0`, `k(X^n ⋆ τ)=0` — i.e. drop counterterms whose tree is itself an
integral or a nonconstant-polynomial multiple of a smaller tree (tex line 47 of digest). We will
flag both bases (raw `G⁻` and admissible `G⁻_ad`) but default to `G⁻_ad`.

### 3.8 (v3, optional) the canonical BHZ character

If we ever want the *values* of `c_τ`: build the extraction coproduct `δ` (tex `D⁻` formula,
line 5539+), the **negative twisted antipode** `S'_-` recursion
`S'_-(P_-τ) = −M_-(S'_-⊗Id)(δτ − P_-τ ⊗ ●^{0,|τ|})` (tex line 121), and the character
`k^ζ = h^ζ ∘ S'_-` with `h^ζ(τ) = 𝔼[Π^ζ τ](0)`. **`h^ζ` needs Wick/Gaussian expectations of
the specific noise — that is the analytic/probabilistic part we do not implement.** The plugin
would, at most, return the *symbolic forest structure* (which constants are tied / which trees
diverge), leaving the numbers symbolic. It also needs the richer **extended-decoration** tree
type (red nodes with `𝔬 : N^red → ℤ[β₀]`), distinct from the plain MVP tree (adversarial-review
correction). Keep it a separate, opt-in module.

---

## 4. SPDE input representation (the design the user asked us to think hard about)

A small **Python DSL on top of SymPy**. The user writes the equation almost as they would on
paper; the package *derives* the rule from it (it does not ask the user to specify trees/rules).

```python
from counterterms import Unknown, Noise, Parabolic, SPDE, kappa
from sympy import Function, Derivative, Symbol

# --- symbols -------------------------------------------------------------
u   = Unknown("u")                       # scalar unknown (systems: pass several / vector)
xi  = Noise("xi", regularity=-1 - kappa) # |Ξ| = regularity (Hölder exponent), kappa>0 infinitesimal
x   = u.space_coords                      # (x1,) here; d inferred from `Parabolic(dim=...)`
f, g = Function("f"), Function("g")       # abstract smooth coefficients (or concrete exprs)

# --- the equation --------------------------------------------------------
eq = SPDE(
    operator   = Parabolic(dim=1, mass=1),          # L = ∂_t − Δ + 1; scaling (2,1); Schauder gain +2
    unknown    = u,
    rhs        = f(u)*xi + g(u)*Derivative(u, x[0])**2,
    noises     = [xi],
)

fam = eq.renormalize()        # -> RenormalizedEquation
print(fam.latex())            # the full family with free constants c_τ
fam.counterterms              # [Counterterm(tree, homogeneity, S, c_symbol, vector_field), ...]
```

**Exact data the user supplies:**

| field | type | meaning |
|---|---|---|
| `operator` | `Parabolic(dim, mass=…)` | the differential operator; fixes `d`, scaling `𝔰`, Schauder gain `m=2`. Whitelisted to 2nd-order parabolic (the proven theory); fields `m, 𝔰` are carried generally for forward-compat but validated. |
| `unknown(s)` | `Unknown` (or list) | the solution field(s); for systems each is a sector with its own operator `I^{(a)}`. |
| `noises` | `[Noise(name, regularity)]` | each noise's Hölder regularity (the `β_j`); `\|Ξ_j\| = β_j`. |
| `rhs` | SymPy expr | the nonlinearity in `u`, its spatial derivatives, and the noise(s). Must be **affine in each noise** and **≤ quadratic in `∂u`** (Assumption D2). |
| `rule` *(optional)* | `Rule` | explicit override when auto-derivation can't be certified. |

**Parsing the `rhs`:** split into the per-noise multipliers `f_j(u)` (coefficient of each noise,
checked affine) and the noise-free part `g(u,∂u)`; read the rule from the monomials; build the
base elementary differentials `F(∘_j*) = f_j`, `F(●*) = g`.

**Validation / explicit out-of-scope rejection list** (constructor raises clear errors —
adversarial-review correction):

- `β₀ ≤ −2` → reject, or route through a da Prato–Debussche pre-pass (v3) introducing one noise
  symbol per noise power (this is how the paper handles Φ⁴₃, tex line 2034). **Φ⁴₂/Φ⁴₃/sine-Gordon
  are NOT handled by the direct rule reader.**
- nonlinearity not affine in the noise (general `star`-products of noise) → reject.
- `g` more than quadratic in `∂u` (outside Assumption D2 the BCCH coherence theorem does not
  hold; the family would be unsound) → reject.
- a singular derivative factor with `|n|_𝔰 > 1` (e.g. `∂_xx u`, `∂_t u` in the nonlinearity) →
  reject.
- quasilinear `L(u)·u`, non-parabolic / higher-order `L` → reject (warn that Schauder is unproven).
- subcriticality check fails → reject.

---

## 5. Output representation

`RenormalizedEquation`:

- `.base` — the original equation as a SymPy `Eq`.
- `.counterterms` — list of `Counterterm`:
  - `tree` (DecoratedTree, with a pretty/LaTeX form),
  - `homogeneity` (ordered-ring value),
  - `symmetry_factor` `S(τ)`,
  - `constant` — free SymPy `Symbol` `c_τ` (`= k(τ)/S(τ)`),
  - `vector_field` — `F(τ*)(u,∂u)` as a SymPy expression.
- `.as_sympy()` — the assembled family RHS (a single SymPy `Eq` with the `c_τ` free).
- `.latex()` — rendered family, grouped by homogeneity (mirrors the paper's tables).
- `.constants` — the chart of `G⁻_ad` (one free real per negative tree).

The "family" *is* this object read with `{c_τ}` free over ℝ. Document clearly that pinning them
to the canonical BPHZ values is a separate probabilistic step (v3 / out of scope).

---

## 6. Reuse vs build (research conclusion)

**No existing tool does the end-to-end job** (SPDE → renormalized BCCH/BHZ equation). Verified
across GitHub/GitLab/Zenodo and the whole RS literature (Bruned, Chandra, Chevyrev, Hairer,
Bailleul, Hoshino, Zambotti, Manchon, Foissy, Otto–Linares–Tempelmayr): theory papers only, no
companion code. So the differentiating core is genuinely greenfield. The substrate, however, is
well-trodden.

| Concern | Decision | What |
|---|---|---|
| Symbolic differentiation (Υ-map) + final rendering | **REUSE** | **SymPy** (hard dependency). Mirrors how `BSeries.jl` plugs in SymPy/Symbolics. Keep tree combinatorics *out* of SymPy (dict-of-trees / free module); convert to SymPy only at the leaves and at output. |
| Decorated-tree datatype + canonical form + `S(τ)` | **BUILD** | A focused pure-Python module (~a few hundred lines). It must carry edge-type + edge/node/extended decorations that no library has. |
| Rule → subcriticality-bounded enumeration | **BUILD** | No library does this; spec = BHZ arXiv:1610.08468 + the paper §9. |
| Elementary-differential recursion | **BUILD** (port idea) | Port the *pattern* from `BSeries.jl`'s `elementary_differential` / `modified_equation`. |
| Coproducts `Δ, Δ⁻` + twisted antipode (v3) | **BUILD** | The BHZ extraction/contraction coproducts are *not* the Connes–Kreimer/Grossman–Larson coproducts; must be written from the paper. |

**Oracles / design references (NOT runtime dependencies):**

- **`kauri`** (Python, Apache-2.0, maintained 2026) — rooted trees + Connes–Kreimer coproduct +
  elementary differentials + SymPy backend. Use as a *cross-check oracle* for the undecorated
  CK coproduct and elementary differentials; possibly cannibalize small pieces. Not forked
  wholesale (lacks decorations/edge-types/rule/BHZ — the hard parts).
- **SageMath `sage.combinat`** (`FreePreLieAlgebra`, `LabelledRootedTree`, `GrossmanLarson`) —
  a *correctness oracle* for the pre-Lie grafting product and tree canonicalization. **Not a
  runtime dep** (GPL, not `uv`/pip-friendly; conflicts with the user's `uv` environment).
- **`RootedTrees.jl` / `BSeries.jl`** (Julia, MIT) — architectural blueprint for elementary
  differentials, the modified/modified-equation pattern, and two-interacting-Hopf-algebra
  bookkeeping (`partition_forest`/`all_splittings` mirror the BHZ extraction shape). Reference,
  not dependency.

**Stack decision:** **Python + SymPy, pip/`uv`-installable, pure-Python combinatorial core.**
Rationale: (i) matches the user's `uv`/Python environment; (ii) Sage's reusable bits are either
things we must specialize anyway (BHZ coproducts ≠ CK/GL) or small enough to build; (iii) SymPy
is exactly right for the Υ-map and rendering. This is the lazy *and* correct choice — minimal
deps, reuse the CAS instead of reinventing it, build only the novel core. *(This is the main
decision I'd like confirmed at approval — see §11.)*

---

## 7. Architecture & modules (lean)

```
counterterms/
  trees.py        # DecoratedTree, canonical form, S(τ), Homogeneity (ordered ring ℚ+ℚ·κ + β₀)
  rule.py         # SPDE-RHS → Rule; subcriticality check; SC-tree enumeration with cutoff
  elementary.py   # jet vars, ∂_p, D_i, base F(b*), the F(τ*) recursion (SymPy-backed)
  spde.py         # the input DSL: Unknown, Noise, Parabolic, SPDE; rhs parsing + validation
  renormalize.py  # orchestration: SPDE → RenormalizedEquation (assemble counterterms)
  bhz.py          # (v3, opt-in) extended-decoration trees, δ/Δ⁻, twisted antipode S'_-, character
tests/
  test_golden_gkpz.py   # the 5 counterterms + the KPZ homogeneity table
  test_symmetry.py      # S(τ) vs automorphism counts (incl. factor-2 cases)
  test_homogeneity.py   # ordered-ring comparisons at κ-multiples
references/         # the paper (already present)
notes/             # this plan
```

~5 core modules; the split tracks the math, not speculative layering. `bhz.py` is opt-in.

---

## 8. Phased plan

**MVP (v1) — delivers the stated goal for the core class.**
Scope: scalar unknown, single noise, 2nd-order parabolic `L`, `β₀ ∈ (−2,−1)`,
`g` ≤ quadratic in `∂u`, `|n|_𝔰 ≤ 1`. Covers **gKPZ, KPZ, gPAM, PAM**.
Pipeline: DSL parse + validate → rule → `SC` enumeration of `𝓑_{<0}` → `S(τ)` → `F(τ*)` →
assemble family with free `c_τ`. **No coproducts/antipode.** Golden-tested against the paper's
gKPZ example (the exact 5 counterterms) and KPZ homogeneity table.

**v2 — generality.**
Systems (equation identity on **edge type** `𝔗_e`; vector/matrix-valued `F`-morphism with
cross-component partials `∂_{u_j}F_i`), multiple noises (one `∘`-type per noise, base map
`{●:g, ∘_j:f_j}`), general scaling/operator order carried in the homogeneity arithmetic (with
scope warnings outside the proven parabolic regime), `G⁻` vs `G⁻_ad` toggle.

**v3 — canonical structure (opt-in, partly out of scope).**
The extended-decoration tree type, coproducts `Δ, Δ⁻` and their cointeraction (unit-tested via
coassociativity), the twisted antipode `S'_-`, the symbolic BHZ forest structure. da Prato–
Debussche pre-pass for `β₀ ≤ −2` (Φ⁴₂, Φ⁴₃, sine-Gordon). The *numerical* BPHZ constants need
Wick/Gaussian expectations — **explicitly out of scope** unless a noise-cumulant module is added.

---

## 9. Validation — golden tests from the paper

1. **gKPZ, d=1, `(∂_t−Δ+1)u = f(u)ζ + g(u)(∂_xu)²`, `ζ∈C^{−1−κ}`** (tex 5996–6012). The family
   must reproduce *exactly* (these are the acceptance test for the MVP):
   - `∘`            → `k(∘) f(u)`
   - `∘` with `n=(0,1)` → `k(∘1) f'(u) ∂_x u`
   - `I_{(0,1)}(∘)` under `●` → `2 k(τ) f(u) g(u) ∂_x u`   *(note the factor 2)*
   - `I(∘)` under `∘` → `k(τ) f(u) f'(u)`
   - two `I_{(0,1)}(∘)` under `●` → `k(τ) f²(u) g(u)`
2. **KPZ homogeneity table** (tex 6028–6063): the count of strongly-conforming trees per
   homogeneity row (`β₀, 2β₀+2, 3β₀+4, β₀+1, 4β₀+6, 2β₀+3`).
3. **`S(τ)`** vs hand-computed automorphism counts on the above, incl. the `S=2` and the
   `(i=j)` doubling `F((I_{e_i}Ξ ⋆ I_{e_i}Ξ)*) = 2 g₂^{ii} f²` (tex line 4349).
4. **`F(τ*)` spot checks**: `F((Ξ ⋆ I(Ξ))*) = f f'`, `F((I_{e_i}Ξ ⋆ I_{e_j}Ξ)*) = g₂^{ij} f²`.

---

## 10. Correctness traps (collected from the adversarial review)

1. `|Ξ| = β₀` (the regularity itself); **no off-by-2 α convention**.
2. Homogeneity comparisons in the **ordered ring** `ℚ + ℚ·κ` (+ `β₀`), never floats — critical
   trees sit at `−kκ`.
3. **Edge type** carries equation/component identity for systems; node type stays `{●,∘_j,red}`.
4. `𝓑_{<0}` **includes bare primitives** `∘^n, ●^n` (the `k(∘)`, `k(∘1)` terms).
5. The Υ-map differentiates **all** slots of `g` (function *and* derivative arguments), not just
   `∂/∂u`.
6. Coefficient is `k(τ)/S(τ)`, **never `k(τ)` alone**.
7. Subcriticality uses the **naive** homogeneity `|τ|′` (the `𝔬`-decoration is 0 in this phase).
8. Canonical tree isomorphism is load-bearing — get it wrong and the whole algebra is silently
   inconsistent.
9. Φ⁴₃ etc. (`β₀ ≤ −2`) are **out of the direct rule reader** — reject or DPD-preprocess.
10. Keep the regularity structure `𝒯` and the renormalization structure `𝒰` distinct (two
    gradings on the same trees) — only relevant once v3 coproducts exist.

---

## 11. Open decisions for you (at approval)

1. **Stack** — confirm **Python + SymPy, pure-Python core** (recommended, matches your `uv`
   env), vs. building on **SageMath** (richer primitives, heavier/GPL/non-`uv`), vs. **Julia**
   (best tree substrate, smaller ecosystem here). *Recommendation: Python + SymPy.*
2. **MVP scope** — confirm the v1 class (scalar, single noise, parabolic, gKPZ/KPZ/gPAM/PAM with
   free constants). Systems & multi-noise in v2. OK to defer?
3. **Canonical constants** — confirm that computing the *numerical* BPHZ values (needs Gaussian/
   Wick expectations) is out of scope; the package emits free symbolic constants. OK?
4. **Package name** — resolved: the package is `counterterms`.

---

## 12. References

- Bailleul & Hoshino, *A tourist's guide to regularity structures and singular stochastic PDEs*,
  arXiv:2006.03524 v3 — `references/tourist_guide.{pdf,tex}`. Key tex anchors: rule/decorated
  trees §9.1 (5302), coproducts §9.2 (5539), examples §9.4 (5953), multi-pre-Lie §6.1 (3702),
  renormalized equation Thm `ThmRenormPDEs` (4901), Υ-map (4559), BHZ character §7 (4972).
- Bruned, Hairer, Zambotti, *Algebraic renormalisation of regularity structures*, arXiv:1610.08468
  (the rule/coproduct formalism; completion Prop 5.21).
- Bruned, Chandra, Chevyrev, Hairer, *Renormalising SPDEs in regularity structures*,
  arXiv:1711.10239 (BCCH; the renormalized-equation formula this package emits).
- Bailleul & Bruned, *Renormalised singular stochastic PDEs* (multi-pre-Lie presentation).
- Reference code (oracles only): `kauri` (PyPI), SageMath `sage.combinat`, `BSeries.jl` /
  `RootedTrees.jl`.
