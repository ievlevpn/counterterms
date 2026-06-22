# `counterterms` — architecture

High-level structure of a **modular, extensible, general** library for regularity structures.
Companion to `notes/initial_plan.md` (which fixes the mathematics, conventions, and golden tests).
This document fixes the *module boundaries and interfaces*; the plan's pipeline and math
conventions are unchanged and authoritative.

**Design decisions locked (from approval):**
- **Basis-generic algebra, trees the one concrete basis.** The algebra/Hopf machinery is written
  over a `Symbol` interface; decorated trees are the only implementation now; multi-indices can be
  added later without touching the algebra.
- **Probabilistic layer is a deferred seam.** The library emits *free symbolic* renormalization
  constants. A `NoiseLaw`/character interface exists as a boundary; no Wick engine is built yet.
- **Full algebraic-renormalization machinery is first-class, built after the equation pipeline.**
  Coproducts `Δ, Δ⁻`, twisted antipode, BHZ character, and the regularity/renormalization
  structures are core modules — implemented in a later phase, but their boundaries are fixed now.

---

## 1. The spine: three orthogonal separations

Modularity here is not arbitrary file-splitting; it follows three seams that are *intrinsic to
the mathematics*. Get these three right and everything else is mechanical.

1. **Basis-agnostic algebra ⟂ basis-specific combinatorics.**
   The Hopf-algebra theory — graded connected bialgebra, coproduct application, convolution of
   characters, antipode recursion, comodule coaction — *does not know what a symbol is*. It works
   identically for decorated trees and for multi-indices. So the algebra core depends only on a
   `Symbol` protocol and treats a coproduct as *data* (a linear map it is handed). The tree module
   provides the concrete symbol and **builds** the coproducts; the algebra core consumes them.

2. **Equation-independent structure ⟂ equation-specific rule.**
   The regularity-structure machinery (gradings, coproducts, characters) is generic. The *only*
   thing a particular SPDE injects is a **rule** (which symbols exist) and a **nonlinearity**
   (how trees map back to PDE terms). So `structures/` is equation-agnostic, parameterized by a
   `Rule`; the equation layer produces the rule.

3. **Symbolic/combinatorial ⟂ analytic/probabilistic.**
   Everything the library computes is combinatorial. The analytic model (`Π, Γ`, reconstruction,
   Schauder) and the *numerical* renormalization-constant values (Wick expectations) live
   **outside** the library, reached only through narrow seams (`NoiseLaw`, structure export). The
   library never does analysis; it stops cleanly at the symbolic boundary.

The single object that carries all "generality axes" across these seams is the **`Signature`**
(§3.2): the typed vocabulary the whole library is parametric over. Scalar-vs-system, one-vs-many
noises, parabolic-vs-other operators, the dimension — all are *data in a Signature*, not branches
in the code.

---

## 2. The layered stack

Dependencies point strictly downward; no upward or cyclic imports.

```
            ┌───────────────────────────────────────────────────────────────┐
 render/    │ LaTeX · unicode · SymPy export · structure/model export (handoff)│
            ├───────────────────────────────────────────────────────────────┤
 api/       │ renormalize(spde) · build_structure(spde) · homogeneity_table … │
            ├───────────────────────────────────────────────────────────────┤
 renorm/    │ Υ-map (nonlinearity→F(τ*)) · RenormalizationScheme/Character ·  │
            │ RenormalizedEquation assembler                                  │
            ├───────────────────────────────────────────────────────────────┤
 structures/│ RegularityStructure (T,T⁺,Δ,Δ⁺) · RenormalizationStructure     │  ← Phase 3
            │ (U,U⁻,δ,δ⁻) · twisted antipode S'₋ · BHZ character (symbolic)   │
            ├───────────────────────────────────────────────────────────────┤
 equation/  │ DSL (Unknown·Noise·Operator·SPDE) · IntegrationOperator · Rule  │
            │ · subcriticality · tree Generator                               │
            ├───────────────────────────────────────────────────────────────┤
 trees/     │ DecoratedTree(:Symbol) · canonical form · S(τ) · grafting/⋆/    │
            │ I_p/contraction/cuts · coproduct builders (Δ,Δ⁻ as data)        │
            ├───────────────────────────────────────────────────────────────┤
 core/      │ Homogeneity (ordered ring) · Scaling · Signature · FreeModule/  │
            │ Tensor · Symbol protocol · generic Hopf machinery               │
            └───────────────────────────────────────────────────────────────┘
   seams (out of the library):  NoiseLaw (probabilistic, deferred) ·
   analytic model (never implemented here) · symbolic backend (SymPy, isolated in one adapter)
```

---

## 3. Core abstractions (the seams), with interfaces

Each is given its responsibility, an interface sketch, the v1 concrete implementation, and the
named extension(s) that justify it as a seam (so no abstraction is speculative).

### 3.1 `Homogeneity` and `Scaling`  — `core/homogeneity.py`, `core/scaling.py`

Homogeneities are exact elements of the ordered ring `ℚ ⊕ (⊕_j ℚ·β_j) ⊕ ℚ·κ` — a rational part,
an integer coefficient per noise regularity `β_j`, and a κ-coefficient (κ>0 infinitesimal).
Comparison resolves the `β_j` to their user-given (rational-part, κ-part) values, then compares
lexicographically: rational part first, κ-coefficient as the infinitesimal tiebreak. **Never
floats** (critical trees sit at homogeneity `−kκ`).

```python
@dataclass(frozen=True)
class Homogeneity:
    rational: Fraction
    noise_coeffs: Mapping[NoiseId, int]   # Σ_j m_j β_j
    kappa: Fraction = 0
    def __add__(self, other): ...
    def resolve(self, sig: "Signature") -> tuple[Fraction, Fraction]:  # (standard, κ) parts
    def is_negative(self, sig) -> bool: ...
    def __lt__/__eq__ ...                  # via resolve(); total order

@dataclass(frozen=True)
class Scaling:
    weights: tuple[int, ...]               # 𝔰, e.g. (2,1,...,1)
    def scaled(self, k: MultiIndex) -> int # |k|_𝔰
```
*Seam justification:* `Scaling` is varied the moment a non-parabolic operator appears; the
multi-noise `noise_coeffs` map is needed for systems with several noises.

### 3.2 `Signature` — `core/signature.py`  ★ the linchpin of generality

The typed vocabulary the whole library is parametric over. Encodes *every* generality axis as
data. Trees, rules, homogeneities and structures are all defined *over a Signature*.

```python
@dataclass(frozen=True)
class Signature:
    dim: int                                   # spatial dimension d (spacetime = 1+d)
    scaling: Scaling
    components: tuple[Component, ...]           # the unknowns/sectors (scalar ⇒ length 1)
    operators: tuple[IntegrationOperator, ...]  # 𝔗_e: planting/integration ops, one+ per component
    noises: tuple[Noise, ...]                   # 𝔗_n noise types, each with regularity β_j
    # node types are derived: {● per component} ∪ {∘_j per noise} ∪ {red}
```
- **scalar vs system** ⇒ `len(components)`; component identity rides on the **edge type**
  (operator), never the node type.
- **single vs multiple noises** ⇒ `len(noises)`.
- **operator order / scaling** ⇒ the `operators` carry their own order and scaling.

Adding a feature = building a richer `Signature`; the algorithms do not branch on these.

### 3.3 `Symbol` protocol  — `core/symbol.py`

The atom of the graded free module. Deliberately minimal: the *generic* algebra needs only a
grading and a canonical identity. Everything structural (grafting, contraction, coproducts) is
implemented by the concrete basis and handed to the algebra as maps.

```python
class Symbol(Protocol):
    def homogeneity(self, sig: Signature) -> Homogeneity: ...
    def canonical_key(self) -> Hashable: ...     # equality / hashing / dedup
```
*Seam justification:* a second basis (multi-indices) is a documented goal; the algebra core (3.5)
is written purely against this protocol so that basis swap needs no algebra rewrite.

### 3.4 `FreeModule` / `LinearCombination` / `Tensor` — `core/module.py`

Exact free `ℚ`- (or `ℚ[constants]`-) module over `Symbol` keys; the substrate for every algebraic
object. Tree combinatorics live here, **not** in SymPy (SymPy auto-flattening makes large
tree-sums slow/awkward — convert to SymPy only at the leaves of the Υ-map and at render time).

```python
class LinearCombination(Generic[S]):   # dict {canonical_key: (Symbol, coeff∈Fraction)}
    def __add__, __mul__(scalar), items(), map(f), ...
class Tensor:                          # elements of M ⊗ M (and M ⊗ M ⊗ M)
class LinearMap:                       # S -> LinearCombination, extended linearly
```

### 3.5 Generic Hopf machinery — `core/hopf.py`  (Phase 3)

The basis-agnostic heart. Given a `Coproduct` (a `LinearMap` to a `Tensor`) it provides everything
the renormalization theory needs, knowing nothing about trees:

```python
class Coproduct(LinearMap): ...                       # Δ, Δ⁺, Δ⁻, δ — all just this type
class Character:                                       # algebra morphism Symbol-span → R (or symbols)
    def convolve(self, other, coproduct) -> Character  # (k₁*k₂) = (k₁⊗k₂)∘Δ
    def inverse(self, antipode) -> Character
def antipode(coproduct, grading) -> LinearMap          # connected-graded recursion
def comodule_action(character, coaction) -> LinearMap  # k̃ = (k⊗Id)∘δ
def twisted_negative_antipode(delta, P_minus) -> LinearMap   # S'₋ recursion (BPHZ)
```
*Seam justification:* this is exactly the code that is identical for trees and multi-indices, and
that is reused for all four coproducts. Writing it once, generically, is the payoff of the `Symbol`
seam.

### 3.6 `DecoratedTree` and tree operations — `trees/`

The one concrete `Symbol`. `tree.py`: the datatype (node type+decoration, edges with
type+decoration, children as a multiset), **canonicalization** (the load-bearing correctness
concern), and `symmetry_factor` `S(τ)`. `operations.py`: the multi-pre-Lie grafting, the
`⋆`-product/joining, `I_p` planting, node relabel, contraction `τ/φ`, cuts. `coproducts.py`
(Phase 3): functions that *build* `Δ, Δ⁺, Δ⁻, δ` as `Coproduct` objects over trees from the
extraction/recentering formulas.

```python
@dataclass(frozen=True)
class DecoratedTree:                 # τ = b^n ⋆ ⨉ᵢ I_{p_i}(τ_i)
    node_type: NodeType              # ● / ∘_j / red●^{0,α}
    node_dec: MultiIndex             # n
    children: tuple[Edge, ...]       # each Edge = (operator∈𝔗_e, edge_dec p, subtree), canonically sorted
    # implements Symbol: homogeneity() per §3.1 rules; canonical_key() = the canonical form
```

### 3.7 `IntegrationOperator` / `Kernel` — `equation/operator.py`

A planting/integration operator with its homogeneity gain and identity. Replaces the hardcoded
`+2`.

```python
@dataclass(frozen=True)
class IntegrationOperator:
    label: OperatorId            # ties an edge to a component/equation (𝔗_e)
    order: int                   # Schauder gain m (=2 for the heat kernel)
    scaling: Scaling
    def gain(self, p: MultiIndex) -> int:   # m − |p|_𝔰  ⇒  |I_p τ| = |τ| + gain(p)
```

### 3.8 `Rule` and `Generator` — `equation/rule.py`, `equation/generate.py`

`Rule`: node type → allowed child-edge multisets; `subcritical()`, `complete()` (BHZ closure),
`conforming(tree)`. Built by `rule_from_equation(spde)` *or* supplied explicitly (escape hatch
when auto-derivation can't be certified normal+complete). `Generator`: saturates strongly-
conforming trees up to a homogeneity cutoff (terminates iff subcritical) ⇒ `𝓑_{<0}` (including the
bare primitives `∘^n, ●^n`).

### 3.9 Nonlinearity / the Υ-map — `renorm/nonlinearity.py`

The bridge from trees back to PDE terms, and the only place SymPy enters the math. Holds the base
map `node_type → F(b*)` (`F(∘_j*)=f_j`, `F(●*)=g`, `F(red*)=0`) and the recursion
`F(τ*) = (Πᵢ F(τ_i*))·(D^n Πᵢ ∂_{p_i})F(b*)`, with `∂_p = ∂/∂u_p`, `D_i = Σ_k u_{k+e_i}∂_k`,
`∂_{p_i}` applied before `D^n`. Differentiates **all** slots of `g`. SymPy is isolated behind a
thin `backend` adapter so it could later be swapped (SymEngine/Symbolics) for speed — but not
abstracted prematurely.

### 3.10 `RenormalizationScheme` / `Character` / `NoiseLaw` — `renorm/scheme.py`

```python
class RenormalizationScheme(Protocol):
    def character(self, structure) -> Character: ...
class FreeConstants(RenormalizationScheme):    # default: c_τ are free SymPy symbols
class BPHZ(RenormalizationScheme):             # Phase 4: k^ζ = h^ζ ∘ S'₋ ; needs a NoiseLaw
class NoiseLaw(Protocol):                      # DEFERRED SEAM — not implemented in v1
    def expectation(self, symbol) -> Value: ...   # Gaussian/Wick; powers canonical constant values
```
*Seam justification:* the family is *by definition* the orbit under choices of character; making
`Character`/`Scheme` first-class is intrinsic, not speculative. `NoiseLaw` is named but unbuilt —
the one place we leave a boundary without an implementation, by explicit decision.

### 3.11 Output and rendering — `renorm/equation.py`, `render/`

`RenormalizedEquation` (base PDE + per-tree `Counterterm{tree, homogeneity, S(τ), constant c_τ,
vector_field F(τ*)}`, `.as_sympy()`, `.constants`). `render/` provides `Renderer` implementations
(LaTeX grouped by homogeneity à la the paper tables, unicode/pretty, SymPy). `structures/` can
**export** a built structure (trees, homogeneities, coproducts) for an external analytic/numerical
consumer — the clean handoff at the symbolic boundary.

---

## 4. Build pipeline (data flow)

```
SPDE (DSL)
  │  parse + validate (scope rejection list)              [equation/dsl]
  ▼
Signature + Nonlinearity                                  [core/signature, renorm/nonlinearity]
  │  rule_from_equation  → subcriticality → completion    [equation/rule]
  ▼
Rule ── Generator ──► 𝓑_{<0} (decorated trees)            [equation/generate, trees/]
  │                         │  homogeneity, S(τ)          [core, trees]
  │                         ▼
  │                    Υ-map F(τ*)                         [renorm/nonlinearity]
  │  (Phase 3) build Δ,Δ⁻,structures, S'₋, BHZ char       [trees/coproducts, structures/, core/hopf]
  ▼                         ▼
RenormalizationScheme ─► Character  ─► constants c_τ       [renorm/scheme]
  ▼
RenormalizedEquation ──► Renderer ──► LaTeX / SymPy        [renorm/equation, render/]
```

The Phase-3 structure construction is a *side branch*: the renormalized-equation pipeline reaches
the output via `FreeConstants` without ever building a coproduct. The coproducts/structures are
needed only for the canonical (BPHZ) character — exactly the deferred path.

---

## 5. Module layout

```
counterterms/
  core/
    homogeneity.py     scaling.py     signature.py
    symbol.py          module.py      hopf.py            # hopf.py: Phase 3
  trees/
    tree.py            operations.py  coproducts.py      # coproducts.py: Phase 3
  equation/
    dsl.py             operator.py    rule.py            generate.py
  structures/                                            # Phase 3
    regularity.py      renormalization.py
  renorm/
    nonlinearity.py    scheme.py      equation.py        backend.py   # backend.py: SymPy adapter
  render/
    latex.py           pretty.py
  api.py
tests/
  test_golden_gkpz.py  test_symmetry.py  test_homogeneity.py
  test_coproducts.py   # Phase 3: coassociativity/cointeraction as invariants
```

---

## 6. Extension cookbook (the proof that it's extensible)

Each is "edit/implement *only* the named place; nothing else changes":

- **New differential operator / scaling** → add an `IntegrationOperator(order=m, scaling=𝔰)` to the
  `Signature`. Homogeneity uses `m − |p|_𝔰` automatically. (Validation warns outside the proven
  parabolic regime.)
- **Systems / vector unknowns** → add a `Component` + its planting operator `I^{(a)}` per equation
  to the `Signature`; nonlinearity becomes a tuple `(f_a, g_a)`. Tree/rule/algebra code is
  unchanged (it is Signature-parametric); component identity rides the edge type.
- **Multiple / multi-component noises** → add `Noise(β_j)` entries (and derived `∘_j` node types)
  to the `Signature`; the Υ base map gains `f_j`.
- **A new symbol basis (multi-indices)** → implement the `Symbol` protocol + the coproduct
  builders for the new basis under a `multiindex/` package. `core/` (module, hopf, characters) and
  `renorm/scheme` work unchanged.
- **A new renormalization scheme** (e.g. Bruned preparation maps) → implement
  `RenormalizationScheme.character(structure)`.
- **Canonical BPHZ values** → implement a `NoiseLaw` (Gaussian Wick contractions) and the `BPHZ`
  scheme `k^ζ = h^ζ∘S'₋`. Everything above the seam is already in place.
- **A new renderer / export target** → implement a `Renderer`.
- **Analytic / numerical handoff** → consume `structures/`'s export; the library deliberately does
  not cross into analysis.

---

## 7. Phasing (boundaries fixed now; no empty stubs — each module is populated at first real use)

- **Phase 1 — equation pipeline.** `core/{homogeneity,scaling,signature,symbol,module}`,
  `trees/{tree,operations}`, `equation/*`, `renorm/{nonlinearity,scheme(FreeConstants),equation,
  backend}`, `render/`, `api.renormalize`. Delivers SPDE → family of renormalized equations with
  free constants, golden-tested on gKPZ. (This is the plan's MVP, now seated in the general
  skeleton.)
- **Phase 2 — generality.** Systems, multiple noises, general operators/scaling — almost entirely
  `Signature` enrichment + validation; little new structure, by design.
- **Phase 3 — Hopf core.** `core/hopf`, `trees/coproducts`, `structures/*`, twisted antipode, BHZ
  character (symbolic). Coassociativity/cointeraction become test invariants. Makes it a genuine
  regularity-structures library.
- **Phase 4 — deferred seams (optional).** `NoiseLaw` + BPHZ numerical values; the `multiindex/`
  basis; richer analytic export.

---

## 8. Relationship to `initial_plan.md`

The plan remains authoritative for **the mathematics** (pipeline §3, conventions §10, golden tests
§9, scope §8). This document refines the plan's §6 (reuse) and §7 (architecture) into the general,
modular structure above. Reuse decision is unchanged: **SymPy** (isolated in `renorm/backend.py`)
for the Υ-map and rendering; `kauri`/SageMath/`BSeries.jl` as correctness oracles, not runtime
dependencies. Stack: Python + SymPy, pip/`uv`-installable.
