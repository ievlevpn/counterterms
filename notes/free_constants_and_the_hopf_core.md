# What the Hopf core is for, and in what sense `k(τ)` is "free"

Discussion notes: why Phases 1–2 emit *free* constants `k(τ)` without ever building the
coproducts, what the Hopf algebra actually buys, and how "free" coexists with `k(τ)` being a
multiplicative character.

---

## 1. Why we got the renormalized *equations* without the Hopf core

BCCH (`ThmRenormPDEs`) gives a **closed formula** for the family of renormalized equations:
which trees appear (`|τ|<0`), their coefficient shape `k(τ)/S(τ)`, and the nonlinearity
`F(τ*)`. That's pure tree combinatorics + the Υ-map. No coproduct touches it **as long as the
constants `k(τ)` stay free symbols**.

The Hopf core is needed only for what lies *beyond* that symbolic family. It is two Hopf
algebras doing two different jobs:

**Positive coproduct `Δ⁺` → structure group `G⁺` (the analytic side).**
"Recentering": how a Taylor-like jet expansion around `x` re-expands around `y`. Defines the
structure group, hence what a *model* and *modelled distributions* are, and the fixed-point
map that actually *solves* the PDE. We never built it because we only emit the equation,
never solve it or touch a model.

**Negative coproduct `Δ⁻` + twisted antipode `𝒜⁻` → renormalization group `G⁻` (the
probabilistic side).** Pins the *numerical values* of `k(τ)`: `Δ⁻` extracts divergent
subtrees, the twisted antipode builds the counterterm with correct signs (BPHZ/BHZ
character), and you contract against the **noise covariance** (Wick) to get the actual
divergent constants. Exactly the "canonical constant values" that are out of project scope —
they need probabilistic input. We sidestep by leaving `k(τ)` free.

| What you want | Needs Hopf? |
|---|---|
| Family of renormalized equations, free constants | **No** — BCCH shortcut (Phases 1–2, done) |
| Define/solve the SPDE (models, Γ, fixed point) | **Yes** — `Δ⁺` / `G⁺` |
| Compute *actual* constant values from noise stats | **Yes** — `Δ⁻` + `𝒜⁻` / `G⁻` |

Deep point: BCCH's formula is itself *derived* from `Δ⁻` + the twisted antipode — the Hopf
algebra is upstream of everything. The project consumes the *theorem's output* (a closed
combinatorial formula) instead of re-deriving it. Phase 3 (`core/hopf`, `trees/coproducts`,
twisted antipode, BHZ character) is where it comes back.

---

## 2. "Free" constants — aren't they constrained by the coproducts?

The constraint runs the *other* direction from what one might guess.

**At the level of the family of equations, the constants genuinely are free.**
`G⁻` is the **character group of `H⁻`**, and `H⁻` is the *free* commutative algebra generated
by the connected negative-homogeneity trees. A character is an algebra morphism `H⁻ → ℝ`,
determined by — and free on — its values on the generators = the connected trees = exactly
the counterterms in the equation. So:

- every assignment `k(τ) ∈ ℝ` (one per connected `τ∈𝓑⁻`) gives a **valid** renormalized
  equation, and
- every admissible renormalization arises this way.

There is **no algebraic relation among the `k(τ)`** that the coproduct imposes (a theorem,
resting on the Hopf structure: the renormalized model is still a model). Its content is "the
freedom is total," not "the freedom is constrained."

**What the coproduct actually encodes — not relations on the values, but three other things:**

| Hopf object | What it does | Constrains `k(τ)` values? |
|---|---|---|
| `Δ⁻` (negative coproduct) | the **group law** on `G⁻` — how to *compose* two renormalizations (convolution) | No — multiplication on the group, not a condition on one element |
| `M = (ℓ⊗id)Δ⁻` | how a chosen character `ℓ` **acts on the model**, consistently (maps models→models) | No — consistency of the *action*, automatic for any `ℓ` |
| `Δ⁻` + twisted antipode + noise covariance | computes the **one canonical (BPHZ) point** `ℓ = g` | This is where "not any value" lives |

So "free" should read "**free up to the canonical choice**":
- *Free* family: `k(τ)` independent real parameters. No coproduct needed — why Phases 1–2 work.
- *Canonical* values: one distinguished point `g(τ)` per tree, fixed by noise statistics,
  computed by `Δ⁻` + twisted antipode + Wick. Not a relation tying the `k`'s together — a
  specific number per tree.

Two real caveats, **neither from the coproduct**:
1. *Which* trees carry a constant is constrained — by homogeneity (`|τ|<0`),
   rule-conformance, and `Υ(τ)≠0` (e.g. `F(red*)=0` drops some). Upstream combinatorics,
   already handled.
2. Multiplicativity `k(τ₁·τ₂)=k(τ₁)k(τ₂)` *is* a Hopf constraint — but only on **forests**.
   The equation is indexed by **connected** trees only, so it never bites the listed
   constants. (Expanded below.)

---

## 3. But `k(τ)` IS a character — isn't it multiplicative?

Yes — `k = ℓ` is a character and **is** multiplicative. Multiplicativity and "free" are not in
tension; they are the *same fact* from two sides, because of **where** multiplicativity bites.

`H⁻` is the **free commutative algebra** on the connected trees. Its product is the forest
product (disjoint union); its monomials are forests; its generators are the **connected**
trees. For a free/polynomial algebra:

> A multiplicative map out of a free commutative algebra is **uniquely and freely determined
> by arbitrary values on the generators.**

That is the definition of "free." Hence:

- Multiplicativity is the constraint `ℓ(τ₁·τ₂) = ℓ(τ₁)·ℓ(τ₂)` — relating the value on a
  **forest** (a product) to the values on its connected pieces.
- It says **nothing** about the value on a single **connected** tree, which is a *generator* —
  irreducible under the forest product, not a nontrivial product.

And the equation's counterterms are indexed by **connected** trees only. Forests never appear
as counterterms. So every constraint multiplicativity imposes lands on objects (forests) the
equation never displays.

```
ℓ(  τ₁τ₂  ) = ℓ(τ₁)ℓ(τ₂)     ← multiplicativity, on a FOREST (never a counterterm)
ℓ(   τ    ) = free            ← τ connected = generator = a counterterm
```

"Is it a character?" → yes. "Then isn't it constrained?" → only on forests, which don't show
up. The character being multiplicative is *precisely the mechanism* by which "specify free
values on the connected negative trees and you've specified the entire renormalization" holds.
Free-on-generators **is** the character.

**Where multiplicativity actually does work:** inside the machinery, not the parametrization.
`Δ⁻` extracts *subforests*, and applying `ℓ` to an extracted forest uses
`ℓ(forest)=∏ ℓ(connected pieces)`. So multiplicativity is load-bearing for the action on the
model `M=(ℓ⊗id)Δ⁻` and for the canonical values — but as a *parametrization of the family of
equations*, the independent data is still exactly `{ℓ(τ) : τ connected, |τ|<0}`.

Two normalizations ride along for free (not relations among values): `ℓ(𝟙)=1` (counit) and
`ℓ` supported on `|τ|<0`. Mild support conditions, handled upstream by the rule/homogeneity,
not the coproduct.

---

## TL;DR

- The renormalized **equations** (free constants) need no Hopf core — BCCH hands us a closed
  formula. That's Phases 1–2.
- The Hopf core is for (a) the **structure group `G⁺`** — defining/solving the PDE via models,
  and (b) the **renormalization group `G⁻` + twisted antipode** — fixing the *canonical*
  constant values from noise statistics. That's Phase 3.
- `k(τ)` is a multiplicative character, yet free, because `H⁻` is free on the connected trees
  and the equation only ever evaluates the character on those generators. Multiplicativity
  constrains only forests, which are never counterterms.
- "Free" = free up to the single canonical (BPHZ) point that `Δ⁻` + twisted antipode + Wick
  selects.
