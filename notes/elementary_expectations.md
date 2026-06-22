# Elementary expectations `h(σ)` and which ones we can write as integrals

**Scope.** This note pins down, against the paper (arXiv:2006.03524 v3, cited by `.tex` line),
what the elementary expectations `h(σ) = 𝔼[Π^ζ σ](0)` in the canonical (BPHZ) character are, on
*which* trees `σ` they act, and exactly which of those we can render as an explicit divergent
integral vs. must leave symbolic. It documents the validity domain enforced by
`renorm.scheme.expectation` and the reasoning behind it. (Companion to `notes/phase4_plan.md` B1
and `ROADMAP.md` O3.5.)

## The object

The canonical (BPHZ) renormalization character is (tex 5053–5060)

```
h^ζ(τ) := 𝔼[Π^ζ τ](0)   (a character on ℝ[U]),     k^ζ := h^ζ ∘ S'₋ .
```

- `S'₋ : U⁻ → ℝ[U]` is the **negative twisted antipode** (tex 5028–5036): an algebra morphism
  into the polynomial algebra over the model space `U`. `S'₋(τ)` is a linear combination of
  **forests** (products) of trees in `U`.
- `h^ζ` is multiplicative, so `k^ζ(τ) = Σ coeff · ∏_{σ∈forest} h^ζ(σ)` — a polynomial in the
  elementary expectations `h^ζ(σ)`. This is exactly what `structures.canonical_character`
  returns symbolically (each `h^ζ(σ)` a free symbol `h_i`), parity-reduced.

So the `σ` that get an `h`-symbol are **whatever `S'₋` emits** — and that includes trees with
polynomial `X^n` node-decorations *and* red contraction nodes (see below). The symbolic character
is correct regardless; the only question is which `h(σ)` we can *write down as an integral*.

## What `Π^ζ` is — and the *one* thing that breaks the naive integral

The canonical model `Π^ζ` is defined **explicitly** (tex 4000–4008) — multiplicative for the
`⋆`-product, with

```
Π^ζ(∘^n)(x) = x^n ζ(x),    Π^ζ(●^n)(x) = Π^ζ(●^{n,α})(x) = x^n,    Π^ζ(I_p τ) = ∂^p K(Π^ζ τ).
```

For a **bare** tree — *every node decoration `n = 0`* — this collapses to the kernel⊗noise
integral, and `𝔼[·](0)` is exactly the Wick-pairing sum (Isserlis):

```
h(σ) = Σ_{pairings P} ∫ ∏_{edges} ∂^p K · ∏_{(i,j)∈P} C_ε(z_i − z_j)  d z   (root at 0).
```

**Red contraction nodes are *not* a problem** — this corrects an earlier wrong reading of this
note. `Π^ζ(●^{n,α})(x) = x^n` is **independent of the extended decoration `α`** (tex 4003–4004),
so a red node with `n = 0` is just `Π^ζ = 1`: the plain "1"-vertex that `expectation` already
builds (a vertex with no noise, connected by kernels). E.g. `h(●·𝓘ₓ[●]) = ∫ ∂^{(0,1)}K(−z) dz`
is the *correct* value (a pure-kernel total derivative, = 0 — convergent, not divergent). The
`𝔬`-decoration `𝔬 : N^red → ℤ[β₀]` (tex 5495–5502) is pure grading: it adds to the *extended*
homogeneity (tex 5502) but the *naive* homogeneity ignores it (tex 5506), only `𝔬=0` trees enter
the fixed point (tex 5512), and it has *"no dynamical meaning"* (tex 7409) — including for `Π^ζ`.

The **only** thing that breaks the naive integral is a **non-zero node decoration `n`**, on *any*
node (bullet, noise, or red), because `Π^ζ` then carries an `x^n` factor:

1. **Root `X^n` ⇒ `h = 0`.** `Π^ζ(X^n·τ)(x) = x^n Π^ζτ(x)`; at the base point `x=0` this is
   `0^n = 0` (tex 1809, and tex 5083 spells out `h^ζ(X^n⋆τ)=0`). E.g. `Xₓ·●·𝓘ₓ[Ξ]²` is truly **0**.
2. **Internal `X^n` ⇒ polynomial factor `z^n`** in the integrand — not captured by the naive
   kernel⊗noise build.

`expectation` refuses non-bare (decorated) trees for this reason; red nodes pass through and are
computed correctly.

## The validity domain (what `expectation` enforces)

`renorm.scheme.expectation(σ)` is the **complete** `𝔼[Π^ζσ](0)` **iff `σ` is bare**
(`is_bare`: **every node decoration is zero** — *red color is irrelevant*). Evaluated in order:

1. **Parity ⇒ 0.** Odd count of any noise type ⇒ `𝔼` of an odd number of centered Gaussians is
   `0`. Always valid: red nodes carry no noise (`Π^ζ(●^{0,α})=1`), and the deterministic
   kernel/`x^n` factors don't change parity. So e.g. `Xₓ·Ξ` returns `0` rather than raising.
2. **Non-zero node decoration ⇒ refuse** (`Π^ζ` carries `x^n`; value isn't the naive integral).
3. **Bare ⇒ the integral** — with red nodes built as plain `1`-vertices.

The render (`render/report.py`, `render/latex.py`) shows the explicit `h_ε(σ) = ∫ … d z` for bare
`σ` (red nodes included) and leaves a decorated `σ` as the symbol `h_ε(σ)` with a note. The
symbolic `k^ζ = h^ζ∘S'₋` is untouched and correct; only the *decorated* `h`'s lack an integral.

Note a no-noise bare tree (e.g. `●·𝓘ₓ[●]`) gets a **pure-kernel** integral with *no* covariance
(`∫ ∂^pK d z`) — correct (often a convergent total derivative, = 0), so "every integral has a
`C_ε`" is *not* an invariant any more.

Independent noises factor over noise types: vertices pair only with same-type vertices (no
cross-noise covariance, since `𝔼[ξη]=0`), and any type with odd count ⇒ 0. (Fix in the same pass:
the previous code paired all noises indiscriminately.)

## Are red nodes present? — Yes, and they are computed (not deferred)

Verified on the current code (distinct surviving `h`-arguments per equation), classified by the
*corrected* criterion (red is bare iff its decorations are 0):

| equation | bare (incl. red, get an integral) | `X^n`-decorated (left symbolic) |
|---|---|---|
| gKPZ @ β₀=−1 | 2 | 0 |
| gPAM d=2 | 1 | 0 |
| multinoise | 2 | 0 |
| **KPZ** | 13 | 2 |

So **KPZ** has red contraction nodes in the majority of its surviving `h`-arguments (e.g.
`●·𝓘ₓ[●]`), forced theoretically: `S'₋` keeps the *contracted* right-factors of `δ` (tex 5004,
5028–5034). **They get correct integrals** — `Π^ζ(●^{n,α})=x^n` is `α`-independent (tex 4003).
Only the two genuinely `X^n`-decorated survivors are left symbolic.

## The kernel `K` (a correctness note)

`K` in these integrals is the **singular kernel** of `L⁻¹` — Hairer's `K` in `K̄ = K + R`, the
diagonal-singular, compactly-supported part of the Green's function `K̄` (tex 2105). It *explodes
on the diagonal* (tex 5683), which is exactly why the `ε→0` limits diverge and renormalization is
needed. `K` is **not** the Green's function (that is `K̄ = K + R`). `∂^p K` appears on a derivative
edge `I_p`.

## Convergence (answering "isn't it too singular to converge?")

For the **ε-regularized** noise (`C_ε` smooth) each bare `h_ε(σ)` is *finite*. As `ε→0`,
`C_ε → δ` and `K` explodes on the diagonal, so the divergent `σ` (negative homogeneity — the
counterterm trees) diverge. The individual `h_ε(σ)` blow up; the BPHZ combination `k_τ = h(S'₋τ)`
is the finite object after subtraction. Evaluating that `ε→0` limit is the Track-B2 wall.

## Residual / follow-up

`canonical_character` still keeps `h`-symbols whose true value is `0` — both root-`X^n` trees
(by tex 1809) and no-noise total-derivative trees like `●·𝓘ₓ[●]` (`∫∂^pK = 0`) — exactly as it
pre-reduces odd-parity trees. This is **symbolically correct** (those `h_i` equal 0 downstream),
just not fully reduced. A "root-`X^n` ⇒ 0" reduction alongside the parity reduction (and, with the
explicit pure-kernel integrals, possibly recognising vanishing total derivatives) would tighten
the displayed constants — a clean follow-up, not a bug.

## Tests / anchors

- `tests/test_scheme.py::test_expectation_red_ok_refuses_only_decorations` — a red node with `n=0`
  is bare and **computed** (pure-kernel integral); only an `X^n` decoration is refused.
- `tests/test_scheme.py::test_independent_noises_no_cross_pairing` — the factor-over-noise-types fix.
- `tests/test_render.py::test_canonical_shows_epsilon_regularized_integrals` — the bare-tree
  integrals render with `ξ_ε`, `C_ε`, `∫`, `∂^p K`.
