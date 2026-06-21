# regstruct

Symbolic **renormalized equations for singular SPDEs**, following Bailleul & Hoshino,
*"A tourist's guide to regularity structures and singular stochastic PDEs"*
([arXiv:2006.03524](https://arxiv.org/abs/2006.03524)).

You give it a subcritical singular SPDE; it gives you the **family of renormalized equations**
(the BCCH / `ThmRenormPDEs` formula): the original PDE plus one tree-indexed counterterm per
negative-homogeneity decorated tree, with free renormalization constants.

```
(∂_t − Δ + 1) u⁽ᵏ⁾ = f(u⁽ᵏ⁾)ζ + g(u⁽ᵏ⁾,∂u⁽ᵏ⁾) + Σ_{τ∈𝓑, |τ|<0} (k_τ / S(τ)) F(τ*)(u⁽ᵏ⁾,∂u⁽ᵏ⁾)
```

The free constants `k_τ` parametrise the renormalization group `G⁻_ad`. Their *canonical* (BPHZ)
numerical values need probabilistic input (Gaussian/Wick expectations) and are out of scope — the
library emits the constants symbolically.

## Install

Uses [`uv`](https://docs.astral.sh/uv/) (SymPy is the only runtime dependency).

```sh
uv sync
uv run pytest        # 7 tests
```

## Example — the generalized KPZ equation

`(∂_t − Δ + 1)u = f(u)ζ + g(u)(∂_x u)²` in `d = 1` with `ζ ∈ C^{−1−κ}`
(tourist_guide.tex 5996–6012):

```python
from sympy import Function, Derivative, Rational
from regstruct import Unknown, Noise, Parabolic, SPDE, kappa

u  = Unknown("u", dim=1)
xi = Noise("xi", regularity=Rational(-1) - kappa)     # ζ ∈ C^{−1−κ}
f, g = Function("f"), Function("g")

spde = SPDE(
    operator = Parabolic(dim=1, mass=1),              # ∂_t − Δ + 1
    unknown  = u,
    noises   = [xi],
    rhs      = f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0])**2,
)

print(spde.renormalize().summary())
```

```
Renormalized equation: 5 counterterm(s)
  |τ|=  -1 - κ  S=1  + k_0*f(u(t, x1))
  |τ|=     -2κ  S=1  + k_2*f(u(t, x1))*Derivative(f(u(t, x1)), u(t, x1))
  |τ|=     -2κ  S=2  + k_4*f(u(t, x1))**2*g(u(t, x1))
  |τ|=      -κ  S=1  + k_1*Derivative(f(u(t, x1)), u(t, x1))*Derivative(u(t, x1), x1)
  |τ|=      -κ  S=1  + 2*k_3*f(u(t, x1))*g(u(t, x1))*Derivative(u(t, x1), x1)
```

These are exactly the five counterterms of the paper. `spde.renormalize()` returns a
`RenormalizedEquation`; each entry of `.counterterms` carries the tree, its homogeneity `|τ|`,
the symmetry factor `S(τ)`, the elementary differential `F(τ*)`, and the free constant `k_τ`.

## Scope

**In scope (Phase 1):** scalar unknown, single noise, second-order parabolic operator,
`β₀ ∈ (−2,−1)`, nonlinearity `f(u)ζ + g(u,∂u)` with `g` at most quadratic in `∂u`
(Assumption D2), `|p|_𝔰 ≤ 1`. Covers **gKPZ, KPZ, gPAM, PAM**.

**Rejected with clear errors:** `β₀ ≤ −2` (Φ⁴₂/Φ⁴₃/sine-Gordon — need a da Prato–Debussche
pre-pass), noise nonlinearities not affine in the noise, `g` more than quadratic in `∂u`,
singular derivative factors with `|p|_𝔰 > 1`, quasilinear / non-parabolic operators.

## Status

Phase 1 (`SPDE → renormalized equations`, free constants) is implemented and golden-tested.
Phase 2 (systems, multiple noises, general operators) and Phase 3 (coproducts `Δ/Δ⁻`, the
regularity/renormalization structures, the BHZ character) are designed but not yet built.
See [`notes/architecture.md`](notes/architecture.md) §7 for the phasing and
[`CHANGELOG.md`](CHANGELOG.md).

## Documentation

- [`notes/initial_plan.md`](notes/initial_plan.md) — the mathematics: pipeline, conventions,
  scope, golden tests, and the reuse analysis.
- [`notes/architecture.md`](notes/architecture.md) — the module structure: layered stack, the
  `Signature`-parametric design, interfaces, and the extension cookbook.
- [`references/`](references/) — the source paper (PDF + full LaTeX). Cited by tex line number.

## Validation

The paper is the only oracle (no reference implementation exists). The test suite checks the
gKPZ example's exact five counterterms, the KPZ homogeneity table, gPAM in `d = 2`, the symmetry
factors (including the `S = 2` and factor-2 cases), and the scope rejections.
