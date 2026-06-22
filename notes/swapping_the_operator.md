# Swapping the linear operator `∂_t − Δ + 1`

Short version: the operator enters the renormalization machinery through **exactly three
numbers**, and the literal PDE — the `∂_t`, the `Δ`, the `+1` — never appears anywhere
downstream. Only the *scaling symbol* of the Green's function matters. So "swapping the
operator" = handing the pipeline a new operator object with three fields.

## What actually gets read

`build_context` (`equation/dsl.py:145`) only ever touches three attributes of
`spde.operator`:

| field | what it is | where it bites |
|---|---|---|
| `scaling` | the scaling 𝔰 = weights of (t, x₁,…,x_d), e.g. `(2,1,…,1)` | `|X^n| = |n|_𝔰`, `|p|_𝔰`, the `>1` derivative-scope check |
| `order` | the Schauder smoothing order *m* of `L⁻¹` | `|I_p τ| = |τ| + (m − |p|_𝔰)` (`signature.py:45`) |
| `label` | edge-type tag (`"I"`) | bookkeeping only |

The `mass=0` argument and the symbolic form of `L` are **decorative** — `Parabolic.__init__`
computes `scaling`, `order`, `label` and nothing reads the rest. That's the whole math
content of "the operator": homogeneities are set by the *principal scaling symbol* of `L⁻¹`,
and lower-order terms (the `+1`, drift, etc.) are subcritical perturbations that don't change
a single tree homogeneity.

## How to swap it

**Done.** `equation/dsl.py` now has a generic `Operator(dim, scaling, order, label, symbol,
latex)` — the "general way" — with `Parabolic` and `FractionalHeat` as one-line subclasses:

```python
class FractionalHeat(Operator):     # ∂_t + (−Δ)^σ
    def __init__(self, dim, sigma):
        order = 2 * sigma           # time weight = Schauder gain = 2σ
        super().__init__(dim, Scaling(tuple([order] + [1]*dim)), order, "I",
                         symbol=f"∂_t + (−Δ)^{sigma}", latex=rf"\partial_t + (-\Delta)^{{{sigma}}}")
```

The invariant you must preserve: **`order` = the Schauder gain of `L⁻¹` in the metric
`scaling`** (i.e. `L⁻¹` maps 𝒞^α → 𝒞^{α+order} in the 𝔰-scaled Hölder scale). Get that pair
right and every `S(τ)`, `Υ`-map, and counterterm comes out correct automatically. See
`examples/06_fractional_heat.py` — gPAM at β₀=−1−κ gives **3** counterterms under heat and **5**
under σ=¾ (order 3/2), with the homogeneities shifting by the order.

**Rendering:** the renderer prints `op.symbol` / `op.latex` (display-only — the engine never
reads them); without them it falls back to a generic `L`. So a new operator class must set those
two strings to print nicely.

## What you can / can't reach this way

- **Fractional / higher-order parabolic** (`∂_t + (−Δ)^σ`, `∂_t + Δ²`): yes — just `scaling`
  + `order`. For `Δ²`: `scaling=(4,1,…,1)`, `order=4`.
- **Purely elliptic** (no time, e.g. `−Δ + 1` on ℝ^d): the *math* is fine
  (`scaling=(1,…,1)`, `order=2`), but the **DSL hardcodes a time coord** — `width = 1 + u.dim`,
  `Parabolic` prepends a `2` to the scaling, and `_to_jet` assumes spacetime. So an elliptic
  op needs a small DSL change (drop the implicit time axis), not just a new operator class.
- **Anisotropic scaling** (different weights per axis): yes, fully data-driven via
  `Scaling(...)`.
- **Quasilinear / non-parabolic / `L` depending on `u`**: out of scope by design — and
  subcriticality must still hold (more negative `order` vs. noise regularity), which isn't
  auto-checked beyond the existing `β₀∈(−2,0)` and `|p|_𝔰≤1` guards.

## Not yet done

- The elliptic (no-time) case needs the DSL width change, not just an operator swap (~10-line
  change in `dsl.py`): drop the implicit time axis in `width`, `Parabolic`, `_to_jet`.
- ~~Generic `Operator`~~ — **done** (see above); `Parabolic`/`FractionalHeat` are subclasses.
- Caveat: fractional/high order **enlarges the negative-tree set fast**. gPAM stays small, but
  the gKPZ `(∂u)²` term at order 3/2 explodes enumeration (effectively non-terminating in
  practice). Subcriticality still holds; it's a combinatorial blow-up, not an error. Keep
  fractional examples to mild nonlinearities until generation is memoized.
