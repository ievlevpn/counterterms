# Examples — a quickstart guide

Runnable, self-contained scripts. Each prints labelled output:

```sh
uv run python -u examples/01_renormalized_equation.py
```

| # | Script | What it shows | Key API |
|---|--------|---------------|---------|
| 01 | `01_renormalized_equation.py` | **The headline use case** — an SPDE in, the BCCH counterterm family out (gKPZ). | `SPDE(...).renormalize()` |
| 02 | `02_trees_and_structure.py` | The divergent trees, the graded regularity structure `T`, and the recentering coproduct `Δ`. | `build_regularity_structure`, `render.shorthand` |
| 03 | `03_systems_and_noises.py` | Coupled **systems** (shared constants) and **multiple noises**. | `SPDE(equations=[...])` |
| 04 | `04_daprato_phi4.py` | The **da Prato–Debussche lift** — renormalize a supercritical Φ⁴₃ (β₀≤−order). | `daprato_lift` |
| 05 | `05_bhz_and_export.py` | The **BHZ character** (twisted antipode), the canonical-parity rule, and the JSON **export**. | `build_renormalization`, `render.structure_json` |

## The 60-second version

```python
from sympy import Derivative, Function, Rational
from regstruct import SPDE, Noise, Parabolic, Unknown, kappa

u  = Unknown("u", dim=1)
xi = Noise("xi", regularity=Rational(-1) - kappa)        # ζ ∈ C^{−1−κ}
f, g = Function("f"), Function("g")

spde = SPDE(operator=Parabolic(dim=1, mass=1), unknown=u, noises=[xi],
            rhs=f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0])**2)

eq = spde.renormalize()
print(eq.summary())          # the renormalized family
for ct in eq.counterterms:   # each: tree, |τ|, S(τ), F(τ*), free constant k_τ
    print(ct.homogeneity, ct.symmetry_factor, ct.elem_diff, ct.coefficient)
```

```
[u]  5 counterterm(s)
   |τ|=  -1 - κ  S=1  + k_0*f(u)
   |τ|=     -2κ  S=2  + k_3*f(u)**2*g(u)
   |τ|=     -2κ  S=1  + k_4*f(u)*Derivative(f(u), u)
   |τ|=      -κ  S=1  + k_1*u_x*Derivative(f(u), u)
   |τ|=      -κ  S=1  + 2*k_2*u_x*f(u)*g(u)
```

The free constants `k_τ` parametrise the renormalization group; their canonical (numeric)
values need probabilistic input (out of scope — see the README "What it does *not* do").

See the top-level [`README.md`](../README.md) for the full capability → module map.
