# Getting started

## Install

The project uses [`uv`](https://docs.astral.sh/uv/). SymPy is the only runtime dependency.

```sh
git clone <this repository>
cd regularity_structures_symbolic
uv sync
uv run pytest          # 187 tests, ~10 s — all green means the goldens reproduce
```

Always run Python through uv:

```sh
uv run python -u your_script.py
```

## First renormalization: generalized KPZ

The headline use case, start to finish:

```python
from sympy import Derivative, Function, Rational
from counterterms import SPDE, Noise, Parabolic, Unknown, kappa

u  = Unknown("u", dim=1)                          # one unknown, one space dimension
xi = Noise("xi", regularity=Rational(-1) - kappa) # noise in C^{-1-κ}
f, g = Function("f"), Function("g")

# generalized KPZ:  (∂_t − Δ + 1) u = f(u) ξ + g(u) (∂_x u)²
spde = SPDE(
    operator=Parabolic(dim=1, mass=1),
    unknown=u,
    noises=[xi],
    rhs=f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2,
)

eq = spde.renormalize()        # -> RenormalizedEquation
print(eq.summary())
```

The summary lists the five gKPZ counterterms — for each divergent tree \(\tau\): its
homogeneity \(|\tau|\), symmetry factor \(S(\tau)\), elementary differential \(F(\tau^*)\), and
free constant \(k_\tau\). This is exactly the family displayed in the source paper
(tex 6004–6012), reproduced term for term, including the asymmetric factor 2.

Programmatic access:

```python
for ct in eq.counterterms:
    print(ct.homogeneity, ct.symmetry_factor, ct.elem_diff, ct.coefficient)
    # ct.coefficient == ct.constant / ct.symmetry_factor  (the k_τ/S(τ) convention)
```

## Reading the noise regularity

`kappa` is the project’s positive infinitesimal: homogeneities live in the ordered ring
\(\mathbb{Q} \oplus \mathbb{Q}\kappa\), and `Rational(-1) - kappa` means
\(\beta_0 = -1 - \kappa\), i.e. “Hölder regularity just below \(-1\)”. Typical values:

| Noise | Regularity to pass |
|---|---|
| space white noise, \(d=1\) | `Rational(-1,2) - kappa` |
| space white noise, \(d=2\) | `Rational(-1) - kappa` |
| spacetime white noise, \(d=1\) | `Rational(-3,2) - kappa` |
| spacetime white noise, \(d=3\) | `Rational(-5,2) - kappa` (supercritical — use the [lift](guide/daprato.md)) |

The convention is direct: `regularity` **is** \(|\Xi| = \beta_0\), not “\(\alpha - 2\)” or any
other shifted convention. See [Conventions](conventions.md#1-xi-beta_0-directly) — this is the
single most dangerous off-by-two in the subject and the engine pins it with tests.

## Rendered reports

```python
eq.save("gkpz", outdir="output")     # text + markdown + json + LaTeX (and PDF if latexmk exists)
print(eq.report(canonical=True))     # adds the BPHZ constants k_τ = h(S'_- τ)
print(eq.report(reduced=True))       # …fully reduced for a symmetric noise
```

See [Reports & output](guide/reports.md) for the formats and the meaning of
`canonical` / `reduced` / `symmetric`.

## Runnable examples

Six commented scripts under `examples/` cover the whole surface; each is a
one-command demo:

```sh
uv run python -u examples/01_renormalized_equation.py   # the quickstart above
uv run python -u examples/02_trees_and_structure.py     # trees, (T,T⁺), recentering Δ
uv run python -u examples/03_systems_and_noises.py      # coupled systems, several noises
uv run python -u examples/04_daprato_phi4.py            # Φ⁴₃ via the da Prato–Debussche lift
uv run python -u examples/05_bhz_and_export.py          # twisted antipode, BHZ character, JSON export
uv run python -u examples/06_fractional_heat.py         # swapping the linear operator
```
