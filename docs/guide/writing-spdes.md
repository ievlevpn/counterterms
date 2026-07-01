# Writing an SPDE

The input DSL is four classes — `Unknown`, `Noise`, an `Operator`, and `SPDE` — plus plain SymPy
for the right-hand side.

## The pieces

### Unknown

```python
u = Unknown("u", dim=1)
```

A solution component in `dim` space dimensions. It exposes the SymPy objects you build the
right-hand side from:

- `u.field` — the function \(u(t, x_1, \dots)\);
- `u.t`, `u.x` — the time symbol and the list of space symbols;
- write \(\partial_{x_1} u\) as `Derivative(u.field, u.x[0])`.

All components of a system must share the same `dim`.

### Noise

```python
from counterterms import kappa
xi = Noise("xi", regularity=Rational(-3, 2) - kappa)
```

`regularity` is the noise’s Hölder regularity \(\beta_0\) **directly**
(see [Conventions §1](../conventions.md#1-xi-beta_0-directly)). It may be a rational, or a
rational minus a multiple of `kappa` (the positive infinitesimal); white noises always carry a
`- kappa`. Use `xi.symbol` in the right-hand side.

### Operator

```python
Parabolic(dim=1)                          # ∂_t − Δ            (the default choice)
Parabolic(dim=1, mass=1)                  # ∂_t − Δ + 1
Parabolic(dim=2, order=4)                 # order-4 parabolic — warns, see below
FractionalHeat(dim=1, sigma=Rational(3,4))  # ∂_t + (−Δ)^σ     — warns, see below
```

The engine reads **exactly two numbers** from the operator: the metric scaling \(\mathfrak{s}\)
and the Schauder order \(m\) of \(L^{-1}\) (it maps \(\mathcal{C}^\alpha \to
\mathcal{C}^{\alpha+m}\)). Everything downstream — homogeneities, \(S(\tau)\), the Υ-map — is
computed from that pair; the display strings are cosmetic. You can subclass `Operator` to model
any linear operator whose inverse has a Schauder estimate in a scaled Hölder scale.

!!! warning "Order ≠ 2 runs with a warning"
    The regularity-structure theory behind the output is proven for 2nd-order parabolic
    operators. For other orders the engine computes the combinatorics anyway (the tree
    enumeration is complete for any order) and warns that the analytic theory is unverified —
    see [Scope & limitations](../scope.md#operator-order-2). Treat order ≠ 2 output as
    exploratory.

### SPDE

Scalar form:

```python
spde = SPDE(operator=op, unknown=u, noises=[xi], rhs=...)
```

System form — a list of `(unknown, operator, rhs)` triples:

```python
spde = SPDE(noises=[xi], equations=[
    (u, op, a(v.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0])**2),
    (v, op, b(u.field) * xi.symbol),
])
```

Then:

```python
eq = spde.renormalize()      # -> RenormalizedEquation
```

## The right-hand side

Plain SymPy, with three kinds of atoms: noise symbols (`xi.symbol`), the fields (`u.field`,
`v.field`), and their first spatial derivatives (`Derivative(u.field, u.x[i])`). Coefficient
functions are undeclared SymPy `Function`s of the fields — `f(u.field)`,
`g(u.field, v.field)`, etc.

What the parser accepts (and the theory covers):

- **affine in the noise**: `f(u)·ξ + (noise-free part)` — a `ξ²` or `f(ξ)` is rejected;
- **at most quadratic in \(\partial u\)** (Assumption D2), counting *total* gradient degree —
  `g(u)·(∂ₓu)²` is fine, `(∂ₓu)³` is rejected;
- polynomial or symbolic-function dependence on the fields themselves is unrestricted.

## Systems and multiple noises

```python
xi  = Noise("xi",  regularity=Rational(-1) - kappa)
eta = Noise("eta", regularity=Rational(-1) - kappa)
res = SPDE(operator=op, unknown=u, noises=[xi, eta],
           rhs=f(u.field) * xi.symbol + h(u.field) * eta.symbol).renormalize()
```

- Noises are treated as **independent**; with several noises, subcriticality binds through the
  *worst* regularity, and Wick pairing (in the canonical constants) is within-type only.
- In systems, one tree can produce counterterms in several equations; the free constant is
  **shared** — `res.per_component[a]` lists each equation’s counterterms, and the same
  `k_i` symbol may appear in several components (that is the theorem, not a bug).
- `res.n_components`, `res.per_component[comp]`, `res.counterterm_rhs(comp)` give per-equation
  access.

## What gets rejected, and why

The engine enforces its scope with explicit errors rather than wrong output:

| Input | Response |
|---|---|
| Supercritical \(\beta_0 \le -\text{order}\), polynomial additive noise (\(\Phi^4_2, \Phi^4_3\)) | `ValueError` pointing to [`daprato_lift`](daprato.md) |
| Supercritical, not liftable (sine-Gordon — needs Wick exponentials) | rejected |
| Nonlinearity not affine in the noise | rejected |
| \(g\) more than quadratic in \(\partial u\) | rejected (Assumption D2) |
| Derivative factors with \(|p|_\mathfrak{s} > 1\) (e.g. \(\partial_t u\) or \(\partial_x^2 u\) in the nonlinearity) | rejected |
| Quasilinear / non-parabolic \(L\) | not expressible in the DSL |
