"""Quickstart — derive the renormalized equation for a singular SPDE.

The headline use case: SPDE in → the BCCH counterterm family out.
Run:  uv run python -u examples/01_renormalized_equation.py
"""
from sympy import Derivative, Function, Rational

from counterterms import SPDE, Noise, Parabolic, Unknown, kappa

# generalized KPZ:  (∂_t − Δ + 1)u = f(u)ζ + g(u)(∂_x u)²,  d=1,  ζ ∈ C^{−1−κ}
u = Unknown("u", dim=1)
xi = Noise("xi", regularity=Rational(-1) - kappa)
f, g = Function("f"), Function("g")
spde = SPDE(operator=Parabolic(dim=1, mass=1), unknown=u, noises=[xi],
            rhs=f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2)

eq = spde.renormalize()                          # -> RenormalizedEquation
print(eq.summary())                              # the family of counterterms

print("\nProgrammatic access (each counterterm):")
for ct in eq.counterterms:
    print(f"  |τ|={ct.homogeneity}  S(τ)={ct.symmetry_factor}  "
          f"F(τ*)={ct.elem_diff}   coeff={ct.coefficient}")
