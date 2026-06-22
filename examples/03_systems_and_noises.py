"""Systems (coupled, shared constants) and multiple noises.
Run:  uv run python -u examples/03_systems_and_noises.py
"""
from sympy import Derivative, Function, Rational

from counterterms import SPDE, Noise, Parabolic, Unknown, kappa

a, b, g = Function("a"), Function("b"), Function("g")
u, v = Unknown("u", 1), Unknown("v", 1)
xi = Noise("xi", regularity=Rational(-1) - kappa)
op = Parabolic(dim=1, mass=1)

# coupled system sharing one noise:  u̇ = a(v)ξ + g(u)(∂u)²,   v̇ = b(u)ξ
res = SPDE(noises=[xi], equations=[
    (u, op, a(v.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2),
    (v, op, b(u.field) * xi.symbol),
]).renormalize()
print(f"coupled system: {res.n_components} components")
for comp in range(res.n_components):
    print(f"  equation {comp}: {len(res.per_component[comp])} counterterms")
