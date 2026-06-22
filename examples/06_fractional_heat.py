"""Swapping the linear operator — fractional heat ``∂_t + (−Δ)^σ``.

The operator enters the engine through exactly **(scaling, order)** — the literal
``∂_t``/``Δ`` never appear downstream (see notes/swapping_the_operator.md).
`FractionalHeat(dim, sigma)` sets both: time weight ``2σ`` and Schauder order ``2σ``.
Here we renormalize the *same* multiplicative equation under the heat operator and a
fractional one and watch the counterterms shift with the Schauder order.

(σ≠1 ⇒ order≠2, so the engine warns: the analytic theory is proven only for 2nd-order
parabolic L. The combinatorics are computed regardless.)

Run:  uv run python -u examples/06_fractional_heat.py
"""
from sympy import Function, Rational

from counterterms import SPDE, FractionalHeat, Noise, Parabolic, Unknown, kappa

# gPAM-type:  L u = g(u) ξ,   d=1,   ξ ∈ C^{−1−κ}
# (Heat → 3 counterterms; fractional σ=3/4 → 5, with shifted homogeneities.)
u = Unknown("u", dim=1)
xi = Noise("xi", regularity=Rational(-1) - kappa)
g = Function("g")
rhs = g(u.field) * xi.symbol

for op in [Parabolic(dim=1), FractionalHeat(dim=1, sigma=Rational(3, 4))]:
    spde = SPDE(operator=op, unknown=u, noises=[xi], rhs=rhs)
    eq = spde.renormalize()
    print(f"\n=== {op.symbol}   (order m={op.order}, scaling 𝔰={op.scaling.weights}) ===")
    print(eq.summary())
