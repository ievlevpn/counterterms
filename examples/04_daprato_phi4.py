"""da Prato–Debussche lift — renormalize a supercritical Φ⁴ equation (β₀ ≤ −order).
Run:  uv run python -u examples/04_daprato_phi4.py
"""
from sympy import Rational

from counterterms import SPDE, Noise, Parabolic, Unknown, daprato_lift, kappa

# Φ⁴₃ :  (∂_t − Δ)u = −u³ + ξ,  d=3,  β₀ = −5/2 − κ  (supercritical → rejected directly)
u = Unknown("u", 3); xi = Noise("xi", regularity=Rational(-5, 2) - kappa)
phi43 = SPDE(noises=[xi], operator=Parabolic(dim=3), unknown=u, rhs=xi.symbol - u.field ** 3)

try:
    phi43.renormalize()
except ValueError as e:
    print("direct renormalize rejects it:", str(e)[:55], "...")

lifted = daprato_lift(phi43)                      # u = X + v; Wick-power noises :X^k:
print("\nlifted v-equation noises:",
      ", ".join(f"{n.name}∈C^{{{n.std}{'' if n.kap == 0 else f'{n.kap:+}κ'}}}" for n in lifted.noises))
res = lifted.renormalize()
print(f"renormalized lifted equation: {len(res.counterterms)} counterterms")
