"""The BHZ character (twisted antipode) and machine-readable structure export.
Run:  uv run python -u examples/05_bhz_and_export.py
"""
import json
from sympy import Derivative, Function, Rational

from regstruct import SPDE, Noise, Parabolic, Unknown, kappa, build_renormalization
from regstruct.render import shorthand, structure_json

u = Unknown("u", 1); xi = Noise("xi", regularity=Rational(-1) - kappa)
f, g = Function("f"), Function("g")
spde = SPDE(operator=Parabolic(dim=1, mass=1), unknown=u, noises=[xi],
            rhs=f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2)

rs = build_renormalization(spde)                  # δ, δ⁻, S'₋, and the characters
print("BHZ constant k(τ) = h(S'₋ τ), symbolic in the elementary expectations h_i:")
for t in rs.divergent:
    print(f"  τ={shorthand(t, rs.sig):<14} k(τ) = {rs.bhz_character(t)}")

print("\nCanonical (Gaussian) parity — which constants vanish (odd # noise vertices):")
for t in rs.divergent:
    k = rs.canonical_character(t)
    print(f"  τ={shorthand(t, rs.sig):<14} {'= 0 (canonically)' if k == 0 else '≠ 0 (a genuine divergence)'}")

doc = json.loads(structure_json(spde))            # reconstructible JSON of the whole structure
print(f"\nstructure_json: schema={doc['schema']}, "
      f"{len(doc['divergent'])} divergent trees, coproducts {list(doc['coproducts'])}")
