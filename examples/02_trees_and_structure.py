"""The divergent trees, the regularity structure (T, T⁺), and the coproducts.
Run:  uv run python -u examples/02_trees_and_structure.py
"""
from sympy import Derivative, Function, Rational

from regstruct import SPDE, Noise, Parabolic, Unknown, kappa, build_regularity_structure
from regstruct.render import shorthand

u = Unknown("u", 1); xi = Noise("xi", regularity=Rational(-1) - kappa)
f, g = Function("f"), Function("g")
spde = SPDE(operator=Parabolic(dim=1, mass=1), unknown=u, noises=[xi],
            rhs=f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2)

rs = build_regularity_structure(spde)            # the γ-truncated (T, T⁺)
print("Model-space basis T, graded by homogeneity |τ|:")
for h, trees in sorted(rs.grades().items(), key=lambda kv: kv[0]._key()):
    print(f"  |τ|={str(h):>8}:  " + ", ".join(shorthand(t, rs.sig) for t in trees))

print(f"\nDivergent subspace (counterterm carriers): {len(rs.divergent)} trees")
print("\nRecentering Δ of the noise ∘ (Δ : T → T ⊗ T⁺):")
circ = next(t for t in rs.model_basis if t.nodes() == 1 and t.node_type == "xi")
for (left, right), c in rs.recentering(circ).items():
    print(f"  {c} · {shorthand(left, rs.sig)} ⊗ {shorthand(right, rs.sig)}")
