"""The SPDE → counterterm-tree pipeline mechanics, in two halves.

* **Generation** (`generate_counterterms`): the derived structural rule (which child
  edges a node admits, with the derivative-slot cap = degree in ∂u), the
  subcriticality bound (only finitely many ``|τ|<0`` trees), and `𝓑_{<0}` including
  the bare primitives (convention #6).
* **The Υ-map** (`elem_diff`, convention #5): base cases ``F(∘)=f``, ``F(●)=g``; the
  child derivatives ``∂_{(c,p)}`` applied **before** ``D^n``; and the child equation
  index ``c`` riding on the edge — how systems couple (convention #7).
"""
from sympy import Derivative, simplify

from regstruct import jet
from regstruct.equation.dsl import build_context
from regstruct.equation.generate import generate_counterterms
from regstruct.renorm.nonlinearity import elem_diff
from regstruct.trees.tree import tree

from tests.conftest import gkpz, gpam2, system

U0, U1 = jet(0, (0, 0)), jet(0, (0, 1))
CIRC = tree("xi", (0, 0))


# --------------------------------------------------------------------------- #
# generation: the rule, criticality, and 𝓑_{<0}
# --------------------------------------------------------------------------- #

def test_generation_counts_and_negativity():
    # gKPZ has exactly five divergent trees; gPAM(d=2) four (tourist_guide.tex 6012).
    sig, _b, _u = build_context(gkpz())
    trees = generate_counterterms(sig)
    assert len(trees) == 5
    assert all(t.homogeneity(sig).is_negative() for t in trees)
    assert len(generate_counterterms(build_context(gpam2())[0])) == 4


def test_bare_primitive_is_generated():
    # convention #6: 𝓑_{<0} contains the bare noise ∘ (the classical constant counterterm).
    sig, _b, _u = build_context(gkpz())
    assert CIRC in generate_counterterms(sig)


def test_gradient_cap_from_assumption_d2():
    # g is quadratic in ∂u, so the derivative edge I_{(0,1)} is capped at multiplicity 2:
    # no generated tree plants three gradient children.
    sig, _b, _u = build_context(gkpz())
    for t in generate_counterterms(sig):
        grad_edges = sum(1 for (_c, p, _s) in t.children if p == (0, 1))
        assert grad_edges <= 2


# --------------------------------------------------------------------------- #
# the Υ-map  F_a(τ*)
# --------------------------------------------------------------------------- #

def test_upsilon_base_cases():
    sig, base, _u = build_context(gkpz())
    # F(∘*) = f(u),  F(●*) = g(u,∂u) = g(u)·(∂ₓu)²
    assert elem_diff(CIRC, 0, base, sig) == base[0]["xi"]
    assert elem_diff(tree("bullet", (0, 0)), 0, base, sig) == base[0]["bullet"]


def test_upsilon_node_decoration_is_total_derivative():
    # F((∘ with X_x)*) = D_x f(u) = u₁ · f'(u)   (the gKPZ ∘1 counterterm)
    sig, base, _u = build_context(gkpz())
    f_u0 = base[0]["xi"]
    expected = U1 * Derivative(f_u0, U0)
    assert simplify(elem_diff(tree("xi", (0, 1)), 0, base, sig) - expected) == 0


def test_upsilon_child_derivative_hits_the_gradient_slot():
    # F((●—I_{(0,1)}—∘)*) = ∂_{u₁}(g·u₁²) · f = 2 g u₁ f  — the factor-2 gKPZ term.
    sig, base, _u = build_context(gkpz())
    t = tree("bullet", (0, 0), [(0, (0, 1), CIRC)])
    g_u0 = base[0]["bullet"].subs(U1, 1)        # g(u0)  (strip the u₁² for the expectation)
    expected = 2 * g_u0 * U1 * base[0]["xi"]
    assert simplify(elem_diff(t, 0, base, sig) - expected) == 0


def test_upsilon_threads_the_edge_component_in_a_system():
    # convention #7: the SAME bare tree ∘ gives the noise coefficient of the OUTPUT
    # equation — a(v) in eq u (comp 0), b(u) in eq v (comp 1).
    sig, base, _u = build_context(system())
    assert simplify(elem_diff(CIRC, 0, base, sig) - base[0]["xi"]) == 0
    assert simplify(elem_diff(CIRC, 1, base, sig) - base[1]["xi"]) == 0
    assert base[0]["xi"] != base[1]["xi"]        # a(v) ≠ b(u): genuinely coupled
