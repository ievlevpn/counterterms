"""The generic Hopf machinery (`core/hopf.py`) wired to the two tree Hopf algebras —
``U⁻`` (forests, coproduct ``δ⁻``) and ``T⁺`` (blue trees, coproduct ``Δ⁺``).

The point of the layer is that the *same* generic code serves both bases; the
load-bearing check is the **antipode axiom** ``S⋆id = η∘ε`` on each, plus the
convolution-unit and comodule-counit laws.
"""
from fractions import Fraction

from counterterms.core.hopf import antipode, comodule_action, convolve, counit
from counterterms.equation.dsl import build_context
from counterterms.equation.generate import generate_counterterms
from counterterms.trees.coproducts import _delta_group_forest, delta_minus, delta_plus
from counterterms.trees.tree import tree

from tests.conftest import gkpz

SIG = build_context(gkpz())[0]
W = SIG.width
TREES = generate_counterterms(SIG)


# --- U⁻: the algebra of forests, coproduct δ⁻ ------------------------------- #
def _mul_forest(a, b):
    return tuple(sorted(a + b, key=lambda t: t._sortkey()))


_UNIT_FOREST = ()


def _delta_forest(f):
    return _delta_group_forest(f, SIG)


# --- T⁺: blue-rooted trees, coproduct Δ⁺ ------------------------------------ #
_UNIT_BLUE = tree("bullet", (0,) * W, (), color="blue")


def _mul_blue(a, b):
    nd = tuple(x + y for x, y in zip(a.node_dec, b.node_dec))     # X^j·X^k = X^{j+k}
    return tree("bullet", nd, a.children + b.children, color="blue")


def _delta_blue(b):
    return delta_plus(b, SIG, project_left=True)


def _conv_with_id(S, coproduct, mul):
    """``(S⋆id)(x) = m(S⊗id)Δx`` as a linear combination (for the antipode axiom)."""
    def go(x):
        out: dict = {}
        for (l, r), c in coproduct(x).items():
            for k, v in S(l).items():
                key = mul(k, r)
                out[key] = out.get(key, Fraction(0)) + c * v
        return {k: v for k, v in out.items() if v}
    return go


def test_antipode_axiom_on_U_minus():
    # S⋆id = η∘ε on U⁻ (δ⁻): every generator forest convolves to ε.
    S = antipode(_delta_forest, _mul_forest, _UNIT_FOREST)
    conv = _conv_with_id(S, _delta_forest, _mul_forest)
    assert conv(_UNIT_FOREST) == {_UNIT_FOREST: Fraction(1)}
    for t in TREES:
        assert conv((t,)) == {}                       # ε(x)=0 for x≠𝟙


def test_antipode_axiom_on_T_plus():
    # The structure-group antipode on T⁺ (Δ⁺) — the SAME generic code, new basis.
    S = antipode(_delta_blue, _mul_blue, _UNIT_BLUE)
    conv = _conv_with_id(S, _delta_blue, _mul_blue)
    tplus = {r for t in TREES for (_l, r) in delta_plus(t, SIG)}
    assert conv(_UNIT_BLUE) == {_UNIT_BLUE: Fraction(1)}
    for b in tplus:
        assert conv(b) == ({_UNIT_BLUE: Fraction(1)} if b == _UNIT_BLUE else {})


def test_counit_is_the_convolution_unit():
    # ε⋆f = f for any character f (toy character: integer part of |τ|_std on a forest).
    eps = counit(_UNIT_FOREST)
    f = lambda forest: sum(int(t.homogeneity(SIG).std) for t in forest)
    fstar = convolve(eps, f, _delta_forest)
    for t in TREES:
        assert fstar((t,)) == f((t,))


def test_comodule_action_with_counit_is_identity():
    # k̃ = (k⊗Id)∘δ with k = ε⁻ recovers the model: (ε⁻⊗Id)δ = Id.
    eps_minus = lambda left_forest: 1 if left_forest == () else 0
    kt = comodule_action(eps_minus, lambda t: delta_minus(t, SIG))
    for t in TREES:
        assert kt(t) == {t: 1}
