"""Conventions #3 and #4 — the `DecoratedTree` invariants everything else rests on.

* **#4 canonical isomorphism**: trees are equal/hash-equal iff isomorphic as
  decorated rooted trees (children are an unordered multiset). Equality, set dedup,
  dict keys and `S(τ)` all depend on this.
* **#3 symmetry factor** ``S(τ) = n! · Πⱼ S(σⱼ)^{mⱼ} · mⱼ!`` (tourist_guide.tex 3982).

Plus the homogeneity recursion (naive vs extended, the `o`-decoration on red nodes).
"""
from counterterms.core.homogeneity import Homogeneity
from counterterms.equation.dsl import build_context
from counterterms.trees.tree import red_node, tree

from tests.conftest import gkpz

SIG = build_context(gkpz())[0]
CIRC = tree("xi", (0, 0))                      # ∘  (the noise)


# --------------------------------------------------------------------------- #
# #4 canonical isomorphism
# --------------------------------------------------------------------------- #

def test_child_order_does_not_matter():
    # Same multiset of children built in two orders ⇒ the same tree.
    t1 = tree("bullet", (0, 0), [(0, (0, 1), CIRC), (0, (0, 0), CIRC)])
    t2 = tree("bullet", (0, 0), [(0, (0, 0), CIRC), (0, (0, 1), CIRC)])
    assert t1 == t2 and hash(t1) == hash(t2)


def test_distinct_decorations_are_distinct_trees():
    assert tree("xi", (0, 0)) != tree("xi", (0, 1))                 # node decoration
    assert tree("xi", (0, 0)) != tree("bullet", (0, 0))            # node type
    assert red_node(Homogeneity(-1), width=2) != red_node(Homogeneity(-2), width=2)  # o
    assert CIRC != tree("xi", (0, 0), (), color="red")            # colour


def test_isomorphic_trees_dedup_in_a_set():
    a = tree("bullet", (0, 0), [(0, (0, 1), CIRC), (0, (0, 1), CIRC)])
    b = tree("bullet", (0, 0), [(0, (0, 1), CIRC), (0, (0, 1), CIRC)])
    assert len({a, b}) == 1


def test_node_count():
    assert CIRC.nodes() == 1
    assert tree("xi", (0, 0), [(0, (0, 0), CIRC)]).nodes() == 2


# --------------------------------------------------------------------------- #
# #3 symmetry factor
# --------------------------------------------------------------------------- #

def test_symmetry_factor_node_decoration_factorial():
    assert tree("xi", (0, 0)).symmetry_factor() == 1
    assert tree("xi", (0, 2)).symmetry_factor() == 2            # n! = 2!
    assert tree("xi", (0, 3)).symmetry_factor() == 6            # 3!


def test_symmetry_factor_repeated_children_factorial():
    # two identical child edges ⇒ factor m! = 2! (this is the S=2 gKPZ tree)
    two_same = tree("bullet", (0, 0), [(0, (0, 1), CIRC), (0, (0, 1), CIRC)])
    assert two_same.symmetry_factor() == 2
    # two *distinct* child edges ⇒ no symmetry, S = 1
    two_diff = tree("bullet", (0, 0), [(0, (0, 0), CIRC), (0, (0, 1), CIRC)])
    assert two_diff.symmetry_factor() == 1


def test_symmetry_factor_is_multiplicative_over_subtrees():
    # S(τ) carries the subtree symmetry: child = the S=2 tree, taken twice.
    sub = tree("bullet", (0, 0), [(0, (0, 1), CIRC), (0, (0, 1), CIRC)])   # S(sub)=2
    parent = tree("bullet", (0, 0), [(0, (0, 0), sub), (0, (0, 0), sub)])  # S=2²·2!
    assert parent.symmetry_factor() == 2 ** 2 * 2


# --------------------------------------------------------------------------- #
# homogeneity recursion (naive vs extended)
# --------------------------------------------------------------------------- #

def test_naive_homogeneity_sums_nodes_and_edge_gains():
    assert CIRC.homogeneity(SIG) == Homogeneity(-1, -1)                    # |∘| = β₀
    # ∘ —I₀— ∘ : β₀ + (gain 2) + β₀ = 0 − 2κ
    chain = tree("xi", (0, 0), [(0, (0, 0), CIRC)])
    assert chain.homogeneity(SIG) == Homogeneity(0, -2)
    # node decoration adds |n|_𝔰
    assert tree("xi", (0, 1)).homogeneity(SIG) == Homogeneity(-1, -1) + Homogeneity(1)


def test_extended_homogeneity_adds_the_o_decoration():
    red = red_node(Homogeneity(-1, -1), width=2)
    assert red.homogeneity(SIG) == Homogeneity(0)                 # naive: red counts as 0
    assert red.extended_homogeneity(SIG) == Homogeneity(-1, -1)   # extended: + o
    # a black tree planting a red child: naive ignores o, extended includes it
    t = tree("bullet", (0, 0), [(0, (0, 0), red)])
    assert t.homogeneity(SIG) == Homogeneity(2)                   # 0 + gain 2 + 0
    assert t.extended_homogeneity(SIG) == Homogeneity(2) + Homogeneity(-1, -1)
