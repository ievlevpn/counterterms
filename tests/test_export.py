"""Machine-readable structure export (`render/export.py`): canonical tree round-trip
and a reconstructible JSON document of the full structure (signature, graded `T`,
divergent trees, the coproducts as tensor sums, the BHZ character)."""
import json
from fractions import Fraction

from regstruct.core.homogeneity import Homogeneity
from regstruct.render.export import (
    export_structure, structure_json, tree_from_dict, tree_to_dict)
from regstruct.trees.coproducts import delta_plus
from regstruct.trees.tree import red_node, tree

from tests.conftest import gkpz

CIRC = tree("xi", (0, 0))


def test_tree_roundtrip_is_exact():
    # canonical reconstruction for every flavour of node (decoration, colour, o, nesting)
    trees = [
        CIRC,
        tree("xi", (0, 1)),
        red_node(Homogeneity(-1, -1), width=2),
        tree("bullet", (0, 0), [(0, (0, 1), CIRC), (0, (0, 0), CIRC)]),
        tree("bullet", (0, 0), [(0, (0, 0), red_node(Homogeneity(Fraction(-3, 2), -1), width=2))]),
    ]
    for t in trees:
        assert tree_from_dict(tree_to_dict(t)) == t
    # survives a JSON string round-trip too
    t = trees[-1]
    assert tree_from_dict(json.loads(json.dumps(tree_to_dict(t)))) == t


def test_structure_document_shape():
    doc = export_structure(gkpz())
    assert doc["schema"] == "regstruct/structure/v1"
    assert doc["signature"]["dim"] == 1 and doc["signature"]["scaling"] == [2, 1]
    assert "xi" in doc["signature"]["noises"]
    assert len(doc["divergent"]) == 5                       # gKPZ counterterms
    assert {"recentering", "extraction"} <= doc["coproducts"].keys()
    assert len(doc["bhz_character"]) == 5
    # exact rationals are strings, never floats
    assert doc["divergent"][0]["homogeneity"] == [str(v) for v in doc["divergent"][0]["homogeneity"]]


def test_structure_json_is_valid_and_parses_back():
    doc = json.loads(structure_json(gkpz()))
    assert doc["schema"] == "regstruct/structure/v1"


def test_exported_coproduct_reconstructs_to_the_real_one():
    # take a recentering entry from the export, rebuild the trees, and check it equals
    # the in-memory Δ of that tree (coeffs + tree⊗tree pairs).
    from regstruct.equation.dsl import build_context
    sig = build_context(gkpz())[0]
    doc = export_structure(gkpz())
    entry = next(e for e in doc["coproducts"]["recentering"]
                 if tree_from_dict(e["tree"]) == CIRC)
    rebuilt = {(tree_from_dict(term["left"]), tree_from_dict(term["right"])): Fraction(term["coeff"])
               for term in entry["terms"]}
    assert rebuilt == delta_plus(CIRC, sig)
