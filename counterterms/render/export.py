"""Machine-readable export of the full built structure — interop glue.

Complements the human-oriented ``render('json')`` (which uses lossy *shorthand*
strings) with a **versioned, reconstructible** JSON document for an external analytic
or numerical consumer:

* trees are encoded **canonically** and round-trip exactly (`tree_to_dict` /
  `tree_from_dict`);
* the coproducts ``Δ`` (recentering) and ``δ⁻`` (extraction, group) are emitted as
  explicit ``tree ⊗ tree`` sums;
* the graded model basis ``T``, the divergent trees with ``S(τ)`` and ``F(τ*)``, the
  signature/rule, and the symbolic BHZ character are all included.

Exact rationals are serialized as strings (``"−1/2"``), never floats.
"""
from __future__ import annotations

import json
from fractions import Fraction
from typing import TYPE_CHECKING, Iterable

from ..core.homogeneity import Homogeneity
from ..trees.tree import DecoratedTree, tree as _build_tree

if TYPE_CHECKING:
    from ..equation.dsl import SPDE

SCHEMA = "counterterms/structure/v1"


# --------------------------------------------------------------------------- #
# canonical, round-trippable tree encoding
# --------------------------------------------------------------------------- #

def tree_to_dict(t: DecoratedTree) -> dict:
    """Encode a `DecoratedTree` as a plain dict (canonical, fully reconstructible)."""
    return {
        "type": t.node_type,
        "dec": list(t.node_dec),
        "color": t.color,
        "o": [str(t.o.std), str(t.o.kap)],
        "children": [[c, list(p), tree_to_dict(sub)] for (c, p, sub) in t.children],
    }


def tree_from_dict(d: dict) -> DecoratedTree:
    """Inverse of `tree_to_dict` (rebuilds the canonical tree)."""
    o = Homogeneity(Fraction(d["o"][0]), Fraction(d["o"][1]))
    children = [(c, tuple(p), tree_from_dict(sub)) for (c, p, sub) in d["children"]]
    return _build_tree(d["type"], tuple(d["dec"]), children, color=d["color"], o=o)


def _hom(h: Homogeneity) -> list:
    return [str(h.std), str(h.kap)]


def _forest(f: Iterable[DecoratedTree]) -> list:
    return [tree_to_dict(t) for t in f]


# --------------------------------------------------------------------------- #
# the structure document
# --------------------------------------------------------------------------- #

def export_structure(spde: SPDE, gamma: Fraction = Fraction(1)) -> dict:
    """Build `spde`'s regularity/renormalization structure and return it as a
    JSON-serializable dict (schema ``counterterms/structure/v1``)."""
    from ..structures import build_regularity_structure, build_renormalization
    from ..trees.coproducts import delta_minus_group, delta_plus
    from ..renorm.nonlinearity import elem_diff

    from ..equation.dsl import build_context

    rs = build_regularity_structure(spde, gamma)
    rn = build_renormalization(spde)
    sig = rs.sig
    _sig, base_nl, _u = build_context(spde)        # per-equation nonlinearities for F(τ*)

    def _delta_plus_terms(t: DecoratedTree) -> list:
        return [{"coeff": str(c), "left": tree_to_dict(l), "right": tree_to_dict(r)}
                for (l, r), c in delta_plus(t, sig).items()]

    def _delta_minus_terms(t: DecoratedTree) -> list:
        return [{"coeff": str(c), "left": _forest(l), "right": _forest(r)}
                for (l, r), c in delta_minus_group(t, sig).items()]

    divergent = []
    for t in rs.divergent:
        divergent.append({
            "tree": tree_to_dict(t),
            "homogeneity": _hom(t.homogeneity(sig)),
            "symmetry": t.symmetry_factor(),
            "elem_diff": str(elem_diff(t, 0, base_nl, sig)),
        })

    return {
        "schema": SCHEMA,
        "signature": {
            "dim": sig.dim,
            "scaling": list(sig.scaling.weights),
            "comp_order": list(sig.comp_order),
            "n_components": sig.n_components,
            "noises": {name: _hom(h) for name, h in sig.noise_homog.items()},
            "node_types": list(sig.node_types),
            "rule": {b: [[c, list(p), cap] for (c, p, cap) in rules]
                     for b, rules in sig.allowed.items()},
        },
        "model_basis": [
            {"homogeneity": _hom(h), "trees": [tree_to_dict(t) for t in trees]}
            for h, trees in sorted(rs.grades().items(), key=lambda kv: kv[0]._key())
        ],
        "divergent": divergent,
        "coproducts": {
            "recentering": [{"tree": tree_to_dict(t), "terms": _delta_plus_terms(t)}
                            for t in rs.model_basis],
            "extraction": [{"tree": tree_to_dict(t), "terms": _delta_minus_terms(t)}
                           for t in rs.divergent],
        },
        "bhz_character": [
            {"tree": tree_to_dict(t), "constant": str(rn.bhz_character(t))}
            for t in rn.divergent
        ],
    }


def structure_json(spde: SPDE, gamma: Fraction = Fraction(1), indent: int = 2) -> str:
    """`export_structure` as a JSON string."""
    return json.dumps(export_structure(spde, gamma), indent=indent, ensure_ascii=False)
