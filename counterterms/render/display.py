"""Notebook display: the algebra's results, the equation and the canonical family as pictures.

The coproducts and the twisted antipode return plain dicts (`trees.coproducts.TensorSum` etc.),
which a notebook prints as an unreadable dict of dataclasses.  `TreeSum` is that same dict —
still a ``dict``, equal to the raw result, so nothing downstream changes — carrying the
`Signature` and the shape of its keys, so it can print as a formal sum (``str``) and draw
itself (``_repr_html_``, trees via `render.svg`).  The user-facing entry points
(`Playground`, `RenormalizationStructure`, `RegularityStructure`) wrap their results in it; the
core functions in `trees.coproducts` are untouched.

Everything here emits **inline styles only** (no ``<style>``/``<script>``), so it is safe to
drop into a notebook cell; the full interactive page is `render.html`, embedded by
`explore_html` in an ``<iframe>``.
"""
from __future__ import annotations

from fractions import Fraction
from html import escape
from typing import TYPE_CHECKING, Mapping

from .svg import svg
from .tree import shorthand

if TYPE_CHECKING:
    from ..core.signature import Signature
    from ..renorm.equation import CanonicalFamily, RenormalizedEquation
    from ..trees.tree import DecoratedTree

Forest = tuple["DecoratedTree", ...]

_MAX_TERMS = 40
_MUTED = "opacity:.65"
_TERM = ("display:inline-flex;align-items:center;gap:5px;padding:3px 7px;"
         "border:1px solid rgba(128,128,128,.4);border-radius:6px")


class Picture(str):
    """An HTML fragment that a notebook renders (and that is still a plain ``str``)."""

    def _repr_html_(self) -> str:
        return str(self)


def coeff_text(c: Fraction) -> str:
    """A coefficient as a prefix: ``""`` for 1, ``"−"`` for −1, else ``"3/2·"``."""
    if c == 1:
        return ""
    if c == -1:
        return "−"
    return str(c).replace("-", "−") + "·"


def forest_html(forest: Forest, sig: Signature) -> str:
    """A forest (product of trees) drawn side by side; the empty forest is ``𝟙``."""
    return "".join(svg(t, sig) for t in forest) or "𝟙"


def _forest_text(forest: Forest, sig: Signature) -> str:
    if not forest:
        return "𝟙"
    if len(forest) == 1:
        return shorthand(forest[0], sig)
    return "".join(f"({shorthand(t, sig)})" for t in forest)


class TreeSum(dict):
    """A formal sum of (tensor products of) trees/forests with rational coefficients — the
    result of a coproduct or antipode — that prints and draws itself.

    ``legs`` describes a key: ``("forest",)`` when the key *is* a forest (`ForestSum`, the
    twisted antipode), otherwise one entry per tensor leg, each ``"tree"`` or ``"forest"``
    (``δ``: forest ⊗ tree; ``δ⁻``: forest ⊗ forest; ``Δ``/``Δ⁺``: tree ⊗ tree)."""

    def __init__(self, data: Mapping[object, Fraction], sig: Signature,
                 legs: tuple[str, ...]) -> None:
        super().__init__(data)
        self.sig = sig
        self.legs = legs

    def _legs_of(self, key: object) -> list[Forest]:
        parts = (key,) if len(self.legs) == 1 else key
        return [p if kind == "forest" else (p,) for kind, p in zip(self.legs, parts)]

    def _ordered(self) -> list[tuple[list[Forest], Fraction]]:
        items = [(self._legs_of(k), c) for k, c in self.items()]
        return sorted(items, key=lambda lc: [(len(f), sum(t.nodes() for t in f)) for f in lc[0]])

    def __str__(self) -> str:
        if not self:
            return "0"
        out = ""
        for legs, c in self._ordered():
            body = " ⊗ ".join(_forest_text(f, self.sig) for f in legs)
            if not out:
                out = coeff_text(c) + body
            else:
                out += (" − " + coeff_text(-c) if c < 0 else " + " + coeff_text(c)) + body
        return out

    __repr__ = __str__

    def _repr_html_(self) -> str:
        if not self:
            return "0"
        terms = self._ordered()
        cells = []
        for legs, c in terms[:_MAX_TERMS]:
            body = f'<span style="{_MUTED}">⊗</span>'.join(
                f"<span>{forest_html(f, self.sig)}</span>" for f in legs)
            cells.append(f'<span class="rs-term" style="{_TERM}">{escape(coeff_text(c))}{body}</span>')
        more = (f'<div style="{_MUTED};font-size:90%">… and {len(terms) - _MAX_TERMS} more terms '
                f"(the object is the full dict of {len(terms)})</div>") if len(terms) > _MAX_TERMS else ""
        plus = f'<span style="{_MUTED}">+</span>'
        return ('<div style="display:flex;flex-wrap:wrap;gap:8px 5px;align-items:center">'
                + plus.join(cells) + "</div>" + more)


def tree_picture(t: DecoratedTree, sig: Signature, scale: float = 1.4) -> Picture:
    """``t`` drawn with ``sig`` — unlike a lone tree's own display, this knows the equation's
    noise types and components."""
    return Picture(svg(t, sig, scale=scale))


def summary_html(eq: RenormalizedEquation) -> str:
    """A static drawn summary of ``eq`` for a notebook cell: the free-constant family and
    every divergent tree with its constant, most singular first."""
    from .html import _INERT, _op_html, mathml
    from .report import _sorted_trees, const_map
    import sympy

    sig, kmap = eq.sig, const_map(eq)
    out = ["<div>"]
    for a, (u, op, _rhs) in enumerate(eq.spde.equations):
        # one MathML piece per summand, in tree order — a single SymPy sum would reshuffle them
        cts = sorted(eq.per_component.get(a, []), key=lambda c: c.tree.homogeneity(sig)._key())
        rhs = " + ".join([mathml(eq.original_rhs(a))] + [mathml(eq._pretty(c.term)) for c in cts])
        out.append(f'<div style="margin:4px 0 10px;line-height:2">{_op_html(op)} '
                   f"{mathml(sympy.Symbol(u.name))} = {rhs}</div>")
    out.append('<div style="display:flex;flex-wrap:wrap;gap:10px;align-items:flex-end">')
    card = ("margin:0;padding:6px 8px 3px;border:1px solid rgba(128,128,128,.4);border-radius:7px;"
            "text-align:center;font-size:85%")
    for t in _sorted_trees(eq):
        inert = t not in kmap
        label = _INERT if inert else mathml(kmap[t])
        out.append(f'<figure style="{card}{";opacity:.5" if inert else ""}">{svg(t, sig)}'
                   f"<figcaption>{label}<br>"
                   f'<span style="{_MUTED}">|τ| = {escape(str(t.homogeneity(sig)))}</span>'
                   "</figcaption></figure>")
    out.append("</div></div>")
    return "".join(out)


def explore_html(eq: RenormalizedEquation, canonical: bool = False, reduced: bool = False,
                 symmetric: bool = True, height: int = 640) -> Picture:
    """The interactive explorer (`render.html.html_document`) embedded in an ``<iframe>``, so
    its style and script stay isolated from the notebook page."""
    from .html import html_document
    doc = html_document(eq, canonical, reduced, symmetric)
    return Picture(f'<iframe srcdoc="{escape(doc, quote=True)}" style="width:100%;'
                   f'height:{height}px;border:1px solid rgba(128,128,128,.4);border-radius:6px">'
                   "</iframe>")


def canonical_html(fam: CanonicalFamily) -> str:
    """A `CanonicalFamily` as a table: each tree with its canonical constant, then the legend
    of elementary expectations ``h_i = h_ε(σ_i)`` with ``σ_i`` drawn."""
    import sympy

    from .html import _rate_html, mathml

    sig = fam.eq.sig
    cell = "padding:4px 12px 4px 0;vertical-align:middle;text-align:left"
    rows = []
    for t, v in fam.constants.items():
        rate = "" if v == 0 else "~ " + _rate_html(t.homogeneity(sig))
        rows.append(f'<tr class="rs-const"><td style="{cell}">{svg(t, sig)}</td>'
                    f'<td style="{cell}">{mathml(sympy.Eq(fam.free_constants[t], v, evaluate=False))}'
                    f'</td><td style="{cell};{_MUTED}">{rate}</td></tr>')
    for sym, sigma in fam.legend.items():
        rows.append(f'<tr class="rs-h"><td style="{cell}">{svg(sigma, sig)}</td>'
                    f'<td style="{cell}" colspan="2">{mathml(sym)} = h<sub>ε</sub>(σ), '
                    f'<span style="{_MUTED}">σ the tree on the left</span></td></tr>')
    tag = ("canonical (BPHZ) constants" + (", reduced" if fam.reduced else "")
           + (" — symmetric noise" if fam.reduced and fam.symmetric else ""))
    return (f'<div><div style="{_MUTED};margin-bottom:4px">{tag}</div>'
            f'<table style="border-collapse:collapse">{"".join(rows)}</table></div>')
