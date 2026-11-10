"""A fourth rendering of a ``DecoratedTree``: a standalone inline ``<svg>``.

Same glyph convention as the LaTeX ``forest`` drawer (`render.tree.forest`, tourist_guide.tex
329–337) — noise = open circle, integration node = filled dot, kernel edge = solid line,
derivative kernel = dotted line, red contraction node = red square (with its ``o``), blue
``T⁺`` root = blue dot — but with no LaTeX: it shows in a notebook (`DecoratedTree._repr_html_`)
and inside the HTML explorer (`render.html`).  Root at the bottom, children above.

Node ids (``data-id``) are the pre-order index over the canonical children — the numbering of
`trees.coproducts._explode`, hence of the ``z_i`` in `renorm.scheme.expectation`; ``pairing``
draws a dashed arc per Wick-paired couple of noise nodes.

Black elements use ``currentColor`` so the picture follows the surrounding text colour.
"""
from __future__ import annotations

from html import escape
from typing import TYPE_CHECKING, Mapping

from .tree import _display_sig, _noises, _poly_text, coord_names

if TYPE_CHECKING:
    from ..core.homogeneity import MultiIndex
    from ..core.signature import Signature
    from ..trees.tree import DecoratedTree

_SLOT, _STEP, _MARGIN, _FONT = 30.0, 34.0, 8.0, 8.0
_R_DOT, _R_NOISE, _HALF_SQ = 3.2, 4.2, 3.5
_RED, _BLUE, _GREY = "#c0392b", "#2c6fbb", "#9aa0a6"
_R_HALO = 7.5
_CUT = ("#e69f00", "#56b4e9", "#009e73", "#cc79a7")       # one colour per extracted block
_NOISE_FILL = ("none", _GREY, "#4a7bd0")        # matches forest's noise / noisegray / noiseblue


class _Canvas:
    """Collects SVG elements and the bounding box of everything drawn."""

    def __init__(self) -> None:
        self.under: list[str] = []          # cut shading, drawn beneath everything
        self.edges: list[str] = []
        self.shapes: list[str] = []
        self.x0 = self.y0 = float("inf")
        self.x1 = self.y1 = float("-inf")

    def grow(self, xa: float, ya: float, xb: float, yb: float) -> None:
        self.x0, self.y0 = min(self.x0, xa, xb), min(self.y0, ya, yb)
        self.x1, self.y1 = max(self.x1, xa, xb), max(self.y1, ya, yb)

    def text(self, x: float, y: float, s: str, anchor: str = "start", cls: str = "lbl") -> None:
        w = 0.6 * _FONT * len(s)
        left = {"start": x, "middle": x - w / 2, "end": x - w}[anchor]
        self.grow(left, y - _FONT, left + w, y + 2)
        self.shapes.append(
            f'<text class="{cls}" x="{x:.1f}" y="{y:.1f}" font-size="{_FONT:g}" '
            f'text-anchor="{anchor}" fill="currentColor">{escape(s)}</text>')


def _radius(t: DecoratedTree) -> float:
    if t.color == "red":
        return _HALF_SQ
    return _R_DOT if t.node_type == "bullet" else _R_NOISE


def _edge_label(c: int, p: MultiIndex, sig: Signature, coords: tuple[str, ...]) -> str:
    # same conditions as the forest drawer's `_edge_opts`: component for systems, direction
    # letters for a derivative kernel in d > 1.
    parts = []
    if sig.n_components > 1:
        parts.append(f"({c})")
    if any(p) and sig.dim > 1:
        parts.append("".join(coords[i] for i, e in enumerate(p) for _ in range(e)))
    return ",".join(parts)


def svg(t: DecoratedTree, sig: Signature | None = None, *,
        pairing: tuple[tuple[int, int], ...] = (), scale: float = 1.0,
        mark: Mapping[int, int] | None = None, number: bool = False) -> str:
    """Draw ``t`` as a standalone ``<svg>`` string (see the module docstring).

    ``mark`` shades a cut: ``{node id: block index}`` puts a coloured halo behind each listed
    node and a coloured band along every edge joining two nodes of the same block (an
    extracted subtree of `trees.coproducts.cuts`).  ``number`` writes each node's id beside it."""
    if sig is None:
        sig = _display_sig(t)
    coords = coord_names(sig.width)
    noises = _noises(sig)

    # ---- layout: pre-order ids; leaves take consecutive slots, parents the mean of children
    nodes: list[DecoratedTree] = []
    parent: dict[int, tuple[int, int, MultiIndex]] = {}       # child id -> (parent id, comp, p)
    depth: dict[int, int] = {}
    xs: dict[int, float] = {}
    next_slot = 0

    def place(sub: DecoratedTree, d: int) -> int:
        nonlocal next_slot
        i = len(nodes)
        nodes.append(sub)
        depth[i] = d
        kids = []
        for (c, p, child) in sub.children:
            k = place(child, d + 1)
            parent[k] = (i, c, p)
            kids.append(k)
        if kids:
            xs[i] = sum(xs[k] for k in kids) / len(kids)
        else:
            xs[i] = next_slot * _SLOT
            next_slot += 1
        return i

    place(t, 0)
    top = max(depth.values())
    ys = {i: (top - depth[i]) * _STEP for i in depth}         # root lowest (SVG y grows down)

    cv = _Canvas()
    mark = mark or {}

    # ---- cut shading: a halo per marked node, a band per edge internal to one block
    for i, g in mark.items():
        colour = _CUT[g % len(_CUT)]
        cv.under.append(f'<circle class="cut-node cut-{g}" cx="{xs[i]:.1f}" cy="{ys[i]:.1f}" '
                        f'r="{_R_HALO:g}" fill="{colour}" fill-opacity="0.45"/>')
        cv.grow(xs[i] - _R_HALO, ys[i] - _R_HALO, xs[i] + _R_HALO, ys[i] + _R_HALO)
    for k, (i, _c, _p) in parent.items():
        if i in mark and mark.get(k) == mark[i]:
            cv.under.append(
                f'<line class="cut-edge cut-{mark[i]}" x1="{xs[i]:.1f}" y1="{ys[i]:.1f}" '
                f'x2="{xs[k]:.1f}" y2="{ys[k]:.1f}" stroke="{_CUT[mark[i] % len(_CUT)]}" '
                f'stroke-opacity="0.45" stroke-width="{_R_HALO:g}" stroke-linecap="round"/>')

    # ---- edges (trimmed to the node outlines, so an open noise circle stays empty)
    for k, (i, c, p) in parent.items():
        dx, dy = xs[k] - xs[i], ys[k] - ys[i]
        length = (dx * dx + dy * dy) ** 0.5
        ux, uy = dx / length, dy / length
        ri, rk = _radius(nodes[i]), _radius(nodes[k])
        xa, ya, xb, yb = xs[i] + ux * ri, ys[i] + uy * ri, xs[k] - ux * rk, ys[k] - uy * rk
        deriv = any(p)
        cv.edges.append(
            f'<line class="{"e e-d" if deriv else "e"}" x1="{xa:.1f}" y1="{ya:.1f}" '
            f'x2="{xb:.1f}" y2="{yb:.1f}" stroke="currentColor" stroke-width="1.1"'
            + (' stroke-dasharray="1.5 2.5" stroke-linecap="round"' if deriv else "") + "/>")
        cv.grow(xa, ya, xb, yb)
        lab = _edge_label(c, p, sig, coords)
        if lab:
            # on the outer side of the edge, so sibling edges do not run through it
            mx, my = (xs[i] + xs[k]) / 2, (ys[i] + ys[k]) / 2 + 3
            if xs[k] < xs[i]:
                cv.text(mx - 4, my, lab, "end")
            else:
                cv.text(mx + 4, my, lab)

    # ---- nodes
    red_leaves = 0
    for i, n in enumerate(nodes):
        x, y = xs[i], ys[i]
        if n.color == "red":
            h = _HALF_SQ
            cv.shapes.append(f'<rect class="n n-red" data-id="{i}" x="{x - h:.1f}" y="{y - h:.1f}" '
                             f'width="{2 * h:g}" height="{2 * h:g}" fill="{_RED}"/>')
            cv.grow(x - h, y - h, x + h, y + h)
            if n.o.std != 0 or n.o.kap != 0:
                if n.children:                    # edges leave upward: label beside the node,
                    if i in parent and x <= xs[parent[i][0]]:    # on the side away from siblings
                        cv.text(x - h - 3, y + 3, f"o={n.o}", "end")
                    else:
                        cv.text(x + h + 3, y + 3, f"o={n.o}")
                else:                             # leaf: above; stagger so neighbours don't collide
                    cv.text(x, y - h - 3 - 9 * (red_leaves % 2), f"o={n.o}", "middle")
                    red_leaves += 1
        elif n.color == "blue":
            cv.shapes.append(f'<circle class="n n-blue" data-id="{i}" cx="{x:.1f}" cy="{y:.1f}" '
                             f'r="{_R_DOT:g}" fill="{_BLUE}"/>')
            cv.grow(x - _R_DOT, y - _R_DOT, x + _R_DOT, y + _R_DOT)
        elif n.node_type == "bullet":
            cv.shapes.append(f'<circle class="n n-bullet" data-id="{i}" cx="{x:.1f}" cy="{y:.1f}" '
                             f'r="{_R_DOT:g}" fill="currentColor"/>')
            cv.grow(x - _R_DOT, y - _R_DOT, x + _R_DOT, y + _R_DOT)
        else:
            j = noises.index(n.node_type)
            fill = _NOISE_FILL[j] if j < len(_NOISE_FILL) else "none"
            cv.shapes.append(
                f'<circle class="n n-noise n-noise-{j}" data-id="{i}" cx="{x:.1f}" cy="{y:.1f}" '
                f'r="{_R_NOISE:g}" fill="{fill}" stroke="currentColor" stroke-width="1.1"/>')
            cv.grow(x - _R_NOISE - 1, y - _R_NOISE - 1, x + _R_NOISE + 1, y + _R_NOISE + 1)
            if j >= len(_NOISE_FILL):                        # too many noises to colour: name it
                cv.text(x, y - _R_NOISE - 4, "Ξ_" + n.node_type, "middle")   # above: clear of edges
        if number:
            cv.text(x - _radius(n) - 3, y + 3, str(i), "end", "num")
        poly = _poly_text(n.node_dec, coords)
        if poly:
            # a leaf has free space above; an inner node's edges leave upward, so go below
            cv.text(x + _radius(n) + 2, y + 8 if n.children else y - 3, poly)

    # ---- Wick pairings: a dashed arc — bulging above two side-by-side noises, or out to the
    # right of two stacked ones (same x), where an upward bulge would lie on the kernel edge
    for (a, b) in pairing:
        off = _R_NOISE + 1
        if abs(xs[a] - xs[b]) < 1:
            x = xs[a] + off
            cx, cy = x + 20 + abs(ys[a] - ys[b]) * 0.15, (ys[a] + ys[b]) / 2
            d = f"M {x:.1f} {ys[a]:.1f} Q {cx:.1f} {cy:.1f} {x:.1f} {ys[b]:.1f}"
        else:
            cx, cy = (xs[a] + xs[b]) / 2, min(ys[a], ys[b]) - 16 - abs(xs[a] - xs[b]) * 0.15
            d = (f"M {xs[a]:.1f} {ys[a] - off:.1f} Q {cx:.1f} {cy:.1f} "
                 f"{xs[b]:.1f} {ys[b] - off:.1f}")
        cv.shapes.append(f'<path class="pair" d="{d}" fill="none" stroke="{_RED}" '
                         f'stroke-width="1" stroke-dasharray="3 2"/>')
        cv.grow(min(xs[a], xs[b], cx), min(ys[a], ys[b], cy), max(xs[a], xs[b], cx),
                max(ys[a], ys[b], cy))

    vx, vy = cv.x0 - _MARGIN, cv.y0 - _MARGIN
    vw, vh = cv.x1 - cv.x0 + 2 * _MARGIN, cv.y1 - cv.y0 + 2 * _MARGIN
    return (f'<svg xmlns="http://www.w3.org/2000/svg" class="rs-tree" '
            f'viewBox="{vx:.1f} {vy:.1f} {vw:.1f} {vh:.1f}" '
            f'width="{vw * scale:.1f}" height="{vh * scale:.1f}" '
            f'style="vertical-align:middle">'
            + "".join(cv.under) + "".join(cv.edges) + "".join(cv.shapes) + "</svg>")
