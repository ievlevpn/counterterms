"""Three faithful renderings of one ``DecoratedTree`` (notes/output.md §3).

- ``shorthand``  — one-line product notation, e.g. ``●·𝓘ₓ[Ξ]²``
- ``ascii_art``  — 2D terminal drawing, root at top (the ``tree(1)`` style)
- ``forest``     — a LaTeX ``forest`` snippet (the paper's upward rooted trees)

All three read off the same datatype (node type + node decoration + child edges
``(component, p, subtree)``); coordinate names come from the unknown's scaling
width, so ``p=(0,1)`` prints as an *x*-derivative, not index soup.
"""
from __future__ import annotations

from collections import Counter
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..core.homogeneity import Homogeneity, MultiIndex
    from ..core.signature import Signature
    from ..trees.tree import DecoratedTree

# Unicode super/subscripts.  Only t and x have clean subscript glyphs; other
# coords fall back to ``_name`` (ponytail: covers d≤1 prettily, the common case).
_SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
_SUB_LETTER = {"t": "ₜ", "x": "ₓ"}

_COORDS = ("t", "x", "y", "z")


def coord_names(width: int) -> tuple[str, ...]:
    return tuple(_COORDS[i] if i < len(_COORDS) else f"x{i}" for i in range(width))


def _sub(name: str) -> str:
    return _SUB_LETTER.get(name, "_" + name)


def _noises(sig: Signature) -> list[str]:
    return [nt for nt in sig.node_types if nt != "bullet"]


def _node_glyph(node_type: str, sig: Signature) -> str:
    if node_type == "bullet":
        return "●"
    return "Ξ" if len(_noises(sig)) == 1 else "Ξ_" + node_type


def _poly_text(n: MultiIndex, coords: tuple[str, ...]) -> str:
    parts = []
    for i, e in enumerate(n):
        if e:
            s = "X" + _sub(coords[i])
            if e > 1:
                s += str(e).translate(_SUP)
            parts.append(s)
    return "".join(parts)


def _node_token(t: DecoratedTree, sig: Signature, coords: tuple[str, ...]) -> str:
    glyph = _node_glyph(t.node_type, sig)
    poly = _poly_text(t.node_dec, coords)
    return f"{poly}·{glyph}" if poly else glyph


def edge_glyph_text(c: int, p: MultiIndex, sig: Signature, coords: tuple[str, ...]) -> str:
    s = "𝓘"
    if sig.n_components > 1:
        s += str(c).translate(_SUP)               # component of the planting kernel
    s += "".join(_sub(coords[i]) for i, e in enumerate(p) for _ in range(e))
    return s


def _grouped(
    children: tuple[tuple[int, MultiIndex, DecoratedTree], ...],
) -> object:
    """Distinct child edges with multiplicity, preserving canonical order."""
    return Counter(children).items()


# --------------------------------------------------------------------------- #
# shorthand
# --------------------------------------------------------------------------- #

def shorthand(t: DecoratedTree, sig: Signature, coords: tuple[str, ...] | None = None) -> str:
    if coords is None:
        coords = coord_names(sig.width)
    factors = [_node_token(t, sig, coords)]
    for (c, p, sub), m in _grouped(t.children):
        e = f"{edge_glyph_text(c, p, sig, coords)}[{shorthand(sub, sig, coords)}]"
        if m > 1:
            e += str(m).translate(_SUP)
        factors.append(e)
    return "·".join(factors)


# --------------------------------------------------------------------------- #
# 2D terminal drawing  (fixed 4-char indent; edge embedded in the child token)
# --------------------------------------------------------------------------- #

def ascii_art(t: DecoratedTree, sig: Signature) -> str:
    coords = coord_names(sig.width)
    lines = [_node_token(t, sig, coords)]
    _walk(t, sig, coords, "", lines)
    return "\n".join(lines)


def _walk(
    t: DecoratedTree,
    sig: Signature,
    coords: tuple[str, ...],
    prefix: str,
    lines: list[str],
) -> None:
    kids = list(t.children)                        # expanded multiset, canonically sorted
    for i, (c, p, sub) in enumerate(kids):
        last = i == len(kids) - 1
        branch = "└── " if last else "├── "
        edge = edge_glyph_text(c, p, sig, coords)
        lines.append(prefix + branch + edge + "▸" + _node_token(sub, sig, coords))
        _walk(sub, sig, coords, prefix + ("    " if last else "│   "), lines)


# --------------------------------------------------------------------------- #
# LaTeX forest  — the paper's glyph convention (tourist_guide.tex 329–337):
#   noise = open circle ○ , integration node = filled dot ● , kernel edge = solid
#   line, derivative kernel ∂I = dotted line.  Node decorations / coordinate /
#   component disambiguation appear as small labels only when needed.
# --------------------------------------------------------------------------- #

# Fill styles per noise index (matches the paper's noise / noisegray / noiseblue).
_NOISE_STYLE = ("noise", "noisegray", "noiseblue")


def forest(t: DecoratedTree, sig: Signature) -> str:
    # Wrap in $\vcenter{\hbox{…}}$: gives the picture a correct hbox width (so a table
    # column sizes to the tree) and centers it on the math axis.  A bare forest with a
    # coordinate `baseline` reports a wrong width and overflows its cell.  Usable in
    # text mode (the $…$) and inside surrounding math.
    coords = coord_names(sig.width)
    return ("$\\vcenter{\\hbox{\\begin{forest} rstree\n"
            + _forest_node(t, sig, coords, 1) + "\n\\end{forest}}}$")


def _poly_latex(n: MultiIndex, coords: tuple[str, ...]) -> str:
    parts = []
    for i, e in enumerate(n):
        if e:
            s = "X_{" + coords[i] + "}"
            if e > 1:
                s += "^{" + str(e) + "}"
            parts.append(s)
    return "".join(parts)


def _o_latex(o: Homogeneity) -> str:
    return str(o).replace("κ", r"\kappa")


def _node_opts(t: DecoratedTree, sig: Signature, coords: tuple[str, ...]) -> list[str]:
    # Phase-3 extended decoration: red = contraction node (square, carries o),
    # blue = T⁺ root.  Black nodes keep the paper's circle/dot convention.
    if t.color == "red":
        opts = ["redvertex"]
        if t.o.std != 0 or t.o.kap != 0:
            opts.append("label={[font=\\tiny,fill=white,fill opacity=0.7,text opacity=1,"
                        "inner sep=0.8pt]left:$o{=}" + _o_latex(t.o) + "$}")
    elif t.color == "blue":
        opts = ["bluevertex"]
    elif t.node_type == "bullet":
        opts = ["vertex"]
    else:
        noises = _noises(sig)
        i = noises.index(t.node_type)
        opts = [_NOISE_STYLE[i] if i < len(_NOISE_STYLE) else "noise"]
        if i >= len(_NOISE_STYLE):                         # too many noises to colour
            opts.append("label={[font=\\tiny]below:$\\Xi_{" + t.node_type + "}$}")
    poly = _poly_latex(t.node_dec, coords)
    if poly:
        opts.append("label={[font=\\tiny,inner sep=1.5pt]above right:$" + poly + "$}")
    return opts


def _edge_opts(c: int, p: MultiIndex, sig: Signature, coords: tuple[str, ...]) -> list[str]:
    deriv = any(p)
    opts = ["edge={densely dotted,thick}"] if deriv else []
    labels = []
    if sig.n_components > 1:                               # which kernel/component
        labels.append("(" + str(c) + ")")
    if deriv and sig.dim > 1:                              # which derivative direction
        labels.append("".join(coords[i] for i, e in enumerate(p) for _ in range(e)))
    if labels:
        opts.append("edge label={node[midway,fill=white,inner sep=1pt,font=\\tiny]"
                    "{$" + ",".join(labels) + "$}}")
    return opts


def _forest_node(
    t: DecoratedTree,
    sig: Signature,
    coords: tuple[str, ...],
    depth: int,
    edge_opts: list[str] | None = None,
) -> str:
    pad = "  " * depth
    opts = _node_opts(t, sig, coords) + (edge_opts or [])
    head = pad + "[{}, " + ", ".join(opts)
    if not t.children:
        return head + "]"
    out = [head]
    for (c, p, sub) in t.children:
        out.append(_forest_node(sub, sig, coords, depth + 1, _edge_opts(c, p, sig, coords)))
    out.append(pad + "]")
    return "\n".join(out)
