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

# Unicode super/subscripts.  Only t and x have clean subscript glyphs; other
# coords fall back to ``_name`` (ponytail: covers d≤1 prettily, the common case).
_SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
_SUB_LETTER = {"t": "ₜ", "x": "ₓ"}

_COORDS = ("t", "x", "y", "z")


def coord_names(width: int) -> tuple[str, ...]:
    return tuple(_COORDS[i] if i < len(_COORDS) else f"x{i}" for i in range(width))


def _sub(name: str) -> str:
    return _SUB_LETTER.get(name, "_" + name)


def _noises(sig) -> list[str]:
    return [nt for nt in sig.node_types if nt != "bullet"]


def _node_glyph(node_type: str, sig) -> str:
    if node_type == "bullet":
        return "●"
    return "Ξ" if len(_noises(sig)) == 1 else "Ξ_" + node_type


def _poly_text(n, coords) -> str:
    parts = []
    for i, e in enumerate(n):
        if e:
            s = "X" + _sub(coords[i])
            if e > 1:
                s += str(e).translate(_SUP)
            parts.append(s)
    return "".join(parts)


def _node_token(t, sig, coords) -> str:
    glyph = _node_glyph(t.node_type, sig)
    poly = _poly_text(t.node_dec, coords)
    return f"{poly}·{glyph}" if poly else glyph


def edge_glyph_text(c, p, sig, coords) -> str:
    s = "𝓘"
    if sig.n_components > 1:
        s += str(c).translate(_SUP)               # component of the planting kernel
    s += "".join(_sub(coords[i]) for i, e in enumerate(p) for _ in range(e))
    return s


def _grouped(children):
    """Distinct child edges with multiplicity, preserving canonical order."""
    return Counter(children).items()


# --------------------------------------------------------------------------- #
# shorthand
# --------------------------------------------------------------------------- #

def shorthand(t, sig, coords=None) -> str:
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

def ascii_art(t, sig) -> str:
    coords = coord_names(sig.width)
    lines = [_node_token(t, sig, coords)]
    _walk(t, sig, coords, "", lines)
    return "\n".join(lines)


def _walk(t, sig, coords, prefix, lines):
    kids = list(t.children)                        # expanded multiset, canonically sorted
    for i, (c, p, sub) in enumerate(kids):
        last = i == len(kids) - 1
        branch = "└── " if last else "├── "
        edge = edge_glyph_text(c, p, sig, coords)
        lines.append(prefix + branch + edge + "▸" + _node_token(sub, sig, coords))
        _walk(sub, sig, coords, prefix + ("    " if last else "│   "), lines)


# --------------------------------------------------------------------------- #
# LaTeX forest
# --------------------------------------------------------------------------- #

def forest(t, sig) -> str:
    coords = coord_names(sig.width)
    return "\\begin{forest} rstree\n" + _forest_node(t, sig, coords, 1) + "\n\\end{forest}"


def _node_latex(t, sig, coords) -> str:
    if t.node_type == "bullet":
        sym = "\\bullet"
    else:
        sym = "\\Xi" if len(_noises(sig)) == 1 else "\\Xi_{" + t.node_type + "}"
    poly = _poly_latex(t.node_dec, coords)
    return (poly + " " + sym) if poly else sym


def _poly_latex(n, coords) -> str:
    parts = []
    for i, e in enumerate(n):
        if e:
            s = "X_{" + coords[i] + "}"
            if e > 1:
                s += "^{" + str(e) + "}"
            parts.append(s)
    return "".join(parts)


def _edge_latex(c, p, sig, coords) -> str:
    s = "\\mathcal{I}"
    if sig.n_components > 1:
        s += "^{(" + str(c) + ")}"
    deriv = "".join(coords[i] for i, e in enumerate(p) for _ in range(e))
    if deriv:
        s += "_{" + deriv + "}"
    return s


def _forest_node(t, sig, coords, depth, edge=None) -> str:
    pad = "  " * depth
    line = pad + "[{$" + _node_latex(t, sig, coords) + "$}"
    if edge is not None:
        # forest's `edge label` needs a TikZ node, not bare text; white fill so the
        # kernel label sits cleanly on the edge.
        line += ("," + " edge label={node[midway,fill=white,inner sep=1pt,"
                 "font=\\scriptsize]{$" + edge + "$}}")
    if not t.children:
        return line + "]"
    out = [line]
    for (c, p, sub) in t.children:
        out.append(_forest_node(sub, sig, coords, depth + 1, _edge_latex(c, p, sig, coords)))
    out.append(pad + "]")
    return "\n".join(out)
