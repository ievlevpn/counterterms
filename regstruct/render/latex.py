"""Standalone, ``pdflatex``-ready LaTeX of the renormalized equation (notes/output.md §4–6).

Clean single-column layout in Palatino (``newpx``).  Trees are drawn in the
paper's glyph convention (open circle = noise, filled dot = integration node,
solid edge = ``I``, dotted edge = derivative ``∂I``; see ``render/tree.py``),
laid out in a ``booktabs`` table.  Every formula is printed by SymPy (``flatex``
gives ``f'(u)`` prime notation); we only assemble the document.
"""
from __future__ import annotations

import sympy

from .report import const_map, elem_map, flatex, hom_latex, op_latex, _sorted_trees
from .tree import forest

# Node/edge styles match tourist_guide.tex 329–337 (noise / noisegray / noiseblue,
# K = solid kernel, DK = dotted derivative kernel).
_PREAMBLE = r"""\documentclass[11pt]{article}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb}
\let\Bbbk\relax                         % avoid newpxmath/amssymb \Bbbk clash
\usepackage{newpxtext,newpxmath}        % Palatino text + matching math
\usepackage[margin=1.1in]{geometry}
\usepackage{microtype}
\usepackage{booktabs}
\usepackage{forest}
\tikzset{
  noise/.style={circle, draw=black, fill=white, inner sep=0pt, minimum size=4.4pt},
  noisegray/.style={circle, draw=black, fill=black!35, inner sep=0pt, minimum size=4.4pt},
  noiseblue/.style={circle, draw=black, fill=blue!45, inner sep=0pt, minimum size=4.4pt},
  vertex/.style={circle, draw=black, fill=black, inner sep=0pt, minimum size=3.6pt},
}
\forestset{rstree/.style={for tree={grow'=north, parent anchor=north,
  child anchor=south, s sep=9pt, l sep=11pt, edge={semithick}}}}
\setlength{\parindent}{0pt}
\begin{document}
\thispagestyle{empty}"""

_NOISE_FILL = ("white", "black!35", "blue!45")     # matches tree.py _NOISE_STYLE


def _dot(fill, size):
    return (rf"\tikz[baseline=-.5ex]\node[circle,draw,fill={fill},inner sep=0pt,"
            rf"minimum size={size}pt]{{}};")


def _legend(eq) -> str:
    parts = []
    noises = eq.spde.noises
    for i, nz in enumerate(noises):
        fill = _NOISE_FILL[i] if i < len(_NOISE_FILL) else "white"
        label = "noise" if len(noises) == 1 else f"${sympy.latex(nz.symbol)}$"
        parts.append(f"{_dot(fill, 4.4)}\\,= {label}")
    parts.append(_dot("black", 3.6) + r"\,= integration node")
    parts.append(r"solid edge $=\mathcal I$,\; dotted edge $=\partial\mathcal I$ (derivative)")
    return r"{\footnotesize Legend:\; " + r",\quad ".join(parts) + ".}"


def latex_document(eq) -> str:
    sig = eq.sig
    cmap = const_map(eq)
    emap = elem_map(eq)
    scalar = sig.n_components == 1
    P = [_PREAMBLE]

    P.append(r"\begin{center}{\Large\bfseries Renormalized equation}\end{center}")
    P.append(r"\medskip")

    # the SPDE(s)
    P.append(r"\begin{align*}")
    for a, (u, op, _r) in enumerate(eq.spde.equations):
        P.append(f"  {op_latex(op)}\\, {u.name} &= {flatex(eq.original_rhs(a))} \\\\")
    P.append(r"\end{align*}")

    noi = ", ".join(f"${sympy.latex(nz.symbol)} \\in C^{{{hom_latex(nz.homogeneity)}}}$"
                    for nz in eq.spde.noises)
    P.append(f"with {noi}, parabolic scaling "
             f"$\\mathfrak s = ({', '.join(map(str, sig.scaling.weights))})$, "
             f"spatial dimension $d = {sig.dim}$.")
    if any(o != 2 for o in sig.comp_order):
        P.append(r"\\[2pt] \emph{Non-standard operator order: the theory is proven "
                 r"for 2nd-order parabolic $L$.}")

    # divergent trees
    P.append(r"\section*{Divergent trees}")
    P.append(r"\begin{center}")
    P.append(r"\renewcommand{\arraystretch}{1.2}")
    P.append(r"\begin{tabular}{" + ("c c c c l" if scalar else "c c c c") + "}")
    P.append(r"\toprule")
    P.append(r"$\tau$ & $|\tau|$ & $S(\tau)$ & $k_\tau$"
             + (r" & $F(\tau^*)$ \\" if scalar else r" \\"))
    P.append(r"\midrule")
    for t in _sorted_trees(eq):
        k = cmap.get(t)
        kc = sympy.latex(k) if k is not None else r"\text{--}"
        row = (f"{forest(t, sig)} & ${hom_latex(t.homogeneity(sig))}$ "
               f"& ${t.symmetry_factor()}$ & ${kc}$")
        if scalar:
            fl = flatex(eq._pretty(emap[t])) if k is not None else "0"
            row += f" & ${fl}$"
        P.append(row + r" \\[4pt]")
    P.append(r"\bottomrule")
    P.append(r"\end{tabular}")
    P.append(r"\end{center}")
    P.append(_legend(eq))

    # renormalized family
    P.append(r"\section*{Renormalized family}")
    for a in range(eq.n_components):
        u, op = eq.unknowns[a], eq.spde.equations[a][1]
        P.append(r"\begin{align*}")
        P.append(f"  {op_latex(op)}\\, {u.name} &= {flatex(eq.original_rhs(a))} \\\\")
        for ct in sorted(eq.per_component[a], key=lambda c: c.homogeneity._key()):
            term = sympy.expand(ct.coefficient * eq._pretty(ct.elem_diff))
            P.append(r"  &\quad + " + flatex(term) + r" \\")
        P.append(r"\end{align*}")

    P.append(r"\bigskip")
    P.append(r"{\itshape The algebraic data ($\Delta\tau$, $\Delta^-\tau$, $S'_-\tau$, "
             r"symbolic BHZ character) is available via \texttt{structures.py}; rendering it "
             r"here and the canonical numeric constants (Phase~4) are pending.}")
    P.append(r"\end{document}")
    return "\n".join(P)
