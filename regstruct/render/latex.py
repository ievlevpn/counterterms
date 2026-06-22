"""Standalone, ``pdflatex``-ready LaTeX of the renormalized equation (notes/output.md §4–6).

Clean single-column layout in Palatino (``newpx``).  Trees are drawn in the
paper's glyph convention (open circle = noise, filled dot = integration node,
solid edge = ``I``, dotted edge = derivative ``∂I``; see ``render/tree.py``),
laid out in a ``booktabs`` table.  Every formula is printed by SymPy (``flatex``
gives ``f'(u)`` prime notation); we only assemble the document.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import sympy

from .report import (canonical_data, const_map, elem_map, flatex,
                     hom_latex, op_latex, _sorted_trees)
from .tree import forest

if TYPE_CHECKING:
    from ..renorm.equation import RenormalizedEquation

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
\usepackage{longtable}                   % the trees table may span pages
\usepackage{forest}
\tikzset{
  noise/.style={circle, draw=black, fill=white, inner sep=0pt, minimum size=4.4pt},
  noisegray/.style={circle, draw=black, fill=black!35, inner sep=0pt, minimum size=4.4pt},
  noiseblue/.style={circle, draw=black, fill=blue!45, inner sep=0pt, minimum size=4.4pt},
  vertex/.style={circle, draw=black, fill=black, inner sep=0pt, minimum size=3.6pt},
  redvertex/.style={rectangle, draw=red!70!black, fill=red!30, inner sep=0pt, minimum size=4pt},
  bluevertex/.style={circle, draw=blue!60!black, fill=blue!20, inner sep=0pt, minimum size=3.6pt},
}
\forestset{rstree/.style={for tree={grow'=north, parent anchor=north,
  child anchor=south, s sep=14pt, l sep=14pt, edge={semithick}}}}
\setlength{\parindent}{0pt}
\begin{document}
\thispagestyle{empty}"""

_NOISE_FILL = ("white", "black!35", "blue!45")     # matches tree.py _NOISE_STYLE


def _dot(fill: str, size: float) -> str:
    return (rf"\tikz[baseline=-.5ex]\node[circle,draw,fill={fill},inner sep=0pt,"
            rf"minimum size={size}pt]{{}};")


def _legend(eq: RenormalizedEquation) -> str:
    parts = []
    noises = eq.spde.noises
    for i, nz in enumerate(noises):
        fill = _NOISE_FILL[i] if i < len(_NOISE_FILL) else "white"
        label = "noise" if len(noises) == 1 else f"${sympy.latex(nz.symbol)}$"
        parts.append(f"{_dot(fill, 4.4)}\\,= {label}")
    parts.append(_dot("black", 3.6) + r"\,= integration node")
    parts.append(r"solid edge $=\mathcal I$,\; dotted edge $=\partial\mathcal I$ (derivative)")
    return r"{\footnotesize Legend:\; " + r",\quad ".join(parts) + ".}"


def latex_document(eq: RenormalizedEquation, canonical: bool = False) -> str:
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

    # divergent trees — longtable so many-tree examples (KPZ, …) break across pages
    P.append(r"\section*{Divergent trees}")
    P.append(r"\renewcommand{\arraystretch}{1.2}")
    P.append(r"\begin{longtable}{@{}c c c c l@{}}" if scalar
             else r"\begin{longtable}{@{}c c c c@{}}")
    P.append(r"\toprule")
    P.append(r"$\tau$ & $|\tau|$ & $S(\tau)$ & $k_\tau$"
             + (r" & $F(\tau^*)$ \\" if scalar else r" \\"))
    P.append(r"\midrule\endhead")
    P.append(r"\bottomrule\endfoot")
    for t in _sorted_trees(eq):
        k = cmap.get(t)
        kc = sympy.latex(k) if k is not None else r"\text{--}"
        row = (f"{forest(t, sig)} & ${hom_latex(t.homogeneity(sig))}$ "
               f"& ${t.symmetry_factor()}$ & ${kc}$")
        if scalar:
            fl = flatex(eq._pretty(emap[t])) if k is not None else "0"
            row += f" & ${fl}$"
        P.append(row + r" \\[4pt]")
    P.append(r"\end{longtable}")
    P.append(r"\begin{center}")
    P.append(_legend(eq))
    P.append(r"\end{center}")

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

    # canonical (BPHZ) renormalization — parity-reduced, symbolic in the expectations h(σ)
    if canonical:
        rows, legend = canonical_data(eq)
        nzero = sum(1 for _t, _k, v in rows if v == 0)
        P.append(r"\section*{Canonical (BPHZ) renormalization}")
        P.append(r"Each free constant at its canonical value $k_\tau = h(S'_-\tau)$ for a "
                 r"centered Gaussian noise.  Mean-zero parity makes trees with an odd number of "
                 f"noise vertices vanish ({nzero} of {len(rows)} below); the survivors are "
                 r"polynomials in the elementary expectations $h(\sigma)=\mathbb E[\Pi\sigma](0)$ "
                 r"(numeric $h(\sigma)$: Phase~4).")
        P.append(r"\begin{align*}")
        for t, k_free, v in rows:
            rhs = (r"0 && \text{(vanishes: odd noise parity)}" if v == 0
                   else flatex(v))
            P.append(f"  {sympy.latex(k_free)} &= {rhs} \\\\")
        P.append(r"\end{align*}")
        P.append(r"giving the canonically renormalized equation(s)")
        canon = {t: v for t, _k, v in rows}
        for a, (u, op, _r) in enumerate(eq.spde.equations):
            P.append(r"\begin{align*}")
            P.append(f"  {op_latex(op)}\\, {u.name} &= {flatex(eq.original_rhs(a))} \\\\")
            for ct in sorted(eq.per_component[a], key=lambda c: c.homogeneity._key()):
                v = canon.get(ct.tree)
                if v is not None and v != 0:        # one surviving counterterm per line
                    term = sympy.expand(sympy.Rational(1, ct.symmetry_factor)
                                        * v * eq._pretty(ct.elem_diff))
                    P.append(r"  &\quad + " + flatex(term) + r" \\")
            P.append(r"\end{align*}")
        if legend:
            P.append(r"where each surviving elementary expectation is $h(\sigma)$ for $\sigma$:")
            P.append(r"\begin{center}")
            for sym, tr in legend:
                note = (r"\quad{\footnotesize(contraction node, $o=" + hom_latex(tr.o) + "$)}"
                        if tr.color == "red" else "")
                P.append(f"${sympy.latex(sym)} = h\\bigl({{}}$ {forest(tr, sig)} $\\bigr)$"
                         f"{note}\\par\\bigskip")
            P.append(r"\end{center}")
    else:
        P.append(r"\bigskip")
        P.append(r"{\itshape Pass \texttt{canonical=True} for the canonical (BPHZ) "
                 r"renormalization $k_\tau=h(S'_-\tau)$ of a centered Gaussian noise "
                 r"(parity-reduced; many constants vanish). Numeric $h(\sigma)$ values need a "
                 r"noise law (Phase~4).}")
    P.append(r"\end{document}")
    return "\n".join(P)
