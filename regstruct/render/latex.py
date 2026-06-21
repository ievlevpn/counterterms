"""Standalone, ``pdflatex``-ready LaTeX of the renormalized equation (notes/output.md §4–6).

Reuses ``render/tree.py``'s ``forest`` for the trees and ``sympy.latex`` for every
formula; we only wrap them in a document with the ``forest`` style preamble.
"""
from __future__ import annotations

import sympy

from .report import const_map, elem_map, hom_latex, op_latex, _sorted_trees
from .tree import coord_names, forest

_PREAMBLE = r"""\documentclass{article}
\usepackage{amsmath,amssymb}
\usepackage{forest}
\forestset{rstree/.style={for tree={grow'=north, parent anchor=north,
  child anchor=south, s sep=16pt, l sep=12pt, inner sep=1pt}}}
\begin{document}
\section*{Renormalized equation}"""


def latex_document(eq) -> str:
    sig = eq.sig
    cmap = const_map(eq)
    emap = elem_map(eq)
    scalar = sig.n_components == 1
    P = [_PREAMBLE]

    P.append(r"\begin{align*}")
    for a, (u, op, _r) in enumerate(eq.spde.equations):
        P.append(f"  {op_latex(op)}\\, {u.name} &= {sympy.latex(eq.original_rhs(a))} \\\\")
    P.append(r"\end{align*}")

    P.append(r"\noindent")
    for nz in eq.spde.noises:
        P.append(f"${sympy.latex(nz.symbol)} \\in C^{{{hom_latex(nz.homogeneity)}}}$\\quad")
    P.append("")

    # ponytail: forest-in-tabular is fragile; one centred block per tree instead.
    P.append(r"\subsection*{Divergent trees}")
    for t in _sorted_trees(eq):
        k = cmap.get(t)
        kc = sympy.latex(k) if k is not None else r"\text{(none)}"
        data = [f"|\\tau| = {hom_latex(t.homogeneity(sig))}",
                f"S(\\tau) = {t.symmetry_factor()}",
                f"k_\\tau = {kc}"]
        if scalar:
            fl = sympy.latex(eq._pretty(emap[t])) if k is not None else "0"
            data.append(f"F(\\tau^*) = {fl}")
        P.append(r"\begin{center}")
        P.append(forest(t, sig))
        P.append(r"\\[2pt] $" + r",\quad ".join(data) + r"$")
        P.append(r"\end{center}")

    P.append(r"\subsection*{Renormalized family}")
    for a in range(eq.n_components):
        u, op = eq.unknowns[a], eq.spde.equations[a][1]
        P.append(r"\begin{align*}")
        P.append(f"  {op_latex(op)}\\, {u.name} &= {sympy.latex(eq.original_rhs(a))} \\\\")
        for ct in sorted(eq.per_component[a], key=lambda c: c.homogeneity._key()):
            term = sympy.expand(ct.coefficient * eq._pretty(ct.elem_diff))
            P.append(r"  &\quad + " + sympy.latex(term) + r" \\")
        P.append(r"\end{align*}")

    P.append(r"\subsection*{Algebra (Phase 3) / canonical values (Phase 4)}")
    P.append(r"$\Delta\tau,\ \Delta^-\tau,\ S'_-\tau$, BHZ character and numeric "
             r"$c_\tau$ are not yet computed (see ROADMAP O3/O4).")
    P.append(r"\end{document}")
    return "\n".join(P)
