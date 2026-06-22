"""Assemble the full renormalized-equation report (notes/output.md §6).

``render(eq, fmt)`` dispatches to text / markdown / json / latex.  Every formula
is printed by SymPy; we own only the tree drawings (``render/tree.py``) and this
assembler.  Phase-3/4 information (coproducts, canonical constant values) shows as
a labelled placeholder until those phases land.
"""
from __future__ import annotations

import json
from typing import TYPE_CHECKING, Callable

import sympy
from sympy.core.function import AppliedUndef
from sympy.printing.latex import LatexPrinter
from sympy.printing.str import StrPrinter

from ..renorm.scheme import WHITE_NOISE, Expectation, expectation, is_bare
from .tree import ascii_art, coord_names, edge_glyph_text, shorthand, _node_glyph

if TYPE_CHECKING:
    from ..core.homogeneity import Homogeneity
    from ..renorm.equation import RenormalizedEquation
    from ..trees.tree import DecoratedTree


# --------------------------------------------------------------------------- #
# prime notation:  Derivative(f(u), u)  →  f'(u)   (f''(u), f'''(u), f^{(4)}(u))
# Only single-argument applied functions differentiated in that argument; partial
# derivatives of multi-arg g(u,∂u) fall back to SymPy's default ∂-notation.
# --------------------------------------------------------------------------- #

def _prime(
    expr: sympy.Expr,
    base_print: Callable[[sympy.Expr], str],
    sup: Callable[[int], str],
) -> str | None:
    f = expr.expr
    if isinstance(f, AppliedUndef) and len(f.args) == 1:
        (var,) = f.args
        if all(v == var for v, _ in expr.variable_count):
            n = sum(int(c) for _, c in expr.variable_count)
            mark = "'" * n if n <= 3 else sup(n)
            return f.func.__name__ + mark + "(" + base_print(var) + ")"
    return None


class _PrimeLatex(LatexPrinter):
    def _print_Derivative(self, expr: sympy.Expr) -> str:
        return (_prime(expr, self._print, lambda n: f"^{{({n})}}")
                or super()._print_Derivative(expr))


class _PrimeStr(StrPrinter):
    def _print_Derivative(self, expr: sympy.Expr) -> str:
        return (_prime(expr, self._print, lambda n: f"^({n})")
                or super()._print_Derivative(expr))


def flatex(expr: sympy.Expr) -> str:
    return _PrimeLatex().doprint(expr)


def fstr(expr: sympy.Expr) -> str:
    return _PrimeStr().doprint(expr)


def render(eq: RenormalizedEquation, fmt: str = "text", canonical: bool = False) -> str:
    # `canonical` adds the Phase-3 BHZ section (k_τ = h(S'₋τ)); off by default
    # because the twisted antipode blows up for deep trees (e.g. KPZ).
    if fmt == "text":
        return text_report(eq, canonical)
    if fmt == "markdown":
        return markdown_report(eq, canonical)
    if fmt == "json":
        return json_report(eq, canonical)
    if fmt == "latex":
        from .latex import latex_document
        return latex_document(eq, canonical)
    raise ValueError(f"unknown render format {fmt!r} (text/markdown/json/latex)")


# --------------------------------------------------------------------------- #
# shared helpers
# --------------------------------------------------------------------------- #

def operator_str(op: object) -> str:
    # operators carry their own display symbol (Operator.symbol); fall back to L
    return f"({getattr(op, 'symbol', 'L')})"


def op_latex(op: object) -> str:
    return f"({getattr(op, 'latex', 'L')})"


def hom_latex(h: Homogeneity) -> str:
    return str(h).replace("κ", r"\kappa")


def const_map(eq: RenormalizedEquation) -> dict:
    """tree -> free constant k_τ (constants are shared across components)."""
    m = {}
    for cts in eq.per_component.values():
        for ct in cts:
            m[ct.tree] = ct.constant
    return m


def elem_map(eq: RenormalizedEquation, comp: int = 0) -> dict:
    return {ct.tree: ct.elem_diff for ct in eq.per_component.get(comp, [])}


def _sorted_trees(eq: RenormalizedEquation) -> list[DecoratedTree]:
    return sorted(eq.all_trees, key=lambda t: t.homogeneity(eq.sig)._key())


# --------------------------------------------------------------------------- #
# Phase-3 bridge: the canonical (BHZ) renormalization, symbolic in h(σ).
# --------------------------------------------------------------------------- #

def canonical_data(eq: RenormalizedEquation) -> tuple[list, list]:
    """Build the renormalization structure and return ``(rows, legend)``:

    * ``rows = [(tree, free_const k_τ, value)]`` over the counterterm trees (homogeneity
      order), where ``value`` is the **canonical (BPHZ) constant** ``k_τ = h(S'₋τ)`` for a
      centered-Gaussian noise — *parity-reduced*, so trees with an odd number of noise
      vertices vanish (``value == 0``) and the survivors collapse to short polynomials in
      the even-noise ``h``-symbols;
    * ``legend = [(h_symbol, σ)]`` for each surviving elementary expectation ``h(σ)``.

    The parity rule both sharpens the output (gKPZ: 3 of 5 vanish) and tames it
    (KPZ: ~90-char expressions instead of 144 ``h``-symbols)."""
    from ..structures import build_renormalization
    rs = build_renormalization(eq.spde)
    cmap = const_map(eq)
    rows = [(t, cmap[t], rs.canonical_character(t)) for t in _sorted_trees(eq) if t in cmap]
    # keep only the h-symbols that survive, and renumber them h0, h1, … contiguously
    used = sorted(set().union(*(v.free_symbols for _t, _k, v in rows)) if rows else set(),
                  key=lambda s: int(str(s)[1:]))
    relabel = {s: sympy.Symbol(f"h{i}") for i, s in enumerate(used)}
    rows = [(t, k, v.xreplace(relabel)) for t, k, v in rows]
    legend = sorted(((relabel[sym], tr) for tr, sym in rs._h.items() if sym in relabel),
                    key=lambda x: int(str(x[0])[1:]))
    return rows, legend


def _h_annotation(tr: DecoratedTree) -> str:
    return f"  [o={tr.o}]" if tr.color == "red" else ""


# --------------------------------------------------------------------------- #
# the elementary integrals h_ε(σ) = 𝔼[Π^ε σ](0) — the objects Track B2 must
# evaluate.  Already built (unevaluated) by `renorm.scheme.expectation`: a sum
# over Wick pairings of  ∏ ∂^p K (tree edges) · ∏ C_ε (matched noise pairs),
# integrated over the internal vertices (root at 0).  We only typeset them.  K is
# the (abstract) singular kernel of L⁻¹ — Hairer's K in K̄ = K + R, the
# diagonal-singular, compactly-supported part of the Green's function K̄ (it
# *explodes on the diagonal*, which is why these integrals diverge; tex 2105, 5683).
# C_ε is the regularized covariance.
# --------------------------------------------------------------------------- #

def _kernel_p(name: str) -> tuple[int, ...]:
    """``"K_0_1"`` → ``(0, 1)`` (the edge's derivative multi-index p)."""
    return tuple(int(x) for x in name[2:].split("_"))


def _factor_latex(f: sympy.Expr) -> str:
    """One integrand factor → LaTeX: ``K_p`` → ``∂^p K(·)``, covariance → ``C_ε(·)``."""
    if isinstance(f, sympy.Pow):
        return _factor_latex(f.base) + "^{" + sympy.latex(f.exp) + "}"
    name, arg = f.func.__name__, sympy.latex(f.args[0])
    if name.startswith("K_"):
        p = _kernel_p(name)
        k = "K" if not any(p) else r"\partial^{(" + ",".join(map(str, p)) + ")}K"
        return k + r"\!\left(" + arg + r"\right)"
    return r"C_\varepsilon\!\left(" + arg + r"\right)"         # regularized covariance


def integral_latex_block(sym: str, e: Expectation, width: int,
                         max_units: float = 6.5, per_line: int = 4) -> str:
    """LaTeX display for ``sym = h_ε(σ)``.  Short ⇒ one ``\\[…\\]`` line; long ⇒ a ``multline*``
    broken at factor boundaries (``per_line`` factors per line) so a big integrand (many
    kernel/covariance factors, or a sum over Wick pairings) does not overflow the text width.

    We break on an estimate of *rendered* width — each factor ``∂^pK``/``C_ε`` is ~1 unit, each
    integration variable in the measure ~0.4, and ``+`` between terms ~1 — not the LaTeX *source*
    length (which the verbose ``\\partial^{(0,1)}K\\!\\left(…\\right)`` markup would inflate)."""
    if e.is_zero:
        return r"\[" + sym + r" = 0\]"
    terms = []                                                # (domain, [factors], measure, nvar)
    for integrand, intvars in e.terms:
        facs = [_factor_latex(f) for f in sympy.Mul.make_args(integrand)]
        if intvars:
            dom = r"\int_{(\mathbb R^{%d})^{%d}}" % (width, len(intvars))
            meas = r"\," + r"\,".join(r"\mathrm d" + sympy.latex(v) for v in intvars)
        else:
            dom = meas = ""
        terms.append((dom, facs, meas, len(intvars)))
    units = sum(len(f) + 0.4 * nv for _d, f, _m, nv in terms) + 1.0 * (len(terms) - 1)
    if units <= max_units:
        oneline = " + ".join(" ".join([d, *f]).strip() + m for d, f, m, _ in terms)
        return r"\[" + sym + " = " + oneline + r"\]"
    lines = []                                                # break the integrand across lines
    for ti, (dom, facs, meas, _nv) in enumerate(terms):
        chunks = [facs[i:i + per_line] for i in range(0, len(facs), per_line)] or [[]]
        for ci, chunk in enumerate(chunks):
            seg = (sym + " = ") if (ti == 0 and ci == 0) else ("+ " if ci == 0 else "")
            if ci == 0 and dom:
                seg += dom + " "
            seg += " ".join(chunk)
            if ci == len(chunks) - 1:
                seg += meas
            lines.append(seg)
    return "\\begin{multline*}\n" + " \\\\\n".join(lines) + "\n\\end{multline*}"


def expectation_str(e: Expectation, width: int) -> str:
    """Plain-text counterpart of `expectation_latex` (for the txt/markdown reports)."""
    if e.is_zero:
        return "0"

    def fac(f: sympy.Expr) -> str:
        if isinstance(f, sympy.Pow):
            return fac(f.base) + f"^{f.exp}"
        name, arg = f.func.__name__, str(f.args[0])
        if name.startswith("K_"):
            p = _kernel_p(name)
            k = "K" if not any(p) else "∂^(" + ",".join(map(str, p)) + ")K"
            return k + f"({arg})"
        return f"C_ε({arg})"

    parts = []
    for integrand, intvars in e.terms:
        body = "·".join(fac(f) for f in sympy.Mul.make_args(integrand))
        parts.append(f"∫ {body} " + " ".join("d" + str(v) for v in intvars)
                     if intvars else body)
    return " + ".join(parts)


def canonical_family_rhs(eq: RenormalizedEquation, comp: int, rows: list):
    """The canonically renormalized RHS of component ``comp``: the original equation plus
    each surviving counterterm with its free constant replaced by the canonical (BPHZ)
    value ``(k_τ/S(τ))·F(τ*)``; vanishing (odd-parity) trees drop out."""
    canon = {t: v for t, _k, v in rows}
    expr = eq.original_rhs(comp)
    for ct in eq.per_component[comp]:
        v = canon.get(ct.tree, sympy.Integer(0))
        if v != 0:
            expr += sympy.Rational(1, ct.symmetry_factor) * v * eq._pretty(ct.elem_diff)
    return sympy.expand(expr)


# --------------------------------------------------------------------------- #
# text
# --------------------------------------------------------------------------- #

def _sec(title: str) -> str:
    return f"\n── {title} " + "─" * max(2, 64 - len(title))


def text_report(eq: RenormalizedEquation, canonical: bool = False) -> str:
    sig = eq.sig
    coords = coord_names(sig.width)
    cmap = const_map(eq)
    emap = elem_map(eq)
    trees = _sorted_trees(eq)
    out = ["RENORMALIZED EQUATION", "=" * 21]

    out.append(_sec("EQUATION"))
    for a, (u, op, _r) in enumerate(eq.spde.equations):
        out.append(f"  {operator_str(op)} {u.name} = {fstr(eq.original_rhs(a))}")

    out.append(_sec("DOMAIN & NOISE"))
    out.append(f"  d = {sig.dim}   scaling 𝔰 = {tuple(sig.scaling.weights)}   "
               f"spacetime 1+{sig.dim}")
    if any(o != 2 for o in sig.comp_order):
        out.append(f"  operator order(s) = {sig.comp_order}  "
                   "(⚠ proven theory is 2nd-order parabolic)")
    for nz in eq.spde.noises:
        out.append(f"  {nz.name} ∈ C^({nz.homogeneity})   (β = {nz.homogeneity})")
    out.append("  subcritical: yes — tree generation terminated (finite 𝓑_<0)")

    out.append(_sec("RULE (derived from the nonlinearity)"))
    for b in sig.node_types:
        edges = [f"{edge_glyph_text(c, p, sig, coords)}[·] "
                 f"{'any' if cap is None else f'≤{cap}'}"
                 for (c, p, cap) in sig.allowed[b]]
        out.append(f"  {_node_glyph(b, sig)} → " + (", ".join(edges) or "(leaf)"))

    ndrop = sum(1 for t in trees if t not in cmap)
    out.append(_sec(f"DIVERGENT TREES  (𝓑_<0: {len(trees)} trees, "
                    f"{len(cmap)} counterterms, {ndrop} dropped)"))
    for t in trees:
        k = cmap.get(t)
        tag = str(k) if k is not None else "F(τ*)=0 — no counterterm"
        out.append(f"  τ = {shorthand(t, sig, coords)}      "
                   f"|τ| = {t.homogeneity(sig)}   S = {t.symmetry_factor()}   [{tag}]")
        out += ["      " + ln for ln in ascii_art(t, sig).splitlines()]
        if k is not None and sig.n_components == 1:
            out.append(f"      F(τ*) = {fstr(eq._pretty(emap[t]))}")
        out.append("")

    out.append(_sec("RENORMALIZED FAMILY"))
    for a in range(eq.n_components):
        u, op = eq.unknowns[a], eq.spde.equations[a][1]
        out.append(f"  {operator_str(op)} {u.name} = {fstr(eq.original_rhs(a))}")
        for ct in sorted(eq.per_component[a], key=lambda c: c.homogeneity._key()):
            term = sympy.expand(ct.coefficient * eq._pretty(ct.elem_diff))
            out.append(f"      + {fstr(term)}      "
                       f"[τ={shorthand(ct.tree, sig, coords)}, |τ|={ct.homogeneity}]")
        out.append("")

    if canonical:
        rows, legend = canonical_data(eq)
        nzero = sum(1 for _t, _k, v in rows if v == 0)
        out.append(_sec("CANONICAL (BPHZ) RENORMALIZATION  k_τ = h(S'₋τ)"))
        out.append("  Each free constant at its canonical value for a centered Gaussian noise.")
        out.append(f"  Mean-zero parity ⇒ trees with odd noise count vanish "
                   f"({nzero}/{len(rows)} here); survivors are polynomials in the elementary")
        out.append("  expectations hᵢ ≡ h_ε(σ) = 𝔼[Π^ε σ](0) of the ε-regularized noise ξ_ε.")
        out.append("  These constants are ε-dependent and diverge as ε→0; the integrals below")
        out.append("  are left unevaluated.")
        for t, k_free, v in rows:
            rhs = "0   (vanishes — odd noise parity)" if v == 0 else fstr(v)
            out.append(f"  {k_free} = {rhs}      [τ={shorthand(t, sig, coords)}]")
        out.append("")
        out.append("  Canonically renormalized equation:")
        for a, (u, op, _r) in enumerate(eq.spde.equations):
            out.append(f"    {operator_str(op)} {u.name} = "
                       f"{fstr(canonical_family_rhs(eq, a, rows))}")
        if legend:
            out.append("  where, for the ε-regularized noise ξ_ε (covariance C_ε), each")
            out.append("  elementary expectation, with K the singular kernel of L⁻¹.  Those")
            out.append("  carrying a polynomial factor X^n are left symbolic:")
            for sym, tr in legend:
                out.append(f"    {sym} = h_ε({shorthand(tr, sig, coords)}){_h_annotation(tr)}")
                if is_bare(tr):
                    out.append(f"       = {expectation_str(expectation(tr, sig, WHITE_NOISE), sig.width)}")
                else:
                    out.append("       = (left symbolic — σ carries a polynomial factor X^n)")
    else:
        out.append(_sec("CANONICAL (BPHZ) RENORMALIZATION"))
        out.append("  k_τ = h(S'₋τ) for centered Gaussian noise (parity-reduced; many vanish):")
        out.append("  pass canonical=True.  The h(σ) values are left symbolic.")
    return "\n".join(out)


# --------------------------------------------------------------------------- #
# markdown
# --------------------------------------------------------------------------- #

def _md_cell(s: object) -> str:
    """Escape a value for a GFM table cell: a literal ``|`` ends the cell even inside
    backticks/`$`, so it must be backslash-escaped."""
    return str(s).replace("|", "\\|")


def markdown_report(eq: RenormalizedEquation, canonical: bool = False) -> str:
    sig = eq.sig
    coords = coord_names(sig.width)
    cmap = const_map(eq)
    emap = elem_map(eq)
    trees = _sorted_trees(eq)
    scalar = sig.n_components == 1
    L = ["# Renormalized equation", "", "## Equation", ""]
    for a, (u, op, _r) in enumerate(eq.spde.equations):
        L.append(f"`{operator_str(op)} {u.name} = {fstr(eq.original_rhs(a))}`")
        L.append("")

    L += ["## Domain & noise", ""]
    L.append(f"- `d = {sig.dim}`, scaling `𝔰 = {tuple(sig.scaling.weights)}`")
    for nz in eq.spde.noises:
        L.append(f"- noise `{nz.name} ∈ C^({nz.homogeneity})`")
    L.append("")

    L += ["## Divergent trees", ""]
    if scalar:
        L += ["| τ | \\|τ\\| | S(τ) | constant | F(τ*) |", "|---|---|---|---|---|"]
    else:
        L += ["| τ | \\|τ\\| | S(τ) | constant |", "|---|---|---|---|"]
    for t in trees:
        k = cmap.get(t)
        kc = str(k) if k is not None else "— (F=0)"
        row = f"| `{_md_cell(shorthand(t, sig, coords))}` | `{_md_cell(t.homogeneity(sig))}` | " \
              f"{t.symmetry_factor()} | `{_md_cell(kc)}` |"
        if scalar:
            f = fstr(eq._pretty(emap[t])) if k is not None else "0"
            row += f" `{_md_cell(f)}` |"
        L.append(row)

    L += ["", "## Renormalized family", ""]
    for a in range(eq.n_components):
        u, op = eq.unknowns[a], eq.spde.equations[a][1]
        L.append("```")
        L.append(f"{operator_str(op)} {u.name} = {fstr(eq.original_rhs(a))}")
        for ct in sorted(eq.per_component[a], key=lambda c: c.homogeneity._key()):
            term = sympy.expand(ct.coefficient * eq._pretty(ct.elem_diff))
            L.append(f"   + {fstr(term)}   [τ={shorthand(ct.tree, sig, coords)}]")
        L.append("```")

    L += ["", "## Tree drawings", ""]
    for t in trees:
        L.append("```")
        L.append(ascii_art(t, sig))
        L.append("```")

    if canonical:
        rows, legend = canonical_data(eq)
        nzero = sum(1 for _t, _k, v in rows if v == 0)
        L += ["", "## Canonical (BPHZ) renormalization", "",
              "Each free constant at its canonical value for a centered Gaussian noise. "
              f"Mean-zero parity ⇒ odd-noise trees vanish ({nzero}/{len(rows)} here); "
              "survivors are polynomials in the elementary expectations "
              "`hᵢ ≡ h_ε(σ) = 𝔼[Π^ε σ](0)` of the **ε-regularized noise ξ_ε**. These "
              "constants are ε-dependent and **diverge as ε→0**; the integrals below are "
              "left unevaluated.", ""]
        for t, k_free, v in rows:
            rhs = "0  *(vanishes — odd noise parity)*" if v == 0 else f"`{fstr(v)}`"
            L.append(f"- `{k_free}` = {rhs}   (τ = `{shorthand(t, sig, coords)}`)")
        L += ["", "Canonically renormalized equation:", "", "```"]
        for a, (u, op, _r) in enumerate(eq.spde.equations):
            L.append(f"{operator_str(op)} {u.name} = {fstr(canonical_family_rhs(eq, a, rows))}")
        L.append("```")
        if legend:
            L += ["", "where, for the ε-regularized noise `ξ_ε` (covariance `C_ε`) and the "
                  "singular kernel `K` of `L⁻¹`, each elementary expectation is given below; those "
                  "carrying a polynomial factor `X^n` are left symbolic:", ""]
            for sym, tr in legend:
                head = f"- `{sym} = h_ε({shorthand(tr, sig, coords)})`{_h_annotation(tr)}"
                if is_bare(tr):
                    L.append(head + f"  `= {expectation_str(expectation(tr, sig, WHITE_NOISE), sig.width)}`")
                else:
                    L.append(head + "  *(left symbolic — σ carries a polynomial factor `X^n`)*")
    return "\n".join(L)


# --------------------------------------------------------------------------- #
# json
# --------------------------------------------------------------------------- #

def json_report(eq: RenormalizedEquation, canonical: bool = False) -> str:
    sig = eq.sig
    coords = coord_names(sig.width)
    cmap = const_map(eq)
    emap = elem_map(eq)
    scalar = sig.n_components == 1

    def hom_d(h: Homogeneity) -> dict:
        return {"std": str(h.std), "kap": str(h.kap), "str": str(h)}

    trees = []
    for t in _sorted_trees(eq):
        k = cmap.get(t)
        trees.append({
            "shorthand": shorthand(t, sig, coords),
            "homogeneity": hom_d(t.homogeneity(sig)),
            "S": t.symmetry_factor(),
            "nodes": t.nodes(),
            "contributes": k is not None,
            "constant": str(k) if k is not None else None,
            "F_latex": (flatex(eq._pretty(emap[t]))
                        if (k is not None and scalar) else None),
        })

    data = {
        "equations": [{"unknown": u.name, "lhs": operator_str(op),
                       "rhs": fstr(eq.original_rhs(a))}
                      for a, (u, op, _r) in enumerate(eq.spde.equations)],
        "domain": {"dim": sig.dim, "scaling": list(sig.scaling.weights)},
        "noises": [{"name": nz.name, "homogeneity": hom_d(nz.homogeneity)}
                   for nz in eq.spde.noises],
        "rule": {b: [{"component": c, "p": list(p), "cap": cap}
                     for (c, p, cap) in sig.allowed[b]] for b in sig.node_types},
        "n_components": eq.n_components,
        "trees": trees,
        "family_latex": {eq.unknowns[a].name: flatex(eq.counterterm_rhs(a))
                         for a in range(eq.n_components)},
    }
    if canonical:
        rows, legend = canonical_data(eq)
        data["canonical_bphz"] = [{"tree": shorthand(t, sig, coords), "constant": str(kf),
                                   "value": str(v), "value_latex": flatex(v),
                                   "vanishes": v == 0}
                                  for t, kf, v in rows]
        data["h_legend"] = [{"symbol": str(sym), "tree": shorthand(tr, sig, coords),
                             "contracted": tr.color == "red", "o": str(tr.o)}
                            for sym, tr in legend]
        data["canonical_family_latex"] = {
            eq.unknowns[a].name: flatex(canonical_family_rhs(eq, a, rows))
            for a in range(eq.n_components)}
    return json.dumps(data, indent=2, ensure_ascii=False)
