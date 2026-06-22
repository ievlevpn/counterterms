"""Assemble the full renormalized-equation report (notes/output.md §6).

``render(eq, fmt)`` dispatches to text / markdown / json / latex.  Every formula
is printed by SymPy; we own only the tree drawings (``render/tree.py``) and this
assembler.  Phase-3/4 information (coproducts, canonical constant values) shows as
a labelled placeholder until those phases land.
"""
from __future__ import annotations

import json

import sympy
from sympy.core.function import AppliedUndef
from sympy.printing.latex import LatexPrinter
from sympy.printing.str import StrPrinter

from .tree import ascii_art, coord_names, edge_glyph_text, shorthand, _node_glyph


# --------------------------------------------------------------------------- #
# prime notation:  Derivative(f(u), u)  →  f'(u)   (f''(u), f'''(u), f^{(4)}(u))
# Only single-argument applied functions differentiated in that argument; partial
# derivatives of multi-arg g(u,∂u) fall back to SymPy's default ∂-notation.
# --------------------------------------------------------------------------- #

def _prime(expr, base_print, sup):
    f = expr.expr
    if isinstance(f, AppliedUndef) and len(f.args) == 1:
        (var,) = f.args
        if all(v == var for v, _ in expr.variable_count):
            n = sum(int(c) for _, c in expr.variable_count)
            mark = "'" * n if n <= 3 else sup(n)
            return f.func.__name__ + mark + "(" + base_print(var) + ")"
    return None


class _PrimeLatex(LatexPrinter):
    def _print_Derivative(self, expr):
        return (_prime(expr, self._print, lambda n: f"^{{({n})}}")
                or super()._print_Derivative(expr))


class _PrimeStr(StrPrinter):
    def _print_Derivative(self, expr):
        return (_prime(expr, self._print, lambda n: f"^({n})")
                or super()._print_Derivative(expr))


def flatex(expr) -> str:
    return _PrimeLatex().doprint(expr)


def fstr(expr) -> str:
    return _PrimeStr().doprint(expr)


def render(eq, fmt: str = "text", canonical: bool = False) -> str:
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

def operator_str(op) -> str:
    s = "∂_t − Δ"
    if getattr(op, "mass", 0):
        s += f" + {op.mass}"
    return f"({s})"


def op_latex(op) -> str:
    s = r"\partial_t - \Delta"
    if getattr(op, "mass", 0):
        s += f" + {op.mass}"
    return f"({s})"


def hom_latex(h) -> str:
    return str(h).replace("κ", r"\kappa")


def const_map(eq) -> dict:
    """tree -> free constant k_τ (constants are shared across components)."""
    m = {}
    for cts in eq.per_component.values():
        for ct in cts:
            m[ct.tree] = ct.constant
    return m


def elem_map(eq, comp: int = 0) -> dict:
    return {ct.tree: ct.elem_diff for ct in eq.per_component.get(comp, [])}


def _sorted_trees(eq):
    return sorted(eq.all_trees, key=lambda t: t.homogeneity(eq.sig)._key())


# --------------------------------------------------------------------------- #
# Phase-3 bridge: the canonical (BHZ) renormalization, symbolic in h(σ).
# --------------------------------------------------------------------------- #

def canonical_data(eq):
    """Build the renormalization structure and return
    ``(rows, legend)`` where ``rows = [(tree, free_const k_τ, bhz_expr)]`` over the
    counterterm trees (homogeneity order) and ``legend = [(h_symbol, σ)]`` for every
    elementary expectation ``h(σ)`` appearing.  Heavy for deep trees (KPZ)."""
    from ..structures import build_renormalization
    rs = build_renormalization(eq.spde)
    cmap = const_map(eq)
    rows = [(t, cmap[t], rs.bhz_character(t)) for t in _sorted_trees(eq) if t in cmap]
    legend = sorted(((sym, tr) for tr, sym in rs._h.items()),
                    key=lambda x: int(str(x[0])[1:]))
    return rows, legend


def _h_annotation(tr) -> str:
    return f"  [contracted node, o={tr.o}]" if tr.color == "red" else ""


# --------------------------------------------------------------------------- #
# text
# --------------------------------------------------------------------------- #

def _sec(title: str) -> str:
    return f"\n── {title} " + "─" * max(2, 64 - len(title))


def text_report(eq, canonical: bool = False) -> str:
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
        out.append(_sec("CANONICAL (BHZ) RENORMALIZATION  k_τ = h(S'₋τ)"))
        out.append("  Each free constant, set to its canonical (BPHZ) value — symbolic in the")
        out.append("  elementary expectations h(σ) = 𝔼[Π σ](0) (numeric values: Phase 4).")
        for t, k_free, k_bhz in rows:
            out.append(f"  {k_free} = {fstr(k_bhz)}      [τ={shorthand(t, sig, coords)}]")
        out.append("  where")
        for sym, tr in legend:
            out.append(f"    {sym} = h({shorthand(tr, sig, coords)}){_h_annotation(tr)}")
    else:
        out.append(_sec("CANONICAL (BHZ) RENORMALIZATION"))
        out.append("  k_τ = h(S'₋τ), symbolic in h(σ): pass canonical=True "
                   "(heavy for deep trees, e.g. KPZ).")
        out.append(_sec("CANONICAL VALUES (Phase 4 — needs a NoiseLaw)"))
        out.append("  numeric h(σ) = 𝔼[Π σ](0) via Wick — see ROADMAP O4.")
    return "\n".join(out)


# --------------------------------------------------------------------------- #
# markdown
# --------------------------------------------------------------------------- #

def markdown_report(eq, canonical: bool = False) -> str:
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
        row = f"| `{shorthand(t, sig, coords)}` | `{t.homogeneity(sig)}` | " \
              f"{t.symmetry_factor()} | `{kc}` |"
        if scalar:
            f = fstr(eq._pretty(emap[t])) if k is not None else "0"
            row += f" `{f}` |"
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
        L += ["", "## Canonical (BHZ) renormalization", "",
              "Each free constant at its canonical (BPHZ) value, symbolic in the elementary "
              "expectations `h(σ) = 𝔼[Π σ](0)` (numeric values: Phase 4).", ""]
        for t, k_free, k_bhz in rows:
            L.append(f"- `{k_free} = {fstr(k_bhz)}`   (τ = `{shorthand(t, sig, coords)}`)")
        L += ["", "where", ""]
        for sym, tr in legend:
            L.append(f"- `{sym} = h({shorthand(tr, sig, coords)})`{_h_annotation(tr)}")
    return "\n".join(L)


# --------------------------------------------------------------------------- #
# json
# --------------------------------------------------------------------------- #

def json_report(eq, canonical: bool = False) -> str:
    sig = eq.sig
    coords = coord_names(sig.width)
    cmap = const_map(eq)
    emap = elem_map(eq)
    scalar = sig.n_components == 1

    def hom_d(h):
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
        data["bhz"] = [{"tree": shorthand(t, sig, coords), "constant": str(kf),
                        "value": str(kb), "value_latex": flatex(kb)}
                       for t, kf, kb in rows]
        data["h_legend"] = [{"symbol": str(sym), "tree": shorthand(tr, sig, coords),
                             "contracted": tr.color == "red", "o": str(tr.o)}
                            for sym, tr in legend]
    return json.dumps(data, indent=2, ensure_ascii=False)
