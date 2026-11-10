"""The HTML explorer: one self-contained page for a renormalized equation.

``html_document(eq)`` (= ``eq.render("html")``) returns a single HTML file that opens offline:
the renormalized family, a gallery of **every** divergent tree (drawn by `render.svg`, grouped
by homogeneity, most singular first) and one detail panel per tree.  Hovering a counterterm in
the equation lights up its tree and vice versa; clicking selects a tree.

Everything is rendered here, in Python; the inline script only toggles classes (and the page
reads fine without it — every panel simply shows).  No network: formulas are SymPy
presentation MathML (native in current browsers), trees are inline SVG, style and script are
inline.  This module only *reads* the engine's objects.  See notes/visual_layer_design.md.
"""
from __future__ import annotations

import re
from fractions import Fraction
from html import escape
from typing import TYPE_CHECKING
from xml.dom.minidom import Element

import sympy
from sympy.core.function import AppliedUndef
from sympy.printing.mathml import MathMLPresentationPrinter

from .report import _prime, _sorted_trees, const_map, operator_str
from .svg import svg
from .tree import shorthand

if TYPE_CHECKING:
    from ..core.homogeneity import Homogeneity
    from ..core.signature import Signature
    from ..renorm.equation import CanonicalFamily, Counterterm, RenormalizedEquation
    from ..trees.tree import DecoratedTree

_INERT = "F(τ*) = 0 — no counterterm"


# --------------------------------------------------------------------------- #
# formulas: presentation MathML, with f'(u) instead of SymPy's d/du fraction
# --------------------------------------------------------------------------- #

class _PrimeMathML(MathMLPresentationPrinter):
    def _print_Derivative(self, expr: sympy.Expr) -> Element:
        # same rule as the text/LaTeX reports (`report._prime`): a one-argument applied
        # function differentiated in its argument prints as f'(u), f''(u), f^(4)(u).
        text = _prime(expr, lambda v: "\x00", lambda n: f"^({n})")
        if text is None:
            return super()._print_Derivative(expr)
        head = text.split("(\x00)")[0]                       # e.g.  f''   or  f^(4)
        (var,) = expr.expr.args
        name = expr.expr.func.__name__
        row = self.dom.createElement("mrow")
        mi = self.dom.createElement("mi")                    # single-letter <mi> stays italic
        mi.appendChild(self.dom.createTextNode(name))
        row.appendChild(mi)
        for tok, child in ((head[len(name):].replace("'", "′"), None), ("(", None),
                           ("", self._print(var)), (")", None)):
            if child is None:
                mo = self.dom.createElement("mo")
                mo.appendChild(self.dom.createTextNode(tok))
                row.appendChild(mo)
            else:
                row.appendChild(child)
        return row


    def parenthesize(self, item: sympy.Expr, level: int, strict: bool = False) -> Element:
        # f'(u) is atomic — no brackets around it inside a product (SymPy would add them)
        if isinstance(item, sympy.Derivative) and _prime(item, str, str) is not None:
            return self._print(item)
        return super().parenthesize(item, level, strict)


def _escape_names(expr: sympy.Expr) -> sympy.Expr:
    # SymPy's MathML printer writes text nodes *unescaped* (so its own &dd; entities survive),
    # hence a user name such as "a<b" would break the page: escape names before printing.
    def bad(name: str) -> bool:
        return escape(name) != name

    subs = {s: sympy.Symbol(escape(s.name)) for s in expr.free_symbols if bad(s.name)}
    expr = expr.xreplace(subs)
    return expr.replace(lambda e: isinstance(e, AppliedUndef) and bad(e.func.__name__),
                        lambda e: sympy.Function(escape(e.func.__name__))(*e.args))


def mathml(expr: sympy.Expr) -> str:
    """``expr`` as inline presentation MathML (``<math>…</math>``)."""
    return "<math>" + _PrimeMathML().doprint(_escape_names(sympy.sympify(expr))) + "</math>"


def _op_html(op: object) -> str:
    """The operator's display symbol with ``∂_t`` typeset as a subscript."""
    return re.sub(r"_(\w)", r"<sub>\1</sub>", escape(operator_str(op)))


def _rate_html(h: Homogeneity) -> str:
    """`scheme.divergence_str` with a real superscript: ``ε^(−1)`` → ε⁻¹-style HTML."""
    return "log(1/ε)" if h.std == 0 else f"ε<sup>{escape(str(h.std)).replace('-', '−')}</sup>"


# --------------------------------------------------------------------------- #
# style and script (inline; the script only toggles classes)
# --------------------------------------------------------------------------- #

_CSS = """
:root{--bg:#fff;--fg:#1c1c1e;--mut:#6b6b70;--line:#d9d9de;--card:#f6f6f8;--hl:#fff3c4;--sel:#2c6fbb}
@media(prefers-color-scheme:dark){:root{--bg:#1c1c1e;--fg:#e6e6e8;--mut:#9a9aa0;--line:#3a3a3f;
--card:#26262a;--hl:#4a4320;--sel:#6ea8e6}}
*{box-sizing:border-box}
body{margin:0 auto;padding:24px 20px 60px;max-width:1040px;background:var(--bg);color:var(--fg);
font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
h1{font-size:20px;margin:0 0 4px}h2{font-size:13px;margin:28px 0 8px;color:var(--mut);
text-transform:uppercase;letter-spacing:.06em;font-weight:600}
.meta{color:var(--mut);font-size:13px;margin:0 0 14px}
.eq{padding:12px 14px;border:1px solid var(--line);border-radius:8px;margin:8px 0;line-height:2.1;
overflow-x:auto}
.nb{white-space:nowrap}
.ct{padding:3px 5px;border-radius:5px;cursor:pointer;white-space:nowrap}
.ct.hl,.card.hl{background:var(--hl)}
.canon{color:var(--mut);font-size:14px}
.gallery{display:flex;flex-wrap:wrap;gap:14px 22px;align-items:flex-end}
.grade{display:flex;flex-direction:column;gap:6px}
.grade>.hom{font-size:12px;color:var(--mut);border-bottom:1px solid var(--line);padding-bottom:2px}
.grade>.row{display:flex;gap:8px;align-items:flex-end;flex-wrap:wrap}
.card{margin:0;padding:8px 8px 4px;border:1px solid var(--line);border-radius:8px;
background:var(--card);cursor:pointer;text-align:center;min-width:52px;
display:flex;flex-direction:column;justify-content:flex-end;align-items:center;gap:4px}
.card .k{font-size:12px}.card.inert{opacity:.5}.card.inert .k{font-size:10px;max-width:86px}
.card.sel{border-color:var(--sel);box-shadow:0 0 0 1px var(--sel)}
.detail{border:1px solid var(--line);border-radius:8px;padding:16px;margin:12px 0}
.js .detail{display:none}.js .detail.sel{display:block}
.detail .top{display:flex;gap:28px;align-items:center;flex-wrap:wrap}
.facts{display:grid;grid-template-columns:max-content 1fr;gap:3px 14px;font-size:14px;margin:0}
.facts dt{color:var(--mut)}.facts dd{margin:0}
.sh{font-size:12px;color:var(--mut);margin-top:4px;text-align:center}
.coproduct,.canonical{margin-top:16px;border-top:1px solid var(--line);padding-top:10px}
.terms{display:flex;flex-wrap:wrap;gap:10px 6px;align-items:center}
.term{display:inline-flex;align-items:center;gap:5px;padding:4px 8px;border:1px solid var(--line);
border-radius:6px}
.op{color:var(--mut);padding:0 2px}
.pairs{display:flex;flex-wrap:wrap;gap:6px 16px;align-items:center;margin:4px 0 10px}
h3{font-size:13px;font-weight:600;margin:0 0 8px;color:var(--mut)}
details{margin-top:8px}summary{cursor:pointer;color:var(--mut);font-size:13px}
.note{color:var(--mut);font-size:13px}
"""

_JS = """
(function(){
  document.body.classList.add('js');
  var cards=[].slice.call(document.querySelectorAll('.card'));
  function all(i){return document.querySelectorAll('[data-tree="'+i+'"]');}
  function select(i){
    cards.forEach(function(c){c.classList.toggle('sel',c.dataset.tree===i);});
    [].forEach.call(document.querySelectorAll('.detail'),function(d){
      d.classList.toggle('sel',d.id==='tree-'+i);});
  }
  [].forEach.call(document.querySelectorAll('.card,.ct'),function(el){
    var i=el.dataset.tree;
    el.addEventListener('mouseenter',function(){[].forEach.call(all(i),function(x){x.classList.add('hl');});});
    el.addEventListener('mouseleave',function(){[].forEach.call(all(i),function(x){x.classList.remove('hl');});});
    el.addEventListener('click',function(){select(i);});
  });
  document.addEventListener('keydown',function(e){
    if(e.key!=='ArrowLeft'&&e.key!=='ArrowRight')return;
    var k=cards.findIndex(function(c){return c.classList.contains('sel');});
    k=Math.max(0,Math.min(cards.length-1,k+(e.key==='ArrowRight'?1:-1)));
    select(cards[k].dataset.tree);e.preventDefault();
  });
  if(cards.length)select(cards[0].dataset.tree);
})();
"""


# --------------------------------------------------------------------------- #
# the page
# --------------------------------------------------------------------------- #

def _detail_extras(eq: RenormalizedEquation, t: DecoratedTree,
                   fam: CanonicalFamily | None) -> str:
    """Further sections of a tree's detail panel (the coproduct, the canonical value)."""
    canon = _canonical_html(eq, t, fam) if fam is not None and t in fam.constants else ""
    return canon + _coproduct_html(t, eq.sig)


def _canonical_html(eq: RenormalizedEquation, t: DecoratedTree, fam: CanonicalFamily) -> str:
    """The canonical (BPHZ) value of ``k_τ`` and, for each elementary expectation ``h_i`` in
    it, the tree ``σ_i`` with one diagram per Wick pairing (`scheme.pairings`) — the terms of
    the integral ``h_ε(σ_i)``.  Non-bare ``σ_i`` have no such integral here and say so."""
    from ..renorm.scheme import is_bare, pairings

    sig, value, k = eq.sig, fam.constants[t], fam.free_constants[t]
    head = ("reduced, symmetric noise" if fam.reduced and fam.symmetric
            else "reduced" if fam.reduced else "centered Gaussian noise")
    out = [f'<div class="canonical"><h3>Canonical (BPHZ) value — {head}</h3>'
           f"<p>{mathml(sympy.Eq(k, value, evaluate=False))}"
           + (' <span class="note">— vanishes canonically</span>' if value == 0 else "") + "</p>"]
    for sym, sigma in fam.legend.items():
        if sym not in value.free_symbols:
            continue
        out.append(f'<div class="pairs"><span>{mathml(sym)} = h<sub>ε</sub>( {svg(sigma, sig)} )'
                   "</span>")
        if is_bare(sigma):
            ps = pairings(sigma, sig)
            out.append(f'<span class="note">= sum over {len(ps)} Wick pairing'
                       f'{"s" if len(ps) != 1 else ""}:</span>'
                       + "".join(svg(sigma, sig, pairing=m) for m in ps))
        else:
            out.append('<span class="note">left symbolic — needs the full canonical model '
                       "(polynomial decoration Xⁿ)</span>")
        out.append("</div>")
    out.append("</div>")
    return "".join(out)


_SHOWN_TERMS = 12


def _coproduct_html(t: DecoratedTree, sig: Signature) -> str:
    """The extraction–contraction coaction ``δτ = Σ c · (extracted forest) ⊗ (contracted tree)``
    (`trees.coproducts.delta_minus`, tex 5636), every term drawn.  A term carries no unique
    provenance on τ — different cuts giving the same pair are summed — so the pair itself is
    what is shown.  Terms beyond the first `_SHOWN_TERMS` fold into a ``<details>``."""
    from ..trees.coproducts import delta_minus

    def term(forest: tuple[DecoratedTree, ...], right: DecoratedTree, c: Fraction) -> str:
        coeff = "" if c == 1 else ("−" if c == -1 else escape(str(c)).replace("-", "−") + " ·")
        left = "".join(svg(f, sig) for f in forest) or "𝟙"
        return (f'<span class="term">{coeff}<span>{left}</span><span class="op">⊗</span>'
                f"{svg(right, sig)}</span>")

    items = sorted(delta_minus(t, sig).items(),
                   key=lambda kv: (len(kv[0][0]), sum(f.nodes() for f in kv[0][0]),
                                   kv[0][1].nodes()))
    terms = [term(forest, right, c) for (forest, right), c in items]
    shown, rest = terms[:_SHOWN_TERMS], terms[_SHOWN_TERMS:]
    more = (f'<details><summary>{len(rest)} more terms</summary><div class="terms">'
            + "".join(rest) + "</div></details>") if rest else ""
    return ('<div class="coproduct"><h3>Extraction coproduct  δτ = Σ (extracted forest) ⊗ '
            f'(contracted tree) — {len(terms)} terms</h3><div class="terms">'
            + "".join(shown) + "</div>" + more + "</div>")


def html_document(eq: RenormalizedEquation, canonical: bool = False,
                  reduced: bool = False, symmetric: bool = True) -> str:
    """The explorer page for ``eq`` (see the module docstring).  ``canonical``/``reduced``
    add the canonical (BPHZ) values — and, like every other format, raise for a
    da Prato–Debussche-lifted equation, whose noises are not Gaussian."""
    sig = eq.sig
    fam = eq.canonical(reduced=reduced, symmetric=symmetric) if (canonical or reduced) else None
    trees = _sorted_trees(eq)
    index = {t: i for i, t in enumerate(trees)}
    kmap = const_map(eq)
    by_tree: dict[DecoratedTree, list[tuple[int, Counterterm]]] = {}
    for a, cts in eq.per_component.items():
        for ct in cts:
            by_tree.setdefault(ct.tree, []).append((a, ct))

    names = ", ".join(u.name for u in eq.unknowns)
    out = ['<!doctype html>\n<html lang="en"><head><meta charset="utf-8">'
           '<meta name="viewport" content="width=device-width,initial-scale=1">'
           f"<title>Renormalized equation — {escape(names)}</title>"
           f"<style>{_CSS}</style></head><body>"]

    # ---- header
    noises = "; ".join(f"{escape(n.name)} ∈ C<sup>{escape(str(n.homogeneity))}</sup>" for n in eq.spde.noises)
    n_ct = len(kmap)
    out.append(f"<h1>Renormalized equation</h1><p class=\"meta\">d = {sig.dim} · noise "
               f"{noises} · {len(trees)} divergent trees, {n_ct} counterterms"
               + (f", {len(trees) - n_ct} with F(τ*) = 0" if len(trees) != n_ct else "") + "</p>")

    # ---- the family: free constants, one linked summand per (tree, component)
    for a, (u, op, _rhs) in enumerate(eq.spde.equations):
        lhs = f"{_op_html(op)} {mathml(sympy.Symbol(u.name))} ="
        parts = [f"{lhs} {mathml(eq.original_rhs(a))}"]
        for ct in sorted(eq.per_component.get(a, []), key=lambda c: index[c.tree]):
            parts.append(f'<span class="nb"><span class="op">+</span><span class="ct" '
                         f'data-tree="{index[ct.tree]}" data-comp="{a}">'
                         f"{mathml(eq._pretty(ct.term))}</span></span>")
        out.append('<div class="eq">' + " ".join(parts) + "</div>")
        if fam is not None:
            out.append(f'<p class="canon">canonically renormalized: {lhs} '
                       f"{mathml(eq._pretty(fam.rhs(a)))}</p>")
    if not fam and sig.non_gaussian:
        out.append('<p class="note">Canonical (BPHZ) values are not available for this equation: '
                   f"{escape(', '.join(sorted(sig.non_gaussian)))} are Wick powers of one "
                   "underlying noise, not independent Gaussians.</p>")

    # ---- gallery, grouped by homogeneity (most singular first)
    out.append('<h2>Divergent trees</h2><div class="gallery">')
    current = None
    for i, t in enumerate(trees):
        hom = str(t.homogeneity(sig))
        if hom != current:
            if current is not None:
                out.append("</div></div>")
            out.append(f'<div class="grade" data-hom="{escape(hom)}"><div class="hom">|τ| = '
                       f'{escape(hom)}</div><div class="row">')
            current = hom
        inert = t not in kmap
        label = _INERT if inert else mathml(kmap[t])
        out.append(f'<article class="card{" inert" if inert else ""}" data-tree="{i}">'
                   f'{svg(t, sig)}<div class="k">{label}</div></article>')
    if current is not None:
        out.append("</div></div>")
    out.append("</div>")

    # ---- one detail panel per tree
    out.append("<h2>Selected tree</h2>")
    for i, t in enumerate(trees):
        facts = [("|τ|", escape(str(t.homogeneity(sig)))), ("S(τ)", str(t.symmetry_factor()))]
        if t in kmap:
            facts.append(("constant", mathml(kmap[t])))
            facts.append(("diverges like", _rate_html(t.homogeneity(sig))))
            for a, ct in by_tree[t]:
                where = f" in the {escape(eq.unknowns[a].name)}-equation" if sig.n_components > 1 else ""
                facts.append((f"F(τ*){where}", mathml(eq._pretty(ct.elem_diff))))
        else:
            facts.append(("", _INERT))
        dl = "".join(f"<dt>{k}</dt><dd>{v}</dd>" for k, v in facts)
        out.append(f'<section class="detail" id="tree-{i}"><div class="top"><div>'
                   f'{svg(t, sig, scale=2)}<div class="sh">{escape(shorthand(t, sig))}</div></div>'
                   f'<dl class="facts">{dl}</dl></div>{_detail_extras(eq, t, fam)}</section>')

    out.append(f"<script>{_JS}</script></body></html>\n")
    return "".join(out)
