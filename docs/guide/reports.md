# Reports & output

A `RenormalizedEquation` renders itself in four formats, with three orthogonal views.

## Quick access

```python
eq = spde.renormalize()

print(eq.summary())                  # terse: one line per counterterm
print(eq.report())                   # full text report
print(eq.render("markdown"))         # "text" | "markdown" | "json" | "latex"
eq.save("kpz", outdir="output")      # writes kpz.txt, kpz.md, kpz.json, kpz.tex (+ .pdf if latexmk is installed)
```

Every format contains: the parsed equation, each divergent tree \(\tau\) **drawn** (in the
paper’s convention — \(\circ\) noise, \(\bullet\) integration node, dotted edge = derivative
kernel) with its homogeneity \(|\tau|\), symmetry factor \(S(\tau)\), free constant \(k_\tau\)
and elementary differential \(F(\tau^*)\), and the assembled renormalized family.

Programmatic access, bypassing the renderer:

```python
eq.counterterms          # flat list of Counterterm(tree, homogeneity, symmetry_factor, elem_diff, constant)
eq.per_component[0]      # the same, per equation of a system
eq.counterterm_rhs(0)    # the assembled Σ k_τ/S(τ)·F(τ*) as one SymPy expression
eq.all_trees             # the raw divergent-tree tuple, before Υ-zero trees are dropped
```

## Canonical and reduced views

Three keyword flags on `report` / `render` / `save` / `to_json` control how much of the
BPHZ theory is folded into the constants:

```python
eq.report()                                  # free constants k_τ only (the default)
eq.report(canonical=True)                    # + the BPHZ section: k_τ = h(S'_- τ)
eq.report(reduced=True)                      # + all exact identities folded in
eq.report(reduced=True, symmetric=False)     # …for a noise NOT symmetric under x → −x
```

**Default (free constants).** The family with free \(k_\tau\) — correct for *every*
renormalization character, no probabilistic input assumed. This is the theorem-level output.

**`canonical=True`.** Adds, per tree, the canonical (BPHZ) constant as an exact twisted-antipode
combination \(k_\tau = h(S'_-\tau)\) of elementary-expectation symbols \(h(\sigma)\), each
\(h(\sigma)\) spelled out as an \(\varepsilon\)-regularized Wick integral. Constants that
provably vanish (Gaussian parity, root \(X^n\), pure-kernel total derivative) or duplicate
another are **marked but left in place** — the display stays valid for any centered Gaussian
noise. (This section is off by default because the twisted antipode grows quickly on deep trees —
KPZ at \(\beta_0 = -\tfrac32\) produces \(\sim\)1400-term antipode forests.)

**`reduced=True`** (implies the canonical section). The same constants with the exact identities
**substituted**: provable zeros dropped, duplicates merged — plus, when `symmetric=True`, the
**spatial-reflection identity** (odd total spatial-derivative order on the kernels ⇒
\(h(\sigma) = 0\)). For KPZ this collapses the canonical family to Hairer’s single diverging
constant: \(\partial_t u = \Delta u + (\partial_x u)^2 - C + \xi\).

!!! warning "`symmetric` is an assumption about *your* noise"
    The reflection identity is valid only for a noise whose law is invariant under
    \(x \to -x\) (white noise, symmetric mollifications). It defaults to `True`. If your noise
    is anisotropic, pass `symmetric=False`: the drift-type counterterms then genuinely survive,
    and the report says so explicitly (`reduction_assumes_symmetric_noise: false` in the JSON).
    The other reductions (parity, root \(X^n\), total derivative) are noise-independent and are
    applied in the reduced view regardless.

Side-by-side example (KPZ at \(\beta_0 = -\tfrac32 - \kappa\)):
[canonical view (PDF)](../kpz_canonical.pdf) — 8 constants, zeros marked;
[reduced view (PDF)](../kpz_reduced.pdf) — one constant, Hairer’s \(C_\varepsilon\).

## The JSON format

`eq.to_json(...)` / `eq.save(...)` emit a machine-readable document: the parsed equation,
per-component counterterms (tree as a nested dict, homogeneity as exact
\((\text{std}, \text{kap})\) rationals, \(S\), \(F(\tau^*)\) as a SymPy string, constant name),
and — with `canonical`/`reduced` — the character polynomials and the reduction flags. For the
full algebraic structure (coproducts, antipode, group), see
[`structure_json`](algebra.md#json-export) instead.

## Drawing individual trees

```python
from counterterms.render import shorthand, ascii_art, forest

shorthand(t, sig)    # one-line:  ●(I'[∘])²  style
ascii_art(t, sig)    # multi-line terminal drawing
forest(t, sig)       # LaTeX (forest package) code, as used in the PDF reports
```

!!! note "Cosmetic nondeterminism"
    The RULE line ordering in text/LaTeX reports depends on Python’s hash seed; content is
    otherwise identical run to run. Diff JSON, not text, if you need stable comparisons.
