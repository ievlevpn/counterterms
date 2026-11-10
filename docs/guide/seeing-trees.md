# Seeing the trees

Two renderers show the objects as pictures rather than text: an SVG drawer for a single tree,
and a self-contained HTML explorer for a whole renormalized equation. Neither needs LaTeX, a
network connection or any dependency beyond SymPy.

## Trees as SVG

```python
from counterterms.render import svg

eq = spde.renormalize()            # the gKPZ example from "Getting started"
svg(eq.all_trees[0], eq.sig)       # a standalone "<svg …>…</svg>" string
```

The five gKPZ trees, most singular first:

<p><svg xmlns="http://www.w3.org/2000/svg" class="rs-tree" viewBox="-13.2 -13.2 26.4 26.4" width="37.0" height="37.0" style="vertical-align:middle"><circle class="n n-noise n-noise-0" data-id="0" cx="0.0" cy="0.0" r="4.2" fill="none" stroke="currentColor" stroke-width="1.1"/></svg> &nbsp; <svg xmlns="http://www.w3.org/2000/svg" class="rs-tree" viewBox="-13.2 -13.2 56.4 58.4" width="79.0" height="81.8" style="vertical-align:middle"><line class="e e-d" x1="13.7" y1="31.1" x2="1.7" y2="3.8" stroke="currentColor" stroke-width="1.1" stroke-dasharray="1.5 2.5" stroke-linecap="round"/><line class="e e-d" x1="16.3" y1="31.1" x2="28.3" y2="3.8" stroke="currentColor" stroke-width="1.1" stroke-dasharray="1.5 2.5" stroke-linecap="round"/><circle class="n n-bullet" data-id="0" cx="15.0" cy="34.0" r="3.2" fill="currentColor"/><circle class="n n-noise n-noise-0" data-id="1" cx="0.0" cy="0.0" r="4.2" fill="none" stroke="currentColor" stroke-width="1.1"/><circle class="n n-noise n-noise-0" data-id="2" cx="30.0" cy="0.0" r="4.2" fill="none" stroke="currentColor" stroke-width="1.1"/></svg> &nbsp; <svg xmlns="http://www.w3.org/2000/svg" class="rs-tree" viewBox="-13.2 -13.2 26.4 60.4" width="37.0" height="84.6" style="vertical-align:middle"><line class="e" x1="0.0" y1="29.8" x2="0.0" y2="4.2" stroke="currentColor" stroke-width="1.1"/><circle class="n n-noise n-noise-0" data-id="0" cx="0.0" cy="34.0" r="4.2" fill="none" stroke="currentColor" stroke-width="1.1"/><circle class="n n-noise n-noise-0" data-id="1" cx="0.0" cy="0.0" r="4.2" fill="none" stroke="currentColor" stroke-width="1.1"/></svg> &nbsp; <svg xmlns="http://www.w3.org/2000/svg" class="rs-tree" viewBox="-13.2 -19.0 37.0 32.2" width="51.8" height="45.1" style="vertical-align:middle"><circle class="n n-noise n-noise-0" data-id="0" cx="0.0" cy="0.0" r="4.2" fill="none" stroke="currentColor" stroke-width="1.1"/><text class="lbl" x="6.2" y="-3.0" font-size="8" text-anchor="start" fill="currentColor">Xₓ</text></svg> &nbsp; <svg xmlns="http://www.w3.org/2000/svg" class="rs-tree" viewBox="-13.2 -13.2 26.4 58.4" width="37.0" height="81.8" style="vertical-align:middle"><line class="e e-d" x1="0.0" y1="30.8" x2="0.0" y2="4.2" stroke="currentColor" stroke-width="1.1" stroke-dasharray="1.5 2.5" stroke-linecap="round"/><circle class="n n-bullet" data-id="0" cx="0.0" cy="34.0" r="3.2" fill="currentColor"/><circle class="n n-noise n-noise-0" data-id="1" cx="0.0" cy="0.0" r="4.2" fill="none" stroke="currentColor" stroke-width="1.1"/></svg></p>

The drawing follows the paper's convention (the same one the LaTeX `forest` output uses), with
the root at the bottom:

| Object | Drawing |
|---|---|
| integration node \(\bullet\) | small filled dot |
| noise \(\circ\) | open circle; a second and third noise type are filled grey and blue, further ones are labelled by name |
| kernel edge \(\mathcal{I}\) | solid line |
| derivative kernel \(\mathcal{I}_p\), \(p \neq 0\) | dotted line; the direction is labelled in \(d > 1\) |
| polynomial decoration \(X^n\) | small label beside the node |
| contraction node (extended decoration) | red square, labelled `o=…` |
| \(T^+\) root | blue dot |
| component of the kernel (systems) | `(c)` beside the edge |

Black elements use the surrounding text colour, so the pictures follow a dark theme.

## In a notebook

Leave any of these as the last expression of a cell and it draws itself:

| Expression | Shows |
|---|---|
| a tree `t` | its picture |
| `eq` (a `RenormalizedEquation`) | the family and every divergent tree with its constant |
| `eq.canonical()`, `eq.bphz()` | each tree beside its canonical constant, then the trees \(\sigma_i\) of the \(h_i\) |
| `pg.delta_minus(t)`, `pg.delta_plus(t)`, `pg.antipode(t)` | the sum, term by term, as drawn trees |
| `rn.coaction(t)`, `rn.coproduct(t)`, `rn.twisted_antipode(t)`, `rs.recentering(t)`, `rs.structure_coproduct(t)` | the same |
| `eq.explore()`, `eq.explore(reduced=True)` | the interactive explorer below, embedded in the cell |
| `eq.show(t)`, `pg.show(t)` | `t` drawn with that equation's signature |

The coproduct and antipode results are still ordinary dicts, equal to what the functions in
`counterterms.trees.coproducts` return; they only gained a display. `print(...)` gives the same
sum as text, for example `Ξ ⊗ ●·𝓘ₓ[▪(o=-3/2 - κ)]·𝓘ₓ[Ξ]`. A sum longer than 40 terms draws
the first 40 and says how many are left.

Use `eq.show(t)` rather than a bare `t` when the equation has several noises or components. A
tree on its own does not know them, so `Ξ_eta` alone draws exactly like `Ξ_xi`.

A notebook saved and reopened without being run may show these outputs blank, because
notebooks strip drawings from output they do not trust; re-run the cell.

`svg(t, sig, pairing=((i, j), …))` adds a dashed arc between noise nodes `i` and `j`. Node
numbers are the pre-order index over the tree's children (root = 0); they are the same
numbers as the integration variables \(z_i\) in `renorm.scheme.expectation`, and
`renorm.scheme.pairings(t, sig)` lists every Wick pairing in that form:

```python
from counterterms import playground
from counterterms.renorm.scheme import pairings

pg = playground("-3/2 - kappa")
t = pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"))     # ●·𝓘ₓ[Ξ]²
pairings(t, pg.sig)                                  # [((1, 2),)]
svg(t, pg.sig, pairing=((1, 2),))
```

<p><svg xmlns="http://www.w3.org/2000/svg" class="rs-tree" viewBox="-13.2 -28.5 56.4 73.7" width="79.0" height="103.2" style="vertical-align:middle"><line class="e e-d" x1="13.7" y1="31.1" x2="1.7" y2="3.8" stroke="currentColor" stroke-width="1.1" stroke-dasharray="1.5 2.5" stroke-linecap="round"/><line class="e e-d" x1="16.3" y1="31.1" x2="28.3" y2="3.8" stroke="currentColor" stroke-width="1.1" stroke-dasharray="1.5 2.5" stroke-linecap="round"/><circle class="n n-bullet" data-id="0" cx="15.0" cy="34.0" r="3.2" fill="currentColor"/><circle class="n n-noise n-noise-0" data-id="1" cx="0.0" cy="0.0" r="4.2" fill="none" stroke="currentColor" stroke-width="1.1"/><circle class="n n-noise n-noise-0" data-id="2" cx="30.0" cy="0.0" r="4.2" fill="none" stroke="currentColor" stroke-width="1.1"/><path class="pair" d="M 0.0 -5.2 Q 15.0 -20.5 30.0 -5.2" fill="none" stroke="#c0392b" stroke-width="1" stroke-dasharray="3 2"/></svg></p>

## The explorer

```python
eq.render("html")                                    # the page, as a string
eq.save("kpz", formats=("html",), pdf=False)         # writes output/kpz.html
eq.save("kpz", formats=("html",), pdf=False, reduced=True)   # with the canonical values
```

It is one file that opens offline in a browser. [Open the KPZ example](../kpz_explorer.html)
(reduced view).

- **The equation** is the family with free constants. Each counterterm is linked to its tree:
  hover either one and both light up; click to select.
- **The gallery** shows *every* divergent tree, grouped by homogeneity with the most singular
  first. A tree whose \(F(\tau^*)\) vanishes is dimmed and marked "no counterterm".
- **The detail panel** of the selected tree gives \(|\tau|\), \(S(\tau)\), its constant,
  the divergence rate and \(F(\tau^*)\), then the extraction coproduct \(\delta\tau\)
  drawn term by term as *(extracted forest) ⊗ (contracted tree)*.
- **With `canonical=True` or `reduced=True`** the panel also shows the canonical value of
  \(k_\tau\) and, for each elementary expectation \(h_\varepsilon(\sigma)\) in it, one
  diagram per Wick pairing — the terms of its integral.
- Left and right arrow keys step through the trees. With scripts disabled every panel is
  simply shown.

## Limits

- `"html"` is not among `save`'s default formats; ask for it explicitly.
- A coproduct term is shown as the pair it is. The extracted part is not highlighted on the
  original tree, because different cuts that give the same pair are summed into one term.
- Only the first 12 coproduct terms are shown inline; the rest unfold on click.
- Non-bare \(\sigma\) (with a polynomial decoration) have no pairing diagrams — their
  expectation needs the full canonical model.
- For a [da Prato–Debussche-lifted](daprato.md) equation the page renders without canonical
  values and says why; asking for them raises, as in every other format.
- Formulas are MathML, which current Chrome, Firefox and Safari render natively.
