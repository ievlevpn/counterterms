# API reference

Everything public is importable from the top-level `counterterms` package (the DSL, the five
pipeline entry points) or from the documented submodules. Full type annotations throughout; the
docstrings cite the paper by tex line.

## Top level

```python
from counterterms import (
    SPDE, Unknown, Noise, Operator, Parabolic, FractionalHeat, kappa, jet,
    renormalize, daprato_lift,
    build_renormalization, build_regularity_structure, build_renormalization_group,
    playground, Playground,
)
```

| Callable | Returns | Purpose |
|---|---|---|
| `renormalize(spde)` (= `spde.renormalize()`) | `RenormalizedEquation` | the family of renormalized equations |
| `daprato_lift(spde)` | `SPDE` | da Prato–Debussche change of variables for supercritical polynomial additive-noise equations |
| `build_renormalization(spde)` | `RenormalizationStructure` | \(\delta, \delta^-, S'_-\), BHZ / canonical characters |
| `build_regularity_structure(spde, gamma=1)` | `RegularityStructure` | \(\gamma\)-truncated \((T, T^+)\), graded basis, \(\Delta, \Delta^+\) |
| `build_renormalization_group(spde)` | `RenormalizationGroup` | \(G^-\): characters, convolution, inverse |

## The DSL (`counterterms.equation.dsl`)

| Class | Signature | Notes |
|---|---|---|
| `Unknown` | `Unknown(name, dim)` | exposes `.field`, `.t`, `.x`, `.coords` |
| `Noise` | `Noise(name, regularity, chaos=1)` | `regularity` **is** \(\beta_0\) (may carry `- kappa`; a string such as `"-3/2 - kappa"` is accepted); exposes `.symbol`, `.homogeneity`. `chaos=k>1` marks a Wick power \(:\!X^k\!:\) (set by `daprato_lift`) |
| `Noise.white` | `Noise.white(name, dim, order=2)` | spacetime white noise, \(\beta_0 = -(\text{order}+\text{dim})/2 - \kappa\) |
| `Noise.spatial_white` | `Noise.spatial_white(name, dim)` | time-independent white noise, \(\beta_0 = -\text{dim}/2 - \kappa\) |
| `Parabolic` | `Parabolic(dim, mass=0, order=2)` | \(\partial_t - \Delta\ (+\,m)\) |
| `FractionalHeat` | `FractionalHeat(dim, sigma)` | \(\partial_t + (-\Delta)^\sigma\); order \(2\sigma\) |
| `Operator` | `Operator(dim, scaling, order, label, symbol, latex)` | base class — the engine reads only `(scaling, order)` |
| `SPDE` | `SPDE(noises, operator=…, unknown=…, rhs=…)` or `SPDE(noises, equations=[(u, op, rhs), …])` | scalar or system |
| `kappa` | SymPy symbol | the positive infinitesimal in regularities |

Lower-level (used by the algebra layer): `build_context(spde) -> (sig, base, unknowns)` parses
the SPDE into a `Signature` (the parametric vocabulary all algorithms thread) plus the
nonlinearity base; `check_subcritical(sig)` is the rule-based subcriticality test (automatic in
`build_context`).

## `RenormalizedEquation` (`counterterms.renorm.equation`)

| Member | Type / signature | Meaning |
|---|---|---|
| `.counterterms` | `list[Counterterm]` | flat list over all components |
| `.per_component` | `dict[int, list[Counterterm]]` | per equation of a system |
| `.n_components` | `int` | |
| `.all_trees` | `tuple[DecoratedTree, ...]` | the raw divergent set (before Υ-zero drops) |
| `.counterterm_rhs(comp=0)` | `sympy.Expr` | assembled \(\sum k_\tau/S(\tau)\, F(\tau^*)\) |
| `.original_rhs(comp=0)` | `sympy.Expr` | the parsed input right-hand side |
| `.summary()` | `str` | one line per counterterm |
| `.canonical(reduced=False, symmetric=True)` | `CanonicalFamily` | the canonical (BPHZ) constants as SymPy objects |
| `.bphz(symmetric=True)` | `CanonicalFamily` | sugar for `canonical(reduced=True)` |
| `.report(canonical=False, reduced=False, symmetric=True)` | `str` | full text report |
| `.render(fmt="text", canonical=…, reduced=…, symmetric=…)` | `str` | `fmt ∈ {"text","markdown","json","latex","html"}` |
| `.to_json(…)` | `str` | machine-readable report |
| `.latex_document(…)` | `str` | standalone LaTeX |
| `.save(stem="equation", outdir="output", …)` | paths | writes all formats (+ PDF if `latexmk` exists) |

**`Counterterm`** fields: `.tree`, `.homogeneity`, `.symmetry_factor`, `.elem_diff`
(\(F(\tau^*)\)), `.constant` (\(k_\tau\)); properties `.coefficient` (\(k_\tau/S(\tau)\)),
`.term` (the full summand) and `.divergence` (the leading \(\varepsilon \to 0\) rate of the
canonical constant: \(\varepsilon^{|\tau|_0}\), or \(\log(1/\varepsilon)\) at \(|\tau|_0 = 0\)).

**`CanonicalFamily`** (from `eq.canonical()` / `eq.bphz()`): `.constants` (`{τ: value of k_τ}`, a
polynomial in the elementary expectations \(h_i\)), `.free_constants` (`{τ: k_τ}`), `.legend`
(`{h_i: σ_i}`), `.substitution` (`{k_τ: value}`, for `.xreplace`), `.rhs(comp=0)` (the canonically
renormalized right-hand side) and `.expectations()` (`{h_i: Expectation}` — the Wick integral of
each *bare* \(\sigma_i\)).

!!! warning "Gaussian noises only"
    The canonical paths (`canonical`, `bphz`, `report(canonical=True)`, `canonical_character`,
    `expectation`) assume independent centered Gaussian noises. They raise `NotImplementedError`
    on a `daprato_lift`-ed equation, whose Wick-power noises are correlated and non-Gaussian.
    The free-constant family and the raw `bhz_character` remain available.

## Pictures (`counterterms.render`)

| Callable | Returns | Purpose |
|---|---|---|
| `svg(t, sig=None, *, pairing=(), scale=1.0)` | `str` | one tree as a standalone inline `<svg>`; `pairing=((i, j), …)` adds Wick arcs between noise nodes (pre-order ids) |
| `html_document(eq, canonical=False, reduced=False, symmetric=True)` | `str` | the single-file explorer (= `eq.render("html", …)`) |
| `renorm.scheme.pairings(tree, sig)` | `list[tuple[tuple[int, int], ...]]` | every global Wick pairing of a bare tree — the terms of `expectation` |

`DecoratedTree._repr_html_` returns `svg(self)`, so a tree displays as its picture in a
notebook. See [Seeing the trees](guide/seeing-trees.md).

## The playground (`counterterms.playground`)

`playground(regularity="-3/2 - kappa", dim=1, order=2, noises=None)` builds a `Playground`: a
minimal `Signature` with an empty rule, for hand-built trees (no tree *generation* — that needs
an `SPDE`). Pass `noises={"xi": …, "eta": …}` for several noises.

| Member | Purpose |
|---|---|
| `.Xi`, `.noise(name=None)` | the noise leaf \(\Xi\) |
| `.I(sub, p="", comp=0)` | an integration edge \(\mathcal{I}_p[\text{sub}]\); `p` as coordinate letters (`"x"`, `"xx"`) or a tuple |
| `.tree(*edges, X="")` | a \(\bullet\) node with those edges and polynomial decoration `X` |
| `.homogeneity(t)`, `.delta_minus(t)`, `.delta_plus(t)`, `.antipode(t)` | \(|\tau|\), \(\delta^-\), \(\Delta^+\), \(S'_-\) |
| `.character(t)`, `.expectation(t)` | the canonical constant \(h(S'_-\tau)\) and the Wick integral \(h_\varepsilon(\tau)\) |

```python
from counterterms import playground

pg = playground("-3/2 - kappa")                        # KPZ noise, d = 1
t = pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"))        # ●·𝓘ₓ[Ξ]²
pg.homogeneity(t), t.symmetry_factor(), pg.character(t), pg.expectation(t)
```

## Trees and the rule

| Callable | Module | Purpose |
|---|---|---|
| `DecoratedTree`, `tree(...)`, `red_node(...)` | `trees/tree.py` | the tree data type; canonical form, `.symmetry_factor()`, `.homogeneity(sig)` |
| `generate_counterterms(sig)` | `equation/generate.py` | the divergent set \(\mathcal{B}_{<0}\) |
| `generate_trees(sig)` | `equation/generate.py` | the full (bounded) rule-conforming pool |
| `elem_diff(t, comp, base, sig)` | `renorm/nonlinearity.py` | the elementary differential \(F_\text{comp}(\tau^*)\) (Υ-map) |
| `Homogeneity`, `Scaling` | `core/homogeneity.py` | the ordered ring \(\mathbb{Q} \oplus \mathbb{Q}\kappa\) |
| `jet(comp, k)`, `is_jet`, `jet_parts` | `core/jets.py` | jet variables \(u^c_k\) |
| `Signature` | `core/signature.py` | node/edge vocabulary, homogeneity data, rule caps |

## The algebra (`counterterms.trees.coproducts`, `counterterms.structures`)

| Callable | Purpose |
|---|---|
| `delta_plus(t, sig[, project_left])` | recentering \(\Delta : T \to T \otimes T^+\) / structure coproduct |
| `delta_minus(t, sig)` | extraction–contraction \(\delta\) |
| `delta_minus_group(t, sig)` | \(\delta^- : U^- \to U^- \otimes U^-\) |
| `twisted_antipode(t, sig)` | \(S'_- : U^- \to \mathbb{R}[U]\) |
| `RenormalizationStructure` | `.coaction`, `.coproduct`, `.twisted_antipode`, `.h_symbol`, `.bhz_character(t)`, `.canonical_character(t)`, `.divergent` |
| `RegularityStructure` | `.model_basis`, `.grades()`, `.homogeneities()`, `.divergent`, `.recentering(t)`, `.structure_coproduct(b)`, `.structure_antipode()`, `.positive_basis()` |
| `RenormalizationGroup` | `.character(values)`, `.convolve(f, g)`, `.inverse(f)`, `.identity()`; `.admissible()` raises `NotImplementedError` |
| `convolve`, `antipode`, `comodule_action` (`core/hopf.py`) | generic, basis-agnostic Hopf operations |

## The BPHZ scheme (`counterterms.renorm.scheme`)

| Callable | Purpose |
|---|---|
| `expectation(tree, sig, law=WHITE_NOISE)` | \(h(\sigma)\) as an explicit Wick-pairing integral (`Expectation`; `.is_zero`, `str(...)`) |
| `NoiseLaw`, `WHITE_NOISE` | the (symbolic) noise law; covariance display symbol |
| `zero_note(tree, sig)` | which provable identity zeroes \(h(\sigma)\), if any (root \(X^n\) / parity / total derivative) |
| `spatial_reflection_zero(tree, sig, symmetric=True)` | the reflection identity (gated) |
| `expectation_key(tree)` | duplicate detection (\(\alpha\)-independence of extended decorations) |
| `FreeConstants`, `BPHZ` | `RenormalizationScheme` implementations; `BPHZ.numeric_character` raises until an `IntegralEvaluator` is supplied |

## Rendering (`counterterms.render`)

| Callable | Purpose |
|---|---|
| `shorthand(t, sig)` | one-line tree notation |
| `ascii_art(t, sig)` | terminal drawing |
| `forest(t, sig)` | LaTeX `forest` code |
| `render(eq, fmt, …)` | the report engine behind `eq.render` |
| `structure_json(spde)`, `export_structure`, `tree_to_dict` | JSON export of the full algebraic structure |

## Reading the source

The layered architecture (one direction of dependency: `core` → `trees` → `equation` →
`renorm` → `structures` → `render`) is documented in `notes/architecture.md`; a guided reading
order is in `ENTRYPOINTS.md`. Non-obvious mathematical choices carry a tex-line citation at the
point of use.
