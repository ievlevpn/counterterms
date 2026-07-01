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
| `Noise` | `Noise(name, regularity)` | `regularity` **is** \(\beta_0\) (may carry `- kappa`); exposes `.symbol`, `.homogeneity` |
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
| `.report(canonical=False, reduced=False, symmetric=True)` | `str` | full text report |
| `.render(fmt="text", canonical=…, reduced=…, symmetric=…)` | `str` | `fmt ∈ {"text","markdown","json","latex"}` |
| `.to_json(…)` | `str` | machine-readable report |
| `.latex_document(…)` | `str` | standalone LaTeX |
| `.save(stem="equation", outdir="output", …)` | paths | writes all formats (+ PDF if `latexmk` exists) |

**`Counterterm`** fields: `.tree`, `.homogeneity`, `.symmetry_factor`, `.elem_diff`
(\(F(\tau^*)\)), `.constant` (\(k_\tau\)); properties `.coefficient` (\(k_\tau/S(\tau)\)) and
`.term` (the full summand).

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
