# counterterms

A **symbolic engine for the renormalization of singular SPDEs**. You describe a subcritical
singular stochastic PDE; it returns the **family of renormalized equations** — the original PDE
plus the tree-indexed counterterms the theory prescribes — following Bailleul & Hoshino,
*"A tourist's guide to regularity structures and singular stochastic PDEs"*
([arXiv:2006.03524](https://arxiv.org/abs/2006.03524)), which packages the BHZ/BCCH algebra.

> ⚠️ **Personal research library — no guarantees.** A one-person research project, not production
> software. It is validated against the paper where possible (see [Validation](#validation)), but it
> may be wrong, incomplete, or break on unfamiliar inputs. No warranty, no stability promises, no
> support. Check the output against the mathematics before trusting it.

## It is symbolic

The engine computes *what* to renormalize and the *structure* of the renormalization — which
counterterms appear, with which tree combinatorics, symmetry factors, and elementary differentials.
It works in exact symbolic arithmetic (SymPy + an ordered ring $\mathbb{Q}\oplus\mathbb{Q}\kappa$
for homogeneities).

It does **not** compute the *numeric values* of the renormalization constants (those require
Gaussian/Wick integrals, deferred), and it is **not** a solver: no model construction, no analytic
estimates, no convergence, no time-stepping. The free constants are the complete symbolic answer —
a solution of a singular SPDE *is* the family indexed by the renormalization group.

## What it handles

**Equations** (subcritical, semilinear, parabolic):

- **KPZ** and **generalized KPZ**, $f(u)\zeta + g(u)(\partial_x u)^2$
- **PAM** and **generalized PAM**, $f(u)\zeta$
- **coupled systems** of equations, and **several independent noises**
- **$\Phi^4_2$, $\Phi^4_3$** — supercritical, handled via the da Prato–Debussche lift

**Generality:** scalar or systems; one or many noises; 2nd-order parabolic operator (general
operator order is allowed with a warning); the noise regularity $\beta_0$ is user-supplied and
subcriticality ($\beta_0 > -\text{order}$) is checked, so the counterterm set is finite.

**Assumptions** (enforced — out-of-scope inputs are rejected with a clear error): the nonlinearity
is affine in the noise; $g$ is at most quadratic in $\partial u$ (Assumption D2); derivative factors
satisfy $|p|_\mathfrak{s} \le 1$.

## What you get

- the **renormalized family** — the original PDE plus one counterterm $\frac{k_\tau}{S(\tau)}F(\tau^*)$
  per negative-homogeneity tree, with free constants $k_\tau$;
- the **canonical (BPHZ) character** $k_\tau = h(S'_-\tau)$ — symbolic, as polynomials in the
  elementary expectations $h(\sigma)$, with the centered-Gaussian parity rule applied;
- a **reduced** view that folds in the exact identities (vanishing & duplicate constants, and — for a
  spatially-symmetric noise — the reflection identity);
- typeset **reports** in text / Markdown / JSON / LaTeX→PDF, with the trees drawn;
- the underlying **algebra**: the regularity structure $(T, T^+)$, the extraction/recentering
  coproducts, the twisted antipode, the BHZ character, the renormalization group $G^-$, and a JSON
  export of the whole structure.

## Install

[`uv`](https://docs.astral.sh/uv/); SymPy is the only runtime dependency.

```sh
uv sync
uv run pytest          # 185 tests, ~8s
```

## Quickstart

```python
from sympy import Derivative, Function, Rational
from counterterms import SPDE, Noise, Parabolic, Unknown, kappa

u  = Unknown("u", dim=1)
xi = Noise("xi", regularity=Rational(-1) - kappa)              # noise regularity −1 − κ
f, g = Function("f"), Function("g")

# generalized KPZ:  (∂t − Δ + 1) u = f(u) ξ + g(u) (∂x u)²
spde = SPDE(operator=Parabolic(dim=1, mass=1), unknown=u, noises=[xi],
            rhs=f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0])**2)

print(spde.renormalize().summary())     # the five gKPZ counterterms, each with |τ|, S(τ), F(τ*)
```

The renormalized equation it produces is

$$(\partial_t - \Delta + 1)\,u = f(u)\zeta + g(u,\partial u) + \sum_{\tau \in \mathcal{B},\,|\tau|<0} \frac{k_\tau}{S(\tau)}\,F(\tau^*)(u,\partial u).$$

The same `spde` also yields the **canonical (BPHZ) character** — the symbolic constant
$k_\tau = h(S'_-\tau)$ the twisted antipode prescribes (odd-noise trees vanish by parity; the
$h(\sigma)$ stay free symbols, since their Wick-integral values are out of scope):

```python
from counterterms import build_renormalization

rs = build_renormalization(spde)
for t in rs.divergent:
    print(rs.canonical_character(t))    # exact antipode combination in the h-values
```

More runnable scripts in [`examples/`](examples/) — start with
`uv run python -u examples/01_renormalized_equation.py`.

## Example output

`eq.save()` writes the report as text / Markdown / JSON and a typeset LaTeX → PDF. Both reports below
are for **KPZ**, $(\partial_t - \Delta) u = (\partial_x u)^2 + \xi$ at $\beta_0 = -3/2 - \kappa$:
the parsed equation, every divergent tree $\tau$ (drawn in the paper's convention — $\circ$ noise,
$\bullet$ integration node, dotted = derivative kernel) with its homogeneity $|\tau|$, symmetry
factor $S(\tau)$, free constant $k_\tau$ and elementary differential $F(\tau^*)$, the assembled
renormalized family, and the canonical (BPHZ) constants $k_\tau = h(S'_-\tau)$.

- **[`docs/kpz_canonical.pdf`](docs/kpz_canonical.pdf)** (`canonical=True`) — the BPHZ constants for a
  *general* centered-Gaussian noise. Each $k_\tau$ is a polynomial in the elementary expectations
  $h(\sigma)$; entries that provably vanish ($=0$) or duplicate another ($=h_j$) are *marked* but
  left in place. Correct for any noise (a non-symmetric mollifier genuinely keeps a $\partial_x u$
  drift counterterm).
- **[`docs/kpz_reduced.pdf`](docs/kpz_reduced.pdf)** (`reduced=True`) — the same constants *fully
  reduced* for a spatially-symmetric (e.g. white) noise: provable zeros and o-duplicates folded in,
  **plus** the spatial-reflection identity $x\to-x$, collapsing KPZ to a **single diverging
  constant**, $\partial_t u = \Delta u + (\partial_x u)^2 - C + \xi$ (Hairer's result). The reduction
  states its symmetric-noise assumption and is never claimed for an anisotropic noise.

---

## Under the hood

Everything public is re-exported from the top-level `counterterms` package. The pipeline is
`SPDE (DSL) → structural rule → decorated trees with |τ|<0 → S(τ) + F(τ*) → assembled family`,
optionally composed with the renormalization coproducts and characters.

### Input & the renormalized equation
| You want… | Call | Module |
|---|---|---|
| Write the SPDE (DSL) | `Unknown`, `Noise`, `Parabolic`, `SPDE` | `equation/dsl.py` |
| Derive the renormalized family | `SPDE(...).renormalize()` → `RenormalizedEquation` | `api.py`, `renorm/equation.py` |
| Access each counterterm (tree, $\lvert\tau\rvert$, $S(\tau)$, $F(\tau^*)$, $k_\tau$) | `eq.counterterms`, `eq.per_component` | `renorm/equation.py` |
| Render a report (text / markdown / json / latex; `canonical`, `reduced`) | `eq.summary()`, `eq.render(fmt, …)`, `eq.save(…)` | `render/report.py`, `render/latex.py` |
| $\Phi^4_2 / \Phi^4_3$ (supercritical, $\beta_0 \le -\text{order}$) | `daprato_lift(spde).renormalize()` | `equation/daprato.py` |

### The trees and the rule
| You want… | Call | Module |
|---|---|---|
| The divergent trees $\mathcal{B}_{<0}$ | `generate_counterterms(sig)` | `equation/generate.py` |
| The structural rule from the nonlinearity | `build_context(spde)` → `(sig, base, unknowns)` | `equation/dsl.py` |
| Subcriticality check ($\beta_0 > -\text{order}$) | `check_subcritical(sig)` (auto in `build_context`) | `equation/rule.py` |
| Decorated trees: canonical form, $S(\tau)$, homogeneity | `DecoratedTree`, `tree`, `red_node` | `trees/tree.py` |
| The elementary differential $F(\tau^*)$ ($\Upsilon$-map) | `elem_diff(t, comp, base, sig)` | `renorm/nonlinearity.py` |
| Draw a single tree (shorthand / ascii / `forest`) | `shorthand`, `ascii_art`, `forest` | `render/tree.py` |

### The algebraic structure (regularity & renormalization)
| You want… | Call | Module |
|---|---|---|
| The regularity structure $(T, T^+)$, graded basis | `build_regularity_structure(spde)` | `structures.py` |
| Recentering $\Delta : T \to T \otimes T^+$, structure coproduct $\Delta^+$ | `delta_plus(t, sig[, project_left])` | `trees/coproducts.py` |
| Extraction–contraction $\delta$ / $\delta^-$ (and the cointeraction) | `delta_minus`, `delta_minus_group` | `trees/coproducts.py` |
| Negative twisted antipode $S'_-$ | `twisted_antipode(t, sig)` | `trees/coproducts.py` |
| Renormalization structure + symbolic **BHZ character** $k = h \circ S'_-$ | `build_renormalization(spde)` → `.bhz_character`, `.canonical_character` | `structures.py` |
| The renormalization **group** $G^-$ (convolution, antipode inverse) | `build_renormalization_group(spde)` | `structures.py` |
| Generic, basis-agnostic Hopf ops (convolve / antipode / comodule) | `convolve`, `antipode`, `comodule_action` | `core/hopf.py` |
| Wick-pairing structure of the canonical constants (symbolic) | `expectation`, `NoiseLaw`, `BPHZ`, `FreeConstants` | `renorm/scheme.py` |
| JSON export of the whole structure | `structure_json(spde)`, `export_structure`, `tree_to_dict` | `render/export.py` |

### Foundations
| You want… | Call | Module |
|---|---|---|
| The ordered homogeneity ring $\mathbb{Q} \oplus \mathbb{Q}\cdot\kappa$ | `Homogeneity`, `kappa`, `Scaling` | `core/homogeneity.py` |
| Jet variables $u^c_k$ | `jet`, `is_jet`, `jet_parts` | `core/jets.py` |
| The `Signature` (the parametric vocabulary everything threads) | `Signature` | `core/signature.py` |
| The `Symbol` protocol (basis seam) | `Symbol` | `core/symbol.py` |

### Out of scope (rejected with clear errors)
Non-polynomial supercritical equations (e.g. sine-Gordon, which needs Wick exponentials); noise
nonlinearities not affine in the noise; $g$ more than quadratic in $\partial u$; singular derivative
factors $|p|_\mathfrak{s} > 1$; quasilinear or non-parabolic operators. And — the analysis/probability
wall — no numeric constant values (the Wick integrals are emitted symbolically by `renorm/scheme.py`
but not evaluated), no model construction, estimates, convergence, or solving.

## Status

Phases 1–3 complete and green; Phase 4 partially built.

- **Phase 1–2** — the renormalized family from an `SPDE`, for scalar/systems, one/many noises, general operator order.
- **Phase 3** — the coproducts (the cointeraction holds **including $\beta_0 = -3/2$**), `RegularityStructure` $(T, T^+)$, the generic `core/hopf` layer, subcriticality, the twisted antipode + BHZ character, the group $G^-$.
- **Phase 4 (partial)** — `daprato_lift` ($\Phi^4_{2/3}$), the canonical-character *symbolic* half + Wick parity, the reduced view, the JSON export, and the seams for the unbuilt analytic pieces (evaluating the divergent integrals).

See [`ROADMAP.md`](ROADMAP.md) and [`CHANGELOG.md`](CHANGELOG.md).

## Validation

The paper is the only oracle (no reference implementation exists). The suite reproduces the gKPZ
example's exact five counterterms ($\beta_0 = -1$) **and** the full strongly-conforming-tree table at
$\beta_0 = -3/2$ (43 trees, tex 6024–6163) — matching the latter caught and fixed a real
tree-generation bug. The renormalized KPZ / PAM / gPAM / $\Phi^4$ equations are cross-checked against
the published results (Hairer; BCCH). Plus an independent symmetry-factor cross-check and benchmark
counts. See [`notes/validation.md`](notes/validation.md).

## Documentation

- [`examples/`](examples/) — runnable quickstart scripts (start here).
- [`ENTRYPOINTS.md`](ENTRYPOINTS.md) — a guided reading order through the source.
- [`notes/use_cases.md`](notes/use_cases.md) — what it can and cannot solve, honestly.
- [`notes/initial_plan.md`](notes/initial_plan.md) — the mathematics (pipeline, conventions, scope).
- [`notes/architecture.md`](notes/architecture.md) — the module structure and design.
- [`notes/cointeraction_singular_noise.md`](notes/cointeraction_singular_noise.md),
  [`notes/phase4_plan.md`](notes/phase4_plan.md), [`notes/validation.md`](notes/validation.md) — deep dives.
- [`references/`](references/) — pointer to the source paper ([arXiv:2006.03524](https://arxiv.org/abs/2006.03524)), cited by line number throughout.

## License

[MIT](LICENSE) © 2026 Pavel Ievlev. The Bailleul–Hoshino paper is *not* included (it is the
authors' copyright) — see [`references/`](references/) for the arXiv link.
