# The algebraic layer

Everything in this page is optional for the headline use case — the renormalized family needs no
Hopf algebra. But the full BHZ structure is built, property-tested, and exposed, for when you
want the objects themselves: coproducts, the twisted antipode, characters, the group.

## The regularity structure \((T, T^+)\)

```python
from counterterms import build_regularity_structure

rs = build_regularity_structure(spde, gamma=1)   # γ-truncated
```

- `rs.model_basis` — the decorated trees spanning the model space \(T\) up to homogeneity
  \(\gamma\); `rs.grades()` groups them by homogeneity, `rs.homogeneities()` lists the grading.
- `rs.divergent` — the negative subspace (the counterterm carriers).
- `rs.positive_basis()` — the trees generating \(T^+\).
- `rs.recentering(t)` — the coproduct \(\Delta : T \to T \otimes T^+\) (positive
  renormalization: Taylor recentring of the model). Returned as an exact tensor sum
  `{(left, right): coefficient}`.
- `rs.structure_coproduct(b)` — \(\Delta^+ : T^+ \to T^+ \otimes T^+\), and
  `rs.structure_antipode()` — the antipode of the structure group’s Hopf algebra.

!!! note
    The polynomial sector \(\{X^k\}\) is carried implicitly as node decorations, not as
    standalone basis vectors — consumers counting graded dimensions should account for this.

## The renormalization structure and the BHZ character

```python
from counterterms import build_renormalization

rn = build_renormalization(spde)
```

- `rn.coaction(t)` — \(\delta : U \to U^- \otimes U\), the extraction–contraction coaction
  (cut out divergent subforests; the contracted node keeps an extended decoration; Taylor
  \(\mathfrak{e}'\) terms recentre the cut edges).
- `rn.coproduct(t)` — \(\delta^- : U^- \to U^- \otimes U^-\).
- `rn.twisted_antipode(t)` — the negative twisted antipode \(S'_-\), the recursive
  Bogoliubov-style subtraction *(tex 5034)*. Note \(S'_-\) is *not* a Hopf antipode — it maps
  into the free algebra \(\mathbb{R}[U]\).
- `rn.h_symbol(σ)` — the elementary-expectation symbol \(h(\sigma)\).
- `rn.bhz_character(t)` — the exact polynomial \(k_\tau = h(S'_- t)\) in the \(h\)-symbols.
- `rn.canonical_character(t)` — the same with the centered-Gaussian parity rule applied
  (odd-noise trees ⇒ 0). Symbols that survive here are *not* guaranteed nonzero — the further
  provable reductions live in the [reduced report view](reports.md#canonical-and-reduced-views).

The **cointeraction** between \(\delta\) and the recentering \(\Delta\) *(tex 5717)* — the
compatibility that makes renormalization commute with recentring — is verified by the test suite
over the whole equation corpus, **including the singular \(\beta_0 = -\tfrac32\)** case, where a
naive reading of the paper’s recentring actually fails (see
`notes/cointeraction_singular_noise.md` for that story).

## The renormalization group \(G^-\)

```python
from counterterms import build_renormalization_group

G = build_renormalization_group(spde)
chi = G.character({tree: value, ...})   # a character from tree values (multiplicative on forests)
G.convolve(chi1, chi2)                  # the group product (χ₁ ⊗ χ₂) ∘ δ⁻
G.inverse(chi)                          # χ ∘ S  (antipode inverse)
G.identity()                            # the counit
```

Group axioms (associativity, unit, inverse) are property-tested, including on multi-component
forests. `G.admissible()` — the analytic subgroup \(G^-_{\mathrm{ad}}\) of characters compatible
with the kernel’s vanishing moments — is deliberately unbuilt and raises `NotImplementedError`
(it needs the \(\Pi\)-map, an analytic object).

## Wick structure of the canonical constants

```python
from counterterms.renorm.scheme import expectation, NoiseLaw, WHITE_NOISE

exp = expectation(tree, sig, law=WHITE_NOISE)
str(exp)     # the ε-regularized integral: heat kernels ∂^p K against covariances C(z_i − z_j)
```

`expectation` writes \(h(\sigma) = \mathbb{E}[\Pi^\zeta\sigma](0)\) as an explicit Isserlis /
Wick-pairing sum: one \(\partial^p K\) per kernel edge, one covariance factor per matched noise
pair, within-noise-type pairing only, odd counts ⇒ 0. The integrals are emitted, **not**
evaluated — evaluation is the deliberately unbuilt analytic half
(`BPHZ.numeric_character` raises until an `IntegralEvaluator` exists).

## JSON export

```python
from counterterms.render import structure_json

doc = structure_json(spde)   # one self-contained document
```

Exports the whole structure machine-readably: the basis trees (as nested dicts, exact
homogeneities), the divergent set, the coproducts \(\Delta, \delta^-\) as tensor sums, antipode
values, and the character polynomials — enough to reconstruct or cross-check the algebra in
another system.
