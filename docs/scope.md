# Scope & limitations

The project’s policy is *explicit errors over silent wrong output*, and honest documentation of
what is unbuilt or imperfect. This page is the complete list, including known defects.

## In scope

- Scalar equations **and coupled systems**; one **or several** independent noises.
- Second-order parabolic \(L\) (the proven case); other orders run with a warning (see below).
- Subcritical noise, \(\beta_0 > -\text{order}\), checked rule-based on every input.
- Nonlinearities affine in the noise; \(g\) at most quadratic in \(\partial u\) (Assumption D2);
  derivative factors with \(|p|_\mathfrak{s} \le 1\).
- Supercritical **polynomial additive-noise** equations (\(\Phi^4_2\), \(\Phi^4_3\)) via the
  [da Prato–Debussche lift](guide/daprato.md).

## Rejected, by design

- Supercritical equations not liftable here (sine-Gordon — needs Wick exponentials).
- Nonlinearities not affine in the noise (\(\xi^2\), \(f(\xi)\)).
- \(g\) more than quadratic in \(\partial u\); factors \(\partial_t u\) or \(\partial^2_x u\)
  in the nonlinearity (\(|p|_\mathfrak{s} > 1\)).
- Quasilinear or non-parabolic operators (not expressible in the DSL).

## The analytic wall (deliberately unbuilt)

The engine is symbolic. It stops, on purpose, exactly where analysis/probability starts:

- **No numeric constant values.** The canonical constants are exact polynomials in
  elementary-expectation symbols \(h(\sigma)\), each spelled out as an
  \(\varepsilon\)-regularized Wick integral — never evaluated. `BPHZ.numeric_character` raises
  until an `IntegralEvaluator` is plugged in.
- **No models, estimates, convergence, or solving.**
- **\(G^-_{\mathrm{ad}}\)** (the admissible subgroup — needs kernel moments and the \(\Pi\)-map)
  raises `NotImplementedError`.
- **\(\Phi^4\) canonical values are structure-only**: the lift treats Wick powers as independent
  noises, so the free family is exact but canonical constants ignore the true correlations among
  \(:\!X^k\!:\).

## Known sharp edges

None affects the golden-tested regime (order 2, the validated equations).

### Operator order > 2

Schauder/admissibility is proven only for 2nd-order parabolic \(L\); the engine warns and
computes the combinatorics anyway. The tree *enumeration* is complete for any order (the
node-decoration cap scales with \(-\beta_0\), and the generator does not prune on intermediate
homogeneity sums — both were audit findings, fixed with regression tests in
`tests/test_rule.py`), but treat order ≠ 2 output as exploratory: the underlying analytic
theory is unverified there.

### Asymmetric systems and the raw tree basis

Gradient-degree budgets are shared across components (a `max` over equations), so a system in
which only one equation has a gradient nonlinearity can over-generate trees in the **raw** basis
(`all_trees`, and hence spurious always-zero entries in the canonical constant *list*). The
renormalized equations themselves are unaffected — the spurious trees have \(F(\tau^*) = 0\) and
are dropped.

### Smaller notes

- **Multi-noise covariances share one display symbol**: with several noises of *different* laws,
  all covariances render as the single symbol \(C\) in the Wick integrands. The pairing
  combinatorics (within-type only) are correct; only the label conflates. The API currently
  cannot express distinct laws anyway.
- **Subcriticality compares standard parts only**: a hypothetical noise at
  \(\beta_0 = -\text{order} + \kappa\) (infinitesimally subcritical) is rejected although the
  ordered-ring comparison would accept it. One-sided — it never wrongly *accepts*; no physical
  noise sits there.
- **`symmetric=True` is the default** on the reduced report view. It is an assumption about your
  noise (reflection invariance), stated in-band in every report; pass `symmetric=False` for
  anisotropic noise. In \(d \ge 2\) the reflection reduction uses total spatial parity only —
  conservative: it never wrongly zeroes, but misses per-axis reductions.
- **Report text ordering is hash-seed dependent** (the RULE lines); JSON is the stable format.

## Why trust any of it

Because the scope above is enforced by tests, the golden cases are re-derived from the paper
rather than recorded from the engine, and the codebase has been audited — including
adversarially, by independent re-derivation of the core formulas. See
[Validation](validation.md).
