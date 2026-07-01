# Validation

There is no reference implementation of this pipeline anywhere — the source paper is the only
oracle. The validation strategy is therefore: **pin the paper’s worked examples as golden tests,
property-test the algebraic laws, and cross-check the classical equations against the published
literature.** `uv run pytest` runs all of it (187 tests, ~10 s).

## Tier 1 — the paper’s own worked examples (exact)

**The gKPZ five counterterms** *(tex 6004–6012)*. For
\((\partial_t - \Delta + 1)u = f(u)\zeta + g(u)(\partial_x u)^2\) at \(\beta_0 = -1-\kappa\),
the paper displays the renormalized equation explicitly. The test asserts exact multiset
equality of \((|\tau|, S(\tau), F(\tau^*))\) triples. This is a *stringent* check of the
symmetry factor and Υ-map simultaneously, because the factor 2 behaves differently in two terms:
it **survives** in \(2k\,f g\,\partial_x u\) (one derivative edge, \(S=1\)) and is **cancelled**
in \(k\,f^2 g\) (two identical edges, \(S=2\)) — both must be individually right for the display
to reproduce, and it does, term for term.

**The 43-tree table** *(tex 6028–6163)*. At \(\beta_0 = -\tfrac32 - \kappa\) the paper tabulates
every strongly-conforming gKPZ tree in six homogeneity rows, with
\((1, 2, 6, 2, 23, 9)\) trees per row. The engine reproduces the table row for row — and
building this test caught (and fixed) a genuine tree-generation undercount, which is what golden
tests are for.

**The \(\Phi^4_3\) lift** *(tex 2026–2034)*: the remainder equation
\(-v^3 - 3v^2X - 3vX^2 - X^3\) with \(X^k \in \mathcal{C}^{(-k/2)^-}\), line for line.

**Coproduct examples** *(tex 6168–6205)*: the paper’s worked \(\delta^-\) expansions, including
the term with combinatorial coefficient 2, reproduce exactly.

These goldens were **hand-derived from the paper**, not recorded from engine output — the
distinction between a validation and a tautology. (An independent audit re-derived them from the
tex and confirmed.)

## Tier 2 — the classical renormalized equations (literature, structural)

Running the engine on the standard models and applying the appropriate exact reductions
reproduces the published results:

| Model | Published result | Engine after reduction |
|---|---|---|
| KPZ | one diverging constant \(-C_\varepsilon\) (Hairer) | 8 counterterms → parity + reflection → **one constant** |
| PAM (d=2) | \(u(\xi - C_\varepsilon)\) | 4 → root-\(X^n\) + parity → \(u \cdot \text{const}\) |
| gPAM (d=2) | \(-C_\varepsilon f(u)f'(u)\) (BCCH) | same, exactly \(f f'\) |
| \(\Phi^4_2\) | mass renormalization \(+3C\varphi\) | mass term \(\propto v\) survives |
| \(\Phi^4_3\) | two diverging mass constants | two mass constants \(\propto v\); gradient terms vanish |

## Property tests — the algebraic laws

The Phase-3 layer is tested by verifying the *laws*, parametrized over a corpus of six equations
including the singular \(\beta_0 = -\tfrac32\):

- counits, coassociativity of \(\Delta^+\) and \(\delta^-\), the comodule conditions;
- the **cointeraction** between extraction and recentering *(tex 5717)* — uncapped, over all
  divergent KPZ trees;
- the twisted antipode’s characterizing (Dyson–Salam) relation *(tex 5034)*, checked numerically
  in exact rationals;
- the group axioms of \(G^-\) (associativity, unit, inverse), including on multi-component
  forests;
- homogeneity stability of every coproduct.

## Independent cross-checks

- \(S(\tau)\) is verified against a **brute-force automorphism count** (permutations, not the
  product formula) — so the golden \(S\) values don’t test the formula against itself.
- The Υ-map’s non-commuting order of operations (\(\partial_p\) before \(D^n\)) is tested on a
  tree where the wrong order gives a *different* answer.
- The reflection identity in the reduced view emerges, independently, with the paper’s own
  admissibility statement \(h(X^n\tau) = 0\) *(tex 5083)* as a by-product — an unplanned
  consistency success.
- The repository has been through code-quality and mathematical audits, the latter with
  independent adversarial re-derivation of every load-bearing formula (see `notes/math_audit.md`
  and `notes/validation.md` in the repository).

## What validation does *not* cover

Anything behind the analytic wall (numeric constants, models, convergence) — there is nothing to
validate because nothing is computed. The remaining sharp edges listed in
[Scope & limitations](scope.md#known-sharp-edges) sit outside the test corpus's reach; the
generation defects an audit found there (order > 2 enumeration, partial-sum pruning) have been
fixed and are now pinned by regression tests in `tests/test_rule.py`.
