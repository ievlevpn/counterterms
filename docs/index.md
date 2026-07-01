# counterterms

A **symbolic engine for the renormalization of singular stochastic PDEs**. You describe a
subcritical singular SPDE; it returns the **family of renormalized equations** — the original
PDE plus the tree-indexed counterterms the theory prescribes:

$$(\partial_t - \Delta + 1)\,u \;=\; f(u)\,\zeta + g(u,\partial u)
\;+\; \sum_{\tau \in \mathcal{B},\, |\tau|<0} \frac{k_\tau}{S(\tau)}\, F(\tau^*)(u,\partial u).$$

It implements the pipeline of Bailleul & Hoshino, *“A tourist’s guide to regularity structures
and singular stochastic PDEs”* ([arXiv:2006.03524](https://arxiv.org/abs/2006.03524)), which
packages the BHZ/BCCH renormalization algebra. Every formula in the codebase is traceable to a
line of that paper.

!!! warning "Personal research library — no guarantees"
    A one-person research project, not production software. It is validated against the paper
    where possible (see [Validation](validation.md)), but it may be wrong, incomplete, or break
    on unfamiliar inputs. No warranty, no stability promises. Check the output against the
    mathematics before trusting it. Known gaps are listed honestly in
    [Scope & limitations](scope.md).

## What “symbolic” means here

The engine computes *what* to renormalize and the *structure* of the renormalization — which
counterterms appear, with which tree combinatorics, symmetry factors \(S(\tau)\), and elementary
differentials \(F(\tau^*)\). It works in exact symbolic arithmetic: SymPy for the functional
expressions, and an ordered ring \(\mathbb{Q} \oplus \mathbb{Q}\kappa\) for homogeneities (no
floats anywhere near a critical tree).

It does **not**:

- compute the *numeric values* of the renormalization constants (those are divergent Gaussian /
  Wick integrals; the engine emits them symbolically but does not evaluate them),
- construct models, prove estimates, or solve the equation.

That is not a hedge — it is the point. In the theory, a “solution” of a singular SPDE *is* the
family of renormalized equations indexed by the renormalization group; the free constants
\(k_\tau\) are the complete symbolic answer.

## What it handles

| Equation | Input | Result |
|---|---|---|
| KPZ, generalized KPZ | \(f(u)\zeta + g(u)(\partial_x u)^2\), \(\beta_0 \in (-2, 0)\) | 5 counterterms at \(\beta_0=-1\); the paper’s full 43-tree table at \(\beta_0=-\tfrac32\) |
| PAM, gPAM | \(f(u)\zeta\), \(d = 2\) | the \(C\,f(u)f'(u)\) counterterm |
| Coupled systems | several equations sharing noises | shared constants across components |
| Multiple noises | independent \(\xi, \eta, \dots\) | per-type Wick pairing |
| \(\Phi^4_2\), \(\Phi^4_3\) | supercritical additive noise | via the [da Prato–Debussche lift](guide/daprato.md) |

**Enforced assumptions** (out-of-scope input is rejected with a clear error): the nonlinearity is
affine in the noise, \(g\) is at most quadratic in \(\partial u\) (Assumption D2), derivative
factors satisfy \(|p|_\mathfrak{s} \le 1\), the operator is second-order parabolic (other orders
run with a warning — [read this first](scope.md#operator-order-2)).

## What you get

- the **renormalized family** — one counterterm \(\frac{k_\tau}{S(\tau)} F(\tau^*)\) per
  negative-homogeneity decorated tree, with free symbolic constants \(k_\tau\);
- the **canonical (BPHZ) character** \(k_\tau = h(S'_-\tau)\), symbolic in the elementary
  expectations \(h(\sigma)\), with the centered-Gaussian parity rule applied;
- a **reduced** view that folds in the exact identities (vanishing and duplicate constants,
  plus — for a spatially symmetric noise — the reflection identity), collapsing e.g. KPZ to
  Hairer’s single constant;
- typeset **reports** in text / Markdown / JSON / LaTeX → PDF, with the trees drawn;
- the underlying **algebra**: the regularity structure \((T, T^+)\), the extraction and
  recentering coproducts, the twisted antipode, the BHZ character, the renormalization group
  \(G^-\), and a JSON export of the whole structure.

Example output: [KPZ, canonical view (PDF)](kpz_canonical.pdf) ·
[KPZ, reduced view (PDF)](kpz_reduced.pdf).

## Where to go next

- [Getting started](getting-started.md) — install and run your first renormalization in five lines.
- [The mathematics](mathematics.md) — what the pipeline computes and why, explained from scratch.
- [User guide](guide/writing-spdes.md) — the DSL, reports, the lift, and the Hopf-algebra layer.
- [API reference](api.md) — every public entry point.
