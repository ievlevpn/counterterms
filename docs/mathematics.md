# The mathematics

This page explains, from scratch, what the engine computes and why. It follows Bailleul &
Hoshino’s *tourist guide* ([arXiv:2006.03524](https://arxiv.org/abs/2006.03524)); parenthetical
references like *(tex 4337)* are line numbers in the paper’s LaTeX source, which is how the
codebase cites its oracle throughout.

## 1. Why singular SPDEs need renormalization

Consider a semilinear parabolic equation driven by a very rough noise \(\zeta\):

$$(\partial_t - \Delta + 1)\, u = f(u)\,\zeta + g(u, \partial u).$$

If \(\zeta\) has Hölder regularity \(\beta_0 < 0\) (e.g. spacetime white noise in one space
dimension has \(\beta_0 = -\tfrac32 - \kappa\)), the solution \(u\) is only as regular as the
equation allows, and products like \(f(u)\zeta\) or \((\partial_x u)^2\) are **ill-defined**:
the sum of the regularities of the factors is negative, so classical (Young / paraproduct)
multiplication fails.

The renormalization strategy: mollify the noise (\(\zeta_\varepsilon = \zeta * \rho_\varepsilon\)),
solve the now-classical equation, and try to take \(\varepsilon \to 0\). The limit does not exist —
but it does exist if you first subtract **diverging counterterms** from the equation. The theory
of regularity structures (Hairer; Bruned–Hairer–Zambotti; Bruned–Chandra–Chevyrev–Hairer)
identifies *exactly which* counterterms, and this engine computes them symbolically.

The answer has a remarkable shape: one counterterm per **decorated tree** of negative
homogeneity,

$$(\partial_t - \Delta + 1)\, u^{(k)} = f(u^{(k)})\zeta + g(u^{(k)}, \partial u^{(k)})
 + \sum_{\tau \in \mathcal{B},\, |\tau| < 0} \frac{k(\tau)}{S(\tau)}\, F(\tau^*)(u^{(k)}, \partial u^{(k)}),$$

where the constants \(k(\tau)\) are free — each choice of them is one member of the **family of
renormalized equations** *(the paper’s Theorem “ThmRenormPDEs”, the BCCH formula)*. This display
is the engine’s output.

## 2. Subcriticality

The whole machine only works when the noise is not *too* rough relative to the smoothing of the
operator: **subcriticality** requires \(\beta_0 > -\text{order}\) (for the heat operator,
\(\beta_0 > -2\)) *(tex 5485)*. Subcriticality guarantees that iterating the equation gains
regularity, so only **finitely many** trees have negative homogeneity and the counterterm sum is
finite. The engine checks this rule-based condition on every input and rejects supercritical
equations (for polynomial additive-noise cases like \(\Phi^4_3\), it offers the
[da Prato–Debussche lift](guide/daprato.md) instead).

## 3. Decorated trees: the vocabulary of the expansion

Iterating the mollified equation (Picard iteration) produces a sum of explicit stochastic
integrals; each integral is encoded by a **decorated tree** *(tex 3963–3970)*:

- **\(\circ\) (noise node)** — an occurrence of the noise \(\zeta\). Multiple independent noises
  get distinct node types \(\circ_j\).
- **\(\bullet\) (integration node)** — a point where the nonlinearity \(g\) is evaluated.
- **Edges \(\mathcal{I}_p\)** — a convolution with the kernel \(\partial^p K\) of
  \(L^{-1}\); the multi-index \(p\) records derivatives hitting the kernel. For systems, the edge
  also carries **which equation’s** kernel it is — component identity lives on the *edge type*,
  not the node.
- **Node decorations \(X^n\)** — polynomial factors \(x^n\) from Taylor expansion (recentring).

For example, in KPZ the tree \(\bullet\!\!\Rightarrow\!\!\circ\,\circ\) (an integration node with
two derivative-kernel edges to noises) encodes \(\big(\partial_x K * \zeta\big)^2\) — the famous
ill-defined square.

Two trees are the same iff they are isomorphic as decorated combinatorial trees; the engine uses
a canonical form for equality, dictionary keys, and symmetry counting.

## 4. Homogeneity: which trees diverge

Each tree gets a **homogeneity** \(|\tau|\) — the parabolic Hölder regularity of the stochastic
object it encodes — computed by three rules:

$$|\Xi| = \beta_0, \qquad
|\mathcal{I}_p \tau| = |\tau| + \text{order} - |p|_\mathfrak{s}, \qquad
|X^n| = |n|_\mathfrak{s},$$

multiplicativity over subtrees, with the **parabolic scaling** \(\mathfrak{s} = (2, 1, \dots, 1)\)
(time counts double) *(tex 5276–5283, 2213, 1153)*. Kernel edges *gain* the Schauder order (2 for
the heat operator) and *lose* \(|p|_\mathfrak{s}\) per derivative.

Homogeneities are exact elements of the ordered ring \(\mathbb{Q} \oplus \mathbb{Q}\kappa\), where
\(\kappa > 0\) is an infinitesimal: white noise sits *just below* the rational threshold, and
critical trees sit at homogeneity \(-k\kappa\) — negative, but only by an infinitesimal. Using
floats here would silently misclassify them; the engine never does.

**A tree needs a counterterm exactly when \(|\tau| < 0\)** — then the associated stochastic
object diverges as \(\varepsilon \to 0\) and must be recentred by a constant. The set
\(\mathcal{B}_{<0}\) of such trees is what the engine enumerates. It includes the bare noises
\(\circ^n\) themselves (in KPZ these produce the \(k(\circ)\) and \(k(\circ_1)\) terms); the bare
\(\bullet\) has homogeneity \(0\) and is excluded.

## 5. The rule: which trees can occur at all

Not every decorated tree appears in the expansion of a given equation — the nonlinearity dictates
the possible branching. From the right-hand side the engine derives a **rule** *(tex 5306–5340)*:
e.g. for gKPZ, an integration node may carry at most one noise edge (the nonlinearity is affine in
\(\zeta\)) and at most two derivative-kernel edges (\(g\) is quadratic in \(\partial u\) —
Assumption D2, with the bound on the *total* gradient degree per node). Trees conforming to the
rule (“strongly conforming”, *tex 4555*) with negative homogeneity are exactly the counterterm
carriers.

Validation anchor: for gKPZ at \(\beta_0 = -\tfrac32 - \kappa\) the paper tabulates all 43
strongly-conforming trees in six homogeneity rows *(tex 6028–6163)*; the engine reproduces the
table row for row.

## 6. The coefficient: \(S(\tau)\) and the elementary differential \(F(\tau^*)\)

Each divergent tree contributes the counterterm \(\dfrac{k(\tau)}{S(\tau)} F(\tau^*)(u, \partial u)\):

**The symmetry factor** \(S(\tau)\) counts the tree’s internal symmetries *(tex 3982)*:

$$S(\tau) = n! \; \prod_j S(\sigma_j)^{m_j}\, m_j!$$

over the distinct child branches \(\sigma_j\) with multiplicities \(m_j\) (and the root
decoration’s multi-index factorial \(n!\)). Dividing by \(S(\tau)\), never by anything else, is a
load-bearing convention *(tex 4915)*.

**The elementary differential** \(F(\tau^*)\) — the “Υ-map” — turns the tree back into a function
of the solution *(tex 4337)*:

$$F(\tau^*) = \Big(\prod_i F(\tau_i^*)\Big)\cdot
\Big(D^n \prod_i \partial_{p_i}\Big) F(b^*),$$

recursively over the branches: each child branch contributes its own factor, and the base symbol
of the node (\(F(\circ_j^*) = f_j\), \(F(\bullet^*) = g\)) is differentiated once in the argument
slot matching each child’s edge decoration \(p_i\), then \(n\) more times by the total derivative
\(D_i = \sum_k u_{k+e_i} \partial_{u_k}\) *(tex 4316)*. Two subtleties the engine gets right and
tests:

- the slot derivatives \(\partial_{p_i}\) are applied **before** \(D^n\) — they do not commute;
- \(D^n\) is a genuine iterated total derivative (Faà di Bruno applies).

A worked example of why both \(S\) and \(\Upsilon\) must be individually right: in gKPZ, the tree
with **one** derivative edge to a noise gives
\(\Upsilon = \partial_{u_1}[g\,u_1^2]\cdot f = 2fg\,\partial_x u\) with \(S = 1\) — the paper’s
explicit \(2k(\tau_3)\); the tree with **two** identical derivative edges gives
\(\Upsilon = \partial^2_{u_1}[g\,u_1^2]\cdot f^2 = 2f^2 g\) with \(S = 2\), so the factor 2
*cancels* — the paper’s plain \(k(\tau_5)\) *(tex 6004–6011)*. The engine reproduces both.

## 7. Free constants, the renormalization group, and BPHZ

The family with free \(k(\tau)\) is the complete structural answer: the set of admissible
renormalizations forms a group \(G^-\) (characters on the negative trees under a convolution
product), and each group element gives one renormalized equation. The engine builds \(G^-\)
explicitly — see [the algebraic layer](guide/algebra.md).

Among all choices, the **BPHZ (canonical) character** is the distinguished one: it recentres each
divergent object by its expectation at the origin. Algebraically *(tex 5034)*:

$$k^\zeta = h^\zeta \circ S'_-,$$

where \(h^\zeta(\sigma) = \mathbb{E}\big[\Pi^\zeta \sigma\big](0)\) are the **elementary
expectations** (divergent Gaussian integrals — heat kernels against noise covariances) and
\(S'_-\) is the **negative twisted antipode**, a recursive Bogoliubov-style subtraction operator
on the tree Hopf algebra. The engine computes \(S'_-\) exactly and emits each canonical constant
as an exact polynomial in the symbols \(h(\sigma)\), with the integrals written out in
\(\varepsilon\)-regularized form — it does **not** evaluate them (that is the analytic wall; see
[Scope](scope.md)).

Several \(h(\sigma)\) vanish or coincide for provable, purely structural reasons, and the engine
knows these identities:

- **Gaussian parity**: an odd number of noise vertices ⇒ \(h = 0\) (Isserlis / Wick).
- **Root polynomial decoration**: \(h(X^n \tau) = 0\) for \(n \neq 0\) *(tex 5083)*.
- **Pure-kernel total derivative**: a noiseless tree with a derivative kernel leaf integrates to 0.
- **Spatial reflection** (only for a spatially symmetric noise): odd total spatial derivative
  order on the kernels ⇒ \(h = 0\). This is the identity that collapses KPZ’s canonical family to
  Hairer’s **single** diverging constant.

The first three are noise-independent and always safe; the reflection identity is gated behind an
explicit `symmetric` flag — see [Reports & output](guide/reports.md#canonical-and-reduced-views).

## 8. The full algebraic picture (Phase 3)

Behind the pipeline sits the two-Hopf-algebra structure of BHZ, all of which the engine builds and
property-tests over its equation corpus:

- the **recentering coproduct** \(\Delta : T \to T \otimes T^+\) (positive renormalization —
  Taylor recentring of the model), and the structure group;
- the **extraction–contraction coproduct** \(\delta^- \) (negative renormalization — cutting out
  divergent subforests, with the contracted node inheriting an extended decoration and the Taylor
  \(\mathfrak{e}'\) recentring terms);
- their **cointeraction** compatibility *(tex 5717)* — verified over the whole corpus **including
  the singular \(\beta_0 = -\tfrac32\)** case;
- the **twisted antipode** \(S'_-\) and the **renormalization group** \(G^-\).

None of this is needed to *state* the renormalized family (the free constants are), but it is
what makes the canonical character computable and the theory self-consistent, and it is exposed
programmatically — see [the algebraic layer](guide/algebra.md).

## Further reading

- The source paper, cited by tex line throughout: [arXiv:2006.03524](https://arxiv.org/abs/2006.03524).
- `notes/initial_plan.md` — the project’s authoritative mathematical conventions.
- `notes/free_constants_and_the_hopf_core.md` — why free constants are the right deliverable.
- `notes/elementary_expectations.md` — the \(h(\sigma)\) integrals and their reductions.
