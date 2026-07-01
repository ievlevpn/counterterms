# Conventions

Regularity structures are a minefield of near-equivalent conventions, and getting one wrong is
*silent and fatal* — the output looks plausible and is wrong. These are the project’s
non-negotiable choices, each pinned by tests. If you interoperate with the engine (consume its
JSON, compare against another implementation, extend it), read this page first.

## 1. \(|\Xi| = \beta_0\), directly

The `regularity` you pass to `Noise` **is** the homogeneity of the noise symbol: pass
\(-\tfrac32 - \kappa\) for 1-d spacetime white noise, and the engine uses exactly that
*(tex 5276–5283)*. Some texts parametrize noise by the regularity \(\alpha\) of its *integrated*
object, related by \(\beta_0 = \alpha - 2\); mixing the two conventions is an off-by-two that
shifts every homogeneity in every table. The engine has no such shift anywhere:

$$|\Xi| = \beta_0, \qquad |\mathcal{I}_p \tau| = |\tau| + \text{order} - |p|_\mathfrak{s},
\qquad |X^n| = |n|_\mathfrak{s}.$$

## 2. Homogeneities live in \(\mathbb{Q} \oplus \mathbb{Q}\kappa\), never floats

\(\kappa\) is a formal positive infinitesimal; a homogeneity is a pair
\((\text{std}, \text{kap}) \in \mathbb{Q}^2\) ordered lexicographically. Critical trees sit at
\(-k\kappa\): standard part zero, infinitesimally negative. They **are** in
\(\mathcal{B}_{<0}\) (they diverge logarithmically), and exact zero is **not**. Any float
representation would round them the wrong way; the engine compares only in the ordered ring.

## 3. The coefficient is \(k(\tau)/S(\tau)\), never \(k(\tau)\) alone

The BCCH formula divides each constant by the tree’s symmetry factor
\(S(\tau) = n!\,\prod_j S(\sigma_j)^{m_j} m_j!\) *(tex 3982, 4915)*. The paper’s displayed
coefficients absorb this division (which is why gKPZ shows a bare \(k(\tau_5)\) even though
\(\Upsilon(\tau_5^*) = 2f^2g\)); the engine keeps \(k\), \(S\), and \(F(\tau^*)\) separate and
exposes the ratio as `Counterterm.coefficient`.

## 4. Canonical tree isomorphism is load-bearing

Tree equality, dict keys, deduplication, and \(S(\tau)\) all reduce to a recursive canonical key:

```
(node_type, n, sorted (edge_type, p, canonical(child)))
```

Two decorated trees are equal iff their keys are equal. If you construct trees by hand, never
compare them any other way.

## 5. The Υ-map differentiates *all* slots, in the right order

\(F(\tau^*) = (\prod_i F(\tau_i^*)) \cdot (D^n \prod_i \partial_{p_i}) F(b^*)\) with base cases
\(F(\circ_j^*) = f_j\), \(F(\bullet^*) = g\), \(F(\text{red}^*) = 0\) *(tex 4337)*. The slot
derivatives \(\partial_{p_i}\) act on **every** argument of \(g\) — the function slot *and* the
derivative slots — and are applied **before** the total derivatives
\(D_i = \sum_k u_{k + e_i}\partial_{u_k}\); the two operations do not commute, and the order is
tested.

## 6. \(\mathcal{B}_{<0}\) includes the bare noises

The counterterm sum runs over *all* negative trees, including the primitive symbols
\(\circ^n\) — in KPZ these are the \(k(\circ)\) (constant shift) and \(k(\circ_1)\) (drift)
terms that a naive “only branched trees diverge” reading would drop. The bare \(\bullet\)
(homogeneity 0) is excluded.

## 7. For systems, component identity rides on the *edge*

In a coupled system, node types stay \(\{\bullet, \circ_j, \text{red}\}\); *which equation* a
kernel integrates against is recorded on the **edge type** \(\mathfrak{T}_e\) *(tex 3826–3827)*.
Consequently one tree can contribute counterterms to several components, and the constant
\(k(\tau)\) is **shared** across them — the engine emits one symbol per tree, not per
(tree, component) pair.

## 8. Parabolic scaling: time counts double

\(\mathfrak{s} = (\text{order}, 1, \dots, 1)\), so \(|n|_\mathfrak{s} = \text{order}\cdot n_0 +
n_1 + \dots + n_d\). A time derivative on a kernel edge costs the full operator order; the
engine’s `Scaling` carries this and every homogeneity computation threads through it.

---

Each convention above is enforced by at least one unit test in `tests/test_homogeneity.py`,
`tests/test_trees.py`, or `tests/test_pipeline.py`, and the two paper goldens
(see [Validation](validation.md)) would catch a violation of any of them.
