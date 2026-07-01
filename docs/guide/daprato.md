# The da Prato–Debussche lift

## The problem

The tree machinery requires **subcriticality**: \(\beta_0 > -\text{order}\). The dynamical
\(\Phi^4\) models violate it —

$$(\partial_t - \Delta)\,\varphi = -\varphi^3 + \xi,$$

with spacetime white noise has \(\beta_0 = -2-\kappa\) in \(d=2\) and \(\beta_0 = -\tfrac52-\kappa\)
in \(d=3\), both \(\le -2\). Feeding these to `renormalize()` raises a `ValueError` naming the fix.

## The change of variables

da Prato–Debussche: split off the rough part. Let \(X = L^{-1}\xi\) (the stochastic convolution)
and solve for the **remainder** \(v = u - X\). Substituting \(u = X + v\) *(tex 2026–2034)*:

$$(\partial_t - \Delta)\, v = -v^3 - 3v^2 X - 3v X^2 - X^3.$$

The powers \(X^k\) are ill-defined for white noise, but their **Wick versions** \(:\!X^k\!:\) are
honest random distributions of regularity \((k \cdot \alpha_X)^-\), where \(\alpha_X = \beta_0 +
\text{order}\) is the Schauder gain. The remainder equation is driven by the *family of noises*
\(:\!X\!:, :\!X^2\!:, :\!X^3\!:\), and its worst noise is regular enough that the equation is
**subcritical** — the standard machinery applies.

## Usage

```python
from counterterms import SPDE, Noise, Parabolic, Unknown, daprato_lift, kappa
from sympy import Rational

u  = Unknown("u", 3)
xi = Noise("xi", regularity=Rational(-5, 2) - kappa)     # Φ⁴₃
phi43 = SPDE(noises=[xi], operator=Parabolic(dim=3), unknown=u,
             rhs=xi.symbol - u.field ** 3)

lifted = daprato_lift(phi43)      # the v-equation, with noises X1, X2, X3 = :X:, :X²:, :X³:
eq = lifted.renormalize()         # the renormalized remainder equation
```

For \(\Phi^4_3\) the lift produces the remainder right-hand side
\(-v^3 - 3X_1 v^2 - 3X_2 v - X_3\) with \(X_k \in \mathcal{C}^{-k/2 - k\kappa}\) — matching the
paper’s display line for line — and the renormalized family reproduces the literature’s
structure: the divergences are **mass terms** (\(\propto v\), two independent constants in
\(d=3\)), the gradient-type terms vanish.

## What the lift does and does not claim

- **Applies to**: supercritical equations with **additive** noise and **polynomial**
  nonlinearity — the class where finitely many Wick powers absorb the roughness. Non-polynomial
  supercritical equations (sine-Gordon needs Wick *exponentials*) are rejected.
- **Wick powers as independent symbols.** The lifted noises \(X_1, X_2, X_3\) are treated as
  independent noises with the right regularities. The **free-constant family is exactly right**;
  but the *canonical* constants for \(\Phi^4\) would need the true joint law of the
  \(:\!X^k\!:\) (one underlying Gaussian), which the symbolic scheme does not model. Compare
  canonical \(\Phi^4\) output by structure, not by value.
- If the lift is insufficient (the remainder still supercritical), the downstream subcriticality
  check rejects it — you never get silent wrong output.
