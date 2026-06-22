"""da Prato–Debussche pre-pass — lift a supercritical additive-noise polynomial SPDE
into a subcritical one the engine can renormalise.

For ``(∂_t − L)u = ξ + P(u)`` with ``P`` a polynomial and ``β₀ ≤ −order``
(supercritical, so the direct rule reader rejects it), write ``u = X + v`` with
``X = K∗ξ`` the stochastic convolution — regularity ``α_X = β₀ + order``.  The noise
cancels (``(∂_t − L)X = ξ``), leaving ``(∂_t − L)v = P(X + v)``.  Taylor-expand in ``X``
and replace each power ``X^k`` (k≥1) by a **Wick power** ``:X^k:`` — a fresh noise of
regularity ``k·α_X`` (the Wick renormalisation tames the otherwise ill-defined powers):

    (∂_t − L)v = Σ_{k≥0} (P^{(k)}(v)/k!) :X^k:,      :X^0: = 1  (deterministic part P(v)).

The roughest new noise is ``:X^{deg P}:`` at ``deg·α_X``; one lift suffices when that
exceeds ``−order`` (e.g. Φ⁴₃: ``α_X=−1/2``, so ``:X³:∈C^{−3/2}>−2``).  The result is a
subcritical multi-noise SPDE for ``v`` — feed it to ``.renormalize()``.

**Scope.** Single equation, single *additive* noise (coefficient 1), polynomial
nonlinearity in ``u`` with no derivatives.  Multiplicative-noise and non-polynomial
nonlinearities (e.g. sine-Gordon, which needs Wick *exponentials*) are out of scope here.
The lifted Wick noises are correlated (all built from one ``ξ``); that is irrelevant to
the counterterm *structure* this produces, and matters only for the *canonical values*
(Phase 4 Track B).
"""
from __future__ import annotations

import sympy

from .dsl import SPDE, Noise, Unknown, kappa


def daprato_lift(spde: SPDE) -> SPDE:
    """Return the da Prato–Debussche-lifted SPDE for the remainder ``v = u − X``."""
    if len(spde.equations) != 1 or len(spde.noises) != 1:
        raise ValueError("da Prato–Debussche lift: single equation, single noise only")
    (u, op, rhs) = spde.equations[0]
    xi = spde.noises[0]
    rhs = sympy.expand(sympy.sympify(rhs))

    if rhs.coeff(xi.symbol, 1) != 1 or rhs.coeff(xi.symbol, 2) != 0:
        raise ValueError("da Prato–Debussche lift: needs additive noise "
                         f"(ξ with coefficient 1), got rhs {rhs}")
    P = sympy.expand(rhs - xi.symbol)                  # the noise-free part P(u)
    if P.has(sympy.Derivative):
        raise ValueError("da Prato–Debussche lift: polynomial nonlinearity in u only "
                         "(no ∂u)")
    if not P.is_polynomial(u.field):
        raise ValueError("da Prato–Debussche lift: nonlinearity must be polynomial in u "
                         "(sine-Gordon etc. need Wick exponentials — out of scope)")

    # X = K∗ξ gains `order` regularity: α_X = β₀ + order.
    a_std = xi.std + op.order
    a_kap = xi.kap
    v = Unknown("v", u.dim)
    X = sympy.Symbol("_X_")
    expanded = sympy.expand(P.subs(u.field, v.field + X))

    rhs_v = expanded.coeff(X, 0)                        # X⁰ → deterministic P(v)
    noises = []
    deg = int(sympy.degree(expanded, X)) if expanded.has(X) else 0
    for k in range(1, deg + 1):
        ck = expanded.coeff(X, k)
        if ck == 0:
            continue
        reg = k * (sympy.Rational(a_std.numerator, a_std.denominator)
                   + sympy.Rational(a_kap.numerator, a_kap.denominator) * kappa)
        nk = Noise(f"X{k}", regularity=reg)             # :X^k:, regularity k·α_X
        noises.append(nk)
        rhs_v += ck * nk.symbol

    return SPDE(noises=noises, operator=op, unknown=v, rhs=rhs_v)
