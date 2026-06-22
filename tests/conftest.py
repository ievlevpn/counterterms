"""Shared test corpus and fixtures.

A small **corpus of in-scope SPDEs** that, between them, exercise every generality
axis of the library:

    name        components  noises  dim  β₀        gradient nl   notes
    ----------  ----------  ------  ---  --------  ------------  -----------------------
    gkpz        1           1       1    −1−κ      yes (g·(∂u)²) the running example
    kpz         1           1       1    −3/2−κ    yes           the singular case
    gpam1       1           1       1    −1−κ      no            pure f(u)ξ
    gpam2       1           1       2    −1−κ      no            two space dimensions
    system      2 (coupled) 1       1    −1−κ      yes (in u)    shared noise, coupling
    multinoise  1           2       1    −1−κ      no            two independent noises

The invariant tests (`test_coproducts.py`, …) parametrize over this corpus via the
`ctx` fixture, so the algebra is checked **systematically across the whole input
class**, not on a single example.  End-to-end golden tests call the factories
directly and `.renormalize()`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import pytest
from sympy import Derivative, Function, Rational

from counterterms import Noise, Parabolic, SPDE, Unknown, kappa
from counterterms.equation.dsl import build_context
from counterterms.equation.generate import generate_counterterms

f, g, h = Function("f"), Function("g"), Function("h")
a_coef, b_coef = Function("a"), Function("b")


def gkpz(reg=Rational(-1) - kappa) -> SPDE:
    """Generalised KPZ — scalar, one noise, d=1, ``f(u)ξ + g(u)(∂ₓu)²``."""
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=reg)
    return SPDE(noises=[xi], operator=Parabolic(dim=1, mass=1), unknown=u,
                rhs=f(u.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2)


def kpz() -> SPDE:
    """Classical KPZ — d=1 space-time white noise (β₀=−3/2−κ), ``ξ + (∂ₓu)²``."""
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=Rational(-3, 2) - kappa)
    return SPDE(noises=[xi], operator=Parabolic(dim=1), unknown=u,
                rhs=xi.symbol + Derivative(u.field, u.x[0]) ** 2)


def gpam1() -> SPDE:
    """gPAM in d=1 — ``f(u)ξ`` (no gradient nonlinearity)."""
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    return SPDE(noises=[xi], operator=Parabolic(dim=1), unknown=u, rhs=f(u.field) * xi.symbol)


def gpam2() -> SPDE:
    """gPAM in d=2 — ``f(u)ξ``, space white noise."""
    u = Unknown("u", 2)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    return SPDE(noises=[xi], operator=Parabolic(dim=2), unknown=u, rhs=f(u.field) * xi.symbol)


def system() -> SPDE:
    """Coupled 2-component system sharing one noise: ``u̇=a(v)ξ+g(u)(∂ₓu)²``, ``v̇=b(u)ξ``."""
    u, v = Unknown("u", 1), Unknown("v", 1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    op = Parabolic(dim=1, mass=1)
    return SPDE(noises=[xi], equations=[
        (u, op, a_coef(v.field) * xi.symbol + g(u.field) * Derivative(u.field, u.x[0]) ** 2),
        (v, op, b_coef(u.field) * xi.symbol),
    ])


def multinoise() -> SPDE:
    """Scalar with two independent noises: ``f(u)ξ + h(u)η``."""
    u = Unknown("u", 1)
    xi = Noise("xi", regularity=Rational(-1) - kappa)
    eta = Noise("eta", regularity=Rational(-1) - kappa)
    return SPDE(noises=[xi, eta], operator=Parabolic(dim=1, mass=1), unknown=u,
                rhs=f(u.field) * xi.symbol + h(u.field) * eta.symbol)


CORPUS = {"gkpz": gkpz, "kpz": kpz, "gpam1": gpam1, "gpam2": gpam2,
          "system": system, "multinoise": multinoise}


@dataclass
class Ctx:
    """A built SPDE context: the `Signature`, per-equation base nonlinearities,
    unknowns, and (lazily, cached) the generated divergent-tree set."""

    name: str
    sig: object
    base: dict
    unknowns: list
    _trees: list = field(default=None, repr=False)

    def trees(self, max_nodes: int = 10**9):
        """The divergent trees |τ|<0, optionally capped by node count (for speed)."""
        if self._trees is None:
            self._trees = generate_counterterms(self.sig)
        return [t for t in self._trees if t.nodes() <= max_nodes]


@pytest.fixture(params=list(CORPUS), scope="session")
def ctx(request) -> Ctx:
    """One built context per corpus equation (session-cached). Parametrized: a test
    using `ctx` runs once for each SPDE in the corpus."""
    sig, base, unknowns = build_context(CORPUS[request.param]())
    return Ctx(request.param, sig, base, unknowns)
