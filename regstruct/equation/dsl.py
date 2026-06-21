"""User-facing DSL and the parser SPDE → (Signature, base nonlinearities).

The user writes ``L u = rhs`` with SymPy, tagging the unknown / noise / operator.
The package *derives* the structural rule from the monomials of ``rhs`` (mirroring
tourist_guide.tex 5306–5340); the user never writes trees or rules.

Scope is enforced here with explicit errors (Assumption D2 etc.):
affine-in-noise, g at most quadratic in ∂u, |p|_𝔰 ≤ 1, β₀ ∈ (−2,−1).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction

import sympy

from ..core.homogeneity import Homogeneity, Scaling
from ..core.jets import is_jet, jet, jet_index
from ..core.signature import Signature

kappa = sympy.Symbol("kappa", positive=True)


def _frac(x) -> Fraction:
    r = sympy.Rational(x)
    return Fraction(int(r.p), int(r.q))


def _split_kappa(expr) -> tuple[Fraction, Fraction]:
    expr = sympy.expand(sympy.sympify(expr))
    if expr.has(kappa) and sympy.Poly(expr, kappa).degree() > 1:
        raise ValueError("regularity must be affine in κ")
    return _frac(expr.coeff(kappa, 0)), _frac(expr.coeff(kappa, 1))


class Unknown:
    def __init__(self, name: str, dim: int):
        self.name = name
        self.dim = dim
        self.t = sympy.Symbol("t")
        self.x = [sympy.Symbol(f"x{i + 1}") for i in range(dim)]
        self.coords = (self.t, *self.x)
        self.field = sympy.Function(name)(*self.coords)


class Noise:
    def __init__(self, name: str, regularity):
        self.name = name
        self.symbol = sympy.Symbol(name)
        self.std, self.kap = _split_kappa(regularity)
        self.homogeneity = Homogeneity(self.std, self.kap)


class Parabolic:
    """Second-order parabolic operator ``∂_t − Δ (+ mass)``; the proven theory."""

    def __init__(self, dim: int, mass=0):
        self.dim = dim
        self.mass = mass
        self.scaling = Scaling(tuple([2] + [1] * dim))
        self.order = 2
        self.label = "I"


@dataclass
class SPDE:
    operator: Parabolic
    unknown: Unknown
    noises: list
    rhs: object

    def renormalize(self):
        from ..api import renormalize
        return renormalize(self)


# --------------------------------------------------------------------------- #
# parsing
# --------------------------------------------------------------------------- #

def _deriv_index(d: sympy.Derivative, coords) -> tuple[int, ...]:
    counts = {c: 0 for c in coords}
    for var, cnt in d.variable_count:
        counts[var] += int(cnt)
    return tuple(counts[c] for c in coords)


def _to_jet(expr, u: Unknown):
    expr = sympy.sympify(expr)
    for d in list(expr.atoms(sympy.Derivative)):
        if d.expr == u.field:
            k = _deriv_index(d, u.coords)
            expr = expr.xreplace({d: jet(k)})
    expr = expr.xreplace({u.field: jet((0,) * len(u.coords))})
    return expr


def build_context(spde: SPDE):
    u, noises = spde.unknown, spde.noises
    rhs = sympy.expand(sympy.sympify(spde.rhs))
    width = 1 + u.dim
    zero = (0,) * width

    base: dict[str, object] = {}
    g = rhs
    for nz in noises:
        if sympy.expand(rhs).coeff(nz.symbol, 2) != 0:
            raise ValueError(f"nonlinearity must be affine in noise '{nz.name}'")
        coeff = rhs.coeff(nz.symbol, 1)
        base[nz.name] = _to_jet(coeff, u)
        g = g - coeff * nz.symbol

    g = sympy.expand(g)
    for nz in noises:
        if g.has(nz.symbol):
            raise ValueError("noise-free part still depends on a noise (not affine?)")
    base["bullet"] = _to_jet(g, u)

    # scope checks + structural rule derivation
    for nz in noises:
        derivs = [s for s in base[nz.name].free_symbols if is_jet(s) and any(jet_index(s))]
        if derivs:
            raise ValueError(f"noise coefficient '{nz.name}' must depend on u only (derivative noise out of scope)")
        if not (-2 < nz.std < 0):
            raise ValueError(f"noise '{nz.name}' regularity {nz.std} out of MVP scope β₀∈(−2,0); "
                             "β₀≤−2 needs a da Prato–Debussche pre-pass")

    node_types = tuple(nz.name for nz in noises) + ("bullet",)
    allowed: dict[str, tuple] = {}
    for b in node_types:
        fb = base[b]
        rules = [(spde.operator.label, zero, None)]  # I_0 from the smooth coefficient
        derivs = [s for s in fb.free_symbols if is_jet(s) and any(jet_index(s))]
        for s in derivs:
            idx = jet_index(s)
            if spde.operator.scaling.scaled(idx) > 1:
                raise ValueError(f"singular derivative factor with |p|_𝔰>1 out of scope: {s}")
            rules.append((spde.operator.label, idx, int(sympy.degree(fb, s))))
        if derivs:
            totdeg = sympy.Poly(fb, *derivs).total_degree()
            if totdeg > 2:
                raise ValueError("g must be at most quadratic in ∂u (Assumption D2)")
        allowed[b] = tuple(rules)

    sig = Signature(
        dim=u.dim,
        scaling=spde.operator.scaling,
        op_order=spde.operator.order,
        op_label=spde.operator.label,
        noise_homog={nz.name: nz.homogeneity for nz in noises},
        node_types=node_types,
        allowed=allowed,
    )
    return sig, base, width
