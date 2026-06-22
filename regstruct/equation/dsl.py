"""User-facing DSL and the parser SPDE → (Signature, per-equation nonlinearities).

The user writes one or more equations ``L_a u_a = rhs_a`` with SymPy, tagging the
unknown(s) / noise(s) / operator(s).  The package *derives* the structural rule
from the monomials of each ``rhs_a`` (mirroring tourist_guide.tex 5306–5340); the
user never writes trees or rules.  Scalar = one component; systems = several,
sharing spacetime coordinates and coupling through the nonlinearities.

Scope is enforced here with explicit errors (Assumption D2 etc.): affine-in-noise,
``g`` at most quadratic in ∂u, ``|p|_𝔰 ≤ 1``, ``β₀ ∈ (−2,0)``.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

import sympy

from ..core.homogeneity import Homogeneity, Scaling
from ..core.jets import is_jet, jet, jet_parts
from ..core.signature import Signature
from .rule import check_subcritical

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
    """A solution component. Components share spacetime coords (same `dim`)."""

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
    """``∂_t − Δ (+ mass)`` by default (order 2). `order` is carried into the
    homogeneity arithmetic; only order 2 is the analytically proven regime."""

    def __init__(self, dim: int, mass=0, order: int = 2):
        self.dim = dim
        self.mass = mass
        self.order = order
        self.scaling = Scaling(tuple([order] + [1] * dim))
        self.label = "I"
        if order != 2:
            import warnings
            warnings.warn(
                f"Schauder/admissibility is proven only for 2nd-order parabolic L; "
                f"homogeneities are computed for order={order} but the regularity-"
                f"structure theory is unverified there.", stacklevel=2)


@dataclass
class SPDE:
    equations: list   # list of (Unknown, Parabolic, rhs)
    noises: list

    def __init__(self, noises, operator=None, unknown=None, rhs=None, equations=None):
        if equations is None:
            equations = [(unknown, operator, rhs)]
        self.equations = equations
        self.noises = noises

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


def _to_jet(expr, field_to_comp, coords):
    expr = sympy.sympify(expr)
    for d in list(expr.atoms(sympy.Derivative)):
        if d.expr in field_to_comp:
            comp = field_to_comp[d.expr]
            expr = expr.xreplace({d: jet(comp, _deriv_index(d, coords))})
    subs = {field: jet(comp, (0,) * len(coords)) for field, comp in field_to_comp.items()}
    return expr.xreplace(subs)


def build_context(spde: SPDE):
    equations = spde.equations
    noises = spde.noises
    ncomp = len(equations)
    coords = equations[0][0].coords
    width = len(coords)
    field_to_comp = {eqn[0].field: a for a, eqn in enumerate(equations)}
    scaling = equations[0][1].scaling          # global scaling
    comp_order = tuple(op.order for (_u, op, _r) in equations)

    # per-equation base nonlinearities, in jet variables
    base: dict[int, dict[str, object]] = {}
    for a, (_u, _op, rhs_a) in enumerate(equations):
        rhs = sympy.expand(sympy.sympify(rhs_a))
        ba: dict[str, object] = {}
        g = rhs
        for nz in noises:
            if sympy.expand(rhs).coeff(nz.symbol, 2) != 0:
                raise ValueError(f"equation {a}: must be affine in noise '{nz.name}'")
            coeff = rhs.coeff(nz.symbol, 1)
            ba[nz.name] = _to_jet(coeff, field_to_comp, coords)
            g = g - coeff * nz.symbol
        g = sympy.expand(g)
        for nz in noises:
            if g.has(nz.symbol):
                raise ValueError(f"equation {a}: noise-free part still has '{nz.name}'")
        ba["bullet"] = _to_jet(g, field_to_comp, coords)
        base[a] = ba

    # scope checks (the β₀ lower bound is enforced rule-wise by check_subcritical below)
    for nz in noises:
        if not nz.homogeneity.is_negative():       # ordered-ring: −κ (std 0, kap<0) IS singular
            raise ValueError(f"noise '{nz.name}' regularity {nz.homogeneity} is not singular "
                             "(this package renormalises noises of negative regularity)")
    for a in range(ncomp):
        for nz in noises:
            if any(is_jet(s) and any(jet_parts(s)[1]) for s in base[a][nz.name].free_symbols):
                raise ValueError(f"equation {a}: noise coefficient '{nz.name}' must depend on "
                                 "u only (derivative noise out of scope)")
        gjets = [s for s in base[a]["bullet"].free_symbols if is_jet(s) and any(jet_parts(s)[1])]
        for s in gjets:
            if scaling.scaled(jet_parts(s)[1]) > 1:
                raise ValueError(f"equation {a}: singular derivative factor |p|_𝔰>1: {s}")
        if gjets and sympy.Poly(base[a]["bullet"], *gjets).total_degree() > 2:
            raise ValueError(f"equation {a}: g must be at most quadratic in ∂u (Assumption D2)")

    node_types = ("bullet",) + tuple(nz.name for nz in noises)

    # structural rule: per node type, the child edges (component, p) the nonlinearity
    # depends on, unioned over equations; cap = degree for derivative slots, None for fields.
    allowed: dict[str, tuple] = {}
    for b in node_types:
        caps: dict[tuple[int, tuple], "int | None"] = {}
        for a in range(ncomp):
            fb = base[a][b]
            for s in fb.free_symbols:
                if not is_jet(s):
                    continue
                comp, p = jet_parts(s)
                if any(p):
                    d = int(sympy.degree(fb, s))
                    caps[(comp, p)] = max(caps.get((comp, p), 0), d)
                else:
                    caps[(comp, p)] = None
        allowed[b] = tuple((comp, p, cap) for (comp, p), cap in caps.items())

    sig = Signature(
        dim=equations[0][0].dim,
        scaling=scaling,
        n_components=ncomp,
        comp_order=comp_order,
        noise_homog={nz.name: nz.homogeneity for nz in noises},
        node_types=node_types,
        allowed=allowed,
    )
    check_subcritical(sig)        # rule must be subcritical, else 𝓑_{<0} is infinite
    unknowns = [eqn[0] for eqn in equations]
    return sig, base, unknowns
