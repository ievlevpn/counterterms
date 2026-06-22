"""Generic, basis-agnostic Hopf machinery (architecture.md §3.5).

The renormalisation theory's *algebra* — character convolution, the connected-graded
antipode, the comodule action — does not know what a symbol is.  It works over linear
combinations on any hashable keys, given the structural maps the concrete basis hands
in:

    coproduct : key -> {(left, right): coeff}      # a coalgebra / bialgebra Δ
    mul       : (key, key) -> key                  # the algebra product on basis keys
    unit      : key                                # the algebra unit 𝟙

`DecoratedTree` forests are the one concrete basis today; a multi-index basis (Phase 4)
would reuse this module unchanged.  The same functions serve all four coproducts
(`δ, δ⁻, Δ, Δ⁺`): e.g. `antipode` gives the structure-group antipode on `T⁺` (via `Δ⁺`)
and on `U⁻` (via `δ⁻`) with no per-basis code.

Coefficients are exact `Fraction`s for the antipode; characters may be scalar- or
symbol-valued (the convolution and comodule action only `+`/`*` the values).
"""
from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from typing import Callable, Hashable


def convolve(
    f: Callable[[Hashable], object],
    g: Callable[[Hashable], object],
    coproduct: Callable[[Hashable], dict],
) -> Callable[[Hashable], object]:
    """Character convolution ``(f⋆g)(x) = (f⊗g)Δx = Σ f(x⁽¹⁾)·g(x⁽²⁾)``.

    `f`, `g` are characters (key -> scalar/symbol); returns the character `f⋆g`.
    This is the group law on the character group of the Hopf algebra.
    """
    def conv(x: Hashable) -> object:
        return sum((c * f(l) * g(r) for (l, r), c in coproduct(x).items()), 0)
    return conv


def counit(unit: Hashable) -> Callable[[Hashable], object]:
    """The counit character ``ε``: 1 on the unit, 0 on every other basis key."""
    return lambda x: 1 if x == unit else 0


def antipode(
    coproduct: Callable[[Hashable], dict],
    mul: Callable[[Hashable, Hashable], Hashable],
    unit: Hashable,
) -> Callable[[Hashable], dict]:
    """The connected-graded antipode ``S`` (the convolution inverse of the identity,
    ``S⋆id = η∘ε``), as a memoised linear map ``key -> {key: Fraction}``.

    Recursion (Sweedler, on the reduced coproduct): ``S(𝟙)=𝟙`` and for ``x≠𝟙``

        S(x) = −x − Σ_{(x), x⁽¹⁾≠𝟙≠x⁽²⁾} S(x⁽¹⁾)·x⁽²⁾,

    which terminates because the reduced left factors have strictly lower degree
    (connectedness).
    """
    memo: dict = {}

    def S(x: Hashable) -> dict:
        if x == unit:
            return {unit: Fraction(1)}
        if x in memo:
            return memo[x]
        out: dict = defaultdict(Fraction, {x: Fraction(-1)})
        for (l, r), c in coproduct(x).items():
            if l == unit or r == unit:          # reduced coproduct only
                continue
            for fl, cl in S(l).items():
                out[mul(fl, r)] -= c * Fraction(cl)
        memo[x] = {k: v for k, v in out.items() if v}
        return memo[x]

    return S


def comodule_action(
    character: Callable[[Hashable], object],
    coaction: Callable[[Hashable], dict],
) -> Callable[[Hashable], dict]:
    """The comodule action ``k̃ = (k⊗Id)∘δ`` of a character `k` through a coaction
    ``δ : M → C ⊗ M`` — maps ``x ↦ Σ k(x^{(left)})·x^{(right)}`` (a linear combination
    in M).  This is how a renormalisation character acts on the model.
    """
    def act(x: Hashable) -> dict:
        out: dict = defaultdict(int)
        for (left, right), c in coaction(x).items():
            out[right] += c * character(left)
        return {k: v for k, v in out.items() if v}
    return act
