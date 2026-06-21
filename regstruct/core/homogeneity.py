"""Scaled degrees and the ordered ring of homogeneities.

A homogeneity is an element ``std + kap·κ`` of ``ℚ ⊕ ℚ·κ`` where ``κ`` is a
*positive infinitesimal* (the place-holder for the small regularity loss in the
noise).  Comparison is lexicographic: standard part first, then the κ-coefficient
(κ>0).  We never use floats — critical trees sit exactly at homogeneity ``−kκ``
(tourist_guide.tex 6066, 6140), and a float β₀ would mis-classify them.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

MultiIndex = tuple[int, ...]


@dataclass(frozen=True)
class Scaling:
    """The parabolic scaling 𝔰 (e.g. ``(2,1,...,1)``) and its scaled degree."""

    weights: tuple[int, ...]

    def scaled(self, k: MultiIndex) -> int:
        return sum(w * ki for w, ki in zip(self.weights, k))


@dataclass(frozen=True)
class Homogeneity:
    std: Fraction
    kap: Fraction = Fraction(0)

    def __post_init__(self):
        object.__setattr__(self, "std", Fraction(self.std))
        object.__setattr__(self, "kap", Fraction(self.kap))

    def __add__(self, other: "Homogeneity") -> "Homogeneity":
        return Homogeneity(self.std + other.std, self.kap + other.kap)

    def __neg__(self) -> "Homogeneity":
        return Homogeneity(-self.std, -self.kap)

    def is_negative(self) -> bool:
        return self.std < 0 or (self.std == 0 and self.kap < 0)

    def _key(self):
        return (self.std, self.kap)

    def __lt__(self, other: "Homogeneity") -> bool:
        return self._key() < other._key()

    def __str__(self) -> str:
        parts = []
        if self.std != 0 or self.kap == 0:
            parts.append(str(self.std))
        if self.kap != 0:
            sign = "+" if self.kap > 0 else "-"
            mag = abs(self.kap)
            coeff = "" if mag == 1 else str(mag)
            term = f"{coeff}κ"
            parts.append(f"{sign} {term}" if parts else (f"-{term}" if self.kap < 0 else term))
        return " ".join(parts) if parts else "0"

    __repr__ = __str__
