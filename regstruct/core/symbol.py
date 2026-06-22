"""The ``Symbol`` seam — the atom of the graded free module the whole algebra is written
against (architecture §3.3).

``DecoratedTree`` is the one concrete implementation today; a *second basis* (the
multi-index / Linares–Otto–Tempelmayr approach, Phase 4 A3) only has to satisfy this
protocol — `core/hopf` and the coproduct-level code never touch a tree directly, they
require just a homogeneity and a canonical key.

[SOCKET — protocol only. The multi-index implementation is not built; see
notes/phase4_plan.md A3.]
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Hashable, Protocol, runtime_checkable

from .homogeneity import Homogeneity

if TYPE_CHECKING:
    from .signature import Signature


@runtime_checkable
class Symbol(Protocol):
    """A basis atom: it knows its homogeneity and a canonical key for equality/hashing."""

    def homogeneity(self, sig: Signature) -> Homogeneity:
        ...

    def canonical_key(self) -> Hashable:
        ...
