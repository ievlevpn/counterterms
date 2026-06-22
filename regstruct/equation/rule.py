"""Subcriticality of the structural rule (BHZ §5.5 / tourist_guide.tex 5302–5340).

A rule is **subcritical** iff only finitely many strongly-conforming trees sit below
any homogeneity γ — so tree generation terminates and ``𝓑_{<0}`` is finite.  The
obstruction is an *uncapped* child slot (a field slot ``∂^p u^c`` whose nonlinearity
dependence has no multiplicity bound): a node may carry unboundedly many such
children, so the family stays finite only if each one *strictly raises* the standard
homogeneity, even over the most singular subtree — a bare noise ``β₀``.  Hence

    subcritical  ⟺  for every uncapped edge,  (m_c − |p|_𝔰) + min_ℓ β₀^ℓ  >  0.

For 2nd-order parabolic (``m=2``, field slot ``|p|=0``) this is ``β₀ > −2``; ``β₀ ≤ −2``
(Φ⁴₃, sine-Gordon) is supercritical and needs a da Prato–Debussche pre-pass.  Higher
operator order ``m`` relaxes the threshold to ``β₀ > −m`` — the reason this is a
rule-based check, not a hardcoded ``−2``.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..core.signature import Signature


def check_subcritical(sig: Signature) -> None:
    """Raise ``ValueError`` if the rule carried by `sig` is not subcritical."""
    if not sig.noise_homog:
        return
    beta_min = min(h.std for h in sig.noise_homog.values())   # most singular subtree
    for node_type, rules in sig.allowed.items():
        for (comp, p, cap) in rules:
            if cap is not None:                  # capped (derivative) slot: bounded width
                continue
            gain = sig.edge_gain(comp, p).std    # m_c − |p|_𝔰
            if gain + beta_min <= 0:
                raise ValueError(
                    f"rule is not subcritical: the unbounded field slot "
                    f"I^({comp})_{tuple(p)} (gain {gain}) over the most singular "
                    f"subtree (β₀={beta_min}) contributes {gain + beta_min} ≤ 0, so "
                    f"infinitely many trees lie below any homogeneity. "
                    f"β₀ ≤ −order is supercritical — needs a da Prato–Debussche pre-pass.")
