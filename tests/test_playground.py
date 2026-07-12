"""The standalone algebra playground: hand-built trees must agree with the pipeline's
canonical forms, and the Hopf layer must be reachable without an SPDE."""
from fractions import Fraction

from counterterms import playground


def test_playground_matches_pipeline_trees():
    from tests.conftest import kpz
    pg = playground("-3/2 - kappa")            # KPZ noise, d=1, parabolic

    t = pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"))          # ●·𝓘ₓ[Ξ]²
    assert str(t) == "●·𝓘ₓ[Ξ]²"
    h = pg.homogeneity(t)
    assert (h.std, h.kap) == (Fraction(-1), Fraction(-2))    # |τ| = −1 − 2κ
    assert t.symmetry_factor() == 2

    # the SAME object (canonical iso) as the pipeline's tree — equality and hashing
    eq = kpz().renormalize()
    assert t in set(eq.all_trees)

    # decorated variant: Xₓ·●·𝓘ₓ[Ξ]²
    td = pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"), X="x")
    assert td in set(eq.all_trees) and td.homogeneity(pg.sig).std == Fraction(0)


def test_playground_algebra_one_liners():
    import sympy
    pg = playground("-3/2 - kappa")
    t = pg.tree(pg.I(pg.Xi, "x"), pg.I(pg.Xi, "x"))

    assert pg.delta_minus(t).items()                          # coaction runs
    assert pg.antipode(t).items()                             # twisted antipode runs
    k = pg.character(t)                                       # k(τ) = h(S'₋τ), parity-reduced
    assert k != 0 and all(str(s).startswith("h") for s in k.free_symbols)
    assert pg.character(t) == k                               # h-symbols stable across calls
    e = pg.expectation(t)
    assert len(e.terms) == 1                                  # one Wick pairing of Ξ·Ξ

    # odd noise count ⇒ character and expectation vanish by parity
    single = pg.tree(pg.I(pg.Xi, "x"))
    assert pg.character(single) == 0 and pg.expectation(single).is_zero

    # multi-noise: named leaves required, cross-parity kills Ξ_ξ·𝓘[Ξ_η]
    pg2 = playground(noises={"xi": "-3/2 - kappa", "eta": "-1 - kappa"})
    mixed = pg2.tree(pg2.I(pg2.noise("xi")), pg2.I(pg2.noise("eta")))
    assert pg2.expectation(mixed).is_zero
    try:
        pg2.noise()
        assert False, "ambiguous noise must raise"
    except ValueError:
        pass
