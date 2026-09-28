"""Continuous Taylor remainder and same-family profile connection."""
import pytest
pytest.importorskip('flint')
from flint import arb, ctx

from recursive_horizons.nsc_ks_ball_operator import (
    AnalyticRadiusFamily, time_operator_enclosure, slab_ball, spatial_potential_jets,
)
from recursive_horizons.nsc_ks_ball_trajectory import polynomial_at
from test_nsc_ks_difference_envelope import inputs


def test_dyadic_boundary_ball_contains_the_actual_final_slab():
    with ctx.workprec(120):
        left, right = arb(1), arb(1.0005)
        enclosure = slab_ball(left, right)
        assert enclosure >= 1
        assert enclosure.contains(left) and enclosure.contains(right)


def test_time_polynomial_remainders_contain_independent_point_jets():
    model = AnalyticRadiusFamily(inputs()[2])
    with ctx.workprec(120):
        first, last = arb(1.005), arb(1.0045)
        polynomials, errors = time_operator_enclosure(model, first, last, 2., degree=6)
        for point in (arb(0), arb(1)/4, arb(1)/2, arb(1)):
            exact = model.time_series(first+(last-first)*point, 2., 0, 4)
            for name, coeff in polynomials.items():
                difference = polynomial_at(coeff,point)-exact[name][0]
                # Ball radii of both independently evaluated representations
                # are retained when comparing with the Taylor remainder.
                bound = errors[name]+difference.real.rad()+difference.imag.rad()
                assert difference.abs_upper() <= 2*bound


def test_spatial_powers_retain_the_actual_amplitude_and_mixed_profile():
    model = AnalyticRadiusFamily(inputs()[2])
    with ctx.workprec(120):
        z = arb(float(model.metric.directions[0].w.center))
        direct = model.profile_series(z, 2)
        jets = spatial_potential_jets(model,z,[(1,0),(0,1),(1,1)],2)
        assert jets[(1,0)][0].overlaps(direct[0][0])
        assert jets[(0,1)][0].overlaps(direct[1][0])
        assert jets[(1,1)][1].overlaps(direct[0][1]*direct[1][0]+direct[0][0]*direct[1][1])
