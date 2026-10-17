"""Pure frame/unit and diagnostic checks; no physical evolution."""
import math

import pytest

from recursive_horizons.nsc_atlas_endpoint_comparison import (
    jordan_radius, refinement, scaled_radius, ssprk3_diagnostic_increment,
)


def test_frame_conversion_uses_the_invariant_trace():
    # Independent canonical parametrization of the exact frame change.
    radius, alpha, f = .5, .03, 1.4
    q = math.sqrt(alpha/6)*f
    c = radius*math.sinh(q)/math.sqrt(alpha/6)
    assert jordan_radius(radius,c,alpha) == pytest.approx(radius*math.cosh(q),rel=2e-15)
    assert jordan_radius(radius,c,0.) == radius
    assert jordan_radius(radius,-c,alpha) == jordan_radius(radius,c,alpha)


def test_transport_quadrature_matches_a_nonlinear_rk_update():
    # y'=y^2 makes the predictor and accepted-right slopes different.
    y0, h = 1., .1
    predictor = y0+h*y0*y0
    midpoint = .75*y0+.25*(predictor+h*predictor**2)
    accepted = y0/3+2*(midpoint+h*midpoint**2)/3
    slopes = [y0*y0, predictor**2, midpoint**2, accepted**2]
    increment = ssprk3_diagnostic_increment(h,slopes)
    assert increment == pytest.approx(accepted-y0,abs=3e-16)
    accepted_right_simpson = h*(slopes[0]+4*slopes[2]+slopes[3])/6
    assert abs(increment-accepted_right_simpson) > 1e-5


def test_observable_range_covers_all_grids_and_scales_linearly():
    radii = [.548, .549, .5492]
    result = scaled_radius(radii)
    assert result['value_metres'] == pytest.approx(5492.)
    assert result['resolution_half_width_metres'] == pytest.approx(12.)
    lo,hi = result['resolution_range_metres']
    assert all(lo <= v <= hi for v in result['grid_values_metres'])
    doubled = scaled_radius(radii,length_metres=20000.)
    assert doubled['value_metres'] == 2*result['value_metres']
    assert doubled['resolution_half_width_metres'] == 2*result['resolution_half_width_metres']
    assert result['observational_validation'] is False
    assert result['uncertainty']['rigorous_error_enclosure'] is False


def test_flat_or_worsening_gaps_are_not_reported_as_convergence():
    assert refinement([1.,.5,.2])['strictly_decreasing'] is True
    assert refinement([1.,.5,.6])['strictly_decreasing'] is False
    assert refinement([0.,0.,0.])['strictly_decreasing'] is False
    assert refinement([0.,0.,0.])['fine_over_coarse'] == [None,None]


@pytest.mark.parametrize('bad',[float('nan'),float('inf'),True,-1.])
def test_invalid_values_cannot_make_a_prediction(bad):
    with pytest.raises(ValueError): scaled_radius([.5,.5,bad])
    with pytest.raises(ValueError): refinement([1.,.5,bad])
