"""Check incoming sensitivity against independent constrained-cap evaluations."""
import math

import numpy as np
import pytest

from diagnose_nsc_atlas_endpoint_projection import incoming_offset_response
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_reference import fixed_cap_reference
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


def test_incoming_change_matches_two_independent_cap_evaluations():
    # C3 piecewise quartic crosses the single outer trial element, so the
    # background time jet has a nonzero finite wall-projection defect.
    q=np.array([-.8,0.,.8]);positive=np.maximum(q,0)
    background=RetainedCapBackground(q,.1*positive**4,.4*positive**3,1.2*positive**2)
    budget=Budget(maximum_evaluations=5_000_000,cpu_seconds=30)
    basis=FixedCapBasis([0,.05,.85],budget=budget)
    zeros=np.zeros(basis.dimension);r=.45;Z=math.atanh(r);alpha=.001
    jet=background.reduced_jets(.2,Z,budget=budget)
    incoming=(jet['p']+jet['h'])/r
    def evaluate(value):
        row=fixed_cap_reference(basis,zeros,zeros,Z,alpha=alpha,T_x=1.,radius_rate=0.,
            incoming_w_plus=value,budget=budget,background=background,background_time=.2,
            background_panel_sides=True,rtol=2.5e-14,atol=1e-15)
        w=row['wall']
        trace=w['A']**2*math.exp(row['lapse'])*(w['p_weak']*w['h_weak']-w['p_trace']*w['h_trace'])/2
        return row,trace
    first,t0=evaluate(incoming)
    delta=.17
    second,t1=evaluate(incoming+delta/r)
    assert abs(first['background']['wall_time_projection_gap'])>1e-6
    w=first['wall']
    response=incoming_offset_response(A=w['A'],lapse=first['lapse'],p_weak=w['p_weak'],
        time_wall_gap=first['background']['wall_time_projection_gap'],delta_incoming_sum=delta,alpha=alpha)
    assert second['wall']['p_weak']==pytest.approx(w['p_weak'],abs=2e-13)
    assert second['wall']['h_weak']-w['h_weak']==pytest.approx(delta,abs=2e-13)
    assert t1-t0==pytest.approx(response['normalized_trace_change'],abs=2e-12)
    assert second['explicit_time_remainder']-first['explicit_time_remainder']==pytest.approx(
        response['normalized_affine_change'],abs=2e-12)
    assert alpha*(t1+second['explicit_time_remainder']-t0-first['explicit_time_remainder'])==pytest.approx(
        response['physical_mass_gap_change'],abs=2e-14)


def test_affine_response_adds_and_has_a_zero_offset_limit():
    kwargs=dict(A=.7,lapse=.01,p_weak=-3.,time_wall_gap=2.,alpha=.0001)
    whole=incoming_offset_response(**kwargs,delta_incoming_sum=.5)
    parts=[incoming_offset_response(**kwargs,delta_incoming_sum=v) for v in (.2,.3)]
    for key in whole: assert whole[key]==pytest.approx(sum(p[key] for p in parts))
    assert all(v==0 for v in incoming_offset_response(**kwargs,delta_incoming_sum=0).values())


def test_invalid_metric_cannot_be_used_for_boundary_sensitivity():
    with pytest.raises(ValueError):
        incoming_offset_response(A=-.1,lapse=0,p_weak=1,time_wall_gap=1,delta_incoming_sum=1,alpha=.1)
