import math
import sympy as s
import pytest
from recursive_horizons.nsc_atlas_volume_routes import calculation


def test_endpoint_data_leave_interpolating_barrier_free():
    z=s.symbols('z');f=z**4*(1-z)**4
    assert all(s.diff(f,z,n).subs(z,p)==0 for n in range(3) for p in (0,1))
    assert f.subs(z,s.Rational(1,2))>0


def test_all_five_attractors_and_same_action_decay():
    r=calculation()
    assert not r['physical_selection_derived']
    for row in r['frontier']:
        assert row['endpoint_distance_to_own_minimum']<3e-5
        assert abs(row['endpoint']['energy_ledger_residual'])<1e-8
        assert row['cross_candidate_bounce_to_12_2'] is None
        assert row['entropy_over_reference']<=1
    h=r['same_action_Hawking_Moss']
    assert h['exponent_times_g']==pytest.approx(403.7817410357,rel=1e-10)
    assert 0<h['negative_curvature_over_H_top_squared']<4
    assert h['exponent_times_g']==pytest.approx(24*math.pi**2*(1/h['false_vacuum']-1/h['barrier']))
