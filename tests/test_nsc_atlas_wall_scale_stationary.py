import pytest
from check_nsc_atlas_wall_scale_stationary import calculation


def test_inherited_tension_stationary_points_are_not_static_stable_solutions():
    for row in calculation()['walls']:
        assert abs(row['slope'])<2e-13
        assert row['second_derivative']==pytest.approx(row['independent_stationary_identity'],abs=2e-13)
        assert row['second_derivative']<0
        assert abs(row['potential_at_stationary_radius'])>.01
        assert row['original_velocity']>0
