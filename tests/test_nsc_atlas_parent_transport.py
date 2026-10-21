import math
import pytest

from recursive_horizons.nsc_atlas_scalar_preparation import ScalarPreparation
from recursive_horizons.nsc_atlas_parent_transport import (
    RadialCapTransport, de_sitter_embedding_control, independent_child_arrival)


@pytest.fixture(scope="module")
def transport():
    model = ScalarPreparation()
    solution, _ = model.solve(tolerance=1e-9)
    return RadialCapTransport.from_preparation(model, solution)


def test_transport_recovers_independent_de_sitter_embedding():
    h, wall_angle, source_angle = .3, .4, .6
    a = lambda r: math.sin(h*r)/h
    ap = lambda r: math.cos(h*r)
    t = RadialCapTransport(wall_angle/h, (math.pi-wall_angle)/h,
                           a, a, ap, h, h*h)
    depth = (source_angle-wall_angle)/h
    for gamma in (1.4, 2., 5.):
        row = t.massive(depth, gamma)
        control = de_sitter_embedding_control(h, source_angle, gamma)
        assert row["arrival_time"] == pytest.approx(control["arrival_time"], abs=2e-9)
        assert row["gamma_received"] == pytest.approx(control["gamma_received"], abs=2e-10)
    assert t.central_threshold(depth) == pytest.approx(1/math.cos(source_angle), abs=2e-10)
    null = t.null(depth)
    assert null["one_plus_redshift"] == pytest.approx(1/math.cos(source_angle), abs=2e-10)


def test_child_trajectory_is_independent_of_analytic_inversion(transport):
    for gamma in (1.5, 2., 5.):
        row = transport.massive(.3, gamma)
        direct = independent_child_arrival(transport.child_h, row["radial_charge"], row["log_label"])
        assert direct == pytest.approx(row["arrival_time"], abs=2e-9)


def test_central_horizon_does_not_mean_no_child_entry(transport):
    failed = transport.massive(.3, 1.1)
    assert failed["status"] == "enters_child_but_misses_central_observer"
    assert "arrival_time" not in failed
    threshold = transport.central_threshold(.3)
    assert 1.29 < threshold < 1.30
    assert transport.massive(.3, threshold*1.001)["central_arrival"]
    assert not transport.massive(.3, threshold*.999)["central_arrival"]
    assert not transport.massive(.3, 2., eta=.5)["central_arrival"]


def test_null_limit_and_received_energy_not_envelope_time(transport):
    row = transport.massive(.3, 1e4)
    light = transport.null(.3)
    assert row["arrival_time"] == pytest.approx(light["arrival_time"], rel=1e-7)
    assert row["energy_ratio"] == pytest.approx(1/light["one_plus_redshift"], rel=1e-7)
    primary = transport.massive(.3, 2.)
    assert primary["gamma_received"]**2*(1-primary["speed_received"]**2) == pytest.approx(1.)
    assert primary["speed_received"] > .7  # arrival does not establish cold matter


def test_parent_collar_has_a_finite_null_horizon(transport):
    collar = transport.null_collar()
    d = collar["depth"]
    assert 0 < d < transport.parent_length
    assert .77 < collar["parent_slice_volume_fraction"] < .78
    assert transport.null(d-1e-4)["central_arrival"]
    assert not transport.null(d+1e-4)["central_arrival"]
    assert transport.null_collar(eta=2.)["depth"] == 0


def test_parent_turning_point_is_checked_before_arrival_formula(transport):
    # A source behind the parent's scale-factor maximum can turn on the way in.
    row = transport.massive(transport.parent_length-.2, 1.01)
    assert row["status"] == "turning_or_separatrix_before_child_cone"


def test_common_length_rescaling_preserves_dimensionless_observables(transport):
    t, scale = transport, 3.7
    scaled = RadialCapTransport(t.child_length*scale, t.parent_length*scale,
        lambda r: scale*t.child_a(r/scale), lambda r: scale*t.parent_a(r/scale),
        lambda r: t.parent_ap(r/scale), t.child_h/scale, t.pole_h_squared/scale**2)
    before, after = t.massive(.3, 2.), scaled.massive(.3*scale, 2.)
    assert after["arrival_time"]/scale == pytest.approx(before["arrival_time"], abs=2e-8)
    assert after["energy_ratio"] == pytest.approx(before["energy_ratio"], abs=2e-9)


def test_interface_has_no_impulsive_test_body_energy_jump(transport):
    # Continuous induced metric and no direct membrane interaction: a continuous
    # geodesic tangent has the same wall-frame gamma=K/R on both sides.
    t = transport
    R_i, R_o = t.child_a(t.child_length), t.parent_a(t.parent_length)
    K = t.massive(.3, 2.)["radial_charge"]
    assert K/R_i == pytest.approx(K/R_o, abs=1e-11)


def test_source_domain_rejects_nonphysical_inputs(transport):
    for depth, gamma in ((-.1, 2.), (.3, 1.), (.3, float("nan"))):
        with pytest.raises(ValueError):
            transport.massive(depth, gamma)
