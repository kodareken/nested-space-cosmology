import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_causal_reference import (
    causal_reference_record, initial_slice_curve, parent_event_from_U,
    reference_wall_events, schwarzschild_slice_time, signed_state,
)


def test_imported_branch_has_no_turning_point_and_keeps_signs():
    r = causal_reference_record()["classification"]
    assert r["q_minimum"] > .125
    assert 0 < r["beta_plus_zero_over_ell"] < .8 < 1 < 1.5
    for x in (.02, .5, .8, 1., 1.5, 6.):
        s = signed_state(x)
        assert s["beta_minus"] > 0
        assert s["beta_minus"]-s["beta_plus"] == pytest.approx(2*x)
        assert s["beta_plus"]**2-s["Rdot"]**2 == pytest.approx(s["f_plus"], abs=1e-10)


def test_horizon_crossing_regular_future_and_unit_timelike():
    a = reference_wall_events(np.array([.02, .5, .8-1e-8, .8, .8+1e-8, 1., 1.5, 6.]))
    assert np.all(a["ell_dU_dtau"] > 0)
    assert np.all(a["ell_dV_dtau"] > 0)
    assert np.max(abs(a["metric_norm"]+1)) < 1e-9
    assert a["U"][3] == 0
    assert np.all(a["V"] < 0)
    assert a["U"][2] < 0 < a["U"][4]


def test_null_embedding_derivatives_independent_of_metric_formula():
    x = np.array([.3, .6, .8, 1.2, 2., 4.])
    h = 1e-5
    lower, upper, a = reference_wall_events(x-h), reference_wall_events(x+h), reference_wall_events(x)
    for coordinate in ("U", "V"):
        finite = (upper[coordinate]-lower[coordinate])/(2*h)*signed_state(x)["Rdot"]
        assert np.allclose(finite, a["ell_d"+coordinate+"_dtau"], rtol=2e-7, atol=1e-7)


def test_v4_slice_anchor_is_preserved():
    curve = initial_slice_curve()
    event = reference_wall_events(np.array([1.5]))
    assert curve["U"][0] == pytest.approx(event["U"][0])
    assert curve["V"][0] == pytest.approx(event["V"][0])
    assert schwarzschild_slice_time(1.) == 0
    assert schwarzschild_slice_time(1.5) < 0
    assert min(abs(x)+abs(t) for x, t in zip(curve["X"], curve["T"])) == 0


def test_early_wall_signal_can_reach_parent_but_late_one_cannot():
    early = reference_wall_events(np.array([.5]))
    p = parent_event_from_U(float(early["U"][0]))
    assert p["V"] > 0 > early["V"][0]
    assert p["U"] == early["U"][0]
    assert p["U"]*p["V"] == pytest.approx(-math.exp(2))
    late = reference_wall_events(np.array([1.5]))
    with pytest.raises(ValueError, match="U < 0"):
        parent_event_from_U(float(late["U"][0]))


@pytest.mark.parametrize("radii", [[.01], [6.01], [float("nan")], [], [1.]*301])
def test_finite_numeric_domain_enforced(radii):
    with pytest.raises(ValueError):
        reference_wall_events(np.array(radii))


def test_record_scope_and_numerical_controls():
    r = causal_reference_record()
    assert r["scope"] == "isolated spherical reference completion only"
    assert not r["classification"]["regular_parent_formation_established"]
    assert r["checks"]["wall_metric_norm_max_error"] < 1e-9
    assert r["checks"]["quadrature_32_64_coordinate_max_difference"] < 1e-10
