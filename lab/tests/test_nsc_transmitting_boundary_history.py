"""Missing action data must never become a passing stationarity residual."""
import numpy as np
import pytest

from recursive_horizons.nsc_transmitting_boundary_history import (
    EndpointVariation, covariance_freedom_probe, match_endpoint_variations,
    compose_full_extended_gate,
)


def test_absent_derivative_stays_unevaluated_even_for_zero_local_gradient():
    local = EndpointVariation("same-history", ((0., 0.),)*4)
    result = match_endpoint_variations(local, None)
    assert result["status"] == "OPEN"
    assert result["maximum_absolute"] is None
    assert result["mismatch_components"] is None


def test_endpoint_pullback_is_required_before_addition():
    local = EndpointVariation("KS-temporal", ((1., -1.),)*4)
    other = EndpointVariation("compact-reflecting", ((-1., 1.),)*4)
    with pytest.raises(ValueError, match="same endpoint basis"):
        match_endpoint_variations(local, other)


def test_covectors_add_action_signs_without_an_extra_normal_flip():
    local = EndpointVariation("control", ((1., -2.),)*4)
    other = EndpointVariation("control", ((-0.5, 1.),)*4)
    result = match_endpoint_variations(local, other)
    assert result["status"] == "FAIL"
    assert result["mismatch_components"]["N"] == [0.5, -1.0]


def test_nonfinite_endpoint_or_tolerance_rejected():
    with pytest.raises(ValueError):
        EndpointVariation("control", ((float("nan"), 0.),)*4)
    with pytest.raises(ValueError):
        match_endpoint_variations(EndpointVariation("control", ((0., 0.),)*4), None, tolerance=float("nan"))


def test_covariance_probe_distinguishes_stabilizer_from_state_change():
    inputs = {"sample_offsets": np.array([0, 2]), "frequency": np.array([1., 2.]),
              "covariance_seed": np.array([np.eye(2)*0.2, np.eye(2)*0.8])}
    probe = covariance_freedom_probe(inputs)
    assert probe["maximum_isometry_tangent_residual"] == 0
    assert probe["maximum_state_tangent"] == pytest.approx(0.6)
    inputs["covariance_seed"][1] = inputs["covariance_seed"][0]
    assert covariance_freedom_probe(inputs)["maximum_state_tangent"] == 0


def test_endpoint_numeric_match_cannot_by_itself_prove_existence():
    base = {"components_expose_residuals": True, "killed_homogeneous_regression_pass": True}
    result = compose_full_extended_gate(base, {"status": "PASS", "selection_residual": {}}, {"match": {"status": "PASS"}})
    assert result["decision"] == "OPEN"
    assert result["blocking_residuals"] == ["assembled history stationarity residuals for g_star,C_star"]
    assert not result["extended_existence_claimed"]
    assert not result["metric_timestep_started"]
