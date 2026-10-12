"""Finite observer paths and analytic-domain guards for the fixed wall."""

import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_causal_observers import (
    MAX_RADIUS,
    MIN_RADIUS,
    causal_observers_record,
    centre_emission_to_wall,
    de_sitter_wall_events,
    de_sitter_wall_state,
    require_retained_time,
    wall_to_centre_intercept,
)


@pytest.fixture(scope="module")
def record():
    return causal_observers_record()


def test_regular_boost_survives_both_horizons_and_normal_sign_change():
    radii = np.array([MIN_RADIUS, .5, (4/25)**(1/3), .8, 1., 1.5, MAX_RADIUS])
    data = de_sitter_wall_state(radii)
    assert np.all(np.isfinite(data["gamma"]))
    assert np.all(data["gamma"] >= 1)
    assert np.all(abs(data["v"]) < 1)
    np.testing.assert_allclose(data["gamma"], data["gamma_rationalized"], atol=3e-14)
    np.testing.assert_allclose(data["gamma"]**2-data["gamma_v"]**2, 1, atol=5e-13)
    np.testing.assert_allclose(data["Rdot"], radii*data["gamma"]+data["gamma_v"], atol=1e-13)
    np.testing.assert_allclose(data["beta_minus"], data["gamma"]+radii*data["gamma_v"], atol=1e-13)
    assert data["beta_plus"][1] > 0
    assert data["beta_plus"][3] < 0


def test_unsquared_junction_and_both_bulk_constraints():
    radii = np.linspace(MIN_RADIUS, MAX_RADIUS, 101)
    data = de_sitter_wall_state(radii)
    np.testing.assert_allclose(data["beta_minus"]-data["beta_plus"], 2*radii, atol=1e-14)
    for side in ("minus", "plus"):
        np.testing.assert_allclose(data[f"beta_{side}"]**2,
                                   data["Rdot"]**2+data[f"f_{side}"], atol=1e-12)


def test_anchored_clocks_increase_with_radius_and_follow_input_order():
    radii = np.array([3., .05, 1., 1.5, .8])
    events = de_sitter_wall_events(radii)
    assert events["t_flat"][3] == events["tau_wall"][3] == 0
    ordering = np.argsort(radii)
    assert np.all(np.diff(events["t_flat"][ordering]) > 0)
    assert np.all(np.diff(events["tau_wall"][ordering]) > 0)
    np.testing.assert_allclose(events["r_wall"], -radii*events["eta"])


def test_clock_derivatives_are_proper_time_and_flat_boost():
    centre, step = 1., 1e-5
    event = de_sitter_wall_events([centre-step, centre, centre+step])
    assert (event["tau_wall"][2]-event["tau_wall"][0])/(2*step) == pytest.approx(
        1/event["Rdot"][1], rel=2e-9)
    assert (event["t_flat"][2]-event["t_flat"][0])/(2*step) == pytest.approx(
        event["gamma"][1]/event["Rdot"][1], rel=2e-9)


@pytest.mark.parametrize("ell", [.01, 7., 100.])
def test_units_rescale_clocks_and_coordinates_not_causal_outcomes(ell, record):
    scaled = causal_observers_record(ell)
    assert scaled["access"] == record["access"]
    for name, event in record["named_events"].items():
        assert scaled["named_events"][name] == pytest.approx(event, abs=3e-12)
    for key, values in record["trajectory"]["samples"].items():
        np.testing.assert_allclose(scaled["trajectory"]["samples"][key], values, atol=3e-12)
    base = de_sitter_wall_events([.5, .8, 1.5])
    changed = de_sitter_wall_events(ell*np.array([.5, .8, 1.5]), ell)
    for key in ("t_flat", "tau_wall", "eta", "r_wall"):
        np.testing.assert_allclose(changed[key], ell*base[key], atol=3e-12*ell)


@pytest.mark.parametrize("radius", [0., MIN_RADIUS-1e-8, MAX_RADIUS+1e-8, math.inf, math.nan])
def test_radius_domain_exit_is_explicit(radius):
    with pytest.raises(ValueError, match="radius outside"):
        de_sitter_wall_state(radius)


@pytest.mark.parametrize("ell", [0., -1., math.nan, math.inf])
def test_invalid_scale_is_rejected(ell):
    with pytest.raises(ValueError, match="length scale"):
        de_sitter_wall_events([1.5], ell)


@pytest.mark.parametrize("radii,order", [([], 32), ([[1.5]], 32), ([1.5], 3), ([1.5], 65), ([1.5], 32.5)])
def test_invalid_quadrature_request_is_rejected(radii, order):
    with pytest.raises(ValueError):
        de_sitter_wall_events(radii, order=order)


def test_centre_reception_distinguishes_finite_horizon_and_beyond():
    early = wall_to_centre_intercept(.5, -1.)
    horizon = wall_to_centre_intercept(1., -1.)
    outside = wall_to_centre_intercept(1.5, 0.)
    assert early["status"] == "FINITE_CENTRE_RECEPTION"
    assert early["centre_flat_time"] == pytest.approx(-1+math.log(2))
    assert horizon["status"] == "ONLY_AT_INFINITE_CENTRE_TIME"
    assert horizon["centre_flat_time"] is None
    assert horizon["eta_intercept"] == 0
    assert outside["status"] == "BEYOND_FUTURE_CONFORMAL_BOUNDARY"
    assert outside["centre_flat_time"] is None
    assert outside["eta_intercept"] == .5


def test_historical_centre_emission_is_not_anchor_emission(record):
    event = centre_emission_to_wall(1.5, 0.)
    assert event["centre_flat_time"] == pytest.approx(-math.log(2.5))
    assert event["eta_intercept"] == -2.5
    assert record["named_events"]["C0"]["eta_over_ell"] == -1
    assert record["named_events"]["C_emits_to_W0"]["eta_over_ell"] == -2.5


def test_null_paths_have_correct_orientation_and_remain_inside_cap(record):
    samples = record["trajectory"]["samples"]
    eta, wall_r = np.array(samples["eta_over_ell"]), np.array(samples["r_wall_over_ell"])
    assert np.all(np.diff(eta-wall_r) > 0)
    assert np.all(np.diff(eta+wall_r) > 0)
    for name, receiver in (("C_emits_to_early_wall", "W_early"),
                           ("C_emits_to_W0", "W0"),
                           ("C_escape_threshold", "W_Schwarzschild_horizon")):
        start, end = record["named_events"][name], record["named_events"][receiver]
        start_eta, end_eta = start["eta_over_ell"], end["eta_over_ell"]
        assert end_eta > start_eta
        assert end_eta-start_eta == pytest.approx(end["r_wall_over_ell"])
        mask = (eta > start_eta) & (eta < end_eta)
        assert np.all(eta[mask]-start_eta < wall_r[mask])
    early, reception = record["named_events"]["W_early"], record["named_events"]["C_receives_early_wall"]
    arrival_eta = reception["eta_over_ell"]
    assert arrival_eta > early["eta_over_ell"]
    assert arrival_eta == pytest.approx(early["eta_over_ell"]+early["r_wall_over_ell"])
    mask = (eta > early["eta_over_ell"]) & (eta < arrival_eta)
    assert np.all(arrival_eta-eta[mask] < wall_r[mask])


def test_oldest_positive_control_is_inside_recorded_history(record):
    beginning, end = record["trajectory"]["t_flat_over_ell_domain"]
    earliest_emission = record["named_events"]["C_emits_to_early_wall"]["t_flat_over_ell"]
    require_retained_time(earliest_emission, beginning, end)
    assert earliest_emission-beginning > .019
    with pytest.raises(ValueError, match="outside retained"):
        require_retained_time(beginning-1e-10, beginning, end)
    with pytest.raises(ValueError, match="outside retained"):
        require_retained_time(end+1e-10, beginning, end)
    require_retained_time(beginning, beginning, end)
    require_retained_time(end, beginning, end)


def test_parent_escape_threshold_is_retained_and_strict(record):
    wall = record["named_events"]["W_Schwarzschild_horizon"]
    centre = record["named_events"]["C_escape_threshold"]
    threshold = record["parent_escape_threshold"]
    assert centre["t_flat_over_ell"] == pytest.approx(wall["t_flat_over_ell"]-math.log(1.8))
    assert centre["eta_over_ell"] == pytest.approx(wall["eta_over_ell"]-wall["r_wall_over_ell"])
    assert threshold["matched_path_retained"]
    assert threshold["parent_interpretation_requires_companion_reference"]
    assert not threshold["upper_endpoint_admits_finite_parent_reception"]
    assert record["independent_checks"]["escape_emission_after_retained_start_over_ell"] > .31


def test_future_exclusion_has_analytic_bounds_not_a_grid_inference(record):
    proof = record["future_proof"]
    anchor = 1.5
    assert proof["analytic_not_sample_inference"]
    assert proof["beta_minus_lower_bound"] > 1
    # These constant bounds hold for every x>=anchor, since inverse powers
    # decrease. They certify beta'>0 and q<4x^2 without sampling infinity.
    assert 3/4-2/(5*anchor**3) > 0
    assert 25/16+3/(10*anchor**3)+1/(25*anchor**6) < 4
    assert proof["wall_comoving_radius_lower_bound_over_ell"] > proof["C0_remaining_null_radius_over_ell"]
    assert not record["trajectory"]["finite_plot_is_complete_future"]
    # A different velocity near HR=1 warns against inferring outgoing
    # no-access merely from being at or outside the de Sitter horizon.
    assert de_sitter_wall_state(1.)["v"] < 0
    assert de_sitter_wall_state(1.5)["v"] > 0


def test_independent_residuals_and_finite_budget(record):
    checks = record["independent_checks"]
    for name, value in checks.items():
        if "max_error" in name or "null_error" in name:
            assert value < 1e-11, name
    assert record["admission"]["quadrature_evaluations"] < 20000
    assert record["admission"]["cpu_seconds_max"] == 60
    assert record["theory"]["fixed_within_test"]
    assert not record["admission"]["priority_claim"]
    assert not checks["quadrature_comparison_is_rigorous_error_bound"]
    assert any("field transmission" in line for line in record["limits"])
