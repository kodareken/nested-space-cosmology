"""Conditional initial-lapse bound. Does not step a state or rewrite sealed bytes."""
import copy
import json
import math

import numpy as np
import pytest

from recursive_horizons.nsc_spherical_cauchy_weak import (
    ERROR_RECORD_PATH,
    ERROR_SCHEMA,
    RECORD_PATH,
    V5_JSON,
    V5_NPZ,
    action_ry,
    bandlimited_matter_density,
    dual_norm_E_parts,
    energy_norm,
    fft_periodic_derivative,
    gauge_constants,
    hessian_quadratic,
    initial_lapse_binding_report,
    inverse_square_log_divergence,
    lapse_existence_theorem,
    polynomial_degree_caps,
    geometric_constraint_transport,
    matter_shift_to_g,
    nearby_exact_initial_state,
    _nearby_from_bounds,
    regeneration_window_coverage,
    source_constraint_normalization,
    saved_lapse_existence,
    saved_shift_constraint_bound,
    write_initial_lapse_record,
)
from recursive_horizons.nsc_spherical_feedback_action import feedback_Z, lapse_constraint
from recursive_horizons.nsc_spherical_coupling import source_from_columns
from recursive_horizons.nsc_spherical_galerkin_coupling import (
    blank_state,
    build_grid,
    load_physical_columns,
    prolong_state,
)

SEALED_SOURCE_REF = "5f10ecd365843d1616e50eb16a20d7acd8377e2c"


def test_degree_caps_follow_the_declared_bands():
    caps = polynomial_degree_caps(512, 511)
    assert caps["degree_r"] == 255
    assert caps["degree_radius_products"] == 510
    assert caps["degree_source"] == 511
    assert caps["degree_P"] == 511
    assert caps["highest_antiperiodic_mode"] == 255.5
    with pytest.raises(ValueError):
        polynomial_degree_caps(512, 512)


def test_highest_antiperiodic_product_stops_at_degree_nf_minus_one():
    fermions = 32
    length = 8.0
    quadrature = 128
    nodes = np.arange(quadrature) * (length / quadrature)
    plus = fermions / 2 - 0.5
    minus = -plus
    product = np.conjugate(np.exp(2j * np.pi * plus * nodes / length)) * np.exp(
        2j * np.pi * minus * nodes / length
    )
    coefficients = np.fft.fft(product) / quadrature
    modes = np.fft.fftfreq(quadrature) * quadrature
    amplitude = np.abs(coefficients)
    occupied = amplitude > 1e-10
    assert int(np.max(np.abs(modes[occupied]))) == fermions - 1
    assert float(np.sum(amplitude[np.abs(modes) > fermions - 1] ** 2)) < 1e-20


def test_bandlimited_density_matches_the_owner_source():
    grid = build_grid(32, quadrature=128)
    phi0, phi1, _preparation, _problems = load_physical_columns(32)
    state = blank_state(grid, phi0, phi1)
    fine = prolong_state(grid, state)
    owner = source_from_columns(grid.fine, fine)
    rho_owner = np.asarray(owner["force_L"] / grid.dx_q, dtype=float)
    rho_band = bandlimited_matter_density(phi0, phi1, grid.length, grid.nq, float(state.Q[0]))
    assert float(np.max(np.abs(rho_band - rho_owner))) < 1e-8
    coefficients = np.fft.fft(rho_band) / grid.nq
    modes = np.fft.fftfreq(grid.nq) * grid.nq
    assert float(np.sum(np.abs(coefficients[np.abs(modes) > grid.nf - 1]) ** 2)) < 1e-16


def test_smooth_numerator_matches_the_y_residual_and_the_convex_chain():
    grid = build_grid(32, quadrature=128)
    system = grid.fine
    constants = gauge_constants(system)
    coordinate = system.xi
    length = grid.length
    amplitude = 1e-4
    y_star = np.full(system.points, 2.0)
    capital_g = np.full(system.points, constants["Q"] ** 2 * y_star ** 4)
    y = y_star + amplitude * np.cos(2.0 * math.pi * coordinate / length)
    residual = action_ry(system, y, capital_g, constants)
    radius = y ** 2
    slope = fft_periodic_derivative(radius, length, 1)
    second = fft_periodic_derivative(radius, length, 2)
    numerator = 2.0 * radius * second - slope ** 2 - constants["Q"] ** 2 * radius ** 2 + capital_g
    assert np.max(np.abs(numerator + y ** 3 * residual)) < 1e-8
    dual = dual_norm_E_parts(system, residual, system.points)["full"]
    residual_l2 = math.sqrt(float(system.dx * np.sum(residual ** 2)))
    assert dual <= residual_l2 * (1.0 + 1e-8)
    mu = min(4.0, constants["Q"] ** 2)
    error = energy_norm(system, y - y_star)
    assert error <= residual_l2 / mu * (1.0 + 1e-6)
    test = np.cos(2.0 * math.pi * coordinate / length)
    hessian = hessian_quadratic(system, y, capital_g, test, constants)
    energy = energy_norm(system, test) ** 2
    assert hessian >= mu * energy * (1.0 - 1e-8)


def test_barrier_log_divergence_and_strict_convexity():
    # Equality case: y'=1 on an arc of length x, so ∫(y')^2 = x and y(x)^2 = x^2.
    arc = 0.5
    derivative_l2_squared = arc * 1.0 ** 2
    assert arc ** 2 == derivative_l2_squared * arc
    lower = inverse_square_log_divergence(1.0, 1e-9, 0.25)
    assert lower > math.log(1e8)
    value = 2.0
    step = 1e-4
    second = ((value + step) ** -2 - 2.0 * value ** -2 + (value - step) ** -2) / step ** 2
    assert abs(second - 6.0 * value ** -4) < 1e-6
    left, right, weight = 1.0, 3.0, 0.3
    midpoint = (1.0 - weight) * left + weight * right
    assert midpoint ** -2 < (1.0 - weight) * left ** -2 + weight * right ** -2
    refused = lapse_existence_theorem(length=8.0, q_squared=0.06, g_lower=0.0)
    assert refused["minimizer_exists"] is False
    proved = lapse_existence_theorem(length=8.0, q_squared=0.06, g_lower=0.05)
    assert proved["unique_positive_critical_point"] is True
    assert proved["smooth_critical_point"] is True
    assert proved["evolution_error_bound"] is None
    assert any("Fatou" in step for step in proved["steps"])
    assert any("Cauchy-Schwarz" in step for step in proved["steps"])


def test_saved_existence_and_shift_correction_keep_sealed_bytes():
    error_before = ERROR_RECORD_PATH.read_bytes()
    weak_before = RECORD_PATH.read_bytes()
    v5_before = V5_NPZ.read_bytes()
    v5_json_before = V5_JSON.read_bytes()
    record = json.loads(error_before)
    assert record["schema"] == ERROR_SCHEMA
    assert record["minimizer_existence_proved"] is None
    assert record["evolution_error_bound"] is None
    assert record["geometry_error_upper_bound"] < 5e-6
    existence = saved_lapse_existence()
    assert existence["hypotheses_hold"] is True
    assert existence["unique_positive_critical_point"] is True
    assert existence["smooth_critical_point"] is True
    assert existence["historical_field_left_null"] is True
    assert existence["recomputed_numerator"] is False
    assert existence["g_lower"] == record["G_lower"]
    assert existence["evolution_error_bound"] is None
    shift = saved_shift_constraint_bound()
    assert shift["cpu_seconds"] < 30.0
    assert shift["symmetry_defect_max"] == 0.0
    assert shift["mean_deleted"] is False
    assert shift["source_modified"] is False
    assert shift["mean_exactly_zero"] is False
    assert shift["mean_contains_zero"] is False
    assert shift["mean_j_lower"] > 0.0
    assert shift["current_infinity_upper"] < 1e-12
    assert shift["resolves_initial_shift_against_areal_bar"] is True
    assert shift["recomputed_lapse_numerator"] is False
    assert shift["evolution_error_bound"] is None
    correction = shift["correction"]
    assert correction["admissible"] is True
    assert correction["source_modified"] is False
    assert correction["mean_deleted"] is False
    assert correction["continuum_shift_residual_after_correction"] == 0.0
    assert correction["delta_C_infinity_upper"] < 1e-24
    assert abs(correction["lambda_upper"]) < 1e-12
    assert correction["critical_point_mean_leftover_upper"] < 1e-18
    area = record["settings"]["A"]
    gauge = record["settings"]["Q"]
    zero = lapse_constraint(0, 0, 0, gauge, 2.0, 0, 0.1, -0.2, 0, 0, 0, area, -0.001, 0.01, 1.0)
    moved = lapse_constraint(0, 0.02, 0, gauge, 2.0, 0, 0.1, -0.2, 0, 0, 0, area, -0.001, 0.01, 1.0)
    assert abs((moved - zero) - (0.02 ** 2) / (4 * feedback_Z(area) * gauge)) < 1e-12
    binding = initial_lapse_binding_report(record, source_ref=SEALED_SOURCE_REF)
    assert binding["ok"] is True
    assert binding["effect_comparison_ok"] is True
    assert binding["generator_limits"]
    assert binding["generator_limits"][0]["role"] == "generator"
    mutated = copy.deepcopy(record)
    mutated["settings"]["occupations"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    assert initial_lapse_binding_report(mutated, source_ref=SEALED_SOURCE_REF)["ok"] is False
    source = copy.deepcopy(record)
    source["inputs"]["source_hashes"]["src/recursive_horizons/nsc_spherical_coupling.py"] = "0" * 64
    assert initial_lapse_binding_report(source, source_ref=SEALED_SOURCE_REF)["ok"] is False
    episode = copy.deepcopy(record)
    episode["effect_comparison"]["sha256"] = "0" * 64
    episode_report = initial_lapse_binding_report(episode, source_ref=SEALED_SOURCE_REF)
    assert episode_report["ok"] is True
    assert episode_report["effect_comparison_ok"] is False
    with pytest.raises(FileExistsError):
        write_initial_lapse_record(record, RECORD_PATH)
    with pytest.raises(FileExistsError):
        write_initial_lapse_record(record, V5_JSON)
    with pytest.raises(FileExistsError):
        write_initial_lapse_record(record, ERROR_RECORD_PATH)
    assert ERROR_RECORD_PATH.read_bytes() == error_before
    assert RECORD_PATH.read_bytes() == weak_before
    assert V5_NPZ.read_bytes() == v5_before
    assert V5_JSON.read_bytes() == v5_json_before


def test_sealed_initial_certificate_requires_explicit_historical_source_domain():
    record = json.loads(ERROR_RECORD_PATH.read_text())
    current = initial_lapse_binding_report(record)
    assert current["ok"] is False
    assert "nsc_spherical_galerkin_coupling.py" in current["reason"]
    historical = initial_lapse_binding_report(record, source_ref=SEALED_SOURCE_REF)
    assert historical["ok"] is True
    assert set(historical["source_origins"].values()) == {SEALED_SOURCE_REF}
    missing = initial_lapse_binding_report(record, source_ref="0" * 40)
    assert missing["ok"] is False


def test_lapse_increment_is_the_stated_addition_to_g():
    area = 0.04501936182826115
    gauge = 0.25132493326913924
    zee = feedback_Z(area)
    momentum = 0.02
    base = lapse_constraint(
        0, 0, 0, gauge, 2.0, 0, 0.1, -0.2, 0, 0, 0, area, -0.001, 0.01, 1.0,
    )
    moved = lapse_constraint(
        0, momentum, 0, gauge, 2.0, 0, 0.1, -0.2, 0, 0, 0, area, -0.001, 0.01, 1.0,
    )
    delta = moved - base
    assert abs(delta - momentum ** 2 / (4.0 * zee * gauge)) < 1e-12
    added = matter_shift_to_g(momentum ** 2, area, zee)
    assert abs(delta - (8.0 * math.pi * area / gauge) * added) < 1e-12
    assert zee < 0.0 and added < 0.0


def test_nearby_state_refuses_an_uncontrolled_critical_point_minimum():
    sealed = {
        "G_lower": 0.05, "y_min": 0.1, "y_max": 2.0,
        "pointwise_radius_error_upper_bound": 1.0,
        "geometry_error_upper_bound": 1.0,
        "settings": {"Q": 0.25, "A": 0.045},
    }
    shift = {"mean_j_lower": 1e-17, "mean_j_upper": 2e-17}
    with pytest.raises(ValueError, match="critical-point minimum"):
        _nearby_from_bounds(shift, sealed, 0.2, 0.3, 0.5, 8.0, 2e-15)


def test_nearby_initial_state_and_transport_do_not_reset_records():
    error_before = ERROR_RECORD_PATH.read_bytes()
    weak_before = RECORD_PATH.read_bytes()
    v5_before = V5_NPZ.read_bytes()
    nearby = nearby_exact_initial_state()
    assert nearby["status"] == "EXACT_NEARBY_INITIAL_STATE"
    assert nearby["not_actual_reset"] is True
    assert nearby["actual_saved_pr"] == 0.0
    assert nearby["correction_installed_in_records"] is False
    assert nearby["source_unchanged"] is True
    assert nearby["D0_lower"] > 0.0
    assert nearby["H_plus_lower"] > 0.0
    assert nearby["H_minus_upper"] < 0.0
    assert nearby["constraints_C_and_D_exact"] is True
    assert nearby["g_lambda_lower"] > 0.0
    assert nearby["critical_y_lower"] < json.loads(error_before)["y_min"]
    assert nearby["critical_y_upper"] > json.loads(error_before)["y_max"]
    assert "outward interval" in nearby["arithmetic"]
    assert nearby["evolution_error_bound"] is None
    transport = geometric_constraint_transport()
    assert transport["geometric_identity"] is True
    assert transport["momentum_equations_match_euler"] is True
    assert transport["l2_energy_of_C_and_D_closes"] is False
    assert transport["l2_energy_of_C_and_D_over_Q_closes"] is False
    assert transport["D_t_contains_C_x"] is False
    assert transport["includes_dirac_force"] is False
    assert transport["controls_physical_observables"] is False
    assert "F_Qx" in transport["force_Q_extra_D"]
    radius, radial, spacing = 2.0, 0.3, 0.01
    force_l, force_beta = 0.4, -0.2
    mapped = source_constraint_normalization(force_l, force_beta, radius, radial, spacing)
    sphere = 4.0 * math.pi * radius ** 4
    coordinate_l = force_l / spacing
    coordinate_beta = force_beta / spacing
    assert abs(mapped["physical_density"] - coordinate_l / (sphere * radial)) < 1e-12
    assert abs(mapped["physical_current"] - (-coordinate_beta) / (sphere * radial ** 2)) < 1e-12
    assert abs(mapped["normal_energy_sample"] - force_l / radius) < 1e-12
    assert abs(mapped["current_volume_sample"] - abs(force_beta) / (radius * radial)) < 1e-12
    window = regeneration_window_coverage()
    assert window["missing_observable_primitive"] == transport["missing_observable_primitive"]
    assert window["state_error_bound"] is None
    assert window["constraint_norm_is_state_accuracy"] is False
    assert "one_percent" not in json.dumps(window)
    assert len(window["cases"]) == 4
    for case in window["cases"]:
        assert case["time_start"] == 0.05 and case["time_end"] == 0.085
        assert case["is_trajectory_error"] is False
        assert case["state_error_bound"] is None
        assert case["full_hamilton_max"] > 0.0
        assert case["normal_energy_residual_majorant"] > 0.0
        assert case["smearing_over_field_exchange"] is not None
        assert case["ward_max"] < 1e-15
        assert case["source_gap_max"] == 0.0
    assert ERROR_RECORD_PATH.read_bytes() == error_before
    assert RECORD_PATH.read_bytes() == weak_before
    assert V5_NPZ.read_bytes() == v5_before
