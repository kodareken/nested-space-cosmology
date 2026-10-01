"""Independent checks for the weak spherical initial residual. No evolution."""
import copy
import hashlib
import math
import os
from pathlib import Path

import pytest

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np
import sympy as sp

from recursive_horizons.nsc_spherical_cauchy_weak import (
    RECORD_PATH,
    SCHEMA,
    V5_JSON,
    V5_NPZ,
    action_functional,
    action_ry,
    action_total,
    algebraic_total,
    bracket_G,
    coefficient_c,
    coercivity_witness,
    compare_owner_and_action,
    dual_norm_E_parts,
    gauge_constants,
    hessian_quadratic,
    manufactured_controls,
    observe_discrete_inequality,
    owner_total,
    represented_dual_matrix,
    scientific_mismatches,
    verify_saved,
    weak_pairing,
    write_record,
)
from recursive_horizons.nsc_spherical_coupling import _radius_residual
from recursive_horizons.nsc_spherical_feedback_action import lapse_constraint
from recursive_horizons.nsc_spherical_galerkin_coupling import build_grid

_SEALED_SHA256 = "ce5e4345f75b187661f965f96c680a3076331e2df709845c22251edb99b3f148"
_REPLAY = None
_MODULE = Path(__file__).resolve().parents[1] / "src" / "recursive_horizons" / "nsc_spherical_cauchy_weak.py"


def _replay():
    global _REPLAY
    if _REPLAY is None:
        _REPLAY = verify_saved(source_ref="5f10ecd365843d1616e50eb16a20d7acd8377e2c")
    return _REPLAY


def _integrand_J(system, y, G, Q):
    slope = system.derivative @ y
    return float(system.dx * np.sum(2.0 * slope ** 2 + 0.5 * Q ** 2 * y ** 2 + G / (2.0 * y ** 2)))


def test_helper_does_not_own_production_evolution():
    text = _MODULE.read_text()
    assert "run_window" not in text
    assert "solve_initial_radius" not in text
    assert "best_full_source_radius" not in text


def test_sympy_lapse_matches_the_y_reduction():
    area, weyl, gauge, flux, q_value, rho, coordinate = sp.symbols(
        "A C_W C_F Phi Q rho x", real=True, nonzero=True
    )
    y = sp.Function("y")(coordinate)
    radius = y ** 2
    geometric = lapse_constraint(
        0, 0, 0, q_value, radius, 0,
        sp.diff(radius, coordinate), sp.diff(radius, coordinate, 2),
        0, 0, 0, area, weyl, gauge, flux,
    )
    rmag2 = gauge * flux ** 2 / (4 * area)
    G = q_value ** 2 * rmag2 + q_value * rho / (8 * sp.pi * area)
    Ry = -4 * sp.diff(y, coordinate, 2) + q_value ** 2 * y - G / y ** 3
    claimed = -(8 * sp.pi * area / q_value) * y ** 3 * Ry
    assert sp.simplify(sp.together(geometric + rho - claimed)) == 0


def test_owner_matches_algebra_and_smooth_action_not_aliased_action():
    grid = build_grid(32, quadrature=128)
    system = grid.fine
    coordinate = system.xi
    y = 1.7 + 0.08 * np.cos(2.0 * np.pi * coordinate / grid.length)
    y = y + 0.02 * np.cos(4.0 * np.pi * coordinate / grid.length)
    rho = 0.15 + 0.05 * np.cos(2.0 * np.pi * coordinate / grid.length)
    constants = gauge_constants(system)
    G = bracket_G(rho, constants)
    owner = owner_total(system, y ** 2, rho)
    direct = _radius_residual(system, y ** 2, rho)
    algebraic, _G = algebraic_total(system, y ** 2, rho)
    claimed = action_total(y, action_ry(system, y, G, constants), constants)
    assert np.max(np.abs(owner - direct)) == 0.0
    assert np.max(np.abs(owner - algebraic)) < 1e-8
    assert np.max(np.abs(owner - claimed)) < 1e-6
    report = compare_owner_and_action(system, y, rho)
    assert report["defect_explains_action_gap_max"] < 1e-8

    rough = 1.4 + 0.2 * np.cos(2.0 * np.pi * 40 * coordinate / grid.length)
    rough_owner = _radius_residual(system, rough ** 2, rho)
    rough_algebraic, rough_G = algebraic_total(system, rough ** 2, rho)
    rough_claimed = action_total(rough, action_ry(system, rough, rough_G, constants), constants)
    assert np.max(np.abs(rough_owner - rough_algebraic)) < 1e-8
    assert np.max(np.abs(rough_owner - rough_claimed)) > 1.0
    rough_report = compare_owner_and_action(system, rough, rho)
    assert rough_report["scaled_product_rule_defect_max"] > 1.0
    assert rough_report["defect_explains_action_gap_max"] < 1e-8


def test_positive_action_is_strictly_convex_and_hessian_matches():
    grid = build_grid(32, quadrature=128)
    system = grid.fine
    constants = gauge_constants(system)
    coordinate = system.xi
    y = 1.8 + 0.05 * np.cos(2.0 * np.pi * coordinate / grid.length)
    G = 0.2 + 0.04 * np.cos(2.0 * np.pi * coordinate / grid.length)
    other = 2.3 + 0.03 * np.sin(2.0 * np.pi * coordinate / grid.length)
    q_value = constants["Q"]
    mine = _integrand_J(system, y, G, q_value)
    assert abs(mine - action_functional(system, y, G, constants)) < 1e-10
    midpoint = _integrand_J(system, 0.5 * (y + other), G, q_value)
    average = 0.5 * (
        _integrand_J(system, y, G, q_value) + _integrand_J(system, other, G, q_value)
    )
    assert midpoint < average - 1e-6
    direction = 0.03 * np.cos(2.0 * np.pi * coordinate / grid.length)
    step = 1e-5
    second = (
        _integrand_J(system, y + step * direction, G, q_value)
        - 2.0 * _integrand_J(system, y, G, q_value)
        + _integrand_J(system, y - step * direction, G, q_value)
    ) / step ** 2
    hessian = hessian_quadratic(system, y, G, direction, constants)
    hand_hessian = float(system.dx * np.sum(
        4.0 * (system.derivative @ direction) ** 2
        + (q_value ** 2 + 3.0 * G / y ** 4) * direction ** 2
    ))
    assert abs(hand_hessian - hessian) < 1e-12
    assert hessian > 0.0
    assert abs(second - hessian) < 2e-3 * max(1.0, abs(hessian))
    values = coefficient_c(y, G, constants)
    assert np.min(values) > 0.0
    slope = system.derivative @ direction
    lower = min(4.0, float(np.min(values))) * float(
        system.dx * np.sum(slope ** 2 + direction ** 2)
    )
    assert hessian + 1e-12 >= lower
    first = (
        _integrand_J(system, y + step * direction, G, q_value)
        - _integrand_J(system, y - step * direction, G, q_value)
    ) / (2.0 * step)
    assert abs(first - weak_pairing(system, y, G, direction, constants)) < 1e-6


def test_unknown_rho_independence_does_not_establish_mu():
    grid = build_grid(16, quadrature=32)
    system = grid.fine
    constants = gauge_constants(system)
    y = np.full(system.points, 1.5)
    G = np.full(system.points, 0.3)
    unknown = coercivity_witness(system, y, G)
    affirmed = coercivity_witness(system, y, G, rho_independent=True)
    denied = coercivity_witness(system, y, G, rho_independent=False)
    assert unknown["rho_independent_of_y"] is None
    assert unknown["strict_convexity_assumption"] is False
    assert unknown["uniform_mu"] is None
    assert unknown["error_inequality_applicable"] is False
    assert unknown["total_bound"] is None
    assert affirmed["strict_convexity_assumption"] is True
    assert affirmed["uniform_mu"] == min(4.0, constants["Q"] ** 2)
    assert affirmed["uniform_mu_requires_y_star"] is False
    assert affirmed["error_inequality_applicable"] is False
    assert affirmed["total_bound"] is None
    assert affirmed["continuum_positivity_enclosed"] is False
    assert affirmed["rounding_enclosed"] is False
    assert denied["strict_convexity_assumption"] is False
    assert denied["uniform_mu"] is None
    assert denied["error_inequality_applicable"] is False


def test_negative_G_removes_the_coercivity_assumption():
    grid = build_grid(32, quadrature=64)
    system = grid.fine
    constants = gauge_constants(system)
    y = np.ones(system.points)
    G = np.full(system.points, -1.0)
    coefficient = constants["Q"] ** 2 + 3.0 * G / y ** 4
    assert np.max(coefficient) < 0.0
    hessian = hessian_quadratic(system, y, G, np.ones(system.points), constants)
    hand = float(system.dx * np.sum(coefficient))
    assert abs(hessian - hand) < 1e-10
    assert hessian < 0.0
    controls = manufactured_controls()
    assert controls["negative_G"]["strict_convexity_assumption"] is False
    assert controls["negative_G"]["error_inequality_applicable"] is False
    assert controls["negative_hessian_on_constants"] < 0.0
    assert controls["positive_G"]["strict_convexity_assumption"] is True
    assert controls["positive_G"]["error_inequality_applicable"] is False
    inequality = controls["discrete_segment_inequality"]
    assert inequality["holds"] is True
    assert inequality["total_certified"] is False
    assert inequality["is_continuum_certificate"] is False
    assert inequality["upper_bound"] >= inequality["error_norm_E"]


def test_omitted_mode_is_invisible_to_the_represented_residual():
    grid = build_grid(16, quadrature=32)
    system = grid.fine
    coordinate = grid.xi_q
    density = np.cos(2.0 * np.pi * coordinate / grid.length)
    fourier = dual_norm_E_parts(system, density, grid.ng // 2)
    gram = -(system.derivative @ system.derivative) + np.eye(system.points)
    solved = np.linalg.solve(gram, density)
    dense = math.sqrt(float(system.dx * np.dot(density, solved)))
    assert abs(fourier["full"] - dense) < 1e-8
    matrix = represented_dual_matrix(grid, density)
    assert abs(matrix["dual_E"] - fourier["represented"]) < 1e-8

    controls = manufactured_controls()
    omitted = controls["omitted_unresolved_mode"]
    assert omitted["total_certified"] is False
    assert omitted["represented_weak_certifies_total"] is False
    assert omitted["geometry_error_upper_bound"] is None
    assert omitted["tail_bound"] is None
    # G/y^3 of a pure high mode has an O(amplitude^2) mean, so the represented
    # piece is not machine zero. It remains far below the omitted energy.
    assert omitted["represented_piece_is_nonlinear_image"] is True
    assert omitted["represented_weak_dual_E"] < 1e-3
    assert omitted["unresolved_strong_max"] > 1.0
    assert omitted["omitted_energy_norm_E"] > 100.0 * omitted["represented_weak_dual_E"]
    assert omitted["full_strong_max"] > 100.0 * omitted["represented_pull_max"]
    current = controls["constructor_constant_current"]
    assert current["current_unchanged"] is True
    assert current["source_mean_subtracted"] is False
    assert abs(current["residual_mean"] - 0.3) < 1e-10
    assert current["residual_minus_mean_max"] < 1e-9
    assert current["p_Q_max"] < 1e-8


def test_saved_v5_assessment_keeps_partial_terms_and_original_bytes():
    npz_before = V5_NPZ.read_bytes()
    json_before = V5_JSON.read_bytes()
    record_before = RECORD_PATH.read_bytes()
    assert hashlib.sha256(record_before).hexdigest() == _SEALED_SHA256
    report = _replay()
    record = report["sealed"]
    assert V5_NPZ.read_bytes() == npz_before
    assert V5_JSON.read_bytes() == json_before
    assert RECORD_PATH.read_bytes() == record_before
    assert report["wrote"] is False
    assert report["record_bytes_unchanged"] is True
    assert report["v5_bytes_unchanged"] is True
    assert report["head_equality_required"] is False
    assert report["fabricated_hash"] is False
    assert report["historical_checkpoint_head"] == "9a9090a"
    assert report["current_checkpoint_head"] != "9a9090a"
    assert report["limits"]
    assert all(item["role"] == "generator" for item in report["limits"])
    assert all(item["fabricated_hash"] is False for item in report["limits"])
    assert all(item["current_hash"] != item["declared_hash"] for item in report["limits"])
    assert record["v5_bytes_unchanged"] is True
    assert RECORD_PATH.is_file()
    assert record["schema"] == "NSC-SPHERICAL-CAUCHY-WEAK-v1"
    assert record["evolution_performed"] is False
    assert record["total_certified"] is False
    assert record["represented_weak_certifies_total"] is False
    assert record["geometry_error_upper_bound"] is None
    assert record["tolerance_used_as_proof"] is False
    assert record["cpu_budget_exceeded"] is False
    assert record["cpu_seconds"] < 60.0
    assert record["checkpoint_head"] == "9a9090a"
    gaps = record["replay_gaps"]
    assert abs(gaps["nf256_frozen_minus_saved_frozen_operator"]) < 1e-8
    assert abs(gaps["nf256_source_minus_saved_full_C"]) < 1e-8
    assert abs(gaps["nf256_doubled_minus_saved_nq2048"]) < 1e-8
    assert abs(gaps["nf512_source_minus_saved_full_C"]) < 1e-8
    assert abs(gaps["nf256_constructor_current_minus_saved_current_max"]) < 1e-12

    for name, strong_floor in (("nf256", 1e-4), ("nf512", 1e-6)):
        seed = record[name]
        assert seed["evolution_performed"] is False
        assert seed["radius_newton_rerun"] is False
        assert seed["source_columns_unchanged"] is True
        assert seed["ansatz"]["holds"] is True
        assert seed["rho_independent_assumption"] is True
        assert seed["total_certified"] is False
        assert seed["geometry_error_upper_bound"] is None
        source = seed["source_rho_assessment"]
        assert source["total_certified"] is False
        assert source["represented_weak_certifies_total"] is False
        assert source["tail_bound"] is None
        assert source["partial_terms"]["tail_bound"] is None
        represented = source["represented"]["weak_dual_E_matrix"]
        full_strong = source["full_strong"]["owner_max"]
        unresolved = source["unresolved"]["owner_max"]
        roundoff = source["roundoff"]["algebraic_identity_gap_max"]
        for number in (
            represented,
            full_strong,
            unresolved,
            roundoff,
            source["full_strong"]["weak_dual_E"],
            source["product_rule"]["scaled_defect_max"],
            source["unresolved"]["highest_mode_amplitude"],
            seed["constructor_current"]["current_mean"],
            seed["constructor_current"]["held_out_momentum_max"],
        ):
            assert math.isfinite(number)
        assert full_strong > strong_floor
        assert unresolved > 0.5 * full_strong
        assert represented < 1e-6
        assert represented < 1e-3 * full_strong
        assert math.isfinite(source["roundoff"]["unit_roundoff_times_condition"])
        assert roundoff < 1e-6
        assert roundoff * 10.0 < full_strong
        assert source["represented"]["matrix_fourier_gap"] < 1e-12
        quadrature = seed["quadrature"]
        assert quadrature["computed"] is True
        assert quadrature["tail_bound"] is None
        assert quadrature["difference_is_not_a_tail_bound"] is True
        assert math.isfinite(quadrature["difference_of_full_strong_maxima"])
        assert math.isfinite(quadrature["doubled_nq"]["full_strong_owner_max"])
        current = seed["constructor_current"]
        assert current["source_mean_subtracted"] is False
        assert current["radius_newton_rerun"] is False
        assert current["tolerance_is_not_a_proof"] is True
        assert current["current_unchanged"] is True
    nf512 = record["nf512"]["source_rho_assessment"]
    assert nf512["product_rule"]["scaled_defect_max"] > nf512["full_strong"]["owner_max"]
    assert (
        record["nf512"]["quadrature"]["doubled_nq"]["full_strong_owner_max"]
        > nf512["full_strong"]["owner_max"]
    )
    control = record["polynomial_control"]
    assert control["comparison_count"] == 4
    assert control["temporary_script_required"] is False
    assert control["rounding_enclosure"] is None
    assert control["source_rho_modified"] is False
    assert control["current_modified"] is False
    indicator = control["conditioning_indicator"]
    assert abs(indicator["nf512_dense_nq2048"] - 9.981455056262689e-06) < 1e-15
    assert abs(indicator["nf512_dense_nq4096"] - 4.539313647278211e-05) < 1e-12
    assert abs(indicator["nf512_product_nq2048"] - 5.068517579063356e-06) < 1e-12
    assert abs(indicator["nf512_product_nq4096"] - 5.400936559883778e-06) < 1e-12
    for case in control["comparisons"]:
        assert case["alias_guard"]["alias_free"] is True
        assert case["alias_guard"]["degree_P"] == case["nf"] - 1
        assert case["source_rho_modified"] is False
        assert case["rounding_enclosure"] is None
    mu = record["nf512"]["source_rho_assessment"]["coercivity"]["uniform_mu"]
    assert mu == min(4.0, record["nf512"]["source_rho_assessment"]["Q"] ** 2)
    assert record["nf512"]["source_rho_assessment"]["coercivity"]["uniform_mu_requires_y_star"] is False
    assert record["nf512"]["source_rho_assessment"]["coercivity"]["total_bound"] is None
    assert record["nf512"]["geometry_error_upper_bound"] is None
    assert RECORD_PATH.read_bytes() == record_before
    assert hashlib.sha256(V5_JSON.read_bytes()).hexdigest() == record["v5_json_sha256_before"]
    assert hashlib.sha256(V5_NPZ.read_bytes()).hexdigest() == record["v5_npz_sha256_before"]


def test_volatile_execution_metadata_is_separate_from_scientific_settings():
    report = _replay()
    fresh = report["measured"]
    sealed_bytes = RECORD_PATH.read_bytes()
    volatile = copy.deepcopy(report["sealed"])
    volatile["checkpoint_head"] = "e62f804"
    volatile["cpu_seconds"] = 0.25
    volatile["wall_seconds"] = 0.5
    volatile["cpu_budget_exceeded"] = True
    assert scientific_mismatches(volatile, fresh) == []

    independence = copy.deepcopy(report["sealed"])
    independence["nf512"]["rho_independent_assumption"] = False
    assert "nf512.rho_independent_assumption" in scientific_mismatches(independence, fresh)

    mu = copy.deepcopy(report["sealed"])
    mu["nf512"]["source_rho_assessment"]["coercivity"]["uniform_mu"] = 1.0
    assert "nf512.source_rho_assessment.coercivity.uniform_mu" in scientific_mismatches(mu, fresh)

    profile = copy.deepcopy(report["sealed"])
    profile["saved_v5_reference"]["nf512_nq2048_full_C"] = 1.0
    assert "saved_v5_reference.nf512_nq2048_full_C" in scientific_mismatches(profile, fresh)
    assert RECORD_PATH.read_bytes() == sealed_bytes
    assert hashlib.sha256(sealed_bytes).hexdigest() == _SEALED_SHA256


def test_fresh_write_refuses_the_sealed_record(tmp_path):
    sealed_before = RECORD_PATH.read_bytes()
    json_before = V5_JSON.read_bytes()
    npz_before = V5_NPZ.read_bytes()
    payload = {"schema": SCHEMA, "checkpoint_head": "e62f804", "total_certified": False}
    destination = write_record(payload, tmp_path / "fresh-weak.json")
    written = destination.read_text()
    assert '"checkpoint_head": "e62f804"' in written
    assert destination.resolve() != RECORD_PATH.resolve()
    with pytest.raises(FileExistsError):
        write_record(payload, RECORD_PATH)
    assert RECORD_PATH.read_bytes() == sealed_before
    assert V5_JSON.read_bytes() == json_before
    assert V5_NPZ.read_bytes() == npz_before
