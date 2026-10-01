"""Independent checks for the weak initial-geometry residual.

The identities are recomputed from the lapse constraint and from the
first variation of

    J(y) = ∫ [2 (y')^2 + Q^2 y^2 / 2 + G / (2 y^2)] dx,

with r = y^2, y > 0, and

    G = Q^2 r_mag^2 + Q rho / (8 pi A),
    Ry = -4 y'' + Q^2 y - G / y^3,
    C = -(8 pi A / Q) y^3 Ry.

A green result is a discrete identity or a conditional indicator. It is
not a continuum enclosure and it does not replay the v5 diagnostic.
Declared-input binding is checked separately and does not remeasure.
"""
from __future__ import annotations

import copy
import json
import math
import os
from functools import lru_cache

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

from recursive_horizons.nsc_spherical_cauchy_data import shift_momentum, shift_residual
from recursive_horizons.nsc_spherical_cauchy_weak import (
    RECORD_PATH,
    action_functional,
    action_ry,
    action_total,
    assess_radius,
    binding_report,
    bracket_G,
    coefficient_c,
    coercivity_witness,
    constructor_current,
    energy_norm,
    gauge_constants,
    hessian_quadratic,
    observe_discrete_inequality,
    owner_total,
    represented_dual_matrix,
    weak_pairing,
)
from recursive_horizons.nsc_spherical_coupling import (
    CauchyState,
    hamilton_constraint,
    magnetic_radius_square,
    source_from_columns,
)
from recursive_horizons.nsc_spherical_feedback_action import lapse_constraint
from recursive_horizons.nsc_spherical_galerkin_coupling import (
    blank_state,
    build_grid,
    manufactured_columns,
    prolong_state,
)


@lru_cache(maxsize=1)
def _grid():
    return build_grid(16, quadrature=64)


def _probe(system, radius):
    points = int(system.points)
    gauge = float(system.calibration["b0"] / system.calibration["a0"])
    return CauchyState(
        Q=np.full(points, gauge),
        r=np.asarray(radius, dtype=float),
        chi=np.zeros(points),
        p_Q=np.zeros(points),
        p_r=np.zeros(points),
        p_chi=np.zeros(points),
        phi0=np.zeros((points, 1), dtype=np.complex128),
        phi1=np.zeros((points, 1), dtype=np.complex128),
    )


def _original_constraint(system, radius, rho):
    return hamilton_constraint(system, _probe(system, radius)) + np.asarray(rho, dtype=float)


def _hand_G(rho, constants):
    rho = np.asarray(rho, dtype=float)
    return constants["Q"] ** 2 * constants["rmag2"] + (
        constants["Q"] * rho / constants["eight_pi_A"]
    )


def _hand_algebraic(system, radius, rho, constants):
    """Owner lapse residual with D(F), before the y product rule."""
    radius = np.asarray(radius, dtype=float)
    derivative = system.derivative
    radius_derivative = derivative @ radius
    radius_square_second = derivative @ (derivative @ (radius * radius))
    G = _hand_G(rho, constants)
    total = (constants["eight_pi_A"] / constants["Q"]) * (
        radius_square_second
        - 3.0 * radius_derivative ** 2
        - constants["Q"] ** 2 * radius ** 2
        + G
    )
    return total, G


def _hand_action(system, y, G, constants):
    """C = -(8 pi A / Q) y^3 Ry after r = y^2."""
    y = np.asarray(y, dtype=float)
    second = system.derivative @ (system.derivative @ y)
    Ry = -4.0 * second + constants["Q"] ** 2 * y - np.asarray(G, dtype=float) / y ** 3
    total = -(constants["eight_pi_A"] / constants["Q"]) * y ** 3 * Ry
    return total, Ry


def _scaled_product_defect(system, y, constants):
    y = np.asarray(y, dtype=float)
    derivative = system.derivative
    radius = y * y
    radius_derivative = derivative @ radius
    radius_square_second = derivative @ (derivative @ (radius * radius))
    y_second = derivative @ (derivative @ y)
    defect = radius_square_second - 3.0 * radius_derivative ** 2 - 4.0 * y ** 3 * y_second
    return (constants["eight_pi_A"] / constants["Q"]) * defect


def _roundoff_gap(left, right):
    difference = np.max(np.abs(np.asarray(left) - np.asarray(right)))
    scale = max(1.0, float(np.max(np.abs(left))), float(np.max(np.abs(right))))
    return difference, difference / scale


def _mode_energy(amplitude, mode, length):
    """Exact periodic energy norm of amplitude * sin(2 pi mode x / L)."""
    if mode == 0:
        return abs(amplitude) * math.sqrt(length)
    wavenumber = 2.0 * math.pi * mode / length
    return abs(amplitude) * math.sqrt((wavenumber ** 2 + 1.0) * length / 2.0)


def test_sqrt_reduction_matches_the_original_constraint_and_keeps_the_source():
    area, weyl, gauge, flux, q_symbol, rho_symbol, coordinate = sp.symbols(
        "A C_W C_F Phi Q rho x", real=True, nonzero=True
    )
    y_symbol = sp.Function("y")(coordinate)
    radius_symbol = y_symbol ** 2
    geometric = lapse_constraint(
        0, 0, 0, q_symbol, radius_symbol, 0,
        sp.diff(radius_symbol, coordinate), sp.diff(radius_symbol, coordinate, 2),
        0, 0, 0, area, weyl, gauge, flux,
    )
    rmag2 = gauge * flux ** 2 / (4 * area)
    G_symbol = q_symbol ** 2 * rmag2 + q_symbol * rho_symbol / (8 * sp.pi * area)
    Ry_symbol = (
        -4 * sp.diff(y_symbol, coordinate, 2)
        + q_symbol ** 2 * y_symbol
        - G_symbol / y_symbol ** 3
    )
    claimed = -(8 * sp.pi * area / q_symbol) * y_symbol ** 3 * Ry_symbol
    flipped = -claimed
    assert sp.simplify(sp.together(geometric + rho_symbol - claimed)) == 0
    assert sp.simplify(sp.together(geometric + rho_symbol - flipped)) != 0

    grid = _grid()
    system = grid.fine
    constants = gauge_constants(system)
    assert abs(
        constants["rmag2"] - magnetic_radius_square(system.coefficients)
    ) <= 1e-15 * max(1.0, abs(constants["rmag2"]))
    assert abs(constants["eight_pi_A"] - 8.0 * math.pi * constants["A"]) == 0.0
    coordinate = system.xi
    smooth = (
        1.7
        + 0.08 * np.cos(2.0 * math.pi * coordinate / grid.length)
        + 0.02 * np.cos(4.0 * math.pi * coordinate / grid.length)
    )
    rho = 0.15 + 0.05 * np.cos(2.0 * math.pi * coordinate / grid.length)
    original = _original_constraint(system, smooth ** 2, rho)
    algebraic, G = _hand_algebraic(system, smooth ** 2, rho, constants)
    action, _Ry = _hand_action(system, smooth, G, constants)
    helper_owner = owner_total(system, smooth ** 2, rho)
    helper_action = action_total(smooth, action_ry(system, smooth, G, constants), constants)
    for candidate in (algebraic, action, helper_owner, helper_action):
        _difference, relative = _roundoff_gap(original, candidate)
        assert relative < 1e-9
    defect = _scaled_product_defect(system, smooth, constants)
    explained, explained_relative = _roundoff_gap(original - action, defect)
    assert explained_relative < 1e-9

    rough_mode = 20
    rough = 1.4 + 0.2 * np.cos(2.0 * math.pi * rough_mode * coordinate / grid.length)
    rough_rho = np.full(system.points, 0.2)
    rough_original = _original_constraint(system, rough ** 2, rough_rho)
    rough_G = _hand_G(rough_rho, constants)
    rough_action, _rough_Ry = _hand_action(system, rough, rough_G, constants)
    rough_gap = float(np.max(np.abs(rough_original - rough_action)))
    rough_defect = _scaled_product_defect(system, rough, constants)
    defect_gap, _defect_relative = _roundoff_gap(rough_original - rough_action, rough_defect)
    assert rough_gap > 1.0
    assert defect_gap < 1e-9 * max(1.0, rough_gap)

    rho_constant = np.full(system.points, 0.4)
    vacuum_radius = np.full(system.points, math.sqrt(constants["rmag2"]))
    vacuum_gap, _vacuum_relative = _roundoff_gap(
        _original_constraint(system, vacuum_radius, np.zeros(system.points)),
        np.zeros(system.points),
    )
    assert vacuum_gap < 1e-9
    kept = _original_constraint(system, vacuum_radius, rho_constant)
    kept_gap, _kept_relative = _roundoff_gap(kept, rho_constant)
    assert kept_gap < 1e-9
    assert float(np.max(np.abs(kept))) > 0.3
    G_constant = _hand_G(rho_constant, constants)
    critical_y = (G_constant / constants["Q"] ** 2) ** 0.25
    balanced = _original_constraint(system, np.full(system.points, critical_y ** 2), rho_constant)
    balanced_gap, _balanced_relative = _roundoff_gap(balanced, np.zeros(system.points))
    assert balanced_gap < 1e-9
    one_step = _original_constraint(system, np.ones(system.points), rho_constant)
    assert float(np.max(np.abs(one_step))) > 0.3


def test_coercive_norm_is_an_indicator_until_the_segment_is_controlled():
    grid = _grid()
    system = grid.fine
    constants = gauge_constants(system)
    coordinate = system.xi
    y = 1.8 + 0.05 * np.cos(2.0 * math.pi * coordinate / grid.length)
    other = 2.2 + 0.04 * np.sin(2.0 * math.pi * coordinate / grid.length)
    G = 0.2 + 0.03 * np.cos(2.0 * math.pi * coordinate / grid.length)
    direction = 0.03 * np.cos(2.0 * math.pi * coordinate / grid.length)

    def integrand(profile):
        slope = system.derivative @ profile
        return float(system.dx * np.sum(
            2.0 * slope ** 2
            + 0.5 * constants["Q"] ** 2 * profile ** 2
            + G / (2.0 * profile ** 2)
        ))

    assert abs(integrand(y) - action_functional(system, y, G, constants)) < 1e-12
    average = 0.5 * (integrand(y) + integrand(other))
    assert integrand(0.5 * (y + other)) < average
    step = 1e-3
    first = (integrand(y + step * direction) - integrand(y - step * direction)) / (2.0 * step)
    second = (
        integrand(y + step * direction)
        - 2.0 * integrand(y)
        + integrand(y - step * direction)
    ) / step ** 2
    pairing = weak_pairing(system, y, G, direction, constants)
    hessian = hessian_quadratic(system, y, G, direction, constants)
    assert abs(first - pairing) < 1e-6 * max(1.0, abs(pairing))
    assert abs(second - hessian) < 1e-5 * max(1.0, abs(hessian))
    c_min = float(np.min(coefficient_c(y, G, constants)))
    assert c_min > 0.0
    slope = system.derivative @ direction
    energy_square = float(system.dx * np.sum(slope ** 2 + direction ** 2))
    assert hessian >= min(4.0, c_min) * energy_square - 1e-12

    negative_G = np.full(system.points, -(constants["Q"] ** 2 + 1.0))
    constant_direction = np.ones(system.points)
    negative_hessian = hessian_quadratic(
        system, np.ones(system.points), negative_G, constant_direction, constants
    )
    assert negative_hessian < 0.0
    high = np.cos(2.0 * math.pi * 8.0 * coordinate / grid.length)
    assert hessian_quadratic(system, np.ones(system.points), negative_G, high, constants) > 0.0
    refused = coercivity_witness(system, np.ones(system.points), negative_G, rho_independent=True)
    assert refused["positive_G"] is False
    assert refused["strict_convexity_assumption"] is False
    assert refused["error_inequality_applicable"] is False
    blocked = observe_discrete_inequality(
        system, np.ones(system.points), np.full(system.points, 1.1), negative_G, -1.0,
    )
    assert blocked["holds"] is False
    assert blocked["upper_bound"] is None
    assert blocked["total_certified"] is False
    assert blocked["is_continuum_certificate"] is False

    small_negative = np.full(system.points, -1e-4)
    sampled = coercivity_witness(
        system, np.full(system.points, 1.5), small_negative, rho_independent=True,
    )
    assert sampled["c_min"] > 0.0
    assert sampled["positive_G"] is False
    assert sampled["strict_convexity_assumption"] is False
    assert sampled["error_inequality_applicable"] is False

    aliased = 0.2 + np.cos(2.0 * math.pi * system.points * coordinate / grid.length)
    aliased_witness = coercivity_witness(
        system, np.full(system.points, 1.2), aliased, rho_independent=True,
    )
    analytic_c_min = constants["Q"] ** 2 + 3.0 * (0.2 - 1.0) / 1.2 ** 4
    assert float(np.min(aliased)) > 1.0
    assert analytic_c_min < 0.0
    assert aliased_witness["positive_G"] is True
    assert aliased_witness["frozen_alpha"] > 0.0
    assert aliased_witness["frozen_alpha"] > analytic_c_min
    assert aliased_witness["error_inequality_applicable"] is False

    y_star = np.full(system.points, 2.0)
    G_star = np.full(system.points, constants["Q"] ** 2 * y_star ** 4)
    displacement = 1e-3
    nearby = y_star + displacement
    c_seg = constants["Q"] ** 2 + 3.0 * float(G_star[0]) / float(np.max(nearby)) ** 4
    sharp = observe_discrete_inequality(system, nearby, y_star, G_star, c_seg)
    assert sharp["holds"] is True
    assert sharp["upper_bound"] >= sharp["error_norm_E"]
    assert sharp["upper_bound"] < 1.01 * sharp["error_norm_E"]
    assert sharp["total_certified"] is False
    assert sharp["is_continuum_certificate"] is False
    overstated = observe_discrete_inequality(system, nearby, y_star, G_star, 2.0 * c_seg)
    assert overstated["holds"] is False
    assert overstated["total_certified"] is False
    assert overstated["upper_bound"] < sharp["error_norm_E"]


def test_unknown_source_independence_does_not_establish_convexity():
    grid = _grid()
    system = grid.fine
    y = np.full(system.points, 1.5)
    G = np.full(system.points, 0.3)
    unknown = coercivity_witness(system, y, G)
    affirmed = coercivity_witness(system, y, G, rho_independent=True)
    denied = coercivity_witness(system, y, G, rho_independent=False)
    assert unknown["rho_independent_of_y"] is None
    assert unknown["error_inequality_applicable"] is False
    assert unknown["strict_convexity_assumption"] is False, (
        "omitted rho_independent must not set the convexity hypothesis; "
        f"rho_independent_of_y={unknown['rho_independent_of_y']!r}, "
        f"strict_convexity_assumption={unknown['strict_convexity_assumption']!r}"
    )
    assert affirmed["strict_convexity_assumption"] is True
    assert affirmed["error_inequality_applicable"] is False
    assert denied["strict_convexity_assumption"] is False
    assert denied["error_inequality_applicable"] is False


def test_hidden_high_frequency_residual_is_not_certified_zero():
    grid = _grid()
    system = grid.fine
    constants = gauge_constants(system)
    coordinate = system.xi
    mode = int(grid.ng // 2 + 5)
    assert mode < grid.nq // 2
    amplitude = 0.02
    y_star = np.full(system.points, 2.0)
    G_star = np.full(system.points, constants["Q"] ** 2 * y_star ** 4)
    omitted = amplitude * np.cos(2.0 * math.pi * mode * coordinate / grid.length)
    y = y_star + omitted
    c_seg = constants["Q"] ** 2 + 3.0 * float(G_star[0]) / float(np.max(y)) ** 4
    alpha = min(4.0, c_seg)
    Ry = action_ry(system, y, G_star, constants)
    represented = represented_dual_matrix(grid, Ry)
    assert represented["dual_E"] is not None
    false_upper = represented["dual_E"] / alpha
    error = energy_norm(system, y - y_star)
    assert false_upper < error / 100.0
    discrete = observe_discrete_inequality(system, y, y_star, G_star, c_seg)
    assert discrete["holds"] is True
    assert discrete["upper_bound"] >= discrete["error_norm_E"]
    assert discrete["upper_bound"] > 100.0 * false_upper
    assert discrete["total_certified"] is False
    assert discrete["is_continuum_certificate"] is False
    inflated = observe_discrete_inequality(system, y, y_star, G_star, 10.0 * c_seg)
    assert inflated["holds"] is True
    assert inflated["total_certified"] is False
    assert inflated["is_continuum_certificate"] is False

    source_mode = 0.25 * np.cos(2.0 * math.pi * mode * coordinate / grid.length)
    radius = np.ones(grid.ng)
    assessed = assess_radius(grid, radius, source_mode, rho_independent=True)
    original = _original_constraint(system, np.ones(system.points), source_mode)
    source_gap, _source_relative = _roundoff_gap(original, source_mode)
    assert source_gap < 1e-9
    assert assessed["represented"]["owner_pull_max"] < 1e-6 * assessed["full_strong"]["owner_max"]
    assert assessed["unresolved"]["owner_max"] > 0.9 * assessed["full_strong"]["owner_max"]
    assert assessed["represented"]["weak_dual_E_matrix"] < 1e-6 * assessed["full_strong"]["owner_max"]
    assert assessed["total_certified"] is False
    assert assessed["represented_weak_certifies_total"] is False
    assert assessed["geometry_error_upper_bound"] is None
    assert assessed["tail_bound"] is None
    assert assessed["partial_terms"]["tail_bound"] is None
    assert assessed["coercivity"]["error_inequality_applicable"] is False

    hidden_mode = int(system.points)
    hidden_amplitude = 0.05
    continuum = _mode_energy(hidden_amplitude, hidden_mode, grid.length)
    aligned = hidden_amplitude * np.sin(2.0 * math.pi * hidden_mode * coordinate / grid.length)
    quarter = hidden_amplitude * np.sin(
        2.0 * math.pi * hidden_mode * (coordinate + 0.25 * system.dx) / grid.length
    )
    assert float(np.max(np.abs(aligned))) < 1e-12
    assert float(np.max(np.abs(quarter))) > 0.99 * hidden_amplitude
    assert continuum > 1.0
    invisible = observe_discrete_inequality(
        system, y_star + aligned, y_star, G_star, c_seg,
    )
    assert invisible["holds"] is True
    assert invisible["upper_bound"] < 1e-6 * continuum
    assert invisible["total_certified"] is False
    assert invisible["is_continuum_certificate"] is False

    doubled_mode = 2 * int(grid.nq)
    doubled_energy = _mode_energy(hidden_amplitude, doubled_mode, grid.length)
    sample_maxima = []
    for points in (int(grid.nq), 2 * int(grid.nq)):
        samples = np.arange(points, dtype=float) * (grid.length / points)
        values = hidden_amplitude * np.sin(2.0 * math.pi * doubled_mode * samples / grid.length)
        sample_maxima.append(float(np.max(np.abs(values))))
    assert max(sample_maxima) < 1e-12
    refinement_gap = abs(sample_maxima[0] - sample_maxima[1])
    # Both aligned grids miss this sine, so the refinement gap is below the continuum energy.
    assert refinement_gap < 1e-12
    assert refinement_gap < doubled_energy
    eighth = (np.arange(grid.nq, dtype=float) + 0.125) * system.dx
    visible = hidden_amplitude * np.sin(2.0 * math.pi * doubled_mode * eighth / grid.length)
    assert float(np.max(np.abs(visible))) > 0.99 * hidden_amplitude
    assert doubled_energy > continuum


def test_mean_current_remains_and_the_source_array_is_unchanged():
    grid = _grid()
    current = np.full(grid.nq, 0.3)
    snapshot = current.copy()
    momentum, info = shift_momentum(grid, current)
    residual = shift_residual(grid, momentum, current)
    assert np.max(np.abs(current - snapshot)) == 0.0
    assert info["source_mean_subtracted"] is False
    assert abs(float(np.mean(residual)) - 0.3) < 1e-12
    assert float(np.max(np.abs(residual - 0.3))) < 1e-12
    assert float(info["p_Q_max"]) < 1e-12

    phi0, phi1 = manufactured_columns(grid.nf, seed=3)
    state = blank_state(grid, phi0, phi1)
    built = constructor_current(grid, state)
    fine = prolong_state(grid, state)
    source = source_from_columns(grid.fine, fine)
    independent_current = np.asarray(source["force_beta"] / grid.dx_q, dtype=float)
    independent_momentum, _independent_info = shift_momentum(grid, independent_current)
    independent_residual = shift_residual(grid, independent_momentum, independent_current)
    assert abs(float(np.mean(independent_current))) > 1e-2
    assert abs(float(np.mean(independent_residual)) - float(np.mean(independent_current))) < 1e-12
    assert abs(built["current_mean"] - float(np.mean(independent_current))) < 1e-12
    assert built["source_mean_subtracted"] is False
    assert built["projected_momentum_mean_gap"] < 1e-12
    assert built["radius_newton_rerun"] is False
    assert built["current_unchanged"] is True


def test_declared_inputs_reject_source_and_setting_changes_without_a_fake_hash(tmp_path):
    sealed_before = RECORD_PATH.read_bytes()
    saved = json.loads(sealed_before)
    source_ref = "5f10ecd365843d1616e50eb16a20d7acd8377e2c"
    assert binding_report(saved)["ok"] is False
    report = binding_report(saved, source_ref=source_ref)
    assert report["ok"] is True
    assert report["head_equality_required"] is False
    assert report["fabricated_hash"] is False
    assert report["scientific_sources_match"] is True
    assert report["historical_checkpoint_head"] == "9a9090a"
    generator_limits = [item for item in report["limits"] if item["role"] == "generator"]
    assert len(generator_limits) == 3
    assert all(item["fabricated_hash"] is False for item in generator_limits)
    assert all(item["current_hash"] is not None for item in generator_limits)
    assert all(item["current_hash"] != item["declared_hash"] for item in generator_limits)
    assert all(item["historical_bytes_available"] is False for item in generator_limits)

    shifted = copy.deepcopy(saved)
    shifted["checkpoint_head"] = "e62f804"
    shifted["cpu_seconds"] = 0.01
    shifted["wall_seconds"] = 0.02
    shifted_report = binding_report(shifted, source_ref=source_ref)
    assert shifted_report["ok"] is True
    assert shifted_report["head_equality_required"] is False
    assert shifted_report["historical_checkpoint_head"] == "e62f804"

    key = "src/recursive_horizons/nsc_spherical_coupling.py"
    mismatched = copy.deepcopy(saved)
    mismatched["source_hashes"][key] = "0" * 64
    mismatch_report = binding_report(mismatched, source_ref=source_ref)
    assert mismatch_report["ok"] is False
    entry = next(item for item in mismatch_report["blocking"] if item["path"] == key)
    assert entry["fabricated_hash"] is False
    assert entry["current_hash"] is not None
    assert entry["current_hash"] != entry["declared_hash"]
    assert entry["matches_declared"] is False

    absent = copy.deepcopy(saved)
    absent["source_hashes"].pop(key)
    absent_report = binding_report(absent, source_ref=source_ref)
    assert absent_report["ok"] is False
    absent_entry = next(item for item in absent_report["blocking"] if item["path"] == key)
    assert absent_entry["declared_hash"] is None
    assert absent_entry["current_hash"] is not None
    assert absent_entry["fabricated_hash"] is False
    assert absent_entry["matches_declared"] is False

    missing = binding_report(saved, source_root=tmp_path)
    assert missing["ok"] is False
    assert missing["fabricated_hash"] is False
    missing_sources = [
        item for item in missing["blocking"] if item["role"] == "scientific_source"
    ]
    assert missing_sources
    assert all(item["available"] is False for item in missing_sources)
    assert all(item["current_hash"] is None for item in missing_sources)
    assert all(item["declared_hash"] is not None for item in missing_sources)
    assert all(item["fabricated_hash"] is False for item in missing_sources)
    assert all(item["current_hash"] != item["declared_hash"] for item in missing_sources)

    historical_missing = binding_report(saved, source_ref="0" * 40)
    assert historical_missing["ok"] is False
    assert any(item["historical_bytes_available"] is False for item in historical_missing["blocking"])

    schema = copy.deepcopy(saved)
    schema["schema"] = "NSC-SPHERICAL-CAUCHY-WEAK-v0"
    assert binding_report(schema)["ok"] is False
    opened = copy.deepcopy(saved)
    opened["v5_bytes_unchanged"] = False
    assert binding_report(opened)["ok"] is False
    assert RECORD_PATH.read_bytes() == sealed_before
