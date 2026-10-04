"""Manufactured checks for the neutral RN packet and its common-action source.

No production trajectory and no scientific record is written. Tolerances are
the programme's 0.01 manufactured-check bound except where an identity is exact.
"""
from __future__ import annotations

import os

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_name, "1")

import numpy as np
import pytest

from recursive_horizons.nsc_covariant_operator import SIGMA1, SIGMA2
from recursive_horizons.nsc_rn_source import (
    ENERGY_TARGETS,
    SHELL_MULTIPLICITY,
    action_reference,
    characteristic_transform,
    future_horizon_initial,
    magnetic_radius,
    mode_quadratic_form,
    ode_residual,
    packet_request,
    prepare_radial_packet,
    radial_sbp,
    scattering_mode,
    select_packet_filling,
    shell_multiplicity,
    source_evaluate,
    vacuum_chart,
    weighted_gram,
)


def test_characteristic_basis_and_shell_count():
    transform = characteristic_transform()
    diagonal = transform.conj().T @ SIGMA2 @ transform
    angular = transform.conj().T @ SIGMA1 @ transform
    assert np.allclose(diagonal, np.diag([1.0, -1.0]), atol=1e-12)
    assert np.allclose(angular, SIGMA2, atol=1e-12)
    assert shell_multiplicity(1) == 4
    with pytest.raises(ValueError):
        shell_multiplicity(0)
    with pytest.raises(ValueError):
        shell_multiplicity(-1)


def test_sbp_closure_is_nonperiodic():
    operator = radial_sbp(21, 2.25, 8.0)
    assert operator["closure_residual"] < 1e-12
    assert operator["periodic"] is False
    assert operator["weights"][0] == pytest.approx(operator["step"] / 2.0)


def test_rn_lapse_is_derived_from_the_first_integral():
    model = action_reference(1.04)
    r_m = magnetic_radius(model)
    assert r_m == pytest.approx(1.0)
    assert model.horizons["r_minus"] != pytest.approx(r_m)
    radius = np.linspace(float(model.horizons["r_plus"]) + 0.2, 12.0, 30)
    chart = vacuum_chart(model, radius)
    assert chart["lapse_derived"] is True
    assert chart["gauge_fix_imposed"] is False
    assert np.allclose(chart["q_adm"], 1.0)
    assert np.allclose(chart["conformal_Q"], 1.0 / radius)
    assert np.allclose(chart["lapse"]**2 - chart["shift"]**2, chart["redshift"], atol=1e-12)
    assert np.allclose(chart["lapse"], 1.0, atol=1e-12)
    assert np.allclose(
        chart["shift"]**2, 2.0 * model.mass / radius - model.magnetic_r2 / radius**2, atol=1e-12,
    )
    with pytest.raises(ValueError):
        action_reference(r_m)
    with pytest.raises(ValueError):
        vacuum_chart(model, np.array([0.5 * float(model.horizons["r_minus"])]))


def test_future_horizon_root_and_mode_match_the_static_ode():
    model = action_reference(1.04)
    energy = 1.0
    horizon = float(model.horizons["r_plus"])
    initial, chart = future_horizon_initial(model, energy)
    angular = chart["lapse"] / horizon
    denominator = 1j * energy - 0.5 * chart["speed_plus_derivative"]
    assert initial[0] / initial[1] == pytest.approx(angular / denominator, rel=0.0, abs=1e-12)
    with pytest.raises(ValueError):
        future_horizon_initial(model, energy, radius=horizon + 0.2)
    grid = np.linspace(horizon + 0.3, 10.0, 81)
    solved = scattering_mode(grid, energy, model)
    assert solved["absorbing_matrix_projection"] is False
    assert solved["nodes_dropped"] == 0
    assert solved["interior_nodes"] == 0
    assert ode_residual(model, solved["mode"], grid, energy) < 0.01
    both = np.linspace(
        0.5 * (model.horizons["r_minus"] + horizon), horizon + 2.0, 33,
    )
    crossed = scattering_mode(both, energy, model)
    assert crossed["interior_nodes"] > 0
    assert crossed["nodes_dropped"] == 0
    assert crossed["node_count"] == both.size
    weights = radial_sbp(grid.size, grid[0], grid[-1])["weights"]
    norm = np.sqrt(np.real(np.sum(weights * np.sum(np.abs(solved["mode"])**2, axis=0))))
    mode = solved["mode"] / norm
    operator = radial_sbp(grid.size, grid[0], grid[-1])
    chart = vacuum_chart(model, grid)
    formed = mode_quadratic_form(
        mode, chart["lapse"], chart["shift"], chart["q_adm"],
        chart["radius"], operator["weights"], operator["derivative"],
    )
    assert formed["energy"] == pytest.approx(energy, rel=0.01, abs=0.01)
    current = formed["probability_current"]
    assert (np.max(current) - np.min(current)) / np.max(np.abs(current)) < 0.01
    assert np.max(chart["speed_minus"]) < 0.0
    assert np.min(chart["speed_plus"]) > 0.0


def test_packet_filling_uses_the_measured_energy():
    newton = 0.4
    filling, actual = select_packet_filling(0.5, 2.0, 1.0e-3, newton)
    assert newton * SHELL_MULTIPLICITY * filling * 0.5 / 2.0 == pytest.approx(1.0e-3)
    assert actual == pytest.approx(1.0e-3)
    with pytest.raises(ValueError):
        select_packet_filling(0.5, 2.0, 0.2, newton)
    request = packet_request(1.5)
    assert request["center"] == pytest.approx(18.0)
    assert request["width"] == pytest.approx(3.0)
    assert request["omega"] == pytest.approx(1.0 / 1.5)
    assert request["r_m"] == pytest.approx(1.5)


def test_caller_grid_is_not_resampled():
    model = action_reference(1.04)
    r_m = magnetic_radius(model)
    request = packet_request(r_m)
    left = 0.5 * (float(model.horizons["r_minus"]) + float(model.horizons["r_plus"]))
    right = request["center"] + 4.0 * request["width"]
    operator = radial_sbp(33, left, right)
    with pytest.raises(ValueError, match="resample"):
        prepare_radial_packet(
            mass_over_rm=1.04,
            radii=operator["grid"],
            weights=operator["weights"],
            model=model,
        )


@pytest.fixture(scope="module")
def packet():
    model = action_reference(1.04)
    r_m = magnetic_radius(model)
    request = packet_request(r_m)
    left = 0.5 * (float(model.horizons["r_minus"]) + float(model.horizons["r_plus"]))
    right = request["center"] + 4.0 * request["width"]
    operator = radial_sbp(97, left, right)
    prepared = prepare_radial_packet(
        mass_over_rm=1.04,
        energy_target=1.0e-3,
        rank=5,
        model=model,
        radii=operator["grid"],
        weights=operator["weights"],
        derivative=operator["derivative"],
    )
    return prepared


def test_packet_positivity_gram_and_matched_basis(packet):
    phi = packet["phi"]
    assert phi.shape == (2, packet["background"]["grid"].size, 1)
    assert packet["spectral_phi"].shape[2] == 5
    assert packet["absorbing_matrix_projection"] is False
    assert packet["periodic_identification"] is False
    assert packet["background"]["grid_binding"] == "caller"
    assert np.array_equal(packet["background"]["grid"], packet["background"]["radius"])
    assert packet["spectrum"]["negative_frequency_contamination"] == 0.0
    assert np.min(packet["spectrum"]["requested"]) > 0.0
    assert np.min(packet["spectrum"]["quadratic_form"]) > 0.0
    assert packet["packet_mode_energy"] > 0.0
    assert 0.0 < packet["packet_filling"] < 1.0
    assert packet["occupations"].shape == (1,)
    assert 0.0 < packet["occupations"][0] < 1.0
    gram = weighted_gram(packet["spectral_phi"], packet["background"]["weights"])
    assert np.allclose(gram, packet["car"]["gram"], atol=1e-12)
    assert np.allclose(np.real(np.diag(gram)), 1.0, atol=1e-10)
    assert np.allclose(gram, gram.conj().T, atol=1e-10)
    assert packet["car"]["multiplicity_in_car"] is False
    car_eigenvalues = np.linalg.eigvalsh(packet["car"]["physical"])
    assert int(np.count_nonzero(np.abs(car_eigenvalues) > 1e-12)) == 1
    assert np.max(np.abs(car_eigenvalues)) < 1.0
    requested = packet["spectrum"]["requested"]
    measured = packet["spectrum"]["quadratic_form"]
    carrier = int(np.argmin(np.abs(requested - packet["radial_moments"]["requested_omega"])))
    for index, energy in enumerate(requested):
        residual = ode_residual(
            packet["background"]["model"], packet["spectral_phi"][:, :, index],
            packet["background"]["grid"], energy,
        )
        assert residual < 0.01
        # Edge samples use the second-order SBP energy. The carrier is the tight identity.
        assert measured[index] == pytest.approx(energy, rel=0.02, abs=0.03)
    assert measured[carrier] == pytest.approx(requested[carrier], rel=0.01, abs=0.01)
    ratios = packet["near_horizon_tail"]["horizon_plus_over_minus"]
    assert np.all(np.isfinite(ratios))
    assert np.all(np.isfinite(packet["near_horizon_tail"]["packet_plus"]))


def test_angular_trace_is_four_once_and_energy_hits_both_targets(packet):
    covariance = packet["angular_covariance"]
    assert covariance.shape == (4, 4)
    assert np.allclose(covariance, np.eye(4) * packet["packet_filling"])
    assert np.trace(covariance) == pytest.approx(4.0 * packet["packet_filling"])
    assert packet["multiplicity"] == 4
    assert packet["killing_energy"] == pytest.approx(
        4.0 * packet["packet_filling"] * packet["packet_mode_energy"],
    )
    assert packet["actual_energy_target"] == pytest.approx(1.0e-3, rel=1e-12, abs=1e-15)
    assert packet["G_N"] != pytest.approx(1.0)
    measured = packet["spectrum"]["quadratic_form"]
    sigma = packet["spectrum"]["spectral_width"]
    omega = packet["radial_moments"]["requested_omega"]
    requested = packet["spectrum"]["requested"]
    gaussian = np.exp(-0.25 * ((requested - omega) / sigma) ** 2)
    coefficient = packet["spectral_coefficient_matrix"]
    alpha = coefficient[:, 0]
    scale = alpha[0] / gaussian[0]
    assert np.allclose(alpha, scale * gaussian, rtol=1e-10, atol=1e-14)
    assert measured.shape == requested.shape
    other = prepare_radial_packet(
        mass_over_rm=1.04,
        energy_target=1.0e-2,
        rank=5,
        model=packet["background"]["model"],
        radii=packet["background"]["grid"],
        weights=packet["background"]["weights"],
        derivative=packet["background"]["derivative"],
    )
    assert other["actual_energy_target"] == pytest.approx(1.0e-2, rel=1e-12, abs=1e-15)
    assert set(ENERGY_TARGETS) == {1.0e-3, 1.0e-2}


def test_packet_is_centered_on_the_requested_exterior_window(packet):
    moments = packet["radial_moments"]
    assert moments["center"] == pytest.approx(moments["requested_center"], abs=0.5)
    assert moments["width_rms"] == pytest.approx(moments["requested_width"], rel=0.45, abs=0.2)
    assert moments["width_is_fit"] is False
    tail = np.sum(np.abs(packet["near_horizon_tail"]["packet_plus"])**2
                  + np.abs(packet["near_horizon_tail"]["packet_minus"])**2)
    peak = np.max(np.sum(np.abs(packet["phi"][:, :, 0])**2, axis=0))
    assert tail < peak


def test_source_matches_variation_trace_and_packet_energy(packet):
    background = packet["background"]
    evaluated = source_evaluate(
        packet["phi"], packet["occupations"], background["lapse"], background["shift"],
        background["radial_density"], background["radius"], background["weights"],
        background["derivative"],
    )
    assert evaluated["neutral"] is True
    assert evaluated["magnetic_flux_fixed"] is True
    assert evaluated["weyl_variation_included"] is False
    assert evaluated["dense_covariance"] is False
    assert evaluated["gauge_fix_imposed"] is False
    assert evaluated["coherent_rank_one_column"] is True
    assert np.all(evaluated["electric_current"] == 0.0)
    assert evaluated["multiplicity"] == 4
    assert evaluated["killing_energy"] == pytest.approx(packet["killing_energy"], rel=1e-10, abs=1e-12)
    assert np.max(np.abs(evaluated["massless_trace_residual"])) < 1e-10
    rho = evaluated["physical"]["rho"]
    force_n = evaluated["forces"]["N"]
    volume = background["radial_density"] * background["radius"]**2
    assert np.max(np.abs(rho * 4.0 * np.pi * volume - force_n)) < 1e-10
    weights = background["weights"]
    base = evaluated["killing_energy"]
    for name, field in (("N", "lapse"), ("beta", "shift")):
        delta = 0.05
        changed = background[field] + delta
        arguments = dict(
            lapse=background["lapse"], shift=background["shift"],
            radial_density=background["radial_density"], radius=background["radius"],
        )
        arguments[field] = changed
        moved = source_evaluate(
            packet["phi"], packet["occupations"], arguments["lapse"], arguments["shift"],
            arguments["radial_density"], arguments["radius"], weights,
            background["derivative"],
        )
        contraction = float(np.sum(weights * evaluated["forces"][name])) * delta
        assert moved["killing_energy"] - base == pytest.approx(contraction, rel=1e-9, abs=1e-12)


def test_radius_force_is_taken_before_the_one_over_r_gauge():
    operator = radial_sbp(33, 0.0, 2.0)
    grid = operator["grid"]
    radius = np.full_like(grid, 4.0)
    lapse = np.full_like(grid, 1.2)
    shift = np.zeros_like(grid)
    density = np.full_like(grid, 0.7)
    profile = np.zeros((2, grid.size), dtype=complex)
    profile[0] = np.exp(1j * grid)
    profile[1] = 0.3 * np.exp(-1j * 0.5 * grid)
    phi = profile[:, :, None]
    occupations = np.array([0.02])
    evaluated = source_evaluate(
        phi, occupations, lapse, shift, density, radius,
        operator["weights"], operator["derivative"],
    )
    weights = operator["weights"]
    step = 1.0e-4

    def energy(next_radius, next_density):
        return source_evaluate(
            phi, occupations, lapse, shift, next_density, next_radius,
            weights, operator["derivative"],
        )["killing_energy"]

    partial_radius = (energy(radius + step, density) - energy(radius - step, density)) / (2.0 * step)
    gauge_up = 1.0 / (radius + step)
    gauge_down = 1.0 / (radius - step)
    directional = (energy(radius + step, gauge_up) - energy(radius - step, gauge_down)) / (2.0 * step)
    force_radius = float(np.sum(weights * evaluated["forces"]["r"]))
    assert partial_radius == pytest.approx(force_radius, rel=1e-6, abs=1e-8)
    assert abs(directional - force_radius) > 1e-6


def test_manufactured_cylinder_pressure_and_noether_identity():
    points = 241
    length = 8.0
    operator = radial_sbp(points, 0.0, length)
    grid = operator["grid"]
    kappa = 1
    radius_value = 5.0
    wave = 1.0
    angular = kappa / radius_value
    energy = float(np.sqrt(wave * wave + angular * angular))
    ratio = 1j * (energy - wave) / angular
    mode = np.zeros((2, points), dtype=complex)
    phase = np.exp(1j * wave * grid)
    mode[0] = phase
    mode[1] = ratio * phase
    weights = operator["weights"]
    norm = np.sqrt(np.real(np.sum(weights * np.sum(np.abs(mode)**2, axis=0))))
    mode = mode / norm
    lapse = np.ones(points)
    shift = np.zeros(points)
    density = np.ones(points)
    radius = np.full(points, radius_value)
    formed = mode_quadratic_form(
        mode, lapse, shift, density, radius, weights, operator["derivative"],
    )
    assert formed["energy"] == pytest.approx(energy, rel=0.01, abs=0.01)
    current = formed["probability_current"]
    assert (np.max(current) - np.min(current)) / np.max(np.abs(current)) < 0.01
    evaluated = source_evaluate(
        mode[:, :, None], np.array([0.04]), lapse, shift, density, radius,
        weights, operator["derivative"],
    )
    assert evaluated["killing_energy"] == pytest.approx(4.0 * 0.04 * formed["energy"], rel=1e-8, abs=1e-8)
    assert np.max(np.abs(evaluated["massless_trace_residual"])) < 1e-9
    assert np.max(np.abs(evaluated["angular_pressure"])) > 0.0
