"""Independent algebra for the Fourier–Galerkin coupling.

These checks recompute pullback weights, symplectic factors, the lifted
work sum, the held-out Hamilton residual, the chart velocity, and the
frozen link. A green result is not a physical pass. The saved v2 verdict
stays FAIL_HELD_OUT_CONSTRAINT.
"""
from __future__ import annotations

import json
from functools import lru_cache

import numpy as np

from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_spherical_coupling import (
    CauchyState,
    field_energy,
    geometric_rates,
    hamilton_constraint,
    source_from_columns,
    total_energy,
)

# Hand-set window link. Compared with the Dirac block, not used as a source.
FROZEN_B = np.array(
    [[0.25, 1j / 7], [1.0 / 9.0, 1.0 / 6.0]],
    dtype=np.complex128,
)


def _hand_periodic(ng, nq, length):
    modes = np.arange(-(ng // 2), ng // 2 + 1, dtype=float)
    coarse = np.arange(ng, dtype=float) * (length / ng)
    fine = np.arange(nq, dtype=float) * (length / nq)
    analysis = np.exp(-2j * np.pi * coarse[:, None] * modes[None, :] / length) / ng
    synthesis = np.exp(2j * np.pi * fine[:, None] * modes[None, :] / length)
    return np.real(synthesis @ analysis.T)


def _hand_antiperiodic(nf, nq, length):
    modes = (np.arange(nf) - nf // 2).astype(float) + 0.5
    coarse = np.arange(nf, dtype=float) * (length / nf)
    fine = np.arange(nq, dtype=float) * (length / nq)
    analysis = np.exp(-2j * np.pi * coarse[:, None] * modes[None, :] / length) / nf
    synthesis = np.exp(2j * np.pi * fine[:, None] * modes[None, :] / length)
    return synthesis @ analysis.T


def _pull(grid, values):
    return (grid.ng / float(grid.nq)) * grid.A_g.T @ np.asarray(values, dtype=float)


def _energy(grid, state):
    return total_energy(grid.fine, galerkin.prolong_state(grid, state))


def _field_energy(grid, state):
    return field_energy(grid.fine, galerkin.prolong_state(grid, state))


def _displace(state, name, direction, scale):
    trial = state.copy()
    setattr(trial, name, getattr(trial, name) + scale * direction)
    return trial


def _directional(grid, state, name, direction, eps, energy_fn):
    forward = energy_fn(grid, _displace(state, name, direction, eps))
    backward = energy_fn(grid, _displace(state, name, direction, -eps))
    return (forward - backward) / (2 * eps)


def _manufactured(fermions=32, quadrature=128, seed=2):
    grid = galerkin.build_grid(fermions, quadrature=quadrature)
    phi0, phi1 = galerkin.manufactured_columns(fermions, seed=3)
    state = galerkin.blank_state(grid, phi0, phi1)
    state.r = 2 + 0.05 * np.cos(2 * np.pi * grid.xi_g / grid.length)
    state.Q = state.Q + 0.15 * np.cos(2 * np.pi * (grid.ng // 2) * grid.xi_g / grid.length)
    generator = np.random.default_rng(seed)
    state.p_Q = 0.01 * generator.normal(size=grid.ng)
    state.p_r = 0.01 * generator.normal(size=grid.ng)
    state.p_chi = 0.01 * generator.normal(size=grid.ng)
    state.chi = 0.01 * np.sin(4 * np.pi * grid.xi_g / grid.length)
    return grid, state


def _hand_rates(grid, state):
    fine = galerkin.prolong_state(grid, state)
    q_dot, r_dot, chi_dot, p_q_dot, p_r_dot, p_chi_dot = geometric_rates(grid.fine, fine)
    source = source_from_columns(grid.fine, fine)
    matter = source["force_Q"] / grid.dx_q
    adjoint = grid.U_f.conj().T
    return {
        "Q": _pull(grid, q_dot),
        "r": _pull(grid, r_dot),
        "chi": _pull(grid, chi_dot),
        "p_Q": _pull(grid, p_q_dot - matter),
        "p_Q_without_matter": _pull(grid, p_q_dot),
        "p_r": _pull(grid, p_r_dot),
        "p_chi": _pull(grid, p_chi_dot),
        "phi0": adjoint @ (-1j * source["image0"]),
        "phi1": adjoint @ (-1j * source["image1"]),
        "phi0_with_multiplicity": adjoint @ (-1j * grid.fine.multiplicity * source["image0"]),
        "force_Q": source["force_Q"],
        "q_dot_fine": q_dot,
        "r_dot_fine": r_dot,
        "fine": fine,
    }


def _radius_probe(grid, radius, rho):
    fine_radius = galerkin.prolong_geometry(grid, radius)
    q0 = coupling.CALIBRATION["b0"] / coupling.CALIBRATION["a0"]
    probe = CauchyState(
        Q=np.full(grid.nq, q0),
        r=fine_radius,
        chi=np.zeros(grid.nq),
        p_Q=np.zeros(grid.nq),
        p_r=np.zeros(grid.nq),
        p_chi=np.zeros(grid.nq),
        phi0=np.zeros((grid.nq, 1)),
        phi1=np.zeros((grid.nq, 1)),
    )
    residual = hamilton_constraint(grid.fine, probe) + np.asarray(rho, dtype=float)
    projected = _pull(grid, residual)
    held = residual - grid.A_g @ projected
    return residual, projected, held


def _link_blocks(grid, state):
    fine = galerkin.prolong_state(grid, state)
    source = source_from_columns(grid.fine, fine)
    matrix = fine.phi0.conj().T @ source["image0"] + fine.phi1.conj().T @ source["image1"]
    blocks = []
    for region in range(coupling.PACKET_COUNT - 1):
        block = matrix[2 * region : 2 * region + 2, 2 * region + 2 : 2 * region + 4]
        target = coupling.OMEGA ** region * FROZEN_B
        blocks.append(np.linalg.norm(block - target))
    return blocks


def _velocity_split(grid, state):
    fine = galerkin.prolong_state(grid, state)
    _q_dot, r_dot, _chi, _pq, _pr, _pc = geometric_rates(grid.fine, fine)
    radius_derivative = grid.derivative @ fine.r
    denominator = fine.r * grid.fine.length_density
    chart = (r_dot - grid.fine.shift * radius_derivative) / denominator
    lifted = (grid.A_g @ _pull(grid, r_dot) - grid.fine.shift * radius_derivative) / denominator
    return chart, lifted


@lru_cache(maxsize=1)
def _physical_initial():
    phi0, phi1, preparation, problems = galerkin.load_physical_columns(64)
    grid = galerkin.build_grid(64, quadrature=256)
    solved, info = galerkin.solve_initial_radius(grid, galerkin.blank_state(grid, phi0, phi1))
    return grid, solved, info, preparation, tuple(problems)


def test_variational_pullback_uses_the_nodal_isometry_and_half_density_columns():
    length = galerkin.PERIOD
    nf, nq = 16, 64
    ng = nf - 1
    grid = galerkin.build_grid(nf, quadrature=nq)
    geometry = _hand_periodic(ng, nq, length)
    fermions = _hand_antiperiodic(nf, nq, length)
    assert np.allclose(grid.A_g, geometry)
    assert np.allclose(grid.A_f, fermions)
    weight = ng / float(nq)
    columns = np.sqrt(nf / float(nq)) * fermions
    assert np.allclose(grid.U_f, columns)
    assert abs(grid.weight - weight) < 1e-15
    assert abs(grid.dx_g - grid.dx_q / grid.weight) < 1e-12
    assert np.allclose(weight * geometry.T @ geometry, np.eye(ng), atol=1e-10)
    assert np.linalg.norm((nq / float(ng)) * geometry.T @ geometry - np.eye(ng)) > 1
    generator = np.random.default_rng(4)
    sample = generator.normal(size=ng)
    dx_g = length / ng
    dx_q = length / nq
    lifted_mass = dx_q * np.sum((geometry @ sample) ** 2)
    coarse_mass = dx_g * np.sum(sample ** 2)
    assert abs(lifted_mass - coarse_mass) < 1e-8 * coarse_mass
    probe = generator.normal(size=nq)
    pull = weight * geometry.T @ probe
    assert abs(weight * np.dot(geometry @ sample, probe) - np.dot(sample, pull)) < 1e-8
    raw = generator.normal(size=nf) + 1j * generator.normal(size=nf)
    coefficient = raw / np.linalg.norm(raw)
    dx_f = length / nf
    half_density = (fermions @ (coefficient / np.sqrt(dx_f))) * np.sqrt(dx_q)
    assert np.allclose(half_density, columns @ coefficient)
    assert abs(np.vdot(columns @ coefficient, columns @ coefficient) - 1) < 1e-10
    inverted = np.sqrt(nq / float(nf)) * fermions @ coefficient
    assert abs(np.vdot(inverted, inverted) - 1) > 0.5


def test_symplectic_weight_is_coarse_spacing_and_columns_are_not_scaled_by_4kappa():
    grid, state = _manufactured()
    hand = _hand_rates(grid, state)
    rate = galerkin.rates(grid, state)
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1"):
        assert np.allclose(getattr(rate, name), hand[name], atol=1e-10)
    assert grid.fine.multiplicity == 4 * int(grid.fine.kappa) == 4
    column_gap = np.linalg.norm(rate.phi0 - hand["phi0_with_multiplicity"])
    assert column_gap > 0.5 * np.linalg.norm(rate.phi0)
    generator = np.random.default_rng(11)
    pairs = {
        "Q": ("p_Q", -grid.dx_g),
        "r": ("p_r", -grid.dx_g),
        "chi": ("p_chi", -grid.dx_g),
        "p_Q": ("Q", grid.dx_g),
        "p_r": ("r", grid.dx_g),
        "p_chi": ("chi", grid.dx_g),
    }
    absolute = []
    wrong = []
    scales = []
    for name, (velocity_name, factor) in pairs.items():
        direction = generator.normal(size=grid.ng)
        direction /= np.linalg.norm(direction)
        numeric = _directional(grid, state, name, direction, 1e-6, _energy)
        velocity = hand[velocity_name]
        analytic = factor * float(np.dot(velocity, direction))
        swapped = (grid.dx_q / grid.dx_g) * analytic
        absolute.append(abs(numeric - analytic))
        wrong.append(abs(numeric - swapped))
        scales.append(abs(analytic))
    scale = max(1.0, max(scales))
    assert max(absolute) / scale < galerkin.TOL_FD_RELATIVE
    assert min(wrong) / scale > 1e-3
    direction = generator.normal(size=grid.ng)
    direction /= np.linalg.norm(direction)
    numeric_force = _directional(grid, state, "Q", direction, 1e-6, _field_energy)
    coarse_gradient = grid.A_g.T @ hand["force_Q"]
    weighted_gradient = grid.weight * coarse_gradient
    force_scale = max(1.0, abs(float(np.dot(coarse_gradient, direction))))
    assert abs(numeric_force - float(np.dot(coarse_gradient, direction))) / force_scale < 1e-6
    assert abs(numeric_force - float(np.dot(weighted_gradient, direction))) / force_scale > 1e-2


def test_lifted_fieldwork_matches_the_energy_slope_and_the_fine_rate_does_not():
    grid, state = _manufactured()
    hand = _hand_rates(grid, state)
    dt = 1e-6

    def slope(momentum):
        def at(sign):
            trial = CauchyState(
                state.Q + sign * dt * hand["Q"],
                state.r + sign * dt * hand["r"],
                state.chi + sign * dt * hand["chi"],
                state.p_Q + sign * dt * momentum,
                state.p_r + sign * dt * hand["p_r"],
                state.p_chi + sign * dt * hand["p_chi"],
                state.phi0 + sign * dt * hand["phi0"],
                state.phi1 + sign * dt * hand["phi1"],
            )
            return _energy(grid, trial)

        return (at(1) - at(-1)) / (2 * dt)

    intact = slope(hand["p_Q"])
    mutated = slope(hand["p_Q_without_matter"])
    lifted_q = grid.A_g @ hand["Q"]
    fieldwork = float(np.sum(hand["force_Q"] * lifted_q))
    unprojected = float(np.sum(hand["force_Q"] * hand["q_dot_fine"]))
    spaced = float(grid.dx_q * np.sum(hand["force_Q"] * lifted_q))
    lifted_error = abs((mutated - intact) - fieldwork)
    unprojected_error = abs((mutated - intact) - unprojected)
    spaced_error = abs((mutated - intact) - spaced)
    coarse, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    assert abs(coarse.fieldwork_power - fieldwork) < 1e-8
    assert abs(bundle["unprojected_fieldwork"] - unprojected) < 1e-8
    assert lifted_error < 1e-6
    assert unprojected_error > 1.0
    assert spaced_error > 1.0
    assert abs(fieldwork - unprojected) > 1.0


def test_projected_initial_residual_is_not_the_full_constraint():
    grid, state = _manufactured()
    mode = grid.ng // 2
    rho = 5.0 * np.cos(2 * np.pi * mode * grid.xi_q / grid.length)
    solved, info = galerkin.solve_initial_radius(grid, state, rho=rho)
    residual, projected, held = _radius_probe(grid, solved.r, rho)
    full_max = float(np.max(np.abs(residual)))
    projected_max = float(np.max(np.abs(projected)))
    held_max = float(np.max(np.abs(held)))
    assert info["converged"] is True
    assert projected_max < galerkin.TOL_INITIAL_CONSTRAINT
    assert full_max > galerkin.TOL_CONSTRAINT_VALIDATION
    assert held_max > 0.99 * full_max
    assert abs(info["diagnostics"]["projected_hamilton_max"] - projected_max) < 1e-8
    assert abs(info["diagnostics"]["full_hamilton_max"] - full_max) < 1e-8
    low = 0.02 * np.cos(2 * np.pi * grid.xi_q / grid.length)
    solved_low, info_low = galerkin.solve_initial_radius(grid, state, rho=low)
    low_residual, low_projected, _low_held = _radius_probe(grid, solved_low.r, low)
    assert info_low["converged"] is True
    assert float(np.max(np.abs(low_projected))) < galerkin.TOL_INITIAL_CONSTRAINT
    assert float(np.max(np.abs(low_residual))) < galerkin.TOL_INITIAL_CONSTRAINT
    assert galerkin.RENEWAL is False


def test_lifted_normal_velocity_is_not_the_chart_expression():
    grid, state = _manufactured()
    state.p_Q = np.zeros(grid.ng)
    state.p_r = np.zeros(grid.ng)
    state.p_chi = np.zeros(grid.ng)
    state.chi = np.zeros(grid.ng)
    state.r = 2 + 0.3 * np.cos(2 * np.pi * (grid.ng // 2) * grid.xi_g / grid.length)
    chart, lifted = _velocity_split(grid, state)
    reported = galerkin.normal_velocities(grid, state)
    assert float(np.max(np.abs(chart))) < 1e-10
    assert float(np.max(np.abs(lifted))) > 1e-4
    assert abs(reported["chart_proper_max"] - float(np.max(np.abs(chart)))) < 1e-12
    assert abs(reported["lifted_proper_max"] - float(np.max(np.abs(lifted)))) < 1e-12
    assert reported["lifted_proper_max"] > 100.0 * max(reported["chart_proper_max"], 1e-16)


def test_physical_slice_keeps_the_held_out_failure_and_a_different_link():
    grid, solved, info, preparation, problems = _physical_initial()
    assert problems == ()
    assert preparation["phase_on_minus_spinor_column"] is True
    assert preparation["phase_applied_to_odd_lobe"] is False
    fine = galerkin.prolong_state(grid, solved)
    source = source_from_columns(grid.fine, fine)
    rho = source["force_L"] / grid.dx_q
    residual = hamilton_constraint(grid.fine, fine) + rho
    projected = _pull(grid, residual)
    held = residual - grid.A_g @ projected
    full_max = float(np.max(np.abs(residual)))
    projected_max = float(np.max(np.abs(projected)))
    held_max = float(np.max(np.abs(held)))
    assert info["converged"] is True
    assert info["filtered_v1_radius"] is False
    assert projected_max < galerkin.TOL_INITIAL_CONSTRAINT
    assert full_max > galerkin.TOL_CONSTRAINT_VALIDATION
    assert held_max > 0.99 * full_max
    assert np.max(np.abs(solved.p_r)) == 0
    chart, lifted = _velocity_split(grid, solved)
    assert float(np.max(np.abs(chart))) < 1e-10
    assert float(np.max(np.abs(lifted))) > 1e-6
    assert float(np.max(np.abs(lifted))) > 100.0 * float(np.max(np.abs(chart)) + 1e-16)
    gaps = _link_blocks(grid, solved)
    assert max(gaps) > 1.0
    assert all(np.isfinite(gap) for gap in gaps)
    report = galerkin.physical_link_report(grid, solved)
    assert abs(report["max_difference_from_old_B"] - max(gaps)) < 1e-8
    assert report["B_phys_equals_old_B"] is False
    assert report["copied_from_old_B"] is False
    assert report["old_circulation_transferred"] is False
    varied_radius = solved.copy()
    varied_radius.r = solved.r + 0.05 * np.cos(2 * np.pi * grid.xi_g / grid.length)
    radius_shift = max(
        abs(left - right)
        for left, right in zip(gaps, _link_blocks(grid, varied_radius), strict=True)
    )
    varied_density = solved.copy()
    varied_density.Q = solved.Q * (1.0 + 0.2 * np.cos(2 * np.pi * grid.xi_g / grid.length))
    density_shift = max(
        abs(left - right)
        for left, right in zip(gaps, _link_blocks(grid, varied_density), strict=True)
    )
    assert radius_shift < 1e-8
    assert density_shift > 1e-2
    record = json.loads(galerkin.CONTROL_RECORD.read_text(encoding="utf-8"))
    assert record["verdict"] == "FAIL_HELD_OUT_CONSTRAINT"
    assert record["renewal"] is False
    assert galerkin.RENEWAL is False
