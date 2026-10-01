"""Independent controls for the spherical feedback episode.

The shared normal-observer split is checked on a small manufactured chart,
not by replaying T=0.05. Proper radial velocity is (r_dot - beta r_x) / (r L).
Observer K_perp divides that velocity by the areal radius, so the clock is
the lapse N = r L. Coordinate work is F_Q Q_dot. Observer pressure work is
-N V (p_r K_r + 2 p_perp K_perp) dx. A refinement movement is one percent
of the claimed physical change. A strong-constraint discrepancy is a separate
report. It does not reclassify a resolved physical row unless a derived
sensitivity says that discrepancy can move the row by more than one percent
of the row's own change.

These checks do not start the production episode and do not retune floors.
"""
from __future__ import annotations

import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np

import derive_nsc_spherical_feedback_episode as episode
from recursive_horizons import nsc_regional_energy_exchange as regional
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_spherical_coupling import (
    source_from_columns,
    hamilton_constraint,
    rates,
)

# Hand-set window link. It is a comparison target, not a source term.
FROZEN_B = np.array([[0.25, 1j / 7.0], [1.0 / 9.0, 1.0 / 6.0]], dtype=np.complex128)

# Scales read from the saved episode record's per-observable rows.
# proper_velocity_max agreed under the nf256 comparison.
PROPER_MAX_CHANGE = 0.08815932373520075
PROPER_MAX_SPACE_MOVEMENT = 3.2511467587814646e-06
# held-out Hamilton change fails one percent of itself by a 7.5e-8 movement.
HELD_OUT_CHANGE = 5.530472490607403e-06
HELD_OUT_TIME_MOVEMENT = 7.548507730883449e-08
# Coordinate-energy change sits just above its floor and disagrees with itself.
ENERGY_CHANGE = 3.442300311462532e-08
ENERGY_TIME_MOVEMENT = 3.3211144057077036e-08

BLANKET_PHYSICAL_FAILURE = "EPISODE_MEASURED_EFFECT_UNRESOLVED"


def _grid_view(system):
    return type("GridView", (), {"fine": system})()


def _proper_velocity(r_dot, radius, radius_x, length, shift):
    """(r_dot - beta r_x) / (r L). N = r L is the normal observer's lapse."""
    return (r_dot - shift * radius_x) / (radius * length)


def _observer_expansions(state, rate, system):
    """K_r and K_perp per unit proper time of the normal observer."""
    radius = np.asarray(state.r, dtype=float)
    radial = np.asarray(state.Q, dtype=float)
    length = np.asarray(system.length_density, dtype=float)
    shift = np.asarray(system.shift, dtype=float)
    lapse = radius * length
    areal = radius * radial
    derivative = system.derivative
    radius_x = derivative @ radius
    q_dot = rate.r * radial + radius * rate.Q
    k_r = (q_dot - derivative @ (shift * areal)) / (lapse * areal)
    k_perp = (rate.r - shift * radius_x) / (lapse * radius)
    return k_r, k_perp


def _stresses(state, system, source):
    radius = np.asarray(state.r, dtype=float)
    radial = np.asarray(state.Q, dtype=float)
    length = np.asarray(system.length_density, dtype=float)
    spacing = float(system.dx)
    sphere = 4.0 * np.pi * radius ** 4
    force_l = source["force_L"]
    force_q = source["force_Q"]
    return {
        "rho": force_l / (sphere * radial * spacing),
        "p_r": -force_q / (sphere * length * spacing),
        "p_perp": (length * force_l + radial * force_q) / (2.0 * sphere * length * radial * spacing),
        "volume": 4.0 * np.pi * (radius * radial) * radius ** 2,
    }


def _coordinate_channels(state, rate, system, window):
    """Window flux and coordinate work from the nodal balance, not from the episode assembler."""
    source = source_from_columns(system, state)
    length = np.asarray(system.length_density, dtype=float)
    shift = np.asarray(system.shift, dtype=float)
    radial = np.asarray(state.Q, dtype=float)
    energy = length * source["force_L"] + shift * source["force_beta"]
    phi_shift = -shift * energy
    phi_proper = float(system.multiplicity) * (length / radial) ** 2 * source["Pmom"]
    phi_cross = shift * radial * source["force_Q"]
    work = source["force_Q"] * rate.Q
    derivative = system.derivative @ window
    return {
        "energy": energy,
        "shift_transport": float(np.sum(derivative * phi_shift)),
        "proper_normal_flux": float(np.sum(derivative * phi_proper)),
        "shift_pressure_cross": float(np.sum(derivative * phi_cross)),
        "coordinate_metric_work": float(np.sum(window * work)),
        "source": source,
        "phi_proper": phi_proper,
    }


def _pressure_work(state, rate, system, source):
    length = np.asarray(system.length_density, dtype=float)
    lapse = np.asarray(state.r, dtype=float) * length
    stresses = _stresses(state, system, source)
    k_r, k_perp = _observer_expansions(state, rate, system)
    return -lapse * stresses["volume"] * (
        stresses["p_r"] * k_r + 2.0 * stresses["p_perp"] * k_perp
    ) * float(system.dx)


def _difference_quotient(samples, step):
    forward, backward = samples
    return (forward - backward) / (2.0 * step)


def _series(key, change):
    return {"time": np.array([0.0, 0.05]), key: np.array([0.0, change])}


def _row(name, key, floor, change, movement):
    primary = _series(key, change)
    other = _series(key, change - movement)
    return episode.effect_row(name, key, floor, primary, other, episode.DT, episode.DURATION)


def _sensitivity_reason(verdict):
    """A derived reason is an explicit map from one diagnostic onto one effect.

    The verdict string is not that map. A reason has to name the effect, the
    diagnostic movement, and the one-percent threshold it was compared with.
    """
    if not isinstance(verdict, dict):
        return None
    reason = verdict.get("sensitivity")
    if not isinstance(reason, dict):
        return None
    required = ("effect", "diagnostic", "diagnostic_movement", "effect_scale", "threshold", "blocks")
    if any(key not in reason for key in required):
        return None
    if abs(float(reason["threshold"]) - 0.01 * abs(float(reason["effect_scale"]))) > 1e-15:
        return None
    return reason


def _blanket_physical_failure(verdict):
    if isinstance(verdict, dict):
        label = verdict.get("verdict", verdict.get("physical_verdict"))
        if verdict.get("effect_resolved") is False and verdict.get("physical_summary") == "resolved":
            return _sensitivity_reason(verdict) is None
        return label == BLANKET_PHYSICAL_FAILURE and _sensitivity_reason(verdict) is None
    return verdict == BLANKET_PHYSICAL_FAILURE


def _manufactured_galerkin():
    grid = galerkin.build_grid(32, quadrature=128)
    phi0, phi1 = galerkin.manufactured_columns(32, seed=3)
    state = galerkin.blank_state(grid, phi0, phi1)
    coordinate = grid.xi_g
    state.r = 4.0 + 0.05 * np.cos(2 * np.pi * coordinate / grid.length)
    state.Q = state.Q * (1.0 + 0.02 * np.cos(2 * np.pi * coordinate / grid.length))
    state.chi = 0.01 * np.sin(2 * np.pi * coordinate / grid.length)
    state.p_r = 0.02 * np.sin(2 * np.pi * coordinate / grid.length)
    state.p_Q = 0.01 * np.cos(2 * np.pi * coordinate / grid.length)
    state.p_chi = 0.005 * np.cos(4 * np.pi * coordinate / grid.length)
    return grid, state


def test_proper_motion_uses_the_lapse_and_keeps_the_radial_sign():
    system, state = regional.band_limited_state(64)
    rate = rates(system, state, include_matter_force=True)
    radius_x = system.derivative @ state.r
    expected = _proper_velocity(rate.r, state.r, radius_x, system.length_density, system.shift)
    wrong_sign = (system.shift * radius_x - rate.r) / (state.r * system.length_density)
    missing_radius = (rate.r - system.shift * radius_x) / system.length_density
    motion = episode.signed_proper_motion(_grid_view(system), state, rate, rate.r)
    assert np.max(np.abs(motion["proper"] - expected)) < 1e-12
    assert np.max(np.abs(expected - wrong_sign)) > 1e-3
    assert np.max(np.abs(expected - missing_radius)) > 1e-3
    assert np.max(np.abs(motion["proper"] + motion["shift_piece"] - motion["coordinate_over_N"])) < 1e-12
    k_r, k_perp = _observer_expansions(state, rate, system)
    assert np.max(np.abs(k_perp - expected / state.r)) < 1e-12
    coordinate_expansion = (rate.r - system.shift * radius_x) / state.r
    assert np.max(np.abs(k_perp - coordinate_expansion)) > 1e-3
    assert np.max(np.abs(k_r)) > 1e-6

    zero = state.copy()
    zero.p_Q = np.zeros_like(state.p_Q)
    zero.p_r = np.zeros_like(state.p_r)
    zero.p_chi = np.zeros_like(state.p_chi)
    zero_rate = rates(system, zero, include_matter_force=True)
    zero_motion = episode.signed_proper_motion(_grid_view(system), zero, zero_rate, zero_rate.r)
    assert np.max(np.abs(zero_motion["proper"])) < 1e-12
    _zero_kr, zero_k_perp = _observer_expansions(zero, zero_rate, system)
    assert np.max(np.abs(zero_k_perp)) < 1e-12


def test_coordinate_work_closes_the_window_and_is_not_observer_pressure():
    system, state = regional.band_limited_state(64)
    rate = rates(system, state, include_matter_force=True)
    _coordinate, windows, partition = regional.fixed_window(
        system.length, system.points, episode.WINDOW_CENTERS, episode.WINDOW_HALF_WIDTH
    )
    window = windows[0]
    channels = _coordinate_channels(state, rate, system, window)
    ledger = regional.matter_ledger(system, state)
    terms = regional.proper_balance_terms(system, state, rate, ledger)
    terms = dict(terms)
    terms["rate"] = rate
    geometry = regional.geometric_ledger(system, state, rate)
    row = episode.window_row(system, window, ledger, geometry, terms)
    assert abs(row["shift_transport"] - channels["shift_transport"]) < 1e-12
    assert abs(row["proper_normal_flux"] - channels["proper_normal_flux"]) < 1e-12
    assert abs(row["shift_pressure_cross"] - channels["shift_pressure_cross"]) < 1e-12
    assert abs(row["coordinate_metric_work"] - channels["coordinate_metric_work"]) < 1e-12
    pressure = _pressure_work(state, rate, system, channels["source"])
    assert np.max(np.abs(pressure - terms["proper_pressure_work"])) < 1e-12
    assert abs(row["coordinate_metric_work"] - float(np.sum(window * pressure))) > 1e-2

    step = 1e-6

    def matter_energy(sample):
        source = source_from_columns(system, sample)
        return system.length_density * source["force_L"] + system.shift * source["force_beta"]

    matter_slope = _difference_quotient(
        (
            matter_energy(regional._combine(state, rate, step)),
            matter_energy(regional._combine(state, rate, -step)),
        ),
        step,
    )
    matter_prediction = (
        channels["shift_transport"]
        + channels["proper_normal_flux"]
        + channels["shift_pressure_cross"]
        + channels["coordinate_metric_work"]
    )
    matter_gap = abs(float(np.sum(window * matter_slope)) - matter_prediction)
    assert matter_gap < 1e-8
    pressure_substituted = matter_prediction - channels["coordinate_metric_work"] + float(np.sum(window * pressure))
    assert abs(float(np.sum(window * matter_slope)) - pressure_substituted) > 1e-2

    def total_energy(sample):
        pieces = regional._geometry_pieces(system, sample)["density"]
        return system.dx * pieces + matter_energy(sample)

    total_slope = _difference_quotient(
        (
            total_energy(regional._combine(state, rate, step)),
            total_energy(regional._combine(state, rate, -step)),
        ),
        step,
    )
    geometric_divergence = row["geometric_shift_flux"] + row["geometric_proper_flux"]
    matter_divergence = (
        row["shift_transport"] + row["proper_normal_flux"] + row["shift_pressure_cross"]
    )
    assert abs(float(np.sum(window * total_slope)) - (matter_divergence - geometric_divergence)) < 1e-7
    assert np.max(np.abs(partition - 1.0)) < 1e-12
    flux_sum = 0.0
    for piece in windows:
        piece_channels = _coordinate_channels(state, rate, system, piece)
        flux_sum += piece_channels["proper_normal_flux"]
    assert abs(flux_sum) < 1e-8

    zero = state.copy()
    zero.p_Q = np.zeros_like(state.p_Q)
    zero.p_r = np.zeros_like(state.p_r)
    zero.p_chi = np.zeros_like(state.p_chi)
    zero_rate = rates(system, zero, include_matter_force=True)
    zero_source = source_from_columns(system, zero)
    work_l1 = float(np.sum(np.abs(zero_source["force_Q"] * zero_rate.Q)))
    pressure_l1 = float(np.sum(np.abs(_pressure_work(zero, zero_rate, system, zero_source))))
    assert work_l1 > 1.0
    assert work_l1 > 100.0 * max(pressure_l1, 1e-12)


def test_galerkin_stage_source_binding_projection_and_geometry_link():
    grid, state = _manufactured_galerkin()
    sample, rate, fields = episode.observe(grid, state)
    fine = galerkin.prolong_state(grid, state)
    radius_x = grid.fine.derivative @ fine.r
    lifted_r = galerkin.prolong_geometry(grid, rate.r)
    expected = _proper_velocity(
        lifted_r, fine.r, radius_x, grid.fine.length_density, grid.fine.shift
    )
    assert np.max(np.abs(fields["proper"] - expected)) < 1e-12
    coarse, bundle = galerkin.compose_fine_hamiltonian(grid, state, include_matter_force=True)
    lifted = episode.lifted_rate(grid, coarse, bundle)
    _k_r, k_perp = _observer_expansions(fine, lifted, grid.fine)
    assert np.max(np.abs(fields["K_perp"] - k_perp)) < 1e-12
    assert np.allclose(fields["rho"], bundle["source"]["force_L"] / grid.dx_q)
    assert sample["mean_subtracted"] is False
    assert sample["rho_mean_removed_before_constraint"] is False
    assert sample["frozen_source_used"] is False

    residual = hamilton_constraint(grid.fine, fine) + bundle["source"]["force_L"] / grid.dx_q
    pulled = (grid.ng / float(grid.nq)) * grid.A_g.T @ residual
    defect = residual - grid.A_g @ pulled
    assert abs(float(np.max(np.abs(defect))) - sample["held_out_hamilton_max"]) < 1e-12
    assert abs(float(np.max(np.abs(residual))) - sample["full_hamilton_max"]) < 1e-12
    assert sample["held_out_hamilton_max"] > 1.0

    windows, _partition = episode.windows_for(grid)
    fine_rate = rates(grid.fine, fine, include_matter_force=True)
    algebraic = _coordinate_channels(fine, lifted, grid.fine, windows[0])
    unprojected = _coordinate_channels(fine, fine_rate, grid.fine, windows[0])
    step = 1e-6

    def matter_energy(sample_state):
        source = source_from_columns(grid.fine, sample_state)
        return grid.fine.length_density * source["force_L"] + grid.fine.shift * source["force_beta"]

    lifted_slope = _difference_quotient(
        (
            matter_energy(regional._combine(fine, lifted, step)),
            matter_energy(regional._combine(fine, lifted, -step)),
        ),
        step,
    )
    lifted_prediction = (
        algebraic["shift_transport"]
        + algebraic["proper_normal_flux"]
        + algebraic["shift_pressure_cross"]
        + algebraic["coordinate_metric_work"]
    )
    lifted_gap = abs(float(np.sum(windows[0] * lifted_slope)) - lifted_prediction)
    unprojected_prediction = (
        unprojected["shift_transport"]
        + unprojected["proper_normal_flux"]
        + unprojected["shift_pressure_cross"]
        + unprojected["coordinate_metric_work"]
    )
    unprojected_slope = _difference_quotient(
        (
            matter_energy(regional._combine(fine, fine_rate, step)),
            matter_energy(regional._combine(fine, fine_rate, -step)),
        ),
        step,
    )
    unprojected_gap = abs(float(np.sum(windows[0] * unprojected_slope)) - unprojected_prediction)
    assert lifted_gap > 1e-2
    assert unprojected_gap < 1e-6
    assert abs(sum(sample["window_proper_flux"])) < 1e-8
    assert abs(sum(sample["window_coordinate_work"]) - sample["lifted_fieldwork"]) < 1e-8

    calls = []
    original_compose = galerkin.compose_fine_hamiltonian
    original_rates = galerkin.rates

    def wrapped(grid_argument, state_argument, include_matter_force=True):
        assert include_matter_force is True
        coarse_rate, bundle_now = original_compose(grid_argument, state_argument, include_matter_force)
        calls.append((
            float(np.max(np.abs(bundle_now["source"]["force_L"]))),
            float(np.max(np.abs(state_argument.r))),
            float(np.max(np.abs(
                bundle_now["fine_state"].r - galerkin.prolong_state(grid_argument, state_argument).r
            ))),
        ))
        return coarse_rate, bundle_now

    galerkin.compose_fine_hamiltonian = wrapped
    galerkin.rates = lambda grid_argument, state_argument, include_matter_force=True: wrapped(
        grid_argument, state_argument, include_matter_force
    )[0]
    try:
        galerkin.rk4_step(grid, state, 1e-4)
    finally:
        galerkin.compose_fine_hamiltonian = original_compose
        galerkin.rates = original_rates
    assert len(calls) == 4
    assert len({round(item[0], 10) for item in calls}) == 4
    assert max(item[1] for item in calls) > min(item[1] for item in calls)
    assert max(item[2] for item in calls) < 1e-12

    evolved = episode.evolve_episode(
        grid, state.copy(), 1e-4, 1e-4, lambda: 0.0, budget=30.0, label="independent"
    )
    assert evolved["completed"] is True
    assert evolved["stop_reason"] is None
    first, second = evolved["samples"]
    trap = 0.5 * 1e-4 * (first["lifted_fieldwork"] + second["lifted_fieldwork"])
    assert abs(second["coordinate_work_integral"] - trap) < 1e-12
    proper_trap = 0.5 * 1e-4 * (first["proper_pressure_power"] + second["proper_pressure_power"])
    assert abs(second["proper_pressure_work_integral"] - proper_trap) < 1e-12
    assert abs(second["coordinate_work_integral"] - second["proper_pressure_work_integral"]) > 1e-6
    decision = episode.initial_tolerance_decision(second["full_hamilton_max"], second["full_momentum_max"])
    assert decision["historical_initial_1e-8_passed"] is False
    assert decision["abort_evolution"] is False
    assert decision["continue_evolution"] is True
    assert evolved["steps_completed"] == 1

    link = _dirac_link(fine, bundle["source"])
    varied = state.copy()
    varied.Q = state.Q * 1.2
    _varied_rate, varied_bundle = galerkin.compose_fine_hamiltonian(grid, varied, include_matter_force=True)
    varied_link = _dirac_link(varied_bundle["fine_state"], varied_bundle["source"])
    assert np.linalg.norm(link - FROZEN_B) > 1.0
    assert np.linalg.norm(varied_link - link) > 0.1
    assert "frozen_B_embedding_claimed" in episode.blank_record()
    assert episode.blank_record()["frozen_B_embedding_claimed"] is False
    assert all("link" not in name and not name.startswith("B_") for name, _key, _floor in episode.PHYSICAL_EFFECTS)

    saved_grid, _saved_state, metadata = episode.load_v5_state("nf512")
    expected_binding = _live_binding(saved_grid)
    saved_binding = {
        str(name): float(value)
        for name, value in zip(metadata["binding_names"], metadata["binding_values"])
    }
    assert saved_binding.keys() == expected_binding.keys()
    for name, value in expected_binding.items():
        assert saved_binding[name] == value
    assert metadata["occupation_gap"] == 0.0
    assert metadata["quadrature_radius_gap"] == 0.0
    assert saved_grid.fine.coefficients["C_W"] != 0.0


def _dirac_link(fine, source):
    matrix = fine.phi0.conj().T @ source["image0"] + fine.phi1.conj().T @ source["image1"]
    return matrix[0:2, 2:4]


def _live_binding(grid):
    binding = {
        "coefficient:" + name: float(grid.fine.coefficients[name])
        for name in ("A", "C_W", "C_F", "C_E", "C_box", "flux", "V_rel")
    }
    for name, value in grid.fine.calibration.items():
        binding["calibration:" + name] = float(value)
    return binding


def test_recorded_effect_scales_stay_per_observable():
    proper = _row(
        "proper_velocity_max", "proper_max", 1e-8, PROPER_MAX_CHANGE, PROPER_MAX_SPACE_MOVEMENT
    )
    held_out = _row(
        "held_out_complement_hamilton", "held_out_hamilton_max", 1e-8,
        HELD_OUT_CHANGE, HELD_OUT_TIME_MOVEMENT,
    )
    energy = _row(
        "coordinate_energy", "energy", 1e-8, ENERGY_CHANGE, ENERGY_TIME_MOVEMENT,
    )
    noisy_mean = _row("proper_velocity_mean", "proper_mean", 1e-8, 2e-8, 6e-8)
    assert proper["status"] == "resolved"
    assert proper["movement"] <= 0.01 * proper["effect_scale"]
    assert held_out["status"] == "unresolved"
    assert held_out["movement_over_effect"] > 0.01
    assert held_out["movement"] < 0.01 * proper["effect_scale"]
    assert energy["status"] == "unresolved"
    assert energy["movement"] < 0.01 * proper["effect_scale"]
    assert noisy_mean["status"] == "unresolved"
    assert proper["status"] == "resolved"


def test_diagnostic_discrepancy_does_not_nullify_a_resolved_effect_without_sensitivity():
    proper = _row(
        "proper_velocity_max", "proper_max", 1e-8, PROPER_MAX_CHANGE, PROPER_MAX_SPACE_MOVEMENT
    )
    held_out = _row(
        "held_out_complement_hamilton", "held_out_hamilton_max", 1e-8,
        HELD_OUT_CHANGE, HELD_OUT_TIME_MOVEMENT,
    )
    assert proper["status"] == "resolved"
    assert held_out["status"] == "unresolved"
    assert held_out["movement"] < 0.01 * proper["effect_scale"]
    verdict = episode.verdict_from("resolved", "unresolved", True, True, False, False, False)
    reason = _sensitivity_reason(verdict)
    assert not _blanket_physical_failure(verdict), (
        "VERDICT_DIAGNOSTIC_VETO: diagnostic_summary 'unresolved' produced "
        f"{verdict!r} for a resolved physical effect whose scale is {proper['effect_scale']}. "
        f"The held-out movement {held_out['movement']} is below one percent of that effect "
        "and no derived sensitivity was attached."
    )
    if reason is not None:
        assert reason["blocks"] is True
        assert reason["diagnostic_movement"] > reason["threshold"]


def test_noisy_row_does_not_erase_a_stable_physical_effect_without_sensitivity():
    proper = _row(
        "proper_velocity_max", "proper_max", 1e-8, PROPER_MAX_CHANGE, PROPER_MAX_SPACE_MOVEMENT
    )
    noisy_mean = _row("proper_velocity_mean", "proper_mean", 1e-8, 2e-8, 6e-8)
    assert proper["status"] == "resolved"
    assert noisy_mean["status"] == "unresolved"
    assert noisy_mean["movement"] < 0.01 * proper["effect_scale"]
    collapsed = episode.summarize_statuses([proper["status"], noisy_mean["status"]])
    verdict = episode.verdict_from(collapsed, "resolved", True, True, False, False, False)
    assert not _blanket_physical_failure(verdict), (
        "VERDICT_NOISY_ROW_VETO: summarize_statuses collapsed a resolved effect and a noisy "
        f"mean into {collapsed!r}, and verdict_from returned {verdict!r} without a derived "
        "sensitivity. The stable row stays resolved on its own one-percent test."
    )
