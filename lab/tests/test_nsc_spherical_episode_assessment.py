"""Metric and ledger assessment of a stored coupled-transfer episode. No evolution."""
import os

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

from recursive_horizons.nsc_conformal_adm_source import HALF_DENSITY
from recursive_horizons.nsc_spherical_coupling import periodic_derivative
from recursive_horizons.nsc_spherical_episode_assessment import (
    ACCEPTED_JET_ORIGINS,
    CANONICAL_FRAME,
    DRIVER_JOB,
    INDICATOR_KEYS,
    INPUT_SCHEMA,
    OWNERS,
    PROXY_THRESHOLDS_NOT_APPLIED,
    REGENERATION_NPZ,
    _case_episode,
    assess,
    assess_saved_episode,
    packet_geometry,
    direct_rh_grid,
    inspect_saved_episode,
    ricci_scalar_from_christoffel,
    spectral_dx,
    weyl_actual,
    weyl_proxy_from_chi,
)
from recursive_horizons.nsc_spherical_feedback_action import (
    areal_curvature_scalar,
    weyl_square_from_rh,
)

PERIOD = 8.0
ALPHA = 0.2
EDGES = np.array([[0.0, 2.0], [2.0, 4.0], [4.0, 6.0], [6.0, 8.0]])
OCCUPATIONS = np.array([0.05, 0.05, 0.05, 0.05, 0.7, 0.1])


def _bump(x, center, sigma=0.35):
    delta = (x - center + 0.5 * PERIOD) % PERIOD - 0.5 * PERIOD
    return np.exp(-0.5 * (delta / sigma) ** 2)


def _budget(energies, times, pressure, lapse):
    rate = np.gradient(energies, times, axis=0)
    flux = rate - pressure - lapse
    return flux


def make_episode(
    *,
    n=32,
    times=(0.0, 0.5, 1.0),
    centers=None,
    energies=None,
    pressure=None,
    lapse=None,
    chi=0.0,
    occupations=OCCUPATIONS,
    label="primary",
    origin="finite_ode_derivative",
):
    times = np.asarray(times, dtype=float)
    n_t = int(times.size)
    x = np.arange(n, dtype=float) * (PERIOD / n)
    if centers is None:
        centers = np.full(n_t, 1.0)
    centers = np.asarray(centers, dtype=float)
    weight = np.stack([_bump(x, center) for center in centers])
    radial = np.exp(ALPHA * times)[:, None] * np.ones((1, n))
    radial_t = ALPHA * radial
    radial_tt = ALPHA ** 2 * radial
    step = 1e-6
    if energies is None:
        energies = np.tile(np.array([2.0, 0.2, 0.2, 0.2]), (n_t, 1))
    energies = np.asarray(energies, dtype=float)
    if pressure is None:
        pressure = np.zeros_like(energies)
        pressure[:, 0] = 0.3
    if lapse is None:
        lapse = np.zeros_like(energies)
        lapse[:, 0] = 0.2
    pressure = np.asarray(pressure, dtype=float)
    lapse = np.asarray(lapse, dtype=float)
    chi_values = np.full((n_t, n), float(chi)) if np.ndim(chi) == 0 else np.asarray(chi, dtype=float)
    return {
        "driver_job": DRIVER_JOB,
        "binding": {
            "action_owner": OWNERS["action"],
            "geometry_owner": OWNERS["geometry"],
            "observer_owner": OWNERS["observer"],
            "energy_ledger_owner": OWNERS["energy_ledger"],
            "source_id": "coupled-transfer-source",
            "phi_bound": True,
            "driver_job": DRIVER_JOB,
        },
        "chart": {
            "period": PERIOD,
            "x": x,
            "times": times,
            "L": np.ones(n),
            "beta": np.zeros(n),
            "Q": radial,
            "r": np.full(n, 2.0),
            "chi": chi_values,
            "gauge_hold": True,
        },
        "time_jet": {
            "origin": origin,
            "projection_kept": True,
            "projection": "galerkin_stage_projection",
            "Q_dot": radial_t,
            "Q_dot_rate": radial_tt,
            "neighbor": {
                "dt": step,
                "Q_dot_next": radial_t + step * radial_tt,
                "Q_next": radial + step * radial_t,
                "projection_kept": True,
            },
        },
        "ledger": {
            "source_id": "coupled-transfer-source",
            "flux_is_budget_term": True,
            "normal_energy": energies,
            "flux_budget_term": _budget(energies, times, pressure, lapse),
            "pressure_work": pressure,
            "lapse_exchange": lapse,
            "window_edges": EDGES.copy(),
            "packet_weight_kind": "canonical_half_density",
            "packet_weight": weight,
            "phi_basis_occupations": None if occupations is None else np.asarray(occupations, dtype=float),
        },
        "comparison": {"label": label},
    }


def renewed_energies(times):
    times = np.asarray(times, dtype=float)
    energies = np.zeros((times.size, 4))
    energies[:, 0] = 2.0 - 1.5 * times
    energies[:, 1] = 0.2
    energies[:, 2] = 0.2 + 1.5 * times
    energies[:, 3] = 0.2
    pressure = np.zeros_like(energies)
    lapse = np.zeros_like(energies)
    pressure[:, 0] = 0.3
    pressure[:, 2] = 0.3
    lapse[:, 0] = 0.2
    lapse[:, 2] = 0.2
    return energies, pressure, lapse


def _fd_metric_jets(lapse, radial, shift, t0, x0, step=1e-4):
    def block(t, x):
        g_tt = lapse(t, x) ** 2 - radial(t, x) ** 2 * shift(t, x) ** 2
        g_tx = -(radial(t, x) ** 2) * shift(t, x)
        g_xx = -(radial(t, x) ** 2)
        return np.array([[g_tt, g_tx], [g_tx, g_xx]], dtype=float)

    base = block(t0, x0)
    first = np.zeros((2, 2, 2))
    second = np.zeros((2, 2, 2, 2))
    steps = ((step, 0.0), (0.0, step))
    for axis, (dt, dx) in enumerate(steps):
        first[:, :, axis] = (block(t0 + dt, x0 + dx) - block(t0 - dt, x0 - dx)) / (2.0 * step)
    for left, (dt_left, dx_left) in enumerate(steps):
        for right, (dt_right, dx_right) in enumerate(steps):
            if left == right:
                second[left, :, :, right] = (
                    block(t0 + dt_left, x0 + dx_left)
                    - 2.0 * base
                    + block(t0 - dt_left, x0 - dx_left)
                ) / step**2
            else:
                second[left, :, :, right] = (
                    block(t0 + dt_left + dt_right, x0 + dx_left + dx_right)
                    - block(t0 + dt_left - dt_right, x0 + dx_left - dx_right)
                    - block(t0 - dt_left + dt_right, x0 - dx_left + dx_right)
                    + block(t0 - dt_left - dt_right, x0 - dx_left - dx_right)
                ) / (4.0 * step**2)
    return base, first, second


def test_christoffel_matches_direct_formula_on_analytic_charts():
    alpha = 0.2
    time = 0.4
    radial = np.exp(alpha * time)
    metric = np.array([[1.0, 0.0], [0.0, -radial**2]])
    first = np.zeros((2, 2, 2))
    second = np.zeros((2, 2, 2, 2))
    first[1, 1, 0] = -2.0 * alpha * radial**2
    second[0, 1, 1, 0] = -4.0 * alpha**2 * radial**2
    assert np.isclose(ricci_scalar_from_christoffel(metric, first, second), -2.0 * alpha**2)

    flat = np.diag([1.0, -1.0])
    assert ricci_scalar_from_christoffel(flat, np.zeros((2, 2, 2)), np.zeros((2, 2, 2, 2))) == 0.0

    def compare(lapse, radial_fn, shift, expected):
        base, first_jets, second_jets = _fd_metric_jets(lapse, radial_fn, shift, 0.3, 1.7)
        measured = ricci_scalar_from_christoffel(base, first_jets, second_jets)
        assert np.isclose(measured, expected, atol=1e-6, rtol=1e-6)

    compare(lambda t, x: 1.0, lambda t, x: np.exp(alpha * t), lambda t, x: 0.0, -2.0 * alpha**2)
    wave = 2.0 * np.pi / PERIOD
    amplitude = 0.2

    def static_lapse(t, x):
        return 2.0 + amplitude * np.sin(wave * x)

    x0 = 1.7
    lapse_xx = -amplitude * wave**2 * np.sin(wave * x0)
    compare(static_lapse, lambda t, x: 1.0, lambda t, x: 0.0, 2.0 * lapse_xx / static_lapse(0.3, x0))
    shift_amplitude = 0.3

    def static_shift(t, x):
        return shift_amplitude * np.sin(wave * x)

    compare(
        lambda t, x: 1.0,
        lambda t, x: 1.0,
        static_shift,
        -2.0 * shift_amplitude**2 * wave**2 * np.cos(2.0 * wave * x0),
    )

    symbol_t, symbol_x = sp.symbols("t x", real=True)
    owned = areal_curvature_scalar(
        1,
        sp.exp(sp.Rational(1, 5) * symbol_t),
        1,
        0,
        symbol_t,
        symbol_x,
    )
    assert abs(float(owned) + float(2 * sp.Rational(1, 5) ** 2)) < 1e-12
    beta = sp.Rational(3, 10) * sp.sin(2 * sp.pi * symbol_x / 8)
    owned_shift = areal_curvature_scalar(1, 1, 1, beta, symbol_t, symbol_x)
    expected_shift = -2 * sp.Rational(3, 10) ** 2 * (2 * sp.pi / 8) ** 2 * sp.cos(4 * sp.pi * symbol_x / 8)
    for sample in (sp.Rational(1, 2), sp.Rational(3, 2), sp.Integer(0)):
        assert abs(complex((owned_shift - expected_shift).subs(symbol_x, sample))) < 1e-10


def test_spectral_static_shift_and_time_varying_metrics():
    n = 64
    x = np.arange(n) * (PERIOD / n)
    wave = 2.0 * np.pi / PERIOD
    lapse = 2.0 + 0.2 * np.sin(wave * x)
    zeros = np.zeros(n)
    ones = np.ones(n)
    static = direct_rh_grid(
        lapse, ones, zeros, zeros, zeros, PERIOD, L_dot=zeros, beta_dot=zeros,
    )
    expected_static = 2.0 * (-0.2 * wave**2 * np.sin(wave * x)) / lapse
    assert np.max(np.abs(static["R_h"] - expected_static)) < 1e-10

    shift = 0.3 * np.sin(wave * x)
    shifted = direct_rh_grid(
        ones, ones, shift, zeros, zeros, PERIOD, L_dot=zeros, beta_dot=zeros,
    )
    expected_shift = -2.0 * 0.3**2 * wave**2 * np.cos(2.0 * wave * x)
    assert np.max(np.abs(shifted["R_h"] - expected_shift)) < 1e-10

    radial = np.exp(ALPHA * 0.4) * ones
    varying = direct_rh_grid(
        ones, radial, zeros, ALPHA * radial, ALPHA**2 * radial, PERIOD,
        L_dot=zeros, beta_dot=zeros,
    )
    assert np.max(np.abs(varying["R_h"] + 2.0 * ALPHA**2)) < 1e-12
    wrong_zero_rate = direct_rh_grid(
        ones, radial, zeros, ALPHA * radial, np.zeros(n), PERIOD,
        L_dot=zeros, beta_dot=zeros,
    )
    assert np.max(np.abs(wrong_zero_rate["R_h"])) < 1e-12
    assert np.max(np.abs(wrong_zero_rate["R_h"] - varying["R_h"])) > 1e-2

    matrix = periodic_derivative(n, PERIOD)
    samples = np.sin(wave * x)
    assert np.max(np.abs(matrix @ samples - spectral_dx(samples, PERIOD))) < 1e-12


def test_off_shell_chi_and_el_substitution_rejected():
    episode = make_episode(chi=0.0)
    other = make_episode(n=64, chi=0.0, label="refinement")
    result = assess(episode, comparison_episode=other)
    assert result["completed"] is True
    assert result["metric"]["proxy_gap_max"] > 0.05
    assert result["proxy_agreement_is_metric_verification"] is False
    assert result["el_substitution_accepted"] is False
    assert abs(result["metric"]["R_h_mean"] + 2.0 * ALPHA**2) < 1e-12
    proxy = weyl_proxy_from_chi(0.0, 2.0)
    actual = weyl_actual(-2.0 * ALPHA**2, 2.0)
    assert proxy == 0.0
    assert actual > 0.0
    assert np.isclose(float(weyl_square_from_rh(-2.0 * ALPHA**2, 2.0)), float(actual))
    assert np.isclose(float(weyl_square_from_rh(0.0 + 2.0, 2.0)), 0.0)

    on_shell = make_episode(chi=-2.0 - 2.0 * ALPHA**2)
    on_shell_other = make_episode(n=64, chi=-2.0 - 2.0 * ALPHA**2)
    agreed = assess(on_shell, comparison_episode=on_shell_other)
    assert agreed["metric"]["proxy_gap_max"] < 1e-12
    assert agreed["proxy_agreement_is_metric_verification"] is False
    assert agreed["metric"]["metric_from_realized_jet"] is True

    substituted = make_episode()
    substituted["time_jet"]["origin"] = "el_p_chi_dot"
    substituted["time_jet"]["substitute_chi_for_curvature"] = True
    rejected = assess(substituted, comparison_episode=other)
    assert rejected["completed"] is False
    assert rejected["metric_completed"] is False
    assert rejected["reasons"] == ["el_substitution_rejected"]
    assert "metric" not in rejected
    assert "el_p_chi_dot" not in ACCEPTED_JET_ORIGINS


def test_flux_is_not_the_content_rate_and_labels_do_not_move_the_arc():
    episode = make_episode(label="ledger-a")
    other = make_episode(n=64, times=np.linspace(0.0, 1.0, 5), label="ledger-b")
    result = assess(episode, comparison_episode=other)
    assert result["completed"] is True
    assert result["ledger"]["content_rate_end"] == 0.0
    assert result["ledger"]["flux_budget_end"] == -0.5
    assert result["ledger"]["flux_is_content_rate"] is False
    assert result["flux_is_content_rate"] is False
    assert result["ledger"]["budget_residual_max"] < 1e-10
    assert abs(result["ledger"]["work_integral"] - 0.5) < 1e-12
    assert result["maintained_structure"] is True
    assert result["renewed_structure"] is False
    assert result["goal_met"] is True
    assert abs(result["packet"]["arc_end"] - 1.0) < 1e-8
    assert result["global_Q_or_chi_max_used"] is False
    assert result["packet"]["Q_max_diagnostic"] > 1.2

    relabeled = make_episode(label="renamed-input")
    relabeled_other = make_episode(n=64, label="renamed-pair")
    renamed = assess(relabeled, comparison_episode=relabeled_other)
    assert renamed["goal_met"] is result["goal_met"]
    assert renamed["input_label"] == "renamed-input"
    assert abs(renamed["packet"]["arc_end"] - result["packet"]["arc_end"]) < 1e-12

    permuted = make_episode()
    order = [3, 2, 1, 0]
    for key in ("normal_energy", "flux_budget_term", "pressure_work", "lapse_exchange"):
        permuted["ledger"][key] = permuted["ledger"][key][:, order]
    permuted["ledger"]["window_edges"] = permuted["ledger"]["window_edges"][order]
    permuted_other = make_episode(n=64)
    for key in ("normal_energy", "flux_budget_term", "pressure_work", "lapse_exchange"):
        permuted_other["ledger"][key] = permuted_other["ledger"][key][:, order]
    permuted_other["ledger"]["window_edges"] = permuted_other["ledger"]["window_edges"][order]
    moved = assess(permuted, comparison_episode=permuted_other)
    assert moved["goal_met"] is True
    assert abs(moved["packet"]["region_midpoint"] - 1.0) < 1e-12
    assert moved["packet"]["basis_peak"] == 4


def test_basis_occupations_are_not_the_spatial_window():
    uniform = make_episode()
    uniform["ledger"]["packet_weight"] = np.ones_like(uniform["ledger"]["packet_weight"])
    uniform_other = make_episode(n=64)
    uniform_other["ledger"]["packet_weight"] = np.ones_like(uniform_other["ledger"]["packet_weight"])
    spread = assess(uniform, comparison_episode=uniform_other)
    assert spread["completed"] is True
    assert spread["localized"] is False
    assert spread["goal_met"] is False
    assert spread["spatial_arc_uses_basis"] is False
    assert spread["packet"]["basis_peak"] == 4

    peaked = make_episode()
    peaked["ledger"]["phi_basis_occupations"] = OCCUPATIONS[::-1]
    peaked_other = make_episode(n=64)
    peaked_other["ledger"]["phi_basis_occupations"] = OCCUPATIONS[::-1]
    spatial = assess(peaked, comparison_episode=peaked_other)
    assert spatial["localized"] is True
    assert spatial["packet"]["basis_peak"] == 1
    assert abs(spatial["packet"]["arc_end"] - 1.0) < 1e-8
    assert spatial["maintained_structure"] is True


def test_renewal_follows_the_packet_arc_and_two_cuts_do_not():
    times = np.linspace(0.0, 1.0, 3)
    energies, pressure, lapse = renewed_energies(times)
    episode = make_episode(
        times=times,
        centers=1.0 + 4.0 * times,
        energies=energies,
        pressure=pressure,
        lapse=lapse,
    )
    fine_times = np.linspace(0.0, 1.0, 5)
    fine_energies, fine_pressure, fine_lapse = renewed_energies(fine_times)
    other = make_episode(
        n=64,
        times=fine_times,
        centers=1.0 + 4.0 * fine_times,
        energies=fine_energies,
        pressure=fine_pressure,
        lapse=fine_lapse,
    )
    result = assess(episode, comparison_episode=other)
    reverse = assess(other, comparison_episode=episode)
    assert result["completed"] is True
    assert result["renewed_structure"] is True
    assert result["maintained_structure"] is False
    assert result["goal_met"] is True
    assert result["packet"]["arrived_on_initial_bridge"] is True
    assert abs(result["packet"]["arc_start"] - 1.0) < 1e-8
    assert abs(result["packet"]["arc_end"] - 5.0) < 1e-8
    assert result["ledger"]["flux_budget_end"] != result["ledger"]["content_rate_end"]
    assert reverse["renewed_structure"] is True
    assert reverse["goal_met"] is True
    assert result["actual_bound"] is None
    assert result["comparison_uncertainty"]["actual_bound"] is None
    assert result["refinement_indicator"]["arc_location"] < 1e-6
    assert set(result["refinement_indicator"]) == set(INDICATOR_KEYS)

    cut = make_episode(times=np.array([0.0, 1.0]))
    cut_other = make_episode(n=64, times=np.array([0.0, 1.0]))
    split = assess(cut, comparison_episode=cut_other)
    assert split["completed"] is True
    assert split["arbitrary_two_cut"] is True
    assert split["meaningful_transfer_episode"] is False
    assert split["goal_met"] is False
    assert PROXY_THRESHOLDS_NOT_APPLIED["reversal_at_least"] == 0.1
    assert split["programme_thresholds_not_applied"]["leader_share_at_least"] == 0.5


def test_missing_rate_or_source_rejects_without_zero_fill():
    episode = make_episode()
    other = make_episode(n=64)
    del episode["time_jet"]["Q_dot_rate"]
    missing_rate = assess(episode, comparison_episode=other)
    assert missing_rate["completed"] is False
    assert missing_rate["metric_completed"] is False
    assert "missing_Q_dot_rate" in missing_rate["reasons"]
    assert "metric" not in missing_rate

    unbound = make_episode()
    unbound["binding"]["source_id"] = ""
    unbound["ledger"]["source_id"] = ""
    missing_source = assess(unbound, comparison_episode=other)
    assert missing_source["metric_completed"] is False
    assert "missing_source_binding" in missing_source["reasons"]

    mismatched = make_episode()
    mismatched["ledger"]["source_id"] = "other-source"
    mismatch = assess(mismatched, comparison_episode=other)
    assert "source_id_mismatch" in mismatch["reasons"]
    assert mismatch["completed"] is False

    ambiguous = make_episode()
    ambiguous["ledger"]["flux_is_budget_term"] = False
    assert "flux_sign_ambiguous" in assess(ambiguous, comparison_episode=other)["reasons"]

    disagreed = make_episode()
    disagreed["time_jet"]["neighbor"]["Q_dot_next"] = disagreed["time_jet"]["Q_dot"]
    disagreed_result = assess(disagreed, comparison_episode=other)
    assert disagreed_result["reasons"] == ["neighbor_increment_disagrees"]
    assert "metric" not in disagreed_result

    no_comparison = assess(make_episode())
    assert no_comparison["metric_completed"] is True
    assert no_comparison["completed"] is False
    assert "missing_comparison_uncertainty" in no_comparison["reasons"]
    assert abs(no_comparison["metric"]["R_h_mean"] + 2.0 * ALPHA**2) < 1e-12


def test_positive_chart_and_saved_payload():
    broken = make_episode()
    broken["chart"]["Q"] = -np.abs(broken["chart"]["Q"])
    chart = assess(broken, comparison_episode=make_episode(n=64))
    assert chart["completed"] is False
    assert "positive_chart_failed" in chart["reasons"]
    assert chart["metric_completed"] is False

    forced = make_episode()
    forced["added_force"] = True
    assert assess(forced)["reasons"] == ["manufactured_force_rejected"]

    saved = inspect_saved_episode()
    assert saved["present"] is True
    assert saved["completed"] is False
    assert saved["metric_completed"] is False
    assert saved["missing_rates_set_to_zero"] is False
    assert saved["coarse_frame_difference_used_as_rate"] is False
    assert saved["stored_frame_weyl_is_chi_proxy"] is True
    assert "frame Q_dot" in saved["missing_primitive"]
    assert "frame Q_dot_rate" in saved["missing_primitive"]
    assert "neighbor increment of Q_dot" in saved["missing_primitive"]
    assert "frame p_chi" in saved["missing_primitive"]
    assert saved["available"]["endpoint_Q_dot_shape"] == (511,)
    assert saved["available"]["frame_Q_shape"] == (11, 2048)
    assert saved["available"]["window_normal_shape"] == (101, 4)
    assert saved["available"]["endpoint_rate_matches_frame_nodes"] is False
    assert saved["available"]["frame_Q_shape"] is not None
    assert "Q" in saved["available"]["frame_fields"]
    assert "weyl_C2" in saved["available"]["frame_fields"]
    assert saved["available"]["window_proper_flux"] is True
    assert saved["available"]["window_lapse_work"] is True
    assert saved["closed_regime"] is False
    assert saved["actual_bound"] is None
    assert INPUT_SCHEMA["time_jet"]["origin"].startswith("realized_increment")
    assert INPUT_SCHEMA["binding"]["driver_job"] == DRIVER_JOB
    assert "canonical_half_density" in INPUT_SCHEMA["ledger"]["packet_weight_kind"]


def _static_jet(episode, radial):
    radial = np.asarray(radial, dtype=float)
    shape = episode["chart"]["Q"].shape
    episode["chart"]["Q"] = np.broadcast_to(radial, shape).copy()
    episode["time_jet"]["Q_dot"] = np.zeros(shape)
    episode["time_jet"]["Q_dot_rate"] = np.zeros(shape)
    episode["time_jet"]["neighbor"]["Q_dot_next"] = np.zeros(shape)
    episode["time_jet"]["neighbor"]["Q_next"] = episode["chart"]["Q"].copy()
    return episode


def test_canonical_probability_does_not_pick_up_an_extra_q():
    assert HALF_DENSITY == CANONICAL_FRAME
    uniform = make_episode()
    varied = make_episode()
    profile = 1.0 + 0.35 * np.sin(2.0 * np.pi * varied["chart"]["x"] / PERIOD)
    _static_jet(varied, profile)
    fine = make_episode(n=64)
    fine_profile = 1.0 + 0.35 * np.sin(2.0 * np.pi * fine["chart"]["x"] / PERIOD)
    _static_jet(fine, fine_profile)
    uniform_result = assess(uniform, comparison_episode=make_episode(n=64))
    varied_result = assess(varied, comparison_episode=fine)
    assert uniform_result["packet"]["measure_is_mode_probability"] is True
    assert varied_result["packet"]["q_multiplied_into_measure"] is False
    assert varied_result["packet"]["shell_energy_is_mode_probability"] is False
    assert abs(
        uniform_result["packet"]["measure_mass"] - varied_result["packet"]["measure_mass"]
    ) < 1e-12
    assert abs(
        uniform_result["packet"]["occupation_concentration"]
        - varied_result["packet"]["occupation_concentration"]
    ) < 1e-12
    weight = uniform["ledger"]["packet_weight"][0]
    q_uniform = 2.0 * uniform["chart"]["Q"][0]
    q_varied = 2.0 * varied["chart"]["Q"][0]
    assert abs(float(np.sum(weight * q_uniform) - np.sum(weight * q_varied))) > 1e-3
    phase_uniform = packet_geometry(
        uniform["chart"]["x"], PERIOD, weight, q_uniform,
    )["proper_phase"]
    phase_varied = packet_geometry(
        varied["chart"]["x"], PERIOD, weight, q_varied,
    )["proper_phase"]
    assert abs(phase_uniform - phase_varied) > 1e-3

    bare = make_episode()
    del bare["ledger"]["packet_weight_kind"]
    assert "packet_measure_rejected" in assess(bare, comparison_episode=make_episode(n=64))["reasons"]

    shell = make_episode()
    shell["ledger"]["packet_weight_kind"] = "positive_shell_energy"
    shell_other = make_episode(n=64)
    shell_other["ledger"]["packet_weight_kind"] = "positive_shell_energy"
    shell_result = assess(shell, comparison_episode=shell_other)
    assert shell_result["packet"]["measure_kind"] == "positive_shell_energy"
    assert shell_result["packet"]["measure_is_mode_probability"] is False
    assert shell_result["packet"]["q_multiplied_into_measure"] is False


def test_opposite_clumps_stay_localized_when_the_first_harmonic_cancels():
    def clumps(n, times):
        episode = make_episode(n=n, times=times)
        coordinate = episode["chart"]["x"]
        weight = _bump(coordinate, 1.0) + _bump(coordinate, 5.0)
        episode["ledger"]["packet_weight"] = np.stack([weight for _ in times])
        return episode

    times = np.linspace(0.0, 1.0, 3)
    result = assess(clumps(32, times), comparison_episode=clumps(64, times))
    assert result["completed"] is True
    assert result["packet"]["resultant"] < 1e-6
    assert result["packet"]["second_resultant"] > 0.5
    assert result["packet"]["occupation_excess"] > 1.0
    assert result["localized"] is True
    assert result["homogeneous"] is False
    assert result["single_harmonic_resultant_is_localization"] is False
    assert result["packet"]["location_kind"] == "proper_arc_second_harmonic_axis"
    assert result["packet"]["region_count"] == 2
    assert result["packet"]["q_multiplied_into_measure"] is False


def test_a_proven_bound_is_not_the_refinement_indicator():
    episode = make_episode()
    other = make_episode(n=64)
    episode["mathematical_bound"] = {"name": "supplied_domination", "value": 1.0, "proven": True}
    other["mathematical_bound"] = dict(episode["mathematical_bound"])
    result = assess(episode, comparison_episode=other)
    assert result["closed_regime"] is True
    assert result["actual_bound"] == 1.0
    assert result["refinement_indicator"]["occupation_concentration"] != result["actual_bound"]
    plain = assess(make_episode(), comparison_episode=make_episode(n=64))
    assert plain["closed_regime"] is False
    assert plain["actual_bound"] is None
    assert plain["bound_status"] == "not_supplied"


def test_saved_regeneration_episode_is_maintained_throughflow():
    import hashlib

    before = hashlib.sha256(REGENERATION_NPZ.read_bytes()).hexdigest()
    result = assess_saved_episode()
    after = hashlib.sha256(REGENERATION_NPZ.read_bytes()).hexdigest()
    assert before == after
    assert result["completed"] is True
    assert result["bound_status"] == "propagated_initial_bound_null"
    assert result["actual_bound"] is None
    assert result["closed_regime"] is False
    assert result["indicator_is_not_a_curvature_bound"] is True
    assert result["geometric_event"]["entered_packet_edge"] is True
    assert result["geometric_event"]["bridge_arc_x_start"] == 3.96875
    assert result["el_substitution_accepted"] is False
    primary = result["cases"]["nf512_dtmax_0_00025"]
    assert primary["goal_met"] is True
    assert primary["maintained_structure"] is True
    assert primary["renewed_structure"] is False
    assert primary["packet"]["q_multiplied_into_measure"] is False
    assert primary["packet"]["arc_displacement"] < primary["packet"]["location_width"]
    assert primary["ledger"]["budget_residual_max"] < 1e-5
    assert primary["metric"]["proxy_gap_max"] > 1e-3
    assert primary["record"]["el_p_chi_dot_used"] is False
    assert primary["record"]["gram_admissible"] is True
    assert primary["record"]["positive_chart_series"] is True
    for case in result["cases"].values():
        assert case["completed"] is True
        assert case["maintained_structure"] is True
        for mass in case["record"]["fermion_occupation"]:
            assert abs(mass - 3.0) < 1e-8
        for mass in case["record"]["quadrature_occupation"]:
            assert abs(mass - 3.0) < 1e-8
    try:
        _case_episode({}, "nf512_dtmax_0_00025_", {}, np.ones(6), "nsc-regeneration-episode-v1")
    except KeyError as exc:
        assert "frame_" in str(exc)
    else:
        raise AssertionError("a missing jet was accepted")
