"""Geometry-graded window: translation, leakage, and record replay."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from recursive_horizons.nsc_geometry_graded_window import (
    CPU_LIMIT_S,
    RECORD_PATH,
    _anticommutator_closed,
    _spinor_image,
    inspect_lowdin,
    translation_algebra,
    verify_saved,
)
from recursive_horizons.nsc_spherical_coupling import CALIBRATION, KAPPA, OMEGA

_MODULE = Path(__file__).resolve().parents[1] / "src" / "recursive_horizons" / "nsc_geometry_graded_window.py"


def _matrix(record):
    return np.array(record["real"], dtype=float) + 1j * np.array(record["imag"], dtype=float)


def test_translation_quadrature_and_direct_operator_agree():
    report = translation_algebra()
    assert report["quadrature_relative"] < 1e-10
    assert report["direct_operator_relative"] < 1e-8
    assert report["direct_versus_closed_relative"] < 1e-8

    omega = float(OMEGA)
    log_omega = float(np.log(omega))
    b0 = float(CALIBRATION["b0"])
    beta0 = float(CALIBRATION["beta0"])
    q0 = b0 / float(CALIBRATION["a0"])
    nodes, weights = leggauss(64)
    left, right = 0.2, 0.9
    x = left + 0.5 * (right - left) * (nodes + 1.0)
    w = 0.5 * (right - left) * weights
    center, width = 0.55, 0.07
    z = (x - center) / width
    upper = np.exp(-(z ** 2))
    lower = 0.15 * (x - center) * upper
    d_upper = upper * (-2.0 * z / width)
    d_lower = 0.15 * upper + (x - center) * 0.15 * d_upper
    derivative_of = {id(upper): d_upper, id(lower): d_lower}

    def image_at(coordinate):
        length = b0 * omega ** coordinate
        radial = np.full(coordinate.shape, q0)
        shift = beta0 * omega ** coordinate

        def anti(values, samples):
            return _anticommutator_closed(values, samples, derivative_of[id(samples)], log_omega)

        return _spinor_image(length, radial, shift, KAPPA, upper, lower, anti)

    here = image_at(x)
    moved = image_at(x + 1.0)
    residual = np.sum(w * (np.abs(moved[0] - omega * here[0]) ** 2 + np.abs(moved[1] - omega * here[1]) ** 2))
    norm = np.sum(w * (np.abs(omega * here[0]) ** 2 + np.abs(omega * here[1]) ** 2))
    assert residual / norm < 1e-18


def test_lowdin_support_and_parity_are_measured():
    report = inspect_lowdin(256)
    assert report["problems"] == []
    assert report["phase_on_minus_spinor_column"] is True
    assert report["phase_applied_to_odd_lobe"] is False
    assert report["gram_offdiag_max"] < 1e-12
    assert report["lowdin_deviation_max"] < 1e-12
    assert report["bridge_packet_mass"] < 1e-18
    assert max(report["support_outside_naive"]) < 1e-18
    assert abs(report["parity_even_odd_overlap"]) < 1e-14
    assert report["reconstructed_envelope_max"] < 1e-10


def test_record_separates_leakage_corner_and_refinement():
    saved = json.loads(RECORD_PATH.read_text())
    comparison = saved["comparison"]
    leak = comparison["projection_leakage_max_nf512"]
    discrete = comparison["discrete_max_abs_nf256_vs_nf512"]
    corner = comparison["corner_frobenius_nf512"]
    link = comparison["link_frobenius_nf512"]
    assert leak > 0.5
    assert discrete < 1e-6
    assert corner < 1e-10
    assert link > 1.0
    assert discrete < leak
    assert corner < discrete or corner < 1e-12
    assert saved["closed_six_mode"] is False
    assert saved["nf512"]["direct_operator_gap"] < 1e-10
    assert saved["nf256"]["direct_operator_gap"] < 1e-10
    assert saved["continuum"]["versus_nf512_max"] < 1e-6


def test_normalized_B_and_block_coefficient_are_different():
    saved = json.loads(RECORD_PATH.read_text())
    omega = saved["omega"]
    local = _matrix(saved["derived"]["H"])
    link = _matrix(saved["derived"]["B01"])
    normalized = _matrix(saved["derived"]["normalized_B"])
    assert np.allclose(normalized, link / np.sqrt(omega))
    assert saved["algebra"]["assembly_max"] < 1e-9
    assert saved["algebra"]["schur_matches_resolvent"] < 1e-8
    assert saved["algebra"]["gamma_matches_schur"] < 1e-8
    assert saved["algebra"]["omitted_omega_factor"] > 1e-3
    assert saved["algebra"]["normalized_B_substituted_as_block_coefficient"] > 1e-3
    assert saved["algebra"]["reversed_link_on_terminal_block"] > 1e-3
    assert saved["algebra"]["congruence_link_max"] < 1e-9
    example = saved["historical_rational_example"]
    assert example["fitted"] is False
    assert example["example_checks_pass"] is True
    assert example["H_gap"] < 1e-8
    assert example["B_gap"] > 1.0


def test_saved_final_Q_breaks_kinetic_grading_and_keeps_the_action():
    saved = json.loads(RECORD_PATH.read_text())
    final = saved["feedback_T_0_05"]
    assert final["evolution_rerun"] is False
    assert abs(final["attained_T"] - 0.05) < 1e-12
    assert final["saved_episode_renewal"] is False
    initial_packets = final["initial_packets_on_final_Q"]
    assert initial_packets["kinetic_diagonal_scale_1"] > 1e-3
    assert initial_packets["mass_diagonal_scale_1"] < 1e-8
    assert initial_packets["shift_diagonal_scale_1"] < 1e-8
    assert saved["nf512"]["grading"]["diagonal_scale_1"] < 1e-9
    action = final["action"]
    assert action["owner"] == "nsc_spherical_feedback_action.first_order_density"
    assert action["functional_edited"] is False
    assert action["density_changed"] is True
    assert final["q0_over_Q_max"] > final["q0_over_Q_min"] + 0.05


def test_scope_flags_and_no_stepper():
    text = _MODULE.read_text()
    assert "rk4_step" not in text
    assert "evolve_episode" not in text
    saved = json.loads(RECORD_PATH.read_text())
    assert saved["schema"] == "NSC-GEOMETRY-GRADED-WINDOW-v1"
    assert saved["status"] == "INITIAL_GRADED_COMPRESSION_NOT_INVARIANT"
    assert saved["new_evolution"] is False
    assert saved["incoming_gate_changed"] is False
    assert saved["infinite_nest_required"] is False
    assert saved["lambdacdm_required"] is False
    assert saved["renewal_established"] is False
    assert saved["universal_physical_fractality"] is False
    assert saved["geometry_grading_is_supplied_input"] is True
    assert saved["cpu_seconds"] <= CPU_LIMIT_S
    assert saved["gauge_nf512"]["length_matches_supplied_gauge"] < 1e-12
    assert saved["gauge_nf512"]["profile_equals_x_on_packet_arc"] < 1e-12


def test_record_replay_is_deterministic():
    report = verify_saved()
    assert report["closed_six_mode"] is False
    assert report["cpu_seconds"] <= CPU_LIMIT_S
    assert report["projection_leakage_max"] > 0.5
    assert report["discrete_max_abs"] < 1e-6


@pytest.mark.parametrize("fermions", [256, 512])
def test_both_resolutions_are_present(fermions):
    saved = json.loads(RECORD_PATH.read_text())
    block = saved[f"nf{fermions}"]
    assert block["projection_leakage_max"] > 0.5
    assert block["grading"]["corner_max"] < 1e-7
    assert block["bridge_image_fraction"] < 1e-8
