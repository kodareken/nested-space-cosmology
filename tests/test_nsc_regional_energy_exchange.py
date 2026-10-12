"""Regional ledger controls. These do not accept a continuum limit or a renewal."""
import inspect

import numpy as np

from recursive_horizons import nsc_regional_energy_exchange as regional
from recursive_horizons.nsc_spherical_coupling import field_energy


def test_module_does_not_import_a_frozen_link_current():
    source = inspect.getsource(regional)
    assert "finite_window" not in source
    assert "frozen_onsite" not in source
    assert "1/7" not in source
    assert regional.OLD_CIRCULATION_TRANSFERRED is False
    assert regional.FROZEN_B_EMBEDDING_CLAIMED is False
    assert regional.CONTINUUM_LIMIT_CLAIMED is False
    assert regional.RENEWAL is False
    assert regional.VACUUM_BRANCH_INCLUDED is False


def test_band_limited_channels_close_and_are_independently_visible():
    report = regional.smooth_controls(64)
    assert abs(report["matter_closure_error"]) < regional.TOL_SMOOTH_WINDOW
    assert abs(report["transported_closure_error"]) < regional.TOL_SMOOTH_WINDOW
    assert abs(report["total_closure_error"]) < regional.TOL_SMOOTH_WINDOW
    assert abs(report["global_sum_error"]) < regional.TOL_SUM
    assert abs(report["partition_energy_gap"]) < regional.TOL_SUM
    assert abs(report["field_sum_error"]) < regional.TOL_SUM
    assert abs(report["gravity_sum_error"]) < regional.TOL_SUM
    assert abs(report["quasilocal_gap"]) < regional.TOL_SMOOTH_WINDOW
    assert report["proper_closure_max"] < regional.TOL_PROPER
    assert report["shell_gap_max"] < regional.TOL_SUM
    assert report["ward_max"] < regional.TOL_WARD
    assert report["kernel_gap_max"] < regional.TOL_SUM
    for name in ("shift_pressure_cross", "shift_transport", "coordinate_metric_work"):
        assert report["omitted_channel_error"][name] > 0.5 * abs(report["channels"][name])
        assert report["omitted_channel_error"][name] > 1e-3
    assert report["omitted_channel_error"]["shift_piece_of_Qdot"] > 1e-3
    zero = report["zero_expansion"]
    assert zero["K_r_max"] < 1e-5 and zero["K_perp_max"] < 1e-8
    assert zero["coordinate_work_l1"] > 100.0 * max(zero["proper_pressure_work_l1"], 1e-12)
    # The smooth state is not on the constraint surface, so the regional
    # measured energy is not the constraint smearing.
    assert abs(report["measured_energy"] - report["constraint_smearing"]) > 1e-3
    assert abs(report["matter_energy"]) > 1.0
    assert abs(report["gravity_energy"]) > 1.0


def test_zero_shift_removes_shift_and_mode_channels():
    system, state = regional.band_limited_state(48)
    system.shift = np.zeros(system.points)
    rate = regional.rates(system, state, include_matter_force=True)
    window = 0.5 * (1.0 + np.sin(2 * np.pi * system.xi / system.length))
    _ledger, channels = regional.matter_window_prediction(system, state, rate, window)
    assert abs(channels["shift_transport"]) < 1e-12
    assert abs(channels["shift_pressure_cross"]) < 1e-12
    slope = regional.matter_window_slope(system, state, rate, window)
    assert abs(slope - channels["prediction"]) < regional.TOL_SMOOTH_WINDOW
    assert abs(channels["proper_normal_flux"] + channels["coordinate_metric_work"]) > 1e-3


def test_radius_at_fixed_conformal_data_does_not_change_matter_energy():
    system, state = regional.band_limited_state(32)
    before = field_energy(system, state)
    state.r = state.r * 1.2
    after = field_energy(system, state)
    assert abs(after - before) < 1e-12
    ledger = regional.matter_ledger(system, state)
    # Chain rule at fixed canonical columns: N F_N + q F_q + r F_r = 0.
    assert float(np.max(np.abs(ledger["hamiltonian_chain_ward"]))) < regional.TOL_WARD
    stresses = ledger["stresses"]
    force_l = ledger["source"]["force_L"]
    assert np.max(np.abs(stresses["normal_energy_nodal"] - force_l / state.r)) < 1e-12
    assert np.max(np.abs(stresses["shell_gap"])) < regional.TOL_SUM
    # L*FL is lapse-weighted shell energy, larger than the shell by the lapse factor.
    assert np.max(np.abs(
        ledger["lapse_weighted_integrand"] - (state.r * system.length_density) * stresses["normal_energy_nodal"]
    )) < 1e-8


def test_packet_alias_leakage_is_explicit_and_falls_under_refinement():
    coarse = regional.packet_alias_residual(32)
    fine = regional.packet_alias_residual(64)
    assert coarse["max_abs"] > 1.0
    assert fine["max_abs"] < 0.5 * coarse["max_abs"]


def test_galerkin_lift_is_extra_exchange_beyond_the_fine_identity():
    report = regional.galerkin_leakage(32, 128)
    assert abs(report["unprojected_closure_error"]) < regional.TOL_UNPROJECTED
    assert report["proper_unprojected_max"] < regional.TOL_PROPER
    assert report["proper_lifted_max"] > 10.0 * report["proper_unprojected_max"]
    assert abs(report["projection_leakage"]) > 10 * max(abs(report["unprojected_closure_error"]), 1e-12)
    assert report["frozen_link_copied"] is False
    assert report["uses_manufactured_columns"] is True
