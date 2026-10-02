"""Native period, AP translation, own constraints and bounded extent interface."""
import os
for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(name, "1")
from concurrent.futures import Future
import json
import numpy as np
import pytest
from derive_nsc_discovery_extent import main as cli_main
from recursive_horizons import nsc_discovery_extent as extent


@pytest.fixture(scope="module")
def pairs():
    return {L: extent.build_pair(L, int(8*L)) for L in (8., 12., 16.)}


def test_native_variable_period_operators_and_nominal_resolution(pairs):
    for L, pair in pairs.items():
        grid = pair.grid
        assert grid.nf == int(8*L) and grid.ng == grid.nf-1 and grid.nq == 4*grid.nf
        assert grid.dx_f == pytest.approx(1/8)
        assert grid.dx_q == pytest.approx(1/32)
        k = 2*np.pi/L
        field = np.cos(k*grid.xi_q)
        assert np.max(np.abs(grid.derivative@field + k*np.sin(k*grid.xi_q))) < 1e-10
        ap = np.exp(1j*np.pi*grid.xi_q/L)
        assert np.max(np.abs(grid.momentum@ap - np.pi/L*ap)) < 1e-10
        assert grid.isometry_columns < 1e-10 and grid.isometry_geometry < 1e-10
        assert pair.parent_interval == (2, 6) and pair.child_interval == (3, 5)
        assert pair.clock_locations == (3, 4, 5)
    assert pairs[8.].grid.dx_g != pairs[12.].grid.dx_g
    config = extent.settings()
    assert [(L, int(L*config["points_per_unit"])) for L in config["periods"]] == [(8, 256), (12, 384)]


def test_translation_is_exact_ap_and_source_force_covariance_is_audited(pairs):
    pair = pairs[8.]
    assert extent.translation_audit(pair)["passed"]
    assert extent.translation_audit(pair)["source_carrier_k"] == extent.coupling.CARRIER_K
    for L, current in pairs.items():
        shifted = extent.translation.antiperiodic_shift(current.source_columns[:current.grid.nf], L, L)
        assert np.max(np.abs(shifted+current.source_phi0)) < 1e-12
        assert current.source_metadata["source_support_closures"] == [[2, 3], [3, 5], [5, 6]]
        assert current.source_metadata["physical_widths"] == [1, 2, 1]
        assert current.source_metadata["source_dilated"] is False
        assert np.max(np.abs(current.geometry_map.T@current.geometry_map-np.eye(current.grid.ng))) < 1e-9
    # The matched fermion lattice preserves the source's local sampled shape.
    for name in ("source_phi0", "source_phi1"):
        assert np.max(np.abs(getattr(pairs[8.], name)[:48]-getattr(pairs[12.], name)[:48])) < 1e-12


@pytest.fixture(scope="module")
def solved(pairs):
    results = {}
    for L in (8., 12.):
        state, report = extent.family.solve_prepared(pairs[L])
        results[L] = (state, report)
    return results


def test_each_extent_has_own_dense_source_solve_and_protected_mismatch(pairs, solved):
    samples = {}
    for L, (state, report) in solved.items():
        assert report["converged"]
        assert report["imposed_radius"] is False
        assert report["filtered_v1_radius"] is False
        assert report["shift_residual_max"] >= abs(report["shift_residual_mean"])
        assert np.array_equal(state.phi0, pairs[L].source_phi0)
        assert not np.shares_memory(state.phi0, pairs[L].source_phi0)
        row, samples[L] = extent.observe(pairs[L], state, 0., np.zeros(3), protected=True)
        assert row["length"] == L
        assert row["child_interval"] == [3, 5]
        assert row["tagged_source_columns"] == [2, 3]
        assert row["actual_tides_finite"]
        assert row["chi_substituted_for_curvature"] is False
        assert row["return_is_holding"] is False
        assert row["current_mean_deleted"] is False
        assert row["CAR_min"] >= -1e-8 and row["CAR_max"] <= 1+1e-8
    matching = extent.protected_matching(samples[8.], samples[12.])
    assert set(("r", "N", "q", "r_x", "Q_x", "r_xx", "Q_xx", "R4", "rho")).issubset(matching["gaps"])
    assert matching["gaps"]["r"]["absolute_max"] > 1e-6
    assert matching["coupled_pure_causal_attribution"] == "ambiguous_due_to_changed_initial_local_data"
    assert matching["not_a_continuation_veto"]


def test_checkpoint_preserves_period_not_nf_only(tmp_path, pairs, solved):
    pair = pairs[12.]; state = solved[12.][0]
    arrays = extent.arrays_from(pair, state)
    record = dict(extent.pair_record(pair), case_id="L12_test", input_sha256="test", time=0., steps=0, mode="coupled")
    extent.commit(tmp_path/"checkpoint", record, arrays)
    loaded, saved = extent.load(tmp_path/"checkpoint", "L12_test")
    restored = extent.pair_from_arrays(saved, loaded)
    assert restored.grid.length == 12 and restored.grid.nf == 96 and restored.grid.dx_f == pytest.approx(1/8)
    assert np.array_equal(restored.source_columns, pair.source_columns)
    with pytest.raises(ValueError, match="hashed length"):
        extent.pair_from_arrays(saved, dict(loaded, length=8.))


def test_step_restriction_and_short_adaptation_use_owned_rk4(pairs, solved):
    pair = pairs[12.]; state = solved[12.][0]
    dt, info = extent.episode.step_restriction(pair, state, 1.)
    assert info["geometric_restriction_applied"]
    assert dt <= extent.episode.RK4_HALF_STABILITY/info["geometric"]["omega"]
    fft_pair, fft_state, _ = extent.episode.resolve_pair(pair, state, backend="fft")
    expected = extent.episode.evolving_step(fft_pair, fft_state, 0.002)
    result = extent.advance(pair, state, current_time=0., target=0.002, step_cap=0.002,
                            cadence=0.1, clocks=np.zeros(3), mode="coupled", cpu_allowance=100.)
    assert result["status"] == "target_reached"
    for name in extent.FIELDS:
        assert np.array_equal(getattr(result["state"], name), getattr(expected, name))
    assert result["normal_clocks"][1] > 0
    frozen = extent.advance(pair, state, current_time=0., target=0.002, step_cap=0.002,
                            cadence=0.1, clocks=np.zeros(3), mode="frozen_geometry", cpu_allowance=100.)
    for name in extent.model.GEOMETRY_NAMES+extent.model.MOMENTUM_NAMES:
        assert np.array_equal(getattr(frozen["state"], name), getattr(state, name))


def test_immutable_short_prepare_prefix_and_one_pool(tmp_path):
    directory = tmp_path/"extent-preview"
    prepared = extent.prepare(directory, periods=(8, 12), points_per_unit=8,
                               prefix=0.003, stations=(0.006, 0.009), step_cap=0.003, cadence=0.003)
    assert len(prepared["cases"]) == 2 and not prepared["failures"]
    prefix = extent.materialize(directory)
    assert len(prefix["cases"]) == 4
    for L in (8, 12):
        coupled, a = extent.load(directory, f"L{L}_coupled", 0)
        frozen, b = extent.load(directory, f"L{L}_frozen_geometry", 0)
        assert coupled["source_pins"] == frozen["source_pins"]
        for key in extent.FIELDS+("normal_clocks",):
            assert np.array_equal(a[key], b[key])
        assert coupled["own_source_handoff"] and frozen["own_source_handoff"]
        assert coupled["initial_state_called_during_continuation"] is False
        assert coupled["observations"][0]["mode"] == "coupled"
        assert frozen["observations"][0]["mode"] == "frozen_geometry"
        assert coupled["observations"][0]["worldline_actual_tides"] != frozen["observations"][0]["worldline_actual_tides"]
    class Pool:
        def __init__(self, max_workers):
            assert max_workers == 4
        def __enter__(self):
            return self
        def __exit__(self, *args):
            return False
        def submit(self, fn, payload):
            future = Future()
            future.set_result(fn(payload))
            return future
    result = extent.run(directory, executor=Pool)
    assert result["executor_pools"] == 1 and result["workers"] <= 6
    assert all(case["status"] == "target_reached" and case["time"] <= 0.01 for case in result["results"])
    assert result["echo_proof"] is False and result["return_alone_is_holding"] is False
    assert result["assessment"]["same_extent_feedback_comparisons"]
    assert result["assessment"]["extent_discriminates"] is None
    assert extent.check(directory)["ok"]
    with pytest.raises(FileExistsError):
        extent.materialize(directory)
    with pytest.raises(FileExistsError):
        extent.run(directory, executor=Pool)


def test_preview_production_boundary_and_period16_admission(tmp_path, capsys):
    assert cli_main([]) == 0
    assert "preview" in capsys.readouterr().out
    assert not list(tmp_path.iterdir())
    with pytest.raises(PermissionError, match="period16"):
        extent.settings(periods=(16,))
    confirmed = extent.settings(periods=(16,), period16_admission={"extent_discriminates": True})
    assert confirmed["periods"] == [16.]
    with pytest.raises(PermissionError):
        extent.advance(None, None, current_time=0., target=0.3, step_cap=.001, cadence=.1,
                       clocks=np.zeros(3), mode="coupled", cpu_allowance=100.)


def test_budget_stop_commits_unchanged_last_admissible_state(tmp_path, pairs, solved):
    pair = pairs[12.]; state = solved[12.][0]
    config = extent.settings(periods=(12,), points_per_unit=8, prefix=.003, stations=(.006,),
                             step_cap=.003, cadence=.003)
    directory = tmp_path/"budget-checkpoint"
    record = dict(extent.pair_record(pair), case_id="L12_coupled", input_sha256=config["inputs_sha256"],
                  time=.003, steps=0, mode="coupled", status="handoff", observations=[])
    original = extent.commit(directory, record, extent.arrays_from(pair, state))
    result = extent.execute_case({"directory": str(directory), "case_id": "L12_coupled", "inputs": config,
                                  "cpu_allowance": 0., "production": False, "max_steps": None})
    assert result["status"] == "budget_stop"
    stopped, arrays = extent.load(directory, "L12_coupled")
    assert stopped["ordinal"] == 1
    assert stopped["source_pins"] == original["source_pins"]
    assert stopped["time"] == .003
    assert stopped["steps"] == 0
    assert np.array_equal(arrays["r"], state.r)
    assert np.array_equal(arrays["phi0"], state.phi0)
    assert stopped["cumulative_child_cpu_seconds"] > 0


def test_readonly_legacy_audit_compares_actual_controller_rates_and_time(tmp_path, pairs):
    old = extent.model.build_pair(64, quadrature=256)
    state = extent.model.encode_state(old, extent.galerkin.blank_state(old.grid, old.source_phi0, old.source_phi1))
    moved, shifted = extent.translation.translate_pair_state(old, state, 2.)
    candidate = extent.model.encode_state(pairs[8.], extent.model.reconstruct_state(moved, shifted))
    arrays = extent.episode.arrays_from_state(state, basis_arrays={"W": old.geometry_map,
        "source_phi0": old.source_phi0, "source_phi1": old.source_phi1,
        "observer_columns": old.reference_columns, "source_weights": old.weights},
        clocks={"rates": [1., 1., 1.], "normal_clocks": [0., 0., 0.]})
    record = dict(extent.pair_record(old), case_id="legacy_frozen", nf=64, coordinate_time=.003,
                  control_mode="frozen_geometry", momentum_representation=extent.episode.CANONICAL_PI)
    directory = tmp_path/"legacy"
    committed = extent.episode.commit_checkpoint(directory, record, arrays)
    before = extent.episode.file_sha256(directory/committed["npz"])
    result = extent.audit_existing_period8(directory, "legacy_frozen", pairs[8.], candidate, candidate_time=.003)
    assert result["admissible"]
    assert result["legacy_mode"] == result["candidate_mode"] == "frozen_geometry"
    assert max(result["physical_rate_gaps"].values()) < 1e-8
    assert result["regional_geometric_means"]["child_proper_length"]["absolute_gap"] < 1e-8
    wrong_time = extent.audit_existing_period8(directory, "legacy_frozen", pairs[8.], candidate, candidate_time=.004)
    assert not wrong_time["admissible"]
    assert extent.episode.file_sha256(directory/committed["npz"]) == before
