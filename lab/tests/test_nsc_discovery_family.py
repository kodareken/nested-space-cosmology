"""Source-family preparation, handoff, and regime class. No production trajectory."""
import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import json
from concurrent.futures import Future
from pathlib import Path

import numpy as np
import pytest

from derive_nsc_discovery_family import main as cli_main
from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_discovery_family as family
from recursive_horizons import nsc_discovery_prediction as prediction
from recursive_horizons import nsc_nested_parent_child as model
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


_SOURCE = Path(family.__file__).read_text()


def test_contract_is_one_pool_of_six_without_a_percent_gate_or_a_new_integrator():
    spec = family.specification()
    assert spec["schema"] == "NSC-DISCOVERY-SOURCE-FAMILY-v1"
    assert spec["exploration_cases"] == 6
    assert spec["threads_per_case"] == 1
    assert spec["workers"] == 6
    assert spec["coordinator_pools"] == 1
    assert spec["aggregate_cpu_budget_seconds"] == 21600.0
    assert spec["chunk_limit_bytes"] == 64 * 1024 ** 2
    assert spec["occupation_trace"] == 3.0
    assert spec["scalar_fate_percent_gate"] is None
    assert spec["universal_one_percent_gate"] is False
    assert spec["forced_regeneration"] is False
    assert spec["saved_t03_radius_copied"] is False
    assert spec["third_dynamical_framework"] is False
    assert spec["constructor"] == "nsc_nested_parent_child.initial_state"
    assert spec["evolution"] == "nsc_discovery_prediction.propagate_baseline"
    assert "ProcessPoolExecutor" not in _SOURCE
    assert "def rk4_step" not in _SOURCE
    assert "propagate_baseline" in _SOURCE


def test_endpoints_uniform_weights_and_width_projection():
    grid = galerkin.build_grid(64, gauge="conformal")
    audited = family.audited_chart(grid)
    assert audited["A"] == coupling.locked_coefficients()["A"]
    assert audited["C_W"] == coupling.locked_coefficients()["C_W"]
    assert audited["flux"] == coupling.locked_coefficients()["flux"]
    assert audited["kappa"] == coupling.KAPPA
    assert np.array_equal(family.occupation_weights(1.0), coupling.OCCUPATIONS)
    owned0, owned1, _owned = model._separated_columns(grid)
    for scale in family.WIDTH_SCALES:
        phi0, phi1, meta = family.owned_separated_columns(grid, scale)
        width = family.half_width(scale)
        assert meta["supports"] == [[0.0, 2.0 - width], [2.0 - width, 2.0 + width], [2.0 + width, 4.0]]
        assert meta["explicit_orthogonalization"] == family.EXPLICIT_ORTHOGONALIZATION
        assert meta["child_included_in_lowdin"] is False
        assert meta["child_frame_altered"] is False
        assert meta["child_frame_silently_altered"] is False
        assert meta["mean_current_deleted"] is False
        assert meta["gram_max"] < 1e-8
        assert meta["overlap_after_max"] < 1e-8
        assert len(meta["tails"]) == 6
        assert meta["tail_max"] >= 0.0
        assert meta["exact_spatial_support"] is False
        assert np.isfinite(meta["mean_current"])
        reference0, reference1, _preparation = coupling.prepare_rank6(grid.xi_f, grid.dx_f)
        if scale == 1.0:
            assert np.array_equal(phi0[:, 2:4], reference0[:, 2:4])
            assert np.array_equal(phi1[:, 2:4], reference1[:, 2:4])
            assert meta["child_frame_origin"] == "owned_original_middle"
            assert meta["child_pair_lowdin_explicit"] is False
            assert np.array_equal(phi0, owned0)
            assert np.array_equal(phi1, owned1)
        else:
            assert meta["child_frame_origin"] == "declared_width_scaled_lobes"
            assert not np.array_equal(phi0[:, 2:4], reference0[:, 2:4])
        for imbalance in family.IMBALANCES:
            weights = family.occupation_weights(imbalance)
            eigenvalues = prediction.assert_admissible_source(phi0, phi1, weights)
            assert abs(float(np.sum(weights)) - 3.0) < 1e-12
            assert np.min(eigenvalues) >= -1e-8
            assert np.max(eigenvalues) <= 1.0 + 1e-8
            assert family.regional_occupation_gap(weights) == pytest.approx(imbalance)


def test_frozen_portable_basis_is_complete():
    for fermions in (128, 256):
        frozen = family.load_frozen_basis(fermions)
        degree = fermions - 1
        assert frozen["W"].shape == (degree, degree)
        assert frozen["complete"] is True
        assert frozen["orthogonality_max"] < 1e-8
        assert frozen["basis"] == "frozen_portable_baseline"


def test_regime_class_uses_physical_signs_and_one_discriminating_gap():
    baseline = family.physical_signals(
        length_change=-0.2, width_change=0.2, retained_fraction=0.9, packet_retention=(0.2, 0.9, 0.4),
    )
    tiny = family.physical_signals(
        length_change=-1.0e-4, width_change=1.0e-4, retained_fraction=0.4, packet_retention=(0.2, 0.4, 0.3),
    )
    huge = family.physical_signals(
        length_change=-0.8, width_change=0.8, retained_fraction=0.95, packet_retention=(0.1, 0.95, 0.2),
    )
    opposite = family.physical_signals(
        length_change=0.2, width_change=-0.2, retained_fraction=0.9, packet_retention=(0.9, 0.2, 0.4),
    )
    assert family.fate_class(tiny, baseline) == "same_physical_signs"
    assert family.fate_class(huge, baseline) == "same_physical_signs"
    assert family.fate_class(opposite, baseline) == "distinct_physical_signs"
    assert family.fate_class(baseline, baseline, is_baseline=True) == "baseline"
    assert family.regime_observable(baseline)["name"] == "child_source_retained_fraction"
    assert family.regime_observable(baseline)["percent_gate"] is None
    assert family.regime_observable(baseline)["value"] == pytest.approx(0.9)
    assert family.SCALAR_FATE_PERCENT_GATE is None
    gap_plus = family.regional_occupation_gap(family.occupation_weights(1.0))
    gap_minus = family.regional_occupation_gap(family.occupation_weights(-1.0))
    gap_zero = family.regional_occupation_gap(family.occupation_weights(0.0))
    assert gap_plus == pytest.approx(1.0)
    assert gap_minus == pytest.approx(-1.0)
    assert gap_zero == pytest.approx(0.0)
    assert gap_plus != gap_minus
    missing = family.physical_signals(
        length_change=-0.2, width_change=None, retained_fraction=0.9, packet_retention=(0.2, 0.9, 0.4),
    )
    assert family.fate_class(missing, baseline) == "unresolved"
    assert family.select_confirmation([
        {"case_id": "baseline", "imbalance": 1.0, "width_scale": 1.0, "fate": "baseline"},
        {"case_id": "open", "imbalance": 0.0, "width_scale": 1.0, "fate": "unresolved"},
        {"case_id": "moved", "imbalance": -1.0, "width_scale": 1.2, "fate": "distinct_physical_signs"},
    ])["cases"] == [{
        "imbalance": -1.0, "width_scale": 1.2, "parent_case": "moved", "nf": 256,
        "fate": "distinct_physical_signs", "forced": False,
    }]


def test_stalled_corner_is_the_owned_homotopy_blocker_and_not_a_copied_radius():
    with pytest.raises(ValueError, match="positive homotopy stalled above measured radius-correction floor"):
        family.prepare_member(-1.0, 1.2, nf=64, basis="owned_geometry_frame")


@pytest.fixture(scope="module")
def prepared_geometries():
    baseline = family.prepare_member(1.0, 1.0, nf=64, basis="owned_geometry_frame")
    weights = family.prepare_member(-1.0, 1.0, nf=64, basis="owned_geometry_frame")
    shape = family.prepare_member(1.0, 0.8, nf=64, basis="owned_geometry_frame")
    return baseline, weights, shape


def test_every_prepared_geometry_solves_its_own_source(prepared_geometries):
    baseline, weights, shape = prepared_geometries
    for item in (baseline, weights, shape):
        report = item["report"]
        state = item["state"]
        pair = item["pair"]
        assert report["converged"] is True
        assert report["imposed_radius"] is False
        assert report["filtered_v1_radius"] is False
        assert report["saved_t03_radius_copied"] is False
        assert "retained" in report["shift_method"]
        assert np.isfinite(report["shift_residual_mean"])
        assert np.isfinite(report["solver_diagnostics"]["full_hamilton_max"])
        assert np.isfinite(report["solver_diagnostics"]["held_out_hamilton_max"])
        assert item["initial"]["audited"]["A"] == coupling.locked_coefficients()["A"]
        assert item["initial"]["occupation_trace"] == pytest.approx(3.0)
        assert np.array_equal(state.phi0, np.asarray(pair.source_phi0))
        assert np.array_equal(state.phi1, np.asarray(pair.source_phi1))
        assert state.phi0 is not pair.source_phi0
        assert not np.shares_memory(state.phi0, np.asarray(pair.source_phi0))
        assert not np.shares_memory(state.phi1, np.asarray(pair.source_phi1))
        assert np.max(np.abs(state.r - 1.0)) > 1e-8
    assert np.array_equal(baseline["state"].phi0, np.asarray(weights["pair"].source_phi0))
    assert np.max(np.abs(baseline["state"].r - weights["state"].r)) > 1e-6
    assert not np.array_equal(baseline["state"].phi0, np.asarray(shape["pair"].source_phi0))
    assert np.max(np.abs(baseline["state"].r - shape["state"].r)) > 1e-6
    assert baseline["initial"]["regional_occupation_gap"] == pytest.approx(1.0)
    assert weights["initial"]["regional_occupation_gap"] == pytest.approx(-1.0)
    assert shape["initial"]["width_scale"] == pytest.approx(0.8)
    measured = family.measure_consumers(
        baseline["pair"], baseline["state"], [0.0, 0.0, 0.0], 0.0,
    )
    if measured["regions_ready"]:
        assert measured["child_retained_fraction"] is not None
    else:
        assert measured["blockers"]
        assert measured["blockers"][0]["consumer"] == "regions"
        assert ":" in measured["blockers"][0]["error"]


def test_continued_output_is_the_propagated_state_not_a_reset(tmp_path, prepared_geometries):
    prepared = prepared_geometries[0]
    propagated, clocks, identity = family.evolve_to_handoff(
        prepared["pair"], prepared["state"], duration=0.005, step_cap=0.005, production=False,
    )
    assert identity["handoff_phi_equals_source"] is False
    assert identity["handoff_is_propagate_baseline_output"] is True
    assert identity["saved_t03_radius_copied"] is False
    assert abs(identity["coordinate_time"] - 0.005) < 1e-8
    evolved = propagated["state"]
    assert not np.array_equal(evolved.phi0, np.asarray(prepared["pair"].source_phi0))
    assert episode.state_sha256(evolved) == identity["state_sha256"]
    directory = tmp_path / "handoff"
    manifest = family.commit_campaign(directory, [{
        "pair": prepared["pair"],
        "initial": prepared["initial"],
        "basis_pins": prepared["basis_pins"],
        "propagated": propagated,
        "clocks": clocks,
        "identity": identity,
        "step_cap": 0.005,
    }])
    assert manifest["executor_pools"] == 0
    assert manifest["chunk_limit_bytes"] == 64 * 1024 ** 2
    assert manifest["initial_cpu_seconds"] > 0.0
    assert manifest["cpu_budget_seconds"] < family.AGGREGATE_CPU_BUDGET_SECONDS
    chunk = directory / manifest["cases"][0]["npz"]
    assert chunk.stat().st_size < 64 * 1024 ** 2
    case_id = prepared["initial"]["case_id"]
    _record, arrays = episode.load_checkpoint(directory, case_id, 0)
    assert np.array_equal(arrays["phi0"], evolved.phi0)
    assert np.array_equal(arrays["r"], evolved.r)
    assert not np.array_equal(arrays["phi0"], arrays["source_phi0"])
    seen = {}

    def capture(pair, current, dt):
        seen["sha"] = episode.state_sha256(current)
        seen["equals_source"] = bool(np.array_equal(current.phi0, np.asarray(pair.source_phi0)))
        return episode.evolving_step(pair, current, dt)

    result = episode.execute_case(
        directory, case_id, cpu_allowance=1.0e9, max_steps=1, stepper=capture, backend="dense",
    )
    assert seen["sha"] == identity["state_sha256"]
    assert seen["equals_source"] is False
    assert result["steps"] == 1
    reloaded, reloaded_arrays = episode.load_checkpoint(directory, case_id, 0)
    assert reloaded["arrays_sha256"] == manifest["cases"][0]["arrays_sha256"]
    assert np.array_equal(reloaded_arrays["phi0"], evolved.phi0)
    _later, later_arrays = episode.load_checkpoint(directory, case_id)
    assert not np.array_equal(later_arrays["phi0"], evolved.phi0)
    assert not np.array_equal(later_arrays["phi0"], arrays["source_phi0"])


def test_one_pool_budget_and_separate_output_directories(tmp_path):
    explore = tmp_path / "explore"
    confirm = tmp_path / "confirm"
    left = family.prepare_specification(explore)
    right = family.prepare_confirmation(confirm, explore)
    assert Path(left["output"]) != Path(right["output"])
    assert len(left["cases"]) == 6
    assert len(left["family"]) == 15
    assert left["forecast"]["threads_per_case"] == 1
    assert left["forecast"]["coordinator_pools"] == 1
    assert left["forecast"]["chunk_limit_bytes"] == 64 * 1024 ** 2
    assert left["forecast"]["aggregate_budget_seconds"] == 21600.0
    assert left["forecast"]["dense_initial_solves"] == 6
    assert left["forecast"]["dense_initial_solves_in_rate_sample"] is False
    assert right["forced_regeneration"] is False
    assert right["selected"]["cases"] == []
    assert right["selected"]["reason"] == "no distinct physical regime to confirm"
    assert not (confirm / "manifest.json").exists()
    directory = tmp_path / "pool"
    directory.mkdir()
    (directory / "manifest.json").write_text(json.dumps({
        "schema": episode.SCHEMA,
        "stage": 0,
        "status": "PREPARED",
        "cases": [{"case_id": f"case-{index}"} for index in range(6)],
        "initial_cpu_seconds": 10.0,
        "chunk_limit_bytes": family.CHUNK_LIMIT_BYTES,
        "evolved": False,
    }))
    (directory / "source-family.json").write_text(json.dumps(episode.jsonable({
        **family.specification(),
        "output": str(directory),
        "materialized": True,
        "cases": [f"case-{index}" for index in range(6)],
        "executor_pools": 0,
    })))

    class Pool:
        def __init__(self, max_workers):
            self.max_workers = max_workers

        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

        def submit(self, fn, payload):
            future = Future()
            future.set_result({
                "case_id": payload["case_id"],
                "status": "budget_stop",
                "child_cpu_seconds": 1.0,
                "coordinate_time": 0.3,
                "steps": 0,
                "stations_reached": [],
                "stop": {"stop_class": "budget_stop"},
                "chunk": None,
                "stability_certificate": False,
                "state_clamped": False,
                "physical_instability_claimed": False,
            })
            return future

    result = family.run(directory, executor=Pool, cpu_budget_seconds=21600.0)
    assert result["executor_pools"] == 1
    assert result["coordinator_pools"] == 1
    assert result["workers"] == 6
    assert result["initial_cpu_seconds"] == pytest.approx(10.0)
    assert result["cpu_budget_seconds"] == pytest.approx(21590.0)
    assert result["aggregate_cpu_seconds"] == pytest.approx(16.0 + result["assessment_cpu_seconds"])
    assert family.remaining_budget(21600.0, 10.0) == pytest.approx(21590.0)
    with pytest.raises(RuntimeError):
        family.remaining_budget(21600.0, 21600.1)
    with pytest.raises(FileNotFoundError):
        family.run(explore)


def test_cli_catalog_prepare_confirm_and_run_guard(tmp_path, capsys):
    assert cli_main(["--catalog"]) == 0
    catalog = capsys.readouterr().out
    assert "members 15" in catalog
    assert "cases 6" in catalog
    assert "trace 3.0" in catalog
    assert "budget 21600.0" in catalog
    assert "chunk 67108864" in catalog
    assert "pools 1" in catalog
    explore = tmp_path / "explore"
    confirm = tmp_path / "confirm"
    assert cli_main(["--prepare", "--output", str(explore)]) == 0
    assert cli_main(["--confirm", "--output", str(confirm), "--exploration", str(explore)]) == 0
    assert cli_main(["--check", "--output", str(explore)]) == 0
    assert cli_main(["--run", "--output", str(explore)]) == 2
    error = capsys.readouterr().err
    assert "continuation manifest is missing" in error


def test_full_catalog_can_be_admitted_without_opening_another_pool(tmp_path):
    record = family.prepare_specification(tmp_path / "all", nf=64, all_members=True)
    assert len(record["cases"]) == 15
    assert record["nf"] == 64
    assert record["forecast"]["workers"] == 6
    assert record["forecast"]["coordinator_pools"] == 1
    successor = episode.LAB / "results" / "development" / "nsc-discovery-family-v1"
    assert episode.assert_campaign_output(successor) == successor.resolve()


def test_materialization_keeps_failed_source_and_successful_actual_states(tmp_path, monkeypatch):
    directory = tmp_path / "materialized"
    record = family.prepare_specification(directory, nf=64)
    actual_prepare = family.prepare_member

    def one_failed_member(imbalance, width_scale, **kwargs):
        if imbalance == 0.0 and width_scale == 1.0:
            raise ValueError("deliberate dense-solve blocker")
        return actual_prepare(imbalance, width_scale, **kwargs)

    monkeypatch.setattr(family, "prepare_member", one_failed_member)
    manifest = family.materialize_directory(directory, duration=0.003, step_cap=0.003)
    assert 1 <= len(manifest["cases"]) <= 5
    assert len(manifest["cases"]) + len(manifest["preparation_failures"]) == 6
    failure = next(item for item in manifest["preparation_failures"] if item["case_id"] == family.case_id(0, 1, 64))
    assert failure["case_id"] == family.case_id(0, 1, 64)
    assert "deliberate dense-solve blocker" in failure["error"]
    assert failure["non_existence_claimed"] is False
    assert manifest["preparation_cpu_seconds"] >= manifest["initial_cpu_seconds"] + manifest["handoff_cpu_seconds"]
    assert manifest["handoff_cpu_seconds"] > 0
    assert manifest["cpu_budget_seconds"] == pytest.approx(21600 - manifest["preparation_cpu_seconds"])
    assert family.check(directory)["ok"]
    stored = json.loads((directory / family.FAMILY_RECORD).read_text())
    assert stored["initial_data"][0]["sample"]["actual_geometry"]["finite"]
    assert stored["initial_data"][0]["shift_residual_max"] >= abs(stored["initial_data"][0]["shift_residual_mean"])
    rows = family.assess_regimes(directory, stored, manifest)
    assert len(rows) == 6
    baseline = next(row for row in rows if row.get("imbalance") == 1 and row.get("width_scale") == 1)
    assert baseline["fate"] == "baseline"
    assert baseline["signals"]["actual_geometry_growth"]["R4_l2"] is not None
    assert np.asarray(baseline["final_sample"]["source_packet_window_fractions"]).shape == (3, 3)
    assert baseline["final_sample"]["metric_positive"]
    assert any(row["fate"] == "unresolved" for row in rows)


def test_confirmation_includes_matching_resolution_baseline(tmp_path):
    exploration = tmp_path / "exploration"
    family.prepare_specification(exploration)
    path = exploration / family.FAMILY_RECORD
    record = json.loads(path.read_text())
    record["regimes"] = [{"case_id": family.case_id(-1, 0.8, 128), "imbalance": -1,
                           "width_scale": 0.8, "fate": "distinct_physical_signs"}]
    path.write_text(json.dumps(record))
    result = family.prepare_confirmation(tmp_path / "confirmation", exploration)
    assert result["selected"]["cases"][0]["fate"] == "baseline_comparator"
    assert result["cases"] == [family.case_id(1, 1, 256), family.case_id(-1, 0.8, 256)]
    assert result["selected"]["cases"][1]["parent_case"] == family.case_id(-1, 0.8, 128)


def test_classifier_rejects_nonfinite_retention_and_invalid_budget():
    signals = family.physical_signals(length_change=-1, width_change=1,
                                     retained_fraction=float("nan"), packet_retention=(0.1, 0.2, 0.3))
    assert signals["unresolved"]
    with pytest.raises(ValueError):
        family.remaining_budget(float("nan"), 0)
