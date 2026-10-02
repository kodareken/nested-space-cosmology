"""Fixed-source sector constraints, canonical rates, actual curvature and forks."""
import json
from pathlib import Path
import tempfile

import numpy as np
import pytest
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_initial_sector as sector
from recursive_horizons import nsc_discovery_dynamic_preparation as prep
from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_nested_parent_child as nested
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin


@pytest.fixture(scope="module")
def prepared():
    # Production pins need the future root commit; this unit fixture does
    # not claim a producing revision or create a production record.
    patch = pytest.MonkeyPatch()
    patch.setattr(prep, "_git_hashes", lambda commit, hashes: commit)
    with threadpool_limits(limits=1):
        result = sector.build_record(execute=True, producer_commit="unit-frozen-producer")
    yield result
    patch.undo()


def test_three_new_sectors_preserve_frozen_source_and_observer(prepared):
    report, arrays = prepared
    assert report["status"] == "FINITE_INITIAL_SECTOR_DATA"
    assert set(report["cases"]) == {"h_minus", "h_plus", "chi_minus2"}
    pair, old, _ = prep.load_case(sector.PREPARED, "uniform")
    assert np.array_equal(arrays["source_phi0"], old.phi0)
    assert np.array_equal(arrays["source_phi1"], old.phi1)
    assert np.array_equal(arrays["observer_columns"], pair.reference_columns)
    assert np.array_equal(arrays["W"], pair.geometry_map)
    for name in sector.SECTORS:
        np.testing.assert_array_equal(arrays[name + "_phi0"], old.phi0)
        np.testing.assert_array_equal(arrays[name + "_phi1"], old.phi1)
        measured = report["cases"][name]
        assert measured["gram_max"] < 1e-12
        assert measured["CAR_min"] >= -1e-12
        assert measured["CAR_max"] <= 1 + 1e-12
        assert measured["source_trace"] == pytest.approx(3., abs=1e-12)
        assert measured["multiplicity"] == 4
        assert not measured["current_mean_deleted"]
        assert not measured["field_frame_changed"]
        assert not measured["observer_changed"]
        assert not measured["coefficients_changed"]
    assert report["old_comparison"]["rerun"] is False
    assert all(row["producing_commit"] == "67ab5700229216612e7383dc23d0ae1ce6e479c0"
               for row in report["old_comparison"]["cases"])
    assert all(row["initial_state_sha256"] == report["baseline_prepared_state_sha256"]
               for row in report["old_comparison"]["cases"])


def test_null_kinetic_family_uses_canonical_pi_not_nodal_momenta(prepared):
    report, arrays = prepared
    pair, old, _ = prep.load_case(sector.PREPARED, "uniform")
    nodal = nested.reconstruct_state(pair, old)
    a = -8 * np.pi * pair.grid.fine.A
    f = -8 * np.pi * pair.grid.fine.C_W / 3
    for name, h in (("h_minus", -.5), ("h_plus", .5)):
        assert np.array_equal(arrays[name + "_nodal_r"], nodal.r)
        np.testing.assert_allclose(arrays[name + "_nodal_p_chi"], 2 * f * h)
        np.testing.assert_allclose(arrays[name + "_nodal_p_r"], 2 * a * nodal.r * h)
        expected = pair.grid.dx_g * pair.geometry_map.T @ arrays[name + "_nodal_p_chi"]
        np.testing.assert_allclose(arrays[name + "_p_chi"], expected, atol=1e-15)
        assert not np.allclose(arrays[name + "_p_chi"], arrays[name + "_nodal_p_chi"])
        measured = report["cases"][name]
        assert measured["full_lapse_max"] < 1e-6
        assert measured["full_shift_max"] < 1e-6
        assert abs(measured["source_current_mean"]) < 1e-11
        assert measured["Pi_max"] < 1e-10
        assert abs(measured["gravity_kinetic_budget"]["total"]) < 1e-15
        assert measured["Qdot"]["min"] == pytest.approx(h * np.mean(nodal.Q), abs=2e-11)
        assert measured["rdot"]["max_abs"] < 1e-11
        assert measured["chidot"]["max_abs"] < 1e-9
        expected_K = h / (np.mean(nodal.r) * np.mean(nodal.Q))
        assert measured["K_radial"]["min"] == pytest.approx(expected_K, abs=2e-10)
        assert measured["K_angular"]["max_abs"] < 1e-11


def test_chi_minus2_radius_is_selected_from_its_own_actual_constraint(prepared):
    report, arrays = prepared
    old_r = arrays["h_plus_nodal_r"]
    new_r = arrays["chi_minus2_nodal_r"]
    assert np.linalg.norm(new_r - old_r) > .001
    np.testing.assert_allclose(arrays["chi_minus2_nodal_chi"], -2.)
    for key in nested.MOMENTUM_NAMES:
        assert np.all(arrays["chi_minus2_nodal_" + key] == 0)
    measured = report["cases"]["chi_minus2"]
    assert measured["radius_preparation"]["performed"]
    assert not measured["radius_preparation"]["old_radius_copied"]
    assert measured["full_lapse_max"] < 1e-6
    assert measured["full_shift_max"] < 1e-6
    pair, _old, _ = prep.load_case(sector.PREPARED, "uniform")
    alpha = -4 * np.pi * pair.grid.fine.C_W / 3
    expected_shift = -alpha / (2 * np.pi * pair.grid.fine.A)
    assert np.mean(new_r ** 2 - old_r ** 2) == pytest.approx(expected_shift, abs=2e-9)
    assert measured["momentum_rates"]["p_chi"]["max_abs"] < 1e-10
    assert measured["momentum_rates"]["p_r"]["max_abs"] > 1
    assert measured["momentum_rates"]["p_Q"]["max_abs"] > 1
    assert measured["actual_Qtt"]["max_abs"] < 1e-8
    assert measured["actual_rtt"]["max"] < 0


def test_actual_curvature_and_source_energy_distinguish_the_sector_changes(prepared):
    report, _arrays = prepared
    baseline = report["baseline_initial_measurement"]
    for name in ("h_minus", "h_plus"):
        measured = report["cases"][name]
        assert measured["actual_metric_curvature"]["R_h"]["min"] == pytest.approx(2., abs=2e-7)
        assert measured["actual_metric_curvature"]["owned_W"]["max_abs"] < 1e-12
        assert measured["K_radial"]["max_abs"] > .1
    auxiliary = report["cases"]["chi_minus2"]
    assert auxiliary["actual_metric_curvature"]["R_h"]["max_abs"] < 2e-7
    expected_W = 4 / (3 * auxiliary["r"]["min"] ** 4)
    assert auxiliary["actual_metric_curvature"]["owned_W"]["max_abs"] == pytest.approx(expected_W, rel=2e-7)
    for measured in report["cases"].values():
        assert not measured["chi_substituted_for_curvature"]
        assert measured["source_coordinate_energy"] == pytest.approx(baseline["source_coordinate_energy"], abs=2e-10)
        assert not measured["future_rates_fixed"]


def test_initial_h_is_not_an_ongoing_fixed_rate(prepared):
    _report, arrays = prepared
    pair, _old, _ = prep.load_case(sector.PREPARED, "uniform")
    state = nested.NestedState(*(arrays["h_plus_" + name] for name in nested.STATE_NAMES))
    with threadpool_limits(limits=1):
        advanced = nested.rk4_step(pair, state, 1e-4)
    assert np.linalg.norm(advanced.p_chi - state.p_chi) > 1e-9
    assert not np.array_equal(advanced.phi0, state.phi0)
    assert not np.array_equal(advanced.chi, state.chi)


def test_creation_only_record_and_exact_six_cap_forks(prepared):
    report, arrays = prepared
    with tempfile.TemporaryDirectory(dir="/tmp", prefix="nsc-initial-sector-") as directory:
        record = Path(directory) / "preparation"
        output = Path(directory) / "episode"
        sector.write_record(report, arrays, record)
        assert sector.check_record(record)["ok"]
        with pytest.raises(FileExistsError):
            sector.write_record(report, arrays, record)
        manifest = sector.prepare_episode(record, output, execute=True, producer_commit="unit-frozen-producer")
        assert len(manifest["cases"]) == 6
        assert manifest["stations"] == [.3, 1., 3.]
        assert manifest["cpu_budget_seconds"] == 300.
        assert manifest["forecast_factor"] == 1.5
        assert manifest["evolved"] is False
        for name in sector.SECTORS:
            left, la = episode.load_checkpoint(output, f"nf128_{name}_dt0.001")
            right, ra = episode.load_checkpoint(output, f"nf128_{name}_dt0.0005")
            assert left["source_pins"]["initial_state_sha256"] == right["source_pins"]["initial_state_sha256"]
            for key in la:
                np.testing.assert_array_equal(la[key], ra[key])
            assert left["future_rates_fixed"] is False
            assert left["initial_state_called"] is False
        assert sector.assess_episode(output)["old_baseline_rerun"] is False
        with pytest.raises(FileExistsError):
            sector.prepare_episode(record, output, execute=True, producer_commit="unit-frozen-producer")


def test_new_production_requires_freeze_and_bounds():
    with pytest.raises(ValueError, match="frozen"):
        sector.build_record(execute=True)
    with pytest.raises(ValueError, match="CPU"):
        sector.build_record(cpu_limit=31.)
    with pytest.raises(ValueError, match="300 CPU"):
        sector.run_episode("unused", cpu_budget=301.)
    pair, base, _ = prep.load_case(sector.PREPARED, "uniform")
    with pytest.raises(ValueError, match="declared sector"):
        sector.prepare_uniform_sector(pair, base, "h_large")
    bad = base.copy(); bad.Q += .01
    with pytest.raises(ValueError, match="uniform"):
        sector.prepare_uniform_sector(pair, bad, "h_plus")


def test_existing_worker_one_step_keeps_initial_h_free_and_retains_pins(prepared):
    report, arrays = prepared
    with tempfile.TemporaryDirectory(dir="/tmp", prefix="nsc-sector-step-") as directory:
        initial = Path(directory) / "initial"
        output = Path(directory) / "episode"
        sector.write_record(report, arrays, initial)
        sector.prepare_episode(initial, output, execute=True, producer_commit="unit-frozen-producer")
        identifier = "nf128_h_plus_dt0.0005"
        first, before = episode.load_checkpoint(output, identifier)
        with threadpool_limits(limits=1), prep.episode_adapter():
            result = prep._dynamic_case_worker({"directory": str(output), "case_id": identifier,
                "cpu_budget_seconds": 300., "forecast_factor": 1.5, "memory_limit_bytes": 2**30,
                "max_steps": 1, "backend": "fft"})
        assert result["status"] == "step_limit"
        assert result["steps"] == 1
        last, after = episode.load_checkpoint(output, identifier)
        assert last["coordinate_time"] > 0
        assert np.linalg.norm(after["pi_chi"] - before["pi_chi"]) > 1e-9
        assert not np.array_equal(after["phi0"], before["phi0"])
        assert np.array_equal(after["source_phi0"], before["source_phi0"])
        assert last["source_pins"]["initial_state_sha256"] == first["source_pins"]["initial_state_sha256"]
        assert last["future_rates_fixed"] is False
