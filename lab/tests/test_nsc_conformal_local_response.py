"""Conformal operator and bounded initial-window reduction feasibility."""
from copy import deepcopy
import json

import numpy as np
import pytest

from recursive_horizons import nsc_conformal_local_response as response
from recursive_horizons.nsc_coupled_local_response import sha256_file
from recursive_horizons.nsc_spherical_coupling import apply_dirac
from recursive_horizons.nsc_spherical_galerkin_coupling import build_grid


def test_conformal_matrix_matches_the_owned_fine_dirac_action_without_copy_factor():
    nf, nq, length = 16, 64, 8.
    coordinate = np.arange(nf - 1) * length / (nf - 1)
    base = .25 + .02 * np.sin(2 * np.pi * coordinate / length)
    inputs = {"nf": nf, "nq": nq, "length": length,
              "times": np.array([.025, .03, .035]), "Q": np.array([base, base * 1.01, base * 1.02])}
    schedule = response.ConformalSchedule(inputs)
    rng = np.random.default_rng(173)
    columns = rng.normal(size=(2 * nf, 6)) + 1j * rng.normal(size=(2 * nf, 6))
    grid = build_grid(nf, quadrature=nq, length=length, gauge="conformal")
    fine0, fine1 = grid.U_f @ columns[:nf], grid.U_f @ columns[nf:]
    q = schedule.radial(.0275)
    image0, image1 = apply_dirac(fine0, fine1, q, q, np.zeros(nq), grid.fine.kappa, grid.fine.momentum)
    expected = np.vstack((grid.U_f.conj().T @ image0, grid.U_f.conj().T @ image1))
    matrix = schedule(.0275)
    assert np.max(abs(matrix @ columns - expected)) < 2e-12
    assert np.max(abs(matrix - matrix.conj().T)) < 1e-12
    assert response.action_gate(schedule, columns, .0275) < 1e-12
    # The field generator is the representative block, without M=4 kappa.
    assert np.linalg.norm(grid.fine.multiplicity * matrix @ columns - expected) > 1.
    np.testing.assert_allclose(matrix, .5 * (schedule(.025) + schedule(.03)), rtol=0, atol=1e-12)
    assert len(schedule.cache) <= 2
    with pytest.raises(ValueError, match="saved conformal segment"):
        schedule(.05)


def test_actual_correlated_segment_keeps_original_observer_and_weights():
    inputs = response.load_final_episode()
    blocks = response.initial_blocks(inputs)
    assert blocks["C_AE_frobenius"] == pytest.approx(.176498626520663, abs=1e-12)
    assert blocks["covariance_admissible"] is True
    assert not np.array_equal(inputs["observer"], inputs["columns"][:, :2])
    assert blocks["observer_gram_gap"] < 1e-13
    assert inputs["times"][0] == pytest.approx(.025)
    assert inputs["times"][-1] == pytest.approx(.05)
    assert inputs["regional_exchange_resolved"] is False
    with pytest.raises(RuntimeError, match="feasibility only"):
        response.execute(inputs=inputs)


@pytest.mark.parametrize("field,changed,error", [
    ("gauge", "prescribed", "schema_or_gauge"),
    ("verdict", "PARTIAL_DIAGNOSTIC_CHECKPOINT", "unfinished_episode"),
    ("source_changed", True, "changed_preparation"),
    ("bound_sources_unchanged_during_run", False, "source_binding"),
])
def test_readiness_rejects_changed_or_partial_episode(field, changed, error):
    record = json.loads(response.EPISODE_JSON.read_text())
    assert response.episode_readiness_errors(record) == []
    altered = deepcopy(record)
    altered[field] = changed
    assert error in response.episode_readiness_errors(altered)


def test_readiness_holds_reduction_on_failed_chart_or_constraint_assessment():
    record = json.loads(response.EPISODE_JSON.read_text())
    altered = deepcopy(record)
    altered["results"][response.PRIMARY]["positive_chart"] = False
    altered["results"][response.PRIMARY]["constraint_end_margin"] = -1e-8
    errors = response.episode_readiness_errors(altered)
    assert "chart_or_completion:" + response.PRIMARY in errors
    assert "sampled_constraint_budget:" + response.PRIMARY in errors


def test_saved_probe_is_feasibility_only_and_bound_to_unchanged_sources():
    record = json.loads(response.PROBE_JSON.read_text())
    assert record["status"] == "MEASURED_REDUCTION_FEASIBILITY"
    assert record["programme_complete"] is False
    assert record["production_reduction_executed"] is False
    assert record["full_segment_comparison_executed"] is False
    assert record["regional_throughflow_resolved"] is False
    assert record["time_window"] == [.025, .035]
    assert record["cpu_seconds"] < record["cpu_budget_seconds"] <= 60
    assert record["payload_bytes"] == response.PROBE_JSON.stat().st_size + response.PROBE_NPZ.stat().st_size
    assert record["payload_bytes"] < 64 * 1024 * 1024
    assert record["source_code_sha256"] == sha256_file(response.__file__)
    assert record["payload_sha256"] == sha256_file(response.PROBE_NPZ)
    assert record["allocation"]["time_indexed_exterior_propagator_bytes"] == 0
    assert record["allocation"]["history_shape"] == [3, 1022, 6]
    for path in (response.EPISODE_JSON, response.EPISODE_NPZ, response.V5_NPZ):
        assert record["source_bindings"][path.name] == sha256_file(path)
    with np.load(response.PROBE_NPZ, allow_pickle=False) as saved:
        actual_error = float(np.max(abs(saved["occupation_retained"] - saved["occupation_full"])))
        assert actual_error == record["occupation_error"]["max_abs"]
        assert actual_error < .01 * record["occupation_change"]
        for name in ("memory", "outside_drive"):
            movement = float(np.max(abs(saved["occupation_" + name + "_off"] - saved["occupation_full"])))
            assert movement == record["controls"][name]["movement"]
            assert movement > 20 * actual_error
        cross = saved["covariance_retained"] - saved["covariance_without_cross"]
        movement = float(np.max(abs(np.diagonal(cross, axis1=1, axis2=2).real)))
        assert movement == pytest.approx(record["controls"]["initial_cross"]["movement"], abs=1e-15)
        assert movement > 20 * actual_error
