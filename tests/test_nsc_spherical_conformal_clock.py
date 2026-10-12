"""Exact frozen-clock comparison binds the original observer and endpoint."""
import json
import numpy as np

import derive_nsc_spherical_conformal_episode as episode
import derive_nsc_spherical_conformal_continuation as continuation
import derive_nsc_spherical_conformal_clock as clock


def test_exact_clock_comparison_binds_source_and_removes_time_interpolation():
    record = json.loads(clock.OUT.read_text())
    assert record["source_unchanged"] is True
    assert record["source_bindings"]["episode_v2"] == episode.sha256(continuation.OUT)
    assert record["source_bindings"]["payload_v2"] == episode.sha256(continuation.NPZ)
    assert record["frozen_time_interpolation_used"] is False
    assert record["coupled_time_interpolation_used"] is False
    assert record["trajectory_evolved"] is False
    assert record["observable_error_certified"] is False
    weights = np.asarray(episode.galerkin.OCCUPATIONS)
    with np.load(continuation.NPZ, allow_pickle=False) as payload:
        for spec in episode.RUN_PLAN:
            row = record["results"][spec["name"]]
            assert abs(row["frozen_coordinate_time"] * row["original_clock_rate"] - row["proper_clock"]) < 1e-14
            assert row["frozen_coordinate_time"] <= row["coupled_coordinate_time"]
            observer = np.vstack((payload[f"nf{spec['nf']}_original_T0_phi0"], payload[f"nf{spec['nf']}_original_T0_phi1"]))[:, :2]
            columns = np.vstack((payload[spec["name"] + "_final_phi0"], payload[spec["name"] + "_final_phi1"]))
            projected = observer.conj().T @ columns
            occupation = np.sum(abs(projected) ** 2 * weights[None, :], axis=1)
            assert np.max(abs(occupation - row["coupled_occupation"])) < 1e-14
            assert np.max(abs(np.asarray(row["coupled_occupation"]) - np.asarray(row["frozen_occupation"]) - np.asarray(row["difference"]))) < 1e-14
    assert np.max(record["time_fraction_of_effect"]) < .01
    assert np.max(record["space_fraction_of_effect"]) < .01
