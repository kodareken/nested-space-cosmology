"""Independent checks for the sealed discovery-episode consumer.

The campaign files are read only. Writes stay in a temporary directory.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

import assess_nsc_discovery_episode as assessment

EPISODE = assessment.DEFAULT_EPISODE
FINE = "nf256_coupled_dt0.0005"


def test_missing_defect_stays_unrecorded_and_crossing_is_not_a_gate():
    assert assessment.recorded_defect({"child_normal_energy": 1.0}) is None
    assert assessment.recorded_defect({"normal_energy_projection_defect": 0.0}) == 0.0
    times = [0.0, 1.0, 2.0, 3.0, 4.0]
    gaps = [0.0, 0.2, 1.5, 0.4, 2.0]
    steps = [None, 1.0, 1.0, 1.0, 1.0]
    assert assessment.first_gap_exceeding_step(times, gaps, steps) == 2.0
    assert assessment.relative_gap(34240.3916208148, 29030.988792300224) == pytest.approx(0.15214203406913634)


def test_nodal_trace_uses_period_and_weights():
    weights = np.array([0.75, 0.75, 0.5, 0.5, 0.25, 0.25])
    phi0 = np.zeros((8, 6), dtype=complex)
    phi1 = np.zeros((8, 6), dtype=complex)
    phi0[1, 0] = 1.0
    trace = assessment.occupation_trace(phi0, phi1, weights)
    assert trace == pytest.approx(0.75)
    # x_i = i on this 8-node period, so node 1 lies in the child window [1, 3).
    assert assessment.nodal_window(phi0, phi1, weights, 1.0, 3.0) == pytest.approx(0.75)
    assert assessment.nodal_window(phi0, phi1, weights, 0.0, 1.0) == pytest.approx(0.0)


def test_reversed_rows_still_match_the_station_time():
    rows = [
        {"time": 3.0, "sample_kind": "station", "value": "late"},
        {"time": 0.3, "sample_kind": "start", "value": "handoff"},
        {"time": 1.0, "sample_kind": "station", "value": "middle"},
    ]
    assert assessment.row_at(list(reversed(rows)), 1.0)["value"] == "middle"
    assert assessment.row_at(rows, 3.0)["sample_kind"] == "station"
    with pytest.raises(RuntimeError, match="expected one observation"):
        assessment.row_at(rows + [rows[0]], 3.0)
    assert [(left["time"], right["time"]) for left, right in assessment.match_rows(rows, list(reversed(rows)))] == [(0.3, 0.3), (1., 1.), (3., 3.)]
    with pytest.raises(RuntimeError, match="counts differ"):
        assessment.match_rows(rows, rows[:-1])
    with pytest.raises(RuntimeError, match="expected one observation"):
        assessment.match_rows(rows, [rows[0], rows[1], {"time": 1.1}])


def test_initial_subtracted_control_effect_uses_distinct_initials():
    effects = assessment._control_effect({"R": 10.}, {"R": 13.}, {"R": 1.}, {"R": 3.}, ["R"])
    assert effects["R"]["initial_subtracted_effect"] == 1.
    assert effects["R"]["coupled_change"] == 3.
    assert effects["R"]["frozen_change"] == 2.


def test_periodic_integral_resolves_a_mode_and_constant():
    coordinates = np.arange(32) * assessment.PERIOD / 32
    values = 2. + np.cos(2. * np.pi * coordinates / assessment.PERIOD)
    assert assessment.periodic_integral(values, 0., 2.) == pytest.approx(4. + 4. / np.pi)
    assert assessment.periodic_integral(values, 0., 8.) == pytest.approx(16.)


def test_sealed_batch_authentication_and_physical_anchors():
    npz = EPISODE / f"{FINE}-000002.npz"
    chunk = json.loads((EPISODE / f"{FINE}-000002.json").read_text())
    digest_before = hashlib.sha256(npz.read_bytes()).hexdigest()
    report = assessment.assess(EPISODE, assessment.DEFAULT_ATLAS)
    digest_after = hashlib.sha256(npz.read_bytes()).hexdigest()

    assert digest_before == digest_after == chunk["arrays_sha256"]
    assert report["ok"] is True
    assert report["problems"] == []
    assert report["hashes_unchanged"] is True
    assert report["sealed_bytes_rewritten"] is False
    assert report["producing_commit"] == assessment.PRODUCING_COMMIT
    assert report["n_cases"] == 6
    assert report["n_chunks"] == 18
    assert report["cpu_seconds"] == pytest.approx(277.185634)
    assert report["consumer_sha256"] == hashlib.sha256(Path(assessment.__file__).read_bytes()).hexdigest()
    assert report["read_hashes"]["lab/results/development/nsc-discovery-episode-v1/nf256_coupled_dt0.0005-000002.npz"] == digest_before

    rows = [json.loads(line) for line in (EPISODE / f"{FINE}-observations.jsonl").read_text().splitlines()]
    matched = [row for row in reversed(rows) if abs(row["time"] - chunk["coordinate_time"]) <= 1e-8]
    assert len(matched) == 1
    assert matched[0]["sample_kind"] == "station"
    for key, value in chunk["channel_sample"].items():
        assert matched[0][key] == value

    late = report["cases"][FINE]["stations"]["3.0"]
    assert late["normal_energy_projection_defect"] is None
    assert late["boundary_incoming"] is None
    assert late["boundary_outgoing"] is None
    assert late["energy_accounts"]["energy_Q"] is None
    assert late["energy_accounts"]["energy_field"] == pytest.approx(11.989766083058845)
    assert late["Q_ddot_max"] is None
    assert late["chart_Q_min"] == pytest.approx(2.287805170599047e-05)
    assert late["chart_r_min"] == pytest.approx(20.170607955843877)
    assert late["occupation_trace"] == pytest.approx(3.0, abs=1e-8)
    assert late["child_probability"] == pytest.approx(0.0030091976292384495)
    for name in ("child_probability", "parent_probability", "child_proper_length", "parent_proper_length"):
        assert late["independent_metric_content"][name] == pytest.approx(late[name], abs=1e-10)
    assert late["independent_metric_content"]["clock_rates"] == pytest.approx(late["clock_rates"], abs=1e-10)
    assert report["cases"]["nf256_frozen_geometry_dt0.0005"]["stations"]["3.0"]["Q_ddot_max"] == 0.0

    claims = report["claims"]
    assert claims["normal_energy_projection_defect"] is None
    assert claims["defect_filled_with_zero"] is False
    assert claims["shell_closure_claimed"] is False
    assert claims["incoming_outgoing_split"] == "needed"
    assert claims["one_percent_gate_installed"] is False
    assert claims["analytic_qddot_closes_T3_spatial_gap"] is False
    assert claims["coordinate_station_is_common_proper_time"] is False
    assert claims["stability_theorem"] is False
    assert claims["continuous_evolution_error_bound"] is None
    assert claims["tidal_response_computed_by_this_consumer"] is False

    outcome = report["physical_outcome"]
    assert outcome["T1_R_h_max_space_relative_gap"] == pytest.approx(0.00023250991918969066)
    assert outcome["T3_R_h_max_space_relative_gap"] == pytest.approx(0.15214203406913634)
    assert outcome["T3_weyl_max_space_relative_gap"] == pytest.approx(0.2814915944647413)
    assert outcome["T3_Q_min_space_relative_gap"] == pytest.approx(0.0005690455776077915)
    assert outcome["T3_r_min_space_relative_gap"] == pytest.approx(9.764710626273575e-05)
    assert outcome["T3_child_length_space_relative_gap"] == pytest.approx(0.00046149137189203167)
    assert outcome["T3_R_h_max_time_relative_gap"] == pytest.approx(9.083446646513002e-07)
    assert outcome["late_R_h_gap_exceeds_fine_step_at"] == pytest.approx(2.9)
    assert outcome["late_weyl_gap_exceeds_fine_step_at"] == pytest.approx(2.65)
    assert outcome["areal_r_min_coupled_T3"] == pytest.approx(20.170607955843877)
    assert outcome["proper_length_contraction_factor"] == pytest.approx(0.0004296185078083378)
    control = report["coupled_vs_frozen"]["stations"]["256"]["3.0"]
    assert control["child_probability_over_trace"]["initial_subtracted_effect"] == pytest.approx(-0.0031415689552289905)
    assert control["child_proper_length"]["initial_subtracted_effect"] == pytest.approx(-2.374197705694078)
    assert control["child_probability_per_proper_length"]["initial_subtracted_effect"] == pytest.approx(2.943693543355158)
    assert control["R_h_max"]["initial_subtracted_effect"] == pytest.approx(29015.465414171358)
    assert control["R_h_max"]["coupled_initial"] != control["R_h_max"]["frozen_initial"]

    coupled_clock = report["clocks"][FINE]["1.0->3.0"]
    frozen_clock = report["clocks"]["nf256_frozen_geometry_dt0.0005"]["1.0->3.0"]
    assert coupled_clock["coordinate_dt"] == pytest.approx(2.0)
    assert coupled_clock["proper_dt"][0] == pytest.approx(0.24295814836255603)
    assert frozen_clock["coordinate_dt"] == pytest.approx(2.0)
    assert frozen_clock["proper_dt"][0] == pytest.approx(2.417449033714278)
    assert coupled_clock["observer_coordinates"] == [1., 2., 3.]
    assert "not matched proper times" in coupled_clock["scope"]
    assert report["balance"]["max_work_integral_abs_discrepancy"] == 0.0
    assert report["balance"]["shell_closure_claimed"] is False
    assert report["atlas_handoff"]["normal_energy_projection_defect"] is None
    assert report["atlas_handoff"]["child_balance_residual_without_recorded_defect"] == pytest.approx(0.011538817113910853)


def test_changed_observation_bytes_fail_batch_binding_without_input_writes(monkeypatch):
    original = assessment.read_hashed
    path = EPISODE / f"{FINE}-observations.jsonl"
    before = path.read_bytes()

    def changed(candidate):
        data, digest = original(candidate)
        if candidate == path:
            # Valid JSON with identical rows still differs from the frozen bytes.
            data += b"\n"
            digest = assessment.sha256_bytes(data)
        return data, digest

    monkeypatch.setattr(assessment, "read_hashed", changed)
    report = assessment.assess()
    assert report["ok"] is False
    assert any("does not match the authenticated batch closure" in item for item in report["problems"])
    assert path.read_bytes() == before


def test_write_is_creation_only_and_stays_outside_the_episode(tmp_path):
    destination = tmp_path / "record" / "assessment.json"
    assert assessment.main(["--write", str(destination), "--episode", str(EPISODE)]) == 0
    written = json.loads(destination.read_text())
    assert written["ok"] is True
    assert written["producing_commit"] == assessment.PRODUCING_COMMIT
    with pytest.raises(FileExistsError):
        assessment.main(["--write", str(destination), "--episode", str(EPISODE)])
    forbidden = EPISODE / "assessment-should-not-exist.json"
    with pytest.raises(PermissionError):
        assessment.main(["--write", str(forbidden), "--episode", str(EPISODE)])
    assert not forbidden.exists()
