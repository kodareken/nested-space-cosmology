"""Read-only seals for the parent-episode assessment. No evolution and no new records."""
from pathlib import Path

import assess_nsc_discovery_parent_episode as assessment

DEFAULT_DIGEST = "79dbc6a7bca6010c160a92784dd5d46f57e56ba29b79e6e9d77eb5b79f992bdf"
LAB = Path(assessment.episode.LAB)


def test_default_successor_preserves_measurement_and_authenticates_current_source():
    result = assessment.report([LAB/path for path in assessment.DEFAULTS])
    assert result["schema"] == assessment.SCHEMA
    assert result["predecessor"]["content_sha256"] == DEFAULT_DIGEST
    assert result["measurement_unchanged_from_predecessor"] is True
    owner = str(Path(assessment.__file__).resolve())
    assert result["measurement_source_hashes"][owner] == assessment.file_hash(owner)
    assert result["comparison"]["autonomous_renewal"] is False
    assert result["comparison"]["late_curvature_is_a_continuum_value"] is False


def test_responsive_controls_recompute_arrays_and_agree():
    directories = [LAB/path for path in assessment.RESPONSIVE]
    first = assessment.responsive_report(directories)
    second = assessment.responsive_report(directories)
    assert first["content_sha256"] == second["content_sha256"]
    assert first["header_copied_into_measurement"] is False
    assert first["instability_claimed"] is False
    assert first["new_threshold"] is None
    assert first["late_curvature_is_an_instability"] is False
    strong = first["attribution"]["strong_versus_magnetic_plus_balanced_handoff"]
    assert strong["k"]["authoritative"] == "parent_k"
    assert strong["k"]["weak"] == 0.002621213673327129
    assert strong["k"]["strong"] == 0.45388971484879426
    assert strong["k"]["legacy_common_k"]["authoritative"] is False
    assert strong["source_force_alone"] is False
    assert strong["initial_Q_gap"] == 0.0
    assert strong["initial_weight_gap"] > 0.2
    assert strong["initial_pi_r_gap"] > 0.2
    frozen = first["attribution"]["frozen_versus_magnetic_plus_parent_heavy_handoff"]
    assert frozen["same_canonical_handoff"] is True
    assert frozen["initial_state_gap"] == 0.0
    assert frozen["geometry_change_at_T3"] == {"Q": 0.0, "r": 0.0}
    campaign = next(item for item in first["campaigns"] if item["directory"].endswith("strong-t3-v1"))
    stale = [audit for audit in campaign["header_audits"] if audit["header_stale"]]
    assert {audit["coordinate_time"] for audit in stale} == {1.25, 1.5, 2.25}
    assert all(audit["header_time"] == 1.0 and audit["header_used_as_measurement"] is False for audit in stale)
    plus = next(row for row in campaign["stations"] if row["sign_name"] == "plus" and abs(row["time"]-1.25) < 1e-8)
    audit = plus["checkpoint_header_audit"]
    assert audit["used_as_measurement"] is False
    assert abs(audit["length_header_minus_recomputed"]) > 0.05
    assert plus["proper_lengths"]["child"] == audit["recomputed_child_proper_length"]
    assert {row["time"] for row in first["selected_strong_observations"]} == {1.0, 1.25, 1.5, 2.25, 3.0}
