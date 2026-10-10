"""Keep local dynamics, existence families and physical attachments distinct."""
from copy import deepcopy
import json

import pytest

from test_nsc_cosmological_atlas_v2 import producer as predecessor


def producer():
    predecessor()  # establishes the repository's existing producer import path
    import derive_nsc_cosmological_atlas_v3
    return derive_nsc_cosmological_atlas_v3


def test_successor_preserves_all_previous_scenarios_and_fixed_evidence():
    tool = producer()
    old, _, _ = tool.previous.build_record()
    record, moving, _ = tool.build_record()
    assert record["scenarios"][:4] == old["scenarios"]
    assert record["fitted_observations"] == []
    assert moving["vacuum_wall_control"]["embedded_in_shared_parent"] is False
    assert moving["radiation_flux_screen"]["status"] == "MOTION_ALONE_INCOMPATIBLE_WITH_FIXED_SOURCE"
    assert len(record["scenarios"]) == 6


def test_signed_local_solution_is_neither_descent_nor_sibling_channel():
    tool = producer()
    record, moving, _ = tool.build_record()
    scenario = record["scenarios"][4]
    assert scenario["max_descent_depth"] == 0
    assert [e["kind"] for e in scenario["links"]] == ["attachment"]
    assert all(e["minimum_delay"] is None for e in scenario["links"])
    assert tool.first.atlas.earliest_arrival(scenario, "schwarzschild", "de_sitter", 0, clock="tau") is None
    wall = moving["vacuum_wall_control"]
    assert wall["orientation"]["epsilon_plus"] == -1
    assert wall["clock"]["cross_universe_age_established"] is False
    assert wall["child_formation_established"] is False
    assert wall["sibling_signalling_established"] is False
    assert wall["stability_established"] is False


def test_gluing_family_does_not_certify_numeric_separation_or_wall_data():
    tool = producer()
    record, _, gluing = tool.build_record()
    assert gluing["status"] == "IMPORTED_EXISTENCE_FAMILY"
    assert gluing["separation"]["numeric_threshold"] is None
    assert gluing["separation"]["BL_d12_certified"] is False
    assert gluing["moving_wall_attached"] is False
    assert gluing["explicit_metric_constructed_here"] is False
    assert gluing["unchanged_brill_lindquist"] is False
    assert "delta-layer" in gluing["non_time_symmetric_successor"]["missing"]
    parent = record["scenarios"][5]
    assert parent["max_descent_depth"] == 0
    assert all(e["status"] == "imported_existence_family" for e in parent["links"])
    legacy = deepcopy(record)
    legacy["schema"] = tool.first.atlas.SCHEMA
    assert tool.first.atlas.validate_atlas(legacy)


def test_source_bound_replay_and_creation_only_output(tmp_path):
    tool = producer()
    output = tmp_path / "atlas-v3"
    assert tool.create(output)["outputs"] == 8
    assert tool.check(output)["outputs"] == 8
    manifest = json.loads((output / "manifest.json").read_bytes())
    for name in ("nsc_atlas_moving_boundaries.py", "nsc-atlas-parent-gluing.md", "atlas-working-draft-v3.md"):
        assert any(path.endswith(name) for path in manifest["source_bindings"])
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(output)
    original = (output / "moving-boundary.json").read_bytes()
    (output / "moving-boundary.json").write_bytes(original + b"\n")
    with pytest.raises(ValueError, match="artifact binding mismatch: moving-boundary.json"):
        tool.check(output)
