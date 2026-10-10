"""Initial-data attachment cannot silently become descent or communication."""
from copy import deepcopy
import json

import pytest

from test_nsc_cosmological_atlas_v3 import producer as predecessor


def producer():
    predecessor()
    import derive_nsc_cosmological_atlas_v4
    return derive_nsc_cosmological_atlas_v4


def test_preserves_old_routes_and_marks_historical_interface_status():
    tool = producer()
    old, _, _ = tool.previous.build_record()
    record, _, _, _ = tool.build_record()
    assert record["scenarios"][:6] == old["scenarios"]
    assert len(record["scenarios"]) == 7
    assert record["fitted_observations"] == []
    assert "moving_vacuum_wall" in record["historical_scenario_status"]
    legacy = deepcopy(record)
    legacy["schema"] = tool.first.atlas.SCHEMA
    assert tool.first.atlas.validate_atlas(legacy)


def test_initial_attachment_does_not_promote_causal_path_or_child():
    tool = producer()
    record, connection, assembly, observer = tool.build_record()
    scenario = record["scenarios"][-1]
    assert scenario["max_descent_depth"] == 0
    assert all(e["kind"] == "attachment" for e in scenario["links"])
    assert tool.first.atlas.earliest_arrival(scenario, "parent", "interior_a", 0, clock="tau") is None
    assert not assembly["conclusions"]["formation_history"]
    assert not assembly["conclusions"]["global_shell_evolution"]
    assert not assembly["conclusions"]["perturbation_stability"]
    assert connection["theory"]["fixed_source_within_test"]
    assert not observer["window"]["upper_endpoint_included"]


def test_end_only_hypotheses_and_second_end_survive():
    tool = producer()
    _, connection, assembly, _ = tool.build_record()
    assert assembly["bodies"][0]["mass_over_ell"] != assembly["bodies"][1]["mass_over_ell"]
    assert assembly["hypotheses"]["selected_inputs_smooth_vacuum_AE"]
    assert not assembly["hypotheses"]["whole_body_smooth_corollary_used"]
    assert assembly["topology"]["total_AE_ends"] == 2
    assert not assembly["topology"]["second_cosmological_child_established"]
    assert assembly["separation"]["numeric_threshold"] is None
    assert not assembly["separation"]["BL_d12_certified"]
    broken = deepcopy(connection)
    broken["parent_end"]["exact_vacuum_overlap_available"] = False
    with pytest.raises(ValueError):
        tool.parent_assembly_record(broken)


def test_bound_replay_and_creation_only_output(tmp_path):
    tool = producer()
    output = tmp_path / "atlas-v4"
    assert tool.create(output)["outputs"] == 9
    assert tool.check(output)["outputs"] == 9
    manifest = json.loads((output / "manifest.json").read_bytes())
    for name in ("nsc_atlas_connection_slice.py", "nsc-atlas-connection-gluing.md", "atlas-working-draft-v4.md"):
        assert any(p.endswith(name) for p in manifest["source_bindings"])
    with pytest.raises(FileExistsError):
        tool.create(output)
    artifact = output / "causal-window.json"
    artifact.write_bytes(artifact.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.check(output)
