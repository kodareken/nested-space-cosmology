"""Prevent imported spatial data and a scoped rejection becoming stronger claims."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
import sys

import pytest


def producer():
    scripts = Path(__file__).resolve().parents[1] / "scripts"
    if str(scripts) not in sys.path:
        sys.path.insert(0, str(scripts))
    spec = importlib.util.spec_from_file_location("atlas_v2_producer", scripts / "derive_nsc_cosmological_atlas_v2.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_shared_parent_is_one_spatial_metric_without_formed_children():
    tool = producer()
    record, parent, _ = tool.build_record()
    regions = {r["id"]: r for r in record["scenarios"][0]["regions"]}
    assert record["schema"] == tool.SCHEMA
    assert parent["child_formation_established"] is False
    assert parent["direct_spherical_child_attachment_established"] is False
    assert parent["adm_masses"]["puncture_ends"] == [1.0625, 1.5625]
    assert regions["parent"]["geometry"]["globally_spherical"] is False
    assert regions["hole_a"]["geometry"]["exact_schwarzschild_patch"] is False
    assert regions["hole_a"]["clock"]["kind"] == regions["hole_b"]["clock"]["kind"]
    assert record["fitted_observations"] == []


def test_incompatible_wall_does_not_erase_unrelated_branches():
    tool = producer()
    record, _, result = tool.build_record()
    assert result["status"] == "INCOMPATIBLE_FIXED_SOURCE"
    assessment = record["assessments"]["fixed_source_junction"]
    assert assessment["all_child_attachments_excluded"] is False
    assert assessment["embedded_in_shared_parent"] is False
    assert record["scenarios"][1:] == tool.first.atlas.build_atlas()["scenarios"][1:]
    necks = [e for e in record["scenarios"][0]["links"] if e["id"].startswith("neck-")]
    assert len(necks) == 2 and all(e["status"] == "unresolved" for e in necks)


def test_spatial_connections_cannot_supply_causal_delays_or_reset_depth():
    tool = producer()
    record, _, _ = tool.build_record()
    scenario = record["scenarios"][0]
    assert tool.first.atlas.earliest_arrival(scenario, "child_a", "child_b", 0, clock="tau") is None
    assert all(e["minimum_delay"] is None for e in scenario["links"])
    legacy_shape = deepcopy(record)
    legacy_shape["schema"] = tool.first.atlas.SCHEMA
    legacy_shape["scenarios"][0]["max_descent_depth"] = 1
    with pytest.raises(ValueError, match="cutoff"):
        tool.first.atlas.validate_atlas(legacy_shape)


def test_readonly_replay_overwrite_refusal_and_artifact_tampering(tmp_path):
    tool = producer()
    output = tmp_path / "atlas-v2"
    result = tool.create(output)
    assert result["new_physical_mechanism"] is False
    checked = tool.check(output)
    assert checked["outputs"] == 12
    manifest = json.loads((output / "manifest.json").read_text())
    for name in ("nsc_atlas_parent_initial_data.py", "nsc_atlas_junction_sources.py"):
        assert any(p.endswith(name) for p in manifest["source_bindings"])
    before = (output / "manifest.json").read_bytes()
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(output)
    assert (output / "manifest.json").read_bytes() == before
    with (output / "junction.json").open("ab") as handle:
        handle.write(b"\n")
    with pytest.raises(ValueError, match="binding mismatch: junction.json"):
        tool.check(output)
