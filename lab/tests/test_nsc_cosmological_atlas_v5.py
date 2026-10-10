from copy import deepcopy
import json

import pytest

from test_nsc_cosmological_atlas_v4 import producer as predecessor


def producer():
    predecessor()
    import derive_nsc_cosmological_atlas_v5
    return derive_nsc_cosmological_atlas_v5


def test_preserves_other_routes_and_many_body_unknown():
    tool = producer()
    old, _, _, _ = tool.previous.build_record()
    record, causal, access = tool.build_record()
    assert record["scenarios"][:7] == old["scenarios"]
    assert len(record["scenarios"]) == 8
    assert record["fitted_observations"] == []
    assert not record["assessments"]["observer_access"]["global_many_body_future_solved"]
    assert next(row for row in access["rows"] if row["id"] == "siblings")["verdict"] == "unresolved"
    assert causal["reference"]["scope"] == "isolated spherical reference completion only"


def test_event_access_does_not_become_an_unconditional_message_edge():
    tool = producer()
    record, _, access = tool.build_record()
    scenario = record["scenarios"][-1]
    assert scenario["max_descent_depth"] == 0
    assert all(e["kind"] == "causal" and e["minimum_delay"] is None for e in scenario["links"])
    assert tool.first.atlas.earliest_arrival(scenario, "C", "P", 0., clock="tau") is None
    rows = {r["id"]: r for r in access["rows"]}
    assert rows["early-C-P"]["verdict"] == "allowed geometric path"
    assert rows["C0-P"]["verdict"] == "excluded in reference"
    assert not access["initial_fields_or_messages_prescribed"]
    assert not access["signal_amplitudes_calculated"]
    assert not access["centre_escape_threshold"]["upper_endpoint_admits_finite_parent_reception"]


def test_independent_wall_clocks_match_at_shared_events():
    tool = producer()
    _, causal, _ = tool.build_record()
    interior = causal["interior"]["named_events"]["W_early"]
    exterior = causal["reference"]["events"]["early_wall"]
    assert interior["tau_wall_over_ell"] == pytest.approx(exterior["tau_wall_over_ell"], abs=1e-12)
    legacy, _, _ = tool.build_record()
    legacy = deepcopy(legacy)
    legacy["schema"] = tool.first.atlas.SCHEMA
    assert tool.first.atlas.validate_atlas(legacy)


def test_creation_only_source_bound_diagram_and_access_table(tmp_path):
    tool = producer()
    output = tmp_path / "atlas-v5"
    assert tool.create(output)["outputs"] == 9
    assert tool.check(output)["outputs"] == 9
    m = json.loads((output / "manifest.json").read_bytes())
    assert "spacetime.png" in m["outputs"]
    assert "observer-access.md" in m["outputs"]
    assert any(p.endswith("atlas-working-draft-v5.md") for p in m["source_bindings"])
    with pytest.raises(FileExistsError):
        tool.create(output)
    p = output / "observer-access.json"
    p.write_bytes(p.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.check(output)
