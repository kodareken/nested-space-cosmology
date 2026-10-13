"""Successor integration preserves provenance and the scope of the finite response."""
from copy import deepcopy
import json
import xml.etree.ElementTree as ET

import pytest

def producer():
    import derive_nsc_cosmological_atlas_v8
    return derive_nsc_cosmological_atlas_v8


@pytest.fixture(scope="module")
def records():
    tool = producer()
    return tool, tool.build_record()


def test_evaluated_response_preserves_all_ten_predecessor_scenarios(records):
    tool, (atlas, programme, _, response) = records
    previous, _, _ = tool.previous.build_record()
    assert len(atlas["scenarios"]) == 11
    assert atlas["scenarios"][:-1] == previous["scenarios"]
    scenario = atlas["scenarios"][-1]
    assert scenario["theory"] == previous["scenarios"][-1]["theory"]
    assert atlas["fitted_observations"] == []
    assert response["schema"] == tool.response_model.SCHEMA
    assert tool.first.atlas.validate_programme(programme)
    before = {row["id"]: row for row in previous["programme"]["calculation_graph"]}
    changed = {row["id"] for row in programme["calculation_graph"] if row != before[row["id"]]}
    assert changed == {"wall-trace", "observer-response"}


def test_evaluation_does_not_promote_excluded_or_unsolved_physics(records):
    tool, (atlas, programme, _, _) = records
    scenario = atlas["scenarios"][-1]
    assessment = atlas["assessments"][scenario["id"]]
    assert assessment["wall_waveform_evaluated"] and assessment["observer_waveform_evaluated"]
    assert not assessment["right_parent_access_changed"]
    assert not assessment["permanent_memory_established"]
    assert not assessment["theory_or_preparation_retuned"]
    assert tool.first.atlas.earliest_arrival(scenario, "P", "C", 0., clock="tau") is None
    assert not scenario["connection_records"][0]["remote_preparation_established"]
    assert not scenario["connection_records"][1]["physical_inheritance_established"]
    assert not scenario["connection_records"][2]["equivalence_established"]
    assert programme["completed_calculation"]["future_calculations_require_admission"]
    assert "not evaluated" not in json.dumps(scenario)


def test_waveform_figure_and_narrative_expose_diagnostic_and_domain_limits(records):
    tool, (_, _, prepared, response) = records
    outputs = tool.waveform_figure(response)
    assert set(outputs) == {"prepared-response.svg", "prepared-response.png"}
    assert ET.fromstring(outputs["prepared-response.svg"]).tag.endswith("svg")
    assert outputs["prepared-response.png"].startswith(b"\x89PNG\r\n\x1a\n")
    text = tool.narrative(prepared, response)
    assert "diagnostic, not a second physical" in text
    assert "does not establish permanent memory" in text
    assert "no zero tail or continuation" in text
    assert all(f"| {row['source_coordinate']:.4g} |" in text
               for row in response["analytic_tail"]["independent_controls"])


def test_creation_replay_tampering_and_source_drift(records, tmp_path, monkeypatch):
    tool, record = records
    # Physics/refinement is tested by its owner; reuse this one fixed evaluation
    # while testing artifact reproduction and rejection paths independently.
    monkeypatch.setattr(tool, "build_record", lambda: deepcopy(record))
    output = tmp_path / "atlas-v8"
    assert tool.create(output)["outputs"] == 7
    assert tool.check(output)["outputs"] == 7
    manifest = json.loads((output / "manifest.json").read_bytes())
    assert tool.OWNER in manifest["source_bindings"]
    assert "src/recursive_horizons/nsc_atlas_prepared_response.py" in manifest["source_bindings"]
    assert "paper/nested-space-cosmology.md" not in manifest["source_bindings"]
    assert "paper/nested-space-cosmology.pdf" not in manifest["source_bindings"]
    with pytest.raises(FileExistsError):
        tool.create(output)
    path = output / "prepared-response.json"
    original = path.read_bytes()
    path.write_bytes(original + b"\n")
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.check(output)
    path.write_bytes(original)
    bindings = tool.source_bindings()
    monkeypatch.setattr(tool, "source_bindings", lambda: {**bindings, "changed.py": "0"*64})
    with pytest.raises(ValueError, match="source binding mismatch"):
        tool.check(output)
