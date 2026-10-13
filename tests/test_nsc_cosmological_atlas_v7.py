"""The finite IV preparation advances one connection without replacing other routes."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import pytest

def producer():
    import derive_nsc_cosmological_atlas_v7
    return derive_nsc_cosmological_atlas_v7


def test_preparation_successor_keeps_all_nine_earlier_scenarios_and_open_questions():
    tool = producer()
    before, _ = tool.previous.build_record()
    after, programme, prepared = tool.build_record()
    assert after["scenarios"][:-1] == before["scenarios"]
    assert len(after["scenarios"]) == 10
    assert after["fitted_observations"] == []
    assert after["programme"] == programme
    assert tool.first.atlas.validate_programme(programme)
    assert {q["id"] for q in programme["questions"]} == set(tool.first.atlas.PROGRAMME_THEMES)
    assert programme["new_compute_admitted"] is False
    assert prepared["operator"]["analytic_norm_bound"] == pytest.approx(.13882099, abs=1e-8)


def test_prepared_operator_does_not_assert_a_waveform_or_reopen_parent_access():
    tool = producer()
    record, _, _ = tool.build_record()
    scenario = record["scenarios"][-1]
    findings = record["assessments"][scenario["id"]]
    assert findings["transfer_operator_certified"]
    assert not findings["right_parent_access_changed"]
    assert not findings["wall_waveform_evaluated"]
    assert not findings["observer_waveform_evaluated"]
    assert not findings["global_preparation_or_quantum_state_solved"]
    assert scenario["theory"]["curvature_coupling"] == 1/6
    assert not scenario["theory"]["backreaction"]
    current = json.dumps(scenario)
    assert "field-transfer.json" not in current
    assert "R/ell=0.45" in current
    assert "fixed incoming history, including zero outside its pulse" not in scenario["assumptions"]
    assert "operator" in scenario["links"][0]["evidence"]
    assert not any(link["from"] == "P" for link in scenario["links"])
    assert tool.first.atlas.earliest_arrival(scenario, "P", "C", 0., clock="tau") is None
    claims = scenario["connection_records"]
    assert not claims[0]["remote_preparation_established"]
    assert not claims[1]["physical_inheritance_established"]
    assert not claims[2]["equivalence_established"]


def test_preparation_figures_plot_domain_input_and_research_dependencies_only():
    tool = producer()
    record, programme, prepared = tool.build_record()
    outputs = tool.preparation_figure(prepared) | tool.dependency_figure(programme)
    assert set(outputs) == {"preparation-domain.png", "preparation-domain.svg", "research-dependencies.png", "research-dependencies.svg"}
    for name, data in outputs.items():
        if name.endswith(".svg"):
            assert ET.fromstring(data).tag.endswith("svg")
        else:
            assert data.startswith(b"\x89PNG\r\n\x1a\n")
    assert "waveform has not been evaluated" in tool.narrative(prepared)
    assert "causal exclusion is unchanged" in tool.narrative(prepared)


def test_successor_is_creation_only_source_bound_and_replayable(tmp_path, monkeypatch):
    tool = producer()
    output = tmp_path / "atlas-v7"
    assert tool.create(output)["outputs"] == 8
    assert tool.check(output)["outputs"] == 8
    manifest = json.loads((output / "manifest.json").read_bytes())
    assert "paper/nested-space-cosmology.md" not in manifest["source_bindings"]
    assert "paper/nested-space-cosmology.pdf" not in manifest["source_bindings"]
    assert "docs/nsc-atlas-field-preparation.md" in manifest["source_bindings"]
    assert all(not Path(path).is_absolute() and "\\" not in path for path in manifest["source_bindings"])
    with pytest.raises(FileExistsError):
        tool.create(output)
    payload = output / "preparation-access.json"
    data = payload.read_bytes()
    payload.write_bytes(data + b"\n")
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.check(output)
    payload.write_bytes(data)
    before = tool.source_bindings()
    monkeypatch.setattr(tool, "source_bindings", lambda: {**before, "scripts/changed.py": "0" * 64})
    with pytest.raises(ValueError, match="source binding mismatch"):
        tool.check(output)
