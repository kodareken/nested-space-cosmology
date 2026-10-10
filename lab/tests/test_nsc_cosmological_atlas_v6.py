import json

import pytest

from test_nsc_cosmological_atlas_v5 import producer as predecessor


def producer():
    predecessor()
    import derive_nsc_cosmological_atlas_v6
    return derive_nsc_cosmological_atlas_v6


def test_field_successor_keeps_the_prior_routes_and_causal_exclusions():
    tool = producer()
    before, _, _ = tool.previous.build_record()
    after, response = tool.build_record()
    assert after["scenarios"][:-1] == before["scenarios"]
    assert len(after["scenarios"]) == 9
    assert after["fitted_observations"] == []
    assert response["incoming_state"]["prescribed_characteristics"] == 1
    assert not after["assessments"]["early_conformal_field_response"]["remote_parent_source_solved"]


def test_field_response_does_not_invent_a_global_state_or_unconditional_delay():
    tool = producer()
    record, response = tool.build_record()
    scenario = record["scenarios"][-1]
    assert not scenario["theory"]["backreaction"]
    assert not scenario["theory"]["added_surface_interaction"]
    assert "not supplied" in next(r for r in scenario["regions"] if r["id"] == "P")["state"]
    assert all(e["clock"] is None and e["minimum_delay"] is None for e in scenario["links"])
    assert tool.first.atlas.earliest_arrival(scenario, "W", "C", 0., clock="tau") is None
    assert not response["metrics"]["squared_amplitude_is_energy"]
    assert response["support"]["centre_arrival_t_over_ell"][1] < 0


def test_response_bundle_is_source_bound_creation_only_and_replayable(tmp_path):
    tool = producer()
    output = tmp_path/"atlas-v6"
    assert tool.create(output)["outputs"] == 7
    assert tool.check(output)["outputs"] == 7
    manifest = json.loads((output/"manifest.json").read_bytes())
    assert "field-response.png" in manifest["outputs"]
    assert "field-transfer.json" in manifest["outputs"]
    assert any(p.endswith("atlas-working-draft-v6.md") for p in manifest["source_bindings"])
    with pytest.raises(FileExistsError):
        tool.create(output)
    payload = output/"field-transfer.json"
    payload.write_bytes(payload.read_bytes()+b"\n")
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.check(output)
