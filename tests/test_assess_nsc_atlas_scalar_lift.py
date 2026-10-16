"""Independent polynomial and saved-data provenance checks; no physical solve."""
from copy import deepcopy
import json

import numpy as np
import pytest

import assess_nsc_atlas_scalar_lift as tool


def samples():
    polynomial = np.polynomial.Polynomial([.2, -.3, .7, -.1, .2, -.9, .2, .4])
    derivative = polynomial.deriv()
    # A narrow, offset U domain tests the normalization and slope rescaling.
    def sample(t):
        return {"t": t, "U": -.301+.023*t,
                "b": float(polynomial(t)), "b_U": float(derivative(t)/.023)}
    training = [sample(t) for t in tool.TRAINING_TIMES]
    withheld = [{"t": row["t"], "U": row["U"], "b_reference": row["b"],
                 "b_U_reference": row["b_U"], "b_defect": .1}
                for row in map(sample, tool.WITHHELD_TIMES)]
    return training, withheld


def fixture_record():
    training, withheld = samples()
    return {"comparison": "lift", "solutions": [{"knot_count": 17,
        "knot_samples": training, "validation_samples": withheld}],
        "diagnostics": {"interpolation_target": {"absolute_b_target": 1e-7}}}


@pytest.mark.parametrize("count", tool.NODE_COUNTS)
def test_degree_seven_reproduction_and_normalized_derivative(count):
    training, withheld = samples()
    for query in withheld:
        row = tool.local_prediction(training, query, count)
        assert row["b_interpolated"] == pytest.approx(query["b_reference"], abs=4e-12)
        assert row["cell_width"]*row["b_U_interpolated"] == pytest.approx(
            row["cell_width"]*query["b_U_reference"], abs=4e-12)
        assert row["extrapolated"] is False


def test_source_blocks_cap_available_nodes_and_never_train_on_holdouts():
    record = fixture_record()
    assessment = tool.assess(record)
    for method in assessment["methods"]:
        for row in method["withheld_samples"]:
            lo, hi = tool.BLOCKS[row["block"]]
            assert all(lo <= t <= hi for t in row["training_times"])
            assert not set(row["training_times"]) & set(tool.WITHHELD_TIMES)
            assert row["nodes"] == min(method["requested_nodes"], 9 if row["block"] == 1 else 5)
        assert sum(block["withheld_points"] for block in method["blocks"]) == 16
    # Alter validation values and slopes: predictions and chosen stencils stay fixed.
    changed = deepcopy(record)
    for row in changed["solutions"][0]["validation_samples"]:
        row["b_reference"] += 100.
        row["b_U_reference"] -= 1000.
    other = tool.assess(changed)
    for left, right in zip(assessment["methods"], other["methods"]):
        for a, b in zip(left["withheld_samples"], right["withheld_samples"]):
            for key in ("b_interpolated", "b_U_interpolated", "training_indices"):
                assert a[key] == b[key]


def test_cell_width_is_not_stencil_normalization_width():
    training, withheld = samples()
    row = tool.local_prediction(training, withheld[7], 9)
    assert row["normalization_width"] == pytest.approx(8*row["cell_width"])


@pytest.mark.parametrize("mutation", ("holdout", "edge", "extrapolate", "mismatched_cell", "unordered", "nonfinite"))
def test_reject_invalid_training_or_query(mutation):
    training, withheld = samples()
    query = dict(withheld[1])
    if mutation == "holdout":
        training[1]["t"] = tool.WITHHELD_TIMES[0]
    elif mutation == "edge":
        query["t"] = .25
    elif mutation == "extrapolate":
        query["U"] = 0.
    elif mutation == "mismatched_cell":
        query["U"] = withheld[2]["U"]
    elif mutation == "unordered":
        training[1]["U"] = training[0]["U"]
    else:
        training[1]["b_U"] = float("nan")
    with pytest.raises(ValueError):
        tool.local_prediction(training, query, 4)


def test_source_authentication_is_reused_and_rejection_propagates(tmp_path, monkeypatch):
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    directory = tmp_path/tool.INPUT
    directory.mkdir(parents=True)
    raw = tool.encoded(fixture_record())
    (directory/"inputs.json").write_bytes(raw)
    manifest = {"outputs": {"inputs.json": {"sha256": tool.digest(raw), "bytes": len(raw)}},
                "source_bindings": {"retained_source.py": "saved-hash"}}
    (directory/"manifest.json").write_bytes(tool.encoded(manifest))
    calls = []
    def authenticate(path):
        calls.append(path)
        return {"comparison": "lift", "physical_solves": 0}
    monkeypatch.setattr(tool.predecessor, "check", authenticate)
    record, binding = tool.input_snapshot()
    assert calls == [directory]
    assert record == fixture_record()
    assert binding["source_validation"]["source_bindings_sha256"] == tool.digest(tool.encoded(manifest["source_bindings"]))
    def reject(_):
        raise ValueError("predecessor source authentication failed")
    monkeypatch.setattr(tool.predecessor, "check", reject)
    with pytest.raises(ValueError, match="predecessor"):
        tool.input_snapshot()


@pytest.mark.parametrize("mutation", ("prediction", "scope", "runtime"))
def test_creation_only_and_authenticated_postprocessing_recheck(tmp_path, monkeypatch, mutation):
    monkeypatch.setattr(tool, "input_snapshot", lambda: (fixture_record(), {"authenticated": True}))
    monkeypatch.setattr(tool, "source_bindings", lambda: {"script-and-tests": "test-hash"})
    directory = tmp_path/"assessment"
    created = tool.create(directory)
    assert created["execution"]["scalar_kernel_values"] == 0
    assert tool.check(directory)["physical_solves"] == 0
    with pytest.raises(FileExistsError):
        tool.create(directory)
    record = json.loads((directory/"assessment.json").read_bytes())
    if mutation == "prediction":
        record["assessment"]["methods"][0]["withheld_samples"][0]["b_interpolated"] += 1.
    elif mutation == "scope":
        record["scope"] = "unsupported physical result"
    else:
        record["runtime"]["python"] = "unrecorded"
    raw = tool.encoded(record)
    (directory/"assessment.json").write_bytes(raw)
    manifest = json.loads((directory/"manifest.json").read_bytes())
    manifest["outputs"]["assessment.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
    (directory/"manifest.json").write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError, match="deterministic|runtime"):
        tool.check(directory)


def test_cpu_budget_is_enforced(monkeypatch):
    monkeypatch.setattr(tool.time, "process_time", lambda: 31.)
    with pytest.raises(RuntimeError, match="30-process-CPU"):
        tool._budget(0.)


def test_recomputed_assessment_accepts_ulp_but_preserves_integer_metadata(tmp_path, monkeypatch):
    monkeypatch.setattr(tool, "input_snapshot", lambda: (fixture_record(), {"authenticated": True}))
    monkeypatch.setattr(tool, "source_bindings", lambda: {"script-and-tests": "test-hash"})
    directory = tmp_path/"assessment"
    tool.create(directory)
    payload = directory/"assessment.json"
    manifest_path = directory/"manifest.json"
    record = json.loads(payload.read_bytes())
    row = record["assessment"]["methods"][0]["withheld_samples"][0]
    row["b_interpolated"] = float(np.nextafter(row["b_interpolated"], np.inf))
    def save():
        raw = tool.encoded(record)
        payload.write_bytes(raw)
        manifest = json.loads(manifest_path.read_bytes())
        manifest["outputs"]["assessment.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
        manifest_path.write_bytes(tool.encoded(manifest))
    save()
    assert tool.check(directory)["physical_solves"] == 0
    row["nodes"] = float(row["nodes"])
    save()
    with pytest.raises(ValueError, match="deterministic"):
        tool.check(directory)
