"""Calibration receipts reject drift without evolving or rewriting history."""
import json
import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import derive_nsc_ks_evaluator_accuracy_v2 as A


def test_exclusive_record_preserves_existing_bytes_and_foreign_staging(tmp_path):
    path = tmp_path / "record.json"
    A.write_once(path, {"verdict": "OPEN"})
    original = path.read_bytes()
    with pytest.raises(FileExistsError):
        A.write_once(path, {"verdict": "EXISTENCE"})
    assert path.read_bytes() == original
    stage = path.with_name(path.name + f".{os.getpid()}.tmp")
    stage.write_bytes(b"foreign checkpoint")
    with pytest.raises(FileExistsError):
        A.write_once(path, {"verdict": "OPEN"})
    assert stage.read_bytes() == b"foreign checkpoint"


def fixture(tmp_path, monkeypatch):
    monkeypatch.setattr(A, "ROOT", tmp_path)
    output, directory = A.paths("control")
    source = tmp_path / "source.py"
    source.write_bytes(b"source version one\n")
    binding = {
        "profile_identity": "a" * 64, "history_path": "history.json",
        "target": [0.12, 0.18], "solver": {"integrator": "cf4", "max_step": 0.001},
        "family_keys": [[14, 1]], "source_hashes": {"source.py": A.digest(source)},
    }
    binding_path = directory / "binding.json"
    A.write_once(binding_path, binding)
    payload = directory / "family.npz"
    with payload.open("wb") as stream:
        np.savez(stream, matter_change=np.array([[1e-4, -2e-4], [3e-5, 1e-4]]))
    result = {
        **binding, "schema": A.SCHEMA, "complete": True,
        "binding": A.descriptor(binding_path),
        "families": [{"family": [14, 1], "payload": A.descriptor(payload),
                      "binding_sha256": A.digest(binding_path),
                      "matter_change_maxima": [1e-4, 2e-4]}],
    }
    A.write_once(output, result)
    return output, binding_path, payload, source


def test_replay_authenticates_source_and_reconstructs_saved_maxima(tmp_path, monkeypatch):
    _output, _binding, _payload, source = fixture(tmp_path, monkeypatch)
    assert A.load("control")["families"][0]["matter_change_maxima"] == [1e-4, 2e-4]
    source.write_bytes(b"source version two\n")
    with pytest.raises(ValueError, match="source changed"):
        A.load("control")


@pytest.mark.parametrize("target", ["binding", "payload"])
def test_replay_rejects_changed_descendant(tmp_path, monkeypatch, target):
    _output, binding, payload, _source = fixture(tmp_path, monkeypatch)
    path = binding if target == "binding" else payload
    path.write_bytes(path.read_bytes() + b" ")
    with pytest.raises(ValueError, match="dependency changed"):
        A.load("control")


def test_replay_rejects_header_and_family_coverage_changes(tmp_path, monkeypatch):
    output, _binding, _payload, _source = fixture(tmp_path, monkeypatch)
    saved = json.loads(output.read_text())
    changed = {**saved, "solver": {"integrator": "dop853", "max_step": 0.001}}
    output.write_text(json.dumps(changed))
    with pytest.raises(ValueError, match="header/binding differs"):
        A.load("control")
    saved["families"][0]["family"] = [32, -1]
    output.write_text(json.dumps(saved))
    with pytest.raises(ValueError, match="coverage differs"):
        A.load("control")


def test_comparison_refuses_different_history_before_field_work(monkeypatch):
    records = {
        "left": {"profile_identity": "a"},
        "right": {"profile_identity": "b"},
    }
    monkeypatch.setattr(A, "load", lambda name: records[name])
    with pytest.raises(ValueError, match="physical input or code"):
        A.compare("left", "right")
