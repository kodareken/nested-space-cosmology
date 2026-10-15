"""Exercise evidence preservation separately from the scientific identities."""
import json
from pathlib import Path
import shutil

import pytest

import record_nsc_atlas_cap_wall_response as tool


@pytest.fixture
def isolated_record(tmp_path, monkeypatch):
    original = tool.ROOT
    root = tmp_path / "repository"
    for relative in tool.SOURCES:
        target = root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(original / relative, target)
    monkeypatch.setattr(tool, "ROOT", root)
    # This fixture tests the new record's contract. The real predecessor
    # chain is authenticated by the production creation/check commands.
    monkeypatch.setattr(tool, "authenticate_predecessor", lambda: {"fixture": True})
    monkeypatch.setattr(tool, "algebra", tool.expected_certificate)
    return root / tool.RELATIVE_OUTPUT


def test_create_check_and_preserve_existing_record(isolated_record, monkeypatch):
    assert tool.create(isolated_record)["outputs"] == 2
    before = {p.name: p.read_bytes() for p in isolated_record.iterdir()}
    monkeypatch.setattr(tool, "algebra", lambda: pytest.fail("unexpected symbolic replay"))
    assert tool.check(isolated_record)["numerical_integrations"] == 0
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(isolated_record)
    assert before == {p.name: p.read_bytes() for p in isolated_record.iterdir()}


def test_source_and_output_drift_detected(isolated_record):
    tool.create(isolated_record)
    source = tool.ROOT / tool.OWNER
    before = source.read_bytes()
    source.write_bytes(before + b"\n")
    with pytest.raises(ValueError, match="source binding mismatch"):
        tool.check(isolated_record)
    source.write_bytes(before)
    output = isolated_record / "cap-wall.json"
    output.write_bytes(output.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.check(isolated_record)


@pytest.mark.parametrize("alteration", ["declaration", "certificate", "model"])
def test_rehashed_alteration_rejected(isolated_record, alteration):
    tool.create(isolated_record)
    name = "model-specification.json" if alteration == "model" else "cap-wall.json"
    path = isolated_record / name
    record = json.loads(path.read_bytes())
    if alteration == "model":
        record["id"] = "different-model"
    elif alteration == "certificate":
        record["exact_algebra"]["zero_residuals"] = {}
    else:
        record["status"] = "FULL_EVOLUTION_SOLVED"
    raw = tool.canonical_bytes(record)
    path.write_bytes(raw)
    manifest_path = isolated_record / "manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["outputs"][name] = {"sha256": tool.sha(raw), "bytes": len(raw)}
    manifest_path.write_bytes(tool.canonical_bytes(manifest))
    with pytest.raises(ValueError, match="declaration mismatch|model specification mismatch"):
        tool.check(isolated_record)


def test_explicit_algebra_authenticates_first(isolated_record, monkeypatch):
    tool.create(isolated_record)
    assert tool.check_algebra(isolated_record)["replayed"] == isolated_record.name
    path = isolated_record / "model-specification.json"
    path.write_bytes(path.read_bytes() + b"tamper")
    monkeypatch.setattr(tool, "algebra", lambda: pytest.fail("replay before authentication"))
    with pytest.raises(ValueError, match="artifact binding mismatch"):
        tool.check_algebra(isolated_record)


def test_bad_certificate_or_budget_does_not_create_evidence(isolated_record, monkeypatch):
    with monkeypatch.context() as patch:
        patch.setattr(tool, "algebra", lambda: {})
        with pytest.raises(ValueError, match="certificate mismatch"):
            tool.create(isolated_record)
    assert not isolated_record.exists()
    ticks = iter((0.0, 61.0))
    monkeypatch.setattr(tool.time, "process_time", lambda: next(ticks))
    with pytest.raises(RuntimeError, match="process CPU budget"):
        tool.create(isolated_record)
    assert not isolated_record.exists()
