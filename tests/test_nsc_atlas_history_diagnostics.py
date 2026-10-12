"""Failure diagnostics preserve original checkers and expose JSON differences."""
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import diagnose_history_replay as diagnostics


def test_json_diagnostics_use_original_numeric_tolerance_and_paths():
    stored = {"budget": {"rhs_calls": [10, 20]}, "checks": {"small": 1e-10},
              "samples": [1.0, 2.0], "flag": False}
    regenerated = {"budget": {"rhs_calls": [11, 20]}, "checks": {"small": 1.2e-10},
                   "samples": [1.0+1e-11, 2.0], "flag": False}
    result = diagnostics.json_mismatches(stored, regenerated)
    assert result["comparison"] == {"relative": 1e-10, "absolute": 1e-12}
    assert result["mismatch_count"] == 2
    assert [row["path"] for row in result["reported"]] == ["$.budget.rhs_calls[0]", "$.checks.small"]
    assert result["reported"][1]["absolute_difference"] > result["reported"][1]["allowed_difference"]


def test_json_diagnostics_are_bounded_and_include_structural_differences():
    result = diagnostics.json_mismatches([0]*30, [1]*30, limit=3)
    assert result["mismatch_count"] == 30
    assert len(result["reported"]) == 3 and result["truncated"]
    assert diagnostics.json_mismatches({"a": 1}, {"b": 1})["reported"][0]["reason"] == "key_inventory"
    assert diagnostics.json_mismatches([], [1])["reported"][0]["reason"] == "list_length"
    assert diagnostics.json_mismatches("[]", [1])["reported"][0]["reason"] == "type"
    assert diagnostics.json_mismatches(False, True)["mismatch_count"] == 1


@pytest.mark.parametrize("atlas_version", [1, 5, 6, 7, 8])
def test_explicit_original_context_selection_preserves_legacy_route(tmp_path, monkeypatch, atlas_version):
    replay = diagnostics.replay
    monkeypatch.setattr(replay, "ROOT", tmp_path)
    monkeypatch.setattr(replay, "CONTEXT", tmp_path / "legacy")
    monkeypatch.setattr(replay, "V7_CONTEXT", tmp_path / "v7")
    monkeypatch.setattr(replay, "V8_CONTEXT", tmp_path / "v8")
    calls = []
    for number, name in ((0, "prepare"), (7, "prepare_v7"), (8, "prepare_v8")):
        def prepare(*args, number=number):
            calls.append((number, args))
            return {"source_commit": f"commit-{number}"}
        monkeypatch.setattr(replay, name, prepare)
    selected = diagnostics.target(atlas_version)
    suffix = "" if atlas_version == 1 else f"_v{atlas_version}"
    assert selected["source"] == f"scripts/derive_nsc_cosmological_atlas{suffix}.py"
    if atlas_version < 7:
        assert calls == [(0, ())]
        assert selected["owner"] == tmp_path / "legacy/lab"
        assert selected["bundle"] == tmp_path / "legacy" / replay.BUNDLES[atlas_version-1]
    else:
        context = tmp_path / f"v{atlas_version}"
        assert calls == [(atlas_version, (context,))]
        assert selected["owner"] == context
        assert selected["bundle"] == context / f"results/development/nsc-cosmological-atlas-v{atlas_version}"
        assert selected["environment"]["PYTHONNOUSERSITE"] == "1"
        assert selected["environment"]["GIT_ALTERNATE_OBJECT_DIRECTORIES"] == str(tmp_path / "archive/source-history/objects")


@pytest.mark.parametrize("returncode", [0, 1])
def test_diagnosis_is_creation_only_keeps_outputs_and_never_accepts(tmp_path, monkeypatch, returncode):
    original = tmp_path / "original"
    original.mkdir()
    (original / "manifest.json").write_text(json.dumps({
        "outputs": {"prepared-response.json": {}, "plot.png": {}},
        "numerical_replay_tolerance": {"relative": 1e-10, "absolute": 1e-12},
    }), encoding="utf-8")
    (original / "prepared-response.json").write_text('{"checks":{"small":1e-10}}', encoding="utf-8")
    (original / "plot.png").write_bytes(b"original figure")
    chosen = {"source_commit": "original-commit", "bundle": original,
              "source": "scripts/original-producer.py", "owner": tmp_path / "original-context",
              "environment": {"PYTHONPATH": "original-context/src"}}
    monkeypatch.setattr(diagnostics, "target", lambda version: chosen)
    monkeypatch.setattr(diagnostics, "runtime_metadata", lambda: {"processor": "test CPU"})
    calls = []
    def run(command, **kwargs):
        calls.append((command, kwargs))
        output = Path(command[-1])
        output.mkdir()
        (output / "prepared-response.json").write_text('{"checks":{"small":2e-10}}', encoding="utf-8")
        if not returncode:
            (output / "plot.png").write_bytes(b"different figure")
        return SimpleNamespace(returncode=returncode, stdout="producer output", stderr="producer diagnostic")
    monkeypatch.setattr(diagnostics.subprocess, "run", run)
    destination = tmp_path / "diagnostic"
    report = diagnostics.diagnose(8, destination)
    assert report["accepted"] is False and report["returncode"] == returncode
    assert report["source_commit"] == "original-commit"
    assert report["outputs"]["prepared-response.json"]["json_comparison"]["reported"][0]["path"] == "$.checks.small"
    assert (destination / "runtime.json").is_file()
    assert json.loads((destination / "runtime.json").read_text())["processor"] == "test CPU"
    assert (original / "prepared-response.json").read_text() == '{"checks":{"small":1e-10}}'
    if returncode:
        assert report["outputs"]["plot.png"]["missing_regenerated_output"]
    else:
        assert report["outputs"]["plot.png"]["byte_equal"] is False
    assert calls[0][0] == [diagnostics.sys.executable, "scripts/original-producer.py", "--output", str(destination / "regenerated")]
    assert calls[0][1]["cwd"] == chosen["owner"]
    assert calls[0][1]["env"] == chosen["environment"]
    with pytest.raises(FileExistsError):
        diagnostics.diagnose(8, destination)
    assert len(calls) == 1


def test_original_context_failure_never_starts_producer(tmp_path, monkeypatch):
    def fail(version):
        raise ValueError("historical Git object unavailable")
    monkeypatch.setattr(diagnostics, "target", fail)
    monkeypatch.setattr(diagnostics.subprocess, "run", lambda *a, **k: pytest.fail("must not start"))
    with pytest.raises(ValueError, match="object unavailable"):
        diagnostics.diagnose(8, tmp_path / "diagnostic")
    assert not (tmp_path / "diagnostic").exists()
