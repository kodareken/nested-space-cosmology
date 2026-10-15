"""Original bulk evidence is explicit, isolated and never repaired in place."""
import io
import json
import sys
import tarfile
from types import ModuleType, SimpleNamespace

import pytest

import replay_history as replay


@pytest.fixture
def bulk_archive(tmp_path, monkeypatch):
    root, context = tmp_path / "live", tmp_path / "original"
    files = {"src/recursive_horizons/__init__.py": b"original package",
             "src/recursive_horizons/nsc_atlas_prepared_response.py": b"original field",
             "scripts/record_nsc_atlas_bulk_response.py": b"original recorder",
             "scripts/run_python.py": b"original launcher", "pyproject.toml": b"original packages"}
    common = {name: replay.sha(raw) for name, raw in files.items()}
    for kind, bundle in zip(replay.BULK_V1_COMPARISONS, replay.BULK_V1_BUNDLES):
        source = "docs/" + kind + ".md"
        files[source] = kind.encode()
        raw = json.dumps({"comparison": kind, "runtime": {"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"}}).encode()
        files[bundle + "/response.json"] = raw
        files[bundle + "/manifest.json"] = json.dumps({
            "source_bindings": {**common, source: replay.sha(files[source])},
            "outputs": {"response.json": {"sha256": replay.sha(raw), "bytes": len(raw)}}}).encode()
    for name, raw in files.items():
        path = root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
    calls = []

    def git(*args):
        calls.append(args)
        if args[0] == "show":
            commit, name = args[1].split(":", 1)
            assert commit == replay.BULK_V1_COMMIT
            return files[name]
        assert args[:3] == ("archive", "--format=tar", replay.BULK_V1_COMMIT)
        stream = io.BytesIO()
        with tarfile.open(fileobj=stream, mode="w") as archive:
            for name, raw in files.items():
                member = tarfile.TarInfo(name)
                member.size = len(raw)
                archive.addfile(member, io.BytesIO(raw))
        return stream.getvalue()

    monkeypatch.setattr(replay, "ROOT", root)
    monkeypatch.setattr(replay, "git", git)
    return SimpleNamespace(root=root, context=context, files=files, calls=calls)


def test_bulk_union_uses_one_exact_snapshot_and_receipt(bulk_archive):
    row = bulk_archive
    owner = row.root / "src/recursive_horizons/nsc_atlas_prepared_response.py"
    owner.write_bytes(b"new live API")
    result = replay.prepare_bulk_v1(row.context)
    assert result["source_commit"] == "3470152e176611849140eec3b02bafede405d147"
    assert result["evidence_authenticated"] == 12
    assert result["scientific_replay"] is False
    receipt = json.loads((row.context / "replay-context.json").read_bytes())
    assert receipt["files"] == {name: replay.sha(raw) for name, raw in row.files.items()}
    archives = [call for call in row.calls if call[0] == "archive"]
    assert len(archives) == 1
    assert {"src", "pyproject.toml", "scripts/run_python.py", *row.files.keys()} <= set(archives[0][3:])
    assert replay.prepare_bulk_v1(row.context) == result
    assert owner.read_bytes() == b"new live API"


@pytest.mark.parametrize("changed", ["live-output", "live-manifest", "context-source", "receipt", "source-binding"])
def test_bulk_drift_fails_without_replacement(bulk_archive, changed):
    row = bulk_archive
    replay.prepare_bulk_v1(row.context)
    bundle = replay.BULK_V1_BUNDLES[-1]
    if changed == "source-binding":
        name = bundle + "/manifest.json"
        manifest = json.loads(row.files[name])
        manifest["source_bindings"]["docs/time.md"] = "0" * 64
        row.files[name] = json.dumps(manifest).encode()
        (row.root / name).write_bytes(row.files[name])
        context = row.context.parent / "fresh-context"
        with pytest.raises(ValueError, match="original source does not match"):
            replay.prepare_bulk_v1(context)
        assert not (context / "replay-context.json").exists()
        return
    path = {"live-output": row.root / bundle / "response.json",
            "live-manifest": row.root / bundle / "manifest.json",
            "context-source": row.context / "docs/metric.md",
            "receipt": row.context / "replay-context.json"}[changed]
    path.write_bytes(b"altered")
    with pytest.raises(ValueError, match="evidence changed|context changed|receipt changed"):
        replay.prepare_bulk_v1(row.context)
    assert path.read_bytes() == b"altered"


@pytest.mark.parametrize("failed", [False, True])
def test_bulk_authentication_runs_only_original_checks(bulk_archive, monkeypatch, failed):
    row, calls = bulk_archive, []
    monkeypatch.setenv("PYTHONPATH", "live/src" + replay.os.pathsep + "deps/site-packages")
    def run(command, **options):
        calls.append((command, options))
        return SimpleNamespace(returncode=int(failed), stdout='{"checks":[]}', stderr="original source mismatch" if failed else "")
    monkeypatch.setattr(replay.subprocess, "run", run)
    if failed:
        with pytest.raises(ValueError, match="original source mismatch"):
            replay.check_bulk_v1(row.context)
    else:
        result = replay.check_bulk_v1(row.context)
        assert result["metric_solves"] == 0 and not result["scientific_replay"]
    assert len(calls) == 1
    command, options = calls[0]
    assert command[:2] == [replay.sys.executable, "-c"]
    assert "original.check" in command[2] and "original.compare" not in command[2]
    assert tuple(command[3:]) == replay.BULK_V1_BUNDLES
    assert options["cwd"] == row.context
    assert "live/src" not in options["env"]["PYTHONPATH"]
    assert str(row.context / "src") in options["env"]["PYTHONPATH"]
    assert options["env"]["PYTHONNOUSERSITE"] == "1"


def test_bulk_numerical_replay_is_explicit_and_profile_bound(bulk_archive, monkeypatch):
    row, calls = bulk_archive, []
    replay.prepare_bulk_v1(row.context)
    monkeypatch.setattr(replay, "check_bulk_v1", lambda context: {"scientific_replay": False})
    monkeypatch.setattr(replay.platform, "python_version", lambda: "3.14.8")
    monkeypatch.setattr(replay, "version", lambda name: "wrong")
    monkeypatch.setattr(replay, "_bulk_command", lambda *args: calls.append(args) or {"checks": []})
    with pytest.raises(ValueError, match="check-bulk-v1.*no fallback"):
        replay.run_bulk_v1(row.context)
    assert calls == []
    monkeypatch.setattr(replay, "version", lambda name: {"numpy": "2.5.1", "scipy": "1.17.1"}[name])
    result = replay.run_bulk_v1(row.context)
    assert result["scientific_replay"] and result["metric_solves"] == 12
    assert "original.compare" in calls[0][1]
    assert "original.validate_record" in calls[0][1]
    assert "saved.pop('execution')" in calls[0][1]


def test_explicit_replay_compares_science_but_not_cpu_timing(tmp_path, monkeypatch):
    """Execute the orchestration with a synthetic producer, never physical code."""
    bundle = tmp_path / "synthetic"
    bundle.mkdir()
    (bundle / "response.json").write_text(json.dumps({"comparison": "metric", "coefficient": .2,
                                                      "execution": {"cpu_seconds": 1.}}))
    original = ModuleType("record_nsc_atlas_bulk_response")
    original.check = lambda directory: {"checked": str(directory)}
    original.validate_record = lambda value, kind: True
    original.compare = lambda kind: {"comparison": kind, "coefficient": .2,
                                      "execution": {"cpu_seconds": 2.}}
    monkeypatch.setitem(sys.modules, original.__name__, original)
    monkeypatch.setattr(replay, "check_bulk_v1", lambda context: {})
    monkeypatch.setattr(replay, "require_bulk_v1_profile", lambda context: {})
    code = []
    monkeypatch.setattr(replay, "_bulk_command", lambda context, text: code.append(text) or {})
    replay.run_bulk_v1(tmp_path)
    monkeypatch.setattr(sys, "argv", ["synthetic", str(bundle)])
    exec(compile(code[0], "bulk-replay-test", "exec"), {})
    original.compare = lambda kind: {"comparison": kind, "coefficient": .3,
                                      "execution": {"cpu_seconds": 2.}}
    with pytest.raises(ValueError, match="numerical replay mismatch"):
        exec(compile(code[0], "bulk-replay-test", "exec"), {})


def test_publication_selects_named_original_context_without_fallback(monkeypatch):
    import check_publication
    calls = []
    monkeypatch.setattr(replay, "check_bulk_v1", lambda: calls.append("bulk-v1"))
    errors = []
    check_publication.check_bulk_response_history(errors)
    assert errors == [] and calls == ["bulk-v1"]
    def fail():
        raise ValueError("changed original context")
    monkeypatch.setattr(replay, "check_bulk_v1", fail)
    check_publication.check_bulk_response_history(errors)
    assert errors == ["bulk response original-context evidence: changed original context"]
