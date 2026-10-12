import io
import json
import tarfile
from pathlib import Path
from types import SimpleNamespace

import pytest

import replay_history as replay


def payload(name, data=b"original"):
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode="w") as archive:
        member = tarfile.TarInfo(name)
        member.size = len(data)
        archive.addfile(member, io.BytesIO(data))
    return stream.getvalue()


def test_context_drift_is_never_silently_repaired(tmp_path):
    archive = payload("scripts/old.py")
    replay.install_archive(archive, tmp_path)
    (tmp_path / "scripts/old.py").write_bytes(b"changed")
    with pytest.raises(ValueError, match="context changed"):
        replay.install_archive(archive, tmp_path)


def test_current_evidence_tampering_cannot_be_hidden_by_git(tmp_path, monkeypatch):
    current = tmp_path / "current"
    current.write_bytes(b"changed")
    monkeypatch.setattr(replay, "current_artifact", lambda name: current)
    with pytest.raises(ValueError, match="evidence changed"):
        replay.install_archive(payload("lab/results/evidence.json"), tmp_path / "context")


def test_corrupt_archive_and_escaping_member_fail(tmp_path):
    with pytest.raises(tarfile.ReadError):
        replay.install_archive(b"corrupt", tmp_path)
    with pytest.raises(ValueError, match="escapes"):
        replay.install_archive(payload("../escape"), tmp_path)


def test_missing_git_object_has_no_live_source_fallback(monkeypatch):
    class Missing:
        returncode = 128
        stderr = b"missing object"
    monkeypatch.setattr(replay.subprocess, "run", lambda *a, **k: Missing())
    with pytest.raises(ValueError, match="object unavailable"):
        replay.git("archive", replay.BASELINE)


def test_historical_imports_exclude_live_sources(tmp_path, monkeypatch):
    monkeypatch.setenv("PYTHONPATH", str(Path("live/src")) + replay.os.pathsep + "deps/site-packages")
    env = replay.environment(tmp_path, True)
    assert "live/src" not in env["PYTHONPATH"]
    assert str(tmp_path / "lab/src") in env["PYTHONPATH"]


@pytest.fixture
def v7_archive(tmp_path, monkeypatch):
    """Small original-commit archive with an owner allowed to evolve live."""
    root = tmp_path / "current"
    context = tmp_path / "replay-v7"
    files = {"src/recursive_horizons/__init__.py": b"original package",
             "scripts/derive_nsc_cosmological_atlas_v7.py": b"original checker",
             "docs/v7-owner.md": b"original owner",
             "pyproject.toml": b"original dependency profile",
             replay.V7_BUNDLE + "/atlas.json": b'{"original":true}\n'}
    manifest = {"source_bindings": {name: replay.sha(data) for name, data in files.items()
                                    if name.startswith(("src/", "scripts/", "docs/"))},
                "outputs": {"atlas.json": {"sha256": replay.sha(files[replay.V7_BUNDLE + "/atlas.json"])}},
                "runtime": {"python": "3.14.8", "matplotlib": "3.11.1"}}
    files[replay.V7_BUNDLE + "/manifest.json"] = json.dumps(manifest).encode()
    stream = io.BytesIO()
    with tarfile.open(fileobj=stream, mode="w") as archive:
        for name, data in files.items():
            member = tarfile.TarInfo(name)
            member.size = len(data)
            archive.addfile(member, io.BytesIO(data))
            path = root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
    calls = []
    def original_git(*args):
        calls.append(args)
        if args[0] == "show":
            assert args[1] == replay.V7_COMMIT + ":" + replay.V7_BUNDLE + "/manifest.json"
            return files[replay.V7_BUNDLE + "/manifest.json"]
        assert args[:3] == ("archive", "--format=tar", replay.V7_COMMIT)
        return stream.getvalue()
    monkeypatch.setattr(replay, "ROOT", root)
    monkeypatch.setattr(replay, "git", original_git)
    return root, context, files, calls


def test_v7_preparation_uses_fixed_commit_despite_live_owner_changes(v7_archive):
    root, context, files, calls = v7_archive
    (root / "docs/v7-owner.md").write_bytes(b"current scientific successor")
    result = replay.prepare_v7(context)
    assert result["source_commit"] == "7fee894ad5a4da2ba378f724e45e926401f82670"
    assert not result["scientific_replay"]
    assert result["evidence_authenticated"] == 2
    assert (context / "docs/v7-owner.md").read_bytes() == files["docs/v7-owner.md"]
    assert replay.prepare_v7(context) == result
    receipt = json.loads((context / "replay-context.json").read_bytes())
    assert receipt["files"]["docs/v7-owner.md"] == replay.sha(files["docs/v7-owner.md"])
    assert all(call[0] in {"show", "archive"} for call in calls)


def test_v7_retained_evidence_and_frozen_context_cannot_be_replaced(v7_archive):
    root, context, files, _ = v7_archive
    replay.prepare_v7(context)
    artifact = root / replay.V7_BUNDLE / "atlas.json"
    artifact.write_bytes(b"changed evidence")
    with pytest.raises(ValueError, match="evidence changed"):
        replay.prepare_v7(context)
    artifact.write_bytes(files[replay.V7_BUNDLE + "/atlas.json"])
    (context / "docs/v7-owner.md").write_bytes(b"changed historical owner")
    with pytest.raises(ValueError, match="context changed"):
        replay.prepare_v7(context)


def test_v7_changed_receipt_fails_without_repair(v7_archive):
    _, context, _, _ = v7_archive
    replay.prepare_v7(context)
    receipt = context / "replay-context.json"
    receipt.write_bytes(b"altered receipt")
    with pytest.raises(ValueError, match="receipt changed"):
        replay.prepare_v7(context)
    assert receipt.read_bytes() == b"altered receipt"


@pytest.mark.parametrize("system,matplotlib", [("Linux", "3.11.1"), ("Windows", "3.11.2")])
def test_v7_replay_requires_its_own_profile_without_fallback(v7_archive, monkeypatch, system, matplotlib):
    _, context, _, _ = v7_archive
    monkeypatch.setattr(replay.platform, "system", lambda: system)
    monkeypatch.setattr(replay.platform, "python_version", lambda: "3.14.8")
    monkeypatch.setattr(replay, "version", lambda name: matplotlib)
    monkeypatch.setattr(replay.subprocess, "run", lambda *a, **k: pytest.fail("checker must not start"))
    with pytest.raises(ValueError, match="Use prepare-v7.*no fallback"):
        replay.run_v7(context)
    assert (context / replay.V7_BUNDLE / "atlas.json").is_file()


@pytest.mark.parametrize("returncode", [0, 1])
def test_v7_runs_original_checker_with_exact_comparison_and_isolated_imports(v7_archive, monkeypatch, returncode):
    root, context, _, _ = v7_archive
    monkeypatch.setattr(replay.platform, "system", lambda: "Windows")
    monkeypatch.setattr(replay.platform, "python_version", lambda: "3.14.8")
    monkeypatch.setattr(replay, "version", lambda name: "3.11.1")
    monkeypatch.setenv("PYTHONPATH", "live/src" + replay.os.pathsep + "deps/site-packages")
    calls = []
    def run(command, **options):
        calls.append((command, options))
        return SimpleNamespace(returncode=returncode, stdout='{"checked":"v7","outputs":8}',
                               stderr="atlas replay mismatch: preparation-domain.png" if returncode else "")
    monkeypatch.setattr(replay.subprocess, "run", run)
    if returncode:
        with pytest.raises(ValueError, match="atlas replay mismatch: preparation-domain.png"):
            replay.run_v7(context)
    else:
        result = replay.run_v7(context)
        assert result["scientific_replay"]
        assert result["checks"][0]["result"]["outputs"] == 8
    assert len(calls) == 1
    command, options = calls[0]
    assert command == [replay.sys.executable, "scripts/derive_nsc_cosmological_atlas_v7.py", "--check"]
    assert options["cwd"] == context
    assert "live/src" not in options["env"]["PYTHONPATH"]
    assert str(context / "src") in options["env"]["PYTHONPATH"]
    assert options["env"]["PYTHONNOUSERSITE"] == "1"
    assert options["env"]["GIT_WORK_TREE"] == str(context)
    assert options["env"]["GIT_ALTERNATE_OBJECT_DIRECTORIES"] == str(root / "archive/source-history/objects")


def test_legacy_all_keeps_original_six_atlases_and_profile(tmp_path, monkeypatch):
    monkeypatch.setattr(replay, "prepare", lambda context: {"source_commit": replay.BASELINE})
    monkeypatch.setattr(replay, "prepare_v7", lambda *a: pytest.fail("v7 must be explicitly selected"))
    calls = []
    def run(command, **options):
        calls.append((command, options))
        return SimpleNamespace(returncode=0, stdout="{}", stderr="")
    monkeypatch.setattr(replay.subprocess, "run", run)
    result = replay.run("all", tmp_path)
    assert result["source_commit"] == replay.BASELINE
    assert len(calls) == 9
    assert all("v7" not in str(command) for command, _ in calls)
    assert all(options["cwd"] == tmp_path / "lab" for _, options in calls[:6])
