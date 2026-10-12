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


@pytest.fixture(params=(7, 8))
def root_atlas_archive(tmp_path, monkeypatch, request):
    """Small original-commit archive with an owner allowed to evolve live."""
    label = "v" + str(request.param)
    commit = getattr(replay, label.upper() + "_COMMIT")
    bundle = getattr(replay, label.upper() + "_BUNDLE")
    root = tmp_path / "current"
    context = tmp_path / ("replay-" + label)
    files = {"src/recursive_horizons/__init__.py": b"original package",
             "scripts/derive_nsc_cosmological_atlas_" + label + ".py": b"original checker",
             "docs/atlas-owner.md": b"original owner",
             "pyproject.toml": b"original dependency profile",
             bundle + "/atlas.json": b'{"original":true}\n'}
    manifest = {"source_bindings": {name: replay.sha(data) for name, data in files.items()
                                    if name.startswith(("src/", "scripts/", "docs/"))},
                "outputs": {"atlas.json": {"sha256": replay.sha(files[bundle + "/atlas.json"])}},
                "runtime": {"python": "3.14.8", "matplotlib": "3.11.1"}}
    files[bundle + "/manifest.json"] = json.dumps(manifest).encode()
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
            assert args[1] == commit + ":" + bundle + "/manifest.json"
            return files[bundle + "/manifest.json"]
        assert args[:3] == ("archive", "--format=tar", commit)
        return stream.getvalue()
    monkeypatch.setattr(replay, "ROOT", root)
    monkeypatch.setattr(replay, "git", original_git)
    return SimpleNamespace(root=root, context=context, files=files, calls=calls,
                           bundle=bundle, label=label, commit=commit,
                           prepare=getattr(replay, "prepare_" + label),
                           run=getattr(replay, "run_" + label))


def test_root_atlas_preparation_uses_fixed_commit_despite_live_owner_changes(root_atlas_archive):
    row = root_atlas_archive
    (row.root / "docs/atlas-owner.md").write_bytes(b"current scientific successor")
    result = row.prepare(row.context)
    assert result["source_commit"] == row.commit
    assert not result["scientific_replay"]
    assert result["evidence_authenticated"] == 2
    assert (row.context / "docs/atlas-owner.md").read_bytes() == row.files["docs/atlas-owner.md"]
    assert row.prepare(row.context) == result
    receipt = json.loads((row.context / "replay-context.json").read_bytes())
    assert receipt["files"]["docs/atlas-owner.md"] == replay.sha(row.files["docs/atlas-owner.md"])
    assert all(call[0] in {"show", "archive"} for call in row.calls)


def test_root_atlas_retained_evidence_and_frozen_context_cannot_be_replaced(root_atlas_archive):
    row = root_atlas_archive
    row.prepare(row.context)
    artifact = row.root / row.bundle / "atlas.json"
    artifact.write_bytes(b"changed evidence")
    with pytest.raises(ValueError, match="evidence changed"):
        row.prepare(row.context)
    artifact.write_bytes(row.files[row.bundle + "/atlas.json"])
    (row.context / "docs/atlas-owner.md").write_bytes(b"changed historical owner")
    with pytest.raises(ValueError, match="context changed"):
        row.prepare(row.context)


def test_root_atlas_changed_receipt_fails_without_repair(root_atlas_archive):
    row = root_atlas_archive
    row.prepare(row.context)
    receipt = row.context / "replay-context.json"
    receipt.write_bytes(b"altered receipt")
    with pytest.raises(ValueError, match="receipt changed"):
        row.prepare(row.context)
    assert receipt.read_bytes() == b"altered receipt"


def test_root_atlas_original_manifest_source_binding_is_checked(root_atlas_archive):
    row = root_atlas_archive
    path = row.bundle + "/manifest.json"
    manifest = json.loads(row.files[path])
    manifest["source_bindings"]["docs/atlas-owner.md"] = "0" * 64
    row.files[path] = json.dumps(manifest).encode()
    with pytest.raises(ValueError, match="original source does not match its manifest"):
        row.prepare(row.context)
    assert not (row.context / "replay-context.json").exists()


@pytest.mark.parametrize("system,matplotlib", [("Linux", "3.11.1"), ("Windows", "3.11.2")])
def test_root_atlas_replay_requires_its_own_profile_without_fallback(root_atlas_archive, monkeypatch, system, matplotlib):
    row = root_atlas_archive
    monkeypatch.setattr(replay.platform, "system", lambda: system)
    monkeypatch.setattr(replay.platform, "python_version", lambda: "3.14.8")
    monkeypatch.setattr(replay, "version", lambda name: matplotlib)
    monkeypatch.setattr(replay.subprocess, "run", lambda *a, **k: pytest.fail("checker must not start"))
    with pytest.raises(ValueError, match="Use prepare-" + row.label + ".*no fallback"):
        row.run(row.context)
    assert (row.context / row.bundle / "atlas.json").is_file()


@pytest.mark.parametrize("returncode", [0, 1])
def test_root_atlas_runs_original_checker_with_exact_comparison_and_isolated_imports(root_atlas_archive, monkeypatch, returncode):
    row = root_atlas_archive
    monkeypatch.setattr(replay.platform, "system", lambda: "Windows")
    monkeypatch.setattr(replay.platform, "python_version", lambda: "3.14.8")
    monkeypatch.setattr(replay, "version", lambda name: "3.11.1")
    monkeypatch.setenv("PYTHONPATH", "live/src" + replay.os.pathsep + "deps/site-packages")
    calls = []
    def run(command, **options):
        calls.append((command, options))
        return SimpleNamespace(returncode=returncode, stdout='{"checked":"original","outputs":8}',
                               stderr="atlas replay mismatch: preparation-domain.png" if returncode else "")
    monkeypatch.setattr(replay.subprocess, "run", run)
    if returncode:
        with pytest.raises(ValueError, match="atlas replay mismatch: preparation-domain.png"):
            row.run(row.context)
    else:
        result = row.run(row.context)
        assert result["scientific_replay"]
        assert result["checks"][0]["result"]["outputs"] == 8
    assert len(calls) == 1
    command, options = calls[0]
    assert command == [replay.sys.executable, "scripts/derive_nsc_cosmological_atlas_" + row.label + ".py", "--check"]
    assert options["cwd"] == row.context
    assert "live/src" not in options["env"]["PYTHONPATH"]
    assert str(row.context / "src") in options["env"]["PYTHONPATH"]
    assert options["env"]["PYTHONNOUSERSITE"] == "1"
    assert options["env"]["GIT_WORK_TREE"] == str(row.context)
    assert options["env"]["GIT_ALTERNATE_OBJECT_DIRECTORIES"] == str(row.root / "archive/source-history/objects")


def test_legacy_all_keeps_original_six_atlases_and_profile(tmp_path, monkeypatch):
    monkeypatch.setattr(replay, "prepare", lambda context: {"source_commit": replay.BASELINE})
    monkeypatch.setattr(replay, "prepare_v7", lambda *a: pytest.fail("v7 must be explicitly selected"))
    monkeypatch.setattr(replay, "prepare_v8", lambda *a: pytest.fail("v8 must be explicitly selected"))
    calls = []
    def run(command, **options):
        calls.append((command, options))
        return SimpleNamespace(returncode=0, stdout="{}", stderr="")
    monkeypatch.setattr(replay.subprocess, "run", run)
    result = replay.run("all", tmp_path)
    assert result["source_commit"] == replay.BASELINE
    assert len(calls) == 9
    assert all("v7" not in str(command) and "v8" not in str(command) for command, _ in calls)
    assert all(options["cwd"] == tmp_path / "lab" for _, options in calls[:6])
