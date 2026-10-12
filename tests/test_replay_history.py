import io
import tarfile
from pathlib import Path

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
