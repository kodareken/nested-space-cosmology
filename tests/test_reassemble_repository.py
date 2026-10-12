"""Recovery controls for the root consolidation, independent of its inventory."""
import hashlib
import importlib.util
from pathlib import Path
import zipfile
from types import SimpleNamespace

import pytest


spec = importlib.util.spec_from_file_location(
    "reassemble_repository", Path(__file__).resolve().parents[1] / "scripts/reassemble_repository.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def record(action="move"):
    return {"rows": [{"source": "old/source.bin", "destination": "new/source.bin",
                       "sha256": hashlib.sha256(b"research").hexdigest(),
                       "action": action, "origin": "preview", "original_path": "source.bin"}]}


def test_paths_cannot_escape_workspace(tmp_path):
    with pytest.raises(ValueError, match="leaves"):
        module.confined(tmp_path, "../outside")
    with pytest.raises(ValueError, match="leaves"):
        module.confined(tmp_path, ".")


def test_apply_moves_only_verified_bytes_and_is_resumable(tmp_path):
    source = tmp_path / "old/source.bin"
    source.parent.mkdir()
    source.write_bytes(b"research")
    assert module.apply(record(), tmp_path) == {"move": 1}
    assert not source.exists()
    assert (tmp_path / "new/source.bin").read_bytes() == b"research"
    assert module.apply(record(), tmp_path) == {"already_preserved": 1}
    assert module.check(record(), tmp_path)["total"] == 1


def test_destination_collision_does_not_overwrite_or_move(tmp_path):
    source = tmp_path / "old/source.bin"
    destination = tmp_path / "new/source.bin"
    source.parent.mkdir()
    destination.parent.mkdir()
    source.write_bytes(b"research")
    destination.write_bytes(b"different research")
    with pytest.raises(ValueError, match="destination drift"):
        module.apply(record(), tmp_path)
    assert source.read_bytes() == b"research"
    assert destination.read_bytes() == b"different research"


def test_source_drift_does_not_move(tmp_path):
    source = tmp_path / "old/source.bin"
    source.parent.mkdir()
    source.write_bytes(b"changed")
    with pytest.raises(ValueError, match="source drift"):
        module.apply(record(), tmp_path)
    assert source.exists()
    assert not (tmp_path / "new/source.bin").exists()


def test_copy_preserves_original_and_missing_archive_is_rejected(tmp_path):
    source = tmp_path / "old/source.bin"
    source.parent.mkdir()
    source.write_bytes(b"research")
    module.apply(record("copy"), tmp_path)
    assert source.read_bytes() == b"research"
    (tmp_path / "new/source.bin").unlink()
    with pytest.raises(ValueError, match="original bytes unavailable"):
        module.check(record(), tmp_path)


def test_current_corruption_never_falls_back_to_original(tmp_path, monkeypatch):
    destination = tmp_path / "new/source.bin"
    destination.parent.mkdir()
    destination.write_bytes(b"corrupt")
    monkeypatch.setattr(module, "original_bytes_available", lambda *args: True)
    with pytest.raises(ValueError, match="original bytes unavailable"):
        module.check(record(), tmp_path)


def test_explicit_successor_keeps_original_hash_and_checks_both(tmp_path, monkeypatch):
    destination = tmp_path / "new/source.bin"
    destination.parent.mkdir()
    destination.write_bytes(b"deliberate successor")
    manifest = record()
    monkeypatch.setattr(module, "original_bytes_available", lambda *args: True)
    module.record_successors(manifest, ["new/source.bin"], tmp_path)
    assert manifest["rows"][0]["sha256"] == hashlib.sha256(b"research").hexdigest()
    assert module.check(manifest, tmp_path)["declared_originals_verified"] == 1
    monkeypatch.setattr(module, "original_bytes_available", lambda *args: False)
    with pytest.raises(ValueError, match="original:new/source.bin|original:source.bin"):
        module.check(manifest, tmp_path)


def packed_record(tmp_path, data=b"research", member="new/source.bin"):
    manifest = record()
    archive = tmp_path / "recovery.zip"
    with zipfile.ZipFile(archive, "w") as packed:
        packed.writestr(member, data)
    manifest["rows"][0].update(storage_archive="recovery.zip", archive_member="new/source.bin")
    manifest["storage_archives"] = [{"path": "recovery.zip", "sha256": module.digest(archive),
                                     "bytes": archive.stat().st_size, "members": 1}]
    return manifest


def test_packed_bytes_support_absent_local_ignored_copy(tmp_path):
    manifest = packed_record(tmp_path)
    assert module.check(manifest, tmp_path)["packed_only_destinations"] == 1
    (tmp_path / "recovery.zip").unlink()
    with pytest.raises(ValueError, match="storage archive mismatch"):
        module.check(manifest, tmp_path)


def test_missing_packed_member_fails_even_with_authenticated_zip(tmp_path):
    with pytest.raises(ValueError, match="member inventory mismatch"):
        module.check(packed_record(tmp_path, member="different.bin"), tmp_path)


def test_corrupt_packed_member_fails_even_with_authenticated_zip(tmp_path):
    with pytest.raises(ValueError, match="member mismatch"):
        module.check(packed_record(tmp_path, data=b"changed"), tmp_path)


def test_corrupt_current_copy_is_not_hidden_by_good_packed_original(tmp_path):
    manifest = packed_record(tmp_path)
    destination = tmp_path / "new/source.bin"
    destination.parent.mkdir()
    destination.write_bytes(b"changed")
    with pytest.raises(ValueError, match="original bytes unavailable"):
        module.check(manifest, tmp_path)


def test_runtime_archive_is_deterministic_across_directories(tmp_path, monkeypatch):
    monkeypatch.setattr(module.subprocess, "run", lambda *args, **kwargs:
                        SimpleNamespace(returncode=1, stdout=b"", stderr=b""))
    hashes = []
    for label in ("first", "second"):
        root = tmp_path / label
        source = root / "new/source.bin"
        source.parent.mkdir(parents=True)
        source.write_bytes(b"research")
        manifest = record()
        manifest["rows"][0]["disposition"] = "runtime_preserved"
        hashes.append(module.pack_local_storage(manifest, root)["sha256"])
        assert module.check(manifest, root)["packed_members_verified"] == 1
    assert hashes[0] == hashes[1]
