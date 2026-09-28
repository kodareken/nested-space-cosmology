"""Protect the migration boundary without running scientific campaigns."""
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import check_lab_snapshot as snapshot
import lab as launcher


def isolated_cache(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    shutil.copytree(ROOT / "lab/.source-history", tmp_path / "lab/.source-history")
    return tmp_path


def test_source_replay_without_old_repository(tmp_path):
    assert snapshot.check_objects(isolated_cache(tmp_path)) == 18


def test_corrupt_source_object_rejected(tmp_path):
    root = isolated_cache(tmp_path)
    path = next((root / "lab/.source-history/objects").glob("*/*"))
    path.write_bytes(b"corrupt")
    with pytest.raises(Exception):
        snapshot.check_objects(root)


def test_missing_source_object_rejected(tmp_path):
    root = isolated_cache(tmp_path)
    next((root / "lab/.source-history/objects").glob("*/*")).unlink()
    with pytest.raises(FileNotFoundError):
        snapshot.check_objects(root)


def test_environment_selects_lab_only():
    env = launcher.environment()
    assert env["PYTHONPATH"].split(__import__("os").pathsep) == [
        str(ROOT / "lab/src"), str(ROOT / "lab/scripts")]
    assert env["GIT_ALTERNATE_OBJECT_DIRECTORIES"] == str(ROOT / "lab/.source-history/objects")


def test_snapshot_rejects_modified_payload(tmp_path):
    (tmp_path / "lab").mkdir()
    (tmp_path / "docs").mkdir()
    (tmp_path / "lab/data.json").write_bytes(b"modified")
    manifest = {"schema": "NSC-CONSOLIDATED-LAB-v1", "files": [
        {"path": "lab/data.json", "bytes": 8, "sha256": hashlib.sha256(b"original").hexdigest()}],
        "intentionally_omitted": []}
    (tmp_path / "docs/lab-snapshot.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="snapshot mismatch"):
        snapshot.check_snapshot(tmp_path)


def test_snapshot_rejects_escaping_path(tmp_path):
    with pytest.raises(ValueError, match="invalid snapshot path"):
        snapshot.safe_path(tmp_path, "../outside")
