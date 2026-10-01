"""Protect the migration boundary without running scientific campaigns."""
import hashlib
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
import update_lab_snapshot as updater

# Consolidation-era cache. Later content carriers may follow these rows.
# Identity here is the Merkle object and the pinned file hash, not a count.
ORIGINAL_SOURCE_OBJECTS = (
    ("09f3faa364ca88b8ff482c413022ce38f40bef0d", "commit", 339, "f17ca2a8a430a11161f005be9dac2178f391565c884d3a4475abcee7db6f1f98"),
    ("6be0057f68cc1ca3fc021b0a2865ad87ea79f5c1", "tree", 847, "62fa898de8369496dd1ef1df07113bee809954000a74e71dbe7534176d16ec08"),
    ("22993e882ceb7032db8b233fb71231ee3b226316", "tree", 45, "08ee451400f015721b2ef91e55dcfc0f94eadab268d1618bc4fd6ea4bc93a6a7"),
    ("61286d365dbda4cf1e19b8f154d6f2a298bf92de", "tree", 9458, "2ef8afced2e0123b7a1245674fd78af1f58456751fcd544c3ea509e0539ee0a4"),
    ("a5a85acef3d98a8d6ff66f928a30725856a1ac6c", "blob", 7350, "d6f31c307705b49371242f04c4cd44ae326593b068904cd3d71506751b0a42d9"),
    ("750f08abd50ef461600baf9dd319d1ac38b50ca3", "commit", 317, "6b9db11f0a897997fd1532d39e5062584dafc601376b7f89eb3b591c2192ebb8"),
    ("157a55cdf6ceaf2f69350914046c54fe059855e1", "tree", 882, "8532e5f38ada1e2c1312d753781d953850512c5ba7661034c7fdde06b5c64f87"),
    ("f7b6d43cae388b9e87a586e1bc87ee62327ab96d", "blob", 1262, "935a7d3b693692d5e91d3710934ec68af4e3a532dfba6e905fc5985b89b356c1"),
    ("c7dc0f98cacb28afd030b089bdb9e1369288120c", "commit", 334, "b08c61f8d20dde6afb8eea725ba7ed4ad3d012ee509004e4b1593298316fbadb"),
    ("201a906a41d70ca033fb19f4bf7f1483f1fa2552", "tree", 902, "529fb8a47bc619d56cc636a3c3a1392020bc22c4973edf80bff6b8167636e9cc"),
    ("d778878f61150cfa0dd2cd3b4024aea276001d46", "tree", 45, "0d60fb8d81e087c2e4c67afdd49c86fcf9660dcad1d329ae65d4365e465d7abb"),
    ("4dc21108424b29558c9092c703184543721f5248", "tree", 11006, "db08f2b7c71150db61d7f4f08b3c9760e15fc05e484e5007355430525c0341ff"),
    ("31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227", "blob", 10229, "42cc1f09dd49d9c3ee0eeab22fe3b778f0562729e05414e82e0ef19658c82971"),
    ("e62d28f3230ae6ddb2712df71eb7c61a8ed65ad6", "commit", 341, "6b51fb912f469919df458c2a14b97c8a5e2ce69ff7f8929126cfae19a98d38ef"),
    ("d482613c0af362b6b7440788b335b136e24b8b4f", "tree", 847, "a123fa73f023a45a743a9a7d0712d28038e090a7118cccb0d746b948e0867e5a"),
    ("f2556bf248ef23192a4c864cb816090cb8641931", "tree", 45, "675dfff3c855ac4228ad398fbf28962e7b169d2fb770f0fb5d38b6ff1f38d23b"),
    ("826473ae9a2ff7a3da07718738ae48be7badcd68", "tree", 9798, "c360866cd2d48b03a62f39b9f376a619c24596cd2169862764988ceb2e813e99"),
    ("9efa09bd82b7003391890f1a19828de04148173f", "blob", 9693, "dfb09f21ba206c4b7272d98dc3270ae902f6ade26c2b203596723c0e1b218cfd"),
)
ORIGINAL_PINS = (
    ("09f3faa364ca88b8ff482c413022ce38f40bef0d", "src/recursive_horizons/nsc_local_incoming_family.py", "a5a85acef3d98a8d6ff66f928a30725856a1ac6c", "9d2b9fc57e2057d6c34d8943384428341b314370021bcd4156af88320c3de04b"),
    ("750f08abd50ef461600baf9dd319d1ac38b50ca3", "pyproject.toml", "f7b6d43cae388b9e87a586e1bc87ee62327ab96d", "b5bfe97f7225d19d6036e05d3199376885d4d626db05d49013495c94ec3e400e"),
    ("c7dc0f98cacb28afd030b089bdb9e1369288120c", "src/recursive_horizons/nsc_ks_energy_propagator.py", "31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227", "5f75023efc8c2b637fc4983b414cab50f15f90f8abc1a24dbc1d23ecef9ecf94"),
    ("e62d28f3230ae6ddb2712df71eb7c61a8ed65ad6", "src/recursive_horizons/nsc_ks_difference_envelope.py", "9efa09bd82b7003391890f1a19828de04148173f", "77c89c8e720ca97336ccdee441a645b5c73472cb6f15b09e8120f68ab53db196"),
)
RECOVERED_CONTENT_PINS = (
    ("594a11caeee760d172c4a4f73a7615b05cf30282", "src/recursive_horizons/nsc_spherical_null_expansion.py", "ac014452bed434cff03de0fcd52379c187098247", "1aa47f0b3d680621d23ac9329e0d9e1e38070b7df42fc2ec48848ea367fd82d5"),
    ("594a11caeee760d172c4a4f73a7615b05cf30282", "tests/test_nsc_local_boundary_independent.py", "2c586ad3fc31fe9e6a8fb4dba97fe8323bbec2f8", "68c2fa1559d705c3aae2d3b64d693af72311802e51cb5914432e1422063b4ea9"),
)


def isolated_cache(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    shutil.copytree(ROOT / "lab/.source-history", tmp_path / "lab/.source-history")
    return tmp_path


def test_source_replay_without_old_repository(tmp_path):
    root = isolated_cache(tmp_path)
    manifest = json.loads((root / "lab/.source-history/manifest.json").read_text())
    objects = manifest["objects"]
    pins = manifest["pins"]
    count = snapshot.check_objects(root)
    disk = {
        path.parent.name + path.name
        for path in (root / "lab/.source-history/objects").glob("*/*")
        if path.is_file()
    }
    assert count == len(objects)
    assert {row["oid"] for row in objects} == disk
    assert len(ORIGINAL_SOURCE_OBJECTS) == 18
    assert len(ORIGINAL_PINS) == 4
    assert [
        (row["oid"], row["type"], row["bytes"], row["sha256"]) for row in objects[:18]
    ] == list(ORIGINAL_SOURCE_OBJECTS)
    assert [
        (row["commit"], row["path"], row["blob"], row["sha256"]) for row in pins[:4]
    ] == list(ORIGINAL_PINS)
    pin_rows = [
        (row["commit"], row["path"], row["blob"], row["sha256"]) for row in pins
    ]
    historical_commits = {pin[0] for pin in ORIGINAL_PINS}
    for recovered in RECOVERED_CONTENT_PINS:
        assert recovered in pin_rows
        assert recovered[0] not in historical_commits
        assert recovered[2] in {row["oid"] for row in objects}
        assert recovered[0] in {row["oid"] for row in objects}


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


def test_successor_update_keeps_original_identity(tmp_path):
    (tmp_path / "docs").mkdir()
    (tmp_path / "lab").mkdir()
    (tmp_path / "lab/value.txt").write_text("new")
    target = tmp_path / "docs/lab-snapshot.json"
    original = {"path": "lab/value.txt", "sha256": "old", "bytes": 3,
                "original_sha256": "original", "original_git_blob": "blob"}
    target.write_text(json.dumps({"files": [original], "intentionally_omitted": []}))
    updater.update(tmp_path, ["lab/value.txt"], "parent")
    row = json.loads(target.read_text())["files"][0]
    assert row["original_sha256"] == "original"
    assert row["original_git_blob"] == "blob"
    assert row["sha256"] == hashlib.sha256(b"new").hexdigest()
    assert row["changed_from_commit"] == "parent"


def test_successor_cannot_silently_import_omitted_payload(tmp_path):
    (tmp_path / "docs").mkdir()
    target = tmp_path / "docs/lab-snapshot.json"
    target.write_text(json.dumps({"files": [], "intentionally_omitted": ["old.npz"]}))
    before = target.read_bytes()
    with pytest.raises(ValueError, match="current lab files"):
        updater.update(tmp_path, ["lab/old.npz"], "parent")
    assert target.read_bytes() == before
