from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
from io import BytesIO
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import tomllib
import unittest
from unittest import mock
import zipfile

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto18_production_inputs import (
    Proto18ProductionArchiveError,
    Proto18ProductionPathError,
    Proto18ProductionSelectorError,
    load_source_specs,
    read_source_archive,
    verify_production_corpus,
)
from recursive_horizons.fgc.evolution import proto18_production_inputs as inputs


def config() -> dict[str, object]:
    with (ROOT / "configs/fgc/fgc-1-pro18-frz1.toml").open("rb") as handle:
        return tomllib.load(handle)


class Proto18ProductionInputsTests(unittest.TestCase):
    def test_actual_two_container_corpus_recomputes_all_six_restart_identities(self) -> None:
        corpus = verify_production_corpus(ROOT, config())
        self.assertEqual(tuple(corpus.archives), ("PROTO12", "RSP2"))
        self.assertEqual(
            tuple(corpus.members),
            ("RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"),
        )
        self.assertEqual(
            corpus.archives["PROTO12"].archive_snapshot.sha256,
            "c784d4706911029883d14e1763c2d36c5dbe0e619a4e50416cc6d666d5b7ea17",
        )
        self.assertEqual(
            corpus.archives["RSP2"].archive_snapshot.sha256,
            "000546dcf3726c7882e80b7e70b64a5d2e060e2e0f7a567e174f714432ae7f34",
        )
        for source, archive in corpus.archives.items():
            self.assertEqual(archive.event_snapshot.sha256, archive.embedded_event_log_sha256)
            # This is the single-byte-image property: no parser receives a
            # fresh file read after the retained, hashed source snapshot.
            self.assertEqual(sha256(archive.archive_snapshot.payload).hexdigest(), archive.archive_snapshot.sha256)
            self.assertEqual(sha256(archive.event_snapshot.payload).hexdigest(), archive.event_snapshot.sha256)
            self.assertTrue(archive.metadata["terminal"], source)
        for key, member in corpus.members.items():
            self.assertEqual(member.physical_state_sha256, member.selector.state_sha256, key)
            self.assertEqual(member.restart_payload_sha256, member.selector.restart_sha256, key)
            self.assertEqual(set(member.array_sha256), {
                "u", "p", "q", "tracer_positions", "tracer_proper_times",
                "event_proper_times", "event_fields",
            })
            self.assertFalse(member.arrays["u"].flags.writeable, key)
            with self.assertRaises(ValueError):
                member.arrays["u"][0, 0] = 0.0
            with self.assertRaises(ValueError):
                member.arrays["u"].flags.writeable = True
            # Source metadata is part of the admitted evidence; consumers
            # cannot be allowed to rewrite it after hash classification.
            with self.assertRaises(TypeError):
                member.metadata["time"] = 0.0

    def test_selector_and_corpus_substitution_fail_before_raw_admission(self) -> None:
        value = config()
        bad = deepcopy(value)
        bad["selectors"]["member"][0]["source"] = "PROTO14"
        with self.assertRaises(Proto18ProductionSelectorError):
            load_source_specs(bad)
        bad = deepcopy(value)
        bad["selectors"]["member"][0]["array_prefix"] = bad["selectors"]["member"][1]["array_prefix"]
        with self.assertRaises(Proto18ProductionSelectorError):
            load_source_specs(bad)
        bad = deepcopy(value)
        bad["source_map"]["shared_containers"] = ["runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz"]
        with self.assertRaises(Proto18ProductionSelectorError):
            load_source_specs(bad)

    def test_unsafe_and_hostile_source_images_stop_before_numpy_load(self) -> None:
        value = config()
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "runs/fgc-2-sf1/proto12/calibration"
            target.mkdir(parents=True)
            (target / "events.jsonl").write_text("{}\n")
            # Retain the exact frozen path.  PRO18 must not gain an escape
            # hatch where a hostile image is tested under a substitute path.
            bad_archive = target / "latest-checkpoint.npz"
            with zipfile.ZipFile(bad_archive, "w") as archive:
                archive.writestr("unexpected.npy", b"not-an-npy")
            specs = load_source_specs(value)
            with self.assertRaises(Proto18ProductionArchiveError):
                read_source_archive(root, specs[0])
            payload = target / "payload.npz"
            payload.write_bytes(bad_archive.read_bytes())
            bad_archive.unlink()
            bad_archive.symlink_to(payload)
            with self.assertRaises(Proto18ProductionPathError):
                read_source_archive(root, specs[0])

    def test_same_inode_same_size_mutation_during_snapshot_is_not_admitted(self) -> None:
        value = config()
        specs = load_source_specs(value)
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "runs/fgc-2-sf1/proto12/calibration"
            target.mkdir(parents=True)
            source = target / "latest-checkpoint.npz"
            source.write_bytes(b"x" * 4096)
            (target / "events.jsonl").write_text("{}\n")
            original_read = inputs.os.read
            mutated = False

            def mutate_after_read(descriptor: int, count: int) -> bytes:
                nonlocal mutated
                payload = original_read(descriptor, count)
                if not mutated:
                    mutated = True
                    # Same descriptor, inode, and byte count: only the
                    # snapshot's time/identity guard can expose this change.
                    with source.open("r+b") as writer:
                        writer.write(b"y")
                        writer.flush()
                        os.fsync(writer.fileno())
                return payload

            with mock.patch.object(inputs.os, "read", side_effect=mutate_after_read):
                with self.assertRaises(Proto18ProductionPathError):
                    read_source_archive(root, specs[0])
            self.assertTrue(mutated)

    def test_symlinked_repository_root_is_rejected_before_traversal(self) -> None:
        with TemporaryDirectory() as temporary:
            base = Path(temporary)
            actual = base / "actual"
            target = actual / "runs/fgc-2-sf1/proto12/calibration"
            target.mkdir(parents=True)
            (target / "latest-checkpoint.npz").write_bytes(b"x")
            linked = base / "linked"
            linked.symlink_to(actual, target_is_directory=True)
            with self.assertRaises(Proto18ProductionPathError):
                inputs._regular_snapshot(
                    linked,
                    "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
                    label="root-symlink",
                )

    def test_metadata_overflow_final_hash_npy_version_and_zip_flags_fail_closed(self) -> None:
        with self.assertRaises(Proto18ProductionArchiveError):
            inputs._legacy_json(b'{"finite":1e999}', label="overflow")

        version_two = BytesIO()
        np.lib.format.write_array(version_two, np.zeros((2, 6), dtype="<f8"), version=(2, 0))
        with self.assertRaises(Proto18ProductionArchiveError):
            inputs._load_array(version_two.getvalue(), expected_shape=(2, 6), label="v2")

        value = config()
        bad_hash = deepcopy(value)
        bad_hash["selectors"]["member"][0]["state_sha256"] = "0" * 64
        bad_hash["selectors"]["member"][0]["physical_state_sha256"] = "0" * 64
        with self.assertRaises(Proto18ProductionSelectorError):
            verify_production_corpus(ROOT, bad_hash)

        spec = load_source_specs(value)[0]
        archive = read_source_archive(ROOT, spec)
        payload = bytearray(archive.archive_snapshot.payload)
        local = payload.find(b"PK\x03\x04")
        central = payload.find(b"PK\x01\x02")
        self.assertGreaterEqual(local, 0); self.assertGreaterEqual(central, 0)
        payload[local + 6] |= 1
        payload[central + 8] |= 1
        hostile = inputs.SourceSnapshot(
            archive.archive_snapshot.relative_path, bytes(payload),
            sha256(payload).hexdigest(), archive.archive_snapshot.device,
            archive.archive_snapshot.inode, archive.archive_snapshot.size,
            archive.archive_snapshot.mtime_ns, archive.archive_snapshot.ctime_ns,
        )
        with self.assertRaises(Proto18ProductionArchiveError):
            inputs._archive_members(hostile, spec)

    def test_adapter_never_opens_or_names_future_proto17_roots(self) -> None:
        observed: list[str] = []
        original_open = inputs.os.open

        def record_open(
            path: object, flags: int, *args: object, **kwargs: object,
        ) -> int:
            observed.append(os.fspath(path))
            return original_open(path, flags, *args, **kwargs)

        with mock.patch.object(inputs.os, "open", side_effect=record_open):
            corpus = verify_production_corpus(ROOT, config())
        self.assertEqual(len(corpus.members), 6)
        self.assertTrue(observed)
        self.assertFalse(any("/proto17/" in path or path.endswith("/proto17") for path in observed), observed)


if __name__ == "__main__":
    unittest.main()
