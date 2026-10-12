"""Real-store proof that SID3 preflight is byte-for-byte read-only."""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
SOURCE = str(ROOT / "src")
if SOURCE not in sys.path:
    sys.path.insert(0, SOURCE)

from scripts import run_fgc_pro19_event1 as runner  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
)
from recursive_horizons.fgc.evolution.proto19_progression_contract import (  # noqa: E402
    construct_first_event,
)
from recursive_horizons.fgc.evolution import proto19_resume_authority as sid3  # noqa: E402


LIVE_STORE = ROOT / sid3.STORE_PATH
EXACT_GENERATION8_CHECKPOINT = LIVE_STORE / (
    f"checkpoints/00000000000000000008-{sid3.GENERATION8_CHECKPOINT_SHA256}.json"
)
SID3_PROCESS_MARKER = "FGC-PRO19-SID3-EVENT1"


@dataclass(frozen=True, slots=True)
class _TreeEntry:
    kind: str
    mode: int
    device: int
    inode: int
    link_count: int
    uid: int
    gid: int
    size: int
    modified_ns: int
    changed_ns: int
    flags: int
    content_sha256: str | None


def _stat_identity(metadata: os.stat_result) -> tuple[int, ...]:
    return (
        metadata.st_mode,
        metadata.st_dev,
        metadata.st_ino,
        metadata.st_nlink,
        metadata.st_uid,
        metadata.st_gid,
        metadata.st_size,
        metadata.st_mtime_ns,
        metadata.st_ctime_ns,
        int(getattr(metadata, "st_flags", 0)),
    )


def _regular_file_sha256(path: Path, before: os.stat_result) -> str:
    descriptor = os.open(
        path,
        os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        active = os.fstat(descriptor)
        if _stat_identity(active) != _stat_identity(before):
            raise AssertionError(f"tree leaf raced before read: {path}")
        digest = sha256()
        total = 0
        while block := os.read(descriptor, 1024 * 1024):
            digest.update(block)
            total += len(block)
        after = os.fstat(descriptor)
        if _stat_identity(after) != _stat_identity(active) or total != after.st_size:
            raise AssertionError(f"tree leaf changed during read: {path}")
        return digest.hexdigest()
    finally:
        os.close(descriptor)


def _entry(path: Path, metadata: os.stat_result, *, kind: str) -> _TreeEntry:
    digest = _regular_file_sha256(path, metadata) if kind == "file" else None
    return _TreeEntry(
        kind=kind,
        mode=metadata.st_mode,
        device=metadata.st_dev,
        inode=metadata.st_ino,
        link_count=metadata.st_nlink,
        uid=metadata.st_uid,
        gid=metadata.st_gid,
        size=metadata.st_size,
        modified_ns=metadata.st_mtime_ns,
        changed_ns=metadata.st_ctime_ns,
        flags=int(getattr(metadata, "st_flags", 0)),
        content_sha256=digest,
    )


def _nofollow_tree(root: Path) -> dict[str, _TreeEntry]:
    """Capture every directory and regular leaf without following links."""

    root_metadata = root.stat(follow_symlinks=False)
    if stat.S_ISLNK(root_metadata.st_mode) or not stat.S_ISDIR(root_metadata.st_mode):
        raise AssertionError("copied campaign root is not a direct directory")
    observed = {".": _entry(root, root_metadata, kind="directory")}
    pending = [root]
    while pending:
        directory = pending.pop()
        with os.scandir(directory) as iterator:
            entries = sorted(iterator, key=lambda item: item.name)
        for child in entries:
            metadata = child.stat(follow_symlinks=False)
            path = Path(child.path)
            relative = path.relative_to(root).as_posix()
            if stat.S_ISLNK(metadata.st_mode):
                raise AssertionError(f"copied campaign contains a symlink: {relative}")
            if stat.S_ISDIR(metadata.st_mode):
                observed[relative] = _entry(path, metadata, kind="directory")
                pending.append(path)
            elif stat.S_ISREG(metadata.st_mode):
                observed[relative] = _entry(path, metadata, kind="file")
            else:
                raise AssertionError(
                    f"copied campaign contains a non-file object: {relative}"
                )
    return dict(sorted(observed.items()))


def _authority() -> SimpleNamespace:
    evidence = {
        path: (ROOT / path).read_bytes() for path in sid3.PROGRESSION_EVIDENCE_PATHS
    }
    plan = construct_first_event(
        Path("."),
        authorization_commit=sid3.ORIGINAL_AUTHORIZATION_COMMIT,
        evidence_bytes=evidence,
    )
    static_inputs = {
        path: (ROOT / path).read_bytes() for path in sid3.MANDATORY_STATIC_INPUT_PATHS
    }
    return SimpleNamespace(
        store_path="calibration",
        progression_plan=plan,
        static_input_bytes=static_inputs,
        implementation_commit="a" * 40,
        authority_commit="c" * 40,
        closure_sha256="d" * 64,
        config_sha256="e" * 64,
        result_sha256="f" * 64,
        checkpoint_generation=sid3.GENERATION,
        checkpoint_sha256=sid3.GENERATION8_CHECKPOINT_SHA256,
        journal_sequence=sid3.GENERATION8_JOURNAL_SEQUENCE,
        journal_tip_sha256=sid3.GENERATION8_JOURNAL_SHA256,
        original_authorization_commit=sid3.ORIGINAL_AUTHORIZATION_COMMIT,
        original_plan_sha256=sid3.ORIGINAL_PLAN_SHA256,
        campaign_id=sid3.CAMPAIGN_ID,
        recovery_member_key=sid3.GENERATION8_MEMBER_KEY,
        recovery_cursor_sha256=sid3.GENERATION8_CURSOR_SHA256,
        recovery_pending_cap_hex=sid3.GENERATION8_PENDING_CAP_HEX,
    )


def _marker_processes() -> tuple[str, ...]:
    completed = subprocess.run(
        ("/bin/ps", "-axo", "pid=,command="),
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return tuple(
        line for line in completed.stdout.splitlines() if SID3_PROCESS_MARKER in line
    )


class Proto19SID3RealStorePreflightTests(unittest.TestCase):
    def test_real_generation8_preflight_does_not_mutate_store(self) -> None:
        if not EXACT_GENERATION8_CHECKPOINT.is_file():
            self.skipTest(
                "live exact generation-eight SID3 fixture is absent; "
                f"expected {EXACT_GENERATION8_CHECKPOINT}"
            )

        live = HLT16CampaignStore(LIVE_STORE).authenticated_snapshot()
        if (
            live.checkpoint.generation != sid3.GENERATION
            or live.checkpoint.sha256 != sid3.GENERATION8_CHECKPOINT_SHA256
        ):
            self.skipTest(
                "live SID3 campaign no longer ends at the exact generation-eight fixture"
            )

        with TemporaryDirectory() as directory:
            temporary_root = Path(directory)
            copied_store = temporary_root / "calibration"
            shutil.copytree(LIVE_STORE, copied_store, symlinks=True)
            before = _nofollow_tree(copied_store)
            authority = _authority()
            self.assertEqual(_marker_processes(), ())

            with (
                patch.object(
                    runner,
                    "_validated_sid3_authority",
                    return_value=authority,
                ),
                patch.object(
                    runner,
                    "_acquire_writer",
                    wraps=runner._acquire_writer,
                ) as acquire_writer,
                patch.object(
                    subprocess,
                    "Popen",
                    wraps=subprocess.Popen,
                ) as process_spawn,
            ):
                result = runner.sid3_resume_preflight(
                    temporary_root,
                    authority=authority,
                    store_root=copied_store,
                )

            after = _nofollow_tree(copied_store)
            self.assertEqual(_marker_processes(), ())
            self.assertEqual(after, before)
            acquire_writer.assert_not_called()
            process_spawn.assert_not_called()

            self.assertEqual(result["checkpoint_generation"], sid3.GENERATION)
            self.assertEqual(
                result["checkpoint_sha256"], sid3.GENERATION8_CHECKPOINT_SHA256
            )
            self.assertEqual(
                result["journal_tip_sha256"], sid3.GENERATION8_JOURNAL_SHA256
            )
            self.assertTrue(result["six_member_restore_passed"])
            self.assertTrue(result["safe_to_resume_trajectory"])
            self.assertFalse(result["writer_lease_acquired"])
            self.assertFalse(result["output_created"])
            self.assertFalse(result["state_advanced"])
            self.assertFalse(result["candidate_branch_opened"])

            status = HLT16CampaignStore(copied_store).inspect_recovery()
            self.assertFalse(status.active_write)
            self.assertIsNone(status.writer_state)
            self.assertEqual(status.state, "clean_checkpoint")

            for prefix in ("checkpoints/", "journal/", "payloads/"):
                self.assertEqual(
                    {
                        key: value
                        for key, value in after.items()
                        if key.startswith(prefix)
                    },
                    {
                        key: value
                        for key, value in before.items()
                        if key.startswith(prefix)
                    },
                )


if __name__ == "__main__":
    unittest.main()
