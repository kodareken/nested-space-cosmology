#!/usr/bin/env python3
"""Inspect or continue the authorized RCV3 TDG8 event-23 projection.

Status is read-only.  Run consumes only the authenticated persisted RCV3
projection, rebuilds the immutable original event plan, restores persisted
generation-one state, and advances no farther than event 23 -> 24.  It has no
generation-zero, bootstrap, candidate, or calibration route.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
import os
from pathlib import Path
import re
import stat
import sys
from typing import Any, Mapping, NoReturn


ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = str(ROOT / "src")
sys.path[:] = [entry for entry in sys.path if entry != _SOURCE_ROOT]
sys.path.insert(0, str(ROOT))
from scripts import run_fgc_pro19_event1 as event1  # noqa: E402

sys.path[:] = [entry for entry in sys.path if entry not in {_SOURCE_ROOT, str(ROOT)}]
sys.path[:0] = [_SOURCE_ROOT, str(ROOT)]

from recursive_horizons.fgc.evolution import hlt16_campaign_recovery as recovery  # noqa: E402
from recursive_horizons.fgc.evolution import hlt16_campaign_runtime as runtime  # noqa: E402
from recursive_horizons.fgc.evolution import tdg8_rcv3_execution_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
    HLT16CampaignStoreError,
)
from recursive_horizons.fgc.evolution.proto17_pure_construction import MEMBER_KEYS  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_progression_contract import (  # noqa: E402
    construct_first_event,
)
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    build_static_gr0_shells,
)


RUNNER_ID = "FGC-1-TDG8-RCV3-EVENT-RUNNER"
SCHEMA = "FGC-1-TDG8-RCV3-event-runner-v1"
EXPECTED_DESTINATION_PATH = Path("runs/fgc-2-sf1/tdg8-rcv3/calibration")
EXPECTED_PROJECTION_ID = "FGC-1-TDG8-RCV3-PROJECTION-1"
ORIGINAL_EXECUTION_COMMIT = "6df67967e6e0bc63eca4dc6951a60d3787368f26"
ORIGINAL_PLAN_SHA256 = "e649e477b98411a3051f9946f2a9248b875d2b598818b8983ce06caa969d2126"
EXPECTED_PROTOCOL = "FGC-2-SF1-PROTO18"
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")
_MAX_DETAIL = 480


class TDG8RCV3EventRunnerError(RuntimeError):
    """A fail-closed error with a stable owning layer and bounded detail."""

    def __init__(self, owner: str, code: str, detail: str) -> None:
        self.owner = owner
        self.code = code
        self.detail = _sanitize(detail)
        super().__init__(f"{owner}/{code}: {self.detail}")

    def as_record(self) -> dict[str, object]:
        return {
            "schema": SCHEMA,
            "runner_id": RUNNER_ID,
            "state": "invalid_implementation_or_nonconverged_run",
            "error_owner": self.owner,
            "error_code": self.code,
            "error_detail": self.detail,
            "candidate_branch_opened": False,
            "calibration_result_earned": False,
            "physical_result_earned": False,
        }


def _sanitize(value: object) -> str:
    text = " ".join(str(value).replace(str(ROOT), "<repo>").split())
    text = text.encode("ascii", "backslashreplace").decode("ascii")
    return text[:_MAX_DETAIL] or "no further detail"


def _fail(owner: str, code: str, detail: object) -> NoReturn:
    raise TDG8RCV3EventRunnerError(owner, code, str(detail))


def _from_exception(owner: str, code: str, exc: BaseException) -> NoReturn:
    _fail(owner, code, f"{type(exc).__name__}: {exc}")


def _fixed_root(root: Path) -> Path:
    supplied = Path(os.path.abspath(os.fspath(root)))
    try:
        metadata = supplied.lstat()
    except OSError as exc:
        _from_exception("provenance", "repository_absent", exc)
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _fail("provenance", "repository_unsafe", "repository root is not a real directory")
    if supplied.resolve() != ROOT.resolve():
        _fail("authority", "repository_scope", "runner requires its fixed repository")
    return ROOT.resolve()


def _read_authority_leaf(
    root: Path,
    relative: str,
    *,
    label: str,
    maximum_bytes: int,
) -> bytes:
    """Read one fixed repository leaf without following a redirected leaf."""

    candidate = Path(relative)
    if (
        candidate.is_absolute()
        or not candidate.parts
        or any(part in {"", ".", ".."} for part in candidate.parts)
    ):
        raise ValueError(f"{label} path is unsafe")
    path = root.joinpath(candidate)
    if path.resolve(strict=False) != path:
        raise ValueError(f"{label} path is redirected")
    try:
        before = path.lstat()
    except OSError as exc:
        raise OSError(f"{label} is unavailable") from exc
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_size > maximum_bytes
    ):
        raise ValueError(f"{label} is not a bounded regular file")

    flags = os.O_RDONLY
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    try:
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or (opened.st_dev, opened.st_ino) != (before.st_dev, before.st_ino)
            or opened.st_size != before.st_size
            or opened.st_size > maximum_bytes
        ):
            raise ValueError(f"{label} changed before read")
        chunks: list[bytes] = []
        remaining = maximum_bytes + 1
        while remaining:
            chunk = os.read(descriptor, min(remaining, 1024 * 1024))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        raw = b"".join(chunks)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if len(raw) > maximum_bytes:
        raise ValueError(f"{label} exceeds its byte bound")
    if (
        (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
        != (opened.st_dev, opened.st_ino, opened.st_size, opened.st_mtime_ns)
        or len(raw) != opened.st_size
    ):
        raise ValueError(f"{label} changed during read")
    return raw


def _destination(root: Path) -> Path:
    if (
        Path(authority.DESTINATION_PATH) != EXPECTED_DESTINATION_PATH
        or authority.PROJECTION_ID != EXPECTED_PROJECTION_ID
    ):
        _fail("authority", "destination_scope", "RCV3 projected store identity differs")
    path = root / EXPECTED_DESTINATION_PATH
    if path.resolve(strict=False) != path:
        _fail("provenance", "destination_redirected", "RCV3 campaign destination is redirected")
    try:
        metadata = path.lstat()
    except OSError as exc:
        _from_exception("provenance", "destination_absent", exc)
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _fail("provenance", "destination_unsafe", "RCV3 campaign destination is unsafe")
    return path


def _execution_authority(
    root: Path,
    commit: str,
) -> authority.TDG8RCV3ExecutionAuthority:
    if not isinstance(commit, str) or not _COMMIT.fullmatch(commit):
        _fail("authority", "commit_malformed", "authority commit is not a full Git SHA")
    try:
        config_raw = _read_authority_leaf(
            root,
            authority.CONFIG_PATH,
            label="RCV3 execution config",
            maximum_bytes=1024 * 1024,
        )
        result_raw = _read_authority_leaf(
            root,
            authority.RESULT_PATH,
            label="RCV3 execution result",
            maximum_bytes=4 * 1024 * 1024,
        )
        receipt = authority.authorize_execution(root, config_raw, result_raw, commit)
    except TDG8RCV3EventRunnerError:
        raise
    except Exception as exc:
        _from_exception("authority", "committed_image_rejected", exc)
    if not isinstance(receipt, authority.TDG8RCV3ExecutionAuthority):
        _fail("authority", "typed_scope_rejected", "RCV3 execution authority type differs")
    try:
        observed = (
            receipt.authority_commit,
            receipt.pref1_authority_commit,
            receipt.pref1_result_sha256,
            receipt.projection_id,
            receipt.receipt_sha256,
            receipt.original_execution_commit,
            receipt.original_plan_sha256,
            receipt.campaign_id,
            receipt.destination_path,
            receipt.entry_checkpoint_sha256,
            receipt.entry_journal_sha256,
            receipt.temporal_retry_cap,
            receipt.minimum_macro_step_hex,
            receipt.candidate_branches_forbidden,
        )
    except AttributeError as exc:
        _from_exception("authority", "typed_scope_rejected", exc)
    expected = (
        commit,
        authority.PREF1_AUTHORITY_COMMIT,
        authority.PREF1_RESULT_SHA256,
        EXPECTED_PROJECTION_ID,
        authority.RECEIPT_SHA256,
        ORIGINAL_EXECUTION_COMMIT,
        ORIGINAL_PLAN_SHA256,
        authority.CAMPAIGN_ID,
        EXPECTED_DESTINATION_PATH.as_posix(),
        authority.ENTRY_CHECKPOINT_SHA256,
        authority.ENTRY_JOURNAL_SHA256,
        authority.TEMPORAL_RETRY_CAP,
        authority.MINIMUM_MACRO_STEP_HEX,
        True,
    )
    if (
        authority.ORIGINAL_EXECUTION_COMMIT != ORIGINAL_EXECUTION_COMMIT
        or authority.ORIGINAL_PLAN_SHA256 != ORIGINAL_PLAN_SHA256
        or observed != expected
    ):
        _fail("authority", "typed_scope_rejected", "typed RCV3 execution authority differs")
    return receipt


def _original_plan(root: Path) -> Any:
    try:
        plan = replace(
            construct_first_event(
                root,
                authorization_commit=ORIGINAL_EXECUTION_COMMIT,
            ),
            campaign_id=authority.CAMPAIGN_ID,
        )
    except Exception as exc:
        _from_exception("authority", "original_plan_unavailable", exc)
    if (
        plan.sha256 != ORIGINAL_PLAN_SHA256
        or plan.authorization_commit != ORIGINAL_EXECUTION_COMMIT
        or plan.protocol_artifact_id != EXPECTED_PROTOCOL
        or plan.campaign_id != authority.CAMPAIGN_ID
        or plan.branch != "GR-0"
        or plan.amplitude != "3"
        or plan.common_event_index != authority.EVENT
        or plan.start_time.rational != "23/16"
        or plan.target_time.rational != "3/2"
        or tuple(member.member_key for member in plan.members) != MEMBER_KEYS
    ):
        _fail("authority", "original_plan_rejected", "original immutable plan differs")
    return plan


def _inspect_boundary(
    root: Path,
    receipt: authority.TDG8RCV3ExecutionAuthority,
) -> authority.RCV3BoundedBoundary:
    try:
        boundary = authority.inspect_descendant(root, receipt)
    except Exception as exc:
        _from_exception("provenance", "descendant_rejected", exc)
    if not isinstance(boundary, authority.RCV3BoundedBoundary):
        _fail("authority", "boundary_type", "RCV3 descendant boundary type differs")
    try:
        event = boundary.event
        target = dict(boundary.target)
        disposition = boundary.disposition
    except (AttributeError, TypeError, ValueError) as exc:
        _from_exception("authority", "boundary_shape", exc)
    if event == authority.EVENT:
        if target != dict(authority.TARGET):
            _fail("scope", "boundary_event", "event-23 boundary target differs")
    elif event == authority.SUCCESSOR_EVENT:
        if target != dict(authority.SUCCESSOR_TARGET) or disposition != "event_complete":
            _fail("scope", "boundary_event", "event-24 boundary is not exact completion")
    else:
        _fail("scope", "boundary_event", "RCV3 boundary crossed its single-event scope")
    return boundary


def _summary(state: str, **extra: object) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "runner_id": RUNNER_ID,
        "state": state,
        "projection_id": EXPECTED_PROJECTION_ID,
        "candidate_branch_opened": False,
        "calibration_result_earned": False,
        "physical_result_earned": False,
        **extra,
    }


def status(root: Path, *, authority_commit: str) -> dict[str, Any]:
    repository = _fixed_root(root)
    receipt = _execution_authority(repository, authority_commit)
    _original_plan(repository)
    boundary = _inspect_boundary(repository, receipt)
    state = (
        "same_event_terminal"
        if boundary.terminal
        else "event_one_complete"
        if boundary.event == authority.SUCCESSOR_EVENT
        else "temporal_retry_pending"
        if boundary.temporal_pending_members
        else "publication_recovery_pending"
        if boundary.suffix_classification != "clean"
        else "same_event_ready"
    )
    return _summary(
        state,
        authority_commit=receipt.authority_commit,
        original_execution_commit=receipt.original_execution_commit,
        original_plan_sha256=receipt.original_plan_sha256,
        boundary=asdict(boundary),
    )


def _validate_run_result(result: object) -> Mapping[str, Any]:
    if not isinstance(result, Mapping):
        _fail("runtime", "result_shape", "same-event runner returned no result mapping")
    if (
        result.get("candidate_branch_opened") is not False
        or result.get("calibration_result_earned") is not False
        or result.get("physical_result_earned") is not False
        or result.get("authorization_commit") != ORIGINAL_EXECUTION_COMMIT
        or result.get("campaign_id") != authority.CAMPAIGN_ID
        or result.get("plan_sha256") != ORIGINAL_PLAN_SHA256
        or result.get("event") not in {authority.EVENT, authority.SUCCESSOR_EVENT}
    ):
        _fail("scope", "result_promotion", "same-event result crossed RCV3 scope")
    if result.get("event") == authority.SUCCESSOR_EVENT:
        if (
            result.get("state") != "event_one_complete"
            or result.get("target") != authority.SUCCESSOR_TARGET
        ):
            _fail("scope", "event_boundary", "event-24 result is not exact completion")
    elif (
        result.get("state") not in {"scientific_terminal", "invalid_terminal"}
        or result.get("target") != authority.TARGET
    ):
        _fail("scope", "event_boundary", "event-23 result is not a typed terminal")
    return result


def _require_runtime_checkpoint_scope(checkpoint: object) -> None:
    """Close the inspection/recovery race before the writer can run.

    The store's writer lease is itself exact-checkpoint bound.  This additional
    check ensures that a lawful descendant observed after the read-only AUTH1
    inspection is still inside event 23, or is already the exact event-24
    completion boundary, before `_advance_authenticated_event` may acquire it.
    """
    try:
        event = checkpoint.event
        target = dict(checkpoint.target)
        disposition = checkpoint.disposition
    except (AttributeError, TypeError, ValueError) as exc:
        _from_exception("scope", "runtime_checkpoint_shape", exc)
    if event == authority.EVENT:
        if target != dict(authority.TARGET):
            _fail("scope", "runtime_checkpoint_event", "recovered event-23 target differs")
        return
    if event == authority.SUCCESSOR_EVENT:
        if target != dict(authority.SUCCESSOR_TARGET) or disposition != "event_complete":
            _fail(
                "scope", "runtime_checkpoint_event",
                "recovered event-24 checkpoint is not exact completion",
            )
        return
    _fail("scope", "runtime_checkpoint_event", "recovered checkpoint crossed one-event authority")


def run(root: Path, *, authority_commit: str) -> dict[str, Any]:
    repository = _fixed_root(root)
    receipt = _execution_authority(repository, authority_commit)
    plan = _original_plan(repository)
    boundary = _inspect_boundary(repository, receipt)
    if boundary.terminal and boundary.terminal_lock_present:
        return _summary(
            "same_event_terminal",
            authority_commit=receipt.authority_commit,
            original_execution_commit=receipt.original_execution_commit,
            original_plan_sha256=receipt.original_plan_sha256,
            boundary=asdict(boundary),
            output_created=False,
        )
    if boundary.event == authority.SUCCESSOR_EVENT:
        return _summary(
            "event_one_complete",
            authority_commit=receipt.authority_commit,
            original_execution_commit=receipt.original_execution_commit,
            original_plan_sha256=receipt.original_plan_sha256,
            boundary=asdict(boundary),
            output_created=False,
        )
    if boundary.safe_to_restart is not True:
        _fail(
            "writer",
            "restart_not_authorized",
            (
                "descendant is not safely recoverable: "
                f"writer={boundary.writer_state!r}, "
                f"suffix={boundary.suffix_classification!r}"
            ),
        )

    destination = _destination(repository)
    try:
        store = HLT16CampaignStore(destination)
    except Exception as exc:
        _from_exception("provenance", "campaign_store_rejected", exc)
    try:
        checkpoint = runtime.recover_persisted_generation_one(plan, store)
    except Exception as exc:
        _from_exception("provenance", "persisted_restore_rejected", exc)
    _require_runtime_checkpoint_scope(checkpoint)
    try:
        templates = build_static_gr0_shells(repository)
    except Exception as exc:
        _from_exception("runtime", "gr0_shell_factory_rejected", exc)
    if tuple(templates) != MEMBER_KEYS:
        _fail("scope", "gr0_shell_order", "static GR-0 shell order differs")
    try:
        result = event1._advance_authenticated_event(plan, store, checkpoint, templates)
    except recovery.HLT16RecoveryError as exc:
        _from_exception("recovery", "metadata_reconciliation_rejected", exc)
    except HLT16CampaignStoreError as exc:
        _from_exception("provenance", "campaign_store_rejected", exc)
    except event1.Proto19Event1RunnerError as exc:
        text = str(exc).lower()
        owner = "writer" if "writer" in text or "lease" in text else "runtime"
        _from_exception(owner, "same_event_resume_rejected", exc)
    except Exception as exc:
        _from_exception("runtime", "unexpected_same_event_failure", exc)
    checked = _validate_run_result(result)
    return _summary(
        str(checked["state"]),
        authority_commit=receipt.authority_commit,
        original_execution_commit=receipt.original_execution_commit,
        original_plan_sha256=receipt.original_plan_sha256,
        result=dict(checked),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority-commit", required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--status", action="store_true")
    mode.add_argument("--run", action="store_true")
    args = parser.parse_args(argv)
    try:
        result = (
            status(ROOT, authority_commit=args.authority_commit)
            if args.status
            else run(ROOT, authority_commit=args.authority_commit)
        )
    except TDG8RCV3EventRunnerError as exc:
        print(
            json.dumps(exc.as_record(), sort_keys=True, separators=(",", ":")),
            file=sys.stderr,
        )
        return 2
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
