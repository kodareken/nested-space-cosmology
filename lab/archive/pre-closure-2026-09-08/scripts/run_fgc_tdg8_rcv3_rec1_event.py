#!/usr/bin/env python3
"""Reconcile the exact RCV3 suffix, then continue only event 23 -> 24.

Status never acquires a lease.  Run is a two-phase wrapper: it closes only the
bound generation-nine seq11/seq12 suffix, proves the exact predicted
generation-ten checkpoint, retires the recovery lease, and only then invokes
the existing one-event engine.  It has no GEN0, bootstrap, candidate, or
promotion path.
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
from recursive_horizons.fgc.evolution import tdg8_rcv3_rec1_authority as authority  # noqa: E402
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


RUNNER_ID = "FGC-1-TDG8-RCV3-REC1-EVENT-RUNNER"
SCHEMA = "FGC-1-TDG8-RCV3-REC1-event-runner-v1"
EXPECTED_DESTINATION_PATH = Path("runs/fgc-2-sf1/tdg8-rcv3/calibration")
ORIGINAL_EXECUTION_COMMIT = authority.ORIGINAL_EXECUTION_COMMIT
ORIGINAL_PLAN_SHA256 = authority.ORIGINAL_PLAN_SHA256
EXPECTED_PROTOCOL = authority.TARGET_PROTOCOL
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")
_MAX_DETAIL = 480


class TDG8RCV3REC1RunnerError(RuntimeError):
    def __init__(self, owner: str, code: str, detail: str) -> None:
        self.owner, self.code, self.detail = owner, code, _sanitize(detail)
        super().__init__(f"{owner}/{code}: {self.detail}")

    def as_record(self) -> dict[str, object]:
        return {
            "schema": SCHEMA, "runner_id": RUNNER_ID,
            "state": "invalid_implementation_or_nonconverged_run",
            "error_owner": self.owner, "error_code": self.code,
            "error_detail": self.detail, "candidate_branch_opened": False,
            "calibration_result_earned": False, "physical_result_earned": False,
        }


def _sanitize(value: object) -> str:
    text = " ".join(str(value).replace(str(ROOT), "<repo>").split())
    return text.encode("ascii", "backslashreplace").decode("ascii")[:_MAX_DETAIL] or "no detail"


def _fail(owner: str, code: str, detail: object) -> NoReturn:
    raise TDG8RCV3REC1RunnerError(owner, code, str(detail))


def _from_exception(owner: str, code: str, exc: BaseException) -> NoReturn:
    _fail(owner, code, f"{type(exc).__name__}: {exc}")


def _fixed_root(root: Path) -> Path:
    supplied = Path(os.path.abspath(os.fspath(root)))
    try:
        metadata = supplied.lstat()
    except OSError as exc:
        _from_exception("provenance", "repository_absent", exc)
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _fail("provenance", "repository_unsafe", "repository is not a real directory")
    if supplied.resolve() != ROOT.resolve():
        _fail("authority", "repository_scope", "runner requires its fixed repository")
    return ROOT.resolve()


def _read_leaf(root: Path, relative: str, maximum: int) -> bytes:
    candidate = Path(relative)
    if candidate.is_absolute() or any(part in {"", ".", ".."} for part in candidate.parts):
        raise ValueError("authority leaf path is unsafe")
    path = root / candidate
    if path.resolve(strict=False) != path:
        raise ValueError("authority leaf is redirected")
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_size > maximum:
        raise ValueError("authority leaf is unsafe")
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        active = os.fstat(descriptor)
        identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
        if identity != (active.st_dev, active.st_ino, active.st_size, active.st_mtime_ns):
            raise ValueError("authority leaf raced")
        raw = b""
        while len(raw) <= maximum:
            block = os.read(descriptor, min(1 << 20, maximum + 1 - len(raw)))
            if not block:
                break
            raw += block
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    try:
        path_after = path.lstat()
    except OSError as exc:
        raise ValueError("authority leaf disappeared after read") from exc
    if (
        len(raw) != before.st_size or len(raw) > maximum
        or identity != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
        or (path_after.st_dev, path_after.st_ino, path_after.st_size, path_after.st_mtime_ns) != identity
    ):
        raise ValueError("authority leaf size differs")
    return raw


def _execution_authority(root: Path, commit: str) -> authority.TDG8RCV3REC1Authority:
    if not isinstance(commit, str) or not _COMMIT.fullmatch(commit):
        _fail("authority", "commit_malformed", "authority commit is not a full SHA")
    try:
        config = _read_leaf(root, authority.CONFIG_PATH, 1024 * 1024)
        result = _read_leaf(root, authority.RESULT_PATH, 4 * 1024 * 1024)
        receipt = authority.authorize_execution(root, config, result, commit)
    except TDG8RCV3REC1RunnerError:
        raise
    except Exception as exc:
        _from_exception("authority", "committed_image_rejected", exc)
    if (
        not isinstance(receipt, authority.TDG8RCV3REC1Authority)
        or receipt.authority_commit != commit
        or receipt.prior_auth1_commit != authority.PRIOR_AUTH1_COMMIT
        or receipt.normalization_fix_commit != authority.NORMALIZATION_FIX_COMMIT
        or receipt.original_execution_commit != ORIGINAL_EXECUTION_COMMIT
        or receipt.original_plan_sha256 != ORIGINAL_PLAN_SHA256
        or receipt.campaign_id != authority.CAMPAIGN_ID
        or receipt.destination_path != EXPECTED_DESTINATION_PATH.as_posix()
        or receipt.predicted_checkpoint_sha256 != authority.PREDICTED_CHECKPOINT_SHA256
        or not receipt.candidate_branches_forbidden
    ):
        _fail("authority", "typed_scope_rejected", "REC1 typed authority differs")
    return receipt


def _original_plan(root: Path) -> Any:
    try:
        plan = replace(
            construct_first_event(root, authorization_commit=ORIGINAL_EXECUTION_COMMIT),
            campaign_id=authority.CAMPAIGN_ID,
        )
    except Exception as exc:
        _from_exception("authority", "original_plan_unavailable", exc)
    if (
        plan.sha256 != ORIGINAL_PLAN_SHA256 or plan.authorization_commit != ORIGINAL_EXECUTION_COMMIT
        or plan.protocol_artifact_id != EXPECTED_PROTOCOL or plan.campaign_id != authority.CAMPAIGN_ID
        or plan.branch != "GR-0" or plan.amplitude != "3"
        or plan.common_event_index != authority.EVENT
        or plan.start_time.rational != "23/16" or plan.target_time.rational != "3/2"
        or tuple(member.member_key for member in plan.members) != MEMBER_KEYS
    ):
        _fail("authority", "original_plan_rejected", "immutable GR-0 plan differs")
    return plan


def _destination(root: Path) -> Path:
    if Path(authority.DESTINATION_PATH) != EXPECTED_DESTINATION_PATH:
        _fail("authority", "destination_scope", "REC1 destination differs")
    path = root / EXPECTED_DESTINATION_PATH
    if path.resolve(strict=False) != path:
        _fail("provenance", "destination_redirected", "campaign destination is redirected")
    metadata = path.lstat()
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _fail("provenance", "destination_unsafe", "campaign destination is unsafe")
    return path


def _summary(state: str, **extra: object) -> dict[str, Any]:
    return {
        "schema": SCHEMA, "runner_id": RUNNER_ID, "state": state,
        "projection_id": authority.PROJECTION_ID,
        "candidate_branch_opened": False, "calibration_result_earned": False,
        "physical_result_earned": False, **extra,
    }


def _boundary(root: Path, receipt: authority.TDG8RCV3REC1Authority) -> authority.REC1Boundary:
    try:
        return authority.inspect_boundary(root, receipt)
    except Exception as exc:
        _from_exception("provenance", "boundary_rejected", exc)


def status(root: Path, *, authority_commit: str) -> dict[str, Any]:
    repository = _fixed_root(root)
    receipt = _execution_authority(repository, authority_commit)
    _original_plan(repository)
    boundary = _boundary(repository, receipt)
    return _summary(
        boundary.phase.value, authority_commit=receipt.authority_commit,
        original_execution_commit=receipt.original_execution_commit,
        original_plan_sha256=receipt.original_plan_sha256, boundary=asdict(boundary),
    )


def _checkpoint(plan: Any, store: HLT16CampaignStore) -> Any:
    try:
        checkpoint = runtime.recover_persisted_generation_one(plan, store)
    except Exception as exc:
        _from_exception("provenance", "persisted_restore_rejected", exc)
    if (
        checkpoint.event not in {authority.EVENT, authority.SUCCESSOR_EVENT}
        or checkpoint.authorization_commit != ORIGINAL_EXECUTION_COMMIT
        or checkpoint.plan_sha256 != ORIGINAL_PLAN_SHA256
    ):
        _fail("scope", "checkpoint_scope", "persisted checkpoint escaped REC1 scope")
    return checkpoint


def _require_generation10_barrier(
    destination: Path,
    receipt: authority.TDG8RCV3REC1Authority,
    *,
    owned_writer_pid: int | None = None,
) -> authority.REC1Boundary:
    try:
        boundary = authority.inspect_store(
            destination, receipt, owned_writer_pid=owned_writer_pid,
        )
    except Exception as exc:
        _from_exception("recovery", "predicted_generation10_rejected", exc)
    if (
        boundary.checkpoint_generation != authority.PREDICTED_GENERATION
        or boundary.checkpoint_sha256 != authority.PREDICTED_CHECKPOINT_SHA256
        or not boundary.exact_generation10_ancestor
        or not boundary.physical_state_unchanged_at_recovery
    ):
        _fail("recovery", "predicted_generation10_rejected", "exact barrier is absent")
    return boundary


def _take_over(store: HLT16CampaignStore, checkpoint: Any) -> str:
    try:
        return store.take_over_stale_writer(checkpoint, host=os.uname().nodename)
    except HLT16CampaignStoreError as exc:
        _from_exception("writer", "stale_takeover_rejected", exc)


def _release(store: HLT16CampaignStore, token: str) -> None:
    try:
        store.release_active_writer(
            owner_token=token, host=os.uname().nodename, pid=os.getpid(),
        )
    except HLT16CampaignStoreError as exc:
        _from_exception("writer", "recovery_handoff_rejected", exc)


def _reconcile_to_clean_generation10(
    plan: Any,
    store: HLT16CampaignStore,
    destination: Path,
    receipt: authority.TDG8RCV3REC1Authority,
    boundary: authority.REC1Boundary,
) -> authority.REC1Boundary:
    """Finish only the exact recovery edge and return a lease-free gen10."""
    checkpoint = _checkpoint(plan, store)
    if boundary.phase in {
        authority.REC1Phase.EXACT_PRE_RECOVERY,
        authority.REC1Phase.GEN9_RECONCILE_PENDING,
    }:
        if checkpoint.generation != authority.ENTRY_CHECKPOINT_GENERATION:
            _fail("recovery", "entry_moved", "generation-nine entry moved before takeover")
        token = _take_over(store, checkpoint)
        try:
            decision = recovery.reconcile(store)
            if (
                decision.disposition != "reconciled_checkpoint"
                or decision.checkpoint_sha256 != authority.PREDICTED_CHECKPOINT_SHA256
            ):
                _fail("recovery", "reconciliation_result", "suffix did not yield predicted gen10")
            barrier = _require_generation10_barrier(
                destination, receipt, owned_writer_pid=os.getpid(),
            )
            if barrier.phase is not authority.REC1Phase.GEN10_HANDOFF_PENDING:
                _fail("recovery", "handoff_phase", "predicted gen10 lacks recovery lease")
            _release(store, token)
        except BaseException:
            # The lease deliberately survives an interrupted or rejected edge.
            # Its exact bytes make the next invocation classify one crash cut.
            raise
    elif boundary.phase is authority.REC1Phase.GEN10_HANDOFF_PENDING:
        if checkpoint.sha256 != authority.PREDICTED_CHECKPOINT_SHA256:
            _fail("recovery", "handoff_checkpoint", "handoff does not bind predicted gen10")
        token = _take_over(store, checkpoint)
        _require_generation10_barrier(destination, receipt, owned_writer_pid=os.getpid())
        _release(store, token)
    elif boundary.phase is authority.REC1Phase.EXACT_CLEAN_GEN10:
        return boundary
    else:
        _fail("scope", "reconcile_phase", f"cannot reconcile phase {boundary.phase.value}")
    clean = _require_generation10_barrier(destination, receipt)
    if clean.phase is not authority.REC1Phase.EXACT_CLEAN_GEN10 or clean.active_write:
        _fail("recovery", "clean_handoff_absent", "generation-ten lease did not close")
    return clean


def _validate_event_result(result: object) -> Mapping[str, Any]:
    if not isinstance(result, Mapping):
        _fail("runtime", "result_shape", "event engine returned no mapping")
    if (
        result.get("candidate_branch_opened") is not False
        or result.get("calibration_result_earned") is not False
        or result.get("physical_result_earned") is not False
        or result.get("authorization_commit") != ORIGINAL_EXECUTION_COMMIT
        or result.get("campaign_id") != authority.CAMPAIGN_ID
        or result.get("plan_sha256") != ORIGINAL_PLAN_SHA256
        or result.get("event") not in {authority.EVENT, authority.SUCCESSOR_EVENT}
    ):
        _fail("scope", "result_promotion", "event result crossed REC1 scope")
    if result.get("event") == authority.SUCCESSOR_EVENT:
        if result.get("state") != "event_one_complete" or result.get("target") != authority.SUCCESSOR_TARGET:
            _fail("scope", "event_boundary", "event-24 result differs")
    elif result.get("state") not in {"scientific_terminal", "invalid_terminal"}:
        _fail("scope", "event_boundary", "event-23 result is not typed terminal")
    return result


def run(root: Path, *, authority_commit: str) -> dict[str, Any]:
    repository = _fixed_root(root)
    receipt = _execution_authority(repository, authority_commit)
    plan = _original_plan(repository)
    boundary = _boundary(repository, receipt)
    if boundary.phase in {authority.REC1Phase.EVENT24_COMPLETE, authority.REC1Phase.TYPED_TERMINAL}:
        return _summary(
            boundary.phase.value, authority_commit=receipt.authority_commit,
            boundary=asdict(boundary), output_created=False,
        )
    destination = _destination(repository)
    store = HLT16CampaignStore(destination)
    if boundary.phase in {
        authority.REC1Phase.EXACT_PRE_RECOVERY,
        authority.REC1Phase.GEN9_RECONCILE_PENDING,
        authority.REC1Phase.GEN10_HANDOFF_PENDING,
        authority.REC1Phase.EXACT_CLEAN_GEN10,
    }:
        boundary = _reconcile_to_clean_generation10(plan, store, destination, receipt, boundary)
    elif boundary.phase is not authority.REC1Phase.LAWFUL_EVENT23_DESCENDANT:
        _fail("scope", "entry_phase", f"unrecognized REC1 phase {boundary.phase.value}")

    # This is the sole call into the PDE-capable engine.  The exact generation
    # ten barrier above dominates every path reaching it.
    checkpoint = _checkpoint(plan, store)
    if checkpoint.generation < authority.PREDICTED_GENERATION:
        _fail("recovery", "pre_barrier_engine", "event engine cannot see generation nine")
    if not boundary.exact_generation10_ancestor:
        boundary = _require_generation10_barrier(destination, receipt)
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
        _from_exception("runtime", "same_event_resume_rejected", exc)
    except Exception as exc:
        _from_exception("runtime", "unexpected_same_event_failure", exc)
    checked = _validate_event_result(result)
    return _summary(
        str(checked["state"]), authority_commit=receipt.authority_commit,
        original_execution_commit=receipt.original_execution_commit,
        original_plan_sha256=receipt.original_plan_sha256, result=dict(checked),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority-commit", required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--status", action="store_true")
    mode.add_argument("--run", action="store_true")
    args = parser.parse_args(argv)
    try:
        result = status(ROOT, authority_commit=args.authority_commit) if args.status else run(
            ROOT, authority_commit=args.authority_commit,
        )
    except TDG8RCV3REC1RunnerError as exc:
        print(json.dumps(exc.as_record(), sort_keys=True, separators=(",", ":")), file=sys.stderr)
        return 2
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
