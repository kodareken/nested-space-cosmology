#!/usr/bin/env python3
"""Inspect or continue the authorized bounded TDG8 retry chain.

Status is read-only. Run authenticates the committed RCV2 image, restores only
persisted GR-0 state, reconciles recognized metadata/publication suffixes, and
continues no farther than event 23 -> 24. It cannot import GEN0 or candidate
branches.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
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
from recursive_horizons.fgc.evolution import tdg8_bounded_retry_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution import tdg8_successor_runtime as tdg8  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
    HLT16CampaignStoreError,
)
from recursive_horizons.fgc.evolution.proto17_pure_construction import MEMBER_KEYS  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    build_static_gr0_shells,
)


RUNNER_ID = "FGC-1-TDG8-RCV2-RUNNER"
SCHEMA = "FGC-1-TDG8-RCV2-runner-v1"
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")
_MAX_DETAIL = 480


class TDG8BoundedRetryRunnerError(RuntimeError):
    """A fail-closed error with a stable owner and useful bounded detail."""

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
    raise TDG8BoundedRetryRunnerError(owner, code, str(detail))


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


def _destination(root: Path) -> Path:
    path = root / authority.DESTINATION_PATH
    if path.resolve(strict=False) != path:
        _fail("provenance", "destination_redirected", "campaign destination is redirected")
    try:
        metadata = path.lstat()
    except OSError as exc:
        _from_exception("provenance", "destination_absent", exc)
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _fail("provenance", "destination_unsafe", "campaign destination is unsafe")
    return path


def _execution_authority(
    root: Path, commit: str,
) -> authority.TDG8BoundedRetryExecutionAuthority:
    if not isinstance(commit, str) or not _COMMIT.fullmatch(commit):
        _fail("authority", "commit_malformed", "authority commit is not a full Git SHA")
    try:
        config = authority.rcv1._read_leaf(
            root, authority.CONFIG_PATH, "RCV2 execution config", 1024 * 1024,
        )
        result = authority.rcv1._read_leaf(
            root, authority.RESULT_PATH, "RCV2 execution result", 4 * 1024 * 1024,
        )
        receipt = authority.authorize_execution(root, config, result, commit)
    except Exception as exc:
        _from_exception("authority", "committed_image_rejected", exc)
    if (
        receipt.authority_commit != commit
        or receipt.original_execution_commit != authority.ORIGINAL_EXECUTION_COMMIT
        or receipt.original_plan_sha256 != authority.ORIGINAL_PLAN_SHA256
        or receipt.campaign_id != authority.CAMPAIGN_ID
        or receipt.destination_path != authority.DESTINATION_PATH
        or receipt.entry_checkpoint_sha256 != authority.ENTRY_CHECKPOINT_SHA256
        or receipt.entry_rejection_sha256 != authority.ENTRY_REJECTION_SHA256
        or receipt.entry_transition_sha256 != authority.ENTRY_TRANSITION_SHA256
        or receipt.first_recovery_checkpoint_sha256
        != authority.FIRST_RECOVERY_CHECKPOINT_SHA256
        or receipt.temporal_retry_cap != authority.TEMPORAL_RETRY_CAP
        or receipt.minimum_macro_step_hex != authority.MINIMUM_MACRO_STEP_HEX
        or not receipt.candidate_branches_forbidden
    ):
        _fail("authority", "typed_scope_rejected", "typed execution authority differs")
    return receipt


def _original_plan(root: Path) -> Any:
    try:
        plan = tdg8.construct_successor_plan(
            root, authority.ORIGINAL_EXECUTION_COMMIT,
        )
    except Exception as exc:
        _from_exception("authority", "original_plan_unavailable", exc)
    if (
        plan.sha256 != authority.ORIGINAL_PLAN_SHA256
        or plan.authorization_commit != authority.ORIGINAL_EXECUTION_COMMIT
        or plan.protocol_artifact_id != authority.TARGET_PROTOCOL
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


def _summary(state: str, **extra: object) -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "runner_id": RUNNER_ID,
        "state": state,
        "candidate_branch_opened": False,
        "calibration_result_earned": False,
        "physical_result_earned": False,
        **extra,
    }


def status(root: Path, *, authority_commit: str) -> dict[str, Any]:
    repository = _fixed_root(root)
    receipt = _execution_authority(repository, authority_commit)
    _original_plan(repository)
    try:
        boundary = authority.inspect_descendant(repository, receipt)
    except authority.TDG8BoundedRetryAuthorityError as exc:
        _from_exception("provenance", "descendant_rejected", exc)
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
        or result.get("campaign_id") != authority.CAMPAIGN_ID
        or result.get("plan_sha256") != authority.ORIGINAL_PLAN_SHA256
        or result.get("event") not in {authority.EVENT, authority.SUCCESSOR_EVENT}
    ):
        _fail("scope", "result_promotion", "same-event result crossed RCV2 scope")
    if result.get("event") == authority.SUCCESSOR_EVENT and (
        result.get("state") != "event_one_complete"
        or result.get("target") != authority.SUCCESSOR_TARGET
    ):
        _fail("scope", "event_boundary", "event-24 result is not the exact completion boundary")
    if result.get("event") == authority.EVENT and result.get("state") not in {
        "scientific_terminal", "invalid_terminal",
    }:
        _fail("scope", "event_boundary", "event-23 result is not a typed terminal")
    return result


def run(root: Path, *, authority_commit: str) -> dict[str, Any]:
    repository = _fixed_root(root)
    receipt = _execution_authority(repository, authority_commit)
    plan = _original_plan(repository)
    store = HLT16CampaignStore(_destination(repository))
    try:
        boundary = authority.inspect_descendant(repository, receipt)
    except authority.TDG8BoundedRetryAuthorityError as exc:
        _from_exception("provenance", "descendant_rejected", exc)
    if boundary.terminal and boundary.terminal_lock_present:
        return _summary(
            "same_event_terminal",
            authority_commit=receipt.authority_commit,
            original_plan_sha256=receipt.original_plan_sha256,
            boundary=asdict(boundary),
            output_created=False,
        )
    if boundary.event == authority.SUCCESSOR_EVENT:
        return _summary(
            "event_one_complete",
            authority_commit=receipt.authority_commit,
            original_plan_sha256=receipt.original_plan_sha256,
            boundary=asdict(boundary),
            output_created=False,
        )
    if not boundary.recoverable_under_rcv2:
        _fail(
            "writer",
            "restart_not_authorized",
            f"descendant is not safely recoverable: writer={boundary.writer_state!r}, suffix={boundary.suffix_classification!r}",
        )
    try:
        checkpoint = runtime.recover_persisted_generation_one(plan, store)
    except Exception as exc:
        _from_exception("provenance", "persisted_restore_rejected", exc)
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
        text = str(exc)
        owner = "writer" if "writer" in text.lower() or "lease" in text.lower() else "runtime"
        _from_exception(owner, "same_event_resume_rejected", exc)
    except Exception as exc:
        _from_exception("runtime", "unexpected_same_event_failure", exc)
    checked = _validate_run_result(result)
    return _summary(
        str(checked["state"]),
        authority_commit=receipt.authority_commit,
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
    except TDG8BoundedRetryRunnerError as exc:
        print(json.dumps(exc.as_record(), sort_keys=True, separators=(",", ":")), file=sys.stderr)
        return 2
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
