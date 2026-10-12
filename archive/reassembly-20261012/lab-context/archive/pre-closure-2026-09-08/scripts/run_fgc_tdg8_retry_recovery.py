#!/usr/bin/env python3
"""Inspect, recover, or resume the exact TDG8 retry boundary.

``status`` is read-only.  ``recover`` may only turn the authenticated
generation-seven TDG6 suffix into its unique metadata-only generation-eight
checkpoint.  ``resume`` starts from that exact recovered ancestor and advances
only the original GR-0 event.  No mode can import a candidate branch, reopen
generation-zero arrays, or re-observe destination absence.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import stat
import sys
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = str(ROOT / "src")
sys.path[:] = [entry for entry in sys.path if entry != _SOURCE_ROOT]
sys.path.insert(0, str(ROOT))
from scripts import run_fgc_pro19_event1 as event1  # noqa: E402

sys.path[:] = [entry for entry in sys.path if entry not in {_SOURCE_ROOT, str(ROOT)}]
sys.path[:0] = [_SOURCE_ROOT, str(ROOT)]

from recursive_horizons.fgc.evolution import hlt16_campaign_recovery as recovery  # noqa: E402
from recursive_horizons.fgc.evolution import hlt16_campaign_runtime as runtime  # noqa: E402
from recursive_horizons.fgc.evolution import tdg8_retry_recovery_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution import tdg8_successor_runtime as tdg8  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
    HLT16CampaignStoreError,
)
from recursive_horizons.fgc.evolution.hlt16_state_store import (  # noqa: E402
    HLT16StateStoreError,
    read_nofollow,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto17_pure_construction import (  # noqa: E402
    MEMBER_KEYS,
)
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    build_static_gr0_shells,
)


RUNNER_ID = "FGC-1-TDG8-RCV1-RUNNER"
SCHEMA = "FGC-1-TDG8-RCV1-runner-v1"
_COMMIT = re.compile(r"[0-9a-f]{40}\Z")


class TDG8RetryRecoveryRunnerError(RuntimeError):
    """The exact recovery/resume authority or campaign lineage differs."""


def _fixed_root(root: Path) -> Path:
    supplied = Path(os.path.abspath(os.fspath(root)))
    try:
        metadata = supplied.lstat()
    except OSError as exc:
        raise TDG8RetryRecoveryRunnerError("repository root is absent") from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise TDG8RetryRecoveryRunnerError("repository root is unsafe")
    if supplied.resolve() != ROOT.resolve():
        raise TDG8RetryRecoveryRunnerError("RCV1 runner requires its fixed repository")
    return ROOT.resolve()


def _destination(root: Path) -> Path:
    path = root / authority.DESTINATION_PATH
    if path.resolve(strict=False) != path:
        raise TDG8RetryRecoveryRunnerError("RCV1 destination is redirected")
    try:
        metadata = path.lstat()
    except OSError as exc:
        raise TDG8RetryRecoveryRunnerError("RCV1 destination is absent") from exc
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        raise TDG8RetryRecoveryRunnerError("RCV1 destination is unsafe")
    return path


def _execution_authority(
    root: Path, commit: str,
) -> authority.TDG8RetryRecoveryExecutionAuthority:
    if not isinstance(commit, str) or not _COMMIT.fullmatch(commit):
        raise TDG8RetryRecoveryRunnerError("authority commit is malformed")
    try:
        config = authority._read_leaf(
            root, authority.CONFIG_PATH, "RCV1 execution config", 1024 * 1024
        )
        result = authority._read_leaf(
            root, authority.RESULT_PATH, "RCV1 execution result", 4 * 1024 * 1024
        )
        receipt = authority.authorize_execution(root, config, result, commit)
    except authority.TDG8RetryRecoveryAuthorityError as exc:
        raise TDG8RetryRecoveryRunnerError("RCV1 execution authority failed") from exc
    if (
        receipt.authority_commit != commit
        or receipt.original_execution_commit != authority.ORIGINAL_EXECUTION_COMMIT
        or receipt.original_plan_sha256 != authority.ORIGINAL_PLAN_SHA256
        or receipt.campaign_id != authority.CAMPAIGN_ID
        or receipt.destination_path != authority.DESTINATION_PATH
        or receipt.anchor_checkpoint_sha256 != authority.ANCHOR_CHECKPOINT_SHA256
        or receipt.anchor_journal_sha256 != authority.ANCHOR_JOURNAL_SHA256
        or receipt.expected_checkpoint_sha256 != authority.EXPECTED_CHECKPOINT_SHA256
        or receipt.expected_journal_sha256 != authority.EXPECTED_JOURNAL_SHA256
        or not receipt.candidate_branches_forbidden
    ):
        raise TDG8RetryRecoveryRunnerError("RCV1 typed authority scope differs")
    return receipt


def _original_plan(root: Path) -> Any:
    try:
        plan = tdg8.construct_successor_plan(root, authority.ORIGINAL_EXECUTION_COMMIT)
    except Exception as exc:
        raise TDG8RetryRecoveryRunnerError("original TDG8 plan cannot be reconstructed") from exc
    if (
        plan.sha256 != authority.ORIGINAL_PLAN_SHA256
        or plan.authorization_commit != authority.ORIGINAL_EXECUTION_COMMIT
        or plan.protocol_artifact_id != authority.TARGET_PROTOCOL
        or plan.campaign_id != authority.CAMPAIGN_ID
        or plan.branch != "GR-0"
        or plan.amplitude != "3"
        or plan.common_event_index != 23
        or plan.start_time.rational != "23/16"
        or plan.target_time.rational != "3/2"
        or tuple(member.member_key for member in plan.members) != MEMBER_KEYS
    ):
        raise TDG8RetryRecoveryRunnerError("original TDG8 plan differs")
    return plan


def _expected_generation_eight(store: HLT16CampaignStore, *, require_no_writer: bool) -> Any:
    try:
        snapshot = store.authenticated_snapshot()
        status = store.inspect_recovery()
    except HLT16CampaignStoreError as exc:
        raise TDG8RetryRecoveryRunnerError("recovered store does not validate") from exc
    checkpoint = snapshot.checkpoint
    member = checkpoint.members[authority.MEMBER_KEY]
    if (
        checkpoint.authorization_commit != authority.ORIGINAL_EXECUTION_COMMIT
        or checkpoint.plan_sha256 != authority.ORIGINAL_PLAN_SHA256
        or checkpoint.campaign_id != authority.CAMPAIGN_ID
        or checkpoint.protocol != authority.TARGET_PROTOCOL
        or checkpoint.event != 23
        or checkpoint.target
        != {"rational": "3/2", "binary64_hex": "0x1.8000000000000p+0"}
        or checkpoint.generation != authority.EXPECTED_CHECKPOINT_GENERATION
        or checkpoint.sha256 != authority.EXPECTED_CHECKPOINT_SHA256
        or checkpoint.journal_sequence != authority.EXPECTED_JOURNAL_SEQUENCE
        or checkpoint.journal_tip_sha256 != authority.EXPECTED_JOURNAL_SHA256
        or checkpoint.disposition != "nonterminal"
        or snapshot.suffix_classification != "clean"
        or snapshot.suffix_records
        or snapshot.staging_paths
        or snapshot.orphan_descriptor_sha256 is not None
        or snapshot.orphan_payload_semantic_sha256 is not None
        or snapshot.terminal_lock_present
        or member.cursor["mode"] != "RETRY_PENDING"
        or member.pending_owner != "temporal"
        or member.pending_cap_hex != authority.EXPECTED_PENDING_CAP_HEX
    ):
        raise TDG8RetryRecoveryRunnerError("exact recovered generation eight differs")
    if require_no_writer and (
        status.active_write or status.writer_state is not None or status.state != "clean_checkpoint"
    ):
        raise TDG8RetryRecoveryRunnerError("recovered checkpoint still has a writer")
    try:
        persisted = store.load_state(member.descriptor_sha256)
        physical = array_content_sha256(
            persisted.arrays["u"], persisted.arrays["p"], persisted.arrays["q"]
        )
    except (HLT16CampaignStoreError, KeyError, TypeError) as exc:
        raise TDG8RetryRecoveryRunnerError("recovered physical state cannot be restored") from exc
    if physical != authority.ANCHOR_PHYSICAL_STATE_SHA256:
        raise TDG8RetryRecoveryRunnerError("recovery advanced or replaced the physical state")
    return checkpoint


def _require_generation_eight_ancestor(store: HLT16CampaignStore) -> None:
    relative = (
        f"checkpoints/{authority.EXPECTED_CHECKPOINT_GENERATION:020d}-"
        f"{authority.EXPECTED_CHECKPOINT_SHA256}.json"
    )
    try:
        raw = read_nofollow(store.root, relative, "RCV1 generation-eight ancestor")
        value = json.loads(raw.decode("utf-8"))
    except (HLT16StateStoreError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TDG8RetryRecoveryRunnerError("exact generation-eight ancestor is absent") from exc
    if (
        not isinstance(value, Mapping)
        or value.get("checkpoint_sha256") != authority.EXPECTED_CHECKPOINT_SHA256
        or value.get("generation") != authority.EXPECTED_CHECKPOINT_GENERATION
        or value.get("journal_tip_sha256") != authority.EXPECTED_JOURNAL_SHA256
        or value.get("authorization_commit") != authority.ORIGINAL_EXECUTION_COMMIT
        or value.get("plan_sha256") != authority.ORIGINAL_PLAN_SHA256
    ):
        raise TDG8RetryRecoveryRunnerError("generation-eight ancestor bytes differ")


def _require_same_event_checkpoint(checkpoint: Any) -> None:
    """Reject every lineage position outside the authorized event-23 edge."""

    if checkpoint.event == 23:
        if checkpoint.target != {
            "rational": "3/2",
            "binary64_hex": "0x1.8000000000000p+0",
        }:
            raise TDG8RetryRecoveryRunnerError("RCV1 event-23 target differs")
        return
    if checkpoint.event == 24:
        if (
            checkpoint.disposition != "event_complete"
            or checkpoint.target
            != {
                "rational": "25/16",
                "binary64_hex": "0x1.9000000000000p+0",
            }
        ):
            raise TDG8RetryRecoveryRunnerError(
                "RCV1 completed-event boundary differs"
            )
        return
    raise TDG8RetryRecoveryRunnerError("RCV1 lineage escaped the one-event scope")


def _summary(state: str, **extra: object) -> dict[str, Any]:
    return {
        "schema": SCHEMA, "runner_id": RUNNER_ID, "state": state,
        "candidate_branch_opened": False, "calibration_result_earned": False,
        "physical_result_earned": False, **extra,
    }


def status(root: Path, *, authority_commit: str) -> dict[str, Any]:
    repository = _fixed_root(root)
    receipt = _execution_authority(repository, authority_commit)
    store = HLT16CampaignStore(_destination(repository))
    try:
        snapshot = store.authenticated_snapshot()
        recovery_status = store.inspect_recovery()
    except HLT16CampaignStoreError as exc:
        raise TDG8RetryRecoveryRunnerError("RCV1 status cannot validate the store") from exc
    checkpoint = snapshot.checkpoint
    _require_same_event_checkpoint(checkpoint)
    if checkpoint.sha256 == authority.ANCHOR_CHECKPOINT_SHA256:
        try:
            boundary = authority.inspect_recovery_boundary(
                repository,
                permitted_writer_states=frozenset(
                    {"stale_verified_lock", "live_local_writer"}
                ),
            )
        except authority.TDG8RetryRecoveryAuthorityError as exc:
            raise TDG8RetryRecoveryRunnerError(
                "RCV1 status anchor differs"
            ) from exc
        state = str(boundary["state"])
    elif checkpoint.sha256 == authority.EXPECTED_CHECKPOINT_SHA256:
        if recovery_status.active_write:
            try:
                authority.inspect_recovery_boundary(
                    repository,
                    permitted_writer_states=frozenset(
                        {"stale_verified_lock", "live_local_writer"}
                    ),
                )
            except authority.TDG8RetryRecoveryAuthorityError as exc:
                raise TDG8RetryRecoveryRunnerError(
                    "RCV1 recovered writer boundary differs"
                ) from exc
            state = "generation8_writer_cleanup_pending"
        else:
            _expected_generation_eight(store, require_no_writer=True)
            state = "recovered_retry_ready"
    else:
        _require_generation_eight_ancestor(store)
        if checkpoint.authorization_commit != authority.ORIGINAL_EXECUTION_COMMIT or checkpoint.plan_sha256 != authority.ORIGINAL_PLAN_SHA256:
            raise TDG8RetryRecoveryRunnerError("post-recovery lineage changed authority")
        state = "same_event_descendant"
    return _summary(
        state,
        authority_commit=receipt.authority_commit,
        original_plan_sha256=receipt.original_plan_sha256,
        checkpoint_generation=checkpoint.generation,
        checkpoint_sha256=checkpoint.sha256,
        suffix_kinds=list(snapshot.suffix_kinds),
        recovery=asdict(recovery_status),
    )


def recover(root: Path, *, authority_commit: str) -> dict[str, Any]:
    repository = _fixed_root(root)
    receipt = _execution_authority(repository, authority_commit)
    _original_plan(repository)
    destination = _destination(repository)
    store = HLT16CampaignStore(destination)
    try:
        boundary = authority.inspect_recovery_boundary(repository)
    except authority.TDG8RetryRecoveryAuthorityError as exc:
        raise TDG8RetryRecoveryRunnerError("RCV1 recovery anchor differs") from exc
    try:
        snapshot = store.authenticated_snapshot()
    except HLT16CampaignStoreError as exc:
        raise TDG8RetryRecoveryRunnerError(
            "RCV1 recovery snapshot differs"
        ) from exc
    predecessor = snapshot.checkpoint
    try:
        token = event1._acquire_writer(store, predecessor)
    except Exception as exc:
        raise TDG8RetryRecoveryRunnerError(
            "RCV1 stale-writer takeover failed"
        ) from exc
    clean_release = False
    try:
        # Lease takeover changes no campaign fact.  Re-read the complete anchor
        # and unique preview under the live capability immediately before the
        # sole metadata publication.
        owned_boundary = authority.inspect_recovery_boundary(
            repository, permitted_writer_states=frozenset({"live_local_writer"})
        )
        if owned_boundary["state"] != boundary["state"]:
            raise TDG8RetryRecoveryRunnerError(
                "RCV1 recovery boundary moved during writer takeover"
            )
        if boundary["state"] in {
            "generation7_rejection_only",
            "generation7_transition_published",
        }:
            decision = recovery.reconcile(store)
            if (
                decision.disposition != "reconciled_checkpoint"
                or not decision.published
                or decision.checkpoint_sha256
                != authority.EXPECTED_CHECKPOINT_SHA256
            ):
                raise TDG8RetryRecoveryRunnerError(
                    "RCV1 recovery published a different edge"
                )
        elif boundary["state"] != "generation8_checkpoint_published":
            raise TDG8RetryRecoveryRunnerError(
                "RCV1 recovery phase is not authorized"
            )
        # Once the recovery owner reports the exact checkpoint address, a
        # clean store must not retain our live lease even when the independent
        # post-publication identity check itself rejects.  The release helper
        # still refuses every incomplete or ambiguous suffix.
        clean_release = True
        checkpoint = _expected_generation_eight(store, require_no_writer=False)
    except (
        recovery.HLT16RecoveryError,
        authority.TDG8RetryRecoveryAuthorityError,
    ) as exc:
        raise TDG8RetryRecoveryRunnerError("RCV1 metadata recovery failed") from exc
    finally:
        if clean_release and not event1._release_writer_if_closed(store, owner_token=token):
            raise TDG8RetryRecoveryRunnerError("RCV1 recovery writer did not close")
    _expected_generation_eight(store, require_no_writer=True)
    return _summary(
        "metadata_recovery_complete",
        authority_commit=receipt.authority_commit,
        checkpoint_generation=checkpoint.generation,
        checkpoint_sha256=checkpoint.sha256,
        journal_sequence=checkpoint.journal_sequence,
        journal_sha256=checkpoint.journal_tip_sha256,
        physical_state_advanced=False,
    )


def resume(root: Path, *, authority_commit: str) -> dict[str, Any]:
    repository = _fixed_root(root)
    receipt = _execution_authority(repository, authority_commit)
    plan = _original_plan(repository)
    store = HLT16CampaignStore(_destination(repository))
    try:
        snapshot = store.authenticated_snapshot()
    except HLT16CampaignStoreError as exc:
        raise TDG8RetryRecoveryRunnerError("RCV1 resume snapshot differs") from exc
    _require_same_event_checkpoint(snapshot.checkpoint)
    if snapshot.checkpoint.sha256 == authority.ANCHOR_CHECKPOINT_SHA256:
        raise TDG8RetryRecoveryRunnerError("metadata recovery must complete before resume")
    _require_generation_eight_ancestor(store)
    if snapshot.checkpoint.sha256 == authority.EXPECTED_CHECKPOINT_SHA256:
        _expected_generation_eight(store, require_no_writer=True)
    try:
        checkpoint = runtime.recover_persisted_generation_one(plan, store)
        templates = build_static_gr0_shells(repository)
    except Exception as exc:
        raise TDG8RetryRecoveryRunnerError("RCV1 persisted resume preflight failed") from exc
    if tuple(templates) != MEMBER_KEYS:
        raise TDG8RetryRecoveryRunnerError("RCV1 static GR-0 shell order differs")
    try:
        result = event1._advance_authenticated_event(
            plan, store, checkpoint, templates
        )
    except Exception as exc:
        raise TDG8RetryRecoveryRunnerError(
            "RCV1 same-event resume failed"
        ) from exc
    if (
        not isinstance(result, Mapping)
        or result.get("candidate_branch_opened") is not False
        or result.get("calibration_result_earned") is not False
        or result.get("physical_result_earned") is not False
        or result.get("campaign_id") != authority.CAMPAIGN_ID
        or result.get("plan_sha256") != authority.ORIGINAL_PLAN_SHA256
    ):
        raise TDG8RetryRecoveryRunnerError("RCV1 resume crossed its result scope")
    return _summary(
        str(result.get("state")),
        authority_commit=receipt.authority_commit,
        original_plan_sha256=receipt.original_plan_sha256,
        result=dict(result),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority-commit", required=True)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--status", action="store_true")
    mode.add_argument("--recover", action="store_true")
    mode.add_argument("--resume", action="store_true")
    args = parser.parse_args(argv)
    try:
        if args.status:
            result = status(ROOT, authority_commit=args.authority_commit)
        elif args.recover:
            result = recover(ROOT, authority_commit=args.authority_commit)
        else:
            result = resume(ROOT, authority_commit=args.authority_commit)
    except TDG8RetryRecoveryRunnerError as exc:
        print(str(exc), file=sys.stderr)
        return 2
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
