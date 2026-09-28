#!/usr/bin/env python3
"""Run exactly the first authenticated PROTO18 GR-0 common event.

The command is intentionally smaller than a calibration campaign.  It binds
one committed external image, reconstructs only the six frozen GR-0 members,
bridges the authenticated generation-zero store, and advances one member at a
time until all six reach ``3/2``.  Every accepted macro-step is its own durable
checkpoint.  The command stops after committing common event 24; it cannot
open SGB-L, FGC-QR, activation, DEF1, or any retained-EFT path.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from typing import TYPE_CHECKING, Any, Callable, Mapping


ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = str(ROOT / "src")
if _SOURCE_ROOT not in sys.path:
    # Direct historical inspection still needs the src-layout package root.
    # The SID3 bootstrap installs the stricter [src, repository] guarded path
    # before importing this module; never duplicate or reorder that authority.
    sys.path.insert(0, _SOURCE_ROOT)

from recursive_horizons.fgc.evolution import hlt16_campaign_recovery as recovery  # noqa: E402
from recursive_horizons.fgc.evolution import hlt16_campaign_runtime as runtime  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_schema import CampaignCheckpoint  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
    HLT16CampaignStoreError,
    _stage_target,
    _validated_stage,
    _walk_tree,
)
from recursive_horizons.fgc.evolution.hlt16_state_store import (  # noqa: E402
    SEALED_GEN0_CHECKPOINT_SHA256,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto17_pure_construction import MEMBER_KEYS  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_launch_authority import (  # noqa: E402
    LaunchAuthorityReceipt,
    Proto19LaunchAuthorityError,
    authorize_first_event,
)
from recursive_horizons.fgc.evolution.proto19_cfl_continuation_authority import (  # noqa: E402
    CFLContinuationAuthorityError,
    CFLContinuationReceipt,
    authorize_cfl_continuation,
)
from recursive_horizons.fgc.evolution import proto19_sid1_authority as sid1  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    STATIC_INPUT_PATHS,
    build_static_gr0_shells,
)
from recursive_horizons.fgc.evolution.proto12_runtime import (  # noqa: E402
    Proto12GR0EvolutionOperator,
)
from recursive_horizons.fgc.evolution.proto14_runtime import (  # noqa: E402
    Proto14RunMember,
)

if TYPE_CHECKING:
    from recursive_horizons.fgc.evolution.proto19_progression_inputs import (  # noqa: E402
        ReconstructedGR0MemberSet,
    )


RUNNER_ID = "FGC-1-PRO19-EVENT1-RUNNER"
DEFAULT_MANIFEST = Path("configs/fgc/fgc-1-pro19-launch-authority.toml")
DEFAULT_STORE = Path("runs/fgc-2-sf1/proto17/calibration")


class Proto19Event1RunnerError(RuntimeError):
    """The authenticated first event cannot progress without ambiguity."""


@dataclass(frozen=True, slots=True)
class SID3ResumePreflightContext:
    """One read-only, in-process handoff into the already authenticated event."""

    authority: Any
    store: HLT16CampaignStore
    checkpoint: CampaignCheckpoint
    templates: Mapping[str, object]
    suffix_classification: str
    recovery_state: str


def reconstruct_gr0_members(root: Path) -> "ReconstructedGR0MemberSet":
    """Import the retired raw-source bridge only on the historical GEN0 path."""
    from recursive_horizons.fgc.evolution.proto19_progression_inputs import (
        reconstruct_gr0_members as _reconstruct,
    )

    return _reconstruct(root)


def _canonical(value: object) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    )


def _git_head(root: Path) -> str:
    result = subprocess.run(
        ("git", "-C", str(root), "rev-parse", "HEAD"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        text=True,
    )
    value = result.stdout.strip()
    if result.returncode or len(value) != 40:
        raise Proto19Event1RunnerError("repository HEAD cannot be resolved")
    return value


def _accepted_time(checkpoint: CampaignCheckpoint, key: str) -> float:
    try:
        return float.fromhex(
            str(
                checkpoint.members[key].cursor["accepted_boundary_time"]["binary64_hex"]
            )
        )
    except (KeyError, TypeError, ValueError) as exc:
        raise Proto19Event1RunnerError("member accepted-time identity differs") from exc


def _target(checkpoint: CampaignCheckpoint) -> float:
    try:
        return float.fromhex(str(checkpoint.target["binary64_hex"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise Proto19Event1RunnerError("checkpoint target identity differs") from exc


def _member_to_advance(checkpoint: CampaignCheckpoint) -> str | None:
    """Choose the sole canonical member; never select from an observed outcome."""
    target = _target(checkpoint)
    pending = [
        key
        for key in MEMBER_KEYS
        if checkpoint.members[key].pending_owner is not None
        or checkpoint.members[key].cursor.get("mode") == "RETRY_PENDING"
    ]
    if len(pending) > 1:
        raise Proto19Event1RunnerError("more than one member owns a pending retry")
    if pending:
        key = pending[0]
        if _accepted_time(checkpoint, key) >= target:
            raise Proto19Event1RunnerError("target member retains a pending retry")
        return key
    for key in MEMBER_KEYS:
        accepted = _accepted_time(checkpoint, key)
        if accepted > target:
            raise Proto19Event1RunnerError("member advanced beyond the frozen target")
        if accepted < target:
            return key
        if checkpoint.members[key].cursor.get("mode") != "FRESH_READY":
            raise Proto19Event1RunnerError("target member is not fresh")
    return None


def _summary(
    checkpoint: CampaignCheckpoint,
    *,
    state: str,
    record_sha256s: tuple[str, ...] = (),
) -> dict[str, Any]:
    return {
        "schema": "FGC-1-PRO19-event1-run-status-v1",
        "runner_id": RUNNER_ID,
        "state": state,
        "authorization_commit": checkpoint.authorization_commit,
        "plan_sha256": checkpoint.plan_sha256,
        "campaign_id": checkpoint.campaign_id,
        "checkpoint_generation": checkpoint.generation,
        "checkpoint_sha256": checkpoint.sha256,
        "event": checkpoint.event,
        "target": dict(checkpoint.target),
        "member_times": {
            key: dict(checkpoint.members[key].cursor["accepted_boundary_time"])
            for key in MEMBER_KEYS
        },
        "member_modes": {
            key: checkpoint.members[key].cursor["mode"] for key in MEMBER_KEYS
        },
        "journal_sequence": checkpoint.journal_sequence,
        "journal_tip_sha256": checkpoint.journal_tip_sha256,
        "record_sha256s": list(record_sha256s),
        "candidate_branch_opened": False,
        "calibration_result_earned": False,
        "physical_result_earned": False,
    }


def authorize_event_authority(
    root: Path,
    *,
    authorization_commit: str,
    manifest_path: Path,
) -> LaunchAuthorityReceipt:
    """Close the committed image without importing any historical arrays."""
    repository = Path(root).resolve()
    try:
        relative = manifest_path.resolve().relative_to(repository).as_posix()
    except ValueError as exc:
        raise Proto19Event1RunnerError(
            "launch manifest is outside the repository"
        ) from exc
    try:
        receipt = authorize_first_event(
            repository,
            authorization_commit=authorization_commit,
            manifest_path=relative,
        )
    except (Proto19LaunchAuthorityError, ValueError, OSError) as exc:
        raise Proto19Event1RunnerError("event-one authority differs") from exc
    if receipt.progression_plan.branch != "GR-0":
        raise Proto19Event1RunnerError(
            "event-one authority opened an unauthorized branch"
        )
    return receipt


def _continuation_store_anchor(store: HLT16CampaignStore) -> dict[str, object]:
    """Independently derive the exact pre-repair recovery boundary.

    This intentionally duplicates the compact reproducer's extraction rather
    than importing it.  Agreement between the two paths is part of the
    continuation authority; neither path can mutate the store.
    """
    try:
        snapshot = store.authenticated_snapshot()
        status = store.inspect_recovery()
        checkpoint = snapshot.checkpoint
        member = checkpoint.members["RK4-2049"]
        ledger = member.ledger
        return {
            "authorization_commit": checkpoint.authorization_commit,
            "plan_sha256": checkpoint.plan_sha256,
            "campaign_id": checkpoint.campaign_id,
            "branch": "GR-0",
            "amplitude": "3",
            "event": checkpoint.event,
            "target_rational": checkpoint.target.get("rational"),
            "checkpoint_generation": checkpoint.generation,
            "checkpoint_sha256": checkpoint.sha256,
            "journal_sequence": checkpoint.journal_sequence,
            "journal_tip_sha256": checkpoint.journal_tip_sha256,
            "disposition": checkpoint.disposition,
            "suffix_classification": snapshot.suffix_classification,
            "suffix_kinds": list(snapshot.suffix_kinds),
            "staging_paths": list(snapshot.staging_paths),
            "terminal_lock_present": snapshot.terminal_lock_present,
            "writer_state": status.writer_state,
            "active_write": status.active_write,
            "member_key": "RK4-2049",
            "member_descriptor_sha256": member.descriptor_sha256,
            "member_accepted_time_hex": member.cursor["accepted_boundary_time"][
                "binary64_hex"
            ],
            "member_mode": member.cursor["mode"],
            "source_retry_total": member.source_total,
            "cfl_retry_total": member.cfl_total,
            "temporal_retry_total": ledger["cumulative_temporal_retry_count"],
            "pending_owner": member.pending_owner or "none",
            "pending_cap_hex": member.pending_cap_hex or "none",
            "durable_cfl_rejection_exists": (
                member.cfl_total != 0 or "cfl_rejection" in snapshot.suffix_kinds
            ),
            "durable_terminal_exists": (
                snapshot.terminal_lock_present
                or checkpoint.disposition in {"scientific_terminal", "invalid_terminal"}
            ),
        }
    except (HLT16CampaignStoreError, KeyError, TypeError, ValueError) as exc:
        raise Proto19Event1RunnerError(
            "continuation store anchor cannot be derived"
        ) from exc


def authorize_continuation_event_authority(
    root: Path,
    *,
    continuation_commit: str,
    manifest_path: Path,
    store: HLT16CampaignStore,
) -> CFLContinuationReceipt:
    """Bind the corrected image to only the exact generation-six edge."""
    repository = Path(root).resolve()
    try:
        relative = manifest_path.resolve().relative_to(repository).as_posix()
    except ValueError as exc:
        raise Proto19Event1RunnerError(
            "launch manifest is outside the repository"
        ) from exc
    try:
        receipt = authorize_cfl_continuation(
            repository,
            continuation_commit=continuation_commit,
            launch_manifest_path=relative,
            store_anchor=_continuation_store_anchor(store),
        )
    except (CFLContinuationAuthorityError, ValueError, OSError) as exc:
        raise Proto19Event1RunnerError(
            "event-one continuation authority differs"
        ) from exc
    if receipt.progression_plan.branch != "GR-0":
        raise Proto19Event1RunnerError(
            "event-one continuation opened an unauthorized branch"
        )
    return receipt


def authorize_sid1_event_authority(
    root: Path,
    *,
    continuation_commit: str,
    manifest_path: Path,
    store_anchor: Mapping[str, Any],
) -> sid1.SID1RecoveryReceipt:
    """Bind the corrected image to only the frozen generation-seven edge."""
    repository = Path(root).resolve()
    try:
        relative = manifest_path.resolve().relative_to(repository).as_posix()
    except ValueError as exc:
        raise Proto19Event1RunnerError(
            "launch manifest is outside the repository"
        ) from exc
    try:
        receipt = sid1.authorize_sid1_recovery(
            repository,
            continuation_commit=continuation_commit,
            launch_manifest_path=relative,
            store_anchor=store_anchor,
        )
    except (sid1.SID1AuthorityError, ValueError, OSError) as exc:
        raise Proto19Event1RunnerError("SID1 recovery authority differs") from exc
    if (
        receipt.progression_plan.branch != "GR-0"
        or receipt.progression_plan.amplitude != "3"
        or receipt.campaign_id != sid1.CAMPAIGN_ID
        or receipt.original_plan_sha256 != sid1.PLAN_SHA256
    ):
        raise Proto19Event1RunnerError("SID1 authority opened an unauthorized branch")
    return receipt


def _require_exact_sid1_successor(
    checkpoint: CampaignCheckpoint,
    receipt: sid1.SID1RecoveryReceipt,
) -> None:
    """Require every frozen metadata fact of the no-advance successor."""
    member = checkpoint.members.get(sid1.MEMBER_KEY)
    if member is None:
        raise Proto19Event1RunnerError("SID1 successor omits its owned member")
    cursor = member.cursor
    if (
        checkpoint.generation != 8
        or checkpoint.sha256 != receipt.expected_checkpoint_sha256
        or checkpoint.parent_sha256 != receipt.recovery_checkpoint_sha256
        or checkpoint.journal_sequence != 8
        or checkpoint.journal_tip_sha256 != receipt.expected_transition_sha256
        or checkpoint.disposition != "nonterminal"
        or checkpoint.terminal is not None
        or checkpoint.event != 23
        or checkpoint.target.get("rational") != "3/2"
        or member.descriptor_sha256 != receipt.descriptor_sha256
        or cursor.get("cursor_chain_sha256") != receipt.expected_cursor_sha256
        or cursor.get("accepted_state_sha256") != receipt.descriptor_sha256
        or cursor.get("accepted_boundary_time", {}).get("binary64_hex")
        != sid1.MEMBER_TIME_HEX
        or cursor.get("mode") != "RETRY_PENDING"
        or member.cfl_current != 1
        or member.cfl_total != 1
        or member.ledger.get("cumulative_temporal_retry_count") != 1
        or member.pending_owner != "temporal"
        or member.pending_cap_hex != receipt.expected_pending_cap_hex
    ):
        raise Proto19Event1RunnerError("SID1 successor differs from its frozen edge")


def _require_exact_sid1_preview(
    preview: recovery.TDG6RecoveryPreview,
    receipt: sid1.SID1RecoveryReceipt,
) -> None:
    """Cross-check the pure recovery derivation against the prospective freeze."""
    if (
        preview.predecessor_checkpoint_sha256 != receipt.recovery_checkpoint_sha256
        or preview.evolution_state_sha256 != receipt.evolution_state_sha256
        or preview.physical_state_advanced
        or len(preview.records) != 2
        or preview.records[0].get("kind") != "tdg6_rejection"
        or preview.records[0].get("record_sha256") != receipt.suffix_sha256
        or preview.records[1].get("kind") != "cursor_transition"
        or preview.records[1].get("record_sha256") != receipt.expected_transition_sha256
        or preview.records[1].get("payload", {}).get("successor_cursor_sha256")
        != receipt.expected_cursor_sha256
    ):
        raise Proto19Event1RunnerError("SID1 recovery preview differs")
    _require_exact_sid1_successor(preview.checkpoint, receipt)


def _sid1_preflight_context(
    root: Path,
    *,
    continuation_commit: str,
    manifest_path: Path,
    store_root: Path,
    permitted_writer_states: frozenset[str],
) -> tuple[
    HLT16CampaignStore,
    sid1.SID1RecoveryReceipt,
    recovery.TDG6RecoveryPreview,
]:
    """Read and cross-bind the current SID1 edge without publishing it."""
    repository = Path(root).resolve()
    store = HLT16CampaignStore(store_root)
    try:
        anchor = sid1.derive_live_anchor(repository)
        receipt = authorize_sid1_event_authority(
            repository,
            continuation_commit=continuation_commit,
            manifest_path=manifest_path,
            store_anchor=anchor,
        )
        status = store.inspect_recovery()
        preview = recovery.preview_tdg6_recovery(store)
    except (
        sid1.SID1AuthorityError,
        recovery.HLT16RecoveryError,
        HLT16CampaignStoreError,
        KeyError,
        TypeError,
        ValueError,
    ) as exc:
        raise Proto19Event1RunnerError("SID1 preflight evidence differs") from exc
    if (
        not status.active_write
        or status.writer_state not in permitted_writer_states
        or status.state != "recognized_uncheckpointed_journal_suffix"
    ):
        raise Proto19Event1RunnerError("SID1 writer/recovery state differs")
    _require_exact_sid1_preview(preview, receipt)
    return store, receipt, preview


def _require_exact_sid1_post_recovery(
    store: HLT16CampaignStore,
    receipt: sid1.SID1RecoveryReceipt,
    *,
    expected_owner_token: str | None,
) -> CampaignCheckpoint:
    """Independently authenticate the exact clean generation-eight boundary."""
    try:
        snapshot = store.authenticated_snapshot()
        status = store.inspect_recovery()
    except HLT16CampaignStoreError as exc:
        raise Proto19Event1RunnerError(
            "SID1 recovery checkpoint does not validate"
        ) from exc
    if (
        snapshot.suffix_classification != "clean"
        or snapshot.suffix_records
        or snapshot.staging_paths
        or snapshot.orphan_descriptor_sha256 is not None
        or snapshot.orphan_payload_semantic_sha256 is not None
        or snapshot.terminal_lock_present
    ):
        raise Proto19Event1RunnerError("SID1 post-recovery store is not clean")
    if expected_owner_token is None:
        if status.active_write or status.writer_state is not None:
            raise Proto19Event1RunnerError("SID1 post-recovery writer did not retire")
    else:
        try:
            writer = store._active_writer()
        except HLT16CampaignStoreError as exc:
            raise Proto19Event1RunnerError("SID1 recovery writer differs") from exc
        if (
            not status.active_write
            or status.writer_state != "live_local_writer"
            or writer is None
            or writer.get("owner_token") != expected_owner_token
            or writer.get("host") != os.uname().nodename
            or writer.get("pid") != os.getpid()
        ):
            raise Proto19Event1RunnerError("SID1 recovery writer ownership differs")
    _require_exact_sid1_successor(snapshot.checkpoint, receipt)
    member = snapshot.checkpoint.members[sid1.MEMBER_KEY]
    try:
        persisted = store.load_state(member.descriptor_sha256)
        physical = array_content_sha256(
            persisted.arrays["u"], persisted.arrays["p"], persisted.arrays["q"]
        )
    except (HLT16CampaignStoreError, KeyError, ValueError) as exc:
        raise Proto19Event1RunnerError("SID1 persisted physical state differs") from exc
    if physical != receipt.evolution_state_sha256:
        raise Proto19Event1RunnerError("SID1 persisted physical identity advanced")
    return snapshot.checkpoint


def _require_exact_continuation_checkpoint(
    checkpoint: CampaignCheckpoint,
    receipt: CFLContinuationReceipt,
) -> None:
    """Reject any store movement between anchor capture and continuation use."""
    member = checkpoint.members.get(receipt.recovery_member_key)
    if (
        checkpoint.generation != receipt.recovery_checkpoint_generation
        or checkpoint.sha256 != receipt.recovery_checkpoint_sha256
        or checkpoint.journal_sequence != receipt.recovery_journal_sequence
        or checkpoint.journal_tip_sha256 != receipt.recovery_journal_tip_sha256
        or member is None
        or member.descriptor_sha256 != receipt.recovery_member_descriptor_sha256
    ):
        raise Proto19Event1RunnerError(
            "continuation checkpoint moved after authority capture"
        )


def authorize_and_reconstruct(
    root: Path,
    *,
    authorization_commit: str,
    manifest_path: Path,
    reconstruct: Callable[[Path], ReconstructedGR0MemberSet] = reconstruct_gr0_members,
) -> tuple[LaunchAuthorityReceipt, ReconstructedGR0MemberSet]:
    """Authenticate and perform the one permitted GEN0 raw import."""
    receipt = authorize_event_authority(
        root,
        authorization_commit=authorization_commit,
        manifest_path=manifest_path,
    )
    try:
        members = reconstruct(Path(root).resolve())
    except (ValueError, OSError) as exc:
        raise Proto19Event1RunnerError("GEN0 reconstruction failed") from exc
    if (
        members.branch != "GR-0"
        or members.amplitude != "3"
        or tuple(members.members) != MEMBER_KEYS
        or receipt.progression_plan.branch != "GR-0"
    ):
        raise Proto19Event1RunnerError(
            "event-one reconstruction opened an unauthorized branch"
        )
    return receipt, members


def _acquire_writer(store: HLT16CampaignStore, checkpoint: CampaignCheckpoint) -> str:
    """Acquire a clean lease or take over only a same-host verified-dead lease."""
    host = os.uname().nodename
    try:
        store.recover_interrupted_writer_acquisition(checkpoint, host=host)
        status = store.recover()
        if not status.active_write:
            if status.state != "clean_checkpoint":
                raise Proto19Event1RunnerError(
                    f"campaign requires recovery before writer acquisition: {status.state}"
                )
            return store.acquire_active_writer(checkpoint, host=host)
        if status.writer_state != "stale_verified_lock":
            raise Proto19Event1RunnerError(
                f"campaign writer is not safely recoverable: {status.writer_state or status.state}"
            )
        return store.take_over_stale_writer(checkpoint, host=host)
    except HLT16CampaignStoreError as exc:
        raise Proto19Event1RunnerError("campaign writer authority differs") from exc


def _require_exact_gen0_store(store: HLT16CampaignStore) -> None:
    """Refuse every non-GEN0 suffix before the sole raw-source reimport.

    ``runtime_checkpoint_present`` deliberately detects a published HLT16
    checkpoint, but a crash can leave an earlier payload, descriptor, journal,
    lock, or staging leaf.  Those cases must be invalid rather than an excuse
    to reopen the historical source arrays.  This is read-only and mirrors the
    status command's exact eight-leaf boundary without importing it.
    """
    try:
        source = store._legacy_inventory(repair_stages=False)
        leaves, _directories = _walk_tree(store.root)
    except (HLT16CampaignStoreError, OSError, ValueError) as exc:
        raise Proto19Event1RunnerError(
            "store is not an exact authenticated GEN0 boundary"
        ) from exc
    expected = {
        "receipts/generation-zero.json",
        f"checkpoints/00000000000000000000-{SEALED_GEN0_CHECKPOINT_SHA256}.json",
        *(f"states/{address}.json" for address in source.values()),
    }
    if set(leaves) != expected:
        raise Proto19Event1RunnerError(
            "store has a post-GEN0 runtime, staging, lock, or foreign suffix"
        )


_SHA256_RE = re.compile(r"[0-9a-f]{64}")
_JOURNAL_RE = re.compile(r"00000000000000000000-[0-9a-f]{64}\.journal")


def _recoverable_bootstrap_prefix(store: HLT16CampaignStore) -> bool:
    """Recognize the only pre-generation-one replay prefix.

    The bridge itself publishes six immutable GEN0-derived payload/descriptor
    pairs, one genesis journal record, and then its first checkpoint.  Before
    that checkpoint exists, replay is lawful only when the tree contains that
    bounded task-owned prefix.  A one-link stage is rejected here before the
    permitted GEN0 import because its creator cannot be proved.  Content
    equality for finals and the sole replayable two-link stage is completed by
    ``initialize_or_recover_generation_one`` after that import: its
    content-addressed no-replace publications reject any byte that is not the
    exact reconstructed bridge object.

    Returns ``False`` only when a complete non-GEN0 checkpoint exists, which
    sends the caller to persisted-only recovery.  Every other malformed or
    foreign pre-generation tree raises.
    """
    try:
        source = store._legacy_inventory(repair_stages=False)
        leaves, directories = _walk_tree(store.root)
    except (HLT16CampaignStoreError, OSError, ValueError) as exc:
        raise Proto19Event1RunnerError(
            "store does not retain authenticated GEN0 authority"
        ) from exc
    base = {
        "receipts/generation-zero.json",
        f"checkpoints/00000000000000000000-{SEALED_GEN0_CHECKPOINT_SHA256}.json",
        *(f"states/{address}.json" for address in source.values()),
    }
    known_directories = {
        "checkpoints",
        "receipts",
        "states",
        "payloads",
        "journal",
        "locks",
    }
    if not {"checkpoints", "receipts", "states"}.issubset(directories) or not set(
        directories
    ).issubset(known_directories):
        raise Proto19Event1RunnerError("bootstrap prefix directory inventory differs")

    extras = set(leaves) - base
    if any(
        item.startswith("checkpoints/") and not Path(item).name.startswith(".")
        for item in extras
    ):
        # A complete first checkpoint is the irreversible handoff to
        # persisted-only recovery.  Its full verifier owns all later suffixes.
        return False

    counts = {"payloads": 0, "states": 0, "journal": 0, "stages": 0}

    def permitted_final(directory: str, name: str) -> bool:
        if directory == "payloads":
            return bool(re.fullmatch(r"[0-9a-f]{64}\.npz", name))
        if directory == "states":
            return bool(re.fullmatch(r"[0-9a-f]{64}\.json", name))
        if directory == "journal":
            return bool(_JOURNAL_RE.fullmatch(name))
        if directory == "checkpoints":
            return bool(re.fullmatch(r"00000000000000000001-[0-9a-f]{64}\.json", name))
        return False

    for relative in extras:
        path = Path(relative)
        if len(path.parts) != 2:
            raise Proto19Event1RunnerError("bootstrap prefix leaf differs")
        directory, name = path.parts
        if directory == "locks":
            if name != "bootstrap.guard":
                raise Proto19Event1RunnerError(
                    "bootstrap prefix has writer, terminal, or retired lease"
                )
            continue
        staged = _stage_target(name)
        if staged is not None:
            final, _digest = staged
            if directory not in {
                "payloads",
                "states",
                "journal",
                "checkpoints",
            } or not permitted_final(directory, final):
                raise Proto19Event1RunnerError("bootstrap publication stage differs")
            try:
                validated = _validated_stage(store.root, relative)
            except HLT16CampaignStoreError as exc:
                raise Proto19Event1RunnerError(
                    "bootstrap publication stage cannot be proved"
                ) from exc
            if validated.link_count != 2:
                # Do not reopen the historical source arrays for a tree whose
                # sole prepublication inode has no independently proved owner.
                # Only link(stage, final) creates the exact replay fact used by
                # the campaign store.
                raise Proto19Event1RunnerError(
                    "unlinked bootstrap stage ownership is unprovable"
                )
            counts["stages"] += 1
            continue
        if not permitted_final(directory, name):
            raise Proto19Event1RunnerError(
                "bootstrap prefix has a foreign runtime leaf"
            )
        counts[directory] += 1

    if (
        counts["payloads"] > 6
        or counts["states"] > 6
        or counts["journal"] > 1
        or counts["stages"] > 1
    ):
        raise Proto19Event1RunnerError("bootstrap prefix publication counts differ")
    return True


def _release_writer_if_closed(
    store: HLT16CampaignStore,
    *,
    owner_token: str,
) -> bool:
    """Release only beside a complete checkpoint with no uncommitted suffix."""
    try:
        snapshot = store.authenticated_snapshot()
        if (
            snapshot.suffix_classification != "clean"
            or snapshot.staging_paths
            or snapshot.orphan_descriptor_sha256 is not None
            or snapshot.orphan_payload_semantic_sha256 is not None
        ):
            return False
        store.release_active_writer(
            owner_token=owner_token,
            host=os.uname().nodename,
            pid=os.getpid(),
        )
        return True
    except HLT16CampaignStoreError:
        return False


def _reconcile_before_attempt(
    store: HLT16CampaignStore,
    checkpoint: CampaignCheckpoint,
) -> tuple[CampaignCheckpoint, str | None]:
    """Close a proved suffix or return one payload identity that must be rerun."""
    snapshot = store.authenticated_snapshot()
    if snapshot.checkpoint.sha256 != checkpoint.sha256:
        checkpoint = snapshot.checkpoint
    if snapshot.staging_paths:
        # ``recover`` may remove only a proven two-link after-link stage.  It
        # never removes a one-link prepublication stage.  Reinspect after that
        # bounded cleanup before deciding whether a numerical replay is still
        # necessary.
        store.recover()
        snapshot = store.authenticated_snapshot()
        checkpoint = snapshot.checkpoint
        if not snapshot.staging_paths:
            if snapshot.suffix_classification == "payload_orphan":
                return checkpoint, snapshot.orphan_payload_semantic_sha256
            decision = recovery.reconcile(store)
            if decision.disposition == "invalid_ambiguous_suffix":
                raise Proto19Event1RunnerError(
                    "post-link publication suffix is authenticated but ambiguous"
                )
            if decision.disposition == "rerun_from_checkpoint_required":
                replay = store.authenticated_snapshot()
                return checkpoint, replay.orphan_payload_semantic_sha256
            if decision.disposition == "publication_replay_required":
                raise Proto19Event1RunnerError(
                    "post-link publication stage remained after bounded cleanup"
                )
            return store.authenticated_snapshot().checkpoint, None
        # Repeating the same deterministic numerical/publication edge is the
        # only operation allowed to complete or replace a one-link stage.
        # A target payload hash, when visible, additionally constrains replay.
        semantic: str | None = None
        if len(snapshot.staging_paths) != 1:
            raise Proto19Event1RunnerError("more than one publication stage is visible")
        path = Path(snapshot.staging_paths[0])
        marker = ".hlt16-stage-"
        final = path.name[1:].split(marker, 1)[0]
        if path.parent.as_posix() == "payloads" and final.endswith(".npz"):
            semantic = final[:-4]
        return checkpoint, semantic
    if snapshot.suffix_classification == "payload_orphan":
        return checkpoint, snapshot.orphan_payload_semantic_sha256
    decision = recovery.reconcile(store)
    if decision.disposition == "invalid_ambiguous_suffix":
        raise Proto19Event1RunnerError("campaign suffix is authenticated but ambiguous")
    if decision.disposition in {
        "rerun_from_checkpoint_required",
        "publication_replay_required",
    }:
        raise Proto19Event1RunnerError(
            "campaign recovery did not identify one replay edge"
        )
    return store.authenticated_snapshot().checkpoint, None


def _advance_authenticated_event(
    plan: Any,
    store: HLT16CampaignStore,
    checkpoint: CampaignCheckpoint,
    templates: Mapping[str, object],
    *,
    attempt_runner: Callable[..., Any] | None = None,
) -> dict[str, Any]:
    """Advance one already authenticated persisted event in canonical order."""
    token = _acquire_writer(store, checkpoint)
    capability = store.active_writer_capability(checkpoint)
    clean_release = False
    try:
        while True:
            checkpoint, expected_orphan = _reconcile_before_attempt(store, checkpoint)
            if checkpoint.disposition in {"scientific_terminal", "invalid_terminal"}:
                if not store.authenticated_snapshot().terminal_lock_present:
                    store.publish_terminal_lock(checkpoint, capability=capability)
                clean_release = True
                return _summary(checkpoint, state=checkpoint.disposition)
            if checkpoint.event == 24:
                if checkpoint.disposition != "event_complete":
                    raise Proto19Event1RunnerError(
                        "event 24 lacks its exact event-complete checkpoint"
                    )
                clean_release = True
                return _summary(checkpoint, state="event_one_complete")
            key = _member_to_advance(checkpoint)
            if key is None:
                event = runtime.commit_common_event(checkpoint, store)
                checkpoint = event.checkpoint
                clean_release = True
                return _summary(
                    checkpoint,
                    state="event_one_complete",
                    record_sha256s=event.record_sha256s,
                )
            member = copy.deepcopy(templates[key])
            kwargs: dict[str, Any] = {}
            if attempt_runner is not None:
                kwargs["runner"] = attempt_runner
            attempt = runtime.execute_one_attempt(
                plan,
                store,
                checkpoint,
                member,
                key=key,
                expected_orphan_semantic_sha256=expected_orphan,
                **kwargs,
            )
            checkpoint = attempt.checkpoint
            if checkpoint.disposition in {"scientific_terminal", "invalid_terminal"}:
                clean_release = True
                return _summary(
                    checkpoint,
                    state=checkpoint.disposition,
                    record_sha256s=attempt.record_sha256s,
                )
    finally:
        if clean_release:
            if not _release_writer_if_closed(store, owner_token=token):
                raise Proto19Event1RunnerError(
                    "completed checkpoint exists but its writer lease did not close cleanly"
                )
        # On an interruption or ambiguous exception the lease deliberately
        # remains.  When this PID exits it becomes a verified stale recovery
        # marker bound to the authenticated plan/checkpoint lineage.


def run_first_event(
    root: Path,
    *,
    authorization_commit: str,
    manifest_path: Path,
    store_root: Path,
    reconstruct: Callable[[Path], ReconstructedGR0MemberSet] = reconstruct_gr0_members,
    static_shells: Callable[[Path], Mapping[str, object]] = build_static_gr0_shells,
    attempt_runner: Callable[..., Any] | None = None,
) -> dict[str, Any]:
    """Advance exactly event 23 to event 24, recoverably and in fixed order."""
    store = HLT16CampaignStore(store_root)
    # Once HLT16 has any post-GEN0 checkpoint leaf, even a malformed one, the
    # historical raw source import is permanently unavailable to this runner.
    # A positive inventory fact is authenticated immediately; it is never used
    # to pick a more convenient recovery history.
    if not _recoverable_bootstrap_prefix(store):
        receipt = authorize_event_authority(
            root,
            authorization_commit=authorization_commit,
            manifest_path=manifest_path,
        )
        plan = receipt.progression_plan
        checkpoint = runtime.recover_persisted_generation_one(plan, store)
        try:
            templates = static_shells(Path(root).resolve())
        except (ValueError, OSError) as exc:
            raise Proto19Event1RunnerError(
                "persisted GR-0 shell factory failed"
            ) from exc
        if tuple(templates) != MEMBER_KEYS:
            raise Proto19Event1RunnerError("persisted GR-0 shell order differs")
    else:
        receipt, reconstructed = authorize_and_reconstruct(
            root,
            authorization_commit=authorization_commit,
            manifest_path=manifest_path,
            reconstruct=reconstruct,
        )
        plan = receipt.progression_plan
        checkpoint = runtime.initialize_or_recover_generation_one(
            plan, store, reconstructed
        )
        templates = reconstructed.members
    return _advance_authenticated_event(
        plan,
        store,
        checkpoint,
        templates,
        attempt_runner=attempt_runner,
    )


def continuation_preflight(
    root: Path,
    *,
    continuation_commit: str,
    manifest_path: Path,
    store_root: Path,
) -> dict[str, Any]:
    """Authenticate the corrected image and exact checkpoint without mutation."""
    store = HLT16CampaignStore(store_root)
    if _recoverable_bootstrap_prefix(store):
        raise Proto19Event1RunnerError(
            "continuation requires a persisted post-generation-zero checkpoint"
        )
    receipt = authorize_continuation_event_authority(
        root,
        continuation_commit=continuation_commit,
        manifest_path=manifest_path,
        store=store,
    )
    checkpoint = runtime.recover_persisted_generation_one(
        receipt.progression_plan, store
    )
    _require_exact_continuation_checkpoint(checkpoint, receipt)
    return {
        "schema": "FGC-1-PRO19-event1-continuation-preflight-v1",
        "runner_id": RUNNER_ID,
        "continuation_commit": receipt.continuation_commit,
        "original_authorization_commit": receipt.original_authorization_commit,
        "original_plan_sha256": receipt.original_plan_sha256,
        "campaign_id": receipt.campaign_id,
        "checkpoint_generation": checkpoint.generation,
        "checkpoint_sha256": checkpoint.sha256,
        "journal_sequence": checkpoint.journal_sequence,
        "journal_tip_sha256": checkpoint.journal_tip_sha256,
        "active_target": dict(checkpoint.target),
        "safe_to_resume": True,
        "output_created": False,
        "state_advanced": False,
        "candidate_branch_opened": False,
        "physical_result_earned": False,
    }


def resume_corrected_first_event(
    root: Path,
    *,
    continuation_commit: str,
    manifest_path: Path,
    store_root: Path,
    static_shells: Callable[[Path], Mapping[str, object]] = build_static_gr0_shells,
    attempt_runner: Callable[..., Any] | None = None,
) -> dict[str, Any]:
    """Resume only the exact generation-six campaign under corrected code."""
    store = HLT16CampaignStore(store_root)
    if _recoverable_bootstrap_prefix(store):
        raise Proto19Event1RunnerError(
            "continuation cannot reopen generation-zero raw inputs"
        )
    receipt = authorize_continuation_event_authority(
        root,
        continuation_commit=continuation_commit,
        manifest_path=manifest_path,
        store=store,
    )
    plan = receipt.progression_plan
    checkpoint = runtime.recover_persisted_generation_one(plan, store)
    _require_exact_continuation_checkpoint(checkpoint, receipt)
    try:
        templates = static_shells(Path(root).resolve())
    except (ValueError, OSError) as exc:
        raise Proto19Event1RunnerError("persisted GR-0 shell factory failed") from exc
    if tuple(templates) != MEMBER_KEYS:
        raise Proto19Event1RunnerError("persisted GR-0 shell order differs")
    result = _advance_authenticated_event(
        plan,
        store,
        checkpoint,
        templates,
        attempt_runner=attempt_runner,
    )
    return {
        **result,
        "continuation_commit": receipt.continuation_commit,
        "continuation_config_sha256": receipt.continuation_config_sha256,
        "continuation_result_sha256": receipt.continuation_result_sha256,
        "resumed_from_checkpoint_sha256": receipt.recovery_checkpoint_sha256,
        "earlier_stop_physical_inference": False,
    }


def sid1_preflight(
    root: Path,
    *,
    continuation_commit: str,
    manifest_path: Path,
    store_root: Path,
) -> dict[str, Any]:
    """Prove the exact recovery successor without a writer or publication."""
    _store, receipt, preview = _sid1_preflight_context(
        root,
        continuation_commit=continuation_commit,
        manifest_path=manifest_path,
        store_root=store_root,
        permitted_writer_states=frozenset({"stale_verified_lock"}),
    )
    return {
        "schema": "FGC-1-PRO19-SID1-preflight-v1",
        "runner_id": RUNNER_ID,
        "continuation_commit": receipt.continuation_commit,
        "manifest_sha256": receipt.manifest_sha256,
        "campaign_id": receipt.campaign_id,
        "predecessor_checkpoint_generation": receipt.recovery_checkpoint_generation,
        "predecessor_checkpoint_sha256": receipt.recovery_checkpoint_sha256,
        "suffix_sequence": receipt.suffix_sequence,
        "suffix_sha256": receipt.suffix_sha256,
        "expected_transition_sha256": receipt.expected_transition_sha256,
        "expected_checkpoint_generation": preview.checkpoint.generation,
        "expected_checkpoint_sha256": receipt.expected_checkpoint_sha256,
        "expected_cursor_sha256": receipt.expected_cursor_sha256,
        "expected_pending_cap_hex": receipt.expected_pending_cap_hex,
        "safe_to_recover": True,
        "safe_to_resume_trajectory": False,
        "output_created": False,
        "state_advanced": False,
        "physical_state_advanced": False,
        "rejected_proposal_replayed": False,
        "candidate_branch_opened": False,
        "physical_result_earned": False,
    }


def sid1_recover(
    root: Path,
    *,
    continuation_commit: str,
    manifest_path: Path,
    store_root: Path,
) -> dict[str, Any]:
    """Publish only SID1's exact no-advance cursor/ledger successor."""
    store, receipt, _preview = _sid1_preflight_context(
        root,
        continuation_commit=continuation_commit,
        manifest_path=manifest_path,
        store_root=store_root,
        permitted_writer_states=frozenset({"stale_verified_lock"}),
    )
    predecessor = store.authenticated_snapshot().checkpoint
    token = _acquire_writer(store, predecessor)
    clean_release = False
    try:
        # The lease takeover changes only lock ownership.  Re-read all frozen
        # checkpoint/journal/payload facts before the single publication.
        post_store, post_receipt, post_preview = _sid1_preflight_context(
            root,
            continuation_commit=continuation_commit,
            manifest_path=manifest_path,
            store_root=store_root,
            permitted_writer_states=frozenset({"live_local_writer"}),
        )
        if post_store.root.resolve() != store.root.resolve():
            raise Proto19Event1RunnerError("SID1 store identity changed after lease")
        if post_receipt != receipt:
            raise Proto19Event1RunnerError("SID1 authority changed after lease")
        _require_exact_sid1_preview(post_preview, receipt)
        try:
            decision = recovery.reconcile(store)
        except recovery.HLT16RecoveryError as exc:
            raise Proto19Event1RunnerError("SID1 recovery publication failed") from exc
        if (
            decision.disposition != "reconciled_checkpoint"
            or not decision.published
            or decision.checkpoint_sha256 != receipt.expected_checkpoint_sha256
        ):
            raise Proto19Event1RunnerError("SID1 recovery published a different edge")
        checkpoint = _require_exact_sid1_post_recovery(
            store, receipt, expected_owner_token=token
        )
        clean_release = True
    finally:
        if clean_release:
            if not _release_writer_if_closed(store, owner_token=token):
                raise Proto19Event1RunnerError(
                    "SID1 checkpoint exists but its writer lease did not close cleanly"
                )
    # Re-open the store read-only after lease retirement so the returned claim
    # is about a restartable boundary, not merely the writer's in-memory view.
    checkpoint = _require_exact_sid1_post_recovery(
        HLT16CampaignStore(store_root), receipt, expected_owner_token=None
    )
    return {
        "schema": "FGC-1-PRO19-SID1-recovery-v1",
        "runner_id": RUNNER_ID,
        "continuation_commit": receipt.continuation_commit,
        "campaign_id": receipt.campaign_id,
        "checkpoint_generation": checkpoint.generation,
        "checkpoint_sha256": checkpoint.sha256,
        "journal_sequence": checkpoint.journal_sequence,
        "journal_tip_sha256": checkpoint.journal_tip_sha256,
        "member_cursor_sha256": checkpoint.members[sid1.MEMBER_KEY].cursor[
            "cursor_chain_sha256"
        ],
        "recovery_checkpoint_published": True,
        "safe_to_restart": True,
        "safe_to_resume_trajectory": False,
        "physical_state_advanced": False,
        "rejected_proposal_replayed": False,
        "candidate_branch_opened": False,
        "calibration_result_earned": False,
        "physical_result_earned": False,
    }


def sid1_resume(
    root: Path,
    *,
    continuation_commit: str,
    manifest_path: Path,
    store_root: Path,
    static_shells: Callable[[Path], Mapping[str, object]] = build_static_gr0_shells,
    attempt_runner: Callable[..., Any] | None = None,
) -> dict[str, Any]:
    """Remain closed until a separate resume authority and preflight exist."""
    del (
        root,
        continuation_commit,
        manifest_path,
        store_root,
        static_shells,
        attempt_runner,
    )
    raise Proto19Event1RunnerError(
        "SID1 numerical resume remains closed after SID2; a separately committed "
        "resume authority and read-only preflight are required"
    )


def _validated_sid3_authority(authority: Any) -> Any:
    """Require the exact immutable receipt produced inside the import guard."""
    from recursive_horizons.fgc.evolution import proto19_resume_authority as sid3

    if not isinstance(authority, sid3.ResumeAuthorityReceipt):
        raise Proto19Event1RunnerError("SID3 resume authority receipt type differs")
    plan = authority.progression_plan
    if (
        authority.store_path != sid3.STORE_PATH
        or authority.checkpoint_generation != sid3.GENERATION
        or authority.checkpoint_sha256 != sid3.GENERATION8_CHECKPOINT_SHA256
        or authority.journal_sequence != sid3.GENERATION8_JOURNAL_SEQUENCE
        or authority.journal_tip_sha256 != sid3.GENERATION8_JOURNAL_SHA256
        or authority.original_authorization_commit != sid3.ORIGINAL_AUTHORIZATION_COMMIT
        or authority.original_plan_sha256 != sid3.ORIGINAL_PLAN_SHA256
        or authority.campaign_id != sid3.CAMPAIGN_ID
        or plan.authorization_commit != sid3.ORIGINAL_AUTHORIZATION_COMMIT
        or plan.sha256 != sid3.ORIGINAL_PLAN_SHA256
        or plan.campaign_id != sid3.CAMPAIGN_ID
        or plan.branch != "GR-0"
        or plan.amplitude != "3"
        or plan.common_event_index != 23
        or plan.target_time.rational != "3/2"
        or set(authority.static_input_bytes) != set(sid3.MANDATORY_STATIC_INPUT_PATHS)
    ):
        raise Proto19Event1RunnerError("SID3 resume authority scope differs")
    return authority


def _require_sid3_gr0_templates(templates: Mapping[str, object]) -> None:
    """Reject any shell that is not the exact GR-0 runtime type and operator."""
    if tuple(templates) != MEMBER_KEYS:
        raise Proto19Event1RunnerError("SID3 GR-0 template order differs")
    for key in MEMBER_KEYS:
        member = templates[key]
        if (
            not isinstance(member, Proto14RunMember)
            or member.key != key
            or member.amplitude != "3"
            or not isinstance(member.operator, Proto12GR0EvolutionOperator)
        ):
            raise Proto19Event1RunnerError(
                f"SID3 candidate or template identity differs: {key}"
            )


def _require_sid3_generation8_ancestor(
    store: HLT16CampaignStore,
    authority: Any,
) -> tuple[CampaignCheckpoint, Any, Any]:
    """Bind the current append-only lineage to SID2's exact generation eight."""
    try:
        snapshot = store.authenticated_snapshot()
        status = store.inspect_recovery()
        leaves, _directories = _walk_tree(store.root)
        checkpoint = runtime.recover_persisted_generation_one(
            authority.progression_plan, store
        )
    except Exception as exc:
        raise Proto19Event1RunnerError(
            "SID3 persisted generation-eight lineage does not validate"
        ) from exc
    expected_checkpoint = (
        f"checkpoints/{authority.checkpoint_generation:020d}-"
        f"{authority.checkpoint_sha256}.json"
    )
    expected_journal = (
        f"journal/{authority.journal_sequence:020d}-"
        f"{authority.journal_tip_sha256}.journal"
    )
    if (
        snapshot.checkpoint.sha256 != checkpoint.sha256
        or expected_checkpoint not in leaves
        or expected_journal not in leaves
        or checkpoint.generation < authority.checkpoint_generation
        or checkpoint.journal_sequence < authority.journal_sequence
        or checkpoint.authorization_commit != authority.original_authorization_commit
        or checkpoint.plan_sha256 != authority.original_plan_sha256
        or checkpoint.campaign_id != authority.campaign_id
        or tuple(checkpoint.members) != MEMBER_KEYS
    ):
        raise Proto19Event1RunnerError("SID3 generation-eight ancestor differs")
    return checkpoint, snapshot, status


def _sid3_resume_preflight_context(
    root: Path,
    *,
    authority: Any,
    store_root: Path,
    permit_recoverable_suffix: bool,
) -> SID3ResumePreflightContext:
    """Restore every current member in memory without acquiring a writer."""
    authority = _validated_sid3_authority(authority)
    repository = Path(root).resolve()
    expected_store = (repository / authority.store_path).resolve()
    if Path(store_root).resolve() != expected_store:
        raise Proto19Event1RunnerError("SID3 calibration store path differs")
    store = HLT16CampaignStore(expected_store)
    checkpoint, snapshot, status = _require_sid3_generation8_ancestor(store, authority)
    if (
        checkpoint.disposition != "nonterminal"
        or checkpoint.terminal is not None
        or checkpoint.event != 23
        or checkpoint.target.get("rational") != "3/2"
        or checkpoint.target.get("binary64_hex") != (3 / 2).hex()
        or snapshot.terminal_lock_present
        or status.terminal
    ):
        raise Proto19Event1RunnerError(
            "SID3 event is terminal, complete, or retargeted"
        )

    if permit_recoverable_suffix:
        if status.active_write and status.writer_state != "stale_verified_lock":
            raise Proto19Event1RunnerError("SID3 writer is live or not recoverable")
    elif (
        status.active_write
        or status.writer_state is not None
        or snapshot.suffix_classification != "clean"
        or snapshot.suffix_records
        or snapshot.staging_paths
        or snapshot.orphan_descriptor_sha256 is not None
        or snapshot.orphan_payload_semantic_sha256 is not None
    ):
        raise Proto19Event1RunnerError(
            "SID3 read-only preflight requires a clean writer-free checkpoint"
        )

    if checkpoint.generation == authority.checkpoint_generation:
        from recursive_horizons.fgc.evolution import proto19_resume_authority as sid3

        recovered = checkpoint.members[authority.recovery_member_key]
        if (
            recovered.cursor.get("cursor_chain_sha256")
            != authority.recovery_cursor_sha256
            or recovered.cursor.get("mode") != "RETRY_PENDING"
            or recovered.cursor.get("accepted_boundary_time", {}).get("binary64_hex")
            != sid3.GENERATION8_ACCEPTED_TIME_HEX
            or recovered.pending_owner != "temporal"
            or recovered.pending_cap_hex != authority.recovery_pending_cap_hex
            or recovered.cfl_current != 1
            or recovered.cfl_total != 1
            or recovered.ledger.get("cumulative_temporal_retry_count") != 1
        ):
            raise Proto19Event1RunnerError("SID3 generation-eight retry state differs")

    try:
        factory_inputs = {
            path: authority.static_input_bytes[path] for path in STATIC_INPUT_PATHS
        }
        templates = build_static_gr0_shells(static_input_bytes=factory_inputs)
    except (KeyError, TypeError, ValueError, OSError) as exc:
        raise Proto19Event1RunnerError(
            "SID3 captured GR-0 template construction failed"
        ) from exc
    _require_sid3_gr0_templates(templates)

    # Decode and overlay every selected payload now.  This proves the complete
    # six-member restart in memory while preserving the fresh templates that
    # the immediately following numerical operation will consume.
    try:
        for key in MEMBER_KEYS:
            restored = runtime.restore_member_with_overlay(
                store, checkpoint, copy.deepcopy(templates[key]), key=key
            )
            expected_time = checkpoint.members[key].cursor["accepted_boundary_time"][
                "binary64_hex"
            ]
            if restored.time.hex() != expected_time:
                raise Proto19Event1RunnerError(
                    f"SID3 restored member time differs: {key}"
                )
    except Proto19Event1RunnerError:
        raise
    except Exception as exc:
        raise Proto19Event1RunnerError(
            "SID3 persisted six-member restore failed"
        ) from exc
    _member_to_advance(checkpoint)
    return SID3ResumePreflightContext(
        authority=authority,
        store=store,
        checkpoint=checkpoint,
        templates=templates,
        suffix_classification=snapshot.suffix_classification,
        recovery_state=status.state,
    )


def sid3_resume_preflight(
    root: Path,
    *,
    authority: Any,
    store_root: Path,
) -> dict[str, Any]:
    """Authenticate the live generation-eight boundary without any mutation."""
    context = _sid3_resume_preflight_context(
        root,
        authority=authority,
        store_root=store_root,
        permit_recoverable_suffix=False,
    )
    checkpoint = context.checkpoint
    return {
        "schema": "FGC-1-PRO19-SID3-resume-preflight-v1",
        "runner_id": RUNNER_ID,
        "implementation_commit": authority.implementation_commit,
        "authority_commit": authority.authority_commit,
        "execution_closure_sha256": authority.closure_sha256,
        "authority_config_sha256": authority.config_sha256,
        "authority_result_sha256": authority.result_sha256,
        "campaign_id": checkpoint.campaign_id,
        "checkpoint_generation": checkpoint.generation,
        "checkpoint_sha256": checkpoint.sha256,
        "journal_sequence": checkpoint.journal_sequence,
        "journal_tip_sha256": checkpoint.journal_tip_sha256,
        "event": checkpoint.event,
        "target": dict(checkpoint.target),
        "suffix_classification": context.suffix_classification,
        "recovery_state": context.recovery_state,
        "member_times": {
            key: dict(checkpoint.members[key].cursor["accepted_boundary_time"])
            for key in MEMBER_KEYS
        },
        "member_modes": {
            key: checkpoint.members[key].cursor["mode"] for key in MEMBER_KEYS
        },
        "six_member_restore_passed": True,
        "writer_lease_acquired": False,
        "output_created": False,
        "state_advanced": False,
        "safe_to_resume_trajectory": True,
        "candidate_branch_opened": False,
        "shared_candidate_capable_definitions_loaded": True,
        "candidate_runtime_configuration_state_or_outcome_opened": False,
        "calibration_result_earned": False,
        "physical_result_earned": False,
    }


def sid3_resume(
    root: Path,
    *,
    authority: Any,
    store_root: Path,
) -> dict[str, Any]:
    """Recheck SID3 in process, then advance only the original GR-0 event."""
    context = _sid3_resume_preflight_context(
        root,
        authority=authority,
        store_root=store_root,
        permit_recoverable_suffix=True,
    )
    _require_sid3_gr0_templates(context.templates)
    result = _advance_authenticated_event(
        context.authority.progression_plan,
        context.store,
        context.checkpoint,
        context.templates,
    )
    return {
        **result,
        "implementation_commit": authority.implementation_commit,
        "authority_commit": authority.authority_commit,
        "execution_closure_sha256": authority.closure_sha256,
        "resumed_from_generation": context.checkpoint.generation,
        "resumed_from_checkpoint_sha256": context.checkpoint.sha256,
        "generation8_ancestor_sha256": authority.checkpoint_sha256,
        "candidate_branch_opened": False,
        "shared_candidate_capable_definitions_loaded": True,
        "candidate_runtime_configuration_state_or_outcome_opened": False,
        "calibration_result_earned": False,
        "physical_result_earned": False,
    }


def preflight(
    root: Path,
    *,
    authorization_commit: str,
    manifest_path: Path,
    store_root: Path,
) -> dict[str, Any]:
    # Preflight is allowed to make the single historical GEN0 import only
    # before *any* HLT16 runtime artifact exists.  It is not a shortcut around
    # persisted-payload restart after a bridge or interrupted publication.
    _require_exact_gen0_store(HLT16CampaignStore(store_root))
    receipt, reconstructed = authorize_and_reconstruct(
        root,
        authorization_commit=authorization_commit,
        manifest_path=manifest_path,
    )
    return {
        "schema": "FGC-1-PRO19-event1-preflight-v1",
        "runner_id": RUNNER_ID,
        "authorization_commit": receipt.authorization_commit,
        "manifest_sha256": receipt.manifest_sha256,
        "plan": receipt.progression_plan.mapping(),
        "plan_sha256": receipt.progression_plan.sha256,
        "member_keys": list(reconstructed.members),
        "branch": reconstructed.branch,
        "amplitude": reconstructed.amplitude,
        "output_created": False,
        "state_advanced": False,
        "candidate_branch_opened": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    parser.add_argument("--authorization-commit", required=False)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--store-root", type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--resume-preflight", action="store_true")
    mode.add_argument("--resume", action="store_true")
    mode.add_argument("--sid1-preflight", action="store_true")
    mode.add_argument("--sid1-recover", action="store_true")
    mode.add_argument("--sid1-resume", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = (args.manifest or (root / DEFAULT_MANIFEST)).resolve()
    store_root = (args.store_root or (root / DEFAULT_STORE)).resolve()
    try:
        canonical_root = ROOT.resolve()
        canonical_manifest = (canonical_root / DEFAULT_MANIFEST).resolve()
        canonical_store = (canonical_root / DEFAULT_STORE).resolve()
        if (
            root != canonical_root
            or manifest != canonical_manifest
            or store_root != canonical_store
        ):
            raise Proto19Event1RunnerError(
                "production event-one command requires the canonical repository, manifest, and calibration store"
            )
        commit = args.authorization_commit or _git_head(root)
        if args.preflight:
            result = preflight(
                root,
                authorization_commit=commit,
                manifest_path=manifest,
                store_root=store_root,
            )
        elif args.run:
            result = run_first_event(
                root,
                authorization_commit=commit,
                manifest_path=manifest,
                store_root=store_root,
            )
        elif args.resume_preflight:
            result = continuation_preflight(
                root,
                continuation_commit=commit,
                manifest_path=manifest,
                store_root=store_root,
            )
        elif args.resume:
            result = resume_corrected_first_event(
                root,
                continuation_commit=commit,
                manifest_path=manifest,
                store_root=store_root,
            )
        elif args.sid1_preflight:
            result = sid1_preflight(
                root,
                continuation_commit=commit,
                manifest_path=manifest,
                store_root=store_root,
            )
        elif args.sid1_recover:
            result = sid1_recover(
                root,
                continuation_commit=commit,
                manifest_path=manifest,
                store_root=store_root,
            )
        else:
            result = sid1_resume(
                root,
                continuation_commit=commit,
                manifest_path=manifest,
                store_root=store_root,
            )
    except (Proto19Event1RunnerError, HLT16CampaignStoreError, ValueError) as exc:
        print(
            _canonical(
                {
                    "schema": "FGC-1-PRO19-event1-run-status-v1",
                    "runner_id": RUNNER_ID,
                    "state": "invalid_implementation_or_nonconverged_run",
                    "detail": str(exc),
                    "candidate_branch_opened": False,
                    "physical_result_earned": False,
                }
            )
        )
        return 2
    print(_canonical(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
