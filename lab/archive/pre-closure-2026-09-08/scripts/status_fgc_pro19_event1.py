#!/usr/bin/env python3
"""Read-only status for the first authenticated PRO19 GR-0 event.

This command never imports an evolved payload, acquires a writer lease,
repairs a publication stage, or advances a trajectory.  It distinguishes the
sealed PROTO17 generation-zero predecessor from a later HLT16 checkpoint so a
human can tell whether it is safe to *start* or safe to *restart*.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Iterable, Mapping


ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = str(ROOT / "src")
if _SOURCE_ROOT not in sys.path:
    sys.path.insert(0, _SOURCE_ROOT)

from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
    HLT16CampaignStoreError,
    _walk_tree,
)
from recursive_horizons.fgc.evolution.hlt16_state_store import (  # noqa: E402
    SEALED_GEN0_CHECKPOINT_SHA256,
    read_nofollow,
)
from recursive_horizons.fgc.evolution.proto17_pure_construction import (  # noqa: E402
    MEMBER_KEYS,
    validate_checkpoint as validate_proto17_checkpoint,
)
from recursive_horizons.fgc.evolution.proto19_progression_freeze import (  # noqa: E402
    Proto19FreezeError,
    build_freeze,
)


FREEZE_RELATIVE = Path("configs/fgc/fgc-1-pro19-frz1.toml")
STORE_RELATIVE = Path("runs/fgc-2-sf1/proto17/calibration")
RUNNER_BASENAME = "run_fgc_pro19_event1.py"
SID3_PROCESS_MARKER = "FGC-PRO19-SID3-EVENT1"
RUNNER_MARKERS = (RUNNER_BASENAME, SID3_PROCESS_MARKER)
ROOT_JOURNAL_SHA256 = "0" * 64


class Proto19StatusError(ValueError):
    """The static authority or live store cannot support a status claim."""


def _canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)


def _freeze(root: Path) -> Mapping[str, Any]:
    try:
        raw = (root / FREEZE_RELATIVE).read_text("utf-8")
        return build_freeze(root, tomllib.loads(raw))
    except (OSError, tomllib.TOMLDecodeError, Proto19FreezeError) as exc:
        raise Proto19StatusError("PROTO19 freeze does not validate") from exc


def _process_lines() -> tuple[str, ...]:
    """Return process rows without invoking a shell or changing process state."""
    try:
        completed = subprocess.run(
            ["/bin/ps", "-axo", "pid=,command="],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise Proto19StatusError("runner process status cannot be observed") from exc
    return tuple(line for line in completed.stdout.splitlines() if line.strip())


def _runner_status(lines: Iterable[str]) -> dict[str, Any]:
    """Conservative process classification: a textual match means active."""
    matches: list[dict[str, Any]] = []
    for raw in lines:
        if not any(marker in raw for marker in RUNNER_MARKERS):
            continue
        parts = raw.strip().split(None, 1)
        pid: int | None
        try:
            pid = int(parts[0])
        except (IndexError, ValueError):
            pid = None
        matches.append({"pid": pid, "command": parts[1] if len(parts) == 2 else raw.strip()})
    return {
        "active": bool(matches),
        "classification": "active_runner_detected" if matches else "no_runner_detected",
        "matches": matches,
    }


def _time(value: Mapping[str, Any]) -> dict[str, str]:
    rational = value.get("rational")
    binary = value.get("binary64_hex")
    if not isinstance(rational, str) or not isinstance(binary, str):
        raise Proto19StatusError("time identity differs")
    return {"rational": rational, "binary64_hex": binary}


def _member_status(cursors: Mapping[str, Mapping[str, Any]]) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    if set(cursors) != set(MEMBER_KEYS):
        raise Proto19StatusError("member cursor inventory differs")
    members: dict[str, Any] = {}
    targets: set[str] = set()
    target_values: dict[str, dict[str, str]] = {}
    pending: dict[str, Any] = {}
    for key in MEMBER_KEYS:
        cursor = cursors[key]
        accepted = _time(cursor["accepted_boundary_time"])
        target = _time(cursor["event_target_time"])
        mode = cursor.get("mode")
        if not isinstance(mode, str):
            raise Proto19StatusError("member cursor mode differs")
        members[key] = {"accepted_time": accepted, "mode": mode}
        targets.add(_canonical(target))
        target_values[_canonical(target)] = target
        if mode == "RETRY_PENDING":
            pending[key] = {"mode": mode, "retry_successor": cursor.get("retry_successor_payload_or_none")}
    if len(targets) != 1:
        raise Proto19StatusError("member target times differ")
    return members, target_values[next(iter(targets))], pending


def _legacy_generation_zero(root: Path, store: HLT16CampaignStore, runner: Mapping[str, Any]) -> dict[str, Any]:
    """Validate the immutable eight-leaf predecessor without opening arrays."""
    try:
        source = store._legacy_inventory(repair_stages=False)
        leaves, _directories = _walk_tree(store.root)
        expected = {
            "receipts/generation-zero.json",
            f"checkpoints/00000000000000000000-{SEALED_GEN0_CHECKPOINT_SHA256}.json",
            *(f"states/{address}.json" for address in source.values()),
        }
        if set(leaves) != expected:
            raise Proto19StatusError("legacy GEN0 store has runtime, staging, or foreign suffixes")
        checkpoint_raw = read_nofollow(
            store.root,
            f"checkpoints/00000000000000000000-{SEALED_GEN0_CHECKPOINT_SHA256}.json",
            "GEN0 checkpoint",
        )
        checkpoint = validate_proto17_checkpoint(json.loads(checkpoint_raw.decode("ascii")))
        cursors = checkpoint.get("cursors")
        if not isinstance(cursors, Mapping):
            raise Proto19StatusError("GEN0 cursors are absent")
        members, target, pending = _member_status(cursors)  # type: ignore[arg-type]
    except (HLT16CampaignStoreError, OSError, UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        if isinstance(exc, Proto19StatusError):
            raise
        raise Proto19StatusError("legacy GEN0 store does not validate") from exc

    safe = not bool(runner["active"]) and not pending
    return {
        "schema": "FGC-1-PRO19-event1-status-v1",
        "state": "verified_generation_zero_safe_to_start" if safe else "verified_generation_zero_not_startable",
        "authorization": {
            "protocol_artifact_id": "FGC-2-SF1-PROTO18",
            "freeze_artifact_id": "FGC-1-PRO19-FRZ1",
            "sealed_foundation_commit": "142733cc2a26407e7150962879fb6a760001ea18",
        },
        "checkpoint": {"generation": 0, "sha256": SEALED_GEN0_CHECKPOINT_SHA256},
        "members": members,
        "active_target": target,
        "pending_retry": pending,
        "terminal_lock": False,
        "last_complete_journal": {"sequence": -1, "sha256": ROOT_JOURNAL_SHA256},
        "active_write_lock": False,
        "writer_state": None,
        "runner_process": dict(runner),
        "safe_to_start": safe,
        "safe_to_restart": False,
    }


def _hlt16_status(store: HLT16CampaignStore, runner: Mapping[str, Any]) -> dict[str, Any]:
    try:
        snapshot = store.authenticated_snapshot()
        writer = store._active_writer()
        writer_state = store._writer_liveness(writer) if writer is not None else None
        members: dict[str, Any] = {}
        pending: dict[str, Any] = {}
        for key in MEMBER_KEYS:
            member = snapshot.checkpoint.members[key]
            cursor = member.cursor
            mode = cursor.get("mode")
            members[key] = {"accepted_time": _time(cursor["accepted_boundary_time"]), "mode": mode}
            if mode == "RETRY_PENDING" or member.pending_owner is not None:
                pending[key] = {
                    "mode": mode,
                    "owner": member.pending_owner,
                    "cap_hex": member.pending_cap_hex,
                    "retry_successor": cursor.get("retry_successor_payload_or_none"),
                }
        checkpoint = snapshot.checkpoint
        terminal = snapshot.terminal_lock_present
        terminal_checkpoint = checkpoint.disposition in {"scientific_terminal", "invalid_terminal"}
        blocked_by_suffix = snapshot.suffix_classification != "clean" or bool(snapshot.staging_paths)
        safe = not blocked_by_suffix and writer is None and not runner["active"] and not terminal and not terminal_checkpoint
        if terminal:
            state = "terminal"
        elif terminal_checkpoint:
            state = "terminal_checkpoint_unlocked"
        elif blocked_by_suffix:
            state = f"recovery_required:{snapshot.suffix_classification}"
        elif writer is not None:
            state = "active_write_lock"
        elif runner["active"]:
            state = "runner_process_detected_without_store_lock"
        elif pending:
            state = "checkpoint_ready_with_pending_retry"
        else:
            state = "clean_checkpoint"
        return {
            "schema": "FGC-1-PRO19-event1-status-v1",
            "state": state,
            "authorization": {
                "authorization_commit": checkpoint.authorization_commit,
                "plan_sha256": checkpoint.plan_sha256,
                "campaign_id": checkpoint.campaign_id,
                "protocol": checkpoint.protocol,
            },
            "checkpoint": {"generation": checkpoint.generation, "sha256": checkpoint.sha256},
            "members": members,
            "active_target": _time(checkpoint.target),
            "pending_retry": pending,
            "terminal_lock": terminal,
            "terminal_checkpoint": terminal_checkpoint,
            "last_complete_journal": {"sequence": checkpoint.journal_sequence, "sha256": checkpoint.journal_tip_sha256},
            "suffix": {"classification": snapshot.suffix_classification, "staging_paths": list(snapshot.staging_paths)},
            "active_write_lock": writer is not None,
            "writer_state": writer_state,
            "runner_process": dict(runner),
            "safe_to_start": False,
            "safe_to_restart": safe,
        }
    except (HLT16CampaignStoreError, KeyError, TypeError, ValueError) as exc:
        if isinstance(exc, Proto19StatusError):
            raise
        raise Proto19StatusError("HLT16 campaign store does not validate") from exc


def collect_status(
    root: Path = ROOT,
    *,
    process_lines: Iterable[str] | None = None,
    store_root: Path | None = None,
) -> dict[str, Any]:
    """Return one read-only, fail-closed snapshot of the authorized store."""
    root = Path(root)
    _freeze(root)
    runner = _runner_status(_process_lines() if process_lines is None else process_lines)
    store = HLT16CampaignStore(root / STORE_RELATIVE if store_root is None else Path(store_root))
    try:
        snapshot = store.authenticated_snapshot()
    except HLT16CampaignStoreError as exc:
        # GEN0 has no HLT16 bridge by design.  Its eight leaves are checked by
        # the explicit legacy path below; every other HLT16 validation error is
        # reported there only if that exact predecessor remains intact.
        if "bridge checkpoint absent" not in str(exc):
            try:
                return _legacy_generation_zero(root, store, runner)
            except Proto19StatusError:
                raise Proto19StatusError("campaign store does not validate") from exc
        return _legacy_generation_zero(root, store, runner)
    del snapshot
    return _hlt16_status(store, runner)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    try:
        result = collect_status(args.root)
    except Proto19StatusError as exc:
        print(_canonical({
            "schema": "FGC-1-PRO19-event1-status-v1",
            "state": "invalid_status_observation",
            "safe_to_start": False,
            "safe_to_restart": False,
            "detail": str(exc),
        }))
        return 2
    print(_canonical(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
