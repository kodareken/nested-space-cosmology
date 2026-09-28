#!/usr/bin/env python3
"""Build or verify the premise-only PRO19 CFL1 continuation certificate."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
)
from recursive_horizons.fgc.evolution.proto19_cfl_continuation_authority import (  # noqa: E402
    CONFIG_PATH,
    RESULT_PATH,
    build_cfl1_result,
)


STORE = Path("runs/fgc-2-sf1/proto17/calibration")


def canonical(value: object) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            indent=2,
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("ascii")


def derive_store_anchor(repository: Path) -> dict[str, object]:
    """Derive the exact nonmutating generation-six recovery facts."""
    store = HLT16CampaignStore(repository / STORE)
    snapshot = store.authenticated_snapshot()
    # A certificate verifier must observe a crash suffix, not repair it.
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


def expected_result(repository: Path) -> bytes:
    config = (repository / CONFIG_PATH).read_bytes()
    return canonical(build_cfl1_result(config, derive_store_anchor(repository)))


def write_result(repository: Path) -> None:
    destination = repository / RESULT_PATH
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
    temporary.write_bytes(expected_result(repository))
    os.replace(temporary, destination)


def verify_result(repository: Path) -> None:
    if (repository / RESULT_PATH).read_bytes() != expected_result(repository):
        raise SystemExit("PRO19 CFL1 continuation result differs")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        write_result(ROOT)
    else:
        verify_result(ROOT)
    print(
        json.dumps(
            {
                "artifact_id": "FGC-1-PRO19-CFL1-FRZ1",
                "verified": bool(args.verify),
                "state_advanced": False,
                "candidate_branch_opened": False,
                "physical_result_earned": False,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
