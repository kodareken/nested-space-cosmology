#!/usr/bin/env python3
"""Fast byte-identical continuation of the RSRC3 SSPRK3-16385 chain."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
SOURCE = str(ROOT / "src")
if SOURCE not in sys.path:
    sys.path.insert(0, SOURCE)
ROOT_TEXT = str(ROOT)
if ROOT_TEXT not in sys.path:
    sys.path.insert(0, ROOT_TEXT)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution.pro20_rsrc3_profile import (  # noqa: E402
    CAMPAIGN_ID,
    activate_rsrc3,
)


activate_rsrc3()

from recursive_horizons.fgc.evolution.pro20_rsrc1_isolation import (  # noqa: E402
    decode_bundle,
    encode_attempt_request,
    launch_child_process,
    reduce_completed_child,
)
from recursive_horizons.fgc.evolution.pro20_rsrc1_runtime import (  # noqa: E402
    _attempt_request,
)
from recursive_horizons.fgc.evolution.protocol_v19 import TERMINAL_KINDS  # noqa: E402
from scripts import run_fgc_pro20_rsrc3_parallel as rsrc3  # noqa: E402


ARTIFACT_ID = "FGC-1-PRO20-EV1-RSRC4-FAST1"
SCHEMA = "FGC-1-PRO20-EV1-RSRC4-fast-continuation-v1"
MEMBER_KEY = "SSPRK3-16385"
SOURCE_RESPONSE = (
    "runs/fgc-2-sf1/pro20-rsrc3-event1/event/members/SSPRK3-16385/"
    "attempt-000027-response.json"
)
SOURCE_RESPONSE_SHA256 = (
    "5596cff66b8189e5cf793f6e5514b2a286d8f2c35676419cdcda534273c40801"
)
OUTPUT_CONTAINER = "runs/fgc-2-sf1/pro20-rsrc4-fast1"
OUTPUT_NAMESPACE = f"{OUTPUT_CONTAINER}/event"
CHILD_SCRIPT = "scripts/run_fgc_pro20_rsrc4_fast_child.py"
RUNNER_PATH = "scripts/run_fgc_pro20_rsrc4_fast_continuation.py"
IMPLEMENTATION_PATHS = (
    *rsrc3.IMPLEMENTATION_PATHS,
    CHILD_SCRIPT,
    RUNNER_PATH,
)


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _git(*arguments: str) -> str:
    completed = subprocess.run(
        ("git", *arguments),
        cwd=ROOT,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
        timeout=30.0,
    )
    if completed.returncode != 0:
        raise RuntimeError("RSRC4 Git authority query failed")
    return completed.stdout.strip()


def _implementation_hash(commit: str | None = None) -> str:
    rows: list[tuple[str, str]] = []
    for relative in IMPLEMENTATION_PATHS:
        if commit is None:
            raw = (ROOT / relative).read_bytes()
        else:
            completed = subprocess.run(
                ("git", "show", f"{commit}:{relative}"),
                cwd=ROOT,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            if completed.returncode != 0:
                raise RuntimeError("RSRC4 committed implementation is incomplete")
            raw = completed.stdout
        rows.append((relative, _sha(raw)))
    return _sha("".join(f"{path}\0{digest}\n" for path, digest in rows).encode("ascii"))


def _source_bundle():
    raw = (ROOT / SOURCE_RESPONSE).read_bytes()
    if _sha(raw) != SOURCE_RESPONSE_SHA256:
        raise RuntimeError("RSRC4 source response hash differs")
    value = json.loads(raw.decode("ascii"))
    if value.get("member_key") != MEMBER_KEY or value.get("stop_kind") != "none":
        raise RuntimeError("RSRC4 source response is not an accepted member response")
    bundle = decode_bundle(value["successor_bundle"])
    if bundle.cursor_mapping.get("accepted_time_hex") != "0x1.75d55555554c0p+0":
        raise RuntimeError("RSRC4 source accepted time differs")
    return bundle


def status(authority_commit: str) -> dict[str, object]:
    if len(authority_commit) != 40 or any(c not in "0123456789abcdef" for c in authority_commit):
        raise RuntimeError("RSRC4 authority commit is not exact")
    if _git("rev-parse", "HEAD") != authority_commit:
        raise RuntimeError("RSRC4 authority is not live HEAD")
    if _git("status", "--porcelain=v1", "--untracked-files=all"):
        raise RuntimeError("RSRC4 authority requires a clean worktree")
    live = _implementation_hash()
    committed = _implementation_hash(authority_commit)
    if live != committed:
        raise RuntimeError("RSRC4 live implementation differs from authority")
    _source_bundle()
    if (ROOT / OUTPUT_CONTAINER).exists():
        raise RuntimeError("RSRC4 one-shot namespace already exists")
    return {
        "artifact_id": ARTIFACT_ID,
        "authority_commit": authority_commit,
        "campaign_id": CAMPAIGN_ID,
        "implementation_sha256": live,
        "member_key": MEMBER_KEY,
        "output_namespace": OUTPUT_NAMESPACE,
        "source_response": SOURCE_RESPONSE,
        "source_response_sha256": SOURCE_RESPONSE_SHA256,
        "source_accepted_generation": 28,
        "byte_identical_fast_decimal": True,
        "safe_to_run": True,
        "physics_claimed": False,
    }


def run(authority_commit: str) -> dict[str, object]:
    preflight = status(authority_commit)
    bundle = _source_bundle()
    receipt, source_authority = rsrc3._source_authority(
        str(preflight["implementation_sha256"])
    )
    container = ROOT / OUTPUT_CONTAINER
    container.mkdir()
    event = ROOT / OUTPUT_NAMESPACE
    event.mkdir()
    rsrc3._write_exclusive(
        event / "root.json",
        canonical_json_bytes({**preflight, "schema": SCHEMA, "terminal": False}),
    )
    records: list[dict[str, object]] = []
    started = time.monotonic()
    command = (sys.executable, "-I", "-B", os.fspath(ROOT / CHILD_SCRIPT))
    for index in range(4096):
        accepted_time = str(bundle.cursor_mapping["accepted_time_hex"])
        if accepted_time == float(3.0 / 2.0).hex():
            terminal = {
                "artifact_id": ARTIFACT_ID,
                "schema": SCHEMA,
                "authority_commit": authority_commit,
                "campaign_id": CAMPAIGN_ID,
                "member_key": MEMBER_KEY,
                "disposition": "member_common_event_complete",
                "accepted_generation": bundle.accepted_generation,
                "accepted_time_hex": accepted_time,
                "physical_state_sha256": bundle.physical_state_sha256,
                "attempts": records,
                "terminal": True,
                "wall_seconds_hex": (time.monotonic() - started).hex(),
                "physics_claimed": False,
            }
            rsrc3._write_exclusive(
                event / "terminal.json", canonical_json_bytes(terminal)
            )
            return terminal
        request = _attempt_request(
            MEMBER_KEY,
            bundle,
            receipt=receipt,
            authority=source_authority,
        )
        completed = launch_child_process(
            command,
            repository_root=ROOT,
            request_raw=encode_attempt_request(request),
        )
        classification, metrics = reduce_completed_child(request, completed)
        response_sha256 = rsrc3._write_exclusive(
            event / f"attempt-{index:06d}-response.json", completed.stdout
        )
        bundle = classification.successor_bundle
        records.append(
            {
                "attempt": index,
                "kind": classification.kind,
                "request_sha256": request.request_sha256,
                "response_sha256": response_sha256,
                "accepted_generation": bundle.accepted_generation,
                "accepted_time_hex": bundle.cursor_mapping["accepted_time_hex"],
                "child_peak_rss_bytes": metrics.peak_rss_bytes,
                "child_wall_seconds_hex": metrics.wall_seconds.hex(),
            }
        )
        if classification.kind in TERMINAL_KINDS:
            terminal = {
                "artifact_id": ARTIFACT_ID,
                "schema": SCHEMA,
                "authority_commit": authority_commit,
                "campaign_id": CAMPAIGN_ID,
                "member_key": MEMBER_KEY,
                "disposition": classification.kind,
                "accepted_generation": bundle.accepted_generation,
                "accepted_time_hex": bundle.cursor_mapping["accepted_time_hex"],
                "physical_state_sha256": bundle.physical_state_sha256,
                "attempts": records,
                "terminal_evidence": classification.terminal_evidence,
                "terminal": True,
                "physics_claimed": False,
            }
            rsrc3._write_exclusive(
                event / "terminal.json", canonical_json_bytes(terminal)
            )
            return terminal
    raise RuntimeError("RSRC4 attempt ceiling exceeded")


def _print(value: object) -> None:
    sys.stdout.write(json.dumps(value, sort_keys=True, indent=2) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=ARTIFACT_ID)
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--authority-commit", required=True)
    args = parser.parse_args(argv)
    if args.status == args.run:
        _print({"error": "select exactly one of --status or --run"})
        return 2
    try:
        _print(status(args.authority_commit) if args.status else run(args.authority_commit))
    except Exception as error:
        _print({"error": f"{type(error).__name__}: {error}", "safe_to_run": False})
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
