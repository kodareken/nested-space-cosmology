#!/usr/bin/env python3
"""Six-way RSRC3 development common-event calculation.

Every method/resolution member owns an independent append-only response spool
and launches the already-qualified fresh child for each scheduled attempt.  No
member consumes another member's state.  The parent publishes one terminal only
after all six worker chains return.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from hashlib import sha256
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import time
import tomllib


ROOT = Path(__file__).resolve().parents[1]
SOURCE = str(ROOT / "src")
if SOURCE not in sys.path:
    sys.path.insert(0, SOURCE)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution.pro20_rsrc3_profile import (  # noqa: E402
    ARTIFACT_ID,
    CAMPAIGN_ID,
    FREEZE_ARTIFACT_ID,
    PRODUCTION_CONTAINER,
    PRODUCTION_NAMESPACE,
    activate_rsrc3,
)


activate_rsrc3()

from recursive_horizons.fgc.evolution.pro20_origin import MEMBER_KEYS  # noqa: E402
from recursive_horizons.fgc.evolution.pro20_rsrc1_protocol import (  # noqa: E402
    build_authority_receipt,
)
from recursive_horizons.fgc.evolution.pro20_rsrc1_runtime import (  # noqa: E402
    _attempt_request,
    _launch_attempt,
    _launch_seed,
    _seed_request,
)
from recursive_horizons.fgc.evolution.pro20_rsrc1_isolation import (  # noqa: E402
    reduce_completed_child,
)
from recursive_horizons.fgc.evolution.pro20_rsrc1_seed import (  # noqa: E402
    reduce_completed_seed_child,
)
from recursive_horizons.fgc.evolution.protocol_v19 import TERMINAL_KINDS  # noqa: E402


SCHEMA = "FGC-1-PRO20-EV1-RSRC3-parallel-event-v1"
RSRC2_CONFIG = "configs/fgc/fgc-1-pro20-ev1-rsrc2-frz1.toml"
MAX_ATTEMPTS_PER_MEMBER = 4096
MAX_WORKERS = len(MEMBER_KEYS)
IMPLEMENTATION_PATHS = (
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_attempt.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_isolation.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_member.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_protocol.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_runtime.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc1_seed.py",
    "src/recursive_horizons/fgc/evolution/pro20_rsrc3_profile.py",
    "scripts/run_fgc_pro20_rsrc3_child.py",
    "scripts/run_fgc_pro20_rsrc3_parallel.py",
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
        check=False,
        text=True,
        timeout=30.0,
    )
    if completed.returncode != 0:
        raise RuntimeError("RSRC3 Git authority query failed")
    return completed.stdout.strip()


def _implementation_hashes(commit: str | None = None) -> tuple[dict[str, str], str]:
    hashes: dict[str, str] = {}
    for relative in IMPLEMENTATION_PATHS:
        raw = (
            subprocess.run(
                ("git", "show", f"{commit}:{relative}"),
                cwd=ROOT,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            ).stdout
            if commit is not None
            else (ROOT / relative).read_bytes()
        )
        hashes[relative] = _sha(raw)
    aggregate = _sha(
        "".join(f"{path}\0{hashes[path]}\n" for path in IMPLEMENTATION_PATHS).encode(
            "ascii"
        )
    )
    return hashes, aggregate


def _source_authority(implementation_sha256: str) -> tuple[object, dict[str, object]]:
    config_raw = (ROOT / RSRC2_CONFIG).read_bytes()
    config = tomllib.loads(config_raw.decode("utf-8"))
    source_rows = config["source_configuration"]
    source_configuration = {
        str(row["member_key"]): str(row["sha256"]) for row in source_rows
    }
    specification = {
        "artifact_id": FREEZE_ARTIFACT_ID,
        "campaign_id": CAMPAIGN_ID,
        "implementation_sha256": implementation_sha256,
        "max_workers": MAX_WORKERS,
        "output_namespace": PRODUCTION_NAMESPACE,
        "source_closure_sha256": config["source_closure_sha256"],
        "origin_capture_sha256": config["origin_capture_sha256"],
    }
    config_sha256 = _sha(canonical_json_bytes(specification))
    runner_sha256 = _sha((ROOT / "scripts/run_fgc_pro20_rsrc3_parallel.py").read_bytes())
    environment_sha256 = _sha(canonical_json_bytes(config["environment"]))
    receipt = build_authority_receipt(
        campaign_id=CAMPAIGN_ID,
        authority_sha256=runner_sha256,
        implementation_sha256=implementation_sha256,
        config_sha256=config_sha256,
        source_sha256=str(config["source_closure_sha256"]),
        origin_sha256=str(config["origin_capture_sha256"]),
        environment_sha256=environment_sha256,
    )
    authority = {
        "historical_source_tree": {
            "generation1_bridge_content_id": config["historical"][
                "generation1_bridge_content_id"
            ]
        },
        "source_configuration": source_configuration,
    }
    return receipt, authority


def _namespace_absent() -> bool:
    target = ROOT.joinpath(*PRODUCTION_CONTAINER.split("/"))
    try:
        target.lstat()
    except FileNotFoundError:
        return True
    return False


def status(authority_commit: str) -> dict[str, object]:
    if len(authority_commit) != 40 or any(c not in "0123456789abcdef" for c in authority_commit):
        raise RuntimeError("RSRC3 authority commit is not exact")
    head = _git("rev-parse", "HEAD")
    if head != authority_commit:
        raise RuntimeError("RSRC3 authority is not live HEAD")
    if _git("status", "--porcelain=v1", "--untracked-files=all"):
        raise RuntimeError("RSRC3 authority requires a clean worktree")
    live_hashes, live_aggregate = _implementation_hashes()
    committed_hashes, committed_aggregate = _implementation_hashes(authority_commit)
    if live_hashes != committed_hashes or live_aggregate != committed_aggregate:
        raise RuntimeError("RSRC3 live implementation differs from authority")
    if not _namespace_absent():
        raise RuntimeError("RSRC3 one-shot namespace already exists")
    receipt, _authority = _source_authority(live_aggregate)
    return {
        "artifact_id": FREEZE_ARTIFACT_ID,
        "authority_commit": authority_commit,
        "campaign_id": CAMPAIGN_ID,
        "implementation_sha256": live_aggregate,
        "max_workers": MAX_WORKERS,
        "member_keys": list(MEMBER_KEYS),
        "output_namespace": PRODUCTION_NAMESPACE,
        "receipt_sha256": receipt.sha256,
        "safe_to_run": True,
        "parallel_member_chains": True,
        "campaign_execution_authorized": False,
        "physics_claimed": False,
    }


def _mkdir(path: Path) -> None:
    path.mkdir(mode=0o755)
    info = path.lstat()
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
        raise RuntimeError("RSRC3 output directory is unsafe")


def _write_exclusive(path: Path, raw: bytes) -> str:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(path, flags, 0o644)
    try:
        view = memoryview(raw)
        while view:
            written = os.write(fd, view)
            if written <= 0:
                raise OSError("short RSRC3 evidence write")
            view = view[written:]
        os.fsync(fd)
    finally:
        os.close(fd)
    return _sha(raw)


def _accepted_time_hex(bundle) -> str:
    value = bundle.cursor_mapping.get("accepted_time_hex")
    if type(value) is not str:
        raise RuntimeError("RSRC3 bundle omitted accepted time")
    return value


def _run_member(member_key: str, receipt, authority: dict[str, object], event: Path) -> dict[str, object]:
    member_dir = event / "members" / member_key
    _mkdir(member_dir)
    request = _seed_request(member_key, receipt=receipt, authority=authority)
    completed = _launch_seed(ROOT, request)
    outcome = reduce_completed_seed_child(request, completed)
    seed_sha256 = _write_exclusive(member_dir / "seed-response.json", completed.stdout)
    if not outcome.successful or outcome.bundle is None:
        result = {
            "member_key": member_key,
            "disposition": "seed_resource_stop",
            "resource_evidence": outcome.resource_evidence,
            "seed_response_sha256": seed_sha256,
            "terminal": True,
        }
        _write_exclusive(member_dir / "terminal.json", canonical_json_bytes(result))
        return result
    bundle = outcome.bundle
    records: list[dict[str, object]] = []
    for index in range(MAX_ATTEMPTS_PER_MEMBER):
        if _accepted_time_hex(bundle) == float(3.0 / 2.0).hex():
            result = {
                "member_key": member_key,
                "disposition": "member_common_event_complete",
                "accepted_generation": bundle.accepted_generation,
                "accepted_time_hex": _accepted_time_hex(bundle),
                "physical_state_sha256": bundle.physical_state_sha256,
                "seed_response_sha256": seed_sha256,
                "attempts": records,
                "terminal": True,
            }
            _write_exclusive(member_dir / "terminal.json", canonical_json_bytes(result))
            return result
        request = _attempt_request(
            member_key,
            bundle,
            receipt=receipt,
            authority=authority,
        )
        completed = _launch_attempt(ROOT, request)
        classification, metrics = reduce_completed_child(request, completed)
        response_sha256 = _write_exclusive(
            member_dir / f"attempt-{index:06d}-response.json",
            completed.stdout,
        )
        bundle = classification.successor_bundle
        records.append(
            {
                "attempt": index,
                "kind": classification.kind,
                "request_sha256": request.request_sha256,
                "response_sha256": response_sha256,
                "accepted_generation": bundle.accepted_generation,
                "accepted_time_hex": _accepted_time_hex(bundle),
                "child_peak_rss_bytes": metrics.peak_rss_bytes,
                "child_wall_seconds_hex": metrics.wall_seconds.hex(),
            }
        )
        if classification.kind in TERMINAL_KINDS:
            result = {
                "member_key": member_key,
                "disposition": classification.kind,
                "accepted_generation": bundle.accepted_generation,
                "accepted_time_hex": _accepted_time_hex(bundle),
                "physical_state_sha256": bundle.physical_state_sha256,
                "seed_response_sha256": seed_sha256,
                "attempts": records,
                "terminal_evidence": classification.terminal_evidence,
                "terminal": True,
            }
            _write_exclusive(member_dir / "terminal.json", canonical_json_bytes(result))
            return result
    raise RuntimeError(f"RSRC3 member attempt ceiling exceeded: {member_key}")


def run(authority_commit: str) -> dict[str, object]:
    preflight = status(authority_commit)
    receipt, authority = _source_authority(str(preflight["implementation_sha256"]))
    container = ROOT.joinpath(*PRODUCTION_CONTAINER.split("/"))
    _mkdir(container)
    event = ROOT.joinpath(*PRODUCTION_NAMESPACE.split("/"))
    _mkdir(event)
    _mkdir(event / "members")
    started = time.monotonic()
    root_record = {
        **preflight,
        "schema": SCHEMA,
        "started_monotonic_hex": started.hex(),
        "terminal": False,
    }
    _write_exclusive(event / "root.json", canonical_json_bytes(root_record))
    results: dict[str, dict[str, object]] = {}
    errors: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=MAX_WORKERS, thread_name_prefix="rsrc3") as pool:
        futures = {
            pool.submit(_run_member, key, receipt, authority, event): key
            for key in MEMBER_KEYS
        }
        for future in as_completed(futures):
            key = futures[future]
            try:
                results[key] = future.result()
            except Exception as error:
                errors[key] = f"{type(error).__name__}: {error}"
    ordered = {key: results[key] for key in MEMBER_KEYS if key in results}
    complete = not errors and all(
        ordered[key]["disposition"] == "member_common_event_complete"
        for key in MEMBER_KEYS
    )
    terminal = {
        "artifact_id": ARTIFACT_ID,
        "schema": SCHEMA,
        "authority_commit": authority_commit,
        "campaign_id": CAMPAIGN_ID,
        "members": ordered,
        "errors": errors,
        "first_event_complete": complete,
        "disposition": "six_member_common_event_completed" if complete else "parallel_member_terminal",
        "terminal": True,
        "wall_seconds_hex": (time.monotonic() - started).hex(),
        "physics_claimed": False,
    }
    _write_exclusive(event / "terminal.json", canonical_json_bytes(terminal))
    return terminal


def _print(value: object) -> None:
    sys.stdout.write(json.dumps(value, sort_keys=True, indent=2) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=FREEZE_ARTIFACT_ID)
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
