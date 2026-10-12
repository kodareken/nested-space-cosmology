#!/usr/bin/env python3
"""Result-first RK4-4097 continuation to the first raw trapped GR slice."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import time
import tomllib

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution.pro20_rsrc3_profile import (  # noqa: E402
    CAMPAIGN_ID,
    activate_rsrc3,
)


activate_rsrc3()

from recursive_horizons.fgc.evolution import tdg11_c1r1_enclosure as c1r1_enclosure  # noqa: E402
from recursive_horizons.fgc.evolution import tdg11_c1r1_ring as c1r1_ring  # noqa: E402
from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    gr0_semidiscrete_constraint_snapshot,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    radial_null_observables,
)
from recursive_horizons.fgc.evolution.hlt17_imp1_cursor import (  # noqa: E402
    HLT17_GENESIS_DIGEST,
    encode_hlt17_cursor,
    restore_hlt17_cursor,
)
from recursive_horizons.fgc.evolution.pro20_origin import (  # noqa: E402
    capture_pro20_historical_origin,
)
from recursive_horizons.fgc.evolution.pro20_rsrc1_attempt import (  # noqa: E402
    execute_scheduled_attempt,
    restore_checkpoint_from_bundle,
)
from recursive_horizons.fgc.evolution.pro20_rsrc1_isolation import (  # noqa: E402
    _load_static_inputs,
    decode_bundle,
    encode_bundle,
)
from recursive_horizons.fgc.evolution.pro20_rsrc1_member import (  # noqa: E402
    build_rsrc1_member,
)
from recursive_horizons.fgc.evolution.protocol_v19 import TERMINAL_KINDS  # noqa: E402
from scripts import run_fgc_pro20_rsrc3_parallel as rsrc3  # noqa: E402


ARTIFACT_ID = "FGC-1-GR0-TRAP-DEV1"
SCHEMA = "FGC-1-GR0-TRAP-DEV1-v1"
MEMBER_KEY = "RK4-4097"
SOURCE_RESPONSE = (
    "runs/fgc-2-sf1/pro20-rsrc3-event1/event/members/RK4-4097/"
    "attempt-000044-response.json"
)
SOURCE_RESPONSE_SHA256 = (
    "debb56a8f2eae1aad34227762f7b919ee2088e2e53a94b27448be50255a288e1"
)
OUTPUT_CONTAINER = "runs/fgc-2-sf1/gr0-trap-dev1"
OUTPUT_NAMESPACE = f"{OUTPUT_CONTAINER}/event"
TARGET_TIME = 32.0
MAX_ATTEMPTS = 4096
RUNNER_PATH = "scripts/run_fgc_gr0_trapped_dev1.py"
IMPLEMENTATION_PATHS = (
    *rsrc3.IMPLEMENTATION_PATHS,
    RUNNER_PATH,
)


_chunked_decimal = c1r1_ring.int_to_decimal_text


def _fast_decimal(value: int) -> str:
    if type(value) is not int:
        raise TypeError("integer text requires an int")
    try:
        return str(value)
    except ValueError:
        return _chunked_decimal(value)


c1r1_ring.int_to_decimal_text = _fast_decimal
c1r1_enclosure.int_to_decimal_text = _fast_decimal


def _sha(raw: bytes) -> str:
    return sha256(raw).hexdigest()


def _git(*arguments: str) -> str:
    completed = subprocess.run(
        ("git", *arguments), cwd=ROOT, capture_output=True, text=True, check=False
    )
    if completed.returncode != 0:
        raise RuntimeError("TRAP-DEV1 Git authority query failed")
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
                capture_output=True,
                check=False,
            )
            if completed.returncode != 0:
                raise RuntimeError("TRAP-DEV1 committed implementation is incomplete")
            raw = completed.stdout
        rows.append((relative, _sha(raw)))
    return _sha("".join(f"{path}\0{digest}\n" for path, digest in rows).encode("ascii"))


def _complete_terminal(path: str) -> bool:
    target = ROOT / path
    if not target.is_file():
        return False
    value = json.loads(target.read_text())
    return value.get("disposition") == "member_common_event_complete"


def _six_member_common_event_ready() -> bool:
    members = "runs/fgc-2-sf1/pro20-rsrc3-event1/event/members"
    fixed = (
        "RK4-2049",
        "RK4-4097",
        "RK4-8193",
        "SSPRK3-4097",
        "SSPRK3-8193",
    )
    if not all(_complete_terminal(f"{members}/{key}/terminal.json") for key in fixed):
        return False
    return _complete_terminal(f"{members}/SSPRK3-16385/terminal.json") or _complete_terminal(
        "runs/fgc-2-sf1/pro20-rsrc4-fast1/event/terminal.json"
    )


def _source_bundle():
    raw = (ROOT / SOURCE_RESPONSE).read_bytes()
    if _sha(raw) != SOURCE_RESPONSE_SHA256:
        raise RuntimeError("TRAP-DEV1 source response hash differs")
    value = json.loads(raw.decode("ascii"))
    bundle = decode_bundle(value["successor_bundle"])
    if bundle.cursor_mapping.get("accepted_time_hex") != float(3.0 / 2.0).hex():
        raise RuntimeError("TRAP-DEV1 source is not the common-event boundary")
    return bundle


def _retarget(checkpoint) -> None:
    previous = checkpoint.cursor
    mapping = encode_hlt17_cursor(previous, recompute_tip=False)
    mapping["event_target_hex"] = TARGET_TIME.hex()
    mapping["kernel_transition_parent_sha256"] = previous.kernel_transition_digest
    mapping["kernel_transition_digest"] = HLT17_GENESIS_DIGEST
    mapping["kernel_transition_digest"] = _sha(canonical_json_bytes(mapping))
    checkpoint.cursor = restore_hlt17_cursor(mapping)
    checkpoint.agree()


def _measure(checkpoint, previous_mass: np.ndarray | None, previous_time: float | None):
    state = checkpoint.member.state
    observables = radial_null_observables(state)
    positive_radius = state.u[:, 3] > 0.0
    trapped = positive_radius & (observables.theta_plus < 0.0) & (
        observables.theta_minus < 0.0
    )
    indices = np.flatnonzero(trapped)
    sample: dict[str, object] = {
        "accepted_generation": checkpoint.accepted_generation,
        "coordinate_time_hex": checkpoint.member.time.hex(),
        "minimum_theta_plus_hex": float(np.min(observables.theta_plus[positive_radius])).hex(),
        "minimum_theta_minus_hex": float(np.min(observables.theta_minus[positive_radius])).hex(),
        "maximum_compactness_hex": float(np.max(observables.compactness)).hex(),
        "maximum_misner_sharp_mass_hex": float(np.max(observables.misner_sharp_mass)).hex(),
        "raw_trapped": bool(indices.size),
        "first_trapped_index": None if not indices.size else int(indices[0]),
        "last_trapped_index": None if not indices.size else int(indices[-1]),
    }
    if (
        previous_mass is not None
        and previous_time is not None
        and checkpoint.member.time > previous_time
    ):
        dt = checkpoint.member.time - previous_time
        sample["mass_profile_time_difference_infinity_hex"] = float(
            np.max(np.abs(observables.misner_sharp_mass - previous_mass)) / dt
        ).hex()
    return sample, observables, indices


def status(authority_commit: str) -> dict[str, object]:
    if len(authority_commit) != 40 or any(c not in "0123456789abcdef" for c in authority_commit):
        raise RuntimeError("TRAP-DEV1 authority commit is not exact")
    if _git("rev-parse", "HEAD") != authority_commit:
        raise RuntimeError("TRAP-DEV1 authority is not live HEAD")
    if _git("status", "--porcelain=v1", "--untracked-files=all"):
        raise RuntimeError("TRAP-DEV1 authority requires a clean worktree")
    live = _implementation_hash()
    if live != _implementation_hash(authority_commit):
        raise RuntimeError("TRAP-DEV1 live implementation differs from authority")
    _source_bundle()
    if not _six_member_common_event_ready():
        raise RuntimeError("TRAP-DEV1 requires the six-member common event")
    if (ROOT / OUTPUT_CONTAINER).exists():
        raise RuntimeError("TRAP-DEV1 one-shot namespace already exists")
    return {
        "artifact_id": ARTIFACT_ID,
        "authority_commit": authority_commit,
        "campaign_id": CAMPAIGN_ID,
        "implementation_sha256": live,
        "member_key": MEMBER_KEY,
        "source_response_sha256": SOURCE_RESPONSE_SHA256,
        "source_time_hex": float(3.0 / 2.0).hex(),
        "retarget_time_hex": TARGET_TIME.hex(),
        "output_namespace": OUTPUT_NAMESPACE,
        "six_member_common_event": True,
        "safe_to_run": True,
        "raw_trapped_only": True,
        "physics_claimed": False,
    }


def run(authority_commit: str) -> dict[str, object]:
    preflight = status(authority_commit)
    config = tomllib.loads(
        (ROOT / "configs/fgc/fgc-1-pro20-ev1-rsrc2-frz1.toml").read_text()
    )
    record = build_rsrc1_member(
        capture_pro20_historical_origin(ROOT),
        member_key=MEMBER_KEY,
        static_input_bytes=_load_static_inputs(ROOT),
        source_closure_sha256=config["source_closure_sha256"],
    )
    checkpoint = record.checkpoint
    restore_checkpoint_from_bundle(checkpoint, _source_bundle())
    _retarget(checkpoint)
    container = ROOT / OUTPUT_CONTAINER
    container.mkdir()
    event = ROOT / OUTPUT_NAMESPACE
    event.mkdir()
    rsrc3._write_exclusive(
        event / "root.json",
        canonical_json_bytes({**preflight, "schema": SCHEMA, "terminal": False}),
    )
    attempts: list[dict[str, object]] = []
    previous_mass: np.ndarray | None = None
    previous_time: float | None = None
    for index in range(MAX_ATTEMPTS):
        sample, observables, trapped_indices = _measure(
            checkpoint, previous_mass, previous_time
        )
        rsrc3._write_exclusive(
            event
            / (
                f"sample-{index:06d}-generation-"
                f"{checkpoint.accepted_generation:06d}.json"
            ),
            canonical_json_bytes(sample),
        )
        if trapped_indices.size:
            constraint = gr0_semidiscrete_constraint_snapshot(
                checkpoint.member.state,
                checkpoint.member.initial.grid,
                coordinate_time=checkpoint.member.time,
                diagnostic_spatial_order=4,
            )
            rhs = checkpoint.member.operator(
                checkpoint.member.time, checkpoint.member.state
            )
            terminal = {
                "artifact_id": ARTIFACT_ID,
                "schema": SCHEMA,
                "authority_commit": authority_commit,
                "campaign_id": CAMPAIGN_ID,
                "member_key": MEMBER_KEY,
                "disposition": "first_raw_trapped_slice",
                "sample": sample,
                "theta_plus_hex": [float(item).hex() for item in observables.theta_plus],
                "theta_minus_hex": [float(item).hex() for item in observables.theta_minus],
                "compactness_hex": [float(item).hex() for item in observables.compactness],
                "misner_sharp_mass_hex": [
                    float(item).hex() for item in observables.misner_sharp_mass
                ],
                "constraint_component_names": list(constraint.raw_norms.component_names),
                "constraint_component_infinity_hex": [
                    float(item).hex() for item in constraint.raw_norms.component_infinity
                ],
                "constraint_global_infinity_hex": constraint.raw_norms.global_infinity.hex(),
                "diagnostics": dict(rhs.diagnostics),
                "bundle": encode_bundle(checkpoint.generation_bundle()),
                "attempts": attempts,
                "terminal": True,
                "error_separated_trapped_claimed": False,
                "physics_claimed": False,
            }
            rsrc3._write_exclusive(
                event / "terminal.json", canonical_json_bytes(terminal)
            )
            return terminal
        if previous_time is None or checkpoint.member.time > previous_time:
            previous_mass = np.array(observables.misner_sharp_mass, copy=True)
            previous_time = checkpoint.member.time
        classification = execute_scheduled_attempt(checkpoint)
        attempts.append(
            {
                "attempt": index,
                "kind": classification.kind,
                "accepted_state_advanced": classification.accepted_state_advanced,
                "accepted_generation": checkpoint.accepted_generation,
                "accepted_time_hex": checkpoint.member.time.hex(),
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
                "terminal_evidence": classification.terminal_evidence,
                "attempts": attempts,
                "terminal": True,
                "physics_claimed": False,
            }
            rsrc3._write_exclusive(
                event / "terminal.json", canonical_json_bytes(terminal)
            )
            return terminal
    raise RuntimeError("TRAP-DEV1 attempt ceiling exceeded")


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
