#!/usr/bin/env python3
"""Guarded PRO20-EV1 CLI. Default/status is read-only; no physical source.

A live first-event run requires an exact ``--authority-commit`` together with
``--run``. Authority/delta validation is the later freeze seam, not a boolean
flag. Status never constructs the origin, source, or production namespace.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
_SOURCE_ROOT = str(ROOT / "src")
if _SOURCE_ROOT not in sys.path:
    sys.path.insert(0, _SOURCE_ROOT)

from recursive_horizons.fgc.evolution.pro20_ev1_runtime import (  # noqa: E402
    PRO20EV1RuntimeError,
    implementation_identity,
    plan_status,
    require_authority_delta,
    run_first_event,
)


ARTIFACT_ID = "FGC-1-PRO20-EV1-FRZ1"


def _print(payload: object) -> None:
    sys.stdout.write(
        json.dumps(payload, sort_keys=True, indent=2, ensure_ascii=True) + "\n"
    )


def status() -> dict[str, object]:
    return plan_status()


def committed_image_status(authority_commit: str) -> dict[str, object]:
    """Read-only freeze preflight. Never builds origin, runtime, or source."""

    implementation = implementation_identity(ROOT)
    receipt = require_authority_delta(
        repository_root=ROOT,
        authority_commit=authority_commit,
        implementation_sha256=str(implementation["sha256"]),
    )
    payload = dict(receipt)
    payload["safe_to_run"] = True
    payload["physical_source_constructed"] = False
    payload["physical_source_qualification_executed"] = False
    payload["campaign_execution_authorized"] = False
    payload["calibration_eligible"] = False
    return payload


def _cli_error_payload(error: PRO20EV1RuntimeError) -> dict[str, object]:
    payload: dict[str, object] = {
        "error": str(error),
        "outcome": error.outcome,
        "published": error.published,
        "phase": error.phase,
        "physical_source_constructed": bool(error.physical_source_constructed),
        "safe_to_run": False,
        "campaign_execution_authorized": False,
        "calibration_eligible": False,
    }
    if error.abandon_writer:
        payload["unclosed_session"] = True
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="FGC-1-PRO20-EV1 plan/status or one guarded first-event run"
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="read-only plan or committed-image preflight; never builds source",
    )
    parser.add_argument(
        "--run",
        action="store_true",
        help="one-shot first-event run; requires --authority-commit",
    )
    parser.add_argument(
        "--authority-commit",
        help="exact lowercase freeze commit; required with --run",
    )
    arguments = parser.parse_args(argv)
    has_commit = arguments.authority_commit is not None
    if arguments.status and arguments.run:
        _print(
            {
                "error": "--status never writes a namespace or constructs source",
                "physical_source_constructed": False,
                "safe_to_run": False,
            }
        )
        return 2
    if arguments.status:
        if not has_commit:
            _print(status())
            return 0
        try:
            _print(committed_image_status(arguments.authority_commit))
        except PRO20EV1RuntimeError as error:
            payload = _cli_error_payload(error)
            payload["physical_source_constructed"] = False
            _print(payload)
            return 1
        return 0
    if arguments.run:
        if not has_commit:
            _print(
                {
                    "error": "run refuses absent authority",
                    "physical_source_constructed": False,
                    "safe_to_run": False,
                    "campaign_execution_authorized": False,
                }
            )
            return 2
        try:
            _print(
                run_first_event(
                    ROOT,
                    authority_commit=arguments.authority_commit,
                )
            )
        except PRO20EV1RuntimeError as error:
            _print(_cli_error_payload(error))
            return 1
        return 0
    if has_commit:
        try:
            _print(committed_image_status(arguments.authority_commit))
        except PRO20EV1RuntimeError as error:
            payload = _cli_error_payload(error)
            payload["physical_source_constructed"] = False
            _print(payload)
            return 1
        return 0
    _print(status())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
