#!/usr/bin/env python3
"""Guarded one-shot RSRC2 common-event runner and read-only status CLI."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
SOURCE = str(ROOT / "src")
if SOURCE not in sys.path:
    sys.path.insert(0, SOURCE)

from recursive_horizons.fgc.evolution.pro20_rsrc2_profile import (  # noqa: E402
    FREEZE_ARTIFACT_ID,
    activate_rsrc2,
)


activate_rsrc2()

from recursive_horizons.fgc.evolution.pro20_rsrc1_runtime import (  # noqa: E402
    PRO20RSRC1RuntimeError,
    implementation_identity,
    plan_status,
    require_authority_delta,
    run_first_event,
)
from recursive_horizons.fgc.evolution.pro20_rsrc1_store import (  # noqa: E402
    PRO20RSRC1CampaignStore,
)


def _print(payload: object) -> None:
    sys.stdout.write(
        json.dumps(payload, sort_keys=True, indent=2, ensure_ascii=True) + "\n"
    )


def committed_image_status(authority_commit: str) -> dict[str, object]:
    implementation = implementation_identity(ROOT)
    receipt = require_authority_delta(
        repository_root=ROOT,
        authority_commit=authority_commit,
        implementation_sha256=str(implementation["implementation_sha256"]),
    )
    payload = dict(receipt)
    payload["safe_to_run"] = True
    payload["physical_source_constructed"] = False
    payload["campaign_execution_authorized"] = False
    payload["calibration_eligible"] = False
    return payload


def terminal_receipt() -> dict[str, object]:
    view = PRO20RSRC1CampaignStore.authenticate_read_only(ROOT)
    last = view.generations[-1]
    return {
        "generation_count": len(view.generations),
        "store_generation": last.store_generation,
        "kind": last.kind,
        "disposition": last.disposition,
        "terminal": last.terminal,
        "first_event_complete": last.first_event_complete,
        "checkpoint_sha256": last.checkpoint_sha256,
        "journal_sha256": last.journal_sha256,
        "unclosed_session": view.unclosed_session,
    }


def _fresh_terminal_receipt() -> dict[str, object]:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    completed = subprocess.run(
        (sys.executable, "-I", "-B", os.fspath(Path(__file__).resolve()), "--authenticate-terminal"),
        cwd=ROOT,
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
        env=environment,
        timeout=1800.0,
    )
    if completed.returncode != 0 or completed.stderr:
        raise RuntimeError("fresh RSRC2 terminal authentication failed")
    value = json.loads(completed.stdout.decode("ascii"))
    if type(value) is not dict:
        raise RuntimeError("fresh RSRC2 terminal receipt is malformed")
    return value


def _cli_error_payload(error: PRO20RSRC1RuntimeError) -> dict[str, object]:
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
    parser = argparse.ArgumentParser(description=FREEZE_ARTIFACT_ID)
    parser.add_argument("--status", action="store_true")
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--authenticate-terminal", action="store_true")
    parser.add_argument("--authority-commit")
    arguments = parser.parse_args(argv)
    selected = sum(
        bool(item)
        for item in (arguments.status, arguments.run, arguments.authenticate_terminal)
    )
    if selected > 1:
        _print({"error": "select exactly one RSRC2 operation", "safe_to_run": False})
        return 2
    if arguments.authenticate_terminal:
        try:
            _print(terminal_receipt())
        except Exception as error:
            _print({"error": str(error)})
            return 1
        return 0
    if arguments.status:
        if arguments.authority_commit is None:
            _print(plan_status())
            return 0
        try:
            _print(committed_image_status(arguments.authority_commit))
        except PRO20RSRC1RuntimeError as error:
            _print(_cli_error_payload(error))
            return 1
        return 0
    if arguments.run:
        if arguments.authority_commit is None:
            _print({"error": "run refuses absent authority", "safe_to_run": False})
            return 2
        try:
            result = dict(
                run_first_event(ROOT, authority_commit=arguments.authority_commit)
            )
            receipt = _fresh_terminal_receipt()
            for name in ("checkpoint_sha256", "journal_sha256", "terminal", "first_event_complete", "disposition"):
                if result[name] != receipt[name]:
                    raise RuntimeError(f"fresh terminal {name} differs")
            result["fresh_terminal_authentication"] = receipt
            _print(result)
        except PRO20RSRC1RuntimeError as error:
            _print(_cli_error_payload(error))
            return 1
        return 0
    _print(plan_status())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
