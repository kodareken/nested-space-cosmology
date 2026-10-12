#!/usr/bin/env python3
"""Guarded SRCQ1 CLI. Default is plan/status; no physical source is invoked.

Real qualification requires an exact ``--authority-commit`` and an absent
``--output-directory``. Authority/delta validation is the coordinator freeze-B
seam, not a boolean flag. This script does not retune a cap, adopt an
endpoint, or publish a campaign store.
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

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.evolution.hlt17_srcq1 import (  # noqa: E402
    HLT17SRCQ1Error,
    implementation_identity,
    observe_environment,
    plan_status,
    qualify_and_publish,
    read_static_factory_inputs,
    require_authority_delta,
)
from recursive_horizons.fgc.evolution.pro20_origin import (  # noqa: E402
    capture_pro20_historical_origin,
)


def _print(payload: object) -> None:
    sys.stdout.write(
        json.dumps(payload, sort_keys=True, indent=2, ensure_ascii=True) + "\n"
    )


def status() -> dict[str, object]:
    return plan_status()


def qualify(*, authority_commit: str, output_directory: Path) -> dict[str, object]:
    implementation = implementation_identity(ROOT)
    authority = require_authority_delta(
        repository_root=ROOT,
        authority_commit=authority_commit,
        implementation_sha256=str(implementation["sha256"]),
    )
    origin = capture_pro20_historical_origin(ROOT)
    static_inputs = read_static_factory_inputs(ROOT)
    closure = authority.get("source_closure_sha256", "0" * 64)
    if type(closure) is not str or len(closure) != 64 or closure.lower() != closure:
        raise HLT17SRCQ1Error("authority receipt source_closure_sha256 differs")
    return qualify_and_publish(
        origin,
        repository_root=ROOT,
        output_directory=output_directory,
        authority_commit=authority_commit,
        static_input_bytes=static_inputs,
        source_closure_sha256=closure,
        environment=observe_environment(),
        expected_origin_sha256=origin.capture_sha256,
        authority_delta_validator=lambda **kwargs: authority,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="FGC-1-HLT17-SRCQ1 plan/status or one guarded qualification"
    )
    parser.add_argument(
        "--authority-commit",
        help="exact lowercase freeze-B commit; required with --output-directory",
    )
    parser.add_argument(
        "--output-directory",
        help="absent namespace for the atomic qualification record",
    )
    arguments = parser.parse_args(argv)
    has_commit = arguments.authority_commit is not None
    has_output = arguments.output_directory is not None
    if not has_commit and not has_output:
        _print(status())
        return 0
    if not has_commit or not has_output:
        _print(
            {
                "error": "real qualification requires both --authority-commit and --output-directory",
                "physical_six_member_qualification_executed": False,
            }
        )
        return 2
    try:
        result = qualify(
            authority_commit=arguments.authority_commit,
            output_directory=Path(arguments.output_directory),
        )
    except HLT17SRCQ1Error as error:
        payload = {
            "error": str(error),
            "outcome": error.outcome,
            "published": error.published,
            "physical_six_member_qualification_executed": False,
        }
        _print(payload)
        return 1
    _print(
        {
            "published": True,
            "payload_sha256": result["publication"]["payload_sha256"],
            "output": str(arguments.output_directory),
        }
    )
    canonical_json_bytes(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
