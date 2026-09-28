#!/usr/bin/env python3
"""Guarded SRCQ1-REC1 CLI. Default is plan/status; no physical source is invoked.

Real recovery qualification requires an exact ``--authority-commit`` and an
absent ``--output-directory``. Authority/delta validation is the coordinator
freeze-B seam, not a boolean flag. This script does not remeasure RK members,
retune a cap, adopt an endpoint, or publish a campaign store.
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

from recursive_horizons.evidence_io import (  # noqa: E402
    canonical_json_bytes,
    read_regular_file,
)
from recursive_horizons.fgc.evolution.hlt17_srcq1_rec1 import (  # noqa: E402
    ATTEMPT2_COMPACT_RELATIVE,
    ATTEMPT2_RAW_LEAF_RELATIVE,
    HLT17SRCQ1REC1AuthoritySeamError,
    HLT17SRCQ1REC1Error,
    implementation_identity,
    observe_environment,
    plan_status,
    qualify_and_publish,
    require_authority_delta,
)
from recursive_horizons.fgc.evolution.hlt17_srcq1 import (  # noqa: E402
    read_static_factory_inputs,
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


def committed_image_status(authority_commit: str) -> dict[str, object]:
    """Read-only freeze-B preflight. Never builds origin, runtime, or source."""

    implementation = implementation_identity(ROOT)
    receipt = require_authority_delta(
        repository_root=ROOT,
        authority_commit=authority_commit,
        implementation_sha256=str(implementation["sha256"]),
    )
    payload = dict(receipt)
    payload["safe_to_run"] = True
    payload["physical_source_qualification_executed"] = False
    payload["physical_source_qualification_attempted"] = False
    payload["rk_members_remeasured"] = False
    return payload


def _cli_error_payload(error: HLT17SRCQ1REC1Error) -> dict[str, object]:
    attempted = bool(error.physical_work_started)
    payload: dict[str, object] = {
        "error": str(error),
        "outcome": error.outcome,
        "published": error.published,
        "phase": error.phase,
        "physical_source_qualification_attempted": attempted,
        "rk_members_remeasured": False,
    }
    if error.published and error.outcome == "inconclusive":
        payload["publication_uncertain"] = True
    if not attempted:
        payload["physical_source_qualification_executed"] = False
    return payload


def qualify(*, authority_commit: str, output_directory: Path) -> dict[str, object]:
    implementation = implementation_identity(ROOT)
    authority = require_authority_delta(
        repository_root=ROOT,
        authority_commit=authority_commit,
        implementation_sha256=str(implementation["sha256"]),
    )
    closure = str(authority["source_closure_sha256"])
    environment = authority.get("environment")
    expected_implementation = authority.get("implementation_identity")
    source_configurations = authority.get("source_configuration")
    if (
        not isinstance(environment, dict)
        or not isinstance(expected_implementation, dict)
        or not isinstance(source_configurations, dict)
    ):
        raise HLT17SRCQ1REC1AuthoritySeamError(
            "authority receipt omits pinned implementation/source-config/environment"
        )
    origin_capture = authority.get("origin_capture_sha256")
    expected_origin = origin_capture if type(origin_capture) is str else None
    return qualify_and_publish(
        None,
        repository_root=ROOT,
        output_directory=output_directory,
        authority_commit=authority_commit,
        static_input_bytes=None,
        source_closure_sha256=closure,
        environment=environment,
        attempt2_compact_bytes=None,
        attempt2_raw_bytes=None,
        authority_delta_validator=lambda **kwargs: authority,
        origin_loader=lambda: capture_pro20_historical_origin(ROOT),
        static_input_loader=lambda: read_static_factory_inputs(ROOT),
        environment_loader=observe_environment,
        attempt2_compact_loader=lambda: read_regular_file(ROOT, ATTEMPT2_COMPACT_RELATIVE),
        attempt2_raw_loader=lambda: read_regular_file(ROOT, ATTEMPT2_RAW_LEAF_RELATIVE),
        expected_origin_sha256=expected_origin,
        expected_environment=environment,
        expected_implementation=expected_implementation,
        expected_source_configurations=source_configurations,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="FGC-1-HLT17-SRCQ1-REC1 plan/status or one guarded recovery"
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="read-only plan or committed-image preflight; never builds source",
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
    if arguments.status:
        if has_output:
            _print(
                {
                    "error": "--status never writes an output namespace or builds source",
                    "physical_source_qualification_executed": False,
                    "rk_members_remeasured": False,
                    "safe_to_run": False,
                }
            )
            return 2
        if not has_commit:
            _print(status())
            return 0
        try:
            _print(committed_image_status(arguments.authority_commit))
        except HLT17SRCQ1REC1Error as error:
            payload = _cli_error_payload(error)
            payload["safe_to_run"] = False
            payload["physical_source_qualification_executed"] = False
            _print(payload)
            return 1
        return 0
    if not has_commit and not has_output:
        _print(status())
        return 0
    if not has_commit or not has_output:
        _print(
            {
                "error": "real recovery requires both --authority-commit and --output-directory",
                "physical_source_qualification_executed": False,
                "rk_members_remeasured": False,
            }
        )
        return 2
    try:
        result = qualify(
            authority_commit=arguments.authority_commit,
            output_directory=Path(arguments.output_directory),
        )
    except HLT17SRCQ1REC1Error as error:
        _print(_cli_error_payload(error))
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
