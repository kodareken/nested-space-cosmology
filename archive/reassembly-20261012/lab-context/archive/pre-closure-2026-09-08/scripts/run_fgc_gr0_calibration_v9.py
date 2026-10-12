#!/usr/bin/env python3
"""Execute the HLT7-authorized fresh PROTO9 GR-0 calibration campaign.

The complete evolution loop is the immutable PROTO8 runner.  PROTO9 binds that
engine to a new authorization and enforces the exact ``2049, 4097, 8193``
common-event ladder.  No equation, source solver, threshold, method, retry,
transaction, boundary rule, observable, or candidate branch is changed.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import json
from pathlib import Path
from typing import Any, Iterator, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_hlt7_mon7 as hlt7  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_gr0_calibration_v8 as proto8_runner  # noqa: E402
from recursive_horizons.fgc.evolution.proto9_runtime import (  # noqa: E402
    PROTO9_POINT_COUNTS,
    validate_proto9_common_event_assessment,
)


DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-cal6-run1.toml"
DEFAULT_AUTHORIZATION = REPOSITORY / "results/fgc-1-hlt7-mon7.json"
RUNNER_ID = "FGC-1-CAL6-RUN1-RUNNER"

_INHERITED_VALIDATE_AUTHORIZATION = proto8_runner._validate_authorization
_INHERITED_EVENT_ASSESSMENT = proto8_runner._event_assessment
_INHERITED_RUNNER_ID = proto8_runner.RUNNER_ID


def _validate_authorization(
    plan_path: Path,
    authorization_path: Path,
    *,
    require_fresh_namespace: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = hlt7.validate_run_plan(plan_path)
    authorization = hlt7.load_canonical_result(authorization_path)
    gates = authorization.get("gate_status", {})
    if (
        authorization.get("artifact_id") != hlt7.ARTIFACT_ID
        or gates.get("PROTO9_successor_runtime_implemented") is not True
        or gates.get("PROTO9_fresh_GR0_dynamic_calibration_authorized") is not True
        or gates.get("PROTO9_resolved_holdout_manifest_authorized") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or gates.get("retained_EFT_evolution_authorized") is not False
        or gates.get("physical_transition_claim_authorized") is not False
    ):
        raise ValueError("HLT7 authorization is absent, incomplete, or over-broad")
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if frozen["run_plan_sha256"] != inherited._sha(plan_path):
        raise ValueError("CAL6 run-plan hash differs from HLT7")
    if frozen.get("point_counts") != list(PROTO9_POINT_COUNTS):
        raise ValueError("HLT7 does not authorize the exact PROTO9 ladder")
    for relative, expected in authorization["implementation_sha256"].items():
        current = REPOSITORY / relative
        if not current.is_file() or inherited._sha(current) != expected:
            raise ValueError(f"authorized implementation drifted: {relative}")
    if require_fresh_namespace:
        reproduced = hlt7.record(hlt7.DEFAULT_CONFIG)
        if inherited._serial(reproduced) != authorization:
            raise ValueError("HLT7 authorization does not reproduce before launch")
    return plan, authorization


def _event_assessment(
    plan: Mapping[str, Any],
    members: Mapping[str, Any],
    coordinate_time: float,
) -> dict[str, Any]:
    """Run the inherited compositor and enforce the PROTO9 ladder."""

    assessment = _INHERITED_EVENT_ASSESSMENT(plan, members, coordinate_time)
    validate_proto9_common_event_assessment(assessment["primary_common_event"])
    validate_proto9_common_event_assessment(assessment["comparator_common_event"])
    counts = tuple(
        sorted(
            {
                member.point_count
                for member in members.values()
                if member.method_label == "RK4"
            }
        )
    )
    if counts != PROTO9_POINT_COUNTS:
        raise ValueError("PROTO9 runtime members do not match the frozen ladder")
    answer = dict(assessment)
    answer["PROTO9_resolution_ladder_enforced"] = True
    answer["PROTO9_new_8193_member_present"] = True
    return answer


@contextmanager
def _bound_inherited_engine() -> Iterator[None]:
    """Bind the immutable PROTO8 loop to the explicit PROTO9 contract.

    The predecessor file is hash-bound by HLT7 and remains byte-for-byte
    unchanged.  These three process-local bindings are restored even if the
    campaign raises or is interrupted.
    """

    if (
        proto8_runner._validate_authorization is not _INHERITED_VALIDATE_AUTHORIZATION
        or proto8_runner._event_assessment is not _INHERITED_EVENT_ASSESSMENT
        or proto8_runner.RUNNER_ID != _INHERITED_RUNNER_ID
    ):
        raise RuntimeError("PROTO8 runner bindings were already changed")
    proto8_runner._validate_authorization = _validate_authorization
    proto8_runner._event_assessment = _event_assessment
    proto8_runner.RUNNER_ID = RUNNER_ID
    try:
        yield
    finally:
        proto8_runner._validate_authorization = _INHERITED_VALIDATE_AUTHORIZATION
        proto8_runner._event_assessment = _INHERITED_EVENT_ASSESSMENT
        proto8_runner.RUNNER_ID = _INHERITED_RUNNER_ID


def run_campaign(
    *,
    plan_path: Path,
    authorization_path: Path,
    resume: bool,
) -> dict[str, Any]:
    """Run one authorized campaign through the unchanged predecessor loop."""

    with _bound_inherited_engine():
        result = proto8_runner.run_campaign(
            plan_path=plan_path,
            authorization_path=authorization_path,
            resume=resume,
        )
    if result.get("runner_id") != RUNNER_ID:
        raise RuntimeError("PROTO9 campaign result lost its runner identity")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--authorization", type=Path, default=DEFAULT_AUTHORIZATION)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--check-authorization-only", action="store_true")
    args = parser.parse_args()
    plan_path = args.plan.resolve()
    authorization_path = args.authorization.resolve()
    if args.check_authorization_only:
        inherited._require_clean_tracked_worktree()
        plan, authorization = _validate_authorization(
            plan_path,
            authorization_path,
            require_fresh_namespace=True,
        )
        print(
            json.dumps(
                {
                    "authorized": True,
                    "campaign_id": authorization["artifact_payload"][
                        "frozen_run_plan"
                    ]["campaign_id"],
                    "ordered_amplitudes": plan["candidate_selection"][
                        "ordered_amplitudes"
                    ],
                    "point_counts": list(PROTO9_POINT_COUNTS),
                    "output_created": False,
                },
                sort_keys=True,
            )
        )
        return
    result = run_campaign(
        plan_path=plan_path,
        authorization_path=authorization_path,
        resume=args.resume,
    )
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
