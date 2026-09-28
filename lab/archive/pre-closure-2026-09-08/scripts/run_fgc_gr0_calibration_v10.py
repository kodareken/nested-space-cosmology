#!/usr/bin/env python3
"""Execute the HLT8-authorized fresh PROTO10 GR-0 calibration campaign.

The complete evolution, source, transaction, retry, boundary, temporal, and
trapped-sign loop is the immutable PROTO8 engine.  PROTO10 replaces only the
committed spatial-spectral interpretation: every event retains the raw PROTO9
assessment and freshly recomputes the complete conditioning guard set.
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

from scripts import reproduce_fgc_hlt8_mon8 as hlt8  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_gr0_calibration_v8 as proto8_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v9 as proto9_runner  # noqa: E402
from recursive_horizons.fgc.evolution.proto10_runtime import (  # noqa: E402
    proto10_gr0_common_event,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    proto5_calibration_trapped_assessment,
)
from recursive_horizons.fgc.evolution.proto9_runtime import (  # noqa: E402
    PROTO9_POINT_COUNTS,
)


DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-cal7-run1.toml"
DEFAULT_AUTHORIZATION = REPOSITORY / "results/fgc-1-hlt8-mon8.json"
RUNNER_ID = "FGC-1-CAL7-RUN1-RUNNER"

_INHERITED_VALIDATE_AUTHORIZATION = proto8_runner._validate_authorization
_INHERITED_EVENT_ASSESSMENT = proto9_runner._event_assessment
_INHERITED_RUNNER_ID = proto8_runner.RUNNER_ID


def _validate_authorization(
    plan_path: Path,
    authorization_path: Path,
    *,
    require_fresh_namespace: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = hlt8.validate_run_plan(plan_path)
    authorization = hlt8.load_canonical_result(authorization_path)
    gates = authorization.get("gate_status", {})
    if (
        authorization.get("artifact_id") != hlt8.ARTIFACT_ID
        or gates.get("PROTO10_successor_runtime_implemented") is not True
        or gates.get("PROTO10_fresh_GR0_dynamic_calibration_authorized") is not True
        or gates.get("PROTO10_resolved_holdout_manifest_authorized") is not False
        or gates.get("classical_spherical_diagnostic_authorized") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or gates.get("retained_EFT_evolution_authorized") is not False
        or gates.get("physical_transition_claim_authorized") is not False
    ):
        raise ValueError("HLT8 authorization is absent, incomplete, or over-broad")
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if frozen["run_plan_sha256"] != inherited._sha(plan_path):
        raise ValueError("CAL7 run-plan hash differs from HLT8")
    if frozen.get("point_counts") != list(PROTO9_POINT_COUNTS):
        raise ValueError("HLT8 does not authorize the exact PROTO10 ladder")
    if frozen.get("guarded_PROTO10_t0_passed_count") != 4:
        raise ValueError("HLT8 did not admit all four frozen t0 cases")
    for relative, expected in authorization["implementation_sha256"].items():
        current = REPOSITORY / relative
        if not current.is_file() or inherited._sha(current) != expected:
            raise ValueError(f"authorized implementation drifted: {relative}")
    if require_fresh_namespace:
        reproduced = hlt8.record(hlt8.DEFAULT_CONFIG)
        if inherited._serial(reproduced) != authorization:
            raise ValueError("HLT8 authorization does not reproduce before launch")
    return plan, authorization


def _method_members(members: Mapping[str, Any], method: str) -> tuple[Any, ...]:
    records = tuple(
        sorted(
            (item for item in members.values() if item.method_label == method),
            key=lambda item: item.point_count,
        )
    )
    if len(records) != 3 or tuple(item.point_count for item in records) != PROTO9_POINT_COUNTS:
        raise ValueError(f"PROTO10 {method} members do not match the frozen ladder")
    return records


def _event_assessment(
    plan: Mapping[str, Any],
    members: Mapping[str, Any],
    coordinate_time: float,
) -> dict[str, Any]:
    """Retain raw PROTO9 evidence and freshly apply every PROTO10 guard."""

    base = _INHERITED_EVENT_ASSESSMENT(plan, members, coordinate_time)
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    primary = _method_members(members, "RK4")
    comparator = _method_members(members, "SSPRK3")

    def compose(records: tuple[Any, ...], method: str):
        return proto10_gr0_common_event(
            [item.state for item in records],
            [item.initial.grid for item in records],
            accepted_stage_counts=[
                item.transaction.state.accepted_stage_count for item in records
            ],
            method=method,
            coordinate_time=coordinate_time,
            cutoff=float(inherited.Q(physical["cutoff_Lambda"])),
            measurement_radius_maximum=float(
                inherited.Q(physical["measurement_radius_maximum"])
            ),
            taper_fraction=float(
                inherited.Q(numerics["proper_spectral_taper_fraction"])
            ),
            fixed_outer_rows=numerics["fixed_outer_rows"],
        )

    primary_common = compose(primary, "RK4")
    comparator_common = compose(comparator, "SSPRK3")
    if (
        inherited._serial(primary_common.raw_PROTO9_common_event)
        != inherited._serial(base["primary_common_event"])
        or inherited._serial(comparator_common.raw_PROTO9_common_event)
        != inherited._serial(base["comparator_common_event"])
    ):
        raise RuntimeError("PROTO10 runtime changed the raw PROTO9 assessment")

    trapped = proto5_calibration_trapped_assessment(
        primary_states=[item.state for item in primary],
        comparator_states=[item.state for item in comparator],
        grids=[item.initial.grid for item in primary],
        primary_common_event=primary_common,
        comparator_common_event=comparator_common,
        measurement_radius_maximum=float(
            inherited.Q(physical["measurement_radius_maximum"])
        ),
        minimum_observed_order=float(
            inherited.Q(numerics["minimum_constraint_finest_pair_order"])
        ),
        positive_margin_factor=float(
            plan["candidate_selection"][
                "trapped_sign_margin_over_combined_error_factor"
            ]
        ),
    )
    temporal = base["temporal_spectral_admission"]
    answer = dict(base)
    answer.update(
        {
            "primary_raw_PROTO9_common_event": base["primary_common_event"],
            "comparator_raw_PROTO9_common_event": base[
                "comparator_common_event"
            ],
            "raw_PROTO9_trapped_assessment": base["trapped_assessment"],
            "primary_common_event": primary_common,
            "comparator_common_event": comparator_common,
            "trapped_assessment": trapped,
            "qualified_trapped_common_event": (
                trapped.trapped_sign_passed and temporal["admission_passed"]
            ),
            "PROTO10_resolution_ladder_enforced": True,
            "PROTO10_raw_PROTO9_evidence_retained": True,
            "PROTO10_all_conditioning_guards_recomputed_on_fresh_event": True,
            "PROTO10_round_trip_is_a_continuum_error_bound": False,
            "PROTO10_diagnostic_saturation_is_physical_resolution": False,
        }
    )
    return answer


@contextmanager
def _bound_inherited_engine() -> Iterator[None]:
    """Bind and always restore the immutable PROTO8 campaign loop."""

    if (
        proto8_runner._validate_authorization is not _INHERITED_VALIDATE_AUTHORIZATION
        or proto8_runner._event_assessment is not proto9_runner._INHERITED_EVENT_ASSESSMENT
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
        proto8_runner._event_assessment = proto9_runner._INHERITED_EVENT_ASSESSMENT
        proto8_runner.RUNNER_ID = _INHERITED_RUNNER_ID


def run_campaign(
    *,
    plan_path: Path,
    authorization_path: Path,
    resume: bool,
) -> dict[str, Any]:
    """Run one HLT8-authorized campaign through the unchanged base loop."""

    with _bound_inherited_engine():
        result = proto8_runner.run_campaign(
            plan_path=plan_path,
            authorization_path=authorization_path,
            resume=resume,
        )
    if result.get("runner_id") != RUNNER_ID:
        raise RuntimeError("PROTO10 campaign result lost its runner identity")
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
