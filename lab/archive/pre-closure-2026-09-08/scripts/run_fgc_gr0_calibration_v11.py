#!/usr/bin/env python3
"""Execute one HLT9-authorized fresh PROTO11 GR-0 calibration.

The campaign loop, transaction, retries, time integrators, boundary ledger,
temporal spectra, trapped-sign test, and PROTO10 spatial-spectral rules remain
the immutable predecessor machinery.  PROTO11 binds only four process-local
surfaces: authorization, initial member construction, committed common-event
constraints, and runner identity.  The independently evolved interior ``q``
field is never reprojected after initialization.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import replace
import json
from pathlib import Path
from typing import Any, Iterator, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_hlt9_mon9 as hlt9  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_gr0_calibration_v6 as proto6_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v7 as proto7_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v8 as proto8_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v10 as proto10_runner  # noqa: E402
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    make_gr0_center_boundary_projector,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    SBPFirstDerivative,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    Proto11GR0EvolutionOperator,
    project_gr0_reference_balanced_state,
    proto11_gr0_common_event,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    proto5_calibration_trapped_assessment,
)
from recursive_horizons.fgc.evolution.proto9_runtime import (  # noqa: E402
    PROTO9_POINT_COUNTS,
)


DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-cal8-run1.toml"
DEFAULT_AUTHORIZATION = REPOSITORY / "results/fgc-1-hlt9-mon9.json"
RUNNER_ID = "FGC-1-CAL8-RUN1-RUNNER"

_INHERITED_VALIDATE_AUTHORIZATION = proto8_runner._validate_authorization
_INHERITED_EVENT_ASSESSMENT = proto8_runner._event_assessment
_INHERITED_BUILD_MEMBERS = proto8_runner._build_members
_INHERITED_PROTO6_BUILD_MEMBERS = proto6_runner._build_members
_INHERITED_RUNNER_ID = proto8_runner.RUNNER_ID


def _validate_authorization(
    plan_path: Path,
    authorization_path: Path,
    *,
    require_fresh_namespace: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = hlt9.validate_run_plan(plan_path)
    authorization = hlt9.load_canonical_result(authorization_path)
    gates = authorization.get("gate_status", {})
    if (
        authorization.get("artifact_id") != hlt9.ARTIFACT_ID
        or gates.get("PROTO11_successor_runtime_implemented") is not True
        or gates.get("PROTO11_fresh_GR0_dynamic_calibration_authorized")
        is not True
        or gates.get("PROTO11_resolved_holdout_manifest_authorized") is not False
        or gates.get("classical_spherical_diagnostic_authorized") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or gates.get("retained_EFT_evolution_authorized") is not False
        or gates.get("physical_transition_claim_authorized") is not False
    ):
        raise ValueError("HLT9 authorization is absent, incomplete, or over-broad")
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if frozen["run_plan_sha256"] != inherited._sha(plan_path):
        raise ValueError("CAL8 run-plan hash differs from HLT9")
    if frozen.get("point_counts") != list(PROTO9_POINT_COUNTS):
        raise ValueError("HLT9 does not authorize the exact PROTO11 ladder")
    if frozen.get("PROTO11_t0_passed_count") != 4:
        raise ValueError("HLT9 did not admit all four frozen t0 cases")
    for relative, expected in authorization["implementation_sha256"].items():
        current = REPOSITORY / relative
        if not current.is_file() or inherited._sha(current) != expected:
            raise ValueError(f"authorized implementation drifted: {relative}")
    if require_fresh_namespace:
        reproduced = hlt9.record(hlt9.DEFAULT_CONFIG)
        if inherited._serial(reproduced) != authorization:
            raise ValueError("HLT9 authorization does not reproduce before launch")
    return plan, authorization


def _input_records(
    authorization: Mapping[str, Any], amplitude: str
) -> dict[tuple[str, int], Mapping[str, Any]]:
    records = authorization["artifact_payload"]["frozen_run_inputs"]
    selected = {
        (item["method"], item["point_count"]): item
        for item in records
        if item["amplitude"] == amplitude
    }
    if len(selected) != 6:
        raise ValueError("HLT9 does not contain six inputs for this amplitude")
    return selected


def _legacy_authorization(
    authorization: Mapping[str, Any], amplitude: str
) -> dict[str, Any]:
    records = []
    for item in _input_records(authorization, amplitude).values():
        legacy = dict(item)
        legacy["projected_state_sha256"] = item[
            "predecessor_HLT8_projected_state_sha256"
        ]
        legacy["expanded_run_config_sha256"] = item[
            "predecessor_HLT8_expanded_run_config_sha256"
        ]
        records.append(legacy)
    return {"artifact_payload": {"frozen_run_inputs": records}}


def _build_members(
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
    amplitude: str,
) -> dict[str, proto7_runner.Proto7RunMember]:
    """Build fresh PROTO11 members without changing the inherited loop."""

    base = _INHERITED_PROTO6_BUILD_MEMBERS(
        plan,
        _legacy_authorization(authorization, amplitude),
        amplitude,
    )
    frozen = _input_records(authorization, amplitude)
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    thresholds = plan["universal_thresholds"]
    answer: dict[str, proto7_runner.Proto7RunMember] = {}
    for key, proto6_member in base.items():
        member = proto7_runner._as_proto7_member(proto6_member)
        expected = frozen[(member.method_label, member.point_count)]
        state = project_gr0_reference_balanced_state(
            member.initial,
            spatial_order=member.spatial_order,
        )
        if array_content_sha256(state.u, state.p, state.q) != expected[
            "projected_state_sha256"
        ]:
            raise ValueError("runtime PROTO11 state differs from HLT9 input")
        projected_initial = replace(member.initial, state=state)
        operator = Proto11GR0EvolutionOperator(
            projected_initial.grid,
            spatial_order=member.spatial_order,
            ko_dissipation=float(inherited.Q(numerics["ko_dissipation"])),
            raw_tolerance=float(
                inherited.Q(thresholds["source_residual_infinity_max"])
            ),
            kinetic_condition_maximum=float(
                inherited.Q(thresholds["kinetic_condition_number_max"])
            ),
        )
        initial_rhs = operator(0.0, state)
        if initial_rhs.diagnostics["source_raw_gate_passed"] is not True:
            raise ValueError("runtime PROTO11 source precheck failed")
        if (
            initial_rhs.diagnostics["PROTO11_interior_q_reprojected"] is not False
            or initial_rhs.diagnostics[
                "PROTO11_exact_reference_equilibrium_applied"
            ]
            is not False
        ):
            raise ValueError("nontrivial PROTO11 input used a forbidden map path")
        if (
            member.transaction.causal_state.previous_speed_upper
            != initial_rhs.diagnostics["coordinate_speed_upper"]
        ):
            raise RuntimeError("PROTO11 changed the initial characteristic speed")
        tracers = inherited.NormalFlowTracers.create(
            minimum=float(
                inherited.Q(numerics["normal_flow_tracer_radius_minimum"])
            ),
            maximum=float(
                inherited.Q(numerics["normal_flow_tracer_radius_maximum"])
            ),
            spacing=float(inherited.Q(numerics["normal_flow_tracer_spacing"])),
            state=state,
            coordinates=projected_initial.grid.coordinates,
            cutoff=float(inherited.Q(physical["cutoff_Lambda"])),
            outer_radius=float(inherited.Q(physical["outer_radius"])),
        )
        rebuilt = replace(
            member,
            input_hash=expected["expanded_run_config_sha256"],
            initial=projected_initial,
            state=state,
            operator=operator,
            projector=make_gr0_center_boundary_projector(
                projected_initial,
                fixed_outer_rows=numerics["fixed_outer_rows"],
            ),
            tracers=tracers,
        )
        answer[key] = rebuilt
    if len(answer) != 6:
        raise RuntimeError("PROTO11 GR-0 amplitude requires six run members")
    return answer


def _method_members(
    members: Mapping[str, Any], method: str
) -> tuple[Any, ...]:
    records = tuple(
        sorted(
            (item for item in members.values() if item.method_label == method),
            key=lambda item: item.point_count,
        )
    )
    if (
        len(records) != 3
        or tuple(item.point_count for item in records) != PROTO9_POINT_COUNTS
    ):
        raise ValueError(f"PROTO11 {method} members differ from the frozen ladder")
    return records


def _event_assessment(
    plan: Mapping[str, Any],
    members: Mapping[str, Any],
    coordinate_time: float,
) -> dict[str, Any]:
    """Retain PROTO10 spectra and replace only the operative constraints."""

    base = proto10_runner._event_assessment(plan, members, coordinate_time)
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    primary = _method_members(members, "RK4")
    comparator = _method_members(members, "SSPRK3")

    def compose(records: tuple[Any, ...], method: str):
        return proto11_gr0_common_event(
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
        inherited._serial(primary_common.legacy_PROTO10_common_event)
        != inherited._serial(base["primary_common_event"])
        or inherited._serial(comparator_common.legacy_PROTO10_common_event)
        != inherited._serial(base["comparator_common_event"])
    ):
        raise RuntimeError("PROTO11 changed the inherited PROTO10 spectral event")
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
            "primary_raw_PROTO10_common_event": base["primary_common_event"],
            "comparator_raw_PROTO10_common_event": base[
                "comparator_common_event"
            ],
            "raw_PROTO10_trapped_assessment": base["trapped_assessment"],
            "primary_common_event": primary_common,
            "comparator_common_event": comparator_common,
            "trapped_assessment": trapped,
            "qualified_trapped_common_event": (
                trapped.trapped_sign_passed and temporal["admission_passed"]
            ),
            "PROTO11_resolution_ladder_enforced": True,
            "PROTO11_reference_state_map_applied": True,
            "PROTO11_common_event_constraints_use_reference_map": True,
            "PROTO11_interior_q_reprojected": False,
            "PROTO11_complete_PROTO10_spectral_evidence_retained": True,
        }
    )
    return answer


@contextmanager
def _bound_inherited_engine() -> Iterator[None]:
    """Bind and always restore the immutable campaign machinery."""

    if (
        proto8_runner._validate_authorization is not _INHERITED_VALIDATE_AUTHORIZATION
        or proto8_runner._event_assessment is not _INHERITED_EVENT_ASSESSMENT
        or proto8_runner._build_members is not _INHERITED_BUILD_MEMBERS
        or proto6_runner._build_members is not _INHERITED_PROTO6_BUILD_MEMBERS
        or proto8_runner.RUNNER_ID != _INHERITED_RUNNER_ID
    ):
        raise RuntimeError("inherited runner bindings were already changed")
    proto8_runner._validate_authorization = _validate_authorization
    proto8_runner._event_assessment = _event_assessment
    proto8_runner._build_members = _build_members
    proto6_runner._build_members = _build_members
    proto8_runner.RUNNER_ID = RUNNER_ID
    try:
        yield
    finally:
        proto8_runner._validate_authorization = _INHERITED_VALIDATE_AUTHORIZATION
        proto8_runner._event_assessment = _INHERITED_EVENT_ASSESSMENT
        proto8_runner._build_members = _INHERITED_BUILD_MEMBERS
        proto6_runner._build_members = _INHERITED_PROTO6_BUILD_MEMBERS
        proto8_runner.RUNNER_ID = _INHERITED_RUNNER_ID


def run_campaign(
    *,
    plan_path: Path,
    authorization_path: Path,
    resume: bool,
) -> dict[str, Any]:
    """Run one exact HLT9-authorized campaign through the inherited loop."""

    with _bound_inherited_engine():
        result = proto8_runner.run_campaign(
            plan_path=plan_path,
            authorization_path=authorization_path,
            resume=resume,
        )
    if result.get("runner_id") != RUNNER_ID:
        raise RuntimeError("PROTO11 campaign result lost its runner identity")
    result.update(
        {
            "PROTO11_reference_state_map_enabled": True,
            "PROTO11_interior_q_reprojection_enabled": False,
            "PROTO11_complete_PROTO10_spectral_contract_retained": True,
        }
    )
    result_path = (
        REPOSITORY
        / hlt9.validate_run_plan(plan_path)["provenance"]["output_root"]
        / "campaign-result.json"
    )
    observed = inherited._load_manifest(result_path)
    for key in (
        "PROTO11_reference_state_map_enabled",
        "PROTO11_interior_q_reprojection_enabled",
        "PROTO11_complete_PROTO10_spectral_contract_retained",
    ):
        if key in observed:
            raise RuntimeError("inherited campaign unexpectedly owns PROTO11 result fields")
    inherited._atomic_replace(result_path, inherited._canonical(result))
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
                    "reference_state_map": "PROTO11",
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
