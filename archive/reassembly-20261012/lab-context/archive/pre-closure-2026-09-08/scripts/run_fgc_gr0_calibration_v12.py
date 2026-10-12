#!/usr/bin/env python3
"""Execute one HLT10-authorized fresh PROTO12 GR-0 calibration.

The inherited campaign loop, integrators, transaction, retries, boundary
ledger, temporal spectra, and trapped-sign calculation remain unchanged.
PROTO12 binds only the authorization, frozen member construction, SRC4 source
operator, pairwise common-event classification, and runner identity.
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

from scripts import reproduce_fgc_hlt10_mon10 as hlt10  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_gr0_calibration_v6 as proto6_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v7 as proto7_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v8 as proto8_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v11 as proto11_runner  # noqa: E402
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    make_gr0_center_boundary_projector,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    project_gr0_reference_balanced_state,
)
from recursive_horizons.fgc.evolution.proto12_runtime import (  # noqa: E402
    PROTO12_POINT_COUNTS,
    Proto12GR0EvolutionOperator,
    proto12_gr0_common_event,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    proto5_calibration_trapped_assessment,
)


DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-cal9-run1.toml"
DEFAULT_AUTHORIZATION = REPOSITORY / "results/fgc-1-hlt10-mon10.json"
RUNNER_ID = "FGC-1-CAL9-RUN1-RUNNER"

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
    plan = hlt10.validate_run_plan(plan_path)
    authorization = hlt10.load_canonical_result(authorization_path)
    gates = authorization.get("gate_status", {})
    if (
        authorization.get("artifact_id") != hlt10.ARTIFACT_ID
        or gates.get("PROTO12_successor_runtime_implemented") is not True
        or gates.get("PROTO12_fresh_GR0_dynamic_calibration_authorized")
        is not True
        or gates.get("PROTO12_resolved_holdout_manifest_authorized") is not False
        or gates.get("classical_spherical_diagnostic_authorized") is not False
        or gates.get("SGBL_execution_authorized") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or gates.get("retained_EFT_evolution_authorized") is not False
        or gates.get("physical_transition_claim_authorized") is not False
    ):
        raise ValueError("HLT10 authorization is absent, incomplete, or over-broad")
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if frozen["run_plan_sha256"] != inherited._sha(plan_path):
        raise ValueError("CAL9 run-plan hash differs from HLT10")
    if frozen.get("point_counts") != list(PROTO12_POINT_COUNTS):
        raise ValueError("HLT10 does not authorize the exact PROTO12 ladder")
    if (
        frozen.get("PROTO12_t0_passed_count") != 4
        or frozen.get("SRC4_source_precheck_passed_count") != 12
    ):
        raise ValueError("HLT10 did not admit every frozen t0 premise")
    for relative, expected in authorization["implementation_sha256"].items():
        current = REPOSITORY / relative
        if not current.is_file() or inherited._sha(current) != expected:
            raise ValueError(f"authorized implementation drifted: {relative}")
    if require_fresh_namespace:
        reproduced = hlt10.record(hlt10.DEFAULT_CONFIG)
        if inherited._serial(reproduced) != authorization:
            raise ValueError("HLT10 authorization does not reproduce before launch")
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
        raise ValueError("HLT10 does not contain six inputs for this amplitude")
    return selected


def _legacy_hlt9_authorization(
    authorization: Mapping[str, Any], amplitude: str
) -> dict[str, Any]:
    records = []
    for item in _input_records(authorization, amplitude).values():
        legacy = dict(item)
        legacy["projected_state_sha256"] = item[
            "predecessor_HLT9_projected_state_sha256"
        ]
        legacy["expanded_run_config_sha256"] = item[
            "predecessor_HLT9_expanded_run_config_sha256"
        ]
        legacy["plan_sha256"] = item["predecessor_HLT9_plan_sha256"]
        legacy["protocol_sha256"] = item["predecessor_HLT9_protocol_sha256"]
        legacy["PROTO11_resolution_role"] = item[
            "predecessor_HLT9_resolution_role"
        ]
        records.append(legacy)
    return {"artifact_payload": {"frozen_run_inputs": records}}


def _build_members(
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
    amplitude: str,
) -> dict[str, proto7_runner.Proto7RunMember]:
    """Build fresh PROTO12 members without changing the inherited loop."""

    base = proto11_runner._build_members(
        plan,
        _legacy_hlt9_authorization(authorization, amplitude),
        amplitude,
    )
    frozen = _input_records(authorization, amplitude)
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    thresholds = plan["universal_thresholds"]
    answer: dict[str, proto7_runner.Proto7RunMember] = {}
    for key, member in base.items():
        expected = frozen[(member.method_label, member.point_count)]
        state = project_gr0_reference_balanced_state(
            member.initial,
            spatial_order=member.spatial_order,
        )
        if array_content_sha256(state.u, state.p, state.q) != expected[
            "projected_state_sha256"
        ]:
            raise ValueError("runtime PROTO12 state differs from HLT10 input")
        projected_initial = replace(member.initial, state=state)
        operator = Proto12GR0EvolutionOperator(
            projected_initial.grid,
            spatial_order=member.spatial_order,
            ko_dissipation=float(inherited.Q(numerics["ko_dissipation"])),
            raw_tolerance=float(
                inherited.Q(thresholds["source_residual_infinity_max"])
            ),
            kinetic_condition_maximum=float(
                inherited.Q(thresholds["kinetic_condition_number_max"])
            ),
            maximum_refinement_iterations=thresholds["source_iteration_max"],
            point_batch_size=numerics["source_point_batch_size"],
        )
        initial_rhs = operator(0.0, state)
        if (
            initial_rhs.diagnostics["source_raw_gate_passed"] is not True
            or initial_rhs.diagnostics[
                "SRC4_tensor_contracted_reference_source"
            ]
            is not True
            or initial_rhs.diagnostics["PROTO11_interior_q_reprojected"]
            is not False
            or initial_rhs.diagnostics[
                "PROTO11_exact_reference_equilibrium_applied"
            ]
            is not False
        ):
            raise ValueError("runtime PROTO12 source precheck failed")
        if (
            member.transaction.causal_state.previous_speed_upper
            != initial_rhs.diagnostics["coordinate_speed_upper"]
        ):
            raise RuntimeError("PROTO12 changed the initial characteristic speed")
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
        answer[key] = replace(
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
    if len(answer) != 6:
        raise RuntimeError("PROTO12 GR-0 amplitude requires six run members")
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
        or tuple(item.point_count for item in records) != PROTO12_POINT_COUNTS
    ):
        raise ValueError(f"PROTO12 {method} members differ from the frozen ladder")
    return records


def _event_assessment(
    plan: Mapping[str, Any],
    members: Mapping[str, Any],
    coordinate_time: float,
) -> dict[str, Any]:
    """Retain PROTO11 evidence and replace only pairwise interpretation."""

    base = proto11_runner._event_assessment(plan, members, coordinate_time)
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    primary = _method_members(members, "RK4")
    comparator = _method_members(members, "SSPRK3")

    def compose(records: tuple[Any, ...], method: str):
        return proto12_gr0_common_event(
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
        inherited._serial(primary_common.legacy_PROTO11_common_event)
        != inherited._serial(base["primary_common_event"])
        or inherited._serial(comparator_common.legacy_PROTO11_common_event)
        != inherited._serial(base["comparator_common_event"])
    ):
        raise RuntimeError("PROTO12 changed the inherited PROTO11 common event")
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
            "primary_legacy_PROTO11_common_event": base["primary_common_event"],
            "comparator_legacy_PROTO11_common_event": base[
                "comparator_common_event"
            ],
            "legacy_PROTO11_trapped_assessment": base["trapped_assessment"],
            "primary_common_event": primary_common,
            "comparator_common_event": comparator_common,
            "trapped_assessment": trapped,
            "qualified_trapped_common_event": (
                trapped.trapped_sign_passed and temporal["admission_passed"]
            ),
            "PROTO12_resolution_ladder_enforced": True,
            "PROTO12_SRC4_source_backend_applied": True,
            "PROTO12_pairwise_spectral_classifier_applied": True,
            "PROTO12_all_raw_ratios_retained": True,
            "PROTO11_reference_state_map_applied": True,
            "PROTO11_common_event_constraints_retained": True,
            "PROTO11_interior_q_reprojected": False,
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
    """Run one exact HLT10-authorized campaign through the inherited loop."""

    with _bound_inherited_engine():
        result = proto8_runner.run_campaign(
            plan_path=plan_path,
            authorization_path=authorization_path,
            resume=resume,
        )
    if result.get("runner_id") != RUNNER_ID:
        raise RuntimeError("PROTO12 campaign result lost its runner identity")
    result.update(
        {
            "PROTO12_SRC4_source_backend_enabled": True,
            "PROTO12_pairwise_spectral_classifier_enabled": True,
            "PROTO12_all_raw_ratios_retained": True,
            "PROTO11_reference_state_map_enabled": True,
            "PROTO11_interior_q_reprojection_enabled": False,
        }
    )
    result_path = (
        REPOSITORY
        / hlt10.validate_run_plan(plan_path)["provenance"]["output_root"]
        / "campaign-result.json"
    )
    observed = inherited._load_manifest(result_path)
    for key in (
        "PROTO12_SRC4_source_backend_enabled",
        "PROTO12_pairwise_spectral_classifier_enabled",
        "PROTO12_all_raw_ratios_retained",
        "PROTO11_reference_state_map_enabled",
        "PROTO11_interior_q_reprojection_enabled",
    ):
        if key in observed:
            raise RuntimeError("inherited campaign unexpectedly owns PROTO12 fields")
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
                    "point_counts": list(PROTO12_POINT_COUNTS),
                    "source_backend": "SRC4",
                    "spectral_classifier": "PROTO12_pairwise",
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
