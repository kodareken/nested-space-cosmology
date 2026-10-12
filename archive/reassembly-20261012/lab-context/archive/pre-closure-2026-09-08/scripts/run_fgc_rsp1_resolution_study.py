#!/usr/bin/env python3
"""Run the frozen RSP1 amplitude-three resolution-spectrum study.

The runner advances one GR-0 amplitude, two independently frozen methods, and
the exact ``4097 -> 8193 -> 16385`` ladder to the single synchronized endpoint
``t=1/16``.  It uses the existing transaction, rollback, source-retry, CFL,
boundary, and checkpoint machinery.  It does not select a calibration
amplitude, inspect a candidate branch, or classify collapse.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import time as wall_time
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_rsp1_frz1 as freeze  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_gr0_calibration_v6 as proto6_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v7 as proto7_runner  # noqa: E402
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
    make_gr0_center_boundary_projector,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    project_gr0_reference_balanced_state,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeMonitorState,
)
from recursive_horizons.fgc.evolution.proto6_runtime import (  # noqa: E402
    restore_common_event_snapshots,
)
from recursive_horizons.fgc.evolution.proto7_runtime import (  # noqa: E402
    Proto7TerminalStop,
)
from recursive_horizons.fgc.evolution.rsp1_resolution_runtime import (  # noqa: E402
    RSP1GR0EvolutionOperator,
    RSP1_POINT_COUNTS,
    rsp1_gr0_common_event,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = inherited.Q
DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-rsp1-run1.toml"
DEFAULT_AUTHORIZATION = REPOSITORY / "results/fgc-1-rsp1-frz1.json"
RUNNER_ID = "FGC-1-RSP1-RUN1-RUNNER"


def _validate_authorization(
    plan_path: Path,
    authorization_path: Path,
    *,
    require_fresh_namespace: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = freeze.validate_run_plan(plan_path)
    authorization = freeze.load_canonical_result(authorization_path)
    historical_namespace = authorization.get("artifact_payload", {}).get(
        "namespace_precondition"
    )
    if inherited._serial(
        freeze.reproduce(
            freeze.DEFAULT_CONFIG,
            historical_namespace_evidence=historical_namespace,
        )
    ) != authorization:
        raise ValueError("RSP1 authorization does not reproduce")
    gates = authorization.get("gate_status", {})
    if (
        authorization.get("artifact_id") != freeze.ARTIFACT_ID
        or gates.get("RSP1_runtime_implemented") is not True
        or gates.get("RSP1_execution_authorized") is not True
        or gates.get("amplitude_three_spectral_veto_cleared") is not False
        or gates.get("PROTO12_frozen") is not False
        or gates.get("classical_spherical_diagnostic_authorized") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or gates.get("retained_EFT_evolution_authorized") is not False
        or gates.get("physical_transition_claim_authorized") is not False
    ):
        raise ValueError("RSP1 authorization is absent, incomplete, or over-broad")
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if frozen["run_plan_sha256"] != inherited._sha(plan_path):
        raise ValueError("RSP1 run-plan hash differs from the authorization")
    if frozen["point_counts"] != list(RSP1_POINT_COUNTS):
        raise ValueError("RSP1 authorization has another resolution ladder")
    for relative, expected in authorization["implementation_sha256"].items():
        current = REPOSITORY / relative
        if not current.is_file() or inherited._sha(current) != expected:
            raise ValueError(f"authorized RSP1 implementation drifted: {relative}")
    if require_fresh_namespace:
        output = REPOSITORY / plan["provenance"]["output_root"]
        if output.exists():
            raise ValueError("fresh RSP1 output root already exists")
    return plan, authorization


def _input_records(
    authorization: Mapping[str, Any],
) -> dict[tuple[str, int], Mapping[str, Any]]:
    records = authorization["artifact_payload"]["frozen_run_inputs"]
    selected = {
        (item["method"], item["point_count"]): item for item in records
    }
    expected = {
        (method, point_count)
        for method in ("RK4", "SSPRK3")
        for point_count in RSP1_POINT_COUNTS
    }
    if set(selected) != expected:
        raise ValueError("RSP1 authorization does not contain exactly six inputs")
    return selected


def _build_members(
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
) -> dict[str, proto7_runner.Proto7RunMember]:
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    thresholds = plan["universal_thresholds"]
    frozen = _input_records(authorization)
    parameters = PulseParameters(
        chi_amplitude=float(Q(physical["chi_amplitude"])),
        center=float(Q(physical["chi_center"])),
        half_width=float(Q(physical["chi_half_width"])),
        phi_amplitude=float(Q(physical["phi_seed_amplitude"])),
        planck_mass=float(Q(physical["planck_mass"])),
        scalar_mass=float(Q(physical["scalar_mass"])),
        quartic_coupling=float(Q(physical["quartic_coupling"])),
    )
    answer: dict[str, proto7_runner.Proto7RunMember] = {}
    for method in (plan["primary_method"], plan["comparator_method"]):
        method_label = method["method_label"]
        order = method["spatial_order"]
        for point_count in numerics["resolutions"]:
            expected = frozen[(method_label, point_count)]
            initial = construct_gr0_grid_initial_data(
                parameters,
                point_count=point_count,
                outer_radius=float(Q(physical["outer_radius"])),
                constraint_method=method["constraint_solve_method"],
                diagnostic_spatial_order=order,
            )
            state = project_gr0_reference_balanced_state(
                initial,
                spatial_order=order,
            )
            if array_content_sha256(state.u, state.p, state.q) != expected[
                "projected_state_sha256"
            ]:
                raise ValueError("runtime RSP1 state differs from the frozen input")
            projected_initial = replace(initial, state=state)
            operator = RSP1GR0EvolutionOperator(
                projected_initial.grid,
                spatial_order=order,
                ko_dissipation=float(Q(numerics["ko_dissipation"])),
                raw_tolerance=float(
                    Q(thresholds["source_residual_infinity_max"])
                ),
                kinetic_condition_maximum=float(
                    Q(thresholds["kinetic_condition_number_max"])
                ),
                point_batch_size=numerics["source_point_batch_size"],
            )
            initial_rhs = operator(0.0, state)
            if initial_rhs.diagnostics["source_raw_gate_passed"] is not True:
                raise ValueError("runtime RSP1 initial source precheck failed")
            derivative = SBPFirstDerivative(initial.grid, order)
            transaction = proto6_runner.GR0RuntimeStageTransaction(
                thresholds=proto6_runner.GR0UniversalThresholds(
                    source_residual_maximum=float(
                        Q(thresholds["source_residual_infinity_max"])
                    ),
                    source_iteration_maximum=thresholds["source_iteration_max"],
                    kinetic_condition_maximum=float(
                        Q(thresholds["kinetic_condition_number_max"])
                    ),
                ),
                causal_state=CausalBudgetState(
                    previous_speed_upper=float(
                        initial_rhs.diagnostics["coordinate_speed_upper"]
                    )
                ),
                boundary_geometry=BoundaryGeometry(
                    float(Q(physical["outer_radius"])),
                    float(Q(physical["measurement_radius_maximum"])),
                    float(Q(thresholds["minimum_boundary_causal_buffer"])),
                    derivative.stencil_reach_intervals * initial.grid.spacing,
                ),
                grid_spacing=initial.grid.spacing,
                cfl_maximum=float(Q(numerics["cfl_maximum"])),
                hat_normal_factor=float(Q(thresholds["hat_normal_factor"])),
            )
            tracers = inherited.NormalFlowTracers.create(
                minimum=float(
                    Q(numerics["normal_flow_tracer_radius_minimum"])
                ),
                maximum=float(
                    Q(numerics["normal_flow_tracer_radius_maximum"])
                ),
                spacing=float(Q(numerics["normal_flow_tracer_spacing"])),
                state=state,
                coordinates=projected_initial.grid.coordinates,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                outer_radius=float(Q(physical["outer_radius"])),
            )
            member = proto7_runner.Proto7RunMember(
                amplitude="3",
                method_label=method_label,
                integrator_id=method["integrator_id"],
                spatial_order=order,
                point_count=point_count,
                input_hash=expected["expanded_run_config_sha256"],
                initial=projected_initial,
                state=state,
                operator=operator,
                projector=make_gr0_center_boundary_projector(
                    projected_initial,
                    fixed_outer_rows=numerics["fixed_outer_rows"],
                ),
                transaction=transaction,
                tracers=tracers,
            )
            if member.key in answer:
                raise RuntimeError("duplicate RSP1 run member")
            answer[member.key] = member
    if len(answer) != 6:
        raise RuntimeError("RSP1 requires exactly six run members")
    return answer


def _method_members(
    members: Mapping[str, proto7_runner.Proto7RunMember], method: str
) -> tuple[proto7_runner.Proto7RunMember, ...]:
    records = tuple(
        sorted(
            (item for item in members.values() if item.method_label == method),
            key=lambda item: item.point_count,
        )
    )
    if tuple(item.point_count for item in records) != RSP1_POINT_COUNTS:
        raise ValueError(f"RSP1 {method} members differ from the frozen ladder")
    return records


def _event_assessment(
    plan: Mapping[str, Any],
    members: Mapping[str, proto7_runner.Proto7RunMember],
    coordinate_time: float,
) -> dict[str, Any]:
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    primary = _method_members(members, "RK4")
    comparator = _method_members(members, "SSPRK3")
    for member in (*primary, *comparator):
        if np.float64(member.time).tobytes() != np.float64(
            coordinate_time
        ).tobytes():
            raise ValueError("RSP1 common-event times are not bitwise aligned")

    def compose(records: tuple[proto7_runner.Proto7RunMember, ...], method: str):
        return rsp1_gr0_common_event(
            [item.state for item in records],
            [item.initial.grid for item in records],
            accepted_stage_counts=[
                item.transaction.state.accepted_stage_count for item in records
            ],
            method=method,
            coordinate_time=coordinate_time,
            cutoff=float(Q(physical["cutoff_Lambda"])),
            measurement_radius_maximum=float(
                Q(physical["measurement_radius_maximum"])
            ),
            taper_fraction=float(Q(numerics["proper_spectral_taper_fraction"])),
            fixed_outer_rows=numerics["fixed_outer_rows"],
            point_batch_size=numerics["source_point_batch_size"],
        )

    primary_common = compose(primary, "RK4")
    comparator_common = compose(comparator, "SSPRK3")
    return {
        "coordinate_time": coordinate_time,
        "primary_common_event": primary_common,
        "comparator_common_event": comparator_common,
        "endpoint_question_passed": (
            primary_common.study_admission_passed
            and comparator_common.study_admission_passed
        ),
        "collapse_or_trapped_outcome_classified": False,
        "member_diagnostics": {
            item.key: {
                "step_index": item.step_index,
                "transaction_serial": item.transaction_serial,
                "accepted_stage_count": item.transaction.state.accepted_stage_count,
                "CFL_retry_count": item.CFL_retry_count,
                "source_retry_count": item.source_retry_count,
                "accumulated_characteristic_distance": (
                    item.transaction.causal_state.accumulated_characteristic_distance
                ),
                "remaining_boundary_margin": (
                    item.transaction.boundary_geometry.outer_radius
                    - item.transaction.boundary_geometry.measurement_radius
                    - item.transaction.boundary_geometry.sbp_stencil_reach
                    - item.transaction.causal_state.accumulated_characteristic_distance
                    - item.transaction.boundary_geometry.minimum_causal_buffer
                ),
                "state_sha256": array_content_sha256(
                    item.state.u, item.state.p, item.state.q
                ),
            }
            for item in members.values()
        },
    }


def _manifest(
    plan_path: Path,
    authorization_path: Path,
    authorization: Mapping[str, Any],
) -> dict[str, Any]:
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    return {
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "study_id": frozen["study_id"],
        "git_commit": inherited._git_head(),
        "plan_path": plan_path.relative_to(REPOSITORY).as_posix(),
        "plan_sha256": inherited._sha(plan_path),
        "authorization_path": authorization_path.relative_to(
            REPOSITORY
        ).as_posix(),
        "authorization_sha256": inherited._sha(authorization_path),
        "input_manifest_sha256": frozen["input_manifest_sha256"],
        "point_counts": list(RSP1_POINT_COUNTS),
        "methods": ["RK4", "SSPRK3"],
        "amplitude": "3",
        "implementation_sha256": authorization["implementation_sha256"],
        "created_unix_time": wall_time.time(),
        "outcome_fields_present_at_creation": False,
    }


def _checkpoint_metadata(
    *,
    manifest: Mapping[str, Any],
    event_index: int,
    event_log_path: Path,
    members: Mapping[str, proto7_runner.Proto7RunMember],
    terminal: bool,
    terminal_result: Mapping[str, Any] | None,
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "campaign_id": manifest["study_id"],
        "study_id": manifest["study_id"],
        "amplitude": "3",
        "amplitude_index": 0,
        "completed_common_event_index": event_index,
        "consecutive_qualified_events": 0,
        "completed_amplitude_records": [],
        "event_log": inherited._event_log_identity(
            inherited._event_log_bytes(event_log_path)
        ),
        "terminal": terminal,
        "terminal_result": terminal_result,
        "members": {
            key: proto6_runner._member_metadata(item)
            for key, item in members.items()
        },
    }


def _restore_checkpoint(
    path: Path,
    *,
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, proto7_runner.Proto7RunMember], bytes]:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(bytes(archive["metadata_utf8"]).decode("utf-8"))
        event_log_payload = bytes(archive["event_log_utf8"])
        if metadata.get("event_log") != inherited._event_log_identity(
            event_log_payload
        ):
            raise ValueError("RSP1 checkpoint event-log evidence is inconsistent")
        members = _build_members(plan, authorization)
        for key, member in members.items():
            prefix = key.replace("-", "_")
            stored = metadata["members"][key]
            if stored["input_hash"] != member.input_hash:
                raise ValueError("RSP1 checkpoint input hash differs")
            member.state = EvolutionState(
                archive[f"{prefix}_u"],
                archive[f"{prefix}_p"],
                archive[f"{prefix}_q"],
            )
            member.time = float(stored["time"])
            member.step_index = int(stored["step_index"])
            member.transaction_serial = int(stored["transaction_serial"])
            member.CFL_retry_count = int(stored["CFL_retry_count"])
            member.source_retry_count = int(stored["source_retry_count"])
            member.transaction.state = GR0RuntimeMonitorState(
                **stored["runtime_monitor_state"]
            )
            member.transaction.causal_state = CausalBudgetState(
                **stored["causal_state"]
            )
            member.tracers.positions = archive[
                f"{prefix}_tracer_positions"
            ].copy()
            member.tracers.proper_times = archive[
                f"{prefix}_tracer_proper_times"
            ].copy()
            member.tracers.event_proper_times = [
                row.copy() for row in archive[f"{prefix}_event_proper_times"]
            ]
            member.tracers.event_fields = [
                row.copy() for row in archive[f"{prefix}_event_fields"]
            ]
    return metadata, members, event_log_payload


def _validate_resume_manifest(
    manifest: Mapping[str, Any],
    plan_path: Path,
    authorization_path: Path,
    authorization: Mapping[str, Any],
) -> None:
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if (
        manifest.get("runner_id") != RUNNER_ID
        or manifest.get("study_id") != frozen["study_id"]
        or manifest.get("plan_sha256") != inherited._sha(plan_path)
        or manifest.get("authorization_sha256")
        != inherited._sha(authorization_path)
        or manifest.get("input_manifest_sha256")
        != frozen["input_manifest_sha256"]
    ):
        raise ValueError("RSP1 resume manifest differs from the authorization")


def _restored_common_boundary(
    members: Mapping[str, proto7_runner.Proto7RunMember],
    snapshots: Mapping[str, Mapping[str, Any]],
) -> bool:
    for key, member in members.items():
        before = snapshots[key]
        if (
            array_content_sha256(member.state.u, member.state.p, member.state.q)
            != array_content_sha256(
                before["state"].u, before["state"].p, before["state"].q
            )
            or member.time != before["time"]
            or member.step_index != before["step_index"]
            or member.transaction_serial != before["transaction_serial"]
            or member.CFL_retry_count != before["CFL_retry_count"]
            or member.source_retry_count != before["source_retry_count"]
            or member.transaction.state != before["monitor_state"]
            or member.transaction.causal_state != before["causal_state"]
            or not np.array_equal(
                member.tracers.positions, before["tracer_positions"]
            )
            or not np.array_equal(
                member.tracers.proper_times, before["tracer_proper_times"]
            )
            or len(member.tracers.event_proper_times)
            != len(before["event_proper_times"])
            or len(member.tracers.event_fields) != len(before["event_fields"])
            or any(
                not np.array_equal(left, right)
                for left, right in zip(
                    member.tracers.event_proper_times,
                    before["event_proper_times"],
                    strict=True,
                )
            )
            or any(
                not np.array_equal(left, right)
                for left, right in zip(
                    member.tracers.event_fields,
                    before["event_fields"],
                    strict=True,
                )
            )
        ):
            return False
    return True


def run_study(
    *,
    plan_path: Path,
    authorization_path: Path,
    resume: bool,
) -> dict[str, Any]:
    inherited._require_clean_tracked_worktree()
    plan, authorization = _validate_authorization(
        plan_path,
        authorization_path,
        require_fresh_namespace=not resume,
    )
    output_root = REPOSITORY / plan["provenance"]["output_root"]
    manifest_path = REPOSITORY / plan["provenance"]["campaign_manifest"]
    event_log = output_root / plan["provenance"]["append_only_event_log_name"]
    checkpoint_path = output_root / "latest-checkpoint.npz"
    result_path = output_root / "study-result.json"
    resume_recovery: dict[str, Any] | None = None

    if resume:
        if not output_root.is_dir() or not manifest_path.is_file():
            raise ValueError("RSP1 resume requires the existing study namespace")
        if result_path.exists():
            raise ValueError("RSP1 study is already terminal")
        manifest = inherited._load_manifest(manifest_path)
        _validate_resume_manifest(
            manifest, plan_path, authorization_path, authorization
        )
        if not checkpoint_path.is_file():
            raise ValueError("RSP1 resume requires the atomic checkpoint")
        checkpoint, members, checkpoint_event_log = _restore_checkpoint(
            checkpoint_path,
            plan=plan,
            authorization=authorization,
        )
        if checkpoint["study_id"] != manifest["study_id"]:
            raise ValueError("RSP1 checkpoint study identity differs")
        if checkpoint["terminal"] or checkpoint["completed_common_event_index"] != 0:
            raise ValueError("RSP1 resume boundary is not the frozen t0 checkpoint")
        resume_recovery = inherited._recover_event_log_from_checkpoint(
            event_log, checkpoint_event_log
        )
    else:
        output_root.mkdir(parents=True, exist_ok=False)
        manifest = _manifest(plan_path, authorization_path, authorization)
        inherited._write_exclusive(manifest_path, inherited._canonical(manifest))
        members = _build_members(plan, authorization)
        initial = _event_assessment(plan, members, 0.0)
        inherited._append_fsync(
            event_log,
            inherited._canonical_line(
                {
                    "event_type": "initial_premise_event",
                    "amplitude": "3",
                    "event_index": 0,
                    "assessment": initial,
                    "t0_target_ratio_is_not_the_evolved_endpoint_outcome": True,
                }
            ),
        )
        checkpoint = _checkpoint_metadata(
            manifest=manifest,
            event_index=0,
            event_log_path=event_log,
            members=members,
            terminal=False,
            terminal_result=None,
        )
        inherited._write_checkpoint(
            checkpoint_path,
            metadata=checkpoint,
            members=members,
            event_log_path=event_log,
        )

    target_time = float(Q(plan["numerics"]["final_coordinate_time"]))
    retry_factor = float(Q(plan["numerics"]["retry_factor"]))
    maximum_cfl = plan["numerics"]["maximum_CFL_retries_per_step"]
    maximum_source = plan["numerics"]["maximum_source_retries_per_step"]
    minimum_step = float(Q(plan["numerics"]["minimum_step_size"]))
    snapshots = {key: member.snapshot() for key, member in members.items()}
    active_member: proto7_runner.Proto7RunMember | None = None
    stop: dict[str, Any] | None = None
    endpoint: dict[str, Any] | None = None
    started = wall_time.monotonic()

    def rejected_trial_sink(evidence: Mapping[str, Any]) -> None:
        inherited._append_fsync(event_log, inherited._canonical_line(evidence))

    try:
        for active_member in members.values():
            print(
                f"RSP1 advancing {active_member.key} to t={target_time}",
                flush=True,
            )
            active_member.advance_to(
                target_time,
                retry_factor=retry_factor,
                maximum_CFL_retries=maximum_cfl,
                maximum_source_retries=maximum_source,
                minimum_step_size=minimum_step,
                rejected_trial_sink=rejected_trial_sink,
            )
            print(
                f"RSP1 completed {active_member.key}: "
                f"steps={active_member.step_index} "
                f"source_retries={active_member.source_retry_count}",
                flush=True,
            )
        for member in members.values():
            member.append_common_event()
        endpoint = _event_assessment(plan, members, target_time)
    except Proto7TerminalStop as error:
        restore_common_event_snapshots(members, snapshots)
        if not _restored_common_boundary(members, snapshots):
            raise RuntimeError("RSP1 terminal rollback did not restore t0")
        stop = {
            "classification": "stopped_before_resolution_endpoint",
            "reason": error.reason,
            "failures": list(error.failures),
            "complete_failure_evidence": asdict(error.evidence),
            "member": None if active_member is None else active_member.key,
            "cross_member_state_rolled_back": True,
            "cross_member_rollback_verified": True,
        }
    except proto7_runner.Proto7SourceRetryExhausted as error:
        restore_common_event_snapshots(members, snapshots)
        if not _restored_common_boundary(members, snapshots):
            raise RuntimeError("RSP1 source rollback did not restore t0")
        stop = {
            "classification": "stopped_before_resolution_endpoint",
            "reason": error.reason,
            "complete_failure_evidence": error.evidence,
            "member": None if active_member is None else active_member.key,
            "cross_member_state_rolled_back": True,
            "cross_member_rollback_verified": True,
        }
    except Exception as error:
        restore_common_event_snapshots(members, snapshots)
        if not _restored_common_boundary(members, snapshots):
            raise RuntimeError("RSP1 invalid-run rollback did not restore t0")
        stop = {
            "classification": "invalid_implementation_or_nonconverged_run",
            "reason": inherited._source_failure_reason(error),
            "error_type": type(error).__name__,
            "error_message": str(error),
            "member": None if active_member is None else active_member.key,
            "cross_member_state_rolled_back": True,
            "cross_member_rollback_verified": True,
        }

    if endpoint is not None:
        passed = endpoint["endpoint_question_passed"]
        classification = (
            plan["endpoint_classification"]["positive"]
            if passed
            else plan["endpoint_classification"]["negative"]
        )
        event_index = 1
        inherited._append_fsync(
            event_log,
            inherited._canonical_line(
                {
                    "event_type": "evolved_resolution_endpoint",
                    "amplitude": "3",
                    "event_index": 1,
                    "assessment": endpoint,
                    "classification": classification,
                }
            ),
        )
    else:
        passed = False
        classification = stop["classification"]
        event_index = 0
        inherited._append_fsync(
            event_log,
            inherited._canonical_line(
                {
                    "event_type": "terminal_study_stop",
                    "amplitude": "3",
                    "event_index": 0,
                    "stop": stop,
                }
            ),
        )

    result = {
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "study_id": manifest["study_id"],
        "classification": classification,
        "amplitude": "3",
        "target_coordinate_time": target_time,
        "point_counts": list(RSP1_POINT_COUNTS),
        "endpoint_assessment": endpoint,
        "stop": stop,
        "elapsed_wall_seconds": wall_time.monotonic() - started,
        "event_log": event_log.relative_to(REPOSITORY).as_posix(),
        "event_log_sha256": inherited._sha(event_log),
        "resume_recovery": resume_recovery,
        "amplitude_three_spectral_veto_cleared_for_successor_design": passed,
        "fresh_GR0_calibration_completed": False,
        "SGBL_outcome_read": False,
        "FGCQR_outcome_read": False,
        "holdout_execution_authorized": False,
        "retained_EFT_evolution_authorized": False,
        "physical_transition_claim_authorized": False,
        "mechanism_question_answered": False,
    }
    inherited._write_exclusive(result_path, inherited._canonical(result))
    checkpoint = _checkpoint_metadata(
        manifest=manifest,
        event_index=event_index,
        event_log_path=event_log,
        members=members,
        terminal=True,
        terminal_result={
            "classification": classification,
            "amplitude_three_spectral_veto_cleared_for_successor_design": passed,
        },
    )
    inherited._write_checkpoint(
        checkpoint_path,
        metadata=checkpoint,
        members=members,
        event_log_path=event_log,
    )
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
                    "study_id": authorization["artifact_payload"][
                        "frozen_run_plan"
                    ]["study_id"],
                    "amplitude": plan["scope"]["amplitude"],
                    "point_counts": list(RSP1_POINT_COUNTS),
                    "output_created": False,
                },
                sort_keys=True,
            )
        )
        return
    result = run_study(
        plan_path=plan_path,
        authorization_path=authorization_path,
        resume=args.resume,
    )
    print(json.dumps(inherited._serial(result), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
