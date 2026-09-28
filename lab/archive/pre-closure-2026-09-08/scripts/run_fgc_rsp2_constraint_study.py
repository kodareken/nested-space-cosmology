#!/usr/bin/env python3
"""Run the FGC-1-RSP2-FRZ1 one-member late-event constraint study."""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from hashlib import sha256
import json
from pathlib import Path
import sys
import time as wall_time
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_rsp2_frz1 as freeze  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_gr0_calibration_v6 as proto6_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v7 as proto7_runner  # noqa: E402
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    make_gr0_center_boundary_projector,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeMonitorState,
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.proto7_runtime import (  # noqa: E402
    Proto7TerminalStop,
)
from recursive_horizons.fgc.evolution.rsp2_checkpoint import (  # noqa: E402
    load_rsp2_cal9_predecessors,
)
from recursive_horizons.fgc.evolution.rsp2_constraint_runtime import (  # noqa: E402
    rsp2_constraint_assessment,
)


DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-rsp2-run1.toml"
DEFAULT_AUTHORIZATION = REPOSITORY / "results/fgc-1-rsp2-frz1.json"
RUNNER_ID = "FGC-1-RSP2-RUN1-RUNNER"


def _typed_runtime_stop(
    *,
    reason: str,
    member: str,
    evidence: Mapping[str, Any],
    failures: tuple[str, ...] = (),
) -> dict[str, Any]:
    """Return a fail-closed premise stop that cannot imply persistence."""

    return {
        "classification": "stopped_before_RSP2_endpoint",
        "reason": reason,
        "failures": list(failures),
        "complete_failure_evidence": dict(evidence),
        "member": member,
        "runtime_stop_is_not_target_persistence": True,
    }


def _invalid_runtime_stop(error: Exception, *, member: str) -> dict[str, Any]:
    """Return an invalid-run record without converting it into science."""

    return {
        "classification": "invalid_implementation_or_nonconverged_run",
        "reason": inherited._source_failure_reason(error),
        "error_type": type(error).__name__,
        "error_message": str(error),
        "member": member,
        "invalid_run_is_not_target_persistence": True,
    }


def _validate_authorization(
    plan_path: Path,
    authorization_path: Path,
    *,
    require_fresh_namespace: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = freeze.validate_run_plan(plan_path)
    authorization = freeze.load_canonical_result(authorization_path)
    gates = authorization.get("gate_status", {})
    frozen = authorization.get("artifact_payload", {}).get("frozen_run_plan", {})
    if (
        authorization.get("artifact_id") != freeze.ARTIFACT_ID
        or gates.get("RSP2_runtime_implemented") is not True
        or gates.get("RSP2_execution_authorized") is not True
        or gates.get("RSP2_outcome_read") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or gates.get("retained_EFT_evolution_authorized") is not False
        or frozen.get("run_plan_sha256") != freeze._sha(plan_path)
        or frozen.get("combined_point_counts") != [4097, 8193, 16385]
        or frozen.get("minimum_finest_pair_order") != 1.5
        or frozen.get("only_the_new_16385_member_is_advanced") is not True
    ):
        raise ValueError("RSP2 authorization is absent, incomplete, or over-broad")
    for ledger_name in (
        "source_config_sha256",
        "predecessor_sha256",
        "implementation_sha256",
    ):
        ledger = authorization.get(ledger_name)
        if not isinstance(ledger, dict) or not ledger:
            raise ValueError(f"authorized RSP2 {ledger_name} is absent")
        for relative, expected in ledger.items():
            path = REPOSITORY / relative
            if not path.is_file() or freeze._sha(path) != expected:
                raise ValueError(f"authorized RSP2 file drifted: {relative}")
    owner = REPOSITORY / authorization.get("derivation_document", "")
    if (
        not owner.is_file()
        or freeze._sha(owner) != authorization.get("derivation_document_sha256")
    ):
        raise ValueError("authorized RSP2 owner document drifted")
    if require_fresh_namespace:
        output_root = REPOSITORY / plan["provenance"]["output_root"]
        if output_root.exists():
            raise ValueError("fresh RSP2 output root already exists")
    return plan, authorization


def _build_member(
    plan: Mapping[str, Any], authorization: Mapping[str, Any]
) -> proto7_runner.Proto7RunMember:
    initial, state = freeze.build_rsp2_initial_state(plan)
    frozen_input = authorization["artifact_payload"]["new_run_input"]
    state_hash = array_content_sha256(state.u, state.p, state.q)
    if state_hash != frozen_input["projected_state_sha256"]:
        raise ValueError("RSP2 runtime input differs from the frozen 16385 state")
    operator = freeze.rsp2_operator(plan, initial.grid)
    initial_rhs = operator(0.0, state)
    diagnostics = initial_rhs.diagnostics
    if (
        diagnostics.get("source_raw_gate_passed") is not True
        or diagnostics.get("SRC4_tensor_contracted_reference_source") is not True
        or diagnostics.get("PROTO11_interior_q_reprojected") is not False
    ):
        raise ValueError("RSP2 runtime source precheck failed")

    numerics = plan["numerics"]
    physical = plan["physical_inputs"]
    thresholds = plan["universal_thresholds"]
    derivative = SBPFirstDerivative(initial.grid, 2)
    transaction = GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(
            source_residual_maximum=float(
                freeze.Q(thresholds["source_residual_infinity_max"])
            ),
            source_iteration_maximum=thresholds["source_iteration_max"],
            kinetic_condition_maximum=float(
                freeze.Q(thresholds["kinetic_condition_number_max"])
            ),
        ),
        causal_state=CausalBudgetState(
            previous_speed_upper=float(diagnostics["coordinate_speed_upper"])
        ),
        boundary_geometry=BoundaryGeometry(
            float(freeze.Q(physical["outer_radius"])),
            float(freeze.Q(physical["measurement_radius_maximum"])),
            float(freeze.Q(thresholds["minimum_boundary_causal_buffer"])),
            derivative.stencil_reach_intervals * initial.grid.spacing,
        ),
        grid_spacing=initial.grid.spacing,
        cfl_maximum=float(freeze.Q(numerics["cfl_maximum"])),
        hat_normal_factor=float(freeze.Q(thresholds["hat_normal_factor"])),
    )
    tracers = inherited.NormalFlowTracers.create(
        minimum=float(freeze.Q(numerics["normal_flow_tracer_radius_minimum"])),
        maximum=float(freeze.Q(numerics["normal_flow_tracer_radius_maximum"])),
        spacing=float(freeze.Q(numerics["normal_flow_tracer_spacing"])),
        state=state,
        coordinates=initial.grid.coordinates,
        cutoff=float(freeze.Q(physical["cutoff_Lambda"])),
        outer_radius=float(freeze.Q(physical["outer_radius"])),
    )
    projected_initial = replace(initial, state=state)
    return proto7_runner.Proto7RunMember(
        amplitude="3",
        method_label="SSPRK3",
        integrator_id=plan["method"]["integrator_id"],
        spatial_order=2,
        point_count=16385,
        input_hash=frozen_input["expanded_run_config_sha256"],
        initial=projected_initial,
        state=state,
        operator=operator,
        projector=make_gr0_center_boundary_projector(
            projected_initial, fixed_outer_rows=numerics["fixed_outer_rows"]
        ),
        transaction=transaction,
        tracers=tracers,
    )


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
        "plan_sha256": freeze._sha(plan_path),
        "authorization_path": authorization_path.relative_to(REPOSITORY).as_posix(),
        "authorization_sha256": freeze._sha(authorization_path),
        "input_manifest_sha256": frozen["input_manifest_sha256"],
        "predecessor_checkpoint_sha256": authorization["predecessor_sha256"][
            "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz"
        ],
        "amplitude": "3",
        "method": "SSPRK3",
        "new_point_count": 16385,
        "implementation_sha256": authorization["implementation_sha256"],
        "created_unix_time": wall_time.time(),
        "outcome_fields_present_at_creation": False,
    }


def _member_diagnostics(member: proto7_runner.Proto7RunMember) -> dict[str, Any]:
    boundary = member.transaction.boundary_geometry
    causal = member.transaction.causal_state
    return {
        "time": member.time,
        "step_index": member.step_index,
        "transaction_serial": member.transaction_serial,
        "accepted_stage_count": member.transaction.state.accepted_stage_count,
        "CFL_retry_count": member.CFL_retry_count,
        "source_retry_count": member.source_retry_count,
        "accumulated_characteristic_distance": causal.accumulated_characteristic_distance,
        "remaining_boundary_margin": (
            boundary.outer_radius
            - boundary.measurement_radius
            - boundary.sbp_stencil_reach
            - causal.accumulated_characteristic_distance
            - boundary.minimum_causal_buffer
        ),
        "state_sha256": array_content_sha256(
            member.state.u, member.state.p, member.state.q
        ),
    }


def _checkpoint_metadata(
    *,
    manifest: Mapping[str, Any],
    boundary_index: int,
    event_log_path: Path,
    member: proto7_runner.Proto7RunMember,
    terminal: bool,
    terminal_result: Mapping[str, Any] | None,
) -> dict[str, Any]:
    member_metadata = proto6_runner._member_metadata(member)
    member_metadata["state_sha256"] = array_content_sha256(
        member.state.u, member.state.p, member.state.q
    )
    return {
        "schema_version": 1,
        "campaign_id": manifest["study_id"],
        "study_id": manifest["study_id"],
        "manifest_sha256": sha256(freeze._canonical(manifest)).hexdigest(),
        "amplitude": "3",
        "amplitude_index": 0,
        "completed_common_event_index": boundary_index,
        "consecutive_qualified_events": 0,
        "completed_amplitude_records": [],
        "event_log": inherited._event_log_identity(
            inherited._event_log_bytes(event_log_path)
        ),
        "terminal": terminal,
        "terminal_result": terminal_result,
        "members": {member.key: member_metadata},
    }


def _restore_checkpoint(
    path: Path,
    *,
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
    manifest: Mapping[str, Any],
) -> tuple[dict[str, Any], proto7_runner.Proto7RunMember, bytes]:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(bytes(archive["metadata_utf8"]).decode("utf-8"))
        event_log_payload = bytes(archive["event_log_utf8"])
        members = metadata.get("members")
        boundary_index = metadata.get("completed_common_event_index")
        if (
            metadata.get("schema_version") != 1
            or metadata.get("campaign_id") != manifest.get("study_id")
            or metadata.get("study_id") != manifest.get("study_id")
            or metadata.get("manifest_sha256")
            != sha256(freeze._canonical(manifest)).hexdigest()
            or metadata.get("amplitude") != "3"
            or metadata.get("amplitude_index") != 0
            or isinstance(boundary_index, bool)
            or not isinstance(boundary_index, int)
            or not 0 <= boundary_index <= 23
            or not isinstance(members, dict)
            or set(members) != {"SSPRK3-16385"}
            or metadata.get("event_log")
            != inherited._event_log_identity(event_log_payload)
        ):
            raise ValueError("RSP2 checkpoint event-log evidence is inconsistent")
        member = _build_member(plan, authorization)
        key = member.key
        prefix = key.replace("-", "_")
        stored = members[key]
        if stored["input_hash"] != member.input_hash:
            raise ValueError("RSP2 checkpoint input hash differs")
        restored_state = EvolutionState(
            archive[f"{prefix}_u"].copy(),
            archive[f"{prefix}_p"].copy(),
            archive[f"{prefix}_q"].copy(),
        )
        if (
            array_content_sha256(
                restored_state.u, restored_state.p, restored_state.q
            )
            != stored.get("state_sha256")
        ):
            raise ValueError("RSP2 checkpoint state hash differs")
        member.state = restored_state
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
        member.tracers.positions = archive[f"{prefix}_tracer_positions"].copy()
        member.tracers.proper_times = archive[
            f"{prefix}_tracer_proper_times"
        ].copy()
        member.tracers.event_proper_times = [
            row.copy() for row in archive[f"{prefix}_event_proper_times"]
        ]
        member.tracers.event_fields = [
            row.copy() for row in archive[f"{prefix}_event_fields"]
        ]
    return metadata, member, event_log_payload


def _validate_resume_manifest(
    manifest: Mapping[str, Any],
    plan_path: Path,
    authorization_path: Path,
    authorization: Mapping[str, Any],
) -> None:
    expected = authorization["artifact_payload"]["frozen_run_plan"]
    if (
        manifest.get("runner_id") != RUNNER_ID
        or manifest.get("study_id") != expected["study_id"]
        or manifest.get("plan_sha256") != freeze._sha(plan_path)
        or manifest.get("authorization_sha256") != freeze._sha(authorization_path)
        or manifest.get("input_manifest_sha256")
        != expected["input_manifest_sha256"]
        or manifest.get("git_commit") != inherited._git_head()
        or manifest.get("implementation_sha256")
        != authorization.get("implementation_sha256")
        or manifest.get("predecessor_checkpoint_sha256")
        != authorization.get("predecessor_sha256", {}).get(
            "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz"
        )
    ):
        raise ValueError("resume manifest differs from the authorized RSP2 study")


def _endpoint_assessment(
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
    member: proto7_runner.Proto7RunMember,
):
    checkpoint = REPOSITORY / plan["provenance"]["CAL9_checkpoint"]
    predecessors = load_rsp2_cal9_predecessors(
        checkpoint,
        expected_sha256=authorization["predecessor_sha256"][
            plan["provenance"]["CAL9_checkpoint"]
        ],
        outer_radius=float(freeze.Q(plan["physical_inputs"]["outer_radius"])),
    )
    return rsp2_constraint_assessment(
        [*predecessors.states, member.state],
        [*predecessors.grids, member.initial.grid],
        accepted_stage_counts=[
            *predecessors.accepted_stage_counts,
            member.transaction.state.accepted_stage_count,
        ],
        coordinate_time=member.time,
        fixed_outer_rows=plan["numerics"]["fixed_outer_rows"],
        coarsest_guard_maximum=float(
            freeze.Q(plan["method"]["coarsest_constraint_guard"])
        ),
        finest_guard_maximum=float(
            freeze.Q(plan["method"]["finest_constraint_guard"])
        ),
        minimum_finest_pair_order=float(
            freeze.Q(plan["numerics"]["minimum_constraint_finest_pair_order"])
        ),
        planck_mass=float(freeze.Q(plan["physical_inputs"]["planck_mass"])),
        scalar_mass=float(freeze.Q(plan["physical_inputs"]["scalar_mass"])),
        quartic_coupling=float(
            freeze.Q(plan["physical_inputs"]["quartic_coupling"])
        ),
        length_unit=float(freeze.Q(plan["physical_inputs"]["length_unit_L0"])),
    )


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
    manifest_path = REPOSITORY / plan["provenance"]["study_manifest"]
    event_log = output_root / plan["provenance"]["append_only_event_log_name"]
    checkpoint_path = output_root / plan["provenance"]["checkpoint_name"]
    result_path = output_root / plan["provenance"]["result_name"]
    resume_recovery: dict[str, Any] | None = None

    if resume:
        if not output_root.is_dir() or not manifest_path.is_file():
            raise ValueError("RSP2 resume requires the existing study namespace")
        if result_path.exists():
            raise ValueError("RSP2 study result already exists; resume is unnecessary")
        manifest = inherited._load_manifest(manifest_path)
        _validate_resume_manifest(
            manifest, plan_path, authorization_path, authorization
        )
        if not checkpoint_path.is_file():
            raise ValueError("RSP2 resume requires the latest atomic checkpoint")
        checkpoint, member, checkpoint_event_log = _restore_checkpoint(
            checkpoint_path,
            plan=plan,
            authorization=authorization,
            manifest=manifest,
        )
        if checkpoint["study_id"] != manifest["study_id"]:
            raise ValueError("RSP2 checkpoint study identity differs")
        if checkpoint.get("terminal") is True:
            raise ValueError("RSP2 terminal checkpoint cannot resume")
        boundary_index = int(checkpoint["completed_common_event_index"])
        resume_recovery = inherited._recover_event_log_from_checkpoint(
            event_log, checkpoint_event_log
        )
        if checkpoint.get("terminal") is True:
            recovered = checkpoint.get("terminal_result")
            if (
                not isinstance(recovered, dict)
                or recovered.get("study_id") != manifest["study_id"]
                or recovered.get("event_log_sha256") != freeze._sha(event_log)
                or recovered.get("classification")
                not in {
                    "completed_target_and_complete_constraint_pass",
                    "completed_target_cleared_other_constraint_obstruction",
                    "completed_target_order_persistent_below_threshold",
                    "completed_target_persistent_with_other_constraint_obstruction",
                    "completed_target_inconclusive_premise_failure",
                    "stopped_before_RSP2_endpoint",
                    "invalid_implementation_or_nonconverged_run",
                }
            ):
                raise ValueError("RSP2 terminal checkpoint result is inconsistent")
            inherited._write_exclusive(result_path, freeze._canonical(recovered))
            return freeze._serial(recovered)
        interval = float(freeze.Q(plan["numerics"]["progress_checkpoint_interval"]))
        expected_resume_time = boundary_index * interval
        if np.float64(member.time).tobytes() != np.float64(
            expected_resume_time
        ).tobytes():
            raise ValueError("RSP2 nonterminal checkpoint is off its frozen boundary")
    else:
        # Reconstruct and source-check the frozen input before creating any raw
        # namespace.  A failed prelaunch premise therefore leaves no orphaned
        # run root that could be mistaken for a partially executed study.
        member = _build_member(plan, authorization)
        output_root.mkdir(parents=True, exist_ok=False)
        manifest = _manifest(plan_path, authorization_path, authorization)
        inherited._write_exclusive(manifest_path, freeze._canonical(manifest))
        boundary_index = 0
        inherited._append_fsync(
            event_log,
            inherited._canonical_line(
                {
                    "event_type": "initial_premise_event",
                    "boundary_index": 0,
                    "coordinate_time": 0.0,
                    "member": _member_diagnostics(member),
                    "trajectory_advanced": False,
                }
            ),
        )
        checkpoint = _checkpoint_metadata(
            manifest=manifest,
            boundary_index=0,
            event_log_path=event_log,
            member=member,
            terminal=False,
            terminal_result=None,
        )
        inherited._write_checkpoint(
            checkpoint_path,
            metadata=checkpoint,
            members={member.key: member},
            event_log_path=event_log,
        )

    target_time = float(freeze.Q(plan["numerics"]["final_coordinate_time"]))
    interval = float(freeze.Q(plan["numerics"]["progress_checkpoint_interval"]))
    retry_factor = float(freeze.Q(plan["numerics"]["retry_factor"]))
    maximum_cfl = plan["numerics"]["maximum_CFL_retries_per_step"]
    maximum_source = plan["numerics"]["maximum_source_retries_per_step"]
    minimum_step = float(freeze.Q(plan["numerics"]["minimum_step_size"]))
    started = wall_time.monotonic()
    endpoint = None
    stop = None

    def rejected_trial_sink(evidence: Mapping[str, Any]) -> None:
        inherited._append_fsync(event_log, inherited._canonical_line(evidence))

    try:
        final_boundary = int(round(target_time / interval))
        for next_index in range(boundary_index + 1, final_boundary + 1):
            next_time = next_index * interval
            print(
                f"RSP2 advancing {member.key} to t={next_time:.10g} "
                f"({next_index}/{final_boundary})",
                flush=True,
            )
            member.advance_to(
                next_time,
                retry_factor=retry_factor,
                maximum_CFL_retries=maximum_cfl,
                maximum_source_retries=maximum_source,
                minimum_step_size=minimum_step,
                rejected_trial_sink=rejected_trial_sink,
            )
            member.append_common_event()
            boundary_index = next_index
            inherited._append_fsync(
                event_log,
                inherited._canonical_line(
                    {
                        "event_type": "accepted_progress_boundary",
                        "boundary_index": boundary_index,
                        "coordinate_time": member.time,
                        "member": _member_diagnostics(member),
                    }
                ),
            )
            checkpoint = _checkpoint_metadata(
                manifest=manifest,
                boundary_index=boundary_index,
                event_log_path=event_log,
                member=member,
                terminal=False,
                terminal_result=None,
            )
            inherited._write_checkpoint(
                checkpoint_path,
                metadata=checkpoint,
                members={member.key: member},
                event_log_path=event_log,
            )
        endpoint = _endpoint_assessment(plan, authorization, member)
    except Proto7TerminalStop as error:
        stop = _typed_runtime_stop(
            reason=error.reason,
            member=member.key,
            evidence=asdict(error.evidence),
            failures=tuple(error.failures),
        )
    except proto7_runner.Proto7SourceRetryExhausted as error:
        stop = _typed_runtime_stop(
            reason=error.reason,
            member=member.key,
            evidence=error.evidence,
        )
    except Exception as error:
        stop = _invalid_runtime_stop(error, member=member.key)

    if endpoint is not None:
        classification = endpoint.classification
        inherited._append_fsync(
            event_log,
            inherited._canonical_line(
                {
                    "event_type": "evolved_RSP2_endpoint",
                    "boundary_index": boundary_index,
                    "coordinate_time": member.time,
                    "classification": classification,
                    "assessment": freeze._serial(endpoint),
                }
            ),
        )
    else:
        classification = stop["classification"]
        inherited._append_fsync(
            event_log,
            inherited._canonical_line(
                {
                    "event_type": "terminal_study_stop",
                    "boundary_index": boundary_index,
                    "coordinate_time": member.time,
                    "stop": stop,
                }
            ),
        )

    result = freeze._serial({
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "study_id": manifest["study_id"],
        "classification": classification,
        "amplitude": "3",
        "method": "SSPRK3",
        "target_coordinate_time": target_time,
        "combined_point_counts": [4097, 8193, 16385],
        "endpoint_assessment": endpoint,
        "stop": stop,
        "new_member": _member_diagnostics(member),
        "elapsed_wall_seconds": wall_time.monotonic() - started,
        "event_log": event_log.relative_to(REPOSITORY).as_posix(),
        "event_log_sha256": freeze._sha(event_log),
        "resume_recovery": resume_recovery,
        "RSP2_target_order_cleared": bool(
            endpoint is not None and endpoint.target_order_passed
        ),
        "RSP2_complete_constraint_admission_passed": bool(
            endpoint is not None and endpoint.complete_constraint_admission_passed
        ),
        "PROTO13_frozen": False,
        "fresh_GR0_dynamic_calibration_completed": False,
        "SGBL_outcome_read": False,
        "FGCQR_outcome_read": False,
        "holdout_execution_authorized": False,
        "retained_EFT_evolution_authorized": False,
        "physical_transition_claim_authorized": False,
        "mechanism_question_answered": False,
    })
    checkpoint = _checkpoint_metadata(
        manifest=manifest,
        boundary_index=boundary_index,
        event_log_path=event_log,
        member=member,
        terminal=True,
        terminal_result=result,
    )
    inherited._write_checkpoint(
        checkpoint_path,
        metadata=checkpoint,
        members={member.key: member},
        event_log_path=event_log,
    )
    inherited._write_exclusive(result_path, freeze._canonical(result))
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
            plan_path, authorization_path, require_fresh_namespace=True
        )
        print(
            json.dumps(
                {
                    "authorized": True,
                    "study_id": authorization["artifact_payload"][
                        "frozen_run_plan"
                    ]["study_id"],
                    "new_point_count": plan["numerics"]["new_resolution"],
                    "target_coordinate_time": plan["numerics"][
                        "final_coordinate_time"
                    ],
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
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
