#!/usr/bin/env python3
"""Execute the HLT6-authorized fresh PROTO8 GR-0 calibration campaign.

Evolution, source solving, proposal ownership, rollback, checkpoints, temporal
spectra, trapped-sign errors, thresholds, methods, grids, and physical inputs
are inherited unchanged from the PROTO7 runner.  This successor replaces only
the committed common-event constraint and spatial-spectral compositor frozen
by PROTO8.  It cannot open SGB-L or FGC-QR outcomes.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import time as wall_time
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_hlt6_mon6 as hlt6  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_gr0_calibration_v6 as proto6_runner  # noqa: E402
from scripts import run_fgc_gr0_calibration_v7 as proto7_runner  # noqa: E402
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_SPECTRAL_FIELD_ORDER,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeStop,
    proto5_calibration_trapped_assessment,
    proto5_temporal_spectral_admission,
)
from recursive_horizons.fgc.evolution.proto6_runtime import (  # noqa: E402
    restore_common_event_snapshots,
)
from recursive_horizons.fgc.evolution.proto7_runtime import (  # noqa: E402
    Proto7TerminalStop,
)
from recursive_horizons.fgc.evolution.proto8_runtime import (  # noqa: E402
    proto8_gr0_common_event,
)


Q = inherited.Q
DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-cal5-run1.toml"
DEFAULT_AUTHORIZATION = REPOSITORY / "results/fgc-1-hlt6-mon6.json"
RUNNER_ID = "FGC-1-CAL5-RUN1-RUNNER"
Proto8RunMember = proto7_runner.Proto7RunMember
Proto8SourceRetryExhausted = proto7_runner.Proto7SourceRetryExhausted


def _validate_authorization(
    plan_path: Path,
    authorization_path: Path,
    *,
    require_fresh_namespace: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = hlt6.validate_run_plan(plan_path)
    authorization = hlt6.load_canonical_result(authorization_path)
    gates = authorization.get("gate_status", {})
    if (
        authorization.get("artifact_id") != hlt6.ARTIFACT_ID
        or gates.get("PROTO8_successor_runtime_compositor_implemented") is not True
        or gates.get("PROTO8_fresh_GR0_dynamic_calibration_authorized") is not True
        or gates.get("PROTO8_resolved_holdout_manifest_authorized") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or gates.get("retained_EFT_evolution_authorized") is not False
    ):
        raise ValueError("HLT6 authorization is absent, incomplete, or over-broad")
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if frozen["run_plan_sha256"] != inherited._sha(plan_path):
        raise ValueError("CAL5 run-plan hash differs from HLT6")
    for relative, expected in authorization["implementation_sha256"].items():
        current = REPOSITORY / relative
        if not current.is_file() or inherited._sha(current) != expected:
            raise ValueError(f"authorized implementation drifted: {relative}")
    if require_fresh_namespace:
        reproduced = hlt6.record(hlt6.DEFAULT_CONFIG)
        if inherited._serial(reproduced) != authorization:
            raise ValueError("HLT6 authorization does not reproduce before launch")
    return plan, authorization


def _build_members(
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
    amplitude: str,
) -> dict[str, Proto8RunMember]:
    members = proto7_runner._build_members(plan, authorization, amplitude)
    if len(members) != 6:
        raise RuntimeError("PROTO8 GR-0 amplitude requires six run members")
    return members


def _checkpoint_metadata(
    *,
    manifest: Mapping[str, Any],
    amplitude: str,
    amplitude_index: int,
    event_index: int,
    consecutive_qualified_events: int,
    amplitude_records: list[Mapping[str, Any]],
    event_log_path: Path,
    members: Mapping[str, Proto8RunMember],
    terminal: bool,
) -> dict[str, Any]:
    return proto7_runner._checkpoint_metadata(
        manifest=manifest,
        amplitude=amplitude,
        amplitude_index=amplitude_index,
        event_index=event_index,
        consecutive_qualified_events=consecutive_qualified_events,
        amplitude_records=amplitude_records,
        event_log_path=event_log_path,
        members=members,
        terminal=terminal,
    )


def _restore_checkpoint(
    path: Path,
    *,
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Proto8RunMember], bytes]:
    return proto7_runner._restore_checkpoint(
        path,
        plan=plan,
        authorization=authorization,
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
        "campaign_id": frozen["campaign_id"],
        "git_commit": inherited._git_head(),
        "plan_path": plan_path.relative_to(REPOSITORY).as_posix(),
        "plan_sha256": inherited._sha(plan_path),
        "authorization_path": authorization_path.relative_to(REPOSITORY).as_posix(),
        "authorization_sha256": inherited._sha(authorization_path),
        "input_manifest_sha256": frozen["input_manifest_sha256"],
        "ordered_amplitudes": frozen["ordered_amplitudes"],
        "implementation_sha256": authorization["implementation_sha256"],
        "created_unix_time": wall_time.time(),
        "outcome_fields_present_at_creation": False,
    }


def _validate_resume_manifest(
    manifest: Mapping[str, Any],
    plan_path: Path,
    authorization_path: Path,
    authorization: Mapping[str, Any],
) -> None:
    expected = authorization["artifact_payload"]["frozen_run_plan"]
    if (
        manifest.get("runner_id") != RUNNER_ID
        or manifest.get("campaign_id") != expected["campaign_id"]
        or manifest.get("plan_sha256") != inherited._sha(plan_path)
        or manifest.get("authorization_sha256") != inherited._sha(authorization_path)
        or manifest.get("input_manifest_sha256") != expected["input_manifest_sha256"]
    ):
        raise ValueError("resume manifest differs from the authorized PROTO8 campaign")


def _method_members(
    members: Mapping[str, Proto8RunMember], method: str
) -> tuple[Proto8RunMember, ...]:
    records = tuple(
        sorted(
            (item for item in members.values() if item.method_label == method),
            key=lambda item: item.point_count,
        )
    )
    if len(records) != 3:
        raise RuntimeError(f"method {method} does not own three members")
    return records


def _event_assessment(
    plan: Mapping[str, Any],
    members: Mapping[str, Proto8RunMember],
    coordinate_time: float,
) -> dict[str, Any]:
    """Evaluate only after all six members committed to the same event."""

    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    primary = _method_members(members, "RK4")
    comparator = _method_members(members, "SSPRK3")
    for member in (*primary, *comparator):
        if np.float64(member.time).tobytes() != np.float64(coordinate_time).tobytes():
            raise ValueError("PROTO8 common-event member time is not bitwise aligned")

    primary_common = proto8_gr0_common_event(
        [item.state for item in primary],
        [item.initial.grid for item in primary],
        accepted_stage_counts=[
            item.transaction.state.accepted_stage_count for item in primary
        ],
        method="RK4",
        coordinate_time=coordinate_time,
        cutoff=float(Q(physical["cutoff_Lambda"])),
        measurement_radius_maximum=float(Q(physical["measurement_radius_maximum"])),
        taper_fraction=float(Q(numerics["proper_spectral_taper_fraction"])),
        fixed_outer_rows=numerics["fixed_outer_rows"],
    )
    comparator_common = proto8_gr0_common_event(
        [item.state for item in comparator],
        [item.initial.grid for item in comparator],
        accepted_stage_counts=[
            item.transaction.state.accepted_stage_count for item in comparator
        ],
        method="SSPRK3",
        coordinate_time=coordinate_time,
        cutoff=float(Q(physical["cutoff_Lambda"])),
        measurement_radius_maximum=float(Q(physical["measurement_radius_maximum"])),
        taper_fraction=float(Q(numerics["proper_spectral_taper_fraction"])),
        fixed_outer_rows=numerics["fixed_outer_rows"],
    )
    trapped = proto5_calibration_trapped_assessment(
        primary_states=[item.state for item in primary],
        comparator_states=[item.state for item in comparator],
        grids=[item.initial.grid for item in primary],
        primary_common_event=primary_common,
        comparator_common_event=comparator_common,
        measurement_radius_maximum=float(Q(physical["measurement_radius_maximum"])),
        minimum_observed_order=float(Q(numerics["minimum_constraint_finest_pair_order"])),
        positive_margin_factor=float(
            plan["candidate_selection"]["trapped_sign_margin_over_combined_error_factor"]
        ),
    )
    temporal: dict[str, Any] = {
        "available": False,
        "required_minimum_samples": numerics["temporal_history_minimum_samples"],
        "admission_passed": False,
    }
    sample_count = len(primary[0].tracers.event_fields)
    if sample_count >= numerics["temporal_history_minimum_samples"]:
        method_temporal = {}
        for method, records in (("RK4", primary), ("SSPRK3", comparator)):
            method_temporal[method] = proto5_temporal_spectral_admission(
                point_counts=[item.point_count for item in records],
                proper_times_by_resolution=[
                    np.asarray(item.tracers.event_proper_times) for item in records
                ],
                field_histories_by_resolution=[
                    np.asarray(item.tracers.event_fields) for item in records
                ],
                field_names=PROTO4_SPECTRAL_FIELD_ORDER,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                taper_fraction=float(Q(numerics["proper_spectral_taper_fraction"])),
            )
        temporal = {
            "available": True,
            "sample_count": sample_count,
            "methods": method_temporal,
            "admission_passed": all(
                item.admission_passed for item in method_temporal.values()
            ),
        }
    qualified = trapped.trapped_sign_passed and temporal["admission_passed"]
    return {
        "coordinate_time": coordinate_time,
        "primary_common_event": primary_common,
        "comparator_common_event": comparator_common,
        "temporal_spectral_admission": temporal,
        "trapped_assessment": trapped,
        "qualified_trapped_common_event": qualified,
        "member_diagnostics": {
            item.key: {
                "step_index": item.step_index,
                "transaction_serial": item.transaction_serial,
                "accepted_stage_count": item.transaction.state.accepted_stage_count,
                "CFL_retry_count": item.CFL_retry_count,
                "source_retry_count": item.source_retry_count,
                "accumulated_characteristic_distance": item.transaction.causal_state.accumulated_characteristic_distance,
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
                "minimum_tracer_radius": float(np.min(item.tracers.positions)),
                "maximum_tracer_radius": float(np.max(item.tracers.positions)),
            }
            for item in members.values()
        },
    }


def run_campaign(
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
    result_path = output_root / "campaign-result.json"
    resume_recovery: dict[str, Any] | None = None

    if resume:
        if not output_root.is_dir() or not manifest_path.is_file():
            raise ValueError("resume requires the existing authorized campaign")
        if result_path.exists():
            raise ValueError("campaign is already terminal; overwrite is forbidden")
        manifest = inherited._load_manifest(manifest_path)
        _validate_resume_manifest(manifest, plan_path, authorization_path, authorization)
        if not checkpoint_path.is_file():
            raise ValueError("resume requires the atomic latest checkpoint")
        checkpoint, members, checkpoint_event_log = _restore_checkpoint(
            checkpoint_path,
            plan=plan,
            authorization=authorization,
        )
        if checkpoint["campaign_id"] != manifest["campaign_id"]:
            raise ValueError("checkpoint campaign identity differs")
        amplitude_index = int(checkpoint["amplitude_index"])
        event_index = int(checkpoint["completed_common_event_index"])
        consecutive = int(checkpoint["consecutive_qualified_events"])
        amplitude_records = list(checkpoint["completed_amplitude_records"])
        resume_recovery = inherited._recover_event_log_from_checkpoint(
            event_log, checkpoint_event_log
        )
    else:
        if output_root.exists():
            raise ValueError("fresh PROTO8 calibration output root already exists")
        output_root.mkdir(parents=True, exist_ok=False)
        manifest = _manifest(plan_path, authorization_path, authorization)
        inherited._write_exclusive(manifest_path, inherited._canonical(manifest))
        amplitude_index = 0
        event_index = 0
        consecutive = 0
        amplitude_records: list[dict[str, Any]] = []
        amplitude = plan["candidate_selection"]["ordered_amplitudes"][0]
        members = _build_members(plan, authorization, amplitude)
        initial = _event_assessment(plan, members, 0.0)
        inherited._append_fsync(
            event_log,
            inherited._canonical_line(
                {
                    "event_type": "common_event",
                    "amplitude": amplitude,
                    "event_index": 0,
                    "assessment": initial,
                }
            ),
        )
        checkpoint = _checkpoint_metadata(
            manifest=manifest,
            amplitude=amplitude,
            amplitude_index=amplitude_index,
            event_index=0,
            consecutive_qualified_events=0,
            amplitude_records=amplitude_records,
            event_log_path=event_log,
            members=members,
            terminal=False,
        )
        proto6_runner._write_checkpoint(
            checkpoint_path,
            metadata=checkpoint,
            members=members,
            event_log_path=event_log,
        )

    amplitudes = plan["candidate_selection"]["ordered_amplitudes"]
    output_interval = float(Q(plan["numerics"]["common_output_interval"]))
    checkpoint_interval = float(Q(plan["numerics"]["checkpoint_interval"]))
    checkpoint_ratio = checkpoint_interval / output_interval
    checkpoint_every_events = int(round(checkpoint_ratio))
    if (
        checkpoint_every_events < 1
        or abs(checkpoint_ratio - checkpoint_every_events) > 1.0e-12
    ):
        raise ValueError("checkpoint interval must be an integer number of events")
    final_time = float(Q(plan["numerics"]["final_coordinate_time"]))
    final_event_index = int(round(final_time / output_interval))
    retry_factor = float(Q(plan["numerics"]["retry_factor"]))
    maximum_CFL_retries = plan["numerics"]["maximum_CFL_retries_per_step"]
    maximum_source_retries = plan["numerics"]["maximum_source_retries_per_step"]
    minimum_step = float(Q(plan["numerics"]["minimum_step_size"]))
    required_consecutive = plan["candidate_selection"][
        "minimum_consecutive_trapped_common_events"
    ]
    selected_amplitude: str | None = None
    terminal_classification: str | None = None
    campaign_started = wall_time.monotonic()

    while amplitude_index < len(amplitudes):
        amplitude = amplitudes[amplitude_index]
        if any(item.amplitude != amplitude for item in members.values()):
            raise RuntimeError("checkpoint members belong to another amplitude")
        candidate_stop: dict[str, Any] | None = None
        while event_index < final_event_index:
            target_index = event_index + 1
            target_time = target_index * output_interval
            snapshots = {key: item.snapshot() for key, item in members.items()}
            active_member: Proto8RunMember | None = None

            def rejected_trial_sink(evidence: Mapping[str, Any]) -> None:
                inherited._append_fsync(event_log, inherited._canonical_line(evidence))

            try:
                for active_member in members.values():
                    active_member.advance_to(
                        target_time,
                        retry_factor=retry_factor,
                        maximum_CFL_retries=maximum_CFL_retries,
                        maximum_source_retries=maximum_source_retries,
                        minimum_step_size=minimum_step,
                        rejected_trial_sink=rejected_trial_sink,
                    )
                for member in members.values():
                    member.append_common_event()
                assessment = _event_assessment(plan, members, target_time)
            except Proto7TerminalStop as stop:
                evidence = asdict(stop.evidence)
                restore_common_event_snapshots(members, snapshots)
                restored = proto7_runner._restored_common_boundary(members, snapshots)
                if not restored:
                    raise RuntimeError("PROTO8 terminal rollback did not restore common event")
                candidate_stop = {
                    "classification": "scientific_universal_runtime_stop",
                    "reason": stop.reason,
                    "failures": list(stop.failures),
                    "complete_failure_evidence": evidence,
                    "member": None if active_member is None else active_member.key,
                    "method": None if active_member is None else active_member.method_label,
                    "point_count": None if active_member is None else active_member.point_count,
                    "target_common_event_index": target_index,
                    "last_fully_completed_common_event_index": event_index,
                    "cross_member_state_rolled_back": True,
                    "cross_member_rollback_verified": True,
                }
                break
            except Proto8SourceRetryExhausted as stop:
                restore_common_event_snapshots(members, snapshots)
                restored = proto7_runner._restored_common_boundary(members, snapshots)
                if not restored:
                    raise RuntimeError("PROTO8 exhaustion rollback did not restore common event")
                candidate_stop = {
                    "classification": "scientific_source_retry_exhausted",
                    "reason": stop.reason,
                    "complete_failure_evidence": stop.evidence,
                    "member": None if active_member is None else active_member.key,
                    "method": None if active_member is None else active_member.method_label,
                    "point_count": None if active_member is None else active_member.point_count,
                    "target_common_event_index": target_index,
                    "last_fully_completed_common_event_index": event_index,
                    "cross_member_state_rolled_back": True,
                    "cross_member_rollback_verified": True,
                }
                break
            except GR0RuntimeStop as stop:
                restore_common_event_snapshots(members, snapshots)
                candidate_stop = {
                    "classification": "invalid_implementation_or_nonconverged_run",
                    "reason": "terminal_stop_without_PROTO7_complete_evidence",
                    "inherited_reason": stop.reason,
                    "failures": list(stop.failures),
                    "target_common_event_index": target_index,
                    "last_fully_completed_common_event_index": event_index,
                    "cross_member_state_rolled_back": True,
                }
                break
            except BaseException as error:
                if isinstance(error, KeyboardInterrupt):
                    raise
                restore_common_event_snapshots(members, snapshots)
                reason = inherited._source_failure_reason(error)
                candidate_stop = {
                    "classification": (
                        "scientific_source_or_measurement_stop"
                        if reason != "invalid_implementation_or_nonconverged_run"
                        else reason
                    ),
                    "reason": reason,
                    "error_type": type(error).__name__,
                    "error_message": str(error),
                    "target_common_event_index": target_index,
                    "last_fully_completed_common_event_index": event_index,
                    "cross_member_state_rolled_back": True,
                }
                break

            event_index = target_index
            consecutive = (
                consecutive + 1 if assessment["qualified_trapped_common_event"] else 0
            )
            inherited._append_fsync(
                event_log,
                inherited._canonical_line(
                    {
                        "event_type": "common_event",
                        "amplitude": amplitude,
                        "event_index": event_index,
                        "assessment": assessment,
                        "consecutive_qualified_trapped_common_events": consecutive,
                    }
                ),
            )
            if event_index % checkpoint_every_events == 0:
                checkpoint = _checkpoint_metadata(
                    manifest=manifest,
                    amplitude=amplitude,
                    amplitude_index=amplitude_index,
                    event_index=event_index,
                    consecutive_qualified_events=consecutive,
                    amplitude_records=amplitude_records,
                    event_log_path=event_log,
                    members=members,
                    terminal=False,
                )
                proto6_runner._write_checkpoint(
                    checkpoint_path,
                    metadata=checkpoint,
                    members=members,
                    event_log_path=event_log,
                )
            if event_index % 8 == 0 or consecutive > 0:
                trapped = assessment["trapped_assessment"]
                print(
                    f"amplitude={amplitude} event={event_index}/{final_event_index} "
                    f"t={target_time:.6g} trapped_score="
                    f"{trapped.primary_trapped_scores[-1]:.6g}/"
                    f"{trapped.comparator_trapped_scores[-1]:.6g} "
                    f"qualified_run={consecutive}/{required_consecutive}",
                    flush=True,
                )
            if consecutive >= required_consecutive:
                selected_amplitude = amplitude
                terminal_classification = "GR0_calibration_completed_selected_amplitude"
                break

            common_pass = (
                assessment["primary_common_event"].admission_passed
                and assessment["comparator_common_event"].admission_passed
            )
            temporal = assessment["temporal_spectral_admission"]
            if not common_pass or (
                temporal["available"] and not temporal["admission_passed"]
            ):
                candidate_stop = {
                    "classification": "common_event_constraint_or_spectral_stop",
                    "reason": (
                        "common_event_constraint_or_spatial_spectral_admission"
                        if not common_pass
                        else "causal_past_temporal_spectral_admission"
                    ),
                    "event_index": event_index,
                    "coordinate_time": target_time,
                }
                break

        amplitude_records.append(
            {
                "amplitude": amplitude,
                "last_completed_common_event_index": event_index,
                "last_completed_coordinate_time": event_index * output_interval,
                "consecutive_qualified_trapped_common_events": consecutive,
                "eligible": selected_amplitude == amplitude,
                "stop": candidate_stop,
                "member_source_retry_counts": {
                    key: item.source_retry_count for key, item in members.items()
                },
                "member_final_state_hashes": {
                    key: array_content_sha256(item.state.u, item.state.p, item.state.q)
                    for key, item in members.items()
                },
            }
        )
        if selected_amplitude is not None:
            break
        if (
            candidate_stop is not None
            and candidate_stop["classification"]
            == "invalid_implementation_or_nonconverged_run"
        ):
            terminal_classification = "invalid_implementation_or_nonconverged_run"
            break
        amplitude_index += 1
        if amplitude_index >= len(amplitudes):
            terminal_classification = "calibration_failed_no_eligible_GR0_case"
            break
        amplitude = amplitudes[amplitude_index]
        event_index = 0
        consecutive = 0
        members = _build_members(plan, authorization, amplitude)
        initial = _event_assessment(plan, members, 0.0)
        inherited._append_fsync(
            event_log,
            inherited._canonical_line(
                {
                    "event_type": "common_event",
                    "amplitude": amplitude,
                    "event_index": 0,
                    "assessment": initial,
                }
            ),
        )
        checkpoint = _checkpoint_metadata(
            manifest=manifest,
            amplitude=amplitude,
            amplitude_index=amplitude_index,
            event_index=0,
            consecutive_qualified_events=0,
            amplitude_records=amplitude_records,
            event_log_path=event_log,
            members=members,
            terminal=False,
        )
        proto6_runner._write_checkpoint(
            checkpoint_path,
            metadata=checkpoint,
            members=members,
            event_log_path=event_log,
        )

    if terminal_classification is None:
        terminal_classification = "calibration_failed_no_eligible_GR0_case"
    result = {
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "campaign_id": manifest["campaign_id"],
        "classification": terminal_classification,
        "selected_amplitude": selected_amplitude,
        "amplitude_records": amplitude_records,
        "event_log": event_log.relative_to(REPOSITORY).as_posix(),
        "event_log_sha256": inherited._sha(event_log),
        "elapsed_wall_seconds": wall_time.monotonic() - campaign_started,
        "resume_recovery": resume_recovery,
        "cross_member_common_event_rollback_enabled": True,
        "complete_PROTO7_failure_observability_enabled": True,
        "PROTO8_common_event_ownership_enabled": True,
        "FGCQR_outcome_read": False,
        "SGBL_outcome_read": False,
        "holdout_execution_authorized": False,
        "retained_EFT_evolution_authorized": False,
        "mechanism_question_answered": False,
    }
    inherited._write_exclusive(result_path, inherited._canonical(result))
    checkpoint = _checkpoint_metadata(
        manifest=manifest,
        amplitude=amplitudes[min(amplitude_index, len(amplitudes) - 1)],
        amplitude_index=min(amplitude_index, len(amplitudes) - 1),
        event_index=event_index,
        consecutive_qualified_events=consecutive,
        amplitude_records=amplitude_records,
        event_log_path=event_log,
        members=members,
        terminal=True,
    )
    proto6_runner._write_checkpoint(
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
                    "campaign_id": authorization["artifact_payload"]["frozen_run_plan"]["campaign_id"],
                    "ordered_amplitudes": plan["candidate_selection"]["ordered_amplitudes"],
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
