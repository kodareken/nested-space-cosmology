#!/usr/bin/env python3
"""Execute the HLT4-authorized fresh PROTO6 GR-0 calibration campaign.

The large, unchanged campaign mechanics are imported from the immutable
PROTO5 runner.  This successor owns only the PROTO6 delta: diagnostic exposure
of the same affine source root, accepted-state prechecks, source-only internal
trial rejection, exact rollback, a public retry ledger, and the fresh PROTO6
namespace.  It remains a GR-0 calibration runner and cannot open SGB-L or
FGC-QR outcomes.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, replace
from pathlib import Path
import json
import time as wall_time
from typing import Any, Callable, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_hlt4_mon4 as hlt4  # noqa: E402
from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    project_gr0_semidiscrete_state,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    GR0GridInitialData,
    construct_gr0_grid_initial_data,
    make_gr0_center_boundary_projector,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    CFLRetryRequired,
    GR0RuntimeMonitorState,
    GR0RuntimeStageTransaction,
    GR0RuntimeStop,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.proto6_runtime import (  # noqa: E402
    Proto6GR0EvolutionOperator,
    attempt_proto6_step,
    restore_common_event_snapshots,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = inherited.Q
DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-cal3-run1.toml"
DEFAULT_AUTHORIZATION = REPOSITORY / "results/fgc-1-hlt4-mon4.json"
RUNNER_ID = "FGC-1-CAL3-RUN1-RUNNER"


class Proto6SourceRetryExhausted(RuntimeError):
    """Frozen source-retry budget or minimum step was reached."""

    reason = "newton_residual_limit"

    def __init__(self, message: str, evidence: Mapping[str, Any]) -> None:
        super().__init__(message)
        self.evidence = dict(evidence)


@dataclass
class Proto6RunMember:
    amplitude: str
    method_label: str
    integrator_id: str
    spatial_order: int
    point_count: int
    input_hash: str
    initial: GR0GridInitialData
    state: EvolutionState
    operator: Proto6GR0EvolutionOperator
    projector: Any
    transaction: GR0RuntimeStageTransaction
    tracers: inherited.NormalFlowTracers
    time: float = 0.0
    step_index: int = 0
    transaction_serial: int = 0
    CFL_retry_count: int = 0
    source_retry_count: int = 0

    @property
    def key(self) -> str:
        return f"{self.method_label}-{self.point_count}"

    def snapshot(self) -> dict[str, Any]:
        return {
            "state": self.state,
            "time": self.time,
            "step_index": self.step_index,
            "transaction_serial": self.transaction_serial,
            "CFL_retry_count": self.CFL_retry_count,
            "source_retry_count": self.source_retry_count,
            "monitor_state": self.transaction.state,
            "causal_state": self.transaction.causal_state,
            "tracer_positions": self.tracers.positions.copy(),
            "tracer_proper_times": self.tracers.proper_times.copy(),
            "event_proper_times": [item.copy() for item in self.tracers.event_proper_times],
            "event_fields": [item.copy() for item in self.tracers.event_fields],
        }

    def restore(self, snapshot: Mapping[str, Any]) -> None:
        self.state = snapshot["state"]
        self.time = float(snapshot["time"])
        self.step_index = int(snapshot["step_index"])
        self.transaction_serial = int(snapshot["transaction_serial"])
        self.CFL_retry_count = int(snapshot["CFL_retry_count"])
        self.source_retry_count = int(snapshot["source_retry_count"])
        self.transaction.state = snapshot["monitor_state"]
        self.transaction.causal_state = snapshot["causal_state"]
        self.tracers.positions = snapshot["tracer_positions"].copy()
        self.tracers.proper_times = snapshot["tracer_proper_times"].copy()
        self.tracers.event_proper_times = [
            item.copy() for item in snapshot["event_proper_times"]
        ]
        self.tracers.event_fields = [
            item.copy() for item in snapshot["event_fields"]
        ]

    def advance_to(
        self,
        target_time: float,
        *,
        retry_factor: float,
        maximum_CFL_retries: int,
        maximum_source_retries: int,
        minimum_step_size: float,
        rejected_trial_sink: Callable[[Mapping[str, Any]], None],
    ) -> None:
        tolerance = 4096.0 * np.finfo(np.float64).eps * max(1.0, target_time)
        while self.time < target_time - tolerance:
            previous_speed = max(
                self.transaction.causal_state.previous_speed_upper,
                np.finfo(np.float64).tiny,
            )
            proposed_size = min(
                target_time - self.time,
                self.transaction.cfl_maximum
                * self.initial.grid.spacing
                / previous_speed,
            )
            cfl_retries = 0
            source_retries = 0
            last_source_evidence: dict[str, Any] | None = None
            while True:
                if proposed_size < minimum_step_size:
                    if last_source_evidence is not None:
                        raise Proto6SourceRetryExhausted(
                            "PROTO6 source retry fell below the frozen minimum step",
                            last_source_evidence,
                        )
                    raise RuntimeError("adaptive step fell below the frozen minimum")
                accepted_hash = array_content_sha256(
                    self.state.u, self.state.p, self.state.q
                )
                monitor_before = self.transaction.state
                causal_before = self.transaction.causal_state

                def preview(proposal):
                    return self.tracers.preview_advance(
                        old_state=self.state,
                        new_state=proposal.candidate_state,
                        coordinates=self.initial.grid.coordinates,
                        step_size=proposal.final_time - self.time,
                    )

                try:
                    attempt = attempt_proto6_step(
                        method=self.integrator_id,
                        time=self.time,
                        step_size=proposed_size,
                        state=self.state,
                        rhs=self.operator,
                        projector=self.projector,
                        transaction=self.transaction,
                        previous_step_index=self.step_index,
                        previous_transaction_serial=self.transaction_serial,
                        preaccept=preview,
                    )
                except CFLRetryRequired:
                    cfl_retries += 1
                    self.CFL_retry_count += 1
                    if cfl_retries > maximum_CFL_retries:
                        raise RuntimeError("frozen maximum CFL retries exceeded")
                    proposed_size *= retry_factor
                    continue

                if attempt.retry is not None:
                    if (
                        array_content_sha256(self.state.u, self.state.p, self.state.q)
                        != accepted_hash
                        or self.transaction.state != monitor_before
                        or self.transaction.causal_state != causal_before
                    ):
                        raise RuntimeError("PROTO6 rejected trial changed accepted state")
                    source_retries += 1
                    self.source_retry_count += 1
                    last_source_evidence = {
                        **asdict(attempt.retry),
                        "event_type": "rejected_internal_source_trial",
                        "amplitude": self.amplitude,
                        "member": self.key,
                        "retry_count_for_accepted_step": source_retries,
                        "cumulative_member_source_retry_count": self.source_retry_count,
                        "target_common_event_time": target_time,
                    }
                    rejected_trial_sink(last_source_evidence)
                    if source_retries > maximum_source_retries:
                        raise Proto6SourceRetryExhausted(
                            "PROTO6 source retry exceeded the frozen retry count",
                            last_source_evidence,
                        )
                    proposed_size *= retry_factor
                    continue

                accepted = attempt.accepted
                if accepted is None:
                    raise RuntimeError("PROTO6 attempt produced no accepted step")
                self.tracers.commit_advance(*attempt.preaccept_payload)
                self.state = accepted.state
                self.time = accepted.time
                self.step_index = accepted.step_index
                self.transaction_serial = accepted.transaction_serial
                break
        if np.float64(self.time).tobytes() != np.float64(target_time).tobytes():
            if abs(self.time - target_time) > tolerance:
                raise RuntimeError("member did not land on the common event")
            self.time = target_time

    def append_common_event(self) -> None:
        self.tracers.append_common_event(
            self.state,
            self.initial.grid.coordinates,
        )


def _validate_authorization(
    plan_path: Path,
    authorization_path: Path,
    *,
    require_fresh_namespace: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = hlt4._validate_run_plan(plan_path)
    authorization = hlt4.load_canonical_result(authorization_path)
    gates = authorization.get("gate_status", {})
    if (
        authorization.get("artifact_id") != hlt4.ARTIFACT_ID
        or gates.get("PROTO6_successor_runtime_compositor_implemented") is not True
        or gates.get("PROTO6_fresh_GR0_dynamic_calibration_authorized") is not True
        or gates.get("FGCQR_holdout_execution_authorized") is not False
    ):
        raise ValueError("HLT4 authorization is absent, incomplete, or over-broad")
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if frozen["run_plan_sha256"] != inherited._sha(plan_path):
        raise ValueError("CAL3 run-plan hash differs from HLT4")
    for relative, expected in authorization["implementation_sha256"].items():
        current = REPOSITORY / relative
        if not current.is_file() or inherited._sha(current) != expected:
            raise ValueError(f"authorized implementation drifted: {relative}")
    if require_fresh_namespace:
        reproduced = hlt4.record(
            REPOSITORY / "configs/fgc/fgc-1-hlt4-mon4.toml"
        )
        if inherited._serial(reproduced) != authorization:
            raise ValueError("HLT4 authorization does not reproduce before launch")
    return plan, authorization


def _member_inputs(
    authorization: Mapping[str, Any], amplitude: str
) -> dict[tuple[str, int], Mapping[str, Any]]:
    records = authorization["artifact_payload"]["frozen_run_inputs"]
    selected = {
        (item["method"], item["point_count"]): item
        for item in records
        if item["amplitude"] == amplitude
    }
    if len(selected) != 6:
        raise ValueError("HLT4 does not contain six inputs for this amplitude")
    return selected


def _build_members(
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
    amplitude: str,
) -> dict[str, Proto6RunMember]:
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    thresholds = plan["universal_thresholds"]
    frozen = _member_inputs(authorization, amplitude)
    parameters = PulseParameters(
        chi_amplitude=float(Q(amplitude)),
        center=float(Q(physical["chi_center"])),
        half_width=float(Q(physical["chi_half_width"])),
        phi_amplitude=float(Q(physical["phi_seed_amplitude"])),
        planck_mass=float(Q(physical["planck_mass"])),
        scalar_mass=float(Q(physical["scalar_mass"])),
        quartic_coupling=float(Q(physical["quartic_coupling"])),
    )
    answer: dict[str, Proto6RunMember] = {}
    for method in (plan["primary_method"], plan["comparator_method"]):
        method_label = method["method_label"]
        order = method["spatial_order"]
        for point_count in numerics["resolutions"]:
            input_record = frozen[(method_label, point_count)]
            initial = construct_gr0_grid_initial_data(
                parameters,
                point_count=point_count,
                outer_radius=float(Q(physical["outer_radius"])),
                constraint_method=method["constraint_solve_method"],
                diagnostic_spatial_order=order,
            )
            state = project_gr0_semidiscrete_state(initial, spatial_order=order)
            if array_content_sha256(state.u, state.p, state.q) != input_record[
                "projected_state_sha256"
            ]:
                raise ValueError("runtime projected state differs from HLT4 input")
            projected_initial = replace(initial, state=state)
            operator = Proto6GR0EvolutionOperator(
                initial.grid,
                spatial_order=order,
                ko_dissipation=float(Q(numerics["ko_dissipation"])),
                raw_tolerance=float(Q(thresholds["source_residual_infinity_max"])),
                kinetic_condition_maximum=float(
                    Q(thresholds["kinetic_condition_number_max"])
                ),
            )
            initial_rhs = operator(0.0, state)
            derivative = SBPFirstDerivative(initial.grid, order)
            transaction = GR0RuntimeStageTransaction(
                thresholds=GR0UniversalThresholds(
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
                minimum=float(Q(numerics["normal_flow_tracer_radius_minimum"])),
                maximum=float(Q(numerics["normal_flow_tracer_radius_maximum"])),
                spacing=float(Q(numerics["normal_flow_tracer_spacing"])),
                state=state,
                coordinates=initial.grid.coordinates,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                outer_radius=float(Q(physical["outer_radius"])),
            )
            member = Proto6RunMember(
                amplitude=amplitude,
                method_label=method_label,
                integrator_id=method["integrator_id"],
                spatial_order=order,
                point_count=point_count,
                input_hash=input_record["expanded_run_config_sha256"],
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
                raise RuntimeError("duplicate PROTO6 GR-0 run member")
            answer[member.key] = member
    if len(answer) != 6:
        raise RuntimeError("PROTO6 GR-0 amplitude requires six run members")
    return answer


def _member_metadata(member: Proto6RunMember) -> dict[str, Any]:
    return {
        "time": member.time,
        "step_index": member.step_index,
        "transaction_serial": member.transaction_serial,
        "CFL_retry_count": member.CFL_retry_count,
        "source_retry_count": member.source_retry_count,
        "runtime_monitor_state": asdict(member.transaction.state),
        "causal_state": asdict(member.transaction.causal_state),
        "input_hash": member.input_hash,
    }


def _checkpoint_metadata(
    *,
    manifest: Mapping[str, Any],
    amplitude: str,
    amplitude_index: int,
    event_index: int,
    consecutive_qualified_events: int,
    amplitude_records: Sequence[Mapping[str, Any]],
    event_log_path: Path,
    members: Mapping[str, Proto6RunMember],
    terminal: bool,
) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "campaign_id": manifest["campaign_id"],
        "amplitude": amplitude,
        "amplitude_index": amplitude_index,
        "completed_common_event_index": event_index,
        "consecutive_qualified_events": consecutive_qualified_events,
        "completed_amplitude_records": list(amplitude_records),
        "event_log": inherited._event_log_identity(
            inherited._event_log_bytes(event_log_path)
        ),
        "terminal": terminal,
        "members": {key: _member_metadata(item) for key, item in members.items()},
    }


def _restore_checkpoint(
    path: Path,
    *,
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Proto6RunMember], bytes]:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(bytes(archive["metadata_utf8"]).decode("utf-8"))
        event_log_payload = bytes(archive["event_log_utf8"])
        if metadata.get("event_log") != inherited._event_log_identity(event_log_payload):
            raise ValueError("checkpoint event-log evidence is inconsistent")
        amplitude = metadata["amplitude"]
        members = _build_members(plan, authorization, amplitude)
        for key, member in members.items():
            prefix = key.replace("-", "_")
            stored = metadata["members"][key]
            if stored["input_hash"] != member.input_hash:
                raise ValueError("checkpoint input hash differs from HLT4")
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
    return metadata, members, event_log_payload


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
        raise ValueError("resume manifest differs from the authorized PROTO6 campaign")


def _write_checkpoint(
    path: Path,
    *,
    metadata: Mapping[str, Any],
    members: Mapping[str, Proto6RunMember],
    event_log_path: Path,
) -> None:
    inherited._write_checkpoint(
        path,
        metadata=metadata,
        members=members,
        event_log_path=event_log_path,
    )


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
        _validate_resume_manifest(
            manifest, plan_path, authorization_path, authorization
        )
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
            event_log,
            checkpoint_event_log,
        )
    else:
        if output_root.exists():
            raise ValueError("fresh PROTO6 calibration output root already exists")
        output_root.mkdir(parents=True, exist_ok=False)
        manifest = _manifest(plan_path, authorization_path, authorization)
        inherited._write_exclusive(manifest_path, inherited._canonical(manifest))
        amplitude_index = 0
        event_index = 0
        consecutive = 0
        amplitude_records: list[dict[str, Any]] = []
        amplitude = plan["candidate_selection"]["ordered_amplitudes"][0]
        members = _build_members(plan, authorization, amplitude)
        initial = inherited._event_assessment(plan, members, 0.0)
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
        _write_checkpoint(
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

            def rejected_trial_sink(evidence: Mapping[str, Any]) -> None:
                inherited._append_fsync(
                    event_log,
                    inherited._canonical_line(evidence),
                )

            try:
                for member in members.values():
                    member.advance_to(
                        target_time,
                        retry_factor=retry_factor,
                        maximum_CFL_retries=maximum_CFL_retries,
                        maximum_source_retries=maximum_source_retries,
                        minimum_step_size=minimum_step,
                        rejected_trial_sink=rejected_trial_sink,
                    )
                for member in members.values():
                    member.append_common_event()
                assessment = inherited._event_assessment(plan, members, target_time)
            except GR0RuntimeStop as stop:
                restore_common_event_snapshots(members, snapshots)
                candidate_stop = {
                    "classification": "scientific_universal_runtime_stop",
                    "reason": stop.reason,
                    "failures": list(stop.failures),
                    "target_common_event_index": target_index,
                    "last_fully_completed_common_event_index": event_index,
                    "cross_member_state_rolled_back": True,
                }
                break
            except Proto6SourceRetryExhausted as stop:
                restore_common_event_snapshots(members, snapshots)
                candidate_stop = {
                    "classification": "scientific_source_retry_exhausted",
                    "reason": stop.reason,
                    "evidence": stop.evidence,
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
            consecutive = consecutive + 1 if assessment[
                "qualified_trapped_common_event"
            ] else 0
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
                _write_checkpoint(
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
            and candidate_stop["reason"] == "invalid_implementation_or_nonconverged_run"
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
        initial = inherited._event_assessment(plan, members, 0.0)
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
        _write_checkpoint(
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
    _write_checkpoint(
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
                    "campaign_id": authorization["artifact_payload"][
                        "frozen_run_plan"
                    ]["campaign_id"],
                    "ordered_amplitudes": plan["candidate_selection"][
                        "ordered_amplitudes"
                    ],
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
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
