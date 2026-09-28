#!/usr/bin/env python3
"""Execute the HLT3-authorized fresh PROTO5 GR-0 calibration campaign.

This is an evidence-producing runner, not a reproducer for a tracked result.
It writes only beneath the ignored PROTO5 calibration namespace after
verifying the canonical HLT3 authorization and an otherwise clean tracked
worktree.  The resulting raw campaign must later be reduced by a separate
CAL2 certificate before any holdout manifest can be resolved.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass, replace
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
import os
from pathlib import Path
import subprocess
import sys
import time as wall_time
from typing import Any, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_hlt3_mon3 as hlt3  # noqa: E402
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    project_gr0_semidiscrete_state,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    GR0EvolutionOperator,
    GR0GridInitialData,
    construct_gr0_grid_initial_data,
    make_gr0_center_boundary_projector,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    accept_step,
    array_content_sha256,
    propose_step,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_SPECTRAL_FIELD_ORDER,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    CFLRetryRequired,
    GR0RuntimeMonitorState,
    GR0RuntimeStageTransaction,
    GR0RuntimeStop,
    GR0UniversalThresholds,
    proto5_calibration_trapped_assessment,
    proto5_event_log_recovery_plan,
    proto5_gr0_common_event,
    proto5_temporal_spectral_admission,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = Fraction
DEFAULT_PLAN = REPOSITORY / "configs/fgc/fgc-1-cal2-run1.toml"
DEFAULT_AUTHORIZATION = REPOSITORY / "results/fgc-1-hlt3-mon3.json"
RUNNER_ID = "FGC-1-CAL2-RUN1-RUNNER"


def _serial(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _serial(value.tolist())
    if isinstance(value, np.generic):
        return _serial(value.item())
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if hasattr(value, "__dataclass_fields__"):
        return _serial(asdict(value))
    if isinstance(value, float):
        if not isfinite(value):
            raise ValueError("nonfinite value cannot enter run evidence")
        return value
    if value is None or isinstance(value, (str, bool, int)):
        return value
    raise TypeError(f"unsupported run-evidence type {type(value).__name__}")


def _canonical_line(value: Any) -> bytes:
    return (
        json.dumps(
            _serial(value),
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _canonical(value: Any) -> bytes:
    return (
        json.dumps(
            _serial(value),
            sort_keys=True,
            indent=2,
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _write_exclusive(path: Path, payload: bytes) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
    except BaseException:
        try:
            path.unlink(missing_ok=True)
        finally:
            raise


def _atomic_replace(path: Path, payload: bytes) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    _write_exclusive(temporary, payload)
    os.replace(temporary, path)


def _append_fsync(path: Path, payload: bytes) -> None:
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    with os.fdopen(descriptor, "ab") as handle:
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())


def _event_log_identity(payload: bytes) -> dict[str, Any]:
    if payload and not payload.endswith(b"\n"):
        raise ValueError("event log must end at a complete canonical JSON line")
    return {
        "sha256": sha256(payload).hexdigest(),
        "byte_count": len(payload),
        "line_count": payload.count(b"\n"),
    }


def _event_log_bytes(path: Path) -> bytes:
    if not path.is_file():
        raise ValueError("campaign event log is absent")
    payload = path.read_bytes()
    _event_log_identity(payload)
    return payload


def _write_exclusive_or_match(path: Path, payload: bytes) -> None:
    try:
        _write_exclusive(path, payload)
    except FileExistsError:
        if path.read_bytes() != payload:
            raise ValueError(f"existing recovery evidence differs: {path}")


def _recover_event_log_from_checkpoint(
    path: Path,
    checkpoint_payload: bytes,
) -> dict[str, Any]:
    """Restore the append-only view to the last atomic checkpoint.

    The checkpoint owns accepted state. A process can be interrupted after an
    event line is fsynced but before the new checkpoint is atomically installed.
    Preserve that uncommitted tail in a hash-named recovery sidecar, then restore
    the exact checkpoint-owned bytes. Divergent histories fail closed.
    """

    checkpoint_identity = _event_log_identity(checkpoint_payload)
    current = _event_log_bytes(path) if path.exists() else None
    plan = proto5_event_log_recovery_plan(current, checkpoint_payload)
    if plan.action == "event_log_already_matches_checkpoint":
        return {
            "action": plan.action,
            "checkpoint_event_log": checkpoint_identity,
        }
    if plan.orphaned_tail is not None:
        orphaned_tail = plan.orphaned_tail
        recovery = path.with_name(
            "orphaned-event-tail-"
            f"{sha256(orphaned_tail).hexdigest()}.jsonl"
        )
        _write_exclusive_or_match(recovery, orphaned_tail)
        _atomic_replace(path, plan.restored_payload)
        return {
            "action": plan.action,
            "checkpoint_event_log": checkpoint_identity,
            "orphaned_tail_path": recovery.name,
            "orphaned_tail": _event_log_identity(orphaned_tail),
        }
    _atomic_replace(path, plan.restored_payload)
    return {
        "action": plan.action,
        "checkpoint_event_log": checkpoint_identity,
    }


def _git_head() -> str:
    return subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPOSITORY,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()


def _require_clean_tracked_worktree() -> None:
    status = subprocess.run(
        ["git", "status", "--porcelain", "--untracked-files=no"],
        cwd=REPOSITORY,
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    if status.strip():
        raise ValueError(
            "fresh calibration requires a clean tracked worktree so HLT3 hashes own the run"
        )


def _validate_authorization(
    plan_path: Path,
    authorization_path: Path,
    *,
    require_fresh_namespace: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    plan = hlt3._validate_run_plan(plan_path)  # same exact owner as HLT3
    authorization = hlt3.load_canonical_result(authorization_path)
    gates = authorization.get("gate_status", {})
    if (
        authorization.get("artifact_id") != hlt3.ARTIFACT_ID
        or gates.get("PROTO5_successor_runtime_compositor_implemented") is not True
        or gates.get("PROTO5_fresh_GR0_dynamic_calibration_authorized") is not True
        or gates.get("FGCQR_holdout_execution_authorized") is not False
    ):
        raise ValueError("HLT3 authorization is absent, incomplete, or over-broad")
    frozen = authorization["artifact_payload"]["frozen_run_plan"]
    if frozen["run_plan_sha256"] != _sha(plan_path):
        raise ValueError("CAL2 run-plan hash differs from HLT3")
    for relative, expected in authorization["implementation_sha256"].items():
        current = REPOSITORY / relative
        if not current.is_file() or _sha(current) != expected:
            raise ValueError(f"authorized implementation drifted: {relative}")
    if require_fresh_namespace:
        freshly_reproduced = hlt3.record(
            REPOSITORY / "configs/fgc/fgc-1-hlt3-mon3.toml"
        )
        if _serial(freshly_reproduced) != authorization:
            raise ValueError("HLT3 authorization does not reproduce before launch")
    return plan, authorization


def _interpolate_columns(
    coordinates: np.ndarray,
    values: np.ndarray,
    positions: np.ndarray,
) -> np.ndarray:
    if (
        coordinates.ndim != 1
        or values.ndim != 2
        or values.shape[0] != coordinates.size
        or positions.ndim != 1
    ):
        raise ValueError("tracer interpolation shapes differ")
    return np.column_stack(
        [np.interp(positions, coordinates, values[:, field]) for field in range(values.shape[1])]
    )


@dataclass
class NormalFlowTracers:
    labels: np.ndarray
    positions: np.ndarray
    proper_times: np.ndarray
    event_proper_times: list[np.ndarray]
    event_fields: list[np.ndarray]
    cutoff: float
    outer_radius: float

    @classmethod
    def create(
        cls,
        *,
        minimum: float,
        maximum: float,
        spacing: float,
        state: EvolutionState,
        coordinates: np.ndarray,
        cutoff: float,
        outer_radius: float,
    ) -> "NormalFlowTracers":
        count_float = (maximum - minimum) / spacing
        count = int(round(count_float))
        if abs(count_float - count) > 1.0e-12:
            raise ValueError("tracer interval is not aligned")
        labels = np.linspace(minimum, maximum, count + 1, dtype=np.float64)
        answer = cls(
            labels=labels,
            positions=labels.copy(),
            proper_times=np.zeros(labels.size, dtype=np.float64),
            event_proper_times=[],
            event_fields=[],
            cutoff=cutoff,
            outer_radius=outer_radius,
        )
        answer.append_common_event(state, coordinates)
        return answer

    def _metric_samples(
        self,
        state: EvolutionState,
        coordinates: np.ndarray,
        positions: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        fields = _interpolate_columns(coordinates, state.u, positions)
        lapse = fields[:, 0]
        shift = fields[:, 1]
        if np.any(lapse <= 0.0):
            raise ValueError("normal-flow tracer encountered nonpositive lapse")
        return lapse, shift

    def preview_advance(
        self,
        *,
        old_state: EvolutionState,
        new_state: EvolutionState,
        coordinates: np.ndarray,
        step_size: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        old_lapse, old_shift = self._metric_samples(
            old_state, coordinates, self.positions
        )
        old_velocity = -old_shift
        predictor = self.positions + step_size * old_velocity
        if np.any(predictor <= 0.0) or np.any(predictor >= self.outer_radius):
            raise ValueError("normal-flow tracer predictor left the numerical domain")
        new_lapse, new_shift = self._metric_samples(
            new_state, coordinates, predictor
        )
        new_velocity = -new_shift
        positions = self.positions + 0.5 * step_size * (
            old_velocity + new_velocity
        )
        proper = self.proper_times + 0.5 * step_size * (
            old_lapse + new_lapse
        )
        if (
            np.any(positions <= 0.0)
            or np.any(positions >= self.outer_radius)
            or np.any(proper <= self.proper_times)
        ):
            raise ValueError("normal-flow tracer update lost its timelike domain")
        return positions, proper

    def commit_advance(
        self,
        positions: np.ndarray,
        proper_times: np.ndarray,
    ) -> None:
        if (
            positions.shape != self.positions.shape
            or proper_times.shape != self.proper_times.shape
        ):
            raise ValueError("normal-flow tracer commit shapes differ")
        self.positions = np.asarray(positions, dtype=np.float64).copy()
        self.proper_times = np.asarray(proper_times, dtype=np.float64).copy()

    def field_sample(
        self,
        state: EvolutionState,
        coordinates: np.ndarray,
    ) -> np.ndarray:
        fields = _interpolate_columns(coordinates, state.u, self.positions)
        radius = self.positions
        answer = np.column_stack(
            (
                fields[:, 0] - 1.0,
                fields[:, 1] / radius,
                fields[:, 2] - 1.0,
                fields[:, 3] / radius - 1.0,
                fields[:, 4] / self.cutoff,
                fields[:, 5] / self.cutoff,
            )
        )
        if not np.all(np.isfinite(answer)):
            raise ValueError("normal-flow tracer field sample became nonfinite")
        return answer

    def append_common_event(
        self,
        state: EvolutionState,
        coordinates: np.ndarray,
    ) -> None:
        self.event_proper_times.append(self.proper_times.copy())
        self.event_fields.append(self.field_sample(state, coordinates))


@dataclass
class RunMember:
    amplitude: str
    method_label: str
    integrator_id: str
    spatial_order: int
    point_count: int
    input_hash: str
    initial: GR0GridInitialData
    state: EvolutionState
    operator: GR0EvolutionOperator
    projector: Any
    transaction: GR0RuntimeStageTransaction
    tracers: NormalFlowTracers
    time: float = 0.0
    step_index: int = 0
    transaction_serial: int = 0
    CFL_retry_count: int = 0

    @property
    def key(self) -> str:
        return f"{self.method_label}-{self.point_count}"

    def advance_to(
        self,
        target_time: float,
        *,
        retry_factor: float,
        maximum_retries: int,
        minimum_step_size: float,
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
            retries = 0
            while True:
                if proposed_size < minimum_step_size:
                    raise RuntimeError("adaptive step fell below the frozen minimum")
                old_state = self.state
                try:
                    proposal = propose_step(
                        method=self.integrator_id,
                        time=self.time,
                        step_size=proposed_size,
                        state=self.state,
                        rhs=self.operator,
                        projector=self.projector,
                    )
                    tracer_candidate = self.tracers.preview_advance(
                        old_state=old_state,
                        new_state=proposal.candidate_state,
                        coordinates=self.initial.grid.coordinates,
                        step_size=proposal.final_time - self.time,
                    )
                    accepted = accept_step(
                        proposal,
                        previous_step_index=self.step_index,
                        previous_transaction_serial=self.transaction_serial,
                        guards=(self.transaction,),
                    )
                except CFLRetryRequired:
                    retries += 1
                    self.CFL_retry_count += 1
                    if retries > maximum_retries:
                        raise RuntimeError("frozen maximum CFL retries exceeded")
                    proposed_size *= retry_factor
                    continue
                break
            self.tracers.commit_advance(*tracer_candidate)
            self.state = accepted.state
            self.time = accepted.time
            self.step_index = accepted.step_index
            self.transaction_serial = accepted.transaction_serial
        if np.float64(self.time).tobytes() != np.float64(target_time).tobytes():
            if abs(self.time - target_time) > tolerance:
                raise RuntimeError("member did not land on the common event")
            self.time = target_time

    def append_common_event(self) -> None:
        self.tracers.append_common_event(
            self.state,
            self.initial.grid.coordinates,
        )


def _source_failure_reason(error: BaseException) -> str:
    text = str(error)
    if "alpha, lambda, and R positive" in text:
        return "nonpositive_metric_factor_before_source_evaluation"
    if "kinetic condition limit" in text or "kinetic block is singular" in text:
        return "kinetic_condition_limit"
    if "residual tolerance" in text:
        return "newton_residual_limit"
    if "normal-flow tracer" in text:
        return "stopped_measurement_tracer_contract"
    return "invalid_implementation_or_nonconverged_run"


def _member_inputs(
    authorization: Mapping[str, Any],
    amplitude: str,
) -> dict[tuple[str, int], Mapping[str, Any]]:
    records = authorization["artifact_payload"]["frozen_run_inputs"]
    selected = {
        (item["method"], item["point_count"]): item
        for item in records
        if item["amplitude"] == amplitude
    }
    if len(selected) != 6:
        raise ValueError("HLT3 does not contain six inputs for this amplitude")
    return selected


def _build_members(
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
    amplitude: str,
) -> dict[str, RunMember]:
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
    answer: dict[str, RunMember] = {}
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
                raise ValueError("runtime projected state differs from HLT3 input")
            projected_initial = replace(initial, state=state)
            operator = GR0EvolutionOperator(
                initial.grid,
                spatial_order=order,
                ko_dissipation=float(Q(numerics["ko_dissipation"])),
                residual_tolerance=float(
                    Q(thresholds["source_residual_infinity_max"])
                ),
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
            tracers = NormalFlowTracers.create(
                minimum=float(Q(numerics["normal_flow_tracer_radius_minimum"])),
                maximum=float(Q(numerics["normal_flow_tracer_radius_maximum"])),
                spacing=float(Q(numerics["normal_flow_tracer_spacing"])),
                state=state,
                coordinates=initial.grid.coordinates,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                outer_radius=float(Q(physical["outer_radius"])),
            )
            member = RunMember(
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
                raise RuntimeError("duplicate GR-0 run member")
            answer[member.key] = member
    if len(answer) != 6:
        raise RuntimeError("GR-0 amplitude requires six run members")
    return answer


def _method_members(
    members: Mapping[str, RunMember], method: str
) -> tuple[RunMember, ...]:
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
    members: Mapping[str, RunMember],
    coordinate_time: float,
) -> dict[str, Any]:
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    primary = _method_members(members, "RK4")
    comparator = _method_members(members, "SSPRK3")
    primary_common = proto5_gr0_common_event(
        [item.state for item in primary],
        [item.initial.grid for item in primary],
        method="RK4",
        coordinate_time=coordinate_time,
        cutoff=float(Q(physical["cutoff_Lambda"])),
        measurement_radius_maximum=float(
            Q(physical["measurement_radius_maximum"])
        ),
        taper_fraction=float(Q(numerics["proper_spectral_taper_fraction"])),
    )
    comparator_common = proto5_gr0_common_event(
        [item.state for item in comparator],
        [item.initial.grid for item in comparator],
        method="SSPRK3",
        coordinate_time=coordinate_time,
        cutoff=float(Q(physical["cutoff_Lambda"])),
        measurement_radius_maximum=float(
            Q(physical["measurement_radius_maximum"])
        ),
        taper_fraction=float(Q(numerics["proper_spectral_taper_fraction"])),
    )
    trapped = proto5_calibration_trapped_assessment(
        primary_states=[item.state for item in primary],
        comparator_states=[item.state for item in comparator],
        grids=[item.initial.grid for item in primary],
        primary_common_event=primary_common,
        comparator_common_event=comparator_common,
        measurement_radius_maximum=float(
            Q(physical["measurement_radius_maximum"])
        ),
        minimum_observed_order=float(
            Q(numerics["minimum_constraint_finest_pair_order"])
        ),
        positive_margin_factor=float(
            plan["candidate_selection"][
                "trapped_sign_margin_over_combined_error_factor"
            ]
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
            assessed = proto5_temporal_spectral_admission(
                point_counts=[item.point_count for item in records],
                proper_times_by_resolution=[
                    np.asarray(item.tracers.event_proper_times) for item in records
                ],
                field_histories_by_resolution=[
                    np.asarray(item.tracers.event_fields) for item in records
                ],
                field_names=PROTO4_SPECTRAL_FIELD_ORDER,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                taper_fraction=float(
                    Q(numerics["proper_spectral_taper_fraction"])
                ),
            )
            method_temporal[method] = assessed
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


def _checkpoint_arrays(
    *,
    metadata: Mapping[str, Any],
    members: Mapping[str, RunMember],
    event_log_payload: bytes,
) -> dict[str, np.ndarray]:
    arrays: dict[str, np.ndarray] = {
        "metadata_utf8": np.frombuffer(_canonical(metadata), dtype=np.uint8).copy(),
        "event_log_utf8": np.frombuffer(event_log_payload, dtype=np.uint8).copy(),
    }
    for key, member in members.items():
        prefix = key.replace("-", "_")
        arrays[f"{prefix}_u"] = member.state.u
        arrays[f"{prefix}_p"] = member.state.p
        arrays[f"{prefix}_q"] = member.state.q
        arrays[f"{prefix}_tracer_positions"] = member.tracers.positions
        arrays[f"{prefix}_tracer_proper_times"] = member.tracers.proper_times
        arrays[f"{prefix}_event_proper_times"] = np.asarray(
            member.tracers.event_proper_times
        )
        arrays[f"{prefix}_event_fields"] = np.asarray(member.tracers.event_fields)
    return arrays


def _write_checkpoint(
    path: Path,
    *,
    metadata: Mapping[str, Any],
    members: Mapping[str, RunMember],
    event_log_path: Path,
) -> None:
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}.npz")
    if temporary.exists():
        raise FileExistsError(temporary)
    event_log_payload = _event_log_bytes(event_log_path)
    if metadata.get("event_log") != _event_log_identity(event_log_payload):
        raise ValueError("checkpoint metadata does not match its event log")
    arrays = _checkpoint_arrays(
        metadata=metadata,
        members=members,
        event_log_payload=event_log_payload,
    )
    with temporary.open("xb") as handle:
        np.savez_compressed(handle, **arrays)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def _member_metadata(member: RunMember) -> dict[str, Any]:
    return {
        "time": member.time,
        "step_index": member.step_index,
        "transaction_serial": member.transaction_serial,
        "CFL_retry_count": member.CFL_retry_count,
        "runtime_monitor_state": asdict(member.transaction.state),
        "causal_state": asdict(member.transaction.causal_state),
        "input_hash": member.input_hash,
    }


def _restore_checkpoint(
    path: Path,
    *,
    plan: Mapping[str, Any],
    authorization: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, RunMember], bytes]:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(bytes(archive["metadata_utf8"]).decode("utf-8"))
        event_log_payload = bytes(archive["event_log_utf8"])
        if metadata.get("event_log") != _event_log_identity(event_log_payload):
            raise ValueError("checkpoint event-log evidence is internally inconsistent")
        amplitude = metadata["amplitude"]
        members = _build_members(plan, authorization, amplitude)
        for key, member in members.items():
            prefix = key.replace("-", "_")
            stored = metadata["members"][key]
            if stored["input_hash"] != member.input_hash:
                raise ValueError("checkpoint input hash differs from HLT3")
            member.state = EvolutionState(
                archive[f"{prefix}_u"],
                archive[f"{prefix}_p"],
                archive[f"{prefix}_q"],
            )
            member.time = float(stored["time"])
            member.step_index = int(stored["step_index"])
            member.transaction_serial = int(stored["transaction_serial"])
            member.CFL_retry_count = int(stored["CFL_retry_count"])
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
                row.copy()
                for row in archive[f"{prefix}_event_proper_times"]
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
        "git_commit": _git_head(),
        "plan_path": plan_path.relative_to(REPOSITORY).as_posix(),
        "plan_sha256": _sha(plan_path),
        "authorization_path": authorization_path.relative_to(
            REPOSITORY
        ).as_posix(),
        "authorization_sha256": _sha(authorization_path),
        "input_manifest_sha256": frozen["input_manifest_sha256"],
        "ordered_amplitudes": frozen["ordered_amplitudes"],
        "implementation_sha256": authorization["implementation_sha256"],
        "created_unix_time": wall_time.time(),
        "outcome_fields_present_at_creation": False,
    }


def _load_manifest(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict) or path.read_bytes() != _canonical(value):
        raise ValueError("campaign manifest is not canonical")
    return value


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
        or manifest.get("plan_sha256") != _sha(plan_path)
        or manifest.get("authorization_sha256") != _sha(authorization_path)
        or manifest.get("input_manifest_sha256")
        != expected["input_manifest_sha256"]
    ):
        raise ValueError("resume manifest differs from the authorized campaign")


def _campaign_checkpoint_metadata(
    *,
    manifest: Mapping[str, Any],
    amplitude: str,
    amplitude_index: int,
    event_index: int,
    consecutive_qualified_events: int,
    amplitude_records: Sequence[Mapping[str, Any]],
    event_log_path: Path,
    members: Mapping[str, RunMember],
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
        "event_log": _event_log_identity(_event_log_bytes(event_log_path)),
        "terminal": terminal,
        "members": {key: _member_metadata(item) for key, item in members.items()},
    }


def run_campaign(
    *,
    plan_path: Path,
    authorization_path: Path,
    resume: bool,
) -> dict[str, Any]:
    _require_clean_tracked_worktree()
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
        manifest = _load_manifest(manifest_path)
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
        resume_recovery = _recover_event_log_from_checkpoint(
            event_log,
            checkpoint_event_log,
        )
    else:
        if output_root.exists():
            raise ValueError("fresh calibration output root already exists")
        output_root.mkdir(parents=True, exist_ok=False)
        manifest = _manifest(plan_path, authorization_path, authorization)
        _write_exclusive(manifest_path, _canonical(manifest))
        amplitude_index = 0
        event_index = 0
        consecutive = 0
        amplitude_records: list[dict[str, Any]] = []
        amplitude = plan["candidate_selection"]["ordered_amplitudes"][0]
        members = _build_members(plan, authorization, amplitude)
        initial = _event_assessment(plan, members, 0.0)
        _append_fsync(
            event_log,
            _canonical_line(
                {
                    "event_type": "common_event",
                    "amplitude": amplitude,
                    "event_index": 0,
                    "assessment": initial,
                }
            ),
        )
        checkpoint = _campaign_checkpoint_metadata(
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
    maximum_retries = plan["numerics"]["maximum_CFL_retries_per_step"]
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
            try:
                for member in members.values():
                    member.advance_to(
                        target_time,
                        retry_factor=retry_factor,
                        maximum_retries=maximum_retries,
                        minimum_step_size=minimum_step,
                    )
                for member in members.values():
                    member.append_common_event()
                assessment = _event_assessment(plan, members, target_time)
            except GR0RuntimeStop as stop:
                candidate_stop = {
                    "classification": "scientific_universal_runtime_stop",
                    "reason": stop.reason,
                    "failures": list(stop.failures),
                    "target_common_event_index": target_index,
                    "last_fully_completed_common_event_index": event_index,
                }
                break
            except BaseException as error:
                if isinstance(error, KeyboardInterrupt):
                    # Members advance independently between common events.  The
                    # preceding atomic checkpoint is therefore the only valid
                    # cross-member state if an interrupt lands mid-event.
                    raise
                reason = _source_failure_reason(error)
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
                }
                break

            event_index = target_index
            if assessment["qualified_trapped_common_event"]:
                consecutive += 1
            else:
                consecutive = 0
            event_record = {
                "event_type": "common_event",
                "amplitude": amplitude,
                "event_index": event_index,
                "assessment": assessment,
                "consecutive_qualified_trapped_common_events": consecutive,
            }
            _append_fsync(event_log, _canonical_line(event_record))
            if event_index % checkpoint_every_events == 0:
                checkpoint = _campaign_checkpoint_metadata(
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
                "member_final_state_hashes": {
                    key: array_content_sha256(item.state.u, item.state.p, item.state.q)
                    for key, item in members.items()
                },
            }
        )
        if selected_amplitude is not None:
            break
        if candidate_stop is not None and candidate_stop["reason"] == "invalid_implementation_or_nonconverged_run":
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
        _append_fsync(
            event_log,
            _canonical_line(
                {
                    "event_type": "common_event",
                    "amplitude": amplitude,
                    "event_index": 0,
                    "assessment": initial,
                }
            ),
        )
        checkpoint = _campaign_checkpoint_metadata(
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
    event_log_sha256 = _sha(event_log)
    result = {
        "schema_version": 1,
        "runner_id": RUNNER_ID,
        "campaign_id": manifest["campaign_id"],
        "classification": terminal_classification,
        "selected_amplitude": selected_amplitude,
        "amplitude_records": amplitude_records,
        "event_log": event_log.relative_to(REPOSITORY).as_posix(),
        "event_log_sha256": event_log_sha256,
        "elapsed_wall_seconds": wall_time.monotonic() - campaign_started,
        "resume_recovery": resume_recovery,
        "FGCQR_outcome_read": False,
        "SGBL_outcome_read": False,
        "holdout_execution_authorized": False,
        "retained_EFT_evolution_authorized": False,
        "mechanism_question_answered": False,
    }
    _write_exclusive(result_path, _canonical(result))
    checkpoint = _campaign_checkpoint_metadata(
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
    parser.add_argument(
        "--check-authorization-only",
        action="store_true",
        help="validate HLT3 and the clean fresh namespace without creating output",
    )
    args = parser.parse_args()
    plan_path = args.plan.resolve()
    authorization_path = args.authorization.resolve()
    if args.check_authorization_only:
        _require_clean_tracked_worktree()
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
    print(json.dumps(_serial(result), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
