"""Outcome-blind numerical controls for the NUM1 validation gate.

The controls are deliberately standard problems with independently known
answers.  They validate numerical mechanics without reading an FGC-QR
trajectory:

* exact Minkowski-state preservation for all six ADM fields;
* an even, centre-regular manufactured wave with a declared forcing;
* the odd ``psi=r*chi`` GR-0 spherical scalar control;
* reduction-constraint convergence;
* matched-spacing outer-boundary isolation; and
* bitwise checkpoint/restart equivalence.

The production source remains the unredefined REF1 adapter.  These controls do
not substitute their wave equations for ACT1/VAR1.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from math import exp, log, pi
from typing import Callable, Mapping, Sequence

import numpy as np

from .health_monitor import FIELD_ORDER
from .health_monitor import (
    ClassicalHealthMonitor,
    ClassicalHealthStop,
    HealthThresholds,
    StageHealthSnapshot,
)
from .boundary_domain import BoundaryGeometry, CausalBudgetState
from .numerical_engine import (
    COMPARATOR_METHOD,
    METHODS,
    PRIMARY_METHOD,
    EvolutionCheckpoint,
    EvolutionRHS,
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
    array_content_sha256,
    evolve_fixed_steps,
    parity_project,
    propose_step,
)
from .runtime_transaction import RuntimeStageTransaction


@dataclass(frozen=True, slots=True)
class WaveControl:
    control_id: str
    parity: int
    wavenumber: float
    manufactured: bool

    def __post_init__(self) -> None:
        if not self.control_id:
            raise ValueError("wave control requires an identifier")
        if self.parity not in {-1, 1}:
            raise ValueError("wave parity must be +/-1")
        if self.wavenumber <= 0.0:
            raise ValueError("wave wavenumber must be positive")

    def exact(self, time: float, radii: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        k = self.wavenumber
        if self.manufactured:
            amplitude = np.exp(-time)
            u = amplitude * np.cos(k * radii)
            p = -u
            q = -k * amplitude * np.sin(k * radii)
        else:
            u = np.sin(k * radii) * np.cos(k * time)
            p = -k * np.sin(k * radii) * np.sin(k * time)
            q = k * np.cos(k * radii) * np.cos(k * time)
        return u[:, None], p[:, None], q[:, None]

    def forcing(self, time: float, radii: np.ndarray) -> np.ndarray:
        if not self.manufactured:
            return np.zeros((radii.size, 1), dtype=np.float64)
        u, _p, _q = self.exact(time, radii)
        return (1.0 + self.wavenumber**2) * u


@dataclass(frozen=True, slots=True)
class WaveRunResult:
    control_id: str
    method: str
    point_count: int
    spacing: float
    step_size: float
    step_count: int
    final_time: float
    solution_l2_error: float
    solution_infinity_error: float
    reduction_constraint_l2: float
    reduction_constraint_infinity: float
    relative_energy_drift: float
    endpoint_sha256: str
    endpoint: EvolutionState


def _method_order(method: str) -> int:
    if method == PRIMARY_METHOD:
        return 4
    if method == COMPARATOR_METHOD:
        return 2
    raise ValueError("unknown numerical method")


def _wave_projector(
    control: WaveControl,
    grid: UniformRadialGrid,
) -> Callable[[float, EvolutionState], EvolutionState]:
    radii = grid.coordinates
    u_parity = (control.parity,)
    q_parity = (-control.parity,)

    def project(time: float, state: EvolutionState) -> EvolutionState:
        u = parity_project(state.u, u_parity)
        p = parity_project(state.p, u_parity)
        q = parity_project(state.q, q_parity)
        exact_u, exact_p, exact_q = control.exact(time, radii)
        # The outer control is analytic boundary data.  The inner point uses
        # parity rather than an independently imposed Dirichlet value.
        u[-1] = exact_u[-1]
        p[-1] = exact_p[-1]
        q[-1] = exact_q[-1]
        return EvolutionState(u, p, q)

    return project


def _wave_rhs(
    control: WaveControl,
    grid: UniformRadialGrid,
    operator: SBPFirstDerivative,
) -> Callable[[float, EvolutionState], EvolutionRHS]:
    radii = grid.coordinates
    u_parity = (control.parity,)
    q_parity = (-control.parity,)

    def rhs(time: float, state: EvolutionState) -> EvolutionRHS:
        dp = operator.differentiate(state.q, center_parities=q_parity) + control.forcing(time, radii)
        dq = operator.differentiate(state.p, center_parities=u_parity)
        reduction = state.q - operator.differentiate(state.u, center_parities=u_parity)
        return EvolutionRHS(
            du=state.p,
            dp=dp,
            dq=dq,
            diagnostics={
                "reduction_constraint_infinity": float(np.max(np.abs(reduction))),
            },
        )

    return rhs


def run_wave_control(
    control: WaveControl,
    *,
    method: str,
    point_count: int,
    final_time: float,
    cfl: float = 1.0 / 8.0,
) -> WaveRunResult:
    """Run one manufactured or GR-0 scalar control to a fixed time."""

    if method not in METHODS:
        raise ValueError("unknown numerical method")
    grid = UniformRadialGrid(0.0, pi, point_count)
    operator = SBPFirstDerivative(grid, _method_order(method))
    radii = grid.coordinates
    exact_u, exact_p, exact_q = control.exact(0.0, radii)
    initial = EvolutionState(exact_u, exact_p, exact_q)
    nominal = cfl * grid.spacing
    step_count = max(1, int(np.ceil(final_time / nominal)))
    step_size = final_time / step_count
    retained = evolve_fixed_steps(
        method=method,
        initial_state=initial,
        initial_time=0.0,
        step_size=step_size,
        step_count=step_count,
        rhs=_wave_rhs(control, grid, operator),
        projector=_wave_projector(control, grid),
        retain_every=step_count,
    )
    endpoint = retained[-1].state
    target_u, target_p, target_q = control.exact(final_time, radii)
    error = endpoint.u - target_u
    weights = operator.norm_weights[:, None]
    l2 = float(np.sqrt(np.sum(weights * error**2)))
    infinity = float(np.max(np.abs(error)))
    reduction = endpoint.q - operator.differentiate(
        endpoint.u, center_parities=(control.parity,)
    )
    reduction_l2 = float(np.sqrt(np.sum(weights * reduction**2)))
    reduction_infinity = float(np.max(np.abs(reduction)))

    def energy(p: np.ndarray, q: np.ndarray) -> float:
        return 0.5 * float(np.sum(weights * (p**2 + q**2)))

    initial_energy = energy(exact_p, exact_q)
    final_energy = energy(endpoint.p, endpoint.q)
    drift = abs(final_energy - initial_energy) / max(initial_energy, 1.0e-300)
    return WaveRunResult(
        control_id=control.control_id,
        method=method,
        point_count=point_count,
        spacing=grid.spacing,
        step_size=step_size,
        step_count=step_count,
        final_time=final_time,
        solution_l2_error=l2,
        solution_infinity_error=infinity,
        reduction_constraint_l2=reduction_l2,
        reduction_constraint_infinity=reduction_infinity,
        relative_energy_drift=drift,
        endpoint_sha256=array_content_sha256(endpoint.u, endpoint.p, endpoint.q),
        endpoint=endpoint,
    )


def observed_orders(errors: Sequence[float]) -> tuple[float, ...]:
    if len(errors) < 2 or any(value <= 0.0 for value in errors):
        raise ValueError("observed order requires at least two positive errors")
    return tuple(log(coarse / fine, 2.0) for coarse, fine in zip(errors, errors[1:]))


def convergence_control(
    control: WaveControl,
    *,
    method: str,
    point_counts: Sequence[int],
    final_time: float,
) -> Mapping[str, object]:
    runs = tuple(
        run_wave_control(
            control,
            method=method,
            point_count=count,
            final_time=final_time,
        )
        for count in point_counts
    )
    solution_orders = observed_orders([item.solution_l2_error for item in runs])
    constraint_orders = observed_orders(
        [item.reduction_constraint_l2 for item in runs]
    )
    return {
        "control_id": control.control_id,
        "method": method,
        "runs": runs,
        "solution_orders": solution_orders,
        "constraint_orders": constraint_orders,
        "minimum_solution_order": min(solution_orders),
        "minimum_constraint_order": min(constraint_orders),
    }


def minkowski_preservation_control(
    *,
    method: str,
    point_count: int = 65,
    step_count: int = 16,
) -> Mapping[str, object]:
    grid = UniformRadialGrid(0.0, 4.0, point_count)
    radii = grid.coordinates
    field_count = len(FIELD_ORDER)
    u = np.zeros((point_count, field_count), dtype=np.float64)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = radii
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    initial = EvolutionState(u, p, q)
    parities = (1, -1, 1, -1, 1, 1)
    q_parities = tuple(-value for value in parities)

    def projector(_time: float, state: EvolutionState) -> EvolutionState:
        return EvolutionState(
            parity_project(state.u, parities),
            parity_project(state.p, parities),
            parity_project(state.q, q_parities),
        )

    def rhs(_time: float, state: EvolutionState) -> EvolutionRHS:
        zero = np.zeros_like(state.u)
        return EvolutionRHS(state.p, zero, zero, {"control": "MINKOWSKI-EXACT"})

    final = evolve_fixed_steps(
        method=method,
        initial_state=initial,
        initial_time=0.0,
        step_size=1.0 / 128.0,
        step_count=step_count,
        rhs=rhs,
        projector=projector,
        retain_every=step_count,
    )[-1].state
    drift = max(
        float(np.max(np.abs(final.u - initial.u))),
        float(np.max(np.abs(final.p - initial.p))),
        float(np.max(np.abs(final.q - initial.q))),
    )
    return {
        "method": method,
        "maximum_binary64_drift": drift,
        "bitwise_preserved": array_content_sha256(final.u, final.p, final.q)
        == array_content_sha256(initial.u, initial.p, initial.q),
    }


def checkpoint_restart_control(
    *,
    method: str,
    point_count: int = 65,
    total_steps: int = 16,
) -> Mapping[str, object]:
    if total_steps % 2:
        raise ValueError("restart control requires an even step count")
    control = WaveControl("GR0-RESTART", -1, 1.0, False)
    grid = UniformRadialGrid(0.0, pi, point_count)
    operator = SBPFirstDerivative(grid, _method_order(method))
    initial = EvolutionState(*control.exact(0.0, grid.coordinates))
    step_size = (pi / 32.0) / total_steps
    rhs = _wave_rhs(control, grid, operator)
    projector = _wave_projector(control, grid)
    continuous = evolve_fixed_steps(
        method=method,
        initial_state=initial,
        initial_time=0.0,
        step_size=step_size,
        step_count=total_steps,
        rhs=rhs,
        projector=projector,
        retain_every=total_steps,
    )[-1]
    first = evolve_fixed_steps(
        method=method,
        initial_state=initial,
        initial_time=0.0,
        step_size=step_size,
        step_count=total_steps // 2,
        rhs=rhs,
        projector=projector,
        retain_every=total_steps // 2,
    )[-1]
    checkpoint = EvolutionCheckpoint(
        time=first.time,
        step_index=first.step_index,
        transaction_serial=first.transaction_serial,
        state=first.state,
        metadata={"control_id": control.control_id, "method": method},
    )
    restored = EvolutionCheckpoint.from_bytes(checkpoint.canonical_bytes())
    second = evolve_fixed_steps(
        method=method,
        initial_state=restored.state,
        initial_time=restored.time,
        step_size=step_size,
        step_count=total_steps // 2,
        rhs=rhs,
        projector=projector,
        retain_every=total_steps // 2,
    )[-1]
    equal = all(
        np.array_equal(getattr(continuous.state, name), getattr(second.state, name))
        for name in ("u", "p", "q")
    )
    return {
        "method": method,
        "checkpoint_sha256": checkpoint.sha256,
        "canonical_roundtrip_bitwise": restored.canonical_bytes() == checkpoint.canonical_bytes(),
        "continuous_and_restart_endpoint_bitwise_equal": equal,
        "endpoint_sha256": array_content_sha256(
            second.state.u, second.state.p, second.state.q
        ),
    }


def _compact_bump(radius: np.ndarray, *, center: float, half_width: float) -> np.ndarray:
    x = (radius - center) / half_width
    answer = np.zeros_like(radius)
    mask = np.abs(x) < 1.0
    answer[mask] = np.exp(1.0 - 1.0 / (1.0 - x[mask] ** 2))
    return answer


def boundary_isolation_control(*, method: str) -> Mapping[str, object]:
    """Compare aligned domains before their continuum boundary can arrive."""

    spacing = 1.0 / 32.0
    domains = (6.0, 8.0)
    endpoints = []
    for outer in domains:
        count = int(round(outer / spacing)) + 1
        grid = UniformRadialGrid(0.0, outer, count)
        operator = SBPFirstDerivative(grid, _method_order(method))
        radii = grid.coordinates
        u = _compact_bump(radii, center=1.5, half_width=0.5)[:, None]
        p = np.zeros_like(u)
        q = operator.differentiate(u, center_parities=(1,))
        initial = EvolutionState(u, p, q)

        def projector(_time: float, state: EvolutionState) -> EvolutionState:
            projected_u = parity_project(state.u, (1,))
            projected_p = parity_project(state.p, (1,))
            projected_q = parity_project(state.q, (-1,))
            projected_u[-1] = 0.0
            projected_p[-1] = 0.0
            projected_q[-1] = 0.0
            return EvolutionState(projected_u, projected_p, projected_q)

        def rhs(_time: float, state: EvolutionState) -> EvolutionRHS:
            return EvolutionRHS(
                state.p,
                operator.differentiate(state.q, center_parities=(-1,)),
                operator.differentiate(state.p, center_parities=(1,)),
                {},
            )

        step_size = spacing / 8.0
        step_count = 32
        endpoint = evolve_fixed_steps(
            method=method,
            initial_state=initial,
            initial_time=0.0,
            step_size=step_size,
            step_count=step_count,
            rhs=rhs,
            projector=projector,
            retain_every=step_count,
        )[-1].state
        endpoints.append((grid, endpoint))
    measurement_points = int(round(2.5 / spacing)) + 1
    differences = tuple(
        float(
            np.max(
                np.abs(
                    getattr(endpoints[0][1], name)[:measurement_points]
                    - getattr(endpoints[1][1], name)[:measurement_points]
                )
            )
        )
        for name in ("u", "p", "q")
    )
    return {
        "method": method,
        "spacing_preserved": endpoints[0][0].spacing == endpoints[1][0].spacing,
        "measurement_radius": 2.5,
        "outer_radii": domains,
        "maximum_measurement_difference": max(differences),
        "component_differences": dict(zip(("u", "p", "q"), differences, strict=True)),
    }


def sbp_and_dissipation_controls() -> Mapping[str, object]:
    records = []
    for order in (4, 2):
        grid = UniformRadialGrid(0.0, 1.0, 65)
        operator = SBPFirstDerivative(grid, order)
        high = ((-1.0) ** np.arange(grid.point_count))[:, None]
        # Zero the outer filter closure and test only the actual filtered zone.
        dissipation = operator.kreiss_oliger_dissipation(
            high, coefficient=1.0 / 64.0
        )
        energy_rate = float(
            np.sum(operator.norm_weights[:, None] * high * dissipation)
        )
        records.append(
            {
                "order": order,
                "sbp_identity_residual": operator.sbp_identity_residual(),
                "high_frequency_energy_rate": energy_rate,
                "dissipation_nonpositive": energy_rate < 0.0,
                "stencil_reach_intervals": operator.stencil_reach_intervals,
            }
        )
    return {"records": tuple(records)}


def _healthy_runtime_snapshot(
    thresholds: HealthThresholds,
    *,
    time: float,
    stage_index: int,
    boundary_margin: float,
) -> StageHealthSnapshot:
    """Construct a strictly interior synthetic monitor point for NUM1 only."""

    return StageHealthSnapshot(
        time=time,
        stage_index=stage_index,
        minimum_lapse=1.0,
        minimum_radial_metric=1.0,
        minimum_areal_radius_away_from_center=0.01,
        minimum_hat_lorentzian_margin=0.5,
        newton_residual_infinity=0.0,
        newton_iterations=min(2, thresholds.newton_iteration_max),
        newton_residual_monotonic=True,
        kinetic_condition_infinity=min(100.0, thresholds.kinetic_condition_max / 2.0),
        normalized_branch_displacement=thresholds.branch_displacement_max / 4.0,
        acceleration_infinity=thresholds.acceleration_max / 4.0,
        companion_deformation_operator_2=thresholds.companion_deformation_max / 4.0,
        action_energy_deformation_operator_2=thresholds.action_energy_deformation_max / 4.0,
        kinetic_deformation_operator_2=thresholds.kinetic_deformation_max / 4.0,
        effective_planck_ratio=max(1.0, 2.0 * thresholds.effective_planck_ratio_min),
        characteristic_imaginary_part=0.0,
        eigenframe_condition_number=min(100.0, thresholds.eigenframe_condition_max / 2.0),
        symmetrizer_coercivity=max(0.5, 2.0 * thresholds.symmetrizer_coercivity_min),
        normalized_constraint_infinity=0.0,
        observed_constraint_convergence_order=thresholds.constraint_convergence_order_min + 0.5,
        phi_over_cutoff=thresholds.phi_over_cutoff_max / 4.0,
        proper_frequency_over_cutoff=thresholds.proper_frequency_over_cutoff_max / 4.0,
        proper_wavenumber_over_cutoff=thresholds.proper_wavenumber_over_cutoff_max / 4.0,
        curvature_scale_squared_over_cutoff_squared=(
            thresholds.curvature_scale_squared_over_cutoff_squared_max / 4.0
        ),
        temporal_alias_band_occupied=False,
        radial_alias_band_occupied=False,
        boundary_margin_over_required_buffer=boundary_margin,
    )


def runtime_transaction_control(
    thresholds: HealthThresholds,
) -> Mapping[str, object]:
    """Exercise atomic HLT1/BND2 composition for both frozen integrators.

    This is a synthetic transaction control.  Its zero right-hand side is not
    an FGC solution and none of its monitor values are trajectory evidence.
    """

    if not isinstance(thresholds, HealthThresholds):
        raise TypeError("runtime transaction control requires HealthThresholds")
    initial = EvolutionState(
        np.zeros((9, len(FIELD_ORDER))),
        np.zeros((9, len(FIELD_ORDER))),
        np.zeros((9, len(FIELD_ORDER))),
    )

    def rhs(_time: float, state: EvolutionState) -> EvolutionRHS:
        zero = np.zeros(state.shape)
        return EvolutionRHS(zero, zero, zero, {"control": "SYNTHETIC-TRANSACTION"})

    records = []
    for method in METHODS:
        proposal = propose_step(
            method=method,
            time=0.0,
            step_size=0.125,
            state=initial,
            rhs=rhs,
        )

        def healthy_factory(stage, serial, margin):
            return _healthy_runtime_snapshot(
                thresholds,
                time=stage.time,
                stage_index=serial,
                boundary_margin=margin,
            )

        healthy = RuntimeStageTransaction(
            health_monitor=ClassicalHealthMonitor(thresholds),
            causal_state=CausalBudgetState(previous_speed_upper=1.0),
            boundary_geometry=BoundaryGeometry(96.0, 24.0, 16.0, 0.375),
            speed_evaluator=lambda _stage: 1.25,
            snapshot_factory=healthy_factory,
        )
        receipt = dict(healthy(proposal))

        def late_failure_factory(stage, serial, margin):
            snapshot = healthy_factory(stage, serial, margin)
            if stage.stage_name == "candidate_endpoint":
                return replace(snapshot, minimum_lapse=0.0)
            return snapshot

        failure = RuntimeStageTransaction(
            health_monitor=ClassicalHealthMonitor(thresholds),
            causal_state=CausalBudgetState(previous_speed_upper=1.0),
            boundary_geometry=BoundaryGeometry(96.0, 24.0, 16.0, 0.375),
            speed_evaluator=lambda _stage: 1.25,
            snapshot_factory=late_failure_factory,
        )
        causal_before = failure.causal_state
        try:
            failure(proposal)
        except ClassicalHealthStop as stop:
            late_reason = stop.reason
        else:
            raise ValueError("late synthetic runtime failure was accepted")
        rolled_back = (
            failure.health_monitor.state.accepted_stage_count == 0
            and failure.health_monitor.state.last_accepted_stage_index == -1
            and failure.causal_state == causal_before
        )
        try:
            failure(proposal)
        except ClassicalHealthStop as repeated:
            immutable = repeated.reason == late_reason == "nonpositive_lapse"
        else:
            immutable = False
        if not rolled_back or not immutable:
            raise ValueError("runtime transaction failed atomic rollback control")

        boundary = RuntimeStageTransaction(
            health_monitor=ClassicalHealthMonitor(thresholds),
            causal_state=CausalBudgetState(previous_speed_upper=1.0),
            boundary_geometry=BoundaryGeometry(0.75, 0.2, 0.4, 0.1),
            speed_evaluator=lambda _stage: 1.0,
            snapshot_factory=healthy_factory,
        )
        boundary_before = boundary.causal_state
        try:
            boundary(proposal)
        except ClassicalHealthStop as stop:
            boundary_reason = stop.reason
        else:
            raise ValueError("synthetic boundary-exhaustion proposal was accepted")
        boundary_atomic = (
            boundary_reason == "boundary_causal_buffer"
            and boundary.health_monitor.state.accepted_stage_count == 0
            and boundary.causal_state == boundary_before
        )
        if not boundary_atomic:
            raise ValueError("runtime boundary transaction failed atomic control")
        records.append(
            {
                "method": method,
                "stage_count_including_candidate_endpoint": len(proposal.stages),
                "healthy_receipt": receipt,
                "late_failure_reason": late_reason,
                "late_failure_rolled_back_all_internal_stages": rolled_back,
                "first_failure_immutable": immutable,
                "boundary_failure_reason": boundary_reason,
                "boundary_failure_rolled_back_causal_and_health_state": boundary_atomic,
            }
        )
    return {"records": tuple(records), "all_transactions_passed": True}


__all__ = [
    "WaveControl",
    "WaveRunResult",
    "boundary_isolation_control",
    "checkpoint_restart_control",
    "convergence_control",
    "minkowski_preservation_control",
    "observed_orders",
    "run_wave_control",
    "runtime_transaction_control",
    "sbp_and_dissipation_controls",
]
