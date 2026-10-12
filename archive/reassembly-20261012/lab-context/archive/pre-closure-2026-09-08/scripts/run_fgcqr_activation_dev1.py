#!/usr/bin/env python3
"""Run one result-first FGC-QR activation trajectory on the revised action."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from hashlib import sha256
import json
import multiprocessing
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.action import ActionParameters, ModelID  # noqa: E402
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    radial_null_observables,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionRHS,
    EvolutionState,
    PRIMARY_METHOD,
    SBPFirstDerivative,
    accept_step,
    array_content_sha256,
    propose_step,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    exact_spherical_minkowski_reference,
    reference_balanced_spatial_derivatives,
)
from recursive_horizons.fgc.evolution.static_initial_admission import (  # noqa: E402
    construct_fgcqr_grid_initial_data,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
    BatchScalar,
    regular_center_acceleration_limit,
    solve_grid_accelerations,
)
from recursive_horizons.fgc.initial_data_family import (  # noqa: E402
    FGCQRActionParameters,
    InitialDataParameters,
)
from recursive_horizons.fgc import spherical_reduction as spherical_reduction  # noqa: E402


ARTIFACT_ID = "FGC-1-FGCQR-ACTIVATION-A6042-DEV1"
SCHEMA = "FGC-1-FGCQR-ACTIVATION-A6042-DEV1-v1"
OUTPUT = ROOT / "runs/fgc-2-sf1/fgcqr-activation-a6042-dev1.json"
HISTORY = ROOT / "runs/fgc-2-sf1/fgcqr-activation-a6042-dev1-history.npz"
ACTION = FGCQRActionParameters(
    planck_mass=2.0,
    beta=-128.0,
    scalar_mass=0.02,
    quartic_coupling=256.0,
    eta=16384.0,
)
NUMERICAL_ACTION = ActionParameters(
    model_id=ModelID.FGC_QR,
    planck_mass=ACTION.planck_mass,
    scalar_mass=ACTION.scalar_mass,
    quartic_coupling=ACTION.quartic_coupling,
    pulse_width=2.0,
    scalar_field=0.0,
    ricci_coupling=ACTION.beta,
    quadratic_gb_coupling=ACTION.eta,
)
Q_CENTER_PARITIES = -ADM_CENTER_PARITIES
_ORIGINAL_GB = spherical_reduction._gb


def _fast_batch_gb(riemann, ricci, scalar, inverse):
    """Same four-index contraction as RED1, using ndarray axes for BatchScalar."""

    if not isinstance(scalar, BatchScalar):
        return _ORIGINAL_GB(riemann, ricci, scalar, inverse)
    shape = scalar.data.shape
    gi = np.empty(shape + (4, 4), dtype=np.float64)
    ric = np.empty(shape + (4, 4), dtype=np.float64)
    ri = np.empty(shape + (4, 4, 4, 4), dtype=np.float64)
    for a in range(4):
        for b in range(4):
            gi[..., a, b] = np.broadcast_to(BatchScalar.coerce(inverse[a][b]).data, shape)
            ric[..., a, b] = np.broadcast_to(BatchScalar.coerce(ricci[a][b]).data, shape)
            for c in range(4):
                for d in range(4):
                    ri[..., a, b, c, d] = np.broadcast_to(
                        BatchScalar.coerce(riemann[a][b][c][d]).data, shape
                    )
    raised = np.einsum("...ae,...abcd->...ebcd", gi, ri, optimize=True)
    raised = np.einsum("...bf,...ebcd->...efcd", gi, raised, optimize=True)
    raised = np.einsum("...cg,...efcd->...efgd", gi, raised, optimize=True)
    raised = np.einsum("...dh,...efgd->...efgh", gi, raised, optimize=True)
    r2 = np.einsum("...abcd,...abcd->...", raised, ri, optimize=True)
    ricci2 = np.einsum("...ac,...bd,...ab,...cd->...", gi, gi, ric, ric, optimize=True)
    return BatchScalar(r2 - 4.0 * ricci2 + scalar.data * scalar.data)


spherical_reduction._gb = _fast_batch_gb


def _verify_fast_batch_gb() -> float:
    rng = np.random.default_rng(20260903)
    shape = (1, 1)
    inverse = [
        [BatchScalar(rng.normal(size=shape)) for _ in range(4)] for _ in range(4)
    ]
    ricci = [
        [BatchScalar(rng.normal(size=shape)) for _ in range(4)] for _ in range(4)
    ]
    riemann = [
        [
            [
                [BatchScalar(rng.normal(size=shape)) for _ in range(4)]
                for _ in range(4)
            ]
            for _ in range(4)
        ]
        for _ in range(4)
    ]
    scalar = BatchScalar(rng.normal(size=shape))
    expected = _ORIGINAL_GB(riemann, ricci, scalar, inverse).data
    actual = _fast_batch_gb(riemann, ricci, scalar, inverse).data
    error = float(np.max(np.abs(actual - expected)))
    scale = max(1.0, float(np.max(np.abs(expected))))
    if error > 4096.0 * np.finfo(np.float64).eps * scale:
        raise RuntimeError("fast BatchScalar GB contraction differs from RED1")
    return error


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _projector(initial):
    reference = initial.state

    def project(_time: float, state: EvolutionState) -> EvolutionState:
        u = state.u.copy()
        p = state.p.copy()
        q = state.q.copy()
        u[0, ADM_CENTER_PARITIES == -1] = 0.0
        p[0, ADM_CENTER_PARITIES == -1] = 0.0
        q[0, Q_CENTER_PARITIES == -1] = 0.0
        q[0, 3] = u[0, 2]
        u[-4:] = reference.u[-4:]
        p[-4:] = reference.p[-4:]
        q[-4:] = reference.q[-4:]
        return EvolutionState(u, p, q)

    return project


def _solve_source_chunk(payload):
    u, p, q, p_r, q_r, radii = payload
    return solve_grid_accelerations(
        u,
        p,
        q,
        p_r,
        q_r,
        radii,
        action=NUMERICAL_ACTION,
        residual_tolerance=1.0e-10,
        condition_number_maximum=1.0e10,
        branch_displacement_maximum=64.0,
    )


class FGCQROperator:
    def __init__(
        self,
        grid,
        *,
        spatial_order: int,
        ko_dissipation: float,
        pool: ProcessPoolExecutor,
        workers: int,
    ) -> None:
        self.grid = grid
        self.derivative = SBPFirstDerivative(grid, spatial_order)
        self.ko_dissipation = ko_dissipation
        self.pool = pool
        self.workers = workers
        self.last_physical_acceleration = None
        self.last_p_r = None
        self.last_q_r = None
        self.last_diagnostics = None

    def __call__(self, _time: float, state: EvolutionState) -> EvolutionRHS:
        p_r, q_r, differentiated_u = reference_balanced_spatial_derivatives(
            state, self.derivative
        )
        chunks = [
            chunk
            for chunk in np.array_split(np.arange(1, self.grid.point_count), self.workers)
            if chunk.size
        ]
        payloads = [
            (
                state.u[chunk],
                state.p[chunk],
                state.q[chunk],
                p_r[chunk],
                q_r[chunk],
                self.grid.coordinates[chunk],
            )
            for chunk in chunks
        ]
        sources = list(self.pool.map(_solve_source_chunk, payloads))
        acceleration = np.empty_like(state.p)
        acceleration[1:] = np.concatenate(
            [source.accelerations for source in sources], axis=0
        )
        center = regular_center_acceleration_limit(acceleration[1:5])
        acceleration[0] = center.acceleration
        dissipation = self.derivative.kreiss_oliger_dissipation(
            state.p,
            coefficient=self.ko_dissipation,
            center_parities=ADM_CENTER_PARITIES,
        )
        effective_planck = ACTION.planck_mass**2 + ACTION.beta * state.u[:, 4] ** 2
        speed = np.abs(state.u[:, 1]) + 1.2 * state.u[:, 0] / state.u[:, 2]
        diagnostics = {
            "source_residual_infinity": max(source.residual_infinity for source in sources),
            "source_initial_residual_infinity": max(
                source.initial_residual_infinity for source in sources
            ),
            "kinetic_condition_infinity": max(
                source.condition_infinity_maximum for source in sources
            ),
            "branch_displacement_infinity": max(
                source.branch_displacement_infinity for source in sources
            ),
            "source_iterations": max(source.iterations for source in sources),
            "source_affine_verified": all(
                source.affine_in_accelerations_verified for source in sources
            ),
            "center_acceleration_estimator_infinity": center.estimator_infinity,
            "reduction_constraint_infinity": float(
                np.max(np.abs(state.q - differentiated_u), initial=0.0)
            ),
            "coordinate_speed_upper": float(np.max(speed)),
            "minimum_effective_planck_coefficient": float(np.min(effective_planck)),
        }
        self.last_physical_acceleration = acceleration.copy()
        self.last_p_r = p_r.copy()
        self.last_q_r = q_r.copy()
        self.last_diagnostics = diagnostics
        return EvolutionRHS(state.p, acceleration + dissipation, p_r, diagnostics)


def _sample(time: float, state: EvolutionState, operator: FGCQROperator, seed_peak: float):
    operator(time, state)
    observables = radial_null_observables(state)
    radii = operator.grid.coordinates
    measured = (radii > 0.0) & (radii <= 24.0)
    trapped = measured & (observables.theta_plus < 0.0) & (
        observables.theta_minus < 0.0
    )
    phi_peak = float(np.max(np.abs(state.u[measured, 4]), initial=0.0))
    return {
        "time": time,
        "time_hex": time.hex(),
        "activation_factor": phi_peak / seed_peak,
        "phi_peak": phi_peak,
        "phi_peak_radius": float(
            radii[np.argmax(np.where(measured, np.abs(state.u[:, 4]), -np.inf))]
        ),
        "trapped": bool(np.any(trapped)),
        "trapped_point_count": int(np.sum(trapped)),
        "minimum_theta_plus": float(np.min(observables.theta_plus[measured])),
        "minimum_theta_minus": float(np.min(observables.theta_minus[measured])),
        "maximum_compactness": float(np.max(observables.compactness[measured])),
        **operator.last_diagnostics,
    }


def run(
    *,
    point_count: int,
    target_time: float,
    stop_factor: float,
    workers: int,
    output_tag: str,
) -> dict[str, object]:
    output = (
        OUTPUT
        if output_tag == "dev1"
        else ROOT / f"runs/fgc-2-sf1/fgcqr-activation-a6042-{output_tag}.json"
    )
    history_path = (
        HISTORY
        if output_tag == "dev1"
        else ROOT / f"runs/fgc-2-sf1/fgcqr-activation-a6042-{output_tag}-history.npz"
    )
    artifact_id = (
        ARTIFACT_ID
        if output_tag == "dev1"
        else f"FGC-1-FGCQR-ACTIVATION-A6042-{output_tag.upper()}"
    )
    schema = f"{artifact_id}-v1"
    if output.exists() or history_path.exists():
        raise RuntimeError(f"activation {output_tag} one-shot output already exists")
    fast_gb_comparison_error = _verify_fast_batch_gb()
    parameters = InitialDataParameters(
        chi_amplitude=6.042,
        chi_half_width=2.0,
        phi_amplitude=1.0 / 131072.0,
        center=12.0,
        phi_half_width=2.0,
        outer_radius=128.0,
        action=ACTION,
    )
    initial = construct_fgcqr_grid_initial_data(
        parameters,
        point_count=point_count,
        constraint_method="RK4",
        diagnostic_spatial_order=4,
    )
    derivative = SBPFirstDerivative(initial.grid, 4)
    reference = exact_spherical_minkowski_reference(initial.grid)
    q = derivative.differentiate(
        initial.state.u - reference.u,
        center_parities=ADM_CENTER_PARITIES,
    ) + reference.q
    state = EvolutionState(initial.state.u, initial.state.p, q)
    initial = replace(initial, state=state)
    projector = _projector(initial)
    seed_peak = float(np.max(np.abs(state.u[:, 4])))
    time = 0.0
    step_index = 0
    serial = 0
    retries = 0
    samples = []
    history = {name: [] for name in ("times", "u", "p", "q", "a", "p_r", "q_r")}

    def retain() -> dict[str, object]:
        item = _sample(time, state, operator, seed_peak)
        samples.append(item)
        history["times"].append(time)
        history["u"].append(state.u.copy())
        history["p"].append(state.p.copy())
        history["q"].append(state.q.copy())
        history["a"].append(operator.last_physical_acceleration.copy())
        history["p_r"].append(operator.last_p_r.copy())
        history["q_r"].append(operator.last_q_r.copy())
        return item

    terminal_kind = "target_time_without_resolved_activation"
    error = None
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=multiprocessing.get_context("fork"),
    ) as pool:
        operator = FGCQROperator(
            initial.grid,
            spatial_order=4,
            ko_dissipation=1.0 / 64.0,
            pool=pool,
            workers=workers,
        )
        retain()
        while time < target_time - 1.0e-14:
            speed = max(float(operator.last_diagnostics["coordinate_speed_upper"]), 1.0e-12)
            width = min(target_time - time, 0.125 * initial.grid.spacing / speed)
            local_retries = 0
            while True:
                try:
                    proposal = propose_step(
                        method=PRIMARY_METHOD,
                        time=time,
                        step_size=width,
                        state=state,
                        rhs=operator,
                        projector=projector,
                    )
                    candidate = proposal.candidate_state
                    effective_planck = ACTION.planck_mass**2 + ACTION.beta * candidate.u[:, 4] ** 2
                    if (
                        np.any(candidate.u[:, 0] <= 0.0)
                        or np.any(candidate.u[:, 2] <= 0.0)
                        or np.any(effective_planck <= 0.0)
                    ):
                        raise ValueError("candidate crossed a positive metric or F boundary")
                    accepted = accept_step(
                        proposal,
                        previous_step_index=step_index,
                        previous_transaction_serial=serial,
                    )
                    break
                except (ArithmeticError, ValueError) as exc:
                    local_retries += 1
                    retries += 1
                    width *= 0.5
                    if local_retries > 20 or width < 2.0**-30:
                        error = {"type": type(exc).__name__, "message": str(exc)}
                        terminal_kind = "typed_numerical_or_health_stop"
                        break
            if error is not None:
                break
            state = accepted.state
            time = accepted.time
            step_index = accepted.step_index
            serial = accepted.transaction_serial
            item = retain()
            if item["activation_factor"] >= stop_factor:
                terminal_kind = "resolved_activation_target_reached"
                break

    np.savez_compressed(
        history_path,
        **{name: np.asarray(values) for name, values in history.items()},
        radii=initial.grid.coordinates,
    )
    result = {
        "artifact_id": artifact_id,
        "schema": schema,
        "classification": terminal_kind,
        "action": {
            "planck_mass": ACTION.planck_mass,
            "mu": ACTION.scalar_mass,
            "g4": ACTION.quartic_coupling,
            "beta": ACTION.beta,
            "eta": ACTION.eta,
            "alpha": 0.0,
        },
        "initial_data": {
            "amplitude": 6.042,
            "point_count": point_count,
            "constraint_residual_infinity": initial.physical_constraint_residual_infinity,
            "minimum_effective_planck_coefficient": initial.minimum_effective_planck_coefficient,
            "initially_untrapped": initial.no_initial_trapped_sphere,
            "state_sha256": array_content_sha256(state.u, state.p, state.q),
        },
        "execution": {
            "method": "RK4",
            "spatial_order": 4,
            "target_time": target_time,
            "stop_activation_factor": stop_factor,
            "final_time": time,
            "accepted_steps": step_index,
            "retries": retries,
            "source_workers": workers,
            "fast_batch_GB_RED1_comparison_error": fast_gb_comparison_error,
            "error": error,
        },
        "samples": samples,
        "history": {
            "path": str(history_path.relative_to(ROOT)),
            "sha256": _sha(history_path),
            "sample_count": len(samples),
        },
        "resolved_activation": terminal_kind == "resolved_activation_target_reached",
        "terminal": True,
        "nonclaims": {
            "two_method_three_resolution_confirmation": False,
            "metric_null_defocusing_demonstrated": False,
            "robust_neighborhood_demonstrated": False,
        },
    }
    output.write_bytes(canonical_json_bytes(result))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--point-count", type=int, default=513)
    parser.add_argument("--target-time", type=float, default=1.0)
    parser.add_argument("--stop-factor", type=float, default=2.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output-tag", default="dev1")
    args = parser.parse_args()
    result = run(
        point_count=args.point_count,
        target_time=args.target_time,
        stop_factor=args.stop_factor,
        workers=args.workers,
        output_tag=args.output_tag,
    )
    print(
        json.dumps(
            {
                "classification": result["classification"],
                "execution": result["execution"],
                "final_sample": result["samples"][-1],
                "history": result["history"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
