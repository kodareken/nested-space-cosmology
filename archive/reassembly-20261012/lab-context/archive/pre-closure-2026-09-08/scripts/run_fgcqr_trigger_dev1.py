#!/usr/bin/env python3
"""Measure the fixed FGC-QR curvature trigger on the trapped GR development state."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from fractions import Fraction
from hashlib import sha256
import json
import multiprocessing
from pathlib import Path
import sys
import tomllib

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402
from recursive_horizons.fgc.action import (  # noqa: E402
    ActionParameters,
    ModelID,
    linear_effective_mass_squared,
)
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.floating_jet import (  # noqa: E402
    FloatJet2,
    scalar_primal,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
    make_gr0_center_boundary_projector,
    radial_null_observables,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    SBPFirstDerivative,
    accept_step,
    array_content_sha256,
    propose_step,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    CFLRetryRequired,
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    project_gr0_reference_balanced_state,
    reference_balanced_spatial_derivatives,
)
from recursive_horizons.fgc.evolution.proto12_runtime import (  # noqa: E402
    Proto12GR0EvolutionOperator,
    diagnose_gr0_proto12_accelerations_chunked,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    regular_center_acceleration_limit,
)
from recursive_horizons.fgc.initial_data_preflight import PulseParameters  # noqa: E402
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    _gb,
    _ricci,
    _scalar,
    _state_from_adm_pg_jets,
    _warped_2plus2_metric_connection,
    direct_4d_curvature,
    warped_2plus2_curvature,
)


ARTIFACT_ID = "FGC-1-FGCQR-TRIGGER-A6042-DEV1"
SCHEMA = "FGC-1-FGCQR-TRIGGER-A6042-DEV1-v1"
AMPLITUDE = 6.042
TARGET_TIME = 0.02
MEASUREMENT_RADIUS = 24.0
OUTPUT_DEV1 = ROOT / "runs/fgc-2-sf1/fgcqr-trigger-a6042-dev1.json"
OUTPUT_DEV2 = ROOT / "runs/fgc-2-sf1/fgcqr-trigger-a6042-dev2.json"
TRAPPED_ORIGIN = ROOT / "results/fgc-1-gr0-trap-a6042-dev1.json"
EXPECTED_TRAPPED_ORIGIN_SHA256 = (
    "5334e859a1835e1b220fb1861ddc2bb206737193b2aa9f7919eeeff35a2229b0"
)


def _q(value: object) -> float:
    return float(Fraction(value))


def _sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _build_trapped_state(point_count: int):
    config = tomllib.loads((ROOT / "configs/fgc/fgc-1-cal9-run1.toml").read_text())
    physical = config["physical_inputs"]
    numerics = config["numerics"]
    thresholds = config["universal_thresholds"]
    method = config["primary_method"]
    order = int(method["spatial_order"])
    initial = construct_gr0_grid_initial_data(
        PulseParameters(
            chi_amplitude=AMPLITUDE,
            center=_q(physical["chi_center"]),
            half_width=_q(physical["chi_half_width"]),
            phi_amplitude=_q(physical["phi_seed_amplitude"]),
            planck_mass=_q(physical["planck_mass"]),
            scalar_mass=_q(physical["scalar_mass"]),
            quartic_coupling=_q(physical["quartic_coupling"]),
        ),
        point_count=point_count,
        outer_radius=_q(physical["outer_radius"]),
        constraint_method=str(method["constraint_solve_method"]),
        diagnostic_spatial_order=order,
    )
    state = project_gr0_reference_balanced_state(initial, spatial_order=order)
    initial = replace(initial, state=state)
    operator = Proto12GR0EvolutionOperator(
        initial.grid,
        spatial_order=order,
        ko_dissipation=_q(numerics["ko_dissipation"]),
        raw_tolerance=_q(thresholds["source_residual_infinity_max"]),
        kinetic_condition_maximum=_q(thresholds["kinetic_condition_number_max"]),
        maximum_refinement_iterations=int(thresholds["source_iteration_max"]),
        point_batch_size=int(numerics["source_point_batch_size"]),
    )
    initial_rhs = operator(0.0, state)
    derivative = SBPFirstDerivative(initial.grid, order)
    transaction = GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(
            source_residual_maximum=_q(thresholds["source_residual_infinity_max"]),
            source_iteration_maximum=int(thresholds["source_iteration_max"]),
            kinetic_condition_maximum=_q(thresholds["kinetic_condition_number_max"]),
        ),
        causal_state=CausalBudgetState(
            previous_speed_upper=float(initial_rhs.diagnostics["coordinate_speed_upper"])
        ),
        boundary_geometry=BoundaryGeometry(
            _q(physical["outer_radius"]),
            _q(physical["measurement_radius_maximum"]),
            _q(thresholds["minimum_boundary_causal_buffer"]),
            derivative.stencil_reach_intervals * initial.grid.spacing,
        ),
        grid_spacing=initial.grid.spacing,
        cfl_maximum=_q(numerics["cfl_maximum"]),
        hat_normal_factor=_q(thresholds["hat_normal_factor"]),
    )
    projector = make_gr0_center_boundary_projector(
        initial, fixed_outer_rows=int(numerics["fixed_outer_rows"])
    )
    time = 0.0
    step_index = 0
    serial = 0
    cfl_retries = 0
    while time < TARGET_TIME - 1.0e-14:
        speed = max(transaction.causal_state.previous_speed_upper, np.finfo(np.float64).tiny)
        width = min(
            TARGET_TIME - time,
            _q(numerics["cfl_maximum"]) * initial.grid.spacing / speed,
        )
        while True:
            try:
                proposal = propose_step(
                    method=str(method["integrator_id"]),
                    time=time,
                    step_size=width,
                    state=state,
                    rhs=operator,
                    projector=projector,
                )
                accepted = accept_step(
                    proposal,
                    previous_step_index=step_index,
                    previous_transaction_serial=serial,
                    guards=(transaction,),
                )
                break
            except CFLRetryRequired:
                cfl_retries += 1
                width *= _q(numerics["retry_factor"])
        state = accepted.state
        time = accepted.time
        step_index = accepted.step_index
        serial = accepted.transaction_serial
    p_r, q_r, _ = reference_balanced_spatial_derivatives(state, derivative)
    rhs = operator(time, state)
    return initial.grid, state, rhs.dp, p_r, q_r, step_index, cfl_retries


def _curvature_input(index: int, radii, u, p, q, acceleration, p_r, q_r):
    jets = {
        name: FloatJet2(
            u[index, field],
            p[index, field],
            q[index, field],
            acceleration[index, field],
            p_r[index, field],
            q_r[index, field],
        )
        for field, name in enumerate(
            ("alpha", "shift", "lambda", "areal_radius", "phi", "chi")
        )
    }
    fixture = {
        "model_id": "GR-0",
        "action_parameters": {
            "planck_mass": 2,
            "scalar_mass": 3,
            "quartic_coupling": Fraction(1, 2),
            "beta": 0,
            "alpha_gb": 0,
            "eta": 0,
        },
    }
    return index, float(radii[index]), _state_from_adm_pg_jets(fixture, jets)


def _warped_scalars(item):
    index, radius, state = item
    riemann = warped_2plus2_curvature(state)
    _metric, inverse, _connection = _warped_2plus2_metric_connection(state)
    ricci = _ricci(riemann, inverse)
    scalar = _scalar(ricci, inverse)
    gb = _gb(riemann, ricci, scalar, inverse)
    ricci_value = scalar_primal(scalar)
    gb_value = scalar_primal(gb)
    action = ActionParameters(
        model_id=ModelID.FGC_QR,
        planck_mass=2.0,
        scalar_mass=3.0,
        quartic_coupling=0.5,
        pulse_width=2.0,
        scalar_field=0.0,
        ricci_coupling=-0.25,
        linear_gb_coupling=0.0,
        quadratic_gb_coupling=0.5,
    )
    mass_squared = linear_effective_mass_squared(
        action, ricci_scalar=ricci_value, gauss_bonnet=gb_value, phi=0.0
    )
    return index, radius, ricci_value, gb_value, mass_squared


def _direct_scalars(item):
    index, radius, state = item
    riemann, inverse, _connection = direct_4d_curvature(state)
    ricci = _ricci(riemann, inverse)
    scalar = _scalar(ricci, inverse)
    gb = _gb(riemann, ricci, scalar, inverse)
    return index, radius, scalar_primal(scalar), scalar_primal(gb)


def run(
    *, point_count: int, workers: int, acceleration_source: str = "semidiscrete"
) -> dict[str, object]:
    if acceleration_source not in {"semidiscrete", "physical"}:
        raise ValueError("acceleration_source must be semidiscrete or physical")
    output = OUTPUT_DEV1 if acceleration_source == "semidiscrete" else OUTPUT_DEV2
    artifact_id = ARTIFACT_ID if acceleration_source == "semidiscrete" else ARTIFACT_ID.replace("DEV1", "DEV2")
    schema = SCHEMA if acceleration_source == "semidiscrete" else SCHEMA.replace("DEV1", "DEV2")
    if output.exists():
        raise RuntimeError(f"one-shot output already exists: {output}")
    if _sha(TRAPPED_ORIGIN) != EXPECTED_TRAPPED_ORIGIN_SHA256:
        raise RuntimeError("tracked trapped-origin compact result differs")
    grid, state, acceleration, p_r, q_r, steps, retries = _build_trapped_state(point_count)
    if acceleration_source == "physical":
        source = diagnose_gr0_proto12_accelerations_chunked(
            state.u[1:],
            state.p[1:],
            state.q[1:],
            p_r[1:],
            q_r[1:],
            grid.coordinates[1:],
            raw_tolerance=1.0e-12,
            condition_number_maximum=1.0e10,
            maximum_refinement_iterations=16,
            point_batch_size=2048,
        )
        physical = np.empty_like(acceleration)
        physical[1:] = source.accelerations
        physical[0] = regular_center_acceleration_limit(physical[1:5]).acceleration
        acceleration = physical
    observables = radial_null_observables(state)
    radii = grid.coordinates
    measured = np.flatnonzero((radii > 0.0) & (radii <= MEASUREMENT_RADIUS))
    inputs = [
        _curvature_input(index, radii, state.u, state.p, state.q, acceleration, p_r, q_r)
        for index in measured
    ]
    context = multiprocessing.get_context("fork")
    with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
        rows = list(pool.map(_warped_scalars, inputs, chunksize=max(1, len(inputs) // (8 * workers))))
    values = np.asarray([[row[2], row[3], row[4]] for row in rows], dtype=np.float64)
    effective_mass = values[:, 2]
    trapped = (observables.theta_plus[measured] < 0.0) & (
        observables.theta_minus[measured] < 0.0
    )
    minimum_position = int(np.argmin(effective_mass))
    minimum_index = int(measured[minimum_position])
    trapped_positions = np.flatnonzero(trapped)
    if not trapped_positions.size:
        raise RuntimeError("reconstructed target slice is not trapped")
    minimum_trapped_position = int(trapped_positions[np.argmin(effective_mass[trapped])])
    minimum_trapped_index = int(measured[minimum_trapped_position])
    check_positions = sorted(
        {
            minimum_position,
            minimum_trapped_position,
            int(trapped_positions[0]),
            int(trapped_positions[-1]),
        }
    )
    direct = [_direct_scalars(inputs[position]) for position in check_positions]
    direct_by_index = {row[0]: row for row in direct}
    route_checks = []
    for position in check_positions:
        warped = rows[position]
        exact = direct_by_index[warped[0]]
        route_checks.append(
            {
                "index": int(warped[0]),
                "radius": float(warped[1]),
                "R_warped": float(warped[2]),
                "R_direct": float(exact[2]),
                "R_absolute_difference": abs(float(warped[2]) - float(exact[2])),
                "GB_warped": float(warped[3]),
                "GB_direct": float(exact[3]),
                "GB_absolute_difference": abs(float(warped[3]) - float(exact[3])),
            }
        )
    negative = effective_mass < 0.0
    exterior = ~trapped
    result = {
        "artifact_id": artifact_id,
        "schema": schema,
        "source_trapped_origin_sha256": EXPECTED_TRAPPED_ORIGIN_SHA256,
        "background": {
            "branch": "GR-0",
            "amplitude": AMPLITUDE,
            "method": "RK4",
            "point_count": point_count,
            "coordinate_time_hex": TARGET_TIME.hex(),
            "accepted_steps": steps,
            "cfl_retries": retries,
            "state_sha256": array_content_sha256(state.u, state.p, state.q),
            "acceleration_source": (
                "PROTO12_SRC4_physical_source_without_KO_dissipation"
                if acceleration_source == "physical"
                else "PROTO12_semidiscrete_rhs_including_KO_dissipation"
            ),
        },
        "candidate_action": {
            "branch": "FGC-QR",
            "planck_mass": 2.0,
            "mu": 3.0,
            "g4": 0.5,
            "beta": -0.25,
            "eta": 0.5,
            "alpha": 0.0,
            "linear_trigger_formula": "m_lin_squared=9+R/4-GB/8",
        },
        "measurement": {
            "radius_maximum": MEASUREMENT_RADIUS,
            "measured_point_count": len(rows),
            "trapped_point_count": int(np.sum(trapped)),
            "negative_m_lin_squared_point_count": int(np.sum(negative)),
            "negative_m_lin_squared_trapped_point_count": int(np.sum(negative & trapped)),
            "negative_m_lin_squared_nontrapped_point_count": int(np.sum(negative & exterior)),
            "minimum_m_lin_squared": float(effective_mass[minimum_position]),
            "minimum_m_lin_squared_index": minimum_index,
            "minimum_m_lin_squared_radius": float(radii[minimum_index]),
            "minimum_trapped_m_lin_squared": float(effective_mass[minimum_trapped_position]),
            "minimum_trapped_m_lin_squared_index": minimum_trapped_index,
            "minimum_trapped_m_lin_squared_radius": float(radii[minimum_trapped_index]),
            "minimum_R": float(np.min(values[:, 0])),
            "maximum_R": float(np.max(values[:, 0])),
            "minimum_GB": float(np.min(values[:, 1])),
            "maximum_GB": float(np.max(values[:, 1])),
            "theta_plus_at_minimum_trapped_m_lin_squared": float(
                observables.theta_plus[minimum_trapped_index]
            ),
            "theta_minus_at_minimum_trapped_m_lin_squared": float(
                observables.theta_minus[minimum_trapped_index]
            ),
            "route_checks": route_checks,
        },
        "classification": (
            "fixed_action_curvature_trigger_present_in_trapped_region"
            if np.any(negative & trapped)
            else "fixed_action_curvature_trigger_absent_in_trapped_region"
        ),
        "terminal": True,
        "nonclaims": {
            "nonlinear_FGCQR_trajectory_executed": False,
            "regulator_activation_observed": False,
            "metric_null_defocusing_demonstrated": False,
            "child_domain_demonstrated": False,
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(canonical_json_bytes(result))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--point-count", type=int, default=4097)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--acceleration-source",
        choices=("semidiscrete", "physical"),
        default="semidiscrete",
    )
    args = parser.parse_args()
    result = run(
        point_count=args.point_count,
        workers=args.workers,
        acceleration_source=args.acceleration_source,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
