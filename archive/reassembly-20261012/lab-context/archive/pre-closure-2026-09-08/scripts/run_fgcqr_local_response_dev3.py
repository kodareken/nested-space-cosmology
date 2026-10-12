#!/usr/bin/env python3
"""Evaluate the complete FGC-QR source response at the trapped GR cells."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
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
from recursive_horizons.fgc.evolution.floating_jet import (  # noqa: E402
    FloatJet2,
    scalar_primal,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    radial_null_observables,
)
from recursive_horizons.fgc.evolution.nonlinear_source import (  # noqa: E402
    AccelerationSolveStop,
    NonlinearSolverConfig,
    float_ref1_residual,
    solve_accelerations,
)
from recursive_horizons.fgc.evolution.proto12_runtime import (  # noqa: E402
    diagnose_gr0_proto12_accelerations_chunked,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    regular_center_acceleration_limit,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.sgb1_ctl1_source import (  # noqa: E402
    HAT_NORMAL_FACTOR,
    REFERENCE_RADIAL_MINIMUM,
    TILDE_NORMAL_FACTOR,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    _state_from_adm_pg_jets,
)
from scripts.run_fgcqr_trigger_dev1 import (  # noqa: E402
    AMPLITUDE,
    TARGET_TIME,
    _build_trapped_state,
)


ARTIFACT_ID = "FGC-1-FGCQR-LOCAL-RESPONSE-A6042-DEV3"
SCHEMA = "FGC-1-FGCQR-LOCAL-RESPONSE-A6042-DEV3-v1"
OUTPUT = ROOT / "runs/fgc-2-sf1/fgcqr-local-response-a6042-dev3.json"
ACTION_CASES = (
    {
        "label": "original",
        "mu": 3.0,
        "beta": -0.25,
        "eta": 0.5,
        "g4": 0.5,
    },
    {
        "label": "mass_revision",
        "mu": 0.02,
        "beta": -0.25,
        "eta": 0.5,
        "g4": 0.5,
    },
    {
        "label": "curvature_and_stabilizer_revision",
        "mu": 0.02,
        "beta": -128.0,
        "eta": 16384.0,
        "g4": 256.0,
    },
)
FIELD_NAMES = ("alpha", "shift", "lambda", "areal_radius", "phi", "chi")
REFERENCE = flat_spherical_annulus_reference(
    radial_domain_minimum=REFERENCE_RADIAL_MINIMUM
)


def _fgc_state(case):
    (
        index,
        radius,
        action,
        values,
        dt,
        dr,
        dtt,
        dtr,
        drr,
    ) = case
    jets = {
        name: FloatJet2(
            values[column],
            dt[column],
            dr[column],
            dtt[column],
            dtr[column],
            drr[column],
        )
        for column, name in enumerate(FIELD_NAMES)
    }
    state = _state_from_adm_pg_jets(
        {
            "model_id": "FGC-QR",
            "action_parameters": {
                "planck_mass": 2,
                "scalar_mass": Fraction(str(action["mu"])),
                "quartic_coupling": Fraction(str(action["g4"])),
                "beta": Fraction(str(action["beta"])),
                "alpha_gb": 0,
                "eta": Fraction(str(action["eta"])),
            },
        },
        jets,
    )
    base_acceleration = np.asarray(
        [
            scalar_primal(state.h_tt.dtt),
            scalar_primal(state.h_tr.dtt),
            scalar_primal(state.h_rr.dtt),
            scalar_primal(state.areal_radius.dtt),
            scalar_primal(state.phi.dtt),
            scalar_primal(state.chi.dtt),
        ],
        dtype=np.float64,
    )
    return index, radius, action, state, base_acceleration


def _solve_case(case):
    index, radius, action, state, gr_acceleration = _fgc_state(case)
    config = NonlinearSolverConfig(
        residual_tolerance=1.0e-11,
        maximum_iterations=20,
        condition_number_maximum=1.0e12,
        normalized_branch_displacement_maximum=16.0,
        parameter_half_width=1.0 / 65536.0,
        acceleration_half_width=16.0,
        maximum_backtracks=32,
        minimum_step_fraction=2.0**-32,
    )
    before = float_ref1_residual(
        state,
        gr_acceleration,
        reference=REFERENCE,
        coordinate_radius=Fraction(radius),
        tilde_normal_factor=TILDE_NORMAL_FACTOR,
        hat_normal_factor=HAT_NORMAL_FACTOR,
    )
    try:
        solved = solve_accelerations(
            state,
            reference=REFERENCE,
            coordinate_radius=Fraction(radius),
            tilde_normal_factor=TILDE_NORMAL_FACTOR,
            hat_normal_factor=HAT_NORMAL_FACTOR,
            branch_center=state,
            warm_start=gr_acceleration,
            config=config,
        )
    except AccelerationSolveStop as stop:
        return {
            "index": index,
            "radius": radius,
            "action": action,
            "converged": False,
            "stop_reason": stop.reason,
            "stop_outcome": stop.outcome_label,
            "stop_diagnostics": stop.diagnostics,
            "GR_acceleration": gr_acceleration.tolist(),
            "FGC_residual_at_GR_acceleration": before.tolist(),
            "FGC_residual_at_GR_acceleration_infinity": float(
                np.linalg.norm(before, ord=np.inf)
            ),
        }
    acceleration = np.asarray(solved.accelerations, dtype=np.float64)
    return {
        "index": index,
        "radius": radius,
        "action": action,
        "converged": solved.converged,
        "iterations": solved.iterations,
        "jacobian_condition_infinity": solved.jacobian_condition_infinity,
        "branch_displacement_infinity": solved.branch_displacement_infinity,
        "GR_acceleration": gr_acceleration.tolist(),
        "FGC_acceleration": acceleration.tolist(),
        "FGC_minus_GR_acceleration": (acceleration - gr_acceleration).tolist(),
        "phi": float(case[3][4]),
        "phi_t": float(case[4][4]),
        "phi_r": float(case[5][4]),
        "GR_phi_acceleration": float(gr_acceleration[4]),
        "FGC_phi_acceleration": float(acceleration[4]),
        "signed_phi_growth_drive": float(case[3][4] * acceleration[4]),
        "FGC_residual_at_GR_acceleration": before.tolist(),
        "FGC_residual_at_GR_acceleration_infinity": float(
            np.linalg.norm(before, ord=np.inf)
        ),
        "solved_residual": list(solved.residual_vector),
        "solved_residual_infinity": solved.residual_infinity,
    }


def run(*, point_count: int, workers: int) -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError(f"one-shot output already exists: {OUTPUT}")
    grid, state, semidiscrete_acceleration, p_r, q_r, steps, retries = _build_trapped_state(
        point_count
    )
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
    acceleration = np.empty_like(semidiscrete_acceleration)
    acceleration[1:] = source.accelerations
    acceleration[0] = regular_center_acceleration_limit(
        acceleration[1:5]
    ).acceleration
    observables = radial_null_observables(state)
    radii = grid.coordinates
    trapped = np.flatnonzero(
        (radii > 0.0)
        & (radii <= 24.0)
        & (observables.theta_plus < 0.0)
        & (observables.theta_minus < 0.0)
    )
    if not trapped.size:
        raise RuntimeError("reconstructed state contains no trapped cell")
    cases = []
    for action in ACTION_CASES:
        for index in trapped:
            cases.append(
                (
                    int(index),
                    float(radii[index]),
                    action,
                    state.u[index].copy(),
                    state.p[index].copy(),
                    state.q[index].copy(),
                    acceleration[index].copy(),
                    p_r[index].copy(),
                    q_r[index].copy(),
                )
            )
    with ProcessPoolExecutor(
        max_workers=min(workers, len(cases)),
        mp_context=multiprocessing.get_context("fork"),
    ) as pool:
        responses = list(pool.map(_solve_case, cases))
    responses.sort(key=lambda item: (item["action"]["label"], item["index"]))
    by_label = {
        action["label"]: [
            row for row in responses if row["action"]["label"] == action["label"]
        ]
        for action in ACTION_CASES
    }
    selected = by_label["curvature_and_stabilizer_revision"]
    selected_converged = [row for row in selected if row["converged"]]
    result = {
        "artifact_id": ARTIFACT_ID,
        "schema": SCHEMA,
        "background": {
            "branch": "GR-0",
            "amplitude": AMPLITUDE,
            "method": "RK4",
            "point_count": point_count,
            "coordinate_time_hex": TARGET_TIME.hex(),
            "accepted_steps": steps,
            "cfl_retries": retries,
            "trapped_indices": [int(value) for value in trapped],
            "trapped_radii": [float(radii[value]) for value in trapped],
            "theta_plus": [float(observables.theta_plus[value]) for value in trapped],
            "theta_minus": [float(observables.theta_minus[value]) for value in trapped],
        },
        "action_cases": list(ACTION_CASES),
        "selected_action_design": {
            "planck_mass": 2.0,
            "alpha": 0.0,
            "selection": "beta_and_eta_chosen_to_make_both_measured_trapped_cells_order_one_tachyonic_and_g4_chosen_to_saturate_below_F_zero",
            "target_linear_mass_squared_order": -1.0,
            "effective_planck_zero_field_magnitude": float(np.sqrt(4.0 / 128.0)),
            "largest_estimated_quartic_saturation_field_magnitude": float(
                np.sqrt(1.5 / 256.0)
            ),
        },
        "acceleration_chart": "base_metric_dtt_after_ADM_to_metric_jet_map",
        "GR_acceleration_source": "PROTO12_SRC4_physical_source_without_KO_dissipation",
        "responses": responses,
        "classification": (
            "equation_directed_action_has_positive_complete_source_growth_drive"
            if selected_converged
            and all(row["signed_phi_growth_drive"] > 0.0 for row in selected_converged)
            else "equation_directed_action_did_not_establish_positive_complete_source_growth_drive"
        ),
        "summary": {
            label: {
                "converged_cell_count": sum(bool(row["converged"]) for row in rows),
                "positive_growth_cell_count": sum(
                    bool(row["converged"]) and row["signed_phi_growth_drive"] > 0.0
                    for row in rows
                ),
                "maximum_phi_acceleration": max(
                    (row["FGC_phi_acceleration"] for row in rows if row["converged"]),
                    default=None,
                ),
                "maximum_condition": max(
                    (row["jacobian_condition_infinity"] for row in rows if row["converged"]),
                    default=None,
                ),
            }
            for label, rows in by_label.items()
        },
        "terminal": True,
        "nonclaims": {
            "finite_time_activation_observed": False,
            "nonlinear_grid_trajectory_executed": False,
            "metric_null_defocusing_demonstrated": False,
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_bytes(canonical_json_bytes(result))
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--point-count", type=int, default=4097)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    result = run(point_count=args.point_count, workers=args.workers)
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
