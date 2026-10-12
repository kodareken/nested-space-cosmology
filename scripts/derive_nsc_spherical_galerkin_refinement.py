#!/usr/bin/env python3
"""Matched spatial refinement of the unchanged variational spherical coupling.

Reuses ``nsc_spherical_galerkin_coupling`` as it stands. Writes only
``results/development/nsc-spherical-coupling-refinement-v3.json``. Does not
filter, project, reset, relax a tolerance, extend the window to T=0.05, or
rewrite the v1 or v2 records.

The declared lines stay separate: full initial residuals at 1e-8, and the
short-window full residuals at 1e-3. A projected residual near 1e-10 is not
a pass of the full initial line.
"""
from __future__ import annotations

import hashlib
import json
import os
import signal
import time
import traceback
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np

from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results" / "development" / "nsc-spherical-coupling-refinement-v3.json"
V1 = LAB / "results" / "development" / "nsc-spherical-coupling-control-v1.json"
V2 = LAB / "results" / "development" / "nsc-spherical-coupling-control-v2.json"
DRIVER = Path(__file__).resolve()

EXPECTED = {
    "galerkin": "85d1f3dbbd82a38fe14d2405ad9782af6de68095052ff522ac6602e5d50a14ab",
    "tests": "e928b0d9b9d1162e3f665681adf1d4c9a3ff00c2941189d20bdda937c4b51a19",
    "coupling": "64e066b9b3210121461f902e748d4fc8a9cfc5978af2af47ab4b12fd8797e4f5",
    "v1": "88bea1c962478a2ec57e0c0f3844ab6ee9237d148dc1863f10122ee8fa5a9b60",
    "v2": "26030900805ebe27553e1b8d1cd490c6ca1a79c4631b65dad5483db3c27b370c",
}

NF = 256
NQ = 1024
DURATION = 0.005
DT = 0.0005
NF512_LIMIT_S = 120.0

CHECKPOINT_KEYS = (
    "time",
    "energy",
    "full_hamilton_max",
    "full_momentum_max",
    "projected_hamilton_max",
    "projected_momentum_max",
    "held_out_hamilton_max",
    "held_out_momentum_max",
    "lifted_proper_max",
    "chart_proper_max",
    "lifted_fieldwork",
    "unitarity",
    "positive_r",
    "positive_Q",
    "r_min_quadrature",
    "Q_min_quadrature",
)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def save(record):
    OUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUT.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(galerkin._jsonable(record), indent=2) + "\n")
    os.replace(temporary, OUT)


def hashes():
    return {
        "galerkin": sha256(galerkin._MODULE_PATH),
        "tests": sha256(galerkin._TEST_PATH),
        "coupling": sha256(galerkin._COUPLING_PATH),
        "action": sha256(galerkin._ACTION_PATH),
        "source": sha256(galerkin._SOURCE_PATH),
        "v1": sha256(V1),
        "v2": sha256(V2),
        "driver": sha256(DRIVER),
    }


def assert_predecessors(found):
    mismatched = [
        name for name in ("galerkin", "tests", "coupling", "v1", "v2")
        if found[name] != EXPECTED[name]
    ]
    if mismatched:
        raise RuntimeError("predecessor bytes differ from the recorded hashes: " + ", ".join(mismatched))


def v2_reference():
    payload = json.loads(V2.read_text())
    initial = payload["runs"]["initial_nf128_nq512"]
    window = payload["runs"]["nf128_nq512_T0.005"]
    return {
        "record": "results/development/nsc-spherical-coupling-control-v2.json",
        "nf": 128,
        "ng": 127,
        "nq": 512,
        "initial_projected_hamilton_max": initial["diagnostics"]["projected_hamilton_max"],
        "initial_full_hamilton_max": initial["diagnostics"]["full_hamilton_max"],
        "initial_full_momentum_max": initial["diagnostics"]["full_momentum_max"],
        "window_full_hamilton_max": window["full_hamilton_max"],
        "window_full_momentum_max": window["full_momentum_max"],
        "window_projected_hamilton_max": window["projected_hamilton_max"],
        "cited_window_full_hamilton": 0.03427,
        "cited_window_full_momentum": 0.00623,
        "projected_initial_does_not_pass_full_initial_criterion": True,
    }


def timed(timings, name):
    class _Timer:
        def __enter__(self):
            self.wall = time.perf_counter()
            self.cpu = time.process_time()
            return self

        def __exit__(self, *_args):
            timings[name] = {
                "wall_seconds": time.perf_counter() - self.wall,
                "cpu_seconds": time.process_time() - self.cpu,
            }

    return _Timer()


def preparation_summary(loaded):
    preparation = loaded[2]
    return {
        "problems": list(loaded[3]),
        "phase_applied_to_odd_lobe": preparation.get("phase_applied_to_odd_lobe"),
        "phase_on_minus_spinor_column": preparation.get("phase_on_minus_spinor_column"),
        "column_orthonormal_defect": preparation.get("column_orthonormal_defect"),
        "complement_occupation": preparation.get("complement_occupation"),
        "vacuum_identification": preparation.get("vacuum_identification"),
        "sea_subtracted": preparation.get("sea_subtracted"),
    }


def admissibility(grid, state, evolution=None):
    fine = galerkin.prolong_state(grid, state)
    names = ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")
    finite = all(np.isfinite(getattr(fine, name)).all() for name in names)
    report = {
        "positive_r": bool(np.min(fine.r) > 0),
        "positive_Q": bool(np.min(fine.Q) > 0),
        "finite_state": bool(finite),
        "r_min": float(np.min(fine.r)),
        "r_max": float(np.max(fine.r)),
        "Q_min": float(np.min(fine.Q)),
        "Q_max": float(np.max(fine.Q)),
    }
    if evolution is None:
        report["admissible"] = bool(report["positive_r"] and report["positive_Q"] and report["finite_state"])
        return report
    samples = evolution["samples"]
    report.update({
        "window_positive_r": all(bool(sample["positive_r"]) for sample in samples),
        "window_positive_Q": all(bool(sample["positive_Q"]) for sample in samples),
        "window_r_min": float(min(sample["r_min_quadrature"] for sample in samples)),
        "window_Q_min": float(min(sample["Q_min_quadrature"] for sample in samples)),
        "window_finite_energy": all(np.isfinite(sample["energy"]) for sample in samples),
        "unitarity_max": evolution.get("unitarity_max"),
        "chart_stopped": bool(evolution["stopped"]),
        "chart_reason": evolution["reason"],
    })
    report["admissible"] = bool(
        report["window_positive_r"]
        and report["window_positive_Q"]
        and report["window_finite_energy"]
        and report["finite_state"]
        and not evolution["stopped"]
    )
    return report


def over_tolerance(items, tolerance):
    failed = []
    for name, value in items:
        if value is None or not (value <= tolerance):
            failed.append({
                "name": name,
                "value": None if value is None else float(value),
                "tolerance": float(tolerance),
                "ratio_to_tolerance": None if value is None or not tolerance else float(value / tolerance),
            })
    return failed


def initial_axis(info, doubled):
    tolerance = galerkin.TOL_INITIAL_CONSTRAINT
    diagnostics = info.get("diagnostics") or {}
    projected_items = (
        ("projected_hamilton_max", diagnostics.get("projected_hamilton_max")),
        ("projected_momentum_max", diagnostics.get("projected_momentum_max")),
    )
    full_items = (
        ("full_hamilton_max", diagnostics.get("full_hamilton_max")),
        ("full_momentum_max", diagnostics.get("full_momentum_max")),
        ("held_out_hamilton_max", diagnostics.get("held_out_hamilton_max")),
        ("held_out_momentum_max", diagnostics.get("held_out_momentum_max")),
    )
    doubled_items = ()
    if doubled is not None:
        doubled_items = (
            ("doubled_full_hamilton_max", doubled["full_hamilton_max"]),
            ("doubled_full_momentum_max", doubled["full_momentum_max"]),
            ("doubled_held_out_hamilton_max", doubled["held_out_hamilton_max"]),
            ("doubled_held_out_momentum_max", doubled["held_out_momentum_max"]),
        )
    local_pass = bool(info.get("converged") and galerkin._local_checks_pass(info))
    doubled_pass = doubled is not None and (
        doubled["full_hamilton_max"] <= tolerance and doubled["full_momentum_max"] <= tolerance
    )
    projected_within = info.get("converged") and not over_tolerance(projected_items, tolerance)
    failures = []
    if not info.get("converged"):
        failures.append({
            "name": "initial_solve",
            "value": info.get("residual_max"),
            "tolerance": galerkin.TOL_NEWTON,
            "blocker": info.get("blocker"),
        })
    failures.extend(over_tolerance(full_items, tolerance))
    failures.extend(over_tolerance(projected_items, tolerance))
    failures.extend(over_tolerance(doubled_items[:2], tolerance))
    return {
        "tolerance": tolerance,
        "counts_projected_success_as_full_pass": False,
        "projected_within_tolerance": bool(projected_within),
        "full_local_pass": local_pass,
        "doubled_quadrature_full_pass": bool(doubled_pass),
        "pass": bool(local_pass and doubled_pass),
        "failures": failures,
        "full_items_reported": {name: value for name, value in full_items + doubled_items},
    }


def peak(samples, key):
    return float(max(sample[key] for sample in samples))


def window_view(result):
    samples = result["samples"]
    tolerance = galerkin.TOL_CONSTRAINT_VALIDATION
    energy_initial = float(samples[0]["energy"])
    energy_drift = result.get("energy_drift")
    if energy_drift is None:
        energy_drift = float(max(abs(sample["energy"] - energy_initial) for sample in samples))
    view = {
        "stopped": result["stopped"],
        "reason": result["reason"],
        "time": result["time"],
        "dt": result["dt"],
        "steps": result["steps"],
        "success_interval_completed": result["success"],
        "full_hamilton_max": peak(samples, "full_hamilton_max"),
        "full_momentum_max": peak(samples, "full_momentum_max"),
        "projected_hamilton_max": peak(samples, "projected_hamilton_max"),
        "projected_momentum_max": peak(samples, "projected_momentum_max"),
        "held_out_hamilton_max": peak(samples, "held_out_hamilton_max"),
        "held_out_momentum_max": peak(samples, "held_out_momentum_max"),
        "lifted_proper_max": peak(samples, "lifted_proper_max"),
        "chart_proper_max": peak(samples, "chart_proper_max"),
        "lifted_proper_initial": float(samples[0]["lifted_proper_max"]),
        "chart_proper_initial": float(samples[0]["chart_proper_max"]),
        "energy_initial": energy_initial,
        "energy_final": float(samples[-1]["energy"]),
        "energy_drift": energy_drift,
        "fieldwork_integral": result["fieldwork_integral"],
        "unitarity_max": result.get("unitarity_max"),
        "number_drift": result.get("number_drift"),
        "checkpoints": [{key: sample[key] for key in CHECKPOINT_KEYS} for sample in samples],
    }
    items = (
        ("full_hamilton_max", view["full_hamilton_max"]),
        ("full_momentum_max", view["full_momentum_max"]),
        ("held_out_hamilton_max", view["held_out_hamilton_max"]),
    )
    failures = over_tolerance(items, tolerance)
    if result["stopped"] or not result["success"]:
        failures.insert(0, {
            "name": "window_completion",
            "value": result["time"],
            "tolerance": DURATION,
            "blocker": result["reason"],
        })
    view["tolerance"] = tolerance
    view["pass"] = bool(galerkin._window_passes(result))
    view["failures"] = failures
    scale = max(1.0, abs(view["energy_initial"]))
    view["energy_relative_drift"] = float(view["energy_drift"] / scale) if view["energy_drift"] is not None else None
    view["energy_within_existing_relative_tolerance"] = bool(
        view["energy_drift"] is not None and view["energy_drift"] <= galerkin.TOL_ENERGY_RELATIVE * scale
    )
    return view


def ratio(new, old):
    if old == 0:
        return None
    return float(new / old)


def improvement(reference, initial_full, window_dt):
    pairs = {
        "initial_full_hamilton": (
            reference["initial_full_hamilton_max"],
            initial_full,
        ),
        "window_full_hamilton": (
            reference["window_full_hamilton_max"],
            window_dt["full_hamilton_max"],
        ),
        "window_full_momentum": (
            reference["window_full_momentum_max"],
            window_dt["full_momentum_max"],
        ),
    }
    report = {}
    for name, (previous, current) in pairs.items():
        report[name] = {
            "v2_nf128": previous,
            "nf256": current,
            "ratio_256_over_128": ratio(current, previous),
            "decreased": bool(current < previous),
        }
    return report


def classify_spatial(previous_ratio, next_ratio, next_value, tolerance):
    if next_value <= tolerance:
        return "nf512_initial_full_hamilton_meets_1e-8"
    if next_ratio is None:
        return "not_comparable"
    if next_ratio > 0.7:
        return "reduction_stalled_on_this_doubling"
    if previous_ratio is not None and next_ratio <= max(2.0 * previous_ratio, 0.25):
        return "still_decreasing_at_a_comparable_rate"
    if next_ratio < 0.5:
        return "still_decreasing_slower_than_the_previous_doubling"
    return "decreasing_without_a_stable_rate"


def solve_initial(fermions, quadrature, timings, prefix):
    with timed(timings, prefix + "_columns"):
        loaded = galerkin.load_physical_columns(fermions)
    with timed(timings, prefix + "_grid"):
        grid = galerkin.build_grid(fermions, quadrature=quadrature)
    state = galerkin.blank_state(grid, loaded[0], loaded[1])
    with timed(timings, prefix + "_newton"):
        solved, info = galerkin.solve_initial_radius(grid, state)
    return loaded, grid, solved, info


def estimate_nf512(timings):
    grid_1024 = timings["nf256_grid"]["wall_seconds"]
    grid_2048 = timings["nf256_doubled_grid"]["wall_seconds"]
    newton = timings["nf256_newton"]["wall_seconds"]
    evolve = timings.get("nf256_evolve_dt", {}).get("wall_seconds", 0.0)
    observed = grid_2048 / max(grid_1024, 1e-9)
    scale = max(8.0, observed)
    # Same nq=2048 operator build as the doubled grid, plus a larger odd subspace.
    grid_estimate = 1.25 * grid_2048
    newton_estimate = scale * newton
    evolve_estimate = evolve * (2048 / 1024) ** 2
    total = grid_estimate + newton_estimate + evolve_estimate
    return {
        "estimated_seconds": total,
        "grid_seconds": grid_estimate,
        "newton_seconds": newton_estimate,
        "one_window_seconds": evolve_estimate,
        "nq_cost_ratio_observed": observed,
        "nq_cost_ratio_used": scale,
        "basis": (
            "nq^3 operator build measured from nq=1024 to nq=2048, "
            "with the arithmetic floor 8, plus one T=0.005 matvec window scaled by nq^2"
        ),
        "cheap_under_120s": bool(total < NF512_LIMIT_S),
        "includes_dt_half": False,
        "includes_quadrupled_quadrature": False,
    }


def run_nf512(timings):
    loaded, grid, solved, info = solve_initial(512, 2048, timings, "nf512")
    diagnostics = info.get("diagnostics")
    velocities = info.get("velocities")
    energy0 = None
    balance = None
    evolution = None
    window = None
    if info.get("converged"):
        with timed(timings, "nf512_observables"):
            energy0 = float(galerkin.energy(grid, solved))
            balance = galerkin.work_balance(grid, solved)
        with timed(timings, "nf512_evolve_dt"):
            evolution = galerkin.evolve(grid, solved, DURATION, DT)
        window = window_view(evolution)
    return {
        "nf": 512,
        "ng": 511,
        "nq": 2048,
        "purpose": "one matched comparison of the full initial residual and one T=0.005 window",
        "dt_half_not_repeated": True,
        "quadrupled_quadrature_not_run": True,
        "preparation": preparation_summary(loaded),
        "converged": bool(info.get("converged")),
        "blocker": info.get("blocker"),
        "r_min_coarse": info.get("r_min_coarse"),
        "r_max_coarse": info.get("r_max_coarse"),
        "bracket_positive": info.get("bracket_positive"),
        "diagnostics": diagnostics,
        "velocities": velocities,
        "lifted_proper_max": None if not velocities else velocities["lifted_proper_max"],
        "chart_proper_max": None if not velocities else velocities["chart_proper_max"],
        "initial_admissibility": admissibility(grid, solved) if info.get("converged") else None,
        "initial_axis": initial_axis(info, None),
        "energy": energy0,
        "work_balance": balance,
        "window_dt": window,
        "isometry_geometry": grid.isometry_geometry,
        "isometry_columns": grid.isometry_columns,
    }


class _Nf512Timeout(Exception):
    pass


def _alarm(_signum, _frame):
    raise _Nf512Timeout("nf512 comparison exceeded 120s")


def seal(record, started, cpu, before):
    after = hashes()
    record["hashes_after"] = after
    record["v1_bytes_preserved"] = before["v1"] == after["v1"] == EXPECTED["v1"]
    record["v2_bytes_preserved"] = before["v2"] == after["v2"] == EXPECTED["v2"]
    record["module_unchanged"] = before["galerkin"] == after["galerkin"] == EXPECTED["galerkin"]
    record["runtime_seconds"] = time.perf_counter() - started
    record["cpu_seconds"] = time.process_time() - cpu
    record.setdefault("renewal", False)
    return record


def campaign():
    started = time.perf_counter()
    cpu = time.process_time()
    before = hashes()
    assert_predecessors(before)
    reference = v2_reference()
    timings = {}
    record = {
        "schema": "NSC-SPHERICAL-COUPLING-REFINEMENT-v3",
        "status": "PROVISIONAL",
        "renewal": False,
        "global_regeneration": False,
        "absolute_vacuum_claim": False,
        "module_unchanged": True,
        "filtered_v1_radius": False,
        "post_step_filter": False,
        "constraint_projection": False,
        "prescribed_radius": False,
        "threshold_relaxation": False,
        "time_extension_T_0_05": False,
        "v1_bytes_preserved": True,
        "v2_bytes_preserved": True,
        "representation": "variational Fourier-Galerkin matched spatial refinement",
        "question": (
            "Does matched refinement ng=nf-1, nq=4*nf bring the full quadrature "
            "constraints under the declared initial and window tolerances?"
        ),
        "phi_dot": "U_f† (-i H_fine U_f Phi)",
        "thread_limits": {
            "OMP_NUM_THREADS": os.environ["OMP_NUM_THREADS"],
            "OPENBLAS_NUM_THREADS": os.environ["OPENBLAS_NUM_THREADS"],
            "MKL_NUM_THREADS": os.environ["MKL_NUM_THREADS"],
            "NUMEXPR_NUM_THREADS": os.environ["NUMEXPR_NUM_THREADS"],
            "VECLIB_MAXIMUM_THREADS": os.environ["VECLIB_MAXIMUM_THREADS"],
        },
        "tolerances": {
            "initial_full": galerkin.TOL_INITIAL_CONSTRAINT,
            "window_full": galerkin.TOL_CONSTRAINT_VALIDATION,
            "energy_relative": galerkin.TOL_ENERGY_RELATIVE,
            "newton": galerkin.TOL_NEWTON,
            "projected_residual_satisfies_full_initial_criterion": False,
        },
        "hashes_before": before,
        "v2_reference": reference,
        "costs": timings,
    }
    save(record)

    loaded, grid, solved, info = solve_initial(NF, NQ, timings, "nf256")
    record["preparation"] = preparation_summary(loaded)
    if grid.ng != 255 or grid.nf != 256 or grid.nq != 1024:
        raise RuntimeError(f"unexpected grid {(grid.ng, grid.nf, grid.nq)}")
    omega, nyquist = (None, None)
    if info.get("converged"):
        omega, nyquist = galerkin.subspace_frequency(grid, solved)
    suggested, _omega_check, _nyquist_check = (
        galerkin.stable_timestep(grid, solved) if info.get("converged") else (None, None, None)
    )
    initial = {
        "nf": NF,
        "ng": grid.ng,
        "nq": NQ,
        "doubled_nq": 2 * NQ,
        "converged": bool(info.get("converged")),
        "blocker": info.get("blocker"),
        "newton_steps": info.get("steps"),
        "rho_independent_of_r_max": info.get("rho_independent_of_r_max"),
        "r_min_coarse": info.get("r_min_coarse"),
        "r_max_coarse": info.get("r_max_coarse"),
        "bracket_min": info.get("bracket_min"),
        "bracket_positive": info.get("bracket_positive"),
        "imposed_radius": False,
        "filtered_v1_radius": False,
        "diagnostics": info.get("diagnostics"),
        "velocities": info.get("velocities"),
        "lifted_proper_max": None if not info.get("velocities") else info["velocities"]["lifted_proper_max"],
        "chart_proper_max": None if not info.get("velocities") else info["velocities"]["chart_proper_max"],
        "isometry_geometry": grid.isometry_geometry,
        "isometry_columns": grid.isometry_columns,
        "omega_subspace": omega,
        "omega_quadrature_nyquist": nyquist,
        "module_stable_dt": suggested,
        "requested_dt": DT,
        "requested_dt_half": 0.5 * DT,
    }
    record["nf256_initial_solve"] = initial
    save(record)
    if loaded[3]:
        record["verdict"] = "BLOCKED_UNCORRECTED_PREPARATION"
        record["exact_failure"] = list(loaded[3])
        return seal(record, started, cpu, before)
    if not info.get("converged"):
        record["initial_axis"] = initial_axis(info, None)
        record["verdict"] = "FAIL_INITIAL_SOLVE"
        record["exact_failure"] = record["initial_axis"]["failures"]
        return seal(record, started, cpu, before)

    with timed(timings, "nf256_initial_observables"):
        initial["admissibility"] = admissibility(grid, solved)
        initial["energy"] = float(galerkin.energy(grid, solved))
        initial["work_balance"] = galerkin.work_balance(grid, solved)
    with timed(timings, "nf256_doubled_grid"):
        doubled_grid = galerkin.build_grid(NF, quadrature=2 * NQ)
    with timed(timings, "nf256_doubled_residual"):
        doubled = galerkin.constraint_diagnostics(doubled_grid, solved)
    initial["doubled_quadrature"] = doubled
    initial["doubled_minus_base_full_hamilton"] = float(
        doubled["full_hamilton_max"] - info["diagnostics"]["full_hamilton_max"]
    )
    record["initial_axis"] = initial_axis(info, doubled)
    save(record)

    with timed(timings, "nf256_evolve_dt"):
        primary = galerkin.evolve(grid, solved, DURATION, DT)
    with timed(timings, "nf256_evolve_dt_half"):
        halved = galerkin.evolve(grid, solved, DURATION, 0.5 * DT)
    window_dt = window_view(primary)
    window_half = window_view(halved)
    window_dt["admissibility"] = admissibility(grid, primary["state"], primary)
    window_half["admissibility"] = admissibility(grid, halved["state"], halved)
    window_pass = bool(window_dt["pass"] and window_half["pass"])
    record["window_axis"] = {
        "tolerance": galerkin.TOL_CONSTRAINT_VALIDATION,
        "duration": DURATION,
        "pass": window_pass,
        "dt": window_dt,
        "dt_half": window_half,
        "failures": [{"dt": DT, **item} for item in window_dt["failures"]]
        + [{"dt": 0.5 * DT, **item} for item in window_half["failures"]],
    }
    record["improvement_vs_v2_nf128"] = improvement(
        reference,
        info["diagnostics"]["full_hamilton_max"],
        window_dt,
    )
    save(record)

    initial_pass = bool(record["initial_axis"]["pass"])
    exact = []
    if not initial_pass:
        exact.extend({"axis": "initial", **item} for item in record["initial_axis"]["failures"])
    if not window_pass:
        exact.extend({"axis": "window", **item} for item in record["window_axis"]["failures"])
    chart_stop = primary["reason"] or halved["reason"]
    if chart_stop is not None:
        verdict = "STOP_POSITIVE_CHART"
    elif not initial_pass and not window_pass:
        verdict = "FAIL_INITIAL_AND_WINDOW"
    elif not initial_pass:
        verdict = "FAIL_INITIAL_FULL_CONSTRAINT"
    elif not window_pass:
        verdict = "FAIL_WINDOW_CONSTRAINT"
    else:
        verdict = "PASS_DECLARED_TOLERANCES_ON_SHORT_WINDOW"
    record["verdict"] = verdict
    record["exact_failure"] = exact
    record["renewal"] = False
    save(record)

    failed = verdict != "PASS_DECLARED_TOLERANCES_ON_SHORT_WINDOW"
    comparison = {"ran": False}
    if not failed:
        comparison["reason"] = "nf256 met both declared constraint tolerances; nf512 was not required"
    else:
        comparison["estimate"] = estimate_nf512(timings)
        comparison["reason_required"] = verdict
        if not comparison["estimate"]["cheap_under_120s"]:
            comparison["reason"] = (
                "nf512 comparison estimated at "
                f"{comparison['estimate']['estimated_seconds']:.3f}s, at or above 120s; not run"
            )
        else:
            signal.signal(signal.SIGALRM, _alarm)
            signal.alarm(int(NF512_LIMIT_S))
            try:
                comparison["result"] = run_nf512(timings)
                comparison["ran"] = True
                comparison["reason"] = "estimate was below 120s"
            except _Nf512Timeout as error:
                comparison["ran"] = False
                comparison["reason"] = str(error)
            finally:
                signal.alarm(0)
            result = comparison.get("result")
            if result and result.get("diagnostics"):
                full_512 = result["diagnostics"]["full_hamilton_max"]
                full_256 = info["diagnostics"]["full_hamilton_max"]
                full_128 = reference["initial_full_hamilton_max"]
                previous = ratio(full_256, full_128)
                following = ratio(full_512, full_256)
                comparison["spatial_reading"] = {
                    "initial_full_hamilton_nf128": full_128,
                    "initial_full_hamilton_nf256": full_256,
                    "initial_full_hamilton_nf512": full_512,
                    "ratio_256_over_128": previous,
                    "ratio_512_over_256": following,
                    "classification": classify_spatial(
                        previous, following, full_512, galerkin.TOL_INITIAL_CONSTRAINT,
                    ),
                }
                if result.get("window_dt"):
                    comparison["window_reading"] = {
                        "full_hamilton_max": result["window_dt"]["full_hamilton_max"],
                        "full_momentum_max": result["window_dt"]["full_momentum_max"],
                        "projected_hamilton_max": result["window_dt"]["projected_hamilton_max"],
                        "projected_momentum_max": result["window_dt"]["projected_momentum_max"],
                        "held_out_hamilton_max": result["window_dt"]["held_out_hamilton_max"],
                        "lifted_proper_max": result["window_dt"]["lifted_proper_max"],
                        "pass": result["window_dt"]["pass"],
                        "dt_half_not_repeated": True,
                    }
    record["nf512_comparison"] = comparison
    record["unresolved_physical_embedding"] = {
        "frozen_six_mode_generator_embedded": False,
        "link_B_mapping_closed": False,
        "common_action_closure": False,
        "counted_action": "Gamma_one",
        "vacuum_matched_ctp_used": False,
        "incoming_gate": "OPEN",
        "incoming_gate_campaign": "paused",
        "renewal": False,
        "statement": (
            "Matched refinement tests the existing variational subspace only. "
            "Exact physical embedding of the frozen six-mode generator remains incomplete: "
            "site copy into this chart is still obstructed, link B is unmatched, and there is "
            "no common-action closure. The counted action remains Gamma_one. "
            "The vacuum-matched CTP branch was not used. The incoming gate stays OPEN and its "
            "campaign stays paused. A short window is not a continuum limit and not a renewal."
        ),
    }
    return seal(record, started, cpu, before)


def main():
    record = {
        "schema": "NSC-SPHERICAL-COUPLING-REFINEMENT-v3",
        "status": "PROVISIONAL",
        "renewal": False,
    }
    try:
        record = campaign()
    except Exception as error:
        record["verdict"] = record.get("verdict") or "RUN_FAILED"
        record["run_error"] = {
            "type": type(error).__name__,
            "message": str(error),
            "traceback": traceback.format_exc(),
        }
        save(record)
        raise
    save(record)
    preserved = hashes()
    if preserved["v1"] != EXPECTED["v1"] or preserved["v2"] != EXPECTED["v2"]:
        raise RuntimeError("v1 or v2 bytes changed after the record write")
    if preserved["galerkin"] != EXPECTED["galerkin"]:
        raise RuntimeError("Galerkin module bytes changed during the run")
    summary = {
        "verdict": record.get("verdict"),
        "initial_pass": (record.get("initial_axis") or {}).get("pass"),
        "window_pass": (record.get("window_axis") or {}).get("pass"),
        "exact_failure": record.get("exact_failure"),
        "nf512": {
            "ran": (record.get("nf512_comparison") or {}).get("ran"),
            "reason": (record.get("nf512_comparison") or {}).get("reason"),
            "spatial_reading": (record.get("nf512_comparison") or {}).get("spatial_reading"),
        },
        "runtime_seconds": record.get("runtime_seconds"),
        "cpu_seconds": record.get("cpu_seconds"),
        "record_sha256": sha256(OUT),
        "v1_sha256": preserved["v1"],
        "v2_sha256": preserved["v2"],
        "galerkin_sha256": preserved["galerkin"],
    }
    print(json.dumps(galerkin._jsonable(summary), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
