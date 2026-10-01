#!/usr/bin/env python3
"""Diagnostic measurement of the full-source nf256 state.

Constructs the same best full-source radius as /tmp/nf256_compensate_probe.py:
prolong the converged nf128 radius, Newton the full nf256 source, and keep the
iterate with the smallest projected residual. A projected residual near 1e-9 is
recorded as the precision limit of that iteration. It is not a pass of the
Newton 1e-10 line or of the full initial 1e-8 line, and it does not stop the
measurement.

Writes only the v5 JSON and NPZ. Does not filter, retune, extend to T=0.05,
or claim a renewal, a continuum limit, or a frozen-B embedding.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import time
import traceback
from pathlib import Path

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np

try:
    from threadpoolctl import threadpool_limits

    threadpool_limits(limits=1).__enter__()
    THREADPOOL_LIMIT = 1
except Exception as error:
    THREADPOOL_LIMIT = "not applied: " + type(error).__name__

from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results" / "development" / "nsc-spherical-coupling-refinement-v5.json"
NPZ = LAB / "results" / "development" / "nsc-spherical-coupling-refinement-v5.npz"
DRIVER = Path(__file__).resolve()
V1 = LAB / "results" / "development" / "nsc-spherical-coupling-control-v1.json"
V2 = LAB / "results" / "development" / "nsc-spherical-coupling-control-v2.json"
V3_DRIVER = LAB / "scripts" / "derive_nsc_spherical_galerkin_refinement.py"
V3 = LAB / "results" / "development" / "nsc-spherical-coupling-refinement-v3.json"
V4_DRIVER = LAB / "scripts" / "derive_nsc_spherical_galerkin_refinement_v4.py"
V4 = LAB / "results" / "development" / "nsc-spherical-coupling-refinement-v4.json"
LOCKED = LAB / "results" / "development" / "nsc-subgap-history-response.json"
PROBE = Path("/tmp/nf256_compensate_probe.py")
SEED_PROBE = Path("/tmp/nf256_seed_probe.py")

DURATION = 0.005
DT = 0.0005
DT_HALF = 0.00025
NF512_BUDGET_S = 120.0
# Printed v4 archive of the compensate probe. Comparison only, not a measurement.
ARCHIVE_FULL_C = 0.001768137526894975
ARCHIVE_PROJECTED_OPERATOR = 9.608390314276184e-10

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

BINDING_COEFFICIENTS = ("A", "C_W", "C_F", "C_E", "C_box", "flux", "V_rel")
BINDING_CALIBRATION = (
    "Mweight_executed",
    "a0",
    "beta0",
    "b0",
    "phase",
    "kappa",
    "Omega",
    "k",
    "ell",
    "length",
)


def sha256(path):
    digest = hashlib.sha256()
    file_path = Path(path)
    if not file_path.is_file():
        return None
    with file_path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def source_hashes():
    paths = {
        "driver": DRIVER,
        "galerkin": galerkin._MODULE_PATH,
        "tests": galerkin._TEST_PATH,
        "coupling": galerkin._COUPLING_PATH,
        "action": galerkin._ACTION_PATH,
        "source": galerkin._SOURCE_PATH,
        "locked_record": LOCKED,
        "v1": V1,
        "v2": V2,
        "v3_driver": V3_DRIVER,
        "v3_record": V3,
        "v4_driver": V4_DRIVER,
        "v4_record": V4,
        "compensate_probe": PROBE,
        "seed_probe": SEED_PROBE,
        "npz": NPZ,
    }
    return {name: sha256(path) for name, path in paths.items()}


def jsonable(value):
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if not math.isfinite(number):
            return {"nonfinite": True, "repr": repr(number)}
        return number
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, complex):
        return {"real": jsonable(value.real), "imag": jsonable(value.imag)}
    if value is None or isinstance(value, str):
        return value
    return str(value)


def save(record):
    OUT.parent.mkdir(parents=True, exist_ok=True)
    temporary = OUT.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(jsonable(record), indent=2) + "\n")
    os.replace(temporary, OUT)


def write_npz(arrays):
    NPZ.parent.mkdir(parents=True, exist_ok=True)
    temporary = NPZ.with_name(NPZ.name + ".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    os.replace(temporary, NPZ)


class Clock:
    def __init__(self, timings, name):
        self.timings = timings
        self.name = name

    def __enter__(self):
        self.wall = time.perf_counter()
        self.cpu = time.process_time()
        print("start", self.name, flush=True)
        return self

    def __exit__(self, *_args):
        self.timings[self.name] = {
            "wall_seconds": time.perf_counter() - self.wall,
            "cpu_seconds": time.process_time() - self.cpu,
        }
        print("done", self.name, self.timings[self.name], flush=True)


def embed_radius(radius, ng_coarse, ng_fine, length):
    modes_c = galerkin.geometry_modes(ng_coarse)
    modes_f = galerkin.geometry_modes(ng_fine)
    coarse = np.arange(ng_coarse, dtype=float) * (length / ng_coarse)
    fine = np.arange(ng_fine, dtype=float) * (length / ng_fine)
    analysis = np.exp(-2j * np.pi * coarse[:, None] * modes_c[None, :] / length) / ng_coarse
    coefficients = analysis.T @ np.asarray(radius, dtype=float)
    coeff_f = np.zeros(ng_fine, dtype=complex)
    lookup = {float(mode): index for index, mode in enumerate(modes_f)}
    for coefficient, mode in zip(coefficients, modes_c, strict=True):
        coeff_f[lookup[float(mode)]] = coefficient
    synthesis = np.exp(2j * np.pi * fine[:, None] * modes_f[None, :] / length)
    nodal = synthesis @ coeff_f
    high = np.abs(modes_f) > (ng_coarse // 2)
    return (
        np.real(nodal).astype(np.float64, copy=True),
        float(np.max(np.abs(np.imag(nodal)))),
        float(np.max(np.abs(coeff_f[high]))) if np.any(high) else 0.0,
    )


def best_full_source_radius(grid, radius, rho, rounds=6):
    """Compensate-probe acceptance: full step, keep the smallest projected residual."""
    history = []
    best_residual = None
    best_radius = np.asarray(radius, dtype=np.float64).copy()
    current = best_radius.copy()
    for iteration in range(rounds):
        projected, full = galerkin.projected_radius_operator(grid, current, rho)
        residual_max = float(np.max(np.abs(projected)))
        full_max = float(np.max(np.abs(full)))
        entry = {
            "iteration": iteration,
            "projected": residual_max,
            "full_operator": full_max,
        }
        if best_residual is None or residual_max < best_residual:
            best_residual = residual_max
            best_radius = current.copy()
            entry["retained_as_best"] = True
        else:
            entry["retained_as_best"] = False
        if residual_max < galerkin.TOL_NEWTON:
            entry["stop"] = "newton_tolerance"
            history.append(entry)
            break
        jacobian = galerkin._projected_radius_jacobian(grid, current)
        try:
            delta = np.linalg.solve(jacobian, -projected)
        except np.linalg.LinAlgError as error:
            entry["stop"] = "singular_jacobian"
            entry["error"] = str(error)
            history.append(entry)
            break
        delta_max = float(np.max(np.abs(delta)))
        entry["delta_max"] = delta_max
        trial = current + delta
        if not np.isfinite(trial).all() or np.min(galerkin.prolong_geometry(grid, trial)) <= 0.0:
            entry["stop"] = "left_positive_chart"
            history.append(entry)
            break
        trial_projected, _full_trial = galerkin.projected_radius_operator(grid, trial, rho)
        trial_max = float(np.max(np.abs(trial_projected)))
        entry["trial_projected"] = trial_max
        if trial_max < residual_max:
            current = trial
            entry["accepted"] = True
            history.append(entry)
        else:
            entry["accepted"] = False
            entry["stop"] = "no_decrease"
            history.append(entry)
            break
        if delta_max < 1e-12:
            break
    return best_radius, float(best_residual), history


def full_source_rho(grid, state):
    fine = galerkin.prolong_state(grid, state)
    source = galerkin.source_from_columns(grid.fine, fine)
    rho = source["force_L"] / grid.dx_q
    varied = state.copy()
    varied.r = state.r * 1.7
    other = galerkin.source_from_columns(grid.fine, galerkin.prolong_state(grid, varied))
    variation = float(np.max(np.abs(other["force_L"] - source["force_L"])))
    return rho, variation, float(np.max(np.abs(rho)))


def binding_arrays(grid):
    names = []
    values = []
    coefficients = grid.fine.coefficients
    calibration = grid.fine.calibration
    for name in BINDING_COEFFICIENTS:
        names.append("coefficient:" + name)
        values.append(float(coefficients[name]))
    for name in BINDING_CALIBRATION:
        names.append("calibration:" + name)
        values.append(float(calibration[name]))
    return np.array(names), np.asarray(values, dtype=np.float64)


def state_arrays(grid, state, rho, *, prefix, projected_residual, full_operator):
    names, values = binding_arrays(grid)
    return {
        prefix + "nf": np.int32(grid.nf),
        prefix + "ng": np.int32(grid.ng),
        prefix + "nq": np.int32(grid.nq),
        prefix + "length": np.float64(grid.length),
        prefix + "radius": np.asarray(state.r, dtype=np.float64),
        prefix + "geometry_Q": np.asarray(state.Q, dtype=np.float64),
        prefix + "geometry_r": np.asarray(state.r, dtype=np.float64),
        prefix + "geometry_chi": np.asarray(state.chi, dtype=np.float64),
        prefix + "geometry_p_Q": np.asarray(state.p_Q, dtype=np.float64),
        prefix + "geometry_p_r": np.asarray(state.p_r, dtype=np.float64),
        prefix + "geometry_p_chi": np.asarray(state.p_chi, dtype=np.float64),
        prefix + "geometry_modes": np.asarray(grid.modes_g, dtype=np.float64),
        prefix + "quadrature_radius": np.asarray(galerkin.prolong_geometry(grid, state.r), dtype=np.float64),
        prefix + "columns_phi0": np.asarray(state.phi0),
        prefix + "columns_phi1": np.asarray(state.phi1),
        prefix + "occupations": np.asarray(grid.fine.occupations, dtype=np.float64),
        prefix + "binding_names": names,
        prefix + "binding_values": values,
        prefix + "rho": np.asarray(rho, dtype=np.float64),
        prefix + "projected_residual": np.float64(projected_residual),
        prefix + "full_operator": np.float64(full_operator),
    }


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


def occupation_report(state):
    values = np.sort(np.real(galerkin.occupation_eigenvalues(state.phi0, state.phi1)))
    declared = np.sort(np.asarray(galerkin.OCCUPATIONS, dtype=float))
    leading = values[-declared.size :]
    return {
        "declared": declared.tolist(),
        "leading_eigenvalues": leading.tolist(),
        "declared_gap": float(np.max(np.abs(leading - declared))),
        "tail_max": float(np.max(np.abs(values[: -declared.size]))) if values.size > declared.size else 0.0,
    }


def admissibility(grid, state, evolution=None):
    fine = galerkin.prolong_state(grid, state)
    names = ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")
    finite = all(np.isfinite(getattr(fine, name)).all() for name in names)
    report = {
        "positive_r": bool(np.min(fine.r) > 0.0),
        "positive_Q": bool(np.min(fine.Q) > 0.0),
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
        "window_finite_energy": all(math.isfinite(float(sample["energy"])) for sample in samples),
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


def constraint_view(diagnostics, velocities):
    full_c = float(diagnostics["full_hamilton_max"])
    full_d = float(diagnostics["full_momentum_max"])
    return {
        "full_C_hamilton_max": full_c,
        "full_D_momentum_max": full_d,
        "projected_C_hamilton_max": float(diagnostics["projected_hamilton_max"]),
        "projected_D_momentum_max": float(diagnostics["projected_momentum_max"]),
        "held_out_C_hamilton_max": float(diagnostics["held_out_hamilton_max"]),
        "held_out_D_momentum_max": float(diagnostics["held_out_momentum_max"]),
        "positive_r": bool(diagnostics["positive_r"]),
        "positive_Q": bool(diagnostics["positive_Q"]),
        "r_min_quadrature": float(diagnostics["r_min_quadrature"]),
        "Q_min_quadrature": float(diagnostics["Q_min_quadrature"]),
        "hamilton_unresolved_fraction": float(diagnostics["hamilton_modes"]["unresolved_fraction"]),
        "hamilton_unresolved_max_mode": int(diagnostics["hamilton_modes"]["unresolved_max_mode"]),
        "momentum_unresolved_fraction": float(diagnostics["momentum_modes"]["unresolved_fraction"]),
        "lifted_normal_velocity_max": float(velocities["lifted_proper_max"]),
        "projected_normal_velocity_max": float(velocities["projected_proper_max"]),
        "chart_proper_velocity_max": float(velocities["chart_proper_max"]),
        "coordinate_velocity_max": float(velocities["coordinate_max"]),
        "shift_piece_max": float(velocities["shift_piece_max"]),
        "full_initial_1e-8_passed": False,
        "full_C_leq_1e-8": bool(full_c <= galerkin.TOL_INITIAL_CONSTRAINT),
        "full_D_leq_1e-8": bool(full_d <= galerkin.TOL_INITIAL_CONSTRAINT),
        "diagnostics": diagnostics,
        "velocities": velocities,
    }


def measure_constraints(grid, state):
    diagnostics = galerkin.constraint_diagnostics(grid, state)
    velocities = galerkin.normal_velocities(grid, state)
    return constraint_view(diagnostics, velocities)


def peak(samples, key):
    return float(max(sample[key] for sample in samples))


def window_report(grid, result):
    samples = result["samples"]
    energy_initial = float(samples[0]["energy"])
    if result.get("energy_drift") is None:
        energy_drift = float(max(abs(float(sample["energy"]) - energy_initial) for sample in samples))
    else:
        energy_drift = float(result["energy_drift"])
    scale = max(1.0, abs(energy_initial))
    full_c = peak(samples, "full_hamilton_max")
    full_d = peak(samples, "full_momentum_max")
    held_c = peak(samples, "held_out_hamilton_max")
    tolerance = float(galerkin.TOL_CONSTRAINT_VALIDATION)
    within = bool(
        result["success"]
        and not result["stopped"]
        and full_c <= tolerance
        and full_d <= tolerance
        and held_c <= tolerance
    )
    return {
        "diagnostic_state": True,
        "labeled_as_window_pass": False,
        "measured_within_window_1e-3": within,
        "stopped": bool(result["stopped"]),
        "reason": result["reason"],
        "time": float(result["time"]),
        "requested_duration": DURATION,
        "dt": float(result["dt"]),
        "steps": int(result["steps"]),
        "success_interval_completed": bool(result["success"]),
        "full_C_hamilton_max": full_c,
        "full_D_momentum_max": full_d,
        "projected_C_hamilton_max": peak(samples, "projected_hamilton_max"),
        "projected_D_momentum_max": peak(samples, "projected_momentum_max"),
        "held_out_C_hamilton_max": held_c,
        "held_out_D_momentum_max": peak(samples, "held_out_momentum_max"),
        "lifted_normal_velocity_max": peak(samples, "lifted_proper_max"),
        "chart_proper_velocity_max": peak(samples, "chart_proper_max"),
        "energy_initial": energy_initial,
        "energy_final": float(samples[-1]["energy"]),
        "energy_drift": energy_drift,
        "energy_relative_drift": float(energy_drift / scale),
        "energy_within_existing_relative_1e-6": bool(
            energy_drift <= galerkin.TOL_ENERGY_RELATIVE * scale
        ),
        "work_fieldwork_integral": float(result["fieldwork_integral"]),
        "unitarity_max": result.get("unitarity_max"),
        "number_drift": result.get("number_drift"),
        "positivity": {
            "all_positive_r": all(bool(sample["positive_r"]) for sample in samples),
            "all_positive_Q": all(bool(sample["positive_Q"]) for sample in samples),
            "r_min": float(min(sample["r_min_quadrature"] for sample in samples)),
            "Q_min": float(min(sample["Q_min_quadrature"] for sample in samples)),
        },
        "admissibility": admissibility(grid, result["state"], result),
        "checkpoints": [{key: sample[key] for key in CHECKPOINT_KEYS} for sample in samples],
    }


def run_window(grid, state, dt):
    try:
        result = galerkin.evolve(grid, state.copy(), DURATION, dt)
    except Exception as error:
        return {
            "diagnostic_state": True,
            "labeled_as_window_pass": False,
            "failed": True,
            "error": {
                "type": type(error).__name__,
                "message": str(error),
                "traceback": traceback.format_exc(),
            },
        }
    report = window_report(grid, result)
    report["failed"] = False
    return report


def safe_call(function, *args):
    try:
        return function(*args), None
    except Exception as error:
        return None, {
            "type": type(error).__name__,
            "message": str(error),
            "traceback": traceback.format_exc(),
        }


def blank_record():
    return {
        "schema": "NSC-SPHERICAL-COUPLING-REFINEMENT-v5",
        "status": "DIAGNOSTIC",
        "verdict": "RUNNING",
        "pass": False,
        "renewal": False,
        "continuum_limit_claimed": False,
        "frozen_B_embedding_claimed": False,
        "time_extension_T_0_05": False,
        "counts_projected_success_as_full_pass": False,
        "newton_1e-10_passed": False,
        "full_initial_1e-8_passed": False,
        "module_unchanged_during_run": None,
        "historical_v1_v4_unchanged": None,
        "state_label": "DIAGNOSTIC",
        "question": (
            "On the best full-source nf256 iterate, what are the full and projected "
            "constraints and the lifted normal velocity at nq=1024 and nq=2048, and "
            "what are the full C and D maxima, positivity, admissibility, energy drift, "
            "and work through T=0.005 at dt=0.0005 and dt=0.00025?"
        ),
        "construction": (
            "Prolong the converged nf128 radius onto ng=255 and Newton the full nf256 "
            "source with the compensate-probe rule: one full step, stop when the step "
            "does not decrease the projected residual, keep the best iterate."
        ),
        "phi_dot": "U_f† (-i H_fine U_f Phi)",
        "duration": DURATION,
        "dt": DT,
        "dt_half": DT_HALF,
        "thread_limits": {
            "OMP_NUM_THREADS": os.environ["OMP_NUM_THREADS"],
            "OPENBLAS_NUM_THREADS": os.environ["OPENBLAS_NUM_THREADS"],
            "MKL_NUM_THREADS": os.environ["MKL_NUM_THREADS"],
            "NUMEXPR_NUM_THREADS": os.environ["NUMEXPR_NUM_THREADS"],
            "VECLIB_MAXIMUM_THREADS": os.environ["VECLIB_MAXIMUM_THREADS"],
            "threadpoolctl": THREADPOOL_LIMIT,
        },
        "platform": {
            "np_longdouble_equals_float64": bool(np.dtype(np.longdouble) == np.dtype(np.float64)),
            "np_longdouble_eps": float(np.finfo(np.longdouble).eps),
        },
        "tolerances_read_from_source_not_retuned": {
            "initial_full": float(galerkin.TOL_INITIAL_CONSTRAINT),
            "window_full": float(galerkin.TOL_CONSTRAINT_VALIDATION),
            "energy_relative": float(galerkin.TOL_ENERGY_RELATIVE),
            "newton_internal": float(galerkin.TOL_NEWTON),
            "projected_residual_satisfies_full_initial_criterion": False,
            "newton_1e-10_passed": False,
            "full_initial_1e-8_passed": False,
        },
        "precision_limit": {
            "internal_newton_residual_near_1e-9_aborts_diagnostics": False,
            "statement": (
                "The internal projected residual near 1e-9 is a reported precision "
                "limit of this Newton iteration. It is not labeled as a pass of "
                "Newton 1e-10 or of the full initial residual 1e-8."
            ),
        },
        "open_claims": {
            "continuum_limit": "open",
            "frozen_B_embedding": "open",
            "renewal": False,
            "vacuum_matched_ctp_used": False,
            "incoming_gate": "OPEN",
            "incoming_gate_campaign": "paused",
        },
        "costs": {},
        "npz_persisted_before_optional": False,
    }


def campaign(record, arrays):
    timings = record["costs"]
    before = source_hashes()
    record["hashes_before"] = before
    save(record)

    with Clock(timings, "nf128_solve"):
        loaded128 = galerkin.load_physical_columns(128)
        grid128 = galerkin.build_grid(128, quadrature=512)
        solved128, info128 = galerkin.solve_initial_radius(
            grid128, galerkin.blank_state(grid128, loaded128[0], loaded128[1])
        )
    record["nf128_seed"] = {
        "nf": 128,
        "nq": 512,
        "preparation": preparation_summary(loaded128),
        "converged": bool(info128.get("converged")),
        "residual_max": info128.get("residual_max"),
        "r_min_coarse": info128.get("r_min_coarse"),
        "r_max_coarse": info128.get("r_max_coarse"),
        "blocker": info128.get("blocker"),
    }
    save(record)
    if loaded128[3] or not info128.get("converged"):
        record["verdict"] = "BLOCKED_NF128_SEED"
        record["exact_failure"] = record["nf128_seed"]
        return record

    with Clock(timings, "nf256_grid_and_source"):
        loaded = galerkin.load_physical_columns(256)
        grid = galerkin.build_grid(256, quadrature=1024)
        state = galerkin.blank_state(grid, loaded[0], loaded[1])
        rho, rho_variation, rho_max = full_source_rho(grid, state)
    if grid.ng != 255 or grid.nf != 256 or grid.nq != 1024:
        raise RuntimeError(f"unexpected nf256 grid {(grid.ng, grid.nf, grid.nq)}")
    radius, embed_imag, high_mode = embed_radius(solved128.r, 127, 255, grid.length)
    seed_projected, seed_full = galerkin.projected_radius_operator(grid, radius, rho)
    with Clock(timings, "nf256_newton"):
        radius, best_projected, history = best_full_source_radius(grid, radius, rho)
    state.r = radius
    projected_now, full_now = galerkin.projected_radius_operator(grid, state.r, rho)
    record["preparation"] = preparation_summary(loaded)
    record["nf256_construction"] = {
        "nf": 256,
        "ng": 255,
        "nq": 1024,
        "state_label": "DIAGNOSTIC",
        "embed_imaginary_max": embed_imag,
        "high_mode_coefficient_max": high_mode,
        "nodal_r_min": float(np.min(radius)),
        "nodal_r_max": float(np.max(radius)),
        "rho_max": rho_max,
        "rho_independent_of_r_max": rho_variation,
        "seed_projected_operator_max": float(np.max(np.abs(seed_projected))),
        "seed_full_operator_max": float(np.max(np.abs(seed_full))),
        "best_projected_operator_max": float(best_projected),
        "best_full_operator_max": float(np.max(np.abs(full_now))),
        "newton_history": history,
        "newton_1e-10_passed": False,
        "projected_residual_leq_1e-10": bool(best_projected <= galerkin.TOL_NEWTON),
        "precision_limit_not_an_abort": True,
        "preparation_problems": list(loaded[3]),
    }
    record["occupations"] = occupation_report(state)
    arrays.update(state_arrays(
        grid,
        state,
        rho,
        prefix="",
        projected_residual=best_projected,
        full_operator=float(np.max(np.abs(full_now))),
    ))
    arrays["schema"] = np.array("NSC-SPHERICAL-COUPLING-REFINEMENT-v5")
    arrays["state_label"] = np.array("DIAGNOSTIC_FULL_SOURCE_BEST")
    with Clock(timings, "persist_nf256_npz"):
        write_npz(arrays)
    record["npz_persisted_before_optional"] = True
    record["npz_sha256_after_state"] = sha256(NPZ)
    save(record)
    if loaded[3]:
        record["verdict"] = "BLOCKED_UNCORRECTED_PREPARATION"
        record["exact_failure"] = list(loaded[3])
        return record

    with Clock(timings, "nf256_nq1024_constraints"):
        base = measure_constraints(grid, state)
    base["frozen_rho_projected_operator_max"] = float(np.max(np.abs(projected_now)))
    base["frozen_rho_full_operator_max"] = float(np.max(np.abs(full_now)))
    base["difference_from_v4_printed_full_C"] = float(base["full_C_hamilton_max"] - ARCHIVE_FULL_C)
    base["difference_from_v4_printed_projected_operator"] = float(
        best_projected - ARCHIVE_PROJECTED_OPERATOR
    )
    base["admissibility"] = admissibility(grid, state)
    with Clock(timings, "nf256_energy_and_work"):
        energy0, energy_error = safe_call(galerkin.energy, grid, state)
        balance, balance_error = safe_call(galerkin.work_balance, grid, state)
    base["energy"] = None if energy0 is None else float(energy0)
    base["energy_error"] = energy_error
    base["work_balance"] = balance
    base["work_balance_error"] = balance_error
    record["nq1024"] = base
    save(record)

    with Clock(timings, "nf256_nq2048_grid"):
        doubled = galerkin.build_grid(256, quadrature=2048)
    if doubled.ng != 255 or doubled.nq != 2048:
        raise RuntimeError(f"unexpected doubled grid {(doubled.ng, doubled.nf, doubled.nq)}")
    with Clock(timings, "nf256_nq2048_constraints"):
        record["nq2048"] = measure_constraints(doubled, state)
    record["nq2048"]["same_coarse_state_as_nq1024"] = True
    record["quadrature_comparison"] = {
        "full_C_nq1024": record["nq1024"]["full_C_hamilton_max"],
        "full_C_nq2048": record["nq2048"]["full_C_hamilton_max"],
        "full_C_2048_minus_1024": float(
            record["nq2048"]["full_C_hamilton_max"] - record["nq1024"]["full_C_hamilton_max"]
        ),
        "full_D_nq1024": record["nq1024"]["full_D_momentum_max"],
        "full_D_nq2048": record["nq2048"]["full_D_momentum_max"],
        "lifted_normal_velocity_nq1024": record["nq1024"]["lifted_normal_velocity_max"],
        "lifted_normal_velocity_nq2048": record["nq2048"]["lifted_normal_velocity_max"],
    }
    save(record)

    with Clock(timings, "nf256_evolve_dt"):
        record["window_dt_0_0005"] = run_window(grid, state, DT)
    save(record)
    with Clock(timings, "nf256_evolve_dt_half"):
        record["window_dt_0_00025"] = run_window(grid, state, DT_HALF)
    save(record)

    v2_payload = json.loads(V2.read_text())
    v2_full_c = float(v2_payload["runs"]["initial_nf128_nq512"]["diagnostics"]["full_hamilton_max"])
    full_c = record["nq1024"]["full_C_hamilton_max"]
    record["comparison_with_recorded_v2_initial_full_C"] = {
        "v2_nf128_full_C": v2_full_c,
        "nf256_nq1024_full_C": full_c,
        "ratio_256_over_128": float(full_c / v2_full_c) if v2_full_c else None,
    }
    print(
        "nf256 full C",
        full_c,
        "full D",
        record["nq1024"]["full_D_momentum_max"],
        "lifted",
        record["nq1024"]["lifted_normal_velocity_max"],
        "window C",
        record["window_dt_0_0005"].get("full_C_hamilton_max"),
        record["window_dt_0_00025"].get("full_C_hamilton_max"),
        flush=True,
    )
    record["verdict"] = "DIAGNOSTIC_MEASURED"
    save(record)
    maybe_nf512(record, arrays, radius)
    return record


def maybe_nf512(record, arrays, radius256):
    timings = record["costs"]
    newton = timings.get("nf256_newton", {}).get("cpu_seconds", math.inf)
    grid_2048 = timings.get("nf256_nq2048_grid", {}).get("cpu_seconds", math.inf)
    evolve = timings.get("nf256_evolve_dt", {}).get("cpu_seconds", math.inf)
    evolve_half = timings.get("nf256_evolve_dt_half", {}).get("cpu_seconds", math.inf)
    estimate = 1.25 * (grid_2048 + 8.0 * newton + 4.0 * (evolve + evolve_half))
    decision = {
        "budget_cpu_seconds": NF512_BUDGET_S,
        "estimated_additional_cpu_seconds": estimate,
        "safety_factor": 1.25,
        "scaling": "nq^3 on the Newton Jacobian from 1024 to 2048, nq^2 on the two short windows",
        "components_cpu_seconds": {
            "measured_nq2048_grid": grid_2048,
            "measured_nf256_newton": newton,
            "scaled_nf512_newton": 8.0 * newton,
            "measured_both_windows": evolve + evolve_half,
            "scaled_both_windows": 4.0 * (evolve + evolve_half),
        },
    }
    if not math.isfinite(estimate) or estimate >= NF512_BUDGET_S:
        decision["ran"] = False
        decision["reason"] = (
            f"estimated additional numerical CPU {estimate} is not under {NF512_BUDGET_S}"
        )
        record["nf512"] = decision
        save(record)
        return
    decision["ran"] = True
    decision["reason"] = "estimate was under 120 CPU seconds"
    record["nf512"] = decision
    started = time.process_time()

    def spent():
        return time.process_time() - started

    decision["stages_completed"] = []
    try:
        with Clock(timings, "nf512_grid_and_source"):
            loaded = galerkin.load_physical_columns(512)
            grid = galerkin.build_grid(512, quadrature=2048)
            blank = galerkin.blank_state(grid, loaded[0], loaded[1])
            rho, rho_variation, rho_max = full_source_rho(grid, blank)
        decision["stages_completed"].append("grid_and_source")
        if grid.nf != 512 or grid.ng != 511 or grid.nq != 2048:
            raise RuntimeError(f"unexpected nf512 grid {(grid.ng, grid.nf, grid.nq)}")
        radius, embed_imag, high_mode = embed_radius(radius256, 255, 511, grid.length)
        seed_projected, seed_full = galerkin.projected_radius_operator(grid, radius, rho)
        history = []
        best_residual = None
        best_radius = radius.copy()
        current = radius.copy()
        with Clock(timings, "nf512_newton"):
            for iteration in range(6):
                if spent() >= NF512_BUDGET_S and iteration:
                    history.append({"iteration": iteration, "stop": "cpu_budget"})
                    break
                projected, full = galerkin.projected_radius_operator(grid, current, rho)
                residual_max = float(np.max(np.abs(projected)))
                full_max = float(np.max(np.abs(full)))
                entry = {
                    "iteration": iteration,
                    "projected": residual_max,
                    "full_operator": full_max,
                }
                if best_residual is None or residual_max < best_residual:
                    best_residual = residual_max
                    best_radius = current.copy()
                    entry["retained_as_best"] = True
                else:
                    entry["retained_as_best"] = False
                if residual_max < galerkin.TOL_NEWTON:
                    entry["stop"] = "newton_tolerance"
                    history.append(entry)
                    break
                if spent() >= NF512_BUDGET_S:
                    entry["stop"] = "cpu_budget_before_jacobian"
                    history.append(entry)
                    break
                jacobian = galerkin._projected_radius_jacobian(grid, current)
                try:
                    delta = np.linalg.solve(jacobian, -projected)
                except np.linalg.LinAlgError as error:
                    entry["stop"] = "singular_jacobian"
                    entry["error"] = str(error)
                    history.append(entry)
                    break
                delta_max = float(np.max(np.abs(delta)))
                entry["delta_max"] = delta_max
                trial = current + delta
                if not np.isfinite(trial).all() or np.min(galerkin.prolong_geometry(grid, trial)) <= 0.0:
                    entry["stop"] = "left_positive_chart"
                    history.append(entry)
                    break
                trial_projected, _full_trial = galerkin.projected_radius_operator(grid, trial, rho)
                trial_max = float(np.max(np.abs(trial_projected)))
                entry["trial_projected"] = trial_max
                if trial_max < residual_max:
                    current = trial
                    entry["accepted"] = True
                    history.append(entry)
                else:
                    entry["accepted"] = False
                    entry["stop"] = "no_decrease"
                    history.append(entry)
                    break
                if delta_max < 1e-12:
                    break
        state = blank
        state.r = best_radius
        projected_now, full_now = galerkin.projected_radius_operator(grid, state.r, rho)
        arrays.update(state_arrays(
            grid,
            state,
            rho,
            prefix="nf512_",
            projected_residual=float(best_residual),
            full_operator=float(np.max(np.abs(full_now))),
        ))
        with Clock(timings, "persist_nf512_npz"):
            write_npz(arrays)
        decision["stages_completed"].append("persisted")
        decision["preparation"] = preparation_summary(loaded)
        decision["construction"] = {
            "nf": 512,
            "ng": 511,
            "nq": 2048,
            "state_label": "DIAGNOSTIC",
            "embed_imaginary_max": embed_imag,
            "high_mode_coefficient_max": high_mode,
            "nodal_r_min": float(np.min(best_radius)),
            "nodal_r_max": float(np.max(best_radius)),
            "rho_max": rho_max,
            "rho_independent_of_r_max": rho_variation,
            "seed_projected_operator_max": float(np.max(np.abs(seed_projected))),
            "seed_full_operator_max": float(np.max(np.abs(seed_full))),
            "best_projected_operator_max": float(best_residual),
            "best_full_operator_max": float(np.max(np.abs(full_now))),
            "newton_history": history,
            "newton_1e-10_passed": False,
            "projected_residual_leq_1e-10": bool(best_residual <= galerkin.TOL_NEWTON),
            "precision_limit_not_an_abort": True,
        }
        decision["occupations"] = occupation_report(state)
        save(record)
        with Clock(timings, "nf512_constraints"):
            decision["nq2048"] = measure_constraints(grid, state)
        decision["nq2048"]["admissibility"] = admissibility(grid, state)
        decision["stages_completed"].append("constraints")
        full_256 = record["nq1024"]["full_C_hamilton_max"]
        full_512 = decision["nq2048"]["full_C_hamilton_max"]
        decision["spatial_comparison"] = {
            "nf256_nq1024_full_C": full_256,
            "nf256_nq2048_full_C": record["nq2048"]["full_C_hamilton_max"],
            "nf512_nq2048_full_C": full_512,
            "ratio_512_over_256": float(full_512 / full_256) if full_256 else None,
            "nf256_lifted_normal_velocity": record["nq1024"]["lifted_normal_velocity_max"],
            "nf512_lifted_normal_velocity": decision["nq2048"]["lifted_normal_velocity_max"],
        }
        save(record)
        if spent() < NF512_BUDGET_S:
            with Clock(timings, "nf512_energy_and_work"):
                energy0, energy_error = safe_call(galerkin.energy, grid, state)
                balance, balance_error = safe_call(galerkin.work_balance, grid, state)
            decision["energy"] = None if energy0 is None else float(energy0)
            decision["energy_error"] = energy_error
            decision["work_balance"] = balance
            decision["work_balance_error"] = balance_error
            decision["stages_completed"].append("energy_and_work")
        else:
            decision["energy_and_work_skipped"] = "cpu budget reached after constraints"
        if spent() < NF512_BUDGET_S:
            with Clock(timings, "nf512_evolve_dt"):
                decision["window_dt_0_0005"] = run_window(grid, state, DT)
            decision["stages_completed"].append("window_dt")
            save(record)
        else:
            decision["window_dt_0_0005"] = {"skipped": "cpu budget reached"}
        if spent() < NF512_BUDGET_S:
            with Clock(timings, "nf512_evolve_dt_half"):
                decision["window_dt_0_00025"] = run_window(grid, state, DT_HALF)
            decision["stages_completed"].append("window_dt_half")
        else:
            decision["window_dt_0_00025"] = {"skipped": "cpu budget reached"}
        decision["additional_cpu_seconds"] = spent()
    except Exception as error:
        decision["failed"] = True
        decision["error"] = {
            "type": type(error).__name__,
            "message": str(error),
            "traceback": traceback.format_exc(),
        }
        decision["additional_cpu_seconds"] = spent()
        try:
            write_npz(arrays)
        except Exception as write_error:
            decision["npz_rewrite_error"] = str(write_error)
    record["nf512"] = decision
    save(record)


def seal(record, started, cpu, before):
    after = source_hashes()
    record["hashes_after"] = after
    watched = ("v1", "v2", "v3_driver", "v3_record", "v4_driver", "v4_record", "galerkin", "coupling", "action", "source")
    record["historical_v1_v4_unchanged"] = all(before[name] == after[name] for name in watched if name != "galerkin")
    record["module_unchanged_during_run"] = before["galerkin"] == after["galerkin"]
    record["v4_is_printed_output_archive_only"] = True
    record["runtime_seconds"] = time.perf_counter() - started
    record["cpu_seconds"] = time.process_time() - cpu
    record["renewal"] = False
    record["pass"] = False
    record["newton_1e-10_passed"] = False
    record["full_initial_1e-8_passed"] = False
    record["continuum_limit_claimed"] = False
    record["frozen_B_embedding_claimed"] = False
    record["time_extension_T_0_05"] = False
    return record


def main():
    started = time.perf_counter()
    cpu = time.process_time()
    record = blank_record()
    arrays = {}
    before = None
    try:
        before = source_hashes()
        record = campaign(record, arrays)
        record = seal(record, started, cpu, record.get("hashes_before") or before)
    except Exception as error:
        record["verdict"] = record.get("verdict") if record.get("verdict") not in (None, "RUNNING") else "RUN_FAILED"
        if record.get("nq1024") or record.get("npz_persisted_before_optional"):
            record["verdict"] = "DIAGNOSTIC_PARTIAL"
        record["run_error"] = {
            "type": type(error).__name__,
            "message": str(error),
            "traceback": traceback.format_exc(),
        }
        if arrays:
            try:
                write_npz(arrays)
                record["npz_persisted_before_optional"] = True
            except Exception as write_error:
                record["npz_rewrite_error"] = str(write_error)
        record = seal(record, started, cpu, record.get("hashes_before") or before or source_hashes())
        save(record)
        print(json.dumps(jsonable({
            "verdict": record.get("verdict"),
            "run_error": record.get("run_error", {}).get("message"),
            "npz": str(NPZ),
            "json": str(OUT),
        }), indent=2), flush=True)
        return 1
    save(record)
    nf512 = record.get("nf512") or {}
    summary = {
        "verdict": record.get("verdict"),
        "pass": False,
        "newton_1e-10_passed": False,
        "full_initial_1e-8_passed": False,
        "renewal": False,
        "nf256_full_C_nq1024": (record.get("nq1024") or {}).get("full_C_hamilton_max"),
        "nf256_full_D_nq1024": (record.get("nq1024") or {}).get("full_D_momentum_max"),
        "nf256_lifted_nq1024": (record.get("nq1024") or {}).get("lifted_normal_velocity_max"),
        "nf256_full_C_nq2048": (record.get("nq2048") or {}).get("full_C_hamilton_max"),
        "nf256_lifted_nq2048": (record.get("nq2048") or {}).get("lifted_normal_velocity_max"),
        "window_dt_full_C": (record.get("window_dt_0_0005") or {}).get("full_C_hamilton_max"),
        "window_dt_full_D": (record.get("window_dt_0_0005") or {}).get("full_D_momentum_max"),
        "window_dt_energy_drift": (record.get("window_dt_0_0005") or {}).get("energy_drift"),
        "window_dt_work": (record.get("window_dt_0_0005") or {}).get("work_fieldwork_integral"),
        "window_dt_admissible": ((record.get("window_dt_0_0005") or {}).get("admissibility") or {}).get("admissible"),
        "window_half_full_C": (record.get("window_dt_0_00025") or {}).get("full_C_hamilton_max"),
        "window_half_full_D": (record.get("window_dt_0_00025") or {}).get("full_D_momentum_max"),
        "window_half_energy_drift": (record.get("window_dt_0_00025") or {}).get("energy_drift"),
        "window_half_work": (record.get("window_dt_0_00025") or {}).get("work_fieldwork_integral"),
        "window_half_admissible": ((record.get("window_dt_0_00025") or {}).get("admissibility") or {}).get("admissible"),
        "quadrature_comparison": record.get("quadrature_comparison"),
        "nf512_ran": nf512.get("ran"),
        "nf512_reason": nf512.get("reason"),
        "nf512_full_C": (nf512.get("nq2048") or {}).get("full_C_hamilton_max"),
        "nf512_spatial": nf512.get("spatial_comparison"),
        "npz_persisted_before_optional": record.get("npz_persisted_before_optional"),
        "cpu_seconds": record.get("cpu_seconds"),
        "runtime_seconds": record.get("runtime_seconds"),
        "json_sha256": sha256(OUT),
        "npz_sha256": sha256(NPZ),
    }
    print(json.dumps(jsonable(summary), indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
