#!/usr/bin/env python3
"""One replayable episode of the autonomous spherical Galerkin model.

The production evolution is this driver. It loads the saved self-contained v5
preparation, calls the existing RK4 step, and writes one JSON record and one
NPZ payload. It does not retune floors after the run, project out a mean, or
add a reset, a prescribed trajectory, or a restoring force.

Each RK4 stage evaluates the current source, geometry, and Dirac image through
``nsc_spherical_galerkin_coupling.rk4_step``. The historical initial 1e-8 line
is recorded and does not stop the episode. T* D = D is reused and is not
recomputed. The incoming gate is not a prerequisite.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import sys
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

from recursive_horizons import nsc_regional_energy_exchange as regional
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_spherical_coupling import (
    CauchyRate,
    CauchyState,
    PositiveChartExit,
    field_energy,
    gravity_energy,
)
from recursive_horizons.nsc_spherical_feedback_action import weyl_square_from_rh

LAB = Path(__file__).resolve().parents[1]
ROOT = LAB.parent
OUT = LAB / "results" / "development" / "nsc-spherical-feedback-episode-v1.json"
NPZ = LAB / "results" / "development" / "nsc-spherical-feedback-episode-v1.npz"
DRIVER = Path(__file__).resolve()
V5_JSON = LAB / "results" / "development" / "nsc-spherical-coupling-refinement-v5.json"
V5_NPZ = LAB / "results" / "development" / "nsc-spherical-coupling-refinement-v5.npz"
REGIONAL_JSON = LAB / "results" / "development" / "nsc-regional-energy-exchange-v1.json"

SCHEMA = "NSC-SPHERICAL-FEEDBACK-EPISODE-v1"
DURATION = 0.05
DT = 0.0005
DT_HALF = 0.00025
CPU_BUDGET_S = 300.0
FRAME = 0.005
PAYLOAD_LIMIT = 64 * 1024 * 1024
ONE_PERCENT = 0.01
BUDGET_SAFETY = 1.15

WINDOW_CENTERS = (1.0, 3.0, 5.0, 7.0)
WINDOW_HALF_WIDTH = 2.0
WINDOW_NAMES = ("packet_0_2", "packet_2_4", "packet_4_6", "complement_6_8")

INTENDED_RUNS = (
    {"name": "nf512_dt_0_0005", "resolution": "nf512", "dt": DT, "role": "primary"},
    {"name": "nf512_dt_0_00025", "resolution": "nf512", "dt": DT_HALF, "role": "time"},
    {"name": "nf256_dt_0_0005", "resolution": "nf256", "dt": DT, "role": "space"},
    {"name": "nf256_dt_0_00025", "resolution": "nf256", "dt": DT_HALF, "role": "supporting_time"},
)

# Absolute floors declared before the episode. A smaller change is not claimed
# as an effect. Floors are not pass lines and are not retuned after the run.
PHYSICAL_EFFECTS = (
    ("proper_velocity_max", "proper_max", 1e-8),
    ("proper_velocity_min", "proper_min", 1e-8),
    ("proper_velocity_mean", "proper_mean", 1e-8),
    ("observer_K_perp_max", "K_perp_max", 1e-8),
    ("observer_K_perp_min", "K_perp_min", 1e-8),
    ("observer_K_r_max", "K_r_max", 1e-8),
    ("observer_K_r_min", "K_r_min", 1e-8),
    ("areal_radius_min", "r_min", 1e-8),
    ("areal_radius_max", "r_max", 1e-8),
    ("conformal_Q_min", "Q_min", 1e-8),
    ("conformal_Q_max", "Q_max", 1e-8),
    ("chi_min", "chi_min", 1e-8),
    ("chi_max", "chi_max", 1e-8),
    ("weyl_C2_max", "c2_max", 1e-12),
    ("normal_shell_packet_0_2", "window_normal_0", 1e-6),
    ("normal_shell_packet_2_4", "window_normal_1", 1e-6),
    ("normal_shell_packet_4_6", "window_normal_2", 1e-6),
    ("normal_shell_complement_6_8", "window_normal_3", 1e-6),
    ("proper_boundary_flux_packet_0_2", "window_proper_flux_0", 1e-6),
    ("proper_boundary_flux_packet_2_4", "window_proper_flux_1", 1e-6),
    ("proper_boundary_flux_packet_4_6", "window_proper_flux_2", 1e-6),
    ("proper_boundary_flux_complement_6_8", "window_proper_flux_3", 1e-6),
    ("quasilocal_boundary_packet_0_2", "window_quasilocal_0", 1e-6),
    ("quasilocal_boundary_packet_2_4", "window_quasilocal_1", 1e-6),
    ("quasilocal_boundary_packet_4_6", "window_quasilocal_2", 1e-6),
    ("quasilocal_boundary_complement_6_8", "window_quasilocal_3", 1e-6),
    ("coordinate_work_integral", "coordinate_work_integral", 1e-6),
    ("proper_pressure_work_integral", "proper_pressure_work_integral", 1e-6),
    ("shared_mode_radius_l2", "mode_r", 1e-8),
    ("shared_mode_Q_l2", "mode_Q", 1e-8),
    ("shared_mode_chi_l2", "mode_chi", 1e-8),
)

DIAGNOSTIC_EFFECTS = (
    ("full_quadrature_hamilton", "full_hamilton_max", 1e-8),
    ("geometry_band_projected_hamilton", "projected_hamilton_max", 1e-10),
    ("held_out_complement_hamilton", "held_out_hamilton_max", 1e-8),
    ("full_quadrature_momentum", "full_momentum_max", 1e-10),
    ("geometry_band_projected_momentum", "projected_momentum_max", 1e-12),
    ("held_out_complement_momentum", "held_out_momentum_max", 1e-10),
    ("coordinate_energy", "energy", 1e-8),
)

SERIES_KEYS = (
    "time",
    "energy",
    "field_energy",
    "gravity_energy",
    "full_hamilton_max",
    "projected_hamilton_max",
    "held_out_hamilton_max",
    "full_momentum_max",
    "projected_momentum_max",
    "held_out_momentum_max",
    "mean_rho",
    "mean_current",
    "mean_K",
    "mean_Pmom",
    "proper_min",
    "proper_max",
    "proper_mean",
    "lifted_proper_max",
    "unprojected_proper_max",
    "chart_proper_max",
    "coordinate_over_N_max",
    "shift_piece_max",
    "r_dot_max",
    "K_r_min",
    "K_r_max",
    "K_perp_min",
    "K_perp_max",
    "chi_min",
    "chi_max",
    "c2_max",
    "r_min",
    "r_max",
    "Q_min",
    "Q_max",
    "lifted_fieldwork",
    "proper_pressure_power",
    "momentum_lapse_power",
    "coordinate_work_integral",
    "proper_pressure_work_integral",
    "unitarity",
    "number",
    "ward_max",
    "shell_gap_max",
    "source_gap_max",
    "partition_gap",
)

WINDOW_SERIES = (
    "window_normal",
    "window_matter",
    "window_gravity",
    "window_coordinate",
    "window_shift_transport",
    "window_proper_flux",
    "window_shift_pressure",
    "window_coordinate_work",
    "window_geometric_shift",
    "window_geometric_proper",
    "window_proper_work",
    "window_lapse_work",
    "window_quasilocal",
)

WITNESS_KEYS = (
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
    "r_min",
    "Q_min",
)

FROZEN_PATHS = {
    "v5_json": V5_JSON,
    "v5_npz": V5_NPZ,
    "regional_json": REGIONAL_JSON,
    "v1_json": LAB / "results" / "development" / "nsc-spherical-coupling-control-v1.json",
    "v2_json": LAB / "results" / "development" / "nsc-spherical-coupling-control-v2.json",
    "v3_json": LAB / "results" / "development" / "nsc-spherical-coupling-refinement-v3.json",
    "galerkin": galerkin._MODULE_PATH,
    "regional_module": regional._MODULE_PATH,
    "feedback_action": ROOT / "lab" / "src" / "recursive_horizons" / "nsc_spherical_feedback_action.py",
    "coupling": galerkin._COUPLING_PATH,
    "conformal_source": galerkin._SOURCE_PATH,
    "manuscript_md": ROOT / "paper" / "nested-space-cosmology.md",
    "manuscript_pdf": ROOT / "paper" / "nested-space-cosmology.pdf",
    "companion_pdf": ROOT / "paper" / "local-incoming-gate-draft.pdf",
}


def limit_threads():
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["NUMEXPR_NUM_THREADS"] = "1"
    os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
    try:
        from threadpoolctl import threadpool_limits as limits

        limits(limits=1).__enter__()
        return 1
    except Exception as error:
        return "not applied: " + type(error).__name__


def sha256(path):
    file_path = Path(path)
    if not file_path.is_file():
        return None
    digest = hashlib.sha256()
    with file_path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def frozen_hashes():
    return {name: sha256(path) for name, path in FROZEN_PATHS.items()}


def jsonable(value):
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": jsonable(np.real(value)), "imag": jsonable(np.imag(value))}
        return jsonable(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, complex):
        return {"real": jsonable(value.real), "imag": jsonable(value.imag)}
    if value is None or isinstance(value, str):
        return value
    return str(value)


def write_json(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(jsonable(record), indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def write_npz(path, arrays):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    os.replace(temporary, path)


def budget_allows(elapsed, step_cost, budget=CPU_BUDGET_S, safety=BUDGET_SAFETY):
    """True when another measured step still fits in the remaining CPU budget."""
    elapsed = float(elapsed)
    step_cost = float(step_cost)
    budget = float(budget)
    if elapsed >= budget:
        return False
    if step_cost <= 0.0:
        return True
    return elapsed + safety * step_cost < budget


def classify_movement(effect_scale, movement, floor):
    """Compare one refinement movement with one percent of the episode change.

    ``below_noise_floor`` is not a resolved physical effect. ``unresolved`` is
    kept when the movement exceeds one percent of a real effect, or when the
    runs disagree above the floor while the primary change does not.
    """
    effect_scale = abs(float(effect_scale))
    movement = abs(float(movement))
    floor = abs(float(floor))
    if effect_scale < floor and movement < floor:
        return "below_noise_floor"
    if effect_scale >= floor and movement <= ONE_PERCENT * effect_scale:
        return "resolved"
    return "unresolved"


def summarize_statuses(statuses):
    statuses = list(statuses)
    if not statuses or any(status == "incomplete" for status in statuses):
        return "incomplete"
    if any(status == "unresolved" for status in statuses):
        return "unresolved"
    if any(status == "resolved" for status in statuses):
        return "resolved"
    if all(status == "below_noise_floor" for status in statuses):
        return "no_measurable_evolution"
    return "incomplete"


def balance_against_exchange(balance_change, exchange):
    """Compare a near-zero balance drift with the observed energy exchange.

    One percent of the drift itself is not the scale. The exchange is the
    field-energy change, whose size is the integrated coordinate work.
    """
    balance = abs(float(balance_change))
    scale = abs(float(exchange))
    ratio = None if scale == 0.0 else balance / scale
    return {
        "balance_change": float(balance_change),
        "exchange": float(exchange),
        "balance_over_exchange": ratio,
        "within_one_percent_of_exchange": bool(ratio is not None and ratio <= ONE_PERCENT),
        "comparison": "balance error against observed energy exchange",
        "not_one_percent_of_the_balance_error": True,
    }


def common_comparison_time(time_a, time_b, dt):
    if time_a is None or time_b is None or len(time_a) == 0 or len(time_b) == 0:
        return None
    limit = min(float(time_a[-1]), float(time_b[-1]))
    steps = int(math.floor(limit / float(dt) + 1e-9))
    if steps <= 0:
        return None
    return float(steps * float(dt))


def value_at(times, values, time):
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    if times.size == 0:
        raise KeyError("empty series")
    index = int(np.argmin(np.abs(times - float(time))))
    if abs(float(times[index]) - float(time)) > 1e-9:
        raise KeyError(time)
    return float(values[index])


def initial_tolerance_decision(full_hamilton_max, full_momentum_max, threshold=None):
    threshold = float(galerkin.TOL_INITIAL_CONSTRAINT if threshold is None else threshold)
    passed = bool(float(full_hamilton_max) <= threshold and float(full_momentum_max) <= threshold)
    return {
        "threshold": threshold,
        "full_hamilton_max": float(full_hamilton_max),
        "full_momentum_max": float(full_momentum_max),
        "historical_initial_1e-8_passed": passed,
        "role": "diagnostic",
        "abort_evolution": False,
        "continue_evolution": True,
    }


def _prefix(label):
    return "" if label in ("", "nf256") else label + "_"


def load_v5_state(label):
    """Load one saved v5 state. Does not solve a radius and does not write."""
    prefix = _prefix(label)
    with np.load(V5_NPZ, allow_pickle=False) as data:
        fermions = int(data[prefix + "nf"])
        quadrature = int(data[prefix + "nq"])
        length = float(data[prefix + "length"])
        state = CauchyState(
            Q=np.array(data[prefix + "geometry_Q"], dtype=float, copy=True),
            r=np.array(data[prefix + "geometry_r"], dtype=float, copy=True),
            chi=np.array(data[prefix + "geometry_chi"], dtype=float, copy=True),
            p_Q=np.array(data[prefix + "geometry_p_Q"], dtype=float, copy=True),
            p_r=np.array(data[prefix + "geometry_p_r"], dtype=float, copy=True),
            p_chi=np.array(data[prefix + "geometry_p_chi"], dtype=float, copy=True),
            phi0=np.array(data[prefix + "columns_phi0"], copy=True),
            phi1=np.array(data[prefix + "columns_phi1"], copy=True),
        )
        recorded_radius = np.array(data[prefix + "quadrature_radius"], dtype=float, copy=True)
        occupations = np.array(data[prefix + "occupations"], dtype=float, copy=True)
        binding_names = np.array(data[prefix + "binding_names"])
        binding_values = np.array(data[prefix + "binding_values"], dtype=float, copy=True)
        saved_rho = np.array(data[prefix + "rho"], dtype=float, copy=True)
        geometry_modes = np.array(data[prefix + "geometry_modes"], dtype=float, copy=True)
    grid = galerkin.build_grid(fermions, quadrature=quadrature, length=length)
    if state.r.shape != (grid.ng,) or state.phi0.shape != (grid.nf, 6):
        raise ValueError(f"saved {label} arrays do not match the owned grid")
    if not np.array_equal(geometry_modes, grid.modes_g):
        raise ValueError(f"saved {label} geometry modes differ from the owned grid")
    metadata = {
        "label": "nf256" if prefix == "" else "nf512",
        "nf": int(grid.nf),
        "ng": int(grid.ng),
        "nq": int(grid.nq),
        "length": float(grid.length),
        "recorded_radius": recorded_radius,
        "occupations": occupations,
        "binding_names": binding_names,
        "binding_values": binding_values,
        "saved_rho": saved_rho,
        "v5_sha256": sha256(V5_NPZ),
        "occupation_gap": float(np.max(np.abs(occupations - grid.fine.occupations))),
        "quadrature_radius_gap": float(np.max(np.abs(
            galerkin.prolong_geometry(grid, state.r) - recorded_radius
        ))),
    }
    return grid, state, metadata


def windows_for(grid):
    _coordinate, windows, partition = regional.fixed_window(
        grid.length, grid.nq, WINDOW_CENTERS, WINDOW_HALF_WIDTH
    )
    return windows, partition


def _smear(weights, values):
    return float(np.sum(np.asarray(weights) * np.asarray(values)))


def lifted_rate(grid, coarse, bundle):
    source = bundle["source"]
    return CauchyRate(
        Q=galerkin.prolong_geometry(grid, coarse.Q),
        r=galerkin.prolong_geometry(grid, coarse.r),
        chi=galerkin.prolong_geometry(grid, coarse.chi),
        p_Q=galerkin.prolong_geometry(grid, coarse.p_Q),
        p_r=galerkin.prolong_geometry(grid, coarse.p_r),
        p_chi=galerkin.prolong_geometry(grid, coarse.p_chi),
        phi0=galerkin.prolong_columns(grid, coarse.phi0),
        phi1=galerkin.prolong_columns(grid, coarse.phi1),
        force_L=source["force_L"],
        force_Q=source["force_Q"],
        force_beta=source["force_beta"],
        fieldwork_power=float(coarse.fieldwork_power),
    )


def signed_proper_motion(grid, fine, lifted, unprojected_r):
    """Split coordinate radial velocity into shift drag and proper motion.

    Lifted proper velocity is (r_dot_lifted - beta r_x) / N, the same normal
    change recorded by the v5 diagnostic. Observer K_perp divides by an extra
    areal radius. Neither one is F_Q Q_dot.
    """
    system = grid.fine
    radius_x = system.derivative @ fine.r
    denominator = fine.r * system.length_density
    proper = (lifted.r - system.shift * radius_x) / denominator
    coordinate = lifted.r / denominator
    shift_piece = system.shift * radius_x / denominator
    unprojected = (unprojected_r - system.shift * radius_x) / denominator
    return {
        "proper": proper,
        "coordinate_over_N": coordinate,
        "shift_piece": shift_piece,
        "unprojected_proper": unprojected,
        "r_dot": lifted.r,
        "identity_gap": float(np.max(np.abs(proper + shift_piece - coordinate))),
    }


def curvature_fields(fine):
    """chi is the chart variable with the owned identity chi = R_h - 2.

    The four-dimensional scalar used here is the owned product
    C^2 = (R_h - 2)^2 / (3 r^4). The 4D Riemann tensor is not recomputed.
    """
    chi = np.asarray(fine.chi, dtype=float)
    weyl = np.asarray(weyl_square_from_rh(chi + 2.0, fine.r), dtype=float)
    identity = chi ** 2 / (3.0 * np.asarray(fine.r, dtype=float) ** 4)
    return {
        "chi": chi,
        "weyl_C2": weyl,
        "identity_gap": float(np.max(np.abs(weyl - identity))),
        "chi_equals_Rh_minus_2": True,
        "owner": "nsc_spherical_feedback_action.weyl_square_from_rh",
        "riemann_tensor_recomputed": False,
    }


def window_row(system, window, ledger, geometry, terms):
    derivative = system.derivative @ window
    work = regional.coordinate_metric_work(ledger, terms["rate"])
    columns = []
    for column in ledger["columns"]:
        columns.append({
            "occupation": float(column["occupation"]),
            "shift_transport": _smear(derivative, column["shift"]),
            "proper_normal_flux": _smear(derivative, column["proper"]),
            "shift_pressure_cross": _smear(derivative, column["shift_pressure"]),
            "lapse_weighted_energy": _smear(window, column["lapse_weighted"]),
            "hamiltonian_energy": _smear(window, column["hamiltonian"]),
        })
    return {
        "normal_energy": _smear(window, ledger["normal_energy_nodal"]),
        "lapse_weighted_energy": _smear(window, ledger["lapse_weighted_integrand"]),
        "matter_energy": _smear(window, ledger["hamiltonian_integrand"]),
        "gravity_energy": _smear(window, geometry["nodal_energy"]),
        "measured_energy": _smear(window, geometry["nodal_energy"] + ledger["hamiltonian_integrand"]),
        "shift_transport": _smear(derivative, ledger["flux_shift"]),
        "proper_normal_flux": _smear(derivative, ledger["flux_proper"]),
        "shift_pressure_cross": _smear(derivative, ledger["flux_shift_pressure"]),
        "coordinate_metric_work": _smear(window, work),
        "geometric_shift_flux": _smear(derivative, geometry["nodal_shift_flux"]),
        "geometric_proper_flux": _smear(derivative, geometry["nodal_proper_flux"]),
        "proper_pressure_work": _smear(window, terms["proper_pressure_work"]),
        "momentum_lapse_work": _smear(window, terms["momentum_lapse_work"]),
        "quasilocal_boundary": -_smear(derivative, system.dx * geometry["quasilocal_density"]),
        "mode_columns": columns,
    }


def observe(grid, state, windows=None, partition=None):
    """One fresh source, geometry, and Hamiltonian evaluation of this state."""
    coarse, bundle = galerkin.compose_fine_hamiltonian(grid, state, include_matter_force=True)
    fine = bundle["fine_state"]
    source = bundle["source"]
    lifted = lifted_rate(grid, coarse, bundle)
    system = grid.fine
    diagnostics = galerkin.constraint_diagnostics(grid, state, source, fine)
    velocities = galerkin.normal_velocities(grid, state, coarse, bundle)
    motion = signed_proper_motion(grid, fine, lifted, bundle["unprojected_r"])
    curvature = curvature_fields(fine)
    ledger = regional.matter_ledger(system, fine)
    terms = regional.proper_balance_terms(system, fine, lifted, ledger)
    terms = dict(terms)
    terms["rate"] = lifted
    geometry = regional.geometric_ledger(system, fine, lifted)
    if windows is None or partition is None:
        windows, partition = windows_for(grid)
    rows = [window_row(system, window, ledger, geometry, terms) for window in windows]
    for row, name in zip(rows, WINDOW_NAMES, strict=True):
        row["name"] = name
    rho = source["force_L"] / grid.dx_q
    current = source["force_beta"] / grid.dx_q
    images = (source["image0"], source["image1"])
    field = float(field_energy(system, fine, images))
    gravity = float(gravity_energy(system, fine))
    gram = fine.phi0.conj().T @ fine.phi0 + fine.phi1.conj().T @ fine.phi1
    number = float(np.sum(
        system.occupations * np.sum(np.abs(fine.phi0) ** 2 + np.abs(fine.phi1) ** 2, axis=0)
    ))
    kernel_gap = max(float(value) for value in ledger["kernel_gap"].values())
    source_gap = float(np.max(np.abs(ledger["source"]["force_L"] - source["force_L"])))
    sample = {
        "time": 0.0,
        "coordinate_work_integral": 0.0,
        "proper_pressure_work_integral": 0.0,
        "energy": field + gravity,
        "field_energy": field,
        "gravity_energy": gravity,
        "full_hamilton_max": float(diagnostics["full_hamilton_max"]),
        "projected_hamilton_max": float(diagnostics["projected_hamilton_max"]),
        "held_out_hamilton_max": float(diagnostics["held_out_hamilton_max"]),
        "full_momentum_max": float(diagnostics["full_momentum_max"]),
        "projected_momentum_max": float(diagnostics["projected_momentum_max"]),
        "held_out_momentum_max": float(diagnostics["held_out_momentum_max"]),
        "mean_rho": float(np.mean(rho)),
        "mean_current": float(np.mean(current)),
        "mean_K": float(np.mean(source["K"])),
        "mean_Pmom": float(np.mean(source["Pmom"])),
        "mean_subtracted": False,
        "proper_min": float(np.min(motion["proper"])),
        "proper_max": float(np.max(motion["proper"])),
        "proper_mean": float(np.mean(motion["proper"])),
        "lifted_proper_max": float(velocities["lifted_proper_max"]),
        "unprojected_proper_max": float(np.max(np.abs(motion["unprojected_proper"]))),
        "chart_proper_max": float(velocities["chart_proper_max"]),
        "coordinate_over_N_max": float(np.max(np.abs(motion["coordinate_over_N"]))),
        "shift_piece_max": float(np.max(np.abs(motion["shift_piece"]))),
        "r_dot_max": float(np.max(np.abs(motion["r_dot"]))),
        "proper_shift_identity_gap": motion["identity_gap"],
        "signed_versus_owner_gap": float(np.max(np.abs(motion["proper"])) - velocities["lifted_proper_max"]),
        "K_r_min": float(np.min(terms["K_r"])),
        "K_r_max": float(np.max(terms["K_r"])),
        "K_perp_min": float(np.min(terms["K_perp"])),
        "K_perp_max": float(np.max(terms["K_perp"])),
        "chi_min": float(np.min(curvature["chi"])),
        "chi_max": float(np.max(curvature["chi"])),
        "c2_max": float(np.max(np.abs(curvature["weyl_C2"]))),
        "c2_identity_gap": curvature["identity_gap"],
        "r_min": float(diagnostics["r_min_quadrature"]),
        "r_max": float(np.max(fine.r)),
        "Q_min": float(diagnostics["Q_min_quadrature"]),
        "Q_max": float(np.max(fine.Q)),
        "lifted_fieldwork": float(coarse.fieldwork_power),
        "proper_pressure_power": float(np.sum(terms["proper_pressure_work"])),
        "momentum_lapse_power": float(np.sum(terms["momentum_lapse_work"])),
        "coordinate_work_l1": float(np.sum(np.abs(regional.coordinate_metric_work(ledger, lifted)))),
        "proper_pressure_work_l1": float(np.sum(np.abs(terms["proper_pressure_work"]))),
        "unitarity": float(np.max(np.abs(gram - np.eye(6)))),
        "number": number,
        "ward_max": float(np.max(np.abs(ledger["hamiltonian_chain_ward"]))),
        "shell_gap_max": float(np.max(np.abs(ledger["stresses"]["shell_gap"]))),
        "kernel_gap_max": kernel_gap,
        "source_gap_max": source_gap,
        "partition_gap": float(np.max(np.abs(partition - 1.0))),
        "positive_r": bool(diagnostics["positive_r"]),
        "positive_Q": bool(diagnostics["positive_Q"]),
        "mode_projectors": {
            "full_quadrature_hamilton_max": float(diagnostics["full_hamilton_max"]),
            "geometry_band_projected_hamilton_max": float(diagnostics["projected_hamilton_max"]),
            "held_out_complement_hamilton_max": float(diagnostics["held_out_hamilton_max"]),
            "full_quadrature_momentum_max": float(diagnostics["full_momentum_max"]),
            "geometry_band_projected_momentum_max": float(diagnostics["projected_momentum_max"]),
            "held_out_complement_momentum_max": float(diagnostics["held_out_momentum_max"]),
            "lifted_proper_velocity": "subspace pullback of r_dot, then (r_dot - beta r_x) / N",
            "unprojected_proper_velocity": "fine geometric r_dot before the geometry pullback",
            "chart_proper_velocity": "momenta expression from the owned chart, not the integrated r_dot",
        },
        "windows": rows,
        "window_normal": [row["normal_energy"] for row in rows],
        "window_matter": [row["matter_energy"] for row in rows],
        "window_gravity": [row["gravity_energy"] for row in rows],
        "window_coordinate": [row["measured_energy"] for row in rows],
        "window_shift_transport": [row["shift_transport"] for row in rows],
        "window_proper_flux": [row["proper_normal_flux"] for row in rows],
        "window_shift_pressure": [row["shift_pressure_cross"] for row in rows],
        "window_coordinate_work": [row["coordinate_metric_work"] for row in rows],
        "window_geometric_shift": [row["geometric_shift_flux"] for row in rows],
        "window_geometric_proper": [row["geometric_proper_flux"] for row in rows],
        "window_proper_work": [row["proper_pressure_work"] for row in rows],
        "window_lapse_work": [row["momentum_lapse_work"] for row in rows],
        "window_quasilocal": [row["quasilocal_boundary"] for row in rows],
        "rho_mean_removed_before_constraint": False,
        "frozen_source_used": False,
    }
    for index in range(4):
        sample[f"window_normal_{index}"] = sample["window_normal"][index]
        sample[f"window_proper_flux_{index}"] = sample["window_proper_flux"][index]
        sample[f"window_quasilocal_{index}"] = sample["window_quasilocal"][index]
    fields = {
        "r": np.asarray(fine.r, dtype=float),
        "Q": np.asarray(fine.Q, dtype=float),
        "chi": np.asarray(curvature["chi"], dtype=float),
        "proper": np.asarray(motion["proper"], dtype=float),
        "K_r": np.asarray(terms["K_r"], dtype=float),
        "K_perp": np.asarray(terms["K_perp"], dtype=float),
        "rho": np.asarray(rho, dtype=float),
        "current": np.asarray(current, dtype=float),
        "weyl_C2": np.asarray(curvature["weyl_C2"], dtype=float),
    }
    return sample, coarse, fields


def modal_coefficients(grid, nodal):
    modes = grid.modes_g
    analysis = np.exp(-2j * np.pi * grid.xi_g[:, None] * modes[None, :] / grid.length) / grid.ng
    return analysis.T @ np.asarray(nodal, dtype=float)


def _frame_time(time, final):
    if final or time <= 1e-15:
        return True
    nearest = round(float(time) / FRAME) * FRAME
    return abs(float(time) - nearest) <= 1e-10


def evolve_episode(grid, state, duration, dt, elapsed, budget=CPU_BUDGET_S, last_step_cost=0.0, label="episode"):
    """Integrate with the owned RK4 step until the duration or the CPU budget."""
    steps = int(np.ceil(float(duration) / float(dt) - 1e-12))
    dt_used = float(duration) / steps
    windows, partition = windows_for(grid)
    current = state.copy()
    samples = []
    frames = []
    modes = {"time": [], "r": [], "Q": [], "chi": []}
    coordinate_integral = 0.0
    proper_integral = 0.0
    previous_power = None
    previous_proper = None
    stop_reason = None
    chart_stop = False
    completed_steps = 0
    step_cost = float(last_step_cost)
    state_after_one_step = None
    rate_after_one_step = None
    last_rate = None
    cpu_start = time.process_time()
    try:
        sample, last_rate, fields = observe(grid, current, windows, partition)
        sample["time"] = 0.0
        sample["coordinate_work_integral"] = 0.0
        sample["proper_pressure_work_integral"] = 0.0
        samples.append(sample)
        previous_power = sample["lifted_fieldwork"]
        previous_proper = sample["proper_pressure_power"]
        _store_frame(grid, current, sample, fields, last_rate, frames, modes, final=False)
        for step in range(steps):
            now = float(elapsed())
            if not budget_allows(now, step_cost if step else 0.0, budget):
                stop_reason = "CPU_BUDGET_300S"
                break
            step_started = time.process_time()
            try:
                candidate = galerkin.rk4_step(grid, current, dt_used)
                sample, rate, fields = observe(grid, candidate, windows, partition)
            except PositiveChartExit as exit_chart:
                stop_reason = str(exit_chart.reason)
                chart_stop = True
                break
            current = candidate
            last_rate = rate
            step_time = (step + 1) * dt_used
            coordinate_integral += 0.5 * dt_used * (previous_power + sample["lifted_fieldwork"])
            proper_integral += 0.5 * dt_used * (previous_proper + sample["proper_pressure_power"])
            previous_power = sample["lifted_fieldwork"]
            previous_proper = sample["proper_pressure_power"]
            sample["time"] = float(step_time)
            sample["coordinate_work_integral"] = float(coordinate_integral)
            sample["proper_pressure_work_integral"] = float(proper_integral)
            samples.append(sample)
            final = step + 1 == steps
            if _frame_time(step_time, final):
                _store_frame(grid, current, sample, fields, rate, frames, modes, final=final)
            if step == 0:
                state_after_one_step = current.copy()
                rate_after_one_step = rate
            completed_steps = step + 1
            step_cost = time.process_time() - step_started
            if step == 0 or completed_steps % 10 == 0 or final:
                print(
                    label,
                    "step",
                    completed_steps,
                    "T",
                    sample["time"],
                    "cpu",
                    round(time.process_time() - cpu_start, 3),
                    "proper_max",
                    sample["proper_max"],
                    "full_H",
                    sample["full_hamilton_max"],
                    flush=True,
                )
    except Exception as error:
        if not samples:
            raise
        stop_reason = "CODE_EXCEPTION"
        error_text = traceback.format_exc()
        print("exception", label, error, flush=True)
    else:
        error_text = None
    attained = 0.0 if not samples else float(samples[-1]["time"])
    return {
        "samples": samples,
        "state": current,
        "state_after_one_step": state_after_one_step,
        "rate_after_one_step": rate_after_one_step,
        "last_rate": last_rate,
        "frames": frames,
        "modes": modes,
        "attained_time": attained,
        "requested_time": float(duration),
        "dt": dt_used,
        "steps_requested": steps,
        "steps_completed": completed_steps,
        "completed": bool(completed_steps == steps and stop_reason is None),
        "stop_reason": stop_reason,
        "chart_stop": chart_stop,
        "cpu_seconds": float(time.process_time() - cpu_start),
        "step_cost": float(step_cost),
        "last_nq": int(grid.nq),
        "error": error_text,
        "positive_r": all(bool(sample["positive_r"]) for sample in samples),
        "positive_Q": all(bool(sample["positive_Q"]) for sample in samples),
    }


def _store_frame(grid, state, sample, fields, rate, frames, modes, final):
    frames.append({
        "time": float(sample["time"]),
        "final": bool(final),
        "fields": fields,
        "state": state.copy(),
        "rate": rate,
    })
    modes["time"].append(float(sample["time"]))
    modes["r"].append(modal_coefficients(grid, state.r))
    modes["Q"].append(modal_coefficients(grid, state.Q))
    modes["chi"].append(modal_coefficients(grid, state.chi))


def _series_map(samples):
    data = {key: np.asarray([float(sample[key]) for sample in samples], dtype=float) for key in SERIES_KEYS}
    for name in WINDOW_SERIES:
        data[name] = np.asarray([sample[name] for sample in samples], dtype=float)
    for index in range(4):
        data[f"window_normal_{index}"] = data["window_normal"][:, index]
        data[f"window_proper_flux_{index}"] = data["window_proper_flux"][:, index]
        data[f"window_quasilocal_{index}"] = data["window_quasilocal"][:, index]
    return data


def _checkpoint(sample, include_columns):
    kept = {
        key: sample[key]
        for key in (
            "time",
            "energy",
            "field_energy",
            "gravity_energy",
            "full_hamilton_max",
            "projected_hamilton_max",
            "held_out_hamilton_max",
            "full_momentum_max",
            "projected_momentum_max",
            "held_out_momentum_max",
            "mean_rho",
            "mean_current",
            "proper_min",
            "proper_max",
            "proper_mean",
            "lifted_proper_max",
            "chart_proper_max",
            "coordinate_over_N_max",
            "K_r_min",
            "K_r_max",
            "K_perp_min",
            "K_perp_max",
            "chi_min",
            "chi_max",
            "c2_max",
            "r_min",
            "r_max",
            "Q_min",
            "Q_max",
            "lifted_fieldwork",
            "proper_pressure_power",
            "coordinate_work_integral",
            "proper_pressure_work_integral",
            "unitarity",
            "positive_r",
            "positive_Q",
            "window_normal",
            "window_proper_flux",
            "window_coordinate_work",
            "window_proper_work",
            "window_quasilocal",
        )
    }
    if include_columns:
        kept["windows"] = []
        for row in sample["windows"]:
            copy = {key: row[key] for key in row if key != "mode_columns"}
            copy["mode_columns"] = row["mode_columns"]
            kept["windows"].append(copy)
    return kept


def json_checkpoints(samples):
    kept = []
    last = len(samples) - 1
    for index, sample in enumerate(samples):
        final = index == last
        if index == 0 or final or _frame_time(sample["time"], False):
            kept.append(_checkpoint(sample, include_columns=(index == 0 or final)))
    return kept


def effect_row(name, key, floor, primary, other, dt, requested):
    row = {
        "name": name,
        "series": key,
        "floor": float(floor),
        "requested_T": float(requested),
        "rule": "movement of episode changes <= 0.01 * primary episode change",
        "comparison_time": None,
        "primary_initial": None,
        "primary_change": None,
        "other_change": None,
        "effect_scale": None,
        "movement": None,
        "movement_over_effect": None,
        "terminal_difference": None,
        "status": "incomplete",
        "overlap_status": None,
    }
    if primary is None or other is None or key not in primary or key not in other:
        row["missing_series"] = key
        return row
    primary_attained = float(primary["time"][-1])
    other_attained = float(other["time"][-1])
    row["primary_attained_T"] = primary_attained
    row["other_attained_T"] = other_attained
    comparison = common_comparison_time(primary["time"], other["time"], dt)
    row["comparison_time"] = comparison
    if comparison is None:
        return row
    try:
        p0 = value_at(primary["time"], primary[key], 0.0)
        pT = value_at(primary["time"], primary[key], comparison)
        o0 = value_at(other["time"], other[key], 0.0)
        oT = value_at(other["time"], other[key], comparison)
    except KeyError:
        return row
    primary_change = pT - p0
    other_change = oT - o0
    effect = abs(primary_change)
    movement = abs(primary_change - other_change)
    overlap = classify_movement(effect, movement, floor)
    finished = primary_attained + 1e-12 >= requested and other_attained + 1e-12 >= requested
    row.update({
        "primary_initial": p0,
        "primary_at_comparison": pT,
        "other_initial": o0,
        "other_at_comparison": oT,
        "primary_change": primary_change,
        "other_change": other_change,
        "effect_scale": effect,
        "movement": movement,
        "movement_over_effect": None if effect == 0.0 else movement / effect,
        "terminal_difference": abs(pT - oT),
        "overlap_status": overlap,
        "status": overlap if finished else "incomplete",
    })
    return row


def _mode_array(modes, name):
    if not modes["time"]:
        return np.zeros(0), np.zeros((0, 0), dtype=complex)
    return np.asarray(modes["time"], dtype=float), np.stack(modes[name])


def shared_indices(modes_primary, modes_other):
    lookup = {float(mode): index for index, mode in enumerate(modes_other)}
    left = []
    right = []
    for index, mode in enumerate(modes_primary):
        key = float(mode)
        if key in lookup:
            left.append(index)
            right.append(lookup[key])
    return np.asarray(left, dtype=int), np.asarray(right, dtype=int)


def modal_effect_row(name, key, floor, primary_modes, other_modes, primary_grid_modes, other_grid_modes, requested):
    row = {
        "name": name,
        "series": key,
        "floor": float(floor),
        "requested_T": float(requested),
        "rule": "shared-mode episode-change movement <= 0.01 * primary change",
        "status": "incomplete",
        "effect_scale": None,
        "movement": None,
    }
    if primary_modes is None or other_modes is None:
        return row
    time_p, coeff_p = _mode_array(primary_modes, key)
    time_o, coeff_o = _mode_array(other_modes, key)
    if time_p.size == 0 or time_o.size == 0:
        return row
    if float(time_p[-1]) + 1e-12 < requested or float(time_o[-1]) + 1e-12 < requested:
        row["primary_attained_T"] = float(time_p[-1])
        row["other_attained_T"] = float(time_o[-1])
        return row
    left, right = shared_indices(primary_grid_modes, other_grid_modes)
    delta_p = coeff_p[-1, left] - coeff_p[0, left]
    delta_o = coeff_o[-1, right] - coeff_o[0, right]
    effect = float(np.linalg.norm(delta_p))
    movement = float(np.linalg.norm(delta_p - delta_o))
    row.update({
        "comparison_time": float(time_p[-1]),
        "shared_mode_count": int(left.size),
        "initial_shared_l2": float(np.linalg.norm(coeff_p[0, left] - coeff_o[0, right])),
        "effect_scale": effect,
        "movement": movement,
        "movement_over_effect": None if effect == 0.0 else movement / effect,
        "primary_change": effect,
        "status": classify_movement(effect, movement, floor),
        "overlap_status": classify_movement(effect, movement, floor),
    })
    return row


def build_effect_table(runs, grids):
    primary = runs.get("nf512_dt_0_0005")
    time_other = runs.get("nf512_dt_0_00025")
    space_other = runs.get("nf256_dt_0_0005")
    physical = []
    diagnostic = []
    for name, key, floor in PHYSICAL_EFFECTS:
        if key.startswith("mode_"):
            continue
        physical.append({
            "time": effect_row(name, key, floor, primary, time_other, DT, DURATION),
            "space": effect_row(name, key, floor, primary, space_other, DT, DURATION),
        })
    for name, key, floor in DIAGNOSTIC_EFFECTS:
        diagnostic.append({
            "time": effect_row(name, key, floor, primary, time_other, DT, DURATION),
            "space": effect_row(name, key, floor, primary, space_other, DT, DURATION),
        })
    if primary is not None and "nf512" in grids:
        for name, key, floor in PHYSICAL_EFFECTS:
            if not key.startswith("mode_"):
                continue
            mode_key = key.split("_", 1)[1]
            physical.append({
                "time": modal_effect_row(
                    name, mode_key, floor,
                    primary.get("modes"), None if time_other is None else time_other.get("modes"),
                    grids["nf512"].modes_g, grids["nf512"].modes_g, DURATION,
                ),
                "space": modal_effect_row(
                    name, mode_key, floor,
                    primary.get("modes"), None if space_other is None else space_other.get("modes"),
                    grids["nf512"].modes_g,
                    grids["nf256"].modes_g if "nf256" in grids else np.zeros(0),
                    DURATION,
                ),
            })
    return physical, diagnostic


def _statuses(table):
    statuses = []
    for item in table:
        statuses.append(item["time"]["status"])
        statuses.append(item["space"]["status"])
    return statuses


def verdict_from(physical_summary, diagnostic_summary, runs_complete, witness_ok, chart_stop, exception, budget_stop):
    """Label the episode without letting a constraint rollup veto physical rows.

    The one-percent test belongs to each claimed physical observable. A
    constraint residual, including a near-zero total-energy drift, stays in
    the diagnostic report. It becomes a veto only through an explicit
    sensitivity map, which this classifier does not invent. A constraint-
    consistent solution stays false unless a real bound is supplied.
    """
    assessment = {
        "physical_summary": physical_summary,
        "diagnostic_summary": diagnostic_summary,
        "constraint_consistent_solution": False,
        "constraint_solution_requires_real_bound": True,
        "continuum_error": "unknown",
        "sensitivity": None,
        "diagnostic_self_scale_is_physical_veto": False,
    }
    if exception:
        assessment.update(verdict="FAILED_EXCEPTION", effect_resolved=False)
        return assessment
    if chart_stop:
        assessment.update(verdict="STOP_POSITIVE_CHART", effect_resolved=False)
        return assessment
    if not witness_ok:
        assessment.update(verdict="SAME_MODEL_WITNESS_FAILED", effect_resolved=False)
        return assessment
    if not runs_complete:
        assessment.update(
            verdict="PARTIAL_CPU_BUDGET" if budget_stop else "PARTIAL_INCOMPLETE",
            effect_resolved=False,
        )
        return assessment
    if physical_summary == "resolved":
        assessment.update(
            verdict="MEASURED_FEEDBACK_UNRESOLVED_CONSTRAINT_CONTROL",
            effect_resolved=True,
            constraint_control="unresolved",
        )
        return assessment
    if physical_summary == "no_measurable_evolution":
        assessment.update(
            verdict="EPISODE_MEASURED_NO_MEASURABLE_EVOLUTION",
            effect_resolved=False,
            constraint_control="unresolved",
        )
        return assessment
    assessment.update(
        verdict="MEASURED_FEEDBACK_PER_OBSERVABLE",
        effect_resolved=None,
        per_observable_indicators_stand=True,
        constraint_control="unresolved",
    )
    return assessment


def same_model_witness(samples):
    """Compare the opening of this episode with the saved v5 nf512 window."""
    payload = json.loads(V5_JSON.read_text())
    checkpoints = payload["nf512"]["window_dt_0_0005"]["checkpoints"]
    gaps = {key: 0.0 for key in WITNESS_KEYS}
    compared = 0
    for checkpoint in checkpoints:
        time_value = float(checkpoint["time"])
        if time_value > DURATION + 1e-12:
            continue
        try:
            mine = next(sample for sample in samples if abs(sample["time"] - time_value) <= 1e-12)
        except StopIteration:
            continue
        compared += 1
        for key in WITNESS_KEYS:
            saved_key = "r_min_quadrature" if key == "r_min" else "Q_min_quadrature" if key == "Q_min" else key
            gaps[key] = max(gaps[key], abs(float(mine[key]) - float(checkpoint[saved_key])))
    recorded_work = float(payload["nf512"]["window_dt_0_0005"]["work_fieldwork_integral"])
    mine_at = next((sample for sample in samples if abs(sample["time"] - 0.005) <= 1e-12), None)
    work_gap = None if mine_at is None else abs(float(mine_at["coordinate_work_integral"]) - recorded_work)
    worst = max(gaps.values()) if compared else None
    attained = float(samples[-1]["time"]) if samples else 0.0
    full_window = attained + 1e-12 >= 0.005
    matches = bool(
        compared > 0
        and worst is not None
        and worst <= 1e-8
        and (work_gap is None or work_gap <= 1e-8)
        and (not full_window or (compared == 11 and work_gap is not None and work_gap <= 1e-8))
    )
    return {
        "compared_checkpoints": compared,
        "expected_checkpoints": 11,
        "attained_T": attained,
        "max_abs_gap": worst,
        "gaps": gaps,
        "work_integral_gap_at_T_0_005": work_gap,
        "matches_saved_v5_window": matches,
    }


def _state_arrays(prefix, state):
    return {
        prefix + "Q": np.asarray(state.Q),
        prefix + "r": np.asarray(state.r),
        prefix + "chi": np.asarray(state.chi),
        prefix + "p_Q": np.asarray(state.p_Q),
        prefix + "p_r": np.asarray(state.p_r),
        prefix + "p_chi": np.asarray(state.p_chi),
        prefix + "phi0": np.asarray(state.phi0),
        prefix + "phi1": np.asarray(state.phi1),
    }


def _rate_arrays(prefix, rate):
    if rate is None:
        return {}
    return {
        prefix + "Q_dot": np.asarray(rate.Q),
        prefix + "r_dot": np.asarray(rate.r),
        prefix + "chi_dot": np.asarray(rate.chi),
        prefix + "p_Q_dot": np.asarray(rate.p_Q),
        prefix + "p_r_dot": np.asarray(rate.p_r),
        prefix + "p_chi_dot": np.asarray(rate.p_chi),
        prefix + "phi0_dot": np.asarray(rate.phi0),
        prefix + "phi1_dot": np.asarray(rate.phi1),
    }


def _series_arrays(prefix, samples):
    data = _series_map(samples)
    arrays = {prefix + key: value for key, value in data.items()}
    return arrays


def _frame_arrays(prefix, frames):
    if not frames:
        return {}
    arrays = {prefix + "frame_time": np.asarray([frame["time"] for frame in frames], dtype=float)}
    for name in ("r", "Q", "chi", "proper", "K_r", "K_perp", "rho", "current", "weyl_C2"):
        arrays[prefix + "frame_" + name] = np.stack([frame["fields"][name] for frame in frames])
    return arrays


def _mode_arrays(prefix, modes):
    if not modes["time"]:
        return {}
    arrays = {prefix + "mode_time": np.asarray(modes["time"], dtype=float)}
    for name in ("r", "Q", "chi"):
        arrays[prefix + "mode_" + name] = np.stack(modes[name])
    return arrays


def run_summary(result):
    samples = result["samples"]
    initial = samples[0]
    final = samples[-1]
    return {
        "attained_T": result["attained_time"],
        "requested_T": result["requested_time"],
        "completed": result["completed"],
        "stop_reason": result["stop_reason"],
        "chart_stop": result["chart_stop"],
        "dt": result["dt"],
        "steps_requested": result["steps_requested"],
        "steps_completed": result["steps_completed"],
        "cpu_seconds": result["cpu_seconds"],
        "positive_r": result["positive_r"],
        "positive_Q": result["positive_Q"],
        "initial": _checkpoint(initial, include_columns=True),
        "final": _checkpoint(final, include_columns=True),
        "checkpoints": json_checkpoints(samples),
        "proper_max_change": float(final["proper_max"] - initial["proper_max"]),
        "chi_max_change": float(final["chi_max"] - initial["chi_max"]),
        "r_min_change": float(final["r_min"] - initial["r_min"]),
        "coordinate_work_integral": float(final["coordinate_work_integral"]),
        "proper_pressure_work_integral": float(final["proper_pressure_work_integral"]),
        "full_hamilton_initial": float(initial["full_hamilton_max"]),
        "full_hamilton_final": float(final["full_hamilton_max"]),
        "full_momentum_initial": float(initial["full_momentum_max"]),
        "full_momentum_final": float(final["full_momentum_max"]),
        "mean_rho_initial": float(initial["mean_rho"]),
        "mean_rho_final": float(final["mean_rho"]),
        "mean_current_initial": float(initial["mean_current"]),
        "mean_current_final": float(final["mean_current"]),
        "mean_subtracted": False,
    }


def _binding_record(metadata):
    return {
        "names": [str(name) for name in metadata["binding_names"]],
        "values": [float(value) for value in metadata["binding_values"]],
        "occupation_gap": metadata["occupation_gap"],
        "quadrature_radius_gap": metadata["quadrature_radius_gap"],
        "v5_npz_sha256": metadata["v5_sha256"],
    }


def _setup_record(grid, state, metadata, sample):
    return {
        "nf": metadata["nf"],
        "ng": metadata["ng"],
        "nq": metadata["nq"],
        "length": metadata["length"],
        "preparation_source": "saved v5 npz state, not a new Newton solve",
        "frozen_rho_used_in_evolution": False,
        "binding": _binding_record(metadata),
        "historical_initial_1e-8": initial_tolerance_decision(
            sample["full_hamilton_max"], sample["full_momentum_max"]
        ),
        "initial_observation": _checkpoint(sample, include_columns=True),
        "state_copied_from_v5": True,
        "phi_shape": [int(state.phi0.shape[0]), int(state.phi0.shape[1])],
        "grid_nf": int(grid.nf),
    }


def saved_rho_gap(grid, state, metadata):
    _sample, _rate, fields = observe(grid, state)
    return float(np.max(np.abs(fields["rho"] - metadata["saved_rho"]))), _sample


def blank_record():
    return {
        "schema": SCHEMA,
        "status": "RUNNING",
        "verdict": "RUNNING",
        "renewal": False,
        "continuum_limit_claimed": False,
        "frozen_B_embedding_claimed": False,
        "embedding_match_required": False,
        "vacuum_branch_included": False,
        "incoming_gate": "OPEN",
        "incoming_gate_campaign": "paused",
        "incoming_gate_prerequisite": False,
        "lambda_cdm_conflict_claimed": False,
        "inheritance_T_star_recomputed": False,
        "inheritance_statement": "T* D = D is reused. This episode does not compute an infinite nest.",
        "question": (
            "From the saved v5 autonomous state, what evolves through T=0.05, "
            "and do the episode changes agree within one percent under dt halving "
            "at nf512 and under comparison with nf256?"
        ),
        "model": {
            "driver": "derive_nsc_spherical_feedback_episode.py",
            "step": "nsc_spherical_galerkin_coupling.rk4_step",
            "phi_dot": "U_f† (-i H_fine U_f Phi)",
            "source_geometry_H_every_rk4_stage": True,
            "hamiltonian_cached_across_stages": False,
            "include_matter_force": True,
            "post_step_filter": False,
            "constraint_projection": False,
            "prescribed_radius": False,
            "state_reset": False,
            "restoring_force": False,
            "static_gauge": "L and beta stay the calibration samples",
            "evolving": ["Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1"],
            "mean_source_subtracted": False,
            "mean_current_subtracted": False,
            "curvature": "chi = R_h - 2; C^2 from weyl_square_from_rh; Riemann tensor not recomputed",
        },
        "budget": {
            "cpu_seconds": CPU_BUDGET_S,
            "includes_all_intended_runs": True,
            "unbudgeted_refinement": False,
            "intended_runs": [dict(item) for item in INTENDED_RUNS],
        },
        "resolution_rule": (
            "effect_scale is the absolute episode change on nf512 at dt=0.0005. "
            "movement is the absolute difference of the episode changes. "
            "resolved means the requested T was reached and movement <= 0.01 * effect_scale, "
            "with effect_scale at or above the predeclared floor. "
            "A change below the floor on both runs is below_noise_floor, not a resolved effect. "
            "Disagreement above the floor stays unresolved."
        ),
        "floors": {
            "physical": {name: floor for name, _key, floor in PHYSICAL_EFFECTS},
            "diagnostic": {name: floor for name, _key, floor in DIAGNOSTIC_EFFECTS},
        },
        "thread_limits": {
            "OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"),
            "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
            "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
            "NUMEXPR_NUM_THREADS": os.environ.get("NUMEXPR_NUM_THREADS"),
            "VECLIB_MAXIMUM_THREADS": os.environ.get("VECLIB_MAXIMUM_THREADS"),
            "threadpoolctl": THREADPOOL_LIMIT,
        },
        "runs": {},
        "costs": {},
    }


def run_episode():
    limit_threads()
    if OUT.exists() or NPZ.exists():
        raise FileExistsError(f"refusing to overwrite episode evidence: {OUT.name} or {NPZ.name}")
    started_wall = time.perf_counter()
    started_cpu = time.process_time()
    elapsed = lambda: time.process_time() - started_cpu
    record = blank_record()
    record["hashes_before"] = frozen_hashes()
    record["hashes_before"]["driver"] = sha256(DRIVER)
    arrays = {
        "schema": np.array(SCHEMA),
        "duration": np.float64(DURATION),
        "dt": np.float64(DT),
        "dt_half": np.float64(DT_HALF),
    }
    exception_text = None
    results = {}
    series = {}
    grids = {}
    last_step_cost = 0.0
    last_nq = 0
    try:
        loaded = {}
        for label in ("nf512", "nf256"):
            if not budget_allows(elapsed(), 0.0):
                record["setup_stop"] = "CPU_BUDGET_300S"
                break
            grid, state, metadata = load_v5_state(label)
            rho_gap, sample = saved_rho_gap(grid, state, metadata)
            sample["time"] = 0.0
            sample["coordinate_work_integral"] = 0.0
            sample["proper_pressure_work_integral"] = 0.0
            loaded[label] = (grid, state, metadata, sample, rho_gap)
            grids[label] = grid
            arrays.update(_state_arrays(label + "_initial_", state))
            arrays[label + "_binding_names"] = np.array([str(item) for item in metadata["binding_names"]])
            arrays[label + "_binding_values"] = np.asarray(metadata["binding_values"], dtype=float)
            arrays[label + "_occupations"] = np.asarray(metadata["occupations"], dtype=float)
            record["setup_" + label] = _setup_record(grid, state, metadata, sample)
            record["setup_" + label]["saved_rho_gap"] = rho_gap
            print("loaded", label, "rho_gap", rho_gap, "full_H", sample["full_hamilton_max"], flush=True)
        write_json(OUT, record)
        for spec in INTENDED_RUNS:
            name = spec["name"]
            if spec["resolution"] not in loaded:
                record["runs"][name] = {"completed": False, "stop_reason": "SETUP_MISSING", "attained_T": 0.0}
                break
            grid, state, _metadata, _sample, _gap = loaded[spec["resolution"]]
            estimated = last_step_cost * (grid.nq / last_nq) ** 2 if last_nq else last_step_cost
            if not budget_allows(elapsed(), estimated):
                record["runs"][name] = {
                    "completed": False,
                    "stop_reason": "CPU_BUDGET_300S",
                    "attained_T": 0.0,
                    "cpu_seconds_at_refusal": float(elapsed()),
                    "estimated_step_cpu": float(estimated),
                }
                continue
            print("start", name, "elapsed_cpu", elapsed(), flush=True)
            try:
                result = evolve_episode(
                    grid, state.copy(), DURATION, spec["dt"], elapsed,
                    budget=CPU_BUDGET_S, last_step_cost=last_step_cost, label=name,
                )
            except Exception:
                exception_text = traceback.format_exc()
                record["runs"][name] = {
                    "completed": False,
                    "stop_reason": "CODE_EXCEPTION",
                    "attained_T": 0.0,
                    "error": exception_text,
                }
                break
            last_step_cost = result["step_cost"]
            last_nq = int(result["last_nq"])
            results[name] = result
            packed = _series_map(result["samples"])
            packed["modes"] = result["modes"]
            series[name] = packed
            summary = run_summary(result)
            if result["error"]:
                summary["error"] = result["error"]
                exception_text = result["error"]
            record["runs"][name] = summary
            record["costs"][name] = {"cpu_seconds": result["cpu_seconds"]}
            arrays.update(_series_arrays(name + "_", result["samples"]))
            arrays.update(_frame_arrays(name + "_", result["frames"]))
            arrays.update(_mode_arrays(name + "_", result["modes"]))
            arrays.update(_state_arrays(name + "_final_", result["state"]))
            arrays.update(_rate_arrays(name + "_final_rate_", result["last_rate"]))
            if name == "nf512_dt_0_0005" and result["state_after_one_step"] is not None:
                arrays.update(_state_arrays("replay_", result["state_after_one_step"]))
                arrays.update(_rate_arrays("replay_rate_", result["rate_after_one_step"]))
                arrays["replay_time"] = np.float64(result["dt"])
                record["replay"] = {
                    "run": name,
                    "time": result["dt"],
                    "stored": True,
                }
                record["same_model_witness"] = same_model_witness(result["samples"])
                print(
                    "witness",
                    record["same_model_witness"]["matches_saved_v5_window"],
                    "gap",
                    record["same_model_witness"]["max_abs_gap"],
                    "work_gap",
                    record["same_model_witness"]["work_integral_gap_at_T_0_005"],
                    flush=True,
                )
            write_json(OUT, record)
            write_npz(NPZ, arrays)
            print("done", name, "T", result["attained_time"], "reason", result["stop_reason"], flush=True)
            if result["stop_reason"] == "CODE_EXCEPTION":
                break
    except Exception:
        exception_text = traceback.format_exc()
        record["campaign_error"] = exception_text
    physical, diagnostic = build_effect_table(series, grids)
    record["physical_effects"] = physical
    record["constraint_diagnostics"] = diagnostic
    physical_summary = summarize_statuses(_statuses(physical))
    diagnostic_summary = summarize_statuses(_statuses(diagnostic))
    runs_complete = all(
        bool(record["runs"].get(spec["name"], {}).get("completed")) for spec in INTENDED_RUNS
    )
    chart_stop = any(bool(record["runs"].get(spec["name"], {}).get("chart_stop")) for spec in INTENDED_RUNS)
    budget_stop = any(
        record["runs"].get(spec["name"], {}).get("stop_reason") == "CPU_BUDGET_300S" for spec in INTENDED_RUNS
    ) or record.get("setup_stop") == "CPU_BUDGET_300S"
    witness = record.get("same_model_witness") or {}
    if "nf512_dt_0_0005" not in series:
        witness_ok = True
    else:
        witness_ok = bool(witness.get("matches_saved_v5_window"))
    exception = exception_text is not None or any(
        record["runs"].get(spec["name"], {}).get("stop_reason") == "CODE_EXCEPTION" for spec in INTENDED_RUNS
    )
    assessment = verdict_from(
        physical_summary, diagnostic_summary, runs_complete, witness_ok, chart_stop, exception, budget_stop
    )
    record["physical_effect_summary"] = physical_summary
    record["constraint_diagnostic_summary"] = diagnostic_summary
    record["runs_complete"] = runs_complete
    record["assessment"] = assessment
    record["verdict"] = assessment["verdict"]
    record["effect_resolved"] = assessment["effect_resolved"] is True
    record["constraint_consistent_solution"] = False
    record["renewal"] = False
    label = assessment["verdict"]
    record["status"] = "PARTIAL" if label.startswith("PARTIAL") or label.startswith("FAILED") or label.startswith("STOP") or label.startswith("SAME_MODEL") else "MEASURED"
    record["answer"] = _answer(record, physical, diagnostic)
    record["cpu_seconds"] = float(elapsed())
    record["wall_seconds"] = float(time.perf_counter() - started_wall)
    record["answer"]["cpu_seconds"] = record["cpu_seconds"]
    record["hashes_after"] = frozen_hashes()
    record["hashes_after"]["driver"] = sha256(DRIVER)
    record["frozen_bytes_unchanged"] = _frozen_same(record)
    write_json(OUT, record)
    arrays["verdict"] = np.array(record["verdict"])
    write_npz(NPZ, arrays)
    payload = OUT.stat().st_size + NPZ.stat().st_size
    record["payload_bytes"] = int(payload)
    record["payload_within_64MiB"] = bool(payload <= PAYLOAD_LIMIT)
    record["npz_sha256"] = sha256(NPZ)
    record["json_sha256_without_self_hash"] = None
    write_json(OUT, record)
    if payload > PAYLOAD_LIMIT:
        record["verdict"] = "FAILED_PAYLOAD_LIMIT"
        record["effect_resolved"] = False
        write_json(OUT, record)
    print("verdict", record["verdict"], "cpu", record["cpu_seconds"], "payload", payload, flush=True)
    return record


def _frozen_same(record):
    before = dict(record["hashes_before"])
    after = dict(record["hashes_after"])
    before.pop("driver", None)
    after.pop("driver", None)
    return before == after


def _answer(record, physical, diagnostic):
    evolving = []
    unresolved = []
    resolved = []
    quiet = []
    for item in physical:
        for axis_name in ("time", "space"):
            row = item[axis_name]
            label = row["name"] + ":" + axis_name
            status = row["status"]
            if status == "resolved":
                resolved.append(label)
            elif status == "unresolved":
                unresolved.append(label)
            elif status == "below_noise_floor":
                quiet.append(label)
            if row.get("effect_scale") is not None and row["effect_scale"] >= row["floor"]:
                evolving.append({
                    "name": label,
                    "effect_scale": row["effect_scale"],
                    "movement": row["movement"],
                    "status": status,
                })
    constraint_unresolved = [
        row["name"] + ":" + axis
        for item in diagnostic
        for axis, row in item.items()
        if row["status"] == "unresolved"
    ]
    attained = {
        name: record["runs"].get(name, {}).get("attained_T")
        for name in (spec["name"] for spec in INTENDED_RUNS)
    }
    return {
        "attained_T": attained,
        "cpu_seconds": None,
        "physical_effect_summary": record.get("physical_effect_summary"),
        "constraint_diagnostic_summary": record.get("constraint_diagnostic_summary"),
        "effect_resolved": record.get("effect_resolved"),
        "constraint_consistent_solution": False,
        "verdict": record.get("verdict"),
        "renewal": False,
        "what_changes_above_the_floor": evolving,
        "resolved": resolved,
        "unresolved": unresolved,
        "below_noise_floor": quiet,
        "constraint_changes_unresolved": constraint_unresolved,
        "same_model_witness": (record.get("same_model_witness") or {}).get("matches_saved_v5_window"),
        "historical_initial_1e-8_aborted_the_run": False,
    }


def refinement_indicators(physical_effects):
    """One time indicator and one space indicator for each claimed physical effect."""
    indicators = []
    for item in physical_effects:
        time_row = item["time"]
        space_row = item["space"]
        indicators.append({
            "name": time_row["name"],
            "time_status": time_row["status"],
            "space_status": space_row["status"],
            "time_effect_scale": time_row.get("effect_scale"),
            "time_movement": time_row.get("movement"),
            "time_movement_over_effect": time_row.get("movement_over_effect"),
            "space_effect_scale": space_row.get("effect_scale"),
            "space_movement": space_row.get("movement"),
            "space_movement_over_effect": space_row.get("movement_over_effect"),
            "meets_one_percent": time_row["status"] in ("resolved", "below_noise_floor") and space_row["status"] in ("resolved", "below_noise_floor"),
        })
    return indicators


def energy_balance_from_run(run):
    """Field and gravity changes against the exchange already stored on the run."""
    initial = run["initial"]
    final = run["final"]
    field_change = float(final["field_energy"]) - float(initial["field_energy"])
    gravity_change = float(final["gravity_energy"]) - float(initial["gravity_energy"])
    total_change = float(final["energy"]) - float(initial["energy"])
    work = float(final["coordinate_work_integral"])
    return {
        "field_energy_initial": float(initial["field_energy"]),
        "field_energy_final": float(final["field_energy"]),
        "gravity_energy_initial": float(initial["gravity_energy"]),
        "gravity_energy_final": float(final["gravity_energy"]),
        "field_energy_change": field_change,
        "gravity_energy_change": gravity_change,
        "field_plus_gravity_change": field_change + gravity_change,
        "total_energy_change": total_change,
        "coordinate_work_integral": work,
        "proper_pressure_work_integral": float(final["proper_pressure_work_integral"]),
        "balance_versus_field_exchange": balance_against_exchange(total_change, field_change),
        "balance_versus_coordinate_work": balance_against_exchange(total_change, work),
    }


def _constraint_control(record, balance):
    primary = record["runs"]["nf512_dt_0_0005"]
    coarse = record["runs"].get("nf256_dt_0_0005", {})
    residuals = []
    for item in record.get("constraint_diagnostics", []):
        for axis in ("time", "space"):
            row = item[axis]
            residuals.append({
                "name": row["name"],
                "axis": axis,
                "self_scale_status": row["status"],
                "effect_scale": row.get("effect_scale"),
                "movement": row.get("movement"),
                "movement_over_effect": row.get("movement_over_effect"),
                "primary_initial": row.get("primary_initial"),
                "primary_at_comparison": row.get("primary_at_comparison"),
                "other_at_comparison": row.get("other_at_comparison"),
                "vetoes_physical_effects": False,
            })
    return {
        "consistent_solution_validated": False,
        "continuum_error": "unknown",
        "full_projected_and_held_out_residuals_retained": True,
        "spatial_constraint_gap_retained": True,
        "full_hamilton_initial": primary.get("full_hamilton_initial"),
        "full_hamilton_final": primary.get("full_hamilton_final"),
        "full_momentum_initial": primary.get("full_momentum_initial"),
        "full_momentum_final": primary.get("full_momentum_final"),
        "nf256_full_hamilton_initial": coarse.get("full_hamilton_initial"),
        "nf256_full_hamilton_final": coarse.get("full_hamilton_final"),
        "nf256_full_momentum_final": coarse.get("full_momentum_final"),
        "diagnostic_rows": residuals,
        "energy_balance": balance,
        "statement": (
            "Full, projected, and held-out residuals remain, and the nf256 "
            "Hamilton residual stays far above the nf512 residual. Continuum "
            "error is unknown. No real bound validates a constraint-consistent "
            "solution. A one-percent test of those residual changes, or of the "
            "near-zero total-energy drift, does not veto the physical effects."
        ),
    }


def reassess_draft():
    """Rewrite assessment metadata from the saved four-run record.

    Does not integrate, does not rebuild states, and does not edit the v5 files.
    The uncommitted v1 JSON and NPZ are still drafts, so this correction is
    allowed before their first commit.
    """
    if not OUT.is_file() or not NPZ.is_file():
        raise FileNotFoundError("draft episode evidence is missing")
    record = json.loads(OUT.read_text())
    previous_verdict = record.get("verdict")
    campaign_driver = (record.get("hashes_after") or {}).get("driver")
    indicators = refinement_indicators(record["physical_effects"])
    balance = {
        name: energy_balance_from_run(run)
        for name, run in record["runs"].items()
        if run.get("initial") and run.get("final")
    }
    for item in record.get("constraint_diagnostics", []):
        for axis in ("time", "space"):
            item[axis]["vetoes_physical_effects"] = False
            item[axis]["self_scale_is_not_the_physical_test"] = True
    physical_summary = summarize_statuses(_statuses(record["physical_effects"]))
    previous_diagnostic = record.get("constraint_diagnostic_summary")
    assessment = verdict_from(physical_summary, previous_diagnostic, True, True, False, False, False)
    record["physical_effect_summary"] = physical_summary
    record["constraint_diagnostic_summary"] = "reported_not_a_physical_veto"
    record["assessment"] = assessment
    record["verdict"] = assessment["verdict"]
    record["effect_resolved"] = assessment["effect_resolved"] is True
    record["constraint_consistent_solution"] = False
    record["refinement_indicators"] = indicators
    record["energy_balance"] = balance
    record["constraint_control"] = _constraint_control(record, balance.get("nf512_dt_0_0005"))
    record["answer"]["effect_resolved"] = record["effect_resolved"]
    record["answer"]["verdict"] = record["verdict"]
    record["answer"]["constraint_consistent_solution"] = False
    record["answer"]["physical_effect_summary"] = physical_summary
    record["answer"]["constraint_diagnostic_summary"] = record["constraint_diagnostic_summary"]
    record["answer"]["unresolved"] = [
        row["name"] + ":" + axis
        for item in record["physical_effects"]
        for axis, row in (("time", item["time"]), ("space", item["space"]))
        if row["status"] == "unresolved"
    ]
    record["status"] = "MEASURED"
    record["assessment_revision"] = {
        "revised_without_rerunning_evolution": True,
        "numerical_values_recomputed": False,
        "historical_v5_json_changed": False,
        "historical_v5_npz_changed": False,
        "previous_verdict": previous_verdict,
        "campaign_driver_sha256": campaign_driver,
        "reason": (
            "The draft applied the one-percent target to constraint-diagnostic "
            "changes, including a total-energy drift of about 3.44e-8 beside an "
            "energy exchange of about 0.4123, and that rollup vetoed the physical "
            "effects. The one-percent target applies to claimed physical effects. "
            "Constraint residuals stay reported and do not by themselves unresolve "
            "those effects. This revision rewrites the assessment from the saved "
            "four-run arrays before the first commit of the draft."
        ),
    }
    record["hashes_after"].pop("episode_tests", None)
    record["hashes_after"]["driver"] = sha256(DRIVER)
    record["assessment_revision"]["reassessed_driver_sha256"] = record["hashes_after"]["driver"]
    record["assessment_revision"]["reassessed_test_sha256"] = sha256(
        LAB / "tests" / "test_nsc_spherical_feedback_episode.py"
    )
    with np.load(NPZ, allow_pickle=False) as stored:
        arrays = {key: np.array(stored[key]) for key in stored.files}
    numerical = {
        key: np.array(value, copy=True)
        for key, value in arrays.items()
        if key != "verdict"
    }
    arrays["verdict"] = np.array(record["verdict"])
    write_npz(NPZ, arrays)
    with np.load(NPZ, allow_pickle=False) as checked:
        for key, value in numerical.items():
            if not np.array_equal(checked[key], value):
                raise RuntimeError(f"reassessment changed numerical array {key}")
        if str(checked["verdict"]) != record["verdict"]:
            raise RuntimeError("verdict metadata was not stored")
    record["npz_sha256"] = sha256(NPZ)
    record["payload_bytes"] = int(OUT.stat().st_size + NPZ.stat().st_size)
    write_json(OUT, record)
    record["payload_bytes"] = int(OUT.stat().st_size + NPZ.stat().st_size)
    record["payload_within_64MiB"] = bool(record["payload_bytes"] <= PAYLOAD_LIMIT)
    write_json(OUT, record)
    return record


def verify_saved(replay=True):
    """Replay the stored first step and check the record against its payload."""
    limit_threads()
    if not OUT.is_file() or not NPZ.is_file():
        raise FileNotFoundError("episode evidence is missing")
    record = json.loads(OUT.read_text())
    payload = OUT.stat().st_size + NPZ.stat().st_size
    if payload > PAYLOAD_LIMIT:
        raise AssertionError(f"payload {payload} exceeds 64MiB")
    if record["schema"] != SCHEMA:
        raise AssertionError("schema mismatch")
    if record["renewal"] is not False:
        raise AssertionError("renewal must stay false")
    if record["model"]["mean_source_subtracted"] is not False:
        raise AssertionError("mean source was subtracted")
    if record["model"]["source_geometry_H_every_rk4_stage"] is not True:
        raise AssertionError("stage recomputation was not declared")
    if record["incoming_gate_prerequisite"] is not False:
        raise AssertionError("incoming gate was treated as a prerequisite")
    if not _frozen_same(record):
        raise AssertionError("frozen bytes changed during the episode")
    current_frozen = frozen_hashes()
    for name, digest in current_frozen.items():
        if record["hashes_after"].get(name) != digest:
            raise AssertionError(f"frozen hash drifted after the run: {name}")
    with np.load(NPZ, allow_pickle=False) as episode:
        with np.load(V5_NPZ, allow_pickle=False) as saved:
            for episode_key, saved_key in (
                ("nf512_initial_r", "nf512_geometry_r"),
                ("nf512_initial_Q", "nf512_geometry_Q"),
                ("nf512_initial_chi", "nf512_geometry_chi"),
                ("nf512_initial_phi0", "nf512_columns_phi0"),
                ("nf512_initial_phi1", "nf512_columns_phi1"),
                ("nf256_initial_r", "geometry_r"),
                ("nf256_initial_phi0", "columns_phi0"),
            ):
                if not np.array_equal(episode[episode_key], saved[saved_key]):
                    raise AssertionError(f"setup does not copy {saved_key}")
        if replay:
            if "replay_r" not in episode:
                raise AssertionError("replay state was not stored")
            grid, state, _metadata = load_v5_state("nf512")
            stepped = galerkin.rk4_step(grid, state, DT)
            for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1"):
                if not np.allclose(getattr(stepped, name), episode["replay_" + name], rtol=0.0, atol=1e-12):
                    raise AssertionError(f"replay mismatch in {name}")
            if np.array_equal(stepped.Q, state.Q):
                raise AssertionError("first step left Q unchanged")
        primary = "nf512_dt_0_0005_"
        if primary + "time" in episode and record["runs"].get("nf512_dt_0_0005"):
            attained = float(episode[primary + "time"][-1])
            reported = float(record["runs"]["nf512_dt_0_0005"]["attained_T"])
            if abs(attained - reported) > 1e-12:
                raise AssertionError("attained T does not match the series")
            for item in record.get("physical_effects", []):
                for axis in ("time", "space"):
                    row = item[axis]
                    if row["series"].startswith("mode_"):
                        continue
                    if row["status"] == "incomplete" or row.get("effect_scale") is None:
                        continue
                    fresh = classify_movement(row["effect_scale"], row["movement"], row["floor"])
                    if row["status"] != fresh:
                        raise AssertionError(f"stored status disagrees with the rule for {row['name']}")
                    change = abs(row["primary_change"])
                    if abs(change - row["effect_scale"]) > 1e-12:
                        raise AssertionError("effect scale is not the episode change")
    decision = record["setup_nf512"]["historical_initial_1e-8"]
    if decision["abort_evolution"] is not False or decision["continue_evolution"] is not True:
        raise AssertionError("historical 1e-8 flag aborted or blocked the episode")
    return record


def main(argv=None):
    arguments = list(sys.argv[1:] if argv is None else argv)
    if "--reassess" in arguments:
        record = reassess_draft()
        print(json.dumps({
            "verdict": record["verdict"],
            "effect_resolved": record["effect_resolved"],
            "constraint_consistent_solution": record["constraint_consistent_solution"],
            "physical_indicators": len(record["refinement_indicators"]),
        }, indent=2))
        return 0
    if "--check" in arguments:
        record = verify_saved(replay=True)
        print(json.dumps({
            "verdict": record["verdict"],
            "effect_resolved": record["effect_resolved"],
            "attained_T": record["answer"]["attained_T"],
            "cpu_seconds": record["cpu_seconds"],
            "payload_bytes": record["payload_bytes"],
        }, indent=2))
        return 0
    run_episode()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
