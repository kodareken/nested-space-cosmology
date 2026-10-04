#!/usr/bin/env python3
"""Sole campaign CLI for the classical RN radial parent.

Default mode is read-only preview. Calibration measures a source-free PG
reconstruction, its rates, its metric jets and one RK4 step against the
reference RN. The pass mark is 1e-4 on those measurements. A boolean cannot
grant admission. Neutral evolution calls the PG and source owners when their
callables exist, after that immutable binding and the CPU admission. A missing
callable is an adapter gap and is not filled with a manufactured state.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import subprocess

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import numpy as np

from recursive_horizons import nsc_rn_observables as observables
from recursive_horizons import nsc_rn_reference as reference

SCHEMA = "NSC-RN-PARENT-v2"
API_VERSION = "nsc-rn-parent-driver-v2"
LAB = Path(__file__).resolve().parents[1]
ROOT = LAB.parent
MASS_OVER_RM = (1.002, 1.01, 1.04)
ENERGY_FRACTIONS = (1e-3, 1e-2)
STATION_TIMES_OVER_RM = (5.0, 10.0, 20.0, 40.0)
ROUTE_OUTER_OVER_RM = 32.0
CONTROL_OUTER_OVER_RM = 64.0
ROUTE_FLOOR = 32
CONTROL_FLOOR = 64
CALIBRATION_TOLERANCE = observables.CALIBRATION_TOLERANCE
MAGNETIC_Q = 1.0
V4 = 0.0
CPU_BUDGET_SECONDS = 6.0 * 3600.0
CASE_PROCESSES = 6
THREADS_PER_CASE = 1
MEMORY_BYTES = 8 * 1024 ** 3
CHUNK_BYTES = 64 * 1024 * 1024
ADMISSION_FACTOR = 1.5
STEP_CAP_OVER_RM = 0.001
FIRST_PAIR = {"mass_over_rm": 1.04, "energy_fraction": 1.0e-3}
PILOT_RANK = 9
FIXED_CONTROL = "initial_sourced_metric_frozen"
TIME_JET_PRIMITIVE = (
    "beta_tr is the nodal radial derivative of the constraint shift rate; "
    "N_tr and second time derivatives cancel in this Riemann contraction "
    "and are not invented"
)
PRESERVED_NF256 = "lab/results/development/nsc-discovery-parent-cut-confirmation-v1"
OWNED = (
    "src/recursive_horizons/nsc_rn_observables.py",
    "scripts/derive_nsc_rn_parent.py",
    "docs/nsc-rn-parent-episode.md",
)
DEPENDENCIES = (
    "src/recursive_horizons/nsc_rn_reference.py",
    "src/recursive_horizons/nsc_rn_pg.py",
    "src/recursive_horizons/nsc_rn_source.py",
    "src/recursive_horizons/nsc_rn_observables.py",
)
PG_CALLS = (
    "build_grid",
    "sourcefree_rates",
    "metric_jets",
    "rk4_step",
)
SOURCE_CALLS = ("prepare_radial_packet",)
MEASUREMENT_CHANNELS = (
    "reconstruction_lapse",
    "reconstruction_shift",
    "reconstruction_shift_r",
    "reconstruction_shift_rr",
    "rate_max",
    "rk4_stationarity_lapse",
    "rk4_stationarity_shift",
    "charged_mass",
    "R4",
    "Ricci2",
    "trapping_minus_rn_r_plus",
)
IGNORED_ADMISSION_FLAGS = (
    "passed",
    "ok",
    "root_scientific_calibration_passed",
    "within_manufactured_tolerance",
)


class CampaignBlocker(RuntimeError):
    def __init__(self, blocker, **fields):
        super().__init__(blocker)
        self.payload = {
            "schema": SCHEMA,
            "ok": False,
            "blocker": blocker,
            "evolved": False,
            "evolution_called": False,
            "bytes_written": 0,
            "production_trajectory": False,
            "mock_physical_output": False,
            "root_review_overrides": False,
            **fields,
        }


def file_sha256(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def array_sha256(array):
    contiguous = np.ascontiguousarray(array)
    hasher = hashlib.sha256()
    hasher.update(str(contiguous.dtype).encode())
    hasher.update(str(contiguous.shape).encode())
    hasher.update(contiguous.tobytes())
    return hasher.hexdigest()


def configuration():
    return {
        "schema": SCHEMA,
        "api_version": API_VERSION,
        "chart": reference.ACTION["line_element"],
        "matter": {
            "experiment": 1,
            "kappa": reference.KAPPA,
            "multiplicity": reference.MULTIPLICITY,
            "multiplicity_counted_once": True,
            "electric_current": reference.ELECTRIC_CURRENT,
            "magnetic_q": MAGNETIC_Q,
            "V4": V4,
            "filled_sea": False,
            "magnetic_flux_fixed": True,
            "killing_frequency": "positive_exterior",
            "angular_coupling": "kappa/r",
            "four_dimensional_mass_from_kappa": False,
            "weyl_pole": False,
        },
        "cases": {
            "mass_over_rm": MASS_OVER_RM,
            "energy_fraction": ENERGY_FRACTIONS,
            "station_times_over_rm": STATION_TIMES_OVER_RM,
            "count": len(MASS_OVER_RM) * len(ENERGY_FRACTIONS),
            "r_m": "magnetic radius sqrt(P^2)",
            "initial_P2": "r_m^2",
            "r_m_is_inner_horizon": False,
        },
        "grid": {
            "route_outer_over_rm": ROUTE_OUTER_OVER_RM,
            "control_outer_over_rm": CONTROL_OUTER_OVER_RM,
            "route_point_floor": ROUTE_FLOOR,
            "control_point_floor": CONTROL_FLOOR,
            "spacing": observables.FORMULAS["spacing"],
            "throat": "(r_plus - r_minus)/16",
            "derivative_ladder": False,
            "cfl": observables.FORMULAS["full_cfl"],
        },
        "calibration_tolerance": CALIBRATION_TOLERANCE,
        "resources": {
            "cpu_budget_seconds": CPU_BUDGET_SECONDS,
            "case_processes": CASE_PROCESSES,
            "threads_per_case": THREADS_PER_CASE,
            "memory_bytes": MEMORY_BYTES,
            "chunk_bytes": CHUNK_BYTES,
            "admission_factor": ADMISSION_FACTOR,
        },
        "phases": {
            "default": "preview",
            "read_only_default": True,
            "calibration": "sourcefree PG reconstruction, rates, jets, one RK4 step",
            "physical_run_this_followup": False,
        },
        "admission": "measured residuals at 1e-4 plus source binding",
        "root_review_overrides": False,
        "writer": "creation_only_after_frozen_root_commit",
        "hard_acceptance_gate": False,
        "periodic_runs": False,
    }


def _load_module(name):
    try:
        return importlib.import_module(f"recursive_horizons.{name}")
    except ImportError:
        return None


def _missing_callables(module, names, owner):
    if module is None:
        return [f"{owner} (module missing)"]
    return [
        f"{owner}.{name}" for name in names
        if not callable(getattr(module, name, None))
    ]


def adapter_gap():
    pg = _load_module("nsc_rn_pg")
    source = _load_module("nsc_rn_source")
    return _missing_callables(pg, PG_CALLS, "nsc_rn_pg") + _missing_callables(
        source, SOURCE_CALLS, "nsc_rn_source",
    )


def magnetic_reference(mass_over_rm, r_m=1.0):
    """RNReference with P² = r_m² and fixed magnetic q = 1. r_m is not r₋."""
    ratio = float(mass_over_rm)
    r_m = float(r_m)
    if not math.isfinite(ratio) or ratio <= 1.0 or r_m <= 0.0:
        raise CampaignBlocker("M/r_m must exceed 1 at a positive magnetic radius")
    model = observables.reference_model(ratio * r_m, r_m)
    if abs(model.magnetic_r2 - r_m * r_m) > 1e-9 * max(1.0, r_m * r_m):
        raise CampaignBlocker("reference magnetic radius did not realize P^2 = r_m^2")
    if not model.horizons["black_hole"]:
        raise CampaignBlocker("magnetic family has no horizon")
    return model


def _member_grid(model, r_m, outer_over_rm, floor):
    horizons = model.horizons
    r_minus = float(horizons["r_minus"])
    r_plus = float(horizons["r_plus"])
    gap = r_plus - r_minus
    excision = 0.5 * (r_minus + r_plus)
    if not bool(np.asarray(model.pure_outflow(excision))):
        raise CampaignBlocker("midpoint excision is not pure outflow")
    outer = float(outer_over_rm) * float(r_m)
    child = (r_plus + gap / 16.0, r_plus + 0.5 * gap)
    parent = (child[1], outer)
    placed = reference.placement(
        model, parent=parent, child=child, excision=excision, outer="absorbing",
    )
    if not placed["supported"]:
        raise CampaignBlocker("reference placement refused the fixed cuts: " + "; ".join(placed["reasons"]))
    scales = _load_module("nsc_rn_source").packet_resolution_scales(model, rank=PILOT_RANK)
    spacing = observables.throat_spacing(
        gap, lambda_min=scales["lambda_min"], sigma=scales["sigma"],
    )
    points = observables.route_points(outer - excision, spacing["spacing"], floor)
    return {
        "r_minus": r_minus,
        "r_plus": r_plus,
        "horizon_gap": gap,
        "excision": excision,
        "child": child,
        "parent": parent,
        "outer": outer,
        "same_areal_radius": False,
        "r_m_identified_with_r_minus": False,
        "pure_outflow": True,
        "placement_supported": True,
        "spacing": spacing["spacing"],
        "pending_scales": spacing["pending_scales"],
        "lambda_min": scales["lambda_min"],
        "sigma": scales["sigma"],
        "packet_center": scales["center"],
        "packet_width": scales["width"],
        "e_max": scales["e_max"],
        "points": points,
        "point_floor": int(floor),
        "derivative_ladder": False,
    }


def case_geometry(mass_over_rm, r_m=1.0):
    model = magnetic_reference(mass_over_rm, r_m)
    route = _member_grid(model, r_m, ROUTE_OUTER_OVER_RM, ROUTE_FLOOR)
    control = _member_grid(model, r_m, CONTROL_OUTER_OVER_RM, CONTROL_FLOOR)
    if route["child"] != control["child"] or route["excision"] != control["excision"]:
        raise CampaignBlocker("child and excision cuts changed between route and control")
    return {
        "mass_over_rm": float(mass_over_rm),
        "r_m": float(r_m),
        "mass": float(model.mass),
        "magnetic_r2": float(model.magnetic_r2),
        "magnetic_q": MAGNETIC_Q,
        "V4": V4,
        "filled_sea": False,
        "A": float(model.A),
        "C_F": float(model.C_F),
        "flux": float(model.flux),
        "r_minus": route["r_minus"],
        "r_plus": route["r_plus"],
        "horizon_gap": route["horizon_gap"],
        "r_m_identified_with_r_minus": False,
        "excision": route["excision"],
        "child": route["child"],
        "route": route,
        "control": control,
        "observers": observables.OBSERVERS,
        "observers_declared_once": True,
        "model": model,
    }


def case_table(r_m=1.0):
    rows = []
    for ratio in MASS_OVER_RM:
        geometry = case_geometry(ratio, r_m)
        geometry.pop("model")
        for fraction in ENERGY_FRACTIONS:
            rows.append({
                **geometry,
                "energy_fraction": float(fraction),
                "station_times_over_rm": STATION_TIMES_OVER_RM,
                "launched": False,
            })
    if len(rows) != 6:
        raise CampaignBlocker("production case table is not the six declared cases")
    return rows


def producers():
    pins = {}
    for relative in OWNED + DEPENDENCIES:
        path = LAB / relative
        pins[str(Path("lab") / relative)] = file_sha256(path) if path.is_file() else None
    return pins


def dependency_status():
    status = {}
    for name in ("nsc_rn_reference", "nsc_rn_pg", "nsc_rn_source"):
        status[name] = "present" if _load_module(name) is not None else "pending"
    status["nsc_rn_observables"] = "present"
    status["adapter_gap"] = adapter_gap()
    return status


def _git_blob_sha(commit, relative):
    if len(commit) != 40 or any(character not in "0123456789abcdef" for character in commit):
        raise CampaignBlocker("frozen root commit must be a full hexadecimal SHA")
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts:
        raise CampaignBlocker("producer path must stay inside the repository")
    result = subprocess.run(
        ["git", "cat-file", "blob", f"{commit}:{path.as_posix()}"],
        cwd=ROOT, capture_output=True, check=True, timeout=5,
    )
    return hashlib.sha256(result.stdout).hexdigest()


def frozen_blob_match(commit, pins):
    pending = [path for path, digest in pins.items() if digest is None]
    if pending:
        raise CampaignBlocker("pending dependency: " + ", ".join(pending), mode="writer")
    try:
        matched = all(_git_blob_sha(commit, path) == digest for path, digest in pins.items())
    except CampaignBlocker:
        raise
    except (OSError, subprocess.SubprocessError):
        return False
    return matched


def preserved_nf256():
    return {
        "path": PRESERVED_NF256,
        "mode": "read_only",
        "opened": False,
        "new_periodic_run": False,
    }


def preview():
    return {
        "schema": SCHEMA,
        "mode": "preview",
        "pure": True,
        "read_only": True,
        "configuration": configuration(),
        "cases": case_table(),
        "dependencies": dependency_status(),
        "producers": producers(),
        "preserved_nf256": preserved_nf256(),
        "evolved": False,
        "evolution_called": False,
        "bytes_written": 0,
        "production_trajectory": False,
        "processes_launched": 0,
        "mock_physical_output": False,
        "root_review_overrides": False,
    }


def _max_abs(left, right):
    return float(np.max(np.abs(np.asarray(left, dtype=float) - np.asarray(right, dtype=float))))


def _rate_max(rates):
    def magnitude(item):
        array = np.asarray(item)
        if array.size == 0:
            return None
        return np.abs(array)

    if isinstance(rates, dict):
        values = [magnitude(item) for item in rates.values()]
        values = [item for item in values if item is not None]
        if not values:
            raise CampaignBlocker("sourcefree_rates returned no numeric samples")
        return float(max(np.max(item) for item in values))
    array = magnitude(rates)
    if array is None or not np.all(np.isfinite(array)):
        raise CampaignBlocker("sourcefree_rates returned an empty sample")
    return float(np.max(array))


def _reconstruct_sourcefree(pg, geometry):
    """Owned PG reconstruction. The wrapper's reference call uses a different signature."""
    route = geometry["route"]
    mass = geometry["mass"]
    charge = geometry["r_m"]
    points = int(route["points"])
    outer = route["outer"]
    try:
        return pg.sourcefree_rates(
            mass, r_m=charge, points=points, r_out=outer, killing_frequency=1.0,
        )[:3]
    except (TypeError, AttributeError) as error:
        raise CampaignBlocker(
            "adapter gap: nsc_rn_pg.sourcefree_rates is not callable on the current reference: " + str(error),
            adapter_gap=["nsc_rn_pg.sourcefree_rates", str(error)],
            mock_physical_output=False,
        ) from error


def _measure_case(pg, geometry):
    """One source-free PG reconstruction and one RK4 step against the reference."""
    model = geometry["model"]
    try:
        grid, state, rates = _reconstruct_sourcefree(pg, geometry)
    except (TypeError, pg.ChartAdmissionError, observables.ObservableError) as error:
        raise CampaignBlocker(
            "source-free PG reconstruction was not admitted: " + str(error),
            adapter_gap=["nsc_rn_pg.sourcefree_rates"],
            mock_physical_output=False,
        ) from error
    radii = np.asarray(grid["radius"], dtype=float)
    analytic_sample = reference.sourcefree_static_rn(
        radii, model.mass, A=model.A, C_F=model.C_F, flux=model.flux,
    )
    analytic = analytic_sample["metric_jets"]
    jets = rates["jets"]
    artifact = observables.endpoint_second_derivative_bias(jets["shift_rr"], analytic["beta_rr"])
    dt = 0.5 * float(rates["cfl_dt"])
    try:
        advanced = pg.rk4_step(state, grid, dt)
        advanced_rates = pg.stage_rates(advanced, grid)
    except TypeError as error:
        raise CampaignBlocker(
            "adapter gap: rk4_step signature mismatch: " + str(error),
            adapter_gap=["nsc_rn_pg.rk4_step"],
            mock_physical_output=False,
        ) from error
    located = observables.trapping_horizon(radii, rates["lapse"], rates["shift"])
    trapping_gap = None if located["trapping_horizon"] is None else abs(
        located["trapping_horizon"] - float(model.horizons["r_plus"])
    )
    curvature = observables.curvature_from_pg_jets(
        radii, rates["lapse"], rates["shift"],
        lapse_r=jets["lapse_r"], shift_r=jets["shift_r"],
        lapse_rr=jets["lapse_rr"], shift_rr=jets["shift_rr"],
        lapse_t=rates["lapse_rate"], shift_t=rates["shift_rate"],
        shift_tr=rates["time_jets"]["shift_tr"],
    )
    analytic_invariants = analytic_sample["invariants"]
    if curvature["status"] != "stationary_radial_jets":
        r4_gap = None
        ricci_gap = None
    else:
        r4_gap = _max_abs(curvature["R4"], analytic_invariants["R4"])
        ricci_gap = _max_abs(curvature["Ricci2"], analytic_invariants["Ricci2"])
    rate_values = {
        "phi_rate": rates["phi_rate"],
        "lapse_rate": rates["lapse_rate"],
        "shift_rate": rates["shift_rate"],
        "mass_rate": rates["mass_rate"],
    }
    residuals = {
        "reconstruction_lapse": _max_abs(rates["lapse"], analytic["N"]),
        "reconstruction_shift": _max_abs(rates["shift"], analytic["beta"]),
        "reconstruction_shift_r": _max_abs(jets["shift_r"], analytic["beta_r"]),
        "reconstruction_shift_rr": _max_abs(jets["shift_rr"], analytic["beta_rr"]),
        "rate_max": _rate_max(rate_values),
        "rk4_stationarity_lapse": _max_abs(advanced_rates["lapse"], analytic["N"]),
        "rk4_stationarity_shift": _max_abs(advanced_rates["shift"], analytic["beta"]),
        "charged_mass": _max_abs(rates["mass"], analytic_sample["charged_mass"]),
        "R4": r4_gap,
        "Ricci2": ricci_gap,
        "trapping_minus_rn_r_plus": trapping_gap,
    }
    return {
        "mass_over_rm": geometry["mass_over_rm"],
        "r_m": geometry["r_m"],
        "r_minus": geometry["r_minus"],
        "r_m_identified_with_r_minus": False,
        "points": int(radii.size),
        "dt": dt,
        "endpoint_second_derivative": artifact["artifact"],
        "measured_residuals": residuals,
        "lapse_sha256": array_sha256(rates["lapse"]),
        "shift_sha256": array_sha256(rates["shift"]),
        "injected_analytic_jets": False,
        "measurement_origin": "nsc_rn_pg.sourcefree_rates+rk4_step",
    }


def _worst(rows):
    worst = {}
    for channel in MEASUREMENT_CHANNELS:
        values = [row["measured_residuals"][channel] for row in rows]
        if any(value is None for value in values):
            worst[channel] = None
        else:
            worst[channel] = max(values)
    return worst


def calibrate(*, r_m=1.0):
    """Source-free numerical calibration. Analytic jets alone do not pass."""
    missing = _missing_callables(_load_module("nsc_rn_pg"), PG_CALLS, "nsc_rn_pg")
    families = [case_geometry(ratio, r_m) for ratio in MASS_OVER_RM]
    public = [{key: value for key, value in row.items() if key != "model"} for row in families]
    if missing:
        return {
            "schema": SCHEMA,
            "mode": "calibration",
            "ok": False,
            "passed": False,
            "blocker": "adapter gap: " + ", ".join(missing),
            "adapter_gap": missing,
            "calibration_executed": False,
            "injected_analytic_jets": False,
            "measured_residuals": None,
            "cases": public,
            "tolerance": CALIBRATION_TOLERANCE,
            "source_binding": producers(),
            "evolved": False,
            "evolution_called": False,
            "bytes_written": 0,
            "production_trajectory": False,
            "mock_physical_output": False,
            "root_review_overrides": False,
        }
    pg = _load_module("nsc_rn_pg")
    try:
        measured = [_measure_case(pg, row) for row in families]
    except CampaignBlocker as blocked:
        report = dict(blocked.payload)
        report.update({
            "mode": "calibration",
            "calibration_executed": False,
            "injected_analytic_jets": False,
            "measured_residuals": None,
            "cases": public,
            "tolerance": CALIBRATION_TOLERANCE,
            "source_binding": producers(),
        })
        return report
    residuals = _worst(measured)
    residual_ok = all(
        value is not None and value <= CALIBRATION_TOLERANCE for value in residuals.values()
    )
    endpoint_ok = all(row["endpoint_second_derivative"] is None for row in measured)
    passed = residual_ok and endpoint_ok
    reasons = []
    if not residual_ok:
        reasons.append("measured residual exceeds 1e-4")
    if not endpoint_ok:
        reasons.append("naive D@D endpoint jet is half the consistent second derivative")
    payload = json.dumps(
        [{"residuals": row["measured_residuals"], "lapse": row["lapse_sha256"], "shift": row["shift_sha256"]}
         for row in measured],
        sort_keys=True,
    ).encode()
    return {
        "schema": SCHEMA,
        "mode": "calibration",
        "ok": passed,
        "passed": passed,
        "blocker": None if passed else "; ".join(reasons),
        "adapter_gap": [],
        "calibration_executed": True,
        "injected_analytic_jets": False,
        "measured_residuals": residuals,
        "measurement_sha256": hashlib.sha256(payload).hexdigest(),
        "measurement_origin": "nsc_rn_pg.sourcefree_rates+rk4_step",
        "cases": measured,
        "tolerance": CALIBRATION_TOLERANCE,
        "source_binding": producers(),
        "evolved": True,
        "evolution_called": True,
        "rk4_steps_per_case": 1,
        "production_trajectory": False,
        "bytes_written": 0,
        "mock_physical_output": False,
        "root_review_overrides": False,
        "ignored_flags": IGNORED_ADMISSION_FLAGS,
    }


def admission_decision(record, pins=None):
    """Residuals and the producer binding decide. Root-review booleans do not."""
    pins = producers() if pins is None else pins
    reasons = []
    if not isinstance(record, dict) or record.get("schema") != SCHEMA or record.get("mode") != "calibration":
        reasons.append("calibration record is missing or has the wrong schema")
        record = record if isinstance(record, dict) else {}
    if record.get("adapter_gap"):
        reasons.append("adapter gap")
    if record.get("endpoint_second_derivative") == "D_at_D_half_actual":
        reasons.append("naive D@D endpoint jet")
    if record.get("calibration_executed") is not True:
        reasons.append("calibration was not executed")
    if not record.get("measurement_sha256") or record.get("measurement_origin") != "nsc_rn_pg.sourcefree_rates+rk4_step":
        reasons.append("immutable measurement payload is missing")
    if any(value is None for value in pins.values()):
        reasons.append("pending dependency")
    if record.get("source_binding") != pins:
        reasons.append("source binding does not match the current producers")
    residuals = record.get("measured_residuals")
    if not isinstance(residuals, dict):
        reasons.append("measured residuals are missing")
    else:
        for channel in MEASUREMENT_CHANNELS:
            if channel not in residuals or residuals[channel] is None:
                reasons.append("missing measurement " + channel)
                continue
            try:
                value = abs(float(residuals[channel]))
            except (TypeError, ValueError):
                reasons.append("nonnumeric measurement " + channel)
                continue
            if value > CALIBRATION_TOLERANCE:
                reasons.append(channel + " exceeds 1e-4")
    return {
        "admitted": not reasons,
        "blocker": None if not reasons else "; ".join(reasons),
        "tolerance": CALIBRATION_TOLERANCE,
        "root_review_overrides": False,
        "boolean_pass_ignored": True,
        "ignored_flags": IGNORED_ADMISSION_FLAGS,
    }


def admit(step_seconds, steps, *, remaining, factor=ADMISSION_FACTOR):
    step_seconds = float(step_seconds)
    steps = int(steps)
    remaining = float(remaining)
    factor = float(factor)
    if min(step_seconds, remaining, factor) < 0.0 or steps < 0:
        raise CampaignBlocker("admission timing must be nonnegative")
    forecast = factor * step_seconds * steps
    return {
        "forecast_cpu_seconds": forecast,
        "remaining_cpu_seconds": remaining,
        "admission_factor": factor,
        "admitted": bool(forecast <= remaining or math.isclose(forecast, remaining, rel_tol=0.0, abs_tol=1e-12)),
    }


def _sealed(path):
    resolved = Path(path).expanduser().resolve()
    for prefix in (LAB / "results", ROOT / "results"):
        prefix = prefix.resolve()
        if resolved == prefix or prefix in resolved.parents:
            return True
    return False


def write_stage(directory, stage, record, arrays=None, *, blob_match):
    """Creation-only checkpoint. Refuses sealed results and unmatched commits."""
    if blob_match is not True:
        raise CampaignBlocker("creation requires a frozen root commit", mode=stage)
    if stage not in ("prepare", "dependency") and not str(stage).startswith("station-"):
        raise CampaignBlocker("unknown checkpoint stage", mode=stage)
    directory = Path(directory).expanduser().resolve()
    if _sealed(directory):
        raise CampaignBlocker("existing outputs are sealed", mode=stage)
    json_path = directory / f"{stage}.json"
    npz_path = directory / f"{stage}.npz"
    if json_path.exists() or npz_path.exists():
        raise CampaignBlocker("checkpoint stage already exists", mode=stage)
    payload = {} if arrays is None else {key: np.asarray(value) for key, value in arrays.items()}
    blob = _npz_bytes(payload)
    body = dict(record)
    body.update({
        "schema": SCHEMA,
        "mode": stage,
        "creation_only": True,
        "evolved": False,
        "production_trajectory": False,
        "payload_sha256": hashlib.sha256(blob).hexdigest(),
        "array_sha256": {key: array_sha256(value) for key, value in payload.items()},
        "payload_bytes": len(blob),
    })
    encoded = (json.dumps(body, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    if len(blob) > CHUNK_BYTES or len(blob) + len(encoded) > CHUNK_BYTES:
        raise CampaignBlocker("checkpoint exceeds 64 MiB", mode=stage)
    directory.mkdir(parents=True, exist_ok=True)
    for path, data in ((npz_path, blob), (json_path, encoded)):
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o444)
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data)
    body["bytes_written"] = len(blob) + len(encoded)
    return body


def _npz_bytes(payload):
    import io
    stream = io.BytesIO()
    np.savez_compressed(stream, **payload)
    return stream.getvalue()


def prepare(directory, *, producer_commit, calibration_record, inputs=()):
    """Creation-only plan after a measured calibration binding. Does not evolve."""
    decision = admission_decision(calibration_record)
    if not decision["admitted"]:
        raise CampaignBlocker(decision["blocker"], mode="prepare", admission=decision)
    if not producer_commit:
        raise CampaignBlocker("creation requires a frozen root commit", mode="prepare")
    pins = producers()
    if not frozen_blob_match(producer_commit, pins):
        raise CampaignBlocker("producer bytes are not the frozen root commit", mode="prepare")
    record = {
        "producing_commit": producer_commit,
        "producers": pins,
        "inputs": bind_inputs(inputs) if inputs else {},
        "source_binding": calibration_record.get("source_binding"),
        "measured_residuals": calibration_record.get("measured_residuals"),
        "cases": case_table(),
        "neutral_initial_data": "not_launched",
        "evolution_called": False,
        "physical_run_this_followup": False,
        "root_review_overrides": False,
    }
    return write_stage(directory, "prepare", record, {}, blob_match=True)


def matched_step(cfl_dt, r_m):
    """March step. The cap is 0.001 r_m. Half the full CFL is used only when it is smaller."""
    cap = float(STEP_CAP_OVER_RM) * float(r_m)
    cfl_half = 0.5 * float(cfl_dt)
    step = min(cfl_half, cap)
    return {
        "dt": step,
        "cfl_dt": float(cfl_dt),
        "cfl_half_step": cfl_half,
        "step_cap": cap,
        "step_cap_over_rm": STEP_CAP_OVER_RM,
        "limited_by": "step_cap" if cap < cfl_half else "half_cfl",
    }


def freeze_sourced_metric(rates):
    """Hold the initial sourced lapse and shift. This is not a second coupled evolution."""
    return {
        "choice": FIXED_CONTROL,
        "reconstructs_metric": False,
        "lapse": np.array(rates["lapse"], dtype=float, copy=True),
        "shift": np.array(rates["shift"], dtype=float, copy=True),
    }


def _clone_state(state):
    import copy
    return copy.deepcopy(state)


def fixed_rk4(state, grid, metric, dt):
    """Dirac RK4 on a frozen lapse and shift. stage_rates is not called."""
    pg = _load_module("nsc_rn_pg")
    observables_mod = _load_module("nsc_rn_observables")
    step = float(dt)
    cfl = observables_mod.full_cfl_dt(metric["lapse"], metric["shift"], grid["spacing"])
    if step > cfl * (1.0 + 1e-12):
        raise CampaignBlocker("fixed Dirac step exceeds the frozen-metric CFL", mode="march")
    lapse = metric["lapse"]
    shift = metric["shift"]
    kappa = state["kappa"]

    def sample(phi):
        phi_rate, sat = pg.dirac_rates(grid, phi, lapse, shift, kappa=kappa)
        boundary = pg._probability_boundary(phi, lapse, shift)
        return phi_rate, float(sat), boundary

    phi0 = state["phi"]
    k1, s1, b1 = sample(phi0)
    k2, s2, b2 = sample(phi0 + 0.5 * step * k1)
    k3, s3, b3 = sample(phi0 + 0.5 * step * k2)
    k4, s4, b4 = sample(phi0 + step * k3)
    phi = phi0 + (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    weights = (1.0, 2.0, 2.0, 1.0)
    sat = sum(weight * value for weight, value in zip(weights, (s1, s2, s3, s4))) / 6.0
    excision = sum(weight * item["excision_outflow"] for weight, item in zip(weights, (b1, b2, b3, b4))) / 6.0
    outer = sum(weight * item["outer_outflow"] for weight, item in zip(weights, (b1, b2, b3, b4))) / 6.0
    metric_tt = float(lapse[-1] ** 2 - shift[-1] ** 2)
    killing_rate = math.sqrt(metric_tt) if metric_tt > 0.0 else 0.0
    updated = _clone_state(state)
    updated["phi"] = phi
    updated["inner_mass"] = float(state["inner_mass"])
    updated["clocks"] = {
        "t": state["clocks"]["t"] + step,
        "killing": state["clocks"]["killing"] + step * killing_rate,
        "normal": state["clocks"]["normal"] + step * float(lapse[0]),
    }
    updated["stocks"] = {
        "excision_flux": float(state["stocks"]["excision_flux"]),
        "outer_flux": float(state["stocks"]["outer_flux"]),
        "sat_debit": float(state["stocks"]["sat_debit"]) + step * sat,
    }
    updated["probability_stocks"] = {
        "remaining": pg._probability_norm(phi, grid["weights"]),
        "excision_outflow": float(state["probability_stocks"]["excision_outflow"]) + step * excision,
        "outer_outflow": float(state["probability_stocks"]["outer_outflow"]) + step * outer,
    }
    return updated


def state_arrays(coupled, fixed, metric):
    """Scalar stocks and both spinors. The derivative matrix is not stored."""
    return {
        "coupled_phi": np.asarray(coupled["phi"]),
        "fixed_phi": np.asarray(fixed["phi"]),
        "occupations": np.asarray(coupled["occupations"], dtype=float),
        "coupled_inner_mass": np.array([coupled["inner_mass"]], dtype=float),
        "fixed_inner_mass": np.array([fixed["inner_mass"]], dtype=float),
        "coupled_t": np.array([coupled["clocks"]["t"]], dtype=float),
        "coupled_killing": np.array([coupled["clocks"]["killing"]], dtype=float),
        "coupled_normal": np.array([coupled["clocks"]["normal"]], dtype=float),
        "fixed_t": np.array([fixed["clocks"]["t"]], dtype=float),
        "fixed_killing": np.array([fixed["clocks"]["killing"]], dtype=float),
        "fixed_normal": np.array([fixed["clocks"]["normal"]], dtype=float),
        "coupled_excision_flux": np.array([coupled["stocks"]["excision_flux"]], dtype=float),
        "coupled_outer_flux": np.array([coupled["stocks"]["outer_flux"]], dtype=float),
        "coupled_sat": np.array([coupled["stocks"]["sat_debit"]], dtype=float),
        "fixed_excision_flux": np.array([fixed["stocks"]["excision_flux"]], dtype=float),
        "fixed_outer_flux": np.array([fixed["stocks"]["outer_flux"]], dtype=float),
        "fixed_sat": np.array([fixed["stocks"]["sat_debit"]], dtype=float),
        "coupled_prob_remaining": np.array([coupled["probability_stocks"]["remaining"]], dtype=float),
        "coupled_prob_excision": np.array([coupled["probability_stocks"]["excision_outflow"]], dtype=float),
        "coupled_prob_outer": np.array([coupled["probability_stocks"]["outer_outflow"]], dtype=float),
        "fixed_prob_remaining": np.array([fixed["probability_stocks"]["remaining"]], dtype=float),
        "fixed_prob_excision": np.array([fixed["probability_stocks"]["excision_outflow"]], dtype=float),
        "fixed_prob_outer": np.array([fixed["probability_stocks"]["outer_outflow"]], dtype=float),
        "frozen_lapse": np.asarray(metric["lapse"], dtype=float),
        "frozen_shift": np.asarray(metric["shift"], dtype=float),
    }


def apply_checkpoint(state, arrays, arm):
    """Replace spinor, clocks, stocks and mass. Occupations stay the stored weights."""
    updated = _clone_state(state)
    prefix = "coupled" if arm == "coupled" else "fixed"
    updated["phi"] = np.array(arrays[prefix + "_phi"], copy=True)
    updated["occupations"] = np.array(arrays["occupations"], dtype=float, copy=True)
    updated["inner_mass"] = float(arrays[prefix + "_inner_mass"][0])
    updated["clocks"] = {
        "t": float(arrays[prefix + "_t"][0]),
        "killing": float(arrays[prefix + "_killing"][0]),
        "normal": float(arrays[prefix + "_normal"][0]),
    }
    updated["stocks"] = {
        "excision_flux": float(arrays[prefix + "_excision_flux"][0]),
        "outer_flux": float(arrays[prefix + "_outer_flux"][0]),
        "sat_debit": float(arrays[prefix + "_sat"][0]),
    }
    updated["probability_stocks"] = {
        "remaining": float(arrays[prefix + "_prob_remaining"][0]),
        "excision_outflow": float(arrays[prefix + "_prob_excision"][0]),
        "outer_outflow": float(arrays[prefix + "_prob_outer"][0]),
    }
    return updated


def admit_prepared_packet(packet):
    """Inward current, moments, energy and rank-one CAR of one prepared column."""
    car = np.asarray(packet["car"]["physical"])
    hermitian = 0.5 * (car + car.conj().T)
    eigenvalues = np.linalg.eigvalsh(hermitian)
    moments = packet["radial_moments"]
    inward = packet["inward_flux"]
    return {
        "inward_at_center": bool(inward["inward_at_center"]),
        "center_probability_current": float(inward["center_probability_current"]),
        "integrated_negative_current": float(inward["integrated_negative_current"]),
        "mean_radius": float(moments["center"]),
        "width_rms": float(moments["width_rms"]),
        "requested_center": float(moments["requested_center"]),
        "requested_width": float(moments["requested_width"]),
        "packet_mode_energy": float(packet["packet_mode_energy"]),
        "killing_energy": float(packet["killing_energy"]),
        "filling": float(packet["packet_filling"]),
        "car_rank": int(np.count_nonzero(np.abs(eigenvalues) > 1.0e-8)),
        "phi_shape": [int(item) for item in np.asarray(packet["phi"]).shape],
        "requested_frequencies": [float(item) for item in packet["spectrum"]["requested"]],
    }


def sat_norm_ledger(raw_sat, probability, initial_norm):
    """Raw SAT stock versus the probability loss.

    The historical field stocks.sat_debit stores -a_minus |chi_minus|^2.
    The SBP identity 2 Re(chi, rate) makes the norm loss twice that coefficient.
    The raw number is not copied into a second energy or into the mass flux.
    """
    raw = float(raw_sat)
    actual = 2.0 * raw
    remaining = float(probability["remaining"])
    excision = float(probability["excision_outflow"])
    outer = float(probability["outer_outflow"])
    return {
        "historical_field": "stocks.sat_debit",
        "raw_coefficient": "-a_minus |chi_minus|^2",
        "raw_sat_coefficient": raw,
        "norm_factor": 2.0,
        "actual_probability_loss": actual,
        "probability_closure_residual": remaining + excision + outer + actual - float(initial_norm),
        "coefficient_duplicated": False,
        "included_in_mass_flux": False,
    }


def station_observation(state, grid, rates, packet=None, initial_norm=None):
    """Stocks, trapping radius and charge-corrected mass. Sourced curvature stays null."""
    observables_mod = _load_module("nsc_rn_observables")
    curvature = observables_mod.curvature_from_pg_jets(
        grid["radius"], rates["lapse"], rates["shift"],
        lapse_r=rates["jets"]["lapse_r"], shift_r=rates["jets"]["shift_r"],
        lapse_rr=rates["jets"]["lapse_rr"], shift_rr=rates["jets"]["shift_rr"],
        lapse_t=rates["lapse_rate"], shift_t=rates["shift_rate"],
        shift_tr=rates["time_jets"]["shift_tr"],
    )
    located = observables_mod.trapping_horizon(grid["radius"], rates["lapse"], rates["shift"])
    charged = observables_mod.charged_mass(
        grid["radius"], rates["lapse"], rates["shift"], grid["magnetic_r2"],
    )
    unresolved = curvature["R4"] is None
    if initial_norm is None:
        initial_norm = state["probability_stocks"]["remaining"]
    report = {
        "probability_stocks": dict(state["probability_stocks"]),
        "mass_flux_stocks": {
            "excision_flux": state["stocks"]["excision_flux"],
            "outer_flux": state["stocks"]["outer_flux"],
        },
        "sat_norm_stock": sat_norm_ledger(
            state["stocks"]["sat_debit"], state["probability_stocks"], initial_norm,
        ),
        "stocks_are_separate": True,
        "trapping_horizon": located["trapping_horizon"],
        "charged_mass_inner": float(np.asarray(charged).reshape(-1)[0]),
        "charged_mass_outer": float(np.asarray(charged).reshape(-1)[-1]),
        "curvature_status": curvature["status"],
        "R4": None if unresolved else curvature["R4"],
        "Ricci2": None if unresolved else curvature["Ricci2"],
        "required_time_jet_primitive": TIME_JET_PRIMITIVE if unresolved else None,
        "analytic_curvature_injected": False,
        "source_force_each_rk_stage": True,
        "clock_t": state["clocks"]["t"],
    }
    source = _load_module("nsc_rn_source")
    account = source.energy_account(
        state["phi"], state["occupations"], rates["lapse"], rates["shift"],
        grid["radius"], grid["weights"], grid["derivative"],
        phi_rate=rates["phi_rate"], lapse_t=rates["lapse_rate"], shift_t=rates["shift_rate"],
        lapse_r=rates["jets"]["lapse_r"], shift_r=rates["jets"]["shift_r"],
    )
    mass_from_energy = float(grid["G_N"]) * account["e_H_density"] / rates["lapse"]
    report["energy"] = {
        "E_coordinate": account["E_coordinate"],
        "E_normal": account["E_normal"],
        "metric_power": account["metric_power"],
        "boundary_energy_flux": account["boundary_energy_flux"],
        "boundary_re_A_H": account["boundary_re_A_H"],
        "sat_energy_rate": account["sat_energy_rate"],
        "phi_motion_rate": account["phi_motion_rate"],
        "predicted_coordinate_rate": account["predicted_coordinate_rate"],
        "probability_sat_raw": account["probability_sat_raw"],
        "probability_sat_is_energy": False,
        "shell_occupation_factor": account["shell_occupation_factor"],
        "multiplicity_applied_once": True,
        "normal_continuity_residual": account["normal_continuity_residual"],
        "mass_slope_residual": float(np.max(np.abs(rates["mass_slope"] - mass_from_energy))),
        "field_energy_not_added_to_mass": True,
    }
    if packet is not None:
        report["packet"] = admit_prepared_packet(packet)
    return report


def march_pair(coupled, fixed, grid, metric, *, dt, steps=None, station_time=None,
               cpu_seconds=None):
    """Advance the coupled RK4 and the frozen-metric Dirac RK4 together.

    The coupled arm calls rk4_step, which rebuilds the source at every stage.
    The fixed arm does not. cpu_seconds stops the loop before the next step
    when the factor-1.5 projection would pass the cap. The first step is timed
    and is not refused in advance.
    """
    import time
    pg = _load_module("nsc_rn_pg")
    spent = 0.0
    taken = 0
    last = None
    target = None if station_time is None else float(station_time)
    while True:
        if steps is not None and taken >= int(steps):
            break
        if target is not None and coupled["clocks"]["t"] >= target - 1e-15:
            break
        if last is not None and cpu_seconds is not None and ADMISSION_FACTOR * (spent + last) > float(cpu_seconds):
            return {
                "coupled": coupled, "fixed": fixed, "steps": taken, "cpu_seconds": spent,
                "stopped": "cpu_budget", "dt": float(dt),
            }
        current = float(dt)
        if target is not None:
            current = min(current, target - coupled["clocks"]["t"])
        if current <= 0.0:
            break
        started = time.perf_counter()
        coupled = pg.rk4_step(coupled, grid, current)
        fixed = fixed_rk4(fixed, grid, metric, current)
        last = time.perf_counter() - started
        spent += last
        taken += 1
    return {
        "coupled": coupled, "fixed": fixed, "steps": taken, "cpu_seconds": spent,
        "stopped": "station" if target is not None else "steps", "dt": float(dt),
    }


def run(directory, *, producer_commit, calibration_record, probe_step_seconds=None,
        mass_over_rm=None, energy_fraction=None, station_limit=None, cpu_seconds=None):
    """Matched-pair march after the measurement binding and the root freeze."""
    missing = adapter_gap()
    decision = admission_decision(calibration_record)
    if missing:
        raise CampaignBlocker(
            "adapter gap: " + ", ".join(missing),
            mode="run",
            adapter_gap=missing,
            admission=decision,
            evolution_called=False,
            mock_physical_output=False,
        )
    if not decision["admitted"]:
        raise CampaignBlocker(decision["blocker"], mode="run", admission=decision, evolution_called=False)
    ratio = FIRST_PAIR["mass_over_rm"] if mass_over_rm is None else float(mass_over_rm)
    fraction = FIRST_PAIR["energy_fraction"] if energy_fraction is None else float(energy_fraction)
    limit = 5.0 if station_limit is None else float(station_limit)
    budget = CPU_BUDGET_SECONDS if cpu_seconds is None else float(cpu_seconds)
    if ratio not in MASS_OVER_RM or fraction not in ENERGY_FRACTIONS:
        raise CampaignBlocker("case is outside the magnetic family", mode="run", evolution_called=False)
    if limit not in STATION_TIMES_OVER_RM and limit != 5.0:
        raise CampaignBlocker("station limit is not 5, 10, 20, or 40", mode="run", evolution_called=False)
    if not producer_commit or not frozen_blob_match(producer_commit, producers()):
        raise CampaignBlocker("creation requires a frozen root commit", mode="run", evolution_called=False)
    if directory is None:
        raise CampaignBlocker("run needs an existing prepare checkpoint", mode="run", evolution_called=False)
    pg = _load_module("nsc_rn_pg")
    source = _load_module("nsc_rn_source")
    geometry = case_geometry(ratio)
    try:
        grid = pg.build_grid(
            geometry["mass"], r_m=geometry["r_m"], points=geometry["route"]["points"],
            r_out=geometry["route"]["outer"],
        )
        packet = source.prepare_radial_packet(
            mass_over_rm=ratio, energy_target=fraction, rank=PILOT_RANK, model=grid["model"],
            radii=grid["radius"], weights=grid["weights"], derivative=grid["derivative"],
        )
        phi = np.asarray(packet["phi"])
        if phi.shape != (2, grid["points"], 1) or not np.array_equal(
            np.asarray(packet["background"]["grid"]), grid["radius"],
        ):
            raise CampaignBlocker(
                "adapter gap: source packet is not the rank-1 column on the PG grid",
                mode="run", evolution_called=False, mock_physical_output=False,
            )
        coupled = pg.make_state(
            grid, phi, inner_mass=geometry["mass"],
            killing_frequency=float(packet["radial_moments"]["requested_omega"]),
            occupations=np.asarray(packet["occupations"], dtype=float),
        )
        initial = pg.stage_rates(coupled, grid)
        metric = freeze_sourced_metric(initial)
        fixed = _clone_state(coupled)
        initial_norm = float(coupled["probability_stocks"]["remaining"])
        step = matched_step(initial["cfl_dt"], geometry["r_m"])
        if probe_step_seconds is None:
            import time
            started = time.perf_counter()
            pg.rk4_step(_clone_state(coupled), grid, step["dt"])
            probe_step_seconds = time.perf_counter() - started
        steps_needed = int(math.ceil(limit * geometry["r_m"] / step["dt"]))
        forecast = admit(probe_step_seconds, steps_needed, remaining=budget)
        if not forecast["admitted"]:
            raise CampaignBlocker(
                "timing admission refused the station march", mode="run",
                admission=forecast, evolution_called=False,
            )
        stations = [item for item in STATION_TIMES_OVER_RM if item <= limit]
        reached = []
        bytes_written = 0
        spent = 0.0
        for station in stations:
            marched = march_pair(
                coupled, fixed, grid, metric,
                dt=step["dt"], station_time=station * geometry["r_m"],
                cpu_seconds=budget - spent,
            )
            coupled = marched["coupled"]
            fixed = marched["fixed"]
            spent += marched["cpu_seconds"]
            if marched["stopped"] == "cpu_budget":
                break
            rates = pg.stage_rates(coupled, grid)
            observation = station_observation(
                coupled, grid, rates, packet, initial_norm=initial_norm,
            )
            record = {
                "station_over_rm": station,
                "fixed_control": FIXED_CONTROL,
                "mass_over_rm": ratio,
                "energy_fraction": fraction,
                "observation": _plain(observation),
                "step": step,
                "producing_commit": producer_commit,
            }
            written = write_stage(
                directory, f"station-{station:g}", record,
                state_arrays(coupled, fixed, metric), blob_match=True,
            )
            bytes_written += written["bytes_written"]
            reached.append(station)
    except CampaignBlocker:
        raise
    except (TypeError, KeyError, IndexError, ValueError) as error:
        raise CampaignBlocker(
            "adapter gap: neutral PG contract: " + str(error),
            mode="run", adapter_gap=["nsc_rn_pg/nsc_rn_source contract"],
            evolution_called=False, mock_physical_output=False,
        ) from error
    return {
        "schema": SCHEMA,
        "mode": "run",
        "ok": True,
        "evolution_called": True,
        "mass_over_rm": ratio,
        "energy_fraction": fraction,
        "stations_reached": reached,
        "fixed_control": FIXED_CONTROL,
        "step": step,
        "production_trajectory": False,
        "bytes_written": bytes_written,
        "mock_physical_output": False,
        "root_review_overrides": False,
        "admission": decision,
        "timing": forecast,
        "cpu_seconds": spent,
    }


def write_dependency_manifest(directory, *, producer_commit, inputs):
    pins = producers()
    if not producer_commit or not frozen_blob_match(producer_commit, pins):
        raise CampaignBlocker("creation requires a frozen root commit", mode="dependency")
    if not inputs:
        raise CampaignBlocker("final writer requires the frozen input dependency chain", mode="dependency")
    record = {
        "producing_commit": producer_commit,
        "producers": pins,
        "inputs": bind_inputs(inputs),
        "dependency_order": [
            "nsc_rn_reference", "nsc_rn_pg", "nsc_rn_source",
            "nsc_rn_observables", "derive_nsc_rn_parent",
        ],
        "writer": "creation_only_after_frozen_root_commit",
        "evolution_called": False,
    }
    return write_stage(directory, "dependency", record, {}, blob_match=True)


def bind_inputs(paths):
    bound = {}
    for raw in paths:
        path = Path(raw).expanduser().resolve()
        if not path.is_file():
            raise CampaignBlocker(f"frozen input missing: {raw}", mode="dependency")
        bound[str(path)] = file_sha256(path)
    return bound


def _checkpoint_names(directory):
    names = []
    for path in sorted(Path(directory).glob("*.json")):
        stage = path.stem
        if stage in ("prepare", "dependency") or stage.startswith("station-"):
            names.append(stage)
    return names


def resume(directory):
    """Read an existing stage. Do not evolve and do not create one."""
    directory = Path(directory).expanduser().resolve()
    found = []
    arrays = {}
    for stage in _checkpoint_names(directory):
        json_path = directory / f"{stage}.json"
        npz_path = directory / f"{stage}.npz"
        if not json_path.is_file() or not npz_path.is_file():
            raise CampaignBlocker("checkpoint stage is incomplete", mode="resume")
        record = json.loads(json_path.read_text())
        if file_sha256(npz_path) != record.get("payload_sha256"):
            raise CampaignBlocker("checkpoint payload hash does not match", mode="resume")
        with np.load(npz_path, allow_pickle=False) as payload:
            loaded = {key: np.array(payload[key]) for key in payload.files}
        found.append({
            "stage": stage,
            "evolved": record.get("evolved"),
            "creation_only": record.get("creation_only"),
            "station_over_rm": record.get("station_over_rm"),
        })
        arrays[stage] = loaded
    if not found:
        raise CampaignBlocker("checkpoint not started", mode="resume")
    return {
        "schema": SCHEMA,
        "mode": "resume",
        "ok": True,
        "stages": found,
        "arrays": arrays,
        "evolved": False,
        "evolution_called": False,
        "bytes_written": 0,
        "production_trajectory": False,
        "mock_physical_output": False,
    }


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    return value


def forecast_first_pair():
    """CPU and memory forecast for M/r_m=1.04, ε=1e-3, up to T/r_m=5.

    Times one coupled RK4 step on a rank-one column. Does not march to the
    station and does not write.
    """
    import time
    geometry = case_geometry(1.04)
    pg = _load_module("nsc_rn_pg")
    route = geometry["route"]
    grid, state, rates = _load_module("nsc_rn_pg").sourcefree_rates(
        geometry["mass"], r_m=geometry["r_m"], points=int(route["points"]),
        r_out=route["outer"], killing_frequency=1.0,
    )[:3]
    phi = np.zeros((2, grid["points"], 1), dtype=complex)
    phi[0, :, 0] = 1.0e-3
    coupled = pg.make_state(
        grid, phi, occupations=np.array([1.0e-4]), inner_mass=float(geometry["mass"]),
        killing_frequency=1.0,
    )
    sourced = pg.stage_rates(coupled, grid)
    metric = freeze_sourced_metric(sourced)
    step = matched_step(sourced["cfl_dt"], geometry["r_m"])
    fixed = _clone_state(coupled)
    started = time.perf_counter()
    pg.rk4_step(_clone_state(coupled), grid, step["dt"])
    coupled_seconds = time.perf_counter() - started
    started = time.perf_counter()
    fixed_rk4(fixed, grid, metric, step["dt"])
    fixed_seconds = time.perf_counter() - started
    pair_seconds = coupled_seconds + fixed_seconds
    station = 5.0 * float(geometry["r_m"])
    steps = int(math.ceil(station / step["dt"]))
    forecast = ADMISSION_FACTOR * pair_seconds * steps
    derivative = np.asarray(grid["derivative"])
    memory_bytes = int(phi.nbytes + grid["radius"].nbytes + grid["weights"].nbytes + derivative.nbytes)
    return {
        "schema": SCHEMA,
        "mode": "forecast",
        "mass_over_rm": 1.04,
        "energy_fraction": 1.0e-3,
        "r_m": geometry["r_m"],
        "points": int(grid["points"]),
        "spacing": route["spacing"],
        "lambda_min": route["lambda_min"],
        "sigma": route["sigma"],
        "horizon_gap": geometry["horizon_gap"],
        "dt": step["dt"],
        "step": step,
        "cfl_half_step_not_used": step["limited_by"] == "step_cap",
        "station_time": station,
        "station_over_rm": 5.0,
        "steps_to_first_station": steps,
        "measured_coupled_rk4_seconds": coupled_seconds,
        "measured_fixed_rk4_seconds": fixed_seconds,
        "measured_pair_step_seconds": pair_seconds,
        "admission_factor": ADMISSION_FACTOR,
        "forecast_cpu_seconds": forecast,
        "forecast_cpu_seconds_by_station": {
            str(item): ADMISSION_FACTOR * pair_seconds * int(math.ceil(item * geometry["r_m"] / step["dt"]))
            for item in STATION_TIMES_OVER_RM
        },
        "cpu_budget_seconds": CPU_BUDGET_SECONDS,
        "within_cpu_budget": forecast <= CPU_BUDGET_SECONDS,
        "memory_bytes": memory_bytes,
        "memory_budget_bytes": 8 * 1024 ** 3,
        "within_memory_budget": memory_bytes <= 8 * 1024 ** 3,
        "chunk_bytes": CHUNK_BYTES,
        "later_stations_over_rm": STATION_TIMES_OVER_RM[1:],
        "fixed_control": FIXED_CONTROL,
        "fixed_background": "initial sourced metric held frozen; Dirac rates only",
        "coupled_state": "same initial phi; rk4_step rebuilds the source every stage",
        "admitted": False,
        "root_freeze_required": True,
        "probe_rk4_called": True,
        "evolution_called": False,
        "production_trajectory": False,
        "bytes_written": 0,
        "stations_reached": (),
    }


def _assessment_context_hashes(directory):
    """Hash retrospective context. Absence is recorded and is not a commit."""
    found = {}
    for name in (
        "summary.json", "initial-packet.json", "producer-hashes-before.json", "initial-step.json",
    ):
        path = directory / name
        found[name] = file_sha256(path) if path.is_file() else None
    receipt = None
    for candidate in (directory.parent / "receipt.json", directory.parent / "receipt", directory / "receipt.json"):
        if candidate.is_file():
            receipt = file_sha256(candidate)
            found[str(candidate)] = receipt
            break
    found["parent_receipt"] = receipt
    wrapper = None
    for candidate in (
        directory / "wrapper.py",
        directory.parent / "nsc_rn_pilot_march.py",
        directory.parent / "wrapper.py",
    ):
        if candidate.is_file():
            wrapper = file_sha256(candidate)
            found["original_wrapper"] = wrapper
            found[str(candidate)] = wrapper
            break
    if wrapper is None:
        found["original_wrapper"] = None
    return found


def _full_interval(rows, stride, value_key, rate_of):
    """Net trapezoid over a partition of the whole record.

    The largest single-interval gap is recorded and is not the full-run closure.
    """
    predicted = 0.0
    gaps = []
    for index in range(0, len(rows) - stride, stride):
        left = rows[index]
        right = rows[index + stride]
        step = right["t"] - left["t"]
        piece = 0.5 * step * (rate_of(left) + rate_of(right))
        predicted += piece
        gaps.append(right[value_key] - left[value_key] - piece)
    delta = rows[-1][value_key] - rows[0][value_key]
    return {
        "full_delta": delta,
        "physical_integral": predicted,
        "full_run_gap": delta - predicted,
        "interval_max_abs_gap": float(max(abs(item) for item in gaps)),
        "interval_max_is_not_full_run_closure": True,
    }


def assess_archive(directory):
    """Read saved states, reconstruct each once, and return scalar balances.

    No time step is taken. The fixed arm uses the archived initial sourced
    metric. The coupled arm uses one constraint reconstruction of that snapshot.
    """
    import platform
    import sys
    import time
    directory = Path(directory).expanduser().resolve()
    if not directory.is_dir():
        raise CampaignBlocker("archive path is not a directory", mode="assess")
    snapshots = sorted(directory.glob("t*.npz"), key=lambda item: float(item.stem[1:]))
    if len(snapshots) < 2:
        raise CampaignBlocker("archive does not contain a time sequence", mode="assess")
    cpu_started = time.process_time()
    wall_started = time.perf_counter()
    pg = _load_module("nsc_rn_pg")
    source = _load_module("nsc_rn_source")
    observables_mod = _load_module("nsc_rn_observables")
    geometry = case_geometry(FIRST_PAIR["mass_over_rm"])
    grid = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"], points=geometry["route"]["points"],
        r_out=geometry["route"]["outer"],
    )
    archive_hashes = {}
    context_hashes = _assessment_context_hashes(directory)
    sidecar_authenticated = []
    stations = []
    for path in snapshots:
        digest = file_sha256(path)
        archive_hashes[path.name] = digest
        sidecar = path.with_suffix(".json")
        if sidecar.is_file():
            archive_hashes[sidecar.name] = file_sha256(sidecar)
            claimed = json.loads(sidecar.read_text()).get("payload_sha256")
            if claimed is not None and claimed != digest:
                raise CampaignBlocker(
                    "archived npz hash does not match sidecar payload_sha256: " + path.name,
                    mode="assess",
                )
            if claimed is not None:
                sidecar_authenticated.append(path.name)
        with np.load(path, allow_pickle=False) as payload:
            arrays = {key: np.array(payload[key]) for key in payload.files}
        if not np.array_equal(arrays["radius"], grid["radius"]):
            raise CampaignBlocker("saved radius is not the first-pair route grid", mode="assess")
        template = pg.make_state(
            grid, arrays["coupled_phi"], occupations=arrays["occupations"],
            inner_mass=float(arrays["coupled_inner_mass"][0]), killing_frequency=1.0,
        )
        coupled = apply_checkpoint(template, arrays, "coupled")
        fixed = apply_checkpoint(template, arrays, "fixed")
        rates = pg.stage_rates(coupled, grid)
        frozen_lapse = arrays["frozen_lapse"]
        frozen_shift = arrays["frozen_shift"]
        fixed_rate, _probability_sat = pg.dirac_rates(grid, fixed["phi"], frozen_lapse, frozen_shift)
        fixed_jets = pg.metric_jets(grid["radius"], frozen_lapse, frozen_shift)
        coupled_energy = source.energy_account(
            coupled["phi"], coupled["occupations"], rates["lapse"], rates["shift"],
            grid["radius"], grid["weights"], grid["derivative"],
            phi_rate=rates["phi_rate"], lapse_t=rates["lapse_rate"], shift_t=rates["shift_rate"],
            lapse_r=rates["jets"]["lapse_r"], shift_r=rates["jets"]["shift_r"],
        )
        fixed_energy = source.energy_account(
            fixed["phi"], fixed["occupations"], frozen_lapse, frozen_shift,
            grid["radius"], grid["weights"], grid["derivative"],
            phi_rate=fixed_rate, lapse_t=np.zeros(grid["points"]), shift_t=np.zeros(grid["points"]),
            lapse_r=fixed_jets["lapse_r"], shift_r=fixed_jets["shift_r"],
        )
        normal = source.normal_energy_ledger(
            coupled["phi"], coupled["occupations"], rates["lapse"], rates["shift"],
            grid["radius"], grid["weights"], grid["derivative"],
            phi_rate=rates["phi_rate"], lapse_t=rates["lapse_rate"], shift_t=rates["shift_rate"],
            lapse_r=rates["jets"]["lapse_r"], shift_r=rates["jets"]["shift_r"],
        )
        curvature = observables_mod.curvature_from_pg_jets(
            grid["radius"], rates["lapse"], rates["shift"],
            lapse_r=rates["jets"]["lapse_r"], shift_r=rates["jets"]["shift_r"],
            lapse_rr=rates["jets"]["lapse_rr"], shift_rr=rates["jets"]["shift_rr"],
            lapse_t=rates["lapse_rate"], shift_t=rates["shift_rate"],
            shift_tr=rates["time_jets"]["shift_tr"],
        )
        vacuum = np.asarray(grid["model"].invariants(grid["radius"])["Ricci2"], dtype=float)
        ricci = np.asarray(curvature["Ricci2"], dtype=float)
        scalar = np.asarray(curvature["R4"], dtype=float)
        curvature_mask = np.ones(scalar.size, dtype=bool)
        curvature_mask[:10] = False
        curvature_mask[-10:] = False
        mass_from_energy = float(grid["G_N"]) * coupled_energy["e_H_density"] / rates["lapse"]
        stations.append({
            "tag": path.stem,
            "t": float(coupled["clocks"]["t"]),
            "t_over_rm": float(coupled["clocks"]["t"] / grid["r_m"]),
            "E_coordinate": coupled_energy["E_coordinate"],
            "E_normal": normal["E_normal"],
            "fixed_E_coordinate": fixed_energy["E_coordinate"],
            "fixed_E_normal": fixed_energy["E_normal"],
            "metric_power": coupled_energy["metric_power"],
            "boundary_energy_flux": coupled_energy["boundary_energy_flux"],
            "sat_energy_rate": coupled_energy["sat_energy_rate"],
            "phi_motion_rate": coupled_energy["phi_motion_rate"],
            "predicted_coordinate_rate": coupled_energy["predicted_coordinate_rate"],
            "J_N_left": normal["J_N_left"],
            "J_N_right": normal["J_N_right"],
            "pressure_integral": normal["pressure_integral"],
            "lapse_integral": normal["lapse_integral"],
            "sat_normal_energy": normal["sat_normal_energy"],
            "normal_global_residual": normal["global_residual"],
            "normal_point_max": normal["full_residual_max"],
            "normal_point_index": normal["full_residual_index"],
            "normal_point_radius": normal["full_residual_radius"],
            "normal_sbp_interior_max": normal["bulk_residual_max_outside_stencil"],
            "inner_bulk_residual": normal["bulk_residual_max"],
            "product_rule_beta": normal["product_rule_beta_max"],
            "R4_max_abs": float(np.max(np.abs(scalar))),
            "R4_bulk_10_max_abs": float(np.max(np.abs(scalar[curvature_mask]))),
            "Ricci2_max": float(np.max(ricci)),
            "Ricci2_vacuum_departure_max": float(np.max(np.abs(ricci - vacuum))),
            "Ricci2_vacuum_bulk_10_departure_max": float(
                np.max(np.abs(ricci[curvature_mask] - vacuum[curvature_mask]))
            ),
            "mass_slope_residual": float(np.max(np.abs(rates["mass_slope"] - mass_from_energy))),
            "outer_lapse_rate": float(rates["lapse_rate"][-1]),
            "occupation": float(coupled["occupations"][0]),
        })
    def normal_rate(row):
        return row["J_N_left"] - row["J_N_right"] + row["pressure_integral"] + row["lapse_integral"] + row["sat_normal_energy"]
    normal_quarter = _full_interval(stations, 1, "E_normal", normal_rate)
    normal_half = _full_interval(stations, 2, "E_normal", normal_rate)
    coordinate_quarter = _full_interval(stations, 1, "E_coordinate", lambda row: row["predicted_coordinate_rate"])
    coordinate_half = _full_interval(stations, 2, "E_coordinate", lambda row: row["predicted_coordinate_rate"])
    point = max(stations, key=lambda row: row["normal_point_max"])
    inner = max(stations, key=lambda row: row["inner_bulk_residual"])
    commutator = max(stations, key=lambda row: row["product_rule_beta"])
    return {
        "schema": SCHEMA,
        "mode": "assess",
        "ok": True,
        "evolution_called": False,
        "dynamics_loop": False,
        "bytes_written": 0,
        "production_trajectory": False,
        "mock_physical_output": False,
        "archive": str(directory),
        "frames": len(stations),
        "cpu_seconds": time.process_time() - cpu_started,
        "wall_seconds": time.perf_counter() - wall_started,
        "masks": {
            "normal_continuity": {
                "nodes_excluded_each_end": 1,
                "nodes_excluded_total": 2,
                "boundary_width": 1,
                "reason": "SBP endpoints",
            },
            "curvature_bulk": {
                "nodes_excluded_each_end": 10,
                "nodes_excluded_total": 20,
                "boundary_width": 10,
                "reason": "ten nodes at each endpoint",
            },
            "full_unmasked_max_retained": True,
            "cells_forced_to_zero": False,
        },
        "interval_max_is_not_full_run_closure": True,
        "lineage": {
            "producing_commit": None,
            "old89_produced_this_run": False,
            "producer_hashes_before_label": "retrospective_context_not_a_producing_commit",
            "sidecar_payloads_authenticated": sidecar_authenticated,
        },
        "context_hashes": context_hashes,
        "highlights": {
            "normal_point_max": point,
            "inner_bulk": inner,
            "product_rule": commutator,
            "R4_max_abs": float(max(row["R4_max_abs"] for row in stations)),
            "Ricci2_vacuum_departure_max": float(max(row["Ricci2_vacuum_departure_max"] for row in stations)),
            "mass_slope_residual_max": float(max(row["mass_slope_residual"] for row in stations)),
        },
        "balances": {
            "normal_dt_0.25": normal_quarter,
            "normal_dt_0.5": normal_half,
            "coordinate_dt_0.25": coordinate_quarter,
            "coordinate_dt_0.5": coordinate_half,
        },
        "stations": stations,
        "archive_hashes": archive_hashes,
        "producers": producers(),
        "environment": {
            "python": sys.version,
            "numpy": np.__version__,
            "platform": platform.platform(),
            "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
            "openblas_num_threads": os.environ.get("OPENBLAS_NUM_THREADS"),
        },
    }


def write_assessment(destination, report, *, producer_commit):
    """Creation-only JSON. Requires a frozen 40-character commit matching these bytes."""
    if not producer_commit or not frozen_blob_match(producer_commit, producers()):
        raise CampaignBlocker("creation requires a frozen root commit", mode="assess")
    destination = Path(destination).expanduser().resolve()
    if _sealed(destination):
        raise CampaignBlocker("existing outputs are sealed", mode="assess")
    if destination.exists():
        raise CampaignBlocker("assessment output already exists", mode="assess")
    body = dict(report)
    body["producing_commit"] = producer_commit
    body["creation_only"] = True
    encoded = (json.dumps(_plain(body), indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    if len(encoded) > CHUNK_BYTES:
        raise CampaignBlocker("assessment exceeds 64 MiB", mode="assess")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("xb") as handle:
        handle.write(encoded)
    return {
        "bytes_written": len(encoded),
        "path": str(destination),
        "sha256": hashlib.sha256(encoded).hexdigest(),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--calibrate", action="store_true")
    modes.add_argument("--prepare", action="store_true")
    modes.add_argument("--run", action="store_true")
    modes.add_argument("--resume", action="store_true")
    modes.add_argument("--assess-archive", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--producer-commit")
    parser.add_argument("--calibration-record", type=Path)
    parser.add_argument("--input", action="append", default=[], dest="inputs")
    parser.add_argument("--probe-step-seconds", type=float)
    parser.add_argument("--mass-over-rm", type=float)
    parser.add_argument("--energy-fraction", type=float)
    parser.add_argument("--station-limit", type=float)
    parser.add_argument("--cpu-seconds", type=float)
    parser.add_argument("--assessment-output", type=Path)
    args = parser.parse_args(argv)
    try:
        if args.assess_archive is not None:
            report = assess_archive(args.assess_archive)
            if args.assessment_output is not None:
                written = write_assessment(
                    args.assessment_output, report, producer_commit=args.producer_commit,
                )
                report = dict(report)
                report["bytes_written"] = written["bytes_written"]
                report["assessment_sha256"] = written["sha256"]
                report["assessment_path"] = written["path"]
        elif args.calibrate:
            report = calibrate()
        elif args.prepare or args.run or args.resume:
            if args.output is None:
                raise CampaignBlocker("checkpoint path is required")
            if args.resume:
                report = resume(args.output)
            else:
                if args.calibration_record is None:
                    record = None
                else:
                    record = json.loads(Path(args.calibration_record).read_text())
                if args.prepare:
                    report = prepare(
                        args.output, producer_commit=args.producer_commit,
                        calibration_record=record, inputs=tuple(args.inputs),
                    )
                else:
                    report = run(
                        args.output, producer_commit=args.producer_commit,
                        calibration_record=record, probe_step_seconds=args.probe_step_seconds,
                        mass_over_rm=args.mass_over_rm, energy_fraction=args.energy_fraction,
                        station_limit=args.station_limit, cpu_seconds=args.cpu_seconds,
                    )
        else:
            report = preview()
        code = 0
    except (CampaignBlocker, observables.ObservableError) as blocked:
        report = blocked.payload if isinstance(blocked, CampaignBlocker) else {
            "schema": SCHEMA,
            "ok": False,
            "blocker": str(blocked),
            "evolved": False,
            "evolution_called": False,
            "bytes_written": 0,
            "production_trajectory": False,
            "mock_physical_output": False,
            "root_review_overrides": False,
        }
        code = 0
    print(json.dumps(_plain(report), indent=2, sort_keys=True, allow_nan=False))
    return code


if __name__ == "__main__":
    raise SystemExit(main())
