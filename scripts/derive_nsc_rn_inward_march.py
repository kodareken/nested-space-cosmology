#!/usr/bin/env python3
"""Inward-march campaign: extend the RN radial grid to the inner horizon.

The recorded family stops at the horizon midpoint (the classical inner
mass), which sits OUTSIDE the RN inner horizon in every case. The
gradient-bubble note therefore says regularity of the interior support
cannot be decided on that domain. This runner extends the areal grid
inward from the midpoint down to just above r_minus, keeping the spacing,
the outer boundary, the source packet, and the coupled/frozen two-arm
march unchanged, and marches the same neutral shell so the force slopes
are measured on a domain that reaches the inner horizon.

The new inner edge stays strictly above r_minus: the chart requires pure
outflow there, and just above r_minus the outgoing speed is still
negative (trapped), so the extension is admissible. Below r_minus the
outgoing speed is positive and the excision would emit; that region is a
new construction and is not opened here.

Writer. Creation-only per case directory; a failure file is written with
a structured blocker on any exception, and no station is committed from a
partial step. Producing commit is recorded from git and the runner pins
its own bytes.

Run one case per process:

    python scripts/lab.py scripts/derive_nsc_rn_inward_march.py --mu 1.04 --eps 0.001
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Dict

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_name, "1")

import numpy as np

LAB = Path(__file__).resolve().parents[1]
for _path in (LAB / "src", LAB / "scripts"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from recursive_horizons import nsc_rn_observables as observables  # noqa: E402
from recursive_horizons import nsc_rn_pg as pg  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402
import derive_nsc_rn_parent as driver  # noqa: E402

SCHEMA = "NSC-RN-INWARD-MARCH-v1"
CAMPAIGN = LAB / "results/development/nsc-rn-inward-march-v1"
STATION_TIMES = (5.0, 10.0, 15.0, 20.0)
HORIZON_TOL_FRACTION = 0.05  # inner edge sits spacing * this above r_minus
RANK = driver.PILOT_RANK


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _producing_commit() -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=LAB, capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()
    except Exception:
        return None


def _exclusive_bytes(path: Path, data: bytes) -> int:
    descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o444)
    with os.fdopen(descriptor, "wb") as handle:
        handle.write(data)
    return len(data)


def to_jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        array = np.asarray(value)
        if array.size > 8:
            summary: Dict[str, Any] = {
                "shape": [int(item) for item in array.shape],
                "dtype": str(array.dtype),
                "sha256": driver.array_sha256(array),
            }
            if not np.iscomplexobj(array) and array.size:
                summary["max_abs"] = float(np.max(np.abs(array)))
            return summary
        return driver._plain(array)
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [to_jsonable(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    return value


def extended_grid(mass_over_rm: float) -> Dict[str, Any]:
    """Production grid whose inner edge sits just above r_minus."""
    geometry = driver.case_geometry(mass_over_rm)
    base = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"],
        points=geometry["route"]["points"], r_out=geometry["route"]["outer"],
    )
    r_minus = geometry["r_minus"]
    r0 = float(base["radius"][0])
    spacing = float(base["spacing"])
    added = int(math.ceil((r0 - r_minus) / spacing))
    inner = r0 - added * spacing
    if inner <= r_minus:
        inner = r_minus + HORIZON_TOL_FRACTION * spacing
        added = int(math.floor((r0 - inner) / spacing))
        inner = r0 - added * spacing
    if inner <= r_minus or added < 1:
        raise RuntimeError("inward extension did not open at least one cell above r_minus")
    lower = np.linspace(inner, r0 - spacing, added)
    radius = np.concatenate([lower, np.asarray(base["radius"], dtype=float)])
    weights, derivative, q_matrix = pg.sbp_operator(int(radius.size), spacing, layout="csr")
    grid = dict(base)
    grid.update({
        "radius": radius,
        "weights": weights,
        "derivative": derivative,
        "q_matrix": q_matrix,
        "spacing": spacing,
        "points": int(radius.size),
    })
    grid["inner_edge"] = float(inner)
    grid["r_minus"] = float(r_minus)
    grid["inner_edge_gap"] = float(inner - r_minus)
    grid["route_points"] = grid["points"]
    del grid["cuts"]  # cuts referred to the old excision; keep the march honest
    return grid, geometry


def stress_report(phi, occupations, lapse, shift, radius, weights, derivative):
    evaluated = source.source_evaluate(
        phi, occupations, lapse, shift, np.ones_like(radius), radius, weights, derivative,
    )
    rho = np.asarray(evaluated["physical"]["rho"], dtype=float)
    radial = np.asarray(evaluated["physical"]["radial_pressure"], dtype=float)
    angular = np.asarray(evaluated["physical"]["angular_pressure"], dtype=float)
    return rho, radial, angular


def station_arrays(coupled, fixed, metric, grid, rates, rho, radial, angular,
                   observation) -> Dict[str, np.ndarray]:
    arrays = driver.state_arrays(coupled, fixed, metric)
    arrays.update({
        "radius": np.asarray(grid["radius"], dtype=float),
        "coupled_mass": np.asarray(rates["mass"], dtype=float),
        "rho": rho,
        "radial_pressure": radial,
        "angular_pressure": angular,
    })
    if observation.get("R4") is not None and observation.get("Ricci2") is not None:
        arrays.update({
            "R4": np.asarray(observation["R4"], dtype=float),
            "Ricci2": np.asarray(observation["Ricci2"], dtype=float),
        })
    return arrays


def write_payload(directory: Path, name: str, record: Dict[str, Any],
                  arrays: Dict[str, np.ndarray]) -> Dict[str, Any]:
    buffer = io.BytesIO()
    np.savez(buffer, **arrays)
    payload_bytes = buffer.getvalue()
    payload_sha = _sha256_bytes(payload_bytes)
    record = dict(record)
    record.update({
        "schema": SCHEMA,
        "creation_only": True,
        "payload_sha256": payload_sha,
        "payload_bytes": len(payload_bytes),
    })
    json_path = directory / f"{name}.json"
    npz_path = directory / f"{name}.npz"
    if json_path.exists() or npz_path.exists():
        raise RuntimeError(f"refusing to overwrite sealed station {name}")
    _exclusive_bytes(npz_path, payload_bytes)
    try:
        _exclusive_bytes(json_path, json.dumps(record, indent=1, sort_keys=True).encode("utf-8") + b"\n")
    except FileExistsError:
        return {"written": npz_path.name, "json_refused": True, "bytes_written": len(payload_bytes)}
    return {"written": name, "bytes_written": len(payload_bytes)}


def run_case(mu: float, eps: float, campaign: Path | None = None,
             stations: tuple[float, ...] | None = None) -> int:
    campaign = CAMPAIGN if campaign is None else campaign
    stations = STATION_TIMES if stations is None else stations
    case_name = f"mu-{mu:.3f}-eps-{eps:.3f}"
    case_dir = campaign / case_name
    runner_sha = _sha256_file(Path(__file__))
    produce_commit = _producing_commit()
    geometry = None
    coupled = fixed = None
    clock_now = 0.0
    stage = "extension"
    try:
        campaign.mkdir(parents=True, exist_ok=True)
        grid, geometry = extended_grid(mu)
        if (case_dir / "prepare.json").exists():
            print(f"case already prepared: {case_dir}", file=sys.stderr)
            return 1
        case_dir.mkdir()
        packet = source.prepare_radial_packet(
            mass_over_rm=mu, energy_target=eps, rank=RANK, model=grid["model"],
            radii=grid["radius"], weights=grid["weights"], derivative=grid["derivative"],
        )
        phi = np.asarray(packet["phi"])
        coupled = pg.make_state(
            grid, phi, inner_mass=geometry["mass"],
            killing_frequency=float(packet["radial_moments"]["requested_omega"]),
            occupations=np.asarray(packet["occupations"], dtype=float),
        )
        initial = pg.stage_rates(coupled, grid)
        metric = driver.freeze_sourced_metric(initial)
        fixed = driver._clone_state(coupled)
        initial_norm = float(coupled["probability_stocks"]["remaining"])
        step = driver.matched_step(initial["cfl_dt"], geometry["r_m"])

        prepare = {
            "schema": SCHEMA,
            "creation_only": True,
            "mode": "inward_march",
            "mass_over_rm": mu,
            "energy_fraction": eps,
            "r_minus": geometry["r_minus"],
            "inner_edge": grid["inner_edge"],
            "inner_edge_gap": grid["inner_edge_gap"],
            "grid_points": grid["points"],
            "spacing": grid["spacing"],
            "outer": grid["r_out"],
            "packet_center": packet["radial_moments"]["center"],
            "packet_width": packet["radial_moments"]["width_rms"],
            "station_times": [float(s) for s in stations],
            "step_cap_over_rm": step["step_cap_over_rm"],
            "producing_commit": produce_commit,
            "runner_sha256": runner_sha,
            "evolution_called": False,
        }
        _exclusive_bytes(case_dir / "prepare.json",
                         json.dumps(prepare, indent=1, sort_keys=True).encode("utf-8") + b"\n")

        stage = "march"
        for station in stations:
            marched = driver.march_pair(
                coupled, fixed, grid, metric,
                dt=step["dt"], station_time=station * geometry["r_m"],
                cpu_seconds=None,
            )
            coupled = marched["coupled"]
            fixed = marched["fixed"]
            clock_now = float(coupled["clocks"]["t"])
            rates = pg.stage_rates(coupled, grid)
            observation = driver.station_observation(
                coupled, grid, rates, packet=packet, initial_norm=initial_norm,
            )
            rho, radial, angular = stress_report(
                coupled["phi"], coupled["occupations"], rates["lapse"], rates["shift"],
                grid["radius"], grid["weights"], grid["derivative"],
            )
            arrays = station_arrays(coupled, fixed, metric, grid, rates, rho, radial, angular, observation)
            record = {
                "station_over_rm": float(station),
                "clock_t_over_rm": clock_now / geometry["r_m"],
                "mass_over_rm": mu,
                "energy_fraction": eps,
                "r_minus": geometry["r_minus"],
                "inner_edge": grid["inner_edge"],
                "inner_edge_gap": grid["inner_edge_gap"],
                "points": grid["points"],
                "fixed_control": driver.FIXED_CONTROL,
                "closure_residual": observation["sat_norm_stock"]["probability_closure_residual"],
                "charged_mass_inner": observation["charged_mass_inner"],
                "trapping_horizon": observation["trapping_horizon"],
                "producing_commit": produce_commit,
                "runner_sha256": runner_sha,
            }
            result = write_payload(case_dir, f"station-{station:g}", record, arrays)
            print(f"  {case_name} station {station:g} clock {clock_now:.3f} -> {result['written']}", flush=True)
    except Exception as error:  # noqa: BLE001 -- recorded as a structured blocker
        if case_dir is not None and case_dir.exists():
            failure = {
                "schema": SCHEMA,
                "blocker": getattr(error, "blocker", None) or str(error),
                "kind": type(error).__name__,
                "stage": stage,
                "clock_t": clock_now,
                "mass_over_rm": mu,
                "energy_fraction": eps,
                "runner_sha256": runner_sha,
                "producing_commit": _producing_commit(),
                "traceback_tail": traceback.format_exc()[-4000:],
            }
            path = case_dir / "station-failure.json"
            if not path.exists():
                path.write_text(json.dumps(failure, indent=1, sort_keys=True), encoding="utf-8")
        print(f"{case_name} blocker at stage {stage}: {error}", file=sys.stderr)
        return 1
    print(f"{case_name} complete: {len(stations)} stations")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mu", type=float, required=True)
    parser.add_argument("--eps", type=float, required=True)
    parser.add_argument("--campaign", default=str(CAMPAIGN.relative_to(LAB)),
                        help="creation-only campaign directory")
    parser.add_argument("--stations", type=float, nargs="+", default=list(STATION_TIMES),
                        help="station times over r_m, in increasing order")
    args = parser.parse_args()
    if args.mu not in driver.MASS_OVER_RM or args.eps not in driver.ENERGY_FRACTIONS:
        print("mu and eps must belong to the magnetic family", file=sys.stderr)
        return 1
    stations = tuple(sorted(float(s) for s in args.stations))
    if not stations or any(s <= 0.0 for s in stations):
        print("stations must be positive and increasing", file=sys.stderr)
        return 1
    return run_case(args.mu, args.eps, campaign=(LAB / args.campaign).resolve(), stations=stations)


if __name__ == "__main__":
    raise SystemExit(main())