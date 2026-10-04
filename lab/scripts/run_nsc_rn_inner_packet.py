#!/usr/bin/env python3
"""March a single RN radial case with the source packet released inside the collar.

Why this exists. The recorded neutral family places the source packet at
``PACKET_CENTER_OVER_RM = 12.0``, which is a construction choice. Measured
by the mass-shell reduction, that packet falls as a coherent shell at about
``-1.44`` per unit of ``t/r_m`` and crosses the grid's inner edge at
``t = 7.0`` in all six cases. By the end of every recorded run, between 61%
and 99% of the accreted mass lies *below* ``r = 1.04``. The accumulation,
which is the part of the picture that matters, happens outside the domain.

Putting the same packet lower changes what the computation can see. The
packet is still the same construction, same width, same spectrum; only its
radius is changed. The shell then forms inside the collar and falls a
shorter distance, so its arrival at the inner edge, and the mass profile it
leaves behind, are resolved rather than truncated.

This is a different construction from the recorded family. It is not a
replay of it and its numbers are not comparable with the recorded family
except through quantities that depend only on the construction, not on where
the packet was released.

Not a sealed record. Writes to a creation-only directory and refuses to
overwrite. No time step is taken by the reduction driver that reads this.

Run from the repository root:

    python scripts/lab.py scripts/run_nsc_rn_inner_packet.py --center 3.0
"""
from __future__ import annotations

import os

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import argparse
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

LAB = Path(__file__).resolve().parents[1]
for _path in (str(LAB / "src"), str(LAB / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

import derive_nsc_rn_parent as driver  # noqa: E402
from recursive_horizons import nsc_rn_pg as pg  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402

OUT = LAB / "results/development/nsc-rn-inner-packet-v1"

#: Station spacing in units of t/r_m. The recorded family used 0.5; the
#: inner packet crosses the collar in about 1.4, so finer stations are
#: worth the files.
STATION_STEP = 0.1
#: Hard cap on the march, in units of t/r_m.
TARGET_OVER_RM = 3.0


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _station_name(t_over_rm: float) -> str:
    """Recorded convention: the name encodes t/r_m at 0.1 resolution.

    ``lab.py`` runs with the working directory at ``lab/``, so the same
    helper is used here to keep new stations interchangeable with the
    recorded family.
    """
    deci = int(round(float(t_over_rm) * 10.0))
    if deci % 10 == 0:
        return f"station-{deci // 10}"
    return f"station-{deci // 10}.{deci % 10}"


def _write_station(path: Path, coupled, fixed, grid, rates) -> Dict[str, Any]:
    """Write the arrays the mass-shell reduction and its controls consume.

    Curvature is the geometry itself, so it is recorded rather than
    reconstructed later. It comes from the evolved lapse and shift and
    their jets, not from the vacuum reference, which is why this is the
    same-state quantity the curvature provenance note asks for. It stays
    null when the required time-jet primitive is unavailable, and that is
    recorded rather than hidden.
    """
    phi = np.asarray(coupled["phi"], dtype=np.complex128)
    observation = driver.station_observation(coupled, grid, rates)
    arrays = {
        "coupled_phi": phi,
        "coupled_mass": np.asarray(rates["mass"], dtype=np.float64),
        "coupled_inner_mass": np.asarray([float(coupled["inner_mass"])], dtype=np.float64),
        "coupled_t": np.asarray([float(coupled["clocks"]["t"])], dtype=np.float64),
        "coupled_killing": np.asarray([float(coupled["clocks"]["killing"])], dtype=np.float64),
        "coupled_normal": np.asarray([float(coupled["clocks"]["normal"])], dtype=np.float64),
        "occupations": np.asarray(coupled["occupations"], dtype=np.float64),
        "radius": np.asarray(grid["radius"], dtype=np.float64),
        "frozen_lapse": np.asarray(rates["lapse"], dtype=np.float64),
        "frozen_shift": np.asarray(rates["shift"], dtype=np.float64),
    }
    if observation.get("R4") is not None and observation.get("Ricci2") is not None:
        arrays["R4"] = np.asarray(observation["R4"], dtype=np.float64)
        arrays["Ricci2"] = np.asarray(observation["Ricci2"], dtype=np.float64)
    np.savez_compressed(path, **arrays)
    return {
        "curvature_status": observation.get("curvature_status"),
        "R4_recorded": observation.get("R4") is not None,
        "trapping_horizon": observation.get("trapping_horizon"),
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mu", type=float, default=1.04)
    parser.add_argument("--eps", type=float, default=1.0e-3)
    parser.add_argument("--center", type=float, default=3.0,
                        help="packet centre in units of r_m")
    parser.add_argument("--target", type=float, default=TARGET_OVER_RM,
                        help="march target in units of t/r_m")
    parser.add_argument("--points", type=int, default=None)
    parser.add_argument("--excision", type=float, default=None,
                        help="inner edge of the grid, overriding the horizon midpoint. "
                             "Must lie strictly between r_minus and r_plus.")
    parser.add_argument("--out", default="results/development/nsc-rn-inner-packet-v1",
                        help="output directory, relative to the lab root")
    parser.add_argument("--cpu-cap", type=float, default=1800.0)
    args = parser.parse_args(argv)

    # ``lab.py`` runs with the working directory at ``lab/``, so a relative
    # --out is resolved against the lab root rather than the process CWD.
    out = Path(args.out)
    out = out if out.is_absolute() else (LAB / out)
    if out.exists() and any(out.iterdir()):
        print(f"refusing to overwrite non-empty record: {out}", file=sys.stderr)
        return 1
    out.mkdir(parents=True, exist_ok=True)

    geometry = driver.case_geometry(args.mu)
    points = int(args.points or geometry["route"]["points"])
    grid = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"], points=points,
        r_out=geometry["route"]["outer"], excision=args.excision,
    )

    r_inner = float(grid["radius"][0])
    r_outer = float(grid["radius"][-1])
    r_minus = float(geometry["r_minus"]) if "r_minus" in geometry else None
    if r_minus is not None and not r_minus < r_inner:
        raise SystemExit("the grid's inner edge is not above the inner horizon")

    # The whole point of this construction is that the packet forms inside
    # the collar. Refuse anything that would put it back outside, and refuse
    # a centre whose *left edge* already sits below the inner boundary. The
    # left edge is the binding condition, not the centre: the packet is two
    # r_m wide, so a centre of 2.0 puts half of it outside the collar and
    # the shell is born pinned to the boundary with no fall at all. That
    # case was run once and is kept in the record as a boundary, but it is
    # not a measurement of a shell falling.
    width = float(source.PACKET_WIDTH_OVER_RM)
    left_edge = args.center - 0.5 * width
    if not r_inner < args.center < r_outer:
        raise SystemExit(f"packet centre {args.center} is not strictly inside the collar "
                         f"[{r_inner}, {r_outer}]")
    if left_edge <= r_inner:
        raise SystemExit(f"packet left edge {left_edge:.4f} is at or below the inner edge "
                         f"{r_inner}; with a width of {width} r_m the shell would be born "
                         f"pinned to the boundary and no fall would be measured")
    if r_minus is not None and left_edge <= r_minus:
        raise SystemExit("packet left edge is at or inside the inner horizon")

    started = time.perf_counter()
    packet = source.prepare_radial_packet(
        mass_over_rm=args.mu, energy_target=args.eps, rank=driver.PILOT_RANK,
        model=grid["model"], radii=grid["radius"], weights=grid["weights"],
        derivative=grid["derivative"], center=args.center * geometry["r_m"],
    )
    prepare_seconds = time.perf_counter() - started
    admission = driver.admit_prepared_packet(packet)

    phi = np.asarray(packet["phi"])
    coupled = pg.make_state(
        grid, phi, inner_mass=geometry["mass"],
        killing_frequency=float(packet["radial_moments"]["requested_omega"]),
        occupations=np.asarray(packet["occupations"], dtype=float),
    )
    initial_rates = pg.stage_rates(coupled, grid)
    metric = driver.freeze_sourced_metric(initial_rates)
    fixed = driver._clone_state(coupled)
    if not np.array_equal(coupled["phi"], fixed["phi"]):
        raise SystemExit("the frozen control did not start from the same phi")

    step = driver.matched_step(initial_rates["cfl_dt"], geometry["r_m"])
    if step["dt"] > 0.001 * geometry["r_m"] + 1e-15 or step["dt"] > step["cfl_dt"]:
        raise SystemExit("step exceeds the cap or the full CFL")

    stations: List[Dict[str, Any]] = []
    first_report = _write_station(out / _station_name(0.0), coupled, fixed, grid, initial_rates)
    stations.append({"t_over_rm": 0.0,
                     "shell_inner_mass": float(initial_rates["mass"][0]),
                     "shell_outer_mass": float(initial_rates["mass"][-1]),
                     **first_report})

    cpu = 0.0
    wall0 = time.perf_counter()
    stop_reason = "not_started"
    target = STATION_STEP
    while target <= args.target + 1e-12:
        if cpu >= args.cpu_cap:
            stop_reason = "cpu_budget"
            break
        remaining = args.cpu_cap - cpu
        chunk_time = target * geometry["r_m"]
        marched = driver.march_pair(
            coupled, fixed, grid, metric, dt=step["dt"],
            station_time=chunk_time,
            cpu_seconds=remaining * driver.ADMISSION_FACTOR,
        )
        coupled = marched["coupled"]
        fixed = marched["fixed"]
        cpu += float(marched["cpu_seconds"])
        rates = pg.stage_rates(coupled, grid)
        mass = np.asarray(rates["mass"], dtype=float)
        t_over = float(coupled["clocks"]["t"]) / geometry["r_m"]
        report = _write_station(out / _station_name(t_over), coupled, fixed, grid, rates)
        stations.append({
            "t_over_rm": t_over,
            "shell_inner_mass": float(mass[0]),
            "shell_outer_mass": float(mass[-1]),
            "stopped": marched["stopped"],
            **report,
        })
        print(f"t/rm={t_over:.4f} mass[0]={mass[0]:.7f} cpu={cpu:.1f}s "
              f"stop={marched['stopped']}", flush=True)
        if marched["stopped"] == "cpu_budget":
            stop_reason = "cpu_budget"
            break
        target += STATION_STEP
    else:
        stop_reason = "target"

    summary = {
        "schema": "NSC-RN-INNER-PACKET-v1",
        "creation_only": True,
        "not_a_sealed_record": True,
        "producing_commit": _git_head(),
        "producer_sha256": _sha256_file(Path(__file__).resolve()),
        "mu": args.mu,
        "eps": args.eps,
        "packet_center_over_rm": args.center,
        "packet_width_over_rm": source.PACKET_WIDTH_OVER_RM,
        "recorded_family_packet_center_over_rm": source.PACKET_CENTER_OVER_RM,
        "points": int(grid["points"]),
        "r_inner": r_inner,
        "r_outer": r_outer,
        "r_minus": r_minus,
        "excision_default": bool(grid["excision_is_default"]),
        "excision_requested": args.excision,
        "horizon_midpoint": 0.5 * (r_minus + float(geometry["route"]["r_plus"]))
        if "r_plus" in geometry else None,
        "target_over_rm": args.target,
        "station_step_over_rm": STATION_STEP,
        "dt": float(step["dt"]),
        "prepare_seconds": prepare_seconds,
        "cpu_seconds": cpu,
        "wall_seconds": time.perf_counter() - wall0,
        "stop_reason": stop_reason,
        "stations": len(stations),
        "admission": admission,
        "recorded_stations": stations,
        "why_this_is_a_different_construction": (
            "The recorded family releases the packet at 12 r_m, from outside the "
            "collar it is about to cross. This one releases it at the centre above, "
            "inside. Same width, same spectrum, same chart. The mass profile left "
            "behind is a different problem and is not comparable with the recorded "
            "family except through construction-only quantities."
        ),
        "scope_note": (
            "A construction check on one finite marched collar. It does not measure "
            "the actual universe."
        ),
    }
    (out / "summary.json").write_text(
        json.dumps(summary, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {out}: {len(stations)} stations, stop={stop_reason}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())