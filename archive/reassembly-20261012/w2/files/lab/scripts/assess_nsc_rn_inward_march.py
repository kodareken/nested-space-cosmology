"""Assess the inward-march stations: decide the Hayward criterion where it lives.

Reads the stations written by ``derive_nsc_rn_inward_march`` (one cell
above the RN inner horizon) and reports the compressive-force log-slope
on the actually evaluated inner region. Creation-only record.

Run from the repository root:

    python scripts/lab.py scripts/assess_nsc_rn_inward_march.py
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np

LAB = Path(__file__).resolve().parents[1]
for _p in (LAB / "src", LAB / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from recursive_horizons import nsc_rn_cone_map as cone  # noqa: E402
from recursive_horizons import nsc_rn_gradient_bubble as bubble  # noqa: E402
from recursive_horizons import nsc_rn_inward as inward  # noqa: E402

CAMPAIGN = LAB / "results/development/nsc-rn-inward-march-v1"
RECORD = LAB / "results/development/nsc-rn-inward-march-assessment-v1.json"
CASE_PATTERN = re.compile(r"^mu-(?P<mu>[0-9.]+)-eps-(?P<eps>[0-9.]+)$")
STATION_PATTERN = re.compile(r"^station-(?P<t>[0-9.]+)\.npz$")


def _producing_commit() -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=LAB, capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()
    except Exception:
        return None


def _station_time(path: Path) -> float:
    match = STATION_PATTERN.match(path.name)
    return float(match.group("t")) if match else float("inf")


def _station_row(station: Path, mu: float, r_minus: float) -> Dict[str, Any]:
    with np.load(station) as payload:
        station_t = float(np.asarray(payload["coupled_t"])[0])
        grid = cone.station_grid(dict(payload))
        decision = inward.inner_region_decision(payload["radius"], payload["radial_pressure"], r_minus)
        assessment = bubble.gradient_bubble_assessment(payload, {"R4": payload["R4"]} if "R4" in payload.files else None,
                                                       mass_over_rm=mu)
        trap = cone.cone_profile(grid, cone.coupled_metric(dict(payload), grid)["geometry"])["trap_radius"]
    return {
        "station": station.stem,
        "t_over_rm": station_t,
        "throat_slope_inner": decision["throat_slope_inner"],
        "inner_region_slope_min": decision["inner_region_slope_min"],
        "meets_hayward_on_evaluated_domain": decision["meets_hayward_on_evaluated_domain"],
        "gap_to_horizon_cells": decision["gap_in_cells"],
        "trap_radius": trap,
        "gradient_bubble_verdict": assessment["verdict"],
        "wall_count": assessment["profile"]["wall_count"],
    }


def main() -> int:
    if not CAMPAIGN.is_dir():
        print(f"missing campaign directory: {CAMPAIGN}", file=sys.stderr)
        return 1
    if RECORD.exists():
        print(f"refusing to overwrite sealed record: {RECORD}", file=sys.stderr)
        return 1

    cases = []
    for case_dir in sorted(CAMPAIGN.iterdir()):
        match = CASE_PATTERN.match(case_dir.name)
        if not case_dir.is_dir() or match is None:
            continue
        mu = float(match.group("mu"))
        prepare = json.loads((case_dir / "prepare.json").read_text())
        r_minus = float(prepare["r_minus"])
        stations = sorted(
            (p for p in case_dir.glob("station-*.npz") if STATION_PATTERN.match(p.name)),
            key=_station_time,
        )
        rows = [_station_row(s, mu, r_minus) for s in stations]
        last = rows[-1] if rows else None
        cases.append({
            "case": case_dir.name,
            "mass_over_rm": mu,
            "r_minus": r_minus,
            "inner_edge": prepare["inner_edge"],
            "station_count": len(rows),
            "meets_at_last_station": last["meets_hayward_on_evaluated_domain"] if last else None,
            "throat_slope_last": last["throat_slope_inner"] if last else None,
            "inner_region_slope_min_last": last["inner_region_slope_min"] if last else None,
            "stations": rows,
        })
    cases.sort(key=lambda c: c["mass_over_rm"])

    record = {
        "schema": "NSC-RN-INWARD-MARCH-ASSESSMENT-v1",
        "api_version": inward.API_VERSION,
        "creation_only": True,
        "producing_commit": _producing_commit(),
        "campaign": CAMPAIGN.name,
        "criterion": (
            "Hayward gr-qc/0506126 (2005): a static region supported by finite "
            "density and pressure behaves as r^g with g >= -1 inward. Decided on "
            "the evaluated cells one cell above r_minus, no further extrapolation."
        ),
        "coverage": {
            "cases": len(cases),
            "total_stations": sum(c["station_count"] for c in cases),
        },
        "summary": {
            "meets_hayward_at_last_station": {
                c["case"]: c["meets_at_last_station"] for c in cases
            },
            "throat_slope_last_station": {
                c["case"]: c["throat_slope_last"] for c in cases
            },
            "inner_region_slope_min_last": {
                c["case"]: c["inner_region_slope_min_last"] for c in cases
            },
        },
        "cases": cases,
    }

    RECORD.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {RECORD.relative_to(LAB)}")
    for case in cases:
        print(
            f"  {case['case']:22s} stations={case['station_count']:3d} "
            f"slope_last={case['throat_slope_last']} meets={case['meets_at_last_station']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())