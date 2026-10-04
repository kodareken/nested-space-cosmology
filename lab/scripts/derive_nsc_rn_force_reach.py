"""Measure the effective reach of the two RN radial forces across the family.

Answers one question put by Douglas: does the larger force dominate a
larger region and reach further, with nothing pushing back, and does the
opposing force simply cancel away once other zones appear.

Reads every recorded station already written by
``nsc-rn-neutral-family-v1``. No time step is taken, nothing is re-run,
and no existing file is rewritten.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_force_reach.py
"""
from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_rn_force_reach as reach  # noqa: E402
from recursive_horizons import nsc_rn_gradient_bubble as bubble  # noqa: E402

FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"
RECORD = LAB / "results/development/nsc-rn-force-reach-v1.json"
CASE_PATTERN = re.compile(r"^mu-(?P<mu>[0-9.]+)-eps-(?P<eps>[0-9.]+)$")
STATION_PATTERN = re.compile(r"^station-(?P<t>[0-9.]+)\.npz$")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _station_time(path: Path) -> float:
    match = STATION_PATTERN.match(path.name)
    return float(match.group("t")) if match else float("inf")


def _producing_commit() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _reduce_case(case_dir: Path) -> Dict[str, Any]:
    match = CASE_PATTERN.match(case_dir.name)
    stations = sorted((p for p in case_dir.glob("station-*.npz")
                       if STATION_PATTERN.match(p.name)), key=_station_time)
    rows: List[Dict[str, Any]] = []
    points = None
    for station in stations:
        with np.load(station) as payload:
            if points is None:
                points = int(payload["radius"].shape[0])
            profile = reach.force_reach(payload)
            slope = bubble._log_slope(payload["radius"],
                                      np.abs(payload["radial_pressure"]))
            throat_slope = float(np.nanmin(slope[1:31]))
            clock = (float(payload["coupled_t"][0])
                     if "coupled_t" in payload.files else _station_time(station))
        rows.append({
            "station": station.stem,
            "throat_slope_inner": throat_slope,
            "t_over_rm": clock,
            "reach_inner": profile["reach_inner"],
            "reach_outer": profile["reach_outer"],
            "reach_longer": profile["reach_longer"],
            "larger_force": profile["larger_force"],
            "dominant_fraction_inner": profile["dominant_fraction_inner"],
            "crossover_count": profile["crossover_count"],
            "verdict": reach.force_reach_verdict(profile),
        })
    verdicts = sorted({row["verdict"] for row in rows})
    last = rows[-1] if rows else None
    return {
        "case": case_dir.name,
        "mass_over_rm": float(match.group("mu")) if match else None,
        "epsilon": float(match.group("eps")) if match else None,
        "points": points,
        "station_count": len(rows),
        "reached_t_over_rm": last["t_over_rm"] if last else None,
        "verdicts": verdicts,
        "larger_force_stations": len([r for r in rows
                                      if r["larger_force"] == "radial"]),
        "final": last,
        "stations": rows,
    }


def main() -> int:
    if not FAMILY.is_dir():
        print(f"missing family directory: {FAMILY}", file=sys.stderr)
        return 1
    if RECORD.exists():
        print(f"refusing to overwrite sealed record: {RECORD}", file=sys.stderr)
        return 1

    cases = [_reduce_case(p) for p in sorted(FAMILY.iterdir())
             if p.is_dir() and CASE_PATTERN.match(p.name)]

    finals = [c["final"] for c in cases if c["final"]]
    longer = [row["reach_longer"] for row in finals]
    verdicts = sorted({v for c in cases for v in c["verdicts"]})
    supporting = len([c for c in cases if c["verdicts"] == ["larger_reaches_further"]])
    pairs = [(c["final"]["throat_slope_inner"], c["final"]["reach_inner"])
             for c in cases if c["final"]]
    correlation = (float(np.corrcoef(np.array(pairs)[:, 0], np.array(pairs)[:, 1])[0, 1])
                   if len(pairs) > 2 else None)

    record = {
        "schema": "NSC-RN-FORCE-REACH-v1",
        "api_version": reach.API_VERSION,
        "creation_only": True,
        "producing_commit": _producing_commit(),
        "family_source": FAMILY.name,
        "family_source_sha256": _sha256_file(FAMILY / "prepare.json"),
        "question": (
            "Does the larger of the two forces dominate a larger region and reach "
            "further, with nothing pushing back, and does the opposing force cancel "
            "away once other zones appear?"
        ),
        "domain": {
            "chart": "areal radius",
            "r_outer": 32.0,
            "observer": "areal chart of the recorded coupled arm",
            "state": "recorded station payload, no time step taken",
            "compared_quantity": (
                "radius at which each force stays above a declared fraction of its "
                "own peak, taken at the last crossing"
            ),
        },
        "scope_note": (
            "A construction check on recorded finite states. It does not measure the "
            "actual universe, and a reach here is a property of one computed collar on "
            "one grid, not a force law."
        ),
        "definition": {
            "reach_floor": reach.REACH_FLOOR,
            "last_crossing": (
                "reach is the last node at which the magnitude is still above the floor, "
                "so near-throat ringing cannot make a force look like it died early"
            ),
            "supported_prefix": (
                "the node count inside which the force never returns above the floor, so a "
                "later dip does not shorten the reported reach"
            ),
            "dominance": "fraction of grid nodes where the radial magnitude is larger",
        },
        "coverage": {
            "cases": len(cases),
            "total_stations": sum(c["station_count"] for c in cases),
            "cases_supporting_larger_reaches_further": supporting,
        },
        "verdicts": verdicts,
        "reach_longer_counts": {value: longer.count(value) for value in sorted(set(longer))},
        "summary": {
            "reach_inner_min": min((row["reach_inner"] for row in finals), default=None),
            "reach_inner_max": max((row["reach_inner"] for row in finals), default=None),
            "reach_outer_min": min((row["reach_outer"] for row in finals), default=None),
            "reach_outer_max": max((row["reach_outer"] for row in finals), default=None),
            "dominant_fraction_inner_min": min(
                (row["dominant_fraction_inner"] for row in finals), default=None),
        },
        "compression_relation": {
            "question": (
                "Douglas asked whether a strongly compressed field is a tiny, strongly "
                "localised one rather than a long-reach one. Does the dominant force's "
                "reach shorten as the throat compression steepens?"
            ),
            "throat_slope_vs_reach_inner_pearson_r": correlation,
            "case_pairs": [{"case": c["case"],
                            "throat_slope_inner": c["final"]["throat_slope_inner"],
                            "reach_inner": c["final"]["reach_inner"]}
                           for c in cases if c["final"]],
            "reading": (
                "A strongly negative correlation means steeper compression shortens the "
                "dominant force's reach, which is the 'strong but local' reading. Six "
                "points, one family, no error bar: it is a direction, not a law."
            ),
        },
        "cases": cases,
    }

    RECORD.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {RECORD.relative_to(LAB)}")
    print(f"  cases           : {record['coverage']['cases']}")
    print(f"  stations        : {record['coverage']['total_stations']}")
    print(f"  verdicts        : {verdicts}")
    print(f"  reach_inner     : {record['summary']['reach_inner_min']:.3f}"
          f" .. {record['summary']['reach_inner_max']:.3f}")
    print(f"  reach_outer     : {record['summary']['reach_outer_min']:.3f}"
          f" .. {record['summary']['reach_outer_max']:.3f}")
    print(f"  radial dominant : {record['summary']['dominant_fraction_inner_min']:.4f} of nodes (min)")
    if correlation is not None:
        print(f"  slope vs reach r: {correlation:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())