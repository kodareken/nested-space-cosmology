"""Drive the two-force gradient balance across the recorded RN neutral family.

Reads every station already written by ``nsc-rn-neutral-family-v1`` and
reduces each one to a finite two-force balance record. No time step is
taken, no station is re-run, and no existing file is rewritten.

The family is incomplete on disk: three cases stopped on the CPU
projection and two stopped early, and no case wrote a ``summary.json``.
This driver therefore reports coverage honestly per case instead of
describing a family that does not exist.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_gradient_bubble.py
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

from recursive_horizons import nsc_rn_gradient_bubble as bubble  # noqa: E402

FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"
RECORD = LAB / "results/development/nsc-rn-gradient-bubble-v1.json"
CASE_PATTERN = re.compile(r"^mu-(?P<mu>[0-9.]+)-eps-(?P<eps>[0-9.]+)$")
STATION_PATTERN = re.compile(r"^station-(?P<t>[0-9.]+)\.npz$")


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


def _station_time(path: Path) -> float:
    match = STATION_PATTERN.match(path.name)
    return float(match.group("t")) if match else float("inf")


def _reduce_case(case_dir: Path) -> Dict[str, Any]:
    match = CASE_PATTERN.match(case_dir.name)
    stations = sorted((p for p in case_dir.glob("station-*.npz") if STATION_PATTERN.match(p.name)),
                      key=_station_time)
    failure_path = case_dir / "station-failure.json"
    failure = json.loads(failure_path.read_text()) if failure_path.exists() else None

    rows: List[Dict[str, Any]] = []
    points = None
    for station in stations:
        with np.load(station) as payload:
            if points is None:
                points = int(payload["radius"].shape[0])
            assessment = bubble.gradient_bubble_assessment(payload, {"R4": payload["R4"]})
            clock = float(payload["coupled_t"][0]) if "coupled_t" in payload.files else _station_time(station)
        profile = assessment["profile"]
        rows.append({
            "station": station.stem,
            "t_over_rm": clock,
            "verdict": assessment["verdict"],
            "dead_end_flags": assessment["dead_end_flags"],
            "balance_at_throat": profile["balance_at_throat"],
            "balance_min": profile["balance_min"],
            "balance_at_mid": profile["balance_at_mid"],
            "r_decouple": profile["r_decouple"],
            "throat_slope_inner": profile["throat_slope_inner"],
            "force_ratio_integrated": profile["force_ratio_integrated"],
            "curvature_force_ratio": assessment["curvature_force_ratio"],
        })

    verdicts = sorted({row["verdict"] for row in rows})
    flags = sorted({flag for row in rows for flag in row["dead_end_flags"]})
    last = rows[-1] if rows else None
    return {
        "case": case_dir.name,
        "mass_over_rm": float(match.group("mu")) if match else None,
        "epsilon": float(match.group("eps")) if match else None,
        "points": points,
        "station_count": len(rows),
        "reached_t_over_rm": last["t_over_rm"] if last else None,
        "verdicts": verdicts,
        "dead_end_flags_observed": flags,
        "has_summary_json": (case_dir / "summary.json").exists(),
        "failure_blocker": failure.get("blocker") if failure else None,
        "failure_clock_t": failure.get("clock_t") if failure else None,
        "stations": rows,
    }


def main() -> int:
    if not FAMILY.is_dir():
        print(f"missing family directory: {FAMILY}", file=sys.stderr)
        return 1
    if RECORD.exists():
        print(f"refusing to overwrite sealed record: {RECORD}", file=sys.stderr)
        return 1

    cases = [_reduce_case(p) for p in sorted(FAMILY.iterdir()) if p.is_dir() and CASE_PATTERN.match(p.name)]

    all_flags = sorted({flag for case in cases for flag in case["dead_end_flags_observed"]})
    verdicts = sorted({v for case in cases for v in case["verdicts"]})
    balanced = [c for c in cases if "two_force_balance" in c["verdicts"]]
    throats = [c["stations"][-1]["balance_at_throat"] for c in cases if c["stations"]]
    ratios = [c["stations"][-1]["force_ratio_integrated"] for c in cases if c["stations"]]
    decouple = [c["stations"][-1]["r_decouple"] for c in cases if c["stations"]]

    record = {
        "schema": "NSC-RN-GRADIENT-BUBBLE-v1",
        "api_version": bubble.API_VERSION,
        "creation_only": True,
        "producing_commit": _producing_commit(),
        "family_source": FAMILY.name,
        "family_source_sha256": _sha256_file(FAMILY / "prepare.json"),
        "domain": {
            "chart": "areal radius",
            "r_inner": "fixed inner mass at grid[0]",
            "r_outer": 32.0,
            "observer": "areal chart of the recorded coupled arm",
            "state": "recorded station payload, no time step taken",
            "compared_quantity": "balance of inherited radial and generated angular pressure",
        },
        "scope_note": (
            "A construction check on recorded finite states. It does not measure the "
            "actual universe, and a bounded balance here is not evidence about cosmology."
        ),
        "coverage": {
            "cases": len(cases),
            "cases_with_summary_json": sum(1 for c in cases if c["has_summary_json"]),
            "cases_showing_two_force_balance": len(balanced),
            "total_stations": sum(c["station_count"] for c in cases),
            "family_complete": all(c["has_summary_json"] for c in cases),
        },
        "verdicts": verdicts,
        "dead_end_flags_observed": all_flags,
        "summary": {
            "balance_at_throat_min": min(throats) if throats else None,
            "balance_at_throat_max": max(throats) if throats else None,
            "integrated_outer_over_inner_min": min(ratios) if ratios else None,
            "integrated_outer_over_inner_max": max(ratios) if ratios else None,
            "r_decouple_min": min(decouple) if decouple else None,
            "r_decouple_max": max(decouple) if decouple else None,
        },
        "limits": {
            "throat_slope_limit": bubble.THROAT_SLOPE_LIMIT,
            "single_force_balance": bubble.SINGLE_FORCE_BALANCE,
            "exterior_decoupled": bubble.EXTERIOR_DECOUPLED,
            "curvature_force_ratio": bubble.CURVATURE_FORCE_RATIO,
            "throat_fraction": bubble.THROAT_FRACTION,
        },
        "cases": cases,
    }

    RECORD.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {RECORD.relative_to(LAB)}")
    print(f"  cases            : {record['coverage']['cases']}")
    print(f"  stations reduced : {record['coverage']['total_stations']}")
    print(f"  family complete  : {record['coverage']['family_complete']}")
    print(f"  verdicts         : {verdicts}")
    print(f"  flags observed   : {all_flags}")
    print(f"  balance@throat   : {record['summary']['balance_at_throat_min']:.4f}"
          f" .. {record['summary']['balance_at_throat_max']:.4f}")
    print(f"  F_out/F_in       : {record['summary']['integrated_outer_over_inner_min']:.4f}"
          f" .. {record['summary']['integrated_outer_over_inner_max']:.4f}")
    print(f"  r_decouple       : {record['summary']['r_decouple_min']:.3f}"
          f" .. {record['summary']['r_decouple_max']:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())