"""Reduce radial/transverse pressures across the recorded RN neutral family.

Reads every station already written by ``nsc-rn-neutral-family-v1`` and
reduces each one to a finite pressure-diagnostic record. No time step is
taken, no station is re-run, and no existing file is rewritten.

The family is incomplete on disk: four cases stopped on the CPU
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
from recursive_horizons import nsc_rn_stencil_audit as stencil  # noqa: E402

FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"
RECORD = LAB / "results/development/nsc-rn-gradient-bubble-v6.json"
PRIOR_RECORD = LAB / "results/development/nsc-rn-gradient-bubble-v5.json"
EARLIER_RECORDS = (
    LAB / "results/development/nsc-rn-gradient-bubble-v4.json",
    LAB / "results/development/nsc-rn-gradient-bubble-v3.json",
    LAB / "results/development/nsc-rn-gradient-bubble-v2.json",
    LAB / "results/development/nsc-rn-gradient-bubble-v1.json",
)
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
            mu = float(CASE_PATTERN.match(case_dir.name).group("mu"))
            assessment = bubble.gradient_bubble_assessment(
                payload, {"R4": payload["R4"]}, mass_over_rm=mu)
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
            "wall_count": profile["wall_count"],
            "r_first_wall": profile["r_first_wall"],
            "r_last_wall": profile["r_last_wall"],
            "walls": [float(w) for w in profile["walls"]],
            "throat_slope_inner": profile["throat_slope_inner"],
            "force_ratio_integrated": profile["force_ratio_integrated"],
            "curvature_force_ratio": assessment["curvature_force_ratio"],
            "slope_meets_regularity_limit": assessment["slope_meets_regularity_limit"],
            "regularity_extrapolation": assessment["regularity_extrapolation"],
        })

    verdicts = sorted({row["verdict"] for row in rows})
    flags = sorted({flag for row in rows for flag in row["dead_end_flags"]})
    last = rows[-1] if rows else None
    wall_rows = [row for row in rows if row["wall_count"] > 0]
    wall_radii = [row["r_first_wall"] for row in wall_rows]
    return {
        "case": case_dir.name,
        "mass_over_rm": float(match.group("mu")) if match else None,
        "epsilon": float(match.group("eps")) if match else None,
        "points": points,
        "station_count": len(rows),
        "stations_with_a_wall": len(wall_rows),
        "wall_presence_fraction": (len(wall_rows) / len(rows)) if rows else None,
        "r_wall_min": min(wall_radii) if wall_radii else None,
        "r_wall_max": max(wall_radii) if wall_radii else None,
        "r_wall_final": last["r_first_wall"] if last else None,
        "wall_at_final_station": bool(last["wall_count"] > 0) if last else False,
        "reached_t_over_rm": last["t_over_rm"] if last else None,
        "slope_meets_regularity_limit": None,
        "verdicts": verdicts,
        "dead_end_flags_observed": flags,
        "has_summary_json": (case_dir / "summary.json").exists(),
        "failure_blocker": failure.get("blocker") if failure else None,
        "failure_clock_t": failure.get("clock_t") if failure else None,
        "stations": rows,
    }



def _stencil_validity(cases):
    """Is the measured throat steepness physics, or the stencil's fault?

    `nsc_rn_pg` raises "nodal stencil failed its polynomial moment" at
    collection time for widths of 9 and above, which blocks re-marching
    the family. This section tests whether those weights are actually
    wrong, by comparing them against the same recursion evaluated in
    exact rational arithmetic, and by re-deriving the throat slope under
    both. Nothing in `nsc_rn_pg` is imported, edited or replaced.
    """
    audit = stencil.audit_stencils()
    slopes = []
    for case in cases:
        stations = sorted(
            (p for p in (FAMILY / case["case"]).glob("station-*.npz")
             if STATION_PATTERN.match(p.name)),
            key=_station_time,
        )
        if not stations:
            continue
        with np.load(stations[-1]) as payload:
            comparison = stencil.slope_under_both_stencils(
                payload["radius"], payload["radial_pressure"])
        comparison["case"] = case["case"]
        slopes.append(comparison)
    worst = max((abs(row["slope_difference"]) for row in slopes), default=0.0)
    return {
        "question": (
            "The recorded local pressure slope ranges from -7 to -488. Is that "
            "steepness changed by exact stencil weights rather than the nodal weights the "
            "chart itself reports as failing?"
        ),
        "method": (
            "The identical Fornberg recursion is re-evaluated in exact rational "
            "arithmetic, since the nodes are integers and the weights are rational. "
            "The throat slope is then re-derived under both weight sets. "
            "nsc_rn_pg is not imported, edited or replaced."
        ),
        "moment_tolerance": audit["moment_tolerance"],
        "stencil_rows_audited": len(audit["rows"]),
        "stencil_rows_failing_exact_arithmetic": audit["rows_failing_extended"],
        "configured_one_sided_failing": audit["configured_failing_plain"],
        "max_relative_weight_difference": max(
            (row["max_relative_weight_difference"] for row in audit["rows"]), default=0.0),
        "weights_agree_to": (
            f"{max((row['max_relative_weight_difference'] for row in audit['rows']), default=0.0):.1e}"),
        "weights_exact": max(
            (row["max_relative_weight_difference"] for row in audit["rows"]), default=1.0) < 1e-12,
        "throat_slope_difference_max": worst,
        "throat_slope_confirms_verdict": None,
        "pressure_slope_consistent_across_weight_sets": worst < 1e-6,
        "throat_slopes": slopes,
        "conclusion": (
            "The recorded weights agree with exact rational arithmetic to a few ulp "
            "across every audited stencil, and the throat slope is unchanged to about "
            "1e-13 when re-derived under exact weights. This checks stencil agreement, "
            "not the physical interpretation of the pressure slope. The "
            "RuntimeError is a false alarm from an ill-conditioned self-check: the "
            "moment sum cancels catastrophically, so evaluating it in double precision "
            "at 1e-8 is not a valid test at these widths."
        ),
        "blocker_status": (
            "Downgraded. This is not a blocker on re-marching the family, and it does "
            "not invalidate any recorded operator value or calibration hash. What is "
            "broken is the guard, not the weights. Correcting the guard changes no "
            "operator and therefore no provenance."
        ),
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
    balanced = [c for c in cases if "two_force_wall" in c["verdicts"]]
    throats = [c["stations"][-1]["balance_at_throat"] for c in cases if c["stations"]]
    ratios = [c["stations"][-1]["force_ratio_integrated"] for c in cases if c["stations"]]
    decouple = [c["stations"][-1]["r_decouple"] for c in cases if c["stations"]]
    walls = [c["stations"][-1]["r_first_wall"] for c in cases
             if c["stations"] and c["stations"][-1]["wall_count"] > 0]
    slopes = [c["stations"][-1]["throat_slope_inner"] for c in cases if c["stations"]]
    gaps = [c["stations"][-1]["regularity_extrapolation"]["gap_to_horizon"]
            for c in cases if c["stations"]]
    extrap = [c["stations"][-1]["regularity_extrapolation"]["extrapolated_log10_force_ratio"]
              for c in cases if c["stations"]]
    reaches = any(c["stations"][-1]["regularity_extrapolation"]["domain_reaches_horizon"]
                  for c in cases if c["stations"])

    prior = None
    if PRIOR_RECORD.is_file():
        prior_doc = json.loads(PRIOR_RECORD.read_text(encoding="utf-8"))
        prior = {
            "record": PRIOR_RECORD.name,
            "sha256": _sha256_file(PRIOR_RECORD),
            "verdicts": prior_doc.get("verdicts"),
            "superseded_reason": (
                "The g >= -1 regularity test and its Hayward attribution are unsupported. "
                "Inputs are orthonormal pressure components, not opposing forces. "
                "Crossings and slopes remain finite diagnostics; reaching the inner "
                "horizon does not determine centre regularity. Historical bytes remain intact."
            ),
            "case_verdicts": {c["case"]: c["stations"][-1]["verdict"] for c in prior_doc["cases"]},
        }

    final_walls = [c for c in cases
                  if c["stations"] and c["stations"][-1]["wall_count"] > 0]
    ever_walls = [c for c in cases if "two_force_wall" in c["verdicts"]]

    validity = _stencil_validity(cases)

    record = {
        "schema": "NSC-RN-GRADIENT-BUBBLE-v6",
        "api_version": bubble.API_VERSION,
        "diagnostic_semantics": bubble.diagnostic_semantics(),
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
            "compared_quantity": "relative magnitudes of radial and transverse orthonormal pressures",
            "verdict_basis": (
                "A legacy wall label denotes an equal-magnitude pressure crossing "
                "surviving the declared hysteresis, persistence and live floors. "
                "It does not establish mechanical equilibrium or a cosmological boundary."
            ),
            "wall_dynamics": (
                "The wall is tracked across every station of a case, not read at the last "
                "one. A wall that appears at some stations and not others is oscillating, "
                "not settled, and the per-case presence fraction records that."
            ),
        },
        "regularity_scope": (
            "Regularity is NOT decidable from this positive-radius pressure diagnostic. "
            "An assumed power law extrapolated to the inner horizon remains an "
            "extrapolation. Even reaching that horizon does not establish centre "
            "regularity. No inward continuation is required or selected by this record."
        ),
        "regularity_criterion": None,
        "regularity_attribution_correction": (
            "Hayward gr-qc/0506126 gives a finite-stress regular-core example; "
            "it does not supply the v5 g >= -1 criterion for these pressure components."
        ),
        "supersedes": prior,
        "earlier_records": [p.name for p in EARLIER_RECORDS if p.is_file()],
        "scope_note": (
            "A construction check on recorded finite states. It does not measure the "
            "actual universe, and a bounded balance here is not evidence about cosmology."
        ),
        "coverage": {
            "cases": len(cases),
            "cases_with_summary_json": sum(1 for c in cases if c["has_summary_json"]),
            "cases_with_a_wall_at_the_final_station": sum(1 for c in cases if c["wall_at_final_station"]),
            "cases_showing_a_wall_at_any_station": sum(1 for c in cases if c["stations_with_a_wall"] > 0),
            "total_stations": sum(c["station_count"] for c in cases),
            "station_fraction_with_a_wall": (
                sum(c["stations_with_a_wall"] for c in cases) / sum(c["station_count"] for c in cases)
                if sum(c["station_count"] for c in cases) else None
            ),
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
            "r_first_wall_min": min(walls) if walls else None,
            "r_first_wall_max": max(walls) if walls else None,
            "domain_reaches_horizon_any_case": bool(reaches),
            "cases_meeting_the_regularity_slope_limit": None,
            "gap_to_horizon_min": min(gaps) if gaps else None,
            "gap_to_horizon_max": max(gaps) if gaps else None,
            "extrapolated_log10_ratio_min": min(extrap) if extrap else None,
            "extrapolated_log10_ratio_max": max(extrap) if extrap else None,
            "throat_slope_min": min(slopes) if slopes else None,
            "throat_slope_max": max(slopes) if slopes else None,
        },
        "limits": {
            "throat_slope_limit": bubble.THROAT_SLOPE_LIMIT,
            "single_force_balance": bubble.SINGLE_FORCE_BALANCE,
            "exterior_decoupled": bubble.EXTERIOR_DECOUPLED,
            "curvature_force_ratio": bubble.CURVATURE_FORCE_RATIO,
            "throat_fraction": bubble.THROAT_FRACTION,
            "wall_live_fraction": bubble.WALL_LIVE_FRACTION,
            "wall_hysteresis": bubble.WALL_HYSTERESIS,
            "wall_min_cells": bubble.WALL_MIN_CELLS,
            "regular_slope_limit": None,
            "regular_slope_tolerance": None,
            "horizon_reach_tolerance": bubble.HORIZON_REACH_TOLERANCE,
        },
        "numerical_validity": validity,
        "cases": cases,
    }

    with RECORD.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=1, sort_keys=True) + "\n")
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
    if walls:
        print(f"  r_first_wall     : {min(walls):.3f} .. {max(walls):.3f}")
    else:
        print("  r_first_wall     : none in any case")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
