"""Drive the cone map across the recorded RN neutral family.

Reads every station already written by ``nsc-rn-neutral-family-v1`` and
reduces each one to its cone facts: the local speed (one), the chart
speeds, the ray frequency dilution along the ride, the chart's own
trapped boundary, and the redshift at the inner edge. No time step is
taken, no station is re-run, and no existing file is rewritten.

The family is incomplete on disk: no case wrote a ``summary.json``.
Coverage is therefore reported honestly per case, exactly like the
gradient-bubble driver, instead of describing a family that does not
exist.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_cone_map.py

The same reduction runs over any recorded family layout (each case directory
holding ``mu-*-eps-*`` with ``station-*.npz``), so the inward-march campaign can
be reduced on the domain that reaches the inner horizon:

    python scripts/lab.py scripts/derive_nsc_rn_cone_map.py \\
        --family results/development/nsc-rn-inward-march-v1 \\
        --record results/development/nsc-rn-cone-map-v2.json
"""
from __future__ import annotations

import argparse
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

from recursive_horizons import nsc_rn_cone_map as cone  # noqa: E402

FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"
RECORD = LAB / "results/development/nsc-rn-cone-map-v1.json"
CASE_PATTERN = re.compile(r"^mu-(?P<mu>[0-9.]+)-eps-(?P<eps>[0-9.]+)$")
STATION_PATTERN = re.compile(r"^station-(?P<t>[0-9.]+)\.npz$")


def _family_digest(family: Path) -> tuple[str, str]:
    """Digest of the recorded inputs a reduction stands on.

    The neutral family keeps one root ``prepare.json``. A campaign layout keeps
    one per case instead, so the digest folds the case preparation markers in
    a fixed order. No file is rewritten either way.
    """
    root = family / "prepare.json"
    if root.is_file():
        return "prepare.json", _sha256_file(root)
    parts = sorted(p for p in family.glob("*/prepare.json"))
    if not parts:
        raise FileNotFoundError(f"no prepare.json under {family}")
    digest = hashlib.sha256()
    for path in parts:
        digest.update(path.parent.name.encode("utf-8"))
        digest.update(_sha256_file(path).encode("utf-8"))
    return "case prepare.json set", digest.hexdigest()


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


def _station_row(station: Path, mu: float) -> Dict[str, Any]:
    with np.load(station) as payload:
        assessment = cone.station_assessment(payload, mass_over_rm=mu)
    profile = assessment["profile"]
    return {
        "station": station.stem,
        "t_over_rm": assessment["t_over_rm"],
        "verdict": assessment["verdict"],
        "flags": assessment["flags"],
        "mass_residual": assessment["mass_residual"],
        "lapse_at_inner": profile["lapse_at_inner"],
        "lapse_min": profile["lapse_min"],
        "redshift_at_excision": profile["redshift_at_excision"],
        "chart_speeds_at_inner": profile["chart_speeds_at_inner"],
        "chart_speeds_at_outer": profile["chart_speeds_at_outer"],
        "trap_radius": profile["trap_radius"],
        "trap_minus_reference": assessment["horizon_comparison"]["trap_minus_reference"],
        "lapse_log_slope_at_throat": profile["lapse_log_slope_at_throat"],
        "crossing_coordinate_time_full": float(
            profile["crossing_coordinate_time"][0]
        ),
        "light_delay_full": float(profile["light_delay_vs_flat"][0]),
    }


def _array_list(values: Any) -> List[float]:
    return [float(item) for item in np.asarray(values, dtype=np.float64)]


def _last_station_arrays(station: Path, mu: float) -> Dict[str, Any]:
    """The final station's profile arrays, so the integrals are re-checkable."""
    with np.load(station) as payload:
        assessment = cone.station_assessment(payload, mass_over_rm=mu)
    profile = assessment["profile"]
    return {
        "station": station.stem,
        "lapse": _array_list(profile["lapse"]),
        "beta": _array_list(profile["beta"]),
        "v_out": _array_list(profile["v_out"]),
        "v_in": _array_list(profile["v_in"]),
        "crossing_coordinate_time": _array_list(profile["crossing_coordinate_time"]),
        "light_delay_vs_flat": _array_list(profile["light_delay_vs_flat"]),
    }


def _reduce_case(case_dir: Path) -> Dict[str, Any]:
    match = CASE_PATTERN.match(case_dir.name)
    mu = float(match.group("mu"))
    epsilon = float(match.group("eps"))
    stations = sorted(
        (p for p in case_dir.glob("station-*.npz") if STATION_PATTERN.match(p.name)),
        key=_station_time,
    )
    rows = [_station_row(s, mu) for s in stations]
    failure_path = case_dir / "station-failure.json"
    failure = json.loads(failure_path.read_text()) if failure_path.exists() else None
    last = rows[-1] if rows else None
    return {
        "case": case_dir.name,
        "mass_over_rm": mu,
        "epsilon": epsilon,
        "station_count": len(rows),
        "reached_t_over_rm": last["t_over_rm"] if last else None,
        "stations": rows,
        "last_station_arrays": _last_station_arrays(stations[-1], mu) if stations else None,
        "redshift_at_excision_first": rows[0]["redshift_at_excision"] if rows else None,
        "redshift_at_excision_last": last["redshift_at_excision"] if last else None,
        "trap_radius_first": rows[0]["trap_radius"] if rows else None,
        "trap_radius_last": last["trap_radius"] if last else None,
        "trap_minus_reference_last": last["trap_minus_reference"] if last else None,
        "lapse_log_slope_last": last["lapse_log_slope_at_throat"] if last else None,
        "has_summary_json": (case_dir / "summary.json").exists(),
        "failure_blocker": failure.get("blocker") if failure else None,
        "failure_clock_t": failure.get("clock_t") if failure else None,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", default=str(FAMILY.relative_to(LAB)))
    parser.add_argument("--record", default=str(RECORD.relative_to(LAB)))
    args = parser.parse_args()
    family = (LAB / args.family).resolve()
    record_path = (LAB / args.record).resolve()

    if not family.is_dir():
        print(f"missing family directory: {family}", file=sys.stderr)
        return 1
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    cases = [
        _reduce_case(p) for p in sorted(family.iterdir())
        if p.is_dir() and CASE_PATTERN.match(p.name)
    ]
    if not cases:
        print(f"no case directories under {family}", file=sys.stderr)
        return 1
    cases.sort(key=lambda c: (c["mass_over_rm"], c["epsilon"]))
    source_kind, source_digest = _family_digest(family)

    flags = sorted({flag for c in cases for r in c["stations"] for flag in r["flags"]})
    red_first = {c["case"]: c["redshift_at_excision_first"] for c in cases}
    red_last = {c["case"]: c["redshift_at_excision_last"] for c in cases}
    trap_first = {c["case"]: c["trap_radius_first"] for c in cases}
    trap_last = {c["case"]: c["trap_radius_last"] for c in cases}
    residues = [r["mass_residual"] for c in cases for r in c["stations"]]

    record = {
        "schema": "NSC-RN-CONE-MAP-v1",
        "api_version": cone.API_VERSION,
        "creation_only": True,
        "producing_commit": _producing_commit(),
        "family_source": family.name,
        "family_source_kind": source_kind,
        "family_source_sha256": source_digest,
        "domain": {
            "chart": "areal radius",
            "observer": "chart of the recorded coupled arm",
            "state": "recorded station payload, no time step taken",
            "compared_quantity": (
                "light cone of the state's own metric: local speed, chart "
                "speeds, ray frequency dilution, trapped boundary, redshift "
                "at the inner edge"
            ),
            "metric_reconstruction": (
                "lapse and shift are reconstructed with "
                "nsc_rn_pg.reconstruct_geometry at weight 1 from the stored "
                "coupled_phi, occupations and coupled_inner_mass, the same "
                "route the march used; a station whose mass does not match "
                "the stored coupled_mass below 1e-9 is refused"
            ),
        },
        "reading": (
            "Light always rides its own space at unit speed. The chart value of "
            "the cone is -beta +/- N and changes with radius because the "
            "recorded lapse and shift carry the room's stress; the ray "
            "frequency dilution 1/N(r) is how much the value of light has "
            "changed after riding the space. This is a reading of measured "
            "numbers on evaluated grids, not a law and not a cosmology claim."
        ),
        "coverage": {
            "cases": len(cases),
            "cases_with_summary_json": sum(1 for c in cases if c["has_summary_json"]),
            "total_stations": sum(c["station_count"] for c in cases),
            "family_complete": all(c["has_summary_json"] for c in cases),
        },
        "flags_observed": flags,
        "mass_reconstruction_residual_max": max(residues) if residues else None,
        "mass_reconstruction_tolerance": cone.MASS_RECONSTRUCTION_TOLERANCE,
        "summary": {
            "redshift_at_excision_first_station": red_first,
            "redshift_at_excision_last_station": red_last,
            "trap_radius_first_station": trap_first,
            "trap_radius_last_station": trap_last,
            "trap_minus_reference_last_station": {
                c["case"]: c["trap_minus_reference_last"] for c in cases
            },
            "lapse_log_slope_last_station": {
                c["case"]: c["lapse_log_slope_last"] for c in cases
            },
        },
        "cases": cases,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  cases          : {record['coverage']['cases']}")
    print(f"  stations       : {record['coverage']['total_stations']}")
    print(f"  family complete: {record['coverage']['family_complete']}")
    print(f"  flags observed : {flags}")
    print(f"  mass residual  : {record['mass_reconstruction_residual_max']:.3e}")
    for case in cases:
        print(
            f"  {case['case']:22s} stations={case['station_count']:3d} "
            f"redshift@inner={case['redshift_at_excision_last']:.6f} "
            f"trap-r_plus={case['trap_minus_reference_last']:+.5f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())