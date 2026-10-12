"""Write the two-ends room-quality record across the recorded RN family.

Reads every recorded station and reduces its two ends (the angular
ground ``N/r`` and the CFL resolution ``(N+|beta|)/spacing``) to one
dimensionless per-room number. Creation only.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_two_ends.py
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_rn_two_ends as ends  # noqa: E402

FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"
RECORD = LAB / "results/development/nsc-rn-two-ends-v1.json"
CASE_PATTERN = re.compile(r"^mu-(?P<mu>[0-9.]+)-eps-(?P<eps>[0-9.]+)$")
STATION_PATTERN = re.compile(r"^station-(?P<t>[0-9.]+)\.npz$")


def _family_digest(family: Path) -> tuple[str, str]:
    """Digest of the recorded inputs a reduction stands on.

    The neutral family keeps one root ``prepare.json``; a campaign layout keeps
    one per case. The digest folds whichever layout is present, in a fixed
    order, and rewrites nothing.
    """
    root = family / "prepare.json"
    if root.is_file():
        return "prepare.json", _sha256_file(root)
    parts = sorted(family.glob("*/prepare.json"))
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
        assessment = ends.station_two_ends(payload, mass_over_rm=mu)
    return {
        "station": station.stem,
        "t_over_rm": assessment["t_over_rm"],
        "trap_radius": assessment["trap_radius"],
        "room_quality_inner": assessment["room_quality_inner"],
        "room_quality_at_r_plus": assessment["room_quality_at_r_plus"],
        "room_quality_mid": assessment["room_quality_mid"],
        "room_quality_outer": assessment["room_quality_outer"],
        "quality_ratio_inner_over_outer": assessment["quality_ratio_inner_over_outer"],
        "quality_ratio_inner_over_at_r_plus": assessment["quality_ratio_inner_over_at_r_plus"],
        "dimensionless_invariant_is_per_room": assessment["dimensionless_invariant_is_per_room"],
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

    cases = []
    for case_dir in sorted(family.iterdir()):
        match = CASE_PATTERN.match(case_dir.name)
        if not case_dir.is_dir() or match is None:
            continue
        mu = float(match.group("mu"))
        stations = sorted(
            (p for p in case_dir.glob("station-*.npz") if STATION_PATTERN.match(p.name)),
            key=_station_time,
        )
        rows = [_station_row(s, mu) for s in stations]
        last = rows[-1]
        cases.append({
            "case": case_dir.name,
            "mass_over_rm": mu,
            "epsilon": float(match.group("eps")),
            "station_count": len(rows),
            "has_summary_json": (case_dir / "summary.json").exists(),
            "room_quality_inner_last": last["room_quality_inner"],
            "room_quality_outer_last": last["room_quality_outer"],
            "quality_ratio_inner_over_outer_last": last["quality_ratio_inner_over_outer"],
            "stations": rows,
        })
    cases.sort(key=lambda c: (c["mass_over_rm"], c["epsilon"]))
    source_kind, source_digest = _family_digest(family)

    record = {
        "schema": "NSC-RN-TWO-ENDS-v1",
        "api_version": ends.API_VERSION,
        "creation_only": True,
        "producing_commit": _producing_commit(),
        "family_source": family.name,
        "family_source_kind": source_kind,
        "family_source_sha256": source_digest,
        "definitions": {
            "omega_zero": "N/r, the angular ground of the room",
            "omega_cutoff": "(N + |beta|)/spacing, the CFL resolution of the room",
            "room_quality": "omega_zero / omega_cutoff, dimensionless per-room number",
            "scale_invariant_expression": "N^2 spacing / (r (N + |beta|))",
        },
        "reading": (
            "The bottom and top ends of a room are both set by the local space. "
            "For the recorded states the room quality is many orders of magnitude "
            "larger at the inner edge than at the outer boundary in this flat-slice "
            "gauge, so the recorded states hold a sharper room at the bottom than at "
            "the top. This is a reading of measured numbers, not a derivation of "
            "Planck's constant or of absolute zero."
        ),
        "coverage": {
            "cases": len(cases),
            "total_stations": sum(c["station_count"] for c in cases),
            "family_complete": all(c["has_summary_json"] for c in cases),
        },
        "summary": {
            "quality_ratio_inner_over_outer_last": {
                c["case"]: c["quality_ratio_inner_over_outer_last"] for c in cases
            },
            "room_quality_inner_last": {
                c["case"]: c["room_quality_inner_last"] for c in cases
            },
        },
        "cases": cases,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  cases     : {record['coverage']['cases']}")
    print(f"  stations  : {record['coverage']['total_stations']}")
    for case in cases:
        print(
            f"  {case['case']:22s} Q_inner={case['room_quality_inner_last']:.6e} "
            f"Q_ratio_inner/outer={case['quality_ratio_inner_over_outer_last']:.4e}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())