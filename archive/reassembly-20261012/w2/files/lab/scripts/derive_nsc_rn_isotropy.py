"""Reduce the recorded states to an isotropy trajectory and register the test.

No time step is taken. Every station of every supplied layout is read, its
two-force balance profile recomputed with the owning module, and reduced to the
single departure-from-isotropic scalar. The registered pass/fail of the
not-yet-run relaxation test is written into the record, so the condition is on
disk before the run exists.

    python scripts/lab.py scripts/derive_nsc_rn_isotropy.py \\
        --family results/development/nsc-rn-neutral-family-v1 \\
        --family results/development/nsc-rn-inward-march-v1 \\
        --record results/development/nsc-rn-isotropy-v1.json
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
from recursive_horizons import nsc_rn_gradient_bubble as bubble  # noqa: E402
from recursive_horizons import nsc_rn_isotropy as isotropy  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402

FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"
RECORD = LAB / "results/development/nsc-rn-isotropy-v1.json"
CASE_PATTERN = re.compile(r"^mu-(?P<mu>[0-9.]+)-eps-(?P<eps>[0-9.]+)$")
STATION_PATTERN = re.compile(r"^station-(?P<t>[0-9.]+)\.npz$")

ARMS = {
    "coupled": (
        "the recorded coupled arm: pressures as written by the march, with the "
        "coupled lapse and shift"
    ),
    "fixed": (
        "the frozen-metric control arm: pressures recomputed from the stored "
        "fixed_phi on its own frozen lapse and shift, so the no-feedback "
        "control carries the same statistic"
    ),
}


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


def _station_balance(payload, arm: str):
    """Balance profile and clock of one station, for the requested arm."""
    if arm == "coupled":
        profile = bubble.gradient_bubble_profile(payload)
        return np.asarray(profile["balance"]), float(np.asarray(payload["coupled_t"])[0])
    for key in ("fixed_phi", "frozen_lapse", "frozen_shift", "occupations"):
        if key not in payload.files:
            raise KeyError(f"station payload has no {key}; the fixed arm needs it")
    grid = cone.station_grid(dict(payload))
    evaluated = source.source_evaluate(
        np.asarray(payload["fixed_phi"]), np.asarray(payload["occupations"]),
        np.asarray(payload["frozen_lapse"]), np.asarray(payload["frozen_shift"]),
        np.ones_like(grid["radius"]), grid["radius"], grid["weights"], grid["derivative"],
    )
    radial = np.asarray(evaluated["physical"]["radial_pressure"], dtype=float)
    angular = np.asarray(evaluated["physical"]["angular_pressure"], dtype=float)
    balance = isotropy.balance_from_pressures(radial, angular)
    return balance, float(np.asarray(payload["fixed_t"])[0])


def _reduce_case(case_dir: Path, arm: str) -> Dict[str, Any]:
    match = CASE_PATTERN.match(case_dir.name)
    mu = float(match.group("mu"))
    stations = sorted(
        (p for p in case_dir.glob("station-*.npz") if STATION_PATTERN.match(p.name)),
        key=_station_time,
    )
    rows: List[Dict[str, Any]] = []
    for station in stations:
        with np.load(station) as payload:
            balance, clock = _station_balance(payload, arm)
        rows.append(isotropy.station_row(station.stem, clock, np.asarray(balance)))
    verdict = isotropy.relaxation_verdict(rows) if rows else None
    return {
        "case": case_dir.name,
        "mass_over_rm": mu,
        "epsilon": float(match.group("eps")),
        "station_count": len(rows),
        "has_summary_json": (case_dir / "summary.json").exists(),
        "isotropy_index_first": rows[0]["isotropy_index"] if rows else None,
        "isotropy_index_last": rows[-1]["isotropy_index"] if rows else None,
        "isotropy_index_max": max((r["isotropy_index"] for r in rows), default=None),
        "relaxation_registration": verdict,
        "stations": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", action="append", default=None,
                        help="layout directory; repeat for several layouts")
    parser.add_argument("--record", default=str(RECORD.relative_to(LAB)))
    parser.add_argument("--arm", choices=sorted(ARMS), default="coupled")
    args = parser.parse_args()
    families = [LAB / f for f in (args.family or [str(FAMILY.relative_to(LAB))])]
    record_path = (LAB / args.record).resolve()

    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    cases: List[Dict[str, Any]] = []
    sources: List[Dict[str, Any]] = []
    for family in families:
        family = family.resolve()
        if not family.is_dir():
            print(f"missing family directory: {family}", file=sys.stderr)
            return 1
        found = [
            _reduce_case(p, args.arm) for p in sorted(family.iterdir())
            if p.is_dir() and CASE_PATTERN.match(p.name)
        ]
        for case in found:
            case["layout"] = family.name
        cases.extend(found)
        root = family / "prepare.json"
        if root.is_file():
            sources.append({"layout": family.name, "kind": "prepare.json",
                            "sha256": _sha256_file(root)})
        else:
            digest = hashlib.sha256()
            for path in sorted(family.glob("*/prepare.json")):
                digest.update(path.parent.name.encode("utf-8"))
                digest.update(_sha256_file(path).encode("utf-8"))
            sources.append({"layout": family.name, "kind": "case prepare.json set",
                            "sha256": digest.hexdigest()})
    cases.sort(key=lambda c: (c["layout"], c["mass_over_rm"], c["epsilon"]))

    record = {
        "schema": "NSC-RN-ISOTROPY-v1",
        "api_version": isotropy.API_VERSION,
        "creation_only": True,
        "producing_commit": _producing_commit(),
        "arm": args.arm,
        "arm_meaning": ARMS[args.arm],
        "sources": sources,
        "registration": {
            "isotropic_band": isotropy.ISOTROPIC_BAND,
            "decay_fraction": isotropy.RELAXATION_DECAY_FRACTION,
            "refinement_multiple": isotropy.REFINEMENT_MULTIPLE,
            "written": "before the relaxation run exists",
        },
        "limits": [
            "the index measures the matter balance, not the metric anisotropy",
            "every recorded station is driven: the shell is still falling, so a "
            "trend here is a statement about the driven phase only",
        ],
        "coverage": {
            "layouts": [s["layout"] for s in sources],
            "cases": len(cases),
            "total_stations": sum(c["station_count"] for c in cases),
        },
        "summary": {
            f"{c['layout']}/{c['case']}": {
                "first": c["isotropy_index_first"],
                "max": c["isotropy_index_max"],
                "last": c["isotropy_index_last"],
                "registered_verdict": (c["relaxation_registration"] or {}).get("verdict"),
            } for c in cases
        },
        "cases": cases,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  layouts : {record['coverage']['layouts']}")
    print(f"  cases   : {record['coverage']['cases']}  stations: {record['coverage']['total_stations']}")
    for case in cases:
        reg = case["relaxation_registration"] or {}
        print(
            f"  {case['layout'][:18]:18s} {case['case']:22s} "
            f"I_first={case['isotropy_index_first']:.4f} I_max={case['isotropy_index_max']:.4f} "
            f"I_last={case['isotropy_index_last']:.4f} -> {reg.get('verdict')}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())