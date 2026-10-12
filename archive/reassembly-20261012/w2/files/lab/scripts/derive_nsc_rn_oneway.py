"""Drive the one-way enclosure measure, bounded, across the recorded family.

Builds the instantaneous Dirac rate operator of the final recorded
station of each case and reports the direction split of its eigenmodes
across the state's own trapped boundary. No time step is taken and no
existing file is rewritten.

Cost bound. The dense spectrum of a 1738-order advective operator is
expensive and, for some recorded spectra, the LAPACK eigensolve does not
converge promptly. Each station is therefore evaluated in its own worker
subprocess under a hard wall-clock guard; a station that does not return
in time is recorded with the blocker `eigen_timeout` and no partial
spectrum is claimed. The sample is the final station of each case, where
the shell has sunk farthest inward, which is the decisive instant for
whether the inside can emit outward.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_oneway.py

Any recorded family layout works; the inward-march campaign is reduced with

    python scripts/lab.py scripts/derive_nsc_rn_oneway.py \\
        --family results/development/nsc-rn-inward-march-v1 \\
        --record results/development/nsc-rn-oneway-v3.json
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

from recursive_horizons import nsc_rn_oneway as oneway  # noqa: E402

FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"
RECORD_V1 = LAB / "results/development/nsc-rn-oneway-v1.json"
RECORD = LAB / "results/development/nsc-rn-oneway-v2.json"
CASE_PATTERN = re.compile(r"^mu-(?P<mu>[0-9.]+)-eps-(?P<eps>[0-9.]+)$")
STATION_PATTERN = re.compile(r"^station-(?P<t>[0-9.]+)\.npz$")
WORKER_TIMEOUT_SECONDS = 240.0

WORKER = r"""
import json, sys, numpy as np
sys.path.insert(0, r"{src}")
from recursive_horizons import nsc_rn_oneway as oneway
with np.load(r"{station}") as payload:
    out = oneway.station_oneway(payload, mass_over_rm={mu})
band = [{{"band_low": b["band"][0], "band_high": b["band"][1],
         "interior_max_direction": b["interior_max_direction"],
         "exterior_min_direction": b["exterior_min_direction"]}} for b in out["bands"]]
row = {{
    "station": "{stem}",
    "t_over_rm": out["t_over_rm"],
    "verdict": out["verdict"],
    "flags": out["flags"],
    "trap_radius": out["trap_radius"],
    "mode_count": out["mode_count"],
    "interior_mode_count": out["interior_mode_count"],
    "interior_max_direction": out["interior_max_direction"],
    "exterior_min_direction": out["exterior_min_direction"],
    "spectral_radius": out["eigenvalue_spectral_radius"],
    "max_real_eigenvalue": out["max_real_eigenvalue"],
    "bands": band,
}}
print(json.dumps(row, sort_keys=True))
"""


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


def _last_station(case_dir: Path) -> Path | None:
    stations = sorted(
        (p for p in case_dir.glob("station-*.npz") if STATION_PATTERN.match(p.name)),
        key=_station_time,
    )
    return stations[-1] if stations else None


def _family_digest(family: Path) -> tuple[str, str]:
    """Digest of the recorded inputs a reduction stands on.

    The neutral family keeps one root ``prepare.json``; a campaign layout keeps
    one per case. The digest folds whichever layout is present and rewrites
    nothing.
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


def _evaluate(station: Path, mu: float) -> Dict[str, Any]:
    code = WORKER.format(
        src=str(LAB / "src"), station=str(station.resolve()), mu=mu, stem=station.stem,
    )
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=LAB, capture_output=True, text=True,
        timeout=WORKER_TIMEOUT_SECONDS,
    )
    if result.returncode != 0:
        return {
            "station": station.stem,
            "t_over_rm": _station_time(station),
            "verdict": "evaluator_error",
            "flags": ["evaluator_error"],
            "blocker": (result.stderr or "worker failed").strip()[-600:],
        }
    return json.loads(result.stdout.strip().splitlines()[-1])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family", default=str(FAMILY.relative_to(LAB)))
    parser.add_argument("--record", default=str(RECORD.relative_to(LAB)))
    parser.add_argument("--supersedes", default=str(RECORD_V1.relative_to(LAB)),
                        help="sealed record this reduction follows")
    parser.add_argument("--supersedes-reason", default=(
        "v1 measured channel weight with speed magnitudes, which counts the "
        "chi_plus content of a trapped mode as outward weight even though both "
        "characteristics are negative inside the trap; the eigen data are "
        "recomputed with the signed current."
    ))
    args = parser.parse_args()
    family = (LAB / args.family).resolve()
    record_path = (LAB / args.record).resolve()
    prior = (LAB / args.supersedes).resolve()

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
        station = _last_station(case_dir)
        if station is None:
            continue
        with np.load(station) as payload:
            witness = oneway.state_witness(dict(payload))
        try:
            row = _evaluate(station, mu)
        except subprocess.TimeoutExpired:
            row = {
                "station": station.stem,
                "t_over_rm": _station_time(station),
                "verdict": "eigen_timeout",
                "flags": ["eigen_timeout"],
                "blocker": f"dense spectrum exceeded {WORKER_TIMEOUT_SECONDS:g}s",
            }
        cases.append({
            "case": case_dir.name,
            "mass_over_rm": mu,
            "epsilon": float(match.group("eps")),
            "state_witness": witness,
            "final_station": row,
        })
    cases.sort(key=lambda c: (c["mass_over_rm"], c["epsilon"]))
    source_kind, source_digest = _family_digest(family)

    supersedes = {
        "record": prior.name,
        "sha256": _sha256_file(prior) if prior.is_file() else None,
        "superseded_reason": args.supersedes_reason,
    }

    flags = sorted({f for c in cases for f in c["final_station"].get("flags", [])})
    cardinal = sum(1 for c in cases
                   if "interior_cardinal_outward" in c["final_station"].get("flags", []))
    timed_out = sum(1 for c in cases
                    if "eigen_timeout" in c["final_station"].get("flags", []))

    record = {
        "schema": "NSC-RN-ONEWAY-v1",
        "api_version": oneway.API_VERSION,
        "creation_only": True,
        "producing_commit": _producing_commit(),
        "family_source": family.name,
        "family_source_kind": source_kind,
        "family_source_sha256": source_digest,
        "sample": "final recorded station of each case",
        "worker_timeout_seconds": WORKER_TIMEOUT_SECONDS,
        "domain": {
            "chart": "areal radius",
            "operator": (
                "instantaneous Dirac rate matrix of nsc_rn_pg.dirac_rates on the "
                "reconstructed coupled lapse and shift; identical to the march route"
            ),
            "state": "recorded station payload, no time step taken",
            "compared_quantity": (
                "transport-direction statistic S of every eigenmode, split at the "
                "state's own measured trap radius (nsc-rn-cone-map)"
            ),
            "one_way_definition": (
                "no interior-supported mode carries outward probability weight "
                "above the direction floor, while exterior modes fall inward"
            ),
        },
        "limits": {"direction_floor": oneway.DIRECTION_FLOOR,
                   "interior_fraction": oneway.INTERIOR_FRACTION},
        "coverage": {"cases": len(cases), "final_station_sample": True},
        "supersedes": supersedes,
        "flags_observed": flags,
        "summary": {
            "stations_with_interior_cardinal_outward": cardinal,
            "stations_with_eigen_timeout": timed_out,
            "cases_strictly_oneway": sum(
                1 for c in cases
                if "interior_cardinal_outward" not in c["final_station"].get("flags", [])
                and "eigen_timeout" not in c["final_station"].get("flags", [])),
        },
        "reading": (
            "If no interior mode carries outward weight, the recorded operator's "
            "instantaneous passage across the measured trap is one way in the "
            "direction set by the 2005 reading. This is a spectral statement about "
            "the recorded state, not the actual universe."
        ),
        "cases": cases,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  cases        : {record['coverage']['cases']}")
    print(f"  cardinal swap: {cardinal}")
    print(f"  eigen timeout: {timed_out}")
    print(f"  flags        : {flags}")
    for case in cases:
        row = case["final_station"]
        print(
            f"  {case['case']:22s} {row.get('station'):14s} "
            f"maxS_in={row.get('interior_max_direction')} minS_out={row.get('exterior_min_direction')}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())