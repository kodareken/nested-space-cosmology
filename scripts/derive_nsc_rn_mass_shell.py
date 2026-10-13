"""Reduce the recorded RN radial family to the motion of its mass shell.

Reads every station written by ``nsc-rn-neutral-family-v1``, locates the
mass shell in each, and fits its initial radius and inward speed per case.
No time step is taken, nothing is re-run, and no existing file is
rewritten.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_mass_shell.py
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

from recursive_horizons import nsc_rn_mass_shell as shell  # noqa: E402

FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"
RECORD = LAB / "results/development/nsc-rn-mass-shell-v6.json"
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
    mass_over_rm = float(match.group("mu")) if match else None
    epsilon = float(match.group("eps")) if match else None
    stations = sorted((p for p in case_dir.glob("station-*.npz")
                       if STATION_PATTERN.match(p.name)), key=_station_time)

    rows: List[Dict[str, Any]] = []
    points = None
    for station in stations:
        with np.load(station) as payload:
            if points is None:
                points = int(payload["radius"].shape[0])
            state = shell.shell_state(payload, base_mass=mass_over_rm)
            clock = (float(payload["coupled_t"][0])
                     if "coupled_t" in payload.files else _station_time(station))
        state["t_over_rm"] = clock
        rows.append(state)

    motion = shell.shell_motion(rows)
    last = rows[-1] if rows else {}
    return {
        "case": case_dir.name,
        "mass_over_rm": mass_over_rm,
        "epsilon": epsilon,
        "points": points,
        "station_count": len(rows),
        "reached_t_over_rm": last.get("t_over_rm"),
        "motion": motion,
        "final_shell_radius": last.get("shell_radius"),
        "final_fraction_below_inner_edge": (
            last["gain_below_inner_edge"] / last["total_gain"]
            if last.get("total_gain") else None),
        "stations": rows,
    }


def main() -> int:
    if RECORD.exists():
        print(f"refusing to overwrite sealed record: {RECORD}", file=sys.stderr)
        return 1
    if not FAMILY.is_dir():
        print(f"missing family directory: {FAMILY}; restore it from Zenodo DOI 10.5281/zenodo.23090015", file=sys.stderr)
        return 1

    cases = [_reduce_case(p) for p in sorted(FAMILY.iterdir())
             if p.is_dir() and CASE_PATTERN.match(p.name)]

    starts = [c["motion"]["intercept_radius"] for c in cases if c["motion"].get("fit_available")]
    rates: Dict[float, list] = {}
    all_rates: list = []
    for case in cases:
        if not case["motion"].get("fit_available"):
            continue
        rate = case["motion"]["infall_rate_mean"]
        rates.setdefault(case["mass_over_rm"], []).append(rate)
        all_rates.append(rate)
    amps: Dict[float, list] = {}
    for case in cases:
        if case["stations"]:
            amps.setdefault(case["epsilon"], []).append(case["stations"][0]["shell_amplitude"])

    # The question the family is asked: do the two knobs set two different
    # things, or do they both move the mass geometry together?
    spread_over_mu = max(all_rates) - min(all_rates) if all_rates else None
    means_over_mu = {str(k): sum(v) / len(v) for k, v in sorted(rates.items())}
    mu_means = list(means_over_mu.values())
    mu_spread = (max(mu_means) - min(mu_means)) if mu_means else None

    eps_keys = sorted(amps)
    eps_ratio = (sum(amps[eps_keys[-1]]) / len(amps[eps_keys[-1]])
                 / (sum(amps[eps_keys[0]]) / len(amps[eps_keys[0]]))) if len(eps_keys) >= 2 else None

    record = {
        "schema": "NSC-RN-MASS-SHELL-v1",
        "api_version": shell.API_VERSION,
        "creation_only": True,
        "producing_commit": _producing_commit(),
        "family_source": FAMILY.name,
        "family_source_sha256": _sha256_file(FAMILY / "prepare.json"),
        "what_this_measures": (
            "The geometry itself. In the PG chart the geometry is fixed by M(r), so the "
            "shape of dM/dr is the shape of the object. The force diagnostics in this "
            "line measure what the two pressures do; this measures where the mass is."
        ),
        "domain": {
            "chart": "areal radius",
            "r_outer": 32.0,
            "observer": "areal chart of the recorded coupled arm",
            "state": "recorded station payload, no time step taken",
            "compared_quantity": "radius of the smoothed dM/dr peak and its inward drift",
        },
        "scope_note": (
            "A construction check on recorded finite states. It does not measure the "
            "actual universe. The shell leaves the computed domain at the inner edge in "
            "every case, so what happens to it afterwards is not measured here."
        ),
        "method": {
            "shell_density": "centred dM/dr of the recorded coupled_mass",
            "smooth_window": shell.SMOOTH_WINDOW,
            "search_limit": shell.SEARCH_LIMIT,
            "feature_continuity": shell.FEATURE_CONTINUITY,
            "why_smoothed": (
                "The oscillation audit measured ringing of one to three cells over much "
                "of this profile. A median window of fifteen cells sits far above that and "
                "well below the shell's own width."
            ),
        },
        "coverage": {
            "cases": len(cases),
            "total_stations": sum(c["station_count"] for c in cases),
            "cases_with_a_fit": sum(1 for c in cases if c["motion"].get("fit_available")),
        },
        "separation": {
            "reading": (
                "Two knobs, and they do not act on the same quantity. The mass ratio sets "
                "the force balance found earlier in this line; in the mass geometry it "
                "changes almost nothing. The source strength scales the shell amplitude "
                "exactly, decade for decade. So the family separates a force question from "
                "a mass question rather than sweeping one dial."
            ),
            "infall_rate_all": sorted(all_rates),
            "infall_rate_spread_all": spread_over_mu,
            "infall_rate_mean_by_mass_ratio": means_over_mu,
            "infall_rate_spread_of_means_over_mu": mu_spread,
            "infall_rate_spread_relative": (mu_spread / abs(sum(mu_means) / len(mu_means))
                                            if mu_means and sum(mu_means) else None),
            "source_scaling_ratio_between_largest_and_smallest_epsilon": eps_ratio,
            "shell_start_radius_min": min(starts) if starts else None,
            "shell_start_radius_max": max(starts) if starts else None,
            "what_this_does_not_show": (
                "It does not show that the mass ratio is irrelevant. The mass ratio moves "
                "the radial and angular pressures, and therefore the balance wall and the "
                "force reach, both of which were measured and both of which do respond. "
                "What is insensitve to it here is the mass geometry alone."
            ),
        },
        "cases": cases,
    }

    RECORD.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {RECORD.relative_to(LAB)}")
    print(f"  cases            : {record['coverage']['cases']}")
    print(f"  stations         : {record['coverage']['total_stations']}")
    print(f"  falling phase    : {sum(c['motion']['stations_falling'] for c in cases)} of"
          f" {record['coverage']['total_stations']} stations")
    print(f"  edge arrival t   : "
          f"{sorted({c['motion']['inner_edge_arrival_t'] for c in cases})}")
    print(f"  shell start r    : {record['separation']['shell_start_radius_min']:.3f}"
          f" .. {record['separation']['shell_start_radius_max']:.3f}")
    print(f"  infall rate      : {record['separation']['infall_rate_all'][0]:.4f}"
          f" .. {record['separation']['infall_rate_all'][-1]:.4f}"
          f"   (spread {record['separation']['infall_rate_spread_all']:.4f})")
    for key, value in means_over_mu.items():
        print(f"    mean rate mu={key}: {value:.4f}")
    print(f"  eps scaling ratio: {eps_ratio:.3f}  (expected 10 for a linear source)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
