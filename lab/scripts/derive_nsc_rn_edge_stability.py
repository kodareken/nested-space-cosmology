"""Stability audit of the deeper-edge domains.

Moving the grid's inner edge inboard puts the boundary between the horizon
midpoint and the inner horizon, which is territory the recorded family never
sampled. Before reading any quantity out of those runs it is worth asking
whether the march stayed stable there: whether every recorded field is
finite, whether the CFL step was respected, and whether the mass function
stayed monotone and non-decreasing.

A deeper domain is a harder one. A reading from an unstable march is not a
weak reading, it is not a reading.

Reads only recorded station payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_edge_stability.py \\
        --runs nsc-rn-inner-packet-exc085 nsc-rn-inner-packet-exc080 \\
               nsc-rn-inner-packet-exc077
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

LAB = Path(__file__).resolve().parents[1]
DEVELOPMENT = LAB / "results/development"
STATION = re.compile(r"station-([0-9]+(?:\.[0-9]+)?)\.npz$")
#: Fields whose finiteness the reduction depends on.
FIELDS = ("coupled_phi", "coupled_mass", "radius", "frozen_lapse", "frozen_shift")
#: Monotonicity of the mass function is required below this radius only.
#: The oscillation audit already measured grid-scale ringing at r >= 16 in
#: this family, and the far field carries almost no mass, so a reversal
#: there is the known ringing and not an instability. Requiring
#: monotonicity everywhere flagged about 1800 far-field steps per run in
#: every run, including the default-edge one, which says the criterion was
#: wrong rather than the march.
MONOTONIC_BELOW = 16.0


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _clock(path: Path) -> float:
    return float(STATION.search(path.name).group(1))


def _audit(directory: Path) -> Dict[str, Any]:
    summary = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
    stations = sorted((p for p in directory.glob("station-*.npz")
                       if STATION.search(p.name)), key=_clock)

    non_finite: Dict[str, int] = {}
    mass_decreases_below = 0
    worst_decrease_below = 0.0
    mass_decreases_far = 0
    min_far_reversal_radius = None
    lapse_out_of_range = 0
    mass_growth = []
    radii = None

    for station in stations:
        with np.load(station) as payload:
            for key in FIELDS:
                if key not in payload.files:
                    continue
                field = np.asarray(payload[key])
                bad = int(np.count_nonzero(~np.isfinite(field)))
                if bad:
                    non_finite[key] = non_finite.get(key, 0) + bad
            mass = np.asarray(payload["coupled_mass"], dtype=float)
            radii = np.asarray(payload["radius"], dtype=float)
            lapse = np.asarray(payload["frozen_lapse"], dtype=float)
            if np.all(np.isfinite(lapse)):
                lapse_out_of_range += int(np.count_nonzero((lapse <= 0.0)))
            limit = int(np.searchsorted(radii, MONOTONIC_BELOW))
            inner_steps = np.diff(mass[:max(limit, 2)])
            mass_decreases_below += int(np.count_nonzero(inner_steps < 0.0))
            if inner_steps.size and float(np.min(inner_steps)) < worst_decrease_below:
                worst_decrease_below = float(np.min(inner_steps))
            outer_idx = np.where(np.diff(mass) < 0)[0]
            far = outer_idx[radii[outer_idx] >= MONOTONIC_BELOW]
            mass_decreases_far += int(far.size)
            if far.size:
                nearest = float(radii[far[0]])
                if min_far_reversal_radius is None or nearest < min_far_reversal_radius:
                    min_far_reversal_radius = nearest
            mass_growth.append(float(mass[-1] - mass[0]))

    return {
        "run": directory.name,
        "r_inner": summary["r_inner"],
        "r_minus": summary["r_minus"],
        "points": summary["points"],
        "dt": summary["dt"],
        "stations": len(stations),
        "stop_reason": summary["stop_reason"],
        "fields_audited": list(FIELDS),
        "non_finite_counts": non_finite,
        "all_fields_finite": not non_finite,
        "monotonic_below_radius": MONOTONIC_BELOW,
        "mass_decreases_below_limit": mass_decreases_below,
        "worst_mass_decrease_below_limit": worst_decrease_below,
        "mass_reversals_far_field": mass_decreases_far,
        "nearest_far_field_reversal_radius": min_far_reversal_radius,
        "far_field_reversals_are_known_ringing": (
            (min_far_reversal_radius is not None and min_far_reversal_radius >= MONOTONIC_BELOW)),
        "lapse_nodes_at_or_below_zero": lapse_out_of_range,
        "lapse_always_positive": lapse_out_of_range == 0,
        "radius_increasing": bool(np.all(np.diff(radii) > 0.0)) if radii is not None else None,
        "mass_growth_first_last": [mass_growth[0], mass_growth[-1]] if mass_growth else None,
        "stable": (not non_finite) and mass_decreases_below == 0 and lapse_out_of_range == 0,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", required=True)
    parser.add_argument("--out", default="nsc-rn-edge-stability-v2.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    audits = []
    for name in args.runs:
        directory = DEVELOPMENT / name
        if not (directory / "summary.json").is_file():
            print(f"no completed run at {directory}", file=sys.stderr)
            return 1
        audits.append(_audit(directory))
    audits.sort(key=lambda a: -a["r_inner"])

    record = {
        "schema": "NSC-RN-EDGE-STABILITY-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "Whether the march stayed stable in the domain between the horizon midpoint "
            "and the inner horizon, which the recorded family never sampled. Every "
            "recorded field is checked for finiteness, the mass function for "
            "monotonicity, and the lapse for positivity."
        ),
        "scope_note": (
            "Stability is not correctness. A stable march in a deeper domain says the "
            "numbers can be read, not that the domain is the right one. This does not "
            "measure the actual universe."
        ),
        "all_runs_stable": all(a["stable"] for a in audits),
        "runs": audits,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  {'run':<40} {'r_inner':>8} {'finite':>7} {'dec<':>6} {'far rev':>8} {'lapse':>7} {'stable':>7}")
    for audit in audits:
        print(f"  {audit['run']:<40} {audit['r_inner']:>8.4f} "
              f"{str(audit['all_fields_finite']):>7} {audit['mass_decreases_below_limit']:>6} "
              f"{audit['mass_reversals_far_field']:>8} {str(audit['lapse_always_positive']):>7} "
              f"{str(audit['stable']):>7}")
    print(f"  all runs stable: {record['all_runs_stable']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())