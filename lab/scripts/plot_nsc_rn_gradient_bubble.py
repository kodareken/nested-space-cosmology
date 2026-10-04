"""Render the two-force geometry and the regularity verdict.

Three panels, all from recorded states and the v4 reduction:

1. the two force profiles on the areal radius, with the balance wall
   marked, which is the shape of the nested object;
2. the wall radius against time for every case, which is the oscillation
   and the trend with exterior strength;
3. the throat slope against the sourced regularity limit, which is the
   negative result.

Run from the repository root:

    python scripts/lab.py scripts/plot_nsc_rn_gradient_bubble.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_rn_gradient_bubble as bubble  # noqa: E402

FAMILY = LAB / "results/development/nsc-rn-neutral-family-v1"
RECORD = LAB / "results/development/nsc-rn-gradient-bubble-v4.json"
OUT = LAB / "results/development/nsc-rn-gradient-bubble-v4.png"
CASE_PATTERN = re.compile(r"^mu-(?P<mu>[0-9.]+)-eps-(?P<eps>[0-9.]+)$")
STATION_PATTERN = re.compile(r"^station-(?P<t>[0-9.]+)\.npz$")


def _station_time(path: Path) -> float:
    match = STATION_PATTERN.match(path.name)
    return float(match.group("t")) if match else float("inf")


def _load_record() -> Dict:
    return json.loads(RECORD.read_text(encoding="utf-8"))


def _last_station(case_dir: Path) -> Path:
    stations = sorted((p for p in case_dir.glob("station-*.npz")
                       if STATION_PATTERN.match(p.name)), key=_station_time)
    return stations[-1]


def panel_forces(ax, case_dirs: Dict[str, Path]) -> None:
    """The two force profiles on the areal radius, at the last station."""
    for case, directory in sorted(case_dirs.items()):
        with np.load(_last_station(directory)) as payload:
            radius = payload["radius"]
            inner = np.abs(payload["radial_pressure"])
            outer = np.abs(payload["angular_pressure"])
            profile = bubble.gradient_bubble_profile(payload)
        strong = CASE_PATTERN.match(case).group("mu")
        ax.loglog(radius, inner, color="#1b3a6b", lw=1.6,
                  label=f"radial (inherited)  $\\mu$={strong}")
        ax.loglog(radius, outer, color="#b3452c", lw=1.4, ls="--",
                  label=f"angular (generated) $\\mu$={strong}")
        walls = profile["walls"]
        if walls.size:
            for wall in walls:
                ax.axvline(wall, color="#444444", lw=0.6, alpha=0.35)
    ax.set_xlabel("areal radius $r$   [$r_m$]")
    ax.set_ylabel("force density  [log]")
    ax.set_title("The two forces, and the wall where they are equal\n"
                 "(vertical lines: balance walls; dark = radial, red = angular)",
                 fontsize=9)
    ax.legend(fontsize=6, ncol=2, frameon=False)
    ax.grid(alpha=0.25, which="both")


def panel_wall(ax, record: Dict) -> None:
    """Wall radius against time for every case."""
    for case in record["cases"]:
        times = []
        radii = []
        for row in case["stations"]:
            match = STATION_PATTERN.match(row["station"] + ".npz")
            times.append(float(match.group("t")) if match else 0.0)
            radii.append(row["r_first_wall"] if row["wall_count"] else np.nan)
        mu = case["mass_over_rm"]
        ax.plot(times, radii, marker="o", ms=3, lw=1.2,
                label=f"$\\mu$={mu}, $\\varepsilon$={case['epsilon']}")
    ax.set_xlabel("time  $t/r_m$")
    ax.set_ylabel("balance wall radius $r_{\\mathrm{wall}}$")
    ax.set_title("The wall oscillates, and stronger compression dissolves it\n"
                 "(gaps are stations with no wall)", fontsize=9)
    ax.legend(fontsize=6, frameon=False, ncol=2)
    ax.grid(alpha=0.25)


def panel_regularity(ax, record: Dict) -> None:
    """Throat slope against the sourced regularity limit."""
    names = []
    slopes = []
    validity = {row["case"]: row for row in record["numerical_validity"]["throat_slopes"]}
    for case in record["cases"]:
        names.append(f"$\\mu$={case['mass_over_rm']}\n$\\varepsilon$={case['epsilon']}")
        slopes.append(abs(validity[case["case"]]["slope_exact"]))
    limit = abs(record["limits"]["regular_slope_limit"])
    positions = np.arange(len(names))
    ax.bar(positions, slopes, color="#b3452c", alpha=0.85)
    ax.axhline(limit, color="#1b3a6b", lw=2.0)
    ax.annotate("Hayward regularity limit  $g=-1$\n(gr-qc/0506126, 2005)",
                xy=(len(names) - 0.5, limit), xytext=(-6, 10),
                textcoords="offset points", ha="right", fontsize=7, color="#1b3a6b")
    ax.set_yscale("log")
    ax.set_xticks(positions)
    ax.set_xticklabels(names, fontsize=6)
    ax.set_ylabel("$|$ throat log-slope $|$")
    ax.set_title("Every case is far steeper than the sourced limit:\n"
                 "no case has a regular interior", fontsize=9)
    ax.grid(alpha=0.25, axis="y")


def main() -> int:
    if not RECORD.is_file():
        print(f"missing record: {RECORD}", file=sys.stderr)
        return 1
    record = _load_record()
    case_dirs = {
        p.name: p for p in sorted(FAMILY.iterdir())
        if p.is_dir() and CASE_PATTERN.match(p.name)
    }

    figure, axes = plt.subplots(1, 3, figsize=(16.5, 4.9))
    panel_forces(axes[0], case_dirs)
    panel_wall(axes[1], record)
    panel_regularity(axes[2], record)
    figure.suptitle(
        "RN radial parent: the nested geometry, the oscillating wall, and the regularity failure",
        fontsize=11)
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(OUT, dpi=140)
    print(f"wrote {OUT.relative_to(LAB)}  ({OUT.stat().st_size} bytes)")
    print(f"  cases plotted   : {len(case_dirs)}")
    print(f"  station fraction with a wall: "
          f"{record['coverage']['station_fraction_with_a_wall']:.3f}")
    print(f"  cases with a regular interior: "
          f"{record['summary']['cases_with_a_regular_interior']} of {record['coverage']['cases']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())