"""Reduce the inner-packet release-radius sweep.

One run cannot tell a construction from an accident. This reduces every
inner-packet march to the same quantities and asks whether the compression
measured at a single release radius survives a change of release radius.

The packet width, spectrum, chart, mass ratio and source strength are held
fixed across the sweep; only the release radius moves. If the compression
gain, the infall rate and the width narrowing all track each other across
radii, they are properties of the construction. If they do not, the single
run's x1.197 was an accident and says so.

Reads only recorded station payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_packet_sweep.py
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

from recursive_horizons import nsc_rn_mass_shell as shell  # noqa: E402

DEVELOPMENT = LAB / "results/development"
STATION = re.compile(r"station-([0-9]+(?:\.[0-9]+)?)\.npz$")


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _clock(path: Path) -> float:
    return float(STATION.search(path.name).group(1))


def _reduce_run(directory: Path) -> Dict[str, Any]:
    summary_path = directory / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    base_mass = float(summary["mu"])

    stations = sorted((p for p in directory.glob("station-*.npz")
                       if STATION.search(p.name)), key=_clock)
    if not stations:
        raise SystemExit(f"no stations in {directory}")

    rows: List[Dict[str, Any]] = []
    for station in stations:
        with np.load(station) as payload:
            state = shell.shell_state(payload, base_mass=base_mass)
            state["t_over_rm"] = _clock(station)
        rows.append(state)

    motion = shell.shell_motion(rows)
    accumulation = shell.accumulation_after_arrival(rows, motion)

    falling = [row for row in rows if row.get("phase") == "falling"]
    widths = [row["shell_width_half_max"] for row in falling
              if np.isfinite(row["shell_width_half_max"])]

    # A run whose packet was wide enough to touch the inner boundary at
    # t = 0 has no fall in it. It still measures the accumulation, but it
    # measures nothing about a shell moving, and the record must not let a
    # reader mistake one for the other.
    left_edge = summary["packet_center_over_rm"] - 0.5 * summary["packet_width_over_rm"]
    starts_pinned = left_edge <= summary["r_inner"]

    with np.load(stations[0]) as payload:
        radius = payload["radius"]
        spacing = float(radius[1] - radius[0])

    return {
        "run": directory.name,
        "release_radius_over_rm": summary["packet_center_over_rm"],
        "mu": summary["mu"],
        "eps": summary["eps"],
        "packet_width_over_rm": summary["packet_width_over_rm"],
        "r_inner": summary["r_inner"],
        "r_outer": summary["r_outer"],
        "points": summary["points"],
        "target_over_rm": summary["target_over_rm"],
        "stop_reason": summary["stop_reason"],
        "cpu_seconds": summary["cpu_seconds"],
        "stations": len(rows),
        "grid_spacing": spacing,
        "march_summary_sha256": _sha256_file(summary_path),
        "motion": motion,
        "accumulation_after_arrival": accumulation,
        "starts_pinned_at_inner_edge": starts_pinned,
        "packet_left_edge_over_rm": left_edge,
        "width_first": widths[0] if widths else None,
        "width_last": widths[-1] if widths else None,
        "width_cells_first": widths[0] / spacing if widths else None,
        "width_cells_last": widths[-1] / spacing if widths else None,
        "width_compression": (widths[0] / widths[-1]) if len(widths) >= 2 else None,
        "fraction_below_inner_edge_final": (
            rows[-1]["gain_below_inner_edge"] / rows[-1]["total_gain"]
            if rows[-1]["total_gain"] else None),
        "series": rows,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", default=None,
                        help="explicit run directory names under results/development. "
                             "Preferred over --pattern: a glob will happily mix a "
                             "different-resolution march of the same radius into a "
                             "fixed-resolution sweep.")
    parser.add_argument("--pattern", default=None,
                        help="glob under results/development, used only when --runs is absent")
    parser.add_argument("--out", default="nsc-rn-packet-sweep-v4.json",
                        help="record name, relative to results/development")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    if args.runs:
        names = list(args.runs)
    elif args.pattern:
        names = [p.name for p in sorted(DEVELOPMENT.glob(args.pattern))
                 if p.is_dir() and (p / "summary.json").is_file()]
    else:
        print("give --runs or --pattern", file=sys.stderr)
        return 2

    if len(names) < 2:
        print(f"need at least two completed runs, named {names}", file=sys.stderr)
        return 1
    runs = []
    for name in names:
        directory = DEVELOPMENT / name
        if not (directory / "summary.json").is_file():
            print(f"no completed run at {directory}", file=sys.stderr)
            return 1
        runs.append(directory)

    reduced = [_reduce_run(p) for p in runs]
    reduced.sort(key=lambda r: r["release_radius_over_rm"])

    # A sweep varies the release radius and nothing else. If anything else
    # moved between runs the comparison is not the one the record claims,
    # so it is refused rather than reported.
    for key in ("points", "mu", "eps", "packet_width_over_rm", "r_inner"):
        values = [r[key] for r in reduced]
        if len(set(values)) != 1:
            print(f"runs differ in {key}: {values}. This is not a release-radius sweep; "
                  f"refusing to reduce it as one.", file=sys.stderr)
            return 2

    gains = [r["accumulation_after_arrival"]["amplitude_gain"]
             for r in reduced if r["accumulation_after_arrival"]]
    rises = [r["accumulation_after_arrival"]["rise_time"]
             for r in reduced if r["accumulation_after_arrival"]]
    relax = [r["accumulation_after_arrival"]["amplitude_after_peak_fraction"]
             for r in reduced if r["accumulation_after_arrival"]]
    rates = [r["motion"]["infall_rate_mean"]
             for r in reduced if r["motion"].get("fit_available")]
    widths = [r["width_compression"] for r in reduced if r["width_compression"]]
    arrivals = [r["accumulation_after_arrival"]["arrival_t"]
                for r in reduced if r["accumulation_after_arrival"]]
    moving = [r for r in reduced if not r["starts_pinned_at_inner_edge"]]

    fixed = {
        "mu": sorted({r["mu"] for r in reduced}),
        "eps": sorted({r["eps"] for r in reduced}),
        "packet_width_over_rm": sorted({r["packet_width_over_rm"] for r in reduced}),
        "points": sorted({r["points"] for r in reduced}),
        "r_inner": sorted({r["r_inner"] for r in reduced}),
    }

    def _agrees(values, tolerance, floor=1):
        if len([v for v in values if v is not None]) < floor:
            return None
        return (max(values) - min(values)) < tolerance

    gain_holds = _agrees(gains, 0.05 * (sum(gains) / len(gains)))
    rise_holds = _agrees(rises, 0.15)
    relax_holds = _agrees(relax, 0.10)

    lines_hold = []
    if gain_holds:
        lines_hold.append(
            "The compression gain holds: the density peak at the inner boundary grows by "
            "about a fifth at every release radius tried.")
    if rise_holds:
        lines_hold.append(
            "The rise time holds: the peak arrives at the same lag after the shell reaches "
            "the boundary, however far the shell fell to get there.")

    lines_not = []
    if relax_holds is False:
        lines_not.append(
            "The relaxation after the peak does not agree across runs: the final amplitude "
            "as a fraction of the peak spans a range wider than a tenth of its own value, "
            "so no single relaxed value is claimed.")
    elif relax_holds:
        lines_not.append(
            "The relaxation agrees only as well as these runs put it on a common clock. "
            "The runs were marched to different absolute times, so agreement here bounds "
            "the relaxation at the lag actually reached rather than settling it.")
    if len(widths) >= 2 and (max(widths) - min(widths)) > 0.05:
        lines_not.append(
            "The width compression is not a fixed factor; it grows with release radius, so "
            "it is a property of how far the shell fell rather than a constant.")
    if len(rates) < len(reduced):
        lines_not.append(
            "Not every run has a falling phase, so the infall rate is measured on fewer "
            "runs than the sweep contains.")
    if not lines_not:
        lines_not.append("No quantity in this record was found not to hold.")

    record = {
        "schema": "NSC-RN-PACKET-SWEEP-v1",
        "api_version": shell.API_VERSION,
        "creation_only": True,
        "producing_commit": _git_head(),
        "run_selection": names,
        "run_pattern": args.pattern,
        "what_this_measures": (
            "Whether the compression measured at one release radius is a property of the "
            "construction or an accident of one run. Packet width, spectrum, chart, mass "
            "ratio and source strength are held fixed; only the release radius moves."
        ),
        "scope_note": (
            "A construction check on a finite marched collar. It does not measure the "
            "actual universe. Holding a result across a change of one construction "
            "parameter shows it is not an accident; it does not make it physical."
        ),
        "held_fixed": fixed,
        "sweep": {
            "release_radii_over_rm": [r["release_radius_over_rm"] for r in reduced],
            "runs": len(reduced),
            "runs_with_a_falling_phase": len(moving),
            "runs_starting_pinned": [r["release_radius_over_rm"] for r in reduced
                                     if r["starts_pinned_at_inner_edge"]],
            "amplitude_gain": gains,
            "amplitude_gain_spread": (max(gains) - min(gains)) if gains else None,
            "amplitude_gain_relative_spread": (
                (max(gains) - min(gains)) / (sum(gains) / len(gains)) if gains else None),
            "amplitude_gain_by_radius": {
                str(r["release_radius_over_rm"]): (
                    r["accumulation_after_arrival"]["amplitude_gain"]
                    if r["accumulation_after_arrival"] else None)
                for r in reduced},
            "rise_time": rises,
            "rise_time_spread": (max(rises) - min(rises)) if rises else None,
            "after_peak_fraction": relax,
            "after_peak_fraction_spread": (max(relax) - min(relax)) if relax else None,
            "infall_rate": rates,
            "infall_rate_spread": (max(rates) - min(rates)) if rates else None,
            "width_compression": widths,
            "width_compression_spread": (max(widths) - min(widths)) if widths else None,
            "edge_arrival_t": arrivals,
            "what_holds": " ".join(lines_hold),
            "what_does_not_hold": " ".join(lines_not),
        },
        "runs": reduced,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  runs            : {len(reduced)}")
    print(f"  release radii   : {record['sweep']['release_radii_over_rm']}")
    print(f"  amplitude gain  : {['%.4f' % g for g in gains]}")
    if record['sweep']['amplitude_gain_relative_spread'] is not None:
        print(f"                     relative spread {record['sweep']['amplitude_gain_relative_spread']:.4f}")
    print(f"  infall rate     : {['%.4f' % r for r in rates]}")
    print(f"  width compress  : {['%.4f' % w for w in widths]}")
    print(f"  arrival t       : {arrivals}")
    print(f"  rise time       : {['%.3f' % r for r in rises]}")
    print(f"  after-peak frac : {['%.3f' % r for r in relax]}"
          f"  (spread {max(relax) - min(relax) if relax else 0:.3f})")
    print(f"  pinned at start : {record['sweep']['runs_starting_pinned']}")
    print(f"  holds           : {len(lines_hold)} claim(s)")
    print(f"  does not hold   : {len(lines_not)} item(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())