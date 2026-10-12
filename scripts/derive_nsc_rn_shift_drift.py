"""Reduce the whole causal story to one radius-varying quantity.

Everything measured in this line about causality reduces to a single
number's drift with radius. The lapse is pinned at one, the two
characteristic speeds sum to twice it, and the only thing in the metric
that changes with radius is the shift:

    speed_plus  =  N - beta
    speed_minus = -(N + beta)
    speed_plus - speed_minus = 2N = 2        exactly, at every radius

So light rides the same amount of space per unit time always, which is
what makes its local speed a constant. What varies is how far the space
has pulled the two directions apart, and that asymmetry is ``beta``.

The outgoing leg goes negative when ``beta`` exceeds the lapse. That is a
coordinate statement about the reference chart. The causal statement is
that the evolved metric has an outermost simple zero of its outgoing
speed, and this checks the two agree, which is a consistency test between
the reference geometry and the marched one rather than a restatement of
either.

Reads the reference chart and, for the comparison, recorded station
payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_shift_drift.py
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
for _path in (str(LAB / "src"), str(LAB / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from recursive_horizons import nsc_rn_observables as obs  # noqa: E402
from recursive_horizons import nsc_rn_pg as pg  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402
import derive_nsc_rn_parent as driver  # noqa: E402

DEVELOPMENT = LAB / "results/development"
STATION = re.compile(r"station-([0-9]+(?:\.[0-9]+)?)\.npz$")

#: The light constant: speed_plus - speed_minus must equal twice the lapse.
LIGHT_CONSTANT_TOLERANCE = 1.0e-9


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _zero_of_outgoing(radius: np.ndarray, outgoing: np.ndarray):
    for index in range(1, radius.size):
        left, right = outgoing[index - 1], outgoing[index]
        if left == 0.0:
            return float(radius[index - 1])
        if left < 0.0 <= right:
            span = right - left
            fraction = 0.0 if span == 0.0 else -left / span
            return float(radius[index - 1]
                         + fraction * (radius[index] - radius[index - 1]))
    return None


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mu", type=float, default=1.04)
    parser.add_argument("--samples", type=int, default=5000)
    parser.add_argument("--run", default="nsc-rn-inner-packet-exc07555-long",
                        help="a recorded run, for the evolved-metric comparison")
    parser.add_argument("--out", default="nsc-rn-shift-drift-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    geometry = driver.case_geometry(args.mu)
    model = pg.action_model(geometry["mass"], geometry["r_m"])
    r_minus = float(model.horizons["r_minus"])
    r_plus = float(model.horizons["r_plus"])
    gap = r_plus - r_minus

    grid = np.linspace(r_minus + 1.0e-3, 32.0, args.samples)
    chart = source.vacuum_chart(model, grid)
    outgoing = np.asarray(chart["speed_plus"], dtype=float)
    ingoing = np.asarray(chart["speed_minus"], dtype=float)
    lapse = np.asarray(chart["lapse"], dtype=float)
    shift = np.asarray(chart["shift"], dtype=float)

    constant = outgoing - ingoing
    reference_zero = _zero_of_outgoing(grid, outgoing)

    # The evolved metric, from the last recorded station.
    evolved_zero = None
    evolved_spacing = None
    directory = DEVELOPMENT / args.run
    stations = sorted((p for p in directory.glob("station-*.npz")
                       if STATION.search(p.name)),
                      key=lambda p: float(STATION.search(p.name).group(1)))
    if stations:
        with np.load(stations[-1]) as payload:
            radius = np.asarray(payload["radius"], dtype=float)
            evolved_zero = obs.trapping_horizon(
                radius,
                np.asarray(payload["frozen_lapse"], dtype=float),
                np.asarray(payload["frozen_shift"], dtype=float))["trapping_horizon"]
            evolved_spacing = float(radius[1] - radius[0])

    record = {
        "schema": "NSC-RN-SHIFT-DRIFT-v1",
        "api_version": "nsc-rn-shift-drift-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "That the causal structure of this metric reduces to one radius-varying "
            "quantity. The lapse is pinned, the two characteristic speeds sum to twice it, "
            "and the shift is the only thing that changes with radius. The horizon is "
            "where the shift has drifted far enough for the outgoing leg to go negative."
        ),
        "light_constant": {
            "identity": "speed_plus - speed_minus = 2N",
            "lapse_min": float(lapse.min()),
            "lapse_max": float(lapse.max()),
            "constant_min": float(constant.min()),
            "constant_max": float(constant.max()),
            "holds": bool(float(np.max(np.abs(constant - 2.0 * lapse)))
                         <= LIGHT_CONSTANT_TOLERANCE),
            "meaning": (
                "Light rides the same amount of space per unit time at every radius, "
                "which is what makes its local speed a constant. Only the split varies."
            ),
        },
        "shift_drift": {
            "beta_at_inner": float(shift[0]),
            "beta_at_outer": float(shift[-1]),
            "beta_monotone_decreasing": bool(np.all(np.diff(shift) < 0.0)),
            # beta = sqrt(2M/r - Q^2/r^2) is extremal where its radicand is,
            # at r = Q^2/M. On the default grid the domain starts above that,
            # so beta falls throughout. Only a grid reaching below it sees
            # the turn, and the turn is a property of the charge rather
            # than of anything dynamical.
            "beta_max_radius": float(grid[int(np.argmax(shift))]),
            "beta_max": float(shift.max()),
            "beta_min_radius": float(grid[int(np.argmin(shift))]),
            "beta_min": float(shift.min()),
            "beta_turn_is_above_default_edge": bool(
                float(grid[int(np.argmax(shift))]) >= 0.5 * (r_minus + r_plus)),
            "beta_where_outgoing_vanishes":
                float(shift[int(np.argmin(np.abs(outgoing)))])
                if reference_zero else None,
            "beta_is_the_only_radius_varying_quantity": (
                "the lapse is identically one and the speed sum is identically two; "
                "the shift carries all of the spatial variation"
            ),
        },
        "chart": {
            "mu": args.mu,
            "r_minus": r_minus,
            "r_plus": r_plus,
            "child_start": r_plus + gap / 16.0,
            "child_end": r_plus + gap / 2.0,
            "samples": int(args.samples),
        },
        "reference_zero_of_outgoing": reference_zero,
        "evolved_trapping_horizon": evolved_zero,
        "evolved_grid_spacing": evolved_spacing,
        "reference_and_evolved_agree": (
            None if (reference_zero is None or evolved_zero is None)
            else bool(abs(reference_zero - evolved_zero) <= evolved_spacing)),
        "agreement_in_cells": (
            None if (reference_zero is None or evolved_zero is None or not evolved_spacing)
            else abs(reference_zero - evolved_zero) / evolved_spacing),
        "scope_note": (
            "A construction check on one reference chart and one recorded run. It does "
            "not measure the actual universe. The identity is exact for this family and "
            "says nothing about any other; the comparison between the reference zero and "
            "the evolved trapping horizon is a consistency check, not a derivation of one "
            "from the other."
        ),
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    constant_block = record["light_constant"]
    print(f"  light constant holds : {constant_block['holds']}"
          f"   (lapse {constant_block['lapse_min']:.9f}..{constant_block['lapse_max']:.9f})")
    drift = record["shift_drift"]
    print(f"  beta                 : {drift['beta_at_inner']:.5f} -> "
          f"{drift['beta_at_outer']:.5f}  monotone={drift['beta_monotone_decreasing']}")
    print(f"  beta turns at r      : {drift['beta_max_radius']:.4f}"
          f" (max {drift['beta_max']:.5f});"
          f" above the default edge: {drift['beta_turn_is_above_default_edge']}")
    print(f"  reference zero       : {reference_zero:.6f}"
          f"   (r_plus = {r_plus:.6f})")
    print(f"  evolved horizon      : {evolved_zero}")
    print(f"  agree within a cell  : {record['reference_and_evolved_agree']}"
          f"   ({record['agreement_in_cells']:.3f} cells)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())