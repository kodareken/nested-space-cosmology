"""Map the region below the inner horizon, which the chart refuses to sample.

The map exposed that the geometry divides into three regions rather than
two: a parent, a trapped region between the two horizons, and a region
below the inner horizon. The third is where a further nest would have to
be, and `vacuum_chart` refuses to sample it with the message "static modes
stay strictly above the inner horizon".

That guard is worth examining rather than accepting. The shift is
``sqrt(2M/r - Q^2/r^2)``, which is real wherever ``r > Q^2/(2M)``. For
``mu = 1.04`` that is ``0.4808``, well below ``r_minus = 0.7543``. So the
real domain of the geometry reaches further in than the chart's guard
allows, and the horizons are where the shift passes one, not where it
stops being real.

This calls the model's own functions directly, without the chart wrapper,
to establish whether the third region is actually describable or whether
something breaks there. It reads the reference model only: no station
payloads, no time step, no marching, and no production guard is changed.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_below_horizon.py
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

import numpy as np

LAB = Path(__file__).resolve().parents[1]
for _path in (str(LAB / "src"), str(LAB / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from recursive_horizons import nsc_rn_pg as pg  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402
import derive_nsc_rn_parent as driver  # noqa: E402

DEVELOPMENT = LAB / "results/development"


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def _entry(mu: float, samples: int) -> Dict[str, Any]:
    geometry = driver.case_geometry(mu)
    model = pg.action_model(geometry["mass"], geometry["r_m"])
    r_minus = float(model.horizons["r_minus"])
    r_plus = float(model.horizons["r_plus"])
    mass = float(model.mass)
    charge2 = float(model.magnetic_r2)
    real_limit = charge2 / (2.0 * mass)

    # Sample the whole real domain of the shift, including the part the
    # chart's guard excludes.
    grid = np.linspace(real_limit * (1.0 + 1e-6), 8.0, samples)
    radicand = 2.0 * mass / grid - charge2 / grid ** 2
    beta = np.sqrt(np.clip(radicand, 0.0, None))
    lapse = np.sqrt(np.clip(1.0 - 2.0 * mass / grid + charge2 / grid ** 2 + radicand,
                            0.0, None))
    outgoing = lapse - beta
    ingoing = -(lapse + beta)

    # Where does the outgoing speed vanish, going inward as well as out?
    crossings: List[float] = []
    for i in range(1, grid.size):
        if outgoing[i - 1] * outgoing[i] < 0.0:
            span = outgoing[i] - outgoing[i - 1]
            crossings.append(float(grid[i - 1]
                                   + (-outgoing[i - 1]) / span
                                   * (grid[i] - grid[i - 1])))

    # Does the model itself evaluate there, without the chart wrapper?
    probes: Dict[str, Any] = {}
    for probe in (0.9 * r_minus, 0.7 * r_minus, 1.05 * real_limit):
        try:
            jets = model.pg_metric_jets(np.array([probe]))
            out_speed, in_speed = model.characteristic_speeds(np.array([probe]))
            probes[f"{probe:.6f}"] = {
                "jets_ok": True,
                "lapse": float(np.asarray(jets["N"])[0]),
                "shift": float(np.asarray(jets["beta"])[0]),
                "outgoing": float(out_speed[0]),
                "ingoing": float(in_speed[0]),
                "all_finite": bool(np.isfinite(np.asarray(jets["N"])[0])
                                   and np.isfinite(np.asarray(jets["beta"])[0])),
            }
        except Exception as exc:  # noqa: BLE001 - the point is to see if it raises
            probes[f"{probe:.6f}"] = {"jets_ok": False, "error": str(exc)}

    return {
        "mu": mu,
        "mass": mass,
        "charge_squared": charge2,
        "r_minus": r_minus,
        "r_plus": r_plus,
        "shift_real_domain_limit": real_limit,
        "chart_guard_limit": r_minus,
        "guard_is_beyond_the_real_domain": bool(real_limit < r_minus),
        "guard_reaches_beyond_the_geometry": bool(
            r_minus - real_limit > 0.1 * (r_plus - r_minus)),
        "samples": samples,
        "outgoing_zeroes": crossings,
        "zero_count": len(crossings),
        "trapped_between_the_two_zeroes": bool(
            len(crossings) == 2 and outgoing[0] > 0.0),
        "positive_below_inner_horizon": bool(
            float(outgoing[-1]) > 0.0),
        "model_probes": probes,
        "model_evaluates_below_inner_horizon": all(
            probe.get("jets_ok") for probe in probes.values()),
        "lapse_is_one_throughout": bool(
            abs(float(lapse.min()) - 1.0) < 1e-9 and abs(float(lapse.max()) - 1.0) < 1e-9),
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mu", nargs="+", type=float,
                        default=[1.002, 1.010, 1.040])
    parser.add_argument("--samples", type=int, default=20000)
    parser.add_argument("--out", default="nsc-rn-below-horizon-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    entries = [_entry(mu, args.samples) for mu in sorted(args.mu)]

    # What the chart refuses, stated in the chart's own terms.
    guard_message = None
    try:
        geometry = driver.case_geometry(args.mu[0])
        model = pg.action_model(geometry["mass"], geometry["r_m"])
        source.vacuum_chart(model, np.array([float(model.horizons["r_minus"]) * 0.9]))
    except Exception as exc:  # noqa: BLE001
        guard_message = str(exc)

    record = {
        "schema": "NSC-RN-BELOW-HORIZON-v1",
        "api_version": "nsc-rn-below-horizon-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "Whether the region below the inner horizon, which the chart refuses to "
            "sample, is actually describable by the model. The shift is real wherever "
            "r > Q^2/(2M), which is further in than r_minus, so the guard may be a "
            "convention rather than a property of the geometry."
        ),
        "chart_guard_message": guard_message,
        "production_guard_changed": False,
        "why_no_guard_was_relaxed": (
            "Nothing here relaxes a production guard. Whether the third region is "
            "reachable in a march is a separate question with a separate answer, and "
            "changing vacuum_chart on the strength of a reading would be the wrong order "
            "of operations."
        ),
        "scope_note": (
            "A construction check on the reference model. It does not measure the actual "
            "universe. Showing the third region is describable says nothing about whether "
            "it contains anything, and nothing about whether the nesting iterates."
        ),
        "guard_is_beyond_the_real_domain_everywhere": all(
            entry["guard_is_beyond_the_real_domain"] for entry in entries),
        "model_evaluates_below_inner_horizon_everywhere": all(
            entry["model_evaluates_below_inner_horizon"] for entry in entries),
        "two_zeroes_everywhere": all(entry["zero_count"] == 2 for entry in entries),
        "lapse_is_one_everywhere": all(entry["lapse_is_one_throughout"] for entry in entries),
        "entries": entries,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  chart guard says        : {guard_message}")
    for entry in entries:
        print(f"  {entry['run'] if 'run' in entry else 'mu=%.3f' % entry['mu']}"
              f"  guard {entry['chart_guard_limit']:.4f}"
              f" vs real limit {entry['shift_real_domain_limit']:.4f}"
              f"   zeroes {entry['zero_count']}"
              f"   model ok below: {entry['model_evaluates_below_inner_horizon']}")
        print(f"      outgoing zero crossings: "
              f"{['%.5f' % c for c in entry['outgoing_zeroes']]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())