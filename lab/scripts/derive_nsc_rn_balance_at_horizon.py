"""Does the force balance coincide with the causal threshold?

The claim under test: gravity has to overcome the force of the space holding
it back, and at the point where it does, the object becomes a black hole. If
that is right then the mechanical threshold and the causal threshold are the
same radius, and the balance value there is exactly one half.

The two quantities are already computed in this line and have never been set
against each other:

* the **balance** ``B(r) = |F_in| / (|F_in| + |F_out|)``, from the radial and
  angular pressures of the sourced matter, where one half means the inward
  and outward forces are equal;
* the **trapping horizon**, the outermost simple zero of the outgoing speed
  in the evolved metric, which is where nothing inside can escape.

This evaluates ``B`` at the located trapping radius in every recorded case,
and reports the balance's own crossing radius alongside it so the two can be
compared directly rather than through a summary.

Reads only recorded station payloads. No time step is taken.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_balance_at_horizon.py
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
import derive_nsc_rn_parent as driver  # noqa: E402

DEVELOPMENT = LAB / "results/development"
STATION = re.compile(r"station-([0-9]+(?:\.[0-9]+)?)\.npz$")


def _git_head() -> str | None:
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=LAB,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def balance(radial: np.ndarray, angular: np.ndarray) -> np.ndarray:
    inward = np.abs(np.asarray(radial, dtype=float))
    outward = np.abs(np.asarray(angular, dtype=float))
    total = inward + outward
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(total > 0.0, inward / total, np.nan)


def crossings(radius: np.ndarray, values: np.ndarray) -> List[float]:
    out: List[float] = []
    for i in range(1, radius.size):
        a, b = values[i - 1], values[i]
        if not (np.isfinite(a) and np.isfinite(b)):
            continue
        if (a - 0.5) * (b - 0.5) < 0.0:
            span = b - a
            out.append(float(radius[i - 1] + (0.5 - a) / span * (radius[i] - radius[i - 1])))
    return out


def _reduce(directory: Path, mu: float) -> Dict[str, Any]:
    geometry = driver.case_geometry(mu)
    model = pg.action_model(geometry["mass"], geometry["r_m"])
    r_plus = float(model.horizons["r_plus"])

    stations = sorted((p for p in directory.glob("station-*.npz")
                       if STATION.search(p.name)),
                      key=lambda p: float(STATION.search(p.name).group(1)))
    rows: List[Dict[str, Any]] = []
    for station in stations:
        with np.load(station) as payload:
            if "radial_pressure" not in payload.files:
                return {"run": directory.name, "assessable": False,
                        "reason": "station carries no pressure channels"}
            radius = np.asarray(payload["radius"], dtype=float)
            values = balance(payload["radial_pressure"], payload["angular_pressure"])
            located = obs.trapping_horizon(
                radius,
                np.asarray(payload["frozen_lapse"], dtype=float),
                np.asarray(payload["frozen_shift"], dtype=float))
            trapping = located["trapping_horizon"]
            clock = float(payload["coupled_t"][0])

        at_horizon = None
        if trapping is not None and trapping >= radius[0]:
            at_horizon = float(np.interp(trapping, radius, values))
        near = crossings(radius, values)
        rows.append({
            "t_over_rm": clock,
            "trapping_horizon": trapping,
            "balance_at_trapping_horizon": at_horizon,
            "balance_min": float(np.nanmin(values)),
            "balance_max": float(np.nanmax(values)),
            "balance_crossings": near,
            "nearest_crossing_to_trapping": (min(near, key=lambda x: abs(x - trapping))
                                             if near and trapping is not None else None),
            "trapping_inside_domain": bool(trapping is not None and trapping >= radius[0]),
        })

    usable = [r for r in rows if r["balance_at_trapping_horizon"] is not None]
    values_at = [r["balance_at_trapping_horizon"] for r in usable]
    return {
        "run": directory.name,
        "mu": mu,
        "r_plus": r_plus,
        "stations": len(rows),
        "stations_assessable": len(usable),
        "assessable": bool(usable),
        "balance_at_trapping_horizon": values_at,
        "balance_mean_at_trapping_horizon": (float(np.mean(values_at)) if values_at else None),
        "balance_spread_at_trapping_horizon": (max(values_at) - min(values_at)
                                              if values_at else None),
        "distance_to_half_mean": (abs(float(np.mean(values_at)) - 0.5)
                                  if values_at else None),
        "trapping_horizon_constant": len({round(r["trapping_horizon"], 9)
                                          for r in rows
                                          if r["trapping_horizon"] is not None}) == 1,
        "rows": rows,
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="+", required=True,
                        help="RUN:MU pairs under results/development")
    parser.add_argument("--out", default="nsc-rn-balance-at-horizon-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    runs = []
    for pair in args.runs:
        name, _, mu = pair.partition(":")
        directory = DEVELOPMENT / name
        if not directory.is_dir():
            print(f"no run at {directory}", file=sys.stderr)
            return 1
        runs.append(_reduce(directory, float(mu or 1.04)))

    assessable = [r for r in runs if r.get("assessable")]
    means = [r["balance_mean_at_trapping_horizon"] for r in assessable]
    record = {
        "schema": "NSC-RN-BALANCE-AT-HORIZON-v1",
        "api_version": "nsc-rn-balance-at-horizon-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "claim_under_test": (
            "Gravity has to overcome the force of the space holding it back, and at the "
            "point where it does the object becomes a black hole. If that is right, the "
            "mechanical threshold and the causal threshold are the same radius and the "
            "balance there is exactly one half."
        ),
        "balance_definition": "B(r) = |radial_pressure| / (|radial_pressure| + |angular_pressure|)",
        "why_half_matters": (
            "B = 0.5 is where the inward and outward forces are equal, which is the "
            "mechanical statement of the claim. The trapping horizon is the causal "
            "statement. If they agree, the two are the same event."
        ),
        "scope_note": (
            "A construction check on recorded finite states. It does not measure the "
            "actual universe. Agreement here would show the two criteria coincide in "
            "this construction, which is not evidence that black holes in general form at "
            "a force balance."
        ),
        "runs_assessable": len(assessable),
        "runs_total": len(runs),
        "balance_mean_over_runs": means,
        "balance_mean_of_means": float(np.mean(means)) if means else None,
        "distance_of_mean_from_half": (abs(float(np.mean(means)) - 0.5)
                                       if means else None),
        "coincides_with_half": bool(means and abs(float(np.mean(means)) - 0.5) < 0.05),
        "runs": runs,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    print(f"  assessable runs : {len(assessable)} of {len(runs)}")
    for run in assessable:
        print(f"  {run['run']:<44} B(trapping) = "
              f"{run['balance_mean_at_trapping_horizon']:.4f}"
              f"  spread {run['balance_spread_at_trapping_horizon']:.4f}")
    print(f"  mean of means   : {record['balance_mean_of_means']:.4f}")
    print(f"  coincides with 0.5: {record['coincides_with_half']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())