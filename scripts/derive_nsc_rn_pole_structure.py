"""Why the integrator stops at the inner horizon.

The last commit described the obstruction below the inner horizon as
stiffness. That was an assertion. This measures it.

The characteristic right-hand side divides by the two speeds:

    d_plus  = (...)/speed_plus
    d_minus = (...)/speed_minus

and ``speed_plus`` vanishes wherever the shift equals the lapse, which is
at **both** horizons, not only at ``r_plus``. So the right-hand side carries
a simple pole at ``r_minus`` as well. The code clears the pole at ``r_plus``
by starting the branch from a horizon jet; it does not clear the one at
``r_minus``, and so integrates straight into it.

This locates both poles from the chart, measures the speed's approach to
zero on each side, and shows that the inward branch is heading for the second
one. It calls the model's own functions and reads no station payloads.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_pole_structure.py
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


def _zeroes(radius: np.ndarray, values: np.ndarray) -> List[float]:
    out = []
    for i in range(1, radius.size):
        if values[i - 1] * values[i] < 0.0:
            span = values[i] - values[i - 1]
            out.append(float(radius[i - 1]
                             + (-values[i - 1]) / span * (radius[i] - radius[i - 1])))
    return out


def _entry(mu: float, samples: int) -> Dict[str, Any]:
    geometry = driver.case_geometry(mu)
    model = pg.action_model(geometry["mass"], geometry["r_m"])
    r_minus = float(model.horizons["r_minus"])
    r_plus = float(model.horizons["r_plus"])
    real_limit = float(model.magnetic_r2) / (2.0 * float(model.mass))

    grid = np.linspace(real_limit * (1.0 + 1e-9), 12.0, samples)
    outgoing = np.asarray(source.vacuum_chart(
        model, grid, allow_below_inner_horizon=True)["speed_plus"], dtype=float)
    ingoing = np.asarray(source.vacuum_chart(
        model, grid, allow_below_inner_horizon=True)["speed_minus"], dtype=float)

    out_zeroes = _zeroes(grid, outgoing)
    in_zeroes = _zeroes(grid, ingoing)

    # How fast the denominator collapses on each side of r_minus, which is
    # what sets how hard the integrator is pushed.
    probes = []
    for delta in (1e-1, 1e-2, 1e-3, 1e-4):
        probes.append({
            "offset_above_r_minus": delta,
            "speed_plus_just_above": float(np.interp(r_minus + delta, grid, outgoing)),
            "reciprocal": float(1.0 / np.interp(r_minus + delta, grid, outgoing)),
        })

    return {
        "mu": mu,
        "r_minus": r_minus,
        "r_plus": r_plus,
        "shift_real_limit": real_limit,
        "samples": samples,
        "outgoing_zeroes": out_zeroes,
        "outgoing_zero_count": len(out_zeroes),
        "ingoing_zeroes": in_zeroes,
        "ingoing_never_vanishes": len(in_zeroes) == 0,
        "both_poles_are_on_the_outgoing_leg": bool(
            len(out_zeroes) == 2 and len(in_zeroes) == 0),
        "second_pole_at_r_minus": any(abs(z - r_minus) < 1e-3 for z in out_zeroes),
        "first_pole_at_r_plus": any(abs(z - r_plus) < 1e-3 for z in out_zeroes),
        "approach_to_r_minus": probes,
        "reciprocal_at_one_thousandth_above": probes[2]["reciprocal"],
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mu", nargs="+", type=float, default=[1.002, 1.010, 1.040])
    parser.add_argument("--samples", type=int, default=40000)
    parser.add_argument("--out", default="nsc-rn-pole-structure-v1.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    entries = [_entry(mu, args.samples) for mu in sorted(args.mu)]
    record = {
        "schema": "NSC-RN-POLE-STRUCTURE-v1",
        "api_version": "nsc-rn-pole-structure-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "Where the characteristic right-hand side is singular. It divides by the two "
            "speeds, and the outgoing speed vanishes at both horizons, so the equation "
            "carries a pole the inward branch integrates into."
        ),
        "mechanism": (
            "The code clears the pole at r_plus by starting the branch from a horizon jet "
            "rather than by dividing at the horizon. It does not do the same at r_minus, "
            "and that is where the inward integration stops."
        ),
        "correction": (
            "The previous commit called this stiffness. It is not stiffness in the solver's "
            "sense; it is a simple pole in the right-hand side at a known radius, and the "
            "solver reports it as a step-size failure because the solution is unbounded "
            "there."
        ),
        "scope_note": (
            "A construction check on the reference chart. It does not measure the actual "
            "universe. Locating the pole says where the obstruction is, not what is on "
            "the other side of it."
        ),
        "two_poles_everywhere": all(e["both_poles_are_on_the_outgoing_leg"]
                                    for e in entries),
        "second_pole_at_r_minus_everywhere": all(e["second_pole_at_r_minus"]
                                                 for e in entries),
        "ingoing_never_vanishes": all(e["ingoing_never_vanishes"] for e in entries),
        "entries": entries,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    for entry in entries:
        print(f"  mu={entry['mu']:.3f}  r_minus={entry['r_minus']:.6f}  "
              f"r_plus={entry['r_plus']:.6f}")
        print(f"     outgoing zero crossings: {['%.6f' % z for z in entry['outgoing_zeroes']]}"
              f"   ingoing: {len(entry['ingoing_zeroes'])}")
        print(f"     1/speed_plus just above r_minus:")
        for probe in entry["approach_to_r_minus"]:
            print(f"        +{probe['offset_above_r_minus']:.0e}  speed {probe['speed_plus_just_above']:+.3e}"
                  f"   reciprocal {probe['reciprocal']:+.4e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())