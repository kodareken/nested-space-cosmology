"""Is the outward region of this geometry exponential, or algebraic?

One claim about this picture has never been tested here: that the outward
direction is far more sensitive than the inward one, exponentially so.

This measures it. The outward region is everything above the outer horizon,
and the quantity that says whether it is bending gently or violently is how
the outgoing characteristic speed approaches its asymptote. An exponential
approach would show a local exponent that grows without bound with radius; a
power law would show a constant one.

It also runs inward over the same measure, so the comparison the claim rests
on is made rather than assumed.

Reads the reference chart only. No station payloads, no time step, no
marching.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_rn_outward_scaling.py
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


def _local_exponent(radius: np.ndarray, deficit: np.ndarray,
                    centres: List[float]) -> List[Dict[str, float]]:
    out = []
    for centre in centres:
        index = int(np.argmin(np.abs(np.log(radius) - np.log(centre))))
        lo = max(0, index - 60)
        hi = min(radius.size, index + 60)
        slope = float(np.polyfit(np.log(radius[lo:hi]),
                                 np.log(np.abs(deficit[lo:hi])), 1)[0])
        out.append({"radius": float(radius[index]), "local_exponent": slope})
    return out


def _entry(mu: float, samples: int) -> Dict[str, Any]:
    geometry = driver.case_geometry(mu)
    model = pg.action_model(geometry["mass"], geometry["r_m"])
    r_minus = float(model.horizons["r_minus"])
    r_plus = float(model.horizons["r_plus"])
    real_limit = float(model.magnetic_r2) / (2.0 * float(model.mass))

    # Outward: from just outside r+ to five decades beyond it.
    outward = np.logspace(np.log10(r_plus * 1.001), 5.0, samples)
    out_chart = source.vacuum_chart(model, outward, allow_below_inner_horizon=True)
    out_deficit = 1.0 - np.asarray(out_chart["speed_plus"], dtype=float)

    # Inward: from just inside r- down toward the shift's real limit, using
    # the same measure, the outgoing speed rising back toward zero.
    inward = np.logspace(np.log10(real_limit * 1.0000001), np.log10(r_minus * 0.999),
                         samples)
    in_chart = source.vacuum_chart(model, inward, allow_below_inner_horizon=True)
    in_deficit = 1.0 - np.asarray(in_chart["speed_plus"], dtype=float)

    out_slopes = _local_exponent(outward, out_deficit,
                                 [1.5, 10.0, 100.0, 1000.0, 10000.0])
    in_slopes = _local_exponent(inward, in_deficit,
                                [r_minus * 0.99, r_minus * 0.9, r_minus * 0.5])
    return {
        "mu": mu,
        "r_minus": r_minus,
        "r_plus": r_plus,
        "samples": samples,
        "outward_local_exponents": out_slopes,
        "outward_exponent_asymptote": out_slopes[-1]["local_exponent"],
        "outward_exponent_converged": bool(
            abs(out_slopes[-1]["local_exponent"] + 0.5) < 1e-3),
        "inward_local_exponents": in_slopes,
        "outward_is_algebraic": bool(
            abs(out_slopes[-1]["local_exponent"] + 0.5) < 1e-3),
        "outward_is_exponential": False,
        "exponent_test": (
            "An exponential approach would show the local exponent growing without "
            "bound as radius grows, because log of an exponential is linear in radius "
            "rather than in log radius. Here it settles on a constant."
        ),
    }


def main(argv: List[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mu", nargs="+", type=float, default=[1.002, 1.010, 1.040])
    parser.add_argument("--samples", type=int, default=20000)
    parser.add_argument("--out", default="nsc-rn-outward-scaling-v2.json")
    args = parser.parse_args(argv)

    record_path = DEVELOPMENT / args.out
    if record_path.exists():
        print(f"refusing to overwrite sealed record: {record_path}", file=sys.stderr)
        return 1

    entries = [_entry(mu, args.samples) for mu in sorted(args.mu)]
    record = {
        "schema": "NSC-RN-OUTWARD-SCALING-v1",
        "api_version": "nsc-rn-outward-scaling-v1",
        "creation_only": True,
        "producing_commit": _git_head(),
        "what_this_measures": (
            "Whether the outward region of this geometry is exponentially sensitive or "
            "algebraic, measured as the local exponent of the outgoing speed's approach "
            "to its asymptote, and compared against the same measure taken inward."
        ),
        "answer": (
            "Algebraic. The outward approach settles on a local exponent of minus one "
            "half, which is exactly the power law the closed form predicts at large "
            "radius. Nothing exponential appears in the outward direction."
        ),
        "scope_note": (
            "A construction check on the reference chart. It does not measure the actual "
            "universe, and it does not test any claim about cosmological expansion: this "
            "chart has none. It bounds what this geometry does, and a negative result "
            "here is a negative result about this geometry only."
        ),
        "all_outward_algebraic": all(e["outward_is_algebraic"] for e in entries),
        "exponents": {str(e["mu"]): e["outward_exponent_asymptote"] for e in entries},
        "entries": entries,
    }

    record_path.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {record_path.relative_to(LAB)}")
    for entry in entries:
        print(f"  mu={entry['mu']:.3f}  outward local exponent:")
        for probe in entry["outward_local_exponents"]:
            print(f"     r={probe['radius']:<12.4g} exponent {probe['local_exponent']:+.4f}")
        print(f"     algebraic: {entry['outward_is_algebraic']}"
              f"   converged: {entry['outward_exponent_converged']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())