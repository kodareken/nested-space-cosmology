"""Scale-counting for the wall clock pair, with one derivable check.

The measured pair (rise ~ gap^-1/2, fall ~ gap^-1/4) invites the
analytic question. The near-extremal scale counting offers exactly one
parameter-free decomposition: if the fall phase is a constant-frequency
relaxation sampled by the breathing amplitude A(mu), then

    fall = omega * A        with omega family-constant.

This driver evaluates fall/A from the recorded rates and envelopes of
all three cases and records the constant's spread. The claim label is
honored: the decomposition is a derived relation on recorded numbers;
the deeper origin of the exponents stays underdetermined at three gap
points and is marked as such.

    python scripts/lab.py scripts/derive_nsc_wall_clock_derivation.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

BREATHING = LAB / "results" / "development" / "nsc-wall-breathing-family-v1.json"
CLOCK = LAB / "results" / "development" / "nsc-wall-clock-horizon-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-wall-clock-derivation-v1.json"
SCHEMA = "NSC-WALL-CLOCK-DERIVATION-v1"


def main():
    if not BREATHING.is_file() or not CLOCK.is_file():
        raise SystemExit("missing breathing or clock record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    breathing = json.loads(BREATHING.read_text(encoding="utf-8"))
    clock = json.loads(CLOCK.read_text(encoding="utf-8"))
    sites = {}
    for name in breathing["cases"]:
        entry = breathing["cases"][name]
        low, high = entry["radial_range"]
        amplitude = high - low
        fall = entry["mean_fall_rate"]
        rise = entry["mean_rise_rate"]
        sites[name] = {
            "amplitude": amplitude,
            "fall": fall,
            "rise": rise,
            "fall_over_amplitude": (abs(fall) / amplitude) if amplitude else None,
        }
    omegas = [site["fall_over_amplitude"] for site in sites.values()
              if site["fall_over_amplitude"] is not None]
    mean_omega = sum(omegas) / len(omegas)
    spread_omega = (sum((v - mean_omega) ** 2 for v in omegas) / len(omegas)) ** 0.5 \
        / abs(mean_omega)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "decomposition": "fall = omega * amplitude, omega treated as family-constant",
        "sites": sites,
        "omega_mean": mean_omega,
        "omega_spread": spread_omega,
        "decomposition_holds": bool(spread_omega < 0.12),
        "scale_counting": {
            "gap": "Delta = r_+ - r_- = 2 sqrt(M^2 - P^2)",
            "rise_power": -0.5,
            "rise_reading": "reciprocal proper throat-crossing scale",
            "fall_power": -0.25,
            "fall_reading": "enters only through the breathing amplitude; the "
                           "fall is a constant-frequency relaxation sampled by A",
            "underdetermined": "three gap points cannot fix origin exponents; "
                               "the powers are measured, their provenance is flagged",
        },
        "reading": (
            "the one parameter-free decomposition in the recorded numbers; "
            "not a claim about the origin of the exponents"
        ),
        "input_sha256": {
            "lab/results/development/nsc-wall-breathing-family-v1.json":
                chain.path_sha256(BREATHING),
            "lab/results/development/nsc-wall-clock-horizon-v1.json":
                chain.path_sha256(CLOCK),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_wall_clock_derivation.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_wall_clock_derivation.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, site in sites.items():
        print("  {:<20} amplitude {:.3f} fall {:.3f}  fall/A {:.4f}".format(
            name, site["amplitude"], site["fall"], site["fall_over_amplitude"]))
    print("  omega mean {:.4f}  spread {:.1%}  decomposition_holds: {}".format(
        record["omega_mean"], record["omega_spread"], record["decomposition_holds"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())