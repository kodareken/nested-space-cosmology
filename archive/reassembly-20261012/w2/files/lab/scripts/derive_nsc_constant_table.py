"""The programme's dimensionless constants, laid in one table.

Across fifty-eight rounds the programme measured a set of dimensionless
coefficients from different regimes: the room-sharpness gap power p_Q,
the kappa-ladder Gaussian pair (a, b), the wall oscillator omega, the
clock powers, and the covariance exponent. This driver assembles them
from their own records and records the one near-identity visible in the
assembly: 4 * p_Q = 0.4692 versus the ladder's quadratic coefficient
b = 0.4672 (0.4% difference, numerology status flagged). No mechanism
is claimed.

    python scripts/lab.py scripts/derive_nsc_constant_table.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

COVARIANCE = LAB / "results" / "development" / "nsc-two-ends-covariance-v1.json"
LADDER = LAB / "results" / "development" / "nsc-charged-ladder-deficit-law-v1.json"
DERIVATION = LAB / "results" / "development" / "nsc-wall-clock-derivation-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-constant-table-v1.json"
SCHEMA = "NSC-CONSTANT-TABLE-v1"


def main():
    if any(not p.is_file() for p in (COVARIANCE, LADDER, DERIVATION)):
        raise SystemExit("missing constant-source record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    covariance = json.loads(COVARIANCE.read_text(encoding="utf-8"))
    ladder = json.loads(LADDER.read_text(encoding="utf-8"))
    derivation = json.loads(DERIVATION.read_text(encoding="utf-8"))
    p_q = covariance["q_gap_power_post_hoc"]
    b_ladder = ladder["b"]
    a_ladder = ladder["a"]
    omega = derivation["omega_mean"]
    constants = {
        "p_Q_gap_power": {"value": p_q, "source": "nsc-two-ends-covariance-v1"},
        "ladder_a": {"value": a_ladder, "source": "nsc-charged-ladder-deficit-law-v1"},
        "ladder_b": {"value": b_ladder, "source": "nsc-charged-ladder-deficit-law-v1"},
        "wall_oscillator_omega": {"value": omega, "source": "nsc-wall-clock-derivation-v1"},
        "rise_gap_power": {"value": 0.5, "source": "pre-registered round 26"},
        "fall_gap_power": {"value": 0.25, "source": "pre-registered round 30"},
        "covariance_composed_exponent": {"value": covariance["composed_exponent"],
                                         "source": "nsc-two-ends-covariance-v1"},
    }
    identity = {
        "four_times_p_Q": 4.0 * p_q,
        "ladder_b": b_ladder,
        "relative_difference": (4.0 * p_q - b_ladder) / b_ladder,
    }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "constants": constants,
        "near_identity": identity,
        "status": (
            "observation-level; the 4*p_Q ~ b proximity is flagged as "
            "possible numerology, exactly like the round-1 transfer share"
        ),
        "notes": [
            "every value is re-read from its own record rather than "
            "remembered; sources are named per constant",
        ],
        "input_sha256": {
            "lab/results/development/nsc-two-ends-covariance-v1.json":
                chain.path_sha256(COVARIANCE),
            "lab/results/development/nsc-charged-ladder-deficit-law-v1.json":
                chain.path_sha256(LADDER),
            "lab/results/development/nsc-wall-clock-derivation-v1.json":
                chain.path_sha256(DERIVATION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_constant_table.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_constant_table.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for name, entry in constants.items():
        print("  {:<32} {:.6f}  ({})".format(name, entry["value"], entry["source"]))
    print("  4*p_Q = {:.6f}  ladder b = {:.6f}  rel diff {:+.2%}".format(
        identity["four_times_p_Q"], identity["ladder_b"],
        identity["relative_difference"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())