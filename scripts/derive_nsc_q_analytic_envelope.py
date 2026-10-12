"""Analytic envelope of the room-stretch power.

Round-67 showed p_Q is the room stretch 32/r_inner. The analytic
derivative of that stretch follows: with the inner edge one cell above
r_minus, the local gap power is x/mu with x = sqrt(mu^2 - 1), so the
chord power measured across any pair must lie inside the envelope
[x/mu over the family]. This driver evaluates the envelope, verifies
the measured chord sits inside it, and checks the implied range for
4*p_Q against the ladder's b (the round-58 identity): 4*p over the
envelope is [0.25, 1.10], and the recorded b = 0.4672 must lie inside.

    python scripts/lab.py scripts/derive_nsc_q_analytic_envelope.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

LADDER = LAB / "results" / "development" / "nsc-charged-ladder-deficit-law-v1.json"
COVARIANCE = LAB / "results" / "development" / "nsc-two-ends-covariance-v1.json"
DECOMPOSITION = LAB / "results" / "development" / "nsc-q-decomposition-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-q-analytic-envelope-v1.json"
SCHEMA = "NSC-Q-ANALYTIC-ENVELOPE-v1"
MASSES = (1.002, 1.010, 1.040)


def main():
    if any(not p.is_file() for p in (LADDER, COVARIANCE, DECOMPOSITION)):
        raise SystemExit("missing input record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    ladder = json.loads(LADDER.read_text(encoding="utf-8"))
    covariance = json.loads(COVARIANCE.read_text(encoding="utf-8"))
    decomposition = json.loads(DECOMPOSITION.read_text(encoding="utf-8"))
    envelope = []
    for mass in MASSES:
        x = math.sqrt(mass * mass - 1.0)
        envelope.append({"mu": mass, "x": x, "x_over_mu": x / mass})
    low = envelope[0]["x_over_mu"]
    high = envelope[-1]["x_over_mu"]
    measured = covariance["q_gap_power_post_hoc"]
    b = ladder["b"]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "envelope": envelope,
        "envelope_range": [low, high],
        "measured_chord_p_q": measured,
        "chord_inside_envelope": bool(low <= measured <= high),
        "four_p_envelope": [4.0 * low, 4.0 * high],
        "ladder_b": b,
        "b_inside_four_p_envelope": bool(4.0 * low <= b <= 4.0 * high),
        "reading": (
            "the local power of the stretch is x/mu; any chord over the "
            "family lands inside its range, and the identity 4*p_Q = b is "
            "achievable only if b sits inside the four-fold envelope"
        ),
        "notes": [
            "the precise chord depends on the march's per-case inner-edge "
            "construction; the envelope only bounds it",
            "this is constraint satisfaction, not a derivation of b",
        ],
        "input_sha256": {
            "lab/results/development/nsc-charged-ladder-deficit-law-v1.json":
                chain.path_sha256(LADDER),
            "lab/results/development/nsc-two-ends-covariance-v1.json":
                chain.path_sha256(COVARIANCE),
            "lab/results/development/nsc-q-decomposition-v1.json":
                chain.path_sha256(DECOMPOSITION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_q_analytic_envelope.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_q_analytic_envelope.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for entry in envelope:
        print("  mu {:<5} x {:.4f}  x/mu {:.4f}".format(entry["mu"], entry["x"],
                                                       entry["x_over_mu"]))
    print("  envelope [{:.4f}, {:.4f}]  measured chord {:.4f} inside {}".format(
        low, high, measured, record["chord_inside_envelope"]))
    print("  4p envelope [{:.4f}, {:.4f}]  ladder b {:.4f} inside {}".format(
        4.0 * low, 4.0 * high, b, record["b_inside_four_p_envelope"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())