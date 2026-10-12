"""The measured chain fee meets the Schroeder skeleton.

The matryoshka record anchors the NSC chain to the small-Schroeder
skeleton with per-link coefficient 1/Omega (Omega = 4). This driver
joins that skeleton to the first-principles fee programme: with path
branching b = 3 (the recorded Schroeder asymptotic), a chain of
per-level fee f sums to Sigma = 1/(1 - 3f). It evaluates Sigma for the
three recorded fees - the structural resolvent 1/Omega = 0.2517, the
round-10 lapse tangent 0.2610, and the round-1 transfer 0.0446 - and
holds each against the Omega = 4 anchor. Parameter-free: no fit.

    python scripts/lab.py scripts/derive_nsc_matryoshka_fee_join.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

MATRYOSHKA = LAB / "results" / "development" / "nsc-matryoshka-comparison-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-matryoshka-fee-join-v1.json"
SCHEMA = "NSC-MATRYOSHKA-FEE-JOIN-v1"
BRANCHING = 3.0
OMEGA_ANCHOR = 4.0
FEES = {
    "resolvent_1_over_omega": 0.2517,
    "lapse_tangent_f": 0.2610,
    "transfer_fee": 0.0446,
}


def main():
    if not MATRYOSHKA.is_file():
        raise SystemExit("missing matryoshka record: " + str(MATRYOSHKA))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    matryoshka = json.loads(MATRYOSHKA.read_text(encoding="utf-8"))
    sites = {}
    for name, fee in FEES.items():
        per_level_growth = BRANCHING * fee
        if per_level_growth < 1.0:
            sigma = 1.0 / (1.0 - per_level_growth)
        else:
            sigma = None
        sites[name] = {
            "fee": fee,
            "per_level_growth": per_level_growth,
            "sigma": sigma,
            "sigma_distance_to_omega4": None if sigma is None else sigma - OMEGA_ANCHOR,
            "sigma_relative": None if sigma is None else sigma / OMEGA_ANCHOR,
        }
    schroder = matryoshka["comparison"]["schroder_computed"]
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "branching": BRANCHING,
        "schroder_counts": schroder,
        "sites": sites,
        "reading": (
            "with the recorded Schroeder branching b=3, a per-level fee f "
            "sums to 1/(1-3f); the Omega=4 resolvent anchor of the "
            "matryoshka record is the reference"
        ),
        "notes": [
            "branching 3 is the recorded Schroeder asymptotic, not a fit",
            "the lapse tangent fee 0.2610 and the transfer fee 0.0446 are "
            "recorded first-principles measurements",
        ],
        "input_sha256": {
            "lab/results/development/nsc-matryoshka-comparison-v1.json":
                chain.path_sha256(MATRYOSHKA),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_matryoshka_fee_join.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_matryoshka_fee_join.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  Schroeder counts:", schroder, " branching 3")
    for name, site in sites.items():
        sigma_text = "diverges" if site["sigma"] is None else "{:.4f}".format(site["sigma"])
        dist_text = "-" if site["sigma_distance_to_omega4"] is None \
            else "{:+.4f} ({:+.1%})".format(site["sigma_distance_to_omega4"], site["sigma_relative"] - 1.0)
        print("  {:<22} growth {:.4f} Sigma {} vs Omega=4: {}".format(
            name, site["per_level_growth"], sigma_text, dist_text))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())