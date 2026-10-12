"""Close the nested-influence chain from recorded NSC bytes. No evolution.

Inputs are existing, hash-bound records:
    - results/development/scale-binding.json        (Omega, transfer)
    - results/development/nsc-rn-neutral-continuation-v1/summary.json
                                                    (measured mass fee)

For each candidate per-level coefficient the self-similar chain is summed
and the local fraction is compared against the pre-registered band around
the BBN-anchored baryon share (0.049). The record is creation-only and
refuses to overwrite an existing output.

    python scripts/lab.py scripts/derive_nsc_influence_chain.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

SCALE_BINDING = LAB / "results" / "development" / "scale-binding.json"
CONTINUATION = (
    LAB / "results" / "development" / "nsc-rn-neutral-continuation-v1" / "summary.json"
)
OUTPUT = LAB / "results" / "development" / "nsc-influence-chain-v1.json"
SCHEMA = "NSC-INFLUENCE-CHAIN-v1"


def _load_json(path):
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def main():
    for path in (SCALE_BINDING, CONTINUATION):
        if not path.is_file():
            raise SystemExit("missing input record: " + str(path))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))

    scale_record = _load_json(SCALE_BINDING)
    continuation_record = _load_json(CONTINUATION)

    scale = chain.scale_binding_fees(scale_record)
    mass = chain.continuation_mass_fee(continuation_record)

    fees = {
        "structural_resolvent_one_over_omega": scale["structural_one_over_omega"],
        "recorded_transfer": scale["recorded_transfer"],
        "measured_mass_fee": mass["measured_relative_mass_renormalisation"],
    }
    rows = chain.candidate_table(fees)
    proximity = chain.transfer_baryon_proximity(fees)
    implied = {
        "from_recorded_transfer": chain.matrix_coefficient(fees["recorded_transfer"]),
        "from_measured_mass_fee": chain.matrix_coefficient(fees["measured_mass_fee"]),
    }

    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "band": [chain.BAND_LO, chain.BAND_HI],
        "target_local_fraction": chain.TARGET_LOCAL_FRACTION,
        "target_anchor": "BBN-anchored baryon share Omega_b ~ 0.049, Planck 2018 frame",
        "candidates": rows,
        "transfer_baryon_proximity": proximity,
        "implied_layer_ratio": implied,
        "route_ratio": chain.ROUTE_RATIO,
        "route_ratio_source": "ROUTE_OUTER_OVER_RM = 32.0 in lab/scripts/derive_nsc_rn_parent.py",
        "measured_mass_fee_detail": mass,
        "scale_fee_detail": scale,
        "input_sha256": {
            str(SCALE_BINDING.relative_to(LAB)): chain.path_sha256(SCALE_BINDING),
            str(CONTINUATION.relative_to(LAB)): chain.path_sha256(CONTINUATION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_influence_chain.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_influence_chain.py"
            ),
        },
        "notes": [
            "No coefficient is fitted anywhere; the three candidate fees are "
            "extracted from frozen recorded bytes and the band is pre-registered.",
            "Proposal, not observation: the chain closure is a model comparison, "
            "not a cosmology measurement.",
            "The recorded-transfer proximity to 0.049 is an observation about two "
            "fixed numbers and is flagged as possible numerology.",
        ],
    }
    encoded = json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    OUTPUT.write_text(encoded, encoding="utf-8")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print(
            "  {:<38} f={: .6e}  local={: .6f}  vertex={}".format(
                row["candidate"], row["f"], row["local_fraction"], row["verdict"]
            )
        )
    print(
        "  recorded_transfer vs baryon target: offset {:+g} ({:+.2%})".format(
            proximity["absolute_offset"], proximity["relative_offset"]
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())