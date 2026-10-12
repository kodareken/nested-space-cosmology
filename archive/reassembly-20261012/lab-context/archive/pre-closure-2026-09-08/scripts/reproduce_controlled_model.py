#!/usr/bin/env python3
"""Write the versioned machine record for the v0.11 controlled model."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.evidence import reproduce_controlled_model  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPOSITORY / "results" / "controlled-model.json",
        help="JSON output path",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    record = reproduce_controlled_model()
    payload = json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(output)

    core = record["covariant_core"]
    information = record["information_map"]
    perturbations = record["linear_perturbations"]
    domain = record["domain_synthesis"]
    discriminator = record["independent_discriminator"]
    gates = record["literal_model_gates"]
    print(f"wrote {output}")
    print(
        "CCT-1 maximum Friedmann residual: "
        f"{core['max_abs_friedmann_residual']:.3e}"
    )
    print(
        "IM-1 isometry error: "
        f"{information['audit']['isometry_error']:.3e} "
        "(finite code toy; not gravitational factorization)"
    )
    print(
        "tensor final determinant error: "
        f"{perturbations['tensor_refinement']['final_determinant_error']:.3e}"
    )
    print(
        "DST-1 compensated spectator endpoint residual: "
        f"{domain['fixed_h_spectator_endpoint']['endpoint_residual']:.3e} "
        "(effective-fluid benchmark; not dark-sector identification)"
    )
    print(
        "Kottler balance acceleration residual: "
        f"{domain['kottler_balance']['acceleration_at_balance_m_s^-2']:.3e} m/s^2 "
        "(unstable weak-field sign change; not a domain boundary)"
    )
    print(
        "conditional globally closed-branch discriminator: "
        f"{discriminator['prediction']} (independent sign given the added branch "
        "postulate; not confirmation)"
    )
    print(
        "literal-model gate controls/rejections: "
        f"{gates['classification']} "
        "(this controlled record retains published-summary gates; the separate "
        "BAO artifact evaluates the real distance likelihood)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
