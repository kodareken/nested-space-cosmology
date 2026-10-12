#!/usr/bin/env python3
"""Write the versioned, machine-readable core calculation record."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons import CosmologyInputs, reproduce_core  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPOSITORY / "results" / "core-identities.json",
        help="JSON output path",
    )
    parser.add_argument("--h0", type=float, default=67.4)
    parser.add_argument("--omega-m", type=float, default=0.315)
    parser.add_argument("--omega-lambda", type=float, default=0.685)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    inputs = CosmologyInputs(args.h0, args.omega_m, args.omega_lambda)
    record = reproduce_core(inputs)
    payload = json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n"

    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(output)

    standard = record["standard_cosmology"]
    conditional = record["conditional_horizon_saturation"]
    print(f"wrote {output}")
    print(f"Hubble-sphere compactness identity: {standard['hubble_sphere_compactness']:.16g}")
    print(f"Lambda conversion: {standard['lambda_m^-2']:.8e} m^-2")
    print(f"Matter+Lambda age benchmark: {standard['matter_plus_lambda_age_Gyr']:.5f} Gyr")
    print(
        "H-SAT parent mass inference: "
        f"{conditional['parent_mass_kg']:.8e} kg "
        "(conditional; not an independent prediction)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
