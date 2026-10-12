#!/usr/bin/env python3
"""Write the bounded, real DESI DR2 BAO-only likelihood profile artifact."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.bao import reproduce_desi_dr2_bao  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPOSITORY / "results" / "desi-dr2-bao-profile.json")
    parser.add_argument("--simpson-intervals", type=int, default=256)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    record = reproduce_desi_dr2_bao(intervals=args.simpson_intervals)
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(output)
    best = record["best_fit"]
    print(f"wrote {output}")
    print(f"BAO-only best fit: w={best['w']:.6f}, Omega_m={best['omega_m']:.6f}, H0*r_d={best['H0rd_km_s']:.3f} km/s, chi2={best['chi2']:.6f}")
    print("Scope: real Gaussian BAO distance likelihood only; not CMB/SN, full shape, or external-origin evidence.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
