#!/usr/bin/env python3
"""Write or replay the geometry-graded window of the initial three-packet compression.

The measurement reads the saved T=0.05 feedback payload. It does not evolve a
state and it does not rewrite that payload.
"""
from __future__ import annotations

import argparse

from recursive_horizons.nsc_geometry_graded_window import measure, verify_saved, write_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        report = verify_saved()
        print(
            report["status"],
            "cpu",
            report["cpu_seconds"],
            "leak",
            report["projection_leakage_max"],
            "discrete",
            report["discrete_max_abs"],
            "closed_six_mode",
            report["closed_six_mode"],
        )
        return 0
    payload = measure()
    path = write_record(payload)
    comparison = payload["comparison"]
    print(path)
    print(
        payload["status"],
        "cpu",
        payload["cpu_seconds"],
        "leak",
        comparison["projection_leakage_max_nf512"],
        "discrete",
        comparison["discrete_max_abs_nf256_vs_nf512"],
        "corner",
        comparison["corner_frobenius_nf512"],
        "H_gap",
        payload["historical_rational_example"]["H_gap"],
        "B_gap",
        payload["historical_rational_example"]["B_gap"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
