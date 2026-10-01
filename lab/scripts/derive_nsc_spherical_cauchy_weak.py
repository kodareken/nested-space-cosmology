#!/usr/bin/env python3
"""Write the weak initial-residual record from the saved v5 seed.

Does not evolve, retune a tolerance, or rewrite the v5 JSON or NPZ.
"""
from recursive_horizons.nsc_spherical_cauchy_weak import build_record


def main():
    record = build_record()
    print(
        record["status"],
        "cpu",
        record["cpu_seconds"],
        "total_certified",
        record["total_certified"],
        "nf256_strong",
        record["nf256"]["source_rho_assessment"]["full_strong"]["owner_max"],
        "nf512_strong",
        record["nf512"]["source_rho_assessment"]["full_strong"]["owner_max"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
