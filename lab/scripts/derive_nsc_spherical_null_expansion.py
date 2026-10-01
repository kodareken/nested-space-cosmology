#!/usr/bin/env python3
"""Write the spherical null-expansion record from the saved T=0.05 episode.

Does not evolve, recompute inheritance, or rewrite the episode JSON or NPZ.
"""
from recursive_horizons.nsc_spherical_null_expansion import build_record, write_record


def main():
    record = build_record()
    path = write_record(record)
    print(
        record["status"],
        "cpu",
        record["cpu_seconds"],
        "bytes",
        path.stat().st_size,
        "margin",
        record["margin"]["absolute_expansion_margin"],
        "shared_disagreements",
        record["comparisons"]["shared_node_sign_disagreements_max"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
