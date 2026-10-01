#!/usr/bin/env python3
"""Replay or write the weak initial-residual record.

``--check`` reads the sealed record, remeasures scientific fields, and does
not write. A fresh record is written only to an explicit ``--output`` path,
which cannot be the sealed JSON. The v5 JSON and NPZ are read and not
rewritten.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from recursive_horizons.nsc_spherical_cauchy_weak import (
    build_record,
    verify_saved,
    write_record,
)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    if args.output is not None and args.check:
        parser.error("use either --check or --output")
    if args.output is not None:
        record = build_record()
        path = write_record(record, args.output)
        print(path)
        print(
            record["status"],
            "cpu",
            record["cpu_seconds"],
            "checkpoint_head",
            record["checkpoint_head"],
            "total_certified",
            record["total_certified"],
        )
        return 0
    report = verify_saved()
    print(
        report["status"],
        "cpu",
        report["cpu_seconds"],
        "historical_head",
        report["historical_checkpoint_head"],
        "current_head",
        report["current_checkpoint_head"],
        "head_equality_required",
        report["head_equality_required"],
        "v5_bytes_unchanged",
        report["v5_bytes_unchanged"],
        "record_bytes_unchanged",
        report["record_bytes_unchanged"],
        "generator_limits",
        len(report["limits"]),
        "total_certified",
        report["sealed"]["total_certified"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
