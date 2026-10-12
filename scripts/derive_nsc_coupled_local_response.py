#!/usr/bin/env python3
"""Conditional consumer of stored Galerkin geometry.

The default domain measures one streamed probe on the T=0.05 episode, then
keeps only the prefix that fits in the process-time budget. ``--successor``
reads the completed regeneration episode on ``[0.05, 0.085]`` and does not
evolve that geometry again. ``--check`` rereads the successor record and
does not write.
"""
from __future__ import annotations

import argparse
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

from recursive_horizons.nsc_coupled_local_response import execute, execute_successor, verify_successor


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--domain", choices=("pilot", "stored-frames"), default="pilot")
    parser.add_argument("--case", default="nf512_dt_0_0005")
    parser.add_argument("--budget-s", type=float, default=None)
    parser.add_argument("--confirm-budget-s", type=float, default=1800.0)
    parser.add_argument("--substeps", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--successor", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        record = verify_successor()
        gap = record["phases"]["conditional_versus_autonomous"]["full"]["max_abs"]
        drive = record["phases"]["full_controls"]["omissions"]["outside_drive"]["occupation_error_against_reference"]["max_abs"]
        cross = record["phases"]["full_controls"]["omissions"]["initial_cross"]
        print(f"status {record['status']}")
        print(f"cpu_seconds {record.get('cpu_seconds_this_process')}")
        print(f"conditional_versus_autonomous {gap}")
        print(f"drive_omission {drive}")
        print(f"cross_active {cross.get('active')} cross_occupation {cross.get('occupation_diagonal_movement')}")
        print(f"method_error_within_one_percent {record.get('method_error_within_one_percent')}")
        print("check readonly")
        return 0
    if args.successor:
        record = execute_successor(
            budget_s=600.0 if args.budget_s is None else args.budget_s,
            confirm_budget_s=args.confirm_budget_s,
            resume=args.resume,
        )
        full = record.get("phases", {}).get("full_window") or {}
        gap = (record.get("phases", {}).get("conditional_versus_autonomous") or {}).get("full") or {}
        print(f"status {record['status']}")
        print(f"full_window {full.get('t_start')} {full.get('t_end')} substeps {full.get('substeps')}")
        print(f"conditional_versus_autonomous {gap.get('max_abs')}")
        print(f"cpu_seconds {record.get('cpu_seconds_this_process')}")
        print(f"payload_bytes {record.get('payload_bytes')} within_64MiB {record.get('payload_within_64MiB')}")
        if record.get("status") == "REJECTED_REQUEST_MISMATCH":
            print("resume rejected a different successor binding")
            return 2
        return 0 if record.get("status") == "MEASURED_SUCCESSOR_RESPONSE" else 1
    record = execute(
        domain=args.domain,
        case=args.case,
        budget_s=60.0 if args.budget_s is None else args.budget_s,
        resume=args.resume,
        substeps=args.substeps,
    )
    headline = record.get("headline") or {}
    error = (headline.get("occupation_error") or {}).get("relative_to_signal")
    full = (record.get("phases") or {}).get("full_window") or {}
    full_error = ((full.get("comparison") or {}).get("occupation_error") or {}).get("relative_to_signal")
    print(f"status {record['status']}")
    print(f"pilot_window {headline.get('t_start')} {headline.get('t_end')} steps {headline.get('steps')} occupation_relative {error}")
    print(f"pilot_cpu_seconds {record.get('cpu_seconds_this_process')}")
    if full:
        print(
            f"full_window {full.get('t_start')} {full.get('t_end')} "
            f"substeps {full.get('substeps')} occupation_relative {full_error} "
            f"cpu_seconds {record.get('full_window_cpu_seconds')}"
        )
        saved = (full.get("saved_autonomous") or {}).get("conditional_full_versus_saved") or {}
        print(f"conditional_versus_saved_relative_to_change {saved.get('relative_to_change')}")
    print(f"payload_bytes {record.get('payload_bytes')} within_64MiB {record.get('payload_within_64MiB')}")
    if record.get("status") == "REJECTED_REQUEST_MISMATCH":
        print("resume rejected a different domain, substeps, source, or settings")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
