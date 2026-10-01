#!/usr/bin/env python3
"""Pilot consumer of the stored T=0.05 Galerkin geometry.

The default domain measures one streamed probe, then keeps only the prefix
that fits in the process-time budget. It does not regenerate the episode and
does not launch the feedback-episode driver. ``--domain stored-frames`` is the
full stored-frame comparison for a later, larger budget.
"""
from __future__ import annotations

import argparse
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

from recursive_horizons.nsc_coupled_local_response import execute


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--domain", choices=("pilot", "stored-frames"), default="pilot")
    parser.add_argument("--case", default="nf512_dt_0_0005")
    parser.add_argument("--budget-s", type=float, default=60.0)
    parser.add_argument("--substeps", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    record = execute(
        domain=args.domain,
        case=args.case,
        budget_s=args.budget_s,
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
