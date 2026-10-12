"""Measure what a grid of N points costs in memory on this machine.

The recorded family stopped silently at 3921 points on the largest grids,
with memory scaling as N squared, and the natural inference was a resource
ceiling. That inference was never tested. The machine here reports far more
memory than the recorded family's per-case figure of 205 MB at 869 points
would suggest for the depths now wanted.

This builds a real grid, packet and state at a sequence of point counts and
reports resident memory, so the next step is chosen from a measurement rather
than from an inference about a machine that may no longer be the one
involved.

Nothing is marched. This only builds and discards.
"""
from __future__ import annotations

import ctypes
import gc
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, List

LAB = Path(__file__).resolve().parents[1]
for _path in (str(LAB / "src"), str(LAB / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

import numpy as np  # noqa: E402

import derive_nsc_rn_parent as driver  # noqa: E402
from recursive_horizons import nsc_rn_pg as pg  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402

COUNTS = [869, 2000, 4000, 6000, 8000]
R_OUT = 32.0
EXCISION = 0.7555


class _MemoryStatus(ctypes.Structure):
    _fields_ = [
        ("dwLength", ctypes.c_ulong),
        ("dwMemoryLoad", ctypes.c_ulong),
        ("ullTotalPhys", ctypes.c_ulonglong),
        ("ullAvailPhys", ctypes.c_ulonglong),
        ("ullTotalPageFile", ctypes.c_ulonglong),
        ("ullAvailPageFile", ctypes.c_ulonglong),
        ("ullTotalVirtual", ctypes.c_ulonglong),
        ("ullAvailVirtual", ctypes.c_ulonglong),
        ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
    ]


def _status() -> _MemoryStatus:
    status = _MemoryStatus()
    status.dwLength = ctypes.sizeof(_MemoryStatus)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
    return status


def available_bytes() -> int:
    return int(_status().ullAvailPhys)


def committed_bytes() -> int:
    """Process commit charge: total page file less what is still available."""
    status = _status()
    return int(status.ullTotalPageFile - status.ullAvailPageFile)


def _measure(points: int) -> Dict[str, Any]:
    gc.collect()
    before = committed_bytes()
    started = time.perf_counter()
    geometry = driver.case_geometry(1.04)
    grid = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"], points=points,
        r_out=R_OUT, excision=EXCISION,
    )
    packet = source.prepare_radial_packet(
        mass_over_rm=1.04, energy_target=1.0e-3, rank=driver.PILOT_RANK,
        model=grid["model"], radii=grid["radius"], weights=grid["weights"],
        derivative=grid["derivative"], center=3.0 * geometry["r_m"],
    )
    state = pg.make_state(
        grid, np.asarray(packet["phi"]), inner_mass=geometry["mass"],
        killing_frequency=float(packet["radial_moments"]["requested_omega"]),
        occupations=np.asarray(packet["occupations"], dtype=float),
    )
    rates = pg.stage_rates(state, grid)
    elapsed = time.perf_counter() - started
    after = committed_bytes()
    spacing = float(grid["spacing"])
    derivative = grid["derivative"]
    nnz = int(derivative.nnz) if hasattr(derivative, "nnz") else int(derivative.size)
    result = {
        "points": points,
        "r_inner": float(grid["radius"][0]),
        "spacing": spacing,
        "gap_in_cells": (grid["radius"][0] - grid["r_minus"]) / spacing,
        "commit_delta_bytes": int(after - before),
        "commit_total_bytes": int(after),
        "operator_nnz": nnz,
        "operator_bytes": nnz * 8,
        "mass_finite": bool(np.all(np.isfinite(rates["mass"]))),
        "build_seconds": elapsed,
        "available_after_bytes": available_bytes(),
    }
    del state, rates, packet, grid
    gc.collect()
    return result


def main() -> int:
    status = _status()
    print(f"machine: {status.ullTotalPhys / 2**30:.1f} GB total, "
          f"{status.ullAvailPhys / 2**30:.1f} GB available, "
          f"{os.cpu_count()} logical cores")
    print()
    print("%8s %12s %12s %14s %14s %9s" %
          ("points", "spacing", "gap cells", "commit delta", "commit total", "build s"))
    rows: List[Dict[str, Any]] = []
    for points in COUNTS:
        if available_bytes() < 8 * 2 ** 30 and points > 4000:
            print("%8d skipped: only %.1f GB available"
                  % (points, available_bytes() / 2 ** 30))
            continue
        try:
            row = _measure(points)
        except MemoryError:
            print("%8d MemoryError" % points)
            gc.collect()
            continue
        rows.append(row)
        print("%8d %12.6f %12.2f %12.1f GB %12.1f GB %9.1f"
              % (row["points"], row["spacing"], row["gap_in_cells"],
                 row["commit_delta_bytes"] / 2 ** 30,
                 row["commit_total_bytes"] / 2 ** 30, row["build_seconds"]))

    out = {
        "schema": "NSC-RN-GRID-MEMORY-v1",
        "what_this_measures": (
            "What a grid of N points costs on this machine, at a fixed inner edge of "
            "0.7555 and outer radius 32. Nothing is marched. The point is to choose the "
            "next step from a measurement rather than from an inference about a resource "
            "ceiling that was never tested."
        ),
        "machine": {
            "total_gb": status.ullTotalPhys / 2 ** 30,
            "available_gb": status.ullAvailPhys / 2 ** 30,
            "logical_cores": os.cpu_count(),
            "measurement": "GlobalMemoryStatusEx; psutil is not installed here",
        },
        "r_out": R_OUT,
        "excision": EXCISION,
        "rows": rows,
    }
    path = LAB / "results/development/nsc-rn-grid-memory-v1.json"
    path.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print()
    print("wrote %s" % path.relative_to(LAB))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())