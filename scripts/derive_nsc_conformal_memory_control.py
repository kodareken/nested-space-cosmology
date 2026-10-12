#!/usr/bin/env python3
"""Bounded independent reference for the existing conformal memory omission."""
import os

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

from recursive_horizons.nsc_conformal_memory_control import execute

if __name__ == "__main__":
    record = execute()
    print(record["status"], "CPU", record["cpu_seconds"])
    print("memory-off error", record["streamed_memory_off_error_against_independent_reference"],
          "fraction", record["error_fraction_of_memory_effect"],
          "midpoint indicator", record["reference_substeps_2_to_4_indicator"])
