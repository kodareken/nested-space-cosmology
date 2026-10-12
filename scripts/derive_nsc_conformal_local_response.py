#!/usr/bin/env python3
"""Read the finalized conformal episode; never integrate geometry."""
import argparse
import os

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")

from recursive_horizons.nsc_conformal_local_response import probe

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget-s", type=float, default=60.0)
    args = parser.parse_args()
    result = probe(args.budget_s)
    print(result["status"])
    print("occupation change", result["occupation_change"], "reduction error", result["occupation_error"]["max_abs"])
    print("probe omissions", {name: item["movement"] for name, item in result["controls"].items()})
    print("CPU", result["cpu_seconds"], "payload", result["payload_bytes"])
