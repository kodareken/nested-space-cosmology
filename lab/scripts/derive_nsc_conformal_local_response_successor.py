#!/usr/bin/env python3
"""Production local reduction on the finalized conformal .2->.3 segment."""
import argparse
import os

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"

from recursive_horizons.nsc_conformal_local_response_successor import execute, verify

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-budget-s", type=float, default=300.)
    parser.add_argument("--confirmation-budget-s", type=float, default=1800.)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    record = verify() if args.check else execute(args.pilot_budget_s, args.confirmation_budget_s)
    print(record["status"], "CPU", record.get("cpu_seconds"))
    print("occupation change", record["occupation_effect"], "error", record["occupation_error"]["max_abs"])
    print("controls", record["controls"])
    print("all numerical movements within one percent", record["all_numerical_movements_within_one_percent"])
