#!/usr/bin/env python3
"""Derive the initial source-selected scale and measure the withheld family.

The run prints the analytic comparison and the held-out prediction. It does
not call the initial-radius solver and it does not write a production record.
"""
from __future__ import annotations

import argparse
import json

from recursive_horizons.nsc_discovery_scale import public_summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solve", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--output", default=None, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.solve or args.output:
        raise SystemExit("this derivation does not solve the initial radius or write a record")
    summary = public_summary()
    held = summary["heldout"]["members"][0]
    print(json.dumps({
        "schema": summary["schema"],
        "equation": summary["equation"],
        "uniform_formula": summary["uniform_formula"],
        "conditional_bound": summary["conditional_bound"],
        "selection": summary["selection"],
        "dynamical_selection": summary["dynamical_selection"],
        "initial_solver_called": summary["initial_solver_called"],
        "production_record_written": summary["production_record_written"],
        "one_percent_pass_declared": summary["one_percent_pass_declared"],
        "uniform_residual_max": summary["analytic"]["uniform_residual_max"],
        "uniform_r_flat2": summary["analytic"]["uniform_r_flat2"],
        "variable_identity_gap": summary["analytic"]["variable_owner_comparison"]["identity_gap"],
        "integrated_identity_gap": summary["analytic"]["integrated_identity_gap"],
        "max_principle_holds": summary["analytic"]["max_principle_holds"],
        "algebraic_owner_gap": summary["projection"]["algebraic_owner_gap"],
        "projected_pull_gap": summary["projection"]["projected_pull_gap"],
        "reviewer_unit_coefficient_gap": summary["projection"]["reviewer_unit_coefficient_gap"],
        "areal_radius_weight": summary["coordinate"]["areal_radius_weight"],
        "held_out_carrier_nu": summary["heldout"]["held_out_carrier_nu"],
        "held_out_radius": held["prediction"]["radius"],
        "held_out_clock_rate": held["prediction"]["clock_rate"],
        "held_out_proper_ratio": held["prediction"]["proper_ratio"],
        "held_out_conditional_relative_gap": held["prediction"]["conditional_relative_gap"],
        "held_out_owner_identity_gap": held["owner_comparison"]["identity_gap"],
        "held_out_owner_relative_identity_gap": held["owner_comparison"]["relative_identity_gap"],
        "held_out_owner_residual_max": held["owner_comparison"]["owner_residual_max"],
        "held_out_source_positivity": held["source_positivity"],
        "phase_rho_integral_gap": summary["heldout"]["phase_rho_integral_gap"],
        "kinetic_sum_gap": summary["heldout"]["kinetic_sum_gap"],
    }, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
