#!/usr/bin/env python3
"""Run the full-interval L0 A4 attempt. Do not publish a null remainder."""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from recursive_horizons.nsc_ks_full_interval_uv import (  # noqa: E402
    full_interval_uv_report,
    validate_full_interval_uv_report,
)


def calculate(root=None):
    report = full_interval_uv_report(ROOT if root is None else root)
    validate_full_interval_uv_report(report)
    return report


def publish(report):
    if report.get('numerical_C_M') is None or report.get('induced_N_beta_remainder') is None:
        raise SystemExit(
            'no finite C_M or induced N/beta remainder; exclusive record not published')
    raise SystemExit('publication of a finite remainder is not implemented in this attempt')


def _summary(report):
    return {
        'schema': report['schema'],
        'history': report['history_identity'][:16],
        'full_interval_covered': report['full_interval_covered'],
        'domain_failure': report['domain']['domain_failure'],
        'major_A2_contains_zero': report['major_A2_contains_zero'],
        'reentry_upper': report['reentry']['reentry_majorant']['binary64_upper'],
        'pointwise_reentry_contains_zero': report['pointwise_counterexample'][
            'pointwise_reentry_contains_zero'],
        'absolute_common_carrier_integral_upper': report['carrier'][
            'absolute_common_carrier_integral_upper'],
        'numerical_C_M': report['numerical_C_M'],
        'induced_N_beta_remainder': report['induced_N_beta_remainder'],
        'unconstructed_term': report['unconstructed_term'],
        'cpu_seconds': report['cpu_seconds'],
        'rss_mib': report['rss_mib'],
        'physical_local_gate': report['physical_local_gate'],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--record', action='store_true')
    args = parser.parse_args(argv)
    report = calculate()
    print(json.dumps(_summary(report), indent=2, sort_keys=True))
    if args.record:
        publish(report)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
