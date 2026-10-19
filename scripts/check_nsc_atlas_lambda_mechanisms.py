"""Recompute the strongest presently specified nesting/Lambda comparison.

Print a compact summary only; full arrays and result files are unnecessary.
"""
from __future__ import annotations

import json

from recursive_horizons.nsc_atlas_gravitational_reduction import paired_nesting_vacuum_test


if __name__ == '__main__':
    print(json.dumps(paired_nesting_vacuum_test(), indent=2, allow_nan=False))
