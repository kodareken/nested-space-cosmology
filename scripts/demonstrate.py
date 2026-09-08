#!/usr/bin/env python3
"""Recompute three compact demonstrations and state their physical scope."""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

from reproduce_public_results import load_manifest, validate_checkout

ROOT = Path(__file__).resolve().parents[1]
CASES = (
    ('scripts/check_nsc_boundary_response.py',
     'Finite radial boundary maps: independently computed spatial response; physical Lorentzian matching remains open.'),
    ('scripts/check_nsc_vacuum_work.py',
     'Vacuum response: prescribed geometric work produces finite Dirac excitations; the geometry work source remains unsolved.'),
    ('scripts/check_nsc_observable_bridge.py',
     'Spinor bridge: exact identification requirements and a kinematic scalar link; particle and cosmology closure remains open.'),
)


def main() -> int:
    manifest = load_manifest()
    validate_checkout(manifest)
    env = os.environ.copy()
    env.update({key: '1' for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS',
                                    'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')})
    env['PYTHONHASHSEED'] = '0'
    for generator, scope in CASES:
        subprocess.run([sys.executable, str(ROOT / generator), '--check'],
                       cwd=ROOT, env=env, check=True)
        print(scope, flush=True)
    print(f'Three scoped demonstrations passed within the {manifest["result_count"]}-record release.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
