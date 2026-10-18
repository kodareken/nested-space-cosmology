"""Evaluate the initial Dirac time-generator and moving-wall compatibility."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from recursive_horizons.nsc_atlas_dirac_evolution import reference_wall_time_report

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('src/recursive_horizons/nsc_atlas_dirac_evolution.py',
    'src/recursive_horizons/nsc_atlas_dirac_slice.py',
    'src/recursive_horizons/nsc_atlas_connection_slice.py',
    'src/recursive_horizons/nsc_atlas_moving_boundaries.py',
    'tests/test_nsc_atlas_dirac_evolution.py', 'scripts/measure_nsc_atlas_dirac_wall_time.py')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if not args.output.parent.is_dir():
        raise ValueError('output parent must exist')
    started = time.process_time()
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    bindings = {}
    for name in SOURCES:
        raw = (ROOT/name).read_bytes()
        if subprocess.check_output(['git', 'show', f'{commit}:{name}'], cwd=ROOT) != raw:
            raise ValueError('freeze the wall time-generator method before evaluation: '+name)
        bindings[name] = hashlib.sha256(raw).hexdigest()
    reference = reference_wall_time_report()
    if max(reference['checks'].values()) > 1e-12:
        raise RuntimeError('Dirac moving-wall identity check failed')
    result = dict(schema='nsc-atlas-dirac-wall-time-v1',
        status='initial_time_generator_and_moving_domain_compatible',
        source_commit=commit, source_sha256=bindings, reference=reference,
        runtime_cpu_seconds=time.process_time()-started,
        limits=dict(CPU_seconds=30, summary_bytes=8000))
    if result['runtime_cpu_seconds'] > 30:
        raise RuntimeError('wall time-generator CPU budget exhausted')
    raw = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'
    if len(raw.encode()) > 8000:
        raise RuntimeError('wall time-generator summary exceeds byte budget')
    args.output.write_text(raw, encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=result['status'], wall=reference['wall'],
        checks=reference['checks'], CPU_seconds=result['runtime_cpu_seconds'])))


if __name__ == '__main__':
    main()
