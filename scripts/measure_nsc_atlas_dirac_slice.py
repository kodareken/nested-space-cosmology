"""Record a quiet source-anchored joined slice; no quantum slip is inferred."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from recursive_horizons.nsc_atlas_dirac_slice import reference_slice_report

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('src/recursive_horizons/nsc_atlas_dirac_slice.py',
    'src/recursive_horizons/nsc_atlas_connection_slice.py',
    'src/recursive_horizons/nsc_atlas_moving_boundaries.py',
    'tests/test_nsc_atlas_dirac_slice.py', 'scripts/measure_nsc_atlas_dirac_slice.py')


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
            raise ValueError('freeze the joined-slice method before evaluation: '+name)
        bindings[name] = hashlib.sha256(raw).hexdigest()
    report = reference_slice_report()
    if max(report['checks'].values()) > 1e-10:
        raise RuntimeError('joined-slice geometry check failed')
    if report['counts']['geometric_points'] > 40000:
        raise RuntimeError('geometric evaluation budget exhausted')
    result = dict(schema='nsc-atlas-dirac-source-slice-v1',
        status='joined_quiet_initial_surface_and_intrinsic_operator_only',
        source_commit=commit, source_sha256=bindings, reference=report,
        runtime_cpu_seconds=time.process_time()-started,
        limits=dict(CPU_seconds=60, geometric_points=40000, summary_bytes=10000))
    if result['runtime_cpu_seconds'] > 60:
        raise RuntimeError('joined-slice CPU budget exhausted')
    raw = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'
    if len(raw.encode()) > 10000:
        raise RuntimeError('joined-slice summary exceeds byte budget')
    args.output.write_text(raw, encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=result['status'], geometry=report['geometry'],
        intrinsic_spatial_Dirac=report['intrinsic_spatial_Dirac'],
        maximum_check_error=max(report['checks'].values()),
        CPU_seconds=result['runtime_cpu_seconds'])))


if __name__ == '__main__':
    main()
