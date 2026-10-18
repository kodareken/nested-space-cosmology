"""Preserve the fixed quiet time-chart reference, including a coordinate-domain exit."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from recursive_horizons.nsc_atlas_dirac_time_chart import reference_chart_report

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('src/recursive_horizons/nsc_atlas_dirac_time_chart.py',
    'src/recursive_horizons/nsc_atlas_dirac_slice.py',
    'src/recursive_horizons/nsc_atlas_dirac_evolution.py',
    'src/recursive_horizons/nsc_atlas_connection_slice.py',
    'src/recursive_horizons/nsc_atlas_moving_boundaries.py',
    'src/recursive_horizons/nsc_atlas_cap_wall_response.py',
    'src/recursive_horizons/nsc_atlas_causal_reference.py',
    'tests/test_nsc_atlas_dirac_time_chart.py', 'scripts/measure_nsc_atlas_dirac_time_chart.py')


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
            raise ValueError('freeze the coordinate-chart method before evaluation: '+name)
        bindings[name] = hashlib.sha256(raw).hexdigest()
    failure = None
    try:
        reference = reference_chart_report()
        if reference['counts']['quadrature_points'] > 200000:
            raise RuntimeError('fixed chart quadrature budget exhausted')
        if max(reference['checks'].values()) > 1e-8:
            raise RuntimeError('fixed chart identity check failed')
        status = ('coordinate_domain_exit_recorded' if reference['domain']['status'] == 'exited'
                  else 'fixed_chart_samples_checked')
    except (RuntimeError, ValueError) as error:
        failure = str(error)
        reference = None
        status = 'chart_identity_or_budget_failure'
    elapsed = time.process_time()-started
    if elapsed > 120:
        failure = 'fixed chart process-CPU budget exhausted'
        status = 'chart_identity_or_budget_failure'
    result = dict(schema='nsc-atlas-dirac-time-chart-v1', status=status,
        source_commit=commit, source_sha256=bindings, reference=reference,
        failure=failure, runtime_cpu_seconds=elapsed,
        limits=dict(CPU_seconds=120, quadrature_points=200000, summary_bytes=12000))
    raw = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'
    if len(raw.encode()) > 12000:
        raise RuntimeError('fixed chart summary exceeds byte budget')
    args.output.write_text(raw, encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=status, domain=None if reference is None else reference['domain'],
        checks=None if reference is None else reference['checks'], CPU_seconds=elapsed)))
    if failure is not None:
        raise RuntimeError(failure)


if __name__ == '__main__':
    main()
