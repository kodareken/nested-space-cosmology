"""Run the frozen shared-scalar action control once and retain a compact summary."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_collar_field import reference_report

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('src/recursive_horizons/nsc_atlas_collar_field.py',
    'src/recursive_horizons/nsc_atlas_cap_wall_response.py',
    'src/recursive_horizons/nsc_atlas_causal_reference.py',
    'src/recursive_horizons/nsc_atlas_conformal_action.py',
    'tests/test_nsc_atlas_collar_field.py', 'scripts/measure_nsc_atlas_collar_field.py')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    if not args.output.parent.is_dir():
        raise ValueError('output parent must exist')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    bindings = {}
    for name in SOURCES:
        raw = (ROOT/name).read_bytes()
        if subprocess.check_output(['git', 'show', f'{commit}:{name}'], cwd=ROOT) != raw:
            raise ValueError('freeze the shared-action method before evaluation: '+name)
        bindings[name] = hashlib.sha256(raw).hexdigest()
    started = time.process_time()
    failure, reference = None, None
    try:
        with threadpool_limits(limits=1):
            reference = reference_report()
        if reference['status'] != 'shared_action_manufactured_control_passed':
            failure = 'shared-action field/energy control exceeded 1e-6'
    except (RuntimeError, ValueError) as error:
        failure = str(error)
    elapsed = time.process_time()-started
    if elapsed > 180:
        failure = 'shared-action reference process-CPU budget exhausted'
    status = 'failed' if failure else reference['status']
    result = dict(schema='nsc-atlas-collar-field-v1', status=status,
        source_commit=commit, source_sha256=bindings, reference=reference,
        failure=failure, runtime_cpu_seconds=elapsed,
        limits=dict(CPU_seconds=180, field_rhs=6000, background_rhs=400000, summary_bytes=16000))
    raw = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'
    if len(raw.encode()) > 16000:
        raise RuntimeError('shared-action summary exceeds byte budget')
    args.output.write_text(raw, encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=status, CPU_seconds=elapsed,
        checks=None if reference is None else reference['normalized_errors'], failure=failure)))
    if failure:
        raise RuntimeError(failure)


if __name__ == '__main__':
    main()
