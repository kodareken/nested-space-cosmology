"""Record conformal Dirac matching and actual one-sided radial curvature."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from recursive_horizons.nsc_atlas_dirac_conformal import reference_conformal_report
from recursive_horizons.nsc_atlas_dirac_trace import reference_bulk_trace_report

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ('src/recursive_horizons/nsc_atlas_dirac_conformal.py',
    'src/recursive_horizons/nsc_atlas_dirac_probe.py',
    'src/recursive_horizons/nsc_atlas_cap_wall_response.py',
    'src/recursive_horizons/nsc_atlas_causal_reference.py',
    'src/recursive_horizons/nsc_atlas_dirac_trace.py',
    'tests/test_nsc_atlas_dirac_trace.py',
    'tests/test_nsc_atlas_dirac_conformal.py', 'scripts/measure_nsc_atlas_dirac_conformal.py')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--bulk-trace', action='store_true',
        help='evaluate the local one-loop bulk trace in the original physical metric')
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
            raise ValueError('freeze the conformal Dirac method before evaluation: '+name)
        bindings[name] = hashlib.sha256(raw).hexdigest()
    reference = reference_bulk_trace_report() if args.bulk_trace else reference_conformal_report()
    if args.bulk_trace:
        errors = list(reference['checks'].values())
    else:
        interface = reference['conformal_interface']
        errors = [abs(interface[key]) for key in ('conformal_trace_K_jump',
            'conformal_temporal_K_jump', 'conformal_angular_K_jump',
            'conformal_normal_spinor_derivative_jump')]
        errors.append(abs(reference['radial_curvature_jump']-reference['radial_curvature_jump_control']))
    if max(errors) > 1e-12:
        raise RuntimeError('conformal material-interface check failed')
    schema = 'nsc-atlas-dirac-bulk-trace-v1' if args.bulk_trace else 'nsc-atlas-dirac-conformal-interface-v1'
    status = ('local_one_loop_bulk_trace_only' if args.bulk_trace
              else 'classical_interface_jets_cancel_bulk_curvature_change_remains')
    result = dict(schema=schema, status=status,
        source_commit=commit, source_sha256=bindings, reference=reference,
        runtime_cpu_seconds=time.process_time()-started,
        limits=dict(CPU_seconds=30, summary_bytes=8000))
    if result['runtime_cpu_seconds'] > 30:
        raise RuntimeError('conformal Dirac CPU budget exhausted')
    raw = json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+'\n'
    if len(raw.encode()) > 8000:
        raise RuntimeError('conformal Dirac summary exceeds byte budget')
    args.output.write_text(raw, encoding='utf-8', newline='\n')
    print(json.dumps(dict(status=result['status'],
        result=reference['one_sided'] if args.bulk_trace else reference['conformal_interface'],
        maximum_identity_error=max(errors), CPU_seconds=result['runtime_cpu_seconds'])))


if __name__ == '__main__':
    main()
