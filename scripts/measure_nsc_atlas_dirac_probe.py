"""Evaluate geometry-fixed Dirac inputs; a nested gravitational slip is not supplied."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

from recursive_horizons.nsc_atlas_dirac_probe import reference_probe_report

ROOT=Path(__file__).resolve().parents[1]
SOURCES=('src/recursive_horizons/nsc_atlas_dirac_probe.py',
    'src/recursive_horizons/nsc_atlas_cap_wall_response.py',
    'src/recursive_horizons/nsc_atlas_causal_reference.py',
    'tests/test_nsc_atlas_dirac_probe.py','scripts/measure_nsc_atlas_dirac_probe.py')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    if not args.output.parent.is_dir(): raise ValueError('output parent must exist')
    started=time.process_time()
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bindings={}
    for name in SOURCES:
        raw=(ROOT/name).read_bytes()
        if subprocess.check_output(['git','show',f'{commit}:{name}'],cwd=ROOT)!=raw:
            raise ValueError('freeze the reference Dirac method before evaluation: '+name)
        bindings[name]=hashlib.sha256(raw).hexdigest()
    result=dict(schema='nsc-atlas-dirac-probe-reference-v1',
        status='geometry_inputs_computed_quantum_slip_underdetermined',
        source_commit=commit,source_sha256=bindings,reference=reference_probe_report(),
        earlier_flat_slip_is_nested_prediction=False,
        counts=dict(algebraic_reference_points=1,loop_polarizations=0,evolutions=0,new_grids=0),
        runtime_cpu_seconds=time.process_time()-started,
        limits=dict(CPU_seconds=60,summary_bytes=8000))
    if result['runtime_cpu_seconds']>60: raise RuntimeError('reference CPU budget exhausted')
    raw=json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n'
    if len(raw.encode())>8000: raise RuntimeError('reference summary exceeds byte budget')
    args.output.write_text(raw,encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=result['status'],reference=result['reference'],CPU_seconds=result['runtime_cpu_seconds'])))


if __name__=='__main__': main()
