#!/usr/bin/env python3
"""One authenticated interval pass for the first onset geometric coefficients."""
import argparse
from fractions import Fraction as Q
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys
import time

from check_nsc_atlas_onset_endpoint_bound import INPUTS as SCALAR_INPUTS
from check_nsc_atlas_onset_endpoint_bound import OWNERS as SCALAR_OWNERS
from check_nsc_atlas_onset_endpoint_bound import display, plain, read_bound
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_onset_geometry_bound import onset_geometry_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT = Path(__file__).resolve().parents[1]
INPUTS = {**SCALAR_INPUTS,
    'results/summaries/nsc-atlas-onset-endpoint-bound-v1.json':
        'c03772142d1852b80c292225e39d5a6f5bf7bf1b4ab6725870a48fb4127ada41'}
OWNERS = (*SCALAR_OWNERS,
    'scripts/check_nsc_atlas_onset_geometry_bound.py',
    'src/recursive_horizons/nsc_atlas_onset_geometry_bound.py')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('require an existing parent and a new output path')
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    bindings = original_source_graph(ROOT, OWNERS, commit)
    for name, digest in bindings.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
            raise ValueError('commit calculation source before execution: '+name)
    source, gauge, endpoint = [read_bound(name, digest) for name, digest in INPUTS.items()]
    raw_input = source['input']
    original = read_bound(raw_input['path'], raw_input['sha256'])
    raw_identity = {key:raw_input[key] for key in ('path','sha256')}
    if (source['status'] != 'certified_source_representation_cell'
            or gauge['status'] != 'certified_exact_reference_gauge_cell'
            or endpoint['status'] != 'certified_leading_wall_endpoint'
            or source['cell'] != gauge['cell']):
        raise ValueError('require original source, reference and leading-field certificates')
    for name in ('left', 'right', 'radius_lower', 'radius_upper'):
        if Q(gauge['certificate']['inputs'][name]) != Q(endpoint['parameters'][name]):
            raise ValueError('endpoint/reference domain mismatch')
    if endpoint['inputs'] != [dict(path=p,sha256=h) for p,h in SCALAR_INPUTS.items()]+[raw_identity]:
        raise ValueError('leading-field input identities changed')
    budget = Budget(maximum_evaluations=200000, cpu_seconds=30.)
    report = dict(schema='nsc-atlas-onset-geometry-bound-v1',status='incomplete',
        source_commit=commit,source_bindings=bindings,
        inputs=[dict(path=p,sha256=h) for p,h in INPUTS.items()]+[raw_identity],
        runtime=dict(python=sys.version.split()[0],python_flint=version('python-flint'),platform=sys.platform),
        parameters=endpoint['parameters'],
        limits=dict(charged_values=200000,cpu_seconds=30.,precision_bits=160,interval_panels=128,
                    physical_events=1))
    try:
        result = onset_geometry_bound(gauge['certificate'],original['geometry'],endpoint['certificate'],
                                      budget=budget,panels=128)
        report.update(status='certified_first_geometric_endpoint',certificate=plain(result),
                      display={key:display(value) for key,value in result['outputs'].items()})
        budget.check()
    except (ValueError,RuntimeError,KeyError,TypeError,OverflowError) as error:
        report.update(status='stopped',failure=str(error))
    report['execution'] = dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        interval_panels=128 if report['status']=='certified_first_geometric_endpoint' else None,
        source_ode_solves=0,field_evolutions=0,geometric_ode_solves=0,root_searches=0,
        retries=0,raw_history_saved=False)
    report['regenerate'] = ('python scripts/run_python.py scripts/check_nsc_atlas_onset_geometry_bound.py'
                            ' --output <new-output-path>')
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',
                          encoding='utf-8',newline='\n')
    print(json.dumps({key:report.get(key) for key in ('status','display','execution','failure')}))
    return 0 if report['status']=='certified_first_geometric_endpoint' else 1


if __name__ == '__main__':
    raise SystemExit(main())
