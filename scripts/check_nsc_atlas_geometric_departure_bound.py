#!/usr/bin/env python3
"""Consume the original uniform tube once; no physical field evaluation."""
import argparse
from fractions import Fraction as Q
import hashlib
from importlib.metadata import version
import json
import math
from pathlib import Path
import subprocess
import sys
import time

from check_nsc_atlas_onset_endpoint_bound import plain,read_bound
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_geometric_departure_bound import geometric_departure_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT=Path(__file__).resolve().parents[1]
INPUTS={
 'results/summaries/nsc-atlas-source-forcing-bound-v1.json':'2ccf4ad7326291b79500f4d78a0f858907a80901c74b771f8b22fd8cb735c611',
 'results/summaries/nsc-atlas-reference-gauge-bound-v1.json':'e8a8400cec0f6ae0bd831723696ebf1d41c38bdf0b052301ffab837465f1d216',
 'results/summaries/nsc-atlas-coupled-tube-bound-v1.json':'192caf213cd4b67961f815203687a87748a9668c4a97d2d77f9ee2f5d78a8dab'}
OWNERS=('scripts/check_nsc_atlas_geometric_departure_bound.py',
 'src/recursive_horizons/nsc_atlas_geometric_departure_bound.py',
 'src/recursive_horizons/nsc_atlas_source_forcing_bound.py',
 'src/recursive_horizons/nsc_source_binding_graph.py','scripts/run_python.py','pyproject.toml')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('require a new output path in an existing directory')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bindings=original_source_graph(ROOT,OWNERS,commit)
    for name,digest in bindings.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('commit calculation source before execution: '+name)
    source,gauge,tube=[read_bound(name,digest) for name,digest in INPUTS.items()]
    authenticated=dict(INPUTS)
    for row in [*tube['inputs'].values(),source['input']]:
        read_bound(row['path'],row['sha256']);authenticated[row['path']]=row['sha256']
    if (source['status']!='certified_actual_source_forcing'
            or gauge['status']!='certified_exact_reference_gauge_cell'
            or not tube['certificate']['strict'] or tube['uniform_in_coupling'] is not True):
        raise ValueError('require original uniform source/reference/strict-tube certificates')
    for key in ('left','right'):
        if Q(source['parameters'][key])!=Q(gauge['certificate']['inputs'][key]):
            raise ValueError('source/reference domain mismatch')
    budget=Budget(maximum_evaluations=32768,cpu_seconds=10.)
    report=dict(schema='nsc-atlas-geometric-departure-bound-v1',status='incomplete',
        source_commit=commit,source_bindings=bindings,
        inputs=[dict(path=n,sha256=h) for n,h in authenticated.items()],
        runtime=dict(python=sys.version.split()[0],python_flint=version('python-flint'),platform=sys.platform),
        limits=dict(charged_values=32768,cpu_seconds=10.,precision_bits=160,
                    affine_directions=4,matrix_dimension=3),
        cell={k:source['parameters'][k] for k in ('left','right','alpha_max')})
    try:
        result=geometric_departure_bound(source['certificate'],gauge['certificate'],tube['certificate'],budget=budget)
        report.update(status='certified_amplitude_geometric_departure',certificate=plain(result),
            display=dict(weighted_norm_upper=math.nextafter(float(result['weighted_norm']),math.inf),
                coefficients={k:0. if v==0 else math.nextafter(float(v),math.inf) for k,v in result['coefficients'].items()},
                endpoint_bounds={k:0. if v==0 else math.nextafter(float(v),math.inf) for k,v in result['endpoint_bounds'].items()}))
        budget.check()
    except (ValueError,RuntimeError,KeyError,TypeError,OverflowError) as error:
        report.update(status='stopped',failure=str(error))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_queries=0,source_integrations=0,field_evolutions=0,time_steps=0,root_searches=0,
        retries=0,raw_history_saved=False)
    report['regenerate']=('python scripts/run_python.py scripts/check_nsc_atlas_geometric_departure_bound.py'
                          ' --output <new-output-path>')
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report.get(k) for k in ('status','display','execution','failure')}))
    return 0 if report['status']=='certified_amplitude_geometric_departure' else 1


if __name__=='__main__':
    raise SystemExit(main())
