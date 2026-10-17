#!/usr/bin/env python3
"""One radial cap comparison from authenticated onset certificates."""
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
from recursive_horizons.nsc_atlas_cap_departure_bound import cap_departure_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT=Path(__file__).resolve().parents[1]
INPUTS={
 'results/summaries/nsc-atlas-coupled-tube-bound-v1.json':'192caf213cd4b67961f815203687a87748a9668c4a97d2d77f9ee2f5d78a8dab',
 'results/summaries/nsc-atlas-scalar-departure-bound-v1.json':'fb52782856689b88ce2d6bee893085352eda9c380372ec178d2aa670196ef7f5'}
OWNERS=('scripts/check_nsc_atlas_cap_departure_bound.py',
 'src/recursive_horizons/nsc_atlas_cap_departure_bound.py',
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
    tube,scalar=[read_bound(name,digest) for name,digest in INPUTS.items()]
    authenticated=dict(INPUTS)
    for row in [*tube['inputs'].values(),*scalar['inputs']]:
        read_bound(row['path'],row['sha256']);authenticated[row['path']]=row['sha256']
    if (not tube['certificate']['strict'] or tube['uniform_in_coupling'] is not True
            or scalar['status']!='certified_amplitude_scalar_departure'):
        raise ValueError('require original strict tube and scalar comparison')
    if Q(scalar['cell']['right'])-Q(scalar['cell']['left'])!=Q(tube['certificate']['duration']):
        raise ValueError('scalar/cap original domain mismatch')
    budget=Budget(maximum_evaluations=12288,cpu_seconds=10.)
    report=dict(schema='nsc-atlas-cap-departure-bound-v1',status='incomplete',
        source_commit=commit,source_bindings=bindings,
        inputs=[dict(path=n,sha256=h) for n,h in authenticated.items()],
        runtime=dict(python=sys.version.split()[0],python_flint=version('python-flint'),platform=sys.platform),
        limits=dict(charged_values=12288,cpu_seconds=10.,precision_bits=160,applications=1),
        cell=scalar['cell'])
    try:
        result=cap_departure_bound(tube['certificate'],scalar['certificate'],budget=budget)
        outward=lambda v:0. if v==0 else math.nextafter(float(v),math.inf)
        report.update(status='certified_amplitude_cap_departure',certificate=plain(result),
            display={section:{k:outward(v) for k,v in result[section].items()} for section in
                ('coefficients','maximum_amplitude_differences','maximum_amplitude_second_order')})
        budget.check()
    except (ValueError,RuntimeError,KeyError,TypeError,OverflowError) as error:
        report.update(status='stopped',failure=str(error))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_queries=0,source_integrations=0,field_evolutions=0,time_steps=0,root_searches=0,
        retries=0,raw_history_saved=False)
    report['regenerate']=('python scripts/run_python.py scripts/check_nsc_atlas_cap_departure_bound.py'
                          ' --output <new-output-path>')
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report.get(k) for k in ('status','display','execution','failure')}))
    return 0 if report['status']=='certified_amplitude_cap_departure' else 1


if __name__=='__main__':
    raise SystemExit(main())
