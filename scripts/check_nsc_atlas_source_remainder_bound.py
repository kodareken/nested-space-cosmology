#!/usr/bin/env python3
"""One exact rational source remainder application, without reintegration."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

from check_nsc_atlas_onset_endpoint_bound import plain,read_bound
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_source_remainder_bound import source_remainder_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT=Path(__file__).resolve().parents[1]
INPUTS={'results/summaries/nsc-atlas-source-forcing-bound-v1.json':
        '2ccf4ad7326291b79500f4d78a0f858907a80901c74b771f8b22fd8cb735c611'}
OWNERS=('scripts/check_nsc_atlas_source_remainder_bound.py',
 'src/recursive_horizons/nsc_atlas_source_remainder_bound.py',
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
    source=read_bound(*next(iter(INPUTS.items())))
    row=source['input'];original=read_bound(row['path'],row['sha256'])
    authenticated={**INPUTS,row['path']:row['sha256']}
    if (source['status']!='certified_actual_source_forcing'
            or original['status']!='certified_source_representation_cell'):
        raise ValueError('require original certified source and representation cell')
    for key in ('left','right'):
        if Q(source['parameters'][key])!=Q(original['cell'][key]):
            raise ValueError('source/representation cell mismatch')
    budget=Budget(maximum_evaluations=8192,cpu_seconds=5.)
    report=dict(schema='nsc-atlas-source-remainder-bound-v1',status='incomplete',
        source_commit=commit,source_bindings=bindings,
        inputs=[dict(path=n,sha256=h) for n,h in authenticated.items()],
        runtime=dict(python=sys.version.split()[0],platform=sys.platform),
        limits=dict(charged_values=8192,cpu_seconds=5.,arithmetic='exact rational',applications=1),
        cell={k:source['parameters'][k] for k in ('left','right','alpha_max')})
    try:
        result=source_remainder_bound(source['parameters'],source['certificate'],budget=budget)
        outward=lambda v:0. if v==0 else math.nextafter(float(v),math.inf)
        report.update(status='certified_quadratic_source_remainder',certificate=plain(result),
            display={section:{k:outward(v) for k,v in result[section].items()} for section in
                ('coefficients','maximum_amplitude_remainders')})
        budget.check()
    except (ValueError,RuntimeError,KeyError,TypeError,OverflowError) as error:
        report.update(status='stopped',failure=str(error))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_queries=0,source_integrations=0,field_evolutions=0,time_steps=0,root_searches=0,
        retries=0,raw_history_saved=False)
    report['regenerate']=('python scripts/run_python.py scripts/check_nsc_atlas_source_remainder_bound.py'
                          ' --output <new-output-path>')
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report.get(k) for k in ('status','display','execution','failure')}))
    return 0 if report['status']=='certified_quadratic_source_remainder' else 1


if __name__=='__main__':
    raise SystemExit(main())
