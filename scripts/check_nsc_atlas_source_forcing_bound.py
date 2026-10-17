#!/usr/bin/env python3
"""Convert existing pulse enclosures into actual source-forcing bounds."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_source_forcing_bound import source_forcing_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT=Path(__file__).resolve().parents[1]
INPUT='results/summaries/nsc-atlas-source-cell-bound-v2.json'
INPUT_SHA='707bb8d76359edaf2a26228657a219933ebd164edcf89eb87742f9e66bec2bdd'
OWNERS=('scripts/check_nsc_atlas_source_forcing_bound.py',
        'src/recursive_horizons/nsc_atlas_source_forcing_bound.py',
        'src/recursive_horizons/nsc_atlas_source_cell_bound.py',
        'src/recursive_horizons/nsc_atlas_characteristic_initialization.py',
        'src/recursive_horizons/nsc_source_binding_graph.py','scripts/run_python.py','pyproject.toml')


def plain(value):
    if isinstance(value,Q):return str(value)
    if isinstance(value,dict):return {k:plain(v) for k,v in value.items()}
    return value


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('require an existing parent and a new output path')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bindings=original_source_graph(ROOT,OWNERS,commit)
    for name,expected in bindings.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('freeze changed producer before execution: '+name)
    raw=(ROOT/INPUT).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=INPUT_SHA:raise ValueError('original source certificate changed')
    source=json.loads(raw);cert=source['certificate']
    if (source['status']!='certified_source_representation_cell'
            or not cert['tube']['exact_continuation_proved']
            or any(Q(v)!=0 for v in source['cell']['normalized_exact_binary_powers'][:2])):
        raise ValueError('source does not provide the original quiet-initial interval')
    inputs={k:Q(source['cell'][v]) for k,v in (('left','left'),('right','right'),('alpha_max','alpha'))}
    for k in ('g','g_x','g_xx'):
        row=cert['analytic_source_bounds'][k]['upper_dyadic']
        inputs[k]=Q(int(row['mantissa']))*Q(2)**row['exponent']
    budget=Budget(maximum_evaluations=4096,cpu_seconds=5.)
    report=dict(schema='nsc-atlas-source-forcing-bound-v1',source_commit=commit,source_bindings=bindings,
                input=dict(path=INPUT,sha256=INPUT_SHA),parameters=plain(inputs),status='incomplete')
    try:
        certificate=source_forcing_bound(**inputs,budget=budget)
        report.update(status='certified_actual_source_forcing',certificate=plain(certificate),
                      display={k:(0. if v==0 else math.nextafter(float(v),math.inf)) for k,v in certificate['forcing'].items()})
    except (ValueError,RuntimeError,KeyError,TypeError,OverflowError) as error:
        report.update(status='stopped',failure=str(error))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_queries=0,source_integrations=0,time_steps=0,root_searches=0,retries=0,raw_history_saved=False)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report.get(k) for k in ('status','display','execution','failure')}))
    return 0 if report['status']=='certified_actual_source_forcing' else 1


if __name__=='__main__':raise SystemExit(main())
