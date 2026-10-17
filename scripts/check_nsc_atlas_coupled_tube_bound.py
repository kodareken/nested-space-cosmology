#!/usr/bin/env python3
"""Evaluate one predeclared conditional sourced-tube image, without evolution."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import time

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_coupled_tube_bound import DECLARED_BOX,coupled_vacuum_tube_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT=Path(__file__).resolve().parents[1]
INPUTS={
 'source':('results/summaries/nsc-atlas-source-forcing-bound-v1.json',
           '2ccf4ad7326291b79500f4d78a0f858907a80901c74b771f8b22fd8cb735c611'),
 'gauge':('results/summaries/nsc-atlas-reference-gauge-bound-v1.json',
          'e8a8400cec0f6ae0bd831723696ebf1d41c38bdf0b052301ffab837465f1d216'),
 'initial':('results/summaries/nsc-atlas-reference-initial-bound-v1.json',
            'b1d29d8a2c78518889706adbb8ab3293edda596725949347d2830c187f9ffde9')}
OWNERS=('scripts/check_nsc_atlas_coupled_tube_bound.py',
        'src/recursive_horizons/nsc_atlas_coupled_tube_bound.py',
        'src/recursive_horizons/nsc_atlas_cap_characteristics.py',
        'src/recursive_horizons/nsc_atlas_nonlinear_boundary.py',
        'src/recursive_horizons/nsc_atlas_nonlinear_bulk.py',
        'src/recursive_horizons/nsc_source_binding_graph.py','scripts/run_python.py','pyproject.toml')


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('require an existing parent and a new output path')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bindings=original_source_graph(ROOT,OWNERS,commit)
    for name,expected in bindings.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('freeze changed producer before execution: '+name)
    saved={}
    for key,(name,expected) in INPUTS.items():
        raw=(ROOT/name).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=expected:raise ValueError('input identity changed: '+name)
        saved[key]=json.loads(raw)
    source,gauge,initial=(saved[k] for k in ('source','gauge','initial'))
    if (source['status']!='certified_actual_source_forcing'
            or gauge['status']!='certified_exact_reference_gauge_cell'
            or initial['status']!='certified_quiet_exact_reference_initial_event'
            or source['parameters']['left']!=initial['certificate']['left']
            or Q(source['parameters']['left'])!=Q(gauge['cell']['left'])
            or Q(source['parameters']['right'])!=Q(gauge['cell']['right'])):
        raise ValueError('source, gauge and initial certificates do not describe the same cell')
    original=initial['inputs']['onset'];raw=(ROOT/original['path']).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=original['sha256']:
        raise ValueError('original collar input changed')
    zstar=Q(json.loads(raw)['geometry']['exact_cap_knots'][1])
    budget=Budget(maximum_evaluations=20000,cpu_seconds=10.)
    report=dict(schema='nsc-atlas-coupled-tube-bound-v1',source_commit=commit,source_bindings=bindings,
                inputs={k:dict(path=v[0],sha256=v[1]) for k,v in INPUTS.items()},
                original_collar_input=original,status='incomplete')
    try:
        certificate=coupled_vacuum_tube_bound(source['certificate'],gauge['certificate'],
                                              z_star=zstar,box=DECLARED_BOX,budget=budget)
        report.update(status=('strict_conditional_bootstrap_image' if certificate['strict']
                              else 'unclosed_conditional_bootstrap_image'),certificate=certificate,
                      existence_proved=False,uniform_in_coupling=True)
    except (ValueError,RuntimeError,KeyError,TypeError,OverflowError) as error:
        report.update(status='stopped',failure=str(error))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_queries=0,source_integrations=0,time_steps=0,root_searches=0,retries=0,raw_history_saved=False)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report.get(k) for k in ('status','execution','failure')}))
    if 'certificate' in report:
        print(json.dumps(dict(ratios={k:v['upper_float'] for k,v in certificate['ratios'].items()},
                              margins={k:v['lower_float'] for k,v in certificate['margins'].items()})))
    return 0 if report['status']=='strict_conditional_bootstrap_image' else 1


if __name__=='__main__':raise SystemExit(main())
