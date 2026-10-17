#!/usr/bin/env python3
"""One finite-amplitude material-radius frame conversion, without replay."""
import argparse
import hashlib
from importlib.metadata import version
import json
from fractions import Fraction as Q
from pathlib import Path
import subprocess
import sys
import time

from check_nsc_atlas_onset_endpoint_bound import plain,read_bound
from check_nsc_atlas_geometric_remainder_bound import outward
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_jordan_radius_bound import jordan_radius_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT=Path(__file__).resolve().parents[1]
INPUTS={
 'results/summaries/nsc-atlas-coupled-tube-bound-v1.json':'192caf213cd4b67961f815203687a87748a9668c4a97d2d77f9ee2f5d78a8dab',
 'results/summaries/nsc-atlas-geometric-departure-bound-v1.json':'1c6fab9b40ba7fb623dffca80ec8e06b859531ccec3497361649269e580afd62',
 'results/summaries/nsc-atlas-scalar-departure-bound-v1.json':'fb52782856689b88ce2d6bee893085352eda9c380372ec178d2aa670196ef7f5',
 'results/summaries/nsc-atlas-geometric-remainder-bound-v2.json':'b7416d98046cd73b4b76d984de4ea634726e28700f8375b835d6afc62ffe86b2',
 'results/summaries/nsc-atlas-onset-geometry-bound-v1.json':'a8f84fd4ede226068134a3e7b05379524d0c5f683e116318a6cc3c35e101df6c'}
OWNERS=('scripts/check_nsc_atlas_jordan_radius_bound.py',
 'src/recursive_horizons/nsc_atlas_jordan_radius_bound.py',
 'src/recursive_horizons/nsc_source_binding_graph.py','scripts/run_python.py','pyproject.toml')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():raise ValueError('new output path required')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bindings=original_source_graph(ROOT,OWNERS,commit)
    for n,h in bindings.items():
        if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=h:raise ValueError('freeze producer before execution: '+n)
    records=[read_bound(n,h) for n,h in INPUTS.items()]
    tube,geometry,scalar,remainder,first=records
    expected=('strict_conditional_bootstrap_image','certified_amplitude_geometric_departure',
        'certified_amplitude_scalar_departure','certified_geometric_quadratic_remainder','certified_first_geometric_endpoint')
    if tuple(r['status'] for r in records)!=expected or tube['uniform_in_coupling'] is not True:
        raise ValueError('certified original onset chain required')
    authenticated=dict(INPUTS)
    for record in records:
        rows=record.get('inputs',[])
        if isinstance(rows,dict):rows=list(rows.values())
        for row in rows:
            n,h=row['path'],row['sha256']
            if n in authenticated:
                if authenticated[n]!=h:raise ValueError('inconsistent input: '+n)
                continue
            with (ROOT/n).open('rb') as stream:
                if hashlib.file_digest(stream,'sha256').hexdigest()!=h:raise ValueError('changed input: '+n)
            authenticated[n]=h
    for k in ('left','right'):
        if Q(first['parameters'][k])!=Q(remainder['cell'][k]):raise ValueError('original endpoint cell mismatch')
    budget=Budget(maximum_evaluations=12288,cpu_seconds=10.)
    report=dict(schema='nsc-atlas-jordan-radius-bound-v1',status='incomplete',source_commit=commit,
        source_bindings=bindings,inputs=[dict(path=n,sha256=h) for n,h in authenticated.items()],
        cell=remainder['cell'],runtime=dict(python=sys.version.split()[0],python_flint=version('python-flint'),platform=sys.platform),
        limits=dict(charged_values=12288,cpu_seconds=10.,precision_bits=160,applications=1))
    try:
        result=jordan_radius_bound(*(r['certificate'] for r in records),budget=budget)
        endpoints={}
        for name,row in result['endpoint'].items():
            endpoints[name]=dict(strict_on_full_interval=row['strict_on_full_interval'],
                strict_sign_alpha_max=outward(row['strict_sign_alpha_max'],True),
                relative_error_upper=outward(row['relative_error_at_maximum']),
                normalized_interval=dict(lower=outward(row['normalized_lower_at_maximum'],True),upper=outward(row['normalized_upper_at_maximum'])),
                shift_interval=dict(lower=outward(row['shift_lower_at_maximum'],True),upper=outward(row['shift_upper_at_maximum'])))
        report.update(status='certified_jordan_radius_remainder',certificate=plain(result),
            display=dict(first_jordan=dict(lower=outward(result['first_jordan']['lower'],True),upper=outward(result['first_jordan']['upper'])),
                frame_half_hessian=outward(result['frame']['half_hessian']),
                conformal_factor_lower=outward(result['frame']['conformal_factor_lower'],True),
                coefficients={k:outward(v) for k,v in result['coefficients'].items()},
                maximum_amplitude_remainders={k:outward(v) for k,v in result['maximum_amplitude_remainders'].items()},
                endpoint=endpoints,common_strict_sign_alpha_max=outward(result['common_strict_sign_alpha_max'],True),
                strict_endpoint_sign_on_full_interval=result['strict_endpoint_sign_on_full_interval']))
        budget.check()
    except (ValueError,RuntimeError,KeyError,TypeError,OverflowError) as error:
        report.update(status='stopped',failure=str(error))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_queries=0,field_evolutions=0,integrations=0,scans=0,retries=0,raw_history_saved=False)
    report['regenerate']='python scripts/run_python.py scripts/check_nsc_atlas_jordan_radius_bound.py --output <new-output-path>'
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report.get(k) for k in ('status','display','execution','failure')}))
    return 0 if report['status']=='certified_jordan_radius_remainder' else 1


if __name__=='__main__':raise SystemExit(main())
