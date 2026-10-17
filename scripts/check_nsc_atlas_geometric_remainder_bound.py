#!/usr/bin/env python3
"""One interval-Hessian and positive-matrix remainder application."""
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
from recursive_horizons.nsc_atlas_geometric_remainder_bound import geometric_remainder_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT=Path(__file__).resolve().parents[1]
INPUTS={
 'results/summaries/nsc-atlas-source-forcing-bound-v1.json':'2ccf4ad7326291b79500f4d78a0f858907a80901c74b771f8b22fd8cb735c611',
 'results/summaries/nsc-atlas-reference-gauge-bound-v1.json':'e8a8400cec0f6ae0bd831723696ebf1d41c38bdf0b052301ffab837465f1d216',
 'results/summaries/nsc-atlas-coupled-tube-bound-v1.json':'192caf213cd4b67961f815203687a87748a9668c4a97d2d77f9ee2f5d78a8dab',
 'results/summaries/nsc-atlas-geometric-departure-bound-v1.json':'1c6fab9b40ba7fb623dffca80ec8e06b859531ccec3497361649269e580afd62',
 'results/summaries/nsc-atlas-scalar-departure-bound-v1.json':'fb52782856689b88ce2d6bee893085352eda9c380372ec178d2aa670196ef7f5',
 'results/summaries/nsc-atlas-cap-departure-bound-v1.json':'61190f1d3029ccda920c0a6783d9e9d682408e43a3078a816672a182502e44d2',
 'results/summaries/nsc-atlas-source-remainder-bound-v1.json':'56d8c15170909d2e8fdd65e96bcfcc9fe2f082b733deb58bf706fbd3b38be721',
 'results/summaries/nsc-atlas-onset-geometry-bound-v1.json':'a8f84fd4ede226068134a3e7b05379524d0c5f683e116318a6cc3c35e101df6c'}
OWNERS=('scripts/check_nsc_atlas_geometric_remainder_bound.py',
 'src/recursive_horizons/nsc_atlas_geometric_remainder_bound.py',
 'src/recursive_horizons/nsc_source_binding_graph.py','scripts/run_python.py','pyproject.toml')


def outward(value,lower=False):
    return 0. if value==0 else math.nextafter(float(value),-math.inf if lower else math.inf)


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
    records=[read_bound(name,digest) for name,digest in INPUTS.items()]
    source,gauge,tube,geometry,scalar,cap,source2,first=records
    authenticated=dict(INPUTS)
    for record in records:
        rows=record.get('inputs',[])
        if isinstance(rows,dict):rows=list(rows.values())
        if 'input' in record:rows=[*rows,record['input']]
        for row in rows:
            name,digest=row['path'],row['sha256']
            if name in authenticated:
                if authenticated[name]!=digest:raise ValueError('inconsistent original input hash: '+name)
                continue
            # Authenticate historical payloads without parsing or sampling
            # their arrays. The calculation consumes only compact receipts.
            with (ROOT/name).open('rb') as stream:
                if hashlib.file_digest(stream,'sha256').hexdigest()!=digest:
                    raise ValueError('original input changed: '+name)
            authenticated[name]=digest
    expected=('certified_actual_source_forcing','certified_exact_reference_gauge_cell',
        'strict_conditional_bootstrap_image','certified_amplitude_geometric_departure',
        'certified_amplitude_scalar_departure','certified_amplitude_cap_departure',
        'certified_quadratic_source_remainder','certified_first_geometric_endpoint')
    if tuple(row['status'] for row in records)!=expected or tube['uniform_in_coupling'] is not True:
        raise ValueError('require original chain of certified onset inputs')
    for key in ('left','right'):
        values=(source['parameters'][key],gauge['certificate']['inputs'][key],first['parameters'][key])
        if len(set(map(Q,values)))!=1:raise ValueError('first response/source cell mismatch')
    budget=Budget(maximum_evaluations=49152,cpu_seconds=15.)
    report=dict(schema='nsc-atlas-geometric-remainder-bound-v1',status='incomplete',
        source_commit=commit,source_bindings=bindings,
        inputs=[dict(path=n,sha256=h) for n,h in authenticated.items()],
        runtime=dict(python=sys.version.split()[0],python_flint=version('python-flint'),platform=sys.platform),
        limits=dict(charged_values=49152,cpu_seconds=15.,precision_bits=160,applications=1,
                    interval_passes=1,hessian_maps=4,matrix_dimension=3),
        cell={k:source['parameters'][k] for k in ('left','right','alpha_max')})
    try:
        result=geometric_remainder_bound(*(record['certificate'] for record in records),budget=budget)
        endpoint={}
        for name,row in result['endpoint'].items():
            endpoint[name]=dict(strict_on_full_interval=row['strict_on_full_interval'],
                strict_sign_alpha_max=outward(row['strict_sign_alpha_max'],True),
                relative_error_upper=outward(row['relative_error_at_maximum']),
                normalized_interval=dict(lower=outward(row['normalized_lower_at_maximum'],True),
                                         upper=outward(row['normalized_upper_at_maximum'])),
                shift_interval=dict(lower=outward(row['shift_lower_at_maximum'],True),
                                    upper=outward(row['shift_upper_at_maximum'])))
        report.update(status='certified_geometric_quadratic_remainder',certificate=plain(result),
            display=dict(coefficients={k:outward(v) for k,v in result['coefficients'].items()},
                maximum_amplitude_remainders={k:outward(v) for k,v in result['maximum_amplitude_remainders'].items()},
                endpoint=endpoint,common_strict_sign_alpha_max=outward(result['common_strict_sign_alpha_max'],True),
                strict_endpoint_sign_on_full_interval=result['strict_endpoint_sign_on_full_interval']))
        budget.check()
    except (ValueError,RuntimeError,KeyError,TypeError,OverflowError) as error:
        report.update(status='stopped',failure=str(error))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_queries=0,source_integrations=0,field_evolutions=0,time_steps=0,root_searches=0,
        retries=0,raw_history_saved=False)
    report['regenerate']=('python scripts/run_python.py scripts/check_nsc_atlas_geometric_remainder_bound.py'
                          ' --output <new-output-path>')
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report.get(k) for k in ('status','display','execution','failure')}))
    return 0 if report['status']=='certified_geometric_quadratic_remainder' else 1


if __name__=='__main__':
    raise SystemExit(main())
