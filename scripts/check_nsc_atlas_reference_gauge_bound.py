#!/usr/bin/env python3
"""Certify the reference wall tube on the already recorded first source cell."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import time

import flint
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_reference_gauge_bound import reference_gauge_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph

ROOT = Path(__file__).resolve().parents[1]
INPUT = 'results/summaries/nsc-atlas-source-cell-bound-v2.json'
INPUT_SHA = '707bb8d76359edaf2a26228657a219933ebd164edcf89eb87742f9e66bec2bdd'
OWNERS = ('scripts/check_nsc_atlas_reference_gauge_bound.py',
          'src/recursive_horizons/nsc_atlas_reference_gauge_bound.py',
          'src/recursive_horizons/nsc_source_binding_graph.py',
          'scripts/run_python.py','pyproject.toml')


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('require an existing parent and a new summary path')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bindings=original_source_graph(ROOT,OWNERS,commit)
    for name,expected in bindings.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('freeze changed producer before execution: '+name)
    raw=(ROOT/INPUT).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=INPUT_SHA:
        raise ValueError('source certificate identity changed')
    source=json.loads(raw)
    if source['status']!='certified_source_representation_cell':
        raise ValueError('original cell certificate not successful')
    budget=Budget(maximum_evaluations=100000,cpu_seconds=20.)
    report=dict(schema='nsc-atlas-reference-gauge-bound-v1',source_commit=commit,
                source_bindings=bindings,input=dict(path=INPUT,sha256=INPUT_SHA),
                cell=source['cell'],status='incomplete',runtime=dict(python_flint=flint.__version__))
    try:
        certificate=reference_gauge_bound(Q(source['cell']['left']),Q(source['cell']['right']),
                                         Q(474,1000),Q(477,1000),budget=budget)
        report.update(status='certified_exact_reference_gauge_cell',certificate=certificate)
    except (ValueError,RuntimeError,KeyError,TypeError,OverflowError) as error:
        report.update(status='stopped',failure=str(error))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_queries=0,source_integrations=0,time_steps=0,root_searches=0,retries=0,raw_history_saved=False)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({key:report.get(key) for key in ('status','execution','failure')}))
    return 0 if report['status']=='certified_exact_reference_gauge_cell' else 1


if __name__=='__main__':
    raise SystemExit(main())
