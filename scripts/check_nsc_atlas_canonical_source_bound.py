#!/usr/bin/env python3
"""Carry the existing source-cell certificate through the canonical frame."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import subprocess
import time

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_transmission_bound import canonical_source_error_bound
from recursive_horizons.nsc_source_binding_graph import original_source_graph


ROOT = Path(__file__).resolve().parents[1]
INPUT = 'results/summaries/nsc-atlas-source-cell-bound-v2.json'
INPUT_SHA = '707bb8d76359edaf2a26228657a219933ebd164edcf89eb87742f9e66bec2bdd'
OWNERS = ('scripts/check_nsc_atlas_canonical_source_bound.py',
          'src/recursive_horizons/nsc_atlas_transmission_bound.py',
          'src/recursive_horizons/nsc_source_binding_graph.py',
          'src/recursive_horizons/nsc_atlas_characteristic_initialization.py',
          'scripts/run_python.py', 'pyproject.toml')


def dyadic(row):
    value = row['upper_dyadic']
    exponent = value['exponent']
    if type(exponent) is not int:
        raise ValueError('invalid dyadic exponent')
    result = Q(int(value['mantissa']))*Q(2)**exponent
    if result < 0:
        raise ValueError('negative error envelope')
    return result


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('require an existing parent and a new summary path')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    bindings=original_source_graph(ROOT,OWNERS,commit)
    for name, expected in bindings.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('freeze changed producer before execution: '+name)
    raw=(ROOT/INPUT).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=INPUT_SHA:
        raise ValueError('source certificate identity changed')
    budget=Budget(maximum_evaluations=4096,cpu_seconds=5.)
    report=dict(schema='nsc-atlas-canonical-source-bound-v1',source_commit=commit,
                source_bindings=bindings,input=dict(path=INPUT,sha256=INPUT_SHA),status='incomplete')
    try:
        budget.check(512)
        original=json.loads(raw); certificate=original['certificate']; tube=certificate['tube']
        if (original['status']!='certified_source_representation_cell'
                or not tube['exact_continuation_proved'] or not tube['polynomial_inside']
                or not certificate['right_owner']['supplied']):
            raise ValueError('source certificate does not supply the required common tube')
        F=Q(tube['F_lower']); rJ=Q(tube['rJ_lower']); vJ=Q(tube['absolute_rJ_x_upper'])
        if not 0<F<=1 or rJ<=0 or vJ<0:
            raise ValueError('source tube is not in the declared positive frame')
        alpha=Q(original['cell']['alpha'])
        if alpha<0:
            raise ValueError('negative source coupling')
        g=dyadic(certificate['analytic_source_bounds']['g'])
        gx=dyadic(certificate['analytic_source_bounds']['g_x'])
        radius=F*rJ
        common=dict(alpha=alpha,radius_floor=radius,u_bound=g/radius,
                    c_x_bound=gx,r_x_bound=vJ/F+alpha*g*gx/(6*radius))
        errors={k:dyadic(certificate['frame_error'][k]) for k in ('r','c','r_x','c_x')}
        out=canonical_source_error_bound(**common,errors=errors,budget=budget)
        report.update(status='certified_canonical_source_representation',
                      cell=original['cell'],common_tube={k:str(v) for k,v in common.items()},
                      input_errors={k:str(v) for k,v in errors.items()},
                      bound={k:str(v) for k,v in out.items()},
                      display={k:float(out[k]) for k in ('value','first')},
                      scope='source representation only; floating evaluation and coupled propagation additional')
    except (ValueError,RuntimeError,KeyError,TypeError,OverflowError) as error:
        report.update(status='stopped',failure=str(error))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_queries=0,source_integrations=0,time_steps=0,retries=0,raw_history_saved=False)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps({k:report.get(k) for k in ('status','display','execution','failure')}))
    return 0 if report['status']=='certified_canonical_source_representation' else 1


if __name__=='__main__':
    raise SystemExit(main())
