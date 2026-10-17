#!/usr/bin/env python3
"""Check radial quadrature at one authenticated failed stage; no time step."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

import run_nsc_atlas_fixed_endpoint as endpoint
from recursive_horizons import nsc_atlas_cap_fixed_reference as fixed

CASES=((32,2e-13,2e-15),(64,2e-13,2e-15),(128,2e-13,2e-15),
       (64,2e-14,1e-15))
CONFIRMATION_CASES=((32,2.5e-14,1e-15),(64,2.5e-14,1e-15))


def read_seed(receipt,root=endpoint.ROOT):
    binding=receipt['failed_stage_seed']
    path=(root/binding['path']).resolve()
    if not path.is_relative_to((root/'.build/results').resolve()):
        raise ValueError('stage seed must stay under the ignored result directory')
    raw=path.read_bytes()
    if len(raw)!=binding['bytes'] or hashlib.sha256(raw).hexdigest()!=binding['sha256']:
        raise ValueError('stage seed bytes changed')
    state={}
    with np.load(path,allow_pickle=False) as saved:
        if (saved['level'].item()!=receipt['level'] or
            saved['source_commit'].item()!=receipt['source_commit'] or
            saved['x'].item()!=binding['source_x']):
            raise ValueError('stage seed metadata changed')
        x=float(saved['x'].item())
        for key in saved.files:
            if '.' in key:
                group,name=key.split('.',1)
                state.setdefault(group,{})[name]=saved[key].copy()
    return state,x


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--receipt',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--confirm-forward-accuracy',action='store_true',
                        help='the admitted 32/64-point check of the selected tighter tolerance')
    args=parser.parse_args(argv)
    cases=CONFIRMATION_CASES if args.confirm_forward_accuracy else CASES
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('choose a new output in an existing directory')
    raw=args.receipt.read_bytes();receipt=json.loads(raw)
    if receipt['status']!='stopped':
        raise ValueError('require a stopped evolution receipt')
    state,x=read_seed(receipt)
    commit,bindings=endpoint.frozen_sources()
    budget=endpoint.Budget(maximum_evaluations=40_000_000,cpu_seconds=60.)
    report=dict(schema='nsc-atlas-endpoint-radial-diagnostic-v1',status='incomplete',
        source_commit=commit,source_bindings=bindings,
        input_receipt=dict(path=args.receipt.resolve().relative_to(endpoint.ROOT).as_posix(),
                           bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()),
        stage_seed=receipt['failed_stage_seed'],x=x,level=receipt['level'],
        physical_time_steps=0,gradient_tolerance=2e-8,
        comparison_kind='selected_forward_accuracy' if args.confirm_forward_accuracy else 'error_attribution',
        limits=dict(cpu_seconds=60.,charged_values=40_000_000,cases=len(cases)),cases=[])
    original=fixed._radial_variations
    radial_results=[]
    try:
        endpoint.authenticate(endpoint.ROOT,budget)
        source_raw=(endpoint.ROOT/endpoint.producer.BASE_INPUT).read_bytes()
        if hashlib.sha256(source_raw).hexdigest()!=receipt['input']['sha256']:
            raise ValueError('physical source changed')
        base=json.loads(source_raw);geometry=base['geometry'];data=base['initial_data']['finite']
        rows={k:np.asarray(v,dtype=float) for k,v in data['exact_rows'].items()}
        spec=endpoint.level_spec(receipt['level'])
        basis=endpoint.FixedCapBasis(endpoint.cap_knots(receipt['level']),budget=budget)
        background,offset=endpoint.prepare_history(budget,report.setdefault('preparation',{}),
                                                  enrich_ordinary_pieces=True)
        for key in ('input_sha256','original_input_sha256'):
            if report['preparation'][key]!=receipt['preparation'][key]:
                raise ValueError('prepared history changed')
        snapshot=endpoint.SourceSnapshot(geometry,alpha=data['alpha'],x_nodes=rows['x'],
                                         source_rows=rows,budget=budget)
        gauge=endpoint.ReferenceWallGauge(geometry,order=64,budget=budget)
        for order,rtol,atol in cases:
            case=dict(adjoint_order=order,rtol=rtol,atol=atol,
                      effective_rtol=max(rtol,100*np.finfo(float).eps))
            captured=[];started=time.process_time();before=budget.evaluations
            def routed(*a,**kw):
                kw.update(adjoint_order=order)
                result=original(*a,**kw)
                captured.append(result)
                return result
            fixed._radial_variations=routed
            try:
                stage=endpoint.coupled_fixed_stage(state,basis=basis,x=x,
                    eta=np.linspace(0.,1.,spec['exterior_cells']+1),source=snapshot.sample(x),
                    gauge=gauge.sample(x),alpha=data['alpha'],tube=endpoint.TUBE,budget=budget,
                    background=background,background_clock_offset=offset,rtol=rtol,atol=atol)
                cap=stage['connection']['cap']
                case.update(status='passed_existing_gradient_guard',mass=cap['mass'],lapse=cap['lapse'],
                    gradient_discrepancy=cap['gradient_discrepancy'],
                    projected_gradient_discrepancy=cap['projected_gradient_discrepancy'],
                    field_sites=cap['field_sites'],
                    physical_gaps=endpoint.physical_gaps(stage['connection']['residual']))
            except (ValueError,RuntimeError,FloatingPointError,OverflowError) as error:
                case.update(status='stopped',reason=str(error))
                budget.check()
            finally:
                fixed._radial_variations=original
            radial_results.append(captured[-1] if captured else None)
            case.update(cpu_seconds=time.process_time()-started,charged_values=budget.evaluations-before)
            report['cases'].append(case)
        report['pair_differences']=[]
        for i,j in ((0,1),(1,2),(1,3)):
            if j>=len(radial_results):
                continue
            a,b=radial_results[i],radial_results[j]
            if a is not None and b is not None:
                report['pair_differences'].append(dict(cases=[i,j],
                    **{key:float(np.max(np.abs(a[key]-b[key]))) for key in
                       ('mass','lapse','gradient','adjoint_gradient')}))
        report['status']='completed_one_stage_comparison'
    except (ValueError,RuntimeError,FloatingPointError,OverflowError) as error:
        report.update(status='stopped',reason=str(error))
    finally:
        fixed._radial_variations=original
    report['execution']=dict(cpu_seconds=time.process_time()-budget.started,
                              charged_values=budget.evaluations)
    report=endpoint.plain(report)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',
                          encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=report['status'],cases=report['cases'],execution=report['execution'])))
    return 0 if report['status']=='completed_one_stage_comparison' else 1


if __name__=='__main__':
    raise SystemExit(main())
