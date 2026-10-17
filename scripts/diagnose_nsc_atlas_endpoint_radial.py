#!/usr/bin/env python3
"""Check radial quadrature at one authenticated failed stage; no time step."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

import run_nsc_atlas_fixed_endpoint as endpoint
from recursive_horizons import nsc_atlas_cap_fixed_reference as fixed
from recursive_horizons import nsc_atlas_cap_variational_reference as radial

CASES=((32,2e-13,2e-15),(64,2e-13,2e-15),(128,2e-13,2e-15),
       (64,2e-14,1e-15))
CONFIRMATION_CASES=((32,2.5e-14,1e-15),(64,2.5e-14,1e-15))
PANEL_CASES=((32,2.5e-14,1e-15,.5),(32,2.5e-14,1e-15,.25),
             (64,2.5e-14,1e-15,.5))
ONE_SIDED_CASES=((32,2.5e-14,1e-15,None),(32,2.5e-14,1e-15,.25),
                 (64,2.5e-14,1e-15,None))


def inward_panel_point(point,span):
    """Numerical one-sided limits for a split radial integration panel.

    This diagnostic changes only points within 64 coordinate ulps of an
    endpoint. Refuse a panel too small to resolve that displacement.
    """
    lo,hi=map(float,span)
    margin=64*np.spacing(max(1.,abs(lo),abs(hi)))
    if not lo < hi or margin > 1e-6*(hi-lo):
        raise ValueError('one-sided radial endpoint is not resolved')
    if point <= lo+margin:
        return lo+margin
    if point >= hi-margin:
        return hi-margin
    return point


def panel_step_options(span, options, fraction):
    """Bound the forward step inside each existing smooth radial panel."""
    options=dict(options)
    if fraction is not None:
        if not np.isfinite(fraction) or not 0 < fraction <= 1:
            raise ValueError('radial panel fraction must be in (0, 1]')
        maximum=(float(span[1])-float(span[0]))*fraction
        if not np.isfinite(maximum) or maximum <= 0:
            raise ValueError('positive radial panel required')
        options['max_step']=min(options.get('max_step',np.inf),maximum)
        if 'first_step' in options:
            options['first_step']=min(options['first_step'],options['max_step'])
    return options


_WALL_TRACES=('f','p_trace','h_trace','q_trace','velocity','A','B')


def _safe_ratio(difference, reference):
    if reference==0.:
        return None
    return difference/abs(reference)


def physical_field_comparison(basis,reference,other,budget):
    """Galerkin L2 of one frozen stage's cap rates, plus wall-trace differences.

    ``squared_norm`` integrates the reconstructed field against the model's
    physical-z coordinate measure dz, not proper volume or energy. Diagonal Gram scaling
    stays inside that integral, so a coefficient infinity norm is not reused
    as a physical difference. This reports one stage. It does not decide
    equivalence and it does not replace the coordinate-rate gate.
    """
    left,right=reference['connection']['cap'],other['connection']['cap']
    gram=basis.gram(left['wall_coordinate'],budget=budget)
    def galerkin(coefficients):
        return math.sqrt(gram.squared_norm(np.asarray(coefficients,dtype=float),budget=budget))
    fields={}
    for key in ('chi','canonical_density'):
        original=np.asarray(reference['rates']['cap'][key],dtype=float)
        delta=original-np.asarray(other['rates']['cap'][key],dtype=float)
        reference_norm=galerkin(original)
        difference=galerkin(delta)
        fields[key]=dict(reference_galerkin_L2=reference_norm,difference_galerkin_L2=difference,
                         difference_over_reference=_safe_ratio(difference,reference_norm),
                         coordinate_infinity=float(np.max(np.abs(delta))))
    traces={}
    for key in _WALL_TRACES:
        value=float(left['wall'][key]);delta=abs(value-float(right['wall'][key]))
        traces[key]=dict(reference=value,absolute_difference=delta,
                         difference_over_reference=_safe_ratio(delta,value))
    flux=None
    if 'residual' in reference['connection'] and 'residual' in other['connection']:
        flux=abs(float(reference['connection']['residual']['conformal_wall_flux'])
                 -float(other['connection']['residual']['conformal_wall_flux']))
    return dict(gram_condition=float(gram.condition),
                wall_coordinate_absolute_difference=abs(float(left['wall_coordinate'])
                                                        -float(right['wall_coordinate'])),
                cap_rate_fields=fields,wall_traces=traces,conformal_wall_flux_absolute_difference=flux,
                scope='one frozen stage; Galerkin L2 of reconstructed cap rates on the reference wall Gram; not trajectory equivalence and not a coordinate-rate gate')


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
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument('--confirm-forward-accuracy',action='store_true',
                        help='the admitted 32/64-point check of the selected tighter tolerance')
    mode.add_argument('--panel-step-check',action='store_true',
                      help='half/quarter-panel forward steps at the selected tolerance')
    mode.add_argument('--one-sided-panel-check',action='store_true',
                      help='test inward endpoint limits at existing radial formula breaks')
    mode.add_argument('--confirm-panel-polynomials',action='store_true',
                      help='32/64-point check of explicit panel polynomial ownership at unchanged coordinates')
    mode.add_argument('--compare-compiled-coefficients',action='store_true',
                      help='one frozen stage: reference, first compiled call, warm compiled call; no time steps')
    args=parser.parse_args(argv)
    cases=((CONFIRMATION_CASES[0],)*3 if args.compare_compiled_coefficients else
           ONE_SIDED_CASES if args.one_sided_panel_check else
           PANEL_CASES if args.panel_step_check else
           CONFIRMATION_CASES if (args.confirm_forward_accuracy or args.confirm_panel_polynomials) else CASES)
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('choose a new output in an existing directory')
    raw=args.receipt.read_bytes();receipt=json.loads(raw)
    if receipt['status']!='stopped':
        raise ValueError('require a stopped evolution receipt')
    state,x=read_seed(receipt)
    commit,bindings=endpoint.frozen_sources()
    diagnostic_path=Path(__file__).resolve()
    bindings[diagnostic_path.relative_to(endpoint.ROOT).as_posix()]=hashlib.sha256(
        diagnostic_path.read_bytes()).hexdigest()
    budget=endpoint.Budget(maximum_evaluations=40_000_000,cpu_seconds=60.)
    report=dict(schema='nsc-atlas-endpoint-radial-diagnostic-v1',status='incomplete',
        source_commit=commit,source_bindings=bindings,
        input_receipt=dict(path=args.receipt.resolve().relative_to(endpoint.ROOT).as_posix(),
                           bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()),
        stage_seed=receipt['failed_stage_seed'],x=x,level=receipt['level'],
        physical_time_steps=0,gradient_tolerance=2e-8,
        comparison_kind=('compiled_coefficients' if args.compare_compiled_coefficients else
                         'explicit_panel_polynomials' if args.confirm_panel_polynomials else
                         'one_sided_panel_limits' if args.one_sided_panel_check else
                         'forward_panel_step' if args.panel_step_check else
                         'selected_forward_accuracy' if args.confirm_forward_accuracy else 'error_attribution'),
        limits=dict(cpu_seconds=60.,charged_values=40_000_000,cases=len(cases)),cases=[])
    original=fixed._radial_variations
    original_solve=radial.solve_ivp
    radial_results=[];stage_results=[]
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
        for case_index,specification in enumerate(cases):
            order,rtol,atol=specification[:3]
            fraction=specification[3] if len(specification)==4 else None
            case=dict(adjoint_order=order,rtol=rtol,atol=atol,
                      effective_rtol=max(rtol,100*np.finfo(float).eps))
            if fraction is not None:
                case['maximum_radial_step_panel_fraction']=fraction
            panel_sides=args.confirm_panel_polynomials or args.compare_compiled_coefficients
            compiled=args.compare_compiled_coefficients and case_index>0
            if panel_sides:
                case['background_panel_sides']=True
            if args.compare_compiled_coefficients:
                case['coefficient_backend']=('reference','compiled_first','compiled_warm')[case_index]
            captured=[];started=time.process_time();started_wall=time.perf_counter();before=budget.evaluations
            stage=None
            panels=[];shifts=[]
            def solve(fun,span,initial,**kw):
                def one_sided(point,values):
                    interior=inward_panel_point(point,span)
                    if interior!=point:
                        shifts.append(abs(interior-point))
                    return fun(interior,values)
                result=original_solve(one_sided if args.one_sided_panel_check else fun,
                    span,initial,**panel_step_options(span,kw,fraction))
                budget.check(result.y.size)
                panels.append(dict(nfev=result.nfev,accepted_steps=len(result.t)-1,
                    maximum_internal_gradient=float(np.max(np.abs(result.y[2:])))))
                return result
            def routed(*a,**kw):
                kw.update(adjoint_order=order)
                result=original(*a,**kw)
                captured.append(result)
                return result
            fixed._radial_variations=routed
            radial.solve_ivp=solve
            try:
                stage=endpoint.coupled_fixed_stage(state,basis=basis,x=x,
                    eta=np.linspace(0.,1.,spec['exterior_cells']+1),source=snapshot.sample(x),
                    gauge=gauge.sample(x),alpha=data['alpha'],tube=endpoint.TUBE,budget=budget,
                    background=background,background_clock_offset=offset,rtol=rtol,atol=atol,
                    background_panel_sides=panel_sides,compiled_coefficients=compiled)
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
                radial.solve_ivp=original_solve
            case['forward_panels']=dict(count=len(panels),nfev=sum(p['nfev'] for p in panels),
                accepted_steps=sum(p['accepted_steps'] for p in panels),
                maximum_internal_gradient=max((p['maximum_internal_gradient'] for p in panels),default=0.))
            if args.one_sided_panel_check:
                case['endpoint_limits']=dict(shifted_samples=len(shifts),
                    maximum_coordinate_displacement=max(shifts,default=0.),
                    coordinate_ulps=64,largest_allowed_panel_fraction=1e-6)
            radial_results.append(captured[-1] if captured else None)
            stage_results.append(stage)
            case.update(cpu_seconds=time.process_time()-started,wall_seconds=time.perf_counter()-started_wall,
                        charged_values=budget.evaluations-before)
            report['cases'].append(case)
        report['pair_differences']=[]
        pairs=((0,1),(0,2)) if (args.panel_step_check or args.one_sided_panel_check or args.compare_compiled_coefficients) else ((0,1),(1,2),(1,3))
        for i,j in pairs:
            if j>=len(radial_results):
                continue
            a,b=radial_results[i],radial_results[j]
            if a is not None and b is not None:
                report['pair_differences'].append(dict(cases=[i,j],
                    **{key:float(np.max(np.abs(a[key]-b[key]))) for key in
                       ('mass','lapse','gradient','adjoint_gradient')}))
        if args.compare_compiled_coefficients:
            comparison=[]
            for j in (1,2):
                a,b=stage_results[0],stage_results[j]
                if a is None or b is None:
                    comparison.append(dict(cases=[0,j],passed=False,reason='a stage stopped'))
                    continue
                rates={group+'.'+key:float(np.max(np.abs(value-b['rates'][group][key])))
                       for group,values in a['rates'].items() for key,value in values.items()}
                passed=all(np.allclose(value,b['rates'][group][key],rtol=1e-11,atol=1e-10)
                           for group,values in a['rates'].items() for key,value in values.items())
                scalars={key:abs(a['connection']['cap'][key]-b['connection']['cap'][key])
                         for key in ('mass','lapse')}
                passed=passed and all(value<=1e-11 for value in scalars.values())
                item=dict(cases=[0,j],passed=bool(passed),maximum_rate_differences=rates,
                          scalar_differences=scalars,
                          coordinate_gate='strict coefficient infinity norm and mass/lapse; unchanged')
                try:
                    item['physical_field_comparison']=physical_field_comparison(basis,a,b,budget)
                except (ValueError,RuntimeError,FloatingPointError,OverflowError) as error:
                    item['physical_field_comparison']=dict(status='unavailable',reason=str(error))
                comparison.append(item)
            report['compiled_equivalence']=dict(rate_atol=1e-10,rate_rtol=1e-11,
                mass_lapse_atol=1e-11,comparisons=comparison,
                passed=all(item['passed'] for item in comparison),
                scope='one frozen stage; coordinate-rate failure stays visible; Galerkin fields are not a trajectory equivalence or a replacement gate')
        report['status']='completed_one_stage_comparison'
    except (ValueError,RuntimeError,FloatingPointError,OverflowError) as error:
        report.update(status='stopped',reason=str(error))
    finally:
        fixed._radial_variations=original
        radial.solve_ivp=original_solve
    report['execution']=dict(cpu_seconds=time.process_time()-budget.started,
                              charged_values=budget.evaluations)
    report=endpoint.plain(report)
    args.output.write_text(json.dumps(report,indent=2,sort_keys=True,allow_nan=False)+'\n',
                          encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=report['status'],cases=report['cases'],execution=report['execution'])))
    return 0 if report['status']=='completed_one_stage_comparison' else 1


if __name__=='__main__':
    raise SystemExit(main())
