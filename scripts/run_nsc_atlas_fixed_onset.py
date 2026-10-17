#!/usr/bin/env python3
"""One retained-source onset window; summaries/hashes, never raw histories."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

from diagnose_nsc_atlas_boundary_rows import authenticate,producer,plain,constraint_terms
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_trace import SourceSnapshot
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_linear_bound import frozen_cap_bound
from recursive_horizons.nsc_atlas_fixed_cap_stage import quiet_fixed_state,coupled_fixed_stage
from recursive_horizons.nsc_atlas_fixed_cap_time import (
    advance_fixed_ssprk3,fixed_time_residual,fixed_derived_residual,
)
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget,solve_prepared_field
from recursive_horizons.nsc_atlas_cap_field import CapLeadingField
from recursive_horizons.nsc_atlas_background_profile import build_prepared_background
from recursive_horizons.nsc_source_binding_graph import original_source_graph
from check_nsc_atlas_prepared_history import array_hash

ROOT=Path(__file__).resolve().parents[1]
ENRICHED_REFERENCE='results/summaries/nsc-atlas-prepared-endpoint-enriched-v1.json'


def digest_state(value):
    return hashlib.sha256(json.dumps(plain(value),sort_keys=True,separators=(',',':'),
                                    allow_nan=False).encode()).hexdigest()


def frozen_sources():
    """Pin imported local owners and their AST dependencies before computing."""
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    owners=set()
    for module in tuple(sys.modules.values()):
        name=getattr(module,'__file__',None)
        if name is None: continue
        path=Path(name).resolve()
        if path.is_relative_to(ROOT) and path.suffix=='.py':
            relative=path.relative_to(ROOT)
            if relative.parts[0] in ('src','scripts'): owners.add(relative.as_posix())
    owners.update(('pyproject.toml','scripts/run_python.py'))
    bindings=original_source_graph(ROOT,tuple(sorted(owners)),commit)
    for name,expected in bindings.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise ValueError('commit the calculation source before execution: '+name)
    return commit,bindings


def prepare_history(budget,report,*,enrich_ordinary_pieces=False):
    """One preparation phase, charged into the original caller budget once."""
    budget.check()
    remaining=budget.maximum_evaluations-budget.evaluations
    if remaining<=0: raise RuntimeError('no remaining preparation value budget')
    preparation=PreparedFieldBudget(cpu_seconds=min(60.,budget.cpu_seconds),
                                    kernel_values_max=min(10000000,remaining))
    preparation.started=budget.started
    report.update(source_integrations_attempted=1,source_integrations_completed=0,
                  phase_limits=dict(charged_values=preparation.kernel_values_max,
                                    cpu_seconds_from_run_start=preparation.cpu_seconds))
    try:
        field=solve_prepared_field(nodes=257,geometry_knots=257,budget=preparation)
        report['source_integrations_completed']=1
        cap=CapLeadingField(field,derivative_order=4)
        result=(build_prepared_background(cap,enrich_ordinary_pieces=True) if enrich_ordinary_pieces
                else build_prepared_background(cap))
        background=result['background']
        report.update(inputs=dict(field_nodes=257,geometry_knots=257,derivative_order=4,
                                  quiet_tolerances=[0.,0.,0.]),
                      input_sha256=array_hash(background.inputs),onset=result['onset'],
                      diagnostics=result['diagnostics'],clock_offset=cap.domain['Ta'],
                      numerical_error_enclosure=False)
        if enrich_ordinary_pieces:
            report.update(original_input_sha256=array_hash(result['original_history_inputs']),
                          enrichment=result['enrichment'])
        return background,float(cap.domain['Ta'])
    finally:
        report['charged_values']=preparation.evaluations
        report['cpu_seconds_from_run_start']=time.process_time()-budget.started
        budget.check(preparation.evaluations)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--retained-history',action='store_true',
                        help='same full interval and schedule, with the prepared background and departures')
    parser.add_argument('--enrich-history',action='store_true',help='use the checked degree-seven ordinary history pieces')
    args=parser.parse_args()
    if args.enrich_history and not args.retained_history:
        parser.error('--enrich-history requires --retained-history')
    if args.output.exists(): raise ValueError('preserve the existing result; choose a new admitted output')
    if not args.output.parent.is_dir(): raise ValueError('output parent must already exist')
    commit,bindings=frozen_sources()
    values,cpu=(2000000000,900.) if args.retained_history else (600000000,180.)
    budget=Budget(maximum_evaluations=values,cpu_seconds=cpu)
    schema=('nsc-atlas-fixed-enriched-history-onset-v1' if args.enrich_history else
            'nsc-atlas-fixed-history-onset-v1' if args.retained_history else 'nsc-atlas-fixed-onset-v1')
    report=dict(schema=schema,
        status='incomplete',runtime=producer.runtime_versions(),source_commit=commit,
        scope='retained-source numerical common evolution; sampled defects, no continuum error enclosure or observational prediction',
        limits=dict(charged_values=values,cpu_seconds=cpu,accepted_steps=1024),source_bindings=bindings)
    phase='input authentication';calls=accepted=0;maxima={};minimum_A=None;last_x=None
    history=hashlib.sha256();state=None;background=None;clock_offset=0.
    def maximum(key,value):
        maxima[key]=max(maxima.get(key,0.),float(np.max(np.abs(value))))
    try:
        _,_,_,identity=authenticate(ROOT,budget)
        raw=(ROOT/producer.BASE_INPUT).read_bytes()
        if hashlib.sha256(raw).hexdigest()!=identity['source_bindings'][producer.BASE_INPUT]:
            raise ValueError('source input changed')
        base=json.loads(raw);geometry=base['geometry'];data=base['initial_data']['finite']
        rows={k:np.asarray(v,dtype=float)[::2] for k,v in data['exact_rows'].items()}
        x=rows['x'];front=geometry['support_x'][0];j=int(np.searchsorted(x,front)-1)
        if j!=35 or len(x)!=129: raise ValueError('retained onset interval changed')
        left,right=float(x[j]),float(x[j+1]);alpha=float(data['alpha'])
        snapshot=SourceSnapshot(geometry,alpha=alpha,x_nodes=x,source_rows=rows,budget=budget)
        gauge=ReferenceWallGauge(geometry,order=64,budget=budget)
        if args.retained_history:
            phase='one prepared history and transfer';report['preparation']={}
            background,clock_offset=prepare_history(budget,report['preparation'],enrich_ordinary_pieces=args.enrich_history)
            if args.enrich_history:
                phase='checked enriched input identity'
                reference_raw=(ROOT/ENRICHED_REFERENCE).read_bytes()
                pinned=subprocess.check_output(['git','show',commit+':'+ENRICHED_REFERENCE],cwd=ROOT)
                if reference_raw!=pinned: raise ValueError('enriched reference differs from committed receipt')
                checked=json.loads(reference_raw)
                if (checked['schema']!='nsc-atlas-prepared-endpoint-enriched-v1'
                        or checked['status']!='completed_reference_comparison'
                        or set(checked['acceptance'])!={'history_material_rate','conformal_wall_gap',
                                                       'prepared_prefix_gap','velocity_channel','density_force_L2'}
                        or not all(item['passed'] for item in checked['acceptance'].values())
                        or any(report['preparation'][key]!=checked['preparation'][key]
                               for key in ('input_sha256','original_input_sha256'))):
                    raise ValueError('enriched history differs from its checked preparation')
                report['enriched_reference']=dict(path=ENRICHED_REFERENCE,sha256=hashlib.sha256(reference_raw).hexdigest())
        reference=gauge.sample(np.array([left,right]));Z0,Z1=np.arctanh(reference['r'])
        dZ=Q(float(Z1))-Q(float(Z0));dT=Q(float(reference['T'][1]))-Q(float(reference['T'][0]))
        a=max(4*dZ,2*dT);knots=[Q(0),Q(float(Z0))-a,Q(float(Z0))+a]
        if not 0<knots[1]<Q(float(Z0))<knots[2]: raise ValueError('invalid onset geometry scale')
        basis=FixedCapBasis(knots,budget=budget)
        phase='initial frozen operator'
        bound=frozen_cap_bound(basis,Q(float(Z0)),budget=budget)
        L=bound['integer_bound'];tx_ceiling=math.ceil(float(np.max(reference['T_x'])))
        pieces=[(left,front),(front,right)]
        counts=[math.ceil((b-a)*16*L*tx_ceiling) for a,b in pieces]
        report.update(input=dict(path=producer.BASE_INPUT,sha256=hashlib.sha256(raw).hexdigest(),
            operator_sha256=identity['sha256'],manifest_sha256=identity['manifest_sha256'],
            authenticated_input_files=104,source_knots=129,alpha=alpha),
            geometry=dict(x_left=left,scalar_front=front,x_right=right,
                reference_radius=reference['r'],reference_T=reference['T'],
                reference_delta_Z=float(dZ),reference_delta_T=float(dT),
                reference_incoming_label_advance=float(dT+dZ),
                exact_cap_knots=[str(v) for v in knots],cap_cells=2,cap_dimension=basis.dimension,exterior_cells=64),
            schedule=dict(piece_steps=counts,initial_linear_bound_squared=str(bound['bound_squared']),
                initial_linear_integer_bound=L,T_x_ceiling=tx_ceiling,maximum_initial_scaled_step=1/16,
                scope='initial frozen alpha=0 time scale; moving coupled stability not established'))
        if sum(counts)>1024: raise ValueError('derived time schedule exceeds admitted steps')
        print(json.dumps(dict(event='prepared',steps=sum(counts),L=L,T_x_ceiling=tx_ceiling)),flush=True)
        eta=np.linspace(0.,1.,65);state=quiet_fixed_state(basis,exterior_cells=64)
        patch=np.r_[np.arange(28,37),np.arange(56,64)]
        report['selected_exterior_constraint_nodes']=patch.tolist()
        last_x=left

        def stage(value,*,x,budget):
            nonlocal calls,minimum_A
            calls+=1
            source=snapshot.sample(x);ref=gauge.sample(x)
            if background is None and x<=front and (any(np.any(v!=0.) for v in value['cap'].values())
                             or np.any(value['departure']['c']!=0.)):
                raise ValueError('scalar zero-source recovery failed before the front')
            row=coupled_fixed_stage(value,basis=basis,x=x,eta=eta,source=source,gauge=ref,
                alpha=alpha,tube=dict(r=[.2,.8],L=[-2.,2.],c=[-4.,4.]),budget=budget,
                background=background,background_clock_offset=clock_offset)
            connection=row['connection'];cap=connection['cap'];tx=connection['wall_rates']['T']
            if not 0<tx<=tx_ceiling: raise ValueError('actual cap clock exits declared ceiling')
            floor=cap['metric_certificate']['A_lower_float_down']
            minimum_A=floor if minimum_A is None else min(minimum_A,floor)
            for key,v in connection['residual'].items(): maximum('connection.'+key,v)
            for group in ('wall_kinematic_defect','shared_wall_update_defect','wall_incoming_jet_defect'):
                for key,v in row[group].items(): maximum(group+'.'+key,v)
            for key in ('source_constraint','source_h_xx_defect','source_constraint_jordan',
                        'source_constraint_from_curvature_defect'):
                maximum('source.'+key,source[key])
            maximum('source.stage_minus_snapshot',row['source_constraint']-source['source_constraint'])
            maximum('cap.weak_balance_defect',cap['weak_balance_defect'])
            maximum('cap.projected_gradient_wall_velocity',cap['projected_gradient_discrepancy']['wall_velocity'])
            if background is not None:
                maximum('cap.explicit_time_remainder',cap['explicit_time_remainder'])
                maximum('cap.remainder_identity_defect',cap['remainder_identity_defect'])
                maximum('cap.explicit_mass_clock_derivative',cap['gradient']['T'])
                maximum('cap.background_radial_intervals',cap['background']['radial_intervals'])
                if x<=front:
                    maximum('quiet.cap_wall_scalar',row['physical_wall']['c'])
                    maximum('quiet.exterior_scalar_departure',value['departure']['c'])
                    for key,v in value['cap'].items(): maximum('quiet.cap_coefficients.'+key,v)
            # Reuse the existing independent pure-U constraint algebra on
            # its original centre/wall patches, without a history replay.
            vac=row['exterior_prefix']['vacuum']
            offset=float(source['r'])-float(vac['r'][0])
            delta_r=np.r_[0.,value['departure']['r'],row['physical_wall']['r']-vac['r'][-1]-offset]
            second=(delta_r[2:]-2*delta_r[1:-1]+delta_r[:-2])*64**2
            pure=constraint_terms({k:row['exterior'][k][patch] for k in ('r','c')},
                {k:vac[k][patch] for k in ('r_Ubar','r_Ubar_Ubar','log_C_Ubar')},
                {k:row['departure_gradients'][k][patch] for k in ('r','L','c')},
                second[patch-1],float(ref['W']),alpha,budget)
            maximum('exterior.selected_pure_U_constraint',pure['E_Ubar'])
            maximum('identity.selected_pure_U_decomposition',pure['algebraic_identity_defect'])
            for key in ('r','L','c'):
                defect=row['wall_kinematic_defect'][key]+ref['W_x']*row['wall_incoming_jet_defect'][key]
                maximum('identity.wall.'+key,defect)
                if abs(defect)>2e-11: raise ValueError('wall jet/kinematic identity failed')
            R=row['physical_wall']['r'];rootD=math.sqrt(float(connection['hybrid']['D']))
            defect=connection['residual']['conformal_wall_flux']-R*rootD*connection['residual']['scalar_tangent']
            maximum('identity.conformal_scalar',defect)
            if abs(defect)>2e-11: raise ValueError('conformal scalar tangent identity failed')
            row['onset_reference']=ref
            return row

        history.update(digest_state(dict(x=left,state=state)).encode())
        first_row=final_row=None
        if background is not None:
            # This starts at the original left, not a presumed continuation
            # of a skipped numerical quiet prefix. The whole cap is quiet.
            initial_label=float(reference['T'][0]+clock_offset+Z0)
            onset=float(background.inputs['labels'][background.onset['index']])
            report['initial_history_quiet_margin']=onset-initial_label
            if not initial_label<onset: raise ValueError('initial cap history is not entirely quiet')
        for piece,(pa,pb) in enumerate(pieces):
            nodes=np.linspace(pa,pb,counts[piece]+1)
            for index in range(counts[piece]):
                xa,xb=float(nodes[index]),float(nodes[index+1]);xm=xa+(xb-xa)/2
                phase='common step'
                step=advance_fixed_ssprk3(state,xa,xm,xb,basis=basis,exterior_cells=64,stage=stage,budget=budget)
                if first_row is None:
                    first_row=step['stages'][0]
                    report['initial']=dict(wall=first_row['physical_wall'],
                        maximum_absolute_rate=max(float(np.max(abs(v))) for g in first_row['rates'].values() for v in g.values()),
                        source_constraint=first_row['source_constraint'])
                    if report['initial']['maximum_absolute_rate']>1e-10:
                        raise ValueError('initial quiet value/first-jet recovery failed')
                phase='continuous midpoint residual'
                dense=fixed_time_residual(step,xm,basis=basis,exterior_cells=64,stage=stage,budget=budget)
                dr=dense['residual'];physical=dense['stage'];R=physical['physical_wall']['r']
                Z=float(np.arctanh(R));gram=basis.gram(Z,budget=budget)
                for key in ('chi','canonical_density'):
                    maximum('time_L2.'+key,math.sqrt(gram.squared_norm(dr['cap'][key],budget=budget)))
                for group in ('departure','wall'):
                    for key,v in dr[group].items(): maximum('time_component.'+group+'.'+key,v)
                derived=fixed_derived_residual(dr,physical,basis=basis,budget=budget,background=background)
                for key,v in derived.items(): maximum('time_derived.'+key,v)
                for item in step['stages']:
                    maximum('exterior_characteristic_courant',(xb-xa)*item['exterior_characteristic_courant_rate'])
                state=step['accepted_state'];final_row=step['stages'][3];last_x=xb;accepted+=1
                history.update(digest_state(dict(x=xb,state=state)).encode())
                if accepted%16==0:
                    print(json.dumps(dict(event='progress',accepted=accepted,x=xb,
                        charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started)),flush=True)
                budget.check()
        phase='endpoint summary'
        cap=final_row['connection']['cap'];Z=float(np.arctanh(final_row['physical_wall']['r']))
        gram=basis.gram(Z,budget=budget)
        report['final']=dict(wall=final_row['physical_wall'],wall_departures=state['wall'],
            state=state,wall_rates=final_row['wall_rates'],
            cap_wall_diagnostics=cap['wall'],
            projected_gradient_discrepancy=cap['projected_gradient_discrepancy'],
            normalized_cap_mass=cap['mass'],cap_lapse=cap['lapse'],
            cap_L2={k:math.sqrt(gram.squared_norm(v,budget=budget)) for k,v in state['cap'].items()},
            exterior_departure_max={k:float(np.max(abs(v))) for k,v in state['departure'].items()},
            connection_residuals=final_row['connection']['residual'],
            source_constraint=final_row['source_constraint'],
            state_sha256=digest_state(state))
        if background is not None:
            report['final']['background']=dict(time=cap['background']['time'],
                explicit_time_remainder=cap['explicit_time_remainder'],
                explicit_mass_clock_derivative=cap['gradient']['T'])
        if calls!=5*accepted: raise ValueError('common stage accounting mismatch')
        if args.enrich_history:
            measured=maxima['connection.conformal_wall_flux']
            report['scalar_boundary_decision']=dict(maximum_conformal_wall_flux=measured,
                limit=1e-10,passed=measured<=1e-10,
                scope='sampled scalar boundary decision only; no continuum or complete junction error enclosure')
        budget.check();report['status']='completed_numerical_interval'
    except (ValueError,RuntimeError,FloatingPointError,OverflowError,KeyError,TypeError) as error:
        report['status']='stopped';report['failure']=dict(phase=phase,reason=str(error))
    report.update(maximum_absolute=maxima,minimum_stage_metric_certificate=minimum_A,
        accepted_history_sha256=history.hexdigest(),
        execution=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
            common_stage_calls=calls,accepted_steps=accepted,last_accepted_x=last_x,
            source_integrations=report.get('preparation',{}).get('source_integrations_completed',0),
            source_integrations_attempted=report.get('preparation',{}).get('source_integrations_attempted',0),
            retries=0,raw_history_saved=False,final_state_saved='final' in report))
    args.output.write_text(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False)+'\n',
                           encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=report['status'],execution=report['execution'],
                         failure=report.get('failure'))),flush=True)
    return 0 if report['status']=='completed_numerical_interval' else 1


if __name__=='__main__': raise SystemExit(main())
