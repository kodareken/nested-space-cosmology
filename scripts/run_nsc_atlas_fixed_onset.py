#!/usr/bin/env python3
"""One retained-source onset window; summaries/hashes, never raw histories."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from diagnose_nsc_atlas_boundary_rows import authenticate,producer,plain,constraint_terms
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_trace import SourceSnapshot
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_linear_bound import frozen_cap_bound
from recursive_horizons.nsc_atlas_fixed_cap_stage import quiet_fixed_state,coupled_fixed_stage
from recursive_horizons.nsc_atlas_fixed_cap_time import advance_fixed_ssprk3,fixed_time_residual

ROOT=Path(__file__).resolve().parents[1]
SOURCE_NAMES=('nsc_atlas_characteristic_initialization','nsc_atlas_characteristic_trace',
    'nsc_atlas_reference_wall_gauge','nsc_atlas_cap_fixed_basis','nsc_atlas_cap_linear_bound',
    'nsc_atlas_cap_fixed_bound','nsc_atlas_cap_fixed_reference','nsc_atlas_cap_fixed_hybrid',
    'nsc_atlas_fixed_cap_stage','nsc_atlas_fixed_cap_time','nsc_atlas_vacuum_departure',
    'nsc_atlas_departure_gradient','nsc_atlas_nonlinear_boundary')


def digest_state(value):
    return hashlib.sha256(json.dumps(plain(value),sort_keys=True,separators=(',',':'),
                                    allow_nan=False).encode()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    budget=Budget(maximum_evaluations=600000000,cpu_seconds=180.)
    report=dict(schema='nsc-atlas-fixed-onset-v1',status='incomplete',runtime=producer.runtime_versions(),
        scope='retained-source numerical common evolution; sampled defects, no continuum error enclosure or observational prediction',
        limits=dict(charged_values=600000000,cpu_seconds=180.,accepted_steps=1024),
        source_bindings={str(Path('src/recursive_horizons')/(name+'.py')).replace('\\','/'):
            hashlib.sha256((ROOT/'src/recursive_horizons'/str(name+'.py')).read_bytes()).hexdigest()
            for name in SOURCE_NAMES})
    report['source_bindings']['scripts/run_nsc_atlas_fixed_onset.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    phase='input authentication';calls=accepted=0;maxima={};minimum_A=None;last_x=None
    history=hashlib.sha256();state=None
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
            if x<=front and (any(np.any(v!=0.) for v in value['cap'].values())
                             or np.any(value['departure']['c']!=0.)):
                raise ValueError('scalar zero-source recovery failed before the front')
            row=coupled_fixed_stage(value,basis=basis,x=x,eta=eta,source=source,gauge=ref,
                alpha=alpha,tube=dict(r=[.2,.8],L=[-2.,2.],c=[-4.,4.]),budget=budget)
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
                # Differential of the derived scalar applied to the time
                # residual, with the same reference rate on both paths.
                local=gram.wall;idx=local['indices'];coeff=dense['state']['cap']['chi']
                chi=float(local['value']@coeff[idx]);chiz=float(local['first']@coeff[idx]);f=chi/R
                k=math.sqrt(alpha/6);u=math.sinh(k*f)/k;root=math.cosh(k*f)
                radius_residual=float(dr['wall']['r'])
                scalar_residual=root*(float(local['value']@dr['cap']['chi'][idx])
                    +chiz*radius_residual/(1-R*R))+(u-f*root)*radius_residual
                maximum('time_derived.wall_c',scalar_residual)
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
            normalized_cap_mass=cap['mass'],cap_lapse=cap['lapse'],
            cap_L2={k:math.sqrt(gram.squared_norm(v,budget=budget)) for k,v in state['cap'].items()},
            exterior_departure_max={k:float(np.max(abs(v))) for k,v in state['departure'].items()},
            connection_residuals=final_row['connection']['residual'],
            source_constraint=final_row['source_constraint'],
            state_sha256=digest_state(state))
        if calls!=5*accepted: raise ValueError('common stage accounting mismatch')
        budget.check();report['status']='completed_numerical_interval'
    except (ValueError,RuntimeError,FloatingPointError,OverflowError,KeyError,TypeError) as error:
        report['status']='stopped';report['failure']=dict(phase=phase,reason=str(error))
    report.update(maximum_absolute=maxima,minimum_stage_metric_certificate=minimum_A,
        accepted_history_sha256=history.hexdigest(),
        execution=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
            common_stage_calls=calls,accepted_steps=accepted,last_accepted_x=last_x,
            source_integrations=0,retries=0,raw_history_saved=False))
    args.output.write_text(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False)+'\n',
                           encoding='utf-8',newline='\n')
    print(json.dumps(dict(status=report['status'],execution=report['execution'],
                         failure=report.get('failure'))),flush=True)
    return 0 if report['status']=='completed_numerical_interval' else 1


if __name__=='__main__': raise SystemExit(main())
