#!/usr/bin/env python3
"""One leading prepared-field endpoint and its existing common-stage ports."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from run_nsc_atlas_fixed_onset import ROOT,frozen_sources,plain,producer
from check_nsc_atlas_prepared_history import array_hash
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget,source_jet
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget,solve_prepared_field
from recursive_horizons.nsc_atlas_cap_field import CapLeadingField
from recursive_horizons.nsc_atlas_background_profile import build_prepared_background
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_wedge_field import WedgeField
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_fixed_cap_stage import quiet_fixed_state,coupled_fixed_stage

INPUT='results/summaries/nsc-atlas-fixed-history-onset-v1.json'
PREPARED='results/summaries/nsc-atlas-prepared-history-v2.json'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists(): raise ValueError('preserve the existing comparison')
    if not args.output.parent.is_dir(): raise ValueError('output parent must exist')
    commit,bindings=frozen_sources()
    raw=(ROOT/INPUT).read_bytes();saved=json.loads(raw)
    prepared_raw=(ROOT/PREPARED).read_bytes();reference=json.loads(prepared_raw)
    original=(ROOT/producer.BASE_INPUT).read_bytes()
    if hashlib.sha256(original).hexdigest()!=saved['input']['sha256']:
        raise ValueError('original source geometry changed')
    if saved['status']!='completed_numerical_interval': raise ValueError('require completed endpoint')
    geometry=json.loads(original)['geometry']
    budget=Budget(maximum_evaluations=10000000,cpu_seconds=60.)
    report=dict(schema='nsc-atlas-prepared-endpoint-connection-v1',source_commit=commit,
        source_bindings=bindings,input=dict(path=INPUT,sha256=hashlib.sha256(raw).hexdigest()),
        prepared_reference=dict(path=PREPARED,sha256=hashlib.sha256(prepared_raw).hexdigest()),
        original_input=dict(path=producer.BASE_INPUT,sha256=hashlib.sha256(original).hexdigest()),
        limits=dict(charged_values=10000000,cpu_seconds=60.),
        scope='one alpha-zero reference event from prepared leading fields and zero cap departures; not the evolved finite-alpha state')
    phase='reference endpoint';preparation=None;source_completed=False
    try:
        x=saved['geometry']['x_right']
        gauge=ReferenceWallGauge(geometry,order=64,budget=budget).sample(x)
        preparation=PreparedFieldBudget(cpu_seconds=60.,kernel_values_max=10000000-budget.evaluations)
        preparation.started=budget.started
        phase='single preparation and endpoint jets'
        try:
            field=solve_prepared_field(nodes=257,geometry_knots=257,budget=preparation)
            source_completed=True
            if any(field.domain[key]!=geometry[key] for key in ('U0','Ua','Ub','Va','Vb')):
                raise ValueError('prepared and coupled null domains differ')
            cap=CapLeadingField(field,derivative_order=4)
            transfer=build_prepared_background(cap);background=transfer['background']
            identity=array_hash(background.inputs)
            if not identity==reference['input_sha256']==saved['preparation']['input_sha256']:
                raise ValueError('prepared history identity differs')
            V=float(gauge['V']);Vx=V*(-x/(.8*(.8-x)))
            t=(V-field.domain['Va'])/(field.domain['Vb']-field.domain['Va'])
            R=float(gauge['r']);A0=1-R*R;offset=cap.domain['Ta']
            q=float(gauge['T'])+offset+float(np.arctanh(R))
            minus=float(gauge['T'])+offset-float(np.arctanh(R))
            if minus>background.inputs['labels'][background.onset['index']]:
                raise ValueError('reference endpoint no longer has a quiet reflected foot')
            q_direct=float(field.geometry['arrival'](V))
            qx=float(gauge['T_x']+gauge['r_x_total']/A0)
            direct=field.wall(np.array([t]))
            direct_profile=cap.profile(np.array([q_direct,q]))
            retained=background.profile.jet(np.array([q_direct,q]),budget=preparation)
            eta=np.linspace(0.,1.,65)
            U=field.domain['U0']+eta[1:-1]*float(gauge['W'])/(-3.2/geometry['U0'])
            exterior=WedgeField(field).bulk(U,t,order=64)
            p,px,_=source_jet(x,geometry);g=float(x*p);gx=float(p+x*px)
            report['preparation']=dict(input_sha256=identity,charged_values=preparation.evaluations,
                field_nodes=257,geometry_knots=257,clock_offset=offset)
            Jx=float(direct['wall_dchi_dV_over_epsilon'][0])*Vx
            report['clock']=dict(x=x,V=V,V_x=Vx,radius=R,reference_T=float(gauge['T']),
                history_label=q,prepared_label=q_direct,label_difference=q-q_direct,
                quiet_reflected_label=minus,
                history_label_rate=qx,prepared_label_rate=float(field.geometry['arrival'](V,1))*Vx,
                wall_U_difference=float(gauge['U'])-float(field.geometry['wall'](V)))
            report['field']=dict(prepared_wall_value=float(direct['wall_chi_over_epsilon'][0]),
                prepared_wall_material_rate=Jx,direct_profiles=direct_profile,retained_profiles=retained,
                derivative_conventions='profile derivatives use q; material rates use x; J*V_x uses the prepared zero chi_U wall condition',
                retained_material_rate=float(retained['F_prime'][1])*qx,
                retained_minus_prepared_rate=float(retained['F_prime'][1])*qx-Jx,
                source_value_identity=g-float(direct['source_chi_over_epsilon'][0]),
                interior_values_sha256=array_hash({'U':U,'chi':exterior['chi']}))
        finally:
            budget.check(preparation.evaluations)
        phase='one leading common stage'
        basis=FixedCapBasis([Q(v) for v in saved['geometry']['exact_cap_knots']],budget=budget)
        state=quiet_fixed_state(basis,exterior_cells=64)
        state['departure']['c']=exterior['chi']-g
        source=dict(r=x,log_C=0.,c=g,r_x=1.,log_C_x=0.,c_x=gx,r_xx=0.)
        row=coupled_fixed_stage(state,basis=basis,x=x,eta=eta,source=source,gauge=gauge,alpha=0.,
            tube=dict(r=[.2,.8],L=[-2.,2.],c=[-4.,4.]),budget=budget,
            background=background,background_clock_offset=offset)
        connection=row['connection'];caprow=connection['cap'];rates=connection['wall_rates']
        phase='endpoint comparison'
        gram=basis.gram(Q(float(np.arctanh(R))),budget=budget)
        trace=np.zeros(basis.dimension);trace[gram.wall['indices']]=gram.wall['value']
        port=-rates['T']*A0*R*caprow['wall']['incoming_discrepancy']*gram.solve(trace,budget=budget)
        density_error=caprow['rates']['canonical_density']-port
        report.update(status='completed_reference_comparison',common_wall=row['physical_wall'],
            wall_rates=rates,cap_wall=caprow['wall'],residual=connection['residual'],
            projected_gradient_discrepancy=caprow['projected_gradient_discrepancy'],
            prepared_minus_prefix_rate=Jx-connection['outgoing']['c_V'],
            prefix_outgoing=connection['outgoing'],
            prefix_derivative_coordinate='x (the common stage names this null coordinate V)',
            cap_c_rate=connection['cap_c_x'],hybrid_c_rate=rates['c'],
            cap_rate_minus_retained_rate=connection['cap_c_x']-report['field']['retained_material_rate'],
            cap_departure_rate_L2={key:float(np.sqrt(gram.squared_norm(value,budget=budget)))
                for key,value in caprow['rates'].items()},
            incoming_boundary_force_L2=float(np.sqrt(gram.squared_norm(port,budget=budget))),
            density_rate_minus_boundary_force_L2=float(np.sqrt(gram.squared_norm(density_error,budget=budget))),
            exterior_departure_max=float(np.max(abs(state['departure']['c']))),
            numerical_error_enclosure=False)
        budget.check()
    except (ValueError,RuntimeError,FloatingPointError,OverflowError,KeyError,TypeError) as error:
        report.update(status='stopped',failure=dict(phase=phase,reason=str(error)))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_integrations=int(source_completed),retries=0,time_steps=0,raw_history_saved=False)
    args.output.write_text(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False)+'\n',
                           encoding='utf-8',newline='\n')
    keys=('status','execution','failure','clock','field','residual','prepared_minus_prefix_rate','cap_rate_minus_retained_rate')
    print(json.dumps(plain({k:report.get(k) for k in keys})),flush=True)
    return 0 if report['status']=='completed_reference_comparison' else 1


if __name__=='__main__': raise SystemExit(main())
