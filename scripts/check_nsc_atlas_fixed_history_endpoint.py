#!/usr/bin/env python3
"""One prepared-history zero-coupling velocity control; no time evolution."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from run_nsc_atlas_fixed_onset import ROOT,frozen_sources,prepare_history,plain
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_reference import fixed_cap_reference

INPUT='results/summaries/nsc-atlas-fixed-history-onset-v1.json'
PREPARED='results/summaries/nsc-atlas-prepared-history-v2.json'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    if args.output.exists(): raise ValueError('preserve the existing diagnostic')
    if not args.output.parent.is_dir(): raise ValueError('output parent must exist')
    commit,bindings=frozen_sources()
    raw=(ROOT/INPUT).read_bytes();saved=json.loads(raw)
    prepared_raw=(ROOT/PREPARED).read_bytes();reference=json.loads(prepared_raw)
    if saved['status']!='completed_numerical_interval': raise ValueError('require completed reference interval')
    budget=Budget(maximum_evaluations=10000000,cpu_seconds=60.)
    report=dict(schema='nsc-atlas-fixed-history-endpoint-velocity-v1',source_commit=commit,
        source_bindings=bindings,input=dict(path=INPUT,sha256=hashlib.sha256(raw).hexdigest()),
        prepared_reference=dict(path=PREPARED,sha256=hashlib.sha256(prepared_raw).hexdigest()),
        limits=dict(charged_values=10000000,cpu_seconds=60.),preparation={},
        scope='one zero-coupling, zero-departure limit control at the recorded physical-history endpoint; not the finite-alpha endpoint state')
    phase='history preparation'
    try:
        background,_=prepare_history(budget,report['preparation'])
        if not (report['preparation']['input_sha256']==reference['input_sha256']
                ==saved['preparation']['input_sha256']):
            raise ValueError('prepared history identity differs')
        R=saved['final']['wall']['r'];T=saved['final']['background']['time']
        basis=FixedCapBasis([Q(v) for v in saved['geometry']['exact_cap_knots']],budget=budget)
        Z=Q(float(np.arctanh(R)));zero=np.zeros(basis.dimension)
        report['point']=dict(radius=R,background_time=T,alpha=0.,zero_departures=True,
            T_x=1.,R_x=0.,incoming_w_plus=0.,exact_cap_knots=saved['geometry']['exact_cap_knots'])
        phase='one radial reference'
        cap=fixed_cap_reference(basis,zero,zero,Z,alpha=0.,T_x=1.,radius_rate=0.,incoming_w_plus=0.,
                                background=background,background_time=T,budget=budget)
        phase='independent projected wall contractions'
        gram=basis.gram(Z,budget=budget);local=gram.wall;idx=local['indices'];e=local['value']
        a=cap['background']['projected_time_chi'];mu=np.exp(cap['lapse'])
        contractions={}
        for name,key in (('forward','gradient'),('adjoint','adjoint_gradient')):
            residual=2*mu*gram.solve(cap[key]['canonical_density'],budget=budget)-a
            contractions[name]=dict(wall=float(e@residual[idx]),
                L2=float(np.sqrt(gram.squared_norm(residual,budget=budget))))
        report.update(status='completed_limit_diagnostic',expected_velocity_discrepancy=0.,
            measured_velocity_discrepancy=cap['wall']['velocity_discrepancy'],
            projected_contractions=contractions,
            contraction_identity_defect=contractions['forward']['wall']-cap['wall']['velocity_discrepancy'],
            projected_gradient_discrepancy=cap['projected_gradient_discrepancy'],
            mass=cap['mass'],lapse=cap['lapse'],wall=cap['wall'],
            radial_intervals=cap['background']['radial_intervals'],
            scaled_gram_condition=cap['scaled_gram_condition'],numerical_error_enclosure=False)
        budget.check()
    except (ValueError,RuntimeError,FloatingPointError,OverflowError,KeyError,TypeError) as error:
        report.update(status='stopped',failure=dict(phase=phase,reason=str(error)))
    report['execution']=dict(charged_values=budget.evaluations,cpu_seconds=time.process_time()-budget.started,
        source_integrations=report['preparation'].get('source_integrations_completed',0),
        retries=0,time_steps=0,raw_history_saved=False)
    args.output.write_text(json.dumps(plain(report),indent=2,sort_keys=True,allow_nan=False)+'\n',
                           encoding='utf-8',newline='\n')
    print(json.dumps(plain({k:report.get(k) for k in ('status','execution','failure',
        'measured_velocity_discrepancy','projected_contractions','contraction_identity_defect')})),flush=True)
    return 0 if report['status']=='completed_limit_diagnostic' else 1


if __name__=='__main__': raise SystemExit(main())
