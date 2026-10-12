#!/usr/bin/env python3
"""One time-box high-k moment control, with the full axial cell enclosed."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_local_prepared_response as P
from recursive_horizons.nsc_radius_transfer_sums import load_radius_transfer_inputs,evaluate_transfer_catalog
from recursive_horizons.nsc_ks_uniform_axial_geometry import uniform_axial_geometry_jets
from recursive_horizons.nsc_ks_high_radius_derivatives import restore_packed
from recursive_horizons.nsc_scaled_reference_projector import scaled_reference_projector
from recursive_horizons.nsc_global_projector_bounds import global_projector_bounds
from recursive_horizons.nsc_weyl_integrated_moments import integrated_defect_moments
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper
from recursive_horizons.nsc_ks_profile_identity import profile_identity

OUTPUT=ROOT/'results/development/nsc-weyl-integrated-moment-control.json'
RADIUS=ROOT/'results/development/nsc-ks-radius-covariance-bounds.json'
FINE=ROOT/'results/development/nsc-ks-fine-trajectory.json'
BITS=160
CAP=180.
SPLITS=(256,1024)
OWNED=('scripts/derive_nsc_weyl_integrated_moment_control.py','docs/nsc-weyl-integrated-moments.md',
       'src/recursive_horizons/nsc_weyl_integrated_moments.py','src/recursive_horizons/nsc_ks_uniform_axial_geometry.py',
       'src/recursive_horizons/nsc_global_projector_bounds.py','src/recursive_horizons/nsc_scaled_reference_projector.py',
       'src/recursive_horizons/nsc_radius_transfer_sums.py')


def digest(path):
    p=Path(path);return sha256((p if p.is_absolute() else ROOT/p).read_bytes()).hexdigest()


def compute():
    radius=json.loads(RADIUS.read_text());fine=json.loads(FINE.read_text())
    for record in (radius,fine):
        for category in ('source_hashes','input_hashes'):
            for path,expected in record.get(category,{}).items():
                if digest(path)!=expected:raise ValueError('bound input changed: '+path)
    if digest(fine['payload']['path'])!=fine['payload']['sha256']:raise ValueError('saved history changed')
    with np.load(ROOT/fine['payload']['path'],allow_pickle=False) as f:meta=json.loads(f['metadata_json'].tobytes())
    family=P.family(P.ALPHA)
    if profile_identity(family)!=radius['history_identity']:raise ValueError('history differs from radius certificate')
    with ctx.workprec(BITS):
        inputs=load_radius_transfer_inputs(root=ROOT,bits=BITS)
        W=[restore_packed(v) for v in radius['geometry']['W_fourier_l1']]
        geometry=uniform_axial_geometry_jets(family,arb((259,-8),(1,-20)),W,[arb(0)]*10,
            period_left=meta['period_origin'],period_length=Q(205,512),
            axial_profile_identity=radius['geometry']['axial_profile_identity'],bits=BITS)
        m,ell=meta['parameters']['mass'],meta['parameters']['angular']
        global_bounds=global_projector_bounds(mass=m,angular=ell,geometry=geometry,bits=BITS)
        rows=[]
        for K in SPLITS:
            midpoint=(arb(2)/K).man_exp()
            mu=arb(midpoint,midpoint)
            scaled={s:scaled_reference_projector(mu=mu,sign=s,mass=m,angular=ell,geometry=geometry,bits=BITS) for s in (-1,1)}
            catalog=evaluate_transfer_catalog(inputs,K,bits=BITS)
            for cut in (False,True):
                result=integrated_defect_moments(scaled,global_bounds,geometry.derivative('inv_a'),Q(205,512),catalog,K,
                    low_momentum_auxiliary=cut,uniform_axial_coverage=True,bits=BITS)
                result['canonical_momentum_split']=K
                result['time_box_integral_uppers']=[exact_upper(restored_upper(result[key])/2**19)
                    for key in ('zeroth_moment_density_upper','first_moment_density_upper')]
                rows.append(result)
        return {'schema':'NSC-WEYL-INTEGRATED-MOMENT-CONTROL-v1','accountable_author':'Douglas Ek',
            'status':'OPEN: high-k defect moments on one time box; no physical UV/gate certificate',
            'rows':rows,'history_identity':geometry.history_identity,'source_family_control':'14_1',
            'geometry_box':{'rho_center':'259/256','rho_radius':'1/1048576',
                'whole_axial_cell':True,'period_left_exact':str(Q(meta['period_origin'])),'period_length_exact':'205/512'},
            'global_P4_operator_upper':exact_upper(global_bounds.derivative_bound(4)),
            'global_P4_zz_operator_upper':exact_upper(global_bounds.derivative_bound(4,z=2)),
            'source_hashes':{p:digest(p) for p in OWNED},
            'input_hashes':{str(p.relative_to(ROOT)):digest(p) for p in (RADIUS,FINE,ROOT/fine['payload']['path'])},
            'scope':{'physical_local_gate':'OPEN','whole_time_domain_covered':False,
                'physical_subtraction_changed':False,'source_energy_tail_bound':None,
                'initial_covariance_error_bound':None,'low_canonical_momentum_residual_bound':None,
                'source_value_error_bound':None,'new_field_or_source_runs':0,'metric_timestep':False},
            'reproducer':'python scripts/derive_nsc_weyl_integrated_moment_control.py --check'}


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--record',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args()
    def stop(*_):raise TimeoutError('two-split coefficient control exceeded180 CPU seconds')
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CAP)
    start=time.process_time()
    try:result=compute()
    finally:signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    cpu=time.process_time()-start
    if args.record:
        if OUTPUT.exists():raise FileExistsError('integrated coefficient record exists; use --check')
        result['runtime']={'CPU_seconds':cpu,'CPU_cap':CAP}
        OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    else:
        old=json.loads(OUTPUT.read_text());old.pop('runtime')
        if result!=old:raise ValueError('integrated moment replay differs')
    with ctx.workprec(BITS):display=[{'K':r['canonical_momentum_split'],'auxiliary':r['auxiliary'],
        'moment_densities':[float(restored_upper(r[k])) for k in ('zeroth_moment_density_upper','first_moment_density_upper')],
        'one_time_box_integrals':[float(restored_upper(v)) for v in r['time_box_integral_uppers']]} for r in result['rows']]
    print(json.dumps({'status':result['status'],'comparisons':display,'CPU_seconds':cpu},indent=2))
