#!/usr/bin/env python3
"""One bounded actual-history coefficient box, not a whole-history UV certificate."""
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
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_spacetime_geometry_jets import spacetime_geometry_jets
from recursive_horizons.nsc_scaled_reference_projector import scaled_reference_projector
from recursive_horizons.nsc_weyl_small_transfer_bound import small_transfer_coefficients,evaluate_coefficient_row
from recursive_horizons.nsc_ks_high_radius_derivatives import restore_packed
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper

OUTPUT=ROOT/'results/development/nsc-weyl-small-transfer-control.json'
RADIUS=ROOT/'results/development/nsc-ks-radius-covariance-bounds.json'
FINE=ROOT/'results/development/nsc-ks-fine-trajectory.json'
BITS=160
CAP=90.
K=256
OWNED=('scripts/derive_nsc_weyl_small_transfer_control.py','docs/nsc-weyl-small-transfer-bound.md',
       'src/recursive_horizons/nsc_weyl_small_transfer_bound.py','src/recursive_horizons/nsc_scaled_reference_projector.py',
       'src/recursive_horizons/nsc_ks_spacetime_geometry_jets.py','src/recursive_horizons/nsc_weyl_commutator_remainder.py')


def digest(path):
    p=Path(path);return sha256((p if p.is_absolute() else ROOT/p).read_bytes()).hexdigest()


def compute():
    radius=json.loads(RADIUS.read_text());fine=json.loads(FINE.read_text())
    for record in (radius,fine):
        for category in ('source_hashes','input_hashes'):
            for path,expected in record.get(category,{}).items():
                if digest(path)!=expected:raise ValueError('coefficient input changed: '+path)
    if digest(fine['payload']['path'])!=fine['payload']['sha256']:raise ValueError('history payload changed')
    with np.load(ROOT/fine['payload']['path'],allow_pickle=False) as f:meta=json.loads(f['metadata_json'].tobytes())
    family=P.family(P.ALPHA)
    if profile_identity(family)!=meta['history_identity'] or radius['history_identity']!=meta['history_identity']:
        raise ValueError('coefficient box must use the same saved history')
    center=family.directions[0].w.center
    with ctx.workprec(BITS):
        rho=arb((259,-8),(1,-20))
        z=arb(arb(float(center))+arb(3)/64,arb((1,-20)))
        mu=arb((1,-8),(1,-8))
        geometry=spacetime_geometry_jets(family,rho,z,order=9,bits=BITS)
        q=[restore_packed(v) for v in radius['geometry']['q_difference_fourier_l1']]
        rows=[]
        for sign in (-1,1):
            jets=scaled_reference_projector(mu=mu,sign=sign,mass=meta['parameters']['mass'],
                angular=meta['parameters']['angular'],geometry=geometry,bits=BITS)
            row=small_transfer_coefficients(jets,geometry.derivative('inv_a'),q,K,axial_order=2,bits=BITS)
            row['at_split_rho_upper']=[evaluate_coefficient_row(r,K,bits=BITS) for r in row['rows']]
            rows.append(row)
        result={'schema':'NSC-WEYL-SMALL-TRANSFER-CONTROL-v1','accountable_author':'Douglas Ek',
            'status':'PASS: one finite-history coefficient box; full UV/local gate OPEN',
            'canonical_momentum_cutoff':K,'rows_by_sign':rows,
            'box':{'rho_center':'259/256','rho_radius':'1/1048576',
                'z_center_exact':str(Q(float(center))+Q(3,64)),'z_radius_nominal':'1/1048576',
                'z_ball_display':str(z),'mu_interval_exact':['0','1/128']},
            'geometry_history_identity':geometry.history_identity,
            'source_family_control':'14_1','parameters':meta['parameters'],
            'scope':{'physical_local_gate':'OPEN','whole_preparation_slab_covered':False,
                'whole_axial_support_covered':False,'large_transfer_bound':None,
                'time_integral_bound':None,'initial_covariance_error_bound':None,
                'physical_source_energy_tail_bound':None,'new_field_or_source_runs':0,
                'metric_timestep':False},
            'source_hashes':{p:digest(p) for p in OWNED},
            'input_hashes':{str(p.relative_to(ROOT)):digest(p) for p in (RADIUS,FINE,ROOT/fine['payload']['path'])},
            'reproducer':'python scripts/derive_nsc_weyl_small_transfer_control.py --check'}
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--record',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args()
    def stop(*_):raise TimeoutError('coefficient control exceeded90 CPU seconds')
    previous=signal.signal(signal.SIGPROF,stop);signal.setitimer(signal.ITIMER_PROF,CAP)
    start=time.process_time()
    try:result=compute()
    finally:signal.setitimer(signal.ITIMER_PROF,0);signal.signal(signal.SIGPROF,previous)
    cpu=time.process_time()-start
    if args.record:
        if OUTPUT.exists():raise FileExistsError('coefficient record exists; use --check')
        result['runtime']={'CPU_seconds':cpu,'CPU_cap':CAP}
        OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    else:
        saved=json.loads(OUTPUT.read_text());saved.pop('runtime')
        if result!=saved:raise ValueError('coefficient-box replay differs')
    with ctx.workprec(BITS):display=[[float(restored_upper(v)) for v in r['at_split_rho_upper']] for r in result['rows_by_sign']]
    print(json.dumps({'status':result['status'],'rho_defect_bounds_at_k256_by_sign':display,'CPU_seconds':cpu},indent=2))
