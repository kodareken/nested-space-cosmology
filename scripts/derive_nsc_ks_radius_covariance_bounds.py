#!/usr/bin/env python3
"""Bind the true reciprocal Fourier moments to covariance transport constants."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb,ctx
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_ks_high_radius_derivatives import bound_pure_radius_difference,restore_packed
from recursive_horizons.nsc_fourier_covariance_moments import pure_radius_potential_integrals
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper

OUTPUT=ROOT/'results/development/nsc-ks-radius-covariance-bounds.json'
OWNED=('scripts/derive_nsc_ks_radius_covariance_bounds.py','docs/nsc-ks-radius-covariance-bounds.md',
       'src/recursive_horizons/nsc_ks_high_radius_derivatives.py','src/recursive_horizons/nsc_fourier_covariance_moments.py',
       'src/recursive_horizons/nsc_bloch_trace_bound.py','src/recursive_horizons/nsc_ks_profile_fourier_bound.py')
RECORDS=('results/development/nsc-ks-fine-profile-table.json','results/development/nsc-ks-radius-bound.json',
         'results/development/nsc-ks-source-inventory.json','results/development/nsc-ks-fine-trajectory.json')
BITS=160


def digest(path):
    p=Path(path);return sha256((p if p.is_absolute() else ROOT/p).read_bytes()).hexdigest()


def require_normal_support(geometry,history):
    sigma=Q(geometry['normal_support']['exact_rational'])
    directions=history['history_description']['directions']
    if not directions:raise ValueError('explicit normal windows required')
    for direction in directions:
        inner,outer=(Q(float.fromhex(direction[k])) for k in ('normal_inner','normal_outer'))
        if not 0<inner<outer<=sigma:raise ValueError('history exceeds the certified normal support')


def compute():
    records=[json.loads((ROOT/p).read_text()) for p in RECORDS]
    for record in records:
        for category in ('source_hashes','input_hashes'):
            for p,h in record.get(category,{}).items():
                if digest(p)!=h:raise ValueError('input dependency changed: '+p)
        if 'payload' in record:
            if digest(record['payload']['path'])!=record['payload']['sha256']:raise ValueError('input payload changed')
    with np.load(ROOT/records[2]['payload']['path'],allow_pickle=False) as f:source=json.loads(f['metadata_json'].tobytes())
    with np.load(ROOT/records[3]['payload']['path'],allow_pickle=False) as f:history=json.loads(f['metadata_json'].tobytes())
    with ctx.workprec(BITS):
        geometry=bound_pure_radius_difference(root=ROOT,max_derivative=9,bits=BITS)
        if geometry['axial_profile_identity']!=history['axial_profile_identity']:raise ValueError('radius and evolved profile identities differ')
        require_normal_support(geometry,history)
        if history['rho_up']>float(Q(records[1]['bounds']['rho_upper']['exact_rational'])):raise ValueError('history outside geometry slab')
        q0,q1=(restore_packed(geometry['q_difference_fourier_l1'][j]) for j in (0,1))
        rows=[];growth_max=arb(0);shift_max=arb(0)
        for channel in source['channels']:
            potential=pure_radius_potential_integrals(q0,q1,absolute_angular=channel['angular_eigenvalue'],
                axial_lower=Q(records[1]['bounds']['axial_lower']['exact_rational']),rho_up=history['rho_up'],bits=BITS)
            K0=restored_upper(potential['potential_zero_integral_upper'])
            K1=restored_upper(potential['potential_first_integral_upper'])
            growth=(2*K0).exp();shift=growth*K1
            growth_max=arb.max(growth_max,growth);shift_max=arb.max(shift_max,shift)
            rows.append({'group':channel['index'],'angular_abs':channel['angular_eigenvalue'],
                'compact_mass':channel['compact_mass'],'actual_angular_signs':1 if channel['angular_eigenvalue']==0 else 2,
                'potential_integrals':potential,'moment_comparison_diagonal_upper':exact_upper(growth),
                'zeroth_to_first_moment_upper':exact_upper(shift),
                'initial_moments':None,'true_defect_moments':None})
        return {'schema':'NSC-KS-RADIUS-COVARIANCE-BOUNDS-v1','accountable_author':'Douglas Ek',
            'status':'PASS: nonlinear spatial geometry and covariance-growth components; UV/local gate OPEN',
            'geometry':geometry,'groups':rows,'actual_angular_families':sum(r['actual_angular_signs'] for r in rows),
            'maximum_moment_diagonal_upper':exact_upper(growth_max),'maximum_moment_off_diagonal_upper':exact_upper(shift_max),
            'history_identity':history['history_identity'],'source_inventory_unchanged':True,
            'source_hashes':{p:digest(p) for p in OWNED},'input_hashes':{p:digest(p) for p in RECORDS},
            'scope':{'physical_local_gate':'OPEN','UV_error_bound':None,'source_error_bound':None,
                'actual_projector_derivative_bound':None,'time_integrated_defect_bound':None,
                'initial_C_up_minus_P_up_moments':None,'new_state_or_action':False,
                'new_field_or_source_runs':0,'metric_timestep':False},
            'reproducer':'python scripts/derive_nsc_ks_radius_covariance_bounds.py --check'}


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--record',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args();result=compute()
    if args.record:
        if OUTPUT.exists():raise FileExistsError('radius covariance record exists; use --check')
        OUTPUT.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')
    elif result!=json.loads(OUTPUT.read_text()):raise ValueError('radius covariance bound replay differs')
    with ctx.workprec(BITS):
        display={'status':result['status'],'families':result['actual_angular_families'],
            'q_Fourier_moments_0_1':[float(restore_packed(v)) for v in result['geometry']['q_difference_fourier_l1'][:2]],
            'moment_diagonal_upper':float(restored_upper(result['maximum_moment_diagonal_upper'])),
            'moment_off_diagonal_upper':float(restored_upper(result['maximum_moment_off_diagonal_upper']))}
    print(json.dumps(display,indent=2))
