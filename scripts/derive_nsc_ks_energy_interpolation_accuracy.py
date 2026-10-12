#!/usr/bin/env python3
"""Bound the interpolation component of the completed retained-source control."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np
from flint import arb,ctx,fmpq
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'));sys.path.insert(0,str(ROOT/'scripts'))
import derive_nsc_ks_retained_response_sum as S
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper,complex_ball
from recursive_horizons.nsc_ks_energy_node_bound import actual_node_remainder_factor,multiply_derivative_bound
from recursive_horizons.nsc_ks_energy_interpolation_bound import energy_derivative_majorants,weighted_A_up_frobenius,field_remainders_from_basis
from recursive_horizons.nsc_ks_finite_matter_error import finite_matter_error

OUTPUT=ROOT/'results/development/nsc-ks-energy-interpolation-accuracy.json'
RADIUS=ROOT/'results/development/nsc-ks-radius-bound.json'
BITS=160
OWNED=('scripts/derive_nsc_ks_energy_interpolation_accuracy.py','docs/nsc-ks-energy-interpolation-accuracy.md',
       'src/recursive_horizons/nsc_ks_energy_node_bound.py','src/recursive_horizons/nsc_ks_energy_interpolation_bound.py',
       'src/recursive_horizons/nsc_ks_finite_matter_error.py')


def rational(value):
    q=Fraction(value);return arb(fmpq(q.numerator,q.denominator))


def gamma_from_blocks(blocks):
    row_max,column_max=arb(0),arb(0)
    for block in blocks:
        for row in block:
            row_max=arb.max(row_max,sum((complex_ball(v).abs_upper() for v in row),arb(0)))
        for column in block.T:
            column_max=arb.max(column_max,sum((complex_ball(v).abs_upper() for v in column),arb(0)))
    # Retain any stored roundoff asymmetry instead of assuming exact Hermiticity.
    return (row_max*column_max).sqrt().upper()


def integrals(radius,mass,angular,rho_up,normal_support):
    b=radius['bounds'];get=lambda name:rational(b[name]['exact_rational'])
    if not arb(1)<arb(rho_up)<=get('rho_upper'):raise ValueError('preparation outside the proven geometry slab')
    alpha=arb(float(radius['control']['history_amplitude']))
    if not alpha>0 or any(Fraction(v['exact_rational'])!=0 for v in b['U_sup_orders_0_1_2']):
        raise ValueError('this component record is bound to the single nonzero w-control')
    W,Wz=(rational(b['w_sup_orders_0_1_2'][j]['exact_rational']) for j in (0,1))
    a,r=get('axial_lower'),get('radius_lower')
    duration=arb(rho_up)-1;sigma=arb(normal_support);ell=abs(arb(angular))
    if not arb(0)<sigma<=arb(.03):raise ValueError('normal support differs from the geometry bound')
    common=duration*ell/a
    D=duration/a**2;K0=duration/a*(abs(arb(mass))+ell/r)
    K1=common*sigma*Wz/r**2
    J0=common*sigma*W/(alpha*r**2)
    J1=common*sigma*Wz/(alpha*r**2)*(1+2*sigma*W/r)
    return tuple(v.upper() for v in (K0,D,K1,J0,J1)),a,r


def compute():
    if not S.OUTPUT.exists():raise ValueError('complete finite inventory sum required')
    summed=json.loads(S.OUTPUT.read_text())
    if not summed['all_retained_nonzero_angular_families_included']:raise ValueError('partial inventory')
    op,meta,coeff,baseline,grouped,signature=S.context()
    if summed['source_hashes']!=signature:raise ValueError('source sum dependencies changed')
    radius=json.loads(RADIUS.read_text())
    for mapping in ('source_hashes','input_hashes'):
        for p,h in radius[mapping].items():
            if S.digest(p)!=h:raise ValueError('radius-bound owner changed')
    results=[];total=[arb(0),arb(0)]
    with ctx.workprec(BITS):
        for key in sorted(grouped):
            record,_,_=S.read_family(key,op,meta,coeff,grouped[key],signature)
            with np.load(ROOT/record['payload']['path'],allow_pickle=False) as f:a={k:f[k] for k in f.files}
            mass,ell=grouped[key][0].mass,grouped[key][0].angular
            bounds,amin,rmin=integrals(radius,mass,ell,op.rho_up,max(op.binding.outer_radii))
            node=actual_node_remainder_factor(a['operator/energies'],record['operator_energy_interval'],bits=BITS)
            majorants=energy_derivative_majorants(*bounds,len(a['operator/energies']),bits=BITS)
            errors={name:multiply_derivative_bound(node,majorants['dE_n_'+name+'_upper'],bits=BITS)
                    for name in ('W','Wz','Y','Yz')}
            norm2=arb(0);positive_gamma=arb(0);negative_gamma=arb(0);max_energy=arb(0)
            for panel in grouped[key]:
                positive_gamma=arb.max(positive_gamma,gamma_from_blocks(panel.covariance))
                negative_gamma=arb.max(negative_gamma,gamma_from_blocks(panel.negative_partner(meta['config']).covariance))
                max_energy=arb.max(max_energy,arb(float(np.max(panel.energies))))
                for part in panel.batches(S.P.BATCH_ROWS):
                    # Same repeated physical fiber weights as FixedSourcePreparation.
                    src=S.FixedSourcePreparation.from_signed_blocks(part.energies,part.weights,part.covariance)
                    norm=weighted_A_up_frobenius(a['upstream/'+part.name],src.column_weights,bits=BITS)
                    norm2+=norm**2
            norm=norm2.sqrt().upper()
            field=field_remainders_from_basis(errors['W'],errors['Wz'],norm,max_energy,
                Y_remainder=errors['Y'],Yz_remainder=errors['Yz'],bits=BITS)
            growth=bounds[0].exp();Fnorm=arb(2).sqrt()*growth*norm
            Fznorm=arb(2).sqrt()*growth*(bounds[2]+max_energy)*norm
            channel=meta['channels'][key[0]];mu=arb(channel['copy_count']*channel['degeneracy'])/2
            contributions=[]
            for gamma in (positive_gamma,negative_gamma):
                err=finite_matter_error(Fnorm,Fznorm,restored_upper(field['weighted_F_remainder_upper']),
                    restored_upper(field['weighted_F_z_remainder_upper']),gamma,mass=mass,absolute_angular=abs(ell),
                    axial_lower=amin,radius_lower=rmin,multiplicity=mu,bits=BITS)
                contributions.append(err)
                for j,name in enumerate(('N','beta')):total[j]+=restored_upper(err[name])
            results.append({'positive_family':list(key),'operator_digest':record['operator_digest'],
                'integral_bounds':dict(zip(('K0','D','K1','J0','J1'),map(exact_upper,bounds))),
                'actual_node_remainder_factor':node,'row_operator_remainder':errors,
                'weighted_A_up_norm_upper':exact_upper(norm),'field_interpolation_remainder':field,
                'separate_signed_matter_error_uppers':contributions})
        return {'schema':'NSC-KS-ENERGY-INTERPOLATION-ACCURACY-v1','accountable_author':'Douglas Ek',
            'status':'PASS: interpolation remainder component; node evolution and physical gate OPEN',
            'families':results,'total_N_beta_interpolation_error_upper':[exact_upper(v) for v in total],
            'scope':{'actual_binary_nodes_used':True,'unchanged_original_source':True,
                'node_evolution_error_bound':None,'interpolation_arithmetic_error_bound':None,
                'source_preparation_error_bound':None,'source_quadrature_error_bound':None,
                'changed_history_UV_bound':None,'physical_local_gate':'OPEN','constraint_root_claimed':False,
                'metric_timestep':False,'new_field_or_source_solves':0},
            'source_hashes':{p:S.digest(p) for p in OWNED},
            'input_hashes':{str(S.OUTPUT.relative_to(ROOT)):S.digest(S.OUTPUT),str(RADIUS.relative_to(ROOT)):S.digest(RADIUS)},
            'reproducer':'python scripts/derive_nsc_ks_energy_interpolation_accuracy.py --check'}


if __name__=='__main__':
    p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group(required=True)
    g.add_argument('--run',action='store_true');g.add_argument('--check',action='store_true')
    args=p.parse_args()
    if args.run:
        if OUTPUT.exists():raise FileExistsError('component record exists; use --check')
        r=compute();OUTPUT.write_text(json.dumps(r,sort_keys=True,indent=2,allow_nan=False)+'\n')
    else:
        r=compute()
        if r!=json.loads(OUTPUT.read_text()):raise ValueError('interpolation component replay differs')
    with ctx.workprec(BITS):values=[float(restored_upper(v)) for v in r['total_N_beta_interpolation_error_upper']]
    print(json.dumps({'status':r['status'],'families':len(r['families']),
        'N_beta_interpolation_component_upper':values},indent=2))
