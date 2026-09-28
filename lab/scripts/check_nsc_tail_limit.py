#!/usr/bin/env python3
"""Verify tail convergence criteria without inserting a cosmological fraction."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_tail_limit import motif_hoppings,coefficients,finite_response,weyl_disk,record_family
from recursive_horizons.nsc_regulated import recursive_response,direct_chain
from check_nsc_scale_closure import compare
OUTPUT=ROOT/'results/nsc-5-tail-limit.json'


def calculate():
    n=sp.symbols('n',positive=True,integer=True)
    gamma=(n+1)/n
    assert sp.simplify(gamma.subs(n,n+1)-(2-1/gamma))==0
    critical=[{'depth':k,'inverse_response':float(gamma.subs(n,k)),'error_to_limit':1/k}
              for k in (1,2,4,8,16,32,64)]
    families=[record_family(omega) for omega in (1.3,4.,50.)]
    assert families[0]['endpoint_class']==families[1]['endpoint_class']=='limit_point'
    assert families[2]['endpoint_class']=='limit_circle'
    assert families[0]['rows'][-1]['log10_diameter'] < -150
    assert families[2]['rows'][-1]['diameter_numeric']>.005

    motif,base=motif_hoppings(16)
    matrix_controls=[]
    for omega in (1.3,4.,50.):
        for depth in (1,2,4):
            with mp.workdps(60):
                projected=complex(finite_response(coefficients(base,omega,depth),mp.mpc('.3','.4')))
            block=recursive_response(motif['H'],motif['B'],.3+.4j,omega,depth)
            direct=direct_chain(motif['H'],motif['B'],.3+.4j,omega,depth)
            from_block=np.linalg.inv(block)[0,0];from_direct=np.linalg.inv(direct)[0,0]
            errors=[abs(projected-from_block),abs(projected-from_direct)]
            assert max(errors)<1e-12
            matrix_controls.append({'Omega':omega,'depth':depth,'block_error':float(errors[0]),'direct_error':float(errors[1])})

    # Independent direct small-matrix circle check and Gram identity.
    small=[]
    with mp.workdps(80):
        for a_raw in ([1.],[.7,1.2],[1.1,.8,1.4]):
            a=[mp.mpf(v) for v in a_raw];z=mp.mpc('.6','.8')
            center,radius,(p,s,c)=weyl_disk(a,z)
            gram=abs(p*s-abs(c)**2-mp.im(c)/mp.im(z))
            assert gram<mp.mpf('1e-60')
            mat=np.diag(a_raw,1)+np.diag(a_raw,-1)
            direct=np.linalg.inv(complex(z)*np.eye(len(mat))-mat)[0,0]
            error=abs(direct-complex(finite_response(a,z)))
            assert error<1e-14
            small.append({'couplings':a_raw,'gram_residual':float(gram),'direct_error':float(error),
                          'radius':float(radius)})
    resolution=[]
    for points in (16,32,64,128):
        m,_=motif_hoppings(points)
        resolution.append({'intervals':points,'q_lattice':m['zero_energy_transfer'],
                           'q_continuum':m['continuum_zero_transfer'],
                           'q_error':abs(m['zero_energy_transfer']-m['continuum_zero_transfer']),
                           'critical_Omega_lattice':1/m['zero_energy_transfer']})
    high=record_family(1.3,depths=(16,),precision=260)['rows'][0]
    precision_difference=abs(high['log10_diameter']-families[0]['rows'][4]['log10_diameter'])
    assert precision_difference<1e-12
    return {'schema':'nsc-tail-limit-v1','artifact_id':'NSC-5-TAIL-LIMIT',
            'classification':'exact_endpoint_and_Weyl_disk_tail_criteria_for_the_declared_smooth_spatial_chain_not_cosmological_digits',
            'source_hashes':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in (
                'src/recursive_horizons/nsc_tail_limit.py','scripts/check_nsc_tail_limit.py',
                'src/recursive_horizons/nsc_smooth_geometry.py','src/recursive_horizons/nsc_regulated.py')},
            'conventions':{'geometry':'smooth_motif(a=1,R=2,kappa=1), 16 intervals per motif',
                           'Hilbert_space':'unweighted ell² of the interleaved node/edge chain',
                           'scaling':'H_n=Omega^n H; B_n=Omega^n B; one common dimensional energy',
                           'resolvent':'G00(z)=<e0,(z-J)^-1 e0>, z=0.3+0.4i in parent units',
                           'Jacobi_sign_gauge':'unitary sign change of the interleaved path; first-site projection unchanged'},
            'endpoint_theorem':{'q':'product_i[(1/h-w_i/2)/(1/h+w_i/2)]',
                                'two_zero_energy_motif_multipliers':'(-1)^m q and (-1)^m/(Omega*q)',
                                'limit_point_condition':'Omega*q<=1 (q<1); equality is nondecaying odd solution',
                                'limit_circle_condition':'Omega*q>1; both independent endpoint solutions are square summable',
                                'left_boundary_caution':'odd formal solution need not satisfy regular left boundary but is required to classify infinity',
                                'scope':'fixed lattice and stipulated scaling/Hilbert measure; not a continuum or full NSC theorem'},
            'tail_bound':{'polynomials':'p_-1=0,p0=1,p_(j+1)=(z*p_j-a_(j-1)*p_(j-1))/a_j',
                          'radius':'R_M=1/(2 Im(z) sum_{j=0}^{M-1}|p_j|²)',
                          'consequence':'any two self-adjoint completions, including a Dirichlet truncation comparison, differ in G00 by at most 2R_M',
                          'passive_convention':'terminal self-energy Im Sigma<=0',
                          'mathematical_bound_exact':True,'displayed_values_interval_certified':False,
                          'bound_is_for_full_matrix_or_stress_derivatives':False,
                          'precision_difference_log10_diameter':precision_difference},
            'families':families,'matrix_and_scalar_controls':matrix_controls,'small_Gram_controls':small,
            'spatial_resolution':resolution,
            'critical_scalar_counterexample':{'Omega':4,'b':2,'K':2,'recursion':'Gamma_(n+1)=2-1/Gamma_n, Gamma_1=2',
                                             'exact_solution':'Gamma_n=(n+1)/n','limit':1,
                                             'successive_correction':'-1/[n(n+1)]','error':'1/n, not O(4^-n)',
                                             'controls':critical},
            'decimal_counterexamples':{'infinite_positive_geometric_sum':'sum_(n>=1)2^-n=1 has terminating decimal',
                                       'nonterminating_without_fractal_geometry':'1/3=0.333... in base ten and 0.1 in base three',
                                       'decimal_expansion_is_a_physical_signature':False},
            'observable_requirement':'derive the stress/observation functional and its error amplification before converting a resolvent bound into a dark fraction or detectability threshold',
            'nonclaims':{'physical_Omega_selected':False,'dark_fraction_95_percent_derived':False,
                         'positive_dark_density_increment_inferred':False,'cosmological_time_attractor_proved':False,
                         'precision_of_physical_constants_exceeds_lattice_accuracy':False,'new_to_world_mathematics_claimed':False},
            'primary_source':'https://www.mat.univie.ac.at/~gerald/ftp/book-jac/jacop.pdf',
            'terminal':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError('refusing to overwrite result')
    record=calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()),record);print('Tail endpoint, convergence and precision requirements reproduced; every field checked.')
    elif args.output:
        with args.output.open('x') as stream:json.dump(record,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
        print(args.output)
    else:print(json.dumps(record,indent=2,allow_nan=False))


if __name__=='__main__':main()
