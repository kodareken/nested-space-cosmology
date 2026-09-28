#!/usr/bin/env python3
"""Compute the repeated-throat gap and geometric-link finite recursion."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_geometric_chain import geometric_motif,bloch_gap,continuum_first_band_edge,continuum_discriminant
from recursive_horizons.nsc_regulated import recursive_response,direct_chain
from check_nsc_scale_closure import compare
OUTPUT=ROOT/'results/nsc-3-geometric-chain.json'


def calculate():
    families=[]
    for radius in (2.,4.,8.):
        continuum=continuum_first_band_edge(radius)
        lattice=[bloch_gap(radius,n) for n in (32,64,128)]
        errors=[abs(v['minimum_sampled_bloch_gap']-continuum['band_edge']) for v in lattice]
        orders=[float(np.log2(errors[i]/errors[i+1])) for i in range(2)]
        zero_trace,zero_det=continuum_discriminant(0.,radius)
        expected_trace=2*np.cosh(2*np.arcsinh(radius))
        assert abs(zero_trace-expected_trace)<1e-8 and abs(zero_det-1)<1e-8
        assert min(orders)>1.8 and max(orders)<2.3 and errors[-1]<.002
        families.append({'radius':radius,'continuum':continuum,'lattice':lattice,
                         'absolute_errors':errors,'observed_orders':orders,
                         'zero_discriminant':zero_trace,'analytic_zero_discriminant':float(expected_trace)})
    motif=geometric_motif(2.,8)
    free=geometric_motif(2.,32,kappa=0.)
    free_zero_gap=float(np.min(np.abs(np.linalg.eigvalsh(free['H']+free['B']+free['B'].T))))
    assert free_zero_gap < 1e-12
    controls=[]
    for omega in (1.,1.3):
        previous=None
        for depth in (2,4,8,16):
            reduced=recursive_response(motif['H'],motif['B'],.3+.4j,omega,depth)
            direct=direct_chain(motif['H'],motif['B'],.3+.4j,omega,depth)
            residual=float(np.linalg.norm(reduced-direct))
            assert residual<1e-10
            controls.append({'Omega':omega,'depth':depth,'Schur_residual':residual,
                             'response_00':[float(reduced[0,0].real),float(reduced[0,0].imag)],
                             'depth_change':None if previous is None else float(np.linalg.norm(reduced-previous))})
            previous=reduced
    return {
        'schema':'nsc-geometric-chain-v1','artifact_id':'NSC-3-GEOMETRIC-CHAIN',
        'classification':'periodic_repetition_of_the_spatial_throat_generates_a_band_gap_without_an_inserted_mass_and_geometric_links_close_finite_recursion',
        'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
            'src/recursive_horizons/nsc_geometric_chain.py','src/recursive_horizons/nsc_regulated.py','scripts/check_nsc_geometric_chain.py')},
        'declared_geometry':{'motif':'rho in [-R,R], w=1/sqrt(1+rho²), repeated periodically',
                             'R_is_input':True,'Omega_is_input':True,'independent_mass_inserted':False,
                             'seams':'w is continuous and periodic; derivative of metric profile changes at seams; no gravitational junction equations solved'},
        'zero_gap_exclusion':{'continuum_transfer_eigenvalues':'exp(plus_or_minus 2 asinh(R))',
                              'argument':'neither multiplier has modulus 1 for R>0; zero is outside every real Bloch fiber; compact phase and continuity give an open spectral gap',
                              'discrete_transfer':'product[(1-h*w/2)/(1+h*w/2)] in (0,1)',
                              'scope':'equal-scale periodic spatial family only'},
        'gap_families':families,
        'removed_angular_potential_control':{'kappa':0.,'zero_phase_gap':free_zero_gap,
                                             'interpretation':'formal control with w removed; kappa=0 is not a spinor-sphere eigenvalue'},
        'geometric_link':{'formula':'B[last_edge,first_node]=(1/h)+w_last/2; all other entries zero',
                          'radius':2.,'intervals':8,'entry':float(motif['B'][-1,0]),
                          'independently_fitted_link':False,'normalized_finite_chain_controls':controls},
        'nonclaims':{'stationary_geometry_or_R_Omega_derived':False,'periodic_spatial_array_is_a_Lorentzian_nested_cosmology':False,
                     'infinite_unequal_scale_tail_solved':False,'particle_species_or_cosmological_constant_predicted':False,
                     'physical_metric_stability_or_conservation_proved':False,'new_to_world_mechanism_established':False},
        'terminal':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__); group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args();result=calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()),result)
        print('Repeated-throat gap and geometric-link recursion reproduced; every field checked.')
    elif args.output:
        with args.output.open('x') as stream: json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write('\n')
        print(args.output)
    else: print(json.dumps(result,indent=2,allow_nan=False))


if __name__=='__main__':main()
