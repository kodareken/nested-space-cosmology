#!/usr/bin/env python3
"""Reproduce causal fixed-background PG Dirac transport and its flux ledger."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_lorentzian import evolve,transport_control,geometry,horizon
from check_nsc_scale_closure import compare
OUTPUT=ROOT/'results/nsc-4-lorentzian-transport.json'


def compact(row):
    sample=np.linspace(0,row['points']-1,65,dtype=int)
    result={k:v for k,v in row.items() if k not in ('field','x','initial')}
    for key in ('incoming_channels','outgoing_channels'):
        result[key]=[list(item) for item in result[key]]
    result['samples']={'rho':row['x'][sample].tolist(),
                       'field_real':row['field'][:,sample].real.tolist(),
                       'field_imag':row['field'][:,sample].imag.tolist()}
    return result


def calculate():
    runs=[evolve(n,final_time=1.2) for n in (801,1601,3201)]
    errors=[float(np.sqrt(a['h']*np.sum(abs(a['field']-b['field'][:,::2])**2)))
            for a,b in zip(runs,runs[1:])]
    assert errors[1]<errors[0]
    for row in runs:
        assert abs(row['probability_balance_residual'])<1e-8
        assert abs(row['child_balance_residual'])<1e-8
        assert row['incoming_channels']==[(1,-1)]
        assert row['boundary_form_residual']<1e-12
    assert runs[-1]['child_probability']+runs[-1]['left_child_outflow']>.9999
    assert runs[-1]['outside_horizon_probability']<1e-10
    tails=[row['outside_characteristic_cone_probability'] for row in runs]
    assert tails[2]<tails[1]<tails[0]
    time_runs=[runs[0],evolve(801,final_time=1.2,cfl=.225),evolve(801,final_time=1.2,cfl=.1125)]
    time_errors=[float(np.sqrt(a['h']*np.sum(abs(a['field']-b['field'])**2)))
                 for a,b in zip(time_runs,time_runs[1:])]
    time_order=float(np.log2(time_errors[0]/time_errors[1]))
    assert 3.95<time_order<4.05
    transport=[]
    for n in (801,1601,3201):
        row=evolve(n,final_time=.4,kappa=0.)
        exact=transport_control(row['x'],.4)/row['initial_normalization']
        error=float(np.sqrt(row['h']*np.sum(abs(row['field'][0]-exact)**2)))
        transport.append({'points':n,'L2_error':error,'other_channel_norm':float(row['h']*np.sum(abs(row['field'][1])**2))})
    assert transport[2]['L2_error']<transport[1]['L2_error']<transport[0]['L2_error']
    outflow=evolve(801,final_time=.35,center=-5,half_width=.5)
    assert outflow['physical_outflow']>.999 and abs(outflow['probability_balance_residual'])<1e-6
    probe=[-8.,0.,horizon(),8.]
    speeds=[]
    for x in probe:
        beta,db,A=geometry(x)
        speeds.append({'rho':x,'A':float(A),'beta':float(beta),'beta_prime':float(db),
                       'speeds':[float(1-beta),float(-1-beta)]})
    return {'schema':'nsc-lorentzian-transport-v1','artifact_id':'NSC-4-LORENTZIAN-TRANSPORT',
            'classification':'causal_massless_Dirac_probe_on_fixed_black_universe_with_correct_characteristic_outflow_and_conserved_Dirac_norm_transfer',
            'source_hashes':{p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in (
                'src/recursive_horizons/nsc_lorentzian.py','scripts/check_nsc_lorentzian.py')},
            'independent_tetrad_record_sha256':hashlib.sha256((ROOT/'results/nsc-4-dirac-tetrad.json').read_bytes()).hexdigest(),
            'model':{'geometry':'r=sqrt(1+rho²), A=1+3rho+3(1+rho²)(atan(rho)-pi/2), beta=sqrt(1-A)',
                     'metric':'ds²=d_tau²-(d_rho+beta*d_tau)²-r²dOmega²',
                     'Hamiltonian':'-i(sigma2-beta*I)d_rho+i*beta_prime/2*I+(kappa/r)*sigma1',
                     'current':'u_dagger*(sigma2-beta*I)*u',
                     'basis':'sigma2 characteristic eigenbasis; normal-slice radial L² measure',
                     'units':'c=hbar=throat_radius=1',
                     'field_content':'one fixed kappa=1 angular sector of a massless Dirac probe',
                     'state':'compact Cauchy packet, not a vacuum or stationary quantum state',
                     'box':[-8.,8.],'initial_packet':{'center':.4,'half_width':.25,'momentum':4.,'characteristic_channel':0},
                     'boundaries':'left trapped cut: no incoming data; right parent cut: homogeneous incoming minus field',
                     'method':'SBP split-form central space, dissipative inflow SAT, RK4 with matching-stage flux quadrature',
                     'cfl':.45},
            'horizon':horizon(),'characteristic_probes':speeds,
            'coupled_runs':[compact(row) for row in runs],
            'spatial_convergence':{'successive_grid_L2_differences':errors,'observed_order':float(np.log2(errors[0]/errors[1])),
                                   'continuum_error_bound_proved':False,
                                   'warning':'narrow-packet profile convergence is pre-asymptotic; integrated charge is better resolved than pointwise phase'},
            'temporal_convergence':{'points':801,'cfl_values':[.45,.225,.1125],
                                    'successive_grid_L2_differences':time_errors,
                                    'observed_order_rounded_2dp':round(time_order,2),
                                    'accepted_order_interval':[3.95,4.05]},
            'independent_half_density_transport':{'kappa':0,'scope':'formal numerical control, not an allowed spinor-sphere eigenvalue',
                                                  'formula':'u(t,x)=sqrt(dY/dx)*u0(Y), Y=backward_characteristic_flow',
                                                  'cases':transport},
            'inner_outflow_control':compact(outflow),
            'nonclaims':{'background_sourced_by_same_fermions':False,'quantum_vacuum_selected':False,
                         'physical_metric_Hessian_verified':False,'positive_Killing_energy_in_trapped_region':False,
                         'total_gravitating_energy_balance_derived':False,'stationary_scale_or_Phi_selected':False,
                         'full_black_to_child_formation_solved':False,'new_to_world_physics_discovery':False},
            'terminal':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument('--check',action='store_true');group.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError('refusing to overwrite result')
    record=calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()),record)
        print('Lorentzian Dirac propagation, flux ledger and characteristic controls reproduced; every field checked.')
    elif args.output:
        with args.output.open('x') as stream:json.dump(record,stream,sort_keys=True,indent=2,allow_nan=False);stream.write('\n')
        print(args.output)
    else:print(json.dumps(record,indent=2,allow_nan=False))


if __name__=='__main__':main()
