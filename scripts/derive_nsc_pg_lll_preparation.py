#!/usr/bin/env python3
"""Bind the existing transparent horizon state to an actual PG covariance."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import expit

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_lll_preparation import PGLLLPreparation, require_full_retained_preparation
from recursive_horizons.nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain, MODE_TO_CURRENT
from recursive_horizons.nsc_transmitting_ctp_resolvent import CTPRetardedMemory
from derive_nsc_transmitting_boundary_binding import compare

OUTPUT='results/development/nsc-pg-lll-preparation.json'
SOURCES=('scripts/derive_nsc_pg_lll_preparation.py','src/recursive_horizons/nsc_pg_lll_preparation.py',
         'tests/test_nsc_pg_lll_preparation.py','docs/nsc-pg-lll-preparation.md')
INPUTS=('results/development/charged-ctp-neck-source.json','results/development/nsc-mode-resolved-cauchy-state.json',
        'results/development/nsc-transmitting-cross-resolvent.json','results/development/nsc-transmitting-ctp-resolvent.json',
        'results/development/nsc-transmitting-boundary-remainder.json',
        'docs/nsc-unruh-state.md','docs/nsc-transmitting-dirac-domain.md',
        'src/recursive_horizons/nsc_lorentzian.py','src/recursive_horizons/nsc_unruh_state.py',
        'src/recursive_horizons/nsc_transmitting_dirac_domain.py','src/recursive_horizons/nsc_transmitting_ctp_resolvent.py',
        'scripts/derive_nsc_transmitting_boundary_binding.py')


def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def pairs(value):
    a=np.asarray(value);return np.stack((a.real,a.imag),axis=-1).tolist()


def calculate():
    charged,state,response,ctp,remainder=[json.loads((ROOT/p).read_text()) for p in INPUTS[:5]]
    config=charged['runs']['base']['config'];locked=response['locked_inputs']
    if config['omega']!=locked['Omega'] or config['magnetic_flux']!=locked['magnetic_flux']:
        raise ValueError('locked state parameters differ')
    payload_paths=(state['payload']['path'],response['payload']['path'])
    for rec in (state,response):
        if sha(rec['payload']['path'])!=rec['payload']['sha256']:raise ValueError('input artifact changed')
    modes=np.load(ROOT/payload_paths[0],allow_pickle=False);old=np.load(ROOT/payload_paths[1],allow_pickle=False)
    p=PGLLLPreparation(config['horizon_rho'],config['surface_gravity'],config['omega'])
    require_full_retained_preparation(0)
    lll=state['channels'][0];start=lll['sample_offset'];stop=start+lll['sample_count']
    frequency=modes['frequency'][start:stop];original=modes['covariance_seed'][start:stop]
    reconstructed=np.array([p.seed_restriction(e) for e in frequency])
    signed=np.r_[-frequency[::-1],frequency]
    eigenvalues=np.array([np.linalg.eigvalsh(p.spectral_covariance(e)) for e in signed])
    # One spatial quadrature refinement of the new preparation only.
    base=p.project(128);fine=p.project(256);c=fine['covariance'];gram=fine['CAR_gram']
    z=complex(*response['parameters']['z'])
    memory=CTPRetardedMemory(old['upper_response'][0],old['lower_response'][0])
    characteristic=p.packet_retarded_resolvent(z,48)
    atlas=fine['atlas'];locations=np.array([-.9,-.2,.3,.8,3.2,3.8]);h=2e-4
    derivative_error=0.
    for sign in (-1,1):
        derivative=(-atlas.coordinate(sign,locations+2*h)+8*atlas.coordinate(sign,locations+h)
                    -8*atlas.coordinate(sign,locations-h)+atlas.coordinate(sign,locations-2*h))/(12*h)
        derivative_error=max(derivative_error,float(np.max(abs(derivative-1/atlas.velocity(sign,locations)))))
    # The seed trace is recovered with the same T half-density/current map.
    beta0=float(np.sqrt(3*np.pi/2));seam=TransmittingDiracSeamDomain(1.,1.,beta0,1.)
    trace,inverse,_=seam.trace_map(np.ones(1))
    amplitude=MODE_TO_CURRENT@np.diag([1/np.sqrt(beta0-1),1/np.sqrt(beta0+1)])
    trace_error=float(np.linalg.norm(amplitude-trace[0]))
    pullback=np.array([inverse[0]@amplitude@cc@amplitude.conj().T@inverse[0].conj().T for cc in reconstructed])
    point_pairs=((.4,-.6),(3.5,-.6),(-2.,5.))
    hermitian_kernel=max(float(np.linalg.norm(p.off_diagonal_kernel(a,b)-p.off_diagonal_kernel(b,a).conj().T)) for a,b in point_pairs)
    # Resolve the Fourier sign/normalization at actual flow separations. The
    # filled-sea distribution is kept analytically, not cut at old omega_max.
    fourier_error=0.
    for beta in (p.beta_out,p.beta_in):
        distance=float(atlas.coordinate(-1,.4)-atlas.coordinate(-1,-.6))
        correction=quad(lambda e:float(expit(-beta*e))*np.sin(e*distance),0.,np.inf,epsabs=2e-13,epsrel=2e-13)[0]
        direct=-1j/(2*np.pi*distance)+1j/np.pi*correction
        analytic=-1j/(2*beta*np.sinh(np.pi*distance/beta))
        fourier_error=max(fourier_error,float(abs(direct-analytic)))
    distance=float(atlas.coordinate(1,3.5)-atlas.coordinate(1,-.6))
    integral=quad(lambda e:np.exp(-p.beta_out*e/2)/(1+np.exp(-p.beta_out*e))*np.cos(e*distance),0.,np.inf,epsabs=2e-13,epsrel=2e-13)[0]
    fourier_error=max(fourier_error,float(abs(-1j/np.pi*integral+1j/(2*p.beta_out*np.cosh(np.pi*distance/p.beta_out)))))
    eig=np.linalg.eigvalsh(c)
    residuals={
        'seed_occupation_reconstruction':float(np.max(abs(reconstructed-original))),
        'T_trace_amplitude':trace_error,
        'T_seed_covariance_pullback':float(np.max(abs(pullback-original))),
        'spectral_CAR_lower_violation':max(0.,-float(eigenvalues.min())),
        'spectral_CAR_upper_violation':max(0.,float(eigenvalues.max())-1.),
        'spatial_probe_CAR_gram':float(np.linalg.norm(gram-np.eye(7))),
        'projected_CAR_lower_violation':max(0.,-float(eig.min())),
        'projected_CAR_upper_violation':max(0.,float(eig.max())-1.),
        'covariance_hermiticity':float(np.linalg.norm(c-c.conj().T)),
        'kernel_hermiticity':hermitian_kernel,
        'spatial_quadrature_refinement':float(np.linalg.norm(c-base['covariance'])),
        'characteristic_response_vs_locked_R':float(np.linalg.norm(-characteristic-memory.G_retarded)),
        'flow_coordinate_derivative':derivative_error,
        'thermal_Fourier_normalization':fourier_error,
        'initial_CTP_CAR_identity':float(np.linalg.norm(fine['greater_equal_time']-fine['lesser_equal_time']+1j*gram)),
    }
    tolerance=3e-11
    if max(residuals.values())>tolerance:raise ArithmeticError(json.dumps(residuals,indent=2))
    if np.linalg.norm(c[:4,4:])<1e-6:raise ArithmeticError('packet/bulk preparation correlations missing')
    return {
        'schema':'NSC-PG-LLL-PREPARATION-v1',
        'status':'LLL physical PG Cauchy preparation PASS; full retained C1b and C2-C4 OPEN',
        'source_hashes':{f:sha(f) for f in SOURCES},
        'input_hashes':{f:sha(f) for f in (*INPUTS,*payload_paths)},
        'locked_inputs':locked,
        'domain':{
            'operator':'same whole-line Z1 PG Dirac operator restricted to m=lambda=0',
            'surface':'one common tau=constant PG Cauchy surface; fixed metric and rho=0 cut',
            'flow_sectors':['sigma2+ exterior rho>rho_h','sigma2+ interior rho<rho_h','sigma2- whole line'],
            'new_particle_sectors':False,
            'flow_map':'dx_s/drho=1/a_s; chi_s(rho)=tilde_chi_s(x_s)/sqrt(abs(a_s))',
            'relative_affine_phase':'same logarithmic finite part across horizon; common translations cancel',
            'full_global_LLL_state_defined':True,
            'massive_channels_bound':False,
        },
        'state_definition':{
            'input':'owned affine-horizon C_H plus inherited incoming Fermi occupation',
            'frequency_order':['outgoing exterior','outgoing interior','incoming'],
            'C_E':'[[f,-i sqrt(f(1-f)),0],[i sqrt(f(1-f)),1-f,0],[0,0,n]]',
            'f':'expit(-2*pi*E/kappa_h)','n':'expit(-2*pi*E/(Omega*kappa_h))',
            'signed_frequencies':'full real E; original thermal law, no finite-frequency CAR truncation',
            'spectral_eigenvalues':'0,1,n(E)',
            'Fourier_sign':'exp(+i E (x-y)); H=-i partial_x, U(t)=exp(-i H t)',
            'contact_term':'I delta(rho-rho_prime)/2 retained in smeared covariance',
            'seed_covariance_read_only_as_verification':True,
            'seed_covariance_copied_to_spatial_basis':False,
            'horizon_rho':p.horizon_rho,'surface_gravity':p.surface_gravity,
            'beta_out':p.beta_out,'beta_in':p.beta_in,
        },
        'projection':{
            'order':['parent spin0','parent spin1','child spin0','child spin1','child_bulk spin0','child_bulk spin1','exterior_bulk sigma2+'],
            'original_J_unchanged':True,
            'additional_probes':'sqrt210*t(1-t)(2t-1) with t=-rho on(-1,0), two spins; sqrt30*t(1-t) with t=rho-3 on(3,4), sigma2+',
            'bulk_probes_are_not_new_source_channels':True,
            'finite_extra_probes_exhaust_complement':False,
            'covariance_real_imag':pairs(c),
            'CAR_gram_real_imag':pairs(gram),
            'eigenvalues':eig.tolist(),
            'parent_child_correlation_norm':float(np.linalg.norm(c[:2,2:4])),
            'J_to_sampled_bulk_correlation_norm':float(np.linalg.norm(c[:4,4:])),
            'J_to_outgoing_exterior_correlation_norm':float(np.linalg.norm(c[:4,6:])),
            'lesser_equal_time_real_imag':pairs(fine['lesser_equal_time']),
            'greater_equal_time_real_imag':pairs(fine['greater_equal_time']),
            'Keldysh_equal_time_real_imag':pairs(fine['Keldysh_equal_time']),
            'spatial_quadrature_points_per_interval':[128,256],
        },
        'verification':{'residuals':residuals,'tolerance':tolerance,'locked_LLL_seed_nodes_compared':len(frequency)},
        'gate':{
            'C1b_LLL':'PASS','C1b_full_retained':'OPEN','prepared_channel_indices':[0],
            'missing_channel_indices':list(range(1,33)),
            'massive_next_owner':'global signed-frequency spin/angular mode reconstruction under the existing horizon/parent preparation',
            'C2':'OPEN: full transmitting history derivatives and KS endpoint pullback not supplied',
            'transmitting_EndpointBranchJets':None,
            'C3':'OPEN','Gamma_rest_boundary_derivative':None,'physical_Weyl_mismatch':None,
            'locked_Weyl_coefficients':remainder['locked_local_components'],
            'Weyl_extraction_residual':remainder['locked_extraction_residual'],
            'preserved_Weyl_value':remainder['preserved_Weyl_value'],
            'C4':'OPEN','physical_Vc_or_stationary_history':None,
        },
        'scope':{'parameters_refitted':False,'new_physical_terms':False,'prior_scientific_files_modified':False,
                 'old_generators_rerun':False,'finite_stress_or_nulls_assigned':False,
                 'metric_timestep_started':False,'physical_history_selected':False,
                 'full_C1_substrate_claimed':False,'Z3':'OUT OF SCOPE'},
        'comparison':response['comparison'],
    }


def main():
    p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group()
    g.add_argument('--check',action='store_true');g.add_argument('--write',action='store_true');a=p.parse_args()
    r=calculate();path=ROOT/OUTPUT
    if a.check:
        compare(json.loads(path.read_text()),r);print('LLL PG preparation and new CTP data verified; full retained C1b OPEN')
    elif a.write:
        if path.exists():raise FileExistsError('immutable result already exists')
        path.write_text(json.dumps(r,indent=2,sort_keys=True)+'\n');print(r['status'],r['verification'])
    else:print(json.dumps(r,indent=2,sort_keys=True))
    return 0


if __name__=='__main__':raise SystemExit(main())
