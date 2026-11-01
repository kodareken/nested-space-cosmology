"""Actual low-mode pump: complete signed Gaussian coincidence benchmark.

This reads existing spatial/onset records and caches. It adds no spectral
cutoff, geometry evolution or occupation fit, and does not compute kinetic
or renormalized stress.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_parent_circulation import Parent
from recursive_horizons.nsc_atlas_valve_neck import NeckInitialData
from recursive_horizons.nsc_atlas_valve_dirac_modes import mode_radial_profiles,initial_parent_mode_pump
from recursive_horizons.nsc_atlas_valve_source_gate import (
    geometric_spinor_candidate,signed_gaussian_coincidence,scalar_coherent_gaussian_control,
)

ROOT=Path(__file__).resolve().parents[1]
SOURCES=('src/recursive_horizons/nsc_atlas_valve_dirac_modes.py',
         'src/recursive_horizons/nsc_atlas_valve_source_gate.py','src/recursive_horizons/nsc_atlas_valve_neck.py',
         'src/recursive_horizons/nsc_atlas_parent_circulation.py','src/recursive_horizons/nsc_atlas_parent_scale_selector.py',
         'src/recursive_horizons/nsc_atlas_inherited_charge.py','scripts/check_nsc_atlas_valve_gaussian_source.py')
INPUTS=('results/summaries/nsc-atlas-valve-dirac-modes-v1.json',
        'results/summaries/nsc-atlas-valve-dirac-onset-v1.json')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calculate():
    start=time.perf_counter();before={p:sha(ROOT/p) for p in SOURCES}
    input_pins={p:sha(ROOT/p) for p in INPUTS}
    spatial,onset=[json.loads((ROOT/p).read_text()) for p in INPUTS]
    metric_record=spatial['defined_geometry']['raw_cache'];metric_path=ROOT/metric_record['path']
    if sha(metric_path)!=metric_record['sha256']:raise RuntimeError('defined initial geometry cache hash mismatch')
    with np.load(metric_path) as raw:
        coefficients=raw['coefficients'].copy();degree=int(raw['degree']);epsilon=float(raw['epsilon']);H=float(raw['H'])
    parent=Parent.retained()
    if parent.H!=H:raise RuntimeError('current parent normalization differs from the recorded mode geometry')
    initial=NeckInitialData(parent,degree,epsilon,coefficients,())
    candidate=geometric_spinor_candidate()
    tau1=np.array([[0,1],[1,0]],complex);tau2=np.array([[0,-1j],[1j,0]],complex)
    tau3=np.diag([1.,-1.]).astype(complex);I2=np.eye(2);I4=np.eye(4)
    gamma0=np.kron(tau3,I2);gamma5=np.kron(tau1,I2)
    rotations=[(I4-1j*a*gamma5)/math.sqrt(2) for a in (1,-1)]
    # Fixed diagnostic points, evaluating existing coefficients, not a new mode grid.
    u=np.linspace(-1.,1.,1001);Omega=initial.values(u)**2/H;channels=[]
    with threadpool_limits(limits=1):
        for item in spatial['levels'][-1]['lowest_by_sign']:
            path=ROOT/item['raw_mode_cache']['path']
            if sha(path)!=item['raw_mode_cache']['sha256']:raise RuntimeError('input spatial mode hash mismatch')
            with np.load(path) as raw:
                mode={k:raw[k].copy() for k in ('degree','k','sign','eigenvalue','coefficients')}
            lam=float(mode['eigenvalue']);mass=float(item['declared_mass']);E=math.hypot(mass,lam)
            if abs(mass-candidate['lowest_mass'])>2e-14:raise RuntimeError('candidate mass differs from recorded mode mass')
            pump=initial_parent_mode_pump(initial,mode,mass,(320,140))
            M=pump['total_pair_Hessian']
            if int(mode['sign'])>0 and abs(M-onset['initial_pump_quadrature'][-1]['total_pair_Hessian'])>2e-13:
                raise RuntimeError('actual positive pump does not match its input onset record')
            F,G=mode_radial_profiles(mode,u);rho=(F*F+G*G)/(4*math.pi*Omega**2)
            a=math.sqrt((E+mass)/(2*E));b=lam/math.sqrt(2*E*(E+mass));Q=np.array([[a,-b],[b,a]])
            C3=Q@(-M*tau2/6)@Q.T;C4=Q@(E*M*tau1/12)@Q.T
            S3=np.kron(C3@tau3,I2);S4=np.kron(C4@tau3,I2)
            channels.append(dict(sign=int(mode['sign']),lam=lam,mass=mass,energy=E,M=M,rho=rho,S3=S3,S4=S4,
                                 cache=item['raw_mode_cache'],
                                 matrix3_error=float(np.max(abs(S3+1j*M*gamma5/6))),
                                 matrix4_error=float(np.max(abs(S4+M*(lam*I4+mass*gamma0@gamma5)/12)))))
        external3=sum(c['rho'][:,None,None]*c['S3'] for c in channels)
        external4=sum(c['rho'][:,None,None]*c['S4'] for c in channels)
        completed3=sum(R@external3@R for R in rotations)
        completed4=sum(R@external4@R for R in rotations)
        partial4=[sum(R@(c['rho'][:,None,None]*c['S4'])@R for R in rotations) for c in channels]
        plus=next(c for c in channels if c['sign']>0);minus=next(c for c in channels if c['sign']<0)
        scalar=np.trace(external4,axis1=1,axis2=2)/4
        nonscalar=external4-scalar[:,None,None]*I4
        norm=lambda a:np.sqrt(np.sum(abs(a)**2,axis=(-2,-1)))
        controls=[]
        for name,C in (('tau1',tau1),('tau2',tau2),('tau3',tau3),('combination',.17*tau1-.29*tau2+.31*tau3)):
            control=signed_gaussian_coincidence(C)
            controls.append(dict(input=name,external_norm=control['external_norm'],complete_norm=control['complete_norm'],
                                 scalar_coefficient=control['scalar_coefficient']))
        # This control is a different local Gaussian state, not an inherited global state.
        coherent=scalar_coherent_gaussian_control()
        compact_coherent={key:coherent[key] for key in (
            'Clifford_dimension','spinor_dimension','chiral_commutator_norm','rank3_trace_coefficient',
            'literal_Wick_three_form','contact_lagrangian_over_kappa',
            'contact_coefficient_per_squared_scalar_coincidence','reference_is_vacuum',
            'global_parent_state_selected','UV_completion_selected','radial_NEC_computed','scope')}
        compact_coherent['covariance_eigenvalue_range']=[float(np.min(coherent['covariance_eigenvalues'])),
                                                        float(np.max(coherent['covariance_eigenvalues']))]
        compact_coherent['illustrative_local_control_c']=.1
        compact_coherent['illustrative_profile_density']=1.
        compact_coherent['not_used_to_fit_or_prepare_the_actual_parent']=True
        K=coherent['covariance'];grading=coherent['chirality'];xi=coherent['chiral_commutator']
        U=(np.eye(len(K))+1j*grading)/math.sqrt(2)
        final=U@K@U.conj().T
        inheritance=dict(formula='Xi_out=orientation*U Xi_in Udagger when Gamma_out U=orientation U Gamma_in.',
            unitary_is_an_algebraic_control_not_the_actual_parent_propagator=True,
            unitary_error=float(np.max(abs(U@U.conj().T-np.eye(len(K))))),
            initial_Xi_norm=float(np.linalg.norm(xi)),
            orientation_controls=[dict(orientation=e,
                Xi_transport_error=float(np.max(abs(final@(e*grading)-(e*grading)@final-e*U@xi@U.conj().T))))
                for e in (1,-1)],
            state_selected_by_operator_inheritance=False,
            scope='Finite graded-unitary matrix identity; neither an actual horizon channel nor initial-state selection.')
    import sympy as sp
    sn,sH,sc,sm,sms,sdot=sp.symbols('n H c m m_S Sdot')
    density=sm*sn-sc*sn**2;pressure=-sc*sn**2;ndot=-3*sH*sn
    ward=sp.simplify(sp.diff(density,sn)*ndot+sp.diff(density,sm)*sms*sdot+3*sH*(density+pressure)-sms*sn*sdot)
    direct_missing=sp.simplify(sp.diff(-sc*sn**2,sn)*ndot+3*sH*(-sc*sn**2+sc*sn**2))
    ward_control=dict(full_residual=str(ward),direct_contact_only_missing_term=str(direct_missing),
        assumptions='Homogeneous algebra sanity: n=N/a^3,rho=m(S)n-c n^2,p=-c n^2,c constant.',
        actual_parent_equation_of_state_derived=False)
    if ward!=0:raise RuntimeError('independent homogeneous Ward control failed')
    Vgamma=2*math.pi**6.5/math.gamma(6.5)*candidate['internal_radius_squared']**6
    factor=math.exp(6*candidate['s_star'])/Vgamma
    norms=dict(spatial_tau3_max_norm=float(np.max(norm(external3))),
               spatial_tau4_nonscalar_max_norm=float(np.max(norm(nonscalar))),
               spatial_tau4_scalar_equator=float(scalar[len(u)//2].real),
               spatial_tau4_scalar_sampled_max_abs=float(np.max(abs(scalar))),
               internal_complete_single_spatial_sign_tau4_max_norm=float(max(np.max(norm(v)) for v in partial4)),
               all_signs_tau3_max_norm=float(np.max(norm(completed3))),
               all_signs_tau4_max_norm=float(np.max(norm(completed4))),
               full256_tau3_max_Frobenius=float(8*factor*np.max(norm(completed3))),
               full256_tau4_max_Frobenius=float(8*factor*np.max(norm(completed4))))
    after={p:sha(ROOT/p) for p in SOURCES}
    if before!=after or input_pins!={p:sha(ROOT/p) for p in INPUTS}:
        raise RuntimeError('source or preserved input record changed during calculation')
    return dict(schema='nsc-atlas-valve-gaussian-source-v1',
                provenance=dict(source_sha256_before=before,source_sha256_after=after,source_unchanged=True,
                                input_record_sha256=input_pins,defined_metric_cache=metric_record,
                                git_parent_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()),
                conventions=dict(C4='Here C_4 is a four-component coincident matrix, not a connected four-point cumulant.',
                    covariance='C>=<psi psi^dagger>; DeltaS4=DeltaC> gamma0. Fixed initial canonical low-mode vacuum P+=(I+tau3)/2 tensorIspin.',
                    Clifford='gamma0=tau3 tensorI2,gamma5=tau1 tensorI2.',
                    mass_rotation='Q_lambda=[[A,-B],[B,A]],A=sqrt((E+m)/(2E)),B=lambda/sqrt(2E(E+m)).',
                    internal_rotation='R_a=(I-i*a*gamma5)/sqrt2,psi_a=R_a chi_a,barpsi_a=barchi_a R_a.',
                    complete_signs='Both spatial Dirac signs and both internal Killing signs; all64 internal labels per sign and both angular spin labels.',
                    internal_coincidence='Each complete internal projector isI64/Vgamma; DeltaS16=e6s/Vgamma sum_a R_a DeltaS4 R_a tensorI64.',
                    normal_order='Complete relative vacuum Wick ordering for the declared finite Gaussian benchmark, not an absolute renormalization prescription.'),
                actual_parent_channels=[dict(sign=c['sign'],spatial_eigenvalue=c['lam'],mass=c['mass'],energy=c['energy'],
                    pair_Hessian=c['M'],matrix_tau3_identity_error=c['matrix3_error'],matrix_tau4_identity_error=c['matrix4_error'],
                    input_mode_cache=c['cache']) for c in channels],
                matrix_coefficients=dict(energy_basis='DeltaC_E=-(M/6)tau2*tau³+(E*M/12)tau1*tau4+...',
                    per_spatial_sign='DeltaS_sigma=-i*rho_sigma*M_sigma*gamma5*tau³/6-rho_sigma*M_sigma*(lambda_sigma I+m gamma0 gamma5)*tau4/12+...',
                    explanation='The individual tau4 covariance is nonscalar; cancellation is checked after actual sign completion.'),
                independent_partner_comparison=dict(sampled_max_density_difference=float(np.max(abs(plus['rho']-minus['rho']))),
                    eigenvalue_sum=plus['lam']+minus['lam'],pair_Hessian_sum=plus['M']+minus['M']),
                raw_matrix_norms=norms,exact_generality_controls=controls,
                conditional_identity=dict(hypotheses=[
                    'Identical positive-mass preparation of the complete internal signs and common symmetry-preserving relative subtraction.',
                    'Spatial partner covariance C_-=tau3 C_+ tau3, unpolarized spin, equal exact partner densities and tracelessDeltaC.',
                    'The symmetrically completed vacuum is chiral invariant; the geometry alone has not selected this complete initial state.'],
                    spatial_sum='rho{DeltaC,tau3}=2rho*c3 I4',internal_sum='R+ I R+ +R- I R-=0',
                    result='Exact coincident covariance/contact cancellation for this balanced same-mode Gaussian benchmark.',
                    caveat='The independent L24 modes have finite residuals, so their small raw cancellation error is retained, not promoted to a physical source.'),
                bare_contact=dict(functional='sum_ABC[tr(GammaABC DeltaS16)tr(Gamma^ABC DeltaS16)-tr(GammaABC DeltaS16 Gamma^ABC DeltaS16)]',
                    physical_U1_anomalous_field_pairing_zero=True,
                    tau6_bare_contact_energy_EL_coefficient_over_g=0.,
                    balanced_closed_sector_bare_contact_energy_EL_over_g=0.,
                    energy_units='Dimensionless contact energyE*L scales asg=kappa_nat/L².',
                    Einstein_source_units='Einstein-normalized contact and induced kinetic response scale together asg²; free fermion stress starts asg.',
                    total_kinetic_stress_or_renormalized_source_computed=False),
                local_chiral_coherent_Gaussian_counterexample=compact_coherent,
                chiral_inheritance_control=inheritance,
                homogeneous_Ward_control=ward_control,
                bare_action_Ward=dict(action='L_D=K-sum_a m_a(S) B_a+C4, with the full matched quartic vertex.',
                    off_shell_stress='T_D,mu nu=Bkin_mu nu+g_mu nu L_D; vary before using field equations.',
                    direct_contact='T_contact,direct_mu nu=g_mu nu C4, hence directT_kk=0.',
                    on_shell='K=sum mB-2C4,L_D,on=-C4,T_D,kk=Bkin_kk evaluated with the interacting equation.',
                    radion_exchange='nabla^mu T_D,mu nu=(-sum m_a,S B_a+partial_S C4)partial_nu S; barepartial_S C4=0.',
                    normalization='m_S=-7m/sqrt84 in the retainedkappa=1 normalization; scalar stress has opposite exchange.',
                    counterterm_scope='The displayed identities are bare. A renormalized source must also vary its counterterms: J_S=sum m_S<B>ren-(delta Sct/delta S)/sqrt(-g), with div<T_D+T_ct>=-J_S*dS and the opposite radion source. No counterterm value is selected here.',
                    perturbative_order='Contact Wick term and induced kinetic/radion response must be retained at the same order.',
                    physical_stress_result_claimed=False),
                verdict=dict(balanced_Gaussian_contact_cancelled=True,all_Gaussian_states_contact_zero=False,
                    nonGaussianity_required_for_every_possible_source=False,
                    chiral_coherent_Gaussian_preparation_selected=False,connected_four_point_computed=False,
                    intermode_and_UV_Hadamard_completion_selected=False,g_selected=False,
                    radial_NEC_computed=False,turnaround_decided=False,baryon_photon_eta=None,
                    scope='Finite balanced same-mode Gaussian benchmark, not the full kinetic stress, whole quantum state or source no-go.'),
                numerical=dict(existing_spatial_degree=24,diagnostic_profile_samples=len(u),
                    new_spectral_or_geometry_evolution_grids=0,elapsed_seconds=time.perf_counter()-start,
                    raw_arrays_committed=False),
                reproduction='python scripts/run_python.py scripts/check_nsc_atlas_valve_gaussian_source.py --output results/summaries/nsc-atlas-valve-gaussian-source-v1.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/summaries/nsc-atlas-valve-gaussian-source-v1.json')
    args=parser.parse_args();result=calculate();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(args.output),norms=result['raw_matrix_norms'],controls=result['exact_generality_controls'],
                          alternative=result['local_chiral_coherent_Gaussian_counterexample']),indent=2))


if __name__=='__main__':main()
