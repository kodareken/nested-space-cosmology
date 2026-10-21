"""Constrained unaveraged parent growth, compact packet and scale dependence."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import sympy as sp
from recursive_horizons.nsc_atlas_parent_perturbations import ParentModes,ParentPacket,fixed_metric_control,jeans_comparison

ROOT=Path(__file__).resolve().parents[1]


def identities():
    a,H,p,x,w,kap,K,l,X,Y,F=sp.symbols('a H p x w kap K l X Y F',nonzero=True)
    hd=K/a**2-kap*p*p/2;pd=-3*H*p-w*w*x
    fd=-H*F+kap*p*X/2
    yd=-3*H*Y-(w*w+l/a**2)*X+4*p*fd-2*w*w*x*F
    vars=(a,H,p,x,X,Y,F);rates=(a*H,hd,pd,p,Y,yd,fd)
    dt=lambda z:sum(sp.diff(z,v)*r for v,r in zip(vars,rates))
    rho=p*Y-p*p*F+w*w*x*X;pressure=rho-2*w*w*x*X
    C=(l-3*K)*F/a**2+kap*(rho+3*H*p*X)/2
    T=sp.symbols('T')
    raw={
        'Hamilton_constraint_propagation':dt(C)+3*H*C,
        'independent_spatial_Einstein_equation':dt(fd)+4*H*fd+(2*hd+3*H*H-K/a**2)*F-kap*pressure/2,
        'perturbed_scalar_energy_conservation':dt(rho)+3*H*(rho+pressure)-3*p*p*fd+l*p*X/a**2,
        'comoving_density_gauge_invariance':(rho+3*H*p*p*T)+3*H*p*(X-p*T)-(rho+3*H*p*X),
        'initial_metric_and_scalar_constraint':(4*(3*K-l)*F/a**2-2*kap*w*w*x*X).subs(X,2*(3*K-l)*F/(kap*a**2*w*w*x)),
    }
    eps,e,q=sp.symbols('eps e q',positive=True)
    raw['two_step_inherited_excitation']=w*w*((sp.sqrt(1-eps)/w)*(2*e))**2/2-2*(1-eps)*e*e
    z=sp.symbols('z');phi=sp.Function('phi')(z)
    # Laplace-Beltrami operator in x=cos(theta) on radial S3 functions.
    theta=sp.symbols('theta',real=True)
    for n in (2,3):
        f=sp.chebyshevu(n,z)/(n+1)
        raw[f'S3_radial_eigenvalue_{n}']=(1-z*z)*sp.diff(f,z,2)-3*z*sp.diff(f,z)+n*(n+2)*f
    results={k:str(sp.simplify(v)) for k,v in raw.items()}
    if any(v!='0' for v in results.values()):raise RuntimeError(results)
    return results


def calculation():
    start=time.process_time();m=ParentModes();p=ParentPacket(m)
    coarse=m.solve();fine=m.solve(clock='efolds',rtol=2e-11,atol=2e-13)
    summaries=[x.summary() for x in (coarse,fine)]
    if (max(x['max_normalized_Hamilton'] for x in summaries)>3e-9
            or max(x['max_density_constraint_gap'] for x in summaries)>2e-8
            or max(x['max_Friedmann'] for x in summaries)>2e-9):
        raise RuntimeError('linear or background gravitational constraints exceed tolerance')
    gap=max(abs(a[k]-b[k]) for a,b in zip(summaries[0]['rows'],summaries[1]['rows'])
        for k in ('density_seed_final_Delta','phase_seed_final_Delta','density_seed_potential_ratio'))
    theta=np.linspace(0,2.85,801);packet_rows=[];last=None
    for cutoff in (64,128,256):
        sol,coeff,initial=p.evolve(cutoff)
        end=p.reconstruct(sol,initial,sol.solution.t[-1],theta)
        first=p.reconstruct(sol,initial,0,theta)
        eta=float(sol.sample([sol.solution.t[-1]])['conformal'][0])
        row=dict(cutoff=cutoff,initial_center=float(first['Delta'][0]),final_center=float(end['Delta'][0]),
            final_center_metric_potential=float(end['Phi'][0]),
            maximum_initial_density_reconstruction_error=float(np.max(abs(first['Delta']-p.profiles(theta)['Delta']))),
            final_maximum_tail_beyond_causal_support=float(np.max(abs(end['Delta'][theta>p.support_angle+eta]))),
            density_difference_from_previous_cutoff=None if last is None else float(np.max(abs(last['Delta']-end['Delta']))),
            initial_lowest_mode_coefficients=coeff[:2].tolist(),
            coefficient_sha256=hashlib.sha256(np.asarray(coeff,dtype='<f8').tobytes()).hexdigest())
        packet_rows.append(row);last=end
    # Independent time formulation for the retained packet; no raw histories
    # or coefficient arrays are stored in the scientific summary.
    alternative,_,initial2=p.evolve(256,clock='efolds',rtol=2e-11,atol=2e-13)
    end2=p.reconstruct(alternative,initial2,alternative.solution.t[-1],theta)
    packet_clock_gap=float(np.max(abs(end2['Delta']-last['Delta'])))
    control=fixed_metric_control(sol)
    control_center=float(initial@control['Delta'])
    sample=sol.sample(np.linspace(0,sol.solution.t[-1],257))
    center_history=initial@sample['Delta'][:,0,:]
    domain=m.atlas.parent_causal_domain()
    paths=('src/recursive_horizons/nsc_atlas_parent_perturbations.py',
        'src/recursive_horizons/nsc_atlas_shared_scalar_slice.py',
        'src/recursive_horizons/nsc_atlas_flux_conversion.py',
        'scripts/check_nsc_atlas_parent_perturbations.py',
        'scripts/plot_nsc_atlas_parent_perturbations.py',
        'tests/test_nsc_atlas_parent_perturbations.py')
    record=dict(schema='nsc-atlas-parent-perturbations-v1',
        result='A regular, constraint-propagating linear Einstein-scalar system and a compact parent-supported perturbation with gravitational density growth and scale-dependent transfer.',
        inherited_background=m.initial,
        equations=dict(metric='ds2=-(1+2Phi)dt2+a2(1-2Phi)gamma_ij dxidxj; curvature K=1, no linear anisotropic stress',
            harmonics='lambda=n(n+2), n>=2; computational S3 extension, not wall spectrum',
            constraint='C=(lambda-3)*Phi/a2+kappa*(Pi*Y-Pi2*Phi+omega2*chi*X+3H*Pi*X)/2; Cdot=-3H*C',
            observable='Delta=(delta_rho+3H*Pi*X)/E=-2*(lambda-3)*Phi/(kappa*a2*E); lensing potential Phi+Psi=2Phi',
            regularity='No division by Pi or by a time-dependent constraint coefficient; all scalar turning points retained.'),
        mode_transfer=summaries,
        transfer_basis=dict(density='Delta0=1, X0=chi0/2, Y0=0; Phi0 fixed by Hamilton constraint',
            phase='Delta0=Phi0=X0=0, Y0=omega*chi0; free scalar velocity remains physical initial data',
            amplitude='Unit linear response columns; finite benchmark multiplies by epsilon=1e-5, without a certified nonlinear remainder.'),
        packet=dict(support_angle=p.support_angle,potential_regular_class='C7; scalar displacement C5',
            prescription='(cos(theta)-cos(theta_b))^8 times a degree-two polynomial orthogonal to 1 and cos(theta) in the weighted cap measure; zero outside; central initial Delta=1, Y=0.',
            polynomial_coefficients_ascending=p.polynomial.coef.tolist(),normalization=p.normalization,
            preserves_wall_initial_data=True,
            whole_future_causal_margin=domain['initial_parent_angle']-p.support_angle-2*domain['future_conformal_distance_upper_bound'],
            finite_endpoint_proper_time=float(sample['proper_time'][-1]),
            truncation_checks=packet_rows,independent_clock_max_density_gap=packet_clock_gap,
            central_growth_factor=packet_rows[-1]['final_center'],
            central_growth_percent=100*(packet_rows[-1]['final_center']-1),
            central_density_history_range=[float(np.min(center_history)),float(np.max(center_history))],
            fixed_metric_control_center=control_center,
            fixed_metric_control_scope=control['scope'],
            central_growth_due_to_metric_response_comparison=packet_rows[-1]['final_center']-control_center,
            scope='A constrained linear perturbation of the actual compact initial atlas, extended by zero near both walls. Its smooth-parent causal development is calculated; no nonlinear finite-amplitude completion, formation rate, observed power spectrum or full perturbed wall future is supplied.'),
        scale_comparison=jeans_comparison(m),
        checks=dict(symbolic_zero_identities=identities(),independent_mode_clock_gap=gap,
            retained_numerical_comparison_tolerance=1e-7,
            meaning='Finite numerical consistency tolerance; spectral differences and sampled residuals are recorded separately, not a rigorous uniform enclosure.'),
        limits='Initial perturbation shape and scalar phase are specified, not derived from inheritance. The sampled longer modes grow weakly; shorter modes oscillate and are preparation dependent. The Jeans comparison uses averaged short-curvature-scale mathematics only; the evolution itself is unaveraged and includes positive spatial curvature. No SI ruler or measured dark-matter abundance is selected.',
        next='Derive a prepared perturbation covariance or formation mechanism, and connect the same-field density and lensing transfer to observers. Check physical scale selection before interpreting this benchmark as small-scale dark matter.',
        process_cpu_seconds=time.process_time()-start,
        source_sha256={x:hashlib.sha256((ROOT/x).read_bytes()).hexdigest() for x in paths},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_parent_perturbations.py')
    if max(gap,packet_clock_gap,packet_rows[-1]['density_difference_from_previous_cutoff'])>1e-7:
        raise RuntimeError('retained perturbation comparisons exceed tolerance')
    if record['process_cpu_seconds']>120:raise RuntimeError('120 CPU-second budget exceeded')
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path)
    args=p.parse_args();encoded=json.dumps(calculation(),indent=2,allow_nan=False)+'\n'
    if len(encoded.encode())>22000:raise RuntimeError('22 KB summary limit exceeded')
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded,end='')
