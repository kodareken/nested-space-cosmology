"""Finite observer-clock transfer and a conditional state-free null relation."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import sympy as sp
from recursive_horizons.nsc_atlas_parent_clocks import ParentClocks,infer_two_amplitudes

ROOT=Path(__file__).resolve().parents[1]


def identities():
    H,p,w,x,lam,a,kap=sp.symbols('H p w x lam a kap',nonzero=True)
    flow=sp.Matrix([[0,1,0],[-w*w-lam/a**2+2*kap*p*p,-3*H,-4*H*p-2*w*w*x],
        [kap*p/2,0,-H]])
    A,B,At,Bt,pt=sp.symbols('A B At Bt pt',nonzero=True)
    # On a radial null curve p_theta^2=B*p_t^2/A. Hamilton equations
    # independently give the logarithmic covariant photon energy rate.
    dpdt=(-pt*pt*At/(2*A*A)+(B*pt*pt/A)*Bt/(2*B*B))/(-pt/A)
    F,Ft=sp.symbols('F Ft')
    ratio=At/(2*A)-Bt/(2*B)
    raw={
        'finite_state_flow_trace':sp.trace(flow)+4*H,
        'null_covariant_energy_transport':dpdt/pt-ratio,
        'perturbed_metric_null_energy':ratio.subs({A:1+2*F,At:2*Ft,B:a*a*(1-2*F),Bt:2*a*a*(H*(1-2*F)-Ft)})+H-2*Ft/(1-4*F*F),
    }
    ae,ao,He,Ho,dto,dx,J,L,Pe,Po,v=sp.symbols('ae ao He Ho dto dx J L Pe Po v')
    dte=ae*(dto/ao+2*J-dx)
    response=He*dte-Ho*dto+Pe-Po-v+2*L
    raw['observer_age_response']=sp.diff(response,dto)-(He*ae/ao-Ho)
    amplitude=sp.symbols('amplitude')
    raw['static_Killing_frequency_sign']=sp.diff(sp.log((1+2*amplitude*Pe)/(1+2*amplitude*Po))/2,amplitude).subs(amplitude,0)-Pe+Po
    raw['recession_Doppler_sign']=sp.diff(-sp.atanh(amplitude*v),amplitude).subs(amplitude,0)+v
    u,v1,w1,z=sp.symbols('u v1 w1 z',nonzero=True)
    r,s=sp.symbols('r s')
    M=sp.Matrix([[u,v1],[w1,z]]);third=sp.Matrix([[r,s]])
    weights=-(third*M.inv())
    raw['first_state_eliminated']=(weights*M+third)[0]
    raw['second_state_eliminated']=(weights*M+third)[1]
    answers={k:str(sp.factor(y)) for k,y in raw.items()}
    if any(v!='0' for v in answers.values()):raise RuntimeError(answers)
    return answers


def calculation():
    start=time.process_time();m=ParentClocks();positions=(.04,.08,.12)
    rows=[m.response(x) for x in positions]
    R=np.array([r['log_frequency_response'] for r in rows]);inverse=np.linalg.inv(R[:2])
    weights=np.r_[-R[2]@inverse,1.]
    fine=ParentClocks(cutoff=256,rtol=2e-11,atol=2e-13)
    refined=[fine.response(x) for x in positions];R2=np.array([r['log_frequency_response'] for r in refined])
    w2=np.r_[-R2[2]@np.linalg.inv(R2[:2]),1.]
    zero=m.finite_metric_clock(.12,[0,0])
    zero_gap=zero['log_frequency_ratio']-math.log(rows[-1]['background_frequency_ratio'])
    comparisons=[(.12,[1.,0.],(.001,.0005)),(.12,[0.,1.],(.001,.0005)),
        (.12,[1.,-.4],(.001,.0005)),(.04,[1.,-.4],(.0005,)),(.08,[1.,-.4],(.0005,))]
    geodesics=[]
    for position,direction,steps in comparisons:
        direction=np.array(direction);target=float(R[positions.index(position)]@direction);derivatives=[]
        for step in steps:
            plus=m.finite_metric_clock(position,step*direction)
            minus=m.finite_metric_clock(position,-step*direction)
            measured=(plus['log_frequency_ratio']-minus['log_frequency_ratio'])/(2*step)
            derivatives.append(dict(step=step,value=measured,error=measured-target))
        geodesics.append(dict(initial_angle=position,direction=direction.tolist(),linear=target,central_differences=derivatives))
    truth=np.array([1e-5,-4e-6]);y=R@truth;recovered=infer_two_amplitudes(rows,y[:2])
    sigma=1e-10;noise_gain=float(np.linalg.norm(weights))
    noise_amplitude=sigma*np.sqrt(np.diag(inverse@inverse.T))
    corrupt=y.copy();corrupt[2]+=2.5e-9
    omitted=[]
    bo=m.background(m.to)
    for row in rows:
        be=m.background(row['background_emission_time'])
        omitted.append((np.array(row['log_frequency_response'])-
            np.array(row['receiver_time_shift'])*(be['H']*be['a']/bo['a']-bo['H'])).tolist())
    pars=m.packet.model.atlas
    prediction=dict(observable='y_i=log[(nu_received/nu_emitted)_i/(a_e/a_o)_i]; first-order log frequency departure',
        reference_length='Lref=a0, initial parent spatial curvature radius; no SI value selected',
        initial_proper_separations_over_Lref=list(positions),
        receiver_proper_age_over_Lref=m.to/m.initial['a'],
        emitter_rule='Radial freely falling test clocks, initially at rest relative to the Cauchy normal, at fixed initial proper separations.',
        receiver_rule='Central freely falling receiver, fixed proper age from the shared initial slice.',
        preparation_family='Existing compact support-angle1.4 potential, power8 compensated profile; independent density and phase amplitudes, fixed shape and couplings.',
        response_matrix=R.tolist(),null_weights=weights.tolist(),
        null_relation='sum_i null_weights_i*y_i=0 at retained linear order; third row withheld after first two calibrate the state',
        coefficient_comparison_tolerance=1e-7,
        observational_validation=False)
    files=('src/recursive_horizons/nsc_atlas_parent_clocks.py',
        'src/recursive_horizons/nsc_atlas_parent_perturbations.py',
        'src/recursive_horizons/nsc_atlas_shared_scalar_slice.py',
        'scripts/check_nsc_atlas_parent_clocks.py','scripts/plot_nsc_atlas_parent_clocks.py',
        'tests/test_nsc_atlas_parent_clocks.py')
    record=dict(schema='nsc-atlas-parent-clocks-v1',
        result='Operational freely falling clock response in the inherited scalar parent; a three-clock relation eliminates both free state amplitudes within the fixed compact-profile family.',
        same_action=dict(kappa=pars.scalar.kappa,epsilon=pars.scalar.epsilon,omega=pars.scalar.omega,
            charge=pars.charge,tension=pars.tension,bare_density=pars.scalar.bare_density,parent_flux=pars.parent_flux),
        state_nonselection=dict(flow_trace='-4H',full_flow_determinant='exp(-4N)>0 for every finite N',
            determinant_at_one_efold=math.exp(-4),
            scope='Finite deterministic evolution maps the two-dimensional constraint space bijectively. It does not select a unique phase or covariance. An asymptotic selection rule, a boundary condition, open-system tracing or a different formation law is a separate physical input.'),
        prediction=prediction,clock_rows=rows,
        calibration=dict(matrix_condition_number=float(np.linalg.cond(R[:2])),
            two_clock_singular_values=np.linalg.svd(R[:2],compute_uv=False).tolist(),
            three_clock_singular_values=np.linalg.svd(R,compute_uv=False).tolist(),
            determinant=float(np.linalg.det(R[:2])),
            basis='alpha is initial central density contrast; beta is the existing unit phase-column amplitude',
            scope='Conditioning depends on the specified amplitude normalization; noise propagation is reported in that same basis.'),
        synthetic_test=dict(scope='Generated check only, not astronomical data, fitted fundamental parameters or observational confirmation.',
            true_amplitudes=truth.tolist(),log_frequency_departures=y.tolist(),
            recovered_from_first_two=recovered.tolist(),predicted_third=float(R[2]@recovered),
            withheld_error=float(R[2]@recovered-y[2]),
            hypothetical_independent_clock_sigma=sigma,amplitude_standard_deviations=noise_amplitude.tolist(),
            null_residual_noise_gain=noise_gain,
            corrupted_third_shift=2.5e-9,corrupted_null_residual=float(weights@corrupt),
            corrupted_standardized_residual=float(weights@corrupt)/(sigma*noise_gain)),
        independent_geodesic_check=dict(zero_source_log_frequency_gap=zero_gap,variations=geodesics,
            scope='Exact radial geodesics and clock/distance integrals in a prescribed finite metric deformation validate its first variation. The finite metric is not asserted to solve the nonlinear Einstein-scalar equations.'),
        coordinate_age_counterfactual=dict(response_matrix=omitted,
            meaning='Dropping the receiver proper-age shift changes the events and these responses; it is not an equivalent observable.'),
        checks=dict(symbolic_zero_identities=identities(),maximum_refined_response_gap=float(np.max(abs(R2-R))),
            maximum_refined_null_weight_gap=float(np.max(abs(w2-weights))),
            maximum_null_relation_coefficient_residual=float(np.max(abs(weights@R)))),
        limits='The compact shape, background, initial proper source separations and clock rule are fixed inputs. Arbitrary spatial profiles, unknown source peculiar velocities or additional matter add independent degrees of freedom and need not obey this three-clock closure. No empirical data, SI ruler, primordial spectrum or nonlinear remainder bound is supplied.',
        next='Use additional operational observables to test the declared shape and source-motion assumptions, and connect a physical formation mechanism to the remaining initial-state family before a cosmological inference.',
        process_cpu_seconds=time.process_time()-start,
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_parent_clocks.py')
    all_errors=[abs(d['error']) for row in geodesics for d in row['central_differences']]
    if max(all_errors+[abs(zero_gap),record['checks']['maximum_refined_response_gap']])>1e-8:
        raise RuntimeError('clock response comparison exceeded 1e-8')
    if record['checks']['maximum_refined_null_weight_gap']>1e-7:raise RuntimeError('null weights unresolved')
    if record['process_cpu_seconds']>120:raise RuntimeError('120 CPU-second budget exceeded')
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path)
    args=p.parse_args();encoded=json.dumps(calculation(),indent=2,allow_nan=False)+'\n'
    if len(encoded.encode())>22000:raise RuntimeError('22 KB summary budget exceeded')
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded,end='')
