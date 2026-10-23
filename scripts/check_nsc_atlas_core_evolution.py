"""One sourced evolution, a diffuse concentration test and a received clock."""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_core_evolution import EvolutionRun,run_evolution
from recursive_horizons.nsc_atlas_core_reservoir import CoreReservoir

ROOT=Path(__file__).resolve().parents[1]


def identities():
    R,H,Hd,al,ar,m,Pi,F,U,k=sp.symbols('R H Hd al ar m Pi F U k',real=True)
    b=sp.symbols('b',positive=True);v=H*R
    rho=(Pi**2+b*b*F*F)/2+U;pr=rho-2*U
    mr=k*R*R*(rho-v*Pi*F)/2
    mt=al*k*R*R*((b*b+v*v)*Pi*F-v*(Pi*Pi+b*b*F*F))/2
    br=2*H*H*R+2*m/R**2-2*mr/R
    bt=2*H*Hd*R**2-2*mt/R
    hd=(b*b*ar-al*(H*H*R+m/R**2+k*R*pr/2))/R
    C=al*(v+b);Cr=ar*(v+b)+al*(H+br/(2*b))
    normal=al*(-H-k*R*Pi*F/2)-b*ar
    raw={'canonical_and_normal_photon':(-Cr+(bt+C*br)/(2*b*b)-normal).subs(Hd,hd)}
    db,dp,bp,cb,dc,M=sp.symbols('db dp bp cb dc M',real=True)
    dpress=((bp+dp)**2-bp**2+b*b*F*F)/2-(2*cb*dc+dc**2)/2
    Q=M/R**2+k*R*dpress/2
    raw['subtracted_lapse']=(1-al)*Hd*R+al*Q-((Q-Hd*R)*(al-1)+Q)
    fixed=mt.subs(Pi,bp+dp)+k*H*R**3*bp*bp/2
    subtracted=al*k*R*R*((b*b+v*v)*(bp+dp)*F-v*(dp*(dp+2*bp)+b*b*F*F))/2-(al-1)*k*H*R**3*bp*bp/2
    raw['subtracted_time_mass_flux']=fixed-subtracted
    out={name:str(sp.simplify(expr)) for name,expr in raw.items()}
    if any(v!='0' for v in out.values()):raise RuntimeError(out)
    return out


def independent_homogeneous(emission_age,receiver=.6):
    """Unaveraged DOP853 cosmology and analytic FLRW null relation."""
    model=CoreReservoir();s=model.reference;a0=s.cap_geometry['parent_maximum_radius']
    def rhs(t,y):
        chi,p,N,H,eta,tau=y
        rates=s.parent_rates(chi,p,N,H);a=a0*math.exp(N);B=1-receiver**2/a**2
        return (*rates[:4],rates[5],math.sqrt(1-(H*receiver)**2/B))
    end=emission_age+1
    sol=solve_ivp(rhs,(0,end),(model.ambient,0.,0.,model.H,0.,0.),method='DOP853',
        rtol=2e-12,atol=(2e-15,2e-12,2e-13,2e-13,2e-14,2e-13),dense_output=True,max_step=.002)
    if not sol.success:raise RuntimeError(sol.message)
    e=sol.sol(emission_age)
    def radius(t):
        y=sol.sol(t);return a0*math.exp(y[2])*math.sin(y[4]-e[4])
    t=brentq(lambda t:radius(t)-receiver,emission_age,end,xtol=1e-13)
    y=sol.sol(t);a=a0*math.exp(y[2]);beta=math.sqrt(1-receiver**2/a**2);v=y[3]*receiver
    ratio=math.exp(e[2]-y[2])*math.sqrt((beta+v)/(beta-v))
    return dict(coordinate_time=t,receiver_proper_age=float(y[5]),received_to_emitted_frequency=ratio,
                rhs_calls=sol.nfev,emission_proper_age=emission_age)


def worker(item):
    label,c=item;r,d=run_evolution(c,retain=True)
    state=np.column_stack([d['radius'],d['delta'],d['Pi']])
    r['final_state_sha256']=hashlib.sha256(np.asarray(state,dtype='<f8').tobytes()).hexdigest()
    history=d['history']
    # Histories are transient plotting inputs; compact checkpoints are retained.
    times=(.5,1.,2.,3.,3.5,4.,4.5,5.)
    r['checkpoints']=[dict(time=float(history[j,0]),central_excitation=float(history[j,11]),
        mass_inside_025=float(history[j,12]),half_mass_radius=float(history[j,-1]) if np.isfinite(history[j,-1]) else None)
        for t in times if t<=c.end for j in [int(np.argmin(abs(history[:,0]-t)))]]
    return label,r,history


def calculation():
    start=time.process_time();wall=time.perf_counter()
    base=EvolutionRun();diffuse=EvolutionRun(kind='diffuse',end=5.,emission_age=4.)
    configs={
        'diffuse_fine':replace(diffuse,spacing=.00025,step=.000025),
        'diffuse_medium':replace(diffuse,spacing=.0005,step=.00005),
        'diffuse_coarse':diffuse,
        'diffuse_time_control':replace(diffuse,spacing=.0005,step=.000025),
        'diffuse_outer_control':replace(diffuse,outer=4.5,spacing=.0005,step=.00005),
        'fixed_geometry_control':replace(diffuse,spacing=.0005,step=.00005,backreaction=False),
        'core_fine':replace(base,spacing=.00025,step=.000025),
        'core_medium':replace(base,spacing=.0005,step=.00005),
        'core_coarse':base,
        'homogeneous_early':replace(base,kind='homogeneous',spacing=.002),
        'homogeneous_late':replace(diffuse,kind='homogeneous',spacing=.002),
    }
    with ProcessPoolExecutor(max_workers=6) as pool:
        outputs=list(pool.map(worker,configs.items()))
    runs={label:r for label,r,_ in outputs};plots={label:h for label,_,h in outputs}
    backgrounds={key:independent_homogeneous(t) for key,t in (('early',0.),('late',4.))}
    predictions={}
    for family,background in (('core','early'),('diffuse','late')):
        ref=backgrounds[background];fine=runs[family+'_fine'];medium=runs[family+'_medium']
        a=fine['arrival'];factor=a['received_to_emitted_frequency']/ref['received_to_emitted_frequency']
        sigma=2*max(abs(a['received_to_emitted_frequency']-medium['arrival']['received_to_emitted_frequency']),
            abs(fine['geodesic_frequency_gap']),
            abs(runs['homogeneous_'+background]['arrival']['received_to_emitted_frequency']-ref['received_to_emitted_frequency']))
        delay=a['receiver_proper_age']-ref['receiver_proper_age']
        delay_error=2*max(abs(a['receiver_proper_age']-medium['arrival']['receiver_proper_age']),
            abs(runs['homogeneous_'+background]['arrival']['receiver_proper_age']-ref['receiver_proper_age']))
        if family=='diffuse':
            for control in ('diffuse_time_control','diffuse_outer_control'):
                sigma=max(sigma,2*abs(runs[control]['arrival']['received_to_emitted_frequency']-medium['arrival']['received_to_emitted_frequency']))
                delay_error=max(delay_error,2*abs(runs[control]['arrival']['receiver_proper_age']-medium['arrival']['receiver_proper_age']))
        predictions[family]=dict(frequency_fractional_departure=factor-1,
            frequency_log_departure=math.log(factor),frequency_ratio_numerical_error=sigma,
            fractional_departure_numerical_error=sigma/ref['received_to_emitted_frequency'],
            receiver_proper_delay=delay,proper_delay_numerical_error=delay_error,
            interpretation='Same initial slice, central proper emission age and held receiver areal radius; error is a finite comparison envelope, not a rigorous error bound or physical/model uncertainty.')
    fine=runs['diffuse_fine'];medium=runs['diffuse_medium'];control=runs['fixed_geometry_control']
    formation=dict(initial_core_mass=fine['initial_mass_at_cuts'][2],
        initial_half_mass_radius=fine['initial_half_mass_radius'],late_half_mass_radius=fine['late_mean_half_mass_radius'],
        radius_ratio=fine['late_mean_half_mass_radius']/fine['initial_half_mass_radius'],
        initial_central_excitation=fine['initial_central_excitation'],late_central_excitation=fine['late_mean_central_excitation'],
        central_excitation_gain=fine['late_mean_central_excitation']/fine['initial_central_excitation'],
        fixed_geometry_late_excitation=control['late_mean_central_excitation'],
        gravity_to_fixed_geometry_ratio=fine['late_mean_central_excitation']/control['late_mean_central_excitation'],
        radius_medium_fine_spread=abs(fine['late_mean_half_mass_radius']-medium['late_mean_half_mass_radius']),
        excitation_medium_fine_spread=abs(fine['late_mean_central_excitation']-medium['late_mean_central_excitation']),
        scope='Finite nonlinear concentration followed by breathing. Neither an asymptotic stationary remnant nor universal formation from arbitrary inherited data is asserted.')
    cpu=time.process_time()-start+sum(r['cpu_seconds'] for r in runs.values())
    files=('src/recursive_horizons/nsc_atlas_core_evolution.py','src/recursive_horizons/nsc_atlas_core_reservoir.py',
        'src/recursive_horizons/nsc_atlas_scalar_core.py','src/recursive_horizons/nsc_atlas_shared_scalar_slice.py',
        'src/recursive_horizons/nsc_atlas_flux_conversion.py','scripts/check_nsc_atlas_core_evolution.py',
        'scripts/plot_nsc_atlas_core_evolution.py','tests/test_nsc_atlas_core_evolution.py',
        'results/summaries/nsc-atlas-core-reservoir-v1.json')
    record=dict(schema='nsc-atlas-core-evolution-v1',
        result='Full finite spherical Einstein-scalar evolution connects a compensated inherited preparation, nonlinear concentration, mass/flux conservation and operational clock signals.',
        symbolic_zero_identities=identities(),runs=runs,formation=formation,observer_predictions=predictions,
        independent_homogeneous=backgrounds,
        units=dict(kappa=1,omega=2000,reference_length='a0=6.62607707610 atlas length units, fixed by initial parent geometry',
            receiver_radius=.6,receiver_radius_over_reference=.6/CoreReservoir().reference.cap_geometry['parent_maximum_radius'],
            vacuum_Lambda_reference_squared=12.3482523989285,
            physical_scale='No SI length or observed scalar mass is selected. Frequency departures are dimensionless; delays scale as the declared atlas length divided by c.'),
        protocol='Central geodesic emits at proper age0 for the core control and4 for the diffuse formation run. Receiver is actively held at areal radius.6, with its clock set to zero on the initial slice. Arrival age is an output, not a fixed coordinate-time comparison.',
        method='Unaveraged constrained Einstein-Klein-Gordon; eighth-order spatial derivatives, fourth-order radial constraints and Lawson RK4 exact frequency rotation. Smooth outer moving mesh is fixed belowR=1 and comoving above initial radius2; outer initial radius4 lies in the protected reference annulus.',
        limits='Fixed frequency/action, supplied initial shape/phase/pulse and compensated mass. Spherical finite-time evidence; no complete membrane future, nonspherical stability, thermal history, primordial abundance or measured cosmological identification. The finite breathing concentration is not an asymptotic formation theorem.',
        next='Test how inherited preparation selects the cloud mass/phase and population, and whether the formed concentrations persist under nonradial disturbances; compare multiple late proper-age signals without retuning the initial state.',
        worker_count=6,process_cpu_seconds_including_workers=cpu,wall_seconds=time.perf_counter()-wall,
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_core_evolution.py --figure')
    if cpu>3600:raise RuntimeError('3600 CPU-second evolution campaign budget exceeded')
    for family in ('core','diffuse'):
        errors=[max(runs[family+'_'+level]['maximum_mass_flux_error']) for level in ('coarse','medium','fine')]
        if not errors[0]>errors[1]>errors[2]:raise RuntimeError('mass-flux refinement does not improve: '+str(errors))
        if errors[2]>3e-9:raise RuntimeError('fine mass-flux mismatch exceeds the admitted threshold')
        if predictions[family]['frequency_ratio_numerical_error']>2e-6:raise RuntimeError('clock frequency error unresolved')
        if abs(runs[family+'_fine']['geodesic_frequency_gap'])>1e-7:raise RuntimeError('independent photon equations disagree')
    if formation['radius_medium_fine_spread']>5e-4:raise RuntimeError('concentration radius unresolved')
    if max(r['maximum_background_constraint'] for r in runs.values())>1e-8:raise RuntimeError('background Einstein constraint unresolved')
    return record,plots


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);p.add_argument('--figure',action='store_true');args=p.parse_args()
    result,plots=calculation();encoded=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if len(encoded.encode())>65000:raise RuntimeError('65KB compact record limit exceeded')
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded,end='')
    if args.figure:
        from plot_nsc_atlas_core_evolution import render
        render(result,plots)
