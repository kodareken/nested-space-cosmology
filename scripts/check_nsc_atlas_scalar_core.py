"""Self-gravitating scalar profile, finite cooling and sourced atlas admission."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import time

import numpy as np
import sympy as sp
from recursive_horizons.nsc_atlas_scalar_core import solve_core,evolve_core,CoolingRun,CoreAtlasEmbedding
from recursive_horizons.nsc_atlas_shared_scalar_slice import SharedScalarSlice
from recursive_horizons.nsc_atlas_flux_conversion import ScalarFlux

ROOT=Path(__file__).resolve().parents[1]


def identities():
    w,epsilon,q,qg,kap,zeta,r,V,S,Sp,M,Mp,U0,v,vp=sp.symbols('w epsilon q qg kap zeta r V S Sp M Mp U0 v vp',positive=True)
    mu=w*sp.sqrt(1-epsilon);phi0=-sp.sqrt(1-epsilon)*qg/w
    U=epsilon*w*w*phi0**2/2+(q+mu*phi0)**2/2
    beta2=1+v*v-kap*U0*r*r/3-2*M/r
    bp=2*v*vp-2*kap*U0*r/3-2*Mp/r+2*M/r**2
    R3=-2*bp/r+2*(1-beta2)/r**2
    raw={
        'frequency_independent_initial_energy':sp.diff(U,w),
        'initial_energy_completed_square':U-epsilon*q*q/2-(1-epsilon)*(q-qg)**2/2,
        'spherical_Hamilton_mass_identity':R3+4*v*vp/r+2*v*v/r**2-2*kap*U0-4*Mp/r**2,
        'momentum_constraint':-2*(-vp/r+v/r**2)+2*(-vp+v/r)/r,
        'geometric_mass_normalization':kap/4*(w*w*(2/kap)*zeta**4)*r*r/(w*zeta)-zeta*(w*zeta*r)**2/(2*w),
    }
    scale=sp.symbols('scale',positive=True)
    # Fixed-mass dilation of the envelope energy: kinetic, gravity and the
    # cosmological quadratic potential scale differently.
    T,W,L=sp.symbols('T W L')
    raw['virial_scaling']=sp.diff(T*scale**2+W*scale+L/scale**2,scale).subs(scale,1)-(2*T+W-2*L)
    ar,br,arr,brr=sp.symbols('ar br arr brr',real=True)
    aa,bb=sp.symbols('aa bb',real=True)
    raw['real_potential_norm_flux']=2*aa*(-brr/2+V*bb)+2*bb*(arr/2-V*aa)+aa*brr-bb*arr
    chi,delta=sp.symbols('chi delta',real=True)
    raw['inherited_phase_interference']=w*w*((chi+delta)**2-(chi-delta)**2)/2-2*w*w*chi*delta
    answer={k:str(sp.simplify(x)) for k,x in raw.items()}
    if any(x!='0' for x in answer.values()):raise RuntimeError(answer)
    return answer


def worker(item):
    label,config=item
    record,data=evolve_core(config,retain=True)
    digest=hashlib.sha256(np.asarray(data['complex_radial_field'],dtype='<c16').tobytes()).hexdigest()
    record['final_field_sha256']=digest
    # Only in-memory plot projections return to the parent; no raw waveform
    # history or radial array is committed.
    scale=record['mass_selected_scale'];x=np.linspace(0,10,241)
    curve=np.interp(x/scale,data['radius'],data['mean_amplitude'])/scale**2
    plotted=dict(radius=x,profile=curve,snapshots=data['snapshots'])
    return label,record,plotted


def calculation():
    start=time.process_time();begin=time.perf_counter()
    ground=solve_core();gs=ground.summary()
    check=solve_core(tolerance=2e-9,outer=20.).summary()
    original=SharedScalarSlice();controlled=SharedScalarSlice(scalar=ScalarFlux(omega=2000.))
    geometry_gap=max(abs(a-b) for a,b in zip(original.densities+original.expansions,
        controlled.densities+controlled.expansions))
    configs={
        'gaussian_coarse':CoolingRun(points=1023,step=.04),
        'gaussian_fine':CoolingRun(),
        'gaussian_outer':CoolingRun(points=3071,outer=96),
        'supergaussian_fine':CoolingRun(shape='supergaussian'),
        'stationary_control':CoolingRun(shape='stationary',points=1023,step=.04,end=40,h=0.,absorber=False),
        'dispersive_control':CoolingRun(points=1023,step=.04,end=60,h=0.,gravity=False,absorber=False),
    }
    workers=min(6,os.cpu_count() or 1)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        outputs=list(pool.map(worker,configs.items()))
    runs={name:r for name,r,_ in outputs};plots={name:p for name,_,p in outputs}
    baseline=runs['gaussian_fine'];coarse=runs['gaussian_coarse'];outer=runs['gaussian_outer']
    # The coherent representative is fixed by the late core mass, not an
    # independent peak-amplitude fit. Residual breathing is reported above.
    formed_weak_amplitude=.1*baseline['mass_selected_scale']
    embedding=CoreAtlasEmbedding(weak_amplitude=formed_weak_amplitude).solve()
    embedding['representative_from_run']='gaussian_fine'
    embedding['representative_rule']='Nodeless coherent representative with scale fixed by the late time-averaged mass inside envelope radius10; not the full oscillating/radiating snapshot.'
    opposite=CoreAtlasEmbedding(weak_amplitude=formed_weak_amplitude,phase_sign=-1).solve()
    zero=CoreAtlasEmbedding(with_core=False).solve()
    convergence={}
    for name,other in (('coarse_to_fine',coarse),('outer_domain',outer)):
        convergence[name]={key:other[key]-baseline[key] for key in
            ('final_mass','mean_core_mass_inside10','mass_selected_scale','late_profile_gap','outward_flux_at20')}
    omega=2000.;zeta=.1;vacuum=controlled.scalar.vacuum(controlled.parent_flux);Hvac=math.sqrt(vacuum/3)
    physical=dict(frequency=omega,weak_amplitude=zeta,vacuum_density=vacuum,
        inherited_vacuum_H=Hvac,weak_field_order=zeta*zeta,
        vacuum_binding_ratio=Hvac/(omega*zeta*zeta),
        original_frequency_binding_ratio=Hvac/(20*zeta*zeta),
        formed_coherent_weak_amplitude=formed_weak_amplitude,
        formed_coherent_binding_ratio=Hvac/(omega*formed_weak_amplitude**2),
        physical_radius50=gs['mass_fraction_radii']['0.5']/(omega*zeta),
        physical_radius99=gs['mass_fraction_radii']['0.99']/(omega*zeta),
        physical_envelope_time=1/(omega*zeta*zeta),
        initial_core_amplitude=math.sqrt(2)*zeta*zeta*gs['central_amplitude'],
        ambient_parent_displacement=controlled.field-controlled.scalar.minimum(controlled.parent_flux),
        core_density_leading=omega*omega*zeta**4*gs['central_amplitude']**2,
        local_formation_duration=configs['gaussian_fine'].end/(omega*zeta*zeta),
        coupling_status='Frequency2000 is an explicitly different shared coupling, chosen for an approximation hierarchy, not inferred or fitted. Original frequency20 and P003 are not reinterpreted.')
    files=('src/recursive_horizons/nsc_atlas_scalar_core.py',
        'src/recursive_horizons/nsc_atlas_shared_scalar_slice.py',
        'src/recursive_horizons/nsc_atlas_flux_conversion.py','src/recursive_horizons/nsc_atlas_flux_transition.py',
        'scripts/check_nsc_atlas_scalar_core.py','scripts/plot_nsc_atlas_scalar_core.py',
        'tests/test_nsc_atlas_scalar_core.py')
    cpu=time.process_time()-start+sum(r['cpu_seconds'] for r in runs.values())
    record=dict(schema='nsc-atlas-scalar-core-v1',
        result='A calculated self-gravitating real-scalar envelope profile, finite local cooling controls and an exactly sourced compact-atlas initial-constraint embedding.',
        frequency_family=dict(exact_identity='U(phi0,q)=epsilon*q^2/2+(1-epsilon)*(q-q_g)^2/2 with phi0=-sqrt(1-epsilon)*q_g/omega',
            maximum_initial_energy_or_wall_expansion_gap=geometry_gap,
            interpretation='Initial stress and geometry do not select the scalar frequency. The stated frequency member is part of the model input.'),
        stationary_profile=gs,stationary_independent_domain=dict(
            central_difference=check['central_amplitude']-gs['central_amplitude'],mass_difference=check['radial_mass_norm']-gs['radial_mass_norm']),
        physical_hierarchy=physical,
        local_envelope_equations=dict(field='chi=sqrt(2/kappa)*zeta^2*Re[psi(rho,tau)*exp(-i*omega*t)]',
            coordinates='rho=omega*zeta*r; tau=omega*zeta^2*t',
            evolution='i psi_tau=-laplacian(psi)/2+V psi-h^2 rho^2 psi/2; laplacian(V)=|psi|^2/2',
            scope='Leading weak-field slow-envelope equations in a local de Sitter vacuum limit; the complex envelope represents one real scalar, not an added charged field. No primordial seed or full evolving-parent formation is asserted.'),
        cooling_runs=runs,finite_comparisons=convergence,atlas_embedding=embedding,
        opposite_phase_embedding=opposite,
        checks=dict(symbolic_zero_identities=identities(),
            zero_core_parent_length_gap=zero['parent_pole_to_wall_length']-controlled.cap_geometry['parent_length'],
            zero_core_outer_wall_speed_change=zero['outer_wall_speed_change'],
            independent_profile_tolerance=2e-7,initial_constraint_tolerance=1e-7,
            mass_balance_tolerance=2e-9,
            note='Envelope discretization checks and physical residual breathing are distinct from relativistic and embedding-history errors.'),
        limits='Two initial cloud shapes and their mass/phase are supplied. Finite cooling demonstrates radial profile selection for those examples, not a global attractor theorem or an origin of primordial structure. The isolated vacuum-limit cooling evolution and the exact sourced initial atlas are separate constructions; their full dynamical matching, ambient-field interaction, membrane evolution, nonspherical stability and observational identification remain open. No quantum reduction, unique particle mass or generated cosmic Lambda is claimed.',
        next='Develop the interacting core and ambient inherited field on the shared evolving geometry, using the computed cross-term and conserved flux as compatibility data; retain the supplied cloud origin and coupling selection as explicit questions.',
        process_cpu_seconds_including_workers=cpu,wall_seconds=time.perf_counter()-begin,worker_count=workers,
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_scalar_core.py')
    if cpu>240:raise RuntimeError('240 CPU-second formation budget exceeded')
    if max(abs(r['mass_balance']) for r in runs.values())>2e-9:raise RuntimeError('cooling mass accounting failed')
    if max(embedding['maximum_sampled_Hamilton_residual'],opposite['maximum_sampled_Hamilton_residual'])>1e-7:
        raise RuntimeError('sourced initial constraints failed')
    if max(runs[n]['late_profile_gap'] for n in ('gaussian_fine','supergaussian_fine'))>.06:
        raise RuntimeError('finite examples do not approach the declared coherent profile class')
    if max(abs(runs[n]['independent_flux_balance']) for n in configs)>3e-5:
        raise RuntimeError('independent outgoing flux check failed')
    if abs(convergence['outer_domain']['mean_core_mass_inside10'])/baseline['mean_core_mass_inside10']>.005:
        raise RuntimeError('outer-boundary sensitivity too large for the core claim')
    if runs['dispersive_control']['analytic_free_field_error_inside20']>1e-7:
        raise RuntimeError('independent free-field control failed')
    return record,plots,ground


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path)
    parser.add_argument('--figure',action='store_true');args=parser.parse_args()
    result,plots,ground=calculation();encoded=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if len(encoded.encode())>26000:raise RuntimeError('26 KB compact record limit exceeded')
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded,end='')
    if args.figure:
        from plot_nsc_atlas_scalar_core import render
        render(result,plots,ground)
