"""Compact shared scalar slice, all-sector obstruction and protected cosmology."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import sympy as sp
from recursive_horizons.nsc_atlas_shared_scalar_slice import SharedScalarSlice,alternative_floor_test

ROOT=Path(__file__).resolve().parents[1]


def identities():
    R,M,U,kappa,v,vp=sp.symbols('R M U kappa v vp',nonzero=True)
    f=1-2*M/R-kappa*U*R**2/3
    beta2=f+v*v
    Rll=sp.diff(f,R)/2+v*vp
    Lr,Lt=-vp,-v/R
    R3=-4*Rll/R+2*(1-beta2)/R**2
    phi,q,mu,m2=sp.symbols('phi q mu m2',real=True)
    omega2=m2+mu*mu
    potential=(m2*phi*phi+(q+mu*phi)**2)/2
    a,H,Pi,chi,rv=sp.symbols('a H Pi chi rv',nonzero=True)
    E=(Pi*Pi+omega2*chi*chi)/2
    Hd=1/a**2-kappa*Pi*Pi/2
    Pid=-3*H*Pi-omega2*chi
    C=3*(H*H+1/a**2)-kappa*(rv+E)
    cd=sp.diff(C,a)*a*H+sp.diff(C,H)*Hd+sp.diff(C,chi)*Pi+sp.diff(C,Pi)*Pid
    S,F0,t=sp.symbols('S F0 t',real=True)
    b=F0+S/2
    raw={
        'bulk_Hamilton':R3+(Lr+2*Lt)**2-Lr*Lr-2*Lt*Lt-2*kappa*U,
        'bulk_momentum_divided_by_beta':-2*(-vp/R+v/R**2)+2*(Lr-Lt)/R,
        'moving_Misner_Sharp_mass':R*(1+v*v-beta2)/2-M-kappa*U*R**3/6,
        'canonical_ground_floor':potential-omega2*(phi+mu*q/omega2)**2/2-m2*q*q/(2*omega2),
        'parent_Friedmann_propagation':cd,
        'parent_scalar_work_balance':sp.diff(E,chi)*Pi+sp.diff(E,Pi)*Pid+3*H*Pi*Pi,
        'weak_corner_inward_flux':(-b*t)*(-S*t/2)-S*b*t*t/2,
        'weak_corner_pressure_drop':S*(-b*t*t/2)+S*b*t*t/2,
    }
    k,sigma=sp.symbols('k sigma',nonzero=True)
    raw['surface_Hamilton']=(-4*(-k*R)/R-2*kappa*sigma).subs(k,kappa*sigma/2)
    z=sp.symbols('z',real=True)
    shape=z*(1-z*z)**4
    raw['lapse_shape_extremum']=sp.diff(shape,z).subs(z,sp.Rational(1,3))
    raw['lapse_shape_bound']=shape.subs(z,sp.Rational(1,3))-sp.Rational(4096,19683)
    results={name:str(sp.factor(expr)) for name,expr in raw.items()}
    if any(x!='0' for x in results.values()):raise RuntimeError(results)
    return results


def calculation():
    start=time.process_time();s=SharedScalarSlice()
    checks=s.checks();observers=[s.observer(j) for j in (0,1)]
    # A raw squared-gradient residual is more informative than the square
    # root of positive floating-point roundoff around this exact zero state.
    for row in observers:row.pop('normal_magnitude',None)
    parent=[s.evolve_parent(),s.evolve_parent(clock='efolds')]
    comparison={key:parent[0][key]-parent[1][key] for key in
        ('proper_time','H','scalar_displacement','scalar_velocity','excitation_density','conformal_distance')}
    residuals=[checks['middle_Hamilton_max'],checks['middle_momentum_max'],
        abs(checks['independent_radius_endpoint_gap']),abs(checks['independent_volume_gap'])]
    residuals.extend(abs(x) for x in comparison.values())
    residuals.extend(p[name] for p in parent for name in ('max_Friedmann_residual','max_scalar_energy_balance_residual'))
    if max(residuals)>1e-9:raise RuntimeError('declared shared-atlas tolerance exceeded')
    paths=('src/recursive_horizons/nsc_atlas_shared_scalar_slice.py',
        'src/recursive_horizons/nsc_atlas_flux_transition.py',
        'src/recursive_horizons/nsc_atlas_flux_conversion.py',
        'src/recursive_horizons/nsc_atlas_geometric_observer.py',
        'src/recursive_horizons/nsc_atlas_scalar_preparation.py',
        'scripts/check_nsc_atlas_shared_scalar_slice.py',
        'scripts/plot_nsc_atlas_shared_scalar_slice.py',
        'tests/test_nsc_atlas_shared_scalar_slice.py')
    record=dict(schema='nsc-atlas-shared-scalar-slice-v1',
        result='Exact compact three-region canonical scalar/membrane constraint data, all-sector alternative-law obstruction, and smooth protected parent evolution.',
        same_action=dict(kappa=s.scalar.kappa,epsilon=s.scalar.epsilon,omega=s.scalar.omega,
            bare_density=s.scalar.bare_density,charge=s.charge,tension=s.tension),
        initial_data=dict(fluxes_inside_outward=s.fluxes,uniform_scalar=s.field,scalar_Pi=0,scalar_normal=0,
            density=s.densities,ground_densities=[s.scalar.vacuum(q) for q in s.fluxes],
            shared_middle_mass=s.middle_mass,radii=[s.inner_radius,s.outer_radius],
            preparation='Specified uniform scalar at grandchild ground value, with zero canonical momentum. Radii and middle mass are supplied initial data, not derived couplings.'),
        wall_jets=s.walls,cap_geometry=s.cap_geometry,proper_geometry=s.lengths_and_volumes,
        geometric_observers=observers,
        conditional_corner_coefficients=[s.tangent_coefficients(j) for j in (0,1)],
        alternative_law_obstruction=alternative_floor_test(),
        parent_causal_domain=s.parent_causal_domain(),parent_one_efold=parent,
        checks=dict(symbolic_zero_identities=identities(),constraints=checks,
            independent_clock_differences=comparison,comparison_tolerance=1e-9),
        limits='The entire constrained initial atlas and the protected smooth parent domain are constructed. A global weak development of both charged walls, scalar corners, inter-wall transmission, formation, stability and clustering are not constructed. The conditional pulse coefficients use separate local wall clocks. No SI ruler or observed abundance is predicted; Lambda*a0^2 is conditional on the declared action and initial geometry.',
        next='Compute scalar and metric perturbations on the evolving protected parent to test growth, propagation scale and observer transfer from the same inherited source; retain full wall evolution as a separate causal coupling problem.',
        process_cpu_seconds=time.process_time()-start,
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_shared_scalar_slice.py')
    if record['process_cpu_seconds']>120:raise RuntimeError('120 CPU-second budget exceeded')
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path)
    args=p.parse_args();encoded=json.dumps(calculation(),indent=2,allow_nan=False)+'\n'
    if len(encoded.encode())>24000:raise RuntimeError('24 KB summary budget exceeded')
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded,end='')
