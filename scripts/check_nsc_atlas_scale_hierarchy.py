"""Fixed-vacuum core hierarchy, exact supply bound and solved reservoir."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import sympy as sp
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_scale_hierarchy import ScaleHierarchy, depletion_integral
from recursive_horizons.nsc_atlas_core_reservoir import CoreReservoir

ROOT=Path(__file__).resolve().parents[1]


def identities():
    x,w,z,k,d,c=sp.symbols('x w z k d c',positive=True)
    N,A,Ap=sp.symbols('N A Ap',real=True)
    R=x/(w*z);M=z*N/w;f=sp.sqrt(2/k)*z*z*A;F=sp.sqrt(2/k)*w*z**3*Ap
    B=1-c*R*R-2*M/R
    original=k*R*R*(B*F*F/4+w*w*(d/w)*f/2+w*w*f*f/4)/z**2
    polynomial=-z**4*x*Ap**2*N+x*x*(A*A+z*z*Ap*Ap)/2+x*x*sp.sqrt(k/2)*d*A/(w*z*z)-c*x**4*Ap*Ap/(2*w*w)
    b=1024*x**5*(1-x)**5;a=sp.symbols('a',positive=True)
    J=sp.integrate(x*x*(2*a*b-a*a*b*b),(x,0,1))
    fixed=sp.factor(J.subs(a,sp.Rational(1,2)))
    P,v,F,B,U=sp.symbols('P v F B U',real=True)
    square=(P*P+B*F*F)/2+U-v*P*F-(U+(P-v*F)**2/2+(B-v*v)*F*F/2)
    return dict(core_constraint_polynomial_residual=str(sp.simplify(original-polynomial)),
        untrapped_supply_square_identity=str(sp.simplify(square)),
        weighted_depletion_integral=str(sp.factor(J)), half_depletion_integral=str(fixed),
        asymptotic_vacuum_radius_mass_ratio=str(sp.Rational(3,2)/fixed))


def worker(job):
    name,w,z,a,rtol=job;t=time.process_time()
    result=ScaleHierarchy(w,z,a).solve(rtol=rtol)
    result['cpu_seconds']=time.process_time()-t
    return name,result


def dimensional_example(row, vacuum_per_m2=1e-52):
    # A declared dimensional reference, not a fitted or predicted Lambda.
    c=299792458.;au=149597870700.;pc=648000/math.pi*au;gm_sun=1.3271244e20
    L=math.sqrt(row['parent_vacuum']/vacuum_per_m2)
    return dict(assumed_Lambda_per_m2=vacuum_per_m2, length_unit_metres=L,
        core_geometric_mass_in_nominal_solar_units=c*c*L*row['core']['mass']/gm_sun,
        core_half_mass_radius_kpc=L*row['core']['half_mass_radius']/(1000*pc),
        core_cutoff_radius_kpc=L*row['core']['cutoff_radius']/(1000*pc),
        reservoir_outer_radius_kpc=L*row['reservoir_outer_radius']/(1000*pc),
        interpretation='Unit conversion of an explicitly supplied omega=1e9, zeta=.001, depletion=.5 initial preparation at the illustrative Lambda=1e-52 m^-2. No observation was fitted, and no actual galaxy or stable halo is identified.',
        constants=dict(c_m_per_s=c,au_metres=au,parsec_metres=pc,nominal_solar_GM_m3_per_s2=gm_sun),
        constants_source='https://www.iau.org/static/resolutions/IAU2015_English.pdf')


def calculation():
    wall=time.perf_counter();cpu=time.process_time();z=.08266217132541569
    jobs=[(f'frequency_{w:g}',w,z,.5,2e-11) for w in (2000.,8000.,32000.,128000.,512000.,2e6,1e8,1e10)]
    jobs += [('small_compactness',1e9,.001,.5,2e-11),('small_compactness_tight',1e9,.001,.5,2e-13),
        ('baseline_tight',2000.,z,.5,2e-13),('depletion_quarter',1e9,.001,.25,2e-11),
        ('depletion_three_quarters',1e9,.001,.75,2e-11)]
    with ProcessPoolExecutor(max_workers=6) as pool:runs=dict(pool.map(worker,jobs))
    symbolic=identities();base=ScaleHierarchy();data,_=base.core_data()
    coeff=np.array(data['coefficients']);k=base.reference.scalar.kappa
    # Necessary ideal supply bound at Rd=2Rc; equality is a quadratic
    # in omega because the exact core mass has three inverse powers.
    rhs=k*base.excitation*(2**3-1)*20**3/(6*z**4)
    threshold=float(np.max(np.roots([coeff[0],coeff[1],coeff[2]-rhs])))
    old=CoreReservoir(pulse_fraction=0.)
    old_solution=old.solve()
    new_at_old,_=base.reservoir(old.reservoir_outer_factor*base.radius,data['mass'],depletion=old_solution.depletion)
    independent_core=old.core.mass_solution(old.H)[0].y[0,-1]
    comparisons=dict(original_core_mass_gap=data['mass']-independent_core,
        original_depletion=old_solution.depletion,
        original_no_pulse_mass_endpoint=float(old_solution.solution.y[0,-1]),
        scaled_reservoir_at_original_depletion=new_at_old['endpoint_mass_ratio'])
    refinement={}
    for name,control in (('small_compactness','small_compactness_tight'),('frequency_2000','baseline_tight')):
        r,t=runs[name],runs[control]
        refinement[name]=dict(core_mass_relative=abs(r['core']['mass']/t['core']['mass']-1),
            reservoir_radius_relative=abs(r['reservoir_outer_radius']/t['reservoir_outer_radius']-1),
            half_mass_radius_relative=abs(r['core']['half_mass_radius']/t['core']['half_mass_radius']-1))
    errors=[]
    if symbolic['core_constraint_polynomial_residual']!='0' or symbolic['untrapped_supply_square_identity']!='0':
        errors.append('symbolic constraint identity failed')
    if abs(comparisons['original_core_mass_gap'])>2e-12:errors.append('original core comparison failed')
    if abs(comparisons['scaled_reservoir_at_original_depletion'])>1e-7:errors.append('original reservoir comparison failed')
    for name,r in runs.items():
        if abs(r['reservoir_constraint']['endpoint_mass_ratio'])>1e-9:errors.append(name+': mass compensation unresolved')
        if abs(r['independent_proper_distance']['endpoint_mass_ratio'])>1e-8:errors.append(name+': independent Hamilton check failed')
        if min(r['core']['minimum_core_beta_squared'],r['reservoir_constraint']['minimum_beta_squared'])<=.98:
            errors.append(name+': positive chart margin failed')
    if max(v for row in refinement.values() for v in row.values())>1e-7:errors.append('tolerance refinement failed')
    files=('src/recursive_horizons/nsc_atlas_scale_hierarchy.py','scripts/check_nsc_atlas_scale_hierarchy.py',
        'src/recursive_horizons/nsc_atlas_scalar_core.py','src/recursive_horizons/nsc_atlas_core_reservoir.py',
        'src/recursive_horizons/nsc_atlas_shared_scalar_slice.py','src/recursive_horizons/nsc_atlas_flux_conversion.py',
        'src/recursive_horizons/nsc_atlas_flux_transition.py','scripts/plot_nsc_atlas_scale_hierarchy.py',
        'tests/test_nsc_atlas_scale_hierarchy.py')
    total=time.process_time()-cpu+sum(r['cpu_seconds'] for r in runs.values())
    return dict(schema='nsc-atlas-scale-hierarchy-v1',
        result='An exact initial-constraint family separates core and supply radii while retaining inherited parent vacuum and original initial wall data.',
        symbolic_identities=symbolic,runs=runs,original_implementation_checks=comparisons,
        fixed_ratio_obstruction=dict(reservoir_to_core_cutoff=2,necessary_frequency_ceiling=threshold,
            scope='Above this frequency, even complete removal of available excitation throughout the annulus cannot compensate the supplied positive core. The bound permits arbitrary depletion profiles with Pi=0 and U>=Uvac; gradient energy only strengthens it.'),
        refinement=refinement,illustrative_physical_conversion=dimensional_example(runs['small_compactness']),
        numerical_validation_errors=errors,process_cpu_seconds_including_workers=total,wall_seconds=time.perf_counter()-wall,
        compute_budget_seconds=120,budget_exhausted=total>120,
        limits='Initial Einstein-scalar constraints and conditional asymptotic hierarchy, not a high-frequency formation/evolution or lens prediction. Frequency, weak amplitude, reservoir bump and depletion are supplied. The ratio18.330588 is specific to half-depletion and the declared bump, not a universal constant or an observed dark-matter fraction.',
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_scale_hierarchy.py --output results/summaries/nsc-atlas-scale-hierarchy-v1.json')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path)
    p.add_argument('--figure',action='store_true');args=p.parse_args()
    result=calculation();text=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if args.output:args.output.write_text(text,encoding='utf-8',newline='\n')
    else:print(text)
    print(json.dumps(dict(errors=result['numerical_validation_errors'],cpu=result['process_cpu_seconds_including_workers'],
        obstruction=result['fixed_ratio_obstruction'],example=result['illustrative_physical_conversion'])))
    if result['numerical_validation_errors']:raise RuntimeError(result['numerical_validation_errors'])
    if result['budget_exhausted']:raise RuntimeError('120 CPU-second budget exhausted; diagnose before more calculations')
    if args.figure:
        from plot_nsc_atlas_scale_hierarchy import render
        render(result)
