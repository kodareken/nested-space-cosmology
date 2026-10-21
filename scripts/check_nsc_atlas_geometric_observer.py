"""Exact wall observer identities, coupled reference, and synthetic calibration."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import time

from scipy.integrate import quad
import sympy as sp
from recursive_horizons.nsc_atlas_geometric_observer import (
    ScalarWallObserver,regular_wall_geometry,regular_bulk_Ruu,
    equivalent_one_wall_example,infer_mixing_and_charge)

ROOT=Path(__file__).resolve().parents[1]


def identities():
    R,v,d,A,B,Ad,Bd,kap=sp.symbols('R v d A B Ad Bd kap',nonzero=True)
    m=R*(1+v*v-B*B*R*R)/2
    md=sp.diff(m,R)*v+sp.diff(m,v)*d+sp.diff(m,B)*Bd
    pn=(-2*d/R-(1+v*v)/R**2+2*A*B+B*B)/kap
    F=2*(v/R*(A-B)-Bd)/kap
    Pi,Phi,U,S,omega,Pc,Lambda,chi=sp.symbols('Pi Phi U S omega Pc Lambda chi',nonzero=True)
    p=(Pi*Pi+Phi*Phi)/2-U
    p1,p2,p3,C=sp.symbols('p1 p2 p3 C',nonzero=True)
    Y=lambda P:-Lambda-C*(P-Pc)**2
    dd=sum(Y(a)/((a-b)*(a-c)) for a,b,c in ((p1,p2,p3),(p2,p3,p1),(p3,p1,p2)))
    eps,e,q,bare=sp.symbols('eps e q bare',nonzero=True)
    vac=lambda f:kap*(bare+eps*f*f/2)
    spacing=vac(q+e)-2*vac(q)+vac(q-e)
    raw={
        'finite_time_mass_balance':md-kap*R*R*(B*R*F-v*pn)/2,
        'vacuum_subtracted_mass_balance':md-Lambda*R*R*v/2-kap*R*R*(B*R*F-v*(pn+Lambda/kap))/2,
        'canonical_scalar_stress_identity':(Pi*Phi)**2-Pi*Pi*(2*(p+U)-Pi*Pi),
        'independent_tidal_parabola':kap*(Pi*Pi-Lambda/kap-omega**2*chi*chi/2)-kap*(S*Pi)**2/S**2
             +Lambda+kap*omega**2*(S*chi)**2/(2*S*S),
        'source_free_second_divided_difference':dd+C,
        'equal_charge_ground_spacing':spacing-kap*eps*e*e,
        'local_tidal_and_inherited_vacuum_relation':kap/(2*(1-eps)*e*e)*spacing/kap**2-eps/(2*(1-eps)),
    }
    phi=sp.symbols('phi',real=True)
    UA=lambda f,g:50*f*f+(g+10*sp.sqrt(3)*f)**2/2
    UB=lambda f,g:sp.Rational(1250,7)*f*f+(g+10*sp.sqrt(sp.Rational(3,7))*f)**2/2
    translation=3*sp.sqrt(3)/100
    diff=sp.expand(UB(phi+translation,sp.sqrt(7)*(q-sp.Rational(6,5)))-UA(phi,q))
    raw['two_sector_potential_equivalence']=diff-3*(q-1)*(q-sp.Rational(3,2))
    x,z=kap*eps*e*e/kap,kap/(2*(kap/(2*(1-eps)*e*e)))
    raw['inheritance_charge_inverse']=x+z-e*e
    raw['inheritance_mixing_inverse']=x/(x+z)-eps
    results={k:str(sp.factor(x)) for k,x in raw.items()}
    if any(x!='0' for x in results.values()):raise RuntimeError(results)
    return results


def synthetic_samples(observer):
    """Algebraic field/curvature samples only; not an evolved spacetime."""
    m=observer.prepared
    rows=[]
    for P,Pi in ((.20,.03),(.22,.04),(.24,.02)):
        chi=(P-observer.ground_pressure)/observer.force_scale
        U=m.child_vacuum+m.scalar.omega**2*chi*chi/2
        rows.append(dict(pressure=P,pressure_dot=observer.force_scale*Pi,
                         Ruu=m.scalar.kappa*(Pi*Pi-U)))
    return rows


def calculation():
    start=time.process_time()
    observer=ScalarWallObserver()
    model=observer.prepared
    sol,meta=model.solve(tolerance=1e-9)
    curvature=regular_bulk_Ruu(sol)
    actual_phi=model.scalar.minimum(model.child_flux)+float(sol.y[2,-1])
    actual_normal=float(sol.y[3,-1])
    actual_U=model.scalar.potential(actual_phi,model.child_flux)
    actual_p=actual_normal**2/2-actual_U
    rows=[]
    for eta in (0.,.5,1.):
        geo=regular_wall_geometry(model,sol,eta)
        row=observer.reconstruct(geo,measured_Ruu=curvature)
        row.update(eta=eta,proper_time=float(sol.y[0,-1])*eta,
            scalar_value_error=row['scalar_value']-actual_phi,
            normal_magnitude_error=row['normal_magnitude']-abs(actual_normal),
            direct_field_pressure_error=row['radial_pressure']-actual_p,
            direct_field_mass_derivative_error=row['geometric_mass_dot']+model.scalar.kappa*geo.radius**2*geo.radius_dot*actual_p/2)
        rows.append(row)
    R0=float(sol.y[0,-1])
    direct_work=quad(lambda t:-model.scalar.kappa*(R0*math.cosh(t/R0))**2
                     *math.sinh(t/R0)*(actual_p+model.child_vacuum)/2,0,R0,epsabs=1e-13)[0]
    mass_increase=rows[-1]['vacuum_subtracted_mass']-rows[0]['vacuum_subtracted_mass']
    samples=synthetic_samples(observer)
    calibrated=observer.calibrate_tidal_pair(samples[0],samples[1])
    prediction=observer.predict_tidal(samples[2]['pressure'],samples[2]['pressure_dot'],calibrated)
    dd=observer.divided_difference(samples)
    corrupted=[dict(s) for s in samples]
    corrupted[2]['Ruu']+=.001
    qvalues=[model.parent_flux,model.child_flux,model.child_flux-model.charge]
    vacua=[model.scalar.kappa*model.scalar.vacuum(q) for q in qvalues]
    delta=vacua[0]-2*vacua[1]+vacua[2]
    first,second,translation=equivalent_one_wall_example()
    alternatives=[]
    for candidate in (first,second):
        obs=ScalarWallObserver(candidate)
        qs=[candidate.parent_flux,candidate.child_flux,candidate.child_flux-candidate.charge]
        ls=[candidate.scalar.kappa*candidate.scalar.vacuum(q) for q in qs]
        D=ls[0]-2*ls[1]+ls[2]
        alternatives.append(dict(epsilon=candidate.scalar.epsilon,charge=candidate.charge,
            parent_flux=candidate.parent_flux,child_flux=candidate.child_flux,mu=candidate.scalar.mu,
            force_scale=obs.force_scale,ground_pressure=obs.ground_pressure,
            tidal_C=obs.tidal_quadratic_coefficient,ground_Lambda=ls,
            vacuum_second_difference=D,next_ground_transition_is_downward=ls[2]<ls[1],
            recovered=infer_mixing_and_charge(obs.tidal_quadratic_coefficient,D,kappa=candidate.scalar.kappa)))
    potential_gaps=[second.scalar.potential(phi+translation,qb)-first.scalar.potential(phi,qa)
        for phi in (actual_phi-.01,actual_phi,actual_phi+.01)
        for qa,qb in ((first.child_flux,second.child_flux),(first.parent_flux,second.parent_flux))]
    errors=[abs(row[name]) for row in rows for name in
        ('scalar_value_error','normal_magnitude_error','tidal_residual','direct_field_pressure_error','direct_field_mass_derivative_error')]
    errors.extend([abs(direct_work-mass_increase),abs(prediction-samples[2]['Ruu']),
                   abs(dd+observer.tidal_quadratic_coefficient)])
    errors.extend(map(abs,potential_gaps))
    if max(errors)>2e-8:raise RuntimeError('observer checks exceed declared numerical comparison tolerance')
    files=('src/recursive_horizons/nsc_atlas_geometric_observer.py',
           'src/recursive_horizons/nsc_atlas_scalar_preparation.py',
           'src/recursive_horizons/nsc_atlas_flux_conversion.py',
           'scripts/check_nsc_atlas_geometric_observer.py',
           'scripts/plot_nsc_atlas_geometric_observer.py',
           'tests/test_nsc_atlas_geometric_observer.py')
    record=dict(schema='nsc-atlas-geometric-observer-v1',
        result='One-sided local scalar reconstruction, independent tidal correlation and exact pressure-work/mass balance for the existing charged canonical scalar action.',
        geometry_inputs='Child-side R,Rdot,Rddot,A,B,Adot,Bdot in wall proper time; common normal points child to parent. Independent bulk Ruu is a test output, not an input to scalar reconstruction.',
        fixed_parameters=dict(kappa=model.scalar.kappa,omega=model.scalar.omega,epsilon=model.scalar.epsilon,
            charge=model.charge,tension=model.tension,parent_flux=model.parent_flux,bare_density=model.scalar.bare_density),
        reconstructed_reference=rows,reference_solver=meta,
        mass_balance=dict(interval_wall_proper_time=[0.,R0],inward_flux=0.,
            vacuum_subtracted_mass_change=mass_increase,independent_excitation_pressure_work=direct_work,
            gap=direct_work-mass_increase,
            conclusion='Mass above the vacuum contribution grows in this evolved coupled reference with zero wall-rest-frame scalar influx; pressure work accounts for the change.'),
        tidal_relation=dict(form='Y=Ruu-kappa*Pdot^2/(mu*e)^2=-Lambda_c-C*(P-Pc)^2',
            C=observer.tidal_quadratic_coefficient,second_divided_difference=-observer.tidal_quadratic_coefficient,
            meaning='Exact necessary local relation on any regular finite interval of a solution; independent of the unknown preparation history, not independent of the fixed action.'),
        synthetic_calibration=dict(scope='Algebraic recovery and rejection test, not observations or new bulk solutions.',
            samples=samples,calibration_indices=[0,1],withheld_index=2,calibration=calibrated,
            predicted_withheld_Ruu=prediction,withheld_error=prediction-samples[2]['Ruu'],
            computed_second_divided_difference=dd,corrupted_withheld_shift=.001,
            corrupted_prediction_error=corrupted[2]['Ruu']-prediction),
        inheritance_potential_identity=dict(fluxes=qvalues,ground_Lambda_values=vacua,second_difference=delta,
            local_C_times_spacing_over_kappa_squared=observer.tidal_quadratic_coefficient*delta/model.scalar.kappa**2,
            exact_relation='C*Delta2Lambda/kappa^2=epsilon/[2(1-epsilon)]',
            scope='Same-action potential-level relation. No composed three-generation scalar spacetime or origin of the couplings is proved; these are not the pure-flux P002 numbers.'),
        one_wall_nonidentifiability=dict(models=alternatives,scalar_translation=translation,
            exact_potential_difference='U_B(phi+3sqrt(3)/100,sqrt(7)*(q-6/5))-U_A(phi,q)=3*(q-1)*(q-3/2)',
            maximum_numerical_two_sector_potential_gap=max(map(abs,potential_gaps)),
            interpretation='The first two charge sectors have identical reduced classical scalar/gravity/wall dynamics after field translation. Same-oriented inheritance predicts different third-sector ground energies. Independent charge probes, additional matter couplings and quantum formation are outside this equivalence.',
            inverse='e^2=D/kappa+kappa/(2C); epsilon=(D/kappa)/e^2',
            limits='Three compatible physical regions are not constructed here. The second model predicts an upward next ground transition, so the ordinary descending-cap preparation cannot simply be copied.'),
        checks=dict(symbolic_zero_identities=identities(),maximum_numerical_comparison_error=max(errors),
                    comparison_tolerance=2e-8),
        limits='Requires normal geometric information, not radius or expansion alone. At Pi=0 this input set identifies only the magnitude of the normal scalar gradient; extra curvature derivatives can add information. Derivative noise and nearby calibration pressures can make inversion ill-conditioned. Local identities are necessary, not a bulk existence theorem, finite-time excited-parent evolution, empirical cosmology, or prediction of the observed vacuum scale.',
        next='Use the invariant flux, pressure-work balance and independent tidal relation to test the coupled excited-parent evolution; freeze calibrations before testing independent observations.',
        process_cpu_seconds=time.process_time()-start,
        source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_geometric_observer.py')
    if record['process_cpu_seconds']>120:raise RuntimeError('observer calculation exceeded 120 CPU seconds')
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    encoded=json.dumps(calculation(),indent=2,allow_nan=False)+'\n'
    if len(encoded.encode())>20000:raise RuntimeError('observer record exceeds 20 KB budget')
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded,end='')
