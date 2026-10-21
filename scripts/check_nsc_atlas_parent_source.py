"""Finite parent source, conformal gravity and scalar/Israel initial matching."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import sympy as sp
from recursive_horizons.nsc_atlas_parent_source import (
    ParentScalarSource, uniform_excited_parent, planar_interface_control, planar_child_energy_integral,
    planar_volume_energies)

ROOT=Path(__file__).resolve().parents[1]


def identities():
    a,ap,W,p,g0,g,U0,U,kap=sp.symbols('a ap W p g0 g U0 U kap',positive=True)
    app=(1-ap**2)/(2*a)-kap*a*(g0*g0/2+U0)/2
    F=kap*(g0*g0-g*g)*W/8+kap*(U0*W-U*W**5)/4
    pp=F-2*ap*p/a
    R=W*W*a
    beta=ap+2*a*p/W
    Rll=(app+2*ap*p/W+2*a*pp/W-2*a*p*p/W**2)/W**2
    hamilton=-4*Rll*R+2*(1-beta**2)-2*kap*(g*g/(2*W**4)+U)*R**2
    k,Q,z=sp.symbols('k Q z',real=True)
    jet=Q*(W**4-1)+2*k*z*(W**2-1)
    jump_p=k*(W-W**3)/2
    normal_jump=(Q+3*k*z+jet-2*z*jump_p/W)/W**4
    scalar_match=normal_jump-3*k*z/W**2-Q
    radius,B,sigma,P,pr=sp.symbols('radius B sigma P pr',nonzero=True,real=True)
    Bo=B-k*radius
    mean=P/sigma-(B+Bo)/radius
    Ai,Ao=mean+k/2,mean-k/2
    accel=lambda bb,aa,pressure: -(1-bb*bb)/(2*radius)-k*radius*pressure/sigma+aa*bb
    wall_acc=accel(Bo,Ao,pr-P)-accel(B,Ai,pr)
    m_jump=radius*(B*B-Bo*Bo)/2-k*radius*radius*(B+Bo)/2
    raw={'conformal_Hamiltonian':hamilton,
         'fixed_tension_spatial_jump':(-k+2*jump_p/W)/W**2+k,
         'scalar_second_normal_matching':scalar_match,
         'two_sided_wall_acceleration':wall_acc,
         'surface_mass_relation':m_jump}
    t,x,S=sp.symbols('t x S',real=True)
    child=-S*(t+x)**2/4
    parent=-S*(t*t+2*t*x-x*x)/4
    flux=(sp.diff(child,t)*sp.diff(child,x)).subs(x,0)
    energy=sp.integrate((sp.diff(child,t)**2+sp.diff(child,x)**2)/2,(x,-t,0))
    parent_energy=sp.integrate((sp.diff(parent,t)**2+sp.diff(parent,x)**2)/2+S*parent,(x,0,t))
    raw.update(tangent_child_wave=sp.diff(child,t,2)-sp.diff(child,x,2),
        tangent_parent_force=sp.diff(parent,t,2)-sp.diff(parent,x,2)+S,
        tangent_interface_value=(child-parent).subs(x,0),
        tangent_interface_normal=sp.diff(child-parent,x).subs(x,0),
        tangent_child_energy_balance=sp.diff(energy,t)-flux,
        tangent_parent_supplies_child_energy=parent_energy+energy,
        tangent_pressure_inflow_relation=-S*child.subs(x,0)-flux)
    result={name:str(sp.factor(value)) for name,value in raw.items()}
    if any(v!='0' for v in result.values()):raise RuntimeError(result)
    return result


def calculation():
    start=time.process_time()
    model=ParentScalarSource()
    background,bg_meta=model.prepared.solve(tolerance=1e-9)
    rows=[]
    primary=None
    for amplitude in (0.,.01,.02):
        src=ParentScalarSource(amplitude=amplitude)
        solved,meta=src.solve(background,tolerance=1e-9)
        row=dict(amplitude=amplitude,solver=meta,wall=src.wall(background,solved),
                 constraints=src.diagnostics(background,solved),energy=src.energy(background,solved),
                 mass_checks=src.mass_checks(background,solved))
        rows.append(row)
        if amplitude==.02:primary=solved
        if max(row['constraints']['weighted_Hamilton_residual'])>1e-6:
            raise RuntimeError('physical Hamiltonian residual exceeds declared tolerance')
        if max(abs(c['gap']) for c in row['mass_checks'])>1e-7:
            raise RuntimeError('integrated mass constraint failed')
        if max(abs(v) for k,v in row['wall'].items() if k.endswith('_gap'))>1e-9:
            raise RuntimeError('initial wall matching failed')
    coarse_bg,_=model.prepared.solve(tolerance=1e-7)
    coarse,_=model.solve(coarse_bg,tolerance=1e-7)
    wc,wf=model.wall(coarse_bg,coarse),rows[-1]['wall']
    ec=model.energy(coarse_bg,coarse)
    spreads={key:abs(wc[key]-wf[key]) for key in ('radius','second_jet','wall_factor')}
    for key in ('wall_areal_acceleration','wall_field_acceleration','outward_flux_slope'):
        spreads[key]=abs(wc['inside'][key]-wf['inside'][key])
    spreads['total_scalar_excess']=abs(ec['total_scalar_excess']-rows[-1]['energy']['total_scalar_excess'])
    if max(spreads.values())>2e-7:raise RuntimeError('coarse/fine source comparison failed')
    frozen,_=model.solve(background,match_scalar_jet=False)
    bad_wall=model.wall(background,frozen)
    frozen_diag=model.diagnostics(background,frozen)
    if abs(bad_wall['scalar_second_jet_gap'])<.3:
        raise RuntimeError('frozen-layer obstruction control did not reproduce')
    uniform=uniform_excited_parent(model.prepared)
    planar=planar_interface_control(.01,force_jump=uniform['field_force_jump'],omega=model.prepared.scalar.omega)
    volume=planar_child_energy_integral(.01,force_jump=uniform['field_force_jump'],omega=model.prepared.scalar.omega)
    volume_fine=planar_child_energy_integral(.01,force_jump=uniform['field_force_jump'],omega=model.prepared.scalar.omega,order=36)
    both_volumes=planar_volume_energies(.01,force_jump=uniform['field_force_jump'],omega=model.prepared.scalar.omega)
    if abs(volume-planar['child_energy_per_area'])>1e-13 or abs(volume-volume_fine)>1e-13:
        raise RuntimeError('retarded volume energy does not match integrated boundary flux')
    if abs(both_volumes['net'])>1e-13:raise RuntimeError('parent loss and child gain do not balance')
    files=('src/recursive_horizons/nsc_atlas_parent_source.py',
           'src/recursive_horizons/nsc_atlas_scalar_preparation.py',
           'src/recursive_horizons/nsc_atlas_flux_conversion.py',
           'scripts/check_nsc_atlas_parent_source.py','tests/test_nsc_atlas_parent_source.py',
           'scripts/plot_nsc_atlas_parent_source.py')
    record=dict(schema='nsc-atlas-parent-source-v1',
        question='Which finite parent preparations satisfy the coupled initial constraints, and what early energy transfer follows under their declared interface regularity?',
        action='Unchanged Einstein/canonical scalar/four-form/constant-tension charged membrane; no extra species or source coupling.',
        inputs=dict(epsilon=model.prepared.scalar.epsilon,omega=model.prepared.scalar.omega,
                    kappa=model.prepared.scalar.kappa,parent_flux=model.prepared.parent_flux,
                    charge=model.prepared.charge,tension=model.prepared.tension,bare_density=0.,
                    packet_center=model.center,packet_half_width=model.half_width,collar_width=model.collar_width,
                    packet='f=A*(1-z^2)^6 for abs(z)<1, z=(parent seed distance from wall-center)/half_width; zero outside.',
                    collar='g=lambda*d^2*(1-d/collar_width)^6/2 for 0<=d<collar_width; lambda is solved, not fitted.',
                    Cauchy_state='K_ij=0, scalar Pi=0; original scalar values and first seed-coordinate derivatives at the wall; compact parent displacement. Physical normal derivatives scale with Psi^-2 and remain continuous.'),
        units='Existing normalized variables, c=1. Proper energy is an integral on the sourced slice; no SI ruler or global ADM energy is selected.',
        solved=rows,background_solver=bg_meta,
        exact_wall_rule='lambda=mu*e*(Psi_W^4-1)+2*k*phi0_n*(Psi_W^2-1).',
        rejection_control=dict(name='Freeze the old scalar second-normal layer after changing the geometry',
            weighted_Hamilton_residual=frozen_diag['weighted_Hamilton_residual'],
            scalar_compatibility_gap=bad_wall['scalar_second_jet_gap'],
            wall_field_acceleration_gap=bad_wall['field_acceleration_gap'],
            interpretation='Passing the Einstein constraint does not give C2 scalar compatibility at the initial corner. A causal finite-energy weak transient is not ruled out.'),
        uniform_parent_preparation=uniform,
        conditional_weak_transfer=dict(
            assumptions='A causal local development with continuous field/normal matching, bounded one-sided geometric coefficients and no extra incoming initial waves; existence of the full Einstein/membrane development is not established here.',
            scalar_wall_displacement='-S*tau^2/4+O(tau^3), S=mu*e, from the retarded tangent equation.',
            inward_flux='S^2*tau^2/4+O(tau^3)',
            child_excitation_per_initial_area='S^2*tau^3/12+O(tau^4)',
            pressure_relation='(P(0)-sigma*K_mean)/F_in tends to 1; same charge/scalar source fixes both coefficients.',
            contribution='Exact round constrained initial data and a nonzero conditional causal transfer coefficient. This allows an initial weak corner instead of treating C2 compatibility as a compulsory physical gate.'),
        planar_massive_control=dict(**planar,independent_volume_energy=volume,
            volume_boundary_energy_gap=volume-planar['child_energy_per_area'],
            quadrature_order_gap=volume_fine-volume,two_sided_volume_energy=both_volumes),
        checks=dict(symbolic_zero_identities=identities(),two_tolerance_spreads=spreads,
                    numerical_comparison_tolerance=2e-7),
        initial_flux=dict(preparation='Smooth compact packet with amplitude .02, distinct from the uniform weak-corner seed.',
                          sign='Positive is child-to-parent; F=-T(n,u).',
                          value=0.,slope=wf['inside']['outward_flux_slope'],
                          conclusion='The compatible preparation has an initially outward flux slope, not demonstrated matter deposition in the child.'),
        validity='Three finite Cauchy examples with second-order smooth interface compatibility; also one exact uniform-parent constrained seed with a weak initial corner and conditional leading transfer coefficients. The planar massive control has its own fixed-background domain.',
        limits='No full coupled finite-time evolving solution, all-order corner compatibility, global existence, stability, finite-time curved deposited energy/error bound, isotropization, clustering, selected source history or observational abundance. Initial geometry differences are not instantaneous causal signals.',
        next='Construct the coupled causal evolution of the uniform excited parent and test the nonzero limiting flux and pressure relation; retain the smooth constrained source as a distinct initial-data control.',
        process_cpu_seconds=time.process_time()-start,
        source_sha256={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce='python scripts/run_python.py scripts/check_nsc_atlas_parent_source.py')
    if record['process_cpu_seconds']>120:raise RuntimeError('source calculation exceeded 120 CPU seconds')
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    encoded=json.dumps(calculation(),indent=2,allow_nan=False)+'\n'
    if len(encoded.encode())>24000:raise RuntimeError('compact source record exceeded 24 KB')
    if args.output:args.output.write_text(encoded,encoding='utf-8',newline='\n')
    else:print(encoded,end='')
