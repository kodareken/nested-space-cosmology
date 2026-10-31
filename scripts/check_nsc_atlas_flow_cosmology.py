"""Frozen P006 flow-to-vacuum comparison and same-state receiving pressure.

No null grid, new geometry or observation score is computed.  V_star is
retained only as the historical background/normalization of this example.
This script supplies no inferred Omega_Lambda(today).
"""
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp

from recursive_horizons.nsc_atlas_horizon_response import solve_source
from recursive_horizons.nsc_atlas_characteristic_handoff import (
    C, Diamond, fields_and_derivatives, geometry, inherited_reference,
)


def receiver_pressure(rtol, atol):
    record=json.loads(Path('results/summaries/nsc-atlas-horizon-response-v1.json').read_text(encoding='utf-8'))
    axis=np.array(record['angular_transfer']['axis_in_n_coordinates'])
    source=solve_source(endpoint_u=2.,rtol=rtol,atol=atol)
    domain=Diamond.inherited(source.parent)
    u=np.array(domain.u0+domain.width);v=np.array(domain.v0+domain.width)
    psi,du,dv=inherited_reference(source,u,v)
    F,Fu,Fv=fields_and_derivatives(source,u,v,psi,du,dv)
    eta,chi,a,_=geometry(source,u,v);a=float(a);chi=float(chi)
    field=float(axis[:2]@F)
    Pi=float(axis[:2]@(Fu+Fv))/a
    radial=float(axis[:2]@(Fv-Fu))/a
    angular=(float(F@F)-field**2)/(a*a*math.sin(chi)**2)
    grad2=radial**2+angular
    mass=source.parent.mass_squared*field**2
    rho=(Pi**2+grad2+mass)/2
    pressure=(Pi**2-grad2/3-mass)/2

    # Independent ambient-S3 gradient, rather than the null-coordinate
    # conversion used above.  X4=cos(chi), Xa=sin(chi)*na for a=1,2,3.
    dimensionless_time=math.asinh(math.tan(float(eta)))
    state=source.state(dimensionless_time)
    X=np.r_[math.sin(chi)*axis,math.cos(chi)]
    ambient=np.array([C[0]*state[2],C[1]*state[4]*X[3],0.,C[1]*state[4]*X[1]])
    grad2_check=(ambient@ambient-(X@ambient)**2)/a**2
    Pi_check=C[0]*state[3]*X[0]+C[1]*state[5]*X[1]*X[3]
    return dict(rtol=rtol,atol=atol,
        proper_event=dimensionless_time/source.parent.H,Ht=dimensionless_time,
        field_over_epsilon=field,Pi_over_epsilon=Pi,
        gradient_squared_over_epsilon_squared=grad2,
        rho_over_epsilon_squared=rho,pressure_over_epsilon_squared=pressure,
        rho_over_epsilon_squared_Vstar=rho/source.parent.vacuum,
        pressure_over_epsilon_squared_Vstar=pressure/source.parent.vacuum,
        w_trace=pressure/rho,
        direct_embedding_gradient_gap=float(grad2_check-grad2),
        direct_embedding_Pi_gap=float(Pi_check-Pi))


def late_pressure_identity():
    A,B,beta,nu,mu=sp.symbols('A B beta nu mu',real=True)
    # phi=a^-3/2(A cos(nu*t)+B sin(nu*t)), beta=3H/2.
    # The displayed averages are coefficients before the common a^-3.
    mean_kinetic=((nu*B-beta*A)**2+(nu*A+beta*B)**2)/2
    mean_mass=mu**2*(A*A+B*B)/2
    pressure=sp.expand((mean_kinetic-mean_mass)/2)
    residual=sp.factor(pressure.subs(nu**2,mu**2-beta**2))
    assert residual==0
    return dict(beta='3H/2',frequency_relation='nu_phys^2+beta^2=mu^2',
        mean_kinetic_coefficient=str(sp.factor(mean_kinetic)),
        mean_mass_coefficient=str(mean_mass),
        phase_averaged_leading_pressure_residual=str(residual),
        density_scale='a^-3',gradient_stress_scale='a^-5',
        phase_averaged_w_limit=0,
        scope='Late massive scalar on the supplied asymptotically de Sitter background. Instantaneous pressure oscillates; this is not a generated vacuum floor.')


def calculation():
    import hashlib
    import subprocess
    from recursive_horizons.nsc_atlas_offset_generator import offset_value
    forecast_path=Path('results/summaries/nsc-atlas-flow-omega-forecast-v1.json')
    forecast=json.loads(forecast_path.read_bytes())
    assert not forecast['comparison_evaluated']
    rows=[receiver_pressure(rtol,atol) for rtol,atol in
          ((2e-12,2e-14),(3e-13,3e-15))]
    assert max(abs(r['direct_embedding_gradient_gap']) for r in rows)<1e-16
    assert max(abs(r['direct_embedding_Pi_gap']) for r in rows)<1e-12
    assert abs(rows[0]['w_trace']-rows[1]['w_trace'])<1e-8
    # Vacuum shifts are invisible to a null stress projection.  This is
    # an independent algebraic identity, not an alternate fitted model.
    c=sp.symbols('c',real=True)
    metric=sp.diag(-1,1,1,1); null=sp.Matrix([1,1,0,0])
    null_vacuum_shift=sp.expand((null.T*(-c*metric)*null)[0])
    assert null_vacuum_shift==0
    late=late_pressure_identity()
    return dict(schema='nsc-flow-omega-evaluation-v1',prediction_id='NSC-P006',
        forecast_sha256=hashlib.sha256(forecast_path.read_bytes()).hexdigest(),
        protocol_commit='057ecd9e',analysis_parent_commit=subprocess.check_output(
            ['git','rev-parse','HEAD'],text=True).strip(),
        Omega_Lambda_today=None,comparison_pass=False,empirical_rejection=False,
        status='prediction_not_identified',
        missing_closure='No constitutive receiving-vacuum/expansion history selecting rho_vac_minus, H_today and event_today. The offset map fixes phases, not that history.',
        no_supplied_vacuum_as_answer=True,
        Lambda_eff_Lref_squared=None,
        added_vacuum_Lambda_Lref_squared_late_leading_field=0,
        late_limit_is_not_today=True,
        null_projection_vacuum_shift=str(null_vacuum_shift),
        receiving_pressure=rows,tighter_w_gap=rows[1]['w_trace']-rows[0]['w_trace'],
        receiving_w_minus_vacuum_value=rows[1]['w_trace']+1,
        late_field=late,
        phase_averaged_late_w_minus_vacuum_value=1,
        child=dict(initial_strict_future_trapping=False,amplitude_interval=[0,1],
            future_child_formation=None,threshold_amplitude=None,child_Lambda_Lref_squared=None,
            child_grandchild_fixed_point=None,
            criterion='theta_plus*theta_minus=-H_surface^2 <= 0 at Kij=0; a future child requires a later trapped surface and an admissible interior continuation.',
            epsilon_one_nonlinear_status='under-resolved, not used for claims'),
        parameters=dict(action_alpha=float(offset_value('action')),
            alpha_status='Derived mode ratio within the supplied action/frequency normalization; identification as the generator offset is assumed.',
            n100_status='Chosen reporting count, not a selected epoch.',
            omega0_20_status='Supplied reference scalar frequency.',
            spatial_dimension_status='3D completion is proved for the specified kernel and two accumulators; selection of those construction rules is not derived.'),
        charge=dict(signed_charge_resolved_flux=0,baryon_photon_eta=None,
            requested_eta_target=6.1e-10,
            reason='The real two-scalar horizon sector has no baryon current, charge-biased term or photon density; crossing is charge blind. A symplectic Green current is not particle charge.'),
        no_new_evolution_campaign=True,
        evidence_sources=['results/summaries/nsc-atlas-parent-circulation-v1.json',
            'results/summaries/nsc-atlas-horizon-response-v1.json',
            'results/summaries/nsc-atlas-characteristic-handoff-v1.json',
            'results/summaries/nsc-atlas-finite-parent-data-v1.json'])


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=calculation()
    raw=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(raw,encoding='utf-8',newline='\n')
        print('P006: prediction not identified; conditional receiving w =',result['receiving_pressure'][1]['w_trace'])
    else:print(raw)


if __name__=='__main__':main()
