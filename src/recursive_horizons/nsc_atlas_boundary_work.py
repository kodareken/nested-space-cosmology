"""Average geometric work required by the finite fluid/gravity continuation.

This closes the fourth-order averaged boundary Ward identity under a stated
homogeneous-strain ansatz. It does not solve the local fourth-order bulk.
"""
import math


def work_coefficient(initial_energy,energy,loss,*,length=4.,nu=.001):
    V=length**3
    return dict(b2=2*(initial_energy-energy)/(3*V*nu),
        db2_dslowtime=2*loss/(3*V*nu),
        proper_H_over_epsilon4=2*loss/(3*V*nu**1.5),
        flat_boundary_Ward_residual_over_epsilon4=-2*loss/(V*nu**2),
        work_term_over_epsilon4=2*loss/(V*nu**2),
        maximum_energy_budget_b2=2*initial_energy/(3*V*nu))


def calculation(flow):
    rows=[]
    for case in flow['cases']:
        E0=case['initial']['energy']
        values=work_coefficient(E0,case['final']['energy'],case['final']['dissipation'])
        values.update(n=case['n'],dt=case['dt'],endpoint=case['endpoint'],
            b2_from_integrated_loss=2*case['integrated_viscous_loss']/(3*64*.001),
            history=[dict(t=row['t'],**work_coefficient(E0,row['energy'],row['dissipation'])) for row in case['history']])
        rows.append(values)
    space=abs(rows[2]['b2']-rows[1]['b2']);time=abs(rows[3]['b2']-rows[2]['b2'])
    return dict(schema='nsc-atlas-boundary-work-v1',
        assumptions=['same finite periodic fluid state and zero force',
                     'vacuum bulk flux through the cutoff, with Brown-York convention 16piG5=1',
                     'background boundary pressure p0=1/sqrt(nu), rho0=0',
                     'homogeneous boundary strain b=1+epsilon^2 b2(t), b2(0)=0',
                     'no additional order-epsilon^2 averaged energy mode or imposed source'],
        identities=dict(
            exact_boundary_law='d rho_bar/d t_hat + 3 H_b*(rho_bar+p_bar)=0 at the retained averaged order',
            clocks='t=epsilon^2 tau; t_hat=sqrt(nu) tau',
            energy='rho_bar=2 epsilon^2 E/(V nu^(3/2)); E_prime=-D',
            residual='Flat cutoff residual/epsilon^4=-2D/(Vnu^2)',
            work='H_b/epsilon^4=b2_prime/sqrt(nu); 3H_b p0/epsilon^4=3b2_prime/nu',
            selected_rate='b2_prime=2D/(3Vnu)',
            integrated='b2(t)=2(E(0)-E(t))/(3Vnu)',
            orientation='D>=0 makes b2 nondecreasing; with zero mean periodic flow D=0 only for u=0.',
            Poincare='D>=2nu(2pi/L)^2 E; finite smooth interval E(t)<=E0 exp[-2nu(2pi/L)^2 t]'),
        cases=rows,
        refinement=dict(last_spatial_change=space,half_step_change=time,
            empirical_indicator=2*max(space,time),
            not_a_certified_bound=True),
        lower_growth_coefficient=rows[-1]['maximum_energy_budget_b2']*(1-math.exp(-2*.001*(2*math.pi/4)**2*.05)),
        scope=dict(average_Ward_identity_closed=True,local_bulk_Einstein_order4_solved=False,
            horizon_regular_extension_solved=False,new_dimension_created=False,
            accelerated_cosmology_predicted=False,SI_scale_selected=False,
            time_origin_derived=False,
            interpretation='A conditional energy-versus-boundary-scale law and a dissipative ordering of states. It does not derive the existence of time or cosmological dark energy.'),
        alternatives='A source flux, nonhomogeneous work, or an additional allowed energy sector would alter this relation; the Ward identity alone does not select the strain ansatz.')
