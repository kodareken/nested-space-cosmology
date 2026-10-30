"""Specified connections from circulation to boundary geometry and horizon flux.

The five-dimensional Rindler lift and four-dimensional rotating-horizon
scalar are distinct attributed sectors. Their equality as an NSC history is
not assumed. This module calculates their common stress/flux interfaces.
"""
import math
from .nsc_atlas_vortex_transfer import reynolds,mean_theta_rate,energy,dissipation


def boundary_coefficients(E,D,nu=.001):
    if E<0 or D<0 or nu<=0:raise ValueError('nonnegative energy/loss and positive viscosity required')
    return dict(integrated_boundary_energy_coefficient=2*E/nu**1.5,
        boundary_energy_slow_time_derivative=-2*D/nu**1.5,
        integrated_Kretschmann_coefficient=-3*D/nu**3,
        scope='Slow-coordinate integrals: boundary energy at order epsilon^2, Lorentzian Riemann-squared at leading order epsilon^4. Not a cosmological Lambda.')


def horizon_flux(omega,winding,angular_velocity,amplitude_norm=1.):
    """Scalar energy and angular momentum flux INTO a Killing horizon."""
    if omega<=0 or amplitude_norm<0:raise ValueError('positive frequency and nonnegative horizon norm required')
    delta=omega-winding*angular_velocity
    E=amplitude_norm*omega*delta;J=amplitude_norm*winding*delta
    return dict(offset=delta,energy=E,angular_momentum=J,
        corotating_energy=E-angular_velocity*J,
        null_stress_integral=amplitude_norm*delta*delta)


def connections():
    nu=.001;rows=[]
    for radius in (.75,1.,1.25):
        R=reynolds(radius,0.,4)
        rows.append(dict(radius=radius,covariance=R.tolist(),
            proper_boundary_energy_coefficient=float(R.trace()/nu**1.5),
            offdiagonal_boundary_stress_coefficient=float(R[0,1]/nu**1.5),
            boundary_momentum_slow_time_rate=float(mean_theta_rate(radius,4)/nu)))
    return dict(
        hydrodynamic_lift=dict(source='Bredberg et al 1101.2451v2, equations 14-25',
            fluid_spatial_dimensions=3,bulk_spacetime_dimensions=5,cutoff=nu,
            scaling='v=epsilon*u(epsilon*x,epsilon^2*tau); P=epsilon^2*p',
            proper_boundary_clock='t_hat=sqrt(nu)*tau',
            offset_tensor='R_ij=<u_i u_j>-<u_i><u_j>',
            boundary_excess='rho_BY=epsilon^2*tr(R)/nu^(3/2) for zero mean; momentum density=epsilon*U_i/nu',
            metric_covariance='delta g_ij=(1-r_bulk/nu)*epsilon^2*R_ij/nu at fixed mean pressure',
            vacuum_error='Einstein tensor vanishes through epsilon^3 for the published metric; higher orders are not bounded by this calculation.',
            Kretschmann='K=epsilon^4*(-3|curl_y u|^2/nu^2)+higher orders at finite bulk radius',
            initial=boundary_coefficients(energy(4),dissipation(4)),local_rows=rows,
            averaged_metric_is_a_closed_vacuum_solution=False,
            inherited_4D_parent_matched=False,
            interpretation='A computed gravitational boundary/curvature response to the same finite fluid state. Not an attachment of the 5D vacuum patch to the 4D positive-vacuum cap.'),
        rotating_horizon=dict(source='Canonical scalar stress and Killing generator; superradiant/synchronization comparison to 1403.2757 and 1501.06570',
            generator='k=xi+Omega_H eta; psi=A exp[i(m phi-omega t)]',
            offset='Delta_H=omega-m Omega_H',
            norm='N=2 integral_horizon |A|^2 dA',
            law='F_E=N omega Delta_H; F_J=N m Delta_H; F_E-Omega_H F_J=N Delta_H^2 >=0',
            normalized_controls=[dict(omega_over_Omega=x,**horizon_flux(x,1,1)) for x in (.5,1.,1.5)],
            positive_frequency_antiparticle_internal_charge_flip='Same omega,m and opposite internal charge have identical stress flux on this neutral rotating background.',
            spectrum_selected=False,scattering_amplification_computed=False,
            child_spacetime_created=False,
            interpretation='Mode-dependent transfer or rotational-energy extraction; not material escaping from inside the horizon and not a new-dimension theorem.'),
        spherical_balance=dict(energy='E_surface=gamma A at fixed enclosed volume, gamma>0 isotropic',
            theorem='A^3>=36*pi*V^2, equality only for a sphere (up to rigid motions)',
            pressure='Delta p=2 gamma/R at stationary spherical balance',
            rotation='For angular momentum J and moment of inertia I, E_rot>=J^2/(2I); compression at fixed J raises this lower bound when I decreases.',
            scope='Isoperimetric minimizer and angular-momentum energy bound, not a count-of-forces rule or proof all compact geometries are spherical.'),
        dimensional_handoff=dict(reduced_coordinate='Areal radius or a generating curve is a reduced description; suppressed angular coordinates remain part of the geometry.',
            required='A regular metric and field transition with conserved flux and declared topology; an actual dimension change needs a different rank/geometry construction.',
            established=False))
