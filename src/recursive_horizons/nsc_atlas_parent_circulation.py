"""Compensated circulation and a differential light record on a smooth S3 parent.

Uses the retained VolumeParent vacuum and canonical real scalar. The coherent
two-harmonic preparation and anchored light guides are declared benchmark data.
Only the leading quadratic vector response is calculated; it is not a complete
finite-amplitude parent/daughter solution or a particle-number asymmetry.
"""
from dataclasses import dataclass
import math

import mpmath as mp
import numpy as np
from scipy.integrate import quad, solve_ivp

from .nsc_atlas_offset_generator import offset_value
from .nsc_atlas_parent_scale_selector import VolumeParent

VOLUME = 2 * math.pi**2
HARMONIC_PRODUCT = 2 * math.sqrt(6) / math.pi**2


@dataclass(frozen=True)
class Parent:
    vacuum: float
    mass_squared: float
    offset: float
    kappa: float = 1.

    @classmethod
    def retained(cls):
        model = VolumeParent()
        minimum = next(row for row in model.vacuum_roots() if row['kind'] == 'minimum')
        with mp.workdps(60):
            alpha = float(offset_value('action'))
        return cls(minimum['vacuum'], 20.**2 * math.exp(-model.d * minimum['s']), alpha)

    @property
    def H(self):
        return math.sqrt(self.kappa * self.vacuum / 3)

    @property
    def reference_frequency(self):
        return math.sqrt(self.mass_squared + 3 * self.H**2)

    def initial(self, eigenvalues=(3., 8.)):
        # <rho_chi>(0) = epsilon^2 * vacuum; equal mode energy, zero velocity.
        q = [math.sqrt(self.vacuum * VOLUME / (self.mass_squared + l * self.H**2))
             for l in eigenvalues]
        # a, a_dot, q1, p1, q2, p2, J, integral beta dt, proper readout, work
        return np.array([1 / self.H, 0., q[0], 0., q[1], 0., 0., 0., 0., 0.])


def proper_rhs(z, parent, eigenvalues=(3., 8.)):
    a, v, q1, p1, q2, p2, J = z[:7]
    h = v / a
    l1, l2 = eigenvalues
    mu2 = parent.mass_squared
    work = (-3 * h * a**3 * (p1*p1 + p2*p2)
            + 3 * h * a**3 * mu2 * (q1*q1 + q2*q2)
            + h * a * (l1*q1*q1 + l2*q2*q2)) / 2
    return np.array([v, parent.H**2 * a, p1, -3*h*p1-(mu2+l1/a**2)*q1,
                     p2, -3*h*p2-(mu2+l2/a**2)*q2,
                     -3*h*J+(l2-l1)*q1*q2/a**2, parent.kappa*J/5, 1., work])


def phase_speed(z, parent, first_eigenvalue=3.):
    a, v, q, p = z[:4]
    w = parent.reference_frequency
    return w * (p*p + 3*v/a*q*p + (parent.mass_squared+first_eigenvalue/a**2)*q*q) / (w*w*q*q+p*p)


def energy(z, parent, eigenvalues=(3., 8.)):
    a, _, q1, p1, q2, p2 = z[:6]
    return a**3 * (p1*p1+p2*p2 + parent.mass_squared*(q1*q1+q2*q2)
                  + (eigenvalues[0]*q1*q1+eigenvalues[1]*q2*q2)/a**2) / 2


def horizon_flux(t, z, parent):
    """Integrated T_kk and Killing-energy flux, divided by epsilon^2.

    North-pole observer, k=partial_t-a^-1 partial_chi, normalized by n.k=-1.
    The static Killing vector on its future horizon is tanh(Ht)*k.
    """
    u = math.tanh(parent.H*t)
    v = 1 / math.cosh(parent.H*t)
    q1,p1,q2,p2=z[2:6]
    A1 = math.sqrt(2)/math.pi*v*(p1-parent.H*u*q1)
    A2 = math.sqrt(12)/math.pi*v*(u*p2+parent.H*(1-2*u*u)*q2)
    null_flux = 4*math.pi/(3*parent.H**2)*(A1*A1+A2*A2)
    return null_flux, u*null_flux


def static_patch_energy(t, z, parent, order=48):
    """Independent volume integral of T(n,K), divided by epsilon^2.

    The null boundary is chi_H=pi/2-atan(sinh(Ht)). Integration is over
    the smooth parent field; an outside continuation already exists here.
    """
    eta=math.atan(math.sinh(parent.H*t))
    x,w=np.polynomial.legendre.leggauss(order)
    chi=(x+1)*(math.pi/2-eta)/2
    weight=w*(math.pi/2-eta)/2
    a=z[0];q1,p1,q2,p2=z[2:6]
    sin,cos=np.sin(chi),np.cos(chi)
    C1,C2=math.sqrt(2)/math.pi,math.sqrt(12)/math.pi
    f=np.array([C1*q1*sin,C2*q2*sin*cos])
    fchi=np.array([C1*q1*cos,C2*q2*(cos*cos-sin*sin)])
    ft=np.array([C1*p1*sin,C2*p2*sin*cos])
    rho=2*math.pi/3*np.sum(ft*ft+(fchi*fchi+2*f*f/(sin*sin))/a**2+parent.mass_squared*f*f,axis=0)
    momentum=4*math.pi/3*np.sum(ft*fchi,axis=0)
    integrand=a**3*sin*sin*(cos*rho-parent.H*math.tanh(parent.H*t)*sin*momentum)
    return float(np.dot(weight,integrand))


def evolve(*, parent=None, clock='phase', passes=300, endpoint=None,
           eigenvalues=(3., 8.), method='DOP853', rtol=3e-12, atol=2e-14):
    """Finite linear field evolution, with an independent phase/proper control.

    For the phase parametrization an event advances Theta by 2*pi*alpha.
    This synchronization with the drawing is a candidate rule, not a derived
    origin of time or a selected cosmic nesting depth.
    """
    parent = parent or Parent.retained()
    if clock not in ('phase', 'proper') or passes <= 0 or passes > 1000:
        raise ValueError('phase/proper clock and 0 < passes <= 1000 required')
    stop = 2*math.pi*parent.offset*passes if clock == 'phase' else endpoint
    if stop is None or stop <= 0:
        raise ValueError('positive finite endpoint required')
    def rhs(_, z):
        result = proper_rhs(z, parent, eigenvalues)
        if clock == 'phase':
            speed = phase_speed(z, parent, eigenvalues[0])
            if speed <= 0:
                raise ValueError('phase stopped being monotone')
            return result / speed
        return result
    sol = solve_ivp(rhs, (0., stop), parent.initial(eigenvalues), method=method,
                    rtol=rtol, atol=atol, dense_output=True, max_step=.12 if clock == 'phase' else .03)
    if not sol.success or sol.nfev > 160000:
        raise RuntimeError('finite field evolution failed or budget exceeded')
    return sol


def initial_constraint(parent=None):
    """Exact coefficients of Psi=1+epsilon^2*u+O(epsilon^4), not a nonlinear solve."""
    p = parent or Parent.retained()
    z = p.initial(); q1, q2 = z[2], z[4]
    A = p.mass_squared / p.H**2
    B1 = q1*q1*2/math.pi**2
    B2 = q2*q2*12/math.pi**2
    B12 = 2*q1*q2*HARMONIC_PRODUCT
    c1 = p.kappa * B1 * (A-1) / 40
    c3 = p.kappa * B2 * (A-4) / 168
    c2 = p.kappa * B2 / 40 + 2*c3/5
    c4 = p.kappa * B12 * (A-2) / 96
    c0 = (-p.kappa*B1/8-2*c1-4*c2) / 3
    bound = abs(c0)+abs(c1)+abs(c2)+abs(c3)/4+abs(c4)/(3*math.sqrt(3))
    return dict(basis=['1', 'X1^2', 'X2^2+X4^2', 'X2^2*X4^2', 'X1*X2*X4'],
                coefficients=[c0,c1,c2,c3,c4],
                mean_u=c0+c1/4+c2/2+c3/24,
                absolute_u_upper_bound=bound,
                scope='Quadratic corrector. Symmetry-restricted IFT gives an exact local data branch, without an amplitude bound.')


def ring_return(parent, radius):
    if not 0 < radius < .25:
        raise ValueError('guided loop must return before the future conformal endpoint')
    return math.asinh(math.tan(2*math.pi*radius)) / parent.H


def differential_coefficient(parent, proper_solution, radius=.125, separation=math.pi/6):
    t = ring_return(parent, radius)
    z = proper_solution.sol(t)
    B = math.sqrt(1-radius*radius)
    coeff = 2*z[0]*HARMONIC_PRODUCT*B*radius*(1-math.cos(separation))*z[7]
    causal_margin = math.pi/2-(2*math.pi*radius+separation/2+math.asin(radius))
    return dict(radius=radius, separation=separation, return_proper=t,
                return_phase_readout=None, integral_beta=z[7],
                delay_over_epsilon_squared=coeff, H_delay_over_epsilon_squared=parent.H*coeff,
                fractional_delay_over_initial_energy_fraction=coeff/t,
                common_comparison_conformal_margin=causal_margin)


def direct_light_control(parent, solution, *, radius=.125, epsilon=.05,
                         separation=math.pi/6, azimuths=8, frame_amplitude=0.):
    """Null integration of the vector-only metric, a control of the quadratic term.

    Other second-order metric sectors cancel in the phase-averaged orientation
    difference analytically. This truncated metric is not a nonlinear Einstein
    solution. No error bar for omitted fourth-order physics is claimed.
    """
    if epsilon <= 0 or azimuths < 4:
        raise ValueError('positive amplitude and >=4 azimuths required')
    phases = 2*math.pi*np.arange(azimuths)/azimuths
    b = np.repeat([math.sqrt(1-radius**2), math.sqrt(1-radius**2)*math.cos(separation)], 2*azimuths)
    sign = np.tile(np.repeat([1.,-1.],azimuths),2)
    phi0 = np.tile(phases,4)
    def rhs(xi, times):
        z = solution.sol(times)
        a, J = z[0], z[6]
        Wphi = HARMONIC_PRODUCT*b*radius**2*(1-np.cos(2*(phi0+sign*xi))/3)
        shift = epsilon**2*a*a*(parent.kappa*J/5*Wphi
                 + frame_amplitude*parent.H*np.cos(parent.H*times)*radius**2)
        return sign*shift + np.sqrt(shift*shift+a*a*radius**2)
    sol = solve_ivp(rhs, (0.,2*math.pi), np.zeros(4*azimuths),method='DOP853',
                    rtol=3e-14, atol=3e-15, max_step=.015)
    if not sol.success or np.max(sol.y[:,-1]) > solution.t[-1]:
        raise RuntimeError('light control left its finite field interval')
    arrivals = sol.y[:,-1].reshape(2,2,azimuths)
    ring_delays = np.mean(arrivals[:,0]-arrivals[:,1], axis=1)
    result = float(ring_delays[0]-ring_delays[1])/epsilon**2
    expected = differential_coefficient(parent,solution,radius,separation)['delay_over_epsilon_squared']
    return dict(radius=radius, epsilon=epsilon, azimuths=azimuths, frame_amplitude=frame_amplitude,
                delay_over_epsilon_squared=result, relative_difference_from_quadratic=result/expected-1,
                rhs_evaluations=sol.nfev)


def calculation():
    p = Parent.retained()
    phase = evolve(parent=p)
    tend = float(phase.y[8,-1])
    proper = evolve(parent=p,clock='proper',endpoint=tend)
    radau = evolve(parent=p,clock='proper',endpoint=tend,method='Radau',rtol=2e-10,atol=2e-13)
    samples = np.linspace(0.,phase.t[-1],601)
    yphase = phase.sol(samples)
    yp = proper.sol(yphase[8])
    zr = proper.sol(np.linspace(0.,tend,1001))
    q1,p1,q2,p2 = zr[2:6]
    Jdirect = p1*q2-q1*p2
    scaleJ = max(abs(zr[6]))
    E0 = energy(proper.y[:,0],p)
    controls = dict(
        phase_proper_max_state_difference=float(np.max(np.abs(yphase-yp))),
        independent_integrator_max_state_difference=float(np.max(np.abs(proper.sol(np.linspace(0.,tend,301))-radau.sol(np.linspace(0.,tend,301))))),
        wronskian_relative_residual=float(np.max(abs(Jdirect-zr[6]))/scaleJ),
        closed_parent_constraint_max_residual=float(np.max(abs(zr[1]**2+1-p.H**2*zr[0]**2))),
        energy_work_relative_residual=float(abs(energy(proper.y[:,-1],p)-E0-proper.y[9,-1])/E0),
        min_phase_speed=float(min(phase_speed(z,p) for z in zr.T)),
        phase_monotonic_sufficient_margin=p.mass_squared-9*p.H**2/4,
        phase_rhs_evaluations=phase.nfev, proper_rhs_evaluations=proper.nfev,radau_rhs_evaluations=radau.nfev)
    rows = [differential_coefficient(p,proper,r) for r in (.125,.0625)]
    for row in rows:
        # Invert the monotone readout rather than give t the role of generator.
        from scipy.optimize import brentq
        theta = brentq(lambda v:phase.sol(v)[8]-row['return_proper'],0.,phase.t[-1])
        row['return_phase_readout']=theta
        row['offset_passes_at_return']=theta/(2*math.pi*p.offset)
        radau_coefficient=differential_coefficient(p,radau,row['radius'])['delay_over_epsilon_squared']
        row['independent_integrator_relative_spread']=abs(radau_coefficient/row['delay_over_epsilon_squared']-1)
    equal = evolve(parent=p,clock='proper',endpoint=tend,eigenvalues=(3.,3.))
    controls['equal_eigenvalue_max_current']=float(max(abs(equal.y[6])))
    light = [direct_light_control(p,proper,epsilon=e) for e in (.1,.05)]
    light.append(direct_light_control(p,proper,epsilon=.05,frame_amplitude=.03))
    emitted,quadrature_error=quad(lambda t:horizon_flux(t,proper.sol(t),p)[1],0.,tend,
                                   epsabs=2e-12,epsrel=2e-11,limit=600)
    EK0=static_patch_energy(0.,proper.y[:,0],p)
    EKend=static_patch_energy(tend,proper.y[:,-1],p)
    EKend96=static_patch_energy(tend,proper.y[:,-1],p,order=96)
    exchange=dict(horizon_areal_radius=1/p.H,
                  horizon_shell_tension=0,initial_Killing_energy_over_epsilon_squared=EK0,
                  final_Killing_energy_over_epsilon_squared=EKend,
                  outward_Killing_energy_over_epsilon_squared=emitted,
                  depleted_initial_patch_fraction=emitted/EK0,
                  energy_balance_relative_residual=abs(EKend+emitted-EK0)/EK0,
                  flux_quadrature_error_over_initial_energy=quadrature_error/EK0,
                  volume_quadrature_relative_spread=abs(EKend96-EKend)/EK0,
                  minimum_sampled_null_flux=min(horizon_flux(t,proper.sol(t),p)[0] for t in np.linspace(0.,tend,601)),
                  region_interpretation='Energy leaves one observer static patch and continues on the same smooth parent. A new child is not constructed.')
    a = 2*math.pi*p.offset
    defect=4*math.sin(a)**2*math.cos(a)*(1-math.cos(a))
    return dict(schema='nsc-atlas-parent-circulation-v1',
        background=dict(spatial_geometry='Smooth round S3; time-symmetric inherited late-vacuum slice',
                        vacuum=p.vacuum,H=p.H,mass_squared=p.mass_squared,kappa=p.kappa,
                        offset=p.offset,finite_proper_endpoint=tend,phase_passes=300,
                        SI_ruler_selected=False),
        preparation=dict(field='One existing real canonical scalar, two coherent spatial harmonics',
                         harmonic_eigenvalues=[3,8],initial_q=[p.initial()[2],p.initial()[4]],
                         initial_velocities=[0,0],initial_mean_energy_over_vacuum='epsilon^2',
                         unique_inheritance_selection=False),
        generator=dict(first_distance_square_defect=defect,
                       drawing_update_is_round_isometry=False,
                       physical_phase_synchronization='One candidate event advances the first mode phase by 2*pi*alpha.',
                       physical_time_origin_proved=False),
        exact_relations=dict(vector_eigenvalue=-7,vector_constraint_eigenvalue=-5,
                             W_norm_squared=16/(3*math.pi**2),
                             current='J=p1*q2-q1*p2; Jdot+3hJ=5q1q2/a^2',
                             vector_response='beta=kappa*J/5',
                             total_rotation_charges=[0]*6,
                             scalar_internal_charge=0,
                             baryon_asymmetry_prediction=None),
        initial_constraint=initial_constraint(p),controls=controls,light_comparisons=rows,
        direct_null_controls=light,horizon_exchange=exchange,
        scope=['Leading epsilon^2 vector Einstein-scalar response on a smooth evolving closed parent',
               'Symmetry-restricted exact small-amplitude constraint-data branch; explicit quadratic corrector',
               'Phase-averaged differential test-guide proper-time record, in model units',
               'Fourth-order metric response, finite-amplitude error bound, selected preparation and child formation remain open'],
        no_observational_fit=True)
