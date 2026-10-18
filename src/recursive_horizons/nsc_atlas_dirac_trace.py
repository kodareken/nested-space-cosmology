"""Local massless Dirac bulk and Euler-surface terms on the actual quiet atlas.

Free four-component Dirac quantization and the Laplace-type heat coefficient
are imported. No singular wall curvature is squared. A trace projection is
not the full stress, a state, gravitational slip or self-consistent Lambda.
"""
from __future__ import annotations

import math

from .nsc_atlas_dirac_probe import quiet_reference_curvatures
from .nsc_atlas_cap_wall_response import reference_coefficients


def dirac_bulk_trace(*, R, Ricci_squared, Riemann_squared):
    """Return <T^mu_mu>/(hbar*c) for a smooth constant-R region, in length^-4.

    For L=D^2=-(nabla^2+E), E=-R/4 and tr I=4,
    tr Omega_mn Omega^mn=-Riemann^2/2. The traced heat numerator is
    5R^2-8Ricci^2-7Riemann^2. The fermionic determinant reverses its sign.
    The total-derivative term is omitted only because Box R=0 here.
    """
    R, Ric, Riem = (float(v) for v in (R, Ricci_squared, Riemann_squared))
    if not all(math.isfinite(v) for v in (R, Ric, Riem)):
        raise ValueError('smooth bulk curvature invariants must be finite')
    return (7*Riem+8*Ric-5*R*R)/(360*(4*math.pi)**2)


def reference_bulk_trace_report(radius=.45, *, ell_m=10000.):
    r, ell = float(radius), float(ell_m)
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError('fiducial physical length must be positive and finite')
    geometry = quiet_reference_curvatures(r)
    coefficients = {}
    for name, row in geometry.items():
        shape = dirac_bulk_trace(R=row['R'], Ricci_squared=row['Ricci_squared'],
                                Riemann_squared=row['Riemann_squared'])
        coefficients[name] = dict(trace_times_ell4_over_hbar_c=shape,
            delta_Lambda_trace_shape=-2*math.pi*shape)
    # Same supplied physical constants and 10 km fiducial used by the scalar
    # trace/P001 scale calculations. h and c exact; G is supplied measured input.
    h_SI, c_SI, G_SI = 6.62607015e-34, 299792458., 6.67430e-11
    lP2 = (h_SI/(2*math.pi))*G_SI/c_SI**3
    zeta = lP2/ell**2
    for row in coefficients.values():
        row['delta_Lambda_trace_times_ell2'] = row['delta_Lambda_trace_shape']*zeta
        row['delta_Lambda_trace_per_m2'] = row['delta_Lambda_trace_times_ell2']/ell**2
    return dict(inputs=dict(Rwall_over_ell=r, H_ell=1., M_over_ell=.4,
            k_wall_ell=2., Dirac_mass=0., four_component_Dirac_fields=1,
            ell_m=ell, ell_is_predicted=False, G_SI=G_SI, h_SI=h_SI, c_SI=c_SI,
            Planck_length_squared_over_ell_squared=zeta,
            independent_new_quantum_state_is_selected=False),
        one_sided=coefficients,
        checks=dict(cap_shape_error=abs(coefficients['cap']['delta_Lambda_trace_shape']-11/(120*math.pi)),
            exterior_shape_error=abs(coefficients['exterior']['delta_Lambda_trace_shape']
                                    +7*.4**2/(60*math.pi*r**6))),
        counts=dict(algebraic_wall_points=1, local_one_loop_bulk_coefficients=2,
            nonlocal_loop_polarizations=0, P001_time_steps=0, Dirac_time_steps=0,
            new_evolution_grids=0),
        nested_Psi_over_Phi=None,
        scope='One smooth-bulk trace term of the declared Dirac probe on the original quiet metric. No interface distribution products, full vacuum stress, state, nonlocal polarization, backreaction or derived cosmological Lambda.',
        uncertainty_scope='Algebraic curvature coefficient conditional on imported free Dirac theory; supplied G and fiducial length, not a total physical uncertainty or a converged full quantum model.')


def euler_chern_scalar(*, Ktau, Ksphere, intrinsic_acceleration, intrinsic_curvature):
    """4*(J-2*Ghat.K) for induced -dtau^2+R(tau)^2*dOmega^2.

    Intrinsic inputs are Rddot/R and (1+Rdot^2)/R^2. All K use ONE
    common normal, not the two outward normals of the disjoint regions.
    """
    t, b, A, B = (float(v) for v in (Ktau, Ksphere,
                                     intrinsic_acceleration, intrinsic_curvature))
    if not all(math.isfinite(v) for v in (t, b, A, B)):
        raise ValueError('finite common-normal and intrinsic wall curvatures required')
    return -8*t*b*b+8*B*t+16*A*b


def quiet_euler_surface(radius):
    """Euler coefficient multiplying invariant proper-normal delta_Sigma.

    Derived by integrating a smooth normal jet, not by multiplying
    distributions. Conditional on bounded conformal-product curvature in
    the declared conformal smoothing class. This is not the full anomaly.
    """
    R = float(radius)
    g = {key: float(value) for key, value in reference_coefficients(R).items()}
    A = g['radial_acceleration']/R
    B = (1+g['v']**2)/R**2
    chern = {}
    for side in ('minus', 'plus'):
        chern[side] = euler_chern_scalar(Ktau=g['beta_'+side+'_x'],
            Ksphere=g['beta_'+side]/R, intrinsic_acceleration=A, intrinsic_curvature=B)
    return dict(euler_surface=chern['minus']-chern['plus'],
        chern_minus=chern['minus'], chern_plus=chern['plus'],
        intrinsic_acceleration=A, intrinsic_curvature=B)


def dirac_euler_surface_trace(euler_surface):
    """Type-A trace/(hbar*c), coefficient of delta_Sigma, in length^-3.

    Excludes Box R and finite interface terms; not a bulk trace or full
    vacuum stress. An arbitrary microscopic wall is not assumed equivalent.
    """
    E = float(euler_surface)
    if not math.isfinite(E):
        raise ValueError('finite Euler surface coefficient required')
    return -11*E/(360*(4*math.pi)**2)


def reference_euler_wall_report(radius=.45, *, ell_m=10000.):
    ell = float(ell_m)
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError('positive finite supplied physical scale required')
    wall = quiet_euler_surface(radius)
    coefficient = dirac_euler_surface_trace(wall['euler_surface'])
    shape = -2*math.pi*coefficient
    h_SI, c_SI, G_SI = 6.62607015e-34, 299792458., 6.67430e-11
    lP2 = h_SI*G_SI/(2*math.pi*c_SI**3)
    zeta = lP2/ell**2
    wall.update(trace_times_ell3_over_hbar_c=coefficient,
        trace_delta_coefficient_per_m3=coefficient/ell**3,
        delta_Lambda_trace_surface_shape=shape,
        ell_times_integrated_delta_Lambda_trace=shape*zeta,
        integrated_delta_Lambda_trace_per_m=shape*zeta/ell)
    return dict(inputs=dict(Rwall_over_ell=float(radius), H_ell=1., k_wall_ell=2.,
            M_over_ell=.4, Dirac_mass=0., four_component_Dirac_fields=1,
            ell_m=ell, ell_is_predicted=False, h_SI=h_SI, G_SI=G_SI, c_SI=c_SI,
            Planck_length_squared_over_ell_squared=zeta,
            normal_orientation='cap to exterior on both sides',
            smoothing_class='smooth kinked conformal factor; conformal product curvature remains bounded'),
        wall=wall, smooth_cap_continuation_euler_surface=0.,
        checks=dict(signed_chern_vs_closed_junction=abs(wall['euler_surface']-56.),
            trace_coefficient_vs_closed_junction=abs(coefficient+77/(720*math.pi**2))),
        excluded_terms=dict(Box_R='scheme dependent; not set to zero at the wall',
            finite_wall_counterterms='not selected', state_dependent_traceless_stress='not calculated'),
        counts=dict(algebraic_wall_points=1, type_A_surface_coefficients=1,
            new_evolution_grids=0, Dirac_time_steps=0, P001_time_steps=0,
            nonlocal_loop_polarizations=0), nested_Psi_over_Phi=None,
        scope='Euler part of the actual material-wall Dirac anomaly in the declared smoothing class. Not full surface stress, bulk Lambda, a microscopic wall model, nonlocal response or a slip prediction.',
        uncertainty_scope='Exact algebraic coefficient within stated wall class; no physical state, full quantum uncertainty or wall-thickness extrapolation is supplied.')
