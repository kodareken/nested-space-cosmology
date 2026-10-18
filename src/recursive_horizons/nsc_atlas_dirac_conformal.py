"""Massless Dirac calculation variables at the quiet pure-tension interface.

The original gravitational metric and observer clocks are unchanged. This
is a local, nonsingular conformal spinor calculation at r>0, not a state,
one-loop Weyl invariance assertion, reflection coefficient or slip value.
"""
from __future__ import annotations

from math import isfinite

from .nsc_atlas_cap_wall_response import reference_coefficients
from .nsc_atlas_dirac_probe import reference_dirac_interface


def radial_product_curvature(*, radius, f, f_r, f_rr, sphere_radius=1.0):
    """Two-dimensional curvature after g_tilde=(L/r)^2*g.

    On the one-sided smooth static radial metric -f dt2+dr2/f,
    R2=-f_rr and box2(log r)=f_r/r-f/r2. No division by f remains;
    the resulting invariant continues through f=0 in regular coordinates.
    """
    r, L = float(radius), float(sphere_radius)
    vals = tuple(float(v) for v in (f, f_r, f_rr))
    if (not all(isfinite(v) for v in (r, L, *vals)) or r <= 0 or L <= 0):
        raise ValueError('require positive finite r,L and finite one-sided metric jets')
    f, fp, fpp = vals
    radial = (-r*r*fpp+2*r*fp-2*f)/(L*L)
    return dict(radial_scalar_curvature=radial,
        four_dimensional_scalar_curvature=radial+2/(L*L))


def conformal_interface_data(radius):
    """Transform the existing reference domain by Omega=ell/r, ell=1.

    Ktilde^a_b=Omega^-1*(K^a_b+n(log Omega)*delta^a_b).
    psi_tilde=Omega^-3/2*psi; ntilde=Omega^-1*n. Return coefficients
    multiplying psi_tilde in its normal derivative jump, not amplitudes.
    """
    old = reference_dirac_interface(radius)
    geometry = reference_coefficients(radius)
    r = float(old['radius'])
    normal_log_r_jump = (old['beta_plus']-old['beta_minus'])/r
    normal_log_Omega_jump = -normal_log_r_jump
    trace_K_jump = r*(old['trace_K_jump']+3*normal_log_Omega_jump)
    spinor_jump = r*(old['squared_Dirac_normal_derivative_jump_coefficient']
                    -1.5*normal_log_Omega_jump)
    return dict(radius=r, Omega=1/r,
        original_normal_log_radius_jump=normal_log_r_jump,
        normal_log_Omega_jump=normal_log_Omega_jump,
        original_trace_K_jump=old['trace_K_jump'],
        conformal_trace_K_jump=trace_K_jump,
        conformal_angular_K_jump=r*(normal_log_r_jump+normal_log_Omega_jump),
        conformal_temporal_K_jump=r*(float(geometry['beta_plus_x']-geometry['beta_minus_x'])
                                    +normal_log_Omega_jump),
        conformal_normal_spinor_derivative_jump=spinor_jump,
        conformal_spinor_is_continuous=True,
        transformation_regular_at_material_wall=True,
        transformation_regular_at_centre=False,
        global_state_selected=False,
        is_quantum_conformal_invariance_claim=False)


def reference_conformal_report(radius=.45):
    r, H, M = float(radius), 1., .4
    interface = conformal_interface_data(r)
    cap = radial_product_curvature(radius=r, f=1-H*H*r*r,
        f_r=-2*H*H*r, f_rr=-2*H*H)
    exterior = radial_product_curvature(radius=r, f=1-2*M/r,
        f_r=2*M/r**2, f_rr=-4*M/r**3)
    return dict(inputs=dict(Rwall_over_ell=r, H_ell=H, M_over_ell=M,
            k_wall_ell=2., sphere_radius_over_ell=1., scalar_amplitude=0., Dirac_mass=0.),
        conformal_interface=interface, cap=cap, exterior=exterior,
        radial_curvature_jump=exterior['radial_scalar_curvature']-cap['radial_scalar_curvature'],
        radial_curvature_jump_control=12*M/r,
        cap_Weyl_squared=0., exterior_Weyl_squared=48*M*M/(r*r),
        smooth_cap_control_radial_curvature_jump=0.,
        smooth_cap_control_interface_derivative_jump=0.,
        counts=dict(algebraic_wall_points=1, Dirac_time_steps=0, P001_time_steps=0,
                    new_evolution_grids=0, loop_polarizations=0),
        nested_Psi_over_Phi=None,
        scope='Massless classical operator and domain in local conformal variables; physical clocks remain in g, and quantum state/anomaly/metric response are not computed.')
