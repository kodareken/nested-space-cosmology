"""One-sided time generator and compatible initial moving-wall spinor jets.

The declared lapse/frame differs from the intrinsic spatial operator. These
local identities supply no full foliation, quantum state or metric loop.
"""
from __future__ import annotations

import numpy as np

from .nsc_atlas_dirac_slice import source_anchor_slice
from .nsc_atlas_moving_boundaries import vacuum_wall_state


ALPHA = np.array([[0., -1j], [1j, 0.]])


def radial_time_coefficients(*, lapse, lapse_l, radial_length, radial_length_l,
                             radius, angular_k=1):
    """H=-i alpha*(c*d_l+c_l/2)+sigma1*N*k/r for chi=sqrt(a)*r*psi.

    The smooth one-sided metric is -N^2 dt^2+a^2 dl^2+r^2 dOmega^2,
    shift zero, c=N/a, with a nonrotating spherical frame. Derivatives are
    coordinate l derivatives. The half-density norm is integral |chi|^2 dl.
    """
    N, Nl, a, al, r = np.broadcast_arrays(*[np.asarray(v, dtype=float) for v in
        (lapse, lapse_l, radial_length, radial_length_l, radius)])
    if (not all(np.all(np.isfinite(v)) for v in (N, Nl, a, al, r))
            or np.any(N <= 0) or np.any(a <= 0) or np.any(r <= 0)
            or isinstance(angular_k, bool) or not isinstance(angular_k, int) or not angular_k):
        raise ValueError('require positive N,a,r, finite derivatives and a nonzero integer k')
    return dict(speed=N/a, derivative_companion=(Nl/a-N*al/a**2)/2,
        angular_coefficient=N*angular_k/r)


def initial_wall_time_data(wall):
    """Required lapse jets for the checked slice, and their transparent trace jump.

    a=1, Kll=0 on both sides initially; Nwall=1 and shift=0. The cap
    clock is its static clock normalized to the wall at this instant.
    The exterior lapse beyond this first jet is not fixed by these data.
    """
    p, r = wall['parameters'], wall['R']
    gamma, gamma_v = wall['gamma'], wall['gamma_v']
    s = gamma_v/gamma
    a_minus = -p.H*p.H*r/wall['Rprime_minus']
    a_jump = -p.k/gamma
    a_plus = a_minus+a_jump
    data = vacuum_wall_state(r, p)
    theta_tau = float(data['Ktau_minus'])-gamma*a_minus
    # [chi_l] = (identity_coefficient*I+alpha_coefficient*alpha)*chi.
    identity_coefficient = p.k*gamma/2
    alpha_coefficient = p.k*gamma_v/2
    return dict(wall_coordinate_speed_over_N=s,
        d_log_lapse_minus=a_minus, d_log_lapse_plus=a_plus,
        d_log_lapse_jump=a_jump, common_boost_rapidity_tau=theta_tau,
        Ktau_minus=float(data['Ktau_minus']), Ktau_plus=float(data['Ktau_plus']),
        reconstructed_Ktau_minus=theta_tau+gamma*a_minus,
        reconstructed_Ktau_plus=theta_tau+gamma*a_plus,
        density_derivative_jump_identity=identity_coefficient,
        density_derivative_jump_alpha=alpha_coefficient,
        density_derivative_jump_alpha_eigenvalues=[identity_coefficient-alpha_coefficient,
                                                  identity_coefficient+alpha_coefficient],
        normal_component_derivative_jump_identity=3*p.k/2,
        normal_component_derivative_jump_alpha=p.k*s/2,
        normal_spin_connection_jump=-p.k*s,
        parallel_frame_normal_derivative_jump=3*p.k/2,
        global_initial_lapse_is_supplied=False, full_time_development_is_supplied=False)


def moving_wall_density_jump(trace, wall):
    """Initial [d_l chi] implied by moving transparency and the time equation."""
    chi = np.asarray(trace, dtype=complex)
    if chi.shape != (2,) or not np.all(np.isfinite(chi)):
        raise ValueError('trace must be a finite complex radial two-column')
    return wall['parameters'].k/2*(wall['gamma']*chi+wall['gamma_v']*(ALPHA @ chi))


def reference_wall_time_report():
    w = source_anchor_slice()
    d = initial_wall_time_data(w)
    trace = np.array([1+.3j, -.7+.2j])
    delta_l = moving_wall_density_jump(trace, w)
    s = d['wall_coordinate_speed_over_N']
    delta_t = -ALPHA @ (delta_l+d['d_log_lapse_jump']*trace/2)
    p, gamma, gamma_v = w['parameters'], w['gamma'], w['gamma_v']
    # Divide the density jets by r and retain both geometric derivative jumps.
    delta_psi_l = delta_l/w['R']+p.k*gamma*trace/w['R']
    delta_psi_t = delta_t/w['R']-p.k*gamma_v*trace/w['R']
    normal = gamma_v*delta_psi_t+gamma*delta_psi_l
    normal += d['normal_spin_connection_jump']/2*(ALPHA @ trace)/w['R']
    checks = dict(
        temporal_Israel_error=abs(d['reconstructed_Ktau_plus']-d['Ktau_plus']),
        moving_trace_error=float(np.max(np.abs(delta_t+s*delta_l))),
        tangential_component_jump=float(np.max(np.abs(gamma*delta_psi_t+gamma_v*delta_psi_l))),
        parallel_normal_jump_error=float(np.max(np.abs(normal-3*p.k/2*trace/w['R']))))
    return dict(inputs=dict(Rwall_over_ell=.45, H_ell=1., k_wall_ell=2., M_over_ell=.4,
            scalar_amplitude=0., N_wall=1., shift=0., a_initial=1., Kll_initial_both_sides=0.),
        wall=d, checks=checks,
        nested_Psi_over_Phi=None,
        scope='One-sided initial Dirac time generator and moving-wall jet compatibility; no complete foliation, quantum state, loop or evolved spinor.',
        counts=dict(algebraic_wall_points=1, Dirac_time_steps=0, P001_time_steps=0,
            new_evolution_grids=0, loop_polarizations=0))
