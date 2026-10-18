"""Quiet source-anchor slice for the joined cap/wall Dirac calculation.

This supplies spatial geometry and its intrinsic Dirac coefficients. It is
not a physical time Hamiltonian, two-point state or one-loop metric response.
The older flat-cap slice at Rwall/ell=1.5 remains in connection_slice.
"""
from __future__ import annotations

import math

import numpy as np

from .nsc_atlas_connection_slice import smooth_step
from .nsc_atlas_moving_boundaries import VacuumWallParameters, vacuum_wall_state


def source_anchor_slice(length_scale=1.0):
    """Fixed quiet model at Rwall/ell=.45, with the cap static slice inside."""
    ell = float(length_scale)
    if not math.isfinite(ell) or ell <= 0:
        raise ValueError('length_scale must be finite and positive')
    p = VacuumWallParameters(H=1/ell, k=2/ell, mass=.4*ell)
    rw = .45*ell
    state = vacuum_wall_state(rw, p)
    speed, bm, bp = (float(state[name]) for name in ('Rdot', 'beta_minus', 'beta_plus'))
    root = math.sqrt(1-(p.H*rw)**2)
    gamma, gamma_v = bm/root, speed/root
    J = gamma*bp-gamma_v*speed
    E = gamma*speed-gamma_v*bp
    a = math.sqrt(1+E*E)
    rmin = 2*p.mass/(a*a)
    if not (gamma >= 1 and gamma_v > 0 and J < 0 and E > 0 and 0 < rmin < rw):
        raise ValueError('source anchor does not have the declared slice orientation')
    qw = -math.acosh(math.sqrt(rw/rmin))
    lw = math.asin(p.H*rw)/p.H
    lmin = lw-rmin/a*(qw+math.sinh(qw)*math.cosh(qw))
    return dict(ell=ell, parameters=p, R=rw, Rdot=speed,
        beta_minus=bm, beta_plus=bp, gamma=gamma, gamma_v=gamma_v,
        Rprime_minus=root, Rprime_plus=J, V_minus=0., V_plus=E,
        minimum_radius=rmin, minimum_distance=lmin, wall_distance=lw,
        q_wall=qw, taper_start=ell, taper_stop=1.2*ell,
        energy_surface=p.sigma*gamma, momentum_surface=p.sigma*gamma_v)


def cap_slice(proper_distance, wall):
    l = np.asarray(proper_distance, dtype=float)
    if not np.all(np.isfinite(l)) or np.any(l < 0) or np.any(l > wall['wall_distance']):
        raise ValueError('cap distance must lie between centre and wall')
    H = wall['parameters'].H
    R = np.sin(H*l)/H
    return dict(l=l, R=R, Rprime=np.cos(H*l), Rsecond=-H*H*R,
        V=np.zeros_like(R), Kll=np.zeros_like(R), Ktheta=np.zeros_like(R),
        scalar_curvature=np.full_like(R, 6*H*H), mass=H*H*R**3/2)


def _vacuum_fields(R, Rp, V, V_R, wall):
    M = wall['parameters'].mass
    Rpp = M/R**2+V*V_R
    Kll, Ktheta = -V_R, -V/R
    Ktheta_l = (-V_R/R+V/R**2)*Rp
    scalar = -4*Rpp/R+2*(1-Rp*Rp)/R**2
    # Outgoing EF: ds2=-f du2-2 du dR. These expressions remain finite at f=0.
    f = 1-2*M/R
    ef = 1/(V+Rp)
    spatial_norm = -f*ef**2+2*ef*Rp
    normal_norm = -f*ef**2-2*ef*V
    orthogonality = f*ef**2+(V-Rp)*ef
    return dict(R=R, Rprime=Rp, Rsecond=Rpp, V=V, V_R=V_R,
        Kll=Kll, Ktheta=Ktheta, scalar_curvature=scalar,
        mass=R/2*(1+V*V-Rp*Rp),
        hamiltonian_residual=scalar+4*Kll*Ktheta+2*Ktheta*Ktheta,
        momentum_residual=2*Rp/R*(Kll-Ktheta)-2*Ktheta_l,
        outgoing_EF_u_l=-ef, outgoing_EF_normal_u=ef,
        embedding_spatial_norm=spatial_norm, embedding_normal_norm=normal_norm,
        embedding_orthogonality=orthogonality)


def constant_rate_bridge(q, wall):
    """Proper-distance regular Schwarzschild segment, including its minimum."""
    q = np.asarray(q, dtype=float)
    rmin, E = wall['minimum_radius'], wall['V_plus']
    qend = math.acosh(math.sqrt(wall['taper_start']/rmin))
    if (not np.all(np.isfinite(q)) or np.any(q < wall['q_wall'])
            or np.any(q > qend)):
        raise ValueError('bridge parameter is outside wall-to-taper domain')
    a = math.sqrt(1+E*E)
    R = rmin*np.cosh(q)**2
    fields = _vacuum_fields(R, a*np.tanh(q), np.full_like(R, E), np.zeros_like(R), wall)
    fields['l'] = wall['minimum_distance']+rmin/a*(q+np.sinh(q)*np.cosh(q))
    return fields


def outgoing_slice(radius, wall):
    """Vacuum taper on the increasing-radius branch, followed by an exact AE end."""
    R = np.asarray(radius, dtype=float)
    if not np.all(np.isfinite(R)) or np.any(R < wall['taper_start']):
        raise ValueError('outgoing radius must start at or beyond the taper')
    width = wall['taper_stop']-wall['taper_start']
    step, derivative = smooth_step((R-wall['taper_start'])/width)
    V = wall['V_plus']*(1-step)
    V_R = -wall['V_plus']*derivative/width
    Rp = np.sqrt(1-2*wall['parameters'].mass/R+V*V)
    return _vacuum_fields(R, Rp, V, V_R, wall)


def spatial_dirac_coefficients(radius, radial_slope, *, angular_k=1):
    """Intrinsic D_Sigma,k=-i sigma2 d/dl+sigma1*k/r, in norm int |u|2 dl.

    u=r*psi is a spatial density rescaling. This is not the optical-cap
    Hamiltonian, the full spacetime squared operator or a state prescription.
    """
    R, Rp = np.broadcast_arrays(np.asarray(radius, dtype=float),
                                np.asarray(radial_slope, dtype=float))
    if (not np.all(np.isfinite(R)) or not np.all(np.isfinite(Rp)) or np.any(R <= 0)
            or isinstance(angular_k, bool) or not isinstance(angular_k, int) or not angular_k):
        raise ValueError('require finite positive radius, finite slope and nonzero integer k')
    W = angular_k/R
    W_l = -angular_k*Rp/R**2
    return dict(W=W, W_l=W_l, V_plus=W*W+W_l, V_minus=W*W-W_l)


def _outgoing_distance(wall, order, right_radius):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    start, stop = wall['taper_start'], wall['taper_stop']
    R = (start+stop)/2+(stop-start)/2*nodes
    taper = (stop-start)/2*float(weights @ (1/outgoing_slice(R, wall)['Rprime']))
    M = wall['parameters'].mass
    def primitive(r):
        return math.sqrt(r*(r-2*M))+2*M*math.acosh(math.sqrt(r/(2*M)))
    return taper+primitive(right_radius)-primitive(stop)


def reference_slice_report():
    wall = source_anchor_slice()
    p, E = wall['parameters'], wall['V_plus']
    rmin, rw = wall['minimum_radius'], wall['R']
    qend = math.acosh(math.sqrt(wall['taper_start']/rmin))
    qh = math.acosh(math.sqrt(2*p.mass/rmin))
    qs = np.concatenate((np.linspace(wall['q_wall'], qend, 65), [0., qh]))
    bridge = constant_rate_bridge(qs, wall)
    taper = outgoing_slice(np.linspace(wall['taper_start'], wall['taper_stop'], 65), wall)
    jump_rp = wall['Rprime_plus']-wall['Rprime_minus']
    jump_Ktheta = -E/rw
    # Independent tensor projection S_ab=-sigma*(g_ab-n_a*n_b), restricted
    # to the slice with delta(n)=delta(l-lw)/gamma.
    g = np.diag([-1., 1.])
    n = np.array([-wall['gamma_v'], wall['gamma']])
    stress = -p.sigma*(g-np.outer(n, n))/wall['gamma']
    checks = dict(
        boost_unit_error=abs(wall['gamma']**2-wall['gamma_v']**2-1),
        exterior_mass_identity_error=abs(wall['Rprime_plus']**2-E*E-(1-2*p.mass/rw)),
        shell_hamiltonian_error=abs(-4*jump_rp/rw-16*math.pi*wall['energy_surface']),
        shell_momentum_error=abs(-2*jump_Ktheta-8*math.pi*wall['momentum_surface']),
        surface_energy_projection_error=abs(stress[0, 0]-wall['energy_surface']),
        surface_momentum_projection_error=abs(-stress[0, 1]-wall['momentum_surface']))
    for name, data in (('bridge', bridge), ('taper', taper)):
        for key, target in (('hamiltonian_residual', 0), ('momentum_residual', 0),
                ('mass', p.mass), ('embedding_spatial_norm', 1),
                ('embedding_normal_norm', -1), ('embedding_orthogonality', 0)):
            checks[name+'_'+key+'_max_error'] = float(np.max(np.abs(data[key]-target)))
    wall_field = constant_rate_bridge(np.array([wall['q_wall']]), wall)
    checks['wall_overlap_error'] = max(abs(float(wall_field['R'][0])-rw),
        abs(float(wall_field['Rprime'][0])-wall['Rprime_plus']),
        abs(float(wall_field['l'][0])-wall['wall_distance']))
    before = constant_rate_bridge(np.array([qend]), wall)
    after = outgoing_slice(np.array([wall['taper_start']]), wall)
    checks['taper_join_max_error'] = max(abs(float(before[key][0]-after[key][0]))
        for key in ('R', 'Rprime', 'Rsecond', 'V', 'V_R', 'Kll', 'Ktheta'))
    dist32 = _outgoing_distance(wall, 32, 4.)
    dist64 = _outgoing_distance(wall, 64, 4.)
    right_l = float(before['l'][0])+dist64
    cap_coeff = spatial_dirac_coefficients(rw, wall['Rprime_minus'])
    ext_coeff = spatial_dirac_coefficients(rw, wall['Rprime_plus'])
    return dict(inputs=dict(ell=1., H_ell=1., k_wall_ell=2., M_over_ell=.4,
            Rwall_over_ell=.45, scalar_amplitude=0., taper_R_over_ell=[1., 1.2],
            right_display_R_over_ell=4.),
        wall={key: wall[key] for key in ('Rdot', 'beta_minus', 'beta_plus', 'gamma',
            'gamma_v', 'Rprime_minus', 'Rprime_plus', 'V_minus', 'V_plus', 'wall_distance')},
        geometry=dict(minimum_R_over_ell=rmin, minimum_l_over_ell=wall['minimum_distance'],
            event_horizon_R_over_ell=2*p.mass, right_display_l_over_ell=right_l,
            radial_order='regular centre, cap, wall, minimum sphere, horizon, taper, one AE end',
            minimum_is_event_horizon=False, analytic_AE_ends=1,
            vacuum_data_beyond_taper='Schwarzschild, K=0, ADM mass .4*ell'),
        intrinsic_spatial_Dirac=dict(angular_k=1,
            W_wall=float(cap_coeff['W']), W_l_minus=float(cap_coeff['W_l']),
            W_l_plus=float(ext_coeff['W_l']),
            W_l_jump=float(ext_coeff['W_l']-cap_coeff['W_l']),
            V_plus_minus=float(cap_coeff['V_plus']), V_plus_plus=float(ext_coeff['V_plus']),
            V_minus_minus=float(cap_coeff['V_minus']), V_minus_plus=float(ext_coeff['V_minus']),
            W_l_jump_control=p.k*wall['gamma']/rw,
            square_has_interface_delta=False, square_has_interface_step=True,
            is_four_dimensional_Dirac_squared=False, is_physical_time_Hamiltonian=False),
        checks=checks,
        numerical_comparison=dict(outgoing_distance_32_64_abs_difference=abs(dist32-dist64),
            is_rigorous_error_enclosure=False),
        counts=dict(geometric_points=int(1+qs.size+65+2+32+64),
            P001_evolution_steps=0, new_P001_grids=0, loop_polarizations=0),
        nested_Psi_over_Phi=None,
        scope='Quiet sourced initial surface and intrinsic spatial operator only; no physical quantum state, metric response or formation claim.')
