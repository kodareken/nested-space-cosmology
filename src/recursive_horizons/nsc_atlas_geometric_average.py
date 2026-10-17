"""Scalar averages of the owned quiet cap/wall/bridge initial slice.

The result is an instantaneous scalar fitting model in a specified normal
frame. It neither averages a four-metric nor selects a cosmological clock.
All dimensional inputs use the existing geometric G=c=1 convention.
"""
from __future__ import annotations

import math
import time

import numpy as np
from scipy.integrate import quad_vec

from .nsc_atlas_connection_slice import vacuum_transition, wall_slice_data


KAPPA = 8*math.pi


def scalar_average_sources(volume, theta_integral, theta_squared_integral,
                           shear_squared_integral, curvature_integral,
                           energy_integral, stress_trace_integral):
    """Average the Hamiltonian and instantaneous geodesic Raychaudhuri law.

    Stress_trace is the spatial stress trace S, not four-dimensional T.
    Curvature must include surface distributions. A dot denotes the proper
    time of the hypersurface normal, with the domain transported along it.
    """
    values = np.asarray((volume, theta_integral, theta_squared_integral,
                         shear_squared_integral, curvature_integral,
                         energy_integral, stress_trace_integral), dtype=float)
    if not np.isfinite(values).all() or volume <= 0:
        raise ValueError('finite integrals and positive volume are required')
    if theta_squared_integral < 0 or shear_squared_integral < 0:
        raise ValueError('squared kinematic integrals must be nonnegative')
    theta = theta_integral/volume
    variance = theta_squared_integral/volume-theta**2
    if variance < -1e-12*max(1., theta_squared_integral/volume, theta**2):
        raise ValueError('the integrals violate the variance inequality')
    Q = 2*variance/3-2*shear_squared_integral/volume
    R = curvature_integral/volume
    kE, kS = KAPPA*energy_integral/volume, KAPPA*stress_trace_integral/volume
    H = theta/3
    k_rho_g, k_p_g = -(Q+R)/2, -Q/2+R/6
    acceleration = -(kE+kS)/2+Q  # 3 a_ddot/a, not a_ddot/a.
    return {
        'H_normal': H, 'Q': Q, 'mean_intrinsic_curvature': R,
        'kappa_mean_energy': kE, 'kappa_mean_pressure': kS/3,
        'kappa_geometric_density': k_rho_g,
        'kappa_geometric_pressure': k_p_g,
        'geometric_trace_projection': (Q-R)/4,
        'geometric_vacuum_condition': 3*Q+R,
        'three_volume_acceleration': acceleration,
        'total_trace_projection': (kE-kS+Q-R)/4,
        'hamiltonian_residual': 3*H**2-kE+(Q+R)/2,
        'scalar_fitting_acceleration_residual': (
            acceleration+(kE+k_rho_g+kS+3*k_p_g)/2),
    }


def wall_slice_stress(wall):
    """Coefficients of delta(l-lw), including delta(eta)=delta(l-lw)/gamma."""
    sigma, gamma = wall['parameters'].sigma, wall['gamma']
    return {'energy': sigma*gamma,
            'radial_stress': sigma*(gamma**2-1)/gamma,
            'angular_stress': -sigma/gamma,
            'spatial_stress_trace': sigma*(gamma**2-3)/gamma}


def _bridge_primitive(u):
    # Integral cosh(u)^6 du, independently usable against radial quadrature.
    return (5*u/16+15*math.sinh(2*u)/64+3*math.sinh(4*u)/64
            +math.sinh(6*u)/192)


def geometric_slice_average(*, length_scale=1., quadrature='gauss64'):
    """One finite domain: centre to the right-end sphere Rcut/ell=4.

    Alternatives are quadrature controls at this same domain, not a family
    of geometries. Budget: 30 CPU seconds and 20,000 sampled spatial points.
    The stress and curvature delta are integrated exactly before averaging.
    """
    wall = wall_slice_data(length_scale)
    ell, params = wall['ell'], wall['parameters']
    Rw, R0, M = wall['R'], wall['taper_stop'], params.mass
    started, evaluations = time.process_time(), 0

    def integrate(function, a, b):
        def bounded(x):
            nonlocal evaluations
            evaluations += np.asarray(x).size
            if evaluations > 20000 or time.process_time()-started > 30:
                raise RuntimeError('slice-average quadrature budget exhausted')
            return function(x)
        if quadrature in ('gauss32', 'gauss64'):
            n = int(quadrature[5:])
            nodes, weights = np.polynomial.legendre.leggauss(n)
            x = (a+b)/2+(b-a)*nodes/2
            return (b-a)/2*np.einsum('i,i...->...', weights, bounded(x))
        if quadrature == 'adaptive':
            # Work in ell units so the tolerance is invariant under unit changes.
            return quad_vec(bounded, a, b, epsabs=1e-11, epsrel=1e-12)[0]
        raise ValueError('choose gauss32, gauss64 or adaptive quadrature')

    def transition(x):
        fields = vacuum_transition(np.asarray(x)*ell, wall)
        R = fields['R']/ell
        V, VR = fields['V'], fields['V_R']*ell
        theta = VR+2*V/R
        shear_squared = (VR-V/R)**2/3
        measure = 4*math.pi*R**2/(-fields['Rprime'])
        return np.stack((measure, measure*theta, measure*theta**2,
                         measure*shear_squared,
                         measure*fields['scalar_curvature']*ell**2), axis=-1)

    trans = integrate(transition, R0/ell, Rw/ell)
    u0 = -math.acosh(math.sqrt(R0/(2*M)))
    ucut = math.acosh(math.sqrt(4*ell/(2*M)))
    bridge_volume = float(integrate(
        lambda u: 4*math.pi*2*(2*M/ell*np.cosh(u)**2)**3, u0, ucut))
    bridge_exact = 64*math.pi*(M/ell)**3*(_bridge_primitive(ucut)-_bridge_primitive(u0))
    cap_volume = 4*math.pi*(Rw/ell)**3/3
    H = params.H*ell
    volume = cap_volume+float(trans[0])+bridge_volume
    Itheta = 3*H*cap_volume+float(trans[1])
    Itheta2 = 9*H**2*cap_volume+float(trans[2])
    Ishear2 = float(trans[3])
    curvature_delta = -16*math.pi*(Rw/ell)*(wall['Rprime_plus']-1.)
    IR = float(trans[4])+curvature_delta
    surface = wall_slice_stress(wall)
    area = 4*math.pi*(Rw/ell)**2
    Ecap = 3*H**2/KAPPA*cap_volume
    Ewall = area*surface['energy']*ell
    Swall = area*surface['spatial_stress_trace']*ell
    sources = scalar_average_sources(volume, Itheta, Itheta2, Ishear2, IR,
                                     Ecap+Ewall, -3*Ecap+Swall)
    if time.process_time()-started > 30:
        raise RuntimeError('slice-average CPU budget exhausted')
    return {
        'schema': 'nsc-atlas-geometric-slice-average-v1',
        'domain': {'Rwall_over_ell': 1.5, 'Rcut_right_over_ell': 4.,
                   'H_ell': H, 'k_ell': params.k*ell, 'M_over_ell': M/ell},
        'clock': 'Instantaneous geodesic hypersurface normal; not an observer light-cone map.',
        'volume_over_ell3': volume,
        'cap_volume_fraction': cap_volume/volume,
        'integrals_in_ell_units': {
            'theta': Itheta, 'theta_squared': Itheta2,
            'shear_squared': Ishear2, 'curvature_bulk': float(trans[4]),
            'curvature_wall_delta': curvature_delta, 'curvature_total': IR,
            'energy_cap': Ecap, 'energy_wall': Ewall,
            'stress_trace_cap': -3*Ecap, 'stress_trace_wall': Swall,
        },
        'normalized_sources': sources,  # H*ell, every other entry *ell^2.
        'checks': {
            'hamiltonian_ell2': sources['hamiltonian_residual'],
            'scalar_fitting_acceleration_ell2': sources['scalar_fitting_acceleration_residual'],
            'bridge_volume_ell3_difference': bridge_volume-bridge_exact,
            'wall_trace_ell': (surface['spatial_stress_trace']-surface['energy']
                               +3*params.sigma/wall['gamma'])*ell,
        },
        'quadrature': quadrature, 'spatial_evaluations': evaluations,
        'scope': 'Scalar fitting coefficients of one fixed sourced initial slice; no four-tensor vacuum or cosmological prediction.',
    }
