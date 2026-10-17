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


def clocked_scalar_average_sources(
        volume, *, lapse_theta_integral, lapse_squared_theta_squared_integral,
        lapse_squared_shear_squared_integral, lapse_squared_curvature_integral,
        lapse_squared_energy_integral, lapse_squared_stress_trace_integral,
        lapse_rate_theta_integral, lapse_laplacian_integral):
    """Scalar fitting in a specified observer clock, with zero shift.

    The domain is transported by the hypersurface normals. ``lapse`` means
    nu = N/N_ref, where d tau_ref = N_ref(t) dt; its rate is at fixed spatial
    coordinate with respect to tau_ref. All integrals use proper spatial
    volume. The last integral is int nu * D^2(nu) dV, including any boundary
    flux. It is not twice that integral, nor a derivative of the mean lapse.

    Einstein equations and the stated foliation give these scalar identities.
    This function does not select that geometry/clock or average a tensor.
    """
    clock_integrals = np.asarray((lapse_rate_theta_integral,
                                  lapse_laplacian_integral), dtype=float)
    if not np.isfinite(clock_integrals).all():
        raise ValueError('finite clock integrals are required')
    row = scalar_average_sources(
        volume, lapse_theta_integral, lapse_squared_theta_squared_integral,
        lapse_squared_shear_squared_integral, lapse_squared_curvature_integral,
        lapse_squared_energy_integral, lapse_squared_stress_trace_integral)
    P = float(clock_integrals.sum()/volume)
    row['H_observer'] = row.pop('H_normal')
    row['lapse_acceleration'] = P
    row['kappa_geometric_pressure'] -= 2*P/3
    row['geometric_trace_projection'] += P/2
    row['geometric_vacuum_condition'] += 2*P
    row['three_volume_acceleration'] += P
    row['total_trace_projection'] += P/2
    return row


def _periodic_derivative(values, coordinate_length, order=1):
    count = values.size
    spectrum = np.fft.fft(values)
    wavenumber = 2*math.pi*np.fft.fftfreq(count, d=coordinate_length/count)
    multiplier = (1j*wavenumber)**order
    return np.fft.ifft(multiplier*spectrum).real


def periodic_spherical_average(
        areal_radius, radial_scale, lapse, radius_rate, scale_rate, lapse_rate,
        radius_acceleration, scale_acceleration, coordinate_length, observer_index):
    """Kinematic average of one periodic spherical metric snapshot.

    The line element is ds^2 = -N^2 dt^2 + a^2 dx^2 + r^2 dOmega^2, with zero
    shift and periodic x. Coordinates t,x are dimensionless, so r, a, N and
    every coordinate rate have length dimension. The clock is
    d tau_o = N[observer_index] dt. This reconstructs the Einstein source from
    the metric, including whatever higher-curvature content that representation
    contains. It is not the Dirac stress, not an on-shell matter check, and
    not the P001 cap-wall parent. A nonzero trace is not a vacuum identification.
    """
    arrays = [np.asarray(values, dtype=float) for values in (
        areal_radius, radial_scale, lapse, radius_rate, scale_rate, lapse_rate,
        radius_acceleration, scale_acceleration)]
    if any(values.ndim != 1 or values.shape != arrays[0].shape for values in arrays):
        raise ValueError('spherical snapshot arrays must be equal-length one-dimensional')
    if isinstance(coordinate_length, bool) or not isinstance(observer_index, (int, np.integer)):
        raise ValueError('coordinate length and observer index must be real and integral')
    length = float(coordinate_length)
    if (not math.isfinite(length) or length <= 0 or any(not np.isfinite(values).all() for values in arrays)
            or np.any(arrays[0] <= 0) or np.any(arrays[1] <= 0) or np.any(arrays[2] <= 0)):
        raise ValueError('positive geometry, finite rates and positive coordinate length are required')
    count = arrays[0].size
    if count < 2 or not 0 <= int(observer_index) < count:
        raise ValueError('observer index must lie on a periodic grid of at least two nodes')
    radius, scale, lapse, radius_rate, scale_rate, lapse_rate, radius_acceleration, scale_acceleration = arrays
    dx = length/count
    measure = 4*math.pi*scale*radius**2
    radius_x = _periodic_derivative(radius, length)
    scale_x = _periodic_derivative(scale, length)
    radius_xx = _periodic_derivative(radius, length, 2)
    observer_lapse = float(lapse[int(observer_index)])
    observer_lapse_rate = float(lapse_rate[int(observer_index)])
    reduced_lapse = lapse/observer_lapse
    reduced_rate = (lapse_rate-reduced_lapse*observer_lapse_rate)/observer_lapse**2
    reduced_x = _periodic_derivative(reduced_lapse, length)
    reduced_xx = _periodic_derivative(reduced_lapse, length, 2)
    theta = (scale_rate/scale+2*radius_rate/radius)/lapse
    shear_squared = (scale_rate/scale-radius_rate/radius)**2/(3*lapse**2)
    curvature = (2/radius**2-4*radius_xx/(scale**2*radius)
                 +4*scale_x*radius_x/(scale**3*radius)-2*radius_x**2/(scale**2*radius**2))
    laplacian = (reduced_xx+(2*radius_x/radius-scale_x/scale)*reduced_x)/scale**2
    lapse_laplacian = laplacian*observer_lapse
    normal_curvature = ((-scale_acceleration/scale-2*radius_acceleration/radius
                         +(lapse_rate/lapse)*(scale_rate/scale+2*radius_rate/radius))/lapse**2
                        +lapse_laplacian/lapse)
    kappa_energy = (curvature+2*theta**2/3-2*shear_squared)/2
    kappa_stress = 2*normal_curvature-kappa_energy
    volume = float(np.sum(measure)*dx)
    def moment(weight):
        return float(np.sum(weight*measure)*dx)
    row = clocked_scalar_average_sources(
        volume, lapse_theta_integral=moment(reduced_lapse*theta),
        lapse_squared_theta_squared_integral=moment(reduced_lapse**2*theta**2),
        lapse_squared_shear_squared_integral=moment(reduced_lapse**2*shear_squared),
        lapse_squared_curvature_integral=moment(reduced_lapse**2*curvature),
        lapse_squared_energy_integral=moment(reduced_lapse**2*kappa_energy/KAPPA),
        lapse_squared_stress_trace_integral=moment(reduced_lapse**2*kappa_stress/KAPPA),
        lapse_rate_theta_integral=moment(reduced_rate*theta),
        lapse_laplacian_integral=moment(reduced_lapse*laplacian))
    volume_rate = float(np.sum(4*math.pi*(scale_rate*radius**2+2*scale*radius*radius_rate))*dx)
    volume_acceleration = float(np.sum(4*math.pi*(
        scale_acceleration*radius**2+4*scale_rate*radius*radius_rate
        +2*scale*(radius_rate**2+radius*radius_acceleration)))*dx)
    direct_acceleration = (volume_acceleration/(volume*observer_lapse**2)
                           -volume_rate*observer_lapse_rate/(volume*observer_lapse**3)
                           -(2/3)*(volume_rate/(volume*observer_lapse))**2)
    curvature_integral = float(np.sum(curvature*measure)*dx)
    curvature_reference = float(8*math.pi*np.sum(scale+radius_x**2/scale)*dx)
    lapse_reference = float(-4*math.pi*np.sum(radius**2/scale*reduced_x**2)*dx)
    row.update({
        'volume': volume,
        'observer_lapse': observer_lapse,
        'direct_volume_acceleration': direct_acceleration,
        'curvature_proper_integral': curvature_integral,
        'curvature_ibp_residual': curvature_integral-curvature_reference,
        'lapse_ibp_residual': moment(reduced_lapse*laplacian)-lapse_reference,
        'scope': ('periodic S1 x S2 Einstein representation of one inherited '
                  'common-action snapshot; not the P001 parent and not a vacuum identification'),
    })
    return row


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
