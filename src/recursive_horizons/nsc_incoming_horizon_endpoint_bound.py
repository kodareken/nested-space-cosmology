"""Finite-collar vacuum formula versus the exact affine-horizon projector.

One group32 E=[24,32] enclosure. No field ODE or radial recurrence is solved.
The frozen source kappa is retained; the profile slope is a geometry quantity.
"""
import math

import mpmath as mp
import numpy as np
import sympy as sp

from .nsc_incoming_vacuum_tail_bound import _precision, _horizon_enclosure, _lo, _hi, _range
from .nsc_incoming_middle_bound import _parameters
from .nsc_incoming_source_quadrature_bound import pack, unpack, float_enclosure
from .nsc_unruh_state import ParentDirac
from .nsc_paired_horizon_preparation import PairedHorizonSeedMap


def frobenius_residual_identity():
    """Reuse the owned finite frame series; verify its two scalar cancellations."""
    j, nu = sp.symbols('j nu', real=True)
    # B_(j-1)/A_j=2*i*j and B_j/A_j=-i/[2*(j+1/2+i*nu)].
    first = sp.simplify(j+sp.I*(2*sp.I*j)/2)
    second = sp.simplify((j+sp.Rational(1, 2)+sp.I*nu)*(-sp.I/(2*(j+sp.Rational(1, 2)+sp.I*nu)))+sp.I/2)
    if first != 0 or second != 0:
        raise ArithmeticError('existing Frobenius coefficient recurrence failed')
    return {'coefficient_residuals': [str(first), str(second)],
            'phase_free_equations': 'a_y+i*z*b/2=R_a; b_y+i*z*a/2+i*nu*b=0',
            'only_finite_series_residual': 'R_a=i*z*b_last/2; terms0..9 imply R_a=O(delta^10)',
            'vacuum_phase_identity': '(exp(i*E*x)*v)(exp(i*E*x)*v)^dagger=v*v^dagger',
            'covariance_residual': 'R*v^dagger+v*R^dagger; norm <=2*norm(v)*norm(R)',
            'normalization_policy': 'finite column outer product retained; no normalization or clipping',
            'unitary_bound_owner': 'docs/nsc-incoming-vacuum-tail-bound.md'}


def _abs_upper(value):
    return mp.iv.mpf(_hi(abs(value)))


def _positive_terms(t, nu, shift, count=10):
    values = [mp.iv.mpf(1)]
    for j in range(1, count):
        values.append(values[-1]*t/(4*j*mp.iv.sqrt((j-1+shift)**2+nu*nu)))
    return values


def _report(value):
    return {'binary_interval': pack(value), 'float_enclosure': float_enclosure(value)}


def endpoint_bound(channel, config, *, precision=80, contextual_lapse_budget=2.10969e-11):
    """Directed uniform mathematical-initializer error; arithmetic stays OPEN.

    Compare the actual finite formula at delta_num with the exact-profile
    affine vacuum at q=q_h_num-delta_num. Source occupations are not touched.
    """
    if channel['index'] != 32 or channel['compact_level'] != 2 or channel['angular_level'] != 9:
        raise ValueError('this bounded endpoint pilot owns group32 only')
    if config['horizon_offset'] != 1e-10 or config['magnetic_flux'] != 4 or precision < 50:
        raise ValueError('unchanged inherited collar, magnetic sector and directed precision required')
    with _precision(precision):
        horizon, derivative, _ = _horizon_enclosure(config['horizon_rho'], config['surface_gravity'], precision)
        k_geometry = derivative/2
        q_exact = mp.iv.pi/2+mp.iv.atan2(horizon, mp.iv.mpf(1))
        r_exact = mp.iv.sqrt(1+horizon*horizon)
        background = ParentDirac(config['horizon_rho'], config['surface_gravity'])
        q_stored = mp.iv.mpf(float(background.horizon_q))
        r_stored = mp.iv.mpf(math.sqrt(1+config['horizon_rho']**2))
        k_source = mp.iv.mpf(config['surface_gravity'])
        delta_stored = mp.iv.mpf(math.exp(math.log(config['horizon_offset'])))
        q_shift = q_exact-q_stored
        delta_physical = delta_stored+q_shift
        delta_max = mp.iv.mpf(max(_hi(delta_stored), _hi(delta_physical)))
        delta_min = mp.iv.mpf(min(_lo(delta_stored), _lo(delta_physical)))
        if _lo(delta_min) <= 0:
            raise ArithmeticError('finite archived collar does not lie inside the exact horizon')

        # Geometry only: exact q profile W=3(pi-q)+1.5 sin(2q)-sin(q)^2.
        q_near = _range(_lo(q_exact)-_hi(delta_max), _hi(q_exact))
        W2 = -6*mp.iv.sin(2*q_near)-2*mp.iv.cos(2*q_near)
        W5 = 48*mp.iv.cos(2*q_near)-16*mp.iv.sin(2*q_near)
        kr = _abs_upper(mp.iv.cos(q_near)/(mp.iv.sin(q_near)**2))
        epsilon_coefficient = _abs_upper(W2)/(4*k_geometry)
        epsilon_max = epsilon_coefficient*delta_max
        if _hi(epsilon_max) >= 1:
            raise ArithmeticError('exact near-horizon W positivity not certified')
        # W/(2*k_geometry*delta)=1+epsilon, |epsilon|<=J*delta.
        inv_geometry = 1/mp.iv.sqrt(2*k_geometry)
        inv_source = 1/mp.iv.sqrt(2*k_source)
        inv_max = inv_geometry/mp.iv.sqrt(1-epsilon_max)
        inv_zero_difference = _abs_upper(inv_geometry-inv_source)
        inv_linear = inv_geometry*epsilon_coefficient/(2*(1-epsilon_max)**mp.iv.mpf('1.5'))
        diagonal_zero = _abs_upper(1/(2*k_geometry)-1/(2*k_source))
        diagonal_linear = epsilon_coefficient/(2*k_geometry*(1-epsilon_max))

        a_in, r_in, mass, ell, factor = _parameters(channel)
        m0, ell0 = mp.iv.mpf(channel['compact_mass']), mp.iv.mpf(channel['angular_eigenvalue'])
        mu0 = mp.iv.sqrt((m0*r_stored)**2+ell0*ell0)
        radius_zero = _abs_upper(r_exact-r_stored)
        radius_max = r_exact+kr*delta_max
        label_zero = _abs_upper(mass-m0)*radius_max+_abs_upper(ell-ell0)
        transverse_zero = mu0*inv_zero_difference+inv_max*(m0*radius_zero+label_zero)
        transverse_linear = mu0*inv_linear+inv_max*mass*kr

        nu_min = mp.iv.mpf(24)/k_source
        t_max = 2*mu0*mu0*delta_max/k_source
        terms_a = _positive_terms(t_max, nu_min, mp.iv.mpf('0.5'))
        terms_b = _positive_terms(t_max, nu_min, mp.iv.mpf('1.5'))
        A = sum(terms_a, mp.iv.mpf(0))
        B = mu0*mp.iv.sqrt(2/k_source)/mp.iv.sqrt(1+4*nu_min*nu_min)*sum(terms_b, mp.iv.mpf(0))
        column_norm_squared = A*A+B*B*delta_max
        sqrt_delta = mp.iv.sqrt(delta_max)
        integral_half = 2*sqrt_delta
        integral_three_halves = 2*delta_max*sqrt_delta/3
        transverse_error = 2*column_norm_squared*(transverse_zero*integral_half+transverse_linear*integral_three_halves)
        diagonal_error = 2*A*B*32*(diagonal_zero*integral_half+diagonal_linear*integral_three_halves)
        residual_endpoint = t_max/(2*mp.iv.sqrt(1+4*nu_min*nu_min))*terms_b[-1]
        truncation_error = 2*mp.iv.sqrt(column_norm_squared)*residual_endpoint/10
        aligned_error = transverse_error+diagonal_error+truncation_error

        # Coordinate bridge: P_frame(delta_physical)-P_frame(delta_stored).
        # Bound its derivative on the tiny segment, without a field solve.
        da = sum((j*term for j, term in enumerate(terms_a)), mp.iv.mpf(0))/delta_max
        db = mu0*mp.iv.sqrt(2/k_source)/mp.iv.sqrt(1+4*nu_min*nu_min)*sum(
            ((j+mp.iv.mpf('0.5'))*term for j, term in enumerate(terms_b)), mp.iv.mpf(0))/mp.iv.sqrt(delta_min)
        projector_derivative = 2*mp.iv.sqrt(column_norm_squared)*mp.iv.sqrt(da*da+db*db)
        coordinate_error = projector_derivative*_abs_upper(q_shift)
        total_error = aligned_error+coordinate_error

        H_upper = mp.iv.sqrt(mass*mass+(ell/r_in)**2+(mp.iv.mpf(32)/a_in)**2)
        # factor already includes the full actual angular pair and negative E.
        density_conversion = 2*factor*8*H_upper
        lapse_conversion = 4*mp.iv.pi*a_in*r_in*r_in*density_conversion
        lapse = lapse_conversion*total_error
        aligned_lapse = lapse_conversion*aligned_error
        coordinate_lapse = lapse_conversion*coordinate_error
        # The old small-delta W routine retains degree4. Its omitted fifth
        # derivative is bounded on the same exact-profile interval. We do not
        # substitute that polynomial into the proof's generator.
        W_tail = _abs_upper(W5)*delta_max**5/mp.factorial(5)
        W_tail_relative = W_tail/(2*k_geometry*delta_max*(1-epsilon_max))
        values = {'exact_horizon_rho': horizon, 'exact_horizon_q': q_exact,
                  'stored_horizon_q': q_stored, 'q_coordinate_shift': q_shift,
                  'exact_geometric_kappa': k_geometry, 'frozen_source_kappa': k_source,
                  'geometric_minus_source_kappa': k_geometry-k_source,
                  'stored_delta': delta_stored, 'true_distance_at_stored_endpoint': delta_physical,
                  'relative_W_correction_upper': epsilon_max,
                  'finite_series_norm_defect_upper': truncation_error,
                  'finite_series_covariance_residual_integral_upper': truncation_error,
                  'geometry_transverse_residual_integral_upper': transverse_error,
                  'geometry_diagonal_residual_integral_upper': diagonal_error,
                  'aligned_coordinate_projector_error_upper': aligned_error,
                  'coordinate_bridge_projector_error_upper': coordinate_error,
                  'total_mathematical_initializer_error_upper': total_error,
                  'aligned_coordinate_lapse_error_upper': aligned_lapse,
                  'coordinate_bridge_lapse_error_upper': coordinate_lapse,
                  'total_mathematical_initializer_lapse_error_upper': lapse,
                  'exact_profile_degree4_W_remainder_upper': W_tail,
                  'exact_profile_degree4_W_relative_remainder_upper': W_tail_relative}
        return {'group': 32, 'energy_interval': [24, 32], 'actual_angular_signs': [1, -1],
                'Frobenius_terms': 10, 'precision': precision,
                'bounds': {key: _report(value) for key, value in values.items()},
                'contextual_remaining_lapse_budget': contextual_lapse_budget,
                'mathematical_endpoint_budget_pass': _hi(lapse) <= mp.mpf(contextual_lapse_budget),
                'aligned_coordinate_budget_pass': _hi(aligned_lapse) <= mp.mpf(contextual_lapse_budget),
                'input_interpretation': {'source_kappa': 'unchanged binary input in frame and physical state; never replaced by profile derivative',
                                         'geometric_kappa': 'derivative of the existing exact profile at its enclosed root; no refit',
                                         'point_comparison': 'q=q_h_stored-delta_stored, hence true delta=delta_stored+q_h_exact-q_h_stored',
                                         'frame': 'mathematical ten-term working_frame formula with inherited scalar arguments; floating evaluation error is separate',
                                         'mass_angular': 'exact defining labels and stored binary labels enclosed together'},
                'scope': {'endpoint_formula_bias_enclosed': True, 'both_signs_included_once': True,
                          'archived_floating_initializer_error_bound': None,
                          'archived_DOP853_propagation_error_bound': None,
                          'physical_thermal_error_included': False,
                          'radial_or_mode_solves': 0, 'metric_timestep': False,
                          'source_covariance_modified': False, 'finite_frame_normalized': False,
                          'actual_source_quadrature_certified': False}}


def phase_free_column(delta, energy, mass, angular, horizon_radius, kappa):
    """High-precision evaluation control of the actual finite frame formula."""
    d, E, m, ell, rh, k = map(mp.mpf, (delta, energy, mass, angular, horizon_radius, kappa))
    nu = E/k; mu = mp.sqrt((m*rh)**2+ell*ell); z = mu*mp.sqrt(2*d/k)
    def series(parameter):
        value, term = mp.mpc(1), mp.mpc(1)
        for j in range(1, 10):
            term *= -z*z/(4*j*(parameter+j-1)); value += term
        return value
    vector = mp.matrix([series(mp.mpf('0.5')+1j*nu),
                        -1j*z/(1+2j*nu)*series(mp.mpf('1.5')+1j*nu)])
    angle = mp.pi/4+mp.atan2(m*rh, ell)/2
    return mp.matrix([mp.exp(-1j*angle)*vector[0], mp.exp(1j*angle)*vector[1]])


def archived_formula_controls(channel, config):
    """Numerical controls only, never a rigorous floating/ODE error claim."""
    preparation = PairedHorizonSeedMap(*(config[key] for key in
        ('horizon_rho', 'surface_gravity', 'omega', 'horizon_offset', 'scattering_tolerance')))
    y = math.log(config['horizon_offset']); delta = math.exp(y)
    difference = norm = 0.
    with mp.workdps(80):
        for sign in (1, -1):
            for energy in (24., 28., 32.):
                actual = preparation.working_frame(y, energy, channel['compact_mass'], sign*channel['angular_eigenvalue'])[:, 0]
                controlled = phase_free_column(delta, energy, channel['compact_mass'], sign*channel['angular_eigenvalue'],
                                               preparation.horizon_radius, config['surface_gravity'])
                v = np.array([complex(x) for x in controlled])
                difference = max(difference, float(np.max(abs(np.outer(actual, actual.conj())-np.outer(v, v.conj())))))
                norm = max(norm, float(abs(np.vdot(actual, actual)-1)))
    return {'sampled_working_frame_covariance_difference': difference,
            'sampled_raw_frame_norm_residual': norm, 'samples': 6,
            'interpretation': 'arithmetic/convention indicators only; not uniform error bounds'}
