"""Geometry-matched affine vacuum initializer, before any field propagation.

This constructs a NEW numerical approximation from a two-term Frobenius
ratio. No archived column is normalized, edited or promoted to an exact mode.
"""
import mpmath as mp
import numpy as np
import sympy as sp

from .nsc_incoming_vacuum_tail_bound import _precision, _horizon_enclosure, _lo, _hi, _range
from .nsc_incoming_middle_bound import _parameters
from .nsc_incoming_source_quadrature_bound import pack, unpack, float_enclosure


def matched_ratio_identity():
    """Exact first two coefficients of the owned canonical Dirac Riccati law."""
    D0, D1 = sp.symbols('D0 D1', real=True)
    K0, Kb0, K1 = sp.symbols('K0 Kbar0 K1')
    w0 = -sp.I*K0/(sp.Rational(1, 2)+2*sp.I*D0)
    w1 = (-sp.I*K1-2*sp.I*D1*w0+sp.I*Kb0*w0**2)/(sp.Rational(3, 2)+2*sp.I*D0)
    residuals = [sp.simplify((sp.Rational(1, 2)+2*sp.I*D0)*w0+sp.I*K0),
                 sp.simplify((sp.Rational(3, 2)+2*sp.I*D0)*w1+sp.I*K1+2*sp.I*D1*w0-sp.I*Kb0*w0**2)]
    if residuals != [0, 0]:
        raise ArithmeticError('geometry-matched Frobenius cancellation failed')
    return {'coefficient_residuals': [str(v) for v in residuals],
            'ratio': 'w=sqrt(delta)*(w0+w1*delta)',
            'Riccati_equation': 'delta*w_delta=-i*sqrt(delta)*K-2i*D*w+i*sqrt(delta)*conj(K)*w^2',
            'residual_order': 'delta^(5/2)',
            'projector_defect_identity': 'norm(P_y+i[G,P])=abs(Riccati_residual)/(1+abs(w)^2)',
            'projector_error_upper': '(2/5)*residual_coefficient_upper*delta0^(5/2)',
            'vacuum_limit': 'P(delta)->diag(1,0)',
            'source_kappa_policy': 'frozen source occupation kappa is not the geometric indicial exponent'}


def _upper(value):
    return mp.iv.mpf(_hi(abs(value)))


def _describe(value):
    return {'binary_interval': pack(value), 'float_enclosure': float_enclosure(value)}


def _complex_describe(value):
    return {'real': _describe(value.real), 'imag': _describe(value.imag)}


def _geometry(channel, config):
    if channel['index'] != 32 or channel['compact_level'] != 2 or channel['angular_level'] != 9:
        raise ValueError('matched initializer pilot owns group32 only')
    if config['horizon_offset'] != 1e-10 or config['magnetic_flux'] != 4:
        raise ValueError('fixed inherited collar and magnetic sector required')
    h, derivative, _ = _horizon_enclosure(config['horizon_rho'], config['surface_gravity'], mp.iv.dps)
    kg = derivative/2
    qh = mp.iv.pi/2+mp.iv.atan2(h, mp.iv.mpf(1))
    rh = mp.iv.sqrt(1+h*h)
    d = mp.iv.mpf(config['horizon_offset'])
    delta = _range(mp.mpf(0), _hi(d))
    q = _range(_lo(qh)-_hi(d), _hi(qh))
    W2 = -6*mp.iv.sin(2*q)-2*mp.iv.cos(2*q)
    W3 = -12*mp.iv.cos(2*q)+4*mp.iv.sin(2*q)
    b0 = 2*kg
    b1 = (-6*mp.iv.sin(2*qh)-2*mp.iv.cos(2*qh))/2
    b2 = -W3/6
    B = b0+b1*delta+b2*delta*delta
    if _lo(B) <= 0:
        raise ArithmeticError('exact-profile W/delta is not positively enclosed')
    # Integral Taylor formulas of B=W/delta give B'=int t W'' dt
    # and B''=-int t² W''' dt, including delta=0 without division by0.
    derivative_B = W2/2
    second_B = -W3/3
    f0 = 1/mp.iv.sqrt(b0)
    f1 = -b1/(2*b0**mp.iv.mpf('1.5'))
    f2 = 3*derivative_B*derivative_B/(8*B**mp.iv.mpf('2.5'))-second_B/(4*B**mp.iv.mpf('1.5'))
    r1 = mp.iv.cos(qh)/(mp.iv.sin(qh)**2)
    r2 = (1/mp.iv.sin(q)+2*mp.iv.cos(q)**2/mp.iv.sin(q)**3)/2
    a, r_in, mass, angular, factor = _parameters(channel)
    return {'horizon': h, 'qh': qh, 'rh': rh, 'kg': kg, 'delta0': d, 'delta': delta,
            'b0': b0, 'b1': b1, 'b2': b2, 'B': B,
            'f0': f0, 'f1': f1, 'f2': f2, 'r1': r1, 'r2': r2,
            'a': a, 'r_in': r_in, 'mass': mass, 'angular': angular, 'factor': factor}


def _coefficients(g, energy, sign):
    if sign not in (1, -1):
        raise ValueError('actual signed angular channel required')
    E, m, ell = energy, g['mass'], sign*g['angular']
    delta, b0, b1, b2, B = (g[k] for k in ('delta', 'b0', 'b1', 'b2', 'B'))
    N0 = -m*g['rh']+mp.iv.j*ell
    N1, N2 = -m*g['r1'], -m*g['r2']
    K0, K1 = N0*g['f0'], N1*g['f0']+N0*g['f1']
    K2 = N2*g['f0']+N1*g['f1']+N0*g['f2']+delta*(N2*g['f1']+N1*g['f2'])+delta*delta*N2*g['f2']
    D0, D1 = E/b0, -E*b1/(b0*b0)
    D2 = E*(b1*b1-b0*b2+b1*b2*delta)/(b0*b0*B)
    w0 = -mp.iv.j*K0/(mp.iv.mpf('0.5')+2*mp.iv.j*D0)
    Kbar0 = mp.iv.mpc(K0.real, -K0.imag)
    w1 = (-mp.iv.j*K1-2*mp.iv.j*D1*w0+mp.iv.j*Kbar0*w0*w0)/(mp.iv.mpf('1.5')+2*mp.iv.j*D0)
    return {'K0': K0, 'K1': K1, 'K2': K2, 'D0': D0, 'D1': D1, 'D2': D2, 'w0': w0, 'w1': w1}


def _residual_coefficient(g, c):
    d = g['delta0']
    K0, K1, K2, D1, D2, w0, w1 = (_upper(c[key]) for key in ('K0', 'K1', 'K2', 'D1', 'D2', 'w0', 'w1'))
    q = w0+d*w1
    return K2+2*(D1*w1+D2*q)+K1*q*q+K0*(2*w0*w1+d*w1*w1)+d*K2*q*q


def matched_initializer_bound(channel, config, *, precision=80, contextual_lapse_budget=2.10969e-11):
    """Uniform mathematical endpoint certificate on both E24–32 channels."""
    if precision < 50:
        raise ValueError('directed precision>=50 required')
    with _precision(precision):
        g = _geometry(channel, config)
        E = _range(mp.mpf(24), mp.mpf(32))
        H = mp.iv.sqrt(g['mass']**2+(g['angular']/g['r_in'])**2+(32/g['a'])**2)
        # One half of the group factor for each actual angular sign.
        conversion_per_sign = 4*mp.iv.pi*g['a']*g['r_in']**2*2*(g['factor']/2)*8*H
        rows, total = [], mp.iv.mpf(0)
        for sign in (1, -1):
            c = _coefficients(g, E, sign)
            R = _residual_coefficient(g, c)
            distance = 2*R*g['delta0']**mp.iv.mpf('2.5')/5
            lapse = conversion_per_sign*distance
            total += lapse
            rows.append({'sign': sign, 'coefficients': {key: _complex_describe(value) for key, value in c.items()},
                         'residual_coefficient_upper': _describe(R), 'projector_error_upper': _describe(distance),
                         'lapse_error_upper': _describe(lapse)})
        return {'group': 32, 'energy_interval': [24, 32], 'precision': precision,
                'source_order': 'new normalized rankone approximation from sqrt(delta)*(w0+w1 delta)',
                'geometry': {key: _describe(g[key]) for key in ('horizon', 'qh', 'rh', 'kg', 'delta0', 'B', 'b2', 'f2', 'r2')},
                'frozen_source_kappa': config['surface_gravity'], 'angular_signs': [1, -1], 'per_sign': rows,
                'total_lapse_error_upper': _describe(total),
                'contextual_remaining_lapse_budget': contextual_lapse_budget,
                'mathematical_endpoint_budget_pass': _hi(total) <= mp.mpf(contextual_lapse_budget),
                'source_covariance_policy': 'vacuum affine limiting projector unchanged; thermal/coherent source covariance never evaluated or changed',
                'scope': {'new_numerical_initializer': True, 'archived_columns_normalized_or_changed': False,
                          'uniform_mathematical_endpoint_error': True,
                          'uniform_binary_representation_error': None,
                          'actual_ODE_or_spectral_quadrature_error': None,
                          'radial_or_mode_solves': 0, 'physical_IV_selected': False, 'metric_timestep': False}}


def represent_initializer(channel, config, energy, sign, *, precision=80):
    """Directed point-energy ratio and a binary representation with its bound.

    The returned binary ratio defines a new normalized mathematical rankone
    projector. No floating matrix normalization is assumed exact.
    """
    if not 24 <= energy <= 32 or sign not in (1, -1):
        raise ValueError('point within E24–32 and actual angular sign required')
    with _precision(precision):
        g = _geometry(channel, config)
        c = _coefficients(g, mp.iv.mpf(energy), sign)
        w = mp.iv.sqrt(g['delta0'])*(c['w0']+c['w1']*g['delta0'])
        real = float((_lo(w.real)+_hi(w.real))/2)
        imag = float((_lo(w.imag)+_hi(w.imag))/2)
        distance = _upper(w-mp.iv.mpc(real, imag))
        # Exact normalized-projector distance =|w-wb|/sqrt((1+|w|²)(1+|wb|²)) <=|w-wb|.
        return {'energy': energy, 'sign': sign, 'ratio_interval': _complex_describe(w),
                'binary_ratio': [real, imag], 'projector_representation_error_upper': _describe(distance),
                'interpretation': 'point-energy representation only; normalization defines the new mathematical projector, not an alteration of archived modes'}
