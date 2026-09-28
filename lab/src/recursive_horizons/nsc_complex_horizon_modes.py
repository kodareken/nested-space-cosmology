"""Holomorphic complex-frequency extension of the published real PG/Jost modes.

The real horizon frame, massive outgoing spinor series, exterior Dirac
Riccati, matched_delta_q interior collar, inverse-r_h exterior collar, and
PG clock/frame conversions are reused. Frequency z is not conjugated inside
the second Frobenius column. No seed covariance, R sewing, metric evolution,
or packet integral is performed here.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import atan2, exp, log, pi, sqrt

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_paired_horizon_preparation import PairedHorizonSeedMap
from .nsc_pg_massive_modes import MassivePGModeResolution


def _frobenius_series(sign, argument_squared, parameter, terms):
    value = 1. + 0j
    term = 1. + 0j
    for j in range(1, terms):
        denominator = parameter + j - 1
        if abs(denominator) < 1e-14:
            raise ValueError('Frobenius series encountered a pole inside the local band')
        term *= sign * argument_squared / (4 * j * denominator)
        value += term
    return value


def _require_channel(kappa, effective, distance):
    if not np.isfinite(kappa) or kappa <= 0 or effective < 0 or distance < 0:
        raise ValueError('positive horizon scale and nonnegative effective channel required')


def _require_frequency(z, kappa, *, upper_half=False):
    z = complex(z)
    kappa = float(kappa)
    if not np.isfinite([z.real, z.imag, kappa]).all() or kappa <= 0:
        raise ValueError('finite frequency and positive surface gravity required')
    if abs(z.imag) >= kappa / 2:
        raise ValueError(
            '|Im z| < kappa/2 Frobenius pole bound required; contour height is not a history duration'
        )
    if z.real < 0:
        raise ValueError('nonnegative real-frequency chart required')
    if upper_half and z.imag < 0:
        raise ValueError('outgoing root is defined for Re z>0, Im z>=0')
    if upper_half and z.real == 0:
        raise ValueError('outgoing massive root requires Re z>0')
    return z, kappa


def analytic_horizon_frame(z, kappa, effective, distance, tortoise, timelike=False, terms=10):
    """Holomorphic Dirac Frobenius frame with conjugate-analytic second column.

    Plus/minus exponentials and series use ``.5 +/- i*z/kappa`` of the same
    z. Second-column signs reduce to ``horizon_frame`` on the real axis.
    """
    z, kappa = _require_frequency(z, kappa)
    _require_channel(kappa, effective, distance)
    if terms < 2:
        raise ValueError('resolved Frobenius truncation required')
    coupling = float(effective) * float(distance)
    sign = -1 if timelike else 1
    nu = z / kappa
    phase_plus = np.exp(1j * z * tortoise)
    phase_minus = np.exp(-1j * z * tortoise)
    a_plus = phase_plus * _frobenius_series(sign, coupling * coupling, .5 + 1j * nu, terms)
    a_minus = phase_minus * _frobenius_series(sign, coupling * coupling, .5 - 1j * nu, terms)
    i_plus = -1j if timelike else 1j
    i_minus = 1j if timelike else -1j
    b_plus = i_plus * coupling / (1 + 2j * nu) * phase_plus * _frobenius_series(
        sign, coupling * coupling, 1.5 + 1j * nu, terms)
    b_minus = i_minus * coupling / (1 - 2j * nu) * phase_minus * _frobenius_series(
        sign, coupling * coupling, 1.5 - 1j * nu, terms)
    top_minus = -b_minus if timelike else b_minus
    return np.array([[a_plus, top_minus], [b_plus, a_minus]], complex)


def outgoing_momentum(z, mass):
    """Outgoing root p=sqrt(z^2-m^2) on Re z>0, Im z>=0."""
    z = complex(z)
    mass = float(mass)
    if z.real <= 0 or z.imag < 0 or mass < 0:
        raise ValueError('outgoing root requires Re z>0, Im z>=0 and nonnegative mass')
    return np.sqrt(z * z - mass * mass)


def complex_outgoing_ratio(z, mass, angular, radius, order=8):
    """Holomorphic massive outgoing spinor series; same recurrence as outgoing_ratio."""
    z = complex(z)
    mass = float(mass)
    angular = float(angular)
    momentum = outgoing_momentum(z, mass)
    if radius <= 0 or order < 1:
        raise ValueError('positive radius and resolved series order required')
    if abs(momentum) < 1e-14:
        raise ValueError('outgoing root vanishes at the massive threshold')
    metric = np.zeros(order + 2)
    metric[0] = 1.
    for j in range((order + 2) // 2):
        power = 2 * j + 1
        if power < len(metric):
            metric[power] = -6 * (-1) ** j / ((2 * j + 1) * (2 * j + 3))
    root = np.zeros_like(metric)
    root[0] = 1.
    for n in range(1, len(root)):
        root[n] = (metric[n] - sum(root[j] * root[n - j] for j in range(1, n))) / 2
    sphere = np.zeros_like(metric)
    sphere[0] = 1.
    coefficient = 1.
    for n in range(1, len(sphere) // 2):
        coefficient *= (-.5 - (n - 1)) / n
        sphere[2 * n] = coefficient
    first = np.zeros_like(metric)
    first[1:] = angular * np.convolve(root, sphere)[:len(metric) - 1]
    second = mass * root
    minus = first - 1j * second
    plus = first + 1j * second
    ratio = np.zeros(order + 1, complex)
    ratio[0] = 1j * mass / (z + momentum) if mass else 0.
    for n in range(1, order + 1):
        derivative = -sum(metric[n - j - 1] * j * ratio[j] for j in range(1, n))
        square = np.convolve(ratio, ratio)
        forcing = 1j * (plus[n] + np.convolve(minus, square)[n])
        ratio[n] = (forcing - derivative) / (2j * momentum)
    powers = radius ** (-np.arange(order + 1, dtype=float))
    value = ratio @ powers
    derivative = -(np.arange(order + 1) * ratio @ powers) / radius
    return complex(value), complex(derivative), ratio, momentum


def _oriented_collar_frame(z, mass, angular, preparation, *, timelike):
    p = preparation
    rh = p.horizon_radius
    effective = float(np.hypot(angular, mass * rh))
    phase = atan2(mass * rh, angular)
    kappa = p.surface_gravity
    collar = p.collar_delta_q
    if timelike:
        signed_offset = -rh * rh * collar + p.horizon_rho * rh * rh * collar ** 2
        distance = sqrt(2 * collar / kappa)
        rotation = np.diag(np.exp(np.array([-1j, 1j]) * (pi / 4 + phase / 2)))
        tortoise = p.background.near_tortoise(signed_offset)
        return rotation @ analytic_horizon_frame(
            z, kappa, effective, distance, tortoise, timelike=True)
    distance = sqrt(2 * collar / kappa) / rh
    rotation = np.diag(np.exp(np.array([-1j, 1j]) * phase / 2))
    tortoise = p.background.near_tortoise(collar)
    return rotation @ analytic_horizon_frame(
        z, kappa, effective, distance, tortoise, timelike=False)


@dataclass(frozen=True)
class ComplexHorizonReflection:
    frequency: complex
    mass: float
    angular: float
    momentum: complex
    reflection: complex
    horizon_coefficients: np.ndarray
    outer_ratio: complex
    outer_radius: float
    residuals: dict


def complex_reflection(z, mass, angular, background, *, radial_collar=1e-10,
                       outer_radius=None, order=8, rtol=2e-13, atol=2e-15):
    """Exterior Jost reflection R(z) with the owned massive Riccati/series.

    ``background`` is the existing ParentDirac chart. The returned reflection
    is the horizon-coefficient ratio, not a current-renormalized magnitude.
    """
    mass = float(mass)
    angular = float(angular)
    z, kappa = _require_frequency(z, background.surface_gravity, upper_half=True)
    if mass < 0 or radial_collar <= 0:
        raise ValueError('nonnegative mass and positive radial collar required')
    if abs(z - mass) < 1e-12 or abs(z + mass) < 1e-12:
        raise ValueError('positive frequency away from massive threshold required')
    momentum = outgoing_momentum(z, mass)
    end = (max(60., 30 * max(1., abs(angular), mass) / max(abs(momentum), .2))
           if outer_radius is None else float(outer_radius))
    for outer_expansions in range(16):
        ratio, derivative, _, momentum = complex_outgoing_ratio(z, mass, angular, end, order)
        metric_a = background.A_from_offset(end - background.horizon_rho)
        v1 = angular * sqrt(metric_a) / sqrt(1 + end * end)
        v2 = mass * sqrt(metric_a)
        asymptotic_residual = abs(
            metric_a * derivative - 1j * (v1 + 1j * v2) + 2j * z * ratio
            - 1j * (v1 - 1j * v2) * ratio * ratio)
        initial_current = float(1 - abs(ratio) ** 2)
        if outer_radius is not None or asymptotic_residual < 3e-13:
            break
        end *= 2
    else:
        raise ArithmeticError('outer mode expansion failed its algebraic accuracy gate')
    if z.imag == 0 and z.real > mass and initial_current <= 0:
        raise ArithmeticError('outer outgoing current is unresolved')
    if z.imag == 0 and z.real < mass and abs(initial_current) > 1e-8:
        raise ArithmeticError('decaying outer expansion is not current-null to its declared accuracy')

    def rhs(y, state):
        offset = exp(y)
        rho = background.horizon_rho + offset
        a = background.A_from_offset(offset)
        first = angular * sqrt(a) / sqrt(1 + rho * rho)
        second = mass * sqrt(a)
        r = state[0]
        return [offset / a * (1j * (first + 1j * second) - 2j * z * r + 1j * (first - 1j * second) * r * r),
                1j * offset / a * (z - (first - 1j * second) * r)]

    start, stop = log(end - background.horizon_rho), log(radial_collar)
    if z.imag == 0 and z.real < mass:
        def phase_rhs(y, state):
            theta = float(state[0].real)
            r = np.exp(1j * theta)
            dr, dlog = rhs(y, np.array([r, state[1]], complex))
            return [float((dr / (1j * r)).real), dlog]
        run = solve_ivp(phase_rhs, (start, stop), np.array([np.angle(ratio), 0j]),
                        method='DOP853', rtol=rtol, atol=atol, max_step=.15)
        if not run.success:
            raise ArithmeticError(run.message)
        run.y[0] = np.exp(1j * run.y[0].real)
    else:
        run = solve_ivp(rhs, (start, stop), np.array([ratio, 0j]),
                        method='DOP853', rtol=rtol, atol=atol, max_step=.15)
        if not run.success:
            raise ArithmeticError(run.message)
    rh = sqrt(1 + background.horizon_rho ** 2)
    phase = atan2(mass * rh, angular)
    orient = np.diag(np.exp(np.array([-1j, 1j]) * phase / 2))
    frame = orient @ analytic_horizon_frame(
        z, kappa, np.hypot(angular, mass * rh),
        sqrt(2 * radial_collar / kappa) / rh, background.near_tortoise(radial_collar))
    coefficients = np.linalg.solve(frame, np.array([1., run.y[0, -1]], complex))
    reflection = coefficients[1] / coefficients[0]
    return ComplexHorizonReflection(
        z, mass, angular, momentum, complex(reflection), coefficients, complex(ratio), float(end),
        {
            'outer_Riccati_residual': float(asymptotic_residual),
            'outer_expansions': outer_expansions,
            'outer_radius': float(end),
            'closed_outer_current': abs(initial_current) if (z.imag == 0 and z.real < mass) else 0.,
            'function_evaluations': int(run.nfev),
        })


def wirtinger_residual(func, z, step, *, order=2):
    """Finite-difference Cauchy-Riemann residual of a holomorphic map."""
    z = complex(z)
    step = float(step)
    if step <= 0 or order not in (2, 4):
        raise ValueError('positive finite-difference step and order 2 or 4 required')

    def eval_at(point):
        return np.asarray(func(point), complex)

    if order == 2:
        dx = (eval_at(z + step) - eval_at(z - step)) / (2 * step)
        dy = (eval_at(z + 1j * step) - eval_at(z - 1j * step)) / (2 * step)
    else:
        dx = (-eval_at(z + 2 * step) + 8 * eval_at(z + step)
              - 8 * eval_at(z - step) + eval_at(z - 2 * step)) / (12 * step)
        dy = (-eval_at(z + 2j * step) + 8 * eval_at(z + 1j * step)
              - 8 * eval_at(z - 1j * step) + eval_at(z - 2j * step)) / (12 * step)
    d_bar = .5 * (dx + 1j * dy)
    d_z = .5 * (dx - 1j * dy)
    holomorphic = float(np.linalg.norm(d_z))
    residual = float(np.linalg.norm(d_bar))
    return {
        'cauchy_riemann': residual,
        'holomorphic_derivative': holomorphic,
        'relative': residual / max(1., holomorphic),
        'step': step,
        'order': order,
    }


def _matrix_ode_residual(sol, y, rhs_matrix, step=2e-4):
    def u(point):
        return sol(point).reshape(2, 2)
    derivative = (-u(y + 2 * step) + 8 * u(y + step) - 8 * u(y - step) + u(y - 2 * step)) / (12 * step)
    return float(np.linalg.norm(derivative - rhs_matrix @ u(y)))


def _transport_horizon_columns(z, mass, angular, resolution, interior_rho, exterior_rho):
    p = resolution.preparation
    inner_y = log(p.background.horizon_q - (pi / 2 + np.arctan(interior_rho)))
    begin = log(p.collar_delta_q)
    if inner_y < begin:
        raise ValueError('interior sample must lie beyond the declared matched numerical collar')
    outer_y = log(exterior_rho - p.horizon_rho)
    if outer_y < log(p.collar_delta_q):
        raise ValueError('exterior sample must lie beyond the declared radial collar')
    step = 2e-4
    frame_h = _oriented_collar_frame(z, mass, angular, p, timelike=True)
    run_i = solve_ivp(
        lambda y, u: (-1j * p.generator(y, z, mass, angular) @ u.reshape(2, 2)).ravel(),
        (begin, inner_y), frame_h.ravel(), method='DOP853',
        rtol=resolution.rtol, atol=resolution.atol, max_step=.1, dense_output=True)
    ext_h = _oriented_collar_frame(z, mass, angular, p, timelike=False)
    run_e = solve_ivp(
        lambda y, u: (resolution.exterior_generator(y, z, mass, angular) @ u.reshape(2, 2)).ravel(),
        (log(p.collar_delta_q), outer_y), ext_h.ravel(), method='DOP853',
        rtol=resolution.rtol, atol=resolution.atol, max_step=.1, dense_output=True)
    if not run_i.success or not run_e.success:
        raise ArithmeticError('complex-frequency horizon-column integration failed')
    inner_ks = run_i.y[:, -1].reshape(2, 2)
    outer_ks = run_e.y[:, -1].reshape(2, 2)
    inner_frame, _ = resolution.frame(interior_rho)
    outer_frame, _ = resolution.frame(exterior_rho)
    clock_i = float(resolution.clock(np.array(interior_rho)))
    clock_e = float(resolution.clock(np.array(exterior_rho)))
    interior = np.exp(1j * z * clock_i) * inner_frame @ inner_ks
    exterior = np.exp(1j * z * clock_e) * outer_frame @ outer_ks
    interior_ode = _matrix_ode_residual(
        run_i.sol, inner_y - 2*step, -1j * p.generator(inner_y - 2*step, z, mass, angular), step)
    exterior_ode = _matrix_ode_residual(
        run_e.sol, outer_y - 2*step, resolution.exterior_generator(outer_y - 2*step, z, mass, angular), step)
    frame_checks = {
        **{f'interior_{k}': v for k, v in resolution.frame_residual(interior_rho, z, mass, angular).items()},
        **{f'exterior_{k}': v for k, v in resolution.frame_residual(exterior_rho, z, mass, angular).items()},
    }
    return interior, exterior, {
        'interior_ode': interior_ode,
        'interior_ode_log_distance': inner_y - 2*step,
        'exterior_ode': exterior_ode,
        'exterior_ode_log_distance': outer_y - 2*step,
        'pg_clock_interior': clock_i,
        'pg_clock_exterior': clock_e,
        'function_evaluations': int(run_i.nfev + run_e.nfev),
        **frame_checks,
    }


@dataclass(frozen=True)
class FinitePGHorizonBasis:
    frequency: complex
    conjugate_frequency: complex
    mass: float
    angular: float
    interior_rho: float
    exterior_rho: float
    interior: np.ndarray
    exterior: np.ndarray
    interior_conjugate: np.ndarray
    exterior_conjugate: np.ndarray
    residuals: dict


def finite_pg_horizon_basis(z, mass, angular, preparation, *, interior_rho=0., exterior_rho=3.,
                            rtol=2e-13, atol=2e-15):
    """Unsewn PG horizon columns at one interior and one exterior chart point.

    Columns are evaluated at z and at conjugate(z) by the holomorphic frame,
    not by conjugating the z-columns. R sewing and covariance are omitted.
    """
    if not isinstance(preparation, PairedHorizonSeedMap):
        raise TypeError('PairedHorizonSeedMap owns the matched collar and PG clock data')
    mass = float(mass)
    angular = float(angular)
    interior_rho = float(interior_rho)
    exterior_rho = float(exterior_rho)
    z, _ = _require_frequency(z, preparation.surface_gravity)
    if mass < 0:
        raise ValueError('nonnegative compact mass required')
    if not (interior_rho < preparation.horizon_rho < exterior_rho):
        raise ValueError('finite samples on both sides of the fixed horizon required')
    resolution = MassivePGModeResolution(preparation, rtol=rtol, atol=atol)
    interior, exterior, residuals_z = _transport_horizon_columns(
        z, mass, angular, resolution, interior_rho, exterior_rho)
    z_bar = z.conjugate()
    if z.imag == 0:
        interior_c, exterior_c, residuals_c = interior, exterior, residuals_z
    else:
        _require_frequency(z_bar, preparation.surface_gravity)
        interior_c, exterior_c, residuals_c = _transport_horizon_columns(
            z_bar, mass, angular, resolution, interior_rho, exterior_rho)
    residuals = {
        'frequency': residuals_z,
        'conjugate_frequency': residuals_c,
        'pg_clock_is_geometric': float(abs(residuals_z['pg_clock_interior']
                                           - residuals_c['pg_clock_interior'])
                                       + abs(residuals_z['pg_clock_exterior']
                                             - residuals_c['pg_clock_exterior'])),
        'naive_interior_conjugation': float(np.linalg.norm(interior_c - interior.conj())),
        'naive_exterior_conjugation': float(np.linalg.norm(exterior_c - exterior.conj())),
    }
    return FinitePGHorizonBasis(
        z, z_bar, mass, angular, interior_rho, exterior_rho,
        interior, exterior, interior_c, exterior_c, residuals)
