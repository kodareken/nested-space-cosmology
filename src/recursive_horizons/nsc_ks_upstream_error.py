"""Directed a-posteriori error for homogeneous KS continuation of stored columns.

The exact generator G(rho) is anti-Hermitian for real a,r, so the 2-norm
evolution operator U is unitary. A continuous interpolant Y of a numerical
continuation on [1, rho_up] therefore satisfies

    e(rho_up) = U(rho_up, 1) e(1) + int_1^{rho_up} U(rho_up, s) R(s) ds,
    R = Y' - G Y,

    ||e(rho_up)||_2 <= exp(K) ( ||e(1)||_2 + int ||R||_2 ds ),

with K the integral of the 2-logarithmic norm. Exact anti-Hermiticity gives
K = 0; an interval remainder of the Hermitian part is retained when it is
proved positive. Stored binary64 rho=1 columns may be treated as exact only
for this numerical-continuation sub-bound. Physical rho=1 source error is a
separate explicit input and is never replaced by zero, by a Pauli identity,
or by a nested-node identity. A refinement movement is an indicator, never
the certificate bound.
"""
from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from math import comb
from pathlib import Path

import numpy as np
from flint import acb, arb, ctx, fmpq
from scipy.integrate import DOP853
from scipy.linalg import expm

from . import nsc_ks_ball_geometry as G
from .nsc_ks_ball_geometry import RHO_MAX, RHO_MIN, background_series
from .nsc_ks_ball_operator import slab_ball
from .nsc_ks_ball_trajectory import complex_ball, exact_upper, power_to_bernstein
from .nsc_ks_reference_error import homogeneous_reference_norms
from .nsc_ks_residual_error import matter_error_bounds, nonnegative
from .nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from .nsc_retarded_radial_response import ks_generator
from .nsc_transmitting_dirac_domain import S1, S2, S3


SCHEMA_CONTINUATION = 'NSC-KS-UPSTREAM-ERROR-PILOT-v1'
SCHEMA_DECISION = 'NSC-KS-LOW-SUBGAP-DECISION-v1'
RHO_UP = Fraction(1159676904047903, 1125899906842624)
RHO_UP_BINARY64 = float(RHO_UP)
AXIAL_LOWER = Fraction(4, 5)
RADIUS_LOWER = Fraction(7, 5)
ELL0_GROUPS = (0, 13, 23)
CONTINUATION_SOLVER = dict(rtol=5e-14, atol=5e-16, max_step=0.00025)
NINE_COMPONENT_ORDER = (
    'field_space_time',
    'changed_history_UV_tail',
    'baseline_low_subgap',
    'upstream',
    'between_node',
    'arithmetic',
    'covered_regions',
    'energy_interpolation',
    'phase_value',
)
UPSTREAM_ALLOCATION = 2e-12
KNOWN_PARTIAL_N = 8.999152980170087e-12
REMAINING_N_AFTER_KNOWN = 2e-11 - KNOWN_PARTIAL_N
IDENTITIES_NOT_ERROR_BOUNDS = (
    'Pauli 4-component subgap covariance reconstruction',
    'negative-source complement I-conj(C_src)',
    'ell=0 Dirac r-independence of the radius-history response',
    'nested-node / nested-gradient algebraic identities',
    'quadrature-node or refinement movement',
)
DIRECTED_FORMULA = (
    '||e(rho_up)||_2 <= exp(K) ( ||e(1)||_2 + int_1^{rho_up} ||Y\'-G Y||_2 d rho ), '
    'K = int mu_2(G) d rho, mu_2(G)=0 for the exact anti-Hermitian KS generator'
)
ASSUMPTIONS = (
    'homogeneous ks_generator on the background chart; no history potential',
    'stored binary64 rho=1 columns exact only for the numerical-continuation sub-bound',
    'physical rho=1 source error is a separate explicit input',
    'DOP853 dense output is the interpolant whose defect is enclosed',
    'refinement movement is an indicator, never the certificate',
    'V5 covered action regions are not mode-column reconstruction bounds',
)


def _nonneg(value, name):
    return nonnegative(value, name)


def _arb_fraction(value):
    q = Fraction(value)
    return arb(fmpq(q.numerator, q.denominator))


def pack_upper(value):
    return exact_upper(value)


def restore_packed_upper(value):
    if not isinstance(value, dict) or set(value) != {'mantissa', 'exponent'}:
        raise ValueError('exact-upper mantissa/exponent required')
    return arb(int(value['mantissa'])) * arb(2) ** int(value['exponent'])


def packed_upper_fraction(value):
    if not isinstance(value, dict) or set(value) != {'mantissa', 'exponent'}:
        raise ValueError('exact-upper mantissa/exponent required')
    exponent = int(value['exponent'])
    mantissa = int(value['mantissa'])
    return Fraction(mantissa * (1 << max(exponent, 0)), 1 << max(-exponent, 0))


def rho_up_is_archived_binary64(value):
    return float(value) == RHO_UP_BINARY64


def propagate_homogeneous_continuation(initial_row_error, residual_integral,
                                       logarithmic_norm_integral=0):
    """Endpoint 2-norm majorant of Y'=G Y from an a-posteriori defect integral.

    K is the integral of mu_2(G). It is not inferred from a refinement
    difference. Missing residual integrals are rejected, not treated as zero.
    """
    e0 = _nonneg(initial_row_error, 'initial row error')
    residual = _nonneg(residual_integral, 'continuous residual integral')
    growth_log = _nonneg(logarithmic_norm_integral, 'logarithmic-norm integral')
    from math import exp as _exp
    # Fraction exp is not available; return exact rational growth via series
    # enclosure only when K=0, otherwise a dyadic exp upper.
    if growth_log == 0:
        growth = Fraction(1)
        endpoint = e0 + residual
    else:
        growth = None
        endpoint = None
    return {
        'initial_row_error': e0,
        'residual_integral': residual,
        'logarithmic_norm_integral': growth_log,
        'growth_upper': growth,
        'endpoint_row_error_upper': endpoint,
        'formula': DIRECTED_FORMULA,
        'refinement_used_as_bound': False,
        'physical_source_error_included': False,
    }


def continuation_growth_and_endpoint(initial_row_error, residual_integral,
                                     logarithmic_norm_integral=0, *, bits=80):
    """Same majorant with a directed exponential enclosure when K>0."""
    e0 = _nonneg(initial_row_error, 'initial row error')
    residual = _nonneg(residual_integral, 'continuous residual integral')
    k = _nonneg(logarithmic_norm_integral, 'logarithmic-norm integral')
    with ctx.workprec(bits):
        growth = _arb_fraction(k).exp()
        endpoint = growth * (_arb_fraction(e0) + _arb_fraction(residual))
        if not growth.is_finite() or not endpoint.is_finite():
            raise ArithmeticError('finite continuation growth enclosure required')
        return {
            'growth_upper': pack_upper(growth),
            'endpoint_row_error_upper': pack_upper(endpoint),
            'formula': DIRECTED_FORMULA,
            'refinement_used_as_bound': False,
            'physical_source_error_included': False,
        }


def combine_upstream_row_error(continuation_row_error, physical_source_row_error):
    """Transport-ready row error. A missing physical input keeps the sum OPEN."""
    continuation = None if continuation_row_error is None else _nonneg(
        continuation_row_error, 'continuation row error')
    if physical_source_row_error is None:
        return {
            'continuation_row_error_upper': continuation,
            'physical_source_row_error_upper': None,
            'combined_row_error_upper': None,
            'combined_status': (
                'OPEN: physical rho=1 source error remains an explicit missing input'),
        }
    physical = _nonneg(physical_source_row_error, 'physical source row error')
    combined = None if continuation is None else continuation + physical
    return {
        'continuation_row_error_upper': continuation,
        'physical_source_row_error_upper': physical,
        'combined_row_error_upper': combined,
        'combined_status': (
            'OPEN: combined numerical-plus-physical upstream row error'
            if combined is None else
            'conditional: continuation plus explicit physical source row error'),
    }


def nine_component_upstream_feed(continuation_n_beta, physical_n_beta=None):
    """How a continuation sub-bound enters the nine-component budget's upstream slot.

    This does not write gate-budget v4. Physical source error remains explicit.
    """
    combined = None
    if continuation_n_beta is not None and physical_n_beta is not None:
        combined = [continuation_n_beta[0] + physical_n_beta[0],
                    continuation_n_beta[1] + physical_n_beta[1]]
    fits_allocation = None
    if continuation_n_beta is not None:
        fits_allocation = bool(continuation_n_beta[0] <= UPSTREAM_ALLOCATION
                               and continuation_n_beta[1] <= UPSTREAM_ALLOCATION)
    return {
        'nine_component_order': list(NINE_COMPONENT_ORDER),
        'budget_slot': 'upstream',
        'allocation_per_component': UPSTREAM_ALLOCATION,
        'continuation_sub_bound_N_beta': continuation_n_beta,
        'physical_source_error_bound_N_beta': physical_n_beta,
        'combined_upstream_bound_N_beta': combined,
        'continuation_fits_upstream_allocation': fits_allocation,
        'known_partial_error_sum_N': KNOWN_PARTIAL_N,
        'remaining_N_after_known_components': REMAINING_N_AFTER_KNOWN,
        'gate_budget_v4_edited': False,
        'combined_status': (
            'OPEN: physical rho=1 source error remains an explicit missing input'
            if physical_n_beta is None else
            'conditional combined upstream N,beta'),
    }


def independent_constant_matrix_control():
    """A-posteriori majorant versus an independently solved constant 2x2 system.

    G = i diag(1,-1) on [0, h], constant interpolant Y_h = Y(0). The residual
    is -G Y0. The exact endpoint error uses scipy.linalg.expm, not this
    formula and not a refinement difference.
    """
    h = Fraction(3, 100)
    y0 = np.array([1.0, 1.0], complex)
    g = 1j * np.diag([1.0, -1.0])
    residual_frobenius = float(np.linalg.norm(g @ y0))
    residual_integral = Fraction(residual_frobenius).limit_denominator(10**12) * h
    # The stored binary64 residual norm is enclosed from above.
    residual_integral += Fraction(1, 10**15)
    bound = continuation_growth_and_endpoint(0, residual_integral, 0)
    exact = y0 - expm(np.asarray(g * float(h), complex)) @ y0
    exact_norm = float(np.linalg.norm(exact))
    bound_float = float(arb(int(bound['endpoint_row_error_upper']['mantissa']))
                        * arb(2)**int(bound['endpoint_row_error_upper']['exponent']))
    finer = y0 - expm(np.asarray(g * float(h), complex)) @ y0
    return {
        'generator': 'i diag(1,-1)',
        'interpolant': 'constant Y_h = Y(0)',
        'duration': str(h),
        'independent_endpoint_error': exact_norm,
        'majorant': bound_float,
        'majorant_contains_independent_error': bool(bound_float >= exact_norm),
        'refinement_movement_is_not_the_bound': True,
        'refinement_indicator_same_as_exact_for_this_linear_system': float(np.linalg.norm(finer)),
        'physical_source_error_included': False,
    }


def dop853_power_coefficients(dense):
    """Power basis in the step fraction x=(rho-rho_old)/h of a DOP853 interpolant."""
    terms = [np.zeros_like(dense.y_old)]
    for i, stage in enumerate(reversed(dense.F)):
        terms[0] = terms[0] + np.asarray(stage)
        zero = np.zeros_like(dense.y_old)
        if i % 2 == 0:
            terms = [zero] + terms
        else:
            lifted = [c.copy() for c in terms] + [zero]
            for k, coeff in enumerate(terms):
                lifted[k + 1] = lifted[k + 1] - coeff
            terms = lifted
    terms[0] = terms[0] + dense.y_old
    return tuple(terms)


def eval_power_polynomial(coefficients, x):
    acc = np.array(coefficients[-1], copy=True, dtype=complex)
    for coeff in reversed(coefficients[:-1]):
        acc = acc * x + coeff
    return acc


def _rho_ball(lo, hi):
    return slab_ball(min(float(lo), float(hi)), max(float(lo), float(hi)))


def _pauli_acb(matrix):
    return [[complex_ball(matrix[i, j]) for j in range(2)] for i in range(2)]


def _scale_mat(scale, matrix):
    return [[scale * matrix[i][j] for j in range(2)] for i in range(2)]


def _add_mat(left, right):
    return [[left[i][j] + right[i][j] for j in range(2)] for i in range(2)]


def enclosed_generator_split(rho_ball, mass, angular, *, bits=80):
    """G(E) = G0 + E GE with enclosed 1/a, 1/r and exact Pauli matrices."""
    with ctx.workprec(bits):
        bg = background_series(rho_ball, 0)
        inv_a, inv_r, inv_a2 = bg.inv_a[0], arb(1) / bg.r[0], bg.inv_a2[0]
        if not (inv_a.is_finite() and inv_r.is_finite() and inv_a2.is_finite()):
            raise ArithmeticError('finite enclosed 1/a and 1/r required')
        imag = acb(0, 1)
        mass_b, ell = arb(float(mass)), arb(float(angular))
        s1, s2, s3 = _pauli_acb(S1), _pauli_acb(S2), _pauli_acb(S3)
        g0 = _add_mat(_scale_mat(-imag * mass_b * inv_a, s1),
                      _scale_mat(imag * ell * inv_a * inv_r, s2))
        ge = _scale_mat(-imag * inv_a2, s3)
        mu = arb(0)
        for i in range(2):
            for j in range(2):
                herm = (g0[i][j] + g0[j][i].conjugate()) / 2
                mu += herm.abs_upper() ** 2
        return g0, ge, mu.sqrt().upper(), inv_a.abs_upper()


def _matvec2(matrix, vec0, vec1):
    return (matrix[0][0] * vec0 + matrix[0][1] * vec1,
            matrix[1][0] * vec0 + matrix[1][1] * vec1)


def _horner_spinor(coeffs, x, off0, off1):
    y0 = coeffs[-1][off0]
    y1 = coeffs[-1][off1]
    for coeff in reversed(coeffs[:-1]):
        y0 = y0 * x + coeff[off0]
        y1 = y1 * x + coeff[off1]
    d0 = coeffs[-1][off0] * (len(coeffs) - 1)
    d1 = coeffs[-1][off1] * (len(coeffs) - 1)
    for k in range(len(coeffs) - 2, 0, -1):
        d0 = d0 * x + coeffs[k][off0] * k
        d1 = d1 * x + coeffs[k][off1] * k
    return y0, y1, d0, d1


def enclose_dense_step_interval(dense, energies, mass, angular, *, panels=4, bits=80):
    """Per-column a-posteriori residual integrals on one DOP853 step."""
    h = float(dense.h)
    if not np.isfinite(h) or h <= 0:
        raise ValueError('positive DOP853 step required')
    energies = np.asarray(energies, float)
    n_energy = len(energies)
    if dense.y_old.shape != (n_energy * 6,):
        raise ValueError('dense output must match energy x 2 x 3 source fibers')
    power = dop853_power_coefficients(dense)
    with ctx.workprec(bits):
        acb_power = tuple([complex_ball(v) for v in np.ravel(coeff)] for coeff in power)
        residual = [arb(0) for _ in range(n_energy * 3)]
        mu_integral = arb(0)
        for panel in range(panels):
            x0, x1 = panel / panels, (panel + 1) / panels
            x_ball = arb(x0).union(arb(x1))
            rho_lo = float(dense.t_old) + h * x0
            rho_hi = float(dense.t_old) + h * x1
            g0, ge, mu, _ = enclosed_generator_split(
                _rho_ball(rho_lo, rho_hi), mass, angular, bits=bits)
            width = arb(h) * arb(x1 - x0)
            mu_integral += mu * width
            for energy_index, energy in enumerate(energies):
                g = _add_mat(g0, _scale_mat(acb(float(energy)), ge))
                base = energy_index * 6
                for column in range(3):
                    y0, y1, dy0, dy1 = _horner_spinor(
                        acb_power, x_ball, base + column, base + 3 + column)
                    gy0, gy1 = _matvec2(g, y0, y1)
                    r0 = dy0 / arb(h) - gy0
                    r1 = dy1 / arb(h) - gy1
                    residual[energy_index * 3 + column] += (
                        width * (r0.abs_upper()**2 + r1.abs_upper()**2).sqrt())
        return residual, mu_integral.upper(), h


def _background_power_polynomials(rho_start, rho_end, degree, *, bits=80):
    """Power polynomials in x on [0,1] plus uniform Taylor remainders."""
    with ctx.workprec(bits):
        start, end = arb(float(rho_start)), arb(float(rho_end))
        if not end > start:
            raise ValueError('increasing homogeneous continuation step required')
        h, center = end - start, (start + end) / 2
        at_center = background_series(center, degree)
        on_interval = background_series(slab_ball(start, end), degree + 1)
        center_inv_ar = at_center.inv_a / at_center.r
        interval_inv_ar = on_interval.inv_a / on_interval.r
        series = {
            'inv_a2': at_center.inv_a2,
            'inv_a': at_center.inv_a,
            'inv_ar': center_inv_ar,
        }
        interval_series = {
            'inv_a2': on_interval.inv_a2,
            'inv_a': on_interval.inv_a,
            'inv_ar': interval_inv_ar,
        }
        polynomials, errors = {}, {}
        for name, value in series.items():
            coefficients = G._coeffs(value, degree)
            polynomials[name] = [sum((
                coefficients[j] * h**j * comb(j, k) * (-arb(1)/2)**(j-k)
                for j in range(k, degree + 1)), arb(0))
                for k in range(degree + 1)]
            derivative = G._coeffs(interval_series[name], degree + 1)[degree + 1]
            errors[name] = (abs(derivative) * (h/2)**(degree + 1)).upper()
        return polynomials, errors, h


def _poly_convolution(left, right):
    result = [acb(0) for _ in range(len(left) + len(right) - 1)]
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            result[i+j] += acb(x) * y
    return result


def _poly_add(target, source, scale=1):
    if len(target) < len(source):
        target.extend(acb(0) for _ in range(len(source) - len(target)))
    for index, value in enumerate(source):
        target[index] += scale * value


def _vector_bernstein_norm(first, second):
    degree = max(len(first), len(second)) - 1
    first = power_to_bernstein(first, degree)
    second = power_to_bernstein(second, degree)
    return max(((a.abs_upper()**2 + b.abs_upper()**2).sqrt().upper()
                for a, b in zip(first, second)), default=arb(0))


def enclose_dense_step(dense, energies, mass, angular, *, panels=4, bits=80):
    """Cancellation-preserving Bernstein enclosure of one DOP853 defect.

    The dense solution is a degree-seven polynomial in step fraction x.
    Background coefficients are Taylor polynomials on the same x interval.
    Their product is formed before taking norms, then Bernstein convexity
    bounds the vector residual. The analytic generator is anti-Hermitian, so
    its exact 2-logarithmic norm is zero. ``panels`` is retained only for API
    compatibility with the superseded natural-interval diagnostic.
    """
    h = float(dense.h)
    if not np.isfinite(h) or h <= 0:
        raise ValueError('positive DOP853 step required')
    if isinstance(panels, bool) or not isinstance(panels, int) or panels < 1:
        raise ValueError('positive legacy panel count required')
    energies = np.asarray(energies, float)
    n_energy = len(energies)
    if dense.y_old.shape != (n_energy * 6,):
        raise ValueError('dense output must match energy x 2 x 3 source fibers')
    power = dop853_power_coefficients(dense)
    with ctx.workprec(bits):
        coefficients = tuple([complex_ball(value) for value in np.ravel(row)]
                             for row in power)
        geometry, errors, h_ball = _background_power_polynomials(
            dense.t_old, dense.t, 8, bits=bits)
        residual_integrals = []
        for energy_index, energy in enumerate(energies):
            for column in range(3):
                base = energy_index * 6
                y0 = [row[base + column] for row in coefficients]
                y1 = [row[base + 3 + column] for row in coefficients]
                dy0 = [(index + 1) * y0[index + 1] / h_ball
                       for index in range(len(y0) - 1)]
                dy1 = [(index + 1) * y1[index + 1] / h_ball
                       for index in range(len(y1) - 1)]
                r0, r1 = list(dy0), list(dy1)
                # G=-i E/a^2 S3 - i m/a S1 + i ell/(a r) S2.
                diagonal0 = _poly_convolution(geometry['inv_a2'], y0)
                diagonal1 = _poly_convolution(geometry['inv_a2'], y1)
                _poly_add(r0, diagonal0, acb(0, float(energy)))
                _poly_add(r1, diagonal1, -acb(0, float(energy)))
                mass0 = _poly_convolution(geometry['inv_a'], y1)
                mass1 = _poly_convolution(geometry['inv_a'], y0)
                angular0 = _poly_convolution(geometry['inv_ar'], y1)
                angular1 = _poly_convolution(geometry['inv_ar'], y0)
                _poly_add(r0, mass0, acb(0, float(mass)))
                _poly_add(r1, mass1, acb(0, float(mass)))
                # i ell S2 sends (y0,y1) to (ell*y1,-ell*y0).
                _poly_add(r0, angular0, -float(angular))
                _poly_add(r1, angular1, float(angular))
                polynomial = _vector_bernstein_norm(r0, r1)
                y_bound = _vector_bernstein_norm(y0, y1)
                generator_remainder = (
                    abs(arb(float(energy))) * errors['inv_a2']
                    + abs(arb(float(mass))) * errors['inv_a']
                    + abs(arb(float(angular))) * errors['inv_ar'])
                residual_integrals.append(
                    (h_ball * (polynomial + generator_remainder * y_bound)).upper())
        return residual_integrals, arb(0), h


def continue_columns(amplitudes_at_one, energies, mass, angular, *,
                     rho_up=RHO_UP_BINARY64, rtol=CONTINUATION_SOLVER['rtol'],
                     atol=CONTINUATION_SOLVER['atol'],
                     max_step=CONTINUATION_SOLVER['max_step'],
                     panels=4, bits=80):
    """Homogeneous continuation with a directed defect majorant of the interpolant."""
    initial = np.asarray(amplitudes_at_one, complex)
    energies = np.asarray(energies, float)
    if initial.shape != (len(energies), 2, 3) or not np.isfinite(initial).all():
        raise ValueError('finite amplitudes_at_one of shape (n, 2, 3) required')
    if not rho_up_is_archived_binary64(rho_up):
        raise ValueError('archived binary64 rho_up required')

    def rhs(rho, state):
        return (ks_generator(rho, energies, mass, angular)
                @ state.reshape(initial.shape)).ravel()

    solver = DOP853(rhs, 1.0, initial.ravel(), float(rho_up),
                    rtol=rtol, atol=atol, max_step=max_step)
    dense_steps = []
    while solver.status == 'running':
        message = solver.step()
        if solver.status == 'failed':
            raise ArithmeticError(message)
        dense_steps.append(solver.dense_output())
    endpoint = solver.y.reshape(initial.shape)
    with ctx.workprec(bits):
        column_integrals = [arb(0) for _ in range(len(energies) * 3)]
        mu = arb(0)
        for dense in dense_steps:
            residual, mu_step, _ = enclose_dense_step(
                dense, energies, mass, angular, panels=panels, bits=bits)
            for index, value in enumerate(residual):
                column_integrals[index] += value
            mu += mu_step
        growth = mu.exp()
        transported_columns = [growth * value for value in column_integrals]
        frobenius = arb(0)
        for value in transported_columns:
            frobenius += value * value
        frobenius = frobenius.sqrt().upper()
    amplitudes = endpoint.transpose(1, 0, 2).reshape(2, -1)
    return {
        'amplitudes_at_rho_up': amplitudes,
        'nfev': int(solver.nfev),
        'accepted_steps': len(dense_steps),
        'column_residual_integrals': [pack_upper(v) for v in column_integrals],
        'transported_column_errors': [pack_upper(v) for v in transported_columns],
        'unweighted_frobenius_error_upper': pack_upper(frobenius),
        'logarithmic_norm_integral': pack_upper(mu),
        'growth_upper': pack_upper(growth),
        'physical_source_error_included': False,
        'refinement_used_as_bound': False,
    }


def weighted_column_errors(unweighted_column_integrals, weights, energies, *, bits=80):
    weights = np.asarray(weights, float)
    energies = np.asarray(energies, float)
    integrals = tuple(unweighted_column_integrals)
    if weights.shape != (len(integrals),) or energies.shape != weights.shape:
        raise ValueError('matching weights, energies and column integrals required')
    with ctx.workprec(bits):
        field = arb(0)
        axial = arb(0)
        for weight, energy, err in zip(weights, energies, integrals):
            error = restore_packed_upper(err) if isinstance(err, dict) else arb(float(err))
            term = arb(float(weight)) * error
            field += term * term
            axial += (arb(abs(float(energy))) * term) ** 2
        return (packed_upper_fraction(pack_upper(field.sqrt().upper())),
                packed_upper_fraction(pack_upper(axial.sqrt().upper())))


def continuation_matter_bounds(field_norm, axial_norm, field_error, axial_error, *,
                               mass, absolute_angular, multiplicity):
    """N,beta from the continuation field error at incoming background lowers.

    Covariance/source reconstruction error is not included. History radius
    below the background lower bound is a separate field-geometry input.
    """
    bound = matter_error_bounds(
        field_norm, axial_norm, field_error, axial_error,
        mass=mass, absolute_angular=absolute_angular,
        axial_lower=AXIAL_LOWER, radius_lower=RADIUS_LOWER,
        multiplicity=multiplicity, source_operator_error=0)
    return {
        'N': bound['N'],
        'beta': bound['beta'],
        'axial_lower': AXIAL_LOWER,
        'radius_lower': RADIUS_LOWER,
        'source_accuracy_included': False,
        'physical_source_error_included': False,
        'scope': bound['scope'],
    }


def _digest(root, path):
    return sha256((Path(root) / path).read_bytes()).hexdigest()


def _method_kind(method):
    if method.startswith('accepted physical low-mode'):
        return 'low_mode_reconstruction'
    if method.startswith('same original middle producer'):
        return 'middle_mode_columns'
    if method.startswith('eight inherited real subgap'):
        return 'subgap_rows'
    return 'other_source_input'


def classify_low_subgap_windows(inventory_record, panels, coverage):
    """Compact coverage table: ell=0 identity, V5 action regions, unbounded reconstruction."""
    windows = []
    ell0 = tuple(inventory_record['radius_response_exactly_zero_groups'])
    if ell0 != list(ELL0_GROUPS) and ell0 != ELL0_GROUPS:
        raise ValueError('ell=0 radius-response groups changed')
    windows.append({
        'group': 0,
        'sign': 1,
        'kind': 'lll_analytic',
        'panel': None,
        'energy_interval': None,
        'radius_history_response': 'exact_zero',
        'baseline_action_region': 'analytic_LLL_owner',
        'source_reconstruction': 'not_applicable',
        'error_bound': None,
        'identity_is_not_a_source_error_bound': True,
        'notes': 'group 0 has no synthetic mode panel; exact ell=0 zero radius-history response',
    })
    panels = tuple(panels)
    for panel in sorted(panels, key=lambda item: (item.group, item.angular_sign, item.name)):
        name, group, sign = panel.name, panel.group, panel.angular_sign
        energies = np.asarray(panel.energies, float)
        joined = coverage[str(group)]['joined']
        joined_left = float(joined[0])
        method = panel.provenance['method']
        kind = _method_kind(method)
        e_min, e_max = float(energies.min()), float(energies.max())
        action_covered = bool(e_min + 1e-12 >= joined_left)
        windows.append({
            'group': group,
            'sign': sign,
            'kind': kind,
            'panel': name,
            'energy_interval': [e_min, e_max],
            'n_rows': int(len(energies)),
            'radius_history_response': (
                'exact_zero' if group in ELL0_GROUPS else 'nontrivial'),
            'baseline_action_region': (
                'covered_v5_joined' if action_covered else 'not_covered'),
            'v5_joined': joined,
            'source_reconstruction': 'unbounded',
            'error_bound': None,
            'identity_is_not_a_source_error_bound': True,
            'method': method,
            'notes': (
                'V5 covers the baseline action on this energy window; mode columns remain unbounded'
                if action_covered else
                'source reconstruction of stored rho=1 columns is unbounded'),
        })
    for group, row in sorted(coverage.items(), key=lambda item: int(item[0])):
        group = int(group)
        windows.append({
            'group': group,
            'sign': row['actual_signs'],
            'kind': 'covered_baseline_action',
            'panel': None,
            'energy_interval': row['joined'],
            'radius_history_response': (
                'exact_zero' if group in ELL0_GROUPS else 'nontrivial'),
            'baseline_action_region': 'covered_v5_joined',
            'source_reconstruction': 'not_a_mode_column_certificate',
            'error_bound': 'v5_covered_spectral_regions_action_error_upper',
            'identity_is_not_a_source_error_bound': True,
            'tail': row.get('tail'),
            'validated_lower_window': row.get('validated_lower_window'),
            'notes': 'existing V5 action certificate; not a KS mode-column reconstruction bound',
        })
    unbounded = [w for w in windows if w.get('source_reconstruction') == 'unbounded']
    return {
        'ell0_zero_response_groups': list(ELL0_GROUPS),
        'identities_that_are_not_error_bounds': list(IDENTITIES_NOT_ERROR_BOUNDS),
        'windows': windows,
        'unbounded_reconstruction_windows': len(unbounded),
        'unbounded_reconstruction_rows': int(sum(w.get('n_rows', 0) for w in unbounded)),
        'remaining_low_subgap_source_error_bound': None,
        'pauli_or_nested_node_relabelled_as_bound': False,
    }


def family14_signed_multiplicity(channel):
    copies, degeneracy = channel['copy_count'], channel['degeneracy']
    if any(isinstance(v, (bool, np.bool_)) or int(v) != v or v <= 0 for v in (copies, degeneracy)):
        raise ValueError('positive integer copies/degeneracy required')
    signs = 1 if channel['angular_eigenvalue'] == 0 else 2
    return copies * degeneracy / signs, signs
