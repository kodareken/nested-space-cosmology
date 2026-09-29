"""Signed one-sided subgap Jost phase transport on a finite exterior band.

Successor of ``nsc_massive_jost_transport_bound``. The normalized defect,
the single original cell width, the composed metric tails and the exterior
half-line initializer are the same objects. The phase derivative is only
required to have a finite lower bound: inward transport may expand.

The predecessor remains unmodified. This module does not replace source
columns, sew the horizon, or close the local gate.
"""
import math

from flint import arb

from .nsc_ks_ball_trajectory import exact_upper, restored_upper
from .nsc_massive_jost_modes import outgoing_ratio, solve_jost
from .nsc_massive_jost_phase_bound import subgap_phase_bound
from .nsc_massive_jost_transport_bound import (
    PhaseTransportSegment,
    _endpoint,
    _metric_value,
    _normalized_defect,
    _series_context,
    _validate_orders,
    _validate_proof_parameters,
    _variable,
    original_dense_solution,
    phase_transport_segments,
)


MISSING_STEP = (
    'this phase owner stops at the declared inner radius; matched horizon '
    'sewing and the physical rho=1 source column are not replaced'
)


def _explicit_real(value, label):
    # Arb(None) is zero: absent proof data must be rejected before conversion.
    if value is None or isinstance(value, bool):
        raise ValueError('explicit ' + label + ' required')
    return arb(value)


def _phi_tail(magnitude, order):
    """Geometric majorant of sum_{n>=order} M^n/(n+1)!."""
    magnitude = _explicit_real(magnitude, 'magnitude')
    if magnitude.is_zero():
        return arb(0)
    if not magnitude >= 0 or not magnitude < order + 2:
        raise ArithmeticError('analytic phi tail is unresolved')
    denominator = arb(1)
    for factor in range(1, order + 2):
        denominator *= factor
    base = magnitude**order / denominator
    return (base / (1 - magnitude / (order + 2))).upper()


def _phi_series(rate):
    """Integrate the exponential series. No division by the rate."""
    rate = _explicit_real(rate, 'rate')
    magnitude = rate.abs_upper()
    if not magnitude.is_finite():
        raise ArithmeticError('analytic phi weight is unresolved')
    order = 1
    tail = arb(1)
    while True:
        if not magnitude < order + 2:
            order *= 2
        else:
            tail = _phi_tail(magnitude, order)
            if tail < arb(2)**(-160):
                break
            order *= 2
        if order > 4096:
            raise ArithmeticError('analytic phi tail is unresolved')
    term = arb(1)
    total = term
    for index in range(order - 1):
        term *= (-rate) / (index + 2)
        total += term
    value = total + arb(0, tail)
    if not value.is_finite():
        raise ArithmeticError('analytic phi weight is unresolved')
    return value


def _phi_atom(point):
    point = _explicit_real(point, 'rate')
    if not point.is_finite():
        raise ArithmeticError('finite transport rate required')
    if not (point > 0 or point < 0):
        return _phi_series(point)
    # (1-exp(-k))/k = -expm1(-k)/k, with k bounded away from 0.
    value = (-(point.neg().expm1())) / point
    if not value.is_finite() or not value > 0:
        raise ArithmeticError('analytic phi weight is unresolved')
    return value


def phi_weight(rate):
    """Enclose φ(k)=∫_0^1 exp(-k t) dt on a real ball.

    φ(0)=1 and φ(k)=(1-exp(-k))/k otherwise. The function is entire and
    strictly decreasing. A ball containing 0 is enclosed by the integrated
    exponential series when it is small, and by its endpoint values when it
    is wide. Neither branch divides by an interval containing 0.
    """
    rate = _explicit_real(rate, 'rate')
    if not rate.is_finite():
        raise ArithmeticError('finite transport rate required')
    if rate > 0 or rate < 0:
        return _phi_atom(rate)
    if rate.abs_upper() < 1:
        return _phi_series(rate)
    return _phi_atom(rate.lower()).union(_phi_atom(rate.upper()))


def _as_nonnegative(value, label):
    value = _explicit_real(value, label)
    if not value.is_finite() or not value >= 0:
        raise ArithmeticError(f'nonnegative finite {label} is unresolved')
    return value


def signed_transport_image(initial, residual, rate):
    """Enclose exp(-k)*ein + eps*φ(k) at the left endpoint of ``rate``.

    Both factors decrease in k, so the upper bound sits at the smaller
    admissible rate. ``k=0`` and a rate ball containing 0 stay finite.
    """
    initial = _as_nonnegative(initial, 'initial error')
    residual = _as_nonnegative(residual, 'cell defect')
    rate = _explicit_real(rate, 'rate')
    if not rate.is_finite():
        raise ArithmeticError('finite transport rate required')
    k = rate.lower()
    image = (-k).exp() * initial + residual * phi_weight(k)
    if not image.is_finite():
        raise ArithmeticError('signed transport bound is unresolved')
    return image


def comparison_value(initial, residual, rate, fraction):
    """Value of the constant-rate majorant at a cell fraction in [0,1].

    u'(x)=exp(-k x)*(eps-k*ein) has constant sign, so the maximum on the
    cell is attained at an endpoint.
    """
    fraction = _explicit_real(fraction, 'comparison fraction')
    if not fraction.is_finite() or not fraction >= 0 or not fraction <= 1:
        raise ValueError('comparison fraction must lie in [0,1]')
    initial = _as_nonnegative(initial, 'initial error')
    residual = _as_nonnegative(residual, 'cell defect')
    k = _explicit_real(rate, 'rate').lower()
    if not k.is_finite():
        raise ArithmeticError('finite transport rate required')
    span = k * fraction
    image = (-span).exp() * initial + residual * fraction * phi_weight(span)
    if not image.is_finite():
        raise ArithmeticError('comparison value is unresolved')
    return image


def _minimum_lower(values):
    acc = values[0]
    for item in values[1:]:
        acc = acc.union(item)
    point = acc.lower()
    if not point.is_finite():
        raise ArithmeticError('finite lower bound required')
    return point


def _maximum_upper(values):
    acc = values[0]
    for item in values[1:]:
        acc = acc.union(item)
    point = acc.upper()
    if not point.is_finite():
        raise ArithmeticError('finite upper bound required')
    return point


def enclose_mixed_phase_transport_cell(segment, energy, mass, angular, horizon_rho, *,
                                       bits=192, degree=12, metric_terms=12,
                                       tube='0.00001', defect_subdivisions=1):
    """Enclose one inward cell with a signed phase-derivative lower bound.

    Dyadic subdivisions change only the Taylor majorant of the normalized
    residual theta_x-h*f_y. The original cell width h is used once.
    """
    if not isinstance(segment, PhaseTransportSegment):
        raise TypeError('PhaseTransportSegment required')
    _validate_orders(degree, metric_terms, bits)
    if (isinstance(defect_subdivisions, bool) or not isinstance(defect_subdivisions, int)
            or defect_subdivisions not in (1, 2, 4, 8, 16)):
        raise ValueError('defect subdivisions must be one of 1,2,4,8,16')
    with _series_context(bits, degree):
        energy, mass, angular, horizon = (_explicit_real(v, n) for v, n in zip(
            (energy, mass, angular, horizon_rho), ('energy', 'mass', 'angular', 'horizon')))
        eta = _explicit_real(tube, 'phase tube')
        if (not all(v.is_finite() for v in (energy, mass, angular, horizon, eta))
                or not eta > 0 or not eta < arb(1) / 4):
            raise ValueError('finite channel parameters and a small phase tube required')
        step = arb(segment.y_end) - arb(segment.y_start)
        if not step < 0:
            raise ValueError('inward phase cell required')
        residuals, rates, metrics = [], [], []
        for index in range(defect_subdivisions):
            left, right = arb(index) / defect_subdivisions, arb(index + 1) / defect_subdivisions
            center_x, half_width = (left + right) / 2, (right - left) / 2
            center, _, _, _, _ = _normalized_defect(
                _variable(center_x, degree), segment, energy, mass, angular, horizon,
                metric_terms, degree)
            interval, metric, theta, radius, _ = _normalized_defect(
                _variable(left.union(right), degree), segment, energy, mass, angular,
                horizon, metric_terms, degree)
            residual = sum((center[k].abs_upper() * half_width**k
                            for k in range(degree)), arb(0))
            residual = (residual + interval[degree].abs_upper() * half_width**degree).upper()
            y_ball = arb(segment.y_start) + step * left.union(right)
            r_ball = horizon + y_ball.exp()
            a_ball = _metric_value(r_ball, metric_terms)
            theta_ball = theta[0] + arb(0, eta.upper())
            bracket = (mass * theta_ball.cos()
                       - angular * theta_ball.sin() / (1 + r_ball * r_ball).sqrt())
            if (not residual.is_finite() or not a_ball > 0 or not metric[0] > 0
                    or not r_ball > 1 or not radius[0] > 1):
                raise ArithmeticError('cell metric is unresolved')
            derivative = 2 * y_ball.exp() / a_ball.sqrt() * bracket
            lam = derivative.lower()
            if not lam.is_finite():
                raise ArithmeticError('cell phase derivative is unresolved')
            residuals.append(residual)
            rates.append(lam)
            metrics.append(a_ball.lower())
        lam = _minimum_lower(rates)
        rate = ((-step) * lam).lower()
        if not rate.is_finite():
            raise ArithmeticError('cell phase derivative is unresolved')
        return {
            'defect_upper': _maximum_upper(residuals),
            'phase_derivative_lower': lam,
            'transport_rate_lower': rate,
            'step': step,
            'tube_upper': eta,
            'metric_lower': _minimum_lower(metrics),
            'defect_subdivisions': defect_subdivisions,
            'normalized_defect': 'theta_x-h*f_y',
            'cell_width': 'exact endpoint difference, inserted once',
        }


def advance_mixed_tube(error, cell, tube):
    """Advance one cell, or raise if the comparison leaves the phase tube.

    The constant-rate majorant peaks at an endpoint, so the initial error
    and the transported end are the only tube checks.
    """
    eta = _explicit_real(tube, 'phase tube')
    if not eta.is_finite() or not eta > 0:
        raise ValueError('positive phase tube required')
    defect = cell['defect_upper']
    if not defect.is_finite() or not defect < eta.lower():
        raise ArithmeticError('normalized cell defect is not inside the tube')
    transported = signed_transport_image(error, defect, cell['transport_rate_lower'])
    peak = arb(error).union(transported).upper()
    if not peak.is_finite() or not peak < eta.lower():
        raise ArithmeticError('phase tube is not invariant under this cell')
    return transported.upper()


def _rate_sign(rate):
    if rate > 0:
        return 'contracting'
    if rate < 0:
        return 'expanding'
    if rate.is_zero():
        return 'zero'
    return 'straddling'


def subgap_mixed_phase_transport_bound(energy, mass, angular, background, *,
                                       inner_radius=8.0, bits=192, degree=12,
                                       metric_terms=12, tube='0.00001', mode=None,
                                       defect_subdivisions=1):
    """Bound |theta_true-theta_interpolant| down to a declared inner radius.

    ``k=-h*lambda`` may be negative, zero, or an enclosure containing zero.
    Failure of the metric enclosure or of tube invariance raises.
    """
    inner_radius, horizon_rho = _validate_proof_parameters(
        degree, metric_terms, bits, inner_radius, background.horizon_rho)
    if mode is None:
        if any(isinstance(value, arb) for value in (energy, mass, angular)):
            raise ValueError('Arb parameter intervals require a solved Jost mode')
        energy_float, mass_float, angular_float = map(float, (energy, mass, angular))
        if not all(math.isfinite(v) for v in (energy_float, mass_float, angular_float)):
            raise ValueError('finite energy, mass and angular channel required')
        if energy_float <= 0 or not mass_float > energy_float:
            raise ValueError('subgap positive mass required')
        mode = solve_jost(background, energy_float, mass_float, angular_float)
    energy_float, mass_float, angular_float = (
        float(mode.energy), float(mode.mass), float(mode.angular))
    if energy_float <= 0 or not mass_float > energy_float:
        raise ValueError('subgap positive mass required')
    for name, ball, number in (
            ('energy', energy, energy_float), ('mass', mass, mass_float),
            ('angular', angular, angular_float)):
        enclosure = _explicit_real(ball, name)
        if not enclosure.contains(arb(number)):
            raise ValueError(f'{name} interval does not contain the Jost label')
    if float(mode.background.horizon_rho) != horizon_rho:
        raise ValueError('mode and proof must use the same radial coordinate origin')
    dense = original_dense_solution(mode.run)
    segments = phase_transport_segments(dense, horizon_rho, inner_radius)
    outer_radius = float(mode.outer_radius)
    _, _, coefficients = outgoing_ratio(
        energy_float, mass_float, angular_float, outer_radius, order=8)
    cell_defects = []
    signs = {'contracting': 0, 'expanding': 0, 'zero': 0, 'straddling': 0}
    with _series_context(bits, degree):
        actual_outer = arb(horizon_rho) + arb(segments[0].y_start).exp()
        initializer = subgap_phase_bound(
            energy, mass, angular, actual_outer, coefficients, segments[0].theta_start,
            bits=bits, degree=degree, metric_terms=metric_terms, tube=tube)
        error = restored_upper(initializer['combined_phase_initial_error_upper'])
        eta = _explicit_real(tube, 'phase tube')
        if not error < eta.lower():
            raise ArithmeticError('initializer already leaves the phase tube')
        min_rate = None
        for index, segment in enumerate(segments):
            cell = enclose_mixed_phase_transport_cell(
                segment, energy, mass, angular, horizon_rho, bits=bits,
                degree=degree, metric_terms=metric_terms, tube=tube,
                defect_subdivisions=defect_subdivisions)
            rate = cell['transport_rate_lower']
            signs[_rate_sign(rate)] += 1
            min_rate = rate if min_rate is None else min_rate.union(rate).lower()
            try:
                error = advance_mixed_tube(error, cell, eta)
            except ArithmeticError as exc:
                raise ArithmeticError(
                    f'{exc} (cell {index} of {len(segments)})') from exc
            cell_defects.append(cell['defect_upper'])
        inner_y = arb(segments[-1].y_end)
        inner_r = arb(horizon_rho) + inner_y.exp()
        inner_offset = inner_r - arb(horizon_rho)
        if not inner_r >= inner_radius or not inner_offset >= arb('1e-4'):
            raise ArithmeticError('certified endpoint is not outside the declared floor')
        max_defect = _maximum_upper(cell_defects)
        return {
            'scope': 'signed subgap phase transport on a finite exterior band',
            'predecessor': 'nsc_massive_jost_transport_bound.subgap_phase_transport_bound',
            'metric': 'A=1-3*((1+r^2)*atan(1/r)-r)',
            'independent_variable': 'y=log(r-horizon_rho)',
            'f_y': '2*exp(y)/A*(sqrt(A)*(ell/sqrt(1+r^2)*cos(theta)+m*sin(theta))-E)',
            'f_y_theta': (
                '2*exp(y)/sqrt(A)*(m*cos(theta)-ell*sin(theta)/sqrt(1+r^2))'),
            'normalized_defect': 'theta_x-h*f_y on the DOP853 cell fraction x',
            'interpolant': (
                'degree-7 Hairer form from y_old, next interpolant y_old, and F[1:]'),
            'transport': (
                'k=-h*lambda for signed lower lambda; '
                'exp(-k)*ein+eps*phi(k), phi=integral_0^1 exp(-k*t) dt'),
            'peak': 'constant-rate comparison attains its maximum at an endpoint',
            'parameter_intervals': {
                name: {'lower': _endpoint(v, True), 'upper': _endpoint(v)}
                for name, v in [
                    ('energy', arb(energy)), ('mass', arb(mass)),
                    ('angular', arb(angular)),
                    ('outer_radius', arb(outer_radius)),
                    ('inner_radius', inner_r)]},
            'bits': bits, 'taylor_degree': degree, 'metric_terms': metric_terms,
            'phase_tube_upper': exact_upper(eta),
            'cells': len(segments),
            'rate_signs': signs,
            'min_transport_rate_lower': {
                'lower': _endpoint(min_rate, True), 'upper': _endpoint(min_rate)},
            'outer_radius': outer_radius,
            'actual_outer_radius_interval': {
                'lower': _endpoint(actual_outer, True), 'upper': _endpoint(actual_outer)},
            'taylor_center': 'subinterval midpoint',
            'defect_subdivisions': defect_subdivisions,
            'actual_initializer_phase_hex': segments[0].theta_start.hex(),
            'inner_radius_certified_lower': _endpoint(inner_r, True),
            'inner_offset_certified_lower': _endpoint(inner_offset, True),
            'inner_y_node_hex': float(segments[-1].y_end).hex(),
            'initializer_phase_error_upper': initializer[
                'combined_phase_initial_error_upper'],
            'max_cell_defect_upper': exact_upper(max_defect),
            'phase_error_inner_upper': exact_upper(error),
            'all_exterior_radii': False,
            'certified_finite_band': True,
            'horizon_sewing_error': None,
            'physical_rho1_source_error': None,
            'source_columns_replaced': False,
            'missing_step': MISSING_STEP,
            'initializer': initializer,
        }
