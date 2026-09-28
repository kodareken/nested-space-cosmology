"""Directed subgap Jost phase transport on a finite exterior band.

The original solve_jost DOP853 interpolants are the numerical objects. Their
producer and the background metric are unchanged. This owner encloses the
normalized cell defect of the degree-7 phase interpolant and contracts that
defect inward. It does not enclose the near-horizon collar, matched sewing
or a physical rho=1 source column.
"""
from contextlib import contextmanager
from dataclasses import dataclass
import math

import numpy as np
from flint import arb, arb_series, ctx
from scipy.integrate._ivp.common import OdeSolution

from .nsc_ks_ball_trajectory import exact_upper, restored_upper
from .nsc_massive_jost_modes import outgoing_ratio, solve_jost
from .nsc_massive_jost_phase_bound import subgap_phase_bound


MISSING_STEP = (
    'near-horizon collar below the certified inner radius, matched horizon '
    'sewing, and interior transport to rho=1'
)


@contextmanager
def _series_context(bits, degree):
    previous = ctx.cap
    try:
        with ctx.workprec(bits):
            ctx.cap = degree + 1
            yield
    finally:
        ctx.cap = previous


def _variable(center, degree):
    return arb_series([center, 1], prec=degree + 1)


def _coeffs(series, degree):
    values = list(series.coeffs())
    if len(values) < degree + 1:
        values.extend(arb(0) for _ in range(degree + 1 - len(values)))
    return values[:degree + 1]


def _series(values, degree):
    return arb_series(_coeffs(values, degree), prec=degree + 1)


def _endpoint(value, lower=False):
    point = value.lower() if lower else value.upper()
    if not point.is_finite():
        raise ArithmeticError('finite directed endpoint required')
    mantissa, exponent = point.man_exp()
    return {'mantissa': str(mantissa), 'exponent': int(exponent)}


def original_dense_solution(run):
    """Return the OdeSolution retained in solve_jost as the closure name dense."""
    sol = getattr(run, 'sol', None)
    code = getattr(sol, '__code__', None)
    closure = getattr(sol, '__closure__', None)
    if (code is None or closure is None or code.co_freevars != ('dense',)
            or len(closure) != 1):
        raise ValueError('solve_jost must retain OdeSolution in the dense closure')
    dense = closure[0].cell_contents
    if not isinstance(dense, OdeSolution):
        raise TypeError('dense closure is not the original SciPy OdeSolution')
    return dense


@dataclass(frozen=True)
class PhaseTransportSegment:
    """One inward DOP853 cell in y=log(r-horizon), Hairer/TrajectorySegment basis."""
    y_start: float
    y_end: float
    theta_start: float
    theta_end: float
    corrections: object

    def __post_init__(self):
        if (not math.isfinite(self.y_start) or not math.isfinite(self.y_end)
                or not math.isfinite(self.theta_start)
                or not math.isfinite(self.theta_end)
                or not self.y_end < self.y_start):
            raise ValueError('inward finite phase cell required')
        object.__setattr__(self, 'corrections',
                           np.asarray(self.corrections, dtype=float).reshape(6))
        if not np.isfinite(self.corrections).all():
            raise ValueError('finite DOP853 correction coefficients required')


def phase_transport_segments(dense, horizon_rho, inner_radius):
    """Anchored degree-7 phase cells from original interpolants, r_end >= inner_radius.

    Endpoints are solver nodes: the current interpolant y_old and the next
    interpolant y_old. The last included cell uses that exact next node, which
    is also y_old+F[0].
    """
    if not isinstance(dense, OdeSolution) or not dense.interpolants:
        raise ValueError('original DOP853 OdeSolution interpolants required')
    horizon_rho = float(horizon_rho)
    inner_radius = float(inner_radius)
    if not math.isfinite(horizon_rho) or not math.isfinite(inner_radius):
        raise ValueError('finite horizon and inner radius required')
    segments = []
    interpolants = dense.interpolants
    for index, interpolant in enumerate(interpolants):
        F = np.asarray(interpolant.F)
        y_old = np.asarray(interpolant.y_old)
        if F.shape != (7, len(y_old)) or y_old.ndim != 1 or len(y_old) < 1:
            raise ValueError('unsupported SciPy DOP853 reconstruction layout')
        if index + 1 < len(interpolants):
            theta_end_state = np.asarray(interpolants[index + 1].y_old)
            if float(interpolant.t) != float(interpolants[index + 1].t_old):
                raise ArithmeticError('DOP853 phase cells are not contiguous')
            y_end = float(interpolants[index + 1].t_old)
            if not np.array_equal(y_old + F[0], theta_end_state):
                raise ArithmeticError('next interpolant y_old is not the exact node')
        else:
            theta_end_state = y_old + F[0]
            y_end = float(interpolant.t)
        if np.any(F[:,0].imag != 0) or y_old[0].imag != 0 or theta_end_state[0].imag != 0:
            raise ValueError('real subgap phase data required')
        r_end = horizon_rho + math.exp(y_end)
        if not math.isfinite(r_end) or r_end < inner_radius:
            break
        segments.append(PhaseTransportSegment(
            float(interpolant.t_old), y_end, float(y_old[0].real),
            float(theta_end_state[0].real), F[1:, 0].real))
    if not segments:
        raise ArithmeticError('no DOP853 cell remains outside the inner radius')
    return tuple(segments)


def _metric_truncated(inv_r, terms, degree):
    value = arb_series([1], prec=degree + 1)
    for j in range(terms):
        value -= arb(6)*(-1)**j*inv_r**(2*j+1)/((2*j+1)*(2*j+3))
    return _series(value, degree)


def _composed_metric_tail(inv_r, terms, degree):
    """Faà-di-Bruno majorant of the omitted 1/r series, as a series in the cell.

    D_k = k-th u-derivative of 6 u^{2J+1}/(1-u) at u=1/r_min, composed with the
    positive Taylor series of |d^j(1/r)/dx^j|/j!. These are not u-derivatives
    copied into the y (or cell-x) basis.
    """
    u_max = inv_r[0].abs_upper()
    if not u_max < 1:
        raise ArithmeticError('metric series coordinate 1/r is not inside the unit disk')
    majorant = arb(6)*_variable(u_max, degree)**(2*terms+1)/(1-_variable(u_max, degree))
    generating = arb_series([arb(majorant[k].abs_upper()) for k in range(degree + 1)],
                            prec=degree + 1)
    increment = arb_series(
        [arb(0)] + [inv_r[k].abs_upper() for k in range(1, degree + 1)],
        prec=degree + 1)
    composed = arb_series([0], prec=degree + 1)
    for k in range(degree, -1, -1):
        composed = composed*increment + generating[k]
    return arb_series([arb(0, composed[k].abs_upper()) for k in range(degree + 1)],
                      prec=degree + 1)


def metric_series(inv_r, terms, degree):
    """Exact truncated 1/r series plus composed derivative-tail majorant."""
    return _metric_truncated(inv_r, terms, degree) + _composed_metric_tail(
        inv_r, terms, degree)


def _hairer_series(theta_start, theta_end, corrections, x, degree):
    value = arb_series([0], prec=degree + 1)
    terms = [arb(theta_end) - arb(theta_start),
             *[arb(float(coefficient)) for coefficient in corrections]]
    for index, term in enumerate(reversed(terms)):
        value = value + term
        value = value*x if index % 2 == 0 else value*(1 - x)
    return _series(value + arb(theta_start), degree)


def _metric_value(radius, terms):
    """Value enclosure of A(r) from the 1/r series; no y-derivative reuse."""
    if not radius > 1:
        raise ArithmeticError('metric series requires r>1')
    u = 1/radius
    value = arb(1)
    for j in range(terms):
        value -= arb(6)*(-1)**j*u**(2*j+1)/((2*j+1)*(2*j+3))
    upper = u.abs_upper()
    if not upper < 1:
        raise ArithmeticError('metric series coordinate 1/r is not inside the unit disk')
    majorant = arb(6)*upper**(2*terms+1)/(1-upper)
    return value + arb(0, majorant.abs_upper())


def _phase_rhs_series(x, segment, energy, mass, angular, horizon, terms, degree):
    y = arb(segment.y_start) + (arb(segment.y_end) - arb(segment.y_start))*x
    radius = arb(horizon) + y.exp()
    inv_r = radius.inv()
    metric = metric_series(inv_r, terms, degree)
    theta = _hairer_series(segment.theta_start, segment.theta_end,
                           segment.corrections, x, degree)
    root = metric.sqrt()
    sphere = (1 + radius*radius).sqrt()
    return (2*y.exp()/metric*(root*(angular/sphere*theta.cos()
                                    + mass*theta.sin()) - energy),
            metric, theta, radius)


def _normalized_defect(x, segment, energy, mass, angular, horizon, terms, degree):
    step = (arb(segment.y_end) - arb(segment.y_start))
    rhs, metric, theta, radius = _phase_rhs_series(
        x, segment, energy, mass, angular, horizon, terms, degree)
    theta_x = _series(theta.derivative(), degree)
    return _series(theta_x - step*rhs, degree), metric, theta, radius, step


def _contraction(initial, residual, rate):
    decay = (-rate).exp()
    quotient = -(-rate).expm1()/rate
    return decay*initial + residual*quotient


def enclose_phase_transport_cell(segment, energy, mass, angular, horizon_rho, *,
                                 bits=192, degree=12, metric_terms=12,
                                 tube='0.00001', defect_subdivisions=1):
    """Enclose a cell using maxima/minima over an exact dyadic partition.

    Subdivision changes the enclosure only. Every subinterval bounds the SAME
    normalized defect theta_x-h*f_y; no second step/partition width is inserted
    into the whole-cell contraction estimate.
    """
    if not isinstance(segment, PhaseTransportSegment):
        raise TypeError('PhaseTransportSegment required')
    _validate_orders(degree, metric_terms, bits)
    if (isinstance(defect_subdivisions, bool) or not isinstance(defect_subdivisions, int)
            or defect_subdivisions not in (1, 2, 4, 8, 16)):
        raise ValueError('defect subdivisions must be one of 1,2,4,8,16')
    with _series_context(bits, degree):
        energy, mass, angular, horizon = map(arb, (energy, mass, angular, horizon_rho))
        eta = arb(tube)
        if (not all(v.is_finite() for v in (energy, mass, angular, horizon, eta))
                or not eta > 0 or not eta < arb(1)/4):
            raise ValueError('finite channel parameters and a small phase tube required')
        step = arb(segment.y_end)-arb(segment.y_start)
        if not step < 0:
            raise ValueError('inward phase cell required')
        residuals, rates, metrics = [], [], []
        for index in range(defect_subdivisions):
            left, right = arb(index)/defect_subdivisions, arb(index+1)/defect_subdivisions
            center_x, half_width = (left+right)/2, (right-left)/2
            center, _, _, _, _ = _normalized_defect(
                _variable(center_x, degree), segment, energy, mass, angular, horizon,
                metric_terms, degree)
            interval, metric, theta, radius, _ = _normalized_defect(
                _variable(left.union(right), degree), segment, energy, mass, angular,
                horizon, metric_terms, degree)
            residual = sum((center[k].abs_upper()*half_width**k
                            for k in range(degree)), arb(0))
            residual = (residual+interval[degree].abs_upper()*half_width**degree).upper()
            y_ball = arb(segment.y_start)+step*left.union(right)
            r_ball = horizon+y_ball.exp()
            a_ball = _metric_value(r_ball, metric_terms)
            theta_ball = theta[0]+arb(0,eta.upper())
            bracket = (mass*theta_ball.cos()
                       -angular*theta_ball.sin()/(1+r_ball*r_ball).sqrt())
            lam = (2*y_ball.exp()/a_ball.sqrt()*bracket).lower()
            if (not residual.is_finite() or not a_ball > 0 or not metric[0] > 0
                    or not r_ball > 1 or not radius[0] > 1 or not lam > 0):
                raise ArithmeticError('cell metric or phase monotonicity is unresolved')
            residuals.append(residual)
            rates.append(lam)
            metrics.append(a_ball.lower())
        residual = max(residuals)
        lam = min(rates)
        rate = ((-step)*lam).lower()
        if not rate > 0:
            raise ArithmeticError('positive contraction rate is unresolved')
        return {'defect_upper': residual, 'phase_derivative_lower': lam,
                'contraction_rate_lower': rate, 'step': step, 'tube_upper': eta,
                'metric_lower': min(metrics), 'defect_subdivisions': defect_subdivisions}


def _validate_orders(degree, metric_terms, bits):
    if (isinstance(degree, bool) or not isinstance(degree, int) or degree < 7
            or isinstance(metric_terms, bool) or not isinstance(metric_terms, int)
            or metric_terms < 1 or 2*metric_terms + 1 < degree
            or isinstance(bits, bool) or not isinstance(bits, int) or bits < 80):
        raise ValueError('valid Taylor degree, metric tail order and precision required')


def _validate_proof_parameters(degree, metric_terms, bits, inner_radius, horizon_rho):
    _validate_orders(degree, metric_terms, bits)
    inner_radius = float(inner_radius)
    horizon_rho = float(horizon_rho)
    if (not math.isfinite(inner_radius) or not math.isfinite(horizon_rho)
            or not inner_radius > 1 or not inner_radius > horizon_rho + 1e-4):
        raise ValueError(
            'certified inner radius must satisfy r>1 and lie at least 1e-4 outside the coordinate horizon')
    return inner_radius, horizon_rho


def subgap_phase_transport_bound(energy, mass, angular, background, *,
                                 inner_radius=8.0, bits=192, degree=12,
                                 metric_terms=12, tube='0.00001', mode=None, defect_subdivisions=1):
    """Bound |theta_true - theta_interpolant| from the solver R down to inner r.

    The independent variable is the original y=log(r-horizon). Forward cell
    fraction x in [0,1] has negative width h, contraction k=-h f_{y,theta}, and
    endpoint bound exp(-k) ein + eps (1-exp(-k))/k inside a validated tube.
    Assumptions that are not proved raise ArithmeticError (fail closed).
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
        enclosure = ball if isinstance(ball, arb) else arb(float(ball))
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
    with _series_context(bits, degree):
        # log/exp are not an exact round trip. Bind the mathematical radius at
        # the stored initial y node, and the ACTUAL solver's initializer angle.
        actual_outer = arb(horizon_rho)+arb(segments[0].y_start).exp()
        initializer = subgap_phase_bound(
            energy, mass, angular, actual_outer, coefficients, segments[0].theta_start,
            bits=bits, degree=degree, metric_terms=metric_terms, tube=tube)
        error = restored_upper(initializer['combined_phase_initial_error_upper'])
        eta = arb(tube)
        if not error < eta.lower():
            raise ArithmeticError('initializer already leaves the phase tube')
        for segment in segments:
            cell = enclose_phase_transport_cell(
                segment, energy, mass, angular, horizon_rho, bits=bits,
                degree=degree, metric_terms=metric_terms, tube=tube,
                defect_subdivisions=defect_subdivisions)
            if not cell['defect_upper'] < eta.lower():
                raise ArithmeticError('normalized cell defect is not inside the tube')
            transported = _contraction(error, cell['defect_upper'],
                                       cell['contraction_rate_lower'])
            peak = error.union(transported).upper()
            if not peak.is_finite() or not peak < eta.lower():
                raise ArithmeticError('phase tube is not invariant under this cell')
            error = transported.upper()
            cell_defects.append(cell['defect_upper'])
        inner_y = arb(segments[-1].y_end)
        inner_r = arb(horizon_rho) + inner_y.exp()
        inner_offset = inner_r - arb(horizon_rho)
        if not inner_r >= inner_radius or not inner_offset >= arb('1e-4'):
            raise ArithmeticError('certified endpoint is not outside the declared floor')
        max_defect = arb(0)
        for item in cell_defects:
            upper = item.abs_upper()
            if upper > max_defect:
                max_defect = upper
        return {
            'scope': 'decaying subgap phase transport on a finite exterior band',
            'metric': 'A=1-3*((1+r^2)*atan(1/r)-r)',
            'independent_variable': 'y=log(r-horizon_rho)',
            'f_y': '2*exp(y)/A*(sqrt(A)*(ell/sqrt(1+r^2)*cos(theta)+m*sin(theta))-E)',
            'normalized_defect': 'theta_x-h*f_y on the DOP853 cell fraction x',
            'interpolant': (
                'degree-7 Hairer form from y_old, next interpolant y_old, and F[1:]'),
            'contraction': 'k=-h*f_y,theta; exp(-k)*ein+eps*(1-exp(-k))/k',
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
            'outer_radius': outer_radius,
            'actual_outer_radius_interval': {'lower': _endpoint(actual_outer, True),
                                             'upper': _endpoint(actual_outer)},
            'taylor_center': 'subinterval midpoint',
            'defect_subdivisions': defect_subdivisions,
            'actual_initializer_phase_hex': segments[0].theta_start.hex(),
            'inner_radius_certified_lower': _endpoint(inner_r, True),
            'inner_offset_certified_lower': _endpoint(inner_offset, True),
            'inner_y_node_hex': float(segments[-1].y_end).hex(),
            'initializer_phase_error_upper': initializer['combined_phase_initial_error_upper'],
            'max_cell_defect_upper': exact_upper(max_defect),
            'phase_error_inner_upper': exact_upper(error),
            'all_exterior_radii': False,
            'certified_finite_band': True,
            'horizon_sewing_error': None,
            'physical_rho1_source_error': None,
            'missing_step': MISSING_STEP,
            'initializer': initializer,
        }
