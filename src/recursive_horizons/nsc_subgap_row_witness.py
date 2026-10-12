"""Replay one archived subgap phase prefix and its Bloch witness.

``subgap_mixed_phase_transport_bound`` reads the original DOP853 solution
only through ``original_dense_solution``. That helper requires ``run.sol`` to
be a closure whose sole free variable is named ``dense`` and whose cell holds
a SciPy ``OdeSolution``. The saved prefix rebuilds that object from the phase
component. The original log-amplitude component is not an input of the phase
proof and is not stored.

This module does not call the Jost solver or the Bloch integrator. Capture of
one original positive row stays in the row-15 driver. Replay validates the
saved prefix and the saved Bloch trace with the existing owners.
"""
import math
from types import SimpleNamespace

import numpy as np
from flint import arb, ctx
from scipy.integrate._ivp.common import OdeSolution
from scipy.integrate._ivp.rk import Dop853DenseOutput

from .nsc_ks_ball_trajectory import exact_upper, restored_upper
from .nsc_massive_jost_mixed_transport import subgap_mixed_phase_transport_bound
from .nsc_massive_jost_transport_bound import (
    original_dense_solution, phase_transport_segments)
from .nsc_metric_horizon_frame import metric_horizon_frame, reflection_from_phase
from .nsc_subgap_source_covariance import (
    BlochSource, initial_bloch, original_covariance_error, validate_bloch)
from .nsc_subgap_upstream_covariance import negative_subgap_covariance_error


PHASE_COLUMNS = (
    'y0', 'y1', 'theta0', 'theta1', 'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6')
PAYLOAD_NAMES = (
    'phase_prefix', 'bloch_trace', 'positive_columns', 'negative_columns',
    'positive_covariance', 'negative_covariance', 'outer_radius',
    'phase_nfev', 'bloch_nfev')
DEGREE = 16
METRIC_TERMS = 48
DEFECT_SUBDIVISIONS = 4
TUBE = '0.00001'
BITS = 192
FRAME_ORDER = 16
BLOCH_DEGREE = 12
BLOCH_Y_START = -18.0
BLOCH_MAX_STEP = 0.025
INNER_REQUEST = 1.01e-4
SOLVER = {'order': 8, 'radial_collar': 1e-10, 'rtol': 2e-13, 'atol': 2e-15}
MISSING_STEP = (
    'unweighted source-covariance errors for one signed row; the quadrature '
    'weight, N/beta contraction, remaining rows and families, and the '
    'physical local gate remain open'
)


class PhasePreparation:
    """Phase proof plus the horizon objects needed to validate a Bloch trace."""

    __slots__ = (
        'bound', 'segments', 'mass_hull', 'angular_hull', 'frame', 'reflection',
        'distance', 'tail', 'initial', 'model', 'target')

    def __init__(self, **kwargs):
        for name in self.__slots__:
            setattr(self, name, kwargs[name])


def requested_inner_radius(horizon_rho):
    """Same inner request as the row-15 mixed-transport pilot."""
    return float(horizon_rho) + INNER_REQUEST


def validate_phase_prefix(prefix):
    """Require finite inward nodes, exact joins, and theta0+F0 == theta1."""
    rows = np.asarray(prefix)
    if (rows.dtype != np.float64 or rows.ndim != 2 or rows.shape[1] != len(PHASE_COLUMNS)
            or rows.shape[0] < 1):
        raise ValueError('phase prefix must be a finite float64 array of 11 columns')
    rows = np.array(rows, dtype=np.float64, copy=True)
    if not np.isfinite(rows).all():
        raise ValueError('phase prefix must be finite')
    if not np.all(rows[:, 1] < rows[:, 0]):
        raise ValueError('phase nodes must step inward')
    if not np.array_equal(rows[:, 2] + rows[:, 4], rows[:, 3]):
        raise ValueError('phase node is not theta0+F0')
    if rows.shape[0] > 1:
        if not np.array_equal(rows[:-1, 1], rows[1:, 0]):
            raise ValueError('phase cells do not join at identical radii')
        if not np.array_equal(rows[:-1, 3], rows[1:, 2]):
            raise ValueError('phase cells do not join at identical angles')
    return rows


def validate_bloch_trace(trace, y_start):
    """Require a finite outward Bloch witness joined to the horizon initializer."""
    rows = np.asarray(trace)
    if rows.dtype != np.float64 or rows.ndim != 2 or rows.shape[1] != 26 or len(rows) < 1:
        raise ValueError('complete Bloch trajectory required')
    rows = np.array(rows, dtype=np.float64, copy=True)
    if not np.isfinite(rows).all():
        raise ValueError('Bloch trajectory must be finite')
    if not np.all(rows[:, 1] > rows[:, 0]):
        raise ValueError('Bloch cells must advance outward')
    if (not np.array_equal(rows[:-1, 1], rows[1:, 0])
            or not np.array_equal(rows[:-1, 5:8], rows[1:, 2:5])):
        raise ValueError('Bloch cells do not join at identical endpoints')
    if rows[0, 0] != float(y_start):
        raise ValueError('trajectory must start at its enclosed horizon initializer')
    return rows


def phase_prefix_from_dense(dense, horizon_rho, inner_radius):
    """Save each transported cell as y0, y1, theta0, theta1 and all seven F's.

    ``phase_transport_segments`` drops F0 because the next node already carries
    it. The prefix keeps F0 so theta0+F0 == theta1 can be checked without the
    solver. Cells inside the inner radius are not part of the phase proof.
    """
    segments = phase_transport_segments(dense, horizon_rho, inner_radius)
    interpolants = getattr(dense, 'interpolants', None)
    if interpolants is None or len(interpolants) < len(segments):
        raise ArithmeticError('DOP853 interpolants ended before the phase cells')
    rows = np.empty((len(segments), len(PHASE_COLUMNS)), dtype=np.float64)
    for index, segment in enumerate(segments):
        interpolant = interpolants[index]
        coefficients = np.asarray(interpolant.F)
        y_old = np.asarray(interpolant.y_old)
        if (coefficients.shape != (7, y_old.shape[0]) or y_old.ndim != 1
                or y_old.shape[0] < 1):
            raise ValueError('unsupported SciPy DOP853 reconstruction layout')
        if np.any(coefficients[:, 0].imag != 0) or y_old[0].imag != 0:
            raise ValueError('real subgap phase data required')
        if (float(interpolant.t_old) != segment.y_start
                or float(y_old[0].real) != segment.theta_start):
            raise ArithmeticError('phase segment is not the next DOP853 cell')
        stored = np.asarray(coefficients[:, 0].real, dtype=np.float64).reshape(7)
        if float(segment.theta_start) + float(stored[0]) != float(segment.theta_end):
            raise ArithmeticError('phase node is not theta0+F0')
        if not np.array_equal(stored[1:], np.asarray(segment.corrections, dtype=np.float64)):
            raise ArithmeticError('saved F coefficients do not match the phase segment')
        rows[index, 0] = segment.y_start
        rows[index, 1] = segment.y_end
        rows[index, 2] = segment.theta_start
        rows[index, 3] = segment.theta_end
        rows[index, 4:] = stored
    return validate_phase_prefix(rows)


def retain_dense_solution(dense):
    """Close over ``dense`` under that exact free name.

    A default argument would hide the name from the closure and
    ``original_dense_solution`` would reject it. The body is not called by
    the phase proof. The omitted amplitude component is not reconstructed.
    """
    def sol(y):
        return dense(y)
    return sol


def ode_solution_from_prefix(prefix):
    """Rebuild the DOP853 phase solution on the saved 11-column prefix."""
    rows = validate_phase_prefix(prefix)
    interpolants = []
    for row in rows:
        coefficients = np.zeros((7, 1), dtype=np.complex128)
        coefficients[:, 0] = row[4:]
        y_old = np.array([row[2]], dtype=np.complex128)
        interpolants.append(Dop853DenseOutput(float(row[0]), float(row[1]), y_old, coefficients))
    nodes = np.empty(len(rows) + 1, dtype=np.float64)
    nodes[0] = rows[0, 0]
    nodes[1:] = rows[:, 1]
    return OdeSolution(nodes, interpolants)


def run_from_prefix(prefix):
    """Return ``run`` whose ``sol`` closure satisfies ``original_dense_solution``."""
    dense = ode_solution_from_prefix(prefix)
    run = SimpleNamespace(sol=retain_dense_solution(dense))
    recovered = original_dense_solution(run)
    if recovered is not dense:
        raise RuntimeError('dense closure did not round-trip')
    return run, dense


def replay_mode(prefix, *, energy, mass, angular, background, outer_radius):
    """Namespace carrying the source labels and the reconstructed solution."""
    run, _dense = run_from_prefix(prefix)
    return SimpleNamespace(
        energy=float(energy), mass=float(mass), angular=float(angular),
        background=background, outer_radius=float(outer_radius), run=run)


def require_same_phase_segments(left, right):
    """Exact endpoint and Hairer-coefficient identity. No tolerance."""
    if len(left) != len(right):
        raise ArithmeticError('phase segment count changed')
    for index, (first, second) in enumerate(zip(left, right)):
        if (first.y_start != second.y_start or first.y_end != second.y_end
                or first.theta_start != second.theta_start
                or first.theta_end != second.theta_end
                or not np.array_equal(first.corrections, second.corrections)):
            raise ArithmeticError(
                f'reconstructed phase segment {index} differs from the original capture')


def require_segments_match_prefix(segments, prefix):
    rows = validate_phase_prefix(prefix)
    if len(segments) != len(rows):
        raise ArithmeticError('phase prefix row count is not the transported cell count')
    for index, (segment, row) in enumerate(zip(segments, rows)):
        if (segment.y_start != row[0] or segment.y_end != row[1]
                or segment.theta_start != row[2] or segment.theta_end != row[3]
                or not np.array_equal(segment.corrections, row[5:])):
            raise ArithmeticError(f'phase segment {index} does not match the saved prefix')


def segments_of_prefix(prefix, horizon_rho, inner_radius):
    run, _dense = run_from_prefix(prefix)
    return phase_transport_segments(
        original_dense_solution(run), float(horizon_rho), float(inner_radius))


def require_outer_node(prefix, outer_radius, horizon_rho):
    node = math.log(float(outer_radius) - float(horizon_rho))
    if not math.isfinite(node) or prefix[0, 0] != node:
        raise ValueError('phase prefix does not start at the source outer radius')


def prepare_phase(prefix, *, energy, mass, angular, background, outer_radius,
                  horizon_rho, kappa, rho_up):
    """Replay the mixed phase owner and the existing horizon frame.

    The mass and angular arguments are the archived positive floats. Their
    hulls are the pilot's unions with pi/2 and sqrt(5). The mode label remains
    the archived float, so a substituted label fails inside the owner.
    """
    prefix = validate_phase_prefix(prefix)
    require_outer_node(prefix, outer_radius, horizon_rho)
    mode = replay_mode(
        prefix, energy=energy, mass=mass, angular=angular, background=background,
        outer_radius=outer_radius)
    inner = requested_inner_radius(horizon_rho)
    segments = phase_transport_segments(
        original_dense_solution(mode.run), horizon_rho, inner)
    require_segments_match_prefix(segments, prefix)
    with ctx.workprec(BITS):
        mass_hull = arb(mass).union(arb.pi() / 2)
        angular_hull = arb(angular).union(arb(5).sqrt())
        bound = subgap_mixed_phase_transport_bound(
            energy, mass_hull, angular_hull, background, mode=mode,
            inner_radius=inner, degree=DEGREE, metric_terms=METRIC_TERMS,
            defect_subdivisions=DEFECT_SUBDIVISIONS, tube=TUBE, bits=BITS)
    if float(segments[-1].y_end).hex() != bound['inner_y_node_hex']:
        raise ArithmeticError('transported node does not match the captured endpoint')
    if int(bound['cells']) != len(segments):
        raise ArithmeticError('transported cell count does not match the saved prefix')
    phase_error = restored_upper(bound['phase_error_inner_upper'])
    with ctx.workprec(BITS):
        rho = arb(horizon_rho) + arb(float(segments[-1].y_end)).exp()
        theta = arb(float(segments[-1].theta_end)) + arb(0, phase_error)
        frame = metric_horizon_frame(
            energy, mass_hull, angular_hull, horizon_rho, order=FRAME_ORDER, bits=BITS)
        reflection, distance, tail = reflection_from_phase(frame, rho, theta)
        initial = initial_bloch(frame, reflection, kappa, BLOCH_Y_START)
        model = BlochSource(frame.q, energy, mass_hull, angular_hull, bits=BITS)
        target = (frame.q - arb.pi() / 2 - arb(rho_up).atan()).log()
    return PhasePreparation(
        bound=bound, segments=segments, mass_hull=mass_hull, angular_hull=angular_hull,
        frame=frame, reflection=reflection, distance=distance, tail=tail,
        initial=initial, model=model, target=target)


def signed_source_errors(positive_columns, positive_covariance, negative_columns,
                         negative_covariance, validated):
    """Recompute each sign. The negative path applies the complement itself."""
    positive = original_covariance_error(
        positive_columns, positive_covariance, validated)
    negative = negative_subgap_covariance_error(
        negative_columns, negative_covariance, validated)
    return positive, negative


def finalize_bloch(prepared, trace, positive_columns, positive_covariance,
                   negative_columns, negative_covariance):
    """Validate a saved trace at rho_up. This does not integrate the trace."""
    rows = validate_bloch_trace(trace, BLOCH_Y_START)
    validated = validate_bloch(
        prepared.model, prepared.initial, rows, prepared.target, degree=BLOCH_DEGREE)
    positive, negative = signed_source_errors(
        positive_columns, positive_covariance, negative_columns, negative_covariance,
        validated)
    return validated, positive, negative


def phase_record(prepared):
    bound = prepared.bound
    return {
        'phase_cells': int(bound['cells']),
        'rate_signs': {key: int(value) for key, value in bound['rate_signs'].items()},
        'min_transport_rate_lower': bound['min_transport_rate_lower'],
        'parameter_intervals': bound['parameter_intervals'],
        'phase_tube_upper': bound['phase_tube_upper'],
        'actual_initializer_phase_hex': bound['actual_initializer_phase_hex'],
        'inner_radius_certified_lower': bound['inner_radius_certified_lower'],
        'inner_offset_certified_lower': bound['inner_offset_certified_lower'],
        'inner_y_node_hex': bound['inner_y_node_hex'],
        'initializer_phase_error_upper': bound['initializer_phase_error_upper'],
        'max_cell_defect_upper': bound['max_cell_defect_upper'],
        'phase_error_inner_upper': bound['phase_error_inner_upper'],
        'reflection_tail_upper': exact_upper(prepared.tail),
        'compact_distance_upper': exact_upper(prepared.distance),
    }


def bloch_record(validated, positive_error, negative_error):
    return {
        'bloch_cells': int(validated['cells']),
        'bloch_error_upper': exact_upper(validated['bloch_error']),
        'initial_bloch_error_upper': exact_upper(validated['initial_error']),
        'normalized_defect_integral_upper': exact_upper(validated['defect_integral']),
        'endpoint_bridge_upper': exact_upper(validated['endpoint_bridge']),
        'positive_covariance_error_upper': exact_upper(positive_error),
        'negative_covariance_error_upper': exact_upper(negative_error),
    }


def require_exact_source(payload, positive_columns, negative_columns,
                         positive_covariance, negative_covariance):
    expected = {
        'positive_columns': positive_columns,
        'negative_columns': negative_columns,
        'positive_covariance': positive_covariance,
        'negative_covariance': negative_covariance,
    }
    for name, original in expected.items():
        value = np.asarray(payload[name])
        if value.shape != original.shape or value.dtype != original.dtype or not np.array_equal(value, original):
            raise ValueError('archived upstream source changed in witness: ' + name)


def _scalar(value, dtype, name):
    array = np.asarray(value)
    if array.dtype != dtype or array.shape != (1,):
        raise ValueError('malformed row witness payload: ' + name)
    if dtype == np.float64 and not np.isfinite(array).all():
        raise ValueError('malformed row witness payload: ' + name)
    return array[0]


def load_witness_arrays(arrays):
    """Check the saved witness inventory before any proof call."""
    if set(arrays) != set(PAYLOAD_NAMES):
        raise ValueError('complete row witness payload required')
    phase = validate_phase_prefix(arrays['phase_prefix'])
    trace = validate_bloch_trace(arrays['bloch_trace'], BLOCH_Y_START)
    for name, shape in (
            ('positive_columns', (2, 3)), ('negative_columns', (2, 3)),
            ('positive_covariance', (3, 3)), ('negative_covariance', (3, 3))):
        value = np.asarray(arrays[name])
        if value.dtype != np.complex128 or value.shape != shape or not np.isfinite(value).all():
            raise ValueError('malformed row witness payload: ' + name)
    outer = float(_scalar(arrays['outer_radius'], np.float64, 'outer_radius'))
    phase_nfev = int(_scalar(arrays['phase_nfev'], np.int64, 'phase_nfev'))
    bloch_nfev = int(_scalar(arrays['bloch_nfev'], np.int64, 'bloch_nfev'))
    if not math.isfinite(outer) or outer <= 0 or phase_nfev < 0 or bloch_nfev < 0:
        raise ValueError('malformed row witness payload: capture scalars')
    return {
        'phase_prefix': phase,
        'bloch_trace': trace,
        'positive_columns': np.array(arrays['positive_columns'], copy=True),
        'negative_columns': np.array(arrays['negative_columns'], copy=True),
        'positive_covariance': np.array(arrays['positive_covariance'], copy=True),
        'negative_covariance': np.array(arrays['negative_covariance'], copy=True),
        'outer_radius': outer,
        'phase_nfev': phase_nfev,
        'bloch_nfev': bloch_nfev,
    }
