"""Homogeneous Bloch transport of one interior vacuum branch.

The horizon frame is built inside this owner from the same positive energy,
mass, angular label and horizon hint as the generator. An unrelated frame
cannot be supplied. VacuumCorrection is not used: its expansion forcing is
not subtracted. Finite occupation stays outside the archived-matrix distance.

This is a reusable per-row method. It is not full source coverage, a tail
remainder, or local-gate closure.
"""
import time

import numpy as np
from flint import arb, ctx
from scipy.integrate import solve_ivp

from .nsc_metric_horizon_frame import metric_horizon_frame
from .nsc_source_occupation_enclosure import occupation_vacuum_distance
from .nsc_subgap_source_covariance import (
    BlochSource, _coefficients, original_covariance_distance_bounds, validate_bloch,
)
from .nsc_vacuum_source_remainder import _real


def _explicit_real(value, name):
    """Reject missing and boolean inputs before arb can turn them into zero."""
    if value is None or isinstance(value, (bool, np.bool_)):
        raise ValueError('explicit real '+name+' required')
    if isinstance(value, (list, tuple, dict, set, np.ndarray)):
        raise ValueError('explicit real '+name+' required')
    try:
        return _real(value, name)
    except (TypeError, ValueError) as exc:
        raise ValueError('explicit real '+name+' required') from exc


def _declared_float(value, name):
    number = _explicit_real(value, name)
    if not number.is_finite():
        raise ValueError('finite '+name+' required')
    out = float(number.mid())
    if not np.isfinite(out):
        raise ValueError('finite '+name+' required')
    return out


def _positive_control(value, name):
    if value is None or isinstance(value, (bool, np.bool_)):
        raise ValueError('explicit positive '+name+' required')
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError('explicit positive '+name+' required') from exc
    if not np.isfinite(number) or number <= 0:
        raise ValueError('explicit positive '+name+' required')
    return number


def _archived_matrix(value, shape, name):
    if value is None or isinstance(value, (bool, np.bool_)):
        raise ValueError('explicit archived '+name+' required')
    if isinstance(value, np.ndarray):
        array = value
    elif isinstance(value, (list, tuple)):
        array = np.array(value)
    else:
        raise ValueError('explicit archived '+name+' required')
    if array.dtype == object or np.issubdtype(array.dtype, np.bool_):
        raise ValueError('explicit numeric archived '+name+' required')
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError('finite archived '+name+' with shape %s required' % (shape,))
    return array


class DirectVacuumBloch:
    """Positive-energy vacuum projector transported by the homogeneous Bloch law.

    Any positive energy inside the frame and metric domains is allowed,
    including E>m. Mass is nonnegative. Angular momentum is signed.
    """

    def __init__(self, energy, mass, angular, horizon_hint, *, bits=192,
                 frame_order=16, metric_terms=48, analytic_radius='0.1'):
        if type(bits) is not int or bits < 80:
            raise ValueError('directed precision >=80 required')
        if type(frame_order) is not int or frame_order < 2:
            raise ValueError('integer frame order >=2 required')
        if type(metric_terms) is not int or metric_terms < 16:
            raise ValueError('metric tail terms >=16 required')
        with ctx.workprec(bits):
            energy = _explicit_real(energy, 'energy')
            mass = _explicit_real(mass, 'mass')
            angular = _explicit_real(angular, 'angular')
            hint = _explicit_real(horizon_hint, 'horizon hint')
            radius = _explicit_real(analytic_radius, 'analytic radius')
            if not energy.is_finite() or not energy > 0:
                raise ValueError('positive finite energy required')
            if not mass.is_finite() or not mass >= 0:
                raise ValueError('nonnegative finite mass required')
            if not angular.is_finite():
                raise ValueError('finite signed angular label required')
            if not hint.is_finite() or not radius.is_finite() or not radius > 0:
                raise ValueError('finite horizon hint and positive analytic radius required')
            frame = metric_horizon_frame(
                energy, mass, angular, hint, order=frame_order,
                analytic_radius=radius, bits=bits)
            self._frame = frame
            self.q = frame.q
            self.energy = frame.energy
            self.mass = mass
            self.angular = angular
            self.bits = bits
            self.metric_terms = metric_terms
            self.numeric_coefficients = np.array(
                [float(term) for term in _coefficients(self.q, metric_terms)])

    def hamiltonian_series(self, y, order):
        return BlochSource.hamiltonian_series(self, y, order)

    def rhs_numeric(self, y, n):
        # Homogeneous cross product only. Do not subtract VacuumCorrection forcing.
        return BlochSource.rhs_numeric(self, y, n)

    def vacuum_initial(self, y_start):
        """Enclosed Bloch vector of interior frame column 0, tail included."""
        with ctx.workprec(self.bits):
            y = _explicit_real(y_start, 'start')
            if not y.is_finite():
                raise ValueError('finite start log-distance required')
            columns, _tail = self._frame.evaluate(y.exp(), interior=True)
            top, bottom = columns[0][0], columns[1][0]
            coherence = top*bottom.conjugate()
            initial = (
                2*coherence.real,
                -2*coherence.imag,
                (top*top.conjugate()).real-(bottom*bottom.conjugate()).real)
            if any(not component.is_finite() for component in initial):
                raise ArithmeticError('nonfinite vacuum projector')
            return initial

    def target_log_distance(self, rho):
        """Log compact distance of rho in the same horizon chart as the generator."""
        with ctx.workprec(self.bits):
            rho = _explicit_real(rho, 'target rho')
            if not rho.is_finite() or not rho >= 1:
                raise ValueError('finite target rho>=1 required')
            # metric_H is proved only for |u|<0.4. The path is monotone in y.
            delta = self.q-arb.pi()/2-rho.atan()
            if not delta > 0 or not delta < arb('0.4'):
                raise ValueError('target rho is outside the homogeneous metric collar')
            return delta.log()


def capture_direct_vacuum(model, y_start, rho, *, max_step=None, cpu_limit=None):
    """Integrate the homogeneous vacuum. An explicit CPU budget is required.

    Timeout raises and does not return a trajectory. A failed solve is not
    repackaged as a complete witness.
    """
    if type(model) is not DirectVacuumBloch:
        raise TypeError('direct homogeneous vacuum model required')
    step = _positive_control(max_step, 'max step')
    budget = _positive_control(cpu_limit, 'CPU budget')
    initial = model.vacuum_initial(y_start)
    target = model.target_log_distance(rho)
    start = _declared_float(y_start, 'start')
    end = float(target.mid())
    if not np.isfinite(end) or not end > start:
        raise ValueError('target rho must lie outward of the finite start')
    numeric = np.array([float(component.mid()) for component in initial])
    started = time.process_time()
    timed_out = False

    def rhs(y, state):
        nonlocal timed_out
        if time.process_time()-started > budget:
            timed_out = True
            raise TimeoutError('direct vacuum Bloch CPU budget exhausted')
        return model.rhs_numeric(y, state)

    try:
        run = solve_ivp(
            rhs, (start, end), numeric, method='DOP853', rtol=2e-13, atol=2e-15,
            max_step=step, dense_output=True)
    except TimeoutError:
        raise
    if timed_out:
        raise TimeoutError('direct vacuum Bloch CPU budget exhausted')
    if not getattr(run, 'success', False) or getattr(run, 'sol', None) is None:
        raise ArithmeticError(getattr(run, 'message', 'direct vacuum capture failed'))
    rows = [
        np.r_[run.t[j], run.t[j+1], run.y[:, j], run.y[:, j+1], polynomial.F[1:].ravel()]
        for j, polynomial in enumerate(run.sol.interpolants)]
    return np.asarray(rows), int(run.nfev)


def validate_direct_vacuum(model, trace, y_start, rho, *, degree=12):
    """Replay a saved trace against the enclosed column-0 initial and target rho.

    The ball-to-float initial discrepancy is retained. Joins and the binary
    start time are required; the metric bridge uses the actual target rho.
    """
    if type(model) is not DirectVacuumBloch:
        raise TypeError('direct homogeneous vacuum model required')
    if trace is None or isinstance(trace, (bool, np.bool_)):
        raise ValueError('explicit Bloch trajectory required')
    if isinstance(trace, np.ndarray) and np.issubdtype(trace.dtype, np.bool_):
        raise ValueError('explicit Bloch trajectory required')
    if np.iscomplexobj(trace):
        raise ValueError('real Bloch trajectory required')
    trace = np.asarray(trace, float)
    if (trace.ndim != 2 or trace.shape[1] != 26 or len(trace) == 0
            or not np.isfinite(trace).all()):
        raise ValueError('complete finite Bloch trajectory required')
    if trace[0, 0] != _declared_float(y_start, 'start'):
        raise ValueError('Bloch trajectory must start at its declared start time')
    initial = model.vacuum_initial(y_start)
    target = model.target_log_distance(rho)
    return validate_bloch(model, initial, trace, target, degree=degree)


def _endpoint(validation):
    if not isinstance(validation, dict) or 'endpoint' not in validation or 'bits' not in validation:
        raise ValueError('validated Bloch endpoint required')
    if type(validation['bits']) is not int:
        raise ValueError('validated precision required')
    error = _explicit_real(validation.get('bloch_error'), 'Bloch error')
    if not error.is_finite() or not error >= 0:
        raise ValueError('nonnegative finite Bloch error required')
    endpoint = validation['endpoint']
    if (endpoint is None or isinstance(endpoint, (bool, np.bool_)) or len(endpoint) != 3
            or any(component is None or isinstance(component, (bool, np.bool_))
                   for component in endpoint)):
        raise ValueError('complete three-component Bloch endpoint required')
    return endpoint


def archived_vacuum_distance(columns, covariance, validation, *, complement=False):
    """Two-sided archived distance. Occupation is not added.

    complement=True compares the actual supplied columns with (nx, -ny, -nz),
    the Bloch image of I-S3 conj(Q) S3. It does not reuse another row's bound.
    """
    if type(complement) is not bool:
        raise ValueError('explicit complement choice required')
    columns = _archived_matrix(columns, (2, 3), 'columns')
    covariance = _archived_matrix(covariance, (3, 3), 'covariance')
    endpoint = _endpoint(validation)
    used = validation
    if complement:
        used = dict(validation)
        # Q_- = I - S3 conj(Q_+) S3. Do not replace this by -n or by the positive bound.
        used['endpoint'] = (endpoint[0], -endpoint[1], -endpoint[2])
    return original_covariance_distance_bounds(columns, covariance, used)


def compare_signed_archives(model, validation, positive_columns, positive_covariance,
                            negative_columns, negative_covariance, *, kappa, omega):
    """Positive and negative archived distances, plus a separate occupation allowance.

    The occupation term is the reviewed max(sqrt(f), n). It is not added here
    and is not a substitute for either archived matrix.
    """
    if type(model) is not DirectVacuumBloch:
        raise TypeError('direct homogeneous vacuum model required')
    if (not isinstance(validation, dict) or validation.get('bits') != model.bits):
        raise ValueError('validation precision must match the vacuum model')
    for value, name in ((kappa, 'kappa'), (omega, 'omega')):
        number = _explicit_real(value, name)
        if not number.is_finite() or not number > 0:
            raise ValueError('positive '+name+' required')
    positive = archived_vacuum_distance(
        positive_columns, positive_covariance, validation, complement=False)
    negative = archived_vacuum_distance(
        negative_columns, negative_covariance, validation, complement=True)
    occupation = occupation_vacuum_distance(
        model.energy, model.mass, kappa, omega, bits=model.bits)
    return {
        'archived_positive': positive,
        'archived_negative': negative,
        'occupation': occupation,
    }
