"""Damped coupled Newton CONTROL for local incoming (w,U) histories.

This is a numerical control interface, not a physical gate. It accepts an
evaluator that returns the KS constraint schema (raw N,beta and the full
retarded history Jacobian) for each LocalIncomingFamily candidate. Geometry
Jacobians may precondition the linear step; they never replace the retarded
state response. Every result remains OPEN. A stop is not NON-EXISTENCE.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from inspect import Parameter, signature
from types import MappingProxyType

import numpy as np

from .nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response,
)
from .nsc_ks_source_envelope import physical_incoming_interval
from .nsc_local_incoming_constraints import (
    ERROR_BUDGET_COMPONENTS, STATIONARITY_TOLERANCE, local_error_budget,
)
from .nsc_local_incoming_family import LocalIncomingFamily
from .nsc_ks_profile_identity import profile_identity


ALLOWED_COEFFICIENT_COUNTS = (8, 16, 32)
NUMERICAL_CONTROL_TOLERANCE = STATIONARITY_TOLERANCE
PRODUCTION_LOCAL_GATE = 'I=S(1)+[.12,.18]'
FIXED_AXIAL_INNER = 0.03
FIXED_AXIAL_OUTER = 0.06
FIXED_NORMAL_INNER = 0.007
FIXED_NORMAL_OUTER = 0.03
STOP_CONVERGED = 'converged_numerical_control'
STOP_LINE_SEARCH = 'failed_line_search'
STOP_RANK = 'rank_or_conditioning'
STOP_BUDGET = 'iteration_budget'
STOP_VERIFICATION = 'resolution_or_verification'
STOP_REASON_PHRASES = MappingProxyType({
    STOP_CONVERGED: 'converged numerical control',
    STOP_LINE_SEARCH: 'failed line search',
    STOP_RANK: 'rank/conditioning',
    STOP_BUDGET: 'iteration budget',
    STOP_VERIFICATION: 'resolution/verification',
})
_FROZEN_DERIVATIVE_FLAGS = (
    'constant_frozen_matter_used_as_physical_state',
    'reference_state_tangent_subtracted',
    'incoming_state_derivative_removed',
    'frozen_incoming_state_derivative',
)
_FORBIDDEN_STATE = frozenset({
    'incoming_C0', 'C0', 'stress_approximant', 'parts', 'node_fields',
    'HistoryBinding', 'frozen_matter', 'matter_dict',
})
_RECTANGULAR_JUSTIFICATION = (
    'The local geometric operators have mixed differential orders '
    '(N ~ w\'\'/U, beta ~ w\'\'\'/U\'), so a unique square collocation solve '
    'is not assumed. The tau layout uses n N-equations, n-3 beta-equations '
    'and 3 declared jet conditions for 2n coefficients. The actual residual '
    'also contains retarded state response, which is not a local ODE of those '
    'orders; dropping beta nodes discards genuine Jacobian information. '
    'Rectangular Gauss-Newton on all n N and n beta nodes plus the same 3 '
    'declared conditions is therefore the better production linearization. '
    'Neither layout certifies the residual between nodes.'
)


def _finite_real(value, name):
    if value is None:
        raise ValueError('explicit ' + name + ' required')
    raw = np.asarray(value)
    if np.iscomplexobj(raw) and np.any(raw.imag != 0):
        raise ValueError('real ' + name + ' required')
    result = np.asarray(raw.real if np.iscomplexobj(raw) else raw, dtype=float)
    if not np.isfinite(result).all():
        raise ValueError('finite ' + name + ' required')
    return result


def _copy_real(value, name):
    result = np.array(_finite_real(value, name), dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _positive_int(value, name, allowed=None):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(name + ' must be a positive integer')
    value = int(value)
    if allowed is not None:
        if value not in allowed:
            raise ValueError(name + ' must be 8, 16 or 32')
    elif value < 1:
        raise ValueError(name + ' must be a positive integer')
    return value


def _require_family(family):
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError('LocalIncomingFamily required; frozen C0 and matter dictionaries are not this control')
    n = len(family.coefficients[0])
    if n not in ALLOWED_COEFFICIENT_COUNTS:
        raise ValueError('coefficient count must be 8, 16 or 32; resolution is not increased automatically')
    expected_interval = physical_incoming_interval()
    if not np.allclose(family.interval, expected_interval, atol=3e-12, rtol=0):
        raise ValueError('incoming interval is the fixed I=S(1)+[.12,.18]')
    if not np.allclose(
            [family.axial_inner, family.axial_outer, family.normal_inner, family.normal_outer],
            [FIXED_AXIAL_INNER, FIXED_AXIAL_OUTER, FIXED_NORMAL_INNER, FIXED_NORMAL_OUTER],
            atol=3e-15, rtol=0):
        raise ValueError('axial cutoff and normal support are fixed solver-domain data')
    return family, n


def _family_key(family):
    coeff = np.asarray(family.coefficients, float)
    return (
        tuple(coeff[0].tolist()), tuple(coeff[1].tolist()),
        float(family.center), float(family.axial_inner), float(family.axial_outer),
        float(family.normal_inner), float(family.normal_outer),
    )


def _nodes_match(actual, expected):
    actual = np.asarray(actual, float)
    expected = np.asarray(expected, float)
    return actual.shape == expected.shape and np.array_equal(actual, expected)


def _snapshot_identity(value):
    """Immutable snapshot. Mutable dict/list aliases cannot mask later mutation."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (str, bytes)):
        return value
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        if not np.isfinite(number):
            raise ValueError('finite identity values required')
        return number
    if value is None:
        return None
    if isinstance(value, Mapping):
        return tuple((str(key), _snapshot_identity(value[key])) for key in sorted(value, key=str))
    if isinstance(value, (list, tuple)):
        return tuple(_snapshot_identity(item) for item in value)
    if isinstance(value, np.ndarray):
        array = np.array(value, copy=True)
        return ('ndarray', array.shape, str(array.dtype), array.tobytes())
    raise TypeError(
        'source/history identity must be a snapshotable value, not '
        + type(value).__name__)


def _error_budget(raw=None):
    if raw is None:
        return local_error_budget()
    if not isinstance(raw, Mapping):
        raise TypeError('error_budget mapping required')
    return local_error_budget(**{
        name: raw[name] if name in raw else None for name in ERROR_BUDGET_COMPONENTS
    })


def _batch_digest_bindings(records, label):
    if not isinstance(records, Mapping):
        raise TypeError(label + ' must be a mapping of production batch bindings')
    if not records:
        raise ValueError(label + ' are empty; every production batch must bind source and preparation')
    missing = []
    items = []
    for key in sorted(records, key=str):
        record = records[key]
        if not isinstance(record, Mapping):
            missing.append(str(key))
            continue
        source = record.get('source_digest')
        preparation = record.get('preparation_digest')
        if source is None or preparation is None:
            missing.append(str(key))
            continue
        items.append((str(key), str(source), str(preparation)))
    if missing:
        raise ValueError(
            'partial missing batch bindings: source_digest and preparation_digest '
            'required for every ' + label + ' (' + ', '.join(missing) + ')')
    return tuple(items)


def _extract_source_identity(raw):
    parts = []
    if raw.get('source_identity') is not None:
        parts.append(('explicit', _snapshot_identity(raw['source_identity'])))
    if raw.get('batch_records') is not None:
        parts.append(('batch_records', _batch_digest_bindings(raw['batch_records'], 'batch_records')))
    elif raw.get('family_records') is not None:
        parts.append(('family_records', _batch_digest_bindings(raw['family_records'], 'family_records')))
    if not parts:
        raise ValueError('upstream/source identity required; a missing source cannot be assumed unchanged')
    if len(parts) == 1:
        return parts[0][1]
    return tuple(parts)


def _history_identity_matches(supplied, expected, family):
    if isinstance(supplied,str):
        return supplied==profile_identity(family)
    try:
        if _snapshot_identity(supplied) == _snapshot_identity(expected):
            return True
    except (TypeError, ValueError):
        pass
    return False


def _declared_history_identity(raw, family):
    supplied = raw.get('history_identity')
    if supplied is None:
        supplied = raw.get('profile_identity')
    if supplied is None:
        raise ValueError(
            'explicit matching history identity required; geometry slots do not prove state evolution')
    expected = _family_key(family)
    if not _history_identity_matches(supplied, expected, family):
        raise ValueError('mismatched history identity for the supplied LocalIncomingFamily')
    return expected


def _require_retarded_declaration(raw):
    scope = raw.get('scope') if isinstance(raw.get('scope'), Mapping) else {}
    for name in _FROZEN_DERIVATIVE_FLAGS:
        if raw.get(name) is True or scope.get(name) is True:
            raise ValueError(
                'scope admits frozen matter or a removed state derivative; '
                'full retarded state derivative is required')
    declared = raw.get('full_retarded_state_derivative')
    if declared is None:
        declared = scope.get('full_retarded_state_derivative')
    if declared is not True:
        raise ValueError(
            'explicit full_retarded_state_derivative declaration required; '
            'this is a production-owner contract, not a physics certificate')


def _open_scope(interval, extra=None):
    scope = {
        'interval': list(interval),
        'production_local_gate': PRODUCTION_LOCAL_GATE,
        'global_matching': False,
        'extended_stationarity': 'OPEN',
        'metric_timestep': False,
        'observations': False,
        'physical_constraint_status': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'certified_constraint_residual': None,
        'full_source_error_bound': None,
        'between_node_remainder_certified': False,
        'constant_frozen_matter_used_as_physical_state': False,
        'optimization_failure_is_nonexistence': False,
        'fresh_source_evolution_required_for_physical_gate': True,
        'unresolved': [
            'baseline low-energy/subgap source accuracy',
            'changed-history finite-energy coverage',
            'changed-history tail',
            'preparation/field/axial-derivative error',
            'coefficient interval and arithmetic',
            'between-node remainder on I',
            'complete source/error certificate',
            'fresh source-fixed evolution of the production evaluator',
        ],
    }
    if extra:
        scope.update(extra)
    return MappingProxyType(scope)


@dataclass(frozen=True)
class DeclaredJetBoundary:
    """Caller-declared w, w', w'' at one point of I.

    These numbers are solver controls that close the mixed-order free data.
    They are not physical matching data and are not selected by the residual.
    """
    point: float
    w: float
    w_z: float
    w_zz: float

    def __post_init__(self):
        point, w, w_z, w_zz = (_finite_real(v, name).item()
                               for v, name in ((self.point, 'boundary point'),
                                               (self.w, 'declared w'),
                                               (self.w_z, "declared w'"),
                                               (self.w_zz, "declared w''")))
        object.__setattr__(self, 'point', point)
        object.__setattr__(self, 'w', w)
        object.__setattr__(self, 'w_z', w_z)
        object.__setattr__(self, 'w_zz', w_zz)

    def description(self):
        return {
            'point': self.point, 'w': self.w, 'w_z': self.w_z, 'w_zz': self.w_zz,
            'role': 'solver_control', 'physical_matching_data': False,
            'specification': "declared Cauchy jet (w, w', w'') at a point of I",
        }


@dataclass(frozen=True)
class HistoryCollocation:
    """Tau or rectangular collocation on the fixed incoming interval.

    Default n=8. n=16 or 32 is an explicit caller choice. The solver never
    raises the resolution by itself. Verification nodes are a separate grid.
    """
    coefficient_count: int = 8
    layout: str = 'tau'

    def __post_init__(self):
        count = _positive_int(self.coefficient_count, 'coefficient_count', ALLOWED_COEFFICIENT_COUNTS)
        layout = str(self.layout)
        if layout not in ('tau', 'rectangular'):
            raise ValueError("collocation layout must be 'tau' or 'rectangular'")
        object.__setattr__(self, 'coefficient_count', count)
        object.__setattr__(self, 'layout', layout)

    def solve_nodes(self, family):
        family, n = _require_family(family)
        if n != self.coefficient_count:
            raise ValueError('collocation coefficient_count must match the family')
        return np.array(family.collocation_nodes(n), dtype=float, copy=True)

    def beta_indices(self):
        n = self.coefficient_count
        if self.layout == 'rectangular':
            return np.arange(n)
        # Tau: n-3 beta nodes. Drop the left endpoint and the two rightmost
        # Chebyshev-Lobatto nodes; the 3 jet conditions replace that data.
        return np.arange(1, n - 2)

    def n_equation_counts(self):
        n = self.coefficient_count
        n_beta = n if self.layout == 'rectangular' else n - 3
        return {'N': n, 'beta': n_beta, 'boundary': 3, 'unknowns': 2 * n}

    def verification_nodes(self, family, count=None):
        family, n = _require_family(family)
        if n != self.coefficient_count:
            raise ValueError('collocation coefficient_count must match the family')
        count = 2 * n + 1 if count is None else _positive_int(count, 'verification node count')
        if count == n:
            raise ValueError('verification nodes must be separate from the solve nodes')
        return np.array(family.collocation_nodes(count), dtype=float, copy=True)

    def description(self, family=None):
        counts = self.n_equation_counts()
        record = {
            'coefficient_count': self.coefficient_count, 'layout': self.layout,
            'equations': counts,
            'tau_beta_rule': None if self.layout == 'rectangular' else (
                'beta collocation uses Chebyshev-Lobatto indices 1:n-2; '
                'drops the left endpoint and the two rightmost nodes'),
            'rectangular_justification': _RECTANGULAR_JUSTIFICATION,
            'automatic_resolution': False,
            'verification_nodes_separate': True,
        }
        if family is not None:
            record['solve_nodes'] = self.solve_nodes(family).tolist()
            record['beta_indices'] = self.beta_indices().tolist()
        return record


@dataclass(frozen=True)
class LocalHistoryNewtonSettings:
    """Numerical-control knobs. No physical source, seam or stress parameters."""

    residual_tolerance: float = NUMERICAL_CONTROL_TOLERANCE
    max_iterations: int = 20
    max_step_norm: float = 1.0
    min_line_search_scale: float = 1 / 64
    line_search_factor: float = 0.5
    armijo: float = 1e-4
    rank_tolerance: float = 1e-14
    max_condition_number: float = 1e16
    geometric_preconditioner: bool = True

    def __post_init__(self):
        residual_tolerance = float(_finite_real(self.residual_tolerance, 'residual_tolerance'))
        max_step_norm = float(_finite_real(self.max_step_norm, 'max_step_norm'))
        min_line_search_scale = float(_finite_real(self.min_line_search_scale, 'min_line_search_scale'))
        line_search_factor = float(_finite_real(self.line_search_factor, 'line_search_factor'))
        armijo = float(_finite_real(self.armijo, 'armijo'))
        rank_tolerance = float(_finite_real(self.rank_tolerance, 'rank_tolerance'))
        max_condition_number = float(_finite_real(self.max_condition_number, 'max_condition_number'))
        if residual_tolerance <= 0 or max_step_norm <= 0 or rank_tolerance <= 0 or max_condition_number <= 1:
            raise ValueError('positive ordered residual/step/rank/condition settings required')
        if not 0 < min_line_search_scale <= 1:
            raise ValueError('min_line_search_scale must be in (0, 1]')
        if not 0 < line_search_factor < 1:
            raise ValueError('line_search_factor must be in (0, 1)')
        if not 0 <= armijo < 1:
            raise ValueError('armijo must be in [0, 1)')
        object.__setattr__(self, 'residual_tolerance', residual_tolerance)
        object.__setattr__(self, 'max_step_norm', max_step_norm)
        object.__setattr__(self, 'min_line_search_scale', min_line_search_scale)
        object.__setattr__(self, 'line_search_factor', line_search_factor)
        object.__setattr__(self, 'armijo', armijo)
        object.__setattr__(self, 'rank_tolerance', rank_tolerance)
        object.__setattr__(self, 'max_condition_number', max_condition_number)
        iterations = self.max_iterations
        if isinstance(iterations, (bool, np.bool_)) or not isinstance(iterations, (int, np.integer)) or int(iterations) < 0:
            raise ValueError('max_iterations must be a nonnegative integer')
        object.__setattr__(self, 'max_iterations', int(iterations))
        object.__setattr__(self, 'geometric_preconditioner', bool(self.geometric_preconditioner))


@dataclass(frozen=True)
class HistoryEvaluatorReceipt:
    """Typed binding of one evaluator result to a candidate history.

    Matches the KSConstraintAccumulator residual/Jacobian schema. Missing
    error-budget components remain None. Physical status is OPEN.
    """
    family: LocalIncomingFamily
    z: object
    action_gradient: object
    history_jacobian: object
    source_identity: object
    history_identity: object
    incoming_slots: object
    incoming_slot_tangents: object
    local_reference_gradient_tangent: object
    error_budget: object
    full_retarded_state_derivative: bool = False
    constraint_order: object = ('N', 'beta')
    physical_constraint_status: str = 'OPEN'
    physical_EXISTENCE_certificate: bool = False
    certified_constraint_residual: object = None
    full_source_error_bound: object = None
    missing_error_budget: object = None

    def __post_init__(self):
        z = _copy_real(self.z, 'incoming z')
        if z.ndim != 1 or len(z) < 3:
            raise ValueError('incoming z must be a one-dimensional node set')
        gradient = _copy_real(self.action_gradient, 'action_gradient')
        jacobian = _copy_real(self.history_jacobian, 'history_jacobian')
        slots = _copy_real(self.incoming_slots, 'incoming_slots')
        slot_tangents = _copy_real(self.incoming_slot_tangents, 'incoming_slot_tangents')
        n = len(self.family.coefficients[0])
        ndir = 2 * n
        nz = len(z)
        if gradient.shape != (nz, 2):
            raise ValueError('raw N,beta action_gradient must have shape (nz, 2)')
        if jacobian.shape != (ndir, nz, 2):
            raise ValueError(
                'full retarded history_jacobian of shape (ndir, nz, 2) required; '
                'geometric polynomials are not a substitute')
        if slots.shape != (nz, 6) or slot_tangents.shape != (ndir, nz, 6):
            raise ValueError('(nz,6) slots and (ndir,nz,6) slot tangents required')
        geom = self.local_reference_gradient_tangent
        if geom is not None:
            geom = _copy_real(geom, 'local_reference_gradient_tangent')
            if geom.shape != (ndir, nz, 2):
                raise ValueError('geometric tangent, if supplied, must have shape (ndir, nz, 2)')
        object.__setattr__(self, 'z', z)
        object.__setattr__(self, 'action_gradient', gradient)
        object.__setattr__(self, 'history_jacobian', jacobian)
        object.__setattr__(self, 'incoming_slots', slots)
        object.__setattr__(self, 'incoming_slot_tangents', slot_tangents)
        object.__setattr__(self, 'local_reference_gradient_tangent', geom)
        object.__setattr__(self, 'source_identity', _snapshot_identity(self.source_identity))
        object.__setattr__(self, 'history_identity', _snapshot_identity(self.history_identity))
        if self.full_retarded_state_derivative is not True:
            raise ValueError('HistoryEvaluatorReceipt requires full_retarded_state_derivative=True')
        object.__setattr__(self, 'error_budget', MappingProxyType(dict(self.error_budget)))
        object.__setattr__(self, 'full_retarded_state_derivative', True)
        object.__setattr__(self, 'constraint_order', ('N', 'beta'))
        object.__setattr__(self, 'physical_constraint_status', 'OPEN')
        object.__setattr__(self, 'physical_EXISTENCE_certificate', False)
        object.__setattr__(self, 'certified_constraint_residual', None)
        object.__setattr__(self, 'full_source_error_bound', None)
        object.__setattr__(self, 'missing_error_budget', None)


def bind_history_evaluator_receipt(raw, family, *, expected_nodes=None,
                                   expected_source_identity=None):
    """Validate one evaluator mapping against the KS residual/Jacobian schema."""
    family, n = _require_family(family)
    if isinstance(raw, HistoryEvaluatorReceipt):
        raw = {
            'z': raw.z, 'action_gradient': raw.action_gradient,
            'history_jacobian': raw.history_jacobian,
            'source_identity': raw.source_identity,
            'history_identity': raw.history_identity,
            'incoming_slots': raw.incoming_slots,
            'incoming_slot_tangents': raw.incoming_slot_tangents,
            'local_reference_gradient_tangent': raw.local_reference_gradient_tangent,
            'error_budget': raw.error_budget,
            'full_retarded_state_derivative': raw.full_retarded_state_derivative,
            'constraint_order': raw.constraint_order,
            'physical_constraint_status': raw.physical_constraint_status,
            'physical_EXISTENCE_certificate': raw.physical_EXISTENCE_certificate,
            'constant_frozen_matter_used_as_physical_state': False,
            'reference_state_tangent_subtracted': False,
        }
    if not isinstance(raw, Mapping):
        raise TypeError('evaluator must return a KS-constraint-schema mapping')
    forbidden = _FORBIDDEN_STATE.intersection(raw)
    if forbidden:
        raise ValueError(
            'frozen C0/matter-dict shortcut is not a physical certificate: '
            + ', '.join(sorted(forbidden)))
    if raw.get('physical_EXISTENCE_certificate') is True:
        raise ValueError('this CONTROL owner does not accept a physical EXISTENCE certificate')
    status = raw.get('physical_constraint_status')
    if status not in (None, 'OPEN'):
        raise ValueError('physical constraint status remains OPEN in this CONTROL owner')
    order = raw.get('constraint_order', ('N', 'beta'))
    if tuple(order) != ('N', 'beta'):
        raise ValueError("constraint_order must be ('N', 'beta')")
    _require_retarded_declaration(raw)
    if 'history_jacobian' not in raw or raw['history_jacobian'] is None:
        raise ValueError('full retarded history_jacobian is required and must be consumed')
    z = _finite_real(raw['z'], 'incoming z')
    if expected_nodes is not None and not _nodes_match(z, expected_nodes):
        raise ValueError('mismatched incoming nodes: evaluator z must match the bound solver nodes')
    ndir = 2 * n
    slots, slot_tangents = compatible_history_slots(family.metric(), z, ndir)
    if 'incoming_slots' in raw and raw['incoming_slots'] is not None:
        supplied = _finite_real(raw['incoming_slots'], 'incoming_slots')
        if supplied.shape != slots.shape or np.max(np.abs(supplied - slots)) > 3e-11:
            raise ValueError('mismatched history: incoming slots do not belong to the supplied family')
    if 'incoming_slot_tangents' in raw and raw['incoming_slot_tangents'] is not None:
        supplied_t = _finite_real(raw['incoming_slot_tangents'], 'incoming_slot_tangents')
        if supplied_t.shape != slot_tangents.shape or np.max(np.abs(supplied_t - slot_tangents)) > 3e-11:
            raise ValueError('mismatched history: slot tangents do not belong to the supplied family')
    source_identity = _extract_source_identity(raw)
    if expected_source_identity is not None and source_identity != expected_source_identity:
        raise ValueError('upstream/source identity changed across history iterations')
    history_identity = _declared_history_identity(raw, family)
    geom = raw.get('local_reference_gradient_tangent')
    coefficients = raw.get('coefficients') or raw.get('surface_coefficients')
    if geom is None and coefficients is not None:
        geom = surface_geometry_response(slots, slot_tangents, coefficients)['action_gradient_tangent']
    elif geom is None and 'local_reference_gradient_change' in raw and coefficients is not None:
        geom = surface_geometry_response(slots, slot_tangents, coefficients)['action_gradient_tangent']
    return HistoryEvaluatorReceipt(
        family=family, z=z, action_gradient=raw['action_gradient'],
        history_jacobian=raw['history_jacobian'], source_identity=source_identity,
        history_identity=history_identity, incoming_slots=slots,
        incoming_slot_tangents=slot_tangents,
        local_reference_gradient_tangent=geom,
        error_budget=_error_budget(raw.get('error_budget')),
        full_retarded_state_derivative=True,
    )


def _invoke_evaluator(evaluator, family, z):
    if hasattr(evaluator, 'evaluate') and callable(evaluator.evaluate):
        fn = evaluator.evaluate
    elif callable(evaluator):
        fn = evaluator
    else:
        raise TypeError('evaluator protocol requires evaluate(family) with optional z')
    try:
        params = signature(fn).parameters
    except (TypeError, ValueError):
        params = {}
    requested = np.array(z, dtype=float, copy=True)
    if 'z' in params:
        return fn(family, z=requested)
    positional = [
        item for item in params.values()
        if item.kind in (Parameter.POSITIONAL_ONLY, Parameter.POSITIONAL_OR_KEYWORD)
    ]
    if len(positional) >= 2:
        return fn(family, requested)
    raise TypeError('evaluator must accept the requested node array z')


class GeometryBoundEvaluatorCache:
    """Cache immutable receipts by geometry identity and node set only."""

    def __init__(self, evaluator):
        self.evaluator = evaluator
        self._store = {}

    def evaluate(self, family, z, expected_source_identity=None):
        requested = np.array(z, dtype=float, copy=True)
        key = (_family_key(family), tuple(requested.tolist()))
        if key not in self._store:
            raw = _invoke_evaluator(self.evaluator, family, requested)
            self._store[key] = bind_history_evaluator_receipt(
                raw, family, expected_nodes=requested,
                expected_source_identity=expected_source_identity)
        receipt = self._store[key]
        if expected_source_identity is not None and receipt.source_identity != expected_source_identity:
            raise ValueError('upstream/source identity changed across history iterations')
        return receipt


def _replace_family(family, coefficients):
    family, n = _require_family(family)
    coeff = np.asarray(coefficients, float)
    if coeff.shape != (2, n):
        coeff = coeff.reshape(2, n)
    return LocalIncomingFamily(
        coeff, center=family.center, axial_inner=family.axial_inner,
        axial_outer=family.axial_outer, normal_inner=family.normal_inner,
        normal_outer=family.normal_outer,
    )


def _boundary_system(family, boundary):
    family, n = _require_family(family)
    interval = family.interval
    if not (interval[0] - 3e-12 <= boundary.point <= interval[1] + 3e-12):
        raise ValueError('declared boundary point must lie in I')
    z = np.array([boundary.point], dtype=float)
    slots, tangents = compatible_history_slots(family.metric(), z, 2 * n)
    value = slots[0, :3]
    declared = np.array([boundary.w, boundary.w_z, boundary.w_zz], dtype=float)
    inner = float(family.axial_inner)
    scale = np.array([1.0, inner, inner * inner], dtype=float)
    residual = (value - declared) * scale
    jacobian = tangents[:, 0, :3].T * scale[:, None]
    w_fn, _ = family.functions
    from_functions = np.array([w_fn(boundary.point, k) for k in range(3)], dtype=float)
    if np.max(np.abs(from_functions - value)) > 3e-11:
        raise ArithmeticError('boundary jet is not the derivative of the same w function')
    return residual, jacobian, scale


def _assemble_rows(block, indices, component):
    """block is (ndir, nz, 2) -> (len(indices), ndir) for one N/beta component."""
    return np.asarray(block, float)[:, indices, component].T


def _control_system(receipt, collocation, boundary):
    n = collocation.coefficient_count
    if len(receipt.z) != n:
        raise ValueError('solve-node evaluator must return exactly n collocation nodes')
    n_index = np.arange(n)
    beta_index = collocation.beta_indices()
    r_n = np.asarray(receipt.action_gradient, float)[n_index, 0]
    r_beta = np.asarray(receipt.action_gradient, float)[beta_index, 1]
    j_n = _assemble_rows(receipt.history_jacobian, n_index, 0)
    j_beta = _assemble_rows(receipt.history_jacobian, beta_index, 1)
    r_bc, j_bc, bc_scale = _boundary_system(receipt.family, boundary)
    residual = np.concatenate([r_n, r_beta, r_bc])
    jacobian = np.vstack([j_n, j_beta, j_bc])
    geom = receipt.local_reference_gradient_tangent
    if geom is None:
        tangents = receipt.incoming_slot_tangents
        geom_n = (tangents[:, n_index, 4] + tangents[:, n_index, 2]).T
        geom_beta = (tangents[:, beta_index, 5] + tangents[:, beta_index, 3]).T
        geom_matrix = np.vstack([geom_n, geom_beta, j_bc])
        geom_source = 'structural_order_map'
        preconditioner = None
    else:
        geom_matrix = np.vstack([
            _assemble_rows(geom, n_index, 0),
            _assemble_rows(geom, beta_index, 1),
            j_bc,
        ])
        geom_source = 'local_reference_gradient_tangent'
        preconditioner = geom_matrix
    discrepancy = float(np.max(np.abs(jacobian - geom_matrix)))
    return {
        'residual': residual, 'jacobian': jacobian, 'preconditioner': preconditioner,
        'n_residual': r_n, 'beta_residual': r_beta, 'boundary_residual': r_bc,
        'beta_indices': beta_index, 'boundary_scale': bc_scale,
        'geometric_source': geom_source, 'full_minus_geometric_max': discrepancy,
        'newton_jacobian_source': 'full_retarded_history_jacobian',
    }


def _linear_step(jacobian, residual, preconditioner, settings):
    J = np.asarray(jacobian, float)
    r = np.asarray(residual, float)
    n_eq, n_unknowns = J.shape
    use_geom = bool(settings.geometric_preconditioner and preconditioner is not None)
    P = np.asarray(preconditioner, float) if use_geom else None
    if P is not None and P.shape != J.shape:
        P = None
        use_geom = False
    row = np.linalg.norm(J, axis=1)
    col = np.linalg.norm(J, axis=0)
    if P is not None:
        col_p = np.linalg.norm(P, axis=0)
        # Geometric column scales are a preconditioner only when they stay
        # finite and do not erase a retarded column.
        if np.all(col_p > 0) and np.all(np.isfinite(col_p)):
            col = np.where(col > 0, np.sqrt(col * col_p), col_p)
    row = np.where(row > 0, row, 1.0)
    col = np.where(col > 0, col, 1.0)
    Js = J / row[:, None] / col[None, :]
    rs = r / row
    u, singular, vt = np.linalg.svd(Js, full_matrices=False)
    smax = float(singular[0]) if singular.size else 0.0
    machine = np.finfo(float).eps * max(n_eq, n_unknowns, 1) * max(smax, 1.0)
    floor = max(settings.rank_tolerance * max(smax, 1.0), machine)
    keep = singular > floor
    rank = int(np.count_nonzero(keep))
    smin = float(singular[-1]) if singular.size else 0.0
    condition = (smax / smin) if smin > 0 else float('inf')
    deficient = rank < n_unknowns or smax == 0.0 or not np.isfinite(smax)
    if np.isfinite(condition) and condition > settings.max_condition_number and rank == n_unknowns:
        deficient = True
    diagnostics = {
        'rank': rank, 'n_unknowns': n_unknowns, 'n_equations': n_eq,
        'singular_values': np.array(singular, copy=True),
        'condition_number': float(condition),
        'row_scales': row, 'column_scales': col,
        'rank_or_conditioning': deficient,
        'usable_range': rank, 'geometric_preconditioner': use_geom,
    }
    if deficient:
        return None, diagnostics
    delta_scaled = vt[keep].T @ ((u[:, keep].T @ (-rs)) / singular[keep])
    delta = delta_scaled / col
    if not np.isfinite(delta).all():
        diagnostics['rank_or_conditioning'] = True
        return None, diagnostics
    return delta, diagnostics


def _clip_step(delta, max_norm):
    delta = np.asarray(delta, float)
    norm = float(np.linalg.norm(delta))
    if norm > max_norm > 0:
        delta = delta * (max_norm / norm)
        norm = float(max_norm)
    return delta, norm


def _control_max(system):
    return float(np.max(np.abs(system['residual']))) if system['residual'].size else 0.0


def _evaluate_receipt(cache, family, z, expected_source_identity):
    return cache.evaluate(family, z, expected_source_identity=expected_source_identity)


def _step_record(iteration, scale, step_norm, radius, residual_max, accepted, reason, extra=None):
    record = {
        'iteration': iteration, 'step_scale': scale, 'step_norm': step_norm,
        'radius_lower_bound': radius, 'control_residual_max': residual_max,
        'accepted': accepted, 'reason': reason,
        'physical_constraint_status': 'OPEN',
        'nonexistence': False,
    }
    if extra:
        record.update(extra)
    return record


def solve_local_history(
        evaluator, initial_family, *, boundary, collocation=None, settings=None,
        verification_nodes=None, expected_source_identity=None):
    """Damped Newton/Gauss-Newton CONTROL on a LocalIncomingFamily.

    The evaluator must supply raw N,beta and the full retarded Jacobian for
    every candidate. Stop reasons are numerical control only. Physical gate
    remains OPEN regardless of the nodal residual.
    """
    if boundary is None:
        raise TypeError(
            "declared w, w', w'' data at a point of I (or an equally explicit "
            'boundary specification) are required; a unique square collocation '
            'solve is not assumed')
    if not isinstance(boundary, DeclaredJetBoundary):
        raise TypeError('DeclaredJetBoundary required')
    family, n = _require_family(initial_family)
    if collocation is None:
        collocation = HistoryCollocation(coefficient_count=n)
    if not isinstance(collocation, HistoryCollocation):
        raise TypeError('HistoryCollocation required')
    if collocation.coefficient_count != n:
        raise ValueError('collocation coefficient_count must match the initial family; n is not raised automatically')
    settings = LocalHistoryNewtonSettings() if settings is None else settings
    if not isinstance(settings, LocalHistoryNewtonSettings):
        raise TypeError('LocalHistoryNewtonSettings required')
    if family.radius_lower_bound() <= 0:
        raise ValueError('initial family radius_lower_bound must be positive before evaluator call')
    solve_nodes = _copy_real(collocation.solve_nodes(family), 'solve nodes')
    if verification_nodes is None:
        verification_nodes = _copy_real(collocation.verification_nodes(family), 'verification nodes')
    else:
        verification_nodes = _copy_real(verification_nodes, 'verification nodes')
        if _nodes_match(verification_nodes, solve_nodes):
            raise ValueError('verification nodes must be separate from the solve nodes')
    cache = GeometryBoundEvaluatorCache(evaluator)
    receipt = _evaluate_receipt(cache, family, solve_nodes, expected_source_identity)
    source_identity = receipt.source_identity
    system = _control_system(receipt, collocation, boundary)
    accepted_steps = []
    rejected_steps = []
    linear = None
    stop_reason = None
    x = np.asarray(family.coefficients, float).ravel().copy()

    def snapshot(reason, iterations):
        verify = _evaluate_receipt(cache, family, verification_nodes, source_identity)
        verify_max = np.max(np.abs(verify.action_gradient), axis=0)
        solve_max = np.max(np.abs(receipt.action_gradient), axis=0)
        control_converged = _control_max(system) <= settings.residual_tolerance
        solve_all_nodes_pass = float(np.max(np.abs(receipt.action_gradient))) <= settings.residual_tolerance
        verification_nodes_pass = float(np.max(np.abs(verify.action_gradient))) <= settings.residual_tolerance
        requested_nodes_pass = control_converged and solve_all_nodes_pass and verification_nodes_pass
        if reason == STOP_CONVERGED and not requested_nodes_pass:
            reason = STOP_VERIFICATION
        numerical_control_converged = reason == STOP_CONVERGED and requested_nodes_pass
        error_budget = _error_budget(receipt.error_budget)
        missing = [name for name in ERROR_BUDGET_COMPONENTS if error_budget[name] is None]
        geom_used = bool(
            settings.geometric_preconditioner and system['preconditioner'] is not None)
        return {
            'constraint_order': ('N', 'beta'),
            'family': family,
            'coefficients': np.array(family.coefficients, dtype=float, copy=True),
            'z': np.array(solve_nodes, copy=True),
            'action_gradient': np.array(receipt.action_gradient, copy=True),
            'history_jacobian': np.array(receipt.history_jacobian, copy=True),
            'incoming_slots': np.array(receipt.incoming_slots, copy=True),
            'incoming_slot_tangents': np.array(receipt.incoming_slot_tangents, copy=True),
            'local_reference_gradient_tangent': (
                None if receipt.local_reference_gradient_tangent is None
                else np.array(receipt.local_reference_gradient_tangent, copy=True)),
            'control_residual': np.array(system['residual'], copy=True),
            'control_residual_max': _control_max(system),
            'solve_node_maximum_absolute_residual': solve_max,
            'boundary_residual': np.array(system['boundary_residual'], copy=True),
            'boundary_control': boundary.description(),
            'collocation': collocation.description(family),
            'newton_jacobian_source': system['newton_jacobian_source'],
            'geometric_jacobian_used_as_newton_matrix': False,
            'geometric_preconditioner': geom_used,
            'geometric_preconditioner_source': system['geometric_source'],
            'full_minus_geometric_jacobian_max': system['full_minus_geometric_max'],
            'full_retarded_jacobian_consumed': True,
            'linear_algebra': None if linear is None else {
                'rank': linear['rank'], 'n_unknowns': linear['n_unknowns'],
                'n_equations': linear['n_equations'],
                'condition_number': linear['condition_number'],
                'rank_or_conditioning': linear['rank_or_conditioning'],
                'column_scales': np.array(linear['column_scales'], copy=True),
                'row_scales': np.array(linear['row_scales'], copy=True),
                'geometric_preconditioner': linear['geometric_preconditioner'],
            },
            'accepted_steps': tuple(accepted_steps),
            'rejected_steps': tuple(rejected_steps),
            'iterations': iterations,
            'stop_reason': reason,
            'stop_reason_phrase': STOP_REASON_PHRASES[reason],
            'control_converged': control_converged,
            'solve_all_nodes_pass': solve_all_nodes_pass,
            'verification_nodes_pass': verification_nodes_pass,
            'numerical_control_converged': numerical_control_converged,
            'residual_tolerance': settings.residual_tolerance,
            'stationarity_tolerance': NUMERICAL_CONTROL_TOLERANCE,
            'source_identity': source_identity,
            'history_identity': receipt.history_identity,
            'radius_lower_bound': family.radius_lower_bound(),
            'verification': {
                'z': np.array(verification_nodes, copy=True),
                'action_gradient': np.array(verify.action_gradient, copy=True),
                'maximum_absolute_residual': verify_max,
                'nodes_separate_from_solve_nodes': True,
                'verification_nodes_pass': verification_nodes_pass,
                'between_node_certified': False,
                'missing_error_certified': False,
                'automatic_resolution_increase': False,
            },
            'error_budget': error_budget,
            'missing_error_budget': None,
            'missing_error_budget_components': tuple(missing),
            'interval': family.interval,
            'scope': _open_scope(family.interval, {
                'stop_reason': reason,
                'control_converged': control_converged,
                'solve_all_nodes_pass': solve_all_nodes_pass,
                'verification_nodes_pass': verification_nodes_pass,
                'numerical_control_converged': numerical_control_converged,
                'collocation_layout': collocation.layout,
                'declared_boundary_is_solver_control': True,
                'full_retarded_state_derivative_declared': True,
            }),
            'physical_constraint_status': 'OPEN',
            'physical_EXISTENCE_certificate': False,
            'physical_NONEXISTENCE_certificate': False,
            'certified_constraint_residual': None,
            'full_source_error_bound': None,
            'optimization_failure_is_nonexistence': False,
            'stop_reason_is_nonexistence': False,
            'fresh_source_evolution_required_for_physical_gate': True,
        }

    if _control_max(system) <= settings.residual_tolerance:
        return snapshot(STOP_CONVERGED, 0)

    for iteration in range(1, settings.max_iterations + 1):
        delta, linear = _linear_step(
            system['jacobian'], system['residual'], system['preconditioner'], settings)
        if delta is None:
            return snapshot(STOP_RANK, iteration - 1)
        delta, step_norm = _clip_step(delta, settings.max_step_norm)
        scale = 1.0
        accepted = False
        current_norm = float(np.linalg.norm(system['residual']))
        while scale >= settings.min_line_search_scale:
            trial_x = x + scale * delta
            try:
                trial_family = _replace_family(family, trial_x)
            except ValueError as error:
                rejected_steps.append(_step_record(
                    iteration, scale, step_norm * scale, None, None, False,
                    'invalid_trial_family', {'detail': str(error)}))
                scale *= settings.line_search_factor
                continue
            radius = trial_family.radius_lower_bound()
            if radius <= 0:
                rejected_steps.append(_step_record(
                    iteration, scale, step_norm * scale, radius, None, False,
                    'nonpositive_radius'))
                scale *= settings.line_search_factor
                continue
            trial_receipt = _evaluate_receipt(cache, trial_family, solve_nodes, source_identity)
            trial_system = _control_system(trial_receipt, collocation, boundary)
            trial_max = _control_max(trial_system)
            trial_norm = float(np.linalg.norm(trial_system['residual']))
            if trial_norm <= current_norm * (1.0 - settings.armijo * scale):
                accepted_steps.append(_step_record(
                    iteration, scale, step_norm * scale, radius, trial_max, True,
                    'sufficient_decrease'))
                family = trial_family
                x = trial_x
                receipt = trial_receipt
                system = trial_system
                accepted = True
                break
            rejected_steps.append(_step_record(
                iteration, scale, step_norm * scale, radius, trial_max, False,
                'residual_increase'))
            scale *= settings.line_search_factor
        if not accepted:
            return snapshot(STOP_LINE_SEARCH, iteration)
        if _control_max(system) <= settings.residual_tolerance:
            return snapshot(STOP_CONVERGED, iteration)
    return snapshot(STOP_BUDGET, settings.max_iterations)
