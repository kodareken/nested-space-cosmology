"""Insert an upstream covariance error into the history-minus-reference matter difference.

For one complete same-energy fiber the homogeneous upstream matrix A is 2x3 and
Q=A C_src A^dagger. The same retarded solution acts on every column, so the
weighted fields are F=K A and F_z=L A. The right inverse
R=A^dagger (A A^dagger)^-1 gives an equivalent source insertion
delta C=R delta Q R^dagger. Original C_src and A are not modified.

If lambda_min(A A^dagger)>=gamma>0 and ||delta Q||_op<=epsilon, then
||delta C||_op<=epsilon/gamma. Energy-diagonal fibers use
eta=max_e epsilon_e/gamma_e. No column phases are aligned. A common phase
cancels in Q; relative phases in a coherent source remain part of Q.
Evolved momentum is never replaced by -E. Column weights already contain
sqrt(dE/(2 pi)) and are not applied again.

The difference error is Tr(delta Q (H[g]-H[ref])). Certified field errors are
added to captured weighted Frobenius uppers before that norm inequality, so a
source error on an uncertain field remains enclosed. The later inhomogeneous
kernel is not treated as pointwise unitary.

Directed Arb enclosures reject a nonpositive Gram eigenvalue, an interval that
straddles an invalid sign, a missing epsilon, and an incomplete fiber. Absent
bounds are rejected before any Arb conversion: arb(None) would silently become
zero. This owner does not load upstream records, does not copy a positive-energy
error onto the negative partner, and does not fill a physical budget. The full
upstream gate stays OPEN.
"""
from collections.abc import Mapping, Sequence

import numpy as np
from flint import acb, arb, ctx

from .nsc_evolved_incoming_constraints import raw_ks_vertex_coefficients, source_column_matter
from .nsc_ks_ball_trajectory import complex_ball, exact_upper, restored_upper
from .nsc_ks_batched_constraints import KSUpstreamBatch
from .nsc_ks_finite_matter_error import finite_matter_error


SOURCE_BOUND_SCOPE = 'caller-supplied covariance bounds; source evidence is authenticated by the caller'
ENDPOINT_FIELD_SLOTS = (
    'reference_norm', 'reference_axial_norm', 'difference_norm', 'difference_axial_norm',
    'reference_error', 'reference_axial_error', 'difference_error', 'difference_axial_error',
    'mass', 'absolute_angular', 'axial_lower', 'radius_lower', 'multiplicity',
)
_NEGATIVE_COPY = (
    'negative energy must use its own covariance bound; '
    'copying the positive bound is forbidden')


def _open_scope():
    return {
        'physical_source_coverage': False,
        'source_bound_scope': SOURCE_BOUND_SCOPE,
        'physical_upstream_budget_component': None,
        'full_physical_budget_populated': False,
        'physical_local_gate': 'OPEN',
        'full_upstream_gate': 'OPEN',
        'quadrature_applied_again': False,
        'evolved_momentum_replaced_by_minus_energy': False,
        'inhomogeneous_kernel_assumed_pointwise_unitary': False,
        'negative_covariance_bound_copied_from_positive': False,
        'gram_identity_assumed': False,
        'history_minus_reference_cancellation': True,
        'source_accuracy_included': False,
        'conditional_source_error_included': True,
        'upstream_records_loaded': False,
    }


def _missing(name):
    raise ValueError('missing ' + name)


def _as_ball(value, name):
    """Restore a directed bound. None is rejected before Arb conversion."""
    if value is None or isinstance(value, bool):
        _missing(name)
    if isinstance(value, Mapping):
        if 'mantissa' in value and 'exponent' in value:
            return restored_upper(value)
        if 'value' in value:
            if value['value'] is None:
                _missing(name)
            return _as_ball(value['value'], name)
        _missing(name)
    if isinstance(value, arb):
        return value
    if isinstance(value, (int, np.integer)):
        return arb(int(value))
    if isinstance(value, float):
        number = float(value)
        if not np.isfinite(number):
            raise ValueError('finite ' + name + ' required')
        return arb(number)
    if isinstance(value, np.floating):
        number = float(value)
        if not np.isfinite(number):
            raise ValueError('finite ' + name + ' required')
        return arb(number)
    _missing(name)


def _nonnegative_upper(value, name):
    number = _as_ball(value, name)
    if not number.is_finite() or not (number >= 0):
        raise ValueError('interval bound straddles an invalid sign: ' + name)
    return number.upper()


def _positive_magnitude(value, name):
    number = _as_ball(value, name)
    if not number.is_finite() or not (number > 0):
        raise ValueError('interval bound straddles an invalid sign: ' + name)
    return number.upper()


def _positive_lower(value, name):
    number = _as_ball(value, name)
    if not number.is_finite() or not (number > 0):
        raise ValueError('interval bound straddles an invalid sign: ' + name)
    return number.lower()


def _pack_lower(value):
    endpoint = value.lower()
    if not endpoint.is_finite() or not (endpoint > 0):
        raise ValueError('rank failure: lambda_min(A A^dagger) is not certified positive')
    packed = exact_upper(endpoint)
    packed['endpoint'] = 'lower'
    return packed


def _require_scope(caller_scope):
    if not isinstance(caller_scope, str) or not caller_scope.strip():
        raise ValueError('explicit partial-source/history scope required')
    return caller_scope.strip()


def _fiber_columns(columns):
    values = np.asarray(columns, complex)
    if values.shape != (2, 3) or not np.isfinite(values).all():
        raise ValueError('complete finite 2x3 same-energy fiber required')
    return values


def homogeneous_gram_lower(columns, *, bits=120):
    """Directed lower endpoint of lambda_min(A A^dagger) for one original 2x3 fiber.

    The eigenvalue is the smaller root of the exact 2x2 Hermitian Gram. An
    enclosure that is not strictly positive is a rank failure, including a
    ball that straddles zero. Column phases are not compared.
    """
    values = _fiber_columns(columns)
    with ctx.workprec(bits):
        matrix = [[complex_ball(entry) for entry in row] for row in values]
        gram_00 = sum((matrix[0][k] * matrix[0][k].conjugate() for k in range(3)), acb(0)).real
        gram_11 = sum((matrix[1][k] * matrix[1][k].conjugate() for k in range(3)), acb(0)).real
        gram_01 = sum((matrix[0][k] * matrix[1][k].conjugate() for k in range(3)), acb(0))
        trace = gram_00 + gram_11
        gap = gram_00 - gram_11
        # Only an upper discriminant is needed for a lower eigenvalue.
        # Magnitude uppers also preserve its known nonnegative sign when
        # cancellation makes an interval product straddle zero.
        discriminant = (gap.abs_upper()**2 + 4*gram_01.abs_upper()**2).upper()
        if not discriminant.is_finite() or not (discriminant >= 0):
            raise ValueError('Gram discriminant enclosure straddles an invalid sign')
        lesser = (trace - discriminant.sqrt()) / 2
        if not lesser.is_finite() or not (lesser > 0):
            raise ValueError('rank failure: lambda_min(A A^dagger) is not certified positive')
        return lesser.lower()


def covariance_insertion_factor(columns, epsilon, *, bits=120):
    """Upper bound eta>=epsilon/gamma. gamma is not replaced by 1."""
    if epsilon is None or isinstance(epsilon, bool):
        _missing('covariance operator error epsilon')
    with ctx.workprec(bits):
        gamma = homogeneous_gram_lower(columns, bits=bits)
        eps = _nonnegative_upper(epsilon, 'covariance operator error epsilon')
        if not gamma > 0:
            raise ValueError('rank failure: lambda_min(A A^dagger) is not certified positive')
        eta = (eps / gamma).upper()
        return {
            'gamma_lower': _pack_lower(gamma),
            'epsilon_upper': exact_upper(eps),
            'eta_upper': exact_upper(eta),
            'gram_identity_assumed': False,
            'ignored_gamma': False,
            'column_phases_compared': False,
        }


def _energy_sign(value, name='energy sign'):
    if isinstance(value, bool) or value is None:
        _missing(name)
    sign = int(value)
    if sign not in (-1, 1) or value != sign:
        raise ValueError('explicit energy sign +1 or -1 required; energies are not collapsed')
    return sign


def energy_diagonal_eta(fibers, *, bits=120, copy_positive_bound=False):
    """eta=max epsilon_e/gamma_e. Each fiber, including a negative energy, supplies epsilon."""
    if copy_positive_bound:
        raise ValueError(_NEGATIVE_COPY)
    if isinstance(fibers, (str, bytes)) or not isinstance(fibers, Sequence) or not fibers:
        raise ValueError('incomplete fibers')
    factors = []
    seen = []
    with ctx.workprec(bits):
        for index, fiber in enumerate(fibers):
            if not isinstance(fiber, Mapping):
                raise ValueError('incomplete fibers')
            epsilon = fiber.get('epsilon', None)
            sign = _energy_sign(fiber.get('energy_sign'))
            if epsilon is None or isinstance(epsilon, bool):
                if sign < 0:
                    raise ValueError(_NEGATIVE_COPY)
                _missing('covariance operator error epsilon')
            columns = _fiber_columns(fiber.get('columns'))
            if 'energy' in fiber and fiber['energy'] is not None:
                energy = float(fiber['energy'])
                if not np.isfinite(energy) or energy == 0 or int(np.sign(energy)) != sign:
                    raise ValueError('complete same-energy fibers required; energies are not collapsed')
                if any(np.sign(previous) == np.sign(energy) and previous == energy for previous in seen):
                    raise ValueError('duplicate energy fiber; do not split one energy across blocks')
                seen.append(energy)
            else:
                _missing('nonzero source energy')
            factor = covariance_insertion_factor(columns, epsilon, bits=bits)
            factor['energy_sign'] = sign
            factor['energy_hex'] = energy.hex()
            factor['fiber_index'] = index
            factors.append(factor)
        eta = restored_upper(factors[0]['eta_upper'])
        for factor in factors[1:]:
            eta = arb.max(eta, restored_upper(factor['eta_upper']))
        if not eta.is_finite() or not (eta >= 0):
            raise ValueError('interval bound straddles an invalid sign: eta')
        return {
            'eta_upper': exact_upper(eta.upper()),
            'fibers': tuple(factors),
            'negative_covariance_bound_copied_from_positive': False,
            'gram_identity_assumed': False,
            'energies_collapsed': False,
        }


def signed_covariance_factors(positive_columns, positive_epsilon, negative_columns,
                              negative_epsilon, *, bits=120, copy_positive_bound=False):
    """Separate factors. The negative epsilon is never defaulted from the positive one."""
    if copy_positive_bound or negative_epsilon is None or positive_epsilon is None:
        raise ValueError(_NEGATIVE_COPY)
    return {
        'positive': covariance_insertion_factor(positive_columns, positive_epsilon, bits=bits),
        'negative': covariance_insertion_factor(negative_columns, negative_epsilon, bits=bits),
        **_open_scope(),
    }


def require_energy_diagonal_preparation(columns, energies, covariance=None):
    """Reject an incomplete fiber, a collapsed sign, or a cross-energy coherence."""
    upstream = np.asarray(columns, complex)
    labels = np.asarray(energies, float)
    if (upstream.ndim != 2 or upstream.shape[0] != 2 or labels.ndim != 1
            or labels.shape != (upstream.shape[1],) or upstream.shape[1] % 3
            or not upstream.shape[1] or not np.isfinite(upstream).all()
            or not np.isfinite(labels).all()):
        raise ValueError('incomplete same-energy fibers: homogeneous upstream A must be 2x3 per energy')
    seen = []
    blocks = upstream.shape[1] // 3
    for index in range(blocks):
        block = labels[3 * index:3 * index + 3]
        if not np.all(block == block[0]) or block[0] == 0:
            raise ValueError('complete same-energy fibers required; energies are not collapsed')
        if any(previous == block[0] for previous in seen):
            raise ValueError('duplicate energy fiber; do not split one energy across blocks')
        seen.append(float(block[0]))
    if covariance is not None:
        source = np.asarray(covariance, complex)
        if source.shape != (upstream.shape[1], upstream.shape[1]) or not np.isfinite(source).all():
            raise ValueError('energy-diagonal covariance required')
        for left in range(blocks):
            for right in range(blocks):
                if left == right:
                    continue
                block = source[3 * left:3 * left + 3, 3 * right:3 * right + 3]
                if np.any(block != 0):
                    raise ValueError('energy-diagonal covariance required; cross-energy coherence is not an insertion fiber')
    return tuple(seen)


def numerical_right_inverse(columns):
    """R=A^dagger (A A^dagger)^-1. Validation only; not a certificate."""
    upstream = np.array(_fiber_columns(columns), dtype=complex, copy=True)
    gram = upstream @ upstream.conj().T
    eigenvalues = np.linalg.eigvalsh(gram)
    if eigenvalues.min() <= 0:
        raise ValueError('rank failure: lambda_min(A A^dagger) is not certified positive')
    return upstream.conj().T @ np.linalg.inv(gram)


def numerical_covariance_insertion(columns, delta_Q):
    """delta C=R delta Q R^dagger. Inputs are not modified. Not a certificate."""
    if delta_Q is None:
        _missing('delta Q')
    perturbation = np.array(delta_Q, dtype=complex, copy=True)
    if perturbation.shape != (2, 2) or not np.isfinite(perturbation).all():
        raise ValueError('finite 2x2 Hermitian delta Q required')
    if np.max(np.abs(perturbation - perturbation.conj().T)) > 1e-10:
        raise ValueError('finite 2x2 Hermitian delta Q required')
    right_inverse = numerical_right_inverse(columns)
    inserted = right_inverse @ perturbation @ right_inverse.conj().T
    return {
        'delta_C': inserted,
        'certificate': False,
        'original_columns_modified': False,
        'original_source_modified': False,
        'inhomogeneous_kernel_assumed_pointwise_unitary': False,
    }


def numerical_retarded_kernels(upstream, columns, axial, *, tolerance=1e-8):
    """Recover K,L from F=KA and F_z=LA. Not a certificate and not a unitary completion."""
    fiber = np.array(_fiber_columns(upstream), dtype=complex, copy=True)
    field = np.array(np.asarray(columns, complex), dtype=complex, copy=True)
    derivative = np.array(np.asarray(axial, complex), dtype=complex, copy=True)
    if (field.ndim != 3 or field.shape[1:] != (2, 3) or derivative.shape != field.shape
            or not np.isfinite(field).all() or not np.isfinite(derivative).all()):
        raise ValueError('weighted (nz,2,3) columns and actual axial derivatives required')
    right_inverse = numerical_right_inverse(fiber)
    kernel = np.einsum('zis,sa->zia', field, right_inverse, optimize=True)
    axial_kernel = np.einsum('zis,sa->zia', derivative, right_inverse, optimize=True)
    residual = np.einsum('zia,as->zis', kernel, fiber, optimize=True) - field
    residual_z = np.einsum('zia,as->zis', axial_kernel, fiber, optimize=True) - derivative
    if max(np.max(np.abs(residual)), np.max(np.abs(residual_z))) > tolerance:
        raise ValueError('columns are not one retarded image of this fiber; K is not a certificate')
    return {
        'K': kernel,
        'L': axial_kernel,
        'reconstruction_residual': float(np.max(np.abs(residual))),
        'axial_reconstruction_residual': float(np.max(np.abs(residual_z))),
        'certificate': False,
        'evolved_momentum_replaced_by_minus_energy': False,
        'inhomogeneous_kernel_assumed_pointwise_unitary': False,
    }


def numerical_action_operators(kernel, axial_kernel, *, mass, angular, axial_scale, radius, multiplicity):
    """H_N and H_beta for the supplied kernels. Validation only; not a certificate.

    V=-m sigma_1+(ell/r) sigma_2,
    H_N=-mu [K^dagger V K+(K^dagger sigma_3 L-L^dagger sigma_3 K)/(2 i a)],
    H_beta=mu (K^dagger L-L^dagger K)/(2 i).
    The same intrinsic a,r are the ones passed to source_column_matter.
    """
    if mass is None or angular is None or axial_scale is None or radius is None or multiplicity is None:
        _missing('action parameters')
    spin = np.asarray(kernel, complex)
    axial = np.asarray(axial_kernel, complex)
    if spin.ndim != 3 or spin.shape[1:] != (2, 2) or axial.shape != spin.shape:
        raise ValueError('matching (nz,2,2) retarded kernels required')
    scalars = [float(mass), float(angular), float(axial_scale), float(radius), float(multiplicity)]
    if not np.isfinite(scalars).all() or scalars[0] < 0 or min(scalars[2], scalars[3], scalars[4]) <= 0:
        raise ValueError('nonnegative mass and positive intrinsic a, r, multiplicity required')
    mass, angular, axial_scale, radius, mu = scalars
    vertices = raw_ks_vertex_coefficients(
        np.array([1.0, 0.0, axial_scale, radius]), mass, angular, envelopes=np.ones(4))
    potential = vertices.multiplication[0]
    axial_vertex = vertices.momentum[0]
    hamiltonian_N = np.empty_like(spin)
    hamiltonian_beta = np.empty_like(spin)
    for index, (current, current_z) in enumerate(zip(spin, axial)):
        adjoint = current.conj().T
        adjoint_z = current_z.conj().T
        radial = (adjoint @ axial_vertex @ current_z - adjoint_z @ axial_vertex @ current) / (2j)
        hamiltonian_N[index] = -mu * (adjoint @ potential @ current + radial)
        hamiltonian_beta[index] = mu * (adjoint @ current_z - adjoint_z @ current) / (2j)
    return {
        'H_N': hamiltonian_N,
        'H_beta': hamiltonian_beta,
        'intrinsic_axial_scale': axial_scale,
        'intrinsic_radius': radius,
        'certificate': False,
        'evolved_momentum_replaced_by_minus_energy': False,
    }


def numerical_operator_trace(covariance, operators):
    """Per-z Tr(Q H). Validation only; not a certificate."""
    if covariance is None:
        _missing('covariance')
    source = np.asarray(covariance, complex)
    if source.shape != (2, 2) or not np.isfinite(source).all():
        raise ValueError('finite 2x2 covariance required')
    rows = []
    for left, right in zip(operators['H_N'], operators['H_beta']):
        rows.append((
            complex(np.trace(source @ left)),
            complex(np.trace(source @ right)),
        ))
    traced = np.asarray(rows, complex)
    return {
        'action_trace': traced,
        'certificate': False,
        'history_minus_reference_cancellation': True,
    }


def compare_kernel_to_source_column_matter(upstream, columns, axial, covariance, *,
                                           mass, angular, axial_scale, radius, multiplicity,
                                           tolerance=1e-8):
    """Compare Tr(Q H) with source_column_matter. Not a certificate."""
    fiber = _fiber_columns(upstream)
    field = np.asarray(columns, complex)
    derivative = np.asarray(axial, complex)
    source = np.asarray(covariance, complex)
    if source.shape != (3, 3):
        raise ValueError('one 3x3 source block required')
    kernels = numerical_retarded_kernels(fiber, field, derivative, tolerance=tolerance)
    operators = numerical_action_operators(
        kernels['K'], kernels['L'], mass=mass, angular=angular,
        axial_scale=axial_scale, radius=radius, multiplicity=multiplicity)
    q_matrix = fiber @ source @ fiber.conj().T
    traced = numerical_operator_trace(q_matrix, operators)['action_trace'].real
    empty = np.zeros((0, *field.shape), complex)
    direct = source_column_matter(
        field, derivative, source, empty, empty, mass=mass, angular=angular,
        axial_scale=axial_scale, radius=radius, multiplicity=multiplicity)['action_gradient']
    discrepancy = float(np.max(np.abs(traced - direct)))
    if discrepancy > tolerance:
        raise ArithmeticError('kernel trace does not match source_column_matter')
    return {
        'discrepancy': discrepancy,
        'certificate': False,
        'evolved_momentum_replaced_by_minus_energy': False,
        'quadrature_applied_again': False,
        'intrinsic_axial_scale': float(axial_scale),
        'intrinsic_radius': float(radius),
    }


def source_difference_insertion_bound(reference_norm, difference_norm, reference_axial_norm,
                                      difference_axial_norm, reference_error, difference_error,
                                      reference_axial_error, difference_axial_error, eta, *,
                                      mass, absolute_angular, axial_lower, radius_lower,
                                      multiplicity, axial_scale=None, radius=None, bits=120):
    """Source-only N,beta difference bound after field errors enlarge the norms.

    density<=eta (2 f d+d^2) and momentum<=eta (f dz+fz d+d dz), with f,fz,d,dz
    replaced by captured upper + certified field-error upper. The vertex factor
    is the existing finite_matter_error contraction. Zero enlarged difference
    norms return a zero source-difference bound even when eta and the reference
    field are positive: the shared source error cancels in history minus reference.
    """
    arguments = (
        (reference_norm, 'reference norm'), (difference_norm, 'difference norm'),
        (reference_axial_norm, 'reference axial norm'), (difference_axial_norm, 'difference axial norm'),
        (reference_error, 'reference field error'), (difference_error, 'difference field error'),
        (reference_axial_error, 'reference axial field error'),
        (difference_axial_error, 'difference axial field error'),
        (eta, 'eta'), (mass, 'mass'), (absolute_angular, 'absolute angular'),
        (axial_lower, 'axial lower'), (radius_lower, 'radius lower'),
        (multiplicity, 'multiplicity'),
    )
    for value, name in arguments:
        if value is None or isinstance(value, bool):
            _missing(name)
    with ctx.workprec(bits):
        reference = _nonnegative_upper(reference_norm, 'reference norm')
        reference_axial = _nonnegative_upper(reference_axial_norm, 'reference axial norm')
        difference = _nonnegative_upper(difference_norm, 'difference norm')
        difference_axial = _nonnegative_upper(difference_axial_norm, 'difference axial norm')
        reference_error_upper = _nonnegative_upper(reference_error, 'reference field error')
        difference_error_upper = _nonnegative_upper(difference_error, 'difference field error')
        reference_axial_error_upper = _nonnegative_upper(
            reference_axial_error, 'reference axial field error')
        difference_axial_error_upper = _nonnegative_upper(
            difference_axial_error, 'difference axial field error')
        eta_upper = _nonnegative_upper(eta, 'eta')
        mass_upper = _nonnegative_upper(mass, 'mass')
        angular_upper = _nonnegative_upper(absolute_angular, 'absolute angular')
        multiplicity_upper = _positive_magnitude(multiplicity, 'multiplicity')
        enlarged_reference = (reference + reference_error_upper).upper()
        enlarged_reference_axial = (reference_axial + reference_axial_error_upper).upper()
        enlarged_difference = (difference + difference_error_upper).upper()
        enlarged_difference_axial = (difference_axial + difference_axial_error_upper).upper()
        axial_bound = _positive_lower(axial_lower, 'axial lower')
        radius_bound = _positive_lower(radius_lower, 'radius lower')
        if axial_scale is not None:
            declared_axial = _positive_magnitude(axial_scale, 'intrinsic axial scale')
            if not (declared_axial >= axial_bound):
                raise ValueError('intrinsic axial scale does not match the declared slice')
        if radius is not None:
            declared_radius = _positive_magnitude(radius, 'intrinsic radius')
            if not (declared_radius >= radius_bound):
                raise ValueError('intrinsic radius does not match the declared slice')
        bound = finite_matter_error(
            enlarged_reference, enlarged_reference_axial, enlarged_difference,
            enlarged_difference_axial, eta_upper, mass=mass_upper,
            absolute_angular=angular_upper, axial_lower=axial_bound,
            radius_lower=radius_bound, multiplicity=multiplicity_upper, bits=bits)
    result = dict(bound)
    result.update(_open_scope())
    result.update({
        'field_errors_added_to_norm_uppers': True,
        'enlarged_reference_norm': exact_upper(enlarged_reference),
        'enlarged_reference_axial_norm': exact_upper(enlarged_reference_axial),
        'enlarged_difference_norm': exact_upper(enlarged_difference),
        'enlarged_difference_axial_norm': exact_upper(enlarged_difference_axial),
        'eta_upper': exact_upper(eta_upper),
        'source_only_difference_bound': True,
        'quadrature_applied_again': False,
    })
    return result


def insert_source_difference(fibers, *, reference_norm, difference_norm, reference_axial_norm,
                             difference_axial_norm, reference_error, difference_error,
                             reference_axial_error, difference_axial_error, mass,
                             absolute_angular, axial_lower, radius_lower, multiplicity,
                             caller_scope, axial_scale=None, radius=None, bits=120,
                             copy_positive_bound=False):
    """Combine per-fiber eta with the enlarged source-difference bound."""
    scope = _require_scope(caller_scope)
    factor = energy_diagonal_eta(fibers, bits=bits, copy_positive_bound=copy_positive_bound)
    # Pass the packed upper through so it is restored at the bound's working precision.
    bound = source_difference_insertion_bound(
        reference_norm, difference_norm, reference_axial_norm, difference_axial_norm,
        reference_error, difference_error, reference_axial_error, difference_axial_error,
        factor['eta_upper'], mass=mass, absolute_angular=absolute_angular,
        axial_lower=axial_lower, radius_lower=radius_lower, multiplicity=multiplicity,
        axial_scale=axial_scale, radius=radius, bits=bits)
    bound['caller_scope'] = scope
    bound['fiber_factors'] = factor['fibers']
    bound['eta_upper'] = factor['eta_upper']
    bound['negative_covariance_bound_copied_from_positive'] = False
    bound['physical_source_coverage'] = False
    bound['full_upstream_gate'] = 'OPEN'
    return bound


def _slot_value(matter, name):
    if not isinstance(matter, Mapping) or name not in matter:
        _missing(name)
    slot = matter[name]
    if slot is None:
        _missing(name)
    if isinstance(slot, Mapping):
        if 'value' in slot:
            if slot['value'] is None:
                _missing(name)
            return slot['value']
        if 'mantissa' not in slot or 'exponent' not in slot:
            _missing(name)
    return slot


def insert_from_endpoint_slots(matter, fibers, *, caller_scope, bits=120,
                               copy_positive_bound=False):
    """Read named norm and error slots. Does not assemble a physical budget.

    ``source_norm`` is not an epsilon. The caller states the partial history
    and source scope; this function still reports the upstream gate OPEN.
    """
    if not isinstance(matter, Mapping):
        raise TypeError('endpoint matter slots required')
    values = {name: _slot_value(matter, name) for name in ENDPOINT_FIELD_SLOTS}
    result = insert_source_difference(
        fibers, caller_scope=caller_scope, bits=bits, copy_positive_bound=copy_positive_bound,
        **values)
    result['endpoint_slot_source_norm_used_as_epsilon'] = False
    result['caller_scope'] = _require_scope(caller_scope)
    result['partial_history_scope'] = True
    result['partial_source_scope'] = True
    return result


def insert_from_upstream_batch(batch, epsilons, norms, *, caller_scope, bits=120,
                               copy_positive_bound=False):
    """Use one KSUpstreamBatch fiber layout. Epsilons stay caller-supplied."""
    if copy_positive_bound or (batch is not None and getattr(batch, 'energy_sign', 1) < 0
                               and (epsilons is None or copy_positive_bound)):
        raise ValueError(_NEGATIVE_COPY)
    if not isinstance(batch, KSUpstreamBatch):
        raise TypeError('KSUpstreamBatch required')
    if epsilons is None:
        if batch.energy_sign < 0:
            raise ValueError(_NEGATIVE_COPY)
        _missing('covariance operator error epsilon')
    if isinstance(epsilons, np.ndarray):
        epsilons = tuple(np.asarray(epsilons).reshape(-1).tolist())
    elif isinstance(epsilons, (str, bytes, int, float, np.floating, np.integer, arb)):
        epsilons = (epsilons,)
    if isinstance(epsilons, (str, bytes)) or not isinstance(epsilons, Sequence) or any(
            item is None or isinstance(item, bool) for item in epsilons):
        if batch.energy_sign < 0:
            raise ValueError(_NEGATIVE_COPY)
        _missing('covariance operator error epsilon')
    upstream = np.asarray(batch.initial_columns, complex)
    labels = np.asarray(batch.source.energies, float)
    require_energy_diagonal_preparation(upstream, labels, batch.source.covariance)
    count = upstream.shape[1] // 3
    if len(epsilons) != count:
        _missing('covariance operator error epsilon')
    fibers = [{
        'columns': upstream[:, 3 * index:3 * index + 3],
        'epsilon': epsilons[index],
        'energy': float(labels[3 * index]),
        'energy_sign': int(batch.energy_sign),
    } for index in range(count)]
    if not isinstance(norms, Mapping):
        raise TypeError('explicit field norm and error slots required')
    with ctx.workprec(bits):
        if (_nonnegative_upper(norms.get('mass'), 'mass') < arb(batch.mass)
                or _nonnegative_upper(norms.get('absolute_angular'), 'absolute angular') < abs(arb(batch.angular))):
            raise ValueError('channel bounds do not cover the bound upstream batch')
    result = insert_source_difference(
        fibers, caller_scope=caller_scope, bits=bits, copy_positive_bound=False,
        reference_norm=norms.get('reference_norm', None),
        difference_norm=norms.get('difference_norm', None),
        reference_axial_norm=norms.get('reference_axial_norm', None),
        difference_axial_norm=norms.get('difference_axial_norm', None),
        reference_error=norms.get('reference_error', None),
        difference_error=norms.get('difference_error', None),
        reference_axial_error=norms.get('reference_axial_error', None),
        difference_axial_error=norms.get('difference_axial_error', None),
        mass=norms.get('mass', None),
        absolute_angular=norms.get('absolute_angular', None),
        axial_lower=norms.get('axial_lower', None),
        radius_lower=norms.get('radius_lower', None),
        multiplicity=norms.get('multiplicity', None),
        axial_scale=norms.get('axial_scale', None),
        radius=norms.get('radius', None))
    result['upstream_panel'] = batch.original_panel
    result['upstream_rows'] = tuple(batch.rows)
    result['energy_sign'] = int(batch.energy_sign)
    result['angular_sign'] = int(batch.angular_sign)
    result['preparation_digest'] = batch.preparation_digest
    result['source_digest'] = batch.source.digest
    result['source_evidence_authenticated'] = False
    result['physical_source_coverage'] = False
    result['preparation_ode_rerun'] = False
    result['archived_columns_modified'] = False
    return result
