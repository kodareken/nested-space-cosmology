"""Numerical incoming-boundary buffer by extending the computational right edge.

The physical source, signed energies, three fibers, C_src and clock remain
the archived group14_1 data. Only rho>1.2 nodes are appended, using the owned
canonical KS generator and the actual restriction/inverse. Finite-grid
causality is not claimed; this is not a new source or metric history.
"""
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.sparse import bmat, csr_matrix, diags
from scipy.sparse.linalg import expm_multiply

from .nsc_common_ks_trace import incoming_KS_state, restrict_resolved_modes
from .nsc_evolved_incoming_constraints import source_column_matter
from .nsc_evolved_incoming_state import rho1_node_index, rho1_restriction_map
from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_lorentzian import geometry, horizon
from .nsc_retarded_radial_response import archived_amplitudes, canonical_map, ks_generator
from .nsc_transmitting_dirac_domain import MODE_TO_CURRENT, TransmittingDiracSeamDomain
from .nsc_transmitting_history_jets import incoming_reference_columns
from .nsc_transmitting_history_modes import characteristic_columns


ORIGINAL_LEFT = -2.
ORIGINAL_RIGHT = 1.2
EXTENDED_RIGHT = 1.8
ORIGINAL_COARSE_POINTS = 801
ORIGINAL_FINE_POINTS = 3201
COARSE_POINTS = 951
FINE_POINTS = 3801
COARSE_SPACING = .004
FINE_SPACING = .001
TIME_STOP = .3
SAMPLE_TIMES = 65
RTOL = 2e-13
ATOL = 2e-15
MAX_STEP = 1 / 2000
MULTIPLICITY = 12.
DRIFT_TARGET = 3e-11
BOUNDARY_STENCIL = 4
CPU_CAP = 60.


def maxabs(value):
    return float(np.max(np.abs(value), initial=0.))


def uniform_grid(left, right, points):
    x = np.linspace(float(left), float(right), int(points))
    if x.ndim != 1 or len(x) < 9 or np.any(np.diff(x) <= 0):
        raise ValueError('ordered uniform buffer grid required')
    return x


def original_grid(points):
    if points not in (ORIGINAL_COARSE_POINTS, ORIGINAL_FINE_POINTS):
        raise ValueError('original prefix has 801 or 3201 nodes')
    return uniform_grid(ORIGINAL_LEFT, ORIGINAL_RIGHT, points)


def extended_grid(points):
    if points not in (COARSE_POINTS, FINE_POINTS):
        raise ValueError('extended buffer has 951 or 3801 nodes')
    x = uniform_grid(ORIGINAL_LEFT, EXTENDED_RIGHT, points)
    # linspace's explicit endpoint differs by one ulp from the corresponding
    # interior node of the longer grid. Preserve authenticated prefix bytes.
    old_points = ORIGINAL_FINE_POINTS if points == FINE_POINTS else ORIGINAL_COARSE_POINTS
    x[:old_points] = original_grid(old_points)
    if float(x[-1]) >= float(horizon()):
        raise ValueError('extended edge must remain strictly below the fixed horizon')
    return x


def nested_prefix(x, phi, points):
    x = np.asarray(x)
    phi = np.asarray(phi)
    if points < 2 or (len(x) - 1) % (points - 1) or phi.shape != (2 * len(x), 12):
        raise ValueError('nested source grid and twelve spin-major columns required')
    stride = (len(x) - 1) // (points - 1)
    return x[::stride], phi.reshape(2, len(x), 12)[:, ::stride].reshape(2 * points, 12)


def prefix_fields(x_old, phi_old, x_new, phi_new):
    n = len(x_old)
    if len(x_new) < n or phi_old.shape != (2 * n, 12) or phi_new.shape != (2 * len(x_new), 12):
        raise ValueError('extended fields must contain the original prefix layout')
    return x_new[:n], phi_new.reshape(2, len(x_new), 12)[:, :n].reshape(2 * n, 12)


def prefixes_equal(x_old, phi_old, x_new, phi_new):
    x_pref, phi_pref = prefix_fields(x_old, phi_old, x_new, phi_new)
    return np.array_equal(x_pref, x_old) and np.array_equal(phi_pref, phi_old)


def pg_from_characteristic(phi, index):
    phi = np.asarray(phi, complex)
    n = phi.shape[0] // 2
    if phi.shape != (2 * n, 12) or not 0 <= index < n:
        raise ValueError('spin-major twelve-column characteristic node required')
    current = MODE_TO_CURRENT @ phi[[index, n + index]]
    return current.reshape(2, 4, 3).transpose(1, 0, 2)


def ks_from_characteristic(x, phi, energies, index):
    """Owned canonical restriction plus the existing clock phase."""
    pg = pg_from_characteristic(phi, index)
    canonical, chart = restrict_resolved_modes(float(x[index]), np.asarray(energies, float), pg)
    archived = archived_amplitudes(x, phi, energies, index)
    if maxabs(canonical - archived) > 3e-11:
        raise ArithmeticError('canonical restriction disagrees with the archived KS map')
    return canonical, chart


def pg_from_ks_amplitudes(rho, energies, amplitudes):
    """Inverse of restrict_resolved_modes: trace map times radius, with the clock removed."""
    E = np.asarray(energies, float)
    A = np.asarray(amplitudes, complex)
    if E.shape != (4,) or A.shape != (4, 2, 3) or not np.isfinite(A).all():
        raise ValueError('four finite three-fiber KS amplitudes required')
    rho = float(rho)
    beta = float(geometry(rho)[0])
    radius = float(np.sqrt(1. + rho * rho))
    domain = TransmittingDiracSeamDomain(lapse=1., radial_scale=1., shift=beta, radius=radius,
                                         surface_id=f'buffer-KS-inverse-rho={rho:g}')
    trace, inverse, _ = domain.trace_map(np.ones(1))
    if maxabs(inverse[0] / radius - canonical_map(rho)) > 3e-11:
        raise ArithmeticError('inverse frame is not the owned canonical map')
    clock = chart_coordinates(rho)[1]
    unclocked = A * np.exp(-1j * E[:, None, None] * clock)
    psi = np.einsum('ab,ebs->eas', trace[0], unclocked)
    return radius * psi


def characteristic_from_ks_path(radii, energies, amplitudes):
    radii = np.asarray(radii, float)
    amplitudes = np.asarray(amplitudes, complex)
    if amplitudes.shape != (len(radii), 4, 2, 3):
        raise ValueError('one KS amplitude block per new radius required')
    fields = np.stack([pg_from_ks_amplitudes(rho, energies, block)
                       for rho, block in zip(radii, amplitudes)], axis=0)
    return characteristic_columns(fields)


def continue_canonical_amplitudes(energies, mass, angular, initial, sample_rho, *,
                                  rho_start=ORIGINAL_RIGHT, rho_stop=EXTENDED_RIGHT,
                                  rtol=RTOL, atol=ATOL, max_step=MAX_STEP):
    """One coupled DOP853 solve of G_E on the new .6 interval; both grids sample it."""
    E = np.asarray(energies, float)
    A0 = np.asarray(initial, complex)
    sample_rho = np.asarray(sample_rho, float)
    if E.shape != (4,) or A0.shape != (4, 2, 3) or not np.isfinite(A0).all():
        raise ValueError('four energies and three source fibers required')
    if sample_rho.ndim != 1 or len(sample_rho) < 1 or np.any(np.diff(sample_rho) <= 0):
        raise ValueError('ordered new radii required')
    if np.any(sample_rho <= rho_start) or float(sample_rho[-1]) > float(rho_stop) + 1e-12:
        raise ValueError('append only new nodes with original_right < rho <= extended_right')
    if abs(float(rho_stop) - float(rho_start) - .6) > 1e-12:
        raise ValueError('the authorized continuation occupies only the new .6 interval')

    def rhs(rho, state):
        G = ks_generator(float(rho), E, mass, angular)
        return np.einsum('eab,ebc->eac', G, state.reshape(4, 2, 3)).ravel()

    run = solve_ivp(rhs, (float(rho_start), float(rho_stop)), np.asarray(A0.ravel(), complex),
                    method='DOP853', rtol=float(rtol), atol=float(atol),
                    max_step=float(max_step), dense_output=True)
    if not run.success:
        raise ArithmeticError(run.message)
    sampled = np.asarray(run.sol(sample_rho).T, complex).reshape((-1, 4, 2, 3))
    return sampled, {
        'function_evaluations': int(run.nfev),
        'accepted_steps': int(len(run.t) - 1),
        'rtol': float(rtol),
        'atol': float(atol),
        'max_step': float(max_step),
        'rho_start': float(rho_start),
        'rho_stop': float(rho_stop),
    }


def append_extension(x_old, phi_old, x_new, extension):
    n_old = len(x_old)
    n_new = len(x_new)
    extension = np.asarray(extension, complex)
    if not np.array_equal(x_new[:n_old], x_old):
        raise ValueError('extended grid must retain the original coordinates byte-for-byte')
    n_add = n_new - n_old
    if extension.shape != (2 * n_add, 12) or phi_old.shape != (2 * n_old, 12):
        raise ValueError('extension must supply only the new spin-major nodes')
    phi = np.empty((2 * n_new, 12), complex)
    phi[:n_old] = phi_old[:n_old]
    phi[n_old:n_new] = extension[:n_add]
    phi[n_new:n_new + n_old] = phi_old[n_old:]
    phi[n_new + n_old:] = extension[n_add:]
    return phi


def fiber_gram_residuals(amplitudes, projector):
    A = np.asarray(amplitudes, complex)
    P = np.asarray(projector, complex)
    if A.ndim == 3:
        A = A[None]
    if A.shape[-3:] != (4, 2, 3) or P.shape != (4, 3, 3):
        raise ValueError('KS amplitudes and the original 4x3x3 projector required')
    adjoint = A.swapaxes(-1, -2).conj()
    gram = A @ adjoint
    closed = A @ P - A
    identity = np.eye(2, dtype=complex)
    return {
        'canonical_gram': maxabs(gram - identity),
        'source_fiber_closure': maxabs(closed),
        'closed_third_column': maxabs(A[..., 2]),
    }


def slice_state_residuals(energies, pg_fields, source, projector, rho):
    state = incoming_KS_state(energies, pg_fields, source, projector, rho=float(rho))
    return {key: float(state['residuals'][key])
            for key in ('source_coisometry', 'closed_source_projection',
                        'source_complement_projector', 'Hermiticity')}


def continuum_incoming_buffer(left=1., right=EXTENDED_RIGHT, duration=TIME_STOP):
    """Lower bound .8/(1+beta(1)) from decreasing beta on [1,1.8]; not a grid theorem."""
    left, right, duration = float(left), float(right), float(duration)
    beta_left = float(geometry(left)[0])
    beta_right = float(geometry(right)[0])
    if not beta_right < beta_left:
        raise ArithmeticError('background beta must decrease on the buffer interval')
    lower = (right - left) / (1. + beta_left)
    travel = quad(lambda rho: 1. / (1. + float(geometry(rho)[0])), left, right,
                  epsabs=2e-13, epsrel=2e-13)[0]
    return {
        'interval': [left, right],
        'horizon': float(horizon()),
        'beta_left': beta_left,
        'beta_right': beta_right,
        'fast_speed_left': 1. + beta_left,
        'travel_time_lower_bound': lower,
        'travel_time_quadrature': float(travel),
        'inherited_duration': duration,
        'causal_margin_lower_bound': lower - duration,
        'finite_grid_exact_causality': False,
        'characteristic_speeds_changed': False,
    }


def sigma_parameters(owner, multiplicity=MULTIPLICITY):
    index = rho1_node_index(owner.x)
    lapse, shift, radial_scale, radius = (float(v) for v in owner.reference[:, index])
    domain = TransmittingDiracSeamDomain(lapse=lapse, radial_scale=radial_scale,
                                         shift=shift, radius=radius,
                                         surface_id='buffer-fixed-rho1')
    return {
        'mass': float(owner.mass),
        'angular': float(owner.angular),
        'axial_scale': float(domain.induced_axial_scale),
        'radius': radius,
        'multiplicity': float(multiplicity),
    }


def initial_stationary_defect(owner, phi, energies, stencil=BOUNDARY_STENCIL):
    forcing = incoming_reference_columns(owner, phi)
    generator, algebra = owner.generator(owner.reference)
    energy = np.repeat(np.asarray(energies, float), 3)
    residual = generator @ phi + forcing + 1j * energy * phi
    n = len(owner.x)
    spins = []
    for spin in range(2):
        block = residual[spin * n:(spin + 1) * n]
        interior = block[stencil:-stencil]
        spins.append({
            'right_boundary': maxabs(block[-1]),
            'left_boundary': maxabs(block[0]),
            'interior_excluding_four_boundary_rows': maxabs(interior),
            'all': maxabs(block),
        })
    return residual, spins, float(algebra['norm_flux_algebra'])


def exact_harmonic_reference(owner, phi, energies, times):
    """Exact source-phase block of evolve_field_jets; F_z is the actual PDE row."""
    times = np.asarray(times, float)
    if (times.ndim != 1 or len(times) < 2 or not np.isfinite(times).all()
            or times[0] != 0 or np.any(np.diff(times) <= 0)
            or not np.allclose(np.diff(times), np.diff(times)[0], rtol=1e-12, atol=1e-15)):
        raise ValueError('uniform increasing sample times starting at the archived zero required')
    energy = np.repeat(np.asarray(energies, float), 3)
    forcing = incoming_reference_columns(owner, phi)
    generator, algebra = owner.generator(owner.reference)
    rows, sources = phi.shape
    augmented = bmat([[generator, csr_matrix(forcing)],
                      [None, diags(-1j * energy)]], format='csr')
    initial = np.vstack((phi, np.eye(sources, dtype=complex)))
    solution = expm_multiply(augmented, initial, start=float(times[0]), stop=float(times[-1]),
                             num=len(times), endpoint=True, traceA=augmented.diagonal().sum())
    index = rho1_node_index(owner.x)
    take = np.array([index, len(owner.x) + index])
    lapse, shift, radial_scale, radius = owner.reference[:, index]
    restriction = rho1_restriction_map(lapse, radial_scale, shift, radius)
    fields = solution[:, take, :]
    columns = np.einsum('ab,tbs->tas', restriction, fields)
    axial = np.array([restriction @ (generator[take] @ state[:rows] + forcing[take] @ state[rows:])
                      for state in solution])
    expected = np.eye(sources)[None] * np.exp(-1j * times[:, None, None] * energy)
    reference = (restriction @ phi[take])[None] * np.exp(-1j * times[:, None, None] * energy)
    return {
        'columns': columns,
        'pde_axial_columns': axial,
        'reference_columns': reference,
        'reference_axial_columns': -1j * energy * reference,
        'phase': solution[:, rows:],
        'expected_phase': expected,
        'flux_algebra': float(algebra['norm_flux_algebra']),
        'restriction': restriction,
        'sample_rows': take,
    }


def weighted_matter(columns, axial, covariance, weights, parameters):
    zero = np.zeros((0, *np.asarray(columns).shape), complex)
    return source_column_matter(np.asarray(columns) * weights, np.asarray(axial) * weights,
                                covariance, zero, zero, **parameters)['action_gradient']


def matter_drift(columns, axial, reference_columns, reference_axial,
                 covariance, weights, parameters, times):
    value = weighted_matter(columns, axial, covariance, weights, parameters)
    reference = weighted_matter(reference_columns, reference_axial, covariance, weights, parameters)
    difference = value - reference
    per_time = np.max(np.abs(difference), axis=1)
    index = int(np.argmax(per_time))
    return {
        'PDE_matter_drift': np.max(np.abs(difference), axis=0).tolist(),
        'PDE_values': value.tolist(),
        'reference_values': reference.tolist(),
        'peak_time': float(times[index]),
        'peak_matter_drift': float(per_time[index]),
        'meets_partial_control_target': bool(np.max(np.abs(difference)) < DRIFT_TARGET),
    }
