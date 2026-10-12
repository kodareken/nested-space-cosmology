"""Reference plus difference coordinates for the same KS Dirac evolution.

The homogeneous reference A and D=X-A evolve in one ODE.  Only D and the
retarded tangents are Fourier differentiated.  No measured stress is removed.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import DOP853
from scipy.sparse.linalg import LinearOperator, expm_multiply

from . import nsc_ks_source_envelope as E


@dataclass(frozen=True)
class KSDifferenceIncoming(E.KSEnvelopeIncoming):
    reference_amplitudes: object
    envelope_difference: object
    envelope_difference_z: object

    def __post_init__(self):
        super().__post_init__()
        for name, shape in (
            ('reference_amplitudes', self.initial_columns.shape),
            ('envelope_difference', self.columns.shape),
            ('envelope_difference_z', self.columns.shape),
        ):
            value = E._finite_complex(getattr(self, name), name)
            if value.shape != shape:
                raise ValueError(f'{name} has the wrong shape')
            object.__setattr__(self, name, value)
        phase = np.exp(-1j * self.source_energies * self.z[:, None])[:, None, :]
        total = self.reference_amplitudes[None] + self.envelope_difference
        if not np.array_equal(self.columns, phase * total):
            raise ValueError('columns must reconstruct the evolved reference plus difference')
        if not np.array_equal(self.axial_columns, phase * (
                self.envelope_difference_z - 1j * self.source_energies * total)):
            raise ValueError('F_z must include the evolved difference derivative')


def _radius_multipliers(axial, angular, r_ref, delta_r, basis):
    """M = correction S2 and K_j = k_coeff_j S2 from the owned radius law."""
    delta_r = np.asarray(delta_r)
    radius = r_ref + delta_r
    # This algebraic form avoids cancellation of two nearly equal reciprocals.
    correction = -(1j * angular / axial) * delta_r / (r_ref * radius)
    basis = np.asarray(basis)
    if basis.size == 0:
        k_coeff = np.zeros((0, delta_r.shape[-1]), complex)
    else:
        if basis.ndim == 1:
            basis = basis[None]
        k_coeff = -(1j * angular / axial) * basis / radius**2
    return correction, k_coeff


def _apply_L(field, G, axial, correction, wave, adjoint):
    """L = Fourier spatial + G + M on a z-field. Spatial adjoint flips ik."""
    sign = -1j if adjoint else 1j
    dz = np.fft.ifft(sign * wave * np.fft.fft(field, axis=-1), axis=-1)
    spatial = np.einsum('ab,...bsz->...asz', E.S3, dz) / axial**2
    if adjoint:
        generated = np.einsum('sab,...asz->...bsz', np.conjugate(G), field)
        potential = np.conjugate(correction) * np.einsum('ab,...bsz->...asz', E.S2, field)
    else:
        generated = np.einsum('sab,...bsz->...asz', G, field)
        potential = correction * np.einsum('ab,...bsz->...asz', E.S2, field)
    return spatial + generated + potential


def _g_ref(G, A, adjoint):
    if adjoint:
        return np.einsum('sab,...as->...bs', np.conjugate(G), A)
    return np.einsum('sab,...bs->...as', G, A)


def _difference_primal_fields(A, D, Y, G, axial, correction, k_coeff, wave):
    A_dot = _g_ref(G, A, False)
    D_dot = _apply_L(D, G, axial, correction, wave, False)
    D_dot = D_dot + correction * np.einsum('ab,...bs->...as', E.S2, A)[..., None]
    Y_dot = _apply_L(Y, G, axial, correction, wave, False)
    if Y.shape[1]:
        spin = np.einsum('ab,...bsz->...asz', E.S2, A[..., None] + D)
        Y_dot = Y_dot + k_coeff[:, None, None, :] * spin[:, None]
    return A_dot, D_dot, Y_dot


def _difference_adjoint_fields(A, D, Y, G, axial, correction, k_coeff, wave):
    A_dot = _g_ref(G, A, True)
    D_dot = _apply_L(D, G, axial, correction, wave, True)
    Y_dot = _apply_L(Y, G, axial, correction, wave, True)
    A_dot = A_dot + np.sum(np.conjugate(correction) * np.einsum('ab,...bsz->...asz', E.S2, D), axis=-1)
    if Y.shape[1]:
        weighted = np.conjugate(k_coeff)[:, None, None, :] * np.einsum('ab,...bsz->...asz', E.S2, Y)
        A_dot = A_dot + weighted.sum(axis=(-4, -1))
        D_dot = D_dot + weighted.sum(axis=-4)
    return A_dot, D_dot, Y_dot


def difference_rhs(A, D, Y, G, axial, angular, r_ref, delta_r, basis, wave):
    """Exact change of coordinates: A'=G A, (A+D)'=L_g(A+D)."""
    A, D, Y = np.asarray(A, complex), np.asarray(D, complex), np.asarray(Y, complex)
    correction, k_coeff = _radius_multipliers(axial, angular, r_ref, delta_r, basis)
    Ad, Dd, Yd = _difference_primal_fields(
        A[None], D[None], Y[None], G, axial, correction, k_coeff, wave)
    return Ad[0], Dd[0], Yd[0]


def difference_rhs_adjoint(A, D, Y, G, axial, angular, r_ref, delta_r, basis, wave):
    """Conjugate transpose of difference_rhs on the discrete (A, D, Y) inner product."""
    A, D, Y = np.asarray(A, complex), np.asarray(D, complex), np.asarray(Y, complex)
    correction, k_coeff = _radius_multipliers(axial, angular, r_ref, delta_r, basis)
    Ad, Dd, Yd = _difference_adjoint_fields(
        A[None], D[None], Y[None], G, axial, correction, k_coeff, wave)
    return Ad[0], Dd[0], Yd[0]


def difference_operator_trace(G, nsrc, nz, ndir, axial, correction, wave):
    """Analytic trace of the augmented vector field. G is not assumed traceless."""
    G = np.asarray(G, complex)
    g_trace = np.trace(G, axis1=-2, axis2=-1)
    if np.ndim(g_trace):
        g_trace = g_trace.sum()
    ndir = int(ndir)
    nz = int(nz)
    nsrc = int(nsrc)
    multiplicity = 1 + nz * (1 + ndir)
    s2_trace = np.trace(E.S2)
    s3_trace = np.trace(E.S3)
    m_trace = s2_trace * np.sum(np.asarray(correction)) * nsrc * (1 + ndir)
    spatial_trace = (s3_trace / axial**2) * (1j * np.sum(np.asarray(wave))) * nsrc * (1 + ndir)
    return g_trace * multiplicity + m_trace + spatial_trace


class KSDifferenceOperator(LinearOperator):
    """Vector field of the augmented (A, D, Y) system at one rho, with adjoint."""

    def __init__(self, G, axial, angular, r_ref, delta_r, basis, wave, evaluations=None):
        G = np.asarray(G, complex)
        if G.ndim != 3 or G.shape[-2:] != (2, 2):
            raise ValueError('G must have shape (nsrc, 2, 2)')
        self.G = G
        self.axial = float(axial)
        self.angular = float(angular)
        self.r_ref = float(r_ref)
        self.delta_r = np.asarray(delta_r, float)
        self.wave = np.asarray(wave, float)
        basis = np.asarray(basis, float)
        self.basis = np.zeros((0, self.delta_r.shape[-1]), float) if basis.size == 0 else np.atleast_2d(basis)
        if self.wave.shape != self.delta_r.shape:
            raise ValueError('wave must match the z-grid')
        self.nsrc = int(G.shape[0])
        self.nz = int(self.delta_r.shape[-1])
        self.ndir = int(self.basis.shape[0])
        self.size_A = 2 * self.nsrc
        self.size_D = self.size_A * self.nz
        self.correction, self.k_coeff = _radius_multipliers(
            self.axial, self.angular, self.r_ref, self.delta_r, self.basis)
        self.trace = difference_operator_trace(
            G, self.nsrc, self.nz, self.ndir, self.axial, self.correction, self.wave)
        self._evaluations = evaluations
        super().__init__(complex, (self.size_A + self.size_D * (1 + self.ndir),) * 2)

    def _unpack(self, vector):
        vector = np.asarray(vector, complex)
        squeeze = vector.ndim == 1
        if squeeze:
            vector = vector[:, None]
        elif vector.ndim != 2:
            raise ValueError('expected a vector or column matrix')
        if vector.shape[0] != self.shape[0]:
            raise ValueError('state does not match the augmented (A, D, Y) length')
        k = vector.shape[1]
        A = np.moveaxis(vector[:self.size_A], -1, 0).reshape(k, 2, self.nsrc)
        D = np.moveaxis(vector[self.size_A:self.size_A + self.size_D], -1, 0).reshape(
            k, 2, self.nsrc, self.nz)
        Y = np.moveaxis(vector[self.size_A + self.size_D:], -1, 0).reshape(
            k, self.ndir, 2, self.nsrc, self.nz)
        return A, D, Y, squeeze

    def _pack(self, A, D, Y, squeeze):
        k = A.shape[0]
        parts = (
            np.moveaxis(np.reshape(A, (k, -1)), 0, -1),
            np.moveaxis(np.reshape(D, (k, -1)), 0, -1),
            np.moveaxis(np.reshape(Y, (k, -1)), 0, -1),
        )
        out = np.concatenate(parts, axis=0)
        return out[:, 0] if squeeze else out

    def _apply(self, vector, adjoint):
        A, D, Y, squeeze = self._unpack(vector)
        if self._evaluations is not None:
            self._evaluations[0] += int(A.shape[0])
        apply = _difference_adjoint_fields if adjoint else _difference_primal_fields
        Ad, Dd, Yd = apply(A, D, Y, self.G, self.axial, self.correction, self.k_coeff, self.wave)
        return self._pack(Ad, Dd, Yd, squeeze)

    def _matvec(self, x):
        return self._apply(x, False)

    def _rmatvec(self, x):
        return self._apply(x, True)

    def _matmat(self, X):
        return self._apply(X, False)

    def _rmatmat(self, X):
        return self._apply(X, True)


# Blanes–Moan (2006) equation 43: two Gauss nodes, two exponentials.
_CF4_A1 = (3 - 2 * np.sqrt(3.)) / 12
_CF4_A2 = (3 + 2 * np.sqrt(3.)) / 12
_CF4_C1 = 0.5 - np.sqrt(3.) / 6
_CF4_C2 = 0.5 + np.sqrt(3.) / 6
# Gauss nodes A1=A(left+c1 h), A2=A(left+c2 h). The 4th-order product, valid for
# signed h, is exp(h(α1 A1+α2 A2)) exp(h(α2 A1+α1 A2)): apply (α2,α1) then (α1,α2).
_CF4_WEIGHTS = ((_CF4_A2, _CF4_A1), (_CF4_A1, _CF4_A2))


def _cf4_gauss_nodes(left, right):
    span = float(right - left)
    return span, (float(left) + _CF4_C1 * span, float(left) + _CF4_C2 * span)


def _cf4_step(op0, op1, state, span, reverse=False):
    """One CF4 step of y' = G(rho) y. span = right-left may be negative."""
    span = float(span)
    if span == 0.0:
        return np.asarray(state, complex)
    pairs = _CF4_WEIGHTS[::-1] if reverse else _CF4_WEIGHTS
    updated = np.asarray(state, complex)
    tr0 = complex(op0.trace)
    tr1 = complex(op1.trace)
    for weight0, weight1 in pairs:
        operator = span * (weight0 * op0 + weight1 * op1)
        traceA = span * (weight0 * tr0 + weight1 * tr1)
        updated = expm_multiply(operator, updated, traceA=traceA)
        if not np.isfinite(updated).all():
            raise ArithmeticError('nonfinite CF4 state')
    return updated


def evolve_ks_difference_envelope(
        source, initial_columns, family, z_grid, target_z, mass, angular, rho_up, *,
        axial_support, rtol, atol, max_step, tangents='all', rho_sigma=E.RHO_SIGMA,
        integrator='dop853', step_control='joint', tangent_rtol=None, tangent_atol=None):
    """Evolve a fixed upstream preparation without differentiating its carrier.

    Array dimensions and the source/history contract are those of
    evolve_ks_source_envelope.  Integration and continuum errors remain OPEN.
    """
    source = E._require_source(source)
    if float(rho_sigma) != E.RHO_SIGMA:
        raise ValueError('incoming surface is the owned rho_sigma=1 slice')
    rho_up, mass, angular = map(float, (rho_up, mass, angular))
    if not np.isfinite([rho_up, mass, angular]).all() or rho_up < E.RHO_UP_MIN or mass < 0:
        raise ValueError('finite mass>=0, angular and rho_up>=1.03 required')
    rtol, atol, max_step = (E._solver_option(v, n) for v, n in
                          zip((rtol, atol, max_step), ('rtol', 'atol', 'max_step')))
    if tangents not in ('all', 'zero'):
        raise ValueError("tangents must be 'all' or 'zero'")
    if integrator not in ('dop853', 'cf4'):
        raise ValueError("integrator must be 'dop853' or 'cf4'")
    if step_control not in ('joint', 'primal'):
        raise ValueError("step_control must be 'joint' or 'primal'")
    if (tangent_rtol is None) != (tangent_atol is None):
        raise ValueError('tangent_rtol and tangent_atol must be supplied together')
    if step_control == 'joint' and tangent_rtol is not None:
        raise ValueError("tangent accuracy control requires step_control='primal'")
    if tangent_rtol is not None:
        tangent_rtol = E._solver_option(tangent_rtol, 'tangent_rtol')
        tangent_atol = E._solver_option(tangent_atol, 'tangent_atol')
    z, dz, length = E._require_uniform_centered_grid(z_grid)
    target = E._finite_real(target_z, 'target z')
    interval = E.physical_incoming_interval()
    support = E._scalar_pair(axial_support, 'axial_support')
    domain_end = float(z[0] + length)
    if target.ndim != 1 or len(target) == 0 or np.any(np.diff(target) <= 0):
        raise ValueError('nonempty increasing target z required')
    if not (support[0] <= interval[0] < interval[1] <= support[1]):
        raise ValueError('physical I must lie inside axial_support')
    if (np.any(target < interval[0] - 1e-12) or np.any(target > interval[1] + 1e-12)
            or interval[0] < z[0] or interval[1] >= domain_end):
        raise ValueError('target z must lie in physical I inside the computational domain')
    local, metric = E._as_metric(family)
    if local is not None and (not np.allclose(local.interval, interval, rtol=0, atol=1e-12)
                             or not np.allclose(local.support, support, rtol=0, atol=1e-12)):
        raise ValueError('local interval/support must match the supplied domain')
    initial = E._finite_complex(initial_columns, 'upstream columns')
    nsrc = source.covariance.shape[0]
    if initial.shape != (2, nsrc):
        raise ValueError('initial columns must have shape (2, nsrc)')
    w, U = E._sample_axial_profiles(metric.directions, z)
    outside = (z < support[0]) | (z > support[1])
    if max(np.max(abs(w[:, outside]), initial=0), np.max(abs(U[:, outside]), initial=0)) > 1e-12:
        raise ValueError('w/U samples must vanish outside axial_support')
    bound, global_bound = E._radius_lower_bound(local, metric, w, U)
    if bound <= 0:
        raise ValueError('coefficient bound does not keep radius positive')
    inners = tuple(float(d.inner_radius) for d in metric.directions)
    outers = tuple(float(d.outer_radius) for d in metric.directions)
    T1 = E.chart_coordinates(1.)[0]
    s_up = E.chart_coordinates(rho_up)[0] - T1
    if any(abs(s_up) < outer for outer in outers):
        raise ValueError('normal support reaches the fixed upstream preparation')
    distance = E.continuum_speed_distance(rho_up)
    padding = (support[0] - float(z[0]), domain_end - support[1])
    if min(padding) <= 2 * distance:
        raise ValueError('computational padding must exceed 2D')
    energies = np.asarray(source.energies)
    amps = np.asarray(metric.amplitudes)
    ndir = len(amps) if tangents == 'all' else 0
    nz, size_A = len(z), initial.size
    size_D = size_A * nz
    wave = 2 * np.pi * np.fft.fftfreq(nz, d=dz)
    y0 = np.concatenate((initial.ravel(), np.zeros((1 + ndir) * size_D, complex)))

    def unpack(state):
        return (state[:size_A].reshape(2, nsrc),
                state[size_A:size_A + size_D].reshape(2, nsrc, nz),
                state[size_A + size_D:].reshape(ndir, 2, nsrc, nz))

    evaluations = [0]
    stage_evaluations = 0

    def stage_operator(rho):
        nonlocal stage_evaluations
        stage_evaluations += 1
        rho = float(rho)
        beta = float(E.geometry(rho)[0])
        axial = np.sqrt(beta**2 - 1)
        if not np.isfinite(axial) or axial <= 0:
            raise ArithmeticError('strictly trapped a>0 required')
        s = E.chart_coordinates(rho)[0] - T1
        chi = np.array([float(E._plateau(np.asarray([s]), i, o)[0])
                        for i, o in zip(inners, outers)])
        basis = chi[:, None] * (s * w + s**3 * U / 6)
        delta_r = amps @ basis
        r_ref = np.sqrt(1 + rho**2)
        if not np.isfinite(delta_r).all() or np.any(r_ref + delta_r <= 0):
            raise ArithmeticError('actual radius must stay finite and positive')
        return KSDifferenceOperator(
            E.ks_generator(rho, energies, mass, angular), axial, angular, r_ref,
            delta_r, basis[:ndir], wave, evaluations=evaluations)

    def rhs(rho, state):
        unpack  # kept in this closure for primal_block_length
        return stage_operator(rho).matvec(np.asarray(state, complex))

    primal_error_norm = None
    tangent_error_norm = None
    if integrator == 'dop853':
        solver_kwargs = dict(rtol=rtol, atol=atol, max_step=max_step)
        if step_control == 'primal':
            from .nsc_ks_primal_step_control import PrimalBlockDOP853
            solution = PrimalBlockDOP853(
                rhs, rho_up, y0, E.RHO_SIGMA, n_primal=size_A + size_D,
                tangent_rtol=tangent_rtol, tangent_atol=tangent_atol,
                tangent_block_size=size_D, **solver_kwargs)
        else:
            solution = DOP853(rhs, rho_up, y0, E.RHO_SIGMA, **solver_kwargs)
        steps = 0
        while solution.status == 'running':
            message = solution.step()
            if solution.status == 'failed':
                raise ArithmeticError(message)
            steps += 1
        final = solution.y
        nfev = int(solution.nfev)
        primal_error_norm = getattr(solution, 'primal_error_norm', None)
        tangent_error_norm = getattr(solution, 'tangent_error_norm', None)
    else:
        final = np.asarray(y0, complex)
        rho = float(rho_up)
        steps = 0
        while rho > float(E.RHO_SIGMA):
            nxt = max(float(E.RHO_SIGMA), rho - max_step)
            span, nodes = _cf4_gauss_nodes(rho, nxt)
            final = _cf4_step(stage_operator(nodes[0]), stage_operator(nodes[1]), final, span)
            rho = nxt
            steps += 1
        nfev = int(evaluations[0])
    A, D, Y = unpack(final)
    D_t = np.moveaxis(E.trigonometric_polynomial(D, z, target), -1, 0)
    Dz_t = np.moveaxis(E.trigonometric_polynomial(D, z, target, derivative=True), -1, 0)
    if ndir:
        Y_t = np.moveaxis(E.trigonometric_polynomial(Y, z, target), -1, 1)
        Yz_t = np.moveaxis(E.trigonometric_polynomial(Y, z, target, derivative=True), -1, 1)
    else:
        Y_t = np.zeros((0, len(target), 2, nsrc), complex)
        Yz_t = np.zeros_like(Y_t)
    total = A[None] + D_t
    phase = np.exp(-1j * energies * target[:, None])[:, None, :]
    prep = E._fixed_preparation_digest(source, initial, mass, angular, rho_up, E.RHO_SIGMA)
    binding = E.KSEnvelopeBinding(z, w, U, tuple(amps), inners, outers, support,
                                  rho_up, E.RHO_SIGMA, rtol, atol, max_step)
    initial_norm = np.linalg.norm(initial) * np.sqrt(nz)
    final_norm = np.linalg.norm(A[:, :, None] + D)
    diagnostics = {
        'field_representation': 'joint homogeneous reference plus difference',
        'function_evaluations': nfev, 'accepted_steps': steps,
        'time_integrator': integrator,
        'step_control': step_control,
        'n_primal': size_A + size_D,
        'tangent_directions': ndir,
        'tangent_rtol': tangent_rtol,
        'tangent_atol': tangent_atol,
        'primal_error_norm': primal_error_norm,
        'tangent_error_norm': tangent_error_norm,
        'stage_coefficient_evaluations': stage_evaluations,
        'analytic_operator_trace': integrator == 'cf4',
        'fixed_preparation_digest': prep, 'continuum_speed_distance': distance,
        'padding': padding, 'computational_count': nz,
        'radius_lower_bound': bound if global_bound else None,
        'global_radius_bound_supplied': global_bound,
        'computational_norm_residual': float(abs(final_norm - initial_norm)),
        'stress_drift_subtracted': False, 'reference_solves_same_generator': True,
        'reference_z_derivative_exactly_zero': True, 'metric_timestep': False,
    }
    return KSDifferenceIncoming(
        target, phase * total, phase[None] * Y_t,
        phase * (Dz_t - 1j * energies * total),
        phase[None] * (Yz_t - 1j * energies * Y_t),
        source.covariance, source.column_weights, source.energies, mass, angular,
        initial, rho_up, binding, prep, E._coverage(nsrc, len(np.unique(energies))),
        diagnostics, A, D_t, Dz_t)
