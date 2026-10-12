"""Adjoint, analytic trace and CF4 order for the augmented difference operator."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm, norm
from scipy.sparse.linalg import LinearOperator, expm_multiply, onenormest

from recursive_horizons import nsc_ks_source_envelope as E
from recursive_horizons.nsc_ks_difference_envelope import (
    KSDifferenceOperator, _cf4_gauss_nodes, _cf4_step, difference_operator_trace,
    difference_rhs, difference_rhs_adjoint, evolve_ks_difference_envelope)
from recursive_horizons.nsc_ks_primal_step_control import evolve_primal_controlled
from tests.test_nsc_ks_difference_envelope import evolve, inputs


def _synthetic(nsrc=2, nz=6, ndir=1, seed=4, traced=False):
    rng = np.random.default_rng(seed)
    G = rng.normal(size=(nsrc, 2, 2)) + 1j * rng.normal(size=(nsrc, 2, 2))
    if traced:
        G = G + np.eye(2)[None] * (0.4 - 0.2j)
    axial, angular, r_ref = 1.2, np.sqrt(5), 1.4
    delta_r = 0.03 * rng.normal(size=nz)
    basis = np.zeros((0, nz)) if ndir == 0 else 0.04 * rng.normal(size=(ndir, nz))
    wave = 2 * np.pi * np.fft.fftfreq(nz, d=0.05)
    return G, axial, angular, r_ref, delta_r, basis, wave


def _pack(A, D, Y):
    return np.concatenate([np.asarray(A, complex).ravel(),
                           np.asarray(D, complex).ravel(),
                           np.asarray(Y, complex).ravel()])


def _random_state(nsrc, nz, ndir, rng):
    A = rng.normal(size=(2, nsrc)) + 1j * rng.normal(size=(2, nsrc))
    D = rng.normal(size=(2, nsrc, nz)) + 1j * rng.normal(size=(2, nsrc, nz))
    Y = np.zeros((0, 2, nsrc, nz), complex) if ndir == 0 else (
        rng.normal(size=(ndir, 2, nsrc, nz)) + 1j * rng.normal(size=(ndir, 2, nsrc, nz)))
    return A, D, Y


class _MatrixOp(LinearOperator):
    def __init__(self, matrix):
        matrix = np.asarray(matrix, complex)
        super().__init__(complex, matrix.shape)
        self._A = matrix
        self.trace = np.trace(matrix)

    def _matvec(self, x):
        return self._A @ x

    def _rmatvec(self, x):
        return self._A.conj().T @ x

    def _matmat(self, X):
        return self._A @ X

    def _rmatmat(self, X):
        return self._A.conj().T @ X


def test_all_ones_probe_cancels_on_shift_but_one_norm_is_two():
    block = 1j * np.array([[1.0, -1.0], [-1.0, 1.0]])
    matrix = np.zeros((4, 4), complex)
    matrix[:2, :2] = block
    ones = np.ones(4, complex)
    bound = float(4 * norm(matrix @ ones, 1) / max(norm(ones, 1), 1e-30))
    assert bound == 0.0
    assert abs(norm(matrix, 1) - 2.0) < 1e-15
    probe = np.array([1.0, 0.0, 0.0, 0.0], complex)

    def fake_rmatvec(vector):
        return bound * np.asarray(vector, complex)

    fake = LinearOperator((4, 4), matvec=lambda v: matrix @ v, rmatvec=fake_rmatvec,
                          dtype=complex)
    true = _MatrixOp(matrix)
    np.testing.assert_array_equal(fake.rmatvec(probe), 0)
    np.testing.assert_allclose(true.rmatvec(probe), matrix.conj().T @ probe, atol=0, rtol=0)
    assert np.max(np.abs(true.rmatvec(probe))) > 0
    op = KSDifferenceOperator(*_synthetic(nsrc=1, nz=5, ndir=1, seed=6))
    dense = op.matmat(np.eye(op.shape[0], dtype=complex))
    exact = float(norm(dense, 1))
    estimated = float(onenormest(op))
    assert estimated > exact / 3
    assert estimated <= exact * 1.05


def test_random_complex_adjoint_identity_all_and_zero_tangents():
    rng = np.random.default_rng(21)
    for ndir in (0, 1, 2):
        args = _synthetic(ndir=ndir, seed=10 + ndir, traced=True)
        op = KSDifferenceOperator(*args)
        nsrc, nz = args[0].shape[0], args[4].shape[0]
        for _ in range(3):
            x = _pack(*_random_state(nsrc, nz, ndir, rng))
            z = _pack(*_random_state(nsrc, nz, ndir, rng))
            Tx = op.matvec(x)
            THz = op.rmatvec(z)
            np.testing.assert_allclose(np.vdot(z, Tx), np.vdot(THz, x), atol=2e-12, rtol=0)
        X = rng.normal(size=(op.shape[0], 3)) + 1j * rng.normal(size=(op.shape[0], 3))
        np.testing.assert_allclose(
            op.matmat(X), np.column_stack([op.matvec(X[:, j]) for j in range(3)]), atol=2e-12, rtol=0)
        np.testing.assert_allclose(
            op.rmatmat(X), np.column_stack([op.rmatvec(X[:, j]) for j in range(3)]), atol=2e-12, rtol=0)
        column = X[:, :1]
        assert op.matvec(column).shape == (op.shape[0], 1)
        assert op.rmatvec(column).shape == (op.shape[0], 1)


def test_adjoint_includes_nonnormal_cross_blocks():
    rng = np.random.default_rng(3)
    args = _synthetic(ndir=2, seed=8)
    op = KSDifferenceOperator(*args)
    nsrc, nz = args[0].shape[0], args[4].shape[0]
    A, D, Y = _random_state(nsrc, nz, 2, rng)
    packed = _pack(A, D, Y)
    dense = op.matmat(np.eye(op.shape[0], dtype=complex))
    gram = dense.conj().T @ dense - dense @ dense.conj().T
    assert norm(gram, 2) > 0.1
    Ad, Dd, Yd = difference_rhs_adjoint(A, D, Y, *args)
    np.testing.assert_allclose(_pack(Ad, Dd, Yd), op.rmatvec(packed), atol=2e-12, rtol=0)
    Pd, Qd, Rd = difference_rhs(A, D, Y, *args)
    np.testing.assert_allclose(_pack(Pd, Qd, Rd), op.matvec(packed), atol=2e-12, rtol=0)
    zA, zD, zY = _random_state(nsrc, nz, 2, rng)
    z = _pack(zA, zD, zY)
    Zd = difference_rhs_adjoint(zA, zD, zY, *args)
    np.testing.assert_allclose(np.vdot(z, _pack(Pd, Qd, Rd)), np.vdot(_pack(*Zd), packed),
                               atol=2e-12, rtol=0)


def test_analytic_trace_uses_G_multiplicities_without_assuming_zero():
    nsrc, nz, ndir = 2, 5, 3
    G = np.zeros((nsrc, 2, 2), complex)
    G[0] = np.diag([1.5 + 0.2j, -0.4])
    G[1] = np.array([[0.3j, 1.0], [0.0, 2.0]])
    g_trace = np.trace(G, axis1=-2, axis2=-1).sum()
    assert abs(g_trace) > 1
    wave = 2 * np.pi * np.fft.fftfreq(nz, d=0.1)
    correction = 0.2 * np.linspace(-1, 1, nz)
    got = difference_operator_trace(G, nsrc, nz, ndir, 1.1, correction, wave)
    expected = g_trace * (1 + nz * (1 + ndir))
    np.testing.assert_allclose(got, expected, atol=1e-13, rtol=0)
    op = KSDifferenceOperator(G, 1.1, 0.5, 1.3, 0.01 * np.ones(nz),
                              0.02 * np.ones((ndir, nz)), wave)
    dense = op.matmat(np.eye(op.shape[0], dtype=complex))
    np.testing.assert_allclose(op.trace, np.trace(dense), atol=1e-10, rtol=0)
    physical = E.ks_generator(1.02, np.array([0.7, -0.4]), np.pi / 2, np.sqrt(5))
    phys = difference_operator_trace(physical, 2, nz, ndir, 1.1, correction, wave)
    np.testing.assert_allclose(phys, 0.0, atol=1e-12, rtol=0)


def test_dense_expm_matches_linear_operator():
    args = _synthetic(nsrc=1, nz=4, ndir=1, seed=2)
    op = KSDifferenceOperator(*args)
    vector = np.linspace(0.2, 1.1, op.shape[0]) + 1j * np.linspace(-0.4, 0.3, op.shape[0])
    dense = op.matmat(np.eye(op.shape[0], dtype=complex))
    dt = 0.03
    np.testing.assert_allclose(
        expm_multiply(dt * op, vector, traceA=dt * op.trace),
        expm(dt * dense) @ vector, atol=2e-12, rtol=0)
    constant = _cf4_step(op, op, vector, dt)
    np.testing.assert_allclose(constant, expm(dt * dense) @ vector, atol=2e-12, rtol=0)


def _pauli_flow(t):
    sx = np.array([[0.0, 1.0], [1.0, 0.0]])
    sz = np.array([[1.0, 0.0], [0.0, -1.0]])
    return 1j * ((1.2 + t) * sz + 0.8 * t * sx)


def _cf4_integrate(fun, y0, t0, t1, steps, reverse=False, absolute_span=False):
    y = np.asarray(y0, complex)
    nodes = np.linspace(t0, t1, steps + 1)
    for left, right in zip(nodes[:-1], nodes[1:]):
        span, stages = _cf4_gauss_nodes(left, right)
        if absolute_span:
            span = abs(span)
        y = _cf4_step(_MatrixOp(fun(stages[0])), _MatrixOp(fun(stages[1])), y, span, reverse=reverse)
    return y


def test_cf4_noncommuting_order_and_backward_integration():
    y0 = np.array([0.4 + 0.1j, -0.2 + 0.3j])
    reference = solve_ivp(
        lambda t, y: _pauli_flow(t) @ y, (0.0, 0.6), y0, method='DOP853',
        rtol=1e-12, atol=1e-14, dense_output=True)
    assert reference.success
    y_end = reference.sol(0.6)
    errors = []
    for steps in (4, 8, 16):
        got = _cf4_integrate(_pauli_flow, y0, 0.0, 0.6, steps)
        errors.append(float(np.max(np.abs(got - y_end))))
    assert errors[1] < errors[0] / 8, errors
    assert errors[2] < errors[1] / 8, errors
    backward = _cf4_integrate(_pauli_flow, y_end, 0.6, 0.0, 16)
    np.testing.assert_allclose(backward, y0, atol=3e-8, rtol=0)
    wrong_sign = _cf4_integrate(_pauli_flow, y_end, 0.6, 0.0, 16, absolute_span=True)
    assert np.max(np.abs(wrong_sign - y0)) > 1e-2
    reversed_errors = []
    for steps in (4, 8, 16):
        got = _cf4_integrate(_pauli_flow, y0, 0.0, 0.6, steps, reverse=True)
        reversed_errors.append(float(np.max(np.abs(got - y_end))))
    assert reversed_errors[1] > reversed_errors[0] / 8
    paper = _cf4_integrate(_pauli_flow, y0, 0.0, 0.6, 8)
    round_trip = _cf4_integrate(_pauli_flow, paper, 0.6, 0.0, 8)
    np.testing.assert_allclose(round_trip, y0, atol=2e-12, rtol=0)


def test_cf4_finite_difference_fields_and_tangents():
    kwargs = dict(integrator='cf4', max_step=0.01, rtol=1e-8, atol=1e-10)
    source, initial, family = inputs()
    grid = E.computational_z_grid(32)
    target = np.linspace(*E.physical_incoming_interval(), 9)
    common = dict(axial_support=E.usual_axial_support(), **kwargs)
    base = evolve_ks_difference_envelope(
        source, initial, family, grid, target, np.pi / 2, np.sqrt(5), 1.03, **common)
    h = 1e-4
    plus = evolve_ks_difference_envelope(
        *inputs(0.002 + h), grid, target, np.pi / 2, np.sqrt(5), 1.03,
        tangents='zero', **common)
    minus = evolve_ks_difference_envelope(
        *inputs(0.002 - h), grid, target, np.pi / 2, np.sqrt(5), 1.03,
        tangents='zero', **common)
    for name, derivative in (('columns', 'column_tangents'), ('axial_columns', 'axial_tangents')):
        fd = (getattr(plus, name) - getattr(minus, name)) / (2 * h)
        np.testing.assert_allclose(fd, getattr(base, derivative)[0], atol=4e-6, rtol=0)
    assert base.diagnostics['time_integrator'] == 'cf4'
    assert base.diagnostics['analytic_operator_trace'] is True
    assert base.diagnostics['stage_coefficient_evaluations'] == 2 * base.diagnostics['accepted_steps']
    assert base.diagnostics['function_evaluations'] > base.diagnostics['stage_coefficient_evaluations']
    assert base.fixed_preparation_digest == plus.fixed_preparation_digest


def test_zero_tangent_directions_do_not_change_primal():
    all_t = evolve(owner=evolve_primal_controlled, tangents='all')
    zero = evolve(owner=evolve_primal_controlled, tangents='zero')
    np.testing.assert_allclose(all_t.reference_amplitudes, zero.reference_amplitudes, atol=1e-14, rtol=0)
    np.testing.assert_allclose(all_t.envelope_difference, zero.envelope_difference, atol=1e-14, rtol=0)
    assert all_t.diagnostics['accepted_steps'] == zero.diagnostics['accepted_steps']
    assert all_t.diagnostics['step_control'] == 'primal'
    assert zero.diagnostics['tangent_directions'] == 0
    assert all_t.diagnostics['tangent_directions'] == 1
