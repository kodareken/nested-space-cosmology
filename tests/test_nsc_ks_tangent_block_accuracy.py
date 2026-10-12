"""Separate retarded directions cannot dilute one another's adaptive error."""
import numpy as np

from recursive_horizons.nsc_ks_primal_step_control import PrimalBlockDOP853


def solve(extra_zeros):
    initial = np.zeros(2 + extra_zeros, complex)
    initial[:2] = 1

    def rhs(_t, y):
        result = np.zeros_like(y)
        result[0] = -y[0]
        result[1] = 100j * y[1]
        return result

    solver = PrimalBlockDOP853(
        rhs, 0, initial, 0.1, n_primal=1, rtol=1e-3, atol=1e-8,
        tangent_rtol=1e-10, tangent_atol=1e-12, tangent_block_size=1)
    times = []
    while solver.status == "running":
        solver.step()
        times.append(solver.t)
    assert solver.status == "finished"
    return np.array(times), solver.y


def test_zero_padding_does_not_relax_nonzero_tangent_accuracy():
    times, value = solve(0)
    padded_times, padded_value = solve(64)
    assert len(times) == len(padded_times)
    # BLAS reduction layout can shift the adaptive times slightly. The
    # workload and local step scales must stay fixed, and the physical
    # result below has its own much tighter comparison.
    np.testing.assert_allclose(
        np.diff(np.r_[0, times]), np.diff(np.r_[0, padded_times]),
        rtol=5e-8, atol=1e-12)
    np.testing.assert_allclose(value, padded_value[:2], rtol=0, atol=2e-12)
    assert abs(value[1] - np.exp(10j)) < 1e-10
    assert np.count_nonzero(padded_value[2:]) == 0
