"""Independent algebraic and ODE controls for correlated difference bounds."""
import numpy as np
import pytest
pytest.importorskip("flint")
from scipy.linalg import expm

from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_difference_error import (
    difference_matter_error, propagate_difference_error)


def number(record, name):
    return float(restored_upper(record[name]))


def test_characteristic_majorant_encloses_independent_positive_coupled_ode():
    # The constant positive system saturates the individual norm inequalities.
    # Its exact exponential is independent of the integrated formula under test.
    kr, k, kz, mint, mzint = 0.2, 0.3, 0.1, 0.04, 0.07
    a0, d0, dz0 = 0.03, 0.02, 0.04
    ra, rd, rdz = 0.01, 0.02, 0.01
    matrix = np.array([[kr, 0, 0, ra], [mint, k, 0, rd],
                       [mzint, kz, k, rdz], [0, 0, 0, 0]])
    exact = expm(matrix) @ np.array([a0, d0, dz0, 1.0])
    bound = propagate_difference_error(
        a0, ra, kr, [d0, dz0], [rd, rdz], k, kz, mint, mzint, 4.0)
    for field, actual in zip(("reference_row_error_upper", "difference_row_error_upper",
                              "difference_axial_envelope_error_upper"), exact[:3]):
        assert number(bound, field) >= actual
    assert bound["source_accuracy_included"] is False
    assert bound["tangent_accuracy_included"] is False


def test_zero_coupling_preserves_exact_zero_difference_despite_reference_error():
    bound = propagate_difference_error(1e-6, 1e-6, 1, [0, 0], [0, 0], 1, 0, 0, 0, 100)
    assert number(bound, "difference_row_error_upper") == 0
    assert number(bound, "difference_axial_envelope_error_upper") == 0
    result = difference_matter_error(3, 0, 7, 0, 1e-6, 0, 2e-6, 0, 2,
        mass=1, absolute_angular=2, axial_lower=0.8, radius_lower=1.3, multiplicity=9)
    assert number(result, "N") == number(result, "beta") == 0


def test_correlated_matter_bound_contains_full_coherent_matrix_perturbations():
    rng = np.random.default_rng(20260925)
    m, ell, scale, radius, mu = 1.7, 2.3, 0.8, 1.3, 7
    s1 = np.array([[0, 1], [1, 0]], complex)
    s3 = np.diag([1, -1])
    potential = m*s3 + ell/radius*s1
    for _ in range(32):
        values = [rng.normal(size=(2, 5))+1j*rng.normal(size=(2, 5)) for i in range(8)]
        a, d, az, dz = values[:4]
        ea, ed, eaz, edz = [v*10**rng.uniform(-8, -2) for v in values[4:]]
        # Dense Hermitian covariance retains off-diagonal coherence.
        base = rng.normal(size=(5, 5))+1j*rng.normal(size=(5, 5))
        cov = (base + base.conj().T)/2
        def stress(a, d, az, dz):
            cross = d@cov@a.conj().T
            density = cross + cross.conj().T + d@cov@d.conj().T
            crossz = dz@cov@a.conj().T + az@cov@d.conj().T + dz@cov@d.conj().T
            momentum = (crossz-crossz.conj().T)/(2j)
            return mu*np.array([np.trace(potential@density)+np.trace(s3@momentum)/scale,
                                np.trace(momentum)]).real
        measured = abs(stress(a+ea, d+ed, az+eaz, dz+edz)-stress(a, d, az, dz))
        norms = [np.linalg.norm(v) for v in (a, d, az, dz, ea, ed, eaz, edz)]
        result = difference_matter_error(*norms, np.linalg.norm(cov, 2),
            mass=m, absolute_angular=ell, axial_lower=scale, radius_lower=radius,
            multiplicity=mu)
        assert np.all(measured <= [number(result, "N"), number(result, "beta")])


@pytest.mark.parametrize("bad", [None, True, -1, float("nan"), float("inf")])
def test_unenclosed_or_invalid_input_is_rejected(bad):
    with pytest.raises(ValueError, match="bound"):
        propagate_difference_error(0, bad, 1, [0, 0], [0, 0], 1, 1, 1, 1, 1)
