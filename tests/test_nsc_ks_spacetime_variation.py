"""Checks of the new spacetime chain rule, not historical generators."""
import sys
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from recursive_horizons.nsc_ks_spacetime_variation import (
    CompactKSHarmonic, reference_metric_jacobian_jet,
    pullback_variation_jet, harmonic_dyson_derivative,
)


def test_chart_chain_rule_includes_axial_phase_and_metric_jacobian_derivative():
    direction = CompactKSHarmonic(.4)
    e = np.array([.3, -.2, .1, .4])
    def field(rho):
        return reference_metric_jacobian_jet(rho)[0]@e*direction.data(rho)[0]
    for rho in (-.6, 0., .6):
        f, fT, fz, _ = direction.data(rho)
        value, derivative = pullback_variation_jet(rho, e*f, e*fT, e*fz)
        h = 1e-5
        finite = (-field(rho+2*h)+8*field(rho+h)-8*field(rho-h)+field(rho-2*h))/(12*h)
        assert np.max(abs(value-field(rho))) < 3e-11
        assert np.max(abs(derivative-finite)) < 3e-8
        # Identifying KS z with PG tau omits a nonzero radial phase term.
        _, missing_phase = pullback_variation_jet(rho, e*f, e*fT, np.zeros(4))
        assert np.linalg.norm(derivative-missing_phase) > 1e-3


def test_dyson_integral_matches_direct_time_quadrature_and_unitary_tangent():
    # Finite algebra check only: no matrix here is assigned to an NSC link.
    E = np.array([-.7, .2, .9])
    M = np.array([[1., .2j, .3], [.1, -.4, .2], [.1j, -.3j, .6]])
    ti, tf, omega = .2, .6, .4
    U, du, generator = harmonic_dyson_derivative(E, M, omega, ti, tf)
    nodes, weights = leggauss(32)
    direct = np.zeros((3, 3), complex)
    for time, weight in zip((ti+tf)/2+(tf-ti)*nodes/2, weights*(tf-ti)/2):
        h1 = .5*(np.exp(-1j*omega*time)*M+np.exp(1j*omega*time)*M.conj().T)
        direct += -1j*weight*np.exp(-1j*E*(tf-time))[:, None]*h1*np.exp(-1j*E*(time-ti))[None, :]
    assert np.max(abs(du-direct)) < 3e-11
    assert np.max(abs(du.conj().T@U+U.conj().T@du)) < 3e-11
    assert np.max(abs(generator-generator.conj().T)) < 3e-11
    with pytest.raises(ValueError, match='ordered caller'):
        harmonic_dyson_derivative(E, M, omega, tf, ti)


def test_no_endpoint_or_history_selection_is_returned_by_the_direction():
    direction = CompactKSHarmonic(.4)
    assert not hasattr(direction, 'physical_history')
    assert not hasattr(direction, 'as_endpoint_branch_jets')
    with pytest.raises(ValueError, match='four finite'):
        pullback_variation_jet(0., [0., 1.], [0., 1.], [0., 1.])
