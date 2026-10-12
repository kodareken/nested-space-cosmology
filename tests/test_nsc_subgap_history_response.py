"""Focused contracts for the new subgap/history connection."""
import numpy as np
import pytest

from recursive_horizons.nsc_subgap_history_response import (
    trapped_subgap_contractions, AnalyticResponsePanel, integrate_subgap_response,
)
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_transmitting_history_jets import relative_branch_ctp_control


def test_signed_source_matches_existing_gaussian_differential():
    E, kappa = 1.2, .23832579963401956
    B = np.array([[[.03, .002+.004j], [.002-.004j, -.01]]]*4)
    R = np.exp(.7j)
    smooth, reflection = trapped_subgap_contractions(E, kappa, B)
    expected = -(smooth.real+2*np.real(R*reflection))/np.pi
    D = np.diag([R, 1.])
    A = np.array([D.conj().T@(-1j*b)@D for b in B])
    kernels = np.zeros((4, 6, 6), complex)
    kernels[:, :2, :2] = A.conj()
    kernels[:, 3:5, 3:5] = A
    C = np.array([source_covariance(e, kappa, 3.973074368754331, np.pi/2) for e in (-E, E)])
    P = np.array([np.diag([1., 1., 0.])]*2)
    rows = relative_branch_ctp_control(np.array([-E, E]), np.ones(2), C, P, kernels)
    assert np.max(abs(expected-np.array([row['derivative'] for row in rows]))) < 3e-14


def test_analytic_duality_is_not_same_complex_energy_adjoint():
    x = np.cos(np.pi*np.arange(9)/8)
    E = 1.25+.25*x
    def exact(z):
        B = np.zeros((4, 2, 2), complex)
        B[:, 0, 0] = z*z
        B[:, 1, 1] = -z
        B[:, 0, 1] = 1j*z
        B[:, 1, 0] = -1j*z
        return B
    panel = AnalyticResponsePanel(E, np.array([exact(e) for e in E]), 1., 1.5)
    z = 1.2+.06j
    assert np.max(abs(panel(z)-exact(z))) < 1e-12
    assert np.max(abs(panel(z)-panel(z.conjugate()).swapaxes(-1, -2).conj())) < 1e-12
    assert np.max(abs(panel(z)-panel(z).swapaxes(-1, -2).conj())) > .1


def test_pole_strip_and_unresolved_panel_rejected():
    E = np.linspace(1., 1.5, 5)
    panel = AnalyticResponsePanel(E, np.zeros((5, 4, 2, 2)), 1., 1.5)
    with pytest.raises(ValueError, match='strip'):
        integrate_subgap_response(panel, .2, lambda z: 1., height=.1)
    with pytest.raises(ValueError, match='samples'):
        AnalyticResponsePanel(E[:2], np.zeros((2, 4, 2, 2)), 1., 1.5)
