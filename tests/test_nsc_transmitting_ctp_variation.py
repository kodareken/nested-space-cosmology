"""Checks of the new action differential, not of the old NSC source records."""
import numpy as np
import pytest
from scipy.linalg import expm

from recursive_horizons.nsc_influence import influence
from recursive_horizons.nsc_transmitting_boundary_history import EndpointVariation
from recursive_horizons.nsc_transmitting_ctp_variation import (
    EndpointBranchJets, ctp_first_variation, endpoint_ctp_pullback,
)


def test_noncommuting_two_branch_differential_agrees_with_existing_action():
    c = np.array([[0.7, 0.1j], [-0.1j, 0.2]])
    x = np.array([[0.3, 0.2j], [-0.2j, -0.6]])
    y = np.array([[-0.1, 0.4], [0.4, 0.2]])
    up, um = expm(0.17j*x), expm(-0.23j*y)
    result = ctp_first_variation(c, up, um, 1j*x@up, -1j*y@um)
    eps = 2e-5
    def action(s):
        return influence(c, expm(1j*s*x)@up, expm(-1j*s*y)@um)["principal_action"]
    numeric = (action(eps)-action(-eps))/(2*eps)
    assert result["derivative"] == pytest.approx(numeric, abs=3e-9)


def test_common_and_difference_variations_are_distinct():
    c = np.diag([0.8, 0.2])
    u = np.eye(2);x = np.diag([1., 0.])
    common = ctp_first_variation(c, u, u, 1j*x, 1j*x)
    relative = ctp_first_variation(c, u, u, 0.5j*x, -0.5j*x)
    assert common["derivative"] == 0
    assert relative["derivative"] == pytest.approx(0.8)


def test_nonunitary_tangent_rejected():
    with pytest.raises(ValueError, match="tangent-compatible"):
        ctp_first_variation(np.eye(2)*0.3, np.eye(2), np.eye(2), np.eye(2), np.zeros((2, 2)))


def test_missing_physical_branch_jets_cannot_be_zero_derivatives():
    basis = EndpointVariation("locked-basis", ((0., 0.),)*4)
    result = endpoint_ctp_pullback(np.eye(2)*0.3, np.eye(2), np.eye(2), basis, None)
    assert result["status"] == "OPEN"
    assert result["gaussian_endpoint_covector"] is None


def test_pullback_uses_difference_coordinate_and_original_endpoint_signs():
    basis = EndpointVariation("control-basis", ((0., 0.),)*4)
    # Test input only: these dimensionless endpoint jets are not an NSC history.
    coefficients = np.arange(1, 9).reshape(4, 2)*np.array([1., -1.])
    plus = 0.5j*coefficients[:, :, None, None]*np.diag([1., 0.])
    jets = EndpointBranchJets(basis.basis_id, plus, -plus, "test-only algebraic jets")
    result = endpoint_ctp_pullback(np.diag([0.8, 0.2]), np.eye(2), np.eye(2), basis, jets)
    assert np.asarray(list(result["gaussian_endpoint_covector"].values())) == pytest.approx(0.8*coefficients)
    assert result["two_sided_physical_match"] is None
    wrong = EndpointBranchJets("different-domain", plus, -plus, "test")
    with pytest.raises(ValueError, match="matching KS"):
        endpoint_ctp_pullback(np.diag([0.8, 0.2]), np.eye(2), np.eye(2), basis, wrong)
