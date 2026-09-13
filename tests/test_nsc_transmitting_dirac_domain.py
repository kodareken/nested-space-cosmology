"""Domain-definition checks for the same-surface trace, not a new dynamics test."""
import numpy as np
import pytest

from recursive_horizons.nsc_causal_common import LinkHistory
from recursive_horizons.nsc_transmitting_dirac_domain import (
    TransmittingDiracSeamDomain, DomainRepresentationError, require_transmitting_domain,
    MODE_TO_CURRENT, I2, S1, S2, S3,
)


def test_rotation_transports_current_and_both_other_axes():
    r = MODE_TO_CURRENT
    assert np.max(abs(r@S3@r.conj().T-S2)) < 1e-14
    assert np.max(abs(r@S2@r.conj().T+S3)) < 1e-14
    assert np.max(abs(r@S1@r.conj().T-S1)) < 1e-14


def test_normal_and_trace_weight_follow_nonunit_lapse_and_radial_scale():
    # Definition-level geometry case, not a new physical NSC solution.
    domain = TransmittingDiracSeamDomain(1.2, 0.8, 2.0, 1.1)
    n = domain.normal_geometry()
    assert n['normal_unit_residual'] < 1e-14
    assert n['normal_tangent_residual'] < 1e-14
    result = domain.evaluate_channel([0.1, 0.7], np.array([0.2*I2, 0.8*I2]))
    assert max(result['residuals'].values()) < 1e-14


def test_null_or_timelike_cut_is_outside_this_domain():
    for shift in (0.5, 1.):
        with pytest.raises(DomainRepresentationError):
            TransmittingDiracSeamDomain(1., 1., shift, 1.)


def test_given_link_matrix_cannot_satisfy_domain_definition():
    with pytest.raises(DomainRepresentationError, match='given B'):
        require_transmitting_domain(LinkHistory(np.zeros((2, 1, 1))))


def test_sewing_does_not_create_independent_room_dofs_or_hopping():
    domain = TransmittingDiracSeamDomain(1., 1., 2., 1.)
    result = domain.evaluate_channel([0.2, 0.3], np.array([0.3*I2, 0.6*I2]))
    assert result['independent_coefficients'] == 4
    assert result['allowed_graph_dimension'] == 4
    assert result['doubled_trace_dimension'] == 8
    with pytest.raises(DomainRepresentationError, match='instantaneous'):
        domain.instantaneous_link_history()
