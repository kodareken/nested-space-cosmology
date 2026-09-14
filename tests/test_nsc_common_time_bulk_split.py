"""Independent bulk support and transmitted operator-domain checks."""
import numpy as np
import pytest

from recursive_horizons.nsc_causal_common import LinkHistory
from recursive_horizons.nsc_common_time_bulk_split import (
    CommonTimeBulkSplit, BulkOperatorDomainError, require_common_time_bulk_split,
)
from recursive_horizons.nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain, I2, S1, S2


def test_unequal_independent_bulk_parts_recombine_without_duplication():
    domain = CommonTimeBulkSplit()
    maps = domain.projectors([-3., -1., 0.5, 1., 3.])
    pp, pc, reorder = maps['P_parent'], maps['P_child'], maps['canonical_reorder']
    assert np.trace(pp) == 6 and np.trace(pc) == 4
    assert np.array_equal(pp@pc, np.zeros((10, 10)))
    assert np.array_equal(pp+pc, np.eye(10))
    assert np.array_equal(reorder@reorder.T, np.eye(10))


def test_given_B_and_complete_trace_pair_are_not_bulk_projectors():
    for wrong in (LinkHistory(np.zeros((1, 1, 1))), TransmittingDiracSeamDomain(1., 1., 2., 1.)):
        with pytest.raises(BulkOperatorDomainError): require_common_time_bulk_split(wrong)


def test_seam_sample_is_not_an_independent_duplicate_bulk_node():
    with pytest.raises(BulkOperatorDomainError):CommonTimeBulkSplit().projectors([-1., 0., 1.])


def test_current_phase_preserves_inherited_massless_PG_operator():
    domain = CommonTimeBulkSplit()
    phase = (I2+1j*S2)/np.sqrt(2)
    local = domain.local_potential(1., 2., 0., 3.)
    assert np.max(abs(phase@local@phase.conj().T-1.5*S1)) < 1e-14
    assert np.max(abs(phase@S2@phase.conj().T-S2)) < 1e-14


def test_T_normal_measure_matches_common_time_flux_with_nonunit_q():
    seam = TransmittingDiracSeamDomain(1.2, 0.8, 2., 1.1)
    w = np.array([0.2, 0.3])
    _, _, gram_T = seam.trace_map(w)
    v = CommonTimeBulkSplit.principal(seam.lapse, seam.radial_scale, seam.shift)
    assert np.max(abs(gram_T+seam.radial_scale*seam.radius**2*w[:,None,None]*v)) < 1e-14


def test_domain_delta_is_not_an_ordinary_link_or_a_covariance_map():
    domain = CommonTimeBulkSplit()
    defect = domain.projection_domain_defect(N=1., q_PG=1., beta=2.)
    assert defect['smallest_singular_value'] == pytest.approx(1.)
    assert defect['largest_singular_value'] == pytest.approx(3.)
    assert defect['ordinary_link_matrix'] is None
    with pytest.raises(BulkOperatorDomainError, match='does not preserve'):
        domain.ordinary_hamiltonian_link()
    with pytest.raises(BulkOperatorDomainError, match='Cauchy-state map'):
        domain.attach_seed_trace_covariance(np.eye(4))
