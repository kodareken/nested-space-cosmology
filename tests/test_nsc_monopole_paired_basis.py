"""New NSC magnetic normalization only; no old angular/source generators."""
import numpy as np
import pytest

from recursive_horizons.nsc_monopole_paired_basis import (
    orthonormality_residuals,ladder_residuals,gauge_residuals,
    charge_map_residuals,paired_sigma3_lock_residuals,paired_mode_overlap,
    require_charged_spinor_eigenspace,
)


def test_magnetic_pair_is_orthonormal_and_has_declared_eigenvalue():
    assert orthonormality_residuals(n_max=2,n_theta=16,n_phi=16)['max_gram_residual']<3e-11
    assert max(ladder_residuals(n_max=2,n_theta=16).values())<3e-11
    assert np.linalg.norm(paired_mode_overlap(4,1)-np.eye(2))<3e-11


def test_bundle_and_signed_radial_phase_are_compatible():
    assert max(gauge_residuals(n_max=2,n_theta=16,n_phi=16).values())<3e-11
    c=charge_map_residuals(n_max=2)
    assert c['stated_q_positive_residual']<3e-11
    assert c['negative_bundle_physical_residual']<3e-11
    assert paired_sigma3_lock_residuals()['max_eta1_sigma3K_lock_residual']<3e-11
    with pytest.raises(ValueError):require_charged_spinor_eigenspace(degeneracy=12)
