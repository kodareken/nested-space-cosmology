"""Focused checks for the new massive representation and information gate."""
import numpy as np
import pytest

from recursive_horizons.nsc_massive_signed_preparation import (
    exact_pairing_identities,constant_single_block_obstruction,map_residuals,
    mixed_pointwise_probe,paired_seed_data_witness,require_physical_pg_preparation,
)


def test_exact_pairing_and_mixed_single_block_obstruction():
    assert exact_pairing_identities()['single_block_exact_kernel_dimension']==0
    assert constant_single_block_obstruction()['kernel_dimension']==0
    for mass,angular in ((0.,np.sqrt(5)),(np.pi/2,0.),(np.pi/2,np.sqrt(5))):
        assert max(map_residuals(mass,angular).values())<3e-11
    p=mixed_pointwise_probe(np.pi/2,np.sqrt(5))
    assert p['pointwise_potential_flip_residual']<3e-11
    assert p['fixed_seed_axis_potential_mismatch']>.1
    assert p['required_basis_connection_norm']>.1


def test_seed_marginal_does_not_determine_angular_partner():
    # Rank-one case: missing Y persists even when CAR forces many Z entries.
    r=paired_seed_data_witness(np.diag([1.,0.]))
    assert r['same_observed_marginal_residual']==0
    assert r['known_charge_paired_marginal_residual']==0
    assert r['witness_CAR_lower_violation']==r['witness_CAR_upper_violation']==0
    assert r['unobserved_negative_same_angular_difference']>1.
    assert r['physical_horizon_preparation_satisfied'] is None
    assert not r['witness_adopted_as_spatial_state']


@pytest.mark.parametrize('input',[
    {'seed_covariance':np.eye(2)/2}, {'lll_map':'flow'}, {'multiplicity':4}, {},
])
def test_fake_physical_preparation_is_rejected(input):
    with pytest.raises(ValueError):require_physical_pg_preparation(**input)
