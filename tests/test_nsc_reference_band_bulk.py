"""Exact band-boundary algebra and the actual retained multiplicity contract."""
import copy
import json
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_reference_band_bulk import (
    S1, S2, S3, exact_leading_coefficients, exact_star_trace_coefficients,
    recursion_power_certificate, angular_pairing_certificate, paired_unweighted_bulk_closure,
)
from recursive_horizons.nsc_incoming_state_moments import incoming_group_factor, RADIUS, AXIAL
from recursive_horizons.nsc_spatial_reference_symbol import SIGMA


def retained_channels():
    root = Path(__file__).resolve().parents[1]
    return json.loads((root/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']


def test_polar_frame_and_all_leading_primitives_have_exact_residuals():
    np.testing.assert_array_equal(np.array([S1.tolist(), S2.tolist(), S3.tolist()], complex), SIGMA)
    result = exact_leading_coefficients()
    for row in result['massive_chart'].values():
        assert all(set(values) == {'0'} for values in row['residuals'].values())
        assert row['unpaired_primitive_generically_nonzero']
        # Lambda=0 removes v and all of its metric/coordinate variations.
        for name in ('push_k_primitive_constant', 'band_k_primitive_constant'):
            expression = sp.sympify(row[name], locals={'beta': sp.Symbol('beta')})
            zeros = {symbol: 0 for symbol in expression.free_symbols
                     if symbol.name in ('v_z', 'delta_v_z')}
            assert expression.subs(zeros) == 0
    for row in result['massless_nonzero_angular_chart'].values():
        assert all(set(values) == {'0'} for values in row.values())
    for row in result['LLL']['residuals'].values():
        assert all(set(values) == {'0'} for values in row.values())
    assert result['LLL']['band_exchange_scope'] == 'each fixed nonzero-k chiral chart of the owned canonical reference'
    assert not result['LLL']['distributional_gap_zero_bridge_proved']


def test_star_trace_identity_and_every_recursion_order_include_shift():
    star = exact_star_trace_coefficients()
    assert [row['k_primitive_coefficient'] for row in star] == ['I', '0', '-I/24', '0']
    assert all(set(row['coefficient_residuals']) == {'0'} for row in star)
    rows = recursion_power_certificate()['rows']
    for n, row in enumerate(rows, 1):
        assert row['F_power'] == -n and row['G_power'] == -n-2
        assert row['P_power'] == -n-1 and row['U_and_deltaU_power'] == -n-1
        assert row['D_power'] == -n-1 and row['effective_h_power'] == -n
        assert row['recursion_residuals'] == [0, 0, 0]
        assert row['frame_recursion_residuals'] == [0, 0, 0, 0, 0]
        assert row['k_boundary_power_before_pairing'] == 1-n
    assert not rows[0]['higher_boundary_decays']
    assert all(row['higher_boundary_decays'] for row in rows[1:])
    for row in exact_leading_coefficients()['massive_chart'].values():
        assert 'delta_beta' in row['push_k_primitive_constant']
        assert 'beta' in row['band_k_primitive_constant']


def test_exact_pair_weights_are_the_owned_63_family_trace():
    channels = retained_channels()
    result = angular_pairing_certificate(channels)
    assert (result['group_count'], result['signed_family_count'], result['paired_group_count']) == (33, 63, 30)
    common_measure = 1/(4*np.pi**2*RADIUS**2*AXIAL)
    for channel, row in zip(channels, result['groups']):
        weight = row['weight_per_sign']
        expected = weight['numerator']/weight['denominator']*common_measure
        # Check the existing executable source-owner contract, rather than
        # treating degeneracy as a separate multiplicity for each sign.
        observed = incoming_group_factor(channel, row['angular_signs'])
        np.testing.assert_allclose(observed, expected, rtol=2e-15, atol=0)
        if channel['angular_eigenvalue']:
            assert row['odd_angular_weight_sum'] == '0'
            with pytest.raises(ValueError, match='sign inventory'):
                incoming_group_factor(channel, [1])
    bad = copy.deepcopy(channels); bad[0]['degeneracy'] = -1
    with pytest.raises(ValueError, match='positive integer'):
        angular_pairing_certificate(bad)


def test_bulk_closure_rejects_finite_cutoff_weighted_or_endpoint_promotion():
    domain = dict(full_momentum_line=True, unweighted=True,
                  compact_test_variations=True, smooth_positive_metric_jets=True)
    result = paired_unweighted_bulk_closure(retained_channels(), **domain)
    assert result['bulk_action_gradient'] == dict(N=0, beta=0, a=0, r=0)
    assert result['scope']['paired_unweighted_full_k_bulk_exchange_closed']
    assert not result['scope']['finite_momentum_window_exchange_zero']
    assert not result['scope']['EndpointBranchJets_completed']
    assert not result['scope']['physical_source_spectral_convergence_proved']
    assert not result['scope']['Gamma_rest_assigned']
    assert not result['scope']['LLL_distributional_gap_zero_bridge_proved']
    assert result['scope']['LLL_trace_convention'] == 'sum of the two owned fixed nonzero-k chiral charts; retain the separate c=4 geometric allocation'
    assert set(result['boundary_accounting'].values()) == {'required separately'}
    for name in domain:
        changed = dict(domain); changed[name] = False
        with pytest.raises(ValueError, match='full unweighted'):
            paired_unweighted_bulk_closure(retained_channels(), **changed)
