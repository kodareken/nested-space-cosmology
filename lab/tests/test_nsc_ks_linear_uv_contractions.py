"""The next diagonal coefficient is necessary for an E-linear observable."""
from recursive_horizons.nsc_ks_linear_uv_contractions import contraction_identity


def test_projector_pair_cancels_both_constraint_coefficients_per_angular_sign():
    result = contraction_identity()
    residuals = result['residuals']
    assert residuals['projector_pair_E_minus3'] == ['0'] * 4
    assert all(residuals[name] == '0' for name in ('N_E_minus1', 'N_E_minus2', 'beta_E_minus1', 'beta_E_minus2'))
    assert result['omitted_diagonal_false_N_coefficient'] != '0'
    assert result['angular_pairing_needed_for_these_zeros'] is False
    assert result['scope']['finite_history_bound'] is False
    assert result['scope']['numeric_tail_constant'] is None
