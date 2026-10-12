"""Finite-history UV order, defect and honest missing-constant controls."""
from recursive_horizons.nsc_ks_finite_history_uv_remainder import (
    CURRENT_HISTORY_IDENTITY,
    QUANTITY_NAMES,
    authenticate_changed_history_uv_bindings,
    changed_history_uv_bindings,
    changed_history_uv_quantities,
    changed_history_uv_successor_report,
    conditional_m4_remainder_majorant,
    current_history_remainder_pilot,
    fermi_thermal_tail_bound,
    first_noncancelling_paired_coefficient,
    manufactured_defect_majorant_control,
    negative_pair_and_omission_controls,
    remainder_recurrence_identities,
    thermal_tail_direct_window_check,
)


def test_m4_recurrence_uses_actual_transport_and_telescopes():
    result = remainder_recurrence_identities()
    assert result['expansion_order'] == 4
    assert result['defect'] == (
        'formal e^{-4} Pi_{-s} L0 A4 identity conditional on constructed '
        'A0..A4; production A3/A4 are not constructed here')
    assert not result['A3_locally_algebraic_in_metric_jets']
    assert not result['major_A2_replaced_by_local_metric_jet']
    assert result['dummy_symbol_L0_A_j']
    assert result['formal_symbolic_telescope']
    assert not result['production_A3_A4_constructed']
    residuals = result['residuals']
    zero_leaves = []
    for key, value in residuals.items():
        if key.endswith(('depends_on_major_A2', 'depends_on_major_Aj',
                         'depends_on_q', 'independent_of_q',
                         'defect_is_Pi_minus_L0_AM_over_eM')):
            continue
        zero_leaves.extend(value if isinstance(value, list) else [value])
    assert set(zero_leaves) == {'0'}


def test_first_noncancelling_pair_is_e_minus3_for_both_constraints():
    result = first_noncancelling_paired_coefficient()
    assert result['first_noncancelling_paired_inverse_energy_order'] == 3
    assert result['leading_paired_e_minus2_cancels_for_equal_mu']
    assert result['e_minus3_even_in_ell']
    assert set(result['involution_residuals'].values()) == {'0'}
    assert any('/N_parity_p3' in key for key in result['involution_residuals'])
    assert any('/beta_parity_p3' in key for key in result['involution_residuals'])
    assert result['numerical_C4'] is None


def test_missing_transport_keeps_c_m_none_and_supplied_inputs_close_primitive():
    missing = conditional_m4_remainder_majorant(None, None, 0.1, 0.2)
    assert missing['numerical_C_M'] is None
    assert len(missing['missing_inputs']) == 2
    complete = conditional_m4_remainder_majorant(
        (1., 2., 3.), (.1, .2, .3), .01, .02)
    assert complete['numerical_C_M'] is not None
    assert complete['numerical_C_M'] >= max(complete['C_M_h2'])
    assert complete['missing_inputs'] == ()


def test_thermal_tail_is_finite_but_not_rounded_to_zero():
    result = fermi_thermal_tail_bound(0., 5**.5, 6., kappa=.2, omega=4.)
    for row in result['cutoffs']:
        assert row['N'] is None and row['beta'] is None
        assert row['N_decimal_upper'] != '0.0'
        assert row['beta_decimal_upper'] != '0.0'
        assert row['log_N'] < -100
        assert row['below_provisional_allocation']
        assert row['underflow_is_not_exact_zero']


def test_manufactured_and_negative_controls_detect_omissions():
    manufactured = manufactured_defect_majorant_control()
    assert manufactured['periodic_manufactured_enclosed']
    assert manufactured['small_matrix_enclosed']
    assert manufactured['omitting_residual_fails_to_enclose']
    assert manufactured['tautological_zero_commutator_residual_fixture']
    assert not manufactured['production_L0_A4_control']
    assert not manufactured['certificate_use']
    negative = negative_pair_and_omission_controls()
    assert negative['equal_mu_e_minus2_max_abs'] == 0
    assert negative['broken_mu_e_minus2_max_abs'] > 0
    assert not negative['broken_mu_cancels']
    assert negative['paired_e_minus3_max_abs'] > 0
    assert negative['omitting_e_minus3_leaves_false_e_minus4_lead']
    assert negative['manufactured_even_parity_fixture']
    assert not negative['production_C4_coefficient_control']
    assert not negative['certificate_use']


def test_current_history_pilot_names_exact_missing_inputs_without_runs():
    result = current_history_remainder_pilot()
    assert result['profile_identity'] == CURRENT_HISTORY_IDENTITY
    assert result['first_noncancelling_paired_order'] == 3
    assert result['numerical_C4'] is None
    assert result['numerical_C_M'] is None
    assert len(result['majorant']['missing_inputs']) == 2
    assert result['field_or_source_evolutions'] == 0
    assert all(row['thermal_below_allocation'] for row in result['provisional_comparison'])
    assert all(row['vacuum_below_allocation'] is None for row in result['provisional_comparison'])


def test_successor_quantities_are_bound_and_remain_open():
    report = changed_history_uv_successor_report()
    assert report['schema'] == 'NSC-KS-CURRENT-UV-REMAINDER-v2'
    assert report['first_noncancelling_paired_order'] == 3
    assert report['leading_paired_e_minus2_cancels_for_equal_mu']
    assert report['leading_e_minus2_cancellation_is_not_a_tail_bound']
    assert report['e_minus3_first_allowed_order_is_not_a_tail_bound']
    assert report['numerical_C4'] is None
    assert report['numerical_C_M'] is None
    assert report['vacuum_integrated_tail_N_beta'] is None
    assert report['physical_local_gate'] == 'OPEN'
    assert report['field_or_source_evolutions'] == 0
    assert tuple(report['quantities']) == QUANTITY_NAMES
    quantities = changed_history_uv_quantities()
    assert quantities['A2']['imag_major_transport_rhs_this_massless_pair'] == 0.0
    assert quantities['A2']['imag_major_A2_numerical_value'] is None
    assert not quantities['A2']['upstream_imag_major_A2_owned']
    assert quantities['A2']['real_major_A2_numerical_value'] is None
    assert quantities['A3']['numerical_value'] is None
    assert not quantities['A3']['A3_locally_algebraic_in_metric_jets']
    assert quantities['n4']['retains_minor_A3']
    assert quantities['C4']['numerical_value'] is None
    assert quantities['C_M']['numerical_value'] is None
    assert quantities['C_M']['Bz_integral'] > 0
    assert quantities['current_history_L0_A4_H2_integrals']['numerical_value'] is None
    assert quantities['upstream_higher_order_H2_remainder']['numerical_value'] is None
    assert not quantities['upstream_higher_order_H2_remainder']['frozen_background_tail_copied']
    assert quantities['transported_majors']['real_major_A2'] is None
    assert quantities['transported_majors']['major_A3'] is None
    assert quantities['transported_majors']['imag_major_A2_this_massless_pair'] is None
    assert quantities['transported_majors'][
        'imag_major_A2_characteristic_integral_this_massless_pair'] is None
    assert quantities['transported_majors'][
        'imag_major_A2_transport_rhs_this_massless_pair'] == 0.0
    assert not quantities['transported_majors']['characteristic_integrals_owned']
    assert not quantities['C4']['history_minus_reference_coefficient_formed']
    assert quantities['vacuum_N_beta_tail']['N'] is None


def test_successor_bindings_and_direct_windows_authenticate():
    bindings = changed_history_uv_bindings()
    authenticate_changed_history_uv_bindings(bindings)
    assert bindings['profile']['identity'] == CURRENT_HISTORY_IDENTITY
    assert bindings['subtraction']['vertices']['numerical_residuals']['N_density'] == 0.0
    windows = thermal_tail_direct_window_check(0.0, 5 ** 0.5, 6.0, kappa=0.2, omega=4.0)
    assert windows['vacuum_integrated_tail_N_beta'] is None
    for row in windows['cutoffs']:
        assert row['N'] is None and row['vacuum_N'] is None
        assert row['window_plus_rest_encloses_infinite_N']
        assert row['rest_is_not_exact_zero']
    thermal = fermi_thermal_tail_bound(0.0, 5 ** 0.5, 6.0, kappa=0.2, omega=4.0)
    assert thermal['cutoffs'][0]['N'] is None
    assert thermal['cutoffs'][0]['N_decimal_upper'] != '0.0'
