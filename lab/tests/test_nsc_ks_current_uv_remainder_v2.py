"""Changed-history UV remainder successor: bindings, nulls and negative controls."""
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import json
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_ks_current_uv_remainder_v2 as D
from recursive_horizons.nsc_ks_finite_history_uv_remainder import (
    CURRENT_HISTORY_IDENTITY, NONE_REQUIRED_QUANTITIES, PILOT_RECORD,
    PILOT_SCHEMA, QUANTITY_NAMES, SUCCESSOR_SCHEMA,
    validate_changed_history_uv_successor,
)


def test_successor_binds_quantities_and_keeps_missing_values_none():
    value = D.calculate()
    assert value['schema'] == SUCCESSOR_SCHEMA
    assert value['status'].startswith('OPEN')
    assert value['first_noncancelling_paired_order'] == 3
    assert value['leading_paired_e_minus2_cancels_for_equal_mu']
    assert value['leading_e_minus2_cancellation_is_not_a_tail_bound']
    assert value['e_minus3_first_allowed_order_is_not_a_tail_bound']
    assert value['numerical_C4'] is None
    assert value['numerical_C_M'] is None
    assert value['vacuum_integrated_tail_N_beta'] is None
    assert value['physical_EXISTENCE_certificate'] is False
    assert value['physical_NONEXISTENCE_certificate'] is False
    assert value['v1_pilot_rewritten'] is False
    remainder = value['remainder']
    assert remainder['identities']['dummy_symbol_L0_A_j']
    assert remainder['identities']['formal_symbolic_telescope']
    assert not remainder['identities']['production_A3_A4_constructed']
    assert remainder['defect']['dummy_symbol_L0_A_j']
    assert remainder['defect']['formal_symbolic_telescope']
    assert not remainder['defect']['production_A3_A4_constructed']
    assert set(remainder['quantities']) == set(QUANTITY_NAMES)
    for name in NONE_REQUIRED_QUANTITIES:
        item = remainder['quantities'][name]
        assert item['numerical_value'] is None
        assert item['status'] == 'OPEN'
        assert item['bindings']['profile'] == CURRENT_HISTORY_IDENTITY
        assert item['bindings']['preparation_restart']
        assert item['bindings']['preparation_cauchy']
        assert item['bindings']['inventory_record']
        assert item['bindings']['inventory_cutoff_bridge']
    vacuum = remainder['quantities']['vacuum_N_beta_tail']
    assert vacuum['N'] is None and vacuum['beta'] is None
    assert remainder['field_or_source_evolutions'] == 0
    D.check(value)


def test_bindings_authenticate_state_law_profile_preparation_inventory_subtraction():
    value = D.calculate()
    bindings = value['remainder']['bindings']
    assert bindings['state_law']['owner'].endswith('nsc_evolved_incoming_state.py')
    assert 'C_Sigma[g]=F[g] C_src F[g]^dagger' in bindings['state_law']['statement']
    assert bindings['profile']['identity'] == CURRENT_HISTORY_IDENTITY
    assert bindings['preparation']['owner'].endswith('nsc_paired_horizon_preparation.py')
    assert bindings['preparation']['physical_source_renormalized'] is False
    assert bindings['inventory']['payload_sha256'] == sha256(
        (ROOT / bindings['inventory']['payload_path']).read_bytes()).hexdigest()
    assert bindings['subtraction']['owner'].endswith('nsc_common_subtracted_ks_source.py')
    assert bindings['subtraction']['background_tail_certificate_transferred'] is False
    residuals = bindings['subtraction']['vertices']['numerical_residuals']
    assert residuals['N_density'] == 0.0
    assert residuals['N_momentum'] == 0.0
    assert residuals['beta_momentum'] == 0.0
    D.check(value)

    forged = deepcopy(value['remainder'])
    forged['bindings']['preparation']['surface_gravity'] += 1e-6
    with pytest.raises(ValueError, match='binding ledger differs'):
        validate_changed_history_uv_successor(forged, ROOT)
    forged_cauchy = deepcopy(value['remainder'])
    forged_cauchy['bindings']['preparation']['cauchy_record_sha256'] = '0' * 64
    with pytest.raises(ValueError, match='binding ledger differs'):
        validate_changed_history_uv_successor(forged_cauchy, ROOT)
    forged_cutoff = deepcopy(value['remainder'])
    forged_cutoff['bindings']['inventory']['cutoff_bridge_sha256'] = 'f' * 64
    with pytest.raises(ValueError, match='binding ledger differs'):
        validate_changed_history_uv_successor(forged_cutoff, ROOT)


def test_manufactured_and_negative_controls_are_retained():
    value = D.calculate()
    manufactured = value['remainder']['manufactured']
    assert manufactured['periodic_manufactured_enclosed']
    assert manufactured['small_matrix_enclosed']
    assert manufactured['omitting_residual_fails_to_enclose']
    assert manufactured['current_history_C_M'] is None
    assert manufactured['tautological_zero_commutator_residual_fixture']
    assert not manufactured['production_L0_A4_control']
    assert not manufactured['certificate_use']
    negative = value['remainder']['negative_controls']
    assert negative['equal_mu_e_minus2_max_abs'] == 0
    assert negative['broken_mu_e_minus2_max_abs'] > 0
    assert not negative['broken_mu_cancels']
    assert negative['paired_e_minus3_max_abs'] > 0
    assert negative['omitting_e_minus3_leaves_false_e_minus4_lead']
    assert negative['leading_e_minus2_cancellation_is_not_a_tail_bound']
    assert negative['manufactured_even_parity_fixture']
    assert not negative['production_C4_coefficient_control']
    assert not negative['certificate_use']


def test_direct_windows_enclose_thermal_tails_without_inventing_vacuum():
    value = D.calculate()
    windows = value['remainder']['direct_windows']
    assert windows['vacuum_integrated_tail_N_beta'] is None
    assert windows['background_tail_certificate_transferred'] is False
    assert windows['vacuum_envelope_used'] is False
    for row in windows['cutoffs']:
        assert row['N'] is None and row['beta'] is None
        assert row['vacuum_N'] is None and row['vacuum_beta'] is None
        assert row['infinite_N_decimal_upper'] != '0.0'
        assert row['window_N_decimal_upper'] != '0.0'
        assert row['rest_N_decimal_upper'] != '0.0'
        assert row['window_plus_rest_encloses_infinite_N']
        assert row['window_plus_rest_encloses_infinite_beta']
        assert row['window'][1] == 2.0 * row['cutoff']


def test_validator_rejects_invented_constants_and_false_certificates():
    value = D.calculate()
    forged_c4 = deepcopy(value['remainder'])
    forged_c4['quantities']['C4']['numerical_value'] = 0.0
    with pytest.raises(ValueError, match='presented as zero'):
        validate_changed_history_uv_successor(forged_c4, ROOT)
    forged_cm = deepcopy(value['remainder'])
    forged_cm['quantities']['C_M']['numerical_value'] = 1.25
    with pytest.raises(ValueError, match='may not be invented'):
        validate_changed_history_uv_successor(forged_cm, ROOT)
    for name in ('A2', 'A3', 'n4', 'transported_majors'):
        forged_quantity = deepcopy(value['remainder'])
        forged_quantity['quantities'][name]['numerical_value'] = 0.0
        with pytest.raises(ValueError):
            validate_changed_history_uv_successor(forged_quantity, ROOT)
    forged_imag = deepcopy(value['remainder'])
    forged_imag['quantities']['A2']['imag_major_A2_numerical_value'] = 0.0
    with pytest.raises(ValueError, match='presented as zero'):
        validate_changed_history_uv_successor(forged_imag, ROOT)
    forged_transport = deepcopy(value['remainder'])
    forged_transport['quantities']['transported_majors'][
        'imag_major_A2_this_massless_pair'] = 0.0
    with pytest.raises(ValueError, match='presented as zero'):
        validate_changed_history_uv_successor(forged_transport, ROOT)
    forged_thermal = deepcopy(value['remainder'])
    forged_thermal['thermal']['cutoffs'][0]['N'] = 0.0
    with pytest.raises(ValueError, match='presented as zero'):
        validate_changed_history_uv_successor(forged_thermal, ROOT)
    forged_formal = deepcopy(value['remainder'])
    forged_formal['identities']['production_A3_A4_constructed'] = True
    with pytest.raises(ValueError, match='formal M=4 telescope'):
        validate_changed_history_uv_successor(forged_formal, ROOT)
    forged_pass = deepcopy(value['remainder'])
    forged_pass['physical_EXISTENCE_certificate'] = True
    with pytest.raises(ValueError, match='PASS or NON_EXISTENCE'):
        validate_changed_history_uv_successor(forged_pass, ROOT)
    forged_tail = deepcopy(value['remainder'])
    forged_tail['certificate_from_leading_e_minus2_cancellation'] = True
    with pytest.raises(ValueError, match='not a tail bound'):
        validate_changed_history_uv_successor(forged_tail, ROOT)
    forged_background = deepcopy(value['remainder'])
    forged_background['quantities']['vacuum_N_beta_tail'][
        'background_tail_certificate_transferred'] = True
    with pytest.raises(ValueError, match='must not transfer'):
        validate_changed_history_uv_successor(forged_background, ROOT)


def test_v1_pilot_bytes_and_schema_are_unchanged():
    value = D.calculate()
    pilot = (ROOT / PILOT_RECORD).read_bytes()
    recorded = json.loads(pilot)
    assert recorded['schema'] == PILOT_SCHEMA
    assert sha256(pilot).hexdigest() == value['v1_pilot_sha256']
    assert recorded['numerical_C4'] is None
    assert recorded['numerical_C_M'] is None
    assert recorded['vacuum_integrated_tail_N_beta'] is None
    assert value['v1_pilot_rewritten'] is False
