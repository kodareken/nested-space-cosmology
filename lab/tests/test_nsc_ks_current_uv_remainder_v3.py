"""Changed-history UV remainder v3: transports, C4 scope and negative controls."""
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import json
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_ks_current_uv_remainder_v3 as D
from recursive_horizons.nsc_ks_current_uv_transport import (
    C4_SCOPE, CURRENT_HISTORY_IDENTITY, NEXT_PRIMITIVE,
    characteristic_transport_identities,
    history_minus_reference_imag_major_A2,
    massless_j3_drops_upstream_imag_major_A2,
    validate_characteristic_transports,
)
from recursive_horizons.nsc_ks_current_uv_report import (
    PILOT_RECORD, PILOT_SCHEMA, SUCCESSOR_SCHEMA, V2_REMAINDER_RECORD,
    V3_SUCCESSOR_SCHEMA, QUANTITY_NAMES, validate_changed_history_uv_v3,
)


def test_transport_identities_use_actual_L0_not_dummy_symbols():
    identities = characteristic_transport_identities()
    assert identities['production_A2_transport_constructed']
    assert not identities['dummy_symbol_L0_A_j']
    assert not identities['production_A3_A4_constructed']
    assert not identities['massless_Im_Ts_major_A2_sets_upstream_datum_to_zero']
    assert not identities['massless_Im_Ts_major_A2_sets_full_envelope_coefficient_to_zero']
    residuals = identities['residuals']
    for source_sign in (1, -1):
        tag = 's%+d' % source_sign
        assert residuals['%s/Im_Ts_major_A2_parent' % tag] == '0'
        assert residuals['%s/Re_Ts_major_A2_parent' % tag] == '0'
        assert residuals['%s/Im_Ts_major_A2_vanishes_for_m_0' % tag] == '0'
        assert residuals['%s/Re_Ts_major_A2_depends_on_q' % tag] == 'depends'
        assert residuals['%s/minor_A3_depends_on_major_A2' % tag] == 'depends'
        assert residuals['%s/Re_Ts_major_A3_massless_is_ell2_h_over_2r2' % tag] == '0'


def test_massless_difference_is_zero_without_zeroing_upstream_or_full_envelope():
    closed = history_minus_reference_imag_major_A2(0.0)
    assert closed['numerical_value'] == 0.0
    assert closed['proved_exact']
    assert closed['scope'] == 'history-minus-reference'
    assert closed['full_envelope_value'] is None
    assert not closed['upstream_imag_major_A2_owned']
    assert closed['transport_rhs_this_massless_pair'] == 0.0
    massive = history_minus_reference_imag_major_A2(1.0)
    assert massive['numerical_value'] is None
    assert massive['full_envelope_value'] is None
    j3 = massless_j3_drops_upstream_imag_major_A2()
    assert j3['massless_J3_independent_of_upstream_imag_major_A2']
    assert not j3['massless_Im_Ts_major_A2_sets_upstream_datum_to_zero']
    assert not j3['full_envelope_imag_major_A2_set_to_zero']


def test_v3_successor_proves_transports_and_keeps_c4_cm_none():
    value, payload = D.calculate()
    assert value['schema'] == V3_SUCCESSOR_SCHEMA
    assert value['status'].startswith('OPEN')
    assert value['first_noncancelling_paired_order'] == 3
    assert value['leading_paired_e_minus2_cancels_for_equal_mu']
    assert value['numerical_C4'] is None
    assert value['numerical_C_M'] is None
    assert value['vacuum_integrated_tail_N_beta'] is None
    assert value['physical_EXISTENCE_certificate'] is False
    assert value['physical_NONEXISTENCE_certificate'] is False
    assert value['v1_pilot_rewritten'] is False
    assert value['v2_record_rewritten'] is False
    assert value['next_primitive'] == NEXT_PRIMITIVE
    proved = value['proved']
    assert proved['history_minus_reference_imag_major_A2_massless'] == 0.0
    assert proved['q_abs_upper'] > 0
    assert proved['delta_q_abs_upper'] > 0
    assert proved['history_minus_reference_real_major_A2_abs_upper'] > 0
    assert proved['history_minus_reference_C4_definition_formed']
    assert proved['numerical_C4'] is None
    remainder = value['remainder']
    assert remainder['quantities']['C4']['history_minus_reference_coefficient_formed']
    assert remainder['quantities']['C4']['definition_scope'] == C4_SCOPE
    assert not remainder['quantities']['C4']['full_envelope_coefficient_formed']
    assert remainder['quantities']['A2']['history_minus_reference_imag_major_A2'] == 0.0
    assert remainder['quantities']['A2']['imag_major_A2_numerical_value'] is None
    assert not remainder['quantities']['A2']['upstream_imag_major_A2_owned']
    assert remainder['identities']['production_A2_transport_constructed']
    assert not remainder['identities']['dummy_symbol_L0_A_j']
    assert not remainder['identities']['production_A3_A4_constructed']
    assert remainder['transport']['profile_identity'] == CURRENT_HISTORY_IDENTITY
    assert sha256(payload).hexdigest() == value['payload']['sha256']
    D.check(value, payload)


def test_bindings_and_profile_identity_are_authenticated():
    value, payload = D.calculate()
    bindings = value['remainder']['bindings']
    assert bindings['state_law']['owner'].endswith('nsc_evolved_incoming_state.py')
    assert 'C_Sigma[g]=F[g] C_src F[g]^dagger' in bindings['state_law']['statement']
    assert bindings['profile']['identity'] == CURRENT_HISTORY_IDENTITY
    assert bindings['preparation']['physical_source_renormalized'] is False
    assert bindings['subtraction']['background_tail_certificate_transferred'] is False
    D.check(value, payload)


def test_negative_invented_a2_a3_n4_c4_cm_are_rejected():
    value, _ = D.calculate()
    for name in ('A2', 'A3', 'n4', 'C4', 'C_M'):
        forged = deepcopy(value['remainder'])
        forged['quantities'][name]['numerical_value'] = 0.0
        with pytest.raises(ValueError):
            validate_changed_history_uv_v3(forged, ROOT)
    forged_full = deepcopy(value['remainder'])
    forged_full['quantities']['A2']['imag_major_A2_numerical_value'] = 0.0
    with pytest.raises(ValueError, match='presented as zero'):
        validate_changed_history_uv_v3(forged_full, ROOT)
    forged_transport = deepcopy(value['remainder'])
    forged_transport['quantities']['transported_majors'][
        'imag_major_A2_this_massless_pair'] = 0.0
    with pytest.raises(ValueError, match='presented as zero'):
        validate_changed_history_uv_v3(forged_transport, ROOT)
    forged_cm = deepcopy(value['remainder'])
    forged_cm['quantities']['C_M']['numerical_value'] = 1.25
    with pytest.raises(ValueError, match='may not be invented'):
        validate_changed_history_uv_v3(forged_cm, ROOT)


def test_missing_upstream_datum_is_not_coerced_to_zero():
    value, _ = D.calculate()
    forged = deepcopy(value['remainder'])
    forged['quantities']['A2']['upstream_imag_major_A2_owned'] = True
    forged['transport']['identities'][
        'massless_Im_Ts_major_A2_sets_upstream_datum_to_zero'] = True
    with pytest.raises(ValueError, match='upstream datum'):
        validate_changed_history_uv_v3(forged, ROOT)
    forged_value = deepcopy(value['remainder'])
    forged_value['quantities']['A2']['imag_major_A2_numerical_value'] = 0.0
    forged_value['quantities']['A2']['upstream_imag_major_A2_owned'] = False
    with pytest.raises(ValueError):
        validate_changed_history_uv_v3(forged_value, ROOT)


def test_wrong_coefficient_scope_is_rejected():
    value, _ = D.calculate()
    forged = deepcopy(value['remainder'])
    forged['quantities']['C4']['definition_scope'] = (
        'full-envelope e^{-3} N,beta coefficient of g, including the reference '
        'affine column')
    forged['quantities']['C4']['full_envelope_coefficient_formed'] = True
    with pytest.raises(ValueError, match='full-envelope'):
        validate_changed_history_uv_v3(forged, ROOT)
    forged_a2 = deepcopy(value['remainder'])
    forged_a2['quantities']['A2']['coefficient_scope'] = 'full-envelope'
    with pytest.raises(ValueError, match='wrong coefficient scope'):
        validate_changed_history_uv_v3(forged_a2, ROOT)
    forged_transport = deepcopy(value['remainder']['transport'])
    forged_transport['imag_major_A2']['scope'] = 'full-envelope'
    with pytest.raises(ValueError, match='wrong coefficient scope'):
        validate_characteristic_transports(forged_transport, ROOT)


def test_transferred_background_tail_and_forged_bindings_are_rejected():
    value, _ = D.calculate()
    forged = deepcopy(value['remainder'])
    forged['quantities']['vacuum_N_beta_tail'][
        'background_tail_certificate_transferred'] = True
    with pytest.raises(ValueError, match='must not transfer'):
        validate_changed_history_uv_v3(forged, ROOT)
    forged_up = deepcopy(value['remainder'])
    forged_up['quantities']['upstream_higher_order_H2_remainder'][
        'frozen_background_tail_copied'] = True
    with pytest.raises(ValueError, match='without a proved bridge'):
        validate_changed_history_uv_v3(forged_up, ROOT)
    forged_bind = deepcopy(value['remainder'])
    forged_bind['bindings']['preparation']['surface_gravity'] += 1e-6
    with pytest.raises(ValueError, match='binding ledger differs'):
        validate_changed_history_uv_v3(forged_bind, ROOT)
    forged_cauchy = deepcopy(value['remainder'])
    forged_cauchy['bindings']['preparation']['cauchy_record_sha256'] = '0' * 64
    with pytest.raises(ValueError, match='binding ledger differs'):
        validate_changed_history_uv_v3(forged_cauchy, ROOT)


def test_manufactured_fixtures_cannot_supply_production_constants():
    value, _ = D.calculate()
    remainder = value['remainder']
    assert remainder['manufactured']['tautological_zero_commutator_residual_fixture']
    assert not remainder['manufactured']['production_L0_A4_control']
    assert remainder['manufactured']['current_history_C_M'] is None
    assert not remainder['negative_controls']['production_C4_coefficient_control']
    forged = deepcopy(remainder)
    forged['quantities']['C_M']['manufactured_fixture_used_as_production_constant'] = True
    with pytest.raises(ValueError, match='manufactured fixture'):
        validate_changed_history_uv_v3(forged, ROOT)
    forged_l0 = deepcopy(remainder)
    forged_l0['quantities']['current_history_L0_A4_H2_integrals'][
        'production_L0_A4_control'] = True
    with pytest.raises(ValueError, match='manufactured fixture'):
        validate_changed_history_uv_v3(forged_l0, ROOT)


def test_historical_v1_and_v2_bytes_are_unchanged():
    value, payload = D.calculate()
    pilot = (ROOT / PILOT_RECORD).read_bytes()
    v2 = (ROOT / V2_REMAINDER_RECORD).read_bytes()
    assert json.loads(pilot)['schema'] == PILOT_SCHEMA
    assert json.loads(v2)['schema'] == SUCCESSOR_SCHEMA
    assert sha256(pilot).hexdigest() == value['v1_pilot_sha256']
    assert sha256(v2).hexdigest() == value['v2_record_sha256']
    assert value['v1_pilot_rewritten'] is False
    assert value['v2_record_rewritten'] is False
    D.check(value, payload)


def test_quantity_names_and_gate_remain_open():
    value, _ = D.calculate()
    remainder = value['remainder']
    assert set(remainder['quantities']) == set(QUANTITY_NAMES)
    assert remainder['physical_local_gate'] == 'OPEN'
    assert 'Gate 2 remains OPEN' in remainder['budget_implication']
    assert remainder['next_primitive'] == NEXT_PRIMITIVE
    for name in QUANTITY_NAMES:
        assert remainder['quantities'][name]['numerical_value'] is None
        assert remainder['quantities'][name]['status'] == 'OPEN'
