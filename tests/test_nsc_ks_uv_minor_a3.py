"""Minor A3 jet bound: L0 parent formulae, second jets, no invented C4."""
from copy import deepcopy
from pathlib import Path
import json
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_ks_uv_minor_a3 as D
from recursive_horizons.nsc_ks_current_uv_transport import (
    CURRENT_HISTORY_IDENTITY,
    characteristic_transport_identities,
)
from recursive_horizons.nsc_ks_finite_history_uv_remainder import (
    FIRST_NONCANCELLING_PAIRED_ORDER,
)
from recursive_horizons.nsc_ks_uv_minor_a3 import (
    MISSING_UPSTREAM_MAJOR_A2,
    SCHEMA,
    minor_a3_identities,
    minor_a3_interval_majorants,
    second_geometry_phase_jet_majorants,
    validate_minor_a3_report,
)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


def zero_family():
    return LocalIncomingFamily(np.zeros((2, 8)))


def test_parent_formulae_use_actual_L0_and_retain_major_A2():
    identities = minor_a3_identities()
    assert not identities['dummy_symbol_L0_A_j']
    assert not identities['A3_locally_algebraic_in_metric_jets']
    assert not identities['production_A3_A4_constructed']
    assert not identities['massless_Im_Ts_major_A2_sets_upstream_datum_to_zero']
    assert identities['leading_paired_e_minus2_cancels_for_equal_mu']
    assert identities['leading_e_minus2_cancellation_is_not_this_value']
    residuals = identities['residuals']
    for source_sign in (1, -1):
        tag = 's%+d' % source_sign
        assert residuals['%s/minor_A2_parent' % tag] == '0'
        assert residuals[
            '%s/minor_A3_is_Tminus_minor_A2_plus_algebraic_major_A2' % tag] == '0'
        assert residuals['%s/minor_A3_depends_on_major_A2' % tag] == 'depends'
        assert residuals['%s/minor_A3_jet_independent_of_major_A2' % tag] == (
            'independent')
        assert residuals['%s/massless_Re_minor_A3_is_s_a_ell_h_over_2r' % tag] == '0'
        assert residuals[
            '%s/massless_Im_minor_A3_retains_R_and_Tminus_minor_A2' % tag] == '0'
        assert residuals['%s/massless_jet_independent_of_h' % tag] == 'independent'
        assert residuals['%s/T_s_q_z_on_shell_is_ell2_r_z_over_r3' % tag] == '0'
    parent = characteristic_transport_identities()
    assert not parent['dummy_symbol_L0_A_j']
    assert not parent['production_A3_A4_constructed']


def test_zero_history_difference_is_exact_zero_without_zeroing_upstream():
    family = zero_family()
    jets = second_geometry_phase_jet_majorants(family, 0, 5**0.5)
    assert jets['deformation_vanishes']
    assert jets['r_z_abs_upper'] == 0
    assert jets['q_z_abs_upper'] == 0
    assert jets['q_reference_z_abs_upper'] == 0
    assert jets['delta_radius_abs_upper'] == 0
    majorants = minor_a3_interval_majorants(family, 0, 5**0.5)
    assert majorants['history_minus_reference_minor_A2_abs_upper'] == 0
    assert majorants['history_minus_reference_Tminus_minor_A2_abs_upper'] == 0
    assert majorants['history_minus_reference_minor_A3_jet_abs_upper'] == 0
    assert majorants['generated_algebraic_major_A2_coupling_abs_upper'] == 0
    assert majorants['upstream_major_A2_difference_kernel_abs_upper'] == 0
    assert majorants['complete_history_minus_reference_minor_A3_abs_upper'] == 0
    assert majorants['complete_full_envelope_minor_A3_abs_upper'] is None
    assert not majorants['upstream_major_A2_owned']
    assert majorants['full_envelope_real_major_A2'] is None
    assert majorants['minor_A2_ref_abs_upper'] > 0
    assert majorants['minor_A3_jet_ref_abs_upper'] > 0


def test_live_history_encloses_jets_and_keeps_complete_minor_a3_open():
    value = D.calculate()
    assert value['schema'] == SCHEMA
    assert value['profile_identity'] == CURRENT_HISTORY_IDENTITY
    assert value['pair']['mass'] == 0.0
    assert value['pair']['absolute_angular'] > 0
    assert value['first_noncancelling_paired_order'] == FIRST_NONCANCELLING_PAIRED_ORDER
    assert value['leading_paired_e_minus2_cancels_for_equal_mu']
    assert value['leading_e_minus2_cancellation_is_not_this_value']
    assert value['e_minus3_first_allowed_order_is_not_a_tail_bound']
    assert value['numerical_C4'] is None
    assert value['numerical_C_M'] is None
    assert value['vacuum_integrated_tail_N_beta'] is None
    assert value['physical_local_gate'] == 'OPEN'
    assert value['missing_primitive'] == MISSING_UPSTREAM_MAJOR_A2
    assert not value['upstream_major_A2_owned']
    assert not value['production_A3_A4_constructed']
    assert not value['v1_pilot_rewritten']
    assert not value['v2_record_rewritten']
    assert not value['v3_record_rewritten']
    majorants = value['majorants']
    assert majorants['complete_history_minus_reference_minor_A3_abs_upper'] is None
    assert majorants['complete_full_envelope_minor_A3_abs_upper'] is None
    assert majorants['full_envelope_real_major_A2'] is None
    for name in (
            'q_z_abs_upper', 'q_zz_abs_upper', 'r_zz_abs_upper',
            'r_rho_rho_abs_upper', 'minor_A2_g_abs_upper',
            'history_minus_reference_minor_A2_abs_upper',
            'minor_A3_jet_g_abs_upper',
            'history_minus_reference_minor_A3_jet_abs_upper',
            'generated_algebraic_major_A2_coupling_abs_upper',
            'upstream_major_A2_difference_kernel_abs_upper'):
        assert majorants[name]['binary64_upper'] > 0
        assert majorants[name]['sign'] == 1
    assert majorants['q_abs_upper']['binary64_upper'] > 0
    validate_minor_a3_report(value, ROOT)


def test_invented_c4_cm_and_upstream_zero_are_rejected():
    value = D.calculate()
    forged = deepcopy(value)
    forged['numerical_C4'] = 0.0
    with pytest.raises(ValueError, match='presented as zero'):
        validate_minor_a3_report(forged, ROOT)
    forged_cm = deepcopy(value)
    forged_cm['numerical_C_M'] = 1.25
    with pytest.raises(ValueError, match='may not be invented'):
        validate_minor_a3_report(forged_cm, ROOT)
    forged_up = deepcopy(value)
    forged_up['upstream_major_A2_owned'] = True
    with pytest.raises(ValueError, match='upstream major A2'):
        validate_minor_a3_report(forged_up, ROOT)
    forged_complete = deepcopy(value)
    forged_complete['majorants'][
        'complete_history_minus_reference_minor_A3_abs_upper'] = 0.0
    with pytest.raises(ValueError, match='presented as zero'):
        validate_minor_a3_report(forged_complete, ROOT)
    forged_gate = deepcopy(value)
    forged_gate['physical_local_gate'] = 'PASS'
    with pytest.raises(ValueError, match='OPEN'):
        validate_minor_a3_report(forged_gate, ROOT)
    forged_jet = deepcopy(value)
    forged_jet['identities']['A3_locally_algebraic_in_metric_jets'] = True
    with pytest.raises(ValueError, match='local metric jet'):
        validate_minor_a3_report(forged_jet, ROOT)


def test_driver_replay_matches_owned_record():
    value = D.calculate()
    D.check(value)
    assert json.loads((ROOT / D.OUTPUT).read_text()) == value


def test_massive_channel_cannot_reuse_massless_difference_coupling():
    from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
    from recursive_horizons.nsc_ks_uv_minor_a3 import minor_a3_interval_majorants
    import numpy as np
    import pytest
    family = zero_family()
    with pytest.raises(ValueError, match="massless"):
        minor_a3_interval_majorants(family, 1, 2)
