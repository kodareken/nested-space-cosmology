"""Generated massless e^{-3} N,beta coefficient: actual vertices, no invented C4."""
from copy import deepcopy
from pathlib import Path
import json
import sys

import numpy as np
import pytest
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_ks_cubic_uv_current as D
from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter
from recursive_horizons.nsc_ks_current_uv_transport import (
    CURRENT_HISTORY_IDENTITY, HISTORY_RECORD, ARCHIVE_RHO_UP,
)
from recursive_horizons.nsc_ks_cubic_uv_current import (
    SCHEMA,
    confirm_action_signs_against_source_column_matter,
    evaluate_massless_cubic_on_history,
    massless_cubic_identities,
    unlinked_R_negative_control,
    validate_cubic_uv_report,
    geometry_jets, transport_massless_cubic_fields, _normal_window_jets,
    reduced_massless_transport_identities,
)
from recursive_horizons.nsc_dirac_source_phase import formal_source_phase_coefficient
from recursive_horizons.nsc_ks_finite_history_uv_remainder import (
    FIRST_NONCANCELLING_PAIRED_ORDER,
)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


def test_massless_reductions_and_actual_pauli_contraction():
    identities = massless_cubic_identities()
    assert not identities['dummy_symbol_L0_A_j']
    assert not identities['A3_locally_algebraic_in_metric_jets']
    assert not identities['R_zeroed_while_freezing_h3_n4']
    assert not identities['invariance_uses_n2_up_equals_zero']
    assert identities['all_linked_C_R_h3_n4_retained']
    assert identities['first_noncancelling_paired_order'] == FIRST_NONCANCELLING_PAIRED_ORDER
    assert identities['action_beta3'] == '+mu/(2 pi) J_I_3'
    assert identities['action_N3'].startswith('-mu/(2 pi)')
    residuals = identities['residuals']
    for source_sign in (1, -1):
        tag = 's%+d' % source_sign
        assert residuals['%s/l_is_minus_i_s_b' % tag] == '0'
        assert residuals['%s/B_reduction' % tag] == '0'
        assert residuals['%s/B_from_L0_A1' % tag] == '0'
        assert residuals['%s/n2_is_s_qz_plus_k' % tag] == '0'
        assert residuals['%s/C_from_recurrence' % tag] == '0'
        assert residuals['%s/C_massless_parent' % tag] == '0'
        assert residuals['%s/Ts_h3' % tag] == '0'
        assert residuals['%s/Ts_R3' % tag] == '0'
        assert residuals['%s/Re_conj_l_C_is_s_b_D' % tag] == '0'
        assert residuals['%s/continuity_e_minus4_is_Ts_n4_source' % tag] == '0'
        assert residuals['%s/J_I_3' % tag] == '0'
        assert residuals['%s/J_I_3_coefficient_of_h' % tag] == '0'
        assert residuals['%s/linked_scalar_R_cancels_in_J_I_3' % tag] == '0'
        assert residuals['%s/unlinked_R_spurious_c_qz' % tag] == '0'
        assert residuals['%s/J_I_3_retains_h3_n4_R_B' % tag] == 'depends'


def test_action_signs_match_source_column_matter_not_old_v3():
    signs = confirm_action_signs_against_source_column_matter()
    assert signs['beta'] == '+mu/(2 pi) I current'
    assert signs['N'] == '-mu/(2 pi) (V density + S3 current / a)'
    assert not signs['old_v3_beta_shorthand_is_authority']
    assert not signs['second_two_pi_inserted']
    for source_sign in (1, -1):
        columns = np.zeros((1, 2, 1), complex)
        columns[0, 0 if source_sign == 1 else 1, 0] = 1
        axial = -2j * source_sign * columns
        empty = np.zeros((0, *columns.shape), complex)
        result = source_column_matter(
            columns, axial, np.eye(1), empty, empty,
            mass=0, angular=0, axial_scale=1, radius=1, multiplicity=1)
        j_i = -2.0 * source_sign
        assert result['action_gradient'][0, 1] == j_i
        assert result['action_gradient'][0, 1] != -j_i or source_sign == 0
        assert signs['source_column_matter_residuals'][
            's%+d/beta_equals_plus_J_I' % source_sign]


def test_reduced_transport_and_surface_contraction_follow_original_vertices():
    proof=reduced_massless_transport_identities()
    assert len(proof['residuals']) == 6
    assert set(proof['residuals'].values()) == {'0'}
    assert proof['required_transport_quantities'] == ['q_z','q_zz','J3']
    assert proof['common_k_cancels_on_Sigma']
    assert proof['directed_quadrature_bound'] is None


def test_unlinked_R_is_the_spurious_control():
    control = unlinked_R_negative_control()
    assert control['unlinked_R_shift_of_J_I_3'] == 'c q_z'
    assert control['linked_scalar_shift_of_J_I_3'] == 0
    assert control['R_may_not_be_zeroed_independently']
    c, qz, rr, rz, q, h3z, n4, b, bz, B, Bz = sp.symbols(
        'c qz R Rz q h3z n4 b bz B Bz', real=True)
    for source_sign in (1, -1):
        current = h3z + rr * qz - q * rz + source_sign * (b * Bz - B * bz) - source_sign * n4
        assert sp.expand(current.subs(rr, rr + c) - current) == c * qz
        linked = current.xreplace({
            rr: rr + c, h3z: h3z + c * qz, n4: n4 + 2 * c * source_sign * qz})
        assert sp.expand(linked - current) == 0


def test_zero_history_difference_vanishes_without_zeroing_R():
    family = LocalIncomingFamily(np.zeros((2, 8)))
    value = evaluate_massless_cubic_on_history(
        family, angular=5**0.5, multiplicity=6.0, rho_steps=8, z_points=11)
    assert value['I_max_abs_N'] == 0.0
    assert value['I_max_abs_beta'] == 0.0
    assert value['R_from_continuity_and_k']
    assert not value['R_zeroed_independently']
    assert value['certified_quadrature_error_bound'] is None
    assert not value['validated']


def test_massive_channel_is_rejected():
    from recursive_horizons.nsc_ks_cubic_uv_current import _require_massless_angular
    with pytest.raises(ValueError, match='massless reduction'):
        _require_massless_angular(1.0, 5**0.5)
    with pytest.raises(ValueError, match='nonzero original angular'):
        _require_massless_angular(0.0, 0.0)


def test_inner_ramp_derivatives_survive_rounded_unit_window():
    inner, outer = .007, .03
    value, first, second = _normal_window_jets(inner+(outer-inner)*.02,inner,outer)
    assert value == 1.0  # rounded logistic value, not an exactly flat point
    assert first < 0 and second < 0


def test_mixed_geometry_jets_match_coordinate_differences():
    family = LocalIncomingFamily(np.asarray(json.loads((ROOT/HISTORY_RECORD).read_text())['history']['coefficients']))
    rho, step = 1.018, 2e-7
    z = family.center+np.array([-.041,-.017,.013,.044])
    j = geometry_jets(family,rho,z)
    zp,zm = geometry_jets(family,rho,z+step),geometry_jets(family,rho,z-step)
    rp,rm = geometry_jets(family,rho+step,z),geometry_jets(family,rho-step,z)
    for target,key,plus,minus in (
        ('r_zzz','r_zz',zp,zm),('r_rhozz','r_rhoz',zp,zm),
        ('r_rho_rhoz','r_rhoz',rp,rm),('r_rho_rho','r_rho',rp,rm),
    ):
        fd = (plus[key]-minus[key])/(2*step)
        assert np.max(abs(j[target]-fd)) < 1e-6*(1+np.max(abs(fd)))


def test_characteristic_phase_matches_independent_existing_gauss_owner():
    family = LocalIncomingFamily(np.asarray(json.loads((ROOT/HISTORY_RECORD).read_text())['history']['coefficients']))
    reference = LocalIncomingFamily(np.zeros((2,32)))
    z = np.linspace(*family.interval,5)
    phase = formal_source_phase_coefficient(family,z,angular=np.sqrt(5),
        rho_up=ARCHIVE_RHO_UP,gauss_nodes=48)
    for index,sign in enumerate((1,-1)):
        changed = transport_massless_cubic_fields(family,z,angular=np.sqrt(5),
            source_sign=sign,rho_steps=96)
        base = transport_massless_cubic_fields(reference,z,angular=np.sqrt(5),
            source_sign=sign,rho_steps=96)
        assert np.max(abs(changed['q']-base['q']-phase.f[index])) < 1e-11
        assert np.max(abs(changed['q_z']-base['q_z']-phase.f_z[index])) < 1e-10


def test_live_history_measures_unvalidated_cubic_coefficient():
    value = D.calculate()
    assert value['schema'] == SCHEMA
    assert value['profile_identity'] == CURRENT_HISTORY_IDENTITY
    assert value['pair']['mass'] == 0.0
    assert value['pair']['absolute_angular'] > 0
    assert value['first_noncancelling_paired_order'] == FIRST_NONCANCELLING_PAIRED_ORDER
    assert value['generated_massless_A3_n4_constructed']
    assert not value['production_A4_constructed']
    assert not value['validated_cubic_coefficient']
    assert value['numerical_C4'] is None
    assert value['numerical_C_M'] is None
    assert value['vacuum_integrated_tail_N_beta'] is None
    assert value['certified_quadrature_error_bound'] is None
    assert value['physical_local_gate'] == 'OPEN'
    assert not value['R_zeroed_while_freezing_h3_n4']
    assert not value['source_columns_changed']
    diagnostic = value['diagnostic_quadrature']
    assert not diagnostic['validated']
    assert not diagnostic['remainder_proof_supplied']
    assert diagnostic['certified_quadrature_error_bound'] is None
    assert diagnostic['I_max_abs_N'] > 0
    assert diagnostic['I_max_abs_beta'] > 0
    assert np.isfinite(diagnostic['I_integral_N'])
    assert np.isfinite(diagnostic['I_integral_beta'])
    assert np.isfinite(diagnostic['refinement_discrepancy_N'])
    assert np.isfinite(diagnostic['refinement_discrepancy_beta'])
    assert diagnostic['ell_odd_max_abs_N'] == 0.0
    assert diagnostic['ell_odd_max_abs_beta'] == 0.0
    assert len(diagnostic['z_I']) == len(diagnostic['action_N3'])
    assert value['zero_history_control']['I_max_abs_N'] == 0.0
    assert value['zero_history_control']['I_max_abs_beta'] == 0.0
    validate_cubic_uv_report(value, ROOT)


def test_invented_c4_cm_and_independent_R_zero_are_rejected():
    value = D.calculate()
    forged_reduction = deepcopy(value)
    forged_reduction['reduced_transport']['residuals']['s+1/Ts_J3'] = '1'
    with pytest.raises(ValueError, match='reduced transport'):
        validate_cubic_uv_report(forged_reduction, ROOT)
    forged = deepcopy(value)
    forged['numerical_C4'] = 0.0
    with pytest.raises(ValueError, match='presented as zero'):
        validate_cubic_uv_report(forged, ROOT)
    forged_cm = deepcopy(value)
    forged_cm['numerical_C_M'] = 1.25
    with pytest.raises(ValueError, match='may not be invented'):
        validate_cubic_uv_report(forged_cm, ROOT)
    forged_valid = deepcopy(value)
    forged_valid['validated_cubic_coefficient'] = True
    with pytest.raises(ValueError, match='remainder-validated'):
        validate_cubic_uv_report(forged_valid, ROOT)
    forged_r = deepcopy(value)
    forged_r['R_zeroed_while_freezing_h3_n4'] = True
    with pytest.raises(ValueError, match='independently'):
        validate_cubic_uv_report(forged_r, ROOT)
    forged_gate = deepcopy(value)
    forged_gate['physical_local_gate'] = 'PASS'
    with pytest.raises(ValueError, match='OPEN'):
        validate_cubic_uv_report(forged_gate, ROOT)
    forged_v3 = deepcopy(value)
    forged_v3['action_signs']['old_v3_beta_shorthand_is_authority'] = True
    with pytest.raises(ValueError, match='shorthand'):
        validate_cubic_uv_report(forged_v3, ROOT)
