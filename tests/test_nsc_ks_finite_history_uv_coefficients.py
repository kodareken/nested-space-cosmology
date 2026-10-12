"""Constructed vacuum-envelope coefficients, not a dummy L0 telescope."""
import json
import time
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_compatible_history_geometry import _plateau
from recursive_horizons.nsc_ks_finite_history_uv_coefficients import (
    ARCHIVE_RHO_UP,
    CURRENT_HISTORY_IDENTITY_PREFIX,
    PHYSICAL_LOCAL_GATE,
    SOURCE_SIGNS,
    V2_RECORD,
    _plateau_jet,
    chebyshev_profile_point_enclosures,
    coefficient_identities,
    confirm_raw_ks_vertices,
    e_minus2_current_bookkeeping,
    e_minus2_sigma_contraction,
    finite_history_uv_coefficient_report,
    integrate_imag_major_A2,
    integrate_imag_major_A2_z,
    named_gaps,
    occupied_basis,
    original_angular_pair_ledger,
    paired_leading_e_minus2_weight,
    radius_partial_jets,
    vacuum_upstream_conditions,
)
from recursive_horizons.nsc_ks_linear_uv_contractions import contraction_identity
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


ROOT = Path(__file__).resolve().parents[1]
CURRENT_HISTORY = ROOT / 'results/development/nsc-ks-gate-history-lm-broyden.json'
ANGULAR = float(np.sqrt(5.0))
MASS = 1.0


def n8_family():
    coeff = np.zeros((2, 8))
    coeff[0, :4] = [0.03, -0.02, 0.01, 0.005]
    coeff[1, :3] = [1.2, -0.4, 0.25]
    return LocalIncomingFamily(coeff)


def zero_family(n=8):
    return LocalIncomingFamily(np.zeros((2, n)))


def current_family():
    if not CURRENT_HISTORY.exists():
        pytest.skip('current 0b0e4ced history record is not present')
    payload = json.loads(CURRENT_HISTORY.read_text())
    return LocalIncomingFamily(np.array(payload['history']['coefficients']))


def test_parent_algebraic_leads_and_L0_telescope_are_exact():
    result = coefficient_identities()
    residuals = result['residuals']
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        assert residuals['%s/A0_is_e_s' % tag] == ['0', '0']
        assert residuals['%s/minor_A1_parent' % tag] == '0'
        assert residuals['%s/minor_A1_from_L0_A0' % tag] == '0'
        assert residuals['%s/T_s_major_A1_parent' % tag] == '0'
        assert residuals['%s/Im_T_s_major_A2_parent' % tag] == '0'
        assert residuals['%s/L0_telescope_M2' % tag] == ['0', '0']
        assert residuals['%s/Im_T_s_major_A2_odd_in_ell' % tag] == '0'
        assert residuals['%s/minor_A3_depends_on_major_A2' % tag] == 'depends'
        assert residuals['%s/Re_major_A1_normalization' % tag] == '0'
    assert residuals['s+1/minor_A1_matches_P1'] == '0'
    assert result['dummy_symbol_L0_A_j'] is False
    assert result['physical_source_columns_replaced'] is False


def test_homogeneous_and_free_controls_kill_the_constructed_imaginary_A2_rhs():
    residuals = coefficient_identities()['residuals']
    for source_sign in SOURCE_SIGNS:
        tag = 's%+d' % source_sign
        assert residuals['%s/Im_T_s_major_A2_vanishes_for_m_0' % tag] == '0'
        assert residuals['%s/Im_T_s_major_A2_vanishes_for_ell_0' % tag] == '0'


def test_e_minus2_bookkeeping_keeps_h_z_from_being_the_full_coefficient():
    current = e_minus2_current_bookkeeping()
    for source_sign in SOURCE_SIGNS:
        terms = current['terms']['s%+d' % source_sign]
        assert terms['h_z_is_the_full_e_minus2_coefficient'] is False
        assert terms['beta_e_minus2_contains_h_z'] is True
        assert terms['beta_e_minus2_contains_A1_interference'] is True
        assert terms['beta_e_minus2_contains_n3_hence_A3'] is True
        assert terms['A3_locally_algebraic_in_metric_jets'] is False
        assert terms['Im_A1dag_A1_z_residual'] == '0'
        assert 'S3/a' in terms['N_momentum_vertex']
        assert terms['beta_momentum_vertex'].startswith('-I')
    linear = contraction_identity()
    assert linear['residuals']['N_E_minus2'] == '0'
    assert linear['scope']['finite_history_bound'] is False
    assert 'not a finite-history' in current['linear_fixed_transfer_e_minus2_zero']


def test_gaps_and_upstream_conditions_stay_explicit():
    gaps = named_gaps()
    assert gaps['numerical_C_M'] is None
    assert gaps['thermal_coherent_source_remainder'] is None
    assert gaps['e_minus3_paired_C4'] is None
    assert gaps['complete_e_minus2_N_beta_coefficient'] is None
    assert gaps['physical_local_gate'] == 'OPEN'
    assert 'leading 1/Lambda' in gaps['finite_history_UV_cancellation']
    upstream = vacuum_upstream_conditions()
    assert upstream['A0'] == 'e_s, real positive major'
    assert upstream['major_A1'] == 0
    assert upstream['physical_source_renormalized'] is False
    report = finite_history_uv_coefficient_report()
    assert report['schema'] == 'NSC-KS-FINITE-HISTORY-UV-COEFFICIENTS-v3'
    assert report['certificate_from_asymptotic_orders_alone'] is False
    assert report['C_M_or_physical_gate_from_leading_cancellation'] is False
    assert report['existing_source_ad759_evaluation_changed'] is False
    assert report['physical_local_gate'] == PHYSICAL_LOCAL_GATE
    assert occupied_basis(1).tolist() == [1.0, 0.0]
    assert occupied_basis(-1).tolist() == [0.0, 1.0]


def test_n3_continuity_reduces_sigma_e_minus2_without_solving_A3():
    result = e_minus2_sigma_contraction()
    residuals = result['residuals']
    assert set(residuals.values()) == {'0'}
    assert result['n3'] == '2 s h_z after T_s transport and normalized z-homogeneous upstream'
    assert result['Sigma_I_momentum_e_minus2'] == '-h_z'
    assert result['action_N2'] == 'mu s delta(h_z) / (2 pi a)'
    assert result['action_beta2'] == '-mu delta(h_z) / (2 pi)'
    assert result['local_delta_r_rho_cancels_in_N'] is True
    assert result['A3_solved_componentwise'] is False
    assert result['n3_eliminated_by_continuity_and_upstream'] is True
    assert result['thermal_occupation_replaced_by_one'] is False
    assert result['C4_or_C_M_bounded'] is False
    assert result['leading_paired_one_over_Lambda_cancels_for_equal_mu'] is True
    vertices = confirm_raw_ks_vertices()
    assert vertices['numerical_residuals']['N_density'] == 0.0
    assert vertices['numerical_residuals']['N_momentum'] == 0.0
    assert vertices['numerical_residuals']['beta_density'] == 0.0
    assert vertices['numerical_residuals']['beta_momentum'] == 0.0


def test_original_pairs_have_equal_mu_and_broken_weights_do_not_cancel():
    pairs = original_angular_pair_ledger()
    assert len(pairs) == 30
    assert all(row['same_cutoff_grid'] for row in pairs)
    assert all(row['source_grid_and_weight_residual'] == 0.0 for row in pairs)
    h = np.array([1.25, -0.5, 0.25])
    mu = pairs[0]['multiplicity_per_signed_family']
    assert np.array_equal(paired_leading_e_minus2_weight(h, -h, mu, mu), np.zeros(3))
    broken = paired_leading_e_minus2_weight(h, -h, mu, mu * 0.5)
    assert not np.allclose(broken, 0.0)
    assert broken == pytest.approx(mu * 0.5 * h)


def test_plateau_jet_matches_owned_window_and_finite_difference():
    inner, outer = 0.007, 0.03
    step = 1e-8
    for normal in (0.0, -0.01, -0.02, -0.029, -0.04):
        value, deriv = _plateau_jet(normal, inner, outer)
        owned = float(_plateau(np.array([normal]), inner, outer)[0])
        assert value == pytest.approx(owned, abs=0.0, rel=0.0)
        finite = (
            _plateau_jet(normal + step, inner, outer)[0]
            - _plateau_jet(normal - step, inner, outer)[0]) / (2.0 * step)
        assert deriv == pytest.approx(finite, rel=1e-6, abs=1e-8)


def test_radius_partials_match_coordinate_finite_differences():
    family = n8_family()
    z = family.collocation_nodes(3)
    rho = 1.015
    jets = radius_partial_jets(family, rho, z)
    step = 1e-8
    plus = radius_partial_jets(family, rho + step, z)
    minus = radius_partial_jets(family, rho - step, z)
    assert jets['r_rho'] == pytest.approx((plus['r'] - minus['r']) / (2.0 * step), rel=1e-6, abs=1e-8)
    plus_z = radius_partial_jets(family, rho, z + step)
    minus_z = radius_partial_jets(family, rho, z - step)
    assert jets['r_z'] == pytest.approx((plus_z['r'] - minus_z['r']) / (2.0 * step), rel=1e-6, abs=1e-8)
    assert jets['r_rhoz'] == pytest.approx(
        (plus_z['r_rho'] - minus_z['r_rho']) / (2.0 * step), rel=1e-6, abs=1e-8)
    sigma = radius_partial_jets(family, 1.0, z)
    assert np.allclose(sigma['r_z'], 0.0)
    assert np.allclose(sigma['r'] - sigma['r_ref'], 0.0)


def test_reference_history_gives_zero_h_z_and_current_history_does_not():
    z = n8_family().collocation_nodes(3)
    reference = integrate_imag_major_A2_z(
        zero_family(), z, mass=MASS, angular=ANGULAR, gauss_nodes=8)
    current = integrate_imag_major_A2_z(
        n8_family(), z, mass=MASS, angular=ANGULAR, gauss_nodes=8)
    assert np.array_equal(reference['h_z'], np.zeros((2, 3)))
    assert np.max(np.abs(current['h_z'])) > 1e-6
    assert current['is_full_e_minus2_coefficient'] is False
    assert current['numerical_C_M'] is None
    assert current['quadrature']['certified_quadrature_error_bound'] is None
    opposite = integrate_imag_major_A2_z(
        n8_family(), z, mass=MASS, angular=-ANGULAR, gauss_nodes=8)
    assert opposite['h_z'] == pytest.approx(-current['h_z'], abs=0.0, rel=0.0)
    skipped = integrate_imag_major_A2_z(
        n8_family(), z, mass=0.0, angular=ANGULAR)
    assert skipped['quadrature']['skipped'] is True
    assert np.array_equal(skipped['h_z'], np.zeros((2, 3)))


def test_h_z_matches_finite_difference_of_the_h_integral():
    family = n8_family()
    z = np.array([family.center])
    step = 1e-6
    sample = np.array([family.center - step, family.center, family.center + step])
    field = integrate_imag_major_A2(
        family, sample, mass=MASS, angular=ANGULAR, gauss_nodes=16)
    slope = integrate_imag_major_A2_z(
        family, z, mass=MASS, angular=ANGULAR, gauss_nodes=16)
    finite = (field['h'][:, 2] - field['h'][:, 0]) / (2.0 * step)
    assert slope['h_z'][:, 0] == pytest.approx(finite, rel=1e-4, abs=1e-8)


def test_current_0b0e4ced_three_node_geometry_pilot():
    family = current_family()
    identity = profile_identity(family)
    assert identity.startswith(CURRENT_HISTORY_IDENTITY_PREFIX)
    z = family.collocation_nodes(3)
    assert z.shape == (3,)
    started = time.process_time()
    result = integrate_imag_major_A2_z(
        family, z, mass=MASS, angular=ANGULAR, rho_up=ARCHIVE_RHO_UP, gauss_nodes=48)
    cpu = time.process_time() - started
    assert result['profile_identity'] == identity
    assert result['h_z'].shape == (2, 3)
    assert np.max(np.abs(result['h_z'])) > 1e-6
    assert result['is_full_e_minus2_coefficient'] is False
    assert result['numerical_C_M'] is None
    assert cpu < 120.0
    finer = integrate_imag_major_A2_z(
        family, z, mass=MASS, angular=ANGULAR, rho_up=ARCHIVE_RHO_UP, gauss_nodes=64)
    assert np.max(np.abs(result['h_z'] - finer['h_z'])) < 1e-9
    reference = integrate_imag_major_A2_z(
        zero_family(32), z, mass=MASS, angular=ANGULAR, rho_up=ARCHIVE_RHO_UP, gauss_nodes=8)
    assert np.array_equal(reference['h_z'], np.zeros((2, 3)))


def test_profile_jets_record_reconstruction_error_across_cutoff_regions():
    pytest.importorskip('flint')
    family = current_family()
    center = family.center
    inner = family.axial_inner
    outer = family.axial_outer
    samples = {
        'interior': np.array([center]),
        'ramp': np.array([center + 0.5 * (inner + outer)]),
        'exterior': np.array([center + outer + 0.02]),
        'collocation': family.collocation_nodes(3),
    }
    seen = set()
    for label, z in samples.items():
        for order in (0, 1):
            rows = chebyshev_profile_point_enclosures(family, z, derivative_order=order)
            for row in rows:
                assert row['profile'] in ('w', 'U')
                for point in row['points']:
                    seen.add(point['region'])
                    assert point['reconstruction_abs_error'] >= 0.0
                    assert point['directed_distance_to_cutoff_enclosure'] >= 0.0
                    if label == 'exterior':
                        assert point['region'] == 'exterior'
                        assert point['numerical_value'] == 0.0
                        assert point['analytic_zero_contained'] is True
                        assert point['polynomial_only_ignores_cutoff'] is True
                    if label == 'ramp':
                        assert point['region'] == 'ramp'
                        assert point['polynomial_only_ignores_cutoff'] is True
                        assert abs(point['numerical_value'] - point['polynomial_only_mid']) > (
                            1e6 * (point['cutoff_enclosure_rad'] + 1e-18))
                    if label == 'interior':
                        assert point['region'] == 'interior'
                        assert point['polynomial_only_ignores_cutoff'] is False
    assert seen == {'interior', 'ramp', 'exterior'}
    # Binary64 evaluation is not an analytic target; containment of that
    # float in the 80-bit ball is recorded and is not a pass criterion.


def test_frozen_v2_record_bytes_are_unchanged():
    path = ROOT / V2_RECORD
    if not path.exists():
        pytest.skip('historical v2 record is not present')
    digest = __import__('hashlib').sha256(path.read_bytes()).hexdigest()
    assert digest == 'de8e67e4cd89eb0331aff4e18d577220b83ddca94c0c6bac6d0c5a0fe32585df'
    payload = json.loads(path.read_text())
    assert payload['schema'] == 'NSC-KS-CURRENT-UV-COEFFICIENTS-v2'
