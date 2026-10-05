"""Full-interval UV attempt: domain refusal, phase-fixed reentry, no invented C_M."""
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_ks_full_interval_uv as derive
from recursive_horizons.nsc_ks_full_interval_uv import (
    UNCONSTRUCTED_TERM,
    full_interval_uv_report,
    validate_full_interval_uv_report,
)
from recursive_horizons.nsc_ks_uv_minor_a3 import minor_a3_identities


@pytest.fixture(scope='module')
def report():
    value = full_interval_uv_report(ROOT)
    validate_full_interval_uv_report(value)
    return value


def test_live_history_full_interval_is_not_enclosed(report):
    assert report['history_identity'].startswith('0b0e4ced')
    assert report['domain']['interval_I'] == 'S(1)+[0.12,0.18]'
    assert report['domain']['dyadic_radius_exceeds_axial_inner']
    assert 'straddles the plateau center' in report['domain']['domain_failure']
    assert report['domain']['full_interval_jets_enclosed'] is False
    assert report['domain']['proper_subinterval_enclosed'] is True
    assert report['domain']['proper_subinterval_is_full_I'] is False
    assert report['domain']['subinterval_b_z_contains_zero_on_sigma'] is True
    assert report['preimage']['both_characteristics_exit_support'] is True
    assert report['preimage']['characteristic_distance_lower'] > report['domain']['axial_outer'] - report['domain']['axial_inner']
    assert report['full_interval_covered'] is False
    assert report['cpu_seconds'] < 180
    assert report['rss_mib'] < 1024


def test_phase_fixed_major_a2_is_not_the_zero_rhs_value(report):
    assert report['major_A2']['contains_zero'] is False
    assert report['opposite_angular_major_A2']['contains_zero'] is False
    assert report['major_A2']['strictly_negative'] is True
    assert report['opposite_angular_major_A2']['strictly_negative'] is True
    assert report['major_A4']['contains_zero'] is False
    assert report['upstream_order'] == 5
    assert report['normalization_residual_contains_zero'] is True
    assert report['A2_set_to_zero_from_vanishing_imaginary_rhs'] is False
    assert report['degree4_pointwise_tail_upper'] > 0
    identities = minor_a3_identities()
    assert identities['residuals']['s+1/minor_A3_depends_on_major_A2'] == 'depends'
    assert identities['massless_Im_Ts_major_A2_sets_upstream_datum_to_zero'] is False


def test_reentry_is_separated_and_is_not_a_remainder(report):
    assert report['reentry']['zero_deformation_kernel'] == 0
    assert report['reentry']['zero_rhs_reentry_majorant'] == 0
    assert report['reentry']['zero_rhs_drops_positive_reentry_majorant'] is True
    assert report['reentry']['reentry_majorant']['binary64_upper'] > 0
    assert report['reentry']['complete_minor_A3'] is None
    point = report['pointwise_counterexample']
    assert point['sample_is_full_interval'] is False
    assert point['sample_ray_meets_I'] is True
    assert point['zero_rhs_removes_separated_term'] is True
    assert point['R_up_times_delta_inv_r']['contains_zero'] is False
    assert point['delta_b']['contains_zero'] is False
    assert point['delta_b']['strictly_positive'] is True
    assert report['numerical_C_M'] is None
    assert report['induced_N_beta_remainder'] is None
    assert report['current_history_L0_A4_H2_integrals'] is None
    assert report['unconstructed_term'] == UNCONSTRUCTED_TERM
    assert report['commutators_are_C_M'] is False
    assert report['Bz_integral'] > 0 and report['Bzz_integral'] > 0
    assert report['loose_absolute_jet_used_as_L0_A4'] is False
    assert report['thermal_coherent_tail_included'] is False
    assert report['physical_local_gate'] == 'OPEN'


def test_carrier_has_the_leading_cutoff_weight_and_is_not_the_difference(report):
    carrier = report['carrier']
    assert carrier['F_z_equals_energy_times_F'] is True
    assert carrier['same_leading_cutoff_weight'] == '1/(2*Ec^2)'
    assert carrier['placeholder_period_used_as_physical_domain'] is False
    assert carrier['this_integral_is_history_minus_reference'] is False
    assert carrier['absolute_common_tail_exceeds_allocation'] is True
    assert carrier['exceeding_allocation_is_not_nonexistence'] is True
    assert carrier['absolute_common_carrier_integral_upper'] > carrier['provisional_allocation']
    assert report['moments']['angular_families'] == 60
    assert report['moments']['leading_cutoff_factor'].startswith('per family 1/(2*Ec^2)')
    assert report['moments']['moment_fold_supplies_C_M'] is False
    assert all(value > 0 for value in report['moments']['leading_moment_uppers'])


def test_record_is_refused_without_a_finite_remainder(report):
    with pytest.raises(SystemExit, match='exclusive record not published'):
        derive.publish(report)
