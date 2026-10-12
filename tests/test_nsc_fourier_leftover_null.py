"""Leftover-null and the compact cokernel stay short of a class certificate."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_s5_is_retracted_and_stall_is_finite_image():
    retract = json.loads(
        (ROOT / 'results/development/nsc-local-incoming-gate-status-retract-s5.json').read_text())
    assert retract['s5_retracted'] is True
    assert retract['class_rewrite_is_not_next'] is True
    assert retract['finite_leftover_is_not_class_identity'] is True
    assert retract['named_search_stall'] == 'finite_wu_geometry_image_cannot_reach_existence'
    assert 's^5' not in (retract['next_named_owner'] or '')
    assert retract['previous_named_search_stall'] == 'declared_wu_image_cannot_reach_existence'
    assert 's^5' in retract['previous_next_named_owner']
    assert retract['physical_local_gate'] == 'OPEN'
    assert retract['physical_NONEXISTENCE_certificate'] is False


def test_fourier_leftover_null_is_not_a_class_identity():
    leftover_null = json.loads(
        (ROOT / 'results/development/nsc-ks-fourier-leftover-null.json').read_text())
    assert leftover_null['not_n16_walk_v'] is True
    assert leftover_null['L_E_not_reevaluated'] is True
    assert leftover_null['leftover_null_norm'] == 1.0
    assert leftover_null['physical_NONEXISTENCE_certificate'] is False
    assert leftover_null['unenclosed_edge_tail'] is None
    assert leftover_null['scope']['missing_tail_blocks_gap'] is True
    assert leftover_null['scope']['finite_span_is_not_the_class'] is True
    assert leftover_null['scope']['new_field_or_source_evolutions'] == 0
    assert leftover_null['selected_mode_count'] in (8, 16)
    assert len(leftover_null['walk']) == 5


def test_continuous_cokernel_is_trivial_and_not_L_or_flux():
    cokernel = json.loads(
        (ROOT / 'results/development/nsc-ks-declared-wu-cokernel.json').read_text())
    assert cokernel['continuous_cokernel_trivial'] is True
    assert cokernel['is_L_or_flux_again'] is False
    assert cokernel['recovered_L_E'] is False
    assert cokernel['recovered_flux_divergence'] is False
    assert cokernel['principal_det_excludes_zero'] is True
    assert cokernel['certified_det_excludes_zero'] is True
    assert cokernel['compact_support_forces_zero_multiplier'] is True
    assert cokernel['physical_NONEXISTENCE_certificate'] is False
    assert cokernel['geometry_killed'] is False
    assert cokernel['unenclosed_edge_tail'] is None


def test_declared_class_gate_stays_open_and_paper_blocked():
    status = json.loads(
        (ROOT / 'results/development/nsc-local-incoming-gate-status-declared-class.json').read_text())
    bounds = json.loads(
        (ROOT / 'results/development/nsc-ks-deciding-bounds-after-cokernel.json').read_text())
    blocked = json.loads(
        (ROOT / 'results/development/nsc-027-blocked-declared-class.json').read_text())
    assert status['physical_local_gate'] == 'OPEN'
    assert status['physical_EXISTENCE_certificate'] is False
    assert status['physical_NONEXISTENCE_certificate'] is False
    assert status['s5_retracted'] is True
    assert status['named_search_stall'] == 'finite_wu_geometry_image_cannot_reach_existence'
    assert 's^5' not in (status['next_named_owner'] or '')
    assert 'no class rewrite' in status['next_named_owner']
    assert bounds['bounds_spent'] is False
    assert bounds['changed_history_UV_tail'] is None
    assert bounds['field_error_bound'] is None
    assert bounds['baseline_low_subgap'] is None
    assert blocked['manuscript_027_written'] is False
    assert blocked['arxiv_metadata_prepared'] is False
    assert blocked['public_checkout_unused'] is True
    assert blocked['named_search_stall'] == status['named_search_stall']
