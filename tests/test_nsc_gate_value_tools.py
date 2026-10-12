#!/usr/bin/env python3
"""Value-only bookkeeping, the Dirac lock, and the declared merit."""
import os
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n64_truncated_iterate as N64
import derive_nsc_ks_gate_step as Step
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from recursive_horizons.nsc_ks_local_constraints import _family_matter
from recursive_horizons.nsc_ks_family_pool import (
    acquire_dirac_lock, echo_item, map_families, other_dirac_runs, release_dirac_lock)
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_source_envelope import (
    computational_z_grid, physical_incoming_interval, usual_axial_support)
from recursive_horizons.nsc_ks_value_evaluator import (
    MERIT_DEFINITION, assemble_residual, geometry_change, matter_change, merit,
    mutated_matter, target_nodes)
from recursive_horizons.nsc_local_incoming_constraints import _channel_record
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from test_nsc_ks_local_constraints import coefficients_for


def _prepared(tangents):
    source = FixedSourcePreparation.from_signed_blocks(
        np.array([0.7]), np.array([1.2]), np.array([[[0.6, 0.04j], [-0.04j, 0.4]]]))
    initial = np.array([[0.6, 0.3j], [0.2j, 0.7]], complex)
    coefficients = np.zeros((2, 8))
    coefficients[0, 1] = 0.002
    coefficients[1, 0] = 0.001
    family = LocalIncomingFamily(coefficients)
    prepared = evolve_ks_difference_envelope(
        source, initial, family, computational_z_grid(16),
        np.linspace(*physical_incoming_interval(), 5), np.pi / 2, np.sqrt(5), 1.03,
        axial_support=usual_axial_support(), rtol=2e-10, atol=2e-12, max_step=0.002,
        tangents=tangents)
    channel = _channel_record('synthetic', {'synthetic': {
        'index': 14, 'compact_mass': np.pi / 2, 'angular_eigenvalue': np.sqrt(5),
        'angular_sign': 1, 'copy_count': 2, 'degeneracy': 12,
    }})
    return family, prepared, channel


def test_sixty_four_coefficients_match_the_script_family():
    low = LocalIncomingFamily(np.zeros((2, 8)))
    again = LocalIncomingFamily(np.zeros((2, 8)))
    assert profile_identity(low, include_normal_window=True) == profile_identity(
        again, include_normal_window=True)
    coefficients = np.zeros((2, 64))
    owned = LocalIncomingFamily(coefficients)
    script = N64.LocalIncomingFamily64(coefficients)
    assert owned.radius_lower_bound() == script.radius_lower_bound()
    assert owned.description()['coefficients'] == script.description()['coefficients']
    assert owned.radius_lower_bound() > 0


def test_target_counts_and_merit_definition():
    family = LocalIncomingFamily(np.zeros((2, 8)))
    for count in (47, 129, 257):
        assert len(target_nodes(family, count)) == count
    assert merit([[1e-4, -2e-4], [3e-6, 1e-8]]) == 2e-4
    assert MERIT_DEFINITION.startswith('m(g)=')


def test_zero_tangents_agree_with_tangent_mode_on_a_synthetic_family():
    _family, zero, channel = _prepared('zero')
    _family, full, _channel = _prepared('all')
    coefficients = coefficients_for()
    assert zero.column_tangents.shape[0] == 0
    change = matter_change(zero, channel, coefficients)
    tangent_change = _family_matter(
        full, channel, coefficients, full.column_tangents.shape[0], np.asarray(full.z))['correction']
    assert change.shape == tangent_change.shape
    assert np.max(np.abs(change - tangent_change)) < 1e-8


def test_mutations_do_not_reproduce_the_honest_matter_change():
    _family, prepared, channel = _prepared('zero')
    coefficients = coefficients_for()
    honest = matter_change(prepared, channel, coefficients)
    for name in ('drop-delta-C', 'drop-coherence', 'double-weights'):
        mutated = mutated_matter(prepared, channel, coefficients, name)
        assert not np.allclose(mutated, honest)
    with pytest.raises(ValueError, match='k=-E'):
        mutated_matter(prepared, channel, coefficients, 'k=-E')


def test_geometry_is_added_once():
    coefficients = np.zeros((2, 8))
    coefficients[1, 0] = 0.02
    family = LocalIncomingFamily(coefficients)
    z = target_nodes(family, 47)
    geometry = geometry_change(family, z, coefficients_for())
    matter = np.zeros_like(geometry)
    edge = np.zeros_like(geometry)
    once = assemble_residual(np.array([0.1, -0.2]), geometry, matter, edge)
    twice = assemble_residual(np.array([0.1, -0.2]), 2 * geometry, matter, edge)
    assert not np.allclose(once, twice)


def test_dirac_lock_refuses_a_live_holder_and_another_driver(tmp_path):
    lines = [f'{os.getpid()} python3 scripts/derive_nsc_ks_gate_value.py --run']
    assert other_dirac_runs(lines, os.getpid()) == []
    assert len(other_dirac_runs(lines, os.getpid() + 1)) == 1
    assert other_dirac_runs(
        ['99 python3 scripts/derive_nsc_ks_n64_primal_norm_trust03.py --check'],
        os.getpid()) == []
    lock = tmp_path/'nsc-ks-gate-dirac.lock'
    held = acquire_dirac_lock(lock, 'test', os.getpid(), lines=[])
    try:
        with pytest.raises(RuntimeError, match='live pid'):
            acquire_dirac_lock(lock, 'other', os.getpid() + 1, lines=[])
    finally:
        release_dirac_lock(held, os.getpid())
    assert not lock.exists()


def test_pool_returns_items_in_submission_order():
    assert map_families(echo_item, (3, 1, 2), workers=2) == (3, 1, 2)


def test_ledger_declares_merit_and_appends_without_rewriting(tmp_path):
    path = tmp_path/'ledger.json'
    Step.declare_ledger(path)
    row = {
        'parent': 'b4b53f0b', 'rule': 'unit-test', 'parameters': {'nodes': 47},
        'step_size': 0.0, 'predicted_merit': None, 'measured_merit': 1e-6,
        'accepted': False,
    }
    Step.append_row(row, path)
    Step.append_row({**row, 'parent': 'second'}, path)
    ledger = Step.load_ledger(path)
    assert ledger['merit_definition'] == MERIT_DEFINITION
    assert [item['parent'] for item in ledger['rows']] == ['b4b53f0b', 'second']
    with pytest.raises(ValueError, match='forbidden'):
        Step.append_row({**row, 'candidate_profile_identity': Step.FORBIDDEN[0]}, path)
