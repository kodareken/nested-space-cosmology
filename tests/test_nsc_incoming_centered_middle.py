"""Retained-order scope and replay checks for the group12 coefficient pilot."""
import importlib.util
import json
from pathlib import Path

import pytest

from recursive_horizons.nsc_incoming_centered_middle import centered_middle_bound
from recursive_horizons.nsc_incoming_projector_energy_bound import projected_energy_bound

ROOT = Path(__file__).resolve().parents[1]


def test_incompatible_domain_is_rejected_before_any_radial_work():
    config = {}; channel = {'index': 14}
    with pytest.raises(ValueError, match='other than14'):
        centered_middle_bound(config, channel)
    with pytest.raises(ValueError, match='fixed128'):
        centered_middle_bound(config, {'index': 12}, intervals=16)


def test_record_replays_saved_coefficients_and_keeps_all_equations_open(monkeypatch):
    record = json.loads((ROOT/'results/development/nsc-incoming-centered-middle.json').read_text())
    spec = importlib.util.spec_from_file_location('centered_middle_replay', ROOT/'scripts/derive_nsc_incoming_centered_middle.py')
    script = importlib.util.module_from_spec(spec); spec.loader.exec_module(script)
    def forbidden(*args, **kwargs): raise AssertionError('replay must not generate radial coefficients')
    monkeypatch.setattr('recursive_horizons.nsc_incoming_centered_middle.centered_middle_bound', forbidden)
    actual = script.replay(ROOT/record['payload']['path'])
    assert actual == record
    assert not actual['scope']['source_values_changed']
    assert not actual['scope']['constraints_solved']
    assert not actual['scope']['physical_NON_EXISTENCE_claimed']
    assert actual['same_numerical_order'] == 16
    assert actual['stationarity_tolerance'] == 3e-11


def test_saved_intersections_improve_valid_natural_coefficients():
    record = json.loads((ROOT/'results/development/nsc-incoming-centered-middle.json').read_text())
    artifact = json.loads((ROOT/record['payload']['path']).read_text())
    control = artifact['controls']['centered']; radial = control['radial']
    channel = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][12]
    assert projected_energy_bound(channel, radial, control['local_cross'], 40., 320.) == control['energy']
    for row in radial['per_sign']:
        assert len(row['I16_to_I32_upper']) == 17
        assert all(0 <= x <= y for x, y in zip(row['I16_to_I32_upper'], row['natural_I16_to_I32_upper']))
