"""Same-source convention, explicit remainders and no-generator replay."""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

import recursive_horizons.nsc_incoming_low_high_source as owner
from recursive_horizons.nsc_incoming_source_tail import _riccati_at_one, paired_source_function


ROOT = Path(__file__).resolve().parents[1]


def test_zero_extended16_source_matches_existing_ad4_owner():
    # Finite algebra control of signs, subtraction and normalization. No
    # physical mode, scattering or radial interval calculation is called.
    channel = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][22]
    with mp.workdps(70):
        mass, ell, E = mp.mpf(channel['compact_mass']), mp.mpf(channel['angular_eigenvalue']), mp.mpf(35)
        reference = paired_source_function(channel)(1/E)
        actual = mp.matrix([0, 0, 0, 0])
        for sign in (1, -1):
            coeff = _riccati_at_one(mass, sign*ell, order=16)+[mp.mpc(0)]*8
            actual += owner.vacuum_source_kernel(E, mass, sign*ell, coeff)
        assert max(abs(actual[i]-reference[j]) for i, j in ((0, 0), (1, 1), (3, 2))) < mp.mpf('1e-60')
        assert actual[2] == 0
        with pytest.raises(ValueError, match='order24'):
            owner.vacuum_source_kernel(E, mass, ell, coeff[:16])
        with pytest.raises(ValueError, match='LOW'):
            owner.vacuum_source_kernel(31, mass, ell, coeff)


def test_original_source_cells_and_single_group_scope():
    owned = owner.authenticated_low_inputs(ROOT)
    assert owned['channel']['index'] == 22
    assert owned['measure_residual'] <= 3e-13
    for sign in (1, -1):
        for points in (24, 32):
            row = owned['selected'][f'{sign}/old{points}']
            assert len(row['energies']) == points
            assert np.all(row['weights'] > 0)
            assert 32 < row['energies'].min() < row['energies'].max() < 40
    wrong = dict(owned['channel'], index=14)
    with pytest.raises(ValueError, match='group22'):
        owner.prepare_group22_defect(owned['config'], wrong)


def test_record_preserves_old_values_and_nonzero_thermal_remainder():
    record = json.loads((ROOT/'results/development/nsc-incoming-low-high-source.json').read_text())
    assert record['group'] == 22 and record['interval'] == [32., 40.]
    assert not record['failures']
    assert record['numerical_Riccati_order'] == record['radial_certificate']['physical_Riccati_order'] == 24
    assert record['radial_certificate']['intervals'] == 128
    assert record['radial_certificate']['centered_remainder_depth'] == 4
    assert np.array_equal(np.array(record['new_order24_vacuum_source'])-record['archived_LOW32_source'],
                          record['new_vacuum_minus_archived_LOW32'])
    assert record['new_order24_vacuum_source'][2] == 0.
    assert all(v > 0 for v in record['omitted_thermal_stress_upper'])
    assert not record['scope']['thermal_source_physically_zero']
    assert not record['scope']['old_finite_offset_accuracy_certified']
    assert not record['scope']['old_difference_identified_as_inverse_noise']
    assert not record['physical_error_budget']['rigorous_numerical_source_quadrature_error_included']
    assert record['physical_error_budget']['stationarity_tolerance'] == 3e-11


def test_replay_does_not_run_interval_or_source_preparation(monkeypatch):
    def forbidden(*args, **kwargs): raise AssertionError('preparation entered during authenticated replay')
    for name in ('trimmed_defect_jets', '_geometry_jets', '_riccati_at_one', 'positive_grid'):
        monkeypatch.setattr(owner, name, forbidden)
    spec = importlib.util.spec_from_file_location('nsc_low_high_replay', ROOT/'scripts/derive_nsc_incoming_low_high_source.py')
    script = importlib.util.module_from_spec(spec); spec.loader.exec_module(script)
    prior = json.loads((ROOT/script.OUTPUT).read_text())
    assert script.make_record(ROOT/prior['payload']['path']) == prior
