"""Exact Horner identity, actual interval coverage and no-generator resume."""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

import recursive_horizons.nsc_incoming_high_band_batch_source as owner
from recursive_horizons.nsc_incoming_window_refinement import window_kernel

ROOT=Path(__file__).resolve().parents[1]


def script_module():
    spec=importlib.util.spec_from_file_location('high_band_source_replay',ROOT/'scripts/derive_nsc_incoming_high_band_batch_source.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def test_exact_horner_identity_and_frozen_kernel_convention():
    assert owner.horner_identity()==['0','0','0']
    with mp.workdps(70):
        c=[mp.mpc(1,j+1)/(j+1) for j in range(24)]
        for kind,E in (('low_vacuum',mp.mpf(35)),('middle_correction',mp.mpf(75))):
            actual=owner.horner_kernel(E,mp.mpf('1.2'),mp.mpf('-2.3'),c,kind,dps=70)
            original=window_kernel(E,mp.mpf('1.2'),mp.mpf('-2.3'),c,kind)
            assert max(abs(a-b) for a,b in zip(actual,original))<mp.mpf('1e-60')
            assert actual[2]==0.


def test_group6_first_and_group13_no_overlapping_low():
    assert owner.GROUPS[0]==6 and len(owner.GROUPS)==28
    assert set(owner.GROUPS)==set(range(1,33))-{12,14,22,32}
    cases=owner.authenticated_cases(str(ROOT),6)
    assert set(cases)=={'low','middle'}
    assert (cases['low']['left'],cases['low']['right'])==(32.,40.)
    early=owner.authenticated_cases(str(ROOT),13)
    assert set(early)=={'middle'}
    assert (early['middle']['left'],early['middle']['right'])==(16.,160.)
    assert early['middle']['signs']==(1,)
    for group in (12,14,22,32):
        with pytest.raises(ValueError,match='remaining groups'):
            owner.authenticated_cases(str(ROOT),group)


def test_complete_record_preserves_signed_contributions_and_scope():
    record=json.loads((ROOT/'results/development/nsc-incoming-high-band-batch-source.json').read_text())
    assert not record['failures'] and len(record['groups'])==28 and len(record['windows'])==55
    assert record['group6_gate_payload']['group']==6
    assert 'group13_low' not in record['windows']
    assert record['physical_order24_error_bound'] is None
    assert record['scope']['new_radial_mode_or_scattering_runs']==0
    for item in record['windows'].values():
        assert not item['failures']
        assert set(item['per_sign'])=={str(s) for s in item['actual_angular_signs']}
        assert np.allclose(np.array(item['old_source'])+item['explicit_source_delta'],item['updated_source_approximant'],rtol=0,atol=3e-22)
        assert all(v>0 for v in item['thermal_scattering_stress_upper'])
        assert item['physical_order24_error_bound'] is None


def test_per_group_resume_and_full_replay_never_prepare(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('preparation entered during authenticated resume/replay')
    for name in ('prepare_group_source','channel_coefficients','cached_grid'):
        monkeypatch.setattr(owner,name,forbidden)
    script=script_module()
    prior=json.loads((ROOT/script.OUTPUT).read_text())
    cached,group6=script.prepare_one(6)
    assert cached==prior['group6_gate_payload'] and not group6['failures']
    assert script.make_record(prior['payloads'])==prior
