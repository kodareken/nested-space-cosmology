"""True retained-group routing, zero-angular factors and immutable pilot reuse."""
import importlib.util
import json
from pathlib import Path

import mpmath as mp
import pytest

from recursive_horizons.nsc_incoming_tail_quadrature_batch import endpoint,signs,prepare_tail,replay_tail,thermal_tail,RetainedTailKernel
from recursive_horizons.nsc_incoming_middle_bound import _parameters
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision,_lo,_hi
from recursive_horizons.nsc_incoming_source_quadrature_bound import unpack

ROOT=Path(__file__).resolve().parents[1]
CHANNELS=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
CONFIG=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
OUTPUT=ROOT/'results/development/nsc-incoming-tail-quadrature-batch.json'


def test_actual_endpoint_and_sign_routing():
    for g in range(1,33):
        assert endpoint(CHANNELS[g])==(320. if g in (10,11,12,31,32) else 160.)
        assert signs(CHANNELS[g])==((1,) if g in (13,23) else (1,-1))
    with pytest.raises(ValueError,match='retained'):endpoint(CHANNELS[0])
    with pytest.raises(ValueError,match='completed group22'):prepare_tail(CHANNELS[22],{})


def test_generic_zero_angular_factor_and_outward_thermal_inventory():
    with _precision(80):
        kernel=RetainedTailKernel(CHANNELS[13],1,[mp.iv.mpc(0)]*16)
        assert kernel.factor._mpi_==_parameters(CHANNELS[13])[4]._mpi_
        with pytest.raises(ValueError,match='actual sign'):RetainedTailKernel(CHANNELS[13],-1,[mp.iv.mpc(0)]*16)
        for g in (1,13,23,32):
            result=thermal_tail(CHANNELS[g],CONFIG)
            assert result['exact_zero_components']==([3] if g in (13,23) else [])
            for i,value in enumerate(result['intervals']):
                interval=unpack(value)
                assert _lo(interval)==_hi(interval)==0 if i in result['exact_zero_components'] else _lo(interval)>0


def test_pilot_payload_replays_identically_without_preparation():
    record=json.loads((ROOT/'results/development/nsc-incoming-tail-quadrature-bound.json').read_text())
    payload=json.loads((ROOT/record['payload']['path']).read_text())
    assert replay_tail(payload['integral'],CHANNELS[22])==record['certified_integral']
    assert thermal_tail(CHANNELS[22],CONFIG)==record['thermal_tail']


def test_batch_saved_only_replay_and_archive_preservation(monkeypatch):
    if not OUTPUT.exists():pytest.skip('post-preparation batch gate')
    spec=importlib.util.spec_from_file_location('tail_batch',ROOT/'scripts/derive_nsc_incoming_tail_quadrature_batch.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    def forbidden(*args,**kwargs):raise AssertionError('replay prepared scientific integral')
    monkeypatch.setattr(module,'prepare_tail',forbidden)
    monkeypatch.setattr(RetainedTailKernel,'__call__',forbidden)
    old=json.loads(OUTPUT.read_text());assert module.replay(old['manifest'])==old
    assert old['group_count']==32 and old['manifest']['new_group_count']==31
    assert old['manifest']['reused_groups']==[22]
    assert old['aggregate']['archived_source_unchanged']
    assert not old['scope']['pilot22_recomputed']
