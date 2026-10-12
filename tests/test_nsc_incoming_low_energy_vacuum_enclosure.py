"""New interval algebra and authenticated error-sum replay, never a second run."""
import gzip
import importlib.util
import json
from pathlib import Path
import mpmath as mp
import pytest
import recursive_horizons.nsc_incoming_low_energy_vacuum_enclosure as M
ROOT=Path(__file__).resolve().parents[1]


def test_exact_frozen_H_interval_update_and_point_rounding():
    with M._precision(M.PRECISION):
        D,x,y=map(mp.mpf,('.7','.2','-.3'));h=mp.iv.mpf(1)/100
        column=(mp.mpc(1),mp.mpc(0))
        after,image,rounding=M.frozen_update(column,D,x,y,h)
        n=mp.sqrt(D*D+x*x+y*y);s=mp.sin(n/100)/n;c=mp.cos(n/100)
        exact=(c+mp.j*s*D,-mp.j*s*(x+mp.j*y))
        for v,box in zip(exact,image):
            assert M._lo(box.real)<=v.real<=M._hi(box.real)
            assert M._lo(box.imag)<=v.imag<=M._hi(box.imag)
        distance=mp.sqrt(sum(abs(a-b)**2 for a,b in zip(after,exact)))
        assert distance<=M._hi(rounding)
        assert abs(sum(abs(v)**2 for v in after)-1)<mp.mpf('1e-55')


def test_whole_cell_norm_and_start_tail_are_outward():
    channel=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels'][14]
    config=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    with M._precision(M.PRECISION):
        G=M.SelectedGenerator(channel,config)
        y=M._lo(mp.iv.ln(mp.iv.mpf(M.DELTA0)));tail=G.tail(y)
        assert M.unpoint(tail['tail_delta_upper'])>=M._hi(mp.iv.exp(M.ivpoint(y)))
        assert M._hi(M.unpack(tail['affine_projector_tail_error_upper']))>0
        left,right=mp.mpf('-2'),mp.mpf('-1.99')
        whole=G.components(M._range(left,right))
        for point in (left,(left+right)/2,right):
            values=G.components(M.ivpoint(point))
            for bound,value in zip(whole,values):
                assert M._lo(bound)<=M._lo(value)<=M._hi(value)<=M._hi(bound)
        # Frozen Hermitian difference norm is exact for a varying diagonal.
        eta=M.cell_defect((mp.iv.mpf([.5,1]),mp.iv.mpf(0),mp.iv.mpf(0)),(mp.mpf('.75'),mp.mpf(0),mp.mpf(0)),mp.iv.mpf(1)/16)
        assert M._hi(eta)>=mp.mpf(1)/64
        with pytest.raises(ValueError):M.SelectedGenerator(channel,config,24.)


def test_receipt_replays_without_generator_or_frozen_updates(monkeypatch):
    spec=importlib.util.spec_from_file_location('low_energy_proof',ROOT/'scripts/derive_nsc_incoming_low_energy_vacuum_enclosure.py')
    script=importlib.util.module_from_spec(spec);spec.loader.exec_module(script)
    def forbidden(*args,**kwargs):raise AssertionError('propagation entered error replay')
    monkeypatch.setattr(M,'propagate',forbidden);monkeypatch.setattr(M,'frozen_update',forbidden)
    monkeypatch.setattr(M.SelectedGenerator,'components',forbidden)
    record=script.check()
    if record['result'] is not None:
        assert record['result']['lossless_step_replay_residual']==0
        assert not record['result']['column_normalized']
        assert record['status'].startswith('PASS')==record['result']['strict_Re_P01_above_one_quarter']
    assert not record['scope']['coherent_source_replaced']
    assert not record['scope']['group32_bound_applied_to_low_energy']
    assert record['algorithm']['steps']==2048


def test_replay_rejects_incomplete_cell_coverage_and_precision():
    saved=json.loads((ROOT/'results/development/nsc-incoming-low-energy-vacuum-enclosure.json').read_text())
    if saved['result'] is None:pytest.skip('budget checkpoint contains no completed certificate')
    payload=json.loads(gzip.decompress((ROOT/saved['payload']['path']).read_bytes()))
    wrong=dict(payload,precision=50)
    with pytest.raises(ValueError,match='fixed-size'):M.replay(wrong)
    wrong=dict(payload,rows=payload['rows'][:-1])
    with pytest.raises(ValueError,match='fixed-size'):M.replay(wrong)
