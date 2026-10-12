"""Checkpoint coverage and directed accumulation; no field runs."""
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest
pytest.importorskip('flint')
from flint import arb,ctx

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import validate_nsc_ks_fine_history as V
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper,restored_upper


def fixture(monkeypatch,tmp_path):
    monkeypatch.setattr(V,'ROOT',tmp_path)
    directory=tmp_path/'proofs';directory.mkdir()
    monkeypatch.setattr(V,'DIRECTORY',directory)
    owner=SimpleNamespace(total=2,signature={'fixture':'yes'},signature_digest='fixture',
        segment=lambda i:SimpleNamespace(rho_start=1.03-.01*i,rho_end=1.02-.01*i))
    with ctx.workprec(120):
        row={'segment':0,'rho_start':1.03,'rho_end':1.02,'CPU_seconds':.1,'signature_digest':'fixture',
             'continuous_residual_integral_upper':[exact_upper(arb(1)/100),exact_upper(arb(1)/50)]}
    (directory/'segment-000000.json').write_text(json.dumps(row))
    return owner,row


def test_partial_coverage_cannot_create_a_full_history_error(monkeypatch,tmp_path):
    owner,row=fixture(monkeypatch,tmp_path)
    rows=V.saved_rows(owner)
    result=V.summary(owner,rows)
    assert result['validated_segments']==1 and not result['all_segments_covered']
    assert result['numerical_field_error_relative_to_saved_upstream'] is None
    assert result['physical_local_gate']=='OPEN'


def test_stale_or_relocated_cell_is_rejected(monkeypatch,tmp_path):
    owner,row=fixture(monkeypatch,tmp_path)
    path=V.DIRECTORY/'segment-000000.json'
    row['rho_end']=1.019
    path.write_text(json.dumps(row))
    with pytest.raises(ValueError,match='another rho cell'):V.saved_rows(owner)
    row['rho_end']=1.02;row['signature_digest']='stale'
    path.write_text(json.dumps(row))
    with pytest.raises(ValueError,match='owner/input changed'):V.saved_rows(owner)
