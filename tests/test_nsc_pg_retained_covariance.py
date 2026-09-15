"""Focused checks for signed assembly, source provenance and CTP normalization."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
sys.path.insert(0,str(ROOT/'src'))
from recursive_horizons.nsc_pg_retained_covariance import positive_panel,signed_covariance
from recursive_horizons.nsc_pg_packet_modes import SIGNED_PACKET_MAP as S


def test_signed_pair_keeps_independent_gram_and_both_angular_families():
    # Algebraic implementation control; none of these is a physical state.
    G=np.diag(np.linspace(.2,.4,8));Go=np.diag(np.linspace(.1,.3,8))
    K=.1*G;Ko=-.2*Go
    r=signed_covariance((G,K),(Go,Ko))
    assert np.allclose(r['car_gram'],G+S@Go.conj()@S,rtol=0,atol=1e-15)
    assert np.allclose(r['covariance'],.5*r['car_gram']+K-S@Ko.conj()@S,rtol=0,atol=1e-15)
    assert np.linalg.norm(r['car_gram']-np.eye(8))>.1
    assert np.allclose(1j*(r['greater']-r['lesser']),r['car_gram'],atol=1e-15)
    with pytest.raises(ValueError):signed_covariance((G,K),(Go,np.full((8,8),np.nan)))


def test_source_domain_rejects_seed_and_noncanonical_covariance():
    F=np.zeros((2,8,3),complex);w=np.ones(2)
    with pytest.raises(ValueError,match='source-fiber'):
        positive_panel(F,w,np.zeros((2,2,2)),np.broadcast_to(np.eye(2),(2,2,2)))
    P=np.broadcast_to(np.eye(3),(2,3,3)).copy()
    with pytest.raises(ValueError,match='CAR interval'):positive_panel(F,w,1.1*P,P)
    with pytest.raises(ValueError,match='accompany'):positive_panel(F,w,P)
    with pytest.raises(ValueError,match='positive finite'):positive_panel(F,-w)


def test_complete_collection_reproduces_every_retained_group_without_old_generators():
    from derive_nsc_pg_retained_covariance import read,assemble
    r=json.loads((ROOT/'results/development/nsc-pg-retained-covariance.json').read_text())
    a=read(ROOT/r['payload']['path']);m=json.loads(a['metadata_json'].tobytes())
    results,checks=assemble(a,m)
    assert {g for g,s in results}==set(range(33))
    assert len(results)==63
    for key,v in results.items():
        C=v['covariance'];G=v['car_gram']
        assert np.linalg.norm(G-np.eye(len(G)),2)<3e-9
        assert np.linalg.eigvalsh(C).min()>-3e-9
        assert np.linalg.eigvalsh(C).max()<1+3e-9
        assert np.linalg.norm(C[:4,4:])>1e-4
        assert np.linalg.norm(1j*(v['greater']-v['lesser'])-G)<3e-14
    assert r['gate']['metric_timestep'] is False
    assert r['gate']['two_sided_Weyl_mismatch'] is None
    assert r['domain']['closed_packet_system'] is False
    assert m['seed_covariance_used_as_input'] is False
