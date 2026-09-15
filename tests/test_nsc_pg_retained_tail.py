"""Tail-only retained massive inputs; not a full covariance."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from derive_nsc_pg_retained_tail import (
    C_TOL,EXCLUDED,G_TOL,OUTPUT,families,load_artifact,make_record,verify_arrays,_run,
)
from recursive_horizons.nsc_pg_packet_modes import SIGNED_PACKET_MAP


def _loaded():
    record=json.loads((ROOT/OUTPUT).read_text())
    arrays,metadata,_=load_artifact(ROOT/record['payload']['path'])
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    return record,arrays,metadata,state


def test_family_enumeration_is_sixty_one_signed_labels():
    state=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    rows=families(state)
    assert [r[0] for r in rows].count(EXCLUDED)==0
    assert {r[1] for r in rows if r[0]==23}=={1}
    nonzero=[r for r in rows if r[0]!=23]
    assert {r[0] for r in nonzero}==set(range(1,33))-{13,23}
    for g in {r[0] for r in nonzero}:
        assert {r[1] for r in nonzero if r[0]==g}=={1,-1}
    assert len(rows)==61


def test_retained_tail_record_is_tail_only_and_covers_all_families():
    record,arrays,metadata,state= _loaded()
    verify_arrays(arrays,metadata,state)
    assert record['gate']['full_C1b']=='OPEN'
    assert record['gate']['tail_inputs'] in ('PASS','OPEN')
    assert record['gate']['full_covariance_assigned'] is False
    assert record['gate']['signed_fold'] is False
    assert record['gate']['A_q_Omega_zeta_Vfull_changed'] is False
    assert record['domain']['families']==61
    assert 13 not in record['coverage']['groups']
    assert record['coverage']['family_count']==61
    restart=json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())
    assert record['locked_inputs']==restart['locked_inputs']
    rebuilt,_=make_record(ROOT/record['payload']['path'],arrays,metadata)
    assert json.dumps(record,sort_keys=True)==json.dumps(rebuilt,sort_keys=True)


def test_positive_energy_matrices_are_not_signed_folded_or_fabricated():
    record,arrays,metadata,_=_loaded()
    S=SIGNED_PACKET_MAP
    assert arrays['centered'].shape==(61,8,8) and arrays['gram'].shape==(61,8,8)
    i=int(np.flatnonzero((arrays['group']==1)&(arrays['sign']==1))[0])
    j=int(np.flatnonzero((arrays['group']==1)&(arrays['sign']==-1))[0])
    K=arrays['centered'][i];folded=K-S@K.conj()@S
    assert np.linalg.norm(K-folded)>1e-18
    assert np.linalg.norm(arrays['centered'][i]-arrays['centered'][j])>1e-18
    for g,L in record['coverage']['selected_lower'].items():
        idx=np.flatnonzero(arrays['group']==int(g))
        assert set(arrays['selected_lower'][idx])=={float(L)}
        c=float(np.sum(arrays['bound_C'][idx]));gsum=float(np.sum(arrays['bound_G'][idx]))
        assert c==pytest.approx(record['coverage']['C_bound_sum'][g],rel=1e-7,abs=0)
        if record['coverage']['all_pass']:
            assert c<=C_TOL and gsum<=G_TOL
    for g in record['groups_requiring_320']:
        assert g in set(map(int,arrays['diag160_group']))
        assert record['coverage']['selected_lower'][str(g)]==320
    assert all(item['group']!=13 for item in record['mode_order_12_vs_16'])


def test_worker_cap_is_two():
    with pytest.raises(ValueError,match='at most 2'):
        _run([],3,'cap')
