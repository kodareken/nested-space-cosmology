"""Signed source error retains the numerical Gram/complement discrepancy."""
import numpy as np
import pytest
from flint import arb, ctx

from recursive_horizons.nsc_subgap_source_covariance import original_covariance_error
from recursive_horizons.nsc_subgap_upstream_covariance import (
    negative_subgap_covariance_error, numerical_gram_error,
)
from recursive_horizons.nsc_transmitting_dirac_domain import S1,S2,S3


def test_bloch_complement_matches_the_owned_antiunitary_map():
    vector=np.array([.3,-.4,.7])
    Q=(np.eye(2)+vector[0]*S1+vector[1]*S2+vector[2]*S3)/2
    negative=(np.eye(2)+vector[0]*S1-vector[1]*S2-vector[2]*S3)/2
    np.testing.assert_allclose(negative,np.eye(2)-S3@Q.conj()@S3,atol=1e-15)


def test_copying_positive_error_would_miss_negative_gram_error():
    A=np.array([[1,0,0],[0,2,0]],complex)
    C=np.diag([1.,0.,0.])
    proof={'bits':160,'endpoint':(arb(0),arb(0),arb(1)),
           'bloch_error':arb(0)}
    positive=original_covariance_error(A,C,proof)
    negative=negative_subgap_covariance_error(S3@A.conj(),np.eye(3)-C.conj(),proof)
    assert positive == 0
    assert negative >= 3
    assert negative < arb('3.00000000000001')
    assert numerical_gram_error(A) == 3


def test_signed_comparison_accounts_for_covariance_arithmetic():
    A=np.array([[1,0,0],[0,1,0]],complex)
    proof={'bits':160,'endpoint':(arb(0),arb(0),arb(1)),
           'bloch_error':arb(0)}
    Cminus=np.diag([1e-7,1.,0.])
    error=negative_subgap_covariance_error(S3@A,Cminus,proof)
    assert error >= arb(float(1e-7))
    assert error < arb('1.00000000001e-7')


def test_gram_domain_and_context():
    before=ctx.prec,ctx.cap
    with pytest.raises(ValueError):numerical_gram_error(np.eye(2))
    numerical_gram_error(np.array([[1,0,0],[0,1,0]],complex),bits=160)
    assert (ctx.prec,ctx.cap)==before


def test_actual_archived_signed_pair_is_bounded_without_new_preparation(monkeypatch):
    import runpy
    from pathlib import Path
    from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
    root=Path(__file__).resolve().parents[1]
    driver=runpy.run_path(str(root/'scripts/derive_nsc_subgap_upstream_covariance.py'))
    record,payload=driver['calculate']()
    assert record['source_signs']==[1,-1]
    assert record['angular_signs']==[1,-1]
    assert record['covered_rows']==[0,1]
    for key in ('positive_covariance_operator_error_upper','negative_covariance_operator_error_upper'):
        assert restored_upper(record[key])<arb('1e-11')
    assert restored_upper(record['numerical_positive_gram_frobenius_error_upper'])>0
    assert record['negative_error_copied_from_positive'] is False
    assert record['preparation_ode_rerun'] is False
    assert record['archived_upstream_columns_replaced'] is False
    assert record['all_source_families'] is False
    assert record['physical_upstream_budget_component'] is None
    assert record['physical_local_gate']=='OPEN'
    assert len(payload)<100_000
    def no_new_solve(*args,**kwargs):
        raise AssertionError('replay must consume the saved trajectory')
    monkeypatch.setitem(driver['calculate'].__globals__,'capture_bloch',no_new_solve)
    replay,replayed_payload=driver['calculate'](saved=(record,payload))
    assert replay==record and replayed_payload==payload
    from copy import deepcopy
    from io import BytesIO
    from hashlib import sha256
    from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
    with np.load(BytesIO(payload),allow_pickle=False) as saved:
        arrays={k:saved[k].copy() for k in saved.files}
    arrays['trace'][0,0]-=.01
    wrong=deterministic_npz_bytes(arrays)
    forged=deepcopy(record);forged['payload']['sha256']=sha256(wrong).hexdigest()
    with pytest.raises(ValueError,match='initializer'):
        driver['calculate'](saved=(forged,wrong))
    with pytest.raises(ValueError,match='payload hash'):
        driver['calculate'](saved=(record,wrong))
    with np.load(BytesIO(payload),allow_pickle=False) as saved:
        arrays={k:saved[k].copy() for k in saved.files}
    arrays['positive_columns'][0,0]+=.001
    wrong=deterministic_npz_bytes(arrays)
    forged=deepcopy(record);forged['payload']['sha256']=sha256(wrong).hexdigest()
    with pytest.raises(ValueError,match='archived upstream source'):
        driver['calculate'](saved=(forged,wrong))
