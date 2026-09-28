"""Full-energy covariance and its non-closed probe/bulk interface."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from derive_nsc_pg_group13_covariance import assemble
from recursive_horizons.nsc_pg_state_covariance import spectral_probe_bulk_blocks,equal_time_ctp_blocks


def test_group13_full_energy_covariance_and_car():
    record=json.loads((ROOT/'results/development/nsc-pg-group13-covariance.json').read_text())
    with np.load(ROOT/record['payload']['path'],allow_pickle=False) as a:arrays={k:a[k].copy() for k in a.files}
    result,residuals=assemble(arrays);C=result['covariance'];G=result['car_gram']
    assert residuals['CAR']<3e-9
    assert np.linalg.eigvalsh(C).min()>0 and np.linalg.eigvalsh(C).max()<1
    assert np.linalg.norm(C[:4,4:])>1e-3
    assert np.linalg.norm(1j*(result['greater']-result['lesser'])-G)<3e-14
    assert np.linalg.norm(result['keldysh']-result['greater']-result['lesser'])<3e-14
    assert record['gate']['full_C1b']=='OPEN' and record['gate']['metric_timestep'] is False


def test_bulk_cross_terms_are_retained_in_field_reconstruction():
    # Numerical interface control, not a chosen NSC covariance completion.
    C=np.diag([.2,.8]);P=np.eye(2)
    F=np.array([[.4,.1j],[.2j,.3]],complex)
    J=np.array([[.5,0],[0,.2]],complex)
    phi=np.array([[.8,.3j],[-.2j,.9]],complex)
    blocks=spectral_probe_bulk_blocks(C,P,F,phi,J)
    restored=J@blocks['probe']@J.conj().T+J@blocks['probe_bulk']+blocks['probe_bulk'].conj().T@J.conj().T+blocks['bulk_diagonal']
    assert np.linalg.norm(restored-phi@C@phi.conj().T)<3e-14
    assert np.linalg.norm(blocks['probe_bulk'])>1e-3
    with pytest.raises(ValueError,match='projector'):
        spectral_probe_bulk_blocks(C,2*P,F,phi,J)


def test_ctp_relation_does_not_replace_a_computed_gram_by_identity():
    C=np.diag([.2,.7]);G=np.diag([1.,.999])
    result=equal_time_ctp_blocks(C,G)
    assert np.array_equal(1j*(result['greater']-result['lesser']),G)
