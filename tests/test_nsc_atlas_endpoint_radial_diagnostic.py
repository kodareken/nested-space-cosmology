"""A one-stage radial diagnosis authenticates its exact local state."""
import hashlib
import numpy as np
import pytest

from diagnose_nsc_atlas_endpoint_radial import read_seed, panel_step_options, inward_panel_point


def fixture(root):
    path=root/'.build/results/stage.npz'
    path.parent.mkdir(parents=True)
    np.savez_compressed(path,x=np.array(.43),level=np.array(0),
        source_commit=np.array('a'*40),**{'cap.chi':np.array([.1,.2]),'wall.r':np.array(.5)})
    raw=path.read_bytes()
    return path,dict(level=0,source_commit='a'*40,failed_stage_seed=dict(
        path=path.relative_to(root).as_posix(),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),
        source_x=.43))


def test_authenticated_seed_recovers_exact_arrays(tmp_path):
    _,receipt=fixture(tmp_path)
    state,x=read_seed(receipt,tmp_path)
    assert x==.43 and float(state['wall']['r'])==.5
    np.testing.assert_array_equal(state['cap']['chi'],[.1,.2])


def test_changed_seed_or_metadata_is_rejected(tmp_path):
    path,receipt=fixture(tmp_path)
    receipt['level']=1
    with pytest.raises(ValueError,match='metadata'):read_seed(receipt,tmp_path)
    receipt['level']=0
    path.write_bytes(path.read_bytes()+b'changed')
    with pytest.raises(ValueError,match='bytes changed'):read_seed(receipt,tmp_path)
    receipt['failed_stage_seed']['path']='outside.npz'
    with pytest.raises(ValueError,match='ignored result directory'):read_seed(receipt,tmp_path)


def test_panel_step_check_preserves_tolerance_and_existing_smaller_step():
    options=dict(rtol=2.5e-14,atol=1e-15,first_step=.4)
    assert panel_step_options((1.,1.5),options,None)==options
    assert panel_step_options((1.,1.5),options,.5)==dict(
        rtol=2.5e-14,atol=1e-15,first_step=.25,max_step=.25)
    assert options['first_step']==.4
    assert panel_step_options((1.,1.5),dict(max_step=.1),.5)['max_step']==.1
    for fraction in (0,-1,1.1,np.inf,np.nan):
        with pytest.raises(ValueError,match='fraction'):
            panel_step_options((1.,1.5),options,fraction)
    with pytest.raises(ValueError,match='positive radial panel'):
        panel_step_options((1.,1.),options,.5)


def test_one_sided_diagnostic_moves_only_resolved_endpoint_samples():
    span=(1.,1.25)
    assert span[0]<inward_panel_point(span[0],span)<1.000000000001
    assert 1.249999999999<inward_panel_point(span[1],span)<span[1]
    assert inward_panel_point(1.125,span)==1.125
    with pytest.raises(ValueError,match='not resolved'):
        inward_panel_point(1.,(1.,np.nextafter(1.,2.)))
