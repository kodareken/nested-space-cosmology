"""Reference shell power, temporal gaps and matrix-record binding controls."""
import json
import math

import numpy as np
import pytest

import diagnose_nsc_atlas_boson_shells as diagnostic


def test_reference_shell_power_and_norm_fractions_are_distinct():
    sector=diagnostic.reconstruct_sector(2);n=len(sector['fields'])
    states=np.zeros((2,2*n));masses=np.array([2.,5.])
    low=np.flatnonzero((sector['fields']==0)&(sector['degrees']==0))[0]
    high=np.flatnonzero((sector['fields']==0)&(sector['degrees']==2))[0]
    chi=np.flatnonzero((sector['fields']==1)&(sector['degrees']==1))[0]
    states[0,n+low]=4.;states[0,n+high]=3.;states[1,chi]=1.
    rows=diagnostic.shell_diagnosis(states,sector,masses)
    first=rows[0]
    assert first['reference_norm']==5.
    assert first['last_degree_power_fraction']==pytest.approx(9/25)
    assert first['last_degree_norm_fraction']==pytest.approx(3/5)
    assert first['power_by_degree']['0']==16.
    assert rows[1]['reference_norm']==0 and rows[1]['last_degree_power_fraction'] is None
    assert rows[3]['reference_norm']==pytest.approx(math.sqrt(5.))
    assert rows[3]['last_degree_power_fraction']==0.
    assert rows[3]['last_two_degrees_power_fraction']==pytest.approx(1.)
    np.testing.assert_allclose(np.linalg.norm(sector['kicks'],axis=1),1.,rtol=0,atol=4e-16)
    assert sector['kicks'][0]@sector['kicks'][1]==pytest.approx(0.,abs=2e-15)


def test_whole_state_gap_uses_the_tighter_state_as_denominator():
    sector=diagnostic.reconstruct_sector(2);states=sector['kicks'].copy()
    rows=diagnostic.relative_state_gaps(states,2*states,sector,np.array([2.,5.]))
    assert rows[0]['gap_over_tighter_reference_norm'] is None
    assert rows[1]['absolute_reference_norm_gap']==pytest.approx(1.)
    assert rows[1]['gap_over_tighter_reference_norm']==pytest.approx(.5)
    assert rows[3]['gap_over_tighter_reference_norm']==pytest.approx(.5)


def make_case(step,*,control=False,background='same'):
    sector=diagnostic.reconstruct_sector(2)
    record=dict(probe_degree=2,max_step=step,decoupled_control=control,
        background_record_sha256=background,background_history_sha256='history',
        background_producer_commit='producer',background_degree=16,background_epsilon=.25,
        background_temporal_tolerances={'rtol':3e-8},source_sha256={'src/model.py':'fixed'})
    return dict(record=record,sector=sector,states=sector['kicks'].copy(),
        summary=dict(coordinate_endpoint=4.,observer_proper_endpoint=4.36))


def test_only_like_time_controls_are_compared_and_duplicates_are_rejected():
    action=dict(sstar=0.,masses_squared=np.array([2.,5.]),source_sha256={})
    result=diagnostic.aggregate([make_case(.02),make_case(.01),make_case(.02,control=True)],action)
    assert len(result['temporal_changes'])==1
    assert result['temporal_changes'][0]['from_max_step']==.02
    assert result['temporal_changes'][0]['to_max_step']==.01
    assert not result['physical_prediction_admitted'] and not result['spatial_convergence_decided']
    with pytest.raises(ValueError,match='duplicate'):
        diagnostic.aggregate([make_case(.02),make_case(.02)],action)
    with pytest.raises(ValueError,match='share one background'):
        diagnostic.aggregate([make_case(.02),make_case(.01,background='different')],action)


def test_cache_sha_metadata_field_order_and_source_binding_are_checked(tmp_path,monkeypatch):
    monkeypatch.setattr(diagnostic,'ROOT',tmp_path)
    monkeypatch.setattr(diagnostic,'source_binding',lambda record:record.get('bound',False))
    sector=diagnostic.reconstruct_sector(2);n=len(sector['fields'])
    record_path=tmp_path/'record.json';matrix_path=tmp_path/'matrix.npz'
    action=dict(sstar=0.,masses_squared=np.array([2.,5.]),source_sha256={'src/model.py':'fixed'})
    record=dict(bound=True,status='endpoint_reached',probe_degree=2,max_step=.02,
        decoupled_control=False,source_sha256=action['source_sha256'],
        background_history_sha256='history',checkpoints=[
            dict(coordinate_t=0.,observer_proper_duration=0.),
            dict(coordinate_t=1.,observer_proper_duration=.9)])
    metadata=dict(probe_degree=2,history_sha256='history',source_sha256=record['source_sha256'])
    def write_cache(*,field_order=None,cache_metadata=None):
        np.savez_compressed(matrix_path,times=np.array([0.,1.]),
            matrices=np.array([np.eye(2*n),np.eye(2*n)]),
            fields=sector['fields'] if field_order is None else field_order,
            metadata=json.dumps(metadata if cache_metadata is None else cache_metadata))
        record['matrix_cache_sha256']=diagnostic.digest(matrix_path)
        record_path.write_text(json.dumps(record),encoding='utf-8')
    write_cache()
    case=diagnostic.load_verified_case(record_path,matrix_path,action)
    assert case['summary']['matrix_binding_verified'] and case['summary']['source_binding_verified']
    record['matrix_cache_sha256']='wrong';record_path.write_text(json.dumps(record),encoding='utf-8')
    with pytest.raises(ValueError,match='SHA'):diagnostic.load_verified_case(record_path,matrix_path,action)
    write_cache(cache_metadata=dict(metadata,history_sha256='other'))
    with pytest.raises(ValueError,match='metadata'):diagnostic.load_verified_case(record_path,matrix_path,action)
    write_cache(field_order=1-sector['fields'])
    with pytest.raises(ValueError,match='field order'):diagnostic.load_verified_case(record_path,matrix_path,action)
    write_cache();record['bound']=False;record_path.write_text(json.dumps(record),encoding='utf-8')
    with pytest.raises(ValueError,match='source-bound'):diagnostic.load_verified_case(record_path,matrix_path,action)
    record['bound']=True;record_path.write_text(json.dumps(record),encoding='utf-8')
    different=dict(action,source_sha256={'src/model.py':'changed'})
    with pytest.raises(ValueError,match='differs from producer'):
        diagnostic.load_verified_case(record_path,matrix_path,different)
