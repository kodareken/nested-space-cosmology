"""Finite response comparisons and a directly solved eliminated oscillator."""
import copy
import json
import math

import numpy as np
import pytest

from summarize_nsc_atlas_boson_response import (
    ONSET,READOUTS,_load_cache,aggregate,digest,memory_matrix,phase_partition,response_difference,
)


def record(degree=2,step=.02,value=1.,*,control=False,status='endpoint_reached'):
    responses=[dict(source_harmonic=source,**{name:source*value for name in READOUTS})
        for source in (1,2)]
    return dict(probe_degree=degree,max_step=step,decoupled_control=control,status=status,
        source_sha256={'src/model.py':'fixed'},background_record_sha256='record',
        background_history_sha256='history',background_producer_commit='producer',
        background_degree=16,background_epsilon=.25,
        background_temporal_tolerances=dict(rtol=3e-8,atol=3e-10,max_step=.05),
        checkpoints=[dict(coordinate_t=0.,observer_proper_duration=0.),
            dict(coordinate_t=4.,observer_proper_duration=4.36,response=responses,
                symplectic=dict(max_abs=1e-13))],
        onset_cubic_response=[dict(source_harmonic=source,**{name:source*value for name in ONSET})
            for source in (1,2)],
        max_checkpoint_symplectic_defect=2e-13)


def test_spatial_time_and_decoupled_comparisons_remain_separate():
    records=[record(6,value=1.1),record(2,value=2),record(4,value=1.2),
        record(4,.01,1.201),record(4,value=0,control=True)]
    result=aggregate(records,bindings=[True]*len(records))
    assert len(result['cases'])==5 and result['all_source_bindings_verified']
    assert len(result['spatial_changes'])==2 and len(result['temporal_changes'])==1
    assert len(result['decoupled_to_coupled_changes'])==1
    change=result['spatial_changes'][1]
    assert (change['from_probe_degree'],change['to_probe_degree'])==(4,6)
    assert change['responses'][0]['signed_change']['local_log_g_KK_per_unit_momentum_kick']==pytest.approx(-.1)
    assert change['onset_cubic_changes'][0]['signed_change']['volume_q_cubic_coefficient_norm']==pytest.approx(-.1)
    assert result['temporal_changes'][0]['responses'][1]['absolute_change']['scalar_reference_response_norm']==pytest.approx(.002)
    assert not result['physical_prediction_admitted'] and not result['entropy_or_noise_selected']
    assert result['cases'][0]['run'] is records[0]


def test_failed_and_unverified_records_are_preserved_without_comparing():
    records=[record(),record(4,status='failed'),record(6)]
    result=aggregate(records,bindings=[True,True,False])
    assert [case['disposition'] for case in result['cases']]==[
        'finite_endpoint_control','failed_attempt','unverified_source_binding']
    assert result['spatial_changes']==[] and not result['all_source_bindings_verified']


@pytest.mark.parametrize('changed',[
    'background_record_sha256','background_history_sha256','background_temporal_tolerances',
])
def test_distinct_backgrounds_cannot_be_merged(changed):
    a,b=record(),record(4)
    b[changed]='other'
    with pytest.raises(ValueError,match='share one background'):
        aggregate([a,b],bindings=[True,True])


def test_duplicates_sources_and_inconsistent_endpoints_are_rejected():
    a,b=record(),record(4)
    with pytest.raises(ValueError,match='duplicate'):
        aggregate([a,copy.deepcopy(a)],bindings=[True,True])
    b['source_sha256']={'src/model.py':'changed'}
    with pytest.raises(ValueError,match='producer source'):
        aggregate([a,b],bindings=[True,True])
    b=record(4);b['checkpoints'][-1]['observer_proper_duration']+=.01
    with pytest.raises(ValueError,match='physical endpoint'):
        aggregate([a,b],bindings=[True,True])


def test_undefined_and_zero_ratios_are_retained():
    a,b=record(value=0),record(4,value=1)
    a['checkpoints'][-1]['response'][1]['volume_to_scalar_reference_response_ratio']=None
    result=response_difference(a,b)
    assert result['responses'][1]['signed_change']['volume_to_scalar_reference_response_ratio'] is None
    assert result['responses'][0]['change_over_previous_absolute_value']['scalar_reference_response_norm'] is None


def oscillator_map(omega_s,omega_c,t):
    frequencies=np.array([omega_s,omega_c]);cos=np.diag(np.cos(frequencies*t))
    sine=np.diag(np.sin(frequencies*t))
    return np.block([[cos,np.diag(1/frequencies)@sine],[-np.diag(frequencies)@sine,cos]])


def generator(coupling):
    return np.block([[np.zeros((2,2)),np.eye(2)],
        [-np.array([[1.7**2,coupling],[coupling,2.3**2]]),np.zeros((2,2))]])


def test_schur_memory_uses_eliminated_volume_oscillator():
    c=.31;t=1.1;s=.27;D=generator(c)
    K=memory_matrix(D,D,oscillator_map(1.7,2.3,t),oscillator_map(1.7,2.3,s),np.array([0,1]))
    expected=np.zeros((2,2));expected[1,0]=c*c*math.sin(1.7*(t-s))/1.7
    np.testing.assert_allclose(K,expected,rtol=2e-15,atol=2e-17)
    zero=memory_matrix(generator(0),D,np.eye(4),np.eye(4),np.array([0,1]))
    np.testing.assert_array_equal(zero,np.zeros((2,2)))


def test_memory_uses_actual_inverse_without_symplectic_repair():
    fields=np.array([0,1]);_,y,_,_=phase_partition(fields)
    earlier=np.eye(4);later=np.eye(4)
    earlier[np.ix_(y,y)]=[[2.,.3],[0.,1.]]
    later[np.ix_(y,y)]=[[1.,.2],[-.1,1.]]
    c=.4;K=memory_matrix(generator(c),generator(c),later,earlier,fields)
    assert K[1,0]==pytest.approx(c*c*.05)
    earlier[np.ix_(y,y)]=0
    with pytest.raises(ValueError,match='singular'):
        memory_matrix(generator(c),generator(c),later,earlier,fields)


def test_cache_metadata_and_hash_binding(tmp_path):
    r=record();times=np.array([0.,4.]);fields=np.array([0,1]);path=tmp_path/'matrix.npz'
    metadata=dict(probe_degree=r['probe_degree'],history_sha256=r['background_history_sha256'],
        source_sha256=r['source_sha256'])
    np.savez_compressed(path,times=times,matrices=np.array([np.eye(4),np.eye(4)]),
        fields=fields,metadata=json.dumps(metadata))
    r['matrix_cache_sha256']=digest(path)
    actual,maps,labels=_load_cache(path,r)
    np.testing.assert_array_equal(actual,times);np.testing.assert_array_equal(labels,fields)
    assert maps.shape==(2,4,4)
    r['matrix_cache_sha256']='incorrect'
    with pytest.raises(ValueError,match='record hash'):_load_cache(path,r)
