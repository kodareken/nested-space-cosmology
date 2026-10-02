"""Finite inherited grandchild: driven composition, clocks and locked prediction."""
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons import nsc_discovery_grandchild as grand
from recursive_horizons import nsc_discovery_prediction as prediction
from recursive_horizons import nsc_discovery_response as response
from recursive_horizons import nsc_nested_parent_child as model


def test_sequential_detail_requires_outer_drive_and_memory():
    rng=np.random.default_rng(901)
    size=20
    basis,_=np.linalg.qr(rng.normal(size=(size,8))+1j*rng.normal(size=(size,8)))
    H=rng.normal(size=(size,size))+1j*rng.normal(size=(size,size));H=H+H.conj().T
    field=rng.normal(size=(size,6))+1j*rng.normal(size=(size,6))
    direct,sequential=grand.initial_routes(field,basis)
    # Both exterior drive and memory are active; total source stays identical.
    outer=sequential['drive'].copy()
    sequential['drive']=.6*outer;sequential['memory']=.4*outer
    complement=direct['drive'].copy()
    direct['drive']=.6*complement;direct['memory']=.4*complement
    image=lambda values:H@values
    dr,sr=grand.route_rates(image,basis,direct,sequential)
    expected=-1j*H@field
    assert np.max(np.abs(grand.reconstruct_route(dr,basis)-expected))<2e-13
    assert np.max(np.abs(grand.reconstruct_route(sr,basis,True)-expected))<2e-13
    _,wrong=grand.route_rates(image,basis,direct,sequential,omit_outer_detail=True)
    assert np.linalg.norm(grand.reconstruct_route(wrong,basis,True)-expected)>1.
    exterior=(basis[:,4:8]@(sequential['f_drive']+sequential['f_memory'])+sequential['drive']+sequential['memory'])
    missing=-1j*basis[:,2:4].conj().T@H@exterior
    assert sr['d_memory']-wrong['d_memory']==pytest.approx(missing,abs=2e-13)


@pytest.fixture(scope='module')
def opening():
    return grand.load_opening(128)


def test_nested_bases_contain_original_observers_and_report_projection_tails(opening):
    pair=opening['pair'];basis,report=grand.observer_basis(pair)
    assert basis.shape==(2*128,8)
    assert report['ranks']==[2,4,8]
    assert report['single_opening_QR'] is True
    assert report['dilated_seed_projection_tail']>0.
    assert not basis.flags.writeable
    child=basis[:,:4];parent=basis[:,:8]
    assert np.max(np.abs(child@child.conj().T@basis[:,:2]-basis[:,:2]))<1e-13
    original_middle=pair.original_columns[:,[2,3]]
    assert np.max(np.abs(child@(child.conj().T@original_middle)-original_middle))<1e-13
    assert np.max(np.abs(parent@(parent.conj().T@pair.original_columns)-pair.original_columns))<1e-13
    assert report['spatial_tails']['grandchild']['exact_compact_support_claimed'] is False
    diagnostic=grand.endpoint_diagnostics(pair,opening['state'],basis)
    assert diagnostic['NG']!=pytest.approx(diagnostic['modal_trace_is_secondary'],abs=.01)
    assert diagnostic['grand_exterior_cross_frobenius']>1e-5
    assert diagnostic['source_CAR_admissible'] is True
    assert max(diagnostic['cross_force_max'].values())>1e-5


def test_opening_clock_alignment_recovers_missing_integral_and_keeps_delta_c(opening):
    clock=opening['clock_alignment']
    assert clock['delta_tau_primary']!=0.
    assert clock['delta_tau_primary']==pytest.approx(clock['delta_tau_secondary'],abs=1e-12)
    assert clock['new_segment_delta_tau']==0.
    assert opening['tangent'].occupations==pytest.approx([1.,1.,0.,0.,-1.,-1.])
    assert opening['held_tau_opening']-opening['tau_opening']!=0.
    # Physical content includes both field and occupation direction at fixed clock.
    raw=response.nested_readout(opening['pair'],opening['state'],opening['tangent'],coordinate_duration=0.)
    source=json.loads((grand.INPUT/'prediction-run.json').read_text())
    case=next(case for case in source['cases'] if case['nf']==128)
    assert raw['delta_child_regional_content_t']==pytest.approx(case['analytic_primary'],abs=1e-11)
    assert raw['delta_child_proper_mean_r_t']==pytest.approx(case['analytic_secondary'],abs=1e-11)


def test_routed_RK4_reconstructs_owned_step_without_source_or_geometry_reset(opening):
    pair=opening['held_pair'];state=opening['held'];basis,_=grand.observer_basis(opening['pair'])
    direct,sequential=grand.initial_routes(grand.columns(state),basis)
    token=grand.episode.state_sha256(state)
    stepped=grand.routed_step(pair,state,direct,sequential,basis,.0005)
    expected=prediction.nonlinear_rk4_step(pair,state,.0005)
    for name in model.STATE_NAMES:
        assert np.array_equal(getattr(stepped['state'],name),getattr(expected['state'],name))
    assert stepped['tau']==expected['tau']
    assert stepped['reconstruction_gap']<1e-12
    assert stepped['stage_source_gap']<1e-10
    assert grand.episode.state_sha256(state)==token


def test_locked_forecast_then_new_equal_tau_measurement_and_immutable_replay(tmp_path,monkeypatch):
    monkeypatch.setattr(model,'initial_state',lambda *args,**kwargs:pytest.fail('initializer used during inherited handoff'))
    destination=tmp_path/'short-grandchild'
    forecast=grand.prepare(destination,128,proper_increment=.001,matched=True)
    assert forecast['forecast_locked_before_measurement'] is True
    assert forecast['new_heldout_future_measured'] is False
    assert forecast['forecast_admission']['probes_use_baseline_only'] is True
    assert not (destination/'measurement.json').exists()
    before=(destination/'forecast.json').read_bytes()
    held_target=forecast['absolute_future_tau_target']-forecast['held_tau_opening']
    assert held_target!=forecast['proper_increment']
    measured=grand.run(destination)
    assert (destination/'forecast.json').read_bytes()==before
    assert measured['aggregate_CPU_seconds']<300.
    for row in measured['results']:
        assert row['status']=='heldout_future_measured'
        assert row['absolute_future_tau']==pytest.approx(forecast['absolute_future_tau_target'],abs=1e-9)
        assert row['event']['uses_linear_clock_correction'] is False
        assert row['measurement']['direct_sequential_Phi_gap']<1e-12
        assert row['measurement']['grand_covariance_composition_gap']<1e-12
        assert row['measurement']['direct_NG']==pytest.approx(row['measurement']['NG'],abs=1e-12)
        assert row['measurement']['sequential_NG']==pytest.approx(row['measurement']['NG'],abs=1e-12)
        assert row['measurement']['direct_clock_rate']==row['measurement']['sequential_clock_rate']
    assert measured['conditional_finite_handoff_only'] is True
    assert measured['universal_Omega_scaling_claimed'] is False
    assert grand.check(destination)['ok'] is True
    with pytest.raises(FileExistsError):grand.prepare(destination,128)
    with pytest.raises(FileExistsError):grand.run(destination)
    payload=destination/'measurement.npz'
    payload.chmod(0o644)
    payload.write_bytes(payload.read_bytes()+b'x')
    with pytest.raises(ValueError,match='payload binding'):grand.check(destination)
