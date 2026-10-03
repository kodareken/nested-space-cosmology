"""Independent equation, source/tangent/event and immutable lock checks."""
import json
from types import SimpleNamespace
import numpy as np
import pytest
from threadpoolctl import threadpool_limits
from recursive_horizons import nsc_discovery_homogeneous_turn as turn


@pytest.fixture(scope='module')
def data():
    return turn.load_preparation()


def test_symmetry_reduction_matches_full_leading_RHS_and_energy(data):
    k=data['constants'];y=data['initial_state'].copy()
    y[:4]=[.35,1.2,.1,-.2]
    phi=turn.reconstruct_columns(y,data['frame_coefficients'],k,32)
    le=turn.leading
    with threadpool_limits(limits=1),le.backend.fft_thread_limit(1):
        pair=le.nested.build_pair(32,coarse_modes=1,child_details=2,source_layout='override',
                                  columns_override=phi,occupations=np.repeat(turn.BASE_C,2))
        pair=le.backend.make_fft_pair(pair)
        state=le.encode(pair,le.State(*(np.full(pair.grid.ng,value) for value in y[:4]),*phi))
        rate=le.decode(pair,le.rates(pair,state));reduced=turn.rhs(0,y,turn.BASE_C,k)
        for index,name in enumerate(('Q','r','p_Q','p_r')):
            np.testing.assert_allclose(getattr(rate,name),reduced[index],rtol=1e-11,atol=1e-11)
        expected=turn.reconstruct_columns(reduced,data['frame_coefficients'],k,32)
        np.testing.assert_allclose(rate.phi0,expected[0],atol=1e-12)
        np.testing.assert_allclose(rate.phi1,expected[1],atol=1e-12)
        energy=le.energy(pair,state);values=turn.quantities(y,turn.BASE_C,k)
        assert energy['field']==pytest.approx(values['field_energy'],abs=1e-11)
        assert energy['gravity']==pytest.approx(values['geometry_energy'],abs=1e-11)


def test_analytic_state_and_source_Jv_and_initial_constraint_derivative(data):
    k=data['constants'];y=data['initial_state'].copy();y[:4]=[.3,1.4,.2,-.1]
    delta=np.random.default_rng(25).normal(size=16)*.03;dc=turn.DIRECTION
    image=turn.jvp(y,delta,turn.BASE_C,dc,k);h=1e-6
    numeric=(turn.rhs(0,y+h*delta,turn.BASE_C+h*dc,k)-turn.rhs(0,y-h*delta,turn.BASE_C-h*dc,k))/(2*h)
    np.testing.assert_allclose(image,numeric,rtol=1e-7,atol=2e-9)
    initial=turn.initial_tangent(turn.BASE_C,dc,k,data['Q0'])
    numeric=(turn.initial_state(turn.BASE_C+h*dc,k,data['Q0'])-turn.initial_state(turn.BASE_C-h*dc,k,data['Q0']))/(2*h)
    np.testing.assert_allclose(initial,numeric,atol=1e-9)
    assert np.count_nonzero(initial)==1 and initial[0]==0 and np.all(initial[4:]==0)


def test_autonomous_event_constraint_energy_transfer_and_reference(data):
    report,arrays=turn.integrate(turn.BASE_C,data['constants'],data['Q0'],direction=turn.DIRECTION)
    assert 0<report['field_energy_zero_time']<report['time']<3
    assert report['p_Q_dot']<0 and report['r_ddot']>0
    assert report['event_energy']['field_energy']<0 and report['submagnetic_turn_requires_negative_field_energy']
    assert abs(report['turn_constraint_radius_squared_gap'])<1e-8
    assert report['r_ddot']==pytest.approx(report['turn_acceleration_from_energy'],abs=1e-9)
    assert report['maximum_total_energy_residual']<1e-8
    assert abs(report['fieldwork_closure_gap'])<1e-8 and abs(report['geometry_transfer_closure_gap'])<1e-8
    assert report['source_trace']==3 and report['spinor_norm_gap']<1e-10
    assert report['Q_positive'] and report['r_positive']
    assert report['final_energy']['field_energy']==pytest.approx(data['fullband_reference']['energy']['field'],abs=1e-8)
    assert report['radius']<data['fullband_reference']['sampled_minimum_radius']
    assert arrays['event_tangent'].shape==(16,)
    with pytest.raises(ValueError):turn.initial_state(np.array([1.01,.5,.25]),data['constants'],data['Q0'])


def test_creation_only_lock_replay_before_independent_measurement(tmp_path,monkeypatch,data):
    assert str(turn.INITIAL)+'.json' in data['input_hashes']
    assert str(turn.REFERENCE)+'.npz' in data['input_hashes']
    before={path:turn.sha(path) for path in data['input_hashes']}
    with pytest.raises(FileNotFoundError):turn.measure(tmp_path)
    turn.prepare(tmp_path);locked=turn.predict(tmp_path)
    assert locked['preparation']['leading_initial_record']==str(turn.INITIAL)+'.json'
    assert locked['locked_before_held_data'] and not locked['held_data_generated']
    assert locked['preparation']['held_pair_occupations']==[.752,.752,.5,.5,.248,.248]
    assert max(abs(row['analytic_difference']) for row in locked['centred_event_checks'])<1e-6
    assert turn.check(tmp_path)['prediction_replayed']
    with pytest.raises(FileExistsError):turn.prepare(tmp_path)
    with pytest.raises(FileExistsError):turn.predict(tmp_path)
    def independent(path,*,cpu_limit_seconds):
        observed=json.loads(path.read_text())
        assert observed['locked_before_held_data'] and observed['mode']=='locked_prediction'
        assert cpu_limit_seconds<=30
        return {'events':[{'source':'held','radius':observed['predicted_held_radius']+1e-6,'time':2.6}]},{'tiny_endpoint':np.array([1.])}
    monkeypatch.setattr(turn.importlib,'import_module',lambda name:SimpleNamespace(measure_locked_prediction=independent))
    measured=turn.measure(tmp_path)
    assert measured['comparisons'][0]['prediction_error']==pytest.approx(1e-6)
    assert turn.check(tmp_path)['measurement_authenticated']
    assert before=={path:turn.sha(path) for path in before}
