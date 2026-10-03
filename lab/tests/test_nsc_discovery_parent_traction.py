"""Reviewed local equations, independent linear checks and immutable evidence order."""
import json
import numpy as np
import pytest
import sympy as sp
from scipy.integrate import quad
from recursive_horizons import nsc_discovery_parent_traction as collar


@pytest.fixture(scope='module')
def constants():return collar.load_inputs()[0]


@pytest.fixture
def frozen_writes(monkeypatch):
    # Unit records use a simulated immutable commit, never mutate real Git.
    monkeypatch.setattr(collar,'producing_commit',lambda pins:'a'*40)
    monkeypatch.setattr(collar.provenance,'resolve_pinned_source_bytes',lambda *args,**kwargs:b'unit frozen source')


def test_symbolic_first_integral_is_constant_for_actual_static_equations():
    r,s,Q,hx,u,v,g,mag,M,w,k,eps=sp.symbols('r s Q hx u v g mag M w k eps',nonzero=True)
    n=u*u+v*v;S=2*u*v;K=eps*n-k*Q*S
    I=Q*Q*(r*r-mag/g)-3*s*s-2*r*s*hx+M*w*K/g
    rxx=Q*Q*(r*r-mag/g)/r-s*s/r-M*w*k*Q*S/(2*g*r)
    hxx=Q*Q-3*rxx/r
    directions=(s,rxx,Q*hx,hxx,eps*v-k*Q*u,k*Q*v-eps*u)
    assert sp.factor(sum(sp.diff(I,key)*rate for key,rate in zip((r,s,Q,hx,u,v),directions)))==0
    B=g*(r*s+r*r*hx)
    derivative=sum(sp.diff(B,key)*rate for key,rate in zip((r,s,Q,hx,u,v),directions))
    assert sp.factor(derivative-mag*Q*Q-M*w*eps*n+g*I)==0


def test_exact_BR_profiles_static_RHS_and_linear_Green_function(constants):
    k=constants;x=.7;Q=1/np.cos(x);u,v=collar.baseline_spinor(x)
    state=np.array([1.,0.,np.log(Q),np.tan(x),u,v,0.,0.,0.])
    image=collar.rhs(x,state,0.,k)
    expected_u=u*(-np.tan(x)+np.cos(x)/(2*(1+np.sin(x))))
    expected_v=v*(-np.tan(x)-np.cos(x)/(2*(1-np.sin(x))))
    np.testing.assert_allclose(image[:6],[0.,0.,np.tan(x),Q*Q,expected_u,expected_v],atol=1e-14)
    assert u*u+v*v==pytest.approx(2*np.cos(x)**2/np.pi,abs=1e-15)
    report,arrays=collar.exact_local(0.,k)
    np.testing.assert_allclose(arrays['state'][:,0],1.,atol=1e-11)
    np.testing.assert_allclose(np.exp(arrays['state'][:,2]),1/np.cos(arrays['x']),atol=2e-10)
    up,vp=collar.baseline_spinor(arrays['x'])
    np.testing.assert_allclose(arrays['state'][:,4],up,atol=2e-10)
    np.testing.assert_allclose(arrays['state'][:,5],vp,atol=2e-10)
    linear=collar.linear_solution(k)
    for cut in (.5,1.):
        index=np.argmin(abs(linear['x']-cut))
        green=quad(lambda s:(np.tan(cut)*(1+s*np.tan(s))-(1+cut*np.tan(cut))*np.tan(s))
                          *(-3*collar.radius_response(s,k)[2]),0.,cut,epsabs=2e-13,epsrel=2e-13)[0]
        assert linear['q'][index]==pytest.approx(green,abs=2e-11)
    assert linear['a'][-1]==pytest.approx(-2.07627439614124,abs=1e-13)
    assert linear['q'][-1]==pytest.approx(5.66623962782446,abs=2e-11)


def test_other_weights_quadratic_remainder_admissibility_and_quadrature(constants):
    k=constants;linear=collar.linear_solution(k);remainders=[]
    for w in (.0003,.0006):
        measured,arrays=collar.exact_local(w,k)
        predicted=collar.linear_readouts(linear,k,w)
        assert abs(collar.first_integral(collar.center_state(w,k),w,k))<2e-15
        assert measured['first_integral_max']<1e-10 and measured['r_min']>0 and measured['Q_min']>0
        assert measured['center_proper_lapse_N']==measured['center_radius']<1
        for name in collar.CUTS:
            value=measured['readouts'][name]
            assert abs(value['static_boundary_identity_residual'])<1e-9
            assert value['covariance_admissible_on_this_collar']
            assert value['finite_collar_covariance_eigenvalue']==pytest.approx(w*value['finite_collar_mode_norm'])
            assert value['partial_source_energy']==pytest.approx(k['M']*k['epsilon']*value['finite_collar_covariance_eigenvalue'])
            assert abs(predicted[name]['linear_static_boundary_identity_residual'])<1e-11
        remainders.append({key:measured['readouts']['parent'][key]-predicted['parent'][key] for key in ('radius','Q','full_collar_boundary_charge')})
    for key in remainders[0]:assert 3.8<remainders[1][key]/remainders[0][key]<4.2
    coarse,_=collar.exact_local(.0006,k,points=129)
    fine,_=collar.exact_local(.0006,k,points=257,rtol=2e-13,atol=1e-14)
    for key in ('mode_norm_quadrature_gap','magnetic_quadrature_gap'):
        assert abs(fine['readouts']['parent'][key])<abs(coarse['readouts']['parent'][key])/10
    with pytest.raises(ValueError):collar.center_state(-.1,k)
    with pytest.raises(ValueError):collar.center_state(100.,k)


def test_prediction_lock_creation_only_replay_and_input_preservation(tmp_path,monkeypatch,frozen_writes):
    _,hashes,_=collar.load_inputs();original={name:collar.sha(name) for name in hashes}
    with pytest.raises(FileNotFoundError):collar.measure(tmp_path)
    prepared=collar.prepare(tmp_path);locked=collar.predict(tmp_path)
    assert prepared['held_w']==.0013 and not prepared['nonlinear_held_generated']
    assert locked['locked_before_held_nonlinear'] and locked['mode']=='locked_linear_prediction'
    assert collar.check(tmp_path)['prediction_replayed']
    with pytest.raises(FileExistsError):collar.prepare(tmp_path)
    with pytest.raises(FileExistsError):collar.predict(tmp_path)
    def verify_lock_then_stop(w,k,**kwargs):
        record=json.loads((tmp_path/'prediction.json').read_text())
        assert record['locked_before_held_nonlinear'] and w==record['held_w']
        assert kwargs['cpu_remaining']<=30
        raise RuntimeError('test stops before held nonlinear solve')
    monkeypatch.setattr(collar,'exact_local',verify_lock_then_stop)
    with pytest.raises(RuntimeError,match='test stops'):collar.measure(tmp_path)
    assert not (tmp_path/'measurement.json').exists()
    assert original=={name:collar.sha(name) for name in hashes}


def test_unfrozen_record_creation_is_rejected_without_writes(tmp_path,monkeypatch):
    monkeypatch.setattr(collar,'producing_commit',lambda pins:None)
    with pytest.raises(ValueError,match='frozen non-None'):collar.prepare(tmp_path)
    assert list(tmp_path.iterdir())==[]


def test_historical_source_change_authenticates_without_current_math_replay(tmp_path,monkeypatch,frozen_writes):
    collar.prepare(tmp_path);collar.predict(tmp_path)
    original=collar.source_hashes();changed=dict(original);changed[next(iter(changed))]='0'*64
    monkeypatch.setattr(collar,'source_hashes',lambda:changed)
    monkeypatch.setattr(collar,'linear_readouts',lambda *args:pytest.fail('historical data replayed with current math'))
    checked=collar.check(tmp_path)
    assert checked['historical_sources_authenticated'] and checked['mode']=='historical_authentication_only'
    assert not checked['numerical_replay'] and not checked['prediction_replayed']


def edit_json(path,update):
    record=json.loads(path.read_text());update(record);path.chmod(0o600)
    path.write_text(json.dumps(record,sort_keys=True))


def test_measurement_rejects_preparation_mutation_before_solver(tmp_path,monkeypatch,frozen_writes):
    first=tmp_path/'preparation_changed';collar.prepare(first);collar.predict(first)
    edit_json(first/'prepare.json',lambda record:record.update(test_mutation=True))
    monkeypatch.setattr(collar,'exact_local',lambda *args,**kwargs:pytest.fail('solve started before preparation authentication'))
    with pytest.raises(ValueError,match='preparation binding failed'):collar.measure(first)
    assert not (first/'measurement.json').exists()


def test_prediction_mutation_during_measurement_is_rejected(tmp_path,monkeypatch,frozen_writes,constants):
    collar.prepare(tmp_path);collar.predict(tmp_path)
    cached,arrays=collar.exact_local(.0003,constants)
    def mutate_then_return_other_weight(*args,**kwargs):
        edit_json(tmp_path/'prediction.json',lambda record:record.update(test_mutation=True))
        return cached,{key:value.copy() for key,value in arrays.items()}
    monkeypatch.setattr(collar,'exact_local',mutate_then_return_other_weight)
    with pytest.raises(ValueError,match='prediction changed during collar measurement'):collar.measure(tmp_path)
    assert not (tmp_path/'measurement.json').exists() and not (tmp_path/'measurement.npz').exists()


def test_missing_cut_grid_rejected_before_solver(constants,monkeypatch):
    monkeypatch.setattr(collar,'solve_ivp',lambda *args,**kwargs:pytest.fail('invalid grid reached solver'))
    with pytest.raises(ValueError,match='every declared cut'):collar.exact_local(.0003,constants,points=256)


def test_small_weight_measurement_replays_primary_gradients_and_center_clock(tmp_path,monkeypatch,frozen_writes):
    monkeypatch.setattr(collar,'HELD_WEIGHT',.0003)
    collar.prepare(tmp_path);locked=collar.predict(tmp_path);collar.measure(tmp_path)
    assert collar.check(tmp_path)['measurement_replayed']
    def coordinated_readout_edit(record):
        record['measurement']['readouts']['parent']['boundary_gradient_r_x']+=1e-4
        record['comparisons']=collar.comparisons(locked,record['measurement'])
    edit_json(tmp_path/'measurement.json',coordinated_readout_edit)
    with pytest.raises(ValueError,match='endpoint/source/clock replay differs'):collar.check(tmp_path)


def test_direct_clock_and_charge_replay_rejects_misreported_small_weight(constants):
    import copy
    report,arrays=collar.exact_local(.0006,constants)
    for key in ('boundary_gradient_h_x','boundary_gradient_Q_x','right_endpoint_charge','left_endpoint_charge','static_boundary_identity_residual'):
        mutated=copy.deepcopy(report);mutated['readouts']['child'][key]+=1e-4
        with pytest.raises(ValueError,match='endpoint/source/clock'):collar._verify_nonlinear(mutated,arrays,constants,.0006)
    mutated=copy.deepcopy(report);mutated['center_proper_frequency']+=1e-4
    with pytest.raises(ValueError,match='center/clock/conservation'):collar._verify_nonlinear(mutated,arrays,constants,.0006)
