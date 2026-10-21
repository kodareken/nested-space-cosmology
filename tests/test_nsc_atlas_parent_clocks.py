import math

import numpy as np
import pytest
from scipy.special import eval_gegenbauer

from recursive_horizons.nsc_atlas_parent_clocks import ParentClocks,infer_two_amplitudes


@pytest.fixture(scope='module')
def clocks():
    return ParentClocks()


@pytest.fixture(scope='module')
def responses(clocks):
    return [clocks.response(x) for x in (.04,.08,.12)]


def test_spatial_derivative_has_an_independent_polynomial_check(clocks):
    for theta in (0.,1e-7,.04,.12):
        Q,Qr=clocks.harmonics(theta)
        expected=-2*math.sin(theta)*eval_gegenbauer(clocks.modes-1,2,math.cos(theta))/(clocks.modes+1)
        np.testing.assert_allclose(Qr,expected,atol=3e-11,rtol=2e-10)
        np.testing.assert_allclose(Q,clocks.packet.radial_harmonics(clocks.modes,[theta])[:,0],atol=2e-12)


def test_clock_protocol_is_causal_and_uses_initial_proper_distance(clocks,responses):
    assert all(0<r['background_emission_time']<clocks.to for r in responses)
    assert all(0<r['background_frequency_ratio']<1 for r in responses)
    assert all(r['initial_coordinate_displacement'][0]<0 for r in responses)
    assert all(r['initial_coordinate_displacement'][1]==0 for r in responses)
    assert max(r['initial_angle'] for r in responses)<clocks.packet.support_angle
    for r in responses:
        assert r['initial_proper_separation']==pytest.approx(clocks.initial['a']*r['initial_angle'])
        np.testing.assert_allclose(sum(map(np.array,r['parts'].values())),r['log_frequency_response'],atol=1e-15)


def test_receiver_proper_age_cannot_be_replaced_by_coordinate_time(clocks,responses):
    r=responses[-1];be=clocks.background(r['background_emission_time']);bo=clocks.background(clocks.to)
    correction=np.array(r['receiver_time_shift'])*(be['H']*be['a']/bo['a']-bo['H'])
    assert abs(correction[0])>.01
    assert abs(correction[1])>.0005


@pytest.mark.parametrize('index',[0,1,2])
def test_finite_metric_geodesics_reproduce_the_linear_observable(clocks,responses,index):
    row=responses[index];position=row['initial_angle']
    quiet=clocks.finite_metric_clock(position,[0,0])
    assert quiet['log_frequency_ratio']==pytest.approx(math.log(row['background_frequency_ratio']),abs=2e-11)
    direction=np.array([1.,-.4]);target=np.array(row['log_frequency_response'])@direction
    errors=[]
    for h in (.001,.0005):
        p=clocks.finite_metric_clock(position,h*direction)
        m=clocks.finite_metric_clock(position,-h*direction)
        derivative=(p['log_frequency_ratio']-m['log_frequency_ratio'])/(2*h)
        errors.append(abs(derivative-target))
    assert max(errors)<2e-9
    assert errors[1]<errors[0]/2


def test_two_clocks_recover_state_and_predict_a_withheld_clock(responses):
    R=np.array([r['log_frequency_response'] for r in responses]);amplitudes=np.array([1e-5,-4e-6])
    y=R@amplitudes
    recovered=infer_two_amplitudes(responses,y[:2])
    np.testing.assert_allclose(recovered,amplitudes,atol=2e-18)
    assert R[2]@recovered==pytest.approx(y[2],abs=2e-21)
    weights=np.r_[-R[2]@np.linalg.inv(R[:2]),1.]
    np.testing.assert_allclose(weights@R,0,atol=1e-16)
    corrupted=y.copy();corrupted[2]+=2.5e-9
    assert (weights@corrupted)/(np.linalg.norm(weights)*1e-10)>5


def test_spatial_and_time_accuracy_preserve_the_clock_relation(responses):
    other=ParentClocks(cutoff=256,rtol=2e-11,atol=2e-13)
    R=np.array([r['log_frequency_response'] for r in responses])
    R2=np.array([other.response(x)['log_frequency_response'] for x in (.04,.08,.12)])
    assert np.max(abs(R2-R))<1e-9
    w=-R[2]@np.linalg.inv(R[:2]);w2=-R2[2]@np.linalg.inv(R2[:2])
    assert np.max(abs(w2-w))<1e-7


def test_out_of_domain_clocks_and_unresolved_inverses_are_rejected(clocks,responses):
    with pytest.raises(ValueError):clocks.response(.5)
    with pytest.raises(ValueError):ParentClocks(observation_efold=1)
    with pytest.raises(ValueError):clocks.finite_metric_clock(.12,[1.,0.])
    with pytest.raises(ValueError):infer_two_amplitudes([responses[0],responses[0]],[0.,0.])
