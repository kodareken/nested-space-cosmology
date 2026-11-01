"""Wrapper unit checks with zero forcing/stub diagnostics, not a physical PDE run."""
import numpy as np
import pytest

from recursive_horizons import nsc_atlas_valve_spherical as candidate


def diagnostic(T,tau,y,yt,yx,yxx,h,H,kappa,sstar,*,ytt=None):
    c=y.shape[1]//2;em=np.exp(-y[1,c]);R=y[0,c]
    Rtt=yxx[0,c] if ytt is None else ytt[0,c]
    return dict(conformal_time=T,center_proper_time=tau,center_radius=R,center_S=y[2,c],
        center_Pi_S=em*yt[2,c],center_theta_plus_over_H=2*em*(yt[0,c]+yx[0,c])/(R*H),
        center_Rproper=em*yt[0,c],center_Rproper_acc=em**2*(Rtt-yt[1,c]*yt[0,c]-yx[1,c]*yx[0,c]),
        Hamiltonian_relative_RMS_term_budget=0.,Hamiltonian_relative_actual_metric_Riemann_scale=0.,
        metric_Einstein_residual_over_source_RMS=0.,central_half_Hamiltonian_RMS_over_H2=0.,
        min_radius=float(np.min(y[0])))


@pytest.fixture
def wrapper(monkeypatch):
    zero=lambda y,yt,yx,sstar,kappa:np.zeros_like(y)
    for name,value in [('conformal_forcing',zero),('conformal_diagnostics',diagnostic),('conformal_metric_diagnostics',diagnostic)]:
        monkeypatch.setattr(candidate,name,value)
    return candidate._conformal_evolve,candidate._conformal_evolve



def initial(asymmetric=False):
    x=np.linspace(-.2,.2,41);y=np.zeros((4,len(x)));v=np.zeros_like(y);dx=np.zeros_like(y);dxx=np.zeros_like(y)
    y[0]=2-.02*x*x;v[0]=-.1;dx[0]=-.04*x;dxx[0]=-.04
    y[2]=.005+.01*x*x;v[2]=.003;dx[2]=.02*x;dxx[2]=.02
    y[3]=.02*x;v[3]=.001*x;dx[3]=.02
    if asymmetric:
        y[0]+=.03*x;dx[0]+=.03;y[1]=.1*x;v[1]=.02+.01*x;dx[1]=.1
        y[3]+=.01;v[3]+=.005
    return dict(x=x,y=y,yt=v,yx=dx,yxx=dxx,step=float(x[1]-x[0]),half_width=.2)


def run(function,prepared,**kwargs):
    return function(prepared,H=.2,kappa=1.,sstar=0.,proper_time=.7,endpoint=.02,
        validity_mode='relative_metric_source',sample_every=.01,**kwargs)


def test_default_path_equals_explicit_reflection_and_exact_toy_update(wrapper):
    new,_=wrapper;a=run(new,initial());b=run(new,initial(),reflection_symmetry=True)
    np.testing.assert_array_equal(b['center_history'],a['center_history'])
    for key in ('y','yt','yx','yxx','ytt'):
        np.testing.assert_array_equal(b['end_profile'][key],a['end_profile'][key])
    assert b['events']==a['events'] and b['status']==a['status']
    assert b['reflection_symmetry'] is True
    p=a['end_profile'];x=p['x'];t=p['T']
    np.testing.assert_allclose(p['y'][0],2-.02*(x*x+t*t)-.1*t,rtol=0,atol=4e-15)
    np.testing.assert_allclose(p['y'][3],.02*x+.001*x*t,rtol=0,atol=3e-16)


def test_asymmetric_supplied_fields_and_jets_are_not_reflected(wrapper):
    new,_=wrapper;prepared=initial(True);before={k:v.copy() for k,v in prepared.items() if isinstance(v,np.ndarray)}
    result=run(new,prepared,reflection_symmetry=False);p=result['end_profile']
    # The zero-forcing toy has an exact linear chi solution.
    expected=.01+.02*p['x']+(.005+.001*p['x'])*p['T']
    np.testing.assert_allclose(p['y'][3],expected,rtol=0,atol=3e-16)
    assert abs(float(np.mean(p['y'][3])))>.009
    for key,value in before.items():np.testing.assert_array_equal(prepared[key],value)
    assert result['reflection_symmetry'] is False


def test_nonstationary_center_is_not_reported_as_a_neck_type_change(wrapper):
    new,_=wrapper;result=run(new,initial(True),reflection_symmetry=False)
    assert all(e['kind']!='spatial_neck_type_change_on_conformal_slices' for e in result['events'])
    assert 'constant X' in result['center_clock_role'] and 'accelerated' in result['center_clock_role']


def test_general_center_theta_and_hessian_roles_include_spatial_terms(wrapper):
    new,_=wrapper;result=run(new,initial(True),reflection_symmetry=False);p=result['end_profile'];c=len(p['x'])//2
    em=np.exp(-p['y'][1,c]);R=p['y'][0,c];history=result['center_history'][-1]
    expected_theta=2*em*(p['yt'][0,c]+p['yx'][0,c])/(R*.2)
    assert history[result['center_columns'].index('center_theta_plus_over_H')]==pytest.approx(expected_theta,abs=2e-15)
    proper_second=em**2*(p['ytt'][0,c]-p['yt'][1,c]*p['yt'][0,c])
    hessian=proper_second-em**2*p['yx'][1,c]*p['yx'][0,c]
    actual=history[result['center_columns'].index('center_Rproper_acc')]
    assert actual==pytest.approx(hessian,abs=2e-15)
    assert abs(actual-proper_second)>1e-4
    assert 'not d2R/dtau2' in result['center_Rproper_acc_role']


@pytest.mark.parametrize('bad',(None,0,1,'False',[]))
def test_reflection_flag_requires_boolean_before_initialization(wrapper,bad):
    new,_=wrapper
    with pytest.raises(ValueError,match='boolean reflection_symmetry'):
        run(new,{},reflection_symmetry=bad)
