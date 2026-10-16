"""Small manufactured caps; no retained state, source solve or trajectory."""
import copy

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_background_cap import background_cap_metric
from recursive_horizons.nsc_atlas_cap_mass_tangent import mass_tangent
from recursive_horizons.nsc_atlas_cap_transport_fields import CapTransportFields
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile, SMALL_RADIUS
from recursive_horizons.nsc_atlas_cap_radial_comparison import (
    compare_cap_radial_maps, _reference_pieces, _reference_map, _native_map,
)

CONTROLS=ControlBudget(maximum_evaluations=1000000,cpu_seconds=40.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('cap radial composition controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def setup(budget,*,radius=.4,alpha=.002,cubic=False,state=None,rates=None,order=8):
    knots=np.array([-1.,0.,1.]); c=float(cubic)
    profile=CharacteristicProfile(knots,c*knots**3,3*c*knots**2,6*c*knots)
    state=state or dict(Ip=np.array([.002,.005,.013,.025]),Iq=np.array([.003,.001,-.0005,-.002]))
    rates=rates or dict(R=.02,T=.3,Ip=np.array([.01,.03,.02,.04]),Iq=np.array([-.01,.02,.04,.01]))
    metric=background_cap_metric(np.linspace(0,1,5),radius,state,T=.1,profile=profile,
                                 alpha=alpha,quadrature_order=order,budget=budget)
    tangent=mass_tangent(metric,rates,budget=budget)
    return metric,rates,tangent


def test_exact_zero_coupling_constant_field_mass_and_full_moving_rate(budget):
    R,p,q,rx,pdot,qdot=.4,.3,.2,.07,.11,-.04
    faces=R*np.linspace(0,1,5)
    state=dict(Ip=p*np.diff(faces**2)/2,Iq=q*np.diff(faces))
    rates=dict(R=rx,T=.2,Ip=(pdot/p+2*rx/R)*state['Ip'],Iq=(qdot/q+rx/R)*state['Iq'])
    metric,rates,tangent=setup(budget,alpha=0.,state=state,rates=rates)
    result=compare_cap_radial_maps(metric,rates,tangent,budget=budget)
    mass=q*q*R**3/6+p*p*(R**3/12-R**5/20)
    rho=R*R*q*q/2+R*R*(1-R*R)*p*p/4
    rate=q*qdot*R**3/3+p*pdot*(R**3/6-R**5/10)+rx*rho
    for key in ('native_mass','reference_mass'):
        np.testing.assert_allclose(result[key],mass,rtol=1e-12,atol=2e-15)
    for key in ('native_rate','reference_rate'):
        np.testing.assert_allclose(result[key],rate,rtol=1e-12,atol=2e-15)
    assert np.all(abs(result['telescope_defects'])<=result['telescope_roundoff_bounds'])
    assert not result['uniform_error_enclosure']


def independent_variation(metric,rates,panels,budget):
    view=CapTransportFields(metric,rates,budget=budget)
    n=material=0.
    for panel in panels:
        for piece in _reference_pieces(panel):
            def coefficients(r):
                if r==0: return view.reference_axis_coefficients(budget=budget)
                row=view.reference_sample(piece['owner'],[r],budget=budget,
                    small_radius_side=piece['small_radius_side'] if r==SMALL_RADIUS else None)
                return {k:float(row[k][0]) for k in ('a','b','a_x','b_x')}
            low,high=piece['left'],piece['right']
            begin=coefficients(low); rho=begin['b']-begin['a']*n
            if low==high:
                material+=(piece['right_rate']-piece['left_rate'])*rho
                continue
            nu=material-piece['left_rate']*rho
            def rhs(r,state):
                budget.check(8); f=coefficients(r); value,partial=state
                return [f['b']-f['a']*value,f['b_x']-f['a_x']*value-f['a']*partial]
            flow=solve_ivp(rhs,(low,high),[n,nu],method='DOP853',rtol=2e-12,atol=2e-14,
                           max_step=(high-low)/4)
            assert flow.success
            n,nu=flow.y[:,-1]; end=coefficients(high)
            material=nu+piece['right_rate']*(end['b']-end['a']*n)
    return n,material


def test_nonconstant_cap_matches_native_tangent_and_independent_reference_variation(budget):
    # Order four leaves a resolvable map error in this small smooth cap;
    # order eight already agrees at roundoff and cannot exercise this check.
    metric,rates,tangent=setup(budget,cubic=True,order=4)
    original=copy.deepcopy(metric['state'])
    result=compare_cap_radial_maps(metric,rates,tangent,budget=budget)
    direct=independent_variation(metric,rates,[r['panel'] for r in result['panels']],budget)
    np.testing.assert_allclose([result['reference_mass'],result['reference_rate']],direct,rtol=2e-10,atol=2e-12)
    for key in original: np.testing.assert_array_equal(metric['state'][key],original[key])
    assert abs(result['native_tangent_defect'])<=result['native_tangent_allowance']
    assert abs(result['error_rate'])>1e-12
    assert result['reference_coefficient_calls']>0


@pytest.mark.parametrize('rx',[.02,-.02])
def test_moving_cell_face_at_fixed_switch_retains_reference_born_pieces(rx,budget):
    rates=dict(R=rx,T=.3,Ip=np.zeros(4),Iq=np.zeros(4))
    metric,rates,tangent=setup(budget,radius=.04,cubic=True,
        state=dict(Ip=np.zeros(4),Iq=np.zeros(4)),rates=rates)
    result=compare_cap_radial_maps(metric,rates,tangent,budget=budget)
    assert sum(p['born'] for row in result['panels'] for p in row['pieces'])==1
    direct=independent_variation(metric,rates,[r['panel'] for r in result['panels']],budget)
    np.testing.assert_allclose([result['reference_mass'],result['reference_rate']],direct,rtol=2e-10,atol=2e-12)


def test_switch_split_has_both_correct_sides_and_born_map_derivatives(budget):
    panel=dict(left=SMALL_RADIUS,right=SMALL_RADIUS,left_rate=1.,right_rate=-1.,owner=0)
    pieces=_reference_pieces(panel)
    assert [p['small_radius_side'] for p in pieces]==['left','right']
    assert [(p['left_rate'],p['right_rate']) for p in pieces]==[(1.,0.),(0.,-1.)]
    class PiecewiseConstant:
        def reference_sample(self,owner,radii,*,budget,small_radius_side):
            budget.check(6)
            b=2. if small_radius_side=='left' else 5.
            return {k:np.array([b if k=='b' else 0.]) for k in ('a','b','a_x','b_x','a_r','b_r')}
    result,_,_= _reference_map(PiecewiseConstant(),panel,budget,2e-11,2e-13)
    assert result==dict(S=1.,T=0.,S_x=0.,T_x=-7.)
    with pytest.raises(ValueError,match='reversed'):
        _reference_pieces(panel|dict(left_rate=-1.,right_rate=1.))


def test_native_partition_is_frozen_and_mismatched_evidence_fails(budget):
    metric,rates,tangent=setup(budget)
    view=CapTransportFields(metric,rates,budget=budget)
    before=view.native_panels(budget=budget)
    metric['y_faces'][1]=.19
    rates['R']=9.
    assert view.native_panels(budget=budget)==before
    # Use an intact manufactured context for mismatched proof controls.
    metric,rates,tangent=setup(budget)
    wrong=copy.deepcopy(tangent); wrong['panels'][0]['owner']=3
    with pytest.raises(ValueError,match='owner'):
        compare_cap_radial_maps(metric,rates,wrong,budget=budget)
    wrong=copy.deepcopy(tangent); wrong['forward']=np.nan
    with pytest.raises(ValueError,match='native tangent'):
        compare_cap_radial_maps(metric,rates,wrong,budget=budget)


def test_native_born_switch_keeps_existing_tangent_domain_guard(budget):
    metric,rates,_=setup(budget)
    view=CapTransportFields(metric,rates,budget=budget)
    panel=dict(left=SMALL_RADIUS,right=SMALL_RADIUS,left_rate=1.,right_rate=-1.,owner=0)
    with pytest.raises(ValueError,match='small-radius formula switch'):
        _native_map(view,panel,metric['quadrature_order'],budget)
