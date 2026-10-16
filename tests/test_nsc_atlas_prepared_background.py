"""Manufactured interface controls; no prepared physical source solve."""
import time
from types import SimpleNamespace

import numpy as np
from numpy.polynomial import Polynomial
from scipy.interpolate import CubicHermiteSpline
import pytest

from recursive_horizons.nsc_atlas_cap_field import CapLeadingField
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget
from recursive_horizons.nsc_atlas_background_profile import build_prepared_background


@pytest.fixture(scope='module')
def budget():
    value=PreparedFieldBudget(cpu_seconds=20.,kernel_values_max=100000)
    yield value
    value.check()
    print('prepared history adapter controls:',value.kernel_values,
          'charged values;',time.process_time()-value.started,'CPU seconds')


class AnalyticPreparedCap(CapLeadingField):
    """Actual source polynomial in an affine clock; no field evolution."""
    def __init__(self,budget):
        self.budget=budget;self.Va=0.;self.Vb=1.;self.dV=1.
        nodes=np.array([0.,.25,.6,1.])
        self.geometry={'arrival':CubicHermiteSpline(nodes,1+2*nodes,np.full(4,2.))}
        self.prepared_field=SimpleNamespace(dense_data={'cuts':np.array([0.,.25,.5,.75,1.])})
        self.domain=dict(qa=1.,qb=3.,qmin=0.,Ta=.1,Tb=.2)
        self.quiet_certificate={'quiet_cap':True}
        self.max_inverse_bracket_width=0.;self.calls=[];self.quiet_noise=0.

    def profile(self,q,*,second=False):
        q=np.array(q,dtype=float);self.calls.append((q.copy(),second))
        V=(q-1)/2;z=Polynomial([-2,4]);P=(Polynomial([1])-z*z)**4
        active=(V>.25)&(V<.75)
        row={k:np.where(active,P.deriv(i)(V)/2**i,0.)
             for i,k in enumerate(('f','f_prime','f_second'))}
        row['f'][q<=1.5]+=self.quiet_noise
        self.budget.check(3*len(q))
        return row


def test_native_knots_derivatives_and_source_coefficient_are_retained(budget,monkeypatch):
    cap=AnalyticPreparedCap(budget)
    def forbidden(*args,**kwargs):raise AssertionError('waveform differences are not an input')
    monkeypatch.setattr(np,'gradient',forbidden)
    before=budget.evaluations
    result=build_prepared_background(cap)
    assert result['charged_evaluations']==budget.evaluations-before
    np.testing.assert_array_equal(result['source_null_labels'],[0.,.25,.5,.6,.75,1.])
    assert len(cap.calls)==1 and cap.calls[0][1] is True
    assert result['onset']['quartic_coefficient']==256.
    assert result['onset']['source_V_roundoff']==0.
    bg=result['background'];sample=cap.profile(bg.inputs['labels'],second=True)
    for key,name in (('values','f'),('first','f_prime'),('second','f_second')):
        np.testing.assert_array_equal(bg.inputs[key],sample[name])
    onset=result['onset']['label']
    assert bg.profile._spline.derivative(4)(onset)==pytest.approx(24*256,rel=2e-14)
    assert np.all(bg.profile._spline.c[:,0]==0.)
    assert bg.profile.quiet_extension
    assert result['physical_input_domain']['qmin']==0.
    assert not result['numerical_error_enclosure']
    assert max(result['diagnostics']['knot_jet_discrepancies'])<2e-10
    with pytest.raises(ValueError,match='supplied finite'):
        bg.profile.jet(np.nextafter(3.,np.inf),budget=budget)


def test_quiet_discrepancy_is_rejected_or_explicitly_retained(budget):
    cap=AnalyticPreparedCap(budget);cap.quiet_noise=1e-8
    with pytest.raises(ValueError,match='quiet-jet discrepancy'):
        build_prepared_background(cap)
    out=build_prepared_background(cap,quiet_tolerances=(2e-8,0.,0.))
    np.testing.assert_array_equal(out['diagnostics']['quiet_jet_discrepancies'],[1e-8,0.,0.])
    assert np.all(out['background'].inputs['values'][:2]==0.)
    assert out['background'].onset['quartic_coefficient']==256.


def test_inputs_are_copied_and_budget_counter_has_one_owner(budget):
    cap=AnalyticPreparedCap(budget);out=build_prepared_background(cap)
    labels=out['background'].inputs['labels'].copy()
    cap.prepared_field.dense_data['cuts'][:]=99.
    cap.calls[0][0][:]=99.
    np.testing.assert_array_equal(out['background'].inputs['labels'],labels)
    before=budget.evaluations;budget.check(11)
    assert budget.evaluations==budget.kernel_values==before+11
    with pytest.raises(AttributeError): budget.evaluations=0


def test_knot_capacity_and_domain_fail_before_profile_query(budget):
    cap=AnalyticPreparedCap(budget);nodes=np.linspace(0.,1.,513)
    cap.geometry['arrival']=CubicHermiteSpline(nodes,1+2*nodes,np.full(len(nodes),2.))
    cap.prepared_field.dense_data['cuts']=np.array([0.,.3,1.])
    with pytest.raises(ValueError,match='exceeds 513'):build_prepared_background(cap)
    assert not cap.calls
    cap=AnalyticPreparedCap(budget);cap.domain['qb']=4.
    with pytest.raises(ValueError,match='finite domain'):build_prepared_background(cap)
    assert not cap.calls
    for values in ((-1.,0.,0.),(1.,2.),(True,0.,0.)):
        with pytest.raises(ValueError):build_prepared_background(AnalyticPreparedCap(budget),quiet_tolerances=values)
    with pytest.raises(ValueError):build_prepared_background(object())


def test_budget_exhaustion_precedes_profile_query(budget):
    cap=AnalyticPreparedCap(budget);limit=budget.kernel_values_max
    budget.kernel_values_max=budget.kernel_values+1
    try:
        with pytest.raises(RuntimeError,match='budget exhausted'):build_prepared_background(cap)
        assert not cap.calls
    finally:budget.kernel_values_max=limit


def test_enrichment_queries_existing_profile_without_changing_knots_or_onset(budget):
    cap=AnalyticPreparedCap(budget);base=build_prepared_background(cap);cap.calls.clear()
    before=budget.evaluations;out=build_prepared_background(cap,enrich_ordinary_pieces=True)
    assert out['charged_evaluations']==budget.evaluations-before
    assert len(cap.calls)==2 and cap.calls[0][1] is True and cap.calls[1][1] is False
    bg=out['background'];extra=bg.interior_jets
    assert out['enrichment']['ordinary_pieces']==3
    for key in ('labels','values','first','second'):
        np.testing.assert_array_equal(bg.inputs[key],base['background'].inputs[key])
        np.testing.assert_array_equal(out['original_history_inputs'][key],base['background'].inputs[key])
    np.testing.assert_array_equal(bg.profile._spline.c[:,:2],base['background'].profile._spline.c[:,:2])
    np.testing.assert_array_equal(extra['labels'],cap.calls[1][0])
    original=cap.profile(extra['labels'])
    computed=bg.profile.jet(extra['labels'],budget=budget)
    np.testing.assert_allclose(computed['F'],original['f'],atol=3e-13,rtol=3e-13)
    np.testing.assert_allclose(computed['F_prime'],original['f_prime'],atol=3e-12,rtol=3e-12)
