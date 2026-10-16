"""Independent energy, mass sensitivities and strong/weak fixed-cap traces."""
from fractions import Fraction as Q

import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_fixed_basis import FixedCapBasis
from recursive_horizons.nsc_atlas_cap_fixed_reference import _delta,_Fields,fixed_cap_reference

CONTROLS=ControlBudget(maximum_evaluations=5_000_000,cpu_seconds=60.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value: yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('fixed mass reference controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def case(budget):
    basis=FixedCapBasis([0,Q(1,4),Q(1,2)],budget=budget)
    c=np.array([1,2,-1,3,-2,1,-1])*.002
    d=np.array([-1,2,3,1,2,-2,1])*.003
    return basis,c,d,Q(3,8)


def reference(basis,c,d,Z,budget,**changes):
    args=dict(alpha=.002,T_x=1.2,radius_rate=.001,incoming_w_plus=.03,budget=budget)
    args.update(changes)
    return fixed_cap_reference(basis,c,d,Z,**args)


def exact_operators(basis,Z,budget):
    n=basis.dimension
    M=[[Q(0) for _ in range(n)] for _ in range(n)]
    K=[[Q(0) for _ in range(n)] for _ in range(n)]
    e=[Q(0)]*n
    for cell in range(basis.cells):
        low,high=basis.exact_knots[cell:cell+2];width=high-low
        cut=min(Q(1),(Z-low)/width);idx=basis.indices(cell);s=basis.exact_shapes(cell)
        for i,gi in enumerate(idx):
            if cell==basis.cells-1: e[int(gi)]=sum(s[p][i]*cut**p for p in range(6))
            for j,gj in enumerate(idx):
                M[gi][gj]+=width*sum(s[p][i]*s[q][j]*cut**(p+q+1)/Q(p+q+1)
                                    for p in range(6) for q in range(6))
                K[gi][gj]+=sum(p*q*s[p][i]*s[q][j]*cut**(p+q-1)/Q(p+q-1)
                               for p in range(1,6) for q in range(1,6))/width
        budget.check(5000)
    return np.array(M,dtype=float),np.array(K,dtype=float),np.array(e,dtype=float)


def test_delta_and_axis_reduced_fields_preserve_their_exact_limits(budget):
    power=tuple(Q(v,7) for v in (2,-3,5,-1,4,6));left,width=Q(1,8),Q(1,4)
    delta=_delta(power,left,width)
    for xi in (Q(i,7) for i in range(7)):
        p=sum(v*xi**j for j,v in enumerate(power))
        derivative=sum(j*power[j]*xi**(j-1)/width for j in range(1,6))
        assert sum(v*xi**j for j,v in enumerate(delta))==(left+width*xi)*derivative-p
    basis=FixedCapBasis([0,Q(1,2)],budget=budget)
    c=np.array([.5,.5,.5]);d=2*c
    fields=_Fields(basis,c,d,budget)
    axis=fields.evaluate(0,0.,budget)
    assert (axis['f'],axis['p'],axis['h'])==(1.,0.,0.)
    tiny=fields.evaluate(0,1e-7,budget);z=tiny['z']
    assert tiny['h']>0
    assert abs(tiny['h']/((2/3)*z*z)-1)<1e-12
    budget.check(2000)


def test_zero_coupling_matches_independently_integrated_energy_and_rates(budget):
    basis,c,d,Z=case(budget);M,K,e=exact_operators(basis,Z,budget)
    value=reference(basis,c,d,Z,budget,alpha=0.)
    R=np.tanh(float(Z));A0=1/np.cosh(float(Z))**2;beta=A0/R
    stiffness=K-beta*np.outer(e,e)
    mass=(c@stiffness@c+d@M@d)/4
    gc=stiffness@c/2;gd=M@d/2
    np.testing.assert_allclose(value['mass'],mass,rtol=2e-10,atol=2e-12)
    np.testing.assert_allclose(value['gradient']['chi'],gc,rtol=2e-9,atol=2e-10)
    np.testing.assert_allclose(value['gradient']['canonical_density'],gd,rtol=2e-9,atol=2e-10)
    np.testing.assert_allclose(value['rates']['chi'],1.2*d,rtol=2e-8,atol=2e-10)
    expected_d=1.2*np.linalg.solve(M,-2*gc+e*(A0*R*.03-e@d))
    np.testing.assert_allclose(value['rates']['canonical_density'],expected_d,rtol=3e-8,atol=3e-9)
    assert abs(value['weak_balance_defect'])<1e-11
    assert abs(value['wall']['momentum_discrepancy'])<1e-10
    assert abs(value['wall']['normal_discrepancy'])>1e-3
    assert abs(value['wall']['incoming_discrepancy']-(value['wall']['normal_discrepancy']+
               value['wall']['momentum_discrepancy'])/R)<1e-12


def test_finite_coefficient_and_shape_derivatives_are_independent_mass_differences(budget):
    basis,c,d,Z=case(budget);c0,d0=c.copy(),d.copy()
    base=reference(basis,c,d,Z,budget)
    vc=np.array([1,-2,1,0,1,-1,2])/8
    vd=np.array([-1,0,2,1,-1,1,0])/8
    step=2e-4
    plus=reference(basis,c+step*vc,d+step*vd,Z,budget)
    minus=reference(basis,c-step*vc,d-step*vd,Z,budget)
    observed=(plus['mass']-minus['mass'])/(2*step)
    predicted=base['gradient']['chi']@vc+base['gradient']['canonical_density']@vd
    assert abs(observed-predicted)<2e-8
    shape_step=Q(1,100000)
    upper=reference(basis,c,d,Z+shape_step,budget)
    lower=reference(basis,c,d,Z-shape_step,budget)
    observed_shape=(upper['mass']-lower['mass'])/(2*float(shape_step))
    assert abs(observed_shape-base['shape_derivative_Z'])<5e-8
    assert abs(base['weak_balance_defect'])<1e-11
    assert base['gradient_discrepancy']<base['gradient_absolute_tolerance']
    expected_lapse=-np.log(np.cosh(np.sqrt(.002/6)*base['centre']['f']))
    assert abs(base['centre']['delta']-expected_lapse)<1e-15
    np.testing.assert_array_equal(c,c0);np.testing.assert_array_equal(d,d0)


def test_propagated_gradient_difference_is_small_in_physical_norm(budget):
    basis,c,d,Z=case(budget);M,_,_=exact_operators(basis,Z,budget)
    gram=basis.gram(Z,budget=budget)
    np.testing.assert_allclose(gram.squared_norm(c,budget=budget),c@M@c,rtol=2e-12,atol=1e-18)
    assert gram.squared_norm(np.zeros(basis.dimension),budget=budget)==0.
    with pytest.raises(ValueError,match='dimension'):
        gram.squared_norm(c[:-1],budget=budget)
    value=reference(basis,c,d,Z,budget)
    discrepancy=value['projected_gradient_discrepancy']
    assert discrepancy['chi_rate_L2']<2e-10
    assert discrepancy['canonical_density_rate_L2']<2e-10
    assert abs(discrepancy['wall_velocity'])<2e-10
    print('projected gradient comparison:',discrepancy)


def test_quiet_total_fields_recover_empty_cap(budget):
    basis=FixedCapBasis([0,Q(1,2)],budget=budget);zero=np.zeros(basis.dimension)
    value=reference(basis,zero,zero,Q(1,4),budget,incoming_w_plus=0.)
    assert value['mass']==value['lapse']==value['mass_rate']==0
    assert value['shape_derivative_Z']==0
    assert all(np.all(v==0) for v in value['rates'].values())
    assert value['wall']['normal_discrepancy']==value['wall']['momentum_discrepancy']==0


def test_bad_domains_and_nontimelike_wall_are_rejected(budget):
    basis,c,d,Z=case(budget)
    for changes in (dict(T_x=0.),dict(alpha=-1.),dict(rtol=1e-16),dict(atol=1.)):
        with pytest.raises(ValueError): reference(basis,c,d,Z,budget,**changes)
    with pytest.raises(ValueError,match='dimensions'): reference(basis,c[:-1],d,Z,budget)
    with pytest.raises(ValueError,match='wall cut'): reference(basis,c,d,Q(1,4),budget)
    with pytest.raises(ValueError,match='not timelike'): reference(basis,c,d,Z,budget,radius_rate=2.)
