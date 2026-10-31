"""Actual dense-polynomial differentiation and consistent null metric jets."""
import math

import numpy as np
import pytest
from scipy.integrate import DOP853

from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
from recursive_horizons.nsc_atlas_null_parent import dense_history_jet,GHNullHistorySampler
from recursive_horizons.nsc_atlas_parent_gh import ParentGH,PAIRS


def analytic_history():
    system=ParentGH(2,0);initial=system.initial();H=system.r;normalization=math.sqrt(2*math.pi**2)
    def rhs(t,state):
        C=system.unpack(state);derivative=np.zeros_like(C)
        for k,(a,b) in enumerate(PAIRS):
            derivative[k]=C[k+10]
            if a==b and a:derivative[k+10,0]=2*H*H*math.cosh(2*H*t)*normalization
        out=np.r_[system.pack(derivative),np.zeros(10)];out[system.dof+8]=1
        return out
    solver=DOP853(rhs,0.,initial,.4,rtol=1e-12,atol=1e-14,max_step=.05)
    history=DenseGHHistory()
    while solver.status=='running':
        solver.step();history.record(solver.t_old,solver.t,solver.dense_output())
    return system,history


def test_actual_dense_metric_polynomial_derivatives_match_analytic_background():
    system,history=analytic_history();H=system.r;normalization=math.sqrt(2*math.pi**2)
    index=PAIRS.index((1,1))
    for t in (.013,.087,.223,.397):
        y,dy,ddy=dense_history_jet(history,t)
        assert np.max(abs(y-history.evaluate(t)))<2e-15
        C,C1,C2=map(system.unpack,(y,dy,ddy))
        assert abs(C[index,0]/normalization-math.cosh(H*t)**2)<1e-12
        assert abs(C1[index,0]/normalization-H*math.sinh(2*H*t))<1e-10
        assert abs(C2[index,0]/normalization-2*H**2*math.cosh(2*H*t))<1e-8
    with pytest.raises(ValueError,match='extrapolation'):dense_history_jet(history,.401)


def test_consistent_metric_curvature_and_scalar_source_are_independent_controls():
    system,history=analytic_history();sampler=GHNullHistorySampler(system,history,physical_interval=(0.,.4))
    X=np.array([[0,0,0,1],[.5,.5,.5,.5]])
    data=sampler(.223,X);g=data['g'];H=system.r
    exact=H**2*(np.einsum('nac,nbd->nabcd',g,g)-np.einsum('nad,nbc->nabcd',g,g))
    assert np.max(abs(data['R_cov']-exact))<2e-9
    assert np.max(abs(data['metric_velocity_consistency']))<2e-10
    assert np.max(abs(data['Ricci_source_residual']))<4e-9
    assert np.max(abs(data['field_velocity_consistency']))==0
    assert sampler.calls==1 and sampler.points_evaluated==2
    with pytest.raises(ValueError,match='physical history'):sampler(.401,X)
