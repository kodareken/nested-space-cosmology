"""Known integrals, rejected stages, dense cache and exact passivity checks."""
import math

import numpy as np
import pytest
from scipy.integrate import DOP853

from recursive_horizons.nsc_atlas_passive_dop853 import PassiveDOP853Quadrature


def test_initial_step_probe_does_not_replace_true_initial_work():
    Q=PassiveDOP853Quadrature(initial=.3)
    def fun(t,y):
        Q.record(math.exp(5000*t))
        return np.ones_like(y)
    solver=DOP853(fun,0.,[0.],.001,rtol=1e-10,atol=1e-12)
    assert solver.nfev==2  # f(t0) and the separate select_initial_step probe.
    Q.initialize(solver);solver.step();Q.accepted_step(solver)
    expected=.3+math.expm1(5000*solver.t)/5000
    assert Q.value==pytest.approx(expected,abs=4e-14)


def test_known_sine_integral_and_dense_values_add_no_rhs_calls():
    Q=PassiveDOP853Quadrature();omega=7.
    def fun(t,y):
        Q.record(math.sin(omega*t))
        return np.array([omega*math.cos(omega*t)])
    solver=DOP853(fun,0.,[0.],1.7,rtol=3e-11,atol=2e-13)
    Q.initialize(solver)
    while solver.status=='running':
        solver.step();Q.accepted_step(solver)
        before=solver.nfev
        solver.dense_output()
        assert solver.nfev==before+3
        times=solver.t_old+np.array([0.,.2,.5,.8,1.])*(solver.t-solver.t_old)
        count=solver.nfev
        actual=Q.dense_value(solver,times)
        exact=(1-np.cos(omega*times))/omega
        np.testing.assert_allclose(actual,exact,atol=2e-11,rtol=0)
        assert solver.nfev==count
        assert Q.dense_value(solver,solver.t)==Q.value
        assert len(Q._records)<=16
    assert Q.value==pytest.approx((1-math.cos(omega*solver.t))/omega,abs=2e-11)


def test_rejected_trials_and_old_dense_stages_are_not_integrated():
    Q=PassiveDOP853Quadrature();omega=15.
    def fun(t,y):
        Q.record(y[0]*y[1])  # Work oscillates at twice the controlled frequency.
        return np.array([y[1],-omega*omega*y[0]])
    solver=DOP853(fun,0.,[1.,0.],1.1,first_step=.9,rtol=3e-12,atol=3e-14)
    Q.initialize(solver);rejected=False;steps=0
    while solver.status=='running':
        before=solver.nfev
        solver.step();rejected|=(solver.nfev-before)>12
        Q.accepted_step(solver)
        if steps%2==0:
            # Some diagnostic dense outputs are not used by the quadrature.
            # Their extra three work records must be ignored at the next step.
            solver.dense_output()
            if steps%4==0:
                midpoint=(solver.t_old+solver.t)/2
                assert Q.dense_value(solver,midpoint)==pytest.approx(
                    .5*(math.cos(omega*midpoint)**2-1),abs=3e-10)
        steps+=1
        assert len(Q._records)<=16
    assert rejected
    assert Q.value==pytest.approx(.5*(math.cos(omega*solver.t)**2-1),abs=3e-10)


def integrate_controller(with_quadrature):
    omega=9.;Q=PassiveDOP853Quadrature() if with_quadrature else None
    def fun(t,y):
        work=y[0]*y[1]
        if Q is not None:Q.record(work)
        return np.array([y[1],-omega*omega*y[0]])
    solver=DOP853(fun,0.,[1.,0.],1.4,first_step=.8,rtol=1e-10,atol=1e-12)
    if Q is not None:Q.initialize(solver)
    steps=[];states=[];counter=0
    while solver.status=='running':
        solver.step()
        if Q is not None:Q.accepted_step(solver)
        steps.append(solver.t);states.append(solver.y.copy())
        if counter%3==0:
            solver.dense_output()
            if Q is not None:Q.dense_value(solver,(solver.t_old+solver.t)/2)
        counter+=1
    return np.array(steps),np.array(states),solver.nfev


def test_quadrature_leaves_all_adaptive_steps_states_and_rhs_counts_bitwise_unchanged():
    baseline=integrate_controller(False);passive=integrate_controller(True)
    np.testing.assert_array_equal(baseline[0],passive[0])
    np.testing.assert_array_equal(baseline[1],passive[1])
    assert baseline[2]==passive[2]


def test_dense_cache_is_current_interval_only_and_can_be_reused_or_rebuilt():
    Q=PassiveDOP853Quadrature()
    def fun(t,y):
        Q.record(math.sin(t))
        return np.array([math.cos(t)])
    solver=DOP853(fun,0.,[0.],.5,max_step=.1,first_step=.1,rtol=1e-10,atol=1e-12)
    Q.initialize(solver);solver.step();Q.accepted_step(solver)
    with pytest.raises(RuntimeError,match='dense_output'):Q.dense_value(solver,.05)
    solver.dense_output();first=Q.dense_value(solver,.05)
    nfev=solver.nfev
    assert Q.dense_value(solver,.05)==first
    assert solver.nfev==nfev
    solver.dense_output()
    assert Q.dense_value(solver,.05)==first
    solver.step();Q.accepted_step(solver)
    with pytest.raises(RuntimeError,match='dense_output'):Q.dense_value(solver,solver.t)
    solver.dense_output()
    with pytest.raises(ValueError,match='last accepted interval'):Q.dense_value(solver,.05)
    with pytest.raises(RuntimeError,match='already recorded'):Q.accepted_step(solver)


def test_backward_steps_use_the_signed_quadrature_interval():
    Q=PassiveDOP853Quadrature(initial=.2)
    def fun(t,y):
        Q.record(math.sin(3*t))
        return np.array([math.cos(3*t)])
    solver=DOP853(fun,1.,[0.],0.,rtol=1e-11,atol=1e-13)
    Q.initialize(solver)
    while solver.status=='running':
        solver.step();Q.accepted_step(solver);solver.dense_output()
        t=(solver.t_old+solver.t)/2
        expected=.2+(math.cos(3)-math.cos(3*t))/3
        assert Q.dense_value(solver,t)==pytest.approx(expected,abs=3e-11)
    assert Q.value==pytest.approx(.2+(math.cos(3)-1)/3,abs=3e-11)
