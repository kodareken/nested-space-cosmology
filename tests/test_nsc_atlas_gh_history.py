"""Dense archives agree with the actual solver and do not steer its state."""
import numpy as np
import pytest
from scipy.integrate import DOP853

from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
from recursive_horizons.nsc_atlas_parent_gh import evolve_parent


def test_polynomial_round_trip_and_no_extrapolation(tmp_path):
    solver=DOP853(lambda t,y:np.array([y[1],-4*y[0]]),0.,[1.,.3],t_bound=.8,max_step=.1,rtol=1e-11,atol=1e-13)
    history=DenseGHHistory();checks=[]
    while solver.status=='running':
        solver.step();dense=solver.dense_output();history.record(solver.t_old,solver.t,dense)
        for x in (.13,.71,1.):
            t=solver.t_old+x*(solver.t-solver.t_old);checks.append((t,dense(t)))
    path=tmp_path/'history.npz';history.save(path,{'example':True});loaded=DenseGHHistory.load(path)
    for t,y in checks:np.testing.assert_array_equal(loaded.evaluate(t),y)
    assert loaded.metadata=={'example':True}
    for t in (-.01,.81,float('nan')):
        with pytest.raises(ValueError):loaded.evaluate(t)
    with pytest.raises(ValueError,match='fresh'):history.save(path,{})


def test_dense_sampler_leaves_nonzero_parent_state_unchanged():
    captured=[];results=[]
    history=DenseGHHistory()
    for enabled in (False,True):
        states=[]
        result=evolve_parent(degree=2,epsilon=.05,threads=1,proper_endpoint=.3,budget_seconds=90,
            state_sink=lambda t,y:states.append(y.copy()),dense_sink=history.record if enabled else None)
        captured.append(states);results.append(result)
    assert results[0]['steps']==results[1]['steps']
    for a,b in zip(*captured):np.testing.assert_array_equal(a,b)
    assert len(history.bounds)==results[1]['steps']
    t=results[1]['snapshots'][-1]['coordinate_t']
    np.testing.assert_allclose(history.evaluate(t)[:-1],captured[1][-1][:-1],atol=0,rtol=0)
