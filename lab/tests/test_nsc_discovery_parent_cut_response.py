"""Fixed-k implicit preparation, AP characteristics and short shared-clock response."""
from dataclasses import replace
import numpy as np
import pytest
from threadpoolctl import threadpool_limits
from recursive_horizons import nsc_discovery_parent_cut_response as cut


@pytest.fixture(scope='module')
def fixture():
    p=cut.parent;r=cut.response;le=cut.leading
    with threadpool_limits(limits=1),p.backend.fft_thread_limit(1):
        grid=r.install_native_fft(p.galerkin.build_grid(16,quadrature=64,gauge='conformal'))
        x=grid.xi_q-grid.length/2;g=np.exp(-x*x/.5);packet=np.exp(-(abs(x)-2.)**2/.2)
        raw=[p._project_spinor(grid,g*np.exp(.2*x),g*np.exp(-.2*x)),p._project_spinor(grid,packet,packet)]
        phi0=np.column_stack([a[0]/np.sqrt(a[2]) for a in raw]);phi1=np.column_stack([a[1]/np.sqrt(a[2]) for a in raw])
        columns=np.vstack((phi0,phi1));reference=np.linalg.qr(columns)[0]
        pair=r.make_holder(grid,weights=np.array([.02,.01]),reference_columns=reference,source_columns=columns)
        a,Z,mag=le.coefficients(pair.grid.fine);radius=np.sqrt(mag/abs(a))
        state=le.encode(pair,le.State(np.ones(grid.ng),np.full(grid.ng,radius),np.zeros(grid.ng),np.zeros(grid.ng),phi0,phi1))
        problem,seed,guess,theta=cut.problem(pair,state,.1,-1)
        state=p._state_from_theta(pair,problem,guess)
        pair,state,_=cut.reprepare(pair,state,.1,-1,pair.weights,cpu_limit=3.)
    return pair,state,.1


def test_analytic_initial_actual_momentum_constraints_match_two_widths(fixture):
    pair,state,k=fixture;before={n:getattr(state,n).tobytes() for n in cut.leading.FIELDS}
    with threadpool_limits(limits=1),cut.parent.backend.fft_thread_limit(1):
        tangent,report=cut.initial_tangent(pair,state,k,width=.001,checks=True,cpu_limit=5.)
    assert report['analytic_rank_usable'] and report['rank']==2*((pair.grid.ng+1)//2)
    assert report['linear_residual_max']<1e-12
    assert len(report['finite_difference_checks'])==2
    assert max(row['momentum_max_gap'] for row in report['finite_difference_checks'])<1e-8
    assert np.array_equal(tangent.occupations,[0.,pair.weights[1]])
    for n in ('Q','r','phi0','phi1'):assert not np.any(getattr(tangent,n))
    assert before=={n:getattr(state,n).tobytes() for n in before}


def test_AP_sampler_and_incoming_projectors_do_not_periodic_wrap(fixture):
    pair,state,k=fixture;nf=pair.grid.nf;dx=pair.grid.length/nf;x=np.arange(nf)*dx
    wave=2*np.pi*1.5/pair.grid.length;phi=np.sqrt(dx)*np.exp(1j*wave*x)[:,None]
    points=np.array([0.,.3,8.,8.3])
    np.testing.assert_allclose(cut.ap_values(phi,pair.grid.length,points)[:,0],np.exp(1j*wave*points),atol=2e-15)
    assert cut.ap_values(phi,8.,[8.])[0,0]==pytest.approx(-cut.ap_values(phi,8.,[0.])[0,0],abs=2e-15)
    # A sigma2-positive wave is incoming at left, outgoing at right.
    field=np.sqrt(dx)*np.exp(1j*wave*x)[:,None]/np.sqrt(2)
    source=np.vstack((field,1j*field));reference=source/np.linalg.norm(source)
    holder=cut.response.make_holder(pair.grid,weights=np.array([.02]),reference_columns=reference,source_columns=reference)
    nodal=cut.leading.decode(pair,state)
    test=cut.leading.encode(holder,cut.leading.State(nodal.Q,nodal.r,nodal.p_Q,nodal.p_r,reference[:nf],reference[nf:]))
    values,_=cut.observables(holder,test)
    assert values['left_Dirac_incoming_normal_current']>0
    assert abs(values['right_Dirac_incoming_normal_current'])<1e-25


def test_short_retarded_tangent_matches_full_RK_and_accumulated_clock(fixture):
    pair,state,k=fixture;h=.001
    with threadpool_limits(limits=1),cut.parent.backend.fft_thread_limit(1):
        tangent,_=cut.initial_tangent(pair,state,k,checks=False)
        up,plus,_=cut.reprepare(pair,state,k,-1,pair.weights+h*tangent.occupations)
        down,minus,_=cut.reprepare(pair,state,k,-1,pair.weights-h*tangent.occupations)
        baseline=state.copy();tau=dtau=tu=td=0.
        for _ in range(3):
            baseline,tangent,clock=cut.response.advance_tangent(pair,baseline,tangent,.001)
            tau+=clock['tau'];dtau+=clock['delta_tau']
            tu+=cut.response._nonlinear_centre_increment(up,plus,.001);td+=cut.response._nonlinear_centre_increment(down,minus,.001)
            plus=cut.leading.rk4_step(up,plus,.001);minus=cut.leading.rk4_step(down,minus,.001)
        for n in cut.leading.FIELDS:
            np.testing.assert_allclose(getattr(tangent,n),(getattr(plus,n)-getattr(minus,n))/(2*h),rtol=2e-5,atol=2e-9)
        assert dtau==pytest.approx((tu-td)/(2*h),abs=1e-10)
        assert tau>clock['tau']
        values,delta=cut.observables(pair,baseline,tangent)
        upper,_=cut.observables(up,plus);lower,_=cut.observables(down,minus)
        for key in values:assert delta[key]==pytest.approx((upper[key]-lower[key])/(2*h),abs=2e-8)
        row=cut.row(pair,baseline,tangent,.003,tau,dtau)
        _,slope=cut.observables(pair,baseline,cut.response.velocity_tangent(pair,baseline))
        for key in values:
            assert row['delta_equal_tau'][key]==pytest.approx(delta[key]-slope[key]*dtau/values['centre_clock_rate'],abs=1e-14)


def test_protected_cut_geometry_incoming_signs_are_distinct_from_child_window(fixture):
    pair,state,k=fixture;rate,bundle=cut.leading.rates(pair,state,return_bundle=True)
    fine=bundle['fine_state'];velocity=cut.leading.prolong(pair.grid,cut.leading.decode(pair,rate));D=bundle['fine_system'].derivative
    values,_=cut.observables(pair,state)
    assert cut.CUTS==(3.,5.) and pair.child_interval==(3.5,4.5)
    ev=lambda a:cut.extent.real_periodic_values(pair.grid,a,cut.CUTS)
    expected=(ev(velocity.r)+np.array([-1.,1.])*ev(D@fine.r))/ev(fine.r*fine.Q)
    assert values['left_r_incoming_normal']==pytest.approx(expected[0],abs=1e-14)
    assert values['right_r_incoming_normal']==pytest.approx(expected[1],abs=1e-14)
