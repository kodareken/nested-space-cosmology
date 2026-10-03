"""Rank-two source normalization, actual constraints, metric and bounded preparation."""
from dataclasses import replace
import hashlib
import json
import time
import numpy as np
import pytest
from threadpoolctl import threadpool_limits
from recursive_horizons import nsc_discovery_parent as parent


@pytest.fixture(scope='module')
def setup():
    with threadpool_limits(limits=1),parent.backend.fft_thread_limit(1):
        table=parent.parent_continuum(cpu_limit=5.)
        grid=parent.backend.make_fft_grid(parent.galerkin.build_grid(32,gauge='conformal'))
        r,Q=parent._band_geometry(grid,table);columns=parent._source_columns(grid,table)
        gram=columns['source_columns'].conj().T@columns['source_columns']
        weights,meta=parent.population_weights(columns['probabilities'],columns['raw_norms'],gram,1)
        grid.fine=replace(grid.fine,occupations=weights)
        pair=parent._holder(grid,weights,columns,meta,table,1)
    return table,pair,r,Q,columns


def test_project_then_normalize_and_distinct_observer(setup):
    table,pair,r,Q,col=setup
    gram=col['source_columns'].conj().T@col['source_columns']
    np.testing.assert_allclose(np.diag(gram),1.,atol=2e-15)
    assert abs(gram[0,1])>1e-5
    np.testing.assert_allclose(col['reference_columns'].conj().T@col['reference_columns'],np.eye(2),atol=2e-15)
    assert not np.array_equal(col['reference_columns'],col['source_columns'])
    base=parent.SEED_WEIGHT*col['raw_norms'][0]*col['unit_child_mass'][0]
    assert base==pytest.approx(parent.SEED_WEIGHT*col['raw_child_mass'],abs=1e-17)
    assert table['collar_radius_gap']<1e-12 and table['collar_Q_gap']<1e-12
    assert table['dirac_extension']['follows_blended_Q_outside_core']
    assert parent.backend.operator_backend(pair.grid)=='fft'


def test_actual_regions_fixed_trace_symmetric_population_intervention(setup):
    _,pair,r,Q,col=setup;gram=col['source_columns'].conj().T@col['source_columns'];values=[]
    for population in parent.POPULATIONS:
        weights,meta=parent.population_weights(col['probabilities'],col['raw_norms'],gram,population)
        a,b=col['probabilities'];delta=np.dot(weights,a-b)/np.dot(weights,a+b)
        assert delta==pytest.approx(meta['target_imbalance'],abs=1e-15)
        assert weights.sum()==pytest.approx(2*parent.SEED_WEIGHT*col['raw_norms'][0],abs=1e-17)
        assert np.min(weights)>=0 and np.max(parent._weighted_eigenvalues(gram,weights))<1
        values.append((weights,meta))
    assert values[0][1]['imbalance']==pytest.approx(-values[2][1]['imbalance'],abs=1e-15)
    assert values[1][1]['imbalance']==pytest.approx(0.,abs=1e-15)
    assert abs(values[1][0][0]-values[1][0][1])>1e-6
    with pytest.raises(ValueError,match='do not straddle'):
        parent.population_weights(np.array([[.5,.6],[.1,.2]]),np.ones(2),np.eye(2),1)


def problem_for(setup,sign=1):
    _,pair,r,Q,col=setup
    seed=parent._offset_seed(pair,r,Q,col['phi0'],col['phi1'],.05)
    problem=parent._MomentumProblem(pair,r,Q,col['phi0'],col['phi1'],sign,0.,seed['g'])
    theta=parent._seed_theta(problem,seed,sign)
    return problem,seed,theta


def test_actual_finite_momentum_Jacobian_physical_metric_and_anchor(setup):
    with threadpool_limits(limits=1),parent.backend.fft_thread_limit(1):
        problem,seed,theta=problem_for(setup)
        direction=np.random.default_rng(11).normal(size=theta.size);direction/=np.linalg.norm(direction)
        jac=problem.jacobian(theta);h=1e-6
        numeric=(problem.residual(theta+h*direction)-problem.residual(theta-h*direction))/(2*h)
        np.testing.assert_allclose(jac@direction,numeric,rtol=2e-6,atol=2e-7)
        np.testing.assert_allclose(problem.whitener.T@problem.metric@problem.whitener,np.eye(theta.size),atol=2e-14)
        p,v=problem.fine_momenta(theta)
        assert problem.anchor==pytest.approx(np.mean(p*p/problem.fine_r**3),abs=1e-14)
        assert abs(problem.residual(theta)[-1])<1e-14
        # The anchor uses the actual projected seed; it is not set to k.
        assert abs(problem.anchor-seed['suggested_k'])>1e-8
        stored=problem.anchor
        theta,report=parent.correct_momenta(problem,theta,parent._cpu_time()+2.,max_accepted=1)
        assert problem.anchor==stored and report['accepted']<=1 and report['physical_metric']


def test_both_signs_preserve_source_bytes_and_actual_source_current(setup):
    with threadpool_limits(limits=1),parent.backend.fft_thread_limit(1):
        plus,seed,tp=problem_for(setup,1);minus,_,tm=problem_for(setup,-1)
        sp=parent._state_from_theta(plus.pair,plus,tp);sm=parent._state_from_theta(minus.pair,minus,tm)
        assert sp.phi0.tobytes()==sm.phi0.tobytes() and sp.phi1.tobytes()==sm.phi1.tobytes()
        np.testing.assert_allclose(sp.p_Q,-sm.p_Q,atol=1e-13)
        np.testing.assert_allclose(sp.p_r,-sm.p_r,atol=1e-13)
        fine=parent.leading.fine_state(plus.pair,sp);system=parent.leading.active_system(plus.grid,fine)
        source=parent.coupling.source_from_columns(system,fine)
        C,D=parent.leading.constraint_arrays(plus.grid,fine,system,source)
        assert np.max(abs(source['force_L']))>0
        assert 'force_beta' in source and np.isfinite(D).all()
        restriction=parent.step_control.scaled_step_admission(plus.pair,sp,.001)
        assert restriction['source_rank']==2 and restriction['local_real_variables']==12


def test_expansion_uses_projected_velocity_and_normal_proper_lapse(setup):
    with threadpool_limits(limits=1),parent.backend.fft_thread_limit(1):
        problem,seed,theta=problem_for(setup)
        state=parent._state_from_theta(problem.pair,problem,theta)
        pair=problem.pair;rate,bundle=parent.leading.rates(pair,state,return_bundle=True)
        _,velocity,_=parent.leading.metric_jets(pair,state,rate,bundle);fine=bundle['fine_state']
        assert np.max(abs(bundle['unprojected_rates'][1]-velocity.r))>1e-8
        coordinate=velocity.r/fine.r;normal=coordinate/(fine.r*fine.Q)
        radial=(coordinate+velocity.Q/fine.Q)/(fine.r*fine.Q)
        diagnosed=parent._diagnose(pair,state,problem)
        for key,expected in (('expansion',normal),('coordinate_angular_rate',coordinate),('normal_radial_rate',radial)):
            for statistic,value in parent._summary(expected).items():
                assert diagnosed[key][statistic]==pytest.approx(value,abs=1e-12)
        center=parent.extent.real_periodic_values(pair.grid,normal,(parent.CENTER,))[0]
        assert diagnosed['expansion']['centre']==pytest.approx(center,abs=1e-12)
        assert diagnosed['raw_unprojected_rates_used_for_acceptance'] is False
        assert diagnosed['rates']['r']['max_abs']==pytest.approx(np.max(abs(velocity.r)),abs=1e-12)
        # At fixed r and P, scaling Q doubles coordinate rdot and N equally.
        doubled=state.copy();doubled.Q*=2
        scaled=parent._diagnose(pair,doubled,problem)
        assert scaled['coordinate_angular_rate']['max_abs']==pytest.approx(2*diagnosed['coordinate_angular_rate']['max_abs'],abs=1e-12)
        for statistic in ('min','max','max_abs','centre'):
            assert scaled['expansion'][statistic]==pytest.approx(diagnosed['expansion'][statistic],abs=1e-12)


def test_bounded_small_constructor_and_creation_only_temp_record(tmp_path,monkeypatch):
    before={str(parent.collar.OUTPUT/name):parent._sha_file(parent.collar.OUTPUT/name) for name in ('measurement.json','measurement.npz')}
    with threadpool_limits(limits=1),parent.backend.fft_thread_limit(1):
        pair,state,report=parent.prepare_parent(32,population=1,cpu_limit=5.)
    assert state.phi0.shape==(32,2) and report['evolved'] is False and report['trajectory'] is False
    assert report['diagnosis'] is not None and report['diagnosis']['timestep']['source_rank']==2
    assert report['correction'].get('accepted_total',report['correction']['accepted'])<=12
    assert report['fallback']['steps']<=parent.FALLBACK_STEPS
    monkeypatch.setattr(parent,'producing_commit',lambda pins:None)
    with pytest.raises(ValueError,match='frozen producer'):parent.prepare(tmp_path,nf=32,cpu_limit=5.)
    assert not (tmp_path/'parent.json').exists()
    # A simulated immutable commit permits bounded unit evidence only.
    monkeypatch.setattr(parent,'producing_commit',lambda pins:'a'*40)
    monkeypatch.setattr(parent.provenance,'resolve_pinned_source_bytes',lambda *args,**kwargs:b'unit source')
    parent.prepare(tmp_path,nf=32,population=1,cpu_limit=5.)
    assert parent.check(tmp_path)['numerical_replay']
    with pytest.raises(FileExistsError):parent.prepare(tmp_path,nf=32,population=1,cpu_limit=5.)
    assert before=={path:parent._sha_file(path) for path in before}


def test_pure_preview_and_shared_budget_common_offset_pilot(tmp_path,monkeypatch):
    monkeypatch.setattr(parent,'parent_continuum',lambda *args,**kwargs:pytest.fail('preview solved spatial ODE'))
    assert parent.preview()['constructor_called'] is False
    monkeypatch.setattr(parent,'producing_commit',lambda pins:'a'*40)
    monkeypatch.setattr(parent.provenance,'resolve_pinned_source_bytes',lambda *args,**kwargs:b'unit source')
    seed_calls=[];case_calls=[]
    def seed(nf,*,cpu_limit,profile=None):
        seed_calls.append((nf,cpu_limit));return {'k':float(nf)/100.,'kmin':.1}
    def case(path,**kwargs):
        case_calls.append(kwargs);path.mkdir(parents=True)
        record={'nf':kwargs['nf'],'k':kwargs['k_override'],'kmin':.1,'diagnosis':{'anchor':kwargs['k_override']+.05},
                'converged':False,'blocker':'unit bounded checkpoint','CPU_seconds':0.}
        (path/'parent.json').write_text(json.dumps(record));return record
    monkeypatch.setattr(parent,'common_k',seed);monkeypatch.setattr(parent,'prepare',case)
    ledger=parent.pilot(tmp_path,cpu_limit=2.)
    assert [value['nf'] for value in case_calls]==[64,128]
    assert all(value['k_override']==1.28 for value in case_calls)
    assert case_calls[0]['cpu_limit']+case_calls[1]['cpu_limit']<=3.9  # second receives the measured remaining budget.
    assert ledger['common_k']==1.28 and ledger['aggregate_CPU_seconds']<=2.
    assert ledger['complete'] and all(row['kinetic_anchor']!=row['k'] for row in ledger['records'])
    with pytest.raises(FileExistsError):parent.pilot(tmp_path,cpu_limit=2.)


def test_scalar_anchor_cannot_be_truncated_by_reported_rank_cutoff():
    # Reproduce the immutable pilot's singular scales without reading/writing it.
    target=.000135234609664;actual=.000154624662736;gap=actual-target
    jac=np.diag(np.r_[np.linspace(2922.,1.,127),.16]);rows=np.ones(128)
    rows[-1]=.16/8.235e-8
    residual=np.zeros(128);residual[-1]=gap
    scaled=jac/rows[:,None];old,oldrank,singular=parent._svd_step(scaled,residual/rows)
    assert oldrank==127 and singular[-1]<1e-10*singular[0]
    assert (jac@old+residual)[-1]==pytest.approx(gap,abs=1e-16)
    hard,rank,_=parent._hard_anchor_step(jac,residual,np.eye(128),rows)
    assert rank==128 and abs((jac@hard+residual)[-1])<1e-18
    np.testing.assert_allclose((jac@hard+residual)[:-1],0.,atol=1e-14)


def test_hard_anchor_retraction_convergence_and_sign_margin_preserve_source(setup):
    with threadpool_limits(limits=1),parent.backend.fft_thread_limit(1):
        problem,seed,theta=problem_for(setup)
        before=(problem.phi0.tobytes(),problem.phi1.tobytes(),problem.pair.weights.tobytes())
        target=problem.anchor;theta[:problem.count]*=np.sqrt(.000154624662736/.000135234609664)
        assert not problem.anchor_status(theta)['satisfied']
        restored=problem.restore_anchor(theta)
        assert problem.anchor_status(restored)['satisfied'] and problem.anchor==target
        residual=problem.residual(theta);jac=problem.jacobian(theta)
        step,rank,_=parent._hard_anchor_step(jac,residual,problem.whitener,problem.row_scales(theta))
        assert abs(jac[-1]@problem.whitener@step+residual[-1])<1e-13*target
        corrected,report=parent.correct_momenta(problem,theta,parent._cpu_time()+2.,max_accepted=1)
        assert report['hard_anchor']['satisfied'] and report['hard_anchor']['relative_error']<=parent.ANCHOR_RELATIVE_TOL
        assert report['initial_anchor']['relative_error']>.1
        assert report['anchor_never_subject_to_relative_SVD_cutoff']
        assert report['final_sign_margin']>=report['required_nonzero_sign_margin']>0
        assert before==(problem.phi0.tobytes(),problem.phi1.tobytes(),problem.pair.weights.tobytes())
        negative=restored.copy();negative[:problem.count]*=-1
        stopped,blocked=parent.correct_momenta(problem,negative,parent._cpu_time()+2.,max_accepted=1)
        assert not blocked['converged'] and blocked['accepted']==0
        assert 'sign margin' in blocked['blocker']


def test_transition_option_is_forwarded_without_new_constructor_run(monkeypatch,tmp_path):
    import derive_nsc_discovery_parent as cli
    calls=[]
    monkeypatch.setattr(parent,'pilot',lambda *args,**kwargs:calls.append(kwargs) or {'unit':True})
    cli.main(['--pilot','--transition-end','1.4','--output',str(parent.OUTPUT/'unit-unwritten')])
    assert calls[0]['profile']=={'initial_radius':'collar','transition_end':1.4}
    assert parent._profile(calls[0]['profile'])['transition_end']<np.pi/2


def test_magnetic_initial_radius_preserves_original_source_weights_Q_and_observer(setup):
    with threadpool_limits(limits=1),parent.backend.fft_thread_limit(1):
        table,pair,r,Q,columns=setup
        magnetic=dict(table,profile=dict(table['profile'],initial_radius='magnetic'))
        mr,mQ=parent._band_geometry(pair.grid,magnetic)
        copied=parent._source_columns(pair.grid,magnetic)
        assert mQ.tobytes()==Q.tobytes()
        assert np.array_equal(mr,np.full(pair.grid.ng,np.sqrt(table['constants']['mag']/table['constants']['g'])))
        for key in ('phi0','phi1','source_columns','reference_columns','raw_norms','probabilities'):
            assert copied[key].tobytes()==columns[key].tobytes()
        gram=columns['source_columns'].conj().T@columns['source_columns']
        for population in parent.POPULATIONS:
            first,_=parent.population_weights(columns['probabilities'],columns['raw_norms'],gram,population)
            second,_=parent.population_weights(copied['probabilities'],copied['raw_norms'],gram,population)
            assert first.tobytes()==second.tobytes()
        targets=parent._geometry_targets(pair.grid,magnetic)
        np.testing.assert_array_equal(targets['radius'],parent._sample_circle(pair.grid,table,'radius'))
        changes=parent._geometry_changes(pair.grid,mr,mQ,targets)
        assert changes['protected_core_radius_relative_max']>0 and changes['protected_core_exactly_preserved'] is False


def test_pilot_explicit_k_and_magnetic_profile_forwarding(tmp_path,monkeypatch):
    import derive_nsc_discovery_parent as cli
    calls=[]
    monkeypatch.setattr(parent,'pilot',lambda *args,**kwargs:calls.append(kwargs) or {'unit':True})
    cli.main(['--pilot','--initial-radius','magnetic','--transition-end','1.4','--k','.002621213673327129',
              '--output',str(parent.OUTPUT/'unit-magnetic-unwritten')])
    assert calls[0]['profile']=={'initial_radius':'magnetic','transition_end':1.4}
    assert calls[0]['k_override']==.002621213673327129


def test_physical_fallback_budget_uses_original_candidate_clock_and_lengths(setup):
    table,pair,r,Q,col=setup;grid=pair.grid
    original=dict(table,sigma=np.array([0.,.5,1.,3.]),radius=np.ones(4),Q=np.ones(4),
                  exterior_radius=1.,profile=dict(table['profile'],exterior_Q=1.))
    targets=parent._geometry_targets(grid,original)
    flat=np.ones(grid.ng)
    baseline=parent._geometry_changes(grid,flat,flat,targets)
    assert baseline['within_budget'] and baseline['relative_budget']==.05
    assert baseline['proper_lengths']==pytest.approx({'child':1.,'collar':2.,'parent':6.,'whole':8.})
    # 6% against the original must fail even if only 3.9% from a changed band seed.
    assert parent._geometry_changes(grid,flat,flat*1.02,targets)['within_budget']
    assert not parent._geometry_changes(grid,flat,flat*1.06,targets)['within_budget']
    compounded=parent._geometry_changes(grid,flat*1.03,flat*1.03,targets)
    assert compounded['radius_relative_max']<.05 and compounded['Q_relative_max']<.05
    assert compounded['child_clock_relative_max']>.05 and not compounded['within_budget']
    huge=parent._geometry_changes(grid,flat,flat*17.,targets)
    assert not huge['within_budget'] and huge['child_clock_relative_max']>15
    allowed=parent._geometry_changes(grid,flat*.98,flat*1.04,targets)
    assert allowed['within_budget'] and allowed['protected_core_N_relative_max']>0
    assert allowed['protected_core_exactly_preserved'] is False


def test_fallback_refuses_large_geometry_moves_and_preserves_source(setup,monkeypatch):
    with threadpool_limits(limits=1),parent.backend.fft_thread_limit(1):
        table,pair,r,Q,col=setup;problem,seed,theta=problem_for(setup)
        targets=parent._geometry_targets(pair.grid,table)
        before=(col['phi0'].tobytes(),col['phi1'].tobytes(),pair.weights.tobytes())
        monkeypatch.setattr(parent,'_svd_step',lambda jac,res:(np.array([1e4,0.,0.,0.]),4,np.ones(4)))
        best,report=parent._fallback(pair,r,Q,col['phi0'],col['phi1'],1,seed['g'],seed['suggested_k'],.05,
                                     parent._cpu_time()+2.,problem.anchor,targets)
        assert report['steps']==0 and report['refused_move_count']>=parent.MAX_BACKTRACKS
        np.testing.assert_allclose(best[0],r,atol=2e-15,rtol=0.)
        np.testing.assert_allclose(best[1],Q,atol=2e-15,rtol=0.)
        assert report['actual_parameter_changes']==[0.,0.,0.,0.]
        assert any(row['kind']=='backtrack' and not row['changes']['within_budget'] for row in report['refused_moves'])
        assert before==(col['phi0'].tobytes(),col['phi1'].tobytes(),pair.weights.tobytes())
