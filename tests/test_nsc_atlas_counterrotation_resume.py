"""Portable checkpoint identity, admission, time and cumulative-work controls."""
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import pytest
from scipy.integrate import DOP853

from recursive_horizons import nsc_atlas_counterrotation_resume as resume
from recursive_horizons.nsc_atlas_passive_dop853 import PassiveDOP853Quadrature


@pytest.fixture(scope='module')
def saved_run(tmp_path_factory):
    root=tmp_path_factory.mktemp('counterrotation-resume-repository')
    package=root/'src/recursive_horizons';package.mkdir(parents=True)
    (root/'scripts').mkdir();folder=root/'.build/fine';folder.mkdir(parents=True)
    source_names=[resume.PRODUCER,
        'src/recursive_horizons/nsc_atlas_parent_gh.py',
        'src/recursive_horizons/nsc_atlas_gh_pointwise.py',
        'src/recursive_horizons/nsc_atlas_gh_history.py',
        'src/recursive_horizons/nsc_atlas_s3_hopf.py']
    for i,name in enumerate(source_names):(root/name).write_bytes(f'# synthetic source {i}\n'.encode())
    subprocess.run(['git','init','-q'],cwd=root,check=True,capture_output=True)
    subprocess.run(['git','-c','core.autocrlf=false','add','.'],cwd=root,check=True,capture_output=True)
    subprocess.run(['git','-c','user.name=NSC test','-c','user.email=nsc-test@example.invalid',
        'commit','-qm','Synthetic checkpoint decoder fixture'],cwd=root,check=True,capture_output=True)
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
    pins={name:resume.sha(root/name) for name in source_names}
    # The historical producer is verified in Git. Its present text may acquire
    # the resume feature; all old package decoder/RHS bytes must remain fixed.
    (root/resume.PRODUCER).write_bytes(b'# synthetic current producer with resume support\n')
    initial=root/'initial.json';initial.write_bytes(b'{"synthetic":true}\n')
    meta=dict(degree=4,epsilon=1.,gauge='lapse_controlled',backend='hopf',quadrature_order=[5,14],
        kick_amplitude=1.,kick_source=2,source_sha256=pins,producer_commit=commit,
        preparation_id='one-mode-counterrotation-v1')
    times=np.array([0.,.101,.205]);proper=np.array([0.,.1,.2])
    states=np.arange(54,dtype=float).reshape(3,18)/100
    states[:,-10:-6]=[0.,2**-.5,0.,2**-.5]
    states[:,-6:-2]=[1.01,.02,-.03,.04]
    states[:,-2]=proper;states[:,-1]=[.2,.4,.6]
    cache=folder/'sparse.npz'
    np.savez_compressed(cache,times=times,states=states,metadata=json.dumps(meta))
    dense=folder/'dense-0000.npz'
    np.savez_compressed(dense,bounds=np.array([[0.,.15],[.15,.23]]),old=states[:2],
        coefficients=np.zeros((2,7,18)),metadata=json.dumps(meta))
    tail=folder/'accepted-endpoint.npz'
    np.savez_compressed(tail,state=states[-1]+.001,coordinate_t=.23,metadata=json.dumps(meta))
    snapshots=[dict(coordinate_t=float(t),observer_proper_duration=float(p),sampled_admitted=True,
        normal_energy_minus_fixed_vacuum=12.+w+.01*(i>0),integrated_excitation_work=w)
        for i,(t,p,w) in enumerate(zip(times,proper,[0.,.4,.7]))]
    record=dict(schema='nsc-atlas-parent-gh-v3',status='incomplete',
        stop_reason='declared formation run budget exhausted',degree=4,epsilon=1.,
        prescribed_GH_gauge='lapse_controlled',spatial_transform='hopf',
        producer_commit=commit,source_sha256=pins,source_bound_to_commit=True,
        producer_source_unchanged_during_run=True,initial_record_unchanged=True,passive_work_valid=True,
        initial_record='initial.json',initial_record_sha256=resume.sha(initial),
        model_inputs=dict(H_reference=.3,vacuum=.27,scalar_mass_squared=400.,volume_kinetic=84.,
            sstar=-.03,kick_amplitude=1.,kick_source=2),snapshots=snapshots,
        state_cache=cache.relative_to(root).as_posix(),state_cache_sha256=resume.sha(cache),
        dense_history_chunks=[dict(path=dense.relative_to(root).as_posix(),sha256=resume.sha(dense),
            coordinate_bounds=[0.,.23],intervals=2)],last_admitted_coordinate_t=.205,
        actual_last_accepted_coordinate_t=.23,actual_last_accepted_proper=.224,
        accepted_endpoint_cache=dict(path=tail.relative_to(root).as_posix(),sha256=resume.sha(tail),
            admitted_at_this_exact_time=False))
    path=folder/'run.json';path.write_text(json.dumps(record))
    return root,path,record,states


def test_full_state_observer_and_work_are_exact_admitted_snapshot(saved_run):
    root,path,record,states=saved_run
    c=resume.load_resume_checkpoint(path,repository_root=root,degree=4)
    np.testing.assert_array_equal(c.state,states[-1])
    assert c.coordinate_t==.205 and c.proper==.2
    assert c.cumulative_work==.7 and c.baseline_excitation_energy==12.
    assert c.provenance['state_array_sha256']==hashlib.sha256(states[-1].tobytes()).hexdigest()
    assert c.provenance['physical_history_prefix'][-1]['usable_coordinate_bounds']==[0.,.205]
    assert c.provenance['preserved_unpromoted_parent_tail']['coordinate_bounds']==[.205,.23]
    assert record['accepted_endpoint_cache']['admitted_at_this_exact_time'] is False


def test_live_admission_failure_and_wrong_grid_are_rejected(saved_run):
    root,path,_,_=saved_run
    with pytest.raises(ValueError,match='terminal budget'):resume._terminal_allowed({'status':'initializing'})
    with pytest.raises(ValueError,match='admission failures'):
        resume._terminal_allowed({'status':'incomplete','stop_reason':'sampled constraint/gauge/observer admission ended; no physical endpoint is inferred'})
    with pytest.raises(ValueError,match='same-grid'):
        resume.load_resume_checkpoint(path,repository_root=root,degree=6)


@pytest.mark.parametrize('broken_name,expected',[
    ('nsc_atlas_parent_gh.py','current package/source differs'),
    ('sparse.npz','sparse cache hash mismatch'),
    ('dense-0000.npz','parent dense history hash failed'),
])
def test_identity_failures_do_not_silently_replay(saved_run,monkeypatch,broken_name,expected):
    root,path,_,_=saved_run;actual=resume.sha
    monkeypatch.setattr(resume,'sha',lambda p:'0'*64 if Path(p).name==broken_name else actual(p))
    with pytest.raises(ValueError,match=expected):
        resume.load_resume_checkpoint(path,repository_root=root,degree=4)


def test_original_baseline_survives_repeated_resume(saved_run):
    root,path,record,_=saved_run
    changed=dict(record,excitation_work_origin_energy=9.)
    other=path.with_name('with-original-baseline.json');other.write_text(json.dumps(changed))
    c=resume.load_resume_checkpoint(other,repository_root=root,degree=4)
    assert c.baseline_excitation_energy==9.


def test_targets_are_after_handoff():
    assert resume.remaining_proper_targets(2.,1.1999999999999997)==[1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.]
    with pytest.raises(ValueError,match='later'):resume.remaining_proper_targets(1.2,1.2)


def test_coordinate_time_and_cumulative_passive_work_continue_without_reset():
    t0=1.2;t1=1.5;normalization=2.;work0=5.;baseline=3.;old_residual=.25
    ledger=PassiveDOP853Quadrature(initial=work0/normalization);calls=[]
    def rhs(t,y):
        calls.append(t);ledger.record(3*t*t/normalization)
        return np.array([t])
    solver=DOP853(rhs,t0,np.array([10.]),t1,rtol=1e-10,atol=1e-12,max_step=.03)
    ledger.initialize(solver)
    while solver.status=='running':solver.step();ledger.accepted_step(solver)
    assert calls[0]==t0
    assert solver.y[0]==pytest.approx(10.+.5*(t1*t1-t0*t0),abs=2e-12)
    expected=work0+t1**3-t0**3;work=ledger.value*normalization
    assert work==pytest.approx(expected,abs=2e-12)
    energy=baseline+old_residual+expected
    assert energy-baseline-work==pytest.approx(old_residual,abs=2e-12)
