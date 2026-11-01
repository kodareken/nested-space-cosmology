"""Same-state directional handoff and geometry controls; no PDE evolution."""
import json
import math
from pathlib import Path
import subprocess

import numpy as np
import pytest
from scipy.integrate import DOP853
from threadpoolctl import threadpool_limits

import evolve_nsc_atlas_valve_directional as driver


@pytest.fixture
def saved(tmp_path,monkeypatch):
    monkeypatch.setattr(driver,'ROOT',tmp_path)
    with threadpool_limits(limits=1):s=driver.DirectionalEvolution(16)
    folder=tmp_path/'.build/parent';folder.mkdir(parents=True)
    names=[driver.RUNNER,'src/recursive_horizons/nsc_atlas_valve_spherical.py',
        'src/recursive_horizons/nsc_atlas_valve_neck.py','src/recursive_horizons/nsc_atlas_passive_dop853.py']
    for i,name in enumerate(names):
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(f'# source{i}\n'.encode())
    subprocess.run(['git','init','-q'],cwd=tmp_path,check=True,capture_output=True)
    subprocess.run(['git','-c','core.autocrlf=false','add','.'],cwd=tmp_path,check=True,capture_output=True)
    subprocess.run(['git','-c','user.name=NSC test','-c','user.email=nsc-test@example.invalid','commit','-qm','Source fixture'],cwd=tmp_path,check=True,capture_output=True)
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=tmp_path,text=True).strip()
    pins={name:driver.sha(tmp_path/name) for name in names}
    (tmp_path/driver.RUNNER).write_bytes(b'# current runner has resume support\n')
    initial=tmp_path/'initial.json';initial.write_bytes(b'{"synthetic":true}\n')
    state=s.initial.copy();state[2]+=.001*s.u;state[7]+=.002
    cache=folder/'sparse.npz';np.savez_compressed(cache,times=[0.,2.],states=np.array([s.initial,state]))
    def cap(energy,material,work,residual):return dict(energy=energy,integrated_material_received=material,
        integrated_pressure_work=work,energy_ledger_residual=residual)
    rows=[dict(proper_time=0.,sampled_admitted=True,minimum={'u':0.},
        fixed_receiving_cap=cap(100.,0.,0.,0.),fixed_donor_cap=cap(200.,0.,0.,0.)),
        dict(proper_time=2.,sampled_admitted=True,minimum={'u':.012},
        fixed_receiving_cap=cap(101.01,.4,.6,.01),fixed_donor_cap=cap(199.82,-.4,.2,.02))]
    record=dict(schema='nsc-atlas-valve-directional-evolution-v1',status='endpoint_reached',stop_reason=None,
        modes=16,reflection_projection=False,source_bound_to_commit=True,source_unchanged=True,
        producer_commit=commit,source_sha256=pins,initial_record='initial.json',initial_record_sha256=driver.sha(initial),
        state_cache=cache.relative_to(tmp_path).as_posix(),state_cache_sha256=driver.sha(cache),
        snapshots=rows,last_admitted_proper=2.,actual_last_accepted_proper=2.)
    path=folder/'run.json';path.write_text(json.dumps(record))
    return s,path,record,state,[{'energy':100.},{'energy':200.}]


def test_exact_all_field_handoff_and_three_cumulative_integrals(saved):
    s,path,record,state,baseline=saved
    c=driver.load_resume(s,path,initial_record=Path('initial.json'),baseline_caps=baseline)
    np.testing.assert_array_equal(c['state'],state)
    assert c['proper_time']==2. and c['previous_minimum_u']==.012
    np.testing.assert_array_equal(c['integrals'],[.4,.6,.2])
    assert c['baselines']==[100.,200.]
    assert c['provenance']['parent_record_sha256']==driver.sha(path)
    assert c['provenance']['state_cache_sha256']==record['state_cache_sha256']


@pytest.mark.parametrize('change,match',[
    ({'status':'running'},'terminal endpoint/budget'),
    ({'status':'incomplete','stop_reason':'sampled ADM admission ended; no physical endpoint inferred'},'admission-failed'),
    ({'modes':24},'same-N'),
    ({'source_unchanged':False},'same-N'),
])
def test_invalid_terminal_or_grid_cannot_resume(saved,change,match):
    s,path,record,_,baseline=saved;record.update(change);path.write_text(json.dumps(record))
    with pytest.raises(ValueError,match=match):driver.load_resume(s,path,initial_record=Path('initial.json'),baseline_caps=baseline)


def test_current_physical_rhs_and_raw_cache_must_match(saved,monkeypatch):
    s,path,_,_,baseline=saved;actual=driver.sha
    monkeypatch.setattr(driver,'sha',lambda p:'0'*64 if Path(p).name=='nsc_atlas_valve_spherical.py' else actual(p))
    with pytest.raises(ValueError,match='current physical package'):
        driver.load_resume(s,path,initial_record=Path('initial.json'),baseline_caps=baseline)
    monkeypatch.setattr(driver,'sha',lambda p:'0'*64 if Path(p).name=='sparse.npz' else actual(p))
    with pytest.raises(ValueError,match='sparse cache identity'):
        driver.load_resume(s,path,initial_record=Path('initial.json'),baseline_caps=baseline)


def test_baselines_and_anchor_survive_a_lost_minimum(saved):
    s,path,record,_,baseline=saved
    record['cap_energy_baselines']=[100.,200.]
    record['snapshots'][-1]['minimum']=None;record['snapshots'][-1]['minimum_anchor_u']=-.075
    path.write_text(json.dumps(record))
    c=driver.load_resume(s,path,initial_record=Path('initial.json'),baseline_caps=baseline)
    assert c['previous_minimum_u']==-.075 and c['baselines']==[100.,200.]
    with pytest.raises(ValueError,match='baselines changed'):
        driver.load_resume(s,path,initial_record=Path('initial.json'),baseline_caps=[{'energy':101.},{'energy':200.}])


def test_three_passive_ledgers_start_at_original_proper_time(saved):
    s,path,_,state,baseline=saved
    c=driver.load_resume(s,path,initial_record=Path('initial.json'),baseline_caps=baseline)
    t0=c['proper_time'];end=2.1;calls=[]
    ledgers=[driver.PassiveDOP853Quadrature(initial=v) for v in c['integrals']]
    def rhs(t,y):
        calls.append(t)
        for ledger,rate in zip(ledgers,[1+t,2*t,3*t*t]):ledger.record(rate)
        return np.zeros_like(y)
    solver=DOP853(rhs,t0,state.ravel(),end,first_step=.002,max_step=.02,rtol=1e-10,atol=1e-12)
    for ledger in ledgers:ledger.initialize(solver)
    while solver.status=='running':
        solver.step()
        for ledger in ledgers:ledger.accepted_step(solver)
    assert calls[0]==2.
    np.testing.assert_array_equal(solver.y,state.ravel())
    expected=c['integrals']+np.array([end-t0+(end**2-t0**2)/2,end**2-t0**2,end**3-t0**3])
    np.testing.assert_allclose([v.value for v in ledgers],expected,atol=3e-14,rtol=0)


def test_minimum_second_derivative_comes_from_actual_spatial_metric(saved):
    s,_,_,_,_=saved;measure=driver.Measurements(s);state=s.initial
    result=measure.minimum(state);assert result is not None
    u=np.array([result['u']]);v=1-u*u;z=s.pieces(state,u);_,yu,yuu=s.jets(state,u)
    bu,cu=yu[:2];qu=bu-2*u*z['y'][1]+v*cu
    Rchichi=-np.sqrt(v)*(z['b']+3*u*bu-v*yuu[0])
    qchi=-np.sqrt(v)*qu
    direct=(Rchichi-qchi/z['q']*z['rr'])/z['q']**2
    assert result['R_ll']==pytest.approx(float(direct[0]),rel=2e-13,abs=2e-14)
    assert result['R_ll']>0


def test_band_check_keeps_last_anchor_when_minimum_disappears(monkeypatch):
    measure=driver.Measurements.__new__(driver.Measurements);seen=[]
    monkeypatch.setattr(measure,'minimum',lambda state,previous:None)
    def trapped(state,anchor):seen.append(anchor);return {'kept':True}
    monkeypatch.setattr(measure,'trapped',trapped)
    minimum,anchor,band=measure.tracked_geometry(None,-.075,2.2)
    assert minimum is None and anchor==-.075 and seen==[-.075]
    assert band['anchor_u']==-.075 and not band['anchor_is_current_minimum']
