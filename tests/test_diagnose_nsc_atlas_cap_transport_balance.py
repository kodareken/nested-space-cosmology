"""Synthetic fixed-case routing, numerical gates, provenance and replay receipts."""
import copy
import json
from pathlib import Path

import pytest

import diagnose_nsc_atlas_cap_transport_balance as tool
import test_nsc_atlas_cap_mass_balance as older
from control_budget import ControlBudget

CONTROLS = ControlBudget(maximum_evaluations=250_000, cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('terminal transport recorder controls:',CONTROLS.budget.evaluations,
          'charged values;',CONTROLS.active_cpu_seconds,'active CPU seconds')


def manufactured_result(old, settings):
    channels = dict.fromkeys(tool.CHANNELS,0.)
    channels['field_Pi'] = old['B_mass']-.02
    channels['interfaces'] = .02
    return dict(channels=channels, total=old['B_mass'], native_balance=old['B_mass'],
        actual_native_balance=old['actual_B_mass'], balance_defect=0., reference_mass_defect=0.,
        reference_rate_defect=0., native_wall=dict(rho=old['radial_density'],J=old['flux_left']),
        maps=dict(error=1e-7,native_mass=.2,reference_mass=.2-1e-7,reference_rate=.3,
            telescope_defects=[0.,0.],telescope_roundoff_bounds=[1e-12,1e-12],
            native_tangent_defect=0.,native_tangent_allowance=1e-12,
            reference_coefficient_calls=30,requested_rtol=settings['rtol'],requested_atol=settings['atol']),
        reference=dict(wall_fields=dict(a=.2),jump_sum_defect=0.,delta=.001,minimum_sampled_A=.7,
            sampled_A_positive=True,jumps=[dict(r=.2)],stages=[dict(piece=1),dict(piece=2)],
            coefficient_calls=40,rhs_calls=40),evaluations=70)


def fixtures(monkeypatch,budget):
    data = older.fixtures(budget)
    seen = older.install(monkeypatch)
    prior = older.calc.analyze(*data,budget)
    for case,original in zip(prior['cases'],data[0]['cases']):
        state = original['snapshots'][128]['state']['cap']
        R,T=.49,.08
        case['orders']['8']['tangent']['N'] = R*R*T+sum(state['Ip'])*sum(state['Iq'])+T**3+8e-5*R
    for key in seen:
        seen[key].clear()
    for name in ('background_cap_metric','background_cap_rhs'):
        monkeypatch.setattr(tool,name,getattr(older.calc,name))
    def tangent(metric,direction,*,budget):
        row = older.calc.mass_tangent(metric,direction,budget=budget)
        row['N'] = float(metric['n_faces'][-1])
        return row
    balance_calls = []
    def balance(metric,direction,tangent,*,budget,rtol,atol):
        assert metric['quadrature_order']==8
        balance_calls.append((metric['alpha'],copy.deepcopy(direction),rtol,atol))
        row = prior['cases'][int(metric['alpha']>0)]['orders']['8']
        budget.check(70)
        return manufactured_result(row,dict(rtol=rtol,atol=atol))
    monkeypatch.setattr(tool,'mass_tangent',tangent)
    monkeypatch.setattr(tool,'cap_transport_balance',balance)
    return (*data[:3],dict(data=prior)),seen,balance_calls


def test_two_states_two_tolerances_fixed_inputs_and_accounting(monkeypatch,budget):
    data,seen,calls = fixtures(monkeypatch,budget)
    original = copy.deepcopy(data); start = budget.evaluations
    result = tool.analyze(*data,budget)
    assert data==original
    assert [len(v) for v in seen.values()]==[2,2,2] and len(calls)==4
    assert [(r[2],r[3]) for r in calls]==[(s['rtol'],s['atol']) for s in tool.SETTINGS]*2
    assert result['evaluated_points']==budget.evaluations-start==sum(result['execution']['charges'].values())
    for key,value in tool.EXPECTED_CALLS.items():
        assert result['execution'][key]==value
    assert all(all(v==0 for v in c['fine_minus_coarse'].values()) for c in result['summary']['cases'])


@pytest.mark.parametrize('kind',['native','reference_mass','reference_rate','jump','balance','sensitivity'])
def test_each_material_failure_stops_without_retry(monkeypatch,budget,kind):
    data,seen,calls = fixtures(monkeypatch,budget)
    original = tool.cap_transport_balance
    def fail(*a,**kw):
        row = original(*a,**kw)
        if kind=='native': row['native_balance']+=.01
        elif kind=='reference_mass': row['reference_mass_defect']=2e-8
        elif kind=='reference_rate': row['reference_rate_defect']=3e-7
        elif kind=='jump': row['reference']['jump_sum_defect']=2e-7
        elif kind=='balance': row['balance_defect']=2e-6
        elif len(calls)==2: row['channels']['field_f']=6e-7
        return row
    monkeypatch.setattr(tool,'cap_transport_balance',fail)
    with pytest.raises(ValueError) as caught:
        tool.analyze(*data,budget)
    assert [len(v) for v in seen.values()]==[1,1,1]
    assert len(calls)==(2 if kind=='sensitivity' else 1)
    receipt = caught.value.transport_failure
    assert receipt['counts']['balance_calls']==len(calls)
    if kind!='native': assert receipt['accuracy']


def test_wall_effect_and_telescope_are_separate_gates():
    old=dict(B_mass=.2,actual_B_mass=.01,radial_density=.3,flux_left=.4)
    result=manufactured_result(old,tool.SETTINGS[0])
    result['reference_mass_defect']=9e-9
    with pytest.raises(ValueError) as caught:
        tool.run_summary(result,dict(R=100.),old)
    assert 'wall_density_effect' in caught.value.accuracy_failure['failures']
    result['reference_mass_defect']=0.
    result['maps']['telescope_defects'][1]=2e-12
    with pytest.raises(ValueError) as caught:
        tool.run_summary(result,dict(R=1.),old)
    assert 'rate_telescope' in caught.value.accuracy_failure['failures']


def synthetic(monkeypatch,root,budget):
    for name in tool.OWN_SOURCES:
        path=root/name; path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(b'synthetic frozen producer')
    calls=[]
    def auth(root,shared,sources):
        assert shared is budget
        calls.append('auth')
        return (1,2,3,4),dict(fixture=True)
    def analyze(*args):
        assert args==(1,2,3,4,budget)
        calls.append('analysis'); budget.check(7)
        return dict(summary=dict(value=.2,gates=tool.GATES),details=[dict(coefficient_calls=9)],
            evaluated_points=7,execution=dict(tool.EXPECTED_CALLS,
            charges=dict(input_and_profile=1,metrics=1,rhs=1,tangents=1,balances=3)))
    monkeypatch.setattr(tool,'authenticate',auth)
    monkeypatch.setattr(tool,'analyze',analyze)
    monkeypatch.setattr(tool,'validate_saved_gates',lambda *a:None)
    return calls


def test_create_reconstruct_check_no_overwrite_and_adaptive_counts(monkeypatch,tmp_path,budget):
    calls=synthetic(monkeypatch,tmp_path,budget); path=tmp_path/'.build/results/audit.json'
    saved=tool.operation(path,root=tmp_path,budget=budget)
    assert calls==['auth','analysis','analysis'] and saved['execution']['evaluated_points']==14
    saved['analysis']['details'][0]['coefficient_calls']+=1
    saved['analysis']['summary']['value']+=1e-10
    path.write_bytes(tool.encoded(saved))
    assert tool.operation(path,checking=True,root=tmp_path,budget=budget)['evaluated_points']==7
    raw=path.read_bytes()
    with pytest.raises(FileExistsError): tool.operation(path,root=tmp_path,budget=budget)
    assert path.read_bytes()==raw


@pytest.mark.parametrize('kind',['science','gate','receipt','producer','counter'])
def test_saved_drift_has_bounded_failure_receipt(monkeypatch,tmp_path,budget,kind):
    calls=synthetic(monkeypatch,tmp_path,budget); path=tmp_path/'.build/results/audit.json'
    saved=tool.operation(path,root=tmp_path,budget=budget)
    if kind=='science': saved['analysis']['summary']['value']+=.001
    elif kind=='gate': saved['analysis']['summary']['gates']=dict(tool.GATES,reference_mass=2e-8)
    elif kind=='receipt': saved['execution']['hash_authentication_passes']=True
    elif kind=='counter': saved['analysis']['details'][0]['coefficient_calls']=True
    else: (tmp_path/tool.OWN_SOURCES[-1]).write_bytes(b'changed producer')
    path.write_bytes(tool.encoded(saved)); raw=path.read_bytes()
    with pytest.raises(ValueError) as caught:
        tool.operation(path,checking=True,root=tmp_path,budget=budget)
    receipt=json.loads(Path(caught.value.failure_receipt).read_bytes())
    assert receipt['retries']==0 and receipt['phase']=='check reconstruction'
    assert path.read_bytes()==raw


def test_mass_summary_pin_and_one_inherited_hash_pass(monkeypatch,tmp_path,budget):
    from test_diagnose_nsc_atlas_cap_mass_direction import pinned_fixture
    pinned_fixture(tmp_path,monkeypatch,budget)
    for name in tool.OWN_SOURCES:
        path=tmp_path/name; path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(b'synthetic frozen producer')
    # The older direction summary and its mass-ledger successor are byte pins.
    prior=tmp_path/tool.previous.PRIOR; prior.parent.mkdir(parents=True,exist_ok=True)
    raw=tool.encoded(dict(data=dict(fixture=True))); prior.write_bytes(raw)
    monkeypatch.setattr(tool.previous,'PRIOR_SHA256',tool.digest(raw))
    path=tmp_path/tool.PRIOR
    raw=tool.encoded(dict(data=dict(fixture=True),producer_bindings={n:dict(sha256=tool.digest((tmp_path/n).read_bytes()))
                                                        for n in tool.previous.OWN_SOURCES}))
    path.write_bytes(raw); monkeypatch.setattr(tool,'PRIOR_SHA256',tool.digest(raw))
    original=tool.previous.authenticate; calls=[]
    def auth(*a): calls.append(1); return original(*a)
    def forbidden(*a,**kw): raise AssertionError('old scientific replay')
    monkeypatch.setattr(tool.previous,'authenticate',auth)
    monkeypatch.setattr(tool.previous,'analysis',forbidden)
    start=budget.evaluations
    data,identity=tool.authenticate(tmp_path,budget,tool.bindings(tmp_path))
    assert data[3]['data']['fixture'] and calls==[1] and budget.evaluations==start
    assert identity['previous_analysis_replays']==0
    path.write_bytes(b'changed summary')
    with pytest.raises(ValueError,match='summary mismatch'):
        tool.authenticate(tmp_path,budget,tool.bindings(tmp_path))


def test_allowance_cannot_be_expanded_or_clock_reset():
    for kwargs in (dict(maximum_evaluations=48_000_001,cpu_seconds=300.),
                   dict(maximum_evaluations=1,cpu_seconds=301.),
                   dict(maximum_evaluations=True,cpu_seconds=10.)):
        with pytest.raises(ValueError): tool.shared_budget(tool.Budget(**kwargs))


def test_saved_gates_rechecked_even_within_replay_allowance(monkeypatch,budget):
    data,_,_ = fixtures(monkeypatch,budget)
    result = tool.analyze(*data,budget)
    tool.validate_saved_gates(result,data[3])
    changed=copy.deepcopy(result)
    changed['details'][0]['runs'][0]['result']['reference_mass_defect']=1.1e-8
    with pytest.raises(ValueError,match='accuracy gate'):
        tool.validate_saved_gates(changed,data[3])
    changed=copy.deepcopy(result)
    changed['summary']['cases'][0]['runs'][0]['defects']['balance']=1e-10
    with pytest.raises(ValueError,match='scalar summary'):
        tool.validate_saved_gates(changed,data[3])
    assert tool.failure_values(dict(value=float('nan')))==dict(value='nan')


@pytest.mark.parametrize('name',['results/summaries/full.json','.build/elsewhere/full.json','outside.json'])
def test_full_output_cannot_be_written_to_tracked_path(tmp_path,budget,name):
    with pytest.raises(ValueError,match='ignored'):
        tool.operation(tmp_path/name,root=tmp_path,budget=budget)
    assert not (tmp_path/name).exists()


def test_summary_selects_only_scalars_and_binds_regeneration(monkeypatch,tmp_path,budget):
    synthetic(monkeypatch,tmp_path,budget)
    for name in (tool.PRIOR,tool.previous.PRIOR,tool.previous.previous.CAP_INPUT):
        p=tmp_path/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b'{}')
    raw=tmp_path/'.build/results/audit.json'
    payload=tool.operation(raw,root=tmp_path,budget=budget)
    summary=tool.summarize(raw,tmp_path/'results/summaries/audit.json',root=tmp_path)
    assert summary['data']==payload['analysis']['summary']
    assert 'details' not in summary['data']
    assert summary['raw_payload']['sha256']==tool.digest(raw.read_bytes())
    assert summary['regeneration']['argv'][-1]=='.build/results/audit.json'
    assert set(summary['producer_bindings'])==set(tool.OWN_SOURCES)
