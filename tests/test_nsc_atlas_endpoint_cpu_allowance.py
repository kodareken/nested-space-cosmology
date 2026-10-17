"""CPU allowance bookkeeping and rollback; no physical or remote execution."""
from copy import deepcopy
import json
import os
from types import SimpleNamespace

import pytest

import amend_nsc_atlas_endpoint_cpu as controller
import run_nsc_atlas_fixed_endpoint as runner


def setup_job(tmp_path, monkeypatch):
    monkeypatch.setattr(controller, 'ROOT', tmp_path)
    monkeypatch.setattr(controller, 'RUNNER', tmp_path/'scripts/run_nsc_atlas_fixed_endpoint.py')
    monkeypatch.setattr(controller.time, 'process_time', lambda: 8.)
    (tmp_path/'results/summaries').mkdir(parents=True)
    budget = SimpleNamespace(cpu_seconds=10000., started=5., evaluations=123,
                             maximum_evaluations=500_000_000_000)
    report = dict(schema='nsc-atlas-fixed-endpoint-v1', status='incomplete',level=0,
        source_commit=controller.V4_SOURCE,
        limits=dict(cpu_seconds=10000., charged_values=500_000_000_000,
                    accepted_steps=2000, budget_reset=False),
        radial_accuracy=dict(rtol=2.5e-14,atol=1e-15,gradient_absolute_tolerance=2e-8),
        retained_state=dict(field=[.2,.3], radius=.51))
    args = SimpleNamespace(level=0,output=tmp_path/'results/summaries/nsc-atlas-fixed-endpoint-level0-v4.json')
    namespace=dict(controller=controller, os=os)
    source='''def main(args,budget,report,accepted,last_x):
    controller.apply_to_current_process(pid=os.getpid(),level=0,
        receipt='results/summaries/amendment.json',provenance={'review':'fixture'})
'''
    exec(compile(source,str(controller.RUNNER),'exec'),namespace)
    return budget,report,args,namespace['main']


def test_extension_changes_only_cpu_ceiling_and_preserves_counters(tmp_path,monkeypatch):
    budget,report,args,main=setup_job(tmp_path,monkeypatch)
    before=deepcopy(report); counters=vars(budget).copy()
    main(args,budget,report,720,.457)
    saved=json.loads((tmp_path/'results/summaries/amendment.json').read_text())
    assert saved==report.pop('cpu_allowance_amendment')
    assert saved['elapsed_cpu_seconds']==3. and saved['charged_values']==123
    assert saved['accepted_steps']==720 and saved['last_accepted_x']==.457
    assert saved['physical_state_changed'] is saved['numerical_settings_changed'] is False
    assert budget.cpu_seconds==report['limits'].pop('cpu_seconds')==20000.
    before['limits'].pop('cpu_seconds')
    assert report==before
    counters['cpu_seconds']=20000.
    assert vars(budget)==counters


def test_second_extension_is_rejected_without_erasing_first(tmp_path,monkeypatch):
    budget,report,args,main=setup_job(tmp_path,monkeypatch)
    main(args,budget,report,720,.457)
    before=deepcopy(report)
    with pytest.raises(ValueError,match='unamended'):
        controller.extend_allowance(budget,report,accepted=720,last_x=.457,provenance={})
    assert report==before and budget.cpu_seconds==20000.


@pytest.mark.parametrize('failure',['expired','source','charged'])
def test_invalid_target_has_no_mutation_or_receipt(tmp_path,monkeypatch,failure):
    budget,report,args,main=setup_job(tmp_path,monkeypatch)
    if failure=='expired':budget.started=-10000.
    elif failure=='source':report['source_commit']='b'*40
    else:budget.evaluations=budget.maximum_evaluations+1
    before=deepcopy(report); counters=vars(budget).copy()
    with pytest.raises(ValueError):main(args,budget,report,720,.457)
    assert report==before and vars(budget)==counters
    assert not (tmp_path/'results/summaries/amendment.json').exists()


def test_failed_receipt_write_rolls_back_all_bookkeeping(tmp_path,monkeypatch):
    budget,report,args,main=setup_job(tmp_path,monkeypatch)
    before=deepcopy(report); counters=vars(budget).copy()
    def fail(_):raise OSError('manufactured receipt failure')
    monkeypatch.setattr(controller.os,'fsync',fail)
    with pytest.raises(OSError,match='manufactured'):main(args,budget,report,720,.457)
    assert report==before and vars(budget)==counters
    assert not (tmp_path/'results/summaries/amendment.json').exists()


def test_fresh_replay_selects_same_finite_cpu_ceiling(tmp_path,monkeypatch):
    class StopBeforeAnyCalculation(Exception):pass
    seen={}
    def capture(**kwargs):
        seen.update(kwargs)
        raise StopBeforeAnyCalculation
    monkeypatch.setattr(runner,'Budget',capture)
    monkeypatch.setattr(runner,'frozen_sources',lambda:('a'*40,{}))
    monkeypatch.setattr(runner.sys,'argv',['runner','--level','2','--extended-cpu-allowance',
                                       '--output',str(tmp_path/'fresh.json')])
    with pytest.raises(StopBeforeAnyCalculation):runner.main()
    assert seen==dict(maximum_evaluations=500_000_000_000,cpu_seconds=80000.)
    assert runner.level_spec(2)['cpu_seconds']==40000.
