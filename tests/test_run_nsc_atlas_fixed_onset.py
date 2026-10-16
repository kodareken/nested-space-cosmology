"""Preparation failure accounting and preservation; no physical solve."""
from types import SimpleNamespace

import numpy as np
import pytest

import run_nsc_atlas_fixed_onset as runner
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


@pytest.mark.parametrize('failure',[None,'source','transfer'])
def test_preparation_charges_attempts_once_without_resetting_the_clock(monkeypatch,failure):
    budget=Budget(maximum_evaluations=1000,cpu_seconds=60.)
    budget.check(3);started=budget.started;report={};seen=[]
    background=SimpleNamespace(inputs={'labels':np.array([-2.,0.])})
    def source(*,nodes,geometry_knots,budget):
        assert nodes==geometry_knots==257
        assert budget.started==started
        assert budget.kernel_values_max==997
        seen.append(budget);budget.check(11)
        if failure=='source': raise ValueError('source failure')
        return SimpleNamespace(budget=budget)
    def cap(field,*,derivative_order):
        assert derivative_order==4
        return SimpleNamespace(budget=field.budget,domain={'Ta':-1.})
    def transfer(cap):
        cap.budget.check(17)
        if failure=='transfer': raise ValueError('transfer failure')
        return dict(background=background,onset={'index':1},diagnostics={})
    monkeypatch.setattr(runner,'solve_prepared_field',source)
    monkeypatch.setattr(runner,'CapLeadingField',cap)
    monkeypatch.setattr(runner,'build_prepared_background',transfer)
    if failure:
        with pytest.raises(ValueError,match=failure+' failure'):
            runner.prepare_history(budget,report)
    else:
        assert runner.prepare_history(budget,report)==(background,-1.)
    charged=11 if failure=='source' else 28
    assert report['charged_values']==charged
    assert budget.evaluations==3+charged and budget.started==started
    assert len(seen)==1 and report['source_integrations_attempted']==1
    assert report['source_integrations_completed']==int(failure!='source')


def test_existing_output_is_preserved_before_any_calculation(tmp_path,monkeypatch):
    path=tmp_path/'saved.json';path.write_bytes(b'preserve this receipt')
    monkeypatch.setattr(runner.sys,'argv',['runner','--retained-history','--output',str(path)])
    def forbidden(): raise AssertionError('source discovery must not run')
    monkeypatch.setattr(runner,'frozen_sources',forbidden)
    with pytest.raises(ValueError,match='preserve the existing result'):
        runner.main()
    assert path.read_bytes()==b'preserve this receipt'
