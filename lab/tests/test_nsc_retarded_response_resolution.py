"""Focused layout/gating checks and replay; no new scientific producers."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('retarded_resolution',ROOT/'scripts/derive_nsc_retarded_response_resolution.py')
R=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(R)


def test_nested_grid_preserves_spin_energy_source_layout():
    x=np.arange(9.);phi=np.arange(18*12).reshape(18,12)
    xx,ff=R.subgrid(x,phi,5)
    assert np.array_equal(xx,x[::2])
    assert np.array_equal(ff[:5],phi[:9:2])
    assert np.array_equal(ff[5:],phi[9::2])
    with pytest.raises(ValueError,match='nested'):R.subgrid(x,phi,6)


def test_comparison_contract_restricts_both_covariance_time_axes():
    def fixture(points,times):
        field=np.ones((2*points,12),complex)
        return {'x':np.linspace(0,1,points),'times':np.linspace(0,1,times),
                'field_tangent':field,'trace_tangent':np.ones((times,2,12),complex),
                'covariance_tangent':np.ones((2*times,2*times),complex)}
    c=fixture(5,3);f=fixture(9,3)
    assert all(v['passes'] for v in R.compare(c,f,'spatial',.25).values())
    c=fixture(9,3);f=fixture(9,5)
    assert all(v['passes'] for v in R.compare(c,f,'temporal',.1).values())
    f['covariance_tangent'][0,0]=2
    assert not R.compare(c,f,'temporal',.1)['covariance']['passes']


def test_outside_support_retains_late_values_as_diagnostic():
    t=np.linspace(0,.3,7);a=np.ones((7,2,12),complex)
    value=R.outside_causal_support(a,t)
    assert value['maximum_entry']==1
    assert 0<value['sampled_energy_fraction']<1
    assert not value['directed_bound'] and not value['acceptance_gate']


def test_saved_replay_binds_sources_and_full_coherent_kernels(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('producer called during saved replay')
    monkeypatch.setattr(R,'reference_fields_on_grid',forbidden)
    monkeypatch.setattr(R,'one_case',forbidden)
    monkeypatch.setattr(R.P,'evolve_prepared_field_jets',forbidden)
    record=R.check()
    assert record['control']['reference_sampling_calls']==1
    assert record['residuals']['archived801_resampling']<3e-11
    assert len(record['control']['cases_completed'])<=3
    assert record['scope']['continuum_error_bound'] is None
    assert not record['scope']['full_CAR_verified']
    with np.load(ROOT/record['payload']['path'],allow_pickle=False) as archive:
        a={k:archive[k] for k in archive.files}
    rootw=np.repeat(np.sqrt(a['source/weights']/(2*np.pi)),3)
    for case in record['control']['cases_completed']:
        c=R.case_arrays(a,case);F=(c['trace']*rootw).reshape(-1,12);dF=(c['trace_tangent']*rootw).reshape(-1,12)
        for key,source in (('covariance_tangent','source'),('projector_tangent','projector')):
            C=R.P.block_diag(*a['source/'+source]);expected=dF@C@F.conj().T+F@C@dF.conj().T
            assert np.max(abs(expected-c[key]))<1e-15
        assert np.max(abs(c['trace_tangent'][0]))==0
        assert np.array_equal(c['x'][c['sample_rows'][0]],1.)
    completed=record['control']['cases_completed']
    if '3201_256' in completed:
        assert all(v['passes'] for v in record['comparison']['spatial'].values())
    if not record['runtime']['cpu_budget_exceeded']:
        assert record['runtime']['CPU_seconds']<=120
