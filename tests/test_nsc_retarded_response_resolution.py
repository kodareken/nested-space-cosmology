"""Focused layout/gating checks and replay; no new scientific producers."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import pytest

from recursive_horizons.nsc_frozen_sources import frozen_source_context
# Original producer bytes are explicit; numerical fixtures remain current local inputs.
FROZEN_PRODUCER_CONTEXT = "329ef809e466a86ed4a290ee32c19485f85add9d"

ROOT=Path(__file__).resolve().parents[1]
from recursive_horizons.nsc_retarded_radial_response import resolution_control as R


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


def test_zero_energy_fraction_cannot_be_a_boolean():
    values=np.zeros((3,2,12),complex);times=np.linspace(0.,1.,3)
    with pytest.raises(ValueError,match='invalid sampled-energy'):
        R.energy_fraction_roundoff(values,times,False,0.)


@frozen_source_context(FROZEN_PRODUCER_CONTEXT)
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


@frozen_source_context(FROZEN_PRODUCER_CONTEXT)
def test_energy_fraction_replay_rejects_change_beyond_reduction_roundoff():
    from copy import deepcopy
    record = R.check()
    with np.load(ROOT / record['payload']['path'], allow_pickle=False) as archive:
        arrays = {k: archive[k] for k in archive.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    case = meta['completed'][0]
    a = R.case_arrays(arrays, case)
    field = 'characteristic_trace_outside_causal_support'
    value = record['per_case'][case][field]['sampled_energy_fraction']
    bound = R.energy_fraction_roundoff(a['characteristic_trace_tangent'], a['times'], value, value)
    changed = deepcopy(record)
    changed['per_case'][case][field]['sampled_energy_fraction'] += 4 * bound
    with pytest.raises(ValueError, match='roundoff envelope'):
        R.validate_replay(record, changed, arrays, meta)
    changed = deepcopy(record)
    changed['per_case'][case][field]['acceptance_gate'] = True
    with pytest.raises(ValueError, match='receipt differs'):
        R.validate_replay(record, changed, arrays, meta)
    assert R.energy_fraction_roundoff(a['characteristic_trace_tangent'], a['times'], 0., 0.) == 0.
