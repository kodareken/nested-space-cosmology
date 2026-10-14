"""No-evolution diagnosis replay and the owned interior symbol discriminator."""
import importlib.util
import json
from pathlib import Path
import numpy as np
import pytest
import sympy as sp
from recursive_horizons.nsc_frozen_sources import frozen_source_context
from recursive_horizons.nsc_retarded_saved_arithmetic import (
    localization_allowances, require_diagnostic_receipt)
ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('spatial_diagnosis',ROOT/'scripts/derive_nsc_retarded_spatial_diagnosis.py')
D=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(D)
ORIGINAL_SOURCE='329ef809e466a86ed4a290ee32c19485f85add9d'


@pytest.fixture(scope='module')
def saved_localization_arrays():
    # The original source scope authenticates provenance only. The pilot check
    # and NPZ load still require the exact current local numerical bytes.
    with frozen_source_context(ORIGINAL_SOURCE):
        receipt=D.pilot.check()
    with np.load(ROOT/receipt['payload']['path'],allow_pickle=False) as payload:
        return {key:payload[key] for key in payload.files}


def drift_maximum_record(case,name,value):
    return {'localization':{'per_case':{case:{name:{'maximum':value}}}}}


@pytest.mark.parametrize('case',['401_64','801_64','801_128'])
@pytest.mark.parametrize('name',['stationary_field_drift','stationary_trace_drift'])
def test_drift_maxima_charge_both_saved_operand_evaluations(saved_localization_arrays,case,name):
    arrays=saved_localization_arrays
    key='field' if name=='stationary_field_drift' else 'trace'
    maximum=float(np.max(abs(arrays[case+'/'+key]-arrays[case+'/stationary_'+key])))
    path=('localization','per_case',case,name,'maximum')
    allowances=localization_allowances(arrays)
    # No locations, norms, source bindings, or physical criteria are added to
    # the allowlist by the six drift-maximum comparison rules.
    drift_paths={p for p in allowances if p[:2]==('localization','per_case')}
    assert drift_paths=={
        ('localization','per_case',c,n,'maximum')
        for c in ('401_64','801_64','801_128')
        for n in ('stationary_field_drift','stationary_trace_drift')}
    actual=drift_maximum_record(case,name,maximum)
    rounded=drift_maximum_record(case,name,float(np.nextafter(maximum,np.inf)))
    differences=require_diagnostic_receipt(actual,rounded,allowances)
    assert len(differences)==1
    assert differences[0]['path']==list(path)
    cap=allowances[path]
    assert np.isfinite(cap) and cap>0
    beyond=float(np.nextafter(maximum+cap,np.inf))
    assert beyond-maximum>cap
    with pytest.raises(ValueError,match='beyond finite operation policy'):
        require_diagnostic_receipt(actual,drift_maximum_record(case,name,beyond),allowances)


@pytest.mark.parametrize('key',['field','trace'])
def test_cancelled_drift_maximum_retains_uncancelled_operand_charge(saved_localization_arrays,key):
    arrays=dict(saved_localization_arrays)
    case='801_128'
    a=np.full(arrays[case+'/'+key].shape,1+2j,dtype=np.complex128)
    b=a-complex(2.**-40,0.)
    arrays[case+'/'+key],arrays[case+'/stationary_'+key]=a,b
    path=('localization','per_case',case,'stationary_'+key+'_drift','maximum')
    cap=localization_allowances(arrays)[path]
    maximum=float(np.max(abs(a-b)))
    # Independent scale check: subtracting nearly equal large operands must
    # not charge only the tiny final drift. This is a host-arithmetic policy,
    # and contributes no physical or continuum error bound.
    assert cap>=np.finfo(float).eps*float(np.max(abs(a)+abs(b)))
    assert cap>1000*np.finfo(float).eps*maximum


@pytest.mark.parametrize('changed',[1,'0.0',False,float('nan'),float('inf')])
def test_drift_allowance_cannot_exempt_numeric_types_or_nonfinite_values(changed):
    path=('localization','per_case','401_64','stationary_field_drift','maximum')
    original=drift_maximum_record('401_64','stationary_field_drift',0.)
    with pytest.raises(ValueError,match='diagnostic receipt'):
        require_diagnostic_receipt(original,drift_maximum_record('401_64','stationary_field_drift',changed),{path:1.})


def test_drift_allowance_does_not_exempt_peak_identity_or_bound_status():
    path=('localization','per_case','401_64','stationary_field_drift','maximum')
    for original,changed in [
        ({'location':{'rho':1.}}, {'location':{'rho':float(np.nextafter(1.,np.inf))}}),
        ({'matching_status':'OPEN'}, {'matching_status':'PASS'}),
        ({'physical_error_bound':None}, {'physical_error_bound':0.}),
        ({'source_hash':'original'}, {'source_hash':'changed'})]:
        with pytest.raises(ValueError,match='diagnostic receipt'):
            require_diagnostic_receipt(original,changed,{path:1.})


def test_owned_fourth_order_group_factor_is_symbol_derivative():
    theta=sp.Symbol('theta',real=True)
    symbol=(8*sp.sin(theta)-sp.sin(2*theta))/6
    assert sp.trigsimp(sp.diff(symbol,theta)-(4*sp.cos(theta)-sp.cos(2*theta))/3)==0


def test_saved_and_local_diagnosis_replays_without_any_evolution(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('propagation/producer entered a diagnosis')
    monkeypatch.setattr(D.pilot,'run',forbidden)
    monkeypatch.setattr(D.pilot,'one_case',forbidden)
    monkeypatch.setattr(D.pilot,'evolve_prepared_field_jets',forbidden)
    monkeypatch.setattr(D.pilot.FourthOrderModePropagator,'evolve',forbidden)
    record,diagnostic=D.replay_recorded_diagnosis(source_context='329ef809e466a86ed4a290ee32c19485f85add9d')
    # Immutable captured diagnostics and a current binary64 recomputation have
    # distinct roles; finite-operation comparisons retain all metadata exactly.
    assert diagnostic['physical_error_bound'] is diagnostic['continuum_error_bound'] is None
    assert diagnostic['matching_status']=='OPEN'
    assert not diagnostic['configured_chart_indicators']['proved_error_bound']
    assert record['scope']['field_propagations']==record['scope']['radial_solves']==0
    assert record['decision']['existing_three_solve_response_status'].startswith('OPEN')
    assert record['residuals']['common_node_forcing_identity']==0
    local=record['local_discriminator']['times']
    assert np.allclose([r['LF_relative_difference'] for r in local],[.5239659969886,.45211274657919764,.4935501286168124],rtol=0,atol=2e-12)
    assert np.allclose([r['grids']['401']['geometry_power_fraction_group_error_above_0_25'] for r in local],
                       [.3371253835272251,.3650634154225146,.3386758253347509],rtol=0,atol=2e-12)
    localization=record['localization']
    assert localization['spatial_field']['unweighted_error_squared_norm_fraction_by_spin'][0]>.999
    assert localization['spatial_field']['spin0_error_fraction_in_rho_0_95_to_1_05']>.90
    assert localization['spatial_trace']['location']['PG_time']==.253125
    assert localization['spatial_field']['unweighted_error_fraction_by_source'][2]==0
    assert not record['decision']['pulse_widening_or_source_alteration_selected']
    assert record['scope']['continuum_error_bound'] is None


def test_conditional_fast_slow_support_identities_and_saved_leakage():
    beta=sp.Symbol('beta',positive=True)
    assert sp.factor(1/(beta+1)-beta/(beta**2-1)+1/(beta**2-1))==0
    assert sp.factor(1/(beta-1)-beta/(beta**2-1)-1/(beta**2-1))==0
    saved=json.loads((ROOT/D.OUTPUT).read_text());support=saved['conditional_causal_support']
    assert np.allclose(support['PG_time_window'],[.05404795482456391,.24595204517543606],rtol=0,atol=2e-14)
    assert support['per_case']['401_64']['location']['PG_time']>support['PG_time_window'][1]
    assert not support['directed_edge_certificate']
    assert not support['discrete_finite_propagation_claimed']
