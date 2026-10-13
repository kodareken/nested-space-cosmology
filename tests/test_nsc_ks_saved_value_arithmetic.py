"""Finite-operation comparison does not alter source or physical contracts."""
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_evolved_incoming_constraints import coefficients_from_records
from recursive_horizons.nsc_ks_saved_value_arithmetic import (
    _Ball, _profile, coefficients_midpoint_roundoff, saved_value_replay_envelope,
    require_saved_value_summaries, require_saved_value_receipt,
)
from recursive_horizons.nsc_ks_value_evaluator import geometry_change, edge_change
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def finite_recipe():
    rawc=json.loads((ROOT/'results/development/nsc-incoming-surface-coefficients.json').read_text())
    rawb=json.loads((ROOT/'results/development/nsc-incoming-surface-regular-branch.json').read_text())
    c=coefficients_from_records(rawc,rawb)
    errors=coefficients_midpoint_roundoff(rawc,rawb)
    params=np.zeros((2,8));params[0,1]=.001;params[1,2]=.003
    family=LocalIncomingFamily(params)
    nodes=family.collocation_nodes(3)
    weights=np.array([107952.,107952.])
    geometry=geometry_change(family,nodes,c)
    edge=edge_change(family,nodes,weights,c['a'],1.03,8)
    arguments=(family,nodes,c,np.array([.2,.3]),geometry,edge,
               [np.zeros((3,2)) for _ in range(60)],weights,1.03,8)
    policy=saved_value_replay_envelope(*arguments,coefficient_errors=errors)
    return arguments,errors,policy


def test_profile_error_covers_owned_clenshaw_derivatives(finite_recipe):
    args,_,_=finite_recipe
    family,nodes=args[:2]
    for row,function in zip(family.coefficients,family.functions,strict=True):
        for order in range(4):
            bounded=_profile(row,family,_Ball(nodes),order)
            actual=np.array([function(node,order) for node in nodes])
            assert np.all(abs(actual-bounded.v)<=bounded.e)


def test_arithmetic_and_unproved_configured_indicators_are_separate(finite_recipe):
    _,_,policy=finite_recipe
    assert np.all(policy['arithmetic']>0)
    assert np.all(policy['coordinate_accuracy_indicator']>=0)
    assert np.all(policy['summary_envelope']>=policy['arithmetic'])
    assert policy['physical_error_bound'] is None
    assert policy['phase_quadrature_error_bound'] is None
    assert policy['physical_constraint_status']=='OPEN'
    assert 'not proved' in policy['scope']


def test_authentic_array_count_and_shape_are_required(finite_recipe):
    args,errors,_=finite_recipe
    changed=list(args);changed[6]=changed[6][:-1]
    with pytest.raises(ValueError,match='60 authenticated'):
        saved_value_replay_envelope(*changed,coefficient_errors=errors)
    changed=list(args);changed[6]=[np.zeros((2,2))]*60
    with pytest.raises(ValueError,match='60 authenticated'):
        saved_value_replay_envelope(*changed,coefficient_errors=errors)


@pytest.mark.parametrize('index', [4,5])
def test_current_part_mutation_cannot_enlarge_its_allowance(finite_recipe,index):
    args,errors,_=finite_recipe
    changed=list(args)
    changed[index]=args[index].copy()
    changed[index][0,0]+=1.
    with pytest.raises(ValueError,match='changed beyond finite operation model'):
        saved_value_replay_envelope(*changed,coefficient_errors=errors)


@pytest.mark.parametrize('dtype', [complex,str,bool])
def test_nonnumeric_or_complex_inputs_cannot_be_coerced(finite_recipe,dtype):
    args,errors,_=finite_recipe
    changed=list(args)
    changed[6]=[a.astype(dtype) for a in args[6]]
    with pytest.raises(ValueError,match='real numeric'):
        saved_value_replay_envelope(*changed,coefficient_errors=errors)


def test_endpoint_cancellation_is_charged_before_midpoint():
    # Reuse the real interval schema, changing only an in-memory endpoint pair.
    rawc=json.loads((ROOT/'results/development/nsc-incoming-surface-coefficients.json').read_text())
    rawb=json.loads((ROOT/'results/development/nsc-incoming-surface-regular-branch.json').read_text())
    rawb['branch']['imported_unchanged']['A0']=[-1e6,1e6]
    errors=coefficients_midpoint_roundoff(rawc,rawb)
    assert errors['A0']>0  # zero midpoint is not a zero rounding allowance


def test_summary_mutation_outside_bound_is_rejected(finite_recipe):
    _,_,policy=finite_recipe
    residual=np.array([[.001,.0003],[.0001,.0004],[.0002,.0001]])
    saved={'merit':.001,'constraint_maxima':[.001,.0004]}
    require_saved_value_summaries(saved,residual,policy['summary_envelope'])
    wrong=deepcopy(saved)
    wrong['constraint_maxima'][0]+=max(1.,2*np.max(policy['summary_envelope']))
    wrong['merit']=wrong['constraint_maxima'][0]
    with pytest.raises(ValueError,match='beyond replay envelope'):
        require_saved_value_summaries(wrong,residual,policy['summary_envelope'])


@pytest.mark.parametrize('mutation', ['source_hashes','physical_error_bound','status','extra','typed_merit'])
def test_nonallowed_receipt_drift_is_rejected(finite_recipe,mutation):
    _,_,policy=finite_recipe
    saved={'merit':.001,'constraint_maxima':[.001,.0004],
           'source_hashes':{'src/owner.py':'a'*64},'physical_error_bound':None,'status':'OPEN'}
    actual=deepcopy(saved)
    if mutation=='source_hashes':actual[mutation]['src/owner.py']='b'*64
    elif mutation=='physical_error_bound':actual[mutation]=0.
    elif mutation=='status':actual[mutation]='CERTIFIED'
    elif mutation=='extra':actual['extra']=True
    else:actual['merit']='0.001'
    with pytest.raises(ValueError):
        require_saved_value_receipt(actual,saved,policy['summary_envelope'])
