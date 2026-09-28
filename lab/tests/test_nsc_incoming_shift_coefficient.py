"""Exact Weyl sign, full local Euler extraction and stored-control replay."""
import importlib.util
import json
from pathlib import Path

import numpy as np

import recursive_horizons.nsc_incoming_shift_coefficient as owner
from recursive_horizons.nsc_spatial_reference_symbol import SymbolJet,star_order,INDEX,SIGMA

ROOT=Path(__file__).resolve().parents[1]


def test_trace_formula_has_the_original_star_order_sign():
    identity=owner.weyl_trace_identity()
    assert identity['first_star_order_residual']=='0' and identity['rank_one_block_residual']=='0'
    z=np.array([1.,2.,0.]);k=np.array([-.3,.7,0.])
    Q=SymbolJet.constant((np.eye(2)+SIGMA[2])/2)
    Q.data[INDEX[0,1,0]]=np.einsum('i,ijk->jk',z,SIGMA)/2
    Q.data[INDEX[0,0,1]]=np.einsum('i,ijk->jk',k,SIGMA)/2
    actual=star_order(Q,Q,1).value[0]
    expected=-np.einsum('i,ijk->jk',np.cross(z,k),SIGMA)/4
    assert np.max(abs(actual-expected))<3e-15


def test_exact_reference_normalization_and_full_line_coefficients():
    reference=owner.reference_shift_expressions();integrated=owner.integrated_reference_expressions()
    assert reference['normalization_residuals']==['0']*4
    assert integrated['compact_reduction_residuals']==['0','0']
    assert integrated['odd_order_full_line_identity']=='0'
    assert integrated['moments']=={'p2_gap5':'2/(3*M**2)','p6_gap11':'4/(63*M**4)',
                                  'p4_gap11':'16/(315*M**6)','p2_gap11':'32/(315*M**8)'}
    # Order3 is an odd integrand, not an invented pointwise zero.
    assert reference['kernels'][3]!=0 and reference['paired_kernels'][3]==0


def test_full2d_local_euler_and_original_density_bridge():
    identity=owner.local_compact_identities()
    assert set(identity['closed_form_residuals'].values())=={'0'}
    assert set(identity['plane_residuals'].values())=={'0'}
    ledger=json.loads((ROOT/'results/development/nsc-incoming-local-constraints.json').read_text())['locked_inputs']
    assert owner.local_density_controls(ledger)['maximum']<3e-13


def test_positive_coefficient_and_zero_angular_reference():
    channels=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    ledger=json.loads((ROOT/'results/development/nsc-incoming-local-constraints.json').read_text())['locked_inputs']
    result=owner.closed_shift_coefficient(channels,ledger)
    lo,hi=result['total_coefficient_interval']
    assert 16.5759756220584<lo<hi<16.5759756220586
    assert result['strictly_positive']
    for row in result['per_group_reference']:
        if channels[row['group']]['angular_eigenvalue']==0:
            for key in ('reference_order2_interval','reference_order4_interval'):
                assert row[key][0]<=0<=row[key][1]


def test_replay_never_runs_symbolic_or_previous_coefficient_producers(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('producer ran during directed coefficient replay')
    for name in ('reference_shift_expressions','integrated_reference_expressions','local_shift_expressions','local_density_controls'):
        monkeypatch.setattr(owner,name,forbidden)
    import recursive_horizons.nsc_incoming_lapse_coefficient as lapse
    import recursive_horizons.nsc_incoming_reference_response as reference
    monkeypatch.setattr(lapse,'closed_lapse_coefficient',forbidden)
    monkeypatch.setattr(reference,'incoming_reference_response',forbidden)
    spec=importlib.util.spec_from_file_location('shift_coefficient_replay',ROOT/'scripts/derive_nsc_incoming_shift_coefficient.py')
    script=importlib.util.module_from_spec(spec);spec.loader.exec_module(script)
    prior=json.loads((ROOT/script.OUTPUT).read_text())
    assert not prior['failures']
    assert not prior['scope']['c_u_proof_rerun'] and not prior['scope']['c_v2_newly_certified']
    assert script.make_record(ROOT/prior['payload']['path'])==prior
