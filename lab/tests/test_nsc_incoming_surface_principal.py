"""Principal identities, actual density controls and saved-cu directed replay."""
import importlib.util
import json
from pathlib import Path

import numpy as np

import recursive_horizons.nsc_incoming_surface_principal as owner
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets,IncomingNormalJetChange
from recursive_horizons.nsc_spatial_reference_symbol import INDEX

ROOT=Path(__file__).resolve().parents[1]


def data():
    channels=json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    ledger=json.loads((ROOT/'results/development/nsc-incoming-local-constraints.json').read_text())['locked_inputs']
    cu=json.loads((ROOT/'results/development/nsc-incoming-lapse-coefficient.json').read_text())['coefficient']['total_coefficient_interval']
    return channels,ledger,cu


def test_full_weyl_shift_response_and_new_principal_moments():
    identity=owner.linear_weyl_response_identity()
    assert all(v=='0' for k,v in identity.items() if k.endswith('_residual'))
    principal=owner.reference_shift_principal()
    assert principal['kernel_identity_residuals']==['0','0']
    assert principal['moment_identity_residuals']==['0','0']
    assert principal['beta_wave_ratio_residual']=='0'
    combined=owner.combined_principal_identities()
    assert combined['reference_reused_cu_relation']=='0'
    assert set(combined['local_reused_cu_relations'].values())=={'0'}
    assert combined['determinant_identity_residual']==combined['elimination_identity_residual']=='0'


def test_full_local_euler_extraction_matches_original_density_controls():
    exact=owner.local_euler_principal_identity()
    assert {v for row in exact.values() for v in row.values()}=={'0'}
    channels,ledger,_=data();directed=owner.directed_d_components(channels,ledger)
    for step in (.5,1.):
        for name,row in owner.local_principal_controls(ledger,step).items():
            for key,value in row.items():
                lo,hi=directed['local_entries'].get(name,{}).get(key+'_interval',[0.,0.])
                assert lo-3e-13<=value<=hi+3e-13


def test_compatible_six_slot_factorials_are_the_owned_derivative_values():
    from math import factorial
    slots=((1,0),(1,1),(1,2),(1,3),(3,0),(3,1))
    values=(.2,.3,.4,.5,.6,.7)
    domain=incoming_cauchy_jets(tuple(IncomingNormalJetChange('r',t,z,v) for (t,z),v in zip(slots,values)))
    for (t,z),value in zip(slots,values):
        change=domain.fields[3].data[INDEX[t,z,0],0,0,0]-domain.baseline[3].data[INDEX[t,z,0],0,0,0]
        assert abs(change-value/(factorial(t)*factorial(z)))<3e-15


def test_directed_matrix_excludes_zero_and_reuses_saved_cu():
    channels,ledger,cu=data();result=owner.directed_principal_matrix(channels,ledger,cu)
    A,B=result['matrix_intervals'][0];d,e=result['matrix_intervals'][1]
    assert A[0]<=cu[0]<cu[1]<=A[1]
    assert A[0]>0 and B[1]<0 and d[1]<0 and e[0]>0
    assert result['determinant_interval'][1]<0
    assert result['direct_matrix_determinant_interval'][1]<0
    assert result['eliminated_w3_coefficient_interval'][1]<0
    assert result['strictly_nonzero']


def test_replay_runs_no_old_or_new_symbolic_producer(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('coefficient producer entered during directed replay')
    for name in ('local_principal_expressions','local_euler_principal_identity','local_principal_controls',
                 'reference_shift_principal','linear_weyl_response_identity','combined_principal_identities'):
        monkeypatch.setattr(owner,name,forbidden)
    import recursive_horizons.nsc_incoming_lapse_coefficient as cu
    import recursive_horizons.nsc_incoming_shift_coefficient as cv
    monkeypatch.setattr(cu,'closed_lapse_coefficient',forbidden);monkeypatch.setattr(cv,'closed_shift_coefficient',forbidden)
    spec=importlib.util.spec_from_file_location('principal_replay',ROOT/'scripts/derive_nsc_incoming_surface_principal.py')
    script=importlib.util.module_from_spec(spec);spec.loader.exec_module(script)
    prior=json.loads((ROOT/script.OUTPUT).read_text())
    assert script.make_record(ROOT/prior['payload']['path'])==prior
    assert not prior['scope']['cu_cv_producers_rerun']
    assert not prior['scope']['spatial_domain_boundary_length_selected']
