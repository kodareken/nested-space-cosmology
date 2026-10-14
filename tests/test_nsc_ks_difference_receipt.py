"""Strict construction and authenticated historical phase replay are distinct."""
from dataclasses import fields, replace
import json
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons import nsc_ks_source_envelope as E
from recursive_horizons.nsc_ks_difference_envelope import KSDifferenceIncoming
from recursive_horizons.nsc_ks_difference_receipt import phase_reconstruction_diagnostic

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = 'results/development/nsc-ks-difference-control.json'
RECEIPT_SHA = '01ac13653c8b787567de920d66c13e11dd6cc343d4c1ff4df6f991df1c463bb0'
CONTEXT = '329ef809e466a86ed4a290ee32c19485f85add9d'
SOURCE_CONTEXTS = {
    'src/recursive_horizons/nsc_local_incoming_family.py':'09f3faa364ca88b8ff482c413022ce38f40bef0d',
    'src/recursive_horizons/nsc_ks_difference_envelope.py':'e62d28f3230ae6ddb2712df71eb7c61a8ed65ad6',
}


@pytest.fixture(scope='module')
def recipe():
    import derive_nsc_ks_difference_control as control
    family=control.R.C.P.family(control.R.C.P.ALPHA)
    name,count,options=control.CASES[-1]
    record=json.loads((ROOT/RECEIPT).read_text())
    with np.load(ROOT/record['payload']['path'],allow_pickle=False) as archive:
        arrays={key:archive[key].copy() for key in archive.files}
    meta=json.loads(arrays['metadata_json'].tobytes())
    mesh=E.computational_z_grid(count)
    w,U=E._sample_axial_profiles(family.directions,mesh)
    binding=E.KSEnvelopeBinding(mesh,w,U,family.amplitudes,
        tuple(d.inner_radius for d in family.directions),tuple(d.outer_radius for d in family.directions),
        E.usual_axial_support(),meta['rho_up'],1.,**options)
    return {'case':name,'binding':binding,'family':family,'source_context':CONTEXT,
            'source_contexts':SOURCE_CONTEXTS}


@pytest.fixture(scope='module')
def authenticated(recipe):
    return KSDifferenceIncoming.from_authenticated_historical_receipt(ROOT,RECEIPT,RECEIPT_SHA,**recipe)


def test_exact_local_arrays_and_separate_phase_indicators(authenticated):
    record=json.loads((ROOT/RECEIPT).read_text())
    with np.load(ROOT/record['payload']['path'],allow_pickle=False) as archive:
        metadata=json.loads(archive['metadata_json'].tobytes())
        original_fingerprint=metadata['cases']['tighter512']['binding']
        for attr,key in [('columns','F'),('column_tangents','dF'),('axial_columns','Fz'),
                         ('axial_tangents','dFz'),('reference_amplitudes','A'),
                         ('envelope_difference','D'),('envelope_difference_z','Dz')]:
            actual=getattr(authenticated,attr)
            expected=archive['tighter512/'+key]
            assert actual.dtype==expected.dtype
            assert actual.tobytes()==expected.tobytes()
            assert not actual.flags.writeable
    diagnostic=authenticated.diagnostics['phase_replay_diagnostic']
    assert diagnostic['physical_field_error_bound'] is None
    assert diagnostic['integration_error_bound'] is None
    assert diagnostic['physical_constraint_status']=='OPEN'
    assert diagnostic['maximum_F_arithmetic_bound']>0
    assert diagnostic['maximum_F_library_accuracy_indicator']>0
    # Original and current identities are separately authenticated. Their
    # hashes can agree when the current profile bytes match the captured ones.
    assert authenticated.diagnostics['original_binding_fingerprint']==original_fingerprint
    binding=authenticated.binding
    current_fingerprint=E._profile_fingerprint(
        binding.computational_z,binding.w_samples,binding.U_samples,binding.amplitudes,
        binding.inner_radii,binding.outer_radii,binding.axial_support,
        binding.rho_up,binding.rho_sigma,binding.rtol,binding.atol,binding.max_step)
    assert current_fingerprint==binding.fingerprint
    assert authenticated.diagnostics['current_binding_fingerprint']==authenticated.binding.fingerprint
    assert authenticated.validate()


def test_default_constructor_stays_exact(authenticated):
    values={field.name:getattr(authenticated,field.name) for field in fields(KSDifferenceIncoming)}
    # Even a one-bit phase change cannot be accepted by normal construction.
    values['columns']=values['columns'].copy()
    phase=np.exp(-1j*values['source_energies']*values['z'][:,None])[:,None]
    total=values['reference_amplitudes'][None]+values['envelope_difference']
    values['columns']=phase*total
    values['columns'][0,0,0]+=1e-12
    with pytest.raises(ValueError,match='columns must reconstruct'):
        KSDifferenceIncoming(**values)


def test_missing_current_payload_never_uses_git(tmp_path,recipe):
    target=tmp_path/RECEIPT
    target.parent.mkdir(parents=True)
    target.write_bytes((ROOT/RECEIPT).read_bytes())
    with pytest.raises(FileNotFoundError,match='current source unavailable'):
        KSDifferenceIncoming.from_authenticated_historical_receipt(tmp_path,RECEIPT,RECEIPT_SHA,**recipe)


def test_changed_local_payload_is_rejected_before_source_or_phase(tmp_path,recipe):
    raw=(ROOT/RECEIPT).read_bytes()
    record=json.loads(raw)
    target=tmp_path/RECEIPT;target.parent.mkdir(parents=True);target.write_bytes(raw)
    payload=tmp_path/record['payload']['path'];payload.parent.mkdir(parents=True,exist_ok=True)
    changed=bytearray((ROOT/record['payload']['path']).read_bytes());changed[-1]^=1
    payload.write_bytes(changed)
    with pytest.raises(RuntimeError,match='current source drift'):
        KSDifferenceIncoming.from_authenticated_historical_receipt(tmp_path,RECEIPT,RECEIPT_SHA,**recipe)


def test_changed_receipt_hash_and_missing_explicit_context_are_rejected(recipe):
    with pytest.raises(RuntimeError,match='current source drift'):
        KSDifferenceIncoming.from_authenticated_historical_receipt(ROOT,RECEIPT,'0'*64,**recipe)
    wrong={**recipe,'source_context':None}
    with pytest.raises(ValueError,match='explicit full'):
        KSDifferenceIncoming.from_authenticated_historical_receipt(ROOT,RECEIPT,RECEIPT_SHA,**wrong)


def test_wrong_source_pin_is_rejected(recipe):
    wrong={**recipe,'source_contexts':{**SOURCE_CONTEXTS,
        'src/recursive_horizons/nsc_ks_difference_envelope.py':CONTEXT}}
    with pytest.raises(RuntimeError,match='frozen source drift'):
        KSDifferenceIncoming.from_authenticated_historical_receipt(ROOT,RECEIPT,RECEIPT_SHA,**wrong)


def test_changed_current_binding_and_case_are_rejected(recipe):
    wrong={**recipe,'binding':replace(recipe['binding'],rtol=1e-6,fingerprint='')}
    with pytest.raises(ValueError,match='solver/domain'):
        KSDifferenceIncoming.from_authenticated_historical_receipt(ROOT,RECEIPT,RECEIPT_SHA,**wrong)
    wrong={**recipe,'case':'foreign-case'}
    with pytest.raises(ValueError,match='completed difference case'):
        KSDifferenceIncoming.from_authenticated_historical_receipt(ROOT,RECEIPT,RECEIPT_SHA,**wrong)


@pytest.mark.parametrize('field',['columns','axial_columns'])
def test_phase_mutation_exceeds_derived_envelope(authenticated,field):
    state=authenticated
    F,Fz=state.columns.copy(),state.axial_columns.copy()
    (F if field=='columns' else Fz)[0,0,0]+=1e-9
    with pytest.raises(ValueError,match='exceeds finite replay'):
        phase_reconstruction_diagnostic(state.z,state.source_energies,state.reference_amplitudes,
            state.envelope_difference,state.envelope_difference_z,F,Fz)


def test_zero_column_mass_has_no_absolute_allowance():
    z=np.array([1.]);energies=np.array([2.]);A=np.zeros((2,1),complex)
    D=np.zeros((1,2,1),complex)
    F=D.copy();F[0,0,0]=np.nextafter(0.,1.)
    with pytest.raises(ValueError,match='exceeds finite replay'):
        phase_reconstruction_diagnostic(z,energies,A,D,D,F,D)


def test_finite_operands_with_overflow_do_not_bypass_the_diagnostic():
    with np.errstate(over='ignore',invalid='ignore'):
        with pytest.raises(ValueError,match='nonfinite phase replay'):
            phase_reconstruction_diagnostic(np.array([1e308]),np.array([1e308]),
                np.ones((2,1),complex),*[np.zeros((1,2,1),complex)]*4)


@pytest.mark.parametrize('dtype',[str,bool,object])
def test_phase_operands_do_not_coerce_invalid_types(dtype):
    with pytest.raises(ValueError,match='numeric phase operands'):
        phase_reconstruction_diagnostic(np.array([1.]),np.array([2.]),
            np.zeros((2,1)).astype(dtype),*[np.zeros((1,2,1),complex)]*4)
