#!/usr/bin/env python3
"""Dirac-free provenance checks for the truncated trail and the ruler."""
import json
import sys
from pathlib import Path

import numpy as np
import pytest
from hashlib import sha256

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_ks_coupled_newton_n32_truncated_iterate as N32
import derive_nsc_ks_coupled_newton_n64_high_modes as N64
import derive_nsc_ks_gate_ruler_drift as Ruler
from recursive_horizons.nsc_frozen_sources import frozen_source_context
from recursive_horizons.nsc_ks_campaign_inputs import (
    ExternalRawInputsUnavailable,descriptor_inventory,raw_input_availability,
    require_campaign_raw_inputs,validate_truncated_linear_step,
    stage0_capture_source_provenance,require_original_campaign_sources)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

SOURCE_CONTEXT='329ef809e466a86ed4a290ee32c19485f85add9d'
CAPTURED_PROPOSAL_SHA256='ffeb0f442eb5e29766687c3b28af73b39341a80e02640f10874f705b1cc0e898'
DESCRIPTOR_CENSUS={
'nsc-ks-coupled-newton-n32-truncated-iterate2':'cd7e75f8056f7d970c2a893f40740283d79b15eb404595f37cbf6a95afe89cc7',
'nsc-ks-coupled-newton-n32-truncated-iterate3':'fe54ac81c805419f933590f5952d62d18c04617bfb4a177d83d55aaf407f4636',
'nsc-ks-coupled-newton-n32-truncated-iterate4':'927c6f294ac1630abbd9eadf114ef9c2cbe7050cf7c72361a2c9b468aad6e128',
'nsc-ks-coupled-newton-n32-truncated-iterate5':'49b7fc55d313668b396fc4ac7a476f679ca6c0041a0ca63b52e9a8439b64938a',
'nsc-ks-coupled-newton-n32-truncated-iterate6':'55b7d689a243e807abccf0e4602bde89617426678916e919ad552000868039dd',
'nsc-ks-coupled-newton-n32-truncated-trial':'4fcebf334b687b0e3a87d67dd49f4b8daedf068daa7bca57c3d0936684e9e022',
'nsc-ks-coupled-newton-n64-high-modes':'c6577c981dd42540d0ea4679e9b0496393c76d07ee78fe2a4cd26c9b0a9a9723',
}

ITERATE6 = 'b4b53f0bba7290cfde14af741c82fda9af1b49a13f7a3b207a7230f39ea106e9'
TRAIL = (
    'nsc-ks-coupled-newton-n32-truncated-trial',
    'nsc-ks-coupled-newton-n32-truncated-iterate2',
    'nsc-ks-coupled-newton-n32-truncated-iterate3',
    'nsc-ks-coupled-newton-n32-truncated-iterate4',
    'nsc-ks-coupled-newton-n32-truncated-iterate5',
    'nsc-ks-coupled-newton-n32-truncated-iterate6',
)


def test_name_rejects_path_characters():
    for name in ('../etc', 'iterate6.json', 'Iterate', 'a/b', ''):
        try:
            N32.paths(name)
        except ValueError:
            continue
        raise AssertionError('path-like name was accepted: '+name)
    directory, output, proposal = N32.paths('iterate6')
    assert directory.name == 'nsc-ks-coupled-newton-n32-truncated-iterate6'
    assert output.parent == proposal.parent


def captured_proposal():
    path=ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6-proposal.json'
    raw=path.read_bytes()
    assert sha256(raw).hexdigest()==CAPTURED_PROPOSAL_SHA256
    return json.loads(raw)


def test_campaign_best_regenerates_iterate6_identity():
    # The immutable captured witness owns its exact binary64 identity.
    # A fresh SVD has its own finite linear-system contract on this runtime.
    parent_path=ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate5.json'
    best_path=ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate3.json'
    parent=json.loads(parent_path.read_text());best=json.loads(best_path.read_text())
    saved=captured_proposal()
    captured=LocalIncomingFamily(np.asarray(saved['candidate_coefficients'],float))
    assert saved['candidate_profile_identity']==ITERATE6
    assert N32.profile_identity(captured,include_normal_window=True)==ITERATE6
    assert saved['condition_cutoff']==1e6 and saved['step_scale']==1.
    assert saved['step_norm']<=1.+1e-12 and captured.radius_lower_bound()>0
    provenance=stage0_capture_source_provenance(ROOT,saved['source_hashes'])
    assert provenance['creating_CLI_authenticated']
    assert provenance['inherited_math_authenticated'] is False
    with frozen_source_context(SOURCE_CONTEXT):
        for path,digest in saved['input_hashes'].items():
            assert N32.digest(path)==digest
    proposal=N32.build_proposal(parent_path,1e6,1.,best_path)
    for path,digest in {**proposal['source_hashes'],**proposal['input_hashes']}.items():
        assert sha256((ROOT/path).read_bytes()).hexdigest()==digest
    assert 'docs/active-code-map.md' in proposal['source_hashes']
    assert 'docs/nsc-ks-coupled-newton-n32-truncated-iterate.md' not in proposal['source_hashes']
    assert proposal['campaign_best_constraint_maxima']==best['constraint_maxima']
    assert proposal['immediate_parent_constraint_maxima']==parent['constraint_maxima']
    assert proposal['kept_modes']==saved['kept_modes']
    assert proposal['condition_cutoff']==saved['condition_cutoff']
    assert proposal['step_scale']==saved['step_scale']
    assert proposal['step_norm']<=1.+1e-12 and proposal['radius_lower_bound']>0
    fresh=LocalIncomingFamily(np.asarray(proposal['candidate_coefficients'],float))
    assert N32.profile_identity(fresh,include_normal_window=True)==proposal['candidate_profile_identity']
    gradient=np.asarray(parent['action_gradient'],float)
    jacobian=np.asarray(parent['history_jacobian'],float)
    delta,kept,_=N32.truncated_delta(gradient,jacobian,1e6)
    diagnostic=validate_truncated_linear_step(gradient,jacobian,delta,1e6,kept)
    assert diagnostic['physical_error_bound'] is None
    assert diagnostic['physical_local_gate']=='OPEN'
    assert N32.D.beats_measured(np.asarray(proposal['predicted_all_node_residual_maxima']),
                               np.asarray(best['constraint_maxima']))
    assert proposal['prediction_is_not_a_nonlinear_residual'] is True
    bad=delta.copy();bad[0]+=1e-5
    with pytest.raises(ValueError,match='retained coefficient|projected linear|right-space'):
        validate_truncated_linear_step(gradient,jacobian,bad,1e6,kept)
    with pytest.raises(ValueError,match='cutoff coverage'):
        validate_truncated_linear_step(gradient,jacobian,delta,1e6,kept+1)


def _external_raw(names,test_name):
    missing=[]
    for name in names:
        directory=ROOT/'results/development/artifacts'/name
        # A missing array never bypasses the exact descriptor integrity census.
        descriptor_inventory(ROOT,directory,expected_digest=DESCRIPTOR_CENSUS[name])
        missing.extend(raw_input_availability(ROOT,directory)['missing'])
    if missing:
        pytest.skip(f'external_data: {test_name} requires original SHA256-bound campaign arrays; '
                    f'{len(missing)} NPZs unavailable (first: {missing[0]}); '
                    'no numeric replay, fetch, regeneration or metadata substitution')


@pytest.fixture
def external_trail(request):
    _external_raw(TRAIL,request.node.name)


@pytest.fixture
def external_iterate6(request):
    _external_raw((TRAIL[-1],),request.node.name)


@pytest.fixture
def external_n64(request):
    _external_raw(('nsc-ks-coupled-newton-n64-high-modes',),request.node.name)


@pytest.fixture
def external_ruler(request):
    _external_raw((TRAIL[-1],'nsc-ks-coupled-newton-n64-high-modes'),request.node.name)



@pytest.mark.external_data
@pytest.mark.slow
def test_truncated_trail_payloads_match_their_records(external_trail):
    for name in TRAIL:
        directory = ROOT/'results/development/artifacts'/name
        records = sorted(directory.glob('family-*.json'))
        assert len(records) == 60
        first = json.loads(records[0].read_text())
        payload = ROOT/first['payload']['path']
        assert N32.digest(payload) == first['payload']['sha256']
        for path in records:
            record = json.loads(path.read_text())
            assert (ROOT/record['payload']['path']).stat().st_size == record['payload']['bytes']


@pytest.mark.external_data
@pytest.mark.slow
def test_iterate6_replay_keeps_the_measured_residual(external_iterate6):
    result = N32.replay_name('iterate6', family_key=(1, -1),source_context=SOURCE_CONTEXT)
    saved = json.loads((ROOT/'results/development/nsc-ks-coupled-newton-n32-truncated-iterate6.json').read_text())
    assert result['replay_dirac_evolution'] is False
    assert result['constraint_maxima'] == saved['constraint_maxima']
    assert result['profile_identity'] == ITERATE6


@pytest.mark.external_data
@pytest.mark.slow
def test_n64_column_replay_keeps_the_rank_table(external_n64):
    result = N64.replay(family_key=(1, -1),source_context=SOURCE_CONTEXT)
    saved = json.loads(N64.OUTPUT.read_text())
    assert result['replay_dirac_evolution'] is False
    assert result['ranks'] == saved['ranks']
    assert result['profile_identity'] == ITERATE6


@pytest.mark.external_data
@pytest.mark.slow
def test_ruler_drift_is_the_saved_matter_sums(external_ruler):
    result = Ruler.compute(source_context=SOURCE_CONTEXT)
    assert result['families_evolved'] == 0
    assert result['physical_EXISTENCE_certificate'] is False
    assert result['drift_subtracted_as_physics'] is False
    n_drift, beta_drift = result['sum_difference_max_abs']
    assert 5.784e-5 < n_drift < 5.785e-5
    assert 1.284e-10 < beta_drift < 1.285e-10
    assert result['per_family_relative']['median'][0] > 1e-5
    assert result['per_family_relative']['maximum'][0] > 1e-3
    saved = json.loads(Ruler.OUTPUT.read_text())
    assert result == saved


def test_descriptor_census_authenticates_metadata_not_missing_arrays():
    expected_keys=None
    for name,digest in DESCRIPTOR_CENSUS.items():
        inventory=descriptor_inventory(ROOT,ROOT/'results/development/artifacts'/name,
                                       expected_digest=digest)
        assert not inventory['raw_array_integrity_claimed']
        assert len(inventory['rows'])==60
        keys={tuple(row['record']['positive_family']) for row in inventory['rows']}
        if expected_keys is None:expected_keys=keys
        assert keys==expected_keys
        if name in (TRAIL[-1],'nsc-ks-coupled-newton-n64-high-modes'):
            assert {row['record']['profile_identity'] for row in inventory['rows']}=={ITERATE6}
        for row in inventory['rows']:
            if 'physical_local_gate' in row['record']:
                assert row['record']['physical_local_gate']=='OPEN'
            if 'source_and_numerical_error_bound' in row['record']:
                assert row['record']['source_and_numerical_error_bound'] is None


def test_captured_history_keeps_known_metric_without_field_evolution():
    saved=captured_proposal()
    family=LocalIncomingFamily(np.asarray(saved['candidate_coefficients'],float))
    assert N32.profile_identity(family,include_normal_window=True)==ITERATE6
    assert family.radius_lower_bound()>0
    metric=N64.next_metric(family)
    assert len(metric.directions)==65
    assert metric.amplitudes==(1.,)+(0.,)*64
    reference=N64.CompatibleIncomingMetric((1.,),(metric.directions[0],))
    rho=np.array([1.,1.01,1.02])
    assert np.array_equal(metric.values(.15,rho),reference.values(.15,rho))


def test_missing_raw_inputs_are_refused_before_numerical_replay():
    unavailable=raw_input_availability(ROOT,ROOT/'results/development/artifacts'/TRAIL[-1])
    if unavailable['all_available']:
        # The absence contract is still covered by the temporary real-file test.
        return
    for loader in (lambda:N32.replay_name('iterate6'),
                   lambda:N64.replay(family_key=(1,-1)),lambda:Ruler.compute()):
        with pytest.raises(ExternalRawInputsUnavailable,match='external raw campaign inputs unavailable'):
            loader()


def test_raw_input_refusal_never_hides_malformed_or_corrupt_present_bytes(tmp_path):
    directory=tmp_path/'results/development/artifacts/example'
    directory.mkdir(parents=True)
    payload=directory/'family-01-+1.npz'
    descriptor=directory/'family-01-+1.json'
    record={'positive_family':[1,1],'profile_identity':'1'*64,'source_hashes':{'scripts/owner.py':'2'*64},
            'payload':{'path':payload.relative_to(tmp_path).as_posix(),'bytes':1,'sha256':'3'*64}}
    descriptor.write_text(json.dumps(record),encoding='utf8')
    with pytest.raises(ExternalRawInputsUnavailable):
        require_campaign_raw_inputs(tmp_path,directory,expected_count=1)
    bad=json.loads(json.dumps(record));bad['payload']['bytes']=True
    descriptor.write_text(json.dumps(bad),encoding='utf8')
    with pytest.raises(ValueError,match='byte count'):
        require_campaign_raw_inputs(tmp_path,directory,expected_count=1)
    bad=json.loads(json.dumps(record));bad['payload']['sha256']='bad'
    descriptor.write_text(json.dumps(bad),encoding='utf8')
    with pytest.raises(ValueError,match='SHA256'):
        require_campaign_raw_inputs(tmp_path,directory,expected_count=1)
    np.savez(payload,sample=np.array([1.]))
    raw=payload.read_bytes()
    record['payload'].update(bytes=len(raw),sha256=sha256(raw).hexdigest())
    descriptor.write_text(json.dumps(record),encoding='utf8')
    assert require_campaign_raw_inputs(tmp_path,directory,expected_count=1)['all_available']
    payload.write_bytes(raw[:-1]+bytes([raw[-1]^1]))
    with pytest.raises(ValueError,match='SHA256 changed'):
        require_campaign_raw_inputs(tmp_path,directory,expected_count=1)


def test_unavailable_inherited_source_is_disclosed_and_strict_replay_refuses():
    saved=captured_proposal()
    provenance=stage0_capture_source_provenance(ROOT,saved['source_hashes'])
    unresolved=[row for row in provenance['bindings'] if not row['authenticated']]
    assert unresolved==[{
        'path':'src/recursive_horizons/nsc_ks_history_evaluator.py',
        'expected_sha256':'86ed312db133f201a3cce11cacc96b62dc89d6c3fe434bdac0f77ca54d023cd6',
        'authenticated':False,
        'reason':'declared inherited history-evaluator bytes unavailable in bounded preserved source history',
    }]
    assert provenance['creating_CLI_authenticated'] is True
    assert provenance['inherited_math_authenticated'] is False
    with pytest.raises(ValueError,match='inherited source unavailable'):
        require_original_campaign_sources(ROOT,{unresolved[0]['path']:unresolved[0]['expected_sha256']},
                                          source_context=SOURCE_CONTEXT)
    bad=dict(saved['source_hashes']);bad[unresolved[0]['path']]='0'*64
    with pytest.raises(ValueError,match='SHA256 changed'):
        stage0_capture_source_provenance(ROOT,bad)
    bad=dict(saved['source_hashes']);bad['scripts/derive_nsc_ks_coupled_newton_n32_truncated_iterate.py']='0'*64
    with pytest.raises(RuntimeError,match='source drift'):
        stage0_capture_source_provenance(ROOT,bad)


def test_finite_svd_control_rejects_coerced_types_shapes_and_nonfinite_operands():
    gradient=np.array([[.1,.2],[.3,.4]])
    jacobian=np.array([[[1.,.5],[.2,.3]],[[.2,.4],[.7,.1]]])
    delta,kept,_=N32.truncated_delta(gradient,jacobian,1e6)
    validate_truncated_linear_step(gradient,jacobian,delta,1e6,kept)
    for bad in (gradient.astype(str),gradient.astype(bool),gradient.astype(complex),
                np.array([[np.nan,.2],[.3,.4]]),np.array([[np.inf,.2],[.3,.4]])):
        with pytest.raises(ValueError,match='real numeric'):
            validate_truncated_linear_step(bad,jacobian,delta,1e6,kept)
    for bad in ('1000000',True,1e6+0j,np.nan,np.inf,[1e6]):
        with pytest.raises(ValueError,match='real numeric'):
            validate_truncated_linear_step(gradient,jacobian,delta,bad,kept)
    for bad in (delta.astype(str),delta.astype(bool),delta.astype(complex),delta[:,None]):
        with pytest.raises(ValueError,match='numeric|layout'):
            validate_truncated_linear_step(gradient,jacobian,bad,1e6,kept)
    for bad in (True,0,len(delta)+1,float(kept)):
        with pytest.raises(ValueError,match='layout'):
            validate_truncated_linear_step(gradient,jacobian,delta,1e6,bad)
