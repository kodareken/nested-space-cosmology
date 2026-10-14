"""Metadata maintenance cannot change science, source captures or array bits."""
import io
import json
import zipfile

import numpy as np
import pytest

import rebind_nsc_prepared_history_receipts as maintenance

PATH='results/development/input.json'
OLD='a'*64;NEW='b'*64
ORIGINAL={PATH:{'sha256':OLD,'bytes':17}}
CURRENT={PATH:{'sha256':NEW,'bytes':17}}


def npy(value):
    stream=io.BytesIO();np.save(stream,value,allow_pickle=False);return stream.getvalue()


def test_signed_zero_arrays_and_executed_source_members_stay_byte_identical():
    meta=json.dumps({'signature':{PATH:OLD},'verification_signature':{PATH:OLD},
                     'execution_signature':{PATH:OLD},'CPU_seconds':1.234567}).encode()
    members={'signed.npy':npy(np.array([-0.,0.,-1.])),
             'executed_sources/producer.py.npy':npy(np.frombuffer(b'original capture\n',np.uint8)),
             'metadata_json.npy':npy(np.frombuffer(meta,np.uint8))}
    source=io.BytesIO()
    with zipfile.ZipFile(source,'w') as archive:
        for name,raw in members.items():archive.writestr(name,raw)
    rebound,changes,preserved=maintenance.rewrite_npz(source.getvalue(),ORIGINAL,CURRENT)
    with zipfile.ZipFile(io.BytesIO(rebound)) as archive:
        for name in preserved:assert archive.read(name)==members[name]
        signed=np.lib.format.read_array(io.BytesIO(archive.read('signed.npy')),allow_pickle=False)
        assert np.signbit(signed[0]) and not np.signbit(signed[1])
    actual=json.loads(maintenance.metadata_json(rebound))
    assert actual['signature'][PATH]==actual['verification_signature'][PATH]==NEW
    assert actual['execution_signature'][PATH]==OLD
    assert actual['CPU_seconds']==1.234567
    assert len(changes)==2


@pytest.mark.parametrize('raw',[
    json.dumps({'input_hashes':{PATH:'c'*64}}).encode(),
    json.dumps({'unqualified':OLD}).encode(),
])
def test_wrong_or_unqualified_original_digest_is_rejected(raw):
    with pytest.raises(ValueError,match='original receipt digest changed|unqualified live digest'):
        maintenance.rewrite_json(raw,ORIGINAL,CURRENT)


def test_scientific_numeric_token_spellings_are_preserved_and_mutations_rejected():
    raw=(' {"input_hashes":{"'+PATH+'":"'+OLD+'"},'
         '"measurement":1.2300e-12,"negative_zero":-0.0,"count":3}\n').encode()
    rebound,changes=maintenance.rewrite_json(raw,ORIGINAL,CURRENT)
    assert b'1.2300e-12' in rebound and b'-0.0' in rebound
    damaged=rebound.replace(b'1.2300e-12',b'1.2400e-12')
    with pytest.raises(ValueError,match='nonmetadata receipt token changed'):
        maintenance.validate_json_changes(raw,damaged,changes)
    changed_type=rebound.replace(b'"count":3',b'"count":true')
    with pytest.raises(ValueError,match='receipt token type changed'):
        maintenance.validate_json_changes(raw,changed_type,changes)


@pytest.mark.parametrize('size',['17',True,18])
def test_payload_size_cannot_be_coerced_or_refresh_wrong_original_identity(size):
    raw=json.dumps({'payload':{'path':PATH,'sha256':OLD,'bytes':size}}).encode()
    with pytest.raises(ValueError,match='original payload size changed'):
        maintenance.rewrite_json(raw,ORIGINAL,CURRENT)


def test_duplicate_keys_cannot_hide_a_changed_receipt():
    with pytest.raises(ValueError,match='ambiguous JSON key'):
        maintenance.json_tokens(b'{"value":1,"value":2}')


def test_numeric_drift_cannot_self_approve_by_refreshing_ledger_hash(tmp_path,monkeypatch):
    """A synthetic authenticated original isolates the verifier's trust boundary."""
    from copy import deepcopy
    seed=b'fixed corrected metadata seed'
    name='results/development/synthetic-live.json'
    original=(' {"input_hashes":{"'+maintenance.SEED+'":"'+OLD+'"},'
              '"measurement":1.2500e-12,"count":3}\n').encode()
    monkeypatch.setattr(maintenance,'SEED_OLD',OLD)
    monkeypatch.setattr(maintenance,'SEED_CURRENT',maintenance.sha(seed))
    monkeypatch.setattr(maintenance,'ORIGINALS',{name:maintenance.sha(original)})
    monkeypatch.setattr(maintenance,'authenticated_originals',lambda root:{name:original})
    outputs,audit=maintenance.build_plan({name:original},seed)
    seed_path=tmp_path/maintenance.SEED;seed_path.parent.mkdir(parents=True);seed_path.write_bytes(seed)
    path=tmp_path/name;path.write_bytes(outputs[name])
    maintenance.verify_audit(tmp_path,audit)
    changed=outputs[name].replace(b'1.2500e-12',b'9.9990e-12')
    assert changed!=outputs[name]
    path.write_bytes(changed)
    refreshed=deepcopy(audit)
    refreshed['rows'][0]['current'].update(sha256=maintenance.sha(changed),bytes=len(changed))
    with pytest.raises(ValueError,match='audit differs from authenticated transformation'):
        maintenance.verify_audit(tmp_path,refreshed)
    with pytest.raises(ValueError,match='current bytes differ from authenticated transformation'):
        maintenance.verify_audit(tmp_path,audit)
