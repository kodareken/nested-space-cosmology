"""Only two exact creator-label replacements can bridge original receipts."""
from copy import deepcopy
import hashlib
from pathlib import Path
import re

import pytest

from recursive_horizons import nsc_threshold18_metadata as metadata
from recursive_horizons.nsc_frozen_sources import (
    frozen_source_context, selected_source_sha256)

ROOT=Path(__file__).resolve().parents[1]
PATHS=tuple(metadata.METADATA_IDENTITIES)


@pytest.fixture(scope='module')
def pairs():
    return {path:metadata._authenticated_pair(
        ROOT,path,identity['original_sha256'])
        for path,identity in metadata.METADATA_IDENTITIES.items()}


@pytest.mark.parametrize('path',PATHS)
def test_actual_local_metadata_has_separate_exact_identities(path,pairs):
    original,current=pairs[path]
    identity=metadata.METADATA_IDENTITIES[path]
    receipt=metadata.authenticate_normalized_metadata(ROOT,path,identity['original_sha256'])
    assert receipt['original_sha256']==hashlib.sha256(original).hexdigest()
    assert receipt['current_sha256']==hashlib.sha256(current).hexdigest()
    assert receipt['original_sha256']!=receipt['current_sha256']
    assert receipt['changed_pointer']=='/mac_integration/lab_root'
    assert receipt['current_label']=='.'
    with frozen_source_context(metadata.AUTHENTICATION_COMMIT):
        assert selected_source_sha256(ROOT,path)==receipt['current_sha256']


@pytest.mark.parametrize('kind',('numeric','digest','different_label','wrong_root'))
@pytest.mark.parametrize('path',PATHS)
def test_refreshed_current_identity_cannot_authorize_other_changes(path,kind,pairs):
    original,current=pairs[path]
    if kind=='numeric':
        changed=current.replace(b'"angular_sign": 1',b'"angular_sign": -1',1)
    elif kind=='digest':
        match=re.search(rb'"([0-9a-f]{64})"',current)
        assert match is not None
        start=match.start(1)
        changed=current[:start]+(b'0' if current[start:start+1]!=b'0' else b'1')+current[start+1:]
    elif kind=='different_label':
        changed=current.replace(b'"panel": "group14/low32_1"',b'"panel": "changed"',1)
    else:
        changed=current.replace(b'"lab_root": "."',b'"lab_root": "./"',1)
    assert changed!=current
    refreshed=deepcopy(metadata.METADATA_IDENTITIES[path])
    refreshed.update(current_sha256=hashlib.sha256(changed).hexdigest(),current_bytes=len(changed))
    with pytest.raises(ValueError,match='beyond the authorized path label'):
        metadata._prove_pair(original,changed,refreshed)


@pytest.mark.parametrize('path',PATHS)
def test_missing_local_json_never_uses_original_as_a_fallback(tmp_path,monkeypatch,path):
    def forbidden(*_):
        pytest.fail('missing local numerical metadata reached Git')
    monkeypatch.setattr(metadata,'_original_metadata_bytes',forbidden)
    with pytest.raises(ValueError,match='local regular file'):
        metadata.authenticate_normalized_metadata(
            tmp_path,path,metadata.METADATA_IDENTITIES[path]['original_sha256'])


def test_other_paths_and_wrong_original_digests_are_rejected():
    with pytest.raises(ValueError,match='exact original path and digest'):
        metadata.authenticate_normalized_metadata(ROOT,PATHS[0],'0'*64)
    with pytest.raises(ValueError,match='exact original path and digest'):
        metadata.authenticate_normalized_metadata(ROOT,'results/other.json','0'*64)


def test_current_and_original_25_file_tree_pins_are_both_authenticated():
    receipt=metadata.authenticate_import_tree(ROOT)
    assert receipt['files']==25
    assert receipt['original_sha256']==metadata.ORIGINAL_IMPORT_TREE_SHA256
    assert receipt['current_sha256']==metadata.CURRENT_IMPORT_TREE_SHA256
    assert receipt['original_sha256']!=receipt['current_sha256']
    assert {row['path'] for row in receipt['metadata_normalizations']}==set(PATHS)
