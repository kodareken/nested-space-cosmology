"""Authenticate the two authorized threshold18 creator path-label changes.

Only the current local JSON is a computational input. Original Git bytes
authenticate the one string-token replacement and the original tree identity;
they are never substituted for missing or changed numerical evidence.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess


AUTHENTICATION_COMMIT = 'f2671a8c4c72add006ba746a555a601a311b9629'
IMPORT_RELATIVE = 'results/development/nsc-threshold18-import-v1'
ORIGINAL_IMPORT_TREE_SHA256 = 'bc7c5c2424810650acbeae20c4bedfc18548737f112725c63702eb3f9920a4c7'
CURRENT_IMPORT_TREE_SHA256 = '593df76bf8f85278069f5d2e26fa2329aac821b83bf1cc10aa914f61155973df'
METADATA_IDENTITIES = {
    IMPORT_RELATIVE+'/evidence/actual18-replay.json': {
        'original_sha256':'1708ac804c931b594b64bf13842ec6d66325dacc1543b7c78f0d879f702b5405',
        'original_bytes':59095,
        'current_sha256':'7c8f61de68c40f6e11e2f3f8285999fa04bc699a1ade74fa8af74fbad90b625b',
        'current_bytes':59050,
    },
    IMPORT_RELATIVE+'/evidence/actual18-report.json': {
        'original_sha256':'b65ccf7752db441c635f4ea6f8b5efd7faed735048cef1d3c3155c1bbad12601',
        'original_bytes':59103,
        'current_sha256':'4714f949fd9c4cbe2bb2383aa4de444169bd64602ee145c2b4c70262b9ce13d9',
        'current_bytes':59058,
    },
}


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _original_metadata_bytes(root, relative):
    """Authentication-only read, restricted to the fixed two metadata paths."""
    if relative not in METADATA_IDENTITIES:
        raise ValueError('metadata normalization path is not authorized')
    environment=os.environ.copy()
    for key in ('GIT_DIR','GIT_WORK_TREE','GIT_INDEX_FILE','GIT_COMMON_DIR'):
        environment.pop(key,None)
    try:
        return subprocess.check_output(
            ['git','--no-pager','-C',str(root),'show','--no-textconv',
             AUTHENTICATION_COMMIT+':'+relative],
            env=environment,stderr=subprocess.PIPE,timeout=30)
    except (OSError,subprocess.SubprocessError) as error:
        raise ValueError('original threshold metadata authentication unavailable') from error


def _prove_pair(original, current, identity):
    """Exact transformation proof, independently of any current-file manifest."""
    if len(current)!=identity['current_bytes'] or _sha(current)!=identity['current_sha256']:
        raise ValueError('current normalized threshold metadata bytes changed')
    if len(original)!=identity['original_bytes'] or _sha(original)!=identity['original_sha256']:
        raise ValueError('original threshold metadata authentication changed')
    try:
        label=json.loads(original)['mac_integration']['lab_root']
    except (ValueError,KeyError,TypeError) as error:
        raise ValueError('original threshold path-label metadata malformed') from error
    if not isinstance(label,str) or label=='.':
        raise ValueError('original threshold path label must be a distinct string')
    token=json.dumps(label).encode('utf-8')
    # Byte comparison preserves numerical spellings, digest values, ordering,
    # whitespace and every other label. No JSON reserialization is accepted.
    if original.count(token)!=1 or original.replace(token,b'"."')!=current:
        raise ValueError('threshold metadata differs beyond the authorized path label')


def _authenticated_pair(root, relative, expected_original_sha256):
    root=Path(root).resolve()
    identity=METADATA_IDENTITIES.get(relative)
    if identity is None or expected_original_sha256!=identity['original_sha256']:
        raise ValueError('metadata normalization requires the exact original path and digest')
    path=root/relative
    if path.is_symlink() or not path.resolve().is_relative_to(root) or not path.is_file():
        raise ValueError('current threshold metadata must be a local regular file')
    current=path.read_bytes()
    if len(current)!=identity['current_bytes'] or _sha(current)!=identity['current_sha256']:
        raise ValueError('current normalized threshold metadata bytes changed')
    original=_original_metadata_bytes(root,relative)
    _prove_pair(original,current,identity)
    return original,current


def authenticate_normalized_metadata(root, relative, expected_original_sha256):
    """Return separate original/current identities after exact transform proof."""
    original,current=_authenticated_pair(root,relative,expected_original_sha256)
    return {'path':relative,'original_sha256':_sha(original),'current_sha256':_sha(current),
            'original_bytes':len(original),'current_bytes':len(current),
            'authentication_commit':AUTHENTICATION_COMMIT,
            'changed_pointer':'/mac_integration/lab_root','current_label':'.',
            'scope':'Single creator path-label token; numerical and digest bytes unchanged.'}


def authenticate_import_tree(root):
    """Authenticate current local leaves and the separately named original tree."""
    root=Path(root).resolve()
    folder=root/IMPORT_RELATIVE
    if folder.is_symlink() or not folder.resolve().is_relative_to(root):
        raise ValueError('threshold tree must stay inside the current repository')
    originals={}
    normalizations=[]
    for relative,identity in METADATA_IDENTITIES.items():
        original,_current=_authenticated_pair(root,relative,identity['original_sha256'])
        originals[relative]=original
        normalizations.append({
            **identity,'path':relative,'authentication_commit':AUTHENTICATION_COMMIT,
            'changed_pointer':'/mac_integration/lab_root','current_label':'.',
            'scope':'Single creator path-label token; numerical and digest bytes unchanged.'})
    files=sorted(path for path in folder.rglob('*') if path.is_file())
    current_digest=hashlib.sha256(); original_digest=hashlib.sha256()
    for path in files:
        if path.is_symlink() or not path.resolve().is_relative_to(root):
            raise ValueError('threshold tree leaf must be a local regular file')
        prefix=path.relative_to(folder).as_posix().encode()+b'\0'
        current=path.read_bytes()
        original=originals.get(path.relative_to(root).as_posix(),current)
        for digest,raw in ((current_digest,current),(original_digest,original)):
            digest.update(prefix);digest.update(raw);digest.update(b'\0')
    if (len(files)!=25 or current_digest.hexdigest()!=CURRENT_IMPORT_TREE_SHA256
            or original_digest.hexdigest()!=ORIGINAL_IMPORT_TREE_SHA256):
        raise ValueError('threshold original or current import tree changed')
    return {'original_sha256':original_digest.hexdigest(),
            'current_sha256':current_digest.hexdigest(),'files':len(files),
            'metadata_normalizations':normalizations}
