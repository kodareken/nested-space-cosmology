"""Exact authentication of the bounded prepared-history metadata rebinding.

This API reports original and current carrier identities separately. It never
changes a current-file digest or supplies historical data to a computation.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path


LEDGER_PATH='results/metadata/prepared-history/receipt-rebind.json'
LEDGER_SHA256='46857b7ba2f51b53038305d3796b4de03d4e7a2551e167e20f9100976cc3ceeb'
BACKEND_PATH='scripts/rebind_nsc_prepared_history_receipts.py'
BACKEND_SHA256='3ba886fe408e5758c3b66beff06e83fb4d67cc07cf377a3f97fcf732d88aefaa'


def _pinned_file(root,relative,expected):
    path=root/relative
    if path.is_symlink() or not path.resolve().is_relative_to(root) or not path.is_file():
        raise ValueError('rebind authentication file must be local: '+relative)
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=expected:
        raise ValueError('rebind authentication file changed: '+relative)
    return path,raw


def _maintenance(root):
    path,_raw=_pinned_file(root,BACKEND_PATH,BACKEND_SHA256)
    spec=importlib.util.spec_from_file_location('nsc_prepared_history_metadata_backend',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def authenticate_rebound_receipt(root,path,expected_original_sha):
    """Prove one allowlisted old/current metadata-only identity relation.

    The fixed ledger and fixed baseline originals authenticate the complete
    topological transformation. Actual current files must equal the rebuilt
    bytes, including unchanged numeric spellings, source captures and NPY
    members. Original Git bytes serve authentication only, never fallback.
    """
    root=Path(root).resolve()
    _ledger,raw=_pinned_file(root,LEDGER_PATH,LEDGER_SHA256)
    backend=_maintenance(root)
    if not isinstance(path,str) or backend.ORIGINALS.get(path)!=expected_original_sha:
        raise ValueError('exact allowlisted original receipt identity required')
    audit=json.loads(raw)
    verified=backend.verify_audit(root,audit)
    row=next(row for row in verified['rows'] if row['path']==path)
    return {'path':path,'original_sha256':row['original']['sha256'],
            'current_sha256':row['current']['sha256'],
            'original_bytes':row['original']['bytes'],'current_bytes':row['current']['bytes'],
            'ledger_sha256':LEDGER_SHA256,'original_source_context':verified['original_source_context'],
            'scope':'Exact bounded metadata rebinding; scientific tokens, execution capture and nonmetadata NPY members unchanged.'}
