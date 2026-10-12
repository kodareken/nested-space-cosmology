"""A scientific replay preserves incomplete operator artifacts for recovery."""
import importlib
from pathlib import Path
import sys

import pytest


@pytest.mark.parametrize('mode', ('assemble', 'run'))
def test_missing_receipt_never_deletes_an_existing_operator(tmp_path, monkeypatch, mode):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/'scripts'))
    module = importlib.import_module('derive_nsc_ks_coupled_retained_control')
    monkeypatch.setattr(module, 'DIRECTORY', tmp_path)
    context = {'archived': {(14, 1): {}}, 'target': (0., 1., 2.)}
    monkeypatch.setattr(module, 'context', lambda: context)
    _, path = module.family_paths((14, 1))
    original = b'owned operator artifact awaiting receipt recovery'
    path.write_bytes(original)
    with pytest.raises(ValueError, match='retained for recovery'):
        if mode == 'assemble':
            module.assemble(context, replay=True)
        else:
            module.run(5., 0)
    assert path.read_bytes() == original
    assert list(tmp_path.iterdir()) == [path]
