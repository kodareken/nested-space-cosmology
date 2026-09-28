"""A later clipped n=16 trial replay preserves incomplete operator artifacts."""
import importlib
from pathlib import Path

import pytest


@pytest.mark.parametrize('mode', ('assemble', 'run'))
def test_missing_receipt_never_deletes_an_existing_operator(tmp_path, monkeypatch, mode):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/'scripts'))
    module = importlib.import_module('derive_nsc_ks_coupled_newton_n16_damped_iterate')
    module.bind()
    monkeypatch.setattr(module.N16, 'DIRECTORY', tmp_path)
    context = {'archived': {(11, 1): {}}, 'target': (0., 1., 2.)}
    monkeypatch.setattr(module.N16, 'context', lambda: context)
    monkeypatch.setattr(module, 'context', lambda: context)
    _, path = module.N16.family_paths((11, 1))
    original = b'owned clipped n=16 iterate operator artifact awaiting receipt recovery'
    path.write_bytes(original)
    with pytest.raises(ValueError, match='retained for recovery'):
        if mode == 'assemble':
            module.assemble(context, replay=True)
        else:
            module.N16.run(5., 0)
    assert path.read_bytes() == original
    assert list(tmp_path.iterdir()) == [path]
