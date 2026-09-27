"""Windows all-field replay with only the known path spelling normalized."""
import os
from pathlib import PureWindowsPath
import pytest

from test_nsc_pg_retained_tail import ROOT, _loaded, make_record, verify_arrays


@pytest.mark.skipif(os.name != 'nt', reason='Linux runs the original exact replay')
def test_windows_tail_replay_preserves_every_scientific_field():
    record, arrays, metadata, state = _loaded()
    verify_arrays(arrays, metadata, state)
    rebuilt, _ = make_record(ROOT / record['payload']['path'], arrays, metadata)
    # Do not normalize numbers, hashes, other paths, or original record bytes.
    actual = PureWindowsPath(rebuilt['payload']['path']).as_posix()
    assert actual == record['payload']['path']
    rebuilt['payload']['path'] = actual
    assert rebuilt == record
