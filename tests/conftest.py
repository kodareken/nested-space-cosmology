"""One root package, with explicit Windows limitations for historical tests."""
import os
from pathlib import Path
import sys
import warnings

import pytest

ROOT = Path(__file__).resolve().parents[1]
for _path in (str(ROOT), str(ROOT / "src"), str(ROOT / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from recursive_horizons import nsc_posix  # noqa: E402

_installed = nsc_posix.install()
if _installed:
    warnings.warn(
        "POSIX stand-ins registered for: " + ", ".join(_installed)
        + ". CPU budget guards do not fire and advisory run locks are not"
        " enforced on this platform; see recursive_horizons.nsc_posix.",
        RuntimeWarning,
        stacklevel=1,
    )


def pytest_collection_modifyitems(items):
    if os.name != 'nt':
        return
    legacy = ('tests/test_nsc_pg_retained_tail.py::'
              'test_retained_tail_record_is_tail_only_and_covers_all_families')
    for item in items:
        if item.nodeid.replace('\\', '/') == legacy:
            item.add_marker(pytest.mark.skip(reason=(
                'Historical exact replay serializes native path separators; '
                'test_portable_tail_replay authenticates every field on Windows. '
                'The unchanged exact test runs on Linux.')))
