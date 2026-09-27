"""Keep authenticated historical tests unchanged while routing portable replay."""
import os
import pytest


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
