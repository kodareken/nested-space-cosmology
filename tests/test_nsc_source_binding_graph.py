"""Original AST source authentication is explicit; outputs stay local."""
from pathlib import Path
from hashlib import sha256

import pytest

from recursive_horizons.nsc_ks_evaluation_binding import (
    LEGACY_REPLAY_MODE, implementation_hashes, resolve_recorded_source_bytes)
from recursive_horizons.nsc_frozen_sources import source_sha256

ROOT = Path(__file__).resolve().parents[1]
CONTEXT = '329ef809e466a86ed4a290ee32c19485f85add9d'


def test_original_graph_authenticates_retired_sources_with_explicit_pin():
    owners = ('scripts/derive_nsc_evolved_incoming_state.py', 'docs/nsc-evolved-incoming-state.md')
    hashes = implementation_hashes(ROOT, owners, source_context=CONTEXT)
    assert hashes[owners[0]] == 'f33ee6bb7854feb424aefcd5063105f6c1222ff082a570c70ff36b841d411939'
    assert hashes[owners[1]] == source_sha256(ROOT, owners[1], context=CONTEXT)
    assert 'src/recursive_horizons/nsc_evolved_incoming_state.py' in hashes
    assert all(not name.startswith('results/') for name in hashes)
    with pytest.raises(ValueError, match='missing binding field'):
        implementation_hashes(ROOT, owners)


def test_original_graph_rejects_output_history_and_abbreviated_context():
    with pytest.raises(ValueError, match='implementation source'):
        implementation_hashes(ROOT, ('results/nsc-3-boundary-response.json',), source_context=CONTEXT)
    with pytest.raises(ValueError, match='full'):
        implementation_hashes(ROOT, ('docs/nsc-evolved-incoming-state.md',), source_context='329ef809')


@pytest.mark.parametrize('name', (
    'docs/numeric-audit.json', 'paper/evidence/receipt.json',
    'paper/figure-evidence/calculation.md',
))
def test_legacy_source_resolution_never_loads_numeric_or_evidence_history(tmp_path, name):
    raw = b'original numerical receipt\n'
    digest = sha256(raw).hexdigest()

    def forbidden(*args, **kwargs):
        raise AssertionError('historical output lookup attempted')

    with pytest.raises(ValueError, match='recorded data hash mismatch'):
        resolve_recorded_source_bytes(
            tmp_path, name, digest, mode=LEGACY_REPLAY_MODE,
            git_commit=CONTEXT, git_show=forbidden, git_log=forbidden)
    current = tmp_path/name
    current.parent.mkdir(parents=True, exist_ok=True)
    current.write_bytes(raw)
    assert resolve_recorded_source_bytes(
        tmp_path, name, digest, mode=LEGACY_REPLAY_MODE,
        git_commit=CONTEXT, git_show=forbidden, git_log=forbidden,
    ) == (raw, 'working-tree')
