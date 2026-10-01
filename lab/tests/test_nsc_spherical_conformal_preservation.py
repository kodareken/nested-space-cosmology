"""Sealed outputs refuse recomputation; explicit historical replay stays read-only."""
from copy import deepcopy
import json

import pytest

import derive_nsc_spherical_conformal_episode as episode
import derive_nsc_spherical_conformal_frames as frames
import derive_nsc_spherical_conformal_continuation as continuation
import derive_nsc_spherical_conformal_transport as transport
import derive_nsc_spherical_conformal_clock as clock

OWNERS = (episode, frames, continuation, transport, clock)


@pytest.mark.parametrize("owner", OWNERS)
def test_existing_output_refuses_before_source_read_or_compute(owner, monkeypatch, tmp_path):
    path = tmp_path / "sealed.json"
    path.write_bytes(b"immutable sentinel")
    monkeypatch.setattr(owner, "OUT", path)
    with pytest.raises(FileExistsError, match="sealed conformal output"):
        owner.run()
    assert path.read_bytes() == b"immutable sentinel"


def test_existing_payload_alone_also_refuses_compute(monkeypatch, tmp_path):
    payload = tmp_path / "sealed.npz"
    payload.write_bytes(b"partial checkpoint")
    monkeypatch.setattr(episode, "OUT", tmp_path / "missing.json")
    monkeypatch.setattr(episode, "NPZ", payload)
    with pytest.raises(FileExistsError):
        episode.run()
    assert payload.read_bytes() == b"partial checkpoint"


@pytest.mark.parametrize("owner", OWNERS)
def test_historical_check_authenticates_sources_without_writing(owner):
    paths = [owner.OUT]
    if hasattr(owner, "NPZ"):
        paths.append(owner.NPZ)
    before = {path: episode.sha256(path) for path in paths}
    report = owner.verify_saved(source_ref=episode.SEALED_SOURCE_REF)
    assert report["wrote"] is False
    assert report["source_ref"] == episode.SEALED_SOURCE_REF
    assert before == {path: episode.sha256(path) for path in paths}


def test_current_requests_do_not_substitute_history_and_replay_rejects_mutations():
    record = json.loads(episode.OUT.read_text())
    paths = episode.source_paths()
    with pytest.raises(ValueError, match="recorded source changed: driver"):
        episode.check_recorded_sources(record["source_hashes_before"], paths)
    bad = deepcopy(record["source_hashes_before"])
    bad["driver"] = "0" * 64
    with pytest.raises(RuntimeError, match="pinned source drift"):
        episode.check_recorded_sources(bad, paths, source_ref=episode.SEALED_SOURCE_REF)
    del bad["driver"]
    with pytest.raises(ValueError, match="inventory"):
        episode.check_recorded_sources(bad, paths, source_ref=episode.SEALED_SOURCE_REF)
    with pytest.raises(RuntimeError, match="unavailable"):
        episode.check_recorded_sources(record["source_hashes_before"], paths, source_ref="0" * 40)
    numerical = deepcopy(record["source_hashes_before"])
    numerical["v5_npz"] = "0" * 64
    with pytest.raises(ValueError, match="v5_npz"):
        episode.check_recorded_sources(numerical, paths, source_ref=episode.SEALED_SOURCE_REF)
    with pytest.raises(ValueError, match="numerical value"):
        episode.check_saved_values(0.03937795, 0.04, "surface effect")
