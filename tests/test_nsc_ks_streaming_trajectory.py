"""Tiny synthetic streaming capture, chunk replay, and tamper/incomplete checks."""
import gc
from hashlib import sha256
import json
from pathlib import Path
import weakref

import numpy as np
import pytest

from recursive_horizons.nsc_ks_difference_envelope import evolve_ks_difference_envelope
from recursive_horizons.nsc_ks_source_envelope import (
    computational_z_grid, physical_incoming_interval, usual_axial_support,
)
from recursive_horizons.nsc_ks_streaming_trajectory import (
    MAX_SEGMENTS_PER_CHUNK, SCHEMA, TrajectoryChunkReader, TrajectoryChunkWriter,
    chunk_filename, stream_ks_trajectory,
)
from recursive_horizons.nsc_ks_trajectory import capture_ks_trajectory
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from test_nsc_ks_difference_envelope import inputs


def tiny():
    args = (*inputs(), computational_z_grid(16), np.linspace(*physical_incoming_interval(), 5),
            1.5, np.sqrt(5), 1.03)
    options = dict(axial_support=usual_axial_support(), rtol=2e-12, atol=2e-14,
                   max_step=.001, tangents='zero')
    return args, options


def test_stream_matches_direct_tiny_evolution_without_patching_globals():
    args, options = tiny()
    source = args[0]
    before = evolve_ks_difference_envelope.__globals__['DOP853']
    original = evolve_ks_difference_envelope(*args, **options)
    observed = []
    prepared, summary = stream_ks_trajectory(*args, on_segment=observed.append, **options)
    captured, segments = capture_ks_trajectory(*args, **options)
    assert evolve_ks_difference_envelope.__globals__['DOP853'] is before
    for name in ('columns', 'axial_columns', 'reference_amplitudes', 'envelope_difference',
                 'source_covariance', 'column_weights'):
        np.testing.assert_array_equal(getattr(original, name), getattr(prepared, name))
        np.testing.assert_array_equal(getattr(original, name), getattr(captured, name))
    np.testing.assert_array_equal(prepared.source_covariance, source.covariance)
    np.testing.assert_array_equal(prepared.column_weights, source.column_weights)
    assert prepared.fixed_preparation_digest == original.fixed_preparation_digest
    assert summary.accepted_steps == original.diagnostics['accepted_steps'] == len(observed)
    assert summary.rho_start == prepared.rho_up and summary.rho_end == 1.
    assert len(observed) == len(segments)
    for streamed, stored in zip(observed, segments):
        assert streamed.rho_start == stored.rho_start and streamed.rho_end == stored.rho_end
        np.testing.assert_array_equal(streamed.start, stored.start)
        np.testing.assert_array_equal(streamed.end, stored.end)
        np.testing.assert_array_equal(streamed.coefficients, stored.coefficients)


def test_stream_does_not_retain_the_trajectory():
    args, options = tiny()
    refs = []

    def observe(segment):
        refs.append(weakref.ref(segment))

    prepared, summary = stream_ks_trajectory(*args, on_segment=observe, **options)
    del prepared
    gc.collect()
    assert summary.accepted_steps == len(refs)
    assert not hasattr(summary, 'segments')
    assert all(reference() is None for reference in refs)


def test_chunk_bytes_are_deterministic_and_shared_endpoints_are_bitwise(tmp_path):
    args, options = tiny()
    first, second = tmp_path / 'a', tmp_path / 'b'
    hashes = []
    for directory in (first, second):
        with TrajectoryChunkWriter(directory) as writer:
            prepared, summary = stream_ks_trajectory(*args, on_segment=writer, **options)
        assert writer.buffered_segments == 0
        reader = TrajectoryChunkReader(directory)
        assert reader.manifest['schema'] == SCHEMA
        assert reader.manifest['complete'] is True
        assert reader.manifest['accepted_steps'] == summary.accepted_steps
        assert reader.manifest['max_segments_per_chunk'] == MAX_SEGMENTS_PER_CHUNK
        hashes.append(tuple(entry['sha256'] for entry in reader.manifest['chunks']))
        previous = None
        total = 0
        for chunk in reader:
            assert 1 <= len(chunk) <= MAX_SEGMENTS_PER_CHUNK
            total += len(chunk)
            if previous is not None:
                assert previous.rho_nodes[-1].tobytes() == chunk.rho_nodes[0].tobytes()
                assert previous.state_nodes[-1].tobytes() == chunk.state_nodes[0].tobytes()
            previous = chunk
        assert total == summary.accepted_steps
        np.testing.assert_array_equal(prepared.source_covariance, args[0].covariance)
        np.testing.assert_array_equal(prepared.column_weights, args[0].column_weights)
    assert hashes[0] == hashes[1] and hashes[0]
    with pytest.raises(FileExistsError, match='already exists'):
        TrajectoryChunkWriter(first)


def test_reader_rejects_tampered_hash_continuity_and_incomplete_runs(tmp_path):
    args, options = tiny()
    complete = tmp_path / 'complete'
    with TrajectoryChunkWriter(complete) as writer:
        stream_ks_trajectory(*args, on_segment=writer, **options)
    interrupted = tmp_path / 'interrupted'
    writer = TrajectoryChunkWriter(interrupted)
    count = 0

    def stop_after_one_chunk(segment):
        nonlocal count
        writer(segment)
        count += 1
        if count == MAX_SEGMENTS_PER_CHUNK:
            raise RuntimeError('interrupt after a complete chunk')

    with pytest.raises(RuntimeError, match='interrupt'):
        stream_ks_trajectory(*args, on_segment=stop_after_one_chunk, **options)
    writer.close(complete=False)
    with pytest.raises(ValueError, match='incomplete'):
        TrajectoryChunkReader(interrupted)
    recovered = TrajectoryChunkReader(interrupted, require_complete=False)
    assert recovered.manifest['complete'] is False
    chunks = list(recovered)
    assert len(chunks) == 1 and len(chunks[0]) == MAX_SEGMENTS_PER_CHUNK
    with pytest.raises(ValueError, match='cannot be labeled complete'):
        writer.close(complete=True)

    tampered = tmp_path / 'tampered'
    _copy_run(complete, tampered)
    target = tampered / chunk_filename(0)
    data = bytearray(target.read_bytes())
    data[-20] ^= 0xFF
    target.write_bytes(bytes(data))
    with pytest.raises(ValueError, match='sha256'):
        list(TrajectoryChunkReader(tampered))

    broken = tmp_path / 'broken'
    _copy_run(complete, broken)
    reader = TrajectoryChunkReader(broken)
    if len(reader) < 2:
        pytest.skip('tiny solve produced a single chunk')
    first = reader.load_chunk(0)
    second = reader.load_chunk(1)
    rho = np.array(second.rho_nodes, copy=True)
    rho[0] = np.nextafter(rho[0], np.inf)
    raw = deterministic_npz_bytes({
        'rho_nodes': rho,
        'state_nodes': second.state_nodes,
        'dense_corrections': second.dense_corrections,
    })
    (broken / chunk_filename(1)).write_bytes(raw)
    manifest = json.loads((broken / 'manifest.json').read_text())
    manifest['chunks'][1]['sha256'] = sha256(raw).hexdigest()
    manifest['chunks'][1]['bytes'] = len(raw)
    (broken / 'manifest.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    with pytest.raises(ValueError, match='bitwise continuous|rho bounds'):
        list(TrajectoryChunkReader(broken))
    # The first chunk still loads by itself after the second is corrupted.
    restored = TrajectoryChunkReader(broken).load_chunk(0)
    np.testing.assert_array_equal(restored.rho_nodes, first.rho_nodes)
    np.testing.assert_array_equal(restored.state_nodes, first.state_nodes)
    np.testing.assert_array_equal(restored.dense_corrections, first.dense_corrections)


def test_manifest_records_the_empty_incomplete_start_and_rejects_false_bounds(tmp_path):
    directory=tmp_path/'stream'
    writer=TrajectoryChunkWriter(directory)
    initial=json.loads((directory/'manifest.json').read_text())
    assert initial['complete'] is False and initial['accepted_steps']==0
    args,options=tiny()
    stream_ks_trajectory(*args,on_segment=writer,**options)
    writer.close(complete=True)
    path=directory/'manifest.json'
    original=json.loads(path.read_text())
    changed=json.loads(path.read_text());changed['rho_end']=.9
    path.write_text(json.dumps(changed))
    with pytest.raises(ValueError,match='global rho bounds'):
        TrajectoryChunkReader(directory)
    changed=original
    changed['chunks'][0]['rho_end']=np.nextafter(changed['chunks'][0]['rho_end'],np.inf).item()
    path.write_text(json.dumps(changed))
    with pytest.raises(ValueError,match='chunk rho bounds'):
        TrajectoryChunkReader(directory).load_chunk(0)


def _copy_run(source, destination):
    destination.mkdir()
    for path in Path(source).iterdir():
        (destination / path.name).write_bytes(path.read_bytes())
