"""Stream accepted KS reconstruction segments without retaining the trajectory.

The owned difference evolution runs from its exact code object in a private
namespace whose DOP853 factory observes dense output.  No process global is
patched.  Chunk files hold at most four complete segments; source/history
bindings and the final prepared state remain caller responsibilities.
"""
from dataclasses import dataclass
from hashlib import sha256
from io import BytesIO
import json
import os
from pathlib import Path
import tempfile
from types import FunctionType

import numpy as np
from scipy.integrate import DOP853

from .nsc_ks_difference_envelope import evolve_ks_difference_envelope
from .nsc_ks_trajectory import TrajectorySegment
from .nsc_mode_resolved_cauchy_state import deterministic_npz_bytes


SCHEMA = 'NSC-KS-STREAMING-TRAJECTORY-v1'
MANIFEST_NAME = 'manifest.json'
MAX_SEGMENTS_PER_CHUNK = 4
CHUNK_ARRAYS = ('dense_corrections', 'rho_nodes', 'state_nodes')
# Documented later control; this module does not execute it.
INTENDED_LATER_CONTROL = {
    'computational_nodes': 2048,
    'computational_length': '205/512',
    'physical_interval': 'S(1)+[.12,.18]',
    'axial_support': 'usual compact support unchanged',
    'rho_up': 'existing source control',
    'rtol': 5e-14,
    'atol': 5e-19,
    'max_step': '1/8192',
    'tangents': 'zero',
}


@dataclass(frozen=True)
class StreamSummary:
    accepted_steps: int
    rho_start: float
    rho_end: float


@dataclass(frozen=True)
class TrajectoryChunk:
    index: int
    rho_nodes: object
    state_nodes: object
    dense_corrections: object

    def __post_init__(self):
        rho = np.ascontiguousarray(self.rho_nodes, dtype=np.float64)
        state = np.ascontiguousarray(self.state_nodes)
        corrections = np.ascontiguousarray(self.dense_corrections)
        _require_chunk_arrays(rho, state, corrections)
        rho.setflags(write=False)
        state.setflags(write=False)
        corrections.setflags(write=False)
        object.__setattr__(self, 'rho_nodes', rho)
        object.__setattr__(self, 'state_nodes', state)
        object.__setattr__(self, 'dense_corrections', corrections)

    def __len__(self):
        return len(self.dense_corrections)

    def segment(self, index):
        if index < 0 or index >= len(self):
            raise IndexError('chunk segment index out of range')
        return TrajectorySegment(
            float(self.rho_nodes[index]), float(self.rho_nodes[index + 1]),
            self.state_nodes[index], self.state_nodes[index + 1],
            self.dense_corrections[index])

    def segments(self):
        return tuple(self.segment(index) for index in range(len(self)))


def stream_ks_trajectory(*args, on_segment, **kwargs):
    """Same upstream/history solve; one callback per accepted step, no segment list."""
    if not callable(on_segment):
        raise TypeError('on_segment callback required')
    count = 0
    first_rho = None
    previous_rho = None
    previous_end = None

    def emit(segment):
        nonlocal count, first_rho, previous_rho, previous_end
        if count == 0:
            first_rho = segment.rho_start
        elif previous_rho != segment.rho_start or not np.array_equal(previous_end, segment.start):
            raise ValueError('captured reconstruction endpoints are not continuous')
        on_segment(segment)
        count += 1
        previous_rho = segment.rho_end
        previous_end = segment.end

    class Recorder(DOP853):
        def step(self):
            message = super().step()
            if self.status != 'failed':
                dense = self.dense_output()
                if dense.F.shape != (7, self.n):
                    raise ValueError('unsupported SciPy DOP853 reconstruction layout')
                emit(TrajectorySegment(self.t_old, self.t, self.y_old, self.y, dense.F[1:]))
            return message

    original = evolve_ks_difference_envelope
    namespace = {**original.__globals__, 'DOP853': Recorder}
    runner = FunctionType(original.__code__, namespace, original.__name__,
                          original.__defaults__, original.__closure__)
    runner.__kwdefaults__ = dict(original.__kwdefaults__ or {})
    prepared = runner(*args, **kwargs)
    if not count or first_rho != prepared.rho_up or previous_rho != 1.:
        raise ValueError('captured trajectory does not cover the prepared slab')
    return prepared, StreamSummary(count, first_rho, previous_rho)


def chunk_filename(index):
    if isinstance(index, bool) or not isinstance(index, int) or index < 0:
        raise ValueError('nonnegative chunk index required')
    return f'chunk-{index:06d}.npz'


class TrajectoryChunkWriter:
    """Deterministic NPZ chunks of at most four segments; atomic incomplete manifest."""

    def __init__(self, directory):
        directory = Path(directory)
        if directory.exists():
            if any(directory.iterdir()):
                raise FileExistsError(f'streaming trajectory run already exists: {directory}')
        else:
            directory.mkdir(parents=True)
        self.directory = directory
        self._pending = []
        self._closed = False
        self._last_rho = None
        self._last_end = None
        self._manifest = {
            'schema': SCHEMA,
            'complete': False,
            'max_segments_per_chunk': MAX_SEGMENTS_PER_CHUNK,
            'accepted_steps': 0,
            'chunk_count': 0,
            'state_size': None,
            'rho_start': None,
            'rho_end': None,
            'chunks': [],
        }
        self._write_manifest()

    @property
    def buffered_segments(self):
        return len(self._pending)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, *_):
        self.close(complete=exc_type is None)
        return False

    def __call__(self, segment):
        if self._closed:
            raise ValueError('cannot append to a closed streaming trajectory')
        if not isinstance(segment, TrajectorySegment):
            raise TypeError('TrajectorySegment required')
        if self._pending:
            previous = self._pending[-1]
            if previous.rho_end != segment.rho_start or not np.array_equal(previous.end, segment.start):
                raise ValueError('captured reconstruction endpoints are not continuous')
        elif self._last_end is not None:
            if self._last_rho != segment.rho_start or not np.array_equal(self._last_end, segment.start):
                raise ValueError('captured reconstruction endpoints are not continuous')
        self._pending.append(segment)
        if len(self._pending) > MAX_SEGMENTS_PER_CHUNK:
            raise ValueError('chunk buffer exceeded the four-segment bound')
        if len(self._pending) == MAX_SEGMENTS_PER_CHUNK:
            self._flush()

    def close(self, *, complete=False):
        if self._closed:
            if complete and not self._manifest['complete']:
                raise ValueError('an interrupted stream cannot be labeled complete')
            return
        self._flush()
        if complete:
            if self._manifest['chunk_count'] < 1:
                raise ValueError('complete stream requires at least one written chunk')
            self._manifest['complete'] = True
        self._closed = True
        self._write_manifest()

    def _flush(self):
        if not self._pending:
            return
        arrays = _chunk_arrays(self._pending)
        raw = deterministic_npz_bytes(arrays)
        index = self._manifest['chunk_count']
        name = chunk_filename(index)
        target = self.directory / name
        if target.exists():
            raise FileExistsError(f'chunk already exists: {target}')
        _replace_bytes(target, raw)
        rho = arrays['rho_nodes']
        state = arrays['state_nodes']
        entry = {
            'index': index,
            'filename': name,
            'sha256': sha256(raw).hexdigest(),
            'bytes': len(raw),
            'segment_count': len(self._pending),
            'rho_start': float(rho[0]),
            'rho_end': float(rho[-1]),
        }
        if self._manifest['state_size'] is None:
            self._manifest['state_size'] = int(state.shape[1])
            self._manifest['rho_start'] = entry['rho_start']
        elif self._manifest['state_size'] != int(state.shape[1]):
            raise ValueError('chunk state size changed during the stream')
        self._manifest['chunks'].append(entry)
        self._manifest['chunk_count'] = index + 1
        self._manifest['accepted_steps'] += len(self._pending)
        self._manifest['rho_end'] = entry['rho_end']
        self._last_rho = self._pending[-1].rho_end
        self._last_end = self._pending[-1].end
        self._pending.clear()
        self._write_manifest()

    def _write_manifest(self):
        _replace_bytes(
            self.directory / MANIFEST_NAME,
            json.dumps(self._manifest, indent=2, sort_keys=True, allow_nan=False).encode() + b'\n')


class TrajectoryChunkReader:
    """Replay one verified chunk, or iterate chunks, without loading the full run."""

    def __init__(self, directory, *, require_complete=True):
        directory = Path(directory)
        manifest = json.loads((directory / MANIFEST_NAME).read_text())
        if manifest.get('schema') != SCHEMA:
            raise ValueError('unsupported streaming trajectory schema')
        if require_complete and not manifest.get('complete'):
            raise ValueError('incomplete streaming trajectory run')
        chunks = list(manifest.get('chunks') or [])
        if manifest.get('complete') and not chunks:
            raise ValueError('complete stream requires stored chunks')
        if manifest.get('max_segments_per_chunk') != MAX_SEGMENTS_PER_CHUNK:
            raise ValueError('manifest chunk bound is not four segments')
        if manifest.get('chunk_count') != len(chunks):
            raise ValueError('manifest chunk_count does not match the chunk list')
        accepted = 0
        for index, entry in enumerate(chunks):
            if entry.get('index') != index or entry.get('filename') != chunk_filename(index):
                raise ValueError('chunk files are not monotonically numbered')
            count = int(entry['segment_count'])
            if count < 1 or count > MAX_SEGMENTS_PER_CHUNK:
                raise ValueError('declared chunk segment count is outside 1..4')
            accepted += count
        if manifest.get('accepted_steps') != accepted:
            raise ValueError('manifest accepted_steps does not match declared chunks')
        if chunks and (manifest.get('rho_start') != chunks[0]['rho_start']
                       or manifest.get('rho_end') != chunks[-1]['rho_end']):
            raise ValueError('manifest global rho bounds do not match its chunks')
        self.directory = directory
        self.manifest = manifest

    def __len__(self):
        return int(self.manifest['chunk_count'])

    def load_chunk(self, index):
        if index < 0 or index >= len(self):
            raise IndexError('chunk index out of range')
        entry = self.manifest['chunks'][index]
        path = self.directory / entry['filename']
        raw = path.read_bytes()
        if len(raw) != int(entry['bytes']) or sha256(raw).hexdigest() != entry['sha256']:
            raise ValueError('chunk sha256 does not match the manifest')
        with np.load(BytesIO(raw), allow_pickle=False) as payload:
            if set(payload.files) != set(CHUNK_ARRAYS):
                raise ValueError('chunk arrays must be rho_nodes, state_nodes, dense_corrections')
            arrays = {name: np.array(payload[name], copy=True) for name in CHUNK_ARRAYS}
        chunk = TrajectoryChunk(index, arrays['rho_nodes'], arrays['state_nodes'],
                                arrays['dense_corrections'])
        if len(chunk) != int(entry['segment_count']):
            raise ValueError('chunk segment count does not match the manifest')
        if self.manifest['state_size'] != int(chunk.state_nodes.shape[1]):
            raise ValueError('chunk state size does not match the manifest')
        if float(chunk.rho_nodes[0]) != entry['rho_start'] or float(chunk.rho_nodes[-1]) != entry['rho_end']:
            raise ValueError('chunk rho bounds do not match the manifest')
        return chunk

    def __iter__(self):
        previous_rho = None
        previous_state = None
        for index in range(len(self)):
            chunk = self.load_chunk(index)
            if previous_rho is not None:
                if not np.array_equal(previous_rho, chunk.rho_nodes[:1]):
                    raise ValueError('chunk rho endpoints are not bitwise continuous')
                if not np.array_equal(previous_state, chunk.state_nodes[0]):
                    raise ValueError('chunk state endpoints are not bitwise continuous')
            previous_rho = np.array(chunk.rho_nodes[-1:], copy=True)
            previous_state = np.array(chunk.state_nodes[-1], copy=True)
            yield chunk


def _chunk_arrays(segments):
    if not segments or len(segments) > MAX_SEGMENTS_PER_CHUNK:
        raise ValueError('chunk must contain 1 to 4 complete segments')
    for previous, following in zip(segments, segments[1:]):
        if previous.rho_end != following.rho_start or not np.array_equal(previous.end, following.start):
            raise ValueError('captured reconstruction endpoints are not continuous')
    rho = np.ascontiguousarray(
        np.array([segments[0].rho_start, *[part.rho_end for part in segments]], dtype=np.float64))
    state = np.ascontiguousarray(np.stack([segments[0].start, *[part.end for part in segments]]))
    corrections = np.ascontiguousarray(np.stack([part.coefficients for part in segments]))
    _require_chunk_arrays(rho, state, corrections)
    return {'rho_nodes': rho, 'state_nodes': state, 'dense_corrections': corrections}


def _require_chunk_arrays(rho, state, corrections):
    if rho.ndim != 1 or rho.dtype != np.float64 or len(rho) < 2:
        raise ValueError('rho_nodes must be a float64 node vector')
    if not np.isfinite(rho).all() or np.any(np.diff(rho) >= 0):
        raise ValueError('rho_nodes must be strictly decreasing and finite')
    segments = len(rho) - 1
    if segments > MAX_SEGMENTS_PER_CHUNK:
        raise ValueError('chunk contains more than four segments')
    if state.ndim != 2 or state.shape[0] != len(rho) or not np.issubdtype(state.dtype, np.complexfloating):
        raise ValueError('state_nodes must match rho nodes as complex endpoint vectors')
    if corrections.shape != (segments, 6, state.shape[1]) or not np.issubdtype(corrections.dtype, np.complexfloating):
        raise ValueError('dense_corrections must have shape (segments, 6, state)')
    if not np.isfinite(state).all() or not np.isfinite(corrections).all():
        raise ValueError('finite chunk state and dense corrections required')


def _replace_bytes(path, data):
    path = Path(path)
    handle, temporary = tempfile.mkstemp(dir=path.parent, prefix=path.name + '.', suffix='.tmp')
    try:
        with os.fdopen(handle, 'wb') as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
