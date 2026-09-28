"""Value-only evaluation identity: structured binding plus a deterministic digest.

Production caches bind source, history, channel, solver settings, nodes and
implementation hashes. A missing or mismatched field is rejected before a
run is marked complete and before pending family jobs are launched. v1
records stay immutable; they are replayed only through an explicit legacy
path that authenticates the recorded hashes, including Git-pinned source
bytes when the working tree has moved on.

Public v2 API for calibration and tests
---------------------------------------
Constants:
    SCHEMA_V1, SCHEMA_V2, BINDING_SCHEMA
    PRODUCTION_MODE, LEGACY_REPLAY_MODE
    BINDING_DIMENSIONS, IMPLEMENTATION_OWNERS, GIT_SOURCE_PREFIXES
    BindingMismatch

Canonicalisation:
    hex_float, array_digest, file_digest, canonicalize, canonical_dumps,
    structured_digest, family_key_text, family_key_tuple

Binding construction:
    solver_binding, coefficient_binding, baseline_binding, target_binding,
    grid_binding, geometry_binding, batch_source_binding, family_source_binding,
    source_identity_of, implementation_hashes, split_hashes, evaluation_binding,
    require_complete_binding, require_family_source, binding_digest,
    family_identity, dimension_payload, compare_bindings

Authentication and records:
    classify_family_record, classify_top_record
    resolve_recorded_source_bytes, authenticate_source_hashes
    authenticate_payload, compare_family_cache, legacy_family_consistent
    payload_descriptor, write_bytes_atomic, write_json_atomic
"""
__all__ = (
    'SCHEMA_V1', 'SCHEMA_V2', 'BINDING_SCHEMA', 'PRODUCTION_MODE',
    'LEGACY_REPLAY_MODE', 'BINDING_DIMENSIONS', 'IMPLEMENTATION_OWNERS',
    'GIT_SOURCE_PREFIXES', 'BindingMismatch', 'hex_float', 'array_digest',
    'file_digest', 'canonicalize', 'canonical_dumps', 'structured_digest',
    'family_key_text', 'family_key_tuple', 'solver_binding',
    'coefficient_binding', 'baseline_binding', 'target_binding', 'grid_binding',
    'geometry_binding', 'batch_source_binding', 'family_source_binding',
    'source_identity_of', 'implementation_hashes', 'split_hashes',
    'evaluation_binding', 'require_complete_binding', 'require_family_source',
    'binding_digest', 'family_identity', 'dimension_payload', 'compare_bindings',
    'classify_family_record', 'classify_top_record',
    'resolve_recorded_source_bytes', 'authenticate_source_hashes',
    'authenticate_payload', 'compare_family_cache', 'legacy_family_consistent',
    'payload_descriptor', 'write_bytes_atomic', 'write_json_atomic',
)
from collections.abc import Mapping, Sequence
import ast
from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath
import subprocess

import numpy as np

from .nsc_evolved_incoming_state import _digest_arrays
from .nsc_ks_profile_identity import profile_description, profile_identity
from .nsc_ks_source_envelope import (
    RHO_SIGMA, _sample_axial_profiles, physical_incoming_interval,
    usual_axial_support)


SCHEMA_V1 = 'NSC-KS-GATE-VALUE-v1'
SCHEMA_V2 = 'NSC-KS-GATE-VALUE-v2'
BINDING_SCHEMA = 'NSC-KS-EVALUATION-BINDING-v2'
PRODUCTION_MODE = 'production'
LEGACY_REPLAY_MODE = 'legacy-replay'
BINDING_DIMENSIONS = (
    'source', 'history', 'channel', 'settings', 'nodes', 'code',
)
IMPLEMENTATION_OWNERS = (
    'scripts/derive_nsc_ks_gate_value.py',
    'docs/nsc-ks-gate-value.md',
    'src/recursive_horizons/nsc_ks_evaluation_binding.py',
    'src/recursive_horizons/nsc_ks_value_evaluator.py',
    'src/recursive_horizons/nsc_ks_family_pool.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_ks_energy_propagator.py',
    'src/recursive_horizons/nsc_ks_primal_step_control.py',
)
GIT_SOURCE_PREFIXES = ('scripts/', 'src/', 'docs/')
_HASH = 64
_REQUIRED_BINDING = (
    'schema', 'profile_identity', 'history_path', 'tangents', 'node_count',
    'target', 'grid', 'geometry', 'solver', 'coefficients', 'baseline',
    'source_families', 'source_identity', 'implementation_hashes',
    'input_hashes', 'source_hashes',
)
_REQUIRED_BATCH = (
    'batch_id', 'panel_name', 'original_panel', 'rows', 'source_digest',
    'preparation_digest', 'covariance_digest', 'covariance_shape',
    'weights_digest', 'energies_digest', 'initial_columns_digest',
    'initial_columns_shape', 'channel',
)
_REQUIRED_CHANNEL = (
    'group', 'angular_sign', 'energy_sign', 'mass', 'angular', 'rho_up',
    'copy_count', 'degeneracy', 'multiplicity', 'compact_mass',
    'angular_eigenvalue',
)
_REQUIRED_SOLVER = (
    'integrator', 'rtol', 'atol', 'max_step', 'degree', 'phase_nodes',
    'grid_nodes', 'length', 'step_control',
)
_REQUIRED_TARGET = ('digest', 'shape', 'values')
_REQUIRED_GRID = ('digest', 'shape', 'nodes', 'length')
_REQUIRED_GEOMETRY = (
    'profile_identity', 'normal_windows', 'axial_support', 'domain',
    'w_digest', 'U_digest', 'w_shape', 'U_shape',
)
_REQUIRED_BASELINE = ('digest', 'shape', 'values')


class BindingMismatch(ValueError):
    """One named binding dimension does not match the declared evaluation."""

    def __init__(self, dimension, message):
        self.dimension = dimension
        super().__init__(f'{dimension} binding mismatch: {message}')


def family_key_text(key):
    group, sign = int(key[0]), int(key[1])
    if group < 0 or group > 99 or sign not in (-1, 1):
        raise ValueError('family key is not a retained angular family')
    return f'{group:02d},{sign:+d}'


def family_key_tuple(value):
    if isinstance(value, str):
        group_text, sign_text = value.split(',')
        return (int(group_text), int(sign_text))
    return (int(value[0]), int(value[1]))


def hex_float(value):
    number = float(value)
    if not np.isfinite(number):
        raise ValueError('finite binding number required')
    return number.hex()


def array_digest(value):
    return _digest_arrays(np.ascontiguousarray(value))


def file_digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def canonicalize(value):
    if isinstance(value, Mapping):
        return {str(key): canonicalize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [canonicalize(item) for item in value]
    if isinstance(value, np.ndarray):
        raise TypeError('arrays must be digested or listed before binding')
    if isinstance(value, (np.floating, float)):
        return hex_float(value)
    if isinstance(value, (np.integer, int)) and not isinstance(value, (bool, np.bool_)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if value is None or isinstance(value, str):
        return value
    raise TypeError(f'unsupported binding type {type(value).__name__}')


def canonical_dumps(value):
    return json.dumps(
        canonicalize(value), sort_keys=True, separators=(',', ':'),
        allow_nan=False)


def structured_digest(value):
    return sha256(canonical_dumps(value).encode()).hexdigest()


def repository_relative(path, root):
    root = Path(root).resolve()
    candidate = Path(path)
    if candidate.is_absolute():
        resolved = candidate.resolve()
        if not resolved.is_relative_to(root):
            raise ValueError('binding path escapes repository')
        candidate = resolved.relative_to(root)
    text = candidate.as_posix()
    parsed = PurePosixPath(text)
    if parsed.is_absolute() or '..' in parsed.parts or '\\' in text or not text:
        raise ValueError('binding path must be canonical and repository-relative')
    return text


def _require_hash(value, name):
    if not isinstance(value, str) or len(value) != _HASH or any(
            char not in '0123456789abcdef' for char in value):
        raise ValueError(name + ' must be a lowercase SHA-256 hex digest')
    return value


def _require_keys(record, required, label):
    if not isinstance(record, Mapping):
        raise ValueError(label + ' must be an object')
    missing = [name for name in required if name not in record or record[name] is None]
    if missing:
        raise ValueError('missing binding field: ' + label + '.' + missing[0])
    return record


def _field(item, name):
    if isinstance(item, Mapping):
        return item[name]
    return getattr(item, name)


def _optional(item, name, default=None):
    if isinstance(item, Mapping):
        return item.get(name, default)
    return getattr(item, name, default)


def _shape(value):
    array = np.asarray(value)
    return [int(dim) for dim in array.shape]


def listed_hex(values):
    array = np.asarray(values, float)
    if array.ndim != 1:
        raise ValueError('one-dimensional hex list required')
    return [hex_float(item) for item in array]


def solver_binding(solver):
    record = _require_keys(dict(solver), _REQUIRED_SOLVER, 'solver')
    return {
        'integrator': str(record['integrator']),
        'rtol': hex_float(record['rtol']),
        'atol': hex_float(record['atol']),
        'max_step': hex_float(record['max_step']),
        'degree': int(record['degree']),
        'phase_nodes': int(record['phase_nodes']),
        'grid_nodes': int(record['grid_nodes']),
        'length': hex_float(record['length']),
        'step_control': str(record['step_control']),
    }


def coefficient_binding(coeff):
    if not isinstance(coeff, Mapping) or not coeff:
        raise ValueError('missing binding field: coefficients')
    bound = {}
    for name in sorted(coeff):
        value = coeff[name]
        if isinstance(value, np.ndarray):
            bound[name] = {
                'digest': array_digest(value),
                'shape': _shape(value),
                'values': [hex_float(item) for item in np.asarray(value, float).ravel()],
            }
        else:
            bound[name] = hex_float(value)
    return bound


def baseline_binding(baseline):
    array = np.asarray(baseline, float)
    if array.size == 0 or not np.isfinite(array).all():
        raise ValueError('finite baseline values required')
    return {
        'digest': array_digest(array),
        'shape': _shape(array),
        'values': [hex_float(item) for item in array.ravel()],
    }


def target_binding(target, node_count):
    array = np.asarray(target, float)
    if array.ndim != 1 or len(array) != int(node_count):
        raise ValueError('target nodes must match the declared node count')
    if not np.isfinite(array).all():
        raise ValueError('finite target nodes required')
    return {
        'digest': array_digest(array),
        'shape': _shape(array),
        'values': listed_hex(array),
    }


def grid_binding(grid, solver):
    array = np.asarray(grid, float)
    if array.ndim != 1 or len(array) != int(solver['grid_nodes']):
        raise ValueError('computational grid must match solver grid_nodes')
    if not np.isfinite(array).all():
        raise ValueError('finite computational grid required')
    return {
        'digest': array_digest(array),
        'shape': _shape(array),
        'nodes': int(solver['grid_nodes']),
        'length': hex_float(solver['length']),
    }


def geometry_binding(family, grid, target, solver):
    metric = family.metric()
    w, u = _sample_axial_profiles(metric.directions, grid)
    support = usual_axial_support()
    domain = physical_incoming_interval()
    identity = profile_identity(family, include_normal_window=True)
    return {
        'profile_identity': identity,
        'profile_description': profile_description(
            family, include_normal_window=True),
        'amplitudes': [hex_float(value) for value in metric.amplitudes],
        'normal_windows': [
            {'inner': hex_float(direction.inner_radius),
             'outer': hex_float(direction.outer_radius)}
            for direction in metric.directions
        ],
        'axial_support': [hex_float(support[0]), hex_float(support[1])],
        'domain': [hex_float(domain[0]), hex_float(domain[1])],
        'rho_sigma': hex_float(RHO_SIGMA),
        'w_digest': array_digest(w),
        'U_digest': array_digest(u),
        'w_shape': _shape(w),
        'U_shape': _shape(u),
        'target_digest': array_digest(np.asarray(target, float)),
        'grid_digest': array_digest(np.asarray(grid, float)),
        'grid_nodes': int(solver['grid_nodes']),
        'length': hex_float(solver['length']),
    }


def _source_arrays(batch):
    source = _optional(batch, 'source')
    if source is None:
        return None
    return (
        np.asarray(_field(source, 'covariance')),
        np.asarray(_field(source, 'column_weights')),
        np.asarray(_field(source, 'energies')),
    )


def _source_digest(batch, arrays):
    source = _optional(batch, 'source')
    recorded = _optional(source, 'digest') if source is not None else None
    if recorded:
        return str(recorded)
    if arrays is not None:
        from .nsc_evolved_incoming_state import FixedSourcePreparation
        return FixedSourcePreparation(arrays[0], arrays[1], arrays[2]).digest
    return str(_field(batch, 'source_digest'))


def _int_field(channel, batch, name):
    value = _optional(channel, name, _optional(batch, name))
    if value is None:
        raise ValueError('missing binding field: channel.' + name)
    return int(value)


def batch_source_binding(batch, channel=None):
    channel = {} if channel is None else dict(channel)
    arrays = _source_arrays(batch)
    if arrays is not None:
        covariance, weights, energies = arrays
        covariance_digest = array_digest(covariance)
        weights_digest = array_digest(weights)
        energies_digest = array_digest(energies)
        covariance_shape = _shape(covariance)
    else:
        covariance_digest = _field(batch, 'covariance_digest')
        weights_digest = _field(batch, 'weights_digest')
        energies_digest = _field(batch, 'energies_digest')
        covariance_shape = [int(v) for v in _field(batch, 'covariance_shape')]
    initial = _optional(batch, 'initial_columns')
    if initial is not None:
        initial = np.asarray(initial)
        initial_digest = array_digest(initial)
        initial_shape = _shape(initial)
    else:
        initial_digest = _field(batch, 'initial_columns_digest')
        initial_shape = [int(v) for v in _field(batch, 'initial_columns_shape')]
    group = _int_field(channel, batch, 'group')
    angular_sign = int(_optional(channel, 'angular_sign', _field(batch, 'angular_sign')))
    energy_sign = int(_optional(channel, 'energy_sign', _field(batch, 'energy_sign')))
    mass = _optional(channel, 'compact_mass', _field(batch, 'mass'))
    angular = _optional(channel, 'angular_eigenvalue')
    if angular is None:
        angular = abs(float(_field(batch, 'angular')))
    signed_angular = _optional(batch, 'angular')
    if signed_angular is None:
        signed_angular = angular_sign * float(angular)
    copy_count = _int_field(channel, batch, 'copy_count')
    degeneracy = _int_field(channel, batch, 'degeneracy')
    n_signs = 1 if float(angular) == 0.0 else 2
    multiplicity = _optional(channel, 'multiplicity', copy_count * degeneracy / n_signs)
    panel_name = str(_field(batch, 'panel_name'))
    rows = [int(v) for v in _field(batch, 'rows')]
    preparation = _optional(batch, 'preparation_digest')
    if preparation is None:
        raise ValueError('missing binding field: batch.preparation_digest')
    return {
        'batch_id': str(_optional(
            batch, 'batch_id',
            f'{panel_name}:{energy_sign:+d}:{rows[0]}:{rows[1]}')),
        'panel_name': panel_name,
        'original_panel': str(_optional(batch, 'original_panel', panel_name)),
        'rows': rows,
        'source_digest': _require_hash(_source_digest(batch, arrays), 'source_digest'),
        'preparation_digest': _require_hash(str(preparation), 'preparation_digest'),
        'covariance_digest': _require_hash(str(covariance_digest), 'covariance_digest'),
        'covariance_shape': covariance_shape,
        'weights_digest': _require_hash(str(weights_digest), 'weights_digest'),
        'energies_digest': _require_hash(str(energies_digest), 'energies_digest'),
        'initial_columns_digest': _require_hash(str(initial_digest), 'initial_columns_digest'),
        'initial_columns_shape': initial_shape,
        'channel': {
            'group': group,
            'angular_sign': angular_sign,
            'energy_sign': energy_sign,
            'mass': hex_float(mass),
            'angular': hex_float(signed_angular),
            'rho_up': hex_float(_field(batch, 'rho_up')),
            'copy_count': copy_count,
            'degeneracy': degeneracy,
            'multiplicity': hex_float(multiplicity),
            'compact_mass': hex_float(mass),
            'angular_eigenvalue': hex_float(angular),
        },
    }


def _source_only(batch):
    return {name: batch[name] for name in (
        'batch_id', 'panel_name', 'original_panel', 'rows', 'source_digest',
        'preparation_digest', 'covariance_digest', 'covariance_shape',
        'weights_digest', 'energies_digest', 'initial_columns_digest',
        'initial_columns_shape')}


def source_identity_of(batches):
    return structured_digest([_source_only(batch) for batch in batches])


def family_source_binding(entries):
    if not isinstance(entries, Sequence) or isinstance(entries, (str, bytes)) or not entries:
        raise ValueError('missing binding field: source family batches')
    batches = []
    for item in entries:
        if isinstance(item, Sequence) and not isinstance(item, (str, bytes, Mapping)):
            if len(item) != 2:
                raise ValueError('each source entry is (batch, channel)')
            batches.append(batch_source_binding(item[0], item[1]))
        else:
            batches.append(batch_source_binding(item))
    batches = sorted(batches, key=lambda row: row['batch_id'])
    ids = [row['batch_id'] for row in batches]
    if len(ids) != len(set(ids)):
        raise ValueError('duplicate source batch ids')
    return {
        'batches': batches,
        'source_identity': source_identity_of(batches),
        'batch_count': len(batches),
    }


def implementation_hashes(root, owners=IMPLEMENTATION_OWNERS):
    root = Path(root)
    pending, visited = list(owners), set()
    hashes = {}

    def add_module(parts):
        if not parts:
            return
        if parts[0] == 'recursive_horizons':
            base = root / 'src'
        elif len(parts) == 1 and parts[0].startswith(('derive_', 'check_', 'validate_')):
            base = root / 'scripts'
        else:
            return
        candidate = base.joinpath(*parts).with_suffix('.py')
        if not candidate.is_file():
            candidate = base.joinpath(*parts, '__init__.py')
        if candidate.is_file():
            pending.append(candidate.relative_to(root).as_posix())

    while pending:
        name = pending.pop()
        if name in visited:
            continue
        visited.add(name)
        path = root / name
        if not path.is_file():
            raise ValueError('missing binding field: implementation hash ' + name)
        if path.suffix != '.py':
            hashes[name] = file_digest(path)
            continue
        tree = ast.parse(path.read_text(), filename=name)
        if path.name == '__init__.py' and all(
                (isinstance(item, ast.Expr) and isinstance(item.value, ast.Constant)
                 and isinstance(item.value.value, str))
                or (isinstance(item, ast.Assign)
                    and all(isinstance(target, ast.Name) and target.id == '__version__'
                            for target in item.targets)
                    and isinstance(item.value, ast.Constant) and isinstance(item.value.value, str))
                for item in tree.body):
            # Publication changes the package version, not these equations.
            # Any executable initializer or other assignment is still bound.
            continue
        hashes[name] = file_digest(path)
        package = tuple(path.parent.relative_to(root / 'src').parts) if name.startswith('src/') else ()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    add_module(tuple(alias.name.split('.')))
            elif isinstance(node, ast.ImportFrom):
                if node.level:
                    if node.level > len(package):
                        continue
                    prefix = package[:len(package)-node.level+1]
                else:
                    prefix = ()
                parts = prefix + tuple(node.module.split('.') if node.module else ())
                add_module(parts)
                for alias in node.names:
                    if alias.name != '*':
                        add_module(parts + tuple(alias.name.split('.')))
    return hashes


def split_hashes(hashes, owners=IMPLEMENTATION_OWNERS):
    if not isinstance(hashes, Mapping) or not hashes:
        raise ValueError('missing binding field: source_hashes')
    missing = [name for name in owners if name not in hashes]
    if missing:
        raise ValueError('missing binding field: implementation hash ' + missing[0])
    impl = {name: _require_hash(str(value), name)
            for name, value in hashes.items()
            if name in owners or name.startswith(('src/', 'scripts/'))}
    inputs = {
        name: _require_hash(str(value), name)
        for name, value in hashes.items() if name not in impl
    }
    if not inputs:
        raise ValueError('missing binding field: input_hashes')
    return impl, inputs


def evaluation_binding(
        *, family=None, target, grid, solver, hashes, history_path,
        profile_identity_value, coefficients, baseline, source_families,
        node_count, tangents='zero', geometry=None):
    if tangents != 'zero':
        raise ValueError('value evaluation binds tangents=zero')
    solver_bound = solver_binding(solver)
    if geometry is None:
        if family is None:
            raise ValueError('missing binding field: geometry')
        geometry = geometry_binding(family, grid, target, solver)
    if geometry.get('profile_identity') != profile_identity_value:
        raise BindingMismatch(
            'history', 'sampled geometry identity differs from the declared profile')
    impl, inputs = split_hashes(hashes)
    if not isinstance(source_families, Mapping) or not source_families:
        raise ValueError('missing binding field: source_families')
    families = {}
    for key, record in source_families.items():
        families[str(key)] = require_family_source(record)
    source_identity = structured_digest({
        key: families[key]['source_identity'] for key in sorted(families)
    })
    binding = {
        'schema': BINDING_SCHEMA,
        'profile_identity': str(profile_identity_value),
        'history_path': str(history_path),
        'tangents': 'zero',
        'node_count': int(node_count),
        'target': target_binding(target, node_count),
        'grid': grid_binding(grid, solver),
        'geometry': geometry,
        'solver': solver_bound,
        'coefficients': coefficient_binding(coefficients),
        'baseline': baseline_binding(baseline),
        'source_families': families,
        'source_identity': source_identity,
        'implementation_hashes': impl,
        'input_hashes': inputs,
        'source_hashes': {**inputs, **impl},
    }
    return require_complete_binding(binding)


def require_family_source(record):
    record = _require_keys(record, ('batches', 'source_identity', 'batch_count'),
                           'source_families')
    batches = record['batches']
    if not isinstance(batches, Sequence) or not batches:
        raise ValueError('missing binding field: source_families.batches')
    cleaned = []
    for batch in batches:
        batch = _require_keys(batch, _REQUIRED_BATCH, 'batch')
        _require_keys(batch['channel'], _REQUIRED_CHANNEL, 'channel')
        cleaned.append(batch)
    expected = source_identity_of(cleaned)
    if expected != record['source_identity']:
        raise BindingMismatch('source', 'family source_identity does not match batches')
    if int(record['batch_count']) != len(cleaned):
        raise BindingMismatch('source', 'batch_count does not match batches')
    return {
        'batches': cleaned,
        'source_identity': record['source_identity'],
        'batch_count': int(record['batch_count']),
    }


def require_complete_binding(binding):
    binding = _require_keys(binding, _REQUIRED_BINDING, 'binding')
    if binding['schema'] != BINDING_SCHEMA:
        raise ValueError('evaluation binding must use ' + BINDING_SCHEMA)
    _require_keys(binding['target'], _REQUIRED_TARGET, 'target')
    _require_keys(binding['grid'], _REQUIRED_GRID, 'grid')
    _require_keys(binding['geometry'], _REQUIRED_GEOMETRY, 'geometry')
    _require_keys(binding['solver'], _REQUIRED_SOLVER, 'solver')
    _require_keys(binding['baseline'], _REQUIRED_BASELINE, 'baseline')
    if not binding['coefficients']:
        raise ValueError('missing binding field: coefficients')
    if not binding['implementation_hashes']:
        raise ValueError('missing binding field: implementation_hashes')
    if not binding['input_hashes']:
        raise ValueError('missing binding field: input_hashes')
    if not binding['source_families']:
        raise ValueError('missing binding field: source_families')
    for key, record in binding['source_families'].items():
        require_family_source(record)
        family_key_tuple(key)
    missing_owners = [name for name in IMPLEMENTATION_OWNERS
                      if name not in binding['implementation_hashes']]
    if missing_owners:
        raise ValueError('missing binding field: implementation hash ' + missing_owners[0])
    return binding


def binding_digest(binding):
    return structured_digest(require_complete_binding(binding))


def family_identity(binding, key):
    text = family_key_text(key) if not isinstance(key, str) else key
    families = binding['source_families']
    if text not in families:
        raise ValueError('missing binding field: source_families.' + text)
    return structured_digest({
        'evaluation_identity': binding_digest(binding),
        'positive_family': list(family_key_tuple(text)),
        'source_identity': families[text]['source_identity'],
    })


def _channel_view(source_families):
    rows = []
    for key in sorted(source_families):
        for batch in source_families[key]['batches']:
            rows.append((key, batch['batch_id'], batch['channel']))
    return rows


def dimension_payload(binding, dimension):
    require_complete_binding(binding)
    if dimension == 'source':
        return {
            'source_identity': binding['source_identity'],
            'families': {
                key: record['source_identity']
                for key, record in binding['source_families'].items()
            },
        }
    if dimension == 'history':
        geometry = binding['geometry']
        return {
            'profile_identity': binding['profile_identity'],
            'history_path': binding['history_path'],
            'w_digest': geometry['w_digest'],
            'U_digest': geometry['U_digest'],
            'normal_windows': geometry['normal_windows'],
            'axial_support': geometry['axial_support'],
            'domain': geometry['domain'],
        }
    if dimension == 'channel':
        return _channel_view(binding['source_families'])
    if dimension == 'settings':
        return binding['solver']
    if dimension == 'nodes':
        return {
            'node_count': binding['node_count'],
            'target': binding['target'],
            'grid': binding['grid'],
        }
    if dimension == 'code':
        return binding['implementation_hashes']
    raise ValueError('unknown binding dimension ' + dimension)


def compare_bindings(expected, actual):
    expected = require_complete_binding(expected)
    actual = require_complete_binding(actual)
    for dimension in BINDING_DIMENSIONS:
        if canonicalize(dimension_payload(expected, dimension)) != canonicalize(
                dimension_payload(actual, dimension)):
            raise BindingMismatch(dimension, 'declared evaluation differs')
    if binding_digest(expected) != binding_digest(actual):
        raise BindingMismatch('source', 'evaluation digest differs')
    return True


def classify_family_record(record):
    if not isinstance(record, Mapping):
        raise ValueError('family record must be an object')
    schema = record.get('schema', SCHEMA_V1)
    if schema == SCHEMA_V2:
        return SCHEMA_V2
    if schema == SCHEMA_V1:
        return SCHEMA_V1
    raise ValueError('unknown value family schema ' + str(schema))


def classify_top_record(record):
    if not isinstance(record, Mapping):
        raise ValueError('value register must be an object')
    schema = record.get('schema', SCHEMA_V1)
    if schema not in (SCHEMA_V1, SCHEMA_V2):
        raise ValueError('unknown value register schema ' + str(schema))
    return schema


def _git_log_commits(root, relative, max_commits, git_log=None):
    if git_log is not None:
        return list(git_log(relative, max_commits))
    output = subprocess.check_output(
        ['git', '-C', str(root), 'log', '-n', str(int(max_commits)),
         '--format=%H', '--', relative],
        text=True, timeout=30)
    return [line.strip() for line in output.splitlines() if line.strip()]


def _git_show(root, commit, relative, git_show=None):
    if git_show is not None:
        return git_show(commit, relative)
    return subprocess.check_output(
        ['git', '-C', str(root), 'show', f'{commit}:{relative}'],
        stderr=subprocess.DEVNULL, timeout=30)


def resolve_recorded_source_bytes(
        root, relative, expected, *, mode, git_commit=None, max_commits=32,
        git_log=None, git_show=None):
    """Return authenticated bytes for one recorded hash.

    Production reads the working tree only. Legacy replay may resolve
    implementation sources from a bounded Git log of that path; numerical
    result artifacts are never walked in Git history.
    """
    root = Path(root).resolve()
    relative = repository_relative(relative, root)
    expected = _require_hash(expected, relative)
    path = root / relative
    if mode == PRODUCTION_MODE:
        if not path.is_file():
            raise ValueError('missing production source: ' + relative)
        raw = path.read_bytes()
        if sha256(raw).hexdigest() != expected:
            raise ValueError('source hash mismatch: ' + relative)
        return raw, 'working-tree'
    if mode != LEGACY_REPLAY_MODE:
        raise ValueError('mode must be production or legacy-replay')
    if path.is_file():
        raw = path.read_bytes()
        if sha256(raw).hexdigest() == expected:
            return raw, 'working-tree'
    if not relative.startswith(GIT_SOURCE_PREFIXES):
        raise ValueError(
            'legacy recorded data hash mismatch: ' + relative
            + ' (working tree is not compared as a silent substitute)')
    if git_commit:
        raw = _git_show(root, git_commit, relative, git_show=git_show)
        if sha256(raw).hexdigest() != expected:
            raise ValueError('pinned git source hash mismatch: ' + relative)
        return raw, git_commit
    for commit in _git_log_commits(root, relative, max_commits, git_log=git_log):
        try:
            raw = _git_show(root, commit, relative, git_show=git_show)
        except (subprocess.CalledProcessError, ValueError, FileNotFoundError):
            continue
        if sha256(raw).hexdigest() == expected:
            return raw, commit
    raise ValueError('historical source bytes not found for ' + relative)


def authenticate_source_hashes(
        hashes, root, *, mode, git_commit=None, git_log=None, git_show=None):
    if not isinstance(hashes, Mapping) or not hashes:
        raise ValueError('recorded source hashes required')
    origins = {}
    for relative, expected in hashes.items():
        _raw, origin = resolve_recorded_source_bytes(
            root, relative, expected, mode=mode, git_commit=git_commit,
            git_log=git_log, git_show=git_show)
        origins[relative] = origin
    return origins


def write_bytes_atomic(path, raw, *, exclusive=False):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f'.{os.getpid()}.tmp')
    created = False
    try:
        with temporary.open('xb') as stream:
            created = True
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        if exclusive:
            os.link(temporary, path)
        else:
            os.replace(temporary, path)
    finally:
        if created and temporary.exists():
            temporary.unlink()
    return path


def write_json_atomic(path, data, *, exclusive=False, indent=None):
    raw = json.dumps(data, sort_keys=True, indent=indent, allow_nan=False) + '\n'
    return write_bytes_atomic(path, raw.encode(), exclusive=exclusive)


def payload_descriptor(path, root):
    path = Path(path)
    relative = repository_relative(path, root)
    raw = path.read_bytes()
    return {'path': relative, 'sha256': sha256(raw).hexdigest(), 'bytes': len(raw)}


def authenticate_payload(record, payload_path, root, *, target=None):
    payload_path = Path(payload_path)
    declared = record.get('payload')
    if not isinstance(declared, Mapping) or 'sha256' not in declared:
        raise ValueError('family record/payload binding is missing')
    raw = payload_path.read_bytes()
    digest = sha256(raw).hexdigest()
    if digest != declared['sha256']:
        raise BindingMismatch('source', 'family payload sha256 changed')
    if declared.get('bytes') is not None and int(declared['bytes']) != len(raw):
        raise BindingMismatch('source', 'family payload size changed')
    declared_path = declared.get('path')
    if declared_path is not None:
        actual = repository_relative(payload_path, root)
        if declared_path != actual:
            raise BindingMismatch('source', 'family payload path changed')
    arrays = {}
    with np.load(payload_path, allow_pickle=False) as handle:
        if 'matter_change' not in handle.files:
            raise ValueError('family payload missing matter_change')
        arrays['matter_change'] = np.array(handle['matter_change'], float, copy=True)
        if 'metadata_json' in handle.files:
            arrays['metadata'] = json.loads(np.asarray(handle['metadata_json']).tobytes())
        if 'z' in handle.files:
            arrays['z'] = np.array(handle['z'], float, copy=True)
    change = arrays['matter_change']
    if change.ndim != 2 or change.shape[1] != 2 or not np.isfinite(change).all():
        raise ValueError('family matter_change must be a finite nodal (N, beta) array')
    if target is not None:
        target = np.asarray(target, float)
        if change.shape[0] != len(target):
            raise BindingMismatch('nodes', 'stored matter_change length differs from target nodes')
        stored_z = arrays.get('z')
        if stored_z is not None and (
                stored_z.shape != target.shape or not np.array_equal(stored_z, target)):
            raise BindingMismatch('nodes', 'stored operator nodes differ from the target grid')
    return arrays


def compare_family_cache(record, expected_binding, key, *, mode, payload_arrays=None):
    schema = classify_family_record(record)
    if mode == PRODUCTION_MODE:
        if schema != SCHEMA_V2:
            raise ValueError(
                'v1 family cache is not a production evaluation; choose a new v2 label')
        expected_binding = require_complete_binding(expected_binding)
        expected_id = family_identity(expected_binding, key)
        if list(record.get('positive_family', [])) != [int(key[0]), int(key[1])]:
            raise BindingMismatch('channel', 'stored family id differs')
        if record.get('family_identity') != expected_id:
            raise BindingMismatch('source', 'stored family identity differs')
        if record.get('evaluation_identity') != binding_digest(expected_binding):
            raise BindingMismatch('source', 'stored evaluation identity differs')
        if record.get('binding_digest') != binding_digest(expected_binding):
            raise BindingMismatch('source', 'record/payload binding digest differs')
        text = family_key_text(key)
        expected_source = expected_binding['source_families'][text]['source_identity']
        if record.get('source_identity') != expected_source:
            raise BindingMismatch('source', 'exact source identity differs')
        if record.get('profile_identity') != expected_binding['profile_identity']:
            raise BindingMismatch('history', 'stored profile identity differs')
        if record.get('tangents') != 'zero':
            raise BindingMismatch('history', 'value family tangents are not zero')
        if int(record.get('node_count', -1)) != int(expected_binding['node_count']):
            raise BindingMismatch('nodes', 'stored node count differs')
        if record.get('solver') is not None and canonicalize(record['solver']) != canonicalize(
                expected_binding['solver']):
            raise BindingMismatch('settings', 'stored solver settings differ')
        if record.get('source_hashes'):
            if record['source_hashes'] != expected_binding['source_hashes']:
                raise BindingMismatch('code', 'stored source hashes differ')
        if payload_arrays is not None:
            metadata = payload_arrays.get('metadata')
            if not isinstance(metadata, Mapping) or payload_arrays.get('z') is None:
                raise BindingMismatch('source', 'v2 payload identity or nodes missing')
            for name in ('family_identity', 'evaluation_identity', 'source_identity', 'binding_digest'):
                if metadata.get(name) != record.get(name):
                    raise BindingMismatch('source', 'payload identity differs: ' + name)
            change = payload_arrays['matter_change']
            expected_shape = [int(expected_binding['node_count']), 2]
            if list(change.shape) != expected_shape:
                raise BindingMismatch('nodes', 'stored matter shape differs')
            if record.get('matter_shape') not in (None, expected_shape):
                raise BindingMismatch('nodes', 'recorded matter shape differs')
        return SCHEMA_V2
    if mode != LEGACY_REPLAY_MODE:
        raise ValueError('mode must be production or legacy-replay')
    if schema != SCHEMA_V1:
        raise ValueError('explicit legacy replay requires a v1 family record')
    if list(record.get('positive_family', [])) != [int(key[0]), int(key[1])]:
        raise ValueError('legacy family id is inconsistent with its filename')
    if record.get('tangents') != 'zero':
        raise ValueError('legacy value family tangents are not zero')
    return SCHEMA_V1


def legacy_family_consistent(top, record):
    if top.get('profile_identity') != record.get('profile_identity'):
        raise ValueError('legacy family profile identity differs from the top record')
    if record.get('tangents') != top.get('tangents', 'zero'):
        raise ValueError('legacy family tangents differ from the top record')
    top_hashes = top.get('source_hashes') or {}
    family_hashes = record.get('source_hashes') or {}
    shared = set(top_hashes) & set(family_hashes)
    if not shared:
        raise ValueError('legacy family has no recorded source hashes in common with the top record')
    for name in sorted(shared):
        if family_hashes[name] != top_hashes[name]:
            raise ValueError('legacy family source hash is inconsistent: ' + name)
    return True
