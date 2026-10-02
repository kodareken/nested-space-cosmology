"""Bounded whole-realization translation on the saved periodic/AP carrier.

Translate the columns of W themselves. Canonical geometry/momentum coordinates
and parent/detail indices stay fixed, so W's spatial interpretation follows the
translated windows. This is a coordinate covariance control, not an ambient-size
or recurrence experiment.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from . import nsc_discovery_episode as episode
from . import nsc_discovery_backend as backend
from . import nsc_discovery_tidal as tidal
from . import nsc_nested_parent_child as model
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin

LAB = Path(__file__).resolve().parents[2]
DEFAULT_EPISODE = LAB / 'results/development/nsc-discovery-episode-v1'
DEFAULT_SHIFT = 5.5
SCHEMA = 'NSC-DISCOVERY-TRANSLATION-v1'
CPU_LIMIT = 30.
PERIOD = 8.


def periodic_shift(values, shift, period=PERIOD):
    """f_shift(x)=f(x-shift), integer Fourier modes along the first axis."""
    values = np.asarray(values)
    count = values.shape[0]
    modes = np.fft.fftfreq(count) * count
    phase = np.exp(-2j * np.pi * modes * float(shift) / period)
    result = np.fft.ifft(np.fft.fft(values, axis=0) * phase.reshape((count,) + (1,) * (values.ndim - 1)), axis=0)
    if not np.iscomplexobj(values):
        if np.max(np.abs(result.imag)) > 1e-10 * max(1., float(np.max(np.abs(result.real)))):
            raise ValueError('translation leaves the real Nyquist carrier; choose a quadrature-lattice offset')
        result = result.real
    return result


def antiperiodic_shift(values, shift, period=PERIOD):
    """Exact half-integer Fourier shift. A full-period translation is minus I."""
    values = np.asarray(values, dtype=complex)
    count = values.shape[0]
    shape = (count,) + (1,) * (values.ndim - 1)
    demodulator = np.exp(-1j * np.pi * np.arange(count) / count).reshape(shape)
    modes = np.fft.fftfreq(count) * count + .5
    phase = np.exp(-2j * np.pi * modes * float(shift) / period).reshape(shape)
    return np.fft.ifft(np.fft.fft(values * demodulator, axis=0) * phase, axis=0) / demodulator


def interval_pieces(interval, shift, period=PERIOD):
    start, end = map(float, interval)
    span = end - start
    if not 0. < span <= period:
        raise ValueError('window must have positive span no larger than the carrier')
    moved = (start + float(shift)) % period
    if moved + span <= period:
        return [(moved, moved + span)]
    return [(moved, period), (0., moved + span - period)]


def translated_integral(grid, values, interval, shift):
    return sum(model.interval_integral(grid, values, piece) for piece in interval_pieces(interval, shift, grid.length))


def _readonly(values):
    values = np.array(values, copy=True)
    values.setflags(write=False)
    return values


def translate_pair_state(pair, state, shift):
    """Move basis functions, sources, references, field and observer worldlines."""
    # The saved quadrature products then translate by an exact node permutation.
    if abs(shift * pair.grid.nf / pair.grid.length - round(shift * pair.grid.nf / pair.grid.length)) > 1e-10:
        raise ValueError('whole-equation control requires an offset on the fermion lattice')
    ap = lambda values: _readonly(antiperiodic_shift(values, shift, pair.grid.length))
    phi0, phi1 = ap(pair.original_phi0), ap(pair.original_phi1)
    source0, source1 = ap(pair.source_phi0), ap(pair.source_phi1)
    moved = replace(pair, geometry_map=_readonly(periodic_shift(pair.geometry_map, shift, pair.grid.length)),
                    original_phi0=phi0, original_phi1=phi1, original_columns=_readonly(np.vstack((phi0, phi1))),
                    source_phi0=source0, source_phi1=source1, source_columns=_readonly(np.vstack((source0, source1))),
                    clock_locations=tuple((x + shift) % pair.grid.length for x in pair.clock_locations),
                    parent_interval=(pair.parent_interval[0] + shift, pair.parent_interval[1] + shift),
                    child_interval=(pair.child_interval[0] + shift, pair.child_interval[1] + shift),
                    geometry_metadata=dict(pair.geometry_metadata, translation_shift=float(shift)),
                    source_metadata=dict(pair.source_metadata, translation_shift=float(shift)))
    changed = state.copy()
    changed.phi0 = antiperiodic_shift(state.phi0, shift, pair.grid.length)
    changed.phi1 = antiperiodic_shift(state.phi1, shift, pair.grid.length)
    return moved, changed


def _gap(actual, expected, *, reference_scale=None):
    actual, expected = np.asarray(actual), np.asarray(expected)
    error = float(np.max(np.abs(actual - expected)))
    scale = max(1., float(np.max(np.abs(expected))))
    if reference_scale is not None:
        scale = max(scale, float(reference_scale))
    return {'absolute_gap': error, 'scale': scale, 'relative_to_scale': error / scale,
            'within_roundoff_tolerance': error <= 2e-9 + 2e-9 * scale}


def _diagnostics(pair, state, jets, shift=0., original_pair=None):
    fine, grid = jets['fine'], pair.grid
    source = jets['bundle']['source']
    radial = fine.r * fine.Q
    mass = np.sum((np.abs(fine.phi0)**2 + np.abs(fine.phi1)**2) * pair.weights[None, :], axis=1) / grid.dx_q
    normal = source['force_L'] / fine.r / grid.dx_q
    original_pair = pair if original_pair is None else original_pair
    windows = {}
    for name, interval in [('child', original_pair.child_interval), ('parent', original_pair.parent_interval)]:
        windows[name] = {
            'pieces': interval_pieces(interval, shift, grid.length),
            'proper_length': translated_integral(grid, radial, interval, shift),
            'probability': translated_integral(grid, mass, interval, shift),
            'normal_energy': translated_integral(grid, normal, interval, shift),
        }
    return {
        'windows': windows,
        'clock_locations': list(pair.clock_locations),
        'clock_rates': model._periodic_values(grid, radial, pair.clock_locations).tolist(),
        'weighted_trace': float(np.sum(mass) * grid.dx_q),
        'reference_amplitudes': pair.original_columns.conj().T @ np.vstack((state.phi0, state.phi1)),
        'source_amplitudes': pair.source_columns.conj().T @ np.vstack((state.phi0, state.phi1)),
    }


def compare_translation(pair, state, shift=DEFAULT_SHIFT, *, step=None):
    before_state = episode.state_sha256(state)
    before_basis = episode.array_sha256(pair.geometry_map)
    before_source = episode.array_sha256(pair.source_columns)
    moved, shifted = translate_pair_state(pair, state, shift)
    original_jets = tidal.analytic_accelerations(pair, state, 'coupled')
    shifted_jets = tidal.analytic_accelerations(moved, shifted, 'coupled')
    checks = {}
    original_nodal = model.reconstruct_state(pair, state)
    shifted_nodal = model.reconstruct_state(moved, shifted)
    original_rate, original_bundle = model.rates(pair, state, return_bundle=True)
    shifted_rate, shifted_bundle = model.rates(moved, shifted, return_bundle=True)
    for name in model.GEOMETRY_NAMES + model.MOMENTUM_NAMES:
        checks['physical_state_' + name] = _gap(getattr(shifted_nodal, name), periodic_shift(getattr(original_nodal, name), shift))
        checks['canonical_rate_' + name] = _gap(getattr(shifted_rate, name), getattr(original_rate, name))
        checks['fine_prolongation_' + name] = _gap(getattr(shifted_jets['fine'], name), periodic_shift(getattr(original_jets['fine'], name), shift))
    for name in ('phi0', 'phi1'):
        checks['field_rate_' + name] = _gap(getattr(shifted_rate, name), antiperiodic_shift(getattr(original_rate, name), shift))
        checks['AP_prolongation_' + name] = _gap(getattr(shifted_jets['fine'], name), antiperiodic_shift(getattr(original_jets['fine'], name), shift))
    checks['basis_orthogonality'] = _gap(moved.geometry_map.T @ moved.geometry_map, np.eye(pair.grid.ng))
    for name in ('force_Q', 'force_L', 'force_beta'):
        checks[name] = _gap(shifted_bundle['source'][name], periodic_shift(original_bundle['source'][name], shift))
    for name in ('Q_t', 'r_t', 'Q_tt', 'r_tt'):
        checks[name] = _gap(shifted_jets[name], periodic_shift(original_jets[name], shift))
    original_constraints = coupling.constraint_residuals(original_bundle['fine_system'], original_bundle['fine_state'], original_bundle['source'])
    shifted_constraints = coupling.constraint_residuals(shifted_bundle['fine_system'], shifted_bundle['fine_state'], shifted_bundle['source'])
    for name in ('hamilton', 'momentum', 'rho', 'current'):
        # Constraints subtract large geometric and matter terms. Compare the
        # covariance gap against those terms, not a near-zero residual alone.
        natural_scale = max(float(np.max(np.abs(original_constraints['rho']))), float(np.max(np.abs(original_constraints['current']))))
        checks['constraint_' + name] = _gap(shifted_constraints[name], periodic_shift(original_constraints[name], shift), reference_scale=natural_scale)
    original_profiles = tidal.profiles_on_grid(pair, original_jets)
    shifted_profiles = tidal.profiles_on_grid(moved, shifted_jets)
    tidal_values = {}
    for name in ('R_h', 'R4', 'R_0101', 'R_0202', 'owned_W'):
        original, changed = original_profiles['tides'][name], shifted_profiles['tides'][name]
        checks['actual_curvature_' + name] = _gap(changed, periodic_shift(original, shift))
        tidal_values[name] = {'original_abs_max': float(np.max(np.abs(original))), 'translated_abs_max': float(np.max(np.abs(changed)))}
    first = _diagnostics(pair, state, original_jets)
    second = _diagnostics(moved, shifted, shifted_jets, shift, pair)
    for name in ('clock_rates', 'weighted_trace', 'reference_amplitudes', 'source_amplitudes'):
        checks['diagnostic_' + name] = _gap(second[name], first[name])
    for window in first['windows']:
        for name in ('proper_length', 'probability', 'normal_energy'):
            checks[window + '_' + name] = _gap(second['windows'][window][name], first['windows'][window][name])
    if step is not None:
        if not 0. < float(step) <= .001:
            raise ValueError('optional control is limited to one RK4 step at most 0.001')
        original_after = model.rk4_step(pair, state, step)
        changed_after = model.rk4_step(moved, shifted, step)
        for name in model.STATE_NAMES:
            expected = antiperiodic_shift(getattr(original_after, name), shift) if name in ('phi0', 'phi1') else getattr(original_after, name)
            checks['one_step_' + name] = _gap(getattr(changed_after, name), expected)
    preserved = (episode.state_sha256(state) == before_state and episode.array_sha256(pair.geometry_map) == before_basis and episode.array_sha256(pair.source_columns) == before_source)
    first.pop('reference_amplitudes'); first.pop('source_amplitudes')
    second.pop('reference_amplitudes'); second.pop('source_amplitudes')
    return {'shift': float(shift), 'checks': checks, 'ok': preserved and all(item['within_roundoff_tolerance'] for item in checks.values()),
            'input_state_basis_source_preserved': preserved, 'original': first, 'translated': second,
            'actual_tidal_scalars': tidal_values, 'optional_RK4_step': step,
            'basis_interpretation': 'W columns translated; canonical coefficients and parent/detail indices preserved',
            'clock_worldlines_shifted': True, 'spinor_boundary': 'antiperiodic half-integer Fourier shift'}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assessment(directory=DEFAULT_EPISODE, *, shift=DEFAULT_SHIFT, step=None):
    directory = Path(directory)
    started = time.process_time()
    producer_paths = (LAB / 'src/recursive_horizons/nsc_discovery_translation.py', LAB / 'scripts/derive_nsc_discovery_translation.py')
    producers_before = tidal.producer_hashes(producer_paths)
    hashes = {}
    cases = []
    for nf in (128, 256):
        case_id = f'nf{nf}_coupled_dt0.0005'
        for ordinal, station in ((0, .3), (1, 1.)):
            if time.process_time() - started >= CPU_LIMIT:
                raise RuntimeError('translation assessment exhausted its 30 CPU-second budget')
            stem = directory / f'{case_id}-{ordinal:06d}'
            for extension in ('.json', '.npz'):
                path = Path(str(stem) + extension)
                hashes[str(path.relative_to(LAB))] = sha256(path)
            record, arrays = episode.load_checkpoint(directory, case_id, ordinal)
            if abs(float(record['coordinate_time']) - station) > 1e-10:
                raise ValueError('saved ordinal is not the requested coordinate station')
            pair = backend.make_fft_pair(episode.pair_from_arrays(arrays, record))
            state = episode.state_from_arrays(arrays, record['momentum_representation'])
            compared = compare_translation(pair, state, shift, step=step)
            compared.update(case_id=case_id, time=station, nf=nf)
            compared['saved_proper_clocks_at_shifted_worldlines'] = arrays['normal_clocks'].tolist()
            compared['clock_history_scope'] = 'same saved accumulated proper clocks assigned to the translated physical worldlines; no new clock integration'
            cases.append(compared)
    preserved = all(sha256(LAB / path) == digest for path, digest in hashes.items())
    producers = tidal.producer_hashes(producer_paths)
    if producers != producers_before:
        raise RuntimeError('translation producer source changed during assessment')
    producers = {str(Path(path).relative_to(LAB)): digest for path, digest in producers.items()}
    elapsed = time.process_time() - started
    if elapsed > CPU_LIMIT:
        raise RuntimeError('translation assessment exceeded its 30 CPU-second budget')
    return {'schema': SCHEMA, 'ok': preserved and all(case['ok'] for case in cases), 'cases': cases,
            'input_hashes': hashes, 'producer_hashes': producers, 'hashes_unchanged': preserved,
            'cpu_seconds': elapsed, 'cpu_limit_seconds': CPU_LIMIT,
            'scope': 'whole-state coordinate/seam covariance on the existing periodic carrier; not an ambient-extent or echo control',
            'physical_mechanism_claimed_from_identity': False, 'continuum_certified': False,
            'trajectory_evolved': step is not None}
