"""Redshift, light focusing and radar clock distance in saved spherical geometry.

The fixed-coordinate probes are normal observers, possibly accelerated.
No homogeneous expansion law, field solve or cosmological parameter fit is used.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from measure_nsc_conformal_geometric_average import (
    ROOT, EPISODE, CURVATURE, EXPECTED_PAYLOAD, CASES, PRIMARY,
    restore_geometry, sha256,
)
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_cosmological_distance import radial_conformal_lightcone
from recursive_horizons.nsc_coupled_local_response import hermite_with_slopes
from recursive_horizons.nsc_spherical_galerkin_coupling import periodic_interpolation


PERIOD = 8.
RECEPTION_TIME = .3
RECEPTION_X = 1.
LOOKBACKS = (.025, .075, .125)
ORDERS = (32, 64)
CONSISTENCY_TOLERANCE = 1e-10


class FourierHermiteMetric:
    """The original real Fourier band and supplied nodal time derivatives."""

    def __init__(self, times, radius, Q, radius_rate, Q_rate, *, length, budget):
        self.times = np.asarray(times, dtype=float)
        fields = [np.asarray(value, dtype=float) for value in (radius, Q, radius_rate, Q_rate)]
        if (any(value.ndim != 2 or value.shape != fields[0].shape for value in fields)
                or fields[0].shape[0] != len(self.times)
                or not all(np.isfinite(value).all() for value in fields)):
            raise ValueError('finite aligned stored metric values and rates required')
        self.count = fields[0].shape[1]
        self.length = float(length)
        self.budget = budget
        if self.count % 2 != 1 or self.length <= 0:
            raise ValueError('the owned geometry band is odd with positive period')
        coefficients = [np.fft.fft(value, axis=1)/self.count for value in fields]

        def pack(a, b):
            return np.concatenate((a.real, a.imag, b.real, b.imag), axis=1)

        self.values = pack(coefficients[0], coefficients[1])
        self.rates = pack(coefficients[2], coefficients[3])
        self.wave = 2*np.pi*np.fft.fftfreq(self.count, d=self.length/self.count)

    def _spatial(self, packed, positions):
        n = self.count
        phase = np.exp(1j*np.atleast_1d(positions)[:, None]*self.wave[None, :])
        r_coeff = packed[:, :n]+1j*packed[:, n:2*n]
        q_coeff = packed[:, 2*n:3*n]+1j*packed[:, 3*n:]
        r = np.sum(r_coeff*phase, axis=1)
        q = np.sum(q_coeff*phase, axis=1)
        if max(float(np.max(abs(r.imag))), float(np.max(abs(q.imag)))) > 1e-10:
            raise ValueError('real periodic metric left its Fourier subspace')
        return r.real, q.real

    def __call__(self, coordinate_time, position):
        marks, locations = np.broadcast_arrays(np.asarray(coordinate_time, dtype=float),
                                               np.asarray(position, dtype=float))
        self.budget.check(4*self.count*marks.size)
        packed = hermite_with_slopes(self.times, self.values, self.rates, marks.ravel())
        radius, Q = self._spatial(packed, np.mod(locations.ravel(), self.length))
        if np.any(radius <= 0) or np.any(Q <= 0):
            raise ValueError('interpolated metric left the positive chart')
        return {'radius': radius.reshape(marks.shape), 'Q': Q.reshape(marks.shape)}


def observer_radar_distances(metric, time_nodes, *, reception_time, reception_x,
                             lookbacks, quadrature_order):
    """Half the same observer's proper round-trip duration to fixed mirrors."""
    largest = 2*max(lookbacks)
    if reception_time-largest < float(time_nodes[0])-1e-12:
        raise ValueError('radar trip leaves the saved metric interval')
    cuts = sorted({0., *(2*float(value) for value in lookbacks),
                   *(reception_time-float(mark) for mark in time_nodes
                     if 0 < reception_time-float(mark) < largest), largest})
    left, right = np.asarray(cuts[:-1]), np.asarray(cuts[1:])
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    marks = (left+right)[:, None]/2+(right-left)[:, None]*nodes[None, :]/2
    geometry = metric(reception_time-marks.ravel(), np.full(marks.size, reception_x))
    lapse = (geometry['radius']*geometry['Q']).reshape(marks.shape)
    integrals = (right-left)/2*(lapse@weights)
    cumulative = np.r_[0., np.cumsum(integrals)]
    indices = {value: index for index, value in enumerate(cuts)}
    return [float(cumulative[indices[2*float(value)]]/2) for value in lookbacks]


def source_bindings():
    names = (
        'src/recursive_horizons/nsc_atlas_cosmological_distance.py',
        'src/recursive_horizons/nsc_coupled_local_response.py',
        'src/recursive_horizons/nsc_spherical_galerkin_coupling.py',
        'scripts/measure_nsc_conformal_geometric_average.py',
        'scripts/measure_nsc_conformal_lightcone.py',
        'tests/test_nsc_atlas_cosmological_distance.py',
        'tests/test_nsc_conformal_lightcone.py',
    )
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    bindings = {}
    for name in names:
        raw = (ROOT/name).read_bytes()
        if raw != subprocess.check_output(['git', 'show', commit+':'+name], cwd=ROOT):
            raise ValueError('light-cone calculation source is not frozen: '+name)
        bindings[name] = sha256(ROOT/name)
    return commit, bindings


def metric_from_payload(payload, case, interpolation, budget):
    times = payload[case+'_frame_times']
    ng, nq = interpolation.shape[1], interpolation.shape[0]
    radius, Q = (payload[case+'_frames_'+key] for key in ('r', 'Q'))
    radius_rate, Q_rate = ((ng/nq)*payload[case+'_frames_rate_'+key]@interpolation
                           for key in ('r', 'Q'))
    return FourierHermiteMetric(times, radius, Q, radius_rate, Q_rate,
                                length=PERIOD, budget=budget)


def normalize_point(point, radar, reference_length):
    return {
        'direction': int(point['direction']), 'coordinate_lookback': float(point['coordinate_lookback']),
        'emission_time': float(point['emission_time']), 'emission_x': float(point['emission_x']),
        'redshift': float(point['redshift']),
        'observer_lapse_over_Lref': float(point['observer_lapse'])/reference_length,
        'emission_lapse_over_Lref': float(point['emission_lapse'])/reference_length,
        'angular_distance_over_Lref': float(point['angular_diameter_distance'])/reference_length,
        'luminosity_distance_over_Lref': float(point['luminosity_distance'])/reference_length,
        'radar_distance_over_Lref': float(radar)/reference_length,
        'luminosity_over_radar': float(point['luminosity_distance'])/radar,
    }


def measure(geometry_path, budget, record):
    commit, bindings = source_bindings()
    episode = json.loads(EPISODE.read_text())
    curvature = json.loads(CURVATURE.read_text())
    if (sha256(geometry_path) != EXPECTED_PAYLOAD
            or episode['payload_sha256'] != EXPECTED_PAYLOAD
            or sha256(EPISODE) != curvature['source_bindings']['episode_v2']
            or bindings['src/recursive_horizons/nsc_spherical_galerkin_coupling.py']
               != curvature['source_bindings']['galerkin']):
        raise ValueError('the original metric or projection is not authenticated')
    record.update({
        'schema': 'nsc-conformal-lightcone-summary-v1', 'status': 'incomplete',
        'source_commit': commit, 'source_sha256': bindings,
        'input': {'episode_sha256': sha256(EPISODE), 'payload_sha256': EXPECTED_PAYLOAD,
                  'payload_bytes': geometry_path.stat().st_size},
        'primary': PRIMARY,
        'protocol': {'reception_time': RECEPTION_TIME, 'reception_x': RECEPTION_X,
                     'coordinate_lookbacks': list(LOOKBACKS), 'directions': [-1, 1],
                     'probe_worldlines': 'fixed-coordinate normals, possibly accelerated',
                     'light': 'transparent photon-conserving geometric optics',
                     'distance_ratio': 'luminosity distance / observer radar distance'},
        'field_solves': 0, 'physical_time_steps': 0, 'rows': [],
        'limits': {'cpu_seconds': 60., 'charged_values': 90_000_000,
                   'dimensionless_quadrature_consistency_tolerance': CONSISTENCY_TOLERANCE},
        'scope': ('Finite causal observations in the separately specified conformal common-action '
                  'metric. No geodesic matter congruence, full sky, cosmological Hubble inference, '
                  'P001 parent or generated Lambda is claimed.'),
    })
    with np.load(geometry_path, allow_pickle=False) as payload:
        maps = {nf: periodic_interpolation(nf-1, 4*nf, PERIOD) for nf in (256, 512)}
        initial_radius = maps[512]@payload[PRIMARY+'_frames_r'][0]
        Lref = float(initial_radius[int(RECEPTION_X/PERIOD*len(initial_radius))])
        record['reference_length'] = {'value_in_original_action_units': Lref,
            'definition': 'initial areal radius of the original observer sphere at T=0.05,x=1'}
        for case in CASES:
            nf = int(case.split('_')[0][2:])
            metric = metric_from_payload(payload, case, maps[nf], budget)
            at_observer = metric(RECEPTION_TIME, RECEPTION_X)
            clock = float(at_observer['radius']*at_observer['Q'])
            original_clock = float(episode['results'][case]['final']['clock_rate'])
            if abs(clock-original_clock)/Lref > CONSISTENCY_TOLERANCE:
                raise ValueError('original observer clock changed')
            calculations = {}
            for order in ORDERS:
                radar = observer_radar_distances(metric, metric.times,
                    reception_time=RECEPTION_TIME, reception_x=RECEPTION_X,
                    lookbacks=LOOKBACKS, quadrature_order=order)
                for direction in (-1, 1):
                    cone = radial_conformal_lightcone(metric, metric.times,
                        reception_time=RECEPTION_TIME, reception_x=RECEPTION_X,
                        lookbacks=LOOKBACKS, direction=direction, quadrature_order=order)
                    calculations[(order, direction)] = [
                        normalize_point(point, clock_distance, Lref)
                        for point, clock_distance in zip(cone, radar)]
            checked_keys = ('redshift', 'angular_distance_over_Lref', 'luminosity_distance_over_Lref',
                            'radar_distance_over_Lref', 'luminosity_over_radar')
            for direction in (-1, 1):
                for coarse, reference in zip(calculations[(32, direction)], calculations[(64, direction)]):
                    gaps = {key: coarse[key]-reference[key] for key in checked_keys}
                    record['rows'].append({'case': case, 'observables': reference,
                                           'order32_minus64': gaps})
                    if max(abs(value) for value in gaps.values()) > CONSISTENCY_TOLERANCE:
                        record.update(status='stopped', reason='quadrature consistency tolerance failed')
                        return record
            budget.check()
    record['primary_rows'] = [row for row in record['rows'] if row['case'] == PRIMARY]
    record['refinement_indicators'] = []
    for direction in (-1, 1):
        for lookback in LOOKBACKS:
            rows = {row['case']: row['observables'] for row in record['rows']
                    if row['observables']['direction'] == direction
                    and row['observables']['coordinate_lookback'] == lookback}
            record['refinement_indicators'].append({
                'direction': direction, 'coordinate_lookback': lookback,
                'coarser_time_minus_primary': {key: rows[CASES[2]][key]-rows[PRIMARY][key] for key in checked_keys},
                'coarser_space_minus_primary': {key: rows[CASES[1]][key]-rows[PRIMARY][key] for key in checked_keys},
                'four_case_maximum_deviation_from_primary': {
                    key: max(abs(row[key]-rows[PRIMARY][key]) for row in rows.values()) for key in checked_keys},
            })
    record['original_finite_constraint_defects'] = {case: {
        key: episode['results'][case]['final'][key] for key in ('full_hamilton_max', 'full_momentum_max')}
        for case in CASES}
    record['uncertainty_scope'] = (
        'Quadrature and original time/space refinement indicators only; no bound on Hermite interpolation, '
        'continuum state, physical model or observed data. No observational dataset was fitted.')
    record['reproduce'] = ('python scripts/run_python.py scripts/measure_nsc_conformal_lightcone.py '
        '--restore-geometry --geometry .build/data/nsc-spherical-conformal-episode-v2.npz '
        '--output results/summaries/nsc-conformal-lightcone-v1.json')
    budget.check()
    record['status'] = 'completed_finite_lightcone'
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--restore-geometry', action='store_true')
    args = parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('choose a new output in an existing directory')
    budget = Budget(maximum_evaluations=90_000_000, cpu_seconds=60.)
    record = {'schema': 'nsc-conformal-lightcone-summary-v1', 'status': 'incomplete',
              'field_solves': 0, 'physical_time_steps': 0}
    try:
        if args.restore_geometry:
            restore_geometry(args.geometry)
        measure(args.geometry, budget, record)
    except (ValueError, RuntimeError, FloatingPointError, OSError) as error:
        record.update(status='stopped', reason=str(error))
    record['execution'] = {'cpu_seconds': time.process_time()-budget.started,
                            'charged_values': budget.evaluations}
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+'\n',
                          encoding='utf-8', newline='\n')
    print(json.dumps({'status': record['status'], 'execution': record['execution'],
                      'reason': record.get('reason')}))
    return 0 if record['status'] == 'completed_finite_lightcone' else 1


if __name__ == '__main__':
    raise SystemExit(main())
