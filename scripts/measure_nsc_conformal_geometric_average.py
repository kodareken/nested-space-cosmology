"""Lapse-weighted averages of the authenticated, already evolved conformal metric.

This is the separate Einstein/Weyl/Dirac realization, not P001's absent parent.
The source NPZ stays in an ignored input cache. No field or geometry evolves.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import struct
import time
import urllib.request
import zlib

import numpy as np

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_geometric_average import periodic_spherical_average
from recursive_horizons.nsc_spherical_coupling import locked_coefficients, PERIOD
from recursive_horizons.nsc_spherical_feedback_action import partial_F, feedback_Z
from recursive_horizons.nsc_spherical_galerkin_coupling import periodic_interpolation


ROOT = Path(__file__).resolve().parents[1]
EPISODE = ROOT/'results/development/nsc-spherical-conformal-episode-v2.json'
CURVATURE = ROOT/'results/development/nsc-spherical-conformal-curvature-v1.json'
EXPECTED_PAYLOAD = '40a2c377e2dbf9af83f38670fc36c7029a73c7ae7367dd6cd95e4e93857aa877'
CASES = ('nf256_dt_0.0010', 'nf256_dt_0.0005', 'nf512_dt_0.0010', 'nf512_dt_0.0005')
PRIMARY = CASES[-1]
TIMES = (.05, .2, .3)
OBSERVATION_X = 1.
IDENTITY_TOLERANCE = 1e-10
ARCHIVE_URL = 'https://zenodo.org/api/records/23090015/files/kodareken/nested-space-cosmology-v0.27.0.zip/content'
ARCHIVE_MEMBER = 'kodareken-nested-space-cosmology-eea83b2/lab/results/development/nsc-spherical-conformal-episode-v2.npz'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def restore_geometry(target):
    """Recover only the immutable ZIP member, bounded to 60 MB of response data."""
    if target.exists():
        if sha256(target) != EXPECTED_PAYLOAD:
            raise ValueError('existing geometry cache differs; refusing to overwrite it')
        return
    archive_bytes, offset = 2_070_378_031, 1_942_517_559
    compressed, uncompressed = 53_041_551, 53_100_722
    charged = 0

    def read_range(start, size):
        nonlocal charged
        if charged+size > 60_000_000:
            raise RuntimeError('archive recovery exceeded its 60 MB transport limit')
        request = urllib.request.Request(ARCHIVE_URL, headers={
            'Range': f'bytes={start}-{start+size-1}', 'Accept-Encoding': 'identity'})
        with urllib.request.urlopen(request, timeout=40) as response:
            content_range = f'bytes {start}-{start+size-1}/{archive_bytes}'
            if response.status != 206 or response.headers.get('Content-Range') != content_range:
                raise ValueError('bounded archive range unavailable; full archive not downloaded')
            raw = response.read(size+1)
        if len(raw) != size:
            raise ValueError('archive range length mismatch')
        charged += len(raw)
        return raw

    header = struct.unpack('<4s5H3L2H', read_range(offset, 30))
    if header[0] != b'PK\x03\x04' or header[3] != 8:
        raise ValueError('archive member header changed')
    name = read_range(offset+30, header[-2]).decode('utf-8')
    if name != ARCHIVE_MEMBER:
        raise ValueError('archive member name changed')
    packed = read_range(offset+30+header[-2]+header[-1], compressed)
    decoder = zlib.decompressobj(-15)
    raw = decoder.decompress(packed, uncompressed+1)
    if (len(raw) != uncompressed or not decoder.eof or decoder.unused_data
            or hashlib.sha256(raw).hexdigest() != EXPECTED_PAYLOAD):
        raise ValueError('recovered geometry differs from the original payload')
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(raw)


def snapshot_jets(interpolation, arrays, index, *, A, C_W):
    """Differentiate the actual projected bilinear velocities, not sample times."""
    ng = interpolation.shape[1]
    nq = interpolation.shape[0]

    def lift(name):
        return interpolation@arrays[name][index]

    def project(value):
        return interpolation@((ng/nq)*(interpolation.T@value))

    radius, Q = lift('r'), lift('Q')
    pr, pc = lift('p_r'), lift('p_chi')
    rt, Qt = arrays['rate_r'][index], arrays['rate_Q'][index]
    prt, pct = arrays['rate_p_r'][index], arrays['rate_p_chi'][index]
    Fr, Fchi = partial_F(radius, A, C_W)
    Fchi, Z = float(Fchi), float(feedback_Z(A))
    Fr_rate = -8*np.pi*A*rt
    first_Q = project(Q*pc/(2*Fchi))
    first_r = project((pr-Fr*pc/Fchi)/(2*Z))
    Qtt = project((Qt*pc+Q*pct)/(2*Fchi))
    rtt = project((prt-(Fr_rate*pc+Fr*pct)/Fchi)/(2*Z))
    radial = radius*Q
    radial_t = rt*Q+radius*Qt
    radial_tt = rtt*Q+2*rt*Qt+radius*Qtt
    return {
        'radius': radius, 'radial_scale': radial, 'lapse': radial,
        'radius_rate': rt, 'radial_scale_rate': radial_t, 'lapse_rate': radial_t,
        'radius_second_rate': rtt, 'radial_scale_second_rate': radial_tt,
    }, {
        'Q_rate': float(np.max(abs(first_Q-Qt))),
        'radius_rate': float(np.max(abs(first_r-rt))),
        'stored_lapse': float(np.max(abs(Q-arrays['L'][index]))),
        'stored_lapse_rate': float(np.max(abs(Qt-arrays['L_dot'][index]))),
    }


def _source_bindings():
    names = (
        'src/recursive_horizons/nsc_atlas_geometric_average.py',
        'src/recursive_horizons/nsc_spherical_feedback_action.py',
        'src/recursive_horizons/nsc_spherical_coupling.py',
        'src/recursive_horizons/nsc_spherical_galerkin_coupling.py',
        'scripts/measure_nsc_conformal_geometric_average.py',
        'tests/test_nsc_atlas_geometric_average.py',
        'tests/test_nsc_conformal_geometric_average.py',
    )
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    bindings = {}
    for name in names:
        raw = (ROOT/name).read_bytes()
        if raw != subprocess.check_output(['git', 'show', commit+':'+name], cwd=ROOT):
            raise ValueError('calculation source is not frozen: '+name)
        bindings[name] = hashlib.sha256(raw).hexdigest()
    return commit, bindings


def dimensionless_average(row, reference_length):
    """Use one independently selected proper ruler for every stored case."""
    inverse_square = (
        'Q', 'mean_intrinsic_curvature', 'kappa_mean_energy', 'kappa_mean_pressure',
        'kappa_geometric_density', 'kappa_geometric_pressure', 'geometric_trace_projection',
        'geometric_vacuum_condition', 'three_volume_acceleration', 'total_trace_projection',
        'hamiltonian_residual', 'scalar_fitting_acceleration_residual',
        'lapse_acceleration', 'direct_volume_acceleration',
    )
    L = reference_length
    normalized = {key: float(row[key])*L**2 for key in inverse_square}
    normalized.update(H_observer=float(row['H_observer'])*L,
                      volume_over_reference_length_cubed=float(row['volume'])/L**3,
                      observer_lapse_over_reference_length=float(row['observer_lapse'])/L)
    checks = {
        'volume_acceleration': (row['three_volume_acceleration']-row['direct_volume_acceleration'])*L**2,
        'curvature_integration_by_parts': row['curvature_ibp_residual']/L,
        'lapse_integration_by_parts': row['lapse_ibp_residual']/L,
        'hamiltonian_metric_definition': row['hamiltonian_residual']*L**2,
    }
    return normalized, {key: float(value) for key, value in checks.items()}


def measure(geometry, budget, record):
    source_commit, bindings = _source_bindings()
    episode = json.loads(EPISODE.read_text())
    old_curvature = json.loads(CURVATURE.read_text())
    old_bindings = old_curvature['source_bindings']
    if (sha256(EPISODE) != old_bindings['episode_v2']
            or sha256(geometry) != EXPECTED_PAYLOAD
            or episode['payload_sha256'] != EXPECTED_PAYLOAD
            or episode['gauge'] != 'conformal'
            or episode['all_target_reached'] is not True
            or episode['source_changed'] is not False
            or episode['midtrajectory_reset'] is not False):
        raise ValueError('the original evolved geometry is not authenticated')
    for key, path in (('galerkin', 'src/recursive_horizons/nsc_spherical_galerkin_coupling.py'),
                      ('action', 'src/recursive_horizons/nsc_spherical_feedback_action.py')):
        if bindings[path] != old_bindings[key]:
            raise ValueError('original action or projection owner changed')
    if PERIOD != 8.:
        raise ValueError('original period changed')
    couplings = locked_coefficients()
    A, C_W = float(couplings['A']), float(couplings['C_W'])
    record.update({
        'schema': 'nsc-conformal-geometric-average-summary-v1',
        'status': 'incomplete', 'source_commit': source_commit, 'source_sha256': bindings,
        'input': {'episode_sha256': sha256(EPISODE), 'payload_sha256': EXPECTED_PAYLOAD,
                  'payload_bytes': geometry.stat().st_size,
                  'archive_url': ARCHIVE_URL, 'archive_member': ARCHIVE_MEMBER},
        'domain': {'coordinate_period': PERIOD, 'times': list(TIMES), 'observer_x': OBSERVATION_X,
                   'topology': 'S1 x S2', 'zero_shift': True},
        'primary': PRIMARY, 'rows': [], 'physical_time_steps': 0, 'field_solves': 0,
        'kinetic_couplings': {'A': A, 'C_W': C_W, 'Z': float(feedback_Z(A)),
                             'F_chi': float(partial_F(1., A, C_W)[1])},
        'scope': ('Lapse-weighted scalar fitting of the separate already evolved '
                  'Einstein/Weyl/Dirac metric. The Einstein source is metric-defined; '
                  'its physical interpretation inherits the saved finite constraint defects. '
                  'No independent matter-constraint certification, '
                  'P001 parent, measured fractal dimension or cosmological Lambda.'),
        'limits': {'cpu_seconds': 30., 'charged_values': 20_000_000,
                   'dimensionless_identity_tolerance': IDENTITY_TOLERANCE},
    })
    with np.load(geometry, allow_pickle=False) as payload:
        maps = {nf: periodic_interpolation(nf-1, 4*nf, PERIOD) for nf in (256, 512)}
        index0 = int(np.argmin(abs(payload[PRIMARY+'_frame_times']-TIMES[0])))
        primary_radius = maps[512]@payload[PRIMARY+'_frames_r'][index0]
        observer_primary = int(OBSERVATION_X/PERIOD*len(primary_radius))
        reference_length = float(primary_radius[observer_primary])
        record['reference_length'] = {
            'value_in_original_action_units': reference_length,
            'definition': 'Areal radius of the original observer sphere x=1 at T=0.05 in the primary case',
            'selected_by_target_observation': False,
        }
        for case in CASES:
            nf = int(case.split('_')[0][2:])
            interpolation = maps[nf]
            times = payload[case+'_frame_times']
            arrays = {name: payload[case+'_frames_'+name] for name in
                      ('Q', 'r', 'p_r', 'p_chi', 'rate_Q', 'rate_r', 'rate_p_r',
                       'rate_p_chi', 'L', 'L_dot')}
            for requested in TIMES:
                index = int(np.argmin(abs(times-requested)))
                if abs(float(times[index])-requested) > 1e-12:
                    raise ValueError('predeclared stored time is missing')
                jets, rate_gaps = snapshot_jets(interpolation, arrays, index, A=A, C_W=C_W)
                budget.check(40*interpolation.shape[0])
                normalized_rate_gaps = dict(rate_gaps)
                normalized_rate_gaps['radius_rate'] /= reference_length
                if max(normalized_rate_gaps.values()) > IDENTITY_TOLERANCE:
                    record.update(status='stopped', reason='saved projected first-rate reconstruction disagrees',
                                  failed_first_rate={'case': case, 'time': requested,
                                                     'dimensionless_differences': normalized_rate_gaps})
                    return record
                result = periodic_spherical_average(
                    areal_radius=jets['radius'], radial_scale=jets['radial_scale'], lapse=jets['lapse'],
                    radius_rate=jets['radius_rate'], scale_rate=jets['radial_scale_rate'],
                    lapse_rate=jets['lapse_rate'], radius_acceleration=jets['radius_second_rate'],
                    scale_acceleration=jets['radial_scale_second_rate'], coordinate_length=PERIOD,
                    observer_index=int(OBSERVATION_X/PERIOD*interpolation.shape[0]))
                source_row = min(episode['series'][case], key=lambda row: abs(row['time']-requested))
                normalized, checks = dimensionless_average(result, reference_length)
                checks['original_observer_clock'] = (result['observer_lapse']-source_row['clock_rate'])/reference_length
                record['rows'].append({
                    'case': case, 'time': requested, 'dimensionless': normalized,
                    'dimensionless_metric_identity_residuals': checks,
                    'first_rate_dimensionless_absolute_differences': normalized_rate_gaps,
                    'original_finite_constraint_defects': {
                        key: source_row[key] for key in ('full_hamilton_max', 'full_momentum_max')},
                })
                if max(abs(value) for value in checks.values()) > IDENTITY_TOLERANCE:
                    record.update(status='stopped', reason='metric identity tolerance failed')
                    return record
                budget.check()
    record['primary_rows'] = [row for row in record['rows'] if row['case'] == PRIMARY]
    comparison_keys = ('geometric_trace_projection', 'geometric_vacuum_condition',
                       'kappa_geometric_density', 'kappa_geometric_pressure', 'H_observer')
    record['refinement_indicators'] = []
    for requested in TIMES:
        rows = {row['case']: row['dimensionless'] for row in record['rows'] if row['time'] == requested}
        record['refinement_indicators'].append({
            'time': requested,
            'coarser_time_minus_primary': {key: rows[CASES[2]][key]-rows[PRIMARY][key] for key in comparison_keys},
            'coarser_space_minus_primary': {key: rows[CASES[1]][key]-rows[PRIMARY][key] for key in comparison_keys},
            'four_case_maximum_deviation_from_primary': {
                key: max(abs(row[key]-rows[PRIMARY][key]) for row in rows.values()) for key in comparison_keys},
        })
    budget.check()
    record['status'] = 'completed_saved_geometry_average'
    record['uncertainty_scope'] = 'Existing time/space refinement differences only; not a certified propagated or cosmological uncertainty.'
    record['reproduce'] = ('python scripts/run_python.py scripts/measure_nsc_conformal_geometric_average.py '
                           '--geometry .build/data/nsc-spherical-conformal-episode-v2.npz '
                           '--output results/summaries/nsc-conformal-geometric-average-v1.json')
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--geometry', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--restore-geometry', action='store_true',
                        help='recover the exact archived NPZ into the specified input cache if missing')
    args = parser.parse_args()
    if args.output.exists() or not args.output.parent.is_dir():
        raise ValueError('choose a new output in an existing directory')
    budget = Budget(maximum_evaluations=20_000_000, cpu_seconds=30.)
    record = {'schema': 'nsc-conformal-geometric-average-summary-v1', 'status': 'incomplete',
              'physical_time_steps': 0, 'field_solves': 0}
    try:
        if args.restore_geometry:
            restore_geometry(args.geometry)
        budget.check()
        record = measure(args.geometry, budget, record)
    except (ValueError, RuntimeError, FloatingPointError, OSError) as error:
        record.update(status='stopped', reason=str(error))
    budget_elapsed = time.process_time()-budget.started
    record['execution'] = {'cpu_seconds': budget_elapsed, 'charged_values': budget.evaluations}
    args.output.write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+'\n',
                          encoding='utf-8', newline='\n')
    print(json.dumps({'status': record['status'], 'execution': record['execution'],
                      'reason': record.get('reason')}))
    return 0 if record['status'] == 'completed_saved_geometry_average' else 1


if __name__ == '__main__':
    raise SystemExit(main())
