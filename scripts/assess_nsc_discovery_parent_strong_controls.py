#!/usr/bin/env python3
"""Read authenticated strong/empty controls; reuse frozen columns without evolution."""
from __future__ import annotations
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
from threadpoolctl import threadpool_limits
import assess_nsc_discovery_parent_episode as saved

LIMIT = 64*1024*1024
DEFAULTS = tuple(saved.episode.LAB/'results/development'/('nsc-discovery-parent-'+name+'-v1')
                 for name in ('strong-t1', 'strong-t3', 'strong-empty-k', 'frozen-parent-heavy'))
EQUAL = ('Q', 'r', 'W', 'phi0', 'phi1', 'source_phi0', 'source_phi1', 'observer_columns')


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def exact_initial(strong, frozen):
    for name in EQUAL:
        if not np.array_equal(strong[name], frozen[name]):
            raise ValueError('frozen reuse initial equality failed: '+name)
    weights = np.asarray(strong['source_weights'])/np.asarray(frozen['source_weights'])
    if weights.shape != (2,) or not np.isfinite(weights).all() or np.any(weights <= 0):
        raise ValueError('invalid column reweight')
    return weights


def reweight(tagged, ratio):
    child = np.asarray(tagged['child'], dtype=float)*ratio
    total = np.asarray(tagged['ambient_total'], dtype=float)*ratio
    return {'child_probability': float(child.sum()), 'probability_total': float(total.sum()),
            'child_fraction': None if total.sum() == 0 else float(child.sum()/total.sum())}


def interpolate(clock, values, target):
    clock, values = np.asarray(clock, float), np.asarray(values, float)
    if len(clock) < 2 or values.shape[0] != len(clock) or not np.isfinite(clock).all() or np.any(np.diff(clock) <= 0):
        raise ValueError('proper clocks must be finite and strictly increasing')
    if not np.isfinite(values).all() or not np.isfinite(target):
        raise ValueError('nonfinite clock interpolation')
    if target < clock[0] or target > clock[-1]:
        return None
    return float(np.interp(target, clock, values))


def anchored_clocks(rows, stations):
    """Positive-rate trapezoids, scaled within each interval to exact saved clocks."""
    t = np.asarray([r['time'] for r in rows]); rate = np.asarray([r['centre_rate'] for r in rows])
    if np.any(np.diff(t) <= 0) or np.any(rate <= 0) or not np.isfinite(rate).all():
        raise ValueError('scalar time/positive clock check failed')
    tau = np.zeros(len(rows)); discrepancy = 0.
    for left, right in zip(stations[:-1], stations[1:]):
        i, j = np.searchsorted(t, [left['time'], right['time']])
        if j >= len(t) or t[i] != left['time'] or t[j] != right['time']:
            raise ValueError('scalar stream lacks a saved clock anchor')
        delta = right['centre_tau']-left['centre_tau']
        if delta <= 0:
            raise ValueError('saved proper clocks are not increasing')
        integral = np.r_[0., np.cumsum(np.diff(t[i:j+1])*(rate[i:j]+rate[i+1:j+1])/2)]
        discrepancy = max(discrepancy, abs(float(integral[-1])-delta))
        tau[i:j+1] = left['centre_tau']+integral*delta/integral[-1]
    if np.any(np.diff(tau) <= 0):
        raise ValueError('scalar stream exceeds saved clock domain')
    return tau, discrepancy


def load(directory, hashes):
    directory = Path(directory).resolve()
    def track(path):
        if path.stat().st_size > LIMIT:
            raise ValueError('input must be <=64 MiB: '+str(path))
        hashes[str(path)] = saved.file_hash(path)
    manifest_path = directory/'manifest.json'; track(manifest_path)
    manifest = json.loads(manifest_path.read_text()); saved.authenticate_producers(manifest, pinned=None)
    if manifest.get('status') in ('RUNNING', 'PREPARED'):
        raise ValueError('control campaign is not completed')
    cases = {}
    for case in manifest['cases']:
        name = case['case_id']; chunks = []
        for path in sorted(directory.glob(name+'-*.json')):
            track(path); track(path.with_suffix('.npz'))
            record, arrays = saved.authenticate_snapshot(directory, path)
            if record['case_id'] != name:
                raise ValueError('snapshot case mismatch')
            saved.authenticate_producers(record, pinned=None)
            for rel in record.get('input_hashes', {}):
                hashes.setdefault(str(saved.episode.LAB/rel), saved.file_hash(saved.episode.LAB/rel))
            chunks.append((record, arrays))
        if not chunks:
            raise ValueError('no authenticated snapshots: '+name)
        path = directory/(name+'-observations.jsonl'); track(path)
        rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
        cases[name] = {'chunks': chunks, 'rows': rows, 'status': manifest['status']}
    return cases


def station(record, arrays):
    pair, state = saved.reconstruct(record, arrays)
    fine = saved.leading.fine_state(pair, state); grid = pair.grid
    index = int(np.argmin(np.abs(np.asarray(pair.clock_locations)-grid.length/2)))
    mass = (abs(fine.phi0)**2+abs(fine.phi1)**2)*pair.weights/grid.dx_q
    child = [float(saved.extent.real_interval_integral(grid, mass[:,j], pair.child_interval)) for j in range(2)]
    total = float(np.sum(mass)*grid.dx_q)
    mode = 'frozen_geometry' if record['control_mode'] == 'frozen_geometry' else 'coupled'
    rate, bundle = saved.leading.rates(pair, state, return_bundle=True, control_mode=mode)
    curvature, _, _ = saved.leading.metric_jets(pair, state, rate, bundle, control_mode=mode)
    return {'time': float(record['coordinate_time']), 'centre_tau': float(arrays['normal_clocks'][index]),
            'parent_k': float(record['parent_k']), 'child_length': float(saved.extent.real_interval_integral(grid, fine.r*fine.Q, pair.child_interval)),
            'centre_radius': float(saved.extent.real_periodic_values(grid, fine.r, [grid.length/2])[0]),
            'Nmin': float(np.min(fine.r*fine.Q)), 'child_probability': sum(child), 'probability_total': total,
            'child_fraction': None if total == 0 else sum(child)/total, 'column_child': child,
            'energy': saved.leading.energy(pair, state),
            'constraints': {k:v for k,v in saved.leading.constraints(pair, state).items() if k in ('raw_C_max','D_max','projected_h_c_max','held_out_h_c_max','observable_error_bound')},
            'projected_R4_max_abs': float(np.max(abs(curvature['tides']['R4']))),
            'snapshot_kind': record['snapshot_kind']}


def scalar_rows(rows, ratio=None):
    result = []
    for row in rows:
        observer = row['external_observer']['result']; locations = np.asarray(observer['proper_clocks']['locations'])
        centre = int(np.argmin(abs(locations-4.)))
        p = observer['channels']['probability']; total = float(p['total'])
        values = reweight(row['tagged_column_probability'], ratio) if ratio is not None else {
            'child_probability': float(p['child']), 'probability_total': total,
            'child_fraction': None if total == 0 else float(p['child'])/total}
        result.append(dict(time=float(row['time']), centre_rate=float(row['normal_clock_rates'][centre]),
                           **({} if ratio is not None else {'energy_total': float(row['energy']['total'])}), **values))
    return result


def report(directories=DEFAULTS):
    # Bind the actual evaluator and all imported numerical owners, not a historical label.
    modules = [saved, saved.backend, saved.episode, saved.extent, saved.leading, saved.observer,
               saved.tidal, saved.galerkin, saved.provenance]
    modules += [m for n,m in sys.modules.items() if n.startswith('recursive_horizons.') and getattr(m, '__file__', None)]
    sources = {str(Path(m.__file__).resolve()): saved.file_hash(m.__file__) for m in modules}
    sources[str(Path(__file__).resolve())] = saved.file_hash(__file__)
    hashes = {}; groups = [load(d, hashes) for d in directories]
    strong, empty, frozen = {}, groups[2], next(iter(groups[3].values()))
    for group in groups[:2]:
        for name, case in group.items():
            target = strong.setdefault(name, {'chunks': [], 'rows': []})
            target['chunks'] += case['chunks']; target['rows'] += case['rows']
    tables, matched, direct = {}, [], []
    with threadpool_limits(limits=1), saved.backend.fft_thread_limit(1):
        strong_initial = next(c for c in strong.values())['chunks'][0]
        frozen_initial = frozen['chunks'][0]; ratio = exact_initial(strong_initial[1], frozen_initial[1])
        for name in ('nf', 'child_interval', 'parent_interval', 'clock_locations'):
            if strong_initial[0][name] != frozen_initial[0][name]:
                raise ValueError('frozen carrier/observer mismatch: '+name)
        pair_s, state_s = saved.reconstruct(*strong_initial); pair_f, state_f = saved.reconstruct(*frozen_initial)
        system_s = saved.leading.active_system(pair_s.grid, saved.leading.fine_state(pair_s, state_s))
        system_f = saved.leading.active_system(pair_f.grid, saved.leading.fine_state(pair_f, state_f))
        for n in ('length_density', 'shift', 'kappa'):
            if not np.array_equal(getattr(system_s,n), getattr(system_f,n)):
                raise ValueError('frozen Dirac generator coefficient differs: '+n)
        if type(system_s.momentum) is not type(system_f.momentum) or any(not np.array_equal(getattr(system_s.momentum,n), getattr(system_f.momentum,n)) for n in ('points','length','_symbol','_twist','_untwist')):
            raise ValueError('frozen Fourier momentum differs')
        rate_s = saved.leading.rates(pair_s, state_s, control_mode='frozen_geometry')
        rate_f = saved.leading.rates(pair_f, state_f, control_mode='frozen_geometry')
        if not all(np.array_equal(getattr(rate_s,n), getattr(rate_f,n)) for n in ('phi0','phi1')):
            raise ValueError('frozen initial Dirac generator action differs')
        for record, arrays in frozen['chunks']:
            for n in ('Q','r','W','pi_Q','pi_r'):
                if not np.array_equal(arrays[n], frozen_initial[1][n]):
                    raise ValueError('frozen geometry/canonical array changed: '+n)
            row = station(record, arrays); t = row['time']
            sample = next(r for r in frozen['rows'] if r['time'] == t)
            if not sample['external_observer']['result']['rates']['geometry_jets_zero']:
                raise ValueError('frozen saved geometry rates are not zero')
            weighted = reweight(sample['tagged_column_probability'], ratio)
            direct_child = float(np.dot(row['column_child'], ratio))
            if not np.isclose(direct_child, weighted['child_probability'], rtol=1e-11, atol=1e-13):
                raise ValueError('scalar reweight disagrees with direct NPZ column integrals')
            direct.append({'time': t, **weighted, 'direct_NPZ_child_gap': direct_child-weighted['child_probability'],
                           'centre_tau': row['centre_tau']})
        fs = scalar_rows(frozen['rows'], ratio); ftau, fgap = anchored_clocks(fs, direct)
        for name, case in strong.items():
            sign = case['chunks'][0][0]['sign_name']; ec = next(c for n,c in empty.items() if '_'+sign+'_' in n)
            s0, e0 = case['chunks'][0], ec['chunks'][0]
            if not np.array_equal(exact_initial(s0[1], frozen_initial[1]), ratio):
                raise ValueError('strong signs have different frozen reweights')
            
            if any(not np.array_equal(s0[1][n], e0[1][n]) for n in ('Q','r','W','source_weights','observer_columns')):
                raise ValueError('empty initial geometry/weights/observer mismatch')
            if s0[0]['parent_k'] != e0[0]['parent_k'] or any(np.any(e0[1][n]) for n in ('phi0','phi1')):
                raise ValueError('empty control is not matched-k source-free')
            bytime = {}
            for record, arrays in case['chunks']:
                t = float(record['coordinate_time'])
                if t in bytime and record['arrays_sha256'] != bytime[t][0]['arrays_sha256']:
                    raise ValueError('strong successor changed its exact handoff')
                bytime[t] = (record, arrays)
            stations = [station(*pair) for _,pair in sorted(bytime.items())]
            empty_stations = [station(*pair) for pair in ec['chunks']]
            tables[sign] = {'occupied': stations, 'same_k_empty': empty_stations,
                            'pi_r_initial_gap': float(np.max(abs(s0[1]['pi_r']-e0[1]['pi_r']))),
                            'force_deletion_on_same_Cauchy_state': False, 'empty_status': ec['status']}
            rows = {r['time']: r for r in scalar_rows(case['rows'])}; ss = [r for _,r in sorted(rows.items())]
            stau, sgap = anchored_clocks(ss, stations)
            es = scalar_rows(ec['rows']); etau, egap = anchored_clocks(es, empty_stations)
            coarse = np.array([i for i,r in enumerate(fs) if abs(r['time']*10-round(r['time']*10)) < 1e-7])
            for row in stations:
                tau = row['centre_tau']; f = interpolate(ftau, [r['child_fraction'] for r in fs], tau)
                fc = interpolate(ftau[coarse], [fs[i]['child_fraction'] for i in coarse], tau)
                matched.append({'sign': sign, 'coordinate_time': row['time'], 'centre_tau': tau,
                    'coupled_fraction': row['child_fraction'], 'coupled_child_length': row['child_length'],
                    'coupled_centre_radius': row['centre_radius'], 'frozen_fraction': f,
                    'frozen_coordinate_time': interpolate(ftau, [r['time'] for r in fs], tau),
                    'coupled_minus_frozen': None if f is None else row['child_fraction']-f,
                    'frozen_0_1_sampling_gap': None if f is None or fc is None else fc-f,
                    'empty_coordinate_time': interpolate(etau, [r['time'] for r in es], tau),
                    'empty_energy_total': interpolate(etau, [r['energy_total'] for r in es], tau),
                    'empty_child_length_station_interpolation': interpolate([r['centre_tau'] for r in empty_stations], [r['child_length'] for r in empty_stations], tau),
                    'empty_centre_radius_station_interpolation': interpolate([r['centre_tau'] for r in empty_stations], [r['centre_radius'] for r in empty_stations], tau)})
            tables[sign]['clock_anchor_trapezoid_indicator'] = {'occupied': sgap, 'empty': egap}
    body = {'schema': 'NSC-DISCOVERY-PARENT-STRONG-CONTROLS-v1', 'stations': tables,
            'frozen_reweight': {'ratios': ratio.tolist(), 'stored_times': direct, 'generator_action_equal': True,
                               'scalar_selected_times': [r for r in fs if any(abs(r['time']-t)<1e-8 for t in (1.,1.25,1.5,2.25,3.))],
                               'clock_anchor_trapezoid_indicator': fgap}, 'matched_centre_clock': matched,
            'input_sha256': hashes, 'measurement_source_sha256': sources,
            'limitations': ['Re-solved empty-source C/D; same initial Q,r and parent_k; canonical momenta differ.',
                'Frozen work, aggregate energy and constraints are not reused under changed weights/momenta.',
                'Clock interpolation and 0.1 sampling are finite indicators without an interpolation error bound.',
                'Source-free matter fractions are undefined (null); no negative filled vacuum is inferred.',
                'Projected leading curvature/constraints are diagnostics; no continuum instability or autonomous renewal proof.'],
            'evolved': False, 'extrapolated': False, 'renewal_asserted': False,
            'next_physical_question': 'Does source-dependent transfer persist at matched proper clocks under a resolved spatial/time refinement?'}
    for path, digest in {**hashes, **sources}.items():
        if saved.file_hash(path) != digest:
            raise ValueError('input or evaluator changed during assessment: '+path)
    body['content_sha256'] = hashlib.sha256(canonical(body).encode()).hexdigest()
    return body


def write_record(path, payload, directories=DEFAULTS):
    path = Path(path).expanduser().resolve()
    if any(path.is_relative_to(Path(d).resolve()) for d in directories):
        raise ValueError('output must be outside sealed input directories')
    encoded = (json.dumps(payload, indent=2, allow_nan=False)+'\n').encode()
    if len(encoded) > LIMIT:
        raise ValueError('output exceeds 64 MiB')
    with path.open('xb') as stream:
        stream.write(encoded)
    path.chmod(0o444)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument('--output', type=Path)
    args = parser.parse_args(argv)
    try:
        if args.output and args.output.exists():
            raise FileExistsError('refusing to overwrite assessment output')
        result = report()
        if args.output:
            write_record(args.output, result)
        print(json.dumps(result, indent=2, allow_nan=False)); return 0
    except (ValueError, RuntimeError, OSError, KeyError) as error:
        print(type(error).__name__+': '+str(error), file=sys.stderr); return 2


if __name__ == '__main__':
    raise SystemExit(main())
