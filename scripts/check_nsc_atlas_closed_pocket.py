"""Reproduce one jointly constrained closed parent-to-compact preparation.

Five independent, spawned runs retain the same declared source:64/96 modes
in both time directions and one64-mode future time-step control. No membrane,
q2 parent, old state cache or selected quantum state is inserted.
"""
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.integrate import DOP853
from scipy.optimize import brentq
from scipy.special import roots_jacobi
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_valve_neck import solve_pocket_data
from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve, gegenbauer_jets
from recursive_horizons.nsc_atlas_passive_dop853 import PassiveDOP853Quadrature
from recursive_horizons.nsc_atlas_sphere_flux_completion import combined_record
from recursive_horizons.nsc_atlas_shape_transition import flat_horizon_pp_record


ROOT = Path(__file__).resolve().parents[1]
SELECTED = frozenset((-25., -20., -10., 0., 10., 20., 25.))
LEDGER_ORDER = ('north_carrier', 'north_radion', 'south_carrier', 'south_radion',
                'north_pressure_work', 'south_pressure_work')
PROTOCOL = (
    dict(name='N64_past', modes=64, endpoint=-25., max_step=.02),
    dict(name='N64_future', modes=64, endpoint=25., max_step=.02),
    dict(name='N96_past', modes=96, endpoint=-25., max_step=.02),
    dict(name='N96_future', modes=96, endpoint=25., max_step=.02),
    dict(name='N64_future_time_control', modes=64, endpoint=25., max_step=.01),
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_pins():
    package = (ROOT/'src/recursive_horizons').resolve()
    files = {Path(__file__).resolve()}
    for module in tuple(sys.modules.values()):
        name = getattr(module, '__file__', None)
        if name:
            path = Path(name).resolve()
            if path.suffix == '.py' and path.is_relative_to(package):
                files.add(path)
    return {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(files)}


def git_source_check(commit, sources):
    mismatches = []
    for name, expected in sources.items():
        old = subprocess.run(['git', 'show', commit+':'+name], cwd=ROOT, capture_output=True)
        if old.returncode or hashlib.sha256(old.stdout).hexdigest() != expected:
            mismatches.append(name)
    return mismatches


class CachedSpherical(SphericalValve):
    """Cache evaluation rules only; the package physical RHS is unchanged."""
    def __init__(self, modes, initial_data):
        super().__init__(modes, initial_data=initial_data, reflection_symmetry=False,
                         budget_seconds=45.)
        self.rules = {}

    def retain(self, points):
        points = np.asarray(points, dtype=float)
        self.rules[(points.shape, points.tobytes())] = gegenbauer_jets(points, self.modes)
        return points

    def jets(self, state, u=None):
        if u is None:
            return super().jets(state)
        u = np.atleast_1d(u)
        rule = self.rules.get((u.shape, u.tobytes()))
        if rule is None:
            return super().jets(state, u)
        coefficients = self.coefficients(state)
        return tuple(coefficients@matrix.T for matrix in rule)


class Measurements:
    """Physical fixed-cap fluxes, pressure work and independent ADM controls."""
    def __init__(self, system):
        self.system = system
        angle = np.arange(1, 402)*math.pi/402
        self.global_u = system.retain(np.cos(angle))
        self.global_w = 4*math.pi**2/402*np.sin(angle)**2
        self.poles = system.retain(np.array([1., 0., -1.]))
        x, w = roots_jacobi(48, .5, 0.)
        self.positive = .75+.25*x
        self.cap_u = system.retain(np.r_[self.positive, -self.positive, .5, -.5])
        self.cap_w = 4*math.pi*.125*w*np.sqrt(1+self.positive)

    def matter(self, state, *, details=False):
        z = self.system.pieces(state, self.cap_u)
        v = 1-self.cap_u**2
        radion = np.sqrt(v)/z['q']*z['y'][6]*z['yu'][4]
        carrier = np.sqrt(v)/z['q']*z['y'][7]*z['yu'][5]
        areas = 4*math.pi*z['radius'][-2:]**2
        rates = np.empty(6)
        rates[:4] = [-areas[0]*carrier[-2], -areas[0]*radion[-2],
                      areas[1]*carrier[-1], areas[1]*radion[-1]]
        caps = {}
        for index, name in enumerate(('north', 'south')):
            sl = slice(index*48, (index+1)*48)
            measure = self.cap_w*z['q'][sl]*z['b'][sl]**2
            rates[4+index] = measure@(z['ka'][sl]*z['pr'][sl]+2*z['kb'][sl]*z['pt'][sl])
            if details:
                energy = float(measure@z['rho'][sl])
                volume = float(sum(measure))
                pressure = float(measure@((z['pr'][sl]+2*z['pt'][sl])/3))
                s = self.system.sstar+z['y'][4, sl]/math.sqrt(84)
                rate = z['y'][6, sl]/math.sqrt(84)
                carrier_density=.5*(z['y'][7,sl]**2+(1-self.cap_u[sl]**2)/z['q'][sl]**2*z['yu'][5,sl]**2
                    +400*np.exp(-12*s)*z['y'][5,sl]**2)
                caps[name] = dict(energy=energy, volume=volume, density=energy/volume,
                                  volume_rate=float(measure@(-z['K'][sl])),
                                  mean_isotropic_pressure=pressure/volume,
                                  carrier_energy=float(measure@carrier_density),
                                  w=pressure/energy, s_range=[float(min(s)), float(max(s))],
                                  s_rate_range=[float(min(rate)), float(max(rate))])
        return (rates, caps) if details else rates

    def area_maximum(self, state):
        """Actual spatial jets/null signs, not curvature supplied by a source."""
        s = self.system
        def radial(coord):
            return float(s.pieces(state, [coord])['rr'][0])
        try:
            root = brentq(radial, -.5, .5, xtol=2e-13)
        except ValueError:
            return dict(found=False, bracket=[-.5, .5],
                        scope='No central sign-changing stationary sphere in the declared bracket.')
        z = s.pieces(state, [root])
        radius = float(z['radius'][0])
        second = -.5*float(z['Rr'][0])*radius
        plus, minus = (float(x[0]) for x in s.expansions(state, [root]))
        if plus < -1e-12 and minus < -1e-12:
            signs = 'future_trapped_cosmological_cut'
        elif plus > 1e-12 and minus > 1e-12:
            signs = 'anti_trapped_cosmological_cut'
        else:
            signs = 'marginal_or_mixed_null_signs'
        return dict(found=True, u=root, radius=radius,
                    radius_normal_rate=-float(z['kb'][0])*radius,
                    R_ll=second, spatial_maximum=second < 0,
                    theta_plus=plus, theta_minus=minus, null_signs=signs,
                    scope='Spatial/null geometry in the global closed cosmology; no black-hole, event-horizon or causal-child classification.')

    def snapshot(self, t, state, integrals, baseline):
        s = self.system
        z = s.pieces(state, self.global_u)
        v = 1-self.global_u**2
        R3 = z['Rr']+2*z['Rt']
        K2 = z['ka']**2+2*z['kb']**2
        ham = R3+z['K']**2-K2-2*s.kappa*z['rho']
        mom = 2*((self.global_u/np.sqrt(v)-np.sqrt(v)*z['yu'][0]/z['b'])*(z['ka']-z['kb'])
                 +np.sqrt(v)*z['yu'][2])/z['q']-s.kappa*z['j']
        scale = s.H**2+abs(R3)+z['K']**2+K2+2*s.kappa*abs(z['rho'])
        weight = self.global_w*z['q']*z['b']**2
        rates, caps = self.matter(state, details=True)
        for index, name in enumerate(('north', 'south')):
            c = caps[name]
            carrier, radion = integrals[2*index:2*index+2]
            work = integrals[4+index]
            residual = c['energy']-baseline[name]['energy']-carrier-radion-work
            c.update(carrier_received=float(carrier), radion_received=float(radion),
                     material_received=float(carrier+radion), pressure_work=float(work),
                     carrier_influx=float(rates[2*index]), radion_influx=float(rates[2*index+1]),
                     ledger_absolute=float(residual), ledger_relative=float(abs(residual)/max(1., abs(c['energy']))))
        poles = s.pieces(state, self.poles)
        row = dict(proper_time=float(t), Hrelative=float(max(abs(ham)/scale)),
                   Mrelative=float(max(abs(mom)/scale)),
                   Hamiltonian_RMS=float(np.sqrt(weight@(ham*ham)/sum(weight))),
                   momentum_RMS=float(np.sqrt(weight@(mom*mom)/sum(weight))),
                   positive_spatial_scale_min=float(min(np.min(z['q']), np.min(z['b']))),
                   angular_expansion_range=[float(min(-z['kb'])), float(max(-z['kb']))],
                   normal_volume_expansion_range=[float(min(-z['K'])), float(max(-z['K']))],
                   pole_order=['north', 'center', 'south'],
                   pole_s=(s.sstar+poles['y'][4]/math.sqrt(84)).tolist(),
                   pole_s_rate=(poles['y'][6]/math.sqrt(84)).tolist(),
                   pole_chi=poles['y'][5].tolist(), pole_Pi_chi=poles['y'][7].tolist(),
                   caps=caps)
        row['sampled_admitted'] = max(row['Hrelative'], row['Mrelative']) <= .02
        if t in SELECTED:
            row['central_area_maximum'] = self.area_maximum(state)
        return row


def evolve_worker(spec, raw_directory):
    """One fresh initialization and physical run; failures remain visible."""
    started = time.perf_counter()
    cpu = time.process_time()
    before = source_pins()
    rows, times, states, integral_rows = [], [], [], []
    data = None
    status, reason = 'incomplete', 'worker did not initialize'
    initial_hash = None
    accepted = 0.
    with threadpool_limits(limits=1):
        try:
            data = solve_pocket_data(40, epsilon=0., carrier_scale=.25, compression=.01,
                                     budget_seconds=15.)
            system = CachedSpherical(spec['modes'], data)
            measure = Measurements(system)
            initial = system.initial.copy()
            initial_hash = hashlib.sha256(initial.tobytes()).hexdigest()
            _, baseline = measure.matter(initial, details=True)
            ledgers = [PassiveDOP853Quadrature() for _ in range(6)]
            def rhs(t, flat):
                if time.perf_counter()-started > 45.:
                    raise TimeoutError('declared45s worker wall allowance exhausted')
                derivative = system.rhs(t, flat)
                for ledger, rate in zip(ledgers, measure.matter(flat.reshape(8, spec['modes']))):
                    ledger.record(rate)
                return derivative
            first = measure.snapshot(0., initial, np.zeros(6), baseline)
            rows.append(first); times.append(0.); states.append(initial); integral_rows.append(np.zeros(6))
            if not first['sampled_admitted']:
                raise RuntimeError('joint initial physical constraints failed the2% guard')
            solver = DOP853(rhs, 0., initial.ravel(), spec['endpoint'],
                            rtol=2e-10, atol=2e-12, max_step=spec['max_step'])
            for ledger in ledgers:
                ledger.initialize(solver)
            direction = 1 if spec['endpoint'] > 0 else -1
            target = .25*direction
            while solver.status == 'running':
                solver.step()
                if solver.status == 'failed':
                    raise RuntimeError('DOP853 failed')
                accepted = float(solver.t)
                for ledger in ledgers:
                    ledger.accepted_step(solver)
                if direction*(solver.t-target) >= 0:
                    dense = solver.dense_output()
                    state = dense(target).reshape(8, spec['modes'])
                    integrals = np.array([ledger.dense_value(solver, target) for ledger in ledgers])
                    row = measure.snapshot(target, state, integrals, baseline)
                    rows.append(row); times.append(target); states.append(state.copy()); integral_rows.append(integrals)
                    if not row['sampled_admitted']:
                        raise RuntimeError('sampled independent ADM2% guard ended')
                    target += .25*direction
            if not times or times[-1] != spec['endpoint']:
                raise RuntimeError('declared endpoint has no retained event state')
            status, reason = 'endpoint_reached', 'declared proper-time endpoint'
        except Exception as exc:
            reason = type(exc).__name__+': '+str(exc)
    after = source_pins()
    changed = [name for name, digest in before.items() if after.get(name) != digest]
    if changed:
        status, reason = 'source_changed', 'source bytes changed during worker: '+', '.join(changed)
    raw = Path(raw_directory)/(spec['name']+'.npz')
    payload = dict(times=np.asarray(times), states=np.asarray(states).reshape((-1, 8, spec['modes'])),
                   ledgers=np.asarray(integral_rows).reshape((-1, 6)),
                   diagnostics_json=np.array(json.dumps(rows, allow_nan=False)))
    if data is not None:
        payload.update(initial_coefficients=data.coefficients, initial_labels=data.labels,
                       initial_S_coefficients=data.S_coefficients, initial_S_labels=data.S_labels)
    np.savez_compressed(raw, **payload)
    caps = ('north', 'south')
    valid = [row for row in rows if row['sampled_admitted']]
    return dict(**spec, status=status, stop_reason=reason,
                rtol=2e-10, atol=2e-12, threads=1, wall_allowance_seconds=45.,
                elapsed_seconds=time.perf_counter()-started, cpu_seconds=time.process_time()-cpu,
                accepted_proper_time=accepted,
                last_retained_proper_time=times[-1] if times else None,
                last_admitted_proper_time=valid[-1]['proper_time'] if valid else None,
                initial_state_sha256=initial_hash,
                initial_constraint=dict(degree=data.degree, labels=data.labels.tolist(),
                    positivity_lower_bound=data.positivity_lower_bound,
                    continuation=data.continuation) if data is not None else None,
                max_constraints=dict(H=max((row['Hrelative'] for row in rows), default=None),
                                     M=max((row['Mrelative'] for row in rows), default=None)),
                max_ledger_relative={name: max((row['caps'][name]['ledger_relative'] for row in rows), default=None) for name in caps},
                selected_events=[row for row in rows if row['proper_time'] in SELECTED],
                endpoint_state=rows[-1] if rows and rows[-1]['proper_time'] not in SELECTED else None,
                raw_cache=raw.relative_to(ROOT).as_posix(), raw_cache_sha256=sha(raw),
                source_sha256=after, source_unchanged=not changed,
                scope='Finite sampled controls; retained unadmitted events remain diagnostic only. No raw accepted tail beyond the last event is archived.')


def _event(run, t):
    return next((row for row in run.get('selected_events', ()) if row['proper_time'] == t), None)


def _interval(first, last):
    if first is None or last is None or not first['sampled_admitted'] or not last['sampled_admitted']:
        return dict(available=False, reason='Both interval events must be retained and admitted')
    out = {}
    for name in ('north', 'south'):
        a, b = first['caps'][name], last['caps'][name]
        carrier = b['carrier_received']-a['carrier_received']
        radion = b['radion_received']-a['radion_received']
        work = b['pressure_work']-a['pressure_work']
        out[name] = dict(carrier=carrier, radion=radion, material=carrier+radion,
                         pressure_work=work, energy_change=b['energy']-a['energy'],
                         energy_identity_residual=b['energy']-a['energy']-carrier-radion-work)
    return dict(available=True, from_proper=first['proper_time'], to_proper=last['proper_time'], caps=out)


def comparison_records(runs):
    by = {run['name']: run for run in runs}
    initials = {}
    for modes in (64, 96):
        past, future = by[f'N{modes}_past'], by[f'N{modes}_future']
        arrays = []
        for run in (past, future):
            if run.get('raw_cache'):
                with np.load(ROOT/run['raw_cache'], allow_pickle=False) as raw:
                    if len(raw['states']):
                        arrays.append(raw['states'][0].copy())
        initials[str(modes)] = dict(available=len(arrays) == 2,
            bitwise_identical=len(arrays) == 2 and np.array_equal(*arrays),
            hashes=[past.get('initial_state_sha256'), future.get('initial_state_sha256')])
    controls = [by['N64_future'], by['N64_future_time_control']]
    arrays = []
    for run in controls:
        if run.get('raw_cache'):
            with np.load(ROOT/run['raw_cache'], allow_pickle=False) as raw:
                if len(raw['states']):
                    arrays.append(raw['states'][0].copy())
    initials['time_control64'] = dict(available=len(arrays) == 2,
        bitwise_identical=len(arrays) == 2 and np.array_equal(*arrays),
        hashes=[run.get('initial_state_sha256') for run in controls])
    full = {str(n): _interval(_event(by[f'N{n}_past'], -25.), _event(by[f'N{n}_future'], 25.)) for n in (64, 96)}
    late = {name: _interval(_event(by[name], 20.), _event(by[name], 25.))
            for name in ('N64_future', 'N96_future', 'N64_future_time_control')}
    differences = {}
    for label, first, second in (
            ('spatial64_vs96', 'N64_future', 'N96_future'),
            ('time_step64', 'N64_future', 'N64_future_time_control')):
        a, b = late[first], late[second]
        differences[label] = (dict(available=True,
            cap_carrier_spread={cap: abs(a['caps'][cap]['carrier']-b['caps'][cap]['carrier']) for cap in ('north', 'south')},
            scope='Empirical differences, not rigorous error bounds')
            if a['available'] and b['available'] else dict(available=False, reason='A required interval failed admission'))
    return dict(initial_states=initials, physical_past_to_future=full,
                compact_interval20_to25=late, numerical_difference_indicators=differences)


def calculation(raw_directory, workers=5):
    started, cpu = time.perf_counter(), time.process_time()
    before = source_pins()
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    mismatches = git_source_check(commit, before)
    if mismatches:
        raise RuntimeError('Commit the actual producer/physical sources before this calculation: '+', '.join(mismatches))
    raw_directory.mkdir(parents=True, exist_ok=False)
    results = {}
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('spawn')) as pool:
        pending = {pool.submit(evolve_worker, dict(spec), str(raw_directory)): spec for spec in PROTOCOL}
        try:
            for future in as_completed(pending, timeout=120.):
                spec = pending[future]
                try:
                    results[spec['name']] = future.result()
                except Exception as exc:
                    results[spec['name']] = dict(**spec, status='worker_failure',
                        stop_reason=type(exc).__name__+': '+str(exc), selected_events=[], endpoint_state=None,
                        source_sha256={}, source_unchanged=False)
        except TimeoutError:
            for future, spec in pending.items():
                if spec['name'] not in results:
                    future.cancel()
                    results[spec['name']] = dict(**spec, status='coordinator_timeout',
                        stop_reason='Coordinator120s allowance ended; no completed result admitted',
                        selected_events=[], endpoint_state=None, source_sha256={}, source_unchanged=False)
    runs = [results[spec['name']] for spec in PROTOCOL]
    after = source_pins()
    sources = dict(after)
    inconsistent = []
    for run in runs:
        for name, digest in run['source_sha256'].items():
            if name in sources and sources[name] != digest:
                inconsistent.append(name)
            sources[name] = digest
    bound_mismatches = git_source_check(commit, sources)
    unchanged = (all(after.get(name) == digest for name, digest in before.items())
                 and all(run['source_unchanged'] for run in runs) and not inconsistent
                 and all(sha(ROOT/name) == digest for name, digest in sources.items()))
    comparison = comparison_records(runs)
    complete = all(run['status'] == 'endpoint_reached' for run in runs)
    identical = all(row['bitwise_identical'] for row in comparison['initial_states'].values())
    return dict(schema='nsc-atlas-closed-pocket-v1', scientific_owner='docs/nsc-atlas-valve-formation.md',
        action=combined_record().action,
        flat_horizon_comparison=flat_horizon_pp_record(),
        preparation=dict(q=1.5, epsilon_radion=0., carrier_scale=.25, compression=.01,
            initial_factory='solve_pocket_data(40,epsilon=0,carrier_scale=.25,compression=.01)',
            initial_s='Instantaneous force barrier at actualchi; no suppliedY2 offset',
            chi='Achi*u', Pi_chi_conformal='-Bchi*sqrt(3/5)*(u²-1/6)',
            Pi_S_conformal='-sqrt(84)*.01', trace_K=0., conformal_TT=0.,
            full_conformal_basis='All GegenbauerU_l,l=0..40; scalarYork+carrierYork cross term retained',
            reflection_projection=False, membrane_inserted=False, flux_change=0.,
            q2_parent_transplanted=False, benchmark_amplitudes_selected=False),
        ledger_convention=dict(order=list(LEDGER_ORDER),
            north='u>=.5; inward fixed-cap material rate=-4piR²j at u=.5',
            south='u<=-.5; inward fixed-cap material rate=+4piR²j at u=-.5',
            j='sqrt(1-u²)/q*(PiS*S_u+PiChi*chi_u), in the increasingchi direction',
            work='integral(Ka*pr+2Kb*pt)dV; separate from material transfer',
            signed='All raw ledgers integrate0→t. Physical−25→+25 equals future endpoint minus backward endpoint.',
            interfaces='Fixed coordinate caps, not horizons or causal child boundaries'),
        protocol=dict(runs=list(PROTOCOL), spatial_validation_order=401, cap_Jacobi_order=48,
            event_spacing=.25, selected_events=sorted(SELECTED), independent_ADM_relative_limit=.02,
            rtol=2e-10, atol=2e-12, worker_wall_seconds=45., coordinator_wall_allowance_seconds=120.),
        runs=runs, comparison=comparison,
        verdict=dict(all_five_runs_completed=complete, shared_initial_states_bitwise_identical=identical,
            finite_history_numerically_admitted=complete and identical and unchanged and not bound_mismatches,
            child_proved=False, offspring_proved=False, causal_separation_proved=False,
            state_selection_proved=False, SI_scale_selected=False,
            scope='One conditional closed parent-to-compact benchmark with regional stress-energy transfer; not a black-hole birth, selected quantum state or complete past.'),
        provenance=dict(producer_commit=commit, source_sha256=sources, source_unchanged=unchanged,
            source_bound_to_commit=not bound_mismatches, source_commit_mismatches=bound_mismatches,
            inconsistent_worker_sources=sorted(set(inconsistent))),
        compute=dict(workers=workers, BLAS_threads_per_worker=1,
            elapsed_seconds=time.perf_counter()-started, coordinator_cpu_seconds=time.process_time()-cpu,
            sum_worker_cpu_seconds=sum(run.get('cpu_seconds', 0.) for run in runs),
            scope='Measured process wall/CPU; worker allowances include fresh constraint initialization and output.'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-directory', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=5, choices=range(1, 6))
    args = parser.parse_args()
    raw = (ROOT/args.raw_directory).resolve()
    output = (ROOT/args.output).resolve()
    if not raw.is_relative_to(ROOT/'.build') or raw.exists():
        raise ValueError('Fresh raw directory under.build required')
    if not output.is_relative_to(ROOT) or output.exists() or output == raw:
        raise ValueError('Fresh in-repository summary path required')
    result = calculation(raw, args.workers)
    result['reproduction'] = subprocess.list2cmdline([
        'python', 'scripts/run_python.py', Path(__file__).relative_to(ROOT).as_posix(),
        '--raw-directory', raw.relative_to(ROOT).as_posix(),
        '--output', output.relative_to(ROOT).as_posix(), '--workers', str(args.workers)])
    def exact(value):
        if isinstance(value,Fraction):return str(value)
        raise TypeError(type(value).__name__)
    content = json.dumps(result, indent=2, allow_nan=False,default=exact)+'\n'
    if len(content.encode()) >= 100_000:
        raise RuntimeError('Compact scientific summary exceeded100KB')
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(content, encoding='utf-8', newline='\n')
    print(json.dumps(dict(output=str(output), verdict=result['verdict'], compute=result['compute']), indent=2))
    return 0 if result['verdict']['finite_history_numerically_admitted'] else 2


if __name__ == '__main__':
    raise SystemExit(main())
