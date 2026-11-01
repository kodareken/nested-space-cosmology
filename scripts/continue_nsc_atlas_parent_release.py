"""Continue the unchanged original parent from its admitted proper-.9 state.

No initializer, new preparation or physical equation is supplied. Absolute
coordinate/proper clocks and the complete Cauchy state are retained. Event
locations use actual accepted DOP853 dense polynomials. No child is inferred.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from numba import set_num_threads
from scipy.fft import set_workers
from scipy.integrate import DOP853
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.check_nsc_atlas_parent_valve import load_original_parent
from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
from recursive_horizons.nsc_atlas_parent_gh import PAIRS
from recursive_horizons.nsc_atlas_parent_scale_selector import VolumeParent
from recursive_horizons.nsc_atlas_s3_frame import point_values


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def json_scalar(value):
    """Preserve numeric types; never stringify unsupported scientific data."""
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError('unsupported JSON value: '+type(value).__name__)


class Budget:
    def __init__(self, function, started, seconds):
        self.function, self.started, self.seconds, self.enabled = function, started, seconds, True

    def __call__(self, t, y):
        if self.enabled and time.perf_counter()-self.started >= self.seconds:
            raise TimeoutError('declared continuation wall-time budget exhausted')
        return self.function(t, y)


def source_identity():
    paths = {Path(__file__).resolve(), ROOT/'scripts/check_nsc_atlas_parent_valve.py'}
    package = ROOT/'src/recursive_horizons'
    for module in tuple(sys.modules.values()):
        path = getattr(module, '__file__', None)
        if path and Path(path).suffix == '.py' and Path(path).resolve().is_relative_to(package):
            paths.add(Path(path).resolve())
    commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    pins = {p.relative_to(ROOT).as_posix(): sha(p) for p in sorted(paths)}
    for name, digest in pins.items():
        archived = subprocess.run(['git', 'show', commit+':'+name], cwd=ROOT, capture_output=True)
        if archived.returncode or hashlib.sha256(archived.stdout).hexdigest() != digest:
            raise ValueError('commit unchanged continuation/decoder sources before running: '+name)
    return commit, pins


def anchor(system, state):
    """Actual first jets at the tracked geodesic; no averaged scalar phase."""
    C = system.unpack(state)
    Y = point_values(state[system.dof:system.dof+4], system.basis.labels)
    value = C@Y
    g = np.empty((4, 4))
    for k, (a, b) in enumerate(PAIRS):
        g[a, b] = g[b, a] = value[k]
    hi = np.linalg.inv(g[1:, 1:]); shift = hi@g[0, 1:]
    N2 = -g[0, 0]+g[0, 1:]@shift
    if N2 <= 0 or min(np.linalg.eigvalsh(g[1:, 1:])) <= 0:
        raise ValueError('actual accepted observer metric is not regular Lorentzian')
    dq = np.vstack([value[22:24], [system.r*Y@(D@C[20:22].T) for D in system.basis.generators]])
    u = state[system.dof+4:system.dof+8]; norm = float(u@g@u)
    if norm >= 0 or u[0] <= 0:
        raise ValueError('actual tracked observer is not future timelike')
    normal = np.r_[1., -shift]/math.sqrt(N2)
    return dict(s=float(system.sstar+value[20]/math.sqrt(system.K)), chi=float(value[21]),
        uS=float(u@dq[:, 0]/math.sqrt(-norm)),
        nS=float(normal@dq[:, 0]), nChi=float(normal@dq[:, 1]),
        observer_unit_norm_error=abs(norm+1))


def snapshot(system, t, state, phase):
    row = system.diagnostic(t, state)
    row.update(anchor(system, state))
    J = system.jets(state)
    shift = np.einsum('nij,nj->ni', np.linalg.inv(J.g[:, 1:, 1:]), J.g[:, 0, 1:])
    lapse = np.sqrt(-J.g[:, 0, 0]+np.einsum('ni,ni->n', shift, J.g[:, 0, 1:]))
    normal_pi = (J.p-np.einsum('ni,nia->na', shift, J.dq))/lapse[:, None]
    row['whole_grid_normal_momentum_ranges'] = {
        name: [float(min(normal_pi[:system.N, i])), float(max(normal_pi[:system.N, i]))]
        for i, name in enumerate(('S', 'chi'))}
    row['passive_sampled_phase'] = phase
    row['observer_sampled_unwrapped_phase'] = phase
    row['numerical_guard_passed'] = (row['metric_valid']
        and row['ADM_H_RMS_over_H_squared'] <= .1 and row['GH_C_RMS_over_H'] <= .1)
    return row


def run(args):
    started = time.perf_counter()
    folder = (ROOT/args.output).resolve()
    if folder.exists() or not folder.is_relative_to(ROOT/'.build'):
        raise ValueError('fresh output directory under .build required')
    set_num_threads(args.threads)
    with threadpool_limits(limits=1), set_workers(2):
        data = load_original_parent(args.degree)
        system, t0 = data['system'], float(data['end'])
        state0 = data['history'].evaluate(t0).copy()
        if abs(state0[system.dof+8]-.9) > 2e-10:
            raise ValueError('restart must preserve the admitted proper-.9 event')
        system.X0 = state0[system.dof:system.dof+4].copy()
        system.started, system.budget = started, math.inf
        commit, pins = source_identity()
        folder.mkdir(parents=True)
        top = next(row['s'] for row in VolumeParent().vacuum_roots() if row['kind']=='maximum')
        fold = -.5*math.log(brentq(lambda x: 225*x**6-560*x+336, .99, 1.))
        phase_row = next(row for row in data['provenance']['auxiliary_phase']['values']
                         if abs(row['reference_observer_proper_duration']-.9) < 2e-10)
        phase = phase_row['sparse_sampled_phase']
        angle = system.observer_angle(state0)
        result = dict(schema='nsc-atlas-parent-release-continuation-v1', degree=args.degree,
            status='running', producer_commit=commit, source_sha256=pins, source_bound_to_commit=True,
            original_parent=data['provenance'], start_coordinate_t=t0, start_proper=.9,
            requested_proper_endpoint=args.proper_endpoint, budget_seconds=args.budget,
            threads=dict(Numba=args.threads, FFT=2, BLAS=1),
            temporal_tolerances=dict(rtol=args.rtol, atol=args.atol, max_step=args.max_step),
            model_inputs=data['record']['model_inputs'], prescribed_GH_gauge=system.gauge,
            spatial_transform=system.backend, epsilon=system.epsilon,
            thresholds=dict(nonterminal_reference_fold_s=fold, terminal_bare_top_s=top,
                terminal_ADM_H_RMS_over_H2=.1, terminal_GH_C_RMS_over_H=.1),
            auxiliary_phase=dict(ODE_dummy=float(state0[-1]), sparse_phase_carry=phase,
                feedback=False, sparse_states_keep_dummy=True),
            start_full_state_sha256=hashlib.sha256(state0.tobytes()).hexdigest(),
            physical_history_prefix=dict(path=data['provenance']['dense_cache'], until_coordinate_t=t0),
            snapshots=[], events=[], child_formation_proved=False, escape_proved=False)
        meta = dict(degree=args.degree, epsilon=system.epsilon, gauge=system.gauge, backend=system.backend,
            quadrature_order=list(system.basis.quadrature_shape[:2]),
            original_requested_quadrature_order=data['history'].metadata.get('quadrature_order'),
            source_sha256=pins, producer_commit=commit, phase_dummy=float(state0[-1]))
        raw_metadata = json.dumps(meta, allow_nan=False, sort_keys=True, default=json_scalar)
        history, saved = DenseGHHistory(), []
        solver = None; last_state, last_t = state0, t0

        def save_record():
            result['elapsed_seconds'] = time.perf_counter()-started
            temp = folder/'run.tmp'
            temp.write_text(json.dumps(result, allow_nan=False, indent=2, default=json_scalar)+'\n', encoding='utf-8', newline='\n')
            temp.replace(folder/'run.json')

        def checkpoint(t, state, sampled_phase, event):
            row = snapshot(system, t, state, sampled_phase)
            row['event'] = event
            # Reject a nonfinite diagnostic before it can poison the terminal
            # JSON; the exact accepted state is still preserved by finally.
            json.dumps(row, allow_nan=False, default=json_scalar)
            result['snapshots'].append(row); saved.append((float(t), state.copy()))
            save_record()
            print(json.dumps(dict(degree=args.degree, proper=row['observer_proper_duration'], event=event,
                s=row['s'], uS=row['uS'], H_res=row['ADM_H_RMS_over_H_squared'], GH_res=row['GH_C_RMS_over_H'])), flush=True)
            if not row['numerical_guard_passed']:
                raise RuntimeError('sampled numerical admission guard ended; no physical endpoint inferred')

        status, reason = 'incomplete', None
        try:
            checkpoint(t0, state0, phase, 'same_Cauchy_restart')
            if anchor(system, state0)['s'] >= top:
                raise ValueError('bare top already reached before the requested continuation')
            budget = Budget(system.rhs, started, args.budget)
            solver = DOP853(budget, t0, state0, t_bound=8/system.r,
                rtol=args.rtol, atol=args.atol, max_step=args.max_step)
            targets = list(np.arange(1., args.proper_endpoint-1e-10, .1))+[args.proper_endpoint]
            seen_positive = anchor(system, state0)['uS'] > 0
            fold_recorded = anchor(system, state0)['s'] >= fold
            while solver.status == 'running':
                message = solver.step()
                if solver.status == 'failed':
                    raise RuntimeError(message or 'DOP853 failed')
                last_t, last_state = float(solver.t), solver.y.copy()
                budget.enabled = False
                try:
                    dense = solver.dense_output()
                finally:
                    budget.enabled = True
                history.record(solver.t_old, solver.t, dense)
                before_phase, before_angle = phase, angle
                angle = system.observer_angle(solver.y)
                if angle is not None and before_angle is not None and phase is not None:
                    phase += math.atan2(math.sin(angle-before_angle), math.cos(angle-before_angle))
                else:
                    phase = None
                def phase_at(t):
                    current = system.observer_angle(dense(t))
                    return (before_phase+math.atan2(math.sin(current-before_angle), math.cos(current-before_angle))
                        if current is not None and before_angle is not None and before_phase is not None else None)
                candidates = []
                probes = np.linspace(solver.t_old, solver.t, 5)
                points = [anchor(system, dense(t)) for t in probes]
                for left, right, a, b in zip(probes[:-1], probes[1:], points[:-1], points[1:]):
                    for name, level, terminal in (('reference_fold_s', fold, False), ('bare_top', top, True)):
                        if name == 'reference_fold_s' and fold_recorded:
                            continue
                        if a['s'] < level <= b['s']:
                            root = brentq(lambda t: anchor(system, dense(t))['s']-level, left, right, xtol=2e-12)
                            candidates.append((root, name, terminal))
                            if not terminal:
                                fold_recorded = True
                    if seen_positive and a['uS'] > 0 >= b['uS']:
                        root = brentq(lambda t: anchor(system, dense(t))['uS'], left, right, xtol=2e-12)
                        candidates.append((root, 'actual_uS_reversal', True))
                    seen_positive = seen_positive or b['uS'] > 0
                while targets and solver.y[system.dof+8] >= targets[0]:
                    target = targets.pop(0)
                    root = brentq(lambda t: dense(t)[system.dof+8]-target, solver.t_old, solver.t, xtol=2e-12)
                    terminal = bool(target == args.proper_endpoint)
                    candidates.append((root, 'proper_endpoint' if terminal else 'proper_checkpoint', terminal))
                for t, event, terminal in sorted(candidates):
                    point = dense(t)
                    result['events'].append(dict(event=event, coordinate_t=float(t),
                        proper=float(point[system.dof+8]), terminal=bool(terminal),
                        accepted_coordinate_bracket=[float(solver.t_old), float(solver.t)]))
                    checkpoint(t, point, phase_at(t), event)
                    if terminal:
                        status, reason = 'event_reached', event
                        result['physical_stop_coordinate_t'] = float(t)
                        break
                if status == 'event_reached':
                    break
            if status == 'incomplete':
                reason = 'coordinate allowance ended before requested proper endpoint'
        except (TimeoutError, RuntimeError, ValueError, TypeError, FloatingPointError) as exc:
            reason = str(exc)
            status = ('budget_stop' if isinstance(exc, TimeoutError) else
                      'export_stop' if isinstance(exc, TypeError) else 'numerical_stop')
        finally:
            if solver is not None:
                last_t, last_state = float(solver.t), solver.y.copy()
            # RAW FIRST: no diagnostic or public JSON call may prevent this.
            cache = folder/'states.npz'
            np.savez_compressed(cache, times=np.array([p[0] for p in saved]),
                states=np.array([p[1] for p in saved]), accepted_endpoint=last_state,
                accepted_endpoint_coordinate_t=last_t, metadata=raw_metadata)
            if history.bounds:
                path = folder/'history.npz'; history.save(path, meta)
                result['dense_history'] = dict(path=path.relative_to(ROOT).as_posix(), sha256=sha(path),
                    coordinate_bounds=[history.bounds[0][0], history.bounds[-1][1]])
            if status != 'event_reached' and (not saved or last_t > saved[-1][0]+2e-12):
                try:
                    endpoint = snapshot(system, last_t, last_state, phase)
                    json.dumps(endpoint, allow_nan=False, default=json_scalar)
                    result['accepted_endpoint_diagnostic'] = endpoint
                except (RuntimeError, ValueError, TypeError, FloatingPointError) as exc:
                    result['accepted_endpoint_diagnostic_failure'] = str(exc)
            controlled = [row['coordinate_t'] for row in result['snapshots'] if row['numerical_guard_passed']]
            endpoint = result.get('accepted_endpoint_diagnostic')
            if endpoint and endpoint['numerical_guard_passed']:
                controlled.append(endpoint['coordinate_t'])
            result.setdefault('physical_stop_coordinate_t', max(controlled) if controlled else t0)
            if 'dense_history' in result:
                result['dense_history']['unpromoted_tail_after'] = result['physical_stop_coordinate_t']
            result.update(status=status, stop_reason=reason, state_cache=cache.relative_to(ROOT).as_posix(),
                state_cache_sha256=sha(cache), actual_last_accepted_coordinate_t=last_t,
                actual_last_accepted_proper=float(last_state[system.dof+8]),
                accepted_endpoint_covered_by_dense=bool(history.bounds and history.bounds[-1][1] == last_t),
                last_numerical_guard_passed_coordinate_t=max(controlled) if controlled else None,
                steps=len(history.bounds), rhs_calls=system.calls,
                producer_source_unchanged_during_run=all(sha(ROOT/name)==digest for name, digest in pins.items()),
                history_scope='Original prefix only through proper.9, new dense history thereafter through physical_stop; any later accepted tail is preserved but not promoted.',
                event_scope='Finite event bracketing with five probes per accepted interval, not a theorem excluding missed substep double crossings. Bare-top crossing alone is not escape or a child.')
            try:
                save_record()
            except (TypeError, ValueError) as exc:
                # A minimal terminal record names already-flushed raw bytes.
                fallback = dict(schema='nsc-atlas-parent-release-export-error-v1', status='export_failure',
                    export_error=str(exc), producer_commit=commit, source_sha256=pins,
                    original_parent=data['provenance'], state_cache=result['state_cache'],
                    state_cache_sha256=result['state_cache_sha256'], dense_history=result.get('dense_history'),
                    actual_last_accepted_coordinate_t=last_t,
                    actual_last_accepted_proper=float(last_state[system.dof+8]))
                (folder/'run.json').write_text(json.dumps(fallback, allow_nan=False, default=json_scalar)+'\n',
                    encoding='utf-8', newline='\n')
                result['status'], result['export_error'] = 'export_failure', str(exc)
        print(json.dumps(dict(status=status, reason=reason, output=str(folder/'run.json'))), flush=True)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--degree', type=int, choices=(20, 24), required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--budget', type=float, default=1800.)
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--max-step', type=float, default=.01)
    parser.add_argument('--proper-endpoint', type=float, default=1.6)
    parser.add_argument('--rtol', type=float, default=3e-9)
    parser.add_argument('--atol', type=float, default=3e-11)
    args = parser.parse_args()
    if (not all(math.isfinite(v) and v > 0 for v in (args.budget, args.max_step, args.rtol, args.atol))
            or not .9 < args.proper_endpoint <= 1.6 or not 1 <= args.threads <= 8):
        raise ValueError('finite positive controls, .9<proper_endpoint<=1.6 and threads1..8 required')
    run(args)


if __name__ == '__main__':
    main()
