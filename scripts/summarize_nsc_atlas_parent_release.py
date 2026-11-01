"""Authenticate terminal same-state runs and export their compact diagnostics.

No producer, initializer, history evaluator or physical evolution is imported.
Run only after both continuation jobs are terminal; no polling or retry occurs.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TERMINAL = {'event_reached', 'budget_stop', 'numerical_stop', 'export_stop', 'incomplete'}
CHECKPOINTS = {'same_Cauchy_restart', 'proper_checkpoint', 'proper_endpoint'}
ROW_KEYS = ('event', 'coordinate_t', 'observer_proper_duration', 's', 'uS', 'nS', 'chi', 'nChi',
    'whole_grid_normal_momentum_ranges', 'ADM_H_RMS_over_H_squared', 'ADM_M_RMS_over_H_squared',
    'GH_C_RMS_over_H', 'scalar_last_two_shell_power_fraction', 'volume_last_two_shell_power_fraction',
    'scalar_acceleration_unresolved_fraction', 'volume_acceleration_unresolved_fraction',
    'observer_unit_norm_error', 'minimum_lapse', 'maximum_lapse', 'metric_valid', 'numerical_guard_passed')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1024*1024), b''):
            digest.update(part)
    return digest.hexdigest()


def strict_json(text):
    def reject(value):
        raise ValueError('nonfinite JSON constant: '+value)
    return json.loads(text, parse_constant=reject)


def json_scalar(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError('unsupported JSON type: '+type(value).__name__)


def local_path(name):
    path = (ROOT/Path(name)).resolve()
    if not path.is_relative_to(ROOT) or not path.is_file():
        raise ValueError('missing or external repository input: '+str(name))
    return path


def checked_file(name, expected):
    path = local_path(name)
    if not re.fullmatch('[0-9a-f]{64}', str(expected)) or sha(path) != expected:
        raise ValueError('raw SHA256 mismatch: '+str(name))
    return path


def potential(s):
    return 7*math.exp(-12*s)-7.5*math.exp(-14*s)+(25/32)*math.exp(-24*s)


def scalar_row(row, sstar, top):
    selected = {key: row[key] for key in ROW_KEYS if key in row}
    s, velocity = float(row['s']), float(row['uS'])
    V = potential(s)
    selected.update(bare_potential=V, mechanical_radion_energy=.5*velocity**2+V,
        bare_mechanical_headroom=potential(top)-.5*velocity**2-V,
        internal_volume_over_vacuum=math.exp(12*(s-sstar)))
    return selected


def read_terminal(path):
    path = local_path(path); record = strict_json(path.read_text(encoding='utf-8'))
    if record.get('schema') != 'nsc-atlas-parent-release-continuation-v1' or record.get('status') not in TERMINAL:
        raise ValueError('complete terminal continuation record required: '+str(path))
    if not record.get('stop_reason') or not record.get('source_bound_to_commit') or not record.get('producer_source_unchanged_during_run'):
        raise ValueError('terminal/source authentication markers missing')
    commit, pins = record['producer_commit'], record['source_sha256']
    if not re.fullmatch('[0-9a-f]{40}', commit) or not pins:
        raise ValueError('invalid producer identity')
    for name, expected in pins.items():
        pure = PurePosixPath(name)
        if pure.is_absolute() or '..' in pure.parts or ':' in name or '\\' in name:
            raise ValueError('invalid archived source path')
        archived = subprocess.run(['git', 'show', commit+':'+name], cwd=ROOT, capture_output=True)
        if archived.returncode or hashlib.sha256(archived.stdout).hexdigest() != expected:
            raise ValueError('archived producer source mismatch: '+name)
    states = checked_file(record['state_cache'], record['state_cache_sha256'])
    with np.load(states, allow_pickle=False) as raw:
        times, saved = raw['times'], raw['states']
        metadata = strict_json(str(raw['metadata']))
        if saved.ndim != 2 or not len(saved) or times.shape != (len(saved),) or not np.all(np.isfinite(saved)):
            raise ValueError('finite first saved Cauchy state required')
        if hashlib.sha256(saved[0].tobytes()).hexdigest() != record['start_full_state_sha256']:
            raise ValueError('restart first full-state SHA mismatch')
        if abs(float(times[0])-record['start_coordinate_t']) > 2e-12 or abs(float(saved[0, -2])-record['start_proper']) > 2e-10:
            raise ValueError('restart coordinate/proper event mismatch')
    if metadata['source_sha256'] != pins or metadata['producer_commit'] != commit or metadata['degree'] != record['degree']:
        raise ValueError('state metadata producer mismatch')
    dense = record.get('dense_history')
    if dense:
        history = checked_file(dense['path'], dense['sha256'])
        with np.load(history, allow_pickle=False) as raw:
            meta = strict_json(str(raw['metadata'])); bounds = raw['bounds']
        if meta['source_sha256'] != pins or meta['producer_commit'] != commit or not len(bounds):
            raise ValueError('history metadata producer mismatch')
        if list(map(float, (bounds[0, 0], bounds[-1, 1]))) != dense['coordinate_bounds']:
            raise ValueError('history coordinate-bound mismatch')
    rows = record['snapshots']
    if len(rows) != len(times) or any(abs(float(r['coordinate_t'])-t) > 2e-12 for r, t in zip(rows, times)):
        raise ValueError('saved-state/snapshot events differ')
    sstar = float(record['model_inputs']['sstar']); top = float(record['thresholds']['terminal_bare_top_s'])
    keys = ('degree', 'status', 'stop_reason', 'producer_commit', 'source_sha256', 'start_coordinate_t',
        'start_proper', 'requested_proper_endpoint', 'physical_stop_coordinate_t', 'last_numerical_guard_passed_coordinate_t',
        'actual_last_accepted_coordinate_t', 'actual_last_accepted_proper', 'accepted_endpoint_covered_by_dense',
        'budget_seconds', 'elapsed_seconds', 'threads', 'temporal_tolerances', 'steps', 'rhs_calls', 'events',
        'auxiliary_phase', 'history_scope', 'event_scope', 'physical_history_prefix', 'original_parent')
    summary = {key: record[key] for key in keys if key in record}
    summary.update(raw_record=dict(path=path.relative_to(ROOT).as_posix(), sha256=sha(path)),
        state_cache=dict(path=record['state_cache'], sha256=record['state_cache_sha256'], first_state_sha256=record['start_full_state_sha256']),
        dense_history=dense, bare_top=dict(s=top, potential=potential(top)),
        scalar_rows=[scalar_row(row, sstar, top) for row in rows],
        accepted_endpoint_diagnostic=(scalar_row(record['accepted_endpoint_diagnostic'], sstar, top)
            if 'accepted_endpoint_diagnostic' in record else None),
        accepted_endpoint_diagnostic_failure=record.get('accepted_endpoint_diagnostic_failure'))
    return summary


def common_checkpoints(first, second):
    out = []
    for left in first['scalar_rows']:
        if left['event'] not in CHECKPOINTS or not left['numerical_guard_passed']:
            continue
        matches = [r for r in second['scalar_rows'] if r['event'] in CHECKPOINTS and r['numerical_guard_passed']
                   and abs(r['observer_proper_duration']-left['observer_proper_duration']) <= 2e-10]
        if len(matches) > 1:
            raise ValueError('ambiguous common proper checkpoint')
        if matches:
            right = matches[0]
            keys = ('s', 'uS', 'nS', 'nChi', 'ADM_H_RMS_over_H_squared', 'ADM_M_RMS_over_H_squared', 'GH_C_RMS_over_H')
            out.append(dict(proper=left['observer_proper_duration'], coordinate_times=[left['coordinate_t'], right['coordinate_t']],
                proper_mismatch=abs(left['observer_proper_duration']-right['observer_proper_duration']),
                absolute_differences={key: abs(right[key]-left[key]) for key in keys}))
    return out


def fixed_surface_check(runs):
    """Check only the predeclared measurement boundary at common proper1.2.

    This is a readout of authenticated states, not a new horizon search.
    """
    if str(ROOT) not in sys.path:
        sys.path.insert(0,str(ROOT))
    from scripts.parent_release_region import source_region
    from recursive_horizons.nsc_atlas_parent_gh import ParentGH
    from recursive_horizons.nsc_atlas_trapped_parent import GHSurfaceGeometry
    from recursive_horizons.nsc_atlas_trapped_surfaces import round_surface,evaluate_surface
    from threadpoolctl import threadpool_limits
    rows=[]
    with threadpool_limits(limits=4):
        for run in runs:
            with np.load(ROOT/run['state_cache']['path'],allow_pickle=False) as raw:
                meta=strict_json(str(raw['metadata'])); times=raw['times']; states=raw['states']
            candidates=[i for i,row in enumerate(run['scalar_rows'])
                if row['numerical_guard_passed'] and abs(row['observer_proper_duration']-1.2)<2e-10]
            if len(candidates)!=1:
                raise ValueError('One guard-passing proper1.2 state is required')
            index=candidates[0]
            system=ParentGH(run['degree'],meta['epsilon'],gauge=meta['gauge'],
                            backend=meta['backend'],quadrature_order=meta['quadrature_order'])
            definition=source_region(system.parent)
            surface=round_surface(definition['sampler_arguments']['anchor'],definition['sampler_arguments']['radius'])
            provider=GHSurfaceGeometry(system,states[index],coordinate_t=float(times[index]),
                provenance=dict(source_binding='SHA-authenticated terminal continuation state',
                    state_cache_sha256=run['state_cache']['sha256']))
            low=evaluate_surface(surface,provider,frame_scale=system.r,quadrature_order=(32,64))
            high=evaluate_surface(surface,provider,frame_scale=system.r,quadrature_order=(40,80))
            rows.append(dict(degree=run['degree'],proper=run['scalar_rows'][index]['observer_proper_duration'],
                coordinate_t=float(times[index]),state_cache=run['state_cache'],
                state_sha256=hashlib.sha256(states[index].tobytes()).hexdigest(),
                region_definition=definition,nominal=low.summary(),higher_angular_order=high.summary(),
                area_difference=abs(low.area-high.area),
                scope='One fixed source-derived boundary only; no optimization, alternative surface search or all-surface nonexistence theorem.'))
    paths=['scripts/summarize_nsc_atlas_parent_release.py','scripts/parent_release_region.py',
        'src/recursive_horizons/nsc_atlas_trapped_surfaces.py','src/recursive_horizons/nsc_atlas_trapped_parent.py']
    return dict(schema='nsc-atlas-parent-release-surface-v1',runs=rows,
        source_sha256={name:sha(ROOT/name) for name in paths},
        physical_child=dict(formation_proved=False,volume=None,energy=None,w=None))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--records', nargs=2, type=Path, default=[ROOT/f'.build/parent-release-continuation-L{n}-v3/run.json' for n in (20, 24)])
    parser.add_argument('--output', type=Path, default=ROOT/'results/summaries/nsc-atlas-parent-release-continuation-v1.json')
    parser.add_argument('--surface-output',type=Path,help='Optional separate fixed-boundary check; does not rewrite the terminal summary')
    args = parser.parse_args(); output = (ROOT/args.output).resolve()
    if args.surface_output:
        output=(ROOT/args.surface_output).resolve()
    if output.exists() or not output.is_relative_to(ROOT/'results/summaries'):
        raise ValueError('fresh summary output under results/summaries required')
    runs = [read_terminal(path) for path in args.records]
    if [r['degree'] for r in runs] != [20, 24]:
        raise ValueError('explicit degree20 then degree24 terminal records required')
    if args.surface_output:
        result=fixed_surface_check(runs)
        output.parent.mkdir(parents=True,exist_ok=True)
        with output.open('x',encoding='utf-8',newline='\n') as stream:
            stream.write(json.dumps(result,indent=2,allow_nan=False,default=json_scalar)+'\n')
        print(json.dumps(dict(output=str(output),fixed_boundary_only=True)))
        return
    result = dict(schema='nsc-atlas-parent-release-continuation-summary-v1', runs=runs,
        shared_proper_checkpoints=common_checkpoints(*runs),
        comparison_scope='Empirical same-proper checkpoint differences only; no off-clock event comparison or error enclosure.',
        interpretation='Reference fold is not a physical transition; mechanical barrier crossing is not escape. Accepted tails stay unpromoted.',
        physical_child=dict(formation_proved=False, volume=None, energy=None, w=None),
        exporter=dict(path=Path(__file__).relative_to(ROOT).as_posix(), sha256=sha(Path(__file__))),
        scope='Terminal authenticated diagnostic export only; no new fields, states, PDE, causal handoff or child claimed.')
    payload = json.dumps(result, indent=2, allow_nan=False, default=json_scalar)+'\n'
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x', encoding='utf-8', newline='\n') as stream:
        stream.write(payload)
    print(json.dumps(dict(output=output.relative_to(ROOT).as_posix(), statuses=[r['status'] for r in runs])))


if __name__ == '__main__':
    main()
