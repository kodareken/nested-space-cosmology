#!/usr/bin/env python3
"""Two retained terminal cap balances; no grid scan, fit or time evolution."""
import argparse
import json
import math
from pathlib import Path
import time
import uuid

import numpy as np

import diagnose_nsc_atlas_cap_mass_balance as previous
from analyze_nsc_atlas_cap_consistency import _array, _numbers, _scalar
from analyze_nsc_atlas_cap_mass_direction import _saved_terminal
from analyze_nsc_atlas_cap_mass_balance import same_values
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_background_cap import background_cap_metric, background_cap_rhs
from recursive_horizons.nsc_atlas_cap_mass_tangent import mass_tangent
from recursive_horizons.nsc_atlas_cap_transport_balance import cap_transport_balance
from validation_source_successors import verify_saved_bindings
from summarize_result import relative_file, create_summary

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = 'nsc-atlas-cap-transport-balance-audit-v1'
SCOPE = 'fixed terminal cap bulk/interface/radial-error attribution; not corrected evolution or closure'
PRIOR = 'results/summaries/nsc-atlas-cap-mass-balance-v1.json'
PRIOR_SHA256 = 'd3d7be8d7b2a124b6e6facbfc8ddfb86ff1e66eceb1cd16553fa7c243c8109e0'
LIMITS = dict(creation=dict(evaluated_points=48_000_000, cpu_seconds=300.),
              check=dict(evaluated_points=24_000_000, cpu_seconds=150.))
SETTINGS = [dict(rtol=2e-10, atol=2e-12), dict(rtol=2e-12, atol=2e-14)]
GATES = dict(reference_mass=1e-8, wall_density_effect=1e-7, reference_rate=2e-7,
             jump_attenuation=1e-7, balance=1e-6, channel_sensitivity=5e-7)
CHANNELS = ('radial_map_error_rate', 'wall_density_conversion', 'wall_flux_conversion',
            'field_f', 'field_Pi', 'field_Phi', 'interfaces')
EXPECTED_CALLS = dict(metric_calls=2, rhs_calls=2, tangent_calls=2, balance_calls=4,
    profile_constructions=1, profile_fits=0, source_solves=0, evolution_calls=0, time_steps=0)
CHARGES = ('input_and_profile', 'metrics', 'rhs', 'tangents', 'balances')
OWN_SOURCES = (*previous.OWN_SOURCES,
    *('src/recursive_horizons/nsc_atlas_cap_'+name+'.py' for name in
      ('transport_fields', 'radial_reference', 'radial_comparison', 'transport_reference', 'transport_balance')),
    *('tests/test_nsc_atlas_cap_'+name+'.py' for name in
      ('transport_fields', 'radial_reference', 'radial_comparison', 'transport_reference', 'transport_balance')),
    'scripts/diagnose_nsc_atlas_cap_transport_balance.py',
    'tests/test_diagnose_nsc_atlas_cap_transport_balance.py',
    'scripts/summarize_result.py')
Budget, encoded, digest, plain = previous.Budget, previous.encoded, previous.digest, previous.plain
runtime_versions = previous.runtime_versions


def shared_budget(value=None, *, checking=False):
    limit = LIMITS['check' if checking else 'creation']
    if value is None:
        value = Budget(maximum_evaluations=limit['evaluated_points'], cpu_seconds=limit['cpu_seconds'])
    if (not isinstance(value, Budget) or type(value.maximum_evaluations) is not int
            or not 0 < value.maximum_evaluations <= limit['evaluated_points']
            or type(value.evaluations) is not int or value.evaluations < 0
            or type(value.cpu_seconds) not in (int, float) or not math.isfinite(value.cpu_seconds)
            or not 0 < value.cpu_seconds <= limit['cpu_seconds']
            or type(value.started) not in (int, float) or not math.isfinite(value.started)
            or not 0 <= value.started <= time.process_time()):
        raise ValueError('invalid shared terminal transport allowance')
    value.check()
    return value


def bindings(root):
    return {name: digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}


def authenticate(root, budget, sources):
    inputs, identity = previous.authenticate(root, budget, sources)
    path = Path(root)/PRIOR
    if path.is_symlink():
        raise ValueError('mass ledger summary must be a regular file')
    raw = path.read_bytes()
    if digest(raw) != PRIOR_SHA256:
        raise ValueError('pinned mass ledger summary mismatch')
    summary = json.loads(raw)
    verify_saved_bindings({k: v['sha256'] for k, v in summary['producer_bindings'].items()},
                          root, expected_paths=previous.OWN_SOURCES)
    budget.check()
    return (*inputs[:3], summary), dict(inherited=identity,
        mass_ledger_summary=dict(path=PRIOR, sha256=PRIOR_SHA256, bytes=len(raw)),
        previous_analysis_replays=0)


def run_summary(result, direction, old):
    """Gate independent routes before exposing an interpretation of channels."""
    maps, ref = result['maps'], result['reference']
    native = [result['native_balance'], result['native_wall']['rho'], result['native_wall']['J']]
    retained = [old['B_mass'], old['radial_density'], old['flux_left']]
    try:
        same_values(native, retained, 'saved native balance/density/flux')
    except ValueError as error:
        error.accuracy_failure = dict(native=native, retained=retained,
            defects=[a-b for a,b in zip(native,retained)])
        raise
    if set(result['channels']) != set(CHANNELS):
        raise ValueError('transport channel inventory differs')
    defects = dict(reference_mass=result['reference_mass_defect'],
        wall_density_effect=direction['R']*ref['wall_fields']['a']*result['reference_mass_defect'],
        reference_rate=result['reference_rate_defect'], jump_attenuation=ref['jump_sum_defect'],
        balance=result['balance_defect'])
    failures = {key: dict(value=float(value), limit=GATES[key]) for key, value in defects.items()
                if not math.isfinite(value) or abs(value) > GATES[key]}
    if len(maps['telescope_defects'])!=2 or len(maps['telescope_roundoff_bounds'])!=2:
        raise ValueError('mass/rate telescope inventory differs')
    for label, value, bound in zip(('mass_telescope', 'rate_telescope'),
                                 maps['telescope_defects'], maps['telescope_roundoff_bounds']):
        if not math.isfinite(value) or not math.isfinite(bound) or bound < 0 or abs(value) > bound:
            failures[label] = dict(value=float(value), limit=float(bound))
    if failures:
        error = ValueError('terminal transport accuracy gate failed')
        error.accuracy_failure = plain(dict(failures=failures, defects=defects, channels=result['channels']))
        raise error
    values = dict(channels=result['channels'], total=result['total'], native_balance=result['native_balance'],
        actual_native_balance=result['actual_native_balance'], defects=defects,
        radial_mass_error=maps['error'], native_mass=maps['native_mass'], reference_mass=maps['reference_mass'],
        reference_rate=maps['reference_rate'], reference_lapse=ref['delta'],
        telescope_defects=maps['telescope_defects'], telescope_roundoff_bounds=maps['telescope_roundoff_bounds'],
        native_tangent_defect=maps['native_tangent_defect'], native_tangent_allowance=maps['native_tangent_allowance'],
        minimum_sampled_A=ref['minimum_sampled_A'], sampled_A_positive=ref['sampled_A_positive'],
        interface_count=len(ref['jumps']), piece_count=len(ref['stages']),
        uniform_error_enclosure=False, accuracy_gates_passed=True)
    encoded(plain(values))
    return values


def case_summary(label, alpha, direction, runs):
    changes = {key: runs[1]['channels'][key]-runs[0]['channels'][key] for key in CHANNELS}
    if any(not math.isfinite(v) or abs(v) > GATES['channel_sensitivity'] for v in changes.values()):
        error = ValueError('terminal transport tolerance sensitivity failed')
        error.accuracy_failure = dict(channel_changes=changes, limit=GATES['channel_sensitivity'])
        raise error
    return dict(amplitude=label, alpha=alpha, direction=direction, runs=runs,
                fine_minus_coarse=changes, tolerance_sensitivity_passed=True)


def analyze(current, cap, cache, prior, budget):
    shared_budget(budget)
    start = budget.evaluations
    charges = dict.fromkeys(CHARGES, 0)
    counts = dict.fromkeys(('metric_calls', 'rhs_calls', 'tangent_calls', 'balance_calls'), 0)
    phase = dict(name='fixed input declarations')
    try:
        if [(v['amplitude'], v['cells'], v['steps']) for v in current['cases']] != [('zero',64,512), ('finite',64,512)]:
            raise ValueError('original zero/finite 64-cell cases required')
        old = cap['analysis']; ledger = prior['data']
        if (old['schema'] != 'nsc-atlas-cap-consistency-analysis-v1'
                or ledger['schema'] != 'nsc-atlas-cap-mass-balance-analysis-v1'
                or ledger['cut'] != 128 or ledger['quadrature_orders'] != [4,8]
                or any([c['amplitude'] for c in d['cases']] != ['zero','finite'] for d in (old,ledger))):
            raise ValueError('retained comparison inventory differs')
        x = _array(cache['x_nodes'], (1025,), 'cache labels', budget)
        if np.any(np.diff(x)<=0) or not np.allclose(np.diff(x), (x[-1]-x[0])/1024, rtol=0, atol=8*np.finfo(float).eps):
            raise ValueError('uniform original cache required')
        profile = CharacteristicProfile(*[_array(current['profile'][k], (65,), 'profile.'+k, budget)
                                        for k in ('labels','values','first','second')])
        if not profile.quiet_extension:
            raise ValueError('original quiet extension required')
        charges['input_and_profile'] = budget.evaluations-start
        details, summaries = [], []
        for case, cap_case, prior_case in zip(current['cases'], old['cases'], ledger['cases']):
            before = budget.evaluations
            label, alpha = case['amplitude'], _scalar(case['alpha'], 'alpha')
            phase = dict(name='terminal input', amplitude=label, cut=128, order=8)
            if ((label=='zero' and alpha!=0) or (label=='finite' and alpha<=0)
                    or cap_case['alpha']!=alpha or prior_case['alpha']!=alpha or len(case['snapshots'])!=129):
                raise ValueError('terminal amplitude/cut declarations differ')
            sample = case['snapshots'][128]
            saved, target = _saved_terminal(cap_case, 8)
            if sample['x']!=x[-1] or saved['x']!=sample['x']:
                raise ValueError('terminal label differs')
            physical = _numbers(sample['physical']['wall'], ('r','T','L','s','c'), 'wall')
            for k in ('r','T','L','s'):
                same_values(physical[k], cache['gauge_rows'][k][-1]+sample['state']['wall'][k], 'wall.'+k)
            state = {k: _array(sample['state']['cap'][k], (64,), 'cap.'+k, budget).copy() for k in ('Ip','Iq')}
            charges['input_and_profile'] += budget.evaluations-before
            phase['name'] = 'base metric'; before = budget.evaluations; counts['metric_calls'] += 1
            metric = background_cap_metric(np.linspace(0,1,65), physical['r'], state, T=physical['T'],
                profile=profile, alpha=alpha, quadrature_order=8, budget=budget)
            charges['metrics'] += budget.evaluations-before
            for key in ('n_faces','delta_faces'):
                same_values(metric[key], saved[key], key)
            for key in ('A','delta','f','chi','c'):
                same_values(metric['wall'][key], saved['wall'][key], 'wall.'+key)
            same_values(metric['wall']['c'], physical['c'], 'physical scalar')
            h = _numbers(target['hybrid'], ('r_s','T_s','Vdot_plus','w_plus'), 'hybrid')
            if h['T_s']<=0 or h['Vdot_plus']<=0:
                raise ValueError('positive stored clocks required')
            direction = dict(R_x=h['r_s']/h['Vdot_plus'], T_x=h['T_s']/h['Vdot_plus'], w_plus=h['w_plus'])
            prior_order = prior_case['orders']['8']
            same_values(list(direction.values()), [prior_order['direction'][k] for k in direction], 'direction')
            phase['name'] = 'base departure RHS'; before = budget.evaluations; counts['rhs_calls'] += 1
            rhs = background_cap_rhs(metric, T_x=direction['T_x'], radius_x=direction['R_x'],
                                     incoming_w_plus=direction['w_plus'], budget=budget)
            charges['rhs'] += budget.evaluations-before
            full_direction = dict(R=direction['R_x'], T=direction['T_x'], **rhs['rates'])
            phase['name'] = 'mass tangent'; before = budget.evaluations; counts['tangent_calls'] += 1
            tangent = mass_tangent(metric, full_direction, budget=budget)
            charges['tangents'] += budget.evaluations-before
            for key in ('N','N_R','N_T','forward','reverse'):
                try:
                    same_values(tangent[key], prior_order['tangent'][key], 'saved tangent.'+key)
                except ValueError as error:
                    error.accuracy_failure = dict(tangent_key=key, actual=tangent[key],
                        retained=prior_order['tangent'][key], defect=tangent[key]-prior_order['tangent'][key])
                    raise
            runs, compact = [], []
            for settings in SETTINGS:
                phase.update(name='full reference balance', settings=settings)
                before = budget.evaluations; counts['balance_calls'] += 1
                result = cap_transport_balance(metric, full_direction, tangent, budget=budget, **settings)
                charges['balances'] += budget.evaluations-before
                phase['name'] = 'saved native ledger and accuracy gates'
                compact.append(dict(settings=settings, **run_summary(result, full_direction, prior_order)))
                runs.append(dict(settings=settings, result=result))
            phase['name'] = 'tolerance sensitivity'
            summaries.append(case_summary(label, alpha, direction, compact))
            details.append(dict(amplitude=label, runs=runs))
        result = plain(dict(summary=dict(cut=128, cells=64, native_order=8, gates=GATES, cases=summaries,
            scope=SCOPE, uniform_error_enclosure=False), details=details,
            execution=dict(EXPECTED_CALLS, charges=charges), evaluated_points=budget.evaluations-start))
        if counts != {k: EXPECTED_CALLS[k] for k in counts}:
            raise ValueError('terminal transport call inventory mismatch')
        check_accounting(result)
        budget.check()
        return result
    except Exception as error:
        error.transport_failure = plain(dict(phase=phase, counts=counts,
            evaluated_points=budget.evaluations-start, accuracy=getattr(error,'accuracy_failure',None)))
        raise


def check_accounting(result):
    execution = result['execution']; cost = result['evaluated_points']; charges = execution['charges']
    if (type(cost) is not int or cost<=0 or cost>LIMITS['check']['evaluated_points']
            or set(execution)!=set(EXPECTED_CALLS)|{'charges'}
            or any(type(execution[k]) is not int or execution[k]!=v for k,v in EXPECTED_CALLS.items())
            or set(charges)!=set(CHARGES) or any(type(v) is not int or v<0 for v in charges.values())
            or sum(charges.values())!=cost):
        raise ValueError('terminal transport accounting mismatch')


def replay_equal(actual, expected):
    """Numerical payload replay: exact structure, bounded arithmetic, adaptive counts recorded separately."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual)!=set(expected):
            return False
        for key, value in expected.items():
            if key in ('settings','gates','requested_rtol','requested_atol'):
                if actual[key]!=value:
                    return False
            elif key in ('evaluations','coefficient_calls','reference_coefficient_calls','rhs_calls'):
                if type(actual[key]) is not int or not 0<=actual[key]<=LIMITS['check']['evaluated_points']:
                    return False
            elif not replay_equal(actual[key], value):
                return False
        return True
    if isinstance(expected, list):
        return isinstance(actual,list) and len(actual)==len(expected) and all(replay_equal(a,b) for a,b in zip(actual,expected))
    if type(expected) is float:
        return (type(actual) is float and math.isfinite(actual) and math.isfinite(expected)
                and math.isclose(actual,expected,rel_tol=2e-12,abs_tol=1e-8))
    return type(actual) is type(expected) and actual==expected


def validate_saved_gates(analysis, prior):
    """A replay allowance never turns a saved failing gate into a passing one."""
    summary = analysis['summary']
    if (summary['gates']!=GATES or summary['cut']!=128 or summary['cells']!=64
            or summary['native_order']!=8 or summary['scope']!=SCOPE
            or summary['uniform_error_enclosure'] is not False
            or [c['amplitude'] for c in analysis['details']]!=['zero','finite']
            or [c['amplitude'] for c in summary['cases']]!=['zero','finite']):
        raise ValueError('saved terminal summary inventory differs')
    for detailed, compact, retained in zip(analysis['details'], summary['cases'], prior['data']['cases']):
        direction = compact['direction']
        same_values([direction[k] for k in ('R_x','T_x','w_plus')],
                    [retained['orders']['8']['direction'][k] for k in ('R_x','T_x','w_plus')], 'saved direction')
        if [r['settings'] for r in detailed['runs']]!=SETTINGS or compact['alpha']!=retained['alpha']:
            raise ValueError('saved terminal tolerance/amplitude declarations differ')
        runs = [dict(settings=r['settings'], **run_summary(r['result'],
                    dict(R=direction['R_x'],T=direction['T_x']),retained['orders']['8'])) for r in detailed['runs']]
        expected = plain(case_summary(compact['amplitude'],compact['alpha'],direction,runs))
        if compact!=expected:
            raise ValueError('saved scalar summary differs from its detailed balances')


def validate(saved, inputs, identity, sources, budget):
    if (saved.get('schema')!=SCHEMA or saved.get('scope')!=SCOPE or saved.get('limits')!=LIMITS
            or saved.get('settings')!=SETTINGS or saved.get('input')!=identity
            or saved.get('runtime')!=runtime_versions() or saved.get('audit_source_bindings')!=sources):
        raise ValueError('terminal transport declaration/input/producer mismatch')
    check_accounting(saved['analysis'])
    validate_saved_gates(saved['analysis'],inputs[3])
    data = saved['execution']
    if (set(data)!={'analysis_evaluations','creation_validation_evaluations','evaluated_points','initial_evaluations',
                   'cpu_seconds','hash_authentication_passes','previous_analysis_replays','time_steps'}
            or any(type(v) is not int or v<0 for k,v in data.items() if k!='cpu_seconds')
            or type(data['cpu_seconds']) not in (float,int) or not 0<=data['cpu_seconds']<=LIMITS['creation']['cpu_seconds']
            or data['analysis_evaluations']!=saved['analysis']['evaluated_points']
            or not 0<data['creation_validation_evaluations']<=LIMITS['check']['evaluated_points']
            or data['evaluated_points']!=data['analysis_evaluations']+data['creation_validation_evaluations']
            or data['evaluated_points']+data['initial_evaluations']>LIMITS['creation']['evaluated_points']
            or (data['hash_authentication_passes'],data['previous_analysis_replays'],data['time_steps'])!=(1,0,0)):
        raise ValueError('invalid terminal transport execution receipt')
    expected = analyze(*inputs, budget)
    # Adaptive work need not be bitwise identical on another platform. Inputs,
    # producer bytes, fixed call counts and each analysis's accounting remain exact.
    if not replay_equal({k:saved['analysis'][k] for k in ('summary','details')},
                        {k:expected[k] for k in ('summary','details')}):
        raise ValueError('terminal transport numerical reconstruction mismatch')
    return expected['evaluated_points']


def operation(path, *, checking=False, root=ROOT, budget=None):
    budget = shared_budget(budget, checking=checking)
    path, relative = relative_file(root,path,must_exist=checking)
    if Path(relative).parts[:2]!=('.build','results') or path.suffix!='.json':
        raise ValueError('full terminal transport data belong under ignored .build/results')
    if not checking and path.exists():
        raise FileExistsError('terminal transport output is creation-only')
    initial = budget.evaluations; phase = 'producer bindings'
    try:
        sources, runtime = bindings(root), runtime_versions()
        phase = 'input authentication'; inputs, identity = authenticate(root,budget,sources)
        if checking:
            phase = 'check reconstruction'; payload = json.loads(path.read_bytes())
            cost = validate(payload,inputs,identity,sources,budget)
        else:
            phase = 'terminal transport analysis'; result = analyze(*inputs,budget)
            cost = result['evaluated_points']
            payload = dict(schema=SCHEMA, scope=SCOPE, limits=LIMITS, settings=SETTINGS, input=identity,
                runtime=runtime, audit_source_bindings=sources, analysis=result,
                execution=dict(analysis_evaluations=cost, creation_validation_evaluations=cost,
                    evaluated_points=2*cost, initial_evaluations=initial, cpu_seconds=time.process_time()-budget.started,
                    hash_authentication_passes=1, previous_analysis_replays=0, time_steps=0))
            phase = 'creation reconstruction'; repeat = validate(payload,inputs,identity,sources,budget)
            payload['execution'].update(creation_validation_evaluations=repeat, evaluated_points=cost+repeat,
                                        cpu_seconds=time.process_time()-budget.started)
            cost += repeat
        if bindings(root)!=sources or runtime_versions()!=runtime:
            raise ValueError('terminal transport producer/runtime changed during calculation')
        if budget.evaluations-initial!=cost:
            raise ValueError('unexpected charges outside terminal analyses')
        budget.check()
        if not checking:
            raw = encoded(payload); budget.check()
    except Exception as error:
        target = Path(root)/'.build/cap-transport-balance-failures'/(uuid.uuid4().hex+'.json')
        target.parent.mkdir(parents=True,exist_ok=True)
        with target.open('xb') as stream:
            stream.write(encoded(failure_values(dict(schema=SCHEMA+'-failed-attempt', phase=phase, exception=type(error).__name__,
                reason=str(error), retries=0, evaluated_points=budget.evaluations-initial,
                cpu_seconds=time.process_time()-budget.started, analysis_failure=getattr(error,'transport_failure',None)))))
        error.failure_receipt = str(target)
        raise
    if checking:
        return dict(checked=str(path), evaluated_points=cost, cpu_seconds=time.process_time()-budget.started,
                    hash_authentication_passes=1, time_steps=0)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('xb') as stream:
        stream.write(raw)
    return payload


def failure_values(value):
    if isinstance(value,dict):
        return {k:failure_values(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):
        return [failure_values(v) for v in value]
    if isinstance(value,float) and not math.isfinite(value):
        return str(value)
    return value


def summarize(path, output, *, root=ROOT):
    _, relative = relative_file(root,path)
    if Path(relative).parts[:2]!=('.build','results'):
        raise ValueError('terminal raw data must be under .build/results')
    return create_summary(path,output,root=root,selection='/analysis/summary',
        argv=['python','scripts/run_python.py','scripts/diagnose_nsc_atlas_cap_transport_balance.py',
              '--output',relative],producer_paths=OWN_SOURCES,
        input_paths=(PRIOR,previous.PRIOR,previous.previous.CAP_INPUT))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--check',action='store_true')
    parser.add_argument('--summary',type=Path,help='create a compact tracked summary after successful creation/check')
    args = parser.parse_args()
    try:
        result = operation(args.output,checking=args.check)
        if args.summary:
            summarize(args.output,args.summary)
    except Exception as error:
        if hasattr(error,'failure_receipt'):
            print(json.dumps(dict(failed=str(error),diagnostic=error.failure_receipt)))
        raise
    print(json.dumps(result if args.check else dict(created=str(args.output),execution=result['execution']),sort_keys=True))


if __name__=='__main__':
    main()
