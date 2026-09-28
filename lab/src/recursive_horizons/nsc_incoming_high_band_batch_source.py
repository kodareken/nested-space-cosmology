"""Resumable source preparation for the28 remaining high-band channel groups."""
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from .nsc_incoming_window_refinement import authenticated_window, window_kernel, window_quadrature
from .nsc_incoming_source_assembly import read_panels
from .nsc_incoming_source_tail import _riccati_at_one, _mp_bloch
from .nsc_incoming_middle_order_correction import stable_delta_bloch, source16_control
from .nsc_incoming_state_moments import incoming_group_factor
from .nsc_incoming_middle_bound import thermal_middle_difference_bound


EXCLUDED = (12, 14, 22, 32)
GROUPS = (6,)+tuple(g for g in range(1, 33) if g not in (*EXCLUDED, 6))
TOLERANCES = {'contraction_replay': 3e-22, 'prefix16': 3e-40,
    'original_measure': 3e-11, 'quadrature_measure': 3e-11,
    'source16_convention': 3e-11, 'horner_vs_frozen_kernel': 3e-60,
    'original_to48': 3e-20, 'quadrature48_to64': 3e-20, 'precision50_to70': 3e-35}


def horner(values, x):
    result = 0
    for coefficient in reversed(values): result = coefficient+x*result
    return result


@lru_cache(maxsize=1)
def horner_identity():
    x = sp.Symbol('x'); c = sp.symbols('c1:25')
    differences = [x*horner(c[:16], x)-sum(c[j]*x**(j+1) for j in range(16)),
        x**17*horner(c[16:], x)-sum(c[j]*x**(j+1) for j in range(16, 24)),
        x*horner(c, x)-sum(c[j]*x**(j+1) for j in range(24))]
    residuals = [str(sp.expand(v)) for v in differences]
    if residuals != ['0']*3: raise ArithmeticError('Horner polynomial identity failed')
    return residuals


@lru_cache(maxsize=2)
def geometry_constants(dps):
    with mp.workdps(dps): return mp.sqrt(3*mp.pi/2-4), mp.sqrt(2), 3*mp.pi/4


@lru_cache(maxsize=16)
def cached_grid(dps, left, right, points):
    with mp.workdps(dps):
        E, w = window_quadrature(left, right, points)
        if any(v <= 0 for v in w): raise ArithmeticError('positive source quadrature required')
        return tuple(E), tuple(w)


@lru_cache(maxsize=128)
def channel_coefficients(dps, mass_value, angular_value):
    with mp.workdps(dps):
        mass, angular = mp.mpf(mass_value), mp.mpf(angular_value)
        c = _riccati_at_one(mass, angular, order=24); c16 = _riccati_at_one(mass, angular, order=16)
        return tuple(c), float(max(abs(a-b) for a, b in zip(c[:16], c16)))


def horner_kernel(energy, mass, angular, coefficients, kind, *, dps):
    E = mp.mpf(energy)
    if len(coefficients) != 24 or E <= 0 or kind not in ('low_vacuum', 'middle_correction'):
        raise ValueError('same explicit order24 positive-energy window required')
    a, r, q = geometry_constants(dps); x = 1/(2*E)
    if kind == 'middle_correction':
        S = x*horner(coefficients[:16], x); delta_S = x**17*horner(coefficients[16:], x)
        delta = stable_delta_bloch(S, delta_S, a)
    else:
        S = x*(angular/r+1j*mass)+x*x*horner(coefficients[1:], x)
        u = a*a*abs(S)**2; y = a*a*(mass*mass+angular*angular/(r*r))/(E*E); root = mp.sqrt(1+y)
        delta = [2*a*mp.im(S)/(1+u)-mass*a/(E*root), -2*a*mp.re(S)/(1+u)+angular*a/(r*E*root),
                 y/(root*(1+root))-2*u/(1+u)]
        terms = _mp_bloch()(q, mass, angular, E)
        for i in range(3): delta[i] -= sum((terms[j][i, 0] for j in range(1, 5)), mp.mpf(0))
    parallel = -E/a*delta[2]
    return mp.matrix([-mass*delta[0]+angular/r*delta[1]+parallel, parallel, mp.mpf(0), angular/(2*r)*delta[1]])


def _group13_middle(root):
    """Special archived group13 interval; never add overlapping LOW[32,40]."""
    root = Path(root); panels, meta = read_panels(root)
    record = json.loads((root/'results/development/nsc-pg-group13-covariance.json').read_text())
    for kind in ('source_hashes', 'input_hashes'):
        for p, expected in record[kind].items():
            if sha256((root/p).read_bytes()).hexdigest() != expected: raise ValueError('group13 source owner/input changed')
    spec = record['payload']
    if sha256((root/spec['path']).read_bytes()).hexdigest() != spec['sha256']: raise ValueError('group13 artifact changed')
    finite = json.loads((root/'results/development/nsc-incoming-spectral-source.json').read_text())
    inventory = json.loads((root/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())
    channel = inventory['channels'][13]; factor = incoming_group_factor(channel, (1,))
    config = json.loads((root/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    if record['energy_account']['fast_mode_window'] != [16., 160.] or record['energy_account']['low'] != [0., 16.]:
        raise ValueError('group13 early middle and nonoverlapping LOW coverage changed')
    with np.load(root/spec['path'], allow_pickle=False) as a:
        E, w = a['fast_energy_24'], a['fast_weight_24']; original_meta = json.loads(a['metadata_json'].tobytes())
    name = 'group13/mid'
    if original_meta['channel'] != channel or factor != finite['groups']['13']['factor']:
        raise ValueError('group13 labels/factor changed')
    if not np.array_equal(E, panels[name+'/energies']) or not np.array_equal(w, panels[name+'/weights']) or np.any(w <= 0):
        raise ValueError('group13 source differs from original positive middle quadrature')
    names = [p['panel'] for p in finite['groups']['13']['signed_inventory']['1']['pieces'] if '/mid' in p['panel']]
    if names != [name]: raise ValueError('actual selected group13 middle panel changed')
    selected = {k: panels[name+'/'+k].copy() for k in ('energies', 'weights', 'kernels', 'arithmetic_indicator', 'vacuum_kernels', 'thermal_kernels')}
    return {'group': 13, 'kind': 'middle_correction', 'left': 16., 'right': 160., 'channel': channel,
        'config': config, 'signs': (1,), 'factor': factor, 'selected': {'1': selected},
        'descriptors': [{'sign': 1, 'panel': name, 'rows': len(E), 'interval': [16., 160.], 'method': meta['panels'][name]['method']}],
        'measure_residual': abs(float(w.sum())-144), 'input_payloads': [finite['payload'], spec]}


@lru_cache(maxsize=28)
def authenticated_cases(root, group):
    if group not in GROUPS: raise ValueError('only the28 remaining groups are owned by this batch')
    if group == 13: return {'middle': _group13_middle(root)}
    finite = json.loads((Path(root)/'results/development/nsc-incoming-spectral-source.json').read_text())
    right = finite['groups'][str(group)]['signed_inventory']['1']['upper_energy_endpoint']
    return {'low': authenticated_window(root, group, 'low_vacuum', 32., 40.),
            'middle': authenticated_window(root, group, 'middle_correction', 40., right)}


def prepare_group_source(root, group):
    """One resolution, two precisions; no radial or physical mode calculation."""
    cases = authenticated_cases(str(Path(root)), group); arrays = {}; windows = {}; exact = {}
    for name, owned in cases.items():
        windows[name] = {'group': group, 'kind': owned['kind'], 'interval': [owned['left'], owned['right']],
            'source_order': 24, 'factor_per_sign': owned['factor'], 'selected_panels': owned['descriptors'],
            'input_payloads': owned['input_payloads'], 'original_measure_residual': owned['measure_residual'], 'reports': {}}
    for dps in (50, 70):
        with mp.workdps(dps):
            channel = next(iter(cases.values()))['channel']; mass = mp.mpf(channel['compact_mass'])
            for sign in next(iter(cases.values()))['signs']:
                angular_float = sign*channel['angular_eigenvalue']; angular = mp.mpf(angular_float)
                c, prefix_error = channel_coefficients(dps, channel['compact_mass'], angular_float)
                for name, owned in cases.items():
                    report = windows[name]; selected = owned['selected'][str(sign)]
                    rules = {'refined64': cached_grid(dps, owned['left'], owned['right'], 64)}
                    if dps == 70:
                        rules['refined48'] = cached_grid(dps, owned['left'], owned['right'], 48)
                        rules['original'] = (tuple(mp.mpf(v) for v in selected['energies']), tuple(mp.mpf(v) for v in selected['weights']))
                    for label, (E, w) in rules.items():
                        kernels = [horner_kernel(e, mass, angular, c, owned['kind'], dps=dps) for e in E]
                        physical = mp.mpf(owned['factor'])*mp.matrix([sum((wi*K[i] for wi, K in zip(w, kernels)), mp.mpf(0)) for i in range(4)])
                        key = f'{name}/{sign}/{label}/dps{dps}'; exact[key] = physical
                        arrays[key+'/energies'] = np.array([float(v) for v in E]); arrays[key+'/weights'] = np.array([float(v) for v in w])
                        arrays[key+'/kernels'] = np.array([[float(v) for v in K] for K in kernels])
                        report['reports'][f'{sign}/{label}/dps{dps}'] = {'physical_value': [float(v) for v in physical],
                            'physical_value_decimal': [mp.nstr(v, dps) for v in physical], 'rows': len(E),
                            'prefix16_residual': prefix_error, 'measure_residual': float(abs(sum(w)-(owned['right']-owned['left'])))}
                    if dps == 70:
                        for key, value in selected.items(): arrays[f'{name}/{sign}/archived/'+key] = value.copy()
                        errors = []
                        for fraction in (mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('0.9')):
                            E = owned['left']+fraction*(owned['right']-owned['left'])
                            new = horner_kernel(E, mass, angular, c, owned['kind'], dps=dps)
                            original = window_kernel(E, mass, angular, c, owned['kind'])
                            errors.append(max(abs(a-b) for a, b in zip(new, original)))
                        report['reports'][str(sign)] = {'horner_vs_frozen_kernel': float(max(errors)),
                            'source16_convention_residual': (float(np.max(abs(source16_control(selected['energies'], mass, angular, c[:16])-selected['vacuum_kernels'])))
                                                            if owned['kind'] == 'middle_correction' else None)}
        print(f'source batch group{group}: precision{dps}, actual signs and windows complete', flush=True)
    with mp.workdps(70):
        for name, owned in cases.items():
            indicators = {}
            for label, first, second in (('original_to48','original/dps70','refined48/dps70'),
                                         ('quadrature48_to64','refined48/dps70','refined64/dps70'),
                                         ('precision50_to70','refined64/dps50','refined64/dps70')):
                indicators[label] = [float(sum(abs(exact[f'{name}/{s}/{first}'][i]-exact[f'{name}/{s}/{second}'][i]) for s in owned['signs'])) for i in range(4)]
            windows[name]['indicators_sum_absolute_signs'] = indicators
    return arrays, {'group': group, 'windows': windows, 'horner_identity': horner_identity()}


def replay_group_source(root, arrays, prepared):
    """Authenticate and contract arrays only; never prepare missing data."""
    group = prepared['group']; cases = authenticated_cases(str(Path(root)), group)
    if set(prepared['windows']) != set(cases) or prepared['horner_identity'] != horner_identity():
        raise ValueError('actual windows or exact Horner identity changed')
    windows = {}; failures = {}
    for name, owned in cases.items():
        report = prepared['windows'][name]
        if report['source_order'] != 24 or report['factor_per_sign'] != owned['factor'] or report['selected_panels'] != owned['descriptors'] or report['input_payloads'] != owned['input_payloads']:
            raise ValueError('prepared source labels, factor or provenance changed')
        old = np.zeros(4); old_thermal = np.zeros(4); value = np.zeros(4); per_sign = {}; residuals = dict.fromkeys(TOLERANCES, 0.)
        residuals['original_measure'] = report['original_measure_residual']
        for sign in owned['signs']:
            selected = owned['selected'][str(sign)]
            for key, original in selected.items():
                if not np.array_equal(arrays[f'{name}/{sign}/archived/'+key], original): raise ValueError('old source array changed')
            old_sign = owned['factor']*np.einsum('n,nv->v', selected['weights'], selected['kernels']); old += old_sign
            if owned['kind'] == 'middle_correction': old_thermal += owned['factor']*np.einsum('n,nv->v', selected['weights'], selected['thermal_kernels'])
            new_sign = None
            for label, dps in (('original',70), ('refined48',70), ('refined64',50), ('refined64',70)):
                key = f'{name}/{sign}/{label}/dps{dps}'; row = report['reports'][f'{sign}/{label}/dps{dps}']
                E, w, K = [arrays[key+'/'+field] for field in ('energies','weights','kernels')]
                if K.shape != (len(E),4) or np.any(w <= 0) or np.any(E <= owned['left']) or np.any(E >= owned['right']) or np.count_nonzero(K[:,2]):
                    raise ValueError('positive owned window and exact vacuum-current identity required')
                if label == 'original' and (not np.array_equal(E, selected['energies']) or not np.array_equal(w, selected['weights'])):
                    raise ValueError('original source quadrature changed')
                replay = owned['factor']*np.einsum('n,nv->v', w, K)
                residuals['contraction_replay'] = max(residuals['contraction_replay'], float(np.max(abs(replay-row['physical_value']))))
                residuals['prefix16'] = max(residuals['prefix16'], row['prefix16_residual'])
                residuals['quadrature_measure'] = max(residuals['quadrature_measure'], row['measure_residual'])
                if label == 'refined64' and dps == 70: new_sign = np.array(row['physical_value']); value += new_sign
            control = report['reports'][str(sign)]
            residuals['horner_vs_frozen_kernel'] = max(residuals['horner_vs_frozen_kernel'], control['horner_vs_frozen_kernel'])
            residuals['source16_convention'] = max(residuals['source16_convention'], control['source16_convention_residual'] or 0.)
            sign_delta = new_sign-old_sign if owned['kind'] == 'low_vacuum' else new_sign
            per_sign[str(sign)] = {'old_source': old_sign.tolist(), 'new_numerical_value': new_sign.tolist(), 'explicit_source_delta': sign_delta.tolist()}
        residuals.update({key:max(v) for key,v in report['indicators_sum_absolute_signs'].items()})
        bad = {key:v for key,v in residuals.items() if not np.isfinite(v) or v>TOLERANCES[key]}
        if bad: failures[name] = bad
        delta = value-old if owned['kind']=='low_vacuum' else value
        updated = value if owned['kind']=='low_vacuum' else old+value
        thermal = thermal_middle_difference_bound(owned['channel'], owned['config'], owned['left'], owned['right'])
        windows[name] = {'group':group, 'kind':owned['kind'], 'interval':[owned['left'],owned['right']], 'source_order':24,
            'factor_per_sign':owned['factor'], 'selected_panels':owned['descriptors'], 'actual_angular_signs':list(owned['signs']),
            'per_sign':per_sign, 'old_source':old.tolist(), 'new_numerical_value':value.tolist(), 'explicit_source_delta':delta.tolist(),
            'updated_source_approximant':updated.tolist(), 'retained_archived_thermal':old_thermal.tolist() if owned['kind']=='middle_correction' else None,
            'thermal_policy':'old thermal retained; its difference is bounded' if owned['kind']=='middle_correction' else 'physical thermal retained as nonzero bounded remainder',
            'thermal_scattering_stress_upper':thermal, 'physical_order24_error_bound':None, 'physical_error_status':'OPEN: independent matching radial certificate required',
            'numerical_indicators_sum_absolute_signs':report['indicators_sum_absolute_signs'], 'residuals':residuals,
            'numerical_tolerances':dict(TOLERANCES), 'failures':bad, 'status':'PASS: numerical source controls; physical certificate separate' if not bad else 'OPEN'}
    total = sum((np.array(v['explicit_source_delta']) for v in windows.values()), np.zeros(4))
    return {'group':group, 'windows':windows, 'explicit_source_delta':total.tolist(), 'failures':failures,
            'status':'PASS: group source controls; physical mode bound OPEN' if not failures else 'OPEN',
            'horner_identity':prepared['horner_identity'], 'new_radial_mode_or_scattering_runs':0}
