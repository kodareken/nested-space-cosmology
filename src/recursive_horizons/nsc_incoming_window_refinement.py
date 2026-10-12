"""Reusable explicit source windows; radial certificates are never regenerated."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from .nsc_incoming_source_assembly import read_panels
from .nsc_incoming_source_tail import _riccati_at_one, _mp_bloch
from .nsc_incoming_state_moments import incoming_group_factor
from .nsc_incoming_middle_order_correction import stable_delta_bloch, source16_control
from .nsc_incoming_projector_energy_bound import projected_energy_bound


KINDS = ('low_vacuum', 'middle_correction')


def authenticated_window(root, group, kind, left, right):
    """Select an actual LOW cell or whole archived middle interval."""
    if kind not in KINDS or group not in range(1, 33) or group in (13, 14) or not 0 < left < right:
        raise ValueError('owned common retained LOW/middle window required')
    root = Path(root); panels, panel_meta = read_panels(root)
    finite_path = 'results/development/nsc-incoming-spectral-source.json'
    retained_path = 'results/development/nsc-pg-retained-covariance.json'
    finite = json.loads((root/finite_path).read_text()); retained = json.loads((root/retained_path).read_text())
    channels = json.loads((root/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    config = json.loads((root/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    for record in (finite, retained):
        for field in ('source_hashes', 'input_hashes', 'numerical_owner_and_input_hashes'):
            for path, expected in record.get(field, {}).items():
                if sha256((root/path).read_bytes()).hexdigest() != expected: raise ValueError('window owner/input changed: '+path)
        item = record['payload']
        if sha256((root/item['path']).read_bytes()).hexdigest() != item['sha256']: raise ValueError('window artifact changed')
    channel = channels[group]; signs = (1,) if channel['angular_eigenvalue'] == 0 else (1, -1)
    factor = incoming_group_factor(channel, signs)
    if factor != finite['groups'][str(group)]['factor']: raise ValueError('actual group factor changed')
    selected = {}; descriptors = []; measure = 0.
    with np.load(root/retained['payload']['path'], allow_pickle=False) as original:
        for sign in signs:
            inventory = finite['groups'][str(group)]['signed_inventory'][str(sign)]
            starts = ('low/', 'low_ref/') if kind == 'low_vacuum' else ('mid/', 'mid_ref/')
            names = [p['panel'] for p in inventory['pieces'] if p['panel'].startswith(starts)]
            if len(names) != 1: raise ValueError('one selected common panel per signed window required')
            name = names[0]; meta = json.loads(original[name+'/metadata_json'].tobytes())
            if meta['channel'] != channel or meta.get('sign', meta.get('angular_sign')) != sign:
                raise ValueError('signed panel labels differ from current source inventory')
            if kind == 'low_vacuum':
                if (float(left), float(right)) not in list(zip(meta['edges'], meta['edges'][1:])):
                    raise ValueError('requested LOW interval must be one actual archived cell')
                raw_E, raw_w = original[name+'/energy'].ravel(), original[name+'/weight'].ravel()
            else:
                if [meta['left'], meta['right']] != [left, right]:
                    raise ValueError('requested middle interval differs from actual archived endpoints')
                raw_E, raw_w = original[name+'/energies'].ravel(), original[name+'/weights'].ravel()
            if not np.array_equal(raw_E, panels[name+'/energies']) or not np.array_equal(raw_w, panels[name+'/weights']):
                raise ValueError('source nodes/weights differ from original archived producer')
            mask = (raw_E > left)&(raw_E < right); E, w = raw_E[mask], raw_w[mask]
            if not len(E) or np.any(w <= 0): raise ValueError('positive existing source quadrature required')
            measure = max(measure, abs(float(w.sum())-(right-left)))
            if measure > 3e-11: raise ValueError('selected source quadrature does not cover the requested interval')
            quantities = ('energies', 'weights', 'kernels', 'arithmetic_indicator')
            if kind == 'middle_correction': quantities += ('vacuum_kernels', 'thermal_kernels')
            selected[str(sign)] = {k: panels[name+'/'+k][mask].copy() for k in quantities}
            descriptors.append({'sign': sign, 'panel': name, 'rows': len(E), 'interval': [float(left), float(right)],
                                'method': panel_meta['panels'][name]['method']})
    if (right-left) % 8 != 0: raise ValueError('owned eight-unit source window partition required')
    return {'group': group, 'kind': kind, 'left': float(left), 'right': float(right), 'channel': channel,
            'config': config, 'signs': signs, 'factor': factor, 'selected': selected,
            'descriptors': descriptors, 'measure_residual': measure,
            'input_payloads': [finite['payload'], retained['payload']]}


def window_kernel(energy, mass, angular, coefficients, kind):
    """Full ad4 vacuum or stable order24-minus16 vacuum, same raw vertices."""
    E = mp.mpf(energy)
    if len(coefficients) != 24 or E <= 0 or kind not in KINDS:
        raise ValueError('explicit order24 and a declared positive-energy source kind required')
    a, r = mp.sqrt(3*mp.pi/2-4), mp.sqrt(2)
    if kind == 'middle_correction':
        S = sum((coefficients[j]/(2*E)**(j+1) for j in range(16)), mp.mpc(0))
        dS = sum((coefficients[j]/(2*E)**(j+1) for j in range(16, 24)), mp.mpc(0))
        delta = stable_delta_bloch(S, dS, a)
    else:
        S = (angular/r+1j*mass)/(2*E)+sum((coefficients[j]/(2*E)**(j+1) for j in range(1, 24)), mp.mpc(0))
        u = a*a*abs(S)**2; y = a*a*(mass*mass+angular*angular/(r*r))/(E*E); root = mp.sqrt(1+y)
        delta = [2*a*mp.im(S)/(1+u)-mass*a/(E*root), -2*a*mp.re(S)/(1+u)+angular*a/(r*E*root),
                 y/(root*(1+root))-2*u/(1+u)]
        terms = _mp_bloch()(3*mp.pi/4, mass, angular, E)
        for i in range(3): delta[i] -= sum((terms[j][i, 0] for j in range(1, 5)), mp.mpf(0))
    parallel = -E/a*delta[2]
    return mp.matrix([-mass*delta[0]+angular/r*delta[1]+parallel, parallel, mp.mpf(0), angular/(2*r)*delta[1]])


def window_quadrature(left, right, points):
    if points not in (48, 64) or not 0 < left < right or (right-left)%8:
        raise ValueError('positive48/64 quadrature on owned eight-unit cells required')
    x, w = mp.gauss_quadrature(points, 'legendre')
    starts = [mp.mpf(left)+8*i for i in range(int((right-left)/8))]
    return [start+4*(xi+1) for start in starts for xi in x], [4*wi for _ in starts for wi in w]


def prepare_window(root, group, kind, left, right):
    """Source algebra only; independent of any radial-certificate producer."""
    owned = authenticated_window(root, group, kind, left, right); arrays = {}; reports = {}; exact = {}
    for dps in (50, 70):
        with mp.workdps(dps):
            grids = {'refined64': window_quadrature(left, right, 64)}
            if dps == 70: grids['refined48'] = window_quadrature(left, right, 48)
            mass = mp.mpf(owned['channel']['compact_mass']); magnitude = mp.mpf(owned['channel']['angular_eigenvalue'])
            factor = mp.mpf(owned['factor'])
            for sign in owned['signs']:
                ell = sign*magnitude; c = _riccati_at_one(mass, ell, order=24); c16 = _riccati_at_one(mass, ell, order=16)
                prefix = float(max(abs(a-b) for a, b in zip(c[:16], c16)))
                selected = owned['selected'][str(sign)]; rules = dict(grids)
                if dps == 70: rules['original'] = ([mp.mpf(v) for v in selected['energies']], [mp.mpf(v) for v in selected['weights']])
                for label, (E, w) in rules.items():
                    kernels = [window_kernel(e, mass, ell, c, kind) for e in E]
                    physical = factor*mp.matrix([sum((wi*K[i] for wi, K in zip(w, kernels)), mp.mpf(0)) for i in range(4)])
                    key = f'{sign}/{label}/dps{dps}'; exact[key] = physical
                    arrays[key+'/energies'] = np.array([float(v) for v in E]); arrays[key+'/weights'] = np.array([float(v) for v in w])
                    arrays[key+'/kernels'] = np.array([[float(v) for v in K] for K in kernels])
                    reports[key] = {'physical_value': [float(v) for v in physical],
                        'physical_value_decimal': [mp.nstr(v, dps) for v in physical], 'rows': len(E),
                        'prefix16_residual': prefix, 'measure_residual': float(abs(sum(w)-(right-left)))}
                if dps == 70:
                    for key, value in selected.items(): arrays[f'{sign}/archived/'+key] = value.copy()
                    reports[str(sign)] = {'old16_source_convention_residual': (float(np.max(abs(
                        source16_control(selected['energies'], mass, ell, c16)-selected['vacuum_kernels'])))
                        if kind == 'middle_correction' else None)}
            print(f'window group{group} {kind} [{left:g},{right:g}] precision{dps}; source algebra only', flush=True)
    with mp.workdps(70):
        comparisons = {}
        for name, x, y in (('original_to48', 'original/dps70', 'refined48/dps70'),
                           ('quadrature48_to64', 'refined48/dps70', 'refined64/dps70'),
                           ('precision50_to70', 'refined64/dps50', 'refined64/dps70')):
            comparisons[name] = [float(sum(abs(exact[f'{s}/{x}'][i]-exact[f'{s}/{y}'][i]) for s in owned['signs'])) for i in range(4)]
    return arrays, {'group': group, 'kind': kind, 'interval': [float(left), float(right)], 'source_order': 24,
        'factor_per_sign': owned['factor'], 'selected_panels': owned['descriptors'], 'reports': reports,
        'indicators_sum_absolute_signs': comparisons, 'original_measure_residual': owned['measure_residual'],
        'input_payloads': owned['input_payloads']}, owned


def bind_window_certificate(channel, radial, cross, left, right):
    """Explicit fresh integration view; original energy-free coefficients kept."""
    if radial['group'] != channel['index'] or radial['physical_Riccati_order'] != 24 or cross['physical_Riccati_order'] != 24:
        raise ValueError('same channel and matching order24 source/certificate required')
    if not 0 < left < right: raise ValueError('positive fresh finite interval required')
    before = deepcopy(radial)
    view = dict(radial, lower=float(right))
    energy = projected_energy_bound(channel, view, cross, float(left), float(right))
    if radial != before: raise ArithmeticError('original radial certificate was mutated')
    return {'original_endpoint_tag': radial['lower'], 'fresh_interval': [float(left), float(right)],
            'consumer_endpoint_tag': view['lower'], 'coefficient_integrals_unchanged': True,
            'new_radial_recurrence': False,
            'reason': 'saved radial defect integrals have no energy dependence; only finite-energy primitives change',
            'energy': energy}
