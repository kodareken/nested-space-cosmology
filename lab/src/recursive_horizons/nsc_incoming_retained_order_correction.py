"""Explicit group12 order24-minus16 source on its owned middle interval."""
from pathlib import Path
import json

import mpmath as mp
import numpy as np

from .nsc_incoming_middle_bound import authenticated_middle_inputs
from .nsc_incoming_middle_order_correction import stable_delta_bloch, source16_control
from .nsc_incoming_source_tail import _riccati_at_one
from .nsc_incoming_state_moments import incoming_group_factor


def owned_middle_source(root, group=12):
    if group != 12:
        raise ValueError('this correction record owns group12 only')
    root = Path(root); owned = authenticated_middle_inputs(root)
    channel = owned['channels'][group]; interval = owned['middle_intervals'][group-1]
    if (interval['left'], interval['right']) != (40., 320.):
        raise ValueError('owned group12 interval must stay[40,320]')
    finite = json.loads((root/'results/development/nsc-incoming-spectral-source.json').read_text())
    names = [row['panel'] for row in interval['signed_selected_panels']]
    with np.load(root/finite['payload']['path'], allow_pickle=False) as f:
        arrays = {k: f[k].copy() for k in f.files if any(k.startswith(n+'/') for n in names)}
    factor = incoming_group_factor(channel, (1, -1))
    if factor != finite['groups'][str(group)]['factor']:
        raise ValueError('same signed source multiplicity required')
    return {**owned, 'channel': channel, 'interval': interval, 'arrays': arrays, 'factor': factor,
            'thermal_scattering_bound': json.loads((root/'results/development/nsc-incoming-middle-bound.json').read_text())['groups'][group-1]['thermal_scattering_stress_error_upper']}


def delta_kernel(energy, mass, angular, coefficients):
    if len(coefficients) != 24 or not 40 <= energy <= 320:
        raise ValueError('explicit24 coefficients on group12[40,320] required')
    E = mp.mpf(energy); a, r = mp.sqrt(3*mp.pi/2-4), mp.sqrt(2)
    before = sum((coefficients[j]/(2*E)**(j+1) for j in range(16)), mp.mpc(0))
    increment = sum((coefficients[j]/(2*E)**(j+1) for j in range(16, 24)), mp.mpc(0))
    db = stable_delta_bloch(before, increment, a); parallel = -E/a*db[2]
    return mp.matrix([-mass*db[0]+angular/r*db[1]+parallel, parallel, 0, angular/(2*r)*db[1]])


def refined_grid():
    x, w = mp.gauss_quadrature(48, 'legendre')
    return ([start+4*(node+1) for start in range(40,320,8) for node in x],
            [4*weight for _ in range(40,320,8) for weight in w])


def prepare_retained_correction(root):
    owned = owned_middle_source(root); channel = owned['channel']
    arrays, reports, integrals = {}, {}, {}
    for precision in (50,70):
        with mp.workdps(precision):
            mass, angular = mp.mpf(channel['compact_mass']), mp.mpf(channel['angular_eigenvalue'])
            factor = mp.mpf(owned['factor']); new_E, new_w = refined_grid()
            for selected in owned['interval']['signed_selected_panels']:
                sign, prefix = selected['sign'], selected['panel']+'/'
                c24 = _riccati_at_one(mass, sign*angular, order=24)
                c16 = _riccati_at_one(mass, sign*angular, order=16)
                prefix_residual = max(abs(x-y) for x,y in zip(c24[:16],c16))
                choices = [('refined48', new_E, new_w)]
                if precision == 70:
                    choices.append(('archived', [mp.mpf(v) for v in owned['arrays'][prefix+'energies']],
                                                [mp.mpf(v) for v in owned['arrays'][prefix+'weights']]))
                for label,E,w in choices:
                    kernels = [delta_kernel(e,mass,sign*angular,c24) for e in E]
                    integral = [factor*sum((wi*k[i] for wi,k in zip(w,kernels)),mp.mpf(0)) for i in range(4)]
                    key = f'{sign}/{label}/dps{precision}'; integrals[key] = integral
                    arrays[key+'/energies'] = np.array([float(e) for e in E])
                    arrays[key+'/weights'] = np.array([float(v) for v in w])
                    arrays[key+'/kernels'] = np.array([[float(v) for v in k] for k in kernels])
                    reports[key] = {'correction': [float(v) for v in integral],
                                    'correction_decimal': [mp.nstr(v, precision) for v in integral],
                                    'coefficient_prefix_residual': float(prefix_residual),
                                    'interval_measure_residual': float(abs(sum(w)-280))}
                if precision == 70:
                    for label in ('energies','weights','kernels','vacuum_kernels','thermal_kernels'):
                        arrays[f'{sign}/original16/'+label] = owned['arrays'][prefix+label].copy()
                    control = source16_control(owned['arrays'][prefix+'energies'],mass,sign*angular,c16)
                    reports[str(sign)] = {'source16_vacuum_convention': float(np.max(abs(control-owned['arrays'][prefix+'vacuum_kernels']))),
                                         'source16_decomposition': float(np.max(abs(owned['arrays'][prefix+'kernels']-owned['arrays'][prefix+'vacuum_kernels']-owned['arrays'][prefix+'thermal_kernels'])))}
        print(f'group12 additive correction precision{precision}; actual interval and source retained',flush=True)
    with mp.workdps(70):
        quad = [sum(abs(integrals[f'{s}/refined48/dps70'][i]-integrals[f'{s}/archived/dps70'][i]) for s in(1,-1)) for i in range(4)]
        prec = [sum(abs(integrals[f'{s}/refined48/dps70'][i]-integrals[f'{s}/refined48/dps50'][i]) for s in(1,-1)) for i in range(4)]
    return arrays, {'group':12, 'interval':owned['interval'], 'old_order':16, 'new_order':24,
                    'factor_per_sign':owned['factor'], 'reports':reports,
                    'quadrature_indicator_not_bound':[float(v) for v in quad],
                    'precision_indicator_not_bound':[float(v) for v in prec],
                    'thermal_scattering_bound':owned['thermal_scattering_bound'],
                    'input_payloads':owned['input_payloads']}
