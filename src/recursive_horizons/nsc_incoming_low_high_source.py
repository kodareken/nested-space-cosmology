"""One explicit order24 source improvement: group22, existing LOW[32,40]."""
import json
from pathlib import Path
from time import perf_counter

import mpmath as mp
import numpy as np

from .nsc_incoming_middle_bound import authenticated_middle_inputs, _parameters
from .nsc_incoming_source_tail import _riccati_at_one, _mp_bloch
from .nsc_incoming_state_moments import incoming_group_factor
from .nsc_incoming_centered_order24 import CenteredOrder24Geometry, trimmed_defect_jets
from .nsc_incoming_defect_taylor_bound import centered_coefficient_enclosure
from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _geometry_jets, _up_float


GROUP, ORDER, LEFT, RIGHT = 22, 24, 32., 40.
SIGNS = (1, -1)


def authenticated_low_inputs(root):
    root = Path(root); inherited = authenticated_middle_inputs(root)
    finite = json.loads((root/'results/development/nsc-incoming-spectral-source.json').read_text())
    with np.load(root/finite['payload']['path'], allow_pickle=False) as data:
        arrays = {k: data[k].copy() for k in data.files
                  if any(k.startswith(f'{family}/{GROUP}_{s}/') for family in ('low', 'low_ref') for s in SIGNS)}
    retained = json.loads((root/'results/development/nsc-pg-retained-covariance.json').read_text())
    with np.load(root/retained['payload']['path'], allow_pickle=False) as data:
        metadata = {f'{family}/{s}': json.loads(data[f'{family}/{GROUP}_{s}/metadata_json'].tobytes())
                    for family in ('low', 'low_ref') for s in SIGNS}
    channel = inherited['channels'][GROUP]; factor = incoming_group_factor(channel, SIGNS)
    if factor != finite['groups'][str(GROUP)]['factor']:
        raise ValueError('unchanged group22 source normalization required')
    selected = {}; measure_error = 0.
    for family, points in (('low', 24), ('low_ref', 32)):
        for sign in SIGNS:
            meta = metadata[f'{family}/{sign}']; prefix = f'{family}/{GROUP}_{sign}/'
            if (meta['channel'] != channel or meta['sign'] != sign or meta['upper'] != RIGHT
                    or (LEFT, RIGHT) not in list(zip(meta['edges'], meta['edges'][1:]))):
                raise ValueError('the actual group22 signed LOW[32,40] cell must be retained')
            E = arrays[prefix+'energies']; mask = (E > LEFT)&(E < RIGHT)
            w = arrays[prefix+'weights'][mask]
            if mask.sum() != points or np.any(w <= 0):
                raise ValueError('the original positive24/32-node LOW cell is required')
            measure_error = max(measure_error, abs(float(w.sum())-(RIGHT-LEFT)))
            if measure_error > 3e-13: raise ValueError('original LOW cell measure differs from its endpoints')
            selected[f'{sign}/old{points}'] = {key: arrays[prefix+key][mask].copy()
                                              for key in ('energies', 'weights', 'kernels', 'arithmetic_indicator')}
    return {'channel': channel, 'config': inherited['config'], 'selected': selected,
            'factor': factor, 'metadata': metadata, 'measure_residual': measure_error,
            'input_payloads': inherited['input_payloads']}


def vacuum_source_kernel(energy, mass, angular, coefficients):
    """Owned ad4 subtraction of the explicit local order24 vacuum projector."""
    E = mp.mpf(energy)
    if len(coefficients) != ORDER or not LEFT <= E <= RIGHT:
        raise ValueError('explicit order24 source on group22 LOW[32,40] required')
    a, r = mp.sqrt(3*mp.pi/2-4), mp.sqrt(2)
    # Preserve the existing leading identity at the active precision.
    series = (angular/r+1j*mass)/(2*E)
    series += sum((coefficients[j]/(2*E)**(j+1) for j in range(1, ORDER)), mp.mpc(0))
    u = a*a*abs(series)**2; y = a*a*(mass*mass+angular*angular/(r*r))/(E*E); root = mp.sqrt(1+y)
    delta = [2*a*mp.im(series)/(1+u)-mass*a/(E*root),
             -2*a*mp.re(series)/(1+u)+angular*a/(r*E*root),
             y/(root*(1+root))-2*u/(1+u)]
    reference = _mp_bloch()(3*mp.pi/4, mass, angular, E)
    for i in range(3): delta[i] -= sum((reference[j][i, 0] for j in range(1, 5)), mp.mpf(0))
    parallel = -E/a*delta[2]
    return mp.matrix([-mass*delta[0]+angular/r*delta[1]+parallel,
                      parallel, mp.mpf(0), angular/(2*r)*delta[1]])


def positive_grid(points):
    if points not in (48, 64): raise ValueError('bounded positive48/64 source quadratures required')
    x, w = mp.gauss_quadrature(points, 'legendre')
    E, weights = [LEFT+4*(v+1) for v in x], [4*v for v in w]
    if any(v <= 0 for v in weights): raise ArithmeticError('positive source quadrature required')
    return E, weights


def prepare_source(root):
    owned = authenticated_low_inputs(root); channel = owned['channel']; arrays = {}; reports = {}; precise = {}
    for precision in (50, 70):
        with mp.workdps(precision):
            grids = {'new64': positive_grid(64)}
            if precision == 70: grids['new48'] = positive_grid(48)
            mass = mp.mpf(channel['compact_mass']); magnitude = mp.mpf(channel['angular_eigenvalue']); factor = mp.mpf(owned['factor'])
            for sign in SIGNS:
                angular = sign*magnitude; coeff = _riccati_at_one(mass, angular, order=24)
                old_coeff = _riccati_at_one(mass, angular, order=16)
                prefix_error = max(abs(a-b) for a, b in zip(coeff[:16], old_coeff))
                sign_grids = dict(grids)
                if precision == 70:
                    for points in (24, 32):
                        old = owned['selected'][f'{sign}/old{points}']
                        sign_grids[f'old{points}'] = ([mp.mpf(v) for v in old['energies']], [mp.mpf(v) for v in old['weights']])
                        for key, value in old.items(): arrays[f'{sign}/archived{points}/'+key] = value
                for label, (E, w) in sign_grids.items():
                    kernels = [vacuum_source_kernel(e, mass, angular, coeff) for e in E]
                    value = factor*mp.matrix([sum((weight*K[i] for weight, K in zip(w, kernels)), mp.mpf(0)) for i in range(4)])
                    key = f'{sign}/{label}/dps{precision}'; precise[key] = value
                    arrays[key+'/energies'] = np.array([float(v) for v in E])
                    arrays[key+'/weights'] = np.array([float(v) for v in w])
                    arrays[key+'/vacuum_kernels'] = np.array([[float(v) for v in K] for K in kernels])
                    reports[key] = {'physical_vacuum_source': [float(v) for v in value],
                        'physical_vacuum_source_decimal': [mp.nstr(v, precision) for v in value],
                        'coefficient16_prefix_residual': float(prefix_error), 'rows': len(E),
                        'quadrature_measure_residual': float(abs(sum(w)-8))}
            print(f'group22 LOW[32,40] local source, both signs, precision{precision}', flush=True)
    with mp.workdps(70):
        indicators = {}
        for name, first, second in (('quadrature32_to48', 'old32/dps70', 'new48/dps70'),
                                     ('quadrature48_to64', 'new48/dps70', 'new64/dps70'),
                                     ('precision50_to70', 'new64/dps50', 'new64/dps70')):
            indicators[name] = [float(sum(abs(precise[f'{s}/{first}'][i]-precise[f'{s}/{second}'][i]) for s in SIGNS)) for i in range(4)]
    return arrays, {'reports': reports, 'indicators_sum_absolute_signs': indicators,
                    'source_factor_per_sign': owned['factor'], 'input_payloads': owned['input_payloads'],
                    'owned_interval_measure_residual': owned['measure_residual']}, owned


def prepare_group22_defect(config, channel, *, progress=None, budget_seconds=900):
    """Wrapper only: frozen trimmed helper, fixed128 cells and Taylor depth4."""
    if channel.get('index') != GROUP or not 0 < budget_seconds <= 1200:
        raise ValueError('one group22 certificate and bounded runtime required')
    geometry = CenteredOrder24Geometry(config, precision=40); start = perf_counter()
    with _precision(40):
        _, _, mass, angular, _ = _parameters(channel)
        totals = {s: [mp.iv.mpf(0) for _ in range(25)] for s in SIGNS}
        for index, (rho, center) in enumerate(geometry.cells):
            if perf_counter()-start > budget_seconds: raise TimeoutError('one fixed group22 interval-certificate budget reached')
            whole_geometry = _geometry_jets(rho, 28); center_geometry = _geometry_jets(center, 28)
            for sign in SIGNS:
                whole = trimmed_defect_jets(*whole_geometry, mass, sign*angular, order=24, depth=4)
                point = trimmed_defect_jets(*center_geometry, mass, sign*angular, order=24, depth=4)
                for j, (w, p) in enumerate(zip(whole, point)):
                    enclosure = centered_coefficient_enclosure(w, p, rho-center, 4)
                    totals[sign][j] += geometry.weight*mp.iv.mpf(_hi(abs(enclosure)))/128
            if progress and (index+1)%16 == 0: progress(index+1)
        return {'group': GROUP, 'lower': RIGHT, 'intervals': 128,
            'physical_Riccati_order': ORDER, 'centered_remainder_depth': 4,
            'first_inverse_power': 24, 'last_inverse_power': 48, 'directed_precision': 40,
            'per_sign': [{'sign': s, 'coefficient_integral_upper': [_up_float(v) for v in totals[s]]} for s in SIGNS],
            'horizon_enclosure': [mp.nstr(_lo(geometry.horizon), 55), mp.nstr(_hi(geometry.horizon), 55)],
            'endpoint_weight_upper': _up_float(geometry.weight), 'initial_projector_difference': 0.,
            'initial_scope': 'same declared affine-horizon vacuum; explicit local order24 projector',
            'finite_offset_archive_accuracy_certified': False,
            'source_and_certificate_order_match': True,
            'bound_method': 'frozen trimmed normalized jets, centered depth4 remainder, fixed128 cells'}
