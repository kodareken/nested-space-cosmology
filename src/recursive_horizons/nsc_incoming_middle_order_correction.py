"""Explicit group14 order24-minus16 vacuum correction; no field evolution."""
from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp
import numpy as np
from numpy.polynomial.legendre import leggauss

from .nsc_incoming_middle_bound import authenticated_middle_inputs
from .nsc_incoming_source_tail import _riccati_at_one
from .nsc_incoming_state_moments import incoming_group_factor, incoming_adiabatic_reference, AXIAL, RADIUS, SIGMA


OLD_ORDER, NEW_ORDER = 16, 24
LEFT, RIGHT = 16., 160.
SIGNS = (1, -1)


def authenticated_group14(root):
    """Read the selected original grids, moments and physical source labels."""
    root = Path(root); inputs = authenticated_middle_inputs(root)
    finite = json.loads((root/'results/development/nsc-incoming-spectral-source.json').read_text())
    path = root/finite['payload']['path']
    with np.load(path, allow_pickle=False) as data:
        arrays = {k: data[k].copy() for k in data.files if k.startswith('group14/mid24_')}
    channel = inputs['channels'][14]; interval = inputs['middle_intervals'][13]
    if interval['group'] != 14 or (interval['left'], interval['right']) != (LEFT, RIGHT):
        raise ValueError('only the existing group14 middle interval is owned')
    factor = incoming_group_factor(channel, SIGNS)
    if factor != finite['groups']['14']['factor']:
        raise ValueError('source group factor changed')
    x, w = leggauss(24)
    expected_E = np.concatenate([start+4*(x+1) for start in np.arange(LEFT, RIGHT, 8.)])
    expected_w = np.tile(4*w, 18)
    grid_residual = 0.
    for sign in SIGNS:
        prefix = f'group14/mid24_{sign}/'
        E, weights = arrays[prefix+'energies'], arrays[prefix+'weights']
        if E.shape != (432,) or weights.shape != E.shape or np.any(weights <= 0):
            raise ValueError('positive original24-point eight-unit panels required')
        grid_residual = max(grid_residual, float(np.max(abs(E-expected_E))), float(np.max(abs(weights-expected_w))))
        if grid_residual > 3e-13:
            raise ValueError('original positive quadrature cells differ from authenticated middle recipe')
    middle = json.loads((root/'results/development/nsc-incoming-middle-bound.json').read_text())
    for field in ('source_hashes', 'input_hashes'):
        for p, expected in middle[field].items():
            if sha256((root/p).read_bytes()).hexdigest() != expected:
                raise ValueError('committed middle bound dependency changed: '+p)
    thermal = middle['groups'][13]['thermal_scattering_stress_error_upper']
    return {'channel': channel, 'config': inputs['config'], 'arrays': arrays, 'factor': factor,
            'interval': interval, 'grid_recipe_residual': grid_residual,
            'thermal_scattering_stress_error_upper': thermal,
            'input_payloads': inputs['input_payloads']}


def stable_delta_bloch(series16, delta_series, axial):
    """P24-P16 in the owned Bloch convention, without unit-trace subtraction."""
    u = axial*axial*abs(series16)**2; D = 1+u
    du = axial*axial*(2*mp.re(mp.conj(series16)*delta_series)+abs(delta_series)**2)
    denominator = D*(D+du)
    return mp.matrix([2*axial*(D*mp.im(delta_series)-mp.im(series16)*du)/denominator,
                     -2*axial*(D*mp.re(delta_series)-mp.re(series16)*du)/denominator,
                     -2*du/denominator])


def correction_kernel(energy, mass, angular, coefficients):
    if len(coefficients) != NEW_ORDER or not LEFT <= energy <= RIGHT:
        raise ValueError('explicit order24 coefficients on the owned finite middle interval required')
    E = mp.mpf(energy); a, r = mp.sqrt(3*mp.pi/2-4), mp.sqrt(2)
    S = sum((coefficients[j]/(2*E)**(j+1) for j in range(OLD_ORDER)), mp.mpc(0))
    dS = sum((coefficients[j]/(2*E)**(j+1) for j in range(OLD_ORDER, NEW_ORDER)), mp.mpc(0))
    db = stable_delta_bloch(S, dS, a); parallel = -E/a*db[2]
    return mp.matrix([-mass*db[0]+angular/r*db[1]+parallel,
                      parallel, mp.mpf(0), angular/(2*r)*db[1]])


def positive_refined_grid():
    """Independent48-point roots/weights at the caller's mpmath precision."""
    x, w = mp.gauss_quadrature(48, 'legendre')
    E = [start+4*(xi+1) for start in range(16, 160, 8) for xi in x]
    weights = [4*wi for _ in range(16, 160, 8) for wi in w]
    if any(v <= 0 for v in weights): raise ArithmeticError('positive refinement weights required')
    return E, weights


def integrate_correction(energies, weights, mass, angular, coefficients):
    if len(energies) != len(weights) or any(w <= 0 for w in weights):
        raise ValueError('matched positive correction quadrature required')
    kernels = [correction_kernel(E, mass, angular, coefficients) for E in energies]
    integral = mp.matrix([sum((w*K[i] for w, K in zip(weights, kernels)), mp.mpf(0)) for i in range(4)])
    return integral, np.array([[float(v) for v in K] for K in kernels])


def source16_control(energies, mass, angular, coefficients16):
    """Independent owned reference/vertex contraction, only to check convention."""
    E = np.asarray(energies, float)
    series = np.array([complex(sum((c/(2*mp.mpf(e))**(j+1) for j, c in enumerate(coefficients16)), mp.mpc(0))) for e in E])
    u = AXIAL**2*abs(series)**2
    bloch = np.stack((2*AXIAL*series.imag/(1+u), -2*AXIAL*series.real/(1+u), 1-2*u/(1+u)), axis=1)
    reference, _ = incoming_adiabatic_reference(E, float(mass), float(angular))
    reference_bloch = np.einsum('nij,vji->nv', reference, SIGMA).real
    db = bloch-reference_bloch; p = -E/AXIAL*db[:, 2]
    return np.stack((-float(mass)*db[:, 0]+float(angular)/RADIUS*db[:, 1]+p,
                     p, np.zeros(len(E)), float(angular)/(2*RADIUS)*db[:, 1]), axis=1)


def prepare_correction(root):
    owned = authenticated_group14(root); channel = owned['channel']; factor = mp.mpf(owned['factor'])
    arrays = {}; reports = {}; exact_integrals = {}
    for precision in (50, 70):
        with mp.workdps(precision):
            mass = mp.mpf(channel['compact_mass']); magnitude = mp.mpf(channel['angular_eigenvalue'])
            refined_E, refined_w = positive_refined_grid()
            for sign in SIGNS:
                angular = sign*magnitude; coefficients = _riccati_at_one(mass, angular, order=NEW_ORDER)
                c16 = _riccati_at_one(mass, angular, order=OLD_ORDER)
                prefix_error = max(abs(a-b) for a, b in zip(coefficients[:OLD_ORDER], c16))
                prefix = f'group14/mid24_{sign}/'; source = owned['arrays']
                original_E = [mp.mpf(v) for v in source[prefix+'energies']]
                original_w = [mp.mpf(v) for v in source[prefix+'weights']]
                choices = [('refined48', refined_E, refined_w)]
                if precision == 70: choices.append(('original24', original_E, original_w))
                for label, E, w in choices:
                    integral, kernels = integrate_correction(E, w, mass, angular, coefficients)
                    physical = factor*integral; key = f'{sign}/{label}/dps{precision}'
                    exact_integrals[key] = physical
                    arrays[key+'/energies'] = np.array([float(v) for v in E])
                    arrays[key+'/weights'] = np.array([float(v) for v in w])
                    arrays[key+'/correction_kernels'] = kernels
                    reports[key] = {'physical_correction': [float(v) for v in physical],
                                    'physical_correction_decimal': [mp.nstr(v, precision) for v in physical],
                                    'positive_weights': True, 'rows': len(E),
                                    'integrated_dE_measure_residual': float(abs(sum(w)-144)),
                                    'coefficient_prefix_residual': float(prefix_error)}
                if precision == 70:
                    control = source16_control(source[prefix+'energies'], mass, angular, c16)
                    reports[str(sign)] = {'order16_vacuum_kernel_vs_archive': float(np.max(abs(control-source[prefix+'vacuum_kernels']))),
                        'source16_decomposition_residual': float(np.max(abs(source[prefix+'kernels']-source[prefix+'vacuum_kernels']-source[prefix+'thermal_kernels']))),
                        'coefficients24_decimal': [[mp.nstr(mp.re(c), 70), mp.nstr(mp.im(c), 70)] for c in coefficients]}
                    for quantity in ('energies', 'weights', 'kernels', 'vacuum_kernels', 'thermal_kernels'):
                        arrays[f'{sign}/archived16/'+quantity] = source[prefix+quantity].copy()
            print(f'group14 order24 correction: precision{precision}, both angular signs; local source algebra only', flush=True)
    with mp.workdps(70):
        grid_indicator = [sum(abs(exact_integrals[f'{s}/refined48/dps70'][i]-exact_integrals[f'{s}/original24/dps70'][i]) for s in SIGNS) for i in range(4)]
        precision_indicator = [sum(abs(exact_integrals[f'{s}/refined48/dps70'][i]-exact_integrals[f'{s}/refined48/dps50'][i]) for s in SIGNS) for i in range(4)]
    return arrays, {'group': 14, 'old_order': OLD_ORDER, 'new_order': NEW_ORDER,
        'group_factor_per_sign': owned['factor'], 'interval': owned['interval'],
        'grid_recipe_residual': owned['grid_recipe_residual'], 'reports': reports,
        'sum_absolute_sign_quadrature_indicator': [float(v) for v in grid_indicator],
        'sum_absolute_sign_precision_indicator': [float(v) for v in precision_indicator],
        'thermal_scattering_stress_error_upper': owned['thermal_scattering_stress_error_upper'],
        'input_payloads': owned['input_payloads']}
