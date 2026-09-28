"""One-sign directed certificates for the remaining retained high bands."""
from functools import lru_cache
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from .nsc_incoming_centered_order24 import CenteredOrder24Geometry, trimmed_defect_jets
from .nsc_incoming_defect_taylor_bound import centered_coefficient_enclosure
from .nsc_incoming_middle_bound import (
    _parameters, authenticated_middle_inputs, thermal_middle_difference_bound, constraint_action_bound,
)
from .nsc_incoming_projector_energy_bound import projected_energy_bound, local_cross_product_coefficients
from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _geometry_jets, _up_float

GROUPS = tuple(g for g in range(1, 33) if g not in (12, 14, 22, 32))


@lru_cache(maxsize=1)
def sign_symmetry_identity():
    """Exact base/induction/defect parity, with arbitrary complex coefficients."""
    L, m, A = sp.symbols('L m A', real=True)
    D, U, S = sp.symbols('D U S')
    base = sp.expand((-L+sp.I*m)+sp.conjugate(L+sp.I*m))
    induction = []
    for n in range(1, 24):
        for j in range(1, n):
            if (-1)**j*(-1)**(n-j) != (-1)**n: raise ArithmeticError('quadratic parity')
        transformed = sp.I*(-1)**n*sp.conjugate(D)-A*sp.conjugate(U)*(-1)**n*sp.conjugate(S)
        expected = (-1)**(n+1)*sp.conjugate(sp.I*D+A*U*S)
        induction.append(str(sp.expand(transformed-expected)))
    defect = []
    for p in range(24, 49):
        transformed = A*sp.conjugate(U)*(-1)**p*sp.conjugate(S)
        expected = (-1)**(p+1)*sp.conjugate(-A*U*S)
        if p == 24:
            transformed -= sp.I*(-1)**24*sp.conjugate(D)
            expected += (-1)**25*sp.conjugate(-sp.I*D)
        defect.append(str(sp.expand(transformed-expected)))
    x, y = sp.symbols('x y', real=True)
    norm = str(sp.expand(((-x)**2+y*y)-(x*x+y*y)))
    if base != 0 or set(induction+defect+[norm]) != {'0'}:
        raise ArithmeticError('angular-sign induction failed')
    return {'base_residual': str(base), 'recurrence_residuals': induction,
            'defect_residuals': defect, 'norm_squared_residual': norm,
            'coefficient_identity': 'c_n(-lambda)=(-1)^n conjugate(c_n(lambda))',
            'defect_identity': 'f_p(-lambda)=(-1)^(p+1) conjugate(f_p(lambda))',
            'physical_occupation_symmetry_assumed': False}


def interval_sign_control():
    """One actual small-cell order24 comparison, not another channel bound."""
    with _precision(40):
        geometry = _geometry_jets(mp.iv.mpf(['1.2', '1.201']), 28)
        plus = trimmed_defect_jets(*geometry, mp.iv.pi/2, mp.iv.sqrt(5))
        minus = trimmed_defect_jets(*geometry, mp.iv.pi/2, -mp.iv.sqrt(5))
        maximum = mp.mpf(0)
        for power, pos, neg in zip(range(24, 49), plus, minus):
            for a, b in zip(pos, neg):
                expected = (-1)**(power+1)*mp.iv.mpc(a.real, -a.imag)
                for actual, wanted in ((b.real, expected.real), (b.imag, expected.imag)):
                    maximum = max(maximum, abs(_lo(actual)-_lo(wanted)), abs(_hi(actual)-_hi(wanted)))
        if maximum != 0: raise ArithmeticError('interval sign transformation differs')
        return {'order': 24, 'depth': 4, 'interval': [1.2, 1.201], 'maximum_endpoint_difference': 0.,
                'all_25_defect_coefficients_and_five_coordinate_jets_checked': True}


def authenticated_batch(root):
    """Actual selected LOW cells and middle intervals; group13 is unsplit."""
    root = Path(root); owned = authenticated_middle_inputs(root)
    finite = json.loads((root/'results/development/nsc-incoming-spectral-source.json').read_text())
    retained = json.loads((root/'results/development/nsc-pg-retained-covariance.json').read_text())
    windows = {}
    with np.load(root/finite['payload']['path'], allow_pickle=False) as src, np.load(root/retained['payload']['path'], allow_pickle=False) as original:
        for group in GROUPS:
            channel = owned['channels'][group]; middle = owned['middle_intervals'][group-1]
            if middle['group'] != group: raise ValueError('middle group order changed')
            right = 320. if group in (10, 11, 31) else 160.
            if (middle['left'], middle['right']) != ((16., 160.) if group == 13 else (40., right)):
                raise ValueError('actual middle endpoints changed')
            rows = [{'kind': 'middle', 'interval': [middle['left'], middle['right']],
                     'panels': middle['signed_selected_panels']}]
            if group != 13:
                low = []
                signs = (1,) if not channel['angular_eigenvalue'] else (1, -1)
                for sign in signs:
                    names = [p['panel'] for p in finite['groups'][str(group)]['signed_inventory'][str(sign)]['pieces']
                             if p['panel'].startswith(('low/', 'low_ref/'))]
                    if len(names) != 1: raise ValueError('one actual LOW panel per sign required')
                    name = names[0]; meta = json.loads(original[name+'/metadata_json'].tobytes())
                    if meta['channel'] != channel or meta['sign'] != sign or (32., 40.) not in list(zip(meta['edges'], meta['edges'][1:])):
                        raise ValueError('actual signed LOW[32,40] cell required')
                    E, w = src[name+'/energies'], src[name+'/weights']
                    if not np.array_equal(E, original[name+'/energy'].ravel()) or not np.array_equal(w, original[name+'/weight'].ravel()):
                        raise ValueError('LOW source nodes differ from authenticated producer')
                    mask = (E > 32)&(E < 40); residual = abs(float(w[mask].sum())-8.)
                    # Metadata's base point count excludes some archived
                    # cell-specific refinements; the authenticated arrays own
                    # the actual selected count and positive measure.
                    if not mask.any() or np.any(w[mask] <= 0) or residual > 3e-11:
                        raise ValueError('positive complete actual LOW cell required')
                    low.append({'sign': sign, 'panel': name, 'rows': int(mask.sum()), 'measure_residual': residual})
                rows.insert(0, {'kind': 'low', 'interval': [32., 40.], 'panels': low})
            windows[str(group)] = rows
    if sum(map(len, windows.values())) != 55: raise ArithmeticError('28 groups and55 disjoint windows required')
    return {'channels': owned['channels'], 'config': owned['config'], 'windows': windows,
            'input_payloads': owned['input_payloads']}


class CachedGeometry:
    """Common128-cell geometry cached once per CPU worker; one sign per group."""
    def __init__(self, config):
        geometry = CenteredOrder24Geometry(config)
        self.weight, self.horizon = geometry.weight, geometry.horizon
        with _precision(40):
            self.cells = [(rho, center, _geometry_jets(rho, 28), _geometry_jets(center, 28))
                          for rho, center in geometry.cells]

    def bound(self, channel, *, progress=None):
        group = channel['index']
        if group not in GROUPS: raise ValueError('completed or unowned channel cannot be rerun')
        sign_symmetry_identity()
        with _precision(40):
            _, _, mass, angular, _ = _parameters(channel)
            totals = [mp.iv.mpf(0) for _ in range(25)]; natural = [mp.iv.mpf(0) for _ in range(25)]
            for index, (rho, center, geometry, point_geometry) in enumerate(self.cells):
                whole = trimmed_defect_jets(*geometry, mass, angular)
                point = trimmed_defect_jets(*point_geometry, mass, angular)
                for j, (w, p) in enumerate(zip(whole, point)):
                    bounded = centered_coefficient_enclosure(w, p, rho-center, 4)
                    totals[j] += self.weight*mp.iv.mpf(_hi(abs(bounded)))/128
                    natural[j] += self.weight*mp.iv.mpf(_hi(abs(w[0])))/128
                if progress and (index+1)%64 == 0: progress(index+1)
            bounds, controls = [_up_float(v) for v in totals], [_up_float(v) for v in natural]
            signs = (1,) if not channel['angular_eigenvalue'] else (1, -1)
            return {'group': group, 'lower': 320. if group in (10, 11, 31) else 160.,
                    'physical_Riccati_order': 24, 'intervals': 128, 'centered_remainder_depth': 4,
                    'first_inverse_power': 24, 'last_inverse_power': 48, 'directed_precision': 40,
                    'per_sign': [{'sign': s, 'coefficient_integral_upper': bounds[:],
                                 'natural_coefficient_integral_upper': controls[:],
                                 'norm_bound_provenance': 'computed positive angular sign' if s == 1 else 'exact defect sign-conjugation symmetry'} for s in signs],
                    'angular_norm_recurrences': 1, 'physical_source_signs_identified': False,
                    'endpoint_weight_upper': _up_float(self.weight),
                    'horizon_enclosure': [mp.nstr(_lo(self.horizon),55), mp.nstr(_hi(self.horizon),55)],
                    'initial_projector_difference': 0.,
                    'initial_scope': 'same declared affine-horizon vacuum; explicit order24 only',
                    'finite_offset_archive_accuracy_certified': False,
                    'explicit_matching_source_corrections_required': True}


def contract_channel(channel, config, windows, radial, cross):
    """Use each disjoint window once; the union is an alternative aggregate."""
    group = channel['index']; result = {}
    specifications = [(r['kind'], *r['interval']) for r in windows]
    if group != 13: specifications.append(('combined', windows[0]['interval'][0], windows[-1]['interval'][1]))
    for name, left, right in specifications:
        view = {**radial, 'lower': right}
        energy = projected_energy_bound(channel, view, cross, left, right)
        thermal = thermal_middle_difference_bound(channel, config, left, right)
        action = constraint_action_bound(thermal)
        with _precision(40):
            total = [_up_float(mp.iv.mpf(energy['lapse_action_error_upper'])+mp.iv.mpf(action[0])), action[1]]
        result[name] = {'interval': [left, right], 'vacuum_energy': energy,
                        'thermal_stress_upper': thermal, 'constraint_action_error_upper': total,
                        'within_component_tolerance': max(total) <= 3e-11,
                        'original_endpoint_tag': radial['lower'], 'consumer_endpoint_tag': right,
                        'radial_coefficients_changed': False}
    aggregate = result['middle' if group == 13 else 'combined']['constraint_action_error_upper']
    excess = max([0.]+[v/n-1 for r in radial['per_sign'] for v,n in zip(r['coefficient_integral_upper'],r['natural_coefficient_integral_upper']) if n])
    additivity = 0.
    if group != 13:
        additivity = max(abs(result['combined']['vacuum_energy'][k]-result['low']['vacuum_energy'][k]-result['middle']['vacuum_energy'][k])
                         for k in ('density_error_upper', 'lapse_action_error_upper'))
    if excess > 3e-15 or additivity > 1e-24: raise ArithmeticError('intersection or interval additivity failed')
    return {'group': group, 'windows': result, 'aggregate_constraint_action_error_upper': aggregate,
            'status': 'PASS conditional on matching order24 source corrections' if max(aggregate) <= 3e-11 else 'OPEN: valid bound exceeds component tolerance',
            'residuals': {'centered_excess_over_natural': excess, 'interval_additivity': additivity},
            'stationarity_tolerance': 3e-11, 'source_windows_count': len(windows),
            'combined_is_alternative_not_additional': True}
