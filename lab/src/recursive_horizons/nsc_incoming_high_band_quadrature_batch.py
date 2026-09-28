"""Direct interval source integrals replacing complete archived high regions.

The frozen ``low_vacuum`` kernel means order24-minus-ad4 DIRECT vacuum. It
is applied to the authenticated LOW-plus-middle union, not a correction.
"""
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path

import mpmath as mp
import numpy as np

from .nsc_incoming_middle_bound import authenticated_middle_inputs, _parameters
from .nsc_incoming_source_assembly import read_panels
from .nsc_incoming_state_moments import incoming_group_factor
from .nsc_incoming_source_quadrature_bound import (
    certify_source_window, gauss_rule, replay_window, pack, unpack, float_enclosure,
)
from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _range


GROUPS = (6,)+tuple(g for g in range(1, 33) if g != 6)


@lru_cache(maxsize=1)
def authenticated_high_inventory(root):
    """Authenticate all old endpoints once; retain only small source summaries."""
    root = Path(root)
    middle = authenticated_middle_inputs(root)
    panels, _ = read_panels(root)
    finite = json.loads((root/'results/development/nsc-incoming-spectral-source.json').read_text())
    retained = json.loads((root/'results/development/nsc-pg-retained-covariance.json').read_text())
    result = {}
    with np.load(root/retained['payload']['path'], allow_pickle=False) as originals:
        for channel in middle['channels'][1:]:
            group = channel['index']; source = finite['groups'][str(group)]
            interval = middle['middle_intervals'][group-1]
            signs = (1,) if channel['angular_eigenvalue'] == 0 else (1, -1)
            left = 16 if group in (13, 14) else 32
            right = int(interval['right'])
            if interval['left'] != (16 if group in (13, 14) else 40):
                raise ValueError('owned middle interval changed')
            factor = incoming_group_factor(channel, signs)
            if factor != source['factor']:
                raise ValueError('signed multiplicity or original source factor changed')
            old = np.zeros(4); descriptors = []; by_sign = {}
            for sign in signs:
                old_sign = np.zeros(4)
                selected = [p for p in interval['signed_selected_panels'] if p['sign'] == sign]
                if len(selected) != 1:
                    raise ValueError('one actual middle panel per signed family required')
                if group not in (13, 14):
                    pieces = source['signed_inventory'][str(sign)]['pieces']
                    lows = [p['panel'] for p in pieces if p['panel'].startswith(('low/', 'low_ref/'))]
                    if len(lows) != 1:
                        raise ValueError('one actual LOW source panel required')
                    name = lows[0]
                    metadata = json.loads(originals[name+'/metadata_json'].tobytes())
                    if (metadata['channel'] != channel or metadata.get('sign', metadata.get('angular_sign')) != sign
                            or (32, 40) not in list(zip(metadata['edges'], metadata['edges'][1:]))):
                        raise ValueError('LOW source labels or actual [32,40] cell changed')
                    E, w = panels[name+'/energies'], panels[name+'/weights']
                    if not np.array_equal(E, originals[name+'/energy'].ravel()) or not np.array_equal(w, originals[name+'/weight'].ravel()):
                        raise ValueError('old LOW source nodes differ from their producer')
                    mask = (E > 32)&(E < 40)
                    if not np.any(mask) or np.any(w[mask] <= 0) or abs(float(w[mask].sum())-8) > 3e-11:
                        raise ValueError('complete positive LOW cell required')
                    old_sign += factor*np.einsum('n,nv->v', w[mask], panels[name+'/kernels'][mask])
                    descriptors.append({'sign': sign, 'panel': name, 'interval': [32, 40], 'rows': int(mask.sum())})
                selection = selected[0]; name = selection['panel']
                old_sign += factor*np.einsum('n,nv->v', panels[name+'/weights'], panels[name+'/kernels'])
                descriptors.append({'sign': sign, 'panel': name, 'interval': [int(selection['left']), right],
                                    'rows': selection['rows']})
                old += old_sign; by_sign[str(sign)] = old_sign.tolist()
            result[group] = {'group': group, 'channel': channel, 'interval': [left, right],
                             'actual_signs': list(signs), 'factor_per_sign': factor,
                             'selected_original_panels': descriptors,
                             'original_high_source': old.tolist(), 'original_high_source_per_sign': by_sign,
                             'input_payloads': middle['input_payloads'],
                             'original_source_scope': 'original finite source arrays, including their high-band thermal insertion; before later numerical refinements'}
    if set(result) != set(range(1, 33)):
        raise ValueError('all32 retained non-LLL groups required')
    return result


@lru_cache(maxsize=1)
def cached_rule():
    return gauss_rule(points=48, precision=80)


def prepare_direct_group(root, group, progress=None):
    """One fixed-resolution new quadrature proof of the complete direct source."""
    if group not in GROUPS:
        raise ValueError('one retained non-LLL group required')
    owned = authenticated_high_inventory(str(Path(root)))[group]
    left, right = owned['interval']
    direct = certify_source_window(owned['channel'], 'low_vacuum', left, right,
                                   rule=cached_rule(), points=48, precision=80, progress=progress)
    return {'purpose': 'direct full order24-minus-ad4 vacuum; replaces complete original high-region source',
            'owned_union': owned, 'direct': direct}


def replay_direct_group(root, prepared):
    """Replay saved directed samples only; no quadrature/source generation."""
    group = prepared['owned_union']['group']
    owned = authenticated_high_inventory(str(Path(root)))[group]
    direct = prepared['direct']
    if prepared['owned_union'] != owned or direct['channel'] != owned['channel']:
        raise ValueError('original high union or channel provenance changed')
    if direct['group'] != group or direct['kind'] != 'low_vacuum' or direct['interval'] != owned['interval'] or direct['source_order'] != 24:
        raise ValueError('complete direct order24 source required, not a correction')
    if direct['points'] != 48 or direct['precision'] != 80:
        raise ValueError('single fixed48-node80-digit batch resolution required')
    left, right = owned['interval']
    expected = {(sign, start, start+8) for sign in owned['actual_signs'] for start in range(left, right, 8)}
    actual = [(row['sign'], row['left'], row['right']) for row in direct['cells']]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError('high-region cells missing, duplicated or overlapping')
    numerical = replay_window(direct)
    old = np.array(owned['original_high_source']); new = np.array(numerical['certified_approximant'])
    return {'group': group, 'interval': owned['interval'], 'purpose': prepared['purpose'],
            'source_order': 24, 'old_baseline_quadrature_carried': False,
            'original_high_source': old.tolist(), 'direct_high_vacuum_source': new.tolist(),
            'replacement_delta_from_original_high_source': (new-old).tolist(),
            'selected_original_panels': owned['selected_original_panels'],
            'actual_signs': owned['actual_signs'], **numerical,
            'thermal_policy': 'old high thermal insertion removed with original high source; physical thermal remains separately bounded, never declared zero'}


def aggregate_groups(groups):
    """Sum exact binary interval endpoints; include final float rounding once."""
    if set(groups) != {str(g) for g in GROUPS}:
        raise ValueError('exactly32 complete direct source receipts required')
    with _precision(80):
        source = [sum((unpack(groups[str(g)]['source_integral_interval'][i]) for g in GROUPS), mp.iv.mpf(0)) for i in range(4)]
        value = [float((_lo(v)+_hi(v))/2) for v in source]
        errors = [mp.iv.mpf(_hi(abs(v-mp.iv.mpf(x)))) for v, x in zip(source, value)]
        a, r = mp.iv.sqrt(3*mp.iv.pi/2-4), mp.iv.sqrt(2)
        lapse = 4*mp.iv.pi*a*r*r*errors[0]
        old = [sum((mp.iv.mpf(groups[str(g)]['original_high_source'][i]) for g in GROUPS), mp.iv.mpf(0)) for i in range(4)]
        old_value = [float((_lo(v)+_hi(v))/2) for v in old]
        return {'direct_high_vacuum_source': value,
                'source_integral_interval': [pack(v) for v in source],
                'source_integral_float_enclosure': [float_enclosure(v) for v in source],
                'stress_error_upper': [float_enclosure(v)[1] for v in errors],
                'lapse_quadrature_arithmetic_error_upper': float_enclosure(lapse)[1],
                'original_high_source': old_value,
                'replacement_delta_from_original_high_source': [n-o for n, o in zip(value, old_value)],
                'stationarity_tolerance': 3e-11,
                'numerical_budget_pass': _hi(lapse) <= mp.mpf('3e-11'),
                'root_composition_instruction': 'replace original finite source high regions once; do not add these values atop old high source or earlier corrections',
                'thermal_policy': 'separate nonzero physical thermal bound required over all replaced regions'}
