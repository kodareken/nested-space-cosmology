"""Finite middle-band bounds from already certified radial defect integrals."""
from hashlib import sha256
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

from .nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi, _range, _up_float


STATIONARITY_TOLERANCE = 3e-11
INPUT_RECORDS = (
    'results/development/nsc-incoming-vacuum-tail-bound.json',
    'results/development/nsc-incoming-spectral-source.json',
    'results/development/nsc-pg-retained-covariance.json',
    'results/development/nsc-pg-group13-covariance.json',
    'results/development/nsc-mode-resolved-cauchy-state.json',
    'results/development/nsc-compact-matched-restart.json',
)


def _payload(root, spec):
    raw = (root/spec['path']).read_bytes()
    if sha256(raw).hexdigest() != spec['sha256']:
        raise ValueError('authenticated middle input changed: '+spec['path'])
    if spec['path'].endswith('.json'): return json.loads(raw)
    with np.load(root/spec['path'], allow_pickle=False) as values:
        return {key: values[key].copy() for key in values.files}


def validate_middle_interval(channel, sign, metadata, source_endpoint):
    """Require the actual signed panel labels and its owned interval."""
    group = channel['index']; expected_left = 16. if group in (13, 14) else 40.
    expected_right = 320. if group in (10, 11, 12, 31, 32) else 160.
    if (metadata['angular_sign'] != sign or metadata['channel'] != channel
            or metadata['left'] != expected_left or metadata['right'] != expected_right
            or source_endpoint != expected_right):
        raise ValueError('actual signed middle-panel labels or endpoints differ from the owned interval')
    return expected_left, expected_right


def authenticated_middle_inputs(root):
    """Read selected middle panels, not a guessed new spectral partition."""
    root = Path(root); records = [json.loads((root/p).read_text()) for p in INPUT_RECORDS]
    bound, finite, retained_record, g13record, inventory, restart = records
    for record in (bound, finite, retained_record, g13record):
        for field in ('source_hashes', 'input_hashes', 'numerical_owner_and_input_hashes'):
            for path, expected in record.get(field, {}).items():
                if sha256((root/path).read_bytes()).hexdigest() != expected:
                    raise ValueError('middle input source or record changed: '+path)
    artifact = _payload(root, bound['payload'])
    if artifact['groups'] != bound['groups']:
        raise ValueError('radial defect record and coefficient artifact differ')
    panels = _payload(root, finite['payload']); retained = _payload(root, retained_record['payload'])
    g13 = _payload(root, g13record['payload'])
    retained_meta = json.loads(retained['metadata_json'].tobytes())
    g14 = _payload(root, retained_meta['group14_input'])
    g13meta = json.loads(g13['metadata_json'].tobytes())
    channels = inventory['channels']; provenance = []; coverage = 0.
    if [c['index'] for c in channels] != list(range(33)):
        raise ValueError('the full owned33-group inventory is required')
    for channel in channels[1:]:
        group = channel['index']; family = finite['groups'][str(group)]; signed = []
        signs = (1,) if channel['angular_eigenvalue'] == 0 else (1, -1)
        if set(family['signed_inventory']) != {str(s) for s in signs}:
            raise ValueError('actual angular-sign inventory changed')
        for sign in signs:
            selected = family['signed_inventory'][str(sign)]
            names = [p['panel'] for p in selected['pieces']
                     if p['panel'].startswith(('mid/', 'mid_ref/', 'group13/mid', 'group14/mid'))]
            if len(names) != 1: raise ValueError('exactly one selected middle panel per signed family required')
            name = names[0]
            if group == 13:
                left, right = g13record['energy_account']['fast_mode_window']
                metadata = {'left': left, 'right': right, 'channel': g13meta['channel'], 'angular_sign': 1,
                            'sources': g13meta['source_hashes']}
                E, w = g13['fast_energy_24'], g13['fast_weight_24']
            else:
                data, prefix = (g14, name.split('/')[1]) if group == 14 else (retained, name)
                metadata = json.loads(data[prefix+'/metadata_json'].tobytes())
                E, w = data[prefix+'/energies'], data[prefix+'/weights']
            for path, expected in metadata['sources'].items():
                if sha256((root/path).read_bytes()).hexdigest() != expected:
                    raise ValueError('selected middle producer changed: '+path)
            left, right = validate_middle_interval(channel, sign, metadata, selected['upper_energy_endpoint'])
            if (not np.array_equal(E, panels[name+'/energies']) or not np.array_equal(w, panels[name+'/weights'])
                    or np.any(E <= left) or np.any(E >= right) or np.any(w <= 0)):
                raise ValueError('selected source nodes differ from authenticated middle producer')
            residual = abs(float(w.sum())-(right-left)); coverage = max(coverage, residual)
            if residual > 3e-11: raise ValueError('middle quadrature measure does not match its metadata endpoints')
            signed.append({'sign': sign, 'panel': name, 'left': left, 'right': right,
                           'rows': len(E), 'endpoint_measure_residual': residual})
        provenance.append({'group': group, 'left': signed[0]['left'], 'right': signed[0]['right'],
                           'signed_selected_panels': signed})
    return {'channels': channels, 'config': restart['scattering_provenance']['config'],
            'radial_groups': bound['groups'], 'middle_intervals': provenance,
            'endpoint_measure_residual': coverage,
            'input_payloads': [bound['payload'], finite['payload'], retained_record['payload'],
                               g13record['payload'], retained_meta['group14_input']]}


def _parameters(channel):
    mass_definition = channel['compact_level']*mp.iv.pi/2
    angular_definition = mp.iv.sqrt(channel['angular_level']*(channel['angular_level']+4))
    mass = _range(min(_lo(mass_definition), mp.mpf(float(channel['compact_mass']))),
                  max(_hi(mass_definition), mp.mpf(float(channel['compact_mass']))))
    ell = _range(min(_lo(angular_definition), mp.mpf(float(channel['angular_eigenvalue']))),
                 max(_hi(angular_definition), mp.mpf(float(channel['angular_eigenvalue']))))
    a, r = mp.iv.sqrt(3*mp.iv.pi/2-4), mp.iv.sqrt(2)
    factor = channel['copy_count']*channel['degeneracy']/(4*mp.iv.pi**2*r*r*a)
    return a, r, mass, ell, factor


def finite_power_primitive(power, moment, left, right):
    """Positive directed integral of E**(moment-power), with power>=16."""
    if power not in range(16, 33) or moment not in (0, 1) or not 0 < left < right:
        raise ValueError('owned positive finite interval and defect powers16..32 required')
    p = power-moment-1
    return (mp.iv.mpf(left)**(-p)-mp.iv.mpf(right)**(-p))/p


def vacuum_middle_bound(channel, radial, left, right, *, precision=40):
    """No recurrence: contract saved positive radial integrals only."""
    if radial['group'] != channel['index'] or radial['lower'] != right:
        raise ValueError('saved defect group and original tail endpoint must match this middle band')
    expected = (1,) if channel['angular_eigenvalue'] == 0 else (1, -1)
    if tuple(p['sign'] for p in radial['per_sign']) != expected:
        raise ValueError('each actual signed radial defect bound required once')
    with _precision(precision):
        a, r, mass, ell, factor = _parameters(channel); factor /= len(expected)
        M = mp.iv.sqrt(mass*mass+(ell/r)**2); total = [mp.iv.mpf(0) for _ in range(4)]
        for row in radial['per_sign']:
            values = row['I16_to_I32_upper']
            if len(values) != 17 or any(not math.isfinite(v) or v < 0 for v in values):
                raise ValueError('positive finite saved radial coefficient bounds required')
            for power, bound in zip(range(16, 33), values):
                j0 = finite_power_primitive(power, 0, left, right)
                j1 = finite_power_primitive(power, 1, left, right)
                common = factor*mp.iv.mpf(bound)/mp.iv.mpf(2)**power
                total[0] += 2*common*(j1/a+M*j0); total[1] += 2*common*j1/a
                total[3] += common*ell*j0/r
        return [_up_float(v) if i != 2 else 0. for i, v in enumerate(total)]


def thermal_middle_difference_bound(channel, config, left, right, *, precision=40):
    """Triangle bound2*b(E) for actual and approximate thermal insertions."""
    if not 0 < left < right: raise ValueError('positive finite middle interval required')
    with _precision(precision):
        a, r, mass, ell, factor = _parameters(channel)
        kappa, omega = mp.iv.mpf(config['surface_gravity']), mp.iv.mpf(config['omega'])
        if _lo(kappa) <= 0 or _lo(omega) <= 0: raise ValueError('positive owned state parameters required')
        i0 = mp.iv.mpf(0); i1 = mp.iv.mpf(0)
        for coefficient, alpha in ((2, mp.iv.pi/kappa), (1, 2*mp.iv.pi/(omega*kappa))):
            base = coefficient*mp.iv.exp(-alpha*left); remaining = mp.iv.exp(-alpha*(right-left))
            i0 += base*(1-remaining)/alpha
            i1 += base*((left+1/alpha)-remaining*(right+1/alpha))/alpha
        M = mp.iv.sqrt(mass*mass+(ell/r)**2)
        # Two departures from respective vacuum projectors, both bounded
        # using the unchanged source law and current-normalized compression.
        total = [2*factor*(i1/a+M*i0), 2*factor*i1/a, 2*factor*i1/a, factor*ell*i0/r]
        return [_up_float(v) for v in total]


def positive_sum(vectors, *, precision=40):
    with _precision(precision):
        sums = [sum((mp.iv.mpf(v[i]) for v in vectors), mp.iv.mpf(0)) for i in range(4)]
        return [_up_float(v) if _hi(v) != 0 else 0. for v in sums]


def constraint_action_bound(stress, *, precision=40):
    """Absolute N,beta conversion of the owned same-surface action map."""
    if len(stress) != 4 or any(not math.isfinite(v) or v < 0 for v in stress):
        raise ValueError('four nonnegative finite stress bounds required')
    with _precision(precision):
        a = mp.iv.sqrt(3*mp.iv.pi/2-4); r2 = mp.iv.mpf(2)
        values = [4*mp.iv.pi*a*r2*stress[0], 4*mp.iv.pi*a*a*r2*stress[2]]
        return [_up_float(v) if _hi(v) else 0. for v in values]
