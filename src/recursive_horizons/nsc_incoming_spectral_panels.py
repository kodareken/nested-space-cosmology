"""Apply the incoming source kernels to authenticated retained spectral panels.

This produces finite state-minus-reference contributions, not a complete
stress tensor. Base/refined quadratures and the eight-node subgap control
remain distinct. Neither packet covariance nor a neck tensor is input.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from .nsc_incoming_state_moments import incoming_state_moments
from .nsc_incoming_high_energy_source import incoming_high_energy_source
from .nsc_paired_horizon_preparation import source_covariance


RECORDS = ('results/development/nsc-pg-spectral-mode-recovery.json',
           'results/development/nsc-pg-retained-covariance.json',
           'results/development/nsc-pg-group13-covariance.json',
           'results/development/nsc-pg-massive-mode-resolution.json',
           'results/development/nsc-mode-resolved-cauchy-state.json',
           'results/development/nsc-compact-matched-restart.json')
NUMERICAL_OWNERS = (
    'src/recursive_horizons/nsc_incoming_spectral_panels.py',
    'src/recursive_horizons/nsc_incoming_state_moments.py',
    'src/recursive_horizons/nsc_incoming_high_energy_source.py',
    'src/recursive_horizons/nsc_pg_archived_high_energy_modes.py',
    'src/recursive_horizons/nsc_pg_high_energy.py',
    'src/recursive_horizons/nsc_pg_fast_packets.py',
    'src/recursive_horizons/nsc_common_ks_trace.py',
    'src/recursive_horizons/nsc_compact_ctp_neck.py',
    'src/recursive_horizons/nsc_paired_horizon_preparation.py',
    'src/recursive_horizons/nsc_lorentzian.py',
    'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
)


def sha(root, path): return hashlib.sha256((Path(root)/path).read_bytes()).hexdigest()
def record(root, path): return json.loads((Path(root)/path).read_text())
def load(path):
    with np.load(path, allow_pickle=False) as a:
        return {k: a[k].copy() for k in a.files}
def authenticated(root, spec):
    if sha(root, spec['path']) != spec['sha256']:
        raise ValueError('input artifact changed: '+spec['path'])
    return load(Path(root)/spec['path'])


def compute_panels(root, progress=None):
    """Evaluate all real-panel kernels plus the inherited subgap base nodes.

    Only the new observable is evaluated. Middle modes are exposed from their
    ORIGINAL producer recipe, avoiding inversion/normalization roundoff in a
    small C-P. No mode, horizon, scattering or packet integrator is called.
    """
    root = Path(root)
    signature = {p: sha(root, p) for p in (*NUMERICAL_OWNERS, *RECORDS)}
    receipts = {p: record(root, p) for p in RECORDS}
    recovery = authenticated(root, receipts[RECORDS[0]]['payload'])
    retained = authenticated(root, receipts[RECORDS[1]]['payload'])
    g13 = authenticated(root, receipts[RECORDS[2]]['payload'])
    mode = authenticated(root, receipts[RECORDS[3]]['payload'])
    retained_meta = json.loads(retained['metadata_json'].tobytes())
    g14 = authenticated(root, retained_meta['group14_input'])
    table = receipts[RECORDS[4]]['channels']
    config = receipts[RECORDS[5]]['scattering_provenance']['config']
    g13meta = json.loads(g13['metadata_json'].tobytes())
    g13middle = {'sources': g13meta['source_hashes'], 'channel': g13meta['channel'],
                 'angular_sign': 1, 'left': 16., 'right': 160., 'config': g13meta['config']}
    arrays = {}; descriptions = {}

    def add(name, group, sign, E, weights, C, P, *, phi=None, middle=None, provenance):
        ch = table[group]; m, ell = ch['compact_mass'], sign*ch['angular_eigenvalue']
        E, weights = np.asarray(E).ravel(), np.asarray(weights).ravel()
        if weights.shape != E.shape or np.any(weights <= 0): raise ValueError('owned positive dE weights required')
        if middle is None:
            result = incoming_state_moments(E, phi, C, P, m, ell)
            diagnostics = result['diagnostics']
            per_energy = result['per_energy_diagnostics']
            method = 'resolved physical columns minus existing ad4 reference'
        else:
            result = incoming_high_energy_source(E, middle, C, repo_root=root)
            diagnostics = {k: float(np.max(v)) for k, v in result['diagnostics'].items()}
            per_energy = result['diagnostics']
            method = 'same archived order16 boundary recipe; stable vacuum/thermal subtraction'
        for key, value in {'energies': E, 'weights': weights, 'kernels': result['kernels'],
                           'arithmetic_indicator': result['roundoff_indicator']}.items():
            arrays[name+'/'+key] = value
        for key in ('raw_kernels', 'raw_roundoff_indicator'):
            if key in result: arrays[name+'/'+key] = result[key]
        arrays[name+'/diagnostics_json'] = np.frombuffer(json.dumps(diagnostics, sort_keys=True).encode(), np.uint8)
        for key, value in per_energy.items(): arrays[name+'/diagnostics/'+key] = np.asarray(value)
        if middle is not None:
            for key in ('vacuum_kernels', 'thermal_kernels', 'direct_kernels'):
                arrays[name+'/'+key] = result[key]
        descriptions[name] = {'group': group, 'sign': sign, 'mass': m, 'angular': ell,
            'rows': len(E), 'energy_range': [float(E.min()), float(E.max())],
            'method': method, 'provenance': provenance,
            'whole_spectrum': False, 'local_restoration_included': False}

    names = sorted(k[:-len('/energies')] for k in recovery if k.endswith('/energies'))
    for i, name in enumerate(names, 1):
        get = lambda k: recovery[name+'/'+k]
        group, sign = int(get('group').item()), int(get('sign').item())
        middle = None
        if name.startswith(('mid/', 'mid_ref/')):
            middle = json.loads(retained[name+'/metadata_json'].tobytes())
        elif name == 'group13/mid': middle = g13middle
        elif name.startswith('group14/mid'):
            middle = json.loads(g14[name.split('/')[1]+'/metadata_json'].tobytes())
        if middle is None and not get('accepted').all(): raise ValueError('unresolved low source field '+name)
        add(name, group, sign, get('energies'), get('weights'), get('source'), get('projector'),
            phi=get('mode_at_one') if middle is None else None, middle=middle,
            provenance='authenticated retained real-panel fields and original middle producer')
        if progress and (i % 12 == 0 or i == len(names)): progress(i, len(names))

    # Existing comparison grids not needed by the inverse are still available
    # from the same analytic middle recipe; no new spectral nodes are chosen.
    def middle_comparison(name, group, sign, E, weights, middle):
        m = table[group]['compact_mass']
        C = np.array([source_covariance(e, config['surface_gravity'], config['omega'], m) for e in E])
        P = np.broadcast_to(np.eye(3), C.shape).copy()
        add(name, group, sign, E, weights, C, P, middle=middle,
            provenance='archived comparison-only middle quadrature, same source law')
    middle_comparison('group13/mid16', 13, 1, g13['fast_energy_16'], g13['fast_weight_16'], g13middle)
    for sign in (1, -1):
        head = f'mid16_{sign}'
        middle_comparison('group14/'+head, 14, sign, g14[head+'/energies'], g14[head+'/weights'],
                          json.loads(g14[head+'/metadata_json'].tobytes()))

    # All 38 compact signed families have eight physical real subgap rows.
    # Their old integrated packet matrices are deliberately not consulted.
    label = mode['label']; at_one = int(np.flatnonzero(mode['rho_interior'] == 1.)[0])
    for group in range(13, 33):
        for sign in ((1,) if table[group]['angular_eigenvalue'] == 0 else (1, -1)):
            ix = np.flatnonzero((label[:, 0] == group) & (label[:, 1] == sign)
                               & (label[:, 5] > 1) & (label[:, 5] < label[:, 3]))
            if len(ix) != 8: raise ValueError('inherited real subgap inventory changed')
            add(f'subgap8/{group}_{sign}', group, sign, label[ix, 5], label[ix, 6],
                mode['source_covariance'][ix], mode['source_projector'][ix],
                phi=mode['interior'][ix, at_one],
                provenance='eight authenticated physical subgap mode rows; coarse source control only')
    if signature != {p: sha(root, p) for p in (*NUMERICAL_OWNERS, *RECORDS)}:
        raise ValueError('numerical owner/input changed during source preparation')
    meta = {'signature': signature, 'panels': descriptions, 'kernel_order': ['rho', 'p_parallel', 'T01', 'p_perp'],
        'input_payloads': [receipts[p]['payload'] for p in RECORDS[:4]]+[retained_meta['group14_input']],
        'geometry': {'rho': 1., 'radius_squared': 2., 'axial_squared': 3*np.pi/2-4},
        'state_changed': False, 'history_selected': False, 'source_integration_scope': 'finite subtracted modes only',
        'old_generators_rerun': False}
    arrays['metadata_json'] = np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)
    return arrays
