"""Finite incoming source from authenticated spectral panels and allocations.

This reader performs quadrature contraction only. It selects no geometry,
solves no mode equation and does not substitute a finite sum for a completed
renormalized source. In particular, subgap source convergence and the omitted
infinite energy/angular/compact remainder remain separate requirements.
"""
import json
from pathlib import Path

import numpy as np

from .nsc_incoming_spectral_panels import load, sha, record
from .nsc_incoming_state_moments import incoming_group_factor, AXIAL, RADIUS
from .nsc_magnetic_light_reference import magnetic_light_spectrum, homogeneous_magnetic_restoration
from .nsc_horizon_source import conformal_stress


PANEL_DIGEST = 'efbce40b2f959c7d0294689661d4e2e7193ade965ece0918c72c5a072ec8ad73'
PANEL_PATH = 'results/development/artifacts/nsc-incoming-spectral-panels.'+PANEL_DIGEST+'.npz'
ORDER = ('rho', 'p_parallel', 'T01', 'p_perp')


def read_panels(root):
    root = Path(root)
    if sha(root, PANEL_PATH) != PANEL_DIGEST:
        raise ValueError('finite source payload digest changed')
    arrays = load(root/PANEL_PATH)
    meta = json.loads(arrays['metadata_json'].tobytes())
    for path, expected in meta['signature'].items():
        if sha(root, path) != expected:
            raise ValueError('numerical owner or input changed: '+path)
    for item in meta['input_payloads']:
        if sha(root, item['path']) != item['sha256']:
            raise ValueError('upstream payload changed: '+item['path'])
    return arrays, meta


def selected_pieces(arrays, descriptions, group, sign, *, refined=True):
    """Return disjoint quadrature pieces; fine panels replace, never add."""
    pieces = []
    def whole(name):
        if name not in descriptions:
            raise ValueError('missing required source panel '+name)
        pieces.append((name, np.ones(len(arrays[name+'/energies']), dtype=bool)))
    if group == 13:
        whole('group13/low')
        whole('group13/mid' if refined else 'group13/mid16')
    elif group == 14:
        name = f'group14/low16_{sign}'
        E = arrays[name+'/energies']
        if refined:
            mass = descriptions[name]['mass']
            mask = ~(((E > .25) & (E < 1.)) | ((E > mass) & (E < 8.)))
            pieces.append((name, mask))
            whole(f'group14/low32_{sign}')
        else:
            whole(name)
        whole(f'group14/mid{24 if refined else 16}_{sign}')
    else:
        for family in ('low', 'mid'):
            fine = f'{family}_ref/{group}_{sign}'
            whole(fine if refined and fine in descriptions else f'{family}/{group}_{sign}')
    if group >= 13:
        whole(f'subgap8/{group}_{sign}')
    if any(descriptions[n]['group'] != group or descriptions[n]['sign'] != sign for n, _ in pieces):
        raise ValueError('source pieces do not belong to the requested signed channel')
    return pieces


def integrate_pieces(arrays, pieces, quantity='kernels'):
    answer = np.zeros(4)
    for name, mask in pieces:
        answer += np.einsum('n,nv->v', arrays[name+'/weights'][mask], arrays[name+'/'+quantity][mask])
    return answer


def light_allocations(config):
    """Same incoming surface, physical LLL counted once, no compact terms."""
    charge = abs(config['magnetic_flux'])
    if charge != 4:
        raise ValueError('only the locked q=4 physical branch is owned')
    kappa, omega = config['surface_gravity'], config['omega']
    tu = charge*kappa*kappa/(48*np.pi); tv = omega*omega*tu
    factor = 1/(4*np.pi*RADIUS**2*AXIAL**2)
    state = factor*np.array([tu+tv, tu+tv, tv-tu, 0.])
    A, Ap, App = 4-3*np.pi/2, 6-3*np.pi/2, 3-3*np.pi/2
    tensor = conformal_stress(A, Ap, App, 0., 0., central_charge=charge)
    # The inherited stationary radial stress, transformed on this actual
    # rho=1 surface; no r=0 or source-selected-neck tensor is copied here.
    uu, uv, vv = tensor
    geometry = factor*np.array([uu-2*uv+vv, uu+2*uv+vv, vv-uu, 0.])
    local = homogeneous_magnetic_restoration(3*np.pi/4, magnetic_light_spectrum(4))
    nonzero = np.array([local['rho'], local['p_parallel'], local['T01'], local['p_sphere']])
    return {'LLL_state': state.tolist(), 'physical_LLL_geometry': geometry.tolist(),
            'nonzero_magnetic_angular_restoration': nonzero.tolist(),
            'sum': (state+geometry+nonzero).tolist(), 'restoration_details': local,
            'compact_complement_included': False, 'LLL_count': 1}


def assemble(root):
    root = Path(root); arrays, meta = read_panels(root)
    table = record(root, 'results/development/nsc-mode-resolved-cauchy-state.json')['channels']
    config = record(root, 'results/development/nsc-compact-matched-restart.json')['scattering_provenance']['config']
    groups = {}; total = np.zeros(4); base_total = np.zeros(4); error_indicator = np.zeros(4)
    weight_residual = 0.; count = 0; maxima = {}; high_context = {}
    for name, desc in meta['panels'].items():
        diag = json.loads(arrays[name+'/diagnostics_json'].tobytes())
        target = high_context if 'middle' in desc['method'] or 'order16' in desc['method'] else maxima
        for key, value in diag.items():
            target[key] = max(target.get(key, 0.), float(value))
    for group in range(1, 33):
        channel = table[group]; signs = (1,) if channel['angular_eigenvalue'] == 0 else (-1, 1)
        factor = incoming_group_factor(channel, signs)
        fine = np.zeros(4); base = np.zeros(4); rounding = np.zeros(4); subgap = np.zeros(4)
        raw_current_delta = 0.; raw_middle_delta = np.zeros(4); signed = {}
        for sign in signs:
            pieces = selected_pieces(arrays, meta['panels'], group, sign)
            coarse = selected_pieces(arrays, meta['panels'], group, sign, refined=False)
            fine += integrate_pieces(arrays, pieces); base += integrate_pieces(arrays, coarse)
            rounding += integrate_pieces(arrays, pieces, 'arithmetic_indicator')
            nodes = sum(int(mask.sum()) for _, mask in pieces); count += nodes
            width = sum(float(arrays[n+'/weights'][mask].sum()) for n, mask in pieces)
            end = 320. if group in (10, 11, 12, 31, 32) else 160.
            weight_residual = max(weight_residual, abs(width-end))
            signed[str(sign)] = {'pieces': [{'panel': n, 'selected_rows': int(m.sum())} for n, m in pieces],
                                 'rows': nodes, 'integrated_dE_measure': width, 'upper_energy_endpoint': end}
            for n, mask in pieces:
                w = arrays[n+'/weights'][mask]
                if n.startswith('subgap8/'):
                    subgap += np.einsum('n,nv->v', w, arrays[n+'/kernels'][mask])
                if n+'/raw_kernels' in arrays:
                    raw_current_delta += float(w @ (arrays[n+'/raw_kernels'][mask, 2]-arrays[n+'/kernels'][mask, 2]))
                if n+'/direct_kernels' in arrays:
                    raw_middle_delta += np.einsum('n,nv->v', w,
                        arrays[n+'/direct_kernels'][mask]-arrays[n+'/kernels'][mask])
        fine *= factor; base *= factor; rounding *= factor; subgap *= factor
        total += fine; base_total += base; error_indicator += rounding
        groups[str(group)] = {'mass': channel['compact_mass'], 'angular_magnitude': channel['angular_eigenvalue'],
            'factor': factor, 'signed_inventory': signed, 'finite_moments': fine.tolist(),
            'base_control_moments': base.tolist(), 'quadrature_change': (fine-base).tolist(),
            'arithmetic_indicator_not_error_bound': rounding.tolist(), 'subgap8_contribution': subgap.tolist(),
            'raw_low_current_minus_stable_current': factor*raw_current_delta,
            'direct_middle_minus_stable_middle': (factor*raw_middle_delta).tolist(),
            'full_source_converged': False}
    allocation = light_allocations(config)
    return {'groups': groups, 'kernel_order': list(ORDER), 'selected_spectral_rows': count,
            'finite_mode_sum': total.tolist(), 'base_mode_sum': base_total.tolist(),
            'net_quadrature_change': (total-base_total).tolist(),
            'sum_absolute_group_quadrature_changes': np.sum([abs(np.array(g['quadrature_change'])) for g in groups.values()], axis=0).tolist(),
            'arithmetic_indicator_not_error_bound': error_indicator.tolist(),
            'light_allocation': allocation,
            'finite_quantum_source_approximant': (total+np.array(allocation['sum'])).tolist(),
            'residuals': {'spectral_interval_measure': weight_residual},
            'low_modal_diagnostics': maxima, 'middle_inherited_approximation_diagnostics': high_context,
            'payload_metadata': meta}
