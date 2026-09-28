"""Reference-action change at fixed incoming intrinsic metric and C0.

Reuse the order-four spatial projector, with its complete local momentum
line and inherited multiplicities. This is a change of subtraction, not a
new physical state, reference-band endpoint result, or constraint root.
"""
import numpy as np

from .nsc_incoming_cauchy_jets import IncomingCauchyJets, _raw_vertex_jets
from .nsc_spatial_reference_symbol import reference_projector


def signed_inventory(channels):
    """Expand the declared group multiplicity once, including the LLL label."""
    rows = []
    for group, channel in enumerate(channels):
        mass = float(channel['compact_mass'])
        angular = float(channel['angular_eigenvalue'])
        copies = float(channel['copy_count'])
        degeneracy = float(channel['degeneracy'])
        if not np.isfinite([mass, angular, copies, degeneracy]).all() or min(mass, angular) < 0 or min(copies, degeneracy) <= 0:
            raise ValueError('finite inherited nonnegative labels and positive multiplicities required')
        signs = (1,) if angular == 0 else (-1, 1)
        for sign in signs:
            rows.append({'group': group, 'mass': mass, 'angular': sign*angular,
                         'multiplicity': copies*degeneracy/len(signs)})
    return rows


def reference_vertex_change(domain, momenta, masses, angular):
    """Ordinary bulk insertion tr[(P_new-P_old) V_B], all four raw fields.

    Complete unweighted phase-space integration makes the Weyl trace cyclic.
    The audited paired band exchange has no compact-test bulk force. Neither
    fact deletes finite-cutoff or physical endpoint terms. Zeroth order is
    identically unchanged on this surface and is checked before exclusion;
    higher orders are subtracted before summation to retain small UV terms.
    The physical LLL is excluded here: its chartwise reference change is zero,
    while its geometric action stays with the existing local owner.
    """
    if not isinstance(domain, IncomingCauchyJets):
        raise ValueError('same-surface incoming Cauchy jets required')
    domain.validate()
    raw = tuple(np.asarray(x) for x in (momenta, masses, angular))
    if any(np.iscomplexobj(x) for x in raw):
        raise ValueError('real momentum and inherited channel labels required')
    k, m, ell = np.broadcast_arrays(*[x.astype(float) for x in raw])
    if k.ndim != 1 or not k.size or not np.isfinite(k+m+ell).all() or np.any(m < 0) or np.any((m == 0) & (ell == 0)):
        raise ValueError('finite non-LLL vector of gapped retained labels required')
    old = reference_projector(*domain.baseline, k, m, ell)
    new = reference_projector(*domain.fields, k, m, ell)
    leading = float(np.max(abs(new['orders'][0]-old['orders'][0])))
    if leading > 3e-12:
        raise ValueError('incoming zeroth reference must be unchanged')
    difference = new['orders'][1:]-old['orders'][1:]
    vertices = np.stack([v.value for v in _raw_vertex_jets(domain.fields, k, m, ell)], axis=1)
    insertion = np.einsum('onij,nbji->onb', difference, vertices)
    formal = max(float(np.max(v['scaled'])) for p in (old, new) for order in p['residuals'] for v in order.values())
    return {'orders': insertion.real, 'leading_projector_change': leading,
            'imaginary_residual': float(np.max(abs(insertion.imag))),
            'formal_projector_residual': formal}


def incoming_reference_response(domain, channels, *, nodes=32, scale_factor=1., chunk_size=16):
    """Integrate the positive reference-action change per dT dz.

    k = a sqrt(m^2+(ell/r)^2) scale_factor tan(theta), 0<theta<pi/2,
    with both k signs. The scale is numerical and covers the same full line
    for every positive value. Compare node counts AND map scales externally;
    agreement is a numerical indicator, not an exact-mode/Hadamard bound.
    """
    if not isinstance(domain, IncomingCauchyJets):
        raise ValueError('same-surface incoming Cauchy jets required')
    domain.validate()
    if isinstance(nodes, bool) or not isinstance(nodes, int) or nodes < 8:
        raise ValueError('at least eight Gauss nodes required')
    if isinstance(chunk_size, bool) or not isinstance(chunk_size, int) or chunk_size < 1:
        raise ValueError('positive integer chunk size required')
    if not np.isfinite(scale_factor) or scale_factor <= 0:
        raise ValueError('positive finite numerical map scale required')
    a = domain.normal_geometry()['intrinsic_a']
    r = domain.normal_geometry()['intrinsic_r']
    x, w = np.polynomial.legendre.leggauss(nodes)
    theta = (x+1)*np.pi/4
    rows = signed_inventory(channels)
    if not rows:
        raise ValueError('declared channel inventory required')
    groups = {}
    samples = []
    maxima = dict(leading_projector_change=0., imaginary_residual=0., formal_projector_residual=0.)
    for row in rows:
        mass, angular = row['mass'], row['angular']
        integral = np.zeros((4, 4))
        if mass or angular:
            scale = scale_factor*a*np.hypot(mass, angular/r)
            positive = scale*np.tan(theta)
            weights = scale*np.pi/4*w/np.cos(theta)**2
            k = np.r_[-positive, positive]
            weights = np.tile(weights, 2)
            insertions = []
            for start in range(0, len(k), chunk_size):
                end = start+chunk_size
                result = reference_vertex_change(domain, k[start:end], mass, angular)
                integral += np.einsum('n,onb->ob', weights[start:end], result['orders'])
                insertions.append(result['orders'])
                for key in maxima:
                    maxima[key] = max(maxima[key], result[key])
            samples.append({**row, 'momenta': k, 'weights': weights,
                            'insertions': np.concatenate(insertions, axis=1)})
        groups.setdefault(str(row['group']), np.zeros((4, 4)))
        groups[str(row['group'])] += row['multiplicity']/(2*np.pi)*integral
    total = sum(groups.values())
    return {'reference_action_gradient_change': total.sum(axis=0),
            'formal_order_gradient_changes': total,
            'group_formal_order_gradient_changes': groups,
            'raw_metric_order': ['N', 'beta', 'a', 'r'], 'formal_orders': [1, 2, 3, 4],
            'maxima': maxima, 'quadrature': {'nodes_per_half_line': nodes, 'scale_factor': scale_factor},
            'changed_normal_entries': domain.changed_normal_entries(),
            'quadrature_samples': samples,
            'scope': {'fixed_intrinsic_surface_and_C0': True,
                      'positive_reference_action_change_only': True,
                      'LLL_chartwise_reference_change': 0,
                      'LLL_geometric_allocation_included': False,
                      'band_or_physical_endpoints_included': False,
                      'physical_mode_error_bounded': False,
                      'constraint_root_or_metric_evolution': False}}
