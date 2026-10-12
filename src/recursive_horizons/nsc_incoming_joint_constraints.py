"""Join owned incoming source and action contributions without declaring a root."""
import numpy as np

from .nsc_incoming_cauchy_jets import IncomingCauchyJets


def source_action_gradient(domain, stress):
    """Convert (rho,p_parallel,T01,p_perp) into raw KS action coefficients."""
    if not isinstance(domain, IncomingCauchyJets):
        raise ValueError('same intrinsic incoming surface required')
    domain.validate()
    tensor = np.asarray(stress)
    if tensor.shape != (4,) or np.iscomplexobj(tensor) or not np.isfinite(tensor).all():
        raise ValueError('four real finite stress moments required')
    N, _, a, r = [float(f.value[0, 0, 0].real) for f in domain.fields]
    rho, parallel, current, sphere = tensor
    return np.array([-4*np.pi*a*r*r*rho, -4*np.pi*a*a*r*r*current,
                     4*np.pi*N*r*r*parallel, 8*np.pi*N*a*r*sphere])


def baseline_matter_source(finite, subgap, tail):
    """Use stored finite moments with existing replacements, excluding geometry.

    Group13's exact closed-channel current is used explicitly; the record's
    raw numerical zero residual is retained separately. Physical LLL state
    is added once. All light geometry is supplied by the local action owner.
    """
    order = ['rho', 'p_parallel', 'T01', 'p_perp']
    if finite['kernel_order'] != order or tail['kernel_order'] != order or subgap['source']['kernel_order'] != order:
        raise ValueError('common stress component order required')
    if subgap['domain']['group'] != 13 or finite['light_allocation']['LLL_count'] != 1:
        raise ValueError('owned group13 replacement and one LLL state required')
    if subgap['source']['exact_physical_subgap_T01'] != 0.:
        raise ValueError('closed subgap current identity required')
    original = np.asarray(finite['groups']['13']['subgap8_contribution'])
    coarse = np.asarray(subgap['source']['coarse8_physical_contribution'])
    if np.max(abs(original-coarse)) > 3e-13:
        raise ValueError('group13 replacement does not match the original panel')
    correction = np.array(subgap['source']['refined_minus_coarse8'], float)
    raw_current_residual = float(correction[2])
    correction[2] = 0.  # Explicit exact closed-channel identity, not thresholding.
    parts = {'finite_nonLLL': np.asarray(finite['finite_mode_sum']),
             'group13_replacement': correction,
             'paired_approximate_tail': np.asarray(tail['total_tail_by_order']['13']),
             'LLL_state': np.asarray(finite['light_allocation']['LLL_state'])}
    if any(v.shape != (4,) or not np.isfinite(v).all() for v in parts.values()):
        raise ValueError('finite four-component owned source contributions required')
    return {'parts': parts, 'stress_approximant': sum(parts.values()),
            'group13_raw_current_zero_residual': raw_current_residual,
            'full_source_error_bound': None,
            'unresolved': ['finite-band mode approximation', 'other compact subgap source refinements',
                           'exact vacuum tail versus archived mode approximation',
                           'retained angular/compact scope versus complete physical inventory']}


def joint_constraint_approximant(domain, matter, local, reference_change):
    """Compose lapse/shift bulk coefficients, preserving the source OPEN gate.

    Fixed C0 means the physical source insertion is unchanged at this surface.
    The subtraction-action change is positive. The local action already owns
    light geometry and compact coefficients. No inferred or fitted force enters.
    """
    if local['constraint_order'] not in (('N', 'beta'), ['N', 'beta']):
        raise ValueError('local lapse/shift order required')
    if reference_change['raw_metric_order'] != ['N', 'beta', 'a', 'r']:
        raise ValueError('all four reference metric directions required')
    if reference_change['changed_normal_entries'] != domain.changed_normal_entries():
        raise ValueError('reference and local domain must share the supplied normal data')
    matter_gradient = source_action_gradient(domain, matter['stress_approximant'])[:2]
    local_gradient = np.asarray(local['local_action_gradient'])
    delta_reference = np.asarray(reference_change['reference_action_gradient_change'])[:2]
    gradient = matter_gradient+local_gradient+delta_reference
    return {'constraint_order': ['N', 'beta'],
            'matter_baseline_action_gradient': matter_gradient,
            'local_action_gradient': local_gradient,
            'reference_action_gradient_change': delta_reference,
            'action_gradient_approximant': gradient, 'force_approximant': -gradient,
            'maximum_absolute_approximant': float(np.max(abs(gradient))),
            'stationarity_tolerance': 3e-11,
            'physical_constraint_status': 'OPEN', 'certified_constraint_residual': None,
            'full_source_error_bound': matter['full_source_error_bound'],
            'normal_data_selected_as_physical': False,
            'extended_stationarity_status': 'OPEN', 'metric_timestep': False}
