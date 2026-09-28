"""Compose explicit, disjoint incoming source refinements and their scope."""
import numpy as np

from .nsc_incoming_joint_constraints import baseline_matter_source, source_action_gradient
from .nsc_incoming_cauchy_jets import incoming_cauchy_jets


def refined_matter_source(finite, subgap, tail, middle14, middle12, low22):
    """Replace only the named numerical approximations, without double counts."""
    expected=['rho','p_parallel','T01','p_perp']
    for record,group in ((middle14,14),(middle12,12),(low22,22)):
        if record['group']!=group or record['kernel_order']!=expected:
            raise ValueError('explicit matching group and stress order required')
    for record in (middle14,middle12):
        if (record['old_order'],record['new_order'])!=(16,24):
            raise ValueError('explicit order24-minus16 correction required')
    if middle14['interval']['left']!=16. or middle14['interval']['right']!=160.:
        raise ValueError('group14 middle interval changed')
    if middle12['interval']['left']!=40. or middle12['interval']['right']!=320.:
        raise ValueError('group12 middle interval changed')
    if low22['interval']!=[32.,40.] or low22['numerical_Riccati_order']!=24:
        raise ValueError('group22 low interval/source order changed')
    previous=baseline_matter_source(finite,subgap,tail)
    updates={
        'group14_middle_16_to_24':np.asarray(middle14['refined48_node_vacuum_correction']),
        'group12_middle_16_to_24':np.asarray(middle12['vacuum_order24_minus16']),
        'group22_low32_40_to_direct24':np.asarray(low22['new_vacuum_minus_archived_LOW32']),
    }
    checks={
        'group14_explicit_difference': np.max(abs(np.asarray(middle14['additively_corrected_middle_approximation'])-np.asarray(middle14['old16_middle_source'])-updates['group14_middle_16_to_24'])),
        'group12_explicit_difference': np.max(abs(np.asarray(middle12['corrected_middle_approximant'])-np.asarray(middle12['old16_middle_approximant'])-updates['group12_middle_16_to_24'])),
        'group22_explicit_difference': np.max(abs(np.asarray(low22['new_order24_vacuum_source'])-np.asarray(low22['archived_LOW32_source'])-updates['group22_low32_40_to_direct24'])),
    }
    if max(checks.values())>3e-22:
        raise ValueError('source update does not match its explicitly recorded original')
    if any(v.shape!=(4,) or not np.isfinite(v).all() for v in updates.values()):
        raise ValueError('finite four-component source increments required')
    increment=sum(updates.values());updated=np.asarray(previous['stress_approximant'])+increment
    return {**previous,'previous_stress_approximant':np.asarray(previous['stress_approximant']),
            'stress_approximant':updated,'refinement_parts':updates,'total_stress_increment':increment,
            'action_gradient_increment':source_action_gradient(incoming_cauchy_jets(),increment),
            'composition_residuals':checks,
            'unresolved':['remaining low/subgap modal accuracy','rigorous source-quadrature error',
                          'remaining middle-group approximation budget','retained angular/compact scope versus complete physical inventory'],
            'full_source_error_bound':None,
            'scope':{'same_target_C0_source_law':True,'old_receipts_overwritten':False,
                     'new_thermal_term_assigned_zero':False,'group22_thermal_remainder_explicit':True,
                     'local_light_geometry_included':False,'source_accuracy_complete':False}}
