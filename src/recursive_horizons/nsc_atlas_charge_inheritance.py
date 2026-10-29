"""Charge transport, distinct from flux discharge and unit normalization.

Only the existing action, its two algebraic vacuum roots, and the stipulated
Parke wall are evaluated. No bounce, field history or scale selector is added.
"""
import math
from fractions import Fraction as F

from .nsc_atlas_flux_transition import FluxLaw
from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_parke_nucleation import thin_wall
from .nsc_atlas_si_hierarchy import A0, units


def volume_charge_ratio(parent_s, child_s, d=12):
    """Regional endpoint comparison in one Einstein-frame reduced action."""
    return math.exp(-d*(child_s-parent_s)/2)


def parke_ruler():
    """Define L_P=R_c without calling the historical curvature preparation.

    The volume well's Lambda is a conversion comparison, not the actual
    vacuum of the auxiliary constant-density Parke interior.
    """
    bubble=thin_wall(Ui=7/32,Uo=21/32,sigma=1/5,kappa=1)
    minimum=next(row for row in VolumeParent().vacuum_roots() if row['kind']=='minimum')
    R=bubble['radius'];square=F(57600,39169)
    return dict(definition='L_P = R_c; a wall areal-radius ruler, not an FRW reference slice',
        inputs=dict(Ui='7/32',Uo='21/32',sigma='1/5',kappa=1),
        L_P_over_L_unit=R,L_P_squared_over_L_unit_squared=str(square),
        late_volume_vacuum=minimum['vacuum'],
        late_volume_Lambda_L_P_squared=minimum['vacuum']*float(square),
        Parke_interior_Lambda_L_P_squared_exact=str(F(7,32)*square),
        Parke_exterior_Lambda_L_P_squared_exact=str(F(21,32)*square),
        log_hierarchy_intercept_eta=math.log(R/5),
        log_hierarchy_intercept_g=math.log(R/2.5),
        old_R_or_M_used=False,FRW_reference_event_selected=False,
        shared_volume_Parke_solution_derived=False,
        scope='The late-volume coefficient is expressed on the Parke ruler only; the source-admission obstruction remains parked.')


def calculation(target_old_log_hierarchy=138.542, target_lambda=1.0904e-52):
    law=FluxLaw();model=VolumeParent();coeff=model.coefficients
    minimum=next(row for row in model.vacuum_roots() if row['kind']=='minimum')
    chain=[{key:value for key,value in law.inherit(1.5,n).items() if key!='Lambda'} for n in range(4)]
    for row in chain:
        row.update(kappa_ratio=1,elementary_charge_ratio=1,eta_ratio_to_parent=1,
                   eta_ratio_to_chain_root=1,screened_scalar_vacuum=row['flux']**2/8)
    p=parke_ruler();u=units()
    # This is explicitly an observation-dependent inverse, never a predictor.
    old_intercept=math.log(A0/5)
    required_log_eta=old_intercept-target_old_log_hierarchy
    required_log_g=2*(math.log(A0/2.5)-target_old_log_hierarchy)
    product=target_lambda*u['ell_e_sigma_m']**2
    parke_target=.5*math.log(p['late_volume_Lambda_L_P_squared']/product)
    return dict(schema='nsc-atlas-charge-inheritance-v1',
        fixed_action=dict(d=12,j=2,Einstein_kappa='kappa_D/V0',
            topological_charge='e_D',Z3='exp(12s)/V0',
            canonical_charge='e_D*sqrt(V0)*exp(-6s)'),
        flux_step=dict(rule='q_child=q_parent-e; e_child=e_parent; kappa_child=kappa_parent',
            r_exact='1',one_line_proof='eta_child/eta_parent=(kappa_child/kappa_parent)*(e_child/e_parent)=1*1=1',
            scope='The implemented constant-coupling flux map; also the two one-sided traces at a finite-energy thin wall of the volume extension.',
            chain=chain),
        volume_extension=dict(regional_ratio='exp[-6*(s_child-s_parent)]',
            wall_continuity='[s]=0: the collar kinetic energy is at least K*(Delta s)^2/(2*epsilon)',
            wall_normal_derivative_jump='K*[normal derivative s]=-6*sigma; derivative, not field, jumps',
            product_rule='product r_k=exp[-6*(s_n-s_root)]',
            repeated_per_level_volume_increment_selected=False,
            regional_r_selected=None,matched_start_to_minimum_ratio=minimum['charge_ratio'],
            matched_start_to_minimum_log_ratio=-6*minimum['s'],minimum_s=minimum['s'],
            matched_start_to_minimum_is_a_nesting_step=False),
        dimensional_reduction=dict(D=16,eta_D='e_D*kappa_D^(4/7)',
            internal_volume='v=V0*exp(12s)/kappa_D^(6/7)',
            exact_ratio='eta_4/eta_D=1/sqrt(v)',
            v_selected=False,repeated_D_to_4_nesting_rule_constructed=False),
        depth=dict(selected_n=None,ancestor_sheet_offset=2,
            ancestor_flux=str(coeff['ancestor_flux']),
            parent_to_zero_flux_steps=3,ancestor_to_zero_flux_steps=5,
            distinction='j is an ancestor flux offset; these counts enumerate flux discharges, not charge-suppression iterations or an observer location.',
            fixed_map_result_for_all_n='ln(eta_n)=ln(eta_root); n*ln(r)=0'),
        root_normalization=dict(selected=False,
            remaining_source='eta_root, equivalently eta_D/sqrt(v_root) in the reduction; not selected by fixed-charge nesting'),
        parke_ruler=p,
        target_diagnostics=dict(old_ruler_target_log=target_old_log_hierarchy,
            required_log_g=required_log_g,required_log_eta=required_log_eta,
            required_eta=math.exp(required_log_eta),
            omitted_factor_correction='ln(eta)=0.5*ln(g)-ln(2), not 0.5*ln(g)',
            observed_Lambda_ell_squared=product,
            Parke_target_log_from_Lambda_obs=parke_target,
            Parke_target_hierarchy_from_Lambda_obs=math.exp(parke_target),
            Parke_target_length_m_from_Lambda_obs=math.sqrt(p['late_volume_Lambda_L_P_squared']/target_lambda),
            ruler_log_shift=math.log(p['L_P_over_L_unit']/A0),
            Parke_target_log_transported_from_rounded_old=target_old_log_hierarchy+math.log(p['L_P_over_L_unit']/A0),
            required_log_eta_on_Parke_ruler=p['log_hierarchy_intercept_eta']-parke_target,
            interpretation='Targets/inversions only. Changing the ruler shifts the target and forward intercept together, without selecting eta.'),
        prediction=dict(ln_eta=None,ln_hierarchy=None,percent_miss=None,
            reason='The implemented flux step has r=1; root eta is absent. The volume extension supplies no fixed regional r or selected depth.'),
        new_evolutions=0,new_fits=0)
