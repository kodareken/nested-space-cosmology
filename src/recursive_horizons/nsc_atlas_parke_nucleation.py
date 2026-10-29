"""Mass-free Parke bubble and its admission into the fixed volume source.

The radius is determined in the stipulated constant-vacuum thin-wall
model. That result is not silently promoted to an instanton of the
nonstationary scalar/volume preparation or to a reheated cosmology.
"""
from fractions import Fraction as F
import math
from scipy.optimize import brentq
from .nsc_atlas_flux_conversion import isolated_bounce
from .nsc_atlas_flux_transition import FluxLaw
from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_inherited_matter import exact_initial_geometry


def thin_wall(Ui=7/32,Uo=21/32,sigma=1/5,kappa=1.):
    if not (Uo>Ui>=0 and sigma>0 and kappa>=0):raise ValueError('downward vacuum step and positive tension required')
    delta=Uo-Ui;k=kappa*sigma/2
    bi=delta/(3*sigma)+kappa*sigma/4
    bo=delta/(3*sigma)-kappa*sigma/4
    if bo<=0:raise ValueError('outside the ordinary positive-cap branch')
    R=1/math.sqrt(kappa*Ui/3+bi*bi)
    betai,betao=bi*R,bo*R
    return dict(radius=R,beta_inside=betai,beta_outside=betao,
        inverse_radius_squared=delta**2/(9*sigma*sigma)+kappa*(Ui+Uo)/6+(kappa*sigma)**2/16,
        angular_residual=betai-betao-k*R,
        inside_constraint=betai**2+kappa*Ui*R*R/3-1,
        outside_constraint=betao**2+kappa*Uo*R*R/3-1,
        inside_Schwarzschild_offset=0.,outside_Schwarzschild_offset=0.,
        vacuum_deficit_term=-kappa*delta*R**3/6,
        redshifted_shell_term=k*R*R*betai,binding_term=-k*k*R**3/2,
        interior_Misner_Sharp_mass=kappa*Ui*R**3/6,
        R_over_false_vacuum_radius=R*math.sqrt(kappa*Uo/3),
        physical_scope='Constant-vacuum, constant-tension, zero-seed-mass thin-wall branch; not a derived (12,2) field profile.')


def geometry_record():
    r=thin_wall();R=r['radius'];k=.1
    # Frozen scalar snapshot represented as constants only for this control.
    existing=isolated_bounce(FluxLaw(bare_lambda=3/128),9/8)
    Ui,Uo=7/32,21/32
    independent_action=2*math.pi**2*.2*R**3+12*math.pi**2*((r['beta_inside']**3-1)/Ui-(r['beta_outside']**3-1)/Uo)
    historical_ruler=1/math.sqrt(float(exact_initial_geometry()['curvature']))
    return dict(schema='nsc-atlas-parke-bubble-v1',inputs=dict(Ui='7/32',Uo='21/32',sigma='1/5',kappa=1),
        bubble=r,R_squared_exact='57600/39169',R_expression='240/sqrt(39169)',
        inverse_curvature_formula='R_c^-2=(DeltaU/(3sigma))^2+kappa*(Ui+Uo)/6+(kappa*sigma)^2/16',
        legacy_comparison=dict(original_R_in_model_units=1.9,R_ratio=R/1.9,
            old_Lref_over_Lunit=historical_ruler,new_R_over_old_Lref=R/historical_ruler,
            old_R_over_old_Lref=1.9/historical_ruler,
            scope='Historical conversion only: the old curvature ruler itself depended on the removed R,M. Neither 1.9 nor R_c was originally measured in L_ref units.'),
        energy_balance_residual=sum(r[key] for key in ('vacuum_deficit_term','redshifted_shell_term','binding_term')),
        action_control=dict(existing_coefficient=existing['action'],independent_Parke_coefficient=independent_action,
                            interpretation='Auxiliary vacuum-bubble action, not a (12,2) nucleation rate.'),
        wall_trajectory='R(tau_wall)=R_c*cosh(tau_wall/R_c), with Rdot(0)=0',
        open_interior='a(t)=sinh(H_i*t)/H_i, H_i^2=Ui/3; wall proper time and interior cosmic time differ',
        open_birth_Omega_k_limit=1.,open_vacuum_Omega_k='sech(H_i*t)^2',
        fixed_original_R_or_M_used_to_determine_radius=False,
        actual_volume_instanton_derived=False)


def source_record():
    model=VolumeParent();top,minimum=model.vacuum_roots()
    Vstar,Vtop=minimum['vacuum'],top['vacuum'];target=3/8
    necessary_exit=Vstar+target
    necessary_s=brentq(lambda s:model.potential(s)-necessary_exit,-.3,minimum['s'],xtol=1e-14)
    m_s=math.sqrt(minimum['canonical_mass_squared']);m_chi=1e9*math.exp(-6*minimum['s'])
    return dict(schema='nsc-atlas-parke-reheating-admission-v1',
        vacuum_admission=dict(prepared_middle_density=7/32,prepared_outer_density=21/32,
            relaxed_middle_density=1/8,relaxed_outer_density=9/32,
            middle_force_dU_d_omega_phi=math.sqrt(3)/4,outer_force_dU_d_omega_phi=math.sqrt(3)/2,
            fixed_volume_minimum=Vstar,fixed_volume_barrier=Vtop,interior_density_minus_minimum=7/32-Vstar,
            endpoints_are_vacua_of_fixed_source=False,
            wall_condition='sigma(s_wall)=.2*exp(-6s_wall); K[normal derivative s]=-6sigma. The constant reference tension alone does not solve this profile.'),
        empty_start=dict(Schwarzschild_offset=0,initial_chi=0,initial_chi_momentum=0,
            zero_matter_is_an_explicit_branch_assumption_not_implied_by_zero_mass_offset=True,
            classical_matter_energy=0,all_time_matter_energy=0,
            identity='E=E_i*exp[-3(N-N_i)-6(s-s_i)]; E_i=0 implies E=0 for every prescribed s(t).',
            reference_energy=target,reference_event_reached=False,
            scope='Homogeneous post-wall fixed-flux classical field/adiabatic matter system without incoming sources; not a claim that all quantum reheating is impossible.'),
        energy_budget=dict(required_total_density_at_reference=necessary_exit,
            literal_child_total_density=7/32,matter_target_over_literal_child=12/7,
            minimum_required_rest_exit_density=necessary_exit,
            necessary_left_branch_s=necessary_s,necessary_s_is_a_CDL_turning_point=False,
            available_excess_from_computed_well_barrier=Vtop-Vstar,
            shortfall_even_at_full_conversion=target-(Vtop-Vstar),
            literal_Ui_excess=7/32-Vstar,parent_Uo_excess=21/32-Vstar,
            domain='The barrier bound applies only to a zero-velocity exit between the computed minimum and barrier with no incoming energy. No unknown CDL exit is replaced by the HM top.',
            parent_energy_caveat='Parent energy could in principle cover the budget, but its transfer into the homogeneous child has not been derived.'),
        perturbative_reheating=dict(volume_small_oscillation_mass=m_s,retained_chi_mass=m_chi,
            pair_threshold=2*m_chi,pair_decay_open=m_s>=2*m_chi,
            scope='Ordinary on-shell one-volume-quantum to two-chi decay near the minimum; nonperturbative and curved-background production are not computed.'),
        downstream=dict(C_over_E0_at_reference=None,Omega_k_observer=None,w_s_redshift=None,
            SN_DESI_refit=None,reason='No compatible (12,2) CDL turning state or E=3/8 event was obtained. A data fit from a reset state would not answer the requested route.'),
        route_dispositions=[
            'Parke radius: fixed for the stipulated vacuum-wall control, but its endpoints are not stationary vacua of the unchanged scalar/volume action.',
            'Regular mass-free interior: M_offset=0 is compatible with wall energy balance; it does not mean zero vacuum or Misner-Sharp energy.',
            'Classical reheating of the retained empty matter field: E=0 is invariant, so E=3/8 is never reached.',
            'Volume energy counted as matter: rolling within the computed well supplies at most .1335536533 above Vstar before expansion losses, below 3/8.',
            'Perturbative decay into the retained massive scalar: the one-to-two channel is below threshold; no alternative reheating law is supplied.'],
        physical_goal_achieved=False,new_observation_fits=0)
