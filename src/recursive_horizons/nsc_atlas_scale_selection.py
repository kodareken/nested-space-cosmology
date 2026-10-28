"""Computed scale-selection ladder for the unchanged inherited parent.

No target value sets a predicted length. Inverse and supplemental quantum
examples carry explicit labels and cannot pass as an inheritance closure.
"""
from fractions import Fraction as F
import math

from .nsc_atlas_inherited_matter import exact_initial_geometry
from .nsc_atlas_flux_conversion import ScalarFlux,isolated_bounce,sector_minimum
from .nsc_atlas_flux_transition import FluxLaw
from .nsc_atlas_shared_scalar_slice import SharedScalarSlice

LIGHT_SPEED=299792458.
PLANCK_H=6.62607015e-34
NEWTON_G=6.67430e-11
NEWTON_G_STANDARD_UNCERTAINTY=1.5e-15
MPC_METRES=1e6*648000/math.pi*149597870700.
TARGET_LAMBDA=1.0904e-52


def calculation(*,target=TARGET_LAMBDA,observer_scale=1.4656003507830409):
    d=exact_initial_geometry();lam=float(d['vacuum']);C=float(d['curvature'])
    a0=1/math.sqrt(C);coefficient=lam/C
    hbar=PLANCK_H/(2*math.pi)
    planck2=hbar*NEWTON_G/LIGHT_SPEED**3
    kappa_natural=8*math.pi*planck2
    kappa_SI=8*math.pi*NEWTON_G/LIGHT_SPEED**4
    assert abs(kappa_natural/(kappa_SI*hbar*LIGHT_SPEED)-1)<3e-16
    eps,q,e,omega=.25,1.5,.5,1e9
    qhat=[F(3,2)-j*F(1,2) for j in range(4)]
    vacua=[F(1,4)*x*x/2 for x in qhat]
    # Formal added S4-period condition: e_phys integral(F4)=2*pi*n.
    # F_hat at the scalar minimum is epsilon*q_hat, NOT q_hat.
    ellhat=math.sqrt(3/lam)
    four_volume=8*math.pi**2*ellhat**4/3
    action_ratio_for_n1=2*math.pi/(e*(eps*q)*four_volume)
    Lunit_n1=math.sqrt(kappa_natural*action_ratio_for_n1)
    Lambda_n1=lam/Lunit_n1**2
    invariant=3*e/q
    assert abs(Lambda_n1*planck2-invariant)<3e-16
    inverse_Lunit=math.sqrt(lam/target)
    inverse_Lref=a0*inverse_Lunit
    # Classical-normalised instanton control, not the interacting scalar-parent rate.
    law=FluxLaw(kappa=1.,bare_lambda=0.,charge=e,tension=.2)
    bounces=[isolated_bounce(law,f) for f in (1.5,1.,.5)]
    atlas=SharedScalarSlice(scalar=ScalarFlux(omega=omega))
    domain=atlas.parent_causal_domain()
    B=domain['initial_parent_angle']
    volume=2*math.pi*a0**3*(B-math.sin(2*B)/2)
    # Conditional number-eigenstate illustration; the prepared classical real
    # scalar is not such a state and supplies no exact conserved particle count.
    Lunit_occupation1=math.sqrt(kappa_natural*omega/(float(d['excitation'])*volume))
    witnesses=[]
    for Lref in (1.,2.):
        Lunit=Lref/a0
        qphys=q/(math.sqrt(kappa_natural)*Lunit)
        ephys=e/(math.sqrt(kappa_natural)*Lunit)
        witnesses.append(dict(arbitrary_reference_metres=Lref,Lambda_m_inverse2=coefficient/Lref**2,
            recovered_dimensionless_charge=math.sqrt(kappa_natural)*Lunit*ephys,
            recovered_dimensionless_flux=math.sqrt(kappa_natural)*Lunit*qphys,
            electric_flux_integer=qphys/ephys,
            interpretation='Counterexample to unique SI identification, not an astronomical prediction.'))
    equality_A=(4/3)**(1/3)
    required_frequency=omega*LIGHT_SPEED/inverse_Lunit
    rows=[
        dict(route='a_flux_charge_quantization',rank=1,
             computed=dict(electric_integer=q/e,scalar_vacuum_ladder=[str(v) for v in vacua],
                           scalar_lattice_minimum=sector_minimum(q,e,kappa=eps)['minimum_vacuum'],
                           scalar_second_difference='1/16',Lambda_Lref_squared=coefficient),
             failure='q/e=3 and the charge steps leave the fundamental SI charge and reference length free.'),
        dict(route='b_inherited_wall_EOS_radius',rank=2,
             computed=dict(radii_in_model_units=[1.,1.9],areal_speeds=[atlas.wall(i)['velocity'] for i in (0,1)]),
             failure='The fixed membranes are evolving initial data, not an inherited stable equilibrium selecting metres; the stable surface-scalar branch has a different action/state.'),
        dict(route='c_recursive_fixed_point_and_level',rank=3,
             computed=dict(vacuum_length_ratios=[1.5,2.],last_vacuum='0 at q=0',
                           flux_fixed_point_residual='-e=-1/2'),
             failure='The supplied map fixes no root length, level or full M/sigma/ell fixed point; its nonzero charge shift has no flux fixed point.'),
        dict(route='d_scalar_clock_epoch',rank=4,
             computed=dict(Lambda_c_squared_over_omega_squared=lam/omega**2,
                           density_ratio_initial=4/3,equality_scale=equality_A,
                           equality_redshift_at_calibrated_observer=observer_scale/equality_A-1),
             failure='The epoch and omega/H ratios are dimensionless; no measured SI scalar clock is inherited.'),
        dict(route='e_compact_four_cycle_quantum_period',rank=5,
             computed=dict(extra_premise='A compact Euclidean S4 period with a real contour, admissible saddle and inherited winding; none is established for the excited parent.',
                           formal_Lambda_lPlanck_squared_times_winding=invariant,winding_1_Lref_metres=a0*Lunit_n1,
                           winding_1_Lambda_m_inverse2=Lambda_n1,winding_1_ratio_to_target=Lambda_n1/target,
                           required_winding_if_target_inserted=Lambda_n1/target),
             failure='The formal period gives Lambda*lP^2=1/n, but neither an admissible parent cycle nor its winding is selected; n=1 gives a Planck-scale value.'),
        dict(route='f_semiclassical_transition_probability',rank=6,
             computed=dict(pure_flux_control_actions=[b['action'] for b in bounces],
                           physical_exponent='S_hat*L_unit^2/kappa_natural',
                           first_exponent_if_target_inserted=bounces[0]['action']*inverse_Lunit**2/kappa_natural),
             failure='The action exponent, determinant/state and parent preparation are not a scale-selection equation; isolated vacuum bounces are not this interacting parent.'),
        dict(route='g_quantum_occupation_or_entropy',rank=7,
             computed=dict(parent_initial_volume_model_units=volume,
                           hypothetical_one_quantum_Lref_metres=a0*Lunit_occupation1,
                           hypothetical_one_quantum_Lambda_m_inverse2=lam/Lunit_occupation1**2),
             failure='An inherited occupation/entropy count and its quantum state are absent; setting one quantum changes the prepared classical source.'),
        dict(route='h_anomaly_running_and_dimensional_transmutation',rank=8,
             computed=dict(existing_classical_added_vacuum_trace=0,
                           generic_quantum_stationary_scale='L^2=-b*kappa_natural/(2*S_hat), only if b and the renormalized action/state are fixed'),
             failure='The compatible renormalized quantum action, state and subtraction conditions do not fix b or an absolute vacuum; old spectral closures are different actions.'),
        dict(route='i_induced_gravity_Casimir_compactification',rank=9,
             computed=dict(stable_gravitational_branch_Lambda_M_squared=.114490614635103,
                           fixed_parent_Lambda_Lref_squared=coefficient),
             failure='Boundary state, compactification modulus or source mass in SI remains supplied; the prior branches do not select the fixed parent ruler.'),
        dict(route='j_infinite_hierarchy_global_self_similarity',rank=10,
             computed=dict(reference_length_relation='ell_n=r^n*ell_0'),
             failure='An inherited ratio or infinite operator composition leaves ell_0 and observer level unselected.'),
        dict(route='k_s_one_spectral_fixed_point',rank=11,
             computed=dict(legacy_gap_threshold=3*math.pi/2,legacy_alternative_vacuum_coefficient=18*math.pi),
             failure='The legacy coefficient dictionary retains a length_anchor and has no solved physical fixed point; its action is not the fixed scalar/four-form parent.')]
    for row in rows:row.update(inherited_SI_scale_selected=False,survives=False)
    return dict(schema='nsc-atlas-SI-scale-ladder-v1',target_Lambda_m_inverse2=target,
        fixed_parent=dict(Lambda_model_unit_squared=lam,Lambda_Lref_squared=coefficient,Lref_over_model_unit=a0),
        SI_conventions=dict(h=PLANCK_H,c=LIGHT_SPEED,G=NEWTON_G,G_standard_uncertainty=NEWTON_G_STANDARD_UNCERTAINTY,
                            lPlanck_metres=math.sqrt(planck2),kappa_natural_metres_squared=kappa_natural,
                            conversion='kappa_natural=(8*pi*G/c^4)*(hbar*c)=8*pi*lPlanck^2'),
        route_ladder=rows,
        scale_nonidentifiability_witnesses=witnesses,
        inverse_target_diagnostic=dict(not_a_prediction=True,model_unit_metres=inverse_Lunit,Lref_metres=inverse_Lref,
            Lref_over_lPlanck=inverse_Lref/math.sqrt(planck2),
            canonical_charge_inverse_metres_squared=e/(math.sqrt(kappa_natural)*inverse_Lunit),
            required_scalar_angular_frequency_s_inverse=required_frequency,
            required_scalar_period_years=2*math.pi/required_frequency/(365.25*86400)),
        conclusion='No audited route selects an inherited SI length from the fixed inputs. This does not reject nesting or exclude a future state/quantum/global selection mechanism.',
        scientific_goal_achieved=False)
