"""Zero/interface tests: late decay, conformal crossover, winding and charge flow.

Signs are Douglas's region/orientation convention, not signs of energy or
automatically particle charges. Model outputs and comparison-only constants
are separate. No conformal source conversion, baryogenesis, or electron
identification is supplied by a notation change.
"""
import math

import numpy as np
import sympy as s

from .nsc_atlas_horizon_response import solve_source
from .nsc_atlas_parent_circulation import VOLUME,energy,proper_rhs
from .nsc_atlas_inherited_charge import harmonic_sector,pair_transport
from .nsc_atlas_s3_polynomials import X,mean


def tail_state(source,u):
    p=source.parent;y=source.modes.sol(u);nu2=source.nu_squared;nu=math.sqrt(nu2)
    b=np.array([3.75,8.75]);q=y[[0,2]];v=y[[1,3]]
    Omega=1/math.cosh(u)
    E=(v*v+(nu2+b*Omega*Omega)*q*q)/2
    amplitude=q*math.cos(nu*u)-v/nu*math.sin(nu*u)
    quadrature=q*math.sin(nu*u)+v/nu*math.cos(nu*u)
    # These are asymptotic phase-quadrature bounds from integrable sech^2 forcing.
    q0=p.initial()[[2,4]];E0=(nu2+b)*q0*q0/2
    phase_tail=b/nu*np.sqrt(2*E0/nu2)*(1-math.tanh(u))
    E_tail=b*E0/nu2*Omega*Omega
    rho_averaged=p.mass_squared*float(sum(E))/(nu2*2*math.pi**2)*Omega**3
    rescaled=math.sqrt(Omega)*q
    rescaled_Omega_derivative=-(v-.5*math.tanh(u)*q)/(math.sqrt(Omega)*math.tanh(u))
    return dict(offset_proper_readout_Ht=u,Omega=Omega,
        amplitude_quadrature=amplitude.tolist(),phase_quadrature=quadrature.tolist(),
        each_quadrature_future_error_upper_bound=phase_tail.tolist(),
        rescaled_oscillator_energy=E.tolist(),energy_future_change_upper_bound=E_tail.tolist(),
        late_phase_average_density_estimate_over_epsilon_squared=rho_averaged,
        late_phase_average_density_estimate_over_epsilon_squared_vacuum=rho_averaged/p.vacuum,
        density_estimate_is_baryon_photon_eta=False,
        H_times_spatial_wavelengths=[2*math.pi/(Omega*math.sqrt(l)) for l in (3,8)],
        conformal_weight_one_mode_values=rescaled.tolist(),
        mode_derivative_with_respect_to_Omega=rescaled_Omega_derivative.tolist())


def reciprocal_fluid(w):
    n=3*(1+w)
    return dict(old_w=w,old_density_exponent=n,new_density_exponent=4-n,new_w=-2/3-w)


def thermodynamic_record(source):
    p=source.parent;rows=[]
    for u in (1.,4.,8.):
        z=source.state(u);a,v,q1,p1,q2,p2=z[:6];h=v/a
        kinetic=p1*p1+p2*p2;mass=p.mass_squared*(q1*q1+q2*q2)
        gradient=(3*q1*q1+8*q2*q2)/a**2
        density=(kinetic+mass+gradient)/(2*VOLUME)
        pressure=(kinetic-mass-gradient/3)/(2*VOLUME)
        volume=VOLUME*a**3;work=-pressure*3*h*volume
        field_energy=energy(z,p);derived_work=proper_rhs(z,p)[9]
        rows.append(dict(Ht=u,field_energy_over_epsilon_squared=field_energy,
            averaged_spatial_density_over_epsilon_squared=density,
            averaged_spatial_pressure_over_epsilon_squared=pressure,
            pressure_work_rate_over_epsilon_squared=work,
            direct_field_work_rate_over_epsilon_squared=derived_work,
            work_identity_absolute_residual=abs(work-derived_work)))
    nu=math.sqrt(source.nu_squared)
    probe_frequency=math.sqrt(p.mass_squared);gamma=3*p.H
    susceptibility=1/complex(p.mass_squared-probe_frequency**2,-gamma*probe_frequency)
    phase_lag=math.atan2(susceptibility.imag,susceptibility.real)
    absorbed=probe_frequency*susceptibility.imag/2
    damping_work=gamma*probe_frequency**2*abs(susceptibility)**2/2
    return dict(samples=rows,
        closed_parent_first_law='dE_chi=-p_chi dV for spatial averages; E=V*rho, V=2*pi^2*a^3',
        retained_patch_first_law='dE_K/dt=-F_K across the moving observer horizon; other region energy definitions require their work/flux terms',
        envelope_rate_per_Ht=1.5,late_angular_rate_per_Ht=nu,
        late_spiral_radius_ratio_per_turn=math.exp(-3*math.pi/nu),
        late_spiral_inverse_amplitude_ratio_per_turn=math.exp(3*math.pi/nu),
        spatial_wavelength_ratio_per_turn=math.exp(2*math.pi/nu),
        phase_plane='(q,(qdot+3h*q/2)/(H*nu)) asymptotically follows an inward logarithmic spiral',
        stereographic_sphere='(2r cosphi,2r sinphi,1-r^2)/(1+r^2)',
        spiral_is_a_mode_state_representation_not_a_derived_horizon=True,
        phase_shift_alone_is_amplitude_decay=False,
        late_retarded_response=dict(
            operator='qddot+3H qdot+mu^2 q=F0 cos(omega t), a weak diagnostic drive on the late operator',
            susceptibility='chi_R(omega)=1/(mu^2-omega^2-3iH omega), Fourier convention exp(-iomega t)',
            probe_frequency=probe_frequency,phase_lag_radians=phase_lag,
            lag_as_proper_duration_readout=phase_lag/probe_frequency,
            mean_absorbed_power_over_F0_squared=absorbed,
            mean_damping_work_over_F0_squared=damping_work,
            power_balance_absolute_residual=abs(absorbed-damping_work),
            isolated_constant_phase_shift_is_this_response_lag=False,
            drive_is_a_probe_not_a_refitted_source=True,
            charge_conjugate_responses_differ_at_zero_internal_potential=False),
        signal_zero_is_zero_energy=False,
        fine_grained_irreversibility_proved_by_coordinate_damping=False)


def spectrum_record():
    return dict(spatial_operator='-Delta_S3 has eigenvalues l(l+2)/a^2 and degeneracy (l+1)^2',
        dispersion='omega_l^2=mu^2+l(l+2)/a^2 for a local slowly varying background',
        lowest_varying_mode='l=1: lambda_spatial=2*pi*a/sqrt(3); the uniform l=0 mode also exists',
        highest_l_selected=False,Planck_length_proved_as_hard_wavelength_cutoff=False,
        thermal_zero_selects_largest_wavelength=False,
        classical_energy_with_cutoff='E=(kB*T)*(L+1)*(L+2)*(2L+3)/6',
        continuum_classical_mode_sum_converges_without_a_state_or_cutoff_rule=False,
        thermal_or_quantum_state_law_derived_from_generator=False,
        interpretation='Geometry fixes the mode spectrum; occupation, ultraviolet completion and physical ruler remain to be derived. A resolution representation does not establish a literal simulator.')


def ccc_record(parent):
    H=parent.H;m2=parent.mass_squared/H**2
    return dict(
        horizon_zero=dict(equation='v=eta+chi=pi/2',normal_norm=0.,causal_type='null',surface_tension=0.,
            signal_flow_balance='z_H=1-H*R; radial null rates in outgoing PG form are dR/dt=H*R+/-1',
            zero_signal_balance_is_zero_energy_density=False,
            endpoint='One north-pole point on the spacelike future conformal boundary'),
        conformal_zero=dict(equation='Omega=cos eta=sech(Ht)=0',normal_norm=-H*H,
            causal_type='spacelike',rescaled_metric='gbar=Omega^2*g=H^-2*(-deta^2+dOmega3^2)',
            same_surface_as_observer_horizon=False),
        reciprocal_comparison=dict(
            rule='a_new=Lx^2/a_old; declared benchmark Lx=1/H, not an action-selected crossover normalization',
            Lx=1/H,H_squared_a_old_a_new=1.,
            conformal_new_scale='a_new(tau)=sin(tau)/H after a separately specified expanding extension',
            proper_new_scale='T=(1-cos(tau))/H; a_new(T)=sqrt(2T/H-T^2)',
            radiation_density_times_a_fourth=3/(parent.kappa*H*H),
            dimensionless_kappa_H_squared_rho_a_fourth=3.,
            old_evolution_Ricci_scalar=12*H*H,new_evolution_Ricci_scalar=0.,
            new_vacuum_from_old_radiation_absent=0.,
            new_radiation_Ricci_squared='12/(H^4*a_new^8)',
            new_radiation_Kretschmann='24/(H^4*a_new^8)',
            physical_bang_regular=False,conformal_cylinder_regular=True,
            new_maximum_scale=1/H,new_proper_age_at_maximum=1/H,
            new_proper_age_at_recollapse=2/H,
            same_action_matter_continuation=False),
        fluid_duality=[reciprocal_fluid(w) for w in (-1.,0.,1/3)],
        actual_scalar=dict(m_squared_over_H_squared=m2,
            conformal_singular_coefficient=m2-2.,
            generic_mode_asymptotic='q/Omega=Omega^(1/2)*[A cos(nu ln Omega)+B sin(nu ln Omega)]+subleading',
            full_field_weight_one_regularity='C0 zero generically, not C1',
            horizon_restricted_extra_node='Y1,Y2 carry an extra Omega along this horizon; its faster tail does not establish full-boundary regularity',
            minimal_mass_fade_alone_makes_conformal_operator_regular=False,
            missing='An action-and-state crossover rule, with a conformally regular operator and inherited state data; background reciprocity alone does not supply it.'),
        regular_operator_comparisons=dict(
            condition='m_bare^2/H^2+12*xi-2=0 on the fixed de Sitter tail',
            minimal_critical_mass_squared=2*H*H,
            actual_mass_squared_over_critical=parent.mass_squared/(2*H*H),
            massless_conformal_xi=1/6,
            regular_real_unit_modes=[dict(l=l,frequency_in_eta=l+1,
                crossover_value=float(round(math.cos((l+1)*math.pi/2))),
                crossover_eta_derivative=float(round(-(l+1)*math.sin((l+1)*math.pi/2)))) for l in (1,2)],
            comparisons_are_actual_massive_state_conversion=False))


def winding_and_charge():
    mode=harmonic_sector();w=mode['charged_mass']
    W=2*s.sqrt(6)/(3*s.pi**2)*s.Matrix([-4*X[1]*X[3],2*X[0]*X[3],0,2*X[0]*X[1]])
    helicity=s.expand(sum(s.LeviCivita(a,b,c,d)*X[a]*W[b]*s.diff(W[d],X[c])
                         for a in range(4) for b in range(4) for c in range(4) for d in range(4)))
    if helicity!=0:raise ArithmeticError('the inherited circulation helicity control changed')
    normK=2*mode['internal_radius_squared']/(mode['d']+1)
    Z0=math.exp((mode['d']+2)*mode['s_star'])*normK/2
    scaled_g2=1/Z0
    return dict(
        inherited_KK_connection=dict(
            generator='2*pi-period internal rotation in the selected (n1,n2) plane',
            mean_Killing_norm_squared_over_Lunit_squared=normK,
            mean_Killing_norm_squared_exact='88/65',
            g_grav_times_gauge_kinetic_coefficient=Z0,
            gauge_coupling_squared_over_g_grav=scaled_g2,
            gauge_coupling_over_sqrt_g_grav=math.sqrt(scaled_g2),
            gravitational_scale_parameter='g_grav=kappa_nat/L_unit^2 remains unselected; arithmetic kappa=1 is not g_grav=1',
            absolute_physical_gauge_coupling_selected=False,
            gauge_coupling_squared_over_kappa_mass_increment=scaled_g2/mode['KK_mass_squared_increment'],
            exact_gauge_gap_ratio='(d+1)/d=13/12',
            volume_log_derivative_gauge_squared=-(mode['d']+2),
            volume_log_derivative_mass_increment_squared=-(mode['d']+2),
            reference_s_zero_to_minimum_squared_quantity_ratio=math.exp(-(mode['d']+2)*mode['s_star']),
            reference_s_zero_to_minimum_gauge_ratio=math.exp(-(mode['d']+2)*mode['s_star']/2),
            ratios_select_the_volume_history=False,
            gauge_derivative='D_mu z=(partial_mu-i q A_mu)z, q=+/-1; A_can=A/g_KK at fixed s',
            kinetic_term='-exp((d+2)s)<K^2> F^2/(8 kappa_nat); retain field-dependent coefficient if s varies',
            neutral_preparation_gauge_invariant_horizon_potential=0.,
            retained_mean_internal_current=0.,
            charge_resolved_cofrequency='omega-q Phi_K; Phi_K is a horizon/reference potential difference or inherited holonomy',
            geometry_alone_selects_nonzero_Phi_K=False,
            identified_as_electromagnetism_or_membrane_charge=False,
            full_nonlinear_finite_KK_truncation_proved=False),
        state_charge=dict(region_convention_is_charge_conjugation=False,
            neutral_state_expected_internal_charge=0.,
            charge_blind_positive_flow_selects_nonzero_mean=False,
            symmetry='[C,P_region]=0; C rho C^-1=rho; C Q C^-1=-Q implies <Q_region>=0',
            physical_baryon_photon_eta=None,
            conditional_eta_if_identified_with_baryons_and_positive_photons=0.,
            target_eta=6.1e-10,conditional_miss_percent=-100.,
            neutral_pair_controls=[dict(squeeze=x,**pair_transport(0.,0.,x)) for x in (.1,1.)],
            same_real_circulation_helicity=str(helicity),
            same_real_circulation_mean_norm_squared=str(mean(W.dot(W))),
            missing_selector='A charge-sensitive inherited connection or non-C-invariant state; baryon/photon currents and conversion must also be specified.'),
        reversed_winding=dict(model_mass=w,spin='scalar',charge='internal KK rotation, not identified as electric or baryon charge',
            comparison_scope='Locally normalized stationary oscillator; expanding-volume charge and energy use their own factors',
            convention='Q=-i(z* zdot-z zdot*); E=|zdot|^2+omega^2|z|^2',
            opposite_winding_Q_per_A_squared=[2*w,-2*w],
            each_positive_E_per_A_squared=2*w*w,
            energy_per_absolute_charge=w,
            fixed_charge_energy_bound='E>=omega*abs(Q); constant-radius circular motion saturates the bound in the stationary normalized mode',
            standing_wave_Q=0.,standing_wave_E_per_A_squared=4*w*w,
            amplitude_sign_flip_preserves_charge=True,complex_conjugation_reverses_charge=True,
            quantum_action_normalization_imported_for_E_equals_hbar_omega=True,
            conjugating_internal_phase_is_not_a_general_time_reversal_identification=True,
            annihilation_or_photons_derived_from_this_free_pair=False),
        two_photon_kinematics=dict(
            supplied_channel='particle+antiparticle at rest in their COM -> exactly two massless quanta',
            each_photon_energy_over_mass_c_squared=1.,cos_opening_angle=-1.,
            annihilation_interaction_and_rate_derived=False,
            electron_mass_spin_and_SI_scale_derived=False))


def pet_comparison_only():
    # 2022 CODATA reference value, not used in the NSC source or crossover.
    electron_MeV=.51099895069;e=1.602176634e-19;h=6.62607015e-34;c=299792458.
    E=electron_MeV*1e6*e;f=E/h
    return dict(reference_only=True,NSC_prediction=False,
        reference_source='https://physics.nist.gov/cuu/pdf/wall_2022.pdf',
        photon_reference_energy_keV=electron_MeV*1000,
        photon_reference_frequency_Hz=f,angular_frequency_rad_s=2*math.pi*f,
        photon_reference_wavelength_m=c/f,
        example_200ps_contains_reference_carrier_cycles=200e-12*f,
        coincidence_relation='x=-(c/2)*(t1-t2), detector 1 at +D and detector 2 at -D',
        recorded_or_calibrated_observables=['electrical-pulse threshold timestamps','calibrated deposited-energy signal or time-over-threshold proxy','reconstructed hit position','coincidence arrival differences'],
        exact_511_and_180_for_every_real_detector_hit=False,
        stored_time_release_is_a_measured_observable=False,
        measured_mass_referenced_phase_is_direct_carrier_cycle_count=False)


def calculation():
    source=solve_source();p=source.parent
    rows=[tail_state(source,u) for u in (4.,6.,8.)]
    q0=p.initial()[[2,4]];b=np.array([3.75,8.75]);nu2=source.nu_squared
    E0=(nu2+b)*q0*q0/2
    lower=E0*nu2/(nu2+b)
    return dict(schema='nsc-atlas-zero-interface-v1',
        convention=dict(plus='energy/source side',minus='counter/receiving side; negative-pressure interpretation is a stress test',
            zero='smooth interface/horizon label; the causal and conformal zero sets are specified separately',
            oscilloscope_zero='q=0 baseline of a signed oscillation; its two lobes can both carry positive energy',
            oscilloscope_baseline_is_already_a_geometric_horizon=False,
            positive_flow='plus -> minus',negative_sign_is_negative_energy=False,
            spatial_description='Three generated spatial directions with offset-parametrized evolution; proper duration is a readout',
            coordinate_description='A covariant representation can use three spatial coordinate arguments plus one offset argument; no fourth spatial axis is added'),
        tail=dict(H=p.H,vacuum_density=p.vacuum,scalar_mass_squared=p.mass_squared,
            rescaled_frequency_nu=math.sqrt(nu2),proper_carrier_frequency_limit=p.H*math.sqrt(nu2),
            scalar_amplitude_exponent_in_a=-1.5,averaged_excitation_density_exponent_in_a=-3.,
            spatial_wavelength_exponent_in_a=1.,conditional_free_radiation_density_exponent_in_a=-4.,
            total_density_limit=p.vacuum,usable_excitation_density_limit=0.,
            tail_sample_offsets_select_a_physical_epoch=False,
            density_estimate_definition='Late adiabatic phase-average estimate using finite-u rescaled oscillator energy; finite-u corrections are not included',
            nonzero_rescaled_energy_lower_bounds=lower.tolist(),samples=rows,
            quantum_Gibbons_Hawking_temperature_if_Euclidean_state=p.H/(2*math.pi)),
        conformal_break=ccc_record(p),charge_and_winding=winding_and_charge(),
        thermodynamics=thermodynamic_record(source),spectrum=spectrum_record(),
        counter_pressure=dict(density=p.vacuum,pressure=-p.vacuum,w=-1.,
            active_gravitational_density=-2*p.vacuum,background_acceleration_over_scale=p.H**2,
            source='Existing positive vacuum minimum, present on both sides of the smooth same-action horizon',
            pressure_sign_derived_from_region_sign=False,
            late_Omega_vacuum=1.,comparison_target_Omega_Lambda=.685,
            if_late_limit_misidentified_as_today_relative_difference_percent=100*(1/.685-1),
            comparison_is_cross_epoch_not_a_current_data_fit=True,
            observed_today_epoch_and_share_derived=False,
            flowing_null_stress_is_pure_vacuum_stress=False),
        detector_comparison=pet_comparison_only(),
        no_observational_fit=True,
        next_connection='Derive an inherited coupling between the spacetime flow orientation and internal winding charge, rather than identifying those signs by notation.')
