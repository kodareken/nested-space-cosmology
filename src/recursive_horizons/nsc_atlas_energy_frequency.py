"""Conditional links between mode energy, pressure, radius and observations.

The occupied surface action, the volume-scalar action and the radiation
benchmark remain distinct. No depletion rate, thermal state or photon
completion is silently added to the fixed cosmological candidate.
"""
import math
from .nsc_atlas_occupied_wall import filled_shell


def flat_equilibrium_ratios(q, reference_q=1000, species=2):
    """Equal bulk pressure, negligible gravity, fixed positive bare tension."""
    if min(q,reference_q)<1:raise ValueError('positive filled levels required')
    a=filled_shell(q,species)['energy_coefficient']
    a0=filled_shell(reference_q,species)['energy_coefficient']
    r=(a/a0)**(1/3)
    return dict(level=q,radius_over_reference=r,total_energy_over_reference=r*r,
        fixed_level_frequency_over_reference=1/r,
        highest_occupied_frequency_over_reference=(q/reference_q)/r,
        highest_frequency_cubed_over_bound=2*q*q/((q+1)*(2*q+1)))


def radiation_temperature_ratio(radius_ratio,energy_ratio):
    """Homogeneous equilibrium radiation of fixed species content only."""
    if min(radius_ratio,energy_ratio)<=0:raise ValueError('positive ratios required')
    return (energy_ratio/radius_ratio**3)**.25


def thermal_log_rate(R,Rdot,E,net_outward_power):
    if min(R,E)<=0:raise ValueError('positive radius and thermal energy required')
    return -Rdot/R-net_outward_power/(4*E)


def clock_redshift(delta_N,delta_s):
    """Conditional null photons and identical clocks scaling as exp(-6s)."""
    return math.exp(delta_N-6*delta_s)


def calculation(occupied_record,volume_record):
    rows={r['level']:r for r in occupied_record['states']}
    before,after=rows[1000],rows[999]
    Rratio=after['radius']/before['radius']
    aratio=after['energy_coefficient']/before['energy_coefficient']
    v=volume_record['averaged_evolution']
    ds=v['s'];N=v['endpoint_N'];frequency=math.exp(-6*ds)
    return dict(schema='nsc-atlas-energy-frequency-v1',
        flat_surface_theorem=dict(
            mode_energy='hbar*c*A_q/R',A_q='N_species*q*(q+1)*(2q+1)/3',
            total_energy='4*pi*sigma0*R^2+hbar*c*A_q/R',
            equilibrium='R^3=hbar*c*A_q/(8*pi*sigma0)',
            total_energy_at_equilibrium='12*pi*sigma0*R^2',
            second_radius_derivative_at_fixed_occupation='24*pi*sigma0 > 0',
            omega_max_cubed_bound='12*pi*sigma0*c^2/(N_species*hbar)',
            bound_ratio_exact='omega_max^3/bound=2*q^2/((q+1)*(2q+1)) < 1',
            frequency_convention='omega is angular frequency; nu_Hz=omega/(2*pi)',
            scope='Filled positive Dirac levels, fixed tension, equal bulk pressures, negligible gravity. This bounds the highest occupied level of this family, not all allowed modes or a thermal temperature.',
            examples=[flat_equilibrium_ratios(q) for q in (999,1000)]),
        retained_gravitating_endpoint_comparison=dict(
            levels=[1000,999],radius_ratio=Rratio,
            geometric_mass_ratio=after['mass']/before['mass'],
            surface_mode_energy_ratio=aratio/Rratio,
            fixed_level_frequency_ratio=1/Rratio,
            highest_occupied_frequency_ratio=(999/1000)/Rratio,
            radius_percent_change=100*(Rratio-1),
            mode_energy_percent_change=100*(aratio/Rratio-1),
            highest_frequency_percent_change=100*((999/1000)/Rratio-1),
            scope='Two retained stationary states under one action. No leakage channel or physical trajectory between them has been constructed.'),
        static_pattern=dict(identity='rho=exp(-i*H*t)*rho0*exp(i*H*t)=rho0 if [H,rho0]=0',
            scope='The filled spectral projector has stationary mean stress; not a claim that every real standing wave has static local energy density.'),
        thermal_benchmark=dict(
            premises='Homogeneous equilibrium radiation, fixed species, E=a_rad*T^4*V and p=E/(3V)',
            first_law='dE/dt+p*dV/dt=P_in-P_out',
            log_temperature_rate='Tdot/T=-Rdot/R-(P_out-P_in)/(4E)',
            integrated_relation='T*R/(T0*R0)=exp(-integral[(P_out-P_in)/(4E)]dt)',
            cooling_while_shrinking='P_out-P_in > 4E*(-Rdot/R)',
            initial_temperature_ceiling_condition='The logarithmic temperature derivative must remain nonpositive for a monotonic global ceiling; compression or heating can violate it.',
            half_radius_half_energy_T_ratio=radiation_temperature_ratio(.5,.5),
            half_radius_one_32_energy_T_ratio=radiation_temperature_ratio(.5,1/32),
            no_flux_adiabatic_half_radius_T_ratio=2.,
            imported_thermal_sector=True,thermal_sector_derived_from_volume_action=False),
        retained_volume_example=dict(
            N=N,s=ds,matter_density_ratio=v['matter_density']/(3/8),
            frequency_ratio=frequency,internal_volume_ratio=math.exp(12*ds),
            identity='rho/rho0=exp(-3N)*(omega/omega0)',
            total_homogeneous_energy_law='d(rho_total)/dt=-3H*(K*sdot^2+E_matter)',
            temperature_defined=False,exterior_leakage_derived=False),
        observer_protocol=dict(
            premises='Transparent null probe on homogeneous metric; identical emission/detection clocks with nu_clock,E proportional exp(-6s)',
            redshift='1+z=exp(DeltaN-6*Delta_s)=(a_E*Omega)_o/(a_E*Omega)_e=a_J,o/a_J,e',
            retained_endpoint_one_plus_z_if_this_protocol_is_used=clock_redshift(N,ds),
            clock_calibrated_contraction_gives_blueshift=True,
            enclosing_wall_radius_is_same_as_internal_scale_factor=False,
            microscopic_photon_detector_completion_derived=False,
            old_likelihoods_recomputed=False),
        conclusion='Radius, spectral state and energy are linked. A universal inflation-temperature ceiling, cooling by deflation alone and a complete shrinking-wall/expanding-interior history do not follow.',
        new_evolutions=0,new_grids=0,new_fits=0)
