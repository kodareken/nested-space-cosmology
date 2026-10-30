"""Three target audit using only retained dimensionless atlas equations.

No SI constants, target-selected depth, state imbalance or component weights
enter a prediction. Conditional controls and inverse diagnostics are explicit.
"""
import math
from fractions import Fraction as F


def room_band(lapse,shift,radius,spacing):
    if not all(math.isfinite(x) for x in (lapse,shift,radius,spacing)) or min(lapse,radius,spacing)<=0:
        raise ValueError('finite metric data and positive lapse/radius/spacing required')
    return (radius/spacing)*(1+abs(shift)/lapse)


def sphere_band(highest_mode,current_over_inflation_radius):
    if type(highest_mode) is not int or highest_mode<1 or current_over_inflation_radius<=0:
        raise ValueError('positive mode label and radius ratio required')
    return highest_mode*current_over_inflation_radius


def orientation_imbalance(n_plus,n_minus):
    if min(n_plus,n_minus)<0 or n_plus+n_minus<=0:raise ValueError('nonempty nonnegative populations required')
    return (n_plus-n_minus)/(n_plus+n_minus)


def outside_energy_share(outside_over_inside_volume):
    v=outside_over_inside_volume
    if v<0:raise ValueError('nonnegative volume ratio required')
    return 3*v/(1+3*v)


def calculation(omega,late_vacuum,occupied_level,targets=(136.84,138.26)):
    logomega=math.log(omega)
    square=F(57600,39169);R=math.sqrt(float(square))
    rows=[]
    for depth in (1,2,100):
        value=depth*logomega
        rows.append(dict(assumed_depth=depth,log_band_if_endpoint_factor_is_one=value,
            log_miss_percent=[100*(value/t-1) for t in targets],
            selected_by_geometry=False,status='Conditional arithmetic only; neither depth nor endpoint factor is selected.'))
    return dict(schema='nsc-atlas-three-targets-v1',
        inputs=dict(omega=omega,late_vacuum_Lunit_squared=late_vacuum,occupied_wall_reference_level=occupied_level,
                    Parke_R_squared=str(square),Parke_R_over_Lunit=R,
                    late_vacuum_L_Parke_squared=late_vacuum*float(square)),
        band=dict(formula='ln(B)=ln(q_max)+ln(R_now/R_infl), for common local c and n_min=1',
            varying_speed_formula='B=q_max*(c_infl/c_now)*(R_now/R_infl)',
            wall_product='omega_max*R_infl/c=q_max; the untruncated n/R spectrum has no largest n',
            ground_wavelength_over_radius='2*pi from lambda=2*pi*c/omega_1 and omega_1=c/R',
            same_radius_occupied_reference_log_band=math.log(occupied_level),
            RN_PG=dict(formula='B=(r/Delta_r)*(1+abs(beta)/N)',N_meaning='lapse, not an integer mode count',
                Delta_r_meaning='numerical spacing in the cited CFL diagnostic',
                refinement_control=[room_band(1.,0.,R,R/32),room_band(1.,0.,R,R/64)],
                refinement_log_change=math.log(2),physical_cutoff_selected=False),
            inheritance=dict(THEORY_Lambda_meaning='frequency scale, not cosmological curvature',
                frequency_step=omega,length_step_at_fixed_c=1/omega,curvature_step_for_equivalent_units=omega**2,
                requested_sqrt_length_step=1/math.sqrt(omega),
                sqrt_step_requires='Reinterpreting Omega as a curvature factor; this is not the frequency-scale definition in THEORY.',
                conditional_rows=rows,
                curvature_interpretation_control=dict(assumed_depth=200,
                    log_band_if_endpoint_factor_is_one=100*logomega,selected_by_geometry=False),
                inverse_depth_interval=[t/logomega for t in targets],
                inverse_depth_interval_if_Omega_is_curvature=[2*t/logomega for t in targets],
                depth_selected=False,candidate_identified_with_Parke_volume_model=False),
            target_log_interval=list(targets),prediction=None,percent_miss=None,
            missing='A state-selected physical endpoint/transition map. The current formulas leave q_max, the inflation/current radius ratio and any cosmological depth independent.'),
        asymmetry=dict(orientation_fraction='A=(n_plus-n_minus)/(n_plus+n_minus)',
            baryon_photon_eta='eta_B=A_b*(n_b+n_antib)/n_gamma',
            phase_alone_changes_occupations=False,
            symmetric_orientation_value=orientation_imbalance(.5,.5),
            alternative_same_operator_population_value=orientation_imbalance(.75,.25),
            target=6.1e-10,prediction=None,percent_miss=None,
            missing='Inherited charged-state preparation and its baryon/photon identification. A relative phase fixes neither the population difference nor photon normalization.',
            scope='The two population examples are controls of the given conjugate-mode algebra, not baryogenesis solutions.'),
        split=dict(formula='Omega_a=kappa*rho_a/(3*H^2), evaluated by one observer from one effective stress tensor',
            auxiliary_outside_fraction='f_out=3*v/(1+3*v), v=V_out/V_in for the stipulated densities',
            volume_controls=[dict(volume_ratio=v,outside_share=outside_energy_share(v)) for v in (1.,1/3)],
            targets=dict(Omega_Lambda=.685,Omega_c=.265,Omega_b=.049),
            predictions=dict(Omega_Lambda=None,Omega_c=None,Omega_b=None),
            percent_miss=dict(Omega_Lambda=None,Omega_c=None,Omega_b=None),
            pure_vacuum_limit_control=dict(Omega_Lambda=1.,Omega_c=0.,Omega_b=0.,
                percent_miss=dict(Omega_Lambda=100*(1/.685-1),Omega_c=-100.,Omega_b=-100.),
                status='Late vacuum-only Friedmann control, not present-day abundance or a baryon-production calculation.'),
            missing='A state-normalized exterior-to-observer stress projection with vacuum, clustering and baryon sectors. A total outside energy share does not supply that decomposition or H.'),
        no_CODATA_or_SI_inputs=True,physical_targets_used_for_selection=False,
        new_field_evolutions=0,new_grids=0,new_likelihood_fits=0,
        stopping_result='All three target predictions remain unselected by the permitted equations; conditional arithmetic is retained without promoting it to prediction.')
