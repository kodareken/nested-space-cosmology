"""SI hierarchy audit: forward identities and explicitly inverse diagnostics.

No value of g or tunnelling-to-ruler identification is selected here.
The source coefficients are retained (12,2) results, not new field runs.
"""
import math
from .nsc_atlas_scale_selection import NEWTON_G,PLANCK_H,LIGHT_SPEED

A0=6.626077076099598
VSTAR=.22098302421260696
B_HM=403.7817410357339
E_HAT=.5
SIGMA_HAT=.2
TARGET_HIERARCHY=1.661269153e60
TARGET_LAMBDA=1.0904e-52


def units():
    kappa=8*math.pi*(PLANCK_H/(2*math.pi))*NEWTON_G/LIGHT_SPEED**3
    ell=E_HAT/SIGMA_HAT*math.sqrt(kappa)
    return dict(kappa_nat_m2=kappa,ell_e_sigma_m=ell)


def conditional_at_g(g):
    """A dimensional completion at a supplied g, not an inherited prediction."""
    if not math.isfinite(g) or g<=0:raise ValueError('positive finite supplied g required')
    u=units();k=u['kappa_nat_m2'];L=math.sqrt(k/g)
    hierarchy=A0*SIGMA_HAT/E_HAT/math.sqrt(g)
    return dict(supplied_g=g,L_unit_m=L,L_ref_m=A0*L,hierarchy=hierarchy,
        log_hierarchy=math.log(A0*SIGMA_HAT/E_HAT)-.5*math.log(g),
        Lambda_star_m_inverse2=VSTAR*g/k,
        Lambda_star_ell_squared=VSTAR*(E_HAT/SIGMA_HAT)**2*g,
        B_HM=B_HM/g,eta_e0=E_HAT*math.sqrt(g),
        canonical_charge=E_HAT/(math.sqrt(k)*L),canonical_tension=SIGMA_HAT/(k*L))


def calculation(*,target_hierarchy=TARGET_HIERARCHY,target_lambda=TARGET_LAMBDA):
    if target_hierarchy<=1 or target_lambda<=0:raise ValueError('positive targets required for inverse diagnostics')
    u=units();logtarget=math.log(target_hierarchy)
    product=target_lambda*u['ell_e_sigma_m']**2
    rows=[]
    for p in (1,2,4):
        inverse_g=B_HM/(p*logtarget)
        actual=conditional_at_g(inverse_g)
        rows.append(dict(imposed_exponential=f'exp(B/{p})',p=p,
            required_B_using_displayed_log=p*138.66,
            required_g_using_displayed_log=B_HM/(p*138.66),
            required_B_using_full_target_log=p*logtarget,required_g_using_full_target_log=inverse_g,
            same_geometry=actual,
            same_geometry_hierarchy_over_target=actual['hierarchy']/target_hierarchy,
            same_geometry_hierarchy_miss_percent=100*(actual['hierarchy']/target_hierarchy-1),
            same_geometry_log_miss_percent=100*(actual['log_hierarchy']/logtarget-1),
            same_geometry_Lambda_over_observed=actual['Lambda_star_ell_squared']/product,
            same_geometry_Lambda_miss_percent=100*(actual['Lambda_star_ell_squared']/product-1),
            status='Target-dependent inversion of an imposed exponential; not a selected coupling or prediction.'))
    factor=A0*SIGMA_HAT/E_HAT
    g_from_hierarchy=(factor/target_hierarchy)**2
    at_hierarchy=conditional_at_g(g_from_hierarchy)
    g_from_lambda=target_lambda*u['kappa_nat_m2']/VSTAR
    return dict(schema='nsc-atlas-si-hierarchy-audit-v1',
        fixed=dict(d=12,j=2,a0=A0,Vstar=VSTAR,HM_coefficient=B_HM,charge_hat=E_HAT,tension_hat=SIGMA_HAT),
        units=u,targets=dict(hierarchy=target_hierarchy,log_hierarchy=logtarget,Lambda_obs_m_inverse2=target_lambda,
                            Lambda_obs_ell_squared=product),
        forward_identities=dict(hierarchy_numerator=factor,log_hierarchy_intercept=math.log(factor),
            Lambda_ell_squared_coefficient=VSTAR*(E_HAT/SIGMA_HAT)**2,
            hierarchy_squared_times_Lambda_ell_squared=A0*A0*VSTAR,
            log_hierarchy='ln(a0*sigma_hat/e_hat)-0.5*ln(g)',
            Lambda_ell_squared='Vstar*(e_hat/sigma_hat)^2*g'),
        dimensional_completions=[conditional_at_g(g) for g in (1.,.25)],
        exponential_inverse_diagnostics=rows,
        target_inverse_consistency=dict(
            g_if_hierarchy_is_inserted=g_from_hierarchy,
            B_if_hierarchy_is_inserted=at_hierarchy['B_HM'],
            Lambda_ell_squared_if_hierarchy_is_inserted=at_hierarchy['Lambda_star_ell_squared'],
            Lambda_miss_percent_if_hierarchy_is_inserted=100*(at_hierarchy['Lambda_star_ell_squared']/product-1),
            g_if_Lambda_obs_is_inserted=g_from_lambda,
            hierarchy_if_Lambda_obs_is_inserted=conditional_at_g(g_from_lambda)['hierarchy'],
            log_hierarchy_if_Lambda_obs_is_inserted=conditional_at_g(g_from_lambda)['log_hierarchy'],
            target_pair_product=target_hierarchy**2*product,
            explanation='The quoted hierarchy was anchored to the old 9/32 vacuum, whose Lambda Lref^2 was 12.3482524. The (12,2) late coefficient is 9.7022370; the two supplied target values cannot both be obtained with it.'),
        tunnelling_dictionary=dict(rate_per_four_volume='Gamma4=A4*exp(-B), [A4]=length^-4',
            unit_expected_event_length='Gamma4^(-1/4)=A4^(-1/4)*exp(B/4)',
            fixed_three_volume_wait='tau=1/(Gamma4*V3)=exp(B)/(A4*V3)',
            half_exponent='Taking an amplitude or square root changes the statistic; no length identification follows by itself.',
            geometric_ruler_identification_derived=False,
            prefactor_length_identified_with_microanchor=False),
        selected_g=None,selected_hierarchy_exponent=None,predicted_log_hierarchy=None,
        predicted_Lambda_ell_squared=None,prediction_miss_percent=None,
        single_missing_numerical_input=dict(name='eta_e0',definition='kappa_nat*e_can at the declared reference slice',
            relation='eta_e0=0.5*sqrt(g), so g=4*eta_e0^2',selected=False),
        stopping_outcome='The geometry does not fix the absolute Planck-normalized membrane charge eta_e0. Stop here; do not continue to another route.',
        new_field_evolutions=0,new_observation_fits=0)
