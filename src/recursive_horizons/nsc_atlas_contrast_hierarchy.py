"""Energy contrast, gravitational length and positive-vacuum curvature scale.

All dimensional completions retain explicit g. Constant-s frame conversion
is distinct from a rolling spacetime horizon. Sector minima below are
conditional endpoints of one fixed-coupling potential, not a joined evolution.
"""
import math
from fractions import Fraction as F
from scipy.optimize import brentq
from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_parke_nucleation import thin_wall
from .nsc_atlas_si_hierarchy import units


def length_data(vacuum, s, g):
    """Lengths in L_unit(0), using LH=1/sqrt(Lambda), at constant s."""
    if vacuum<=0 or g<=0 or not all(math.isfinite(x) for x in (vacuum,s,g)):
        raise ValueError('positive vacuum and positive explicit g required')
    omega=math.exp(-6*s)
    lp=math.sqrt(g/(8*math.pi))
    lh=1/math.sqrt(vacuum)
    return dict(supplied_g=g,Planck_E=lp,Planck_J=omega*lp,
        L_H_E=lh,L_H_J=omega*lh,Lambda_J=vacuum/omega**2,
        kappa_J_over_kappa_E=omega**2,
        log_ratio=.5*math.log(8*math.pi/vacuum)-.5*math.log(g))


def sector_minimum(q, *, inherit_ancestor=True):
    """Relax q and s with UV bare term/curvature fixed, never rematched.

    inherit_ancestor=True implements the declared q_A=q+j*e, j=2 sheet
    prescription. False retains the parent's original sheet as a control.
    Neither constructs the sheet or a cross-wall transition.
    """
    q=F(str(q));base=VolumeParent().coefficients
    uv=base['required_UV_bare'];B=base['B']
    qa=q+1 if inherit_ancestor else base['ancestor_flux']
    A=uv+q*q/8;T=qa*qa/8
    if min(A,B,T)<=0:raise ValueError('declared positive coefficient domain required')
    af,bf,tf=map(float,(A,B,T))
    turn=(7*bf/(72*tf))**.2
    polynomial=lambda x:6*af-7*bf*x+12*tf*x**6
    if polynomial(turn)>=0 or polynomial(2*turn)<=0:
        raise ValueError('the declared finite minimum bracket failed')
    x=brentq(polynomial,turn,2*turn,xtol=5e-15)
    s=-math.log(x)/2
    # Use stationarity to reduce the cancellation in V at its minimum.
    vacuum=x**7*(bf/6-tf*x**5)
    stiffness=144*af*x**6-196*bf*x**7+576*tf*x**12
    ac=5*bf/6*(bf/(6*tf))**.2
    # Exact rational comparison with the zero-vacuum threshold.
    positive=(6*A/(5*B))**5>B/(6*T)
    return dict(q=str(q),ancestor_flux=str(qa),A=str(A),B=str(B),T=str(T),
        UV_bare=str(uv),s_min=s,x_min=x,Lambda_hat_min=vacuum,
        stiffness=stiffness,positive_vacuum_by_exact_comparison=positive,
        critical_A=ac,Jordan_Planck_factor=math.exp(-6*s),
        log_ratio_intercept=(.5*math.log(8*math.pi/vacuum) if positive else None),
        requested_positive_Lambda_length_exists=positive)


def calculation(target_lambda=1.0904e-52):
    u=units();kappa=u['kappa_nat_m2']
    lp=math.sqrt(kappa/(8*math.pi))
    observed=.5*math.log(8*math.pi/(kappa*target_lambda))
    bubble=thin_wall();delta=F(21,32)-F(7,32)
    recovered=3*(F(39169,57600)-(delta/(3*F(1,5))+F(1,5)/4)**2)
    snapshots=[]
    for label,q,rho in [('parent',F(3,2),F(21,32)),('child',F(1),F(7,32)),
                        ('grandchild',F(1,2),F(1,32))]:
        intercept=.5*math.log(8*math.pi/float(rho))
        snapshots.append(dict(region=label,flux=str(q),density=str(rho),
            log_ratio_intercept=intercept,
            log_miss_percent_if_g_is_set_to_one=100*(intercept/observed-1),
            size_ratio_to_observed_if_g_is_set_to_one=math.exp(intercept-observed)))
    relaxed=[dict(flux=str(q),vacuum=str(q*q/8),
                  log_ratio_intercept=.5*math.log(8*math.pi/float(q*q/8)))
             for q in (F(3,2),F(1),F(1,2))]
    sectors=[sector_minimum(q) for q in (1.5,1.,.5)]
    controls=[sector_minimum(q,inherit_ancestor=False) for q in (1.5,1.,.5)]
    minimum=sectors[0]
    fixed_g_control=length_data(minimum['Lambda_hat_min'],minimum['s_min'],1.)
    return dict(schema='nsc-atlas-contrast-hierarchy-v1',
        definitions=dict(L_H='1/sqrt(Lambda); positive-vacuum curvature length',
            de_Sitter_event_horizon='sqrt(3)*L_H',Planck_length_squared='kappa_nat/(8*pi)',
            Planck_scope='Gravitational scale from the Ricci coefficient, not a derived minimum length or a full scalar-mediated Cavendish measurement',
            g='kappa_E/L_unit(0)^2',selected_g=None),
        frame_law=dict(g_J='exp(-12s)*g_E',kappa_E='kappa_D/V0',
            kappa_J='exp(-12s)*kappa_E',Planck_J='exp(-6s)*Planck_E',
            Lambda_J='exp(12s)*Lambda_E',L_H_J='exp(-6s)*L_H_E',
            log_ratio='0.5*ln(8*pi/lambda_hat(s,q))-0.5*ln(g)',
            scope='Constant-s vacua; rolling curvature has derivative terms and the global horizon needs the full metric.'),
        positive_sector_contrast_law=dict(
            Planck_child_over_parent_J='exp(-6*Delta_s)',
            L_H_child_over_parent_J='exp(-6*Delta_s)*sqrt(Lambda_E_parent/Lambda_E_child)',
            log_ratio_gain='0.5*ln(Lambda_E_parent/Lambda_E_child)',
            scope='Positive constant-s vacuum endpoints in the same Einstein action and common absolute normalization; not all endpoints are realized by a transition.'),
        contrast=dict(Ui='7/32',Uo='21/32',sigma='1/5',kappa_hat=1,
            Parke_R_over_L_unit=bubble['radius'],recovered_child_Lambda_hat=str(recovered),
            closed_identity='Lambda_i=3*(R_c^-2-(DeltaU/(3sigma)+kappa*sigma/4)^2)=kappa*Ui',
            exterior_and_tension_select_radius_not_new_Newton_coupling=True,
            volume_endpoint_difference_selected=False,
            missing_relation='A compatible sourced volume profile and its state are needed for Delta s(Ui,Uo,sigma); Parke constant-vacuum data contain no such roll.'),
        observed_benchmark=dict(Lambda_m_inverse2=target_lambda,Planck_length_m=lp,
            curvature_length_m=1/math.sqrt(target_lambda),log_L_H_over_Planck=observed,
            log_de_Sitter_horizon_over_Planck=observed+.5*math.log(3),
            previous_numbers='138.2604 was minus ln(eta); 136.8439 was ln(L_P/ell_e_sigma). Neither is ln(L_H/l_P).'),
        auxiliary_Parke_snapshot_chain=snapshots,
        snapshot_growth=dict(parent_to_child='sqrt(3)',child_to_grandchild='sqrt(7)',
            total_log_growth=.5*math.log(21),
            scope='These are density-scale ratios for the stipulated constant-vacuum controls. Actual parent/child scalar snapshots are not stationary vacua.'),
        relaxed_scalar_without_volume_chain=relaxed,
        scalar_only_total_log_growth=math.log(3),
        fixed_coupling_j2_sectors=sectors,
        fixed_original_sheet_control=controls,
        exact_threshold=dict(x_critical='(B/(6*T))^(1/5)',A_critical='(5*B/6)*(B/(6*T))^(1/5)',
            identity_at_stationary_point='V_min=x^7*(B/6-T*x^5)',
            positive_iff='(6*A/(5*B))^5 > B/(6*T), for the right/minimum stationary root'),
        volume_parent_control=dict(data_if_g_is_set_to_one=fixed_g_control,
            log_miss_percent_if_g_is_set_to_one=100*(fixed_g_control['log_ratio']/observed-1),
            is_prediction=False),
        prediction=dict(log_hierarchy=None,percent_miss=None,
            reason='g and a compatible contrast-to-volume profile are not selected; the fixed-coupling j=2 relaxed child/grandchild minima are negative.'),
        new_field_evolutions=0,new_likelihood_fits=0,
        chain_scope='Conditional homogeneous vacuum sectors with the same UV bare term and internal curvature. No localized ancestor sheet, wall profile or causal transition is constructed.')
