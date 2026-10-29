"""Exact wall projection and the missing HM-to-preparation connection."""
from fractions import Fraction as F
import math
from scipy.optimize import brentq
from .nsc_atlas_inherited_matter import exact_initial_geometry
from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_shared_scalar_slice import SharedScalarSlice


def wall_projection(radius=F(19,10),mass=F(3,25)):
    R,M=F(radius),F(mass);k=F(1,10);Ui,Uo,E0=F(7,32),F(21,32),F(3,8)
    beta=-M/(k*R**2)+(Uo-Ui+3*k*k)*R/(6*k)
    fm=1-2*M/R-Ui*R*R/3
    H2=(beta*beta-fm)/(R*R);C=Uo/3-H2
    assert beta-k*R>0 and H2>0
    return dict(radius=str(R),mass=str(M),beta_minus=float(beta),beta_plus=float(beta-k*R),
                C_exact=str(C),C=float(C),C_over_E0_exact=str(C/E0),C_over_E0=float(C/E0),
                H_squared=float(H2),Omega_k_initial=float(-C/H2))


def calculation():
    reference=wall_projection()
    zero=brentq(lambda r:wall_projection(F(str(r)))['C'],1.9,2.,xtol=1e-14)
    witnesses=[]
    for R in (1.9,zero,2.):
        r=wall_projection(F(str(R)))
        if F(str(R)).denominator>10000:
            r['C_exact']=None;r['C_over_E0_exact']=None
            r['radius_is_numerical_zero_condition']=True
        w=SharedScalarSlice(outer_radius=R).wall(1)
        r['independent_angular_residual']=w['angular_residual']
        r['independent_velocity_squared_error']=w['velocity']**2-r['H_squared']*R*R
        witnesses.append(r)
    m=VolumeParent();top,minimum=m.vacuum_roots()
    beta=-3*top['stiffness']/(m.K*top['vacuum']);growth=(math.sqrt(9+4*beta)-3)/2
    return dict(schema='nsc-atlas-wall-exit-v1',reference_projection=reference,witnesses=witnesses,
        radius_for_C_zero=zero,radius_is_inherited=False,
        density_dictionary=dict(displacement_fluxes=['3/2','1','1/2','0'],ground_vacua=['9/32','1/8','1/32','0'],
            actual_initial_middle='7/32',actual_initial_parent='21/32',parent_excitation='3/8'),
        volume_wall=dict(K=84,jump_at_matched_start='-1/70',
            tension_at_HM_location=.2*math.exp(-6*top['s']),
            normal_derivative_jump_at_HM_location=-.2*math.exp(-6*top['s'])/14,
            missing='A compatible nonuniform volume collar and its gradient stress; the old projection alone is not the full (12,2) junction.'),
        HM=dict(s=top['s'],V=top['vacuum'],vacuum_excitation=0,
            exact_stationary_continuation=True,unique_exit_efolds=None,
            linear_growth_per_efold=growth,efold_change_per_decade_seed=math.log(10)/growth,
            averaged_acceleration_with_loaded_E0=3/112,
            loaded_branch='rightward, toward decompactification; not the stable minimum on the left',
            branch_proof='For s>=s_top, V_s<=0; on v=0 the source force is nonnegative, so a loaded zero-velocity state cannot cross into v<0 and reach the left minimum.',
            energy_if_old_action_preserved=.375*math.exp(-6*top['s']),
            energy_if_old_coherent_field_preserved=.375*math.exp(-12*top['s'])),
        minimum_s=minimum['s'],post_exit_C_over_E0_derived=False,
        post_exit_observer_curvature_derived=False,
        route_dispositions=[
            'Old Israel projection: reproduces C/E0 using supplied R and M; fixed tension and vacua alone do not select them.',
            'Full volume junction: the required normal-volume jump is not supplied by a uniform initial state.',
            'Open CDL child: no compatible (12,2) CDL exit was computed; the open geometry cannot be attached merely by choosing its sign.',
            'Exact HM continuation: remains at the barrier without a specified perturbation/state; no unique rolling e-fold count follows.',
            'HM location plus E0=3/8: is a distinct specified matter preparation and initially accelerates toward the runaway.'])
