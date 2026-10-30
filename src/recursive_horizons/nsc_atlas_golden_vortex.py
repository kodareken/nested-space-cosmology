"""Golden inward maps, source-stability controls and black-hole/BKL interfaces.

Geometry is constructed before a physical identification. Local radial
refinement, global angular sampling, geodesic instability and cosmological
expansion are retained as different readouts. No target selects a parameter.
"""
import math
import numpy as np
from scipy.spatial import SphericalVoronoi
from .nsc_atlas_pen_geometry import persistent_measure,sphere

PHI=(1+math.sqrt(5))/2
GOLDEN_ANGLE=2*math.pi/PHI**2
GOLDEN_PITCH=2*math.log(PHI)/math.pi


def suspend(point,offset):
    """The same set map B_d -> B_(d+1), used twice from a line segment."""
    x=np.asarray(point,dtype=float)
    if abs(offset)>1 or np.dot(x,x)>1+1e-13:raise ValueError('unit ball and offset in [-1,1] required')
    return np.r_[math.sqrt(max(0.,1-offset*offset))*x,offset]


def nested_ball_point(x,t1,t2,level,ratio=PHI**4,delta=GOLDEN_ANGLE):
    y=suspend(suspend([x],t1),t2)
    phase=math.remainder(level*delta,2*math.pi)
    c,s=math.cos(phase),math.sin(phase)
    return np.array([c*y[0]-s*y[1],s*y[0]+c*y[1],y[2]]),-level*math.log(ratio)


def vortex_point(turn,ratio=PHI**4,delta=GOLDEN_ANGLE):
    """Smooth inward curve; return direction and log radius, not an underflow."""
    return sphere(2*math.pi*turn,delta*turn),-turn*math.log(ratio)


def fibonacci_points(count):
    if type(count) is not int or count<4:raise ValueError('at least four points required')
    k=np.arange(count);z=1-2*(k+.5)/count;p=k*GOLDEN_ANGLE
    rho=np.sqrt(1-z*z)
    return np.c_[rho*np.cos(p),rho*np.sin(p),z]


def fibonacci_cover(count):
    p=fibonacci_points(count);vertices=SphericalVoronoi(p,radius=1,center=np.zeros(3)).vertices
    cosine=np.max(vertices@p.T,axis=1)
    actual=float(np.max(np.arccos(np.clip(cosine,-1,1))))
    return dict(points=count,cover_radius=actual,area_lower_bound=math.acos(1-2/count),
        log_angular_resolution=-math.log(actual),
        scope='Finite equal-area Fibonacci sample; its heights depend on count, so successive sets are not a persistent prefix.')


def kerr_control(spin):
    if not 0<=spin<1:raise ValueError('subextremal dimensionless spin required')
    rp=1+math.sqrt(1-spin*spin);A=rp*rp+spin*spin
    rows=[]
    for sign,name in ((1,'prograde'),(-1,'retrograde')):
        r=2*(1+math.cos(2*math.acos(-sign*spin)/3))
        circum=math.sqrt(r*r+spin*spin+2*spin*spin/r)
        exponent=math.sqrt(3)*(r-1)/(2*math.sqrt(r))
        rows.append(dict(orientation=name,photon_radius_over_mass=r,
            BL_photon_over_horizon=r/rp,circumferential_photon_over_horizon=circum/2,
            angular_radial_instability_exponent=exponent,
            log_radial_offset_factor_per_full_orbit=2*math.pi*exponent,
            radial_offset_factor_per_full_orbit=math.exp(2*math.pi*exponent)))
    return dict(spin=spin,horizon_BL_radius_over_mass=rp,horizon_areal_radius_over_mass=math.sqrt(A),
        horizon_area_over_4pi_mass_squared=A,horizon_angular_speed_times_mass=spin/A,
        equatorial_horizon_circumferential_radius_over_mass=2.,photon_orbits=rows)


def kasner(u):
    den=1+u+u*u
    return np.array([-u/den,(1+u)/den,u*(1+u)/den])


def kasner_bounce(p):
    a,b,c=p
    return np.array([-a,b+2*a,c+2*a])/(1+2*a)


def bkl_step(u):
    if u<=1:raise ValueError('nondegenerate u>1 required')
    return u-1 if u>=2 else 1/(u-1)


def contained_prefix(delta,ratio):
    """Analytic bound: convex k*delta+ratio^(-k-1) below first width."""
    first=1/ratio
    K=math.floor(first/delta)-1
    if K<1 or delta+ratio**-2>=first:return None
    margin=first-K*delta
    if margin<=0 or -(K+1)*math.log(ratio)>math.log(margin):return None
    return dict(contained_extra_bands=K,
        area_limit_lower=math.sin(first),
        log_additional_area_upper=-(K+2)*math.log(ratio)-math.log1p(-1/ratio),
        statement='All these bands lie in the first; the remaining area is bounded by the geometric tail. No tiny correction is rounded into an exact equality.')


def calculation(pen,parke,targets=(136.84,138.26)):
    omega=pen['construction']['omega'];delta=pen['construction']['delta'];R=parke['bubble']['radius']
    logphi=math.log(PHI);q=PHI**4;lamG=math.pi**2/(6*math.log(2))
    readings=[('golden full turn',4*logphi),('golden quarter turn',logphi),
              ('golden-angle sampling of planar spiral',GOLDEN_PITCH*GOLDEN_ANGLE),
              ('atlas frequency step',math.log(omega)),('atlas curvature reading',math.log(omega)/2)]
    stability=[]
    for c in pen['retained_frequency_member_controls']:
        measure=persistent_measure(1000,omega,c['delta'])
        row=dict(omega0=c['omega0'],delta=c['delta'],limit=measure,
            percent_miss_dark_energy=100*(measure['undrawn_fraction']/.685-1),
            percent_miss_total_matter=100*(measure['area_fraction']/.315-1),
            contained_prefix=contained_prefix(c['delta'],omega),
            pass_controls=[dict(passes=n,drawn=measure['area_fraction'],undrawn=measure['undrawn_fraction'],
                log_remaining_tail_bound=-(n+1)*math.log(omega)-math.log1p(-1/omega)) for n in (50,100,200,1000)])
        stability.append(row)
    goldarea=persistent_measure(1000,q,GOLDEN_ANGLE)
    for row in (goldarea,):
        row['percent_miss_dark_energy']=100*(row['undrawn_fraction']/.685-1)
        row['percent_miss_total_matter']=100*(row['area_fraction']/.315-1)
    ell=math.sqrt(96/7);m=R**3/(2*ell**2)
    p=kasner(PHI)
    return dict(schema='nsc-atlas-golden-vortex-v1',
        source_stability=stability,
        source_selection='omega0=20 is supplied in the retained parent input; d,j do not derive it. The other two source members are preserved controls.',
        construction=dict(suspension='S(x,t)=(sqrt(1-t^2)*x,t); two independent stored offsets take B1->B2->B3 and boundaries S0->S1->S2',
            inward_map='F_n=Rc*q^(-n)*R_z(n*delta)*S^2(x,t1,t2)',
            vortex='gamma(t)=Rc*q^(-t)*X(2*pi*t,delta*t)',
            exact_similarity='gamma(t+1)=q^(-1)*R_z(delta)*gamma(t)',
            scale_covariance='For round-sphere Dirac operators D_(n+1)=q U D_n U^dagger; normalized pullback is unchanged for every q.',
            scope='A kinematic inward continuation inside Rc. The geometry defines it for every q>1; field equations must select q and the handoff.'),
        golden=dict(phi=PHI,angle=GOLDEN_ANGLE,pitch=GOLDEN_PITCH,
            per_quarter_ratio=PHI,per_full_turn_ratio=q,
            per_golden_angle_sample_ratio=math.exp(GOLDEN_PITCH*GOLDEN_ANGLE),
            local_radial_gap_ratio=1/q,
            angle_miss_from_action_percent=100*(GOLDEN_ANGLE/delta-1),
            full_turn_ratio_miss_from_Omega_percent=100*(q/omega-1),
            phi_cubed_miss_from_Omega_percent=100*(PHI**3/omega-1),
            ln_Omega_over_ln_phi=math.log(omega)/logphi,
            atlas_step_in_golden_full_turns=math.log(omega)/(4*logphi),
            same_full_turn_pitch=False,
            now_positions=[dict(reading=name,log_ratio_per_step=L,position_interval=[t/L for t in targets]) for name,L in readings],
            now_is_position_not_stop=True,physical_observer_event_selected=False),
        golden_area=dict(measure=goldarea,
            definition='Radially project the circle fibers of the same scale-rotation map, retaining the previous geodesic-band footprint rule once per full turn.',
            curve_itself_area=0,
            scope='The smooth vortex, circle-fiber area readout and finite Fibonacci samples are explicit different readouts; neither curve area nor physical stress is inferred from the pen footprint.'),
        fibonacci_samples=[fibonacci_cover(n) for n in (50,100,200,1000)],
        resolution=dict(local_radial='Delta r_n=(1-1/q)*Rc*q^(-n), exactly geometric',
            global_angular='N spherical samples require h>=acos(1-2/N); a spiral alone does not replace the independent angular samples',
            branching='A depth law with N growing proportionally to q^(2n) can provide geometric angular refinement; that population rule is not selected here.'),
        imbalance=dict(counter_rotation_signed_count=0,common_expansion_damping_difference=0,
            real_charge_symmetric_wall_transfer_difference=0,
            identity='Q_comoving proportional a^3*(|A_plus|^2-|A_minus|^2); common damping and conjugate real transmission preserve zero from symmetric data',
            eta_control_if_photon_count_positive=0,eta_target=6.1e-10,eta_control_percent_miss=-100.,
            physical_baryon_photon_eta=None,
            missing='A charge-selective wall/state mechanism and baryon/photon identification are not in the current real scalar/four-form model.'),
        wall=dict(Rc=R,Ui=7/32,Uo=21/32,sigma=.2,
            lapse_inside=1-(7/32)*R*R/3,lapse_outside=1-(21/32)*R*R/3,
            actual_trajectory='R(tau)=Rc*cosh(tau/Rc)',
            kinetic_at_golden_inward_copy=q**-2-1,kinetic_at_Omega_inward_copy=omega**-2-1,
            required_scale_law='Ui,Uo ->q^2*(Ui,Uo), sigma->q*sigma gives Rc->Rc/q for any q; it does not select q.',
            geometry_selected_ratio=None,wall_is_black_hole_horizon=False),
        black_holes=dict(kerr_controls=[kerr_control(a) for a in (0.,.5,.9,.99)],
            entropy='S_BH/kB=A/(4*lP^2); lP remains symbolic',
            Schwarzschild=dict(photon_to_horizon_ratio=1.5,critical_impact_over_mass=3*math.sqrt(3),
                photon_horizon_miss_from_Omega_percent=100*(1.5/omega-1),
                photon_horizon_miss_from_golden_full_percent=100*(1.5/q-1),
                critical_orbit_log_radial_offset_ratio=2*math.pi,
                critical_orbit_radial_offset_ratio=math.exp(2*math.pi),
                radial_factor_miss_from_Omega_percent=100*(math.exp(2*math.pi)/omega-1)),
            scope='Classical photon orbits and near-critical radial offsets, not brightness rings or an established child-universe handoff.'),
        bkl=dict(golden_fixed_point=PHI,era_derivative_at_fixed_point=-PHI**2,
            golden_point_is_attractor=False,kasner_exponents=p.tolist(),
            bounced_exponents=kasner_bounce(p).tolist(),axis_label_cycle_angle=2*math.pi/3,
            axis_cycle_miss_from_golden_angle_percent=100*((2*math.pi/3)/GOLDEN_ANGLE-1),
            volume_law='volume proportional tau; the finite-time metric still has three spatial directions',
            golden_Kretschmann='2/tau^4',
            Gauss_density='1/[ln(2)*(1+x)] on 0<x<1',
            Gauss_mean_log_derivative=lamG,continued_fraction_denominator_factor=math.exp(lamG/2),
            denominator_factor_miss_from_Omega_percent=100*(math.exp(lamG/2)/omega-1),
            golden_denominator_factor=PHI,
            physical_scale_per_era_selected=False,era_count_to_wall_selected=False,
            missing='An initial billiard state, a proper-time/scale map and a compatible curvature cutoff; exponent updates do not create child regions.'),
        core_matching=dict(deSitter_radius=ell,
            continuous_f_Schwarzschild_mass=m,horizon_over_wall_for_that_control=2*m/R,
            boundary_is_inside_black_hole=False,
            invariant_equality_Kasner_time=ell/12**.25,
            scope='Continuous-f and curvature-equality controls only. FMM uses a limiting-curvature spacelike transition, not this timelike Parke turning wall.',
            torsion_comparison='For rho_eff=C_r/a^4-alpha*C_n^2/a^6, flat Lambda0 gives a_min^2=alpha*C_n^2/C_r; spin density and state must be supplied. Microscopic spin is not macroscopic Kerr rotation.',
            mass_inflation='Internal gravitational mass and curvature can grow near an inner horizon while areal radius remains finite; this is not the scale-factor inflation of a child.'),
        no_target_fit=True,no_imported_SI_constants=True,new_spacetime_evolutions=0)
