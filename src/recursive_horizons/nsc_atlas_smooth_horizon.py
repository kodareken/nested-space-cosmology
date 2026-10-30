"""Smooth horizon continuation, charge balance and a sourced drain comparison.

The Parke shell limit, stationary vacuum solution and Vaidya null-source
comparison have separate assumptions. No observed target sets a parameter.
"""
import math
from scipy.optimize import brentq
from .nsc_atlas_parke_nucleation import thin_wall
from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_pen_geometry import action_phase_clock


def metric_f(radius, mass, vacuum):
    if radius <= 0:
        raise ValueError('positive areal radius required')
    return 1-2*mass/radius-vacuum*radius*radius/3


def horizons(mass, vacuum):
    """Simple positive Kottler roots, in the declared Einstein-frame units."""
    if vacuum <= 0 or mass < 0 or 9*mass*mass*vacuum >= 1:
        raise ValueError('nonnegative mass, positive vacuum, sub-Nariai branch required')
    if mass == 0:
        return dict(black_hole=None, cosmological=math.sqrt(3/vacuum))
    f=lambda r:metric_f(r,mass,vacuum)
    middle=1/math.sqrt(vacuum)
    return dict(black_hole=brentq(f,mass,middle,xtol=1e-14),
                cosmological=brentq(f,middle,math.sqrt(3/vacuum),xtol=1e-14))


def infall_at(radius, mass, vacuum, energy=1.):
    """Regular future ingoing radial geodesic in advanced EF coordinates."""
    if energy <= 0:
        raise ValueError('positive Killing energy required')
    f=metric_f(radius,mass,vacuum)
    if energy*energy-f < 0:
        raise ValueError('outside this geodesic segment')
    radial=-math.sqrt(energy*energy-f)
    advanced=1/(energy-radial)
    return dict(dr_dtau=radial,dv_dtau=advanced,
                norm=-f*advanced*advanced+2*advanced*radial)


def horizon_response(radius, vacuum):
    """dR/dm for a simple marginal root; no drain rate is supplied."""
    den=1-vacuum*radius*radius
    if radius <= 0 or abs(den)<1e-13:
        raise ValueError('positive simple root required')
    return 2/den


def calculation():
    ui,uo=7/32,21/32
    bubble=thin_wall();rc=bubble['radius']
    minimum=VolumeParent().vacuum_roots()[1]
    lam=minimum['vacuum'];ell=math.sqrt(3/lam)
    clocks=action_phase_clock(minimum,20.)
    phi=(1+math.sqrt(5))/2
    phase_controls=[]
    for label,delta in (('retained_action_phase',clocks['delta']),('golden_angle',2*math.pi/phi**2)):
        polarization=math.sin(delta)
        phase_controls.append(dict(phase_source=label,delta=delta,charge_per_total_quanta=polarization,
            eta_benchmark_miss_percent=100*(polarization/6.1e-10-1),
            scope='Equal-amplitude degenerate complex-oscillator initial state; this ratio is not baryon/photon eta.'))
    # This inversion is a diagnostic, explicitly NOT a selected mass.
    required_m=rc*(1-lam*rc*rc/3)/2
    roots=horizons(required_m,lam)
    root_rows=[]
    for label,r in roots.items():
        root_rows.append(dict(kind=label,radius=r,f_residual=metric_f(r,required_m,lam),
            dR_dm=horizon_response(r,lam),
            kretschmann=48*required_m**2/r**6+8*lam**2/3,
            crossing=infall_at(r,required_m,lam)))
    tensions=[.2,.02,.002,.0002]
    shell_rows=[dict(sigma=s,radius=thin_wall(sigma=s)['radius'],
        asymptotic_radius=3*s/(uo-ui)) for s in tensions]
    return dict(schema='nsc-atlas-smooth-horizon-v1',
        inputs=dict(Ui=ui,Uo=uo,sigma_reference=.2,kappa=1.,d=12,ancestors=2,
                    eta_comparison=6.1e-10,Omega_Lambda_comparison=.685),
        shell_limit=dict(reference_radius=rc,rows=shell_rows,limit_radius=0.,
            asymptotic='Rc=3sigma/DeltaU + O(sigma^3), fixed DeltaU>0',
            finite_horizon_survives=False,
            Ricci_scalar_jump=4*(uo-ui),
            reference_f_inside=metric_f(rc,0,ui),reference_f_outside=metric_f(rc,0,uo),
            vacuum_horizon_inside=math.sqrt(3/ui),vacuum_horizon_outside=math.sqrt(3/uo),
            scope='Constant-vacuum Parke control; endpoints are the retained nonstationary snapshots, not new vacua of the volume action.'),
        smooth_match=dict(metric='ds^2=-f dv^2+2dv dr+r^2 dOmega^2',
            f='1-2m/r-Lambda r^2/3',determinant='-r^4 sin(theta)^2',
            continuity='[m]=-[Lambda] R^3/6; then [f_prime]=-[Lambda] R',
            conclusion='Smooth common-normalization Kottler matching at finite R forces [Lambda]=[m]=0.',
            domain='Spherical vacuum pieces in one common advanced EF normalization. Null matching uses transverse curvature; no claim about all matter-supported junctions.',
            flux='d(Z3 *F4 + mu chi)=e delta_membrane; without membrane current, q is continuous.',
            missing_source='Keeping nonzero e while setting sigma=0 retains a charged interface; it does not make a source-free horizon.'),
        admitted_vacuum=dict(s_star=minimum['s'],Lambda=lam,horizon_m0=ell,
            horizon_over_old_Rc=ell/rc,Lambda_Rc_squared=lam*rc*rc,
            a_double_dot_over_a=lam/3,H=math.sqrt(lam/3),w=-1.,
            flat_vacuum_only_Omega_Lambda=1.,
            control_miss_percent=100*(1/.685-1),
            scope='Stationary fields in the existing effective Einstein-frame action. m=0 de Sitter is regular; m>0 Kottler is a family, not a derived birth or selected black-hole mass.',
            SI_clock_selected=False),
        radius_inversion_control=dict(required_mass=required_m,mass_selected=False,
            reason='Mass required to force a black-hole horizon to equal old Rc; diagnostic only.',
            roots=root_rows),
        horizon_zero=dict(geodesic_limit='dr/dtau=-E; dv/dtau=1/(2E); norm=-1 at f=0',
            interpretation='The static chart ends, the smooth four-dimensional spacetime continues. No new independent dimension follows from f=0.',
            singularity_distinction='For m>0, K=48m^2/r^6+8Lambda^2/3 diverges at r=0. Crossing the regular horizon has not supplied a central bounce.'),
        charge=dict(geometric_mirror='CPT pairs Q_minus(-t)=-Q_plus(t), leaving the local magnitude free.',
            conservation='Q(V2)-Q(V1)=-integral_boundary J.dSigma + integral_bulk S_B',
            symmetry='C rho C^-1=rho, [C,H(t)]=0, C Q_V C^-1=-Q_V imply <Q_V(t)>=0.',
            symmetry_domain='Internal charge conjugation preserves the chosen region; no charge-selective couplings or initial bias. Global CPT alone is weaker and need not force each side neutral.',
            symmetric_control_eta=0.,control_miss_percent=-100.,physical_eta=None,
            missing='A physically identified charge current, its inherited state/transport law and photon normalization.',
            zero_global_allows_nonzero_regional=True,
            r_time_swap_is_charge_conjugation=False,
            extra_dimension_current_derived=False),
        offset_bridge=dict(
            action='L=(xdot^2+ydot^2-omega^2*(x^2+y^2))/2; psi=(x+i*y)/sqrt(2)',
            state='x=A cos(omega t), y=B cos(omega t-delta)',
            current='Q=x*ydot-y*xdot=A*B*omega*sin(delta)',
            energy='E=omega^2*(A^2+B^2)/2',
            polarization='Q/(E/omega)=2*A*B*sin(delta)/(A^2+B^2)',
            controls=phase_controls,
            expanding_state='x,y acquire a^-3/2 and omega is replaced in phase by nu=sqrt(omega^2-9H^2/4); a^3 Q=A*B*nu*sin(delta) is constant.',
            nondegenerate_obstruction='Qdot=(omega_x^2-omega_y^2)*x*y; the existing volume and scalar modes are not a degenerate U(1) pair.',
            existing_nu_volume=clocks['nu_volume'],existing_nu_scalar=clocks['nu_scalar'],
            scope='Constructed charged comparison, not an added sector silently attributed to the fixed neutral action. The offset specifies initial charge; a horizon has not selected it.'),
        baseline_and_oscillation=dict(
            exact_stress='rho=Vstar+(chidot^2+omega^2 chi^2)/2; p=-Vstar+(chidot^2-omega^2 chi^2)/2',
            phase_ensemble='chi=A a^(-3/2) cos(nu t+theta), uniform theta: <rho_osc>=A^2 omega^2/(2 a^3), <p_osc>=0.',
            readouts='w_vacuum=-1; w_oscillation=0 for the stated constant-H phase ensemble; ordinary fast-time averaging is approximate.',
            scope='Quadratic test oscillator at fixed volume on de Sitter; no self-consistent finite matter cosmology or baryon/dark-matter partition is inferred.',
            amplitude_selected=False),
        historical_material_wall=dict(R='Rc cosh(tau/Rc)',
            R_double_dot_over_R=1/rc**2,R_double_dot_initial=1/rc,
            H_wall='tanh(tau/Rc)/Rc',
            acceleration_over_inner_vacuum=(1/rc**2)/(ui/3),
            scope='Timelike shell proper time; not a null horizon clock. Rdot is an areal derivative, not local three-speed.'),
        drain_comparison=dict(metric='ds^2=-f dv^2+2dv dr+r^2 dOmega^2; m=m(v)',
            null_stress='T_vv=2m_prime/(kappa r^2)',trace=0.,
            root_motion='R_prime=2m_prime/(1-Lambda R^2)',
            signs='m_prime<0 shrinks the black-hole marginal root and enlarges the cosmological root; this advanced patch then needs negative ingoing null energy.',
            flux_rate_selected=False,
            vacuum_added_by_null_flux=0.,
            no_unique_time_acceleration='R_double_prime=2m_double_prime/D+8Lambda R m_prime^2/D^3, D=1-Lambda R^2.',
            scope='Sourced Vaidya-de Sitter comparison. Null fluid is an added effective source, not yet a solution of the retained scalar sector. Marginal/trapping roots are not automatically global event horizons.',
            observational_fit=None),
        dispositions=[
            'Fixed-contrast zero-tension Parke limit: radius tends to zero; retain as one-line completed route.',
            'Smooth same-vacuum horizon: explicit regular continuation and finite geodesic crossing established.',
            'CPT identification: relational opposite sides can be paired, but charge magnitude and transport remain state data.',
            'Moving neutral geometry: symmetric expected charge stays zero even when the horizon moves.',
            'Drain comparison: opposite horizon motion follows from one mass function, while its traceless flux does not generate vacuum Lambda.'],
        physical_eta_or_Omega_prediction=False)
