"""An explicit circle-history construction and its geometric readouts.

The preferred phase increment comes from the existing late-vacuum normal
modes at fixed source coupling. The circle embedding and resolution footprint
are declared kinematic rules; auxiliary phase wheels remain a control.
No target changes an input. Areas and probes use the same persistent drawn
set; geometric residuals are not relabelled baryon charge.
"""
import math
import warnings
import numpy as np
from scipy.integrate import quad, IntegrationWarning


def phase_offset(d=12,j=2):
    """Auxiliary-wheel control, separate from the action-derived phase."""
    ratio=math.sqrt(d/j)
    return 2*math.pi*(ratio-math.floor(ratio))


def action_phase_clock(minimum,omega0):
    """Relative phase of the two linear modes in the existing late vacuum.

    Rescaling either fluctuation by a^(3/2) removes the 3H damping. The
    source's supplied omega0 remains an action input, never fitted here.
    """
    s=minimum['s'];V=minimum['vacuum'];H=math.sqrt(V/3)
    ms2=minimum['canonical_mass_squared'];mc2=omega0**2*math.exp(-12*s)
    ns2=ms2-3*V/4;nc2=mc2-3*V/4
    if min(ns2,nc2)<=0:raise ValueError('two oscillatory linear modes required')
    ratio=math.sqrt(ns2/nc2);delta=2*math.pi*ratio
    return dict(omega0=omega0,s_star=s,H_star=H,nu_volume=math.sqrt(ns2),
        nu_scalar=math.sqrt(nc2),delta=delta,phase_turns=ratio,
        log_expansion_per_scalar_cycle=2*math.pi*H/math.sqrt(nc2),
        formula='delta=2*pi*sqrt((475*x-468)/(32*omega0^2+75*x-84)); 25*x^6-140*x+112=0 at the minimum',
        regime='Linear perturbations about the constant late vacuum; damping removed by a^(3/2). Not the nonlinear formation history.')


def sphere(theta,phi,radius=1.):
    return radius*np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])


def orientations(n,delta):
    if type(n) is not int or not 1<=n<=10000:raise ValueError('finite 1..10000 passes required')
    return np.remainder(np.arange(n,dtype=float)*delta,math.pi)


def angular_gaps(n,delta):
    a=np.sort(orientations(n,delta))
    return np.diff(np.r_[a,a[0]+math.pi])


def current_resolution_area(n,delta,epsilon):
    """Exact band-union formula evaluated without subtracting nearly equal angles."""
    if not 0<epsilon<math.pi/2:raise ValueError('positive sub-hemisphere halfwidth required')
    k=math.sin(epsilon);values=[]
    for gap in angular_gaps(n,delta):
        a=gap/2
        if a<=epsilon:values.append(gap);continue
        sa,ca=math.sin(a),math.cos(a)
        root=math.sqrt(max(0,sa*sa-k*k))
        # Stable evaluation of a-acos(cos(a)/sqrt(1-k^2)).
        difference=math.atan2(ca*k*k/(sa+root),ca*ca+sa*root)
        t=math.acos(min(1.,max(0.,k*ca/(sa*math.sqrt(1-k*k)))))
        values.append(2*difference+2*k*t)
    return math.fsum(values)/math.pi


def longitude_coverage(z,centres,widths):
    """Covered fraction of the longitude circle at |cos(theta)|=z."""
    st=math.sqrt(max(0.,1-z*z))
    if st<=math.sin(max(widths)):return 1.
    intervals=[]
    for centre,width in zip(centres,widths):
        a=math.asin(min(1.,math.sin(width)/st))
        lo,hi=centre-a,centre+a
        if lo<0:intervals.extend(((0.,hi),(lo+math.pi,math.pi)))
        elif hi>math.pi:intervals.extend(((lo,math.pi),(0.,hi-math.pi)))
        else:intervals.append((lo,hi))
    intervals.sort();end=0.;parts=[]
    for lo,hi in intervals:
        if hi>end:parts.append(hi-max(end,lo));end=hi
    return math.fsum(parts)/math.pi


def persistent_measure(n,omega,delta,power_divisor=1,mesh=4096):
    """Area fraction and z-probe stiffness from persistent geodesic pen bands.

    Width of pass k is Omega^(-(k+1)/p) radians. Omitted tail area is bounded
    by its width sum. Monotone Darboux sums give a finite mesh enclosure;
    a separate adaptive quadrature supplies an interior point estimate.
    """
    if power_divisor not in (1,2) or omega<=1 or mesh<64:raise ValueError('declared finite positive refinement required')
    r=omega**(-1/power_divisor);keep=1
    while keep<n and r**(keep+1)/(1-r)>1e-11:keep+=1
    tail=0. if keep==n else r**(keep+1)/(1-r)
    c=orientations(keep,delta);w=[r**(k+1) for k in range(keep)]
    f=lambda z:longitude_coverage(z,c,w)
    z=np.linspace(0,1,mesh+1);vals=np.array([f(x) for x in z])
    # Conservative allowance for finite floating arithmetic, separate from
    # the analytic quadrature-width and omitted-area bounds.
    rounding=1e-8
    area_lo=math.fsum(vals[:-1])/mesh-rounding
    area_hi=math.fsum(vals[1:])/mesh+tail+rounding
    weighted=z*z*vals
    moment_lo=math.fsum(weighted[:-1])/mesh-rounding
    moment_hi=math.fsum(weighted[1:])/mesh+tail+rounding
    zcap=math.cos(w[0]);messages=[]
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always',IntegrationWarning)
        area,ae=quad(f,0,zcap,epsabs=3e-7,epsrel=3e-7,limit=600)
        moment,me=quad(lambda z:z*z*f(z),0,zcap,epsabs=3e-7,epsrel=3e-7,limit=600)
        messages=[str(x.message) for x in caught]
    area+=1-zcap;moment+=(1-zcap**3)/3
    if not area_lo<=area<=area_hi or not moment_lo<=moment<=moment_hi:
        raise ArithmeticError('independent quadrature outside monotone bracket')
    return dict(passes=n,power_divisor=power_divisor,retained_bands=keep,
        area_fraction=area,area_enclosure=[max(0.,area_lo),min(1.,area_hi)],
        undrawn_fraction=1-area,undrawn_enclosure=[max(0.,1-area_hi),min(1.,1-area_lo)],
        z_probe_moment=moment,z_probe_enclosure=[max(0.,moment_lo),min(1/3,moment_hi)],
        omitted_area_bound=tail,rounding_allowance=rounding,quadrature_mesh=mesh,
        quadrature_error_indicators=[ae,me],quadrature_warnings=messages)


def phase_residual(n,delta):
    return abs(math.sin(n*delta/2))/n


def quadratic_irrational_bound(n):
    """Lower bound on |sin(pi*n*sqrt(6))|/n from a nonzero integer norm."""
    floor=math.isqrt(6*n*n)
    nearest=floor+int(24*n*n>(2*floor+1)**2)
    return 2/(n*(nearest+n*math.sqrt(6)))


def calculation(omega,radius,eta_target=6.1e-10,log_targets=(136.84,138.26),*,clock=None):
    auxiliary=clock is None
    delta=phase_offset() if auxiliary else clock['delta']
    logomega=math.log(omega);cases=[]
    for n,p in ((100,1),(200,2),(200,1)):
        ell=n*logomega/p;eps=math.exp(-ell)
        h=float(max(angular_gaps(n,delta))/2)
        residual=phase_residual(n,delta)
        f=current_resolution_area(n,delta,eps)
        cases.append(dict(passes=n,power_divisor=p,prescribed_log_resolution=ell,
            achieved_fill_distance=h,achieved_log_uniform_resolution=-math.log(h),
            achieved_log_miss_percent=[100*(-math.log(h)/t-1) for t in log_targets],
            required_log_circle_count_lower_bound=ell+math.log(math.pi/2),
            dipole_residual=residual,dipole_miss_percent_if_identified_with_eta=100*(residual/eta_target-1),
            residual_lower_bound=quadratic_irrational_bound(n) if auxiliary else None,net_scalar_signed_count=0,
            late_vacuum_log_expansion=(n*clock['log_expansion_per_scalar_cycle'] if clock else None),
            current_width_area=f,log_current_width_area=math.log(f),
            current_width_undrawn_expression=f'1 - {f:.12e}'))
    first=persistent_measure(1,omega,delta)
    second=persistent_measure(2,omega,delta)
    later=persistent_measure(100,omega,delta)
    curved=persistent_measure(200,omega,delta,2)
    f1=math.sin(1/omega);g1=f1/2-f1**3/6
    ratio=second['z_probe_moment']/g1
    tail100=math.exp(-101*logomega)/(1-1/omega)
    for row in (later,curved):
        f=row['area_fraction']
        row['conditional_binary_to_three_assignment']=dict(Omega_Lambda=1-f,Omega_c=0.,Omega_b=f,
            percent_miss=dict(Omega_Lambda=100*((1-f)/.685-1),Omega_c=-100.,Omega_b=100*(f/.049-1)),
            rule='Unresolved area assigned to vacuum; resolved area to baryons; no distinct clustering category. This is a control assignment, not a stress derivation.')
        row['general_three_way_family']='((1-f)*(1-chi), (1-f)*chi, f), with unselected outside clustering fraction chi'
    return dict(schema='nsc-atlas-pen-geometry-v1',
        construction=dict(d=12,j=2,omega=omega,Parke_radius=radius,
            sphere='X(u,v)=R*(sin(u)*cos(v),sin(u)*sin(v),cos(u))',
            phase_rule=('Two auxiliary phase circles have radii sqrt(d),sqrt(j) and equal arc speeds.' if auxiliary else 'Advance by the volume-mode phase during one scalar-mode period of the fixed late-vacuum action.'),
            phase_turns_exact=('sqrt(6)-2' if auxiliary else 'sqrt((475*x-468)/(32*omega0^2+75*x-84))'),
            delta_exact=('2*pi*(sqrt(6)-2)' if auxiliary else '2*pi*sqrt((475*x-468)/(32*omega0^2+75*x-84))'),delta=delta,
            offset_derived_from_fixed_action=not auxiliary,
            phase_to_circle_embedding_is_declared=True,
            clock=clock,
            drawing='Pass k draws the complete great circle at v=k*delta mod pi; independent u draws along that circle.',
            refinement='Pass k has angular resolution footprint Omega^(-(k+1)/p), p=1 or2; initial angular unit is one radian.',
            persistence='Previously drawn bands remain present. Re-resolving all curves at the latest width is a separate coverage diagnostic.'),
        exact_results=dict(dense_sphere_closure=True,finite_ideal_line_area=0,
            persistent_ink_is_open_dense=True,
            persistent_undrawn_area_lower_bound_p1=1-1/(omega-1),
            persistent_complement='Closed, empty interior, but positive area for p=1: total band area <=sum Omega^-k=1/(Omega-1)<1.',
            rank='|dX/du cross dX/dv|=R^2*abs(sin(u)); rank2 away from coordinate poles',
            exact_finite_return=False,
            fill_distance='h_n=max angular gap /2 >=pi/(2n)',
            exponentially_refining_one_circle_per_pass_covers_sphere=False,
            coverage_obstruction='n*Omega^(-n/p)<=p/(e*ln(Omega))<pi/2, so epsilon_n<h_n at every positive n',
            sharp_constructive_refinement='M=ceil(pi/(2*epsilon)) equally spaced meridians, rotated by n*delta, achieve h=pi/(2M)<=epsilon. This requires branching many circles per level and does not select the final level.',
            paired_residual='sum[z_k-z_(k+1)]/(2n)=(1-exp(i*n*delta))/(2n)',
            charge='Each pass has one positive and one negative scalar charge: total zero for every delta.',
            asymmetry_bound=('|sin(pi*n*sqrt6)|/n >= 2/[n*(nearest(n*sqrt6)+n*sqrt6)]' if auxiliary else 'Irrationality follows from the Eisenstein-7 minimum polynomial and its nonconstant rational frequency-square ratio.')),
        target_comparisons=cases,
        persistent_coverage=[later,curved],
        probe=dict(action='E_n(u)=gamma/2 * integral_drawn (u dot X/R)^2 dA/(4*pi*R^2)',
            coupling='g_z(n)=gamma*integral_0^1 z^2*covered_longitude_fraction(z) dz',
            total_inside_plus_outside_z_coupling='gamma/3 at every depth',
            depth_1=first,depth_2=second,ratio_depth_2_over_1=ratio,
            ratio_enclosure=[x/g1 for x in second['z_probe_enclosure']],
            depth_100_over_1=later['z_probe_moment']/g1,
            depth_200_over_100_expression='1 <= ratio <= 1 + bound',
            depth_200_over_100_increment_bound=tail100/later['z_probe_enclosure'][0],
            fundamental_force_identification_derived=False),
        missing_links=[
            ('Offset: derive the auxiliary phase metric/equal-speed state from the parent action.' if auxiliary else 'History: transport the fixed-action linear clocks through formation and nesting; the source frequency is supplied and the phase-to-drawing map is declared.'),
            'History: derive a stopping event or branching circles-per-pass rule; irrational rotation alone selects no finite depth.',
            'Asymmetry: identify a baryon/photon observable with a sourced charge-producing process; a dipole residual is not net charge.',
            'Fractions: derive footprint/energy weights and the outside clustering partition; binary area alone supplies two measures.',
            'Forces: derive the physical probe coupling to this geometric action before identifying fundamental interactions.'],
        no_imported_physical_constants=True,no_target_fit=True,new_spacetime_evolutions=0)
