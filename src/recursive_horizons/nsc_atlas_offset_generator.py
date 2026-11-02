"""One spatial interval and one accumulated offset: explicit metric completions.

Retained history, the distance kernel and metric completion are declared
construction rules. This is not a physical dimension-selection or time-origin
theorem, nor an assertion that a finite drawing fills a volume.
"""
import math
import mpmath as mp
import numpy as np
from scipy.spatial import SphericalVoronoi,cKDTree


def opposite_winding_path(position, contraction, turns=1, radius=1.):
    """A closed forward/reverse loop and its null homotopy, not a time law.

    position is the retained drawing coordinate in [0,2]. At contraction=0
    the first half has degree turns and the second degree -turns. During
    contraction only the combined loop stays closed. The two halves need
    not remain individually closed: their contact is allowed to reconnect.
    """
    if type(turns) is not int or turns == 0:
        raise ValueError('nonzero integer winding required')
    if not math.isfinite(radius) or radius <= 0:
        raise ValueError('positive finite radius required')
    if not math.isfinite(contraction) or not 0 <= contraction <= 1:
        raise ValueError('contraction must lie in [0,1]')
    x = np.asarray(position, dtype=float)
    if not np.all(np.isfinite(x)) or np.any((x < 0) | (x > 2)):
        raise ValueError('drawing coordinate must lie in [0,2]')
    lift = np.minimum(x, 2 - x)
    return radius * np.exp(2j * np.pi * turns * (1 - contraction) * lift)


def winding_variation(contraction, turns=1, radius=1.):
    """Half the integral of |dz/dx|^2 over [0,2], in drawing units.

    This geometric functional has no supplied conversion to physical energy.
    It is the minimum for a fixed circle radius, endpoints and midpoint lift.
    """
    opposite_winding_path(0., contraction, turns, radius)
    return 4 * math.pi**2 * radius**2 * turns**2 * (1 - contraction)**2


def balanced_output_directions(count):
    """Competing balanced channel geometries; count is supplied, not selected."""
    if type(count) is not int or count < 2:
        raise ValueError('at least two output channels required')
    angles = 2 * np.pi * np.arange(count) / count
    return np.c_[np.cos(angles), np.sin(angles), np.zeros(count)]


def winding_transport(position, offset, turns=1, radius=1., closed=False):
    """Exact conservative phase motion, conditional on geometric inertia.

    The new rule is phase_tt=phase_xx in normalized offset units, with
    initial phase equal to the joined winding lift and zero initial velocity.
    closed=True identifies x modulo 2; otherwise the lift is zero outside
    [0,2]. The circle target remains the same in either case. Derivatives
    at the finitely many moving corners are assigned representatives only.
    """
    opposite_winding_path(0., 0., turns, radius)
    if not math.isfinite(offset) or offset < 0:
        raise ValueError('finite nonnegative offset required')
    if type(closed) is not bool:
        raise ValueError('closed must explicitly be a boolean')
    x = np.asarray(position, dtype=float)
    if not np.all(np.isfinite(x)):
        raise ValueError('finite drawing coordinates required')
    t = offset % 2 if closed else offset

    def lift_and_slope(y):
        y = np.remainder(y, 2) if closed else y
        lift = 2*np.pi*turns*np.maximum(0, np.minimum(y, 2-y))
        slope = 2*np.pi*turns*np.where((y > 0) & (y < 1), 1.,
                                     np.where((y > 1) & (y < 2), -1., 0.))
        return lift, slope

    right, right_slope = lift_and_slope(x-t)
    left, left_slope = lift_and_slope(x+t)
    phase = (right+left)/2
    velocity = (left_slope-right_slope)/2
    spatial = (left_slope+right_slope)/2
    kinetic = radius**2*velocity**2/2
    variation = radius**2*spatial**2/2
    return dict(phase=phase, circle=radius*np.exp(1j*phase),
                offset_derivative=velocity, drawing_derivative=spatial,
                kinetic=kinetic, variation=variation, density=kinetic+variation,
                flux=-radius**2*velocity*spatial,
                right_density=radius**2*right_slope**2/4,
                left_density=radius**2*left_slope**2/4)


def winding_transport_ledger(offset, turns=1, radius=1., closed=False):
    """Integrate the piecewise-constant densities at their exact breakpoints."""
    winding_transport(0., offset, turns, radius, closed)
    t = offset % 2 if closed else offset
    if closed:
        cuts = [0., 2.]
        for k in range(-2, 5):
            cuts.extend(x for x in (k-t, k+t) if 0 < x < 2)
    else:
        cuts = [k+sign*t for sign in (-1, 1) for k in (0., 1., 2.)]
    cuts = np.unique(cuts)
    widths = np.diff(cuts)
    values = winding_transport((cuts[:-1]+cuts[1:])/2, offset, turns, radius, closed)
    result = {name:float(np.dot(widths, values[name])) for name in
              ('kinetic','variation','density','right_density','left_density')}
    result['initial_budget'] = winding_variation(0, turns, radius)
    result['conservation_residual'] = result['density']-result['initial_budget']
    if not closed:
        delivered_each = math.pi**2*radius**2*turns**2*min(offset,2.)
        result['received_left_of_initial_interval'] = delivered_each
        result['received_right_of_initial_interval'] = delivered_each
        result['retained_in_initial_interval'] = result['initial_budget']-2*delivered_each
    return result


def winding_transport_calculation():
    """Small reference ledger; channel speed and inertia are explicit premises."""
    cases = []
    for closed in (False, True):
        cases.append(dict(boundary='periodic' if closed else 'open line',
                          history=[dict(offset=t, **winding_transport_ledger(t,closed=closed))
                                   for t in (0., .25, .5, 1., 2.)]))
    return dict(schema='nsc-offset-winding-transport-v1',
        added_premises=['L2 circle-phase inertia with the same metric as squared variation',
                        'zero initial offset velocity',
                        'normalized offset propagation speed one',
                        'open or periodic drawing coordinate, compared separately'],
        equation='phase_tt = phase_xx',
        solution='phase(x,t) = (g(x-t)+g(x+t))/2',
        local_ledger='partial_t[(R^2/2)(phase_t^2+phase_x^2)] + partial_x[-R^2 phase_t phase_x] = 0',
        two_channels='right and left characteristic densities, each half the initial budget',
        cases=cases,
        open_packets_disjoint_after_offset=1.,
        periodic_full_recurrence_offset=2.,
        missing_connection='selection of the inertial rule, initial velocity and receiving/open geometry from the offset generator',
        selected_by_offset_recurrence=False,
        physical_energy_unit_selected=False,
        horizon_or_child_derived=False,
        no_observational_input=True)


def lifted_offset_flow(theta, phi, alpha, offset):
    """An exact affine interpolation of the discrete map on retained lifts.

    No modulo operation is applied: intermediate offsets require the actual
    retained history, since different torus lifts need not agree between steps.
    This interpolation is not uniquely selected by the integer-time map.
    """
    if not math.isfinite(alpha) or not math.isfinite(offset):
        raise ValueError('finite offset and increment required')
    theta, phi = np.broadcast_arrays(np.asarray(theta,dtype=float),np.asarray(phi,dtype=float))
    if not np.all(np.isfinite(theta)) or not np.all(np.isfinite(phi)):
        raise ValueError('finite retained phase lifts required')
    return (theta+alpha*offset,
            phi+offset*theta+alpha*offset*(offset-1)/2)


def offset_driven_winding(position, offset, alpha, coupling_speed, turns=1, radius=1., closed=False):
    """A local coupled extension with the original map as its uniform sector.

    In phase turns: phi_t=theta-alpha/2, theta_t=alpha+c^2*phi_xx.
    c is a declared coupling, not derived from alpha. The geometric circle
    phase is 2*pi*phi; the existing transport solution evolves its relative
    shape. The full background-forced kinetic budget is not called conserved.
    """
    if not math.isfinite(coupling_speed) or coupling_speed < 0:
        raise ValueError('finite nonnegative coupling speed required')
    _, background = lifted_offset_flow(0.,0.,alpha,offset)
    relative = winding_transport(position,coupling_speed*offset,turns,radius,closed)
    phi = relative['phase']/(2*np.pi)+background
    theta = alpha*offset+coupling_speed*relative['offset_derivative']/(2*np.pi)
    return dict(theta_lift=theta,phi_lift=phi,
                circle=radius*np.exp(2j*np.pi*phi),
                relative_phase=relative['phase'],
                relative_offset_derivative=coupling_speed*relative['offset_derivative'],
                drawing_derivative=relative['drawing_derivative'])


def offset_connection_calculation():
    """Compare couplings sharing exactly the same uniform offset history."""
    alpha = float(offset_value('golden'))
    rows = []
    for speed in (0., .5, 1., 2.):
        ledger = winding_transport_ledger(speed)
        rows.append(dict(coupling_speed=speed,update_count=1,
                         fraction_received_outside_initial_interval=
                         (ledger['received_left_of_initial_interval']+
                          ledger['received_right_of_initial_interval'])/ledger['initial_budget'],
                         retained_fraction=ledger['retained_in_initial_interval']/ledger['initial_budget'],
                         shape_variation=ledger['variation'],
                         uniform_theta=alpha,uniform_phi=0.))
    return dict(schema='nsc-offset-contact-connection-v1',
        status='retained isolated-offset diagnostic; active source is enclosing deformation',
        offset_increment=alpha,
        increment_status='retained golden comparison; not newly selected; identities hold for every real alpha',
        exact_integer_history='theta_n=theta_0+n alpha; phi_n=phi_0+n theta_0+alpha n(n-1)/2',
        lifted_flow='theta_t=alpha; phi_t=theta-alpha/2',
        original_pointwise_stationary_winding=dict(
            initial_variation=winding_variation(0.),variation_after_100_updates=winding_variation(0.),
            transported_fraction=0.),
        coupled_extension='theta_t=alpha+c^2 phi_xx; phi_t=theta-alpha/2',
        relative_equation='(phi-background)_tt=c^2 (phi-background)_xx',
        same_history_different_transfer=rows,
        periodic_received_outside_domain_fraction=0.,
        coupling_selected_by_recurrence=False,
        receiving_domain_selected_by_recurrence=False,
        missing_connection='geometric rule fixing neighbouring-point coupling and the receiving domain',
        no_observational_input=True)


def inherited_compression_offset(parent_deformation, child_deformation):
    """Relative logarithmic length compression of corresponding material patches.

    Input histories share reference labels and a comparison frame. Positive
    deformation determinants give their local volume factors. The returned
    measure is derived from those factors; its identification with the
    proposed time/offset is a premise, not an angular calibration or dynamics.
    No isolated alpha, sampling duration or physical stress is supplied here.
    """
    parent = np.asarray(parent_deformation,dtype=float)
    child = np.asarray(child_deformation,dtype=float)
    if (parent.shape != child.shape or parent.ndim != 3 or parent.shape[0] < 1
            or parent.shape[1] != parent.shape[2] or parent.shape[1] < 1):
        raise ValueError('matching nonempty histories of square deformation maps required')
    if not np.all(np.isfinite(parent)) or not np.all(np.isfinite(child)):
        raise ValueError('finite deformation histories required')
    ps,plog = np.linalg.slogdet(parent)
    cs,clog = np.linalg.slogdet(child)
    if np.any(ps <= 0) or np.any(cs <= 0):
        raise ValueError('orientation-preserving nonsingular deformation maps required')
    dimension = parent.shape[1]
    contrast = (plog-plog[0]-clog+clog[0])/dimension
    return dict(compression=contrast,increments=np.diff(contrast),
                parent_log_volume_change=plog-plog[0],
                child_log_volume_change=clog-clog[0],dimension=dimension)


def inherited_offset_calculation():
    """Analytic deformation controls, not a supplied formation history."""
    samples = np.array([0.,.25,.5,1.])
    maps = [np.exp(-rate*samples)[:,None,None]*np.eye(3) for rate in (1.,2.,4.)]
    first = inherited_compression_offset(maps[0],maps[1])
    second = inherited_compression_offset(maps[1],maps[2])
    direct = inherited_compression_offset(maps[0],maps[2])
    common = inherited_compression_offset(maps[1],maps[1])
    return dict(schema='nsc-parent-compression-offset-v1',
        formula='chi_pc = (log(Jp/Jp0)-log(Jc/Jc0))/dimension',
        interpretation='candidate inherited compression coordinate; angular/time identification remains a premise',
        composition='chi_grandparent_child = chi_grandparent_parent + chi_parent_child',
        increment='delta_chi is computed from deformation history, not independently selected',
        normalization='logarithmic length compression; no SI time or angle conversion',
        analytic_controls=dict(sample_parameter=samples.tolist(),
            deformation_rates=[1.,2.,4.],
            grandparent_parent=first['compression'].tolist(),
            parent_child=second['compression'].tolist(),
            grandparent_child=direct['compression'].tolist(),
            composition_residual=(direct['compression']-first['compression']-second['compression']).tolist(),
            grandparent_child_increments=direct['increments'].tolist(),
            common_deformation_relative_offset=common['compression'].tolist()),
        actual_parent_formation_history_computed=False,
        angular_phase_conversion_derived=False,
        pressure_to_deformation_dynamics_derived=False,
        no_observational_input=True)


def intrinsic_phase_cost(deformation, gradient_factor):
    """Intrinsic phase Dirichlet cost under an affine change of actual distances.

    Rows of gradient_factor factor the reference phase-gradient moment:
    M=C.T@C=integral grad(psi) grad(psi).T dV_reference. The factor is
    geometric input from the existing phase, not additional field species.
    Positive determinant deformation A gives Q=det(A)*||C A^-1||^2/2.
    It is a candidate geometric cost, not an assigned physical energy.
    """
    A=np.asarray(deformation,dtype=float);C=np.asarray(gradient_factor,dtype=float)
    if (A.ndim!=2 or A.shape[0]!=A.shape[1] or A.shape[0]<1
            or C.ndim!=2 or C.shape[0]<1 or C.shape[1]!=A.shape[0]):
        raise ValueError('square deformation and matching gradient factor required')
    if not np.all(np.isfinite(A)) or not np.all(np.isfinite(C)):
        raise ValueError('finite geometric inputs required')
    sign,log_volume=np.linalg.slogdet(A)
    if sign<=0:raise ValueError('orientation-preserving nonsingular deformation required')
    transformed=np.linalg.solve(A.T,C.T).T
    if not np.all(np.isfinite(transformed)):
        raise ValueError('gradient exceeds this matrix representation; continue in log-scale principal coordinates')
    largest=float(np.max(np.abs(transformed)))
    if largest==0:
        return dict(cost=0.,log_cost=None,log_volume=float(log_volume),
                    directional_fraction=[0.]*A.shape[0],
                    logarithmic_force_per_cost=None)
    scaled=transformed/largest
    components=np.sum(scaled**2,axis=0);total=float(components.sum())
    fractions=components/total
    log_cost=float(log_volume+2*math.log(largest)+math.log(total/2))
    cost=math.exp(log_cost) if -744<log_cost<709 else None
    return dict(cost=cost,log_cost=log_cost,log_volume=float(log_volume),
                directional_fraction=fractions.tolist(),
                logarithmic_force_per_cost=(2*fractions-1).tolist())


def conservative_neck_shape(offset):
    """Exact constrained shape control; no formed child or parent history.

    Fix volume and A=diag(exp(-s/2),exp(-s/2),exp(s)). For reference
    gradient factor [0,0,1], Q=exp(-2s)/2. With logarithmic-strain metric
    sum_i(d log lambda_i)^2=(3/2)ds^2, unit inertial motion from rest obeys
    s_tt=(2/3)exp(-2s). This is an explicit candidate kinetic rule.
    """
    if not math.isfinite(offset) or offset<0:raise ValueError('finite nonnegative offset required')
    a=math.sqrt(2/3);u=a*offset
    s=float(np.logaddexp(u,-u)-math.log(2))
    velocity=a*math.tanh(u)
    potential=.5*math.exp(-2*s)
    kinetic=.75*velocity**2
    return dict(offset=offset,log_axial_scale=s,log_transverse_scale=-s/2,
                shape_velocity=velocity,phase_cost=potential,
                shape_kinetic=kinetic,total=potential+kinetic,
                log_volume=0.,receiving_child_formed=False)


def intrinsic_compression_calculation():
    """Small exact geometry controls; no selected source, pressure or SI units."""
    scales=(1.,.5,.25)
    dimensional=[]
    for dimension in (1,2,3):
        costs=[intrinsic_phase_cost(scale*np.eye(dimension),np.eye(dimension))['cost']
               for scale in scales]
        dimensional.append(dict(dimension=dimension,scales=list(scales),costs=costs,
                                scale_exponent=dimension-2,
                                isotropic_logarithmic_force_per_cost=(2-dimension)/dimension))
    shapes=[]
    for s in (0.,.5,1.,2.):
        A=np.diag([math.exp(-s/2),math.exp(-s/2),math.exp(s)])
        shapes.append(dict(log_axial_scale=s,
            directional_cost=intrinsic_phase_cost(A,[[0.,0.,1.]])['cost'],
            isotropic_cost=intrinsic_phase_cost(A,np.eye(3))['cost']))
    return dict(schema='nsc-intrinsic-compression-v1',
        premises=['unit target metric for the existing angular phase',
                  'intrinsic squared phase gradient integrated over actual volume',
                  'given reference phase-gradient moment',
                  'affine deformation controls; fixed-volume restriction for shape comparison'],
        cost='Q_g=integral |grad_g psi|^2 dV_g / 2',
        affine_cost='Q(A)=det(A) trace(A^-T M A^-1)/2',
        inherited_isotropic_cost_ratio='Q_child/Q_parent=exp(-(dimension-2)*chi_parent_child), for the same carried reference phase pattern',
        dimensional_controls=dimensional,
        shape_controls=shapes,
        positive_external_compression='For d=3 and uniform scale lambda: d[lambda Q0+P V0 lambda^3]/d lambda = Q0+3 P V0 lambda^2 > 0 for P>0,Q0>0; no finite isotropic minimum.',
        conservative_directional_control=[conservative_neck_shape(t) for t in (0.,1.,2.,4.)],
        no_rotation_torque_from_global_rotation=True,
        actual_parent_gradient_moment_selected=False,
        actual_parent_loading_selected=False,
        formed_child=False,
        no_observational_input=True)


def winding_calculation():
    """Geometry-only reference calculation: no parent action or observed scale."""
    cases = []
    for count in (2, 3, 4):
        directions = balanced_output_directions(count)
        cases.append(dict(channels=count, fraction_each=1/count,
                          normalized_vector_sum=(directions.sum(axis=0)/count).tolist(),
                          adjacent_angle_degrees=360/count))
    return dict(
        schema='nsc-offset-winding-geometry-v1',
        premises=['circle target with unit radius for this dimensionless example',
                  'integer winding of a continuous closed loop',
                  'concatenation joins opposite loops at one contact',
                  'squared-derivative functional in the stated drawing coordinate'],
        pair=dict(degrees=[1, -1], combined_degree=0,
                  variation_each=2*math.pi**2, variation_pair=4*math.pi**2,
                  conjugation_preserves_positive_variation=True),
        contraction=[dict(parameter=h, variation=winding_variation(h))
                     for h in (0., .25, .5, .75, 1.)],
        two_channel_consequence=dict(
            assumptions=['exactly two unit-direction outputs',
                         'positive scalar budgets with total one',
                         'zero weighted direction sum, with no other receiver'],
            fractions=[.5, .5], ratio=1., separation_degrees=180.),
        competing_balanced_channels=cases,
        missing_connection='contact composition/evolution law selecting the channels and transfer',
        physical_energy_unit_selected=False,
        physical_duration_selected=False,
        electron_or_positron_identification_derived=False,
        no_observational_input=True)


def offset_value(kind):
    if kind=='golden':return (3-mp.sqrt(5))/2
    if kind=='action':
        x=mp.findroot(lambda x:25*x**6-140*x+112,(mp.mpf('1.07'),mp.mpf('1.08')))
        return mp.sqrt((475*x-468)/(12800+75*x-84))
    raise ValueError('golden or retained action offset required')


def phase_history(count,kind='golden'):
    if type(count) is not int or count<1:raise ValueError('positive integer count required')
    with mp.workdps(80):
        alpha=offset_value(kind)
        return np.array([[float(mp.frac(n*alpha)),float(mp.frac(n*(n-1)*alpha/2))] for n in range(count)])


def directions(phases):
    theta=2*np.pi*phases[:,0];phi=2*np.pi*phases[:,1]
    return np.c_[np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)]


def euclidean_kernel_distance(r,s,theta,phi,theta2,phi2):
    C=math.cos(theta)*math.cos(theta2)+math.sin(theta)*math.sin(theta2)*math.cos(phi-phi2)
    return math.sqrt(max(0.,r*r+s*s-2*r*s*C))


def closed_parent_point(radial,direction,radius=1.):
    if not 0<=radial<=math.pi*radius:raise ValueError('north-to-south seed interval required')
    u=radial/radius
    return radius*np.r_[math.cos(u),math.sin(u)*np.asarray(direction)]


def sphere_cover(points):
    points=np.asarray(points)
    diagram=SphericalVoronoi(points,radius=1,center=np.zeros(3),threshold=1e-8)
    chord=cKDTree(points).query(diagram.vertices,workers=8)[0]
    return float(np.max(2*np.arcsin(np.clip(chord/2,0,1))))


def finite_record(count,kind):
    phases=phase_history(count,kind);v=directions(phases)
    line_cover=sphere_cover(np.r_[v,-v]);oriented_cover=sphere_cover(v)
    weights=np.abs(np.sin(2*np.pi*phases[:,0]))
    raw=(v.T@v)/count;weighted=(v.T@(v*weights[:,None]))/weights.sum()
    return dict(passes=count,
        line_direction_cover_radians=line_cover,
        ball_Hausdorff_distance_over_radius=math.sin(line_cover),
        oriented_direction_cover_radians=oriented_cover,
        S3_meridian_Hausdorff_distance_over_radius=min(oriented_cover,math.pi/2),
        raw_direction_second_moment=np.diag(raw).tolist(),
        area_weighted_direction_second_moment=np.diag(weighted).tolist(),
        finite_phase_Fourier_moments={f'{k},{l}':float(abs(np.mean(np.exp(2j*np.pi*(k*phases[:,0]+l*phases[:,1])))))
                                     for k,l in ((1,0),(0,1),(1,1),(2,-1))})


def calculation():
    # This is the retained action-based comparison, not a foundation input.
    from .nsc_atlas_parent_scale_selector import VolumeParent
    cases=[]
    for kind in ('action','golden'):
        with mp.workdps(80):alpha=mp.nstr(offset_value(kind),50)
        rows=[finite_record(n,kind) for n in (100,200,1000,4000)]
        cases.append(dict(offset=kind,alpha=alpha,angle=2*math.pi*float(alpha),coverage=rows))
    minimum=VolumeParent().vacuum_roots()[1];Lambda=minimum['vacuum'];R=math.sqrt(3/Lambda)
    return dict(schema='nsc-atlas-offset-generator-v1',
        generator=dict(seed='One continuous spatial interval; n counts retained offset updates, not a supplied proper-time coordinate.',
            rule='theta[n+1]=(theta[n]+alpha) mod1; phi[n+1]=(phi[n]+theta[n]) mod1; theta0=phi0=0',
            closed_form='theta_n=n alpha mod1; phi_n=n(n-1)alpha/2 mod1',
            repeated_operation='Delta theta=alpha; Delta phi=theta: the same discrete accumulation used twice.',
            inverse='theta_old=(theta_new-alpha) mod1; phi_old=(phi_new-theta_old) mod1',
            proof='Every nonzero integer character is a nonconstant polynomial in n with an irrational nonconstant coefficient; Weyl criterion gives joint equidistribution.',
            reversible=True,physical_duration_selected=False),
        metric_completion=dict(
            phase_angles='Theta=2pi theta; Phi=2pi phi',
            cosine_kernel='C_nm=cosTheta_n cosTheta_m+sinTheta_n sinTheta_m cos(Phi_n-Phi_m)',
            disk='Retain the seed interval rotated by Theta_n in one plane; its completion is the disk.',
            ball='d²=r²+s²-2rs C_nm; quotient coincident labels and complete: the Euclidean 3-ball.',
            closed_parent='d=R arccos[cos(r/R)cos(s/R)+sin(r/R)sin(s/R)C_nm], 0<=r,s<=piR; completion is round S3.',
            distinction='S3 has three intrinsic spatial dimensions despite a four-component unit-vector embedding. It is not a five-dimensional spacetime.',
            assumptions=['retain earlier drawing labels','use the stated cross-pass distance kernel','take the metric completion'],
            finite_history_has_nonzero_volume=False,
            travel_only_along_spokes_gives_same_completion=False,
            unique_three_dimension_selection=False),
        measure=dict(raw_phase_limit_second_moment=[.25,.25,.5],
            sphere_area_weight='abs(sinTheta)',area_weighted_limit_second_moment=[1/3]*3,
            ball_volume='r² abs(sinTheta) dr dTheta dPhi',
            closed_volume='R² sin²(r/R) abs(sinTheta) dr dTheta dPhi',
            physical_energy_measure_selected=False),
        cases=cases,
        inherited_parent=dict(Lambda=Lambda,radius=R,scalar_curvature=6/R**2,
            Hamiltonian_residual=6/R**2-2*Lambda,extrinsic_curvature=0,
            first_scalar_Laplacian_eigenvalues=[3/R**2,8/R**2],
            normalized_eigenvalues=[3,8],
            law='Under radius rescaling, R²(-Delta_S3_R) has the same dimensionless spectrum.',
            scope='The generated S3 can carry the existing constant-field, time-symmetric parent data. This is compatible initial geometry, not a selected formation history.'),
        physical_open=['Why the offset dynamics and metric kernel are selected',
                       'Why two accumulation levels followed by nesting are physically selected',
                       'Regional scale/energy transfer, a duration readout and an arrow',
                       'Birth of a daughter and observational abundances'],
        no_observational_fit=True)
