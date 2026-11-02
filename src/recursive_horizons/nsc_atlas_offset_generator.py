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
        bookkeeping_note='chi additivity is an algebraic log-ratio telescoping identity (bookkeeping), not dynamical evidence for the mechanism',
        mechanism_evidence=False,
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


def generator_least_stretch_cost(deformation, gradient_factor, enclosing_load=0., reference_volume=1.):
    """A0 least-stretch cost of 1D generated phase paths under affine deformation.

    Unlike the 3D volume-weighted Dirichlet integral det(A)*||C A^-1||^2/2
    (which scales as lambda^(d-2) and decreases to 0 under d=3 compression),
    a 1D generated curve of reference length L carrying phase change Delta_psi
    has Cauchy-Schwarz minimum arc-length stretch (Delta_psi)^2/(2 L_physical).
    For right-stretch (A^T A)^(1/2) and 1D moment M = C^T C / 2, this gives
    Q_stretch(A) = trace[(A^T A)^(-1/2) C^T C] / 2. Under uniform scaling
    A = lambda I = exp(-chi) I, Q_stretch scales as lambda^-1 = exp(chi) in
    every dimension d >= 1, resisting compression rather than vanishing.
    """
    A=np.asarray(deformation,dtype=float);C=np.asarray(gradient_factor,dtype=float)
    if (A.ndim!=2 or A.shape[0]!=A.shape[1] or A.shape[0]<1
            or C.ndim!=2 or C.shape[0]<1 or C.shape[1]!=A.shape[0]):
        raise ValueError('square deformation and matching gradient factor required')
    if not np.all(np.isfinite(A)) or not np.all(np.isfinite(C)):
        raise ValueError('finite geometric inputs required')
    if not math.isfinite(enclosing_load) or enclosing_load<0 or not math.isfinite(reference_volume) or reference_volume<=0:
        raise ValueError('nonnegative finite enclosing load and positive reference volume required')
    sign,log_volume=np.linalg.slogdet(A)
    if sign<=0:raise ValueError('orientation-preserving nonsingular deformation required')
    vals,vecs=np.linalg.eigh(A.T@A)
    if np.any(vals<=0):raise ValueError('positive definite metric required')
    stretches=np.sqrt(vals)
    projected=C@vecs
    directional_moments=.5*np.sum(projected**2,axis=0)
    stretch_cost=float(np.sum(directional_moments/stretches))
    enclosing_work=float(enclosing_load*reference_volume*math.exp(log_volume))
    return dict(stretch_cost=stretch_cost,enclosing_work=enclosing_work,
                total_cost=stretch_cost+enclosing_work,log_volume=float(log_volume),
                principal_stretches=stretches.tolist(),
                principal_moments=directional_moments.tolist())


def least_resistance_compression_minimum(principal_moments, enclosing_load, reference_volume=1.):
    """Exact minimum of Q_stretch + P V under enclosing compression P > 0.

    For 1D generator moments m_i > 0 along d=3 principal axes:
    E(lambda_1,lambda_2,lambda_3) = sum_i m_i/lambda_i + P V0 lambda_1 lambda_2 lambda_3.
    Stationarity dE/d(lambda_i)=0 gives lambda_i^* = m_i / [(m_1 m_2 m_3)^(1/4) (P V0)^(1/4)].
    If any transverse moment vanishes (1D seed limit m_perp -> 0), the transverse
    radius pinches to 0 (s -> +infinity), triggering the neck-break continuation.
    """
    m=np.asarray(principal_moments,dtype=float)
    if m.shape!=(3,) or not np.all(np.isfinite(m)) or np.any(m<0) or np.all(m==0):
        raise ValueError('three finite nonnegative principal moments with positive sum required')
    if not math.isfinite(enclosing_load) or enclosing_load<=0 or not math.isfinite(reference_volume) or reference_volume<=0:
        raise ValueError('positive finite enclosing load and reference volume required')
    pv=enclosing_load*reference_volume
    if np.any(m==0):
        return dict(outcome='neck_pinch_opening_break',
                    principal_moments=m.tolist(),enclosing_load_work_scale=pv,
                    equilibrium_scales=None,
                    break_reason='zero transverse generator moment cannot hold transverse radius against enclosing compression; continues across neck pinch')
    geom=float(np.prod(m))**.25
    load_factor=pv**.25
    scales=m/(geom*load_factor)
    volume_scale=float(np.prod(scales))
    stretch_cost=float(np.sum(m/scales))
    enclosing_work=float(pv*volume_scale)
    s_neck=float(math.log(scales[2]/math.sqrt(scales[0]*scales[1]))/3)
    outcome='uniform_contraction' if np.allclose(m,m[0]) else 'anisotropic_neck_equilibrium'
    return dict(outcome=outcome,
                principal_moments=m.tolist(),
                enclosing_load_work_scale=pv,
                equilibrium_scales=scales.tolist(),
                equilibrium_volume_ratio=volume_scale,
                inherited_compression_chi=-math.log(volume_scale)/3,
                log_axial_neck_parameter=s_neck,
                axial_to_transverse_aspect_ratio=float(scales[2]/math.sqrt(scales[0]*scales[1])),
                stretch_cost=stretch_cost,
                enclosing_work=enclosing_work,
                total_cost=stretch_cost+enclosing_work,
                virial_ratio_stretch_to_work=stretch_cost/enclosing_work)


def continue_across_neck_break(offset):
    """Apply the NSC break rule to the axial neck pinch a(tau)=exp(-s/2)->0.

    In parent coordinates, conservative_neck_shape drives transverse radius
    a(tau)=[cosh(sqrt(2/3)*tau)]^(-1/2) -> 0 and normalized velocity
    s_dot/sqrt(2/3) = tanh(sqrt(2/3)*tau) -> 1 (hitting 0.999889... by tau=6).
    Changing coordinates at the throat via R_child = a^(-2) = exp(s) = cosh(sqrt(2/3)*tau)
    continues the same conserved trajectory into an expanding child frame obeying
    (R_child_dot / R_child)^2 + 2/(3 R_child^2) = 2/3, with asymptotic expansion
    rate H_child -> sqrt(2/3) and effective child cosmological constant Lambda_child = 2.
    """
    parent=conservative_neck_shape(offset)
    a_rate=math.sqrt(2/3);u=a_rate*offset
    s=parent['log_axial_scale']
    parent_transverse=math.exp(-s/2)
    wall_fraction=math.tanh(u)
    child_scale=math.exp(s)
    child_hubble=parent['shape_velocity']
    child_velocity=child_scale*child_hubble
    friedmann_lhs=child_hubble**2+2/(3*child_scale**2)
    return dict(offset=offset,
                break_type='transverse_throat_pinch_and_velocity_wall',
                parent_frame=dict(
                    log_axial_scale=s,
                    transverse_radius=parent_transverse,
                    phase_cost=parent['phase_cost'],
                    shape_kinetic=parent['shape_kinetic'],
                    velocity_fraction_of_asymptotic_wall=wall_fraction),
                coordinate_change='R_child = a_parent^(-2) = exp(s); chi_child = log(R_child) = s',
                child_frame=dict(
                    scale_factor=child_scale,
                    expansion_log_chi=s,
                    scale_velocity=child_velocity,
                    hubble_rate=child_hubble,
                    asymptotic_hubble_rate=a_rate,
                    effective_child_lambda=3*a_rate**2,
                    friedmann_balance=friedmann_lhs,
                    friedmann_target=2/3,
                    friedmann_residual=friedmann_lhs-2/3,
                    transferred_budget_fraction=parent['shape_kinetic']/parent['total']))


def spiral_decay_to_turn_comparison(alpha=None):
    """Evaluate turning vs compression and compare decay-to-turn ratios without imposing phi.

    Both Q_vol(A) and Q_stretch(A) depend only on symmetric strain and are invariant
    under rigid rotations A -> R A, so pure static compression produces zero torque.
    When unit logarithmic compression lambda(chi)=exp(-chi) is paired with a supplied
    or inherited angular rate omega = 2*pi*alpha per unit chi, z(chi)=exp((-1+2*pi*i*alpha)*chi)
    has decay-per-radian gamma/omega = 1/(2*pi*alpha) and opposite-lobe asymmetry
    eta = (A_+ - |A_-|)/(A_+ + |A_-|) = tanh(1/(4*alpha)).
    """
    phi=(1+math.sqrt(5))/2
    golden_spiral_ratio=2*math.log(phi)/math.pi
    alpha_for_golden_spiral=1/(4*math.log(phi))
    alpha_golden_angle=float(offset_value('golden'))
    alpha_action=float(offset_value('action'))
    chosen=alpha_golden_angle if alpha is None else float(alpha)
    if not math.isfinite(chosen) or chosen<=0:
        raise ValueError('positive finite angular increment required for spiral comparison')

    def metrics(a):
        decay_per_radian=1/(2*math.pi*a)
        lobe_ratio=math.exp(-1/(2*a))
        eta=math.tanh(1/(4*a))
        return dict(alpha=a,
                    decay_per_radian=decay_per_radian,
                    decay_per_full_turn=1/a,
                    opposite_lobe_amplitude_ratio=lobe_ratio,
                    lobe_asymmetry_eta=eta,
                    difference_from_golden_spiral_ratio=decay_per_radian-golden_spiral_ratio)

    quadrupole_omega=math.sqrt(2/3)
    return dict(static_compression_produces_rotation=False,
                reason_no_torque='Q_vol and Q_stretch are SO(3)-invariant functions of symmetric stretch; turning requires rotational/chiral offset input or angular momentum',
                golden_spiral_decay_per_radian=golden_spiral_ratio,
                alpha_required_for_golden_spiral=alpha_for_golden_spiral,
                quadrupole_shape_frequency=quadrupole_omega,
                compression_to_quadrupole_frequency_ratio=1/quadrupole_omega,
                evaluated=metrics(chosen),
                golden_angle_increment=metrics(alpha_golden_angle),
                retained_action_increment=metrics(alpha_action))


def helical_buckling_under_compression(compression_chi, stretch_budget=1., lateral_confinement=1., radius=1.):
    """Derive turning (omega* > 0) by least resistance when axial compression exceeds threshold.

    A 1D seed compressed axially to scale lambda = exp(-chi) that turns laterally
    at angular rate omega around radius r_perp has Euclidean arc-length factor
    l = sqrt(lambda^2 + (r_perp*omega)^2). Its Cauchy-Schwarz least-stretch cost
    plus lateral enclosing confinement is
    E(v_perp; lambda) = Q0 / sqrt(lambda^2 + v_perp^2) + P_perp v_perp^2 / 2,
    where v_perp = r_perp * omega.
    At lambda_crit = (Q0 / P_perp)^(1/3) (chi_crit = -log(Q0/P_perp)/3), the straight
    seed (omega=0) loses transverse stability (d^2E/dv_perp^2 = P_perp - Q0/lambda^3 < 0)
    and buckles into a right- or left-handed helix with exact turning speed
    v_perp^* = sqrt(lambda_crit^2 - lambda^2) and derived turn increment alpha^* = omega^*/(2*pi).
    """
    if (not math.isfinite(compression_chi) or not math.isfinite(stretch_budget) or stretch_budget<=0
            or not math.isfinite(lateral_confinement) or lateral_confinement<=0
            or not math.isfinite(radius) or radius<=0):
        raise ValueError('finite compression and positive stretch_budget, lateral_confinement and radius required')
    lam=math.exp(-compression_chi)
    lam_crit=(stretch_budget/lateral_confinement)**(1/3)
    chi_crit=-math.log(stretch_budget/lateral_confinement)/3
    straight_cost=stretch_budget/lam
    straight_hessian=lateral_confinement-stretch_budget/lam**3
    if lam>=lam_crit:
        return dict(compression_chi=compression_chi,axial_scale=lam,
                    critical_compression_chi=chi_crit,critical_axial_scale=lam_crit,
                    buckled=False,straight_transverse_hessian=straight_hessian,
                    turning_speed=0.,angular_rate_omega=0.,derived_alpha=0.,
                    helical_arc_length=lam,straight_cost=straight_cost,
                    minimum_cost=straight_cost,cost_relief=0.)
    v_perp=math.sqrt(lam_crit**2-lam**2)
    omega=v_perp/radius
    alpha_derived=omega/(2*math.pi)
    min_cost=1.5*lateral_confinement*lam_crit**2-.5*lateral_confinement*lam**2
    return dict(compression_chi=compression_chi,axial_scale=lam,
                critical_compression_chi=chi_crit,critical_axial_scale=lam_crit,
                buckled=True,straight_transverse_hessian=straight_hessian,
                turning_speed=v_perp,angular_rate_omega=omega,derived_alpha=alpha_derived,
                helical_arc_length=lam_crit,straight_cost=straight_cost,
                minimum_cost=min_cost,cost_relief=straight_cost-min_cost)


def damped_horizon_oscillation(offset, damping_rate, turn_frequency, amplitude=1., enclosing_load=1.):
    """Exact two-sided (+/0/-) horizon wave coupled to child geometry expansion.

    Solves q_tt + 2*gamma*q_t + (gamma^2 + omega^2)*q = 0 with q(0)=0, q_t(0)=A0*omega,
    so q(tau) = A0 * exp(-gamma*tau) * sin(omega*tau), coupled to the receiving
    geometry budget by E_geom_t = 2*gamma*q_t^2. Total budget E_wave(tau) + E_geom(tau)
    is identically conserved at E_total = A0^2 * omega^2 / 2.
    Because each + lobe precedes its compensating - lobe by pi/omega, the peak ratio,
    single-cycle area ratio, and infinite cumulative area ratio all yield the same
    invariant asymmetry eta = (I_+ - I_-)/(I_+ + I_-) = tanh(pi*gamma/(2*omega)) > 0,
    while the child's least-resistance holding scale expands to the Goldilocks
    equilibrium radius lambda_Goldilocks = (E_total / (3*P_ext))^(1/4).
    """
    if (not math.isfinite(offset) or offset<0 or not math.isfinite(damping_rate) or damping_rate<0
            or not math.isfinite(turn_frequency) or turn_frequency<=0
            or not math.isfinite(amplitude) or amplitude<=0
            or not math.isfinite(enclosing_load) or enclosing_load<=0):
        raise ValueError('nonnegative finite offset/damping and positive frequency, amplitude and load required')
    g=damping_rate;w=turn_frequency;w0_sq=g*g+w*w
    env=amplitude*math.exp(-g*offset)
    q=env*math.sin(w*offset)
    q_dot=env*(w*math.cos(w*offset)-g*math.sin(w*offset))
    e_total=.5*amplitude**2*w**2
    e_wave=.5*q_dot**2+.5*w0_sq*q**2
    e_geom=max(0.,e_total-e_wave)
    rho=g/w
    exp_half=math.exp(-math.pi*rho)
    first_plus=amplitude*w/w0_sq*(1+exp_half)
    first_minus=first_plus*exp_half
    eta=math.tanh(.5*math.pi*rho)
    if g>0:
        total_plus=first_plus/(1-exp_half**2)
        total_minus=total_plus*exp_half
        net_remainder=amplitude*w/w0_sq
    else:
        total_plus=None;total_minus=None;net_remainder=0.
    child_scale=(e_geom/(3*enclosing_load))**.25 if e_geom>0 else 0.
    goldilocks_scale=(e_total/(3*enclosing_load))**.25
    return dict(offset=offset,displacement=q,velocity=q_dot,envelope=env,
                wave_budget=e_wave,transferred_geometry_budget=e_geom,
                total_budget=e_total,conservation_residual=(e_wave+e_geom)-e_total,
                first_plus_lobe_area=first_plus,first_minus_lobe_area=first_minus,
                cumulative_plus_area=total_plus,cumulative_minus_area=total_minus,
                net_asymmetric_remainder=net_remainder,lobe_asymmetry_eta=eta,
                child_holding_scale=child_scale,
                goldilocks_equilibrium_scale=goldilocks_scale,
                goldilocks_completion_fraction=child_scale/goldilocks_scale)


def two_force_zero_crossing(offset, inward_compression_force, outward_resistance_force, initial_radius=1., initial_velocity=0.):
    """Direct + -> 0 -> - spherical zero-crossing when inward compression wins.

    When inward compression F_in exceeds outward resistance F_out on a sphere of
    signed coordinate R(tau) with R(0)=R0 > 0, the net force F_net = F_in - F_out > 0
    drives R_tt = -F_net < 0:
      R(tau) = R0 + v0*tau - F_net*tau^2/2.
    Comparing the three candidate paths at the horizon R = 0:
      1. Busting back into the parent + sphere (dR > 0, R > 0) requires climbing
         uphill against the winning force F_net > 0 (positive work barrier).
      2. Stopping at R = 0 requires an infinite deceleration impulse to absorb
         kinetic energy K(tau_0) > 0 in zero distance.
      3. Crossing 0 into R < 0 (- region) is smooth (zero impulse, finite R_tt),
         flips the 3D Jacobian determinant det(R/R0 * I_3) = (R/R0)^3 from + to -
         through 0, and converts the inward implosion velocity R_dot < 0 directly
         into the outward expansion velocity r_minus_dot = -R_dot > 0 of a new
         sphere of physical radius r_minus = |R| = -R growing from scratch.
    """
    if (not math.isfinite(offset) or offset<0
            or not math.isfinite(inward_compression_force)
            or not math.isfinite(outward_resistance_force)
            or not math.isfinite(initial_radius) or initial_radius<=0
            or not math.isfinite(initial_velocity) or initial_velocity>0):
        raise ValueError('nonnegative offset, finite forces, positive initial_radius and nonpositive initial_velocity required')
    f_net=inward_compression_force-outward_resistance_force
    if f_net<=0:
        raise ValueError('inward_compression_force must exceed outward_resistance_force for collapse through zero')
    r_signed=initial_radius+initial_velocity*offset-.5*f_net*offset**2
    r_dot=initial_velocity-f_net*offset
    horizon_offset=(initial_velocity+math.sqrt(initial_velocity**2+2*f_net*initial_radius))/f_net
    horizon_speed=math.sqrt(initial_velocity**2+2*f_net*initial_radius)
    horizon_kinetic=.5*horizon_speed**2
    det_ratio=(r_signed/initial_radius)**3
    if r_signed>0:
        region='+'
        physical_radius=r_signed
        radial_expansion_rate=r_dot
    elif r_signed==0:
        region='0'
        physical_radius=0.
        radial_expansion_rate=0.
    else:
        region='-'
        physical_radius=-r_signed
        radial_expansion_rate=-r_dot
    potential=f_net*r_signed
    kinetic=.5*r_dot**2
    initial_total=.5*initial_velocity**2+f_net*initial_radius
    return dict(offset=offset,
                net_compression_force=f_net,
                signed_radius=r_signed,
                signed_velocity=r_dot,
                region=region,
                physical_radius=physical_radius,
                physical_expansion_rate=radial_expansion_rate,
                jacobian_determinant_ratio=det_ratio,
                orientation_sign=1 if r_signed>0 else (0 if r_signed==0 else -1),
                horizon_crossing_offset=horizon_offset,
                horizon_crossing_speed=horizon_speed,
                horizon_kinetic_budget=horizon_kinetic,
                outward_parent_reversal_blocked=True,
                stopping_at_zero_requires_infinite_deceleration=True,
                zero_crossing_impulse_cost=0.,
                unique_finite_cost_continuation='cross 0 into R < 0 (- region), where r_minus = -R expands from 0 with r_minus_dot = -R_dot > 0',
                total_mechanical_budget=kinetic+potential,
                conservation_residual=(kinetic+potential)-initial_total)


def buckling_to_goldilocks_calculation():
    """Connected pipeline: two-force 0-crossing -> helical buckling -> damped +/0/- wave -> Goldilocks zone."""
    t_h=math.sqrt(2.)
    zero_crossing_scan=[two_force_zero_crossing(t,2.,1.) for t in (0.,1.,t_h,2.,3.)]
    buckling_scan=[helical_buckling_under_compression(chi) for chi in (0.,.25,.5,1.,2.)]
    post_buckling=helical_buckling_under_compression(1.)
    omega=2*math.pi
    gamma=.5*math.log(2)
    wave_history=[damped_horizon_oscillation(t,gamma,omega) for t in (0.,.25,.75,1.,2.,4.)]
    return dict(schema='nsc-buckling-damping-goldilocks-v1',
                stage_0_two_force_zero_crossing=dict(
                    rule='R_tt = -(F_in - F_out) < 0; R(tau) crosses + -> 0 -> - smoothly with det(F)=(R/R0)^3 flipping sign',
                    horizon_offset=t_h,
                    scan=zero_crossing_scan),
                stage_1_helical_buckling=dict(
                    rule='E(v_perp; lambda) = Q0/sqrt(lambda^2 + v_perp^2) + P_perp v_perp^2 / 2',
                    critical_condition='lambda_crit = (Q0/P_perp)^(1/3); omega^* = sqrt(lambda_crit^2 - lambda^2)/r_perp for lambda < lambda_crit',
                    scan=buckling_scan,
                    derived_turning_at_chi_1=post_buckling['angular_rate_omega'],
                    derived_alpha_at_chi_1=post_buckling['derived_alpha']),
                stage_2_damped_horizon_asymmetry=dict(
                    equation='q_tt + 2*gamma*q_t + (gamma^2+omega^2)*q = 0; E_geom_t = 2*gamma*q_t^2',
                    invariant_asymmetry_formula='eta = (I_+ - I_-)/(I_+ + I_-) = (I_0^+ - |I_1^-|)/(I_0^+ + |I_1^-|) = tanh(pi*gamma/(2*omega))',
                    net_remainder_formula='Delta_I = I_+ - I_- = A0 * omega / (gamma^2 + omega^2)',
                    history=wave_history),
                stage_3_goldilocks_equilibrium=dict(
                    formula='lambda_child(tau) = (E_geom(tau)/(3 P V0))^(1/4) -> lambda_Goldilocks = (E_total/(3 P V0))^(1/4)',
                    asymptotic_scale=wave_history[-1]['goldilocks_equilibrium_scale'],
                    completion_at_offset_4=wave_history[-1]['goldilocks_completion_fraction']),
                no_observational_input=True)


def intrinsic_compression_calculation():
    """Small exact geometry controls; no selected source, pressure or SI units."""
    scales=(1.,.5,.25)
    dimensional=[]
    for dimension in (1,2,3):
        costs=[intrinsic_phase_cost(scale*np.eye(dimension),np.eye(dimension))['cost']
               for scale in scales]
        stretch_costs=[generator_least_stretch_cost(scale*np.eye(dimension),np.eye(dimension))['stretch_cost']
                       for scale in scales]
        dimensional.append(dict(dimension=dimension,scales=list(scales),costs=costs,
                                generator_least_stretch_costs=stretch_costs,
                                scale_exponent=dimension-2,
                                generator_least_stretch_exponent=-1,
                                isotropic_logarithmic_force_per_cost=(2-dimension)/dimension))
    shapes=[]
    for s in (0.,.5,1.,2.):
        A=np.diag([math.exp(-s/2),math.exp(-s/2),math.exp(s)])
        shapes.append(dict(log_axial_scale=s,
            directional_cost=intrinsic_phase_cost(A,[[0.,0.,1.]])['cost'],
            isotropic_cost=intrinsic_phase_cost(A,np.eye(3))['cost'],
            generator_directional_stretch_cost=generator_least_stretch_cost(A,[[0.,0.,1.]])['stretch_cost'],
            generator_isotropic_stretch_cost=generator_least_stretch_cost(A,np.eye(3))['stretch_cost']))
    outcomes=dict(
        volume_weighted_isotropic_under_load=dict(
            outcome='collapses_to_zero_scale_break',
            equation='d[Q0 lambda + P V0 lambda^3]/d lambda = Q0 + 3 P V0 lambda^2 > 0',
            break_continuation='invert scale lambda_child = lambda^(-1) = exp(chi); parent zero-volume limit becomes child expansion domain (1, infinity)'),
        generator_least_stretch_isotropic_under_load=least_resistance_compression_minimum([.5,.5,.5],1.),
        generator_least_stretch_raw_offset_moment_under_load=least_resistance_compression_minimum([.25,.25,.5],1.),
        generator_least_stretch_1d_seed_under_load=least_resistance_compression_minimum([0.,0.,.5],1.),
        neck_break_continuation=[continue_across_neck_break(t) for t in (0.,1.,2.,4.,6.)],
        turning_and_spiral_check=spiral_decay_to_turn_comparison())
    return dict(schema='nsc-intrinsic-compression-v1',
        premises=['unit target metric for the existing angular phase',
                  'comparison of 3D volume-weighted Dirichlet cost against A0 1D generator least-stretch cost',
                  'given reference phase-gradient moment',
                  'affine deformation controls; enclosing compression load P > 0'],
        cost='Q_g=integral |grad_g psi|^2 dV_g / 2',
        affine_cost='Q(A)=det(A) trace(A^-T M A^-1)/2',
        generator_least_stretch_cost='Q_stretch(A)=trace[(A^T A)^(-1/2) C^T C]/2, derived by Cauchy-Schwarz on 1D generated phase paths',
        cost_justification='Q_g dilutes 1D phase stretch by 3D volume lambda^3 so cost falls as exp(-chi) in d=3; A0 least-stretch of the 1D generated curves scales as exp(+chi) and resists compression',
        inherited_isotropic_cost_ratio='Q_child/Q_parent=exp(-(dimension-2)*chi_parent_child), for the same carried reference phase pattern',
        inherited_generator_stretch_ratio='Q_stretch_child/Q_stretch_parent=exp(+chi_parent_child) in every dimension d',
        dimensional_controls=dimensional,
        shape_controls=shapes,
        positive_external_compression='For d=3 and uniform scale lambda: d[lambda Q0+P V0 lambda^3]/d lambda = Q0+3 P V0 lambda^2 > 0 for P>0,Q0>0; no finite isotropic minimum under volume-weighted Q_g, whereas Q_stretch=Q0/lambda gives unique minimum lambda_*=(Q0/(3 P V0))^(1/4).',
        conservative_directional_control=[conservative_neck_shape(t) for t in (0.,1.,2.,4.)],
        least_resistance_outcomes=outcomes,
        coupling_c_status='c along the 1D drawing coordinate remains supplied; 3D neck shape wave speed sqrt(2/3) follows only after declaring the log-strain kinetic metric.',
        increment_alpha_status='alpha remains supplied (or inherited from external rotational/chiral parent history) because symmetric compression costs are SO(3)-invariant and produce zero torque.',
        chi_bookkeeping_note='chi_gc = chi_gp + chi_pc is an algebraic log-ratio identity (bookkeeping), not dynamical evidence for the mechanism.',
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
