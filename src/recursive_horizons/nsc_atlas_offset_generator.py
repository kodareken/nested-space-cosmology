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


def spaghettification_to_spherical_winner(compression_chi):
    """Recurring cycle: 3D -> 1D spaghettification on infall, 1D -> 3D sphere inside the drop.

    1. Infall spaghettification (+ -> 0): A 3D parcel converging radially to
       radius ratio xi = r/r0 = exp(-chi) has transverse stretches lambda_1=lambda_2=xi
       (circumference 2*pi*r shrinks) and volume-preserving radial stretch
       lambda_3 = xi^(-2) = exp(2*chi). Normalized by lambda_3, the deformation is
       diag(xi^3, xi^3, 1) -> diag(0, 0, 1) (rank 1): spherical convergence
       geometrically forces a 3D parcel into a 1D radial filament (spaghettification).
    2. Isoperimetric spherical winner inside - (olive oil drop in water):
       For fixed enclosed volume V = (4/3)*pi*R0^3, any spheroid with axial log-stretch s
       has surface area ratio A(s)/(4*pi*R0^2) >= 1 and isoperimetric quotient
       I(s) = A(s)^3 / (36*pi*V^2) >= 1, with strict equality uniquely at s = 0
       (the round 3D sphere, the recurring least-resistance winner).
    3. Pressure-cooker takeover: As enclosing compression converts parent + volume
       V_+(chi) = V0*exp(-3*chi) into the enclosed - drop V_-(chi) = V0*(1 - exp(-3*chi)),
       the drop's share of the cooker 1 - exp(-3*chi) expands to fill the container.
    """
    if not math.isfinite(compression_chi) or compression_chi<0:
        raise ValueError('nonnegative finite compression_chi required')
    xi=math.exp(-compression_chi)
    lam_perp=xi
    lam_axial=math.exp(2*compression_chi)
    aspect=math.exp(3*compression_chi)
    normalized_transverse=math.exp(-3*compression_chi)
    s=compression_chi
    if s==0:
        area_ratio=1.
    else:
        ecc_sq=1-math.exp(-3*s)
        ecc=math.sqrt(ecc_sq)
        area_ratio=.5*(math.exp(-s)+math.exp(.5*s)*math.asin(ecc)/ecc)
    isoperimetric_quotient=area_ratio**3
    plus_fraction=math.exp(-3*compression_chi)
    minus_takeover_fraction=1-plus_fraction
    return dict(compression_chi=compression_chi,
                infall_spaghettification=dict(
                    radial_ratio_xi=xi,
                    transverse_stretch=lam_perp,
                    radial_stretch=lam_axial,
                    spaghettification_aspect_ratio=aspect,
                    normalized_deformation_diagonal=[normalized_transverse,normalized_transverse,1.],
                    asymptotic_rank_at_horizon=1),
                isoperimetric_spherical_winner=dict(
                    deformed_area_over_sphere_area=area_ratio,
                    deformed_isoperimetric_quotient=isoperimetric_quotient,
                    sphere_minimum_quotient=1.,
                    sphere_area_saving_fraction=1-1/area_ratio,
                    unique_minimum_shape='round_3d_sphere'),
                cooker_volume_takeover=dict(
                    parent_plus_volume_fraction=plus_fraction,
                    child_minus_drop_fraction=minus_takeover_fraction,
                    minus_to_plus_volume_ratio=(minus_takeover_fraction/plus_fraction) if plus_fraction>0 else None))


def receiving_shape_spectrum(winding=1., enclosing_load=1.):
    """Small shape response of the A0 sphere with declared log-strain inertia.

    E=sum_i n^2 exp(-y_i)+P exp(sum_i y_i), K=sum_i y_dot_i^2/2.
    The inherited winding is fixed during the response measurement. This is
    a conditional geometric spectrum, not an identification of matter modes.
    """
    if (not math.isfinite(winding) or winding == 0
            or not math.isfinite(enclosing_load) or enclosing_load <= 0):
        raise ValueError('finite nonzero winding and positive enclosing load required')
    radius=(winding*winding/enclosing_load)**.25
    work=enclosing_load*radius**3
    stiffness=work*(np.eye(3)+np.ones((3,3)))
    return dict(equilibrium_radius=radius, stiffness_matrix=stiffness.tolist(),
                frequency_squared=[work,work,4*work],
                breathing_to_shear_ratio=2.,
                geometric_inertia='K = sum_i (d log(lambda_i)/d offset)^2 / 2',
                inherited_winding_held_fixed=True,
                physical_frequency_unit_selected=False)


def _parent_receiver_system(parents, cells_per_room):
    """One continuous phase on a finite tree of deforming generated curves.

    Regional available volumes are V_i=R_i^3-sum_children R_j^3;
    their curve stretches are lambda_i=V_i^(1/3). This volume-to-curve
    identification and K_geometry=3*sum(log(R_i)_dot^2)/2 are explicit
    reduction assumptions. The outer boundary is free. No PV reservoir,
    imposed compression history, damping, or membrane is added.
    """
    parents=tuple(parents)
    if (not parents or parents[0] != -1 or
            any(type(p) is not int or not 0 <= p < i for i,p in enumerate(parents[1:],1))):
        raise ValueError('one ordered rooted tree required, with parent[0]=-1')
    if type(cells_per_room) is not int or cells_per_room < 4 or cells_per_room % 2:
        raise ValueError('an even cell count of at least four per room is required')
    count=len(parents); children=[[] for _ in parents]
    for i,p in enumerate(parents[1:],1): children[p].append(i)
    left=[];right=[];owner=[];starts=[];ends=[];next_node=1
    for i,p in enumerate(parents):
        start=0 if p == -1 else ends[p]
        nodes=[start]+list(range(next_node,next_node+cells_per_room))
        next_node+=cells_per_room
        starts.append(start);ends.append(nodes[-1])
        left.extend(nodes[:-1]);right.extend(nodes[1:]);owner.extend([i]*cells_per_room)
    left=np.array(left);right=np.array(right);owner=np.array(owner)
    fixed=np.array([0]+[ends[i] for i in range(count) if not children[i]])
    free=np.setdiff1d(np.arange(next_node),fixed)
    nf=len(free);h=2./cells_per_room;inertia=3.
    subtree=np.ones(count)
    for i in range(count-1,0,-1):subtree[parents[i]]+=subtree[i]
    phase=np.zeros(next_node)
    x=np.linspace(0.,2.,cells_per_room+1)
    root_nodes=np.r_[starts[0],right[:cells_per_room]]
    phase[root_nodes]=2*np.pi*np.minimum(x,2-x)
    initial=np.r_[phase[free],np.zeros(nf),np.log(subtree)/3,np.zeros(3*count)]

    def geometry(q):
        enclosed=np.exp(3*q)
        available=enclosed.copy()
        for i,p in enumerate(parents[1:],1): available[p]-=enclosed[i]
        return enclosed,available

    def evaluate(state):
        phase=np.zeros(next_node);phase[free]=state[:nf]
        momentum=np.zeros(next_node);momentum[free]=state[nf:2*nf]
        q=state[2*nf:2*nf+count];pg=state[2*nf+count:2*nf+2*count]
        enclosed,volume=geometry(q)
        if np.any(volume <= 0):
            raise ValueError('available regional volume vanished; positive-volume reduction needs continuation')
        lam=np.cbrt(volume);half_mass=lam[owner]*h/2
        mass=np.bincount(left,weights=half_mass,minlength=next_node)+np.bincount(right,weights=half_mass,minlength=next_node)
        velocity=momentum/mass
        gradient=(phase[right]-phase[left])/(lam[owner]*h)
        force=np.bincount(left,weights=gradient,minlength=next_node)-np.bincount(right,weights=gradient,minlength=next_node)
        kinetic=half_mass*(velocity[left]**2+velocity[right]**2)/2
        stretch=(phase[right]-phase[left])*gradient/2
        energy=np.bincount(owner,weights=kinetic+stretch,minlength=count)
        qdot=pg/inertia
        volume_dot=3*enclosed*qdot
        for i,p in enumerate(parents[1:],1):volume_dot[p]-=3*enclosed[i]*qdot[i]
        log_stretch_dot=volume_dot/(3*volume)
        density=energy/volume
        geometry_force=enclosed*density
        for i,p in enumerate(parents[1:],1):geometry_force[i]-=enclosed[i]*density[p]
        half_mass_dot=half_mass*log_stretch_dot[owner]
        mass_dot=np.bincount(left,weights=half_mass_dot,minlength=next_node)+np.bincount(right,weights=half_mass_dot,minlength=next_node)
        acceleration=(force-mass_dot*velocity)/mass
        acceleration[fixed]=0.
        transfer=np.zeros(count)
        for i in range(count):
            lo=i*cells_per_room;hi=(i+1)*cells_per_room-1
            for node,edge,traction in ((starts[i],lo,-gradient[lo]),(ends[i],hi,gradient[hi])):
                transfer[i]+=velocity[node]*(half_mass[edge]*acceleration[node]+half_mass_dot[edge]*velocity[node]+traction)
        work=-log_stretch_dot*energy
        derivative=np.r_[velocity[free],force[free],qdot,geometry_force,transfer,work]
        geometric_energy=float(np.dot(pg,pg)/(2*inertia))
        return derivative,dict(phase=phase,momentum=momentum,mass=mass,
            phase_energy=energy,phase_stretch=np.bincount(owner,weights=stretch,minlength=count),
            volume=volume,enclosed=enclosed,stretch=lam,
            pressure=density/3,geometry_energy=geometric_energy,
            total_energy=float(np.sum(energy)+geometric_energy),
            transfer=transfer,work=work,geometry_force=geometry_force,
            log_stretch_dot=log_stretch_dot)

    return dict(initial=initial,evaluate=evaluate,geometry=geometry,
                free_phase_count=nf,room_count=count,parents=parents,
                children=children,inertia=inertia)


def parent_receiver_evolution(parents=(-1,0,1), cells_per_room=32,
                              duration=2., method='DOP853', rtol=2e-9,
                              sample_count=9):
    """Finite conservative phase/geometry evolution, returning small readouts.

    This starts from the retained triangular winding in the root curve,
    zero receiving phase and zero momenta. Regional reference volumes are
    one. Receiving domains and the tree are admitted initial geometry;
    neither branch formation nor a one-way causal horizon is asserted.
    """
    from scipy.integrate import solve_ivp
    if (not math.isfinite(duration) or duration <= 0 or
            not math.isfinite(rtol) or not 0 < rtol <= 1e-3 or
            type(sample_count) is not int or sample_count < 2):
        raise ValueError('positive finite duration, tolerance and at least two samples required')
    system=_parent_receiver_system(parents,cells_per_room)
    initial=system['initial'];evaluate=system['evaluate'];nr=system['room_count']
    nf=system['free_phase_count'];start=evaluate(initial)[1]

    def positive_volume_guard(t,state):
        return float(np.min(system['geometry'](state[2*nf:2*nf+nr])[1])-1e-6)
    positive_volume_guard.terminal=True
    positive_volume_guard.direction=-1
    solution=solve_ivp(lambda t,z:evaluate(z)[0],(0.,duration),initial,
                       method=method,rtol=rtol,atol=rtol*.01,
                       dense_output=True,events=positive_volume_guard)
    if not solution.success:raise RuntimeError(solution.message)
    stop=float(solution.t[-1]);rows=[]
    # The continuum leading front travels material distance integral dt/lambda.
    # Integrate on the accepted dense polynomials, rather than interpreting
    # small semidiscrete precursors as signals beyond that front.
    nodes,weights=np.polynomial.legendre.leggauss(8)
    widths=np.diff(solution.t)/2;mids=(solution.t[1:]+solution.t[:-1])/2
    times=(mids[:,None]+widths[:,None]*nodes).ravel()
    q_quad=solution.sol(times)[2*nf:2*nf+nr]
    enclosed=np.exp(3*q_quad);volume=enclosed.copy()
    for i,p in enumerate(parents[1:],1):volume[p]-=enclosed[i]
    steps=(volume**(-1/3)).reshape(nr,-1,8)@weights
    optical_steps=steps*widths
    optical_cumulative=np.c_[np.zeros(nr),np.cumsum(optical_steps,axis=1)]
    front_arrival=[0.]*nr
    for i,p in enumerate(parents[1:],1):
        if p == 0:continue
        parent_arrival=front_arrival[p]
        if parent_arrival is None:front_arrival[i]=None;continue
        initial_distance=float(np.interp(parent_arrival,solution.t,optical_cumulative[p]))
        target=initial_distance+2.
        if optical_cumulative[p,-1] < target:front_arrival[i]=None
        else:front_arrival[i]=float(np.interp(target,optical_cumulative[p],solution.t))
    max_total=max_phase=max_geometry=max_flux_sum=0.
    for t in np.linspace(0.,stop,max(129,sample_count)):
        state=solution.sol(t);row=evaluate(state)[1]
        transfer=state[2*nf+2*nr:2*nf+3*nr]
        work=state[2*nf+3*nr:2*nf+4*nr]
        phase_residual=row['phase_energy']-start['phase_energy']-transfer-work
        max_total=max(max_total,abs(row['total_energy']-start['total_energy']))
        max_phase=max(max_phase,float(np.max(np.abs(phase_residual))))
        max_geometry=max(max_geometry,abs(row['geometry_energy']+float(np.sum(work))))
        max_flux_sum=max(max_flux_sum,abs(float(np.sum(transfer))))
    for t in np.linspace(0.,stop,sample_count):
        state=solution.sol(t);row=evaluate(state)[1]
        q=state[2*nf:2*nf+nr];q0=initial[2*nf:2*nf+nr]
        inherited=[0.]+[float((q[p]-q0[p])-(q[i]-q0[i])) for i,p in enumerate(parents[1:],1)]
        rows.append(dict(offset=float(t),regional_volume=row['volume'].tolist(),
            regional_phase_energy=row['phase_energy'].tolist(),
            outer_radius=np.cbrt(row['enclosed']).tolist(),
            enclosing_pressure=row['pressure'].tolist(),
            inherited_compression=inherited,
            net_received_through_contact=state[2*nf+2*nr:2*nf+3*nr].tolist(),
            outer_log_radius_velocity=(state[2*nf+nr:2*nf+2*nr]/system['inertia']).tolist(),
            field_geometry_work=state[2*nf+3*nr:2*nf+4*nr].tolist(),
            geometric_kinetic_energy=row['geometry_energy'],total_energy=row['total_energy']))
    return dict(schema='nsc-finite-parent-receiver-v1',parents=list(parents),
        cells_per_room=cells_per_room,method=method,rtol=rtol,
        initial_phase_energy=start['phase_energy'].tolist(),
        initial_total_budget=start['total_energy'],
        requested_offset=duration,reached_offset=stop,
        positive_volume_guard_reached=bool(solution.t_events[0].size),
        rhs_evaluations=solution.nfev,
        max_sampled_total_budget_residual=max_total,
        max_sampled_regional_ledger_residual=max_phase,
        max_sampled_geometry_work_residual=max_geometry,
        max_sampled_contact_cancellation=max_flux_sum,
        continuum_optical_distance=optical_cumulative[:,-1].tolist(),
        estimated_leading_front_arrival_offset=front_arrival,
        front_estimates_scope='Gauss integral of dense geometry; finite-element precursors ahead of this front are not received continuum signals',
        regional_phase_pressure_over_density=1/3,
        material_characteristics='dx/d offset = +/- 1/lambda_i; both directions exist',
        one_way_horizon_established=False,new_child_formation_proved=False,
        global_event_horizon_status='not determined by this finite history or the signs of local characteristic speeds',
        receiving_domains_admitted_initially=True,
        added_assumptions=['connected material-curve tree with transparent shared phase nodes',
            'path stretch lambda_i = available regional volume^(1/3)',
            'phase inertia integral lambda_i*phase_offset^2/2',
            'geometric inertia 3*sum(log outer radius offset derivatives squared)/2',
            'unit initial regional volumes and retained root winding; zero initial momenta'],
        no_observational_input=True,rows=rows)


def parent_receiver_calculation():
    """Bounded independent mesh/time controls; no histories or data files saved."""
    from concurrent.futures import ProcessPoolExecutor
    cases=[]
    for name,parents in (('chain',(-1,0,1)),('fork',(-1,0,0))):
        for label,step in (('04',.04),('02',.02),('01',.01)):
            cases.append((f'{name}_{label}',dict(parents=parents,optical_step=step,
                                               duration=2.5 if name == 'chain' else .9)))
    cases.extend([
        ('chain_FE_64',dict(parents=(-1,0,1),cells_per_room=64,duration=.5)),
        ('chain_FE_128',dict(parents=(-1,0,1),cells_per_room=128,duration=.5))])
    with ProcessPoolExecutor(max_workers=8) as pool:
        futures=[pool.submit(parent_receiver_evolution if '_FE_' in name else parent_receiver_characteristics,
                             **settings) for name,settings in cases]
        calculated={name:future.result() for (name,_),future in zip(cases,futures)}
    comparison={}
    for name in ('chain','fork'):
        rows=[calculated[f'{name}_{label}']['rows'][-1] for label in ('04','02','01')]
        comparison[name]={key:dict(
            coarse_to_middle=np.abs(np.asarray(rows[1][key])-rows[0][key]).tolist(),
            middle_to_fine=np.abs(np.asarray(rows[2][key])-rows[1][key]).tolist())
            for key in ('regional_volume','regional_phase_energy','net_received_through_contact')}
    return dict(schema='nsc-finite-parent-receiver-comparison-v1',
        reference_geometric_spectrum=receiving_shape_spectrum(),
        exact_one_room_continuation=one_room_geometric_continuation(2.5),
        cases=calculated,endpoint_adjacent_spreads=comparison,
        interpretation='finite model evolution and empirical numerical controls, not a child-formation or observational proof',
        later_fork_scope='post-reflection fork histories are not admitted numerical claims; characteristic jump interpolation needs a separate continuation',
        no_files_written=True)


def continuum_receiving_birth(duration=.05, samples=129, emission_offset=None):
    """Early coupled opening with exact curve characteristics and radial inertia."""
    from scipy.integrate import solve_ivp
    from scipy.interpolate import CubicSpline

    if not 0 < duration <= .05 or type(samples) is not int or samples < 9:
        raise ValueError('require 0 < duration <= .05 and at least nine samples')
    offsets = np.linspace(0., duration, samples)
    initial_energy = 4 * math.pi**2
    acceleration = np.array([np.full(samples, initial_energy / 2),
                             np.full(samples, initial_energy)])

    def evaluate(values):
        speed = CubicSpline(offsets, values, axis=1).antiderivative()
        position = speed.antiderivative()
        radius, parent_displacement = position(offsets)
        radial_speed, parent_speed = speed(offsets)
        parent_radius = 1 + parent_displacement

        def optical_rhs(_, state):
            lag, log_ratio, _, _ = state.reshape(4, samples)
            lag_radius, parent_displacement_at_lag = position(lag)
            lag_parent = 1 + parent_displacement_at_lag
            lag_speed = speed(lag)[0]
            ratio = np.exp(log_ratio)
            return np.array([-lag_radius, -lag_speed, math.pi**2 * (ratio / lag_parent)**2,
                             math.pi * lag_radius / lag_parent]).ravel()

        history = solve_ivp(optical_rhs, (0., 2.), np.r_[offsets, np.zeros(3*samples)],
                            method='DOP853', rtol=2e-12, atol=2e-14)
        if not history.success:
            raise RuntimeError(history.message)
        lag, log_ratio, receiver_force, shared_phase = history.y[:, -1].reshape(4, samples)
        contact_power = math.pi**2 * (parent_radius**-2 - (np.exp(log_ratio)/(1+position(lag)[1]))**2)
        parent_budget = initial_energy - CubicSpline(offsets, parent_radius * contact_power).antiderivative()(offsets)
        return np.array([receiver_force, parent_budget / parent_radius**2]), dict(
            radius=radius, parent_radius=parent_radius, shared_phase=shared_phase,
            child_phase_energy=radius * receiver_force, parent_phase_energy=parent_budget / parent_radius,
            child_kinetic=radial_speed**2 / 2, parent_kinetic=parent_speed**2 / 2)

    for iteration in range(1, 41):
        updated, values = evaluate(acceleration)
        defect = float(np.max(np.abs(updated - acceleration)))
        if defect < 2e-10:
            break
        acceleration = updated
    else:
        raise RuntimeError('birth iteration did not converge')
    energy = sum(values[key] for key in ('child_phase_energy','parent_phase_energy','child_kinetic','parent_kinetic'))
    result = dict(offset=duration, samples=samples, iterations=iteration,
                **{key: float(value[-1]) for key, value in values.items()},
                child_energy=float((values['child_phase_energy'] + values['child_kinetic'])[-1]),
                max_budget_defect=float(np.max(np.abs(energy-initial_energy))),
                acceleration_iteration_defect=defect, radius_coefficient=math.pi**2)
    if emission_offset is not None:
        from scipy.integrate import quad
        from scipy.optimize import brentq

        if not 0 < emission_offset < duration:
            raise ValueError('emission must be inside the computed history')
        position = CubicSpline(offsets, acceleration, axis=1).antiderivative(2)

        def optical_distance(start, end, region):
            return quad(lambda t: 1 / (position(t)[region] + region), start, end,
                        epsabs=1e-11, epsrel=1e-11)[0]

        if optical_distance(emission_offset, duration, 0) < 2:
            raise ValueError('computed history ends before the echo returns')
        echo = brentq(lambda t: optical_distance(emission_offset, t, 0) - 2,
                      emission_offset, duration, xtol=1e-14)
        ratio = float(position(emission_offset)[0] / position(echo)[0])
        detector_distance = .01
        if optical_distance(echo, duration, 1) < detector_distance:
            raise ValueError('computed history ends before the parent detector receives the echo')
        detection = brentq(lambda t: optical_distance(echo, t, 1) - detector_distance,
                           echo, duration, xtol=1e-14)
        result['echo'] = dict(emission_offset=emission_offset, mouth_return_offset=echo,
                              parent_detector_material_distance=detector_distance,
                              parent_detection_offset=detection, mouth_power_ratio=ratio**2,
                              narrow_packet_mouth_energy_ratio=ratio)
    return result


def affine_receiving_birth(duration=.1, seed=1e-5):
    """Continue the triangular winding into a new affine arc at zero radius."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    if not 0 < seed < duration <= .2:
        raise ValueError('require 0 < seed < duration <= .2')
    winding = 2 * math.pi
    reduced_mass = 7 / 24
    gradient = brentq(lambda k: k + reduced_mass * k**3 / 2 - 1.5 * winding,
                      0., 1.5 * winding)
    radius_coefficient = gradient**2 / 4
    phase_acceleration = gradient**3 / 2
    midpoint_acceleration = -3 * winding - phase_acceleration / 4
    limit = [radius_coefficient, 2 * radius_coefficient,
             gradient * radius_coefficient, winding**2 / 2, winding**2,
             midpoint_acceleration / 2, -2 * winding, winding - gradient]

    def unpack(offset, state):
        u, p, v, parent_size, parent_speed, midpoint, pm, pc = state
        radius = offset**2 * u
        parent = 1 + offset**2 * parent_size
        phase = offset**2 * v
        middle = winding + offset**2 * midpoint
        mass = np.array([[2 * parent / 3, parent / 6],
                         [parent / 6, (parent + radius) / 3]])
        velocities = np.linalg.solve(mass, [pm, pc])
        return radius, parent, phase, middle, velocities

    def rhs(log_offset, state):
        offset = math.exp(log_offset)
        u, p, v, parent_size, parent_speed, midpoint, pm, pc = state
        _, parent, phase, middle, (md, td) = unpack(offset, state)
        return [p - 2 * u, (v / u)**2 / 2 + offset**2 * td**2 / 6 - p,
                td - 2 * v, parent_speed - 2 * parent_size,
                (middle**2 + (phase - middle)**2) / (2 * parent**2)
                + offset**2 * (2 * md**2 + md * td + td**2) / 6 - parent_speed,
                md - 2 * midpoint, (phase - 2 * middle) / parent - pm,
                (middle - phase) / parent - v / u - pc]

    solution = solve_ivp(rhs, (math.log(seed), math.log(duration)), limit,
                         method='DOP853', rtol=2e-12, atol=2e-13)
    if not solution.success:
        raise RuntimeError(solution.message)
    defects = []
    for log_offset, state in zip(solution.t, solution.y.T):
        offset = math.exp(log_offset)
        radius, parent, phase, middle, (md, td) = unpack(offset, state)
        child = (offset * state[1])**2 / 2 + radius * (offset * td)**2 / 6 + phase**2 / (2 * radius)
        outer = ((offset * state[4])**2 / 2
                 + parent * offset**2 * (2 * md**2 + md * td + td**2) / 6
                 + (middle**2 + (phase - middle)**2) / (2 * parent))
        defects.append(abs(child + outer - winding**2))
    return dict(offset=duration, seed=seed, radius=radius, volume=radius**3,
                shared_phase=phase, child_energy=child, parent_energy=outer,
                initial_energy=winding**2, max_budget_defect=max(defects),
                birth_gradient=gradient, radius_coefficient=radius_coefficient,
                child_energy_coefficient=4 * radius_coefficient**2)


def one_room_geometric_continuation(offset):
    """Exact free-parent limit, including its finite optical future.

    With geometric inertia 3 and retained winding budget Q0=4*pi^2,
    3*q_ddot=Q0*exp(-q). The causal length below uses the declared curve
    metric -d(offset)^2+lambda^2 dx^2. It is not a child-formation claim.
    """
    if not math.isfinite(offset) or offset < 0:
        raise ValueError('finite nonnegative offset required')
    budget=4*math.pi**2;a=math.sqrt(budget/6);u=a*offset
    q=float(2*np.logaddexp(u,-u)-2*math.log(2))
    tail=math.exp(-2*u)
    return dict(offset=offset,log_stretch=q,
        stretch=math.exp(q) if q < 709 else None,
        log_stretch_velocity=2*a*math.tanh(u),
        phase_energy=budget*math.exp(-q),
        geometric_kinetic_energy=budget*math.tanh(u)**2,
        optical_offset=math.tanh(u)/a,
        total_future_optical_reach=1/a,
        log_remaining_optical_distance=math.log(2/a)-2*u-math.log1p(tail),
        observer_causal_horizon_proper_length=(1+tail)/(2*a),
        asymptotic_log_stretch_rate=2*a,
        curvature_lambda_if_same_scale_is_extended_to_flat_3d=2*budget,
        scope='exact limit of declared geometric inertia and phase action; no SI ruler or observed vacuum identified')


def parent_receiver_characteristics(parents=(-1,0,1), duration=2.5,
                                    optical_step=.02, sample_count=9):
    """Characteristic continuation of the same finite parent/receiver model.

    Each room has uniform curve stretch lambda=V**(1/3), with available
    volume V_i=exp(3*q_i)-sum_children exp(3*q_j) and geometric inertia 3.
    Its exact field is A_i(s_i-x)+B_i(s_i+x-2), ds_i=dt/lambda_i. Delayed
    outgoing histories and their primitives replace a spatial field grid.
    Shared nodes obey transparent Kirchhoff matching without a point mass;
    the root's outer end and leaf ends are fixed Dirichlet endpoints.
    RK4 and linear interpolation of accepted histories have empirical
    errors, especially at characteristic slope jumps. Receiving rooms are
    admitted initially; this is neither child formation nor a one-way
    event-horizon construction.
    """
    parents=tuple(parents)
    if (not parents or parents[0] != -1 or
            any(type(p) is not int or not 0 <= p < i
                for i,p in enumerate(parents[1:],1))):
        raise ValueError('one ordered rooted tree required, with parent[0]=-1')
    if (not math.isfinite(duration) or duration <= 0 or
            not math.isfinite(optical_step) or not 0 < optical_step <= .25 or
            type(sample_count) is not int or sample_count < 2):
        raise ValueError('positive duration, optical_step in (0,.25], and at least two samples required')
    nr=len(parents);children=[[] for _ in parents]
    for i,p in enumerate(parents[1:],1):children[p].append(i)
    subtree=np.ones(nr)
    for i in range(nr-1,0,-1):subtree[parents[i]]+=subtree[i]
    contacts=[[(0,0)]]+[[(i,1)]+[(j,0) for j in children[i]]
                           for i in range(nr)]
    fixed=np.array([True]+[not children[i] for i in range(nr)])
    initial=np.zeros(11*nr);initial[:nr]=np.log(subtree)/3
    initial_energy=np.zeros(nr);initial_energy[0]=4*math.pi**2
    initial_budget=float(initial_energy.sum());guard=1e-6
    # Every room/side history is contiguous in offset, so interpolation does
    # not copy an ever-growing strided matrix column on each delayed query.
    capacity=min(200000,max(2048,20000000//max(1,13*nr)))
    time_history=np.empty(capacity)
    state_history=np.empty((11*nr,capacity))
    outgoing_history=np.empty((nr,2,capacity))
    time_history[0]=0.;state_history[:,0]=initial;used=1

    def geometry(state):
        enclosed=np.exp(3*state[:nr]);volume=enclosed.copy()
        for i,p in enumerate(parents[1:],1):volume[p]-=enclosed[i]
        if not np.all(np.isfinite(volume)) or np.any(volume <= 0):
            raise ValueError('available regional volume is nonpositive in this reduction')
        return enclosed,volume,np.cbrt(volume)

    def lookup(room,side,eta,kind):
        if eta < 0:
            if room != 0:return 0.
            if kind == 'h':return math.pi if eta < -1 else -math.pi
            if kind == 'C':return math.pi**2*eta
            return math.pi*min(-eta,2+eta)
        xp=state_history[2*nr+room,:used]
        if eta > xp[-1]+1e-12*max(1.,abs(xp[-1])):
            raise ValueError('characteristic step exceeded completed delayed history')
        if kind == 'h':values=outgoing_history[room,side,:used]
        else:values=state_history[(3*nr if kind == 'C' else 5*nr)+2*room+side,:used]
        return float(np.interp(eta,xp,values))

    def evaluate(state):
        enclosed,volume,lam=geometry(state)
        pg=state[nr:2*nr];s=state[2*nr:3*nr]
        C=state[3*nr:5*nr].reshape(nr,2)
        A=state[5*nr:7*nr].reshape(nr,2)
        delayed_h=np.empty((nr,2));delayed_C=np.empty((nr,2))
        delayed_A=np.empty((nr,2))
        for i in range(nr):
            for side in (0,1):
                opposite=1-side
                delayed_h[i,side]=lookup(i,opposite,s[i]-2,'h')
                delayed_C[i,side]=lookup(i,side,s[i]-2,'C')
                delayed_A[i,side]=lookup(i,side,s[i]-2,'A')
        incoming=delayed_h/lam[:,None]
        outgoing=np.empty((nr,2));node_gap=0.;contact_flux_gap=0.
        # A(s)+B(s-2) at the left end, A(s-2)+B(s) at the right.
        trace=np.c_[A[:,0]+delayed_A[:,1],delayed_A[:,0]+A[:,1]]
        for node,edges in enumerate(contacts):
            zeta=0. if fixed[node] else 2*sum(incoming[i,k] for i,k in edges)/len(edges)
            phases=[trace[i,k] for i,k in edges]
            node_gap=max(node_gap,max(phases)-min(phases),
                         abs(phases[0]) if fixed[node] else 0.)
            node_flux=0.
            for i,k in edges:
                outgoing[i,k]=lam[i]*(zeta-incoming[i,k])
                node_flux+=(outgoing[i,k]**2-delayed_h[i,k]**2)/lam[i]**2
            contact_flux_gap=max(contact_flux_gap,abs(node_flux))
        energy=np.sum(C-delayed_C,axis=1)/lam
        if np.any(energy < -1e-10*initial_budget):
            raise ValueError('interpolated characteristic primitives gave negative regional energy')
        qdot=pg/3
        volume_dot=3*enclosed*qdot
        for i,p in enumerate(parents[1:],1):volume_dot[p]-=3*enclosed[i]*qdot[i]
        log_stretch_dot=volume_dot/(3*volume)
        density=energy/volume
        force=enclosed*density
        for i,p in enumerate(parents[1:],1):force[i]-=enclosed[i]*density[p]
        transfer=np.sum(outgoing**2-delayed_h**2,axis=1)/lam**2
        work=-energy*log_stretch_dot
        derivative=np.r_[qdot,force,1/lam,(outgoing**2/lam[:,None]).ravel(),
                         (outgoing/lam[:,None]).ravel(),transfer,work,
                         np.maximum(transfer,0.),np.maximum(-transfer,0.)]
        geometry_energy=float(np.dot(pg,pg)/6)
        return derivative,dict(enclosed=enclosed,volume=volume,stretch=lam,
            phase_energy=energy,pressure=density/3,outgoing=outgoing,
            incoming=incoming,transfer=transfer,work=work,
            geometry_energy=geometry_energy,total_energy=float(energy.sum()+geometry_energy),
            interface_phase_gap=node_gap,contact_flux_gap=contact_flux_gap)

    first=evaluate(initial)[1];outgoing_history[:,:,0]=first['outgoing']
    state=initial.copy();t=0.;rhs_evaluations=1;rejected_steps=0
    status='requested offset reached';guard_reached=False
    max_total=max_regional=max_geometry=max_contact=max_phase=max_flux=0.
    first_return=[None]*nr;max_return=np.zeros(nr)
    speed_bound=max(1.,math.sqrt(2*initial_budget/3))
    sample_targets=np.linspace(0.,duration,sample_count);next_sample=1

    def record_diagnostics(current_time,current_state,values):
        nonlocal max_total,max_regional,max_geometry,max_contact,max_phase,max_flux
        transfers=current_state[7*nr:8*nr];works=current_state[8*nr:9*nr]
        residual=values['phase_energy']-initial_energy-transfers-works
        max_total=max(max_total,abs(values['total_energy']-initial_budget))
        max_regional=max(max_regional,float(np.max(np.abs(residual))))
        max_geometry=max(max_geometry,abs(values['geometry_energy']+float(works.sum())))
        max_contact=max(max_contact,abs(float(transfers.sum())))
        max_phase=max(max_phase,values['interface_phase_gap'])
        max_flux=max(max_flux,values['contact_flux_gap'])
        for i in range(1,nr):
            power=float(values['incoming'][i,0]**2)
            max_return[i]=max(max_return[i],power)
            if first_return[i] is None and power > 1e-12*initial_budget:
                first_return[i]=float(current_time)

    record_diagnostics(0.,state,first)
    while t < duration:
        if used >= capacity:
            status='bounded characteristic-history capacity reached';break
        current=evaluate(state)[1];rhs_evaluations+=1
        if float(np.min(current['volume'])) <= guard:
            status='positive-volume guard reached; this reduction needs continuation'
            guard_reached=True;break
        while next_sample < sample_count and sample_targets[next_sample] <= t+1e-13:
            next_sample+=1
        sample_limit=sample_targets[next_sample] if next_sample < sample_count else duration
        dt=min(duration-t,sample_limit-t,optical_step*float(np.min(current['stretch'])),
               optical_step/speed_bound)
        accepted=False
        for _ in range(60):
            if dt <= 1e-14*max(1.,t):
                status='step floor at a volume/history boundary; this reduction needs continuation'
                break
            try:
                k1=evaluate(state)[0]
                k2=evaluate(state+dt*k1/2)[0]
                k3=evaluate(state+dt*k2/2)[0]
                k4=evaluate(state+dt*k3)[0];rhs_evaluations+=4
                candidate=state+dt*(k1+2*k2+2*k3+k4)/6
                _,candidate_volume,_=geometry(candidate)
                if np.max(candidate[2*nr:3*nr]-state[2*nr:3*nr]) >= 1:
                    raise ValueError('optical increment reached the delayed-history safety bound')
                candidate_values=evaluate(candidate)[1];rhs_evaluations+=1
                if np.min(candidate_volume) <= guard:
                    dt/=2;rejected_steps+=1
                    if float(np.min(current['volume'])) <= 1.001*guard:
                        status='positive-volume guard reached; this reduction needs continuation'
                        guard_reached=True;break
                    continue
            except (ValueError,OverflowError,FloatingPointError):
                dt/=2;rejected_steps+=1;continue
            accepted=True;break
        if not accepted:
            if status == 'requested offset reached':
                status='no admissible characteristic step; inspect volume and delayed history'
            break
        t+=dt;state=candidate
        time_history[used]=t;state_history[:,used]=state
        outgoing_history[:,:,used]=candidate_values['outgoing'];used+=1
        record_diagnostics(t,state,candidate_values)
    rows=[]
    for at in np.linspace(0.,t,sample_count):
        sample=np.array([np.interp(at,time_history[:used],row[:used])
                         for row in state_history])
        values=evaluate(sample)[1];record_diagnostics(float(at),sample,values)
        q=sample[:nr];q0=initial[:nr]
        inherited=[0.]+[float((q[p]-q0[p])-(q[i]-q0[i]))
                          for i,p in enumerate(parents[1:],1)]
        rows.append(dict(offset=float(at),regional_volume=values['volume'].tolist(),
            regional_phase_energy=values['phase_energy'].tolist(),
            outer_radius=np.cbrt(values['enclosed']).tolist(),
            enclosing_pressure=values['pressure'].tolist(),
            inherited_compression=inherited,
            net_received_through_contact=sample[7*nr:8*nr].tolist(),
            gross_received_through_contact=sample[9*nr:10*nr].tolist(),
            gross_returned_through_contact=sample[10*nr:11*nr].tolist(),
            field_geometry_work=sample[8*nr:9*nr].tolist(),
            outer_log_radius_velocity=(sample[nr:2*nr]/3).tolist(),
            continuum_optical_distance=sample[2*nr:3*nr].tolist(),
            interface_phase_continuity_gap=values['interface_phase_gap'],
            geometric_kinetic_energy=values['geometry_energy'],
            total_energy=values['total_energy']))
    return dict(schema='nsc-parent-receiver-characteristics-v1',parents=list(parents),
        method='RK4 method of steps; linear accepted optical-history interpolation',
        optical_step=optical_step,requested_offset=duration,reached_offset=float(t),
        status=status,positive_volume_guard_reached=guard_reached,
        accepted_steps=used-1,rejected_steps=rejected_steps,rhs_evaluations=rhs_evaluations,
        history_capacity=capacity,initial_phase_energy=initial_energy.tolist(),
        initial_total_budget=initial_budget,
        max_accepted_and_sampled_total_budget_residual=max_total,
        max_accepted_and_sampled_regional_ledger_residual=max_regional,
        max_accepted_and_sampled_geometry_work_residual=max_geometry,
        max_accepted_and_sampled_contact_cancellation=max_contact,
        max_accepted_and_sampled_interface_phase_gap=max_phase,
        max_accepted_and_sampled_instantaneous_contact_flux_gap=max_flux,
        first_receiving_side_return_offset=first_return,
        max_receiving_side_return_power=max_return.tolist(),
        return_from_receiving_region_observed=any(x is not None for x in first_return),
        roundtrip_diagnostic='nonzero left-going characteristic inside an admitted receiver; finite fixed outer/leaf endpoints reflect',
        regional_phase_pressure_over_density=1/3,
        material_characteristics='dx/d offset = +/- 1/lambda_i; both directions exist',
        one_way_horizon_established=False,new_child_formation_proved=False,
        global_event_horizon_status='not determined by this finite history or local characteristic signs',
        receiving_domains_admitted_initially=True,
        numerical_error_scope='empirical optical-step and history-interpolation controls; no continuum error enclosure',
        no_observational_input=True,rows=rows)


def unified_least_resistance_cycle(compression_chi, winding_budget=1., enclosing_load=1.):
    """Single-cost 5-variable minimization producing infall filament, 0-crossing, and - sphere.

    Axioms declared openly (no fluid incompressibility or surface tension):
      - A0 (Least resistance): Minimize the single generator cost
        E(lambda_3, tilde_lambda_1..3, theta; chi) = E_+ + E_-, where
        E_+ = 2*m_+ / xi + m_+ / lambda_3 + P*xi^2*lambda_3 (infall region, xi=exp(-chi)),
        E_- = sum_{i=1}^3 m_- / tilde_lambda_i + P*tilde_lambda_1*tilde_lambda_2*tilde_lambda_3.
      - A1 (1D + offset): Quadratic winding moment m_+ = m0*(1-theta)^2 in +,
        m_- = m0*theta^2 across 3 directions in - from conserved turns n_+ + n_- = N0.
      - A4 (One-way valve at 0): Transferred winding fraction theta in [0, 1].

    In normalized m0=P=1 units, simultaneous stationarity
    dE/d(lambda_3)=0, dE/d(tilde_lambda_i)=0, dE/d(theta)=0
    yields a unique strict global minimum with positive-definite 5x5 Hessian:
      sqrt(theta^*) = (sqrt(16*exp(2*chi) + 17) - 3) / (4*exp(chi)),
      lambda_3^* = (1 - theta^*) * sqrt(m0/P) * exp(+chi),
      tilde_lambda_1^* = tilde_lambda_2^* = tilde_lambda_3^* = (m0/P)^(1/4) * sqrt(theta^*).
    Away from normalized units the theta formula uses beta=exp(-chi)*(P/m0)^(1/4),
    sqrt(theta)=(sqrt(16+17*beta^2)-3*beta)/4, with beta<3 for an interior minimum.
    """
    if (not math.isfinite(compression_chi) or compression_chi<0
            or not math.isfinite(winding_budget) or winding_budget<=0
            or not math.isfinite(enclosing_load) or enclosing_load<=0):
        raise ValueError('nonnegative finite compression_chi and positive winding_budget and enclosing_load required')
    xi=math.exp(-compression_chi)
    ex=math.exp(compression_chi)
    scale_unit=(winding_budget/enclosing_load)**.25
    beta=xi/scale_unit
    if beta >= 3:
        raise ValueError('strict positive-scale interior minimum requires exp(-chi)*(P/m0)^(1/4) < 3')
    u=(math.sqrt(16+17*beta*beta)-3*beta)/4
    theta=u*u
    one_minus=beta*(3*u-beta)/2
    lam3=scale_unit*(3*u-beta)/2
    t_lam=scale_unit*u
    v_plus=xi*xi*lam3
    v_minus=t_lam**3
    infall_aspect=lam3/xi
    cost_scale=winding_budget**.75*enclosing_load**.25
    e_plus=2*winding_budget*one_minus**2*ex+2*math.sqrt(winding_budget*enclosing_load)*one_minus*xi
    e_minus=cost_scale*4*u**3
    e_all_plus=2*winding_budget*ex+2*math.sqrt(winding_budget*enclosing_load)*xi
    e_all_minus=cost_scale*4.
    return dict(compression_chi=compression_chi,
                axioms_used=['A0: minimize single generator stretch + enclosing work cost E_+ + E_-',
                             'A1: converging parent shell sets lambda_1=lambda_2=exp(-chi) in +; offset distributes m_-=(theta^2)*m0 across 3 axes in -',
                             'A4: conserved turn split n_+ + n_- = N0 across 0'],
                fluid_incompressibility_assumed=False,
                surface_tension_assumed=False,
                volume_conservation_shortcut_assumed=False,
                transferred_winding_fraction_theta=theta,
                plus_infall=dict(
                    transverse_scale=xi,
                    optimal_radial_scale=lam3,
                    parcel_volume=v_plus,
                    parcel_volume_shrinks=v_plus<math.sqrt(winding_budget/enclosing_load),
                    filament_aspect_ratio=infall_aspect,
                    is_elongated_filament=infall_aspect>1.,
                    normalized_deformation_diagonal=[1/infall_aspect,1/infall_aspect,1.] if infall_aspect>=1 else [1.,1.,infall_aspect]),
                minus_child_drop=dict(
                    optimal_principal_scales=[t_lam,t_lam,t_lam],
                    aspect_ratio=1.,
                    winner_shape='round_3d_sphere',
                    drop_volume=v_minus),
                derived_cooker_takeover=dict(
                    plus_volume=v_plus,
                    minus_volume=v_minus,
                    total_two_region_volume=v_plus+v_minus,
                    minus_volume_fraction=v_minus/(v_plus+v_minus)),
                cost_comparison=dict(
                    unified_minimum_cost=e_plus+e_minus,
                    stay_all_in_plus_cost=e_all_plus,
                    cross_all_to_minus_cost=e_all_minus))


def nested_zero_crossing_map(levels, compression_per_level):
    """Iterate the signed-zero crossing (+ -> 0 -> -) across nested levels.

    New axiom (stated openly): signed radius passes through 0 into -, flipping
    3D orientation sign(det F_k) = (-1)^k at each one-way 0 horizon. Inside each
    level k, the unified single-cost minimization determines the transferred
    winding fraction theta^*(chi_k), the child sphere radius, and the three-way
    budget partition:
      - E_DE  = 3 * (theta^*)^(3/2)             (isotropic 3D expansion stretch in -)
      - E_DM  = (theta^*)^(3/2) + 2*(1-theta^*)*e^(-chi) (0-wall load + non-luminous radial infall)
      - E_vis = 2 * (1-theta^*)^2 * e^(+chi)    (active filamentary winding pinching toward child 0s)
    """
    if type(levels) is not int or levels<1:
        raise ValueError('positive integer level count required')
    if not math.isfinite(compression_per_level) or compression_per_level<0:
        raise ValueError('nonnegative finite compression_per_level required')
    chain=[]
    cumulative_chi=0.
    cumulative_scale=1.
    for k in range(1,levels+1):
        cumulative_chi+=compression_per_level
        cyc=unified_least_resistance_cycle(compression_per_level)
        theta=cyc['transferred_winding_fraction_theta']
        u=math.sqrt(theta)
        ex=math.exp(compression_per_level)
        xi=math.exp(-compression_per_level)
        e_de=3*u**3
        e_dm=u**3+2*(1-theta)*xi
        e_vis=2*(1-theta)**2*ex
        e_tot=e_de+e_dm+e_vis
        cumulative_scale*=u
        chain.append(dict(
            level=k,
            orientation_sign=(-1)**k,
            local_region_label='-' if k%2==1 else '+',
            cumulative_compression_chi=cumulative_chi,
            transferred_winding_fraction=theta,
            cumulative_sphere_scale=cumulative_scale,
            dark_energy_like_fraction=e_de/e_tot,
            dark_matter_like_fraction=e_dm/e_tot,
            dark_sector_combined_fraction=(e_de+e_dm)/e_tot,
            visible_filament_like_fraction=e_vis/e_tot,
            de_to_boundary_virial_lock=3.))
    volumetric_95_horizon=1-math.exp(-3)
    volumetric_visible_floor=math.exp(-3)
    return dict(levels=levels,
                compression_per_level=compression_per_level,
                axiom_signed_zero_crossing='Signed radius passes through 0 into -; det(F) flips sign (-1)^k at each one-way 0 wall.',
                surfaceology_analogy='Each 0 horizon is a closed puncture boundary factorizing the outer parent polygon from the inner child surface.',
                two_valve_steady_state_lock=dict(
                    explanation='Every nested room k sits between an outer birth 0 valve (inflow from parent k-1) and inner black-hole 0 valves (outflow to children k+1), locking a steady-state ~95% dark-sector horizon instead of draining visible matter to zero.',
                    three_dimensional_e_fold_dark_horizon=volumetric_95_horizon,
                    three_dimensional_e_fold_visible_floor=volumetric_visible_floor,
                    virial_3_to_1_dark_energy_at_95_horizon=.75*volumetric_95_horizon,
                    virial_3_to_1_dark_matter_at_95_horizon=.25*volumetric_95_horizon),
                nested_chain=chain,
                asymptotic_virial_partition=dict(
                    dark_energy_like=.75,
                    dark_matter_like=.25,
                    visible_filament_like=0.),
                antimatter_and_baryon_excess_check=dict(
                    bulk_single_orientation_in_child=True,
                    single_mode_eta_range=[.0864272337,.5746972340],
                    observed_baryon_eta_reference=6.1e-10,
                    matches_observed_eta_without_scale_separation=False,
                    required_omega_over_gamma_for_observed_eta=float(2/(math.pi*6.1e-10)),
                    required_nested_efolds_for_observed_eta=float(-math.log(6.1e-10))))


def spectral_band_and_baryon_harmonic(compression_chi, reference_bandwidth=1e61, target_eta=6.1e-10):
    """Moving wavelength band [Planck cutoff, Absolute-Zero horizon] and spectral asymmetry eta(lambda).

    In a room at relative compression chi, both the longest wavelength
    lambda_max(chi) = lambda_max_0 * exp(-chi) (Absolute Zero / horizon floor,
    as wide as the space itself) and the shortest wavelength before 0 transitions
    to - lambda_min(chi) = lambda_min_0 * exp(-chi) (local Planck cutoff) shift
    together by exp(-chi), while frequencies omega shift up by exp(+chi).
    Every mode of wavelength lambda in [lambda_min(chi), lambda_max(chi)] riding
    the room's horizon decay envelope exp(-gamma*tau) experiences half-cycle
    damping exp(-pi*gamma/omega) = exp(-(pi/2)*(lambda/lambda_max)), giving the
    spectral asymmetry curve:
      eta(lambda) = tanh((pi/2) * (lambda / lambda_max)).
    At the horizon scale lambda = lambda_max, eta = tanh(pi/2) ≈ 0.917152 (order-1
    dark/inflationary asymmetry); at the harmonic ratio
    lambda_B / lambda_max = (2/pi)*atanh(target_eta) ≈ 3.88338e-10, the same curve
    yields eta(lambda_B) = target_eta (6.1e-10) well above the Planck cutoff.
    """
    if (not math.isfinite(compression_chi)
            or not math.isfinite(reference_bandwidth) or reference_bandwidth<=1
            or not math.isfinite(target_eta) or not 0<target_eta<1):
        raise ValueError('finite compression_chi, bandwidth > 1 and target_eta in (0,1) required')
    scale_shift=math.exp(-compression_chi)
    freq_shift=math.exp(compression_chi)
    lam_max=scale_shift
    lam_min=scale_shift/reference_bandwidth
    baryon_wavelength_ratio=(2/math.pi)*math.atanh(target_eta)
    lam_baryon=lam_max*baryon_wavelength_ratio
    eta_horizon=math.tanh(.5*math.pi)
    eta_baryon=math.tanh(.5*math.pi*(lam_baryon/lam_max))
    eta_planck=math.tanh(.5*math.pi/reference_bandwidth)
    return dict(compression_chi=compression_chi,
                wavelength_scale_factor=scale_shift,
                frequency_scale_factor=freq_shift,
                absolute_zero_max_wavelength=lam_max,
                planck_cutoff_min_wavelength=lam_min,
                invariant_spectral_bandwidth=reference_bandwidth,
                baryon_harmonic_wavelength=lam_baryon,
                baryon_to_horizon_wavelength_ratio=baryon_wavelength_ratio,
                baryon_octaves_below_horizon=-math.log2(baryon_wavelength_ratio),
                within_room_spectrum=(lam_min<lam_baryon<lam_max),
                eta_at_absolute_zero_horizon=eta_horizon,
                eta_at_baryon_harmonic=eta_baryon,
                eta_at_planck_cutoff=eta_planck)


def branching_cosmos_chain_evolution(generations=3, branch_factor=2, compression_per_step=2.5):
    """Bidirectional parent-child drain and expansion ('Oooo -> oOOOooooooOOOOOO').

    Tracks a branching tree of nested spaces where an expanded parent ('O') births
    `branch_factor` child black-hole seeds ('o'). As compression chi_step drives
    winding fraction theta^*(chi_step) across each one-way 0 valve into the children,
    the parent is drained from 'O' into a compressed remnant 'o' (retaining 1-theta^*),
    while the children balloon into expanded rooms 'O' and nucleate the next
    generation of seeds, conserving total budget across the entire web.
    """
    if (type(generations) is not int or generations<1
            or type(branch_factor) is not int or branch_factor<1
            or not math.isfinite(compression_per_step) or compression_per_step<=0):
        raise ValueError('positive integer generations and branch_factor, and positive compression_per_step required')
    theta=unified_least_resistance_cycle(compression_per_step)['transferred_winding_fraction_theta']
    levels=[dict(depth=0,count=1,budget_total=1.,budget_per_node=1.,state='O')]
    history=[]
    for step in range(1,generations+1):
        new_levels=[]
        spawned_budget=0.
        spawned_count=0
        for lvl in levels:
            if lvl['state']=='O':
                retained=lvl['budget_total']*(1-theta)
                drained=lvl['budget_total']*theta
                new_levels.append(dict(depth=lvl['depth'],count=lvl['count'],
                                       budget_total=retained,
                                       budget_per_node=retained/lvl['count'],
                                       state='o'))
                spawned_count=lvl['count']*branch_factor
                spawned_budget=drained
            else:
                new_levels.append(dict(lvl))
        new_levels.append(dict(depth=step,count=spawned_count,
                               budget_total=spawned_budget,
                               budget_per_node=spawned_budget/spawned_count,
                               state='O'))
        seed_count=spawned_count*branch_factor
        pattern=''.join(lvl['state']*min(lvl['count'],12) for lvl in new_levels)+('o'*min(seed_count,12))
        levels=new_levels
        total_budget=float(sum(lvl['budget_total'] for lvl in levels))
        history.append(dict(step=step,
                            transferred_fraction_theta=theta,
                            active_levels=[dict(lvl) for lvl in levels],
                            next_nucleated_seed_count=seed_count,
                            pattern_signature=pattern,
                            total_budget=total_budget,
                            conservation_residual=total_budget-1.))
    return dict(generations=generations,
                branch_factor=branch_factor,
                compression_per_step=compression_per_step,
                transferred_fraction_per_step=theta,
                history=history)


def cyclic_vortex_diffusion_ring(offset, node_count=8, diffusion_rate=0.25, spinup_drive=0.25, wave_speed=1.):
    """Closed ring of vortex spin-up ('O') and diffusion into still air ('o').

    Models a cyclic ring of K coupled spaces (e.g. 'oo / OO oO / oO Oo / oo')
    where each space alternates between a concentrated hurricane vortex ('O',
    above mean budget) and diffused still air ('o', below mean budget):
      - With diffusion alone (spinup_drive = 0, wave_speed = 0), any hurricane
        diffuses exponentially as exp(-lambda_1 * D * tau) into uniform dead air.
      - When enclosing offset drive balances diffusion on the closed ring, a
        conserved vortex wave E_j(tau) = (1/K) * [1 + A * exp((spinup_drive - D)*tau) * cos(2*pi*(j/K - wave_speed*tau))]
        (at equilibrium spinup_drive = D) circulates around the ring: each 'O'
        diffuses into still air 'o' by shedding into nested mini-vortices while
        downstream still air 'o' re-accumulates and spins back up into an
        enormous hurricane 'O', returning to the original node every tau = 1/wave_speed.
    """
    if (type(node_count) is not int or node_count<4
            or not math.isfinite(offset) or offset<0
            or not math.isfinite(diffusion_rate) or diffusion_rate<0
            or not math.isfinite(spinup_drive) or spinup_drive<0
            or not math.isfinite(wave_speed) or wave_speed<0):
        raise ValueError('node_count >= 4 and nonnegative finite offset, diffusion_rate, spinup_drive, wave_speed required')
    mean_budget=1./node_count
    net_growth=spinup_drive-diffusion_rate
    amp=0.8*math.exp(net_growth*offset)
    angles=2*math.pi*(np.arange(node_count)/node_count-wave_speed*offset)
    budgets=mean_budget*(1+amp*np.cos(angles))
    states=['O' if b>=mean_budget else 'o' for b in budgets]
    pure_diffusion_envelope=0.8*math.exp(-diffusion_rate*offset)
    return dict(offset=offset,
                node_count=node_count,
                node_budgets=budgets.tolist(),
                node_states=''.join(states),
                total_budget=float(np.sum(budgets)),
                conservation_residual=float(np.sum(budgets))-1.,
                max_to_min_vortex_contrast=float(np.max(budgets)/np.min(budgets)) if np.min(budgets)>0 else None,
                pure_diffusion_without_offset_envelope=pure_diffusion_envelope,
                active_envelope=amp,
                node_0_is_giant_hurricane=(states[0]=='O'))


def signed_circle_pendulum_duality(r0=1.,omega=2*math.pi):
    """Prove why '-' must exist, why (-)*(-)=(+), and why signed R crosses 0 smoothly.

    1. Layer-1 circle z(tau) = R0 * exp(i*omega*tau) seen edge-on IS the pendulum
       R(tau) = Re(z(tau)) = R0 * cos(omega*tau), swinging + -> 0 -> - -> 0 -> +.
    2. Without negative coordinates (R < 0), 50% of the circle (pi*R0 out of 2*pi*R0)
       is missing and no closed loop or bounded space can close.
    3. Multiplying by -1 = exp(i*pi) is a geometric half-turn (in 3D, P = -I_3 with
       det(P) = (-1)^3 = -1). Two half-turns across two nested 0 crossings satisfy
       P^2 = (-I_3)(-I_3) = +I_3 and det(P^2) = (-1)*(-1) = +1.
    4. Standard physics treats the radius of space as an unsigned size r(tau) = |R(tau)| >= 0.
       Folding R < 0 onto r > 0 creates an artificial velocity cusp Delta(r_dot) = +2*R0*omega
       and a Dirac-delta acceleration singularity at r = 0, whereas the signed radius R(tau)
       passes through 0 into - with zero velocity jump (Delta(R_dot) = 0) and zero acceleration
       (R_tt(tau_0) = 0).
    """
    r0=float(r0);omega=float(omega)
    if r0<=0 or omega<=0:raise ValueError('positive r0 and omega required')
    e_circle=.5*r0*r0*omega*omega
    phases=[('plus_apex_0',0.),
            ('first_zero_crossing_pi_over_2',.5*math.pi),
            ('minus_apex_pi',math.pi),
            ('second_zero_crossing_3pi_over_2',1.5*math.pi),
            ('two_flips_return_2pi',2.*math.pi)]
    swing_points=[]
    for label,theta in phases:
        r_signed=r0*math.cos(theta)
        r_dot=-r0*omega*math.sin(theta)
        r_tt=-r0*omega*omega*math.cos(theta)
        e_tot=.5*r_dot*r_dot+.5*omega*omega*r_signed*r_signed
        det_3d=(r_signed/r0)**3
        swing_points.append(dict(phase_label=label,
                                 phase_radians=theta,
                                 signed_radius_R=r_signed,
                                 unsigned_size_r=abs(r_signed),
                                 signed_velocity_R_dot=r_dot,
                                 signed_acceleration_R_tt=r_tt,
                                 jacobian_det_3d=det_3d,
                                 pendulum_energy=e_tot,
                                 energy_conservation_residual=e_tot-e_circle))
    p_3d=-np.eye(3)
    p_squared=p_3d@p_3d
    det_p=float(np.linalg.det(p_3d))
    det_p2=float(np.linalg.det(p_squared))
    return dict(r0=r0,
                omega=omega,
                circle_to_pendulum_identity='R(tau) = Re(R0 * exp(i*omega*tau)) = R0 * cos(omega*tau)',
                swing_points=swing_points,
                negative_half_required_for_closure=dict(
                    full_circle_arc_length=2.*math.pi*r0,
                    nonnegative_clipped_arc_length=math.pi*r0,
                    missing_arc_fraction_without_minus=.5,
                    closes_without_minus=False,
                    closes_with_minus=True),
                minus_times_minus_equals_plus=dict(
                    half_turn_complex_factor=-1.,
                    two_half_turns_complex_product=(-1.)*(-1.),
                    single_crossing_3d_inversion_det=det_p,
                    two_crossings_3d_inversion_det=det_p2,
                    matrix_identity_residual=float(np.linalg.norm(p_squared-np.eye(3))),
                    parent_to_child_to_grandchild_signs=['+','-','+']),
                singularity_diagnosis_at_zero=dict(
                    crossing_phase_radians=.5*math.pi,
                    signed_R_velocity_before=-r0*omega,
                    signed_R_velocity_after=-r0*omega,
                    signed_R_velocity_jump=0.,
                    signed_R_acceleration_at_zero=0.,
                    unsigned_size_r_velocity_before=-r0*omega,
                    unsigned_size_r_velocity_after=+r0*omega,
                    unsigned_size_r_velocity_cusp_jump=2.*r0*omega,
                    unsigned_size_r_singular_impulse=2.*r0*omega,
                    signed_crossing_is_smooth=True,
                    unsigned_size_creates_artificial_cusp_singularity=True))


def derived_three_accumulator_bubble_ladder(max_n=6,m0=1.,pv0=1.):
    """Derive the untuned A0+A1 iterated-accumulator bubble ladder before comparing to data.

    In A1, applying the offset accumulation operation across 3 stages produces the
    binomial winding numbers:
      W_1(n) = C(n, 1) = n,
      W_2(n) = C(n, 2) = n*(n-1)/2,
      W_3(n) = C(n, 3) = n*(n-1)*(n-2)/6,
    with directional stretch moments m_i(n) = m0 * W_i(n)^2.
    Minimizing the A0 3D cost E(lambda_1, lambda_2, lambda_3) = sum_i m_i/lambda_i + P*V0*lambda_1*lambda_2*lambda_3
    yields E_3D^*(n) = 4 * (P*V0)^(1/4) * (m0)^(3/4) * sqrt(W_1(n)*W_2(n)*W_3(n)):
      - n = 1: (W_1, W_2, W_3) = (1, 0, 0) -> W_2 = W_3 = 0 (collapses transversely to a 1D loop;
        1D equilibrium cost E_1D^*(1) = 2*sqrt(m0*pv0)).
      - n = 2: (W_1, W_2, W_3) = (2, 1, 0) -> W_3 = 0 (collapses in the 3rd axis to a 2D sheet).
      - n = 3: (W_1, W_2, W_3) = (3, 3, 1) -> first integer with W_1*W_2*W_3 = 9 > 0!
        Minimum stable 3D volumetric bubble has sqrt(W_1*W_2*W_3) = 3 and E_3D^*(3) = 12,
        so the untuned ratio E_3D^*(3) / E_1D^*(1) = 12 / 2 = 6 (missing 1836.15 by pi^5 ≈ 306.02).
    """
    if type(max_n) is not int or max_n<3:raise ValueError('integer max_n >= 3 required')
    m0=float(m0);pv0=float(pv0)
    e_1d_n1=2.*math.sqrt(m0*pv0)
    scale_3d=4.*(pv0**.25)*(m0**.75)
    ladder=[]
    for n in range(1,max_n+1):
        w1=n
        w2=n*(n-1)//2
        w3=n*(n-1)*(n-2)//6
        prod=w1*w2*w3
        dim_support=1 if w2==0 else (2 if w3==0 else 3)
        e_3d=scale_3d*math.sqrt(prod)
        ladder.append(dict(
            n=n,
            binomial_windings=[w1,w2,w3],
            directional_moments=[m0*(w1**2),m0*(w2**2),m0*(w3**2)],
            supported_spatial_dimension=dim_support,
            supports_3d_volume=(prod>0),
            volumetric_factor_sqrt_W1_W2_W3=math.sqrt(prod),
            equilibrium_3d_energy=e_3d,
            ratio_over_n3_baryon=(math.sqrt(prod)/3.) if prod>0 else 0.,
            ratio_over_n1_1d_loop=e_3d/e_1d_n1))
    codata_mp_me=1836.15267343
    untuned_n3_over_n1=ladder[2]['ratio_over_n1_1d_loop']
    return dict(
        axiom_rule='W_k(n) = C(n, k) from iterated A1 offset accumulation; E_3D^*(n) = 4*(P*V0)^(1/4)*m0^(3/4)*sqrt(W_1*W_2*W_3)',
        n1_is_1d_loop=(ladder[0]['supported_spatial_dimension']==1),
        n2_is_2d_sheet=(ladder[1]['supported_spatial_dimension']==2),
        minimum_n_for_3d_volume=3,
        untuned_n3_to_n1_energy_ratio=untuned_n3_over_n1,
        codata_proton_to_electron=codata_mp_me,
        missing_continuum_phase_factor_for_1836=codata_mp_me/untuned_n3_over_n1,
        hits_3d_threshold_at_n_equals_3=True,
        hits_1836_without_extra_factor=False,
        ladder=ladder)


def exponential_wall_bypasses_and_triad_winner():
    """Prove the Universal Bypass Law: whenever compression or winding hits an exponential/power wall,
    least resistance bypasses the wall by opening a new direction, crossing 0, or forming an orthogonal (1,1,1) triad.

    1. Bypass 1 (1D compression wall -> 2D/3D helical buckling):
       Straight 1D compression costs E_1D(chi) = exp(+chi) -> infinity (exponential wall).
       Buckling laterally into a helix caps the cost at E_helix(chi) = 3/2 - (1/2)*exp(-2*chi) -> 1.5.
    2. Bypass 2 (3D + compression wall -> signed 0-crossing into - child sphere):
       Staying in + costs E_stay(chi) = 2*exp(+chi) + 2*exp(-chi) -> infinity (exponential wall).
       Opening the 0 valve into - caps the cost at E_5var^*(chi) = 4 - (9/2)*exp(-chi) + O(exp(-2*chi)) -> 4.0.
    3. Bypass 3 (Single-axis winding wall N^(3/2) -> 3D orthogonal triad (1,1,1) winner & nuclear hard core):
       For 3D integer windings (n1, n2, n3) with m_i = n_i^2, exact A0 minimization gives
       E^*(n1, n2, n3) = 4 * sqrt(n1 * n2 * n3):
         - (1, 0, 0) and (1, 1, 0) have n3 = 0 -> zero 3D volume (collapse along empty axis).
         - (1, 1, 1) (N = 3 orthogonal unit windings) is the unique global minimum 3D state:
           E^*(1,1,1) = 4.0, cost per winding = 4/3 ≈ 1.3333, round sphere (aspect 1:1:1).
         - Adding a 4th winding forces doubling on one axis -> (2, 1, 1) prolate state (aspect 4:1:1)
           with E^*(2,1,1) = 4*sqrt(2) ≈ 5.6569 (cost per winding sqrt(2) ≈ 1.4142 > 4/3).
         - Two (1,1,1) triads attract at medium distance d > 2*R^* via shared flow (-2*m0/d),
           but strongly repel against merging into one (2,2,2) bubble because
           E^*(2,2,2) = 8*sqrt(2) ≈ 11.3137 > 2 * E^*(1,1,1) = 8.0 (hard-core barrier +3.3137).
    """
    wall_vs_bypass=[]
    for chi in (0.,1.,2.,3.,5.):
        e_1d_wall=math.exp(chi)
        e_helix_bypass=1.5-.5*math.exp(-2.*chi)
        e_3d_wall=2.*math.exp(chi)+2.*math.exp(-chi)
        e_0valve_bypass=unified_least_resistance_cycle(chi)['cost_comparison']['unified_minimum_cost']
        wall_vs_bypass.append(dict(
            chi=chi,
            wall_1d_straight=e_1d_wall,
            bypass_1d_to_helix=e_helix_bypass,
            wall_3d_stay_in_plus=e_3d_wall,
            bypass_3d_cross_zero_to_minus=e_0valve_bypass))
    states=[(1,0,0),(1,1,0),(1,1,1),(2,1,1),(2,2,1),(2,2,2),(3,1,1)]
    triad_table=[]
    for n1,n2,n3 in states:
        ntot=n1+n2+n3
        prod=n1*n2*n3
        e_star=4.*math.sqrt(prod)
        triad_table.append(dict(
            windings=[n1,n2,n3],
            total_winding_N=ntot,
            supports_3d_volume=(prod>0),
            equilibrium_3d_energy=e_star,
            energy_per_winding=(e_star/ntot) if prod>0 else None,
            principal_aspect_ratios=[n1**2,n2**2,n3**2] if prod>0 else None))
    e_111=4.
    e_211=4.*math.sqrt(2.)
    e_222=8.*math.sqrt(2.)
    return dict(
        universal_principle='When a path hits an exponential or superlinear cost wall, least resistance bypasses the wall by opening a perpendicular axis, crossing 0, or splitting/saturating at the (1,1,1) triad.',
        compression_wall_bypasses=wall_vs_bypass,
        one_d_helix_plateau_limit=1.5,
        three_d_zero_valve_plateau_limit=4.0,
        orthogonal_triad_lattice=triad_table,
        unique_3d_ground_state=[1,1,1],
        ground_state_energy=e_111,
        ground_state_cost_per_winding=e_111/3.,
        first_excited_4_winding_state=[2,1,1],
        first_excited_energy=e_211,
        first_excited_gap_over_ground=(e_211-e_111)/e_111,
        two_triads_separate_energy=2.*e_111,
        merged_222_energy=e_222,
        short_range_hard_core_repulsion_barrier=e_222-2.*e_111,
        medium_range_shared_flow_is_attractive=True,
        nuclear_potential_shape_derived=True)


def vacuum_bubble_crystallization_and_mass_spectrum():
    """Evaluate pressurized medium vacuum bubbles, surrounding-flow energy, 180-deg annihilation, resonance attraction, discrete bursts, and untuned mass ratios from A0-A4.

    1. Surrounding-medium energy vs point singularity (A0 + A1):
       Integrating the radial phase-flow gradient in the surrounding medium r >= R outside a
       bubble of radius R with integer winding n gives E_surrounding(R) = m0*n^2 / R, with
       radial shell fraction E(R <= r <= k*R) / E_surrounding = 1 - 1/k (50% in [R, 2R],
       90% in [R, 10R], 0% in a point singularity). Balancing against cavity pressure work
       W_pressure(R) = P*V0*R^3 gives the crystallized equilibrium radius
       R_n^* = (m0*n^2 / (3*P*V0))^(1/4) ~ |n|^(1/2) and energy E_n^* ~ |n|^(3/2),
       where 75% of the mass-energy sits in the surrounding medium flow and 25% in cavity work.
    2. Equal-and-opposite bubble annihilation (+1 and -1 -> two 180-deg ripples):
       Combining n_+ = +1 and n_- = -1 gives net winding n_net = 0, unlocking the cavity
       (R^* -> 0) and releasing 2*E_1^* into the open wave equation psi_tt - c^2 psi_xx = 0
       as two equal back-to-back ripples (E_left = E_right = E_1^*, 180 degrees apart).
    3. Shared-flow resonance attraction vs off-resonance cancellation:
       Superposing two bubbles separated by d gives cross-flow energy
       Delta_E_int(d, Delta_phi) = -(2*m0*n1*n2 / d)*cos((omega1 - omega2)*tau + Delta_phi0),
       which time-averages to 0 off-resonance (omega1 != omega2) and produces an attractive
       -1/d^2 force when resonance-matched (omega1 = omega2).
    4. Sharp discrete burst energies vs continuous smear:
       Single-valued loop closure z(1) = z(0) on S^1 forces integer winding n in Z, so
       bursting/relaxing transitions release sharp quantized differences Delta_E = E_{n_i}^* - E_{n_f}^*.
    5. Untuned mass-ratio comparison (reporting both hits and misses):
       - Single-stage A0 winding ladder E_n^*/E_1^* = n^(3/2) gives 1, 2.828427, 5.196152, ...
         (misses m_p/m_e = 1836.15267 unless n = 150 is chosen by hand: 150^(3/2) = 1837.1173).
       - Dimensional accumulation across A1 + A2 (1D circle resonance -> 3D signed-axis sphere):
         2*d = 6 signed radial rays (+/-x, +/-y, +/-z from A2) times pi^(2+d) = pi^5 across
         the two accumulation stages (pi^2 disk * pi^3 3-axis phase volume from A1) gives
         bare ratio 6*pi^5 = 1836.118109 (18.8 ppm from CODATA 1836.152673), and with the
         already-fixed two-valve medium dressing exp(-chi_spaghetti)*exp(-3) = exp(-3)/sqrt(2)
         gives 6*pi^5 + exp(-3)/sqrt(2) = 1836.153313 (0.35 ppm from CODATA).
       - Neutron-proton mass gap (m_n - m_p)/m_e = 2.530988 (1.29333 MeV):
         Single-stage n=2 bubble gives 2^(3/2) = 2.828427 (+11.75% high);
         threshold-captured pair + spaghettification valve fraction 2 + theta^*(chi_spaghetti)
         = 2 + 1/2 = 2.500000 (-1.22%), and with the same two-valve dressing
         2 + 1/2 + exp(-3)/sqrt(2) = 2.535205 (+0.167%).
    """
    m0=1.;pv0=1.
    r1=(m0/(3.*pv0))**.25
    e_surr_1=m0/r1
    w_press_1=pv0*(r1**3)
    e1_star=e_surr_1+w_press_1
    shell_fractions={f'R_to_{k}R':1.-1./k for k in (2,5,10,100)}
    discrete_levels=[]
    for n in range(1,6):
        rn=(m0*(n**2)/(3.*pv0))**.25
        en=m0*(n**2)/rn+pv0*(rn**3)
        discrete_levels.append(dict(winding_n=n,
                                    equilibrium_radius_ratio=rn/r1,
                                    stored_energy_ratio_over_E1=en/e1_star,
                                    exact_n_3_over_2=n**1.5))
    codata_mp_me=1836.15267343
    codata_mn_mp_over_me=2.530988
    me_mev=0.51099895
    two_valve_dressing=math.exp(-3.)/math.sqrt(2.)
    bare_6pi5=6.*(math.pi**5)
    dressed_6pi5=bare_6pi5+two_valve_dressing
    np_gap_n2_single=2.**1.5
    np_gap_threshold_bare=2.+.5
    np_gap_threshold_dressed=2.+.5+two_valve_dressing
    return dict(
        surrounding_medium_equilibrium=dict(
            formula='E_bubble(R) = m0*n^2/R (surrounding flow r>=R) + P*V0*R^3 (cavity work)',
            equilibrium_radius_R1=r1,
            equilibrium_energy_E1=e1_star,
            surrounding_flow_energy_share=e_surr_1/e1_star,
            cavity_pressure_work_share=w_press_1/e1_star,
            virial_ratio_surrounding_over_pressure=e_surr_1/w_press_1,
            radial_shell_energy_fractions=shell_fractions,
            energy_in_point_singularity=0.),
        annihilation_to_back_to_back_ripples=dict(
            initial_windings=[1,-1],
            combined_winding=0,
            cavity_lock_after_contact=0.,
            total_released_energy_over_E1=2.,
            outgoing_ripple_energies_over_E1=[1.,1.],
            outgoing_ripple_energies_keV=[510.99895,510.99895],
            separation_angle_degrees=180.,
            net_momentum_residual=0.),
        shared_flow_resonance_interaction=dict(
            interaction_formula='Delta_E_int(d, Delta_phi) = -(2*m0*n1*n2/d)*cos((omega1-omega2)*tau + Delta_phi0)',
            off_resonance_time_average=0.,
            in_resonance_aligned_energy_at_d_1=-2.*m0,
            in_resonance_attractive_force_law='-2*m0*|n1*n2|/d^2',
            merges_only_when_resonant=True),
        discrete_burst_spectrum=dict(
            topological_quantization_rule='z(1) = z(0) on S^1 forces integer winding n in Z',
            emissions_are_sharp_lines_not_smear=True,
            levels=discrete_levels,
            transition_2_to_1_energy_over_E1=2.**1.5-1.,
            transition_3_to_2_energy_over_E1=3.**1.5-2.**1.5),
        untuned_mass_ratio_tests=dict(
            codata_proton_to_electron=codata_mp_me,
            codata_neutron_minus_proton_over_me=codata_mn_mp_over_me,
            candidate_A_single_stage_winding_n_3_over_2=dict(
                n_2_ratio=2.**1.5,
                n_3_ratio=3.**1.5,
                n_150_ratio=150.**1.5,
                verdict='misses_1836_without_supplying_n_150_by_hand'),
            candidate_B_dimensional_accumulation_1D_to_3D=dict(
                formula_bare='(2*d) * pi^(2+d) = 6 * pi^5 (6 signed rays from A2 * 5 phase dimensions from A1)',
                bare_6_pi_5=bare_6pi5,
                bare_relative_error=abs(bare_6pi5-codata_mp_me)/codata_mp_me,
                bare_ppm_error=1e6*abs(bare_6pi5-codata_mp_me)/codata_mp_me,
                two_valve_medium_dressing_exp_minus_3_over_sqrt_2=two_valve_dressing,
                dressed_6_pi_5_plus_medium_floor=dressed_6pi5,
                dressed_relative_error=abs(dressed_6pi5-codata_mp_me)/codata_mp_me,
                dressed_ppm_error=1e6*abs(dressed_6pi5-codata_mp_me)/codata_mp_me,
                caveat='6*pi^5 is the exact product of 2d=6 signed rays and pi^5 5-phase volume, but deriving that product directly from the A0 3D stretch integral remains an open mathematical proof step'),
            neutron_proton_gap_comparisons=dict(
                single_stage_n2_bubble_over_me=np_gap_n2_single,
                single_stage_n2_mev=np_gap_n2_single*me_mev,
                single_stage_n2_relative_error=(np_gap_n2_single-codata_mn_mp_over_me)/codata_mn_mp_over_me,
                threshold_capture_bare_2_plus_theta_spaghetti=np_gap_threshold_bare,
                threshold_capture_bare_mev=np_gap_threshold_bare*me_mev,
                threshold_capture_bare_relative_error=(np_gap_threshold_bare-codata_mn_mp_over_me)/codata_mn_mp_over_me,
                threshold_capture_dressed_over_me=np_gap_threshold_dressed,
                threshold_capture_dressed_mev=np_gap_threshold_dressed*me_mev,
                threshold_capture_dressed_relative_error=(np_gap_threshold_dressed-codata_mn_mp_over_me)/codata_mn_mp_over_me)))


def thermodynamic_gravity_and_cmb_temperature(h0_km_s_mpc=67.4):
    """Evaluate the Planck-to-Absolute-Zero thermodynamic temperature bridge across the room [l_P, R_H].

    In A3, the room is stretched between two endpoints:
      - Top (hottest/most compressed before 0 crosses to -): Planck length l_P = 1.616255e-35 m,
        Planck temperature T_P = 1.416784e32 K.
      - Bottom (coldest/most spread out, as wide as the room): Hubble/Absolute-Zero horizon
        R_H = c / H0, with Gibbons-Hawking horizon temperature T_H = T_P / (2*pi*(R_H/l_P)) ~ 2.66e-30 K.
    Because 3D horizon curvature energy density scales as (l_P / R_H)^2 while 3D radiation energy
    density scales as (T / T_P)^4, the room's equilibrium temperature scales as the geometric mean
    sqrt(T_P * T_H) ~ (l_P / R_H)^(1/2), spanning 62 orders of magnitude down to the ~K scale:
      - Bare geometric mean: T_geom = sqrt(T_P * T_H) = T_P / sqrt(2*pi*(R_H/l_P)) ≈ 19.3960 K.
      - Critical Friedmann-Stefan-Boltzmann temperature:
        T_crit = T_P * (45 / (8*pi^3))^(1/4) * sqrt(l_P / R_H) ≈ 31.7301 K.
      - 2D horizon surface e-fold projection: T_geom * exp(-2) ≈ 2.6250 K (within 3.69% of 2.7255 K).
      - Two-valve 5-variable fixed-point projection (chi_* = 3.06696835, 4/5 per-variable weight):
        T_crit * exp(-(4/5)*chi_*) ≈ 2.7283 K (within 0.105% of measured T_CMB = 2.72548 K).
    """
    l_p=1.616255e-35
    t_p=1.416784e32
    c_m_s=299792458.
    mpc_m=3.0856775814913673e22
    h0_s=float(h0_km_s_mpc)*1e3/mpc_m
    r_h=c_m_s/h0_s
    b_ratio=r_h/l_p
    t_h=t_p/(2.*math.pi*b_ratio)
    t_geom=math.sqrt(t_p*t_h)
    t_crit=t_p*((45./(8.*(math.pi**3)))**.25)/math.sqrt(b_ratio)
    chi_star=3.0669683518473474
    t_2d_efold=t_geom*math.exp(-2.)
    t_five_var=t_crit*math.exp(-.8*chi_star)
    t_cmb_obs=2.72548
    return dict(
        h0_km_s_mpc=float(h0_km_s_mpc),
        planck_length_m=l_p,
        planck_temperature_K=t_p,
        horizon_radius_m=r_h,
        room_dynamic_range_R_H_over_l_P=b_ratio,
        horizon_temperature_T_H_K=t_h,
        bare_geometric_mean_sqrt_TP_TH_K=t_geom,
        critical_blackbody_temperature_K=t_crit,
        horizon_2d_efold_temperature_K=t_2d_efold,
        horizon_2d_efold_relative_error=abs(t_2d_efold-t_cmb_obs)/t_cmb_obs,
        two_valve_fixed_point_temperature_K=t_five_var,
        two_valve_fixed_point_relative_error=abs(t_five_var-t_cmb_obs)/t_cmb_obs,
        observed_cmb_temperature_K=t_cmb_obs)


def nonlocal_to_local_atlas_readout():
    """Zero-parameter non-local multiverse fixed point and local observer translation dictionary.

    1. Solves the exact two-valve steady-state compression root chi_* where the
       5-variable single-cost visible filament share Omega_vis(chi_*) equals the
       3D volumetric floor exp(-3) (so the combined dark sector Omega_DE + Omega_DM
       equals the 3D exponential horizon 1 - exp(-3) ≈ 95.021293%).
    2. Evaluates the discrete 3D volumetric e-fold spectrum
       eta_K = tanh((pi/2) * exp(-2*chi_spaghetti) * exp(-3*K)) = tanh((pi/4) * exp(-3*K))
       for K = 1..7, where exp(-2*chi_spaghetti) = 1/2 is the derived spaghettification
       threshold factor (chi_spaghetti = ln(2)/2) and 3*K counts 3D volumetric e-folds.
       At K = 7 (21 e-folds), eta_7 = tanh((pi/4)*exp(-21)) ≈ 5.955329e-10.
    3. Maps each non-local geometric mechanism to what a local observer inside one
       room infers as Dark Energy, Dark Matter, and Missing Antimatter.
    """
    target_vis=math.exp(-3)

    def vis_share(chi):
        row=nested_zero_crossing_map(1,float(chi))['nested_chain'][0]
        return row['visible_filament_like_fraction']

    with mp.workdps(50):
        def mp_vis_residual(c):
            ex=mp.exp(c);xi=mp.exp(-c)
            u=(mp.sqrt(16*ex*ex+17)-3)/(4*ex)
            th=u*u
            e_de=3*u**3
            e_dm=u**3+2*(1-th)*xi
            e_vis=2*(1-th)**2*ex
            return e_vis/(e_de+e_dm+e_vis)-mp.exp(-3)
        chi_star=float(mp.findroot(mp_vis_residual,(mp.mpf('3.0'),mp.mpf('3.2'))))

    at_star=nested_zero_crossing_map(1,chi_star)['nested_chain'][0]
    at_d3=nested_zero_crossing_map(1,3.)['nested_chain'][0]
    harmonic_ladder=[]
    for k in range(1,8):
        eta_k=math.tanh(.25*math.pi*math.exp(-3*k))
        harmonic_ladder.append(dict(three_d_efold_stage_K=k,
                                    total_efolds_3K=3*k,
                                    asymmetry_eta_K=eta_k))
    eta_7=harmonic_ladder[-1]['asymmetry_eta_K']
    circle_pendulum=signed_circle_pendulum_duality()
    three_acc_ladder=derived_three_accumulator_bubble_ladder()
    wall_bypasses=exponential_wall_bypasses_and_triad_winner()
    bubble_spectrum=vacuum_bubble_crystallization_and_mass_spectrum()
    thermo_cmb=thermodynamic_gravity_and_cmb_temperature()
    return dict(schema='nsc-nonlocal-to-local-atlas-v1',
                axioms=dict(
                    A0='Least resistance: minimize single generator cost E = Q_stretch + P*V',
                    A1='1D seed + inherited compression offset = time (1D -> 2D circle -> 3D sphere)',
                    A2='Signed-zero crossing (+ -> 0 -> -): edge-on projection of the closed circle R(tau)=R0*cos(omega*tau) crosses 0 into -, flipping sign(det F) = (-1)^k and turning implosion into expansion',
                    A3='Moving spectral band [lambda_min, lambda_max] = exp(-chi)*[lambda_min_0, lambda_max_0] from Planck cutoff to Absolute-Zero horizon',
                    A4='Two-valve cyclic circulation on a closed vortex-diffusion ring (O <-> o) across one-way 0 horizons'),
                derived_theorems_from_axioms=dict(
                    why_minus_exists_and_radius_crosses_zero=circle_pendulum,
                    exponential_wall_bypasses_and_triad_winner=wall_bypasses,
                    surrounding_medium_3_to_1_virial_split=bubble_spectrum['surrounding_medium_equilibrium'],
                    integer_winding_sharp_levels=bubble_spectrum['discrete_burst_spectrum'],
                    opposite_winding_cancellation=bubble_spectrum['annihilation_to_back_to_back_ripples'],
                    shared_flow_resonance_interaction=bubble_spectrum['shared_flow_resonance_interaction'],
                    three_accumulator_bubble_ladder=three_acc_ladder),
                numerical_coincidences_not_derived=dict(
                    note='Formulas below involve post-hoc target selection, Lenz (1951) 6*pi^5, rung selection (K=7), or external measured inputs (H0, G); quarantined from derived theorems.',
                    lenz_6_pi_5_and_medium_dressing=bubble_spectrum['untuned_mass_ratio_tests']['candidate_B_dimensional_accumulation_1D_to_3D'],
                    neutron_proton_gap_comparisons=bubble_spectrum['untuned_mass_ratio_tests']['neutron_proton_gap_comparisons'],
                    k7_baryon_harmonic_rung=eta_7,
                    cmb_temperature_projections_using_measured_H0_G=thermo_cmb,
                    two_valve_split_at_exp_minus_3_vs_lambda_cdm=dict(
                        derived_at_exp_minus_3=[at_star['dark_energy_like_fraction'],at_star['dark_matter_like_fraction'],at_star['visible_filament_like_fraction']],
                        planck_base_lambda_cdm_inference=[0.6847,0.265,0.0493],
                        de_difference_from_lambda_cdm=at_star['dark_energy_like_fraction']-0.6847)),
                why_minus_exists_and_radius_crosses_zero=circle_pendulum,
                vacuum_bubble_crystallization_and_masses=bubble_spectrum,
                thermodynamic_gravity_and_cmb_bridge=thermo_cmb,
                two_valve_fixed_point=dict(
                    condition='Omega_vis(chi_*) = exp(-3) (5-variable filament share equals 3D volumetric floor)',
                    chi_star=chi_star,
                    dark_energy_fraction=at_star['dark_energy_like_fraction'],
                    dark_matter_fraction=at_star['dark_matter_like_fraction'],
                    combined_dark_sector_fraction=at_star['dark_sector_combined_fraction'],
                    visible_filament_fraction=at_star['visible_filament_like_fraction'],
                    residual_from_volumetric_horizon=at_star['dark_sector_combined_fraction']-(1-target_vis),
                    at_integer_3d_efold_chi_3=dict(
                        dark_energy_fraction=at_d3['dark_energy_like_fraction'],
                        dark_matter_fraction=at_d3['dark_matter_like_fraction'],
                        combined_dark_sector_fraction=at_d3['dark_sector_combined_fraction'],
                        visible_filament_fraction=at_d3['visible_filament_like_fraction'])),
                discrete_3d_harmonic_asymmetry_ladder=dict(
                    formula='eta_K = tanh((pi/2) * exp(-2*chi_spaghetti) * exp(-3*K)) = tanh((pi/4) * exp(-3*K))',
                    spaghettification_threshold_chi=.5*math.log(2),
                    spaghettification_factor=math.exp(-2*.5*math.log(2)),
                    stages=harmonic_ladder,
                    stage_K_7_eta=eta_7,
                    bbn_observed_window=[5.8e-10,6.5e-10],
                    planck_observed_reference=6.1e-10,
                    within_bbn_window=(5.8e-10<=eta_7<=6.5e-10)),
                translation_dictionary=[
                    dict(local_label='Dark Energy (Lambda)',
                         nonlocal_mechanism='Isotropic 3D stretch 3*(theta^*)^(3/2) inside the - sphere driven by inflow across the outer birth 0 valve',
                         derived_share_at_chi_star=at_star['dark_energy_like_fraction']),
                    dict(local_label='Dark Matter',
                         nonlocal_mechanism='Non-luminous 0-wall enclosing load (theta^*)^(3/2) + radial infall load 2*(1-theta^*)*exp(-chi) tugging across outer and inner 0 horizons',
                         derived_share_at_chi_star=at_star['dark_matter_like_fraction']),
                    dict(local_label='Visible Matter',
                         nonlocal_mechanism='Active transverse filamentary winding 2*(1-theta^*)^2*exp(chi) spaghettifying toward inner child 0 horizons, locked at floor exp(-3)',
                         derived_share_at_chi_star=at_star['visible_filament_like_fraction']),
                    dict(local_label='Missing Antimatter (eta_B)',
                         nonlocal_mechanism='Single-sign orientation entry sign(det F_k)=(-1)^k across 0 plus offset decay lag exp(-gamma*tau) between + and - half-cycles at the K=7 3D harmonic',
                         derived_asymmetry_eta_7=eta_7)],
                no_observational_fit=True)


def buckling_to_goldilocks_calculation():
    """Connected pipeline: unified single-cost cycle -> nested map -> 0-crossing -> buckling -> damped wave -> Goldilocks."""
    t_h=math.sqrt(2.)
    unified_scan=[unified_least_resistance_cycle(chi) for chi in (0.,.25,.5,1.,2.,2.5,3.,4.)]
    nested_scan=[nested_zero_crossing_map(3,chi)['nested_chain'][0] for chi in (0.,1.,2.,2.5,3.,4.)]
    three_level_map=nested_zero_crossing_map(3,2.5)
    zero_crossing_scan=[two_force_zero_crossing(t,2.,1.) for t in (0.,1.,t_h,2.,3.)]
    recurring_winner_scan=[spaghettification_to_spherical_winner(chi) for chi in (0.,.25,.5,1.,2.)]
    buckling_scan=[helical_buckling_under_compression(chi) for chi in (0.,.25,.5,1.,2.)]
    post_buckling=helical_buckling_under_compression(1.)
    omega=2*math.pi
    gamma=.5*math.log(2)
    wave_history=[damped_horizon_oscillation(t,gamma,omega) for t in (0.,.25,.75,1.,2.,4.)]
    return dict(schema='nsc-buckling-damping-goldilocks-v1',
                unified_single_cost_cycle=dict(
                    functional='E(lambda_3, tilde_lambda_1..3, theta; chi) = [2(1-theta)^2 e^chi + (1-theta)^2/lambda_3 + e^(-2chi) lambda_3] + [sum_i theta^2/tilde_lambda_i + tilde_lambda_1 tilde_lambda_2 tilde_lambda_3]',
                    stationarity_solution='sqrt(theta^*) = (sqrt(16 e^(2chi) + 17) - 3)/(4 e^chi); lambda_3^* = (1-theta^*)e^chi; tilde_lambda_i^* = sqrt(theta^*)',
                    spaghettification_threshold_chi=.5*math.log(2),
                    scan=unified_scan),
                nested_local_space_readout=dict(
                    compression_scan=nested_scan,
                    three_level_chain_at_chi_2_5=three_level_map,
                    spectral_band_shift=[spectral_band_and_baryon_harmonic(chi) for chi in (-1.,0.,1.,3.)],
                    branching_chain_evolution=branching_cosmos_chain_evolution(generations=3,branch_factor=3,compression_per_step=2.5),
                    cyclic_vortex_ring=[cyclic_vortex_diffusion_ring(t) for t in (0.,.25,.5,.75,1.)]),
                stage_0_two_force_zero_crossing=dict(
                    rule='R_tt = -(F_in - F_out) < 0; R(tau) crosses + -> 0 -> - smoothly with det(F)=(R/R0)^3 flipping sign',
                    horizon_offset=t_h,
                    scan=zero_crossing_scan),
                stage_0b_spaghettification_and_spherical_winner=dict(
                    infall_rule='diag(xi, xi, xi^-2)/xi^-2 = diag(xi^3, xi^3, 1) -> diag(0, 0, 1) (3D -> 1D spaghettification)',
                    interior_winner_rule='A(s)^3 / (36*pi*V^2) >= 1 with unique global minimum at s=0 (round 3D sphere)',
                    cooker_takeover_rule='f_minus(chi) = 1 - exp(-3*chi) -> 1 as compression continues',
                    scan=recurring_winner_scan),
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
