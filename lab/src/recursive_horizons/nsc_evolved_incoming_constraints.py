"""Incoming matter contractions for source-fixed, history-evolved columns.

Momentum acts on the actual axial dependence of the columns. Source-energy
labels are never substituted for outgoing momentum. The old frozen-state
assembly remains an unchanged regression owner.
"""
from hashlib import sha256
import numpy as np

from .nsc_common_subtracted_ks_source import raw_ks_vertex_coefficients
from .nsc_prepared_history_jets import restriction_covariance_tangent


def _adjoint(value):
    return value.swapaxes(-1, -2).conj()


def _finite(value, name, *, real=False):
    if value is None:
        raise ValueError('explicit '+name+' required')
    raw = np.asarray(value)
    if real and np.iscomplexobj(raw) and np.any(raw.imag != 0):
        raise ValueError('real '+name+' required')
    result = np.asarray(raw.real if real else raw, dtype=float if real else complex)
    if not np.isfinite(result).all():
        raise ValueError('finite '+name+' required')
    return result


def source_column_matter(columns, axial_columns, source_covariance,
                         column_tangents, axial_tangents, *, mass, angular,
                         axial_scale, radius, multiplicity, tolerance=3e-11):
    """Raw N,beta matter insertion and its complete fixed-source tangent.

    F,F_z have shape (nz,2,nsrc); dF,dF_z have shape (ndir,nz,2,nsrc).
    Columns already contain sqrt(dE/(2*pi)) quadrature weights ONCE. Csrc
    is the unchanged full source matrix, including its coherences. The
    caller supplies the actual signed-family multiplicity from the ledger.
    This finite-column contraction alone is not a complete physical source.
    """
    F = _finite(columns, 'weighted source columns')
    Fz = _finite(axial_columns, 'actual axial derivatives')
    C = _finite(source_covariance, 'fixed source covariance')
    dF = _finite(column_tangents, 'retarded column tangents')
    dFz = _finite(axial_tangents, 'axial derivatives of column tangents')
    scalars = _finite([mass, angular, axial_scale, radius, multiplicity, tolerance], 'channel/geometry', real=True)
    mass, angular, axial_scale, radius, multiplicity, tolerance = scalars
    if mass < 0 or min(axial_scale, radius, multiplicity, tolerance) <= 0:
        raise ValueError('owned nonnegative mass and positive geometry/multiplicity/tolerance required')
    if F.ndim != 3 or F.shape[1] != 2 or min(F.shape[0], F.shape[2]) < 1 or Fz.shape != F.shape:
        raise ValueError('matching (nz,2,nsrc) columns and axial derivatives required')
    if dF.ndim != 4 or dF.shape[1:] != F.shape or dFz.shape != dF.shape:
        raise ValueError('matching (ndir,nz,2,nsrc) retarded tangents required')
    if C.shape != (F.shape[2], F.shape[2]) or np.max(abs(C-_adjoint(C))) > tolerance:
        raise ValueError('Hermitian full fixed-source matrix required')
    eig = np.linalg.eigvalsh(C)
    if eig.min() < -tolerance or eig.max() > 1+tolerance:
        raise ValueError('fixed source covariance must satisfy its CAR interval')
    density = F@C@_adjoint(F)
    momentum = (Fz@C@_adjoint(F)-F@C@_adjoint(Fz))/(2j)
    density_tangent = restriction_covariance_tangent(F,dF,C,np.zeros_like(C))
    momentum_tangent = (dFz@C@_adjoint(F)+Fz@C@_adjoint(dF)
                        -dF@C@_adjoint(Fz)-F@C@_adjoint(dFz))/(2j)

    vertices = raw_ks_vertex_coefficients(np.array([1.,0.,axial_scale,radius]),
                                         mass,angular,envelopes=np.ones(4))
    def contract(rho, current):
        return -multiplicity*(
            np.einsum('bij,...ji->...b',vertices.multiplication[:2],rho)
            +np.einsum('bij,...ji->...b',vertices.momentum[:2],current))

    raw, tangent = contract(density, momentum), contract(density_tangent, momentum_tangent)
    imag = max(float(np.max(abs(raw.imag))), float(np.max(abs(tangent.imag), initial=0.)))
    if imag > tolerance:
        raise ArithmeticError('non-real symmetric matter insertion exceeds algebra tolerance')
    return {'action_gradient': raw.real, 'action_gradient_tangent': tangent.real,
            'density': density, 'momentum_density': momentum,
            'density_tangent': density_tangent, 'momentum_density_tangent': momentum_tangent,
            'constraint_order': ('N', 'beta'), 'imaginary_residual': imag,
            'source_covariance_tangent': 'zero by unchanged upstream/source law',
            'momentum_prescription': '-i/2*(F_z Csrc Fdagger-F Csrc F_zdagger)',
            'outgoing_momentum_replaced_by_source_energy': False,
            'spectral_and_field_error_bound': None, 'physical_constraint_status': 'OPEN'}


def column_derivatives(columns, tangents, derivative):
    """Apply one explicitly supplied fixed-z differentiation operator.

    Its consistency/accuracy belongs to the supplied operator; no finite
    difference matrix is silently treated as an exact continuum derivative.
    """
    F, dF = _finite(columns, 'columns'), _finite(tangents, 'tangents')
    D = _finite(derivative, 'axial differentiation operator', real=True)
    if F.ndim != 3 or dF.ndim != 4 or dF.shape[1:] != F.shape or D.shape != (len(F), len(F)):
        raise ValueError('fixed axial derivative and column dimensions must agree')
    return np.einsum('ij,jas->ias', D, F), np.einsum('ij,djas->dias', D, dF)


def surface_geometry_response(slots, slot_tangents, coefficients):
    """Reuse compatible-family local/reference changes, excluding baseline.

    Slots are (w,w_z,w_zz,w_zzz,U,U_z). D[0] is deliberately zero: the
    old constant S_N/S_beta is not a geometry coefficient or the new state.
    Coefficient error enclosures remain the caller's recorded evidence.
    """
    s = _finite(slots, 'compatible surface slots', real=True)
    ds = _finite(slot_tangents, 'compatible slot tangents', real=True)
    if s.ndim != 2 or s.shape[1] != 6 or ds.ndim != 3 or ds.shape[1:] != s.shape:
        raise ValueError('(nz,6) slots and (ndir,nz,6) tangents required')
    required = ('a', 'r', 'Hr', 'A0', 'A1', 'd', 'e', 'C')
    a, r, Hr, A0, A1, d, e, c = _finite([coefficients[k] for k in required], 'owned surface coefficients', real=True)
    D = _finite(coefficients['D'], 'nonconstant lapse polynomial', real=True)
    P = _finite(coefficients['F'], 'shift polynomial', real=True)
    if min(a,r) <= 0 or D.shape != (5,) or D[0] != 0 or P.shape != (3,):
        raise ValueError('positive fixed geometry, D[0]=0 and owned polynomial degrees required')
    w, wz, wzz, wzzz, U, Uz = s.T
    A = A0+A1*w
    B = -(2*A+(Hr+w/r)*d)/a**2
    Bw = -(2*A1+d/r)/a**2
    poly = np.polynomial.polynomial
    F = poly.polyval(w, P)
    Dprime = poly.polyval(w, poly.polyder(D))
    Fprime = poly.polyval(w, poly.polyder(P))
    value = np.stack((A*U+B*wzz+c*wz*wz+poly.polyval(w,D),
                      d*Uz+e*wzzz+F*wz), axis=-1)
    jac = np.zeros((len(s),2,6))
    jac[:,0,:] = np.stack((A1*U+Bw*wzz+Dprime,2*c*wz,B,np.zeros_like(w),A,np.zeros_like(w)),axis=-1)
    jac[:,1,:] = np.stack((Fprime*wz,F,np.zeros_like(w),np.full_like(w,e),np.zeros_like(w),np.full_like(w,d)),axis=-1)
    return {'action_gradient_change': value, 'action_gradient_tangent': np.einsum('zbs,dzs->dzb',jac,ds),
            'jacobian': jac, 'constraint_order': ('N','beta'),
            'constant_frozen_matter_insertion_included': False,
            'coefficient_error_bound': None, 'physical_constraint_status': 'OPEN'}


def compatible_history_slots(provider, z, direction_count):
    """Pull the actual radius history into its six owned incoming slots."""
    from .nsc_compatible_history_geometry import CompatibleIncomingMetric
    from .nsc_evolved_incoming_state import CachedCompatibleIncomingMetric, AmplitudeOnlyMetric
    family = provider
    if isinstance(family,AmplitudeOnlyMetric):
        family = family.inner
    if isinstance(family,CachedCompatibleIncomingMetric):
        family = family.original
    if not isinstance(family,CompatibleIncomingMetric):
        raise TypeError('owned compatible history family required')
    z = _finite(z,'incoming coordinates',real=True)
    names = ('w','w_z','w_zz','w_zzz','U','U_z')
    rows = []
    for direction in family.directions:
        values = direction.incoming_slots(z)
        rows.append(np.stack([values[name] for name in names],axis=-1))
    basis = np.asarray(rows).reshape(len(family.directions),len(z),6)
    value = np.einsum('d,dzs->zs',np.asarray(family.amplitudes),basis)
    if direction_count == 0:
        tangent = np.zeros((0,len(z),6))
    elif direction_count == len(family.directions):
        tangent = basis
    else:
        raise ValueError('state and geometry directions do not describe the same history')
    return value,tangent


def evolved_constraint_diagnostic(state, provider, axial_derivative, reference_columns,
                                  reference_axial_columns, baseline_gradient,
                                  coefficients, channel, angular_sign):
    """Connect evolved matter and delta C to both included constraints.

    The reference columns are the fixed computational reference for the
    baseline subtraction; their tangent is ZERO. They are not imposed as the
    physical incoming state. The finite-channel output is always explicitly
    OPEN until full state-correction coverage and errors are supplied.
    """
    from .nsc_evolved_incoming_state import EvolvedIncomingState
    if not isinstance(state,EvolvedIncomingState):
        raise TypeError('evolved incoming state required; frozen C0 and legacy matter dictionaries are regression only')
    state.validate()
    state.require_history(provider,state.history.times,state.history.spatial_grid)
    from .nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain
    N,beta,q,r = state.history.rho1_metric
    intrinsic = TransmittingDiracSeamDomain(lapse=N,radial_scale=q,shift=beta,radius=r)
    if not np.allclose([coefficients['a'],coefficients['r']],
                       [intrinsic.induced_axial_scale,r],atol=3e-12,rtol=0):
        raise ValueError('constraint coefficients and evolved state must use the same intrinsic identification')
    F, dF = state.weighted_columns,state.weighted_column_tangents
    Fz, dFz = column_derivatives(F,dF,axial_derivative)
    ref = _finite(reference_columns,'fixed reference source columns')
    refz = _finite(reference_axial_columns,'fixed reference axial columns')
    if ref.shape != state.columns.shape or refz.shape != ref.shape:
        raise ValueError('reference and evolved source columns must share the source/coordinate inventory')
    angular = float(channel['angular_eigenvalue'])
    signs = (1,) if angular == 0 else (-1,1)
    if angular_sign not in signs:
        raise ValueError('owned angular sign required; no energy-sign folding is performed')
    copies,deg = channel['copy_count'],channel['degeneracy']
    if any(isinstance(v,(bool,np.bool_)) or int(v)!=v or v<=0 for v in (copies,deg)):
        raise ValueError('positive integer copies/degeneracy from the retained ledger required')
    parameters = dict(mass=float(channel['compact_mass']),angular=angular_sign*angular,
                      axial_scale=float(coefficients['a']),radius=float(coefficients['r']),
                      multiplicity=copies*deg/len(signs))
    current = source_column_matter(F,Fz,state.source_covariance,dF,dFz,**parameters)
    zeros = np.zeros((0,*ref.shape),complex)
    reference = source_column_matter(ref*state.column_weights,refz*state.column_weights,
                                    state.source_covariance,zeros,zeros,**parameters)
    slots,slot_tangents = compatible_history_slots(provider,state.z,len(dF))
    geometric = surface_geometry_response(slots,slot_tangents,coefficients)
    baseline = _finite(baseline_gradient,'baseline complete-assembly approximant',real=True)
    if baseline.shape != (2,):
        raise ValueError('baseline N,beta assembly approximant required')
    correction = current['action_gradient']-reference['action_gradient']
    result = baseline+geometric['action_gradient_change']+correction
    tangent = geometric['action_gradient_tangent']+current['action_gradient_tangent']
    return {'constraint_order':('N','beta'),'z':np.array(state.z,copy=True),
            'partial_constraint_diagnostic':result,'history_derivative_diagnostic':tangent,
            'baseline_action_gradient':baseline,
            'local_reference_gradient_change':geometric['action_gradient_change'],
            'local_reference_gradient_tangent':geometric['action_gradient_tangent'],
            'evolved_matter_correction':correction,
            'evolved_matter_tangent':current['action_gradient_tangent'],
            'incoming_slots':slots,'incoming_slot_tangents':slot_tangents,
            'source_binding':state.history.source_digest,'history_binding':state.history.digest,
            'state_law':'C_Sigma[g]=U_g C_up U_gdagger; fixed source, evolved incoming state',
            'reference_state_tangent_subtracted':False,
            'constant_frozen_matter_used_as_physical_state':False,
            'sampled_channel':{'group':channel.get('index'),'angular_sign':angular_sign,
                               'multiplicity':parameters['multiplicity']},
            'physical_constraint_status':'OPEN','certified_constraint_residual':None,
            'full_source_error_bound':None,
            'unresolved':['baseline low-energy/subgap source accuracy',
                          'complete changed-history source-energy and signed-family coverage',
                          'continuum field and axial-derivative error',
                          'coefficient interval propagation for a claimed residual'],
            'spectrum_coverage':dict(state.coverage),
            'stationarity_tolerance':3e-11,'physical_history_selected':False,
            'extended_stationarity':'OPEN','metric_timestep':False}


def evaluate_source_fixed_history(owner, times, provider, initial_fields, source,
                                  axial_derivative, baseline_gradient, coefficients,
                                  channel, angular_sign):
    """Fresh evolution -> evolved state -> both constraints and history tangent.

    This retarded-patch entry accepts the owned fixed source preparation,
    never an incoming C0 or an already frozen matter insertion. Initial and
    incident column tangents are fixed to zero by the unchanged-past class.
    The compact support/past condition remains a physical caller hypothesis;
    the initial grid, channel, and inflow consistency are checked here.
    """
    from .nsc_evolved_incoming_state import FixedSourcePreparation,evolve_incoming_state
    if not isinstance(source,FixedSourcePreparation):
        raise TypeError('fixed upstream/source preparation required, not incoming C0')
    t = _finite(times,'history time grid',real=True)
    phi = _finite(initial_fields,'owned initial source columns')
    x = np.asarray(owner.x)
    if t.ndim!=1 or len(t)<2 or np.any(np.diff(t)<=0):
        raise ValueError('strictly increasing history time grid required')
    if phi.shape!=(2*len(x),len(source.energies)):
        raise ValueError('complete initial columns on the owned spatial/source grid required')
    nsign=1 if float(channel['angular_eigenvalue'])==0 else 2
    if (angular_sign not in ((1,) if nsign==1 else (-1,1))
            or owner.mass!=float(channel['compact_mass'])
            or owner.angular!=angular_sign*float(channel['angular_eigenvalue'])):
        raise ValueError('propagator and constraint vertices must use the same retained channel')
    metric = _finite(provider.values(t[0],x),'initial history metric',real=True)
    if metric.shape!=owner.reference.shape or np.max(abs(metric-owner.reference))>3e-11:
        raise ValueError('retarded patch must start in the unchanged reference geometry')
    directions = _finite(provider.log_directions(t[0],x,metric),'initial metric directions',real=True)
    if directions.ndim!=3 or directions.shape[1:]!=metric.shape or np.max(abs(directions),initial=0.)>3e-11:
        raise ValueError('unchanged-past preparation requires zero initial metric tangents')
    ndir=len(directions)
    initial_tangents=np.zeros((ndir,*phi.shape),complex)
    right=np.array([len(x)-1,2*len(x)-1])
    N,beta,q,_=owner.reference
    speeds=np.array([N[-1]/q[-1]-beta[-1],-N[-1]/q[-1]-beta[-1]])
    sat=-speeds/owner.weights[-1]
    def incoming(time):
        forcing=np.zeros_like(phi)
        forcing[right]=sat[:,None]*phi[right]*np.exp(-1j*source.energies*(time-t[0]))
        return forcing,np.zeros_like(initial_tangents)
    state=evolve_incoming_state(owner,t,provider,phi,initial_tangents,incoming,source)
    index=state.history.rho1_index
    fixed_initial_trace=state.history.restriction@phi[[index,len(x)+index]]
    reference=fixed_initial_trace[None]*np.exp(-1j*(t-t[0])[:,None,None]*source.energies)
    reference_z=-1j*source.energies*reference  # unchanged stationary reference only
    result=evolved_constraint_diagnostic(state,provider,axial_derivative,reference,reference_z,
        baseline_gradient,coefficients,channel,angular_sign)
    result['fresh_history_evolution_executed']=True
    result['source_preparation_derivatives_fixed_to_zero']=True
    result['fixed_preparation_binding']=sha256(
        source.digest.encode()+np.ascontiguousarray(phi).tobytes()
        +np.ascontiguousarray(x).tobytes()
        +repr((owner.mass,owner.angular,float(t[0]))).encode()).hexdigest()
    return state,result
