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

from pathlib import Path as _receipt_Path
from .nsc_frozen_sources import selected_source_sha256 as _receipt_source_sha256

# Retained executable control: constraint_control
# Folded from scripts/derive_nsc_evolved_incoming_constraints.py at 329ef809; original receipts retain their producer identities.
import argparse as _constraint_control_argparse
from hashlib import sha256 as _constraint_control_sha256
import json as _constraint_control_json
from pathlib import Path as _constraint_control_Path
import sys as _constraint_control_sys
import numpy as _constraint_control_np
from scipy.interpolate import make_interp_spline as _constraint_control_make_interp_spline
from recursive_horizons.nsc_posix import ROOT as _constraint_control_ROOT
_constraint_control_sys.path.insert(0, str(_constraint_control_ROOT / 'src'))
_constraint_control_sys.path.insert(0, str(_constraint_control_ROOT / 'scripts'))
from recursive_horizons.nsc_evolved_incoming_constraints import evolved_constraint_diagnostic as _constraint_control_evolved_constraint_diagnostic, source_column_matter as _constraint_control_source_column_matter, column_derivatives as _constraint_control_column_derivatives
from recursive_horizons.nsc_evolved_incoming_state import EvolvedIncomingState as _constraint_control_EvolvedIncomingState, FixedSourcePreparation as _constraint_control_FixedSourcePreparation, HistoryBinding as _constraint_control_HistoryBinding, _coverage as _constraint_control__coverage
from recursive_horizons.nsc_evolved_incoming_state import state_control as _constraint_control_STATE
from derive_nsc_incoming_source_update_v3 import authenticate as _constraint_control_authenticate
_constraint_control_OUTPUT = 'results/development/nsc-evolved-incoming-constraints.json'
_constraint_control_SOURCES = ('src/recursive_horizons/nsc_evolved_incoming_constraints.py', 'scripts/derive_nsc_evolved_incoming_constraints.py', 'tests/test_nsc_evolved_incoming_constraints.py', 'docs/nsc-evolved-incoming-constraints.md')
_constraint_control_RECORDS = ('results/development/nsc-evolved-incoming-state.json', 'results/development/nsc-incoming-source-update-v5.json', 'results/development/nsc-incoming-surface-coefficients.json', 'results/development/nsc-incoming-surface-regular-branch.json', 'results/development/nsc-incoming-fixed-transfer.json')
_constraint_control_OWNERS = ('src/recursive_horizons/nsc_common_subtracted_ks_source.py', 'src/recursive_horizons/nsc_prepared_history_jets.py', 'src/recursive_horizons/nsc_evolved_incoming_state.py', 'scripts/derive_nsc_evolved_incoming_state.py', 'scripts/derive_nsc_incoming_source_update_v3.py', 'results/development/nsc-retarded-compatible-response.json')
_constraint_control_MAX = lambda x: float(_constraint_control_np.max(abs(x), initial=0.0))

def _constraint_control_digest(path):
    candidate = _receipt_Path(path)
    relative = candidate.relative_to(_constraint_control_ROOT).as_posix() if candidate.is_absolute() else candidate.as_posix()
    return _receipt_source_sha256(_constraint_control_ROOT, relative)

def _constraint_control_coefficients_from_records(coeff, branch):
    c, b = (coeff['coefficients'], branch['branch'])
    mid = lambda x: float(_constraint_control_np.mean(x))
    return {'a': mid(b['intrinsic_intervals']['a']), 'r': mid(b['intrinsic_intervals']['r']), 'Hr': mid(b['intrinsic_intervals']['Hr']), 'A0': mid(b['imported_unchanged']['A0']), 'A1': mid(b['intervals']['A1_total']), 'd': mid(b['imported_unchanged']['d']), 'e': mid(b['imported_unchanged']['e']), 'C': mid(c['C_interval']), 'D': _constraint_control_np.r_[0.0, [mid(c['polynomial_intervals']['D'][str(j)]) for j in range(1, 5)]], 'F': _constraint_control_np.array([mid(c['polynomial_intervals']['F'][str(j)]) for j in range(3)])}

def _constraint_control_load_inputs():
    records = [_constraint_control_json.loads((_constraint_control_ROOT / path).read_text()) for path in _constraint_control_RECORDS]
    for record in records:
        _constraint_control_authenticate(record)
    state_record, baseline, coefficient, branch, transfer = records
    if not state_record['status'].startswith('PASS:') or state_record['control']['field_runs'] != 3:
        raise ValueError('completed nonzero-history state control required')
    if baseline['partial_error_budget']['full_source_error_bound'] is not None:
        raise ValueError('this control records the current explicit baseline accuracy hole')
    with _constraint_control_np.load(_constraint_control_ROOT / state_record['payload']['path'], allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    meta = _constraint_control_json.loads(arrays['metadata_json'].tobytes())
    channel = transfer['actual_source_bindings']['channel']
    if channel['index'] != 14 or channel['compact_mass'] != arrays['source/mass'].item() or channel['angular_eigenvalue'] != arrays['source/angular'].item():
        raise ValueError('same authenticated signed group14_1 operator required')
    source = _constraint_control_FixedSourcePreparation(arrays['source/covariance'], arrays['source/column_weights'], arrays['source/energies'])
    return (records, arrays, meta, channel, source, _constraint_control_coefficients_from_records(coefficient, branch))

def _constraint_control_restore_saved_control(arrays, meta, source, name, amplitude):
    """Only called after authentication of the actual evolved-control payload."""
    times, x = (arrays['times'], arrays['x'])
    rho1 = int(_constraint_control_np.flatnonzero(x == 1.0)[0])
    binding = _constraint_control_HistoryBinding(amplitudes=(amplitude,), inner_radii=(0.007,), outer_radii=(0.03,), times=times, spatial_grid=x, rho1_index=rho1, rho1_metric=arrays['center/rho1_metric'], restriction=arrays['center/restriction'], restriction_parameter_derivative=0.0, S1=meta['support']['S1'], source_digest=source.digest)
    fields = arrays[name + '/columns']
    tangents = arrays['center/column_tangents'][None] if name == 'center' else _constraint_control_np.zeros((0, *fields.shape), complex)
    return _constraint_control_EvolvedIncomingState(arrays['center/z'], fields, tangents, source.covariance, source.column_weights, source.energies, binding, _constraint_control__coverage(fields.shape[-1], 4), {'restored_from_authenticated_evolution': True, 'case': name})

def _constraint_control_calculate():
    records, a, meta, channel, source, coeff = _constraint_control_load_inputs()
    state_record, baseline = records[:2]
    z = a['center/z']
    times = a['times']
    Dz = _constraint_control_make_interp_spline(z, _constraint_control_np.eye(len(z)), k=3).derivative()(z)
    x, phi, data, _ = _constraint_control_STATE.RETARDED.inputs()
    if not _constraint_control_np.array_equal(x, a['x']):
        raise ValueError('same original upstream grid required')
    ix = int(_constraint_control_np.flatnonzero(x == 1.0)[0])
    rows = _constraint_control_np.array([ix, len(x) + ix])
    ref0 = a['center/restriction'] @ phi[rows]
    reference = ref0[None] * _constraint_control_np.exp(-1j * times[:, None, None] * source.energies[None, None, :])
    reference_z = -1j * source.energies * reference
    alpha = state_record['control']['amplitude']
    h = state_record['control']['finite_difference_step']
    outputs = {}
    states = {}
    for name, amplitude in (('center', alpha), ('plus', alpha + h), ('minus', alpha - h)):
        state = _constraint_control_restore_saved_control(a, meta, source, name, amplitude)
        states[name] = state
        outputs[name] = _constraint_control_evolved_constraint_diagnostic(state, _constraint_control_STATE.family(amplitude, meta['support']), Dz, reference, reference_z, baseline['baseline']['action_gradient_approximant'], coeff, channel, 1)
    center = outputs['center']
    fd = (outputs['plus']['partial_constraint_diagnostic'] - outputs['minus']['partial_constraint_diagnostic']) / (2 * h)
    residual = _constraint_control_MAX(fd - center['history_derivative_diagnostic'][0])
    state = states['center']
    F = state.weighted_columns
    dF = state.weighted_column_tangents
    Fz, dFz = _constraint_control_column_derivatives(F, dF, Dz)
    params = dict(mass=channel['compact_mass'], angular=channel['angular_eigenvalue'], axial_scale=coeff['a'], radius=coeff['r'], multiplicity=channel['copy_count'] * channel['degeneracy'] / 2)
    actual = _constraint_control_source_column_matter(F, Fz, source.covariance, dF, dFz, **params)
    old = _constraint_control_json.loads((_constraint_control_ROOT / 'results/development/nsc-retarded-compatible-response.json').read_text())
    _constraint_control_authenticate(old)
    with _constraint_control_np.load(_constraint_control_ROOT / old['payload']['path'], allow_pickle=False) as archive:
        zero_columns = archive['801_64/trace']
        if not _constraint_control_np.array_equal(archive['801_64/times'], times) or not _constraint_control_np.array_equal(archive['source/source'], a['source/blocks']) or (not _constraint_control_np.array_equal(archive['source/weights'], a['source/weights'])):
            raise ValueError('saved zero-history control must share the same time/source inventory')
    zeros = _constraint_control_np.zeros((0, *zero_columns.shape), complex)
    zero_weighted = zero_columns * source.column_weights
    zero_z, _ = _constraint_control_column_derivatives(zero_weighted, zeros, Dz)
    zero_matter = _constraint_control_source_column_matter(zero_weighted, zero_z, source.covariance, zeros, zeros, **params)
    reference_matter = _constraint_control_source_column_matter(reference * source.column_weights, reference_z * source.column_weights, source.covariance, zeros, zeros, **params)
    drift = zero_matter['action_gradient'] - reference_matter['action_gradient']
    history_change = actual['action_gradient'] - zero_matter['action_gradient']
    drift_split = {'reference_field_and_derivative_drift_indicator': _constraint_control_np.max(abs(drift), axis=0).tolist(), 'history_change_at_same_discretization': _constraint_control_np.max(abs(history_change), axis=0).tolist(), 'sum_residual': _constraint_control_MAX(drift + history_change - center['evolved_matter_correction']), 'drift_subtracted_from_main_diagnostic': False, 'is_continuum_error_bound': False, 'saved_zero_history_payload': old['payload']}
    wrong = _constraint_control_source_column_matter(F, -1j * source.energies * F, source.covariance, dF, -1j * source.energies * dF, **params)
    diagonal = _constraint_control_source_column_matter(F, Fz, _constraint_control_np.diag(_constraint_control_np.diag(source.covariance)), dF, dFz, **params)
    effects = {'dropped_delta_C': _constraint_control_MAX(center['evolved_matter_tangent']), 'outgoing_k_equals_minus_E': _constraint_control_MAX(actual['action_gradient_tangent'] - wrong['action_gradient_tangent']), 'dropped_source_coherence': _constraint_control_MAX(actual['action_gradient_tangent'] - diagonal['action_gradient_tangent'])}
    rejected = []
    for label, value in (('bare_C0', _constraint_control_np.eye(2) / 2), ('legacy_matter', {'stress_approximant': _constraint_control_np.zeros(4)})):
        try:
            _constraint_control_evolved_constraint_diagnostic(value, _constraint_control_STATE.family(alpha, meta['support']), Dz, reference, reference_z, baseline['baseline']['action_gradient_approximant'], coeff, channel, 1)
        except TypeError:
            rejected.append(label)
        else:
            raise AssertionError('frozen input accepted as evolved state')
    return (records, a, meta, center, residual, effects, rejected, _constraint_control_MAX(Dz @ _constraint_control_np.ones(len(z))), _constraint_control_MAX(_constraint_control_np.einsum('ij,jas->ias', Dz, reference) - reference_z), drift_split)

def _constraint_control_make_record():
    records, arrays, meta, center, residual, effects, rejected, constant_error, ref_derivative_error, drift_split = _constraint_control_calculate()
    derivative_tolerance = 3e-08
    passed = residual < derivative_tolerance and min(effects.values()) > 10 * residual
    return {'schema': 'NSC-EVOLVED-INCOMING-CONSTRAINTS-v1', 'accountable_author': 'Douglas Ek', 'status': 'PASS: evolved-state constraint wiring; full physical constraints OPEN' if passed else 'OPEN: integrated derivative or omitted-term discrimination unresolved', 'state_law': center['state_law'], 'control': {'family': '14_1', 'source_energies': 4, 'amplitude': records[0]['control']['amplitude'], 'finite_difference_step': records[0]['control']['finite_difference_step'], 'field_runs': 0, 'array_only': True, 'axial_derivative': 'derivative of cubic interpolant on saved z nodes', 'source_energy_substitution_used_on_evolved_columns': False, 'reference': 'original stationary canonical columns for fixed computational subtraction only', 'physical_constraint_tolerance': 3e-11, 'derivative_control_tolerance': derivative_tolerance}, 'residuals': {'complete_discrete_constraint_derivative': residual, 'constant_derivative': constant_error, 'reference_axial_derivative_indicator': ref_derivative_error}, 'omitted_term_effects': effects, 'frozen_inputs_rejected': rejected, 'reference_drift_separation': drift_split, 'diagnostics': {'z': center['z'].tolist(), 'partial_constraint_diagnostic': center['partial_constraint_diagnostic'].tolist(), 'history_derivative_diagnostic': center['history_derivative_diagnostic'].tolist(), 'evolved_matter_correction': center['evolved_matter_correction'].tolist(), 'evolved_matter_tangent': center['evolved_matter_tangent'].tolist(), 'maximum_absolute_partial_gradient': _constraint_control_np.max(abs(center['partial_constraint_diagnostic']), axis=0).tolist(), 'sampled_channel': center['sampled_channel']}, 'error_budget': {'baseline_covered_region_error': records[1]['partial_error_budget']['covered_spectral_regions_action_error_upper'], 'baseline_full_source_error': None, 'changed_history_complete_source_error': None, 'field_and_axial_derivative_error': None, 'coefficient_interval_propagation': None}, 'scope': {'fixed_C0_regression_unchanged': True, 'fixed_source_and_evolved_incoming': True, 'retarded_delta_C_included': True, 'old_fixed_C0_roots_relabelled': False, 'reference_counted_once': True, 'reference_state_tangent_subtracted': False, 'partial_channel_control_only': True, 'full_constraint_residual_certified': False, 'physical_constraint_status': 'OPEN', 'physical_history_selected': False, 'metric_timestep': False, 'new_stress_or_action_term': False, 'locked_parameter_refit': False, 'publication': False}, 'source_hashes': {p: _constraint_control_digest(p) for p in _constraint_control_SOURCES}, 'input_hashes': {p: _constraint_control_digest(p) for p in (*_constraint_control_RECORDS, *_constraint_control_OWNERS)}, 'reused_artifacts': [records[0]['payload'], drift_split['saved_zero_history_payload']], 'reproducer': 'python3 scripts/derive_nsc_evolved_incoming_constraints.py --check'}

from types import SimpleNamespace as _constraint_control_Namespace
constraint_control = _constraint_control_Namespace(
    EvolvedIncomingState=_constraint_control_EvolvedIncomingState,
    FixedSourcePreparation=_constraint_control_FixedSourcePreparation,
    HistoryBinding=_constraint_control_HistoryBinding,
    MAX=_constraint_control_MAX,
    OUTPUT=_constraint_control_OUTPUT,
    OWNERS=_constraint_control_OWNERS,
    Path=_constraint_control_Path,
    RECORDS=_constraint_control_RECORDS,
    ROOT=_constraint_control_ROOT,
    SOURCES=_constraint_control_SOURCES,
    STATE=_constraint_control_STATE,
    _coverage=_constraint_control__coverage,
    argparse=_constraint_control_argparse,
    authenticate=_constraint_control_authenticate,
    calculate=_constraint_control_calculate,
    coefficients_from_records=_constraint_control_coefficients_from_records,
    column_derivatives=_constraint_control_column_derivatives,
    digest=_constraint_control_digest,
    evolved_constraint_diagnostic=_constraint_control_evolved_constraint_diagnostic,
    json=_constraint_control_json,
    load_inputs=_constraint_control_load_inputs,
    make_interp_spline=_constraint_control_make_interp_spline,
    make_record=_constraint_control_make_record,
    np=_constraint_control_np,
    restore_saved_control=_constraint_control_restore_saved_control,
    sha256=_constraint_control_sha256,
    source_column_matter=_constraint_control_source_column_matter,
    sys=_constraint_control_sys,
)

coefficients_from_records = _constraint_control_coefficients_from_records

from pathlib import Path
# Fixed-preparation control inherited from the local source-response producer.
# It used at most four solves / 300 CPU seconds, 7601 columns, alpha=.001,
# delta-alpha=.0001, 60/120 time steps and a 1e-11 temporal indicator.
from hashlib import sha256
from types import SimpleNamespace


def _prepared_chart_coordinates(rho):
    from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates
    return chart_coordinates(rho)


def _prepared_coefficients_from_records(*records):
    from recursive_horizons.nsc_evolved_incoming_constraints import coefficients_from_records
    return coefficients_from_records(*records)


import json
from recursive_horizons.nsc_posix import ROOT as _prepared_ROOT

_prepared_OUTPUT = 'results/development/nsc-local-prepared-response.json'

_prepared_ARTIFACT = 'results/development/artifacts/nsc-local-prepared-response.npz'

_prepared_OWNED = ('scripts/derive_nsc_cosmological_atlas_v8.py', 'docs/nsc-atlas-prepared-response.md')

_prepared_INPUTS = ('results/development/nsc-local-reference-accuracy.json', 'results/development/nsc-incoming-source-update-v5.json', 'results/development/nsc-incoming-surface-coefficients.json', 'results/development/nsc-incoming-surface-regular-branch.json', 'results/development/nsc-incoming-fixed-transfer.json', 'src/recursive_horizons/nsc_exact_phase_prepared_state.py', 'src/recursive_horizons/nsc_evolved_incoming_state.py', 'src/recursive_horizons/nsc_evolved_incoming_constraints.py', 'src/recursive_horizons/nsc_local_incoming_family.py', 'src/recursive_horizons/nsc_compatible_history_geometry.py', 'src/recursive_horizons/nsc_transmitting_history_jets.py', 'src/recursive_horizons/nsc_transmitting_history_modes.py', 'scripts/derive_nsc_evolved_incoming_constraints.py')

_prepared_ALPHA = 0.001

_prepared_STEP = 0.0001

_prepared_CAP = 300.0

_prepared_CASES = (('center', _prepared_ALPHA, 60, True), ('plus', _prepared_ALPHA + _prepared_STEP, 60, False), ('minus', _prepared_ALPHA - _prepared_STEP, 60, False), ('time_fine', _prepared_ALPHA, 120, True))

def _prepared_digest(path):
    candidate = _receipt_Path(path)
    relative = candidate.relative_to(_prepared_ROOT).as_posix() if candidate.is_absolute() else candidate.as_posix()
    return _receipt_source_sha256(_prepared_ROOT, relative)

def _prepared_signatures():
    return {p: _prepared_digest(p) for p in (*_prepared_OWNED, *_prepared_INPUTS)}

def _prepared_norm(x):
    return float(np.max(abs(x), initial=0.0))

def _prepared_family(alpha):
    from recursive_horizons.nsc_compatible_history_geometry import CompatibleIncomingMetric, CompatibleRadiusDirection
    from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction
    center = _prepared_chart_coordinates(1.0)[1] + 0.15
    w = LocalAxialFunction((0.0, 1.0), center)
    zero = LocalAxialFunction((0.0,), center)
    return CompatibleIncomingMetric((alpha,), (CompatibleRadiusDirection(w, zero, 0.007, 0.03),))

def _prepared_load_inputs():
    records = [json.loads((_prepared_ROOT / p).read_text()) for p in _prepared_INPUTS[:5]]
    p = records[0]['payload']
    if _prepared_digest(p['path']) != p['sha256']:
        raise ValueError('refined initial column artifact changed')
    with np.load(_prepared_ROOT / p['path'], allow_pickle=False) as a:
        arrays = {k: a[k] for k in a.files}
    channel = records[4]['actual_source_bindings']['channel']
    if channel['compact_mass'] != float(arrays['source/mass'].item()) or channel['angular_eigenvalue'] != float(arrays['source/angular'].item()):
        raise ValueError('source/channel operator changed')
    return (records, arrays, channel)

def _prepared_evaluate(arrays, name, alpha, coeff, baseline, parameters):
    from recursive_horizons.nsc_evolved_incoming_constraints import source_column_matter, surface_geometry_response, compatible_history_slots

    def get(key):
        return arrays[name + '/' + key]
    times = get('times')
    z = times + _prepared_chart_coordinates(1.0)[1]
    F, Fz, dF, dFz = (get(k) for k in ('F', 'Fz', 'dF', 'dFz'))
    weights = arrays['source/column_weights']
    C = arrays['source/covariance']
    E = arrays['source/energies_repeated']
    current = source_column_matter(F * weights, Fz * weights, C, dF * weights, dFz * weights, **parameters)
    ref = arrays['reference_F0'][None] * np.exp(-1j * times[:, None, None] * E)
    refz = -1j * E * ref
    zero = np.zeros((0, *ref.shape), complex)
    reference = source_column_matter(ref * weights, refz * weights, C, zero, zero, **parameters)
    slots, ds = compatible_history_slots(_prepared_family(alpha), z, len(dF))
    geom = surface_geometry_response(slots, ds, coeff)
    value = baseline + geom['action_gradient_change'] + current['action_gradient'] - reference['action_gradient']
    tangent = geom['action_gradient_tangent'] + current['action_gradient_tangent']
    return {'residual': value, 'tangent': tangent, 'matter': current['action_gradient'], 'matter_tangent': current['action_gradient_tangent'], 'matter_correction': current['action_gradient'] - reference['action_gradient']}

def _prepared_analyze(arrays, meta):
    rows = {}
    for name, alpha, steps, _ in _prepared_CASES:
        if name not in meta['completed']:
            continue
        got = _prepared_evaluate(arrays, name, alpha, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        mask = np.arange(steps + 1) >= 2 * steps // 3
        rows[name] = {'sampled_local_constraint_approximant_max': np.max(abs(got['residual'][mask]), axis=0).tolist(), 'sampled_local_matter_correction_max': np.max(abs(got['matter_correction'][mask]), axis=0).tolist(), 'phase_residual': meta['cases'][name]['phase_residual'], 'flux_residual': meta['cases'][name]['flux_residual']}
    checks = {}
    if all((n in rows for n in ('center', 'plus', 'minus'))):
        center = _prepared_evaluate(arrays, 'center', _prepared_ALPHA, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        plus = _prepared_evaluate(arrays, 'plus', _prepared_ALPHA + _prepared_STEP, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        minus = _prepared_evaluate(arrays, 'minus', _prepared_ALPHA - _prepared_STEP, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        fd = (plus['residual'] - minus['residual']) / (2 * _prepared_STEP)
        checks['complete_local_constraint_derivative'] = _prepared_norm((fd - center['tangent'][0])[40:])
        checks['local_matter_derivative'] = _prepared_norm(((plus['matter'] - minus['matter']) / (2 * _prepared_STEP) - center['matter_tangent'][0])[40:])
        checks['dropped_deltaC_effect'] = _prepared_norm(center['matter_tangent'][0, 40:])
        checks['same_fixed_preparation'] = len({meta['cases'][n]['preparation'] for n in ('center', 'plus', 'minus')}) == 1
        checks['derivative_passed'] = checks['complete_local_constraint_derivative'] <= 3e-08 and checks['same_fixed_preparation']
    if 'time_fine' in rows and 'center' in rows:
        coarse = _prepared_evaluate(arrays, 'center', _prepared_ALPHA, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        fine = _prepared_evaluate(arrays, 'time_fine', _prepared_ALPHA, meta['coefficients'], np.asarray(meta['baseline']), meta['parameters'])
        checks['local_time_change'] = np.max(abs(fine['matter'][::2][40:] - coarse['matter'][40:]), axis=0).tolist()
        checks['local_tangent_time_change'] = np.max(abs(fine['matter_tangent'][0, ::2][40:] - coarse['matter_tangent'][0, 40:]), axis=0).tolist()
        checks['time_indicator_target'] = 1e-11
        checks['time_indicator_passed'] = max(checks['local_time_change']) <= 1e-11
    return {'per_case': rows, 'checks': checks}

def _prepared_record(arrays, meta, payload):
    measurements = _prepared_analyze(arrays, meta)
    return {'schema': 'NSC-LOCAL-PREPARED-RESPONSE-v1', 'accountable_author': 'Douglas Ek', 'status': 'OPEN: source-fixed nonzero-history control; full-source local gate uncertified', 'control': {'amplitude': _prepared_ALPHA, 'finite_difference_step': _prepared_STEP, 'runs_max': 4, 'profile': 'w=alpha*Chebyshev_T1*flat axial cutoff, U=0; normal .007/.03', 'family': '14_1', 'points': 7601, 'local_interval': 'S(1)+[.12,.18]', 'derivative_tolerance': 3e-08, 'temporal_indicator_tolerance': 1e-11}, 'measurements': measurements, 'runtime': meta['runtime'], 'scope': {'physical_local_gate': 'OPEN', 'full_source_error': None, 'field_and_derivative_error': None, 'between_node_error': None, 'all_channel_coverage': False, 'source_changed': False, 'stress_drift_subtracted': False, 'frozen_C0_imposed': False, 'metric_timestep': False}, 'source_hashes': {p: meta['signature'][p] for p in _prepared_OWNED}, 'input_hashes': {p: meta['signature'][p] for p in _prepared_INPUTS}, 'reused_artifacts': meta['reused_artifacts'], 'payload': payload, 'reproducer': 'call prepared_response.check() from scripts/derive_nsc_cosmological_atlas_v8.py'}

def _prepared_check():
    saved = json.loads((_prepared_ROOT / _prepared_OUTPUT).read_text())
    _prepared_load_inputs()
    if _prepared_digest(_prepared_ARTIFACT) != saved['payload']['sha256']:
        raise ValueError('control payload changed')
    with np.load(_prepared_ROOT / _prepared_ARTIFACT, allow_pickle=False) as a:
        arrays = {k: a[k] for k in a.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    if meta['signature'] != _prepared_signatures():
        raise ValueError('control source/input changed')
    if _prepared_record(arrays, meta, saved['payload']) != saved:
        raise ValueError('array-only replay differs')
    return saved

prepared_response_control = SimpleNamespace(
    ALPHA=_prepared_ALPHA, STEP=_prepared_STEP, CAP=_prepared_CAP, CASES=_prepared_CASES,
    INPUTS=_prepared_INPUTS, OWNED=_prepared_OWNED, OUTPUT=_prepared_OUTPUT,
    ARTIFACT=_prepared_ARTIFACT, family=_prepared_family, load_inputs=_prepared_load_inputs,
    evaluate=_prepared_evaluate, analyze=_prepared_analyze, check=_prepared_check,
    digest=_prepared_digest, coefficients_from_records=_prepared_coefficients_from_records,
    chart_coordinates=_prepared_chart_coordinates)
