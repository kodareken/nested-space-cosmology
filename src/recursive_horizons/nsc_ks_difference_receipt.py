"""Explicit historical difference-state replay, with current numerical bytes.

Authentication precedes construction. Source pins select read-only Git
bytes; JSON/NPZ witnesses always come from the current root. Original and
current profile fingerprints are authenticated separately. The phase reconstruction
allowance is a finite arithmetic/library diagnostic, not a field error,
integration error or changed scientific acceptance criterion.
"""
from dataclasses import fields
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import re

import numpy as np

from . import nsc_ks_source_envelope as E
from .nsc_frozen_sources import authenticate_source, is_implementation_source, source_sha256
from .nsc_receipt_metadata import authenticate_rebound_receipt


def _gamma(n):
    eps = np.finfo(float).eps
    return n*eps/(1-n*eps)


def _up(value):
    return np.nextafter(np.asarray(value,float)*(1+_gamma(16)),np.inf)


def _numeric(value, *, complex_allowed):
    raw = np.asarray(value)
    allowed = 'fiu'+('c' if complex_allowed else '')
    if raw.dtype.kind not in allowed or not np.isfinite(raw).all():
        raise ValueError('finite numeric phase operands required')
    return raw


def phase_reconstruction_diagnostic(z, energies, A, D, Dz, F, Fz):
    """Two executions of exp(-iEz)*(A+D) and its axial derivative.

    gamma_2*|Ez| charges the phase argument. Complex addition/product chains
    are bounded by uncancelled component 1-norms (gamma_16 for F and
    gamma_32 for Fz). Both replays are charged. Separately, four ulps per
    unit sine/cosine component yield an 8*eps complex phase indicator per
    execution; it is an explicit library-accuracy model, not a proof of
    a physical error. Zero operand mass receives no absolute floor.
    """
    z,energies = (_numeric(v,complex_allowed=False) for v in (z,energies))
    A,D,Dz,F,Fz = (_numeric(v,complex_allowed=True) for v in (A,D,Dz,F,Fz))
    shape = (len(z),2,len(energies))
    if (z.ndim!=1 or energies.ndim!=1 or A.shape!=(2,len(energies))
            or any(v.shape!=shape for v in (D,Dz,F,Fz))):
        raise ValueError('phase operands must retain target, spin and source axes')
    norm1 = lambda value: abs(np.real(value))+abs(np.imag(value))
    mass = norm1(A)[None]+norm1(D)
    axial_mass = norm1(Dz)+abs(energies)[None,None]*mass
    argument = abs(energies[None]*z[:,None])[:,None]
    eps = np.finfo(float).eps
    arithmetic_F = _up(2*(_gamma(2)*argument+2*_gamma(16))*mass*(1+_gamma(32)))
    arithmetic_Fz = _up(2*(_gamma(2)*argument+2*_gamma(32))*axial_mass*(1+_gamma(32)))
    library_F = _up(16*eps*mass*(1+_gamma(32)))
    library_Fz = _up(16*eps*axial_mass*(1+_gamma(32)))
    # Exact zero inputs reconstruct exact zero in both represented chains.
    arithmetic_F = np.where(mass,arithmetic_F,0.)
    library_F = np.where(mass,library_F,0.)
    arithmetic_Fz = np.where(axial_mass,arithmetic_Fz,0.)
    library_Fz = np.where(axial_mass,library_Fz,0.)
    phase = np.exp(-1j*energies*z[:,None])[:,None]
    total = A[None]+D
    reconstructed_F = phase*total
    reconstructed_Fz = phase*(Dz-1j*energies*total)
    difference_F = norm1(F-reconstructed_F)
    difference_Fz = norm1(Fz-reconstructed_Fz)
    limit_F = _up(arithmetic_F+library_F)
    limit_Fz = _up(arithmetic_Fz+library_Fz)
    limit_F = np.where(mass,limit_F,0.)
    limit_Fz = np.where(axial_mass,limit_Fz,0.)
    if any(not np.isfinite(value).all() for value in (
            mass,axial_mass,argument,reconstructed_F,reconstructed_Fz,
            difference_F,difference_Fz,arithmetic_F,arithmetic_Fz,
            library_F,library_Fz,limit_F,limit_Fz)):
        raise ValueError('nonfinite phase replay arithmetic or diagnostic')
    if np.any(difference_F>limit_F) or np.any(difference_Fz>limit_Fz):
        raise ValueError('historical phase reconstruction exceeds finite replay diagnostic')
    return {'maximum_F_difference':float(np.max(difference_F)),
            'maximum_Fz_difference':float(np.max(difference_Fz)),
            'maximum_F_arithmetic_bound':float(np.max(arithmetic_F)),
            'maximum_Fz_arithmetic_bound':float(np.max(arithmetic_Fz)),
            'maximum_F_library_accuracy_indicator':float(np.max(library_F)),
            'maximum_Fz_library_accuracy_indicator':float(np.max(library_Fz)),
            'library_model':'four ulps per unit sine/cosine component; not a proved physical error',
            'physical_field_error_bound':None,'integration_error_bound':None,
            'physical_constraint_status':'OPEN'}


def authenticated_difference_receipt(cls, root, receipt_path, expected_receipt_sha256,
                                     *, case, binding, family, source_context, source_contexts=None):
    if not isinstance(source_context,str) or not re.fullmatch('[0-9a-f]{40}',source_context):
        raise ValueError('explicit full original source context required')
    if not isinstance(binding,E.KSEnvelopeBinding):
        raise TypeError('current KSEnvelopeBinding required')
    root = Path(root).resolve()
    raw = authenticate_source(root,receipt_path,expected_receipt_sha256,context=None)
    record = json.loads(raw)
    if record.get('schema')!='NSC-KS-DIFFERENCE-CONTROL-v1':
        raise ValueError('original difference control receipt required')
    spec = record['payload']
    payload = authenticate_source(root,spec['path'],spec['sha256'],context=None)
    if type(spec.get('bytes')) is not int or spec['bytes']!=len(payload):
        raise ValueError('original difference payload byte count changed')
    with np.load(BytesIO(payload),allow_pickle=False) as loaded:
        arrays = {key:loaded[key].copy() for key in loaded.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    signature = meta['signature']
    if {**record['source_hashes'],**record['input_hashes']}!=signature:
        raise ValueError('original difference source/input receipt changed')
    overrides = {} if source_contexts is None else dict(source_contexts)
    if any(name not in signature or not is_implementation_source(name) for name in overrides):
        raise ValueError('source override must identify a declared implementation source')
    authenticated_contexts = {}
    for name,expected in signature.items():
        historical = is_implementation_source(name)
        context = overrides.get(name,source_context) if historical else None
        if not historical and source_sha256(root,name)!=expected:
            # Preserve this capture's identity while proving the current local
            # carrier differs only by the authorized receipt metadata update.
            authenticate_rebound_receipt(root,name,expected)
        else:
            authenticate_source(root,name,expected,context=context)
        if historical:
            authenticated_contexts[name]=context
    required = {'scripts/derive_nsc_ks_difference_control.py',
                'src/recursive_horizons/nsc_ks_difference_envelope.py',
                'src/recursive_horizons/nsc_ks_source_envelope.py'}
    if not required<=authenticated_contexts.keys():
        raise ValueError('original difference producer source pins required')
    cases = [row for row in record['control']['cases'] if row[0]==case]
    if len(cases)!=1 or case not in meta['completed'] or case not in meta['cases']:
        raise ValueError('original completed difference case required')
    _,points,options = cases[0]
    if (type(points) is not int or points!=len(binding.computational_z)
            or any(getattr(binding,key)!=value for key,value in options.items())
            or binding.rho_up!=meta['rho_up'] or binding.rho_sigma!=1.
            or binding.amplitudes!=(record['control']['history_amplitude'],)
            or binding.inner_radii!=(.007,) or binding.outer_radii!=(.03,)
            or binding.axial_support!=E.usual_axial_support()
            or not np.array_equal(binding.computational_z,E.computational_z_grid(points))):
        raise ValueError('current binding differs from original solver/domain recipe')
    from .nsc_compatible_history_geometry import CompatibleIncomingMetric
    if not isinstance(family,CompatibleIncomingMetric) or len(family.directions)!=1:
        raise ValueError('original one-direction linear profile recipe required')
    direction = family.directions[0]
    if (tuple(family.amplitudes)!=binding.amplitudes
            or direction.w.coefficients!=(0.,1.) or direction.U.coefficients!=(0.,)
            or direction.w.inner!=.03 or direction.w.outer!=.06
            or direction.U.inner!=.03 or direction.U.outer!=.06
            or direction.w.center!=E.chart_coordinates(1.)[1]+.15
            or direction.U.center!=direction.w.center
            or direction.inner_radius!=.007 or direction.outer_radius!=.03):
        raise ValueError('original one-direction linear profile recipe changed')
    expected_w,expected_U = E._sample_axial_profiles(family.directions,binding.computational_z)
    if not np.array_equal(binding.w_samples,expected_w) or not np.array_equal(binding.U_samples,expected_U):
        raise ValueError('current binding samples differ from the declared profile')
    parameters = meta['parameters']
    values = dict(z=arrays['target_z'],columns=arrays[case+'/F'],
                  column_tangents=arrays[case+'/dF'],axial_columns=arrays[case+'/Fz'],
                  axial_tangents=arrays[case+'/dFz'],source_covariance=arrays['source/covariance'],
                  column_weights=arrays['source/column_weights'],source_energies=arrays['source/energies'],
                  mass=parameters['mass'],angular=parameters['angular'],
                  initial_columns=arrays['initial_canonical_columns'],rho_up=meta['rho_up'],
                  binding=binding,fixed_preparation_digest=meta['cases'][case]['preparation'],
                  coverage=None,diagnostics={},reference_amplitudes=arrays[case+'/A'],
                  envelope_difference=arrays[case+'/D'],envelope_difference_z=arrays[case+'/Dz'])
    if set(values)!={field.name for field in fields(cls)}:
        raise ValueError('difference-state storage contract changed')
    # Bypass only the default constructor's bit-equality reconstruction.
    # The complete base storage/source/domain predicates are still invoked.
    state = object.__new__(cls)
    for name,value in values.items():
        object.__setattr__(state,name,value)
    E.KSEnvelopeIncoming.__post_init__(state)
    for name,shape in (('reference_amplitudes',state.initial_columns.shape),
                       ('envelope_difference',state.columns.shape),
                       ('envelope_difference_z',state.columns.shape)):
        value = E._finite_complex(getattr(state,name),name)
        if value.shape!=shape:
            raise ValueError(name+' has the wrong shape')
        object.__setattr__(state,name,value)
    diagnostic = phase_reconstruction_diagnostic(
        state.z,state.source_energies,state.reference_amplitudes,state.envelope_difference,
        state.envelope_difference_z,state.columns,state.axial_columns)
    from types import MappingProxyType
    object.__setattr__(state,'diagnostics',MappingProxyType({
        'historical_receipt_path':str(receipt_path),
        'historical_receipt_sha256':expected_receipt_sha256,
        'historical_payload_sha256':sha256(payload).hexdigest(),
        'original_source_contexts':MappingProxyType(authenticated_contexts),
        'current_source_hashes':MappingProxyType({path:source_sha256(root,path) for path in (
            'src/recursive_horizons/nsc_ks_difference_receipt.py',
            'src/recursive_horizons/nsc_ks_difference_envelope.py',
            'src/recursive_horizons/nsc_ks_source_envelope.py',
            'src/recursive_horizons/nsc_local_incoming_family.py')}),
        'original_binding_fingerprint':meta['cases'][case]['binding'],
        'current_binding_fingerprint':binding.fingerprint,
        'original_profile_samples_retained':False,
        'phase_replay_diagnostic':MappingProxyType(diagnostic)}))
    state.validate()
    return state
