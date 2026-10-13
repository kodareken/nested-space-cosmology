"""Numerical insertion errors for one unchanged finite coherent source table."""
from functools import reduce
from flint import acb,arb,ctx
import numpy as np

from .nsc_ks_ball_trajectory import complex_ball,exact_upper,restored_upper


def source_norm_upper(covariance, *, bits=120):
    """sqrt(||C||_1 ||C||_infinity), retaining all source coherences."""
    matrix=np.asarray(covariance,complex)
    if matrix.ndim!=2 or matrix.shape[0]!=matrix.shape[1] or not matrix.shape[0] or not np.isfinite(matrix).all():
        raise ValueError('finite square coherent source matrix required')
    with ctx.workprec(bits):
        magnitude=[[complex_ball(v).abs_upper() for v in row] for row in matrix]
        row=reduce(arb.max,(sum(values,arb(0)).upper() for values in magnitude))
        column=reduce(arb.max,(sum((magnitude[i][j] for i in range(len(matrix))),arb(0)).upper() for j in range(len(matrix))))
        return (row*column).sqrt().upper()


def endpoint_norms(field,energies):
    """Uniform weighted Frobenius norms of the reconstructed F and actual F_z."""
    if len(energies)!=field.source_count:raise ValueError('matching source energies required')
    with ctx.workprec(field.bits):
        rows=[arb(0),arb(0)];derivative=[arb(0),arb(0)]
        for spin in range(2):
            for k,wave in enumerate(field.waves):
                values=[sum((field.coefficients[p][spin][s][k] for p in range(8)),acb(0)) for s in range(field.source_count)]
                rows[spin]+=sum((v.abs_upper()**2 for v in values),arb(0)).upper().sqrt()
                derivative[spin]+=sum(((acb(0,wave-arb(float(E)))*v).abs_upper()**2
                                      for E,v in zip(energies,values)),arb(0)).upper().sqrt()
        return (sum((v*v for v in rows),arb(0)).sqrt().upper(),
                sum((v*v for v in derivative),arb(0)).sqrt().upper())


def finite_matter_error(field_norm,axial_norm,field_error,axial_error,source_norm,*,
                        mass,absolute_angular,axial_lower,radius_lower,multiplicity,bits=120):
    """Same archived C for exact/reconstructed columns; source error is separate."""
    with ctx.workprec(bits):
        values=list(map(arb,(field_norm,axial_norm,field_error,axial_error,source_norm,
                             mass,absolute_angular,axial_lower,radius_lower,multiplicity)))
        if any(not v.is_finite() or not v>=0 for v in values) or any(not v>0 for v in values[7:]):
            raise ValueError('finite nonnegative bounds and positive geometry/multiplicity required')
        f,z,e,ez,C,m,ell,a,r,mu=values
        density=C*(2*f*e+e*e)
        current=C*(z*e+f*ez+e*ez)
        potential=(m*m+(ell/r)**2).sqrt()
        return {'N':exact_upper(mu*(potential*density+current/a)),
                'beta':exact_upper(mu*current),
                'source_accuracy_included':False,'physical_local_gate':'OPEN'}

from pathlib import Path as _receipt_Path
from .nsc_frozen_sources import selected_source_sha256 as _receipt_source_sha256

# Retained executable control: fine_control
# Folded from scripts/derive_nsc_ks_fine_matter_accuracy.py at 329ef809; original receipts retain their producer identities.
import os as _fine_control_os
import argparse as _fine_control_argparse
from hashlib import sha256 as _fine_control_sha256
import json as _fine_control_json
from pathlib import Path as _fine_control_Path
import sys as _fine_control_sys
import numpy as _fine_control_np
from flint import arb as _fine_control_arb, acb as _fine_control_acb, acb_mat as _fine_control_acb_mat, ctx as _fine_control_ctx
from recursive_horizons.nsc_posix import ROOT as _fine_control_ROOT
_fine_control_sys.path.insert(0, str(_fine_control_ROOT / 'src'))
_fine_control_sys.path.insert(0, str(_fine_control_ROOT / 'scripts'))
import validate_nsc_ks_fine_history as _fine_control_V
import derive_nsc_ks_fine_trajectory as _fine_control_F
from recursive_horizons import nsc_ks_source_envelope as _fine_control_E
from recursive_horizons.nsc_ks_difference_envelope import KSDifferenceIncoming as _fine_control_KSDifferenceIncoming
from recursive_horizons.nsc_ks_ball_trajectory import BallFourierSegment as _fine_control_BallFourierSegment, complex_ball as _fine_control_complex_ball, exact_upper as _fine_control_exact_upper, restored_upper as _fine_control_restored_upper
from recursive_horizons.nsc_ks_reference_error import spectator_reference_error as _fine_control_spectator_reference_error, homogeneous_reference_norms as _fine_control_homogeneous_reference_norms
from recursive_horizons.nsc_ks_finite_matter_error import source_norm_upper as _fine_control_source_norm_upper, endpoint_norms as _fine_control_endpoint_norms, finite_matter_error as _fine_control_finite_matter_error
from recursive_horizons.nsc_ks_local_constraints import homogeneous_reference_amplitudes as _fine_control_homogeneous_reference_amplitudes, stationary_computational_reference as _fine_control_stationary_computational_reference
from recursive_horizons.nsc_ks_matter_difference import source_fixed_matter_difference as _fine_control_source_fixed_matter_difference
from recursive_horizons.nsc_ks_streaming_trajectory import TrajectoryChunkReader as _fine_control_TrajectoryChunkReader
_fine_control_OUTPUT = _fine_control_ROOT / 'results/development/nsc-ks-fine-matter-accuracy.json'
_fine_control_OWNED = ('scripts/derive_nsc_ks_fine_matter_accuracy.py', 'docs/nsc-ks-fine-matter-accuracy.md', 'src/recursive_horizons/nsc_ks_reference_error.py', 'src/recursive_horizons/nsc_ks_finite_matter_error.py', 'src/recursive_horizons/nsc_ks_local_constraints.py', 'src/recursive_horizons/nsc_ks_matter_difference.py', 'src/recursive_horizons/nsc_evolved_incoming_constraints.py', 'src/recursive_horizons/nsc_common_subtracted_ks_source.py')
_fine_control_BITS = 120
_fine_control_CONTROL_TOLERANCE = 1e-11

def _fine_control_digest(path):
    candidate = _receipt_Path(path)
    relative = candidate.relative_to(_fine_control_ROOT).as_posix() if candidate.is_absolute() else candidate.as_posix()
    return _receipt_source_sha256(_fine_control_ROOT, relative)

def _fine_control_restore_state(a, meta, family):
    grid = a['computational_z']
    w, U = _fine_control_E._sample_axial_profiles(family.directions, grid)
    binding = _fine_control_E.KSEnvelopeBinding(grid, w, U, family.amplitudes, tuple((d.inner_radius for d in family.directions)), tuple((d.outer_radius for d in family.directions)), _fine_control_E.usual_axial_support(), meta['rho_up'], 1.0, **meta['options'])
    if binding.fingerprint != meta['history_fingerprint']:
        raise ValueError('fine history binding changed')
    get = lambda key: a['prepared/' + key]
    return _fine_control_KSDifferenceIncoming(a['target_z'], get('columns'), get('column_tangents'), get('axial_columns'), get('axial_tangents'), a['source/covariance'], a['source/column_weights'], a['source/energies'], meta['parameters']['mass'], meta['parameters']['angular'], a['initial_canonical_columns'], meta['rho_up'], binding, meta['preparation'], None, meta['diagnostics'], get('reference_amplitudes'), get('envelope_difference'), get('envelope_difference_z'))

def _fine_control_matrix(array):
    return _fine_control_acb_mat([[_fine_control_complex_ball(v) for v in row] for row in array])

def _fine_control_support_hull(supports):
    if not supports:
        raise ValueError('at least one axial support required')
    low, high = supports[0]
    for left, right in supports[1:]:
        low = _fine_control_arb.min(low, left)
        high = _fine_control_arb.max(high, right)
    return (low, high)

def _fine_control_insertion(Fcol, Fzcol, C, mass, angular, multiplicity):
    D = Fcol * C * Fcol.conjugate().transpose()
    J = (Fzcol * C * Fcol.conjugate().transpose() - Fcol * C * Fzcol.conjugate().transpose()) / _fine_control_acb(0, 2)
    s1 = _fine_control_acb_mat([[0, 1], [1, 0]])
    s2 = _fine_control_acb_mat([[0, _fine_control_acb(0, -1)], [_fine_control_acb(0, 1), 0]])
    s3 = _fine_control_acb_mat([[1, 0], [0, -1]])
    radius = _fine_control_arb(2).sqrt()
    axial = (3 * _fine_control_arb.pi() / 2 - 4).sqrt()
    return (multiplicity * (mass * (s1 * D).trace() - angular / radius * (s2 * D).trace() - (s3 * J).trace() / axial), multiplicity * J.trace())

def _fine_control_node_errors(field, state, reference, reference_z, cache, origin):
    maximum = [_fine_control_arb(0) for _ in range(4)]
    for node, z in enumerate(state.z):
        X = field.value(1, _fine_control_arb(float(z)) - origin)
        Xz = field.value(1, _fine_control_arb(float(z)) - origin, axial_derivative=1)
        errors = [_fine_control_arb(0) for _ in range(4)]
        for spin in range(2):
            for source, (energy, weight) in enumerate(zip(state.source_energies, state.column_weights)):
                phase = _fine_control_acb(0, -_fine_control_arb(float(energy)) * _fine_control_arb(float(z))).exp()
                current = phase * X[spin][source]
                current_z = phase * (Xz[spin][source] - _fine_control_acb(0, _fine_control_arb(float(energy))) * X[spin][source])
                ref = phase * _fine_control_complex_ball(cache[spin, source]) * _fine_control_arb(float(weight))
                refz = -_fine_control_acb(0, _fine_control_arb(float(energy))) * ref
                exact = (current, current_z, ref, refz)
                observed = (state.weighted_columns[node, spin, source], state.weighted_axial_columns[node, spin, source], reference[node, spin, source] * weight, reference_z[node, spin, source] * weight)
                for i, (x, y) in enumerate(zip(exact, observed)):
                    errors[i] += (x - _fine_control_complex_ball(y)).abs_upper() ** 2
        for i, error in enumerate(errors):
            maximum[i] = _fine_control_arb.max(maximum[i], error.upper().sqrt().upper())
    return tuple((v.upper() for v in maximum))

def _fine_control_compute():
    if not _fine_control_V.OUTPUT.exists():
        raise ValueError('all history cells must be validated first')
    validation = _fine_control_json.loads(_fine_control_V.OUTPUT.read_text())
    if not validation['all_segments_covered'] or validation['numerical_field_error_relative_to_saved_upstream'] is None:
        raise ValueError('partial history bounds cannot certify matter accuracy')
    for path, expected in validation['signature'].items():
        if _fine_control_digest(path) != expected:
            raise ValueError('validated history code/input changed')
    for item in validation['segment_proofs']:
        if _fine_control_digest(item['path']) != item['sha256']:
            raise ValueError('segment proof changed')
    fine = _fine_control_json.loads(_fine_control_F.OUTPUT.read_text())
    if _fine_control_digest(fine['payload']['path']) != fine['payload']['sha256']:
        raise ValueError('fine prepared fields changed')
    with _fine_control_np.load(_fine_control_F.FINAL, allow_pickle=False) as f:
        a = {k: f[k] for k in f.files}
    meta = _fine_control_json.loads(a['metadata_json'].tobytes())
    family = _fine_control_F.T.C.R.C.P.family(_fine_control_F.T.C.R.C.P.ALPHA)
    state = _fine_control_restore_state(a, meta, family)
    reader = _fine_control_TrajectoryChunkReader(_fine_control_F.DIRECTORY)
    last = reader.load_chunk(len(reader) - 1)
    segment = last.segment(len(last) - 1)
    cached = _fine_control_homogeneous_reference_amplitudes(state)
    reference, reference_z = _fine_control_stationary_computational_reference(state)
    params = meta['parameters']
    numerical = _fine_control_source_fixed_matter_difference(state, axial_scale=params['axial_scale'], radius=params['radius'], multiplicity=params['multiplicity'], reference='cached')
    with _fine_control_ctx.workprec(_fine_control_BITS):
        field = _fine_control_BallFourierSegment(segment, state.column_weights, len(a['computational_z']), meta['period_length'], bits=_fine_control_BITS)
        norm, z_norm = _fine_control_endpoint_norms(field, state.source_energies)
        field_error = validation['numerical_field_error_relative_to_saved_upstream']
        epsF = _fine_control_restored_upper(field_error['weighted_F_error_upper'])
        epsZ = _fine_control_restored_upper(field_error['weighted_F_z_error_upper'])
        size = state.initial_columns.size
        Dend = segment.end[size:].reshape(2, len(state.source_energies), len(a['computational_z']))
        supports = [(_fine_control_arb(d.w.center) - _fine_control_arb(d.w.outer), _fine_control_arb(d.w.center) + _fine_control_arb(d.w.outer)) for d in family.directions]
        low, high = _fine_control_support_hull(supports)
        origin = _fine_control_arb(float(a['computational_z'][0]))
        ref_error = _fine_control_spectator_reference_error(state.rho_up, origin, field.length, low, high, Dend[:, :, 0], state.column_weights, state.source_energies, _fine_control_restored_upper(field_error['row_field_error_upper']), reference_pair=(state.reference_amplitudes, cached), bits=_fine_control_BITS)
        eval_errors = _fine_control_node_errors(field, state, reference, reference_z, cached, origin)
        refF, refZ = _fine_control_homogeneous_reference_norms(cached, state.column_weights, state.source_energies)
        gamma = _fine_control_source_norm_upper(state.source_covariance)
        parameters = dict(mass=state.mass, absolute_angular=abs(state.angular), axial_lower=_fine_control_arb(4) / 5, radius_lower=_fine_control_arb(7) / 5, multiplicity=params['multiplicity'], bits=_fine_control_BITS)
        current = _fine_control_finite_matter_error(norm + eval_errors[0], z_norm + eval_errors[1], epsF + eval_errors[0], epsZ + eval_errors[1], gamma, **parameters)
        ref = _fine_control_finite_matter_error(refF + eval_errors[2], refZ + eval_errors[3], _fine_control_restored_upper(ref_error['weighted_reference_F_error_upper']) + eval_errors[2], _fine_control_restored_upper(ref_error['weighted_reference_F_z_error_upper']) + eval_errors[3], gamma, **parameters)
        C = _fine_control_matrix(state.source_covariance)
        arithmetic = [_fine_control_arb(0), _fine_control_arb(0)]
        for node in range(len(state.z)):
            cur = _fine_control_insertion(_fine_control_matrix(state.weighted_columns[node]), _fine_control_matrix(state.weighted_axial_columns[node]), C, _fine_control_arb(state.mass), _fine_control_arb(state.angular), _fine_control_arb(params['multiplicity']))
            old = _fine_control_insertion(_fine_control_matrix(reference[node] * state.column_weights), _fine_control_matrix(reference_z[node] * state.column_weights), C, _fine_control_arb(state.mass), _fine_control_arb(state.angular), _fine_control_arb(params['multiplicity']))
            for j in range(2):
                difference = cur[j] - old[j] - _fine_control_arb(float(numerical['action_gradient_change'][node, j]))
                arithmetic[j] = _fine_control_arb.max(arithmetic[j], difference.abs_upper())
        total = [(_fine_control_restored_upper(current[key]) + _fine_control_restored_upper(ref[key]) + arithmetic[j]).upper() for j, key in enumerate(('N', 'beta'))]
        passed = all((v <= _fine_control_arb(_fine_control_CONTROL_TOLERANCE) for v in total))
        result = {'schema': 'NSC-KS-FINE-MATTER-ACCURACY-v1', 'accountable_author': 'Douglas Ek', 'status': ('PASS' if passed else 'OPEN') + ': finite-source numerical matter error; physical gate OPEN', 'numerical_control_pass': passed, 'numerical_control_tolerance': _fine_control_CONTROL_TOLERANCE, 'current_field_error': current, 'reference_field_error': ref, 'reference_error_proof': ref_error, 'node_evaluation_error_upper': [_fine_control_exact_upper(v) for v in eval_errors], 'arithmetic_and_intrinsic_geometry_error_upper': [_fine_control_exact_upper(v) for v in arithmetic], 'total_N_beta_error_upper': [_fine_control_exact_upper(v) for v in total], 'source_matrix_operator_norm_upper': _fine_control_exact_upper(gamma), 'scope': {'source_family': '14_1', 'sampled_incoming_nodes': len(state.z), 'full_source_accuracy': None, 'source_energy_UV_bound': None, 'constraint_root_claimed': False, 'physical_local_gate': 'OPEN', 'metric_timestep': False, 'stress_drift_subtracted': False}, 'source_hashes': {p: _fine_control_digest(p) for p in _fine_control_OWNED}, 'input_hashes': {str(_fine_control_V.OUTPUT.relative_to(_fine_control_ROOT)): _fine_control_digest(_fine_control_V.OUTPUT), str(_fine_control_F.OUTPUT.relative_to(_fine_control_ROOT)): _fine_control_digest(_fine_control_F.OUTPUT)}, 'reproducer': 'python scripts/derive_nsc_ks_fine_matter_accuracy.py --check'}
    return result

from types import SimpleNamespace as _fine_control_Namespace
fine_control = _fine_control_Namespace(
    BITS=_fine_control_BITS,
    BallFourierSegment=_fine_control_BallFourierSegment,
    CONTROL_TOLERANCE=_fine_control_CONTROL_TOLERANCE,
    E=_fine_control_E,
    F=_fine_control_F,
    KSDifferenceIncoming=_fine_control_KSDifferenceIncoming,
    OUTPUT=_fine_control_OUTPUT,
    OWNED=_fine_control_OWNED,
    Path=_fine_control_Path,
    ROOT=_fine_control_ROOT,
    TrajectoryChunkReader=_fine_control_TrajectoryChunkReader,
    V=_fine_control_V,
    acb=_fine_control_acb,
    acb_mat=_fine_control_acb_mat,
    arb=_fine_control_arb,
    argparse=_fine_control_argparse,
    complex_ball=_fine_control_complex_ball,
    compute=_fine_control_compute,
    ctx=_fine_control_ctx,
    digest=_fine_control_digest,
    endpoint_norms=_fine_control_endpoint_norms,
    exact_upper=_fine_control_exact_upper,
    finite_matter_error=_fine_control_finite_matter_error,
    homogeneous_reference_amplitudes=_fine_control_homogeneous_reference_amplitudes,
    homogeneous_reference_norms=_fine_control_homogeneous_reference_norms,
    insertion=_fine_control_insertion,
    json=_fine_control_json,
    matrix=_fine_control_matrix,
    node_errors=_fine_control_node_errors,
    np=_fine_control_np,
    os=_fine_control_os,
    restore_state=_fine_control_restore_state,
    restored_upper=_fine_control_restored_upper,
    sha256=_fine_control_sha256,
    source_fixed_matter_difference=_fine_control_source_fixed_matter_difference,
    source_norm_upper=_fine_control_source_norm_upper,
    spectator_reference_error=_fine_control_spectator_reference_error,
    stationary_computational_reference=_fine_control_stationary_computational_reference,
    support_hull=_fine_control_support_hull,
    sys=_fine_control_sys,
)
