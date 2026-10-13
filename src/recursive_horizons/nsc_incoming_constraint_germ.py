"""Structural two-direction relation for the retained incoming constraints.

This owns a conditional local-germ relation, not selected initial data,
a source-accuracy certificate, a Cauchy-surface solution or metric evolution.
"""
import sympy as sp


def structural_constraint_plane():
    """Grade and reflection of the EXISTING four-derivative bulk owners."""
    u,v=sp.symbols('u v',real=True)
    allowed=[(i,j) for i in range(2) for j in range(3) if 3*i+2*j<=4]
    even=[pair for pair in allowed if pair[1]%2==0]
    odd=[pair for pair in allowed if pair[1]%2==1]
    if set(even)!={(0,0),(1,0),(0,2)} or odd!=[(0,1)]:
        raise ArithmeticError('retained derivative grade produced unexpected monomials')
    N,beta,a,r,m,ell,k=sp.symbols('N beta a r m ell k',real=True,nonzero=True)
    s1=sp.Matrix([[0,1],[1,0]]);s2=sp.Matrix([[0,-sp.I],[sp.I,0]]);s3=sp.diag(1,-1);I=sp.eye(2)
    H=-N*m*s1+N*ell/r*s2+N*k/a*s3-beta*k*I
    reflected=H.subs({ell:-ell,k:-k,beta:-beta},simultaneous=True)
    residual=sp.simplify(reflected-s1*H*s1)
    vertexN=H.diff(N);vertexBeta=H.diff(beta)
    reflectedN=vertexN.subs({ell:-ell,k:-k,beta:-beta},simultaneous=True)
    reflectedBeta=vertexBeta.subs({ell:-ell,k:-k,beta:-beta},simultaneous=True)
    vertex_residuals=[sp.simplify(reflectedN-s1*vertexN*s1),
                      sp.simplify(reflectedBeta+s1*vertexBeta*s1)]
    if residual!=sp.zeros(2) or any(x!=sp.zeros(2) for x in vertex_residuals):
        raise ArithmeticError('Hamiltonian/vertex reflection identity failed')
    SN,Sb,cvv=sp.symbols('S_N S_beta c_v2',real=True)
    cu,cv=sp.symbols('c_u c_v',real=True,nonzero=True)
    EN=SN+cu*u+cvv*v*v;Eb=Sb+cv*v
    formal_v=-Sb/cv;formal_u=-(SN+cvv*formal_v**2)/cu
    substitution={u:formal_u,v:formal_v}
    identities=[sp.simplify(EN.subs(substitution)),sp.simplify(Eb.subs(substitution))]
    if identities!=[0,0]:raise ArithmeticError('conditional algebraic relation failed')
    return {'variables':{'u':'delta r_TTT','v':'delta r_Tz'},'derivative_grades':{'u':3,'v':2},
            'retained_maximum_grade':4,'allowed_lapse_monomials':[list(x) for x in even],
            'allowed_shift_response_monomials':[list(x) for x in odd],
            'constraint_form':{'N':'S_N+c_u*u+c_v2*v^2','beta':'S_beta+c_v*v'},
            'reflection':{'coordinates':'z -> -z; k -> -k; ell -> -ell; beta -> -beta',
                          'matrix':'sigma1','Hamiltonian_residual':'0','vertex_residuals':['0','0'],
                          'state_reflection_required':False,'state_enters_fixed_constants':True},
            'conditional_relation':{'v':'-S_beta/c_v','u':'-(S_N+c_v2*v^2)/c_u',
                                    'requires':['c_u != 0','c_v != 0','finite retained source and reference insertion'],
                                    'substitution_residuals':[str(x) for x in identities]},
            'coefficients_are_action_responses_not_adjustable_couplings':True,
            'physical_initial_data_or_surface_solution':False}

# Retained executable control: germ_control
# Actual API folded from scripts/derive_nsc_incoming_constraint_germ.py at 329ef809; historical source identities remain frozen.
from .nsc_frozen_sources import selected_source_sha256 as _retained_incoming_selected_sha256
def _retained_incoming_source_sha256(root, path):
    from pathlib import Path
    value = Path(path)
    name = value.relative_to(root).as_posix() if value.is_absolute() else path
    return _retained_incoming_selected_sha256(root, name)

import argparse as _germ_control_argparse
import hashlib as _germ_control_hashlib
import json as _germ_control_json
from pathlib import Path as _germ_control_Path
import sys as _germ_control_sys
import numpy as _germ_control_np
from recursive_horizons.nsc_posix import ROOT as _germ_control_ROOT
_germ_control_sys.path.insert(0, str(_germ_control_ROOT / 'src'))
from recursive_horizons.nsc_incoming_constraint_germ import structural_constraint_plane as _germ_control_structural_constraint_plane
_germ_control_OUTPUT = _germ_control_ROOT / 'results/development/nsc-incoming-constraint-germ.json'
_germ_control_SOURCES = ('src/recursive_horizons/nsc_incoming_constraint_germ.py', 'scripts/derive_nsc_incoming_constraint_germ.py', 'tests/test_nsc_incoming_constraint_germ.py', 'docs/nsc-incoming-constraint-germ.md')
_germ_control_INPUTS = ('results/development/nsc-incoming-joint-constraints.json', 'results/development/nsc-incoming-local-constraints.json', 'results/development/nsc-reference-band-bulk.json', 'src/recursive_horizons/nsc_incoming_reference_response.py', 'src/recursive_horizons/nsc_incoming_local_constraints.py', 'src/recursive_horizons/nsc_incoming_cauchy_jets.py', 'src/recursive_horizons/nsc_light_restoration_action.py', 'src/recursive_horizons/nsc_spherical_local_history.py', 'src/recursive_horizons/nsc_spatial_reference_symbol.py', 'src/recursive_horizons/nsc_reference_band_action.py')

def _germ_control_sha(p):
    return _retained_incoming_source_sha256(_germ_control_ROOT, p)

def _germ_control_read(p):
    return _germ_control_json.loads((_germ_control_ROOT / p).read_text())

def _germ_control_signature():
    return {p: _germ_control_sha(p) for p in (*_germ_control_SOURCES, *_germ_control_INPUTS)}

def _germ_control_authenticate():
    for p in _germ_control_INPUTS[:3]:
        data = _germ_control_read(p)
        for key in ('source_hashes', 'input_hashes'):
            for path, digest in data.get(key, {}).items():
                if _germ_control_sha(path) != digest:
                    raise ValueError('response input changed: ' + path)
        if 'payload' in data and _germ_control_sha(data['payload']['path']) != data['payload']['sha256']:
            raise ValueError('response input artifact changed')

def _germ_control_prepare(directory):
    _germ_control_authenticate()
    controls = {'u': _germ_control_json.loads((directory / 'nsc-incoming-independent-lapse-response.json').read_text()), 'mixed': _germ_control_json.loads((directory / 'nsc-incoming-mixed-germ-response.json').read_text())}
    for control in controls.values():
        for p, digest in control['producer_hashes'].items():
            if _germ_control_sha(p) != digest:
                raise ValueError('response producer changed: ' + p)
    raw = (_germ_control_json.dumps({'signature': _germ_control_signature(), 'controls': controls}, indent=2, sort_keys=True) + '\n').encode()
    digest = _germ_control_hashlib.sha256(raw).hexdigest()
    path = _germ_control_ROOT / f'results/development/artifacts/nsc-incoming-constraint-germ.{digest}.json'
    if path.exists() and path.read_bytes() != raw:
        raise ValueError('content-address collision')
    if not path.exists():
        path.write_bytes(raw)
    return path

def _germ_control_reference_replay(reference):
    total = _germ_control_np.zeros((4, 4))
    first_order = 0.0
    for sample in reference['quadrature_samples']:
        insertion = _germ_control_np.asarray(sample['insertions'])
        weights = _germ_control_np.asarray(sample['weights'])
        k = _germ_control_np.asarray(sample['momenta'])
        if not _germ_control_np.array_equal(k[:len(k) // 2], -k[len(k) // 2:]):
            raise ValueError('paired full momentum line required')
        first_order = max(first_order, float(_germ_control_np.max(abs(insertion[0]))))
        total += sample['multiplicity'] / (2 * _germ_control_np.pi) * _germ_control_np.einsum('n,onb->ob', weights, insertion)
    result = total.sum(axis=0)
    error = float(_germ_control_np.max(abs(result - reference['reference_action_gradient_change'])))
    return (result, error, first_order)

def _germ_control_computed_replay(path):
    _germ_control_authenticate()
    payload = _germ_control_json.loads(path.read_text())
    if payload['signature'] != _germ_control_signature():
        raise ValueError('constraint-germ dependency changed')
    controls = payload['controls']
    u = controls['u']
    mixed = controls['mixed']
    old = _germ_control_read(_germ_control_INPUTS[0])
    baseline = _germ_control_read(_germ_control_INPUTS[1])['baseline']
    base = _germ_control_np.asarray(baseline['local_action_gradient'])
    expectations = {'u': {('r', 3, 0): 0.01}, 'mixed': {('r', 3, 0): 0.02, ('r', 1, 1): -0.005}}
    slot_roundoff = 0.0
    actual_changes = {}
    for name, control in controls.items():
        actual = {(row['field'], row['normal_order'], row['spatial_order']): row['delta'] for row in control['probe']}
        if actual.keys() != expectations[name].keys():
            raise ValueError('owned response plane changed')
        slot_roundoff = max(slot_roundoff, max((abs(actual[k] - expectations[name][k]) for k in actual)))
        if slot_roundoff > 3e-15:
            raise ValueError('owned response values changed beyond derivative-slot rounding')
        actual_changes[name] = actual
    u16, err16, first16 = _germ_control_reference_replay(u['references']['16'])
    u24, err24, first24 = _germ_control_reference_replay(u['references']['24'])
    mixref, errmix, firstmix = _germ_control_reference_replay(mixed['reference'])
    response_u = _germ_control_np.asarray(u['local']['local_action_gradient']) - base + u24[:2]
    response_v = _germ_control_np.asarray(old['local_probe']['local_action_gradient']) - base + _germ_control_np.asarray(old['reference_controls']['n24']['reference_action_gradient_change'])[:2]
    response_mixed = _germ_control_np.asarray(mixed['local']['local_action_gradient']) - base + mixref[:2]
    probe_u = actual_changes['u']['r', 3, 0]
    mixed_u = actual_changes['mixed']['r', 3, 0]
    mixed_v = actual_changes['mixed']['r', 1, 1]
    cu = response_u[0] / probe_u
    cv2 = response_v[0] / 0.01 ** 2
    cv = response_v[1] / 0.01
    predicted = _germ_control_np.array([cu * mixed_u + cv2 * mixed_v ** 2, cv * mixed_v])
    difference = response_mixed - predicted
    coefficients = {'c_u': float(cu), 'c_v2': float(cv2), 'c_v': float(cv)}
    residuals = {'normal_slot_roundoff': slot_roundoff, 'reference_node_replay': max(err16, err24, errmix), 'reference16_24_difference': float(_germ_control_np.max(abs(u16 - u24))), 'homogeneous_shift_zero_identity': float(abs(response_u[1])), 'independent_mixed_control': float(_germ_control_np.max(abs(difference))), 'first_reference_order_unchanged': max(first16, first24, firstmix), 'local_resolution_indicator': max(u['local']['maximum_numerical_indicator'], mixed['local']['maximum_numerical_indicator']), 'local_Euler_identity': max(u['local']['Euler_bulk_identity_residual'], mixed['local']['Euler_bulk_identity_residual'])}
    if max(residuals.values()) > 3e-11:
        raise ArithmeticError('independent response verification failed')
    return {'schema': 'NSC-INCOMING-CONSTRAINT-GERM-v1', 'accountable_author': 'Douglas Ek', 'status': 'PASS: structural plane and numerically independent responses; physical root OPEN', 'source_hashes': {p: _germ_control_sha(p) for p in _germ_control_SOURCES}, 'input_hashes': {p: _germ_control_sha(p) for p in _germ_control_INPUTS}, 'payload': {'path': path.relative_to(_germ_control_ROOT).as_posix(), 'sha256': _germ_control_hashlib.sha256(path.read_bytes()).hexdigest()}, 'structure': _germ_control_structural_constraint_plane(), 'numerical_coefficients': coefficients, 'baseline_Jacobian_from_structure': [[float(cu), 0.0], [0.0, float(cv)]], 'numerical_Jacobian_determinant': float(cu * cv), 'rank_evidence': 'two numerically resolved independent action responses; not an interval nonsingularity certificate', 'u_probe_response': response_u.tolist(), 'v_probe_response': response_v.tolist(), 'mixed_control_response': response_mixed.tolist(), 'mixed_control_prediction': predicted.tolist(), 'mixed_control_residual': difference.tolist(), 'residuals': residuals, 'numerical_tolerance': 3e-11, 'local_probe_records': {'u': u['local'], 'mixed': mixed['local']}, 'reference_integrability': {'fixed_grade0_and_grade1_metric_data': True, 'Delta_P0_and_P1': '0', 'first_changed_order': 2, 'reference_difference': 'O(abs(k)^-3)', 'raw_vertices': 'O(abs(k))', 'insertion_tail': 'O(abs(k)^-2)', 'inventory': 'same finite retained33 groups; exact LLL fixed nonzero-k charts', 'owner': 'results/development/nsc-reference-band-bulk.json'}, 'scope': {'matter_source_cancels_in_response_difference': True, 'physical_state_or_couplings_refit': False, 'physical_initial_data_selected': False, 'numerical_root_selected': False, 'source_accuracy_certified': False, 'Cauchy_surface_constraints_solved': False, 'parent_preparation_or_endpoint_completed': False, 'extended_stationarity': 'OPEN', 'metric_evolution': False, 'Gamma_rest_assigned': False}, 'reproducer': 'python3 scripts/derive_nsc_incoming_constraint_germ.py --check'}

def _germ_control_replay(path):
    from .nsc_frozen_sources import selected_source_context
    from .nsc_weighted_receipt_replay import germ_receipt_bounds, require_record
    current=_germ_control_computed_replay(path)
    if selected_source_context() is None:
        return current
    saved=_germ_control_json.loads(_germ_control_OUTPUT.read_text(encoding='utf-8'))
    if saved['payload']!=current['payload']:
        raise ValueError('original germ payload identity changed')
    payload=_germ_control_json.loads(path.read_text(encoding='utf-8'))
    joint=_germ_control_read(_germ_control_INPUTS[0])
    baseline=_germ_control_read(_germ_control_INPUTS[1])['baseline']
    return require_record(saved,current,germ_receipt_bounds(saved,payload,joint,baseline))

from types import SimpleNamespace as _germ_control_Namespace
germ_control = _germ_control_Namespace(
    INPUTS=_germ_control_INPUTS,
    OUTPUT=_germ_control_OUTPUT,
    Path=_germ_control_Path,
    ROOT=_germ_control_ROOT,
    SOURCES=_germ_control_SOURCES,
    argparse=_germ_control_argparse,
    authenticate=_germ_control_authenticate,
    hashlib=_germ_control_hashlib,
    json=_germ_control_json,
    np=_germ_control_np,
    prepare=_germ_control_prepare,
    read=_germ_control_read,
    reference_replay=_germ_control_reference_replay,
    replay=_germ_control_replay,
    sha=_germ_control_sha,
    signature=_germ_control_signature,
    structural_constraint_plane=_germ_control_structural_constraint_plane,
    sys=_germ_control_sys,
)
