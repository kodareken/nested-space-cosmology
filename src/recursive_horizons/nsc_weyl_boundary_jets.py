"""Boundary jets of the already owned NSC homogeneous Weyl action.

Imports the higher-derivative first-variation identity.  This module applies
it to the existing density and decomposes its stored discrete node gradient;
it does not add an action term or solve a metric field equation.
"""
import numpy as np

FIELDS=('lapse','beta','q_ADM','r')


def weyl_density(jets,coefficient):
    """Existing local-owner L_W, with its eight local jets independent.

    Jet order: N,Ndot,a,adot,addot,r,rdot,rddot.  q_ADM=a.
    """
    N,Nd,a,ad,add,r,rd,rdd=np.asarray(jets).T
    A=(add-ad*Nd/N)/(N*N*a)
    B=ad*rd/(N*N*a*r)
    C=(rdd-rd*Nd/N)/(N*N*r)
    D=(1+(rd/N)**2)/(r*r)
    return -(16*np.pi/3)*coefficient*N*a*r*r*(A-B-C+D)**2


def density_partials(jets,coefficient):
    x=np.asarray(jets,float)
    if x.ndim!=2 or x.shape[1]!=8 or not np.isfinite(x).all() or np.any(x[:,(0,2,5)]<=0):
        raise ValueError('finite independent Weyl jets with positive N,a,r required')
    out=np.empty_like(x);h=1e-28
    for j in range(8):
        z=x.astype(complex);z[:,j]+=1j*h
        out[:,j]=weyl_density(z,coefficient).imag/h
    return out


def nodal_decomposition(history,coefficient):
    """Separate bulk, canonical end jets and the owned stencil's defect.

    G=W E+S P+D^T S Q+(B-S)P+D^T(B-S)Q, with
    B=WD+D^TW, P=L_gdot-D L_gddot. This is the exact chain rule for the
    existing trapezoid/np.gradient discretization. B is not silently set to S.
    """
    t=np.asarray(history.time,float);n=len(t)
    if n<7 or np.any(np.diff(t)<=0):raise ValueError('owned increasing KS node grid required')
    D=np.gradient(np.eye(n),t,axis=0,edge_order=2);D2=D@D
    w=np.empty(n);w[0]=(t[1]-t[0])/2;w[-1]=(t[-1]-t[-2])/2;w[1:-1]=(t[2:]-t[:-2])/2
    W=np.diag(w);S=np.diag(np.r_[-1.,np.zeros(n-2),1.]);B=W@D+D.T@W
    N,a,r=map(np.asarray,(history.lapse,history.a_parallel,history.radius))
    Nd,ad,rd=(D@x for x in (N,a,r));add,rdd=D2@a,D2@r
    jets=np.column_stack((N,Nd,a,ad,add,r,rd,rdd));partial=density_partials(jets,coefficient)
    first=np.zeros((4,n));second=first.copy();value=first.copy()
    for field,j0,j1,j2 in ((0,0,1,None),(2,2,3,4),(3,5,6,7)):
        value[field]=partial[:,j0];first[field]=partial[:,j1]
        if j2 is not None:second[field]=partial[:,j2]
    P=first-second@D.T;Q=second
    bulk=(value-P@D.T)@W
    canonical=P@S+Q@S@D
    defect=P@(B-S).T+Q@(B-S).T@D
    chain=value@W+first@W@D+second@W@D2
    reconstructed=bulk+canonical+defect
    # In proper normal variables delta adot=N*a*delta H_a + ... .
    proper_N=P[0]+ad/N*Q[2]+rd/N*Q[3]
    proper_a=P[2]+ad/a*Q[2];proper_r=P[3]+rd/r*Q[3]
    shear= N*a*Q[2]
    k=-(16*np.pi/3)*coefficient
    L=weyl_density(jets,coefficient)
    # Independent explicit highest-jet coefficients from the same density.
    F=((add-ad*Nd/N)/(N*N*a)-ad*rd/(N*N*a*r)
       -(rdd-rd*Nd/N)/(N*N*r)+(1+(rd/N)**2)/(r*r))
    expected_Qa=2*k*r*r*F/N;expected_Qr=-2*k*a*r*F/N
    return {'time':t,'jets':jets,'density':L,'density_partials':partial,
      'gradient':chain,'bulk_component':bulk,'canonical_boundary_component':canonical,
      'stencil_component':defect,'reconstructed_gradient':reconstructed,
      'P':P,'Q':Q,'proper_boundary_N':proper_N,'proper_boundary_a':proper_a,
      'proper_boundary_r':proper_r,'proper_boundary_shear':shear,
      'quadrature_weights':w,'derivative_matrix':D,'boundary_matrix':B,'ideal_boundary_matrix':S,
      'residuals':{'gradient_decomposition':float(np.max(abs(chain-reconstructed))),
          'highest_jet_Qa':float(np.max(abs(Q[2]-expected_Qa))),
          'highest_jet_Qr':float(np.max(abs(Q[3]-expected_Qr))),
          'Weyl_normal_trace':float(np.max(abs(a*Q[2]+r*Q[3]))),
          'proper_boundary_lapse':float(np.max(abs(proper_N))),
          'stencil_SBP_defect_norm':float(np.linalg.norm(B-S,2))}}

# Retained executable control: boundary_control
# Folded from scripts/derive_nsc_weyl_boundary_jets.py at 329ef809; original receipts retain their producer identities.
import hashlib as _boundary_control_hashlib
import json as _boundary_control_json
from pathlib import Path as _boundary_control_Path
import sys as _boundary_control_sys
from types import SimpleNamespace as _boundary_control_SimpleNamespace
import numpy as _boundary_control_np
from recursive_horizons.nsc_posix import ROOT as _boundary_control_ROOT
_boundary_control_sys.path.insert(0, str(_boundary_control_ROOT / 'src'))
from recursive_horizons.nsc_weyl_boundary_jets import nodal_decomposition as _boundary_control_nodal_decomposition, FIELDS as _boundary_control_FIELDS
_boundary_control_OUTPUT = 'results/development/nsc-weyl-boundary-jets.json'
_boundary_control_INPUTS = ('results/development/nsc-general-ks-local-history.json', 'results/development/nsc-landau-cauchy-isometry.json', 'results/development/nsc-transmitting-boundary-remainder.json', 'results/development/compact-boundary-action.json', 'src/recursive_horizons/nsc_general_ks_local_history.py', 'src/recursive_horizons/nsc_boundary_state.py', 'src/recursive_horizons/nsc_causal_common.py', 'src/recursive_horizons/nsc_transmitting_boundary_history.py', 'scripts/derive_nsc_transmitting_boundary_binding.py')
_boundary_control_SOURCES = ('src/recursive_horizons/nsc_weyl_boundary_jets.py', 'scripts/derive_nsc_weyl_boundary_jets.py', 'tests/test_nsc_weyl_boundary_jets.py', 'docs/nsc-weyl-boundary-jets.md')

def _boundary_control_sha(p):
    from .nsc_frozen_sources import selected_source_sha256
    return selected_source_sha256(_boundary_control_ROOT, p)

def _boundary_control_rec(p):
    return _boundary_control_json.loads((_boundary_control_ROOT / p).read_text())

def _boundary_control_calculate():
    local = _boundary_control_rec(_boundary_control_INPUTS[0])
    isometry = _boundary_control_rec(_boundary_control_INPUTS[1])
    binding = _boundary_control_rec(_boundary_control_INPUTS[2])
    payload = isometry['payload']
    if _boundary_control_sha(payload['path']) != payload['sha256']:
        raise ValueError('authenticated history changed')
    with _boundary_control_np.load(_boundary_control_ROOT / payload['path'], allow_pickle=False) as a:
        index = _boundary_control_np.linspace(0, len(a['short_time']) - 1, 33, dtype=int)
        h = _boundary_control_SimpleNamespace(time=a['short_time'][index], lapse=a['short_lapse'][index], shift=a['short_shift'][index], a_parallel=a['short_a'][index], radius=a['short_r'][index])
    result = _boundary_control_nodal_decomposition(h, local['locked_inputs']['C_Weyl'])
    original = _boundary_control_np.array([local['diagnostic_evaluation']['action_gradients']['weyl_bulk'][f] for f in _boundary_control_FIELDS])
    aliases = ('N', 'beta', 'q_ADM', 'r')
    old_endpoints = {f: original[i, [0, -1]].tolist() for i, f in enumerate(aliases)}
    if old_endpoints != binding['locked_local_components']:
        raise ValueError('the shared eight-coefficient object changed')
    residuals = dict(result['residuals'])
    residuals['full_gradient_vs_locked'] = float(_boundary_control_np.max(abs(result['gradient'] - original)))
    residuals['legacy_eight_coefficient_binding'] = 0.0
    tests = {k: v for k, v in residuals.items() if k != 'stencil_SBP_defect_norm'}
    if any((not _boundary_control_np.isfinite(v) or v > 3e-11 for v in tests.values())):
        raise ArithmeticError('boundary-jet reconstruction failed')
    proper = {k: result[k][[0, -1]].tolist() for k in ('proper_boundary_N', 'proper_boundary_a', 'proper_boundary_r', 'proper_boundary_shear')}
    raw = {k: {f: result[k][i, [0, -1]].tolist() for i, f in enumerate(aliases)} for k in ('P', 'Q')}
    parts = {k: {f: result[k][i, [0, -1]].tolist() for i, f in enumerate(aliases)} for k in ('bulk_component', 'canonical_boundary_component', 'stencil_component', 'reconstructed_gradient')}
    record = {'schema': 'NSC-WEYL-BOUNDARY-JETS-v1', 'accountable_author': 'Douglas Ek', 'status': 'PASS: intrinsic boundary-jet extraction and discrete reconstruction; physical interface stationarity OPEN', 'source_hashes': {p: _boundary_control_sha(p) for p in _boundary_control_SOURCES}, 'input_hashes': {p: _boundary_control_sha(p) for p in _boundary_control_INPUTS}, 'history_payload': payload, 'locked_inputs': local['locked_inputs'], 'legacy_endpoint_basis': binding['endpoint_basis'], 'legacy_eight_coefficients': old_endpoints, 'preserved_Weyl_value': binding['preserved_Weyl_value'], 'boundary_jet_basis': {'raw_values': ['N', 'beta', 'a=q_ADM', 'r'], 'independent_coordinate_derivatives': ['adot', 'rdot'], 'orientation': '[boundary form] at final minus initial; P and Q below precede this orientation', 'proper_normal_rates': 'H_a=adot/(N*a), H_r=rdot/(N*r)', 'proper_shear': 'sigma=H_a-H_r', 'normal_lapse_coefficient': 'zero by the same Weyl density identity', 'not_a_Cauchy_state_isometry': True}, 'raw_boundary_coefficients': raw, 'proper_boundary_coefficients': proper, 'endpoint_nodal_decomposition': parts, 'full_nodal_decomposition': {k: result[k].tolist() for k in ('gradient', 'bulk_component', 'canonical_boundary_component', 'stencil_component', 'reconstructed_gradient')}, 'identities': {'first_variation': 'delta S=bulk+[P_N delta N+P_a delta a+P_r delta r+Q_a delta adot+Q_r delta rdot]', 'proper_boundary': '[Pbar_a delta a+Pbar_r delta r+Pi_sigma delta(H_a-H_r)]', 'highest_jet': 'a*Q_a+r*Q_r=0; N*P_N+adot*Q_a+rdot*Q_r=0', 'discrete_gradient': 'G=W E+S P+D^T S Q+(B-S)P+D^T(B-S)Q; B=WD+D^TW', 'reason_for_stencil_term': 'the existing edge_order=2 gradient with trapezoid weights is not an exact endpoint-only summation-by-parts operator'}, 'residuals': residuals, 'tolerance': 3e-11, 'variational_domain': {'fixed_intrinsic_metric_and_normal_rates': 'canonical Weyl boundary form vanishes on these variations; not a self-sourcing claim', 'fixed_intrinsic_metric_only': 'free normal shear still has coefficient Pi_sigma and cannot be discarded', 'free_or_glued_interface': 'requires the actual matching of intrinsic metric, normal shear and their oriented conjugates', 'physical_matching_data_supplied': False, 'time_end_nodes_are_not_the_two_sides_of_rho_zero': True}, 'remaining_action_ownership': {'Gaussian_complement': 'GaussianBoundaryState.overlap_schur factorizes the same Gamma_G; it is not a second added action', 'geometric_terms': 'Euler/boxR endpoints already belong to GeneralKSLocalInducedHistory and are not introduced again', 'Gamma_rest': 'the declared remainder excludes geometric pi but its full transmitting functional and permitted boundary variations are not specified by these records', 'Gamma_rest_assigned_zero': False, 'arbitrary_B_introduced': False}, 'gate': {'B1_full': 'OPEN: physical history/endpoint Cauchy embedding remains distinct from nodal bookkeeping', 'B2': 'OPEN: intrinsic boundary-jet matching and remaining same-action derivative not evaluated', 'physical_two_sided_mismatch': None, 'stationarity': 'OPEN', 'history_selected': False, 'metric_timestep': False, 'new_physical_term': False, 'stress_computed': False, 'Z3': 'OUT OF SCOPE'}, 'comparison': {'fields': 'all', 'float_atol': 3e-13, 'float_rtol': 3e-13, 'exact': 'hashes, original eight coefficients, schemas and scope'}, 'reproducer': 'python3 scripts/derive_nsc_weyl_boundary_jets.py --check'}
    return record

from types import SimpleNamespace as _boundary_control_Namespace
boundary_control = _boundary_control_Namespace(
    FIELDS=_boundary_control_FIELDS,
    INPUTS=_boundary_control_INPUTS,
    OUTPUT=_boundary_control_OUTPUT,
    Path=_boundary_control_Path,
    ROOT=_boundary_control_ROOT,
    SOURCES=_boundary_control_SOURCES,
    SimpleNamespace=_boundary_control_SimpleNamespace,
    calculate=_boundary_control_calculate,
    hashlib=_boundary_control_hashlib,
    json=_boundary_control_json,
    nodal_decomposition=_boundary_control_nodal_decomposition,
    np=_boundary_control_np,
    rec=_boundary_control_rec,
    sha=_boundary_control_sha,
    sys=_boundary_control_sys,
)
