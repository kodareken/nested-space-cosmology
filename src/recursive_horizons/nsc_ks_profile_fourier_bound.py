"""Certified Fourier coefficients of the compact axial potential profiles.

Arb DFT arithmetic is combined with continuous derivative L1 bounds. The
classical alias identity and integration by parts bound coefficients outside
the sampled/retained sets. No empirical Fourier decay fit is used.
"""
from dataclasses import dataclass
from math import factorial

from flint import acb, arb, ctx

from . import nsc_ks_ball_geometry as G


@dataclass(frozen=True)
class ProfileFourierBound:
    key: tuple
    coefficients: object
    derivative_l1_upper: object
    alias_error_upper: object
    uniform_tail_bounds: object
    period_origin: object
    period_length: object
    derivative_order: int
    quadrature_points: int
    retained_index: int
    derivative_cells: int


def alias_and_tail_bounds(derivative_l1, length, derivative_order, quadrature_points,
                          retained_index, max_derivative=2):
    """Bounds with Fourier convention c_k=L^-1 integral f exp(-2 pi i k z/L)."""
    p, M, K = derivative_order, quadrature_points, retained_index
    if any(isinstance(v, bool) or not isinstance(v, int) for v in (p, M, K, max_derivative)):
        raise ValueError('integer derivative and Fourier counts required')
    if not 0 < K < M/2 or p <= max_derivative+1:
        raise ValueError('retained band below Nyquist and integrable derivative tails required')
    L, B = arb(length), arb(derivative_l1)
    if not L > 0 or not B >= 0:
        raise ValueError('positive period and nonnegative derivative bound required')
    scale = B/L*(L/(2*arb.pi()))**p
    alias = (2*scale*p/((p-1)*arb(M-K)**p)).upper()
    tails = tuple((2*scale*(2*arb.pi()/L)**j*arb(K)**(j-p+1)/(p-j-1)).upper()
                  for j in range(max_derivative+1))
    return alias, tails


def _common_window(model):
    profiles = [profile for d in model.metric.directions for profile in (d.w, d.U)]
    windows = {(float(p.center), float(p.inner), float(p.outer)) for p in profiles}
    if len(windows) != 1:
        raise ValueError('the approved common axial cutoff is required')
    return next(iter(windows))


def _derivative_l1(model, keys, order, transition_panels):
    center, inner, outer = map(arb, _common_window(model))
    sums = {key: arb(0) for key in keys}
    intervals = []
    for left, right, count in ((-outer,-inner,transition_panels), (-inner,inner,4),
                                (inner,outer,transition_panels)):
        for i in range(count):
            intervals.append((center+left+(right-left)*i/count,
                              center+left+(right-left)*(i+1)/count, 0))
    cells = 0
    while intervals:
        left, right, depth = intervals.pop()
        try:
            with G._JetWork(order):
                w, U = model.profile_series(left.union(right),order)
                values = {key: G._coeffs(w**key[0]*U**key[1],order)[order]*factorial(order)
                          for key in keys}
        except G.SubdivisionNeeded:
            if depth >= 12:
                raise
            middle = (left+right)/2
            intervals.extend(((left,middle,depth+1),(middle,right,depth+1)))
            continue
        for key, value in values.items():
            sums[key] += (right-left)*abs(value).upper()
        cells += 1
    return {key: value.upper() for key,value in sums.items()}, cells


def enclose_profile_fourier(model, keys, period_origin, period_length, *,
                            derivative_order=8, transition_panels=256,
                            quadrature_points=16384, retained_index=2048,
                            max_derivative=2, bits=120):
    """Return retained coefficient balls and continuous tails for w^p U^q."""
    if (not keys or len(set(keys)) != len(keys) or
            any(len(k)!=2 or min(k)<0 or sum(k)<1 for k in keys)):
        raise ValueError('distinct positive-total-degree profile powers required')
    if transition_panels < 1 or not isinstance(transition_panels,int):
        raise ValueError('positive transition panel count required')
    with ctx.workprec(bits):
        origin, length = arb(period_origin), arb(period_length)
        center, inner, outer = map(arb,_common_window(model))
        if not origin < center-outer or not origin+length > center+outer:
            raise ValueError('compact support must lie strictly inside the numerical period')
        l1, cells = _derivative_l1(model,keys,derivative_order,transition_panels)
        errors = {key:alias_and_tail_bounds(l1[key],length,derivative_order,
                    quadrature_points,retained_index,max_derivative) for key in keys}
        samples = {key:[] for key in keys}
        for node in range(quadrature_points):
            z = origin+length*node/quadrature_points
            if z <= center-outer or z >= center+outer:
                values = {key:arb(0) for key in keys}
            else:
                w, U = model.profile_series(z,0)
                values = {key:w[0]**key[0]*U[0]**key[1] for key in keys}
            for key in keys:
                samples[key].append(acb(values[key]))
        results = {}
        for key in keys:
            transformed = acb.dft(samples[key])
            alias, tails = errors[key]
            error = arb((0,0),alias.man_exp())
            coefficients = tuple(transformed[k % quadrature_points]/quadrature_points+acb(error,error)
                                 for k in range(-retained_index,retained_index+1))
            results[key] = ProfileFourierBound(key,coefficients,l1[key],alias,tails,
                origin,length,derivative_order,quadrature_points,retained_index,cells)
        return results

from pathlib import Path as _receipt_Path
from .nsc_frozen_sources import selected_source_sha256 as _receipt_source_sha256

# Retained executable control: refinement_control
# Folded from scripts/derive_nsc_ks_profile_fourier_refinement.py at 329ef809; original receipts retain their producer identities.
import os as _profile_refinement_control_os
import argparse as _profile_refinement_control_argparse
from hashlib import sha256 as _profile_refinement_control_sha256
import json as _profile_refinement_control_json
from pathlib import Path as _profile_refinement_control_Path
import signal as _profile_refinement_control_signal
import sys as _profile_refinement_control_sys
import time as _profile_refinement_control_time
import numpy as _profile_refinement_control_np
import flint as _profile_refinement_control_flint
from flint import arb as _profile_refinement_control_arb, acb as _profile_refinement_control_acb, ctx as _profile_refinement_control_ctx
from recursive_horizons.nsc_posix import ROOT as _profile_refinement_control_ROOT
_profile_refinement_control_sys.path.insert(0, str(_profile_refinement_control_ROOT / 'src'))
_profile_refinement_control_sys.path.insert(0, str(_profile_refinement_control_ROOT / 'scripts'))
import derive_nsc_ks_profile_fourier_bound as _profile_refinement_control_C
from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily as _profile_refinement_control_AnalyticRadiusFamily
from recursive_horizons.nsc_ks_profile_fourier_bound import enclose_profile_fourier as _profile_refinement_control_enclose_profile_fourier, alias_and_tail_bounds as _profile_refinement_control_alias_and_tail_bounds
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper as _profile_refinement_control_exact_upper, restored_upper as _profile_refinement_control_restored_upper
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes as _profile_refinement_control_deterministic_npz_bytes
_profile_refinement_control_OUTPUT = 'results/development/nsc-ks-profile-fourier-refinement.json'
_profile_refinement_control_PAYLOAD = 'results/development/artifacts/nsc-ks-profile-fourier-refinement.npz'
_profile_refinement_control_OWNED = ('scripts/derive_nsc_ks_profile_fourier_refinement.py', 'docs/nsc-ks-profile-fourier-refinement.md')
_profile_refinement_control_INPUTS = tuple(dict.fromkeys((*_profile_refinement_control_C.OWNED, *_profile_refinement_control_C.INPUTS, _profile_refinement_control_C.OUTPUT)))
_profile_refinement_control_CAP, _profile_refinement_control_BITS = (120.0, 120)
_profile_refinement_control_SETTINGS = {'derivative_order': 16, 'transition_panels': 1024, 'quadrature_points': 65536, 'retained_index': 8192, 'max_derivative': 2, 'bits': _profile_refinement_control_BITS}

def _profile_refinement_control_digest(path):
    candidate = _receipt_Path(path)
    relative = candidate.relative_to(_profile_refinement_control_ROOT).as_posix() if candidate.is_absolute() else candidate.as_posix()
    return _receipt_source_sha256(_profile_refinement_control_ROOT, relative)

def _profile_refinement_control_signatures():
    return {p: _profile_refinement_control_digest(p) for p in (*_profile_refinement_control_OWNED, *_profile_refinement_control_INPUTS)}

def _profile_refinement_control_measure(arrays, meta):
    rows = {}
    with _profile_refinement_control_ctx.workprec(_profile_refinement_control_BITS):
        for key, data in meta['profiles'].items():
            B = _profile_refinement_control_restored_upper(data['derivative_l1_upper'])
            alias, tails = _profile_refinement_control_alias_and_tail_bounds(B, _profile_refinement_control_arb(meta['period_length']), _profile_refinement_control_SETTINGS['derivative_order'], _profile_refinement_control_SETTINGS['quadrature_points'], _profile_refinement_control_SETTINGS['retained_index'], _profile_refinement_control_SETTINGS['max_derivative'])
            if _profile_refinement_control_exact_upper(alias) != data['alias_error_upper'] or [_profile_refinement_control_exact_upper(v) for v in tails] != data['tail_upper']:
                raise ValueError('refined alias/tail replay differs')
            coefficients = [_profile_refinement_control_C.unpack_ball(row) for row in arrays[key + '/coefficients']]
            K = _profile_refinement_control_SETTINGS['retained_index']
            if len(coefficients) != 2 * K + 1 or any((not v.is_finite() for v in coefficients)):
                raise ValueError('refined coefficient inventory is incomplete')
            if not coefficients[K].imag.contains(0) or any((not coefficients[K + k].overlaps(coefficients[K - k].conjugate()) for k in range(1, K + 1))):
                raise ValueError('refined real-profile conjugacy failed')
            old = [_profile_refinement_control_C.unpack_ball(row) for row in arrays[key + '/previous_coefficients']]
            oldK = (len(old) - 1) // 2
            if any((not coefficients[K + k].overlaps(old[oldK + k]) for k in range(-oldK, oldK + 1))):
                raise ValueError('the two certified profile intervals disagree')
            rows[key] = {'retained_coefficients': len(coefficients), 'derivative_cells': data['derivative_cells'], 'tail_sup_orders_0_1_2': [float(v) for v in tails], 'alias_error_upper': float(alias), 'previous_coefficients_compatible': True}
    return rows

def _profile_refinement_control_record(arrays, meta, payload):
    return {'schema': 'NSC-KS-PROFILE-FOURIER-REFINEMENT-v1', 'accountable_author': 'Douglas Ek', 'status': 'PASS: tighter certified axial-profile enclosure; physical local gate OPEN', 'settings': _profile_refinement_control_SETTINGS, 'measurements': _profile_refinement_control_measure(arrays, meta), 'runtime': meta['runtime'], 'scope': {'source_energy_UV_certificate': False, 'continuous_state_error_bound': None, 'physical_local_gate': 'OPEN', 'field_runs': 0, 'source_runs': 0, 'metric_timestep': False, 'geometric_profile_only': True}, 'source_hashes': {p: meta['signature'][p] for p in _profile_refinement_control_OWNED}, 'input_hashes': {p: meta['signature'][p] for p in _profile_refinement_control_INPUTS}, 'reused_artifacts': meta['reused_artifacts'], 'payload': payload, 'reproducer': 'python scripts/derive_nsc_ks_profile_fourier_refinement.py --check'}

def _profile_refinement_control_run():
    if (_profile_refinement_control_ROOT / _profile_refinement_control_OUTPUT).exists() or (_profile_refinement_control_ROOT / _profile_refinement_control_PAYLOAD).exists():
        raise FileExistsError('refinement exists; use --check')
    sig = _profile_refinement_control_signatures()
    previous = _profile_refinement_control_C.check()
    with _profile_refinement_control_np.load(_profile_refinement_control_ROOT / previous['payload']['path'], allow_pickle=False) as f:
        old = {k: f[k] for k in f.files}
    oldmeta = _profile_refinement_control_json.loads(old['metadata_json'].tobytes())
    trajectory = oldmeta['reused_artifact']
    if _profile_refinement_control_digest(trajectory['path']) != trajectory['sha256']:
        raise ValueError('trajectory input changed')
    with _profile_refinement_control_np.load(_profile_refinement_control_ROOT / trajectory['path'], allow_pickle=False) as f:
        a = {k: f[k] for k in f.files}
    prepared, family, _ = _profile_refinement_control_C.C.unpack(a, _profile_refinement_control_json.loads(a['metadata_json'].tobytes()))
    cpu, wall = (_profile_refinement_control_time.process_time(), _profile_refinement_control_time.monotonic())

    def stop(*_):
        raise RuntimeError('geometric Fourier refinement exceeded120 CPU seconds')
    prior = _profile_refinement_control_signal.signal(_profile_refinement_control_signal.SIGPROF, stop)
    _profile_refinement_control_signal.setitimer(_profile_refinement_control_signal.ITIMER_PROF, _profile_refinement_control_CAP)
    try:
        result = _profile_refinement_control_enclose_profile_fourier(_profile_refinement_control_AnalyticRadiusFamily(family), [(n, 0) for n in range(1, 5)], oldmeta['period_origin'], oldmeta['period_length'], **_profile_refinement_control_SETTINGS)
    finally:
        _profile_refinement_control_signal.setitimer(_profile_refinement_control_signal.ITIMER_PROF, 0)
        _profile_refinement_control_signal.signal(_profile_refinement_control_signal.SIGPROF, prior)
    arrays = {}
    profiles = {}
    with _profile_refinement_control_ctx.workprec(_profile_refinement_control_BITS):
        for (p, q), row in result.items():
            key = f'profile_{p}_{q}'
            arrays[key + '/coefficients'] = _profile_refinement_control_np.asarray([_profile_refinement_control_C.pack_ball(v) for v in row.coefficients])
            arrays[key + '/previous_coefficients'] = old[key + '/coefficients']
            profiles[key] = {'powers': [p, q], 'derivative_l1_upper': _profile_refinement_control_exact_upper(row.derivative_l1_upper), 'alias_error_upper': _profile_refinement_control_exact_upper(row.alias_error_upper), 'tail_upper': [_profile_refinement_control_exact_upper(v) for v in row.uniform_tail_bounds], 'derivative_cells': row.derivative_cells}
    meta = {'signature': sig, 'profiles': profiles, 'period_origin': oldmeta['period_origin'], 'period_length': oldmeta['period_length'], 'reused_artifacts': [previous['payload'], trajectory], 'prepared_source_digest': prepared.fixed_preparation_digest, 'python_flint_version': _profile_refinement_control_flint.__version__, 'runtime': {'CPU_seconds': _profile_refinement_control_time.process_time() - cpu, 'wall_seconds': _profile_refinement_control_time.monotonic() - wall, 'CPU_cap': _profile_refinement_control_CAP}}
    if _profile_refinement_control_signatures() != sig:
        raise ValueError('profile refinement source/input changed')
    arrays['metadata_json'] = _profile_refinement_control_np.frombuffer(_profile_refinement_control_json.dumps(meta, sort_keys=True).encode(), _profile_refinement_control_np.uint8)
    raw = _profile_refinement_control_deterministic_npz_bytes(arrays)
    payload = {'path': _profile_refinement_control_PAYLOAD, 'sha256': _profile_refinement_control_sha256(raw).hexdigest(), 'bytes': len(raw)}
    output = _profile_refinement_control_record(arrays, meta, payload)
    (_profile_refinement_control_ROOT / _profile_refinement_control_PAYLOAD).write_bytes(raw)
    (_profile_refinement_control_ROOT / _profile_refinement_control_OUTPUT).write_text(_profile_refinement_control_json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return output

def _profile_refinement_control_check():
    saved = _profile_refinement_control_json.loads((_profile_refinement_control_ROOT / _profile_refinement_control_OUTPUT).read_text())
    _profile_refinement_control_C.check()
    if _profile_refinement_control_digest(_profile_refinement_control_PAYLOAD) != saved['payload']['sha256']:
        raise ValueError('refined Fourier payload changed')
    with _profile_refinement_control_np.load(_profile_refinement_control_ROOT / _profile_refinement_control_PAYLOAD, allow_pickle=False) as f:
        a = {k: f[k] for k in f.files}
    meta = _profile_refinement_control_json.loads(a['metadata_json'].tobytes())
    if meta['signature'] != _profile_refinement_control_signatures():
        raise ValueError('refined Fourier source/input changed')
    if _profile_refinement_control_record(a, meta, saved['payload']) != saved:
        raise ValueError('refined Fourier replay differs')
    return saved

from types import SimpleNamespace as _profile_refinement_control_Namespace
refinement_control = _profile_refinement_control_Namespace(
    AnalyticRadiusFamily=_profile_refinement_control_AnalyticRadiusFamily,
    BITS=_profile_refinement_control_BITS,
    C=_profile_refinement_control_C,
    CAP=_profile_refinement_control_CAP,
    INPUTS=_profile_refinement_control_INPUTS,
    OUTPUT=_profile_refinement_control_OUTPUT,
    OWNED=_profile_refinement_control_OWNED,
    PAYLOAD=_profile_refinement_control_PAYLOAD,
    Path=_profile_refinement_control_Path,
    ROOT=_profile_refinement_control_ROOT,
    SETTINGS=_profile_refinement_control_SETTINGS,
    acb=_profile_refinement_control_acb,
    alias_and_tail_bounds=_profile_refinement_control_alias_and_tail_bounds,
    arb=_profile_refinement_control_arb,
    argparse=_profile_refinement_control_argparse,
    check=_profile_refinement_control_check,
    ctx=_profile_refinement_control_ctx,
    deterministic_npz_bytes=_profile_refinement_control_deterministic_npz_bytes,
    digest=_profile_refinement_control_digest,
    enclose_profile_fourier=_profile_refinement_control_enclose_profile_fourier,
    exact_upper=_profile_refinement_control_exact_upper,
    flint=_profile_refinement_control_flint,
    json=_profile_refinement_control_json,
    measure=_profile_refinement_control_measure,
    np=_profile_refinement_control_np,
    os=_profile_refinement_control_os,
    record=_profile_refinement_control_record,
    restored_upper=_profile_refinement_control_restored_upper,
    run=_profile_refinement_control_run,
    sha256=_profile_refinement_control_sha256,
    signal=_profile_refinement_control_signal,
    signatures=_profile_refinement_control_signatures,
    sys=_profile_refinement_control_sys,
    time=_profile_refinement_control_time,
)

# Retained executable control: fine_profile_control
# Folded from scripts/derive_nsc_ks_fine_profile_table.py at 329ef809; original receipts retain their producer identities.
import os as _fine_profile_control_os
import argparse as _fine_profile_control_argparse
from hashlib import sha256 as _fine_profile_control_sha256
import json as _fine_profile_control_json
from pathlib import Path as _fine_profile_control_Path
import sys as _fine_profile_control_sys
import time as _fine_profile_control_time
import numpy as _fine_profile_control_np
from flint import arb as _fine_profile_control_arb, ctx as _fine_profile_control_ctx
from recursive_horizons.nsc_posix import ROOT as _fine_profile_control_ROOT
_fine_profile_control_sys.path.insert(0, str(_fine_profile_control_ROOT / 'src'))
_fine_profile_control_sys.path.insert(0, str(_fine_profile_control_ROOT / 'scripts'))
import derive_nsc_ks_fine_trajectory as _fine_profile_control_F
import derive_nsc_ks_profile_fourier_bound as _fine_profile_control_P
from recursive_horizons.nsc_ks_profile_fourier_bound import refinement_control as _fine_profile_control_Q
from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily as _fine_profile_control_AnalyticRadiusFamily
from recursive_horizons.nsc_ks_profile_identity import profile_identity as _fine_profile_control_profile_identity, profile_description as _fine_profile_control_profile_description
from recursive_horizons.nsc_ks_profile_fourier_bound import enclose_profile_fourier as _fine_profile_control_enclose_profile_fourier, alias_and_tail_bounds as _fine_profile_control_alias_and_tail_bounds
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper as _fine_profile_control_exact_upper, restored_upper as _fine_profile_control_restored_upper
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes as _fine_profile_control_deterministic_npz_bytes
_fine_profile_control_OUTPUT = _fine_profile_control_ROOT / 'results/development/nsc-ks-fine-profile-table.json'
_fine_profile_control_PAYLOAD = _fine_profile_control_ROOT / 'results/development/artifacts/nsc-ks-fine-profile-table.npz'
_fine_profile_control_OWNED = ('scripts/derive_nsc_ks_fine_profile_table.py', 'docs/nsc-ks-fine-profile-table.md')
_fine_profile_control_INPUTS = tuple(dict.fromkeys((*_fine_profile_control_Q.OWNED, *_fine_profile_control_Q.INPUTS, *_fine_profile_control_F.OWNED, str(_fine_profile_control_F.OUTPUT.relative_to(_fine_profile_control_ROOT)), _fine_profile_control_Q.OUTPUT)))
_fine_profile_control_SETTINGS = _fine_profile_control_Q.SETTINGS

def _fine_profile_control_digest(path):
    candidate = _receipt_Path(path)
    relative = candidate.relative_to(_fine_profile_control_ROOT).as_posix() if candidate.is_absolute() else candidate.as_posix()
    return _receipt_source_sha256(_fine_profile_control_ROOT, relative)

def _fine_profile_control_signatures():
    return {p: _fine_profile_control_digest(p) for p in (*_fine_profile_control_OWNED, *_fine_profile_control_INPUTS)}

def _fine_profile_control_measures(arrays, meta):
    rows = {}
    with _fine_profile_control_ctx.workprec(_fine_profile_control_SETTINGS['bits']):
        for key, proof in meta['profiles'].items():
            alias, tails = _fine_profile_control_alias_and_tail_bounds(_fine_profile_control_restored_upper(proof['derivative_l1_upper']), _fine_profile_control_arb(meta['period_length']), _fine_profile_control_SETTINGS['derivative_order'], _fine_profile_control_SETTINGS['quadrature_points'], _fine_profile_control_SETTINGS['retained_index'], 2)
            if _fine_profile_control_exact_upper(alias) != proof['alias_error_upper'] or [_fine_profile_control_exact_upper(v) for v in tails] != proof['tail_upper']:
                raise ValueError('fine-period Fourier bounds no longer replay')
            coef = [_fine_profile_control_P.unpack_ball(v) for v in arrays[key + '/coefficients']]
            K = _fine_profile_control_SETTINGS['retained_index']
            if len(coef) != 2 * K + 1 or any((not v.is_finite() for v in coef)):
                raise ValueError('incomplete Fourier table')
            if not coef[K].imag.contains(0) or any((not coef[K + k].overlaps(coef[K - k].conjugate()) for k in range(1, K + 1))):
                raise ValueError('real-profile conjugacy failed')
            rows[key] = {'tail_sup_orders_0_1_2': [float(v) for v in tails], 'alias_upper': float(alias)}
    return rows

def _fine_profile_control_record(arrays, meta, payload):
    return {'schema': 'NSC-KS-FINE-PROFILE-TABLE-v1', 'accountable_author': 'Douglas Ek', 'status': 'PASS: fine-period geometric profile enclosure; field and physical gates OPEN', 'settings': _fine_profile_control_SETTINGS, 'measurements': _fine_profile_control_measures(arrays, meta), 'runtime': meta['runtime'], 'period_length_exact': '205/512', 'axial_profile_identity': meta['axial_profile_identity'], 'source_preparation_digest': meta['source_preparation_digest'], 'scope': {'physical_profiles_changed': False, 'physical_local_gate': 'OPEN', 'field_runs': 0, 'source_runs': 0, 'source_energy_UV_certificate': False, 'metric_timestep': False}, 'source_hashes': {p: meta['signature'][p] for p in _fine_profile_control_OWNED}, 'input_hashes': {p: meta['signature'][p] for p in _fine_profile_control_INPUTS}, 'reused_artifacts': meta['reused_artifacts'], 'payload': payload, 'reproducer': 'python scripts/derive_nsc_ks_fine_profile_table.py --check'}

def _fine_profile_control_run():
    if _fine_profile_control_OUTPUT.exists() or _fine_profile_control_PAYLOAD.exists():
        raise FileExistsError('fine profile table exists; use --check')
    sig = _fine_profile_control_signatures()
    fine = _fine_profile_control_json.loads(_fine_profile_control_F.OUTPUT.read_text())
    old = _fine_profile_control_Q.check()
    if not fine['complete_capture']:
        raise ValueError('fine trajectory incomplete')
    if _fine_profile_control_digest(fine['payload']['path']) != fine['payload']['sha256']:
        raise ValueError('fine prepared data changed')
    with _fine_profile_control_np.load(_fine_profile_control_ROOT / fine['payload']['path'], allow_pickle=False) as f:
        fm = _fine_profile_control_json.loads(f['metadata_json'].tobytes())
        grid = f['computational_z']
    with _fine_profile_control_np.load(_fine_profile_control_ROOT / old['payload']['path'], allow_pickle=False) as f:
        oldmeta = _fine_profile_control_json.loads(f['metadata_json'].tobytes())
    family = _fine_profile_control_F.T.C.R.C.P.family(_fine_profile_control_F.T.C.R.C.P.ALPHA)
    identity = _fine_profile_control_profile_identity(family, include_normal_window=False)
    if identity != fine['axial_profile_identity'] or _fine_profile_control_profile_description(family) != fm['history_description']:
        raise ValueError('analytic profile identity differs from the fine trajectory')
    length = float((grid[1] - grid[0]) * len(grid))
    origin = float(grid[0])
    if length != 205 / 512:
        raise ValueError('fine grid period changed')
    cpu, wall = (_fine_profile_control_time.process_time(), _fine_profile_control_time.monotonic())
    result = _fine_profile_control_enclose_profile_fourier(_fine_profile_control_AnalyticRadiusFamily(family), [(n, 0) for n in range(1, 5)], origin, length, **_fine_profile_control_SETTINGS)
    arrays = {}
    proofs = {}
    with _fine_profile_control_ctx.workprec(_fine_profile_control_SETTINGS['bits']):
        for (p, q), item in result.items():
            key = f'profile_{p}_{q}'
            if _fine_profile_control_exact_upper(item.derivative_l1_upper) != oldmeta['profiles'][key]['derivative_l1_upper']:
                raise ValueError('physical derivative bound changed with the numerical period')
            arrays[key + '/coefficients'] = _fine_profile_control_np.asarray([_fine_profile_control_P.pack_ball(v) for v in item.coefficients])
            proofs[key] = {'powers': [p, q], 'derivative_l1_upper': _fine_profile_control_exact_upper(item.derivative_l1_upper), 'alias_error_upper': _fine_profile_control_exact_upper(item.alias_error_upper), 'tail_upper': [_fine_profile_control_exact_upper(v) for v in item.uniform_tail_bounds]}
    meta = {'signature': sig, 'profiles': proofs, 'period_origin': origin, 'period_length': length, 'axial_profile_identity': identity, 'profile_description': _fine_profile_control_profile_description(family, include_normal_window=False), 'source_preparation_digest': fine['source_preparation_digest'], 'reused_artifacts': [fine['payload'], old['payload']], 'runtime': {'CPU_seconds': _fine_profile_control_time.process_time() - cpu, 'wall_seconds': _fine_profile_control_time.monotonic() - wall}}
    if sig != _fine_profile_control_signatures():
        raise ValueError('fine profile owner/input changed')
    arrays['metadata_json'] = _fine_profile_control_np.frombuffer(_fine_profile_control_json.dumps(meta, sort_keys=True).encode(), _fine_profile_control_np.uint8)
    raw = _fine_profile_control_deterministic_npz_bytes(arrays)
    payload = {'path': str(_fine_profile_control_PAYLOAD.relative_to(_fine_profile_control_ROOT)), 'sha256': _fine_profile_control_sha256(raw).hexdigest(), 'bytes': len(raw)}
    output = _fine_profile_control_record(arrays, meta, payload)
    _fine_profile_control_PAYLOAD.write_bytes(raw)
    _fine_profile_control_OUTPUT.write_text(_fine_profile_control_json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return output

def _fine_profile_control_check():
    saved = _fine_profile_control_json.loads(_fine_profile_control_OUTPUT.read_text())
    if _fine_profile_control_digest(_fine_profile_control_PAYLOAD) != saved['payload']['sha256']:
        raise ValueError('fine profile payload changed')
    with _fine_profile_control_np.load(_fine_profile_control_PAYLOAD, allow_pickle=False) as f:
        a = {k: f[k] for k in f.files}
    meta = _fine_profile_control_json.loads(a['metadata_json'].tobytes())
    if meta['signature'] != _fine_profile_control_signatures():
        raise ValueError('fine profile owner/input changed')
    if _fine_profile_control_record(a, meta, saved['payload']) != saved:
        raise ValueError('fine profile replay differs')
    return saved

from types import SimpleNamespace as _fine_profile_control_Namespace
fine_profile_control = _fine_profile_control_Namespace(
    AnalyticRadiusFamily=_fine_profile_control_AnalyticRadiusFamily,
    F=_fine_profile_control_F,
    INPUTS=_fine_profile_control_INPUTS,
    OUTPUT=_fine_profile_control_OUTPUT,
    OWNED=_fine_profile_control_OWNED,
    P=_fine_profile_control_P,
    PAYLOAD=_fine_profile_control_PAYLOAD,
    Path=_fine_profile_control_Path,
    Q=_fine_profile_control_Q,
    ROOT=_fine_profile_control_ROOT,
    SETTINGS=_fine_profile_control_SETTINGS,
    alias_and_tail_bounds=_fine_profile_control_alias_and_tail_bounds,
    arb=_fine_profile_control_arb,
    argparse=_fine_profile_control_argparse,
    check=_fine_profile_control_check,
    ctx=_fine_profile_control_ctx,
    deterministic_npz_bytes=_fine_profile_control_deterministic_npz_bytes,
    digest=_fine_profile_control_digest,
    enclose_profile_fourier=_fine_profile_control_enclose_profile_fourier,
    exact_upper=_fine_profile_control_exact_upper,
    json=_fine_profile_control_json,
    measures=_fine_profile_control_measures,
    np=_fine_profile_control_np,
    os=_fine_profile_control_os,
    profile_description=_fine_profile_control_profile_description,
    profile_identity=_fine_profile_control_profile_identity,
    record=_fine_profile_control_record,
    restored_upper=_fine_profile_control_restored_upper,
    run=_fine_profile_control_run,
    sha256=_fine_profile_control_sha256,
    signatures=_fine_profile_control_signatures,
    sys=_fine_profile_control_sys,
    time=_fine_profile_control_time,
)
