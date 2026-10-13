"""Short-shell Fourier response in the owned unitary KS radial frame.

Real energies label original input fibers and output Fourier arguments.
Neither source occupations nor energy quadrature weights enter this ODE.
No physical preparation accuracy is inferred from its numerical tolerances.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from .nsc_compatible_history_geometry import CompatibleRadiusDirection
from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_lorentzian import geometry
from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain,MODE_TO_CURRENT,I2,S1,S2,S3
from .nsc_transmitting_resolvent import profile


def canonical_map(rho):
    beta=float(geometry(rho)[0]);r=float(np.sqrt(1+rho*rho))
    domain=TransmittingDiracSeamDomain(lapse=1.,radial_scale=1.,shift=beta,radius=r)
    _,inverse,_=domain.trace_map(np.ones(1))
    return inverse[0]/r


def archived_amplitudes(x,phi,energies,index):
    """Spin-major characteristic archive -> unweighted KS source amplitudes."""
    x=np.asarray(x);E=np.asarray(energies);phi=np.asarray(phi)
    if E.shape!=(4,) or phi.shape!=(2*len(x),12) or not 0<=index<len(x):
        raise ValueError('four energies and three archived source columns required')
    current=MODE_TO_CURRENT@phi[[index,len(x)+index]]
    mapped=(canonical_map(float(x[index]))@current).reshape(2,4,3).transpose(1,0,2)
    return mapped*np.exp(1j*E*chart_coordinates(float(x[index]))[1])[:,None,None]


def pulse_fourier(omega,center,halfwidth,order):
    """Integral exp(+i omega z) w(z) dz; inverse measure is d omega/(2 pi)."""
    omega=np.asarray(omega)
    if np.iscomplexobj(omega) or not np.isfinite(omega).all() or halfwidth<=0 or order not in (96,192):
        raise ValueError('real finite Fourier arguments and fixed96/192 rules required')
    nodes,weights=leggauss(order)
    bump=np.array([profile(float(u))[0] for u in nodes])
    return halfwidth*np.exp(1j*omega*center)*np.sum(
        np.exp(1j*omega[...,None]*halfwidth*nodes)*weights*bump,axis=-1)


def ks_generator(rho,energies,mass,angular):
    beta=float(geometry(rho)[0]);a=np.sqrt(beta*beta-1);r=np.sqrt(1+rho*rho)
    E=np.asarray(energies)
    if np.iscomplexobj(E) or not np.isfinite(E).all():raise ValueError('finite real Fourier energies required')
    H=-mass*S1+angular/r*S2-E[...,None,None]/a*S3
    return 1j/a*H


def local_frame_checks(radii,energies,mass,angular):
    """Replay owned frame connection and pure-radius forcing identities."""
    generator=forcing=unitary=0.
    for rho in radii:
        beta,bp,_=map(float,geometry(rho));a=np.sqrt(beta*beta-1);r=np.sqrt(1+rho*rho)
        R=canonical_map(rho);F=np.linalg.inv(R);v=S2-beta*I2
        connection=(-beta*bp*I2-bp*S2)/(2*a*a)
        M=CommonTimeBulkSplit.local_potential(1.,r,mass,angular)
        forcing=max(forcing,float(np.max(abs(R@np.linalg.solve(v,S3)@F-S2/a))))
        for E in energies:
            G=ks_generator(rho,E,mass,angular)
            pg=np.linalg.solve(v,1j*(E*I2-M)+.5*bp*I2)
            mapped=connection-1j*E*beta/(a*a)*I2+F@G@R
            generator=max(generator,float(np.max(abs(pg-mapped))))
            unitary=max(unitary,float(np.max(abs(G+G.conj().T))))
    return {'PG_KS_generator':generator,'radius_forcing_map':forcing,'KS_antiHermitian_generator':unitary}


def short_response(upstream,energies,mass,angular,initial,what,*,rtol,atol):
    """One coupled solve for four input columns and sixteen Fourier responses."""
    E=np.asarray(energies);A0=np.asarray(initial,complex);W=np.asarray(what,complex)
    if E.shape!=(4,) or A0.shape!=(4,2,3) or W.shape!=(4,4) or upstream<=1:
        raise ValueError('four source fibers, sixteen responses and upstream rho>1 required')
    if not all(np.isfinite(v).all() for v in (E,A0,W)):raise ValueError('finite response data required')
    # This constant callback produces exactly s*chi(s), reusing the frozen
    # supplied radius basis rather than implementing another cutoff.
    normal=CompatibleRadiusDirection(lambda z,n:1. if n==0 else 0.,lambda z,n:0.,.007,.03)
    T1=chart_coordinates(1.)[0]
    def rhs(rho,state):
        A=state[:24].reshape(4,2,3);B=state[24:].reshape(4,4,2,3)
        G=ks_generator(rho,E,mass,angular)
        beta=float(geometry(rho)[0]);a=np.sqrt(beta*beta-1);r=np.sqrt(1+rho*rho)
        s=chart_coordinates(float(rho))[0]-T1
        coefficient=-1j*angular*normal.value(s,0.)/(a*r*r)
        dA=np.einsum('iab,ibc->iac',G,A)
        dB=np.einsum('oab,oibc->oiac',G,B)+coefficient*W[:,:,None,None]*np.einsum('ab,ibc->iac',S2,A)[None,:,:,:]
        return np.concatenate((dA.ravel(),dB.ravel()))
    initial_state=np.concatenate((A0.ravel(),np.zeros(96,complex)))
    solution=solve_ivp(rhs,(float(upstream),1.),initial_state,method='DOP853',rtol=rtol,atol=atol)
    if not solution.success:raise ArithmeticError(solution.message)
    end=solution.y[:,-1]
    return end[:24].reshape(4,2,3),end[24:].reshape(4,4,2,3),{'function_evaluations':int(solution.nfev),'accepted_steps':len(solution.t)-1,'rtol':rtol,'atol':atol}

# Retained executable control: resolution_control
# Folded from scripts/derive_nsc_retarded_response_resolution.py at 329ef809; original receipts retain their producer identities.
import os as _resolution_control_os
import argparse as _resolution_control_argparse
import importlib as _resolution_control_importlib
import importlib.util
import json as _resolution_control_json
from pathlib import Path as _resolution_control_Path
import signal as _resolution_control_signal
import time as _resolution_control_time
import numpy as _resolution_control_np
from recursive_horizons.nsc_posix import ROOT as _resolution_control_ROOT
_resolution_control_SPEC = _resolution_control_importlib.util.spec_from_file_location('frozen_retarded_pilot', _resolution_control_ROOT / 'scripts/derive_nsc_retarded_compatible_response.py')
_resolution_control_P = _resolution_control_importlib.util.module_from_spec(_resolution_control_SPEC)
_resolution_control_SPEC.loader.exec_module(_resolution_control_P)
from recursive_horizons.nsc_transmitting_history_modes import reference_fields_on_grid as _resolution_control_reference_fields_on_grid, characteristic_columns as _resolution_control_characteristic_columns
_resolution_control_OUTPUT = 'results/development/nsc-retarded-response-resolution.json'
_resolution_control_ARTIFACT = 'results/development/artifacts/nsc-retarded-response-resolution.npz'
_resolution_control_HISTORICAL_OWNED = ('scripts/derive_nsc_retarded_response_resolution.py', 'tests/test_nsc_retarded_response_resolution.py', 'docs/nsc-retarded-response-resolution.md')
_resolution_control_OWNED = ('src/recursive_horizons/nsc_retarded_radial_response.py', 'tests/test_nsc_retarded_response_resolution.py', 'docs/nsc-normalized-histories.md')
from .nsc_frozen_sources import selected_source_context as _resolution_control_selected_source_context

def _resolution_control_owned_sources():
    return (_resolution_control_HISTORICAL_OWNED if _resolution_control_selected_source_context() is not None
            else _resolution_control_OWNED)

def _resolution_control_reproducer():
    if _resolution_control_selected_source_context() is not None:
        return 'python3 scripts/derive_nsc_retarded_response_resolution.py --check'
    return "python3 -c 'from recursive_horizons.nsc_retarded_radial_response import resolution_control; resolution_control.check()'"

_resolution_control_EXTRA = (_resolution_control_P.OUTPUT, 'src/recursive_horizons/nsc_common_time_bulk_split.py')
_resolution_control_CPU_CAP = 120.0
_resolution_control_RESAMPLING_TOL = 3e-11
_resolution_control_CAUSAL_TRACE_INTERVAL = (0.05404795482456391, 0.24595204517543606)

def _resolution_control_signature():
    return {**_resolution_control_P.signature(), **{p: _resolution_control_P.digest(p) for p in (*_resolution_control_owned_sources(), *_resolution_control_EXTRA)}}

def _resolution_control_subgrid(x, phi, points):
    if points < 2 or (len(x) - 1) % (points - 1) or phi.shape != (2 * len(x), 12):
        raise ValueError('nested source grid and twelve spin-major columns required')
    stride = (len(x) - 1) // (points - 1)
    return (x[::stride], phi.reshape(2, len(x), 12)[:, ::stride].reshape(2 * points, 12))

def _resolution_control_resampling_error(x, phi, oldx, oldphi):
    xx, ff = _resolution_control_subgrid(x, phi, 801)
    if not _resolution_control_np.array_equal(xx, oldx):
        raise ValueError('new reference grid does not exactly nest archived801')
    return _resolution_control_P.MAXABS(ff - oldphi)

def _resolution_control_one_case(points, steps, xmax, phimax, data, guard):
    x, phi = _resolution_control_subgrid(xmax, phimax, points)
    owner = _resolution_control_P.FourthOrderModePropagator(x, data['mass'].item(), data['angular'].item())
    direction = _resolution_control_P.CompatibleRadiusDirection(_resolution_control_P.axial_bump(guard['axial_center']), lambda z, n: 0.0, 0.007, 0.03)
    provider = _resolution_control_P.CachedZeroRadiusProvider(x, direction)
    times = _resolution_control_np.linspace(0.0, 0.3, steps + 1)
    energy = _resolution_control_np.repeat(data['energies'], 3)
    initial_tangents = _resolution_control_np.zeros((1, *phi.shape), complex)
    right = _resolution_control_np.array([points - 1, 2 * points - 1])
    N, beta, q, r = owner.reference
    sat = -_resolution_control_np.array([N[-1] / q[-1] - beta[-1], -N[-1] / q[-1] - beta[-1]]) / owner.weights[-1]

    def incoming(t):
        forcing = _resolution_control_np.zeros_like(phi)
        forcing[right] = sat[:, None] * phi[right] * _resolution_control_np.exp(-1j * energy * t)[None, :]
        return (forcing, _resolution_control_np.zeros_like(initial_tangents))
    ix = _resolution_control_np.flatnonzero(x == 1.0)
    if len(ix) != 1:
        raise ValueError('exact rho1 node required')
    rows = _resolution_control_np.array([ix[0], points + ix[0]])
    _resolution_control_run = _resolution_control_P.evolve_prepared_field_jets(owner, times, provider, phi, initial_tangents, incoming, sample_rows=rows)
    domain = _resolution_control_P.TransmittingDiracSeamDomain(lapse=N[ix[0]], radial_scale=q[ix[0]], shift=beta[ix[0]], radius=r[ix[0]], surface_id='fixed-rho1-retarded-resolution')
    _, inverse, _ = domain.trace_map(_resolution_control_np.ones(1))
    restriction = inverse[0] @ _resolution_control_P.MODE_TO_CURRENT / (r[ix[0]] * _resolution_control_np.sqrt(q[ix[0]]))
    F = _resolution_control_np.einsum('ab,tbc->tac', restriction, _resolution_control_run['sampled_fields'])
    dF = _resolution_control_np.einsum('ab,tbc->tac', restriction, _resolution_control_run['sampled_tangent_fields'][0])
    rootw = _resolution_control_np.repeat(_resolution_control_np.sqrt(data['weights'] / (2 * _resolution_control_np.pi)), 3)
    weightedF = (F * rootw).reshape(-1, 12)
    weighteddF = (dF * rootw).reshape(-1, 12)
    kernels = {}
    for name, key in (('covariance', 'source'), ('projector', 'projector')):
        C = _resolution_control_P.block_diag(*data[key])
        kernels[name + '_tangent'] = _resolution_control_P.restriction_covariance_tangent(weightedF, weighteddF, C, _resolution_control_np.zeros_like(C))
    diagnostics = {k: _resolution_control_run[k] for k in ('norm_flux_residual', 'metric_flux_tangent_residual', 'auxiliary_identity_residual', 'tangent_boundary_trace')}
    diagnostics.update(initial_tangent_norm=_resolution_control_P.MAXABS(initial_tangents), incident_tangent_norm=0.0, source_covariance_tangent_norm=0.0, initial_metric_change=_resolution_control_P.MAXABS(provider.metric - owner.reference), source_closed_column_norm=_resolution_control_P.MAXABS(phi[:, 2::3]), source_covariance_closed_column_norm=_resolution_control_P.MAXABS(data['source'][:, :, 2]), source_covariance_hermiticity=_resolution_control_P.MAXABS(data['source'] - data['source'].swapaxes(-1, -2).conj()), frame_map_parameter_derivative=0.0, stationary_phase_reapplied=False)
    return {'x': x, 'times': times, 'sample_rows': rows, 'field': _resolution_control_run['evolved_field'], 'field_tangent': _resolution_control_run['tangent_fields'][0], 'trace': F, 'trace_tangent': dF, 'characteristic_trace_tangent': _resolution_control_run['sampled_tangent_fields'][0], 'restriction': restriction, **kernels, 'diagnostics_json': _resolution_control_np.frombuffer(_resolution_control_json.dumps(diagnostics, sort_keys=True).encode(), _resolution_control_np.uint8)}

def _resolution_control_compare(c, f, kind, limit):
    rows = {}
    for name, key in (('field', 'field_tangent'), ('trace', 'trace_tangent'), ('covariance', 'covariance_tangent')):
        coarse = c[key]
        fine = f[key]
        if kind == 'spatial' and name == 'field':
            fine = _resolution_control_subgrid(f['x'], fine, len(c['x']))[1]
        if kind == 'temporal':
            stride = (len(f['times']) - 1) // (len(c['times']) - 1)
            if not _resolution_control_np.array_equal(c['times'], f['times'][::stride]):
                raise ValueError('time nodes are not nested')
            if name == 'trace':
                fine = fine[::stride]
            if name == 'covariance':
                fine = fine.reshape(len(f['times']), 2, len(f['times']), 2)[::stride, :, ::stride, :].reshape(coarse.shape)
        error = _resolution_control_P.MAXABS(coarse - fine)
        norm = _resolution_control_P.MAXABS(fine)
        ratio = error / norm if norm else None
        rows[name] = {'max_absolute_difference': error, 'fine_compared_response_norm': norm, 'difference_over_fine_norm': ratio, 'relative_indicator_limit': limit, 'passes': ratio is not None and ratio < limit}
    return rows

def _resolution_control_case_arrays(arrays, name):
    return {k[len(name) + 1:]: v for k, v in arrays.items() if k.startswith(name + '/')}

def _resolution_control_peak(values, times):
    magnitudes = abs(values)
    index = _resolution_control_np.unravel_index(_resolution_control_np.argmax(magnitudes), magnitudes.shape)
    return {'magnitude': float(magnitudes[index]), 'PG_time': float(times[index[0]]), 'source_column': int(index[-1])}

def _resolution_control_outside_causal_support(values, times):
    outside = (times < _resolution_control_CAUSAL_TRACE_INTERVAL[0]) | (times > _resolution_control_CAUSAL_TRACE_INTERVAL[1])
    weights = _resolution_control_np.ones(len(times))
    weights[[0, -1]] = 0.5
    energy = _resolution_control_np.sum(abs(values) ** 2, axis=tuple(range(1, values.ndim)))
    total = float(weights @ energy)
    outer = float(weights[outside] @ energy[outside])
    return {'maximum_entry': _resolution_control_P.MAXABS(values[outside]), 'sampled_energy_fraction': outer / total if total else 0.0, 'interval_estimate': list(_resolution_control_CAUSAL_TRACE_INTERVAL), 'directed_bound': False, 'acceptance_gate': False}

def _resolution_control_analyze(arrays, meta, oldarrays):
    summaries = {}
    comparisons = {}
    historical = {}
    for case in meta['completed']:
        a = _resolution_control_case_arrays(arrays, case)
        x, phi = _resolution_control_subgrid(arrays['reference/x'], arrays['reference/phi'], len(a['x']))
        E = _resolution_control_np.repeat(arrays['source/energies'], 3)
        t = a['times']
        stationary = phi * _resolution_control_np.exp(-1j * E * t[-1])[None, :]
        trace0 = a['restriction'] @ phi[a['sample_rows']]
        stationary_trace = trace0[None, :, :] * _resolution_control_np.exp(-1j * t[:, None, None] * E[None, None, :])
        late = t > meta['support']['PG_time_support'][1]
        slow = a['characteristic_trace_tangent'][:, 0, :]
        summaries[case] = {'response_max_norms': {n: _resolution_control_P.MAXABS(a[k]) for n, k in (('field', 'field_tangent'), ('trace', 'trace_tangent'), ('covariance', 'covariance_tangent'))}, 'covariance_hermiticity_residual': _resolution_control_P.MAXABS(a['covariance_tangent'] - a['covariance_tangent'].conj().T), 'sampled_projector_kernel_tangent_norm': _resolution_control_P.MAXABS(a['projector_tangent']), 'sampled_projector_kernel_hermiticity_residual': _resolution_control_P.MAXABS(a['projector_tangent'] - a['projector_tangent'].conj().T), 'stationary_reference_field_drift': _resolution_control_P.MAXABS(a['field'] - stationary), 'stationary_reference_trace_drift': _resolution_control_P.MAXABS(a['trace'] - stationary_trace), 'slow_characteristic_trace_peak': _resolution_control_peak(slow, t), 'slow_characteristic_trace_after_pulse': _resolution_control_peak(slow[late], t[late]), 'late_window_start': meta['support']['PG_time_support'][1], 'characteristic_trace_outside_causal_support': _resolution_control_outside_causal_support(a['characteristic_trace_tangent'], t), 'slow_trace_outside_causal_support': _resolution_control_outside_causal_support(slow, t), 'diagnostics': _resolution_control_json.loads(a['diagnostics_json'].tobytes())}
    if '3201_128' in meta['completed']:
        comparisons['spatial'] = _resolution_control_compare(_resolution_control_case_arrays(arrays, '1601_128'), _resolution_control_case_arrays(arrays, '3201_128'), 'spatial', 0.25)
    if '3201_256' in meta['completed']:
        comparisons['temporal'] = _resolution_control_compare(_resolution_control_case_arrays(arrays, '3201_128'), _resolution_control_case_arrays(arrays, '3201_256'), 'temporal', 0.1)
    for case in ('1601_128', '3201_128'):
        if case in meta['completed']:
            historical['801_128_to_' + case] = _resolution_control_compare(_resolution_control_case_arrays(oldarrays, '801_128'), _resolution_control_case_arrays(arrays, case), 'spatial', 0.25)
    algebra = max((r['diagnostics'][k] for r in summaries.values() for k in ('norm_flux_residual', 'metric_flux_tangent_residual', 'auxiliary_identity_residual')), default=0.0)
    herm = max((max(r['covariance_hermiticity_residual'], r['sampled_projector_kernel_hermiticity_residual']) for r in summaries.values()), default=0.0)
    passed = len(comparisons) == 2 and all((v['passes'] for comp in comparisons.values() for v in comp.values())) and (algebra < 3e-11) and (herm < 3e-11) and (meta.get('resampling_residual', float('inf')) < _resolution_control_RESAMPLING_TOL) and (not meta['cpu_budget_exceeded'])
    return (summaries, comparisons, historical, passed, algebra, herm)

def _resolution_control_record_from_arrays(arrays, meta, digest, size, oldarrays):
    summaries, comparisons, historical, passed, algebra, herm = _resolution_control_analyze(arrays, meta, oldarrays)
    return {'schema': 'NSC-RETARDED-RESPONSE-RESOLUTION-v1', 'accountable_author': 'Douglas Ek', 'status': 'PASS: sampled resolution indicators only' if passed else 'OPEN: targeted sampled resolution unresolved', 'control': {'family': '14_1', 'amplitude': 0.0, 'cases_completed': meta['completed'], 'support': meta.get('support'), 'resampling_tolerance': _resolution_control_RESAMPLING_TOL, 'spatial_indicator_limit': 0.25, 'temporal_indicator_limit': 0.1, 'algebra_tolerance': 3e-11, 'reference_sampling_calls': meta['reference_sampling_calls'], 'source_law': 'unchanged original signed frequencies, weights and coherent source fibers'}, 'residuals': {'archived801_resampling': meta.get('resampling_residual'), 'maximum_flux_and_auxiliary_algebra': algebra, 'maximum_kernel_hermiticity': herm}, 'per_case': summaries, 'comparison': comparisons, 'historical_context_only': historical, 'runtime': {k: meta[k] for k in ('CPU_seconds', 'wall_seconds', 'cpu_budget_seconds', 'cpu_budget_exceeded', 'case_CPU_seconds', 'resampling_CPU_seconds', 'stop_reason')}, 'scope': {'continuum_error_bound': None, 'full_continuous_C0_match': 'OPEN', 'full_CAR_verified': False, 'sampled_projector_tangent_reported': True, 'stationary_drift_subtracted': False, 'new_horizon_scattering_source_state_solve': False, 'metric_evolution': False, 'physical_constraint_solution': False, 'same_pulse_and_source_as_pilot': True, 'differences_are_indicators': True}, 'source_hashes': {p: meta['signature'][p] for p in _resolution_control_owned_sources()}, 'input_hashes': {p: h for p, h in meta['signature'].items() if p not in _resolution_control_owned_sources()}, 'reused_artifacts': meta['reused_artifacts'], 'payload': {'path': _resolution_control_ARTIFACT, 'sha256': digest, 'bytes': size}, 'reproducer': _resolution_control_reproducer()}

def _resolution_control_publish(arrays, meta, oldarrays):
    encoded = {**arrays, 'metadata_json': _resolution_control_np.frombuffer(_resolution_control_json.dumps(meta, sort_keys=True).encode(), _resolution_control_np.uint8)}
    raw = _resolution_control_P.deterministic_npz_bytes(encoded)
    (_resolution_control_ROOT / _resolution_control_ARTIFACT).write_bytes(raw)
    record = _resolution_control_record_from_arrays(encoded, meta, _resolution_control_P.sha256(raw).hexdigest(), len(raw), oldarrays)
    (_resolution_control_ROOT / _resolution_control_OUTPUT).write_text(_resolution_control_json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + '\n')
    return record

def _resolution_control_run():
    start = _resolution_control_time.process_time()
    wall = _resolution_control_time.monotonic()
    before = _resolution_control_signature()
    arrays = {}
    oldarrays = {}
    meta = {'signature': before, 'completed': [], 'case_CPU_seconds': {}, 'cpu_budget_seconds': _resolution_control_CPU_CAP, 'cpu_budget_exceeded': False, 'reference_sampling_calls': 0, 'resampling_CPU_seconds': None, 'reused_artifacts': [], 'stop_reason': 'checkpoint'}

    def stop(signum, frame):
        raise _resolution_control_P.CPUExceeded('targeted resolution CPU budget exhausted')
    previous = _resolution_control_signal.signal(_resolution_control_signal.SIGPROF, stop)
    _resolution_control_signal.setitimer(_resolution_control_signal.ITIMER_PROF, _resolution_control_CPU_CAP)

    def checkpoint(reason):
        meta.update(CPU_seconds=_resolution_control_time.process_time() - start, wall_seconds=_resolution_control_time.monotonic() - wall, stop_reason=reason)
        if _resolution_control_signature() != before:
            raise ValueError('resolution owner/input changed during bounded run')
        return _resolution_control_publish(arrays, meta, oldarrays)
    try:
        pilot = _resolution_control_P.check()
        oldx, oldphi, data, records = _resolution_control_P.inputs()
        meta['support'] = _resolution_control_P.support_guard()
        with _resolution_control_np.load(_resolution_control_ROOT / pilot['payload']['path'], allow_pickle=False) as a:
            oldarrays = {k: a[k] for k in a.files if k.startswith('801_128/')}
        with _resolution_control_np.load(_resolution_control_ROOT / records[1]['payload']['path'], allow_pickle=False) as a:
            at_zero = a['14_1/at_zero']
        meta['reused_artifacts'] = [pilot['payload'], *[r['payload'] for r in records]]
        arrays.update({'source/' + k: v for k, v in data.items()})
        arrays['source/at_zero'] = at_zero
        x = _resolution_control_np.linspace(-2.0, 1.2, 3201)
        begin = _resolution_control_time.process_time()
        meta['reference_sampling_calls'] = 1
        phi = _resolution_control_characteristic_columns(_resolution_control_reference_fields_on_grid(data['energies'], data['mass'].item(), data['angular'].item(), at_zero, x))
        meta['resampling_CPU_seconds'] = _resolution_control_time.process_time() - begin
        arrays.update({'reference/x': x, 'reference/phi': phi})
        meta['resampling_residual'] = _resolution_control_resampling_error(x, phi, oldx, oldphi)
        if meta['resampling_residual'] >= _resolution_control_RESAMPLING_TOL:
            return checkpoint('OPEN: archived reference mismatch; no propagation')
        print('reference resampling matched: ' + str(meta['resampling_residual']), flush=True)
        for points, steps in ((1601, 128), (3201, 128), (3201, 256)):
            if steps == 256:
                spatial = _resolution_control_compare(_resolution_control_case_arrays(arrays, '1601_128'), _resolution_control_case_arrays(arrays, '3201_128'), 'spatial', 0.25)
                if not all((v['passes'] for v in spatial.values())):
                    return checkpoint('OPEN: spatial indicator failed; temporal solve not authorized')
            case = f'{points}_{steps}'
            begin = _resolution_control_time.process_time()
            result = _resolution_control_one_case(points, steps, x, phi, data, meta['support'])
            arrays.update({case + '/' + k: v for k, v in result.items()})
            meta['completed'].append(case)
            meta['case_CPU_seconds'][case] = _resolution_control_time.process_time() - begin
            checkpoint('completed ' + case)
            print(case + ' checkpoint CPU ' + str(meta['case_CPU_seconds'][case]), flush=True)
        meta['stop_reason'] = 'completed all three prescribed cases'
    except _resolution_control_P.CPUExceeded:
        meta['cpu_budget_exceeded'] = True
        meta['stop_reason'] = 'OPEN: CPU budget exhausted; no automatic extension'
    finally:
        _resolution_control_signal.setitimer(_resolution_control_signal.ITIMER_PROF, 0)
        _resolution_control_signal.signal(_resolution_control_signal.SIGPROF, previous)
    return checkpoint(meta['stop_reason'])

def _resolution_control_energy_fraction_roundoff(values, times, stored, computed):
    """Relative reduction error for a positive sampled-energy ratio.

    The numerator and denominator sum the same nonnegative |value|**2
    entries, with positive endpoint weights.  gamma_k covers complex
    modulus/square, component summation and time reduction; no cancellation
    or absolute near-zero allowance is introduced.  Propagate the two
    positive relative errors through the quotient, then bound two replays.
    This diagnostic is not an acceptance gate or a continuum error bound.
    """
    np = _resolution_control_np
    components = int(np.prod(values.shape[1:]))
    eps = np.finfo(np.float64).eps  # Conservative: twice unit roundoff.
    operations = 8 * components + 2 * len(times) + 16
    gamma = operations * eps / (1 - operations * eps)
    quotient_relative = (2 * gamma + eps * (1 + gamma)) / (1 - gamma)
    if not all(type(v) is float and np.isfinite(v) and 0 <= v <= 1 for v in (stored, computed)):
        raise ValueError('invalid sampled-energy fraction')
    bound = 2 * quotient_relative * max(stored, computed) / (1 - quotient_relative)
    return float(np.nextafter(bound, np.inf)) if bound else 0.0


def _resolution_control_validate_replay(record, replay, arrays, meta):
    from copy import deepcopy
    from .nsc_weighted_receipt_replay import exact_record_equal
    strict = deepcopy(replay)
    for case in meta['completed']:
        a = _resolution_control_case_arrays(arrays, case)
        for field, values in (
            ('characteristic_trace_outside_causal_support', a['characteristic_trace_tangent']),
            ('slow_trace_outside_causal_support', a['characteristic_trace_tangent'][:, 0, :]),
        ):
            stored = record['per_case'][case][field]['sampled_energy_fraction']
            computed = replay['per_case'][case][field]['sampled_energy_fraction']
            if type(stored) is not float or type(computed) is not float:
                raise ValueError('sampled-energy fraction numeric type changed')
            bound = _resolution_control_energy_fraction_roundoff(values, a['times'], stored, computed)
            if abs(stored - computed) > bound:
                raise ValueError('sampled-energy fraction exceeds replay roundoff envelope')
            strict['per_case'][case][field]['sampled_energy_fraction'] = stored
    if not exact_record_equal(strict,record):
        raise ValueError('resolution receipt differs from saved arrays')


def _resolution_control_check():
    record = _resolution_control_json.loads((_resolution_control_ROOT / _resolution_control_OUTPUT).read_text())
    _resolution_control_P.check()
    if _resolution_control_P.digest(_resolution_control_ARTIFACT) != record['payload']['sha256']:
        raise ValueError('resolution artifact changed')
    with _resolution_control_np.load(_resolution_control_ROOT / _resolution_control_ARTIFACT, allow_pickle=False) as a:
        arrays = {k: a[k] for k in a.files}
    meta = _resolution_control_json.loads(arrays['metadata_json'].tobytes())
    if _resolution_control_signature() != meta['signature']:
        raise ValueError('resolution producers/inputs changed')
    oldx, oldphi, data, records = _resolution_control_P.inputs()
    if 'reference/phi' in arrays:
        residual = _resolution_control_resampling_error(arrays['reference/x'], arrays['reference/phi'], oldx, oldphi)
        if residual != meta['resampling_residual']:
            raise ValueError('reference resampling residual changed')
    for k, v in data.items():
        if not _resolution_control_np.array_equal(arrays['source/' + k], v):
            raise ValueError('original source data changed')
    with _resolution_control_np.load(_resolution_control_ROOT / records[1]['payload']['path'], allow_pickle=False) as a:
        if not _resolution_control_np.array_equal(arrays['source/at_zero'], a['14_1/at_zero']):
            raise ValueError('reference initial columns changed')
    pilot = _resolution_control_json.loads((_resolution_control_ROOT / _resolution_control_P.OUTPUT).read_text())
    with _resolution_control_np.load(_resolution_control_ROOT / pilot['payload']['path'], allow_pickle=False) as a:
        oldarrays = {k: a[k] for k in a.files if k.startswith('801_128/')}
    replay = _resolution_control_record_from_arrays(arrays, meta, _resolution_control_P.digest(_resolution_control_ARTIFACT), (_resolution_control_ROOT / _resolution_control_ARTIFACT).stat().st_size, oldarrays)
    _resolution_control_validate_replay(record, replay, arrays, meta)
    return record

from types import SimpleNamespace as _resolution_control_Namespace
resolution_control = _resolution_control_Namespace(
    ARTIFACT=_resolution_control_ARTIFACT,
    CAUSAL_TRACE_INTERVAL=_resolution_control_CAUSAL_TRACE_INTERVAL,
    CPU_CAP=_resolution_control_CPU_CAP,
    EXTRA=_resolution_control_EXTRA,
    OUTPUT=_resolution_control_OUTPUT,
    OWNED=_resolution_control_OWNED,
    P=_resolution_control_P,
    Path=_resolution_control_Path,
    RESAMPLING_TOL=_resolution_control_RESAMPLING_TOL,
    ROOT=_resolution_control_ROOT,
    SPEC=_resolution_control_SPEC,
    analyze=_resolution_control_analyze,
    argparse=_resolution_control_argparse,
    case_arrays=_resolution_control_case_arrays,
    characteristic_columns=_resolution_control_characteristic_columns,
    check=_resolution_control_check,
    energy_fraction_roundoff=_resolution_control_energy_fraction_roundoff,
    compare=_resolution_control_compare,
    importlib=_resolution_control_importlib,
    json=_resolution_control_json,
    np=_resolution_control_np,
    one_case=_resolution_control_one_case,
    os=_resolution_control_os,
    outside_causal_support=_resolution_control_outside_causal_support,
    peak=_resolution_control_peak,
    publish=_resolution_control_publish,
    record_from_arrays=_resolution_control_record_from_arrays,
    reference_fields_on_grid=_resolution_control_reference_fields_on_grid,
    resampling_error=_resolution_control_resampling_error,
    run=_resolution_control_run,
    signal=_resolution_control_signal,
    signature=_resolution_control_signature,
    subgrid=_resolution_control_subgrid,
    time=_resolution_control_time,
    validate_replay=_resolution_control_validate_replay,
)
