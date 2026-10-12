"""Bounded group14 +/- incoming-source pilot with exact source-keyed caches.

This extends the successful group13 LOCAL contraction to the actual mixed
mass/angular vertices. Only missing smooth basis and contour reflection nodes
are evaluated, with the existing physical owners. No historical generator,
history propagator, metric step or additional state is introduced.
"""
from hashlib import sha256
import inspect
import json
from pathlib import Path
import tempfile
from time import perf_counter

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from .nsc_common_ks_trace import restrict_resolved_modes
from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_complex_horizon_modes import complex_reflection
from .nsc_incoming_subgap_source import canonical_analytic_basis, COARSE_PANELS
from .nsc_incoming_state_moments import AXIAL, RADIUS, incoming_adiabatic_reference, incoming_group_factor, incoming_state_moments
from .nsc_lorentzian import geometry
from .nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from .nsc_paired_horizon_preparation import PairedHorizonSeedMap
from .nsc_subgap_history_response import inner_horizon_basis, AnalyticResponsePanel, source_functions
from .nsc_transmitting_dirac_domain import I2, S1, S2, S3


MASS = np.pi/2
MAGNITUDE = np.sqrt(5.)
STAGES = {'base': ((24, 3),), 'refined': ((24, 3), (48, 3)),
          'height': ((24, 3), (48, 3), (48, 4))}
OWNER_PATHS = (
    'src/recursive_horizons/nsc_complex_horizon_modes.py',
    'src/recursive_horizons/nsc_paired_horizon_preparation.py',
    'src/recursive_horizons/nsc_pg_massive_modes.py',
    'src/recursive_horizons/nsc_subgap_history_response.py',
    'src/recursive_horizons/nsc_common_time_bulk_split.py',
    'src/recursive_horizons/nsc_lorentzian.py',
    'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
    'src/recursive_horizons/nsc_unruh_state.py',
    'src/recursive_horizons/nsc_gauge_source.py',
)


def _hash(path): return sha256(Path(path).read_bytes()).hexdigest()
def _load(path):
    with np.load(path, allow_pickle=False) as data: return {k: data[k].copy() for k in data.files}
def _read(root, path): return json.loads((root/path).read_text())


def _authenticated(root, item):
    if _hash(root/item['path']) != item['sha256']: raise ValueError('physical source input changed')
    return _load(root/item['path'])


def recover_unsewn(physical, reflection):
    """Recover (v,w) from authenticated closed-channel (R v,w,0) columns."""
    F, R = np.asarray(physical, complex), np.asarray(reflection, complex)
    if F.shape[-2:] != (2, 3) or F.shape[:-2] != R.shape or not np.isfinite(F).all() or not np.isfinite(R).all():
        raise ValueError('frequency-resolved physical closed-channel columns and reflection required')
    if np.max(abs(abs(R)-1)) > 3e-10 or np.max(abs(F[..., 2])) > 3e-11:
        raise ValueError('retained subgap must have unit reflection and no open infinity column')
    return np.stack((F[..., 0]/R[..., None], F[..., 1]), axis=-1)


def mixed_vertices(energy, angular):
    p = -complex(energy)/AXIAL
    return np.array([-MASS*S1+angular/RADIUS*S2+p*S3, p*S3,
                     -p*I2, angular/(2*RADIUS)*S2])


def mixed_bilinears(energy, angular, upper, dual):
    return np.array([dual.conj().T@V@upper for V in mixed_vertices(energy, angular)])


def mixed_reference_remainder(energies, angular):
    if np.iscomplexobj(energies): raise ValueError('reference subtraction requires real energies')
    P, _ = incoming_adiabatic_reference(energies, MASS, angular)
    V = np.array([mixed_vertices(e, angular) for e in energies])
    return (.5*np.trace(V, axis1=-2, axis2=-1)-np.einsum('nij,nvji->nv', P, V)).real


def _continue_pair(z, angular, upper, dual):
    energies = [complex(z)] if not complex(z).imag else [complex(z), complex(z).conjugate()]
    initial = np.array([upper] if len(energies) == 1 else [upper, dual])
    def rhs(rho, value):
        beta, bp, _ = geometry(rho)
        potential = CommonTimeBulkSplit.local_potential(1., np.sqrt(1+rho*rho), MASS, angular)
        columns = value.reshape(len(energies), 2, 2)
        return np.array([np.linalg.solve(S2-beta*I2, 1j*(E*I2-potential)+.5*bp*I2)@F
                         for E, F in zip(energies, columns)]).ravel()
    run = solve_ivp(rhs, (0., 1.), initial.ravel(), method='DOP853', rtol=2e-13, atol=2e-15)
    if not run.success: raise ArithmeticError(run.message)
    end = run.y[:, -1].reshape(len(energies), 2, 2)
    U, D = end[0], end[0] if len(energies) == 1 else end[1]
    velocity = lambda rho: S2-float(geometry(rho)[0])*I2
    residual = float(np.max(abs(D.conj().T@velocity(1.)@U-dual.conj().T@velocity(0.)@upper)))
    return U, D, residual, int(run.nfev)


def _preparation(config):
    return PairedHorizonSeedMap(*(config[key] for key in
        ('horizon_rho', 'surface_gravity', 'omega', 'horizon_offset', 'scattering_tolerance')))


def _basis_node(z, sign, config):
    start = perf_counter(); p = _preparation(config); angular = sign*MAGNITUDE
    U0 = inner_horizon_basis(z, MASS, angular, p)
    D0 = U0 if not complex(z).imag else inner_horizon_basis(complex(z).conjugate(), MASS, angular, p)
    U, D, residual, evaluations = _continue_pair(z, angular, U0, D0)
    return {'frequency': np.array(z), 'at_zero': U0, 'at_zero_dual': D0,
            'at_one': U, 'at_one_dual': D, 'cross_current': np.array(residual),
            'inner_basis_solves': np.array(1 if not complex(z).imag else 2),
            'short_radial_solves': np.array(1), 'short_radial_function_evaluations': np.array(evaluations),
            'seconds': np.array(perf_counter()-start)}


def _reflection_node(z, sign, config):
    start = perf_counter(); p = _preparation(config)
    value = complex_reflection(z, MASS, sign*MAGNITUDE, p.background, radial_collar=p.collar_delta_q)
    return {'frequency': np.array(z), 'reflection': np.array(value.reflection),
            'diagnostics_json': np.frombuffer(json.dumps(value.residuals, sort_keys=True).encode(), np.uint8),
            'seconds': np.array(perf_counter()-start)}


class PilotBudgetReached(RuntimeError): pass


class ExactNodeCache:
    """Small replayable numerical cache, bound to solver sources and inputs."""
    def __init__(self, root, config, budget_seconds):
        if not np.isfinite(budget_seconds) or not 0 < budget_seconds <= 600:
            raise ValueError('explicit pilot stage budget in (0,600] seconds required')
        self.signature = {'owners': {p: _hash(root/p) for p in OWNER_PATHS}, 'config': config,
            'mass': MASS, 'angular_magnitude': MAGNITUDE,
            'functions': {f.__name__: sha256(inspect.getsource(f).encode()).hexdigest()
                          for f in (_preparation, _continue_pair, _basis_node, _reflection_node)}}
        tag = sha256(json.dumps(self.signature, sort_keys=True).encode()).hexdigest()[:20]
        self.directory = Path(tempfile.gettempdir())/('nsc-incoming-subgap-completion-'+tag)
        self.directory.mkdir(exist_ok=True)
        self.config, self.start, self.budget = config, perf_counter(), budget_seconds
        self.cost = {'new_basis_nodes': 0, 'new_reflection_nodes': 0, 'cache_hits': 0,
                     'inner_basis_solves': 0, 'short_radial_solves': 0,
                     'new_reflection_seconds': 0., 'new_basis_seconds': 0., 'new_reflection_function_evaluations': 0}

    def node(self, kind, z, sign):
        z = complex(z)
        key = sha256((kind+','+str(sign)+','+z.real.hex()+','+z.imag.hex()).encode()).hexdigest()[:24]
        path = self.directory/(key+'.npz')
        if path.exists():
            row = _load(path)
            if row['frequency'].item() != z or json.loads(row['signature_json'].tobytes()) != self.signature:
                raise ValueError('cached node source/configuration differs')
            self.cost['cache_hits'] += 1
            return row
        if perf_counter()-self.start > self.budget:
            raise PilotBudgetReached(json.dumps({'stage_seconds': perf_counter()-self.start, 'cache': str(self.directory),
                                                 'next_kind': kind, 'next_sign': sign, 'next_frequency': [z.real, z.imag], 'cost': self.cost}))
        if kind == 'reflection':
            if self.cost['new_reflection_nodes'] >= 720: raise PilotBudgetReached('720 new-reflection pilot ceiling reached')
            row = _reflection_node(z, sign, self.config)
            self.cost['new_reflection_nodes'] += 1; self.cost['new_reflection_seconds'] += float(row['seconds'])
            self.cost['new_reflection_function_evaluations'] += json.loads(row['diagnostics_json'].tobytes())['function_evaluations']
            if self.cost['new_reflection_nodes'] % 24 == 0:
                print(f"group14 pilot: {self.cost['new_reflection_nodes']} new exact reflections, {perf_counter()-self.start:.1f}s", flush=True)
        elif kind == 'basis':
            row = _basis_node(z, sign, self.config)
            self.cost['new_basis_nodes'] += 1; self.cost['new_basis_seconds'] += float(row['seconds'])
            self.cost['inner_basis_solves'] += int(row['inner_basis_solves']); self.cost['short_radial_solves'] += 1
        else: raise ValueError('only basis and reflection cache nodes are supported')
        row['signature_json'] = np.frombuffer(json.dumps(self.signature, sort_keys=True).encode(), np.uint8)
        path.write_bytes(deterministic_npz_bytes(row))
        return row


def _contour_nodes(kappa, points, fraction):
    x, _ = leggauss(points); h = kappa/fraction
    return [a+(xi+1)*(b-a)/2 for a, b in ((1.+0j, 1.+1j*h), (1.+1j*h, MASS+1j*h), (MASS+1j*h, MASS+0j)) for xi in x]


def integrate_mixed_panel(panel, angular, kappa, reflections, points, fraction):
    x, w = leggauss(points); E = 1+(x+1)*(MASS-1)/2; weights = w*(MASS-1)/2
    smooth = np.zeros(4, complex)
    for e, weight in zip(E, weights):
        d, _ = source_functions(e, kappa); B = panel(e)
        smooth += weight*d*(B[:, 0, 0]-B[:, 1, 1])
    reference = weights@mixed_reference_remainder(E, angular)
    reflected = np.zeros(4, complex); h = kappa/fraction
    for a, b in ((1.+0j, 1.+1j*h), (1.+1j*h, MASS+1j*h), (MASS+1j*h, MASS+0j)):
        for xi, wi in zip(x, w):
            z = a+(xi+1)*(b-a)/2
            if z not in reflections: raise ValueError('exact requested reflection node absent')
            _, s = source_functions(z, kappa)
            reflected += wi*(b-a)/2*(-1j*s*reflections[z]*panel(z)[:, 1, 0])
    return smooth.real+reference+2*reflected.real


def run_group14_pilot(root, *, stage='base', budget_seconds=300):
    """Stage the smallest pilot, then reuse it for requested validation only.

    base:9 smooth nodes,24-point contour. refined:17 nodes and48-point
    contour. height:the second48-point contour height. Both angular signs
    are always retained; other groups are deliberately outside this API.
    """
    if stage not in STAGES: raise ValueError('explicit base/refined/height pilot stage required')
    root = Path(root); clock = perf_counter()
    state = _read(root, 'results/development/nsc-mode-resolved-cauchy-state.json'); channel = state['channels'][14]
    restart = _read(root, 'results/development/nsc-compact-matched-restart.json')
    mode_record = _read(root, 'results/development/nsc-pg-massive-mode-resolution.json'); modes = _authenticated(root, mode_record['payload'])
    if channel['compact_mass'] != MASS or channel['angular_eigenvalue'] != MAGNITUDE:
        raise ValueError('group14 retained mass/angular labels changed')
    if _hash(root/COARSE_PANELS) not in Path(COARSE_PANELS).name: raise ValueError('coarse source payload changed')
    old_source = _load(root/COARSE_PANELS)
    cfg = restart['scattering_provenance']['config']; cache = ExactNodeCache(root, cfg, budget_seconds)
    allE = (1+MASS)/2+(MASS-1)/2*np.cos(np.pi*np.arange(17)/16)
    E = allE[::2] if stage == 'base' else allE
    zcheck = (1+MASS)/2+1j*cfg['surface_gravity']/3
    arrays = {}; per_sign = {}; factor = incoming_group_factor(channel, (-1, 1)); cumulative = {}
    for sign in (-1, 1):
        angular = sign*MAGNITUDE; prefix = str(sign)
        label = modes['label']; ix = np.flatnonzero((label[:, 0] == 14)&(label[:, 1] == sign)&(label[:, 5] > 1)&(label[:, 5] < MASS))
        if len(ix) != 8: raise ValueError('eight physical subgap source rows per sign required')
        E8, w8, R8 = label[ix, 5], label[ix, 6], modes['reflection'][ix, 0]
        recovered0 = recover_unsewn(modes['interior'][ix, 1], R8)
        recovered1 = recover_unsewn(modes['interior'][ix, 2], R8)
        oldF, chart = restrict_resolved_modes(1., E8, np.concatenate((recovered1, np.zeros((8, 2, 1))), axis=2))
        oldF = oldF[..., :2]; S = chart['clock_shift']
        oldB = np.array([mixed_bilinears(e, angular, F, F) for e, F in zip(E8, oldF)])
        formula8 = mixed_reference_remainder(E8, angular)
        for i, (e, B, R) in enumerate(zip(E8, oldB, R8)):
            d, s = source_functions(e, cfg['surface_gravity'])
            formula8[i] += (d*(B[:, 0, 0]-B[:, 1, 1])+2*np.real(-1j*s*R*B[:, 1, 0])).real
        direct = incoming_state_moments(E8, modes['interior'][ix, 2], modes['source_covariance'][ix], modes['source_projector'][ix], MASS, angular)
        coarse = old_source[f'subgap8/14_{sign}/kernels']
        rows = [cache.node('basis', complex(e), sign) for e in E]
        check = cache.node('basis', zcheck, sign)
        F = np.array([canonical_analytic_basis(e, row['at_one'], S) for e, row in zip(E, rows)])
        B = np.array([mixed_bilinears(e, angular, f, f) for e, f in zip(E, F)])
        Fc = canonical_analytic_basis(zcheck, check['at_one'], S)
        Fd = canonical_analytic_basis(zcheck.conjugate(), check['at_one_dual'], S)
        Bc = mixed_bilinears(zcheck, angular, Fc, Fd)
        frequencies = list(dict.fromkeys(z for points, fraction in STAGES[stage] for z in _contour_nodes(cfg['surface_gravity'], points, fraction)))
        rr = [cache.node('reflection', z, sign) for z in frequencies]
        reflection = {z: row['reflection'].item() for z, row in zip(frequencies, rr)}
        panels = {'9': AnalyticResponsePanel(E, B, 1., MASS)} if stage == 'base' else {
            '17': AnalyticResponsePanel(E, B, 1., MASS), '9': AnalyticResponsePanel(E[::2], B[::2], 1., MASS)}
        tables = {f'n{n}_q{points}_h{fraction}': integrate_mixed_panel(panel, angular, cfg['surface_gravity'], reflection, points, fraction)
                  for n, panel in panels.items() for points, fraction in STAGES[stage]}
        selected = 'n9_q24_h3' if stage == 'base' else 'n17_q48_h3'; fine = tables[selected]
        panel = panels['9' if stage == 'base' else '17']
        diag = {
            'recovered_zero_vs_archived_fundamental': float(np.max(abs(recovered0-modes['interior_fundamental'][ix, 1][:, :, [1, 0]]))),
            'recovered_one_vs_archived_fundamental': float(np.max(abs(recovered1-modes['interior_fundamental'][ix, 2][:, :, [1, 0]]))),
            'old_real_basis_unitarity': float(np.max(abs(oldF.swapaxes(-1, -2).conj()@oldF-I2))),
            'new_real_basis_unitarity': float(np.max(abs(F.swapaxes(-1, -2).conj()@F-I2))),
            'new_cross_current': max(float(row['cross_current']) for row in [*rows, check]),
            'complex_cross_Gram': float(np.max(abs(Fd.conj().T@Fc-I2))),
            'complex_bilinear_interpolation': float(np.max(abs(panel(zcheck)-Bc))),
            'real8_bilinear_interpolation': float(np.max(abs(np.array([panel(e) for e in E8])-oldB))),
            'old8_formula_vs_direct': float(np.max(abs(formula8-direct['kernels']))),
            'old8_formula_vs_archived': float(np.max(abs(formula8-coarse))),
            'T01_coherence_overlap': float(np.max(abs(B[:, 2, 1, 0]))),
            'reflection_outer_residual': max(json.loads(row['diagnostics_json'].tobytes())['outer_Riccati_residual'] for row in rr),
        }
        per_sign[prefix] = {'raw_integral': fine.tolist(), 'physical_contribution': (factor*fine).tolist(),
                            'coarse8_physical_contribution': (factor*(w8@coarse)).tolist(),
                            'refined_minus_coarse8': (factor*(fine-w8@coarse)).tolist(), 'diagnostics': diag,
                            'tables': {k: v.tolist() for k, v in tables.items()}}
        arrays.update({prefix+'/energies': E, prefix+'/at_zero': np.array([r['at_zero'] for r in rows]),
                       prefix+'/at_one': np.array([r['at_one'] for r in rows]), prefix+'/canonical': F,
                       prefix+'/bilinears': B, prefix+'/complex_frequency': np.array(zcheck),
                       prefix+'/complex_bilinears': Bc, prefix+'/reflection_frequency': np.array(frequencies),
                       prefix+'/reflection': np.array([reflection[z] for z in frequencies]),
                       prefix+'/coarse_energy': E8, prefix+'/coarse_weight': w8,
                       prefix+'/coarse_bilinears': oldB, prefix+'/coarse_kernels': coarse})
        for key in ('inner_basis_solves', 'short_radial_solves'):
            cumulative[key] = cumulative.get(key, 0)+sum(int(row[key]) for row in [*rows, check])
        cumulative['reflection_nodes'] = cumulative.get('reflection_nodes', 0)+len(rr)
        cumulative['reflection_seconds'] = cumulative.get('reflection_seconds', 0.)+sum(float(r['seconds']) for r in rr)
        cumulative['basis_seconds'] = cumulative.get('basis_seconds', 0.)+sum(float(r['seconds']) for r in [*rows, check])
        print(f'group14 sign {sign}: {selected}, delta={per_sign[prefix]["refined_minus_coarse8"]}', flush=True)
    total = sum(np.array(v['physical_contribution']) for v in per_sign.values())
    coarse = sum(np.array(v['coarse8_physical_contribution']) for v in per_sign.values())
    errors = {}
    if stage != 'base':
        for name, a, b in (('nested9_to17', 'n17_q48_h3', 'n9_q48_h3'),
                           ('contour24_to48', 'n17_q48_h3', 'n17_q24_h3')):
            errors[name] = sum(abs(factor*(np.array(v['tables'][a])-v['tables'][b])) for v in per_sign.values()).tolist()
    if stage == 'height':
        errors['height_kappa3_to4'] = sum(abs(factor*(np.array(v['tables']['n17_q48_h3'])-v['tables']['n17_q48_h4'])) for v in per_sign.values()).tolist()
    report = {'stage': stage, 'group': 14, 'kernel_order': ['rho', 'p_parallel', 'T01', 'p_perp'],
              'group_factor_per_angular_sign': factor, 'per_sign': per_sign,
              'physical_group14_contribution': total.tolist(), 'coarse8_physical_contribution': coarse.tolist(),
              'refined_minus_coarse8': (total-coarse).tolist(), 'error_indicators_sum_absolute_signs': errors,
              'new_stage_cost': cache.cost, 'cumulative_selected_node_cost': cumulative,
              'stage_wall_seconds': perf_counter()-clock, 'cache_directory': str(cache.directory),
              'scope': 'bounded group14 +/- local subgap source; no other groups, state change, history or metric evolution'}
    if stage == 'height':
        metadata = {'schema': 'NSC-INCOMING-GROUP14-SUBGAP-PILOT-v1', 'report': report,
                    'node_signature': cache.signature, 'config': cfg,
                    'source_hashes': {p: _hash(root/p) for p in (*OWNER_PATHS,
                        'src/recursive_horizons/nsc_incoming_subgap_completion.py',
                        'src/recursive_horizons/nsc_incoming_subgap_source.py',
                        'src/recursive_horizons/nsc_incoming_state_moments.py',
                        'src/recursive_horizons/nsc_common_ks_trace.py')},
                    'input_payloads': [mode_record['payload'], {'path': COARSE_PANELS, 'sha256': _hash(root/COARSE_PANELS)}]}
        arrays['metadata_json'] = np.frombuffer(json.dumps(metadata, sort_keys=True).encode(), np.uint8)
        raw = deterministic_npz_bytes(arrays); digest = sha256(raw).hexdigest()
        path = root/f'results/development/artifacts/nsc-incoming-subgap-completion.{digest}.npz'
        if not path.exists(): path.write_bytes(raw)
        report['payload'] = {'path': str(path.relative_to(root)), 'sha256': digest, 'bytes': len(raw)}
    return report
