"""Group-13 incoming local moments from archived unsewn subgap columns.

Only stationary interior continuation from rho=0 to rho=1 is new. The
17 real basis nodes, conjugate-frequency control, normalization control and
all contour reflections are inherited bytes. Complex columns are analytic
bilinear inputs, never Gaussian covariance states. No horizon, scattering or
history propagation is performed, and no other subgap family is inferred.
"""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from .nsc_common_ks_trace import restrict_resolved_modes
from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain, I2, S1, S2, S3
from .nsc_lorentzian import geometry
from .nsc_subgap_history_response import AnalyticResponsePanel, source_functions
from .nsc_incoming_state_moments import (
    AXIAL, RADIUS, incoming_adiabatic_reference, incoming_group_factor, incoming_state_moments,
)
from .nsc_paired_horizon_preparation import source_covariance
from .nsc_mode_resolved_cauchy_state import deterministic_npz_bytes


MASS = np.pi/2
ORDER = ('rho', 'p_parallel', 'T01', 'p_perp')
SUBGAP_RECORD = 'results/development/nsc-subgap-history-response.json'
MODE_RECORD = 'results/development/nsc-pg-massive-mode-resolution.json'
STATE_RECORD = 'results/development/nsc-mode-resolved-cauchy-state.json'
COARSE_PANELS = 'results/development/artifacts/nsc-incoming-spectral-panels.efbce40b2f959c7d0294689661d4e2e7193ade965ece0918c72c5a072ec8ad73.npz'
OWNERS = (
    'src/recursive_horizons/nsc_incoming_subgap_source.py',
    'src/recursive_horizons/nsc_common_ks_trace.py',
    'src/recursive_horizons/nsc_common_time_bulk_split.py',
    'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
    'src/recursive_horizons/nsc_lorentzian.py',
    'src/recursive_horizons/nsc_subgap_history_response.py',
    'src/recursive_horizons/nsc_incoming_state_moments.py',
    'src/recursive_horizons/nsc_paired_horizon_preparation.py',
    'src/recursive_horizons/nsc_compact_ctp_neck.py',
    'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py',
    'src/recursive_horizons/nsc_ks_spacetime_variation.py',
)


def _read(root, path): return json.loads((Path(root)/path).read_text())
def _hash(root, path): return sha256((Path(root)/path).read_bytes()).hexdigest()
def _load(path):
    with np.load(path, allow_pickle=False) as value: return {k: value[k].copy() for k in value.files}


def _payload(root, record):
    item = record['payload']
    if _hash(root, item['path']) != item['sha256']:
        raise ValueError('archived physical input payload changed: '+item['path'])
    return _load(Path(root)/item['path'])


def continue_archived_basis(energy, at_zero, at_zero_dual):
    """One short stationary solve, with a jointly propagated complex dual.

    The generator is the owned homogeneous PG Dirac radial equation. At a
    real frequency the archived dual must equal the upper columns. For a
    complex frequency both separately archived initial bases are evolved in
    one block-diagonal solve. This is not a history or state propagator.
    """
    z = complex(energy); initial = np.asarray(at_zero, complex); dual = np.asarray(at_zero_dual, complex)
    if z.real <= 0 or initial.shape != (2, 2) or dual.shape != (2, 2) or not np.isfinite([initial, dual]).all():
        raise ValueError('positive-frequency archived unsewn columns and their independent dual required')
    real = z.imag == 0
    if real and np.max(abs(initial-dual)) > 3e-13:
        raise ValueError('real archived upper/dual bases disagree')
    frequencies = np.array([z] if real else [z, z.conjugate()])
    bases = np.array([initial] if real else [initial, dual])
    def rhs(rho, value):
        beta, bp, _ = geometry(rho)
        potential = CommonTimeBulkSplit.local_potential(1., np.sqrt(1+rho*rho), MASS, 0.)
        matrices = value.reshape(len(frequencies), 2, 2)
        return np.array([np.linalg.solve(S2-beta*I2, 1j*(E*I2-potential)+.5*bp*I2)@field
                         for E, field in zip(frequencies, matrices)]).ravel()
    run = solve_ivp(rhs, (0., 1.), bases.ravel(), method='DOP853', rtol=2e-13, atol=2e-15)
    if not run.success: raise ArithmeticError(run.message)
    end = run.y[:, -1].reshape(len(frequencies), 2, 2)
    upper, lower = end[0], end[0] if real else end[1]
    velocity = lambda rho: S2-float(geometry(rho)[0])*I2
    initial_current = dual.conj().T@velocity(0.)@initial
    final_current = lower.conj().T@velocity(1.)@upper
    return upper, lower, {'function_evaluations': int(run.nfev), 'short_stationary_solves': 1,
                         'cross_current_residual': float(np.max(abs(final_current-initial_current)))}


def canonical_analytic_basis(energy, pg_basis, clock_shift):
    """Exact same-surface trace inverse, retaining the analytic clock phase.

    This is deliberately separate from the real-frequency Gaussian state API.
    Its matrix is the one used by restrict_resolved_modes at the same rho=1.
    """
    phi = np.asarray(pg_basis, complex)
    if phi.shape != (2, 2) or not np.isfinite(phi).all():
        raise ValueError('two finite unsewn PG columns required')
    beta = float(geometry(1.)[0])
    domain = TransmittingDiracSeamDomain(1., 1., beta, RADIUS,
                                       surface_id='reference-KS-Cauchy-section-rho=1')
    _, inverse, _ = domain.trace_map(np.ones(1))
    return (inverse[0]/RADIUS)@phi*np.exp(1j*complex(energy)*clock_shift)


def local_vertices(energy):
    z = complex(energy); p = -z/AXIAL
    return np.array([-MASS*S1+p*S3, p*S3, -p*I2, np.zeros((2, 2), complex)])


def local_analytic_bilinears(energy, canonical_upper, canonical_dual):
    """B_ab(z)=Phi_a(bar z)^dagger V(z) Phi_b(z), ordered (v,w)."""
    upper, dual = np.asarray(canonical_upper), np.asarray(canonical_dual)
    if upper.shape != (2, 2) or dual.shape != (2, 2):
        raise ValueError('separate two-column analytic basis and dual required')
    return np.array([dual.conj().T@vertex@upper for vertex in local_vertices(energy)])


def real_reference_remainder(energies):
    """Half-identity source term minus the EXISTING Pad4 vertex, real axis only."""
    raw = np.asarray(energies)
    if np.iscomplexobj(raw): raise ValueError('adiabatic source subtraction stays on the real axis')
    reference, _ = incoming_adiabatic_reference(raw, MASS, 0.)
    vertices = np.array([local_vertices(E) for E in raw])
    return (.5*np.trace(vertices, axis1=-2, axis2=-1)
            -np.einsum('nij,nvji->nv', reference, vertices)).real


def integrate_local_subgap(panel, kappa, reflection, *, points=48, height_fraction=3):
    """Positive-E local kernels; no history-action -1/pi factor is applied.

    The reflected coherence alone is continued. The final physical group
    factor is incoming_group_factor, applied once outside this integral.
    ``reflection`` must be an exact archived-node lookup with no fallback.
    """
    if points not in (24, 48) or (points, height_fraction) not in ((24, 3), (48, 3), (48, 4)):
        raise ValueError('only the three archived exact contour quadratures are owned')
    x, w = leggauss(points); left, right = panel.left, panel.right
    E = left+(x+1)*(right-left)/2; weights = w*(right-left)/2
    smooth = np.zeros(4, complex)
    for e, weight in zip(E, weights):
        B = panel(e); d, _ = source_functions(e, kappa)
        smooth += weight*d*(B[:, 0, 0]-B[:, 1, 1])
    reference = weights@real_reference_remainder(E)
    height = kappa/height_fraction; reflected = np.zeros(4, complex)
    for a, b in ((complex(left), left+1j*height),
                 (left+1j*height, right+1j*height), (right+1j*height, complex(right))):
        for xi, wi in zip(x, w):
            z = a+(xi+1)*(b-a)/2
            _, s = source_functions(z, kappa)
            reflected += wi*(b-a)/2*(-1j*s*reflection(z)*panel(z)[:, 1, 0])
    return {'smooth': smooth, 'reference': reference, 'reflected': reflected,
            'integral': smooth.real+reference+2*reflected.real,
            'smooth_imaginary_residual': float(np.max(abs(smooth.imag)))}


def recompose_incoming_subgap(arrays):
    """Replay local contractions and all comparisons from the saved columns."""
    meta = json.loads(arrays['metadata_json'].tobytes()); cfg = meta['config']; kappa = cfg['surface_gravity']
    if meta['group'] != 13 or meta['mass'] != MASS or meta['angular'] != 0:
        raise ValueError('only the retained angular-zero compact group13 is owned')
    E = arrays['energies']; B = arrays['bilinears']
    lookup = dict(zip(arrays['reflection_frequency'], arrays['reflection'])); used = set()
    def reflection(z):
        z = complex(z)
        if z not in lookup: raise ValueError('requested reflection node is absent from authenticated archive')
        used.add(z); return lookup[z]
    tables = {}
    for stride, label in ((1, '17'), (2, '9')):
        panel = AnalyticResponsePanel(E[::stride], B[::stride], 1., MASS)
        for points, fraction in ((24, 3), (48, 3), (48, 4)):
            row = integrate_local_subgap(panel, kappa, reflection, points=points, height_fraction=fraction)
            tables[f'n{label}_q{points}_h{fraction}'] = row
    fine = tables['n17_q48_h3']['integral']; factor = meta['group_factor']
    coarse = arrays['coarse/weights']@arrays['coarse/archived_kernels']
    panel = AnalyticResponsePanel(E, B, 1., MASS)
    expected_complex = arrays['complex_control/bilinears']
    diagnostics = dict(meta['continuation_residuals'])
    diagnostics.update(
        complex_interpolation=float(np.max(abs(panel(arrays['complex_control/frequency'].item())-expected_complex))),
        real8_bilinear_interpolation=float(np.max(abs(np.array([panel(e) for e in arrays['coarse/energies']])-arrays['coarse/bilinears']))),
        real8_formula_vs_archived=float(np.max(abs(arrays['coarse/formula_kernels']-arrays['coarse/archived_kernels']))),
        real8_formula_vs_direct=float(np.max(abs(arrays['coarse/formula_kernels']-arrays['coarse/direct_kernels']))),
        T01_coherence_overlap=float(np.max(abs(B[:, 2, 1, 0]))),
        T01_diagonal_difference=float(np.max(abs(B[:, 2, 0, 0]-B[:, 2, 1, 1]))),
        smooth_imaginary=max(row['smooth_imaginary_residual'] for row in tables.values()),
    )
    error = {
        'nested9_to17': abs(fine-tables['n9_q48_h3']['integral'])*factor,
        'contour24_to48': abs(fine-tables['n17_q24_h3']['integral'])*factor,
        'height_kappa3_to4': abs(fine-tables['n17_q48_h4']['integral'])*factor,
    }
    return {'kernel_order': ORDER, 'raw_integral': fine, 'physical_group13_contribution': factor*fine,
            'coarse8_raw_integral': coarse, 'coarse8_physical_contribution': factor*coarse,
            'refined_minus_coarse8': factor*(fine-coarse), 'error_indicators': error,
            'diagnostics': diagnostics, 'tables': tables, 'reflection_samples_used': len(used),
            'short_stationary_solves': meta['short_stationary_solves'],
            'scope': 'group13 finite [1,pi/2] source panel only; no other group, tail, local restoration or full stress'}


def prepare_incoming_subgap(root):
    """Create one small content-addressed input artifact with at most 19 solves."""
    root = Path(root); source_record = _read(root, SUBGAP_RECORD); old = _payload(root, source_record)
    mode_record = _read(root, MODE_RECORD); mode = _payload(root, mode_record)
    if _hash(root, COARSE_PANELS) not in Path(COARSE_PANELS).name:
        raise ValueError('coarse incoming panel payload digest changed')
    coarse_panels = _load(root/COARSE_PANELS)
    state = _read(root, STATE_RECORD); channel = state['channels'][13]
    cfg = json.loads(old['metadata_json'].tobytes())['config']
    if (channel['index'] != 13 or channel['compact_mass'] != MASS or channel['angular_eigenvalue'] != 0
            or cfg['omega']*MASS <= MASS or len(old['energies']) != 17 or len(old['reflection']) != 360):
        raise ValueError('archived group13 basis/reflection or inherited closed-channel domain changed')
    signature = {p: _hash(root, p) for p in (*OWNERS, SUBGAP_RECORD, MODE_RECORD, STATE_RECORD)}
    at_one, dual_one, diagnostics = [], [], []
    for i, E in enumerate(old['energies']):
        upper, dual, row = continue_archived_basis(E, old[f'node{i}/at_zero'], old[f'node{i}/at_zero_dual'])
        at_one.append(upper); dual_one.append(dual); diagnostics.append(row)
    at_one, dual_one = np.array(at_one), np.array(dual_one)
    padded = np.concatenate((at_one, np.zeros((17, 2, 1))), axis=2)
    canonical, chart = restrict_resolved_modes(1., old['energies'], padded)
    canonical = canonical[..., :2]; S = chart['clock_shift']
    arrays = {'energies': old['energies'], 'at_one': at_one, 'at_one_dual': dual_one,
              'canonical': canonical, 'reflection_frequency': old['reflection_frequency'], 'reflection': old['reflection']}
    arrays['bilinears'] = np.array([local_analytic_bilinears(e, f, f) for e, f in zip(old['energies'], canonical)])
    frame_residual = max(float(np.max(abs(canonical_analytic_basis(e, phi, S)-f)))
                         for e, phi, f in zip(old['energies'], at_one, canonical))
    for name, key in (('complex_control', 'complex_control_frequency'), ('normalization_control', 'normalization_control_frequency')):
        z = old[key].item()
        upper, dual, row = continue_archived_basis(z, old[name+'/at_zero'], old[name+'/at_zero_dual'])
        diagnostics.append(row)
        F, Fd = canonical_analytic_basis(z, upper, S), canonical_analytic_basis(complex(z).conjugate(), dual, S)
        arrays.update({name+'/frequency': np.array(z), name+'/at_one': upper, name+'/at_one_dual': dual,
                       name+'/canonical': F, name+'/canonical_dual': Fd,
                       name+'/bilinears': local_analytic_bilinears(z, F, Fd),
                       name+'/dual_bilinears': local_analytic_bilinears(complex(z).conjugate(), Fd, F)})
    # The old eight REAL source rows are used only for normalization and the
    # coarse-panel comparison, never as a reflection interpolation table.
    labels = mode['label']; indices = np.flatnonzero((labels[:, 0] == 13)&(labels[:, 1] == 1)&(labels[:, 5] > 1)&(labels[:, 5] < MASS))
    if len(indices) != 8 or not np.array_equal(mode['rho_interior'], [-1., 0., 1.]):
        raise ValueError('archived real subgap controls do not match the same incoming slice')
    E8, w8 = labels[indices, 5], labels[indices, 6]
    phi8 = mode['interior'][indices, 2]; unsewn8 = mode['interior_fundamental'][indices, 2][:, :, [1, 0]]
    unsewn_padded = np.concatenate((unsewn8, np.zeros((8, 2, 1))), axis=2)
    F8, _ = restrict_resolved_modes(1., E8, unsewn_padded); F8 = F8[..., :2]
    B8 = np.array([local_analytic_bilinears(e, f, f) for e, f in zip(E8, F8)])
    R8 = mode['reflection'][indices, 0]
    formula8 = real_reference_remainder(E8)
    for i, (e, B, R) in enumerate(zip(E8, B8, R8)):
        d, s = source_functions(e, cfg['surface_gravity'])
        formula8[i] += (d*(B[:, 0, 0]-B[:, 1, 1])+2*np.real(-1j*s*R*B[:, 1, 0])).real
    direct8 = incoming_state_moments(E8, phi8, mode['source_covariance'][indices], mode['source_projector'][indices], MASS, 0.)
    if not np.array_equal(E8, coarse_panels['subgap8/13_1/energies']) or not np.array_equal(w8, coarse_panels['subgap8/13_1/weights']):
        raise ValueError('incoming coarse quadrature rows differ from their physical mode input')
    arrays.update({'coarse/energies': E8, 'coarse/weights': w8, 'coarse/bilinears': B8,
                   'coarse/formula_kernels': formula8, 'coarse/direct_kernels': direct8['kernels'],
                   'coarse/archived_kernels': coarse_panels['subgap8/13_1/kernels']})
    v0 = old['normalization_control/at_zero']; physical0 = old['control/physical_at_zero']
    Rcontrol = np.vdot(v0[:, 0], physical0[:, 0])/np.vdot(v0[:, 0], v0[:, 0])
    reproduced0 = np.column_stack((Rcontrol*v0[:, 0], v0[:, 1], np.zeros(2)))
    gram = canonical.swapaxes(-1, -2).conj()@canonical
    complex_gram = arrays['complex_control/canonical_dual'].conj().T@arrays['complex_control/canonical']
    metadata = {
        'schema': 'NSC-INCOMING-GROUP13-SUBGAP-v1', 'source_hashes': signature,
        'input_payloads': [source_record['payload'], mode_record['payload'],
                           {'path': COARSE_PANELS, 'sha256': _hash(root, COARSE_PANELS)}],
        'group': 13, 'mass': MASS, 'angular': 0., 'config': cfg,
        'group_factor': incoming_group_factor(channel, (1,)), 'clock_shift': S,
        'short_stationary_solves': len(diagnostics), 'function_evaluations': sum(d['function_evaluations'] for d in diagnostics),
        'continuation_residuals': {
            'cross_current': max(d['cross_current_residual'] for d in diagnostics),
            'canonical_frame_vs_real_owner': frame_residual,
            'real_basis_unitarity': float(np.max(abs(gram-I2))),
            'complex_cross_Gram': float(np.max(abs(complex_gram-I2))),
            'complex_bilinear_adjoint': float(np.max(abs(arrays['complex_control/bilinears']-arrays['complex_control/dual_bilinears'].swapaxes(-1, -2).conj()))),
            'normalization_control_at_zero': float(np.max(abs(reproduced0-physical0))),
            'normalization_control_reflection_modulus': float(abs(abs(Rcontrol)-1)),
        },
        'normalization': 'positive-E local kernels; incoming_group_factor counts signed counterpart and copies once',
        'no_history_action_minus_one_over_pi': True, 'horizon_scattering_history_solves': 0,
        'complex_frequency_is_Gaussian_state': False,
    }
    if metadata['short_stationary_solves'] != 19: raise ArithmeticError('stationary continuation budget changed')
    if signature != {p: _hash(root, p) for p in signature}: raise ValueError('source changed during bounded continuation')
    arrays['metadata_json'] = np.frombuffer(json.dumps(metadata, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays); digest = sha256(raw).hexdigest()
    path = root/f'results/development/artifacts/nsc-incoming-subgap-source.{digest}.npz'
    if path.exists() and path.read_bytes() != raw: raise ValueError('content-addressed artifact collision')
    if not path.exists(): path.write_bytes(raw)
    return path, recompose_incoming_subgap(arrays)


def read_prepared_incoming_subgap(root, path):
    path = Path(path); raw = path.read_bytes()
    if sha256(raw).hexdigest() not in path.name: raise ValueError('incoming subgap artifact hash changed')
    arrays = _load(path); meta = json.loads(arrays['metadata_json'].tobytes())
    for p, digest in meta['source_hashes'].items():
        if _hash(root, p) != digest: raise ValueError('incoming subgap dependency changed: '+p)
    for p in meta['input_payloads']:
        if _hash(root, p['path']) != p['sha256']: raise ValueError('incoming subgap input bytes changed')
    return arrays
