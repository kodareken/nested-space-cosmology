"""Read the original shared parent's mode and interface content, without GR runs.

Item1 uses only producer-bound epsilon1 two-mode GH histories. No new
preparation, initial constraint solve, filtered waveform or evolution is made.
The parent-history coordinate time and reference-observer proper-time labels
are distinct. Larger incoming characteristic derivatives of a massive
standing field do not by themselves establish a one-way travelling pulse.
"""
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
import math
import multiprocessing
from pathlib import Path
import subprocess
import time

import numpy as np
from numpy.polynomial import Polynomial
from numpy.polynomial.legendre import leggauss
from scipy.fft import set_workers
from scipy.optimize import brentq, minimize_scalar
from threadpoolctl import threadpool_limits

from recursive_horizons.nsc_atlas_gh_history import DenseGHHistory
from recursive_horizons.nsc_atlas_parent_gh import ParentGH, PAIRS
from recursive_horizons.nsc_atlas_parent_circulation import VOLUME
from recursive_horizons.nsc_atlas_s3_frame import arbitrary_points_values
from recursive_horizons.nsc_atlas_s3_hopf import build_hopf_frame_basis


ROOT = Path(__file__).resolve().parents[1]
CLOCK = {
    'coordinate_t': 'Actual lapse-controlled GH coordinate time; DenseGHHistory.evaluate takes this argument.',
    'reference_observer_proper_duration': 'Elapsed proper duration of one tracked geodesic, initially normal at the prepared scalar-maximum anchor. It labels the GH leaf; it is not global proper time or the interface observers\' clock.',
    'epoch': 't=tau_reference=0 is the declared prepared Cauchy slice, not a derived cosmological origin.',
    'surface': 'Each interface readout lies on the GH coordinate-time leaf selected by the reference geodesic label. Different points have different local tube lapse.',
}
TUBE_LABELS = (0., .02822386, .0538334449, .075, .3, .6, .9)
MODE_KEYS = ('reference_instantaneous', 'reference_quadrature', 'reference_local_mass_quadrature',
             'physical_instantaneous', 'physical_quadrature', 'physical_joint_two_template_quadrature')


def file_sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_neck_reference():
    """Reuse the recorded critical family; never call its constraint solver."""
    path = ROOT/'results/summaries/nsc-atlas-valve-neck.json'
    record = json.loads(path.read_bytes())
    if record.get('schema') != 'nsc-atlas-valve-neck-v1':
        raise ValueError('Recorded one-mode neck reference required')
    row = next(r for r in record['critical_amplitudes'] if r['degree'] == 32)
    return dict(path=path.relative_to(ROOT).as_posix(), sha256=file_sha256(path),
                epsilon_critical=row['epsilon_crit'], psi_critical_at_plane=row['psi_neck'],
                degree=32, full_conformal_residual=row['full_conformal_PDE_RMS'],
                scope='Numerical time-symmetric SO3 one-mode reference; not a certified or state-independent amplitude threshold.')


def verify_parent_record(record):
    """Check original source identity before constructing any decoder."""
    if (record.get('schema') not in ('nsc-atlas-parent-gh-v2', 'nsc-atlas-parent-gh-v3')
            or record.get('epsilon') != 1
            or not record.get('source_bound_to_commit')
            or not record.get('producer_source_unchanged_during_run')
            or (record.get('initial_kick') or {}).get('kick_amplitude', 0.) != 0
            or record.get('spatial_transform') != 'hopf'):
        raise ValueError('Original source-bound, unchanged two-mode/no-kick history required')
    for name, expected in record['source_sha256'].items():
        path = Path(name)
        if path.is_absolute() or '..' in path.parts:
            raise ValueError('Invalid archived source path')
        archived = subprocess.run(['git', 'show', record['producer_commit']+':'+path.as_posix()],
                                  cwd=ROOT, capture_output=True)
        if archived.returncode or hashlib.sha256(archived.stdout).hexdigest() != expected:
            raise ValueError('Original producer Git source mismatch: '+name)
        if file_sha256(ROOT/path) != expected:
            raise ValueError('Current original-parent decoder/source changed: '+name)


def load_original_parent(degree):
    """Read producer-bound dense data; constructor does not initialize/evolve."""
    base = ROOT/f'.build/valve-e1-L{degree}-extended-v1'
    record_path = base.with_suffix('.json')
    record = json.loads(record_path.read_bytes())
    verify_parent_record(record)
    sparse_path = base.with_suffix('.npz')
    dense_path = Path(str(base)+'-history.npz')
    published_path = ROOT/'results/summaries/nsc-atlas-valve-formation-extended-v1.json'
    published = json.loads(published_path.read_bytes())
    published_row = next((r for r in published['evolution_runs'] if r['degree'] == degree), None)
    if published_row is not None:
        if (file_sha256(record_path) != published_row['record_sha256']
                or file_sha256(sparse_path) != published_row['state_cache_sha256']):
            raise ValueError('Published original run/sparse-cache pin failed')
    if file_sha256(dense_path) != record['dense_history']['sha256']:
        raise ValueError('Producer-pinned dense history failed')
    history = DenseGHHistory.load(dense_path)
    if (history.metadata['source_sha256'] != record['source_sha256']
            or history.metadata['producer_commit'] != record['producer_commit']):
        raise ValueError('Dense history source identity differs')
    system = ParentGH(degree, 1., gauge=record['prescribed_GH_gauge'], backend='hopf',
                      quadrature_order=history.metadata.get('quadrature_order'))
    # No system.initial(), system.rhs() or source-preparation function call.
    model = record['model_inputs']
    for name, value in (('H_reference', system.r), ('vacuum', system.parent.vacuum),
                        ('scalar_mass_squared', system.parent.mass_squared),
                        ('volume_kinetic', system.K), ('sstar', system.sstar)):
        if not math.isclose(model[name], value, rel_tol=2e-13, abs_tol=2e-14):
            raise ValueError('Decoder constants differ from original run: '+name)
    with np.load(sparse_path, allow_pickle=False) as cache:
        times, states = cache['times'], cache['states']
        metadata = json.loads(str(cache['metadata']))
        if metadata['source_sha256'] != record['source_sha256']:
            raise ValueError('Sparse cache source identity differs')
        if states.shape != (len(times), system.dof+10) or len(times) != len(record['snapshots']):
            raise ValueError('Original sparse cache shape differs')
        gaps = []
        physical_gaps, observer_gaps, clock_gaps, phase_rows = [], [], [], []
        dense_phase_initial = float(history.old[0][-1])
        for t, state, row in zip(times, states, record['snapshots']):
            if (abs(t-row['coordinate_t']) > 2e-12
                    or abs(state[system.dof+8]-row['observer_proper_duration']) > 2e-10):
                raise ValueError('Original event alignment failed')
            dense_state = history.evaluate(t)
            delta = abs(dense_state-state)
            # The original producer leaves the final ODE component at its
            # initial value (observer_rhs returns0), then writes the sampled
            # unwrapped PHASE into a sparse-copy state before state_sink.
            # It is neither a field nor observer proper time, and work is a
            # separate passive quadrature. Keep every physical/clock guard.
            gap = float(np.max(delta[:-1]))
            if gap > 2e-10*(1+float(np.max(abs(state[:-1])))):
                raise ValueError('Actual dense/sparse physical or observer states disagree')
            physical_gaps.append(float(np.max(delta[:system.dof])))
            observer_gaps.append(float(np.max(delta[system.dof:system.dof+8])))
            clock_gaps.append(float(delta[system.dof+8]))
            if clock_gaps[-1] > 2e-10:
                raise ValueError('Actual dense/sparse reference proper clocks disagree')
            if abs(float(dense_state[-1])-dense_phase_initial) > 2e-12:
                raise ValueError('Dense auxiliary phase is not the producer dummy component')
            reported_phase = row.get('observer_sampled_unwrapped_phase')
            phase_report_gap = (abs(float(state[-1])-reported_phase)
                                if reported_phase is not None else None)
            if phase_report_gap is not None and phase_report_gap > 2e-12:
                raise ValueError('Sparse sampled phase differs from its producer diagnostic')
            phase_rows.append(dict(coordinate_t=float(t),
                reference_observer_proper_duration=float(state[system.dof+8]),
                dense_dummy=float(dense_state[-1]), sparse_sampled_phase=float(state[-1]),
                dense_sparse_auxiliary_gap=float(delta[-1]), sparse_phase_report_gap=phase_report_gap))
            gaps.append(gap)
    endpoint = record['snapshots'][3]
    if not math.isclose(endpoint['observer_proper_duration'], .9, abs_tol=1e-12):
        raise ValueError('Published resolved checkpoint differs')
    return dict(system=system, history=history, record=record, end=endpoint['coordinate_t'],
                provenance=dict(record=record_path.relative_to(ROOT).as_posix(),
                    record_sha256=file_sha256(record_path), sparse_cache=sparse_path.relative_to(ROOT).as_posix(),
                    sparse_sha256=file_sha256(sparse_path), dense_cache=dense_path.relative_to(ROOT).as_posix(),
                    dense_sha256=file_sha256(dense_path), producer_commit=record['producer_commit'],
                    published_record_pin=published_row is not None,
                    published_summary_sha256=file_sha256(published_path), dense_sparse_maximum_gap=max(gaps),
                    dense_sparse_compared_components='All physical coefficients and reference-observer point, four-velocity and proper clock: state[:-1]',
                    dense_sparse_physical_maximum_gap=max(physical_gaps),
                    dense_sparse_observer_point_velocity_maximum_gap=max(observer_gaps),
                    dense_sparse_reference_clock_maximum_gap=max(clock_gaps),
                    auxiliary_phase=dict(index=system.dof+9,
                        role='Producer-owned sampled unwrapped observer phase; NOT physical field, proper clock or work ledger',
                        producer='nsc_atlas_parent_gh.observer_rhs returns0; evolve_parent assigns state[-1]=before_phase+partial only to the sparse-copy state before state_sink',
                        used_by_readout=False, values=phase_rows),
                    source_sha256=record['source_sha256']),
                reference_observer_initial_point=history.evaluate(0.)[system.dof:system.dof+4].tolist(),
                reference_observer_initial_velocity=history.evaluate(0.)[system.dof+4:system.dof+8].tolist())


def neck_normalization(parent, recorded):
    z = parent.mass_squared/parent.H**2
    A = math.sqrt(24/(parent.kappa*(z+3)))
    derivative = recorded['epsilon_critical']*A*parent.H/recorded['psi_critical_at_plane']**2
    return dict(Aneck=A, mass_squared_over_H_squared=z, H_reference=parent.H,
                mass_squared_reference=parent.mass_squared, kappa=parent.kappa,
                epsilon_critical=recorded['epsilon_critical'],
                psi_critical_at_plane=recorded['psi_critical_at_plane'],
                critical_normal_derivative=derivative, critical_abs_win=derivative/2,
                critical_abs_wout=derivative/2, critical_net_flux=0.,
                original_analytical_initial_reference_X1_amplitude=1/math.sqrt(2),
                scope='Critical family has PiS=PiChi=Kij=0,s=sstar,chi=epsilon*Aneck*X1 and its own solved SO3 conformal metric. Its standing in/out traces are equal; content comparison is not parent admission.')


class ParentModeProjection:
    """Reference and physical-measure projections on the original GH leaf."""
    def __init__(self, data, normalization):
        self.system, self.history = data['system'], data['history']
        self.basis = self.system.basis
        self.ell = np.array([label[0] for label in self.basis.labels])
        self.cX = self.basis.project(self.basis.points[:, 0])
        if np.linalg.norm(self.cX[self.ell != 1]) > 1e-11:
            raise ValueError('X1 representation is not the expected pure l1 mode')
        self.cX[self.ell != 1] = 0.
        self.DX = float(self.cX@self.cX)
        self.A = normalization['Aneck']
        self.omega = self.system.parent.reference_frequency
        self.cache = {}

    def evaluate(self, t, basis=None):
        t = float(t)
        if basis is None and t in self.cache:
            return self.cache[t]
        s = self.system
        state = self.history.evaluate(t)
        C = s.unpack(state)
        rule = self.basis if basis is None else basis
        gradient = np.array([s.r*(D@C[21]) for D in s.basis.generators])
        values = rule.evaluate(np.vstack([C[:10], C[20], C[21], C[23], gradient]).T)
        g = values[:, :10]
        h = np.empty((len(values), 3, 3))
        h[:, 0, 0] = g[:, 4]; h[:, 0, 1] = h[:, 1, 0] = g[:, 5]
        h[:, 0, 2] = h[:, 2, 0] = g[:, 6]; h[:, 1, 1] = g[:, 7]
        h[:, 1, 2] = h[:, 2, 1] = g[:, 8]; h[:, 2, 2] = g[:, 9]
        determinant = np.linalg.det(h)
        if (min(determinant) <= 0 or min(h[:, 0, 0]) <= 0
                or min(h[:, 0, 0]*h[:, 1, 1]-h[:, 0, 1]**2) <= 0):
            raise ValueError('Projection metric is not positive definite')
        inverse = np.linalg.inv(h)
        shift = np.einsum('nij,nj->ni', inverse, g[:, 1:4])
        lapse2 = -g[:, 0]+np.einsum('ni,ni->n', g[:, 1:4], shift)
        if min(lapse2) <= 0:
            raise ValueError('Projection slice has nonpositive lapse')
        S, chi, p_coordinate = values[:, 10:13].T
        Pi = (p_coordinate-np.einsum('ni,ni->n', shift, values[:, 13:16]))/np.sqrt(lapse2)
        X1, Q = rule.points[:, 0], rule.points[:, 1]*rule.points[:, 3]
        w = rule.weights
        physical = w*np.sqrt(determinant)/s.r**3
        denominator = float(physical@(X1*X1))
        a_ref = float(C[21]@self.cX/self.DX)
        p_ref = float(w@(Pi*X1)/(VOLUME/4))
        a_physical = float(physical@(chi*X1)/denominator)
        p_physical = float(physical@(Pi*X1)/denominator)
        local_s = s.sstar+S/math.sqrt(s.K)
        mu2 = 400*np.exp(-12*local_s)
        dX = s.r*np.column_stack([rule.points[:, 3], -rule.points[:, 2], rule.points[:, 1]])
        norm_dX = np.einsum('ni,nij,nj->n', dX, inverse, dX)
        omega_physical = math.sqrt(float(physical@(norm_dX+mu2*X1*X1)/denominator))
        gram = np.array([[denominator, physical@(X1*Q)], [physical@(X1*Q), physical@(Q*Q)]])
        fit = np.linalg.solve(gram, [physical@(chi*X1), physical@(chi*Q)])
        fit_pi = np.linalg.solve(gram, [physical@(Pi*X1), physical@(Pi*Q)])
        chi_norm, pi_norm = float(physical@(chi*chi)), float(physical@(Pi*Pi))
        power = float(C[21]@C[21])
        row = dict(coordinate_t=t, reference_observer_proper_duration=float(state[s.dof+8]),
            reference_instantaneous_signed=a_ref/self.A, reference_instantaneous=abs(a_ref)/self.A,
            reference_quadrature=math.hypot(a_ref, p_ref/self.omega)/self.A,
            reference_local_mass_quadrature=math.hypot(a_ref, float(w@(Pi*X1/np.sqrt(mu2))/(VOLUME/4)))/self.A,
            physical_instantaneous_signed=a_physical/self.A, physical_instantaneous=abs(a_physical)/self.A,
            physical_quadrature=math.hypot(a_physical, p_physical/omega_physical)/self.A,
            physical_joint_two_template_X1=float(fit[0]/self.A),
            physical_joint_two_template_quadrature=math.hypot(fit[0], fit_pi[0]/omega_physical)/self.A,
            physical_joint_X2X4_coefficient=float(fit[1]),
            physical_X1_only_relative_L2_residual=math.sqrt(float(physical@((chi-a_physical*X1)**2))/chi_norm),
            physical_two_template_relative_L2_residual=math.sqrt(float(physical@((chi-fit[0]*X1-fit[1]*Q)**2))/chi_norm),
            physical_Pi_two_template_relative_L2_residual=(math.sqrt(float(physical@((Pi-fit_pi[0]*X1-fit_pi[1]*Q)**2))/pi_norm) if pi_norm > 0 else 0.),
            physical_Gram_correlation=float(gram[0, 1]/math.sqrt(gram[0, 0]*gram[1, 1])),
            reference_even_l_power_fraction=float(sum(C[21, self.ell % 2 == 0]**2)/power),
            reference_l1_power_fraction=float(sum(C[21, self.ell == 1]**2)/power),
            physical_frozen_X1_frequency=omega_physical, reference_frequency=self.omega,
            normal_Pi_projection_reference=p_ref, normal_Pi_projection_physical=p_physical,
            physical_volume=float(sum(physical)),
            reference_projection_quadrature_error=float(w@(chi*X1)/(VOLUME/4)-a_ref))
        if basis is None:
            self.cache[t] = row
        return row

    def reference_crests(self, end):
        """Extrema of the actual dense l1 coefficient, not sparse interpolation."""
        s = self.system
        functional = np.zeros(s.dof+10)
        functional[s.offsets[21]:s.offsets[22]] = self.cX[s.indices[21]]/(self.DX*self.A)
        x, one_minus = Polynomial([0., 1.]), Polynomial([1., -1.])
        candidates = [0., end]
        for (start, stop), old, F in zip(self.history.bounds, self.history.old, self.history.coefficients):
            if start > end:
                break
            value = Polynomial([0.])
            for index, coefficient in enumerate((F@functional)[::-1]):
                value += coefficient
                value *= x if index % 2 == 0 else one_minus
            value += old@functional
            for root in value.deriv().roots():
                if abs(root.imag) < 1e-9 and 0 < root.real < 1:
                    t = float(start+(stop-start)*root.real)
                    if t < end:
                        candidates.append(t)
        return np.array(sorted(set(candidates)))

    def physical_instantaneous(self, t):
        state = self.history.evaluate(float(t)); C = self.system.unpack(state)
        values = self.basis.evaluate(C[[4, 5, 6, 7, 8, 9, 21]].T)
        h = np.empty((len(values), 3, 3))
        h[:, 0, 0] = values[:, 0]; h[:, 0, 1] = h[:, 1, 0] = values[:, 1]
        h[:, 0, 2] = h[:, 2, 0] = values[:, 2]; h[:, 1, 1] = values[:, 3]
        h[:, 1, 2] = h[:, 2, 1] = values[:, 4]; h[:, 2, 2] = values[:, 5]
        weight = self.basis.weights*np.sqrt(np.linalg.det(h)); X1 = self.basis.points[:, 0]
        return abs(float(weight@(values[:, 6]*X1)/(weight@(X1*X1)))/self.A)


def interface_rule(order=32):
    """Symmetry-preserving S2 quadrature on F=X1=0, no selected support."""
    z, w = leggauss(order); phi = 2*math.pi*np.arange(2*order)/(2*order)
    zz, pp = np.meshgrid(z, phi, indexing='ij'); ss = np.sqrt(1-zz*zz)
    X = np.stack([np.zeros_like(zz), ss*np.cos(pp), ss*np.sin(pp), zz], axis=-1).reshape(-1, 4)
    theta = np.stack([np.zeros_like(zz), zz*np.cos(pp), zz*np.sin(pp), -ss], axis=-1).reshape(-1, 4)
    azimuth = np.stack([np.zeros_like(zz), -np.sin(pp), np.cos(pp), np.zeros_like(zz)], axis=-1).reshape(-1, 4)
    weight = np.repeat(w, 2*order)*2*math.pi/(2*order)
    x, y, z, w = X.T
    frame = np.stack([np.c_[w, z, -y, -x], np.c_[-z, w, x, -y], np.c_[y, -x, w, -z]], axis=1)
    return X, theta, azimuth, weight, frame


def tube_frame(metric, points, frame_scale):
    """nB=g^-1dF/|dF|; UB is Eulerian normal projected into the tube."""
    h = metric[:, 1:, 1:]
    if (min(h[:, 0, 0]) <= 0 or min(h[:, 0, 0]*h[:, 1, 1]-h[:, 0, 1]**2) <= 0
            or min(np.linalg.det(h)) <= 0):
        raise ValueError('Interface slice metric is not positive definite')
    inverse = np.linalg.inv(metric)
    dF = np.zeros((len(points), 4))
    dF[:, 1:] = frame_scale*np.c_[points[:, 3], -points[:, 2], points[:, 1]]
    norm = np.einsum('ni,nij,nj->n', dF, inverse, dF)
    if min(norm) <= 0 or max(inverse[:, 0, 0]) >= 0:
        raise ValueError('Requested interface is not timelike everywhere; no node filtering allowed')
    nB = np.einsum('nij,nj->ni', inverse, dF)/np.sqrt(norm)[:, None]
    lapse = 1/np.sqrt(-inverse[:, 0, 0])
    shift = np.einsum('nij,nj->ni', np.linalg.inv(metric[:, 1:, 1:]), metric[:, 0, 1:])
    n = np.column_stack([1/lapse, -shift/lapse[:, None]])
    ndot = np.einsum('ni,ni->n', n, dF)/np.sqrt(norm)
    U = (n-ndot[:, None]*nB)/np.sqrt(1+ndot*ndot)[:, None]
    if min(U[:, 0]) <= 0:
        raise ValueError('Tube time orientation is not future directed')
    checks = dict(minimum_spacelike_normal_squared=float(min(norm)),
        U_norm_error=float(max(abs(np.einsum('ni,nij,nj->n', U, metric, U)+1))),
        nB_norm_error=float(max(abs(np.einsum('ni,nij,nj->n', nB, metric, nB)-1))),
        orthogonality_error=float(max(abs(np.einsum('ni,nij,nj->n', U, metric, nB)))))
    if max(checks[key] for key in ('U_norm_error', 'nB_norm_error', 'orthogonality_error')) > 1e-9:
        raise ValueError('Interface unit-frame controls failed')
    return U, nB, 1/U[:, 0], checks


def characteristic_content(U, nB, derivatives):
    """Return BOTH traces and stress flux, without an incoming-only filter."""
    tangent = np.einsum('ni,nia->na', U, derivatives)
    normal = np.einsum('ni,nia->na', nB, derivatives)
    incoming, outgoing = (tangent-normal)/2, (tangent+normal)/2
    return incoming, outgoing, -tangent*normal, tangent, normal


def _field_content(win, wout, j, area_weight, tube_lapse):
    area = float(sum(area_weight)); positive = j > 0
    return dict(area_RMS_win=math.sqrt(float(area_weight@(win*win)/area)),
        area_RMS_wout=math.sqrt(float(area_weight@(wout*wout)/area)),
        peak_abs_win=float(max(abs(win))), peak_abs_wout=float(max(abs(wout))),
        positive_j_area_fraction=float(sum(area_weight[positive])/area),
        area_integral_j=float(area_weight@j),
        area_integral_positive_j=float(area_weight@np.maximum(j, 0)),
        area_integral_negative_j=float(area_weight@np.minimum(j, 0)),
        P_GH=float(area_weight@(tube_lapse*j)),
        P_GH_positive=float(area_weight@(tube_lapse*np.maximum(j, 0))),
        P_GH_negative=float(area_weight@(tube_lapse*np.minimum(j, 0))),
        stress_trace_identity_error=float(max(abs(win*win-wout*wout-j))))


def interface_readouts(data, normalization):
    s, history, end = data['system'], data['history'], data['end']
    X, theta, azimuth, reference_weight, frame = interface_rule()
    Y = arbitrary_points_values(X, s.basis.labels)
    vt = np.einsum('nij,nj->ni', frame, theta)/s.r
    vp = np.einsum('nij,nj->ni', frame, azimuth)/s.r
    raw = dict(points=X, reference_S2_weights=reference_weight)
    snapshots = []
    raw_values = {name: [] for name in ('coordinate_t', 'reference_observer_proper_duration', 'area_weights',
        'tube_U', 'tube_nB', 'tube_lapse', 'canonical_fields', 'U_field_derivatives',
        'nB_field_derivatives', 'w_in', 'w_out', 'j_by_field', 'source_state_sha256')}
    for requested in TUBE_LABELS:
        if requested == 0.:
            t = 0.
        elif requested == .9:
            t = end
        else:
            t = brentq(lambda t: history.evaluate(t)[s.dof+8]-requested, 0., end, xtol=2e-12)
        state = history.evaluate(t); C = s.unpack(state)
        gradient = np.vstack([s.r*(D@C[[20, 21]].T).T for D in s.basis.generators])
        values = Y@np.vstack([C[:10], C[20:22], C[22:24], gradient]).T
        metric = np.empty((len(X), 4, 4))
        for k, (a, b) in enumerate(PAIRS):
            metric[:, a, b] = metric[:, b, a] = values[:, k]
        U, nB, tube_lapse, geometry = tube_frame(metric, X, s.r)
        h = metric[:, 1:, 1:]
        aa = np.einsum('ni,nij,nj->n', vt, h, vt)
        bb = np.einsum('ni,nij,nj->n', vp, h, vp)
        ab = np.einsum('ni,nij,nj->n', vt, h, vp)
        dA = reference_weight*np.sqrt(aa*bb-ab*ab)
        derivatives = np.stack([values[:, 12:14], values[:, 14:16], values[:, 16:18], values[:, 18:20]], axis=1)
        win, wout, j, tangent, normal = characteristic_content(U, nB, derivatives)
        fields = {name: _field_content(win[:, i], wout[:, i], j[:, i], dA, tube_lapse)
                  for i, name in enumerate(('radion', 'carrier'))}
        carrier = fields['carrier']; critical = normalization['critical_abs_win']
        carrier.update(RMS_win_over_critical=carrier['area_RMS_win']/critical,
                       peak_win_over_critical=carrier['peak_abs_win']/critical,
                       field_range_on_interface=[float(min(values[:, 11])), float(max(values[:, 11]))])
        total = j.sum(axis=1); support = total > 0
        positive_area = float(sum(dA[support]))
        snapshot = dict(coordinate_t=float(t), reference_observer_proper_duration=float(state[s.dof+8]),
            requested_reference_observer_label=requested, area=float(sum(dA)), geometry=geometry, fields=fields,
            total_area_integral_j=float(dA@total), total_area_integral_positive_j=float(dA@np.maximum(total, 0)),
            total_area_integral_negative_j=float(dA@np.minimum(total, 0)),
            total_P_GH=float(dA@(tube_lapse*total)),
            total_P_GH_positive=float(dA@(tube_lapse*np.maximum(total, 0))),
            total_P_GH_negative=float(dA@(tube_lapse*np.minimum(total, 0))),
            total_positive_j_area_fraction=positive_area/float(sum(dA)),
            carrier_RMS_win_on_total_positive_j=(math.sqrt(float(dA[support]@(win[support, 1]**2))/positive_area) if positive_area else None),
            carrier_RMS_wout_on_total_positive_j=(math.sqrt(float(dA[support]@(wout[support, 1]**2))/positive_area) if positive_area else None),
            carrier_peak_win_on_total_positive_j=float(max(abs(win[support, 1]))) if positive_area else None,
            carrier_peak_wout_on_total_positive_j=float(max(abs(wout[support, 1]))) if positive_area else None)
        snapshots.append(snapshot)
        payload = dict(coordinate_t=t, reference_observer_proper_duration=float(state[s.dof+8]),
            area_weights=dA, tube_U=U, tube_nB=nB, tube_lapse=tube_lapse,
            canonical_fields=values[:, 10:12], U_field_derivatives=tangent, nB_field_derivatives=normal,
            w_in=win, w_out=wout, j_by_field=j,
            source_state_sha256=hashlib.sha256(state.tobytes()).hexdigest())
        for name, value in payload.items():
            raw_values[name].append(value)
    raw.update({name: np.asarray(values) for name, values in raw_values.items()})
    return snapshots, raw


def _analyse_parent_degree(arguments):
    degree, raw_directory, recorded = arguments
    started = time.perf_counter()
    with threadpool_limits(limits=1), set_workers(2):
        data = load_original_parent(degree)
        norm = neck_normalization(data['system'].parent, recorded)
        projection = ParentModeProjection(data, norm)
        end = data['end']; crests = projection.reference_crests(end)
        times = np.unique(np.r_[np.linspace(0., end, 81), np.asarray(data['history'].bounds)[:, 1], crests])
        times = times[(times >= 0) & (times <= end)]
        samples = [projection.evaluate(t) for t in times]
        maxima = {}
        for key in MODE_KEYS:
            if key == 'reference_instantaneous':
                candidates = [projection.evaluate(t) for t in crests]
            elif key == 'physical_instantaneous':
                candidates = [projection.evaluate(0.), projection.evaluate(end)]
                for i, t in enumerate(crests[1:-1], 1):
                    answer = minimize_scalar(lambda t: -projection.physical_instantaneous(t),
                        bounds=((crests[i-1]+t)/2, (t+crests[i+1])/2), method='bounded', options={'xatol': 2e-9})
                    candidates.append(projection.evaluate(answer.x))
            else:
                best = max(samples, key=lambda row: row[key]); index = int(np.argmin(abs(times-best['coordinate_t'])))
                candidates = [best]
                if 0 < index < len(times)-1:
                    answer = minimize_scalar(lambda t: -projection.evaluate(t)[key],
                        bounds=(times[index-1], times[index+1]), method='bounded', options={'xatol': 2e-9})
                    candidates.append(projection.evaluate(answer.x))
            best = max(candidates, key=lambda row: row[key])
            maxima[key] = dict(value=best[key], coordinate_t=best['coordinate_t'],
                reference_observer_proper_duration=best['reference_observer_proper_duration'],
                even_l_field_power_fraction=best['reference_even_l_power_fraction'],
                physical_two_template_relative_L2_residual=best['physical_two_template_relative_L2_residual'],
                physical_Pi_two_template_relative_L2_residual=best['physical_Pi_two_template_relative_L2_residual'])
        checks = []
        if degree in (20, 24):
            order = data['system'].basis.quadrature_shape[:2]
            refined = build_hopf_frame_basis(degree, (order[0]+4, order[1]+8))
            for t in sorted(set([0.]+[row['coordinate_t'] for row in maxima.values()])):
                first, second = projection.evaluate(t), projection.evaluate(t, refined)
                checks.append(dict(coordinate_t=t, maximum_absolute_projection_change=max(abs(first[k]-second[k]) for k in MODE_KEYS)))
        diagnostic = projection.evaluate(data['record']['snapshots'][-1]['coordinate_t'])
        columns = list(next(iter(projection.cache.values())))
        mode_path = Path(raw_directory)/f'modes-L{degree}.npz'
        series = [row for _, row in sorted(projection.cache.items())]
        np.savez_compressed(mode_path, columns=np.asarray(columns),
            values=np.array([[row[key] for key in columns] for row in series]),
            within_primary_window=np.array([row['coordinate_t'] <= end for row in series]),
            primary_coordinate_interval=np.array([0., end]))
        traces, trace_metadata = [], None
        if degree in (20, 24):
            traces, arrays = interface_readouts(data, norm)
            trace_path = Path(raw_directory)/f'traces-L{degree}.npz'
            np.savez_compressed(trace_path, **arrays)
            trace_metadata = dict(path=trace_path.relative_to(ROOT).as_posix(), sha256=file_sha256(trace_path),
                                  field_order=['canonical_radion_S', 'carrier_chi'], quadrature_order=[32, 64])
    return dict(degree=degree, primary_claim_resolution=degree in (20, 24), normalization=norm,
        interval=dict(reference_observer_proper_duration=[0., .9], coordinate_t=[0., end],
            scope='Published resolved window on20/24;12/16 are coarse comparison diagnostics, not upgraded admission.'),
        initial=projection.evaluate(0.), numerically_located_maxima=maxima, resolved_endpoint=projection.evaluate(end),
        reference_observer=dict(initial_point=data['reference_observer_initial_point'],
            initial_velocity=data['reference_observer_initial_velocity'],
            worldline='Tracked geodesic, initially normal at the scalar-maximum source anchor',
            epoch='Declared prepared Cauchy slice t=tau_reference=0',
            original_endpoint_protocol=data['record']['model_inputs'].get('reference_proper_endpoint_origin')),
        proper1p2_diagnostic_only=diagnostic, projection_quadrature_checks=checks,
        field_basis=dict(mode1='X1, l1', mode2='X2X4, l2 even',
            GH_sector='Exact proper paired-flips(-1,-1), retaining both parity sectors and all admitted harmonics',
            l1_labels=[list(label) for label in data['system'].basis.labels if label[0] == 1]),
        interface_readouts=traces, raw_mode_series=dict(path=mode_path.relative_to(ROOT).as_posix(), sha256=file_sha256(mode_path)),
        raw_interface_arrays=trace_metadata, provenance=data['provenance'],
        snapshot_constraint_indicators=[{key:row[key] for key in
            ('coordinate_t', 'observer_proper_duration', 'ADM_H_RMS_over_H_squared', 'ADM_M_RMS_over_H_squared', 'GH_C_RMS_over_H')}
            for row in data['record']['snapshots']],
        evaluation_count=len(projection.cache), elapsed_seconds=time.perf_counter()-started)


def parent_drive_calculation(raw_directory):
    """Reproduce item1 only. No private helper/readout or prepared state input."""
    started = time.perf_counter()
    raw_directory = (ROOT/Path(raw_directory)).resolve()
    if not raw_directory.is_relative_to(ROOT/'.build'):
        raise ValueError('Raw parent-drive arrays must remain under.build')
    folder = raw_directory/'parent-drive'
    folder.mkdir(parents=True, exist_ok=False)
    recorded = load_neck_reference()
    with ProcessPoolExecutor(max_workers=2, mp_context=multiprocessing.get_context('spawn')) as pool:
        runs = list(pool.map(_analyse_parent_degree, [(degree, str(folder), recorded) for degree in (12, 16, 20, 24)]))
    return dict(schema='nsc-original-parent-drive-item1-v1', clock_definition=CLOCK,
        critical_family_reference=recorded, runs=runs,
        projection_definitions=dict(
            reference='a=<chi,X1>_dOmega/<X1,X1>_dOmega; p uses ACTUAL normalPi=(coordinatep-beta·Dchi)/N.',
            physical='Same templates with dVh=sqrt(deth)dOmega/Href³. The joint2x2Gram control retains X2X4 and reports its fit residual.',
            normalization='Aneck=sqrt(24/[kappa(mu_ref²/Href²+3)]), fixedoriginalParent.retained constants.',
            quadrature='sqrt(a²+(p/omega)²)/Aneck: phase-rotation-invariant for that instantaneous oscillator convention, NOT a conserved coupled-parent amplitude.',
            reference_frequency='sqrt(mu_ref²+3Href²)',
            physical_frequency='Frozen Rayleigh quotientomega²=integral(|DX1|h²+mu(s)²X1²)dVh/integralX1²dVh; X1 is not asserted to be an actual evolving eigenmode.',
            maximum='Passive actualDOP853polynomial evaluation, numerical maximization. Instantaneous physical crests all checked. No certified continuum-in-time error bound.'),
        interface_definitions=dict(F='X1=0 fixed timelike tube; no assumption carrierchi is odd under X1alone',
            normal='nB=g_inverse dF/sqrt(g_inverse(dF,dF)), oriented toward increasingX1',
            U='Future Eulerian normal projected tangent to F0 and normalized; orthogonal to its S2 cuts.',
            traces='win=(Uphi-nBphi)/2,wout=(Uphi+nBphi)/2; j=win²-wout²=-Uphi*nBphi, separately for BOTH canonical fields.',
            support='Actual j_total>0 labels local inward support; no field/trace deleted, no incoming-only initial state created.',
            area_integral='Integral j_B dA on the GHleaf: area integral of local tube-proper-unit flux, NOT dE/dreferenceobserver_tau.',
            coordinate_measure='P_GH=integral N_B*j_B*dA,N_B=1/U_B^0; also reported with positive/negative parts. No global Killing energy or anchor-clock energy derivative is asserted.',
            exact_symmetry='Proper paired flip exchanges F±. Uchi is odd/nBchi even; US is even/nBS odd. Each j changes sign: global net cancels without pointwise flux vanishing.',
            standing_field='Incoming/outgoing derivative RMS is dominated by massive normal-momentum oscillation. Large/equal traces are not evidence of a one-way propagating pulse.'),
        threshold_caveat='The recorded epsilon threshold requires the time-symmetric SO3 one-mode family, zero allnormalmomenta/Kij,s=sstar and its own constrained metric. Neither a projected amplitude nor incoming derivative content admits the actual two-mode parent as that state.',
        no_evolution_or_new_state=True, no_filtered_initial_data=True,
        child_formation_admitted=False, state_selection_proved=False,
        elapsed_seconds=time.perf_counter()-started)


def curvature_invariants(A,B,C,D,F=0.):
    """Independent spherical orthonormal curvature split, signature -+++."""
    ricci=-2*A-4*B+4*C+2*D
    ricci2=(A+2*B)**2+(-A+2*C)**2+2*(-B+C+D)**2-8*F*F
    kretsch=4*(A*A+2*B*B+2*C*C+D*D-4*F*F)
    weyl2=4/3*(A-B+C-D)**2
    return dict(Ricci_scalar=ricci,Ricci_squared=ricci2,Kretschmann=kretsch,Weyl_squared=weyl2,
        decomposition_residual=kretsch-2*ricci2+ricci*ricci/3-weyl2)


def conformal_sections(R,U,sections,power):
    """Diagnostic Omega=R^-power at the same reflection cut, Lunit=1.

    This changes the metric representative, not its physical matter state.
    """
    if not (math.isfinite(R) and R>0 and math.isfinite(U) and math.isfinite(power)):
        raise ValueError('finite radius/rate/power and positive radius required')
    A,B,C,D=(float(sections[k]) for k in ('A','B','C','D'))
    if not math.isclose(D,(1+U*U)/(R*R),rel_tol=2e-10,abs_tol=2e-10):
        raise ValueError('reflection radius-gradient identity required')
    p=power;factor=R**(2*p)
    out=dict(A=factor*(A-p*(B-C+D-R**-2)),B=(1-p)*factor*B,
        C=(1-p)*factor*(C-p*U*U/R**2),D=factor*((1-p)**2*D+p*(2-p)/R**2))
    rho=2*out['C']+out['D'];pr=2*out['B']-out['D'];pt=out['A']+out['B']-out['C']
    return dict(power=p,Omega=R**-p,rescaled_areal_radius=R**(1-p),
        rescaled_normal_radius_rate=(1-p)*U,sectional_components=out,
        invariants=curvature_invariants(**out),
        effective_Einstein_stress=dict(rho=rho,pr=pr,pt=pt,w_iso=(pr+2*pt)/(3*rho) if rho else None,
            radial_NEC=rho+pr,angular_NEC=rho+pt),
        scope='Unphysical metric diagnostic, not a new canonical-field state or an attached child.')


def saved_neck_metric(path,expected):
    """Use actual saved temporal/spatial metric jets, never source repair."""
    path=(ROOT/path).resolve()
    if not path.is_relative_to(ROOT/'.build') or file_sha256(path)!=expected:
        raise ValueError('authenticated saved neck profile required')
    with np.load(path,allow_pickle=False) as raw:
        y,yt,yx,yxx,ytt=(raw[k] for k in ('y','yt','yx','yxx','ytt'))
        x=raw['x'];proper=float(raw['center_proper_time'])
    i=y.shape[1]//2;R,a=y[:2,i];RT,aT=yt[:2,i];RX,aX=yx[:2,i]
    RXX,aXX=yxx[:2,i];RTT,aTT=ytt[:2,i];em=math.exp(-a);em2=em*em
    Rtx=(yt[0,i+1]-yt[0,i-1])/(x[i+1]-x[i-1])
    A=em2*(aXX-aTT);B=-em2*(RTT-aT*RT-aX*RX)/R
    C=-em2*(RXX-aT*RT-aX*RX)/R;D=(1+em2*(RT*RT-RX*RX))/(R*R)
    F=-em2*(Rtx-aX*RT-aT*RX)/R
    Kr=-em*aT;Kt=-em*RT/R;Rll=em2*(RXX-aX*RX)
    sections=dict(A=float(A),B=float(B),C=float(C),D=float(D),F=float(F))
    invariants=curvature_invariants(**sections);U=em*RT
    p1=conformal_sections(float(R),float(U),sections,1.)
    p15=conformal_sections(float(R),float(U),sections,1.5)
    p15['local_volume_expansion']=R**1.5*(-Kr+(2-4.5)*U/R)
    return dict(raw_cache=path.relative_to(ROOT).as_posix(),raw_cache_sha256=expected,
        observer_proper_time=proper,R=float(R),U=float(U),Krad=float(Kr),Ktan=float(Kt),R_ll=float(Rll),
        geometric_mass=float(R*(1+U*U)/2),physical_volume_expansion=float(-Kr-2*Kt),
        sectional_components=sections,invariants={k:float(v) for k,v in invariants.items()},
        Weyl_squared_times_R4=float(invariants['Weyl_squared']*R**4),
        Weyl_squared_times_R6=float(invariants['Weyl_squared']*R**6),
        area_preserving_frame=p1,curvature_normalizing_frame=p15,
        scope='Finite old one-mode neck; actual saved three-level metric jets. No R=0 event or continuation inferred.')


def neck_conformal_calculation():
    paths=[ROOT/'results/summaries/nsc-atlas-valve-restart-v1.json',
           ROOT/'results/summaries/nsc-atlas-valve-conformal-v1.json',
           ROOT/'results/summaries/nsc-atlas-shape-transition-v1.json']
    restart,earlier,shape=(json.loads(p.read_bytes()) for p in paths)
    levels=[]
    for row in restart['levels']:
        levels.append(dict(points=row['points'],original_cut_admitted=row['final_cut_admitted'],
            **saved_neck_metric(row['raw_cache']['path'],row['raw_cache']['sha256'])))
    admitted=[r for r in levels if r['original_cut_admitted']]
    if not admitted:raise ValueError('no admitted old neck endpoint exists')
    cut=max(admitted,key=lambda r:r['points'])
    if abs(cut['sectional_components']['F'])>1e-10:
        raise ValueError('reflection condition needed for the proposed endpoint test')
    trends=[]
    rows=[v[k] for v in earlier['levels'] if v['points']==8001 for k in ('initial','endpoint') if k in v]
    rows += [v['last_actual_cut'] for v in restart['levels'] if v['points']==cut['points']]
    for r in rows:
        R=r['center_radius'];U=r['center_Rproper'];Q=r['center_rho']-r['center_radial_pressure']+r['center_tangential_pressure']
        E=(1+U*U)/(2*R*R)-Q/6
        trends.append(dict(observer_proper_time=r['center_proper_time'],R=R,Weyl_squared=48*E*E,
            Weyl_squared_times_R4=48*E*E*R**4,role='Einstein-source conditional, not an independent metric invariant.'))
    source=rows[-1]
    R=cut['R'];U=cut['U'];rho=source['center_rho'];pr=source['center_radial_pressure'];pt=source['center_tangential_pressure']
    needed=-(1+U*U)/(R*R)
    inside=shape['source_free_transition']['interior'];exit_density=inside['exit_density']
    conditions=dict(reflection=True,spatial_minimum=cut['R_ll']>=0,radial_stretch=cut['Krad']<=0,
        angular_collapse=cut['Ktan']>=0,nonnegative_radial_pressure=pr>=0)
    return dict(saved_levels=levels,finest_admitted=cut,source_inferred_trend=trends,
        current_source=dict(rho=rho,pr=pr,pt=pt,pressure_for_positive_acceleration=needed,
            original_force_requirement_scope='Present reflection sphere with nonzero inward speed; no replacement source inserted.'),
        conformal_convention='g_tilde=Omega^2*g; mixedC invariant, C_tilde^2=Omega^-4*C^2.',
        Omega_tending_zero='Amplifies a nonvanishing Weyl norm. If Omega*R->0, reflection sphere also obeys K_tilde>=4/(Omega*R)^4.',
        conditional_finite_area_obstruction=dict(current_conditions=conditions,
            premise='All current reflection/minimum/stretch/angular-collapse/pr>=0 conditions persist and R actually approaches0.',
            canonical_identity='E=(D-C)/3+kappa*sum(phi_l^2)/6 >=D/3=2m/(3R^3)',
            mass_balance='m_dot=-kappa*pr*R^2*R_dot/2 >=0; m>=m_cut>0',
            bound='C^2>=64*m_cut^2/(3R^6); ifOmega*R->L>0, C_tilde^2>=64*m_cut^2/[3*(Omega*R)^4*R^2]->infinity',
            persistence_proved=False,R_zero_approach_proved=False,
            scope='No finite nondegenerate C2 crossover of this reflection type while those signs persist; not an exclusion of every conformal frame or every gluing law.'),
        conditional_power_blowup=dict(power=1.5,
            premise='Additional positive finite limiting mass and vacuum-dominant m/R^3 sectional curvature.',
            limiting_sections_over_m=['-8','-1/2','2','1/2'],limiting_rho_over_m='9/2',
            limiting_pr_over_m='-3/2',limiting_pt_over_m='-21/2',limiting_w_iso='-5/3',
            limiting_Weyl_squared='48*m^2',limiting_Kretschmann='291*m^2',
            normal_clock='tau_tilde~-lnR/sqrt(2m0), infinite',
            caveat='Stronger Omega->infinity powers can bound curvature; this does not construct a finite CCC crossover or same-action child.'),
        physical_other_side=dict(glued=False,volume=None,energy=None,w=None,expanding=None,
            reason='No same-action conformal matching/extension has been derived. The tested p1.5 representative violates angular canonical NEC.'),
        saved_shape_comparison=dict(exit_density=exit_density,exit_pressure=-exit_density,exit_w=-1.,
            exit_reason='Exact regular pole has zero field velocities.',turn_density=inside['density_at_endpoint'],
            turn_pressure=inside['pressure_at_endpoint'],turn_w=inside['pressure_at_endpoint']/inside['density_at_endpoint'],
            same_epoch=False,matched_to_conformal_other_side=False),
        R_zero_event_reached=False,new_fields_or_states_added=False,
        input_sha256={p.relative_to(ROOT).as_posix():file_sha256(p) for p in paths})


def mcrit_relation():
    """Exact loading fold versus calibrated excitation density, no new state."""
    import sympy as sp
    import mpmath as mp
    x,m=sp.symbols('x m');P=sp.Poly(225*x**6-560*x+336,x)
    equality=sp.Poly(700*x**4-672*x**3-27,x)
    gcd=sp.gcd(P,equality)
    (lo,hi),multiplicity=next((bounds,k) for bounds,k in P.intervals(eps=sp.Rational(1,10**30))
        if bounds[0]>sp.Rational(99,100) and bounds[1]<1)
    A=3*m+sp.Rational(47,8);B=3*m+sp.Rational(51,8);T=sp.Rational(25,32)
    M0=sp.simplify(-2*A+sp.Rational(7,3)*B-4*T)
    Ms=sp.simplify(12*A-sp.Rational(56,3)*B+72*T)
    Mss=sp.simplify(-72*A+sp.Rational(448,3)*B-1296*T)
    E0=(1-sp.Rational(1,4))*(sp.Rational(3,2)-sp.Rational(1,2))**2/2
    with mp.workdps(40):
        xi=(mp.mpf(str(lo.p))/mp.mpf(str(lo.q))+mp.mpf(str(hi.p))/mp.mpf(str(hi.q)))/2
        critical=xi**3*(175*xi-168)/18;ratio=critical/(mp.mpf(3)/8)
        result=dict(Mcrit=float(critical),E0=float(E0),ratio=float(ratio),
            percent_above_E0=float(100*(ratio-1)),percent_E0_below_crit=float(100*(1-1/ratio)),
            fold_x=float(xi),fold_s=float(-mp.log(xi)/2),
            actual_density_at_fold=float(critical*xi**3))
    result.update(exact=dict(fold_polynomial=str(P.as_expr()),equality_polynomial=str(equality.as_expr()),
        gcd=str(gcd.as_expr()),isolating_interval=[str(lo),str(hi)],multiplicity=multiplicity,
        Mcrit='xi^3*(175xi-168)/18',ratio='(4/27)xi^3*(175xi-168)',E0=str(E0),
        matched_A=str(A),matched_B=str(B),matched_M_at0=str(M0),matched_Ms_at0=str(Ms),matched_Mss_at0=str(Mss),
        matched_fold_m=str(sp.solve(Ms,m)[0]),fiducial_Ms_at0=str(Ms.subs(m,E0)),fiducial_Mss_at0=str(Mss.subs(m,E0))),
        identity=False,independent_Parke_prediction=False,
        explanation='The imposed matching gives M(0)=E0 exactly; its nearby adiabatic fold gives the small parameter-dependent margin.',
        units='Geometrized density/loading coefficients with kappa=1,Lunit(0)=1; not geometric mass, total energy or Joules.',
        E0_origin='Original scalar prepared at the grandchild minimum: (1-epsilon)*(qparent-qgrandchild)^2/2. Parke imports this target; it does not produce it.',
        no_new_state_or_fit=True)
    return result


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--raw-directory',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();raw=(ROOT/args.raw_directory).resolve();out=(ROOT/args.output).resolve()
    if not raw.is_relative_to(ROOT/'.build') or raw.exists() or not out.is_relative_to(ROOT) or out.exists():
        raise ValueError('fresh raw folder under.build and fresh in-repository summary required')
    commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    script=Path(__file__).resolve();archived=subprocess.check_output(['git','show',commit+':'+script.relative_to(ROOT).as_posix()],cwd=ROOT)
    if hashlib.sha256(archived).hexdigest()!=file_sha256(script):raise ValueError('commit this readout before calculation')
    raw.mkdir(parents=True);started=time.perf_counter()
    item1=parent_drive_calculation(raw);item2=neck_conformal_calculation();item3=mcrit_relation()
    result=dict(schema='nsc-atlas-parent-valve-test-v1',scientific_owner='docs/nsc-atlas-valve-formation.md',
        parent_trigger=item1,neck_conformal=item2,Mcrit_reference=item3,
        scope=dict(new_fields=0,new_prepared_states=0,GR_or_PDE_runs=0,PDF_growth=0,
            parent_trigger_proved=False,physical_child_glued=False,full_goal_achieved=False),
        provenance=dict(readout_commit=commit,readout_sha256=file_sha256(script),inputs_frozen=True),
        compute=dict(readout_only=True,elapsed_seconds=time.perf_counter()-started))
    (raw/'readout-record.json').write_text(json.dumps(result,allow_nan=False,separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,allow_nan=False,separators=(',',':'))+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(output=str(out),scope=result['scope'],Mcrit_ratio=item3['ratio'],compute=result['compute']),indent=2))


if __name__=='__main__':main()
