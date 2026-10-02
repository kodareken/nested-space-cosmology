"""Bounded static, source-selected control of the existing common action.

A radial eigenfunction seeds the query; it is not a solution of the remaining
constraints. Every stationary residual uses the original discrete rates. The
rank-six covariance commutes with the actual retained Dirac Hamiltonian.
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import scipy
from scipy.linalg import eigh
from scipy.optimize import brentq, least_squares
from threadpoolctl import threadpool_limits

from . import nsc_nested_parent_child as nested
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin

LAB = Path(__file__).resolve().parents[2]
ROOT = LAB.parent
OUTPUT = LAB / "results/development/nsc-discovery-stationary-v1"
SCHEMA = "NSC-DISCOVERY-STATIONARY-v1"
NF = 32
CPU_LIMIT = 120.0
WEIGHTS = np.array([.75, .75, .5, .5, .25, .25])
BLOCKS = ("p_Q", "p_r", "p_chi", "lapse", "shift")
NUMERICAL_TARGET = 1e-7
MAX_BYTES = 2 << 20


def periodic_odd_derivative(points, length):
    """Full odd periodic collocation derivative, with no Nyquist ambiguity."""
    if points < 3 or points % 2 != 1 or length <= 0:
        raise ValueError("positive period and odd grid of at least three required")
    modes = np.fft.fftfreq(points) * points
    fourier = np.exp(2j * np.pi * np.arange(points)[:, None] * modes / points) / np.sqrt(points)
    return (fourier * (2j * np.pi * modes / length)) @ fourier.conj().T


def radial_operator(points=31, length=8., amplitude=.35, scale=0.):
    """L_s=-3 dxx-u_tilde_xx+s^2 exp(2 u_tilde), a radial-only seed."""
    if not np.isfinite(amplitude) or amplitude == 0 or not np.isfinite(scale) or scale < 0:
        raise ValueError("nonconstant finite cosine and nonnegative scale required")
    x = np.arange(points) * length / points
    u = amplitude * np.cos(2 * np.pi * x / length)
    D = periodic_odd_derivative(points, length).real
    D2 = D @ D
    matrix = -3 * D2 + np.diag(-(D2 @ u) + scale**2 * np.exp(2 * u))
    return (matrix + matrix.T) / 2, u, D


def radial_seed(points=31, length=8., amplitude=.35):
    """Unique zero of the strictly increasing lowest eigenvalue; mean(r)=1."""
    lowest = lambda s: float(eigh(radial_operator(points, length, amplitude, s)[0],
                                subset_by_index=(0, 0), eigvals_only=True)[0])
    low = lowest(0.)
    if low >= 0:
        raise RuntimeError("the finite radial seed did not resolve lambda_min(0)<0")
    high = 1.
    while lowest(high) <= 0:
        high *= 2
        if high > 64:
            raise RuntimeError("radial seed failed to bracket its eigenvalue zero")
    scale = brentq(lowest, 0., high, xtol=5e-15)
    matrix, u, D = radial_operator(points, length, amplitude, scale)
    values, vectors = eigh(matrix, subset_by_index=(0, 0))
    r = vectors[:, 0]
    if np.mean(r) < 0:
        r = -r
    if np.min(r) <= 0:
        raise RuntimeError("finite seed ground state is not positive at its nodes")
    r = r / np.mean(r)
    Q = scale * np.exp(u)
    identity = 3 * (D @ D @ r) / r + D @ D @ u - Q**2
    dx = length / points
    return {"Q": Q, "r_shape": r, "logQ": np.log(Q), "u_tilde": u,
            "scale": scale, "lambda_at_zero": low,
            "lambda_at_bracket": lowest(high), "bracket_high": high,
            "lowest_eigenvalue": float(values[0]), "radial_residual_max": maximum(identity),
            "integral_Q_square": float(dx * np.sum(Q**2)),
            "integral_3_logr_x_square": float(3 * dx * np.sum((D @ np.log(r))**2)),
            "integral_r_x_logQ_x": float(dx * np.dot(D @ r, D @ np.log(Q))),
            "minus_integral_Q_square_r": float(-dx * np.dot(Q**2, r)),
            "shape_normalization": "mean(r_shape)=1; physical radius amplitude remains free",
            "full_stationary_balance_claimed": False}


def maximum(values):
    return float(np.max(np.abs(values)))


def make_pair():
    """NF32 only; observer preparation is irrelevant to the stationary source."""
    columns = (np.eye(NF, 6, dtype=complex), np.zeros((NF, 6), dtype=complex))
    return nested.build_pair(NF, coarse_modes=1, child_details=2, source_layout="override",
                             columns_override=columns, occupations=WEIGHTS)


def nodal_state(grid, unknown):
    """All three nodal fields free, including mean(logQ) and mean(logr)."""
    unknown = np.asarray(unknown, dtype=float)
    if unknown.shape != (3 * grid.ng,) or not np.isfinite(unknown).all():
        raise ValueError("stationary coordinates must be three finite ng-vectors")
    u, v, chi = np.split(unknown, 3)
    zeros = np.zeros(grid.ng)
    columns = np.zeros((grid.nf, 6), dtype=complex)
    return coupling.CauchyState(np.exp(u), np.exp(v), chi.copy(), zeros.copy(),
                               zeros.copy(), zeros.copy(), columns.copy(), columns.copy())


def spectral_source(pair, state):
    """Six lowest positive REAL standing modes of the actual U† H_fine U.

    Equal-weight doublets allow arbitrary orthogonal rotations. A numerically
    unresolved degeneracy across unequal weights (including the empty tail)
    is exposed; no crossing is silently declared a differentiable branch.
    """
    H = nested.hamiltonian(pair, nested.encode_state(pair, state))
    real_gap = maximum(H.imag)
    hermitian_gap = maximum(H - H.conj().T)
    scale = max(1., maximum(H))
    if real_gap > 1e-10 * scale or hermitian_gap > 1e-10 * scale:
        raise ValueError("retained static Hamiltonian failed its real symmetric check")
    values, vectors = eigh((H.real + H.real.T) / 2)
    positive = np.flatnonzero(values > 1e-10 * scale)
    if positive.size < 7:
        raise ValueError("retained Hamiltonian has fewer than seven resolved positive modes")
    ids = positive[:6]
    columns = vectors[:, ids].astype(complex)
    state.phi0, state.phi1 = columns[:pair.grid.nf], columns[pair.grid.nf:]
    levels = values[ids]
    covariance = (columns * WEIGHTS) @ columns.conj().T
    separation = np.diff(values[positive[:7]])
    unequal = np.array([False, True, False, True, False, True])
    gap = float(np.min(separation[unequal]))
    return state, {"H": H, "eigenvalues": levels, "all_eigenvalues": values,
                   "columns": columns, "covariance": covariance,
                   "gram_gap": maximum(columns.conj().T @ columns - np.eye(6)),
                   "commutator_gap": maximum(H @ covariance - covariance @ H),
                   "retained_eigen_residual": maximum(H @ columns - columns * levels),
                   "real_hamiltonian_gap": real_gap, "hermitian_gap": hermitian_gap,
                   "unequal_occupation_gap": gap,
                   "unequal_occupation_degeneracy": bool(gap <= 1e-8 * scale),
                   "positive_empty_boundary_gap": float(separation[-1]),
                   "zero_to_occupied_gap": float(values[positive[0]])}


def evaluate(pair, unknown):
    """Raw stationary momenta and retained pre-gauge lapse/shift constraints."""
    state, spectral = spectral_source(pair, nodal_state(pair.grid, unknown))
    rate, bundle = galerkin.compose_fine_hamiltonian(pair.grid, state)
    fine, system, source = bundle["fine_state"], bundle["fine_system"], bundle["source"]
    constraints = coupling.constraint_residuals(system, fine, source)
    raw = {name: np.array(getattr(rate, name), copy=True) for name in BLOCKS[:3]}
    raw["lapse"] = galerkin.pull_geometry(pair.grid, constraints["hamilton"])
    raw["shift"] = galerkin.pull_geometry(pair.grid, constraints["momentum"])
    return state, spectral, raw, rate, bundle, constraints


def initial_unknown(pair, amplitude=.35):
    seed = radial_seed(pair.grid.ng, pair.grid.length, amplitude)
    D = periodic_odd_derivative(pair.grid.ng, pair.grid.length).real
    # Radial normalization is a seed convention only. Both physical means are
    # unconstrained optimization coordinates. Chi is a continuum seed only.
    chi = 2 * (D @ D @ seed["logQ"]) / seed["Q"]**2 - 2
    x = np.concatenate((seed["logQ"], np.log(seed["r_shape"]), chi))
    return x, seed


def residual_scales(raw):
    """Frozen seed RMS scales balance unlike units; all raw units are saved."""
    return {name: max(float(np.sqrt(np.mean(raw[name]**2))), 1e-3) for name in BLOCKS}


def residual_vector(pair, unknown, scales, *, evaluated=None):
    ev = evaluate(pair, unknown) if evaluated is None else evaluated
    phase = np.sin(2 * np.pi * pair.grid.xi_g / pair.grid.length)
    # One translation anchor; this fixes neither mean Q nor radius amplitude.
    anchor = np.sqrt(pair.grid.ng) * float(np.mean(unknown[:pair.grid.ng] * phase))
    return np.concatenate([ev[2][name] / scales[name] for name in BLOCKS] + [[anchor]])


def diagnostics(pair, unknown, evaluated=None):
    state, spectral, raw, rate, bundle, constraints = evaluate(pair, unknown) if evaluated is None else evaluated
    grid, fine, system, source = pair.grid, bundle["fine_state"], bundle["fine_system"], bundle["source"]
    levels = spectral["eigenvalues"]
    images = np.vstack((source["image0"], source["image1"]))
    lifted_columns = np.vstack((fine.phi0, fine.phi1))
    tail = images - lifted_columns * levels
    density = (np.abs(fine.phi0)**2 + np.abs(fine.phi1)**2) / grid.dx_q
    positive_carrier = system.multiplicity * np.sum(density * (WEIGHTS * levels), axis=1) / fine.Q
    source_tail = constraints["rho"] - positive_carrier
    # Hellmann--Feynman check uses the fixed occupation weights on nearby
    # eigensystems, with no extra multiplicity in the representative H.
    direction = np.cos(4 * np.pi * grid.xi_g / grid.length)
    eps = min(1e-6, .1 * float(np.min(fine.Q)) / max(1., maximum(grid.A_g @ direction)))
    def field_eigen_energy(delta):
        varied = state.copy()
        varied.Q = state.Q + delta * direction
        H = nested.hamiltonian(pair, nested.encode_state(pair, varied))
        eigen = eigh((H.real + H.real.T) / 2, eigvals_only=True)
        return float(system.multiplicity * np.dot(WEIGHTS, eigen[eigen > 0][:6]))
    owned = float(np.dot(source["force_Q"] + source["force_L"], grid.A_g @ direction))
    hf = {"step": eps, "owned_source_gradient": owned,
          "unequal_weight_branch_resolved": not spectral["unequal_occupation_degeneracy"]}
    try:
        fd = (field_eigen_energy(eps) - field_eigen_energy(-eps)) / (2 * eps)
        hf.update(performed=True, finite_difference=fd, absolute_gap=abs(fd - owned))
    except (ValueError, np.linalg.LinAlgError, coupling.PositiveChartExit) as error:
        # This auxiliary diagnostic cannot erase the primary residual result.
        hf.update(performed=False, reason=str(error), finite_difference=None, absolute_gap=None)
    D = system.derivative
    radial = 3 * (D @ D @ fine.r) / fine.r + D @ D @ np.log(fine.Q) - fine.Q**2
    report = {"raw_residual_max": {name: maximum(raw[name]) for name in BLOCKS},
              "zero_momentum_geometry_rate_max": {name: maximum(getattr(rate, name)) for name in ("Q", "r", "chi")},
              "full_lapse_residual_max": constraints["hamilton_max"],
              "full_shift_residual_max": constraints["momentum_max"],
              "quadrature_rho_min": float(np.min(constraints["rho"])),
              "positive_full_eigen_carrier_rho_min": float(np.min(positive_carrier)),
              "fine_eigen_tail_max": maximum(tail),
              "pointwise_source_tail_max": maximum(source_tail),
              "rho_positivity_certified": False,
              "quadrature_current_max": constraints["current_max"],
              "Q_min": float(np.min(fine.Q)), "r_min": float(np.min(fine.r)),
              "Q_mean": float(np.mean(fine.Q)), "r_mean": float(np.mean(fine.r)),
              "radial_continuum_identity_sampled_max": maximum(radial),
              "finite_CAR_eigenvalues": np.linalg.eigvalsh(spectral["covariance"]).tolist(),
              "spectral": {key: value for key, value in spectral.items() if not isinstance(value, np.ndarray)},
              "occupied_positive_eigenvalues": levels.tolist(),
              "hellmann_feynman": hf,
              "field_eigen_energy": float(system.multiplicity * np.dot(WEIGHTS, levels)),
              "owned_field_energy": float(coupling.field_energy(system, fine))}
    arrays = {"unknown": unknown.copy(), "eigenvalues": levels.copy(),
              "all_eigenvalues": spectral["all_eigenvalues"].copy(), "H": spectral["H"].copy(),
              "covariance": spectral["covariance"].copy(), "fine_eigen_tail": tail,
              "fine_rho": constraints["rho"], "fine_current": constraints["current"],
              "fine_lapse": constraints["hamilton"], "fine_shift": constraints["momentum"],
              "positive_carrier_rho": positive_carrier, "pointwise_source_tail": source_tail}
    arrays.update({"residual_" + name: raw[name] for name in BLOCKS})
    arrays.update({name: getattr(state, name).copy() for name in nested.STATE_NAMES})
    return report, arrays


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    names = ["src/recursive_horizons/nsc_discovery_stationary.py",
             "scripts/derive_nsc_discovery_stationary.py", "tests/test_nsc_discovery_stationary.py",
             "docs/nsc-discovery-stationary.md", "src/recursive_horizons/nsc_nested_parent_child.py",
             "src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
             "src/recursive_horizons/nsc_spherical_coupling.py",
             "src/recursive_horizons/nsc_spherical_feedback_action.py",
             "src/recursive_horizons/nsc_conformal_adm_source.py",
             "src/recursive_horizons/nsc_covariant_operator.py",
             "src/recursive_horizons/nsc_regulated.py"]
    return {"lab/" + name: sha256(LAB / name) for name in names}


def producer_commit_binding(commit, hashes):
    """Read-only Git pin, accepted only when every producer byte matches."""
    if commit is None:
        return None
    resolved = subprocess.check_output(["git", "rev-parse", "--verify", commit + "^{commit}"], cwd=ROOT, text=True).strip()
    for path, digest in hashes.items():
        data = subprocess.check_output(["git", "show", resolved + ":" + path], cwd=ROOT)
        if hashlib.sha256(data).hexdigest() != digest:
            raise ValueError("producer commit does not contain exact current bytes: " + path)
    return resolved


class BudgetReached(RuntimeError):
    pass


def run_stationary(*, solve=False, cpu_limit=CPU_LIMIT, max_nfev=1000,
                   amplitude=.35, producer_commit=None):
    """One bounded NF32 attempt, no trajectory or production refinement."""
    if not 0 < cpu_limit <= CPU_LIMIT or int(max_nfev) != max_nfev or max_nfev < 1:
        raise ValueError("CPU limit must be in (0,120], with positive integer max_nfev")
    started = time.process_time()
    before = source_hashes()
    commit = producer_commit_binding(producer_commit, before)
    input_hashes = {str(coupling.LOCKED_RECORD.relative_to(ROOT)): sha256(coupling.LOCKED_RECORD)}
    with threadpool_limits(limits=1):
        pair = make_pair()
        initial, seed = initial_unknown(pair, amplitude)
        initial_ev = evaluate(pair, initial)
        scales = residual_scales(initial_ev[2])
        initial_residual = residual_vector(pair, initial, scales, evaluated=initial_ev)
        best = {"unknown": initial.copy(), "cost": float(np.dot(initial_residual, initial_residual)),
                "evaluations": 0, "invalid_chart_evaluations": 0,
                "minimum_sampled_unequal_occupation_gap": initial_ev[1]["unequal_occupation_gap"],
                "unresolved_unequal_occupation_gap_evaluations": 0}
        def objective(unknown):
            if time.process_time() - started >= cpu_limit - min(2., .1 * cpu_limit):
                raise BudgetReached("bounded CPU limit reached")
            best["evaluations"] += 1
            try:
                ev = evaluate(pair, unknown)
                gap = ev[1]["unequal_occupation_gap"]
                best["minimum_sampled_unequal_occupation_gap"] = min(best["minimum_sampled_unequal_occupation_gap"], gap)
                best["unresolved_unequal_occupation_gap_evaluations"] += int(ev[1]["unequal_occupation_degeneracy"])
                residual = residual_vector(pair, unknown, scales, evaluated=ev)
            except coupling.PositiveChartExit:
                best["invalid_chart_evaluations"] += 1
                # Invalid interpolated positive chart is a numerical barrier,
                # never an accepted physical candidate or deleted source.
                fine = galerkin.prolong_state(pair.grid, nodal_state(pair.grid, unknown))
                penalty = 1e4 * (1 + max(0., -float(np.min(fine.Q)), -float(np.min(fine.r))))
                return np.full(5 * pair.grid.ng + 1, penalty)
            cost = float(np.dot(residual, residual))
            if cost < best["cost"]:
                best.update(unknown=unknown.copy(), cost=cost)
            return residual
        optimizer = {"performed": bool(solve), "termination": "seed preview", "success": False}
        if solve:
            ng = pair.grid.ng
            lower = np.r_[np.full(2 * ng, -8.), np.full(ng, -1e4)]
            upper = np.r_[np.full(2 * ng, 4.), np.full(ng, 1e4)]
            if np.any(initial <= lower) or np.any(initial >= upper):
                raise ValueError("radial seed lies outside the declared computational box")
            try:
                result = least_squares(objective, initial, bounds=(lower, upper),
                                       max_nfev=int(max_nfev), ftol=1e-10, xtol=1e-10,
                                       gtol=1e-10, x_scale="jac")
                optimizer.update(termination=result.message, success=bool(result.success),
                                 status=int(result.status), nfev=int(result.nfev),
                                 optimality=float(result.optimality))
            except BudgetReached as error:
                optimizer["termination"] = str(error)
            optimizer["computational_box"] = {"logQ_logr": [-8., 4.], "chi": [-1e4, 1e4]}
        selected = best["unknown"]
        final_ev = evaluate(pair, selected)
        final_residual = residual_vector(pair, selected, scales, evaluated=final_ev)
        measured, arrays = diagnostics(pair, selected, final_ev)
    current_hashes = source_hashes()
    if current_hashes != before or any(sha256(ROOT / p) != digest for p, digest in input_hashes.items()):
        raise RuntimeError("producer bytes or locked coefficients changed during stationary query")
    normalized_max = maximum(final_residual)
    satisfied = normalized_max <= NUMERICAL_TARGET and not measured["spectral"]["unequal_occupation_degeneracy"]
    status = "RADIAL_SEED_PREVIEW" if not solve else ("FINITE_STATIONARY_CANDIDATE" if satisfied else "UNSATISFIED_STATIONARY_RELATIONS")
    arrays["normalized_residual"] = final_residual
    arrays.update({"seed_Q": seed["Q"], "seed_r_shape": seed["r_shape"], "initial_unknown": initial})
    report = {"schema": SCHEMA, "status": status, "nf": NF, "ng": pair.grid.ng, "nq": pair.grid.nq,
              "period": pair.grid.length, "gauge": "conformal", "occupations": WEIGHTS.tolist(),
              "multiplicity": pair.grid.fine.multiplicity, "locked_coefficients": pair.grid.fine.coefficients,
              "source": "six lowest positive retained Dirac eigenmodes; real standing columns, fixed weights",
              "mean_Q_fixed": False, "radius_amplitude_fixed": False, "momenta": "zero nodal densities",
              "derived_static_relations": {
                  "radial_Euler_derivative": "E_r=-8 pi A r^3 R4",
                  "static_R4": "R4=2/r^2 [3 r_xx/r+(logQ)_xx-Q^2]",
                  "radial_equation": "3 r_xx/r+(logQ)_xx-Q^2=0",
                  "integrated_radial_equation": "integral Q^2=3 integral (logr_x)^2",
                  "weighted_integrated_equation": "integral r_x logQ_x=-integral Q^2 r<0",
                  "constant_positive_Q_static_periodic_radius": "excluded by integration; bounded or oscillatory dynamics not excluded"},
              "translation_anchor": "mean(logQ sin(2 pi x/period))=0",
              "seed": {key: value for key, value in seed.items() if not isinstance(value, np.ndarray)},
              "optimizer": dict(optimizer, actual_residual_evaluations=best["evaluations"],
                                invalid_chart_evaluations=best["invalid_chart_evaluations"],
                                minimum_sampled_unequal_occupation_gap=best["minimum_sampled_unequal_occupation_gap"],
                                unresolved_unequal_occupation_gap_evaluations=best["unresolved_unequal_occupation_gap_evaluations"],
                                eigenbranch_policy="ordered six lowest positive modes; no unobserved crossing exclusion"),
              "residual_scales": scales, "normalized_residual_max": normalized_max,
              "numerical_target": NUMERICAL_TARGET, "initial_normalized_residual_max": maximum(initial_residual),
              "measurements": measured, "cpu_seconds": time.process_time() - started,
              "cpu_limit_seconds": cpu_limit, "cpu_budget_scope": "whole query; solver reserves min(2 seconds, 10 percent) for fixed NF32 diagnostics",
              "cpu_budget_exceeded": bool(time.process_time() - started > cpu_limit),
              "source_hashes": current_hashes, "input_hashes": input_hashes,
              "producing_commit": commit, "numpy_version": np.__version__, "scipy_version": scipy.__version__,
              "scope": {"finite_discrete_code_domain": True, "continuum_stationary_solution": False,
                        "long_time_evolution_performed": False, "stability_assessed": False,
                        "nonexistence_inferred_from_optimizer": False,
                        "pointwise_source_positivity_promoted": False},
              "next_unsatisfied_relation": None if satisfied else max(measured["raw_residual_max"], key=lambda k: measured["raw_residual_max"][k] / scales[k])}
    return report, arrays


def output_paths(path=OUTPUT):
    stem = Path(path).expanduser().resolve()
    if stem.suffix in (".json", ".npz"):
        stem = stem.with_suffix("")
    temporary = Path("/tmp").resolve()
    if stem != OUTPUT.resolve() and temporary not in stem.parents:
        raise ValueError("stationary output must use the new v1 stem or a /tmp stem")
    return stem.with_suffix(".json"), stem.with_suffix(".npz")


def validate_new_output(path=OUTPUT):
    paths = output_paths(path)
    for target in paths:
        if target.exists():
            raise FileExistsError("refusing to replace " + str(target))
    return paths


def write_record(report, arrays, path=OUTPUT):
    """Exclusive new JSON+NPZ only; preserve all earlier evidence."""
    json_path, npz_path = validate_new_output(path)
    stream = io.BytesIO()
    np.savez_compressed(stream, **arrays)
    payload = stream.getvalue()
    if len(payload) > MAX_BYTES:
        raise ValueError("stationary payload exceeds two MiB")
    record = dict(report, payload=npz_path.name, payload_sha256=hashlib.sha256(payload).hexdigest(),
                  payload_bytes=len(payload))
    encoded = (json.dumps(record, indent=2, allow_nan=False) + "\n").encode()
    if len(encoded) > MAX_BYTES:
        raise ValueError("stationary JSON exceeds two MiB")
    json_path.parent.mkdir(parents=True, exist_ok=True)
    created = []
    try:
        for target, data in ((npz_path, payload), (json_path, encoded)):
            with target.open("xb") as handle:
                created.append(target)
                handle.write(data)
    except Exception:
        for target in created:
            target.unlink()
        raise
    return record


def check_record(path=OUTPUT):
    """Read-only producer/input/payload hashes and saved algebraic consistency."""
    json_path, npz_path = output_paths(path)
    record = json.loads(json_path.read_text())
    if (record.get("schema") != SCHEMA or record.get("nf") != NF or record.get("ng") != NF - 1
            or record.get("nq") != 4 * NF or record.get("period") != coupling.PERIOD
            or record.get("gauge") != "conformal" or record.get("occupations") != WEIGHTS.tolist()
            or record.get("locked_coefficients") != coupling.locked_coefficients()):
        raise ValueError("stationary record domain does not match")
    if record.get("source_hashes") != source_hashes():
        raise ValueError("stationary producer hashes changed")
    expected_input = {str(coupling.LOCKED_RECORD.relative_to(ROOT)): sha256(coupling.LOCKED_RECORD)}
    if record.get("input_hashes") != expected_input:
        raise ValueError("locked stationary input hash changed")
    producer_commit_binding(record.get("producing_commit"), record["source_hashes"])
    if record.get("payload") != npz_path.name or sha256(npz_path) != record.get("payload_sha256"):
        raise ValueError("stationary payload hash changed")
    if npz_path.stat().st_size != record.get("payload_bytes") or npz_path.stat().st_size > MAX_BYTES:
        raise ValueError("stationary payload size changed")
    with np.load(npz_path, allow_pickle=False) as arrays:
        for name in arrays.files:
            if not np.isfinite(arrays[name]).all():
                raise ValueError("nonfinite stationary array: " + name)
        ng = record["ng"]
        u, v, chi = np.split(arrays["unknown"], 3)
        if u.shape != (ng,) or maximum(arrays["Q"] - np.exp(u)) > 1e-12 or maximum(arrays["r"] - np.exp(v)) > 1e-12 or maximum(arrays["chi"] - chi) > 1e-12:
            raise ValueError("stationary coordinates and saved geometry disagree")
        for name in BLOCKS:
            if maximum(arrays["residual_" + name]) != record["measurements"]["raw_residual_max"][name]:
                raise ValueError("saved stationary residual summary disagrees: " + name)
        normalized = np.concatenate([arrays["residual_" + name] / record["residual_scales"][name]
                                     for name in BLOCKS] + [[np.sqrt(ng) * float(np.mean(
                                         u * np.sin(2 * np.pi * np.arange(ng) / ng)))]])
        if maximum(arrays["normalized_residual"] - normalized) > 1e-12:
            raise ValueError("saved normalized residual and raw blocks disagree")
        if maximum(arrays["normalized_residual"]) != record["normalized_residual_max"]:
            raise ValueError("normalized stationary residual summary disagrees")
        columns = np.vstack((arrays["phi0"], arrays["phi1"]))
        C = (columns * WEIGHTS) @ columns.conj().T
        if maximum(C - arrays["covariance"]) > 1e-12 or maximum(columns.conj().T @ columns - np.eye(6)) > 1e-10:
            raise ValueError("saved stationary CAR columns disagree")
        if np.min(arrays["eigenvalues"]) <= 0 or maximum(arrays["H"] @ columns - columns * arrays["eigenvalues"]) > 1e-8:
            raise ValueError("saved positive retained eigenmodes disagree")
        if float(np.min(arrays["fine_rho"])) != record["measurements"]["quadrature_rho_min"]:
            raise ValueError("saved quadrature source summary disagrees")
    return {"ok": True, "schema": SCHEMA, "status": record["status"], "bytes_written": 0, "solver_executed": False}
