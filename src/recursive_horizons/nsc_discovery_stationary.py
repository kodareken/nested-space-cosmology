"""Bounded static, source-selected control of the existing common action.

A radial eigenfunction seeds the query; it is not a solution of the remaining
constraints. Every stationary residual uses the original discrete rates. The
rank-six covariance commutes with the actual retained Dirac Hamiltonian.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import replace
from functools import lru_cache
import hashlib
import io
import json
from pathlib import Path
import re
import subprocess
import time

import numpy as np
import scipy
from scipy.linalg import eigh, solve
from scipy.optimize import brentq, least_squares
from scipy.special import iv
from threadpoolctl import threadpool_limits

from . import nsc_nested_parent_child as nested
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin

LAB = Path(__file__).resolve().parents[2]
ROOT = LAB.parent
OUTPUT = LAB / "results/development/nsc-discovery-stationary-v2"
LEGACY_OUTPUT = LAB / "results/development/nsc-discovery-stationary-v1"
SCHEMA = "NSC-DISCOVERY-STATIONARY-v2"
LEGACY_SCHEMA = "NSC-DISCOVERY-STATIONARY-v1"
CRITICAL_SCHEMA = "NSC-DISCOVERY-STATIONARY-CRITICAL-v3"
CRITICAL_OUTPUT = LAB / "results/development/nsc-discovery-stationary-critical-v3"
CRITICAL_HARMONIC = 8
CRITICAL_COORDINATES = 7
CRITICAL_SEED_AMPLITUDE = .16415733483383863
CRITICAL_MAX_ITERATIONS = 40
CRITICAL_MAX_TRIALS = 600
BALANCE_NF = 128
BALANCE_HARMONICS = (6, 8, 10, 12)
BALANCE_AMPLITUDES = (.025, .1, .35, .65)
BALANCE_ROOT_EVALUATIONS = 32
MAX_BALANCE_TRIALS = 64
NF = 32
CPU_LIMIT = 120.0
WEIGHTS = np.array([.75, .75, .5, .5, .25, .25])
BLOCKS = ("p_Q", "p_r", "p_chi", "lapse", "shift")
NUMERICAL_TARGET = 1e-7
MAX_BYTES = 16 << 20


def periodic_odd_derivative(points, length):
    """Full odd periodic collocation derivative, with no Nyquist ambiguity."""
    if points < 3 or points % 2 != 1 or length <= 0:
        raise ValueError("positive period and odd grid of at least three required")
    modes = np.fft.fftfreq(points) * points
    fourier = np.exp(2j * np.pi * np.arange(points)[:, None] * modes / points) / np.sqrt(points)
    return (fourier * (2j * np.pi * modes / length)) @ fourier.conj().T


@lru_cache(maxsize=8)
def _odd_operators(points, length):
    D = periodic_odd_derivative(points, length).real
    D2 = D @ D
    for values in (D, D2):
        values.setflags(write=False)
    return D, D2


@lru_cache(maxsize=16)
def _radial_kernel(points, length, harmonic):
    if int(harmonic) != harmonic or not 1 <= harmonic <= points // 2:
        raise ValueError("radial harmonic must be a resolved positive integer")
    x = np.arange(points) * length / points
    cosine = np.cos(2 * np.pi * harmonic * x / length)
    D, D2 = _odd_operators(points, length)
    for values in (cosine, D, D2):
        values.setflags(write=False)
    return cosine, D, D2


def radial_operator(points=31, length=8., amplitude=.35, scale=0., harmonic=1):
    """L_s=-3 dxx-u_tilde_xx+s^2 exp(2 u_tilde), a radial-only seed."""
    if not np.isfinite(amplitude) or amplitude == 0 or not np.isfinite(scale) or scale < 0:
        raise ValueError("nonconstant finite cosine and nonnegative scale required")
    cosine, D, D2 = _radial_kernel(points, float(length), harmonic)
    u = amplitude * cosine
    matrix = -3 * D2 + np.diag(-(D2 @ u) + scale**2 * np.exp(2 * u))
    return (matrix + matrix.T) / 2, u, D


def radial_seed(points=31, length=8., amplitude=.35, harmonic=1):
    """Unique zero of the strictly increasing lowest eigenvalue; mean(r)=1."""
    lowest = lambda s: float(eigh(radial_operator(points, length, amplitude, s, harmonic)[0],
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
    matrix, u, D = radial_operator(points, length, amplitude, scale, harmonic)
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
            "scale": scale, "harmonic": int(harmonic), "amplitude": float(amplitude), "lambda_at_zero": low,
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


def make_pair(nf=NF):
    """Full geometry frame; placeholders never supply a measured source."""
    if nf not in (NF, 64, BALANCE_NF):
        raise ValueError("stationary controls admit nf=32,64,128 only; nf64 is a reduced test domain")
    columns = (np.eye(nf, 6, dtype=complex), np.zeros((nf, 6), dtype=complex))
    return nested.build_pair(nf, coarse_modes=1, child_details=2, source_layout="override",
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
    weights = pair.grid.fine.occupations
    covariance = (columns * weights) @ columns.conj().T
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


def initial_unknown(pair, amplitude=.35, harmonic=1):
    seed = radial_seed(pair.grid.ng, pair.grid.length, amplitude, harmonic)
    D, _ = _odd_operators(pair.grid.ng, pair.grid.length)
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


def virial_diagnostic(pair, state, raw, bundle, constraints, spectral=None):
    """Exact p=0 SBP identity, with adjoint coarse/fine pairings.

    alpha int Q^2 chi^2 - int Q rho - 2 pi C_F flux^2 int Q^2
      = -.5 int r pdot_r - int chi pdot_chi - int Q (C_g+rho).
    This relates residuals; setting its left side to zero is only necessary
    for stationary balance and cannot replace any remaining equation.
    """
    if any(np.any(getattr(state, name) != 0) for name in nested.MOMENTUM_NAMES):
        raise ValueError("static virial diagnostic requires zero momenta")
    grid, fine, system = pair.grid, bundle["fine_state"], bundle["fine_system"]
    alpha = float(coupling.alpha_of(system.C_W))
    auxiliary = float(alpha * grid.dx_q * np.sum(fine.Q**2 * fine.chi**2))
    source = float(grid.dx_q * np.dot(fine.Q, constraints["rho"]))
    magnetic = float(2 * np.pi * system.C_F * system.flux**2 * grid.dx_q * np.sum(fine.Q**2))
    gap = auxiliary - source - magnetic
    fine_terms = {
        "radius": float(-.5 * grid.dx_q * np.dot(fine.r, bundle["unprojected_rates"][4])),
        "chi": float(-grid.dx_q * np.dot(fine.chi, bundle["unprojected_rates"][5])),
        "lapse": float(-grid.dx_q * np.dot(fine.Q, constraints["hamilton"]))}
    coarse_terms = {
        "radius": float(-.5 * grid.dx_g * np.dot(state.r, raw["p_r"])),
        "chi": float(-grid.dx_g * np.dot(state.chi, raw["p_chi"])),
        "lapse": float(-grid.dx_g * np.dot(state.Q, raw["lapse"]))}
    fine_rhs, coarse_rhs = sum(fine_terms.values()), sum(coarse_terms.values())
    owned_energy = float(coupling.field_energy(system, fine))
    eigen_energy = None if spectral is None else float(system.multiplicity * np.dot(system.occupations, spectral["eigenvalues"]))
    return {"auxiliary_term": auxiliary, "source_term": source, "magnetic_term": magnetic,
            "gap": gap, "fine_rhs_terms": fine_terms, "coarse_rhs_terms": coarse_terms,
            "fine_rhs": fine_rhs, "coarse_rhs": coarse_rhs,
            "fine_identity_defect": gap - fine_rhs, "coarse_identity_defect": gap - coarse_rhs,
            "adjoint_pairing_defect": coarse_rhs - fine_rhs,
            "source_owned_field_energy_gap": source - owned_energy,
            "source_spectral_field_energy_gap": None if eigen_energy is None else source - eigen_energy,
            "necessary_only": True, "stationary_balance_claimed": False}


def select_radius_affine(pair, state, spectral=None):
    """Select a^2 from the actual finite pQ rate at r=R and r=2R.

    Q, chi and the spectral covariance stay fixed. The third radius 1.5R
    verifies finite quadratic affinity and source independence separately.
    A negative least-squares root is a branch mismatch, not a positive radius.
    """
    grid = pair.grid
    trials = {}
    for amplitude in (1., 2., 1.5):
        changed = state.copy()
        changed.r = amplitude * state.r
        rate, bundle = galerkin.compose_fine_hamiltonian(grid, changed)
        fine, system = bundle["fine_state"], bundle["fine_system"]
        constraints = coupling.constraint_residuals(system, fine, bundle["source"])
        trials[amplitude] = (changed, rate, bundle, constraints)
    first, second, third = trials[1.], trials[2.], trials[1.5]
    coefficient = (second[1].p_Q - first[1].p_Q) / 3
    remainder = first[1].p_Q - coefficient
    denominator = float(np.dot(coefficient, coefficient))
    square = None if denominator == 0 else float(-np.dot(coefficient, remainder) / denominator)
    expected = remainder + 2.25 * coefficient
    source_gap = max(maximum(trials[a][2]["source"][name] - first[2]["source"][name])
                     for a in (2., 1.5) for name in ("force_L", "force_Q", "force_beta"))
    report = {"least_squares_amplitude_square": square,
              "finite_affinity_verification_amplitude": 1.5,
              "finite_affinity_defect_max": maximum(third[1].p_Q - expected),
              "source_change_max": source_gap,
              "fine_field_eigen_residual_max": None if spectral is None else maximum(
                  np.vstack((first[2]["source"]["image0"], first[2]["source"]["image1"]))
                  - np.vstack((first[2]["fine_state"].phi0, first[2]["fine_state"].phi1)) * spectral["eigenvalues"]),
              "coefficient_max": maximum(coefficient),
              "remainder_max": maximum(remainder), "fit_metric": "unweighted coarse p_Q nodal L2",
              "stationary_balance_claimed": False,
              "branch": "degenerate_radius_coefficient" if square is None else
                        ("nonpositive_radius_square_branch_mismatch" if square <= 0 else "positive_conditional_radius")}
    arrays = {"radius_pQ_coefficient": coefficient, "radius_pQ_remainder": remainder,
              "radius_third_pQ": third[1].p_Q.copy()}
    if square is not None and square > 0:
        amplitude = float(np.sqrt(square))
        changed = state.copy()
        changed.r = amplitude * state.r
        rate, bundle = galerkin.compose_fine_hamiltonian(grid, changed)
        constraints = coupling.constraint_residuals(bundle["fine_system"], bundle["fine_state"], bundle["source"])
        raw = {name: np.array(getattr(rate, name), copy=True) for name in BLOCKS[:3]}
        raw.update(lapse=galerkin.pull_geometry(grid, constraints["hamilton"]),
                   shift=galerkin.pull_geometry(grid, constraints["momentum"]))
        report.update(radius_amplitude=amplitude,
                      raw_residual_max={name: maximum(values) for name, values in raw.items()},
                      fine_momentum_residual_max={name: maximum(bundle["unprojected_rates"][i])
                                                  for name, i in (("p_Q", 3), ("p_r", 4), ("p_chi", 5))},
                      full_lapse_residual_max=constraints["hamilton_max"],
                      full_shift_residual_max=constraints["momentum_max"],
                      virial=virial_diagnostic(pair, changed, raw, bundle, constraints, spectral))
        arrays.update(radius_selected_r=changed.r.copy(), radius_selected_fine_r=bundle["fine_state"].r.copy(),
                      radius_selected_fine_lapse=constraints["hamilton"].copy())
        arrays.update({"radius_selected_residual_" + name: value for name, value in raw.items()})
        arrays.update({"radius_selected_fine_" + name: bundle["unprojected_rates"][index].copy()
                       for name, index in (("p_Q", 3), ("p_r", 4), ("p_chi", 5))})
    return report, arrays


def diagnostics(pair, unknown, evaluated=None, *, hellmann_feynman=True):
    state, spectral, raw, rate, bundle, constraints = evaluate(pair, unknown) if evaluated is None else evaluated
    grid, fine, system, source = pair.grid, bundle["fine_state"], bundle["fine_system"], bundle["source"]
    levels = spectral["eigenvalues"]
    weights = system.occupations
    images = np.vstack((source["image0"], source["image1"]))
    lifted_columns = np.vstack((fine.phi0, fine.phi1))
    tail = images - lifted_columns * levels
    density = (np.abs(fine.phi0)**2 + np.abs(fine.phi1)**2) / grid.dx_q
    positive_carrier = system.multiplicity * np.sum(density * (weights * levels), axis=1) / fine.Q
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
        return float(system.multiplicity * np.dot(weights, eigen[eigen > 0][:6]))
    owned = float(np.dot(source["force_Q"] + source["force_L"], grid.A_g @ direction))
    hf = {"step": eps, "owned_source_gradient": owned,
          "unequal_weight_branch_resolved": not spectral["unequal_occupation_degeneracy"]}
    try:
        if not hellmann_feynman:
            raise ValueError("auxiliary spectral finite difference omitted in bounded selection")
        fd = (field_eigen_energy(eps) - field_eigen_energy(-eps)) / (2 * eps)
        hf.update(performed=True, finite_difference=fd, absolute_gap=abs(fd - owned))
    except (ValueError, np.linalg.LinAlgError, coupling.PositiveChartExit) as error:
        # This auxiliary diagnostic cannot erase the primary residual result.
        hf.update(performed=False, reason=str(error), finite_difference=None, absolute_gap=None)
    D = system.derivative
    radial = 3 * (D @ D @ fine.r) / fine.r + D @ D @ np.log(fine.Q) - fine.Q**2
    fine_momenta = {name: bundle["unprojected_rates"][index]
                    for name, index in (("p_Q", 3), ("p_r", 4), ("p_chi", 5))}
    report = {"virial": virial_diagnostic(pair, state, raw, bundle, constraints, spectral),
              "fine_momentum_residual_max": {name: maximum(values) for name, values in fine_momenta.items()},
              "momentum_projection_defect_max": {
                  name: maximum(values - galerkin.prolong_geometry(grid, raw[name]))
                  for name, values in fine_momenta.items()},
              "fine_chi_seed_relation_defect_max": maximum(fine.Q**2 * (2 + fine.chi) - 2 * (D @ (D @ np.log(fine.Q)))),
              "raw_residual_max": {name: maximum(raw[name]) for name in BLOCKS},
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
              "field_eigen_energy": float(system.multiplicity * np.dot(weights, levels)),
              "owned_field_energy": float(coupling.field_energy(system, fine))}
    arrays = {"unknown": unknown.copy(), "eigenvalues": levels.copy(),
              "all_eigenvalues": spectral["all_eigenvalues"].copy(), "H": spectral["H"].copy(),
              "covariance": spectral["covariance"].copy(), "fine_eigen_tail": tail,
              "fine_rho": constraints["rho"], "fine_current": constraints["current"],
              "fine_lapse": constraints["hamilton"], "fine_shift": constraints["momentum"],
              "positive_carrier_rho": positive_carrier, "pointwise_source_tail": source_tail}
    arrays.update({"fine_" + name: values for name, values in fine_momenta.items()})
    arrays.update({"fine_Q": fine.Q.copy(), "fine_r": fine.r.copy(), "fine_chi": fine.chi.copy()})
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
    report = {"schema": SCHEMA, "mode": "stationary", "status": status, "nf": NF, "ng": pair.grid.ng, "nq": pair.grid.nq,
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


def run_saved_diagnostic(path=LEGACY_OUTPUT, *, producer_commit=None):
    """Re-measure sealed arrays, with no eigenstate selection or optimization."""
    started = time.process_time()
    authenticated = check_record(path)
    json_path, npz_path = output_paths(path)
    old = json.loads(json_path.read_text())
    if old.get("mode") in ("balance_seed", "critical_shape"):
        raise ValueError("saved stationary diagnosis requires a single NF32 stationary iterate")
    before = source_hashes()
    commit = producer_commit_binding(producer_commit, before)
    owned = {"lab/src/recursive_horizons/nsc_discovery_stationary.py",
             "lab/scripts/derive_nsc_discovery_stationary.py", "lab/tests/test_nsc_discovery_stationary.py",
             "lab/docs/nsc-discovery-stationary.md"}
    for name, digest in old["source_hashes"].items():
        if name not in owned and before[name] != digest:
            raise ValueError("saved diagnosis requires unchanged numerical owner: " + name)
    if old["input_hashes"] != {str(coupling.LOCKED_RECORD.relative_to(ROOT)): sha256(coupling.LOCKED_RECORD)}:
        raise ValueError("saved diagnosis requires unchanged locked input")
    with np.load(npz_path, allow_pickle=False) as source:
        stored = {name: source[name].copy() for name in source.files}
    with threadpool_limits(limits=1):
        pair = make_pair()
        state = coupling.CauchyState(*(stored[name] for name in nested.STATE_NAMES))
        rate, bundle = galerkin.compose_fine_hamiltonian(pair.grid, state)
        constraints = coupling.constraint_residuals(bundle["fine_system"], bundle["fine_state"], bundle["source"])
        raw = {name: np.array(getattr(rate, name), copy=True) for name in BLOCKS[:3]}
        raw.update(lapse=galerkin.pull_geometry(pair.grid, constraints["hamilton"]),
                   shift=galerkin.pull_geometry(pair.grid, constraints["momentum"]))
        spectral = dict(old["measurements"]["spectral"], H=stored["H"],
                        eigenvalues=stored["eigenvalues"], all_eigenvalues=stored["all_eigenvalues"],
                        covariance=stored["covariance"], columns=np.vstack((state.phi0, state.phi1)))
        evaluated = state, spectral, raw, rate, bundle, constraints
        measured, arrays = diagnostics(pair, stored["unknown"], evaluated, hellmann_feynman=False)
        normalized = residual_vector(pair, stored["unknown"], old["residual_scales"], evaluated=evaluated)
    arrays["normalized_residual"] = normalized
    for name in ("seed_Q", "seed_r_shape", "initial_unknown"):
        if name in stored:
            arrays[name] = stored[name]
    if source_hashes() != before:
        raise RuntimeError("stationary diagnostic producers changed during read-only measurement")
    report = dict(old)
    for name in ("payload", "payload_sha256", "payload_bytes"):
        report.pop(name, None)
    report.update(schema=SCHEMA, mode="stationary", status="SAVED_STATIONARY_DIAGNOSTIC",
                  source_hashes=before, producing_commit=commit, measurements=measured,
                  normalized_residual_max=maximum(normalized), cpu_seconds=time.process_time() - started,
                  parent_cpu_seconds=old.get("cpu_seconds"), cpu_budget_exceeded=False,
                  optimizer=dict(old["optimizer"], executed_now=False),
                  parent_record={"json": str(json_path.relative_to(ROOT)) if ROOT in json_path.parents else str(json_path),
                                 "json_sha256": sha256(json_path), "payload_sha256": sha256(npz_path),
                                 "producing_commit": old.get("producing_commit"), "status": old["status"]},
                  historical_authentication=authenticated["source_authentication"],
                  next_unsatisfied_relation=max(measured["raw_residual_max"],
                                                key=lambda key: measured["raw_residual_max"][key] / old["residual_scales"][key]))
    return report, arrays


def balance_preflight(cpu_limit=CPU_LIMIT):
    """Declared finite work/memory forecast; no numerical trial is executed."""
    if not 0 < cpu_limit <= CPU_LIMIT:
        raise ValueError("balance selection CPU limit must be in (0,120]")
    # Estimates name the actual dense matrices. Measured cost of the first
    # trial then controls subsequent admission; this is not a timing claim.
    nf, ng, nq = BALANCE_NF, BALANCE_NF - 1, 4 * BALANCE_NF
    return {"nf": nf, "ng": ng, "nq": nq, "harmonics": list(BALANCE_HARMONICS),
            "endpoint_amplitudes": list(BALANCE_AMPLITUDES),
            "endpoint_trials": len(BALANCE_HARMONICS) * len(BALANCE_AMPLITUDES),
            "maximum_actual_trials": MAX_BALANCE_TRIALS,
            "maximum_brent_iterations_per_bracket": BALANCE_ROOT_EVALUATIONS,
            "cpu_limit_seconds": cpu_limit, "reserve_cpu_seconds": min(5., .1 * cpu_limit),
            "dense_operator_bytes_estimate": int(16 * (nq*nq + nq*nf + 4*nf*nf) + 8*nq*ng),
            "retained_eigendecomposition_dimension": 2 * nf,
            "initial_cost_forecast": "uncalibrated; first actual trial forecasts remaining endpoint work",
            "cost_admission": "fixed dense operator; omit further trials when measured forecast exceeds remaining budget",
            "cache": "odd derivative/harmonic kernels and exact (harmonic,amplitude) trial keys",
            "optimization": "none; amplitude Brent only after measured fine virial sign change",
            "radius_selection": "optional actual discrete quadratic p_Q amplitude fit",
            "payload_limit_bytes": MAX_BYTES}


def balance_trial(pair, harmonic, amplitude):
    """One radial seed measured through actual fine Q, chi and spectral source."""
    unknown, seed = initial_unknown(pair, amplitude, harmonic)
    evaluated = evaluate(pair, unknown)
    measured, arrays = diagnostics(pair, unknown, evaluated, hellmann_feynman=False)
    # Compact per-trial carriers preserve every raw measurement. H and C are
    # reproducible from the geometry/columns; dense duplicates would exhaust
    # the bounded payload over a root search. Their measured gaps remain JSON.
    arrays = {name: value for name, value in arrays.items() if name not in ("H", "covariance")}
    arrays.update(seed_Q=seed["Q"], seed_r_shape=seed["r_shape"])
    return {"harmonic": int(harmonic), "amplitude": float(amplitude), "status": "MEASURED",
            "seed": {key: value for key, value in seed.items() if not isinstance(value, np.ndarray)},
            "measurements": measured, "virial_gap": measured["virial"]["gap"]}, arrays


def sign_change_brackets(rows):
    """Only adjacent measured endpoints of the same declared harmonic."""
    brackets = []
    for harmonic in BALANCE_HARMONICS:
        group = sorted((row for row in rows if row["harmonic"] == harmonic and row["status"] == "MEASURED"),
                       key=lambda row: row["amplitude"])
        # A failed intermediate endpoint must not be skipped to manufacture
        # a sign bracket across an unread part of the declared family.
        for first, second in zip(group, group[1:]):
            i = BALANCE_AMPLITUDES.index(first["amplitude"])
            j = BALANCE_AMPLITUDES.index(second["amplitude"])
            if j == i + 1 and first["virial_gap"] * second["virial_gap"] < 0:
                brackets.append((first, second))
    return brackets


def run_balance_seed(*, execute=False, cpu_limit=CPU_LIMIT, radius_select=False,
                     producer_commit=None):
    """One NF128 discriminating selection, capped at 120 aggregate CPU seconds.

    No joint geometry Jacobian solve or trajectory is performed. All measured,
    failed, skipped, nonbracketing and root trials are retained in order.
    """
    plan = balance_preflight(cpu_limit)
    started = time.process_time()
    before = source_hashes()
    commit = producer_commit_binding(producer_commit, before)
    input_hashes = {str(coupling.LOCKED_RECORD.relative_to(ROOT)): sha256(coupling.LOCKED_RECORD)}
    report = {"schema": SCHEMA, "mode": "balance_seed", "status": "BALANCE_SEED_PREFLIGHT",
              "nf": BALANCE_NF, "ng": BALANCE_NF - 1, "nq": 4 * BALANCE_NF,
              "period": coupling.PERIOD, "gauge": "conformal", "occupations": WEIGHTS.tolist(),
              "multiplicity": 4 * coupling.KAPPA, "locked_coefficients": coupling.locked_coefficients(),
              "preflight": plan, "source_hashes": before, "input_hashes": input_hashes,
              "producing_commit": commit, "numpy_version": np.__version__, "scipy_version": scipy.__version__,
              "trials": [], "endpoint_assessment": [], "brackets": [], "roots": [],
              "selected_trial_index": None, "radius_selection": {"performed": False},
              "scope": {"necessary_virial_selection_only": True, "stationary_balance_claimed": False,
                        "joint_jacobian_optimization_performed": False, "long_time_evolution_performed": False,
                        "stability_assessed": False, "nonexistence_inferred_from_nonbracket": False,
                        "pointwise_source_positivity_promoted": False}}
    payload = {}
    if not execute:
        report.update(cpu_seconds=time.process_time() - started, actual_trials=0,
                      cpu_budget_exceeded=False, numerical_selection_executed=False)
        return report, payload
    cache = {}
    trial_arrays = {}
    reserve = plan["reserve_cpu_seconds"]
    cutoff = cpu_limit - reserve
    forecast_cost = None
    payload_forecast_limit = MAX_BALANCE_TRIALS
    def remaining():
        return cutoff - (time.process_time() - started)
    def trial(harmonic, amplitude, kind):
        nonlocal forecast_cost, payload_forecast_limit
        key = (int(harmonic), float(amplitude))
        if key in cache:
            row = report["trials"][cache[key]]
            if row["status"] != "MEASURED":
                raise ValueError(row["reason"])
            return row
        if len(report["trials"]) >= min(MAX_BALANCE_TRIALS, payload_forecast_limit):
            raise BudgetReached("balance selection reached its declared trial/payload admission cap")
        predicted = 0. if forecast_cost is None else forecast_cost
        if remaining() <= predicted:
            raise BudgetReached("measured trial cost exceeds remaining aggregate CPU allowance")
        index = len(report["trials"])
        began = time.process_time()
        try:
            row, arrays = balance_trial(pair, harmonic, amplitude)
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
            row = {"harmonic": int(harmonic), "amplitude": float(amplitude), "status": "FAILED", "reason": str(error)}
            arrays = {}
        cost = time.process_time() - began
        forecast_cost = max(cost, forecast_cost or 0.)
        row.update(index=index, kind=kind, cpu_seconds=cost)
        if arrays:
            per_trial_bytes = sum(np.asarray(values).nbytes for values in arrays.values())
            # Uncompressed upper estimate leaves space for JSON/zip metadata
            # and the optional radius diagnostics. Compressibility is unused.
            payload_forecast_limit = min(payload_forecast_limit,
                                         max(1, (MAX_BYTES - (1 << 20)) // (per_trial_bytes + 8192)))
            row["uncompressed_array_bytes"] = per_trial_bytes
            prefix = f"trial{index:03d}_"
            payload.update({prefix + name: values for name, values in arrays.items()})
            row["array_prefix"] = prefix
            trial_arrays[index] = arrays
        cache[key] = index
        report["trials"].append(row)
        if row["status"] != "MEASURED":
            raise ValueError(row["reason"])
        return row
    with threadpool_limits(limits=1):
        if remaining() <= 0:
            report["status"] = "BALANCE_SEED_BUDGET_EXHAUSTED"
        else:
            pair = make_pair(BALANCE_NF)
            construction_cpu = time.process_time() - started
            report["preflight"]["matrix_and_frame_construction_cpu_seconds"] = construction_cpu
            endpoints = []
            for harmonic in BALANCE_HARMONICS:
                for amplitude in BALANCE_AMPLITUDES:
                    try:
                        row = trial(harmonic, amplitude, "endpoint")
                        endpoints.append(row)
                        report["endpoint_assessment"].append({"harmonic": harmonic, "amplitude": amplitude,
                                                              "status": row["status"], "trial_index": row["index"]})
                    except ValueError as error:
                        failed = report["trials"][-1]
                        endpoints.append(failed)
                        report["endpoint_assessment"].append({"harmonic": harmonic, "amplitude": amplitude,
                                                              "status": "FAILED", "trial_index": failed["index"], "reason": str(error)})
                    except BudgetReached as error:
                        report["endpoint_assessment"].append({"harmonic": harmonic, "amplitude": amplitude,
                                                              "status": "NOT_ADMITTED", "reason": str(error)})
                    if forecast_cost is not None and "measured_first_trial_cpu_seconds" not in report["preflight"]:
                        report["preflight"].update(measured_first_trial_cpu_seconds=forecast_cost,
                                                   forecast_endpoint_cpu_seconds=construction_cpu + plan["endpoint_trials"] * forecast_cost,
                                                   forecast_worst_case_trial_cpu_seconds=construction_cpu + MAX_BALANCE_TRIALS * forecast_cost)
            for first, second in sign_change_brackets(endpoints):
                harmonic = first["harmonic"]
                bracket = {"harmonic": harmonic, "amplitudes": [first["amplitude"], second["amplitude"]],
                           "endpoint_trial_indices": [first["index"], second["index"]],
                           "endpoint_virial_gaps": [first["virial_gap"], second["virial_gap"]]}
                report["brackets"].append(bracket)
                root = dict(bracket, status="UNRESOLVED")
                try:
                    amplitude = brentq(lambda value: trial(harmonic, value, "brent")["virial_gap"],
                                       first["amplitude"], second["amplitude"], xtol=1e-10,
                                       maxiter=BALANCE_ROOT_EVALUATIONS)
                    row = trial(harmonic, amplitude, "brent")
                    root.update(status="VIRIAL_SELECTED_SHAPE", amplitude=amplitude,
                                trial_index=row["index"], virial_gap=row["virial_gap"],
                                stationary_balance_claimed=False)
                except (BudgetReached, ValueError, RuntimeError) as error:
                    root.update(reason=str(error))
                report["roots"].append(root)
            successful = [root for root in report["roots"] if root["status"] == "VIRIAL_SELECTED_SHAPE"]
            measured = [row for row in report["trials"] if row["status"] == "MEASURED"]
            if successful:
                chosen = successful[0]["trial_index"]
                report["status"] = "VIRIAL_SELECTED_SHAPE_ONLY"
                report["selection_policy"] = "first selected harmonic in the predeclared order; every other root retained"
            elif measured:
                chosen = min(measured, key=lambda row: abs(row["virial_gap"]))["index"]
                report["status"] = "NO_COMPLETED_VIRIAL_BRACKET"
                report["selection_policy"] = "smallest measured absolute virial gap for diagnostic reporting only"
            else:
                chosen = None
                report["status"] = "NO_MEASURED_BALANCE_SEED"
            report["selected_trial_index"] = chosen
            if chosen is not None:
                row = report["trials"][chosen]
                report["measurements"] = row["measurements"]
                if radius_select:
                    if remaining() <= (forecast_cost or 0.):
                        report["radius_selection"] = {"performed": False, "reason": "remaining CPU forecast does not admit auxiliary radius fit"}
                    else:
                        arrays = trial_arrays[chosen]
                        state = coupling.CauchyState(*(arrays[name].copy() for name in nested.STATE_NAMES))
                        try:
                            radius_report, radius_arrays = select_radius_affine(pair, state, {"eigenvalues": arrays["eigenvalues"]})
                            report["radius_selection"] = dict(radius_report, performed=True, trial_index=chosen)
                            payload.update(radius_arrays)
                        except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
                            report["radius_selection"] = {"performed": False, "reason": str(error)}
    after = source_hashes()
    if after != before or any(sha256(ROOT / path) != digest for path, digest in input_hashes.items()):
        raise RuntimeError("stationary producer bytes or locked inputs changed during balance selection")
    elapsed = time.process_time() - started
    report.update(cpu_seconds=elapsed, actual_trials=len(report["trials"]), numerical_selection_executed=True,
                  cpu_limit_seconds=cpu_limit, cpu_budget_exceeded=bool(elapsed > cpu_limit),
                  measured_trial_cost_forecast_seconds=forecast_cost,
                  effective_payload_trial_cap=payload_forecast_limit)
    return report, payload


class NonAdmissible(ValueError):
    """A measured trial left the stated source/positive-geometry branch."""
    def __init__(self, reason, **details):
        super().__init__(reason)
        self.details = details


def critical_seed_coefficients(count=CRITICAL_COORDINATES, amplitude=CRITICAL_SEED_AMPLITUDE):
    """Truncated exp(b cos) / I0(b), with its mean exactly one."""
    return 2 * iv(np.arange(1, count+1), amplitude) / iv(0, amplitude)


def critical_context(nf=BALANCE_NF, harmonic=CRITICAL_HARMONIC, count=CRITICAL_COORDINATES):
    """Original fine operators and their exact projected radial kinetic term."""
    pair = make_pair(nf)
    grid = pair.grid
    if not 1 <= harmonic or int(harmonic) != harmonic or harmonic * count > grid.ng // 2:
        raise ValueError("all critical cosine coordinates must fit the odd geometry band")
    cosine = np.cos(2*np.pi*grid.xi_g[:, None]/grid.length * harmonic * np.arange(1, count+1))
    D, U = grid.fine.derivative, grid.A_g
    kinetic = -3 * galerkin.pull_geometry(grid, D @ (D @ U))
    return {"pair": pair, "cosine": cosine, "kinetic": (kinetic+kinetic.T)/2,
            "harmonic": int(harmonic), "count": int(count)}


def critical_shape_evaluate(context, coefficients):
    """Exact finite envelopes and analytic eta gradient; physical eta is fixed.

    The preparation eta changes between candidates, not inside the action
    force. Its gradient is derived from those fixed-eta forces only after
    eliminating the radial eigenvalue and auxiliary critical equations.
    """
    base_pair = context["pair"]
    grid, count = base_pair.grid, context["count"]
    coefficients = np.asarray(coefficients, dtype=float)
    if coefficients.shape != (count,) or not np.isfinite(coefficients).all():
        raise NonAdmissible("critical shape coordinates are not finite")
    S = 1 + context["cosine"] @ coefficients
    fineS = galerkin.prolong_geometry(grid, S)
    if np.min(fineS) <= 0 or np.min(S) <= 0:
        raise NonAdmissible("shape left the actual fine positive chart", fine_S_min=float(np.min(fineS)))
    D, U = grid.fine.derivative, grid.A_g
    logarithmic_jet = D @ ((D @ fineS) / fineS)
    shape_mass = galerkin.pull_geometry(grid, fineS[:, None]**2 * U)
    base_operator = context["kinetic"] - galerkin.pull_geometry(grid, logarithmic_jet[:, None] * U)
    base_operator, shape_mass = (base_operator+base_operator.T)/2, (shape_mass+shape_mass.T)/2
    lowest = lambda scale: float(eigh(base_operator+scale**2*shape_mass,
                                     subset_by_index=(0,0), eigvals_only=True)[0])
    low = lowest(0.)
    if low >= 0:
        raise NonAdmissible("no resolved negative radial eigenvalue at zero scale", lambda_at_zero=low)
    high = 1.
    while lowest(high) <= 0:
        high *= 2
        if high > 64:
            raise NonAdmissible("radial scale bracket exceeded its bounded search")
    scale = brentq(lowest, 0., high, xtol=2e-13)
    operator = base_operator + scale**2 * shape_mass
    eigenvalues, eigenvectors = eigh(operator, subset_by_index=(0,0))
    R = eigenvectors[:, 0]
    if np.mean(R) < 0:
        R = -R
    R /= np.mean(R)
    fineR = galerkin.prolong_geometry(grid, R)
    if np.min(R) <= 0 or np.min(fineR) <= 0:
        raise NonAdmissible("radial ground state is not positive on the fine carrier", fine_R_min=float(np.min(fineR)))
    Q, fineQ = scale*S, scale*fineS
    chi_matrix = scale**2 * shape_mass
    chi_rhs = 2 * galerkin.pull_geometry(grid, logarithmic_jet) - 2 * galerkin.pull_geometry(grid, fineQ**2)
    chi = solve(chi_matrix, chi_rhs, assume_a="pos")
    finechi = galerkin.prolong_geometry(grid, chi)
    zeros, columns = np.zeros(grid.ng), np.zeros((grid.nf,6), dtype=complex)
    state = coupling.CauchyState(Q.copy(), R.copy(), chi.copy(), zeros.copy(), zeros.copy(), zeros.copy(), columns.copy(), columns.copy())
    state, base_spectral = spectral_source(base_pair, state)
    base_energy = float(grid.fine.multiplicity * np.dot(WEIGHTS, base_spectral["eigenvalues"]))
    auxiliary = float(coupling.alpha_of(grid.fine.C_W) * grid.dx_q * np.sum(fineQ**2 * finechi**2))
    magnetic = float(2*np.pi*grid.fine.C_F*grid.fine.flux**2 * grid.dx_q * np.sum(fineQ**2))
    eta = (auxiliary-magnetic) / base_energy
    if not 0 < eta <= 4/3:
        raise NonAdmissible("selected source eta is outside the CAR preparation branch",
                            eta=float(eta), auxiliary_term=auxiliary, magnetic_term=magnetic, base_source_energy=base_energy)
    occupations = eta * WEIGHTS
    fine_system = replace(grid.fine, occupations=occupations)
    pair = replace(base_pair, grid=replace(grid, fine=fine_system), weights=occupations)
    # Reuse the actual eigenvectors; changing weights does not change H(Q).
    spectral = dict(base_spectral, covariance=eta*base_spectral["covariance"],
                    commutator_gap=eta*base_spectral["commutator_gap"])
    first_rate, first_bundle = galerkin.compose_fine_hamiltonian(pair.grid, state)
    double = state.copy()
    double.r *= 2
    double_rate, double_bundle = galerkin.compose_fine_hamiltonian(pair.grid, double)
    PR = (double_rate.p_Q-first_rate.p_Q)/3
    remainder = first_rate.p_Q-PR
    pairing = lambda first, second: float(grid.dx_g * np.dot(first, second))
    denominator = pairing(Q, PR)
    denominator_identity = float(16*np.pi*grid.fine.A*grid.dx_q*np.sum(fineQ**2*fineR**2))
    if denominator <= 0:
        raise NonAdmissible("normal radius denominator is nonpositive", denominator=denominator)
    amplitude_square = -pairing(Q, remainder)/denominator
    dQ0 = scale * context["cosine"]
    mu = grid.dx_g * (dQ0.T @ PR) / denominator
    tangents = dQ0 - Q[:, None] * mu[None, :]
    gradient = grid.dx_g * (tangents.T @ remainder) / base_energy
    if amplitude_square <= 0:
        raise NonAdmissible("normal radius scale has no positive square", eta=float(eta),
                            radius_amplitude_square=float(amplitude_square), gradient=gradient.tolist())
    state.r *= np.sqrt(amplitude_square)
    rate, bundle = galerkin.compose_fine_hamiltonian(pair.grid, state)
    constraints = coupling.constraint_residuals(bundle["fine_system"], bundle["fine_state"], bundle["source"])
    raw = {name: np.array(getattr(rate, name), copy=True) for name in BLOCKS[:3]}
    raw.update(lapse=galerkin.pull_geometry(grid, constraints["hamilton"]),
               shift=galerkin.pull_geometry(grid, constraints["momentum"]))
    unknown = np.r_[np.log(state.Q), np.log(state.r), state.chi]
    source_change = max(maximum(double_bundle["source"][name]-first_bundle["source"][name])
                        for name in ("force_L", "force_Q", "force_beta"))
    return {"coefficients": coefficients.copy(), "eta": float(eta), "gradient": gradient,
            "pair": pair, "state": state, "spectral": spectral, "raw": raw, "rate": rate,
            "bundle": bundle, "constraints": constraints, "unknown": unknown,
            "PR": PR, "remainder": remainder, "tangents": tangents, "mu": mu,
            "R": R.copy(), "scale": float(scale), "radius_amplitude_square": float(amplitude_square),
            "reduction": {"lambda_at_zero": low, "selected_lambda": float(eigenvalues[0]),
                          "radial_operator_residual_max": maximum(operator @ R),
                          "auxiliary_projected_linear_residual_max": maximum(chi_matrix @ chi-chi_rhs),
                          "normal_radius_denominator": denominator,
                          "normal_radius_denominator_identity": denominator_identity,
                          "normal_radius_denominator_defect": denominator-denominator_identity,
                          "tangent_radial_pairing_max": maximum(grid.dx_g * tangents.T @ PR),
                          "normal_Q_force_pairing": pairing(Q, raw["p_Q"]),
                          "fixed_eta_source_change_under_radius_scaling": source_change,
                          "base_source_energy": base_energy, "auxiliary_term": auxiliary,
                          "magnetic_term": magnetic, "fine_S_min": float(np.min(fineS)),
                          "fine_R_min": float(np.min(fineR)), "physical_force_eta_held_fixed": True}}


def critical_gradient_check(context, coefficients, *, evaluator=None, steps=(1e-5, 5e-6)):
    """Two-step finite differences of selected eta, separate from physical forces."""
    evaluate_point = critical_shape_evaluate if evaluator is None else evaluator
    center = evaluate_point(context, coefficients)
    rows = []
    for step in steps:
        finite_difference = []
        for column in range(context["count"]):
            delta = np.zeros(context["count"])
            delta[column] = step
            plus = evaluate_point(context, np.asarray(coefficients)+delta)
            minus = evaluate_point(context, np.asarray(coefficients)-delta)
            finite_difference.append((plus["eta"]-minus["eta"])/(2*step))
        fd = np.asarray(finite_difference)
        rows.append({"step": step, "finite_difference": fd.tolist(),
                     "absolute_gap_max": maximum(fd-center["gradient"])})
    return {"analytic_gradient": center["gradient"].tolist(), "steps": rows,
            "physical_force_eta_held_fixed": True, "eta_differentiated_only_for_preparation_function": True}


def critical_symmetry_diagnostics(core, harmonic):
    """Measure cell translation and sigma1 reflection, without averaging a force."""
    grid, spectral = core["pair"].grid, core["spectral"]
    modes = grid.modes_f
    V = np.exp(2j*np.pi*grid.xi_f[:, None]*modes/grid.length)/np.sqrt(grid.nf)
    T = (V * np.exp(-2j*np.pi*modes/harmonic)) @ V.conj().T
    reflection = V[:, ::-1] @ V.conj().T
    zero = np.zeros_like(T)
    translation = np.block([[T,zero],[zero,T]])
    reflected = np.block([[zero,reflection],[reflection,zero]])
    H, C = spectral["H"], spectral["covariance"]
    return {"cell_translation_covariance_gap": maximum(translation@C@translation.conj().T-C),
            "cell_translation_hamiltonian_gap": maximum(translation@H@translation.conj().T-H),
            "dirac_reflection_covariance_gap": maximum(reflected@C@reflected.conj().T-C),
            "dirac_reflection_hamiltonian_gap": maximum(reflected@H@reflected.conj().T-H),
            "group_average_applied_to_state_or_force": False}


def critical_preflight(cpu_limit=CPU_LIMIT):
    if not 0 < cpu_limit <= CPU_LIMIT:
        raise ValueError("critical query CPU limit must be in (0,120]")
    return {"nf": BALANCE_NF, "ng": BALANCE_NF-1, "nq": 4*BALANCE_NF,
            "harmonic": CRITICAL_HARMONIC, "shape_coordinates": CRITICAL_COORDINATES,
            "initial_coefficients": critical_seed_coefficients().tolist(),
            "seed_parent_amplitude": CRITICAL_SEED_AMPLITUDE,
            "maximum_iterations": CRITICAL_MAX_ITERATIONS, "maximum_actual_trials": CRITICAL_MAX_TRIALS,
            "cpu_limit_seconds": cpu_limit, "reserve_cpu_seconds": min(5., .1*cpu_limit),
            "source_preparation": "c=eta c_base; eta=(W-magnetic)/E_base, 0<eta<=4/3; trace=3 eta",
            "solver": "seven analytic eta gradients; finite differences of these gradients for damped Newton/LM",
            "physical_source_variation": "hold eta fixed in every original common-action force",
            "default_execution": "preflight only; explicit execute-critical required"}


def critical_newton(evaluator, initial, *, maximum_iterations=CRITICAL_MAX_ITERATIONS,
                    target=NUMERICAL_TARGET):
    """Small gradient-root Newton/LM with measured admissibility backtracking."""
    core = evaluator(np.asarray(initial, dtype=float), "initial")
    history = []
    termination = "maximum seven-coordinate iterations reached"
    for iteration in range(maximum_iterations):
        gradient = core["gradient"]
        if maximum(gradient) <= target:
            termination = "seven analytic gradients meet the numerical solver target"
            break
        coefficients = core["coefficients"]
        jacobian = np.empty((gradient.size, gradient.size))
        for column in range(gradient.size):
            delta = np.zeros(gradient.size)
            step = 1e-5 * max(1., abs(coefficients[column]))
            plus = minus = None
            for _ in range(4):
                delta[column] = step
                try:
                    plus = evaluator(coefficients+delta, "gradient_jacobian_plus")
                except NonAdmissible:
                    plus = None
                try:
                    minus = evaluator(coefficients-delta, "gradient_jacobian_minus")
                except NonAdmissible:
                    minus = None
                if plus is not None or minus is not None:
                    break
                step *= .5
            if plus is not None and minus is not None:
                jacobian[:, column] = (plus["gradient"]-minus["gradient"])/(2*step)
            elif plus is not None:
                jacobian[:, column] = (plus["gradient"]-gradient)/step
            elif minus is not None:
                jacobian[:, column] = (gradient-minus["gradient"])/step
            else:
                return core, {"termination": "no admitted local gradient derivative", "iterations": history,
                              "unresolved_coordinate": column, "success": False}
        step, *_ = np.linalg.lstsq(jacobian, -gradient, rcond=None)
        norm = float(np.linalg.norm(step))
        if not np.isfinite(step).all():
            return core, {"termination": "nonfinite seven-coordinate Newton step", "iterations": history, "success": False}
        if norm > 1:
            step /= norm  # Numerical trust step, never clipping Q, eta or a^2.
        baseline = float(np.dot(gradient, gradient))
        accepted = None
        used = None
        method = "Newton"
        for phase in ("Newton", "LM"):
            if phase == "LM":
                method = phase
                damping = 1e-6 * max(1., float(np.linalg.norm(jacobian, 2))**2)
                step = np.linalg.solve(jacobian.T@jacobian+damping*np.eye(gradient.size), -jacobian.T@gradient)
                if np.linalg.norm(step) > 1:
                    step /= np.linalg.norm(step)
            for backtrack in range(12):
                fraction = 2.**(-backtrack)
                try:
                    candidate = evaluator(coefficients+fraction*step, "backtrack_"+phase)
                except NonAdmissible:
                    continue
                if float(np.dot(candidate["gradient"], candidate["gradient"])) < baseline:
                    accepted, used = candidate, fraction
                    break
            if accepted is not None:
                break
        history.append({"iteration": iteration, "gradient_max": maximum(gradient),
                        "gradient_jacobian_symmetry_defect": maximum(jacobian-jacobian.T),
                        "method": method, "accepted_fraction": used})
        if accepted is None:
            termination = "no admitted step decreases the seven-gradient residual"
            break
        core = accepted
    return core, {"termination": termination, "iterations": history,
                  "success": maximum(core["gradient"]) <= target}


def _critical_snapshot(core):
    measured, arrays = diagnostics(core["pair"], core["unknown"],
                                   (core["state"],core["spectral"],core["raw"],core["rate"],
                                    core["bundle"],core["constraints"]), hellmann_feynman=False)
    arrays.update(shape_coefficients=core["coefficients"], eta_gradient=core["gradient"],
                  radial_shape=core["R"], radius_pQ_coefficient=core["PR"],
                  radius_pQ_remainder=core["remainder"], radial_shape_tangents=core["tangents"],
                  radial_shape_mu=core["mu"])
    return measured, arrays


def run_critical_shape(*, execute=False, cpu_limit=CPU_LIMIT, producer_commit=None):
    """True local retained critical query; no NF128 execution by default."""
    plan = critical_preflight(cpu_limit)
    started = time.process_time()
    before = source_hashes()
    commit = producer_commit_binding(producer_commit, before)
    inputs = {str(coupling.LOCKED_RECORD.relative_to(ROOT)): sha256(coupling.LOCKED_RECORD)}
    report = {"schema": CRITICAL_SCHEMA, "mode": "critical_shape", "status": "CRITICAL_SHAPE_PREFLIGHT",
              "nf": BALANCE_NF, "ng": BALANCE_NF-1, "nq": 4*BALANCE_NF,
              "period": coupling.PERIOD, "gauge": "conformal", "preflight": plan,
              "base_occupations": WEIGHTS.tolist(), "occupations": None,
              "source_eta": None, "total_source_trace": None,
              "source_preparation_branch": "eta-scaled standing Gaussian; changed initial covariance selected by static action balance",
              "multiplicity": 4*coupling.KAPPA, "locked_coefficients": coupling.locked_coefficients(),
              "source_hashes": before, "input_hashes": inputs, "producing_commit": commit,
              "numpy_version": np.__version__, "scipy_version": scipy.__version__,
              "trials": [], "actual_trials": 0, "numerical_critical_query_executed": bool(execute),
              "solver": {"performed": False, "termination": "preflight only"},
              "scope": {"finite_retained_common_action_domain": True, "source_preparation_changes_between_candidates": True,
                        "eta_held_fixed_in_physical_force": True, "coefficients_changed": False,
                        "source_pump_or_reset_during_dynamics": False, "dynamics_performed": False,
                        "holding_claim": False, "stability_assessed": False, "continuum_solution_certified": False,
                        "nonexistence_inferred_from_optimizer": False}}
    if not execute:
        report.update(cpu_seconds=time.process_time()-started, cpu_budget_exceeded=False)
        return report, {}
    best = None
    initial_core = None
    cache = OrderedDict()
    largest_trial_cost = 0.
    cutoff = cpu_limit-plan["reserve_cpu_seconds"]
    def oracle(coefficients, kind):
        nonlocal best, largest_trial_cost
        key = tuple(float(value) for value in coefficients)
        if key in cache:
            cache.move_to_end(key)
            return cache[key]
        if len(report["trials"]) >= CRITICAL_MAX_TRIALS:
            raise BudgetReached("critical query reached its finite trial admission cap")
        remaining = cutoff-(time.process_time()-started)
        if remaining <= largest_trial_cost:
            raise BudgetReached("critical query reached aggregate CPU admission allowance")
        began = time.process_time()
        row = {"index": len(report["trials"]), "kind": kind, "coefficients": list(key)}
        try:
            core = critical_shape_evaluate(context, np.asarray(coefficients))
        except (NonAdmissible, np.linalg.LinAlgError, coupling.PositiveChartExit) as error:
            row.update(status="NONADMISSIBLE", reason=str(error), details=getattr(error,"details",{}))
            core = None
        cost = time.process_time()-began
        largest_trial_cost = max(largest_trial_cost,cost)
        row["cpu_seconds"] = cost
        if core is not None:
            score = float(np.dot(core["gradient"],core["gradient"]))
            row.update(status="MEASURED", source_eta=core["eta"], total_source_trace=3*core["eta"],
                       gradient=core["gradient"].tolist(), gradient_max=maximum(core["gradient"]),
                       radius_amplitude_square=core["radius_amplitude_square"],
                       raw_residual_max={name:maximum(core["raw"][name]) for name in BLOCKS},
                       reduction=core["reduction"])
            core["trial_index"] = row["index"]
            if best is None or score < float(np.dot(best["gradient"],best["gradient"])):
                best = core
            cache[key] = core
            if len(cache) > 16:
                cache.popitem(last=False)
        report["trials"].append(row)
        if core is None:
            raise NonAdmissible(row["reason"], **row["details"])
        return core
    with threadpool_limits(limits=1):
        context = critical_context()
        report["preflight"]["matrix_and_frame_cpu_seconds"] = time.process_time()-started
        initial = critical_seed_coefficients()
        try:
            initial_core = oracle(initial,"initial")
            try:
                check = critical_gradient_check(context, initial,
                    evaluator=lambda unused, coefficients: oracle(coefficients,"initial_eta_gradient_check"))
                report["eta_gradient_finite_difference_check"] = dict(check, performed=True)
            except (BudgetReached, NonAdmissible) as error:
                report["eta_gradient_finite_difference_check"] = {"performed":False,"reason":str(error)}
            _, solver = critical_newton(oracle, initial)
            report["solver"] = dict(solver, performed=True)
        except (BudgetReached, NonAdmissible) as error:
            report["solver"] = {"performed":True,"termination":str(error),"success":False}
        if best is None:
            report["status"] = "NO_ADMISSIBLE_CRITICAL_SEED"
            arrays = {}
        else:
            measured, arrays = _critical_snapshot(best)
            report.update(measurements=measured, source_eta=best["eta"], total_source_trace=3*best["eta"],
                          occupations=(best["eta"]*WEIGHTS).tolist(), selected_trial_index=best["trial_index"],
                          analytic_eta_gradient=best["gradient"].tolist(), eta_gradient_max=maximum(best["gradient"]),
                          radius_amplitude_square=best["radius_amplitude_square"], radial_scale=best["scale"],
                          reduction=best["reduction"], symmetry=critical_symmetry_diagnostics(best,CRITICAL_HARMONIC))
            # Verify affinity at an independent radius while eta stays fixed.
            third = best["state"].copy()
            third.r = 1.5*best["R"]
            third_rate, _ = galerkin.compose_fine_hamiltonian(best["pair"].grid,third)
            report["finite_radius_affinity_defect_max"] = maximum(third_rate.p_Q-(2.25*best["PR"]+best["remainder"]))
            arrays["radius_third_pQ"] = third_rate.p_Q.copy()
            reference = initial_core if initial_core is not None else best
            scales = residual_scales(reference["raw"])
            normalized = np.concatenate([best["raw"][name]/scales[name] for name in BLOCKS])
            arrays["normalized_retained_residual"] = normalized
            report.update(residual_scales=scales, normalized_retained_residual_max=maximum(normalized),
                          numerical_solver_target=NUMERICAL_TARGET,
                          initial_source_eta=reference["eta"], initial_eta_gradient=reference["gradient"].tolist())
            initial_measured, initial_arrays = _critical_snapshot(reference)
            report["initial_measurements"] = initial_measured
            arrays.update({"initial_"+name:value for name,value in initial_arrays.items()})
            satisfied = (maximum(best["gradient"]) <= NUMERICAL_TARGET and maximum(normalized) <= NUMERICAL_TARGET
                         and not measured["spectral"]["unequal_occupation_degeneracy"])
            report["status"] = "CRITICAL_RETAINED_CANDIDATE" if satisfied else "UNSATISFIED_LOCAL_CRITICAL_RELATIONS"
    after = source_hashes()
    if after != before or any(sha256(ROOT/path)!=digest for path,digest in inputs.items()):
        raise RuntimeError("critical producer bytes or locked inputs changed during measurement")
    elapsed = time.process_time()-started
    report.update(actual_trials=len(report["trials"]), cpu_seconds=elapsed, cpu_limit_seconds=cpu_limit,
                  cpu_budget_exceeded=bool(elapsed>cpu_limit), largest_measured_trial_cpu_seconds=largest_trial_cost,
                  best_candidate_preserved=best is not None)
    return report, arrays


def output_paths(path=OUTPUT):
    supplied = Path(path).expanduser()
    # scripts/lab.py changes cwd to lab; preserve explicit root-relative
    # lab/... paths used by root commands and the owning note.
    if not supplied.is_absolute() and supplied.parts and supplied.parts[0] == "lab":
        supplied = ROOT / supplied
    stem = supplied.resolve()
    if stem.suffix in (".json", ".npz"):
        stem = stem.with_suffix("")
    temporary = Path("/tmp").resolve()
    successor = (stem.parent == OUTPUT.parent.resolve()
                 and re.fullmatch(r"nsc-discovery-stationary-[A-Za-z0-9][A-Za-z0-9_-]*", stem.name))
    if not successor and temporary not in stem.parents:
        raise ValueError("stationary output must use a stationary-* successor stem or a /tmp stem")
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
        raise ValueError("stationary payload exceeds sixteen MiB")
    record = dict(report, payload=npz_path.name, payload_sha256=hashlib.sha256(payload).hexdigest(),
                  payload_bytes=len(payload))
    encoded = (json.dumps(record, indent=2, allow_nan=False) + "\n").encode()
    if len(encoded) > MAX_BYTES:
        raise ValueError("stationary JSON exceeds sixteen MiB")
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


def _check_saved_measurements(saved, measurements, *, nf, period, weights=WEIGHTS):
    """Algebraic replay shared by a stationary iterate and balance trials."""
    ng, nq = nf - 1, 4 * nf
    u, v, chi = np.split(saved("unknown"), 3)
    if (u.shape != (ng,) or maximum(saved("Q") - np.exp(u)) > 1e-12
            or maximum(saved("r") - np.exp(v)) > 1e-12 or maximum(saved("chi") - chi) > 1e-12):
        raise ValueError("stationary coordinates and saved geometry disagree")
    for name in BLOCKS:
        if maximum(saved("residual_" + name)) != measurements["raw_residual_max"][name]:
            raise ValueError("saved stationary residual summary disagrees: " + name)
    columns = np.vstack((saved("phi0"), saved("phi1")))
    if columns.shape != (2 * nf, 6) or maximum(columns.conj().T @ columns - np.eye(6)) > 1e-10:
        raise ValueError("saved stationary CAR columns disagree")
    levels = saved("eigenvalues")
    if levels.shape != (6,) or np.min(levels) <= 0:
        raise ValueError("saved positive retained eigenvalues disagree")
    if float(np.min(saved("fine_rho"))) != measurements["quadrature_rho_min"]:
        raise ValueError("saved quadrature source summary disagrees")
    if "virial" in measurements:
        fineQ, finer, finechi = (saved(name) for name in ("fine_Q", "fine_r", "fine_chi"))
        if any(values.shape != (nq,) for values in (fineQ, finer, finechi)):
            raise ValueError("saved fine stationary domain disagrees")
        for name in BLOCKS[:3]:
            if maximum(saved("fine_" + name)) != measurements["fine_momentum_residual_max"][name]:
                raise ValueError("saved fine momentum residual summary disagrees")
        coefficients = coupling.locked_coefficients()
        auxiliary = float(coupling.alpha_of(coefficients["C_W"]) * period/nq * np.sum(fineQ**2 * finechi**2))
        source = float(period/nq * np.dot(fineQ, saved("fine_rho")))
        magnetic = float(2 * np.pi * coefficients["C_F"] * coefficients["flux"]**2 * period/nq * np.sum(fineQ**2))
        fine_rhs = float(-.5 * period/nq * np.dot(finer, saved("fine_p_r"))
                         - period/nq * np.dot(finechi, saved("fine_p_chi"))
                         - period/nq * np.dot(fineQ, saved("fine_lapse")))
        coarse_rhs = float(-.5 * period/ng * np.dot(saved("r"), saved("residual_p_r"))
                           - period/ng * np.dot(saved("chi"), saved("residual_p_chi"))
                           - period/ng * np.dot(saved("Q"), saved("residual_lapse")))
        virial = measurements["virial"]
        for name, value in (("auxiliary_term", auxiliary), ("source_term", source), ("magnetic_term", magnetic),
                            ("gap", auxiliary-source-magnetic), ("fine_rhs", fine_rhs), ("coarse_rhs", coarse_rhs)):
            if abs(value - virial[name]) > 1e-11 * max(1., abs(value)):
                raise ValueError("saved virial summary disagrees: " + name)
        energy = float(4 * coupling.KAPPA * np.dot(weights, levels))
        if abs(energy - measurements["field_eigen_energy"]) > 1e-12 * max(1., abs(energy)):
            raise ValueError("saved spectral source energy summary disagrees")


def check_record(path=OUTPUT):
    """Read-only hashes and algebra; sealed producer bytes come from their commit.

    A historical commit pin authenticates the recorded producer set rather
    than requiring today's extended source/docs/tests to equal those bytes.
    Unpinned temporary records still require exact current source hashes.
    """
    json_path, npz_path = output_paths(path)
    record = json.loads(json_path.read_text())
    balance = record.get("mode") == "balance_seed"
    critical = record.get("mode") == "critical_shape"
    nf = BALANCE_NF if balance or critical else NF
    if critical:
        eta = record.get("source_eta")
        occupations_valid = record.get("base_occupations") == WEIGHTS.tolist()
        if eta is None:
            occupations_valid = occupations_valid and record.get("occupations") is None
        else:
            occupations_valid = (occupations_valid and 0 < eta <= 4/3
                                 and record.get("occupations") == (eta*WEIGHTS).tolist()
                                 and abs(record.get("total_source_trace",0)-3*eta) <= 1e-12)
    else:
        occupations_valid = record.get("occupations") == WEIGHTS.tolist()
    if (record.get("schema") not in (LEGACY_SCHEMA, SCHEMA, CRITICAL_SCHEMA) or record.get("nf") != nf
            or record.get("ng") != nf - 1 or record.get("nq") != 4 * nf
            or record.get("period") != coupling.PERIOD or record.get("gauge") != "conformal"
            or not occupations_valid
            or record.get("locked_coefficients") != coupling.locked_coefficients()):
        raise ValueError("stationary record domain does not match")
    if set(record.get("source_hashes", {})) != set(source_hashes()):
        raise ValueError("stationary producer set changed")
    commit = record.get("producing_commit")
    if commit is None:
        if record["source_hashes"] != source_hashes():
            raise ValueError("unpinned stationary producer hashes changed")
        authentication = "current_producer_bytes"
    else:
        producer_commit_binding(commit, record["source_hashes"])
        authentication = "immutable_producing_commit"
    input_path = str(coupling.LOCKED_RECORD.relative_to(ROOT))
    if set(record.get("input_hashes", {})) != {input_path}:
        raise ValueError("locked stationary input set changed")
    if commit is None:
        expected_digest = sha256(coupling.LOCKED_RECORD)
    else:
        data = subprocess.check_output(["git", "show", commit + ":" + input_path], cwd=ROOT)
        expected_digest = hashlib.sha256(data).hexdigest()
    if record["input_hashes"][input_path] != expected_digest:
        raise ValueError("locked stationary input hash changed")
    if "parent_record" in record:
        parent = record["parent_record"]
        parent_path = Path(parent["json"])
        if not parent_path.is_absolute():
            parent_path = ROOT / parent_path
        if sha256(parent_path) != parent["json_sha256"]:
            raise ValueError("sealed parent diagnostic record hash changed")
        parent_record = json.loads(parent_path.read_text())
        if sha256(parent_path.parent / parent_record["payload"]) != parent["payload_sha256"]:
            raise ValueError("sealed parent diagnostic payload hash changed")
    if record.get("payload") != npz_path.name or sha256(npz_path) != record.get("payload_sha256"):
        raise ValueError("stationary payload hash changed")
    if npz_path.stat().st_size != record.get("payload_bytes") or npz_path.stat().st_size > MAX_BYTES:
        raise ValueError("stationary payload size changed")
    with np.load(npz_path, allow_pickle=False) as arrays:
        for name in arrays.files:
            if not np.isfinite(arrays[name]).all():
                raise ValueError("nonfinite stationary array: " + name)
        if critical:
            if record.get("source_eta") is not None:
                weights = np.asarray(record["occupations"])
                saved = lambda name: arrays[name]
                _check_saved_measurements(saved,record["measurements"],nf=nf,period=record["period"],weights=weights)
                columns = np.vstack((arrays["phi0"],arrays["phi1"]))
                covariance = (columns*weights)@columns.conj().T
                if maximum(covariance-arrays["covariance"]) > 1e-11:
                    raise ValueError("critical eta covariance disagrees")
                if maximum(arrays["H"]@columns-columns*arrays["eigenvalues"]) > 1e-8:
                    raise ValueError("critical retained eigenmodes disagree")
                if maximum(arrays["eta_gradient"]-record["analytic_eta_gradient"]) > 1e-12:
                    raise ValueError("critical gradient summary disagrees")
                Ebase = record["reduction"]["base_source_energy"]
                gradient = record["period"]/(nf-1) * arrays["radial_shape_tangents"].T@arrays["radius_pQ_remainder"]/Ebase
                if maximum(gradient-arrays["eta_gradient"]) > 1e-11:
                    raise ValueError("critical analytic gradient and fixed-eta force disagree")
                if maximum(arrays["r"]-np.sqrt(record["radius_amplitude_square"])*arrays["radial_shape"]) > 1e-11:
                    raise ValueError("critical normal radius scale disagrees")
                if maximum(arrays["radius_third_pQ"]-(2.25*arrays["radius_pQ_coefficient"]+arrays["radius_pQ_remainder"])) > 1e-8:
                    raise ValueError("critical finite radius affinity disagrees")
                normalized = np.concatenate([arrays["residual_"+name]/record["residual_scales"][name] for name in BLOCKS])
                if maximum(normalized-arrays["normalized_retained_residual"]) > 1e-12:
                    raise ValueError("critical full retained residual blocks disagree")
                if maximum(normalized) != record["normalized_retained_residual_max"]:
                    raise ValueError("critical retained residual summary disagrees")
                initial_saved = lambda name: arrays["initial_"+name]
                _check_saved_measurements(initial_saved,record["initial_measurements"],nf=nf,period=record["period"],
                                          weights=record["initial_source_eta"]*WEIGHTS)
                selected = record["trials"][record["selected_trial_index"]]
                if maximum(arrays["shape_coefficients"]-selected["coefficients"]) > 1e-12 or selected["source_eta"] != eta:
                    raise ValueError("critical selected trial disagrees")
            elif len(arrays.files) != 0:
                raise ValueError("critical preflight/nonadmitted seed has unexpected payload")
        elif balance:
            for row in record["trials"]:
                if row["status"] != "MEASURED":
                    continue
                prefix = row["array_prefix"]
                saved = lambda name: arrays[prefix + name]
                _check_saved_measurements(saved, row["measurements"], nf=nf, period=record["period"])
                if row["virial_gap"] != row["measurements"]["virial"]["gap"]:
                    raise ValueError("balance trial gap disagrees")
            chosen = record["selected_trial_index"]
            if chosen is not None and record["measurements"] != record["trials"][chosen]["measurements"]:
                raise ValueError("selected balance trial measurements disagree")
            radius = record.get("radius_selection", {})
            if radius.get("performed"):
                coefficient, remainder = arrays["radius_pQ_coefficient"], arrays["radius_pQ_remainder"]
                denominator = float(np.dot(coefficient, coefficient))
                fitted = None if denominator == 0 else float(-np.dot(coefficient, remainder) / denominator)
                if fitted != radius["least_squares_amplitude_square"]:
                    raise ValueError("saved conditional radius fit disagrees")
                if maximum(arrays["radius_third_pQ"] - remainder - 2.25*coefficient) > 1e-8:
                    raise ValueError("saved radius affinity verification disagrees")
                if fitted is not None and fitted > 0:
                    prefix = record["trials"][radius["trial_index"]]["array_prefix"]
                    if maximum(arrays["radius_selected_r"] - radius["radius_amplitude"] * arrays[prefix+"r"]) > 1e-11:
                        raise ValueError("saved conditional radius geometry disagrees")
                    for name in BLOCKS:
                        if maximum(arrays["radius_selected_residual_"+name]) != radius["raw_residual_max"][name]:
                            raise ValueError("saved conditional radius residual summary disagrees")
                    for name in BLOCKS[:3]:
                        if maximum(arrays["radius_selected_fine_"+name]) != radius["fine_momentum_residual_max"][name]:
                            raise ValueError("saved conditional fine radius residual summary disagrees")
            for root in record["roots"]:
                if root["status"] == "VIRIAL_SELECTED_SHAPE":
                    row = record["trials"][root["trial_index"]]
                    if root["amplitude"] != row["amplitude"] or root["virial_gap"] != row["virial_gap"]:
                        raise ValueError("selected virial root disagrees with its measured trial")
        else:
            saved = lambda name: arrays[name]
            _check_saved_measurements(saved, record["measurements"], nf=nf, period=record["period"])
            u = arrays["unknown"][:nf-1]
            normalized = np.concatenate([arrays["residual_" + name] / record["residual_scales"][name]
                                         for name in BLOCKS] + [[np.sqrt(nf-1) * float(np.mean(
                                             u * np.sin(2 * np.pi * np.arange(nf-1) / (nf-1))))]])
            if maximum(arrays["normalized_residual"] - normalized) > 1e-12:
                raise ValueError("saved normalized residual and raw blocks disagree")
            if maximum(arrays["normalized_residual"]) != record["normalized_residual_max"]:
                raise ValueError("normalized stationary residual summary disagrees")
            columns = np.vstack((arrays["phi0"], arrays["phi1"]))
            C = (columns * WEIGHTS) @ columns.conj().T
            if maximum(C - arrays["covariance"]) > 1e-12:
                raise ValueError("saved stationary covariance disagrees")
            if maximum(arrays["H"] @ columns - columns * arrays["eigenvalues"]) > 1e-8:
                raise ValueError("saved positive retained eigenmodes disagree")
    return {"ok": True, "schema": record["schema"], "status": record["status"],
            "source_authentication": authentication, "bytes_written": 0, "solver_executed": False}
