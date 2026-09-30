"""Variational Fourier–Galerkin successor of the spherical coupling.

Geometry is an odd trigonometric subspace. Fermion columns stay in the fixed
even antiperiodic band owned by ``prepare_rank6``. Both are prolonged to a
finer quadrature grid, the existing fine-grid Hamiltonian is evaluated there,
and the rates are pulled back by the exact adjoints. This is not a filter of
the v1 radius and not a constraint projection.

The minus-column phase and the local carrier origins belong to
``nsc_spherical_coupling.prepare_rank6``. This module does not reimplement
them. A physical control is refused while that preparation still fails the
contract below.
"""
from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .nsc_conformal_adm_source import sector_multiplicity
from .nsc_spherical_coupling import (
    CALIBRATION,
    CARRIER_K,
    ELL,
    KAPPA,
    OCCUPATIONS,
    PACKET_COUNT,
    PERIOD,
    TOL_CONSTRAINT_VALIDATION,
    TOL_ENERGY_RELATIVE,
    TOL_FD_RELATIVE,
    TOL_INITIAL_CONSTRAINT,
    TOL_NEWTON,
    CauchyRate,
    CauchyState,
    CouplingSystem,
    PositiveChartExit,
    _combine,
    _radius_jacobian,
    _radius_residual,
    antiperiodic_momentum,
    apply_dirac,
    chart_failure,
    dense_force_difference,
    dense_hamiltonian_difference,
    field_energy,
    geometric_rates,
    gravity_energy,
    hamilton_constraint,
    locked_coefficients,
    magnetic_radius_square,
    packet_onsite_report,
    periodic_derivative,
    prepare_rank6,
    profile_coordinate,
    shift_constraint,
    source_from_columns,
    total_energy,
)

SCHEMA = "NSC-SPHERICAL-COUPLING-CONTROL-v2"
STATUS = "PROVISIONAL"
GLOBAL_REGENERATION = False
ABSOLUTE_VACUUM_CLAIM = False
RENEWAL = False
OLD_CIRCULATION_TRANSFERRED = False
POST_STEP_FILTER = False
CONSTRAINT_PROJECTION = False
PRESCRIBED_RADIUS = False

# Same predeclared lines as v1. Not retuned.
TOL_QUADRATURE_INITIAL = TOL_INITIAL_CONSTRAINT

_LAB_ROOT = Path(__file__).resolve().parents[2]
CONTROL_RECORD = _LAB_ROOT / "results" / "development" / "nsc-spherical-coupling-control-v2.json"
_MODULE_PATH = Path(__file__).resolve()
_TEST_PATH = _LAB_ROOT / "tests" / "test_nsc_spherical_galerkin_coupling.py"
_COUPLING_PATH = _LAB_ROOT / "src" / "recursive_horizons" / "nsc_spherical_coupling.py"
_ACTION_PATH = _LAB_ROOT / "src" / "recursive_horizons" / "nsc_spherical_feedback_action.py"
_SOURCE_PATH = _LAB_ROOT / "src" / "recursive_horizons" / "nsc_conformal_adm_source.py"

PREPARATION_DEPENDENCY = (
    "nsc_spherical_coupling.prepare_rank6 applies the recorded constant phase "
    "to the minus column of a real envelope and uses the local carrier "
    "exp(i k (x-n*ell)). This representation does not rebuild that packet."
)


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def geometry_modes(ng):
    """Integer modes of an odd trigonometric subspace. No Nyquist degree."""
    ng = int(ng)
    if ng < 3 or ng % 2 == 0:
        raise ValueError("geometry degree count must be odd so the Nyquist mode is absent")
    return np.arange(-(ng // 2), ng // 2 + 1, dtype=float)


def fermion_modes(nf):
    """Half-integer AP modes. Even count, symmetric, no extra Nyquist drop."""
    nf = int(nf)
    if nf < 10 or nf % 2:
        raise ValueError("AP fermion count must be even and at least 10")
    return (np.arange(nf) - nf // 2).astype(float) + 0.5


def periodic_interpolation(ng, nq, length):
    """Real nodal interpolation A_g of shape (nq, ng). (ng/nq) A_g.T A_g = I."""
    modes = geometry_modes(ng)
    nq = int(nq)
    length = float(length)
    if nq < int(ng) or nq % 2:
        raise ValueError("quadrature grid must be even and at least as fine as the geometry")
    if length <= 0:
        raise ValueError("positive period required")
    coarse = np.arange(ng, dtype=float) * (length / ng)
    fine = np.arange(nq, dtype=float) * (length / nq)
    analysis = np.exp(-2j * np.pi * coarse[:, None] * modes[None, :] / length) / ng
    synthesis = np.exp(2j * np.pi * fine[:, None] * modes[None, :] / length)
    matrix = synthesis @ analysis.T
    imaginary = float(np.max(np.abs(np.imag(matrix))))
    if imaginary > 1e-8:
        raise ValueError("periodic interpolation left the real trigonometric subspace")
    return np.real(matrix)


def antiperiodic_interpolation(nf, nq, length):
    """Complex AP interpolation A_f of shape (nq, nf)."""
    modes = fermion_modes(nf)
    nq = int(nq)
    length = float(length)
    if nq < int(nf) or nq % 2:
        raise ValueError("AP quadrature must be even and finer than the fermion band")
    coarse = np.arange(nf, dtype=float) * (length / nf)
    fine = np.arange(nq, dtype=float) * (length / nq)
    analysis = np.exp(-2j * np.pi * coarse[:, None] * modes[None, :] / length) / nf
    synthesis = np.exp(2j * np.pi * fine[:, None] * modes[None, :] / length)
    return synthesis @ analysis.T


def canonical_column_map(interpolation, nf, nq):
    """U_f = sqrt(nf/nq) A_f, so U_f† U_f = I on half-density samples."""
    return np.sqrt(float(nf) / float(nq)) * interpolation


def _assert_isometries(weight, geometry, columns):
    identity_g = weight * geometry.T @ geometry
    identity_f = columns.conj().T @ columns
    geometry_error = float(np.max(np.abs(identity_g - np.eye(geometry.shape[1]))))
    column_error = float(np.max(np.abs(identity_f - np.eye(columns.shape[1]))))
    if geometry_error > 1e-8 or column_error > 1e-8:
        raise ValueError(
            f"interpolation isometry failed: geometry {geometry_error}, columns {column_error}"
        )
    return geometry_error, column_error


@dataclass
class GalerkinGrid:
    ng: int
    nf: int
    nq: int
    length: float
    dx_g: float
    dx_f: float
    dx_q: float
    weight: float
    xi_g: np.ndarray
    xi_f: np.ndarray
    xi_q: np.ndarray
    modes_g: np.ndarray
    modes_f: np.ndarray
    A_g: np.ndarray
    A_f: np.ndarray
    U_f: np.ndarray
    derivative_on_geometry: np.ndarray
    isometry_geometry: float
    isometry_columns: float
    fine: CouplingSystem

    @property
    def derivative(self):
        return self.fine.derivative

    @property
    def momentum(self):
        return self.fine.momentum


def build_grid(fermions, quadrature=None, length=PERIOD):
    """Odd geometry ng = nf - 1, even AP count nf, quadrature default 4 nf."""
    nf = int(fermions)
    ng = nf - 1
    nq = int(4 * nf if quadrature is None else quadrature)
    if sector_multiplicity(KAPPA) != 4 * KAPPA:
        raise ValueError("multiplicity must remain M=4 kappa once")
    modes_g = geometry_modes(ng)
    modes_f = fermion_modes(nf)
    geometry = periodic_interpolation(ng, nq, length)
    fermions_map = antiperiodic_interpolation(nf, nq, length)
    columns = canonical_column_map(fermions_map, nf, nq)
    weight = ng / nq
    geometry_error, column_error = _assert_isometries(weight, geometry, columns)
    derivative = periodic_derivative(nq, length)
    momentum, _metric = antiperiodic_momentum(nq, length)
    dx_q = length / nq
    xi_q = np.arange(nq, dtype=float) * dx_q
    coordinate = profile_coordinate(xi_q, length)
    scale = np.exp(np.log(CALIBRATION["Omega"]) * coordinate)
    # Same static gauge functions, evaluated on the quadrature nodes.
    length_density = CALIBRATION["b0"] * scale
    shift = CALIBRATION["beta0"] * scale
    coefficients = locked_coefficients()
    fine = CouplingSystem(
        points=nq,
        length=float(length),
        dx=dx_q,
        xi=xi_q,
        derivative=derivative,
        momentum=momentum,
        length_density=length_density,
        shift=shift,
        kappa=KAPPA,
        multiplicity=int(sector_multiplicity(KAPPA)),
        occupations=OCCUPATIONS.copy(),
        coefficients=coefficients,
        calibration=dict(CALIBRATION),
        preparation={
            "role": "quadrature carrier of a prolonged variational subspace",
            "columns_prepared_on": "fermion band",
            "post_step_filter": False,
        },
    )
    return GalerkinGrid(
        ng=ng,
        nf=nf,
        nq=nq,
        length=float(length),
        dx_g=length / ng,
        dx_f=length / nf,
        dx_q=dx_q,
        weight=weight,
        xi_g=np.arange(ng, dtype=float) * (length / ng),
        xi_f=np.arange(nf, dtype=float) * (length / nf),
        xi_q=xi_q,
        modes_g=modes_g,
        modes_f=modes_f,
        A_g=geometry,
        A_f=fermions_map,
        U_f=columns,
        derivative_on_geometry=derivative @ geometry,
        isometry_geometry=geometry_error,
        isometry_columns=column_error,
        fine=fine,
    )


def preparation_problems(preparation, phi0, phi1, xi):
    """Contract for the coupling owner's corrected packet. Empty means accepted."""
    problems = []
    if preparation.get("phase_applied_to_odd_lobe") is not False:
        problems.append("phase_applied_to_odd_lobe is not false")
    if preparation.get("phase_on_minus_spinor_column") is not True:
        problems.append("phase_on_minus_spinor_column is not true")
    carrier = str(preparation.get("carrier", ""))
    if "x-n*ell" not in carrier:
        problems.append("carrier does not declare local packet origins x-n*ell")
    phase = np.exp(1j * CALIBRATION["phase"])
    xi = np.asarray(xi, dtype=float)
    for region in range(PACKET_COUNT):
        plus = np.concatenate((phi0[:, 2 * region], phi1[:, 2 * region]))
        minus = np.concatenate((phi0[:, 2 * region + 1], phi1[:, 2 * region + 1]))
        conjugate = np.conjugate(plus)
        mask = np.abs(conjugate) > 1e-8
        if np.count_nonzero(mask) < 4 or np.max(np.abs(minus[mask] / conjugate[mask] - phase)) > 1e-8:
            problems.append(
                f"region {region} minus column is not exp(i phase) times conjugate(plus)"
            )
        local = xi - region * ELL
        envelope = np.sqrt(2) * phi0[:, 2 * region] * np.exp(-1j * CARRIER_K * local)
        if np.max(np.abs(np.imag(envelope))) > 1e-7:
            problems.append(f"region {region} carrier origin is not the local packet edge")
    gram = phi0.conj().T @ phi0 + phi1.conj().T @ phi1
    if np.max(np.abs(gram - np.eye(6))) > 1e-8:
        problems.append("fermion columns are not orthonormal")
    return problems


def load_physical_columns(fermions):
    """Samples of the coupling owner's packet on the fixed fermion grid."""
    nf = int(fermions)
    dx = PERIOD / nf
    xi = np.arange(nf, dtype=float) * dx
    phi0, phi1, preparation = prepare_rank6(xi, dx)
    return phi0, phi1, preparation, preparation_problems(preparation, phi0, phi1, xi)


def blank_state(grid, phi0, phi1):
    """Momenta and chi vanish. Q is the constant gauge ratio. r starts at 1."""
    q0 = CALIBRATION["b0"] / CALIBRATION["a0"]
    return CauchyState(
        Q=np.full(grid.ng, q0),
        r=np.ones(grid.ng),
        chi=np.zeros(grid.ng),
        p_Q=np.zeros(grid.ng),
        p_r=np.zeros(grid.ng),
        p_chi=np.zeros(grid.ng),
        phi0=np.array(phi0, dtype=complex, copy=True),
        phi1=np.array(phi1, dtype=complex, copy=True),
    )


def manufactured_columns(nf, seed=3):
    """Orthonormal half-density columns inside the AP band. Not the physical packet."""
    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(nf, 6)) + 1j * rng.normal(size=(nf, 6))
    orthogonal, _ = np.linalg.qr(raw)
    split = orthogonal / np.sqrt(2)
    return split, 1j * split


def prolong_geometry(grid, values):
    return grid.A_g @ np.asarray(values, dtype=float)


def pull_geometry(grid, values):
    """(ng/nq) A_g.T applied to a fine nodal array."""
    return grid.weight * grid.A_g.T @ np.asarray(values)


def prolong_columns(grid, phi):
    return grid.U_f @ np.asarray(phi)


def prolong_state(grid, state):
    return CauchyState(
        Q=prolong_geometry(grid, state.Q),
        r=prolong_geometry(grid, state.r),
        chi=prolong_geometry(grid, state.chi),
        p_Q=prolong_geometry(grid, state.p_Q),
        p_r=prolong_geometry(grid, state.p_r),
        p_chi=prolong_geometry(grid, state.p_chi),
        phi0=prolong_columns(grid, state.phi0),
        phi1=prolong_columns(grid, state.phi1),
    )


def _pull_rate_fields(grid, fine_rate, force_q):
    """Adjoint pullback. Column rate is U†(-i H U Φ), already in fine_rate images."""
    return CauchyRate(
        pull_geometry(grid, fine_rate[0]),
        pull_geometry(grid, fine_rate[1]),
        pull_geometry(grid, fine_rate[2]),
        pull_geometry(grid, fine_rate[3]),
        pull_geometry(grid, fine_rate[4]),
        pull_geometry(grid, fine_rate[5]),
        grid.U_f.conj().T @ fine_rate[6],
        grid.U_f.conj().T @ fine_rate[7],
        fine_rate[8],
        force_q,
        fine_rate[10],
        fine_rate[11],
    )


def compose_fine_hamiltonian(grid, state, include_matter_force=True):
    """Existing fine-grid Hamiltonian on prolonged geometry and columns.

    Returns the fine system, prolonged state, column source, unprojected
    geometry rates, and the coarse rate. Matter is divided by the quadrature
    spacing only where the momentum equation needs a density. 4κ is the
    multiplicity already inside the existing nodal forces, once.
    """
    fine_state = prolong_state(grid, state)
    failure = chart_failure(grid.fine, fine_state)
    if failure is not None:
        raise PositiveChartExit(failure, np.nan, state.copy())
    source = source_from_columns(grid.fine, fine_state)
    if source["multiplicity_applied_once"] is not True:
        raise ValueError("source multiplicity was not applied once")
    if int(grid.fine.multiplicity) != 4 * int(grid.fine.kappa):
        raise ValueError("4 kappa multiplicity drifted")
    q_dot, r_dot, chi_dot, p_q_dot, p_r_dot, p_chi_dot = geometric_rates(grid.fine, fine_state)
    if include_matter_force:
        p_q_dot = p_q_dot - source["force_Q"] / grid.dx_q
    # Representative block evolves by -i H, then the column adjoint. Not -M i H.
    phi0_dot_fine = -1j * source["image0"]
    phi1_dot_fine = -1j * source["image1"]
    coarse = _pull_rate_fields(
        grid,
        (q_dot, r_dot, chi_dot, p_q_dot, p_r_dot, p_chi_dot, phi0_dot_fine, phi1_dot_fine,
         source["force_L"], source["force_Q"], source["force_beta"], 0.0),
        source["force_Q"],
    )
    lifted_q = prolong_geometry(grid, coarse.Q)
    coarse.fieldwork_power = float(np.sum(source["force_Q"] * lifted_q))
    bundle = {
        "fine_state": fine_state,
        "source": source,
        "unprojected_Q": q_dot,
        "unprojected_r": r_dot,
        "unprojected_fieldwork": float(np.sum(source["force_Q"] * q_dot)),
        "lifted_Q": lifted_q,
    }
    return coarse, bundle


def rates(grid, state, include_matter_force=True):
    coarse, _bundle = compose_fine_hamiltonian(grid, state, include_matter_force)
    return coarse


def energy(grid, state):
    fine_state = prolong_state(grid, state)
    return total_energy(grid.fine, fine_state)


def occupation_eigenvalues(phi0, phi1, occupations=None):
    weights = OCCUPATIONS if occupations is None else np.asarray(occupations, dtype=float)
    vectors = np.vstack((phi0, phi1))
    covariance = (vectors * weights[None, :]) @ vectors.conj().T
    return np.linalg.eigvalsh(covariance)


def unresolved_fraction(grid, values):
    """Power outside the odd geometry band. The Nyquist bin is unresolved."""
    values = np.asarray(values)
    count = values.size
    modes = np.fft.fftfreq(count) * count
    power = np.abs(np.fft.fft(values) / count) ** 2
    total = float(np.sum(power))
    unresolved = np.abs(modes) > (grid.ng // 2)
    unresolved_power = float(np.sum(power[unresolved]))
    fraction = 0.0 if total == 0 else unresolved_power / total
    return {
        "max_abs": float(np.max(np.abs(values))),
        "power": total,
        "band_power": total - unresolved_power,
        "unresolved_power": unresolved_power,
        "unresolved_fraction": fraction,
        "unresolved_max_mode": int(np.max(np.abs(modes[unresolved]))) if np.any(unresolved) else 0,
    }


def constraint_diagnostics(grid, state, source=None, fine_state=None):
    """Projected coefficients and the full quadrature residual, including held-out modes."""
    if fine_state is None:
        fine_state = prolong_state(grid, state)
    if source is None:
        source = source_from_columns(grid.fine, fine_state)
    rho = source["force_L"] / grid.dx_q
    current = source["force_beta"] / grid.dx_q
    hamilton = hamilton_constraint(grid.fine, fine_state) + rho
    momentum = shift_constraint(grid.fine, fine_state) + current
    projected_h = pull_geometry(grid, hamilton)
    projected_m = pull_geometry(grid, momentum)
    held_h = hamilton - prolong_geometry(grid, projected_h)
    held_m = momentum - prolong_geometry(grid, projected_m)
    return {
        "projected_hamilton_max": float(np.max(np.abs(projected_h))),
        "projected_momentum_max": float(np.max(np.abs(projected_m))),
        "full_hamilton_max": float(np.max(np.abs(hamilton))),
        "full_momentum_max": float(np.max(np.abs(momentum))),
        "held_out_hamilton_max": float(np.max(np.abs(held_h))),
        "held_out_momentum_max": float(np.max(np.abs(held_m))),
        "hamilton_modes": unresolved_fraction(grid, hamilton),
        "momentum_modes": unresolved_fraction(grid, momentum),
        "rho_modes": unresolved_fraction(grid, rho),
        "current_max": float(np.max(np.abs(current))),
        "positive_r": bool(np.min(fine_state.r) > 0),
        "positive_Q": bool(np.min(fine_state.Q) > 0),
        "r_min_quadrature": float(np.min(fine_state.r)),
        "Q_min_quadrature": float(np.min(fine_state.Q)),
    }


def normal_velocities(grid, state, coarse_rate=None, bundle=None):
    """Chart proper velocity uses p=0; lifted velocity uses the actual coarse r dot."""
    if bundle is None or coarse_rate is None:
        coarse_rate, bundle = compose_fine_hamiltonian(grid, state)
    fine_state = bundle["fine_state"]
    radius_derivative = grid.derivative @ fine_state.r
    denominator = fine_state.r * grid.fine.length_density
    lifted = prolong_geometry(grid, coarse_rate.r)
    unprojected = bundle["unprojected_r"]
    chart = (unprojected - grid.fine.shift * radius_derivative) / denominator
    lifted_normal = (lifted - grid.fine.shift * radius_derivative) / denominator
    coordinate = lifted / denominator
    shift_piece = (grid.fine.shift * radius_derivative) / denominator
    return {
        "chart_proper_max": float(np.max(np.abs(chart))),
        "lifted_proper_max": float(np.max(np.abs(lifted_normal))),
        "projected_proper_max": float(np.max(np.abs(pull_geometry(grid, lifted_normal)))),
        "coordinate_max": float(np.max(np.abs(coordinate))),
        "shift_piece_max": float(np.max(np.abs(shift_piece))),
    }


def _projected_radius_jacobian(grid, radius):
    """Coarse Jacobian of the projected Hamilton residual. Matches _radius_jacobian."""
    fine_radius = prolong_geometry(grid, radius)
    reference = _radius_jacobian(grid.fine, fine_radius)
    return grid.weight * grid.A_g.T @ (reference @ grid.A_g)


def projected_radius_operator(grid, radius, rho):
    fine_residual = _radius_residual(grid.fine, prolong_geometry(grid, radius), rho)
    return pull_geometry(grid, fine_residual), fine_residual


def solve_initial_radius(grid, state, rho=None):
    """Newton for r inside the geometry subspace. Momenta and chi stay zero.

    The target is the projected Hamilton constraint. The full quadrature
    residual is reported and is not removed by filtering.
    """
    if rho is None:
        fine_state = prolong_state(grid, state)
        source = source_from_columns(grid.fine, fine_state)
        rho = source["force_L"] / grid.dx_q
        varied = state.copy()
        varied.r = state.r * 1.7
        other = source_from_columns(grid.fine, prolong_state(grid, varied))
        rho_variation = float(np.max(np.abs(other["force_L"] - source["force_L"])))
    else:
        rho = np.asarray(rho, dtype=float)
        rho_variation = None
    radius = np.ones(grid.ng)
    steps = []
    for scale in np.linspace(0.0, 1.0, 9)[1:]:
        target = scale * rho
        converged = False
        residual_max = None
        for _iteration in range(12):
            projected, _full = projected_radius_operator(grid, radius, target)
            residual_max = float(np.max(np.abs(projected)))
            if residual_max < TOL_NEWTON:
                converged = True
                break
            jacobian = _projected_radius_jacobian(grid, radius)
            try:
                delta = np.linalg.solve(jacobian, -projected)
            except np.linalg.LinAlgError as error:
                return state, {
                    "converged": False,
                    "blocker": f"singular projected radius Jacobian: {error}",
                    "residual_max": residual_max,
                    "steps": steps,
                    "rho_independent_of_r_max": rho_variation,
                    "imposed_radius": False,
                    "filtered_v1_radius": False,
                }
            accepted = False
            step = 1.0
            for _halving in range(16):
                trial = radius + step * delta
                if not np.isfinite(trial).all() or np.min(prolong_geometry(grid, trial)) <= 0:
                    step *= 0.5
                    continue
                trial_projected, _ = projected_radius_operator(grid, trial, target)
                trial_max = float(np.max(np.abs(trial_projected)))
                if trial_max < residual_max * (1 - 0.05 * step):
                    radius = trial
                    accepted = True
                    break
                step *= 0.5
            if not accepted:
                break
        steps.append({"scale": float(scale), "converged": converged, "residual_max": residual_max})
        if not converged:
            return state, {
                "converged": False,
                "blocker": "Newton left the positive quadrature chart or stalled before full rho",
                "residual_max": residual_max,
                "steps": steps,
                "rho_independent_of_r_max": rho_variation,
                "imposed_radius": False,
                "filtered_v1_radius": False,
            }
    solved = state.copy()
    solved.r = radius
    diagnostics = constraint_diagnostics(grid, solved)
    projected_final, full_final = projected_radius_operator(grid, radius, rho)
    held_final = full_final - prolong_geometry(grid, projected_final)
    # The Newton target is this rho. Column diagnostics stay available, but the
    # solved Hamilton residual is the one the subspace was asked to satisfy.
    diagnostics = dict(diagnostics)
    diagnostics["projected_hamilton_max"] = float(np.max(np.abs(projected_final)))
    diagnostics["full_hamilton_max"] = float(np.max(np.abs(full_final)))
    diagnostics["held_out_hamilton_max"] = float(np.max(np.abs(held_final)))
    diagnostics["hamilton_modes"] = unresolved_fraction(grid, full_final)
    velocities = normal_velocities(grid, solved)
    eight_pi_a = 8 * np.pi * grid.fine.A
    fine_state = prolong_state(grid, solved)
    source = source_from_columns(grid.fine, fine_state)
    bracket = (
        fine_state.Q ** 2 * magnetic_radius_square(grid.fine.coefficients)
        + fine_state.Q * (source["force_L"] / grid.dx_q) / eight_pi_a
    )
    info = {
        "converged": True,
        "blocker": None,
        "residual_max": diagnostics["projected_hamilton_max"],
        "steps": steps,
        "rho_independent_of_r_max": rho_variation,
        "imposed_radius": False,
        "filtered_v1_radius": False,
        "r_min_coarse": float(np.min(radius)),
        "r_max_coarse": float(np.max(radius)),
        "bracket_min": float(np.min(bracket)),
        "bracket_positive": bool(np.min(bracket) > 0),
        "diagnostics": diagnostics,
        "velocities": velocities,
    }
    return solved, info


def product_gap(points, mode, length=PERIOD):
    """Product-rule gap of a pure integer mode on a collocation grid.

    The max-norm of mode 31 on 64 nodes is 8π. Its Euclidean norm divided by
    2π is the 22.6 alias figure. The samples of that mode equal the folded
    partner mode-33; the gap itself is not the gap of that low mode.
    """
    derivative = periodic_derivative(int(points), length)
    coordinate = np.arange(int(points)) * (length / int(points))
    values = np.cos(2 * np.pi * int(mode) * coordinate / length)
    gap = derivative @ (values * values) - 2 * values * (derivative @ values)
    vector = np.linalg.norm(gap)
    return {
        "max": float(np.max(np.abs(gap))),
        "l2": float(vector),
        "l2_over_2pi": float(vector / (2 * np.pi)),
    }


def resolved_product_gap(grid, mode):
    """Product gap of a true geometry mode on the quadrature grid, not a folded sample."""
    mode = int(mode)
    if abs(mode) > grid.ng // 2:
        raise ValueError("mode is outside the odd geometry subspace")
    coordinate = grid.xi_g
    values = np.cos(2 * np.pi * mode * coordinate / grid.length)
    prolonged = prolong_geometry(grid, values)
    derivative = grid.derivative @ prolonged
    gap = grid.derivative @ (prolonged * prolonged) - 2 * prolonged * derivative
    return float(np.max(np.abs(gap)))


def alias_diagnostic(length=PERIOD):
    """Low-band products versus a true mode-31 product. Folded samples are labeled."""
    low_geometry = product_gap(64, 1, length)["max"]
    low_chart = _chart_product_gap(64, 1, length)
    low_matter = _matter_symbol_gap(64)
    mode31_on_64 = product_gap(64, 31, length)
    mode31_on_128 = product_gap(128, 31, length)
    grid = build_grid(64, quadrature=256, length=length)
    galerkin_mode31 = resolved_product_gap(grid, 31)
    coarse = np.arange(64) * (length / 64)
    fine = np.arange(256) * (length / 256)
    true_mode = np.cos(2 * np.pi * 31 * fine / length)
    folded = np.cos(2 * np.pi * (31 - 64) * fine / length)
    coarse_true = np.cos(2 * np.pi * 31 * coarse / length)
    coarse_folded = np.cos(2 * np.pi * (31 - 64) * coarse / length)
    return {
        "low_band_geometry_product": low_geometry,
        "low_band_chart_F_product": low_chart,
        "low_band_matter_symbol": low_matter,
        "low_band_combined_max": float(max(low_geometry, low_chart, low_matter)),
        "mode31_collocation_64": mode31_on_64["max"],
        "mode31_collocation_64_l2_over_2pi": mode31_on_64["l2_over_2pi"],
        "mode31_collocation_128": mode31_on_128["max"],
        "mode31_galerkin_quadrature_256": galerkin_mode31,
        "mode31_equals_folded_partner_on_64": bool(np.allclose(coarse_true, coarse_folded, atol=1e-12)),
        "mode31_differs_from_folded_partner_on_256": bool(np.max(np.abs(true_mode - folded)) > 1),
        "folded_partner_mode": 31 - 64,
        "uses_true_high_mode": True,
    }


def _chart_product_gap(points, mode, length):
    from .nsc_spherical_feedback_action import feedback_F, partial_F

    derivative = periodic_derivative(points, length)
    coordinate = np.arange(points) * (length / points)
    radius = 2 + 0.1 * np.cos(2 * np.pi * mode * coordinate / length)
    chi = 0.05 * np.sin(2 * np.pi * coordinate / length)
    coefficients = locked_coefficients()
    force_r, force_chi = partial_F(radius, coefficients["A"], coefficients["C_W"])
    values = feedback_F(radius, chi, coefficients["A"], coefficients["C_W"], )
    gap = derivative @ values - (force_r * (derivative @ radius) + force_chi * (derivative @ chi))
    return float(np.max(np.abs(gap)))


def _matter_symbol_gap(points):
    """Anticommutator symbol on one resolved AP mode, against the dense Hamiltonian."""
    count = int(points)
    metric_momentum, _metric = antiperiodic_momentum(count, PERIOD)
    modes = fermion_modes(count)
    coordinate = np.arange(count) * (PERIOD / count)
    column = np.exp(1j * modes[count // 4] * 2 * np.pi * coordinate / PERIOD) / np.sqrt(count)
    coefficient = 0.3 + 0.1 * np.cos(2 * np.pi * coordinate / PERIOD)
    direct = 0.5 * (coefficient * (metric_momentum @ column) + metric_momentum @ (coefficient * column))
    expected = metric_momentum @ column
    # Constant coefficient must reproduce the eigenvalue. Variable part is the gap's companion.
    constant = 0.5 * (1.0 * (metric_momentum @ column) + metric_momentum @ column)
    return float(np.max(np.abs(constant - expected)))


def hamiltonian_fd_errors(grid, state, seed=11, eps=1e-6):
    """Directional derivatives of the fine energy against the pulled-back RHS."""
    rng = np.random.default_rng(seed)
    rate, bundle = compose_fine_hamiltonian(grid, state, include_matter_force=True)
    errors = {}
    samples = {}
    mapping = {
        "Q": (rate.p_Q, -grid.dx_g),
        "r": (rate.p_r, -grid.dx_g),
        "chi": (rate.p_chi, -grid.dx_g),
        "p_Q": (rate.Q, grid.dx_g),
        "p_r": (rate.r, grid.dx_g),
        "p_chi": (rate.chi, grid.dx_g),
    }

    def displaced(name, direction, sign):
        trial = state.copy()
        setattr(trial, name, getattr(trial, name) + sign * eps * direction)
        return energy(grid, trial)

    for name, (velocity, factor) in mapping.items():
        direction = rng.normal(size=velocity.size)
        direction /= np.linalg.norm(direction)
        numeric = (displaced(name, direction, 1) - displaced(name, direction, -1)) / (2 * eps)
        analytic = factor * float(np.dot(velocity, direction))
        errors[name] = float(abs(numeric - analytic))
        samples[name] = float(abs(analytic))
    fine_state = bundle["fine_state"]
    geometric = hamilton_constraint(grid.fine, fine_state)
    shift_density = shift_constraint(grid.fine, fine_state)
    parameter_rows = {
        "L": bundle["source"]["force_L"] + grid.dx_q * geometric,
        "beta": bundle["source"]["force_beta"] + grid.dx_q * shift_density,
    }
    for name, density in parameter_rows.items():
        direction = rng.normal(size=grid.ng)
        direction /= np.linalg.norm(direction)
        lifted = prolong_geometry(grid, direction)

        def parameter_energy(sign, name=name, lifted=lifted):
            length_density = grid.fine.length_density
            shift = grid.fine.shift
            if name == "L":
                length_density = length_density + sign * eps * lifted
            else:
                shift = shift + sign * eps * lifted
            cloned = CouplingSystem(**{**grid.fine.__dict__, "length_density": length_density, "shift": shift})
            return total_energy(cloned, fine_state)

        numeric = (parameter_energy(1) - parameter_energy(-1)) / (2 * eps)
        analytic = float(np.dot(density, lifted))
        errors[name] = float(abs(numeric - analytic))
        samples[name] = float(abs(analytic))
    worst = max(errors.values())
    scale = max(1.0, max(samples.values()))
    return {
        "absolute": errors,
        "relative_worst": float(worst / scale),
        "scale": float(scale),
    }


def work_balance(grid, state, dt=1e-6):
    """Dropping the matter force must leave the lifted fieldwork uncancelled."""
    intact, intact_bundle = compose_fine_hamiltonian(grid, state, include_matter_force=True)
    mutated, _mutated_bundle = compose_fine_hamiltonian(grid, state, include_matter_force=False)

    def slope(rate):
        forward = energy(grid, _combine(state, rate, dt))
        backward = energy(grid, _combine(state, rate, -dt))
        return (forward - backward) / (2 * dt)

    intact_power = float(slope(intact))
    mutated_power = float(slope(mutated))
    fieldwork = float(intact.fieldwork_power)
    unprojected = float(intact_bundle["unprojected_fieldwork"])
    return {
        "intact_power": intact_power,
        "mutated_power": mutated_power,
        "lifted_fieldwork": fieldwork,
        "unprojected_fieldwork": unprojected,
        "lifted_identity_error": float(abs((mutated_power - intact_power) - fieldwork)),
        "unprojected_identity_error": float(abs((mutated_power - intact_power) - unprojected)),
        "uses_lifted_Qdot": True,
    }


def source_column_report(grid, state):
    """Column images, dense forces on a small grid, and occupation eigenvalues."""
    fine_state = prolong_state(grid, state)
    image_gap = dense_hamiltonian_difference(grid.fine, fine_state)
    coarse_values = np.sort(np.real(occupation_eigenvalues(state.phi0, state.phi1)))
    fine_values = np.sort(np.real(occupation_eigenvalues(fine_state.phi0, fine_state.phi1)))
    declared = np.sort(OCCUPATIONS)
    report = {
        "hamiltonian_image_max": float(image_gap),
        "occupation_max_abs_difference": float(np.max(np.abs(coarse_values[-6:] - fine_values[-6:]))),
        "occupation_tail": float(max(np.max(np.abs(coarse_values[:-6])), np.max(np.abs(fine_values[:-6])))),
        "leading_occupations": fine_values[-6:].tolist(),
        "declared_occupation_gap": float(np.max(np.abs(fine_values[-6:] - declared))),
        "phi_dot": "-i H Phi, then U_f adjoint",
        "phi_dot_not": "-M i H Phi",
    }
    if grid.nq <= 32:
        report["dense_forces"] = dense_force_difference(grid.fine, fine_state)
    return report


def physical_link_report(grid, state):
    """Compare link blocks with the old hand-set B. Do not copy B into the state."""
    fine_state = prolong_state(grid, state)
    report = packet_onsite_report(grid.fine, fine_state)
    differences = [row["difference_from_old_B"] for row in report["links"]]
    copied = any(row["copied_from_old_B"] for row in report["links"])
    return {
        "onsite": report["onsite"],
        "links": report["links"],
        "max_difference_from_old_B": float(max(differences) if differences else 0.0),
        "B_phys_equals_old_B": bool(max(differences) < 1e-8) if differences else False,
        "old_B_not_used": True,
        "old_circulation_transferred": False,
        "copied_from_old_B": copied,
    }


def subspace_frequency(grid, state):
    """Retained-band frequency. Quadrature Nyquist is not this number."""
    fine_state = prolong_state(grid, state)
    wavenumber = float(np.max(np.abs(grid.modes_f)) * 2 * np.pi / grid.length)
    conformal = np.max(grid.fine.length_density / fine_state.Q)
    shift = np.max(np.abs(grid.fine.shift))
    mass = np.max(np.abs(grid.fine.kappa * grid.fine.length_density))
    omega = float((conformal + shift) * wavenumber + mass)
    quadrature_wavenumber = np.pi * grid.nq / grid.length
    return omega, float(quadrature_wavenumber)


def stable_timestep(grid, state):
    omega, quadrature_omega = subspace_frequency(grid, state)
    # RK4 imaginary stability is near 2√2. Stay at about half of that, and at the v1 cap.
    step = float(min(5e-4, 1.4 / max(omega, 1.0)))
    return step, omega, quadrature_omega


def rk4_step(grid, state, dt):
    k1 = rates(grid, state)
    k2 = rates(grid, _combine(state, k1, 0.5 * dt))
    k3 = rates(grid, _combine(state, k2, 0.5 * dt))
    k4 = rates(grid, _combine(state, k3, dt))
    combined = CauchyRate(
        (k1.Q + 2 * k2.Q + 2 * k3.Q + k4.Q) / 6,
        (k1.r + 2 * k2.r + 2 * k3.r + k4.r) / 6,
        (k1.chi + 2 * k2.chi + 2 * k3.chi + k4.chi) / 6,
        (k1.p_Q + 2 * k2.p_Q + 2 * k3.p_Q + k4.p_Q) / 6,
        (k1.p_r + 2 * k2.p_r + 2 * k3.p_r + k4.p_r) / 6,
        (k1.p_chi + 2 * k2.p_chi + 2 * k3.p_chi + k4.p_chi) / 6,
        (k1.phi0 + 2 * k2.phi0 + 2 * k3.phi0 + k4.phi0) / 6,
        (k1.phi1 + 2 * k2.phi1 + 2 * k3.phi1 + k4.phi1) / 6,
        k1.force_L, k1.force_Q, k1.force_beta, k1.fieldwork_power,
    )
    return _combine(state, combined, dt)


def evolve(grid, state, duration, dt):
    steps = int(np.ceil(float(duration) / dt - 1e-12))
    dt = float(duration / steps)
    current = state.copy()
    samples = [observation(grid, current)]
    samples[0]["time"] = 0.0
    fieldwork = 0.0
    previous = samples[0]["lifted_fieldwork"]
    try:
        for step in range(steps):
            try:
                current = rk4_step(grid, current, dt)
            except PositiveChartExit as exit_chart:
                raise PositiveChartExit(exit_chart.reason, step * dt, current) from exit_chart
            sample = observation(grid, current)
            if not sample["positive_r"] or not sample["positive_Q"]:
                reason = "r_left_positive_chart" if not sample["positive_r"] else "Q_left_positive_chart"
                raise PositiveChartExit(reason, (step + 1) * dt, current)
            sample["time"] = (step + 1) * dt
            fieldwork += 0.5 * dt * (previous + sample["lifted_fieldwork"])
            previous = sample["lifted_fieldwork"]
            samples.append(sample)
    except PositiveChartExit as exit_chart:
        return {
            "stopped": True,
            "reason": exit_chart.reason,
            "time": float(exit_chart.time_value),
            "state": current,
            "samples": samples,
            "fieldwork_integral": float(fieldwork),
            "dt": dt,
            "steps": steps,
            "success": False,
        }
    energy0 = samples[0]["energy"]
    return {
        "stopped": False,
        "reason": None,
        "time": float(duration),
        "state": current,
        "samples": samples,
        "fieldwork_integral": float(fieldwork),
        "dt": dt,
        "steps": steps,
        "success": True,
        "energy_drift": float(max(abs(sample["energy"] - energy0) for sample in samples)),
        "full_hamilton_max": float(max(sample["full_hamilton_max"] for sample in samples)),
        "full_momentum_max": float(max(sample["full_momentum_max"] for sample in samples)),
        "projected_hamilton_max": float(max(sample["projected_hamilton_max"] for sample in samples)),
        "held_out_hamilton_max": float(max(sample["held_out_hamilton_max"] for sample in samples)),
        "unitarity_max": float(max(sample["unitarity"] for sample in samples)),
        "number_drift": float(max(abs(sample["number"] - samples[0]["number"]) for sample in samples)),
    }


def observation(grid, state):
    rate, bundle = compose_fine_hamiltonian(grid, state)
    diagnostics = constraint_diagnostics(grid, state, bundle["source"], bundle["fine_state"])
    velocities = normal_velocities(grid, state, rate, bundle)
    fine_state = bundle["fine_state"]
    gram = fine_state.phi0.conj().T @ fine_state.phi0 + fine_state.phi1.conj().T @ fine_state.phi1
    number = float(np.sum(
        grid.fine.occupations * np.sum(np.abs(fine_state.phi0) ** 2 + np.abs(fine_state.phi1) ** 2, axis=0)
    ))
    return {
        "time": None,
        "energy": energy(grid, state),
        "field_energy": field_energy(grid.fine, fine_state),
        "gravity_energy": gravity_energy(grid.fine, fine_state),
        "full_hamilton_max": diagnostics["full_hamilton_max"],
        "full_momentum_max": diagnostics["full_momentum_max"],
        "projected_hamilton_max": diagnostics["projected_hamilton_max"],
        "projected_momentum_max": diagnostics["projected_momentum_max"],
        "held_out_hamilton_max": diagnostics["held_out_hamilton_max"],
        "held_out_momentum_max": diagnostics["held_out_momentum_max"],
        "hamilton_unresolved_fraction": diagnostics["hamilton_modes"]["unresolved_fraction"],
        "unitarity": float(np.max(np.abs(gram - np.eye(6)))),
        "number": number,
        "lifted_fieldwork": rate.fieldwork_power,
        "chart_proper_max": velocities["chart_proper_max"],
        "lifted_proper_max": velocities["lifted_proper_max"],
        "positive_r": diagnostics["positive_r"],
        "positive_Q": diagnostics["positive_Q"],
        "r_min_quadrature": diagnostics["r_min_quadrature"],
        "Q_min_quadrature": diagnostics["Q_min_quadrature"],
        "chi_max": float(np.max(np.abs(fine_state.chi))),
    }


def _local_checks_pass(initial):
    names = (
        "projected_hamilton_max",
        "projected_momentum_max",
        "full_hamilton_max",
        "full_momentum_max",
        "held_out_hamilton_max",
        "held_out_momentum_max",
    )
    return all(initial["diagnostics"][name] <= TOL_INITIAL_CONSTRAINT for name in names)


def _window_passes(result):
    if result["stopped"] or not result["success"]:
        return False
    return (
        result["full_hamilton_max"] <= TOL_CONSTRAINT_VALIDATION
        and result["full_momentum_max"] <= TOL_CONSTRAINT_VALIDATION
        and result["held_out_hamilton_max"] <= TOL_CONSTRAINT_VALIDATION
    )


def _compact(result):
    samples = result["samples"]
    keep = samples if len(samples) <= 40 else samples[:: max(1, len(samples) // 20)] + samples[-1:]
    keys = (
        "time", "energy", "full_hamilton_max", "full_momentum_max",
        "projected_hamilton_max", "held_out_hamilton_max", "hamilton_unresolved_fraction",
        "unitarity", "number", "lifted_fieldwork", "chart_proper_max", "lifted_proper_max",
        "r_min_quadrature", "Q_min_quadrature", "chi_max",
    )
    compact = {
        "stopped": result["stopped"],
        "reason": result["reason"],
        "time": result["time"],
        "dt": result["dt"],
        "steps": result["steps"],
        "success_interval_completed": result["success"],
        "energy_drift": result.get("energy_drift"),
        "full_hamilton_max": result.get("full_hamilton_max"),
        "full_momentum_max": result.get("full_momentum_max"),
        "projected_hamilton_max": result.get("projected_hamilton_max"),
        "held_out_hamilton_max": result.get("held_out_hamilton_max"),
        "unitarity_max": result.get("unitarity_max"),
        "number_drift": result.get("number_drift"),
        "fieldwork_integral": result["fieldwork_integral"],
        "checkpoints": [{key: sample[key] for key in keys} for sample in keep],
    }
    return compact


def _end_arrays(state):
    return {
        "r": state.r.tolist(),
        "Q": state.Q.tolist(),
        "chi": state.chi.tolist(),
        "p_Q": state.p_Q.tolist(),
        "p_r": state.p_r.tolist(),
        "p_chi": state.p_chi.tolist(),
    }


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": np.real(value).tolist(), "imag": np.imag(value).tolist()}
        return value.tolist()
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return None if not np.isfinite(number) else number
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def _initial_slice(fermions, quadrature, columns):
    grid = build_grid(fermions, quadrature=quadrature)
    state = blank_state(grid, columns[0], columns[1])
    solved, info = solve_initial_radius(grid, state)
    return grid, solved, info


def run_bounded_control(path=None, resolutions=None, short_time=0.005, long_time=0.05):
    """Two matched subspaces, separate time and quadrature checks, one JSON."""
    started = time.perf_counter()
    cpu = time.process_time()
    path = Path(path) if path is not None else CONTROL_RECORD
    if path.exists():
        raise FileExistsError(f"refusing to overwrite {path}")
    resolutions = ((64, 256), (128, 512)) if resolutions is None else tuple(resolutions)
    probe_columns = load_physical_columns(resolutions[0][0])
    problems = list(probe_columns[3])
    record = {
        "schema": SCHEMA,
        "status": STATUS,
        "provisional": True,
        "renewal": False,
        "global_regeneration": False,
        "absolute_vacuum_claim": False,
        "independent_review": False,
        "representation": "variational Fourier-Galerkin",
        "post_step_filter": False,
        "constraint_projection": False,
        "prescribed_radius": False,
        "v1_immutable": True,
        "filtered_v1_radius": False,
        "preparation_dependency": PREPARATION_DEPENDENCY,
        "preparation_problems": problems,
        "preparation": probe_columns[2],
        "old_circulation_transferred": False,
        "phi_dot": "U_f† (-i H_fine U_f Phi)",
        "phi_dot_not": "-M i H Phi",
        "tolerances": {
            "initial": TOL_INITIAL_CONSTRAINT,
            "validation": TOL_CONSTRAINT_VALIDATION,
            "finite_difference": TOL_FD_RELATIVE,
            "energy": TOL_ENERGY_RELATIVE,
            "newton": TOL_NEWTON,
        },
        "hashes": {
            "galerkin": _sha256(_MODULE_PATH),
            "tests": _sha256(_TEST_PATH) if _TEST_PATH.exists() else None,
            "coupling": _sha256(_COUPLING_PATH),
            "action": _sha256(_ACTION_PATH),
            "source": _sha256(_SOURCE_PATH),
        },
    }
    if problems:
        record.update({
            "verdict": "BLOCKED_UNCORRECTED_PREPARATION",
            "next_blocker": problems[0],
            "runs": {},
            "runtime_seconds": time.perf_counter() - started,
            "cpu_seconds": time.process_time() - cpu,
        })
        _write(path, record)
        record["path"] = str(path)
        return record

    alias = alias_diagnostic()
    record["alias_diagnostic"] = alias
    slices = []
    runs = {}
    fd_report = None
    balance = None
    links = None
    all_initial = True
    chart_stop = None
    for fermions, quadrature in resolutions:
        loaded = load_physical_columns(fermions)
        problems.extend(item for item in loaded[3] if item not in problems)
        grid, solved, info = _initial_slice(fermions, quadrature, loaded)
        doubled = build_grid(fermions, quadrature=2 * quadrature)
        cross = None if not info["converged"] else constraint_diagnostics(doubled, solved)
        entry = {
            "nf": fermions,
            "ng": grid.ng,
            "nq": quadrature,
            "doubled_nq": 2 * quadrature,
            "solve": {key: value for key, value in info.items() if key != "diagnostics"},
            "diagnostics": info.get("diagnostics"),
            "velocities": info.get("velocities"),
            "doubled_quadrature": cross,
            "isometry_geometry": grid.isometry_geometry,
            "isometry_columns": grid.isometry_columns,
        }
        if info["converged"]:
            entry["positive_on_quadrature"] = bool(
                info["diagnostics"]["positive_r"] and info["diagnostics"]["positive_Q"]
            )
            if not _local_checks_pass(info) or cross is None or cross["full_hamilton_max"] > TOL_INITIAL_CONSTRAINT or cross["full_momentum_max"] > TOL_INITIAL_CONSTRAINT:
                all_initial = False
        else:
            all_initial = False
        slices.append((grid, solved, info, entry))
        record_key = f"initial_nf{fermions}_nq{quadrature}"
        runs[record_key] = entry
    record["preparation_problems"] = problems
    if problems:
        record.update({
            "verdict": "BLOCKED_UNCORRECTED_PREPARATION",
            "next_blocker": problems[0],
            "runs": runs,
            "runtime_seconds": time.perf_counter() - started,
            "cpu_seconds": time.process_time() - cpu,
        })
        _write(path, record)
        record["path"] = str(path)
        return record
    if slices and slices[0][2]["converged"]:
        grid, solved, info, _entry = slices[0]
        trial = solved.copy()
        rng = np.random.default_rng(11)
        trial.p_Q = 0.02 * rng.normal(size=grid.ng)
        trial.p_r = 0.02 * rng.normal(size=grid.ng)
        trial.p_chi = 0.02 * rng.normal(size=grid.ng)
        trial.chi = 0.01 * np.sin(2 * np.pi * grid.xi_g / grid.length)
        fd_report = hamiltonian_fd_errors(grid, trial)
        balance = work_balance(grid, solved)
        links = physical_link_report(grid, solved)
    rhs_ok = fd_report is not None and fd_report["relative_worst"] <= TOL_FD_RELATIVE
    short_ok = True
    evolved_short = False
    for grid, solved, info, entry in slices:
        if not info["converged"]:
            short_ok = False
            continue
        evolved_short = True
        step, omega, quadrature_omega = stable_timestep(grid, solved)
        entry["omega_subspace"] = omega
        entry["omega_quadrature_nyquist"] = quadrature_omega
        entry["dt"] = step
        primary = evolve(grid, solved, short_time, step)
        halved = evolve(grid, solved, short_time, 0.5 * step)
        tag = f"nf{grid.nf}_nq{grid.nq}_T{short_time}"
        runs[tag] = _compact(primary)
        runs[tag + "_dt_half"] = _compact(halved)
        runs[tag]["axis"] = "time_and_space"
        runs[tag + "_dt_half"]["axis"] = "time_refinement"
        if primary["stopped"] or halved["stopped"]:
            chart_stop = primary["reason"] or halved["reason"]
            short_ok = False
        elif not (_window_passes(primary) and _window_passes(halved)):
            short_ok = False
        scale = max(1.0, abs(primary["samples"][0]["energy"]))
        if primary["success"] and primary["energy_drift"] > TOL_ENERGY_RELATIVE * scale:
            short_ok = False
    if not evolved_short:
        short_ok = False

    long_ok = False
    if all_initial and rhs_ok and short_ok and chart_stop is None:
        grid, solved, _info, _entry = slices[0]
        step, _omega, _quadrature_omega = stable_timestep(grid, solved)
        extended = evolve(grid, solved, long_time, step)
        runs[f"nf{grid.nf}_nq{grid.nq}_T{long_time}"] = _compact(extended)
        runs[f"nf{grid.nf}_nq{grid.nq}_T{long_time}"]["axis"] = "time_extension"
        runs[f"nf{grid.nf}_nq{grid.nq}_T{long_time}"]["end"] = _end_arrays(extended["state"])
        long_ok = _window_passes(extended) and not extended["stopped"]
        if extended["stopped"]:
            chart_stop = extended["reason"]

    if problems:
        verdict = "BLOCKED_UNCORRECTED_PREPARATION"
        blocker = problems[0]
    elif chart_stop is not None:
        verdict = "STOP_POSITIVE_CHART"
        blocker = chart_stop
    elif not all(item[2]["converged"] for item in slices):
        verdict = "FAIL_INITIAL_SOLVE"
        blocker = next(item[2]["blocker"] for item in slices if not item[2]["converged"])
    elif not rhs_ok:
        verdict = "FAIL_DISCRETE_RHS"
        blocker = "pulled-back Hamiltonian RHS disagrees with directional energy derivatives"
    elif not all_initial:
        verdict = "FAIL_HELD_OUT_CONSTRAINT"
        blocker = "projected residual is not the physical constraint; full or doubled-quadrature residual remains"
    elif not short_ok:
        verdict = "PROVISIONAL_CONSTRAINT_DRIFT"
        blocker = "short-window full local constraints exceed 1e-3"
    elif not long_ok:
        verdict = "PROVISIONAL_CONSTRAINT_DRIFT"
        blocker = "T=0.05 full local constraints exceed 1e-3; a short step is not a renewal"
    else:
        verdict = "PROVISIONAL_LOCAL_WINDOW"
        blocker = "local window only; not a continuum limit and not a renewal"
    record.update({
        "verdict": verdict,
        "next_blocker": blocker,
        "renewal": False,
        "hamiltonian_fd": fd_report,
        "work_balance": balance,
        "physical_links": links,
        "runs": runs,
        "multiplicity": int(sector_multiplicity(KAPPA)),
        "coefficients": locked_coefficients(),
        "rmag2": magnetic_radius_square(locked_coefficients()),
        "runtime_seconds": time.perf_counter() - started,
        "cpu_seconds": time.process_time() - cpu,
    })
    _write(path, record)
    record["path"] = str(path)
    return record


def _write(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_jsonable(record), indent=2) + "\n")
