"""A nested canonical frame of the existing one-action spherical model.

There is one physical metric and one Gaussian source, not two copied boxes.
The strict coordinate domains are I_C=(1,3) inside I_P=(0,4), on the
period-eight carrier.  A fixed real orthogonal geometry frame has coarse
parent modes and child-localized detail modes; all remaining geometry and
all ambient fermion modes remain dynamical.  Localization after Fourier
interpolation is measured, not asserted to be exact compact support.

For each geometry variable g, a=W.T@g and pi=dx_g*W.T@p are canonical:
dx_g p.dg = pi.da.  The Hamiltonian and rates below are the existing
conformal Galerkin action in these coordinates, including its lapse chain
rule and F_L Ldot work.  No instantaneous spectral inheritance, physical
homothety, running coupling, separate spacetime, or new interface force is
asserted.  ``functional_pullback`` is an exact affine coordinate/clock map.
"""
from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin


GEOMETRY_NAMES = ("Q", "r", "chi")
MOMENTUM_NAMES = ("p_Q", "p_r", "p_chi")
STATE_NAMES = GEOMETRY_NAMES + MOMENTUM_NAMES + ("phi0", "phi1")
PARENT_INTERVAL = (0.0, 4.0)
CHILD_INTERVAL = (1.0, 3.0)
CHILD_INDICES = (2, 3)


@dataclass
class NestedState:
    """Geometry coefficients and canonical pi; full ambient source columns.

    p_Q,p_r,p_chi are pi=dx_g W.T p_nodal, NOT the old nodal momenta.
    phi0/phi1 retain the original nf-by-six half-density representation.
    """

    Q: np.ndarray
    r: np.ndarray
    chi: np.ndarray
    p_Q: np.ndarray
    p_r: np.ndarray
    p_chi: np.ndarray
    phi0: np.ndarray
    phi1: np.ndarray

    def copy(self):
        return NestedState(*(np.array(getattr(self, name), copy=True) for name in STATE_NAMES))


@dataclass
class NestedRate(NestedState):
    fieldwork_power: float
    force_L: np.ndarray
    force_Q: np.ndarray
    force_beta: np.ndarray


@dataclass(frozen=True)
class NestedPair:
    grid: galerkin.GalerkinGrid
    geometry_map: np.ndarray
    geometry_coarse_indices: np.ndarray
    geometry_child_indices: np.ndarray
    geometry_parent_indices: np.ndarray
    original_phi0: np.ndarray
    original_phi1: np.ndarray
    original_columns: np.ndarray
    source_phi0: np.ndarray
    source_phi1: np.ndarray
    source_columns: np.ndarray
    weights: np.ndarray
    source_metadata: dict
    geometry_metadata: dict
    parent_interval: tuple = PARENT_INTERVAL
    child_interval: tuple = CHILD_INTERVAL
    child_indices: tuple = CHILD_INDICES
    clock_locations: tuple = (1.0, 2.0, 3.0)

    @property
    def reference_columns(self):
        """Fixed overlapping observer; distinct from separated source."""
        return self.original_columns


def _readonly(value, dtype=None):
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def _orthogonalize(vector, previous):
    out = np.array(vector, dtype=float, copy=True)
    for _ in range(2):
        for basis in previous:
            out -= basis * np.dot(basis, out)
    norm = np.linalg.norm(out)
    if norm <= 1e-11:
        return None
    return out / norm


def _geometry_frame(grid, coarse_modes=2, child_details=4):
    coarse_modes, child_details = int(coarse_modes), int(child_details)
    if coarse_modes < 0 or child_details < 1 or 1 + 2 * coarse_modes + child_details >= grid.ng:
        raise ValueError("coarse modes and child details must fit the full geometry band")
    x = grid.xi_g
    frame = [np.ones(grid.ng) / np.sqrt(grid.ng)]
    for mode in range(1, coarse_modes + 1):
        angle = 2 * np.pi * mode * x / grid.length
        frame.extend((np.sqrt(2 / grid.ng) * np.cos(angle), np.sqrt(2 / grid.ng) * np.sin(angle)))
    coarse_count = len(frame)
    u = (x - CHILD_INTERVAL[0]) / (CHILD_INTERVAL[1] - CHILD_INTERVAL[0])
    bump = coupling._bump(u)
    # Remove coarse moments INSIDE the child support. Subtracting global
    # coarse vectors from a bump would create nonlocal detail components.
    support_count = int(np.count_nonzero(bump))
    seed_count = min(support_count, coarse_count + child_details + 2)
    if seed_count < coarse_count + child_details:
        raise ValueError("child geometry nullspace is unresolved on this band")
    seeds = []
    for degree in range(seed_count):
        polynomial = np.polynomial.legendre.legval(2 * u - 1, [0.0] * degree + [1.0])
        seed = bump * polynomial
        norm = np.linalg.norm(seed)
        if norm <= 1e-12:
            raise ValueError("child geometry seed is unresolved")
        seeds.append(seed / norm)
    seeds = np.column_stack(seeds)
    coarse_matrix = np.column_stack(frame)
    moment_matrix = coarse_matrix.T @ seeds
    _left, singular, right = np.linalg.svd(moment_matrix, full_matrices=True)
    rank = int(np.sum(singular > 1e-11 * max(1., singular[0])))
    nullspace = right[rank:].T
    localized = seeds @ nullspace
    local_left, local_singular, _local_right = np.linalg.svd(localized, full_matrices=False)
    if local_singular.size < child_details or local_singular[child_details - 1] <= 1e-10:
        raise ValueError("child geometry coarse-moment nullspace is unresolved")
    local_vectors = local_left[:, :child_details]
    # SVD can put roundoff on zero rows. Enforce exactly the declared nodal
    # support, then orthonormalize only WITHIN that support (never globally).
    local_vectors[bump == 0] = 0.
    child_frame = []
    for index in range(child_details):
        candidate = _orthogonalize(local_vectors[:, index], child_frame)
        if candidate is None:
            raise ValueError("child geometry details are linearly dependent")
        child_frame.append(candidate)
    frame.extend(child_frame)
    child_ids = np.arange(coarse_count, len(frame), dtype=int)
    # Complete the same ambient geometry band. No omitted/filter degrees.
    for column in np.eye(grid.ng):
        if len(frame) == grid.ng:
            break
        candidate = _orthogonalize(column, frame)
        if candidate is not None:
            frame.append(candidate)
    if len(frame) != grid.ng:
        raise ValueError("could not complete the full real geometry frame")
    matrix = np.column_stack(frame)
    error = float(np.max(np.abs(matrix.T @ matrix - np.eye(grid.ng))))
    if error > 1e-9:
        raise ValueError("real geometry frame is not orthogonal")
    fine_details = grid.A_g @ matrix[:, child_ids]
    outside = (grid.xi_q < CHILD_INTERVAL[0]) | (grid.xi_q > CHILD_INTERVAL[1])
    leakage = np.sum(fine_details[outside] ** 2, axis=0) / np.sum(fine_details ** 2, axis=0)
    parent_ids = np.setdiff1d(np.arange(grid.ng), child_ids)
    return matrix, np.arange(coarse_count), child_ids, parent_ids, {
        "orthogonality_max": error,
        "coarse_modes": coarse_modes,
        "child_details": child_details,
        "child_seed_nullity": int(nullspace.shape[1]),
        "child_coarse_moment_max": float(np.max(np.abs(coarse_matrix.T @ np.column_stack(child_frame)))),
        "child_nodal_outside_power_fraction": [0.0] * child_details,
        "child_detail_outside_power_fraction": leakage.tolist(),
        "exact_spatial_support": False,
        "all_geometry_degrees_retained": True,
    }


def _separated_columns(grid):
    old0, old1, preparation = coupling.prepare_rank6(grid.xi_f, grid.dx_f)
    old = np.vstack((old0, old1))
    child = old[:, CHILD_INDICES].copy()
    phase = np.exp(1j * coupling.CALIBRATION["phase"])
    norms = coupling._lobe_norms()
    outer = []
    for left, right in ((0.0, 1.0), (3.0, 4.0)):
        lobe_width = (right - left) / 2
        even_u = (grid.xi_f - left) / lobe_width
        odd_u = (grid.xi_f - left - lobe_width) / lobe_width
        even = coupling._bump(even_u) / (np.sqrt(lobe_width) * norms[0])
        odd = (odd_u - .5) * coupling._bump(odd_u) / (np.sqrt(lobe_width) * norms[1])
        even *= np.sqrt(grid.dx_f)
        odd *= np.sqrt(grid.dx_f)
        if min(np.linalg.norm(even), np.linalg.norm(odd)) <= 1e-12:
            raise ValueError("separated source lobe is unresolved")
        envelope = (even / np.linalg.norm(even) + odd / np.linalg.norm(odd)) / np.sqrt(2)
        carrier = np.exp(1j * coupling.CARRIER_K * (grid.xi_f - left)) * envelope
        plus = np.concatenate((carrier, 1j * carrier)) / np.sqrt(2)
        outer.extend((plus, phase * np.conjugate(plus)))
    outer = np.column_stack(outer)
    outer -= child @ (child.conj().T @ outer)
    outer, _ = coupling.lowdin(outer)
    columns = np.column_stack((outer[:, 0:2], child, outer[:, 2:4]))
    return columns[:grid.nf], columns[grid.nf:], {
        "layout": "separated",
        "source_support_closures": [[0., 1.], [1., 3.], [3., 4.]],
        "middle_columns_equal_owned_original": bool(np.array_equal(columns[:, 2:4], child)),
        "historical_source_reused_as_trajectory": False,
        "phase_convention": preparation["minus_columns"],
        "scaled_lobe_half_density": True,
        "outer_projected_against_fixed_child": True,
    }


def build_pair(nf=64, *, quadrature=None, coarse_modes=2, child_details=4,
               source_layout="separated", columns_override=None, occupations=None):
    """Construct frames and source ONLY; no evolution or initial solve.

    ``separated`` is a new explicit preparation. Its owned original middle
    pair remains fixed; the surrounding packets occupy (0,1) and (3,4).
    ``original`` is an explicit comparison, never a silent replay binding.
    """
    if int(nf) != nf or nf < 32 or nf % 2:
        raise ValueError("nested pair needs an even nf of at least 32")
    grid = galerkin.build_grid(int(nf), quadrature=quadrature, gauge="conformal")
    reference0, reference1, _reference_preparation = coupling.prepare_rank6(grid.xi_f, grid.dx_f)
    reference = np.vstack((reference0, reference1))
    if columns_override is not None:
        if source_layout != "override":
            raise ValueError("columns_override requires source_layout='override'")
        phi0, phi1 = (np.array(value, dtype=complex, copy=True) for value in columns_override)
        metadata = {"layout": "override", "caller_declared_preparation": True}
    elif source_layout == "separated":
        phi0, phi1, metadata = _separated_columns(grid)
    elif source_layout == "original":
        phi0, phi1, metadata = galerkin.load_physical_columns(grid.nf)[:3]
        metadata = dict(metadata, layout="original", source_support_closures=[[0., 2.], [1., 3.], [2., 4.]])
    else:
        raise ValueError("source_layout must be separated, original, or explicit override")
    if phi0.shape != (grid.nf, 6) or phi1.shape != (grid.nf, 6):
        raise ValueError("source must have nf-by-six spinor arrays")
    columns = np.vstack((phi0, phi1))
    if not np.isfinite(columns).all():
        raise ValueError("source contains nonfinite columns")
    gram_gap = float(np.max(np.abs(columns.conj().T @ columns - np.eye(6))))
    if gram_gap > 1e-9:
        raise ValueError("source columns must be orthonormal, without repairing the child frame")
    weights = np.array(coupling.OCCUPATIONS if occupations is None else occupations, dtype=float)
    if weights.shape != (6,) or not np.isfinite(weights).all() or np.min(weights) <= 0 or np.max(weights) > 1:
        raise ValueError("six positive CAR occupation weights in (0,1] are required")
    grid.fine = replace(grid.fine, occupations=weights.copy())
    fine0, fine1 = grid.U_f @ phi0, grid.U_f @ phi1
    fine_power = np.abs(fine0) ** 2 + np.abs(fine1) ** 2
    support = metadata.get("source_support_closures", [])
    leakage = []
    for packet, (left, right) in enumerate(support):
        outside = (grid.xi_q < left) | (grid.xi_q > right)
        for column in (2 * packet, 2 * packet + 1):
            leakage.append(float(np.sum(fine_power[outside, column]) / np.sum(fine_power[:, column])))
    metadata = dict(metadata, gram_max=gram_gap, occupations=weights.tolist(),
                    total_occupation=float(np.sum(weights)),
                    fine_outside_support_power_fraction=leakage,
                    exact_spatial_support=False, complement_occupation=0.0,
                    physical_vacuum_identification=False,
                    reference_observer="owned original overlapping rank6, fixed T0",
                    source_equals_reference=bool(np.array_equal(columns, reference)),
                    child_source_equals_reference=bool(np.array_equal(columns[:, 2:4], reference[:, 2:4])))
    matrix, coarse, child, parent, geometry_meta = _geometry_frame(grid, coarse_modes, child_details)
    return NestedPair(grid, _readonly(matrix), _readonly(coarse, int), _readonly(child, int),
                      _readonly(parent, int), _readonly(reference0), _readonly(reference1),
                      _readonly(reference), _readonly(phi0), _readonly(phi1),
                      _readonly(columns), _readonly(weights), metadata, geometry_meta)


def _check(pair, state):
    for name in GEOMETRY_NAMES + MOMENTUM_NAMES:
        array = np.asarray(getattr(state, name))
        if array.shape != (pair.grid.ng,) or np.iscomplexobj(array) or not np.isfinite(array).all():
            raise ValueError(f"{name} must be a finite real ng-vector")
    for name in ("phi0", "phi1"):
        array = np.asarray(getattr(state, name))
        if array.shape != (pair.grid.nf, 6) or not np.isfinite(array).all():
            raise ValueError(f"{name} must be a finite nf-by-six array")


def encode_state(pair, nodal_state):
    """Old nodal state to independent canonical parent/detail coordinates."""
    matrix, spacing = pair.geometry_map, pair.grid.dx_g
    values = [matrix.T @ np.asarray(getattr(nodal_state, name)) for name in GEOMETRY_NAMES]
    values += [spacing * matrix.T @ np.asarray(getattr(nodal_state, name)) for name in MOMENTUM_NAMES]
    values += [np.array(nodal_state.phi0, copy=True), np.array(nodal_state.phi1, copy=True)]
    out = NestedState(*values)
    _check(pair, out)
    return out


def reconstruct_state(pair, state):
    """The ONE physical Cauchy state. Neither metric nor matter is copied."""
    _check(pair, state)
    matrix, spacing = pair.geometry_map, pair.grid.dx_g
    values = [matrix @ getattr(state, name) for name in GEOMETRY_NAMES]
    values += [matrix @ getattr(state, name) / spacing for name in MOMENTUM_NAMES]
    values += [np.array(state.phi0, copy=True), np.array(state.phi1, copy=True)]
    return coupling.CauchyState(*values)


def _pairwise_radius_residual(grid, radius, rho):
    """Same finite SBP polynomial, with an independent reduction ordering.

    This is an arithmetic indicator, not higher precision or an enclosure.
    In particular np.longdouble is only float64 on some supported hosts.
    """
    matvec = lambda matrix, vector: np.sum(matrix * vector[None, :], axis=1)
    fine_radius = matvec(grid.A_g, radius)
    system = grid.fine
    Q = np.full(grid.nq, system.calibration["b0"] / system.calibration["a0"])
    F = coupling.feedback_F(fine_radius, np.zeros(grid.nq), system.A, system.C_W)
    V = coupling.feedback_V(fine_radius, np.zeros(grid.nq), system.A, system.C_W,
                            system.C_F, system.flux)
    derivative = matvec(system.derivative, fine_radius)
    F_x = matvec(system.derivative, F)
    full = float(coupling.feedback_Z(system.A)) * derivative ** 2 / Q - Q * V - \
           2 * matvec(system.derivative, F_x / Q) + rho
    return grid.weight * matvec(grid.A_g.T, full)


def _radius_correction(grid, radius, rho):
    residual, full = galerkin.projected_radius_operator(grid, radius, rho)
    jacobian = galerkin._projected_radius_jacobian(grid, radius)
    correction = np.linalg.solve(jacobian, -residual)
    reordered = _pairwise_radius_residual(grid, radius, rho)
    reordered_correction = np.linalg.solve(jacobian, -reordered)
    noise_correction = np.linalg.solve(jacobian, residual - reordered)
    # ULP floor is in radius units; the other term measures cancellation
    # through this actual Jacobian. Neither is a physical residual tolerance.
    ulp_floor = 8 * float(np.max(np.spacing(np.maximum(1., np.abs(radius)))))
    arithmetic_indicator = float(np.max(np.abs(noise_correction)))
    floor = ulp_floor + arithmetic_indicator
    correction_norm = float(np.max(np.abs(correction)))
    reordered_norm = float(np.max(np.abs(reordered_correction)))
    linear_residual = jacobian @ correction + residual
    denominator = np.linalg.norm(jacobian, np.inf) * correction_norm + np.max(np.abs(residual))
    return residual, full, jacobian, correction, {
        "residual_max": float(np.max(np.abs(residual))),
        "newton_correction_inf": correction_norm,
        "newton_correction_l2": float(np.sqrt(grid.dx_g) * np.linalg.norm(correction)),
        "reordered_correction_inf": reordered_norm,
        "radius_ulp_floor": ulp_floor,
        "arithmetic_correction_indicator": arithmetic_indicator,
        "radius_correction_floor": floor,
        "residual_reduction_difference_inf": float(np.max(np.abs(residual - reordered))),
        "linear_backward_error": float(np.max(np.abs(linear_residual)) / max(denominator, np.finfo(float).tiny)),
        "correction_at_arithmetic_floor": bool(max(correction_norm, reordered_norm) <= floor),
    }


def solve_initial_radius_successor(grid, state):
    """Positive source homotopy stopped by the finite radius correction.

    Uses the OWNED projected residual and Jacobian, without changing them.
    The finite arithmetic stop replaces their historical absolute residual
    line. Full/held-out constraints are observations, never erased. The
    arithmetic reduction comparison is explicitly not a continuum bound.
    """
    if grid.gauge != "conformal":
        raise ValueError("nested successor initializer requires conformal gauge")
    expected_Q = grid.fine.calibration["b0"] / grid.fine.calibration["a0"]
    if not np.array_equal(state.Q, np.full(grid.ng, expected_Q)) or \
       any(np.any(getattr(state, name) != 0) for name in ("chi",) + MOMENTUM_NAMES):
        raise ValueError("radius-only initialization requires owned constant Q, zero chi and momenta")
    fine = galerkin.prolong_state(grid, state)
    source = coupling.source_from_columns(galerkin.active_fine_system(grid, fine), fine)
    rho = source["force_L"] / grid.dx_q
    varied = fine.copy()
    varied.r *= 1.7
    other = coupling.source_from_columns(galerkin.active_fine_system(grid, varied), varied)
    rho_variation = float(np.max(np.abs(other["force_L"] / grid.dx_q - rho)))
    radius = np.ones(grid.ng)
    stages = []
    for scale in np.linspace(0., 1., 9)[1:]:
        target = scale * rho
        converged = False
        status = {}
        for iteration in range(24):
            try:
                residual, _full, jacobian, delta, status = _radius_correction(grid, radius, target)
            except np.linalg.LinAlgError as error:
                return state, {"converged": False, "blocker": f"singular radius Jacobian: {error}",
                               "steps": stages, "method": "source homotopy / radius correction"}
            if status["correction_at_arithmetic_floor"]:
                converged = True
                break
            accepted = False
            factor = 1.
            correction_norm = np.linalg.norm(delta, np.inf)
            for _halving in range(20):
                trial = radius + factor * delta
                if np.isfinite(trial).all() and np.min(galerkin.prolong_geometry(grid, trial)) > 0:
                    trial_residual, _trial_full = galerkin.projected_radius_operator(grid, trial, target)
                    # The raw residual magnifies high-frequency roundoff.
                    # Measure progress in radius units with the CURRENT J,
                    # so a useful low-frequency correction is not rejected.
                    trial_correction = np.linalg.solve(jacobian, -trial_residual)
                    if np.linalg.norm(trial_correction, np.inf) < correction_norm * (1 - .05 * factor):
                        radius = trial
                        accepted = True
                        break
                factor *= .5
            if not accepted:
                break
        stages.append(dict(status, scale=float(scale), iterations=iteration + 1, converged=converged))
        if not converged:
            return state, {"converged": False,
                           "blocker": "positive homotopy stalled above measured radius-correction floor",
                           "steps": stages, "residual_max": status.get("residual_max"),
                           "method": "source homotopy / radius correction",
                           "rho_independent_of_r_max": rho_variation,
                           "imposed_radius": False, "filtered_v1_radius": False}
    solved = state.copy()
    solved.r = radius
    residual, full, jacobian, _correction, status = _radius_correction(grid, radius, rho)
    held = full - galerkin.prolong_geometry(grid, residual)
    diagnostics = dict(galerkin.constraint_diagnostics(grid, solved))
    diagnostics.update(projected_hamilton_max=float(np.max(np.abs(residual))),
                       full_hamilton_max=float(np.max(np.abs(full))),
                       held_out_hamilton_max=float(np.max(np.abs(held))))
    fine = galerkin.prolong_state(grid, solved)
    bracket = fine.Q ** 2 * coupling.magnetic_radius_square(grid.fine.coefficients) + \
              fine.Q * rho / (8 * np.pi * grid.fine.A)
    report = dict(status, converged=True, blocker=None, steps=stages,
                  method="source homotopy / radius correction",
                  jacobian_inf_norm=float(np.linalg.norm(jacobian, np.inf)),
                  jacobian_condition_inf=float(np.linalg.cond(jacobian, np.inf)),
                  rho_independent_of_r_max=rho_variation,
                  imposed_radius=False, filtered_v1_radius=False,
                  r_min_coarse=float(np.min(radius)), r_max_coarse=float(np.max(radius)),
                  bracket_min=float(np.min(bracket)), bracket_positive=bool(np.min(bracket) > 0),
                  diagnostics=diagnostics,
                  approximation_scope="finite projected initial lapse constraint; actual held-out and shift residuals retained",
                  arithmetic_indicator_certified=False, continuum_initial_state_certified=False,
                  stop_rule="both independently ordered Newton corrections <= 8 radius ULP + mapped residual-reduction indicator")
    return solved, report


def initial_state(pair, *, solve_constraints=True):
    """Source-derived initial radius, independently for this source/weights.

    The caller explicitly starts this solve. No old saved state is rebound.
    The owned constructor retains the mean current and its residual domain;
    it does not silently solve a different shift equation or delete that mean.
    """
    nodal = galerkin.blank_state(pair.grid, pair.source_phi0, pair.source_phi1)
    if solve_constraints:
        nodal, report = solve_initial_radius_successor(pair.grid, nodal)
        if not report["converged"]:
            raise ValueError(f"initial source-derived geometry failed: {report.get('blocker')}")
    else:
        report = {"converged": False, "initial_constraints_solved": False,
                  "purpose": "explicit unsolved construction, not a physical episode"}
    encoded = encode_state(pair, nodal)
    reconstructed = reconstruct_state(pair, encoded)
    fine = galerkin.prolong_state(pair.grid, reconstructed)
    active = galerkin.active_fine_system(pair.grid, fine)
    source = coupling.source_from_columns(active, fine)
    report = dict(report, gauge="conformal", initial_radius_method="successor correction-controlled solve of owned projected constraint",
                  shift_method="retained full source and residual; mean current not deleted",
                  source_current_mean=float(np.mean(source["force_beta"] / pair.grid.dx_q)))
    if solve_constraints:
        # The W transform is exact algebraically, but finite arithmetic can
        # change the last bits. Record the actual returned physical state.
        _residual, _full, _jacobian, _delta, actual_correction = _radius_correction(
            pair.grid, reconstructed.r, source["force_L"] / pair.grid.dx_q)
        report["solver_diagnostics"] = report["diagnostics"]
        report["diagnostics"] = galerkin.constraint_diagnostics(pair.grid, reconstructed)
        report["post_canonical_map_radius_correction"] = actual_correction
        report["shift_residual_mean"] = float(np.mean(
            coupling.shift_constraint(active, fine) + source["force_beta"] / pair.grid.dx_q))
    return encoded, report


def rates(pair, state, *, include_matter_force=True, return_bundle=False):
    nodal = reconstruct_state(pair, state)
    old_rate, bundle = galerkin.compose_fine_hamiltonian(pair.grid, nodal, include_matter_force)
    matrix, spacing = pair.geometry_map, pair.grid.dx_g
    values = [matrix.T @ getattr(old_rate, name) for name in GEOMETRY_NAMES]
    values += [spacing * matrix.T @ getattr(old_rate, name) for name in MOMENTUM_NAMES]
    values += [old_rate.phi0, old_rate.phi1]
    result = NestedRate(*values, old_rate.fieldwork_power, old_rate.force_L,
                        old_rate.force_Q, old_rate.force_beta)
    if return_bundle:
        return result, dict(bundle, nodal_state=nodal, nodal_rate=old_rate)
    return result


def _combine(state, rate, factor):
    return NestedState(*(getattr(state, name) + factor * getattr(rate, name) for name in STATE_NAMES))


def rk4_step(pair, state, dt):
    """All stages use their reconstructed full metric/source, with no reset."""
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("a positive finite timestep is required")
    k1 = rates(pair, state)
    k2 = rates(pair, _combine(state, k1, dt / 2))
    k3 = rates(pair, _combine(state, k2, dt / 2))
    k4 = rates(pair, _combine(state, k3, dt))
    out = NestedState(*(getattr(state, name) + dt / 6 *
                        (getattr(k1, name) + 2 * getattr(k2, name) +
                         2 * getattr(k3, name) + getattr(k4, name)) for name in STATE_NAMES))
    # Check the final state as well as every stage; a failure is explicit.
    fine = galerkin.prolong_state(pair.grid, reconstruct_state(pair, out))
    active = galerkin.active_fine_system(pair.grid, fine)
    reason = coupling.chart_failure(active, fine)
    if reason is not None:
        raise coupling.PositiveChartExit(reason, np.nan, out.copy())
    return out


def stable_timestep(pair, state, step_cap=0.001):
    if not np.isfinite(step_cap) or step_cap <= 0:
        raise ValueError("a positive finite step cap is required")
    omega, quadrature_omega = galerkin.subspace_frequency(pair.grid, reconstruct_state(pair, state))
    return min(float(step_cap), 1.4 / max(omega, 1.0)), omega, quadrature_omega


def _active(pair, state):
    fine = galerkin.prolong_state(pair.grid, reconstruct_state(pair, state))
    system = galerkin.active_fine_system(pair.grid, fine)
    reason = coupling.chart_failure(system, fine)
    if reason is not None:
        raise coupling.PositiveChartExit(reason, np.nan, state.copy())
    return fine, system


def apply_hamiltonian(pair, state, columns):
    """Representative H_G on the FULL field band; no angular copy factor."""
    vectors = np.asarray(columns, dtype=complex)
    vector = vectors.ndim == 1
    if vector:
        vectors = vectors[:, None]
    if vectors.ndim != 2 or vectors.shape[0] != 2 * pair.grid.nf or not np.isfinite(vectors).all():
        raise ValueError("Hamiltonian columns must have first dimension 2nf")
    fine, system = _active(pair, state)
    up0 = pair.grid.U_f @ vectors[:pair.grid.nf]
    up1 = pair.grid.U_f @ vectors[pair.grid.nf:]
    image0, image1 = coupling.apply_dirac(up0, up1, system.length_density, fine.Q,
                                       system.shift, system.kappa, system.momentum)
    result = np.vstack((pair.grid.U_f.conj().T @ image0, pair.grid.U_f.conj().T @ image1))
    return result[:, 0] if vector else result


def hamiltonian(pair, state):
    return apply_hamiltonian(pair, state, np.eye(2 * pair.grid.nf, dtype=complex))


def live_blocks(pair, state):
    """Parent contains child; ambient remains in the physical generator."""
    frame = pair.original_columns
    child = frame[:, pair.child_indices]
    detail_indices = [index for index in range(6) if index not in pair.child_indices]
    detail = frame[:, detail_indices]
    image = apply_hamiltonian(pair, state, frame)
    return {
        "parent": frame.conj().T @ image,
        "child": child.conj().T @ image[:, pair.child_indices],
        "parent_detail": detail.conj().T @ image[:, detail_indices],
        "child_parent_link": child.conj().T @ image[:, detail_indices],
        "ambient_image": image - frame @ (frame.conj().T @ image),
        "ambient_retained": True,
    }


def energy(pair, state):
    """Full existing Hamiltonian, once, in canonical hierarchical coordinates."""
    return galerkin.energy(pair.grid, reconstruct_state(pair, state))


def energy_accounting(pair, state):
    fine, system = _active(pair, state)
    columns = np.vstack((state.phi0, state.phi1))
    child = pair.original_columns[:, pair.child_indices]
    ids = [index for index in range(6) if index not in pair.child_indices]
    detail = pair.original_columns[:, ids]
    pieces = {"child": child @ (child.conj().T @ columns),
              "parent_detail": detail @ (detail.conj().T @ columns)}
    pieces["ambient"] = columns - pieces["child"] - pieces["parent_detail"]
    names = tuple(pieces)
    images = apply_hamiltonian(pair, state, np.hstack([pieces[name] for name in names]))
    images = {name: images[:, 6 * index:6 * (index + 1)] for index, name in enumerate(names)}
    factor = system.multiplicity
    weighted = lambda first, second: float(factor * np.dot(pair.weights,
                         np.sum(first.conj() * second, axis=0).real))
    result = {f"field_{name}": weighted(pieces[name], images[name]) for name in names}
    for index, left in enumerate(names):
        for right in names[index + 1:]:
            result[f"field_{left}_{right}_cross"] = 2 * weighted(pieces[left], images[right])
    field = coupling.field_energy(system, fine)
    gravity = coupling.gravity_energy(system, fine)
    pieces_sum = sum(result.values())
    result.update(field=field, gravity=gravity, total=field + gravity,
                  field_closure_error=pieces_sum - field,
                  field_parent=result["field_child"] + result["field_parent_detail"] +
                  result["field_child_parent_detail_cross"],
                  one_source=True, angular_multiplicity_applied_once=True)
    return result


def source_geometry_forces(pair, state):
    """Full cross-complete variational source force in canonical pi units."""
    fine, system = _active(pair, state)
    source = coupling.source_from_columns(system, fine)
    total = source["force_Q"] + source["force_L"]
    canonical_force = pair.geometry_map.T @ (pair.grid.A_g.T @ total)
    return {"Q_energy_gradient": canonical_force,
            "p_Q_source_rate": -canonical_force,
            "parent_p_Q_source_rate": -canonical_force[pair.geometry_parent_indices],
            "child_p_Q_source_rate": -canonical_force[pair.geometry_child_indices],
            "lapse_chain_included": True, "live_links_differentiated": True}


def _periodic_values(grid, values, coordinates):
    array = np.asarray(values)
    if array.shape != (grid.nq,):
        raise ValueError("periodic sampled density must be an nq-vector")
    coefficients = np.fft.fft(array) / grid.nq
    modes = np.fft.fftfreq(grid.nq) * grid.nq
    phase = np.exp(2j * np.pi * np.asarray(coordinates)[:, None] * modes[None, :] / grid.length)
    result = phase @ coefficients
    if np.max(np.abs(result.imag)) > 1e-8 * max(1., np.max(np.abs(result.real))):
        raise ValueError("periodic evaluation did not preserve a real field")
    return result.real


def interval_integral(grid, values, interval):
    """Exact integral of the finite periodic interpolant of this density."""
    left, right = interval
    if not 0 <= left < right <= grid.length:
        raise ValueError("interval must be ordered within the periodic carrier")
    values = np.asarray(values)
    if values.shape != (grid.nq,) or not np.isfinite(values).all():
        raise ValueError("interval density must be a finite nq-vector")
    coefficients = np.fft.fft(values) / grid.nq
    modes = np.fft.fftfreq(grid.nq) * grid.nq
    nonzero = modes != 0
    wave = 2j * np.pi * modes[nonzero] / grid.length
    result = coefficients[0] * (right - left) + np.dot(coefficients[nonzero],
             (np.exp(wave * right) - np.exp(wave * left)) / wave)
    if abs(result.imag) > 1e-8 * max(1., abs(result.real)):
        raise ValueError("interval integral did not preserve a real density")
    return float(result.real)


def metrics(pair, state, rate=None):
    """Both spatial restrictions use the SAME full reconstructed metric."""
    fine, _system = _active(pair, state)
    radial = fine.r * fine.Q
    clocks = _periodic_values(pair.grid, radial, pair.clock_locations)
    parent_length = interval_integral(pair.grid, radial, pair.parent_interval)
    child_length = interval_integral(pair.grid, radial, pair.child_interval)
    annulus_length = parent_length - child_length
    if min(parent_length, child_length, annulus_length) <= 0:
        raise ValueError("full metric has invalid nested proper lengths")
    result = {"r_min": float(np.min(fine.r)), "Q_min": float(np.min(fine.Q)),
              "r_max": float(np.max(fine.r)), "Q_max": float(np.max(fine.Q)),
              "parent_proper_length": parent_length, "child_proper_length": child_length,
              "parent_annulus_proper_length": annulus_length,
              "proper_length_ratio": parent_length / child_length,
              "clock_locations": list(pair.clock_locations), "clock_rates": clocks.tolist(),
              "parent_clock_rate_x1": float(clocks[0]), "child_clock_rate_x2": float(clocks[1]),
              "one_physical_metric": True, "strict_coordinate_inclusion": True,
              "coordinate_scale_ratio": 2.0}
    for name in ("r", "Q"):
        density = radial * getattr(fine, name)
        parent_integral = interval_integral(pair.grid, density, pair.parent_interval)
        child_integral = interval_integral(pair.grid, density, pair.child_interval)
        result[f"child_{name}_proper_mean"] = child_integral / child_length
        result[f"parent_annulus_{name}_proper_mean"] = (parent_integral - child_integral) / annulus_length
    for name in GEOMETRY_NAMES:
        coefficients = getattr(state, name)
        result[f"{name}_child_detail_norm"] = float(np.linalg.norm(coefficients[pair.geometry_child_indices]))
        result[f"{name}_parent_norm"] = float(np.linalg.norm(coefficients[pair.geometry_parent_indices]))
    if rate is not None:
        _check(pair, rate)
        radial_velocity = pair.grid.A_g @ (pair.geometry_map @ rate.r)
        proper_velocity = radial_velocity / radial
        child_velocity = interval_integral(pair.grid, proper_velocity * radial, pair.child_interval)
        parent_velocity = interval_integral(pair.grid, proper_velocity * radial, pair.parent_interval)
        result["child_normal_radial_velocity_mean"] = child_velocity / child_length
        result["parent_annulus_normal_radial_velocity_mean"] = (parent_velocity - child_velocity) / annulus_length
    return result


def functional_pullback(jets, length_scale):
    """Exact inherited functional under x=b+a*xi, T=a*theta, a>0.

    Lhat=a L,Qhat=a Q,rhat=r,chihat=chi,betahat=beta.  The
    first-order density gains a^2; its dtheta dxi measure compensates
    dT dx.  phihat=sqrt(a) phi is the CONTINUUM half-density, not a
    sampled column with an extra sqrt(dx).  pQhat=pQ,prhat=a pr,
    pchihat=a pchi. Physical coefficients and kappa are unchanged.
    Proper clock dtaudtheta=a rL; common-time dtaudT=rL.

    This is coordinate/clock covariance of the same functional. It is not
    physical homothety or equality of the evolving normalized spectra.
    """
    scale = float(length_scale)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("positive finite coordinate scale required")
    powers = {"L": 1, "Q": 1, "L_t": 2, "L_x": 2, "Q_t": 2, "Q_x": 2,
              "L_tt": 3, "L_xx": 3, "Q_tt": 3, "Q_xx": 3,
              "r": 0, "chi": 0, "beta": 0, "r_t": 1, "r_x": 1,
              "chi_t": 1, "chi_x": 1, "beta_t": 1, "beta_x": 1,
              "r_tt": 2, "r_xx": 2, "chi_tt": 2, "chi_xx": 2,
              "p_Q": 0, "p_r": 1, "p_chi": 1,
              "p_Q_x": 1, "p_r_x": 2, "p_chi_x": 2,
              "phi": .5, "phi_t": 1.5, "phi_x": 1.5}
    unknown = set(jets) - set(powers)
    if unknown:
        raise ValueError(f"unsupported affine-pullback jet names: {sorted(unknown)}")
    return {name: np.asarray(value) * scale ** powers[name] for name, value in jets.items()}
