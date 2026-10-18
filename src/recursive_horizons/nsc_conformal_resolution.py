"""Fixed reference resolution of one conformal Cauchy state.

The projector is the canonical dx projector onto a reference eigenspace of
the radial Laplacian. It is not recomputed from later states and it is not
the physical L2(h) projection. Nodal exponentials do not guarantee a
positive prolonged chart.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import eigh

from recursive_horizons.nsc_spherical_coupling import CauchyState


RELATIVE_GATE = 1e-8
_LOG_NAMES = ("a", "b", "chi", "pi_a", "pi_b", "p_chi")


@dataclass(frozen=True)
class ReferenceProjector:
    """Euclidean projector P=U U^T onto the retained reference span."""

    matrix: np.ndarray
    eigenvalues: np.ndarray
    retained: int
    diagnostics: dict
    reference_length: float


def _vector(values, count, name):
    array = np.asarray(values, dtype=float)
    if array.shape != (count,) or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a finite vector of length {count}")
    return array


def _positive_scale(value, name):
    scale = float(value)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return scale


def reference_projector(grid, r_reference, Q_reference, reference_length):
    """Retain mu*L_ref^2<=1 for mu eigenvalues of -Delta_h on the reference slice."""
    length = _positive_scale(reference_length, "reference_length")
    if grid.A_g.shape != (grid.nq, grid.ng) or grid.nq % 2:
        raise ValueError("prolongation must map geometry nodes to an even quadrature grid")
    if np.shape(grid.derivative) != (grid.nq, grid.nq):
        raise ValueError("derivative must be square on the even quadrature grid")
    radius = _vector(r_reference, grid.ng, "r_reference")
    lapse_ratio = _vector(Q_reference, grid.ng, "Q_reference")
    prolonged_radius = grid.A_g @ radius
    prolonged_ratio = grid.A_g @ lapse_ratio
    if np.min(prolonged_radius) <= 0 or np.min(prolonged_ratio) <= 0:
        raise ValueError("prolonged reference radius and Q must be positive")
    slope = grid.derivative @ grid.A_g
    radial_weight = prolonged_radius / prolonged_ratio
    volume_weight = prolonged_radius**3 * prolonged_ratio
    stiffness = grid.dx_q * (slope.T * radial_weight) @ slope
    mass = grid.dx_q * (grid.A_g.T * volume_weight) @ grid.A_g
    stiffness = 0.5 * (stiffness + stiffness.T)
    mass = 0.5 * (mass + mass.T)
    # eigh already checks mass positivity by Cholesky; no second eigensolve.
    try:
        eigenvalues, vectors = eigh(stiffness, mass, check_finite=True)
    except np.linalg.LinAlgError as error:
        raise ValueError("reference mass/eigenproblem is not positive and solvable") from error
    eigenvalues = np.asarray(eigenvalues, dtype=float)
    applied = stiffness @ vectors
    residual = applied - (mass @ vectors) * eigenvalues
    eigen_residual = float(np.linalg.norm(residual) / max(1.0, np.linalg.norm(applied)))
    if eigen_residual > RELATIVE_GATE:
        raise ValueError(f"reference eigen residual {eigen_residual} exceeds {RELATIVE_GATE}")
    scale = float(np.max(np.abs(eigenvalues)))
    floor = 64 * grid.ng * np.finfo(float).eps * scale
    if np.any(eigenvalues < -floor):
        raise ValueError("reference Laplacian has a genuinely negative eigenvalue")
    comparable = eigenvalues.copy()
    # Repair only the known constant kernel. Zeroing every small positive
    # eigenvalue would move physical modes across the cutoff under unit changes.
    constant_overlap = abs(np.sum(vectors,axis=0))/(np.sqrt(grid.ng)*np.linalg.norm(vectors,axis=0))
    constant = int(np.argmax(constant_overlap))
    if abs(eigenvalues[constant]) > floor or constant_overlap[constant] < 1-RELATIVE_GATE:
        raise ValueError("reference eigensystem does not resolve its constant kernel")
    comparable[constant] = 0.0
    if np.any(comparable < 0):
        raise ValueError("nonconstant reference eigenvalue is negative")
    retained_mask = comparable * length**2 <= 1.0
    retained_vectors = vectors[:, retained_mask]
    if retained_vectors.shape[1] == 0:
        matrix = np.zeros((grid.ng, grid.ng))
    else:
        orthogonal, _ = np.linalg.qr(retained_vectors, mode="reduced")
        matrix = orthogonal @ orthogonal.T
    projector_residual = float(np.linalg.norm(matrix @ matrix - matrix) / max(1.0, np.linalg.norm(matrix)))
    symmetry_residual = float(np.linalg.norm(matrix - matrix.T) / max(1.0, np.linalg.norm(matrix)))
    constant_residual = float(np.linalg.norm(matrix@np.ones(grid.ng)-1)/np.sqrt(grid.ng))
    if max(projector_residual, symmetry_residual, constant_residual) > RELATIVE_GATE:
        raise ValueError("reference projector does not resolve its declared span and constant")
    retained_values = eigenvalues[retained_mask]
    excluded_values = eigenvalues[~retained_mask]
    diagnostics = {
        "eigen_residual": eigen_residual,
        "projector_residual": projector_residual,
        "symmetry_residual": symmetry_residual,
        "constant_mode_residual": constant_residual,
        "constant_eigenvalue_before_roundoff_repair": float(eigenvalues[constant]),
        "cutoff_eigenvalue": float(1.0 / length**2),
        "roundoff_floor": float(floor),
        "largest_retained_eigenvalue": None if retained_values.size == 0 else float(np.max(retained_values)),
        "smallest_excluded_eigenvalue": None if excluded_values.size == 0 else float(np.min(excluded_values)),
        "prolonged_radius_min": float(np.min(prolonged_radius)),
        "prolonged_Q_min": float(np.min(prolonged_ratio)),
    }
    return ReferenceProjector(matrix, eigenvalues, int(np.count_nonzero(retained_mask)),
                              diagnostics, length)


def log_state(state):
    """Canonical coordinates a=log r, b=log Q, with chi left unchanged."""
    radius = _vector(state.r, state.r.shape[0], "r")
    ratio = _vector(state.Q, radius.shape[0], "Q")
    if np.min(radius) <= 0 or np.min(ratio) <= 0:
        raise ValueError("log coordinates require positive nodal radius and Q")
    return {
        "a": np.log(radius),
        "b": np.log(ratio),
        "chi": _vector(state.chi, radius.shape[0], "chi"),
        "pi_a": radius * _vector(state.p_r, radius.shape[0], "p_r"),
        "pi_b": ratio * _vector(state.p_Q, radius.shape[0], "p_Q"),
        "p_chi": _vector(state.p_chi, radius.shape[0], "p_chi"),
    }


def from_log_state(values, original_state):
    """Rebuild a Cauchy state. Exponential positivity is only nodal."""
    missing = [name for name in _LOG_NAMES if name not in values]
    if missing:
        raise ValueError("log state is missing " + ",".join(missing))
    count = original_state.r.shape[0]
    coordinate_a = _vector(values["a"], count, "a")
    coordinate_b = _vector(values["b"], count, "b")
    with np.errstate(over="ignore", under="ignore"):
        radius = np.exp(coordinate_a)
        ratio = np.exp(coordinate_b)
    if (not np.isfinite(radius).all() or not np.isfinite(ratio).all()
            or np.min(radius)<=0 or np.min(ratio)<=0):
        raise ValueError("log coordinates must stay finite")
    return CauchyState(
        Q=ratio,
        r=radius,
        chi=_vector(values["chi"], count, "chi"),
        p_Q=_vector(values["pi_b"], count, "pi_b") / ratio,
        p_r=_vector(values["pi_a"], count, "pi_a") / radius,
        p_chi=_vector(values["p_chi"], count, "p_chi"),
        phi0=np.array(original_state.phi0, copy=True),
        phi1=np.array(original_state.phi1, copy=True),
    )


def project_geometry(state, projector_matrix):
    """Apply one fixed projector to every canonical geometry pair."""
    count = state.r.shape[0]
    projector = np.asarray(projector_matrix, dtype=float)
    if projector.shape != (count, count) or not np.isfinite(projector).all():
        raise ValueError("projector must be a finite square geometry matrix")
    logs = log_state(state)
    projected = {name: projector @ logs[name] for name in _LOG_NAMES}
    return from_log_state(projected, state)


def log_rates(state, rate):
    """Push the geometry vector field into the log-canonical coordinates."""
    radius = _vector(state.r, state.r.shape[0], "r")
    ratio = _vector(state.Q, radius.shape[0], "Q")
    if np.min(radius) <= 0 or np.min(ratio) <= 0:
        raise ValueError("rate conversion requires positive nodal radius and Q")
    radius_rate = _vector(rate.r, radius.shape[0], "r_dot")
    ratio_rate = _vector(rate.Q, radius.shape[0], "Q_dot")
    radius_momentum = _vector(state.p_r, radius.shape[0], "p_r")
    ratio_momentum = _vector(state.p_Q, radius.shape[0], "p_Q")
    return {
        "a": radius_rate / radius,
        "b": ratio_rate / ratio,
        "chi": _vector(rate.chi, radius.shape[0], "chi_dot"),
        "pi_a": radius * _vector(rate.p_r, radius.shape[0], "p_r_dot") + radius_momentum * radius_rate,
        "pi_b": ratio * _vector(rate.p_Q, radius.shape[0], "p_Q_dot") + ratio_momentum * ratio_rate,
        "p_chi": _vector(rate.p_chi, radius.shape[0], "p_chi_dot"),
    }


def project_rates(rate_dict, projector_matrix):
    """Apply the same fixed projector to each canonical rate component."""
    missing = [name for name in _LOG_NAMES if name not in rate_dict]
    if missing:
        raise ValueError("canonical rates are missing " + ",".join(missing))
    sample = _vector(rate_dict["a"], np.size(rate_dict["a"]), "a_dot")
    projector = np.asarray(projector_matrix, dtype=float)
    if projector.shape != (sample.shape[0], sample.shape[0]) or not np.isfinite(projector).all():
        raise ValueError("projector must be a finite square geometry matrix")
    return {name: projector @ _vector(rate_dict[name], sample.shape[0], name) for name in _LOG_NAMES}
