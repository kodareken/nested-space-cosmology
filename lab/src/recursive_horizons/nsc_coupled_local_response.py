"""Conditional retained-region consumer of a stored spherical geometry.

The geometry schedule is read from a saved episode. This module does not
regenerate that trajectory and does not launch a campaign driver. The field
operator is the same Fourier–Galerkin Dirac map as ``apply_dirac``: static
lapse and shift, evolved conformal factor ``Q``. The matrix is the
band-limited Fourier multiplication on the fermion prolongation. It is not
the fine-grid identity image.

The retained observer is the initial pair of modes in the first packet,
columns 0 and 1 of the T=0 source, copied with no phase QR. Their weights
stay the original six Gaussian values. Memory, the exterior drive, and the
source weights stay on that preparation. The streamed reducer is called
unchanged. A dense time-indexed exterior propagator is not stored.

``execute_successor`` reads the saved regeneration episode on
``[0.05, 0.085]``. The observer is still the T=0 pair. The initial columns
and the initial ``C_AA``, ``C_EE``, and ``C_AE`` are the transported state
at ``T=0.05``.
"""
from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import LinearOperator, expm_multiply

from .nsc_evolving_reduction import (
    evolve_retained_region,
    hamiltonian_blocks,
    prescribed_six_mode_problem,
)
from .nsc_spherical_coupling import PERIOD, apply_dirac
from .nsc_spherical_galerkin_coupling import (
    build_grid,
    fermion_modes,
    load_physical_columns,
    prolong_geometry,
    pull_geometry,
)

SCHEMA = "NSC-COUPLED-LOCAL-RESPONSE-v1"
BACKEND = "streamed"
PAYLOAD_LIMIT_BYTES = 64 * 1024 * 1024
BUDGET_SAFETY = 1.15
ACTION_GATE = 1e-9
PROLONGATION_GATE = 1e-8
BAND_GATE = 1e-8
_LAB = Path(__file__).resolve().parents[2]
EPISODE_NPZ = _LAB / "results" / "development" / "nsc-spherical-feedback-episode-v1.npz"
DEFAULT_JSON = _LAB / "results" / "development" / "nsc-coupled-local-response-v1.json"
DEFAULT_NPZ = _LAB / "results" / "development" / "nsc-coupled-local-response-v1.npz"

PARTNER_CASE = {
    "nf512_dt_0_0005": "nf512_dt_0_00025",
    "nf512_dt_0_00025": "nf512_dt_0_0005",
    "nf256_dt_0_0005": "nf256_dt_0_00025",
    "nf256_dt_0_00025": "nf256_dt_0_0005",
}

INTERFACE_GAPS = (
    "evolve_retained_region requires a dense Hermitian callback. A LinearOperator is not a supported H(t), so each sample still forms the dense exterior block after the full matrix.",
    "backend='streamed' completes one dense null_space frame and replays the interaction history from the initial time at every new output node. That replay grows quadratically with the number of output nodes and is recomputation, not a stored W.",
    "Passing a full covariance evolves one column per Hilbert-space dimension and stores history with shape (time, N_outside, N). The rank-6 source is evolved with its column weights so that history stays (time, N_outside, 6).",
)

STRESS_GAP = (
    "No effective-action variation is derived. Occupied and empty kernels are not converted into a stress. "
    "A force comparison would need a kernel metric variation; that variation is not computed."
)

NEXT_CLAIM = (
    "On the full stored frame window, with the exterior midpoint refined by substeps=2, "
    "the streamed region-0 occupation and coherence agree with an independent full-band evolution "
    "of the same initial columns and weights to the measured truncation, at both saved resolutions. "
    "The geometry schedule stays the stored frames. The run does not regenerate a trajectory and does not certify a stress."
)


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def jsonable(value):
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": jsonable(np.real(value)), "imag": jsonable(np.imag(value))}
        return jsonable(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if value is None or isinstance(value, str):
        return value
    return str(value)


def write_json(path, record):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(jsonable(record), indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def write_npz(path, arrays):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    temporary.replace(path)


class MeasuredBudget:
    """Process-time budget. An estimate is admitted only with the safety factor."""

    def __init__(self, limit_s, safety=BUDGET_SAFETY):
        self.limit_s = float(limit_s)
        self.safety = float(safety)
        self.started = time.process_time()

    def elapsed(self):
        return float(time.process_time() - self.started)

    def remaining(self):
        return float(self.limit_s - self.elapsed())

    def allows(self, estimate_s):
        estimate_s = float(estimate_s)
        if estimate_s < 0 or not np.isfinite(estimate_s):
            raise ValueError("estimate must be finite and non-negative")
        return self.elapsed() + self.safety * estimate_s < self.limit_s


def quadratic_estimate(probe_seconds, probe_steps, steps, substeps=1, probe_substeps=1):
    """Extrapolate the streamed replay as probe time times (steps/probe steps)^2 times substeps.

    The probe time already includes one frame completion and the Hamiltonian
    builds on that probe. Using the whole probe as the quadratic coefficient
    overestimates the fixed part, which is the conservative budget check.
    """
    if probe_steps < 1 or steps < 1 or substeps < 1 or probe_substeps < 1:
        raise ValueError("step counts must be positive")
    if probe_seconds < 0:
        raise ValueError("probe time must be non-negative")
    ratio = float(steps) / float(probe_steps)
    refinement = float(substeps) / float(probe_substeps)
    return float(probe_seconds) * ratio * ratio * refinement


def plan_remaining(remaining_s, probe_seconds, probe_steps, available_steps, safety=BUDGET_SAFETY):
    """Choose a prefix that still fits after the probe has already been paid."""
    remaining_s = float(remaining_s)
    refine_steps = int(probe_steps) * 2
    refine_cost = quadratic_estimate(probe_seconds, probe_steps, refine_steps)
    omission_cost = 3.0 * float(probe_seconds)
    full_cost = float(probe_seconds)
    overhead = refine_cost + omission_cost + full_cost

    def pack(main_steps, refine, omissions, independent):
        estimate = 0.0
        if main_steps != probe_steps:
            estimate += quadratic_estimate(probe_seconds, probe_steps, main_steps)
        if refine:
            estimate += refine_cost
        if omissions:
            estimate += omission_cost
        if independent and main_steps == probe_steps:
            estimate += full_cost
        elif independent and main_steps != probe_steps:
            estimate += full_cost * (float(main_steps) / float(probe_steps))
        return {
            "main_steps": int(main_steps),
            "refine": bool(refine),
            "omissions": bool(omissions),
            "independent_full": bool(independent),
            "estimate_s": float(estimate),
        }

    if available_steps < probe_steps:
        raise ValueError("available window is shorter than the probe")
    for steps in range(int(available_steps), int(probe_steps), -1):
        proposal = pack(steps, True, True, True)
        if safety * proposal["estimate_s"] < remaining_s:
            return proposal
    for refine, omissions, independent in (
        (True, True, True),
        (True, False, True),
        (False, True, True),
        (False, False, True),
        (False, False, False),
    ):
        proposal = pack(probe_steps, refine, omissions, independent)
        if safety * proposal["estimate_s"] < remaining_s or proposal["estimate_s"] == 0.0:
            return proposal
    return pack(probe_steps, False, False, False)


def _resolution_prefix(case):
    if case not in PARTNER_CASE:
        raise ValueError(f"unknown episode case {case}")
    return "nf512" if case.startswith("nf512") else "nf256"


def load_episode_case(case, path=EPISODE_NPZ):
    """Load stored frames, the initial source, and bindings. Does not evolve."""
    case = str(case)
    resolution = _resolution_prefix(case)
    partner = PARTNER_CASE[case]
    with np.load(path, allow_pickle=False) as data:
        frame_time = np.array(data[case + "_frame_time"], dtype=float, copy=True)
        frame_q = np.array(data[case + "_frame_Q"], dtype=float, copy=True)
        partner_q = np.array(data[partner + "_frame_Q"], dtype=float, copy=True)
        phi0 = np.array(data[resolution + "_initial_phi0"], copy=True)
        phi1 = np.array(data[resolution + "_initial_phi1"], copy=True)
        weights = np.array(data[resolution + "_occupations"], dtype=float, copy=True)
        names = np.array(data[resolution + "_binding_names"])
        values = np.array(data[resolution + "_binding_values"], dtype=float, copy=True)
        coarse_q = np.array(data[resolution + "_initial_Q"], dtype=float, copy=True)
        final_phi0 = np.array(data[case + "_final_phi0"], copy=True)
        final_phi1 = np.array(data[case + "_final_phi1"], copy=True)
        final_phi0_dot = np.array(data[case + "_final_rate_phi0_dot"], copy=True)
        final_phi1_dot = np.array(data[case + "_final_rate_phi1_dot"], copy=True)
        final_q_dot = np.array(data[case + "_final_rate_Q_dot"], dtype=float, copy=True)
        verdict = str(data["verdict"])
    if frame_q.shape[0] != frame_time.shape[0]:
        raise ValueError("stored frame times and Q samples differ in length")
    if partner_q.shape != frame_q.shape:
        raise ValueError("partner geometry does not share the stored frame shape")
    if np.min(frame_q) <= 0:
        raise ValueError("stored conformal factor left the positive chart")
    return {
        "case": case,
        "partner_case": partner,
        "resolution": resolution,
        "frame_time": frame_time,
        "frame_Q": frame_q,
        "phi0": phi0,
        "phi1": phi1,
        "weights": weights,
        "binding_names": [str(item) for item in names.tolist()],
        "binding_values": values,
        "coarse_Q": coarse_q,
        "final_phi0": final_phi0,
        "final_phi1": final_phi1,
        "final_phi0_dot": final_phi0_dot,
        "final_phi1_dot": final_phi1_dot,
        "final_Q_dot": final_q_dot,
        "episode_verdict_label": verdict,
        "episode_verdict_used_as_acceptance": False,
        "time_resolution_geometry_gap": float(np.max(np.abs(frame_q - partner_q))),
        "interpolation": interpolation_scale(frame_time, frame_q),
        "geometry_rate_series_stored": False,
        "final_coarse_Q_dot_norm": float(np.linalg.norm(final_q_dot)),
        "dirac_geometry_argument": "stored fine Q",
        "stored_but_not_dirac_matrix_entries": ["r", "chi", "p_Q", "p_r", "p_chi"],
        "static_gauge": ["length_density", "shift"],
    }


def interpolation_scale(times, values):
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    steps = np.diff(times)
    if times.ndim != 1 or times.size < 3 or np.any(steps <= 0):
        raise ValueError("interpolation scale needs at least three increasing times")
    uniform = bool(np.allclose(steps, steps[0], rtol=0, atol=1e-12 * max(1.0, abs(float(steps[0])))))
    second = values[2:] - 2.0 * values[1:-1] + values[:-2]
    maximum = float(np.max(np.abs(second)))
    return {
        "uniform_stored_frames": uniform,
        "stored_frame_step": float(steps[0]),
        "max_abs_second_difference": maximum,
        "midpoint_linear_scale": maximum / 8.0,
        "meaning": "max |second difference| / 8 on the stored frames; linear midpoint scale for a uniform quadratic",
    }


def interpolate_rows(times, values, query):
    """Piecewise-linear samples. Exact nodes return the stored row."""
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    query = np.atleast_1d(np.asarray(query, dtype=float))
    if values.shape[0] != times.shape[0]:
        raise ValueError("value rows and times differ")
    if np.any(query < times[0] - 1e-12) or np.any(query > times[-1] + 1e-12):
        raise ValueError("query leaves the stored geometry window")
    index = np.searchsorted(times, query, side="right") - 1
    index = np.clip(index, 0, times.size - 2)
    left_time = times[index]
    right_time = times[index + 1]
    width = right_time - left_time
    fraction = (query - left_time) / width
    fraction = np.clip(fraction, 0.0, 1.0)
    left = values[index]
    right = values[index + 1]
    return (1.0 - fraction)[:, None] * left + fraction[:, None] * right


FINISHED_STATUSES = (
    "MEASURED_CONDITIONAL_REDUCTION",
    "BUDGET_EXCEEDED",
)
DECLARED_INTERPOLATION = "piecewise-linear-stored-Q"
CANONICAL_SETTINGS = {
    "backend": BACKEND,
    "interpolation": DECLARED_INTERPOLATION,
    "memory": True,
    "outside_drive": True,
    "coupling": True,
    "drop_cross_covariance": False,
}


def frame_slopes(times, values):
    """Finite-difference slopes of stored rows. These are not saved time derivatives."""
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    slopes = np.empty_like(values)
    steps = np.diff(times)
    slopes[0] = (values[1] - values[0]) / steps[0]
    slopes[-1] = (values[-1] - values[-2]) / steps[-1]
    if times.size > 2:
        span = (times[2:] - times[:-2])[:, None]
        slopes[1:-1] = (values[2:] - values[:-2]) / span
    return slopes


def hermite_rows(times, values, query):
    """Cubic Hermite samples. Stored nodes are returned exactly."""
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    slopes = frame_slopes(times, values)
    query = np.atleast_1d(np.asarray(query, dtype=float))
    if np.any(query < times[0] - 1e-12) or np.any(query > times[-1] + 1e-12):
        raise ValueError("query leaves the stored geometry window")
    index = np.searchsorted(times, query, side="right") - 1
    index = np.clip(index, 0, times.size - 2)
    width = times[index + 1] - times[index]
    fraction = np.clip((query - times[index]) / width, 0.0, 1.0)
    cube = fraction ** 3
    square = fraction ** 2
    h00 = 2.0 * cube - 3.0 * square + 1.0
    h10 = cube - 2.0 * square + fraction
    h01 = -2.0 * cube + 3.0 * square
    h11 = cube - square
    return (
        h00[:, None] * values[index]
        + (h10 * width)[:, None] * slopes[index]
        + h01[:, None] * values[index + 1]
        + (h11 * width)[:, None] * slopes[index + 1]
    )


def interpolation_variant_report(times, values):
    """Compare prescribed interpolants. Saved nodal rates are not a frame series."""
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    midpoints = 0.5 * (times[:-1] + times[1:])
    linear = interpolate_rows(times, values, midpoints)
    hermite = hermite_rows(times, values, midpoints)
    node_gap = float(np.max(np.abs(hermite_rows(times, values, times) - values)))
    return {
        "saved_nodal_Q_rate_series": False,
        "declared_geometry": DECLARED_INTERPOLATION,
        "hermite_slopes": "finite differences of stored Q frames, not saved Q_dot",
        "rate_gap": (
            "Q_dot is stored only as the final coarse sample. "
            "A Hermite curve through the frames is not determined by saved rates."
        ),
        "midpoint_max_abs_linear_versus_hermite": float(np.max(np.abs(linear - hermite))),
        "hermite_node_gap": node_gap,
    }


def saved_rate_versus_frame_secant(grid, frame_time, frame_q, coarse_rate):
    """Compare the only saved coarse Q rate with the last stored-frame secant."""
    step = float(frame_time[-1] - frame_time[-2])
    secant = (frame_q[-1] - frame_q[-2]) / step
    pulled = pull_geometry(grid, secant)
    rate = np.asarray(coarse_rate, dtype=float)
    if pulled.shape != rate.shape:
        raise ValueError("pulled frame secant and saved coarse rate differ in shape")
    absolute = float(np.max(np.abs(pulled - rate)))
    scale = float(np.max(np.abs(rate)))
    return {
        "saved_rate": "final coarse Q_dot only",
        "max_abs": absolute,
        "saved_rate_max_abs": scale,
        "relative_to_saved_rate": None if scale == 0.0 else absolute / scale,
    }


def source_id(phi0, phi1, weights):
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(phi0).view(np.uint8))
    digest.update(np.ascontiguousarray(phi1).view(np.uint8))
    digest.update(np.ascontiguousarray(np.asarray(weights, dtype=np.float64)).tobytes())
    return digest.hexdigest()


def request_identity(domain, case, substeps, source, settings=None):
    chosen = dict(CANONICAL_SETTINGS if settings is None else settings)
    return {
        "domain": str(domain),
        "case": str(case),
        "substeps": int(substeps),
        "source": str(source),
        "backend": str(chosen["backend"]),
        "interpolation": str(chosen["interpolation"]),
        "memory": bool(chosen["memory"]),
        "outside_drive": bool(chosen["outside_drive"]),
        "coupling": bool(chosen["coupling"]),
        "drop_cross_covariance": bool(chosen["drop_cross_covariance"]),
    }


def _legacy_pilot_identity(previous, source):
    return request_identity(
        previous.get("requested_domain"),
        previous.get("requested_case"),
        int(previous.get("requested_substeps", -1)),
        source,
    )


def resume_disposition(previous, identity):
    """Return, reject, or continue a finished record before any recomputation.

    A finished record is returned only when domain, substeps, case, source,
    and the declared settings all match and that request's phase is complete.
    """
    if not isinstance(previous, dict):
        return "continue"
    status = previous.get("status")
    if status not in FINISHED_STATUSES:
        return "continue"
    stored = previous.get("requests", {}).get(identity["domain"])
    if stored is None and identity["domain"] == "pilot":
        stored = _legacy_pilot_identity(previous, identity["source"])
        if (
            previous.get("requested_domain") != identity["domain"]
            or previous.get("requested_case") != identity["case"]
            or int(previous.get("requested_substeps", -1)) != int(identity["substeps"])
        ):
            return "reject"
        if stored != identity:
            return "reject"
        if previous.get("headline", {}).get("complete") is True:
            return "return"
        return "reject"
    if stored is None:
        return "reject"
    if stored != identity:
        return "reject"
    phase_name = "full_window" if identity["domain"] == "stored-frames" else "comparison"
    phase = previous.get("phases", {}).get(phase_name, {})
    if phase.get("complete") is True and int(phase.get("substeps", identity["substeps"])) == int(identity["substeps"]):
        return "return"
    if identity["domain"] == "pilot" and previous.get("headline", {}).get("complete") is True:
        return "return"
    return "continue"


def spinor_local_moments(observer, phi0, phi1, weights):
    stacked = np.concatenate((np.asarray(phi0), np.asarray(phi1)), axis=0)
    projected = np.asarray(observer).conj().T @ stacked
    occupation, coherence = occupations_and_coherence(projected[None, ...], weights)
    return {
        "occupation": [float(item) for item in occupation[0]],
        "coherence": complex(coherence[0]),
    }


def correlated_cross_separation():
    """Active cross omission on the prescribed six-mode state, not the rank-6 source."""
    problem = prescribed_six_mode_problem()
    times = np.linspace(0.0, 1.0, 17)
    kept = evolve_retained_region(
        problem["hamiltonian"],
        problem["local_basis"],
        times,
        covariance=problem["covariance"],
        backend=BACKEND,
    )
    dropped = evolve_retained_region(
        problem["hamiltonian"],
        problem["local_basis"],
        times,
        covariance=problem["covariance"],
        drop_cross_covariance=True,
        backend=BACKEND,
    )
    transformed = kept["frame"].conj().T @ problem["covariance"] @ kept["frame"]
    cross = transformed[:2, 2:]
    separation = float(np.linalg.norm(kept["covariance_total"] - dropped["covariance_total"]))
    return {
        "control": "prescribed_six_mode_problem",
        "cross_norm": float(np.linalg.norm(cross, ord="fro")),
        "covariance_separation": separation,
        "active": bool(separation > 1e-6 and np.linalg.norm(cross) > 1e-6),
        "stores_dense_exterior_propagator": bool(kept["flags"]["stores_dense_exterior_propagator"]),
    }


def multiplication_matrix(factors, values):
    """Galerkin multiplication in the nodal fermion basis, via the Fourier band.

    ``values`` is a fine nodal profile. The product is formed on the fermion
    Fourier columns and rotated by the prolongation change of basis. No
    fine-grid identity and no ``nq × nq`` operator is allocated.
    """
    profile = np.asarray(values, dtype=float)
    fourier = factors["fourier_band"]
    if profile.shape != (fourier.shape[0],):
        raise ValueError("multiplication profile must be sampled on the quadrature grid")
    band = fourier.conj().T @ (profile[:, None] * fourier)
    change = factors["nodal_from_fourier"]
    matrix = change.conj().T @ band @ change
    return 0.5 * (matrix + matrix.conj().T)


def anticommutator_matrix(factors, values):
    multiplied = multiplication_matrix(factors, values)
    momentum = factors["momentum_band"]
    block = 0.5 * (multiplied @ momentum + momentum @ multiplied)
    return 0.5 * (block + block.conj().T)


def fourier_galerkin_factors(grid):
    """Fourier factors for ``U† H_fine U`` on the antiperiodic fermion band."""
    count_f = int(grid.nf)
    count_q = int(grid.nq)
    length = float(grid.length)
    indices = np.arange(count_q) - (count_q // 2)
    fine_keys = np.rint((indices.astype(float) + 0.5) * 2.0).astype(int)
    band_keys = np.rint(fermion_modes(count_f) * 2.0).astype(int)
    position = {int(key): index for index, key in enumerate(fine_keys)}
    try:
        selected = [position[int(key)] for key in band_keys]
    except KeyError as error:
        raise RuntimeError("fermion band is not a subset of the quadrature momenta") from error
    coordinate = indices.astype(float) * (length / count_q)
    momenta = 2.0 * np.pi * (indices[selected].astype(float) + 0.5) / length
    fourier = np.exp(1j * coordinate[:, None] * momenta[None, :]) / np.sqrt(count_q)
    change = fourier.conj().T @ grid.U_f
    reconstructed = fourier @ change
    residual = float(np.linalg.norm(reconstructed - grid.U_f) / np.linalg.norm(grid.U_f))
    if residual > BAND_GATE:
        raise RuntimeError(f"Fourier band does not reproduce the column prolongation: {residual}")
    factors = {
        "nf": count_f,
        "nq": count_q,
        "length": length,
        "kappa": int(grid.fine.kappa),
        "fourier_band": fourier,
        "nodal_from_fourier": change,
        "momenta": np.array(momenta, dtype=float, copy=True),
        "length_density": np.array(grid.fine.length_density, dtype=float, copy=True),
        "shift": np.array(grid.fine.shift, dtype=float, copy=True),
        "band_reconstruction": residual,
        "fine_operator_image": False,
        "fine_identity_image": False,
    }
    momentum = change.conj().T @ (momenta[:, None] * change)
    factors["momentum_band"] = 0.5 * (momentum + momentum.conj().T)
    factors["mass"] = multiplication_matrix(factors, factors["kappa"] * factors["length_density"])
    factors["shift_block"] = anticommutator_matrix(factors, factors["shift"])
    return factors


def galerkin_matrix(factors, radial_density):
    """Dense Dirac matrix in the stacked coarse spinor basis ``[phi0; phi1]``."""
    density = np.asarray(radial_density, dtype=float)
    if density.shape != (factors["nq"],) or not np.isfinite(density).all() or float(np.min(density)) <= 0.0:
        raise ValueError("quadrature radial density must be finite and positive")
    kinetic = anticommutator_matrix(factors, factors["length_density"] / density)
    shift = factors["shift_block"]
    mass = factors["mass"]
    count = int(factors["nf"])
    matrix = np.empty((2 * count, 2 * count), dtype=np.complex128)
    matrix[:count, :count] = -shift
    matrix[:count, count:] = -1j * kinetic + mass
    matrix[count:, :count] = 1j * kinetic + mass
    matrix[count:, count:] = -shift
    defect = float(np.linalg.norm(matrix - matrix.conj().T, ord="fro"))
    scale = max(1.0, float(np.linalg.norm(matrix, ord="fro")))
    if defect > 1e-8 * scale:
        raise RuntimeError(f"Galerkin Dirac matrix is not Hermitian: {defect}")
    return 0.5 * (matrix + matrix.conj().T)


def dirac_pullback(grid, radial_density, phi0, phi1):
    fine0 = grid.U_f @ np.asarray(phi0)
    fine1 = grid.U_f @ np.asarray(phi1)
    image0, image1 = apply_dirac(
        fine0,
        fine1,
        grid.fine.length_density,
        np.asarray(radial_density, dtype=float),
        grid.fine.shift,
        grid.fine.kappa,
        grid.fine.momentum,
    )
    return grid.U_f.conj().T @ image0, grid.U_f.conj().T @ image1


def action_residual(grid, factors, radial_density, phi0, phi1):
    left0, left1 = dirac_pullback(grid, radial_density, phi0, phi1)
    stacked = galerkin_matrix(factors, radial_density) @ np.concatenate(
        (np.asarray(phi0), np.asarray(phi1)),
        axis=0,
    )
    count = int(grid.nf)
    absolute = float(np.linalg.norm(stacked[:count] - left0) + np.linalg.norm(stacked[count:] - left1))
    scale = float(np.linalg.norm(left0) + np.linalg.norm(left1))
    relative = 0.0 if scale == 0.0 else absolute / scale
    return {"relative": relative, "absolute": absolute, "scale": scale}


def region_observer(phi0, phi1, region=0):
    """Fixed columns of one packet. Region 0 is the first two initial modes."""
    phi0 = np.asarray(phi0)
    phi1 = np.asarray(phi1)
    if phi0.shape != phi1.shape or phi0.ndim != 2:
        raise ValueError("spinor columns must share one matrix shape")
    start = 2 * int(region)
    if start < 0 or start + 2 > phi0.shape[1]:
        raise ValueError("region pair is outside the source columns")
    stacked = np.concatenate((phi0, phi1), axis=0).astype(np.complex128, copy=False)
    observer = np.array(stacked[:, start : start + 2], dtype=np.complex128, copy=True)
    columns = np.array(stacked, dtype=np.complex128, copy=True)
    gram = phi0.conj().T @ phi0 + phi1.conj().T @ phi1
    return observer, columns, gram


def source_covariance(columns, weights):
    scale = np.asarray(weights, dtype=float)
    covariance = (np.asarray(columns) * scale) @ np.asarray(columns).conj().T
    return 0.5 * (covariance + covariance.conj().T)


def frame_source_blocks(frame, columns, weights, retained=2):
    coefficients = np.asarray(frame).conj().T @ np.asarray(columns)
    weighted = coefficients * np.asarray(weights, dtype=float)
    parent = weighted[:retained] @ coefficients[:retained].conj().T
    cross = weighted[:retained] @ coefficients[retained:].conj().T
    return {
        "parent": parent,
        "cross": cross,
        "parent_norm": float(np.linalg.norm(parent, ord="fro")),
        "cross_norm": float(np.linalg.norm(cross, ord="fro")),
    }


def independent_radial_density(grid, time):
    """Positive profile for a method control. Not a stored episode frame."""
    coordinate = grid.xi_q
    wave = np.sin(2.0 * np.pi * coordinate / grid.length)
    harmonic = np.cos(4.0 * np.pi * coordinate / grid.length)
    return 0.24 * (1.0 + 0.2 * wave * np.sin(2.0 * float(time)) + 0.05 * harmonic)


def hamiltonian_from_radial(factors, radial):
    """Cached dense ``H(t)``. The cache key is the exact float hex of ``t``."""
    cache = {}
    stats = {"builds": 0, "hits": 0, "matrix_bytes": 0}

    def hamiltonian(time):
        key = float(time).hex()
        cached = cache.get(key)
        if cached is not None:
            stats["hits"] += 1
            return cached
        matrix = galerkin_matrix(factors, radial(float(time)))
        cache[key] = matrix
        stats["builds"] += 1
        stats["matrix_bytes"] = int(matrix.nbytes)
        return matrix

    hamiltonian.cache = cache
    hamiltonian.stats = stats
    return hamiltonian


def independent_control(nf=24):
    """Real rank-6 packet on a small band, with an independent geometry pulse."""
    count = int(nf)
    grid = build_grid(count, quadrature=4 * count, length=PERIOD)
    phi0, phi1, preparation, problems = load_physical_columns(count)
    if problems:
        raise RuntimeError("rank-6 preparation was not accepted: " + "; ".join(problems))
    factors = fourier_galerkin_factors(grid)
    observer, columns, gram = region_observer(phi0, phi1, region=0)
    weights = np.array(grid.fine.occupations, dtype=float, copy=True)
    radial = lambda time: independent_radial_density(grid, time)
    return {
        "grid": grid,
        "factors": factors,
        "phi0": phi0,
        "phi1": phi1,
        "preparation": preparation,
        "observer": observer,
        "columns": columns,
        "gram": gram,
        "weights": weights,
        "hamiltonian": hamiltonian_from_radial(factors, radial),
        "radial": radial,
        "regeneration": False,
        "stored_episode": False,
    }


def stored_radial(case):
    def radial(time, times=case["frame_time"], values=case["frame_Q"]):
        return interpolate_rows(times, values, [float(time)])[0]

    return radial


def signal_error(approximate, reference):
    difference = np.asarray(approximate) - np.asarray(reference)
    absolute = float(np.max(np.abs(difference)))
    scale = float(np.max(np.abs(reference)))
    change = np.asarray(reference) - np.asarray(reference)[0]
    change_scale = float(np.max(np.abs(change)))
    return {
        "max_abs": absolute,
        "signal_max_abs": scale,
        "relative_to_signal": None if scale == 0.0 else absolute / scale,
        "change_from_initial_max_abs": change_scale,
        "relative_to_change": None if change_scale == 0.0 else absolute / change_scale,
    }


def phase_error(approximate, reference):
    reference = np.asarray(reference)
    approximate = np.asarray(approximate)
    peak = float(np.max(np.abs(reference))) if reference.size else 0.0
    if peak == 0.0:
        return {"defined": False, "max_abs_radians": None, "reason": "reference coherence is zero"}
    mask = np.abs(reference) > 1e-8 * peak
    if not np.any(mask):
        return {"defined": False, "max_abs_radians": None, "reason": "reference coherence stays below 1e-8 of its peak"}
    delta = np.angle(approximate[mask] * np.conjugate(reference[mask]))
    return {
        "defined": True,
        "max_abs_radians": float(np.max(np.abs(delta))),
        "samples": int(np.count_nonzero(mask)),
        "reference_peak_abs": peak,
    }


def occupations_and_coherence(projected, weights):
    scale = np.asarray(weights, dtype=float)
    occupation = np.sum((np.abs(projected) ** 2) * scale[None, None, :], axis=2)
    coherence = np.sum(projected[:, 0, :] * np.conjugate(projected[:, 1, :]) * scale[None, :], axis=1)
    return occupation, coherence


def evolve_full_columns(hamiltonian, times, columns, substeps=1):
    """Midpoint exponential of the full ``H(t)`` on columns. No propagator is stored."""
    state = np.array(columns, dtype=np.complex128, copy=True)
    grid = np.asarray(times, dtype=float)
    if state.ndim != 2 or grid.ndim != 1 or grid.size < 2:
        raise ValueError("full evolution needs a time grid and a column matrix")
    history = np.empty((grid.size,) + state.shape, dtype=np.complex128)
    history[0] = state
    refinement = int(substeps)
    if refinement < 1:
        raise ValueError("substeps must be positive")
    for index in range(grid.size - 1):
        width = float(grid[index + 1] - grid[index])
        sample = width / refinement
        for part in range(refinement):
            midpoint = float(grid[index] + (part + 0.5) * sample)
            generator = np.asarray(hamiltonian(midpoint), dtype=np.complex128)
            coefficient = -1j * sample
            state = expm_multiply(
                coefficient * generator,
                state,
                traceA=np.trace(coefficient * generator),
            )
        history[index + 1] = state
    return history


def allocation_report(result):
    allocation = result["allocation"]
    history = result["history"]
    exterior = int(allocation["exterior_dimension"])
    count = int(allocation["output_times"])
    return {
        "backend": allocation["backend"],
        "time_indexed_exterior_propagator_bytes": int(allocation["time_indexed_exterior_propagator_bytes"]),
        "stores_dense_exterior_propagator": bool(result["flags"]["stores_dense_exterior_propagator"]),
        "history_shape": [int(item) for item in history.shape],
        "history_bytes": int(history.nbytes),
        "frame_bytes": int(allocation["frame_bytes"]),
        "frame_completion": allocation["frame_completion"],
        "current_exterior_block_bytes": int(allocation["current_exterior_block_bytes"]),
        "peak_column_workspace_bytes": int(allocation["peak_column_workspace_bytes"]),
        "avoided_dense_propagator_bytes": int(count * exterior * exterior * 16),
        "propagator_storage": allocation["propagator_storage"],
        "output_times": count,
        "exterior_dimension": exterior,
        "retained_dimension": int(allocation.get("retained_dimension", result["local_basis"].shape[1])),
    }


def evolve_streamed(hamiltonian, observer, times, columns, weights, **flags):
    result = evolve_retained_region(
        hamiltonian,
        observer,
        times,
        columns,
        weights=np.array(weights, dtype=float, copy=True),
        backend=BACKEND,
        **flags,
    )
    if not np.array_equal(result["local_basis"], observer):
        raise RuntimeError("fixed observer was not preserved exactly")
    if int(result["allocation"]["time_indexed_exterior_propagator_bytes"]) != 0:
        raise RuntimeError("streamed result stored a time-indexed exterior propagator")
    if result["history"].shape[-1] == result["history"].shape[-2]:
        raise RuntimeError("history has the shape of a dense exterior propagator")
    if flags.get("coupling", True):
        _retained, coupling, _exterior = hamiltonian_blocks(
            hamiltonian(float(np.asarray(times, dtype=float)[0])),
            result["frame"],
            int(observer.shape[1]),
        )
        if not np.allclose(result["coupling_block"][0], coupling, rtol=1e-8, atol=1e-8):
            raise RuntimeError("retained coupling is not the Galerkin H_AE block")
    occupation, coherence = occupations_and_coherence(result["amplitudes"], weights)
    return result, occupation, coherence


def compare_retained(hamiltonian, observer, times, columns, weights, substeps=1, **flags):
    result, occupation, coherence = evolve_streamed(
        hamiltonian,
        observer,
        times,
        columns,
        weights,
        propagator_substeps=int(substeps),
        **flags,
    )
    full = evolve_full_columns(hamiltonian, times, columns, substeps=int(substeps))
    projected = np.einsum("ij,tjk->tik", observer.conj().T, full)
    full_occupation, full_coherence = occupations_and_coherence(projected, weights)
    return {
        "result": result,
        "occupation_retained": occupation,
        "occupation_full": full_occupation,
        "coherence_retained": coherence,
        "coherence_full": full_coherence,
        "occupation_error": signal_error(occupation, full_occupation),
        "coherence_error": signal_error(coherence, full_coherence),
        "phase_error": phase_error(coherence, full_coherence),
        "allocation": allocation_report(result),
        "trapezoid_residual_max": float(result["trapezoid_residual_max"]),
        "history_norm": float(np.linalg.norm(result["history"])),
        "initial_exterior_norm": float(np.linalg.norm(result["initial_exterior"])),
        "observer_preserved": True,
    }


def project_full_occupation(hamiltonian, observer, times, columns, weights, substeps=1):
    full = evolve_full_columns(hamiltonian, times, columns, substeps=substeps)
    projected = np.einsum("ij,tjk->tik", observer.conj().T, full)
    occupation, coherence = occupations_and_coherence(projected, weights)
    return occupation, coherence


def commutator_norm(left, right):
    return float(np.linalg.norm(left @ right - right @ left, ord="fro"))


def verify_small_control(nf=24, seed=5):
    control = independent_control(nf)
    grid = control["grid"]
    factors = control["factors"]
    radial0 = control["radial"](0.0)
    radial1 = control["radial"](0.35)
    rng = np.random.default_rng(seed)
    count = int(grid.nf)
    probe0 = rng.normal(size=(count, 3)) + 1j * rng.normal(size=(count, 3))
    probe1 = rng.normal(size=(count, 3)) + 1j * rng.normal(size=(count, 3))
    residual0 = action_residual(grid, factors, radial0, probe0, probe1)
    residual1 = action_residual(grid, factors, radial1, control["phi0"], control["phi1"])
    leaked = grid.fine.momentum @ grid.U_f
    band_momentum = grid.U_f.conj().T @ leaked
    invariance = float(np.linalg.norm(leaked - grid.U_f @ band_momentum) / np.linalg.norm(leaked))
    matrix0 = galerkin_matrix(factors, radial0)
    matrix1 = galerkin_matrix(factors, radial1)
    covariance = source_covariance(control["columns"], control["weights"])
    if residual0["relative"] > ACTION_GATE or residual1["relative"] > ACTION_GATE:
        raise RuntimeError(f"Fourier Galerkin action failed the apply_dirac gate: {residual0} {residual1}")
    if invariance > 1e-10:
        raise RuntimeError(f"fine momentum does not preserve the fermion band: {invariance}")
    return {
        "complete": True,
        "nf": count,
        "nq": int(grid.nq),
        "action_relative_random": residual0["relative"],
        "action_relative_packet": residual1["relative"],
        "band_invariance": invariance,
        "band_reconstruction": factors["band_reconstruction"],
        "fine_identity_image": False,
        "matrix_shape": [int(matrix0.shape[0]), int(matrix0.shape[1])],
        "commutator_norm": commutator_norm(matrix0, matrix1),
        "covariance_commutator_norm": commutator_norm(matrix0, covariance),
        "fine_operator_image": False,
    }


def verify_stored_operator(grid, factors, case):
    packet = action_residual(grid, factors, case["frame_Q"][0], case["phi0"], case["phi1"])
    final = action_residual(grid, factors, case["frame_Q"][-1], case["final_phi0"], case["final_phi1"])
    pull0, pull1 = dirac_pullback(grid, case["frame_Q"][-1], case["final_phi0"], case["final_phi1"])
    rate_scale = float(np.linalg.norm(case["final_phi0_dot"]) + np.linalg.norm(case["final_phi1_dot"]))
    rate_absolute = float(
        np.linalg.norm(-1j * pull0 - case["final_phi0_dot"])
        + np.linalg.norm(-1j * pull1 - case["final_phi1_dot"])
    )
    prolonged = prolong_geometry(grid, case["coarse_Q"])
    prolong_gap = float(np.max(np.abs(prolonged - case["frame_Q"][0])))
    if packet["relative"] > ACTION_GATE or final["relative"] > ACTION_GATE:
        raise RuntimeError("stored-geometry Galerkin action failed the apply_dirac gate")
    if rate_scale == 0.0 or rate_absolute / rate_scale > 1e-8:
        raise RuntimeError("stored final spinor rate does not match apply_dirac")
    if prolong_gap > PROLONGATION_GATE:
        raise RuntimeError("stored fine Q does not match the prolonged coarse geometry")
    return {
        "complete": True,
        "initial_action_relative": packet["relative"],
        "final_action_relative": final["relative"],
        "final_rate_relative": rate_absolute / rate_scale,
        "initial_prolongation_gap": prolong_gap,
        "fine_identity_image": False,
    }


def build_case_grid(case):
    count_f = int(case["phi0"].shape[0])
    count_q = int(case["frame_Q"].shape[1])
    if case["coarse_Q"].shape != (count_f - 1,):
        raise ValueError("stored coarse geometry does not match the fermion band")
    grid = build_grid(count_f, quadrature=count_q, length=PERIOD)
    factors = fourier_galerkin_factors(grid)
    return grid, factors


def _times_for(case, steps):
    count = int(steps) + 1
    if count > case["frame_time"].size:
        raise ValueError("requested steps exceed the stored frames")
    return np.array(case["frame_time"][:count], dtype=float, copy=True)


def _headline_from_comparison(comparison, times, case, weights):
    change = comparison["occupation_full"] - comparison["occupation_full"][0]
    return {
        "complete": True,
        "times": [float(item) for item in times],
        "t_start": float(times[0]),
        "t_end": float(times[-1]),
        "steps": int(times.size - 1),
        "case": case["case"],
        "occupation_error": comparison["occupation_error"],
        "coherence_error": comparison["coherence_error"],
        "phase_error": comparison["phase_error"],
        "occupation_change_max_abs": float(np.max(np.abs(change))),
        "allocation": comparison["allocation"],
        "trapezoid_residual_max": comparison["trapezoid_residual_max"],
        "history_norm": comparison["history_norm"],
        "initial_exterior_norm": comparison["initial_exterior_norm"],
        "observer_preserved": True,
        "weights": [float(item) for item in weights],
        "backend": BACKEND,
    }


def _effect_record(panel, reference_error):
    records = {}
    threshold = float(reference_error["max_abs"])
    for name, item in panel.items():
        absolute = float(item["occupation_error_against_full"]["max_abs"])
        records[name] = {
            "occupation_error_against_full": item["occupation_error_against_full"],
            "history_norm": item["history_norm"],
            "initial_exterior_norm": item["initial_exterior_norm"],
            "exceeds_reduction_error": bool(absolute > max(threshold, 0.0)),
        }
    active = [name for name, item in records.items() if item["exceeds_reduction_error"]]
    return {"omissions": records, "active": active}


def _load_npz(path):
    if not Path(path).is_file():
        return {}
    with np.load(path, allow_pickle=False) as stored:
        return {key: np.array(stored[key]) for key in stored.files}


def _endpoint_pair(approximate, reference, change):
    absolute = float(np.max(np.abs(np.asarray(approximate) - np.asarray(reference))))
    change_scale = float(np.max(np.abs(change)))
    level = float(np.max(np.abs(reference)))
    return {
        "max_abs": absolute,
        "signal_max_abs": level,
        "relative_to_signal": None if level == 0.0 else absolute / level,
        "change_max_abs": change_scale,
        "relative_to_change": None if change_scale == 0.0 else absolute / change_scale,
    }


def append_full_stored_window(
    case="nf512_dt_0_0005",
    budget_s=600.0,
    resume=False,
    substeps=2,
    output_json=DEFAULT_JSON,
    output_npz=DEFAULT_NPZ,
):
    """Append the stored-frame window without replacing the pilot record."""
    output_json = Path(output_json)
    output_npz = Path(output_npz)
    if not output_json.is_file():
        raise FileNotFoundError("the pilot record must already exist")
    previous = json.loads(output_json.read_text())
    if previous.get("headline", {}).get("t_end") == 0.05 and previous.get("requested_domain") == "pilot":
        raise RuntimeError("pilot headline was already rewritten to T=0.05")
    case_data = load_episode_case(case)
    identity = request_identity(
        "stored-frames",
        case,
        substeps,
        source_id(case_data["phi0"], case_data["phi1"], case_data["weights"]),
    )
    if resume:
        disposition = resume_disposition(previous, identity)
        if disposition == "return":
            previous["resumed_complete"] = True
            return previous
        if disposition == "reject":
            return {
                "schema": SCHEMA,
                "status": "REJECTED_REQUEST_MISMATCH",
                "resumed_complete": False,
                "rejected_identity": identity,
                "preserved_record": True,
                "record_path": str(output_json),
            }
    if previous.get("phases", {}).get("full_window", {}).get("complete") and previous.get("requests", {}).get("stored-frames") == identity:
        previous["resumed_complete"] = True
        return previous
    preserved = {
        "cpu_seconds_this_process": previous.get("cpu_seconds_this_process"),
        "headline": json.loads(json.dumps(previous["headline"])),
        "probe_seconds": previous["phases"]["probe"]["seconds"],
        "comparison_t_end": previous["phases"]["comparison"]["t_end"],
        "cross_norm": previous["initial_source_blocks"]["cross_norm"],
        "omission_active": list(previous["phases"]["omissions"]["active"]),
    }
    previous_arrays = _load_npz(output_npz)
    record = previous
    record["requests"] = dict(record.get("requests") or {})
    record["requests"]["stored-frames"] = identity
    record.setdefault("requests", {})
    if "pilot" not in record["requests"]:
        record["requests"]["pilot"] = request_identity(
            record.get("requested_domain"),
            record.get("requested_case"),
            int(record.get("requested_substeps", 1)),
            identity["source"],
        )
    budget = MeasuredBudget(float(budget_s))

    def checkpoint():
        write_json(output_json, record)

    if not budget.allows(0.0):
        record["full_window_status"] = "BUDGET_EXCEEDED"
        checkpoint()
        return record
    grid, factors = build_case_grid(case_data)
    observer, columns, _gram = region_observer(case_data["phi0"], case_data["phi1"], region=0)
    if not np.array_equal(observer, columns[:, :2]):
        raise RuntimeError("observer is not the initial region-0 pair")
    weights = np.array(case_data["weights"], dtype=float, copy=True)
    times = np.array(case_data["frame_time"], dtype=float, copy=True)
    hamiltonian = hamiltonian_from_radial(factors, stored_radial(case_data))
    estimate = float(record.get("estimate_all_stored_frames_substeps_2_s") or 120.0)
    if int(substeps) != 2:
        estimate = quadratic_estimate(float(record["phases"]["probe"]["seconds"]), 2, times.size - 1, substeps=int(substeps))
    if not budget.allows(estimate):
        record["full_window_status"] = "BUDGET_EXCEEDED"
        record["full_window_stop_reason"] = "stored-frame estimate exceeded the additional budget"
        checkpoint()
        return record
    started = time.process_time()
    comparison = compare_retained(
        hamiltonian,
        observer,
        times,
        columns,
        weights,
        substeps=int(substeps),
    )
    streamed_seconds = float(time.process_time() - started)
    saved_initial = spinor_local_moments(observer, case_data["phi0"], case_data["phi1"], weights)
    saved_final = spinor_local_moments(
        observer,
        case_data["final_phi0"],
        case_data["final_phi1"],
        weights,
    )
    autonomous_change = np.asarray(saved_final["occupation"]) - np.asarray(saved_initial["occupation"])
    full_endpoint = comparison["occupation_full"][-1]
    retained_endpoint = comparison["occupation_retained"][-1]
    variants = interpolation_variant_report(times, case_data["frame_Q"])
    variants["saved_final_rate_versus_last_secant"] = saved_rate_versus_frame_secant(
        grid,
        times,
        case_data["frame_Q"],
        case_data["final_Q_dot"],
    )
    hermite_endpoint = None
    if budget.allows(max(5.0, 0.25 * streamed_seconds)):
        hermite_radial = lambda time: hermite_rows(times, case_data["frame_Q"], [float(time)])[0]
        hermite_hamiltonian = hamiltonian_from_radial(factors, hermite_radial)
        hermite_occupation, hermite_coherence = project_full_occupation(
            hermite_hamiltonian,
            observer,
            times,
            columns,
            weights,
            substeps=int(substeps),
        )
        hermite_endpoint = {
            "occupation_error_against_declared_full": _endpoint_pair(
                hermite_occupation[-1],
                full_endpoint,
                comparison["occupation_full"] - comparison["occupation_full"][0],
            ),
            "coherence_error_against_declared_full": signal_error(
                hermite_coherence,
                comparison["coherence_full"],
            ),
            "role": "indicator of a prescribed finite-difference Hermite Q, not the declared evolution",
        }
    omissions = None
    omission_arrays = {}
    if budget.allows(2.15 * streamed_seconds):
        panel = {}
        for name, flags in (
            ("memory", {"memory": False}),
            ("outside_drive", {"outside_drive": False}),
        ):
            if not budget.allows(streamed_seconds):
                break
            _result, occupation, coherence = evolve_streamed(
                hamiltonian,
                observer,
                times,
                columns,
                weights,
                propagator_substeps=int(substeps),
                **flags,
            )
            panel[name] = {
                "occupation_error_against_full": signal_error(occupation, comparison["occupation_full"]),
                "coherence_error_against_full": signal_error(coherence, comparison["coherence_full"]),
                "history_norm": float(np.linalg.norm(_result["history"])),
                "initial_exterior_norm": float(np.linalg.norm(_result["initial_exterior"])),
            }
            omission_arrays[name] = occupation
        if panel:
            omissions = _effect_record(panel, comparison["occupation_error"])
    cross = correlated_cross_separation()
    field_cross = frame_source_blocks(comparison["result"]["frame"], columns, weights)
    record["phases"]["full_window"] = {
        "complete": True,
        "case": case,
        "substeps": int(substeps),
        "seconds": streamed_seconds,
        "additional_cpu_seconds": budget.elapsed(),
        "t_start": float(times[0]),
        "t_end": float(times[-1]),
        "steps": int(times.size - 1),
        "declared_geometry": DECLARED_INTERPOLATION,
        "comparison": _headline_from_comparison(comparison, times, case_data, weights),
        "saved_autonomous": {
            "field_frames_stored": "initial and final spinors only",
            "initial_occupation": saved_initial["occupation"],
            "final_occupation": saved_final["occupation"],
            "initial_coherence": {
                "real": float(saved_initial["coherence"].real),
                "imag": float(saved_initial["coherence"].imag),
            },
            "final_coherence": {
                "real": float(saved_final["coherence"].real),
                "imag": float(saved_final["coherence"].imag),
            },
            "occupation_change": [float(item) for item in autonomous_change],
            "conditional_full_versus_saved": _endpoint_pair(
                full_endpoint,
                saved_final["occupation"],
                autonomous_change,
            ),
            "conditional_retained_versus_saved": _endpoint_pair(
                retained_endpoint,
                saved_final["occupation"],
                autonomous_change,
            ),
            "volterra_versus_conditional_full": comparison["occupation_error"],
            "volterra_coherence": comparison["coherence_error"],
            "volterra_phase": comparison["phase_error"],
        },
        "interpolation_variants": variants,
        "hermite_indicator": hermite_endpoint,
        "omissions": omissions,
        "field_cross_norm": field_cross["cross_norm"],
        "field_cross_active": bool(field_cross["cross_norm"] > 1e-8),
        "correlated_cross_control": cross,
        "allocation": comparison["allocation"],
        "pilot_window_unchanged": [preserved["headline"]["t_start"], preserved["headline"]["t_end"]],
    }
    record["full_window_cpu_seconds"] = budget.elapsed()
    record["full_window_status"] = "MEASURED_STORED_FRAMES"
    record["stress_claimed"] = False
    if (
        preserved["cpu_seconds_this_process"] != record.get("cpu_seconds_this_process")
        or preserved["headline"] != record["headline"]
        or preserved["probe_seconds"] != record["phases"]["probe"]["seconds"]
        or preserved["comparison_t_end"] != record["phases"]["comparison"]["t_end"]
        or preserved["cross_norm"] != record["initial_source_blocks"]["cross_norm"]
        or preserved["omission_active"] != list(record["phases"]["omissions"]["active"])
    ):
        raise RuntimeError("appending the full window changed the pilot record")
    arrays = dict(previous_arrays)
    arrays["full_times"] = times
    arrays["full_occupation_retained"] = comparison["occupation_retained"]
    arrays["full_occupation_full"] = comparison["occupation_full"]
    arrays["full_coherence_retained"] = comparison["coherence_retained"]
    arrays["full_coherence_full"] = comparison["coherence_full"]
    for name, occupation in omission_arrays.items():
        arrays["full_occupation_" + name + "_off"] = occupation
    for key, value in previous_arrays.items():
        if not np.array_equal(arrays[key], value):
            raise RuntimeError(f"pilot array {key} changed while appending the full window")
    write_npz(output_npz, arrays)
    record["payload_bytes"] = int(output_json.stat().st_size + output_npz.stat().st_size)
    checkpoint()
    record["payload_bytes"] = int(output_json.stat().st_size + output_npz.stat().st_size)
    record["payload_within_64MiB"] = bool(record["payload_bytes"] <= PAYLOAD_LIMIT_BYTES)
    checkpoint()
    if not record["payload_within_64MiB"]:
        raise RuntimeError(f"payload {record['payload_bytes']} exceeds 64MiB")
    if record["headline"]["t_end"] != preserved["comparison_t_end"]:
        raise RuntimeError("pilot end time changed")
    return record


def execute(
    domain="pilot",
    case="nf512_dt_0_0005",
    budget_s=60.0,
    resume=False,
    substeps=1,
    output_json=DEFAULT_JSON,
    output_npz=DEFAULT_NPZ,
):
    """Run the gated consumer. The default domain selects a prefix inside the budget."""
    if domain == "stored-frames":
        return append_full_stored_window(
            case=case,
            budget_s=budget_s,
            resume=resume,
            substeps=substeps,
            output_json=output_json,
            output_npz=output_npz,
        )
    output_json = Path(output_json)
    output_npz = Path(output_npz)
    budget = MeasuredBudget(float(budget_s))
    record = {
        "schema": SCHEMA,
        "status": "RUNNING",
        "renewal": False,
        "regeneration": False,
        "coupled_production_comparison": False,
        "stress_claimed": False,
        "stress_gap": STRESS_GAP,
        "interface_gaps": list(INTERFACE_GAPS),
        "backend": BACKEND,
        "phases": {},
        "budget_s": float(budget_s),
        "safety": BUDGET_SAFETY,
    }
    if resume and output_json.is_file():
        previous = json.loads(output_json.read_text())
        episode_case = load_episode_case(case)
        identity = request_identity(
            domain,
            case,
            substeps,
            source_id(episode_case["phi0"], episode_case["phi1"], episode_case["weights"]),
        )
        disposition = resume_disposition(previous, identity)
        if disposition == "return":
            previous["resumed_complete"] = True
            return previous
        if disposition == "reject":
            return {
                "schema": SCHEMA,
                "status": "REJECTED_REQUEST_MISMATCH",
                "resumed_complete": False,
                "rejected_identity": identity,
                "preserved_record": True,
                "record_path": str(output_json),
            }
        if previous.get("requested_domain") == domain and previous.get("requested_case") == case:
            record = previous
            record["status"] = "RUNNING"
            record["resumed"] = True
    record["requested_case"] = case
    record["requested_domain"] = domain
    record["requested_substeps"] = int(substeps)

    def checkpoint():
        record["cpu_seconds_this_process"] = budget.elapsed()
        write_json(output_json, record)

    checkpoint()
    if not record["phases"].get("verify_small", {}).get("complete"):
        record["phases"]["verify_small"] = verify_small_control()
        checkpoint()
    case_data = load_episode_case(case)
    grid, factors = build_case_grid(case_data)
    if not record["phases"].get("verify_stored", {}).get("complete"):
        record["phases"]["verify_stored"] = verify_stored_operator(grid, factors, case_data)
        checkpoint()
    observer, columns, gram = region_observer(case_data["phi0"], case_data["phi1"], region=0)
    weights = np.array(case_data["weights"], dtype=float, copy=True)
    if not np.array_equal(weights, case_data["weights"]):
        raise RuntimeError("source weights were copied incorrectly")
    gram_defect = float(np.max(np.abs(gram - np.eye(gram.shape[0]))))
    hamiltonian = hamiltonian_from_radial(factors, stored_radial(case_data))
    available = int(case_data["frame_time"].size - 1)
    probe_steps = 2
    probe = None
    if not record["phases"].get("probe", {}).get("complete"):
        if not budget.allows(0.0):
            record["status"] = "BUDGET_EXCEEDED"
            record["stop_reason"] = "no time remained before the cost probe"
            checkpoint()
            return record
        probe_times = _times_for(case_data, probe_steps)
        started = time.process_time()
        probe = compare_retained(hamiltonian, observer, probe_times, columns, weights, substeps=1)
        spent = float(time.process_time() - started)
        record["phases"]["probe"] = {
            "complete": True,
            "seconds": spent,
            "steps": probe_steps,
            "comparison": _headline_from_comparison(probe, probe_times, case_data, weights),
            "hamiltonian_builds": int(hamiltonian.stats["builds"]),
            "hamiltonian_hits": int(hamiltonian.stats["hits"]),
            "matrix_bytes": int(hamiltonian.stats["matrix_bytes"]),
        }
        checkpoint()
    probe_phase = record["phases"]["probe"]
    probe_seconds = float(probe_phase["seconds"])
    if domain == "stored-frames":
        requested = available
        estimate = quadratic_estimate(probe_seconds, probe_steps, requested, substeps=int(substeps))
        # The independent full evolution and three omissions are part of the claim.
        estimate += quadratic_estimate(probe_seconds, probe_steps, requested, substeps=1)
        estimate += 3.0 * quadratic_estimate(probe_seconds, probe_steps, min(requested, 4), substeps=1)
        plan = {
            "main_steps": requested if budget.allows(estimate) else probe_steps,
            "refine": False,
            "omissions": budget.allows(estimate),
            "independent_full": True,
            "estimate_s": estimate,
            "admitted": bool(budget.allows(estimate)),
            "forced_domain": "stored-frames",
        }
        if not plan["admitted"]:
            plan["main_steps"] = probe_steps
            plan["omissions"] = False
    else:
        plan = plan_remaining(
            budget.remaining(),
            probe_seconds,
            probe_steps,
            available,
        )
        plan["admitted"] = True
        plan["forced_domain"] = "pilot"
    record["phases"]["plan"] = plan
    record["estimate_all_stored_frames_s"] = quadratic_estimate(probe_seconds, probe_steps, available, substeps=1)
    record["estimate_all_stored_frames_substeps_2_s"] = quadratic_estimate(
        probe_seconds, probe_steps, available, substeps=2
    )
    record["estimate_hundred_output_nodes_s"] = quadratic_estimate(probe_seconds, probe_steps, 100, substeps=1)
    record["hundred_node_geometry"] = "not stored; extrapolation assumes a Q sample at every output node"
    checkpoint()

    probe_times = _times_for(case_data, probe_steps)
    fresh_probe = probe
    if fresh_probe is None and budget.allows(probe_seconds):
        fresh_probe = compare_retained(hamiltonian, observer, probe_times, columns, weights, substeps=1)
        record["phases"]["probe"]["comparison"] = _headline_from_comparison(
            fresh_probe, probe_times, case_data, weights
        )
    headline = record["phases"]["probe"]["comparison"]
    headline_comparison = fresh_probe
    if plan["main_steps"] != probe_steps and not record["phases"].get("comparison", {}).get("complete"):
        main_times = _times_for(case_data, plan["main_steps"])
        estimate = quadratic_estimate(probe_seconds, probe_steps, plan["main_steps"], substeps=int(substeps))
        if budget.allows(estimate):
            headline_comparison = compare_retained(
                hamiltonian,
                observer,
                main_times,
                columns,
                weights,
                substeps=int(substeps),
            )
            headline = _headline_from_comparison(headline_comparison, main_times, case_data, weights)
            record["phases"]["comparison"] = headline
            checkpoint()
        else:
            plan["main_steps"] = probe_steps
            record["phases"]["plan"] = plan
    elif record["phases"].get("comparison", {}).get("complete"):
        headline = record["phases"]["comparison"]

    arrays = {}
    if headline_comparison is not None and int(headline["steps"]) != probe_steps:
        arrays["occupation_retained"] = headline_comparison["occupation_retained"]
        arrays["occupation_full"] = headline_comparison["occupation_full"]
        arrays["coherence_retained"] = headline_comparison["coherence_retained"]
        arrays["coherence_full"] = headline_comparison["coherence_full"]
    elif fresh_probe is not None:
        arrays["occupation_retained"] = fresh_probe["occupation_retained"]
        arrays["occupation_full"] = fresh_probe["occupation_full"]
        arrays["coherence_retained"] = fresh_probe["coherence_retained"]
        arrays["coherence_full"] = fresh_probe["coherence_full"]
        headline = record["phases"]["probe"]["comparison"]

    if fresh_probe is not None:
        blocks = frame_source_blocks(fresh_probe["result"]["frame"], columns, weights)
        record["initial_source_blocks"] = {
            "parent_norm": blocks["parent_norm"],
            "cross_norm": blocks["cross_norm"],
            "parent_diagonal": [float(item) for item in np.real(np.diag(blocks["parent"]))],
        }
    else:
        record["initial_source_blocks"] = None

    if (
        plan.get("refine")
        and fresh_probe is not None
        and budget.allows(quadratic_estimate(probe_seconds, probe_steps, probe_steps * 2))
    ):
        fine_steps = probe_steps * 2
        fine_times = np.linspace(float(probe_times[0]), float(probe_times[-1]), fine_steps + 1)
        refined = compare_retained(hamiltonian, observer, fine_times, columns, weights, substeps=1)
        coarse_error = signal_error(fresh_probe["occupation_retained"], refined["occupation_full"][::2])
        fine_error = refined["occupation_error"]
        record["phases"]["refinement"] = {
            "complete": True,
            "fine_steps": fine_steps,
            "window": [float(probe_times[0]), float(probe_times[-1])],
            "coarse_against_fine_full": coarse_error,
            "fine_against_fine_full": fine_error,
            "error_decreased": bool(fine_error["max_abs"] < coarse_error["max_abs"]),
        }
        arrays["occupation_refined_retained"] = refined["occupation_retained"]
        arrays["occupation_refined_full"] = refined["occupation_full"]
        checkpoint()

    if plan.get("omissions") and fresh_probe is not None and budget.allows(3.0 * probe_seconds):
        panel = {}
        for name, flags in (
            ("memory", {"memory": False}),
            ("outside_drive", {"outside_drive": False}),
            ("coupling", {"coupling": False}),
        ):
            result, occupation, _coherence = evolve_streamed(
                hamiltonian,
                observer,
                probe_times,
                columns,
                weights,
                **flags,
            )
            panel[name] = {
                "occupation_error_against_full": signal_error(occupation, fresh_probe["occupation_full"]),
                "history_norm": float(np.linalg.norm(result["history"])),
                "initial_exterior_norm": float(np.linalg.norm(result["initial_exterior"])),
            }
            arrays["occupation_" + name + "_off"] = occupation
        record["phases"]["omissions"] = _effect_record(
            panel,
            record["phases"]["probe"]["comparison"]["occupation_error"],
        )
        checkpoint()

    full_window = int(headline["steps"]) == available and domain == "stored-frames" and plan.get("admitted") is True
    record["finite_approximation_domain"] = finite_domain(
        case_data,
        headline,
        factors,
        gram_defect,
        full_window,
    )
    record["headline"] = headline
    record["bindings"] = {
        "names": case_data["binding_names"],
        "values": [float(item) for item in case_data["binding_values"]],
        "time_resolution_geometry_gap": case_data["time_resolution_geometry_gap"],
        "interpolation": case_data["interpolation"],
        "geometry_rate_series_stored": False,
        "final_coarse_Q_dot_norm": case_data["final_coarse_Q_dot_norm"],
        "dirac_geometry_argument": case_data["dirac_geometry_argument"],
        "stored_but_not_dirac_matrix_entries": case_data["stored_but_not_dirac_matrix_entries"],
        "episode_verdict_label": case_data["episode_verdict_label"],
        "episode_verdict_used_as_acceptance": False,
        "episode_npz_sha256": sha256_file(EPISODE_NPZ),
    }
    record["gram_defect"] = gram_defect
    record["next_commands"] = next_commands(case, full_window, int(substeps))
    record["next_command_claim"] = NEXT_CLAIM
    record["hamiltonian_builds"] = int(hamiltonian.stats["builds"])
    record["hamiltonian_hits"] = int(hamiltonian.stats["hits"])
    if domain == "stored-frames" and not plan.get("admitted"):
        record["status"] = "BUDGET_EXCEEDED"
        record["stop_reason"] = "full stored-frame estimate exceeded the measured budget; the probe prefix was retained"
    else:
        record["status"] = "MEASURED_CONDITIONAL_REDUCTION"
        record["stop_reason"] = None
    record["cpu_seconds_this_process"] = budget.elapsed()
    if "occupation_retained" not in arrays:
        record["status"] = "BUDGET_EXCEEDED"
        record["stop_reason"] = "no comparison series was retained inside the budget"
        checkpoint()
        return record
    arrays["times"] = np.asarray(headline["times"], dtype=float)
    if record["phases"].get("refinement", {}).get("complete"):
        arrays["refinement_error_decreased"] = np.bool_(record["phases"]["refinement"]["error_decreased"])
    write_npz(output_npz, arrays)
    record["payload_bytes"] = int(output_json.stat().st_size + output_npz.stat().st_size)
    # json size is not final until this field is written; correct once.
    write_json(output_json, record)
    record["payload_bytes"] = int(output_json.stat().st_size + output_npz.stat().st_size)
    record["payload_within_64MiB"] = bool(record["payload_bytes"] <= PAYLOAD_LIMIT_BYTES)
    write_json(output_json, record)
    if not record["payload_within_64MiB"]:
        raise RuntimeError(f"payload {record['payload_bytes']} exceeds 64MiB")
    return record


def finite_domain(case, headline, factors, gram_defect, full_window):
    return {
        "case": case["case"],
        "fermion_band": int(factors["nf"]),
        "quadrature": int(factors["nq"]),
        "length": float(factors["length"]),
        "spinor_dimension": int(2 * factors["nf"]),
        "observer": "initial mode columns 0 and 1, weights 0.75 and 0.75, no phase QR",
        "source": "six initial episode columns and stored occupation weights",
        "gram_defect": float(gram_defect),
        "time_window": [float(headline["t_start"]), float(headline["t_end"])],
        "output_steps": int(headline["steps"]),
        "stored_frame_count": int(case["frame_time"].size),
        "geometry": "piecewise linear stored fine Q; static calibration lapse and shift",
        "interpolation_scale": case["interpolation"]["midpoint_linear_scale"],
        "time_resolution_geometry_gap": case["time_resolution_geometry_gap"],
        "reference": "independent full-band midpoint column evolution on the same H(Q(t))",
        "reduction": "nsc_evolving_reduction.evolve_retained_region backend='streamed'",
        "full_stored_window": bool(full_window),
        "regeneration": False,
        "production_rk4_nodes": False,
        "stress": False,
    }


def next_commands(case, full_window, substeps):
    commands = []
    if not full_window or int(substeps) < 2:
        commands.append(
            "python scripts/lab.py scripts/derive_nsc_coupled_local_response.py "
            f"--domain stored-frames --case {case} --substeps 2 --budget-s 600"
        )
    partner = PARTNER_CASE[case]
    commands.append(
        "python scripts/lab.py scripts/derive_nsc_coupled_local_response.py "
        f"--domain stored-frames --case {partner} --substeps 2 --budget-s 600"
    )
    return commands


SCHEMA_V2 = "NSC-COUPLED-LOCAL-RESPONSE-v2"
REGENERATION_NPZ = _LAB / "results" / "development" / "nsc-regeneration-episode-v1.npz"
REGENERATION_JSON = _LAB / "results" / "development" / "nsc-regeneration-episode-v1.json"
DEFAULT_JSON_V2 = _LAB / "results" / "development" / "nsc-coupled-local-response-v2.json"
DEFAULT_NPZ_V2 = _LAB / "results" / "development" / "nsc-coupled-local-response-v2.npz"
SUCCESSOR_CASE = "nf512_dtmax_0_00025"
SUCCESSOR_PARTNER = "nf512_dtmax_0_0005"
SUCCESSOR_HANDOFF = "nf512_dt_0_0005"
SUCCESSOR_ORIGIN = "nf512"
SUCCESSOR_T0 = 0.05
SUCCESSOR_T1 = 0.085
SUCCESSOR_PROBE_STEPS = 2
SUCCESSOR_DECLARED_SUBSTEPS = 2
SUCCESSOR_PILOT_BUDGET_S = 600.0
SUCCESSOR_CONFIRM_BUDGET_S = 1800.0
SUCCESSOR_ONE_PERCENT = 0.01
SUCCESSOR_GRAM_GATE = 1e-8
SUCCESSOR_JOB = "634b7687-d576-4638-9dcd-101c8fe7d646"
SUCCESSOR_INTERPOLATION = "piecewise-linear-stored-quad-Q"
HUNDRED_NODE_OUTPUTS = 100


def array_id(values):
    array = np.ascontiguousarray(values)
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode("ascii"))
    digest.update(str(array.shape).encode("ascii"))
    digest.update(array.view(np.uint8))
    return digest.hexdigest()


def fraction_of_effect(movement, effect):
    movement = float(movement)
    effect = float(effect)
    if not np.isfinite(movement) or not np.isfinite(effect):
        raise ValueError("movement and effect must be finite")
    if effect <= 0.0:
        return {
            "movement": movement,
            "effect": effect,
            "fraction": None,
            "within_one_percent": None,
            "reason": "effect is zero",
        }
    fraction = movement / effect
    return {
        "movement": movement,
        "effect": effect,
        "fraction": float(fraction),
        "within_one_percent": bool(fraction <= SUCCESSOR_ONE_PERCENT),
        "reason": None,
    }


def hermite_with_slopes(times, values, slopes, query):
    """Cubic Hermite samples that use supplied nodal slopes.

    Stored nodes are returned exactly. ``slopes`` are time derivatives in
    the same layout as ``values``, not finite differences of those values.
    """
    times = np.asarray(times, dtype=float)
    values = np.asarray(values, dtype=float)
    slopes = np.asarray(slopes, dtype=float)
    query = np.atleast_1d(np.asarray(query, dtype=float))
    if values.shape != slopes.shape or values.shape[0] != times.shape[0]:
        raise ValueError("slopes must share the stored value layout")
    if times.size < 2 or np.any(np.diff(times) <= 0):
        raise ValueError("Hermite samples need increasing times")
    if np.any(query < times[0] - 1e-12) or np.any(query > times[-1] + 1e-12):
        raise ValueError("query leaves the stored geometry window")
    index = np.searchsorted(times, query, side="right") - 1
    index = np.clip(index, 0, times.size - 2)
    width = times[index + 1] - times[index]
    fraction = np.clip((query - times[index]) / width, 0.0, 1.0)
    cube = fraction ** 3
    square = fraction ** 2
    h00 = 2.0 * cube - 3.0 * square + 1.0
    h10 = cube - 2.0 * square + fraction
    h01 = -2.0 * cube + 3.0 * square
    h11 = cube - square
    return (
        h00[:, None] * values[index]
        + (h10 * width)[:, None] * slopes[index]
        + h01[:, None] * values[index + 1]
        + (h11 * width)[:, None] * slopes[index + 1]
    )


def prolonged_nodal_rates(grid, coarse_rates):
    """Prolong coarse Galerkin rates and record the pullback defect.

    ``pull(prolong(rate))`` must reproduce the coarse rate. A nonzero defect
    means the saved rate is not a nodal field in the prolongation range.
    """
    coarse = np.asarray(coarse_rates, dtype=float)
    fine = np.stack([prolong_geometry(grid, row) for row in coarse])
    pulled = np.stack([pull_geometry(grid, row) for row in fine])
    if pulled.shape != coarse.shape:
        raise ValueError("pulled nodal rate does not match the coarse rate")
    gap = float(np.max(np.abs(pulled - coarse)))
    return fine, {
        "projection": "prolong_geometry of the coarse Galerkin Q_dot; pull_geometry checks the range",
        "pullback_max_abs": gap,
        "rate_max_abs": float(np.max(np.abs(coarse))),
        "usable": bool(gap <= 1e-8 * max(1.0, float(np.max(np.abs(coarse))))),
    }


def one_body_blocks(observer, columns, weights):
    """Initial ``C_AA``, ``C_EE``, and ``C_AE`` of one weighted column state.

    The exterior block is not formed as a dense matrix. Its nonzero
    eigenvalues are the spectrum of the weighted exterior column Gram.
    """
    weights = np.asarray(weights, dtype=float)
    columns = np.asarray(columns, dtype=np.complex128)
    observer = np.asarray(observer, dtype=np.complex128)
    if columns.ndim != 2 or observer.ndim != 2 or observer.shape[0] != columns.shape[0]:
        raise ValueError("observer and columns must share the spinor dimension")
    if weights.shape != (columns.shape[1],):
        raise ValueError("weights must match the column count")
    covariance = source_covariance(columns, weights)
    coefficients = observer.conj().T @ columns
    parent = (coefficients * weights) @ coefficients.conj().T
    parent = 0.5 * (parent + parent.conj().T)
    product = observer.conj().T @ covariance
    parent_norm = float(np.linalg.norm(parent, ord="fro"))
    product_norm = float(np.linalg.norm(product, ord="fro"))
    cross = math.sqrt(max(0.0, product_norm * product_norm - parent_norm * parent_norm))
    covariance_norm = float(np.linalg.norm(covariance, ord="fro"))
    exterior_norm = math.sqrt(max(0.0, covariance_norm * covariance_norm - parent_norm * parent_norm - 2.0 * cross * cross))
    scale = np.sqrt(np.clip(weights, 0.0, None))
    weighted = columns * scale
    gram = weighted.conj().T @ weighted
    gram = 0.5 * (gram + gram.conj().T)
    eigenvalues = np.linalg.eigvalsh(gram).real
    parent_eigenvalues = np.linalg.eigvalsh(parent).real
    exterior_columns = columns - observer @ coefficients
    exterior_weighted = exterior_columns * scale
    exterior_gram = exterior_weighted.conj().T @ exterior_weighted
    exterior_gram = 0.5 * (exterior_gram + exterior_gram.conj().T)
    exterior_eigenvalues = np.linalg.eigvalsh(exterior_gram).real
    hermitian = float(np.linalg.norm(covariance - covariance.conj().T, ord="fro"))
    admissible = bool(
        hermitian <= 1e-8
        and float(np.min(eigenvalues)) >= -1e-8
        and float(np.max(eigenvalues)) <= 1.0 + 1e-8
        and float(np.min(parent_eigenvalues)) >= -1e-8
        and float(np.max(parent_eigenvalues)) <= 1.0 + 1e-8
        and float(np.min(exterior_eigenvalues)) >= -1e-8
        and float(np.max(exterior_eigenvalues)) <= 1.0 + 1e-8
    )
    return {
        "parent": parent,
        "parent_eigenvalues": np.array(parent_eigenvalues, dtype=float, copy=True),
        "cross_frobenius": float(cross),
        "exterior_frobenius": float(exterior_norm),
        "covariance_frobenius": covariance_norm,
        "covariance_eigenvalues": np.array(eigenvalues, dtype=float, copy=True),
        "exterior_eigenvalues": np.array(exterior_eigenvalues, dtype=float, copy=True),
        "hermitian_defect": hermitian,
        "exterior_column_norm": float(np.linalg.norm(exterior_columns)),
        "admissible": admissible,
        "dense_exterior_block_stored": False,
    }


def _complex_pair(value):
    number = complex(value)
    return {"real": float(number.real), "imag": float(number.imag)}


def _block_record(blocks):
    parent = blocks["parent"]
    return {
        "C_AA": {
            "real": np.real(parent).tolist(),
            "imag": np.imag(parent).tolist(),
        },
        "C_AA_eigenvalues": [float(item) for item in blocks["parent_eigenvalues"]],
        "C_AA_frobenius": float(np.linalg.norm(parent, ord="fro")),
        "C_AE_frobenius": float(blocks["cross_frobenius"]),
        "C_EE_frobenius": float(blocks["exterior_frobenius"]),
        "covariance_eigenvalues": [float(item) for item in blocks["covariance_eigenvalues"]],
        "exterior_eigenvalues": [float(item) for item in blocks["exterior_eigenvalues"]],
        "hermitian_defect": float(blocks["hermitian_defect"]),
        "exterior_column_norm": float(blocks["exterior_column_norm"]),
        "admissible": bool(blocks["admissible"]),
        "dense_C_EE_stored": False,
        "cross_active_as_a_block": bool(blocks["cross_frobenius"] > 1e-8),
    }


def _sample_series(times, values, query):
    times = np.asarray(times, dtype=float)
    values = np.asarray(values)
    query = np.asarray(query, dtype=float)
    selected = []
    for mark in query:
        index = int(np.argmin(np.abs(times - float(mark))))
        if abs(float(times[index]) - float(mark)) > 1e-12:
            raise ValueError(f"series has no sample at {mark}")
        selected.append(values[index])
    return np.stack(selected)


def load_successor_inputs(case=SUCCESSOR_CASE, regeneration_npz=REGENERATION_NPZ, regeneration_json=REGENERATION_JSON, origin_npz=EPISODE_NPZ):
    """Load the fixed T=0 observer and the transported T=0.05 state.

    No geometry step is taken. The observer is not rebuilt from ``Φ(0.05)``.
    """
    case = str(case)
    if case == SUCCESSOR_CASE:
        partner = SUCCESSOR_PARTNER
    elif case == SUCCESSOR_PARTNER:
        partner = SUCCESSOR_CASE
    else:
        raise ValueError(f"successor case must be an nf512 regeneration run, not {case}")
    resolution = SUCCESSOR_ORIGIN
    with np.load(origin_npz, allow_pickle=False) as origin:
        origin_phi0 = np.array(origin[resolution + "_initial_phi0"], copy=True)
        origin_phi1 = np.array(origin[resolution + "_initial_phi1"], copy=True)
        weights = np.array(origin[resolution + "_occupations"], dtype=float, copy=True)
        names = np.array(origin[resolution + "_binding_names"])
        values = np.array(origin[resolution + "_binding_values"], dtype=float, copy=True)
        handoff_phi0 = np.array(origin[SUCCESSOR_HANDOFF + "_final_phi0"], copy=True)
        handoff_phi1 = np.array(origin[SUCCESSOR_HANDOFF + "_final_phi1"], copy=True)
    with np.load(regeneration_npz, allow_pickle=False) as payload:
        prefix = case + "_"
        other = partner + "_"
        frame_time = np.array(payload[prefix + "frame_time"], dtype=float, copy=True)
        frame_q = np.array(payload[prefix + "frame_quad_Q"], dtype=float, copy=True)
        coarse_q = np.array(payload[prefix + "frame_coarse_Q"], dtype=float, copy=True)
        coarse_rate = np.array(payload[prefix + "frame_Q_dot"], dtype=float, copy=True)
        p_chi = np.array(payload[prefix + "frame_coarse_p_chi"], dtype=float, copy=True)
        p_chi_dot = np.array(payload[prefix + "frame_p_chi_dot"], dtype=float, copy=True)
        increment_q_dot = np.array(payload[prefix + "frame_increment_Q_dot"], dtype=float, copy=True)
        increment_dt = np.array(payload[prefix + "frame_increment_dt"], dtype=float, copy=True)
        indicator_q = np.array(payload[prefix + "frame_indicator_Q_dot"], dtype=float, copy=True)
        phi0 = np.array(payload[prefix + "frame_phi0"], copy=True)
        phi1 = np.array(payload[prefix + "frame_phi1"], copy=True)
        partner_phi0 = np.array(payload[other + "frame_phi0"], copy=True)
        partner_phi1 = np.array(payload[other + "frame_phi1"], copy=True)
        partner_q = np.array(payload[other + "frame_quad_Q"], dtype=float, copy=True)
        partner_time = np.array(payload[other + "frame_time"], dtype=float, copy=True)
        series_time = np.array(payload[prefix + "time"], dtype=float, copy=True)
        proper_clock = np.array(payload[prefix + "proper_clock"], dtype=float, copy=True)
        leader_clock = np.array(payload[prefix + "leader_clock"], dtype=float, copy=True)
        field_energy = np.array(payload[prefix + "field_energy"], dtype=float, copy=True)
        window_normal = np.array(payload[prefix + "window_normal"], dtype=float, copy=True)
        leader = np.array(payload[prefix + "leader"], dtype=float, copy=True)
        leader_share = np.array(payload[prefix + "leader_share"], dtype=float, copy=True)
        gram_gap = np.array(payload[prefix + "gram_gap"], dtype=float, copy=True)
    if not np.array_equal(weights, np.array([0.75, 0.75, 0.5, 0.5, 0.25, 0.25])):
        raise RuntimeError("origin occupations are not the original six Gaussian weights")
    if frame_time.size < 3 or abs(float(frame_time[0]) - SUCCESSOR_T0) > 1e-12:
        raise RuntimeError("successor frames do not start at T=0.05")
    if abs(float(frame_time[-1]) - SUCCESSOR_T1) > 1e-12:
        raise RuntimeError("successor frames do not end at the declared comparison time")
    if phi0.shape[0] != frame_time.size or phi1.shape != phi0.shape:
        raise RuntimeError("stored spinor frames do not match the geometry frames")
    if not np.array_equal(phi0[0], handoff_phi0) or not np.array_equal(phi1[0], handoff_phi1):
        raise RuntimeError("T=0.05 spinor is not the bitwise spherical-episode handoff")
    if not np.array_equal(partner_time, frame_time):
        raise RuntimeError("partner frames do not share the primary times")
    if float(np.min(frame_q)) <= 0.0 or float(np.max(gram_gap)) > SUCCESSOR_GRAM_GATE:
        raise RuntimeError("stored chart or Gram leaves the admitted window")
    observer, origin_columns, origin_gram = region_observer(origin_phi0, origin_phi1, region=0)
    if not np.array_equal(observer, origin_columns[:, :2]):
        raise RuntimeError("observer is not the T=0 region-0 pair")
    transported = np.concatenate((phi0[0], phi1[0]), axis=0)
    if np.array_equal(observer, transported[:, :2]):
        raise RuntimeError("observer was replaced by the transported T=0.05 basis")
    episode = json.loads(Path(regeneration_json).read_text())
    run = episode["runs"][case]
    expansion = run["expansion_frames"][-1]
    anti = expansion["anti_trapped_intervals"][0]
    structure = episode["structures"][case]["slice_ledgers"]
    clock = episode["clocks"][case]
    blocks = one_body_blocks(observer, transported, weights)
    if not blocks["admissible"]:
        raise RuntimeError("initial one-body covariance is not admissible")
    if blocks["cross_frobenius"] <= 1e-8:
        raise RuntimeError("actual C_AE at T=0.05 is inactive; a synthetic cross would be a substitution")
    return {
        "case": case,
        "partner_case": partner,
        "handoff": SUCCESSOR_HANDOFF,
        "requested_job_id": SUCCESSOR_JOB,
        "frame_time": frame_time,
        "frame_Q": frame_q,
        "coarse_Q": coarse_q,
        "coarse_Q_dot": coarse_rate,
        "coarse_p_chi": p_chi,
        "coarse_p_chi_dot": p_chi_dot,
        "increment_Q_dot": increment_q_dot,
        "increment_dt": increment_dt,
        "indicator_Q_dot": indicator_q,
        "phi0_frames": phi0,
        "phi1_frames": phi1,
        "partner_phi0": partner_phi0,
        "partner_phi1": partner_phi1,
        "partner_Q": partner_q,
        "phi0": origin_phi0,
        "phi1": origin_phi1,
        "weights": weights,
        "observer": observer,
        "origin_columns": origin_columns,
        "origin_gram_defect": float(np.max(np.abs(origin_gram - np.eye(origin_gram.shape[0])))),
        "transported_columns": transported,
        "binding_names": [str(item) for item in names.tolist()],
        "binding_values": values,
        "initial_blocks": blocks,
        "proper_clock": _sample_series(series_time, proper_clock, frame_time),
        "leader_clock": _sample_series(series_time, leader_clock, frame_time),
        "field_energy": _sample_series(series_time, field_energy, frame_time),
        "window_normal": _sample_series(series_time, window_normal, frame_time),
        "leader": _sample_series(series_time, leader, frame_time),
        "leader_share": _sample_series(series_time, leader_share, frame_time),
        "series_count": int(series_time.size),
        "clock_protocol": str(clock["protocol"]),
        "episode_stop_reason": str(run["stop_reason"]),
        "episode_attained_T": float(run["attained_T"]),
        "arc_x_start": float(anti["x_start"]),
        "arc_x_end": float(anti["x_end"]),
        "structure_maintained": bool(structure["content"]["maintained"]),
        "structure_renewed": bool(structure["reversal"]["renewed"]),
        "balanced_flow": bool(structure["flux"]["balanced"]),
        "joined_renewal_proxy": bool(episode["goal"]["joined_renewal_proxy"]),
        "old_toy_B_not_used": bool(episode["link_reports"]["512"]["old_B_not_used"]),
        "state_reset": bool(episode["state_reset"]),
        "new_force_added": bool(episode["new_force_added"]),
        "indicator_versus_bound": str(episode["indicator_versus_bound"]),
        "partner_phi_gap": float(np.max(np.abs(phi0 - partner_phi0) + np.abs(phi1 - partner_phi1))),
        "partner_Q_gap": float(np.max(np.abs(frame_q - partner_q))),
        "interpolation": interpolation_scale(frame_time, frame_q),
        "observer_id": array_id(observer),
        "weights_id": array_id(weights),
        "initial_state_id": array_id(transported),
        "geometry_id": array_id(frame_q),
        "episode_sha256": sha256_file(regeneration_npz),
        "episode_json_sha256": sha256_file(regeneration_json),
        "origin_sha256": sha256_file(origin_npz),
        "origin_json_sha256": sha256_file(_LAB / "results" / "development" / "nsc-spherical-feedback-episode-v1.json"),
    }


def successor_identity(inputs, substeps=SUCCESSOR_DECLARED_SUBSTEPS):
    return {
        "domain": "successor",
        "case": str(inputs["case"]),
        "partner_case": str(inputs["partner_case"]),
        "substeps": int(substeps),
        "observer": str(inputs["observer_id"]),
        "weights": str(inputs["weights_id"]),
        "initial_state": str(inputs["initial_state_id"]),
        "geometry": str(inputs["geometry_id"]),
        "episode_npz": str(inputs["episode_sha256"]),
        "origin_npz": str(inputs["origin_sha256"]),
        "backend": BACKEND,
        "interpolation": SUCCESSOR_INTERPOLATION,
        "memory": True,
        "outside_drive": True,
        "coupling": True,
        "drop_cross_covariance": False,
        "observer_protocol": "t0-region-0-columns-0-1-no-rephase",
        "initial_covariance": "actual-transported-phi-0.05",
    }


def successor_resume_disposition(previous, identity):
    """Return a finished successor only when every typed binding matches."""
    if not isinstance(previous, dict):
        return "continue"
    if previous.get("status") != "MEASURED_SUCCESSOR_RESPONSE":
        return "continue"
    if previous.get("schema") != SCHEMA_V2:
        return "reject"
    stored = previous.get("request")
    if stored != identity:
        return "reject"
    full = previous.get("phases", {}).get("full_window", {})
    if full.get("complete") is True and int(full.get("substeps", -1)) == int(identity["substeps"]):
        return "return"
    return "reject"


def autonomous_series(observer, phi0, phi1, weights):
    stacked = np.concatenate((np.asarray(phi0), np.asarray(phi1)), axis=1)
    projected = np.einsum("ij,tjk->tik", np.asarray(observer).conj().T, stacked)
    return occupations_and_coherence(projected, weights)


def _covariance_from_amplitudes(amplitudes, weights):
    scale = np.asarray(weights, dtype=float)
    covariances = np.empty((amplitudes.shape[0], amplitudes.shape[1], amplitudes.shape[1]), dtype=np.complex128)
    for index, sample in enumerate(amplitudes):
        matrix = (sample * scale) @ sample.conj().T
        covariances[index] = 0.5 * (matrix + matrix.conj().T)
    return covariances


def _matrix_series_gap(left, right):
    return float(np.max(np.linalg.norm(np.asarray(left) - np.asarray(right), axis=(1, 2), ord="fro")))


def column_cross_control(hamiltonian, observer, times, columns, weights, substeps=1, baseline=None):
    """Drop the actual initial cross by splitting retained and exterior columns.

    The split is the linear decomposition of the same weighted state. It is
    not a second source. Superposition against the unsplit evolution is the
    numerical error of this control.
    """
    coefficients = np.asarray(observer).conj().T @ np.asarray(columns)
    retained_columns = np.asarray(observer) @ coefficients
    exterior_columns = np.asarray(columns) - retained_columns
    if baseline is None:
        baseline, _occupation, _coherence = evolve_streamed(
            hamiltonian,
            observer,
            times,
            columns,
            weights,
            propagator_substeps=int(substeps),
        )
    retained, _occupation_a, _coherence_a = evolve_streamed(
        hamiltonian,
        observer,
        times,
        retained_columns,
        weights,
        propagator_substeps=int(substeps),
    )
    exterior, _occupation_e, _coherence_e = evolve_streamed(
        hamiltonian,
        observer,
        times,
        exterior_columns,
        weights,
        propagator_substeps=int(substeps),
    )
    combined = retained["amplitudes"] + exterior["amplitudes"]
    superposition = float(np.max(np.abs(combined - baseline["amplitudes"])))
    parent = _covariance_from_amplitudes(retained["amplitudes"], weights)
    child = _covariance_from_amplitudes(exterior["amplitudes"], weights)
    mixed = np.empty_like(parent)
    scale = np.asarray(weights, dtype=float)
    for index in range(parent.shape[0]):
        left = retained["amplitudes"][index]
        right = exterior["amplitudes"][index]
        term = (left * scale) @ right.conj().T + (right * scale) @ left.conj().T
        mixed[index] = 0.5 * (term + term.conj().T)
    total = parent + child + mixed
    dropped = parent + child
    direct = _covariance_from_amplitudes(baseline["amplitudes"], weights)
    assembly = _matrix_series_gap(direct, total)
    separation = _matrix_series_gap(total, dropped)
    diagonal = np.real(np.diagonal(mixed, axis1=1, axis2=2))
    occupation_movement = float(np.max(np.abs(diagonal)))
    numerical = max(superposition, assembly)
    return {
        "covariance_total": total,
        "covariance_without_cross": dropped,
        "cross_covariance": mixed,
        "superposition_max_abs": superposition,
        "assembly_max_frobenius": assembly,
        "separation_max_frobenius": separation,
        "occupation_diagonal_movement": occupation_movement,
        "numerical_error": float(numerical),
        "active": bool(separation > 10.0 * max(numerical, 0.0) and occupation_movement > 10.0 * max(numerical, 0.0)),
        "synthetic_six_mode_substituted": False,
        "stores_dense_exterior_propagator": False,
        "history_width": int(baseline["history"].shape[-1]),
    }


def _effect_scale(reference):
    change = np.asarray(reference) - np.asarray(reference)[0]
    return float(np.max(np.abs(change)))


def _timed_comparison(hamiltonian, observer, times, columns, weights, substeps):
    started = time.process_time()
    result, occupation, coherence = evolve_streamed(
        hamiltonian,
        observer,
        times,
        columns,
        weights,
        propagator_substeps=int(substeps),
    )
    streamed_seconds = float(time.process_time() - started)
    started = time.process_time()
    full = evolve_full_columns(hamiltonian, times, columns, substeps=int(substeps))
    full_seconds = float(time.process_time() - started)
    projected = np.einsum("ij,tjk->tik", observer.conj().T, full)
    full_occupation, full_coherence = occupations_and_coherence(projected, weights)
    comparison = {
        "result": result,
        "occupation_retained": occupation,
        "occupation_full": full_occupation,
        "coherence_retained": coherence,
        "coherence_full": full_coherence,
        "occupation_error": signal_error(occupation, full_occupation),
        "coherence_error": signal_error(coherence, full_coherence),
        "phase_error": phase_error(coherence, full_coherence),
        "allocation": allocation_report(result),
        "trapezoid_residual_max": float(result["trapezoid_residual_max"]),
        "history_norm": float(np.linalg.norm(result["history"])),
        "initial_exterior_norm": float(np.linalg.norm(result["initial_exterior"])),
        "observer_preserved": bool(np.array_equal(result["local_basis"], observer)),
        "streamed_seconds": streamed_seconds,
        "full_seconds": full_seconds,
    }
    return comparison


def _forecast_streamed(probe_seconds, probe_steps, steps, substeps, probe_substeps=1):
    return quadratic_estimate(probe_seconds, probe_steps, steps, substeps=substeps, probe_substeps=probe_substeps)


def _forecast_full(probe_seconds, probe_steps, steps, substeps, probe_substeps=1):
    if probe_steps < 1 or steps < 1 or substeps < 1 or probe_substeps < 1:
        raise ValueError("step counts must be positive")
    return float(probe_seconds) * (float(steps) / float(probe_steps)) * (float(substeps) / float(probe_substeps))


def _indicator_status(movement, effect, numerical):
    report = fraction_of_effect(movement, effect)
    report["numerical_error"] = float(numerical)
    if effect <= max(float(numerical), 0.0):
        report["status"] = "inactive_effect_not_above_numerical_error"
        report["within_one_percent"] = None
        return report
    report["status"] = "resolved" if report["within_one_percent"] else "above_one_percent"
    return report


def _omission_against(occupation, reference, numerical):
    error = signal_error(occupation, reference)
    effect = _effect_scale(reference)
    return {
        "occupation_error_against_reference": error,
        "effect": _indicator_status(error["max_abs"], effect, numerical),
        "exceeds_numerical_error": bool(error["max_abs"] > 10.0 * max(float(numerical), 0.0)),
    }


def _clock_record(inputs, occupation):
    energy = np.asarray(inputs["field_energy"], dtype=float)
    normal = np.asarray(inputs["window_normal"], dtype=float)
    return {
        "protocol": inputs["clock_protocol"],
        "same_observer_as_probability": True,
        "canonical_probability": "observer occupation from the original column weights",
        "normal_energy": "episode window_normal shell content; not the occupation",
        "field_energy": "episode field_energy scalar; not the occupation",
        "proper_clock": [float(item) for item in inputs["proper_clock"]],
        "leader_clock": [float(item) for item in inputs["leader_clock"]],
        "field_energy_samples": [float(item) for item in energy],
        "field_energy_change": float(energy[-1] - energy[0]),
        "window_normal_initial": [float(item) for item in normal[0]],
        "window_normal_final": [float(item) for item in normal[-1]],
        "window_normal_change": [float(item) for item in (normal[-1] - normal[0])],
        "leader_initial": float(inputs["leader"][0]),
        "leader_final": float(inputs["leader"][-1]),
        "leader_share_initial": float(inputs["leader_share"][0]),
        "leader_share_final": float(inputs["leader_share"][-1]),
        "occupation_initial": [float(item) for item in occupation[0]],
        "occupation_final": [float(item) for item in occupation[-1]],
        "identified_with_each_other": False,
    }


def _rate_schedule_report(grid, inputs):
    slopes, projection = prolonged_nodal_rates(grid, inputs["coarse_Q_dot"])
    times = inputs["frame_time"]
    values = inputs["frame_Q"]
    midpoints = 0.5 * (times[:-1] + times[1:])
    linear = interpolate_rows(times, values, midpoints)
    finite = hermite_rows(times, values, midpoints)
    if projection["usable"]:
        rated = hermite_with_slopes(times, values, slopes, midpoints)
        nodes = hermite_with_slopes(times, values, slopes, times)
        node_gap = float(np.max(np.abs(nodes - values)))
        minimum = float(np.min(hermite_with_slopes(
            times,
            values,
            slopes,
            np.linspace(float(times[0]), float(times[-1]), 4 * (times.size - 1) + 1),
        )))
    else:
        rated = None
        node_gap = None
        minimum = None
    secant = (values[-1] - values[-2]) / float(times[-1] - times[-2])
    secant_gap = float(np.max(np.abs(pull_geometry(grid, secant) - inputs["coarse_Q_dot"][-1])))
    return slopes, projection, {
        "declared_geometry": SUCCESSOR_INTERPOLATION,
        "nodal_rate_projection": projection,
        "rate_hermite_node_gap": node_gap,
        "rate_hermite_minimum": minimum,
        "rate_hermite_stays_positive": None if minimum is None else bool(minimum > 0.0),
        "midpoint_linear_versus_rate_hermite": None if rated is None else float(np.max(np.abs(linear - rated))),
        "midpoint_linear_versus_finite_difference_hermite": float(np.max(np.abs(linear - finite))),
        "final_secant_versus_saved_rate_max_abs": secant_gap,
        "saved_rate_max_abs": float(np.max(np.abs(inputs["coarse_Q_dot"][-1]))),
        "increment_dt_max": float(np.max(inputs["increment_dt"])),
        "indicator_max_abs": float(np.nanmax(np.abs(inputs["indicator_Q_dot"]))),
        "indicator_role": inputs["indicator_versus_bound"],
        "p_chi_is_dirac_matrix_entry": False,
        "Q_dot_is_dirac_matrix_entry": False,
    }


def code_identities():
    root = _LAB
    paths = {
        "nsc_coupled_local_response.py": root / "src" / "recursive_horizons" / "nsc_coupled_local_response.py",
        "derive_nsc_coupled_local_response.py": root / "scripts" / "derive_nsc_coupled_local_response.py",
        "test_nsc_coupled_local_response.py": root / "tests" / "test_nsc_coupled_local_response.py",
    }
    return {name: sha256_file(path) for name, path in paths.items()}


def _reject_successor(identity, reason, output_json):
    return {
        "schema": SCHEMA_V2,
        "status": "REJECTED_REQUEST_MISMATCH",
        "resumed_complete": False,
        "preserved_record": True,
        "record_path": str(output_json),
        "rejected_identity": identity,
        "reason": reason,
    }


def execute_successor(
    budget_s=SUCCESSOR_PILOT_BUDGET_S,
    confirm_budget_s=SUCCESSOR_CONFIRM_BUDGET_S,
    resume=False,
    output_json=DEFAULT_JSON_V2,
    output_npz=DEFAULT_NPZ_V2,
):
    """Conditional response on the saved ``[0.05, 0.085]`` regeneration geometry.

    The pilot ceiling is 600 process seconds. A named method error above one
    percent of the measured effect may use the 1800-second confirmation
    ceiling. One hundred output nodes are forecast and not run.
    """
    output_json = Path(output_json)
    output_npz = Path(output_npz)
    if output_json.resolve() == DEFAULT_JSON.resolve() or output_npz.resolve() == DEFAULT_NPZ.resolve():
        raise RuntimeError("successor output must not replace the v1 records")
    budget = MeasuredBudget(float(budget_s))
    confirm_ceiling = float(confirm_budget_s)
    if confirm_ceiling < float(budget_s):
        raise ValueError("confirmation ceiling must cover the pilot ceiling")
    inputs = load_successor_inputs()
    identity = successor_identity(inputs, SUCCESSOR_DECLARED_SUBSTEPS)
    if resume and output_json.is_file():
        previous = json.loads(output_json.read_text())
        disposition = successor_resume_disposition(previous, identity)
        if disposition == "return":
            previous["resumed_complete"] = True
            return previous
        if disposition == "reject":
            return _reject_successor(identity, "finished successor bindings differ", output_json)
    grid, factors = build_case_grid({
        "phi0": inputs["phi0"],
        "frame_Q": inputs["frame_Q"],
        "coarse_Q": inputs["coarse_Q"][0],
    })
    observer = inputs["observer"]
    columns = inputs["transported_columns"]
    weights = np.array(inputs["weights"], dtype=float, copy=True)
    times = np.array(inputs["frame_time"], dtype=float, copy=True)
    if not np.array_equal(observer, region_observer(inputs["phi0"], inputs["phi1"], region=0)[0]):
        raise RuntimeError("observer bytes changed after loading")
    action0 = action_residual(grid, factors, inputs["frame_Q"][0], inputs["phi0_frames"][0], inputs["phi1_frames"][0])
    action1 = action_residual(grid, factors, inputs["frame_Q"][-1], inputs["phi0_frames"][-1], inputs["phi1_frames"][-1])
    prolong_gap = float(np.max(np.abs(prolong_geometry(grid, inputs["coarse_Q"][0]) - inputs["frame_Q"][0])))
    if action0["relative"] > ACTION_GATE or action1["relative"] > ACTION_GATE or prolong_gap > PROLONGATION_GATE:
        raise RuntimeError("stored regeneration geometry failed the Galerkin gate")
    slopes, _projection, geometry_report = _rate_schedule_report(grid, inputs)
    autonomous_occupation, autonomous_coherence = autonomous_series(
        observer, inputs["phi0_frames"], inputs["phi1_frames"], weights
    )
    partner_occupation, partner_coherence = autonomous_series(
        observer, inputs["partner_phi0"], inputs["partner_phi1"], weights
    )
    initial_occupation = np.real(np.diag(inputs["initial_blocks"]["parent"]))
    if float(np.max(np.abs(initial_occupation - autonomous_occupation[0]))) > 1e-12:
        raise RuntimeError("initial occupation is not the diagonal of C_AA")
    record = {
        "schema": SCHEMA_V2,
        "status": "RUNNING",
        "request": identity,
        "requested_job_id": SUCCESSOR_JOB,
        "renewal": False,
        "regeneration": False,
        "global_regeneration": False,
        "old_toy_B_calibrated": False,
        "stress_claimed": False,
        "force_claimed": False,
        "gamma_kernel_variation": False,
        "stress_gap": STRESS_GAP,
        "backend": BACKEND,
        "budget_s": float(budget_s),
        "confirm_budget_s": confirm_ceiling,
        "safety": BUDGET_SAFETY,
        "hundred_node_executed": False,
        "phases": {},
        "admissions": [],
    }

    def checkpoint():
        record["cpu_seconds_this_process"] = budget.elapsed()
        write_json(output_json, record)

    def admit(estimate, label, confirmation=False):
        ceiling = confirm_ceiling if confirmation else float(budget_s)
        allowed = budget.elapsed() + BUDGET_SAFETY * float(estimate) < ceiling
        record["admissions"].append({
            "label": label,
            "estimate_s": float(estimate),
            "admitted": bool(allowed),
            "confirmation": bool(confirmation),
            "elapsed_s": budget.elapsed(),
            "ceiling_s": float(ceiling),
        })
        checkpoint()
        return bool(allowed)

    checkpoint()
    probe_steps = SUCCESSOR_PROBE_STEPS
    probe_times = np.array(times[: probe_steps + 1], dtype=float, copy=True)
    print("successor probe", probe_times[0], probe_times[-1], flush=True)
    radial = stored_radial({"frame_time": times, "frame_Q": inputs["frame_Q"]})
    hamiltonian = hamiltonian_from_radial(factors, radial)
    probe = _timed_comparison(hamiltonian, observer, probe_times, columns, weights, substeps=1)
    if not probe["observer_preserved"] or probe["allocation"]["time_indexed_exterior_propagator_bytes"] != 0:
        raise RuntimeError("probe did not keep the fixed observer and the streamed allocation")
    if probe["allocation"]["history_shape"][-1] != 6:
        raise RuntimeError("probe history is not the six source columns")
    probe_effect = _effect_scale(autonomous_occupation[: probe_steps + 1])
    probe_coherence_effect = _effect_scale(autonomous_coherence[: probe_steps + 1])
    record["phases"]["probe"] = {
        "complete": True,
        "steps": probe_steps,
        "substeps": 1,
        "t_start": float(probe_times[0]),
        "t_end": float(probe_times[-1]),
        "streamed_seconds": probe["streamed_seconds"],
        "full_seconds": probe["full_seconds"],
        "occupation_error": probe["occupation_error"],
        "coherence_error": probe["coherence_error"],
        "occupation_effect": probe_effect,
        "coherence_effect": probe_coherence_effect,
        "allocation": probe["allocation"],
        "hamiltonian_builds": int(hamiltonian.stats["builds"]),
        "hamiltonian_hits": int(hamiltonian.stats["hits"]),
        "matrix_bytes": int(hamiltonian.stats["matrix_bytes"]),
    }
    checkpoint()
    available = int(times.size - 1)
    streamed_probe = float(probe["streamed_seconds"])
    full_probe = float(probe["full_seconds"])
    record["forecast"] = {
        "probe_streamed_s": streamed_probe,
        "probe_full_s": full_probe,
        "full_window_substeps_2_s": (
            _forecast_streamed(streamed_probe, probe_steps, available, 2)
            + _forecast_full(full_probe, probe_steps, available, 2)
        ),
        "full_window_substeps_1_s": (
            _forecast_streamed(streamed_probe, probe_steps, available, 1)
            + _forecast_full(full_probe, probe_steps, available, 1)
        ),
        "hundred_output_nodes_streamed_s": _forecast_streamed(streamed_probe, probe_steps, HUNDRED_NODE_OUTPUTS, 1),
        "hundred_node_geometry": "not stored; the forecast is not a run",
        "rule": "streamed replay scales as steps squared times substeps; full-band evolution scales linearly",
    }
    arrays = {
        "probe_times": probe_times,
        "probe_occupation_retained": probe["occupation_retained"],
        "probe_occupation_full": probe["occupation_full"],
        "probe_coherence_retained": probe["coherence_retained"],
        "probe_coherence_full": probe["coherence_full"],
    }
    method = {}
    fine_steps = probe_steps * 2
    fine_times = np.linspace(float(probe_times[0]), float(probe_times[-1]), fine_steps + 1)
    output_estimate = _forecast_streamed(streamed_probe, probe_steps, fine_steps, 1) + _forecast_full(full_probe, probe_steps, fine_steps, 1)
    if admit(output_estimate, "probe_output_sampling"):
        print("successor output sampling", fine_steps, flush=True)
        refined = _timed_comparison(hamiltonian, observer, fine_times, columns, weights, substeps=1)
        movement = float(np.max(np.abs(probe["occupation_full"] - refined["occupation_full"][::2])))
        coherence_movement = float(np.max(np.abs(probe["coherence_full"] - refined["coherence_full"][::2])))
        method["output_sampling"] = {
            "window": [float(probe_times[0]), float(probe_times[-1])],
            "coarse_nodes": int(probe_times.size),
            "fine_nodes": int(fine_times.size),
            "occupation": _indicator_status(movement, probe_effect, probe["occupation_error"]["max_abs"]),
            "coherence": _indicator_status(coherence_movement, probe_coherence_effect, probe["coherence_error"]["max_abs"]),
        }
        arrays["probe_output_occupation_full"] = refined["occupation_full"]
    substep_estimate = _forecast_streamed(streamed_probe, probe_steps, probe_steps, 2) + _forecast_full(full_probe, probe_steps, probe_steps, 2)
    if admit(substep_estimate, "probe_substeps_2"):
        print("successor probe substeps 2", flush=True)
        refined_sub = _timed_comparison(hamiltonian, observer, probe_times, columns, weights, substeps=2)
        movement = float(np.max(np.abs(probe["occupation_full"] - refined_sub["occupation_full"])))
        coherence_movement = float(np.max(np.abs(probe["coherence_full"] - refined_sub["coherence_full"])))
        method["substeps"] = {
            "window": [float(probe_times[0]), float(probe_times[-1])],
            "occupation": _indicator_status(movement, probe_effect, refined_sub["occupation_error"]["max_abs"]),
            "coherence": _indicator_status(coherence_movement, probe_coherence_effect, refined_sub["coherence_error"]["max_abs"]),
        }
    if geometry_report["rate_hermite_stays_positive"] and admit(_forecast_full(full_probe, probe_steps, probe_steps, 1), "probe_rate_hermite"):
        print("successor probe rate hermite", flush=True)
        rate_radial = lambda time: hermite_with_slopes(times, inputs["frame_Q"], slopes, [float(time)])[0]
        rate_hamiltonian = hamiltonian_from_radial(factors, rate_radial)
        rate_occupation, rate_coherence = project_full_occupation(
            rate_hamiltonian, observer, probe_times, columns, weights, substeps=1
        )
        method["rate_hermite"] = {
            "window": [float(probe_times[0]), float(probe_times[-1])],
            "occupation": _indicator_status(
                float(np.max(np.abs(rate_occupation - probe["occupation_full"]))),
                probe_effect,
                probe["occupation_error"]["max_abs"],
            ),
            "coherence": _indicator_status(
                float(np.max(np.abs(rate_coherence - probe["coherence_full"]))),
                probe_coherence_effect,
                probe["coherence_error"]["max_abs"],
            ),
            "slopes": "prolonged saved coarse Q_dot",
        }
        arrays["probe_rate_hermite_occupation"] = rate_occupation
    timestep_movement = float(np.max(np.abs(partner_occupation - autonomous_occupation)))
    method["timestep_autonomous"] = {
        "cases": [inputs["case"], inputs["partner_case"]],
        "phi_gap": inputs["partner_phi_gap"],
        "Q_gap": inputs["partner_Q_gap"],
        "occupation": _indicator_status(
            timestep_movement,
            _effect_scale(autonomous_occupation),
            float(np.max(inputs["initial_blocks"]["hermitian_defect"])),
        ),
        "conditional_rerun": False,
    }
    record["phases"]["method"] = method
    checkpoint()
    named_error = [
        name
        for name, item in method.items()
        for piece in (item.get("occupation"), item.get("coherence"))
        if isinstance(piece, dict) and piece.get("status") == "above_one_percent"
    ]
    record["named_method_error"] = named_error
    confirmation = bool(named_error)
    if "timestep_autonomous" in named_error and admit(
        _forecast_full(full_probe, probe_steps, probe_steps, 1),
        "partner_conditional_probe",
        confirmation=True,
    ):
        print("successor partner conditional probe", flush=True)
        partner_radial = stored_radial({"frame_time": times, "frame_Q": inputs["partner_Q"]})
        partner_hamiltonian = hamiltonian_from_radial(factors, partner_radial)
        partner_conditional, _partner_coherence = project_full_occupation(
            partner_hamiltonian, observer, probe_times, columns, weights, substeps=1
        )
        method["timestep_conditional_probe"] = {
            "occupation": _indicator_status(
                float(np.max(np.abs(partner_conditional - probe["occupation_full"]))),
                probe_effect,
                probe["occupation_error"]["max_abs"],
            ),
            "conditional_rerun": True,
        }
        record["phases"]["method"] = method
        named_error = [
            name
            for name, item in method.items()
            for piece in (item.get("occupation"), item.get("coherence"))
            if isinstance(piece, dict) and piece.get("status") == "above_one_percent"
        ]
        record["named_method_error"] = named_error
        confirmation = bool(named_error)
        arrays["probe_partner_occupation"] = partner_conditional
    print("successor omissions and actual cross", flush=True)
    panel = {}
    for name, flags in (("memory", {"memory": False}), ("outside_drive", {"outside_drive": False})):
        estimate = streamed_probe
        if not admit(estimate, "probe_" + name, confirmation=confirmation and name in named_error):
            continue
        omitted, occupation, _coherence = evolve_streamed(
            hamiltonian, observer, probe_times, columns, weights, **flags
        )
        panel[name] = _omission_against(occupation, probe["occupation_full"], probe["occupation_error"]["max_abs"])
        panel[name]["history_norm"] = float(np.linalg.norm(omitted["history"]))
        panel[name]["initial_exterior_norm"] = float(np.linalg.norm(omitted["initial_exterior"]))
        arrays["probe_occupation_" + name + "_off"] = occupation
    if admit(2.0 * streamed_probe, "probe_cross", confirmation=confirmation):
        cross = column_cross_control(
            hamiltonian, observer, probe_times, columns, weights, substeps=1, baseline=probe["result"]
        )
        panel["initial_cross"] = {
            "separation_max_frobenius": cross["separation_max_frobenius"],
            "occupation_diagonal_movement": cross["occupation_diagonal_movement"],
            "numerical_error": cross["numerical_error"],
            "active": cross["active"],
            "synthetic_six_mode_substituted": False,
            "effect": _indicator_status(cross["occupation_diagonal_movement"], probe_effect, cross["numerical_error"]),
        }
        arrays["probe_cross_occupation_movement"] = np.real(np.diagonal(cross["cross_covariance"], axis1=1, axis2=2))
    else:
        cross = None
    record["phases"]["probe_controls"] = {
        "complete": True,
        "omissions": panel,
        "active": [name for name, item in panel.items() if item.get("exceeds_numerical_error") or item.get("active")],
    }
    checkpoint()
    full_estimate = record["forecast"]["full_window_substeps_2_s"]
    full_comparison = None
    if admit(full_estimate, "full_window_substeps_2"):
        print("successor full window", times[0], times[-1], flush=True)
        full_comparison = _timed_comparison(
            hamiltonian, observer, times, columns, weights, substeps=SUCCESSOR_DECLARED_SUBSTEPS
        )
        if full_comparison["allocation"]["time_indexed_exterior_propagator_bytes"] != 0:
            raise RuntimeError("full window stored a time-indexed exterior propagator")
        arrays["times"] = times
        arrays["occupation_retained"] = full_comparison["occupation_retained"]
        arrays["occupation_full"] = full_comparison["occupation_full"]
        arrays["occupation_autonomous"] = autonomous_occupation
        arrays["coherence_retained"] = full_comparison["coherence_retained"]
        arrays["coherence_full"] = full_comparison["coherence_full"]
        arrays["coherence_autonomous"] = autonomous_coherence
        record["phases"]["full_window"] = {
            "complete": True,
            "substeps": SUCCESSOR_DECLARED_SUBSTEPS,
            "t_start": float(times[0]),
            "t_end": float(times[-1]),
            "steps": available,
            "streamed_seconds": full_comparison["streamed_seconds"],
            "full_seconds": full_comparison["full_seconds"],
            "occupation_error": full_comparison["occupation_error"],
            "coherence_error": full_comparison["coherence_error"],
            "phase_error": full_comparison["phase_error"],
            "allocation": full_comparison["allocation"],
            "history_norm": full_comparison["history_norm"],
            "initial_exterior_norm": full_comparison["initial_exterior_norm"],
            "trapezoid_residual_max": full_comparison["trapezoid_residual_max"],
        }
    else:
        record["phases"]["full_window"] = {"complete": False, "reason": "forecast exceeded the pilot ceiling"}
    checkpoint()
    if full_comparison is not None:
        effect = _effect_scale(autonomous_occupation)
        coherence_effect = _effect_scale(autonomous_coherence)
        numerical = full_comparison["occupation_error"]["max_abs"]
        coherence_numerical = full_comparison["coherence_error"]["max_abs"]
        record["phases"]["full_window"]["reduction_versus_effect"] = {
            "occupation": _indicator_status(numerical, effect, 0.0),
            "coherence": _indicator_status(coherence_numerical, coherence_effect, 0.0),
        }
        record["phases"]["conditional_versus_autonomous"] = {
            "full": signal_error(full_comparison["occupation_full"], autonomous_occupation),
            "retained": signal_error(full_comparison["occupation_retained"], autonomous_occupation),
            "coherence_full": signal_error(full_comparison["coherence_full"], autonomous_coherence),
            "occupation_effect": effect,
            "coherence_effect": coherence_effect,
            "full_fraction": fraction_of_effect(
                signal_error(full_comparison["occupation_full"], autonomous_occupation)["max_abs"],
                effect,
            ),
        }
        substep_estimate = _forecast_streamed(streamed_probe, probe_steps, available, 1) + _forecast_full(full_probe, probe_steps, available, 1)
        if admit(substep_estimate, "full_substeps_1"):
            print("successor full substeps 1", flush=True)
            coarse_full_occupation, coarse_full_coherence = project_full_occupation(
                hamiltonian, observer, times, columns, weights, substeps=1
            )
            method["full_substeps"] = {
                "occupation": _indicator_status(
                    float(np.max(np.abs(coarse_full_occupation - full_comparison["occupation_full"]))),
                    effect,
                    numerical,
                ),
                "coherence": _indicator_status(
                    float(np.max(np.abs(coarse_full_coherence - full_comparison["coherence_full"]))),
                    coherence_effect,
                    coherence_numerical,
                ),
            }
            arrays["occupation_full_substeps_1"] = coarse_full_occupation
        hermite_estimate = _forecast_full(full_probe, probe_steps, available, SUCCESSOR_DECLARED_SUBSTEPS)
        run_hermite = geometry_report["rate_hermite_stays_positive"] and admit(hermite_estimate, "full_rate_hermite")
        if geometry_report["rate_hermite_stays_positive"] and not run_hermite and confirmation:
            run_hermite = admit(hermite_estimate, "full_rate_hermite_confirmation", confirmation=True)
        if run_hermite:
            print("successor full rate hermite", flush=True)
            rate_radial = lambda time: hermite_with_slopes(times, inputs["frame_Q"], slopes, [float(time)])[0]
            rate_hamiltonian = hamiltonian_from_radial(factors, rate_radial)
            rate_occupation, rate_coherence = project_full_occupation(
                rate_hamiltonian,
                observer,
                times,
                columns,
                weights,
                substeps=SUCCESSOR_DECLARED_SUBSTEPS,
            )
            method["full_rate_hermite"] = {
                "occupation": _indicator_status(
                    float(np.max(np.abs(rate_occupation - full_comparison["occupation_full"]))),
                    effect,
                    numerical,
                ),
                "coherence": _indicator_status(
                    float(np.max(np.abs(rate_coherence - full_comparison["coherence_full"]))),
                    coherence_effect,
                    coherence_numerical,
                ),
                "against_autonomous": signal_error(rate_occupation, autonomous_occupation),
            }
            arrays["occupation_rate_hermite"] = rate_occupation
            arrays["coherence_rate_hermite"] = rate_coherence
        full_panel = {}
        for name, flags in (("memory", {"memory": False}), ("outside_drive", {"outside_drive": False})):
            estimate = _forecast_streamed(streamed_probe, probe_steps, available, SUCCESSOR_DECLARED_SUBSTEPS)
            if not admit(estimate, "full_" + name):
                continue
            print("successor full omission", name, flush=True)
            omitted, occupation, _coherence = evolve_streamed(
                hamiltonian,
                observer,
                times,
                columns,
                weights,
                propagator_substeps=SUCCESSOR_DECLARED_SUBSTEPS,
                **flags,
            )
            full_panel[name] = _omission_against(occupation, full_comparison["occupation_full"], numerical)
            full_panel[name]["history_norm"] = float(np.linalg.norm(omitted["history"]))
            full_panel[name]["initial_exterior_norm"] = float(np.linalg.norm(omitted["initial_exterior"]))
            arrays["occupation_" + name + "_off"] = occupation
        if admit(2.0 * _forecast_streamed(streamed_probe, probe_steps, available, SUCCESSOR_DECLARED_SUBSTEPS), "full_cross"):
            print("successor full cross", flush=True)
            full_cross = column_cross_control(
                hamiltonian,
                observer,
                times,
                columns,
                weights,
                substeps=SUCCESSOR_DECLARED_SUBSTEPS,
                baseline=full_comparison["result"],
            )
            full_panel["initial_cross"] = {
                "separation_max_frobenius": full_cross["separation_max_frobenius"],
                "occupation_diagonal_movement": full_cross["occupation_diagonal_movement"],
                "numerical_error": full_cross["numerical_error"],
                "active": full_cross["active"],
                "synthetic_six_mode_substituted": False,
                "effect": _indicator_status(full_cross["occupation_diagonal_movement"], effect, full_cross["numerical_error"]),
            }
            arrays["covariance_retained"] = full_cross["covariance_total"]
            arrays["covariance_without_cross"] = full_cross["covariance_without_cross"]
            arrays["covariance_eigenvalues"] = np.linalg.eigvalsh(full_cross["covariance_total"]).real
        else:
            full_cross = None
        record["phases"]["full_controls"] = {
            "complete": True,
            "omissions": full_panel,
            "active": [
                name for name, item in full_panel.items() if item.get("exceeds_numerical_error") or item.get("active")
            ],
            "inactive": [
                name for name, item in full_panel.items() if not (item.get("exceeds_numerical_error") or item.get("active"))
            ],
        }
        record["phases"]["method"] = method
    named_error = [
        name
        for name, item in record["phases"].get("method", {}).items()
        for piece in (item.get("occupation"), item.get("coherence"))
        if isinstance(piece, dict) and piece.get("status") == "above_one_percent"
    ]
    record["named_method_error"] = named_error
    record["method_error_within_one_percent"] = not bool(named_error)
    record["finite_domain"] = {
        "case": inputs["case"],
        "partner_case": inputs["partner_case"],
        "fermion_band": int(factors["nf"]),
        "quadrature": int(factors["nq"]),
        "length": float(factors["length"]),
        "spinor_dimension": int(2 * factors["nf"]),
        "observer": "T=0 initial mode columns 0 and 1, no phase QR; not the Phi(0.05) basis",
        "weights": [float(item) for item in weights],
        "initial_state": "bitwise transported Phi(0.05) from the spherical episode handoff",
        "time_window": [float(times[0]), float(times[-1])],
        "stored_frames": int(times.size),
        "geometry": SUCCESSOR_INTERPOLATION,
        "reference": "independent full-band midpoint evolution on H(Q(t)), compared with stored autonomous Phi frames",
        "reduction": "nsc_evolving_reduction.evolve_retained_region backend='streamed'",
        "regeneration_declared": False,
        "production_rk4_rerun": False,
        "stress": False,
    }
    record["initial_blocks"] = _block_record(inputs["initial_blocks"])
    record["geometry"] = geometry_report
    record["clock"] = _clock_record(inputs, autonomous_occupation)
    record["episode_citation"] = {
        "protocol_id": "nsc-regeneration-episode-v1",
        "requested_job_id": SUCCESSOR_JOB,
        "stop_reason": inputs["episode_stop_reason"],
        "attained_T": inputs["episode_attained_T"],
        "arc_x_start": inputs["arc_x_start"],
        "arc_x_end": inputs["arc_x_end"],
        "structure_maintained": inputs["structure_maintained"],
        "balanced_flow": inputs["balanced_flow"],
        "renewed": inputs["structure_renewed"],
        "joined_renewal_proxy": inputs["joined_renewal_proxy"],
        "old_toy_B_not_used": inputs["old_toy_B_not_used"],
        "state_reset": inputs["state_reset"],
        "new_force_added": inputs["new_force_added"],
        "episode_goal_is_not_this_reduction": True,
    }
    record["bindings"] = {
        "names": inputs["binding_names"],
        "values": [float(item) for item in inputs["binding_values"]],
        "observer_id": inputs["observer_id"],
        "weights_id": inputs["weights_id"],
        "initial_state_id": inputs["initial_state_id"],
        "geometry_id": inputs["geometry_id"],
        "episode_npz_sha256": inputs["episode_sha256"],
        "episode_json_sha256": inputs["episode_json_sha256"],
        "origin_npz_sha256": inputs["origin_sha256"],
        "origin_json_sha256": inputs["origin_json_sha256"],
        "handoff_bitwise": True,
        "observer_rephased": False,
    }
    record["operator_gate"] = {
        "initial_action_relative": action0["relative"],
        "final_action_relative": action1["relative"],
        "initial_prolongation_gap": prolong_gap,
        "origin_gram_defect": inputs["origin_gram_defect"],
        "fine_identity_image": False,
    }
    record["immutable_v1"] = {
        "json_sha256": sha256_file(DEFAULT_JSON),
        "npz_sha256": sha256_file(DEFAULT_NPZ),
    }
    record["code_identities"] = code_identities()
    record["hamiltonian_builds"] = int(hamiltonian.stats["builds"])
    record["hamiltonian_hits"] = int(hamiltonian.stats["hits"])
    record["cache_key"] = "float.hex of t; a repeated sample is the same matrix object"
    if record["phases"].get("full_window", {}).get("complete") and record["phases"].get("full_controls", {}).get("complete"):
        record["status"] = "MEASURED_SUCCESSOR_RESPONSE"
        record["stop_reason"] = None
    else:
        record["status"] = "BUDGET_EXCEEDED"
        record["stop_reason"] = "a required full-window phase was not admitted"
    record["cpu_seconds_this_process"] = budget.elapsed()
    record["confirmation_ceiling_used"] = bool(any(item["confirmation"] and item["admitted"] for item in record["admissions"]))
    write_npz(output_npz, arrays)
    record["payload_sha256"] = {"npz": sha256_file(output_npz)}
    record["payload_bytes"] = int(output_json.stat().st_size + output_npz.stat().st_size)
    write_json(output_json, record)
    record["payload_bytes"] = int(output_json.stat().st_size + output_npz.stat().st_size)
    record["payload_within_64MiB"] = bool(record["payload_bytes"] <= PAYLOAD_LIMIT_BYTES)
    write_json(output_json, record)
    if not record["payload_within_64MiB"]:
        raise RuntimeError(f"payload {record['payload_bytes']} exceeds 64MiB")
    print("successor", record["status"], "cpu", record["cpu_seconds_this_process"], flush=True)
    return record


def _same_token(left, right):
    if isinstance(left, bool) or isinstance(right, bool):
        return isinstance(left, bool) and isinstance(right, bool) and bool(left) is bool(right)
    if isinstance(left, (int, float, np.floating)) and not isinstance(left, bool):
        if isinstance(right, (int, float, np.floating)) and not isinstance(right, bool):
            return math.isclose(float(left), float(right), rel_tol=1e-12, abs_tol=1e-12)
    return left == right


def successor_consistency_errors(record, arrays, inputs):
    """Typed checks. This function does not write and does not evolve."""
    errors = []
    if not isinstance(record, dict) or record.get("schema") != SCHEMA_V2:
        return ["schema"]
    identity = successor_identity(inputs, SUCCESSOR_DECLARED_SUBSTEPS)
    if record.get("request") != identity:
        errors.append("request")
    if record.get("renewal") is not False or record.get("global_regeneration") is not False:
        errors.append("regeneration_claim")
    if record.get("old_toy_B_calibrated") is not False or record.get("stress_claimed") is not False:
        errors.append("forbidden_claim")
    if record.get("force_claimed") is not False or record.get("gamma_kernel_variation") is not False:
        errors.append("force_claim")
    if record.get("hundred_node_executed") is not False:
        errors.append("hundred_node")
    blocks = record.get("initial_blocks") or {}
    fresh = _block_record(inputs["initial_blocks"])
    for key in ("C_AE_frobenius", "C_EE_frobenius", "admissible", "cross_active_as_a_block"):
        if not _same_token(blocks.get(key), fresh[key]):
            errors.append("initial_" + key)
    if not fresh["cross_active_as_a_block"]:
        errors.append("actual_cross_inactive")
    allocation = ((record.get("phases") or {}).get("full_window") or {}).get("allocation") or {}
    if allocation and int(allocation.get("time_indexed_exterior_propagator_bytes", -1)) != 0:
        errors.append("dense_propagator")
    if arrays is not None and "occupation_autonomous" in arrays and "occupation_full" in arrays:
        autonomous = np.asarray(arrays["occupation_autonomous"])
        conditional = np.asarray(arrays["occupation_full"])
        retained = np.asarray(arrays["occupation_retained"])
        reported = (record.get("phases") or {}).get("conditional_versus_autonomous") or {}
        full_gap = signal_error(conditional, autonomous)["max_abs"]
        retained_gap = signal_error(retained, autonomous)["max_abs"]
        if not _same_token(reported.get("full", {}).get("max_abs"), full_gap):
            errors.append("autonomous_gap")
        if not _same_token(reported.get("retained", {}).get("max_abs"), retained_gap):
            errors.append("retained_autonomous_gap")
        if "occupation_outside_drive_off" in arrays:
            drive = signal_error(arrays["occupation_outside_drive_off"], conditional)
            stored = (((record.get("phases") or {}).get("full_controls") or {}).get("omissions") or {}).get("outside_drive")
            stored_gap = None if not stored else stored["occupation_error_against_reference"]["max_abs"]
            if not _same_token(stored_gap, drive["max_abs"]):
                errors.append("drive_omission")
        if "covariance_retained" in arrays and "covariance_without_cross" in arrays:
            separation = _matrix_series_gap(arrays["covariance_retained"], arrays["covariance_without_cross"])
            stored_cross = (((record.get("phases") or {}).get("full_controls") or {}).get("omissions") or {}).get("initial_cross")
            if not stored_cross or stored_cross.get("synthetic_six_mode_substituted") is not False:
                errors.append("cross_control")
            elif not _same_token(stored_cross.get("separation_max_frobenius"), separation):
                errors.append("cross_separation")
    identities = record.get("code_identities") or {}
    for name, digest in code_identities().items():
        if identities.get(name) != digest:
            errors.append("code:" + name)
    immutable = record.get("immutable_v1") or {}
    if immutable.get("json_sha256") != sha256_file(DEFAULT_JSON) or immutable.get("npz_sha256") != sha256_file(DEFAULT_NPZ):
        errors.append("immutable_v1")
    return errors


def verify_successor(output_json=DEFAULT_JSON_V2, output_npz=DEFAULT_NPZ_V2):
    """Read the successor. Inputs are hashed and no output file is written."""
    output_json = Path(output_json)
    output_npz = Path(output_npz)
    watched = (
        output_json,
        output_npz,
        DEFAULT_JSON,
        DEFAULT_NPZ,
        REGENERATION_JSON,
        REGENERATION_NPZ,
        EPISODE_NPZ,
    )
    before = {path: sha256_file(path) for path in watched}
    record = json.loads(output_json.read_text())
    with np.load(output_npz, allow_pickle=False) as stored:
        arrays = {key: np.array(stored[key]) for key in stored.files}
    inputs = load_successor_inputs()
    errors = successor_consistency_errors(record, arrays, inputs)
    if record.get("payload_sha256", {}).get("npz") != sha256_file(output_npz):
        errors.append("payload_npz_sha")
    after = {path: sha256_file(path) for path in watched}
    if before != after:
        errors.append("readonly")
    if errors:
        raise RuntimeError("successor check failed: " + ", ".join(errors))
    return record


def reject_linear_operator_hamiltonian():
    """Evidence for the dense-callback gap. The reducer is not modified."""
    dimension = 4

    def hamiltonian(_time):
        return LinearOperator(
            (dimension, dimension),
            matvec=lambda vector: vector,
            dtype=np.complex128,
        )

    observer = np.zeros((dimension, 2), dtype=np.complex128)
    observer[0, 0] = 1.0
    observer[1, 1] = 1.0
    columns = np.eye(dimension, 2, dtype=np.complex128)
    try:
        evolve_retained_region(
            hamiltonian,
            observer,
            np.linspace(0.0, 0.1, 3),
            columns,
            backend=BACKEND,
        )
    except (TypeError, ValueError) as error:
        return type(error).__name__
    raise RuntimeError("a LinearOperator Hamiltonian was accepted")
