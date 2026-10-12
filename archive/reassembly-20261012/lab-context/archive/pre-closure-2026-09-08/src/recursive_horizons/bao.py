"""A small, dependency-free DESI DR2 BAO distance-likelihood reproducer.

The calculation is deliberately narrow: it profiles a flat, late-time,
constant-``w`` background against the public 13-point Gaussian BAO data vector.
It omits radiation (sub-percent over the supplied ``z <= 2.33`` range), early
universe sound-horizon physics, CMB and supernova data, full-shape information,
and every proposed external-domain mechanism.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from math import isfinite, sqrt
from numbers import Real
from pathlib import Path
from typing import Final, Iterable, Sequence


C_KM_S: Final[float] = 299_792.458
ALLOWED_OBSERVABLES: Final[frozenset[str]] = frozenset({"DM_over_rs", "DH_over_rs", "DV_over_rs"})
OMEGA_M_BOUNDS: Final[tuple[float, float]] = (0.05, 0.7)
W_BOUNDS: Final[tuple[float, float]] = (-2.0, 0.0)
PROFILE_68_DELTA_CHI2: Final[float] = 1.0
PROFILE_95_DELTA_CHI2: Final[float] = 3.841458820694124


@dataclass(frozen=True, slots=True)
class BAOPoint:
    """One public BAO measurement in the order used by its covariance matrix."""

    redshift: float
    value: float
    observable: str


@dataclass(frozen=True, slots=True)
class BAODataset:
    """Validated 13-point BAO vector and covariance."""

    points: tuple[BAOPoint, ...]
    covariance: tuple[tuple[float, ...], ...]


@dataclass(frozen=True, slots=True)
class ProfilePoint:
    """The analytic-amplitude profile at fixed ``(Omega_m, w)``."""

    omega_m: float
    w: float
    h0_rd_km_s: float
    amplitude_c_over_h0rd: float
    chi2: float


def _finite(name: str, value: Real) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive(name: str, value: Real) -> float:
    value = _finite(name, value)
    if value <= 0.0:
        raise ValueError(f"{name} must be positive")
    return value


def default_data_directory() -> Path:
    """Return the repository-local public DESI DR2 BAO input directory."""

    return Path(__file__).resolve().parents[2] / "collected-data" / "desi-dr2-bao"


def _cholesky(matrix: Sequence[Sequence[float]]) -> tuple[tuple[float, ...], ...]:
    """Return lower Cholesky factor or reject a non-positive-definite matrix."""

    n = len(matrix)
    if n == 0 or any(len(row) != n for row in matrix):
        raise ValueError("covariance must be a non-empty square matrix")
    lower = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1):
            value = matrix[i][j]
            for k in range(j):
                value -= lower[i][k] * lower[j][k]
            if i == j:
                if not isfinite(value) or value <= 0.0:
                    raise ValueError("covariance must be positive definite")
                lower[i][j] = sqrt(value)
            else:
                lower[i][j] = value / lower[j][j]
    return tuple(tuple(row) for row in lower)


def _solve_cholesky(lower: Sequence[Sequence[float]], vector: Sequence[float]) -> tuple[float, ...]:
    """Solve ``C x=b`` using the supplied lower Cholesky factor of ``C``."""

    n = len(lower)
    if len(vector) != n:
        raise ValueError("linear-solve vector length does not match covariance")
    forward = [0.0] * n
    for i in range(n):
        forward[i] = (vector[i] - sum(lower[i][j] * forward[j] for j in range(i))) / lower[i][i]
    solution = [0.0] * n
    for i in range(n - 1, -1, -1):
        solution[i] = (forward[i] - sum(lower[j][i] * solution[j] for j in range(i + 1, n))) / lower[i][i]
    return tuple(solution)


def _dot(left: Iterable[float], right: Iterable[float]) -> float:
    return sum(a * b for a, b in zip(left, right, strict=True))


def load_bao_dataset(mean_path: Path, covariance_path: Path) -> BAODataset:
    """Parse and strictly validate the DESI public mean vector and covariance."""

    raw_points: list[BAOPoint] = []
    for line_number, line in enumerate(mean_path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        fields = stripped.split()
        if len(fields) != 3:
            raise ValueError(f"mean line {line_number} must have redshift, value, observable")
        redshift = _positive(f"mean redshift at line {line_number}", float(fields[0]))
        value = _positive(f"mean value at line {line_number}", float(fields[1]))
        observable = fields[2]
        if observable not in ALLOWED_OBSERVABLES:
            raise ValueError(f"unsupported observable at line {line_number}: {observable}")
        raw_points.append(BAOPoint(redshift, value, observable))
    if len(raw_points) != 13:
        raise ValueError(f"DESI DR2 BAO mean vector must contain 13 points, found {len(raw_points)}")

    rows: list[tuple[float, ...]] = []
    for line_number, line in enumerate(covariance_path.read_text(encoding="utf-8").splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        row = tuple(_finite(f"covariance entry at line {line_number}", float(item)) for item in stripped.split())
        rows.append(row)
    if len(rows) != 13 or any(len(row) != 13 for row in rows):
        raise ValueError("DESI DR2 BAO covariance must be exactly 13 by 13")
    for i in range(13):
        for j in range(13):
            if abs(rows[i][j] - rows[j][i]) > 1.0e-12:
                raise ValueError("covariance must be symmetric")
    _cholesky(rows)
    return BAODataset(tuple(raw_points), tuple(rows))


def e_squared(redshift: float, omega_m: float, w: float) -> float:
    """Return late-time flat constant-w ``E(z)^2`` with radiation omitted."""

    redshift = _finite("redshift", redshift)
    if redshift < 0.0:
        raise ValueError("redshift must be non-negative")
    omega_m = _finite("omega_m", omega_m)
    w = _finite("w", w)
    if not 0.0 < omega_m < 1.0:
        raise ValueError("omega_m must lie strictly between zero and one")
    one_plus_z = 1.0 + redshift
    result = omega_m * one_plus_z**3 + (1.0 - omega_m) * one_plus_z ** (3.0 * (1.0 + w))
    if not isfinite(result) or result <= 0.0:
        raise ValueError("E(z)^2 must be finite and positive")
    return result


def inverse_e_integral(redshift: float, omega_m: float, w: float, intervals: int = 256) -> float:
    """Integrate ``int_0^z dz'/E(z')`` with deterministic composite Simpson."""

    redshift = _positive("redshift", redshift)
    if not isinstance(intervals, int) or intervals < 2 or intervals % 2:
        raise ValueError("intervals must be an even integer of at least two")
    spacing = redshift / intervals
    total = 0.0
    for index in range(intervals + 1):
        z = index * spacing
        inverse_e = 1.0 / sqrt(e_squared(z, omega_m, w))
        weight = 1.0 if index in (0, intervals) else (4.0 if index % 2 else 2.0)
        total += weight * inverse_e
    result = total * spacing / 3.0
    if not isfinite(result) or result <= 0.0:
        raise ValueError("distance integral must be finite and positive")
    return result


def bao_shape_vector(dataset: BAODataset, omega_m: float, w: float, intervals: int = 256) -> tuple[float, ...]:
    """Return model shapes that multiply ``A=c/(H0*r_d)`` for each datum."""

    integrals: dict[float, float] = {}
    shapes: list[float] = []
    for point in dataset.points:
        if point.redshift not in integrals:
            integrals[point.redshift] = inverse_e_integral(point.redshift, omega_m, w, intervals)
        integral = integrals[point.redshift]
        inverse_e = 1.0 / sqrt(e_squared(point.redshift, omega_m, w))
        if point.observable == "DM_over_rs":
            shape = integral
        elif point.observable == "DH_over_rs":
            shape = inverse_e
        else:
            shape = (point.redshift * integral * integral * inverse_e) ** (1.0 / 3.0)
        shapes.append(_positive("BAO shape", shape))
    return tuple(shapes)


def analytic_amplitude_profile(dataset: BAODataset, omega_m: float, w: float, intervals: int = 256) -> ProfilePoint:
    """Profile the positive nuisance amplitude ``A=c/(H0*r_d)`` analytically."""

    shapes = bao_shape_vector(dataset, omega_m, w, intervals)
    lower = _cholesky(dataset.covariance)
    observations = tuple(point.value for point in dataset.points)
    inverse_covariance_shapes = _solve_cholesky(lower, shapes)
    inverse_covariance_observations = _solve_cholesky(lower, observations)
    denominator = _dot(shapes, inverse_covariance_shapes)
    amplitude = _dot(shapes, inverse_covariance_observations) / denominator
    if not isfinite(amplitude) or amplitude <= 0.0:
        raise ValueError("analytic BAO amplitude must be finite and positive")
    residual = tuple(observation - amplitude * shape for observation, shape in zip(observations, shapes, strict=True))
    chi2 = _dot(residual, _solve_cholesky(lower, residual))
    if not isfinite(chi2) or chi2 < 0.0:
        raise ValueError("profile chi-square must be finite and non-negative")
    return ProfilePoint(omega_m, w, C_KM_S / amplitude, amplitude, chi2)


def golden_minimize(function: object, lower: float, upper: float, tolerance: float = 1.0e-7) -> tuple[float, float]:
    """Deterministically minimize a unimodal scalar function on a closed interval."""

    if not lower < upper or tolerance <= 0.0:
        raise ValueError("golden search requires ordered bounds and positive tolerance")
    evaluator = function  # Keeps the public function dependency-free and duck typed.
    ratio = (sqrt(5.0) - 1.0) / 2.0
    left, right = lower, upper
    x1, x2 = right - ratio * (right - left), left + ratio * (right - left)
    f1, f2 = evaluator(x1), evaluator(x2)  # type: ignore[operator]
    while right - left > tolerance:
        if f1 <= f2:
            right, x2, f2 = x2, x1, f1
            x1 = right - ratio * (right - left)
            f1 = evaluator(x1)  # type: ignore[operator]
        else:
            left, x1, f1 = x1, x2, f2
            x2 = left + ratio * (right - left)
            f2 = evaluator(x2)  # type: ignore[operator]
    position = 0.5 * (left + right)
    return position, evaluator(position)  # type: ignore[operator]


def profile_at_w(dataset: BAODataset, w: float, intervals: int = 256) -> ProfilePoint:
    """Profile ``Omega_m`` in the declared interval at a fixed ``w``."""

    w = _finite("w", w)
    omega_m, _ = golden_minimize(
        lambda value: analytic_amplitude_profile(dataset, value, w, intervals).chi2,
        *OMEGA_M_BOUNDS,
    )
    return analytic_amplitude_profile(dataset, omega_m, w, intervals)


def best_fit_profile(dataset: BAODataset, intervals: int = 256) -> ProfilePoint:
    """Nested profile over declared ``w`` and ``Omega_m`` bounds."""

    w, _ = golden_minimize(lambda value: profile_at_w(dataset, value, intervals).chi2, *W_BOUNDS)
    return profile_at_w(dataset, w, intervals)


def _profile_crossing(dataset: BAODataset, target: float, left: float, right: float, intervals: int) -> float:
    """Bisect a bracketed one-dimensional profile-likelihood crossing."""

    f_left = profile_at_w(dataset, left, intervals).chi2 - target
    f_right = profile_at_w(dataset, right, intervals).chi2 - target
    if f_left * f_right > 0.0:
        raise ValueError("profile interval is not bracketed by declared bounds")
    for _ in range(48):
        midpoint = 0.5 * (left + right)
        f_midpoint = profile_at_w(dataset, midpoint, intervals).chi2 - target
        if f_left * f_midpoint <= 0.0:
            right, f_right = midpoint, f_midpoint
        else:
            left, f_left = midpoint, f_midpoint
    return 0.5 * (left + right)


def profile_interval(dataset: BAODataset, best: ProfilePoint, delta_chi2: float, intervals: int = 2048) -> tuple[float, float]:
    """Return a one-parameter profile interval using a fixed Delta-chi-square."""

    target = best.chi2 + _positive("delta_chi2", delta_chi2)
    return (
        _profile_crossing(dataset, target, W_BOUNDS[0], best.w, intervals),
        _profile_crossing(dataset, target, best.w, W_BOUNDS[1], intervals),
    )


def file_record(path: Path) -> dict[str, object]:
    """Return a small provenance record for a tracked public input."""

    content = path.read_bytes()
    return {"filename": path.name, "bytes": len(content), "sha256": sha256(content).hexdigest()}


def reproduce_desi_dr2_bao(data_directory: Path | None = None, intervals: int = 256) -> dict[str, object]:
    """Compute the intentionally limited DESI DR2 BAO-only profile record."""

    directory = default_data_directory() if data_directory is None else Path(data_directory)
    mean_path = directory / "desi_gaussian_bao_ALL_GCcomb_mean.txt"
    covariance_path = directory / "desi_gaussian_bao_ALL_GCcomb_cov.txt"
    dataset = load_bao_dataset(mean_path, covariance_path)
    best = best_fit_profile(dataset, intervals)
    candidates = {str(w): profile_at_w(dataset, w, intervals) for w in (-1.0, -2.0 / 3.0, -1.0 / 3.0)}
    interval_68 = profile_interval(dataset, best, PROFILE_68_DELTA_CHI2, intervals)
    interval_95 = profile_interval(dataset, best, PROFILE_95_DELTA_CHI2, intervals)
    coarse = best_fit_profile(dataset, intervals // 2)
    input_files = tuple(
        directory / name
        for name in (
            mean_path.name,
            covariance_path.name,
            "chain.margestats",
            "bestfit.minimum.txt",
            "chain.updated.yaml",
        )
    )
    return {
        "schema_version": 1,
        "project_version": "0.11.0",
        "artifact": "DESI_DR2_BAO_only_flat_constant_w_profile",
        "scope": "Real 13-point Gaussian BAO distance likelihood only; not DESI+CMB/SN, not full-shape, not an external-origin inference.",
        "official_likelihood_source": {
            "repository_commit": "b7b8a36e9bccb063081f811f323cada21ab5fbdd",
            "configuration": "chain.updated.yaml",
            "mean_vector": mean_path.name,
            "covariance": covariance_path.name,
        },
        "input_files": [file_record(path) for path in input_files],
        "model": {
            "E_squared": "Omega_m*(1+z)^3 + (1-Omega_m)*(1+z)^(3*(1+w))",
            "radiation": "omitted; sub-percent late-time approximation over z<=2.33",
            "flat": True,
            "constant_w": True,
            "nuisance": "A=c/(H0*r_d), profiled analytically with A>0",
            "c_km_s": C_KM_S,
        },
        "declared_search": {
            "omega_m_bounds": list(OMEGA_M_BOUNDS),
            "w_bounds": list(W_BOUNDS),
            "golden_tolerance": 1.0e-7,
            "simpson_intervals": intervals,
            "one_parameter_delta_chi2": {"68_percent": PROFILE_68_DELTA_CHI2, "95_percent": PROFILE_95_DELTA_CHI2},
        },
        "data_points": [
            {"z": point.redshift, "value": point.value, "observable": point.observable}
            for point in dataset.points
        ],
        "best_fit": _profile_record(best),
        "dof": len(dataset.points) - 3,
        "fixed_w_candidates": {
            label: {**_profile_record(candidate), "delta_chi2_from_best": candidate.chi2 - best.chi2}
            for label, candidate in candidates.items()
        },
        "one_parameter_w_profile_intervals": {
            "68_percent": list(interval_68),
            "95_percent": list(interval_95),
        },
        "official_cross_check": {
            "MAP_from_bestfit_minimum": {"w": -0.91129473, "omega_m": 0.29740107, "H0rd_km_s": 9974.6091, "chi2_BAO": 9.0335199},
            "margestats": {"w_68": "-0.916 +/- 0.078", "omega_m_68": "0.2969 +/- 0.0089", "H0rd_km_s_68": "9984 +160 -180"},
            "interpretation": "Cross-check only: chain.updated.yaml fixes H0 and samples hrdrag, so H0*r_d is effectively free there. The official posterior nevertheless uses its full CAMB radiation/neutrino implementation and stated priors/MCMC sampling; this artifact is a late-time BAO-distance profile.",
        },
        "convergence_refinement": {
            "coarse_simpson_intervals": intervals // 2,
            "fine_simpson_intervals": intervals,
            "abs_best_w_difference": abs(best.w - coarse.w),
            "abs_best_omega_m_difference": abs(best.omega_m - coarse.omega_m),
            "abs_best_chi2_difference": abs(best.chi2 - coarse.chi2),
        },
        "nonclaims": [
            "Not a full DESI+CMB or DESI+supernova likelihood.",
            "Not a DESI full-shape calculation.",
            "Does not determine H0 without an independently modeled sound horizon.",
            "Does not establish an external or boundary origin for dark energy.",
            "Does not apply the separately conserved DEB-1 identity to a flux-fed Q != 0 component.",
        ],
    }


def _profile_record(profile: ProfilePoint) -> dict[str, float]:
    return {
        "w": profile.w,
        "omega_m": profile.omega_m,
        "H0rd_km_s": profile.h0_rd_km_s,
        "amplitude_c_over_H0rd": profile.amplitude_c_over_h0rd,
        "chi2": profile.chi2,
    }
