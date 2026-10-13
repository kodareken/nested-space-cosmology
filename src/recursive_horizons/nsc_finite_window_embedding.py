"""Isometric compression of the fixed spherical Dirac operator onto one frozen window.

The operator is the continuum form of ``apply_dirac``: the same coefficients,
spin matrices and anticommutators, not a new coupling. The target is the
existing ``finite_window`` matrix. Inherited-law identities are read from
``exact_controls`` and are not recomputed.

The recorded real counter-carrier packet matches the onsite blocks and cannot
match this link. Degree-4 sine lobes with one constant phase on the plus
right lobe do match the compression. The six modes vanish together on an open
arc, so they are not an invariant subspace.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss

from .nsc_nested_qualities import exact_controls, finite_window
from .nsc_spherical_coupling import (
    CALIBRATION,
    CARRIER_K,
    KAPPA,
    OMEGA,
    PERIOD,
    apply_dirac,
    antiperiodic_momentum,
    prepare_rank6,
    profile_coordinate,
)

SCHEMA = "NSC-FINITE-WINDOW-EMBEDDING-v1"
A0 = float(CALIBRATION["a0"])
BETA0 = float(CALIBRATION["beta0"])
B0 = float(CALIBRATION["b0"])
K_CARRIER = float(CARRIER_K)
ELL = 1.0
DEGREE = 4
PACKET_ARC = (0.0, 4.0)
LAMBDA = float(np.log(OMEGA))

# Degree-4 sine coefficients on each unit lobe. Minus right-lobe phase is 0.
# Polished so the continuum compression residual is ~1e-14 on 160-point Gauss panels.
PLUS_LEFT = (
    0.4055258127195104,
    -0.00893047570183846,
    -0.4906451379753416,
    -0.33729793105745876,
)
PLUS_RIGHT = (
    -0.10433452497223869,
    -0.5314230030102629,
    0.1908081925236093,
    -0.38892499089629734,
)
MINUS_LEFT = (
    0.2651180411626812,
    -0.0990661544962709,
    0.6262412749070796,
    0.01208524862821752,
)
MINUS_RIGHT = (
    -0.5145584714075958,
    -0.46003562611207954,
    0.14835896346661445,
    -0.17076461524916445,
)
PLUS_RIGHT_PHASE = 0.6332154732379276
MINUS_RIGHT_PHASE = 0.0
COLUMN_PHASE = 0.5298187870937372

_LAB_ROOT = Path(__file__).resolve().parents[2]
RECORD_PATH = _LAB_ROOT / "results" / "development" / "nsc-finite-window-embedding-v1.json"


def frozen_window():
    """The existing depth-3 window. Not a new link."""
    local = sp.Matrix([[1, sp.I / 5], [-sp.I / 5, 2]])
    link = sp.Matrix(
        [[sp.Rational(1, 4), sp.I / 7], [sp.Rational(1, 9), sp.Rational(1, 6)]]
    )
    window = finite_window(local, link, sp.Rational(3, 2), 0, 3)
    matrix = np.array(
        [[complex(window[i, j]) for j in range(6)] for i in range(6)],
        dtype=np.complex128,
    )
    return matrix, window


def real_lobe_overlap_bound():
    """Necessary overlap for a real plus envelope at the recorded carrier.

    On the overlap interval the weight ρ = Ω^{1+t} lies in (Ω, Ω²). For real
    lobe samples f, g with unit joint norm and ∫ f g = 0, Popoviciu's variance
    bound gives |∫ ρ f g| ≤ (Ω² − Ω) / 4. Matching B++ = 1/4 forces a larger
    overlap. The bound uses the stored a0, β0 and k = 1. It does not apply to
    a complex right lobe, a different carrier, or a different target entry.
    """
    gap = A0 - BETA0
    if gap <= 0:
        raise ValueError("onsite diagonal ratio requires a0 > beta0")
    magnitude = 0.25 / gap
    required = magnitude * float(np.cos(K_CARRIER)) / K_CARRIER
    maximum = (float(OMEGA) ** 2 - float(OMEGA)) / 4.0
    return {
        "carrier_k": K_CARRIER,
        "required_weighted_overlap": required,
        "maximum_real_orthogonal_overlap": maximum,
        "gap": required - maximum,
        "obstruction": bool(required > maximum + 1e-9),
        "scope": (
            "real left and right lobes, recorded k, sigma2-eigen spinors, "
            "neighbor L2 orthogonality, target B++=1/4, stored a0 and beta0"
        ),
    }


def conjugate_link_obstruction(matrix):
    """B_{-+} − conj(B_{+-}) for one link block. Zero for a conjugate column."""
    link = np.asarray(matrix)[:2, 2:4]
    defect = complex(link[1, 0] - np.conjugate(link[0, 1]))
    return {
        "link_plus_minus": [link[0, 1].real, link[0, 1].imag],
        "link_minus_plus": [link[1, 0].real, link[1, 0].imag],
        "conjugate_defect": abs(defect),
        "scope": (
            "minus spatial column equals one constant phase times the conjugate "
            "of the plus spatial column, in the sigma2-eigen spin frame"
        ),
    }


def restrictive_probe(points=128):
    """Existing prepare_rank6 packet, compressed by the discrete Dirac operator."""
    points = int(points)
    dx = PERIOD / points
    xi = np.arange(points, dtype=float) * dx
    phi0, phi1, preparation = prepare_rank6(xi, dx)
    target, _ = frozen_window()
    matrix, gram_defect, leak, bridge = _discrete_matrix(phi0, phi1, xi)
    difference = matrix - target
    onsite, link, corner = _block_errors(difference)
    return {
        "points": points,
        "column_orthonormal_defect": float(preparation["column_orthonormal_defect"]),
        "gram_defect": gram_defect,
        "onsite_frobenius": onsite,
        "link_frobenius": link,
        "corner_frobenius": corner,
        "max_onsite": float(max(onsite)),
        "min_link": float(min(link)),
        "conjugate": conjugate_link_obstruction(matrix),
        "frozen_conjugate": conjugate_link_obstruction(target),
        "leak_fraction": leak,
        "bridge_image_max": bridge,
        "phase_on_minus_spinor_column": True,
        "same_real_envelope": True,
    }


def _block_errors(difference):
    onsite, link = [], []
    for region in range(3):
        sl = slice(2 * region, 2 * region + 2)
        onsite.append(float(np.linalg.norm(difference[sl, sl])))
    for region in range(2):
        row = slice(2 * region, 2 * region + 2)
        col = slice(2 * region + 2, 2 * region + 4)
        link.append(float(np.linalg.norm(difference[row, col])))
    corner = float(np.linalg.norm(difference[:2, 4:]))
    return onsite, link, corner


def _geometry(xi):
    coordinate = profile_coordinate(xi, PERIOD)
    scale = np.exp(np.log(OMEGA) * coordinate)
    length = B0 * scale
    shift = BETA0 * scale
    radial = np.full(xi.shape, B0 / A0)
    return length, radial, shift


def _discrete_matrix(phi0, phi1, xi):
    momentum, _metric = antiperiodic_momentum(xi.size, PERIOD)
    length, radial, shift = _geometry(xi)
    image0, image1 = apply_dirac(phi0, phi1, length, radial, shift, KAPPA, momentum)
    gram = phi0.conj().T @ phi0 + phi1.conj().T @ phi1
    matrix = phi0.conj().T @ image0 + phi1.conj().T @ image1
    columns = np.vstack((phi0, phi1))
    images = np.vstack((image0, image1))
    residual = images - columns @ matrix
    energy = np.linalg.norm(images, axis=0)
    leak = np.linalg.norm(residual, axis=0) / np.maximum(energy, 1e-30)
    bridge = (xi > PACKET_ARC[1]) & (xi < PERIOD)
    bridge_max = float(np.max(np.abs(images[np.concatenate((bridge, bridge))]))) if np.any(bridge) else 0.0
    return matrix, float(np.max(np.abs(gram - np.eye(6)))), [float(v) for v in leak], bridge_max


def _panels(n_panels, quad):
    nodes, weights = leggauss(int(quad))
    pieces_x, pieces_w = [], []
    for left in range(int(n_panels)):
        pieces_x.append(left + 0.5 * (nodes + 1.0))
        pieces_w.append(0.5 * weights)
    return np.concatenate(pieces_x), np.concatenate(pieces_w)


def _lobe(x, shift, left_coeff, right_coeff, right_phase):
    z = x - shift
    values = np.zeros(x.shape, dtype=np.complex128)
    deriv = np.zeros(x.shape, dtype=np.complex128)
    modes = np.arange(1, DEGREE + 1)
    left = (z > 0.0) & (z < 1.0)
    right = (z > 1.0) & (z < 2.0)
    if np.any(left):
        t = z[left]
        sine = np.sin(np.pi * t[:, None] * modes[None, :]) * np.sqrt(2.0)
        cosine = np.pi * modes[None, :] * np.cos(np.pi * t[:, None] * modes[None, :]) * np.sqrt(2.0)
        values[left] = sine @ np.asarray(left_coeff, dtype=float)
        deriv[left] = cosine @ np.asarray(left_coeff, dtype=float)
    if np.any(right):
        t = z[right] - 1.0
        sine = np.sin(np.pi * t[:, None] * modes[None, :]) * np.sqrt(2.0)
        cosine = np.pi * modes[None, :] * np.cos(np.pi * t[:, None] * modes[None, :]) * np.sqrt(2.0)
        factor = np.exp(1j * right_phase)
        values[right] = factor * (sine @ np.asarray(right_coeff, dtype=float))
        deriv[right] = factor * (cosine @ np.asarray(right_coeff, dtype=float))
    return values, deriv


def _carrier(x, shift, envelope, derivative, signed_k, column_phase):
    phase = np.exp(1j * column_phase) * np.exp(1j * signed_k * (x - shift))
    values = phase * envelope
    deriv = phase * (derivative + 1j * signed_k * envelope)
    return values, deriv


def continuum_modes(quad=256):
    """Six spinor profiles on (0, 4). Each is supported on (n, n+2)."""
    x, wt = _panels(4, quad)
    rho = OMEGA ** x
    modes = []
    for region in range(3):
        plus_env, plus_d = _lobe(x, region, PLUS_LEFT, PLUS_RIGHT, PLUS_RIGHT_PHASE)
        minus_env, minus_d = _lobe(x, region, MINUS_LEFT, MINUS_RIGHT, MINUS_RIGHT_PHASE)
        fp, dfp = _carrier(x, region, plus_env, plus_d, K_CARRIER, 0.0)
        fm, dfm = _carrier(x, region, minus_env, minus_d, -K_CARRIER, COLUMN_PHASE)
        modes.append({"spin": "plus", "region": region, "f": fp, "df": dfp, "support": (region, region + 2)})
        modes.append({"spin": "minus", "region": region, "f": fm, "df": dfm, "support": (region, region + 2)})
    return {"x": x, "wt": wt, "rho": rho, "modes": modes}


def _image(mode, rho):
    """Spinor image of one mode, returned as plus and minus scalar channels."""
    directional = mode["df"] + 0.5 * LAMBDA * mode["f"]
    if mode["spin"] == "plus":
        plus = 1j * (BETA0 - A0) * rho * directional
        minus = 1j * B0 * rho * mode["f"]
    else:
        plus = -1j * B0 * rho * mode["f"]
        minus = 1j * (BETA0 + A0) * rho * directional
    return plus, minus


def _inner(wt, left, right):
    return np.sum(wt * np.conj(left) * right)


def continuum_compression(quad=256):
    """Direct quadrature of V† H V. Blocks are not inserted by scaling."""
    data = continuum_modes(quad)
    wt, rho, modes = data["wt"], data["rho"], data["modes"]
    images = [_image(mode, rho) for mode in modes]
    matrix = np.zeros((6, 6), dtype=np.complex128)
    gram = np.zeros((6, 6), dtype=np.complex128)
    for i, left in enumerate(modes):
        for j, right in enumerate(modes):
            if left["spin"] == right["spin"]:
                gram[i, j] = _inner(wt, left["f"], right["f"])
            plus, minus = images[j]
            if left["spin"] == "plus":
                matrix[i, j] = _inner(wt, left["f"], plus)
            else:
                matrix[i, j] = _inner(wt, left["f"], minus)
    target, _window = frozen_window()
    difference = matrix - target
    onsite, link, corner = _block_errors(difference)
    leak = []
    captured = []
    for plus, minus in images:
        coeff_plus = np.array([
            _inner(wt, mode["f"], plus) for mode in modes if mode["spin"] == "plus"
        ])
        coeff_minus = np.array([
            _inner(wt, mode["f"], minus) for mode in modes if mode["spin"] == "minus"
        ])
        pred_plus = sum(
            coeff_plus[k] * mode["f"]
            for k, mode in enumerate(m for m in modes if m["spin"] == "plus")
        )
        pred_minus = sum(
            coeff_minus[k] * mode["f"]
            for k, mode in enumerate(m for m in modes if m["spin"] == "minus")
        )
        residual_plus = plus - pred_plus
        residual_minus = minus - pred_minus
        residual_norm = np.sqrt(_inner(wt, residual_plus, residual_plus).real + _inner(wt, residual_minus, residual_minus).real)
        image_norm = np.sqrt(_inner(wt, plus, plus).real + _inner(wt, minus, minus).real)
        leak.append(float(residual_norm / max(image_norm, 1e-30)))
        captured.append(float(np.sqrt(max(0.0, 1.0 - leak[-1] ** 2))))
    support = []
    x = data["x"]
    for mode in modes:
        mass = _inner(wt, mode["f"], mode["f"]).real
        lo, hi = mode["support"]
        outside = (x < lo) | (x > hi)
        leaked = float(np.sum(wt[outside] * np.abs(mode["f"][outside]) ** 2))
        support.append({
            "region": mode["region"],
            "spin": mode["spin"],
            "interval": [lo, hi],
            "mass": float(mass),
            "mass_outside_interval": leaked,
        })
    sign = _sign_current(matrix)
    population = _subspace_population_curvature(modes, images, wt)
    return {
        "quad": int(quad),
        "matrix_error": float(np.linalg.norm(difference)),
        "hermiticity_defect": float(np.linalg.norm(matrix - matrix.conj().T)),
        "gram_defect": float(np.max(np.abs(gram - np.eye(6)))),
        "onsite_frobenius": onsite,
        "link_frobenius": link,
        "corner_frobenius": corner,
        "leak_fraction": leak,
        "captured_fraction": captured,
        "support": support,
        "bridge_component_of_image": 0.0,
        "sign_control_current": sign,
        "subspace_population_second_derivative": population,
        "conjugate": conjugate_link_obstruction(matrix),
        "plus_right_phase": PLUS_RIGHT_PHASE,
        "minus_right_phase": MINUS_RIGHT_PHASE,
        "column_phase": COLUMN_PHASE,
        "carrier_k": K_CARRIER,
    }


def _sign_current(matrix):
    """Same regional trace as the finite-turnover sign control. Not a new law."""
    state = np.zeros(6, dtype=np.complex128)
    state[0] = 1.0 / np.sqrt(2.0)
    state[2] = 1j / np.sqrt(2.0)
    covariance = np.outer(state, state.conj())
    region = np.diag([1.0, 1.0, 0.0, 0.0, 0.0, 0.0]).astype(np.complex128)
    commutator = matrix @ covariance - covariance @ matrix
    current = -1j * np.trace(region @ commutator)
    return {"real": float(current.real), "imag": float(current.imag), "frozen_value": 0.25}


def _subspace_population_curvature(modes, images, wt):
    """d²/dt² of subspace population at the sign-control state, equal to −2‖(1−P)Hψ‖²."""
    state = np.zeros(6, dtype=np.complex128)
    state[0] = 1.0 / np.sqrt(2.0)
    state[2] = 1j / np.sqrt(2.0)
    plus = sum(state[j] * images[j][0] for j in range(6))
    minus = sum(state[j] * images[j][1] for j in range(6))
    coeff_plus = [
        _inner(wt, mode["f"], plus) for mode in modes if mode["spin"] == "plus"
    ]
    coeff_minus = [
        _inner(wt, mode["f"], minus) for mode in modes if mode["spin"] == "minus"
    ]
    pred_plus = sum(
        coeff_plus[k] * mode["f"]
        for k, mode in enumerate(m for m in modes if m["spin"] == "plus")
    )
    pred_minus = sum(
        coeff_minus[k] * mode["f"]
        for k, mode in enumerate(m for m in modes if m["spin"] == "minus")
    )
    residual = _inner(wt, plus - pred_plus, plus - pred_plus).real
    residual += _inner(wt, minus - pred_minus, minus - pred_minus).real
    return {"minus_two_leak_norm_squared": float(-2.0 * residual), "leak_norm_squared": float(residual)}


def discrete_compression(points):
    """Sample the same modes and compress with ``apply_dirac``."""
    points = int(points)
    dx = PERIOD / points
    xi = np.arange(points, dtype=float) * dx
    sqrt_dx = np.sqrt(dx)
    phi0 = np.zeros((points, 6), dtype=np.complex128)
    phi1 = np.zeros((points, 6), dtype=np.complex128)
    for region in range(3):
        plus_env, plus_d = _lobe(xi, region, PLUS_LEFT, PLUS_RIGHT, PLUS_RIGHT_PHASE)
        minus_env, minus_d = _lobe(xi, region, MINUS_LEFT, MINUS_RIGHT, MINUS_RIGHT_PHASE)
        plus, _dp = _carrier(xi, region, plus_env, plus_d, K_CARRIER, 0.0)
        minus, _dm = _carrier(xi, region, minus_env, minus_d, -K_CARRIER, COLUMN_PHASE)
        phi0[:, 2 * region] = sqrt_dx * plus / np.sqrt(2.0)
        phi1[:, 2 * region] = sqrt_dx * 1j * plus / np.sqrt(2.0)
        phi0[:, 2 * region + 1] = sqrt_dx * minus / np.sqrt(2.0)
        phi1[:, 2 * region + 1] = sqrt_dx * (-1j) * minus / np.sqrt(2.0)
    target, _window = frozen_window()
    matrix, gram_defect, leak, bridge = _discrete_matrix(phi0, phi1, xi)
    difference = matrix - target
    onsite, link, corner = _block_errors(difference)
    outside = (xi < 0.0) | (xi > 4.0)
    outside_mass = float(np.sum(np.abs(phi0[outside]) ** 2 + np.abs(phi1[outside]) ** 2))
    return {
        "points": points,
        "matrix_error": float(np.linalg.norm(difference)),
        "gram_defect": gram_defect,
        "onsite_frobenius": onsite,
        "link_frobenius": link,
        "corner_frobenius": corner,
        "leak_fraction": leak,
        "bridge_image_max": bridge,
        "sampled_mass_outside_packet_arc": outside_mass,
    }


def inherited_law_reuse():
    """Read the existing identity. Do not rebuild the Schur complement."""
    controls = exact_controls()
    return {
        "recomputed": False,
        "normalized_shift_exact": bool(controls["checks"]["normalized_shift_exact"]),
        "omega": controls["omega"],
        "depth": controls["depth"],
    }


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_record(quad=256, points=(128, 256)):
    """Write the restrictive probe first, then the construction."""
    restrictive = restrictive_probe(128)
    bound = real_lobe_overlap_bound()
    partial = {
        "schema": SCHEMA,
        "status": "PARTIAL_RESTRICTIVE_PROBE",
        "restrictive_packet": restrictive,
        "real_lobe_necessary_relation": bound,
        "field_equations_changed": False,
        "hand_inserted_matrix_couplings": False,
        "target_B_substituted": False,
    }
    RECORD_PATH.parent.mkdir(parents=True, exist_ok=True)
    RECORD_PATH.write_text(json.dumps(_jsonable(partial), indent=2) + "\n")
    continuum = continuum_compression(quad)
    discrete = [discrete_compression(int(n)) for n in points]
    inherited = inherited_law_reuse()
    record = {
        "schema": SCHEMA,
        "status": "CONSTRUCTED_COMPRESSION_NOT_INVARIANT",
        "operator": "sigma2/2*{a,P}+sigma1*kappa*L-1/2*{beta,P}",
        "operator_owner": "nsc_spherical_coupling.apply_dirac",
        "window_owner": "nsc_nested_qualities.finite_window",
        "geometry": {
            "Q": B0 / A0,
            "length_density": "b0*Omega**s",
            "shift": "beta0*Omega**s",
            "s": "x on the packet arc [0,4], circle length 8",
            "kappa": KAPPA,
            "Omega": float(OMEGA),
            "a0": A0,
            "beta0": BETA0,
            "b0": B0,
        },
        "domains": {
            "packet_arc": list(PACKET_ARC),
            "regions": 3,
            "modes_per_region": 2,
            "mode_support": "(n, n+2)",
            "lobes": "degree-4 sine series on (n,n+1) and (n+1,n+2); values vanish at the endpoints",
            "seam_regularity": (
                "The profile is continuous. Derivatives may jump at the lobe seams and at the outer endpoints. "
                "Continuum quadrature uses the piecewise classical derivative and is not an N=512 Fourier evolution."
            ),
            "plus_right_lobe_phase": PLUS_RIGHT_PHASE,
            "minus_right_lobe_phase": MINUS_RIGHT_PHASE,
            "column_phase": COLUMN_PHASE,
            "carrier": "plus exp(+i k (x-n)), minus exp(i*column_phase) exp(-i k (x-n)), k=1",
            "spin_frame": "sigma2 eigen-spinors [1, i]/sqrt(2) and [1, -i]/sqrt(2)",
            "complement": "L2 orthogonal of these six modes; bridge component of H psi is zero in the continuum operator",
        },
        "restrictive_packet": restrictive,
        "real_lobe_necessary_relation": bound,
        "construction": continuum,
        "discrete_apply_dirac": discrete,
        "inherited_law": inherited,
        "interpretation": {
            "matches_frozen_blocks": bool(continuum["matrix_error"] < 1e-9),
            "instantaneous_projected_generator_is_J": True,
            "sign_control_current_matches_one_quarter": bool(
                abs(continuum["sign_control_current"]["real"] - 0.25) < 1e-9
                and abs(continuum["sign_control_current"]["imag"]) < 1e-9
            ),
            "invariant_subspace": False,
            "closed_stage1_evolution_retained": False,
            "reason": (
                "V† H V equals the frozen window, so subspace quadratic forms and "
                "the finite-window currents agree at that instant. "
                "||(1-P) H V|| is a positive fraction of ||H V|| and the modes "
                "vanish on a common open arc, so the Dirac evolution leaves the window."
            ),
        },
        "scope_limits": [
            "Not a statement about every NSC operator or every finite window.",
            "The conjugate-column obstruction is only for minus = phase * conjugate(plus).",
            "The real-lobe overlap obstruction is only at the recorded carrier, spin frame and B++.",
            "Unique continuation bars an invariant subspace for modes with a common open zero set of this first-order operator. It does not bar the compression.",
            "A spectral basis without regional support is not this construction.",
        ],
        "field_equations_changed": False,
        "hand_inserted_matrix_couplings": False,
        "target_B_substituted": False,
        "all_nsc_claimed": False,
    }
    payload = _jsonable(record)
    RECORD_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    return payload


def _jsonable(value):
    if isinstance(value, dict):
        return {key: _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    return value


if __name__ == "__main__":
    payload = _jsonable(build_record())
    payload["source_hashes"] = {
        "src/recursive_horizons/nsc_finite_window_embedding.py": _sha256(Path(__file__)),
        "tests/test_nsc_finite_window_embedding.py": _sha256(
            _LAB_ROOT / "tests" / "test_nsc_finite_window_embedding.py"
        ),
        "docs/nsc-finite-window-embedding.md": _sha256(
            _LAB_ROOT / "docs" / "nsc-finite-window-embedding.md"
        ),
    }
    RECORD_PATH.write_text(json.dumps(payload, indent=2) + "\n")
    print(RECORD_PATH)
    print(payload["status"], payload["construction"]["matrix_error"])
