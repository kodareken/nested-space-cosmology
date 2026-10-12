"""Initial three-packet compression and its geometry-graded window.

The finite-window algebra is the nested-qualities theorem. The packets, Löwdin
step, gauge and Dirac operator are the active coupling and Galerkin owners.
Blocks are the measured quadratic form. The historical rational link is not a
fitting target. No trajectory is integrated.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
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
from numpy.polynomial.legendre import leggauss

from .nsc_conformal_adm_source import direct_hamiltonian
from .nsc_covariant_operator import CovariantStaticMetric
from .nsc_nested_qualities import exact_controls
from .nsc_spherical_coupling import (
    CALIBRATION,
    CARRIER_K,
    ELL,
    KAPPA,
    OMEGA,
    PACKET_COUNT,
    PERIOD,
    _bump,
    _lobe_norms,
    _sample_lobe,
    apply_dirac,
    lowdin,
    prepare_rank6,
    profile_coordinate,
)
from .nsc_spherical_feedback_action import first_order_density
from .nsc_spherical_galerkin_coupling import (
    blank_state,
    build_grid,
    load_physical_columns,
    preparation_problems,
    prolong_columns,
    prolong_geometry,
    prolong_state,
)

SCHEMA = "NSC-GEOMETRY-GRADED-WINDOW-v1"
STATUS = "INITIAL_GRADED_COMPRESSION_NOT_INVARIANT"
CPU_LIMIT_S = 20.0
DEPTH = 3
BLOCK = 2
RESOLVENT_Z = 1.0 + 2.0j
LEAK_CLOSURE_MAX = 0.1

_LAB = Path(__file__).resolve().parents[2]
RECORD_PATH = _LAB / "results" / "development" / "nsc-geometry-graded-window-v1.json"
EPISODE_JSON = _LAB / "results" / "development" / "nsc-spherical-feedback-episode-v1.json"
EPISODE_NPZ = _LAB / "results" / "development" / "nsc-spherical-feedback-episode-v1.npz"
EPISODE_RUN = "nf512_dt_0_0005"

# Optional historical example. Compared after the compression. Not a target.
_EXAMPLE_H = np.array([[1.0, 1j / 5], [-1j / 5, 2.0]], dtype=np.complex128)
_EXAMPLE_B = np.array(
    [[0.25, 1j / 7], [1.0 / 9.0, 1.0 / 6.0]],
    dtype=np.complex128,
)


def _maxabs(matrix):
    return float(np.max(np.abs(np.asarray(matrix))))


def _frobenius(matrix):
    return float(np.linalg.norm(np.asarray(matrix)))


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _complex_record(matrix):
    array = np.asarray(matrix, dtype=np.complex128)
    return {"real": np.real(array).tolist(), "imag": np.imag(array).tolist()}


def _pair_blocks(matrix):
    """Diagonal, nearest link and corner of a depth-3 two-component window."""
    matrix = np.asarray(matrix, dtype=np.complex128)
    diagonals = [matrix[2 * n:2 * n + 2, 2 * n:2 * n + 2] for n in range(3)]
    links = [matrix[2 * n:2 * n + 2, 2 * n + 2:2 * n + 4] for n in range(2)]
    corner = matrix[:2, 4:6]
    return diagonals, links, corner


def assemble_window(local, link, omega, depth=DEPTH):
    """J_nn = omega^n local and J_n,n+1 = omega^n link. Same block rule as the theorem."""
    local = np.asarray(local, dtype=np.complex128)
    link = np.asarray(link, dtype=np.complex128)
    if local.shape != (BLOCK, BLOCK) or link.shape != (BLOCK, BLOCK):
        raise ValueError("regional blocks are 2 by 2")
    if not np.isfinite(omega) or omega <= 0:
        raise ValueError("positive Omega required")
    size = BLOCK * depth
    window = np.zeros((size, size), dtype=np.complex128)
    for region in range(depth):
        sl = slice(BLOCK * region, BLOCK * region + BLOCK)
        window[sl, sl] = omega ** region * local
        if region + 1 < depth:
            nxt = slice(BLOCK * region + BLOCK, BLOCK * region + 2 * BLOCK)
            window[sl, nxt] = omega ** region * link
            window[nxt, sl] = omega ** region * link.conj().T
    return window


def congruence_normalized(window, omega, depth=DEPTH):
    """Remove Omega^{n/2} from each regional frame. Link becomes Omega^{-1/2} B01."""
    scales = np.repeat(np.power(omega, 0.5 * np.arange(depth)), BLOCK)
    inverse = 1.0 / scales
    matrix = np.asarray(window, dtype=np.complex128)
    return inverse[:, None] * matrix * inverse[None, :]


def schur_top(local, link, omega, z, depth=DEPTH):
    """Terminal elimination. Returns S_0 for a window starting at label 0."""
    local = np.asarray(local, dtype=np.complex128)
    link = np.asarray(link, dtype=np.complex128)
    eye = np.eye(BLOCK, dtype=np.complex128)
    response = z * eye - (omega ** (depth - 1)) * local
    for region in range(depth - 2, -1, -1):
        coupling = (omega ** (2 * region)) * (link @ np.linalg.solve(response, link.conj().T))
        response = z * eye - (omega ** region) * local - coupling
    return response


def gamma_top(local, link, omega, z, depth=DEPTH, include_omega=True):
    """Normalized recurrence. include_omega=False drops the theorem's 1/Omega."""
    local = np.asarray(local, dtype=np.complex128)
    link = np.asarray(link, dtype=np.complex128)
    eye = np.eye(BLOCK, dtype=np.complex128)
    gamma = (z / omega ** (depth - 1)) * eye - local
    weight = (1.0 / omega) if include_omega else 1.0
    for region in range(depth - 2, -1, -1):
        coupling = weight * (link @ np.linalg.solve(gamma, link.conj().T))
        gamma = (z / omega ** region) * eye - local - coupling
    return gamma


def grading_report(matrix, omega):
    """Compare blocks with one H and one B01. Corner is the skipped region pair."""
    diagonals, links, corner = _pair_blocks(matrix)
    local = diagonals[0]
    link = links[0]
    return {
        "H": local,
        "B01": link,
        "diagonal_scale_1": _maxabs(diagonals[1] - omega * local),
        "diagonal_scale_2": _maxabs(diagonals[2] - omega ** 2 * local),
        "link_scale_1": _maxabs(links[1] - omega * link),
        "corner_max": _maxabs(corner),
        "corner_frobenius": _frobenius(corner),
        "link_frobenius": _frobenius(link),
        "onsite_frobenius": _frobenius(local),
        "hermiticity": _maxabs(matrix - matrix.conj().T),
    }


def resolvent_report(local, link, omega, measured):
    """Projected resolvent against the theorem's elimination, on these derived blocks."""
    window = assemble_window(local, link, omega)
    eye = np.eye(window.shape[0], dtype=np.complex128)
    resolvent = np.linalg.inv(RESOLVENT_Z * eye - window)
    top = schur_top(local, link, omega, RESOLVENT_Z)
    gamma = gamma_top(local, link, omega, RESOLVENT_Z, include_omega=True)
    omitted = gamma_top(local, link, omega, RESOLVENT_Z, include_omega=False)
    normalized = link / np.sqrt(omega)
    substituted = schur_top(local, normalized, omega, RESOLVENT_Z)
    congruence = congruence_normalized(measured, omega)
    link_gaps = [
        _maxabs(congruence[2 * n:2 * n + 2, 2 * n + 2:2 * n + 4] - normalized)
        for n in range(2)
    ]
    diagonal_gaps = [
        _maxabs(congruence[2 * n:2 * n + 2, 2 * n:2 * n + 2] - local)
        for n in range(3)
    ]
    # Direction of a non-Hermitian link. Uses the same terminal block as the last Schur step.
    terminal = RESOLVENT_Z * np.eye(BLOCK) - (omega ** (DEPTH - 1)) * local
    forward = link @ np.linalg.solve(terminal, link.conj().T)
    reversed_link = link.conj().T @ np.linalg.solve(terminal, link)
    return {
        "assembly_max": _maxabs(window - measured),
        "schur_matches_resolvent": _maxabs(np.linalg.inv(top) - resolvent[:BLOCK, :BLOCK]),
        "gamma_matches_schur": _maxabs(gamma - top),
        "omitted_omega_factor": _maxabs(omitted - top),
        "normalized_B_substituted_as_block_coefficient": _maxabs(
            np.linalg.inv(substituted) - resolvent[:BLOCK, :BLOCK]
        ),
        "reversed_link_on_terminal_block": _maxabs(forward - reversed_link),
        "congruence_diagonal_max": float(max(diagonal_gaps)),
        "congruence_link_max": float(max(link_gaps)),
        "normalized_B": normalized,
    }


def _fft_derivative(values, spacing):
    count = values.shape[0]
    wavenumber = 2.0 * np.pi * np.fft.fftfreq(count, d=spacing)
    return np.fft.ifft(1j * wavenumber * np.fft.fft(values, axis=0), axis=0)


def _align(values, samples):
    """Broadcast a nodal coefficient across spinor or packet columns."""
    values = np.asarray(values)
    while values.ndim < np.ndim(samples):
        values = values[..., None]
    return values


def _anticommutator_closed(values, samples, derivative, log_omega):
    values = _align(values, samples)
    return -1j * (values * derivative + 0.5 * log_omega * values * samples)


def _anticommutator_direct(values, samples, spacing):
    """{f,P}/2 by the product rule. P is an FFT derivative, not the closed Omega form."""
    values = _align(values, samples)
    momentum = -1j * _fft_derivative(samples, spacing)
    product = -1j * _fft_derivative(values * samples, spacing)
    return 0.5 * (values * momentum + product)


def _spinor_image(length, radial, shift, kappa, upper, lower, anti):
    """Owner block. anti(coefficient, samples) is {coefficient, P}/2 acting on samples."""
    slope = length / radial
    mass = _align(kappa * length, lower)
    image_upper = -1j * anti(slope, lower) + mass * lower - anti(shift, upper)
    image_lower = 1j * anti(slope, upper) + mass * upper - anti(shift, lower)
    return image_upper, image_lower


def _gaussian_spinor(x, center, width):
    z = (x - center) / width
    envelope = np.exp(-(z ** 2))
    upper = envelope
    lower = 0.2 * z * envelope
    d_env = envelope * (-2.0 * z / width)
    d_upper = d_env
    d_lower = 0.2 * (envelope / width + z * d_env)
    return upper, lower, d_upper, d_lower


def translation_algebra():
    """H T_1 = Omega T_1 H when every coefficient scales as Omega^x and Q is constant.

    Quadrature uses the closed anticommutator on Gauss nodes. The direct
    operator uses the FFT product rule on a uniform period and does not insert
    the logarithmic derivative by hand. Both see a translate by one.
    """
    omega = float(OMEGA)
    log_omega = float(np.log(omega))
    a0 = float(CALIBRATION["a0"])
    b0 = float(CALIBRATION["b0"])
    beta0 = float(CALIBRATION["beta0"])
    q0 = b0 / a0
    center = 0.5
    width = 0.08

    nodes, weights = leggauss(96)
    left, right = 0.15, 0.85
    x = left + 0.5 * (right - left) * (nodes + 1.0)
    w = 0.5 * (right - left) * weights
    upper, lower, d_upper, d_lower = _gaussian_spinor(x, center, width)

    def closed_at(coordinate, upper_s, lower_s, d_upper_s, d_lower_s):
        length = b0 * omega ** coordinate
        radial = np.full(coordinate.shape, q0)
        shift = beta0 * omega ** coordinate
        derivative_of = {id(upper_s): d_upper_s, id(lower_s): d_lower_s}

        def anti(values, samples):
            return _anticommutator_closed(values, samples, derivative_of[id(samples)], log_omega)

        return _spinor_image(length, radial, shift, KAPPA, upper_s, lower_s, anti)

    image = closed_at(x, upper, lower, d_upper, d_lower)
    translated = closed_at(x + 1.0, upper, lower, d_upper, d_lower)
    residual = (
        np.abs(translated[0] - omega * image[0]) ** 2
        + np.abs(translated[1] - omega * image[1]) ** 2
    )
    norm = np.abs(omega * image[0]) ** 2 + np.abs(omega * image[1]) ** 2
    quadrature_relative = float(np.sqrt(np.sum(w * residual) / np.sum(w * norm)))

    count = 1024
    length_period = 4.0
    spacing = length_period / count
    grid = np.arange(count) * spacing
    g_upper, g_lower, g_d_upper, g_d_lower = _gaussian_spinor(grid, center, width)
    length = b0 * omega ** grid
    radial = np.full(count, q0)
    shift = beta0 * omega ** grid

    def direct_at(upper_s, lower_s):
        def anti(values, samples):
            return _anticommutator_direct(values, samples, spacing)

        return _spinor_image(length, radial, shift, KAPPA, upper_s, lower_s, anti)

    direct_image = direct_at(g_upper, g_lower)
    derivative_of = {id(g_upper): g_d_upper, id(g_lower): g_d_lower}
    closed_image = _spinor_image(
        length, radial, shift, KAPPA, g_upper, g_lower,
        lambda values, samples: _anticommutator_closed(
            values, samples, derivative_of[id(samples)], log_omega,
        ),
    )
    mask = (grid > 0.2) & (grid < 0.8)
    closed_norm = np.max(np.abs(closed_image[0][mask])) + np.max(np.abs(closed_image[1][mask]))
    operator_gap = (
        np.max(np.abs(direct_image[0][mask] - closed_image[0][mask]))
        + np.max(np.abs(direct_image[1][mask] - closed_image[1][mask]))
    ) / closed_norm
    shift_count = int(round(1.0 / spacing))
    moved_upper = np.roll(g_upper, shift_count)
    moved_lower = np.roll(g_lower, shift_count)
    moved_image = direct_at(moved_upper, moved_lower)
    moved_mask = np.roll(mask, shift_count)
    direct_residual = (
        np.max(np.abs(moved_image[0][moved_mask] - omega * np.roll(direct_image[0], shift_count)[moved_mask]))
        + np.max(np.abs(moved_image[1][moved_mask] - omega * np.roll(direct_image[1], shift_count)[moved_mask]))
    )
    direct_scale = (
        np.max(np.abs(omega * np.roll(direct_image[0], shift_count)[moved_mask]))
        + np.max(np.abs(omega * np.roll(direct_image[1], shift_count)[moved_mask]))
    )
    return {
        "identity": "H T_1 = Omega T_1 H for coefficients proportional to Omega^x and constant Q",
        "quadrature_relative": quadrature_relative,
        "direct_operator_relative": float(direct_residual / direct_scale),
        "direct_versus_closed_relative": float(operator_gap),
        "quadrature_nodes": 96,
        "direct_points": count,
        "direct_period": length_period,
        "gaussian_center": center,
        "gaussian_width": width,
        "kappa": int(KAPPA),
        "Q": q0,
    }


def _bump_derivative(u):
    """Derivative of the owner bump. Tails below exp(-700) are the zero extension."""
    u = np.asarray(u, dtype=float)
    value = np.zeros(u.shape, dtype=float)
    derivative = np.zeros(u.shape, dtype=float)
    mask = (u > 1.0e-6) & (u < 1.0 - 1.0e-6)
    uu = u[mask]
    span = uu * (1.0 - uu)
    keep = span > 1.0 / 700.0
    if not np.any(keep):
        return value, derivative
    host = np.flatnonzero(mask)[keep]
    span = span[keep]
    uu = uu[keep]
    bump = np.exp(-1.0 / span)
    value[host] = bump
    derivative[host] = bump * (1.0 - 2.0 * uu) / span ** 2
    return value, derivative


def _continuum_lobes(x, weights, region, norms):
    even, d_even = _bump_derivative(x - region)
    odd_u = x - (region + 1.0)
    odd_b, odd_db = _bump_derivative(odd_u)
    odd = (odd_u - 0.5) * odd_b
    d_odd = odd_b + (odd_u - 0.5) * odd_db
    even = even / norms[0]
    d_even = d_even / norms[0]
    odd = odd / norms[1]
    d_odd = d_odd / norms[1]
    even_norm = np.sqrt(np.sum(weights * even ** 2))
    odd_norm = np.sqrt(np.sum(weights * odd ** 2))
    even /= even_norm
    d_even /= even_norm
    odd /= odd_norm
    d_odd /= odd_norm
    envelope = (even + odd) / np.sqrt(2.0)
    derivative = (d_even + d_odd) / np.sqrt(2.0)
    return envelope, derivative


def continuum_compression(panels=8, quad=64):
    """Gauss compression on (0, 4). Local derivative, owner lobes, owner Löwdin."""
    omega = float(OMEGA)
    nodes, gauss_w = leggauss(int(quad))
    chunks_x = []
    chunks_w = []
    width = 4.0 / int(panels)
    for panel in range(int(panels)):
        left = panel * width
        chunks_x.append(left + 0.5 * width * (nodes + 1.0))
        chunks_w.append(0.5 * width * gauss_w)
    x = np.concatenate(chunks_x)
    weights = np.concatenate(chunks_w)
    norms = _lobe_norms()
    raw = []
    draw = []
    for region in range(PACKET_COUNT):
        envelope, derivative = _continuum_lobes(x, weights, region, norms)
        raw.append(envelope)
        draw.append(derivative)
    columns = np.sqrt(weights)[:, None] * np.column_stack(raw)
    _orthonormal, gram = lowdin(columns)
    evals, evecs = np.linalg.eigh(0.5 * (gram + gram.conj().T))
    inverse_sqrt = (evecs * (1.0 / np.sqrt(evals))) @ evecs.conj().T
    envelopes = np.column_stack(raw) @ inverse_sqrt
    derivatives = np.column_stack(draw) @ inverse_sqrt
    phase = np.exp(1j * float(CALIBRATION["phase"]))
    spinor = np.zeros((x.size, 2, 6), dtype=np.complex128)
    d_spinor = np.zeros_like(spinor)
    for region in range(PACKET_COUNT):
        carrier = np.exp(1j * CARRIER_K * (x - region * ELL))
        spatial = carrier * envelopes[:, region]
        d_spatial = carrier * (derivatives[:, region] + 1j * CARRIER_K * envelopes[:, region])
        plus_upper = spatial / np.sqrt(2.0)
        plus_lower = 1j * spatial / np.sqrt(2.0)
        spinor[:, 0, 2 * region] = plus_upper
        spinor[:, 1, 2 * region] = plus_lower
        d_spinor[:, 0, 2 * region] = d_spatial / np.sqrt(2.0)
        d_spinor[:, 1, 2 * region] = 1j * d_spatial / np.sqrt(2.0)
        spinor[:, 0, 2 * region + 1] = phase * np.conjugate(plus_upper)
        spinor[:, 1, 2 * region + 1] = phase * np.conjugate(plus_lower)
        d_spinor[:, 0, 2 * region + 1] = phase * np.conjugate(d_spinor[:, 0, 2 * region])
        d_spinor[:, 1, 2 * region + 1] = phase * np.conjugate(d_spinor[:, 1, 2 * region])
    length = float(CALIBRATION["b0"]) * omega ** x
    radial = np.full(x.shape, float(CALIBRATION["b0"]) / float(CALIBRATION["a0"]))
    shift = float(CALIBRATION["beta0"]) * omega ** x
    log_omega = float(np.log(omega))

    def images(kappa):
        upper_s = spinor[:, 0, :]
        lower_s = spinor[:, 1, :]
        derivative_of = {id(upper_s): d_spinor[:, 0, :], id(lower_s): d_spinor[:, 1, :]}

        def anti(values, samples):
            return _anticommutator_closed(values, samples, derivative_of[id(samples)], log_omega)

        upper, lower = _spinor_image(length, radial, shift, kappa, upper_s, lower_s, anti)
        matrix = np.zeros((6, 6), dtype=np.complex128)
        for row in range(6):
            for col in range(6):
                matrix[row, col] = np.sum(
                    weights * (
                        np.conjugate(spinor[:, 0, row]) * upper[:, col]
                        + np.conjugate(spinor[:, 1, row]) * lower[:, col]
                    )
                )
        return matrix

    matrix = images(KAPPA)
    massless = images(0)
    support = []
    for region in range(PACKET_COUNT):
        mass = float(np.sum(weights * np.abs(envelopes[:, region]) ** 2))
        outside = (x < region) | (x > region + 2.0)
        support.append({
            "region": region,
            "mass": mass,
            "mass_outside_naive_support": float(np.sum(weights[outside] * np.abs(envelopes[outside, region]) ** 2)),
        })
    return {
        "quad_panels": int(panels),
        "quad_points_per_panel": int(quad),
        "matrix": matrix,
        "massless_kappa0": massless,
        "gram_offdiag_max": _maxabs(gram - np.eye(3)),
        "lowdin_deviation_max": _maxabs(inverse_sqrt - np.eye(3)),
        "support": support,
        "grading": grading_report(matrix, omega),
        "massless_grading_corner": grading_report(massless, omega)["corner_max"],
        "massless_diagonal_scale_1": grading_report(massless, omega)["diagonal_scale_1"],
    }


def inspect_lowdin(points):
    """Owner lobe samples, owner Löwdin, and where the envelopes actually sit."""
    count = int(points)
    spacing = PERIOD / count
    xi = np.arange(count, dtype=float) * spacing
    norms = _lobe_norms()
    raw = []
    for region in range(PACKET_COUNT):
        even = _sample_lobe(xi, region * ELL, "even", spacing, norms)
        odd = _sample_lobe(xi, region * ELL + ELL, "odd", spacing, norms)
        even = even / np.linalg.norm(even)
        odd = odd / np.linalg.norm(odd)
        raw.append((even + odd) / np.sqrt(2.0))
    raw = np.column_stack(raw)
    orthonormal, gram = lowdin(raw)
    phi0, phi1, preparation = prepare_rank6(xi, spacing)
    problems = preparation_problems(preparation, phi0, phi1, xi)
    envelopes = np.zeros((count, PACKET_COUNT), dtype=float)
    for region in range(PACKET_COUNT):
        local = xi - region * ELL
        carrier = np.sqrt(2.0) * phi0[:, 2 * region] * np.exp(-1j * CARRIER_K * local)
        envelopes[:, region] = np.real(carrier)
    outside = []
    bridge = xi > 4.0
    for region in range(PACKET_COUNT):
        naive = (xi >= region) & (xi <= region + 2.0)
        column = envelopes[:, region]
        outside.append(float(np.sum(column[~naive] ** 2)))
    parity_nodes, parity_weights = leggauss(128)
    parity_u = 0.5 * (parity_nodes + 1.0)
    parity_w = 0.5 * parity_weights
    parity_bump = _bump(parity_u)
    parity_overlap = float(np.sum(parity_w * parity_bump * (parity_u - 0.5) * parity_bump))
    return {
        "points": count,
        "problems": list(problems),
        "gram_offdiag_max": _maxabs(gram - np.diag(np.diag(gram))),
        "lowdin_deviation_max": _maxabs(orthonormal - raw),
        "envelope_imag_after_lowdin": float(preparation["real_envelope_imag_after_lowdin"]),
        "real_envelope_orthonormal_defect": float(preparation["real_envelope_orthonormal_defect"]),
        "column_orthonormal_defect": float(preparation["column_orthonormal_defect"]),
        "phase_on_minus_spinor_column": bool(preparation["phase_on_minus_spinor_column"]),
        "phase_applied_to_odd_lobe": bool(preparation["phase_applied_to_odd_lobe"]),
        "carrier": preparation["carrier"],
        "minus_columns": preparation["minus_columns"],
        "half_density": preparation["half_density"],
        "support_outside_naive": outside,
        "bridge_packet_mass": float(np.sum(envelopes[bridge] ** 2)),
        "reconstructed_envelope_max": _maxabs(envelopes - np.real(orthonormal)),
        "parity_even_odd_overlap": parity_overlap,
        "continuum_lobe_norms": {
            "even": float(preparation["continuum_lobe_norms"]["even"]),
            "odd": float(preparation["continuum_lobe_norms"]["odd"]),
        },
    }


def _images(phi0, phi1, length, radial, shift, kappa, momentum):
    image0, image1 = apply_dirac(phi0, phi1, length, radial, shift, kappa, momentum)
    matrix = phi0.conj().T @ image0 + phi1.conj().T @ image1
    return image0, image1, matrix


def _leakage(phi0, phi1, image0, image1, matrix):
    columns = np.vstack((phi0, phi1))
    images = np.vstack((image0, image1))
    gram = phi0.conj().T @ phi0 + phi1.conj().T @ phi1
    residual = images - columns @ np.linalg.solve(gram, matrix)
    energy = np.linalg.norm(images, axis=0)
    fraction = np.linalg.norm(residual, axis=0) / np.maximum(energy, 1e-30)
    return np.asarray(fraction, dtype=float), gram


def _bridge_fraction(xi, image0, image1):
    bridge = np.asarray(xi) > 4.0
    total = float(np.sum(np.abs(image0) ** 2 + np.abs(image1) ** 2))
    leaked = float(np.sum(np.abs(image0[bridge]) ** 2 + np.abs(image1[bridge]) ** 2))
    return leaked / total


def discrete_compression(fermions, radial=None, phi0=None, phi1=None, radius=None, grid=None, compare_direct=True):
    """Prolonged owner columns and owner Dirac operator on the Galerkin quadrature."""
    if grid is None:
        grid = build_grid(int(fermions))
    elif int(grid.nf) != int(fermions):
        raise ValueError("supplied grid is not this fermion count")
    owned_phi0, owned_phi1, preparation, problems = load_physical_columns(grid.nf)
    if phi0 is None:
        phi0, phi1 = owned_phi0, owned_phi1
    else:
        phi0 = np.asarray(phi0)
        phi1 = np.asarray(phi1)
        if phi0.shape != owned_phi0.shape or phi1.shape != owned_phi1.shape:
            raise ValueError("saved columns are not on this fermion grid")
        problems = preparation_problems(preparation, owned_phi0, owned_phi1, grid.xi_f)
    state = prolong_state(grid, blank_state(grid, phi0, phi1))
    if radial is None:
        radial_fine = state.Q
    else:
        radial = np.asarray(radial, dtype=float)
        if radial.shape != (grid.ng,):
            raise ValueError("saved Q is not on the geometry grid")
        radial_fine = prolong_geometry(grid, radial)
    if radius is None:
        radius_fine = state.r
    else:
        radius = np.asarray(radius, dtype=float)
        if radius.shape != (grid.ng,):
            raise ValueError("saved radius is not on the geometry grid")
        radius_fine = prolong_geometry(grid, radius)
    length = grid.fine.length_density
    shift = grid.fine.shift
    coordinate = profile_coordinate(grid.xi_q, grid.length)
    gauge = {
        "profile_equals_x_on_packet_arc": _maxabs(coordinate[grid.xi_q <= 4.0] - grid.xi_q[grid.xi_q <= 4.0]),
        "length_matches_supplied_gauge": _maxabs(
            length - float(CALIBRATION["b0"]) * float(OMEGA) ** coordinate
        ),
        "shift_matches_supplied_gauge": _maxabs(
            shift - float(CALIBRATION["beta0"]) * float(OMEGA) ** coordinate
        ),
        "initial_Q_span": [float(np.min(state.Q)), float(np.max(state.Q))],
    }
    image0, image1, matrix = _images(
        state.phi0, state.phi1, length, radial_fine, shift, grid.fine.kappa, grid.fine.momentum,
    )
    leak, gram = _leakage(state.phi0, state.phi1, image0, image1, matrix)
    zeros = np.zeros_like(length)
    _k0, _k1, kinetic = _images(
        state.phi0, state.phi1, length, radial_fine, zeros, 0, grid.fine.momentum,
    )
    _s0, _s1, shift_only = _images(
        state.phi0, state.phi1, zeros, radial_fine, shift, 0, grid.fine.momentum,
    )
    mass = matrix - kinetic - shift_only
    if compare_direct:
        metric = CovariantStaticMetric(
            grid.length,
            length * radius_fine,
            radial_fine * radius_fine,
            radius_fine,
            eta=0.5,
        )
        dense = direct_hamiltonian(metric, shift, int(grid.fine.kappa))
        columns = np.vstack((state.phi0, state.phi1))
        direct_gap = _maxabs(dense @ columns - np.vstack((image0, image1)))
    else:
        direct_gap = None
    return {
        "grid": grid,
        "phi0": state.phi0,
        "phi1": state.phi1,
        "radial": radial_fine,
        "matrix": matrix,
        "kinetic": kinetic,
        "mass": mass,
        "shift": shift_only,
        "leakage": leak,
        "gram_defect": _maxabs(gram - np.eye(6)),
        "bridge_image_fraction": _bridge_fraction(grid.xi_q, image0, image1),
        "packet_bridge_mass": float(
            np.sum(np.abs(state.phi0[grid.xi_q > 4.0]) ** 2 + np.abs(state.phi1[grid.xi_q > 4.0]) ** 2)
        ),
        "direct_operator_gap": direct_gap,
        "gauge": gauge,
        "preparation_problems": list(problems),
        "term_sum_gap": _maxabs(kinetic + mass + shift_only - matrix),
    }


def _action_witness(coefficients, q_initial, q_final):
    """Same bulk density at two Q values. The callable is the action owner."""
    lapse = float(CALIBRATION["b0"])
    shift = float(CALIBRATION["beta0"])

    def density(radial):
        # Nonzero r_x makes the -Z L r_x^2 / Q term visible. Same callable either Q.
        return float(first_order_density(
            lapse, 0.0, float(radial), 0.0, 0.0,
            1.0, 0.0, 0.2, 0.0, 0.0, 0.0,
            shift, 0.0,
            coefficients["A"], coefficients["C_W"], coefficients["C_F"], coefficients["flux"],
        ))

    initial = density(q_initial)
    final = density(q_final)
    return {
        "owner": "nsc_spherical_feedback_action.first_order_density",
        "functional_edited": False,
        "sample_Q_initial": float(q_initial),
        "sample_Q_final": float(q_final),
        "density_initial": initial,
        "density_final": final,
        "density_changed": bool(abs(final - initial) > 1e-12),
    }


def _feedback_compression(initial):
    """Saved T=0.05 geometry and columns. Does not call the episode stepper."""
    if not EPISODE_NPZ.is_file() or not EPISODE_JSON.is_file():
        raise FileNotFoundError("saved feedback episode payload is missing")
    episode = json.loads(EPISODE_JSON.read_text())
    with np.load(EPISODE_NPZ, allow_pickle=False) as payload:
        radial = np.array(payload[EPISODE_RUN + "_final_Q"], dtype=float, copy=True)
        radius = np.array(payload[EPISODE_RUN + "_final_r"], dtype=float, copy=True)
        phi0 = np.array(payload[EPISODE_RUN + "_final_phi0"], copy=True)
        phi1 = np.array(payload[EPISODE_RUN + "_final_phi1"], copy=True)
        series_time = float(payload[EPISODE_RUN + "_time"][-1])
        frame_time = float(payload[EPISODE_RUN + "_frame_time"][-1])
        frame_q = np.array(payload[EPISODE_RUN + "_frame_Q"][-1], dtype=float, copy=True)
    if abs(series_time - 0.05) > 1e-12 or abs(frame_time - 0.05) > 1e-12:
        raise ValueError("saved episode does not end at T=0.05")
    if np.min(radius) <= 0 or np.min(radial) <= 0:
        raise ValueError("saved final chart is not positive")
    same_packets = discrete_compression(
        512, radial=radial, grid=initial["grid"], compare_direct=False,
    )
    final_packets = discrete_compression(
        512, radial=radial, phi0=phi0, phi1=phi1, radius=radius,
        grid=initial["grid"], compare_direct=False,
    )
    if _maxabs(same_packets["radial"] - frame_q) > 1e-12:
        raise ValueError("prolonged final Q does not match the saved quadrature frame")
    q0 = float(CALIBRATION["b0"]) / float(CALIBRATION["a0"])
    ratio = q0 / same_packets["radial"]
    coefficients = same_packets["grid"].fine.coefficients
    q_final = float(radial[len(radial) // 2])
    return {
        "payload": "results/development/nsc-spherical-feedback-episode-v1.npz",
        "run": EPISODE_RUN,
        "attained_T": series_time,
        "evolution_rerun": False,
        "saved_episode_renewal": bool(episode["renewal"]),
        "Q_min": float(np.min(radial)),
        "Q_max": float(np.max(radial)),
        "q0_over_Q_min": float(np.min(ratio)),
        "q0_over_Q_max": float(np.max(ratio)),
        "frame_match": _maxabs(same_packets["radial"] - frame_q),
        "initial_packets_on_final_Q": _public_compression(same_packets, float(OMEGA)),
        "final_packets_on_final_Q": _public_compression(final_packets, float(OMEGA)),
        "action": _action_witness(coefficients, q0, q_final),
        "input_sha256": {
            "episode_json": _sha256(EPISODE_JSON),
            "episode_npz": _sha256(EPISODE_NPZ),
        },
    }


def _public_grading(report):
    hidden = {"H", "B01"}
    return {key: value for key, value in report.items() if key not in hidden}


def _public_compression(result, omega):
    full = grading_report(result["matrix"], omega)
    kinetic = grading_report(result["kinetic"], omega)
    mass = grading_report(result["mass"], omega)
    shift = grading_report(result["shift"], omega)
    leak = np.asarray(result["leakage"], dtype=float)
    return {
        "grading": _public_grading(full),
        "H": _complex_record(full["H"]),
        "B01": _complex_record(full["B01"]),
        "kinetic_diagonal_scale_1": kinetic["diagonal_scale_1"],
        "mass_diagonal_scale_1": mass["diagonal_scale_1"],
        "shift_diagonal_scale_1": shift["diagonal_scale_1"],
        "kinetic_link_scale_1": kinetic["link_scale_1"],
        "mass_link_scale_1": mass["link_scale_1"],
        "shift_link_scale_1": shift["link_scale_1"],
        "term_sum_gap": result["term_sum_gap"],
        "projection_leakage": leak.tolist(),
        "projection_leakage_max": float(np.max(leak)),
        "captured_fraction_min": float(np.min(np.sqrt(np.maximum(0.0, 1.0 - leak ** 2)))),
        "gram_defect": result["gram_defect"],
        "bridge_image_fraction": result["bridge_image_fraction"],
        "packet_bridge_mass": result["packet_bridge_mass"],
        "direct_operator_gap": result["direct_operator_gap"],
        "corner_frobenius": full["corner_frobenius"],
    }


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return _complex_record(value)
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def measure():
    """nf=256 and nf=512 compressions, continuum check, and the saved final frame."""
    started = time.process_time()
    omega = float(OMEGA)
    translation = translation_algebra()
    lowdin_report = inspect_lowdin(512)
    continuum = continuum_compression()
    continuum_fine = continuum_compression(panels=8, quad=96)
    discrete = {}
    for fermions in (256, 512):
        discrete[fermions] = discrete_compression(fermions)
    primary = discrete[512]
    graded = grading_report(primary["matrix"], omega)
    algebra = resolvent_report(graded["H"], graded["B01"], omega, primary["matrix"])
    feedback = _feedback_compression(primary)
    historical = exact_controls()
    cpu = float(time.process_time() - started)
    if cpu > CPU_LIMIT_S:
        raise RuntimeError(f"numerical CPU {cpu} exceeds {CPU_LIMIT_S}")
    leak = np.asarray(primary["leakage"], dtype=float)
    discrete_gap = _maxabs(discrete[256]["matrix"] - primary["matrix"])
    continuum_gap = _maxabs(continuum["matrix"] - primary["matrix"])
    quadrature_gap = _maxabs(continuum_fine["matrix"] - continuum["matrix"])
    closed = bool(float(np.max(leak)) <= LEAK_CLOSURE_MAX)
    payload = {
        "schema": SCHEMA,
        "status": STATUS,
        "cpu_seconds": cpu,
        "cpu_limit_seconds": CPU_LIMIT_S,
        "operator": "sigma2 {L/Q, P}/2 + sigma1 kappa L - {beta, P}/2",
        "operator_owner": "nsc_spherical_coupling.apply_dirac",
        "direct_operator_owner": "nsc_conformal_adm_source.direct_hamiltonian",
        "packet_owner": "nsc_spherical_coupling.prepare_rank6",
        "galerkin_owner": "nsc_spherical_galerkin_coupling.build_grid",
        "action_owner": "nsc_spherical_feedback_action.first_order_density",
        "theorem_owner": "docs/nsc-nested-qualities.md",
        "omega": omega,
        "kappa": int(KAPPA),
        "q0": float(CALIBRATION["b0"]) / float(CALIBRATION["a0"]),
        "calibration": {
            "a0": float(CALIBRATION["a0"]),
            "beta0": float(CALIBRATION["beta0"]),
            "b0": float(CALIBRATION["b0"]),
            "phase": float(CALIBRATION["phase"]),
            "k": float(CARRIER_K),
            "ell": float(ELL),
            "period": float(PERIOD),
        },
        "relation": {
            "J_nn": "Omega^n H_derived",
            "J_n_n1": "Omega^n B01_derived",
            "block_coefficient": "B01_derived",
            "congruence_normalized_B": "Omega^{-1/2} B01_derived",
            "published_recurrence_uses": "B01_derived with the factor 1/Omega",
            "normalized_B_is_not_the_block_coefficient": True,
        },
        "translation": translation,
        "lowdin": lowdin_report,
        "gauge_nf512": discrete[512]["gauge"],
        "continuum": {
            "quadrature_gap_64_vs_96": quadrature_gap,
            "versus_nf512_max": continuum_gap,
            "gram_offdiag_max": continuum["gram_offdiag_max"],
            "lowdin_deviation_max": continuum["lowdin_deviation_max"],
            "support": continuum["support"],
            "grading": _public_grading(continuum["grading"]),
            "massless_kappa0_diagonal_scale_1": continuum["massless_diagonal_scale_1"],
            "massless_kappa0_corner": continuum["massless_grading_corner"],
            "massless_H": _complex_record(grading_report(continuum["massless_kappa0"], omega)["H"]),
        },
        "nf256": _public_compression(discrete[256], omega),
        "nf512": _public_compression(primary, omega),
        "derived": {
            "H": _complex_record(graded["H"]),
            "B01": _complex_record(graded["B01"]),
            "normalized_B": _complex_record(algebra["normalized_B"]),
            "from": "nf512 prolonged prepare_rank6 columns, initial constant Q",
        },
        "algebra": {
            key: (_complex_record(value) if isinstance(value, np.ndarray) else value)
            for key, value in algebra.items()
            if key != "normalized_B"
        },
        "historical_rational_example": {
            "fitted": False,
            "role": "optional example in nsc_nested_qualities.exact_controls",
            "example_checks_pass": bool(all(historical["checks"].values())),
            "example_omega": historical["omega"],
            "H_gap": _maxabs(graded["H"] - _EXAMPLE_H),
            "B_gap": _maxabs(graded["B01"] - _EXAMPLE_B),
        },
        "comparison": {
            "discrete_max_abs_nf256_vs_nf512": discrete_gap,
            "discrete_frobenius_nf256_vs_nf512": _frobenius(discrete[256]["matrix"] - primary["matrix"]),
            "continuum_vs_nf512_max": continuum_gap,
            "projection_leakage_max_nf512": float(np.max(leak)),
            "projection_leakage_max_nf256": float(np.max(discrete[256]["leakage"])),
            "corner_frobenius_nf512": graded["corner_frobenius"],
            "corner_frobenius_continuum": continuum["grading"]["corner_frobenius"],
            "link_frobenius_nf512": graded["link_frobenius"],
            "leakage_and_discrete_error_are_distinct": True,
        },
        "feedback_T_0_05": feedback,
        "closed_six_mode": closed,
        "closure_rule": "closed_six_mode requires projection leakage at or below 0.1 on every column",
        "geometry_grading_is_supplied_input": True,
        "universal_physical_fractality": False,
        "incoming_gate_changed": False,
        "infinite_nest_required": False,
        "lambdacdm_required": False,
        "renewal_established": False,
        "new_evolution": False,
        "assumptions": [
            "One angular block of the owner ordering at the calibration kappa.",
            "L and beta are the supplied samples b0*Omega**s and beta0*Omega**s.",
            "s equals x on the packet arc [0, 4]; outside it, s is the owner bridge.",
            "Initial Q is the constant b0/a0. Initial radius is 1 and cancels in L and Q.",
            "Columns are prepare_rank6, prolonged by the Galerkin fermion map.",
            "J is V dagger H V in that orthonormal frame. Occupations do not enter J.",
            "The depth-3 shift covariance checked here is the scaling of these three regions.",
            "The saved T=0.05 frame is read from the episode payload and is not recomputed.",
        ],
        "missing_connection": (
            "The theorem turns these derived blocks into a finite window and a Schur "
            "response. The Dirac image is mostly outside the six columns, so that "
            "window is not the evolution. At the saved T=0.05 frame the same action "
            "functional and the same supplied L and beta gauge remain, while "
            "nonconstant Q removes the Omega^n scaling of the kinetic block. "
            "Renewal, an invariant six-mode structure, and a local response that "
            "stays inside the window are not obtained. The geometric grading is the "
            "supplied initial gauge."
        ),
    }
    return _jsonable(payload)


def write_record(payload, path=RECORD_PATH):
    """Store the measurement and the hashes of the files that define it."""
    record = dict(payload)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    record["source_hashes"] = {
        "src/recursive_horizons/nsc_geometry_graded_window.py": _sha256(Path(__file__)),
        "scripts/derive_nsc_geometry_graded_window.py": _sha256(
            _LAB / "scripts" / "derive_nsc_geometry_graded_window.py"
        ),
        "tests/test_nsc_geometry_graded_window.py": _sha256(
            _LAB / "tests" / "test_nsc_geometry_graded_window.py"
        ),
        "docs/nsc-geometry-graded-window.md": _sha256(
            _LAB / "docs" / "nsc-geometry-graded-window.md"
        ),
    }
    path.write_text(json.dumps(record, indent=2) + "\n")
    return path


def _compare(saved, fresh, path, mismatches):
    if isinstance(fresh, dict):
        if not isinstance(saved, dict) or set(saved) != set(fresh):
            mismatches.append(path or "root")
            return
        for key in fresh:
            _compare(saved[key], fresh[key], f"{path}.{key}" if path else key, mismatches)
        return
    if isinstance(fresh, list):
        if not isinstance(saved, list) or len(saved) != len(fresh):
            mismatches.append(path)
            return
        for index, (left, right) in enumerate(zip(saved, fresh)):
            _compare(left, right, f"{path}[{index}]", mismatches)
        return
    if fresh is None or saved is None:
        if saved is not None or fresh is not None:
            mismatches.append(path)
        return
    if isinstance(fresh, float):
        scale = max(1.0, abs(fresh), abs(float(saved)))
        if abs(float(saved) - fresh) > 1e-9 * scale:
            mismatches.append(path)
        return
    if isinstance(fresh, (bool, int, str)):
        if saved != fresh:
            mismatches.append(path)
        return
    mismatches.append(path)


def verify_saved(path=RECORD_PATH):
    """Recompute the measurement and compare it with the stored record."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    saved = json.loads(path.read_text())
    fresh = measure()
    if saved.get("schema") != SCHEMA or fresh["schema"] != SCHEMA:
        raise AssertionError("schema mismatch")
    if float(saved["cpu_seconds"]) > CPU_LIMIT_S or float(fresh["cpu_seconds"]) > CPU_LIMIT_S:
        raise AssertionError("CPU limit exceeded")
    mismatches = []
    for key, value in fresh.items():
        if key == "cpu_seconds":
            continue
        _compare(saved.get(key), value, key, mismatches)
    hashes = saved.get("source_hashes") or {}
    current_hashes = {
        "src/recursive_horizons/nsc_geometry_graded_window.py": _sha256(Path(__file__)),
        "scripts/derive_nsc_geometry_graded_window.py": _sha256(
            _LAB / "scripts" / "derive_nsc_geometry_graded_window.py"
        ),
        "tests/test_nsc_geometry_graded_window.py": _sha256(
            _LAB / "tests" / "test_nsc_geometry_graded_window.py"
        ),
        "docs/nsc-geometry-graded-window.md": _sha256(
            _LAB / "docs" / "nsc-geometry-graded-window.md"
        ),
    }
    for name, digest in current_hashes.items():
        if hashes.get(name) != digest:
            mismatches.append("source_hashes." + name)
    if mismatches:
        raise AssertionError("record drift: " + ", ".join(mismatches[:12]))
    return {
        "status": fresh["status"],
        "cpu_seconds": fresh["cpu_seconds"],
        "projection_leakage_max": fresh["comparison"]["projection_leakage_max_nf512"],
        "discrete_max_abs": fresh["comparison"]["discrete_max_abs_nf256_vs_nf512"],
        "closed_six_mode": fresh["closed_six_mode"],
    }
