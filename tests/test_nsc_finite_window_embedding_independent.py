"""Independent integrals for the degree-4 finite-window compression.

The matrix residual is a floating comparison of published decimal lobe
coefficients with the rational window. It is not an exact-root theorem.
The exterior complement is the intended dynamics: these tests do not ask
the six modes to be an invariant subspace.
"""
from __future__ import annotations

import os

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

from functools import lru_cache

import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss

from recursive_horizons.nsc_finite_window_embedding import (
    A0,
    B0,
    BETA0,
    COLUMN_PHASE,
    DEGREE,
    K_CARRIER,
    LAMBDA,
    MINUS_LEFT,
    MINUS_RIGHT,
    MINUS_RIGHT_PHASE,
    OMEGA,
    PLUS_LEFT,
    PLUS_RIGHT,
    PLUS_RIGHT_PHASE,
    frozen_window,
)
from recursive_horizons.nsc_nested_qualities import finite_window
from recursive_horizons.nsc_spherical_coupling import (
    KAPPA,
    PERIOD,
    apply_dirac,
    antiperiodic_momentum,
)

PI = np.pi
LOCAL = sp.Matrix([[1, sp.I / 5], [-sp.I / 5, 2]])
LINK = sp.Matrix(
    [[sp.Rational(1, 4), sp.I / 7], [sp.Rational(1, 9), sp.Rational(1, 6)]]
)


def _icos(alpha, beta):
    alpha = complex(alpha)
    beta = float(beta)
    if abs(beta) <= 1e-14 and abs(alpha) <= 1e-14:
        return 1.0 + 0.0j
    if abs(beta) <= 1e-14:
        return (np.exp(alpha) - 1.0) / alpha
    den = alpha * alpha + beta * beta
    return (np.exp(alpha) * (alpha * np.cos(beta) + beta * np.sin(beta)) - alpha) / den


def _isin(alpha, beta):
    alpha = complex(alpha)
    beta = float(beta)
    if abs(beta) <= 1e-14:
        return 0.0j
    den = alpha * alpha + beta * beta
    return (np.exp(alpha) * (alpha * np.sin(beta) - beta * np.cos(beta)) + beta) / den


def _uu(alpha, m, n):
    return _icos(alpha, PI * (m - n)) - _icos(alpha, PI * (m + n))


def _u_up(alpha, m, n):
    return (PI * n) * (_isin(alpha, PI * (m + n)) + _isin(alpha, PI * (m - n)))


def _up_u(alpha, m, n):
    return (PI * m) * (_isin(alpha, PI * (m + n)) - _isin(alpha, PI * (m - n)))


def _up_up(alpha, m, n):
    return (PI * m) * (PI * n) * (
        _icos(alpha, PI * (m - n)) + _icos(alpha, PI * (m + n))
    )


def _phase(lobe):
    return np.exp(1j * lobe["k"] * lobe["delta"])


def _pair(left, right, alpha, kind):
    acc = 0.0j
    for m in range(1, DEGREE + 1):
        for n in range(1, DEGREE + 1):
            kernel = _uu(alpha, m, n) if kind == "value" else _u_up(alpha, m, n)
            acc += np.conjugate(left["coeff"][m - 1]) * right["coeff"][n - 1] * kernel
    return np.conjugate(_phase(left)) * _phase(right) * acc


def _shared(mode_a, mode_b, weighted, deriv_on_b):
    total = 0.0j
    for lobe_a in mode_a:
        for lobe_b in mode_b:
            if lobe_a["origin"] != lobe_b["origin"]:
                continue
            alpha = 1j * (lobe_b["k"] - lobe_a["k"])
            if weighted:
                alpha = alpha + LAMBDA
            if deriv_on_b:
                bare = _pair(lobe_a, lobe_b, alpha, "deriv")
                bare += (1j * lobe_b["k"] + 0.5 * LAMBDA) * _pair(
                    lobe_a, lobe_b, alpha, "value"
                )
            else:
                bare = _pair(lobe_a, lobe_b, alpha, "value")
            if weighted:
                bare *= OMEGA ** lobe_a["origin"]
            total += bare
    return total


def _element(observer, source, spin_o, spin_s):
    if spin_o == "plus" and spin_s == "plus":
        return 1j * (BETA0 - A0) * _shared(observer, source, True, True)
    if spin_o == "plus" and spin_s == "minus":
        return -1j * B0 * _shared(observer, source, True, False)
    if spin_o == "minus" and spin_s == "plus":
        return 1j * B0 * _shared(observer, source, True, False)
    return 1j * (BETA0 + A0) * _shared(observer, source, True, True)


def _lobes(region, left, right, k):
    return (
        {"origin": float(region), "delta": 0.0, "coeff": np.asarray(left, dtype=np.complex128), "k": k, "region": region},
        {"origin": float(region + 1), "delta": 1.0, "coeff": np.asarray(right, dtype=np.complex128), "k": k, "region": region},
    )


def _stored_modes():
    labels = [(region, spin) for region in range(3) for spin in ("plus", "minus")]
    modes = []
    for region, spin in labels:
        if spin == "plus":
            left = np.asarray(PLUS_LEFT, dtype=np.complex128)
            right = np.exp(1j * PLUS_RIGHT_PHASE) * np.asarray(PLUS_RIGHT, dtype=np.complex128)
            k = K_CARRIER
        else:
            left = np.exp(1j * COLUMN_PHASE) * np.asarray(MINUS_LEFT, dtype=np.complex128)
            right = np.exp(1j * (COLUMN_PHASE + MINUS_RIGHT_PHASE)) * np.asarray(
                MINUS_RIGHT, dtype=np.complex128
            )
            k = -K_CARRIER
        modes.append(_lobes(region, left, right, k))
    return labels, modes


def _image_norm_squared(mode, spin):
    total = 0.0
    alpha = 2.0 * LAMBDA
    for lobe in mode:
        coeff = lobe["coeff"]
        gamma = 1j * lobe["k"] + 0.5 * LAMBDA
        deriv = 0.0j
        value = 0.0j
        for m in range(1, DEGREE + 1):
            for n in range(1, DEGREE + 1):
                weight = coeff[m - 1] * np.conjugate(coeff[n - 1])
                deriv += weight * (
                    _up_up(alpha, m, n)
                    + np.conjugate(gamma) * _up_u(alpha, m, n)
                    + gamma * _u_up(alpha, m, n)
                    + (gamma * np.conjugate(gamma)) * _uu(alpha, m, n)
                )
                value += weight * _uu(alpha, m, n)
        scale = (OMEGA ** lobe["origin"]) ** 2
        if spin == "plus":
            total += (A0 - BETA0) ** 2 * scale * deriv.real + (B0 ** 2) * scale * value.real
        else:
            total += (B0 ** 2) * scale * value.real + (A0 + BETA0) ** 2 * scale * deriv.real
    return total


def _image_pair(mode_a, spin_a, mode_b, spin_b):
    """L2 inner product of two continuum images, same spin or opposite."""
    total = 0.0j
    for lobe_a in mode_a:
        for lobe_b in mode_b:
            if lobe_a["origin"] != lobe_b["origin"]:
                continue
            alpha = 2.0 * LAMBDA + 1j * (lobe_b["k"] - lobe_a["k"])
            pre = (
                np.conjugate(_phase(lobe_a))
                * _phase(lobe_b)
                * (OMEGA ** lobe_a["origin"]) ** 2
            )
            gamma_a = 1j * lobe_a["k"] + 0.5 * LAMBDA
            gamma_b = 1j * lobe_b["k"] + 0.5 * LAMBDA
            deriv = value = da_fb = fa_db = 0.0j
            for m in range(1, DEGREE + 1):
                for n in range(1, DEGREE + 1):
                    amp = np.conjugate(lobe_a["coeff"][m - 1]) * lobe_b["coeff"][n - 1]
                    deriv += amp * (
                        _up_up(alpha, m, n)
                        + gamma_b * _up_u(alpha, m, n)
                        + np.conjugate(gamma_a) * _u_up(alpha, m, n)
                        + np.conjugate(gamma_a) * gamma_b * _uu(alpha, m, n)
                    )
                    value += amp * _uu(alpha, m, n)
                    da_fb += amp * (_up_u(alpha, m, n) + np.conjugate(gamma_a) * _uu(alpha, m, n))
                    fa_db += amp * (_u_up(alpha, m, n) + gamma_b * _uu(alpha, m, n))
            deriv *= pre
            value *= pre
            da_fb *= pre
            fa_db *= pre
            if spin_a == spin_b == "plus":
                total += (A0 - BETA0) ** 2 * deriv + (B0 ** 2) * value
            elif spin_a == spin_b == "minus":
                total += (B0 ** 2) * value + (A0 + BETA0) ** 2 * deriv
            elif spin_a == "plus" and spin_b == "minus":
                total += np.conjugate(1j * (BETA0 - A0)) * (-1j * B0) * da_fb
                total += np.conjugate(1j * B0) * (1j * (BETA0 + A0)) * fa_db
            else:
                total += np.conjugate(-1j * B0) * (1j * (BETA0 - A0)) * fa_db
                total += np.conjugate(1j * (BETA0 + A0)) * (1j * B0) * da_fb
    return total


def _rational_window():
    window = finite_window(LOCAL, LINK, sp.Rational(3, 2), 0, 3)
    placed = sp.zeros(6)
    omega = sp.Rational(3, 2)
    for region in range(3):
        block = slice(2 * region, 2 * region + 2)
        placed[block, block] = (omega ** region) * LOCAL
        if region < 2:
            nxt = slice(2 * region + 2, 2 * region + 4)
            placed[block, nxt] = (omega ** region) * LINK
            placed[nxt, block] = (omega ** region) * LINK.H
    return window, placed


def _as_numpy(matrix):
    return np.array(
        [[complex(matrix[i, j]) for j in range(6)] for i in range(6)],
        dtype=np.complex128,
    )


@lru_cache(maxsize=1)
def _compression():
    labels, modes = _stored_modes()
    gram = np.zeros((6, 6), dtype=np.complex128)
    matrix = np.zeros((6, 6), dtype=np.complex128)
    for i, (_, spin_i) in enumerate(labels):
        for j, (_, spin_j) in enumerate(labels):
            if spin_i == spin_j:
                gram[i, j] = _shared(modes[i], modes[j], False, False)
            matrix[i, j] = _element(modes[i], modes[j], spin_i, spin_j)
    return labels, modes, gram, matrix


def _overlap(left_coeff, right_coeff):
    """∫ ρ F G and ∫ ρ F (G' + λ G / 2) on the unit overlap, real lobes."""
    right = np.asarray(right_coeff, dtype=float)
    left = np.asarray(left_coeff, dtype=float)
    mass = slope = 0.0j
    for m in range(1, DEGREE + 1):
        for n in range(1, DEGREE + 1):
            mass += right[m - 1] * left[n - 1] * _uu(LAMBDA, m, n)
            slope += right[m - 1] * left[n - 1] * _u_up(LAMBDA, m, n)
    mass *= OMEGA
    slope *= OMEGA
    return mass, slope + 0.5 * LAMBDA * mass


def _weighted_sine_gram():
    gram = np.zeros((DEGREE, DEGREE))
    for m in range(1, DEGREE + 1):
        for n in range(1, DEGREE + 1):
            gram[m - 1, n - 1] = (OMEGA * _uu(LAMBDA, m, n)).real
    return gram


def _endpoint_derivative(coeff, at_right_end):
    modes = np.arange(1, DEGREE + 1)
    signs = ((-1.0) ** modes) if at_right_end else np.ones(DEGREE)
    return np.sum(np.asarray(coeff, dtype=np.complex128) * np.sqrt(2.0) * PI * modes * signs)


def _lobe_value_change(coeff, k):
    """∫ (s' + i k s) e^{i k t} dt. Zero when the lobe value vanishes at both ends."""
    alpha = 1j * k
    acc = 0.0j
    series = np.asarray(coeff, dtype=np.complex128)
    for n in range(1, DEGREE + 1):
        acc += series[n - 1] * np.sqrt(2.0) * (
            PI * n * _icos(alpha, PI * n) + 1j * k * _isin(alpha, PI * n)
        )
    return acc


def _discrete(points):
    points = int(points)
    dx = PERIOD / points
    xi = np.arange(points, dtype=float) * dx
    sqrt_dx = np.sqrt(dx)
    phi0 = np.zeros((points, 6), dtype=np.complex128)
    phi1 = np.zeros((points, 6), dtype=np.complex128)
    indices = np.arange(1, DEGREE + 1)
    for region in range(3):
        plus = np.zeros(points, dtype=np.complex128)
        for lobe, phase in ((0.0, 0.0), (1.0, PLUS_RIGHT_PHASE)):
            local = xi - region - lobe
            mask = (local > 0.0) & (local < 1.0)
            sine = np.sin(PI * local[mask, None] * indices[None, :]) * np.sqrt(2.0)
            coeff = PLUS_LEFT if lobe == 0.0 else PLUS_RIGHT
            plus[mask] = np.exp(1j * phase) * (sine @ np.asarray(coeff, dtype=float))
        plus *= np.exp(1j * K_CARRIER * (xi - region))
        minus = np.zeros(points, dtype=np.complex128)
        for lobe, phase in ((0.0, 0.0), (1.0, MINUS_RIGHT_PHASE)):
            local = xi - region - lobe
            mask = (local > 0.0) & (local < 1.0)
            sine = np.sin(PI * local[mask, None] * indices[None, :]) * np.sqrt(2.0)
            coeff = MINUS_LEFT if lobe == 0.0 else MINUS_RIGHT
            minus[mask] = np.exp(1j * phase) * (sine @ np.asarray(coeff, dtype=float))
        minus *= np.exp(1j * COLUMN_PHASE) * np.exp(-1j * K_CARRIER * (xi - region))
        phi0[:, 2 * region] = sqrt_dx * plus / np.sqrt(2.0)
        phi1[:, 2 * region] = sqrt_dx * 1j * plus / np.sqrt(2.0)
        phi0[:, 2 * region + 1] = sqrt_dx * minus / np.sqrt(2.0)
        phi1[:, 2 * region + 1] = sqrt_dx * (-1j) * minus / np.sqrt(2.0)
    indices_p = np.arange(points) - points // 2
    momenta = 2.0 * PI * (indices_p + 0.5) / PERIOD
    grid = (np.arange(points) - points // 2) * dx
    fourier = np.exp(1j * grid[:, None] * momenta[None, :]) / np.sqrt(points)
    momentum = (fourier * momenta[None, :]) @ fourier.conj().T
    momentum = (momentum + momentum.conj().T) / 2.0
    scale = OMEGA ** xi
    length = B0 * scale
    shift = BETA0 * scale
    radial = np.full(points, B0 / A0)

    def anticommutator(values, component):
        return 0.5 * (
            values[:, None] * (momentum @ component)
            + momentum @ (values[:, None] * component)
        )

    image0 = (
        -1j * anticommutator(length / radial, phi1)
        + (KAPPA * length)[:, None] * phi1
        - anticommutator(shift, phi0)
    )
    image1 = (
        1j * anticommutator(length / radial, phi0)
        + (KAPPA * length)[:, None] * phi0
        - anticommutator(shift, phi1)
    )
    produced, _metric = antiperiodic_momentum(points, PERIOD)
    reference0, reference1 = apply_dirac(phi0, phi1, length, radial, shift, KAPPA, produced)
    gram = phi0.conj().T @ phi0 + phi1.conj().T @ phi1
    matrix = phi0.conj().T @ image0 + phi1.conj().T @ image1
    bridge = (xi > 4.0) & (xi < PERIOD)
    bridge_max = float(np.max(np.abs(np.concatenate((image0[bridge], image1[bridge])))))
    return {
        "matrix": matrix,
        "gram_defect": float(np.max(np.abs(gram - np.eye(6)))),
        "operator_disagreement": max(
            float(np.linalg.norm(image0 - reference0)),
            float(np.linalg.norm(image1 - reference1)),
            float(np.linalg.norm(momentum - produced)),
        ),
        "bridge_max": bridge_max,
        "outside_mass": float(np.sum(np.abs(phi0[bridge]) ** 2 + np.abs(phi1[bridge]) ** 2)),
    }


def test_antiderivatives_match_an_independent_gauss_panel():
    nodes, weights = leggauss(80)
    t = 0.5 * (nodes + 1.0)
    wt = 0.5 * weights
    alpha = 0.3 + 0.7j
    worst = 0.0
    for m in range(1, DEGREE + 1):
        for n in range(1, DEGREE + 1):
            um = np.sqrt(2.0) * np.sin(PI * m * t)
            un = np.sqrt(2.0) * np.sin(PI * n * t)
            ump = np.sqrt(2.0) * PI * m * np.cos(PI * m * t)
            unp = np.sqrt(2.0) * PI * n * np.cos(PI * n * t)
            factor = np.exp(alpha * t)
            closed = (
                _uu(alpha, m, n),
                _u_up(alpha, m, n),
                _up_u(alpha, m, n),
                _up_up(alpha, m, n),
            )
            sampled = (
                np.sum(wt * factor * um * un),
                np.sum(wt * factor * um * unp),
                np.sum(wt * factor * ump * un),
                np.sum(wt * factor * ump * unp),
            )
            worst = max(worst, max(abs(a - b) for a, b in zip(closed, sampled)))
    assert worst < 1e-11
    assert abs(_uu(0.0, 2, 2) - 1.0) == 0.0
    assert _uu(0.0, 1, 3) == 0.0


def test_closed_form_compression_matches_the_rational_window_numerically():
    """Direct integrals, not the module quadrature. A small residual is not an identity."""
    _labels, _modes, gram, matrix = _compression()
    window, placed = _rational_window()
    assert window == placed
    target = _as_numpy(window)
    built, symbolic = frozen_window()
    assert symbolic == window
    assert np.allclose(built, target)
    difference = matrix - target
    assert 0.0 < np.linalg.norm(difference) < 1e-12
    assert np.max(np.abs(gram - np.eye(6))) < 1e-12
    assert np.linalg.norm(matrix - matrix.conj().T) < 1e-12
    for region in range(3):
        block = slice(2 * region, 2 * region + 2)
        assert np.linalg.norm(difference[block, block]) < 1e-12
    for region in range(2):
        row = slice(2 * region, 2 * region + 2)
        col = slice(2 * region + 2, 2 * region + 4)
        assert np.linalg.norm(difference[row, col]) < 1e-12
    assert matrix[:2, 4:].tolist() == [[0j, 0j], [0j, 0j]]
    assert target[:2, 4:].tolist() == [[0j, 0j], [0j, 0j]]
    tweaked = np.array(PLUS_LEFT, dtype=float)
    tweaked[0] += 1e-6
    left = tweaked
    right = np.exp(1j * PLUS_RIGHT_PHASE) * np.asarray(PLUS_RIGHT, dtype=float)
    moved = _element(_lobes(0, left, right, K_CARRIER), _lobes(1, left, right, K_CARRIER), "plus", "plus")
    assert abs(moved - matrix[0, 2]) > 1e-6


def test_modes_keep_regional_support_and_unit_gram():
    labels, modes, gram, matrix = _compression()
    assert DEGREE == 4
    for index, (region, spin) in enumerate(labels):
        origins = [lobe["origin"] for lobe in modes[index]]
        assert origins == [region, region + 1]
        assert all(origin <= 3.0 for origin in origins)
        mass = sum(float(np.vdot(lobe["coeff"], lobe["coeff"]).real) for lobe in modes[index])
        assert abs(mass - 1.0) < 1e-12
        assert abs(gram[index, index] - 1.0) < 1e-12
    assert abs(np.dot(PLUS_LEFT, PLUS_RIGHT)) < 1e-12
    assert abs(np.dot(MINUS_LEFT, MINUS_RIGHT)) < 1e-12
    # Region 0 and region 2 meet only at an endpoint, so the corner integral is empty.
    assert matrix[0, 4] == 0j and matrix[1, 5] == 0j
    assert MINUS_RIGHT_PHASE == 0.0
    assert abs(np.sin(PLUS_RIGHT_PHASE)) > 0.5


def test_joins_are_continuous_while_the_first_derivative_jumps():
    """Value jump is zero, so the first weak derivative has no boundary delta."""
    assert abs(_endpoint_derivative(PLUS_LEFT, at_right_end=False)) > 1.0
    plus_left = _endpoint_derivative(PLUS_LEFT, at_right_end=True)
    plus_right = _endpoint_derivative(PLUS_RIGHT, at_right_end=False)
    seam = np.exp(1j * K_CARRIER) * (
        np.exp(1j * PLUS_RIGHT_PHASE) * plus_right - plus_left
    )
    assert abs(seam) > 1.0
    assert abs(_lobe_value_change(PLUS_LEFT, K_CARRIER)) < 1e-12
    assert abs(_lobe_value_change(PLUS_RIGHT, K_CARRIER)) < 1e-12
    assert abs(_lobe_value_change(MINUS_LEFT, -K_CARRIER)) < 1e-12
    # sin(π m t) is zero at both integers, for every coefficient. The joined
    # profile is C^0. A jump of the derivative is a delta in the second
    # distributional derivative only, and the Dirac operator is first order.
    for coeff in (PLUS_LEFT, PLUS_RIGHT, MINUS_LEFT, MINUS_RIGHT):
        series = np.asarray(coeff, dtype=float)
        ends = np.sin(PI * np.arange(1, DEGREE + 1))
        assert np.max(np.abs(ends)) < 1e-15
        assert abs(np.dot(series, ends)) < 1e-15


def test_real_lobes_cannot_reach_the_plus_link_and_conjugate_packets_cannot():
    gap = A0 - BETA0
    required = (0.25 / gap) * float(np.cos(K_CARRIER)) / K_CARRIER
    popoviciu = (float(OMEGA) ** 2 - float(OMEGA)) / 4.0
    assert required > popoviciu + 0.01
    eigenvalues = np.linalg.eigvalsh(_weighted_sine_gram())
    degree4_max = (eigenvalues[-1] - eigenvalues[0]) / 4.0
    assert 0.05 < degree4_max < popoviciu
    assert degree4_max < required - 0.05
    rng = np.random.default_rng(7)
    for _ in range(8):
        left = rng.normal(size=DEGREE)
        right = rng.normal(size=DEGREE)
        right -= left * np.dot(left, right) / np.dot(left, left)
        scale = np.sqrt(np.dot(left, left) + np.dot(right, right))
        left /= scale
        right /= scale
        mass, slope = _overlap(left, right)
        modes = [_lobes(region, left, right, K_CARRIER) for region in range(2)]
        element = _element(modes[0], modes[1], "plus", "plus")
        predicted = 1j * (BETA0 - A0) * np.exp(-1j * K_CARRIER) * (
            slope + 1j * K_CARRIER * mass
        )
        assert abs(element - predicted) < 1e-12
        assert abs(mass) <= popoviciu + 1e-12
    bare = [
        _lobes(region, PLUS_LEFT, PLUS_RIGHT, K_CARRIER) for region in range(2)
    ]
    phased_right = np.exp(1j * PLUS_RIGHT_PHASE) * np.asarray(PLUS_RIGHT, dtype=float)
    phased = [
        _lobes(region, PLUS_LEFT, phased_right, K_CARRIER) for region in range(2)
    ]
    assert abs(_element(bare[0], bare[1], "plus", "plus") - 0.25) > 0.05
    assert abs(_element(phased[0], phased[1], "plus", "plus") - 0.25) < 1e-12
    factor = np.exp(1j * 0.4)
    plus = [_lobes(region, PLUS_LEFT, PLUS_RIGHT, K_CARRIER) for region in range(2)]
    minus = [
        _lobes(region, factor * np.asarray(PLUS_LEFT), factor * np.asarray(PLUS_RIGHT), -K_CARRIER)
        for region in range(2)
    ]
    off_plus = _element(plus[0], minus[1], "plus", "minus")
    off_minus = _element(minus[0], plus[1], "minus", "plus")
    assert abs(off_minus - np.conjugate(off_plus)) < 1e-12
    defect = sp.simplify(LINK[1, 0] - sp.conjugate(LINK[0, 1]))
    assert defect == sp.Rational(1, 9) + sp.I / 7
    assert sp.simplify(sp.Abs(defect) - sp.sqrt(130) / 63) == 0
    assert np.max(np.abs(np.asarray(MINUS_LEFT) - np.asarray(PLUS_LEFT))) > 1.0


def test_complement_leak_and_common_zero_bar_invariance():
    labels, modes, gram, matrix = _compression()
    symbol = np.array(
        [[-BETA0, -1j * A0], [1j * A0, -BETA0]],
        dtype=np.complex128,
    )
    assert abs(np.linalg.det(symbol) - (BETA0 ** 2 - A0 ** 2)) < 1e-12
    assert BETA0 ** 2 - A0 ** 2 < 0.0
    leaks = []
    for index, (_region, spin) in enumerate(labels):
        norm2 = _image_norm_squared(modes[index], spin)
        captured = float(np.real(np.vdot(matrix[:, index], gram @ matrix[:, index])))
        leaks.append(np.sqrt(max(0.0, 1.0 - captured / norm2)))
    assert min(leaks) > 0.98
    image_gram = np.zeros((6, 6), dtype=np.complex128)
    spins = [spin for _region, spin in labels]
    for i in range(6):
        for j in range(6):
            image_gram[i, j] = _image_pair(modes[i], spins[i], modes[j], spins[j])
    assert np.linalg.norm(image_gram - image_gram.conj().T) < 1e-10
    values, vectors = np.linalg.eig(_as_numpy(_rational_window()[0]))
    for index in range(6):
        coeff = vectors[:, index]
        coeff = coeff / np.linalg.norm(coeff)
        energy = float(np.real(np.vdot(coeff, image_gram @ coeff)))
        captured = float(np.real(np.vdot(matrix @ coeff, gram @ (matrix @ coeff))))
        assert np.sqrt(max(0.0, 1.0 - captured / energy)) > 0.98
        assert abs(values[index].imag) < 1e-8
    # Every stored lobe lives in [0, 4]. The common zero on (4, 8) is exact.
    assert all(lobe["origin"] + 1.0 <= 4.0 for mode in modes for lobe in mode)


def test_fourier_grid_error_decays_as_n_inverse_square():
    window = _as_numpy(_rational_window()[0])
    measured = {}
    for points in (128, 256, 512):
        row = _discrete(points)
        error = float(np.linalg.norm(row["matrix"] - window))
        measured[points] = error
        assert row["operator_disagreement"] == 0.0
        assert row["gram_defect"] < 1e-12
        assert row["outside_mass"] == 0.0
        assert row["bridge_max"] > 0.5
    assert 3.5 < measured[128] / measured[256] < 4.5
    assert 3.5 < measured[256] / measured[512] < 4.5
    scaled = [measured[points] * points ** 2 for points in (128, 256, 512)]
    assert max(scaled) - min(scaled) < 0.05 * min(scaled)


def test_initial_current_matches_and_closed_window_evolution_does_not():
    labels, modes, gram, matrix = _compression()
    window, _placed = _rational_window()
    sign = sp.zeros(6, 1)
    sign[0, 0] = 1 / sp.sqrt(2)
    sign[2, 0] = sp.I / sp.sqrt(2)
    covariance = sp.simplify(sign * sign.H)
    region = sp.diag(1, 1, 0, 0, 0, 0)
    current = sp.simplify(-sp.I * sp.trace(region * (window * covariance - covariance * window)))
    assert current == sp.Rational(1, 4)
    state = np.zeros(6, dtype=np.complex128)
    state[0] = 1.0 / np.sqrt(2.0)
    state[2] = 1j / np.sqrt(2.0)
    numeric = np.outer(state, state.conj())
    commutator = matrix @ numeric - numeric @ matrix
    observed = -1j * np.trace(np.diag([1.0, 1.0, 0, 0, 0, 0]) @ commutator)
    assert abs(observed.real - 0.25) < 1e-12
    assert abs(observed.imag) < 1e-12
    g00 = _image_pair(modes[0], "plus", modes[0], "plus")
    g22 = _image_pair(modes[2], "plus", modes[2], "plus")
    g02 = _image_pair(modes[0], "plus", modes[2], "plus")
    energy = (
        abs(state[0]) ** 2 * g00
        + abs(state[2]) ** 2 * g22
        + np.conjugate(state[0]) * state[2] * g02
        + np.conjugate(state[2]) * state[0] * np.conjugate(g02)
    )
    captured = float(np.real(np.vdot(matrix @ state, gram @ (matrix @ state))))
    leak2 = energy.real - captured
    population_second = -2.0 * leak2
    assert leak2 > 100.0
    assert population_second < -200.0
    # Closed e^{-i J t} keeps this vector in the window, so its subspace
    # population has second derivative 0. The Dirac image does not.
    assert population_second != 0.0
    stationary = np.eye(6) / 2.0 - np.linalg.matrix_power(_as_numpy(window), 3) / (4.0 * 8.0 ** 3)
    frozen = _as_numpy(window)
    assert np.linalg.norm(frozen @ stationary - stationary @ frozen) < 1e-12
