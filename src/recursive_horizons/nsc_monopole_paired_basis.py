"""Orthonormal charged spinor eigenspaces of the existing NSC magnetic spectrum.

Imported spin_c / monopole-harmonic identities only. The locked background
flux is the integer q already used by the magnetic angular tower; -q is the
conjugate charge bundle, not a new magnetic sector. Occupations and spatial
covariances are not assigned here.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import expm


FIXED_MAGNETIC_FLUX = 4
SIGMA1 = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
SIGMA3 = np.diag([1.0, -1.0]).astype(complex)


def _flux(q) -> int:
    if isinstance(q, bool) or not isinstance(q, (int, np.integer)) or q == 0:
        raise ValueError("nonzero integer magnetic flux required")
    return int(q)


def _level_index(n) -> int:
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 0:
        raise ValueError("nonnegative integer DiracS2 level n required")
    return int(n)


def right_weights(q: int) -> tuple[float, float]:
    q = _flux(q)
    return (q - 1) / 2.0, (q + 1) / 2.0


def magnetic_axis(j: float) -> np.ndarray:
    two = int(round(2.0 * j))
    if abs(2.0 * j - two) > 1e-12 or two < 0:
        raise ValueError("j must be a nonnegative integer or half-integer")
    return 0.5 * np.arange(two, -two - 1, -2, dtype=float)


def dirac_s2_level(q: int, n: int) -> dict:
    """Existing NSC/BKW labels: λ²=n(n+|q|), j=(|q|-1)/2+n, d_n as stored."""
    q = _flux(q)
    n = _level_index(n)
    magnitude = abs(q)
    spin = (magnitude - 1) / 2.0 + n
    lam = float(np.sqrt(n * (n + magnitude))) if n else 0.0
    s_a, s_b = right_weights(q)
    if n == 0:
        surviving = "a" if abs(abs(s_a) - spin) <= 1e-12 else "b"
        degeneracy = magnitude
    else:
        surviving = "pair"
        degeneracy = 2 * (magnitude + 2 * n)
    return {
        "flux": q,
        "n": n,
        "j": float(spin),
        "lam": lam,
        "lam_squared": n * (n + magnitude),
        "s_a": float(s_a),
        "s_b": float(s_b),
        "degeneracy": int(degeneracy),
        "m_axis": tuple(float(m) for m in magnetic_axis(spin)),
        "surviving": surviving,
    }


@lru_cache(maxsize=None)
def _spin_jy(two_j: int) -> tuple[np.ndarray, np.ndarray]:
    axis = magnetic_axis(two_j / 2.0)
    plus = np.zeros((axis.size, axis.size), dtype=complex)
    index = {float(m): i for i, m in enumerate(axis)}
    j = two_j / 2.0
    for i, m in enumerate(axis):
        if m + 1.0 <= j + 1e-12:
            plus[index[float(m + 1.0)], i] = np.sqrt((j - m) * (j + m + 1.0))
    jy = (plus - plus.T.conj()) / (2j)
    return axis, jy


@lru_cache(maxsize=None)
def wigner_little_d(j: float, theta: float) -> tuple[np.ndarray, np.ndarray]:
    """d^j_{m,s}(θ)=⟨j m|exp(-i θ J_y)|j s⟩ with the Condon–Shortley J_y."""
    two = int(round(2.0 * j))
    if abs(2.0 * j - two) > 1e-12:
        raise ValueError("j must be an integer or half-integer")
    axis, jy = _spin_jy(two)
    matrix = expm(-1j * float(theta) * jy)
    imag = float(np.max(np.abs(matrix.imag)))
    if imag > 1e-10:
        raise ValueError("Wigner d-matrix left the real convention")
    return axis, np.real_if_close(matrix, tol=1e5)


def _axis_index(axis: np.ndarray, value: float) -> int:
    two = int(round(2.0 * value))
    if abs(2.0*value-two)>1e-12:
        raise ValueError("magnetic index must be an integer or half-integer")
    for i, m in enumerate(axis):
        if int(round(2.0 * m)) == two:
            return i
    raise ValueError("magnetic index is not in the spin-j axis")


def little_d(j: float, m: float, s: float, theta: float) -> complex:
    axis, matrix = wigner_little_d(j, theta)
    return complex(matrix[_axis_index(axis, m), _axis_index(axis, s)])


def little_d_and_theta_derivative(j: float, m: float, s: float, theta: float) -> tuple[complex, complex]:
    two = int(round(2.0 * j))
    axis, jy = _spin_jy(two)
    unitary = wigner_little_d(j, theta)[1]
    row, col = _axis_index(axis, m), _axis_index(axis, s)
    return complex(unitary[row, col]), complex((unitary @ (-1j * jy))[row, col])


def monopole_harmonic(s: float, j: float, m: float, theta: float, phi: float, patch: str = "N") -> complex:
    """North/south monopole harmonics built from the same little-d matrix.

    Y_N^s(j,m)=√((2j+1)/(4π)) exp(-i(m-s)φ) d^j_{m,s}(θ)
    Y_S^s(j,m)=√((2j+1)/(4π)) exp(-i(m+s)φ) d^j_{m,s}(θ)
    """
    if patch not in ("N", "S"):
        raise ValueError("patch must be N or S")
    if abs(s) - j > 1e-12 or abs(m) - j > 1e-12:
        raise ValueError("monopole harmonic requires |m|,|s| ≤ j")
    pref = np.sqrt((2.0 * j + 1.0) / (4.0 * np.pi))
    azimuth = (m - s) if patch == "N" else (m + s)
    return pref * np.exp(-1j * azimuth * phi) * little_d(j, m, s, theta)


def component_transition(s: float, phi: float) -> complex:
    return np.exp(2j * s * phi)


def spinor_transition(q: int, phi: float) -> np.ndarray:
    """Two-spinor gauge map χ_N = exp(i q φ) diag(exp(-i φ), exp(i φ)) χ_S."""
    q = _flux(q)
    return np.exp(1j * q * phi) * np.diag([np.exp(-1j * phi), np.exp(1j * phi)])


def conjugation_phase(m: float, s: float) -> int:
    """Base-harmonic identity Y_{s,j,m}^* = (-1)^{m-s} Y_{-s,j,-m}."""
    order = m - s
    rounded = int(round(order))
    if abs(order - rounded) > 1e-10:
        raise ValueError("m-s must be an integer for the conjugation phase")
    return -1 if rounded % 2 else 1


def charge_map_phase(q: int, m: float) -> int:
    s_a, _ = right_weights(q)
    return (1 if q>0 else -1)*conjugation_phase(m, s_a)


def charged_eigenspinor(
    q: int,
    n: int,
    m: float,
    theta: float,
    phi: float,
    eta: int | None = None,
    patch: str = "N",
) -> np.ndarray:
    """Orthonormal charged Dirac eigenspinor of one (n,m,η) label.

    n≥1: χ_η=(Y_a, η Y_b)/√2 for q>0. The q<0 bundle is the same pair
    rephased by η, χ_η=(η Y_a, Y_b)/√2. n=0 keeps only the surviving weight.
    """
    level = dirac_s2_level(q, n)
    spin, s_a, s_b = level["j"], level["s_a"], level["s_b"]
    if abs(abs(m) - spin) > 1e-12 and abs(m) > spin:
        raise ValueError("m is outside the j multiplet")
    if n == 0:
        if eta is not None:
            raise ValueError("n=0 has one surviving harmonic; eta is not a label")
        if level["surviving"] == "a":
            return np.array([monopole_harmonic(s_a, spin, m, theta, phi, patch), 0.0], dtype=complex)
        return np.array([0.0, monopole_harmonic(s_b, spin, m, theta, phi, patch)], dtype=complex)
    if eta not in (1, -1):
        raise ValueError("nonzero levels require eta=+1 or eta=-1")
    north = monopole_harmonic(s_a, spin, m, theta, phi, patch)
    south = monopole_harmonic(s_b, spin, m, theta, phi, patch)
    if q > 0:
        return np.array([north, eta * south], dtype=complex) / np.sqrt(2.0)
    return np.array([eta * north, south], dtype=complex) / np.sqrt(2.0)


def apply_sigma1_k(spinor: np.ndarray) -> np.ndarray:
    return SIGMA1 @ np.conjugate(np.asarray(spinor, dtype=complex))


def angular_quadrature(n_theta: int = 48, n_phi: int = 64) -> dict:
    if n_theta < 8 or n_phi < 8:
        raise ValueError("resolved spherical quadrature required")
    nodes, theta_weights = leggauss(n_theta)
    thetas = np.arccos(nodes)
    phis = 2.0 * np.pi * np.arange(n_phi) / n_phi
    phi_weights = np.full(n_phi, 2.0 * np.pi / n_phi)
    return {
        "theta": thetas,
        "phi": phis,
        "theta_weight": theta_weights,
        "phi_weight": phi_weights,
        "sphere_weight": theta_weights[:, None] * phi_weights[None, :],
    }


def require_charged_spinor_eigenspace(**given) -> None:
    blocked = ("degeneracy", "copy_count", "multiplicity", "deg", "copies")
    if not given:
        raise ValueError("an orthonormal charged spinor eigenspace is required")
    if set(given) <= set(blocked) or any(key in blocked for key in given):
        raise ValueError(
            "degeneracy/copy counts do not replace orthonormal charged spinor eigenspaces"
        )
    raise ValueError("an orthonormal charged spinor eigenspace is required")


def _harmonic_on_grid(s, j, m, grid, patch):
    pref = np.sqrt((2.0 * j + 1.0) / (4.0 * np.pi))
    azimuth = (m - s) if patch == "N" else (m + s)
    d_theta = np.array([little_d(j, m, s, theta) for theta in grid["theta"]], dtype=complex)
    phase = np.exp(-1j * azimuth * grid["phi"])
    return pref * d_theta[:, None] * phase[None, :]


def _spinor_grid(q, n, m, eta, grid, patch="N"):
    level=dirac_s2_level(q,n)
    values=np.zeros((len(grid["theta"]),len(grid["phi"]),2),complex)
    if n==0:
        component=0 if level["surviving"]=="a" else 1
        values[:,:,component]=_harmonic_on_grid(level["s_a"] if component==0 else level["s_b"],level["j"],m,grid,patch)
    else:
        a=_harmonic_on_grid(level["s_a"],level["j"],m,grid,patch)
        b=_harmonic_on_grid(level["s_b"],level["j"],m,grid,patch)
        values[:,:,0]=(a if q>0 else eta*a)/np.sqrt(2.)
        values[:,:,1]=(eta*b if q>0 else b)/np.sqrt(2.)
    return values


def _flatten(spinor_grid, sphere_weight):
    weight = np.broadcast_to(sphere_weight[:, :, None], spinor_grid.shape)
    return spinor_grid.reshape(-1), weight.reshape(-1)


def paired_mode_overlap(q: int, n: int, m: float | None = None) -> np.ndarray:
    """Direct magnetic-spinor inner product used by the boundary state map."""
    if n<=0:raise ValueError('nonzero angular eigenvalue pair required')
    level=dirac_s2_level(q,n)
    if m is None:m=level['m_axis'][0]
    grid=angular_quadrature(32,32)
    values=[_flatten(_spinor_grid(q,n,m,eta,grid),grid['sphere_weight'])[0] for eta in (1,-1)]
    weight=np.repeat(grid['sphere_weight'].ravel(),2)
    modes=np.asarray(values)
    return (modes.conj()*weight)@modes.T


def orthonormality_residuals(q: int = FIXED_MAGNETIC_FLUX, n_max: int = 12, **quad) -> dict:
    q = _flux(q)
    grid = angular_quadrature(**quad)
    modes = []
    labels = []
    for n in range(1, n_max + 1):
        level = dirac_s2_level(q, n)
        for eta in (1, -1):
            for m in level["m_axis"]:
                spinor = _spinor_grid(q, n, m, eta, grid)
                vector, weight = _flatten(spinor, grid["sphere_weight"])
                modes.append(vector)
                labels.append((n, eta, m, weight))
    psi = np.stack(modes)
    weight = labels[0][3]
    gram = (psi.conj() * weight) @ psi.T
    identity = np.eye(len(modes))
    off = gram - identity
    return {
        "n_modes": len(modes),
        "max_gram_residual": float(np.max(np.abs(off))),
        "max_diagonal_defect": float(np.max(np.abs(np.diag(gram) - 1.0))),
        "max_off_diagonal": float(np.max(np.abs(off - np.diag(np.diag(off))))),
        "min_diagonal": float(np.min(np.real(np.diag(gram)))),
    }


def _ladder_action(j, m, s_a, s_b, theta):
    d_a, d_a_theta = little_d_and_theta_derivative(j, m, s_a, theta)
    d_b, d_b_theta = little_d_and_theta_derivative(j, m, s_b, theta)
    sine = np.sin(theta)
    cosine = np.cos(theta)
    raise_a = -d_a_theta + (s_a * cosine - m) / sine * d_a
    lower_b = d_b_theta + (s_b * cosine - m) / sine * d_b
    return d_a, d_b, raise_a, lower_b


def ladder_residuals(q: int = FIXED_MAGNETIC_FLUX, n_max: int = 12, n_theta: int = 48) -> dict:
    q = _flux(q)
    nodes, weights = leggauss(n_theta)
    thetas = np.arccos(nodes)
    raise_res = 0.0
    eigen_res = 0.0
    norm_res = 0.0
    for n in range(1, n_max + 1):
        level = dirac_s2_level(q, n)
        j, lam, s_a, s_b = level["j"], level["lam"], level["s_a"], level["s_b"]
        for m in level["m_axis"]:
            for eta in (1, -1):
                k_norm_sq = 0.0
                chi_norm_sq = 0.0
                for theta, weight in zip(thetas, weights):
                    d_a, d_b, raise_a, lower_b = _ladder_action(j, m, s_a, s_b, theta)
                    raise_res = max(raise_res, abs(raise_a - lam * d_b), abs(lower_b - lam * d_a))
                    chi = np.array([d_a, eta * d_b], dtype=complex) / np.sqrt(2.0)
                    acted = np.array([eta * lower_b, raise_a], dtype=complex) / np.sqrt(2.0)
                    eigen_res = max(eigen_res, float(np.linalg.norm(acted - eta * lam * chi)))
                    k_norm_sq += weight * float(np.vdot(acted, acted).real)
                    chi_norm_sq += weight * float(np.vdot(chi, chi).real)
                if chi_norm_sq <= 0.0:
                    raise ValueError("vanishing charged spinor quadrature norm")
                norm_res = max(norm_res, abs(np.sqrt(k_norm_sq / chi_norm_sq) - lam))
    return {
        "max_raise_lower_residual": float(raise_res),
        "max_eta_eigen_residual": float(eigen_res),
        "max_unit_ladder_norm_residual": float(norm_res),
    }


def gauge_residuals(q: int = FIXED_MAGNETIC_FLUX, n_max: int = 12, **quad) -> dict:
    q = _flux(q)
    grid = angular_quadrature(**quad)
    harmonic = 0.0
    spinor = 0.0
    for n in range(1, n_max + 1):
        level = dirac_s2_level(q, n)
        for s in (level["s_a"], level["s_b"]):
            for m in level["m_axis"]:
                for theta in grid["theta"][:: max(1, grid["theta"].size // 8)]:
                    for phi in grid["phi"][:: max(1, grid["phi"].size // 8)]:
                        north = monopole_harmonic(s, level["j"], m, theta, phi, "N")
                        south = monopole_harmonic(s, level["j"], m, theta, phi, "S")
                        harmonic = max(harmonic, abs(north - component_transition(s, phi) * south))
        for eta in (1, -1):
            for m in level["m_axis"][:: max(1, len(level["m_axis"]) // 6)]:
                for theta in grid["theta"][:: max(1, grid["theta"].size // 6)]:
                    for phi in grid["phi"][:: max(1, grid["phi"].size // 6)]:
                        left = charged_eigenspinor(q, n, m, theta, phi, eta=eta, patch="N")
                        right = spinor_transition(q, phi) @ charged_eigenspinor(
                            q, n, m, theta, phi, eta=eta, patch="S"
                        )
                        spinor = max(spinor, float(np.linalg.norm(left - right)))
    return {
        "max_component_transition_residual": float(harmonic),
        "max_spinor_transition_residual": float(spinor),
    }


def conjugation_residuals(q: int = FIXED_MAGNETIC_FLUX, n_max: int = 4) -> dict:
    residual = 0.0
    samples = (
        (0.4, 0.3),
        (1.1, 1.7),
        (2.2, 4.2),
        (2.8, 5.5),
    )
    for n in range(0, n_max + 1):
        for flux in (q, -q):
            level = dirac_s2_level(flux, n)
            weights = (level["s_a"], level["s_b"]) if n else ((level["s_a"],) if level["surviving"] == "a" else (level["s_b"],))
            for s in weights:
                if abs(s) - level["j"] > 1e-12:
                    continue
                for m in level["m_axis"]:
                    for patch in ("N", "S"):
                        for theta, phi in samples:
                            left = np.conjugate(monopole_harmonic(s, level["j"], m, theta, phi, patch))
                            right = conjugation_phase(m, s) * monopole_harmonic(
                                -s, level["j"], -m, theta, phi, patch
                            )
                            residual = max(residual, abs(left - right))
    return {"max_harmonic_conjugation_residual": float(residual)}


def charge_map_residuals(q: int = FIXED_MAGNETIC_FLUX, n_max: int = 12) -> dict:
    """σ1 K χ^q_{η,m} vs the stated phase times the rephased conjugate bundle."""
    samples = ((0.5, 0.4), (1.0, 1.2), (2.0, 3.5), (2.6, 5.1))
    stated = 0.0
    inverse = 0.0
    naive_negative = 0.0
    extra_minus = 0.0
    negative_physical = 0.0
    for n in range(1, n_max + 1):
        level = dirac_s2_level(q, n)
        for m in level["m_axis"]:
            for eta in (1, -1):
                phase = charge_map_phase(q, m)
                for theta, phi in samples:
                    source = charged_eigenspinor(q, n, m, theta, phi, eta=eta)
                    image = apply_sigma1_k(source)
                    target = charged_eigenspinor(-q, n, -m, theta, phi, eta=-eta)
                    stated = max(stated, float(np.linalg.norm(image - phase * target)))
                    back = apply_sigma1_k(target)
                    inverse = max(inverse, float(np.linalg.norm(back - phase * source)))
                    negative = charged_eigenspinor(-q, n, -m, theta, phi, eta=-eta)
                    naive = conjugation_phase(-m, right_weights(-q)[0])
                    mapped_neg = apply_sigma1_k(negative)
                    pos = charged_eigenspinor(q, n, m, theta, phi, eta=eta)
                    negative_physical=max(negative_physical,float(np.linalg.norm(
                        mapped_neg-charge_map_phase(-q,-m)*pos)))
                    naive_negative = max(
                        naive_negative, float(np.linalg.norm(mapped_neg - naive * pos))
                    )
                    extra_minus = max(
                        extra_minus, float(np.linalg.norm(mapped_neg + naive * pos))
                    )
    zero = dirac_s2_level(q, 0)
    zero_res = 0.0
    for m in zero["m_axis"]:
        phase = charge_map_phase(q, m)
        for theta, phi in samples:
            source = charged_eigenspinor(q, 0, m, theta, phi)
            image = apply_sigma1_k(source)
            target = charged_eigenspinor(-q, 0, -m, theta, phi)
            zero_res = max(zero_res, float(np.linalg.norm(image - phase * target)))
    return {
        "stated_q_positive_residual": float(stated),
        "antiunitary_inverse_residual": float(inverse),
        "negative_bundle_physical_residual": float(negative_physical),
        "n0_stated_residual": float(zero_res),
        "naive_negative_q_phase_residual": float(naive_negative),
        "naive_negative_q_extra_minus_residual": float(extra_minus),
        "s_a_plus_s_a_conjugate": float(right_weights(q)[0] + right_weights(-q)[0]),
        "physical_map": (
            "locked under the stated q>0 identity and (σ1 K)^2=I; "
            "the same algebraic phase evaluated at q<0 differs by a minus "
            "because s_a(q)+s_a(-q)=-1"
        ),
    }


def paired_sigma3_lock_residuals(q: int = FIXED_MAGNETIC_FLUX, n: int = 1) -> dict:
    """Angular σ1 K tensored with the existing radial σ3 is η1⊗σ3 K up to phase."""
    level = dirac_s2_level(q, n)
    residual = 0.0
    radial = np.array([0.3 - 0.4j, -1.1 + 0.2j], dtype=complex)
    samples = ((0.7, 0.5), (1.8, 2.4))
    for m in level["m_axis"]:
        for eta in (1, -1):
            phase = charge_map_phase(q, m)
            for theta, phi in samples:
                chi = charged_eigenspinor(q, n, m, theta, phi, eta=eta)
                left = np.kron(apply_sigma1_k(chi), SIGMA3 @ np.conjugate(radial))
                target = charged_eigenspinor(-q, n, -m, theta, phi, eta=-eta)
                right = phase * np.kron(target, SIGMA3 @ np.conjugate(radial))
                residual = max(residual, float(np.linalg.norm(left - right)))
    return {
        "max_eta1_sigma3K_lock_residual": float(residual),
        "radial_factor": "existing σ3 K",
        "angular_factor": "σ1 K on the charged spinor eigenspace",
    }


def zero_mode_residuals(q: int = FIXED_MAGNETIC_FLUX) -> dict:
    level = dirac_s2_level(q, 0)
    missing = 0.0
    samples = ((0.6, 0.2), (2.1, 4.0))
    for m in level["m_axis"]:
        for theta, phi in samples:
            spinor = charged_eigenspinor(q, 0, m, theta, phi)
            if level["surviving"] == "a":
                missing = max(missing, abs(spinor[1]))
            else:
                missing = max(missing, abs(spinor[0]))
    conjugate = dirac_s2_level(-q, 0)
    return {
        "surviving_weight": level["surviving"],
        "conjugate_surviving_weight": conjugate["surviving"],
        "vanishing_component_residual": float(missing),
        "degeneracy": level["degeneracy"],
    }


def verify_locked_q4_basis(n_max: int = 12, n_theta: int = 48, n_phi: int = 64) -> dict:
    if n_max < 1:
        raise ValueError("nonzero DiracS2 levels n=1..n_max required")
    return {
        "flux": FIXED_MAGNETIC_FLUX,
        "n_max": n_max,
        "conventions": {
            "Y_N": "sqrt((2j+1)/(4pi))*exp(-i(m-s)phi)*d^j_{m,s}(theta)",
            "Y_S": "sqrt((2j+1)/(4pi))*exp(-i(m+s)phi)*d^j_{m,s}(theta)",
            "little_d": "d^j_{m,s}(theta)=<j m|exp(-i theta J_y)|j s>",
            "component_transition": "exp(2 i s phi)",
            "spinor_transition": "exp(i q phi)*diag(exp(-i phi), exp(i phi))",
            "chi_q_positive": "(Y_a, eta Y_b)/sqrt(2)",
            "chi_q_negative_n_ge_1": "(eta Y_a, Y_b)/sqrt(2)",
            "charge_map": "sigma1 K chi^q_{eta,m}=(-1)^{m-s_a} chi^{-q}_{-eta,-m}",
            "n0": "only the |s|=j harmonic survives",
            "lambda": "sqrt(n(n+|q|))",
        },
        "zero_mode": zero_mode_residuals(),
        "conjugation": conjugation_residuals(n_max=min(4, n_max)),
        "orthonormality": orthonormality_residuals(n_max=n_max, n_theta=n_theta, n_phi=n_phi),
        "ladder": ladder_residuals(n_max=n_max, n_theta=n_theta),
        "gauge": gauge_residuals(n_max=n_max, n_theta=n_theta, n_phi=n_phi),
        "charge_map": charge_map_residuals(n_max=n_max),
        "paired_lock": paired_sigma3_lock_residuals(),
    }
