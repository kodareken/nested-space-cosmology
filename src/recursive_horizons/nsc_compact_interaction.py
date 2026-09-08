"""Bulk torsion-contact overlaps for the declared free NSC compact interval.

Uses published Einstein--Cartan contact/Clifford identities. The overall
torsion stiffness, boundary action and physical state are not selected here.
The returned vertices omit the common factor kappa4_squared/(32*s_T).
"""
from functools import lru_cache
from itertools import combinations

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

from .nsc_boundary import WARP_COEFFICIENT
from .nsc_spinor_bridge import weyl_matrices


@lru_cache(None)
def clifford_channels():
    frame = weyl_matrices()
    gamma = [np.array(g, dtype=complex) for g in frame["gamma"]]
    gamma5 = np.array(frame["gamma5"], dtype=complex)
    gamma_5d = gamma + [1j * gamma5]
    signs = (1, -1, -1, -1, -1)
    triples = tuple((6 * signs[a] * signs[b] * signs[c],
                     gamma_5d[a] @ gamma_5d[b] @ gamma_5d[c])
                    for a, b, c in combinations(range(5), 3))
    axial = tuple((6 * signs[a], gamma[a] @ gamma5) for a in range(4))
    # Three permutations of the fifth index and both orders of a,b.
    tensor = tuple((6 * signs[a] * signs[b], gamma[a] @ gamma[b] @ gamma5)
                   for a, b in combinations(range(4), 2))
    return {"beta": gamma[0], "triples": triples, "axial": axial, "tensor": tensor}


def profiles(y, level):
    """Dimensionless profiles sqrt(L_star)*f on y=Y/L_star in [-1,1]."""
    if not isinstance(level, (int, np.integer)) or isinstance(level, bool) or level < 0:
        raise ValueError("KK level must be a nonnegative integer")
    y = np.asarray(y, dtype=float)
    if level == 0:
        return np.full_like(y, 1 / np.sqrt(2)), np.zeros_like(y)
    phase = level * np.pi * (y + 1) / 2
    return np.cos(phase), np.sin(phase)


@lru_cache(None)
def quadrature(points=64, warp=WARP_COEFFICIENT):
    y, weights = roots_legendre(points)
    sigma = warp * y * y
    einstein_integral_over_length = float(weights @ np.exp(3 * sigma))
    normalized_contact_weights = einstein_integral_over_length * weights * np.exp(-3 * sigma)
    return y, normalized_contact_weights, einstein_integral_over_length


def profile_overlap(levels, chiralities, *, points=64, warp=WARP_COEFFICIENT):
    """I3 integral e^(-3sigma) f_i f_j f_k f_l dY; dimensionless."""
    if len(levels) != 4 or len(chiralities) != 4 or set(chiralities) - {"L", "R"}:
        raise ValueError("four levels and four L/R labels are required")
    y, weights, _ = quadrature(points, warp)
    product = np.ones_like(y)
    for level, side in zip(levels, chiralities):
        product *= profiles(y, level)[side == "R"]
    return float(weights @ product)


def independent_overlap(levels, chiralities, *, warp=WARP_COEFFICIENT):
    """Adaptive scalar integration, independent of the fixed quadrature."""
    i3 = quad(lambda y: np.exp(3 * warp * y * y), -1, 1, epsabs=1e-13, epsrel=1e-13)[0]
    def integrand(y):
        product = np.exp(-3 * warp * y * y)
        for level, side in zip(levels, chiralities):
            product *= profiles(y, level)[side == "R"]
        return float(product)
    return i3 * quad(integrand, -1, 1, epsabs=1e-13, epsrel=1e-13)[0]


def _embedding(y, level):
    fl, fr = profiles(y, level)
    diagonal = np.stack([fl, fl, fr, fr], axis=-1)
    return diagonal[..., :, None] * np.eye(4)


def vertex(levels, *, points=64, warp=WARP_COEFFICIENT, channels="triples"):
    """Ordered coefficient of bar(psi_i)_a psi_j_b bar(psi_k)_c psi_l_d.

    Computes the full spinor tensor; nonzero raw profile overlaps alone do
    not establish an interaction after fermionic antisymmetrization.
    """
    if len(levels) != 4:
        raise ValueError("four mode labels are required")
    y, weights, _ = quadrature(points, warp)
    data = clifford_channels()
    if channels == "decomposed":
        selected = data["axial"] + data["tensor"]
    elif channels in ("triples", "axial", "tensor"):
        selected = data[channels]
    else:
        raise ValueError("unknown Clifford channel")
    embeddings = [_embedding(y, level) for level in levels]
    beta = data["beta"]
    bars = [beta @ e.conj().transpose(0, 2, 1) @ beta for e in embeddings]
    result = np.zeros((4, 4, 4, 4), dtype=complex)
    for coefficient, matrix in selected:
        first = bars[0] @ matrix @ embeddings[1]
        second = bars[2] @ matrix @ embeddings[3]
        result += coefficient * np.einsum("y,yab,ycd->abcd", weights, first, second)
    return result


def fermionic_vertex(levels, **kwargs):
    """Coefficient projected onto antisymmetric barred and unbarred legs.

    In a full sum over all spin/mode indices, the reordered monomial is
    bar(psi_i)_a bar(psi_k)_c psi_j_b psi_l_d. The minus sign accounts
    for moving psi_j across the second barred field; 1/4 is the projector.
    Output axes remain a,b,c,d for comparison with vertex().
    """
    i, j, k, l = levels
    direct = vertex((i, j, k, l), **kwargs)
    unbarred = vertex((i, l, k, j), **kwargs).transpose(0, 3, 2, 1)
    barred = vertex((k, j, i, l), **kwargs).transpose(2, 1, 0, 3)
    both = vertex((k, l, i, j), **kwargs).transpose(2, 3, 0, 1)
    return (-direct + unbarred + barred - both) / 4
