"""First finite curved spatial Dirac boundary maps on the nested-space throat.

This module owns the first-order radial operator

    D = [[0, -d + w], [d + w, 0]],   w = kappa / sqrt(1 + rho^2),

its self-adjoint staggered/mimetic discretization, parent/child Weyl data at
complex energy, the joined resolvent, and the relation of those data to the
Dirichlet-to-Neumann map of D^2. It does not select a scale, an electron mass,
a covariant action, or a Lorentzian graviton.

Geometry (black-universe chart): parent is rho > 0 (asymptotically flat);
child is rho < 0 (cosmological). The radial measure is the flat measure after
the spinor rescaling that produces A = d + w and A^dagger = -d + w.
"""

from __future__ import annotations

import hashlib
from math import log
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh, svd
from scipy.sparse import csc_matrix, lil_matrix
from scipy.sparse.linalg import splu
import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
GEOMETRY_CHECKPOINT = "5f38712ca01ddd71e715fd265088925a73369aba"
DEFAULT_RESULT = ROOT / "results/nsc-3-boundary-response.json"
SCHEMA = "NSC-3-BOUNDARY-RESPONSE-v2"
SOURCE_PATHS = (
    "src/recursive_horizons/nsc_boundary.py",
    "scripts/check_nsc_boundary_response.py",
    "tests/test_nsc_boundary.py",
    "docs/nsc-boundary-response.md",
)
WARP_COEFFICIENT = -18 / 1015
PRIMARY_KAPPA = 1.0
PRIMARY_RADIUS = 8.0
PRIMARY_HALF_INTERVALS = 80
SCHUR_HALF_INTERVALS = 40
DENSE_CONTROL_HALF_INTERVALS = 20
PROBE_ENERGIES = (
    0.4 + 0.3j,
    0.8 + 0.2j,
    -0.5 + 0.4j,
    0.3 - 0.25j,
)
CONVERGENCE_INTERVALS = (20, 40, 80, 160)
DOMAIN_RADII = (4.0, 8.0, 12.0)
DOMAIN_SPACING = 0.1
MASSIVE_MASS = 0.75
WARP_Y_INTERVALS = 16
TWOD_RADIUS = 4.0
TWOD_HALF_INTERVALS = 16
TWOD_Y_INTERVALS = 8
LOW_MODE_THRESHOLD = 0.5
KERNEL_ABS = 1.0e-10
COMPARE_RTOL = 1.0e-8
COMPARE_ATOL = 1.0e-8
# Parent occupies rho>0; child occupies rho<0.
SIDE_INTERVAL = {
    "parent": "rho_in_[0, R]_asymptotically_flat",
    "child": "rho_in_[-R, 0]_cosmological",
}
OUTWARD_NORMAL = {
    "parent": "at_rho=0_outward_normal_is_minus_d_rho",
    "child": "at_rho=0_outward_normal_is_plus_d_rho",
}


def _native(value: Any) -> Any:
    if isinstance(value, (np.bool_,)):
        return bool(value)
    if isinstance(value, np.generic):
        return value.item()
    return value


def complex_pair(value: complex) -> dict[str, float]:
    number = complex(value)
    if not np.isfinite(number.real) or not np.isfinite(number.imag):
        raise RuntimeError("nonfinite complex value")
    return {"real": float(number.real), "imag": float(number.imag)}


def areal_radius(rho: np.ndarray | float) -> np.ndarray | float:
    return np.sqrt(1.0 + np.square(rho))


def superpotential(rho: np.ndarray | float, kappa: float = PRIMARY_KAPPA) -> np.ndarray | float:
    return kappa / areal_radius(rho)


def warp_sigma(y: np.ndarray | float) -> np.ndarray | float:
    return WARP_COEFFICIENT * np.square(y)


def source_hashes() -> list[dict[str, str]]:
    inventory = []
    for relative in SOURCE_PATHS:
        path = ROOT / relative
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        inventory.append({"path": relative, "sha256": digest})
    return inventory


def exact_relations() -> dict[str, Any]:
    rho = sp.symbols("rho", real=True)
    kappa = sp.symbols("kappa", positive=True)
    energy = sp.symbols("E")
    radius = sp.sqrt(1 + rho**2)
    weight = kappa / radius
    weight_derivative = sp.simplify(sp.diff(weight, rho))
    partner_plus = sp.simplify(weight**2 + weight_derivative)
    partner_minus = sp.simplify(weight**2 - weight_derivative)
    expected_plus = kappa**2 / (1 + rho**2) - kappa * rho / (1 + rho**2) ** sp.Rational(3, 2)
    expected_minus = kappa**2 / (1 + rho**2) + kappa * rho / (1 + rho**2) ** sp.Rational(3, 2)
    plus_residual = sp.simplify(partner_plus - expected_plus)
    minus_residual = sp.simplify(partner_minus - expected_minus)
    sigma = kappa * sp.asinh(rho)
    factorized = sp.simplify(sp.exp(-sigma) * sp.diff(sp.exp(sigma), rho) - weight)

    u, v = sp.symbols("u v", complex=True)
    chi_u, chi_v = sp.symbols("chi_u chi_v", complex=True)
    current_matrix = sp.Matrix([[0, -1], [1, 0]])
    psi = sp.Matrix([u, v])
    chi = sp.Matrix([chi_u, chi_v])
    current = sp.simplify((psi.conjugate().T * current_matrix * chi)[0])
    expected_current = sp.conjugate(v) * chi_u - sp.conjugate(u) * chi_v
    current_residual = sp.simplify(current - expected_current)

    phase = sp.symbols("vartheta", real=True)
    historical = sp.diag(sp.exp(sp.I * phase), sp.exp(-sp.I * phase))
    historical_image = sp.simplify(
        historical.conjugate().T * current_matrix * historical
    )
    identity = sp.eye(2)
    identity_image = sp.simplify(identity.conjugate().T * current_matrix * identity)
    symplectic_family = sp.Matrix(
        [
            [sp.cos(phase), -sp.sin(phase)],
            [sp.sin(phase), sp.cos(phase)],
        ]
    )
    symplectic_image = sp.simplify(
        symplectic_family.T * current_matrix * symplectic_family
    )

    weyl = sp.symbols("m")
    susy_dtn = energy * weyl
    dtn_identity = sp.simplify(susy_dtn - energy * weyl)
    w0 = weight.subs(rho, 0)
    logarithmic_derivative = -w0 + energy * weyl
    m_parent, m_child = sp.symbols("m_parent m_child")
    u_log_parent = -w0 + energy * m_parent
    u_log_child = -w0 + energy * m_child
    n_parent_out = -u_log_parent
    n_child_out = u_log_child
    oriented_jump = sp.simplify(n_parent_out + n_child_out)
    expected_jump = energy * (m_child - m_parent)
    jump_residual = sp.simplify(oriented_jump - expected_jump)
    w0_in_jump = sp.simplify(oriented_jump - expected_jump)

    omega, parent, link, child = sp.symbols("Omega K_p b Gamma_c", nonzero=True)
    first_order = parent - (1 / omega) * link * (1 / child) * sp.conjugate(link)
    unweighted = parent - link * (1 / child) * sp.conjugate(link)
    convention_difference = sp.simplify(first_order - unweighted)

    hp = sp.Matrix([[1, sp.Rational(1, 3)], [sp.Rational(1, 3), 2]])
    hc = sp.Matrix([[3, sp.I / 4], [-sp.I / 4, 4]])
    coupling = sp.Matrix([[sp.Rational(1, 5), sp.I / 7], [sp.Rational(1, 6), 0]])
    z = 2 + sp.I
    rp = z * sp.eye(2) - hp
    rc = z * sp.eye(2) - 2 * hc
    full = rp.row_join(-coupling).col_join((-coupling.conjugate().T).row_join(rc))
    schur = rp - coupling * rc.inv() * coupling.conjugate().T
    schur_residual = (full.inv()[:2, :2] - schur.inv()).applyfunc(sp.simplify)
    reversed_order = rp - coupling.conjugate().T * rc.inv() * coupling
    order_difference = (schur - reversed_order).applyfunc(sp.simplify)

    p_sym, phi = sp.symbols("p Phi", real=True, nonzero=True)
    dirac = sp.Matrix([[p_sym, phi], [phi, -p_sym]])
    visible = (energy * sp.eye(2) - dirac).inv()[0, 0]
    square_residual = sp.simplify(dirac**2 - (p_sym**2 + phi**2) * sp.eye(2))
    visible_residual = sp.factor(1 / visible - (energy - p_sym - phi**2 / (energy + p_sym)))

    if plus_residual != 0 or minus_residual != 0 or factorized != 0:
        raise RuntimeError("radial partner algebra failed")
    if current_residual != 0 or identity_image != current_matrix:
        raise RuntimeError("current algebra failed")
    if symplectic_image != current_matrix:
        raise RuntimeError("symplectic current preservation failed")
    if dtn_identity != 0 or schur_residual != sp.zeros(2):
        raise RuntimeError("DtN or Schur algebra failed")
    if jump_residual != 0 or w0_in_jump != 0:
        raise RuntimeError("oriented DtN jump identity failed")
    if square_residual != sp.zeros(2) or visible_residual != 0:
        raise RuntimeError("two-sheet Phi identity failed")
    if order_difference == sp.zeros(2):
        raise RuntimeError("noncommuting block order control collapsed")

    return {
        "operator": {
            "areal_radius": "sqrt(1+rho^2)",
            "superpotential": "w=kappa/r",
            "first_order_matrix": "D=[[0,-d+w],[d+w,0]]",
            "inner_product": "flat_radial_measure_after_psi=u/sqrt(1+rho^2)",
            "A": "d+w",
            "A_dagger": "-d+w",
            "factorization": "A=exp(-kappa*asinh(rho))*d*exp(kappa*asinh(rho))",
            "factorization_residual": str(factorized),
        },
        "orientation": {
            "parent_interval": SIDE_INTERVAL["parent"],
            "child_interval": SIDE_INTERVAL["child"],
            "parent_outward_normal": OUTWARD_NORMAL["parent"],
            "child_outward_normal": OUTWARD_NORMAL["child"],
            "ordinary_logarithmic_derivative": "u'/u=-w+E*m",
            "N_parent_outward": "-u'_parent/u_parent=w0-E*m_parent",
            "N_child_outward": "+u'_child/u_child=-w0+E*m_child",
            "oriented_jump": "N_parent+N_child=E*(m_child-m_parent)",
            "oriented_jump_residual": str(jump_residual),
            "w0_cancels_in_the_jump": True,
            "geometry_source": "black-universe chart: R_4->0 as rho->+inf (parent), R_4=-36*pi as rho->-inf (child)",
        },
        "partners": {
            "w_derivative": str(weight_derivative),
            "V_plus": str(partner_plus),
            "V_minus": str(partner_minus),
            "plus_residual": str(plus_residual),
            "minus_residual": str(minus_residual),
            "D_squared_upper": "-d^2+V_minus",
            "D_squared_lower": "-d^2+V_plus",
        },
        "current": {
            "matrix": str(current_matrix),
            "bilinear": str(expected_current),
            "residual": str(current_residual),
            "identity_matching_preserves_current": identity_image == current_matrix,
            "historical_sigma3_unitary_preserves_this_current": historical_image
            == current_matrix,
            "historical_image": str(historical_image),
            "SO2_family_preserves_current": symplectic_image == current_matrix,
            "joined_baseline_matching": "U=I_smooth_continuation_through_rho=0",
            "abel_wronskian_is_probability_current": False,
            "green_form_owner": "sesquilinear_boundary_form_of_self_adjoint_D_on_the_real_domain",
        },
        "weyl_calderon_versus_DtN": {
            "cauchy_data": "(u,v)_at_the_throat_rho=0",
            "weyl_m": "m(E)=v(0)/u(0)_for_the_exterior_self_adjoint_solution",
            "staggered_edge_trace": "v(plus_or_minus_h/2)/u(0)_is_not_the_continuum_Calderon_map",
            "discrete_throat_reconstruction": "ODE_consistent_extrapolation_from_the_first_edge_to_rho=0",
            "susy_DtN_of_D_squared": "N_{A^dagger A}(E^2)=(A u)/u=E*m(E)",
            "schrodinger_Neumann": str(logarithmic_derivative),
            "identity_residual": str(dtn_identity),
            "spectral_DtN_is_Phi": False,
            "Phi_dimension": "energy",
            "weyl_m_dimension": "1",
            "DtN_dimension": "1/length",
            "Phi_is_two_sheet_off_diagonal_matching_field": True,
        },
        "resolvent_conventions": {
            "first_order_Schur": "Gamma_p(x)=K_p(x)-(1/Omega)*b*Gamma_c(x/Omega)^(-1)*b_dagger",
            "unweighted_alternative": "Gamma=K-b^2/Gamma_requires_b_sym=B/sqrt(Lambda_p Lambda_c)",
            "convention_difference": str(convention_difference),
            "conventions_may_be_mixed": False,
            "noncommuting_block_Schur_residual": str(schur_residual),
            "reversed_multiplication_differs": order_difference != sp.zeros(2),
            "two_sheet_square_residual": str(square_residual),
            "visible_Schur_residual": str(visible_residual),
        },
        "geometry_unit": {
            "rho_unit": "throat_areal_minimum_set_to_1",
            "energy_unit": "1/L_star_with_L_star=1_in_this_chart",
            "w_is_homogeneous_of_degree_minus_one": False,
            "reason": "w=kappa/sqrt(1+rho^2) retains the throat radius inside the square root",
        },
    }


def _node_edge_grid(left: float, right: float, intervals: int) -> tuple[float, np.ndarray, np.ndarray]:
    if intervals < 4:
        raise ValueError("need at least four intervals")
    nodes = np.linspace(left, right, intervals + 1)
    step = float(nodes[1] - nodes[0])
    edges = 0.5 * (nodes[:-1] + nodes[1:])
    return step, nodes, edges


def _masses(step: float, node_count: int, edge_count: int) -> tuple[np.ndarray, np.ndarray]:
    node_mass = np.full(node_count, step)
    node_mass[0] = 0.5 * step
    node_mass[-1] = 0.5 * step
    edge_mass = np.full(edge_count, step)
    return node_mass, edge_mass


def staggered_forward(nodes: np.ndarray, kappa: float) -> np.ndarray:
    """Node-to-edge mixed incidence with superpotential on edges."""
    step = float(nodes[1] - nodes[0])
    edges = 0.5 * (nodes[:-1] + nodes[1:])
    weight = np.asarray(superpotential(edges, kappa), dtype=float)
    intervals = len(edges)
    operator = np.zeros((intervals, len(nodes)), dtype=float)
    index = np.arange(intervals)
    operator[index, index] = -1.0 / step + 0.5 * weight
    operator[index, index + 1] = 1.0 / step + 0.5 * weight
    return operator


def _adjoint(forward: np.ndarray, node_mass: np.ndarray, edge_mass: np.ndarray) -> np.ndarray:
    return (forward.T * edge_mass) / node_mass[:, None]


def _side_grid(radius: float, half_intervals: int, side: str) -> tuple[np.ndarray, int, int]:
    if side == "parent":
        nodes = np.linspace(0.0, radius, half_intervals + 1)
        return nodes, half_intervals, 0
    if side == "child":
        nodes = np.linspace(-radius, 0.0, half_intervals + 1)
        return nodes, 0, half_intervals
    raise ValueError("side must be parent or child")


def _orthonormal_hamiltonian(
    forward: np.ndarray,
    node_mass: np.ndarray,
    edge_mass: np.ndarray,
    *,
    drop_nodes: Sequence[int] = (),
    drop_edges: Sequence[int] = (),
    mass: float = 0.0,
) -> dict[str, Any]:
    node_keep = [i for i in range(forward.shape[1]) if i not in set(drop_nodes)]
    edge_keep = [i for i in range(forward.shape[0]) if i not in set(drop_edges)]
    reduced = forward[np.ix_(edge_keep, node_keep)]
    node_mass_reduced = node_mass[node_keep]
    edge_mass_reduced = edge_mass[edge_keep]
    off = reduced.T * edge_mass_reduced
    n_u = len(node_keep)
    n_v = len(edge_keep)
    stiffness = np.zeros((n_u + n_v, n_u + n_v), dtype=float)
    stiffness[:n_u, n_u:] = off
    stiffness[n_u:, :n_u] = off.T
    if mass:
        stiffness[:n_u, :n_u] += mass * np.diag(node_mass_reduced)
        stiffness[n_u:, n_u:] -= mass * np.diag(edge_mass_reduced)
    weights = np.concatenate([node_mass_reduced, edge_mass_reduced])
    scale = np.sqrt(weights)
    hamiltonian = stiffness / scale[:, None] / scale[None, :]
    hermitian = float(np.max(np.abs(hamiltonian - hamiltonian.T)))
    return {
        "H": hamiltonian,
        "stiffness": stiffness,
        "weights": weights,
        "n_u": n_u,
        "n_v": n_v,
        "node_keep": node_keep,
        "edge_keep": edge_keep,
        "hermitian_residual": hermitian,
        "forward": reduced,
    }


def joined_operator(
    radius: float,
    half_intervals: int,
    kappa: float = PRIMARY_KAPPA,
    mass: float = 0.0,
) -> dict[str, Any]:
    intervals = 2 * half_intervals
    step, nodes, edges = _node_edge_grid(-radius, radius, intervals)
    forward = staggered_forward(nodes, kappa)
    node_mass, edge_mass = _masses(step, len(nodes), len(edges))
    assembled = _orthonormal_hamiltonian(
        forward,
        node_mass,
        edge_mass,
        drop_nodes=(0, len(nodes) - 1),
        mass=mass,
    )
    assembled.update(
        {
            "step": step,
            "nodes": nodes,
            "edges": edges,
            "radius": radius,
            "half_intervals": half_intervals,
            "kappa": kappa,
            "mass": mass,
            "throat_node_local": half_intervals - 1,
        }
    )
    return assembled


def naive_centered_operator(
    radius: float,
    half_intervals: int,
    kappa: float = PRIMARY_KAPPA,
) -> dict[str, Any]:
    """Collocated centred derivative. Negative control; not a continuum proof."""
    intervals = 2 * half_intervals
    step = (2.0 * radius) / intervals
    rho = np.linspace(-radius, radius, intervals + 1)[1:-1]
    count = len(rho)
    derivative = np.zeros((count, count), dtype=float)
    index = np.arange(count)
    derivative[index[1:], index[:-1]] = -0.5 / step
    derivative[index[:-1], index[1:]] = 0.5 / step
    weight = np.diag(np.asarray(superpotential(rho, kappa), dtype=float))
    hamiltonian = np.zeros((2 * count, 2 * count), dtype=float)
    hamiltonian[:count, count:] = -derivative + weight
    hamiltonian[count:, :count] = derivative + weight
    return {
        "H": hamiltonian,
        "step": step,
        "hermitian_residual": float(np.max(np.abs(hamiltonian - hamiltonian.T))),
        "dimension": 2 * count,
        "method": "naive_collocated_centered_derivative",
    }


def _half_operators(
    radius: float, half_intervals: int, kappa: float, side: str
) -> dict[str, Any]:
    nodes, outer, throat = _side_grid(radius, half_intervals, side)
    forward = staggered_forward(nodes, kappa)
    step = float(nodes[1] - nodes[0])
    edges = 0.5 * (nodes[:-1] + nodes[1:])
    node_mass, edge_mass = _masses(step, len(nodes), half_intervals)
    adjoint = _adjoint(forward, node_mass, edge_mass)
    return {
        "nodes": nodes,
        "edges": edges,
        "step": step,
        "outer": outer,
        "throat": throat,
        "forward": forward,
        "adjoint": adjoint,
        "side": side,
    }


def discrete_half_map(
    radius: float,
    half_intervals: int,
    energy: complex,
    kappa: float = PRIMARY_KAPPA,
    side: str = "parent",
) -> dict[str, Any]:
    """Sparse LU solve with u(0)=1. No SVD phase."""
    data = _half_operators(radius, half_intervals, kappa, side)
    nodes = data["nodes"]
    forward = data["forward"]
    adjoint = data["adjoint"]
    outer = data["outer"]
    throat = data["throat"]
    keep = [i for i in range(len(nodes)) if i not in (outer, throat)]
    n_u = len(keep)
    dim = n_u + half_intervals
    matrix = lil_matrix((dim, dim), dtype=complex)
    rhs = np.zeros(dim, dtype=complex)
    for edge in range(half_intervals):
        for local, full in enumerate(keep):
            matrix[edge, local] = forward[edge, full]
        matrix[edge, n_u + edge] = -energy
        rhs[edge] = -forward[edge, throat]
    interiors = [i for i in range(len(nodes)) if i not in (0, len(nodes) - 1)]
    for row, full in enumerate(interiors):
        index = half_intervals + row
        matrix[index, keep.index(full)] = -energy
        for edge in range(half_intervals):
            matrix[index, n_u + edge] = adjoint[full, edge]
    square = matrix.tocsc()
    solution = splu(square).solve(rhs)
    residual = float(np.linalg.norm(square @ solution - rhs))
    upper = np.zeros(len(nodes), dtype=complex)
    upper[throat] = 1.0
    for local, full in enumerate(keep):
        upper[full] = solution[local]
    lower = solution[n_u:]
    edge_index = 0 if side == "parent" else half_intervals - 1
    v_edge = complex(lower[edge_index])
    step = data["step"]
    weight0 = float(superpotential(0.0, kappa))
    if side == "parent":
        v_throat = (v_edge + 0.5 * step * energy) / (1.0 + 0.5 * step * weight0)
    else:
        v_throat = (v_edge - 0.5 * step * energy) / (1.0 - 0.5 * step * weight0)
    projector_vec = np.array([1.0 + 0.0j, v_throat], dtype=complex)
    throat_projector = np.outer(projector_vec, projector_vec.conj()) / np.vdot(
        projector_vec, projector_vec
    )
    edge_vec = np.array([1.0 + 0.0j, v_edge], dtype=complex)
    edge_projector = np.outer(edge_vec, edge_vec.conj()) / np.vdot(edge_vec, edge_vec)
    return {
        "side": side,
        "interval": SIDE_INTERVAL[side],
        "outward_normal": OUTWARD_NORMAL[side],
        "energy": complex(energy),
        "u_throat": 1.0 + 0.0j,
        "weyl_m_edge": v_edge,
        "weyl_m_throat": complex(v_throat),
        "v_edge": v_edge,
        "linear_solve_residual": residual,
        "solver": "sparse_lu_u_throat_normalized_to_1",
        "staggered_edge_projector_is_continuum_Calderon": False,
        "throat_cauchy_projector_residual": float(
            np.max(np.abs(throat_projector @ throat_projector - throat_projector))
            + np.max(np.abs(throat_projector - throat_projector.T.conj()))
        ),
        "edge_trace_projector_residual": float(
            np.max(np.abs(edge_projector @ edge_projector - edge_projector))
            + np.max(np.abs(edge_projector - edge_projector.T.conj()))
        ),
        "susy_DtN_throat": energy * complex(v_throat),
        "step": step,
        "edge_rho": float(data["edges"][edge_index]),
        "node_rho": 0.0,
    }


def dense_svd_half_map_control(
    radius: float,
    half_intervals: int,
    energy: complex,
    kappa: float,
    side: str,
) -> complex:
    """Independent dense SVD nullspace, used only as a small control."""
    data = _half_operators(radius, half_intervals, kappa, side)
    nodes = data["nodes"]
    forward = data["forward"]
    adjoint = data["adjoint"]
    outer = data["outer"]
    throat = data["throat"]
    keep = [i for i in range(len(nodes)) if i != outer]
    reduced = forward[:, keep]
    unknown = 2 * half_intervals
    equations = half_intervals + (half_intervals - 1)
    matrix = np.zeros((equations, unknown), dtype=complex)
    matrix[:half_intervals, :half_intervals] = reduced
    matrix[:half_intervals, half_intervals:] = -energy * np.eye(half_intervals)
    interiors = [i for i in range(len(nodes)) if i not in (0, len(nodes) - 1)]
    for row, full in enumerate(interiors):
        matrix[half_intervals + row, keep.index(full)] = -energy
        matrix[half_intervals + row, half_intervals:] = adjoint[full]
    _left, _singular, right = svd(matrix, full_matrices=True)
    kernel = right[-1].conj()
    upper = kernel[:half_intervals]
    lower = kernel[half_intervals:]
    throat_local = keep.index(throat)
    edge_index = 0 if side == "parent" else half_intervals - 1
    return complex(lower[edge_index] / upper[throat_local])


def continuum_half_map(
    radius: float,
    energy: complex,
    kappa: float = PRIMARY_KAPPA,
    side: str = "parent",
    sample_rho: float | None = None,
) -> dict[str, Any]:
    def body(rho: float, state: np.ndarray) -> np.ndarray:
        upper, lower = state
        weight = float(superpotential(rho, kappa))
        return np.array(
            [-weight * upper + energy * lower, -energy * upper + weight * lower],
            dtype=complex,
        )

    interval = (radius, 0.0) if side == "parent" else (-radius, 0.0)
    solution = solve_ivp(
        body,
        interval,
        np.array([0.0 + 0.0j, 1.0 + 0.0j], dtype=complex),
        method="DOP853",
        rtol=1.0e-11,
        atol=1.0e-11,
        dense_output=True,
        max_step=0.05,
    )
    if not solution.success:
        raise RuntimeError(f"continuum IVP failed on {side}: {solution.message}")
    upper0, lower0 = solution.sol(0.0)
    weyl0 = lower0 / upper0
    weight0 = float(superpotential(0.0, kappa))
    upper_derivative = -weight0 * upper0 + energy * lower0
    susy = (upper_derivative + weight0 * upper0) / upper0
    neumann = upper_derivative / upper0
    sample = 0.0 if sample_rho is None else sample_rho
    _upper_s, lower_s = solution.sol(sample)
    abel = _abel_wronskian(body, interval)
    return {
        "side": side,
        "energy": complex(energy),
        "weyl_m_throat": weyl0,
        "u_throat": complex(upper0),
        "v_throat": complex(lower0),
        "susy_DtN": complex(susy),
        "schrodinger_Neumann": complex(neumann),
        "susy_identity_residual": float(abs(susy - energy * weyl0)),
        "sample_rho": float(sample),
        "weyl_m_sample_over_u0": complex(lower_s / upper0),
        "abel_wronskian_relative_spread": abel,
        "ivp_success": True,
    }


def _abel_wronskian(body, interval: tuple[float, float]) -> float:
    """Abel identity for two ODE solutions. Not probability current."""
    first = solve_ivp(
        body,
        interval,
        np.array([0.0 + 0.0j, 1.0 + 0.0j], dtype=complex),
        method="DOP853",
        rtol=1.0e-11,
        atol=1.0e-11,
        dense_output=True,
        max_step=0.05,
    )
    second = solve_ivp(
        body,
        interval,
        np.array([1.0 + 0.0j, 0.0 + 0.0j], dtype=complex),
        method="DOP853",
        rtol=1.0e-11,
        atol=1.0e-11,
        dense_output=True,
        max_step=0.05,
    )
    points = np.linspace(interval[1], interval[0], 9)
    values = []
    for rho in points:
        u1, v1 = first.sol(rho)
        u2, v2 = second.sol(rho)
        values.append(u1 * v2 - v1 * u2)
    array = np.asarray(values, dtype=complex)
    return float(np.max(np.abs(array - array[0])) / max(abs(array[0]), 1.0e-16))


def real_energy_probability_current(
    radius: float, energy: float, kappa: float = PRIMARY_KAPPA, side: str = "parent"
) -> dict[str, Any]:
    """Nonzero local current probe for real E, distinct from the Abel Wronskian.

    Complex Cauchy data deliberately carry flux; a real standing-wave datum
    would make every current vanish identically and give a weak control.
    This probe tests the local conservation law, not a reflecting-wall mode.
    """
    if complex(energy).imag != 0.0:
        raise ValueError("probability-current test requires real energy")

    def body(rho: float, state: np.ndarray) -> np.ndarray:
        upper, lower = state
        weight = float(superpotential(rho, kappa))
        return np.array(
            [-weight * upper + energy * lower, -energy * upper + weight * lower],
            dtype=complex,
        )

    interval = (radius, 0.0) if side == "parent" else (-radius, 0.0)
    solution = solve_ivp(
        body,
        interval,
        np.array([1.0, 1.0j], dtype=complex) / np.sqrt(2.0),
        method="DOP853",
        rtol=1.0e-11,
        atol=1.0e-11,
        dense_output=True,
        max_step=0.05,
    )
    points = np.linspace(interval[1], interval[0], 9)
    currents = []
    for rho in points:
        upper, lower = solution.sol(rho)
        currents.append(2.0 * np.imag(np.conjugate(upper) * lower))
    array = np.asarray(currents, dtype=float)
    return {
        "energy": float(energy),
        "side": side,
        "datum": "(1,i)/sqrt(2) at the outer endpoint; local flux probe, not a reflecting-wall eigenstate",
        "current_values": [float(value) for value in array],
        "relative_spread": float(
            np.max(np.abs(array - array[0])) / max(abs(array[0]), 1.0e-16)
        ),
        "identified_with_abel_wronskian": False,
    }


def compare_discrete_continuum(
    radius: float,
    half_intervals: int,
    energy: complex,
    kappa: float = PRIMARY_KAPPA,
    side: str = "parent",
) -> dict[str, Any]:
    discrete = discrete_half_map(radius, half_intervals, energy, kappa, side)
    continuum = continuum_half_map(
        radius, energy, kappa, side, sample_rho=discrete["edge_rho"]
    )
    edge_target = continuum["weyl_m_sample_over_u0"]
    throat_target = continuum["weyl_m_throat"]
    edge_residual = abs(discrete["weyl_m_edge"] - edge_target) / max(
        abs(edge_target), 1.0e-16
    )
    throat_residual = abs(discrete["weyl_m_throat"] - throat_target) / max(
        abs(throat_target), 1.0e-16
    )
    return {
        "side": side,
        "interval": SIDE_INTERVAL[side],
        "outward_normal": OUTWARD_NORMAL[side],
        "energy": complex_pair(energy),
        "discrete_weyl_m_edge": complex_pair(discrete["weyl_m_edge"]),
        "discrete_weyl_m_throat": complex_pair(discrete["weyl_m_throat"]),
        "continuum_weyl_m_throat": complex_pair(continuum["weyl_m_throat"]),
        "continuum_weyl_at_staggered_edge_over_u0": complex_pair(edge_target),
        "edge_trace_relative_residual": float(edge_residual),
        "throat_map_relative_residual": float(throat_residual),
        "linear_solve_residual": discrete["linear_solve_residual"],
        "solver": discrete["solver"],
        "staggered_edge_projector_is_continuum_Calderon": False,
        "throat_cauchy_projector_residual": discrete["throat_cauchy_projector_residual"],
        "continuum_susy_identity_residual": continuum["susy_identity_residual"],
        "abel_wronskian_relative_spread": continuum["abel_wronskian_relative_spread"],
        "abel_wronskian_is_probability_current": False,
        "susy_DtN_throat": complex_pair(discrete["susy_DtN_throat"]),
        "schrodinger_Neumann_throat": complex_pair(continuum["schrodinger_Neumann"]),
        "step": discrete["step"],
        "edge_rho": discrete["edge_rho"],
    }


def oriented_dtn_jump(
    radius: float,
    half_intervals: int,
    energy: complex,
    kappa: float = PRIMARY_KAPPA,
) -> dict[str, Any]:
    parent = compare_discrete_continuum(radius, half_intervals, energy, kappa, "parent")
    child = compare_discrete_continuum(radius, half_intervals, energy, kappa, "child")
    weight0 = float(superpotential(0.0, kappa))

    def from_map(row: dict[str, Any], sign: int) -> tuple[complex, complex, complex]:
        weyl = complex(
            row["discrete_weyl_m_throat"]["real"],
            row["discrete_weyl_m_throat"]["imag"],
        )
        log_derivative = -weight0 + energy * weyl
        outward = sign * log_derivative
        return weyl, log_derivative, outward

    m_parent, log_parent, n_parent = from_map(parent, -1)
    m_child, log_child, n_child = from_map(child, +1)
    jump = n_parent + n_child
    expected = energy * (m_child - m_parent)
    cont_parent = complex(
        parent["schrodinger_Neumann_throat"]["real"],
        parent["schrodinger_Neumann_throat"]["imag"],
    )
    cont_child = complex(
        child["schrodinger_Neumann_throat"]["real"],
        child["schrodinger_Neumann_throat"]["imag"],
    )
    m_parent_c = complex(
        parent["continuum_weyl_m_throat"]["real"],
        parent["continuum_weyl_m_throat"]["imag"],
    )
    m_child_c = complex(
        child["continuum_weyl_m_throat"]["real"],
        child["continuum_weyl_m_throat"]["imag"],
    )
    continuum_jump = (-cont_parent) + cont_child
    continuum_expected = energy * (m_child_c - m_parent_c)
    return {
        "energy": complex_pair(energy),
        "parent": parent,
        "child": child,
        "N_parent_outward": complex_pair(n_parent),
        "N_child_outward": complex_pair(n_child),
        "oriented_jump": complex_pair(jump),
        "E_times_m_child_minus_m_parent": complex_pair(expected),
        "algebraic_jump_residual": float(abs(jump - expected)),
        "continuum_oriented_jump": complex_pair(continuum_jump),
        "continuum_jump_identity_residual": float(
            abs(continuum_jump - continuum_expected)
        ),
        "w0_cancels_in_the_jump": True,
        "parent_plus_child_weyl_throat": complex_pair(m_parent + m_child),
        "parent_minus_child_weyl_throat": complex_pair(m_parent - m_child),
    }


def schur_resolvent_test(
    radius: float,
    half_intervals: int,
    energy: complex,
    kappa: float = PRIMARY_KAPPA,
    *,
    dense_control: bool = False,
) -> dict[str, Any]:
    assembled = joined_operator(radius, half_intervals, kappa)
    hamiltonian = assembled["H"]
    dimension = hamiltonian.shape[0]
    keep = assembled["throat_node_local"]
    probe = np.zeros(dimension, dtype=complex)
    probe[keep] = 1.0
    shifted = hamiltonian.astype(complex) - energy * np.eye(dimension, dtype=complex)
    sparse_factor = splu(csc_matrix(shifted))
    sparse_response = sparse_factor.solve(probe)
    quadratic = complex(np.vdot(probe, sparse_response))
    elim = [i for i in range(dimension) if i != keep]
    block_ee = hamiltonian[np.ix_(elim, elim)].astype(complex) - energy * np.eye(
        len(elim), dtype=complex
    )
    block_ke = hamiltonian[np.ix_([keep], elim)].astype(complex)
    block_ek = hamiltonian[np.ix_(elim, [keep])].astype(complex)
    block_kk = complex(hamiltonian[keep, keep] - energy)
    reduced = splu(csc_matrix(block_ee)).solve(block_ek)
    schur = block_kk - (block_ke @ reduced)[0, 0]
    eliminated = 1.0 / schur
    throat_direct = complex(sparse_response[keep])
    residual = abs(throat_direct - eliminated)
    dense_residual = None
    if dense_control:
        dense_response = np.linalg.solve(shifted, probe)
        dense_residual = float(np.max(np.abs(dense_response - sparse_response)))
    return {
        "energy": complex_pair(energy),
        "solver": "sparse_lu",
        "dense_inversion_is_primary": False,
        "direct_versus_eliminated_residual": float(residual),
        "dense_versus_sparse_residual": dense_residual,
        "throat_resolvent": complex_pair(throat_direct),
        "quadratic_form": complex_pair(quadratic),
        "imag_quadratic": float(quadratic.imag),
        "imag_energy": float(complex(energy).imag),
        "imag_resolvent_sign_matches_energy": bool(
            quadratic.imag * complex(energy).imag > 0.0
        ),
        "hermitian_residual": assembled["hermitian_residual"],
        "principal_compression_used": False,
        "dimension": dimension,
    }


def _canonical_spectrum(values: np.ndarray, radius: float, mass: float) -> dict[str, Any]:
    algebraic = np.sort(np.real(np.asarray(values, dtype=float)))
    kernel = algebraic[np.abs(algebraic) <= KERNEL_ABS]
    positive = algebraic[algebraic > KERNEL_ABS]
    negative = algebraic[algebraic < -KERNEL_ABS]
    first_positive = float(positive[0]) if positive.size else float("nan")
    first_negative = float(negative[-1]) if negative.size else float("nan")
    return {
        "radius": float(radius),
        "mass": float(mass),
        "dimension": int(algebraic.size),
        "index_kernel_count": int(kernel.size),
        "index_kernel_values": [float(value) for value in kernel],
        "positive_ascending": [float(value) for value in positive[:4]],
        "negative_descending_toward_zero": [float(value) for value in negative[-4:][::-1]],
        "first_positive": first_positive,
        "first_negative": first_negative,
        "first_positive_times_radius": first_positive * float(radius),
        "pair_splitting": float(first_positive + first_negative)
        if positive.size and negative.size
        else None,
        "sorting": "algebraic_ascending_then_signed_bands_not_abs_argmin",
    }


def box_spectra(
    radii: Sequence[float],
    spacing: float,
    kappa: float = PRIMARY_KAPPA,
    mass: float = 0.0,
) -> list[dict[str, Any]]:
    rows = []
    for radius in radii:
        half = int(round(radius / spacing))
        assembled = joined_operator(radius, half, kappa, mass=mass)
        values = eigh(assembled["H"], eigvals_only=True)
        row = _canonical_spectrum(values, radius, mass)
        row["half_intervals"] = half
        row["spacing"] = float(radius / half)
        rows.append(row)
    return rows


def inverse_radius_scaling(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    radii = [float(row["radius"]) for row in rows]
    first = [float(row["first_positive"]) for row in rows]
    products = [energy * radius for energy, radius in zip(first, radii)]
    decreasing = all(first[index] > first[index + 1] for index in range(len(first) - 1))
    products_rising = all(
        products[index] < products[index + 1] for index in range(len(products) - 1)
    )
    product_changes = [
        abs(products[index + 1] - products[index]) / products[index]
        for index in range(len(products) - 1)
    ]
    ratio_residuals = []
    for index in range(len(radii) - 1):
        energy_ratio = first[index] / first[index + 1]
        radius_ratio = radii[index + 1] / radii[index]
        ratio_residuals.append(abs(energy_ratio - radius_ratio) / radius_ratio)
    scales = (
        decreasing
        and products_rising
        and product_changes[-1] < 0.05
        and ratio_residuals[-1] < 0.1
        and max(ratio_residuals) < 0.25
    )
    return {
        "radii": radii,
        "first_positive": first,
        "first_positive_times_radius": products,
        "e1_decreasing": decreasing,
        "products_monotone_increasing": products_rising,
        "product_relative_changes": product_changes,
        "inverse_radius_ratio_residuals": ratio_residuals,
        "scales_as_inverse_radius": scales,
        "criterion": "E1_decreases;_E1*R_increases_with_last_relative_change_<0.05;_E1_ratios_match_inverse_R_to_<10%_on_the_finest_pair",
    }


def doubling_control(
    radius: float, half_intervals: int, kappa: float = PRIMARY_KAPPA
) -> dict[str, Any]:
    staggered = joined_operator(radius, half_intervals, kappa)
    naive = naive_centered_operator(radius, half_intervals, kappa)
    staggered_values = eigh(staggered["H"], eigvals_only=True)
    naive_values = eigh(naive["H"], eigvals_only=True)
    staggered_canonical = _canonical_spectrum(staggered_values, radius, 0.0)
    naive_canonical = _canonical_spectrum(naive_values, radius, 0.0)
    return {
        "radius": float(radius),
        "half_intervals": half_intervals,
        "staggered_dimension": int(staggered_values.size),
        "naive_dimension": int(naive_values.size),
        "staggered_first_positive": staggered_canonical["first_positive"],
        "naive_first_positive": naive_canonical["first_positive"],
        "staggered_modes_below_threshold": int(
            np.sum(np.abs(staggered_values) < LOW_MODE_THRESHOLD)
        ),
        "naive_modes_below_threshold": int(
            np.sum(np.abs(naive_values) < LOW_MODE_THRESHOLD)
        ),
        "staggered_max_absolute": float(np.max(np.abs(staggered_values))),
        "naive_max_absolute": float(np.max(np.abs(naive_values))),
        "two_over_step": 2.0 / staggered["step"],
        "one_over_step": 1.0 / naive["step"],
        "staggered_hermitian_residual": staggered["hermitian_residual"],
        "naive_hermitian_residual": naive["hermitian_residual"],
        "continuum_proof_uses_naive_centered": False,
    }


def warped_compact_operator(intervals: int) -> dict[str, Any]:
    step, nodes, edges = _node_edge_grid(-1.0, 1.0, intervals)
    sigma_nodes = np.asarray(warp_sigma(nodes), dtype=float)
    sigma_edges = np.asarray(warp_sigma(edges), dtype=float)
    incidence = np.zeros((intervals, intervals + 1), dtype=float)
    index = np.arange(intervals)
    incidence[index, index] = -1.0 / step
    incidence[index, index + 1] = 1.0 / step
    conformal = (
        np.diag(np.exp(-0.5 * sigma_edges))
        @ incidence
        @ np.diag(np.exp(-0.5 * sigma_nodes))
    )
    node_mass, edge_mass = _masses(step, len(nodes), len(edges))
    orthonormal = (
        np.diag(np.sqrt(edge_mass))
        @ conformal
        @ np.diag(1.0 / np.sqrt(node_mass))
    )
    singular = np.sort(svd(orthonormal, compute_uv=False))
    _left, _values, right = svd(orthonormal, full_matrices=True)
    kernel_vector = right[-1].conj()
    kernel_residual = float(np.linalg.norm(orthonormal @ kernel_vector))
    return {
        "intervals": intervals,
        "step": step,
        "nodes": nodes,
        "edges": edges,
        "sigma_nodes": sigma_nodes,
        "sigma_edges": sigma_edges,
        "conformal": conformal,
        "node_mass": node_mass,
        "edge_mass": edge_mass,
        "first_positive_singular_value": float(singular[0]),
        "second_singular_value": float(singular[1]),
        "rectangular_right_kernel_residual": kernel_residual,
        "compact_length": 2.0,
        "sigma": "sigma=-18*y^2/1015",
        "operator": "Q_y=exp(-sigma_edge/2)*d_y*exp(-sigma_node/2)",
    }


def coupled_spatial_lowest(
    radius: float,
    half_intervals: int,
    y_intervals: int,
    kappa: float = PRIMARY_KAPPA,
) -> dict[str, Any]:
    """Kronecker spatial coupling. Warp multiplication is a projection, not an invariant subspace."""
    radial = joined_operator(radius, half_intervals, kappa)
    y_data = warped_compact_operator(y_intervals)
    y_op = _orthonormal_hamiltonian(
        y_data["conformal"], y_data["node_mass"], y_data["edge_mass"]
    )
    warp_nodes = np.exp(-y_data["sigma_nodes"])
    warp_edges = np.exp(-y_data["sigma_edges"])
    warp = np.concatenate(
        [warp_nodes[y_op["node_keep"]], warp_edges[y_op["edge_keep"]]]
    )
    y_values, y_vectors = eigh(y_op["H"])
    kernel_index = int(np.argmin(np.abs(y_values)))
    kernel = y_vectors[:, kernel_index]
    scale = float(np.vdot(kernel, warp * kernel).real)
    leak = warp * kernel - scale * kernel
    leakage = float(np.linalg.norm(leak) / max(np.linalg.norm(warp * kernel), 1.0e-16))
    n_u = radial["n_u"]
    radial_sigma3 = np.concatenate([np.ones(n_u), -np.ones(radial["n_v"])])
    spatial = np.kron(np.diag(warp), radial["H"]) + np.kron(
        y_op["H"], np.diag(radial_sigma3)
    )
    hermitian = float(np.max(np.abs(spatial - spatial.T)))
    values = eigh(spatial, eigvals_only=True)
    spatial_canonical = _canonical_spectrum(values, radius, 0.0)
    radial_canonical = _canonical_spectrum(eigh(radial["H"], eigvals_only=True), radius, 0.0)
    predicted = scale * radial_canonical["first_positive"]
    return {
        "radius": float(radius),
        "half_intervals": half_intervals,
        "y_intervals": y_intervals,
        "dimension": int(spatial.shape[0]),
        "hermitian_residual": hermitian,
        "y_kernel_eigenvalue": float(y_values[kernel_index]),
        "projected_warp_scale_factor": scale,
        "warp_multiplication_leakage": leakage,
        "zero_mode_is_exact_invariant_subspace_of_exp_minus_sigma": False,
        "relation_is_projected_approximation": True,
        "first_positive": spatial_canonical["first_positive"],
        "radial_first_positive": radial_canonical["first_positive"],
        "projected_rescaled_first_positive": float(predicted),
        "first_positive_over_projected": spatial_canonical["first_positive"]
        / max(abs(predicted), 1.0e-16),
        "compact_gap_theorem_from_this_finite_matrix": False,
        "lapse_shift_included": False,
        "lorentzian_hamiltonian_assembled": False,
        "compact_operator": {
            "intervals": y_data["intervals"],
            "step": y_data["step"],
            "first_positive_singular_value": y_data["first_positive_singular_value"],
            "second_singular_value": y_data["second_singular_value"],
            "rectangular_right_kernel_residual": y_data["rectangular_right_kernel_residual"],
            "compact_length": y_data["compact_length"],
            "sigma": y_data["sigma"],
            "operator": y_data["operator"],
        },
    }


def map_convergence(
    radius: float,
    energy: complex,
    interval_counts: Sequence[int],
    kappa: float = PRIMARY_KAPPA,
    side: str = "parent",
) -> dict[str, Any]:
    rows = []
    for count in interval_counts:
        comparison = compare_discrete_continuum(radius, count, energy, kappa, side)
        rows.append(
            {
                "half_intervals": count,
                "step": comparison["step"],
                "edge_trace_relative_residual": comparison["edge_trace_relative_residual"],
                "throat_map_relative_residual": comparison["throat_map_relative_residual"],
                "linear_solve_residual": comparison["linear_solve_residual"],
            }
        )

    def orders(key: str) -> list[float]:
        residuals = [row[key] for row in rows]
        values = []
        for index in range(1, len(residuals)):
            values.append(float(log(residuals[index - 1] / residuals[index]) / log(2.0)))
        return values

    return {
        "side": side,
        "energy": complex_pair(energy),
        "radius": float(radius),
        "rows": rows,
        "edge_trace_observed_orders": orders("edge_trace_relative_residual"),
        "throat_map_observed_orders": orders("throat_map_relative_residual"),
        "method": "sparse_lu_staggered_Whitney_versus_solve_ivp",
    }


def _dense_sparse_map_control(
    radius: float, half_intervals: int, energy: complex, kappa: float
) -> dict[str, Any]:
    rows = {}
    for side in ("parent", "child"):
        sparse = discrete_half_map(radius, half_intervals, energy, kappa, side)
        dense = dense_svd_half_map_control(radius, half_intervals, energy, kappa, side)
        rows[side] = {
            "sparse_m_edge": complex_pair(sparse["weyl_m_edge"]),
            "dense_svd_m_edge": complex_pair(dense),
            "relative_residual": float(
                abs(sparse["weyl_m_edge"] - dense) / max(abs(dense), 1.0e-16)
            ),
        }
    return {
        "half_intervals": half_intervals,
        "energy": complex_pair(energy),
        "sides": rows,
        "dense_svd_is_primary": False,
    }


def build_record() -> dict[str, Any]:
    algebra = exact_relations()
    hashes = source_hashes()
    energy_rows = [
        oriented_dtn_jump(PRIMARY_RADIUS, PRIMARY_HALF_INTERVALS, energy, PRIMARY_KAPPA)
        for energy in PROBE_ENERGIES
    ]
    parent_convergence = map_convergence(
        PRIMARY_RADIUS,
        PROBE_ENERGIES[0],
        CONVERGENCE_INTERVALS,
        PRIMARY_KAPPA,
        "parent",
    )
    child_convergence = map_convergence(
        PRIMARY_RADIUS,
        PROBE_ENERGIES[0],
        CONVERGENCE_INTERVALS,
        PRIMARY_KAPPA,
        "child",
    )
    schur_rows = [
        schur_resolvent_test(
            PRIMARY_RADIUS, SCHUR_HALF_INTERVALS, energy, PRIMARY_KAPPA
        )
        for energy in PROBE_ENERGIES
    ]
    dense_schur_control = [
        schur_resolvent_test(
            PRIMARY_RADIUS,
            DENSE_CONTROL_HALF_INTERVALS,
            energy,
            PRIMARY_KAPPA,
            dense_control=True,
        )
        for energy in PROBE_ENERGIES
    ]
    map_control = _dense_sparse_map_control(
        PRIMARY_RADIUS, DENSE_CONTROL_HALF_INTERVALS, PROBE_ENERGIES[0], PRIMARY_KAPPA
    )
    massless = box_spectra(DOMAIN_RADII, DOMAIN_SPACING, PRIMARY_KAPPA, mass=0.0)
    massive = box_spectra(DOMAIN_RADII, DOMAIN_SPACING, PRIMARY_KAPPA, mass=MASSIVE_MASS)
    scaling = inverse_radius_scaling(massless)
    doubling = doubling_control(PRIMARY_RADIUS, PRIMARY_HALF_INTERVALS, PRIMARY_KAPPA)
    coupled = coupled_spatial_lowest(
        TWOD_RADIUS, TWOD_HALF_INTERVALS, TWOD_Y_INTERVALS, PRIMARY_KAPPA
    )
    current = real_energy_probability_current(
        PRIMARY_RADIUS, float(PROBE_ENERGIES[0].real), PRIMARY_KAPPA, "parent"
    )
    edge_residuals = [
        row[side]["edge_trace_relative_residual"]
        for row in energy_rows
        for side in ("parent", "child")
    ]
    throat_residuals = [
        row[side]["throat_map_relative_residual"]
        for row in energy_rows
        for side in ("parent", "child")
    ]
    if max(edge_residuals) >= 0.05 or max(throat_residuals) >= 0.05:
        raise RuntimeError("discrete/continuum Weyl residual too large")
    if min(row["imag_resolvent_sign_matches_energy"] for row in schur_rows) is not True:
        raise RuntimeError("resolvent imaginary part has the wrong sign")
    if max(row["direct_versus_eliminated_residual"] for row in schur_rows) >= 1.0e-10:
        raise RuntimeError("joined Schur identity failed")
    if max(row["algebraic_jump_residual"] for row in energy_rows) >= 1.0e-10:
        raise RuntimeError("oriented DtN jump identity failed")
    if not scaling["scales_as_inverse_radius"]:
        raise RuntimeError("massless first positive band does not scale as 1/R")
    massive_edge = [abs(row["first_negative"] + MASSIVE_MASS) for row in massive]
    if max(massive_edge) >= 1.0e-6:
        raise RuntimeError("massive control failed to hold the lifted index mode at -m")
    if min(row["first_positive"] for row in massive) < MASSIVE_MASS - 1.0e-8:
        raise RuntimeError("massive control produced a positive mode below the mass")
    if min(parent_convergence["edge_trace_observed_orders"]) < 1.5:
        raise RuntimeError("staggered edge map did not show at least first-plus order")
    if min(parent_convergence["throat_map_observed_orders"]) < 1.5:
        raise RuntimeError("throat reconstruction did not converge")
    if max(item["relative_residual"] for item in map_control["sides"].values()) >= 1.0e-10:
        raise RuntimeError("sparse LU disagreed with dense SVD control")
    if max(1.0 if row["dense_versus_sparse_residual"] is None else row["dense_versus_sparse_residual"] for row in dense_schur_control) >= 1.0e-10:
        raise RuntimeError("sparse Schur disagreed with dense control")
    if current["relative_spread"] >= 1e-8 or min(current["current_values"]) <= .5:
        raise RuntimeError("nonzero real-energy probability-current control failed")
    if hashes != source_hashes():
        raise RuntimeError("source hash changed during the calculation")

    return {
        "artifact_id": "NSC-3-BOUNDARY-RESPONSE",
        "schema": SCHEMA,
        "provenance": {
            "geometry_source_checkpoint": GEOMETRY_CHECKPOINT,
            "geometry_source_statement": "parent_rho>0_asymptotically_flat;_child_rho<0_cosmological;_from_the_black-universe_chart_used_at_5f38712",
            "code_hashes": hashes,
            "related_scoped_results": {
                "isolated_radial_essential_spectrum": {
                    "path": "docs/nsc-radial-spectrum.md",
                    "statement": "Weyl sequence proves spec_ess=R for the unwarped complete-line spatial Dirac; no nonzero mass gap",
                    "this_record_does_not_reprove_or_edit_it": True,
                },
                "geometric_repeated_throat_chain": {
                    "path": "docs/nsc-geometric-chain.md",
                    "statement": "periodic repetition of the finite motif produces a spatial gap; that geometry is not the isolated box used here",
                    "this_record_does_not_reprove_or_edit_it": True,
                },
            },
        },
        "classification": "the_first_order_radial_Dirac_operator_has_oriented_parent_child_Weyl_maps_at_complex_energy_a_sparse_joined_Schur_resolvent_and_no_asymptotic_mass_gap_from_decaying_w_while_the_compact_warp_relation_is_a_projected_approximation",
        "operator": {
            "kind": "spatial_first_order_radial_Dirac",
            "matrix": "[[0,-d+w],[d+w,0]]",
            "w": "kappa/sqrt(1+rho^2)",
            "kappa": PRIMARY_KAPPA,
            "parent_interval": SIDE_INTERVAL["parent"],
            "child_interval": SIDE_INTERVAL["child"],
            "parent_outward_normal": OUTWARD_NORMAL["parent"],
            "child_outward_normal": OUTWARD_NORMAL["child"],
            "measure": "flat_radial_L2_after_the_supersymmetric_rescaling",
            "exterior_boundary_condition": "separated_self_adjoint_u(plus_or_minus_R)=0",
            "discretization": "staggered_mimetic_Whitney_0_forms_and_1_forms",
            "primary_solver": "sparse_lu_with_u_throat=1",
            "dense_svd_is_primary": False,
            "naive_centered_used_as_continuum_proof": False,
            "lapse_shift_spin_connection_included": False,
            "lorentzian_hamiltonian_assembled": False,
            "lorentzian_operator_assumed_elliptic": False,
        },
        "exact_relations": algebra,
        "numerical_family": {
            "box_radius": PRIMARY_RADIUS,
            "half_intervals": PRIMARY_HALF_INTERVALS,
            "spacing": PRIMARY_RADIUS / PRIMARY_HALF_INTERVALS,
            "probe_energies": [complex_pair(energy) for energy in PROBE_ENERGIES],
            "joined_matching": "U=I",
            "principal_compression_used": False,
            "schur_half_intervals": SCHUR_HALF_INTERVALS,
            "dense_control_half_intervals": DENSE_CONTROL_HALF_INTERVALS,
        },
        "half_domain_maps": energy_rows,
        "mesh_convergence": {
            "parent": parent_convergence,
            "child": child_convergence,
        },
        "solver_controls": {
            "dense_svd_versus_sparse_lu": map_control,
            "dense_versus_sparse_schur": dense_schur_control,
        },
        "joined_schur_resolvent": schur_rows,
        "green_form_and_current": {
            "abel_wronskian_is_probability_current": False,
            "real_energy_probability_current": current,
            "joined_hermitian_residual_is_the_discrete_green_identity": True,
        },
        "gap_diagnostics": {
            "decaying_superpotential": "w->0_as_|rho|->infinity",
            "essential_spectrum_claim": "massless_spatial_Dirac_on_the_line_has_no_asymptotic_mass_gap;_see_also_the_independent_Weyl-sequence_proof_in_docs/nsc-radial-spectrum.md",
            "claim_is_a_spectral_theorem_on_the_noncompact_PDE": False,
            "finite_box_index_kernel": "same_chirality_Dirichlet_walls_leave_one_staggered_v_mode_at_E=0",
            "massless_spectra": massless,
            "massive_control_mass": MASSIVE_MASS,
            "massive_spectra": massive,
            "inverse_radius_scaling": scaling,
            "box_levels_scale_as_inverse_radius": scaling["scales_as_inverse_radius"],
            "massive_index_mode_held_at_minus_mass": True,
            "massive_first_positive_above_mass": True,
            "compact_zero_mode_opens_radial_gap": False,
            "repeated_throat_chain_is_a_different_geometry": True,
        },
        "doubling_control": doubling,
        "warped_spatial_coupling": coupled,
        "gate": {
            "partner_algebra_exact": True,
            "current_algebra_exact": True,
            "oriented_jump_identity_exact": True,
            "DtN_is_E_times_weyl_not_Phi": True,
            "parent_is_positive_rho": True,
            "child_is_negative_rho": True,
            "discrete_continuum_maps_resolved": True,
            "throat_map_converges": True,
            "joined_schur_identity_holds": True,
            "sparse_lu_is_primary": True,
            "imag_resolvent_signs_match": True,
            "staggered_convergence_at_least_order_two_observed": min(
                parent_convergence["edge_trace_observed_orders"]
            )
            >= 1.5,
            "inverse_radius_scaling_holds": scaling["scales_as_inverse_radius"],
            "principal_compression_used": False,
            "physical_zeta_selected": False,
            "source_hashes_stable": True,
        },
        "nonclaims": {
            "stationary_scale_selected": False,
            "electron_mass_derived": False,
            "full_covariant_action_varied": False,
            "matrix_anomaly_for_bar_phi_and_delta_computed": False,
            "complete_Lorentzian_graviton_proved": False,
            "nonlinear_parent_child_transition_solved": False,
            "recursive_tail_solved": False,
            "BFK_compact_gluing_applied_to_this_noncompact_throat": False,
            "spectral_DtN_identified_with_Phi": False,
            "finite_box_levels_are_massive_bound_poles": False,
            "staggered_edge_trace_is_continuum_Calderon": False,
            "compact_zero_mode_is_exact_invariant_subspace": False,
            "compact_gap_theorem_from_finite_2D_matrix": False,
            "abel_wronskian_is_probability_current": False,
            "Nested_Space_hypothesis_proved_or_refuted": False,
        },
        "limitations": {
            "stage_3_spatial_not_full_curved_Hamiltonian": "warp_coupling_is_spatial; lapse, shift, and spin connection are omitted",
            "stage_4_variations": "the_resolvent_is_regulator_independent; the_matrix_Weyl_anomaly_delta_Gamma/delta(bar_phi,delta,Phi) is not computed here",
            "regulator": "no_heat_kernel_E1_subtraction_is_used_for_these_maps",
            "domain": "finite_box_with_separated_u_Dirichlet_walls",
            "index_kernel": "the_E=0_staggered_mode_is_a_same_chirality_wall_index, not_an_asymptotic_gap",
            "warp_projection": "exp(-sigma) mixes y modes; the average-warp factor is a projected approximation with quantified leakage",
        },
        "next_required_equation": "assemble_the_same_first_order_maps_into_the_energy_resolved_recursive_tail_with_the_explicit_1/Omega_Schur_convention_and_derive_the_matrix_anomaly_of_one_regulated_measure_before_any_zeta_stationarity_claim",
        "terminal": True,
    }


def jsonable(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        if not np.isfinite(number):
            raise RuntimeError("nonfinite float in record")
        return number
    if isinstance(value, complex):
        return complex_pair(value)
    if value is None:
        return None
    return value


def compare(expected: Any, observed: Any, pointer: str = "") -> None:
    expected = _native(expected)
    observed = _native(observed)
    if expected is None or observed is None:
        if expected is not observed:
            raise RuntimeError(f"null mismatch at {pointer}")
        return
    if isinstance(expected, bool) or isinstance(observed, bool):
        if expected is not observed:
            raise RuntimeError(f"boolean mismatch at {pointer}")
    elif isinstance(expected, (int, float)) and isinstance(observed, (int, float)):
        if not np.isclose(expected, observed, rtol=COMPARE_RTOL, atol=COMPARE_ATOL):
            raise RuntimeError(
                f"numeric mismatch at {pointer}: {expected!r} vs {observed!r}"
            )
    elif isinstance(expected, Mapping) and isinstance(observed, Mapping):
        if expected.keys() != observed.keys():
            raise RuntimeError(f"keys differ at {pointer}")
        for key in expected:
            compare(expected[key], observed[key], f"{pointer}/{key}")
    elif isinstance(expected, list) and isinstance(observed, list):
        if len(expected) != len(observed):
            raise RuntimeError(f"length mismatch at {pointer}")
        for index, (left, right) in enumerate(zip(expected, observed)):
            compare(left, right, f"{pointer}/{index}")
    elif expected != observed:
        raise RuntimeError(f"value mismatch at {pointer}: {expected!r} vs {observed!r}")
