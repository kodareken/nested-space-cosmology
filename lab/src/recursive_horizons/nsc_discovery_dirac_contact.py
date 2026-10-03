"""Dirac first-order contact consistency for a finite Gaussian source.

The curvature-EFT bulk contact and inverse-metric map are reused. The
scalar/collective spherical specialization is not a Dirac stress. Mean
products and Wick connected bilinears stay separate. The renormalized
four-fermion expectation is not installed.
"""
from __future__ import annotations

import ast
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import sympy as sp
from dataclasses import replace
from scipy.linalg import block_diag, expm
from threadpoolctl import threadpool_limits

from .nsc_conformal_adm_source import (
    block_source,
    conformal_source,
    direct_hamiltonian,
    fixed_gaussian_covariance,
)
from .nsc_covariant_operator import SIGMA1, SIGMA2, smooth_metric
from .nsc_curvature_eft import (
    inverse_metric_shift,
    spherical_contact,
    stress_contact,
    tensor_square,
    tensor_trace,
)
from . import nsc_influence as influence


LAB = Path(__file__).resolve().parents[2]
REPO = LAB.parent
CLI = LAB / "scripts/derive_nsc_discovery_dirac_contact.py"
DEFAULT_OUTPUT = LAB / "results/development/nsc-discovery-dirac-contact-v1.json"
SCHEMA = "NSC-DISCOVERY-DIRAC-CONTACT-v1"
OUTPUT_LIMIT_BYTES = 64 * 1024 * 1024
CPU_LIMIT = 20.0
FINITE_DIFFERENCE_STEP = 1e-6
GRID_POINTS = 12
ANGULAR_LABEL = 1
SHIFT_AMPLITUDE = 0.08
ONE_MODE_OCCUPATION = 0.4
ONE_MODE_LEFT = 1.7
ONE_MODE_RIGHT = -0.6
COPY_COUNT = 4
QUARTIC_EPSILON = 0.35
ETA = sp.diag(1, -1, -1, -1)

MISSING_RELATIONSHIP = (
    "The curvature-EFT contact is the scalar of the Hilbert stress in the "
    "convention T_mu_nu = 2/sqrt|g| delta S_m/delta g^{mu nu}. The owned "
    "orthonormal map is rho = F_N/(4 pi q r^2), T_01 = F_beta/(4 pi q^2 r^2), "
    "p_r = -F_q/(4 pi N r^2), p_perp = -F_r/(8 pi N q r), where F_A = "
    "-delta Gamma/delta A after angular integration. The owned canonical "
    "object is the nodal partial of Tr(C H) on one direct-ordering block, "
    "with M=4 kappa multiplying that mean once. Identifying the partial "
    "with F_A requires partial_A B for B = Gamma_heat - Gamma_canonical_sea, "
    "counted once, and the heat-kernel subtraction of the composite :T T:. "
    "Neither is evaluated. The 2|kappa| angular vertex and the angular "
    "current beyond this radial block are not operators here. Until those "
    "relations are fixed, the renormalized contact expectation is not a number."
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_digest(values):
    values = np.ascontiguousarray(values)
    hasher = hashlib.sha256()
    hasher.update(str(values.dtype).encode())
    hasher.update(str(values.shape).encode())
    hasher.update(values.tobytes())
    return hasher.hexdigest()


def _local_modules(path):
    package = LAB / "src/recursive_horizons"
    names = []
    for node in ast.walk(ast.parse(Path(path).read_text())):
        if isinstance(node, ast.ImportFrom):
            if node.level == 1:
                names.extend([node.module] if node.module else [alias.name for alias in node.names])
            elif node.module and node.module.startswith("recursive_horizons"):
                stem = node.module.removeprefix("recursive_horizons").lstrip(".")
                names.extend([stem] if stem else [alias.name for alias in node.names])
        elif isinstance(node, ast.Import):
            names.extend(alias.name.removeprefix("recursive_horizons.") for alias in node.names
                         if alias.name.startswith("recursive_horizons."))
    return [package / (name.replace(".", "/") + ".py") for name in names]


def computational_paths():
    """Byte closure of this calculation. It is not a campaign replay."""
    pending = [Path(__file__).resolve(), CLI.resolve()]
    visited = set()
    while pending:
        path = pending.pop()
        if path in visited:
            continue
        if not path.is_file():
            raise FileNotFoundError("computational source is unavailable: " + str(path))
        visited.add(path)
        for target in _local_modules(path):
            if target.is_file() and target.resolve() not in visited:
                pending.append(target.resolve())
    return tuple(sorted(visited))


def source_bindings():
    return {path.relative_to(REPO).as_posix(): digest(path) for path in computational_paths()}


def authenticate_git(commit, bindings):
    """Require every declared byte to match one explicit frozen commit."""
    if not isinstance(commit, str) or len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit):
        raise ValueError("an explicit full lowercase frozen Git commit is required")
    requests = [commit + ":" + path for path in bindings]
    completed = subprocess.run(
        ["git", "cat-file", "--batch"], cwd=REPO, input=("\n".join(requests) + "\n").encode(),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    data, offset, objects = completed.stdout, 0, {}
    for path in bindings:
        end = data.index(b"\n", offset)
        header = data[offset:end].decode()
        offset = end + 1
        if header.endswith(" missing"):
            raise ValueError("calculation dependency is absent at the frozen commit: " + path)
        object_id, kind, size = header.split()
        if kind != "blob":
            raise ValueError("calculation dependency is not a Git blob: " + path)
        blob = data[offset:offset + int(size)]
        offset += int(size) + 1
        if hashlib.sha256(blob).hexdigest() != bindings[path]:
            raise ValueError("current calculation bytes differ from the frozen commit: " + path)
        objects[path] = object_id
    return {"commit": commit, "all_declared_bytes_authenticated": True, "blob_ids": objects}


def _symmetric(prefix):
    fields = sp.symbols(prefix + "0:10", real=True)
    matrix, cursor = sp.zeros(4), 0
    for i in range(4):
        for j in range(i, 4):
            matrix[i, j] = matrix[j, i] = fields[cursor]
            cursor += 1
    return matrix


def dirac_orthonormal_stress(rho, radial, angular, current):
    """Covariant orthonormal components in the curvature-EFT contraction frame.

    T_00=rho, T_11=p_r, T_22=T_33=p_perp, T_01=current. These four entries
    are independent. The collective relations are not substituted.
    """
    return sp.Matrix([
        [rho, current, 0, 0],
        [current, radial, 0, 0],
        [0, 0, angular, 0],
        [0, 0, 0, angular],
    ])


def contact_polynomial():
    """Owned bulk contact on an unconstrained Dirac-shaped stress."""
    rho, current, radial, angular = sp.symbols("rho current p_radial p_angular", real=True)
    coefficient_a, c_w, c_r = sp.symbols("A cW cR", real=True, nonzero=True)
    stress = dirac_orthonormal_stress(rho, radial, angular, current)
    square = tensor_square(stress, ETA)
    trace = tensor_trace(stress, ETA)
    contact = stress_contact(stress, ETA, coefficient_a, c_w, c_r)
    split = sp.factor(sp.expand(
        contact + c_w / (2 * coefficient_a**2) * (square - trace**2 / 3)
        + c_r / (4 * coefficient_a**2) * trace**2))
    if split != 0:
        raise RuntimeError("Dirac-shaped stress left the owned contact formula")
    return {
        "contact": str(sp.factor(contact)),
        "trace": str(sp.factor(trace)),
        "trace_equals_rho_minus_p_radial_minus_two_p_angular": str(sp.factor(trace - (rho - radial - 2 * angular))),
        "formula_split_gap": "0",
        "scalar_collective_relations_imposed": False,
    }


def off_shell_identity():
    """The general order-reduction identity on a Dirac-shaped stress."""
    coefficient_a, c_w, c_r = sp.symbols("A cW cR", real=True, nonzero=True)
    rho, current, radial, angular = sp.symbols("rho current p_radial p_angular", real=True)
    einstein = _symmetric("Einstein_residual_")
    stress = dirac_orthonormal_stress(rho, radial, angular, current)
    residual = einstein + stress / (2 * coefficient_a)
    scalar = -tensor_trace(residual, ETA)
    ricci = residual + scalar * ETA / 2
    shift = inverse_metric_shift(ricci, stress, ETA, coefficient_a, c_w, c_r)
    curvature = -2 * c_w * (tensor_square(ricci, ETA) - scalar**2 / 3) - c_r * scalar**2
    action_change = -coefficient_a * sp.trace(einstein * shift)
    gap = sp.factor(sp.expand(curvature + action_change - stress_contact(
        stress, ETA, coefficient_a, c_w, c_r)))
    if gap != 0:
        raise RuntimeError("off-shell Dirac-shaped identity failed")
    return {
        "gap": "0",
        "metric_map": "inverse_metric_shift",
        "boundary_variation_evaluated": False,
        "boundary_variation_discarded": False,
        "euler_box_and_boundary_counted_again": False,
        "euler_shift": "cEuler -> cEuler+cW once, in the inherited static-energy convention",
        "jacobian_set_to_one": False,
        "ghost_initial_data_added": False,
    }


def projection_identity():
    """Solve the owned inverse-metric variation for orthonormal components."""
    lapse, radial_scale, radius, shift = sp.symbols("N q r beta", positive=True)
    rho, current, p_radial, p_angular = sp.symbols("rho j p_radial p_angular", real=True)
    metric = sp.Matrix([
        [lapse**2 - radial_scale**2 * shift**2, -radial_scale**2 * shift],
        [-radial_scale**2 * shift, -radial_scale**2],
    ])
    coframe = sp.Matrix([[lapse, 0], [radial_scale * shift, radial_scale]])
    stress = coframe.T * sp.Matrix([[rho, current], [current, p_radial]]) * coframe
    inverse = sp.simplify(metric.inv())
    volume = 4 * sp.pi * lapse * radial_scale * radius**2

    def variation(field):
        return -volume / 2 * sum(
            stress[i, k] * sp.diff(inverse[i, k], field) for i in range(2) for k in range(2))

    expected = {
        "N": (lapse, 4 * sp.pi * radial_scale * radius**2 * rho),
        "beta": (shift, 4 * sp.pi * radial_scale**2 * radius**2 * current),
        "q": (radial_scale, -4 * sp.pi * lapse * radius**2 * p_radial),
    }
    gaps = {}
    for name, (field, target) in expected.items():
        gap = sp.factor(sp.simplify(variation(field) - target))
        gaps[name] = "0" if gap == 0 else str(gap)
    radial_gap = sp.factor(sp.simplify(
        -volume / 2 * (4 * p_angular / radius) + 8 * sp.pi * lapse * radial_scale * radius * p_angular))
    gaps["r"] = "0" if radial_gap == 0 else str(radial_gap)
    if any(value != "0" for value in gaps.values()):
        raise RuntimeError("orthonormal projection left the owned metric variation")
    return {
        "gaps": gaps,
        "source": "lab/docs/nsc-adm-source-constraints.md and lab/scripts/derive_nsc_adm_source.py",
        "convention": "F_A = -sqrt|g|/(2) T_mu_nu partial g^{mu nu}/partial A",
        "solved": {
            "rho": "F_N/(4*pi*q*r^2)",
            "current": "F_beta/(4*pi*q^2*r^2)",
            "p_radial": "-F_q/(4*pi*N*r^2)",
            "p_angular": "-F_r/(8*pi*N*q*r)",
        },
    }


def orthonormal_from_force_densities(force_n, force_q, force_r, force_beta, lapse, radial_scale, radius):
    """Algebraic image of force densities. It does not decide what F_A is."""
    lapse = np.asarray(lapse, dtype=float)
    radial_scale = np.asarray(radial_scale, dtype=float)
    radius = np.asarray(radius, dtype=float)
    if np.min(lapse) <= 0 or np.min(radial_scale) <= 0 or np.min(radius) <= 0:
        raise ValueError("positive lapse, radial scale and radius are required")
    return {
        "rho": np.asarray(force_n, dtype=float) / (4 * np.pi * radial_scale * radius**2),
        "current": np.asarray(force_beta, dtype=float) / (4 * np.pi * radial_scale**2 * radius**2),
        "p_radial": -np.asarray(force_q, dtype=float) / (4 * np.pi * lapse * radius**2),
        "p_angular": -np.asarray(force_r, dtype=float) / (8 * np.pi * lapse * radial_scale * radius),
    }


def specialization_witness():
    """A pure radial current is outside the collective contact family."""
    coefficient_a, c_w, c_r = sp.symbols("A cW cR", real=True, nonzero=True)
    area_coordinate, field, kinetic = sp.symbols("X chi K", real=True)
    vacuum, charge, flux = sp.symbols("V C_F q_mag", real=True)
    collective = spherical_contact(
        area_coordinate, field, kinetic, coefficient_a, vacuum, charge, flux, c_w, c_r)
    zero_collective = sp.simplify(collective.subs({
        area_coordinate: 2, field: 0, kinetic: 0, coefficient_a: 1,
        vacuum: 0, charge: 1, flux: 0, c_w: 1, c_r: -2}))
    current = dirac_orthonormal_stress(0, 0, 0, 1)
    pure_current = sp.simplify(stress_contact(current, ETA, 1, 1, -2))
    if zero_collective != 0 or pure_current == 0:
        raise RuntimeError("collective specialization witness did not separate")
    return {
        "scalar_collective_specialization_used_as_dirac_stress": False,
        "zero_collective_source_contact": str(zero_collective),
        "pure_radial_current_contact": str(pure_current),
        "difference": str(sp.simplify(pure_current - zero_collective)),
        "reason": "T_01 can be nonzero while the collective kinetic, charge and vacuum scalars vanish",
    }


def _complex_parts(value):
    number = complex(value)
    return {"real": float(number.real), "imag": float(number.imag)}


def one_mode_oracle():
    """Exact one-mode number state, independent of the Wick trace."""
    occupation, left, right = ONE_MODE_OCCUPATION, ONE_MODE_LEFT, ONE_MODE_RIGHT
    density = np.diag([1 - occupation, occupation]).astype(complex)
    number = np.diag([0.0, 1.0]).astype(complex)
    full = complex(np.trace(density @ (left * number) @ (right * number)))
    mean_left = complex(np.trace(density @ (left * number)))
    mean_right = complex(np.trace(density @ (right * number)))
    product = mean_left * mean_right
    wick = complex(influence.connected(
        np.array([[occupation]]), np.array([[left]]), np.array([[right]])))
    connected_exact = full - product
    if abs(connected_exact - occupation * (1 - occupation) * left * right) > 1e-12:
        raise RuntimeError("one-mode connected contraction left n(1-n)")
    if abs(connected_exact - wick) > 1e-12:
        raise RuntimeError("one-mode Wick trace missed the Fock oracle")
    return {
        "occupation": occupation,
        "full_second_moment": _complex_parts(full),
        "product_of_means": _complex_parts(product),
        "connected": _complex_parts(connected_exact),
        "means_entered_the_connected_piece": False,
        "trace_gap": float(abs(connected_exact - wick)),
    }


def _copy_noise(covariance, vertex, copies):
    covariance = np.asarray(covariance, dtype=complex)
    vertex = np.asarray(vertex, dtype=complex)
    one = influence.connected(covariance, vertex, vertex)
    summed = influence.connected(
        block_diag(*([covariance] * copies)), block_diag(*([vertex] * copies)),
        block_diag(*([vertex] * copies)))
    scaled = influence.connected(covariance, copies * vertex, copies * vertex)
    return one, summed, scaled


def multiplicity_noise():
    """Independent copies scale noise by M. A scaled observable scales it by M^2."""
    covariance = np.array([[ONE_MODE_OCCUPATION]], dtype=complex)
    vertex = np.array([[ONE_MODE_LEFT]], dtype=complex)
    one, summed, scaled = _copy_noise(covariance, vertex, COPY_COUNT)
    if abs(summed - COPY_COUNT * one) > 1e-12 or abs(scaled - COPY_COUNT**2 * one) > 1e-12:
        raise RuntimeError("angular-copy noise scaling left M and M^2")
    return {
        "copies": COPY_COUNT,
        "one_copy_connected": _complex_parts(one),
        "independent_product_sum_connected": _complex_parts(summed),
        "single_observable_scaled_by_M_connected": _complex_parts(scaled),
        "product_sum_over_one": float((summed / one).real),
        "scaled_over_one": float((scaled / one).real),
        "declared_default": "independent-copy sum scales by M; M^2 is a different state",
        "noise_set_to_M_squared_without_a_declared_state": False,
    }


def two_mode_wick_oracle():
    covariance = np.array([[0.62, 0.07 + 0.04j], [0.07 - 0.04j, 0.28]], dtype=complex)
    left = np.array([[0.4, 0.25 - 0.1j], [0.25 + 0.1j, -0.3]], dtype=complex)
    right = np.array([[-0.2, 0.15j], [-0.15j, 0.55]], dtype=complex)
    operators = influence.fock_annihilators(2)
    density = influence.gaussian_fock_state(covariance, operators)
    left_operator = influence.second_quantize(left, operators)
    right_operator = influence.second_quantize(right, operators)
    full = complex(np.trace(density @ left_operator @ right_operator))
    product = influence.mean(covariance, left) * influence.mean(covariance, right)
    connected = influence.connected(covariance, left, right)
    gap = full - product - connected
    if abs(gap) > 1e-12:
        raise RuntimeError("two-mode Fock oracle missed Tr C V (I-C) W")
    return {
        "full_second_moment": _complex_parts(full),
        "product_of_means": _complex_parts(product),
        "connected": _complex_parts(connected),
        "gap": float(abs(gap)),
        "covariance": covariance,
        "left": left,
        "right": right,
    }


def _quartic_inputs():
    modes = 3
    index = np.arange(modes)
    raw = np.exp(1j * (index[:, None] * 1.3 - index[None, :] * 0.7)) * (
        0.15 + 0.05 * (index[:, None] + index[None, :]))
    raw = 0.5 * (raw + raw.conj().T)
    values, vectors = np.linalg.eigh(raw)
    del values
    occupation = 0.25 + 0.45 * (index + 1) / (modes + 2)
    covariance = (vectors * occupation) @ vectors.conj().T
    covariance = 0.5 * (covariance + covariance.conj().T)
    vertex_raw = (index[:, None] + 1) + 1j * (index[None, :] - index[:, None])
    vertex = 0.5 * (vertex_raw + vertex_raw.conj().T)
    return covariance, vertex


def quartic_witness():
    """A bilinear square does not stay on the Gaussian family, and its mean does not replace it."""
    covariance, vertex = _quartic_inputs()
    operators = influence.fock_annihilators(3)
    density = influence.gaussian_fock_state(covariance, operators)
    bilinear = influence.second_quantize(vertex, operators)
    quartic = bilinear @ bilinear
    phase = expm(-1j * QUARTIC_EPSILON * quartic)
    evolved = phase @ density @ phase.conj().T

    def one_body(state):
        matrix = np.zeros((3, 3), dtype=complex)
        for i in range(3):
            for j in range(3):
                matrix[i, j] = np.trace(state @ operators[j].conj().T @ operators[i])
        return 0.5 * (matrix + matrix.conj().T)

    evolved_covariance = one_body(evolved)
    rebuilt = influence.gaussian_fock_state(evolved_covariance, operators)
    nongaussian = float(np.linalg.norm(evolved - rebuilt))
    moved = float(np.linalg.norm(evolved - density))
    expectation = complex(np.trace(density @ quartic))
    scalar_generator = expectation * np.eye(evolved.shape[0])
    scalar_phase = expm(-1j * QUARTIC_EPSILON * scalar_generator)
    scalar_state = scalar_phase @ density @ scalar_phase.conj().T
    scalar_distance = float(np.linalg.norm(scalar_state - density))
    mean_bilinear = complex(np.trace(density @ bilinear))
    hartree_phase = expm(-1j * QUARTIC_EPSILON * mean_bilinear * bilinear)
    hartree_state = hartree_phase @ density @ hartree_phase.conj().T
    hartree_covariance = one_body(hartree_state)
    hartree_rebuilt = influence.gaussian_fock_state(hartree_covariance, operators)
    if nongaussian < 1e-8 or moved < 1e-8:
        raise RuntimeError("quartic witness did not leave the Gaussian family")
    if scalar_distance > 1e-12:
        raise RuntimeError("scalar contact expectation moved the density matrix")
    if float(np.linalg.norm(hartree_state - hartree_rebuilt)) > 1e-10:
        raise RuntimeError("mean-field bilinear left the Gaussian family")
    if float(np.linalg.norm(hartree_covariance - evolved_covariance)) < 1e-8:
        raise RuntimeError("mean-field bilinear reproduced the quartic covariance")
    return {
        "modes": 3,
        "epsilon": QUARTIC_EPSILON,
        "generator": "(psi dagger V psi)^2",
        "mean_field_generator": "<psi dagger V psi> (psi dagger V psi)",
        "scalar_generator": "<quartic> I",
        "exact_state_minus_gaussian_rebuild": nongaussian,
        "exact_state_minus_initial": moved,
        "scalar_expectation_moves_state": scalar_distance,
        "hartree_state_minus_its_gaussian": float(np.linalg.norm(hartree_state - hartree_rebuilt)),
        "hartree_covariance_minus_exact": float(np.linalg.norm(hartree_covariance - evolved_covariance)),
        "four_fermion_installed_as_gaussian_mean_force": False,
        "gaussian_preserved_by_the_quartic": False,
        "P_minus_used": False,
    }


def _bump(values, index, delta):
    updated = np.array(values, dtype=float)
    updated[index] = updated[index] + delta
    return updated


def _energy_and_hamiltonian(metric, shift, covariance, kappa):
    hamiltonian = direct_hamiltonian(metric, shift, kappa)
    energy = complex(np.trace(covariance @ hamiltonian))
    if abs(energy.imag) > 1e-10:
        raise RuntimeError("direct-ordering energy left the reals")
    return float(energy.real), hamiltonian


def _varied(metric, shift, kind, index, scale):
    if kind == "lapse":
        return replace(metric, lapse=_bump(metric.lapse, index, scale)), np.array(shift, dtype=float)
    if kind == "Q":
        step = scale * float(metric.sphere_radius[index])
        return replace(metric, radial_scale=_bump(metric.radial_scale, index, step)), np.array(shift, dtype=float)
    if kind == "r":
        return replace(metric, sphere_radius=_bump(metric.sphere_radius, index, scale)), np.array(shift, dtype=float)
    if kind == "shift":
        updated = np.array(shift, dtype=float)
        updated[index] = updated[index] + scale
        return metric, updated
    raise ValueError("geometric variation must be lapse, Q, r or shift")


def _derivative(metric, shift, covariance, kappa, kind, index):
    plus_metric, plus_shift = _varied(metric, shift, kind, index, FINITE_DIFFERENCE_STEP)
    minus_metric, minus_shift = _varied(metric, shift, kind, index, -FINITE_DIFFERENCE_STEP)
    plus_energy, plus_h = _energy_and_hamiltonian(plus_metric, plus_shift, covariance, kappa)
    minus_energy, minus_h = _energy_and_hamiltonian(minus_metric, minus_shift, covariance, kappa)
    step = 2 * FINITE_DIFFERENCE_STEP
    kernel = (plus_h - minus_h) / step
    antisym = float(np.max(np.abs(kernel - kernel.conj().T)))
    kernel = 0.5 * (kernel + kernel.conj().T)
    return (plus_energy - minus_energy) / step, kernel, antisym


def geometric_ledger(metric, shift, covariance):
    """Actual Q, r, lapse and shift variations of the owned one-block trace."""
    owned = block_source(metric, shift, covariance, ANGULAR_LABEL)
    multiplied = conformal_source(metric, shift, covariance, ANGULAR_LABEL)
    names = {"lapse": "N", "Q": "Q", "r": "r", "shift": "beta"}
    analytic = {kind: np.asarray(owned["nodal"][key], dtype=float) for kind, key in names.items()}
    gaps = {}
    kernels = {}
    antisymmetry = {}
    for kind in names:
        derivative = np.empty(metric.points, dtype=float)
        chosen = int(np.argmax(np.abs(analytic[kind])))
        for index in range(metric.points):
            if index == chosen:
                derivative[index], kernels[kind], antisymmetry[kind] = _derivative(
                    metric, shift, covariance, ANGULAR_LABEL, kind, index)
            else:
                derivative[index], _, _ = _derivative(
                    metric, shift, covariance, ANGULAR_LABEL, kind, index)
        gaps[kind] = derivative - analytic[kind]
        if antisymmetry[kind] > 1e-8:
            raise RuntimeError("geometric kernel left the Hermitian direct ordering")
        if np.max(np.abs(gaps[kind])) > 1e-8:
            raise RuntimeError("finite difference left the owned nodal derivative: " + kind)
    projected = orthonormal_from_force_densities(
        analytic["lapse"] / metric.spacing,
        np.asarray(owned["nodal"]["q"], dtype=float) / metric.spacing,
        analytic["r"] / metric.spacing,
        analytic["shift"] / metric.spacing,
        metric.lapse, metric.radial_scale, metric.sphere_radius)
    q_chain = np.asarray(owned["nodal"]["q"], dtype=float) - (
        np.asarray(owned["nodal"]["Q"], dtype=float) / metric.sphere_radius)
    if np.max(np.abs(q_chain)) > 1e-10:
        raise RuntimeError("owned conformal chain left F_q = F_Q / r")
    return owned, multiplied, analytic, gaps, kernels, antisymmetry, projected


def wick_geometric(covariance, kernels, analytic):
    """Means and connected pieces of the varied one-block bilinears."""
    rows = {}
    for kind, kernel in kernels.items():
        mean = influence.mean(covariance, kernel)
        connected = influence.connected(covariance, kernel, kernel)
        product = mean * np.conjugate(mean)
        if abs(mean - analytic[kind][int(np.argmax(np.abs(analytic[kind])))]) > 1e-8:
            raise RuntimeError("kernel trace left the finite-difference derivative")
        rows[kind] = {
            "node": int(np.argmax(np.abs(analytic[kind]))),
            "mean": _complex_parts(mean),
            "product_of_means": _complex_parts(product),
            "connected": _complex_parts(connected),
            "means_added_into_connected": False,
        }
    cross = influence.connected(covariance, kernels["lapse"], kernels["shift"])
    mean_lapse = influence.mean(covariance, kernels["lapse"])
    mean_shift = influence.mean(covariance, kernels["shift"])
    rows["lapse_shift_cross"] = {
        "connected": _complex_parts(cross),
        "product_of_means": _complex_parts(mean_lapse * mean_shift),
        "same_node": rows["lapse"]["node"] == rows["shift"]["node"],
        "cross_contraction_required_to_be_real": False,
        "note": "each kernel is taken at its own largest analytic node",
    }
    return rows


def _sign_map(points):
    unitary = np.kron(SIGMA2, np.eye(points, dtype=complex))
    if np.max(np.abs(unitary.conj().T @ unitary - np.eye(unitary.shape[0]))) > 1e-12:
        raise RuntimeError("sigma2 sign map is not unitary")
    return unitary


def angular_completion(metric, shift, covariance, hamiltonian):
    """Optional block-diagonal sign pair. The map is explicit and not unique."""
    unitary = _sign_map(metric.points)
    negative = direct_hamiltonian(metric, shift, -ANGULAR_LABEL)
    conjugation_gap = float(np.max(np.abs(unitary @ hamiltonian @ unitary.conj().T - negative)))
    if conjugation_gap > 1e-10:
        raise RuntimeError("sigma2 sign map did not send H(kappa) to H(-kappa)")
    wrong = np.kron(SIGMA1, np.eye(metric.points, dtype=complex))
    wrong_gap = float(np.max(np.abs(wrong @ hamiltonian @ wrong.conj().T - negative)))
    if wrong_gap < 1e-6:
        raise RuntimeError("a different Pauli map accidentally reproduced the sign pair")
    generator = expm(1j * hamiltonian / np.linalg.norm(hamiltonian))
    alternate = unitary @ generator
    alternate_gap = float(np.max(np.abs(alternate @ hamiltonian @ alternate.conj().T - negative)))
    completed = block_diag(covariance, unitary @ covariance @ unitary.conj().T)
    alternate_state = block_diag(covariance, alternate @ covariance @ alternate.conj().T)
    if alternate_gap > 1e-8:
        raise RuntimeError("stabilizer image left the negative angular block")
    separation = float(np.linalg.norm(completed - alternate_state))
    if separation < 1e-8:
        raise RuntimeError("sign-pair completion was accidentally unique")
    positive_energy = complex(np.trace(covariance @ hamiltonian))
    completed_energy = positive_energy + complex(np.trace(
        (unitary @ covariance @ unitary.conj().T) @ negative))
    multiplicity = 4 * ANGULAR_LABEL
    return {
        "map": "C plus U C U dagger",
        "U": "sigma2 tensor I_spatial",
        "unique": False,
        "nonuniqueness_witness": "U exp(i H/||H||) conjugates the same Hamiltonian and changes the state",
        "conjugation_gap": conjugation_gap,
        "alternate_conjugation_gap": alternate_gap,
        "sigma1_is_not_this_map_gap": wrong_gap,
        "completed_states_differ_by": separation,
        "completed_energy_over_one_block": float((completed_energy / positive_energy).real),
        "sign_pair_factor": 2,
        "owned_mean_multiplicity": multiplicity,
        "multiplied_again_by_M": False,
        "degeneracy_2_abs_kappa_state_constructed": False,
        "four_dimensional_angular_vertex": None,
        "completed_covariance_shape": [int(completed.shape[0]), int(completed.shape[1])],
        "completed": completed,
        "unitary": unitary,
        "alternate": alternate,
    }


def sea_comparison(hamiltonian, covariance):
    """The negative spectral projector is not the benchmark state."""
    _energies, _vectors, projector = influence.ground_covariance(hamiltonian)
    distance = float(np.linalg.norm(projector - covariance))
    if distance < 1e-8:
        raise RuntimeError("benchmark covariance collapsed onto P_-")
    return {
        "P_minus_frobenius_distance": distance,
        "P_minus_installed_as_physical_vacuum": False,
        "vacuum_branch_evaluated": False,
        "vacuum_branch_owner": "lab/docs/nsc-vacuum-matched-ctp.md",
        "vacuum_branch_definition": "B=Gamma_heat-Gamma_canonical_sea",
    }


def _structures(projected):
    rho, current, radial, angular = sp.symbols("rho current p_radial p_angular", real=True)
    stress = dirac_orthonormal_stress(rho, radial, angular, current)
    square = tensor_square(stress, ETA)
    trace = tensor_trace(stress, ETA)
    evaluate = sp.lambdify(
        (rho, current, radial, angular),
        (square, trace, square - trace**2 / 3, trace**2), "numpy")
    values = evaluate(projected["rho"], projected["current"], projected["p_radial"], projected["p_angular"])
    return {
        "square": np.asarray(values[0], dtype=float),
        "trace": np.asarray(values[1], dtype=float),
        "weyl_structure": np.asarray(values[2], dtype=float),
        "trace_structure": np.asarray(values[3], dtype=float),
    }


def _benchmark_metric():
    metric = smooth_metric(GRID_POINTS, general=True)
    shift = SHIFT_AMPLITUDE * np.cos(2 * np.pi * metric.x / metric.length)
    covariance = fixed_gaussian_covariance(2 * metric.points)
    return metric, shift, covariance


def _summary_array(values):
    array = np.asarray(values, dtype=float)
    return {"max_abs": float(np.max(np.abs(array))), "node0": float(array[0])}


def assemble():
    """Pure finite benchmark. No record is written and no campaign is run."""
    started = time.process_time()
    before = source_bindings()
    with threadpool_limits(limits=1):
        polynomial = contact_polynomial()
        identity = off_shell_identity()
        projection = projection_identity()
        specialization = specialization_witness()
        oracle = one_mode_oracle()
        copies = multiplicity_noise()
        wick = two_mode_wick_oracle()
        quartic = quartic_witness()
        metric, shift, covariance = _benchmark_metric()
        owned, multiplied, analytic, gaps, kernels, antisymmetry, projected = geometric_ledger(
            metric, shift, covariance)
        geometric_wick = wick_geometric(covariance, kernels, analytic)
        structures = _structures(projected)
        algebraic_trace = projected["rho"] - projected["p_radial"] - 2 * projected["p_angular"]
        if np.max(np.abs(algebraic_trace - structures["trace"])) > 1e-10:
            raise RuntimeError("numeric trace left rho - p_r - 2 p_perp")
        measure = (4 * np.pi * metric.lapse * metric.radial_scale
                   * metric.sphere_radius**2 * metric.spacing)
        ward_trace = np.asarray(owned["massless_ward_residual"], dtype=float) / measure
        ward_trace_gap = float(np.max(np.abs(structures["trace"] - ward_trace)))
        if ward_trace_gap > 1e-10:
            raise RuntimeError("raw projected trace left the massless Ward identity")
        hamiltonian = owned["hamiltonian"]
        completion = angular_completion(metric, shift, covariance, hamiltonian)
        sea = sea_comparison(hamiltonian, covariance)
        trace_mean_ratio = multiplied["energy"] / owned["energy"]
        if abs(trace_mean_ratio - multiplied["multiplicity"]) > 1e-12:
            raise RuntimeError("M was not applied exactly once to the mean")
    after = source_bindings()
    if after != before:
        raise RuntimeError("computational source changed during the calculation")
    elapsed = time.process_time() - started
    if elapsed > CPU_LIMIT:
        raise RuntimeError("Dirac contact benchmark exceeded its CPU limit")
    arrays = {
        "covariance": np.ascontiguousarray(covariance),
        "hamiltonian": np.ascontiguousarray(hamiltonian),
        "shift": np.ascontiguousarray(shift),
        "lapse": np.ascontiguousarray(metric.lapse),
        "radial_scale": np.ascontiguousarray(metric.radial_scale),
        "sphere_radius": np.ascontiguousarray(metric.sphere_radius),
        "sign_map": np.ascontiguousarray(completion["unitary"]),
        "alternate_sign_map": np.ascontiguousarray(completion["alternate"]),
        "completed_covariance": np.ascontiguousarray(completion["completed"]),
        "rho": np.ascontiguousarray(projected["rho"]),
        "current": np.ascontiguousarray(projected["current"]),
        "p_radial": np.ascontiguousarray(projected["p_radial"]),
        "p_angular": np.ascontiguousarray(projected["p_angular"]),
        "trace": np.ascontiguousarray(structures["trace"]),
        "weyl_structure": np.ascontiguousarray(structures["weyl_structure"]),
        "trace_structure": np.ascontiguousarray(structures["trace_structure"]),
        "two_mode_covariance": np.ascontiguousarray(wick["covariance"]),
        "two_mode_left": np.ascontiguousarray(wick["left"]),
        "two_mode_right": np.ascontiguousarray(wick["right"]),
    }
    for kind in analytic:
        arrays[f"analytic_{kind}"] = np.ascontiguousarray(analytic[kind])
        arrays[f"finite_difference_gap_{kind}"] = np.ascontiguousarray(gaps[kind])
        arrays[f"kernel_{kind}"] = np.ascontiguousarray(kernels[kind])
    array_hashes = {name: array_digest(values) for name, values in arrays.items()}
    npz_bytes = _npz_bytes(arrays)
    record = {
        "schema": SCHEMA,
        "status": "FINITE_BENCHMARK_CALCULATED_PHYSICAL_DIRAC_CONTACT_OPEN",
        "scope_open": True,
        "blocks_parent_construction": False,
        "missing_relationship": MISSING_RELATIONSHIP,
        "contact_polynomial": polynomial,
        "off_shell_identity": identity,
        "projection_identity": projection,
        "specialization_witness": specialization,
        "one_mode_oracle": oracle,
        "multiplicity_noise": copies,
        "two_mode_wick_oracle": {key: value for key, value in wick.items()
                                 if key not in ("covariance", "left", "right")},
        "quartic_witness": quartic,
        "geometric_variation": {
            "grid_points": GRID_POINTS,
            "angular_label": ANGULAR_LABEL,
            "shift_amplitude": SHIFT_AMPLITUDE,
            "finite_difference_step": FINITE_DIFFERENCE_STEP,
            "owned_function": "block_source",
            "directions": {
                kind: {"max_abs_gap": float(np.max(np.abs(gaps[kind]))),
                       "kernel_antisymmetry": antisymmetry[kind]}
                for kind in analytic},
            "massless_ward_max_abs": float(np.max(np.abs(owned["massless_ward_residual"]))),
            "mass_radius_disagreement_max_abs": float(np.max(np.abs(
                owned["mass_radius_force_disagreement"]))),
            "mean_multiplicity": int(multiplied["multiplicity"]),
            "multiplied_energy_over_one_block": float(trace_mean_ratio),
            "multiplicity_applied_to_noise": False,
        },
        "geometric_wick": geometric_wick,
        "raw_projection": {
            "input": "one-block nodal partial divided by spacing",
            "identified_with_minus_delta_Gamma": False,
            "renormalized": False,
            "composite_subtracted": False,
            "coefficients_chosen": False,
            "extractor_note": "weyl_structure and trace_structure are the two scalars multiplied by -cW/(2 A^2) and -cR/(4 A^2); no physical cW or cR is selected",
            "trace": _summary_array(structures["trace"]),
            "weyl_structure": _summary_array(structures["weyl_structure"]),
            "trace_structure": _summary_array(structures["trace_structure"]),
            "algebraic_trace_gap": float(np.max(np.abs(algebraic_trace - structures["trace"]))),
            "ward_trace_gap": ward_trace_gap,
            "raw_trace_is_the_massless_ward_identity": True,
            "raw_trace_is_not_the_renormalized_dirac_trace": True,
        },
        "angular_completion": {key: value for key, value in completion.items()
                               if key not in ("completed", "unitary", "alternate")},
        "sea_comparison": sea,
        "missing": {
            "four_dimensional_angular_vertex": None,
            "angular_current_map": None,
            "renormalized_contact_expectation": None,
            "vacuum_branch_variation": None,
            "spherical_harmonic_product_state": None,
        },
        "parent_use": (
            "The parent canonical mean force remains M Tr(C H). This benchmark "
            "does not replace it and is not a prerequisite for that evolution."
        ),
        "evolution_ran": False,
        "ultraviolet_campaign_ran": False,
        "new_action_coefficient_selected": False,
        "array_sha256": array_hashes,
        "npz_sha256": hashlib.sha256(npz_bytes).hexdigest(),
        "npz_bytes": len(npz_bytes),
        "source_bindings_before": before,
        "source_bindings_after": after,
        "runtime": {"python": sys.version.split()[0], "numpy": np.__version__, "sympy": sp.__version__},
        "cpu_seconds": elapsed,
        "cpu_limit_seconds": CPU_LIMIT,
    }
    return record, arrays, npz_bytes


def calculate():
    record, _arrays, _npz = assemble()
    return record


def _npz_bytes(arrays):
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **arrays)
    return buffer.getvalue()


def _json_bytes(record):
    return (json.dumps(record, indent=2, allow_nan=False, sort_keys=True) + "\n").encode()


def _exclusive(path, payload):
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY
    descriptor = os.open(path, flags)
    try:
        os.write(descriptor, payload)
    finally:
        os.close(descriptor)


def publish_exclusive(json_path, json_bytes, npz_bytes):
    """Write one JSON record and one NPZ payload, or leave both absent."""
    json_path = Path(json_path)
    if json_path.suffix != ".json":
        raise ValueError("the record path must end with .json")
    npz_path = json_path.with_suffix(".npz")
    total = len(json_bytes) + len(npz_bytes)
    if total > OUTPUT_LIMIT_BYTES:
        raise ValueError("exclusive output exceeds 64 MiB")
    json_path.parent.mkdir(parents=True, exist_ok=True)
    if json_path.exists() or npz_path.exists():
        raise FileExistsError("Dirac contact record exists; use --check or a named successor")
    _exclusive(json_path, json_bytes)
    try:
        _exclusive(npz_path, npz_bytes)
    except Exception:
        json_path.unlink(missing_ok=True)
        raise
    return {"json_bytes": len(json_bytes), "npz_bytes": len(npz_bytes), "total_bytes": total}


def write_record(output, commit):
    """Exclusive frozen record. Root calls this after the source checkpoint."""
    if not isinstance(commit, str) or len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit):
        raise ValueError("an explicit full lowercase frozen Git commit is required")
    path = Path(output).resolve()
    npz_path = path.with_suffix(".npz")
    if path.suffix != ".json":
        raise ValueError("the record path must end with .json")
    if path.exists() or npz_path.exists():
        raise FileExistsError("Dirac contact record exists; use --check or a named successor")
    record, _arrays, npz_bytes = assemble()
    record["frozen_git"] = authenticate_git(commit, record["source_bindings_before"])
    if source_bindings() != record["source_bindings_before"]:
        raise RuntimeError("source changed before the exclusive record write")
    publish_exclusive(path, _json_bytes(record), npz_bytes)
    return record


def same_calculation(saved, fresh):
    if isinstance(fresh, dict):
        if set(saved) != set(fresh):
            raise ValueError("saved Dirac-contact fields changed")
        for name in fresh:
            if name != "cpu_seconds":
                same_calculation(saved[name], fresh[name])
        return
    if isinstance(fresh, list):
        if len(saved) != len(fresh):
            raise ValueError("saved Dirac-contact shape changed")
        for index, item in enumerate(fresh):
            same_calculation(saved[index], item)
        return
    if isinstance(fresh, float):
        if not np.isfinite(saved) or not np.isclose(saved, fresh, rtol=1e-8, atol=1e-12):
            raise ValueError("saved Dirac-contact number changed")
        return
    if saved != fresh:
        raise ValueError("saved Dirac-contact scope changed")


def check_record(output):
    """Replay the benchmark and authenticate the frozen bytes. Nothing is evolved."""
    path = Path(output)
    record = json.loads(path.read_text())
    npz_path = path.with_suffix(".npz")
    if record.get("schema") != SCHEMA:
        raise ValueError("Dirac-contact record has the wrong schema")
    payload = npz_path.read_bytes()
    if len(path.read_bytes()) + len(payload) > OUTPUT_LIMIT_BYTES:
        raise ValueError("Dirac-contact record exceeds 64 MiB")
    if hashlib.sha256(payload).hexdigest() != record["npz_sha256"]:
        raise ValueError("Dirac-contact NPZ bytes changed")
    with np.load(io.BytesIO(payload)) as loaded:
        for name, digest_value in record["array_sha256"].items():
            if array_digest(loaded[name]) != digest_value:
                raise ValueError("Dirac-contact array changed: " + name)
    bindings = source_bindings()
    if bindings != record["source_bindings_before"] or bindings != record["source_bindings_after"]:
        raise ValueError("Dirac-contact producer bytes changed")
    frozen = authenticate_git(record["frozen_git"]["commit"], bindings)
    if frozen != record["frozen_git"]:
        raise ValueError("frozen Git object closure changed")
    fresh, _arrays, _npz = assemble()
    same_calculation({key: value for key, value in record.items() if key != "frozen_git"}, fresh)
    return record
