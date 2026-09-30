"""Direct-order conformal ADM source on the existing Fourier grid.

The continuum spherical Hamiltonian in the canonical half-density frame
u = r sqrt(q) psi is

    H = σ2 {N/q, P}/2 + σ1 N κ/r - I {β, P}/2,    P = -i ∂_x.

Physical fields are N = r L and q = r Q, so the same operator is

    H = σ2 ⊗ {diag(L/Q), P}/2 + σ1 ⊗ diag(κ L) - I ⊗ {diag β, P}/2.

This is that ordering on the existing CovariantStaticMetric momentum
matrix. It is not the historical finite-grid product rule
sqrt(N) {P, 1/q} sqrt(N) retained by nsc_adm_source. That rule is a
different finite ordering. It is not marked false here.

One positive angular block is multiplied by M = 4κ once. That factor
is the existing isotropic copy count; a separate -κ block is not added.
Nodal values are partial derivatives of M Tr(C H), not densities: divide
by the spatial spacing before reading them as densities.

The canonical sea is not subtracted. An absolute renormalized stress
also needs the vacuum branch B = Gamma_heat - Gamma_canonical_sea,
owned by lab/docs/nsc-vacuum-matched-ctp.md, counted once. This module
is not that branch, not a 6-mode embedding, and not a closed source.
Geometry passed to the work ledger is a prescribed control, not a
solution of the source equations.
"""
from __future__ import annotations

import numpy as np
from dataclasses import replace
from scipy.linalg import expm

from .nsc_adm_source import adm_hamiltonian, canonical_force_gradients
from .nsc_covariant_operator import (
    CovariantStaticMetric,
    SIGMA1,
    SIGMA2,
)
from .nsc_influence import _covariance


VACUUM_BRANCH_OWNER = "lab/docs/nsc-vacuum-matched-ctp.md"
VACUUM_BRANCH_DEFINITION = "B=Gamma_heat-Gamma_canonical_sea"
FRAME_MAPPING = "missing"
PRODUCT_RULE_NAME = "finite-grid product rule"
HALF_DENSITY = "u=r sqrt(q) psi"
SOURCE_INTERFACES = ("covariance", "relative_difference")

# Raw Tr((I/2) H) vanishes because this H is traceless on the Fourier grid.
# That identity is not a renormalized stress.
SPIN_IDENTITY_NOTE = (
    "naive one-block spin identity C=I/2 has raw source zero; "
    "this does not mean the renormalized source is zero"
)


def _limits():
    return {
        "absolute_renormalized_stress": False,
        "full_source_closure": False,
        "frame_mapping": FRAME_MAPPING,
        "six_mode_embedding": False,
        "autonomous_metric_evolution": False,
        "vacuum_branch_accounted_in_this_module": False,
        "vacuum_branch_owner": VACUUM_BRANCH_OWNER,
        "vacuum_branch_definition": VACUUM_BRANCH_DEFINITION,
        "sea_accounting": "once, in the vacuum branch, not again in this module",
        "historical_product_rule_marked_false": False,
        "historical_ordering": PRODUCT_RULE_NAME,
        "naive_spin_identity_means_zero_renormalized_source": False,
        "half_density": HALF_DENSITY,
        "geometry_role": "prescribed_control",
    }


def _angular_label(kappa, *, positive):
    if isinstance(kappa, bool) or not isinstance(kappa, (int, np.integer)):
        raise ValueError("integer angular label required")
    kappa = int(kappa)
    if positive:
        if kappa <= 0:
            raise ValueError(
                "positive angular sector required; M=4κ already counts "
                "identical isotropic copies, so -κ must not be added"
            )
        return kappa
    if kappa == 0:
        return 0
    return kappa


def sector_multiplicity(kappa):
    """M=4κ once for one positive sector and its identical isotropic copies.

    The negative angular sign is inside this factor. Do not also trace a
    separate -κ block.
    """
    return 4 * _angular_label(kappa, positive=True)


def conformal_variables(metric):
    """Return L, Q for physical N = r L and q = r Q."""
    if not isinstance(metric, CovariantStaticMetric):
        raise ValueError("existing covariant static metric required")
    radius = metric.sphere_radius
    return metric.lapse / radius, metric.radial_scale / radius


def _shift_vector(metric, shift):
    beta = np.asarray(shift, dtype=float)
    if np.iscomplexobj(shift) and np.max(np.abs(np.imag(np.asarray(shift)))) > 0:
        raise ValueError("physical shift must be real")
    if beta.shape != (metric.points,) or not np.isfinite(beta).all():
        raise ValueError("finite real shift on the existing spatial grid required")
    return beta


def _anticommutator(values, momentum):
    """{diag(values), P}/2 with the existing row/column product."""
    values = np.asarray(values, dtype=float)
    return 0.5 * (values[:, None] * momentum + momentum * values[None, :])


def direct_hamiltonian(metric, shift, kappa=1):
    """One angular block in the direct continuum ordering.

    kappa=0 is the massless angular control, not a physical 4κ sector.
    Negative kappa is accepted only so a caller can see that the coupling
    API does not add it.
    """
    if not isinstance(metric, CovariantStaticMetric):
        raise ValueError("existing covariant static metric required")
    label = _angular_label(kappa, positive=False)
    length_density, radial_density = conformal_variables(metric)
    beta = _shift_vector(metric, shift)
    if np.min(radial_density) <= 0:
        raise ValueError("positive conformal radial density Q required")
    momentum = metric.momentum_matrix
    kinetic = _anticommutator(length_density / radial_density, momentum)
    mass = np.diag(label * length_density)
    shift_term = _anticommutator(beta, momentum)
    return (
        np.kron(SIGMA2, kinetic)
        + np.kron(SIGMA1, mass)
        - np.kron(np.eye(2), shift_term)
    )


def _spin_trace(matrix, gamma, points):
    blocks = np.asarray(matrix, dtype=complex).reshape(2, points, 2, points)
    return np.einsum("ba,aibj->ij", gamma, blocks)


def _prepare_state(state, points, *, state_interface, relative_name):
    if state_interface not in SOURCE_INTERFACES:
        raise ValueError(
            "state_interface must be 'covariance' or 'relative_difference'"
        )
    if state_interface == "covariance":
        if relative_name is not None:
            raise ValueError("covariance interface does not take a relative name")
        matrix = _covariance(state)
        if matrix.shape != (2 * points, 2 * points):
            raise ValueError("covariance and spatial operator dimensions differ")
        return matrix, {"state_interface": "covariance", "relative_name": None,
                        "vacuum_subtracted": False}
    if not isinstance(relative_name, str) or not relative_name.strip():
        raise ValueError(
            "relative_difference requires an explicit name; "
            "the canonical sea is not subtracted silently"
        )
    matrix = np.asarray(state, dtype=complex)
    if matrix.shape != (2 * points, 2 * points) or not np.isfinite(matrix).all():
        raise ValueError("finite relative-difference matrix on this block required")
    if not np.allclose(matrix, matrix.conj().T, rtol=0, atol=1e-11):
        raise ValueError("relative difference must be Hermitian")
    return matrix, {
        "state_interface": "relative_difference",
        "relative_name": relative_name,
        "vacuum_subtracted": False,
    }


def physical_forces_from_conformal(force_l, force_q, force_beta, length_density,
                                   radial_density, radius):
    """Chain rule at fixed covariance, node by node.

    F_N = F_L / r, F_q = F_Q / r,
    F_r = -(L F_L + Q F_Q) / r, and F_β is unchanged.
    """
    arrays = [
        np.asarray(force_l, dtype=float),
        np.asarray(force_q, dtype=float),
        np.asarray(force_beta, dtype=float),
        np.asarray(length_density, dtype=float),
        np.asarray(radial_density, dtype=float),
        np.asarray(radius, dtype=float),
    ]
    if any(array.shape != arrays[0].shape or not np.isfinite(array).all()
           for array in arrays):
        raise ValueError("conformal chain rule needs one finite vector per field")
    if np.min(arrays[5]) <= 0:
        raise ValueError("positive radius required")
    force_l, force_q, force_beta, length_density, radial_density, radius = arrays
    return {
        "N": force_l / radius,
        "q": force_q / radius,
        "r": -(length_density * force_l + radial_density * force_q) / radius,
        "beta": force_beta.copy(),
    }


def _nodal_real(vector):
    return np.asarray(vector, dtype=float).copy()


def block_source(metric, shift, state, kappa=1, *, state_interface="covariance",
                 relative_name=None):
    """Nodal derivatives of one-block Tr(C H). Multiplicity is not applied.

    kappa=0 is allowed and is the massless angular ordering control.
    """
    if not isinstance(metric, CovariantStaticMetric):
        raise ValueError("existing covariant static metric required")
    label = _angular_label(kappa, positive=False)
    matrix, interface = _prepare_state(
        state, metric.points, state_interface=state_interface,
        relative_name=relative_name,
    )
    hamiltonian = direct_hamiltonian(metric, shift, label)
    if matrix.shape != hamiltonian.shape:
        raise ValueError("state and Hamiltonian dimensions differ")
    length_density, radial_density = conformal_variables(metric)
    radius = metric.sphere_radius
    momentum = metric.momentum_matrix
    trace0 = _spin_trace(matrix, np.eye(2), metric.points)
    trace1 = _spin_trace(matrix, SIGMA1, metric.points)
    trace2 = _spin_trace(matrix, SIGMA2, metric.points)
    force_c = 0.5 * np.diag(momentum @ trace2 + trace2 @ momentum).real
    force_l = force_c / radial_density + label * np.diag(trace1).real
    force_q = -force_c * length_density / radial_density**2
    force_beta = -0.5 * np.diag(momentum @ trace0 + trace0 @ momentum).real
    physical = physical_forces_from_conformal(
        force_l, force_q, force_beta, length_density, radial_density, radius,
    )
    # The mass term alone is -∂_r of N κ/r at fixed N. It must agree with
    # the chain-rule radius force; a kinetic ordering error would split them.
    mass_radius = -label * metric.lapse * np.diag(trace1).real / radius**2
    ward = (
        metric.lapse * physical["N"]
        + metric.radial_scale * physical["q"]
        + radius * physical["r"]
    )
    energy = float(np.trace(matrix @ hamiltonian).real)
    return {
        "kappa": label,
        "multiplicity_applied": False,
        "hamiltonian": hamiltonian,
        "energy": energy,
        "nodal": {
            "L": _nodal_real(force_l),
            "Q": _nodal_real(force_q),
            "beta": _nodal_real(force_beta),
            "N": _nodal_real(physical["N"]),
            "q": _nodal_real(physical["q"]),
            "r": _nodal_real(physical["r"]),
        },
        "massless_ward_residual": _nodal_real(ward),
        "mass_radius_force_disagreement": _nodal_real(physical["r"] - mass_radius),
        "state": interface,
        "vacuum_subtracted": False,
        "limits": _limits(),
        "spin_identity_note": SPIN_IDENTITY_NOTE,
    }


def conformal_source(metric, shift, state, kappa=1, *,
                     state_interface="covariance", relative_name=None):
    """Coupling API: M Tr(C H) and its nodal derivatives, M=4κ once.

    state_interface is 'covariance' or an explicitly named
    'relative_difference'. Nothing is subtracted unless that difference
    is the matrix the caller passed. Divide nodal derivatives by the
    spatial spacing before treating them as densities.
    """
    label = _angular_label(kappa, positive=True)
    block = block_source(
        metric, shift, state, label, state_interface=state_interface,
        relative_name=relative_name,
    )
    multiplicity = sector_multiplicity(label)
    nodal = {name: multiplicity * values for name, values in block["nodal"].items()}
    ward = multiplicity * block["massless_ward_residual"]
    return {
        "schema": "NSC-CONFORMAL-ADM-SOURCE-v1",
        "kappa": label,
        "multiplicity": multiplicity,
        "multiplicity_rule": (
            "M=4κ once per positive sector and identical isotropic copies; "
            "no extra ± doubling"
        ),
        "hamiltonian": block["hamiltonian"],
        "hamiltonian_is_one_block": True,
        "energy": multiplicity * block["energy"],
        "sector_energy": block["energy"],
        "spacing": float(metric.spacing),
        "nodal_derivative_convention": "partial derivative of M Tr(C H); divide by spacing for a density",
        "variables": "N=r L, q=r Q",
        "ordering": "direct continuum anticommutator",
        "nodal": nodal,
        "massless_ward_residual": ward,
        "mass_radius_force_disagreement": multiplicity * block["mass_radius_force_disagreement"],
        "state_interface": block["state"]["state_interface"],
        "relative_name": block["state"]["relative_name"],
        "vacuum_subtracted": False,
        "source_interfaces": list(SOURCE_INTERFACES),
        "limits": block["limits"],
        "spin_identity_note": SPIN_IDENTITY_NOTE,
    }


def fixed_gaussian_covariance(dimension):
    """Dense Gaussian covariance from one fixed unitary, not a recorded 6J state.

    Occupations lie in (0, 1). The unitary is exp(i G/||G||) for an
    explicit Hermitian generator, so the draw does not use a random stream.
    """
    if isinstance(dimension, bool) or not isinstance(dimension, (int, np.integer)):
        raise ValueError("covariance dimension must be an integer")
    dimension = int(dimension)
    if dimension < 2:
        raise ValueError("at least a two-dimensional covariance is required")
    index_i, index_j = np.indices((dimension, dimension))
    generator = (
        np.sin((index_i + 1) * (index_j + 2) + 0.3 * (index_i - index_j))
        + 1j * np.cos(0.5 * (index_i + 3) * (index_j + 1))
    )
    generator = 0.5 * (generator + generator.conj().T)
    unitary = expm(1j * generator / np.linalg.norm(generator, ord=2))
    occupation = 0.2 + 0.55 * np.sin((np.arange(dimension) + 1) * 0.7) ** 2
    covariance = (unitary * occupation) @ unitary.conj().T
    covariance = 0.5 * (covariance + covariance.conj().T)
    return _covariance(covariance)


def spin_identity_covariance(points):
    """Naive one-block C=I/2. Not the recorded 6-mode state and not an embedding."""
    if isinstance(points, bool) or not isinstance(points, (int, np.integer)) or int(points) < 1:
        raise ValueError("positive point count required")
    return 0.5 * np.eye(2 * int(points), dtype=complex)


def _radius_factor(metric, amplitude):
    angle = 2 * np.pi * metric.x / metric.length
    return np.exp(amplitude * np.cos(2 * angle))


def _fixed_conformal_radius(metric, factor):
    factor = np.asarray(factor, dtype=float)
    if factor.shape != (metric.points,) or not np.isfinite(factor).all() or np.min(factor) <= 0:
        raise ValueError("positive finite radius factor on the grid required")
    length_density, radial_density = conformal_variables(metric)
    radius = metric.sphere_radius * factor
    return replace(
        metric,
        lapse=length_density * radius,
        radial_scale=radial_density * radius,
        sphere_radius=radius,
    )


def ordering_comparison(metric, shift, covariance, kappa=1, *, amplitude=0.12):
    """Same physical covariance, direct ordering versus the historical product rule.

    The finite-(L, Q) radius change is a control probe. It is not a metric
    update and it does not mark historical product-rule results false.
    Comparison traces are one block, before M=4κ.
    """
    label = _angular_label(kappa, positive=True)
    if not isinstance(metric, CovariantStaticMetric):
        raise ValueError("existing covariant static metric required")
    beta = _shift_vector(metric, shift)
    state = _covariance(covariance)
    if state.shape != (2 * metric.points, 2 * metric.points):
        raise ValueError("covariance and spatial operator dimensions differ")
    factor = _radius_factor(metric, amplitude)
    moved = _fixed_conformal_radius(metric, factor)
    direct_here = direct_hamiltonian(metric, beta, label)
    direct_there = direct_hamiltonian(moved, beta, label)
    product_here = adm_hamiltonian(metric, beta, label)
    product_there = adm_hamiltonian(moved, beta, label)
    direct_block = block_source(metric, beta, state, label)
    historical = canonical_force_gradients(metric, beta, state, label)
    historical_ward = (
        metric.lapse * historical["N"]
        + metric.radial_scale * historical["q"]
        + metric.sphere_radius * historical["r"]
    )
    def _trace_change(hamiltonian_there, hamiltonian_here):
        return float(np.trace(state @ (hamiltonian_there - hamiltonian_here)).real)

    return {
        "ordering_name": PRODUCT_RULE_NAME,
        "historical_results_marked_false": False,
        "same_physical_covariance": True,
        "multiplicity_omitted": True,
        "kappa": label,
        "radius_amplitude": float(amplitude),
        "frame_mapping": FRAME_MAPPING,
        "source_interfaces": list(SOURCE_INTERFACES),
        "direct_fixed_LQ_trace_change": _trace_change(direct_there, direct_here),
        "product_rule_fixed_LQ_trace_change": _trace_change(product_there, product_here),
        "direct_fixed_LQ_hamiltonian_max": float(np.max(np.abs(direct_there - direct_here))),
        "product_rule_fixed_LQ_hamiltonian_frobenius": float(
            np.linalg.norm(product_there - product_here)
        ),
        "direct_ward_l2": float(np.linalg.norm(direct_block["massless_ward_residual"])),
        "product_rule_ward_l2": float(np.linalg.norm(historical_ward)),
        "physical_radius_force_l2": float(np.linalg.norm(direct_block["nodal"]["r"])),
        "historical_radius_force_l2": float(np.linalg.norm(historical["r"])),
        "radius_force_disagreement_l2": float(
            np.linalg.norm(direct_block["nodal"]["r"] - historical["r"])
        ),
        "lapse_force_disagreement_l2": float(
            np.linalg.norm(direct_block["nodal"]["N"] - historical["N"])
        ),
        "shift_force_disagreement_l2": float(
            np.linalg.norm(direct_block["nodal"]["beta"] - historical["beta"])
        ),
        "limits": _limits(),
    }


def _prescribed_knot(metric, shift, time, amplitude):
    """Explicit geometry at one time. The increment is not the source gradient."""
    angle = 2 * np.pi * metric.x / metric.length
    lapse = metric.lapse * np.exp(amplitude * np.sin(time) * np.cos(angle))
    radial = metric.radial_scale * np.exp(
        amplitude * 0.7 * np.sin(time) * np.sin(2 * angle)
    )
    # sin(time) so the supplied metric is recovered at time 0, while r still moves.
    radius = metric.sphere_radius * np.exp(
        amplitude * 0.5 * np.sin(time) * np.cos(2 * angle)
    )
    beta = np.asarray(shift, dtype=float) + amplitude * 0.4 * np.sin(time) * np.cos(angle)
    return replace(metric, lapse=lapse, radial_scale=radial, sphere_radius=radius), beta


def prescribed_work_ledger(metric, shift, covariance, kappa=1, *, times=None,
                           amplitude=0.04, linear_step=1e-6):
    """Unitary evolution of C on a prescribed geometry.

    The canonical identity is d/dt Tr(C H) = Tr(C dH/dt) when
    i dC/dt = [H, C]. Geometry is an input. This does not solve F_A = 0
    and does not include the vacuum branch.
    """
    label = _angular_label(kappa, positive=True)
    if not isinstance(metric, CovariantStaticMetric):
        raise ValueError("existing covariant static metric required")
    beta0 = _shift_vector(metric, shift)
    state = _covariance(covariance)
    if times is None:
        times = (0.0, 0.17, 0.41, 0.73)
    times = tuple(float(item) for item in times)
    if len(times) < 2 or any(not np.isfinite(item) for item in times):
        raise ValueError("at least two finite prescribed times required")
    if np.any(np.diff(times) <= 0):
        raise ValueError("prescribed times must increase")
    if not np.isfinite(amplitude) or abs(amplitude) >= 0.5:
        raise ValueError("prescribed amplitude must keep the metric positive")
    if not np.isfinite(linear_step) or linear_step <= 0 or linear_step > 1e-3:
        raise ValueError("linearization step must be a small positive number")
    multiplicity = sector_multiplicity(label)
    knots = [_prescribed_knot(metric, beta0, time, amplitude) for time in times]
    covariance_now = state.copy()
    unitary_residuals = []
    work_residuals = []
    works = []
    energies = []
    for index, (knot_metric, knot_shift) in enumerate(knots):
        hamiltonian = direct_hamiltonian(knot_metric, knot_shift, label)
        energy = multiplicity * float(np.trace(covariance_now @ hamiltonian).real)
        energies.append(energy)
        if index + 1 == len(knots):
            break
        step = times[index + 1] - times[index]
        unitary = expm(-1j * step * hamiltonian)
        covariance_now = unitary @ covariance_now @ unitary.conj().T
        covariance_now = 0.5 * (covariance_now + covariance_now.conj().T)
        evolved = multiplicity * float(np.trace(covariance_now @ hamiltonian).real)
        unitary_residuals.append(evolved - energy)
        next_metric, next_shift = knots[index + 1]
        next_hamiltonian = direct_hamiltonian(next_metric, next_shift, label)
        work = multiplicity * float(
            np.trace(covariance_now @ (next_hamiltonian - hamiltonian)).real
        )
        after = multiplicity * float(np.trace(covariance_now @ next_hamiltonian).real)
        works.append(work)
        work_residuals.append(after - evolved - work)
    values = np.linalg.eigvalsh(covariance_now)
    start_metric, start_shift = knots[0]
    start_force = conformal_source(start_metric, start_shift, state, label)
    bumped, bumped_shift = _prescribed_knot(metric, beta0, times[0] + linear_step, amplitude)
    base_h = direct_hamiltonian(start_metric, start_shift, label)
    bumped_h = direct_hamiltonian(bumped, bumped_shift, label)
    exact_linear = multiplicity * float(np.trace(state @ (bumped_h - base_h)).real)
    delta = {
        "N": bumped.lapse - start_metric.lapse,
        "q": bumped.radial_scale - start_metric.radial_scale,
        "r": bumped.sphere_radius - start_metric.sphere_radius,
        "beta": bumped_shift - start_shift,
    }
    linear_force = sum(float(np.dot(start_force["nodal"][name], delta[name]))
                       for name in ("N", "q", "r", "beta"))
    probe = np.concatenate([delta[name] for name in ("N", "q", "r", "beta")])
    gradient = np.concatenate([start_force["nodal"][name] for name in ("N", "q", "r", "beta")])
    denom = np.linalg.norm(probe) * np.linalg.norm(gradient)
    alignment = float(np.dot(probe, gradient) / denom) if denom else 0.0
    return {
        "schema": "NSC-CONFORMAL-ADM-WORK-LEDGER-v1",
        "geometry_role": "prescribed_control",
        "autonomous_metric_evolution": False,
        "vacuum_branch_included": False,
        "multiplicity": multiplicity,
        "times": times,
        "energies": energies,
        "works": works,
        "total_work": float(sum(works)),
        "energy_change": float(energies[-1] - energies[0]),
        "unitary_residual_max": float(np.max(np.abs(unitary_residuals))),
        "work_residual_max": float(np.max(np.abs(work_residuals))),
        "ledger_closure": float(energies[-1] - energies[0] - sum(works)),
        "linear_force_work": linear_force,
        "linear_trace_work": exact_linear,
        "linear_residual": float(linear_force - exact_linear),
        "force_alignment_of_prescribed_step": alignment,
        "final_covariance_min": float(np.min(values)),
        "final_covariance_max": float(np.max(values)),
        "limits": _limits(),
    }
