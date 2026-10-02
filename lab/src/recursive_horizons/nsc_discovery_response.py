"""Stage-5 first wave: full-state directional primitive for the conformal rates.

The primitive is the analytic Jacobian-vector product of the existing
conformal Galerkin rate and of its NestedPair image. Geometry, canonical
momenta, both spinor columns and the source occupations are differentiated
together. A field-only or frozen-state derivative is not this map.

The linear operators are the quadrature derivative and the antiperiodic
momentum already stored on the fine system. No new Fourier matrix is built.
Multiplicity stays the existing ``4κ`` factor inside the nodal forces, once.
In the conformal chart both ``force_L`` and ``force_Q`` enter ``p_Q``.
Nested momenta use the owned map ``π = Δx_g Wᵀ p``.

Finite differences are not an analytic substitute. The CTP future effective
stress is specified and not returned as a number: the retarded change of
state is a missing primitive, and the induced spectral force is not
recomputed here.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .nsc_causal_common import SpectralInducedSource
from .nsc_influence import _covariance, connected
from .nsc_spherical_feedback_action import alpha_of, feedback_F, feedback_V, feedback_Z, partial_F
from . import nsc_nested_parent_child as nested
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin


SCHEMA = "NSC-DISCOVERY-RESPONSE-v1"
WAVE = "stage5-first-wave"
METHOD = "analytic_jacobian_vector"
FINITE_DIFFERENCE_USED_AS_ANALYTIC = False
NEWTON_PRODUCTION = False
PRIMARY_READOUT = "child_regional_content"
SECONDARY_READOUT = "child_proper_mean_r"
KERNEL_2X2_IS_SPATIAL_OUTSIDE = False
INDUCED_OWNER = f"{SpectralInducedSource.__module__}.{SpectralInducedSource.__name__}"
JACOBIAN_OWNER = "nsc_spherical_galerkin_coupling._projected_radius_jacobian"
CHILD_CLOCK = 2.0

IMPLEMENTED = (
    "analytic_rate_jacobian_vector",
    "occupation_tangent",
    "force_L_and_force_Q",
    "multiplicity_4kappa_once",
    "canonical_W_pi_over_dx",
    "linearized_radius_constraint_tangent",
    "proper_clock_tangent",
    "matched_tau_correction",
    "full_gaussian_covariance",
    "connected_bilinear",
    "child_regional_content",
    "child_proper_mean_r",
    "initial_modal_cross",
)
MISSING_PRIMITIVES = (
    "ctp_future_effective_stress",
    "retarded_delta_state",
)
BLOCKERS = MISSING_PRIMITIVES + (
    "no_evolution_campaign",
    "induced_coefficient_not_rederived",
    "kernel_2x2_is_not_spatial_outside",
)
CHECKS = (
    "centred two-sided difference of selected full-state directions",
    "small valid manufactured conformal state",
    "one saved initial state, conformal chart, no evolution",
    "connected bilinear on the full column Gaussian",
    "owned projected radius Jacobian, one linear solve",
    "matched proper-time correction on a manufactured clock",
)


def specification():
    """Frozen first-wave contract. Missing primitives stay named."""
    return {
        "schema": SCHEMA,
        "wave": WAVE,
        "method": METHOD,
        "finite_difference_used_as_analytic": FINITE_DIFFERENCE_USED_AS_ANALYTIC,
        "newton_production": NEWTON_PRODUCTION,
        "primary_readout": PRIMARY_READOUT,
        "secondary_readout": SECONDARY_READOUT,
        "kernel_2x2_is_spatial_outside": KERNEL_2X2_IS_SPATIAL_OUTSIDE,
        "induced_owner": INDUCED_OWNER,
        "jacobian_owner": JACOBIAN_OWNER,
        "canonical_momentum": "pi = dx_g W.T p_nodal",
        "multiplicity": "4 kappa once, inside the existing nodal forces",
        "implemented": list(IMPLEMENTED),
        "missing_primitives": list(MISSING_PRIMITIVES),
        "blockers": list(BLOCKERS),
        "checks": list(CHECKS),
        "future_effective_stress": (
            "complementary_force + initial_cross + retarded_delta_state; "
            "sum not formed while retarded_delta_state is missing"
        ),
    }


@dataclass
class StateTangent:
    """Directional increment of one full state, not a second Cauchy state.

    Momenta are nodal when passed to ``rate_jacobian_vector`` and canonical
    ``π`` when passed to ``nested_rate_jacobian_vector``. ``occupations`` is
    ``δc``, not the occupation vector itself.
    """

    Q: np.ndarray
    r: np.ndarray
    chi: np.ndarray
    p_Q: np.ndarray
    p_r: np.ndarray
    p_chi: np.ndarray
    phi0: np.ndarray
    phi1: np.ndarray
    occupations: np.ndarray


@dataclass
class RateTangent:
    """Analytic directional image of one rate. Forces stay on the fine nodes."""

    Q: np.ndarray
    r: np.ndarray
    chi: np.ndarray
    p_Q: np.ndarray
    p_r: np.ndarray
    p_chi: np.ndarray
    phi0: np.ndarray
    phi1: np.ndarray
    force_L: np.ndarray
    force_Q: np.ndarray
    force_beta: np.ndarray
    fieldwork_power: float
    method: str = METHOD


@dataclass(frozen=True)
class PreparedRadiusTangent:
    """One linear solve for ``δr`` on the owned projected radius constraint."""

    available: bool
    delta_r: np.ndarray | None
    missing_primitive: str | None
    linear_residual_max: float | None
    jacobian_owner: str = JACOBIAN_OWNER
    newton_used: bool = False


def zero_tangent(grid):
    """Zero full-state tangent shaped to one Galerkin grid."""
    count = int(grid.fine.occupations.size)
    return StateTangent(
        np.zeros(grid.ng), np.zeros(grid.ng), np.zeros(grid.ng),
        np.zeros(grid.ng), np.zeros(grid.ng), np.zeros(grid.ng),
        np.zeros((grid.nf, count), dtype=complex),
        np.zeros((grid.nf, count), dtype=complex),
        np.zeros(count),
    )


def gaussian_covariance(phi0, phi1, occupations):
    """Full one-particle Gaussian ``C = Φ diag(c) Φ†`` on the column space.

    Complement eigenvalues are the unoccupied directions of this finite
    matrix. Occupation 0 there is the support of ``C``, not a vacuum claim.
    """
    weights = np.asarray(occupations, dtype=float)
    phi0 = np.asarray(phi0, dtype=complex)
    phi1 = np.asarray(phi1, dtype=complex)
    if phi0.shape != phi1.shape or phi0.ndim != 2 or weights.shape != (phi0.shape[1],):
        raise ValueError("covariance needs matching columns and one weight each")
    if not np.isfinite(weights).all() or np.any(weights < 0) or np.any(weights > 1):
        raise ValueError("Gaussian weights must lie in [0, 1]")
    columns = np.vstack((phi0, phi1))
    if not np.isfinite(columns).all():
        raise ValueError("nonfinite columns")
    covariance = (columns * weights) @ columns.conj().T
    checked = _covariance(covariance)
    return checked


def connected_bilinear(covariance, left, right):
    """``Tr(C A (I-C) B)`` from the existing influence owner. Not an induced term."""
    return connected(covariance, left, right)


def matched_tau_correction(delta_observable_t, observable_dot, delta_tau, tau_dot):
    """``δO|τ = δO|t − Ȯ δτ / τ̇`` at one nonzero clock rate."""
    delta_observable_t = float(delta_observable_t)
    observable_dot = float(observable_dot)
    delta_tau = float(delta_tau)
    tau_dot = float(tau_dot)
    if not all(np.isfinite(value) for value in (delta_observable_t, observable_dot, delta_tau, tau_dot)):
        raise ValueError("matched proper time needs finite inputs")
    if tau_dot == 0.0:
        raise ValueError("matched proper time needs a nonzero clock rate")
    return delta_observable_t - observable_dot * delta_tau / tau_dot


def _real_vector(values, size, name):
    array = np.asarray(values, dtype=float)
    if np.iscomplexobj(values):
        raise ValueError(f"{name} tangent must be real")
    if array.shape != (size,) or not np.isfinite(array).all():
        raise ValueError(f"{name} tangent must be a finite vector of length {size}")
    return array


def _columns(values, shape, name):
    array = np.asarray(values, dtype=complex)
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError(f"{name} tangent must be a finite array of shape {shape}")
    return array


def _validate_galerkin(grid, state, tangent):
    if grid.gauge not in ("conformal", "prescribed"):
        raise ValueError("rate tangent requires the prescribed or conformal gauge")
    galerkin._check_state_dimensions(grid, state)
    count = int(grid.fine.occupations.size)
    fields = {
        "Q": _real_vector(tangent.Q, grid.ng, "Q"),
        "r": _real_vector(tangent.r, grid.ng, "r"),
        "chi": _real_vector(tangent.chi, grid.ng, "chi"),
        "p_Q": _real_vector(tangent.p_Q, grid.ng, "p_Q"),
        "p_r": _real_vector(tangent.p_r, grid.ng, "p_r"),
        "p_chi": _real_vector(tangent.p_chi, grid.ng, "p_chi"),
        "phi0": _columns(tangent.phi0, (grid.nf, count), "phi0"),
        "phi1": _columns(tangent.phi1, (grid.nf, count), "phi1"),
        "occupations": _real_vector(tangent.occupations, count, "occupations"),
    }
    return StateTangent(**fields)


def _prolong_tangent(grid, tangent):
    return StateTangent(
        Q=galerkin.prolong_geometry(grid, tangent.Q),
        r=galerkin.prolong_geometry(grid, tangent.r),
        chi=galerkin.prolong_geometry(grid, tangent.chi),
        p_Q=galerkin.prolong_geometry(grid, tangent.p_Q),
        p_r=galerkin.prolong_geometry(grid, tangent.p_r),
        p_chi=galerkin.prolong_geometry(grid, tangent.p_chi),
        phi0=galerkin.prolong_columns(grid, tangent.phi0),
        phi1=galerkin.prolong_columns(grid, tangent.phi1),
        occupations=np.array(tangent.occupations, dtype=float, copy=True),
    )


def _anticommutator(momentum, values, component):
    return 0.5 * (
        values[:, None] * (momentum @ component)
        + momentum @ (values[:, None] * component)
    )


def _images_tangent(momentum, phi0, phi1, dphi0, dphi1, length, dlength, radial, dradial, shift, dshift, kappa):
    """Directional image of ``apply_dirac``. ``a = L/Q`` uses the same product."""
    if np.min(radial) <= 0:
        raise coupling.PositiveChartExit("Q_left_positive_chart", np.nan, None)
    conformal_factor = length / radial
    delta_factor = dlength / radial - length * dradial / radial ** 2
    kinetic0 = _anticommutator(momentum, conformal_factor, phi0)
    kinetic1 = _anticommutator(momentum, conformal_factor, phi1)
    dkinetic0 = _anticommutator(momentum, delta_factor, phi0) + _anticommutator(momentum, conformal_factor, dphi0)
    dkinetic1 = _anticommutator(momentum, delta_factor, phi1) + _anticommutator(momentum, conformal_factor, dphi1)
    shift0 = _anticommutator(momentum, shift, phi0)
    shift1 = _anticommutator(momentum, shift, phi1)
    dshift0 = _anticommutator(momentum, dshift, phi0) + _anticommutator(momentum, shift, dphi0)
    dshift1 = _anticommutator(momentum, dshift, phi1) + _anticommutator(momentum, shift, dphi1)
    mass = kappa * length
    image0 = -1j * kinetic1 + mass[:, None] * phi1 - shift0
    image1 = 1j * kinetic0 + mass[:, None] * phi0 - shift1
    dimage0 = (
        -1j * dkinetic1
        + (kappa * dlength)[:, None] * phi1
        + mass[:, None] * dphi1
        - dshift0
    )
    dimage1 = (
        1j * dkinetic0
        + (kappa * dlength)[:, None] * phi0
        + mass[:, None] * dphi0
        - dshift1
    )
    return image0, image1, dimage0, dimage1


def _moments_tangent(momentum, phi0, phi1, dphi0, dphi1, occupations, delta_occupations):
    """Same nodal moments as ``column_moments``, plus their real directional parts."""
    momentum0 = momentum @ phi0
    momentum1 = momentum @ phi1
    dmomentum0 = momentum @ dphi0
    dmomentum1 = momentum @ dphi1
    kinetic = np.conjugate(phi0) * (-1j * momentum1) + np.conjugate(phi1) * (1j * momentum0)
    dkinetic = (
        np.conjugate(dphi0) * (-1j * momentum1)
        + np.conjugate(phi0) * (-1j * dmomentum1)
        + np.conjugate(dphi1) * (1j * momentum0)
        + np.conjugate(phi1) * (1j * dmomentum0)
    )
    spin = np.conjugate(phi0) * phi1 + np.conjugate(phi1) * phi0
    dspin = (
        np.conjugate(dphi0) * phi1
        + np.conjugate(phi0) * dphi1
        + np.conjugate(dphi1) * phi0
        + np.conjugate(phi1) * dphi0
    )
    current = np.conjugate(phi0) * momentum0 + np.conjugate(phi1) * momentum1
    dcurrent = (
        np.conjugate(dphi0) * momentum0
        + np.conjugate(phi0) * dmomentum0
        + np.conjugate(dphi1) * momentum1
        + np.conjugate(phi1) * dmomentum1
    )
    weights = occupations[None, :]
    dweights = delta_occupations[None, :]
    kinetic_density = np.sum(dkinetic * weights + kinetic * dweights, axis=1).real
    spin_density = np.sum(dspin * weights + spin * dweights, axis=1).real
    current_density = np.sum(dcurrent * weights + current * dweights, axis=1).real
    return kinetic_density, spin_density, current_density


def _force_directional(system, kinetic, dkinetic, spin, dspin, _current, dcurrent, radial, dradial, length, dlength):
    """``M`` once. ``force_L`` uses ``Q``; ``force_Q`` uses both ``L`` and ``Q``."""
    if int(system.multiplicity) != 4 * int(system.kappa):
        raise ValueError("4 kappa multiplicity drifted")
    multiplicity = float(system.multiplicity)
    delta_l = multiplicity * (
        dkinetic / radial - kinetic * dradial / radial ** 2 + system.kappa * dspin
    )
    delta_q = -multiplicity * (
        dlength * kinetic / radial ** 2
        + length * dkinetic / radial ** 2
        - 2 * length * kinetic * dradial / radial ** 3
    )
    delta_beta = -multiplicity * dcurrent
    return delta_l, delta_q, delta_beta


def _geometry_tangents(system, state, tangent, length, dlength, shift, dshift):
    """Differentiate the discrete Hamilton equations. ``D`` stays outside products."""
    derivative = system.derivative
    amplitude = system.A
    weyl = system.C_W
    stiffness = float(feedback_Z(amplitude))
    alpha = float(alpha_of(weyl))
    force_r, force_chi = partial_F(state.r, amplitude, weyl)
    force_chi = float(force_chi)
    delta_force_r = -8 * np.pi * amplitude * tangent.r
    pi = state.p_r - force_r * state.p_chi / force_chi
    delta_pi = (
        tangent.p_r
        - delta_force_r * state.p_chi / force_chi
        - force_r * tangent.p_chi / force_chi
    )
    profile = feedback_F(state.r, state.chi, amplitude, weyl)
    delta_profile = force_r * tangent.r + force_chi * tangent.chi
    potential = feedback_V(state.r, state.chi, amplitude, weyl, system.C_F, system.flux)
    delta_potential = (
        16 * np.pi * amplitude * state.r * tangent.r
        - alpha * (4 + 2 * state.chi) * tangent.chi
    )
    radius = state.Q
    delta_radius = tangent.Q
    radius_x = derivative @ state.r
    delta_radius_x = derivative @ tangent.r
    chi_x = derivative @ state.chi
    delta_chi_x = derivative @ tangent.chi
    momentum_x = derivative @ state.p_Q
    delta_momentum_x = derivative @ tangent.p_Q
    lapse_x = derivative @ length
    delta_lapse_x = derivative @ dlength
    profile_x = derivative @ profile
    delta_profile_x = derivative @ delta_profile
    pi_over = pi / (2 * stiffness * radius)
    delta_pi_over = delta_pi / (2 * stiffness * radius) - pi * delta_radius / (2 * stiffness * radius ** 2)
    velocity = length * pi_over
    delta_velocity = dlength * pi_over + length * delta_pi_over
    delta_q = (
        dlength * state.p_chi / (2 * force_chi)
        + length * tangent.p_chi / (2 * force_chi)
        + derivative @ (dshift * radius + shift * delta_radius)
    )
    delta_r = delta_velocity + dshift * radius_x + shift * delta_radius_x
    delta_chi = (
        dlength * state.p_Q / (2 * force_chi)
        + length * tangent.p_Q / (2 * force_chi)
        - (delta_force_r / force_chi) * velocity
        - (force_r / force_chi) * delta_velocity
        + dshift * chi_x
        + shift * delta_chi_x
    )
    delta_p_q = (
        dlength * pi ** 2 / (4 * stiffness * radius ** 2)
        + length * pi * delta_pi / (2 * stiffness * radius ** 2)
        - length * pi ** 2 * delta_radius / (2 * stiffness * radius ** 3)
        + dlength * stiffness * radius_x ** 2 / radius ** 2
        + 2 * length * stiffness * radius_x * delta_radius_x / radius ** 2
        - 2 * length * stiffness * radius_x ** 2 * delta_radius / radius ** 3
        + dlength * potential
        + length * delta_potential
        + 2 * (delta_profile_x * lapse_x + profile_x * delta_lapse_x) / radius ** 2
        - 4 * profile_x * lapse_x * delta_radius / radius ** 3
        + dshift * momentum_x
        + shift * delta_momentum_x
    )
    pi_radius = (8 * np.pi * amplitude / force_chi) * state.p_chi
    delta_pi_radius = (8 * np.pi * amplitude / force_chi) * tangent.p_chi
    delta_pi_term = (
        -dlength * pi * pi_radius / (2 * stiffness * radius)
        - length * delta_pi * pi_radius / (2 * stiffness * radius)
        - length * pi * delta_pi_radius / (2 * stiffness * radius)
        + length * pi * pi_radius * delta_radius / (2 * stiffness * radius ** 2)
    )
    delta_inner = 2 * stiffness * (
        dlength * radius_x / radius
        + length * delta_radius_x / radius
        - length * radius_x * delta_radius / radius ** 2
    )
    weight = 2 * lapse_x / radius
    delta_weight = 2 * delta_lapse_x / radius - 2 * lapse_x * delta_radius / radius ** 2
    delta_p_r = (
        delta_pi_term
        + derivative @ delta_inner
        + dlength * radius * (16 * np.pi * amplitude * state.r)
        + length * delta_radius * (16 * np.pi * amplitude * state.r)
        + length * radius * (16 * np.pi * amplitude * tangent.r)
        + delta_force_r * (derivative @ weight)
        + force_r * (derivative @ delta_weight)
        + derivative @ (dshift * state.p_r + shift * tangent.p_r)
    )
    delta_p_chi = (
        dlength * radius * (-alpha * (4 + 2 * state.chi))
        + length * delta_radius * (-alpha * (4 + 2 * state.chi))
        + length * radius * (-alpha * 2 * tangent.chi)
        + force_chi * (derivative @ delta_weight)
        + derivative @ (dshift * state.p_chi + shift * tangent.p_chi)
    )
    profile_over = profile_x / radius
    delta_profile_over = delta_profile_x / radius - profile_x * delta_radius / radius ** 2
    delta_hamilton = (
        (tangent.p_Q * state.p_chi + state.p_Q * tangent.p_chi) / (2 * force_chi)
        + pi * delta_pi / (2 * stiffness * radius)
        - pi ** 2 * delta_radius / (4 * stiffness * radius ** 2)
        + 2 * stiffness * radius_x * delta_radius_x / radius
        - stiffness * radius_x ** 2 * delta_radius / radius ** 2
        - delta_radius * potential
        - radius * delta_potential
        - 2 * (derivative @ delta_profile_over)
    )
    return delta_q, delta_r, delta_chi, delta_p_q, delta_p_r, delta_p_chi, delta_hamilton


def _owned_moments(source):
    """Base moments from ``source_from_columns``. Spin is the real part used by the forces."""
    return source["K"], np.real(source["S1"]), source["Pmom"]


def rate_jacobian_vector(grid, state, tangent, *, include_matter_force=True):
    """Analytic ``DR · tangent`` of ``compose_fine_hamiltonian``.

    Occupations are the live system weights. ``tangent.occupations`` is only
    their directional increment. The representative column rate stays ``-iHΦ``;
    multiplicity does not multiply it.
    """
    tangent = _validate_galerkin(grid, state, tangent)
    _coarse, bundle = galerkin.compose_fine_hamiltonian(grid, state, include_matter_force)
    fine_state = bundle["fine_state"]
    fine_system = bundle["fine_system"]
    if bundle["source"]["multiplicity_applied_once"] is not True:
        raise ValueError("source multiplicity was not applied once")
    fine_tangent = _prolong_tangent(grid, tangent)
    if grid.gauge == "conformal":
        length, dlength = fine_state.Q, fine_tangent.Q
        shift = np.zeros(grid.nq)
        dshift = np.zeros(grid.nq)
    else:
        length, dlength = fine_system.length_density, np.zeros(grid.nq)
        shift, dshift = fine_system.shift, np.zeros(grid.nq)
    geometric = _geometry_tangents(
        fine_system, fine_state, fine_tangent, length, dlength, shift, dshift,
    )
    delta_q, delta_r, delta_chi, delta_p_q, delta_p_r, delta_p_chi, delta_hamilton = geometric
    if grid.gauge == "conformal":
        delta_p_q = delta_p_q - delta_hamilton
    occupations = np.asarray(fine_system.occupations, dtype=float)
    kinetic, spin, current = _owned_moments(bundle["source"])
    dkinetic, dspin, dcurrent = _moments_tangent(
        fine_system.momentum,
        fine_state.phi0, fine_state.phi1,
        fine_tangent.phi0, fine_tangent.phi1,
        occupations, fine_tangent.occupations,
    )
    delta_l, delta_q_force, delta_beta = _force_directional(
        fine_system, kinetic, dkinetic, spin, dspin, current, dcurrent,
        fine_state.Q, fine_tangent.Q, length, dlength,
    )
    if include_matter_force:
        delta_p_q = delta_p_q - delta_q_force / grid.dx_q
        if grid.gauge == "conformal":
            delta_p_q = delta_p_q - delta_l / grid.dx_q
    _image0, _image1, dimage0, dimage1 = _images_tangent(
        fine_system.momentum,
        fine_state.phi0, fine_state.phi1,
        fine_tangent.phi0, fine_tangent.phi1,
        length, dlength, fine_state.Q, fine_tangent.Q, shift, dshift, fine_system.kappa,
    )
    # Representative block. Not -M i H.
    delta_phi0 = grid.U_f.conj().T @ (-1j * dimage0)
    delta_phi1 = grid.U_f.conj().T @ (-1j * dimage1)
    delta_q_coarse = galerkin.pull_geometry(grid, delta_q)
    lifted_delta = grid.A_g @ delta_q_coarse
    lapse_delta = lifted_delta if grid.gauge == "conformal" else np.zeros(grid.nq)
    source = bundle["source"]
    fieldwork = float(np.sum(
        delta_q_force * bundle["lifted_Q"]
        + source["force_Q"] * lifted_delta
        + delta_l * bundle["lapse_dot"]
        + source["force_L"] * lapse_delta
    ))
    return RateTangent(
        galerkin.pull_geometry(grid, delta_q),
        galerkin.pull_geometry(grid, delta_r),
        galerkin.pull_geometry(grid, delta_chi),
        galerkin.pull_geometry(grid, delta_p_q),
        galerkin.pull_geometry(grid, delta_p_r),
        galerkin.pull_geometry(grid, delta_p_chi),
        delta_phi0,
        delta_phi1,
        delta_l,
        delta_q_force,
        delta_beta,
        fieldwork,
    )


def _nodal_tangent(pair, tangent):
    nested._check(pair, tangent)
    nodal = nested.reconstruct_state(pair, tangent)
    return StateTangent(
        nodal.Q, nodal.r, nodal.chi, nodal.p_Q, nodal.p_r, nodal.p_chi,
        np.asarray(tangent.phi0, dtype=complex),
        np.asarray(tangent.phi1, dtype=complex),
        _real_vector(tangent.occupations, pair.grid.fine.occupations.size, "occupations"),
    )


def _encode_rate(pair, tangent):
    matrix, spacing = pair.geometry_map, pair.grid.dx_g
    return RateTangent(
        matrix.T @ tangent.Q,
        matrix.T @ tangent.r,
        matrix.T @ tangent.chi,
        spacing * matrix.T @ tangent.p_Q,
        spacing * matrix.T @ tangent.p_r,
        spacing * matrix.T @ tangent.p_chi,
        tangent.phi0,
        tangent.phi1,
        tangent.force_L,
        tangent.force_Q,
        tangent.force_beta,
        tangent.fieldwork_power,
    )


def nested_rate_jacobian_vector(pair, state, tangent, *, include_matter_force=True):
    """Nested image of the same analytic rate tangent.

    ``W`` and ``Δx_g`` are constant, so the chain rule is the owned encode
    map applied to the nodal Jacobian-vector product.
    """
    if pair.grid.gauge != "conformal":
        raise ValueError("nested rate tangent requires the conformal gauge")
    nested._check(pair, state)
    nodal_state = nested.reconstruct_state(pair, state)
    nodal_tangent = _nodal_tangent(pair, tangent)
    nodal_rate = rate_jacobian_vector(
        pair.grid, nodal_state, nodal_tangent, include_matter_force=include_matter_force,
    )
    return _encode_rate(pair, nodal_rate)


def _preparation_domain(grid, state):
    if grid.gauge != "conformal":
        raise ValueError("linearized radius preparation requires the conformal gauge")
    expected = grid.fine.calibration["b0"] / grid.fine.calibration["a0"]
    if float(np.max(np.abs(state.Q - expected))) > 1e-10:
        raise ValueError("linearized radius preparation requires the owned constant Q")
    for name in ("chi", "p_Q", "p_r", "p_chi"):
        if float(np.max(np.abs(getattr(state, name)))) > 1e-10:
            raise ValueError("linearized radius preparation requires zero chi and momenta")


def prepared_radius_tangent(grid, state, tangent):
    """``J δr = -pull(δρ)`` with the owned projected radius Jacobian.

    ``ρ = force_L / Δx`` is the source density at fixed owned ``Q``. The
    base radius is not moved and no Newton loop is run. A singular owned
    Jacobian is a missing primitive, not a finite-difference replacement.
    """
    _preparation_domain(grid, state)
    tangent = _validate_galerkin(grid, state, tangent)
    for name in ("Q", "chi", "p_Q", "p_r", "p_chi"):
        if float(np.max(np.abs(getattr(tangent, name)))) > 0:
            raise ValueError("source preparation tangent keeps Q, chi and momenta fixed")
    fine_state = galerkin.prolong_state(grid, state)
    fine_system = galerkin.active_fine_system(grid, fine_state)
    fine_tangent = _prolong_tangent(grid, tangent)
    occupations = np.asarray(fine_system.occupations, dtype=float)
    kinetic, spin, current = _owned_moments(coupling.source_from_columns(fine_system, fine_state))
    dkinetic, dspin, dcurrent = _moments_tangent(
        fine_system.momentum,
        fine_state.phi0, fine_state.phi1,
        fine_tangent.phi0, fine_tangent.phi1,
        occupations, fine_tangent.occupations,
    )
    # L equals the owned constant Q on this domain, and δL = δQ = 0.
    delta_l, _delta_q, _delta_beta = _force_directional(
        fine_system, kinetic, dkinetic, spin, dspin, current, dcurrent,
        fine_state.Q, fine_tangent.Q, fine_state.Q, fine_tangent.Q,
    )
    delta_rho = delta_l / grid.dx_q
    residual_tangent = galerkin.pull_geometry(grid, delta_rho)
    jacobian = galerkin._projected_radius_jacobian(grid, np.asarray(state.r, dtype=float))
    try:
        delta_r = np.linalg.solve(jacobian, -residual_tangent)
    except np.linalg.LinAlgError:
        return PreparedRadiusTangent(
            available=False, delta_r=None,
            missing_primitive="singular_projected_radius_jacobian",
            linear_residual_max=None,
        )
    if not np.isfinite(delta_r).all():
        return PreparedRadiusTangent(
            available=False, delta_r=None,
            missing_primitive="nonfinite_projected_radius_tangent",
            linear_residual_max=None,
        )
    leftover = jacobian @ delta_r + residual_tangent
    scale = max(1.0, float(np.linalg.norm(residual_tangent, np.inf)))
    return PreparedRadiusTangent(
        available=True,
        delta_r=delta_r,
        missing_primitive=None,
        linear_residual_max=float(np.max(np.abs(leftover)) / scale),
    )


def prepared_nested_radius_tangent(pair, state, tangent):
    """Canonical ``δa_r = Wᵀ δr`` of the same linear solve."""
    nodal_state = nested.reconstruct_state(pair, state)
    nodal_tangent = _nodal_tangent(pair, tangent)
    prepared = prepared_radius_tangent(pair.grid, nodal_state, nodal_tangent)
    if not prepared.available:
        return prepared
    delta_a = pair.geometry_map.T @ prepared.delta_r
    return PreparedRadiusTangent(
        available=True,
        delta_r=delta_a,
        missing_primitive=None,
        linear_residual_max=prepared.linear_residual_max,
    )


def _fine_pair(grid, state):
    fine = galerkin.prolong_state(grid, state)
    return fine, galerkin.active_fine_system(grid, fine)


def _number_density(grid, phi0, phi1, occupations):
    weights = np.asarray(occupations, dtype=float)
    probability = np.sum((np.abs(phi0) ** 2 + np.abs(phi1) ** 2) * weights[None, :], axis=1)
    return probability / grid.dx_q


def _density_tangent(phi0, phi1, dphi0, dphi1, occupations, delta_occupations, spacing):
    weights = occupations[None, :]
    dweights = delta_occupations[None, :]
    base = np.abs(phi0) ** 2 + np.abs(phi1) ** 2
    delta = (
        2 * np.real(np.conjugate(phi0) * dphi0 + np.conjugate(phi1) * dphi1)
    )
    probability = np.sum(delta * weights + base * dweights, axis=1)
    return probability / spacing


def child_regional_content(grid, state, tangent=None, occupations=None):
    """Integral of the one-body density of ``C`` over the child collar ``(1, 3)``.

    This is the primary readout. It is not the trace of the child 2×2 kernel.
    """
    weights = grid.fine.occupations if occupations is None else np.asarray(occupations, dtype=float)
    fine, _system = _fine_pair(grid, state)
    density = _number_density(grid, fine.phi0, fine.phi1, weights)
    content = nested.interval_integral(grid, density, nested.CHILD_INTERVAL)
    if tangent is None:
        return float(content)
    fine_tangent = _prolong_tangent(grid, _validate_galerkin(grid, state, tangent))
    delta_density = _density_tangent(
        fine.phi0, fine.phi1, fine_tangent.phi0, fine_tangent.phi1,
        weights, fine_tangent.occupations, grid.dx_q,
    )
    return float(content), float(nested.interval_integral(grid, delta_density, nested.CHILD_INTERVAL))


def child_proper_mean_r(grid, state, tangent=None):
    """Child proper mean of ``r``, the secondary readout.

    The weight is the same ``r Q`` clock density used by the nested metrics.
    In the conformal chart that density is ``r L``.
    """
    fine, _system = _fine_pair(grid, state)
    radial = fine.r * fine.Q
    denominator = nested.interval_integral(grid, radial, nested.CHILD_INTERVAL)
    numerator = nested.interval_integral(grid, radial * fine.r, nested.CHILD_INTERVAL)
    if denominator <= 0:
        raise ValueError("child proper length is not positive")
    mean = numerator / denominator
    if tangent is None:
        return float(mean)
    validated = _validate_galerkin(grid, state, tangent)
    fine_tangent = _prolong_tangent(grid, validated)
    delta_radial = fine_tangent.r * fine.Q + fine.r * fine_tangent.Q
    delta_denominator = nested.interval_integral(grid, delta_radial, nested.CHILD_INTERVAL)
    delta_numerator = nested.interval_integral(
        grid, delta_radial * fine.r + radial * fine_tangent.r, nested.CHILD_INTERVAL,
    )
    delta_mean = (delta_numerator - mean * delta_denominator) / denominator
    return float(mean), float(delta_mean)


def child_clock_rate(grid, state, tangent=None, coordinate_duration=0.0):
    """Child clock ``τ̇ = (r Q)(x=2)`` and, with a tangent, ``δτ̇`` and ``δτ``.

    ``δτ = δτ̇ Δt`` holds the rate fixed on a declared coordinate interval.
    That product is not an integral along an evolved history. ``Δt = 0``
    gives ``δτ = 0``.
    """
    if not np.isfinite(coordinate_duration):
        raise ValueError("coordinate duration must be finite")
    fine, _system = _fine_pair(grid, state)
    radial = fine.r * fine.Q
    tau_dot = float(nested._periodic_values(grid, radial, (CHILD_CLOCK,))[0])
    if tangent is None:
        return tau_dot
    validated = _validate_galerkin(grid, state, tangent)
    fine_tangent = _prolong_tangent(grid, validated)
    delta_radial = fine_tangent.r * fine.Q + fine.r * fine_tangent.Q
    delta_tau_dot = float(nested._periodic_values(grid, delta_radial, (CHILD_CLOCK,))[0])
    return tau_dot, delta_tau_dot, delta_tau_dot * float(coordinate_duration)


def _velocity_tangent(grid, rate):
    return StateTangent(
        rate.Q, rate.r, rate.chi, rate.p_Q, rate.p_r, rate.p_chi,
        rate.phi0, rate.phi1, np.zeros(grid.fine.occupations.size),
    )


def galerkin_readout(grid, state, tangent, *, coordinate_duration=0.0):
    """Primary child content and secondary proper mean ``r``, at fixed ``t`` and matched ``τ``."""
    tangent = _validate_galerkin(grid, state, tangent)
    content, delta_content = child_regional_content(grid, state, tangent)
    mean_r, delta_mean = child_proper_mean_r(grid, state, tangent)
    tau_dot, delta_tau_dot, delta_tau = child_clock_rate(
        grid, state, tangent, coordinate_duration=coordinate_duration,
    )
    rate = galerkin.rates(grid, state)
    velocity = _velocity_tangent(grid, rate)
    _content, content_dot = child_regional_content(grid, state, velocity)
    _mean, mean_dot = child_proper_mean_r(grid, state, velocity)
    return {
        "primary_readout": PRIMARY_READOUT,
        "secondary_readout": SECONDARY_READOUT,
        "child_regional_content": content,
        "delta_child_regional_content_t": delta_content,
        "child_regional_content_dot": content_dot,
        "delta_child_regional_content_tau": matched_tau_correction(
            delta_content, content_dot, delta_tau, tau_dot,
        ),
        "child_proper_mean_r": mean_r,
        "delta_child_proper_mean_r_t": delta_mean,
        "child_proper_mean_r_dot": mean_dot,
        "delta_child_proper_mean_r_tau": matched_tau_correction(
            delta_mean, mean_dot, delta_tau, tau_dot,
        ),
        "tau_dot": tau_dot,
        "delta_tau_dot": delta_tau_dot,
        "delta_tau": delta_tau,
        "coordinate_duration": float(coordinate_duration),
        "frozen_rate_interval": True,
        "evolved_clock_integral": False,
    }


def nested_readout(pair, state, tangent, *, coordinate_duration=0.0):
    """Same readouts after the owned reconstruction of one nested state."""
    nodal_state = nested.reconstruct_state(pair, state)
    nodal_tangent = _nodal_tangent(pair, tangent)
    return galerkin_readout(
        pair.grid, nodal_state, nodal_tangent, coordinate_duration=coordinate_duration,
    )


def child_modal_kernel(covariance, child_columns):
    """2×2 child-observer block of ``C``.

    This kernel is modal. It is not, by definition, the spatial exterior of
    the child collar, and a nonzero entry is not spatial leakage.
    """
    child = np.asarray(child_columns, dtype=complex)
    covariance = _covariance(covariance)
    if child.ndim != 2 or child.shape[0] != covariance.shape[0] or child.shape[1] != 2:
        raise ValueError("child observer kernel needs two columns in the covariance space")
    if not np.isfinite(child).all():
        raise ValueError("nonfinite child observer")
    gram = child.conj().T @ child
    if np.max(np.abs(gram - np.eye(2))) > 1e-8:
        raise ValueError("child observer columns must be orthonormal")
    kernel = child.conj().T @ covariance @ child
    cross = child.conj().T @ covariance - kernel @ child.conj().T
    return {
        "kernel": kernel,
        "shape": (2, 2),
        "spatial_outside": False,
        "spatial_outside_by_definition": False,
        "initial_cross": cross,
        "initial_cross_frobenius": float(np.linalg.norm(cross)),
        "initial_cross_is_spatial_exterior": False,
    }


def future_effective_stress_specification(covariance, child_columns):
    """CTP future stress, specified and not evaluated.

    The future value needs the complementary-branch force, the initial cross
    and the retarded change of state. This wave evaluates only the initial
    modal cross of the full Gaussian. It does not add an induced coefficient,
    and it does not fill the retarded state by a finite difference.
    """
    modal = child_modal_kernel(covariance, child_columns)
    return {
        "implemented": False,
        "missing_primitive": "ctp_future_effective_stress",
        "effective_stress": None,
        "definition": "complementary_force + initial_cross + retarded_delta_state",
        "parts": {
            "complementary_force": {
                "implemented": False,
                "reason": (
                    "the future complementary CTP branch is not a history in this wave; "
                    "equal-time force_L and force_Q stay inside the rate and are not a substitute"
                ),
            },
            "initial_cross": {
                "implemented": True,
                "frobenius": modal["initial_cross_frobenius"],
                "spatial_exterior": False,
                "meaning": "child observer block of C against the rest of the one-particle space",
            },
            "retarded_delta_state": {
                "implemented": False,
                "missing_primitive": "retarded_delta_state",
                "value": None,
                "reason": "no retarded transport of a full-state tangent is formed in this wave",
            },
        },
        "induced_force_included": False,
        "induced_owner": INDUCED_OWNER,
        "kernel_2x2_is_spatial_outside": False,
        "finite_difference_used_as_analytic": False,
    }
