"""Reusable parent response adapter for the leading six-field chart.

The evolution, the analytic state Jacobian and the RK4 step stay in
``nsc_discovery_leading_einstein``. Source-force and density directions stay
in ``nsc_discovery_response``. The inherited direct and sequential column
split stays in ``nsc_discovery_grandchild.route_rates`` when the observer
widths are the saved ``(2, 2, 4)`` rows. Variable-rank columns do not call
``build_pair`` or the coupled-memory rank-6 guard.

The carrier centre is ``L/2``. Windows are signed-distance collars on that
centre. The primary readout is the child spatial window. The fixed reference
columns are a modal observer and are not that window.

The initial-preparation derivative is a labelled finite-difference
approximation with two step sizes. It is not an analytic preparation
Jacobian. This module does not run a held nonlinear measurement.
"""
from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, replace
import json
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_limits

from . import nsc_discovery_backend as backend
from . import nsc_discovery_coupled_memory as memory
from . import nsc_discovery_episode as episode
from . import nsc_discovery_extent as extent
from . import nsc_discovery_grandchild as grandchild
from . import nsc_discovery_leading_einstein as leading
from . import nsc_discovery_parent_step_control as step_control
from . import nsc_discovery_response as response
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin


SCHEMA = "NSC-DISCOVERY-PARENT-RESPONSE-v1"
METHOD = "leading_analytic_jvp"
PREPARATION_METHOD = "finite_difference_full_preparation"
CARRIER_LENGTH = float(coupling.PERIOD)
CHILD_RADIUS = 0.5
COLLAR_RADIUS = 1.0
PARENT_RADIUS = 3.0
ANNULUS_INNER = 1.2
ANNULUS_OUTER = 3.0
INHERITED_WIDTHS = (2, 2, 4)
CHUNK_LIMIT_BYTES = 64 * 1024 * 1024
DENSE_RADIUS_JACOBIAN_MAX = 31
BASELINE_SEED_AMPLITUDE = 0.0013
PHYSICAL_BINDINGS = ("owner1_prepare_parent", "owner2_cached_episode")
SEALED_DIRECTORY_NAMES = frozenset({
    "nsc-discovery-parent-traction-v1",
    "nsc-discovery-leading-einstein-v1",
    "nsc-discovery-leading-einstein-v2",
    "nsc-discovery-dynamic-preparation-v2",
})
STATE_FIELDS = leading.FIELDS
REAL_FIELDS = leading.REAL_FIELDS
LAB = Path(__file__).resolve().parents[2]
REPO = LAB.parent
DEVELOPMENT_ROOT = LAB / "results" / "development"

IMPLEMENTED = (
    "parent_holder",
    "signed_distance_windows",
    "response_rates",
    "advance_tangent",
    "centre_clock_variation",
    "nonlinear_endpoint_clock_bracket",
    "child_regional_content",
    "proper_localization",
    "fixed_mode_observer",
    "explicit_normalization_derivative",
    "preparation_tangent_finite_difference",
    "hierarchy_rates",
    "full_vs_reduced",
    "prepare_predict_measure_check",
    "owner_callbacks",
    "evolve_prepared_source",
)
NOT_IMPLEMENTED = (
    "analytic_preparation_jacobian",
    "scientific_pass_minted_by_this_adapter",
)
DEPENDENCIES = (
    {
        "owner": "nsc_discovery_leading_einstein",
        "uses": "State, rates, rk4_step, jvp, fine_jvp, encode, decode, constraint_arrays, check_chart",
        "equations_copied": False,
    },
    {
        "owner": "nsc_discovery_response",
        "uses": "matched_tau_correction, _moments_tangent, _force_directional, _images_tangent, _density_tangent",
        "equations_copied": False,
    },
    {
        "owner": "nsc_discovery_grandchild",
        "uses": "initial_routes, reconstruct_route, route_rates for widths (2, 2, 4)",
        "equations_copied": False,
    },
    {
        "owner": "nsc_discovery_coupled_memory",
        "uses": "split interpretation only; _require_rank6 is not called",
        "equations_copied": False,
    },
    {
        "owner": "nsc_spherical_coupling",
        "uses": "apply_dirac and source_from_columns",
        "equations_copied": False,
    },
    {
        "owner": "nsc_discovery_backend",
        "uses": "FFT carrier classes at one worker; make_fft_grid rank-6 guard is not called",
        "equations_copied": False,
    },
    {
        "owner": "nsc_discovery_extent",
        "uses": "real_interval_integral and real_periodic_values",
        "equations_copied": False,
    },
    {
        "owner": "nsc_discovery_parent_step_control",
        "uses": "step_restriction and stable_timestep; admission is not copied",
        "equations_copied": False,
    },
    {
        "owner": "nsc_discovery_parent.prepare_parent",
        "uses": "OWNER1 callback; not called by a default or a manufactured fixture",
        "equations_copied": False,
    },
    {
        "owner": "nsc_discovery_parent_episode",
        "uses": "OWNER2 load_parent_record -> (record, arrays), reconstruct_parent_pair(arrays, record), state_from_arrays",
        "equations_copied": False,
    },
    {
        "owner": "OWNER3 observables",
        "uses": "not imported; spatial windows here are not the fixed-mode observer",
        "equations_copied": False,
    },
)


class BindingUnavailable(RuntimeError):
    """An initial state or input episode is not bound. The stage is open, not a pass."""


class AuthorizationOpen(RuntimeError):
    """A development stage is open until the producer is frozen and the binding is physical."""


@dataclass(frozen=True)
class ParentTangent:
    """Directional increment of the six leading fields plus source weights.

    Momenta are canonical ``π = Δx_g W.T p`` in the same frame as ``leading.State``.
    ``occupations`` is ``δc``, the direction of the live weights, not the weights.
    """

    Q: np.ndarray
    r: np.ndarray
    p_Q: np.ndarray
    p_r: np.ndarray
    phi0: np.ndarray
    phi1: np.ndarray
    occupations: np.ndarray


@dataclass(frozen=True)
class ResponseRate:
    """Directional image of one leading rate, including force and operator pieces."""

    Q: np.ndarray
    r: np.ndarray
    p_Q: np.ndarray
    p_r: np.ndarray
    phi0: np.ndarray
    phi1: np.ndarray
    force_L: np.ndarray
    force_Q: np.ndarray
    force_beta: np.ndarray
    operator_phi0: np.ndarray
    operator_phi1: np.ndarray
    fieldwork_power: float
    method: str = METHOD
    weight_channel_included: bool = True
    dense_propagator_stored: bool = False


@dataclass(frozen=True)
class ParentHolder:
    """ParentPair-compatible carrier. ``W`` may be the identity.

    ``grid.fine.occupations`` is the owned weight vector, installed with
    ``dataclasses.replace``. Columns may have any positive rank.
    """

    grid: galerkin.GalerkinGrid
    geometry_map: np.ndarray
    weights: np.ndarray
    reference_columns: np.ndarray
    source_columns: np.ndarray
    source_metadata: dict
    geometry_metadata: dict
    child_interval: tuple
    parent_interval: tuple
    clock_locations: tuple


def specification():
    """Static contract. Calling it does not read an episode or write a file."""
    return {
        "schema": SCHEMA,
        "method": METHOD,
        "preparation_method": PREPARATION_METHOD,
        "analytic_preparation_jacobian": False,
        "finite_difference_used_as_analytic": False,
        "primary_readout": "child_regional_content",
        "proper_localization": "child_proper_length_and_probability",
        "fixed_mode_observer_is_spatial_window": False,
        "canonical_momentum": "pi = dx_g W.T p",
        "state_fields": list(STATE_FIELDS),
        "centre": "L/2",
        "windows": {
            "child": "|s|<=0.5",
            "protected_collar": "|s|<=1",
            "parent": "|s|<=3",
            "parent_annulus": "1.2<=|s|<=3",
        },
        "one_physical_metric": True,
        "full_geometry_band_active": True,
        "ambient_field_active": True,
        "native_fft_workers": 1,
        "rank6_guard_called": False,
        "build_pair_called": False,
        "dense_propagator_stored": False,
        "production_authorized_by_default": False,
        "baseline_seed_amplitude": BASELINE_SEED_AMPLITUDE,
        "baseline_seed_is_a_measurement_veto": False,
        "held_direction": "sealed record chooses population or gradient before measurement",
        "source_may_be_nonorthogonal": True,
        "modal_observer_must_be_orthonormal": True,
        "implemented": list(IMPLEMENTED),
        "not_implemented": list(NOT_IMPLEMENTED),
        "dependencies": list(DEPENDENCIES),
        "output_limit_bytes": CHUNK_LIMIT_BYTES,
    }


def dependencies():
    """Precise owners this adapter calls. No paper claim is attached."""
    return list(DEPENDENCIES)


def signed_distance(length, coordinates):
    """Periodic offset from ``L/2``, returned in ``(-L/2, L/2]``."""
    period = float(length)
    if not np.isfinite(period) or period <= 0.0:
        raise ValueError("carrier length must be positive")
    centre = period / 2.0
    offset = np.asarray(coordinates, dtype=float) - centre
    return np.mod(offset + period / 2.0, period) - period / 2.0


def regional_map(length=CARRIER_LENGTH):
    """Closed coordinate windows for the signed-distance collars."""
    period = float(length)
    centre = period / 2.0

    def interval(radius):
        left = centre - float(radius)
        right = centre + float(radius)
        if not (0.0 <= left < right <= period):
            raise ValueError("signed-distance window does not fit the carrier")
        return (float(left), float(right))

    annulus = (
        (float(centre - ANNULUS_OUTER), float(centre - ANNULUS_INNER)),
        (float(centre + ANNULUS_INNER), float(centre + ANNULUS_OUTER)),
    )
    if annulus[0][0] < 0.0 or annulus[1][1] > period:
        raise ValueError("parent annulus does not fit the carrier")
    return {
        "length": period,
        "centre": float(centre),
        "child_interval": interval(CHILD_RADIUS),
        "protected_collar_interval": interval(COLLAR_RADIUS),
        "parent_interval": interval(PARENT_RADIUS),
        "parent_annulus_intervals": annulus,
        "child_radius": CHILD_RADIUS,
        "collar_radius": COLLAR_RADIUS,
        "parent_radius": PARENT_RADIUS,
        "annulus_inner": ANNULUS_INNER,
        "annulus_outer": ANNULUS_OUTER,
        "signed_distance": "periodic offset from L/2 in (-L/2, L/2]",
    }


def _readonly(value, dtype=None):
    array = np.array(value, dtype=dtype, copy=True)
    array.setflags(write=False)
    return array


def _real_vector(values, size, name):
    array = np.asarray(values, dtype=float)
    if np.iscomplexobj(values) or array.shape != (size,) or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a finite real vector of length {size}")
    return np.array(array, copy=True)


def _columns(values, shape, name):
    array = np.asarray(values, dtype=complex)
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a finite complex array of shape {shape}")
    return np.array(array, copy=True)


@contextmanager
def _one_thread():
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        yield


def install_native_fft(grid):
    """Install the backend FFT carriers at the current one-worker limit.

    ``make_fft_grid`` refuses every occupation width other than 6. This
    adapter does not call that guard. The operator classes are the backend
    classes, so the derivative and the antiperiodic momentum are not copied
    into a second equation.
    """
    if backend.operator_backend(grid) == "fft":
        return grid
    if not isinstance(grid, galerkin.GalerkinGrid):
        raise TypeError("native FFT installation expects a Galerkin grid")
    derivative = backend.PeriodicDerivative(grid.nq, grid.length)
    momentum = backend.AntiperiodicMomentum(grid.nq, grid.length)
    geometry = backend.GeometryProlongation(grid.ng, grid.nq, grid.length)
    columns = backend.SpinorCarrier(grid.nf, grid.nq, grid.length, canonical=True)
    interpolation = backend.SpinorCarrier(grid.nf, grid.nq, grid.length, canonical=False)
    fine = replace(grid.fine, derivative=derivative, momentum=momentum)
    fft_grid = replace(
        grid,
        A_g=geometry,
        A_f=interpolation,
        U_f=columns,
        derivative_on_geometry=backend.Composed(derivative, geometry),
        fine=fine,
    )
    backend._assert_carrier_isometry(fft_grid)
    return fft_grid


def _install_weights(grid, weights):
    """Owned replacement of the quadrature occupations. No dummy column is added."""
    live = np.array(weights, dtype=float, copy=True)
    if live.ndim != 1 or live.size < 1 or not np.isfinite(live).all():
        raise ValueError("source weights must be a finite nonempty vector")
    if np.min(live) < 0.0 or np.max(live) > 1.0:
        raise ValueError("source weights must lie in [0, 1]")
    return replace(grid, fine=replace(grid.fine, occupations=live))


def _gram(columns):
    return columns.conj().T @ columns


def _require_orthonormal(columns, name, atol=1e-8):
    gram = _gram(columns)
    gap = float(np.max(np.abs(gram - np.eye(columns.shape[1]))))
    if gap > atol:
        raise ValueError(f"{name} columns are not orthonormal")
    return gap


def make_holder(grid, *, weights, reference_columns, source_columns, geometry_map=None,
                source_metadata=None, geometry_metadata=None):
    """Build a holder on one full geometry band and the ambient fermion band."""
    if grid.gauge != "conformal":
        raise ValueError("parent response uses the conformal chart L=Q, beta=0")
    regions = regional_map(grid.length)
    weights = _real_vector(weights, int(np.asarray(weights).shape[0]), "weights")
    if weights.shape != (int(weights.shape[0]),):
        raise ValueError("weights must be one vector")
    rank = int(weights.shape[0])
    stacked = _columns(source_columns, (2 * grid.nf, rank), "source")
    reference = _columns(reference_columns, (2 * grid.nf, int(np.shape(reference_columns)[1])), "reference")
    if reference.shape[1] < 1:
        raise ValueError("the fixed observer needs at least one column")
    source_gap = float(np.max(np.abs(_gram(stacked) - np.eye(rank))))
    reference_gap = _require_orthonormal(reference, "modal observer")
    if geometry_map is None:
        geometry_map = np.eye(grid.ng)
    geometry = _real_vector(np.asarray(geometry_map, dtype=float).reshape(-1), grid.ng * grid.ng, "geometry map")
    geometry = geometry.reshape(grid.ng, grid.ng)
    if geometry.shape != (grid.ng, grid.ng):
        raise ValueError("W must be a full square frame so every geometry degree stays active")
    orthogonality = float(np.max(np.abs(geometry.T @ geometry - np.eye(grid.ng))))
    if orthogonality > 1e-8:
        raise ValueError("W is not a full orthogonal geometry frame")
    weighted = _install_weights(grid, weights)
    fft_grid = install_native_fft(weighted)
    fft_grid = _install_weights(fft_grid, weights)
    metadata = {
        "layout": "variable_rank_antiperiodic",
        "rank": rank,
        "observer_rows": int(reference.shape[1]),
        "source_orthonormal": bool(source_gap <= 1e-8),
        "source_may_be_nonorthogonal": True,
        "gram_distance_from_identity": source_gap,
        "reference_gram_max": reference_gap,
        "modal_observer_orthonormal": True,
        "occupations": weights.tolist(),
        "dummy_positive_columns": False,
        "build_pair_called": False,
        "exact_spatial_support": False,
        "physical_vacuum_identification": False,
        "complement_occupation_is_not_a_vacuum_claim": True,
        "fixed_observer_is_spatial_window": False,
    }
    if source_metadata:
        metadata.update(source_metadata)
    geometry_meta = {
        "W_is_identity": bool(np.array_equal(geometry, np.eye(grid.ng))),
        "orthogonality_max": orthogonality,
        "all_geometry_degrees_retained": True,
        "ambient_field_active": True,
        "native_fft": True,
        "fft_workers": 1,
        **regions,
    }
    if geometry_metadata:
        geometry_meta.update(geometry_metadata)
    holder = ParentHolder(
        grid=fft_grid,
        geometry_map=_readonly(geometry),
        weights=_readonly(weights),
        reference_columns=_readonly(reference),
        source_columns=_readonly(stacked),
        source_metadata=metadata,
        geometry_metadata=geometry_meta,
        child_interval=tuple(regions["child_interval"]),
        parent_interval=tuple(regions["parent_interval"]),
        clock_locations=(float(regions["centre"]),),
    )
    _live_weights(holder)
    return holder


def _live_weights(holder):
    live = np.asarray(holder.grid.fine.occupations, dtype=float)
    owned = np.asarray(holder.weights, dtype=float)
    if live.shape != owned.shape or not np.array_equal(live, owned):
        raise ValueError("grid.fine.occupations must be the owned weight vector")
    return owned


def with_weights(holder, weights):
    """Return a holder whose live occupations are the supplied weights."""
    updated = _install_weights(holder.grid, weights)
    return replace(
        holder,
        grid=updated,
        weights=_readonly(updated.fine.occupations),
        source_metadata=dict(holder.source_metadata, occupations=updated.fine.occupations.tolist(), rank=int(updated.fine.occupations.size)),
    )


def validate_tangent(holder, tangent):
    rank = int(_live_weights(holder).shape[0])
    fields = {
        "Q": _real_vector(tangent.Q, holder.grid.ng, "Q"),
        "r": _real_vector(tangent.r, holder.grid.ng, "r"),
        "p_Q": _real_vector(tangent.p_Q, holder.grid.ng, "p_Q"),
        "p_r": _real_vector(tangent.p_r, holder.grid.ng, "p_r"),
        "phi0": _columns(tangent.phi0, (holder.grid.nf, rank), "phi0"),
        "phi1": _columns(tangent.phi1, (holder.grid.nf, rank), "phi1"),
        "occupations": _real_vector(tangent.occupations, rank, "occupations"),
    }
    return ParentTangent(**fields)


def _as_state(values):
    return leading.State(*(np.array(getattr(values, name), copy=True) for name in STATE_FIELDS))


def _check_state(holder, state):
    rank = int(_live_weights(holder).shape[0])
    for name in REAL_FIELDS:
        _real_vector(getattr(state, name), holder.grid.ng, name)
    _columns(state.phi0, (holder.grid.nf, rank), "phi0")
    _columns(state.phi1, (holder.grid.nf, rank), "phi1")
    return _as_state(state)


def _fine_bundle(holder, state):
    fine = leading.check_chart(holder, state)
    system = leading.active_system(holder.grid, fine)
    return fine, system


def _prolonged_tangent(holder, state, tangent):
    fine, system = _fine_bundle(holder, state)
    lifted = leading.prolong(holder.grid, leading.decode(holder, _as_state(tangent)))
    return fine, system, lifted


def _weight_channel_pi(holder, state, tangent):
    """Canonical ``δπ_Q`` from ``δc`` alone, through the existing nodal forces.

    ``leading.jvp`` holds the occupations fixed. The weight row is this extra
    channel. It uses the same force derivative as the coupled response and the
    same encode as the leading rate. Frozen geometry suppresses it with the
    other geometry rates.
    """
    fine, system, _lifted = _prolonged_tangent(holder, state, tangent)
    zero0 = np.zeros_like(fine.phi0)
    zero1 = np.zeros_like(fine.phi1)
    source = coupling.source_from_columns(system, fine)
    kinetic, spin, current = response._owned_moments(source)
    dkinetic, dspin, dcurrent = response._moments_tangent(
        system.momentum, fine.phi0, fine.phi1, zero0, zero1, system.occupations, tangent.occupations,
    )
    delta_l, delta_q, _delta_beta = response._force_directional(
        system, kinetic, dkinetic, spin, dspin, current, dcurrent,
        fine.Q, np.zeros_like(fine.Q), system.length_density, np.zeros_like(fine.Q),
    )
    delta_p = -(delta_q + delta_l) / holder.grid.dx_q
    nodal = galerkin.pull_geometry(holder.grid, delta_p)
    return holder.grid.dx_g * (holder.geometry_map.T @ nodal)


def _force_and_operator(holder, state, tangent):
    fine, system, lifted = _prolonged_tangent(holder, state, tangent)
    source = coupling.source_from_columns(system, fine)
    kinetic, spin, current = response._owned_moments(source)
    dkinetic, dspin, dcurrent = response._moments_tangent(
        system.momentum, fine.phi0, fine.phi1, lifted.phi0, lifted.phi1,
        system.occupations, tangent.occupations,
    )
    delta_l, delta_q, delta_beta = response._force_directional(
        system, kinetic, dkinetic, spin, dspin, current, dcurrent,
        fine.Q, lifted.Q, system.length_density, lifted.Q,
    )
    _image0, _image1, dimage0, dimage1 = response._images_tangent(
        system.momentum, fine.phi0, fine.phi1, lifted.phi0, lifted.phi1,
        system.length_density, lifted.Q, fine.Q, lifted.Q,
        system.shift, np.zeros_like(system.shift), system.kappa,
    )
    operator0 = holder.grid.U_f.conj().T @ dimage0
    operator1 = holder.grid.U_f.conj().T @ dimage1
    return delta_l, delta_q, delta_beta, operator0, operator1


def response_rates(holder, state, tangent, *, control_mode="coupled"):
    """Analytic leading Jacobian-vector product plus the source-weight channel.

    Geometry, momenta, columns, the Dirac operator and the nodal forces are
    differentiated together. The state piece is ``leading.jvp``. The weight
    piece is the existing force derivative at fixed columns. No second
    Hamiltonian is formed and no dense propagator is stored.
    """
    mode = episode.normalize_control_mode(control_mode)
    state = _check_state(holder, state)
    tangent = validate_tangent(holder, tangent)
    with _one_thread():
        image = leading.jvp(holder, state, _as_state(tangent), control_mode=mode)
        weight_pi = _weight_channel_pi(holder, state, tangent)
        if mode != "frozen_geometry":
            image = replace(image, p_Q=np.array(image.p_Q + weight_pi, copy=True))
        delta_l, delta_q, delta_beta, operator0, operator1 = _force_and_operator(holder, state, tangent)
        _rate, bundle = leading.rates(holder, state, return_bundle=True, control_mode=mode)
    return ResponseRate(
        np.array(image.Q, copy=True), np.array(image.r, copy=True),
        np.array(image.p_Q, copy=True), np.array(image.p_r, copy=True),
        np.array(image.phi0, copy=True), np.array(image.phi1, copy=True),
        delta_l, delta_q, delta_beta, operator0, operator1,
        float(_rate.fieldwork_power),
        weight_channel_included=mode != "frozen_geometry",
    )


def _combine_tangent(tangent, rate, scale):
    values = [getattr(tangent, name) + scale * getattr(rate, name) for name in STATE_FIELDS]
    return ParentTangent(*values, np.array(tangent.occupations, copy=True))


def centre_clock(holder, state, tangent=None):
    """``τ̇ = (r Q)`` at ``x = L/2``, and ``δτ̇`` when a tangent is supplied."""
    with _one_thread():
        fine = leading.check_chart(holder, state)
        proper = fine.r * fine.Q
        centre = float(holder.clock_locations[0])
        if abs(centre - float(holder.grid.length) / 2.0) > 1e-12:
            raise ValueError("the matched clock is the carrier centre")
        tau_dot = float(extent.real_periodic_values(holder.grid, proper, (centre,))[0])
        if tangent is None:
            return tau_dot
        tangent = validate_tangent(holder, tangent)
        lifted = leading.prolong(holder.grid, leading.decode(holder, _as_state(tangent)))
        delta = lifted.r * fine.Q + fine.r * lifted.Q
        delta_tau_dot = float(extent.real_periodic_values(holder.grid, delta, (centre,))[0])
    return tau_dot, delta_tau_dot


def centre_clock_variation(delta_observable_t, observable_dot, delta_tau, tau_dot):
    """``δO|τ = δO|t − Ȯ δτ / τ̇`` from the existing matched-clock owner."""
    return response.matched_tau_correction(delta_observable_t, observable_dot, delta_tau, tau_dot)


def _rk4_average(samples):
    values = [float(sample) for sample in samples]
    if len(values) != 4:
        raise ValueError("RK4 clock quadrature needs four stage samples")
    return (values[0] + 2.0 * values[1] + 2.0 * values[2] + values[3]) / 6.0


def advance_tangent(holder, state, tangent, dt, *, control_mode="coupled"):
    """One shared RK4 step of the state and its retarded tangent.

    The state stages are ``leading.rates``. The tangent stages are
    ``response_rates``. Centre-clock samples from those stages are integrated
    with the same weights. Occupation directions are parameters: their rate is
    zero. This is not a held future measurement.
    """
    step = float(dt)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("a positive finite timestep is required")
    mode = episode.normalize_control_mode(control_mode)
    state = _check_state(holder, state)
    tangent = validate_tangent(holder, tangent)
    stages = [(state, tangent)]
    state_rates = []
    tangent_rates = []
    clocks = []
    with _one_thread():
        for factor in (0.5, 0.5, 1.0, None):
            current, direction = stages[-1]
            state_rate = leading.rates(holder, current, control_mode=mode)
            tangent_rate = response_rates(holder, current, direction, control_mode=mode)
            state_rates.append(state_rate)
            tangent_rates.append(tangent_rate)
            clocks.append(centre_clock(holder, current, direction))
            if factor is not None:
                stages.append((
                    leading.combine(state, state_rate, step * factor),
                    _combine_tangent(tangent, tangent_rate, step * factor),
                ))
        new_state = leading.State(*(
            getattr(state, name) + step * (
                getattr(state_rates[0], name) + 2.0 * getattr(state_rates[1], name)
                + 2.0 * getattr(state_rates[2], name) + getattr(state_rates[3], name)
            ) / 6.0
            for name in STATE_FIELDS
        ))
        if mode == "frozen_geometry":
            for name in REAL_FIELDS:
                setattr(new_state, name, np.array(getattr(state, name), copy=True))
        new_tangent = ParentTangent(*(
            getattr(tangent, name) + step * (
                getattr(tangent_rates[0], name) + 2.0 * getattr(tangent_rates[1], name)
                + 2.0 * getattr(tangent_rates[2], name) + getattr(tangent_rates[3], name)
            ) / 6.0
            for name in STATE_FIELDS
        ), np.array(tangent.occupations, copy=True))
        leading.check_chart(holder, new_state)
    tau_dots = [sample[0] for sample in clocks]
    delta_dots = [sample[1] for sample in clocks]
    report = {
        "dt": step,
        "tau": step * _rk4_average(tau_dots),
        "delta_tau": step * _rk4_average(delta_dots),
        "stage_tau_dot": tau_dots,
        "stage_delta_tau_dot": delta_dots,
        "clock_quadrature": "RK4 stage weights",
        "occupation_rate": "parameter direction held; no occupation ODE",
        "dense_propagator_stored": False,
        "held_nonlinear_measurement": False,
    }
    return new_state, new_tangent, report


def _interval_integral(grid, values, interval):
    return float(extent.real_interval_integral(grid, values, interval))


def _probability_density(holder, fine):
    mass = np.sum((np.abs(fine.phi0) ** 2 + np.abs(fine.phi1) ** 2) * _live_weights(holder), axis=1)
    return mass / holder.grid.dx_q


def child_regional_content(holder, state, tangent=None):
    """Probability in the child window ``|s|<=1/2``. Not the modal kernel."""
    with _one_thread():
        fine = leading.check_chart(holder, state)
        density = _probability_density(holder, fine)
        content = _interval_integral(holder.grid, density, holder.child_interval)
        if tangent is None:
            return content
        tangent = validate_tangent(holder, tangent)
        lifted = leading.prolong(holder.grid, leading.decode(holder, _as_state(tangent)))
        delta_density = response._density_tangent(
            fine.phi0, fine.phi1, lifted.phi0, lifted.phi1,
            _live_weights(holder), tangent.occupations, holder.grid.dx_q,
        )
        delta = _interval_integral(holder.grid, delta_density, holder.child_interval)
    return content, delta


def proper_localization(holder, state, tangent=None):
    """Child proper length and the child share of the carrier probability."""
    with _one_thread():
        fine = leading.check_chart(holder, state)
        proper = fine.r * fine.Q
        density = _probability_density(holder, fine)
        child_proper = _interval_integral(holder.grid, proper, holder.child_interval)
        carrier_probability = _interval_integral(holder.grid, density, (0.0, float(holder.grid.length)))
        child_probability = _interval_integral(holder.grid, density, holder.child_interval)
        collar = holder.geometry_metadata["protected_collar_interval"]
        annulus = holder.geometry_metadata["parent_annulus_intervals"]
        collar_probability = _interval_integral(holder.grid, density, collar)
        annulus_probability = sum(_interval_integral(holder.grid, density, interval) for interval in annulus)
        if tangent is None:
            return {
                "child_proper_length": child_proper,
                "child_probability": child_probability,
                "carrier_probability": carrier_probability,
                "child_probability_fraction": child_probability / carrier_probability,
                "protected_collar_probability": collar_probability,
                "parent_annulus_probability": annulus_probability,
                "spatial_window": True,
                "fixed_mode_observer": False,
            }
        tangent = validate_tangent(holder, tangent)
        lifted = leading.prolong(holder.grid, leading.decode(holder, _as_state(tangent)))
        delta_proper_density = lifted.r * fine.Q + fine.r * lifted.Q
        delta_proper = _interval_integral(holder.grid, delta_proper_density, holder.child_interval)
    return {
        "child_proper_length": child_proper,
        "delta_child_proper_length": delta_proper,
        "child_probability": child_probability,
        "spatial_window": True,
        "fixed_mode_observer": False,
    }


def fixed_mode_observer(holder, state=None):
    """Thin modal kernel of the fixed reference columns.

    ``K = (V† Φ) diag(c) (Φ† V)``. This is not the child spatial integral, and
    the dense one-particle covariance is not formed.
    """
    columns = holder.source_columns if state is None else np.vstack((state.phi0, state.phi1))
    if columns.shape != holder.source_columns.shape:
        raise ValueError("the modal observer uses the live source rank")
    amplitudes = holder.reference_columns.conj().T @ columns
    weights = _live_weights(holder)
    kernel = (amplitudes * weights[None, :]) @ amplitudes.conj().T
    return {
        "kernel": kernel,
        "rows": int(holder.reference_columns.shape[1]),
        "spatial_window": False,
        "dense_covariance_formed": False,
        "fixed_at_construction": state is None,
    }


def velocity_tangent(holder, state, *, control_mode="coupled"):
    """Rate direction with a zero occupation increment."""
    rate = leading.rates(holder, state, control_mode=control_mode)
    return ParentTangent(*(getattr(rate, name) for name in STATE_FIELDS), np.zeros(_live_weights(holder).shape))


def matched_centre_readout(holder, state, tangent, delta_tau, *, control_mode="coupled"):
    """Primary content at fixed coordinate time and at matched centre proper time."""
    content, delta_content = child_regional_content(holder, state, tangent)
    _velocity_content, content_dot = child_regional_content(
        holder, state, velocity_tangent(holder, state, control_mode=control_mode),
    )
    tau_dot = centre_clock(holder, state)
    return {
        "primary_readout": "child_regional_content",
        "child_regional_content": content,
        "delta_child_regional_content_t": delta_content,
        "child_regional_content_dot": content_dot,
        "delta_tau": float(delta_tau),
        "tau_dot": tau_dot,
        "delta_child_regional_content_tau": centre_clock_variation(
            delta_content, content_dot, delta_tau, tau_dot,
        ),
        "clock": "centre",
        "spatial_window_is_fixed_mode_observer": False,
    }


def _nonlinear_centre_increment(holder, state, dt, *, control_mode="coupled"):
    mode = episode.normalize_control_mode(control_mode)
    clocks = []
    current = state
    for factor in (0.5, 0.5, 1.0, None):
        clocks.append(centre_clock(holder, current))
        if factor is not None:
            current = leading.combine(state, leading.rates(holder, current, control_mode=mode), dt * factor)
    return dt * _rk4_average(clocks)


def nonlinear_endpoint_clock_bracket(holder, state, tangent, dt, *, step, control_mode="coupled"):
    """Matched centre bracket from the nonlinear endpoints, not from ``δτ`` of the tangent.

    The base and the perturbed state are each advanced by ``leading.rk4_step``.
    ``δO`` and ``δτ`` are the centred difference of those endpoints. ``Ȯ`` and
    ``τ̇`` are read on the base endpoint. The integrated linear clock is not an
    input.
    """
    probe = float(step)
    if not np.isfinite(probe) or probe == 0.0:
        raise ValueError("the nonlinear clock bracket needs a nonzero probe")
    state = _check_state(holder, state)
    tangent = validate_tangent(holder, tangent)
    with _one_thread():
        base = leading.rk4_step(holder, state, dt, control_mode=control_mode)
        perturbed_holder = with_weights(holder, _live_weights(holder) + probe * tangent.occupations)
        perturbed = leading.rk4_step(
            perturbed_holder, leading.combine(state, _as_state(tangent), probe), dt,
            control_mode=control_mode,
        )
        down_holder = with_weights(holder, _live_weights(holder) - probe * tangent.occupations)
        down = leading.rk4_step(
            down_holder, leading.combine(state, _as_state(tangent), -probe), dt,
            control_mode=control_mode,
        )
        delta_content = (
            child_regional_content(perturbed_holder, perturbed) - child_regional_content(down_holder, down)
        ) / (2.0 * probe)
        delta_tau = (
            _nonlinear_centre_increment(perturbed_holder, leading.combine(state, _as_state(tangent), probe), dt, control_mode=control_mode)
            - _nonlinear_centre_increment(down_holder, leading.combine(state, _as_state(tangent), -probe), dt, control_mode=control_mode)
        ) / (2.0 * probe)
        _content, content_dot = child_regional_content(holder, base, velocity_tangent(holder, base, control_mode=control_mode))
        tau_dot = centre_clock(holder, base)
        bracket = centre_clock_variation(delta_content, content_dot, delta_tau, tau_dot)
    return {
        "delta_child_regional_content_t": float(delta_content),
        "delta_tau": float(delta_tau),
        "tau_dot": float(tau_dot),
        "child_regional_content_dot": float(content_dot),
        "delta_child_regional_content_tau": float(bracket),
        "formula": "delta_O_tau = delta_O_t - Odot delta_tau / tau_dot",
        "source": "nonlinear_endpoint",
        "uses_integrated_linear_delta_tau": False,
        "probe": probe,
        "held_nonlinear_campaign": False,
    }


def explicit_normalization_derivative(phi0, phi1, dphi0, dphi1):
    """Analytic Hermitian projection that keeps an orthonormal frame.

    For ``Q†Q = I`` the unitary tangent is ``δQ = δM − Q herm(Q† δM)``.
    A pure Hermitian gauge mode is removed. A direction orthogonal to the
    frame is retained. This is not a finite-difference substitute.
    """
    stacked = np.vstack((np.asarray(phi0, dtype=complex), np.asarray(phi1, dtype=complex)))
    delta = np.vstack((np.asarray(dphi0, dtype=complex), np.asarray(dphi1, dtype=complex)))
    if stacked.shape != delta.shape or stacked.ndim != 2:
        raise ValueError("normalization derivative needs matching column shapes")
    gap = _require_orthonormal(stacked, "normalization")
    generator = stacked.conj().T @ delta
    hermitian = 0.5 * (generator + generator.conj().T)
    corrected = delta - stacked @ hermitian
    nf = int(np.asarray(phi0).shape[0])
    return corrected[:nf], corrected[nf:], {
        "method": "analytic_hermitian_polar_projection",
        "approximation": False,
        "gram_max": gap,
        "state_derivative_retained": True,
        "hermitian_gauge_removed": True,
    }


def projected_hamilton_constraint(holder, state):
    """Projected ``Q C`` residual from the leading constraint arrays."""
    with _one_thread():
        fine, system = _fine_bundle(holder, state)
        source = coupling.source_from_columns(system, fine)
        density, _shift = leading.constraint_arrays(holder.grid, fine, system, source)
        residual = galerkin.pull_geometry(holder.grid, fine.Q * density)
    return np.array(residual, dtype=float, copy=True)


def _radius_jacobian(holder, state, step):
    size = int(holder.grid.ng)
    columns = []
    for index in range(size):
        direction = np.zeros(size)
        direction[index] = 1.0
        plus = state.copy()
        minus = state.copy()
        plus.r = state.r + step * direction
        minus.r = state.r - step * direction
        columns.append((projected_hamilton_constraint(holder, plus) - projected_hamilton_constraint(holder, minus)) / (2.0 * step))
    return np.column_stack(columns)


def _prepared_probe(holder, state, tangent, scale):
    trial = leading.combine(state, _as_state(replace_radius(tangent)), scale)
    trial.r = np.array(state.r, copy=True)
    probed = with_weights(holder, _live_weights(holder) + scale * tangent.occupations)
    return projected_hamilton_constraint(probed, trial)


def preparation_observables(holder, state):
    """Constraint, momentum constraint, canonical momenta, and column norms.

    Source columns are not forced onto an orthonormal frame. The modal
    observer is a separate object and is not this norm.
    """
    fine, system = _fine_bundle(holder, state)
    source = coupling.source_from_columns(system, fine)
    hamilton, shift = leading.constraint_arrays(holder.grid, fine, system, source)
    stacked = np.vstack((state.phi0, state.phi1))
    return {
        "constraint": np.array(galerkin.pull_geometry(holder.grid, fine.Q * hamilton), copy=True),
        "momentum_constraint": np.array(galerkin.pull_geometry(holder.grid, shift), copy=True),
        "momenta": np.concatenate((np.asarray(state.p_Q, dtype=float), np.asarray(state.p_r, dtype=float))),
        "global_norm": np.linalg.norm(stacked, axis=0),
    }


def _full_preparation_probe(holder, state, tangent, scale):
    trial = leading.combine(state, _as_state(tangent), scale)
    probed = with_weights(holder, _live_weights(holder) + scale * tangent.occupations)
    return preparation_observables(probed, trial)


def _centred_preparation(holder, state, tangent, width):
    plus = _full_preparation_probe(holder, state, tangent, width)
    minus = _full_preparation_probe(holder, state, tangent, -width)
    return {key: (plus[key] - minus[key]) / (2.0 * width) for key in plus}


def replace_radius(tangent, radius=None):
    """Copy a tangent with its radius slot replaced. ``None`` zeroes that slot."""
    radius = np.zeros_like(tangent.r) if radius is None else np.array(radius, dtype=float, copy=True)
    return ParentTangent(tangent.Q, radius, tangent.p_Q, tangent.p_r, tangent.phi0, tangent.phi1, tangent.occupations)


def preparation_tangent(holder, state, tangent, *, step, allow_dense_radius_jacobian=False):
    """Two-step finite difference of the full T0 preparation, not only ``δr``.

    The centred probes move the state, the momenta, the columns and the
    weights. They record the Hamilton constraint, the shift constraint, the
    canonical momenta and the actual column norms. A radius solve remains one
    component. It is not the preparation. The difference is not labelled analytic.
    """
    probe = float(step)
    if not np.isfinite(probe) or probe <= 0.0:
        raise ValueError("preparation probe must be positive")
    if holder.grid.ng > DENSE_RADIUS_JACOBIAN_MAX and not allow_dense_radius_jacobian:
        raise RuntimeError(
            "analytic preparation Jacobian is absent; refusing a dense finite-difference "
            f"radius Jacobian above ng={DENSE_RADIUS_JACOBIAN_MAX}"
        )
    state = _check_state(holder, state)
    tangent = validate_tangent(holder, tangent)
    stacked = np.vstack((state.phi0, state.phi1))
    source_gap = float(np.max(np.abs(_gram(stacked) - np.eye(stacked.shape[1]))))
    if source_gap <= 1e-8:
        dphi0, dphi1, normalization = explicit_normalization_derivative(
            state.phi0, state.phi1, tangent.phi0, tangent.phi1,
        )
    else:
        dphi0 = np.array(tangent.phi0, copy=True)
        dphi1 = np.array(tangent.phi1, copy=True)
        normalization = {
            "method": "source_not_orthonormal",
            "approximation": False,
            "global_norm_in_finite_difference": True,
            "gram_distance_from_identity": source_gap,
        }
    kept = ParentTangent(tangent.Q, tangent.r, tangent.p_Q, tangent.p_r, dphi0, dphi1, tangent.occupations)
    with _one_thread():
        residual_scale = float(np.max(np.abs(projected_hamilton_constraint(holder, state))))
        samples = {}
        solved = {}
        full = {}
        for factor, label in ((1.0, "h"), (0.5, "h_over_2")):
            width = probe * factor
            residual = (_prepared_probe(holder, state, kept, width) - _prepared_probe(holder, state, kept, -width)) / (2.0 * width)
            jacobian = _radius_jacobian(holder, state, width)
            full[label] = _centred_preparation(holder, state, kept, width)
            try:
                delta_r = np.linalg.solve(jacobian, -residual)
            except np.linalg.LinAlgError:
                return {
                    "available": False,
                    "tangent": None,
                    "method": PREPARATION_METHOD,
                    "analytic_preparation_jacobian": False,
                    "finite_difference_used_as_analytic": False,
                    "radius_only": False,
                    "missing_primitive": "singular_finite_difference_radius_jacobian",
                    "normalization": normalization,
                    "full_preparation": {key: full[label][key].tolist() for key in full[label]},
                    "coordinate_time": 0.0,
                    "dropped_state_derivative": False,
                }
            samples[label] = residual
            solved[label] = delta_r
    control = float(np.max(np.abs(solved["h"] - solved["h_over_2"])))
    components = {}
    for key in full["h"]:
        components[key] = {
            "two_step_gap": float(np.max(np.abs(full["h"][key] - full["h_over_2"][key]))),
            "max_abs": float(np.max(np.abs(full["h_over_2"][key]))),
        }
    output = ParentTangent(kept.Q, solved["h_over_2"], kept.p_Q, kept.p_r, kept.phi0, kept.phi1, kept.occupations)
    return {
        "available": True,
        "tangent": output,
        "method": PREPARATION_METHOD,
        "analytic_preparation_jacobian": False,
        "finite_difference_used_as_analytic": False,
        "radius_only": False,
        "full_preparation_components": ("constraint", "momentum_constraint", "momenta", "global_norm"),
        "components": components,
        "normalization": normalization,
        "steps": [probe, probe / 2.0],
        "two_step_radius_gap": control,
        "two_step_residual_gap": float(np.max(np.abs(samples["h"] - samples["h_over_2"]))),
        "two_step_momenta_gap": components["momenta"]["two_step_gap"],
        "two_step_global_norm_gap": components["global_norm"]["two_step_gap"],
        "two_step_constraint_gap": components["constraint"]["two_step_gap"],
        "two_step_momentum_constraint_gap": components["momentum_constraint"]["two_step_gap"],
        "constraint_residual_max": residual_scale,
        "constraint_already_solved": bool(residual_scale <= 1e-8),
        "radius_component": "implicit_linear_solve_plus_full_preparation_difference",
        "other_state_derivatives": "retained",
        "dropped_state_derivative": False,
        "input_radius_tangent_max": float(np.max(np.abs(tangent.r))),
        "coordinate_time": 0.0,
        "held_nonlinear_measurement": False,
    }


def _widths(widths):
    parsed = tuple(int(value) for value in widths)
    if len(parsed) != 3 or min(parsed) < 1:
        raise ValueError("hierarchy widths are three positive row counts")
    return parsed


def _blocks(basis, widths):
    widths = _widths(widths)
    if basis.shape[1] != sum(widths):
        raise ValueError("basis width must equal the grandchild, detail and parent rows")
    grandchild_rows, detail_rows, parent_rows = widths
    return (
        basis[:, :grandchild_rows],
        basis[:, grandchild_rows:grandchild_rows + detail_rows],
        basis[:, grandchild_rows + detail_rows:],
        widths,
    )


def initial_hierarchy(field, basis, widths):
    """Direct and sequential openings. Inherited routes keep the saved widths."""
    field = np.asarray(field, dtype=complex)
    basis = np.asarray(basis, dtype=complex)
    parsed = _widths(widths)
    if parsed == INHERITED_WIDTHS and basis.shape[1] == sum(INHERITED_WIDTHS):
        return grandchild.initial_routes(field, basis)
    grandchild_rows, detail_rows, parent_rows, parsed = _blocks(basis, parsed)
    direct = {
        "a": grandchild_rows.conj().T @ field,
        "drive": field - grandchild_rows @ (grandchild_rows.conj().T @ field),
        "memory": np.zeros_like(field),
    }
    outer = field - basis @ (basis.conj().T @ field)
    sequential = {
        "a": grandchild_rows.conj().T @ field,
        "d_drive": detail_rows.conj().T @ field,
        "d_memory": np.zeros((detail_rows.shape[1], field.shape[1]), dtype=complex),
        "f_drive": parent_rows.conj().T @ field,
        "f_memory": np.zeros((parent_rows.shape[1], field.shape[1]), dtype=complex),
        "drive": outer,
        "memory": np.zeros_like(field),
    }
    return direct, sequential


def reconstruct_hierarchy(route, basis, widths, sequential=False):
    parsed = _widths(widths)
    if parsed == INHERITED_WIDTHS and basis.shape[1] == sum(INHERITED_WIDTHS):
        return grandchild.reconstruct_route(route, basis, sequential=sequential)
    grandchild_rows, detail_rows, parent_rows, _parsed = _blocks(basis, parsed)
    result = grandchild_rows @ route["a"] + route["drive"] + route["memory"]
    if sequential:
        result = (
            result
            + detail_rows @ (route["d_drive"] + route["d_memory"])
            + parent_rows @ (route["f_drive"] + route["f_memory"])
        )
    return result


def extended_hierarchy_rates(image, basis, direct, sequential, widths, *, omit_outer_detail=False):
    """Variable-row form of the inherited streaming split. One image call."""
    grandchild_rows, detail_rows, parent_rows, _parsed = _blocks(basis, widths)
    parts = [
        grandchild_rows @ direct["a"], direct["drive"], direct["memory"],
        grandchild_rows @ sequential["a"],
        detail_rows @ sequential["d_drive"], detail_rows @ sequential["d_memory"],
        parent_rows @ sequential["f_drive"], parent_rows @ sequential["f_memory"],
        sequential["drive"], sequential["memory"],
    ]
    if len({part.shape[1] for part in parts}) != 1:
        raise ValueError("streamed columns must share one source rank")
    images = np.split(image(np.concatenate(parts, axis=1)), 10, axis=1)
    retained, drive, memory, sequential_retained, detail_drive, detail_memory, parent_drive, parent_memory, outer_drive, outer_memory = images

    def complement_grandchild(value):
        return value - grandchild_rows @ (grandchild_rows.conj().T @ value)

    def complement_basis(value):
        return value - basis @ (basis.conj().T @ value)

    first = {
        "a": -1j * grandchild_rows.conj().T @ (retained + drive + memory),
        "drive": -1j * complement_grandchild(drive),
        "memory": -1j * complement_grandchild(memory + retained),
    }
    outside = parent_drive + parent_memory + outer_drive + outer_memory
    detail_source = detail_memory + sequential_retained + (0.0 if omit_outer_detail else outside)
    second = {
        "a": -1j * grandchild_rows.conj().T @ (sequential_retained + detail_drive + detail_memory + outside),
        "d_drive": -1j * detail_rows.conj().T @ detail_drive,
        "d_memory": -1j * detail_rows.conj().T @ detail_source,
        "f_drive": -1j * parent_rows.conj().T @ parent_drive,
        "f_memory": -1j * parent_rows.conj().T @ (
            parent_memory + sequential_retained + detail_drive + detail_memory + outer_drive + outer_memory
        ),
        "drive": -1j * complement_basis(outer_drive),
        "memory": -1j * complement_basis(
            outer_memory + sequential_retained + detail_drive + detail_memory + parent_drive + parent_memory
        ),
    }
    return first, second


def hierarchy_rates(image, basis, direct, sequential, widths, *, omit_outer_detail=False):
    """Stream the direct and sequential splits. No dense exterior propagator."""
    parsed = _widths(widths)
    if parsed == INHERITED_WIDTHS and np.asarray(basis).shape[1] == sum(INHERITED_WIDTHS):
        return grandchild.route_rates(image, basis, direct, sequential, omit_outer_detail=omit_outer_detail)
    return extended_hierarchy_rates(
        image, basis, direct, sequential, parsed, omit_outer_detail=omit_outer_detail,
    )


def column_image(holder, state, columns):
    """One matrix-free Dirac image on the band. The propagator is not stored."""
    columns = np.asarray(columns, dtype=complex)
    if columns.ndim != 2 or columns.shape[0] != 2 * holder.grid.nf:
        raise ValueError("streamed columns must span the fermion band")
    with _one_thread():
        fine, system = _fine_bundle(holder, state)
        nf = holder.grid.nf
        image0, image1 = coupling.apply_dirac(
            holder.grid.U_f @ columns[:nf], holder.grid.U_f @ columns[nf:],
            system.length_density, fine.Q, system.shift, system.kappa, system.momentum,
        )
    return np.vstack((holder.grid.U_f.conj().T @ image0, holder.grid.U_f.conj().T @ image1))


def _thin_cross(amplitudes, exterior, weights):
    return (np.asarray(amplitudes) * np.asarray(weights)[None, :]) @ np.asarray(exterior).conj().T


def split_state(holder, state, basis, widths, *, require_correlation=True):
    """One metric, variable retained rows, and a prepared nonzero cross."""
    state = _check_state(holder, state)
    basis = np.asarray(basis, dtype=complex)
    parsed = _widths(widths)
    if basis.shape != (2 * holder.grid.nf, sum(parsed)):
        raise ValueError("hierarchy basis must cover the fermion band")
    _require_orthonormal(basis, "hierarchy")
    field = np.vstack((state.phi0, state.phi1))
    direct, sequential = initial_hierarchy(field, basis, parsed)
    cross = _thin_cross(direct["a"], direct["drive"], _live_weights(holder))
    correlation = float(np.linalg.norm(cross))
    if require_correlation and correlation <= 1e-12:
        raise ValueError("state correlation control is zero; refusing a product-state placeholder")
    return {
        "state": state,
        "basis": basis,
        "widths": parsed,
        "direct": direct,
        "sequential": sequential,
        "correlation_frobenius": correlation,
        "one_physical_metric": True,
        "retained_rows": int(parsed[0]),
        "source_rank": int(_live_weights(holder).shape[0]),
        "rank6_guard_called": False,
    }


def _forces_of(holder, state, columns):
    nf = holder.grid.nf
    trial = replace(state, phi0=np.array(columns[:nf], copy=True), phi1=np.array(columns[nf:], copy=True))
    fine, system = _fine_bundle(holder, trial)
    source = coupling.source_from_columns(system, fine)
    return source


def _max_abs(left, right):
    return float(np.max(np.abs(np.asarray(left) - np.asarray(right))))


def _route_advance(route, rate, scale):
    return {key: value + scale * rate[key] for key, value in route.items()}


def full_vs_reduced(holder, state, reduced, dt, *, omit_outer_detail=False):
    """One full RK4 step against the streamed direct and sequential reductions.

    Geometry rates come from ``leading.rates`` on the reconstructed source.
    Column rates come from one concatenated Dirac image per stage. Complement
    forces are the existing nodal forces of the exterior pieces. Their cross
    with the retained piece is kept. Memory omission is reported twice:
    conditional on the exterior drive, and as an autonomous retained closure.
    Neither omission is a physical law.
    """
    step = float(dt)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("a positive finite timestep is required")
    basis = reduced["basis"]
    widths = reduced["widths"]
    direct = reduced["direct"]
    sequential = reduced["sequential"]
    opening = np.vstack((state.phi0, state.phi1))
    calls = {"stage": 0, "diagnostic": 0}

    def streamed(columns, current, slot):
        calls[slot] += 1
        return column_image(holder, current, columns)

    with _one_thread():
        direct_gap = _max_abs(reconstruct_hierarchy(direct, basis, widths), opening)
        sequential_gap = _max_abs(
            reconstruct_hierarchy(sequential, basis, widths, sequential=True), opening,
        )
        full_after = leading.rk4_step(holder, state, step)
        stages = [state]
        directs = [direct]
        sequences = [sequential]
        direct_rates = []
        geometry_rates = []
        generator_gap = 0.0
        sequential_gap_max = 0.0
        for factor in (0.5, 0.5, 1.0, None):
            current = stages[-1]
            d_rate, s_rate = hierarchy_rates(
                lambda columns: streamed(columns, current, "stage"),
                basis, directs[-1], sequences[-1], widths, omit_outer_detail=omit_outer_detail,
            )
            full_rate = leading.rates(holder, current)
            full_field = np.vstack((full_rate.phi0, full_rate.phi1))
            generator_gap = max(
                generator_gap,
                _max_abs(reconstruct_hierarchy(d_rate, basis, widths), full_field),
            )
            sequential_gap_max = max(
                sequential_gap_max,
                _max_abs(reconstruct_hierarchy(s_rate, basis, widths, sequential=True), full_field),
            )
            direct_rates.append(d_rate)
            geometry_rates.append(full_rate)
            if factor is not None:
                staged = leading.combine(state, full_rate, step * factor)
                directs.append(_route_advance(direct, d_rate, step * factor))
                sequences.append(_route_advance(sequential, s_rate, step * factor))
                rebuilt = reconstruct_hierarchy(directs[-1], basis, widths)
                staged.phi0 = np.array(rebuilt[:holder.grid.nf], copy=True)
                staged.phi1 = np.array(rebuilt[holder.grid.nf:], copy=True)
                stages.append(staged)
        averaged = {
            key: (direct_rates[0][key] + 2.0 * direct_rates[1][key]
                  + 2.0 * direct_rates[2][key] + direct_rates[3][key]) / 6.0
            for key in direct_rates[0]
        }
        reduced_field = reconstruct_hierarchy(_route_advance(direct, averaged, step), basis, widths)
        opening_image = streamed(opening, state, "diagnostic")
        d_open, _s_open = hierarchy_rates(
            lambda columns: streamed(columns, state, "diagnostic"),
            basis, direct, sequential, widths, omit_outer_detail=omit_outer_detail,
        )
        full_dot = -1j * opening_image
        conditional_rate = dict(d_open)
        conditional_rate["memory"] = np.zeros_like(d_open["memory"])
        conditional_dot = reconstruct_hierarchy(conditional_rate, basis, widths)
        grandchild_rows = basis[:, :widths[0]]
        retained = grandchild_rows @ direct["a"]
        autonomous_amplitudes = -1j * grandchild_rows.conj().T @ streamed(retained, state, "diagnostic")
        autonomous_dot = grandchild_rows @ autonomous_amplitudes
        pieces = {
            "retained": grandchild_rows @ direct["a"],
            "drive": direct["drive"],
            "memory": direct["memory"],
        }
        full_force = _forces_of(holder, state, opening)
        partial = {name: _forces_of(holder, state, columns) for name, columns in pieces.items()}
        complement = _forces_of(holder, state, pieces["drive"] + pieces["memory"])
    force_report = {}
    for name in ("force_L", "force_Q", "force_beta"):
        dropped = partial["retained"][name] + partial["drive"][name] + partial["memory"][name]
        cross = full_force[name] - dropped
        force_report[name] = {
            "cross_max_abs": float(np.max(np.abs(cross))),
            "complement_max_abs": float(np.max(np.abs(complement[name]))),
            "full_max_abs": float(np.max(np.abs(full_force[name]))),
            "geometry_uses_full_source": True,
        }
    def butcher(values):
        return (values[0] + 2.0 * values[1] + 2.0 * values[2] + values[3]) / 6.0

    geometry_gap = max(
        _max_abs(
            getattr(full_after, name),
            getattr(state, name) + step * butcher([getattr(rate, name) for rate in geometry_rates]),
        )
        for name in REAL_FIELDS
    )
    return {
        "dt": step,
        "one_physical_metric": True,
        "retained_rows": int(widths[0]),
        "source_rank": int(state.phi0.shape[1]),
        "opening_direct_gap": direct_gap,
        "opening_sequential_gap": sequential_gap,
        "generator_gap": generator_gap,
        "sequential_generator_gap": sequential_gap_max,
        "endpoint_field_gap": _max_abs(np.vstack((full_after.phi0, full_after.phi1)), reduced_field),
        "endpoint_geometry_gap": geometry_gap,
        "forces": force_report,
        "correlation_frobenius": float(reduced["correlation_frobenius"]),
        "omissions": {
            "conditional": {
                "name": "conditional_memory_omission",
                "autonomous": False,
                "exterior_drive_retained": True,
                "physical_law": False,
                "unitary_global_law": False,
                "claimed_regeneration": False,
                "residual_max_abs": _max_abs(conditional_dot, full_dot),
            },
            "autonomous": {
                "name": "autonomous_retained_closure",
                "autonomous": True,
                "exterior_drive_retained": False,
                "physical_law": False,
                "unitary_global_law": False,
                "claimed_regeneration": False,
                "residual_max_abs": _max_abs(autonomous_dot, full_dot),
            },
        },
        "sequential_omit_outer_detail": {
            "enabled": bool(omit_outer_detail),
            "meaning": "inherited grandchild sequential control; not an autonomous closure and not a rank-6 guard",
        },
        "streaming_calls": int(calls["stage"]),
        "rk4_streaming_calls": int(calls["stage"]),
        "dense_propagator_stored": False,
        "rank6_guard_called": False,
        "memory_owner_guard_not_used": memory.ImplementationGap.__name__,
    }


def _orthonormal_frame(rows, columns, rng):
    raw = rng.normal(size=(rows, columns)) + 1j * rng.normal(size=(rows, columns))
    frame, _r = np.linalg.qr(raw)
    return frame[:, :columns]


def _extend_basis(seed, total, rng):
    seed = np.asarray(seed, dtype=complex)
    if total < seed.shape[1]:
        raise ValueError("basis cannot be shorter than the fixed observer")
    if total == seed.shape[1]:
        return np.array(seed, copy=True)
    raw = _orthonormal_frame(seed.shape[0], total, rng)
    complement = raw - seed @ (seed.conj().T @ raw)
    factor, _r = np.linalg.qr(complement)
    basis = np.concatenate((seed, factor[:, :total - seed.shape[1]]), axis=1)
    _require_orthonormal(basis, "manufactured hierarchy")
    return basis


def manufactured_fixture(*, nf=16, quadrature=32, rank=3, observer_rows=2, widths=(2, 1, 1), seed=11):
    """Small positive chart for tests. Not a physical held run and not a saved episode."""
    if int(rank) == 6:
        raise ValueError("the manufactured fixture does not rebuild the historical six columns")
    parsed = _widths(widths)
    if parsed[0] != int(observer_rows):
        raise ValueError("the fixed observer is the first hierarchy block")
    grid = galerkin.build_grid(int(nf), quadrature=int(quadrature), length=CARRIER_LENGTH, gauge="conformal")
    rng = np.random.default_rng(int(seed))
    source = _orthonormal_frame(2 * grid.nf, int(rank), rng)
    reference = _orthonormal_frame(2 * grid.nf, int(observer_rows), rng)
    weights = np.linspace(0.2, 0.8, int(rank))
    holder = make_holder(grid, weights=weights, reference_columns=reference, source_columns=source,
                         source_metadata={"kind": "manufactured_fixture", "physical_episode": False})
    ratio = float(grid.fine.calibration["b0"] / grid.fine.calibration["a0"])
    nodal = leading.State(
        np.full(grid.ng, ratio), np.ones(grid.ng), np.zeros(grid.ng), np.zeros(grid.ng),
        source[:grid.nf], source[grid.nf:],
    )
    state = leading.encode(holder, nodal)
    angle = 2.0 * np.pi * grid.xi_g / grid.length
    nodal_tangent = leading.State(
        1e-3 * np.cos(angle), 1e-3 * np.sin(angle),
        1e-4 * np.cos(2.0 * angle), -1e-4 * np.sin(2.0 * angle),
        1e-3 * _orthonormal_frame(grid.nf, int(rank), rng),
        1e-3 * _orthonormal_frame(grid.nf, int(rank), rng),
    )
    encoded = leading.encode(holder, nodal_tangent)
    tangent = ParentTangent(
        *(getattr(encoded, name) for name in STATE_FIELDS),
        np.linspace(-0.02, 0.02, int(rank)),
    )
    basis = _extend_basis(reference, sum(parsed), rng)
    binding = {
        "available": True,
        "kind": "manufactured_fixture",
        "physical_episode": False,
        "coordinate_time": 0.0,
        "laboratory_dt": 1.0e-4,
        "holder": holder,
        "state": state,
        "tangent": tangent,
        "basis": basis,
        "widths": parsed,
        "episode_path": None,
    }
    return binding


def _contains(path, root):
    return path == root or root in path.parents


def _sealed_or_input(path, input_paths):
    if any(part in SEALED_DIRECTORY_NAMES for part in path.parts):
        return True
    for raw in input_paths:
        if raw is None:
            continue
        source = Path(raw).expanduser().resolve()
        if path == source:
            return True
        if source.is_dir() and _contains(path, source):
            return True
        if source.is_file() and _contains(source, path):
            return True
    return False


def authorize_stage_output(directory, *, frozen_producer, physical_binding, closure_matches, input_paths=()):
    """Allow a temp stage, or a new development directory after freeze.

    Input and sealed directories are never writable. A development path
    without a frozen producer or a physical binding stays open. It is not a
    blanket refusal of ``lab/results/development``.
    """
    resolved = Path(directory).expanduser().resolve()
    if _sealed_or_input(resolved, input_paths):
        raise PermissionError("refusing to write into an input or sealed directory")
    if not _contains(resolved, REPO):
        return resolved
    if not _contains(resolved, DEVELOPMENT_ROOT.resolve()) or resolved == DEVELOPMENT_ROOT.resolve():
        raise PermissionError("repository stages must be new directories under lab/results/development")
    if not frozen_producer:
        raise AuthorizationOpen("development output stays open until the producer is frozen")
    if not physical_binding:
        raise AuthorizationOpen("development output stays open until a physical source or episode binding exists")
    if not closure_matches:
        raise ValueError("source closure does not match the current producers")
    return resolved


def _directory_size(path):
    total = 0
    if path.exists():
        for child in path.rglob("*"):
            if child.is_file():
                total += child.stat().st_size
    return total


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        raise TypeError("parent-response records do not embed raw arrays")
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if value is None or isinstance(value, str):
        return value
    raise TypeError("unsupported record value " + type(value).__name__)


def _write_exclusive(directory, stem, record, arrays, *, input_paths=()):
    inside = _contains(Path(directory).expanduser().resolve(), REPO)
    destination = authorize_stage_output(
        directory,
        frozen_producer=bool(record.get("frozen_producer")),
        physical_binding=bool(record.get("physical_binding")),
        closure_matches=(record.get("source_closure") == source_closure()) if inside else True,
        input_paths=input_paths,
    )
    destination.mkdir(parents=True, exist_ok=True)
    json_path = destination / f"{stem}.json"
    npz_path = destination / f"{stem}.npz"
    if json_path.exists() or npz_path.exists():
        raise FileExistsError("refusing to overwrite " + stem)
    payload = {key: np.asarray(value) for key, value in arrays.items()}
    array_bytes = sum(value.nbytes for value in payload.values())
    encoded = json.dumps(_jsonable(record), indent=2, sort_keys=True) + "\n"
    if array_bytes + len(encoded.encode()) + _directory_size(destination) > CHUNK_LIMIT_BYTES:
        raise ValueError("exclusive output exceeds 64MiB")
    with npz_path.open("wb") as stream:
        np.savez_compressed(stream, **payload)
    digest = episode.file_sha256(npz_path)
    record = dict(record, payload_sha256=digest, payload_bytes=npz_path.stat().st_size)
    encoded = json.dumps(_jsonable(record), indent=2, sort_keys=True) + "\n"
    if npz_path.stat().st_size + len(encoded.encode()) + _directory_size(destination) - npz_path.stat().st_size > CHUNK_LIMIT_BYTES:
        npz_path.unlink()
        raise ValueError("exclusive output exceeds 64MiB")
    json_path.write_text(encoded)
    return record


def source_closure():
    """Hashes of this adapter and the owners it calls. Not a frozen production commit."""
    import importlib
    owner1 = importlib.import_module("recursive_horizons.nsc_discovery_parent")
    owner2 = importlib.import_module("recursive_horizons.nsc_discovery_parent_episode")
    modules = (leading, response, grandchild, memory, coupling, galerkin, backend, extent, episode, step_control)
    paths = [Path(module.__file__) for module in modules+(owner1, owner2)]
    paths.extend([
        Path(__file__),
        LAB / "scripts" / "derive_nsc_discovery_parent_response.py",
        LAB / "tests" / "test_nsc_discovery_parent_response.py",
        LAB / "docs" / "nsc-discovery-parent-response.md",
    ])
    return {str(path.resolve().relative_to(REPO)): episode.file_sha256(path) for path in paths}


def _binding_arrays(binding):
    holder = binding["holder"]
    state = binding["state"]
    tangent = validate_tangent(holder, binding["tangent"])
    arrays = {name: np.array(getattr(state, name), copy=True) for name in STATE_FIELDS}
    arrays.update({
        "weights": np.array(holder.weights, copy=True),
        "geometry_map": np.array(holder.geometry_map, copy=True),
        "reference_columns": np.array(holder.reference_columns, copy=True),
        "source_columns": np.array(holder.source_columns, copy=True),
        "basis": np.array(binding["basis"], copy=True),
        "occupations_tangent": np.array(tangent.occupations, copy=True),
    })
    for name in STATE_FIELDS:
        arrays["tangent_" + name] = np.array(getattr(tangent, name), copy=True)
    return arrays


def binding_is_physical(binding):
    """True only for an OWNER1 or OWNER2 binding. A manufactured fixture is not physical."""
    return (
        isinstance(binding, dict)
        and binding.get("available") is True
        and binding.get("physical_episode") is True
        and binding.get("kind") in PHYSICAL_BINDINGS
    )


def _validate_held_direction(value):
    """The seeded amplitude 0.0013 is not a veto. The record chooses the new direction."""
    if value is None:
        return None
    if not isinstance(value, dict):
        raise ValueError("a sealed held direction must be a record")
    kind = str(value.get("kind"))
    if kind not in ("population", "gradient"):
        raise ValueError("held direction kind must be population or gradient")
    change = float(value.get("relative_change"))
    if not np.isfinite(change) or change == 0.0:
        raise ValueError("held relative change must be a nonzero finite fraction")
    if value.get("sealed_before_measurement") is not True:
        raise ValueError("held direction must be sealed before the nonlinear measurement")
    baseline = value.get("baseline_amplitude")
    return {
        "kind": kind,
        "relative_change": change,
        "baseline_amplitude": None if baseline is None else float(baseline),
        "sealed_before_measurement": True,
        "baseline_is_measurement_veto": False,
    }


def _input_file_hashes(paths):
    hashes = {}
    for raw in paths:
        path = Path(raw).expanduser().resolve()
        if not path.is_file():
            raise BindingUnavailable("input file is missing; stage open, not a pass: " + str(path))
        hashes[str(path)] = episode.file_sha256(path)
    return hashes


def _require_binding(binding):
    if not isinstance(binding, dict) or binding.get("available") is not True:
        raise BindingUnavailable(
            "initial or input episode binding is unavailable; stage open, not a pass"
        )
    episode_path = binding.get("episode_path")
    if episode_path is not None and not Path(episode_path).is_file():
        raise BindingUnavailable("input episode is not a file; stage open, not a pass")
    for key in ("holder", "state", "tangent", "basis", "widths"):
        if key not in binding:
            raise BindingUnavailable("binding is missing " + key + "; stage open, not a pass")
    return binding


def owner_callbacks():
    """Resolve OWNER1 and OWNER2 callables. This does not call them."""
    import importlib
    owner1 = importlib.import_module("recursive_horizons.nsc_discovery_parent")
    owner2 = importlib.import_module("recursive_horizons.nsc_discovery_parent_episode")
    prepare_parent = getattr(owner1, "prepare_parent", None)
    reconstruct = getattr(owner2, "reconstruct_parent_pair", None)
    loader = getattr(owner2, "load_parent_record", None)
    state_from = getattr(owner2, "state_from_arrays", None)
    if not all(callable(function) for function in (prepare_parent, reconstruct, loader, state_from)):
        raise BindingUnavailable("OWNER1 prepare_parent or OWNER2 cached pair/state is unavailable; stage open")
    return {
        "prepare_parent": prepare_parent,
        "load_parent_record": loader,
        # Compatibility key names the same authenticated public loader; it
        # does not invent a second preparation or bypass record validation.
        "load_parent_preparation": loader,
        "reconstruct_parent_pair": reconstruct,
        "state_from_arrays": state_from,
        "prepare_parent_owner": prepare_parent.__module__,
        "episode_owner": reconstruct.__module__,
        "called": False,
        "loader_return_shape": "(record, arrays)",
        "reconstruction_arguments": "(arrays, record)",
        "state_array_momenta": "pi_Q, pi_r",
        "baseline_seed_amplitude": float(getattr(owner1, "SEED_WEIGHT")),
        "baseline_seed_is_a_measurement_veto": False,
    }


def prepare(directory, binding, *, producer_commit=None):
    """Write the opening record. Defaults do not reach this function.

    A manufactured fixture may be written outside the repository. A new
    ``lab/results/development`` directory also requires a frozen producer and
    a physical OWNER1 or OWNER2 binding. This does not evolve a held arm and
    does not mint a scientific pass.
    """
    binding = _require_binding(binding)
    holder = binding["holder"]
    state = _check_state(holder, binding["state"])
    tangent = validate_tangent(holder, binding["tangent"])
    arrays = _binding_arrays(binding)
    input_paths = tuple(binding.get("input_paths") or ())
    input_hashes = {key: episode.array_sha256(value) for key, value in arrays.items()}
    file_hashes = _input_file_hashes(input_paths)
    regions = regional_map(holder.grid.length)
    physical = binding_is_physical(binding)
    direction = _validate_held_direction(binding.get("sealed_held_direction"))
    record = {
        "schema": SCHEMA,
        "stage": "prepare",
        "kind": str(binding.get("kind", "unspecified")),
        "physical_episode": bool(binding.get("physical_episode", False)),
        "physical_binding": physical,
        "coordinate_time": 0.0,
        "laboratory_dt": float(binding.get("laboratory_dt", 1.0e-4)),
        "step_cap": float(binding.get("step_cap", binding.get("laboratory_dt", 1.0e-4))),
        "nf": int(holder.grid.nf),
        "nq": int(holder.grid.nq),
        "ng": int(holder.grid.ng),
        "length": float(holder.grid.length),
        "widths": [int(value) for value in binding["widths"]],
        "windows": {
            "child_interval": list(holder.child_interval),
            "parent_interval": list(holder.parent_interval),
            "protected_collar_interval": list(regions["protected_collar_interval"]),
            "parent_annulus_intervals": [list(interval) for interval in regions["parent_annulus_intervals"]],
            "centre": regions["centre"],
        },
        "input_hashes": input_hashes,
        "input_file_hashes": file_hashes,
        "input_paths": [str(Path(path).expanduser().resolve()) for path in input_paths],
        "source_closure": source_closure(),
        "producer_commit": producer_commit,
        "frozen_producer": bool(producer_commit),
        "production_authorized": bool(physical and producer_commit),
        "sealed_held_direction": direction,
        "held_direction_sealed": direction is not None,
        "baseline_seed_is_a_measurement_veto": False,
        "scientific_pass": False,
        "analytic_preparation_jacobian": False,
        "child_regional_content": child_regional_content(holder, state),
        "centre_tau_dot": centre_clock(holder, state),
        "correlation_required": True,
        "step_owner": "nsc_discovery_parent_step_control.step_restriction",
    }
    return _write_exclusive(directory, "parent-response-prepare", record, arrays, input_paths=input_paths)


def _load_stage(directory, stem):
    destination = Path(directory).expanduser().resolve()
    if not destination.is_dir():
        raise BindingUnavailable(stem + " record is unavailable; stage open, not a pass")
    json_path = destination / f"{stem}.json"
    npz_path = destination / f"{stem}.npz"
    if not json_path.is_file() or not npz_path.is_file():
        raise BindingUnavailable(stem + " record is unavailable; refusing a placeholder pass")
    record = json.loads(json_path.read_text())
    if episode.file_sha256(npz_path) != record.get("payload_sha256"):
        raise ValueError("payload hash does not match the stage record")
    with np.load(npz_path, allow_pickle=False) as payload:
        arrays = {key: np.array(payload[key], copy=True) for key in payload.files}
    return record, arrays


def _holder_from_prepare(record, arrays):
    grid = galerkin.build_grid(int(record["nf"]), quadrature=int(record["nq"]), length=float(record["length"]), gauge="conformal")
    return make_holder(
        grid,
        weights=arrays["weights"],
        reference_columns=arrays["reference_columns"],
        source_columns=arrays["source_columns"],
        geometry_map=arrays["geometry_map"],
        source_metadata={"kind": record["kind"], "replayed_from_prepare": True},
    )


def _state_tangent_from_arrays(arrays):
    state = leading.State(*(arrays[name] for name in STATE_FIELDS))
    tangent = ParentTangent(*(arrays["tangent_" + name] for name in STATE_FIELDS), arrays["occupations_tangent"])
    return state, tangent


def _assert_inputs_unchanged(record, arrays):
    for key, digest in record["input_hashes"].items():
        if episode.array_sha256(arrays[key]) != digest:
            raise ValueError("prepared input hash changed: " + key)
    for path, digest in record.get("input_file_hashes", {}).items():
        file_path = Path(path)
        if not file_path.is_file() or episode.file_sha256(file_path) != digest:
            raise ValueError("input file changed after prepare: " + path)


def _admitted_step(holder, state, cap):
    dt, admission = step_control.step_restriction(holder, state, cap)
    return float(dt), {
        "dt": float(dt),
        "step_cap": float(admission["step_cap"]),
        "source_rank": int(admission["source_rank"]),
        "restriction_owner": admission["restriction_owner"],
        "stability_certificate": False,
    }


def predict(directory):
    """Lock one tangent forecast from the prepared source. Does not evolve the held arm."""
    record, arrays = _load_stage(directory, "parent-response-prepare")
    if record.get("scientific_pass") is True:
        raise ValueError("prepare record claims a scientific pass")
    holder = _holder_from_prepare(record, arrays)
    state, tangent = _state_tangent_from_arrays(arrays)
    _assert_inputs_unchanged(record, arrays)
    dt, admission = _admitted_step(holder, state, float(record.get("step_cap", record["laboratory_dt"])))
    new_state, new_tangent, clock = advance_tangent(holder, state, tangent, dt)
    matched = matched_centre_readout(holder, new_state, new_tangent, clock["delta_tau"])
    direction = record.get("sealed_held_direction")
    prediction = {
        "schema": SCHEMA,
        "stage": "predict",
        "prepare_payload_sha256": record["payload_sha256"],
        "source_closure": source_closure(),
        "producer_commit": record.get("producer_commit"),
        "frozen_producer": bool(record.get("frozen_producer")),
        "physical_binding": bool(record.get("physical_binding")),
        "production_authorized": bool(record.get("production_authorized")),
        "coordinate_time": dt,
        "admitted_dt": dt,
        "admission": admission,
        "steps": 1,
        "forecast_locked": True,
        "held_arm": False,
        "held_direction_sealed": direction is not None,
        "sealed_held_direction": direction,
        "prediction_before_held_measurement": True,
        "scientific_pass": False,
        "clock": {key: clock[key] for key in ("dt", "tau", "delta_tau", "clock_quadrature")},
        "matched_centre": {key: value for key, value in matched.items() if not isinstance(value, np.ndarray)},
        "baseline_seed_is_a_measurement_veto": False,
        "input_paths": list(record.get("input_paths") or ()),
    }
    payload = {name: np.array(getattr(new_state, name), copy=True) for name in STATE_FIELDS}
    payload.update({"tangent_" + name: np.array(getattr(new_tangent, name), copy=True) for name in STATE_FIELDS})
    payload["occupations_tangent"] = np.array(new_tangent.occupations, copy=True)
    return _write_exclusive(
        directory, "parent-response-predict", prediction, payload,
        input_paths=record.get("input_paths") or (),
    )


def _apply_held_direction(holder, state, tangent, direction):
    change = float(direction["relative_change"])
    if direction["kind"] == "population":
        updated = _live_weights(holder) * (1.0 + change)
        return with_weights(holder, updated), state.copy()
    if direction["kind"] == "gradient":
        return holder, leading.combine(state, _as_state(tangent), change)
    raise ValueError("held direction kind must be population or gradient")


def evolve_prepared_source(holder, state, step_cap):
    """One admitted RK4 step of a fresh prepared source, with its centre clock.

    The timestep comes from ``nsc_discovery_parent_step_control``. This is the
    measurement integrator, not a second stability calculation.
    """
    dt, admission = _admitted_step(holder, state, step_cap)
    advanced = leading.rk4_step(holder, state, dt)
    tau = _nonlinear_centre_increment(holder, state, dt)
    content = child_regional_content(holder, advanced)
    return advanced, {
        "dt": dt,
        "tau": float(tau),
        "child_regional_content": float(content),
        "steps": 1,
        "admission": admission,
        "fresh_prepared_source": True,
        "matched_proper_clock": True,
    }


def _write_measure_record(directory, record, arrays, input_paths):
    if arrays:
        return _write_exclusive(directory, "parent-response-measure", record, arrays, input_paths=input_paths)
    destination = authorize_stage_output(
        directory,
        frozen_producer=bool(record.get("frozen_producer")),
        physical_binding=bool(record.get("physical_binding")),
        closure_matches=True if not _contains(Path(directory).expanduser().resolve(), REPO)
        else record.get("source_closure") == source_closure(),
        input_paths=input_paths,
    )
    path = destination / "parent-response-measure.json"
    if path.exists():
        raise FileExistsError("refusing to overwrite parent-response-measure.json")
    encoded = json.dumps(_jsonable(record), indent=2, sort_keys=True) + "\n"
    if len(encoded.encode()) + _directory_size(destination) > CHUNK_LIMIT_BYTES:
        raise ValueError("exclusive output exceeds 64MiB")
    path.write_text(encoded)
    return record


def measure(directory):
    """Evolve the fresh prepared source when a held direction was sealed first.

    No direction is hardcoded. The baseline amplitude 0.0013 is not a veto.
    A missing seal returns status ``open`` and does not mint a pass. A sealed
    population or gradient direction evolves the prepared source at the matched
    centre clock and keeps the streamed complement forces.
    """
    try:
        forecast, _forecast_arrays = _load_stage(directory, "parent-response-predict")
    except BindingUnavailable as error:
        raise BindingUnavailable("forecast is not locked; stage open, not a pass") from error
    prepared, arrays = _load_stage(directory, "parent-response-prepare")
    _assert_inputs_unchanged(prepared, arrays)
    if forecast.get("prepare_payload_sha256") != prepared.get("payload_sha256"):
        raise ValueError("forecast is not bound to the prepare payload")
    if forecast.get("forecast_locked") is not True:
        raise BindingUnavailable("forecast is not locked; stage open, not a pass")
    input_paths = prepared.get("input_paths") or ()
    common = {
        "schema": SCHEMA,
        "stage": "measure",
        "prediction_payload_sha256": forecast["payload_sha256"],
        "source_closure": source_closure(),
        "producer_commit": prepared.get("producer_commit"),
        "frozen_producer": bool(prepared.get("frozen_producer")),
        "physical_binding": bool(prepared.get("physical_binding")),
        "production_authorized": bool(prepared.get("production_authorized")),
        "scientific_pass": False,
        "baseline_seed_is_a_measurement_veto": False,
    }
    direction = forecast.get("sealed_held_direction")
    if not direction or forecast.get("held_direction_sealed") is not True:
        return _write_measure_record(directory, dict(common, status="open", executed=False,
            reason="held direction is not sealed; measurement stays open"), {}, input_paths)
    holder = _holder_from_prepare(prepared, arrays)
    state, tangent = _state_tangent_from_arrays(arrays)
    cap = float(forecast.get("admitted_dt", prepared.get("step_cap", prepared["laboratory_dt"])))
    baseline, baseline_clock = evolve_prepared_source(holder, state, cap)
    held_holder, held_state = _apply_held_direction(holder, state, tangent, direction)
    held, held_clock = evolve_prepared_source(held_holder, held_state, cap)
    reduced = split_state(holder, state, arrays["basis"], prepared["widths"])
    comparison = full_vs_reduced(holder, state, reduced, baseline_clock["dt"])
    measured = dict(
        common,
        status="measured",
        executed=True,
        held_direction=direction,
        baseline_clock=baseline_clock,
        held_clock=held_clock,
        matched_content_difference=float(held_clock["child_regional_content"] - baseline_clock["child_regional_content"]),
        matched_tau_difference=float(held_clock["tau"] - baseline_clock["tau"]),
        locked_forecast_delta_tau=forecast.get("clock", {}).get("delta_tau"),
        complement_force_max_abs=comparison["forces"]["force_L"]["complement_max_abs"],
        cross_force_max_abs=comparison["forces"]["force_L"]["cross_max_abs"],
        streaming_calls=comparison["rk4_streaming_calls"],
        dense_propagator_stored=False,
        one_physical_metric=True,
        fresh_prepared_source=True,
        physical_held_campaign=bool(prepared.get("physical_binding")),
    )
    payload = {name: np.array(getattr(baseline, name), copy=True) for name in STATE_FIELDS}
    payload.update({"held_" + name: np.array(getattr(held, name), copy=True) for name in STATE_FIELDS})
    return _write_measure_record(directory, measured, payload, input_paths)


def check(directory):
    """Replay hashes and the stage order. An open measurement is not a pass."""
    prepare_record, prepare_arrays = _load_stage(directory, "parent-response-prepare")
    predict_record, _predict_arrays = _load_stage(directory, "parent-response-predict")
    measure_path = Path(directory).expanduser().resolve() / "parent-response-measure.json"
    if not measure_path.is_file():
        raise BindingUnavailable("measure record is unavailable; stage open, not a pass")
    measure_record = json.loads(measure_path.read_text())
    _assert_inputs_unchanged(prepare_record, prepare_arrays)
    if predict_record.get("prepare_payload_sha256") != prepare_record.get("payload_sha256"):
        raise ValueError("prediction is not bound to the prepare payload")
    if predict_record.get("forecast_locked") is not True:
        raise ValueError("prediction was not locked before measurement")
    if measure_record.get("prediction_payload_sha256") != predict_record.get("payload_sha256"):
        raise ValueError("measure record is not bound to the prediction payload")
    status = measure_record.get("status")
    if status not in ("open", "measured"):
        raise ValueError("measure status must be open or measured")
    if status == "open" and measure_record.get("executed") is not False:
        raise ValueError("an open measurement must not claim execution")
    if status == "measured" and measure_record.get("executed") is not True:
        raise ValueError("a measured stage must record the evolution it ran")
    if any(record.get("scientific_pass") is True for record in (prepare_record, predict_record, measure_record)):
        raise ValueError("a stage mints a scientific pass; this adapter does not")
    if _directory_size(Path(directory).expanduser().resolve()) > CHUNK_LIMIT_BYTES:
        raise ValueError("record directory exceeds 64MiB")
    return {
        "schema": SCHEMA,
        "stage": "check",
        "contract_ok": True,
        "scientific_pass": False,
        "measure_status": status,
        "held_measurement_executed": bool(measure_record.get("executed")),
        "forecast_locked_before_measurement": True,
        "inputs_unchanged": True,
        "production_authorized": bool(prepare_record.get("production_authorized")),
        "analytic_preparation_jacobian": False,
        "baseline_seed_is_a_measurement_veto": False,
    }
