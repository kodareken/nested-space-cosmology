"""Exact causal split of the live source, fed back into one common-action geometry.

The retained observer is the fixed original child pair. The source is the full
fermion band. With that fixed frame,

    Phi = V a + e_drive + e_memory
    i a_dot = A a + B (e_drive + e_memory)
    i e_drive_dot = P_E H P_E e_drive,     e_drive(t0) = P_E Phi(t0)
    i e_memory_dot = P_E H P_E e_memory + B† a,   e_memory(t0) = 0

A = V† H V and B e = V† H e, so B† a = P_E H V a. H is the existing Dirac
operator, applied matrix-free. The exterior propagator is not stored.

The reconstructed columns are the source of the existing conformal Galerkin
rates. F_L, F_Q and F_beta are formed before g and the canonical momenta are
advanced, and the multiplicity M = 4 kappa already inside those forces is not
applied again. There is one metric. Child and parent clocks read that metric.

The split is the Volterra decomposition of the same linear law. It is not a
new interaction, not a prescribed-Q schedule, and not a second induced source.
"""
from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import time

import numpy as np

from . import nsc_discovery_backend as backend
from . import nsc_discovery_episode as episode
from . import nsc_nested_parent_child as model
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin


SCHEMA = "NSC-DISCOVERY-COUPLED-MEMORY-CLOSURE-v1"
RUN_SCHEMA = "NSC-DISCOVERY-COUPLED-MEMORY-RUN-v1"
HANDOFF_TIME = 0.3
PRODUCTION_SEGMENT = (0.3, 1.0)
SEGMENT_DURATION = 0.7
SHORT_DURATION = 0.01
MAX_SHORT_DURATION = 0.02
ALLOWED_STEP_CAPS = (0.0005, 0.001, 0.00025)
DEFAULT_STEP_CAP = 0.0005
CHUNK_LIMIT_BYTES = 64 * 1024 * 1024
DEFAULT_CPU_BUDGET_SECONDS = 21600.0
DEFAULT_FORECAST_FACTOR = 1.5
MATCHED_CASES = {
    128: "nf128_baseline_dt0.0005",
    256: "nf256_baseline_dt0.0005",
}
CHILD_COLUMN_INDICES = (2, 3)
SOURCE_RANK = 6
GENERATOR_ATOL = 1e-8
_MODULE = Path(__file__).resolve()
_LAB = _MODULE.parents[2]
_REPO = _MODULE.parents[3]

EXACT_EQUATIONS = {
    "split": "Phi = V a + e_drive + e_memory",
    "retained": "i a_dot = A a + B (e_drive + e_memory)",
    "drive": "i e_drive_dot = P_E H P_E e_drive",
    "drive_initial": "e_drive(t0) = P_E Phi(t0)",
    "memory": "i e_memory_dot = P_E H P_E e_memory + B_dagger a",
    "memory_initial": "e_memory(t0) = 0",
    "blocks": "A = V† H V, B e = V† H e, B† a = P_E H V a",
    "volterra_drive": "e_drive(t) = U_E(t, t0) P_E Phi(t0)",
    "volterra_memory": "e_memory(t) = -i ∫ U_E(t, s) P_E H(s) V a(s) ds",
    "exterior_generator": "i partial_t U_E = P_E H(t) P_E U_E, U_E(t0) = I on the exterior",
    "retained_kernel": "V† H(t) U_E(t, s) P_E H(s) V acts by streaming columns",
    "propagator_stored": False,
    "dense_covariance_stored": False,
    "column_law": "-i H Phi",
    "multiplicity_on_column_law": False,
}

DEPHASING_GAP = {
    "name": "physical_cross_dephasing_preparation",
    "status": "implementation_gap",
    "reason": (
        "A physical cross-dephasing preparation needs an admissible covariance "
        "and a source that remains consistent with the constraints. The live "
        "common-action force is a rank-6 column block. A rank-12 control is not "
        "that block. It is not implemented by truncation, by a hardcoded rank-6 "
        "substitute, or by zeroing the force."
    ),
    "requires_admissible_covariance": True,
    "requires_source_consistent_constraints": True,
    "rank6_block_is_not_a_rank12_control": True,
    "force_zeroed": False,
    "core_modified": False,
    "implemented": False,
    "prerequisite_for_coupled_response": False,
}


class ImplementationGap(RuntimeError):
    """Named missing control. Not a permission to change the force or the core."""


@dataclass
class SplitState:
    """One geometry, one momentum, and the three source pieces. No metric copies."""

    Q: np.ndarray
    r: np.ndarray
    chi: np.ndarray
    p_Q: np.ndarray
    p_r: np.ndarray
    p_chi: np.ndarray
    a: np.ndarray
    e_drive: np.ndarray
    e_memory: np.ndarray

    def copy(self):
        return SplitState(*(np.array(getattr(self, name), copy=True) for name in (
            "Q", "r", "chi", "p_Q", "p_r", "p_chi", "a", "e_drive", "e_memory")))


@dataclass
class SplitRate:
    Q: np.ndarray
    r: np.ndarray
    chi: np.ndarray
    p_Q: np.ndarray
    p_r: np.ndarray
    p_chi: np.ndarray
    a_dot: np.ndarray
    e_drive_dot: np.ndarray
    e_memory_dot: np.ndarray
    force_L: np.ndarray
    force_Q: np.ndarray
    force_beta: np.ndarray
    fieldwork_power: float
    generator_gap: float


def canonical_step_cap(value):
    """Accept the stored cap and the two optional refinements. No other step."""
    number = float(value)
    for cap in ALLOWED_STEP_CAPS:
        if abs(number - cap) <= 1e-15:
            return cap
    raise ValueError("step cap must be 0.0005, 0.001, or 0.00025")


def segment_extent(duration, step_cap):
    """Steps of one coupled segment that starts at the stored T=0.3 and stops by T=1."""
    cap = canonical_step_cap(step_cap)
    span = float(duration)
    steps = int(round(span / cap))
    if steps < 1 or abs(steps * cap - span) > 1e-9:
        raise ValueError("duration must be a positive integer multiple of the step cap")
    final = HANDOFF_TIME + steps * cap
    if final > PRODUCTION_SEGMENT[1] + 1e-12:
        raise ValueError("the coupled segment stops at T=1")
    return {"steps": steps, "dt": cap, "time0": HANDOFF_TIME, "time1": final,
            "duration": steps * cap, "full_segment": abs(span - SEGMENT_DURATION) <= 1e-12}


def cpu_forecast(steps, probe_cpu, *, forecast_factor=DEFAULT_FORECAST_FACTOR,
                 cpu_budget_seconds=DEFAULT_CPU_BUDGET_SECONDS):
    """Linear forecast from one measured coupled step. The factor is the admission margin."""
    if steps < 1 or not np.isfinite(probe_cpu) or probe_cpu < 0:
        raise ValueError("forecast needs a positive step count and a measured probe")
    if not np.isfinite(forecast_factor) or forecast_factor < 1:
        raise ValueError("forecast factor must be at least 1")
    if not np.isfinite(cpu_budget_seconds) or cpu_budget_seconds < 0:
        raise ValueError("CPU budget must be finite and non-negative")
    raw = int(steps) * float(probe_cpu)
    demand = float(forecast_factor) * raw
    return {
        "steps": int(steps),
        "probe_cpu_seconds": float(probe_cpu),
        "forecast_factor": float(forecast_factor),
        "unscaled_cpu_seconds": raw,
        "forecast_cpu_seconds": demand,
        "cpu_budget_seconds": float(cpu_budget_seconds),
        "admitted": bool(demand <= float(cpu_budget_seconds) + 1e-9),
        "checkpoint_limit_bytes": CHUNK_LIMIT_BYTES,
        "executor": "root",
    }


def _require_rank6(columns, weights):
    weight_shape = np.shape(weights)
    column_shape = np.shape(columns)
    if weight_shape != (SOURCE_RANK,) or len(column_shape) != 2 or column_shape[1] != SOURCE_RANK:
        raise ImplementationGap(
            "rank-12 control is an implementation gap named "
            + DEPHASING_GAP["name"]
            + "; refusing to truncate a rank-6 source, to hardcode a substitute block, "
            "or to zero the force"
        )


def fixed_child_observer(pair):
    """Original two child columns, unphased and not replaced by the evolved source."""
    indices = tuple(int(i) for i in pair.child_indices)
    if indices != CHILD_COLUMN_INDICES:
        if any(index < 0 or index >= SOURCE_RANK for index in indices) or len(indices) != 2:
            raise ImplementationGap(
                "observer indices are not the original rank-6 child pair; "
                "geometry-mode indices are not a substitute and the force is left unchanged"
            )
        raise ValueError("this coupling keeps the original child columns (2, 3)")
    observer = np.asarray(pair.original_columns)
    if observer.ndim != 2 or observer.shape[1] != SOURCE_RANK:
        raise ImplementationGap(
            "the fixed observer is not a rank-6 frame; refusing a rank-12 reinterpretation"
        )
    basis = np.array(observer[:, list(indices)], dtype=np.complex128, copy=True)
    gram = basis.conj().T @ basis
    if np.max(np.abs(gram - np.eye(2))) > 1e-8:
        raise ValueError("the fixed child observer is not orthonormal")
    return basis


def _columns_of(state):
    columns = np.vstack((
        np.array(state.phi0, dtype=np.complex128, copy=True),
        np.array(state.phi1, dtype=np.complex128, copy=True),
    ))
    _require_rank6(columns, np.ones(columns.shape[1]))
    return columns


def _spinor(pair, columns):
    nf = int(pair.grid.nf)
    if columns.shape[0] != 2 * nf:
        raise ValueError("source columns must span the full fermion band")
    return columns[:nf], columns[nf:]


def _project_exterior(basis, columns):
    return columns - basis @ (basis.conj().T @ columns)


def reconstruct_columns(pair, split):
    basis = fixed_child_observer(pair)
    return basis @ split.a + split.e_drive + split.e_memory


def nested_from_split(pair, split):
    phi0, phi1 = _spinor(pair, reconstruct_columns(pair, split))
    return model.NestedState(
        np.array(split.Q, dtype=float, copy=True),
        np.array(split.r, dtype=float, copy=True),
        np.array(split.chi, dtype=float, copy=True),
        np.array(split.p_Q, dtype=float, copy=True),
        np.array(split.p_r, dtype=float, copy=True),
        np.array(split.p_chi, dtype=float, copy=True),
        phi0,
        phi1,
    )


def _live_weights(pair):
    weights = np.asarray(pair.weights, dtype=float)
    occupations = np.asarray(pair.grid.fine.occupations, dtype=float)
    _require_rank6(np.zeros((2, SOURCE_RANK), dtype=np.complex128), weights)
    if weights.shape != occupations.shape or not np.array_equal(weights, occupations):
        raise ValueError("source weights and the live occupations must be the same vector")
    if not np.isfinite(weights).all() or np.min(weights) <= 0.0 or np.max(weights) > 1.0:
        raise ValueError("shared source weights must lie in (0, 1]")
    return weights


def initial_split(pair, state):
    """Shared weights. Initial a and e_drive carry the actual cross. Memory starts at 0."""
    weights = _live_weights(pair)
    basis = fixed_child_observer(pair)
    columns = _columns_of(state)
    _require_rank6(columns, weights)
    amplitudes = basis.conj().T @ columns
    drive = columns - basis @ amplitudes
    leakage = float(np.max(np.abs(basis.conj().T @ drive)))
    if leakage > 1e-8:
        raise RuntimeError("initial exterior drive leaked into the child observer")
    memory = np.zeros_like(columns)
    split = SplitState(
        np.array(state.Q, dtype=float, copy=True),
        np.array(state.r, dtype=float, copy=True),
        np.array(state.chi, dtype=float, copy=True),
        np.array(state.p_Q, dtype=float, copy=True),
        np.array(state.p_r, dtype=float, copy=True),
        np.array(state.p_chi, dtype=float, copy=True),
        amplitudes,
        drive,
        memory,
    )
    rebuilt = reconstruct_columns(pair, split)
    gap = float(np.max(np.abs(rebuilt - columns)))
    if gap > 1e-10:
        raise RuntimeError("initial split did not rebuild the full source")
    return split


def _packed_images(pair, nested, basis, amplitudes, drive, memory):
    packed = np.concatenate((basis @ amplitudes, drive, memory), axis=1)
    image = model.apply_hamiltonian(pair, nested, packed)
    width = amplitudes.shape[1]
    return image[:, :width], image[:, width:2 * width], image[:, 2 * width:]


def split_rates(pair, split):
    """Forces from the reconstructed source, then the same g and momenta.

    The geometry derivatives are the existing common-action rates on that
    source. The field derivatives are the exact split of the same H.
    """
    weights = _live_weights(pair)
    del weights
    nested = nested_from_split(pair, split)
    full = model.rates(pair, nested)
    basis = fixed_child_observer(pair)
    drive = _project_exterior(basis, split.e_drive)
    memory = _project_exterior(basis, split.e_memory)
    retained_image, drive_image, memory_image = _packed_images(
        pair, nested, basis, split.a, drive, memory,
    )
    retained_force = _project_exterior(basis, retained_image)
    drive_force = _project_exterior(basis, drive_image)
    memory_force = _project_exterior(basis, memory_image)
    amplitudes_dot = -1j * (basis.conj().T @ (retained_image + drive_image + memory_image))
    drive_dot = -1j * drive_force
    memory_dot = -1j * (memory_force + retained_force)
    split_field = basis @ amplitudes_dot + drive_dot + memory_dot
    full_field = np.vstack((full.phi0, full.phi1))
    gap = float(np.max(np.abs(split_field - full_field)))
    scale = max(1.0, float(np.max(np.abs(full_field))))
    if gap > GENERATOR_ATOL * scale:
        raise RuntimeError("split generator left the common Hamiltonian image")
    multiplicity = int(pair.grid.fine.multiplicity)
    kappa = int(pair.grid.fine.kappa)
    if multiplicity != 4 * kappa:
        raise RuntimeError("multiplicity left M = 4 kappa")
    return SplitRate(
        np.array(full.Q, copy=True),
        np.array(full.r, copy=True),
        np.array(full.chi, copy=True),
        np.array(full.p_Q, copy=True),
        np.array(full.p_r, copy=True),
        np.array(full.p_chi, copy=True),
        amplitudes_dot,
        drive_dot,
        memory_dot,
        np.array(full.force_L, copy=True),
        np.array(full.force_Q, copy=True),
        np.array(full.force_beta, copy=True),
        float(full.fieldwork_power),
        gap,
    )


def _advance(state, rate, factor, *, freeze_geometry=False):
    if freeze_geometry:
        geometry = (state.Q, state.r, state.chi, state.p_Q, state.p_r, state.p_chi)
    else:
        geometry = tuple(
            getattr(state, name) + factor * getattr(rate, name)
            for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi")
        )
    return SplitState(
        *geometry,
        state.a + factor * rate.a_dot,
        state.e_drive + factor * rate.e_drive_dot,
        state.e_memory + factor * rate.e_memory_dot,
    )


def rk4_split(pair, split, dt, *, freeze_geometry=False):
    """Same Butcher weights as the independent full step. Geometry stays one copy."""
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("a positive finite timestep is required")
    k1 = split_rates(pair, split)
    k2 = split_rates(pair, _advance(split, k1, 0.5 * dt, freeze_geometry=freeze_geometry))
    k3 = split_rates(pair, _advance(split, k2, 0.5 * dt, freeze_geometry=freeze_geometry))
    k4 = split_rates(pair, _advance(split, k3, dt, freeze_geometry=freeze_geometry))
    combined = SplitRate(
        (k1.Q + 2 * k2.Q + 2 * k3.Q + k4.Q) / 6,
        (k1.r + 2 * k2.r + 2 * k3.r + k4.r) / 6,
        (k1.chi + 2 * k2.chi + 2 * k3.chi + k4.chi) / 6,
        (k1.p_Q + 2 * k2.p_Q + 2 * k3.p_Q + k4.p_Q) / 6,
        (k1.p_r + 2 * k2.p_r + 2 * k3.p_r + k4.p_r) / 6,
        (k1.p_chi + 2 * k2.p_chi + 2 * k3.p_chi + k4.p_chi) / 6,
        (k1.a_dot + 2 * k2.a_dot + 2 * k3.a_dot + k4.a_dot) / 6,
        (k1.e_drive_dot + 2 * k2.e_drive_dot + 2 * k3.e_drive_dot + k4.e_drive_dot) / 6,
        (k1.e_memory_dot + 2 * k2.e_memory_dot + 2 * k3.e_memory_dot + k4.e_memory_dot) / 6,
        k1.force_L, k1.force_Q, k1.force_beta, k1.fieldwork_power, max(
            k1.generator_gap, k2.generator_gap, k3.generator_gap, k4.generator_gap,
        ),
    )
    out = _advance(split, combined, dt, freeze_geometry=freeze_geometry)
    model._active(pair, nested_from_split(pair, out))
    return out


def independent_equations(pair, split):
    """Recompute each split equation from separate Hamiltonian images."""
    nested = nested_from_split(pair, split)
    basis = fixed_child_observer(pair)
    columns = reconstruct_columns(pair, split)
    drive = _project_exterior(basis, split.e_drive)
    memory = _project_exterior(basis, split.e_memory)
    retained = basis @ split.a
    image_phi = model.apply_hamiltonian(pair, nested, columns)
    image_drive = model.apply_hamiltonian(pair, nested, drive)
    image_memory = model.apply_hamiltonian(pair, nested, memory)
    image_retained = model.apply_hamiltonian(pair, nested, retained)
    amplitudes_dot = -1j * (basis.conj().T @ image_phi)
    drive_dot = -1j * _project_exterior(basis, image_drive)
    memory_dot = -1j * (
        _project_exterior(basis, image_memory) + _project_exterior(basis, image_retained)
    )
    return {
        "a_dot": amplitudes_dot,
        "e_drive_dot": drive_dot,
        "e_memory_dot": memory_dot,
        "phi_dot": -1j * image_phi,
        "dense_exterior_propagator": None,
        "equations": dict(EXACT_EQUATIONS),
    }


def _source_forces(pair, geometry, columns):
    phi0, phi1 = _spinor(pair, columns)
    trial = model.NestedState(
        np.array(geometry.Q, copy=True), np.array(geometry.r, copy=True),
        np.array(geometry.chi, copy=True), np.array(geometry.p_Q, copy=True),
        np.array(geometry.p_r, copy=True), np.array(geometry.p_chi, copy=True),
        phi0, phi1,
    )
    _fine, system = model._active(pair, trial)
    source = coupling.source_from_columns(system, _fine)
    if source["multiplicity_applied_once"] is not True:
        raise RuntimeError("source multiplicity was not applied once")
    if int(system.multiplicity) != 4 * int(system.kappa):
        raise RuntimeError("source multiplicity drifted")
    return source


def cross_force_report(pair, split):
    """Mean force of the full source, with the cross piece kept separate from the modal kernel."""
    basis = fixed_child_observer(pair)
    weights = _live_weights(pair)
    full_columns = reconstruct_columns(pair, split)
    pieces = {
        "retained": basis @ split.a,
        "drive": split.e_drive,
        "memory": split.e_memory,
    }
    full = _source_forces(pair, split, full_columns)
    partial = {name: _source_forces(pair, split, columns) for name, columns in pieces.items()}
    forces = {}
    for name in ("force_L", "force_Q", "force_beta"):
        dropped = partial["retained"][name] + partial["drive"][name] + partial["memory"][name]
        cross = full[name] - dropped
        forces[name] = {
            "cross_max_abs": float(np.max(np.abs(cross))),
            "cross_l2": float(np.linalg.norm(cross)),
            "full_max_abs": float(np.max(np.abs(full[name]))),
            "geometry_uses_full_source": True,
            "cross_omitted_from_geometry": False,
        }
    exterior = split.e_drive + split.e_memory
    modal = (split.a * weights[None, :]) @ split.a.conj().T
    cross_block = (split.a * weights[None, :]) @ exterior.conj().T
    thin = (split.a * weights[None, :]) @ full_columns.conj().T
    thin_gap = float(np.max(np.abs(thin - (modal @ basis.conj().T + cross_block))))
    return {
        "forces": forces,
        "modal_kernel_shape": list(modal.shape),
        "modal_kernel_frobenius": float(np.linalg.norm(modal)),
        "retained_exterior_cross_frobenius": float(np.linalg.norm(cross_block)),
        "thin_row_identity_gap": thin_gap,
        "modal_kernel_is_spatial_stress": False,
        "modal_kernel_used_as_force": False,
        "memory_self_energy_used_as_force": False,
        "shared_weights": True,
        "dense_covariance_formed": False,
        "multiplicity_applied_once": True,
        "full_force_L": full["force_L"],
        "full_force_Q": full["force_Q"],
        "full_force_beta": full["force_beta"],
    }


def car_report(pair, phi0, phi1):
    """Nonzero spectrum of the column Gaussian, without assembling the dense covariance."""
    weights = _live_weights(pair)
    columns = np.vstack((phi0, phi1))
    gram = columns.conj().T @ columns
    scaled = np.sqrt(weights)
    weighted = scaled[:, None] * gram * scaled[None, :]
    weighted = 0.5 * (weighted + weighted.conj().T)
    spectrum = np.linalg.eigvalsh(weighted).real
    return {
        "spectrum_min": float(spectrum[0]),
        "spectrum_max": float(spectrum[-1]),
        "admissible": bool(spectrum[0] >= -1e-8 and spectrum[-1] <= 1.0 + 1e-8),
        "gram_hermitian_gap": float(np.max(np.abs(gram - gram.conj().T))),
        "dense_covariance_formed": False,
    }


def clock_rates(pair, state):
    """Proper-time rates of the one reconstructed metric at x = 1, 2, 3."""
    reading = model.metrics(pair, state)
    if reading["one_physical_metric"] is not True:
        raise RuntimeError("clocks were not read from one metric")
    rates = np.array(reading["clock_rates"], dtype=np.float64)
    if rates.shape != (3,):
        raise RuntimeError("the three fixed observers did not share one clock sample")
    return rates


def _trap(samples, dt):
    """Simpson on an even interval count; trapezoid otherwise. Not a propagator."""
    intervals = len(samples) - 1
    if intervals >= 2 and intervals % 2 == 0:
        total = samples[0] + samples[-1]
        for index, value in enumerate(samples[1:-1], start=1):
            total += (4.0 if index % 2 else 2.0) * value
        return total * dt / 3.0
    total = 0.0j
    for left, right in zip(samples, samples[1:]):
        total += 0.5 * dt * (left + right)
    return total


def _weighted_overlap(left, right, weights):
    return complex(np.sum(weights * np.sum(np.conjugate(left) * right, axis=0)))


def _lifted_direction(pair, direction):
    return pair.grid.A_g @ (pair.geometry_map @ np.asarray(direction, dtype=float))


def same_action_partials(pair, state, *, eps=1e-6):
    """Centred derivatives of the field energy against the existing nodal forces.

    The state is whatever coupled evolution produced. Rank-12 dephasing is not
    consulted. The chart contraction uses F_L + F_Q. The shift contraction is
    the off-chart partial F_beta.
    """
    fine, system = model._active(pair, state)
    source = coupling.source_from_columns(system, fine)
    direction = np.zeros(pair.grid.ng)
    direction[0] = 1.0
    lifted = _lifted_direction(pair, direction)

    def shifted_energy(sign):
        trial = state.copy()
        trial.Q = np.array(state.Q, copy=True) + sign * eps * direction
        trial_fine, trial_system = model._active(pair, trial)
        return coupling.field_energy(trial_system, trial_fine)

    def partial_energy(name, sign):
        length = system.length_density
        shift = system.shift
        if name == "L":
            length = length + sign * eps * lifted
        else:
            shift = shift + sign * eps * lifted
        cloned = replace(system, length_density=np.array(length, copy=True), shift=np.array(shift, copy=True))
        return coupling.field_energy(cloned, fine)

    chart_numeric = (shifted_energy(1.0) - shifted_energy(-1.0)) / (2.0 * eps)
    chart_analytic = float(np.dot(source["force_L"] + source["force_Q"], lifted))
    lapse_numeric = (partial_energy("L", 1.0) - partial_energy("L", -1.0)) / (2.0 * eps)
    lapse_analytic = float(np.dot(source["force_L"], lifted))
    shift_numeric = (partial_energy("beta", 1.0) - partial_energy("beta", -1.0)) / (2.0 * eps)
    shift_analytic = float(np.dot(source["force_beta"], lifted))
    return {
        "chart_numeric": float(chart_numeric),
        "chart_analytic": float(chart_analytic),
        "chart_gap": float(abs(chart_numeric - chart_analytic)),
        "lapse_partial_gap": float(abs(lapse_numeric - lapse_analytic)),
        "shift_partial_gap": float(abs(shift_numeric - shift_analytic)),
        "force_L": source["force_L"],
        "force_Q": source["force_Q"],
        "force_beta": source["force_beta"],
        "dephasing_consulted": False,
    }


def _quadratic_partition(pair, geometry, retained_columns, exterior_columns):
    """Full force, exterior complementary force, and their retained cross.

    Complementary force is the common-action force of the exterior vector.
    The cross is the remainder required by the shared weights. Effective
    outside force is complementary plus that cross.
    """
    full_columns = retained_columns + exterior_columns
    full = _source_forces(pair, geometry, full_columns)
    retained = _source_forces(pair, geometry, retained_columns)
    complementary = _source_forces(pair, geometry, exterior_columns)
    pieces = {}
    for name in ("force_L", "force_Q", "force_beta"):
        cross = full[name] - retained[name] - complementary[name]
        outside = complementary[name] + cross
        pieces[name] = {
            "full": full[name],
            "retained": retained[name],
            "complementary": complementary[name],
            "cross": cross,
            "effective_outside": outside,
        }
    return pieces


def endpoint_force_comparison(pair, full_state, split):
    """Compare force partitions on the evolved full state and the evolved split.

    Memory is the memory the coupled step produced. It is not reset to zero.
    """
    basis = fixed_child_observer(pair)
    columns = np.vstack((
        np.array(full_state.phi0, dtype=np.complex128, copy=True),
        np.array(full_state.phi1, dtype=np.complex128, copy=True),
    ))
    amplitudes = basis.conj().T @ columns
    exterior = columns - basis @ amplitudes
    full_part = _quadratic_partition(pair, full_state, basis @ amplitudes, exterior)
    split_state = nested_from_split(pair, split)
    split_part = _quadratic_partition(
        pair, split_state, basis @ split.a, split.e_drive + split.e_memory,
    )
    partials = same_action_partials(pair, full_state)
    names = ("force_L", "force_Q", "force_beta")
    report = {
        "endpoint_memory_norm": float(np.linalg.norm(split.e_memory)),
        "endpoint_drive_norm": float(np.linalg.norm(split.e_drive)),
        "endpoint_response_uses_evolved_split": True,
        "memory_reset_at_endpoint": False,
        "dephasing_prerequisite": False,
        "endpoint_chart_gap": partials["chart_gap"],
        "endpoint_lapse_partial_gap": partials["lapse_partial_gap"],
        "endpoint_shift_partial_gap": partials["shift_partial_gap"],
        "full_partition": full_part,
        "split_partition": split_part,
    }
    for label in ("full", "complementary", "cross", "effective_outside"):
        report[label + "_force_max_abs" if label != "full" else "endpoint_rebuild_force_max_abs"] = max(
            _max_abs(full_part[name][label], split_part[name][label]) for name in names
        )
    report["effective_outside_force_max_abs"] = report["effective_outside_force_max_abs"]
    return report


def elimination_variation(pair, state, *, eps=1e-6, laboratory_step=2e-4, laboratory_steps=4):
    """Common-action partials, plus the exterior boundary term.

    At fixed reconstructed source the chart derivative of the field energy is
    the contraction of F_L + F_Q, and the off-chart shift derivative is F_beta.
    Those are the forces already passed to the rate.

    After the exterior equation is imposed, a variation that vanishes on the
    initial exterior state integrates to an endpoint term. That term is not
    written into F_L, F_Q, or F_beta. The laboratory that exhibits the endpoint
    term freezes geometry only as an identity check; it is not the coupled law
    and it stores no exterior propagator.
    """
    split = initial_split(pair, state)
    nested = nested_from_split(pair, split)
    partials = same_action_partials(pair, nested, eps=eps)
    chart_numeric = partials["chart_numeric"]
    chart_analytic = partials["chart_analytic"]
    source = {"force_L": partials["force_L"], "force_Q": partials["force_Q"],
              "force_beta": partials["force_beta"]}

    basis = fixed_child_observer(pair)
    weights = _live_weights(pair)
    samples_e = []
    samples_edot = []
    ordered_states = []
    cursor = split
    step = float(laboratory_step)
    count = int(laboratory_steps)
    if count < 2 or step <= 0.0:
        raise ValueError("the identity laboratory needs at least two forward samples")
    duration = step * count
    initial_force = None
    for _index in range(count):
        rate = split_rates(pair, cursor)
        if initial_force is None:
            initial_force = np.array(rate.force_L, copy=True)
        ordered_states.append(cursor)
        samples_e.append(cursor.e_drive + cursor.e_memory)
        samples_edot.append(rate.e_drive_dot + rate.e_memory_dot)
        cursor = rk4_split(pair, cursor, step, freeze_geometry=True)
    final_rate = split_rates(pair, cursor)
    ordered_states.append(cursor)
    samples_e.append(cursor.e_drive + cursor.e_memory)
    samples_edot.append(final_rate.e_drive_dot + final_rate.e_memory_dot)
    eta = _project_exterior(basis, samples_e[-1])
    eta_norm = np.sqrt(max(_weighted_overlap(eta, eta, weights).real, 1e-30))
    eta = eta / eta_norm
    nodes = np.linspace(0.0, duration, count + 1)
    kinetic = []
    adjoint = []
    residual = []
    for time_value, exterior, exterior_dot, sample_state in zip(nodes, samples_e, samples_edot, ordered_states):
        probe = (time_value / duration) * eta
        kinetic.append(_weighted_overlap(probe, 1j * exterior_dot, weights))
        adjoint.append(_weighted_overlap(eta / duration, 1j * exterior, weights))
        residual.append(_weighted_overlap(
            probe, 1j * exterior_dot - _exterior_generator(pair, sample_state, exterior), weights,
        ))
    boundary = _weighted_overlap(eta, 1j * samples_e[-1], weights)
    kinetic_integral = _trap(kinetic, step)
    adjoint_integral = _trap(adjoint, step)
    by_parts_gap = kinetic_integral - (boundary - adjoint_integral)
    residual_integral = _trap(residual, step)
    force_gap = float(np.max(np.abs(initial_force - source["force_L"])))
    return {
        "chart_numeric": partials["chart_numeric"],
        "chart_analytic": partials["chart_analytic"],
        "chart_gap": partials["chart_gap"],
        "lapse_partial_gap": partials["lapse_partial_gap"],
        "shift_partial_gap": partials["shift_partial_gap"],
        "boundary_term": [float(boundary.real), float(boundary.imag)],
        "boundary_abs": float(abs(boundary)),
        "integration_by_parts_gap": float(abs(by_parts_gap)),
        "exterior_residual_integral_abs": float(abs(residual_integral)),
        "force_matches_reconstructed_source": bool(force_gap <= 1e-10),
        "boundary_inserted_into_force": False,
        "force_source_gap": force_gap,
        "second_induced_source": False,
        "laboratory_is_coupled_production": False,
        "laboratory_freezes_geometry": True,
        "dense_exterior_propagator_stored": False,
        "geometry_freeze_is_the_physical_law": False,
    }


def _exterior_generator(pair, split, exterior):
    """P_E H P_E e + P_E H V a at the frozen laboratory geometry."""
    basis = fixed_child_observer(pair)
    nested = nested_from_split(pair, split)
    image_exterior = model.apply_hamiltonian(pair, nested, _project_exterior(basis, exterior))
    image_retained = model.apply_hamiltonian(pair, nested, basis @ split.a)
    return _project_exterior(basis, image_exterior) + _project_exterior(basis, image_retained)


def omission_diagnostic(pair, split, which):
    """Counterfactual projection. Not a unitary full-band law and not a gate."""
    if which not in ("memory", "drive"):
        raise ValueError("omission must name memory or drive")
    rate = split_rates(pair, split)
    basis = fixed_child_observer(pair)
    drive_dot = np.array(rate.e_drive_dot, copy=True)
    memory_dot = np.array(rate.e_memory_dot, copy=True)
    if which == "memory":
        memory_dot[:] = 0.0
    else:
        drive_dot[:] = 0.0
    projected = basis @ rate.a_dot + drive_dot + memory_dot
    columns = reconstruct_columns(pair, split)
    full = -1j * model.apply_hamiltonian(pair, nested_from_split(pair, split), columns)
    residual = float(np.max(np.abs(projected - full)))
    return {
        "which": which,
        "residual_max_abs": residual,
        "unitary_global_law": False,
        "autonomous_regeneration": False,
        "defines_coupled_geometry": False,
        "physical_gate": False,
        "scope": (
            "counterfactual memory or drive projection of the retained description; "
            "the projected vector is not the Hamiltonian image of the full source"
        ),
    }


def _max_abs(left, right):
    return float(np.max(np.abs(np.asarray(left) - np.asarray(right))))


def compare_step(pair, state, dt):
    """One independent full RK4 step against one streamed split step."""
    full_before = state.copy()
    split = initial_split(pair, full_before)
    rate = split_rates(pair, split)
    reference = independent_equations(pair, split)
    equation_gap = max(
        _max_abs(rate.a_dot, reference["a_dot"]),
        _max_abs(rate.e_drive_dot, reference["e_drive_dot"]),
        _max_abs(rate.e_memory_dot, reference["e_memory_dot"]),
    )
    full_after = model.rk4_step(pair, full_before, dt)
    split_after = rk4_split(pair, split, dt)
    nested_after = nested_from_split(pair, split_after)
    crosses = cross_force_report(pair, split)
    return _comparison_record(
        pair, full_before, full_after, nested_after, rate, equation_gap, crosses, steps=1, dt=dt,
    )


def _comparison_record(pair, full_before, full_after, nested_after, rate, equation_gap, crosses,
                       *, steps, dt, clocks_full=None, clocks_split=None, generator_gap=None,
                       energy_full=None, energy_split=None):
    field_gap = max(_max_abs(full_after.phi0, nested_after.phi0), _max_abs(full_after.phi1, nested_after.phi1))
    geometry_gap = max(_max_abs(getattr(full_after, name), getattr(nested_after, name))
                       for name in model.GEOMETRY_NAMES + model.MOMENTUM_NAMES)
    full_rate = model.rates(pair, full_after)
    split_rate = model.rates(pair, nested_after)
    force_gap = max(
        _max_abs(full_rate.force_L, split_rate.force_L),
        _max_abs(full_rate.force_Q, split_rate.force_Q),
        _max_abs(full_rate.force_beta, split_rate.force_beta),
    )
    clock_gap = _max_abs(clock_rates(pair, full_after), clock_rates(pair, nested_after))
    car_full = car_report(pair, full_after.phi0, full_after.phi1)
    car_split = car_report(pair, nested_after.phi0, nested_after.phi1)
    account = model.energy_accounting(pair, nested_after)
    if energy_full is None:
        energy_full = float(model.energy(pair, full_after))
        energy_split = float(model.energy(pair, nested_after))
    clock_integral_gap = None
    if clocks_full is not None:
        clock_integral_gap = _max_abs(clocks_full, clocks_split)
    initial_cross = cross_force_report(pair, initial_split(pair, full_before))
    return {
        "steps": int(steps),
        "dt": float(dt),
        "field_max_abs": field_gap,
        "covariance_spectrum_max_abs": max(
            abs(car_full["spectrum_min"] - car_split["spectrum_min"]),
            abs(car_full["spectrum_max"] - car_split["spectrum_max"]),
        ),
        "force_max_abs": force_gap,
        "geometry_max_abs": geometry_gap,
        "clock_rate_max_abs": clock_gap,
        "clock_integral_max_abs": clock_integral_gap,
        "equation_gap": float(equation_gap),
        "generator_gap": float(rate.generator_gap if generator_gap is None else generator_gap),
        "initial_cross_frobenius": initial_cross["retained_exterior_cross_frobenius"],
        "initial_cross_force_max_abs": max(
            initial_cross["forces"][name]["cross_max_abs"] for name in ("force_L", "force_Q", "force_beta")
        ),
        "cross_retained_in_rate": bool(np.allclose(rate.force_L, crosses["full_force_L"])
                                       and np.allclose(rate.force_Q, crosses["full_force_Q"])
                                       and np.allclose(rate.force_beta, crosses["full_force_beta"])),
        "modal_kernel_used_as_force": False,
        "memory_self_energy_used_as_force": False,
        "one_physical_metric": True,
        "car_full_admissible": car_full["admissible"],
        "car_split_admissible": car_split["admissible"],
        "car_spectrum_split": [car_split["spectrum_min"], car_split["spectrum_max"]],
        "energy_max_abs": abs(energy_full - energy_split),
        "energy_field_closure_abs": abs(float(account["field_closure_error"])),
        "multiplicity_applied_once": bool(account["angular_multiplicity_applied_once"]),
        "omission_used_as_gate": False,
        "dense_covariance_stored": False,
        "dense_exterior_propagator_stored": False,
        "prescribed_q": False,
    }


def compare_window(pair, state, dt, duration, clocks, *, process_budget_seconds=None, sink=None):
    """Coupled window. Clocks advance with the evolved metric, from the stored clock.

    A process budget stops before the next step and returns the last admissible
    state. The short-test helper imposes its own duration cap. This function
    accepts every segment that ends by T=1, including the full 0.7 span.
    """
    extent = segment_extent(duration, dt)
    steps = extent["steps"]
    dt = extent["dt"]
    full = state.copy()
    split = initial_split(pair, full)
    clocks_full = np.array(clocks, dtype=np.float64)
    clocks_split = np.array(clocks, dtype=np.float64)
    if clocks_full.shape != (3,):
        raise ValueError("continuation clocks are the three samples of one metric")
    worst_field = 0.0
    worst_geometry = 0.0
    worst_force = 0.0
    worst_outside = 0.0
    worst_clock = 0.0
    worst_generator = 0.0
    equation_gap = 0.0
    first_rate = None
    first_cross = None
    started = time.process_time()
    status = "completed"
    taken = 0
    for _index in range(steps):
        if process_budget_seconds is not None and time.process_time() - started > float(process_budget_seconds):
            status = "budget_stop"
            break
        rate = split_rates(pair, split)
        reference = independent_equations(pair, split)
        equation_gap = max(equation_gap, _max_abs(rate.a_dot, reference["a_dot"]),
                           _max_abs(rate.e_drive_dot, reference["e_drive_dot"]),
                           _max_abs(rate.e_memory_dot, reference["e_memory_dot"]))
        worst_generator = max(worst_generator, rate.generator_gap)
        if first_rate is None:
            first_rate = rate
            first_cross = cross_force_report(pair, split)
        left_full = clock_rates(pair, full)
        left_split = clock_rates(pair, nested_from_split(pair, split))
        full = model.rk4_step(pair, full, dt)
        split = rk4_split(pair, split, dt)
        nested = nested_from_split(pair, split)
        right_full = clock_rates(pair, full)
        right_split = clock_rates(pair, nested)
        clocks_full = clocks_full + 0.5 * dt * (left_full + right_full)
        clocks_split = clocks_split + 0.5 * dt * (left_split + right_split)
        worst_field = max(worst_field, _max_abs(full.phi0, nested.phi0), _max_abs(full.phi1, nested.phi1))
        worst_geometry = max(worst_geometry, max(
            _max_abs(getattr(full, name), getattr(nested, name))
            for name in model.GEOMETRY_NAMES + model.MOMENTUM_NAMES
        ))
        full_rate = model.rates(pair, full)
        split_rate = model.rates(pair, nested)
        worst_force = max(worst_force, _max_abs(full_rate.force_L, split_rate.force_L),
                          _max_abs(full_rate.force_Q, split_rate.force_Q),
                          _max_abs(full_rate.force_beta, split_rate.force_beta))
        worst_clock = max(worst_clock, _max_abs(right_full, right_split), _max_abs(clocks_full, clocks_split))
        taken += 1
    if first_rate is None:
        halted = {
            "steps": 0, "dt": float(dt), "duration": float(extent["duration"]),
            "status": status, "time0": HANDOFF_TIME, "time1": HANDOFF_TIME,
            "production_steps_taken": 0, "dephasing_prerequisite": False,
        }
        if sink is not None:
            sink.append({"full": full, "split": split, "clocks_full": clocks_full,
                         "clocks_split": clocks_split, "endpoint": None})
        return halted
    nested = nested_from_split(pair, split)
    endpoint = endpoint_force_comparison(pair, full, split)
    worst_outside = endpoint["effective_outside_force_max_abs"]
    record = _comparison_record(
        pair, state, full, nested, first_rate, equation_gap, first_cross,
        steps=taken, dt=dt, clocks_full=clocks_full, clocks_split=clocks_split,
        generator_gap=worst_generator,
        energy_full=float(model.energy(pair, full)),
        energy_split=float(model.energy(pair, nested)),
    )
    record["field_max_abs"] = max(record["field_max_abs"], worst_field)
    record["geometry_max_abs"] = max(record["geometry_max_abs"], worst_geometry)
    record["force_max_abs"] = max(record["force_max_abs"], worst_force)
    record["clock_rate_max_abs"] = max(record["clock_rate_max_abs"], worst_clock)
    record["duration"] = float(taken * dt)
    record["requested_steps"] = int(steps)
    record["status"] = status
    record["time0"] = HANDOFF_TIME
    record["time1"] = HANDOFF_TIME + taken * dt
    record["clocks_follow_evolved_metric"] = True
    record["stored_handoff_rates_reused_as_future_rates"] = False
    record["geometry_motion_max_abs"] = max(
        _max_abs(getattr(full, name), getattr(state, name))
        for name in model.GEOMETRY_NAMES + model.MOMENTUM_NAMES
    )
    record["field_motion_max_abs"] = max(_max_abs(full.phi0, state.phi0), _max_abs(full.phi1, state.phi1))
    record["clock_advance_max_abs"] = _max_abs(clocks_full, np.asarray(clocks, dtype=float))
    record["effective_outside_force_max_abs"] = worst_outside
    record["complementary_force_max_abs"] = endpoint["complementary_force_max_abs"]
    record["endpoint_cross_force_max_abs"] = endpoint["cross_force_max_abs"]
    record["endpoint_memory_norm"] = endpoint["endpoint_memory_norm"]
    record["endpoint_response_uses_evolved_split"] = True
    record["memory_reset_at_endpoint"] = False
    record["dephasing_prerequisite"] = False
    record["endpoint_chart_gap"] = endpoint["endpoint_chart_gap"]
    record["endpoint_lapse_partial_gap"] = endpoint["endpoint_lapse_partial_gap"]
    record["endpoint_shift_partial_gap"] = endpoint["endpoint_shift_partial_gap"]
    if sink is not None:
        sink.append({"full": full, "split": split, "clocks_full": clocks_full,
                     "clocks_split": clocks_split, "endpoint": endpoint})
    return record


def _dependency_paths():
    names = (
        "src/recursive_horizons/nsc_discovery_coupled_memory.py",
        "src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
        "src/recursive_horizons/nsc_nested_parent_child.py",
        "src/recursive_horizons/nsc_spherical_coupling.py",
        "src/recursive_horizons/nsc_discovery_backend.py",
        "src/recursive_horizons/nsc_discovery_episode.py",
        "src/recursive_horizons/nsc_spherical_feedback_action.py",
    )
    return {name: episode.file_sha256(_LAB / name) for name in names}


def load_matched_handoff(fermions, *, backend_name="fft"):
    """Exact stored T=0.3 state, frozen W, and the FFT carrier. No new solve."""
    case_name = MATCHED_CASES[int(fermions)]
    handoff = episode.load_saved_handoff(case_name)
    if int(handoff["nf"]) != int(fermions):
        raise RuntimeError("handoff resolution does not match the requested case")
    source = handoff["source_columns"]
    nf = int(handoff["nf"])
    arrays = {
        "Q": handoff["state"].Q,
        "r": handoff["state"].r,
        "chi": handoff["state"].chi,
        "pi_Q": handoff["state"].p_Q,
        "pi_r": handoff["state"].p_r,
        "pi_chi": handoff["state"].p_chi,
        "phi0": handoff["state"].phi0,
        "phi1": handoff["state"].phi1,
        "W": handoff["W"],
        "source_phi0": source[:nf],
        "source_phi1": source[nf:],
        "observer_columns": handoff["observer_columns"],
        "source_weights": handoff["weights"],
    }
    record = {
        "nf": nf,
        "coarse_indices": handoff["coarse_indices"],
        "child_indices": handoff["child_indices"],
        "parent_indices": handoff["parent_indices"],
        "source_metadata": handoff["source_metadata"],
        "geometry_metadata": handoff["geometry_metadata"],
        "clock_locations": list(handoff["clocks"]["locations"]),
    }
    pair = episode.pair_from_arrays(arrays, record)
    pair, state, info = episode.resolve_pair(pair, handoff["state"], backend=backend_name)
    if info.get("W_changed") or info.get("state_changed") or info.get("initial_state_called"):
        raise RuntimeError("carrier resolution changed W, the stored state, or called initial_state")
    if backend_name == "fft" and info.get("backend") != "fft":
        raise RuntimeError("matrix-free FFT carrier was not applied")
    if tuple(int(i) for i in pair.child_indices) != CHILD_COLUMN_INDICES:
        raise RuntimeError("stored pair did not keep the original child column indices")
    if not np.array_equal(np.asarray(pair.geometry_map), np.asarray(handoff["W"])):
        raise RuntimeError("canonical W was rebuilt")
    clocks = np.array(handoff["clocks"]["normal_clocks"], dtype=np.float64)
    return {
        "case": case_name,
        "nf": nf,
        "pair": pair,
        "state": state,
        "dt": float(handoff["step_cap"]),
        "clocks": clocks,
        "clock_locations": list(handoff["clocks"]["locations"]),
        "pins": dict(handoff["pins"]),
        "time": HANDOFF_TIME,
        "initial_state_called": False,
        "W_rebuilt": False,
        "backend": info.get("backend"),
        "fft_applied": bool(info.get("fft_applied")),
    }


def _array_sha(array):
    return episode.array_sha256(array)


def _cap_token(step_cap):
    return {0.0005: "0.0005", 0.001: "0.001", 0.00025: "0.00025"}[canonical_step_cap(step_cap)]


def _campaign_directory(path, *, create):
    resolved = Path(path).expanduser().resolve()
    if resolved == _REPO or _REPO in resolved.parents:
        development = _REPO / "lab/results/development"
        if resolved.parent != development or not resolved.name.startswith("nsc-discovery-coupled-memory-"):
            raise PermissionError("repository campaigns require a new coupled-memory successor under lab/results/development")
    if create:
        resolved.mkdir(parents=True, exist_ok=True)
        if any(resolved.iterdir()):
            raise FileExistsError("refusing to overwrite an existing campaign directory")
    elif not resolved.is_dir():
        raise FileNotFoundError("campaign directory is missing: " + str(resolved))
    return resolved


def _manifest_path(directory):
    return Path(directory) / "manifest.json"


def _write_manifest(directory, manifest):
    path = _manifest_path(directory)
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    temporary.write_text(json.dumps(_jsonable(manifest), indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _read_manifest(directory):
    path = _manifest_path(directory)
    if not path.is_file():
        raise FileNotFoundError("campaign manifest is missing: " + str(path))
    return json.loads(path.read_text())


def _commit_chunk(directory, case_id, ordinal, arrays, record):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    stored = {}
    for name, value in arrays.items():
        array = np.asarray(value)
        if np.iscomplexobj(array):
            stored[name] = np.ascontiguousarray(array, dtype=np.complex128)
        else:
            stored[name] = np.ascontiguousarray(np.real(array), dtype=np.float64)
    temporary = directory / f".{case_id}-{ordinal:06d}.{os.getpid()}.npz"
    with temporary.open("wb") as stream:
        np.savez(stream, **stored)
    payload = temporary.read_bytes()
    temporary.unlink()
    body = dict(record)
    body.update({
        "schema": RUN_SCHEMA,
        "case_id": case_id,
        "ordinal": int(ordinal),
        "arrays_sha256": hashlib.sha256(payload).hexdigest(),
        "array_sha256": {name: _array_sha(stored[name]) for name in stored},
        "immutable": True,
        "checkpoint_limit_bytes": CHUNK_LIMIT_BYTES,
        "dense_exterior_propagator_stored": False,
        "dense_covariance_stored": False,
    })
    text = json.dumps(_jsonable(body), indent=2, sort_keys=True) + "\n"
    if len(payload) + len(text.encode()) > CHUNK_LIMIT_BYTES:
        raise RuntimeError("checkpoint chunk exceeds 64 MiB")
    npz_path = directory / f"{case_id}-{ordinal:06d}.npz"
    json_path = directory / f"{case_id}-{ordinal:06d}.json"
    if npz_path.exists() or json_path.exists():
        raise FileExistsError("refusing to overwrite immutable chunk " + npz_path.name)
    npz_tmp = directory / f".{npz_path.name}.{os.getpid()}.tmp"
    json_tmp = directory / f".{json_path.name}.{os.getpid()}.tmp"
    npz_tmp.write_bytes(payload)
    json_tmp.write_text(text)
    try:
        os.link(npz_tmp, npz_path)
        os.link(json_tmp, json_path)
    except FileExistsError:
        npz_tmp.unlink(missing_ok=True)
        json_tmp.unlink(missing_ok=True)
        raise FileExistsError("refusing to overwrite immutable chunk " + npz_path.name) from None
    npz_tmp.unlink()
    json_tmp.unlink()
    os.chmod(npz_path, 0o444)
    os.chmod(json_path, 0o444)
    body["npz"] = npz_path.name
    body["json"] = json_path.name
    body["npz_bytes"] = npz_path.stat().st_size
    body["json_bytes"] = json_path.stat().st_size
    return body


def _load_chunk(directory, case_id, ordinal):
    directory = Path(directory)
    json_path = directory / f"{case_id}-{ordinal:06d}.json"
    npz_path = directory / f"{case_id}-{ordinal:06d}.npz"
    record = json.loads(json_path.read_text())
    payload = npz_path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != record["arrays_sha256"]:
        raise ValueError("checkpoint bytes do not match the commit record")
    with np.load(npz_path, allow_pickle=False) as stored:
        arrays = {name: np.array(stored[name], copy=True) for name in stored.files}
    for name, digest in record["array_sha256"].items():
        if _array_sha(arrays[name]) != digest:
            raise ValueError("checkpoint array changed: " + name)
    return record, arrays


def _opening_arrays(loaded):
    state = loaded["state"]
    pair = loaded["pair"]
    return {
        "Q": state.Q, "r": state.r, "chi": state.chi,
        "pi_Q": state.p_Q, "pi_r": state.p_r, "pi_chi": state.p_chi,
        "phi0": state.phi0, "phi1": state.phi1,
        "W": np.asarray(pair.geometry_map),
        "source_phi0": np.asarray(pair.source_phi0),
        "source_phi1": np.asarray(pair.source_phi1),
        "observer_columns": np.asarray(pair.original_columns),
        "source_weights": np.asarray(pair.weights, dtype=float),
        "normal_clocks": np.asarray(loaded["clocks"], dtype=float),
        "clock_rates": clock_rates(pair, state),
    }


def _pair_record(loaded):
    pair = loaded["pair"]
    return {
        "nf": loaded["nf"],
        "coarse_indices": np.asarray(pair.geometry_coarse_indices, dtype=int).tolist(),
        "child_indices": np.asarray(pair.geometry_child_indices, dtype=int).tolist(),
        "parent_indices": np.asarray(pair.geometry_parent_indices, dtype=int).tolist(),
        "source_metadata": dict(pair.source_metadata),
        "geometry_metadata": dict(pair.geometry_metadata),
        "clock_locations": list(loaded["clock_locations"]),
    }


def _rebuild_pair(arrays, geometry_record):
    record = dict(geometry_record)
    record["nf"] = int(geometry_record["nf"])
    pair = episode.pair_from_arrays(arrays, record)
    state = model.NestedState(
        arrays["Q"], arrays["r"], arrays["chi"], arrays["pi_Q"], arrays["pi_r"], arrays["pi_chi"],
        arrays["phi0"], arrays["phi1"],
    )
    pair, state, info = episode.resolve_pair(pair, state, backend="fft")
    if info.get("backend") != "fft" or info.get("W_changed") or info.get("initial_state_called"):
        raise RuntimeError("replay carrier changed the stored pair")
    return pair, state


def probe_coupled_cpu(pair, state, dt):
    """One independent full step and one streamed step, on a copy."""
    started = time.process_time()
    compare_step(pair, state.copy(), float(dt))
    return max(time.process_time() - started, 1e-6)


def prepare_campaign(directory, *, fermions=(128, 256), duration=SEGMENT_DURATION,
                     step_cap=DEFAULT_STEP_CAP, cpu_budget_seconds=DEFAULT_CPU_BUDGET_SECONDS,
                     forecast_factor=DEFAULT_FORECAST_FACTOR):
    """Write the opening JSON+NPZ and the CPU forecast. Does not evolve the segment."""
    directory = _campaign_directory(directory, create=True)
    extent = segment_extent(duration, step_cap)
    remaining = float(cpu_budget_seconds)
    cases = []
    for fermions_count in fermions:
        loaded = load_matched_handoff(int(fermions_count))
        probe = probe_coupled_cpu(loaded["pair"], loaded["state"], extent["dt"])
        forecast = cpu_forecast(
            extent["steps"], probe, forecast_factor=forecast_factor, cpu_budget_seconds=remaining,
        )
        if forecast["admitted"]:
            remaining -= forecast["forecast_cpu_seconds"]
        case_id = f"nf{int(fermions_count)}_dt{_cap_token(extent['dt'])}"
        geometry = _pair_record(loaded)
        body = {
            "case_id": case_id,
            "nf": int(fermions_count),
            "parent_case": loaded["case"],
            "status": "prepared",
            "role": "opening",
            "steps": 0,
            "coordinate_time": HANDOFF_TIME,
            "time0": HANDOFF_TIME,
            "time1": extent["time1"],
            "duration": extent["duration"],
            "full_segment": extent["full_segment"],
            "dt": extent["dt"],
            "step_cap": extent["dt"],
            "pins": {key: loaded["pins"][key] for key in (
                "v1_json_sha256", "v1_npz_sha256", "basis_json_sha256", "basis_npz_sha256",
                "final_state_sha256", "source_columns_sha256", "observer_columns_sha256",
                "W_sha256", "weights_sha256", "phi_sha256",
            )},
            "forecast": forecast,
            "dephasing_prerequisite": False,
            "executor": "root",
            "initial_state_called": False,
            "W_rebuilt": False,
            "child_columns": list(CHILD_COLUMN_INDICES),
            "clock_locations": loaded["clock_locations"],
            "geometry_record": geometry,
            "code_sha256": _dependency_paths(),
        }
        committed = _commit_chunk(directory, case_id, 0, _opening_arrays(loaded), body)
        cases.append({
            "case_id": case_id,
            "nf": int(fermions_count),
            "status": "prepared",
            "opening": committed["npz"],
            "opening_sha256": committed["arrays_sha256"],
            "forecast": forecast,
            "parent_case": loaded["case"],
            "pins": body["pins"],
            "geometry_record": geometry,
        })
    demand = float(sum(case["forecast"]["forecast_cpu_seconds"] for case in cases))
    manifest = {
        "schema": RUN_SCHEMA,
        "production_segment": list(PRODUCTION_SEGMENT),
        "duration": extent["duration"],
        "time0": HANDOFF_TIME,
        "time1": extent["time1"],
        "full_segment": extent["full_segment"],
        "step_cap": extent["dt"],
        "steps": extent["steps"],
        "cpu_budget_seconds": float(cpu_budget_seconds),
        "forecast_factor": float(forecast_factor),
        "forecast_cpu_seconds": demand,
        "admitted": bool(demand <= float(cpu_budget_seconds) + 1e-9),
        "checkpoint_limit_bytes": CHUNK_LIMIT_BYTES,
        "production_evolved": False,
        "t3_in_scope": False,
        "dephasing_prerequisite_for_coupled_response": False,
        "executor": "root",
        "code_sha256": _dependency_paths(),
        "cases": cases,
        "status": "prepared",
    }
    _write_manifest(directory, manifest)
    return manifest


def _state_arrays(state, prefix):
    mapping = {"Q": state.Q, "r": state.r, "chi": state.chi, "pi_Q": state.p_Q,
               "pi_r": state.p_r, "pi_chi": state.p_chi, "phi0": state.phi0, "phi1": state.phi1}
    return {prefix + name: value for name, value in mapping.items()}


def _store_endpoint(opening_arrays, pair, full, split, clocks_full, clocks_split, endpoint):
    stored = {
        "W": opening_arrays["W"],
        "source_phi0": opening_arrays["source_phi0"],
        "source_phi1": opening_arrays["source_phi1"],
        "observer_columns": opening_arrays["observer_columns"],
        "source_weights": opening_arrays["source_weights"],
        "full_normal_clocks": np.asarray(clocks_full, dtype=float),
        "split_normal_clocks": np.asarray(clocks_split, dtype=float),
        "full_clock_rates": clock_rates(pair, full),
        "split_clock_rates": clock_rates(pair, nested_from_split(pair, split)),
    }
    stored.update(_state_arrays(full, "full_"))
    stored.update(_state_arrays(nested_from_split(pair, split), "split_"))
    for label in ("full", "complementary", "cross", "effective_outside"):
        for name in ("force_L", "force_Q", "force_beta"):
            stored[f"full_{label}_{name}"] = endpoint["full_partition"][name][label]
            stored[f"split_{label}_{name}"] = endpoint["split_partition"][name][label]
    return stored


def run_campaign(directory, *, cpu_budget_seconds=None):
    """Evolve a prepared campaign up to its forecast and the CPU budget.

    Root is the intended caller for the full 0.7 segment. A forecast that
    exceeds the budget commits no evolved chunk.
    """
    directory = _campaign_directory(directory, create=False)
    manifest = _read_manifest(directory)
    budget = float(manifest["cpu_budget_seconds"] if cpu_budget_seconds is None else cpu_budget_seconds)
    spent = 0.0
    evolved = []
    for case in manifest["cases"]:
        demand = float(case["forecast"]["forecast_cpu_seconds"])
        case_id = case["case_id"]
        if spent + demand > budget + 1e-9:
            case["status"] = "budget_stop"
            case["production_steps_taken"] = 0
            case["stop"] = "forecast exceeds the remaining CPU budget"
            continue
        opening, arrays = _load_chunk(directory, case_id, 0)
        pair, state = _rebuild_pair(arrays, opening["geometry_record"])
        sink = []
        started = time.process_time()
        comparison = compare_window(
            pair, state, opening["dt"], opening["duration"], arrays["normal_clocks"],
            process_budget_seconds=budget - spent, sink=sink,
        )
        used = time.process_time() - started
        spent += used
        held = sink[-1]
        comparison_public = {key: value for key, value in comparison.items() if not isinstance(value, np.ndarray)}
        if comparison["steps"] == 0 or held["endpoint"] is None:
            case["status"] = comparison["status"]
            case["production_steps_taken"] = 0
            case["cpu_seconds"] = used
            case["comparison"] = comparison_public
            continue
        endpoint_record = dict(opening)
        endpoint_record.update({
            "role": "endpoint",
            "status": comparison["status"],
            "steps": comparison["steps"],
            "coordinate_time": comparison["time1"],
            "comparison": comparison_public,
            "cpu_seconds": used,
            "dephasing_prerequisite": False,
        })
        committed = _commit_chunk(
            directory, case_id, 1,
            _store_endpoint(arrays, pair, held["full"], held["split"], held["clocks_full"],
                            held["clocks_split"], held["endpoint"]),
            endpoint_record,
        )
        case["status"] = comparison["status"]
        case["endpoint"] = committed["npz"]
        case["endpoint_sha256"] = committed["arrays_sha256"]
        case["production_steps_taken"] = comparison["steps"]
        case["coordinate_time"] = comparison["time1"]
        case["cpu_seconds"] = used
        case["comparison"] = comparison_public
        evolved.append(
            comparison["status"] == "completed"
            and abs(float(opening["duration"]) - SEGMENT_DURATION) <= 1e-12
            and abs(comparison["time1"] - PRODUCTION_SEGMENT[1]) <= 1e-12
        )
    manifest["cpu_seconds"] = spent
    manifest["production_evolved"] = bool(evolved) and all(evolved)
    manifest["status"] = "completed" if manifest["production_evolved"] else "incomplete"
    manifest["dephasing_prerequisite_for_coupled_response"] = False
    _write_manifest(directory, manifest)
    return manifest


def _replay_gap(stored, recomputed):
    return float(np.max(np.abs(np.asarray(stored) - np.asarray(recomputed))))


def check_campaign(directory):
    """Reload committed bytes and replay forces from them. Takes no time step."""
    directory = _campaign_directory(directory, create=False)
    manifest = _read_manifest(directory)
    replay_force = 0.0
    replay_clock = 0.0
    replay_partial = 0.0
    chunks = 0
    for case in manifest["cases"]:
        case_id = case["case_id"]
        opening, arrays = _load_chunk(directory, case_id, 0)
        if int(opening["steps"]) != 0:
            raise ValueError("opening chunk records evolved steps")
        for suffix in ("npz", "json"):
            path = directory / f"{case_id}-000000.{suffix}"
            if (path.stat().st_mode & 0o777) != 0o444:
                raise ValueError("opening chunk is not read-only")
        handoff = load_matched_handoff(int(case["nf"]))
        if not np.array_equal(arrays["W"], np.asarray(handoff["pair"].geometry_map)):
            raise ValueError("opening W is not the frozen basis")
        if not np.array_equal(arrays["phi0"], handoff["state"].phi0) or not np.array_equal(arrays["phi1"], handoff["state"].phi1):
            raise ValueError("opening source is not the saved T=0.3 state")
        if _array_sha(arrays["W"]) != case["pins"]["W_sha256"]:
            raise ValueError("opening W hash left the pin")
        if _array_sha(np.vstack((arrays["phi0"], arrays["phi1"]))) != case["pins"]["phi_sha256"]:
            raise ValueError("opening source hash left the pin")
        chunks += 1
        endpoint_path = directory / f"{case_id}-000001.npz"
        if not endpoint_path.is_file():
            continue
        endpoint, final = _load_chunk(directory, case_id, 1)
        if (endpoint_path.stat().st_mode & 0o777) != 0o444:
            raise ValueError("endpoint chunk is not read-only")
        pair, state = _rebuild_pair({
            "Q": final["full_Q"], "r": final["full_r"], "chi": final["full_chi"],
            "pi_Q": final["full_pi_Q"], "pi_r": final["full_pi_r"], "pi_chi": final["full_pi_chi"],
            "phi0": final["full_phi0"], "phi1": final["full_phi1"],
            "W": final["W"], "source_phi0": final["source_phi0"], "source_phi1": final["source_phi1"],
            "observer_columns": final["observer_columns"], "source_weights": final["source_weights"],
        }, opening["geometry_record"])
        partials = same_action_partials(pair, state)
        for name in ("force_L", "force_Q", "force_beta"):
            replay_force = max(replay_force, _replay_gap(final[f"full_full_{name}"], partials[name]))
        replay_clock = max(replay_clock, _replay_gap(final["full_clock_rates"], clock_rates(pair, state)))
        recorded = endpoint["comparison"]
        replay_partial = max(
            replay_partial,
            abs(partials["chart_gap"] - recorded["endpoint_chart_gap"]),
            abs(partials["lapse_partial_gap"] - recorded["endpoint_lapse_partial_gap"]),
            abs(partials["shift_partial_gap"] - recorded["endpoint_shift_partial_gap"]),
        )
        field_gap = max(_replay_gap(final["full_phi0"], final["split_phi0"]),
                        _replay_gap(final["full_phi1"], final["split_phi1"]))
        if field_gap > float(recorded["field_max_abs"]) + 1e-12:
            raise ValueError("stored field gap exceeds the committed comparison")
        geometry_gap = max(
            _replay_gap(final["full_" + name], final["split_" + name])
            for name in ("Q", "r", "chi", "pi_Q", "pi_r", "pi_chi")
        )
        if geometry_gap > float(recorded["geometry_max_abs"]) + 1e-12:
            raise ValueError("stored geometry gap exceeds the committed comparison")
        outside_gap = max(
            _replay_gap(final[f"full_effective_outside_{name}"], final[f"split_effective_outside_{name}"])
            for name in ("force_L", "force_Q", "force_beta")
        )
        if outside_gap > float(recorded["effective_outside_force_max_abs"]) + 1e-12:
            raise ValueError("stored outside-force gap exceeds the committed comparison")
        chunks += 1
    return {
        "stepped": False,
        "chunks": chunks,
        "replay_force_gap": replay_force,
        "replay_clock_gap": replay_clock,
        "replay_partial_gap": replay_partial,
        "production_evolved": bool(manifest.get("production_evolved")),
        "duration": manifest.get("duration"),
        "forecast_cpu_seconds": manifest.get("forecast_cpu_seconds"),
        "cpu_budget_seconds": manifest.get("cpu_budget_seconds"),
        "admitted": manifest.get("admitted"),
        "checkpoint_limit_bytes": manifest.get("checkpoint_limit_bytes"),
        "dephasing_prerequisite_for_coupled_response": False,
    }


def dependency_closure():
    """Pin the stored inputs and the segment configuration. Does not evolve."""
    cases = []
    for fermions in (128, 256):
        loaded = load_matched_handoff(fermions)
        pins = dict(loaded["pins"])
        cases.append({
            "case": loaded["case"],
            "nf": loaded["nf"],
            "time": loaded["time"],
            "dt": loaded["dt"],
            "backend": loaded["backend"],
            "fft_applied": loaded["fft_applied"],
            "initial_state_called": False,
            "W_rebuilt": False,
            "child_columns": list(CHILD_COLUMN_INDICES),
            "clock_locations": loaded["clock_locations"],
            "pins": {key: pins[key] for key in (
                "v1_json_sha256", "v1_npz_sha256", "basis_json_sha256", "basis_npz_sha256",
                "final_state_sha256", "source_columns_sha256", "observer_columns_sha256",
                "W_sha256", "weights_sha256", "phi_sha256",
            )},
        })
    return {
        "schema": SCHEMA,
        "production_segment": list(PRODUCTION_SEGMENT),
        "segment_duration": SEGMENT_DURATION,
        "allowed_step_caps": list(ALLOWED_STEP_CAPS),
        "checkpoint_limit_bytes": CHUNK_LIMIT_BYTES,
        "cpu_budget_seconds": DEFAULT_CPU_BUDGET_SECONDS,
        "forecast_factor": DEFAULT_FORECAST_FACTOR,
        "executor": "root",
        "production_evolved": False,
        "t3_in_scope": False,
        "dephasing_prerequisite_for_coupled_response": False,
        "short_check_duration": SHORT_DURATION,
        "short_check_ran": False,
        "cases": cases,
        "code_sha256": _dependency_paths(),
        "equations": dict(EXACT_EQUATIONS),
        "dephasing": dict(DEPHASING_GAP),
        "omission_is_unitary_global_law": False,
        "omission_is_autonomous_regeneration": False,
        "unphysical_control_gate_required": False,
        "new_interactions": False,
        "dense_covariance_stored": False,
        "dense_exterior_propagator_stored": False,
        "one_physical_metric": True,
        "multiplicity": "M=4 kappa once, inside the existing nodal forces",
    }


def short_coupled_check(fermions, duration=SHORT_DURATION):
    """Bounded comparison used by tests. The segment API is run_campaign."""
    if float(duration) > MAX_SHORT_DURATION:
        raise ValueError("the short check stops at 0.02; the full segment is run_campaign")
    loaded = load_matched_handoff(int(fermions))
    with backend.fft_thread_limit(1):
        window = compare_window(loaded["pair"], loaded["state"], loaded["dt"], float(duration), loaded["clocks"])
        identity = elimination_variation(loaded["pair"], loaded["state"])
        split = initial_split(loaded["pair"], loaded["state"])
        omissions = {name: omission_diagnostic(loaded["pair"], split, name) for name in ("memory", "drive")}
    window["case"] = loaded["case"]
    window["nf"] = loaded["nf"]
    window["time0"] = HANDOFF_TIME
    window["production_evolved"] = False
    window["elimination"] = {
        key: identity[key] for key in (
            "chart_gap", "lapse_partial_gap", "shift_partial_gap", "boundary_abs",
            "integration_by_parts_gap", "exterior_residual_integral_abs",
            "force_matches_reconstructed_source", "force_source_gap", "boundary_inserted_into_force",
            "second_induced_source", "laboratory_is_coupled_production",
        )
    }
    window["omissions"] = {
        name: {key: value[key] for key in (
            "residual_max_abs", "unitary_global_law", "autonomous_regeneration", "physical_gate",
        )}
        for name, value in omissions.items()
    }
    return window


def assert_output_outside_repository(path):
    resolved = Path(path).expanduser().resolve()
    if resolved == _REPO or _REPO in resolved.parents:
        raise PermissionError("refusing to write a closure record inside the repository")
    if resolved.exists():
        raise FileExistsError("refusing to overwrite " + str(resolved))
    return resolved


def _jsonable(value):
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        raise TypeError("closure records do not store raw arrays")
    if isinstance(value, (np.floating, float)):
        return float(value)
    if isinstance(value, (np.integer, int)) and not isinstance(value, bool):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if value is None or isinstance(value, str):
        return value
    raise TypeError("unsupported closure value " + type(value).__name__)


def closure_bytes(record):
    import json
    return json.dumps(_jsonable(record), indent=2, sort_keys=True) + "\n"

