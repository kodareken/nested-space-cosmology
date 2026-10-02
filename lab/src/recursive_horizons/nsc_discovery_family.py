"""Stage-4 source-family regime driver.

The question is whether a change in regional source shape or imbalance
moves contraction, dispersion, or retention off the baseline. The family
is five imbalances and three child-width scales. A bounded exploration
evolves six of those members. It does not open a second dynamical
system: the dense radius solve is ``initial_state``, the march from the
prepared slice to the handoff is ``propagate_baseline``, and the later
stations are the existing episode continuation.

T=0 Phi equals the prepared columns by value and is not an alias.
``nsc_discovery_episode._pins_unchanged`` treats that value equality as
an alias, so this driver does not hand the T=0 slice to that runner.
It evolves the prepared state first and hands the propagated field,
clocks, and source identity onward.
"""
from __future__ import annotations

from dataclasses import replace
import json
import math
from pathlib import Path
import time

import numpy as np

from . import nsc_discovery_episode as episode
from . import nsc_discovery_prediction as prediction
from . import nsc_discovery_regions as regions
from . import nsc_discovery_tidal as tidal
from . import nsc_nested_parent_child as model
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin


SCHEMA = "NSC-DISCOVERY-SOURCE-FAMILY-v1"
QUESTION = (
    "Does varying regional source shape or imbalance change contraction, "
    "dispersion, or retention beyond the baseline?"
)
IMBALANCES = (-1.0, -0.5, 0.0, 0.5, 1.0)
WIDTH_SCALES = (0.8, 1.0, 1.2)
BASELINE_HALF_WIDTH = 1.0
BASELINE_IMBALANCE = 1.0
BASELINE_WIDTH_SCALE = 1.0
OCCUPATION_TRACE = 3.0
EXPLORATION_NF = 128
CONFIRMATION_NF = 256
EXPLORATION_MEMBERS = (
    (1.0, 1.0),
    (0.0, 1.0),
    (-1.0, 1.0),
    (1.0, 0.8),
    (1.0, 1.2),
    (-1.0, 0.8),
)
# The owned positive homotopy stalls on this corner at nf=64 and nf=128.
# The column family still lists it. No substitute radius is constructed.
STALLED_MEMBER = (-1.0, 1.2)
CASE_COUNT = 6
THREADS_PER_CASE = 1
HANDOFF_TIME = float(episode.HANDOFF_TIME)
STATIONS = tuple(episode.DEFAULT_STATIONS)
STEP_CAP = float(episode.MATCHED_STEP_CAP)
NONPRODUCTION_DURATION = float(prediction.NONPRODUCTION_DURATION)
AGGREGATE_CPU_BUDGET_SECONDS = float(episode.DEFAULT_CPU_BUDGET_SECONDS)
FORECAST_FACTOR = float(episode.DEFAULT_FORECAST_FACTOR)
CHUNK_LIMIT_BYTES = int(episode.CHUNK_LIMIT_BYTES)
WORKERS = int(episode.DEFAULT_WORKERS)
SCALAR_FATE_PERCENT_GATE = None
UNIVERSAL_ONE_PERCENT_GATE = False
FAMILY_RECORD = "source-family.json"
EXPLICIT_ORTHOGONALIZATION = "outer_projected_against_child_then_lowdin"
CONSTRUCTOR = "nsc_nested_parent_child.initial_state"
EVOLUTION = "nsc_discovery_prediction.propagate_baseline"
CONTINUATION = "nsc_discovery_episode.run"

_MEMBER_ROLES = {
    (1.0, 1.0): "baseline",
    (0.0, 1.0): "uniform_imbalance",
    (-1.0, 1.0): "endpoint_imbalance",
    (1.0, 0.8): "narrow_shape",
    (1.0, 1.2): "wide_shape",
    (-1.0, 0.8): "crossed_endpoint_and_narrow",
}


def specification():
    """Family contract. Production evolution stays outside this record."""
    return {
        "schema": SCHEMA,
        "question": QUESTION,
        "imbalances": list(IMBALANCES),
        "width_scales": list(WIDTH_SCALES),
        "baseline_imbalance": BASELINE_IMBALANCE,
        "baseline_width_scale": BASELINE_WIDTH_SCALE,
        "occupation_trace": OCCUPATION_TRACE,
        "exploration_nf": EXPLORATION_NF,
        "confirmation_nf": CONFIRMATION_NF,
        "exploration_cases": CASE_COUNT,
        "threads_per_case": THREADS_PER_CASE,
        "workers": WORKERS,
        "coordinator_pools": 1,
        "aggregate_cpu_budget_seconds": AGGREGATE_CPU_BUDGET_SECONDS,
        "budget_covers": "initial_solves_handoff_measurements_and_continuation_child_cpu",
        "chunk_limit_bytes": CHUNK_LIMIT_BYTES,
        "forecast_factor": FORECAST_FACTOR,
        "stations": list(STATIONS),
        "handoff_time": HANDOFF_TIME,
        "step_cap": STEP_CAP,
        "nonproduction_duration": NONPRODUCTION_DURATION,
        "constructor": CONSTRUCTOR,
        "evolution": EVOLUTION,
        "continuation": CONTINUATION,
        "scalar_fate_percent_gate": SCALAR_FATE_PERCENT_GATE,
        "universal_one_percent_gate": UNIVERSAL_ONE_PERCENT_GATE,
        "forced_regeneration": False,
        "saved_t03_radius_copied": False,
        "third_dynamical_framework": False,
    }


def _match(value, choices, label):
    value = float(value)
    for item in choices:
        if abs(value - float(item)) <= 1e-15:
            return float(item)
    raise ValueError(label)


def canonical_imbalance(value):
    return _match(value, IMBALANCES, "imbalances are -1, -0.5, 0, 0.5 and 1")


def canonical_width_scale(value):
    return _match(value, WIDTH_SCALES, "child-width scales are 0.8, 1 and 1.2")


def half_width(width_scale):
    """Child half-width ``w``. Support is ``[2-w, 2+w]``."""
    scale = canonical_width_scale(width_scale)
    width = BASELINE_HALF_WIDTH * scale
    if not 0.0 < width < 2.0:
        raise ValueError("child half-width must leave positive outer annuli inside [0, 4]")
    return float(width)


def source_supports(width_scale):
    width = half_width(width_scale)
    return (
        (0.0, 2.0 - width),
        (2.0 - width, 2.0 + width),
        (2.0 + width, 4.0),
    )


def occupation_weights(imbalance):
    """Six CAR weights. Trace stays 3. Endpoints and the uniform vector lie in ``(0, 1]``."""
    slope = canonical_imbalance(imbalance)
    left = 0.5 + 0.25 * slope
    child = 0.5
    right = 0.5 - 0.25 * slope
    weights = np.array([left, left, child, child, right, right], dtype=float)
    if weights.shape != (6,) or not np.isfinite(weights).all():
        raise ValueError("source weights must be six finite numbers")
    if abs(float(np.sum(weights)) - OCCUPATION_TRACE) > 1e-12:
        raise ValueError("occupation trace left 3")
    if np.any(weights <= 0.0) or np.any(weights > 1.0):
        raise ValueError("source weights left (0, 1]")
    return weights


def regional_occupation_gap(weights):
    """Left-packet occupation minus right-packet occupation. For this family it equals the imbalance."""
    values = np.asarray(weights, dtype=float)
    if values.shape != (6,):
        raise ValueError("occupation gap expects six weights")
    return float(np.sum(values[:2]) - np.sum(values[4:]))


def _number_token(value):
    text = f"{abs(float(value)):.1f}".rstrip("0").rstrip(".")
    text = text.replace(".", "p") or "0"
    return ("m" if float(value) < 0.0 else "") + text


def case_id(imbalance, width_scale, nf):
    return f"nf{int(nf)}_a{_number_token(imbalance)}_w{_number_token(width_scale)}"


def member_role(imbalance, width_scale):
    key = (canonical_imbalance(imbalance), canonical_width_scale(width_scale))
    return _MEMBER_ROLES.get(key, "catalog_only")


def _readonly(value, dtype=None):
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def _packet_pair(xi, spacing, left, right, norms, phase):
    lobe_width = (float(right) - float(left)) / 2.0
    if lobe_width <= 0.0:
        raise ValueError("source lobe width must be positive")
    even_u = (xi - left) / lobe_width
    odd_u = (xi - left - lobe_width) / lobe_width
    even = coupling._bump(even_u) / (math.sqrt(lobe_width) * norms[0])
    odd = (odd_u - 0.5) * coupling._bump(odd_u) / (math.sqrt(lobe_width) * norms[1])
    even = even * math.sqrt(spacing)
    odd = odd * math.sqrt(spacing)
    if min(np.linalg.norm(even), np.linalg.norm(odd)) <= 1e-12:
        raise ValueError("separated source lobe is unresolved")
    envelope = (even / np.linalg.norm(even) + odd / np.linalg.norm(odd)) / math.sqrt(2.0)
    carrier = np.exp(1j * coupling.CARRIER_K * (xi - left)) * envelope
    plus = np.concatenate((carrier, 1j * carrier)) / math.sqrt(2.0)
    minus = phase * np.conjugate(plus)
    return plus, minus


def _tails(grid, phi0, phi1, supports):
    fine0 = grid.U_f @ phi0
    fine1 = grid.U_f @ phi1
    power = np.abs(fine0) ** 2 + np.abs(fine1) ** 2
    leakage = []
    for packet, (left, right) in enumerate(supports):
        outside = (grid.xi_q < left) | (grid.xi_q > right)
        for column in (2 * packet, 2 * packet + 1):
            total = float(np.sum(power[:, column]))
            if total <= 0.0:
                raise ValueError("source column has no fine-grid power")
            leakage.append(float(np.sum(power[outside, column]) / total))
    return leakage


def owned_separated_columns(grid, width_scale=1.0):
    """Bump/lobe packets on ``[0, 2-w]``, ``[2-w, 2+w]`` and ``[2+w, 4]``.

    Outer columns are projected off the child and Löwdin-orthonormalized.
    The child frame is not part of that repair. At unit scale the child
    columns are the owned original middle pair. At the other scales the
    child lobes are an explicit declared frame, then left unchanged.
    """
    scale = canonical_width_scale(width_scale)
    supports = source_supports(scale)
    norms = coupling._lobe_norms()
    phase = np.exp(1j * coupling.CALIBRATION["phase"])
    packets = []
    for left, right in supports:
        plus, minus = _packet_pair(grid.xi_f, grid.dx_f, left, right, norms, phase)
        packets.extend((plus, minus))
    raw = np.column_stack(packets)
    reference0, reference1, _preparation = coupling.prepare_rank6(grid.xi_f, grid.dx_f)
    reference_child = np.vstack((reference0[:, 2:4], reference1[:, 2:4]))
    if scale == BASELINE_WIDTH_SCALE:
        child = np.array(reference_child, copy=True)
        origin = "owned_original_middle"
    else:
        child = np.array(raw[:, 2:4], copy=True)
        origin = "declared_width_scaled_lobes"
    child_gram = child.conj().T @ child
    child_gram_gap = float(np.max(np.abs(child_gram - np.eye(2))))
    child_lowdin = False
    if child_gram_gap > 1e-9:
        if origin == "owned_original_middle":
            raise ValueError("owned original child is not orthonormal; refusing to repair it")
        child, _child_gram_before = coupling.lowdin(child)
        child_lowdin = True
    child_before = np.array(child, copy=True)
    outer = np.column_stack((raw[:, 0:2], raw[:, 4:6]))
    overlap = child.conj().T @ outer
    projected = outer - child @ overlap
    outer_columns, outer_gram = coupling.lowdin(projected)
    if not np.array_equal(child, child_before):
        raise RuntimeError("orthogonalization altered the child frame")
    columns = np.column_stack((outer_columns[:, 0:2], child, outer_columns[:, 2:4]))
    gram = columns.conj().T @ columns
    gram_gap = float(np.max(np.abs(gram - np.eye(6))))
    if gram_gap > 1e-9:
        raise ValueError("source columns must be orthonormal, without repairing the child frame")
    overlap_after = child.conj().T @ np.column_stack((columns[:, 0:2], columns[:, 4:6]))
    phi0 = np.array(columns[:grid.nf], copy=True)
    phi1 = np.array(columns[grid.nf:], copy=True)
    tails = _tails(grid, phi0, phi1, supports)
    mean_current = _probe_current_mean(grid, phi0, phi1)
    metadata = {
        "layout": "owned_separated_columns",
        "width_scale": scale,
        "half_width": half_width(scale),
        "source_support_closures": [list(item) for item in supports],
        "supports": [list(item) for item in supports],
        "child_support": list(supports[1]),
        "outer_annuli": [list(supports[0]), list(supports[2])],
        "child_frame_origin": origin,
        "child_frame_altered": False,
        "child_frame_silently_altered": False,
        "child_pair_lowdin_explicit": child_lowdin,
        "child_gram_before_explicit_lowdin": child_gram_gap,
        "middle_columns_equal_owned_original": bool(np.array_equal(child, reference_child)),
        "explicit_orthogonalization": EXPLICIT_ORTHOGONALIZATION,
        "child_included_in_lowdin": False,
        "outer_projected_against_fixed_child": True,
        "overlap_before_max": float(np.max(np.abs(overlap))),
        "overlap_after_max": float(np.max(np.abs(overlap_after))),
        "outer_gram_before_lowdin_max": float(np.max(np.abs(outer_gram - np.eye(4)))),
        "gram_max": gram_gap,
        "tails": tails,
        "tail_max": float(max(tails)),
        "fine_outside_support_power_fraction": tails,
        "exact_spatial_support": False,
        "mean_current": mean_current,
        "mean_current_deleted": False,
        "historical_source_reused_as_trajectory": False,
        "phase_convention": "exp(i*recorded_phase) times conjugate of plus columns",
        "scaled_lobe_half_density": True,
        "complement_occupation": 0.0,
        "physical_vacuum_identification": False,
    }
    return phi0, phi1, metadata


def _probe_current_mean(grid, phi0, phi1):
    nodal = galerkin.blank_state(grid, phi0, phi1)
    fine = galerkin.prolong_state(grid, nodal)
    active = galerkin.active_fine_system(grid, fine)
    source = coupling.source_from_columns(active, fine)
    return float(np.mean(source["force_beta"] / grid.dx_q))


def audited_chart(grid):
    """The locked ``A``, ``C_W``, flux and kappa. This family does not replace them."""
    coefficients = coupling.locked_coefficients()
    fine = grid.fine
    if float(fine.A) != float(coefficients["A"]) or float(fine.C_W) != float(coefficients["C_W"]):
        raise ValueError("A or C_W left the audited chart")
    if float(fine.flux) != float(coefficients["flux"]) or int(fine.kappa) != int(coupling.KAPPA):
        raise ValueError("flux or kappa left the audited chart")
    return {
        "A": float(coefficients["A"]),
        "C_W": float(coefficients["C_W"]),
        "flux": float(coefficients["flux"]),
        "kappa": int(coupling.KAPPA),
    }


def assert_car(phi0, phi1, weights):
    """CAR stays admissible: weights in ``(0, 1]``, eigenvalues in ``[0, 1]``, trace 3."""
    values = np.asarray(weights, dtype=float)
    if values.shape != (6,) or abs(float(np.sum(values)) - OCCUPATION_TRACE) > 1e-12:
        raise ValueError("CAR occupation trace left 3")
    return prediction.assert_admissible_source(phi0, phi1, values)


def load_frozen_basis(nf):
    """Frozen portable ``W`` for nf 128 or 256. The array is the complete geometry basis."""
    nf = int(nf)
    if nf not in (EXPLORATION_NF, CONFIRMATION_NF):
        raise ValueError("frozen portable baseline W is stored for nf 128 and nf 256")
    basis = json.loads(episode.BASIS_JSON.read_text())
    payload = episode.file_sha256(episode.BASIS_NPZ)
    if payload != basis["payload_sha256"]:
        raise ValueError("frozen replay basis payload does not match its record")
    entry = basis["bases"][str(nf)]
    with np.load(episode.BASIS_NPZ, allow_pickle=False) as frozen:
        matrix = np.array(frozen[f"nf{nf}_W"], copy=True)
        observer = np.array(frozen[f"nf{nf}_reference_columns"], copy=True)
        saved_source = np.array(frozen[f"nf{nf}_source_columns"], copy=True)
    expected = nf - 1
    if matrix.shape != (expected, expected):
        raise ValueError("frozen W is not the complete geometry basis")
    gap = float(np.max(np.abs(matrix.T @ matrix - np.eye(expected))))
    if gap > 1e-8:
        raise ValueError("frozen W is not orthonormal")
    return {
        "nf": nf,
        "W": matrix,
        "observer_columns": observer,
        "saved_source_columns": saved_source,
        "coarse_indices": [int(value) for value in entry["coarse_indices"]],
        "child_indices": [int(value) for value in entry["child_indices"]],
        "parent_indices": [int(value) for value in entry["parent_indices"]],
        "geometry_metadata": dict(entry["geometry"]),
        "orthogonality_max": gap,
        "complete": True,
        "basis": "frozen_portable_baseline",
        "basis_json_sha256": episode.file_sha256(episode.BASIS_JSON),
        "basis_npz_sha256": payload,
        "W_sha256": episode.array_sha256(matrix),
        "observer_sha256": episode.array_sha256(observer),
    }


def _owned_frame(grid):
    matrix, coarse, child, parent, geometry_meta = model._geometry_frame(grid)
    return {
        "nf": int(grid.nf),
        "W": np.array(matrix, copy=True),
        "coarse_indices": [int(value) for value in coarse],
        "child_indices": [int(value) for value in child],
        "parent_indices": [int(value) for value in parent],
        "geometry_metadata": dict(geometry_meta),
        "orthogonality_max": float(geometry_meta["orthogonality_max"]),
        "complete": True,
        "basis": "owned_geometry_frame",
        "frozen_basis_authenticated": False,
    }


def assemble_pair(nf, imbalance, width_scale, *, basis="frozen"):
    """Source columns, explicit observer, and a complete ``W``. No radius solve."""
    nf = int(nf)
    if nf < 32 or nf % 2:
        raise ValueError("source family pair needs an even nf of at least 32")
    weights = occupation_weights(imbalance)
    scale = canonical_width_scale(width_scale)
    grid = galerkin.build_grid(nf, gauge="conformal")
    grid.fine = replace(grid.fine, occupations=weights.copy())
    chart = audited_chart(grid)
    phi0, phi1, metadata = owned_separated_columns(grid, scale)
    eigenvalues = assert_car(phi0, phi1, weights)
    reference0, reference1, _preparation = coupling.prepare_rank6(grid.xi_f, grid.dx_f)
    reference = np.vstack((reference0, reference1))
    columns = np.vstack((phi0, phi1))
    if basis == "frozen":
        frozen = load_frozen_basis(nf)
        matrix = frozen["W"]
        if matrix.shape[0] != grid.ng:
            raise ValueError("frozen W does not match this conformal geometry degree")
        coarse = frozen["coarse_indices"]
        child = frozen["child_indices"]
        parent = frozen["parent_indices"]
        geometry_meta = frozen["geometry_metadata"]
        observer = frozen["observer_columns"]
        if observer.shape != columns.shape:
            raise ValueError("frozen observer frame does not match the source shape")
        basis_name = frozen["basis"]
        basis_pins = {
            "basis_json_sha256": frozen["basis_json_sha256"],
            "basis_npz_sha256": frozen["basis_npz_sha256"],
            "frozen_basis_authenticated": True,
            "W_sha256": frozen["W_sha256"],
        }
    elif basis == "owned_geometry_frame":
        if nf in (EXPLORATION_NF, CONFIRMATION_NF):
            raise ValueError("nf 128 and nf 256 reuse the frozen portable W")
        frame = _owned_frame(grid)
        matrix = frame["W"]
        coarse = frame["coarse_indices"]
        child = frame["child_indices"]
        parent = frame["parent_indices"]
        geometry_meta = frame["geometry_metadata"]
        observer = reference
        basis_name = frame["basis"]
        basis_pins = {"frozen_basis_authenticated": False, "W_sha256": episode.array_sha256(matrix)}
    else:
        raise ValueError("basis must be frozen or owned_geometry_frame")
    metadata = dict(
        metadata,
        occupations=weights.tolist(),
        total_occupation=float(np.sum(weights)),
        occupation_trace=OCCUPATION_TRACE,
        regional_occupation_gap=regional_occupation_gap(weights),
        eigenvalue_min=float(np.min(eigenvalues)),
        eigenvalue_max=float(np.max(eigenvalues)),
        eigenvalues_admissible=True,
        reference_observer="explicit rank6 observer, distinct from the separated source",
        source_equals_reference=bool(np.array_equal(columns, observer)),
        child_source_equals_reference=bool(np.array_equal(columns[:, 2:4], observer[:, 2:4])),
        observer_is_explicit=True,
        basis=basis_name,
        audited=chart,
        saved_t03_radius_copied=False,
    )
    if metadata["source_equals_reference"]:
        raise ValueError("separated source collapsed onto the observer frame")
    pair = model.NestedPair(
        grid,
        _readonly(matrix),
        _readonly(coarse, int),
        _readonly(child, int),
        _readonly(parent, int),
        _readonly(observer[:nf]),
        _readonly(observer[nf:]),
        _readonly(observer),
        _readonly(phi0),
        _readonly(phi1),
        _readonly(columns),
        _readonly(weights),
        metadata,
        dict(geometry_meta),
        clock_locations=(1.0, 2.0, 3.0),
    )
    if tuple(int(index) for index in pair.child_indices) != (2, 3):
        raise RuntimeError("source child columns must stay (2, 3)")
    return pair, basis_pins


def solve_prepared(pair):
    """Source-consistent dense radius solve. The saved T=0.3 radius is not read."""
    if pair.grid.gauge != "conformal":
        raise ValueError("dense initial radius solve requires the conformal chart")
    chart = audited_chart(pair.grid)
    started = time.process_time()
    state, report = model.initial_state(pair, solve_constraints=True)
    elapsed = float(time.process_time() - started)
    if not report.get("converged"):
        blocker = report.get("blocker")
        raise ValueError("source-consistent dense initial radius solve did not converge: " + str(blocker))
    if report.get("imposed_radius") or report.get("filtered_v1_radius"):
        raise ValueError("initial radius was imposed or copied from a saved trajectory")
    if "retained" not in str(report.get("shift_method", "")):
        raise ValueError("full shift residual was not retained")
    source0 = np.asarray(pair.source_phi0)
    source1 = np.asarray(pair.source_phi1)
    if not np.array_equal(state.phi0, source0) or not np.array_equal(state.phi1, source1):
        raise ValueError("T0 Phi is not the prepared source")
    if state.phi0 is pair.source_phi0 or np.shares_memory(state.phi0, source0):
        raise ValueError("T0 Phi aliases the prepared source")
    if state.phi1 is pair.source_phi1 or np.shares_memory(state.phi1, source1):
        raise ValueError("T0 Phi aliases the prepared source")
    nodal = model.reconstruct_state(pair, state)
    fine = galerkin.prolong_state(pair.grid, nodal)
    active = galerkin.active_fine_system(pair.grid, fine)
    source = coupling.source_from_columns(active, fine)
    shift = coupling.shift_constraint(active, fine) + source["force_beta"] / pair.grid.dx_q
    report = dict(report, audited=chart, saved_t03_radius_copied=False,
                  initial_cpu_seconds=float(time.process_time() - started),
                  shift_residual_max=float(np.max(np.abs(shift))),
                  shift_residual_rms=float(np.sqrt(np.mean(np.abs(shift) ** 2))),
                  shift_residual_mean=float(np.mean(shift)))
    return state, report


def initial_record(pair, state, report):
    """Scalar initial data. Arrays stay on the state that is evolved forward."""
    diagnostics = report.get("solver_diagnostics") or {}
    weights = np.asarray(pair.weights, dtype=float)
    metadata = pair.source_metadata
    return {
        "case_id": case_id(regional_occupation_gap(weights), metadata["width_scale"], pair.grid.nf),
        "imbalance": regional_occupation_gap(weights),
        "width_scale": metadata["width_scale"],
        "half_width": metadata["half_width"],
        "nf": int(pair.grid.nf),
        "role": member_role(regional_occupation_gap(weights), metadata["width_scale"]),
        "weights": weights.tolist(),
        "occupation_trace": float(np.sum(weights)),
        "regional_occupation_gap": regional_occupation_gap(weights),
        "supports": metadata["supports"],
        "gram_max": metadata["gram_max"],
        "tails": list(metadata["tails"]),
        "tail_max": metadata["tail_max"],
        "overlap_before_max": metadata["overlap_before_max"],
        "overlap_after_max": metadata["overlap_after_max"],
        "explicit_orthogonalization": metadata["explicit_orthogonalization"],
        "child_frame_origin": metadata["child_frame_origin"],
        "child_frame_altered": False,
        "mean_current": metadata["mean_current"],
        "mean_current_deleted": False,
        "eigenvalue_min": metadata["eigenvalue_min"],
        "eigenvalue_max": metadata["eigenvalue_max"],
        "converged": True,
        "method": report.get("method"),
        "shift_method": report.get("shift_method"),
        "shift_residual_mean": report.get("shift_residual_mean"),
        "shift_residual_max": report.get("shift_residual_max"),
        "shift_residual_rms": report.get("shift_residual_rms"),
        "source_current_mean": report.get("source_current_mean"),
        "imposed_radius": False,
        "filtered_v1_radius": False,
        "saved_t03_radius_copied": False,
        "full_hamilton_max": diagnostics.get("full_hamilton_max"),
        "held_out_hamilton_max": diagnostics.get("held_out_hamilton_max"),
        "projected_hamilton_max": diagnostics.get("projected_hamilton_max"),
        "r_min": float(np.min(state.r)),
        "r_max": float(np.max(state.r)),
        "state_sha256": episode.state_sha256(state),
        "phi_equals_prepared_columns": True,
        "phi_is_alias": False,
        "audited": dict(metadata["audited"]),
        "basis": metadata["basis"],
        "initial_cpu_seconds": report.get("initial_cpu_seconds"),
        "constructor": CONSTRUCTOR,
    }


def prepare_member(imbalance, width_scale, *, nf, basis="frozen"):
    """Build one geometry and solve its own dense initial radius."""
    pair, basis_pins = assemble_pair(nf, imbalance, width_scale, basis=basis)
    state, report = solve_prepared(pair)
    record = initial_record(pair, state, report)
    record["basis_pins"] = {key: value for key, value in basis_pins.items() if key != "W"}
    return {"pair": pair, "state": state, "report": report, "initial": record, "basis_pins": basis_pins}


def evolve_to_handoff(pair, state, *, duration, step_cap, production):
    """Existing production march. The returned field is the handoff, not a reset to T=0."""
    decision = prediction.accepted_duration(duration, production=production)
    if decision["production"] and abs(float(duration) - HANDOFF_TIME) > 1e-12:
        raise ValueError("source-family production handoff is T=0.3")
    started = time.process_time()
    rates_before = np.asarray(model.metrics(pair, state)["clock_rates"], dtype=float)
    propagated = prediction.propagate_baseline(
        pair,
        state,
        duration=float(decision["duration"]),
        step_cap=float(step_cap),
        production=bool(decision["production"]),
        backend_name=prediction.EVOLUTION_BACKEND,
    )
    evolved = propagated["state"]
    carrier = propagated["pair"]
    rates_after = np.asarray(model.metrics(carrier, evolved)["clock_rates"], dtype=float)
    span = float(propagated["coordinate_time"])
    wings = 0.5 * span * (rates_before + rates_after)
    normal = np.array(wings, dtype=float, copy=True)
    normal[1] = float(propagated["tau"])
    source = np.vstack((np.asarray(pair.source_phi0), np.asarray(pair.source_phi1)))
    phi = np.vstack((np.asarray(evolved.phi0), np.asarray(evolved.phi1)))
    same = bool(np.array_equal(phi, source))
    if same:
        raise RuntimeError(
            "propagate_baseline returned Phi equal to the prepared columns. "
            "nsc_discovery_episode._pins_unchanged treats that value equality as an alias "
            "('state Phi is aliased to the source columns'). "
            "Refusing to hand this slice to the old runner."
        )
    clocks = {
        "locations": [1.0, 2.0, 3.0],
        "rates": rates_after,
        "normal_clocks": normal,
        "protocol": {
            "child_x2": "rk4_stages_from_propagate_baseline",
            "wings_x1_x3": "endpoint_trapezoid",
            "frozen_slice_product_used_as_increment": False,
        },
    }
    identity = {
        "evolution": EVOLUTION,
        "production": bool(decision["production"]),
        "target": decision["target"],
        "coordinate_time": span,
        "tau": float(propagated["tau"]),
        "steps": int(propagated["step_count"]),
        "t0_phi_equals_prepared_columns": True,
        "t0_phi_is_alias": False,
        "handoff_phi_equals_source": False,
        "handoff_is_propagate_baseline_output": True,
        "saved_t03_radius_copied": False,
        "state_sha256": episode.state_sha256(evolved),
        "source_columns_sha256": episode.array_sha256(source),
        "phi_sha256": episode.array_sha256(phi),
        "backend": propagated.get("backend"),
        "handoff_cpu_seconds": float(time.process_time() - started),
    }
    return propagated, clocks, identity


def _sign(value):
    if value is None:
        return None
    value = float(value)
    if not math.isfinite(value):
        return None
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def retention_order(fractions):
    """Descending rank of packet retention. No percent cutoff."""
    if fractions is None or len(fractions) != 3:
        return None
    values = [None if item is None else float(item) for item in fractions]
    if any(item is None or not math.isfinite(item) for item in values):
        return None
    order = tuple(int(index) for index in np.argsort(np.array(values, dtype=float))[::-1])
    return order


def physical_signals(*, length_change, width_change, retained_fraction, packet_retention):
    """Signs and a retention rank. Magnitudes are stored and do not choose the class."""
    contraction = _sign(None if length_change is None else -float(length_change))
    dispersion = _sign(width_change)
    order = retention_order(packet_retention)
    child_retention = None if retained_fraction is None else float(retained_fraction)
    unresolved = (contraction is None or dispersion is None or order is None
                  or child_retention is None or not math.isfinite(child_retention))
    return {
        "contraction_sign": contraction,
        "dispersion_sign": dispersion,
        "retention_order": None if order is None else list(order),
        "child_retained_fraction": child_retention,
        "length_change": None if length_change is None else float(length_change),
        "width_change": None if width_change is None else float(width_change),
        "unresolved": bool(unresolved),
        "percent_gate": None,
        "universal_one_percent_gate": False,
    }


def fate_class(signals, reference, *, is_baseline=False):
    """Physical sign pattern. A 1% scalar gate is not a class, and neither is a test count."""
    if is_baseline and not signals.get("unresolved"):
        return "baseline"
    if signals.get("unresolved") or reference is None or reference.get("unresolved"):
        return "unresolved"
    keys = ("contraction_sign", "dispersion_sign", "retention_order", "areal_radius_sign", "geometry_growth_signs", "packet_transfer_signs")
    same = all(signals.get(key) == reference.get(key) for key in keys)
    return "same_physical_signs" if same else "distinct_physical_signs"


def regime_observable(signals):
    """Child-source retention, reported raw."""
    return {
        "name": "child_source_retained_fraction",
        "value": None if signals is None else signals.get("child_retained_fraction"),
        "percent_gate": None,
        "universal_one_percent_gate": False,
        "discriminates_by": "retained fraction and the sign pattern of contraction, dispersion, retention rank",
    }


def measure_consumers(pair, state, clocks, time_value):
    """Regions and tidal read the state. A failure is recorded exactly and is not replaced."""
    blockers = []
    region_row = None
    tidal_row = None
    try:
        region_row = regions.analyze(pair, state, float(time_value), control_mode="coupled")
    except Exception as error:
        blockers.append({"consumer": "regions", "error": type(error).__name__ + ": " + str(error)})
    try:
        jets = tidal.analytic_accelerations(pair, state, "coupled")
        tidal_row = tidal.station_observables(pair, jets, clocks)
    except Exception as error:
        blockers.append({"consumer": "tidal", "error": type(error).__name__ + ": " + str(error)})
    retention = None
    width = None
    if region_row is not None:
        pairs = region_row["localization"]["column_pairs"]
        retention = [item.get("retained_fraction") for item in pairs]
        width = region_row["localization"].get("six_column_total_width")
        if width is None:
            width = region_row["localization"]["originally_child_prepared_pair"].get("proper_width")
    measured = model.metrics(pair, state)
    transfer = None
    if region_row is not None:
        fine, system = model._active(pair, state)
        masses = regions._column_mass(fine, system)
        transfer = []
        for columns in regions.COLUMN_PAIRS:
            mass = np.sum(masses[:, list(columns)], axis=1)
            content = float(np.sum(mass))
            transfer.append([
                None if content <= 0.0 else regions._window_probability(pair.grid, mass, support) / content
                for support in source_supports(pair.source_metadata["width_scale"])
            ])
    return {
        "coordinate_time": float(time_value),
        "child_proper_length": measured["child_proper_length"],
        "child_areal_radius": measured["child_r_proper_mean"],
        "metric_positive": measured["Q_min"] > 0.0 and measured["r_min"] > 0.0,
        "probability_positive": None if region_row is None else region_row["positive_probability"],
        "source_packet_window_fractions": transfer,
        "source_packet_window_fraction_convention": "rows are source pairs; columns are left, child, right windows",
        "actual_geometry": None if tidal_row is None else {
            "finite": tidal_row["finite"],
            "R4_l2": tidal_row["summaries"]["R4"]["l2"],
            "radial_tide_l2": tidal_row["summaries"]["R_0101"]["l2"],
            "angular_tide_l2": tidal_row["summaries"]["R_0202"]["l2"],
            "chi_substituted_for_curvature": False,
        },
        "regions_ready": region_row is not None,
        "tidal_ready": tidal_row is not None,
        "packet_retention": retention,
        "child_retained_fraction": None if not retention else retention[1],
        "packet_width": width,
        "blockers": blockers,
        "source_reset": False,
    }


def signals_from_samples(earlier, later):
    length_change = None
    width_change = None
    if earlier.get("child_proper_length") is not None and later.get("child_proper_length") is not None:
        length_change = float(later["child_proper_length"]) - float(earlier["child_proper_length"])
    if earlier.get("packet_width") is not None and later.get("packet_width") is not None:
        width_change = float(later["packet_width"]) - float(earlier["packet_width"])
    signals = physical_signals(
        length_change=length_change,
        width_change=width_change,
        retained_fraction=later.get("child_retained_fraction"),
        packet_retention=later.get("packet_retention"),
    )
    signals["areal_radius_change"] = (None if earlier.get("child_areal_radius") is None
        or later.get("child_areal_radius") is None else later["child_areal_radius"] - earlier["child_areal_radius"])
    growth = {}
    for name in ("R4_l2", "radial_tide_l2", "angular_tide_l2"):
        first = (earlier.get("actual_geometry") or {}).get(name)
        last = (later.get("actual_geometry") or {}).get(name)
        growth[name] = None if first is None or last is None else float(last) - float(first)
    transfer_first = earlier.get("source_packet_window_fractions")
    transfer_last = later.get("source_packet_window_fractions")
    signals["packet_transfer_change"] = None
    signals["packet_transfer_signs"] = None
    if transfer_first is not None and transfer_last is not None:
        delta = np.asarray(transfer_last, dtype=float) - np.asarray(transfer_first, dtype=float)
        signals["packet_transfer_change"] = delta.tolist()
        signals["packet_transfer_signs"] = [_sign(delta[i, j]) for i in range(3) for j in range(3) if i != j]
    signals["areal_radius_sign"] = _sign(signals["areal_radius_change"])
    signals["actual_geometry_growth"] = growth
    signals["geometry_growth_signs"] = {name: _sign(value) for name, value in growth.items()}
    return signals


def select_confirmation(rows):
    """nf=256 copies of members whose physical signs differ. Unresolved rows are not regenerated."""
    selected = []
    for row in rows or []:
        if row.get("fate") != "distinct_physical_signs":
            continue
        selected.append({
            "imbalance": row["imbalance"],
            "width_scale": row["width_scale"],
            "parent_case": row["case_id"],
            "nf": CONFIRMATION_NF,
            "fate": row["fate"],
            "forced": False,
        })
        if len(selected) == CASE_COUNT:
            break
    return {
        "nf": CONFIRMATION_NF,
        "cases": selected,
        "forced_regeneration": False,
        "reason": None if selected else "no distinct physical regime to confirm",
    }


def answer_question(rows):
    distinct = [row["case_id"] for row in rows or [] if row.get("fate") == "distinct_physical_signs"]
    unresolved = [row["case_id"] for row in rows or [] if row.get("fate") == "unresolved"]
    if not rows:
        return {"answered": False, "changes_beyond_baseline": None, "reason": "no measured regimes"}
    if distinct:
        return {
            "answered": True,
            "changes_beyond_baseline": True,
            "distinct_cases": distinct,
            "percent_gate": None,
        }
    if unresolved:
        return {
            "answered": False,
            "changes_beyond_baseline": None,
            "unresolved_cases": unresolved,
            "reason": "physical signals unresolved; confirmation was not forced",
            "percent_gate": None,
        }
    return {"answered": True, "changes_beyond_baseline": False, "distinct_cases": [], "percent_gate": None}


def campaign_forecast(n_cases=CASE_COUNT, nf=EXPLORATION_NF, step_cap=STEP_CAP):
    """Rate-sample forecast. Dense initial solves are counted and are not inside that sample."""
    n_cases = int(n_cases)
    nf = int(nf)
    step_cap = float(step_cap)
    if n_cases < 1 or step_cap <= 0.0:
        raise ValueError("forecast needs cases and a positive step cap")
    handoff_steps = int(math.ceil(HANDOFF_TIME / step_cap - 1e-12))
    continue_steps = int(math.ceil((STATIONS[-1] - HANDOFF_TIME) / step_cap - 1e-12))
    seconds_per_rate = (0.037 / 12.0) * (nf / 256.0) * (math.log2(nf) / math.log2(256.0))
    evolve = n_cases * handoff_steps * 4 * seconds_per_rate
    child = n_cases * continue_steps * 4 * seconds_per_rate
    forecast_child = FORECAST_FACTOR * child
    return {
        "n_cases": n_cases,
        "nf": nf,
        "step_cap": step_cap,
        "dense_initial_solves": n_cases,
        "dense_initial_solves_in_rate_sample": False,
        "propagate_baseline_cpu_seconds": float(evolve),
        "continuation_raw_cpu_seconds": float(child),
        "continuation_forecast_cpu_seconds": float(forecast_child),
        "forecast_cpu_seconds_without_inits": float(evolve + forecast_child),
        "aggregate_budget_seconds": AGGREGATE_CPU_BUDGET_SECONDS,
        "budget_covers": "initial_solves_handoff_measurements_and_continuation_child_cpu",
        "within_budget_ignoring_unsampled_inits": bool(evolve + forecast_child <= AGGREGATE_CPU_BUDGET_SECONDS),
        "hard_gate": False,
        "chunk_limit_bytes": CHUNK_LIMIT_BYTES,
        "threads_per_case": THREADS_PER_CASE,
        "workers": min(WORKERS, n_cases),
        "coordinator_pools": 1,
    }


def family_catalog(nf=EXPLORATION_NF):
    """All 15 admissible members. This does not solve or evolve."""
    members = []
    for imbalance in IMBALANCES:
        for scale in WIDTH_SCALES:
            weights = occupation_weights(imbalance)
            supports = source_supports(scale)
            members.append({
                "case_id": case_id(imbalance, scale, nf),
                "imbalance": imbalance,
                "width_scale": scale,
                "half_width": half_width(scale),
                "weights": weights.tolist(),
                "occupation_trace": OCCUPATION_TRACE,
                "regional_occupation_gap": regional_occupation_gap(weights),
                "supports": [list(item) for item in supports],
                "role": member_role(imbalance, scale),
                "nf": int(nf),
                "solved": False,
                "saved_t03_radius_copied": False,
                "weights_in_open_unit_interval": True,
            })
    exploration = [case_id(imbalance, scale, nf) for imbalance, scale in EXPLORATION_MEMBERS]
    if len(exploration) != CASE_COUNT or len(set(exploration)) != CASE_COUNT:
        raise RuntimeError("exploration batch is not six distinct cases")
    return members, exploration


def _family_path(directory):
    return Path(directory) / FAMILY_RECORD


def _write_json(path, record):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(episode.jsonable(record), indent=2, allow_nan=False) + "\n")


def _read_family(directory):
    path = _family_path(directory)
    if not path.is_file():
        raise FileNotFoundError("source-family record is missing: " + str(path))
    return json.loads(path.read_text())


def prepare_specification(output, *, nf=EXPLORATION_NF, all_members=False):
    """Write the family and the six-case admission list. Do not solve."""
    directory = episode.assert_campaign_output(output)
    if _family_path(directory).exists() or (directory / "manifest.json").exists():
        raise FileExistsError("source-family output already exists; refusing to overwrite")
    directory.mkdir(parents=True, exist_ok=True)
    members, exploration = family_catalog(nf)
    if all_members:
        exploration = [item["case_id"] for item in members]
    record = dict(specification())
    record.update(
        output=str(directory),
        stage="specification",
        materialized=False,
        evolved=False,
        nf=int(nf),
        family=members,
        cases=exploration,
        forecast=campaign_forecast(len(exploration), int(nf), STEP_CAP),
        all_members=bool(all_members),
        regimes=[],
        question_status=answer_question([]),
        confirmation=None,
        initial_cpu_seconds=0.0,
        child_cpu_seconds=0.0,
        executor_pools=0,
    )
    _write_json(_family_path(directory), record)
    return record


def prepare_confirmation(output, exploration):
    """A separate directory. Selected nf=256 members only. No forced regeneration."""
    destination = episode.assert_campaign_output(output)
    source = Path(exploration).expanduser().resolve()
    if destination == source or source in destination.parents or destination in source.parents:
        raise PermissionError("confirmation output must be a different directory from the exploration")
    if _family_path(destination).exists():
        raise FileExistsError("confirmation output already exists; refusing to overwrite")
    report = _read_family(source)
    selected = select_confirmation(report.get("regimes") or [])
    if selected["cases"]:
        selected["cases"] = selected["cases"][:CASE_COUNT - 1]
        selected["cases"].insert(0, {
            "imbalance": BASELINE_IMBALANCE, "width_scale": BASELINE_WIDTH_SCALE,
            "parent_case": case_id(BASELINE_IMBALANCE, BASELINE_WIDTH_SCALE, EXPLORATION_NF),
            "nf": CONFIRMATION_NF, "fate": "baseline_comparator", "forced": False,
        })
    destination.mkdir(parents=True, exist_ok=True)
    record = dict(specification())
    record.update(
        output=str(destination),
        exploration_directory=str(source),
        stage="confirmation_specification",
        materialized=False,
        evolved=False,
        nf=CONFIRMATION_NF,
        cases=[case_id(item["imbalance"], item["width_scale"], CONFIRMATION_NF) for item in selected["cases"]],
        selected=selected,
        forced_regeneration=False,
        forecast=campaign_forecast(max(1, len(selected["cases"])), CONFIRMATION_NF, STEP_CAP) if selected["cases"] else None,
        regimes=[],
        question_status=report.get("question_status"),
        initial_cpu_seconds=0.0,
        child_cpu_seconds=0.0,
        executor_pools=0,
    )
    _write_json(_family_path(destination), record)
    return record


def _handoff_record(item, stations):
    pair = item["pair"]
    identity = item["identity"]
    initial = item["initial"]
    pins = {
        "source_columns_sha256": identity["source_columns_sha256"],
        "observer_columns_sha256": episode.array_sha256(np.asarray(pair.reference_columns)),
        "W_sha256": item["basis_pins"]["W_sha256"],
        "weights_sha256": episode.array_sha256(np.asarray(pair.weights, dtype=float)),
        "phi_sha256": identity["phi_sha256"],
    }
    if item["basis_pins"].get("frozen_basis_authenticated"):
        pins["basis_json_sha256"] = item["basis_pins"]["basis_json_sha256"]
        pins["basis_npz_sha256"] = item["basis_pins"]["basis_npz_sha256"]
    if pins["phi_sha256"] == pins["source_columns_sha256"]:
        raise RuntimeError(
            "handoff Phi hash equals the source hash. "
            "The episode runner would flag that equality as an alias."
        )
    return {
        "case_id": initial["case_id"],
        "parent_case": initial["case_id"],
        "nf": int(pair.grid.nf),
        "control": "coupled",
        "control_mode": "coupled",
        "geometry": "evolving",
        "step_cap": float(item["step_cap"]),
        "stations": list(stations),
        "coordinate_time": float(identity["coordinate_time"]),
        "steps": 0,
        "stations_reached": [],
        "handoff_time": float(identity["coordinate_time"]),
        "phi": "propagate_baseline output; T0 Phi equalled prepared columns by value and was not an alias",
        "observer": "explicit rank6 observer",
        "momentum_representation": episode.CANONICAL_PI,
        "source_pins": dict(pins),
        "source_pins_before": dict(pins),
        "source_pins_after": dict(pins),
        "source_identity": dict(identity),
        "coarse_indices": [int(value) for value in pair.geometry_coarse_indices],
        "child_indices": [int(value) for value in pair.geometry_child_indices],
        "parent_indices": [int(value) for value in pair.geometry_parent_indices],
        "source_metadata": dict(pair.source_metadata),
        "geometry_metadata": dict(pair.geometry_metadata),
        "clock_locations": [1.0, 2.0, 3.0],
        "clock_protocol": item["clocks"]["protocol"],
        "verify_external_pins": bool(item["basis_pins"].get("frozen_basis_authenticated")),
        "evolution_backend_requested": "auto",
        "initializer": "dense_initial_state_before_propagate_baseline",
        "initial_state_called": False,
        "saved_t03_radius_copied": False,
        "stability_certificate": False,
        "stage": 0,
        "snapshot_kind": "handoff",
        "work_ledger": episode.empty_work_ledger("coupled"),
        "dense_propagator_stored": False,
        "imbalance": initial["imbalance"],
        "width_scale": initial["width_scale"],
    }


def materialize_directory(output, *, duration, step_cap=STEP_CAP, production=False):
    """Solve each admitted member and evolve it to the handoff with ``propagate_baseline``.

    Duration above ``0.01`` requires ``production=True`` and is the root executor's call.
    """
    directory = episode.assert_campaign_output(output)
    family = _read_family(directory)
    if family.get("materialized"):
        raise FileExistsError("source family is already materialized")
    if (directory / "manifest.json").exists():
        raise FileExistsError("continuation manifest already exists")
    if family.get("stage") == "confirmation_specification":
        selected = family.get("selected", {}).get("cases") or []
        if not selected:
            raise ValueError("no distinct physical regime was selected; refusing to regenerate the family")
        members = [(item["imbalance"], item["width_scale"]) for item in selected]
        nf = CONFIRMATION_NF
    else:
        nf = int(family.get("nf", family.get("exploration_nf", EXPLORATION_NF)))
        by_id = {item["case_id"]: item for item in family["family"]}
        members = [(by_id[identifier]["imbalance"], by_id[identifier]["width_scale"])
                   for identifier in family["cases"]]
    if not members or len(members) > len(IMBALANCES) * len(WIDTH_SCALES):
        raise ValueError("campaign must select from the 15-member family")
    basis = "frozen" if nf in (EXPLORATION_NF, CONFIRMATION_NF) else "owned_geometry_frame"
    decision = prediction.accepted_duration(duration, production=production)
    if decision["production"] and abs(float(duration) - HANDOFF_TIME) > 1e-12:
        raise ValueError("source-family production handoff is T=0.3")
    items = []
    failures = []
    started = time.process_time()
    for imbalance, scale in members:
        identifier = case_id(imbalance, scale, nf)
        if time.process_time() - started >= AGGREGATE_CPU_BUDGET_SECONDS:
            failures.append({"case_id": identifier, "imbalance": imbalance, "width_scale": scale,
                             "status": "budget_stop", "fate": "unresolved", "physical_instability_claimed": False})
            continue
        case_started = time.process_time()
        try:
            prepared = prepare_member(imbalance, scale, nf=nf, basis=basis)
            prepared["initial"]["sample"] = measure_consumers(prepared["pair"], prepared["state"], [0.0]*3, 0.0)
            propagated, clocks, identity = evolve_to_handoff(
                prepared["pair"], prepared["state"], duration=duration,
                step_cap=step_cap, production=production,
            )
            item = {
                "pair": prepared["pair"], "initial": prepared["initial"],
                "basis_pins": prepared["basis_pins"], "propagated": propagated,
                "clocks": clocks, "identity": identity, "step_cap": float(step_cap),
            }
            # Preserve each completed handoff before attempting another source.
            item["committed_handoff"] = _commit_handoff(directory, item, STATIONS)
            items.append(item)
        except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as error:
            failures.append({"case_id": identifier, "imbalance": imbalance, "width_scale": scale,
                             "status": "preparation_failed", "fate": "unresolved",
                             "error": type(error).__name__ + ": " + str(error),
                             "cpu_seconds": float(time.process_time() - case_started),
                             "physical_instability_claimed": False, "non_existence_claimed": False})
        family.update(preparation_failures=failures,
                      initial_data=[item["initial"] for item in items],
                      prepared_cases=[item["initial"]["case_id"] for item in items],
                      preparation_cpu_seconds=float(time.process_time() - started))
        _write_json(_family_path(directory), family)
    return commit_campaign(directory, items, preparation_failures=failures,
                           preparation_cpu_seconds=float(time.process_time() - started))


def _commit_handoff(directory, item, stations):
    record = _handoff_record(item, stations)
    arrays = episode.arrays_from_state(
        item["propagated"]["state"], representation=episode.CANONICAL_PI,
        basis_arrays={"W": np.asarray(item["pair"].geometry_map),
            "source_phi0": np.asarray(item["pair"].source_phi0),
            "source_phi1": np.asarray(item["pair"].source_phi1),
            "observer_columns": np.asarray(item["pair"].reference_columns),
            "source_weights": np.asarray(item["pair"].weights, dtype=float)},
        clocks={"rates": item["clocks"]["rates"], "normal_clocks": item["clocks"]["normal_clocks"]},
    )
    return episode.commit_checkpoint(directory, record, arrays, limit=CHUNK_LIMIT_BYTES)


def commit_campaign(directory, items, *, stations=STATIONS, cpu_budget_seconds=AGGREGATE_CPU_BUDGET_SECONDS,
                    preparation_failures=(), preparation_cpu_seconds=None):
    """Write episode handoff chunks from propagated states. One directory, no second framework."""
    directory = episode.assert_campaign_output(directory)
    manifest_path = directory / "manifest.json"
    if manifest_path.exists():
        raise FileExistsError("continuation manifest already exists; refusing to overwrite")
    directory.mkdir(parents=True, exist_ok=True)
    selected = episode.normalize_stations(stations)
    cases = []
    initial_cpu = 0.0
    initials = []
    for item in items:
        if float(item["identity"]["coordinate_time"]) >= selected[0]:
            raise ValueError("handoff time must lie before the first station")
        committed = item.get("committed_handoff") or _commit_handoff(directory, item, selected)
        if int(committed["npz_bytes"]) > CHUNK_LIMIT_BYTES or int(committed["json_bytes"]) > CHUNK_LIMIT_BYTES:
            raise RuntimeError("checkpoint chunk exceeds 64 MiB")
        cases.append({
            "case_id": committed["case_id"],
            "ordinal": committed["ordinal"],
            "npz": committed["npz"],
            "json": committed["json"],
            "arrays_sha256": committed["arrays_sha256"],
            "nf": committed["nf"],
            "coordinate_time": committed["coordinate_time"],
            "step_cap": committed["step_cap"],
        })
        initial_cpu += float(item["initial"]["initial_cpu_seconds"] or 0.0)
        initials.append(item["initial"])
    handoff_cpu = sum(float(item["identity"].get("handoff_cpu_seconds", 0.0)) for item in items)
    preparation_cpu = (initial_cpu + handoff_cpu if preparation_cpu_seconds is None
                       else float(preparation_cpu_seconds))
    manifest = {
        "schema": episode.SCHEMA,
        "stage": 0,
        "status": "PREPARED" if cases else "PREPARATION_FAILED",
        "campaign": "nsc-discovery-source-family",
        "stations": list(selected),
        "cases": cases,
        "workers": min(WORKERS, len(cases)),
        "cpu_budget_seconds": max(0.0, float(cpu_budget_seconds) - preparation_cpu),
        "aggregate_cpu_budget_seconds": float(cpu_budget_seconds),
        "initial_cpu_seconds": initial_cpu,
        "handoff_cpu_seconds": handoff_cpu,
        "preparation_cpu_seconds": preparation_cpu,
        "preparation_failures": list(preparation_failures),
        "cpu_accounting": "preparation_process_time_plus_continuation_child_process_time",
        "forecast_factor": FORECAST_FACTOR,
        "chunk_limit_bytes": CHUNK_LIMIT_BYTES,
        "threads_per_case": THREADS_PER_CASE,
        "initial_state_called_during_continuation": False,
        "evolved": False,
        "stability_certificate": False,
        "zero_residual_gate": False,
        "max_T_forced": False,
        "child_cpu_seconds": 0.0,
        "executor_pools": 0,
        "saved_t03_radius_copied": False,
    }
    _write_json(manifest_path, manifest)
    family_path = _family_path(directory)
    if family_path.is_file():
        family = json.loads(family_path.read_text())
    else:
        family = dict(specification())
        family["output"] = str(directory)
    family.update(
        stage="materialized",
        materialized=True,
        evolved=False,
        initial_data=initials,
        initial_cpu_seconds=initial_cpu,
        handoff_cpu_seconds=handoff_cpu,
        preparation_cpu_seconds=preparation_cpu,
        preparation_failures=list(preparation_failures),
        cases=[item["case_id"] for item in cases],
        continuation_budget_seconds=manifest["cpu_budget_seconds"],
    )
    _write_json(family_path, family)
    return manifest


def remaining_budget(total, initial_cpu_seconds):
    total = float(total)
    spent = float(initial_cpu_seconds)
    if not math.isfinite(total) or not math.isfinite(spent) or total <= 0.0 or spent < 0.0:
        raise ValueError("CPU budget must be positive and initial CPU cannot be negative")
    if spent > total:
        raise RuntimeError("dense initial solves consumed the aggregate CPU budget; continuation pool not opened")
    return total - spent


def run(output, *, workers=WORKERS, cpu_budget_seconds=AGGREGATE_CPU_BUDGET_SECONDS,
        forecast_factor=FORECAST_FACTOR, memory_limit_bytes=episode.DEFAULT_MEMORY_BYTES,
        max_steps=None, backend="auto", executor=None):
    """One episode pool. The budget that remains after the initial solves is the child allowance."""
    directory = episode.assert_campaign_output(output)
    if not (directory / "manifest.json").is_file():
        raise FileNotFoundError(
            "continuation manifest is missing. Materialize the propagate_baseline handoff before run. "
            "A specification directory is not a campaign."
        )
    manifest = episode.read_manifest(directory)
    if len(manifest.get("cases") or []) < 1:
        raise ValueError("no continuation cases; confirmation was not selected")
    if len(manifest["cases"]) > len(IMBALANCES) * len(WIDTH_SCALES):
        raise ValueError("campaign exceeds the 15-member family")
    if manifest.get("status") == "STATION_REACHED" and max_steps is None:
        return dict(manifest, coordinator_pools=1)
    initial_cpu = float(manifest.get("initial_cpu_seconds") or 0.0)
    preparation_cpu = float(manifest.get("preparation_cpu_seconds", initial_cpu))
    prior_child_cpu = float(manifest.get("cumulative_child_cpu_seconds", manifest.get("child_cpu_seconds", 0.0)))
    prior_coordinator_cpu = float(manifest.get("cumulative_coordinator_cpu_seconds", 0.0))
    allowance = remaining_budget(cpu_budget_seconds, preparation_cpu + prior_child_cpu + prior_coordinator_cpu)
    if allowance <= 0.0:
        raise RuntimeError("aggregate CPU budget exhausted; continuation pool not opened")
    pool_workers = min(int(workers), WORKERS, len(manifest["cases"]))
    if pool_workers < 1:
        raise ValueError("workers must be positive")
    result = episode.run(
        directory,
        workers=pool_workers,
        cpu_budget_seconds=allowance,
        forecast_factor=forecast_factor,
        memory_limit_bytes=memory_limit_bytes,
        max_steps=max_steps,
        backend=backend,
        executor=executor,
    )
    if int(result.get("executor_pools") or 0) != 1:
        raise RuntimeError("source-family continuation opened a pool count other than one")
    assessment_started = time.process_time()
    family_path = _family_path(directory)
    family = None
    if family_path.is_file():
        family = json.loads(family_path.read_text())
        family["child_cpu_seconds"] = result.get("child_cpu_seconds")
        family["regimes"] = assess_regimes(directory, family, result)
        family["question_status"] = answer_question(family["regimes"])
        family["executor_pools"] = 1
        family["evolved"] = True
        family["stage"] = "continued"
    assessment_cpu = float(time.process_time() - assessment_started)
    result = dict(result)
    result["initial_cpu_seconds"] = initial_cpu
    result["preparation_cpu_seconds"] = preparation_cpu
    result["cumulative_child_cpu_seconds"] = prior_child_cpu + float(result.get("child_cpu_seconds") or 0.0)
    result["assessment_cpu_seconds"] = assessment_cpu
    result["cumulative_coordinator_cpu_seconds"] = prior_coordinator_cpu + assessment_cpu
    result["aggregate_cpu_seconds"] = (preparation_cpu + result["cumulative_child_cpu_seconds"]
                                       + result["cumulative_coordinator_cpu_seconds"])
    if family is not None:
        family["aggregate_cpu_seconds"] = result["aggregate_cpu_seconds"]
        family["assessment_cpu_seconds"] = assessment_cpu
        _write_json(family_path, family)
    _write_json(directory / "manifest.json", result)
    result["coordinator_pools"] = 1
    return result


def assess_regimes(directory, family, manifest):
    rows = []
    initial_by_id = {item["case_id"]: item for item in family.get("initial_data", [])}
    result_by_id = {item["case_id"]: item for item in manifest.get("results", [])}
    for case in manifest.get("cases", []):
        identifier = case["case_id"]
        initial = initial_by_id.get(identifier, {})
        row = {"case_id": identifier, "imbalance": initial.get("imbalance"),
               "width_scale": initial.get("width_scale"), "fate": "unresolved",
               "renewal_asserted": False, "non_existence_claimed": False}
        try:
            record, arrays = episode.load_checkpoint(directory, identifier)
            pair = episode.pair_from_arrays(arrays, record)
            state = episode.state_from_arrays(arrays, episode.CANONICAL_PI)
            pair, state, _info = episode.resolve_pair(pair, state, backend="fft")
            sample = measure_consumers(pair, state, arrays["normal_clocks"], record["coordinate_time"])
            earlier = initial.get("sample", {})
            signals = signals_from_samples(earlier, sample)
            if (sample["actual_geometry"] is None or earlier.get("actual_geometry") is None
                or not sample["metric_positive"] or sample["probability_positive"] is not True
                or not sample["actual_geometry"].get("finite", False)):
                signals["unresolved"] = True
            row.update(coordinate_time=record["coordinate_time"], initial_sample=earlier,
                       final_sample=sample, signals=signals,
                       stop=result_by_id.get(identifier, {}).get("stop"),
                       stations_reached=record.get("stations_reached", []))
        except (ValueError, RuntimeError, FileNotFoundError) as error:
            row["error"] = type(error).__name__ + ": " + str(error)
        rows.append(row)
    baseline = next((row for row in rows if row["imbalance"] == BASELINE_IMBALANCE
                     and row["width_scale"] == BASELINE_WIDTH_SCALE), None)
    for row in rows:
        # A shortened or failed case does not share the baseline observation time.
        if baseline and row.get("signals") and baseline.get("signals"):
            if abs(row["coordinate_time"] - baseline["coordinate_time"]) <= 1e-10:
                row["fate"] = fate_class(row["signals"], baseline["signals"], is_baseline=row is baseline)
                row["observable"] = regime_observable(row["signals"])
        elif family.get("stage") == "confirmation_specification":
            row["comparison_requires_coarse_baseline"] = True
    rows.extend(dict(item, fate="unresolved") for item in family.get("preparation_failures", []))
    return rows


def check(output):
    """Validate a specification, or the episode chunks once they exist. Does not evolve."""
    directory = episode.assert_campaign_output(output)
    family = _read_family(directory)
    problems = []
    if family.get("schema") != SCHEMA:
        problems.append("unexpected source-family schema")
    if family.get("occupation_trace") != OCCUPATION_TRACE:
        problems.append("occupation trace is not 3")
    if family.get("universal_one_percent_gate") is not False:
        problems.append("scalar fate installed a universal percent gate")
    if family.get("forced_regeneration") is True:
        problems.append("confirmation was forced")
    if int(family.get("coordinator_pools") or 1) != 1 or family.get("executor_pools") not in (0, 1, None):
        problems.append("pool count is not one")
    manifest_path = directory / "manifest.json"
    episode_report = None
    if manifest_path.is_file():
        episode_report = episode.check(directory)
        manifest = episode.read_manifest(directory)
        if int(manifest.get("chunk_limit_bytes") or CHUNK_LIMIT_BYTES) != CHUNK_LIMIT_BYTES:
            problems.append("chunk limit is not 64 MiB")
        if len(manifest.get("cases") or []) > len(IMBALANCES) * len(WIDTH_SCALES):
            problems.append("more than 15 family continuation cases")
    report = {
        "schema": SCHEMA,
        "ok": not problems and (episode_report is None or episode_report["ok"]),
        "materialized": bool(family.get("materialized")),
        "problems": problems,
        "cases": family.get("cases"),
        "episode": None if episode_report is None else {"ok": episode_report["ok"], "checked": episode_report["checked"]},
        "evolved": False,
    }
    if episode_report is not None and not episode_report["ok"]:
        problems.append("episode checkpoint validation failed")
    if problems:
        raise ValueError("source-family check failed: " + "; ".join(problems))
    return report
