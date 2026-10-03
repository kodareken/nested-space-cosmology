"""Parent episode and checkpoint adapter for the leading carrier.

Public reconstruction entry point: ``reconstruct_parent_pair(arrays, record)``.

The adapter does not own a second Hamiltonian. Rates, the analytic Jacobian
action and RK4 are the leading-Einstein functions. The timestep is the
existing physically scaled restriction. Native FFT carriers run at one thread.
Default calls return a plan and do not write, evolve, or launch a pool.
"""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass, replace
import importlib
import inspect
import json
import multiprocessing
import os
from pathlib import Path
import resource
import time

import numpy as np
from threadpoolctl import threadpool_limits

from . import nsc_discovery_backend as backend
from . import nsc_discovery_episode as episode
from . import nsc_discovery_extent as extent
from . import nsc_discovery_leading_einstein as leading
from . import nsc_discovery_leading_step_control as legacy_step
from . import nsc_discovery_parent as prepared_parent
from . import nsc_discovery_parent_step_control as parent_step
from . import nsc_regional_energy_exchange as regional
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin


SCHEMA = "NSC-DISCOVERY-PARENT-EPISODE-v1"
PARENT_RECORD_SCHEMA = prepared_parent.SCHEMA
# The parent preparation owns these indices: delta (+alpha, 0, -alpha).
POPULATION_NAME = {0: "child_heavy", 1: "balanced", 2: "parent_heavy"}
POPULATION_INDEX = {name: index for index, name in POPULATION_NAME.items()}
POPULATION_IDS = tuple(POPULATION_NAME[index] for index in prepared_parent.POPULATIONS)
PARENT_RECORD_ARRAYS = (
    "Q", "r", "pi_Q", "pi_r", "phi0", "phi1", "W", "weights", "source_columns", "reference_columns",
)
SIGNS = (("plus", 1), ("minus", -1))
STATIONS = (1.0, 3.0, 8.0, 16.0, 24.0)
PRIMARY_NF = 128
PRIMARY_STEP_CAP = 0.001
CONFIRM_NF = 256
CONFIRM_STEP_CAP = 0.0005
MAX_WORKERS = episode.DEFAULT_WORKERS
POOL_START_METHOD = "spawn"
CPU_BUDGET_SECONDS = episode.DEFAULT_CPU_BUDGET_SECONDS
FORECAST_FACTOR = episode.DEFAULT_FORECAST_FACTOR
MEMORY_LIMIT_BYTES = episode.DEFAULT_MEMORY_BYTES
CHUNK_LIMIT_BYTES = episode.CHUNK_LIMIT_BYTES
CHILD_ABS_S = 0.5
COLLAR_ABS_S = 1.0
PARENT_ABS_S = 3.0
ANNULUS_INNER_ABS_S = 1.2
ANNULUS_OUTER_ABS_S = 3.0
LAB = episode.LAB
ROOT = episode.REPO
OUTPUT = LAB / "results/development/nsc-discovery-parent-episode-v1"
OWNERS = (
    "lab/src/recursive_horizons/nsc_discovery_parent_episode.py",
    "lab/src/recursive_horizons/nsc_discovery_parent.py",
    "lab/src/recursive_horizons/nsc_discovery_parent_step_control.py",
    "lab/scripts/derive_nsc_discovery_parent_episode.py",
    "lab/tests/test_nsc_discovery_parent_episode.py",
    "lab/docs/nsc-discovery-parent-episode.md",
    "lab/src/recursive_horizons/nsc_discovery_leading_einstein.py",
    "lab/src/recursive_horizons/nsc_discovery_leading_step_control.py",
    "lab/src/recursive_horizons/nsc_discovery_backend.py",
    "lab/src/recursive_horizons/nsc_discovery_episode.py",
    "lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
    "lab/src/recursive_horizons/nsc_spherical_coupling.py",
    "lab/src/recursive_horizons/nsc_regional_energy_exchange.py",
    "lab/src/recursive_horizons/nsc_discovery_extent.py",
    "lab/src/recursive_horizons/nsc_discovery_tidal.py",
)
_ORIGINAL_NORMALIZE = episode.normalize_control_mode
_ORIGINAL_MAKE_FFT_GRID = backend.make_fft_grid
_THREAD_NAMES = (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
)
_worker_scope = None


@dataclass(frozen=True)
class ParentPair:
    """Full-band holder. Column rank is the stored weight length."""

    grid: object
    geometry_map: np.ndarray
    weights: np.ndarray
    reference_columns: np.ndarray
    source_columns: np.ndarray
    source_metadata: dict
    geometry_metadata: dict
    child_interval: tuple
    parent_interval: tuple
    clock_locations: tuple
    source_phi0: np.ndarray
    source_phi1: np.ndarray
    protected_collar: tuple
    parent_annulus: tuple
    column_rank: int
    common_k: float
    covariance_rank: int


def centre(length):
    return 0.5 * float(length)


def signed_distance(x, length):
    """Signed distance from L/2, wrapped into (-L/2, L/2]."""
    length = float(length)
    return (np.asarray(x, dtype=float) - centre(length) + 0.5 * length) % length - 0.5 * length


def region_intervals(length):
    length = float(length)
    origin = centre(length)

    def span(radius):
        left, right = origin - float(radius), origin + float(radius)
        if left < 0.0 or right > length:
            raise ValueError("signed-distance window wraps; the carrier must cover |s|<=3")
        return (float(left), float(right))

    return {
        "centre": float(origin),
        "child_interval": span(CHILD_ABS_S),
        "protected_collar": span(COLLAR_ABS_S),
        "parent_interval": span(PARENT_ABS_S),
        "parent_annulus": (
            (float(origin - ANNULUS_OUTER_ABS_S), float(origin - ANNULUS_INNER_ABS_S)),
            (float(origin + ANNULUS_INNER_ABS_S), float(origin + ANNULUS_OUTER_ABS_S)),
        ),
        "clock_locations": (float(origin), float(origin + COLLAR_ABS_S), float(origin + PARENT_ABS_S)),
    }


def normalize_control_mode(control_mode):
    if control_mode in ("source_free", "source-free"):
        return "source_free"
    return _ORIGINAL_NORMALIZE(control_mode)


def normalize_stations(stations=None):
    values = STATIONS if stations is None else tuple(float(value) for value in stations)
    allowed = set(STATIONS)
    if not values or any(value not in allowed for value in values):
        raise ValueError("stations must be chosen from 1, 3, 8, 16, 24")
    if list(values) != sorted(values) or len(set(values)) != len(values):
        raise ValueError("stations must be strictly increasing")
    return values


def readonly(value, dtype=None):
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def install_occupations(grid, weights):
    """Replace the quadrature occupations with the stored frame weights."""
    weights = np.array(weights, dtype=float, copy=True)
    if weights.ndim != 1 or weights.size < 1 or not np.isfinite(weights).all():
        raise ValueError("frame weights must be a finite rank vector")
    if np.min(weights) < 0.0 or np.max(weights) > 1.0:
        raise ValueError("frame weights must lie in [0,1]")
    return replace(grid, fine=replace(grid.fine, occupations=weights.copy()))


def full_frame(ng, matrix=None):
    if matrix is None:
        matrix = np.eye(int(ng))
    matrix = np.array(matrix, dtype=float, copy=True)
    if matrix.shape != (int(ng), int(ng)) or not np.isfinite(matrix).all():
        raise ValueError("geometry frame must cover the full ng-by-ng band")
    return matrix


def numerical_binding():
    return {
        "mode": parent_step.NUMERICAL_MODE,
        "rates": "nsc_discovery_leading_einstein.rates",
        "jvp": "nsc_discovery_leading_einstein.jvp",
        "integrator": "nsc_discovery_leading_einstein.rk4_step",
        "step_restriction": "nsc_discovery_parent_step_control.step_restriction",
        "raw_coordinate_norm_used": False,
        "fft": "native_one_thread",
        "workers_max": MAX_WORKERS,
        "forecast_factor": FORECAST_FACTOR,
        "memory_limit_bytes": MEMORY_LIMIT_BYTES,
        "chunk_limit_bytes": CHUNK_LIMIT_BYTES,
        "cpu_budget_seconds": CPU_BUDGET_SECONDS,
        "stations": list(STATIONS),
        "source_recomputed_on_every_rk4_stage": True,
        "projector_source_reset": False,
    }


def source_hashes():
    return {name: episode.file_sha256(ROOT / name) for name in OWNERS}


def _summary(values):
    values = np.asarray(values)
    return {"min": float(np.min(values)), "max": float(np.max(values)),
            "max_abs": float(np.max(np.abs(values)))}


def apply_momentum_sign(pi_q, pi_r, sign, *, already_signed=False):
    """Geometric sign of both momenta. Field columns are not an argument.

    A producer state that was already solved at this sign is kept. Otherwise
    the minus sign negates both canonical momenta of one j0-validated field.
    """
    sign = int(sign)
    if sign not in (1, -1):
        raise ValueError("geometric momentum sign must be +1 or -1")
    pi_q = np.array(pi_q, dtype=float, copy=True)
    pi_r = np.array(pi_r, dtype=float, copy=True)
    if not already_signed and sign == -1:
        pi_q *= -1.0
        pi_r *= -1.0
    return pi_q, pi_r


def _lookup_array(population, arrays, key):
    value = population[key]
    if isinstance(value, np.ndarray):
        return np.array(value, copy=True)
    if isinstance(value, str):
        if arrays is None or key and value not in arrays and value not in population:
            if arrays is None or value not in arrays:
                raise ValueError("preparation array is missing: " + value)
        return np.array(arrays[value], copy=True)
    raise ValueError("preparation field " + key + " must be an array or a payload name")


def _population_block(population, arrays, nf):
    if "by_nf" in population:
        block = population["by_nf"][str(int(nf))]
        return block, arrays
    if int(population["nf"]) != int(nf):
        raise ValueError("refusing to resample population arrays onto a different nf")
    return population, arrays


def discover_prepare_parent():
    """Return the native prepare_parent(nf, population, sign, k_override, profile, cpu_limit)."""
    function = getattr(prepared_parent, "prepare_parent", None)
    return function if callable(function) else None


def _population_weights(rank, name):
    rank = int(rank)
    if name not in POPULATION_INDEX:
        raise ValueError("population must be child_heavy, balanced, or parent_heavy")
    if name == "balanced":
        values = np.full(rank, 0.4)
    elif name == "child_heavy":
        values = np.linspace(0.85, 0.25, rank)
    else:
        values = np.linspace(0.25, 0.85, rank)
    return np.asarray(np.clip(values, 0.05, 0.95), dtype=float)


def load_parent_record(path):
    """Read one frozen NSC-DISCOVERY-PARENT-v1 directory. Does not evolve it."""
    record, arrays = prepared_parent._load(path)
    missing = [name for name in PARENT_RECORD_ARRAYS if name not in arrays]
    if missing:
        raise ValueError("parent record is missing arrays: " + ", ".join(missing))
    if int(record["population"]) not in POPULATION_NAME or int(record["sign"]) not in (1, -1):
        raise ValueError("parent record population or geometric sign is not in the producer set")
    return record, arrays


def iter_parent_records(source):
    source = Path(source)
    if (source / "parent.json").is_file():
        yield source
        return
    found = sorted(path for path in source.iterdir() if path.is_dir() and (path / "parent.json").is_file())
    if not found:
        raise ValueError("no NSC-DISCOVERY-PARENT-v1 record at " + str(source))
    yield from found


def analytic_fixture(*, nf=16, rank=2, common_k=1.0):
    """One shared source geometry and three weight contents. Not a production solve."""
    nf, rank = int(nf), int(rank)
    if nf < 10 or nf % 2 or rank < 1:
        raise ValueError("analytic fixture needs an even nf of at least 10 and a positive rank")
    grid = galerkin.build_grid(nf, gauge="conformal")
    phases = []
    for index in range(rank):
        mode = (index - (rank - 1) / 2.0) + 0.5
        phases.append(np.exp(2j * np.pi * mode * grid.xi_f / grid.length))
    raw = np.column_stack(phases)
    orthogonal, _r = np.linalg.qr(np.vstack((raw, 1j * raw)))
    phi0 = orthogonal[:nf].astype(complex, copy=True)
    phi1 = orthogonal[nf:].astype(complex, copy=True)
    q0 = coupling.CALIBRATION["b0"] / coupling.CALIBRATION["a0"]
    pi_q = 0.01 * np.cos(2 * np.pi * grid.xi_g / grid.length)
    pi_r = 0.02 * np.sin(2 * np.pi * grid.xi_g / grid.length)
    populations = []
    for name in POPULATION_IDS:
        populations.append({
            "population": POPULATION_INDEX[name],
            "population_id": name,
            "nf": nf,
            "column_rank": rank,
            "momentum_representation": episode.CANONICAL_PI,
            "momenta_already_signed": False,
            "phi0": phi0,
            "phi1": phi1,
            "weights": _population_weights(rank, name),
            "Q": np.full(grid.ng, q0),
            "r": np.ones(grid.ng),
            "pi_Q": pi_q.copy(),
            "pi_r": pi_r.copy(),
            "W": np.eye(grid.ng),
            "geometry_frame": "identity",
            "observer_columns": np.vstack((phi0, phi1)),
            "analytic_fixture": True,
        })
    return {
        "schema": PARENT_RECORD_SCHEMA,
        "common_k": float(common_k),
        "length": float(grid.length),
        "nf": nf,
        "column_rank": rank,
        "populations": populations,
        "analytic_fixture": True,
        "production": False,
        "same_source_geometry": True,
        "held_out_nonlinear_measurement_before_prediction": False,
    }


def _frame_from_population(block, arrays, ng):
    identity = block.get("geometry_frame") == "identity" or "W" not in block
    if identity:
        if "W" in block:
            matrix = full_frame(ng, _lookup_array(block, arrays, "W"))
            if not np.array_equal(matrix, np.eye(ng)):
                raise ValueError("identity frame does not match the supplied W")
            return matrix, "identity"
        return full_frame(ng), "identity"
    return full_frame(ng, _lookup_array(block, arrays, "W")), "supplied"


def build_parent_pair(grid, *, geometry_map, weights, source_phi0, source_phi1,
                      observer_columns, source_metadata, geometry_metadata, common_k,
                      child_interval=None, parent_interval=None, clock_locations=None):
    intervals = region_intervals(grid.length)
    rank = int(np.asarray(weights).shape[0])
    if source_phi0.shape != (grid.nf, rank) or source_phi1.shape != (grid.nf, rank):
        raise ValueError("source columns must be nf by the weight rank")
    if observer_columns.shape != (2 * grid.nf, rank):
        raise ValueError("observer columns must be 2*nf by the weight rank")
    if geometry_map.shape != (grid.ng, grid.ng):
        raise ValueError("geometry frame dropped part of the band")
    occupied = install_occupations(grid, weights)
    source = np.vstack((source_phi0, source_phi1))
    covariance_rank = source_covariance_rank(source, weights)
    return ParentPair(
        occupied, readonly(geometry_map), readonly(weights), readonly(observer_columns),
        readonly(source), dict(source_metadata), dict(geometry_metadata),
        tuple(child_interval or intervals["child_interval"]),
        tuple(parent_interval or intervals["parent_interval"]),
        tuple(clock_locations or intervals["clock_locations"]),
        readonly(source_phi0), readonly(source_phi1), intervals["protected_collar"],
        intervals["parent_annulus"], rank, float(common_k), covariance_rank,
    )


def source_covariance_rank(columns, weights):
    """Actual rank of Phi diag(c) Phi†, independent of stored slot count."""
    weighted = np.asarray(columns)*np.sqrt(np.asarray(weights))[None, :]
    gram = weighted.conj().T@weighted
    eigenvalues = np.linalg.eigvalsh(gram)
    maximum = float(np.max(abs(eigenvalues)))
    threshold = np.finfo(float).eps*max(weighted.shape)*maximum
    return int(np.count_nonzero(eigenvalues > threshold)) if maximum else 0


def pair_and_state_from_population(population, *, common_k, sign, nf, arrays=None, length=None):
    block, arrays = _population_block(population, arrays, nf)
    grid = galerkin.build_grid(int(nf), gauge="conformal")
    if length is not None and abs(float(length) - grid.length) > 1e-12:
        raise ValueError("preparation length does not match the carrier")
    weights = np.array(_lookup_array(block, arrays, "weights"), dtype=float, copy=True)
    rank = int(block.get("column_rank", weights.size))
    if weights.shape != (rank,):
        raise ValueError("weight rank and declared column rank differ")
    phi0 = np.array(_lookup_array(block, arrays, "phi0"), dtype=complex, copy=True)
    phi1 = np.array(_lookup_array(block, arrays, "phi1"), dtype=complex, copy=True)
    if phi0.shape != (grid.nf, rank) or phi1.shape != (grid.nf, rank):
        raise ValueError("population columns must cover the full AP band at this rank")
    matrix, frame_name = _frame_from_population(block, arrays, grid.ng)
    observer = _lookup_array(block, arrays, "observer_columns") if "observer_columns" in block else np.vstack((phi0, phi1))
    representation = block.get("momentum_representation", episode.CANONICAL_PI)
    if representation == episode.CANONICAL_PI:
        pi_q = np.array(_lookup_array(block, arrays, "pi_Q"), dtype=float, copy=True)
        pi_r = np.array(_lookup_array(block, arrays, "pi_R") if "pi_R" in block else _lookup_array(block, arrays, "pi_r"),
                        dtype=float, copy=True)
    elif representation == "nodal_momentum":
        pi_q = grid.dx_g * (matrix.T @ np.array(_lookup_array(block, arrays, "p_Q"), dtype=float))
        pi_r = grid.dx_g * (matrix.T @ np.array(_lookup_array(block, arrays, "p_r"), dtype=float))
    else:
        raise ValueError("momentum representation must be canonical_pi or nodal_momentum")
    pi_q, pi_r = apply_momentum_sign(pi_q, pi_r, sign, already_signed=bool(block.get("momenta_already_signed", False)))
    geometry = np.array(_lookup_array(block, arrays, "Q"), dtype=float, copy=True)
    radius = np.array(_lookup_array(block, arrays, "r"), dtype=float, copy=True)
    population_id = population.get("population_id") or POPULATION_NAME[int(population["population"])]
    metadata = {
        "population": int(population.get("population", POPULATION_INDEX[population_id])),
        "population_id": population_id,
        "momentum_sign": int(sign),
        "momentum_sign_is_geometric": True,
        "field_conjugated": False,
        "common_k": float(common_k),
        "column_rank": rank,
        "covariance_rank": rank,
        "control_mode": population.get("control_mode", "coupled"),
        "analytic_fixture": bool(population.get("analytic_fixture", False)),
        "production": False if population.get("analytic_fixture") else True,
    }
    geometry_metadata = {
        "frame": frame_name,
        "all_geometry_degrees_retained": True,
        "ambient_field_active": True,
        "geometry_degrees": int(grid.ng),
        "fermion_degrees": int(grid.nf),
        "identity_frame": frame_name == "identity",
        "centre": centre(grid.length),
        "signed_distance": "s=(x-L/2) wrapped into (-L/2, L/2]",
    }
    pair = build_parent_pair(
        grid, geometry_map=matrix, weights=weights, source_phi0=phi0, source_phi1=phi1,
        observer_columns=np.array(observer, dtype=complex, copy=True),
        source_metadata=metadata, geometry_metadata=geometry_metadata, common_k=common_k,
    )
    state = leading.State(geometry, radius, pi_q, pi_r, phi0, phi1)
    return pair, state


def _column_pair(arrays, nf, rank, stacked, first, second):
    if stacked in arrays:
        values = np.array(arrays[stacked], copy=True)
        if values.shape != (2 * nf, rank):
            raise ValueError(stacked + " is not 2*nf by the stored rank")
        return values[:nf], values[nf:]
    return np.array(arrays[first], copy=True), np.array(arrays[second], copy=True)


def _record_windows(record, grid):
    intervals = region_intervals(grid.length)
    nested = record.get("intervals") or {}
    child = tuple(float(value) for value in (record.get("child_interval") or nested.get("child") or intervals["child_interval"]))
    parent = tuple(float(value) for value in (record.get("parent_interval") or nested.get("parent") or intervals["parent_interval"]))
    clocks = tuple(float(value) for value in (record.get("clock_locations") or intervals["clock_locations"]))
    if child != intervals["child_interval"] or parent != intervals["parent_interval"]:
        raise ValueError("cached child or parent window is not the signed-distance window")
    return child, parent, clocks


def reconstruct_parent_pair(arrays, record):
    """Rebuild a ParentPair from cached arrays. No SVD and no source reset."""
    if record.get("projector_source_reset") or record.get("source_reset"):
        raise ValueError("reconstruction refuses a projector source reset")
    nf = int(record["nf"])
    weights = np.array(arrays["weights"] if "weights" in arrays else arrays["source_weights"], dtype=float, copy=True)
    rank = int(record.get("column_rank", weights.size))
    if weights.shape != (rank,):
        raise ValueError("cached column rank does not match the weight vector")
    source_phi0, source_phi1 = _column_pair(arrays, nf, rank, "source_columns", "source_phi0", "source_phi1")
    if "reference_columns" in arrays:
        observer = np.array(arrays["reference_columns"], copy=True)
    elif "observer_columns" in arrays:
        observer = np.array(arrays["observer_columns"], copy=True)
    else:
        raise ValueError("cached observer columns are missing")
    if source_phi0.shape != (nf, rank) or source_phi1.shape != (nf, rank) or observer.shape != (2 * nf, rank):
        raise ValueError("cached source or observer is not the stored AP rank")
    grid = galerkin.build_grid(nf, gauge="conformal")
    matrix = full_frame(grid.ng, arrays["W"])
    if record.get("geometry_frame", record.get("geometry_map")) == "identity" and not np.array_equal(matrix, np.eye(grid.ng)):
        raise ValueError("identity frame record does not match cached W")
    child, parent, clocks = _record_windows(record, grid)
    if "normal_clocks" in arrays and arrays["normal_clocks"].shape != (len(clocks),):
        raise ValueError("cached normal clocks do not match the clock stations")
    common_k = record.get("common_k", record.get("k_common", record.get("k")))
    pair = build_parent_pair(
        grid, geometry_map=matrix, weights=weights, source_phi0=source_phi0, source_phi1=source_phi1,
        observer_columns=observer, source_metadata=dict(record.get("source_metadata") or record.get("weight_metadata") or {}),
        geometry_metadata=dict(record.get("geometry_metadata") or {}), common_k=float(common_k),
        child_interval=child, parent_interval=parent, clock_locations=clocks,
    )
    if "covariance_rank" in record and int(record["covariance_rank"]) != pair.covariance_rank:
        raise ValueError("cached covariance rank does not match the actual source")
    return pair


def adopt_parent_record(record, arrays):
    """Adopt one frozen producer record. Phi and the solved momenta are copied."""
    pair = reconstruct_parent_pair(arrays, record)
    state = leading.State(
        np.array(arrays["Q"], dtype=float, copy=True), np.array(arrays["r"], dtype=float, copy=True),
        np.array(arrays["pi_Q"], dtype=float, copy=True), np.array(arrays["pi_r"], dtype=float, copy=True),
        np.array(arrays["phi0"], dtype=complex, copy=True), np.array(arrays["phi1"], dtype=complex, copy=True),
    )
    return pair, state


def parent_pair_from_arrays(arrays, record):
    return reconstruct_parent_pair(arrays, record)


def state_from_arrays(arrays, representation=episode.CANONICAL_PI):
    return leading.state_from_arrays(arrays, representation)


def arrays_from_parent(pair, state, *, clocks, clock_rates):
    return {
        "Q": np.array(state.Q, copy=True),
        "r": np.array(state.r, copy=True),
        "pi_Q": np.array(state.p_Q, copy=True),
        "pi_r": np.array(state.p_r, copy=True),
        "phi0": np.array(state.phi0, copy=True),
        "phi1": np.array(state.phi1, copy=True),
        "W": np.array(pair.geometry_map, copy=True),
        "source_phi0": np.array(pair.source_phi0, copy=True),
        "source_phi1": np.array(pair.source_phi1, copy=True),
        "observer_columns": np.array(pair.reference_columns, copy=True),
        "source_weights": np.array(pair.weights, copy=True),
        "clock_rates": np.array(clock_rates, dtype=float, copy=True),
        "normal_clocks": np.array(clocks, dtype=float, copy=True),
    }


def initial_ledger(mode):
    ledger = {name: 0.0 for name in episode.WORK_CHANNELS}
    ledger.update(quadrature="trapezoid_coordinate_time", control_mode=normalize_control_mode(mode),
                  dense_propagator_stored=False, stocks_are_coordinate_trapezoid=True)
    return ledger


def parent_pins(arrays, state):
    fields = [episode.array_sha256(np.ascontiguousarray(getattr(state, name))) for name in leading.FIELDS]
    return {
        "source_columns_sha256": episode.array_sha256(np.vstack((arrays["source_phi0"], arrays["source_phi1"]))),
        "observer_columns_sha256": episode.array_sha256(arrays["observer_columns"]),
        "W_sha256": episode.array_sha256(arrays["W"]),
        "weights_sha256": episode.array_sha256(arrays["source_weights"]),
        "phi_sha256": episode.array_sha256(np.vstack((state.phi0, state.phi1))),
        "state_field_sha256": fields,
        "object_alias": any(getattr(state, name) is arrays.get(stored)
                            for name, stored in (("phi0", "source_phi0"), ("phi1", "source_phi1"))),
    }


def pins_unchanged(before, after):
    if not before or not after:
        return
    for key in ("W_sha256", "source_columns_sha256", "observer_columns_sha256", "weights_sha256"):
        if before.get(key) != after.get(key):
            raise ValueError("restart changed a frozen source pin: " + key)
    if after.get("object_alias"):
        raise ValueError("state Phi is aliased to the source columns")


def real_array_names(representation):
    return leading.real_array_names(representation)


def fft_grid_for_parent(grid):
    """Build the native FFT carrier before installing the actual occupation width."""
    weights = np.array(grid.fine.occupations, dtype=float, copy=True)
    if backend.operator_backend(grid) == "fft":
        return replace(grid, fine=replace(grid.fine, occupations=weights))
    carrier = galerkin.build_grid(grid.nf, quadrature=grid.nq, gauge=grid.gauge)
    transformed = _ORIGINAL_MAKE_FFT_GRID(carrier)
    return replace(transformed, fine=replace(transformed.fine, occupations=weights))


_SAVED_MAKE_FFT_PAIR = backend.make_fft_pair


def make_fft_pair(pair):
    if isinstance(pair, ParentPair):
        if backend.operator_backend(pair.grid) == "fft":
            return pair
        return replace(pair, grid=fft_grid_for_parent(pair.grid))
    return _SAVED_MAKE_FFT_PAIR(pair)


def scaled_step_restriction(pair, state, step_cap, control_mode=None):
    """Delegate admission to the shared rank-general owner and pin the source."""
    if control_mode is None:
        control_mode = pair.source_metadata.get("control_mode", "coupled")
    dt, record = parent_step.step_restriction(pair, state, step_cap, control_mode=control_mode)
    record = dict(record)
    record["source_pin"] = {
        "weights_sha256": episode.array_sha256(np.ascontiguousarray(pair.weights)),
        "source_columns_sha256": episode.array_sha256(np.ascontiguousarray(pair.source_columns)),
        "column_rank": int(pair.column_rank),
    }
    record["raw_coordinate_norm_used"] = False
    return dt, record


def source_free_reprepare(pair, state, *, sign=None, deadline=None):
    """Empty the source on a copy and solve momenta for the new C and D.

    Weights and the caller's arrays stay. Zero weights are not the constraint.
    """
    before = {name: np.array(getattr(state, name), copy=True) for name in leading.FIELDS}
    identities = {name: id(getattr(state, name)) for name in leading.FIELDS}
    nodal = leading.decode(pair, state)
    rank = int(np.asarray(pair.weights).shape[0])
    empty0 = np.zeros((pair.grid.nf, rank), dtype=complex)
    empty1 = np.zeros((pair.grid.nf, rank), dtype=complex)
    metadata = dict(pair.source_metadata)
    metadata.update(control_mode="source_free", intended_empty_source=True, field_deleted_inplace=False,
                    weights_preserved=True, weights_used_as_constraint_solution=False)
    fresh = galerkin.build_grid(pair.grid.nf, quadrature=pair.grid.nq, gauge="conformal")
    new_pair = build_parent_pair(
        fresh, geometry_map=np.array(pair.geometry_map, copy=True),
        weights=np.array(pair.weights, dtype=float, copy=True),
        source_phi0=empty0, source_phi1=empty1,
        observer_columns=np.array(pair.reference_columns, copy=True),
        source_metadata=metadata, geometry_metadata=dict(pair.geometry_metadata),
        common_k=pair.common_k, clock_locations=pair.clock_locations,
        child_interval=pair.child_interval, parent_interval=pair.parent_interval,
    )
    sign = int(sign if sign is not None else metadata.get("momentum_sign", 1))
    try:
        seed = prepared_parent._offset_seed(new_pair, nodal.r, nodal.Q, empty0, empty1, 0.05)
        seed["suggested_k"] = float(seed["kmin"] + max(0.05 * float(seed["margin_scale"]), 1e-8))
        problem = prepared_parent._MomentumProblem(
            new_pair, nodal.r, nodal.Q, empty0, empty1, sign, 0.0, seed["g"])
        theta = prepared_parent._seed_theta(problem, seed, sign)
        theta, correction = prepared_parent.correct_momenta(problem, theta, deadline or time.process_time() + 15.0)
        solved = prepared_parent._state_from_theta(new_pair, problem, theta)
        measured = leading.constraints(new_pair, solved)
    finally:
        for name in leading.FIELDS:
            if id(getattr(state, name)) != identities[name] or not np.array_equal(getattr(state, name), before[name]):
                raise RuntimeError("source-free reprepare wrote the caller's field")
    solved_constraints = bool(correction["converged"] and measured["raw_C_max"] <= 1e-6 and measured["D_max"] <= 1e-6)
    report = {
        "source_reprepared": True, "intended_empty_source": True, "field_deleted_inplace": False,
        "column_rank": rank, "covariance_rank": int(new_pair.covariance_rank),
        "source_column_Gram_distance_from_identity": 1.0,
        "weights_preserved": True, "weights_deleted": False, "weights_used_as_constraint_solution": False,
        "converged": bool(correction["converged"]), "constraint_solved": solved_constraints,
        "blocker": correction.get("blocker"), "momentum_sign": sign,
        "leading_constraints": {key: measured[key] for key in ("raw_C_max", "D_max", "h_c_max", "source_current_max")},
    }
    return new_pair, solved, report


def discover_observer_hook():
    try:
        module = importlib.import_module(".nsc_discovery_parent_observer", __package__)
    except ImportError:
        return None
    for name in ("observe_parent", "parent_observation", "observe"):
        hook = getattr(module, name, None)
        if callable(hook):
            return hook
    return None


def call_observer_hook(hook, pair, state, instant, row, control_mode):
    """Call OWNER3 with its declared arguments. The row is not passed as a mode."""
    parameters = inspect.signature(hook).parameters
    names = set(parameters)
    arguments = {}
    mode = normalize_control_mode(control_mode)
    if "control_mode" in names:
        arguments["control_mode"] = mode
    if "observation" in names:
        arguments["observation"] = row
    elif "row" in names:
        arguments["row"] = row
    if any(item.kind == inspect.Parameter.VAR_KEYWORD for item in parameters.values()):
        arguments.setdefault("control_mode", mode)
    return hook(pair, state, instant, **arguments)


def _rank_general_observation(pair, state, instant, mode):
    rate, bundle = leading.rates(pair, state, return_bundle=True, control_mode=mode)
    fine, system = bundle["fine_state"], bundle["fine_system"]
    curvature, velocity, _acceleration = leading.metric_jets(pair, state, rate, bundle, control_mode=mode)
    ledger = regional.matter_ledger(system, fine)
    balance = regional.proper_balance_terms(system, fine, velocity, ledger)
    grid = pair.grid
    integral = lambda values, interval: extent.real_interval_integral(grid, values, interval)
    sample = lambda values, locations: extent.real_periodic_values(grid, values, locations)
    mass = np.sum((np.abs(fine.phi0) ** 2 + np.abs(fine.phi1) ** 2) * pair.weights, axis=1)
    proper = fine.r * fine.Q
    child = pair.child_interval
    boundary = sample(balance["flux_nodal"] / grid.dx_q, child)
    parent_boundary = sample(balance["flux_nodal"] / grid.dx_q, pair.parent_interval)
    child_proper = integral(proper, child)
    return {
        "time": float(instant), "control_mode": mode, "model": "leading_Einstein_EFT_diagnostic",
        "column_rank": int(pair.column_rank),
        "child_proper_length": child_proper,
        "child_probability": integral(mass / grid.dx_q, child),
        "parent_probability": integral(mass / grid.dx_q, pair.parent_interval),
        "coordinate_fieldwork": float(rate.fieldwork_power),
        "pressure_work": float(np.sum(balance["proper_pressure_work"])),
        "lapse_work": float(np.sum(balance["momentum_lapse_work"])),
        "boundary_child": float(boundary[0] - boundary[1]),
        "boundary_parent": float(parent_boundary[0] - parent_boundary[1]),
        "normal_clock_rates": leading.clock_rates(pair, state).tolist(),
        "energy": leading.energy(pair, state),
        "constraints": leading.constraints(pair, state),
        "tides": {name: _summary(curvature["tides"][name]) for name in ("R_0101", "R_0202", "R4", "owned_W")},
        "curvature": {name: _summary(curvature["base"][name]) for name in ("Ricci2", "K")},
        "leading_observe_included": False,
    }


def _annulus_probability(pair, fine):
    mass = np.sum((np.abs(fine.phi0) ** 2 + np.abs(fine.phi1) ** 2) * pair.weights, axis=1) / pair.grid.dx_q
    return float(sum(extent.real_interval_integral(pair.grid, mass, interval) for interval in pair.parent_annulus))


def scalar_observation(pair, state, instant, control_mode="coupled"):
    mode = normalize_control_mode(control_mode)
    rank = int(pair.column_rank)
    row = _rank_general_observation(pair, state, instant, mode)
    fine = leading.fine_state(pair, state)
    system = leading.active_system(pair.grid, fine)
    source = coupling.source_from_columns(system, fine)
    current = source["force_beta"] / pair.grid.dx_q
    density = source["force_L"] / pair.grid.dx_q
    row.update(
        column_rank=rank, covariance_rank=source_covariance_rank(np.vstack((state.phi0, state.phi1)), pair.weights),
        common_k=float(pair.common_k),
        source_column_Gram_distance_from_identity=float(np.max(abs(
            state.phi0.conj().T@state.phi0+state.phi1.conj().T@state.phi1-np.eye(rank)))),
        source_current=_summary(current), source_density=_summary(density),
        source_current_mean=float(np.mean(current)),
        normal_clock_rates=leading.clock_rates(pair, state).tolist(),
        parent_annulus_probability=_annulus_probability(pair, fine),
        protected_collar=list(pair.protected_collar),
        nonlinear_measurement_before_prediction=False,
        population_gradient_held_out=None,
        dense_propagator_stored=False,
        stability_certificate=False,
    )
    hook = discover_observer_hook()
    if hook is None:
        row["external_observer"] = {"available": False}
    else:
        try:
            row["external_observer"] = {
                "available": True,
                "result": call_observer_hook(hook, pair, state, instant, row, mode),
            }
        except Exception as error:
            row["external_observer"] = {
                "available": True,
                "error": type(error).__name__ + ": " + str(error),
            }
    return row


def scalar_observation_row(pair, state, time_value, control_mode="coupled"):
    started = time.process_time()
    full = scalar_observation(pair, state, time_value, control_mode)
    row = {name: full.get(name) for name in episode.WORK_CHANNELS}
    row.update(time=float(time_value), control_mode=normalize_control_mode(control_mode),
               observation_cpu_seconds=time.process_time() - started,
               dense_propagator_stored=False, quadrature="trapezoid_coordinate_time",
               column_rank=full["column_rank"], source_current=full["source_current"],
               source_density=full["source_density"], normal_clock_rates=full["normal_clock_rates"],
               curvature=full.get("curvature"), tides=full.get("tides"),
               constraints=full.get("constraints"), energy=full.get("energy"),
               parent_annulus_probability=full["parent_annulus_probability"],
               nonlinear_measurement_before_prediction=False,
               population_gradient_held_out=None,
               external_observer=full["external_observer"],
               leading_observe_included=full["leading_observe_included"])
    return row


def episode_observe(pair, state, instant, *, control_mode="coupled"):
    return scalar_observation(pair, state, instant, control_mode)


def case_specs(preparation=None, *, confirm=False):
    if preparation is None:
        common_k = None
        populations = [{"population_id": name} for name in POPULATION_IDS]
        nf = CONFIRM_NF if confirm else PRIMARY_NF
        cap = CONFIRM_STEP_CAP if confirm else PRIMARY_STEP_CAP
        analytic = False
    else:
        common_k = preparation.get("common_k", preparation.get("k_common", preparation.get("k")))
        populations = list(preparation["populations"])
        if preparation.get("analytic_fixture"):
            nf = int(preparation["nf"])
            cap = CONFIRM_STEP_CAP if confirm else PRIMARY_STEP_CAP
            analytic = True
        else:
            nf = CONFIRM_NF if confirm else PRIMARY_NF
            cap = CONFIRM_STEP_CAP if confirm else PRIMARY_STEP_CAP
            analytic = False
    specs = []
    for sign_name, sign in SIGNS:
        for population in populations:
            specs.append({
                "case_id": f"nf{nf}_{sign_name}_{population['population_id']}_dt{cap}",
                "nf": int(nf), "population": int(population.get("population", POPULATION_INDEX[population["population_id"]])),
                "population_id": population["population_id"],
                "sign_name": sign_name, "momentum_sign": int(sign), "step_cap": float(cap),
                "control_mode": "coupled", "common_k": common_k, "confirm": bool(confirm),
                "analytic_fixture": analytic, "stations": list(STATIONS),
            })
    return specs


def plan(preparation=None):
    """Pure case list. No files, evolution, or process pool."""
    return {
        "schema": SCHEMA, "status": "PLAN", "pure": True, "evolved": False, "bytes_written": 0,
        "pool_launched": False, "cases": case_specs(preparation, confirm=False),
        "confirmation_cases": case_specs(preparation, confirm=True),
        "stations": list(STATIONS), "numerical_binding": numerical_binding(),
        "owner1_prepare_parent_importable": discover_prepare_parent() is not None,
        "prepare_parent_api": "prepare_parent(nf, population=0|1|2, sign=+1|-1, k_override=None, profile=None, cpu_limit=30) -> (pair, State, report)",
        "accepted_parent_record": {
            "schema": PARENT_RECORD_SCHEMA,
            "files": ["parent.json", "parent.npz"],
            "arrays": list(PARENT_RECORD_ARRAYS),
            "one_record_per": ["nf", "population", "sign"],
            "population_names": {str(index): name for index, name in POPULATION_NAME.items()},
            "momentum_sign": "geometric sign of both p_Q and p_r; phi is not conjugated",
        },
        "held_out_nonlinear_measurement_before_prediction": False,
        "production_campaign": False,
    }


def assert_adapter():
    if episode.SCHEMA != SCHEMA or episode.step_restriction is not scaled_step_restriction:
        raise RuntimeError("parent episode runtime adapter is not active")
    if backend.make_fft_pair is not make_fft_pair:
        raise RuntimeError("parent FFT adapter is not active")


def commit_parent_checkpoint(directory, record, arrays):
    assert_adapter()
    if any(arrays[name] is arrays[source] for name in ("phi0", "phi1") for source in ("source_phi0", "source_phi1")):
        raise ValueError("refusing to store a state field aliased to the source")
    return episode.commit_checkpoint(directory, record, arrays, limit=CHUNK_LIMIT_BYTES)


def build_case_record(spec, pair, state, *, producer_commit, source_binding, input_binding):
    clocks = np.zeros(len(pair.clock_locations))
    rates = leading.clock_rates(pair, state)
    arrays = arrays_from_parent(pair, state, clocks=clocks, clock_rates=rates)
    record = dict(spec)
    record.update(
        schema=SCHEMA, coordinate_time=0.0, steps=0, ordinal=0, status="PREPARED",
        stations=list(STATIONS), control_mode=spec.get("control_mode", "coupled"),
        momentum_representation=episode.CANONICAL_PI, snapshot_kind="handoff",
        column_rank=int(pair.column_rank), covariance_rank=int(pair.covariance_rank),
        common_k=float(pair.common_k), geometry_frame=pair.geometry_metadata.get("frame"),
        child_interval=list(pair.child_interval), parent_interval=list(pair.parent_interval),
        protected_collar=list(pair.protected_collar),
        parent_annulus=[list(item) for item in pair.parent_annulus],
        clock_locations=list(pair.clock_locations),
        source_metadata=dict(pair.source_metadata), geometry_metadata=dict(pair.geometry_metadata),
        work_ledger=initial_ledger(spec.get("control_mode", "coupled")),
        channel_sample=None, channel_time=None, normal_clocks_initial=clocks.tolist(),
        numerical_binding=numerical_binding(), producing_commit=producer_commit,
        source_hashes=source_binding, input_hashes=input_binding,
        projector_source_reset=False, source_reset=False, initial_state_called=False,
        evolved=False, verify_external_pins=False, stability_certificate=False,
        held_out_nonlinear_measurement_before_prediction=False,
    )
    record["source_pins"] = parent_pins(arrays, state)
    return record, arrays


def _directory_files(path):
    path = Path(path)
    if not path.exists():
        return {}
    return {item.name: episode.file_sha256(item) for item in path.iterdir() if item.is_file()}


def prepare(source=None, output=OUTPUT, *, execute=False, producer_commit=None,
            cpu_budget=CPU_BUDGET_SECONDS, confirm=False, preparation=None):
    if not 0 < float(cpu_budget) <= CPU_BUDGET_SECONDS:
        raise ValueError("parent episode budget must be in (0, 21600] CPU seconds")
    report = plan(preparation)
    report.update(status="PREFLIGHT", cpu_budget_seconds=float(cpu_budget), confirm=bool(confirm),
                  output=str(output), execute=False)
    if not execute:
        return report
    if not producer_commit:
        raise ValueError("explicit preparation requires a frozen producer commit")
    if preparation is not None and preparation.get("analytic_fixture"):
        raise ValueError("analytic fixtures are not a production preparation")
    binding = source_hashes()
    try:
        commit = leading.preparation._git_hashes(producer_commit, binding)
    except Exception as error:
        raise ValueError("frozen producer does not match the current source closure") from error
    if source is None:
        raise ValueError("explicit preparation requires one frozen parent record per case")
    records = [(path, *load_parent_record(path)) for path in iter_parent_records(source)]
    if confirm:
        records = [item for item in records if int(item[1]["nf"]) == CONFIRM_NF]
    elif any(int(item[1]["nf"]) == PRIMARY_NF for item in records):
        records = [item for item in records if int(item[1]["nf"]) == PRIMARY_NF]
    destination = episode.assert_campaign_output(output)
    if destination.exists():
        raise FileExistsError("parent episode output already exists")
    started = time.process_time()
    destination.mkdir(parents=True, exist_ok=False)
    inputs = {}
    for path, _record, _arrays in records:
        inputs[str(path / "parent.json")] = episode.file_sha256(path / "parent.json")
        inputs[str(path / "parent.npz")] = episode.file_sha256(path / "parent.npz")
    cases = []
    with runtime_adapter():
        for path, parent_record, parent_arrays in records:
            pair, state = adopt_parent_record(parent_record, parent_arrays)
            sign = int(parent_record["sign"])
            population = int(parent_record["population"])
            nf = int(parent_record["nf"])
            cap = CONFIRM_STEP_CAP if nf == CONFIRM_NF else PRIMARY_STEP_CAP
            name = POPULATION_NAME[population]
            spec = {
                "case_id": f"nf{nf}_{'plus' if sign == 1 else 'minus'}_{name}_dt{cap}",
                "nf": nf, "population": population, "population_id": name,
                "sign_name": "plus" if sign == 1 else "minus", "momentum_sign": sign,
                "step_cap": float(cap), "control_mode": "coupled",
                "common_k": float(parent_record.get("k_common", parent_record["k"])),
                "stations": list(STATIONS), "parent_record": str(path),
                "momenta_already_signed": True, "field_conjugated": False,
            }
            record, arrays = build_case_record(spec, pair, state, producer_commit=commit,
                                               source_binding=binding, input_binding=inputs)
            record["parent_record_schema"] = parent_record.get("schema")
            record["parent_sign"] = sign
            record["parent_population"] = population
            record["parent_k"] = parent_record.get("k")
            committed = commit_parent_checkpoint(destination, record, arrays)
            cases.append({name: committed[name] for name in ("case_id", "ordinal", "npz", "json", "coordinate_time", "steps")})
        used = time.process_time() - started
        if used > float(cpu_budget):
            raise RuntimeError("parent preparation exceeded its CPU budget")
        report.update(status="PREPARED", stage=0, cases=cases, producing_commit=commit,
                      source_hashes=binding, input_hashes=inputs, preparation_cpu_seconds=used,
                      child_cpu_seconds=0.0, evolved=False, pool_launched=False,
                      forecast_factor=FORECAST_FACTOR, memory_limit_bytes=MEMORY_LIMIT_BYTES,
                      chunk_limit_bytes=CHUNK_LIMIT_BYTES, workers_max=MAX_WORKERS)
        episode._write_json(episode._manifest_path(destination), report)
        episode._ledger_update(destination, pilot=used, budget=float(cpu_budget))
    return report


def execute_case(directory, case_id, *, cpu_allowance, forecast_factor=FORECAST_FACTOR,
                 memory_limit_bytes=MEMORY_LIMIT_BYTES, max_steps=None, fft=None, backend=None,
                 diagnostics=False, **kwargs):
    assert_adapter()
    if abs(float(forecast_factor) - FORECAST_FACTOR) > 1e-12:
        raise ValueError("forecast factor is 1.5")
    record, arrays = episode.load_checkpoint(directory, case_id)
    if record.get("projector_source_reset") or record.get("source_reset"):
        raise ValueError("resume refuses a projector source reset")
    if record.get("status") in ("station_reached", "chart_exit"):
        return {"case_id": case_id, "status": record["status"], "child_cpu_seconds": 0.0,
                "coordinate_time": record["coordinate_time"], "steps": record["steps"],
                "stations_reached": record.get("stations_reached", []), "stop": record.get("stop"),
                "chunk": None, "stability_certificate": False}
    pair = reconstruct_parent_pair(arrays, record)
    state = state_from_arrays(arrays)
    pins_unchanged(record.get("source_pins"), parent_pins(arrays, state))
    requested = "fft" if backend is None and fft is None else backend
    pair, state, info = episode.resolve_pair(pair, state, fft=fft, backend=requested or "fft")
    if info.get("state_changed") or info.get("W_changed"):
        raise ValueError("FFT resolve changed the frozen state or frame")
    mode = normalize_control_mode(record.get("control_mode", "coupled"))
    stepper = leading.frozen_geometry_step if mode == "frozen_geometry" else leading.rk4_step
    clocks = np.array(arrays["normal_clocks"], dtype=float, copy=True)
    ledger = dict(record["work_ledger"])
    started = time.process_time()
    result = episode.advance_case(
        pair, state, coordinate_time=record["coordinate_time"], steps=record["steps"],
        stations=record["stations"], step_cap=record["step_cap"],
        geometry_frozen=mode == "frozen_geometry", cpu_allowance=cpu_allowance, spent=0.0,
        forecast_factor=forecast_factor, memory_limit_bytes=memory_limit_bytes, max_steps=max_steps,
        stepper=stepper, normal_clocks=clocks, control_mode=mode, backend=requested or "fft",
        work_ledger=ledger, channel_sample=record.get("channel_sample"),
        channel_time=record.get("channel_time"), **kwargs,
    )
    elapsed = time.process_time() - started
    episode.append_observations(directory, case_id, result["observations"])
    chunk = None
    pins_before = record.get("source_pins")
    for snapshot in result["snapshots"]:
        values = dict(arrays)
        for name, key in (("Q", "Q"), ("r", "r"), ("p_Q", "pi_Q"), ("p_r", "pi_r"),
                          ("phi0", "phi0"), ("phi1", "phi1")):
            values[key] = np.array(getattr(snapshot["state"], name), copy=True)
        values["normal_clocks"] = np.array(snapshot["clocks"], dtype=float, copy=True)
        values["clock_rates"] = np.array(snapshot["clock_rates"], dtype=float, copy=True)
        at_end = abs(snapshot["time"] - result["time"]) <= 1e-12
        updated = dict(record)
        updated.update(
            coordinate_time=snapshot["time"], steps=snapshot["steps"],
            stations_reached=snapshot["stations_reached"], work_ledger=snapshot["ledger"],
            channel_sample=snapshot["channel_sample"], channel_time=snapshot["channel_time"],
            snapshot_kind=snapshot["kind"], status=result["status"] if at_end else "station_retained",
            stop=result["stop"] if at_end else None, evolved=True, control_mode=mode,
            projector_source_reset=False, source_reset=False, initial_state_called=False,
        )
        updated["source_pins_before"] = pins_before
        updated["source_pins_after"] = parent_pins(values, snapshot["state"])
        pins_unchanged(pins_before, updated["source_pins_after"])
        if diagnostics and at_end:
            updated["observation"] = scalar_observation(pair, snapshot["state"], snapshot["time"], mode)
        committed = commit_parent_checkpoint(directory, updated, values)
        pins_before = updated["source_pins_after"]
        chunk = committed["npz"]
        record, arrays, state = updated, values, snapshot["state"]
    return {"case_id": case_id, "status": result["status"], "coordinate_time": result["time"],
            "steps": result["steps"], "stations_reached": result["stations_reached"],
            "stop": result["stop"], "child_cpu_seconds": elapsed, "chunk": chunk,
            "stability_certificate": False, "state_clamped": False,
            "physical_instability_claimed": False, "projector_source_reset": False}


def estimate_case_cpu(pair, state, *, step_cap, stations, current_time, probe_step_cpu, probe_diag_cpu):
    dt, _metadata = scaled_step_restriction(pair, state, step_cap)
    span = max(normalize_stations(stations)[-1] - float(current_time), 0.0)
    steps = int(np.ceil(span / dt)) if dt else 0
    diagnostics = int(np.ceil(span / episode.OBSERVATION_CADENCE)) + len(normalize_stations(stations)) + 1
    return {"dt": dt, "estimated_steps": steps,
            "estimated_cpu_seconds": steps * float(probe_step_cpu) + diagnostics * float(probe_diag_cpu)}


@contextmanager
def runtime_adapter():
    """Install scoped hooks and restore them, including after an exception."""
    saved = {
        "episode": {name: getattr(episode, name) for name in (
            "SCHEMA", "pair_from_arrays", "state_from_arrays", "_real_array_names",
            "normalize_control_mode", "normalize_stations", "step_restriction",
            "evolving_step", "frozen_geometry_step", "_clock_rates", "scalar_observation_row",
            "observe", "execute_case", "estimate_case_cpu",
        )},
        "leading_step": leading.step_restriction,
        "make_fft_pair": backend.make_fft_pair,
    }
    episode.SCHEMA = SCHEMA
    episode.pair_from_arrays = reconstruct_parent_pair
    episode.state_from_arrays = state_from_arrays
    episode._real_array_names = real_array_names
    episode.normalize_control_mode = normalize_control_mode
    episode.normalize_stations = normalize_stations
    episode.step_restriction = scaled_step_restriction
    episode.evolving_step = leading.rk4_step
    episode.frozen_geometry_step = leading.frozen_geometry_step
    episode._clock_rates = leading.clock_rates
    episode.scalar_observation_row = scalar_observation_row
    episode.observe = episode_observe
    episode.execute_case = execute_case
    episode.estimate_case_cpu = estimate_case_cpu
    leading.step_restriction = scaled_step_restriction
    backend.make_fft_pair = make_fft_pair
    try:
        with threadpool_limits(limits=1), backend.fft_thread_limit(1):
            yield
    finally:
        for name, value in saved["episode"].items():
            setattr(episode, name, value)
        leading.step_restriction = saved["leading_step"]
        backend.make_fft_pair = saved["make_fft_pair"]


def _set_one_thread():
    for name in _THREAD_NAMES:
        os.environ[name] = "1"


def worker_initializer():
    """Spawn-safe hooks only. Does not evolve a case."""
    global _worker_scope
    _set_one_thread()
    _worker_scope = runtime_adapter()
    _worker_scope.__enter__()


def worker_identity():
    return {
        "schema": episode.SCHEMA,
        "scaled_step_installed": episode.step_restriction is scaled_step_restriction
        and leading.step_restriction is scaled_step_restriction,
        "raw_step_active": leading.step_restriction is legacy_step.RAW_RESTRICTION,
        "fft_pair_installed": backend.make_fft_pair is make_fft_pair,
        "rates_module": leading.rates.__module__,
        "rk4_module": leading.rk4_step.__module__,
        "trajectory_executed": False,
        "pool_launched_by_initializer": False,
    }


def pool(max_workers):
    if not 1 <= int(max_workers) <= MAX_WORKERS:
        raise ValueError("one pool of at most 6 one-thread workers")
    return ProcessPoolExecutor(
        max_workers=int(max_workers), mp_context=multiprocessing.get_context(POOL_START_METHOD),
        initializer=worker_initializer,
    )


def cpu_usage():
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time() + child.ru_utime + child.ru_stime


def run(output=OUTPUT, *, workers=MAX_WORKERS, cpu_budget=CPU_BUDGET_SECONDS, max_steps=None):
    """Root campaign launch. One pool, native FFT, checkpoint on budget stop."""
    if not 1 <= int(workers) <= MAX_WORKERS:
        raise ValueError("one pool of at most 6 one-thread workers")
    if not 0 < float(cpu_budget) <= CPU_BUDGET_SECONDS:
        raise ValueError("aggregate CPU budget including children must be in (0, 21600]")
    manifest = episode.read_manifest(output)
    if manifest.get("schema") != SCHEMA:
        raise ValueError("not a parent episode checkpoint")
    if manifest.get("numerical_binding", {}).get("raw_coordinate_norm_used"):
        raise ValueError("parent episode refuses the raw coordinate step norm")
    started = cpu_usage()
    result = episode.run(
        output, workers=int(workers), cpu_budget_seconds=float(cpu_budget),
        forecast_factor=FORECAST_FACTOR, memory_limit_bytes=MEMORY_LIMIT_BYTES,
        max_steps=max_steps, backend="fft", executor=pool,
    )
    used = cpu_usage() - started
    result["cpu_seconds_including_children"] = used
    result["cpu_accounting"] = "parent_process_time_plus_children_rusage"
    result["forecast_factor"] = FORECAST_FACTOR
    result["fft"] = "native_one_thread"
    result["executor_pools"] = 1
    result["projector_source_reset"] = False
    result["held_out_nonlinear_measurement_before_prediction"] = False
    if used >= float(cpu_budget) or result.get("status") == "BUDGET_STOP":
        result["status"] = "BUDGET_STOP"
        result["checkpoint_retained"] = True
    episode._write_json(episode._manifest_path(output), result)
    return result


def check(output=OUTPUT):
    """Read bindings, hashes, reconstruction, clocks and ledger stocks."""
    before = _directory_files(output)
    manifest = episode.read_manifest(output)
    if manifest.get("schema") != SCHEMA:
        raise ValueError("not a parent episode checkpoint")
    binding = manifest.get("source_hashes") or {}
    for name, digest in binding.items():
        if episode.file_sha256(ROOT / name) != digest:
            raise ValueError("parent episode source changed: " + name)
    if manifest.get("numerical_binding", {}).get("raw_coordinate_norm_used"):
        raise ValueError("record used the raw coordinate step norm")
    problems = []
    checked = []
    for case in manifest["cases"]:
        ordinals = []
        prefix = case["case_id"] + "-"
        for path in sorted(Path(output).glob(prefix + "*.json")):
            ordinals.append(int(path.name[len(prefix):-len(".json")]))
        if not ordinals:
            problems.append(case["case_id"] + " has no commit")
            continue
        previous = None
        for ordinal in ordinals:
            record, arrays = episode.load_checkpoint(output, case["case_id"], ordinal)
            if record.get("schema") != SCHEMA:
                problems.append(case["case_id"] + " schema drifted")
            if "chi" in arrays or "pi_chi" in arrays:
                problems.append(case["case_id"] + " stored an auxiliary chi field")
            if record.get("projector_source_reset") or record.get("source_reset"):
                problems.append(case["case_id"] + " requests a source reset")
            if record.get("held_out_nonlinear_measurement_before_prediction"):
                problems.append(case["case_id"] + " measured a held-out quantity before prediction")
            pair = reconstruct_parent_pair(arrays, record)
            state = state_from_arrays(arrays)
            if state.phi0.shape[1] != pair.column_rank or pair.weights.shape != (pair.column_rank,):
                problems.append(case["case_id"] + " rank mismatch")
            if pair.geometry_map.shape != (pair.grid.ng, pair.grid.ng):
                problems.append(case["case_id"] + " geometry frame is not the full band")
            if not np.array_equal(state.p_Q, arrays["pi_Q"]) or not np.array_equal(state.p_r, arrays["pi_r"]):
                problems.append(case["case_id"] + " geometry momenta were not restored")
            if not np.array_equal(arrays["normal_clocks"], np.array(record.get("normal_clocks_initial", arrays["normal_clocks"]), dtype=float)) and ordinal == 0:
                problems.append(case["case_id"] + " initial clocks differ from the handoff")
            if "work_ledger" not in record:
                problems.append(case["case_id"] + " is missing ledger stocks")
            pins = parent_pins(arrays, state)
            if previous is not None:
                try:
                    pins_unchanged(previous, pins)
                except ValueError as error:
                    problems.append(str(error))
            previous = pins
            chunk = Path(output) / f"{case['case_id']}-{ordinal:06d}.npz"
            sidecar = Path(output) / f"{case['case_id']}-{ordinal:06d}.json"
            if chunk.stat().st_size + sidecar.stat().st_size > CHUNK_LIMIT_BYTES:
                problems.append(chunk.name + " exceeds 64 MiB")
            checked.append(case["case_id"] + f":{ordinal}")
    after = _directory_files(output)
    if before != after:
        problems.append("check wrote or changed campaign files")
    report = {"schema": SCHEMA, "ok": not problems, "checked": checked, "problems": problems,
              "evolved": False, "bytes_written": 0, "pool_launched": False,
              "stability_certificate": False}
    if problems:
        raise ValueError("parent episode check failed: " + "; ".join(problems))
    return report
