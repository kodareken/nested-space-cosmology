"""Stage 0/2 runner for the first separated-pair continuation batch.

Stage 0 freezes the saved nf=128 and nf=256 matched (dt=0.0005) baseline
states at T=0.3 and the replay basis W. It does not call ``initial_state``
or a new constraint solve. The dense grid remains the initializer. Evolution
asks the carrier resolver for ``auto``, which selects the accepted FFT
carrier in this environment. Stage 2 continues six cases toward stations
1 and 3. A frozen-geometry case keeps the controller's zero geometry jets.

The geometric step limit is a sufficient principal-symbol restriction
from the owned conformal equations. Together with the field CFL it is
not a stability certificate. A positive-chart exit keeps the last
admissible state and the step bracket. It is not clamped and it is not
reported as a physical instability.
"""
from __future__ import annotations

from dataclasses import replace
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
import fcntl
import hashlib
import importlib
import inspect
import json
import math
import os
from pathlib import Path
import resource
import sys
import time

import numpy as np

from . import nsc_evolving_reduction as reduction
from . import nsc_nested_parent_child as model
from . import nsc_regional_energy_exchange as regional
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_episode_assessment as metric
from . import nsc_spherical_feedback_action as action
from . import nsc_spherical_galerkin_coupling as galerkin


SCHEMA = "NSC-DISCOVERY-EPISODE-v1"
HANDOFF_TIME = 0.3
DEFAULT_STATIONS = (1.0, 3.0)
ACCEPTED_LATER_STATIONS = (8.0, 16.0, 24.0)
ALLOWED_STATIONS = DEFAULT_STATIONS + ACCEPTED_LATER_STATIONS
STEP_CAP = 0.001
MATCHED_STEP_CAP = 0.0005
DEFAULT_WORKERS = 6
DEFAULT_CPU_BUDGET_SECONDS = 21600.0
DEFAULT_FORECAST_FACTOR = 1.5
DEFAULT_MEMORY_BYTES = 8 * 1024 ** 3
CHUNK_LIMIT_BYTES = 64 * 1024 ** 2
RK4_HALF_STABILITY = 1.4
OBSERVATION_CADENCE = 0.05
WORK_CHANNELS = (
    "coordinate_fieldwork", "pressure_work", "lapse_work", "boundary_child", "boundary_parent",
)


def _canonical_cap(value):
    value = float(value)
    for cap in (STEP_CAP, MATCHED_STEP_CAP):
        if abs(value - cap) <= 1e-15:
            return cap
    raise ValueError("saved baseline step cap is not 0.001 or 0.0005")
EVOLUTION_BACKEND = "dense_galerkin"
REDUCTION_BACKEND = "dense"
CANONICAL_PI = "canonical_pi"
NODAL_MOMENTUM = "nodal_momentum"

_MODULE = Path(__file__).resolve()
LAB = _MODULE.parents[2]
REPO = _MODULE.parents[3]
V1_JSON = LAB / "results/development/nsc-nested-parent-child-v1.json"
V1_NPZ = LAB / "results/development/nsc-nested-parent-child-v1.npz"
BASIS_JSON = LAB / "results/development/nsc-nested-parent-child-replay-basis-v1.json"
BASIS_NPZ = LAB / "results/development/nsc-nested-parent-child-replay-basis-v1.npz"
BASELINE_CASES = (
    "nf128_baseline_dt0.001",
    "nf128_baseline_dt0.0005",
    "nf256_baseline_dt0.001",
    "nf256_baseline_dt0.0005",
)
_REAL_ARRAYS = ("Q", "r", "chi", "pi_Q", "pi_r", "pi_chi", "source_weights", "clock_rates", "normal_clocks", "W")
_COMPLEX_ARRAYS = ("phi0", "phi1", "source_phi0", "source_phi1", "observer_columns")
_PACKAGE_DIR = Path(model.__file__).resolve().parent

_interface_cache = {}


def _cap_token(step_cap):
    cap = _canonical_cap(step_cap)
    return "0.001" if cap == STEP_CAP else "0.0005"


def first_batch_specs(include_frozen=True):
    """Six cases. Both caps and the frozen control use the dt=0.0005 T=0.3 state."""
    specs = []
    for fermions in (128, 256):
        parent = f"nf{fermions}_baseline_dt0.0005"
        for cap in (STEP_CAP, MATCHED_STEP_CAP):
            specs.append({
                "case_id": f"nf{fermions}_coupled_dt{_cap_token(cap)}",
                "parent_case": parent, "nf": fermions, "control_mode": "coupled",
                "geometry": "evolving", "step_cap": cap,
            })
        if include_frozen:
            specs.append({
                "case_id": f"nf{fermions}_frozen_geometry_dt0.0005",
                "parent_case": parent, "nf": fermions, "control_mode": "frozen_geometry",
                "geometry": "frozen", "step_cap": MATCHED_STEP_CAP,
            })
    return tuple(specs)


def baseline_case_ids(include_frozen=True):
    """Admission order for the first batch."""
    return tuple(spec["case_id"] for spec in first_batch_specs(include_frozen=include_frozen))


def normalize_stations(stations=None):
    if stations is None:
        values = DEFAULT_STATIONS
    else:
        values = tuple(float(value) for value in stations)
    allowed = set(ALLOWED_STATIONS)
    if not values or any(value not in allowed for value in values):
        raise ValueError("stations must be chosen from 1, 3 and the accepted later stations 8, 16, 24")
    if list(values) != sorted(values) or len(set(values)) != len(values):
        raise ValueError("stations must be strictly increasing")
    if values[0] <= HANDOFF_TIME:
        raise ValueError("stations must lie after the T=0.3 handoff")
    return values


def assert_campaign_output(path):
    """Allow a new nsc-discovery successor, or any directory outside the repository.

    Existing sealed records are never overwritten. A repository-wide write
    ban is not the rule. Chunk publication still uses exclusive creates.
    """
    resolved = Path(path).expanduser().resolve()
    sealed = {V1_JSON.resolve(), V1_NPZ.resolve(), BASIS_JSON.resolve(), BASIS_NPZ.resolve()}
    if resolved in sealed or any(parent in sealed for parent in resolved.parents):
        raise PermissionError("refusing to write a sealed scientific record: " + str(resolved))
    if resolved.is_file():
        raise FileExistsError("refusing to overwrite an existing file: " + str(resolved))
    inside = resolved == REPO or REPO in resolved.parents
    if not inside:
        return resolved
    development = (LAB / "results" / "development").resolve()
    if resolved != development and development not in resolved.parents:
        raise PermissionError(
            "repository campaign output must be a new successor under lab/results/development/nsc-discovery-"
        )
    relative = resolved.relative_to(development)
    if not any(part.startswith("nsc-discovery-") for part in relative.parts):
        raise PermissionError(
            "new repository evidence must live under a nsc-discovery- successor name: " + str(resolved)
        )
    return resolved


def assert_output_outside_repository(path):
    """Backward-compatible name. Successor discovery output inside the lab is allowed."""
    return assert_campaign_output(path)


def file_sha256(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def array_sha256(array):
    contiguous = np.ascontiguousarray(array)
    hasher = hashlib.sha256()
    hasher.update(str(contiguous.dtype).encode())
    hasher.update(str(contiguous.shape).encode())
    hasher.update(contiguous.tobytes())
    return hasher.hexdigest()


def state_sha256(state):
    """Same byte recipe as the saved nested-pair driver."""
    hasher = hashlib.sha256()
    for name in model.STATE_NAMES:
        array = np.ascontiguousarray(getattr(state, name))
        hasher.update(name.encode())
        hasher.update(str(array.dtype).encode())
        hasher.update(str(array.shape).encode())
        hasher.update(array.tobytes())
    return hasher.hexdigest()


def jsonable(value):
    if isinstance(value, np.ndarray):
        if np.iscomplexobj(value):
            return {"real": value.real.tolist(), "imag": value.imag.tolist()}
        return value.tolist()
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return value


def resident_peak_bytes():
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    if sys.platform == "darwin":
        return value
    return value * 1024


def reduction_backend_status():
    """The full-state step is the evolution. Dense reduction is selectable."""
    accepted = reduction._backend_name(REDUCTION_BACKEND)
    return {
        "evolution_backend": EVOLUTION_BACKEND,
        "reduction_backend": accepted,
        "reduction_backend_used": False,
        "streamed_reduction_selected": False,
        "dense_exterior_propagator_stored": False,
    }


def _discover_callable(function_name, parameter_names):
    cached = _interface_cache.get(function_name)
    if cached is not None:
        return cached
    notes = []
    found = None
    module_name = None
    for path in sorted(_PACKAGE_DIR.glob("nsc_*.py")):
        if path.resolve() == _MODULE:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except OSError as error:
            notes.append(path.name + ": " + str(error))
            continue
        if ("def " + function_name + "(") not in text:
            continue
        import_name = "recursive_horizons." + path.stem
        try:
            module = importlib.import_module(import_name)
        except Exception as error:
            notes.append(import_name + " import failed: " + type(error).__name__ + ": " + str(error))
            continue
        candidate = getattr(module, function_name, None)
        if not callable(candidate):
            notes.append(import_name + " does not export " + function_name)
            continue
        try:
            signature = inspect.signature(candidate)
        except (TypeError, ValueError) as error:
            notes.append(import_name + " signature unreadable: " + str(error))
            continue
        names = [parameter.name for parameter in signature.parameters.values()
                 if parameter.kind in (parameter.POSITIONAL_ONLY, parameter.POSITIONAL_OR_KEYWORD)]
        if names[:len(parameter_names)] != list(parameter_names):
            notes.append(import_name + " signature " + str(names) + " does not start with " + str(parameter_names))
            continue
        found = candidate
        module_name = import_name
        break
    record = {"function": found, "module": module_name, "notes": notes}
    _interface_cache[function_name] = record
    return record


def resolve_pair(pair, state, *, fft=None, backend=None):
    """Resolve carriers. ``auto`` is the evolution default; the initializer stays dense.

    ``fft=True`` and ``fft=False`` still force the FFT carrier or the dense
    oracle. The saved state bytes and W are not rebuilt. ``initial_state``
    is not called.
    """
    if fft is True:
        requested = "fft"
    elif fft is False:
        requested = "dense"
    elif backend is None:
        requested = "auto"
    else:
        requested = str(backend)
    info = {
        "resolver": "explicit_adapter",
        "fft_requested": requested == "fft",
        "fft_applied": False,
        "backend_requested": requested,
        "backend": requested,
        "initializer": "dense_galerkin_grid",
        "state_changed": False,
        "W_changed": False,
        "initial_state_called": False,
        "notes": [],
    }
    try:
        from .nsc_discovery_backend import resolve_pair as backend_resolve
        from .nsc_discovery_backend import select_backend
    except Exception as error:
        backend_resolve = None
        select_backend = None
        info["notes"].append("backend import failed: " + type(error).__name__ + ": " + str(error))
    if backend_resolve is not None:
        try:
            selected = select_backend(requested)
            resolved = backend_resolve(pair, backend=requested)
        except Exception as error:
            info["notes"].append("resolve_pair failed: " + type(error).__name__ + ": " + str(error))
            return pair, state, info
        if not np.array_equal(np.asarray(resolved.geometry_map), np.asarray(pair.geometry_map)):
            info["W_changed"] = True
            info["notes"].append("resolver changed W; original pair retained")
            return pair, state, info
        info["resolver"] = "recursive_horizons.nsc_discovery_backend"
        info["backend"] = selected
        info["fft_applied"] = selected == "fft"
        info["fft_requested"] = selected == "fft" or requested == "fft"
        return resolved, state, info
    discovered = _discover_callable("resolve_pair", ("pair", "state"))
    info["notes"].extend(discovered["notes"])
    if discovered["function"] is None:
        return pair, state, info
    try:
        returned = discovered["function"](pair, state, fft=requested == "fft")
    except Exception as error:
        info["notes"].append("resolve_pair failed: " + type(error).__name__ + ": " + str(error))
        return pair, state, info
    if not isinstance(returned, tuple) or len(returned) < 2:
        info["notes"].append("pair-state resolver return was not accepted")
        return pair, state, info
    new_pair, new_state = returned[0], returned[1]
    same_state = all(np.array_equal(getattr(new_state, name), getattr(state, name)) for name in model.STATE_NAMES)
    same_w = np.array_equal(np.asarray(new_pair.geometry_map), np.asarray(pair.geometry_map))
    if not same_state or not same_w:
        info["state_changed"] = not same_state
        info["W_changed"] = not same_w
        info["notes"].append("resolver changed the frozen state or W; original bytes retained")
        return pair, state, info
    info["resolver"] = discovered["module"]
    info["fft_applied"] = requested == "fft"
    info["backend"] = "fft" if requested == "fft" else "dense"
    return new_pair, new_state, info


def _readonly(value, dtype=None):
    result = np.array(value, dtype=dtype, copy=True)
    result.setflags(write=False)
    return result


def quadrature_principal_wavenumber(grid):
    """Largest nonzero symbol of the owned periodic derivative.

    That operator sets the Nyquist symbol to zero. The principal part of
    the conformal geometric rates is evaluated with this derivative on the
    quadrature grid before projection.
    """
    max_mode = int(grid.nq) // 2 - 1
    if max_mode < 1:
        raise ValueError("quadrature derivative has no positive frequency")
    return float(2 * np.pi * max_mode / grid.length)


def geometry_band_wavenumber(grid):
    return float(2 * np.pi * np.max(np.abs(grid.modes_g)) / grid.length)


def conformal_principal_speeds(system, fine_state):
    """Coordinate speeds from the owned conformal principal symbol.

    With Z=3a the gauge-fixed second-order principal part on (Q, r, chi)
    is ∂_tt = ∂_xx, so the identity speed is 1. The light-cone speed
    |L|/Q + |β| is retained. If the coefficient identity is not present,
    the step uses the absolute majorant of the uncancelled symbol instead.
    """
    if system is None or fine_state is None:
        raise ValueError("principal speed needs the stage-local conformal system")
    radius = np.asarray(fine_state.r, dtype=float)
    radial = np.asarray(fine_state.Q, dtype=float)
    if np.min(radial) <= 0 or np.min(radius) <= 0:
        raise coupling.PositiveChartExit("nonpositive_chart_for_principal_speed", np.nan, fine_state.copy())
    force_r, f_chi = action.partial_F(1.0, system.A, system.C_W)
    coefficient_a = float(force_r)
    f_value = float(f_chi)
    metric_z = float(action.feedback_Z(system.A))
    identity_gap = abs(coefficient_a - metric_z / 3.0)
    identity_tolerance = 8 * float(np.spacing(max(1.0, abs(metric_z))))
    light = float(np.max(np.abs(system.length_density) / radial + np.abs(system.shift)))
    identity_speed = max(light, 1.0)
    ratio_r = np.max(np.abs(radius) / radial)
    ratio_a = abs(coefficient_a) / max(abs(metric_z), 1e-300)
    row_r = 1.0 + 2.0 * ratio_a * ratio_r
    row_chi = 1.0 + abs(metric_z) * ratio_r / max(3.0 * abs(f_value), 1e-300)
    row_chi += abs(coefficient_a) * ratio_r / max(abs(f_value), 1e-300)
    row_chi += 2.0 * (coefficient_a ** 2) * (ratio_r ** 2) / max(abs(f_value * metric_z), 1e-300)
    majorant_speed = float(np.sqrt(max(1.0, row_r, row_chi)))
    identity_holds = bool(identity_gap <= identity_tolerance and f_value != 0.0 and metric_z != 0.0)
    speed = identity_speed if identity_holds else max(identity_speed, majorant_speed)
    return {
        "speed": float(speed),
        "identity_speed": float(identity_speed),
        "light_cone_speed": light,
        "majorant_speed": majorant_speed,
        "majorant_used": not identity_holds,
        "identity_gap_a_minus_Z_over_3": float(identity_gap),
        "identity_tolerance": identity_tolerance,
        "identity_holds": identity_holds,
        "stability_certificate": False,
    }


def geometric_principal_frequency(pair, state):
    nodal = model.reconstruct_state(pair, state)
    fine = galerkin.prolong_state(pair.grid, nodal)
    system = galerkin.active_fine_system(pair.grid, fine)
    speeds = conformal_principal_speeds(system, fine)
    k_quadrature = quadrature_principal_wavenumber(pair.grid)
    k_band = geometry_band_wavenumber(pair.grid)
    k_used = max(k_quadrature, k_band)
    omega = float(speeds["speed"] * k_used)
    speeds.update(wavenumber_quadrature=k_quadrature, wavenumber_geometry_band=k_band,
                  wavenumber_used=k_used, omega=omega,
                  principal_part="owned conformal gauge-fixed second-order symbol",
                  field_cfl_alone=False, stability_certificate=False)
    return speeds


def field_frequency(pair, state):
    nodal = model.reconstruct_state(pair, state)
    omega, quadrature_omega = galerkin.subspace_frequency(pair.grid, nodal)
    return {"omega": float(omega), "quadrature_nyquist_report": float(quadrature_omega),
            "source": "galerkin.subspace_frequency"}


def step_restriction(pair, state, step_cap):
    """Field CFL and the geometric principal restriction. Not a certificate."""
    cap = float(step_cap)
    if not np.isfinite(cap) or cap <= 0:
        raise ValueError("a positive finite step cap is required")
    field = field_frequency(pair, state)
    geometric = geometric_principal_frequency(pair, state)
    omega = max(float(field["omega"]), float(geometric["omega"]), 1.0)
    dt = min(cap, RK4_HALF_STABILITY / omega)
    field_only = min(cap, RK4_HALF_STABILITY / max(float(field["omega"]), 1.0))
    return dt, {
        "dt": dt,
        "step_cap": cap,
        "field": field,
        "geometric": geometric,
        "combined_omega": omega,
        "field_only_dt": field_only,
        "field_cfl_alone": False,
        "geometric_restriction_applied": True,
        "stability_certificate": False,
        "zero_residual_gate": False,
        "max_T_gate": False,
    }


def _scalar_summary(values):
    array = np.asarray(values, dtype=float)
    return {"min": float(np.min(array)), "max": float(np.max(array)),
            "l2": float(np.sqrt(np.mean(array ** 2)))}


def _actual_curvature(pair, nodal, rate, bundle):
    """Actual metric R_h and Weyl scalar. Chi is not substituted."""
    grid = pair.grid
    fine = bundle["fine_state"]
    _force_r, f_chi = action.partial_F(fine.r, grid.fine.A, grid.fine.C_W)
    f_chi = float(f_chi)
    q_dot = galerkin.prolong_geometry(grid, rate.Q)
    p_chi_dot = galerkin.prolong_geometry(grid, rate.p_chi)
    directional = (q_dot * fine.p_chi + fine.Q * p_chi_dot) / (2 * f_chi)
    q_ddot = galerkin.prolong_geometry(grid, galerkin.pull_geometry(grid, directional))
    zero = np.zeros(grid.nq)
    curvature = metric.direct_rh_grid(
        fine.Q, fine.Q, zero, q_dot, q_ddot, grid.length, L_dot=q_dot, beta_dot=zero,
    )["R_h"]
    weyl = metric.weyl_actual(curvature, fine.r)
    proxy = metric.weyl_proxy_from_chi(fine.chi, fine.r)
    shell_gap = curvature - fine.chi - 2.0
    return {
        "actual_metric_curvature": {
            "R_h": _scalar_summary(curvature),
            "weyl_C2": _scalar_summary(weyl),
            "chi_substituted": False,
            "continuum_curvature_certified": False,
        },
        "chi_shell": {
            "proxy": _scalar_summary(proxy),
            "R_h_minus_chi_minus_2": _scalar_summary(shell_gap),
            "actual_minus_proxy": _scalar_summary(weyl - proxy),
            "substituted_into_curvature": False,
            "meaning": "auxiliary shell chi = R_h - 2 is a comparison, not the curvature",
        },
    }


def normalize_control_mode(control_mode):
    if control_mode in ("frozen", "frozen_geometry"):
        return "frozen_geometry"
    if control_mode in ("coupled", "evolving", None):
        return "coupled"
    raise ValueError("control_mode must be coupled or frozen_geometry")


def actual_control_rate(pair, state, control_mode="coupled"):
    """The rate the controller steps. Frozen geometry zeros Q, r, chi, momenta, and L.

    Field column rates and the source forces stay. Metric coordinate work is
    zero in the frozen mode. This is not the coupled ODE reconstructed after
    the fact.
    """
    mode = normalize_control_mode(control_mode)
    nodal = model.reconstruct_state(pair, state)
    coupled, bundle = galerkin.compose_fine_hamiltonian(pair.grid, nodal)
    if mode == "coupled":
        return nodal, coupled, bundle, mode
    zero_g = np.zeros_like(coupled.Q)
    zero_q = np.zeros(pair.grid.nq, dtype=float)
    rate = coupling.CauchyRate(
        zero_g, zero_g.copy(), zero_g.copy(), zero_g.copy(), zero_g.copy(), zero_g.copy(),
        np.array(coupled.phi0, copy=True), np.array(coupled.phi1, copy=True),
        np.array(coupled.force_L, copy=True), np.array(coupled.force_Q, copy=True),
        np.array(coupled.force_beta, copy=True), 0.0,
    )
    frozen_bundle = dict(bundle)
    frozen_bundle["lifted_Q"] = zero_q
    frozen_bundle["lapse_dot"] = zero_q.copy()
    frozen_bundle["unprojected_Q"] = zero_q.copy()
    frozen_bundle["unprojected_lapse_dot"] = zero_q.copy()
    frozen_bundle["shift_dot"] = zero_q.copy()
    return nodal, rate, frozen_bundle, mode


def _frozen_curvature(pair, fine):
    """Spatial metric jets with the frozen controller's Qdot = Qddot = Ldot = 0."""
    zero = np.zeros(pair.grid.nq, dtype=float)
    curvature = metric.direct_rh_grid(
        fine.Q, fine.Q, zero, zero, zero, pair.grid.length, L_dot=zero, beta_dot=zero,
    )["R_h"]
    weyl = metric.weyl_actual(curvature, fine.r)
    proxy = metric.weyl_proxy_from_chi(fine.chi, fine.r)
    shell_gap = curvature - fine.chi - 2.0
    return {
        "actual_metric_curvature": {
            "R_h": _scalar_summary(curvature),
            "weyl_C2": _scalar_summary(weyl),
            "Q_dot_max": 0.0,
            "Q_ddot_max": 0.0,
            "L_dot_max": 0.0,
            "chi_substituted": False,
            "continuum_curvature_certified": False,
        },
        "chi_shell": {
            "proxy": _scalar_summary(proxy),
            "R_h_minus_chi_minus_2": _scalar_summary(shell_gap),
            "actual_minus_proxy": _scalar_summary(weyl - proxy),
            "substituted_into_curvature": False,
            "meaning": "auxiliary shell chi = R_h - 2 is a comparison, not the curvature",
        },
    }


def stability_diagnostics(pair, state, time, *, control_mode="coupled", rate_bundle=None):
    """Separate curvature, chi-shell, constraints, CAR, and energy/work/chart.

    Frozen mode uses the controller rate. It is not classified as the coupled
    evolution. No single score is formed and nothing here is a certificate.
    """
    mode = normalize_control_mode(control_mode)
    if rate_bundle is None:
        nodal, rate, bundle, mode = actual_control_rate(pair, state, mode)
    else:
        rate, bundle = rate_bundle
        nodal = model.reconstruct_state(pair, state)
    fine = bundle["fine_state"]
    system = bundle["fine_system"]
    chart_reason = coupling.chart_failure(system, fine)
    if mode == "frozen_geometry":
        curvature = _frozen_curvature(pair, fine)
    else:
        curvature = _actual_curvature(pair, nodal, rate, bundle)
        curvature["actual_metric_curvature"]["Q_dot_max"] = float(np.max(np.abs(
            galerkin.prolong_geometry(pair.grid, rate.Q))))
        curvature["actual_metric_curvature"]["L_dot_max"] = float(np.max(np.abs(bundle["lapse_dot"])))
    constraints = galerkin.constraint_diagnostics(pair.grid, nodal, fine_state=fine)
    transport = galerkin.conformal_constraint_transport(pair.grid, nodal, rate, bundle)
    columns = np.vstack((np.array(state.phi0, copy=True), np.array(state.phi1, copy=True)))
    gram = columns.conj().T @ columns
    roots = np.sqrt(np.asarray(pair.weights, dtype=float))
    eigenvalues = np.linalg.eigvalsh(roots[:, None] * gram * roots[None, :])
    accounting = model.energy_accounting(pair, state)
    ledger = regional.matter_ledger(system, fine)
    lifted = coupling.CauchyRate(
        galerkin.prolong_geometry(pair.grid, rate.Q),
        galerkin.prolong_geometry(pair.grid, rate.r),
        galerkin.prolong_geometry(pair.grid, rate.chi),
        galerkin.prolong_geometry(pair.grid, rate.p_Q),
        galerkin.prolong_geometry(pair.grid, rate.p_r),
        galerkin.prolong_geometry(pair.grid, rate.p_chi),
        galerkin.prolong_columns(pair.grid, rate.phi0),
        galerkin.prolong_columns(pair.grid, rate.phi1),
        np.array(rate.force_L, copy=True), np.array(rate.force_Q, copy=True),
        np.array(rate.force_beta, copy=True), float(rate.fieldwork_power),
    )
    balance = regional.proper_balance_terms(system, fine, lifted, ledger)
    metric_work = 0.0 if mode == "frozen_geometry" else float(rate.fieldwork_power)
    report = {
        "time": float(time),
        "control_mode": mode,
        "coupled_stability_classification": False if mode == "frozen_geometry" else None,
        "actual_metric_curvature": curvature["actual_metric_curvature"],
        "chi_shell": curvature["chi_shell"],
        "projected_constraints": {
            "projected_hamilton_max": constraints["projected_hamilton_max"],
            "projected_momentum_max": constraints["projected_momentum_max"],
            "held_out_hamilton_max": constraints["held_out_hamilton_max"],
            "held_out_momentum_max": constraints["held_out_momentum_max"],
            "full_hamilton_max": constraints["full_hamilton_max"],
            "full_momentum_max": constraints["full_momentum_max"],
            "gates_continuation": False,
        },
        "CAR": {
            "meaning": "preserved carrier frame: column Gram and occupation-weighted eigenvalues",
            "gram_gap": float(np.max(np.abs(gram - np.eye(6)))),
            "eigenvalue_min": float(eigenvalues[0]),
            "eigenvalue_max": float(eigenvalues[-1]),
            "weights": np.asarray(pair.weights, dtype=float).tolist(),
            "admissible": bool(eigenvalues[0] >= -1e-8 and eigenvalues[-1] <= 1 + 1e-8),
        },
        "constraint_algebra_residual": {
            "constraint_norm": transport["constraint_norm"],
            "forcing_norm": transport["forcing_norm"],
            "quadrature_forcing_norm": transport["quadrature_forcing_norm"],
            "projection_forcing_norm": transport["projection_forcing_norm"],
            "uses_control_rate": True,
            "continuum_constraint_certified": False,
            "time_integral_certified": False,
        },
        "energy": {
            "total": float(model.energy(pair, state)),
            "field": float(accounting["field"]),
            "gravity": float(accounting["gravity"]),
            "field_closure_error": float(accounting["field_closure_error"]),
            "normal_shell_max": float(np.max(np.abs(ledger["normal_energy_nodal"]))),
        },
        "work": {
            "metric_coordinate_work": metric_work,
            "fieldwork_power": metric_work,
            "spatial_lapse_work_max": float(np.max(np.abs(balance["momentum_lapse_work"]))),
            "source_force_Q_max": float(np.max(np.abs(rate.force_Q))),
            "source_force_L_max": float(np.max(np.abs(rate.force_L))),
            "field_rate_max": float(np.max(np.abs(rate.phi0))),
        },
        "chart": {
            "admissible": chart_reason is None,
            "reason": chart_reason,
            "Q_min": float(np.min(fine.Q)),
            "r_min": float(np.min(fine.r)),
            "clamped": False,
        },
        "collapsed_score": None,
        "stability_certificate": False,
    }
    return report


def _observer_arguments(function, offered):
    """Pass only keyword names the current observer atlas accepts."""
    try:
        parameters = inspect.signature(function).parameters
    except (TypeError, ValueError):
        return {}
    arguments = {}
    for name, value in offered.items():
        if name in parameters:
            arguments[name] = value
    if "control_rate" in parameters and "control_rate" not in arguments:
        arguments["control_rate"] = offered.get("nodal_rate")
    return arguments


def observe(pair, state, time, *, control_mode="coupled"):
    """Call the atlas observer with the actual control rate when it accepts those names."""
    mode = normalize_control_mode(control_mode)
    _nodal, rate, bundle, mode = actual_control_rate(pair, state, mode)
    row = {
        "time": float(time),
        "control_mode": mode,
        "stability": stability_diagnostics(pair, state, time, control_mode=mode, rate_bundle=(rate, bundle)),
        "observer_module": "explicit_adapter",
        "external_observation": None,
        "notes": [],
    }
    discovered = _discover_callable("observe", ("pair", "state", "time"))
    row["notes"] = list(discovered["notes"])
    if discovered["function"] is None:
        return row
    offered = {
        "capture_profiles": False,
        "bundle": bundle,
        "nodal_rate": rate,
        "preparation": "separated_pair_v1",
        "control_mode": mode,
    }
    arguments = _observer_arguments(discovered["function"], offered)
    try:
        row["external_observation"] = discovered["function"](pair, state, float(time), **arguments)
        row["observer_module"] = discovered["module"]
        row["observer_arguments"] = sorted(arguments)
    except TypeError:
        try:
            row["external_observation"] = discovered["function"](pair, state, float(time))
            row["observer_module"] = discovered["module"]
            row["notes"].append("observer was called without the current control keywords")
        except Exception as error:
            row["notes"].append("observe failed: " + type(error).__name__ + ": " + str(error))
    except Exception as error:
        row["notes"].append("observe failed: " + type(error).__name__ + ": " + str(error))
    return row


def _combine(state, rate, factor):
    return model.NestedState(*(getattr(state, name) + factor * getattr(rate, name)
                               for name in model.STATE_NAMES))


def _freeze_geometry_rate(rate):
    zero = np.zeros_like(rate.Q)
    return model.NestedRate(zero, zero, zero, zero, zero, zero, np.array(rate.phi0, copy=True),
                            np.array(rate.phi1, copy=True), 0.0, np.array(rate.force_L, copy=True),
                            np.array(rate.force_Q, copy=True), np.array(rate.force_beta, copy=True))


def _chart_check(pair, state):
    nodal = model.reconstruct_state(pair, state)
    fine = galerkin.prolong_state(pair.grid, nodal)
    system = galerkin.active_fine_system(pair.grid, fine)
    reason = coupling.chart_failure(system, fine)
    if reason is not None:
        raise coupling.PositiveChartExit(reason, np.nan, state.copy())
    return fine, system


def frozen_geometry_step(pair, state, dt):
    """Field RK4 at the fixed T=0.3 geometry. The observer is the pair's."""
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("a positive finite timestep is required")

    def stage(current):
        return _freeze_geometry_rate(model.rates(pair, current))

    k1 = stage(state)
    k2 = stage(_combine(state, k1, dt / 2))
    k3 = stage(_combine(state, k2, dt / 2))
    k4 = stage(_combine(state, k3, dt))
    out = model.NestedState(*(getattr(state, name) + dt / 6 * (
        getattr(k1, name) + 2 * getattr(k2, name) + 2 * getattr(k3, name) + getattr(k4, name))
        for name in model.STATE_NAMES))
    _chart_check(pair, out)
    return out


def evolving_step(pair, state, dt):
    return model.rk4_step(pair, state, dt)


def residual_blocks_continuation(projected_hamilton_max, projected_momentum_max):
    """There is no zero-residual gate. The arguments are recorded only."""
    del projected_hamilton_max, projected_momentum_max
    return False


def _chart_exit_record(reason, time_value, dt):
    start = float(time_value)
    return {
        "stop_class": "chart_exit",
        "reason": str(reason),
        "last_admissible_time": start,
        "event_bracket": [start, start + float(dt)],
        "state_clamped": False,
        "physical_instability_claimed": False,
        "failed_state_retained": False,
        "stability_certificate": False,
    }


def _clock_rates(pair, state):
    """Normal-clock rates dτ/dT = r Q at the fixed observers. A fresh array, not a view."""
    return np.array(model.metrics(pair, state)["clock_rates"], dtype=np.float64, copy=True)


def next_observation_time(time_value, cadence, origin=0.0):
    """Next coordinate time on the cadence grid strictly after ``time_value``."""
    step = float(cadence)
    if not np.isfinite(step) or step <= 0.0:
        return math.inf
    mark = float(time_value)
    index = math.floor((mark - float(origin)) / step + 1e-9) + 1
    return float(origin) + index * step


def empty_work_ledger(control_mode="coupled"):
    ledger = {name: 0.0 for name in WORK_CHANNELS}
    ledger.update(quadrature="trapezoid_coordinate_time",
                  control_mode=normalize_control_mode(control_mode),
                  dense_propagator_stored=False)
    return ledger


def accumulate_work_ledger(ledger, left, right, dt):
    """Trapezoid between two control-rate samples. ``dt`` is the sample spacing."""
    if dt < 0.0:
        raise ValueError("work quadrature requires a forward interval")
    updated = dict(ledger)
    for name in WORK_CHANNELS:
        updated[name] = float(ledger[name]) + 0.5 * float(dt) * (float(left[name]) + float(right[name]))
    updated["quadrature"] = "trapezoid_coordinate_time"
    updated["dense_propagator_stored"] = False
    return updated


def _finite_channel(value, fallback=None):
    if value is None:
        return fallback
    try:
        number = float(value)
    except (TypeError, ValueError):
        return fallback
    if not np.isfinite(number):
        return fallback
    return number


def scalar_observation_row(pair, state, time_value, control_mode="coupled"):
    """Compact physical row. Profiles and propagators are not retained."""
    mode = normalize_control_mode(control_mode)
    started = time.process_time()
    full = observe(pair, state, float(time_value), control_mode=mode)
    elapsed = time.process_time() - started
    external = full.get("external_observation") or {}
    stability = full["stability"]
    child = external.get("child_window") or {}
    parent = external.get("parent_window") or {}
    localization = external.get("localization") or {}
    proper = external.get("proper") or {}
    curvature = external.get("curvature") or {}
    shell = stability["chi_shell"]["R_h_minus_chi_minus_2"]
    aux = _finite_channel(curvature.get("chi_shell_gap_max"), max(abs(shell["min"]), abs(shell["max"])))
    pressure = None
    lapse = None
    if child or parent:
        pressure = float(child.get("proper_pressure_work") or 0.0) + float(parent.get("proper_pressure_work") or 0.0)
        lapse = float(child.get("momentum_lapse_work") or 0.0) + float(parent.get("momentum_lapse_work") or 0.0)
    fieldwork = _finite_channel(
        external.get("coordinate_metric_work_total"), stability["work"]["metric_coordinate_work"],
    )
    row = {
        "time": float(time_value),
        "control_mode": mode,
        "localization_peak_x": localization.get("peak_x"),
        "localization_width": localization.get("proper_width"),
        "child_proper_length": proper.get("child_proper_length"),
        "parent_proper_length": proper.get("parent_proper_length"),
        "parent_annulus_proper_length": proper.get("parent_annulus_proper_length"),
        "child_normal_energy": child.get("normal_energy"),
        "parent_normal_energy": parent.get("normal_energy"),
        "child_probability": child.get("source_probability"),
        "parent_probability": parent.get("source_probability"),
        "coordinate_fieldwork": fieldwork,
        "pressure_work": pressure,
        "lapse_work": lapse,
        "boundary_child": child.get("normal_boundary_flux"),
        "boundary_parent": parent.get("normal_boundary_flux"),
        "R_h_max": _finite_channel(curvature.get("R_h_max"), stability["actual_metric_curvature"]["R_h"]["max"]),
        "weyl_C2_max": _finite_channel(
            curvature.get("weyl_C2_max"), stability["actual_metric_curvature"]["weyl_C2"]["max"],
        ),
        "aux_discrepancy_max": aux,
        "car_gram_gap": stability["CAR"]["gram_gap"],
        "car_eigenvalue_min": stability["CAR"]["eigenvalue_min"],
        "car_eigenvalue_max": stability["CAR"]["eigenvalue_max"],
        "energy_total": stability["energy"]["total"],
        "energy_field": stability["energy"]["field"],
        "energy_gravity": stability["energy"]["gravity"],
        "dense_propagator_stored": False,
        "observation_cpu_seconds": float(elapsed),
        "quadrature": "trapezoid_coordinate_time",
    }
    if mode == "frozen_geometry":
        row["coordinate_fieldwork"] = 0.0
        row["coupled_stability_classification"] = False
    if full.get("notes"):
        row["notes"] = list(full["notes"])
    missing = [name for name in WORK_CHANNELS if row.get(name) is None]
    if missing:
        row["channels_incomplete"] = missing
    return row


def advance_case(pair, state, *, coordinate_time, steps, stations, step_cap, geometry_frozen,
                 cpu_allowance, spent, forecast_factor, memory_limit_bytes, max_steps=None,
                 stepper=None, cpu_per_step=None, rss_bytes=None, resolve=None, fft=None,
                 backend=None, normal_clocks=None, control_mode=None,
                 observation_cadence=OBSERVATION_CADENCE, forced_dt=None, until_time=None,
                 work_ledger=None, channel_sample=None, channel_time=None):
    """Advance one case. Stations and the last admissible state are retained.

    Scalar rows are taken on the coordinate cadence, at stations, and at a
    stop. Work integrals use the trapezoid of those control-rate samples.
    """
    stations = normalize_stations(stations)
    mode = normalize_control_mode(control_mode or ("frozen_geometry" if geometry_frozen else "coupled"))
    if stepper is None:
        stepper = frozen_geometry_step if mode == "frozen_geometry" else evolving_step
    if resolve is None:
        resolve = resolve_pair
    if rss_bytes is None:
        rss_bytes = resident_peak_bytes
    current = state.copy()
    mark = float(coordinate_time)
    taken = int(steps)
    start_steps = taken
    start_time = mark
    used = float(spent)
    step_cpu = 0.0
    observation_cpu = 0.0
    last_cost = None if cpu_per_step is None else float(cpu_per_step)
    reached = []
    restriction = None
    clocks = np.zeros(3, dtype=np.float64) if normal_clocks is None else np.array(normal_clocks, dtype=np.float64, copy=True)
    rate_sample = _clock_rates(pair, current)
    ledger = empty_work_ledger(mode) if work_ledger is None else dict(work_ledger)
    ledger["control_mode"] = mode
    ledger["dense_propagator_stored"] = False
    observations = []
    snapshots = []
    previous_row = None if channel_sample is None else dict(channel_sample)
    previous_time = None if channel_time is None else float(channel_time)
    cadence = float(observation_cadence)

    def remember(kind):
        if snapshots and abs(snapshots[-1]["time"] - mark) <= 1e-12:
            if kind == "station":
                snapshots[-1]["kind"] = "station"
            return
        snapshots.append({
            "kind": kind, "time": float(mark), "state": current.copy(),
            "clocks": clocks.copy(), "clock_rates": rate_sample.copy(),
            "ledger": dict(ledger), "channel_sample": None if previous_row is None else dict(previous_row),
            "channel_time": previous_time, "steps": int(taken), "stations_reached": list(reached),
        })

    def sample(kind):
        nonlocal previous_row, previous_time, ledger, observation_cpu
        try:
            row = scalar_observation_row(pair, current, mark, mode)
        except Exception as error:
            row = {"time": float(mark), "control_mode": mode, "sample_kind": kind,
                   "observation_error": type(error).__name__ + ": " + str(error),
                   "observation_cpu_seconds": 0.0, "dense_propagator_stored": False,
                   "quadrature": "trapezoid_coordinate_time"}
        row["sample_kind"] = kind
        observation_cpu += float(row.get("observation_cpu_seconds") or 0.0)
        if (previous_row is not None and previous_time is not None
                and _channels_ready(previous_row) and _channels_ready(row)):
            ledger = accumulate_work_ledger(ledger, previous_row, row, float(mark) - float(previous_time))
        previous_row = row
        previous_time = float(mark)
        observations.append(row)
        return row

    def finish(status, stop, admissible=None):
        if admissible is not None:
            return _advance_result(
                admissible["state"], admissible["time"], admissible["steps"], status, used, reached,
                restriction, stop, admissible["clocks"], admissible["clock_rates"], snapshots,
                observations, ledger, step_cpu, observation_cpu, previous_row, previous_time,
            )
        retain = status == "chart_exit" or taken != start_steps or abs(mark - start_time) > 1e-12
        if retain:
            remember("last_admissible" if status == "chart_exit" else "event")
        return _advance_result(
            current, mark, taken, status, used, reached, restriction, stop, clocks, rate_sample,
            snapshots, observations, ledger, step_cpu, observation_cpu, previous_row, previous_time,
        )

    def over_memory():
        return int(rss_bytes()) > int(memory_limit_bytes)

    def forecast_blocks():
        if last_cost is None:
            return False
        return used + float(forecast_factor) * last_cost > float(cpu_allowance) + 1e-12

    while mark + 1e-12 < stations[-1]:
        if until_time is not None and mark >= float(until_time) - 1e-12:
            return finish("until_time", {"stop_class": "until_time", "until_time": float(until_time),
                                         "max_T_forced": False})
        if max_steps is not None and taken >= int(max_steps):
            return finish("step_limit", {"stop_class": "step_limit", "max_steps": int(max_steps),
                                         "max_T_forced": False})
        if over_memory():
            return finish("memory_stop", {
                "stop_class": "memory_stop", "memory_limit_bytes": int(memory_limit_bytes),
                "resident_peak_bytes": int(rss_bytes()), "state_clamped": False,
            })
        if forecast_blocks():
            return finish("budget_stop", {
                "stop_class": "budget_stop", "cpu_allowance": float(cpu_allowance),
                "cpu_spent": used, "forecast_factor": float(forecast_factor),
                "next_step_forecast": None if last_cost is None else float(forecast_factor) * last_cost,
            })
        if previous_row is None:
            sample("start")
        target = next(station for station in stations if station > mark + 1e-12)
        cadence_target = next_observation_time(mark, cadence)
        pair, current, _resolve_info = resolve(pair, current, fft=fft, backend=backend)
        if forced_dt is None:
            restriction_dt, restriction = step_restriction(pair, current, step_cap)
        else:
            restriction_dt = float(forced_dt)
            restriction = {"dt": restriction_dt, "forced_dt": True, "stability_certificate": False,
                           "field_cfl_alone": False}
        gaps = [float(restriction_dt), target - mark, cadence_target - mark]
        if until_time is not None:
            gaps.append(float(until_time) - mark)
        positive = [gap for gap in gaps if gap > 1e-12]
        if not positive:
            mark = max(mark, float(target))
            continue
        dt = min(positive)
        snapshot = current.copy()
        snapshot_clocks = clocks.copy()
        snapshot_rates = rate_sample.copy()
        snapshot_ledger = dict(ledger)
        snapshot_row = None if previous_row is None else dict(previous_row)
        snapshot_row_time = previous_time
        started = time.process_time()
        try:
            updated = stepper(pair, current, dt)
        except coupling.PositiveChartExit as exit_chart:
            admissible = {"state": snapshot, "time": mark, "steps": taken,
                          "clocks": snapshot_clocks, "clock_rates": snapshot_rates}
            remember_state = snapshots
            if not remember_state or abs(remember_state[-1]["time"] - mark) > 1e-12:
                snapshots.append({
                    "kind": "last_admissible", "time": float(mark), "state": snapshot,
                    "clocks": snapshot_clocks, "clock_rates": snapshot_rates,
                    "ledger": snapshot_ledger, "channel_sample": snapshot_row,
                    "channel_time": snapshot_row_time, "steps": int(taken),
                    "stations_reached": list(reached),
                })
            return _advance_result(
                snapshot, mark, taken, "chart_exit", used, reached, restriction,
                _chart_exit_record(exit_chart.reason, mark, dt), snapshot_clocks, snapshot_rates,
                snapshots, observations, snapshot_ledger, step_cpu, observation_cpu,
                snapshot_row, snapshot_row_time,
            )
        measured = float(cpu_per_step) if cpu_per_step is not None else time.process_time() - started
        step_cpu += measured
        current = updated
        new_rates = _clock_rates(pair, current)
        clocks = clocks + 0.5 * dt * (rate_sample + new_rates)
        rate_sample = new_rates
        mark += dt
        taken += 1
        used += measured
        last_cost = measured
        on_station = abs(mark - target) <= 1e-12
        on_cadence = abs(mark - cadence_target) <= 1e-12
        if on_station:
            mark = float(target)
            if float(target) not in reached:
                reached.append(float(target))
        if on_cadence or on_station:
            sample("station" if on_station else "cadence")
        if on_station:
            remember("station")
    return finish("station_reached", {"stop_class": "station_reached", "stations": list(stations)})


def _channels_ready(row):
    return all(row.get(name) is not None for name in WORK_CHANNELS)


def _advance_result(state, mark, taken, status, used, reached, restriction, stop, clocks, clock_rates,
                    snapshots, observations, ledger, step_cpu, observation_cpu, channel_sample,
                    channel_time):
    return {"state": state, "time": float(mark), "steps": int(taken), "status": status,
            "cpu_seconds": float(used), "stations_reached": list(reached),
            "restriction": restriction, "stop": stop or {"stop_class": status},
            "normal_clocks": np.asarray(clocks, dtype=np.float64),
            "clock_rates": np.asarray(clock_rates, dtype=np.float64),
            "snapshots": snapshots, "observations": observations, "work_ledger": ledger,
            "step_cpu_seconds": float(step_cpu), "observation_cpu_seconds": float(observation_cpu),
            "channel_sample": channel_sample, "channel_time": channel_time,
            "stability_certificate": False, "state_clamped": False,
            "physical_instability_claimed": False, "dense_propagator_stored": False}


def _manifest_path(directory):
    return Path(directory) / "manifest.json"


def _ledger_path(directory):
    return Path(directory) / "cpu-ledger.json"


def _lock_path(directory):
    return Path(directory) / "cpu-ledger.lock"


def read_manifest(directory):
    path = _manifest_path(directory)
    if not path.is_file():
        raise FileNotFoundError("campaign manifest is missing: " + str(path))
    return json.loads(path.read_text())


def _write_json(path, record):
    payload = json.dumps(jsonable(record), indent=2, allow_nan=False) + "\n"
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(payload)
    os.replace(temporary, path)


def _commit_bytes(directory, case_id, ordinal, suffix, payload, limit):
    directory = Path(directory)
    final = directory / f"{case_id}-{ordinal:06d}{suffix}"
    if final.exists():
        raise FileExistsError("refusing to overwrite immutable chunk " + final.name)
    temporary = directory / f".{final.name}.{os.getpid()}.tmp"
    temporary.write_bytes(payload)
    size = temporary.stat().st_size
    if size > int(limit):
        temporary.unlink()
        raise RuntimeError("checkpoint chunk exceeds 64 MiB limit" if int(limit) == CHUNK_LIMIT_BYTES
                           else "checkpoint chunk exceeds the byte limit")
    try:
        os.link(temporary, final)
    except FileExistsError:
        temporary.unlink(missing_ok=True)
        raise FileExistsError("refusing to overwrite immutable chunk " + final.name) from None
    temporary.unlink()
    os.chmod(final, 0o444)
    return final, size


def _real_array_names(representation):
    common = ("Q", "r", "chi", "source_weights", "clock_rates", "normal_clocks", "W")
    if representation == CANONICAL_PI:
        return ("pi_Q", "pi_r", "pi_chi") + common
    if representation == NODAL_MOMENTUM:
        return ("nodal_p_Q", "nodal_p_r", "nodal_p_chi", "dx_g") + common
    raise ValueError("momentum representation must be canonical_pi or nodal_momentum")


def commit_checkpoint(directory, case, arrays, *, limit=CHUNK_LIMIT_BYTES):
    """Commit one immutable chunk. Continuation reads these bytes back."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    case_id = str(case["case_id"])
    representation = case.get("momentum_representation", CANONICAL_PI)
    ordinal = _next_ordinal(directory, case_id)
    stored = {}
    for name in _real_array_names(representation):
        stored[name] = np.ascontiguousarray(arrays[name], dtype=np.float64)
    for name in _COMPLEX_ARRAYS:
        stored[name] = np.ascontiguousarray(arrays[name], dtype=np.complex128)
    temporary_npz = directory / f".{case_id}-{ordinal:06d}.{os.getpid()}.npz"
    with temporary_npz.open("wb") as stream:
        np.savez(stream, **stored)
    npz_bytes = temporary_npz.read_bytes()
    temporary_npz.unlink()
    digest = hashlib.sha256(npz_bytes).hexdigest()
    record = dict(case)
    record.update({
        "schema": SCHEMA,
        "ordinal": ordinal,
        "arrays_sha256": digest,
        "momentum_representation": case.get("momentum_representation", CANONICAL_PI),
        "immutable": True,
        "stability_certificate": False,
        "array_sha256": {name: array_sha256(stored[name]) for name in stored},
    })
    json_bytes = (json.dumps(jsonable(record), indent=2, allow_nan=False) + "\n").encode()
    if len(npz_bytes) + len(json_bytes) > int(limit):
        raise RuntimeError("checkpoint chunk exceeds the byte limit")
    npz_path, npz_size = _commit_bytes(directory, case_id, ordinal, ".npz", npz_bytes, limit)
    json_path, json_size = _commit_bytes(directory, case_id, ordinal, ".json", json_bytes, limit)
    record["npz"] = npz_path.name
    record["json"] = json_path.name
    record["npz_bytes"] = npz_size
    record["json_bytes"] = json_size
    return record


def _next_ordinal(directory, case_id):
    highest = -1
    for path in Path(directory).glob(case_id + "-*.json"):
        stem = path.name[: -len(".json")]
        try:
            highest = max(highest, int(stem.rsplit("-", 1)[1]))
        except ValueError:
            continue
    return highest + 1


def load_checkpoint(directory, case_id, ordinal=None):
    directory = Path(directory)
    if ordinal is None:
        ordinal = _next_ordinal(directory, case_id) - 1
    if ordinal < 0:
        raise FileNotFoundError("no checkpoint for " + case_id)
    json_path = directory / f"{case_id}-{ordinal:06d}.json"
    npz_path = directory / f"{case_id}-{ordinal:06d}.npz"
    record = json.loads(json_path.read_text())
    digest = file_sha256(npz_path)
    if digest != record["arrays_sha256"]:
        raise ValueError("checkpoint bytes do not match the commit record")
    with np.load(npz_path, allow_pickle=False) as stored:
        arrays = {name: np.array(stored[name], copy=True) for name in stored.files}
    for name, digest_value in record["array_sha256"].items():
        if array_sha256(arrays[name]) != digest_value:
            raise ValueError("checkpoint array changed: " + name)
    return record, arrays


def state_from_arrays(arrays, representation=CANONICAL_PI):
    if representation == CANONICAL_PI:
        momenta = (arrays["pi_Q"], arrays["pi_r"], arrays["pi_chi"])
    elif representation == NODAL_MOMENTUM:
        spacing = float(np.asarray(arrays["dx_g"]).reshape(-1)[0])
        basis = arrays["W"]
        momenta = tuple(spacing * (basis.T @ arrays[name]) for name in ("nodal_p_Q", "nodal_p_r", "nodal_p_chi"))
    else:
        raise ValueError("momentum representation must be canonical_pi or nodal_momentum")
    return model.NestedState(arrays["Q"], arrays["r"], arrays["chi"], *momenta, arrays["phi0"], arrays["phi1"])


def pair_from_arrays(arrays, record):
    """Rebuild the pair from stored W and columns. No SVD and no initial solve."""
    nf = int(record["nf"])
    weights = np.array(arrays["source_weights"], dtype=float, copy=True)
    grid = galerkin.build_grid(nf, gauge="conformal")
    grid.fine = replace(grid.fine, occupations=weights.copy())
    if grid.ng != int(arrays["W"].shape[0]):
        raise ValueError("stored W does not match the conformal geometry degree")
    observer = arrays["observer_columns"]
    source = np.vstack((arrays["source_phi0"], arrays["source_phi1"]))
    return model.NestedPair(
        grid, _readonly(arrays["W"]),
        _readonly(record["coarse_indices"], int), _readonly(record["child_indices"], int),
        _readonly(record["parent_indices"], int),
        _readonly(observer[:nf]), _readonly(observer[nf:]), _readonly(observer),
        _readonly(arrays["source_phi0"]), _readonly(arrays["source_phi1"]), _readonly(source),
        _readonly(weights), dict(record.get("source_metadata") or {}),
        dict(record.get("geometry_metadata") or {}),
        clock_locations=tuple(record.get("clock_locations", (1.0, 2.0, 3.0))),
    )


def arrays_from_state(state, *, basis_arrays, clocks, representation=CANONICAL_PI, nodal=None):
    if representation == CANONICAL_PI:
        payload = {"Q": state.Q, "r": state.r, "chi": state.chi, "pi_Q": state.p_Q,
                   "pi_r": state.p_r, "pi_chi": state.p_chi}
    elif representation == NODAL_MOMENTUM:
        if nodal is None:
            raise ValueError("nodal momentum representation needs the nodal momenta")
        payload = {"Q": state.Q, "r": state.r, "chi": state.chi,
                   "nodal_p_Q": nodal[0], "nodal_p_r": nodal[1], "nodal_p_chi": nodal[2]}
    else:
        raise ValueError("unknown momentum representation")
    payload.update(phi0=state.phi0, phi1=state.phi1, W=basis_arrays["W"],
                   source_phi0=basis_arrays["source_phi0"], source_phi1=basis_arrays["source_phi1"],
                   observer_columns=basis_arrays["observer_columns"],
                   source_weights=basis_arrays["source_weights"],
                   clock_rates=np.asarray(clocks["rates"], dtype=float),
                   normal_clocks=np.asarray(clocks["normal_clocks"], dtype=float))
    return payload


def _saved_records():
    v1 = json.loads(V1_JSON.read_text())
    basis = json.loads(BASIS_JSON.read_text())
    if file_sha256(V1_NPZ) != v1["payload_sha256"]:
        raise ValueError("saved separated-pair payload does not match its record")
    if file_sha256(BASIS_NPZ) != basis["payload_sha256"]:
        raise ValueError("frozen replay basis payload does not match its record")
    return v1, basis


def load_saved_handoff(case_name):
    """Exact T=0.3 baseline state and frozen W. No initial solve and no new W."""
    if case_name not in BASELINE_CASES:
        raise ValueError("handoff cases are the saved nf128/nf256 baselines")
    v1, basis = _saved_records()
    case = v1["results"][case_name]
    if abs(float(case["summary"]["final"]["time"]) - HANDOFF_TIME) > 1e-12:
        raise ValueError("saved baseline does not end at T=0.3")
    step_cap = _canonical_cap(case["step_cap"])
    nf = int(case["nf"])
    entry = basis["bases"][str(nf)]
    with np.load(V1_NPZ, allow_pickle=False) as episode, np.load(BASIS_NPZ, allow_pickle=False) as frozen:
        fields = {name: np.array(episode[case_name + "_" + name][-1], copy=True) for name in model.STATE_NAMES}
        times = np.array(episode[case_name + "_times"], copy=True)
        observer = np.array(episode[case_name + "_reference_columns"], copy=True)
        basis_w = np.array(frozen[f"nf{nf}_W"], copy=True)
        basis_observer = np.array(frozen[f"nf{nf}_reference_columns"], copy=True)
        basis_source = np.array(frozen[f"nf{nf}_source_columns"], copy=True)
    if abs(float(times[-1]) - HANDOFF_TIME) > 1e-12:
        raise ValueError("saved baseline time array does not end at T=0.3")
    if not np.array_equal(observer, basis_observer):
        raise ValueError("saved observer and frozen basis observer differ")
    state = model.NestedState(*(fields[name] for name in model.STATE_NAMES))
    if state_sha256(state) != case["summary"]["final_state_sha256"]:
        raise ValueError("saved baseline final state does not match its pin")
    weights = np.array(case["source"]["occupations"], dtype=float)
    clocks = {
        "locations": list(case["summary"]["final"]["metrics"]["clock_locations"]),
        "rates": list(case["summary"]["final"]["metrics"]["clock_rates"]),
        "normal_clocks": list(case["summary"]["normal_clocks"][-1]),
        "protocol": case["summary"]["clock_protocol"],
    }
    pins = {
        "v1_json_sha256": file_sha256(V1_JSON),
        "v1_npz_sha256": v1["payload_sha256"],
        "basis_json_sha256": file_sha256(BASIS_JSON),
        "basis_npz_sha256": basis["payload_sha256"],
        "final_state_sha256": case["summary"]["final_state_sha256"],
        "source_columns_sha256": array_sha256(basis_source),
        "observer_columns_sha256": array_sha256(basis_observer),
        "W_sha256": array_sha256(basis_w),
        "weights_sha256": array_sha256(weights),
        "phi_sha256": array_sha256(np.vstack((state.phi0, state.phi1))),
    }
    return {"case_name": case_name, "nf": nf, "step_cap": step_cap,
            "state": state, "W": basis_w, "observer_columns": basis_observer,
            "source_columns": basis_source, "weights": weights, "clocks": clocks,
            "coarse_indices": entry["coarse_indices"], "child_indices": entry["child_indices"],
            "parent_indices": entry["parent_indices"], "source_metadata": entry["source"],
            "geometry_metadata": entry["geometry"], "pins": pins,
            "initial_state_called": False, "build_pair_called": False,
            "momentum_representation": CANONICAL_PI}


def _case_record(handoff, *, case_id, control_mode, step_cap, stations):
    mode = normalize_control_mode(control_mode)
    pins = dict(handoff["pins"])
    return {
        "case_id": case_id, "parent_case": handoff["case_name"], "nf": handoff["nf"],
        "control": mode, "control_mode": mode,
        "geometry": "frozen" if mode == "frozen_geometry" else "evolving",
        "step_cap": float(step_cap), "stations": list(stations),
        "coordinate_time": HANDOFF_TIME, "steps": 0, "stations_reached": [],
        "handoff_time": HANDOFF_TIME, "phi": "saved exact Phi(0.3) from the dt=0.0005 baseline",
        "observer": "saved fixed reference columns",
        "momentum_representation": CANONICAL_PI, "source_pins": pins,
        "source_pins_before": pins, "source_pins_after": pins,
        "coarse_indices": handoff["coarse_indices"], "child_indices": handoff["child_indices"],
        "parent_indices": handoff["parent_indices"], "source_metadata": handoff["source_metadata"],
        "geometry_metadata": handoff["geometry_metadata"],
        "clock_locations": list(handoff["clocks"]["locations"]),
        "clock_protocol": handoff["clocks"]["protocol"],
        "verify_external_pins": True, "evolution_backend_requested": "auto",
        "initializer": "dense_saved_state",
        "reduction_backend": REDUCTION_BACKEND, "reduction_backend_used": False,
        "stability_certificate": False, "initial_state_called": False,
        "stage": 0, "snapshot_kind": "handoff", "work_ledger": empty_work_ledger(mode),
        "dense_propagator_stored": False,
    }


def _arrays_for_handoff(handoff):
    source = handoff["source_columns"]
    nf = handoff["nf"]
    return arrays_from_state(handoff["state"], representation=CANONICAL_PI, basis_arrays={
        "W": handoff["W"], "source_phi0": source[:nf], "source_phi1": source[nf:],
        "observer_columns": handoff["observer_columns"], "source_weights": handoff["weights"],
    }, clocks={"rates": handoff["clocks"]["rates"], "normal_clocks": handoff["clocks"]["normal_clocks"]})


def prepare(output, *, stations=None, include_frozen=True):
    """Stage 0. Write exact handoff chunks. Do not evolve."""
    directory = assert_campaign_output(output)
    if _manifest_path(directory).exists():
        raise FileExistsError("campaign output already exists; refusing to overwrite")
    directory.mkdir(parents=True, exist_ok=True)
    selected = normalize_stations(stations)
    cases = []
    handoffs = {}
    for spec in first_batch_specs(include_frozen=include_frozen):
        if spec["parent_case"] not in handoffs:
            handoffs[spec["parent_case"]] = load_saved_handoff(spec["parent_case"])
        handoff = handoffs[spec["parent_case"]]
        record = _case_record(
            handoff, case_id=spec["case_id"], control_mode=spec["control_mode"],
            step_cap=spec["step_cap"], stations=selected,
        )
        committed = commit_checkpoint(directory, record, _arrays_for_handoff(handoff))
        cases.append({key: committed[key] for key in (
            "case_id", "ordinal", "npz", "json", "arrays_sha256", "parent_case", "nf",
            "control", "control_mode", "geometry", "step_cap", "coordinate_time", "steps")})
    manifest = {"schema": SCHEMA, "stage": 0, "status": "PREPARED", "stations": list(selected),
                "cases": cases, "workers": DEFAULT_WORKERS,
                "cpu_budget_seconds": DEFAULT_CPU_BUDGET_SECONDS,
                "cpu_accounting": "sum_of_child_process_time",
                "forecast_factor": DEFAULT_FORECAST_FACTOR,
                "memory_limit_bytes": DEFAULT_MEMORY_BYTES,
                "chunk_limit_bytes": CHUNK_LIMIT_BYTES,
                "evolution_backend": EVOLUTION_BACKEND, **reduction_backend_status(),
                "initial_state_called": False, "evolved": False, "stability_certificate": False,
                "zero_residual_gate": False, "max_T_forced": False,
                "child_cpu_seconds": 0.0}
    _write_json(_manifest_path(directory), manifest)
    return manifest


def _verify_external_pins(record):
    if not record.get("verify_external_pins"):
        return
    pins = record["source_pins"]
    if file_sha256(V1_JSON) != pins["v1_json_sha256"] or file_sha256(V1_NPZ) != pins["v1_npz_sha256"]:
        raise ValueError("saved separated-pair pin changed")
    if file_sha256(BASIS_JSON) != pins["basis_json_sha256"] or file_sha256(BASIS_NPZ) != pins["basis_npz_sha256"]:
        raise ValueError("frozen basis pin changed")


def _live_pins(arrays, state):
    source = np.vstack((np.array(arrays["source_phi0"], copy=True), np.array(arrays["source_phi1"], copy=True)))
    phi = np.vstack((np.array(state.phi0, copy=True), np.array(state.phi1, copy=True)))
    clocks = np.array(arrays["clock_rates"], dtype=float, copy=True)
    weights = np.array(arrays["source_weights"], dtype=float, copy=True)
    return {
        "source_columns_sha256": array_sha256(source),
        "observer_columns_sha256": array_sha256(arrays["observer_columns"]),
        "W_sha256": array_sha256(arrays["W"]),
        "weights_sha256": array_sha256(weights),
        "phi_sha256": array_sha256(phi),
        "state_sha256": state_sha256(state),
        "clock_rates_sha256": array_sha256(clocks),
    }


def _pins_unchanged(before, after):
    if not before or not after:
        return
    for key in ("W_sha256", "source_columns_sha256", "observer_columns_sha256", "weights_sha256"):
        if key in before and key in after and before[key] != after[key]:
            raise ValueError("restart changed a frozen source pin: " + key)
    if after.get("phi_sha256") == after.get("source_columns_sha256"):
        raise ValueError("state Phi is aliased to the source columns")
    if after.get("clock_rates_sha256") == after.get("weights_sha256"):
        raise ValueError("clock rates are aliased to the source weights")


def execute_case(directory, case_id, *, cpu_allowance, forecast_factor=DEFAULT_FORECAST_FACTOR,
                 memory_limit_bytes=DEFAULT_MEMORY_BYTES, max_steps=None, fft=None, backend=None,
                 stepper=None, cpu_per_step=None, rss_bytes=None, diagnostics=False,
                 observation_cadence=OBSERVATION_CADENCE, forced_dt=None, until_time=None):
    """Continue one committed case. The handoff solve is not repeated."""
    record, arrays = load_checkpoint(directory, case_id)
    _verify_external_pins(record)
    _pins_unchanged(record.get("source_pins_before") or record.get("source_pins"),
                    record.get("source_pins_after") or record.get("source_pins"))
    if record.get("initial_state_called"):
        raise RuntimeError("checkpoint claims a new initial solve")
    dense_pair = pair_from_arrays(arrays, record)
    representation = record.get("momentum_representation", CANONICAL_PI)
    if representation != CANONICAL_PI:
        raise ValueError("continuation of this batch stores canonical pi")
    state = state_from_arrays(arrays, CANONICAL_PI)
    pair, state, resolve_info = resolve_pair(dense_pair, state, fft=fft, backend=backend)
    mode = normalize_control_mode(record.get("control_mode") or record.get("geometry"))
    started = time.process_time()
    result = advance_case(
        pair, state, coordinate_time=record["coordinate_time"], steps=record["steps"],
        stations=record["stations"], step_cap=record["step_cap"],
        geometry_frozen=mode == "frozen_geometry", cpu_allowance=cpu_allowance,
        spent=0.0, forecast_factor=forecast_factor, memory_limit_bytes=memory_limit_bytes,
        max_steps=max_steps, stepper=stepper, cpu_per_step=cpu_per_step, rss_bytes=rss_bytes,
        fft=fft, backend=backend, normal_clocks=np.array(arrays["normal_clocks"], dtype=float, copy=True),
        control_mode=mode, observation_cadence=observation_cadence, forced_dt=forced_dt,
        until_time=until_time, work_ledger=record.get("work_ledger"),
        channel_sample=record.get("channel_sample"), channel_time=record.get("channel_time"),
    )
    child_cpu = float(cpu_per_step) * max(result["steps"] - int(record["steps"]), 0) if cpu_per_step is not None else time.process_time() - started
    chunk = None
    append_observations(directory, case_id, result["observations"])
    pins_before = record.get("source_pins_after") or record.get("source_pins")
    for snapshot in result["snapshots"]:
        snap_changed = any(not np.array_equal(getattr(snapshot["state"], name), getattr(state, name))
                           for name in model.STATE_NAMES) or snapshot["time"] != record["coordinate_time"]
        if not snap_changed and result["status"] != "chart_exit":
            continue
        new_record = dict(record)
        new_record.update(coordinate_time=snapshot["time"], steps=snapshot["steps"],
                          stations_reached=snapshot["stations_reached"], stage=2,
                          stop=result["stop"] if abs(snapshot["time"] - result["time"]) <= 1e-12 else None,
                          status=result["status"] if abs(snapshot["time"] - result["time"]) <= 1e-12 else "station_retained",
                          evolved=snap_changed, snapshot_kind=snapshot["kind"],
                          work_ledger=snapshot["ledger"], channel_sample=snapshot["channel_sample"],
                          channel_time=snapshot.get("channel_time"),
                          step_cpu_seconds=result["step_cpu_seconds"],
                          observation_cpu_seconds=result["observation_cpu_seconds"],
                          resolve_pair={key: resolve_info[key] for key in (
                              "resolver", "backend_requested", "backend", "fft_requested",
                              "fft_applied", "state_changed", "W_changed", "initializer")},
                          control_mode=mode, stability_certificate=False, dense_propagator_stored=False)
        new_arrays = dict(arrays)
        new_arrays.update(Q=np.array(snapshot["state"].Q, copy=True),
                          r=np.array(snapshot["state"].r, copy=True),
                          chi=np.array(snapshot["state"].chi, copy=True),
                          pi_Q=np.array(snapshot["state"].p_Q, copy=True),
                          pi_r=np.array(snapshot["state"].p_r, copy=True),
                          pi_chi=np.array(snapshot["state"].p_chi, copy=True),
                          phi0=np.array(snapshot["state"].phi0, copy=True),
                          phi1=np.array(snapshot["state"].phi1, copy=True),
                          clock_rates=np.array(snapshot["clock_rates"], dtype=float, copy=True),
                          normal_clocks=np.array(snapshot["clocks"], dtype=float, copy=True))
        new_record["source_pins_before"] = pins_before
        new_record["source_pins_after"] = _live_pins(new_arrays, snapshot["state"])
        _pins_unchanged(new_record["source_pins_before"], new_record["source_pins_after"])
        if diagnostics and abs(snapshot["time"] - result["time"]) <= 1e-12:
            new_record["stability"] = observe(
                pair, snapshot["state"], snapshot["time"], control_mode=mode,
            )["stability"]
        committed = commit_checkpoint(directory, new_record, new_arrays)
        pins_before = new_record["source_pins_after"]
        chunk = committed["npz"]
        record = new_record
        arrays = new_arrays
        state = snapshot["state"]
    return {"case_id": case_id, "status": result["status"], "child_cpu_seconds": child_cpu,
            "coordinate_time": result["time"], "steps": result["steps"],
            "stations_reached": result["stations_reached"], "stop": result["stop"],
            "chunk": chunk, "stability_certificate": False, "state_clamped": False,
            "physical_instability_claimed": bool(result["stop"].get("physical_instability_claimed", False))}


def _ledger_update(directory, *, pilot=0.0, reserve=0.0, release=0.0, actual=0.0, budget):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    lock_path = _lock_path(directory)
    ledger_path = _ledger_path(directory)
    with lock_path.open("a+") as lock_stream:
        fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX)
        try:
            if ledger_path.exists():
                ledger = json.loads(ledger_path.read_text())
            else:
                ledger = {"spent": 0.0, "reserved": 0.0}
            spent = float(ledger["spent"]) + float(pilot) + float(actual)
            reserved = float(ledger["reserved"]) + float(reserve) - float(release)
            if reserved < -1e-9:
                raise RuntimeError("CPU reservation became negative")
            admitted = True
            if reserve and spent + reserved > float(budget) + 1e-9:
                reserved -= float(reserve)
                admitted = False
            ledger = {"spent": spent, "reserved": max(0.0, reserved), "budget": float(budget),
                      "accounting": "sum_of_child_process_time"}
            _write_json(ledger_path, ledger)
            return admitted, ledger
        finally:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_UN)


def estimate_case_cpu(pair, state, *, step_cap, stations, current_time, probe_step_cpu, probe_diag_cpu):
    dt, _restriction = step_restriction(pair, state, step_cap)
    span = float(normalize_stations(stations)[-1] - float(current_time))
    steps = int(np.ceil(max(span, 0.0) / dt - 1e-12))
    diagnostics = len(normalize_stations(stations)) + 1
    return {"dt": dt, "estimated_steps": steps, "estimated_cpu_seconds": steps * float(probe_step_cpu)
            + diagnostics * float(probe_diag_cpu)}


def _case_worker(payload):
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[name] = "1"
    directory = payload["directory"]
    case_id = payload["case_id"]
    budget = float(payload["cpu_budget_seconds"])
    factor = float(payload["forecast_factor"])
    record, arrays = load_checkpoint(directory, case_id)
    if record.get("status") in ("station_reached", "chart_exit"):
        return {"case_id": case_id, "status": record["status"], "child_cpu_seconds": 0.0,
                "coordinate_time": record["coordinate_time"], "steps": record["steps"],
                "stations_reached": record.get("stations_reached", []), "stop": record.get("stop", {}),
                "chunk": None, "stability_certificate": False, "state_clamped": False,
                "physical_instability_claimed": False}
    dense_pair = pair_from_arrays(arrays, record)
    state = state_from_arrays(arrays, CANONICAL_PI)
    pair, state, _resolved = resolve_pair(
        dense_pair, state, backend=payload.get("backend", "auto"), fft=payload.get("fft"),
    )
    mode = normalize_control_mode(record.get("control_mode") or record.get("geometry"))
    probe_started = time.process_time()
    copy = state.copy()
    stepper = frozen_geometry_step if mode == "frozen_geometry" else evolving_step
    dt, _restriction = step_restriction(pair, copy, record["step_cap"])
    probe_dt = min(dt, 1e-4)
    try:
        stepper(pair, copy, probe_dt)
    except coupling.PositiveChartExit as exit_chart:
        probe_cpu = time.process_time() - probe_started
        _ledger_update(directory, pilot=probe_cpu, budget=budget)
        stop = _chart_exit_record(exit_chart.reason, record["coordinate_time"], probe_dt)
        chart_record = dict(record)
        chart_record.update(stage=2, status="chart_exit", stop=stop, evolved=False,
                            stability_certificate=False, state_clamped=False,
                            physical_instability_claimed=False)
        committed = commit_checkpoint(directory, chart_record, arrays)
        return {"case_id": case_id, "status": "chart_exit", "child_cpu_seconds": probe_cpu,
                "coordinate_time": record["coordinate_time"], "steps": record["steps"],
                "stations_reached": record.get("stations_reached", []),
                "stop": stop, "chunk": committed["npz"], "stability_certificate": False,
                "state_clamped": False, "physical_instability_claimed": False}
    probe_step = time.process_time() - probe_started
    diag_started = time.process_time()
    observe(pair, state, record["coordinate_time"], control_mode=mode)
    probe_diag = time.process_time() - diag_started
    probe_cpu = time.process_time() - probe_started
    estimate = estimate_case_cpu(pair, state, step_cap=record["step_cap"], stations=record["stations"],
                                 current_time=record["coordinate_time"],
                                 probe_step_cpu=max(probe_step, 1e-6), probe_diag_cpu=probe_diag)
    demand = factor * estimate["estimated_cpu_seconds"]
    admitted, _ledger = _ledger_update(directory, pilot=probe_cpu, reserve=demand, budget=budget)
    if not admitted:
        return {"case_id": case_id, "status": "budget_stop", "child_cpu_seconds": probe_cpu,
                "coordinate_time": record["coordinate_time"], "steps": record["steps"],
                "stations_reached": record.get("stations_reached", []),
                "stop": {"stop_class": "budget_stop", "forecast_cpu_seconds": demand,
                         "forecast_factor": factor, "probe_only": True},
                "chunk": None, "stability_certificate": False, "state_clamped": False,
                "physical_instability_claimed": False}
    result = execute_case(directory, case_id, cpu_allowance=demand, forecast_factor=factor,
                          memory_limit_bytes=payload["memory_limit_bytes"],
                          max_steps=payload.get("max_steps"), fft=payload.get("fft"),
                          backend=payload.get("backend", "auto"), diagnostics=True)
    _ledger_update(directory, release=demand, actual=result["child_cpu_seconds"], budget=budget)
    result["child_cpu_seconds"] = float(result["child_cpu_seconds"]) + probe_cpu
    result["forecast_cpu_seconds"] = demand
    return result


def run(output, *, workers=DEFAULT_WORKERS, cpu_budget_seconds=DEFAULT_CPU_BUDGET_SECONDS,
        forecast_factor=DEFAULT_FORECAST_FACTOR, memory_limit_bytes=DEFAULT_MEMORY_BYTES,
        max_steps=None, fft=None, backend=None, executor=None):
    """Stage 2. One process pool. Stop and keep the last admissible checkpoint."""
    directory = assert_campaign_output(output)
    if fft is True:
        requested_backend = "fft"
    elif fft is False:
        requested_backend = "dense"
    elif backend is None:
        requested_backend = "auto"
    else:
        requested_backend = str(backend)
    if int(workers) < 1:
        raise ValueError("workers must be positive")
    if float(forecast_factor) < 1:
        raise ValueError("forecast factor must be at least 1")
    if float(cpu_budget_seconds) <= 0 or int(memory_limit_bytes) <= 0:
        raise ValueError("CPU budget and memory limit must be positive")
    manifest = read_manifest(directory)
    if manifest.get("evolved") and max_steps is None and manifest.get("status") == "STATION_REACHED":
        return manifest
    payloads = [{"directory": str(directory), "case_id": case["case_id"],
                 "cpu_budget_seconds": float(cpu_budget_seconds),
                 "forecast_factor": float(forecast_factor),
                 "memory_limit_bytes": int(memory_limit_bytes),
                 "max_steps": max_steps, "fft": None, "backend": requested_backend}
                for case in manifest["cases"]]
    results = []
    pool_factory = ProcessPoolExecutor if executor is None else executor
    created = 0

    def factory(max_workers):
        nonlocal created
        created += 1
        return pool_factory(max_workers=max_workers)

    with factory(int(workers)) as pool:
        if created != 1:
            raise RuntimeError("the campaign run opened more than one executor pool")
        futures = [pool.submit(_case_worker, payload) for payload in payloads]
        pending = set(futures)
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for future in done:
                results.append(future.result())
    child_cpu = float(sum(item["child_cpu_seconds"] for item in results))
    manifest = read_manifest(directory)
    manifest.update(stage=2, status=_campaign_status(results), results=results,
                    child_cpu_seconds=child_cpu, workers=int(workers),
                    cpu_budget_seconds=float(cpu_budget_seconds),
                    forecast_factor=float(forecast_factor),
                    memory_limit_bytes=int(memory_limit_bytes),
                    cpu_accounting="sum_of_child_process_time", evolved=True,
                    executor_pools=1, stability_certificate=False, max_T_forced=False,
                    zero_residual_gate=False)
    _write_json(_manifest_path(directory), manifest)
    return manifest


def _campaign_status(results):
    statuses = {item["status"] for item in results}
    if statuses == {"station_reached"}:
        return "STATION_REACHED"
    if "chart_exit" in statuses:
        return "CHART_EXIT"
    if "memory_stop" in statuses:
        return "MEMORY_STOP"
    if "budget_stop" in statuses:
        return "BUDGET_STOP"
    if "step_limit" in statuses:
        return "STEP_LIMIT"
    return "INCOMPLETE"


def observation_path(directory, case_id):
    return Path(directory) / f"{case_id}-observations.jsonl"


def append_observations(directory, case_id, rows):
    """Append compact rows. This stream is not a state chunk and holds no propagator."""
    if not rows:
        return None
    path = observation_path(directory, case_id)
    with path.open("a", encoding="utf-8") as stream:
        for row in rows:
            stream.write(json.dumps(jsonable(row), allow_nan=False, separators=(",", ":")) + "\n")
    return path


def read_observations(directory, case_id):
    """Scalar physical rows in write order. Empty when no sample has been stored."""
    path = observation_path(directory, case_id)
    if not path.is_file():
        return []
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.strip():
            rows.append(json.loads(line))
    return rows


def read_station_states(directory, case_id):
    """Full-state commits: handoff, stations, and last-admissible events."""
    directory = Path(directory)
    retained = []
    prefix = str(case_id) + "-"
    for path in sorted(directory.glob(prefix + "*.json")):
        record = json.loads(path.read_text(encoding="utf-8"))
        kind = record.get("snapshot_kind")
        if kind not in ("handoff", "station", "event", "last_admissible"):
            continue
        retained.append({
            "ordinal": record.get("ordinal"),
            "snapshot_kind": kind,
            "coordinate_time": record.get("coordinate_time"),
            "steps": record.get("steps"),
            "stations_reached": record.get("stations_reached", []),
            "work_ledger": record.get("work_ledger"),
            "npz": record.get("npz"),
            "dense_propagator_stored": False,
        })
    return retained


def read_case_ledger(directory, case_id):
    """Accumulated continuation integrals on the latest retained state."""
    record, arrays = load_checkpoint(directory, case_id)
    return {
        "work_ledger": record.get("work_ledger"),
        "normal_clocks": np.array(arrays["normal_clocks"], dtype=float).tolist(),
        "clock_rates": np.array(arrays["clock_rates"], dtype=float).tolist(),
        "coordinate_time": record.get("coordinate_time"),
        "clock_quadrature": "trapezoid_coordinate_time_per_accepted_step",
        "work_quadrature": "trapezoid_coordinate_time_between_control_rate_samples",
        "dense_propagator_stored": False,
    }


def read_physical_audit(directory, case_id):
    """Read API for the retained physical record of one case."""
    rows = read_observations(directory, case_id)
    stations = read_station_states(directory, case_id)
    ledger = read_case_ledger(directory, case_id)
    record, _arrays = load_checkpoint(directory, case_id)
    observation_cpu = float(record.get("observation_cpu_seconds") or 0.0)
    step_cpu = float(record.get("step_cpu_seconds") or 0.0)
    required = (
        "localization_peak_x", "localization_width", "child_proper_length", "parent_proper_length",
        "child_normal_energy", "parent_normal_energy", "coordinate_fieldwork", "pressure_work",
        "lapse_work", "boundary_child", "boundary_parent", "R_h_max", "aux_discrepancy_max",
        "car_gram_gap", "energy_total",
    )
    present = [name for name in required if rows and all(name in row for row in rows)]
    return {
        "case_id": case_id,
        "observation_times": [row["time"] for row in rows],
        "station_times": [item["coordinate_time"] for item in stations],
        "snapshot_kinds": [item["snapshot_kind"] for item in stations],
        "required_evidence_present": present,
        "work_ledger": ledger["work_ledger"],
        "normal_clocks": ledger["normal_clocks"],
        "clock_quadrature": ledger["clock_quadrature"],
        "work_quadrature": ledger["work_quadrature"],
        "dense_propagator_stored": False,
        "overhead": {
            "observation_cpu_seconds": observation_cpu,
            "step_cpu_seconds": step_cpu,
            "observation_over_step": None if step_cpu <= 0.0 else observation_cpu / step_cpu,
        },
    }


def check(output):
    """Validate commits without evolving, writing, or repairing them."""
    directory = assert_campaign_output(output)
    manifest = read_manifest(directory)
    if manifest.get("schema") != SCHEMA:
        raise ValueError("unexpected campaign schema")
    problems = []
    checked = []
    for case in manifest["cases"]:
        ordinals = []
        prefix = case["case_id"] + "-"
        for path in sorted(directory.glob(prefix + "*.json")):
            ordinals.append(int(path.name[len(prefix):-len(".json")]))
        if not ordinals:
            problems.append(case["case_id"] + " has no commit")
            continue
        for ordinal in ordinals:
            json_path = directory / f"{case['case_id']}-{ordinal:06d}.json"
            npz_path = directory / f"{case['case_id']}-{ordinal:06d}.npz"
            record = json.loads(json_path.read_text())
            for path in (json_path, npz_path):
                if not path.is_file():
                    problems.append(path.name + " missing")
                    continue
                mode = path.stat().st_mode & 0o222
                size = path.stat().st_size
                if mode:
                    problems.append(path.name + " is writable")
                if size > CHUNK_LIMIT_BYTES:
                    problems.append(path.name + " exceeds 64 MiB")
            if file_sha256(npz_path) != record.get("arrays_sha256"):
                problems.append(npz_path.name + " hash mismatch")
            with np.load(npz_path, allow_pickle=False) as stored:
                for name, digest in record.get("array_sha256", {}).items():
                    if array_sha256(stored[name]) != digest:
                        problems.append(npz_path.name + ":" + name + " mutated")
            if record.get("momentum_representation") not in (CANONICAL_PI, NODAL_MOMENTUM):
                problems.append(case["case_id"] + " lacks an explicit momentum representation")
            if record.get("stability_certificate") is not False:
                problems.append(case["case_id"] + " claims a stability certificate")
            if record.get("verify_external_pins"):
                try:
                    _verify_external_pins(record)
                    _pins_unchanged(record.get("source_pins_before") or record.get("source_pins"),
                                    record.get("source_pins_after") or record.get("source_pins"))
                except ValueError as error:
                    problems.append(str(error))
            reloaded, reloaded_arrays = load_checkpoint(directory, case["case_id"], ordinal)
            if reloaded["arrays_sha256"] != record["arrays_sha256"]:
                problems.append(case["case_id"] + " reload mismatch")
            del reloaded_arrays
            checked.append(case["case_id"] + f":{ordinal}")
    report = {"schema": SCHEMA, "ok": not problems, "checked": checked, "problems": problems,
              "evolved": False, "bytes_written": 0, "stability_certificate": False}
    if problems:
        raise ValueError("checkpoint check failed: " + "; ".join(problems))
    return report
