"""Stage-1 observer for one nested conformal pair.

``observe(pair, state, time)`` is the compact-row API for campaign worker 4.
Call it on accepted states. Do not call it on Runge--Kutta stages.
``capture_profiles=False`` keeps the row scalar. Profiles are a separate flag.

Spatial windows are the fixed coordinate intervals child ``(1, 3)`` and
parent ``(0, 4)``. They are not geometry-mode projectors and not the
six-column source projectors. In the conformal gauge the null rays have
coordinate speed ``dx/dt = ±1``, so the coordinate crossing times are the
carrier length 8, the parent width 4 and the child width 2. Proper clocks
stay a separate integral. Packet width is not a dispersal time.

A moving coordinate window adds the Leibniz--Reynolds term of an action
density. Contact work is the same observer traction ``-p_r V`` on that
cut. Neither an imposed window velocity nor this row is regeneration
evidence. Live curvature uses
``derive_nsc_spherical_conformal_curvature.actual_q_second_rate`` and
``direct_rh_grid`` after ``reconstruct_state`` through the frozen ``W``.
``chi + 2`` is compared afterward and is not substituted.
"""
from __future__ import annotations

import math

import numpy as np

from . import nsc_nested_parent_child as nested
from . import nsc_regional_energy_exchange as regional
from . import nsc_spherical_episode_assessment as metric
from . import nsc_spherical_galerkin_coupling as galerkin
from .nsc_spherical_coupling import CauchyRate


API_VERSION = "nsc-discovery-observe-v1"
OBSERVATION_ON_RK_STAGE = False
CARRIER_LENGTH = 8.0
PARENT_INTERVAL = (0.0, 4.0)
CHILD_INTERVAL = (1.0, 3.0)
BOUNDARY_ORDER = (0.0, 1.0, 3.0, 4.0)
CHARACTERISTIC_SPEEDS = (-1.0, 1.0)
CROSSING_TIMES = {"carrier": 8.0, "parent": 4.0, "child": 2.0}
PREPARATION_NAMES = (
    "original_conformal_v2",
    "separated_pair_v1",
    "separated_pair_v2",
    "frozen_basis",
    "new_pair",
)

REQUIRED_TOP_KEYS = (
    "api_version",
    "time",
    "capture_profiles",
    "observation_on_rk_stage",
    "gauge",
    "spatial_windows_are_mode_projectors",
    "source_probability_total",
    "source_probability_interpolant",
    "source_probability_per_proper_length_peak",
    "source_probability_per_proper_length_peak_x",
    "normal_energy_total",
    "coordinate_metric_work_total",
    "flux_projection_shift_total",
    "flux_projection_proper_total",
    "flux_projection_shift_pressure_total",
    "nyquist_symbol_max",
    "child_window",
    "parent_window",
    "boundary_flux",
    "curvature",
    "proper",
    "localization",
    "proper_gradients",
    "characteristics",
    "mode_projectors",
    "preparation",
    "control_mode",
    "frozen_geometry_control",
    "normal_balance_slope_evaluated",
    "chi_proxy_substituted",
    "renewal_asserted",
    "continuum_certified",
    "imposed_boundary_trajectory_is_regeneration_evidence",
)
WINDOW_KEYS = (
    "interval",
    "kind",
    "normal_energy",
    "source_probability",
    "source_probability_per_proper_length",
    "proper_pressure_work",
    "momentum_lapse_work",
    "coordinate_metric_work",
    "flux_projection_shift",
    "flux_projection_proper",
    "flux_projection_shift_pressure",
    "normal_boundary_flux",
    "normal_left_outward_flux",
    "normal_right_outward_flux",
    "probability_boundary_flux",
    "boundary_flux_identity_gap",
    "energy_reynolds",
    "probability_reynolds",
    "surface_work",
)
CURVATURE_KEYS = (
    "R_h_min",
    "R_h_max",
    "R_h_l2",
    "weyl_C2_max",
    "weyl_C2_l2",
    "chi_shell_gap_max",
    "Qddot_projection_gap_max",
    "Q_dot_max",
    "Q_ddot_max",
    "frozen_W",
    "auxiliary_substituted",
)
PROPER_KEYS = (
    "parent_proper_length",
    "child_proper_length",
    "parent_annulus_proper_length",
    "proper_length_ratio",
    "child_r_proper_mean",
    "parent_annulus_r_proper_mean",
    "child_Q_proper_mean",
    "parent_annulus_Q_proper_mean",
    "clock_locations",
    "clock_rates",
    "proper_clocks_integrated",
)
LOCALIZATION_KEYS = (
    "peak_x",
    "peak_per_proper_length",
    "arc_location",
    "resultant",
    "proper_width",
    "width_is_dispersal_time",
    "q_multiplied_into_measure",
)
CHARACTERISTIC_KEYS = (
    "coordinate_speed",
    "coordinate_speed_abs",
    "crossing_times",
    "proper_clocks_integrated",
    "dispersal_time_asserted",
    "candidate_tracks",
    "ambiguities",
    "packet_width_used_as_boundary",
)
TRACK_KEYS = (
    "boundary",
    "coordinate_speed",
    "coordinate",
    "kind",
    "measured_boundary",
)
MODE_KEYS = (
    "kind",
    "not_spatial_windows",
    "frame",
    "geometry_child_detail_norm_r",
    "geometry_parent_norm_r",
    "source_child_column_power",
    "source_parent_detail_column_power",
)
CONTROL_KEYS = ("status", "saved_control_applied", "control", "reason")
BOUNDARY_KEYS = (
    "x",
    "normal_energy_flux_density",
    "probability_flux_density",
    "window_outward_normal_flux",
    "window",
    "side",
)
PROFILE_KEYS = (
    "x",
    "source_probability_density",
    "source_probability_per_proper_length",
    "normal_energy_density",
    "normal_energy_flux",
    "probability_flux",
    "normal_interior_source",
    "Q_dot",
    "Q_ddot",
    "R_h",
    "weyl_C2",
    "chi",
    "chi_shell_gap",
    "r",
    "Q",
    "proper_radial_metric",
    "proper_gradient_source_probability",
    "proper_gradient_radius",
)
GRADIENT_KEYS = (
    "source_probability_per_proper_length_max_abs",
    "areal_radius_max_abs",
    "at_boundaries",
)


def observation_due(step_index, *, every=1, on_rk_stage=False):
    """Accepted steps on the caller's stride only. A Runge--Kutta stage is never due."""
    if on_rk_stage or OBSERVATION_ON_RK_STAGE:
        return False
    stride = int(every)
    if stride < 1 or int(step_index) != step_index or step_index < 0:
        raise ValueError("observation stride needs a nonnegative integer step and a positive stride")
    return int(step_index) % stride == 0


def frozen_geometry_policy(preparation):
    """Saved replay basis applies only to the sealed nested-pair records.

    A new pair and the original conformal episode have no frozen ``W``
    control. This does not rebuild a geometry frame.
    """
    name = str(preparation)
    if name not in PREPARATION_NAMES:
        raise ValueError("unknown discovery preparation")
    if name in ("separated_pair_v1", "separated_pair_v2"):
        return {
            "preparation": name,
            "status": "saved_control",
            "control": "nsc-nested-parent-child-replay-basis-v1",
            "reason": "compare this sealed nested-pair record only through the saved replay basis",
        }
    if name == "frozen_basis":
        return {
            "preparation": name,
            "status": "saved_control",
            "control": "nsc-nested-parent-child-replay-basis-v1",
            "reason": "this record is the saved geometry control",
        }
    if name == "original_conformal_v2":
        reason = "the original conformal episode is not a nested pair and has no frozen W control"
    else:
        reason = "a newly built pair has no saved frozen-geometry control"
    return {"preparation": name, "status": "missing", "control": None, "reason": reason}


def reference_weights(count, length, start, end):
    """Same weights as the conformal frame ledger: interpolant integral divided by dx."""
    count, length, start, end = int(count), float(length), float(start), float(end)
    if not 0.0 <= start < end <= length:
        raise ValueError("ordered interval inside one period required")
    wave = 2.0 * np.pi * np.fft.fftfreq(count, d=length / count)
    integrals = np.empty(count, dtype=complex)
    nonzero = wave != 0.0
    integrals[nonzero] = (np.exp(1j * wave[nonzero] * end) - np.exp(1j * wave[nonzero] * start)) / (1j * wave[nonzero])
    integrals[~nonzero] = end - start
    return (np.fft.fft(integrals) / length).real


def reference_window_integral(grid, summand, interval):
    """Reference quadrature of a nodal summand. A full period equals the raw sum."""
    samples = np.asarray(summand, dtype=float)
    if samples.shape != (grid.nq,) or not np.isfinite(samples).all():
        raise ValueError("a reference summand must be a finite nq-vector")
    weights = reference_weights(grid.nq, grid.length, interval[0], interval[1])
    return float(np.dot(weights, samples))


def moving_window_transport(grid, eta, flux, interior_source, interval, boundary_velocity):
    """Leibniz rate of ``E = ∫ eta dx`` on an artificial coordinate cut.

    The Reynolds term is ``eta(b) right_dot - eta(a) left_dot``. Volumetric
    pressure and lapse work stay in ``interior_source``. This cut does not
    add ``-p V``.
    """
    left, right = (float(interval[0]), float(interval[1]))
    left_dot, right_dot = (float(boundary_velocity[0]), float(boundary_velocity[1]))
    if not 0.0 <= left < right <= float(grid.length):
        raise ValueError("moving window must be an ordered interval inside one period")
    if not all(math.isfinite(value) for value in (left_dot, right_dot)):
        raise ValueError("boundary velocity must be finite")
    fields = {}
    for name, values in (("eta", eta), ("flux", flux), ("interior_source", interior_source)):
        array = np.asarray(values, dtype=float)
        if array.shape != (grid.nq,) or not np.isfinite(array).all():
            raise ValueError(f"{name} must be a finite nodal density")
        fields[name] = array
    sampled = {
        name: nested._periodic_values(grid, fields[name], (left, right))
        for name in ("eta", "flux")
    }
    reynolds = float(sampled["eta"][1] * right_dot - sampled["eta"][0] * left_dot)
    boundary_flux = float(sampled["flux"][0] - sampled["flux"][1])
    interior = nested.interval_integral(grid, fields["interior_source"], (left, right))
    return {
        "interval": [left, right],
        "boundary_velocity": [left_dot, right_dot],
        "reynolds": reynolds,
        "boundary_flux": boundary_flux,
        "left_outward_flux": float(-sampled["flux"][0]),
        "right_outward_flux": float(sampled["flux"][1]),
        "interior_source": interior,
        "leibniz_rate": interior + boundary_flux + reynolds,
        "surface_work": 0.0,
        "fixed_window": left_dot == 0.0 and right_dot == 0.0,
    }


def surface_action_work(grid, traction, interval, boundary_velocity):
    """Explicit surface traction on a moving cut. Not added to volumetric ``-p dV``.

    ``traction`` is the caller's surface-action density. The artificial-window
    Reynolds balance does not call this, so a gradient cut is not charged
    ``-p V`` on top of the observer pressure work.
    """
    left, right = (float(interval[0]), float(interval[1]))
    left_dot, right_dot = (float(boundary_velocity[0]), float(boundary_velocity[1]))
    samples = np.asarray(traction, dtype=float)
    if samples.shape != (grid.nq,) or not np.isfinite(samples).all():
        raise ValueError("surface traction must be a finite nodal density")
    if not 0.0 <= left < right <= float(grid.length):
        raise ValueError("surface interval must lie inside one period")
    endpoints = nested._periodic_values(grid, samples, (left, right))
    return float(endpoints[1] * right_dot - endpoints[0] * left_dot)


def candidate_boundary_tracks(time, length=CARRIER_LENGTH):
    """Null rays of ``dx/dt = ±1`` from the fixed cuts. Not measured boundaries."""
    instant = float(time)
    period = float(length)
    if not math.isfinite(instant) or not math.isfinite(period) or period <= 0.0:
        raise ValueError("a finite time and a positive period are required")
    tracks = []
    for boundary in BOUNDARY_ORDER:
        for speed in CHARACTERISTIC_SPEEDS:
            tracks.append({
                "boundary": float(boundary),
                "coordinate_speed": float(speed),
                "coordinate": float((boundary + speed * instant) % period),
                "kind": "candidate_characteristic",
                "measured_boundary": False,
            })
    ambiguities = []
    for track in tracks:
        for other in BOUNDARY_ORDER:
            if other == track["boundary"]:
                continue
            separation = abs((track["coordinate"] - other + 0.5 * period) % period - 0.5 * period)
            if separation <= 1e-8:
                ambiguities.append({
                    "track_boundary": track["boundary"],
                    "coordinate_speed": track["coordinate_speed"],
                    "coordinate": track["coordinate"],
                    "coincides_with": float(other),
                    "kind": "characteristic_image_of_another_fixed_boundary",
                })
    return tracks, ambiguities


def _curvature_helper():
    try:
        import derive_nsc_spherical_conformal_curvature as curvature
    except ImportError as error:
        raise ImportError(
            "live curvature reuses derive_nsc_spherical_conformal_curvature; "
            "import it with lab/scripts on PYTHONPATH"
        ) from error
    return curvature


def _require_pair(pair):
    if pair.grid.gauge != "conformal":
        raise ValueError("discovery observation requires the conformal gauge")
    if pair.geometry_map.flags.writeable:
        raise ValueError("observation requires the frozen read-only geometry map")
    if abs(float(pair.grid.length) - CARRIER_LENGTH) > 1e-12:
        raise ValueError("discovery characteristics assume the period-8 carrier")
    if tuple(float(value) for value in pair.parent_interval) != PARENT_INTERVAL:
        raise ValueError("parent window is the fixed coordinate interval (0, 4)")
    if tuple(float(value) for value in pair.child_interval) != CHILD_INTERVAL:
        raise ValueError("child window is the fixed coordinate interval (1, 3)")


def _bundle_matches(pair, nodal, bundle, nodal_rate):
    fine = bundle["fine_state"]
    expected = galerkin.prolong_state(pair.grid, nodal)
    gap = 0.0
    for name in ("Q", "r", "chi", "p_Q", "p_r", "p_chi"):
        gap = max(gap, float(np.max(np.abs(getattr(fine, name) - getattr(expected, name)))))
    gap = max(gap, float(np.max(np.abs(fine.phi0 - expected.phi0))), float(np.max(np.abs(fine.phi1 - expected.phi1))))
    if gap > 1e-8:
        raise ValueError("bundle is not the frozen-W image of this nested state")
    if np.asarray(nodal_rate.Q).shape != (pair.grid.ng,):
        raise ValueError("nodal_rate must be the coarse Cauchy rate")
    lifted = galerkin.prolong_geometry(pair.grid, nodal_rate.Q)
    if float(np.max(np.abs(lifted - bundle["lifted_Q"]))) > 1e-8:
        raise ValueError("nodal_rate is not the coarse rate that produced this bundle")


def _prepare(pair, state, bundle, nodal_rate):
    _require_pair(pair)
    nested._check(pair, state)
    nodal = nested.reconstruct_state(pair, state)
    if bundle is None and nodal_rate is None:
        nodal_rate, bundle = galerkin.compose_fine_hamiltonian(pair.grid, nodal)
    elif bundle is None or nodal_rate is None:
        raise ValueError("bundle and nodal_rate are supplied together")
    else:
        _bundle_matches(pair, nodal, bundle, nodal_rate)
    return nodal, nodal_rate, bundle


def _lifted_rate(grid, nodal_rate, bundle):
    source = bundle["source"]
    return CauchyRate(
        Q=galerkin.prolong_geometry(grid, nodal_rate.Q),
        r=galerkin.prolong_geometry(grid, nodal_rate.r),
        chi=galerkin.prolong_geometry(grid, nodal_rate.chi),
        p_Q=galerkin.prolong_geometry(grid, nodal_rate.p_Q),
        p_r=galerkin.prolong_geometry(grid, nodal_rate.p_r),
        p_chi=galerkin.prolong_geometry(grid, nodal_rate.p_chi),
        phi0=galerkin.prolong_columns(grid, nodal_rate.phi0),
        phi1=galerkin.prolong_columns(grid, nodal_rate.phi1),
        force_L=source["force_L"],
        force_Q=source["force_Q"],
        force_beta=source["force_beta"],
        fieldwork_power=float(nodal_rate.fieldwork_power),
    )


def _live_curvature(pair, nodal, nodal_rate, bundle, control_mode):
    fine = bundle["fine_state"]
    zeros = np.zeros(pair.grid.nq)
    if control_mode == "frozen_geometry":
        q_dot, q_ddot, before_projection = zeros, zeros, zeros
    else:
        helper = _curvature_helper()
        q_dot, q_ddot, before_projection = helper.actual_q_second_rate(pair.grid, nodal, nodal_rate, bundle)
    direct = metric.direct_rh_grid(
        fine.Q, fine.Q, zeros, q_dot, q_ddot, pair.grid.length, L_dot=q_dot, beta_dot=zeros,
    )
    curvature = np.asarray(direct["R_h"], dtype=float)
    weyl = np.asarray(metric.weyl_actual(curvature, fine.r), dtype=float)
    gap = curvature - np.asarray(fine.chi, dtype=float) - 2.0
    spacing = float(pair.grid.dx_q)
    report = {
        "R_h_min": float(np.min(curvature)),
        "R_h_max": float(np.max(curvature)),
        "R_h_l2": float(np.sqrt(spacing * np.sum(curvature ** 2))),
        "weyl_C2_max": float(np.max(weyl)),
        "weyl_C2_l2": float(np.sqrt(spacing * np.sum(weyl ** 2))),
        "chi_shell_gap_max": float(np.max(np.abs(gap))),
        "Qddot_projection_gap_max": float(np.max(np.abs(q_ddot - before_projection))),
        "Q_dot_max": float(np.max(np.abs(q_dot))),
        "Q_ddot_max": float(np.max(np.abs(q_ddot))),
        "frozen_W": True,
        "auxiliary_substituted": False,
    }
    return report, {"R_h": curvature, "weyl_C2": weyl, "chi_shell_gap": gap, "Q_dot": np.asarray(q_dot, dtype=float), "Q_ddot": np.asarray(q_ddot, dtype=float)}


def _probability(system, fine):
    weights = np.asarray(system.occupations, dtype=float)
    mass = np.sum((np.abs(fine.phi0) ** 2 + np.abs(fine.phi1) ** 2) * weights[None, :], axis=1)
    flux_samples = 2.0 * np.sum(np.imag(np.conjugate(fine.phi0) * fine.phi1) * weights[None, :], axis=1)
    spacing = float(system.dx)
    return mass, mass / spacing, flux_samples / spacing


def _column_power(frame, columns):
    coefficients = np.asarray(frame, dtype=complex).conj().T @ np.asarray(columns, dtype=complex)
    return float(np.sum(np.abs(coefficients) ** 2))


def _mode_projectors(pair, state, measured):
    live = np.vstack((state.phi0, state.phi1))
    child = pair.original_columns[:, pair.child_indices]
    detail_ids = [index for index in range(6) if index not in pair.child_indices]
    detail = pair.original_columns[:, detail_ids]
    return {
        "kind": "canonical_mode_projector",
        "not_spatial_windows": True,
        "frame": "fixed_original_observer",
        "geometry_child_detail_norm_r": float(measured["r_child_detail_norm"]),
        "geometry_parent_norm_r": float(measured["r_parent_norm"]),
        "source_child_column_power": _column_power(child, live),
        "source_parent_detail_column_power": _column_power(detail, live),
    }


def _localization(grid, mass, per_proper, proper_q):
    packet = metric.packet_geometry(grid.xi_q, grid.length, mass, proper_q)
    resultant = float(packet["resultant"])
    if resultant > metric.NUMERICAL_RESULTANT_FLOOR:
        width = float(packet["proper_length"] * math.sqrt(-2.0 * math.log(min(resultant, 1.0))) / (2.0 * math.pi))
    else:
        width = None
    peak_index = int(np.argmax(per_proper))
    return {
        "peak_x": float(grid.xi_q[peak_index]),
        "peak_per_proper_length": float(per_proper[peak_index]),
        "arc_location": None if packet["arc_location"] is None else float(packet["arc_location"]),
        "resultant": resultant,
        "proper_width": width,
        "width_is_dispersal_time": False,
        "q_multiplied_into_measure": False,
    }


def _collocation_flux(grid, summand, interval):
    """Reference boundary flux ``(sample[a] - sample[b]) / dx``. Nyquist stays in the sample."""
    left, right = (float(interval[0]), float(interval[1]))
    spacing = float(grid.dx_q)
    left_index = int(round(left / spacing)) % grid.nq
    right_index = int(round(right / spacing)) % grid.nq
    samples = np.asarray(summand, dtype=float)
    net = float((samples[left_index] - samples[right_index]) / spacing)
    return net, float(-samples[left_index] / spacing), float(samples[right_index] / spacing)


def _flux_identity_gap(grid, summand, interval, net):
    """Derivative-matrix window versus collocation. The Nyquist symbol of the matrix stays zero."""
    weights = reference_weights(grid.nq, grid.length, interval[0], interval[1])
    discrete = float(np.dot(grid.derivative @ weights, np.asarray(summand, dtype=float)))
    return abs(net - discrete)


def _window(grid, interval, *, energy_summand, energy_flux_summand, probability_summand, probability_flux_summand,
            pressure_summand, lapse_summand, coordinate_work, projections, proper_length):
    held = (0.0, 0.0)
    spacing = float(grid.dx_q)
    energy = moving_window_transport(
        grid, energy_summand / spacing, energy_flux_summand / spacing,
        (pressure_summand + lapse_summand) / spacing, interval, held,
    )
    probability = moving_window_transport(
        grid, probability_summand / spacing, probability_flux_summand / spacing,
        np.zeros(grid.nq), interval, held,
    )
    probability_mass = reference_window_integral(grid, probability_summand, interval)
    net, left_out, right_out = _collocation_flux(grid, energy_flux_summand, interval)
    probability_net, _, _ = _collocation_flux(grid, probability_flux_summand, interval)
    return {
        "interval": [float(interval[0]), float(interval[1])],
        "kind": "fixed_coordinate_window",
        "normal_energy": reference_window_integral(grid, energy_summand, interval),
        "source_probability": probability_mass,
        "source_probability_per_proper_length": probability_mass / float(proper_length),
        "proper_pressure_work": reference_window_integral(grid, pressure_summand, interval),
        "momentum_lapse_work": reference_window_integral(grid, lapse_summand, interval),
        "coordinate_metric_work": reference_window_integral(grid, coordinate_work, interval),
        "flux_projection_shift": reference_window_integral(grid, projections["shift"], interval),
        "flux_projection_proper": reference_window_integral(grid, projections["proper"], interval),
        "flux_projection_shift_pressure": reference_window_integral(grid, projections["shift_pressure"], interval),
        "normal_boundary_flux": net,
        "normal_left_outward_flux": left_out,
        "normal_right_outward_flux": right_out,
        "probability_boundary_flux": probability_net,
        "boundary_flux_identity_gap": _flux_identity_gap(grid, energy_flux_summand, interval, net),
        "energy_reynolds": energy["reynolds"],
        "probability_reynolds": probability["reynolds"],
        "surface_work": energy["surface_work"],
    }


def _boundary_catalog(grid, energy_flux, probability_flux):
    energy = nested._periodic_values(grid, energy_flux, BOUNDARY_ORDER)
    probability = nested._periodic_values(grid, probability_flux, BOUNDARY_ORDER)
    rows = []
    for index, position in enumerate(BOUNDARY_ORDER):
        if position in (0.0, 1.0):
            outward = float(-energy[index])
            side = "left"
        else:
            outward = float(energy[index])
            side = "right"
        rows.append({
            "x": float(position),
            "normal_energy_flux_density": float(energy[index]),
            "probability_flux_density": float(probability[index]),
            "window_outward_normal_flux": outward,
            "window": "parent" if position in (0.0, 4.0) else "child",
            "side": side,
        })
    return rows


def _profiles(grid, fields):
    return {name: np.array(fields[name], copy=True) for name in PROFILE_KEYS if name != "x"} | {
        "x": np.array(grid.xi_q, copy=True),
    }


def _reject_arrays(value, path):
    if isinstance(value, np.ndarray):
        raise ValueError(f"compact row contains an array at {path}")
    if isinstance(value, np.generic):
        raise ValueError(f"compact row contains a numpy scalar at {path}")
    if isinstance(value, dict):
        for key, item in value.items():
            _reject_arrays(item, path + "." + str(key))
    elif isinstance(value, (list, tuple)):
        for index, item in enumerate(value):
            _reject_arrays(item, f"{path}[{index}]")
    elif not isinstance(value, (float, int, str, bool, type(None))):
        raise ValueError(f"compact row contains {type(value).__name__} at {path}")


def assert_compact_row(row):
    """Fail if a campaign row drifts from the stage-1 scalar contract."""
    missing = [key for key in REQUIRED_TOP_KEYS if key not in row]
    if missing:
        raise ValueError("compact row is missing " + ", ".join(missing))
    if row["api_version"] != API_VERSION:
        raise ValueError("compact row API version differs from the frozen observer")
    if row["capture_profiles"]:
        profiles = row.get("profiles")
        if not isinstance(profiles, dict) or set(profiles) != set(PROFILE_KEYS):
            raise ValueError("profile capture must contain the frozen profile fields")
        for name in PROFILE_KEYS:
            if not isinstance(profiles[name], np.ndarray):
                raise ValueError("profile fields must be arrays")
        scalar = {key: value for key, value in row.items() if key != "profiles"}
    else:
        if "profiles" in row:
            raise ValueError("profiles were attached without capture_profiles")
        scalar = row
    _reject_arrays(scalar, "row")
    flags = {
        "observation_on_rk_stage": False,
        "spatial_windows_are_mode_projectors": False,
        "normal_balance_slope_evaluated": False,
        "chi_proxy_substituted": False,
        "renewal_asserted": False,
        "continuum_certified": False,
        "imposed_boundary_trajectory_is_regeneration_evidence": False,
    }
    for key, expected in flags.items():
        if row[key] is not expected:
            raise ValueError(f"{key} must be {expected}")
    for name, keys, interval in (
        ("child_window", WINDOW_KEYS, [1.0, 3.0]),
        ("parent_window", WINDOW_KEYS, [0.0, 4.0]),
    ):
        if list(row[name]) != list(keys) and set(row[name]) != set(keys):
            raise ValueError(f"{name} keys differ from the frozen window")
        if row[name]["interval"] != interval or row[name]["kind"] != "fixed_coordinate_window":
            raise ValueError(f"{name} is not the fixed coordinate window")
        if row[name]["energy_reynolds"] != 0.0 or row[name]["probability_reynolds"] != 0.0 or row[name]["surface_work"] != 0.0:
            raise ValueError("the fixed observe window has no Reynolds term or surface work")
    if set(row["curvature"]) != set(CURVATURE_KEYS) or row["curvature"]["auxiliary_substituted"] is not False:
        raise ValueError("curvature block drifted from the frozen jet comparison")
    if row["curvature"]["frozen_W"] is not True:
        raise ValueError("curvature must be converted through the frozen geometry map")
    if set(row["proper"]) != set(PROPER_KEYS) or row["proper"]["proper_clocks_integrated"] is not False:
        raise ValueError("proper clocks are reported as rates and are not integrated here")
    if set(row["localization"]) != set(LOCALIZATION_KEYS) or row["localization"]["width_is_dispersal_time"] is not False:
        raise ValueError("localization width must not be reported as a dispersal time")
    if set(row["characteristics"]) != set(CHARACTERISTIC_KEYS):
        raise ValueError("characteristic block drifted")
    crossings = row["characteristics"]["crossing_times"]
    if crossings != CROSSING_TIMES or row["characteristics"]["coordinate_speed"] != [-1.0, 1.0]:
        raise ValueError("coordinate crossings must stay L8, parent4, child2 at dx/dt = ±1")
    if row["characteristics"]["dispersal_time_asserted"] is not False or row["characteristics"]["packet_width_used_as_boundary"] is not False:
        raise ValueError("packet width is not a boundary and not a dispersal time")
    if len(row["characteristics"]["candidate_tracks"]) != 8:
        raise ValueError("eight candidate characteristic tracks are required")
    for track in row["characteristics"]["candidate_tracks"]:
        if set(track) != set(TRACK_KEYS) or track["measured_boundary"] is not False:
            raise ValueError("a candidate track was marked as a measured boundary")
    if set(row["mode_projectors"]) != set(MODE_KEYS) or row["mode_projectors"]["not_spatial_windows"] is not True:
        raise ValueError("mode projectors must stay distinct from spatial windows")
    if set(row["frozen_geometry_control"]) != set(CONTROL_KEYS):
        raise ValueError("frozen geometry control block drifted")
    if row["control_mode"] not in ("coupled", "frozen_geometry"):
        raise ValueError("control_mode must be coupled or frozen_geometry")
    if row["frozen_geometry_control"]["status"] not in ("missing", "saved_control", "unspecified"):
        raise ValueError("geometry-control status is not one of the preparations")
    if row["control_mode"] == "frozen_geometry" and row["curvature"]["Q_dot_max"] != 0.0:
        raise ValueError("frozen geometry keeps Qdot at zero")
    if row["control_mode"] == "frozen_geometry" and row["coordinate_metric_work_total"] != 0.0:
        raise ValueError("frozen geometry has no pressure metric work")
    if [item["x"] for item in row["boundary_flux"]] != list(BOUNDARY_ORDER):
        raise ValueError("boundary flux must be catalogued at x = 0, 1, 3, 4")
    for item in row["boundary_flux"]:
        if set(item) != set(BOUNDARY_KEYS):
            raise ValueError("boundary flux record drifted")
    if set(row["proper_gradients"]) != set(GRADIENT_KEYS):
        raise ValueError("proper gradient block drifted")
    return row


def _control_report(preparation):
    if preparation is None:
        return {
            "status": "unspecified",
            "saved_control_applied": False,
            "control": None,
            "reason": "no preparation was named; the saved replay basis is not assumed",
        }
    policy = frozen_geometry_policy(preparation)
    return {
        "status": policy["status"],
        "saved_control_applied": False,
        "control": policy["control"],
        "reason": policy["reason"],
    }


def observe(pair, state, time, *, capture_profiles=False, bundle=None, nodal_rate=None,
            preparation=None, control_mode="coupled"):
    """Return one compact scalar row for this nested state at ``time``.

    ``bundle`` and ``nodal_rate``, when both supplied, must already be the
    frozen-``W`` image and the coarse Cauchy rate of this state. ``coupled``
    uses those jets. ``frozen_geometry`` keeps the current fields and sets
    the geometry jets, including ``Qdot`` and ``Qddot``, to zero. Profiles
    are omitted unless ``capture_profiles`` is true.
    """
    if control_mode not in ("coupled", "frozen_geometry"):
        raise ValueError("control_mode must be coupled or frozen_geometry")
    instant = float(time)
    if not math.isfinite(instant):
        raise ValueError("observation time must be finite")
    _nodal, nodal_rate, bundle = _prepare(pair, state, bundle, nodal_rate)
    grid = pair.grid
    fine = bundle["fine_state"]
    system = bundle["fine_system"]
    if min(float(np.min(fine.r)), float(np.min(fine.Q))) <= 0.0:
        raise ValueError("observation requires a positive areal radius and radial metric")
    lifted = _lifted_rate(grid, nodal_rate, bundle)
    if control_mode == "frozen_geometry":
        zeros = np.zeros(grid.nq)
        lifted = CauchyRate(
            zeros, zeros, zeros, zeros, zeros, zeros, lifted.phi0, lifted.phi1,
            lifted.force_L, lifted.force_Q, lifted.force_beta, 0.0,
        )
    ledger = regional.matter_ledger(system, fine)
    terms = regional.proper_balance_terms(system, fine, lifted, ledger)
    mass, probability_density, probability_flux_density = _probability(system, fine)
    probability_flux_summand = probability_flux_density * float(system.dx)
    proper_q = np.asarray(fine.r, dtype=float) * np.asarray(fine.Q, dtype=float)
    per_proper = probability_density / proper_q
    if control_mode == "frozen_geometry":
        coordinate_work = np.zeros(grid.nq)
    else:
        coordinate_work = bundle["source"]["force_Q"] * bundle["lifted_Q"] + bundle["source"]["force_L"] * bundle["lapse_dot"]
    projections = {
        "shift": np.asarray(ledger["flux_shift"], dtype=float),
        "proper": np.asarray(ledger["flux_proper"], dtype=float),
        "shift_pressure": np.asarray(ledger["flux_shift_pressure"], dtype=float),
    }
    carrier = (0.0, float(grid.length))
    measured = nested.metrics(pair, state)
    curvature, curvature_profiles = _live_curvature(pair, _nodal, nodal_rate, bundle, control_mode)
    probability_gradient = (system.derivative @ per_proper) / proper_q
    radius_gradient = (system.derivative @ np.asarray(fine.r, dtype=float)) / proper_q
    gradient_at_cuts = nested._periodic_values(grid, probability_gradient, BOUNDARY_ORDER)
    tracks, ambiguities = candidate_boundary_tracks(instant, grid.length)
    peak_index = int(np.argmax(per_proper))
    nyquist = np.cos(np.pi * np.arange(grid.nq))
    window_arguments = dict(
        energy_summand=np.asarray(ledger["normal_energy_nodal"], dtype=float),
        energy_flux_summand=np.asarray(terms["flux_nodal"], dtype=float),
        probability_summand=np.asarray(mass, dtype=float),
        probability_flux_summand=np.asarray(probability_flux_summand, dtype=float),
        pressure_summand=np.asarray(terms["proper_pressure_work"], dtype=float),
        lapse_summand=np.asarray(terms["momentum_lapse_work"], dtype=float),
        coordinate_work=np.asarray(coordinate_work, dtype=float),
        projections=projections,
    )
    row = {
        "api_version": API_VERSION,
        "time": instant,
        "capture_profiles": bool(capture_profiles),
        "observation_on_rk_stage": False,
        "gauge": "conformal",
        "spatial_windows_are_mode_projectors": False,
        "source_probability_total": float(np.sum(mass)),
        "source_probability_interpolant": reference_window_integral(grid, mass, carrier),
        "source_probability_per_proper_length_peak": float(per_proper[peak_index]),
        "source_probability_per_proper_length_peak_x": float(grid.xi_q[peak_index]),
        "normal_energy_total": reference_window_integral(grid, window_arguments["energy_summand"], carrier),
        "coordinate_metric_work_total": reference_window_integral(grid, window_arguments["coordinate_work"], carrier),
        "flux_projection_shift_total": reference_window_integral(grid, projections["shift"], carrier),
        "flux_projection_proper_total": reference_window_integral(grid, projections["proper"], carrier),
        "flux_projection_shift_pressure_total": reference_window_integral(grid, projections["shift_pressure"], carrier),
        "nyquist_symbol_max": float(np.max(np.abs(grid.derivative @ nyquist))),
        "child_window": _window(grid, CHILD_INTERVAL, proper_length=measured["child_proper_length"], **window_arguments),
        "parent_window": _window(grid, PARENT_INTERVAL, proper_length=measured["parent_proper_length"], **window_arguments),
        "boundary_flux": _boundary_catalog(grid, window_arguments["energy_flux_summand"] / grid.dx_q, probability_flux_density),
        "curvature": curvature,
        "proper": {
            "parent_proper_length": float(measured["parent_proper_length"]),
            "child_proper_length": float(measured["child_proper_length"]),
            "parent_annulus_proper_length": float(measured["parent_annulus_proper_length"]),
            "proper_length_ratio": float(measured["proper_length_ratio"]),
            "child_r_proper_mean": float(measured["child_r_proper_mean"]),
            "parent_annulus_r_proper_mean": float(measured["parent_annulus_r_proper_mean"]),
            "child_Q_proper_mean": float(measured["child_Q_proper_mean"]),
            "parent_annulus_Q_proper_mean": float(measured["parent_annulus_Q_proper_mean"]),
            "clock_locations": [float(value) for value in measured["clock_locations"]],
            "clock_rates": [float(value) for value in measured["clock_rates"]],
            "proper_clocks_integrated": False,
        },
        "localization": _localization(grid, mass, per_proper, proper_q),
        "proper_gradients": {
            "source_probability_per_proper_length_max_abs": float(np.max(np.abs(probability_gradient))),
            "areal_radius_max_abs": float(np.max(np.abs(radius_gradient))),
            "at_boundaries": [
                {"x": float(position), "source_probability_gradient": float(gradient_at_cuts[index])}
                for index, position in enumerate(BOUNDARY_ORDER)
            ],
        },
        "characteristics": {
            "coordinate_speed": [-1.0, 1.0],
            "coordinate_speed_abs": 1.0,
            "crossing_times": dict(CROSSING_TIMES),
            "proper_clocks_integrated": False,
            "dispersal_time_asserted": False,
            "candidate_tracks": tracks,
            "ambiguities": ambiguities,
            "packet_width_used_as_boundary": False,
        },
        "mode_projectors": _mode_projectors(pair, state, measured),
        "preparation": None if preparation is None else str(preparation),
        "control_mode": control_mode,
        "frozen_geometry_control": _control_report(preparation),
        "normal_balance_slope_evaluated": False,
        "chi_proxy_substituted": False,
        "renewal_asserted": False,
        "continuum_certified": False,
        "imposed_boundary_trajectory_is_regeneration_evidence": False,
    }
    if capture_profiles:
        row["profiles"] = _profiles(grid, {
            "source_probability_density": probability_density,
            "source_probability_per_proper_length": per_proper,
            "normal_energy_density": window_arguments["energy_summand"] / grid.dx_q,
            "normal_energy_flux": window_arguments["energy_flux_summand"] / grid.dx_q,
            "probability_flux": probability_flux_density,
            "normal_interior_source": (window_arguments["pressure_summand"] + window_arguments["lapse_summand"]) / grid.dx_q,
            "Q_dot": curvature_profiles["Q_dot"],
            "Q_ddot": curvature_profiles["Q_ddot"],
            "R_h": curvature_profiles["R_h"],
            "weyl_C2": curvature_profiles["weyl_C2"],
            "chi": np.asarray(fine.chi, dtype=float),
            "chi_shell_gap": curvature_profiles["chi_shell_gap"],
            "r": np.asarray(fine.r, dtype=float),
            "Q": np.asarray(fine.Q, dtype=float),
            "proper_radial_metric": proper_q,
            "proper_gradient_source_probability": probability_gradient,
            "proper_gradient_radius": radius_gradient,
        })
    return assert_compact_row(row)
