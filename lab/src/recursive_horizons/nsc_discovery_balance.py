"""RK4-stage normal-energy ledger on the existing nested conformal step.

``model.rk4_step`` and the ODE stay where they are. ``rk4_stage_ledger``
repeats that step's arithmetic. Each stage calls ``rates(..., return_bundle=True)``
once and uses that rate both to advance the state and to sample the balance.
The updated state is the owned Runge--Kutta state.

The accounted windows are the disjoint fixed cuts left ``(0, 1)``, child
``(1, 3)``, right ``(3, 4)`` and ambient ``(4, 8)``. Parent is the sum of the
first three. Global is the sum of all four. Adding the child window to the
overlapping parent window ``(0, 4)`` is recorded and is not used.

The instantaneous balance is the regions helper: collocation flux, pressure,
lapse, the measured pointwise projection defect, and the gap between that
collocation flux and the zero-Nyquist derivative window. Coordinate field
work and those rates are integrated with the RK4 weights, not with samples
spaced by ``0.05``. A gradient measurement cut is not charged an extra
``-p V``. ``M = 4κ`` stays inside the nodal forces once, and each summand is
divided by ``dx`` once.

Proper clocks remain ``dτ/dT = rQ`` at ``x = 1, 2, 3``, accumulated with the
episode trapezoid. The same locations sampled on the four stage states are a
comparison, not a replacement protocol.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np

from . import nsc_discovery_episode as episode
from . import nsc_discovery_observables as observer
from . import nsc_discovery_regions as regions
from . import nsc_nested_parent_child as model
from . import nsc_regional_energy_exchange as regional
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin


API_VERSION = "nsc-discovery-balance-v1"
SCHEMA = "NSC-DISCOVERY-BALANCE-v1"
CHECK_SCHEMA = "NSC-DISCOVERY-BALANCE-CHECK-v1"
CARRIER_LENGTH = 8.0
HANDOFF_TIME = episode.HANDOFF_TIME
TARGET_TIME = 1.0
CPU_BUDGET_SECONDS = 6 * 3600.0
CHUNK_LIMIT_BYTES = 64 * 1024 * 1024
CLOCK_LOCATIONS = (1.0, 2.0, 3.0)
STEP_CAPS = (episode.STEP_CAP, episode.MATCHED_STEP_CAP)
WINDOW_SPECS = (
    ("left", (0.0, 1.0)),
    ("child", (1.0, 3.0)),
    ("right", (3.0, 4.0)),
    ("ambient", (4.0, CARRIER_LENGTH)),
)
DISJOINT_WINDOWS = tuple(name for name, _interval in WINDOW_SPECS)
PARENT_PARTS = ("left", "child", "right")
DIRECT_PARENT_INTERVAL = (0.0, 4.0)
RESOLVED_CUTS = (
    {"x": 0.0, "left_window": "ambient", "right_window": "left"},
    {"x": 1.0, "left_window": "left", "right_window": "child"},
    {"x": 3.0, "left_window": "child", "right_window": "right"},
    {"x": 4.0, "left_window": "right", "right_window": "ambient"},
)
RATE_KEYS = (
    "fixed_boundary_flux",
    "derivative_window_flux",
    "nyquist_gap",
    "pointwise_projection_defect",
    "finite_projection_defect",
    "pressure",
    "lapse",
    "reynolds",
    "slope_integral",
    "balance_rate",
    "coordinate_metric_work",
    "probability_boundary_flux",
    "probability_derivative_window_flux",
    "probability_nyquist_gap",
    "probability_gross_in_lower_bound",
    "probability_gross_out_lower_bound",
    "probability_slope_integral",
    "probability_projection_defect",
    "left_outward_proper",
    "right_outward_proper",
    "left_inward_proper",
    "right_inward_proper",
    "left_outward_probability",
    "right_outward_probability",
    "left_inward_probability",
    "right_inward_probability",
    "gross_cut_magnitude",
    "instantaneous_closure_residual",
)
CUT_KEYS = (
    "proper_normal_flux",
    "sigma2_chiral_probability_flux",
    "left_window_outward_proper",
    "left_window_inward_proper",
    "right_window_outward_proper",
    "right_window_inward_proper",
    "left_window_outward_probability",
    "left_window_inward_probability",
    "right_window_outward_probability",
    "right_window_inward_probability",
    "outward_pair_sum_proper",
    "outward_pair_sum_probability",
)
LAB = Path(__file__).resolve().parents[2]
SEALED_INPUTS = (
    "results/development/nsc-nested-parent-child-v1.json",
    "results/development/nsc-nested-parent-child-v1.npz",
    "results/development/nsc-nested-parent-child-confirmation-v2.json",
    "results/development/nsc-nested-parent-child-confirmation-v2.npz",
    "results/development/nsc-nested-parent-child-replay-basis-v1.json",
    "results/development/nsc-nested-parent-child-replay-basis-v1.npz",
)
EPISODE_MANIFEST = "results/development/nsc-discovery-episode-v1/manifest.json"
PRODUCERS = (
    "src/recursive_horizons/nsc_discovery_balance.py",
    "scripts/derive_nsc_discovery_balance.py",
    "src/recursive_horizons/nsc_discovery_regions.py",
    "src/recursive_horizons/nsc_discovery_observables.py",
    "src/recursive_horizons/nsc_discovery_episode.py",
    "src/recursive_horizons/nsc_nested_parent_child.py",
    "src/recursive_horizons/nsc_regional_energy_exchange.py",
    "src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
    "src/recursive_horizons/nsc_discovery_backend.py",
    "src/recursive_horizons/nsc_spherical_coupling.py",
    "src/recursive_horizons/nsc_spherical_feedback_action.py",
    "src/recursive_horizons/nsc_conformal_adm_source.py",
    "scripts/derive_nsc_nested_parent_child.py",
)
FINE_HANDOFF = "nf256_baseline_dt0.0005"

CLOSURE_IDENTITY = {
    "instantaneous": (
        "d/dt ∫_a^b η dx = (F(a) - F(b))/dx + pressure + lapse "
        "+ pointwise_projection_defect + nyquist_gap + reynolds"
    ),
    "nyquist_gap": "nyquist_gap = -derivative_window_flux - collocation_flux",
    "derivative_window": "interval integral of (∂x F_summand) / dx; the derivative symbol at Nyquist is zero",
    "collocation_flux": "(F(a) - F(b))/dx on the Fourier interpolant",
    "step": (
        "E(t+dt) - E(t) = dt/6 * (f1 + 2*f2 + 2*f3 + f4) "
        "+ quadrature_closure_error"
    ),
    "rk_weight": "dt/6 * (f1 + 2*f2 + 2*f3 + f4)",
    "f_channels": (
        "collocation_flux",
        "pressure",
        "lapse",
        "pointwise_projection_defect",
        "nyquist_gap",
        "reynolds",
    ),
    "windows": {
        "left": [0.0, 1.0],
        "child": [1.0, 3.0],
        "right": [3.0, 4.0],
        "ambient": [4.0, 8.0],
        "parent": "left + child + right",
        "global": "left + child + right + ambient",
    },
    "rejected_overlap": "child (1, 3) + parent (0, 4) double-counts the child and is not accounted",
    "coordinate_fieldwork": "separate channel, same windows and the same RK weights; frozen geometry evaluates to the product with Qdot = Ldot = 0",
    "probability_flux": "occupation-weighted ψ† σ² ψ; the positive part of a net window flux is a gross lower bound and is not the accounted flux",
    "probability_gross": "RK-integrated positive parts at each endpoint are lower bounds on gross counterpropagating probability, not a resolved chiral traffic count",
    "piston_work_on_measurement_cut": False,
    "dx_divisions": 1,
    "multiplicity": "M = 4κ inside the nodal forces once",
    "clocks": "trapezoid of rQ at x = 1, 2, 3 per accepted step; RK stage samples are a comparison",
    "observation_cadence_0_05_used": False,
}


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def jsonable(value):
    if isinstance(value, np.ndarray):
        raise TypeError("balance records stay scalar")
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, (float, int, str, bool)) or value is None:
        return value
    raise TypeError(f"balance record cannot hold {type(value).__name__}")


def producer_hashes(lab=LAB):
    lab = Path(lab)
    return {relative: sha256_file(lab / relative) for relative in PRODUCERS}


def sealed_input_hashes(lab=LAB):
    """Hash frozen handoff bytes. Does not parse trajectory arrays."""
    lab = Path(lab)
    rows = []
    for relative in SEALED_INPUTS:
        path = lab / relative
        if not path.is_file():
            raise FileNotFoundError("missing sealed input " + relative)
        stat = path.stat()
        rows.append({
            "path": relative,
            "sha256": sha256_file(path),
            "bytes": int(stat.st_size),
            "writable": bool(stat.st_mode & 0o222),
            "parsed": False,
        })
    return rows


def sigma2_chiral_summand(phi0, phi1, occupations):
    """Occupation-weighted ``ψ† σ² ψ``.

    ``σ² = [[0, -i], [i, 0]]``, so the column density is
    ``-i φ0* φ1 + i φ1* φ0 = 2 Im(φ0* φ1)``. This is the signed chiral
    probability current sample. It is not the positive part of a net flux.
    """
    phi0 = np.asarray(phi0)
    phi1 = np.asarray(phi1)
    if phi0.shape != phi1.shape:
        raise ValueError("chiral flux needs matching spinor blocks")
    chiral = (-1j * np.conjugate(phi0) * phi1 + 1j * np.conjugate(phi1) * phi0).real
    weights = np.asarray(occupations, dtype=float)
    if chiral.ndim != 2 or weights.shape != (chiral.shape[1],):
        raise ValueError("chiral flux expects nodal columns and one occupation weight each")
    summed = np.sum(chiral * weights[None, :], axis=1)
    if not np.isfinite(summed).all():
        raise ValueError("chiral probability flux is not finite")
    return summed


def acceptance_steps(t0, t1, step_cap):
    """Forward steps of ``step_cap``, with a final remainder that lands on ``t1``."""
    start = float(t0)
    target = float(t1)
    cap = float(step_cap)
    if not all(math.isfinite(value) for value in (start, target, cap)):
        raise ValueError("times and the step cap must be finite")
    if cap <= 0.0:
        raise ValueError("step cap must be positive")
    if target <= start:
        raise ValueError("the ledger segment must advance")
    steps = []
    mark = start
    while mark < target - 1e-15:
        dt = min(cap, target - mark)
        if dt <= 0.0 or not math.isfinite(dt):
            raise RuntimeError("step schedule produced a nonpositive step")
        steps.append(float(dt))
        mark = mark + dt
        if abs(mark - target) <= 1e-12:
            mark = target
        if len(steps) > 1_000_000:
            raise RuntimeError("step schedule did not reach the endpoint")
    if abs(mark - target) > 1e-8:
        raise RuntimeError("step schedule missed the endpoint")
    return steps


def _require_carrier(grid):
    length = float(grid.length)
    spacing = float(grid.dx_q)
    if abs(length - CARRIER_LENGTH) > 1e-12 or spacing <= 0.0:
        raise ValueError("balance windows are the period-8 carrier")
    if getattr(grid, "gauge", "conformal") != "conformal":
        raise ValueError("balance ledger uses the conformal lapse chain")
    return spacing


def _require_clocks(pair):
    locations = tuple(float(value) for value in pair.clock_locations)
    if locations != CLOCK_LOCATIONS:
        raise ValueError("proper clocks stay at x = 1, 2, 3")


def _integrate_summand(grid, summand, interval):
    spacing = float(grid.dx_q)
    return float(model.interval_integral(grid, np.asarray(summand, dtype=float) / spacing, interval))


def _endpoint_flux(grid, summand, interval):
    spacing = float(grid.dx_q)
    ends = model._periodic_values(grid, np.asarray(summand, dtype=float), (float(interval[0]), float(interval[1])))
    return float(ends[0] / spacing), float(ends[1] / spacing), float((ends[0] - ends[1]) / spacing)


def _derivative_window(grid, summand, interval):
    spacing = float(grid.dx_q)
    divergence = np.asarray(grid.derivative @ np.asarray(summand, dtype=float), dtype=float)
    return float(model.interval_integral(grid, divergence / spacing, interval))


def _rk(samples, dt):
    first, second, third, fourth = samples
    return float(dt / 6 * (first + 2 * second + 2 * third + fourth))


def _blank_rates():
    return {key: 0.0 for key in RATE_KEYS}


def _sum_rates(parts):
    return {key: float(sum(part[key] for part in parts)) for key in RATE_KEYS}


def _window_rates(grid, fields, interval):
    """One fixed window. Reynolds is evaluated at zero boundary velocity."""
    piece = regions.moving_normal_energy_rate(
        grid, fields["energy"], fields["flux"], fields["pressure"], fields["lapse"], fields["slope"],
        interval, (0.0, 0.0),
    )
    derivative_flux = _derivative_window(grid, fields["flux"], interval)
    collocation = float(piece["fixed_boundary_flux"])
    nyquist_gap = float(piece["endpoint_identity_gap"])
    gap_from_window = float(-derivative_flux - collocation)
    if abs(nyquist_gap - gap_from_window) > 1e-8:
        raise RuntimeError("Nyquist gap does not match the derivative window")
    probability_left, probability_right, probability_net = _endpoint_flux(grid, fields["probability_flux"], interval)
    probability_derivative = _derivative_window(grid, fields["probability_flux"], interval)
    probability_slope = _integrate_summand(grid, fields.get("probability_slope", np.zeros(grid.nq)), interval)
    proper_left, proper_right, _proper_net = _endpoint_flux(grid, fields["flux"], interval)
    positive_part = float(max(collocation, 0.0))
    gross = float(abs(proper_left) + abs(proper_right))
    if positive_part > gross + 1e-9:
        raise RuntimeError("positive part of the net flux exceeded the gross cut magnitude")
    rates = {
        "fixed_boundary_flux": collocation,
        "derivative_window_flux": derivative_flux,
        "nyquist_gap": nyquist_gap,
        "pointwise_projection_defect": float(piece["pointwise_projection_defect"]),
        "finite_projection_defect": float(piece["finite_projection_defect"]),
        "pressure": float(piece["pressure"]),
        "lapse": float(piece["lapse"]),
        "reynolds": float(piece["reynolds"]),
        "slope_integral": float(piece["slope_integral"]),
        "balance_rate": float(piece["rate"]),
        "coordinate_metric_work": _integrate_summand(grid, fields["coordinate"], interval),
        "probability_boundary_flux": probability_net,
        "probability_derivative_window_flux": probability_derivative,
        "probability_nyquist_gap": float(-probability_derivative - probability_net),
        "probability_gross_in_lower_bound": float(max(probability_left, 0.) + max(-probability_right, 0.)),
        "probability_gross_out_lower_bound": float(max(-probability_left, 0.) + max(probability_right, 0.)),
        "probability_slope_integral": probability_slope,
        "probability_projection_defect": float(probability_slope - probability_net),
        "left_outward_proper": float(-proper_left),
        "right_outward_proper": float(proper_right),
        "left_inward_proper": float(proper_left),
        "right_inward_proper": float(-proper_right),
        "left_outward_probability": float(-probability_left),
        "right_outward_probability": float(probability_right),
        "left_inward_probability": float(probability_left),
        "right_inward_probability": float(-probability_right),
        "positive_part_of_net_flux": positive_part,
        "gross_cut_magnitude": gross,
        "instantaneous_closure_residual": float(piece["closure_residual"]),
    }
    stocks = {
        "normal_energy": _integrate_summand(grid, fields["energy"], interval),
        "probability": _integrate_summand(grid, fields["probability"], interval),
    }
    return {
        "interval": [float(interval[0]), float(interval[1])],
        "kind": "fixed_coordinate_window",
        "boundary_velocity": [0.0, 0.0],
        "rates": rates,
        "stocks": stocks,
        "surface_work": float(piece["surface_work"]),
        "piston_work_included": False,
        "physical_wall": False,
        "dx_divisions": 1,
        "used_in_closure": True,
        "positive_part_used_as_accounted_flux": False,
        "positive_part_net_flux_is_gross_lower_bound": True,
        "regions_closure_residual": float(piece["closure_residual"]),
        "positive_part_of_net_flux": positive_part,
    }


def partition_rates(grid, fields):
    """Disjoint windows, their derived parent and global sums, and the rejected overlap."""
    _require_carrier(grid)
    for name in ("energy", "flux", "pressure", "lapse", "slope", "coordinate", "probability", "probability_flux"):
        array = np.asarray(fields[name], dtype=float)
        if array.shape != (grid.nq,) or not np.isfinite(array).all():
            raise ValueError(f"{name} must be a finite nodal summand")
    windows = {name: _window_rates(grid, fields, interval) for name, interval in WINDOW_SPECS}
    direct_parent = _window_rates(grid, fields, DIRECT_PARENT_INTERVAL)
    direct_parent["kind"] = "overlapping_parent_window"
    direct_parent["used_in_closure"] = False
    derived = _sum_rates([windows[name]["rates"] for name in PARENT_PARTS])
    overlap_gap = float(max(
        abs(derived[key] - direct_parent["rates"][key])
        for key in ("pressure", "lapse", "fixed_boundary_flux", "coordinate_metric_work", "balance_rate")
    ))
    return {
        "windows": windows,
        "direct_parent": direct_parent,
        "derived_parent_matches_direct_gap": overlap_gap,
        "cuts": _cut_rows(grid, fields["flux"], fields["probability_flux"]),
        "owned_fieldwork_power": float(fields["owned_fieldwork_power"]) if "owned_fieldwork_power" in fields else float(np.sum(fields["coordinate"])),
        "window_fieldwork_sum": float(sum(windows[name]["rates"]["coordinate_metric_work"] for name in DISJOINT_WINDOWS)),
    }


def _cut_rows(grid, proper_summand, chiral_summand):
    spacing = float(grid.dx_q)
    positions = [float(cut["x"]) for cut in RESOLVED_CUTS] + [float(grid.length)]
    proper = model._periodic_values(grid, np.asarray(proper_summand, dtype=float), positions)
    chiral = model._periodic_values(grid, np.asarray(chiral_summand, dtype=float), positions)
    branch = float(proper[-1] / spacing)
    origin = float(proper[0] / spacing)
    rows = []
    for index, spec in enumerate(RESOLVED_CUTS):
        proper_flux = float(proper[index] / spacing)
        chiral_flux = float(chiral[index] / spacing)
        left_out = proper_flux
        right_out = -proper_flux
        left_out_p = chiral_flux
        right_out_p = -chiral_flux
        rows.append({
            "x": float(spec["x"]),
            "left_window": spec["left_window"],
            "right_window": spec["right_window"],
            "resolved_by": "fourier_interpolant",
            "proper_normal_flux": proper_flux,
            "sigma2_chiral_probability_flux": chiral_flux,
            "left_window_outward_proper": left_out,
            "left_window_inward_proper": -left_out,
            "right_window_outward_proper": right_out,
            "right_window_inward_proper": -right_out,
            "left_window_outward_probability": left_out_p,
            "left_window_inward_probability": -left_out_p,
            "right_window_outward_probability": right_out_p,
            "right_window_inward_probability": -right_out_p,
            "outward_pair_sum_proper": float(left_out + right_out),
            "outward_pair_sum_probability": float(left_out_p + right_out_p),
        })
    rows[0]["branch_flux_minus_origin_proper"] = float(branch - origin)
    return rows


def _clock_sample(grid, fine):
    radial = np.asarray(fine.r, dtype=float) * np.asarray(fine.Q, dtype=float)
    sampled = model._periodic_values(grid, radial, CLOCK_LOCATIONS)
    return [float(value) for value in sampled]


def _protocol_clock_rates(pair, state):
    return [float(value) for value in model.metrics(pair, state)["clock_rates"]]


def _fields_from_bundle(pair, bundle, control_mode):
    grid = pair.grid
    fine = bundle["fine_state"]
    system = bundle["fine_system"]
    source = bundle["source"]
    if source.get("multiplicity_applied_once") is not True:
        raise ValueError("source multiplicity was not applied once")
    if int(system.multiplicity) != 4 * int(system.kappa):
        raise ValueError("M=4κ multiplicity drifted")
    if float(system.dx) != float(grid.dx_q):
        raise ValueError("fine spacing and carrier spacing differ; refusing a second dx")
    if system.derivative is not grid.derivative:
        raise ValueError("fine derivative is not the carrier Fourier derivative")
    lifted = observer._lifted_rate(grid, bundle["nodal_rate"], bundle)
    if control_mode == "frozen_geometry":
        lifted = regions._freeze_geometry(lifted)
        q_dot = np.asarray(lifted.Q, dtype=float)
        lapse_dot = q_dot
        rate_power = float(lifted.fieldwork_power)
    else:
        q_dot = np.asarray(bundle["lifted_Q"], dtype=float)
        lapse_dot = np.asarray(bundle["lapse_dot"], dtype=float)
        rate_power = float(bundle["nodal_rate"].fieldwork_power)
    coordinate = (
        np.asarray(source["force_Q"], dtype=float) * q_dot
        + np.asarray(source["force_L"], dtype=float) * lapse_dot
    )
    # The full source is already evaluated by this RK stage. Reuse it rather
    # than applying the source operator again inside matter_ledger.
    stresses = regional.observer_stresses(system, fine, source)
    ledger = {"stresses": stresses, "normal_energy_nodal": stresses["normal_energy_nodal"]}
    terms = regional.proper_balance_terms(system, fine, lifted, ledger)
    slope = regions.analytic_normal_energy_slope(system, fine, lifted, source)
    weights = np.asarray(system.occupations, dtype=float)
    mass = np.sum((np.abs(fine.phi0) ** 2 + np.abs(fine.phi1) ** 2) * weights[None, :], axis=1)
    chiral = sigma2_chiral_summand(fine.phi0, fine.phi1, weights)
    fields = {
        "energy": np.asarray(ledger["normal_energy_nodal"], dtype=float),
        "flux": np.asarray(terms["flux_nodal"], dtype=float),
        "pressure": np.asarray(terms["proper_pressure_work"], dtype=float),
        "lapse": np.asarray(terms["momentum_lapse_work"], dtype=float),
        "slope": np.asarray(slope, dtype=float),
        "coordinate": np.asarray(coordinate, dtype=float),
        "probability": np.asarray(mass, dtype=float),
        "probability_flux": np.asarray(chiral, dtype=float),
        "probability_slope": 2. * np.sum(weights[None, :] * np.real(
            np.conjugate(fine.phi0) * lifted.phi0 + np.conjugate(fine.phi1) * lifted.phi1), axis=1),
        "owned_fieldwork_power": float(np.sum(coordinate)),
    }
    return fields, fine, int(system.multiplicity), rate_power


def _stage_sample(pair, state, control_mode):
    rate, bundle = model.rates(pair, state, return_bundle=True)
    if control_mode == "frozen_geometry":
        rate = episode._freeze_geometry_rate(rate)
    fields, fine, multiplicity, rate_power = _fields_from_bundle(pair, bundle, control_mode)
    partition = partition_rates(pair.grid, fields)
    clocks = _clock_sample(pair.grid, fine)
    sample = {
        "partition": partition,
        "clocks": clocks,
        "multiplicity": multiplicity,
        "owned_fieldwork_power": float(fields["owned_fieldwork_power"]),
        "rate_fieldwork_power": float(rate_power),
        "control_rate_fieldwork_power": float(rate.fieldwork_power),
    }
    return rate, sample


def _nyquist_symbol(grid):
    nyquist = np.cos(np.pi * np.arange(grid.nq))
    return float(np.max(np.abs(grid.derivative @ nyquist)))


def _integrate_named(samples, dt, key):
    return _rk([sample[key] for sample in samples], dt)


def _integrated_window(parts, dt, *, kind, interval, used_in_closure):
    integrated = {key: _rk([part["rates"][key] for part in parts], dt) for key in RATE_KEYS}
    window = {
        "interval": list(interval),
        "kind": kind,
        "boundary_velocity": [0.0, 0.0],
        "used_in_closure": used_in_closure,
        "piston_work_included": False,
        "physical_wall": False,
        "surface_work": 0.0,
        "dx_divisions": 1,
        **integrated,
    }
    return _declare_net_bound(window)


def _attach_change(window, start_stock, end_stock):
    window["normal_energy_start"] = float(start_stock["normal_energy"])
    window["normal_energy_end"] = float(end_stock["normal_energy"])
    window["probability_start"] = float(start_stock["probability"])
    window["probability_end"] = float(end_stock["probability"])
    window["measured_change"] = float(end_stock["normal_energy"] - start_stock["normal_energy"])
    window["accounted_change"] = float(
        window["fixed_boundary_flux"] + window["pressure"] + window["lapse"]
        + window["finite_projection_defect"] + window["reynolds"]
    )
    window["quadrature_closure_error"] = float(window["measured_change"] - window["accounted_change"])
    window["probability_measured_change"] = float(end_stock["probability"] - start_stock["probability"])
    window["probability_accounted_change"] = float(window["probability_boundary_flux"] + window["probability_projection_defect"])
    window["probability_quadrature_error"] = float(window["probability_measured_change"] - window["probability_accounted_change"])
    return window


def _declare_net_bound(window):
    """The positive part of the signed net is a lower bound, not the accounted flux."""
    net = float(window["fixed_boundary_flux"])
    gross = float(window["gross_cut_magnitude"])
    positive = float(max(net, 0.0))
    if positive > gross + 1e-8:
        raise RuntimeError("positive part of the net flux exceeded the integrated gross cut magnitude")
    window["positive_part_of_net_flux"] = positive
    window["positive_part_used_as_accounted_flux"] = False
    window["positive_part_net_flux_is_gross_lower_bound"] = True
    return window


def _physical_scale(window):
    return max(
        abs(window["measured_change"]),
        abs(window["accounted_change"]),
        abs(window["fixed_boundary_flux"]),
        abs(window["pressure"]),
        abs(window["lapse"]),
        abs(window["finite_projection_defect"]),
        abs(window["coordinate_metric_work"]),
    )


def _closure_view(windows):
    rows = {}
    for name, window in windows.items():
        scale = _physical_scale(window)
        error = abs(window["quadrature_closure_error"])
        rows[name] = {
            "quadrature_closure_error": window["quadrature_closure_error"],
            "physical_scale": scale,
            "closure_error_over_physical": None if scale == 0.0 else float(error / scale),
            "measured_change": window["measured_change"],
            "accounted_change": window["accounted_change"],
            "finite_projection_defect": window["finite_projection_defect"],
            "nyquist_gap": window["nyquist_gap"],
            "pointwise_projection_defect": window["pointwise_projection_defect"],
        }
    return rows


def _states_bitwise(left, right):
    return all(np.array_equal(getattr(left, name), getattr(right, name)) for name in model.STATE_NAMES)


def _chart_check(pair, state):
    fine = galerkin.prolong_state(pair.grid, model.reconstruct_state(pair, state))
    active = galerkin.active_fine_system(pair.grid, fine)
    reason = coupling.chart_failure(active, fine)
    if reason is not None:
        raise coupling.PositiveChartExit(reason, np.nan, state.copy())
    return fine


def _source_token(pair):
    return episode.array_sha256(np.vstack((np.array(pair.source_phi0, copy=True), np.array(pair.source_phi1, copy=True))))


def energy_stocks(pair, state):
    """Window stocks of ``F_L/r`` and probability. No rate and no step."""
    _require_carrier(pair.grid)
    nodal = model.reconstruct_state(pair, state)
    fine = galerkin.prolong_state(pair.grid, nodal)
    system = galerkin.active_fine_system(pair.grid, fine)
    reason = coupling.chart_failure(system, fine)
    if reason is not None:
        raise coupling.PositiveChartExit(reason, np.nan, state.copy())
    source = coupling.source_from_columns(system, fine)
    weights = np.asarray(system.occupations, dtype=float)
    mass = np.sum((np.abs(fine.phi0) ** 2 + np.abs(fine.phi1) ** 2) * weights[None, :], axis=1)
    energy = np.asarray(source["force_L"] / fine.r, dtype=float)
    stocks = {}
    for name, interval in WINDOW_SPECS:
        stocks[name] = {
            "normal_energy": _integrate_summand(pair.grid, energy, interval),
            "probability": _integrate_summand(pair.grid, mass, interval),
        }
    stocks["direct_parent"] = {
        "normal_energy": _integrate_summand(pair.grid, energy, DIRECT_PARENT_INTERVAL),
        "probability": _integrate_summand(pair.grid, mass, DIRECT_PARENT_INTERVAL),
    }
    return stocks, _clock_sample(pair.grid, fine)


def instantaneous_balance(pair, state, *, bundle=None, nodal_rate=None, control_mode="coupled"):
    """One stage sample. Does not advance the state."""
    mode = episode.normalize_control_mode(control_mode)
    _require_carrier(pair.grid)
    _require_clocks(pair)
    before = episode.state_sha256(state)
    source_before = _source_token(pair)
    if bundle is None and nodal_rate is None:
        _rate, sample = _stage_sample(pair, state, mode)
    elif bundle is None or nodal_rate is None:
        raise ValueError("bundle and nodal_rate are supplied together")
    else:
        bundle = dict(bundle, nodal_rate=nodal_rate)
        fields, fine, multiplicity, rate_power = _fields_from_bundle(pair, bundle, mode)
        sample = {
            "partition": partition_rates(pair.grid, fields),
            "clocks": _clock_sample(pair.grid, fine),
            "multiplicity": multiplicity,
            "owned_fieldwork_power": float(fields["owned_fieldwork_power"]),
            "rate_fieldwork_power": float(rate_power),
            "control_rate_fieldwork_power": 0.0 if mode == "frozen_geometry" else float(rate_power),
        }
    if episode.state_sha256(state) != before or _source_token(pair) != source_before:
        raise RuntimeError("balance sample mutated the state or the source columns")
    partition = sample["partition"]
    windows = {}
    for name, window in partition["windows"].items():
        public = dict(window["rates"])
        public.update(window["stocks"])
        public.update({
            "interval": window["interval"],
            "kind": window["kind"],
            "boundary_velocity": window["boundary_velocity"],
            "surface_work": window["surface_work"],
            "piston_work_included": False,
            "physical_wall": False,
            "dx_divisions": 1,
            "used_in_closure": True,
            "positive_part_used_as_accounted_flux": False,
            "positive_part_net_flux_is_gross_lower_bound": True,
        })
        windows[name] = _declare_net_bound(public)
    parent = _sum_rates([partition["windows"][name]["rates"] for name in PARENT_PARTS])
    parent_stocks = {
        "normal_energy": float(sum(partition["windows"][name]["stocks"]["normal_energy"] for name in PARENT_PARTS)),
        "probability": float(sum(partition["windows"][name]["stocks"]["probability"] for name in PARENT_PARTS)),
    }
    global_rates = _sum_rates([partition["windows"][name]["rates"] for name in DISJOINT_WINDOWS])
    global_stocks = {
        "normal_energy": float(sum(partition["windows"][name]["stocks"]["normal_energy"] for name in DISJOINT_WINDOWS)),
        "probability": float(sum(partition["windows"][name]["stocks"]["probability"] for name in DISJOINT_WINDOWS)),
    }
    windows["parent"] = _declare_net_bound(dict(
        parent, **parent_stocks, interval=[0.0, 4.0], kind="derived_sum",
        boundary_velocity=[0.0, 0.0], used_in_closure=True,
        piston_work_included=False, physical_wall=False, surface_work=0.0, dx_divisions=1,
    ))
    windows["global"] = _declare_net_bound(dict(
        global_rates, **global_stocks, interval=[0.0, CARRIER_LENGTH], kind="derived_sum",
        boundary_velocity=[0.0, 0.0], used_in_closure=True,
        piston_work_included=False, physical_wall=False, surface_work=0.0, dx_divisions=1,
    ))
    child = partition["windows"]["child"]["rates"]
    direct = partition["direct_parent"]["rates"]
    return {
        "api_version": API_VERSION,
        "control_mode": mode,
        "identity": dict(CLOSURE_IDENTITY),
        "windows": windows,
        "cuts": partition["cuts"],
        "direct_parent_window": dict(direct, **partition["direct_parent"]["stocks"],
                                     interval=[0.0, 4.0], kind="overlapping_parent_window",
                                     used_in_closure=False, piston_work_included=False,
                                     physical_wall=False, surface_work=0.0, dx_divisions=1,
                                     positive_part_used_as_accounted_flux=False,
                                     positive_part_net_flux_is_gross_lower_bound=True),
        "rejected_overlap": {
            "definition": "child (1, 3) + parent (0, 4)",
            "pressure": float(child["pressure"] + direct["pressure"]),
            "lapse": float(child["lapse"] + direct["lapse"]),
            "derived_parent_pressure": float(parent["pressure"]),
            "derived_parent_lapse": float(parent["lapse"]),
            "double_count_pressure": float(child["pressure"] + direct["pressure"] - parent["pressure"]),
            "double_count_lapse": float(child["lapse"] + direct["lapse"] - parent["lapse"]),
            "used_in_closure": False,
        },
        "coordinate_fieldwork": {
            "owned_power": sample["owned_fieldwork_power"],
            "rate_power": sample["rate_fieldwork_power"],
            "control_rate_power": sample["control_rate_fieldwork_power"],
            "window_sum": partition["window_fieldwork_sum"],
            "window_minus_owned": float(partition["window_fieldwork_sum"] - sample["owned_fieldwork_power"]),
        },
        "internal_cut_cancellation_residual": float(global_rates["fixed_boundary_flux"]),
        "nyquist_symbol_max": _nyquist_symbol(pair.grid),
        "clock_rates": sample["clocks"],
        "clock_locations": list(CLOCK_LOCATIONS),
        "multiplicity": sample["multiplicity"],
        "multiplicity_applied_once": True,
        "dx_divisions": 1,
        "piston_work_on_measurement_cut": False,
        "gradient_cut_integrated": False,
        "observation_cadence_used": False,
        "child_normal_closure_claimed": False,
        "global_normal_closure_claimed": False,
        "continuum_certified": False,
        "renewal_asserted": False,
    }


def rk4_stage_ledger(pair, state, dt, *, control_mode="coupled"):
    """One owned RK4 step and the RK4-weighted balance on its four stages."""
    mode = episode.normalize_control_mode(control_mode)
    step = float(dt)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("a positive finite timestep is required")
    _require_carrier(pair.grid)
    _require_clocks(pair)
    before = episode.state_sha256(state)
    source_before = _source_token(pair)
    protocol_start = _protocol_clock_rates(pair, state)
    k1, s1 = _stage_sample(pair, state, mode)
    k2, s2 = _stage_sample(pair, model._combine(state, k1, step / 2), mode)
    k3, s3 = _stage_sample(pair, model._combine(state, k2, step / 2), mode)
    k4, s4 = _stage_sample(pair, model._combine(state, k3, step), mode)
    updated = model.NestedState(*(getattr(state, name) + step / 6 * (
        getattr(k1, name) + 2 * getattr(k2, name) + 2 * getattr(k3, name) + getattr(k4, name))
        for name in model.STATE_NAMES))
    _chart_check(pair, updated)
    if episode.state_sha256(state) != before or _source_token(pair) != source_before:
        raise RuntimeError("RK4 ledger mutated the input state or the source columns")
    stages = (s1, s2, s3, s4)
    start_stocks = {name: dict(s1["partition"]["windows"][name]["stocks"]) for name in DISJOINT_WINDOWS}
    start_stocks["direct_parent"] = dict(s1["partition"]["direct_parent"]["stocks"])
    end_stocks, end_clocks = energy_stocks(pair, updated)
    protocol_end = _protocol_clock_rates(pair, updated)
    windows = {}
    for name, _interval in WINDOW_SPECS:
        parts = [stage["partition"]["windows"][name] for stage in stages]
        window = _integrated_window(parts, step, kind="fixed_coordinate_window",
                                    interval=parts[0]["interval"], used_in_closure=True)
        windows[name] = _attach_change(window, start_stocks[name], end_stocks[name])
    parent = _sum_rates([windows[name] for name in PARENT_PARTS])
    parent.update({
        "interval": [0.0, 4.0], "kind": "derived_sum", "boundary_velocity": [0.0, 0.0],
        "used_in_closure": True, "piston_work_included": False, "physical_wall": False,
        "surface_work": 0.0, "dx_divisions": 1, "positive_part_used_as_accounted_flux": False,
        "positive_part_net_flux_is_gross_lower_bound": True,
    })
    parent_start = {
        "normal_energy": float(sum(windows[name]["normal_energy_start"] for name in PARENT_PARTS)),
        "probability": float(sum(windows[name]["probability_start"] for name in PARENT_PARTS)),
    }
    parent_end = {
        "normal_energy": float(sum(end_stocks[name]["normal_energy"] for name in PARENT_PARTS)),
        "probability": float(sum(end_stocks[name]["probability"] for name in PARENT_PARTS)),
    }
    windows["parent"] = _declare_net_bound(_attach_change(parent, parent_start, parent_end))
    global_window = _sum_rates([windows[name] for name in DISJOINT_WINDOWS])
    global_window.update({
        "interval": [0.0, CARRIER_LENGTH], "kind": "derived_sum", "boundary_velocity": [0.0, 0.0],
        "used_in_closure": True, "piston_work_included": False, "physical_wall": False,
        "surface_work": 0.0, "dx_divisions": 1, "positive_part_used_as_accounted_flux": False,
        "positive_part_net_flux_is_gross_lower_bound": True,
    })
    global_start = {
        "normal_energy": float(sum(windows[name]["normal_energy_start"] for name in DISJOINT_WINDOWS)),
        "probability": float(sum(windows[name]["probability_start"] for name in DISJOINT_WINDOWS)),
    }
    global_end = {
        "normal_energy": float(sum(end_stocks[name]["normal_energy"] for name in DISJOINT_WINDOWS)),
        "probability": float(sum(end_stocks[name]["probability"] for name in DISJOINT_WINDOWS)),
    }
    windows["global"] = _declare_net_bound(_attach_change(global_window, global_start, global_end))
    direct_parts = [stage["partition"]["direct_parent"] for stage in stages]
    direct = _integrated_window(direct_parts, step, kind="overlapping_parent_window",
                                interval=[0.0, 4.0], used_in_closure=False)
    direct = _attach_change(direct, start_stocks["direct_parent"], end_stocks["direct_parent"])
    child = windows["child"]
    rejected_pressure = float(child["pressure"] + direct["pressure"])
    rejected_lapse = float(child["lapse"] + direct["lapse"])
    cuts = []
    for index, spec in enumerate(RESOLVED_CUTS):
        row = {"x": spec["x"], "left_window": spec["left_window"], "right_window": spec["right_window"],
               "resolved_by": "fourier_interpolant"}
        for key in CUT_KEYS:
            row[key] = _rk([stage["partition"]["cuts"][index][key] for stage in stages], step)
        cuts.append(row)
    trap = [0.5 * step * (protocol_start[index] + protocol_end[index]) for index in range(3)]
    rk_clocks = [_rk([stage["clocks"][index] for stage in stages], step) for index in range(3)]
    sampler_gap = float(max(abs(s1["clocks"][index] - protocol_start[index]) for index in range(3)))
    fieldwork = _rk([stage["owned_fieldwork_power"] for stage in stages], step)
    window_fieldwork = float(sum(windows[name]["coordinate_metric_work"] for name in DISJOINT_WINDOWS))
    record = {
        "api_version": API_VERSION,
        "control_mode": mode,
        "dt": step,
        "quadrature": "rk4_stage_weights",
        "rk_weight": "dt/6 * (f1 + 2*f2 + 2*f3 + f4)",
        "observation_cadence_used": False,
        "identity": dict(CLOSURE_IDENTITY),
        "windows": windows,
        "cuts": cuts,
        "direct_parent_window": direct,
        "rejected_overlap": {
            "definition": "child (1, 3) + parent (0, 4)",
            "pressure": rejected_pressure,
            "lapse": rejected_lapse,
            "derived_parent_pressure": float(windows["parent"]["pressure"]),
            "derived_parent_lapse": float(windows["parent"]["lapse"]),
            "double_count_pressure": float(rejected_pressure - windows["parent"]["pressure"]),
            "double_count_lapse": float(rejected_lapse - windows["parent"]["lapse"]),
            "used_in_closure": False,
        },
        "coordinate_fieldwork": {
            "owned_power": fieldwork,
            "window_sum": window_fieldwork,
            "window_minus_owned": float(window_fieldwork - fieldwork),
            "frozen_geometry": mode == "frozen_geometry",
        },
        "internal_cut_cancellation_residual": float(windows["global"]["fixed_boundary_flux"]),
        "nyquist_symbol_max": _nyquist_symbol(pair.grid),
        "closure": _closure_view(windows),
        "clocks": {
            "locations": list(CLOCK_LOCATIONS),
            "protocol": "trapezoid_coordinate_time_per_accepted_step",
            "rk_clock_role": "comparison",
            "start_rates": protocol_start,
            "end_rates": protocol_end,
            "stage_rates": [list(stage["clocks"]) for stage in stages],
            "trapezoid_increment": trap,
            "rk_increment": rk_clocks,
            "rk_minus_trapezoid": [float(rk_clocks[index] - trap[index]) for index in range(3)],
            "stage_sampler_minus_protocol_start_max": sampler_gap,
            "end_fine_minus_protocol_max": float(max(abs(end_clocks[index] - protocol_end[index]) for index in range(3))),
        },
        "stage_fieldwork_power": [float(stage["owned_fieldwork_power"]) for stage in stages],
        "stage_rate_fieldwork_power": [float(stage["rate_fieldwork_power"]) for stage in stages],
        "multiplicity": int(s1["multiplicity"]),
        "multiplicity_applied_once": True,
        "dx_divisions": 1,
        "piston_work_on_measurement_cut": False,
        "gradient_cut_integrated": False,
        "source_token": source_before,
        "child_normal_closure_claimed": False,
        "global_normal_closure_claimed": False,
        "continuum_certified": False,
        "renewal_asserted": False,
    }
    return updated, record


def _empty_accumulator():
    windows = {}
    for name, interval in list(WINDOW_SPECS) + [("parent", (0.0, 4.0)), ("global", (0.0, CARRIER_LENGTH))]:
        windows[name] = {key: 0.0 for key in RATE_KEYS}
        windows[name].update({
            "interval": list(interval),
            "measured_change": 0.0,
            "accounted_change": 0.0,
            "quadrature_closure_error": 0.0,
        })
    return {
        "windows": windows,
        "cuts": [{key: 0.0 for key in CUT_KEYS} | {"x": spec["x"], "left_window": spec["left_window"],
                 "right_window": spec["right_window"]} for spec in RESOLVED_CUTS],
        "rejected_pressure": 0.0,
        "rejected_lapse": 0.0,
        "owned_fieldwork": 0.0,
        "window_fieldwork": 0.0,
        "quadrature_error_sum": {name: 0.0 for name in list(DISJOINT_WINDOWS) + ["parent", "global"]},
    }


def _add_step(total, step):
    for name, window in step["windows"].items():
        bucket = total["windows"][name]
        for key in RATE_KEYS:
            bucket[key] += float(window[key])
        bucket["measured_change"] += float(window["measured_change"])
        bucket["accounted_change"] += float(window["accounted_change"])
        bucket["quadrature_closure_error"] += float(window["quadrature_closure_error"])
    for index, cut in enumerate(step["cuts"]):
        for key in CUT_KEYS:
            total["cuts"][index][key] += float(cut[key])
    total["rejected_pressure"] += float(step["rejected_overlap"]["pressure"])
    total["rejected_lapse"] += float(step["rejected_overlap"]["lapse"])
    total["owned_fieldwork"] += float(step["coordinate_fieldwork"]["owned_power"])
    total["window_fieldwork"] += float(step["coordinate_fieldwork"]["window_sum"])
    return total


def integrate_segment(pair, state, t0, t1, step_cap, *, control_mode="coupled",
                      normal_clocks=None, verify_first_step=True,
                      owned_step_restriction=False, cpu_budget_seconds=CPU_BUDGET_SECONDS):
    """Advance ``t0 → t1`` and accumulate the stage ledger. The input state is copied."""
    mode = episode.normalize_control_mode(control_mode)
    import time
    acceptance_steps(t0, t1, step_cap)  # Validate the segment before mutation.
    steps = []
    started = time.process_time()
    mark = float(t0)
    current = state.copy()
    source_before = _source_token(pair)
    start_stocks, _start_clocks = energy_stocks(pair, current)
    protocol = np.zeros(3, dtype=float) if normal_clocks is None else np.array(normal_clocks, dtype=float)
    rk_clocks = np.array(protocol, dtype=float, copy=True)
    if protocol.shape != (3,):
        raise ValueError("three proper clocks are required")
    total = _empty_accumulator()
    bitwise = None
    max_step_ratio = None
    while mark < float(t1) - 1e-12:
        if time.process_time() - started >= float(cpu_budget_seconds):
            raise RuntimeError("six-hour balance CPU budget exhausted before the station")
        restricted = episode.step_restriction(pair, current, step_cap)[0] if owned_step_restriction else step_cap
        dt = min(float(restricted), float(t1) - mark)
        index = len(steps)
        steps.append(dt)
        updated, step = rk4_stage_ledger(pair, current, dt, control_mode=mode)
        if index == 0 and verify_first_step:
            reference_step = episode.frozen_geometry_step if mode == "frozen_geometry" else model.rk4_step
            reference = reference_step(pair, current.copy(), dt)
            bitwise = _states_bitwise(updated, reference)
            if not bitwise:
                raise RuntimeError("stage ledger state does not bit-match the owned RK4 step")
        for clock_index in range(3):
            protocol[clock_index] += float(step["clocks"]["trapezoid_increment"][clock_index])
            rk_clocks[clock_index] += float(step["clocks"]["rk_increment"][clock_index])
        for name, view in step["closure"].items():
            ratio = view["closure_error_over_physical"]
            if ratio is not None and (max_step_ratio is None or ratio > max_step_ratio):
                max_step_ratio = float(ratio)
        total = _add_step(total, step)
        current = updated
        mark += dt
    end_stocks, _end_clocks = energy_stocks(pair, current)
    if _source_token(pair) != source_before:
        raise RuntimeError("segment changed the frozen source columns")
    windows = {}
    for name in list(DISJOINT_WINDOWS) + ["parent", "global"]:
        bucket = total["windows"][name]
        if name in DISJOINT_WINDOWS:
            start = start_stocks[name]
            end = end_stocks[name]
        elif name == "parent":
            start = {
                "normal_energy": float(sum(start_stocks[part]["normal_energy"] for part in PARENT_PARTS)),
                "probability": float(sum(start_stocks[part]["probability"] for part in PARENT_PARTS)),
            }
            end = {
                "normal_energy": float(sum(end_stocks[part]["normal_energy"] for part in PARENT_PARTS)),
                "probability": float(sum(end_stocks[part]["probability"] for part in PARENT_PARTS)),
            }
        else:
            start = {
                "normal_energy": float(sum(start_stocks[part]["normal_energy"] for part in DISJOINT_WINDOWS)),
                "probability": float(sum(start_stocks[part]["probability"] for part in DISJOINT_WINDOWS)),
            }
            end = {
                "normal_energy": float(sum(end_stocks[part]["normal_energy"] for part in DISJOINT_WINDOWS)),
                "probability": float(sum(end_stocks[part]["probability"] for part in DISJOINT_WINDOWS)),
            }
        endpoint_change = float(end["normal_energy"] - start["normal_energy"])
        accounted = float(
            bucket["fixed_boundary_flux"] + bucket["pressure"] + bucket["lapse"]
            + bucket["finite_projection_defect"] + bucket["reynolds"]
        )
        windows[name] = dict(bucket)
        windows[name].update({
            "normal_energy_start": float(start["normal_energy"]),
            "normal_energy_end": float(end["normal_energy"]),
            "probability_start": float(start["probability"]),
            "probability_end": float(end["probability"]),
            "measured_change": endpoint_change,
            "endpoint_measured_change": endpoint_change,
            "accounted_change": accounted,
            "quadrature_closure_error": float(endpoint_change - accounted),
            "telescoped_step_measured_change": float(bucket["measured_change"]),
            "stock_telescoping_gap": float(endpoint_change - bucket["measured_change"]),
            "used_in_closure": True,
            "piston_work_included": False,
            "dx_divisions": 1,
            "kind": "fixed_coordinate_window" if name in DISJOINT_WINDOWS else "derived_sum",
        })
        windows[name] = _declare_net_bound(windows[name])
        windows[name]["probability_measured_change"] = float(end["probability"] - start["probability"])
        windows[name]["probability_accounted_change"] = float(bucket["probability_boundary_flux"] + bucket["probability_projection_defect"])
        windows[name]["probability_quadrature_error"] = float(windows[name]["probability_measured_change"] - windows[name]["probability_accounted_change"])
    closure = _closure_view(windows)
    record = {
        "api_version": API_VERSION,
        "control_mode": mode,
        "t0": float(t0),
        "t1": float(t1),
        "step_cap": float(step_cap),
        "steps": steps,
        "step_count": len(steps),
        "owned_step_restriction": bool(owned_step_restriction),
        "quadrature": "rk4_stage_weights",
        "observation_cadence_used": False,
        "identity": dict(CLOSURE_IDENTITY),
        "windows": windows,
        "cuts": total["cuts"],
        "rejected_overlap": {
            "definition": "child (1, 3) + parent (0, 4)",
            "pressure": float(total["rejected_pressure"]),
            "lapse": float(total["rejected_lapse"]),
            "derived_parent_pressure": float(windows["parent"]["pressure"]),
            "derived_parent_lapse": float(windows["parent"]["lapse"]),
            "double_count_pressure": float(total["rejected_pressure"] - windows["parent"]["pressure"]),
            "double_count_lapse": float(total["rejected_lapse"] - windows["parent"]["lapse"]),
            "used_in_closure": False,
        },
        "coordinate_fieldwork": {
            "owned_power": float(total["owned_fieldwork"]),
            "window_sum": float(total["window_fieldwork"]),
            "window_minus_owned": float(total["window_fieldwork"] - total["owned_fieldwork"]),
            "frozen_geometry": mode == "frozen_geometry",
        },
        "internal_cut_cancellation_residual": float(windows["global"]["fixed_boundary_flux"]),
        "closure": closure,
        "max_step_closure_error_over_physical": max_step_ratio,
        "clocks": {
            "locations": list(CLOCK_LOCATIONS),
            "protocol": "trapezoid_coordinate_time_per_accepted_step",
            "rk_clock_role": "comparison",
            "trapezoid": [float(value) for value in protocol],
            "rk": [float(value) for value in rk_clocks],
            "rk_minus_trapezoid": [float(rk_clocks[index] - protocol[index]) for index in range(3)],
        },
        "first_step_bitwise_owned_rk4": bitwise,
        "source_token": source_before,
        "multiplicity_applied_once": True,
        "dx_divisions": 1,
        "piston_work_on_measurement_cut": False,
        "gradient_cut_integrated": False,
        "child_normal_closure_claimed": False,
        "global_normal_closure_claimed": False,
        "continuum_certified": False,
        "renewal_asserted": False,
        "production_trajectory_stored": False,
    }
    return current, record


def _historical_step_cost(manifest):
    rows = {}
    for case in manifest.get("results", []):
        steps = int(case.get("steps") or 0)
        cpu = case.get("child_cpu_seconds")
        if steps <= 0 or cpu is None:
            continue
        rows[case["case_id"]] = {
            "child_cpu_seconds": float(cpu),
            "steps": steps,
            "seconds_per_step": float(cpu) / steps,
            "includes_observation_samples": True,
        }
    return rows


def production_cost(lab=LAB, *, target_time=TARGET_TIME):
    """Cap-based step counts and historical RK4 CPU to a chosen endpoint.

    Reads the episode manifest and the confirmation forecast. Does not parse
    npz payloads and does not take a step.
    """
    import json
    lab = Path(lab)
    manifest_path = lab / EPISODE_MANIFEST
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {"results": []}
    historical = _historical_step_cost(manifest)
    confirmation_path = lab / "results/development/nsc-nested-parent-child-confirmation-v2.json"
    pilot = None
    if confirmation_path.is_file():
        confirmation = json.loads(confirmation_path.read_text())
        pilot = (confirmation.get("forecasts") or {}).get("nf512_baseline_dt0.0005")
    span = float(target_time - HANDOFF_TIME)
    cases = []
    for nf, case_id, cap in (
        (256, "nf256_coupled_dt0.001", episode.STEP_CAP),
        (256, "nf256_coupled_dt0.0005", episode.MATCHED_STEP_CAP),
        (512, "nf512_baseline_dt0.0005", episode.MATCHED_STEP_CAP),
    ):
        steps = acceptance_steps(HANDOFF_TIME, target_time, cap)
        history = historical.get(case_id)
        if nf == 512 and pilot is not None:
            seconds_per_step = float(pilot["step_cpu_seconds"])
            source = "confirmation_v2_pilot_step"
        elif history is not None:
            seconds_per_step = float(history["seconds_per_step"])
            source = "episode_v1_child_cpu_per_accepted_step"
        else:
            seconds_per_step = None
            source = "unavailable"
        cases.append({
            "nf": nf,
            "historical_case_id": case_id,
            "step_cap": float(cap),
            "step_count": len(steps),
            "span": span,
            "seconds_per_historical_rk4_step": seconds_per_step,
            "estimated_seconds_at_historical_rate": None if seconds_per_step is None else float(len(steps) * seconds_per_step),
            "estimate_source": source,
            "extra_ode_stage_evaluations": 0,
            "endpoint_stock_source_evaluations_per_step": 1,
            "startup_bitwise_reference_steps_per_case": 1,
            "ledger_reuses_stage_rates": True,
            "measured_on_this_call": False,
        })
    two_h = [case for case in cases if case["nf"] == 256]
    spatial = [case for case in cases if case["step_cap"] == episode.MATCHED_STEP_CAP]
    def _sum(rows):
        if any(row["estimated_seconds_at_historical_rate"] is None for row in rows):
            return None
        return float(sum(row["estimated_seconds_at_historical_rate"] for row in rows))
    return {
        "handoff_time": HANDOFF_TIME,
        "target_time": float(target_time),
        "span": span,
        "new_trajectory_measured": False,
        "npz_parsed": False,
        "two_h_estimated_seconds": _sum(two_h),
        "spatial_matched_cap_estimated_seconds": _sum(spatial),
        "full_batch_estimated_seconds_at_historical_rate": _sum(cases),
        "historical_forecast_fraction_of_cpu_budget": None if _sum(cases) is None else _sum(cases) / CPU_BUDGET_SECONDS,
        "cases": cases,
        "cpu_budget_seconds": CPU_BUDGET_SECONDS,
        "chunk_limit_bytes": CHUNK_LIMIT_BYTES,
        "stage_ledger_overhead_measured": False,
        "step_counts_are_cap_based_lower_bounds": True,
        "note": (
            "The historical second is an accepted episode step, including its observation samples. "
            "The ledger uses those same four stage rates and does not add an ODE evaluation. "
            "Window integrals are extra arithmetic on that sample."
        ),
    }


def check_record(lab=LAB, *, target_time=TARGET_TIME):
    """Hash producers and sealed inputs. Writes nothing and takes no step."""
    return {
        "schema": CHECK_SCHEMA,
        "api_version": API_VERSION,
        "identity": dict(CLOSURE_IDENTITY),
        "settings": {
            "api_version": API_VERSION,
            "windows": {name: list(interval) for name, interval in WINDOW_SPECS},
            "parent": "left+child+right",
            "global": "left+child+right+ambient",
            "handoff_time": HANDOFF_TIME,
            "target_time": float(target_time),
            "step_caps": [float(cap) for cap in STEP_CAPS],
            "clock_locations": list(CLOCK_LOCATIONS),
            "clock_protocol": "trapezoid_coordinate_time_per_accepted_step",
            "rk_weight": "dt/6 * (f1 + 2*f2 + 2*f3 + f4)",
            "dx_divisions": 1,
            "piston_work_on_measurement_cut": False,
            "positive_part_net_flux_is_gross_lower_bound": True,
            "probability_flux": "occupation-weighted psi dagger sigma_2 psi",
        },
        "producer_hashes": producer_hashes(lab),
        "sealed_inputs": sealed_input_hashes(lab),
        "sealed_npz_parsed": False,
        "arrays_loaded": False,
        "evolved": False,
        "output_written": False,
        "initial_state_called": False,
        "production_trajectory_stored": False,
        "child_normal_closure_claimed": False,
        "global_normal_closure_claimed": False,
        "continuum_certified": False,
        "renewal_asserted": False,
        "cost": production_cost(lab, target_time=target_time),
        "launch": {
            "check": "python scripts/lab.py scripts/derive_nsc_discovery_balance.py --check",
            "evolve_two_h": (
                "python scripts/lab.py scripts/derive_nsc_discovery_balance.py --run "
                "--write results/development/nsc-discovery-balance-v1/ledger-T1.json"
            ),
            "evolve_spatial": (
                "python scripts/lab.py scripts/derive_nsc_discovery_balance.py --run --spatial "
                "--write results/development/nsc-discovery-balance-v1/ledger-spatial-T1.json"
            ),
            "launched_here": False,
        },
    }


def _pair_from_handoff(handoff):
    source = handoff["source_columns"]
    nf = int(handoff["nf"])
    record = {
        "nf": nf,
        "coarse_indices": handoff["coarse_indices"],
        "child_indices": handoff["child_indices"],
        "parent_indices": handoff["parent_indices"],
        "source_metadata": handoff["source_metadata"],
        "geometry_metadata": handoff["geometry_metadata"],
        "clock_locations": tuple(handoff["clocks"]["locations"]),
    }
    arrays = episode.arrays_from_state(
        handoff["state"],
        representation=episode.CANONICAL_PI,
        basis_arrays={
            "W": handoff["W"],
            "source_phi0": source[:nf],
            "source_phi1": source[nf:],
            "observer_columns": handoff["observer_columns"],
            "source_weights": handoff["weights"],
        },
        clocks={"rates": handoff["clocks"]["rates"], "normal_clocks": handoff["clocks"]["normal_clocks"]},
    )
    pair = episode.pair_from_arrays(arrays, record)
    state = handoff["state"].copy()
    pair, state, info = episode.resolve_pair(pair, state, backend="auto")
    if info.get("state_changed") or info.get("W_changed") or info.get("initial_state_called"):
        raise RuntimeError("resolver changed the frozen handoff")
    return pair, state, info


def _channel_gap(left, right):
    gap = abs(float(left) - float(right))
    scale = max(abs(float(left)), abs(float(right)))
    return {"absolute_gap": gap, "relative_gap": None if scale == 0.0 else float(gap / scale)}


def _compare_ledgers(left, right, label):
    keys = (
        "fixed_boundary_flux", "pressure", "lapse", "finite_projection_defect",
        "nyquist_gap", "pointwise_projection_defect", "coordinate_metric_work",
        "endpoint_measured_change", "quadrature_closure_error",
    )
    windows = {}
    for name in list(DISJOINT_WINDOWS) + ["parent", "global"]:
        windows[name] = {key: _channel_gap(left["windows"][name][key], right["windows"][name][key]) for key in keys}
    return {"label": label, "windows": windows}


def _run_handoff(handoff, step_cap, control_mode, *, target_time=TARGET_TIME, cpu_budget_seconds=CPU_BUDGET_SECONDS):
    import time
    pair, state, info = _pair_from_handoff(handoff)
    source_before = _source_token(pair)
    started = time.process_time()
    _final, ledger = integrate_segment(
        pair, state, HANDOFF_TIME, target_time, step_cap,
        control_mode=control_mode, normal_clocks=handoff["clocks"]["normal_clocks"],
        verify_first_step=True, owned_step_restriction=True,
        cpu_budget_seconds=cpu_budget_seconds,
    )
    elapsed = float(time.process_time() - started)
    if _source_token(pair) != source_before:
        raise RuntimeError("handoff source columns changed")
    return {
        "case_id": f"nf{int(handoff['nf'])}_{control_mode}_dt{float(step_cap):g}_balance",
        "nf": int(handoff["nf"]),
        "parent_case": handoff["case_name"],
        "control_mode": control_mode,
        "step_cap": float(step_cap),
        "handoff_time": HANDOFF_TIME,
        "target_time": float(target_time),
        "cpu_seconds": elapsed,
        "resolver": {key: info.get(key) for key in ("backend", "fft_applied", "state_changed", "W_changed", "initial_state_called")},
        "source_token_unchanged": True,
        "initial_state_called": False,
        "ledger": ledger,
        "initial_state_sha256": episode.state_sha256(state),
        "final_state_sha256": episode.state_sha256(_final),
    }


def evolve_frozen_segment(*, spatial=False, control_mode="coupled", lab=LAB, target_time=TARGET_TIME):
    """Re-evolve frozen ``T=0.3`` handoffs to a chosen coordinate endpoint.

    The default is two step caps on the nf256 handoff. ``spatial=True`` also
    runs the nf512 confirmation handoff at the matched cap. No trajectory
    arrays are retained. Input record hashes are compared before the return.
    """
    mode = episode.normalize_control_mode(control_mode)
    import time
    started = time.process_time()
    producers_before = producer_hashes(lab)
    before = sealed_input_hashes(lab)
    handoff = episode.load_saved_handoff(FINE_HANDOFF)
    runs = [_run_handoff(handoff, cap, mode, target_time=target_time,
                        cpu_budget_seconds=CPU_BUDGET_SECONDS - (time.process_time() - started)) for cap in STEP_CAPS]
    spatial_run = None
    if spatial:
        confirmation = episode.load_confirmation_handoff()
        spatial_run = _run_handoff(confirmation, episode.MATCHED_STEP_CAP, mode, target_time=target_time,
                                  cpu_budget_seconds=CPU_BUDGET_SECONDS - (time.process_time() - started))
        runs.append(spatial_run)
    after = sealed_input_hashes(lab)
    if producer_hashes(lab) != producers_before:
        raise RuntimeError("balance producer source changed during execution")
    if [(row["path"], row["sha256"], row["bytes"]) for row in before] != [
        (row["path"], row["sha256"], row["bytes"]) for row in after
    ]:
        raise RuntimeError("sealed input hashes changed during the ledger segment")
    by_cap = {float(run["step_cap"]): run for run in runs if run["nf"] == 256}
    comparison = {
        "two_h": _compare_ledgers(by_cap[float(episode.MATCHED_STEP_CAP)]["ledger"],
                                  by_cap[float(episode.STEP_CAP)]["ledger"],
                                  "nf256 dt 0.0005 versus dt 0.001"),
        "spatial": None,
    }
    if spatial_run is not None:
        comparison["spatial"] = _compare_ledgers(
            by_cap[float(episode.MATCHED_STEP_CAP)]["ledger"], spatial_run["ledger"],
            "nf256 versus nf512 at dt 0.0005",
        )
    return {
        "schema": SCHEMA,
        "api_version": API_VERSION,
        "identity": dict(CLOSURE_IDENTITY),
        "handoff_time": HANDOFF_TIME,
        "target_time": float(target_time),
        "producer_hashes": producers_before,
        "control_mode": mode,
        "spatial": bool(spatial),
        "cases": runs,
        "comparison": comparison,
        "input_hashes_before": before,
        "input_hashes_after": after,
        "hashes_unchanged": True,
        "cpu_seconds": float(time.process_time() - started),
        "case_integration_cpu_seconds": float(sum(run["cpu_seconds"] for run in runs)),
        "cpu_accounting": "serial process CPU including input binding, handoff preparation, and case integration",
        "evolved": True,
        "initial_state_called": False,
        "production_trajectory_stored": False,
        "creation_only": True,
        "child_normal_closure_claimed": False,
        "global_normal_closure_claimed": False,
        "continuum_certified": False,
        "renewal_asserted": False,
        "piston_work_on_measurement_cut": False,
        "cost_forecast": production_cost(lab, target_time=target_time),
    }


def verify_record(path, *, lab=LAB):
    """Authenticate input/source hashes and replay scalar accounting; no evolution."""
    import json
    path = Path(path)
    record = json.loads(path.read_text())
    problems = []
    if record.get("schema") != SCHEMA:
        problems.append("record is not a completed balance segment")
    if record.get("producer_hashes") != producer_hashes(lab):
        problems.append("current producer hashes differ from the frozen record")
    current = sealed_input_hashes(lab)
    pins = lambda rows: [(row["path"], row["sha256"], row["bytes"]) for row in rows]
    if pins(record.get("input_hashes_before", [])) != pins(current) or pins(record.get("input_hashes_after", [])) != pins(current):
        problems.append("sealed input hashes differ from the record")
    for case in record.get("cases", []):
        windows = case["ledger"]["windows"]
        for name, window in windows.items():
            accounted = sum(window[key] for key in ("fixed_boundary_flux", "pressure", "lapse", "finite_projection_defect", "reynolds"))
            if not np.isclose(accounted, window["accounted_change"], atol=1e-9, rtol=1e-12):
                problems.append(f"{case['case_id']} {name}: accounted channels disagree")
            if not np.isclose(window["normal_energy_end"] - window["normal_energy_start"], accounted + window["quadrature_closure_error"], atol=1e-9, rtol=1e-12):
                problems.append(f"{case['case_id']} {name}: endpoint balance disagrees")
            if not np.isclose(window["finite_projection_defect"], window["pointwise_projection_defect"] + window["nyquist_gap"], atol=1e-9, rtol=1e-12):
                problems.append(f"{case['case_id']} {name}: defect split disagrees")
            if not np.isclose(window["probability_gross_in_lower_bound"] - window["probability_gross_out_lower_bound"], window["probability_boundary_flux"], atol=1e-9, rtol=1e-12):
                problems.append(f"{case['case_id']} {name}: probability traffic disagrees")
        for name, parts in (("parent", PARENT_PARTS), ("global", DISJOINT_WINDOWS)):
            for key in RATE_KEYS:
                if not np.isclose(windows[name][key], sum(windows[part][key] for part in parts), atol=1e-8, rtol=1e-12):
                    problems.append(f"{case['case_id']} {name}: partition sum disagrees on {key}")
        if not case["ledger"]["first_step_bitwise_owned_rk4"]:
            problems.append(f"{case['case_id']}: first step did not bit-match the owned model")
    if not record.get("cases"):
        problems.append("record has no completed cases")
    return {"schema": CHECK_SCHEMA, "ok": not problems, "problems": problems,
            "record_sha256": sha256_file(path), "evolved": False,
            "trajectory_recomputed": False, "scalar_accounting_replayed": True,
            "continuum_certified": False}
