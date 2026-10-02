"""Stage-4 region consumer for one saved nested state.

This module does not step a state, rebuild ``W``, or call ``initial_state``.
It reads a native nested state, or a station chunk written by
``nsc_discovery_episode``, and reports three different localizations of that
same state:

- the six-column packet width of the total probability;
- the originally child-prepared source pair, columns ``(2, 3)``, with its
  own content, proper width, and retained fraction;
- connected gradient contours of the total proper-length density at one
  quarter, one half, and three quarters of the dominant peak.

A gradient contour is a measurement cut. The moving normal-energy rate on a
cut is the Leibniz--Reynolds identity

```
d/dt ∫_a^b η dx = η(b) ḃ − η(a) ȧ + F(a) − F(b) + pressure + lapse
                  + finite projection defect.
```

``η`` and ``F`` are the nodal shell and flux summands divided by ``dx`` once.
The slope of ``η`` is the analytic polarization of ``F_L/r`` owned by
``derive_nsc_nested_parent_child.normal_energy_slope``. The cut does not add
a piston ``−pV``. Column ancestry tags are computational labels, not separate
energy parcels. Modal complement columns are not the spatial exterior.
"""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np

from . import nsc_nested_parent_child as nested
from . import nsc_regional_energy_exchange as regional
from . import nsc_spherical_episode_assessment as metric
from .nsc_discovery_observables import (
    CARRIER_LENGTH,
    CHILD_INTERVAL,
    PARENT_INTERVAL,
    _lifted_rate,
    _prepare,
)
from .nsc_spherical_coupling import CauchyRate


API_VERSION = "nsc-discovery-regions-v1"
SCHEMA = "NSC-DISCOVERY-REGIONS-v1"
CHECK_SCHEMA = "NSC-DISCOVERY-REGIONS-CHECK-v1"
CHILD_SOURCE_COLUMNS = (2, 3)
COLUMN_PAIRS = ((0, 1), (2, 3), (4, 5))
CONTOUR_FRACTIONS = (0.25, 0.5, 0.75)
SHALLOW_SLOPE_FRACTION = 1.0e-3
DOMINANT_PEAK_FRACTION = 0.5
FLAT_CONTRAST = 1.0e-4
RESOLUTION_MODE_KEEP = 0.25
SEPARATED_SUPPORTS = ((0.0, 1.0), (1.0, 3.0), (3.0, 4.0))
PAIR_TAGS = (
    "outer_left_preparation",
    "originally_child_prepared_source_pair",
    "outer_right_preparation",
)
LAB = Path(__file__).resolve().parents[2]
SEALED_INPUTS = (
    "results/development/nsc-spherical-conformal-episode-v2.json",
    "results/development/nsc-spherical-conformal-episode-v2.npz",
    "results/development/nsc-nested-parent-child-v1.json",
    "results/development/nsc-nested-parent-child-v1.npz",
    "results/development/nsc-nested-parent-child-confirmation-v2.json",
    "results/development/nsc-nested-parent-child-confirmation-v2.npz",
    "results/development/nsc-nested-parent-child-replay-basis-v1.json",
    "results/development/nsc-nested-parent-child-replay-basis-v1.npz",
)
PRODUCERS = (
    "src/recursive_horizons/nsc_discovery_regions.py",
    "scripts/derive_nsc_discovery_regions.py",
    "src/recursive_horizons/nsc_discovery_observables.py",
    "src/recursive_horizons/nsc_discovery_episode.py",
    "src/recursive_horizons/nsc_discovery_backend.py",
    "src/recursive_horizons/nsc_nested_parent_child.py",
    "src/recursive_horizons/nsc_regional_energy_exchange.py",
    "src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
    "src/recursive_horizons/nsc_spherical_episode_assessment.py",
    "scripts/derive_nsc_nested_parent_child.py",
)

SETTINGS = {
    "api_version": API_VERSION,
    "contour_fractions": list(CONTOUR_FRACTIONS),
    "shallow_slope_fraction": SHALLOW_SLOPE_FRACTION,
    "dominant_peak_fraction": DOMINANT_PEAK_FRACTION,
    "flat_contrast": FLAT_CONTRAST,
    "resolution_mode_keep": RESOLUTION_MODE_KEEP,
    "child_source_columns": list(CHILD_SOURCE_COLUMNS),
    "child_interval": list(CHILD_INTERVAL),
    "parent_interval": list(PARENT_INTERVAL),
    "carrier_length": CARRIER_LENGTH,
    "dx_convention": "divide each nodal summand by dx once",
    "derivative": "owned periodic Fourier derivative; Nyquist symbol is zero",
    "evaluation": "Fourier interpolant through the nodes",
    "normal_energy_slope": "analytic polarization of F_L/r",
    "piston_work_on_measurement_cut": False,
    "renewal_asserted": False,
    "continuum_certified": False,
    "contour_certified": False,
}

DOMAIN_LIMITS = {
    "evolves_state": False,
    "calls_initial_state": False,
    "rebuilds_geometry_map": False,
    "resets_source_columns": False,
    "rewrites_station_chunks": False,
    "rewrites_sealed_inputs": False,
    "loads_production_trajectories_in_check": False,
    "measurement_cut_is_physical_wall": False,
    "piston_work_included": False,
    "independent_energy_parcel": False,
    "modal_complement_is_spatial_exterior": False,
    "ancestry_tags_are_computational_labels": True,
    "contour_certified": False,
    "continuum_certified": False,
    "renewal_asserted": False,
    "imposed_contour_is_regeneration_evidence": False,
    "output_policy": "creation_only",
    "reader": "station chunks from nsc_discovery_episode",
}

RENEWAL_CLASSIFICATION_REQUIRES = (
    "Saved station states at the handoff and at the later stations under comparison, "
    "for the coupled cases and the frozen-geometry controls, at each prepared resolution. "
    "This consumer reads chunks it is given and does not generate those trajectories.",
    "Contour endpoints and implicit velocities at one quarter, one half, and three quarters "
    "of the dominant peak on each of those stations, together with the resolution-comparison "
    "gaps from this row. A classification needs a pre-registered tolerance for those gaps.",
    "A pre-registered distinction between a maintained packet and a regenerated boundary. "
    "A gradient cut is a measurement surface and carries no wall traction.",
    "The stage-1 null characteristic tracks compared with the measured contour velocities "
    "at the same stations. This row does not turn that comparison into a verdict.",
    "The inherited local response of the realized source inside the measured contour. "
    "That response is a later stage and is not a region integral.",
    "A continuum or refinement certificate. contour_certified, continuum_certified, and "
    "renewal_asserted stay false, and this row does not estimate a state-error bound.",
)


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def jsonable(value):
    if isinstance(value, np.ndarray):
        raise TypeError("region records stay scalar")
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, (float, int, str, bool)) or value is None:
        return value
    raise TypeError(f"region record cannot hold {type(value).__name__}")


def producer_hashes(lab=LAB):
    lab = Path(lab)
    rows = {}
    for relative in PRODUCERS:
        path = lab / relative
        rows[relative] = sha256_file(path)
    return rows


def sealed_input_hashes(lab=LAB):
    """Hash sealed preparation bytes. Does not parse trajectories."""
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


def _state_token(state):
    hasher = hashlib.sha256()
    for name in nested.STATE_NAMES:
        array = np.ascontiguousarray(getattr(state, name))
        hasher.update(name.encode())
        hasher.update(array.tobytes())
    return hasher.hexdigest()


def _nodal(grid, values, name):
    array = np.asarray(values, dtype=float)
    if array.shape != (grid.nq,) or not np.isfinite(array).all():
        raise ValueError(f"{name} must be a finite nodal vector")
    return array


def _polarized_kernel_rates(system, state, rate):
    """Time derivatives of the owned nodal kernels ``K`` and ``S``."""
    momentum0 = system.momentum @ state.phi0
    momentum1 = system.momentum @ state.phi1
    momentum_rate0 = system.momentum @ rate.phi0
    momentum_rate1 = system.momentum @ rate.phi1
    weights = np.asarray(system.occupations, dtype=float)[None, :]
    kinetic = np.sum(weights * (
        np.conjugate(rate.phi0) * (-1j * momentum1)
        + np.conjugate(state.phi0) * (-1j * momentum_rate1)
        + np.conjugate(rate.phi1) * (1j * momentum0)
        + np.conjugate(state.phi1) * (1j * momentum_rate0)
    ), axis=1).real
    sigma = np.sum(weights * (
        np.conjugate(rate.phi0) * state.phi1
        + np.conjugate(state.phi0) * rate.phi1
        + np.conjugate(rate.phi1) * state.phi0
        + np.conjugate(state.phi1) * rate.phi0
    ), axis=1).real
    return kinetic, sigma


def analytic_normal_energy_slope(system, state, rate, source):
    """Nodal ``d/dt (F_L/r)``, the owned polarization.

    Same expression as ``normal_energy_slope`` in
    ``derive_nsc_nested_parent_child.py``. The result has the units of the
    nodal shell summand. Divide by ``dx`` once to obtain ``∂t η``.
    """
    if float(np.min(state.Q)) <= 0.0 or float(np.min(state.r)) <= 0.0:
        raise ValueError("normal-energy slope requires positive Q and r")
    kinetic_rate, sigma_rate = _polarized_kernel_rates(system, state, rate)
    force_rate = float(system.multiplicity) * (
        kinetic_rate / state.Q
        - np.asarray(source["K"], dtype=float) * rate.Q / state.Q ** 2
        + float(system.kappa) * sigma_rate
    )
    return force_rate / state.r - np.asarray(source["force_L"], dtype=float) * rate.r / state.r ** 2


def _freeze_geometry(lifted):
    zeros = np.zeros_like(np.asarray(lifted.Q, dtype=float))
    return CauchyRate(
        zeros, zeros, zeros, zeros, zeros, zeros,
        lifted.phi0, lifted.phi1,
        lifted.force_L, lifted.force_Q, lifted.force_beta, 0.0,
    )


def _sample(grid, values, coordinate):
    return float(nested._periodic_values(grid, values, (float(coordinate),))[0])


def _circular_distance(left, right, length):
    return abs((float(left) - float(right) + 0.5 * length) % length - 0.5 * length)


def moving_normal_energy_rate(grid, energy_summand, flux_summand, pressure_summand,
                               lapse_summand, slope_summand, interval, boundary_velocity):
    """Reynolds rate of the normal shell on one ordered measurement cut.

    Summands use the owned convention: the window integral is the Fourier
    integral of ``summand / dx``. Endpoint flux is ``(F(a) - F(b)) / dx``.
    The finite projection defect contains the pointwise balance residual and
    the gap between that collocation flux and the zero-Nyquist derivative.
    """
    spacing = float(grid.dx_q)
    if spacing <= 0.0:
        raise ValueError("positive node spacing required")
    left, right = (float(interval[0]), float(interval[1]))
    left_dot, right_dot = (float(boundary_velocity[0]), float(boundary_velocity[1]))
    if not 0.0 <= left < right <= float(grid.length):
        raise ValueError("moving window must be an ordered interval inside one period")
    if not all(math.isfinite(value) for value in (left_dot, right_dot)):
        raise ValueError("boundary velocity must be finite")
    fields = {
        "energy": _nodal(grid, energy_summand, "energy"),
        "flux": _nodal(grid, flux_summand, "flux"),
        "pressure": _nodal(grid, pressure_summand, "pressure"),
        "lapse": _nodal(grid, lapse_summand, "lapse"),
        "slope": _nodal(grid, slope_summand, "slope"),
    }
    densities = {name: values / spacing for name, values in fields.items()}
    divergence = np.asarray(grid.derivative @ fields["flux"], dtype=float)
    pointwise = fields["slope"] + divergence - fields["pressure"] - fields["lapse"]
    flux_ends = nested._periodic_values(grid, fields["flux"], (left, right))
    energy_ends = nested._periodic_values(grid, fields["energy"], (left, right))
    fixed_boundary_flux = float((flux_ends[0] - flux_ends[1]) / spacing)
    derivative_flux = nested.interval_integral(grid, divergence / spacing, (left, right))
    # Collocation inflow minus the derivative-matrix divergence, with the sign
    # that makes flux + pressure + lapse + defect reproduce ∫ slope/dx.
    endpoint_gap = float(-derivative_flux - fixed_boundary_flux)
    pointwise_defect = nested.interval_integral(grid, pointwise / spacing, (left, right))
    finite_defect = float(pointwise_defect + endpoint_gap)
    pressure = nested.interval_integral(grid, densities["pressure"], (left, right))
    lapse = nested.interval_integral(grid, densities["lapse"], (left, right))
    reynolds = float(energy_ends[1] / spacing * right_dot - energy_ends[0] / spacing * left_dot)
    slope_integral = nested.interval_integral(grid, densities["slope"], (left, right))
    rate = float(reynolds + fixed_boundary_flux + pressure + lapse + finite_defect)
    return {
        "interval": [left, right],
        "boundary_velocity": [left_dot, right_dot],
        "wrapped": False,
        "reynolds": reynolds,
        "fixed_boundary_flux": fixed_boundary_flux,
        "pressure": float(pressure),
        "lapse": float(lapse),
        "pointwise_projection_defect": float(pointwise_defect),
        "endpoint_identity_gap": endpoint_gap,
        "finite_projection_defect": finite_defect,
        "slope_integral": float(slope_integral),
        "rate": rate,
        "leibniz_from_slope": float(slope_integral + reynolds),
        "closure_residual": float(rate - (slope_integral + reynolds)),
        "surface_work": 0.0,
        "piston_work_included": False,
        "physical_wall": False,
        "dx_divisions": 1,
        "kind": "measurement_cut",
    }


def _add_rates(pieces):
    keys = (
        "reynolds", "fixed_boundary_flux", "pressure", "lapse",
        "pointwise_projection_defect", "endpoint_identity_gap",
        "finite_projection_defect", "slope_integral", "rate",
        "leibniz_from_slope", "closure_residual",
    )
    total = {key: float(sum(piece[key] for piece in pieces)) for key in keys}
    total.update({
        "pieces": [piece["interval"] for piece in pieces],
        "wrapped": True,
        "surface_work": 0.0,
        "piston_work_included": False,
        "physical_wall": False,
        "dx_divisions": 1,
        "kind": "measurement_cut",
    })
    return total


def moving_normal_energy_on_cut(grid, energy_summand, flux_summand, pressure_summand,
                                 lapse_summand, slope_summand, backward, forward,
                                 backward_velocity, forward_velocity, *, wrapped):
    """Rate on the connected cut, splitting a periodic component at ``x=0``.

    The split uses a fixed artificial cut at the branch point. Its two Reynolds
    contributions cancel, so the branch point is not a second physical boundary.
    """
    arguments = (energy_summand, flux_summand, pressure_summand, lapse_summand, slope_summand)
    backward = float(backward) % float(grid.length)
    forward = float(forward) % float(grid.length)
    if backward == forward:
        raise ValueError("measurement cut endpoints coincide")
    if not wrapped and backward < forward:
        return moving_normal_energy_rate(
            grid, *arguments, (backward, forward), (backward_velocity, forward_velocity),
        )
    if forward == 0.0:
        rate = moving_normal_energy_rate(
            grid, *arguments, (backward, float(grid.length)), (backward_velocity, forward_velocity),
        )
        rate["wrapped"] = True
        rate["interval"] = [backward, forward]
        rate["boundary_velocity"] = [float(backward_velocity), float(forward_velocity)]
        return rate
    pieces = []
    if forward > 0.0:
        pieces.append(moving_normal_energy_rate(
            grid, *arguments, (0.0, forward), (0.0, float(forward_velocity)),
        ))
    if backward < float(grid.length):
        pieces.append(moving_normal_energy_rate(
            grid, *arguments, (backward, float(grid.length)), (float(backward_velocity), 0.0),
        ))
    if not pieces:
        raise ValueError("wrapped measurement cut has no interior")
    if len(pieces) == 1:
        pieces[0]["wrapped"] = True
        return pieces[0]
    combined = _add_rates(pieces)
    combined["boundary_velocity"] = [float(backward_velocity), float(forward_velocity)]
    combined["interval"] = [backward, forward]
    return combined


def implicit_level_velocity(time_rate, space_slope, *, peak, spacing):
    """Level-set speed ``-ρ_t / ρ_x``. A shallow slope stays ambiguous."""
    time_rate = float(time_rate)
    space_slope = float(space_slope)
    peak = float(peak)
    spacing = float(spacing)
    if not all(math.isfinite(value) for value in (time_rate, space_slope, peak, spacing)) or spacing <= 0.0:
        raise ValueError("level velocity needs finite rates and a positive spacing")
    threshold = SHALLOW_SLOPE_FRACTION * abs(peak) / spacing
    report = {
        "time_rate": time_rate,
        "space_slope": space_slope,
        "threshold": float(threshold),
        "physical_wall": False,
    }
    if abs(space_slope) < threshold or peak == 0.0:
        report.update(status="ambiguous", velocity=None, division=False)
        return report
    report.update(status="measured", velocity=float(-time_rate / space_slope), division=True)
    return report


def _truncate_fourier(values):
    coefficients = np.fft.fft(np.asarray(values, dtype=float))
    modes = np.fft.fftfreq(values.size) * values.size
    limit = RESOLUTION_MODE_KEEP * values.size
    kept = np.abs(modes) <= limit + 1e-12
    return np.fft.ifft(np.where(kept, coefficients, 0.0)).real


def _local_maxima(values):
    left = np.roll(values, 1)
    right = np.roll(values, -1)
    return np.flatnonzero((values > left) & (values >= right))


def _component_count(values, level):
    above = np.asarray(values, dtype=float) >= float(level)
    if bool(np.all(above)):
        return 1
    previous = np.roll(above, 1)
    return int(np.count_nonzero(above & ~previous))


def _dominant_component(values, level, peak_index):
    count = values.size
    if values[peak_index] < level:
        return {"filled": False, "empty": True, "nodes": []}
    included = np.zeros(count, dtype=bool)
    included[peak_index] = True
    left = peak_index
    for _ in range(count - 1):
        nxt = (left - 1) % count
        if values[nxt] < level:
            break
        included[nxt] = True
        left = nxt
        if nxt == peak_index:
            break
    if int(np.count_nonzero(included)) == count:
        return {"filled": True, "empty": False, "nodes": included}
    right = peak_index
    for _ in range(count - 1):
        nxt = (right + 1) % count
        if values[nxt] < level or included[nxt]:
            break
        included[nxt] = True
        right = nxt
    return {
        "filled": False,
        "empty": False,
        "nodes": included,
        "inside_backward": left,
        "outside_backward": (left - 1) % count,
        "inside_forward": right,
        "outside_forward": (right + 1) % count,
        "wrapped": bool(left > right),
    }


def _refine_crossing(grid, values, level, outside, inside):
    length = float(grid.length)
    delta = (float(inside) - float(outside) + 0.5 * length) % length - 0.5 * length
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        coordinate = (float(outside) + mid * delta) % length
        sample = _sample(grid, values, coordinate)
        if sample >= level:
            hi = mid
        else:
            lo = mid
    return float((float(outside) + 0.5 * (lo + hi) * delta) % length)


def _contour_level(grid, density, density_rate, space_slope, level, peak_index, peak):
    component = _dominant_component(density, level, peak_index)
    report = {
        "level": float(level),
        "filled": bool(component["filled"]),
        "empty": bool(component.get("empty", False)),
        "component_count": _component_count(density, level),
        "wrapped": False,
        "backward": None,
        "forward": None,
        "velocity_backward": None,
        "velocity_forward": None,
        "resolution_gap_backward": None,
        "resolution_gap_forward": None,
        "resolution_status": "unresolved",
        "contour_certified": False,
        "physical_wall": False,
        "kind": "measurement_gradient_cut",
    }
    if component["filled"] or component.get("empty", False):
        return report, component
    backward = _refine_crossing(
        grid, density, level, grid.xi_q[component["outside_backward"]], grid.xi_q[component["inside_backward"]],
    )
    forward = _refine_crossing(
        grid, density, level, grid.xi_q[component["outside_forward"]], grid.xi_q[component["inside_forward"]],
    )
    report["backward"] = backward
    report["forward"] = forward
    # A boundary refined across x=0 lands near the period, so the coordinate
    # order is the wrap test. The nodal index order misses that one cell.
    report["wrapped"] = bool(backward > forward)
    if density_rate is not None:
        report["velocity_backward"] = implicit_level_velocity(
            _sample(grid, density_rate, backward), _sample(grid, space_slope, backward),
            peak=peak, spacing=float(grid.dx_q),
        )
        report["velocity_forward"] = implicit_level_velocity(
            _sample(grid, density_rate, forward), _sample(grid, space_slope, forward),
            peak=peak, spacing=float(grid.dx_q),
        )
    return report, component


def _resolution_gaps(grid, density, fraction, native):
    truncated = _truncate_fourier(density)
    peak = float(np.max(truncated))
    if peak <= 0.0:
        return None, None, "unresolved"
    peak_index = int(np.argmax(truncated))
    comparison, _component = _contour_level(
        grid, truncated, None, None, fraction * peak, peak_index, peak,
    )
    if native["backward"] is None or comparison["backward"] is None:
        return None, None, "unresolved"
    length = float(grid.length)
    return (
        _circular_distance(native["backward"], comparison["backward"], length),
        _circular_distance(native["forward"], comparison["forward"], length),
        "compared",
    )


def density_contours(grid, density, density_rate=None):
    """Dominant-peak contours of one nodal density. Certification stays false."""
    samples = _nodal(grid, density, "density")
    rate = None if density_rate is None else _nodal(grid, density_rate, "density_rate")
    peak = float(np.max(samples))
    floor = float(np.min(samples))
    contrast = (peak - floor) / max(abs(peak), 1.0e-30)
    flat = bool(peak <= 0.0 or contrast < FLAT_CONTRAST)
    peak_index = int(np.argmax(samples))
    slope = np.asarray(grid.derivative @ samples, dtype=float)
    maxima = _local_maxima(samples)
    dominant = [int(index) for index in maxima if samples[index] >= DOMINANT_PEAK_FRACTION * peak] if peak > 0.0 else []
    half_level = DOMINANT_PEAK_FRACTION * peak
    half = _dominant_component(samples, half_level, peak_index) if peak > 0.0 else {"filled": False, "empty": True, "nodes": None}
    if flat or half.get("empty", False):
        merged = 0
        split = False
    elif half.get("filled", False):
        merged = len(dominant)
        split = False
    else:
        inside_half = half["nodes"]
        merged = int(sum(1 for index in dominant if inside_half[index]))
        split = _component_count(samples, half_level) > 1
    flags = {
        "flat": flat,
        "merge": bool(not flat and merged > 1),
        "split": bool(not flat and split),
        "multiple_dominant_peaks": bool(not flat and len(dominant) > 1),
    }
    levels = []
    if not flat:
        for fraction in CONTOUR_FRACTIONS:
            level, _component = _contour_level(
                grid, samples, rate, slope, fraction * peak, peak_index, peak,
            )
            gap_backward, gap_forward, status = _resolution_gaps(grid, samples, fraction, level)
            level.update({
                "fraction": float(fraction),
                "resolution_gap_backward": gap_backward,
                "resolution_gap_forward": gap_forward,
                "resolution_status": status,
            })
            levels.append(level)
    return {
        "peak_x": float(grid.xi_q[peak_index]),
        "peak_value": peak,
        "peak_index": peak_index,
        "contrast": float(contrast),
        "minimum": floor,
        "flags": flags,
        "dominant_peak_count": 0 if flat else len(dominant),
        "levels": levels,
        "contour_certified": False,
        "physical_wall": False,
        "convention": "nodal samples, Fourier interpolant, zero-Nyquist derivative",
    }


def _proper_width(grid, mass, proper_q):
    packet = metric.packet_geometry(grid.xi_q, grid.length, mass, proper_q)
    resultant = float(packet["resultant"])
    if resultant > metric.NUMERICAL_RESULTANT_FLOOR:
        width = float(packet["proper_length"] * math.sqrt(-2.0 * math.log(min(resultant, 1.0))) / (2.0 * math.pi))
    else:
        width = None
    return width, resultant


def _supports(pair):
    closures = (pair.source_metadata or {}).get("source_support_closures")
    if closures is not None and len(closures) == 3:
        return tuple((float(item[0]), float(item[1])) for item in closures)
    return SEPARATED_SUPPORTS


def _column_mass(fine, system):
    weights = np.asarray(system.occupations, dtype=float)
    return (np.abs(fine.phi0) ** 2 + np.abs(fine.phi1) ** 2) * weights[None, :]


def _window_probability(grid, mass, interval):
    return nested.interval_integral(grid, np.asarray(mass, dtype=float) / float(grid.dx_q), interval)


def localization_report(pair, state, fine, system):
    """Separate the total packet, the child source pair, and the spatial windows."""
    if tuple(int(index) for index in pair.child_indices) != CHILD_SOURCE_COLUMNS:
        raise ValueError("child-prepared source pair is columns 2 and 3")
    grid = pair.grid
    spacing = float(grid.dx_q)
    column_mass = _column_mass(fine, system)
    total_mass = np.sum(column_mass, axis=1)
    proper_q = np.asarray(fine.r, dtype=float) * np.asarray(fine.Q, dtype=float)
    per_proper = (total_mass / spacing) / proper_q
    if float(np.min(per_proper)) < -1.0e-10:
        raise ValueError("proper-length probability went negative")
    total_width, total_resultant = _proper_width(grid, total_mass, proper_q)
    supports = _supports(pair)
    measured = nested.metrics(pair, state)
    pairs = []
    for tag, columns, support in zip(PAIR_TAGS, COLUMN_PAIRS, supports):
        mass = np.sum(column_mass[:, list(columns)], axis=1)
        content = float(np.sum(mass))
        width, resultant = _proper_width(grid, mass, proper_q)
        retained = None if content == 0.0 else _window_probability(grid, mass, support) / content
        pairs.append({
            "columns": list(columns),
            "ancestry": tag,
            "ancestry_is_computational_label": True,
            "independent_energy_parcel": False,
            "support": list(support),
            "content": content,
            "proper_width": width,
            "resultant": resultant,
            "retained_fraction": None if retained is None else float(retained),
        })
    child = pairs[1]
    child_probability = _window_probability(grid, total_mass, CHILD_INTERVAL)
    parent_probability = _window_probability(grid, total_mass, PARENT_INTERVAL)
    exterior = float(parent_probability - child_probability)
    modal_child = float(np.sum(column_mass[:, list(CHILD_SOURCE_COLUMNS)]))
    modal_complement = float(np.sum(total_mass) - modal_child)
    return {
        "six_column_content": float(np.sum(total_mass)),
        "six_column_total_width": total_width,
        "six_column_resultant": total_resultant,
        "column_pairs": pairs,
        "originally_child_prepared_pair": child,
        "per_proper_length_peak": float(np.max(per_proper)),
        "per_proper_length_peak_x": float(grid.xi_q[int(np.argmax(per_proper))]),
        "per_proper_length_minimum": float(np.min(per_proper)),
        "positive_probability": bool(float(np.min(per_proper)) >= 0.0 and float(np.sum(total_mass)) > 0.0),
        "child_probability": float(child_probability),
        "parent_probability": float(parent_probability),
        "child_probability_per_proper_length": float(child_probability / measured["child_proper_length"]),
        "parent_probability_per_proper_length": float(parent_probability / measured["parent_proper_length"]),
        "spatial_exterior_probability": exterior,
        "modal_child_content": modal_child,
        "modal_complement_content": modal_complement,
        "modal_complement_is_spatial_exterior": False,
        "q_multiplied_into_measure": False,
        "width_is_dispersal_time": False,
        "mass": total_mass,
        "per_proper": per_proper,
        "proper_q": proper_q,
        "measured": measured,
    }


def _probability_rate(fine, system, lifted):
    weights = np.asarray(system.occupations, dtype=float)[None, :]
    mass_rate = np.sum(weights * 2.0 * (
        np.real(np.conjugate(fine.phi0) * lifted.phi0)
        + np.real(np.conjugate(fine.phi1) * lifted.phi1)
    ), axis=1)
    spacing = float(system.dx)
    proper = np.asarray(fine.r, dtype=float) * np.asarray(fine.Q, dtype=float)
    proper_rate = np.asarray(lifted.r, dtype=float) * np.asarray(fine.Q, dtype=float) + np.asarray(fine.r, dtype=float) * np.asarray(lifted.Q, dtype=float)
    density = _column_mass(fine, system).sum(axis=1) / spacing
    density_rate = mass_rate / spacing
    return density_rate / proper - density * proper_rate / proper ** 2


def _energy_bundle(grid, ledger, terms, slope):
    return {
        "energy": np.asarray(ledger["normal_energy_nodal"], dtype=float),
        "flux": np.asarray(terms["flux_nodal"], dtype=float),
        "pressure": np.asarray(terms["proper_pressure_work"], dtype=float),
        "lapse": np.asarray(terms["momentum_lapse_work"], dtype=float),
        "slope": np.asarray(slope, dtype=float),
    }


def _cut_record(grid, fields, backward, forward, backward_velocity, forward_velocity, *, wrapped, status):
    if status != "measured":
        return {
            "status": status,
            "rate": None,
            "reynolds": None,
            "fixed_boundary_flux": None,
            "pressure": None,
            "lapse": None,
            "finite_projection_defect": None,
            "closure_residual": None,
            "piston_work_included": False,
            "physical_wall": False,
            "kind": "measurement_cut",
        }
    rate = moving_normal_energy_on_cut(
        grid, fields["energy"], fields["flux"], fields["pressure"], fields["lapse"], fields["slope"],
        backward, forward, backward_velocity, forward_velocity, wrapped=wrapped,
    )
    rate["status"] = "measured"
    return rate


def analyze(pair, state, time, *, bundle=None, nodal_rate=None, control_mode="coupled"):
    """Region row for one saved or native state. The state bytes stay in place."""
    if control_mode not in ("coupled", "frozen_geometry"):
        raise ValueError("control_mode must be coupled or frozen_geometry")
    instant = float(time)
    if not math.isfinite(instant):
        raise ValueError("region time must be finite")
    before = _state_token(state)
    _nodal_state, nodal_rate, bundle = _prepare(pair, state, bundle, nodal_rate)
    del _nodal_state
    if _state_token(state) != before:
        raise RuntimeError("region analysis mutated the state")
    grid = pair.grid
    fine = bundle["fine_state"]
    system = bundle["fine_system"]
    if float(system.dx) != float(grid.dx_q):
        raise ValueError("fine spacing and carrier spacing differ; refusing a second dx")
    if system.derivative is not grid.derivative:
        raise ValueError("fine derivative is not the carrier Fourier derivative")
    lifted = _lifted_rate(grid, nodal_rate, bundle)
    if control_mode == "frozen_geometry":
        lifted = _freeze_geometry(lifted)
    located = localization_report(pair, state, fine, system)
    ledger = regional.matter_ledger(system, fine)
    terms = regional.proper_balance_terms(system, fine, lifted, ledger)
    slope = analytic_normal_energy_slope(system, fine, lifted, bundle["source"])
    fields = _energy_bundle(grid, ledger, terms, slope)
    per_proper_rate = _probability_rate(fine, system, lifted)
    contours = density_contours(grid, located["per_proper"], per_proper_rate)
    kinetic_rate, sigma_rate = _polarized_kernel_rates(system, fine, lifted)
    fixed = {
        name: moving_normal_energy_rate(
            grid, fields["energy"], fields["flux"], fields["pressure"], fields["lapse"], fields["slope"],
            interval, (0.0, 0.0),
        )
        for name, interval in (("child", CHILD_INTERVAL), ("parent", PARENT_INTERVAL))
    }
    half = next((level for level in contours["levels"] if level["fraction"] == 0.5), None)
    if half is None or half["backward"] is None:
        contour_rate = _cut_record(grid, fields, None, None, None, None, wrapped=False, status="ambiguous" if not contours["flags"]["flat"] else "flat")
    else:
        velocities = (half["velocity_backward"], half["velocity_forward"])
        if any(item is None or item["status"] != "measured" for item in velocities):
            contour_rate = _cut_record(grid, fields, half["backward"], half["forward"], None, None, wrapped=half["wrapped"], status="ambiguous")
        else:
            contour_rate = _cut_record(
                grid, fields, half["backward"], half["forward"],
                velocities[0]["velocity"], velocities[1]["velocity"],
                wrapped=half["wrapped"], status="measured",
            )
    public_contours = []
    for level in contours["levels"]:
        public_contours.append({key: value for key, value in level.items()})
    closure_values = [fixed["child"]["closure_residual"], fixed["parent"]["closure_residual"]]
    if contour_rate.get("closure_residual") is not None:
        closure_values.append(contour_rate["closure_residual"])
    row = {
        "api_version": API_VERSION,
        "time": instant,
        "control_mode": control_mode,
        "settings": dict(SETTINGS),
        "domain_limits": dict(DOMAIN_LIMITS),
        "renewal_classification_requires": list(RENEWAL_CLASSIFICATION_REQUIRES),
        "renewal_asserted": False,
        "continuum_certified": False,
        "contour_certified": False,
        "source_reset": False,
        "initial_state_called": False,
        "historical_rewrite": False,
        "state_token_unchanged": _state_token(state) == before,
        "positive_probability": located["positive_probability"],
        "probability": {
            "total": located["six_column_content"],
            "per_proper_length_peak": located["per_proper_length_peak"],
            "per_proper_length_peak_x": located["per_proper_length_peak_x"],
            "per_proper_length_minimum": located["per_proper_length_minimum"],
            "child": located["child_probability"],
            "parent": located["parent_probability"],
            "child_per_proper_length": located["child_probability_per_proper_length"],
            "parent_per_proper_length": located["parent_probability_per_proper_length"],
        },
        "localization": {
            "six_column_content": located["six_column_content"],
            "six_column_total_width": located["six_column_total_width"],
            "column_pairs": located["column_pairs"],
            "originally_child_prepared_pair": located["originally_child_prepared_pair"],
            "spatial_exterior_probability": located["spatial_exterior_probability"],
            "modal_child_content": located["modal_child_content"],
            "modal_complement_content": located["modal_complement_content"],
            "modal_complement_is_spatial_exterior": False,
            "q_multiplied_into_measure": False,
            "width_is_dispersal_time": False,
        },
        "contours": {
            "peak_x": contours["peak_x"],
            "peak_value": contours["peak_value"],
            "contrast": contours["contrast"],
            "flags": contours["flags"],
            "dominant_peak_count": contours["dominant_peak_count"],
            "levels": public_contours,
            "contour_certified": False,
            "physical_wall": False,
            "convention": contours["convention"],
        },
        "normal_energy": {
            "child_fixed": fixed["child"],
            "parent_fixed": fixed["parent"],
            "dominant_half_maximum_cut": contour_rate,
            "formula": "reynolds + fixed_boundary_flux + pressure + lapse + finite_projection_defect",
            "analytical_continuity_closure_residual_max": float(np.max(np.abs(closure_values))),
        },
        "rates": {
            "Q_dot_max": float(np.max(np.abs(lifted.Q))),
            "r_dot_max": float(np.max(np.abs(lifted.r))),
            "field_column_rate_max": float(max(np.max(np.abs(lifted.phi0)), np.max(np.abs(lifted.phi1)))),
            "field_source_kernel_rate_max": float(max(np.max(np.abs(kinetic_rate)), np.max(np.abs(sigma_rate)))),
            "geometry_jets_zero": bool(
                float(np.max(np.abs(lifted.Q))) == 0.0 and float(np.max(np.abs(lifted.r))) == 0.0
            ),
        },
    }
    if not row["state_token_unchanged"]:
        raise RuntimeError("region analysis mutated the state")
    return jsonable(row)


def read_station_regions(directory, case_id, *, ordinal=None, snapshot_kinds=("station",)):
    """Read episode station chunks and analyze them. Chunks stay byte-identical."""
    from . import nsc_discovery_episode as episode

    directory = Path(directory)
    listed = episode.read_station_states(directory, case_id)
    chosen = []
    for item in listed:
        if item["snapshot_kind"] not in snapshot_kinds:
            continue
        if ordinal is not None and int(item["ordinal"]) != int(ordinal):
            continue
        chosen.append(item)
    if not chosen:
        raise FileNotFoundError("no matching station chunk")
    stations = []
    for item in chosen:
        stem = f"{case_id}-{int(item['ordinal']):06d}"
        npz_path = directory / f"{stem}.npz"
        json_path = directory / f"{stem}.json"
        before = {"npz": sha256_file(npz_path), "json": sha256_file(json_path)}
        record, arrays = episode.load_checkpoint(directory, case_id, item["ordinal"])
        pair = episode.pair_from_arrays(arrays, record)
        state = episode.state_from_arrays(arrays, record.get("momentum_representation", episode.CANONICAL_PI))
        mode = record.get("control_mode", "coupled")
        report = analyze(pair, state, float(record["coordinate_time"]), control_mode=mode)
        after = {"npz": sha256_file(npz_path), "json": sha256_file(json_path)}
        if before != after:
            raise RuntimeError("station chunk changed while reading")
        same_source = bool(
            np.array_equal(state.phi0, pair.source_phi0) and np.array_equal(state.phi1, pair.source_phi1)
        )
        stations.append({
            "ordinal": int(item["ordinal"]),
            "snapshot_kind": item["snapshot_kind"],
            "coordinate_time": float(record["coordinate_time"]),
            "chunk_sha256": before,
            "source_reset": False,
            "initial_state_called": False,
            "historical_rewrite": False,
            "state_equals_prepared_source": same_source,
            "regions": report,
        })
    return {
        "schema": SCHEMA,
        "case_id": str(case_id),
        "stations": stations,
        "settings": dict(SETTINGS),
        "domain_limits": dict(DOMAIN_LIMITS),
        "renewal_classification_requires": list(RENEWAL_CLASSIFICATION_REQUIRES),
        "producer_hashes": producer_hashes(),
        "sealed_preparations_parsed": False,
        "evolved": False,
        "source_reset": False,
        "historical_rewrite": False,
    }


def check_record(lab=LAB):
    """Hash producers and sealed inputs without parsing trajectory arrays."""
    return {
        "schema": CHECK_SCHEMA,
        "settings": dict(SETTINGS),
        "domain_limits": dict(DOMAIN_LIMITS),
        "renewal_classification_requires": list(RENEWAL_CLASSIFICATION_REQUIRES),
        "producer_hashes": producer_hashes(lab),
        "sealed_inputs": sealed_input_hashes(lab),
        "sealed_preparations_parsed": False,
        "arrays_loaded": False,
        "evolved": False,
        "output_written": False,
        "renewal_asserted": False,
        "contour_certified": False,
    }
