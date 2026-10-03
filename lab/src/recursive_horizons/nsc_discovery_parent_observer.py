"""Leading parent observer on the signed-distance partition.

The measurement cuts are disjoint coordinate regions about the carrier
centre ``L/2``. Child is ``|s| <= 0.5``, the inner collar is
``0.5 < |s| <= 1.2``, the parent annulus is ``1.2 < |s| <= 3``, and the
ambient remainder ``3 < |s| <= L/2`` is the wrapped component. The
protected collar ``|s| <= 1`` overlaps child and inner, so it is not a
fifth accounted window.

Rates, metric jets, and the source come from the leading Einstein
diagnostic that is passed in, or from that same owner when no bundle is
passed. The state flow on the quadrature grid is ``prolong(decode(rate))``
for the geometry and the AP columns. Raw ``unprojected_rates`` are a
projection diagnostic and are not substituted for that flow. A chi-bearing
state or rate is refused. Source-free occupation weights may be zero; the
geometry rate stays the coupled leading rate and the label is preserved.
This module does not repeat the field equations, does not call the
auxiliary chi carrier, and does not add piston work on a measurement cut.
"""
from __future__ import annotations

import hashlib
import math

import numpy as np

from . import nsc_discovery_backend as backend
from . import nsc_discovery_episode as episode
from . import nsc_discovery_extent as extent
from . import nsc_discovery_leading_einstein as leading
from . import nsc_discovery_regions as regions
from . import nsc_regional_energy_exchange as regional


API_VERSION = "nsc-discovery-parent-observer-v1"
STATE_FLOW = "prolong(decode(leading_rate))"
STATE_FLOW_VERSION = "nsc-discovery-parent-observer-projected-flow-v2"
CHILD_RADIUS = 0.5
INNER_RADIUS = 1.2
PARENT_RADIUS = 3.0
PROTECTED_COLLAR_RADIUS = 1.0
REGION_ORDER = ("child", "inner", "parent", "ambient")
CUT_DISTANCES = (
    -PARENT_RADIUS, -INNER_RADIUS, -CHILD_RADIUS,
    CHILD_RADIUS, INNER_RADIUS, PARENT_RADIUS,
)
REQUIRED_PAIR_FIELDS = (
    "grid", "geometry_map", "weights", "reference_columns", "source_columns",
    "source_metadata", "geometry_metadata", "child_interval", "parent_interval",
    "clock_locations",
)
LEDGER_KEYS = (
    "reynolds", "fixed_boundary_flux", "pressure", "lapse",
    "pointwise_projection_defect", "endpoint_identity_gap",
    "finite_projection_defect", "slope_integral", "rate",
    "leibniz_from_slope", "closure_residual",
)
TIDE_CHANNELS = ("R_0101", "R_0202", "R4", "owned_W")
BASE_CHANNELS = ("Ricci2", "K", "normal_tidal_radial", "normal_tidal_angular")
CLOSURE = (
    "child |s|<=0.5; inner 0.5<|s|<=1.2; parent annulus 1.2<|s|<=3; "
    "ambient 3<|s|<=L/2; protected collar |s|<=1 is not a separate window"
)


def centre_of(length):
    length = float(length)
    if not math.isfinite(length) or length <= 0.0:
        raise ValueError("carrier length must be positive")
    return 0.5 * length


def signed_distance(coordinate, length):
    """Signed circular distance from ``L/2``, in ``(-L/2, L/2]``."""
    centre = centre_of(length)
    offset = np.asarray(coordinate, dtype=float) - centre
    return (offset + 0.5 * length) % length - 0.5 * length


def region_of_signed_distance(distance):
    magnitude = abs(float(distance))
    if magnitude <= CHILD_RADIUS:
        return "child"
    if magnitude <= INNER_RADIUS:
        return "inner"
    if magnitude <= PARENT_RADIUS:
        return "parent"
    return "ambient"


def region_partition(length):
    """Disjoint coordinate pieces. Shared cuts have no width."""
    length = float(length)
    centre = centre_of(length)
    if centre <= PARENT_RADIUS:
        raise ValueError("parent annulus does not fit inside one half-period")
    return {
        "centre": centre,
        "length": length,
        "closure": CLOSURE,
        "protected_collar_radius": PROTECTED_COLLAR_RADIUS,
        "protected_collar_accounted_separately": False,
        "accounted_parent_is_annulus": True,
        "overlapping_parent_sum_used": False,
        "child": {"pieces": [(centre - CHILD_RADIUS, centre + CHILD_RADIUS)], "wrapped": False},
        "inner": {
            "pieces": [
                (centre - INNER_RADIUS, centre - CHILD_RADIUS),
                (centre + CHILD_RADIUS, centre + INNER_RADIUS),
            ],
            "wrapped": False,
        },
        "parent": {
            "pieces": [
                (centre - PARENT_RADIUS, centre - INNER_RADIUS),
                (centre + INNER_RADIUS, centre + PARENT_RADIUS),
            ],
            "wrapped": False,
        },
        "ambient": {
            "pieces": [(centre + PARENT_RADIUS, centre - PARENT_RADIUS + length)],
            "wrapped": True,
            "backward": centre + PARENT_RADIUS,
            "forward": centre - PARENT_RADIUS,
        },
    }


def region_masks(grid):
    distance = np.asarray(signed_distance(grid.xi_q, grid.length), dtype=float)
    magnitude = np.abs(distance)
    masks = {
        "signed_distance": distance,
        "child": magnitude <= CHILD_RADIUS,
        "inner": (magnitude > CHILD_RADIUS) & (magnitude <= INNER_RADIUS),
        "parent": (magnitude > INNER_RADIUS) & (magnitude <= PARENT_RADIUS),
        "ambient": magnitude > PARENT_RADIUS,
    }
    covered = np.zeros(grid.nq, dtype=bool)
    for name in REGION_ORDER:
        if np.any(covered & masks[name]):
            raise RuntimeError("signed-distance masks overlap")
        covered |= masks[name]
    if not np.all(covered):
        raise RuntimeError("signed-distance masks do not cover the carrier")
    return masks


def null_crossings(length, time, centre):
    """Coordinate null rays ``dx/dt = ±1`` of the conformal chart ``β = 0``."""
    instant = float(time)
    period = float(length)
    tracks = []
    origins = []
    for distance in CUT_DISTANCES:
        origin = float((centre + distance) % period)
        origins.append(origin)
        for speed in (-1.0, 1.0):
            tracks.append({
                "signed_distance": float(distance),
                "origin": origin,
                "coordinate_speed": float(speed),
                "coordinate": float((origin + speed * instant) % period),
                "kind": "conformal_null_ray",
                "measured_boundary": False,
                "physical_wall": False,
            })
    ambiguities = []
    for track in tracks:
        for other in origins:
            separation = abs((track["coordinate"] - other + 0.5 * period) % period - 0.5 * period)
            own = abs((track["origin"] - other + 0.5 * period) % period - 0.5 * period)
            if separation <= 1.0e-8 and own > 1.0e-8:
                ambiguities.append({
                    "origin": track["origin"],
                    "coordinate_speed": track["coordinate_speed"],
                    "coordinate": track["coordinate"],
                    "coincides_with": float(other),
                    "kind": "null_image_of_another_cut",
                })
    return {
        "coordinate_speed": [-1.0, 1.0],
        "coordinate_speed_abs": 1.0,
        "beta": 0.0,
        "chart": "conformal N=rQ, dx/dt=±1",
        "crossing_time_one_side": {
            "child": 2.0 * CHILD_RADIUS,
            "inner": INNER_RADIUS - CHILD_RADIUS,
            "parent": PARENT_RADIUS - INNER_RADIUS,
            "ambient": 0.5 * period - PARENT_RADIUS,
        },
        "tracks": tracks,
        "ambiguities": ambiguities,
        "proper_clocks_integrated": False,
        "dispersal_time_asserted": False,
        "packet_width_used_as_boundary": False,
    }


def _piece_integral(grid, density, spec):
    if spec["wrapped"]:
        total = 0.0
        forward = float(spec["forward"])
        backward = float(spec["backward"])
        if forward > 0.0:
            total += extent.real_interval_integral(grid, density, (0.0, forward))
        if backward < float(grid.length):
            total += extent.real_interval_integral(grid, density, (backward, float(grid.length)))
        return float(total)
    return float(sum(
        extent.real_interval_integral(grid, density, interval) for interval in spec["pieces"]
    ))


def _fixed_ledger(grid, fields, spec):
    arguments = (
        fields["energy"], fields["flux"], fields["pressure"], fields["lapse"], fields["slope"],
    )
    if spec["wrapped"]:
        rate = regions.moving_normal_energy_on_cut(
            grid, *arguments, spec["backward"], spec["forward"], 0.0, 0.0, wrapped=True,
        )
        rate["status"] = "measured"
        return rate
    pieces = [
        regions.moving_normal_energy_rate(grid, *arguments, interval, (0.0, 0.0))
        for interval in spec["pieces"]
    ]
    if len(pieces) == 1:
        pieces[0]["status"] = "measured"
        return pieces[0]
    total = {key: float(sum(piece[key] for piece in pieces)) for key in LEDGER_KEYS}
    total.update({
        "pieces": [piece["interval"] for piece in pieces],
        "boundary_velocity": [0.0, 0.0],
        "wrapped": False,
        "surface_work": 0.0,
        "piston_work_included": False,
        "physical_wall": False,
        "dx_divisions": 1,
        "kind": "measurement_cut",
        "status": "measured",
    })
    return total


def _summary(values, mask):
    sample = np.asarray(values, dtype=float)[mask]
    if sample.size == 0 or not np.isfinite(sample).all():
        raise ValueError("metric summary needs finite samples in every region")
    return {
        "min": float(np.min(sample)),
        "max": float(np.max(sample)),
        "max_abs": float(np.max(np.abs(sample))),
    }


def _state_token(state):
    hasher = hashlib.sha256()
    for name in leading.FIELDS:
        array = np.ascontiguousarray(getattr(state, name))
        hasher.update(name.encode())
        hasher.update(array.tobytes())
    return hasher.hexdigest()


def _refuse_chi(value, label):
    if hasattr(value, "chi") or hasattr(value, "p_chi"):
        raise ValueError(label + " carries chi; the leading parent observer does not use chi rates")


def _require_pair(pair):
    missing = [name for name in REQUIRED_PAIR_FIELDS if not hasattr(pair, name)]
    if missing:
        raise ValueError("parent holder is missing " + ", ".join(missing))
    if pair.grid.gauge != "conformal":
        raise ValueError("leading parent observer requires the conformal gauge")


def observer_control_mode(control_mode):
    """Accept source-free locally. Do not patch the episode normalizer."""
    if control_mode in ("source_free", "source-free"):
        return "source_free"
    return episode.normalize_control_mode(control_mode)


def geometry_rate_mode(mode):
    """Source-free matter still rides the coupled geometry rate."""
    if mode == "source_free":
        return "coupled"
    return mode


def _weights(pair, state, fine):
    """Nonnegative occupations, including the source-free zero vector.

    A zero weight is an empty column, not a repaired phase. Nonfinite
    columns and weights outside ``[0, 1]`` are refused.
    """
    weights = np.asarray(pair.weights, dtype=float)
    if weights.ndim != 1 or weights.size < 1 or not np.isfinite(weights).all():
        raise ValueError("pair weights must be a finite occupation vector with a declared rank")
    if float(np.min(weights)) < 0.0 or float(np.max(weights)) > 1.0:
        raise ValueError("occupation weights must lie in [0, 1]")
    for name in leading.FIELDS:
        array = np.asarray(getattr(state, name))
        if not np.isfinite(array).all():
            raise ValueError(name + " has a nonfinite phase or sample")
    expected = (pair.grid.nf, int(weights.size))
    if np.shape(state.phi0) != expected or np.shape(state.phi1) != expected:
        raise ValueError("state columns must have shape (nf, weight rank)")
    if state.phi0.dtype.kind != "c" or state.phi1.dtype.kind != "c":
        raise ValueError("AP columns must be complex; an invalid phase is not coerced")
    if np.shape(fine.phi0)[1] != weights.size or np.shape(fine.phi1) != np.shape(fine.phi0):
        raise ValueError("fine columns do not match the actual weights")
    if not np.isfinite(fine.phi0).all() or not np.isfinite(fine.phi1).all():
        raise ValueError("fine columns have a nonfinite phase or sample")
    return weights


def _metadata(value):
    try:
        return regions.jsonable(value)
    except TypeError:
        return {"recorded": False, "reason": "metadata is not a scalar record"}


def _intervals(pair, partition):
    child = partition["child"]["pieces"][0]
    parent_span = (partition["centre"] - PARENT_RADIUS, partition["centre"] + PARENT_RADIUS)
    declared_child = tuple(float(value) for value in pair.child_interval)
    declared_parent = tuple(float(value) for value in pair.parent_interval)
    return {
        "declared_child_interval": [declared_child[0], declared_child[1]],
        "declared_parent_interval": [declared_parent[0], declared_parent[1]],
        "signed_child_interval": [float(child[0]), float(child[1])],
        "signed_parent_interval": [float(parent_span[0]), float(parent_span[1])],
        "matches_signed_distance": bool(
            np.allclose(declared_child, child, rtol=0.0, atol=1.0e-12)
            and np.allclose(declared_parent, parent_span, rtol=0.0, atol=1.0e-12)
        ),
        "declared_intervals_used_as_measurement_cuts": False,
    }


def _bundle_rate(pair, state, rate_mode, bundle):
    """Use a passed leading bundle, or obtain one from the leading owner."""
    if bundle is None:
        with backend.fft_thread_limit(1):
            rate, built = leading.rates(pair, state, return_bundle=True, control_mode=rate_mode)
            jets = leading.metric_jets(pair, state, rate, built, control_mode=rate_mode)
        return rate, built, jets, False
    _refuse_chi(bundle.get("fine_state"), "passed fine state")
    rate = bundle.get("leading_rate")
    if rate is None:
        raise ValueError(
            "passed bundle must carry leading_rate from the leading owner; chi rates are not a fallback"
        )
    _refuse_chi(rate, "passed leading rate")
    for name in ("source", "fine_state", "fine_system"):
        if name not in bundle:
            raise ValueError("passed bundle must carry the leading " + name)
    bundled_mode = bundle.get("control_mode")
    if bundled_mode is not None and episode.normalize_control_mode(bundled_mode) != rate_mode:
        raise ValueError("passed bundle control_mode does not match the geometry rate mode")
    jets = bundle.get("metric_jets")
    if jets is None:
        with backend.fft_thread_limit(1):
            jets = leading.metric_jets(pair, state, rate, bundle, control_mode=rate_mode)
    elif not isinstance(jets, tuple) or len(jets) != 3 or "tides" not in jets[0] or "base" not in jets[0]:
        raise ValueError("passed metric jets are not the leading curvature package")
    return rate, bundle, jets, True


def _component_gap(left, right):
    return float(np.max(np.abs(np.asarray(left) - np.asarray(right))))


def projected_state_rate(pair, rate, bundle=None):
    """Quadrature image of the encoded leading rate.

    ``rk4_step`` adds this encoded rate to the coefficient state. The
    quadrature fields and columns that move are ``prolong(decode(rate))``.
    A raw fine residual left outside that projection is reported and is not
    used as ``rho_t``, pressure, the normal-energy slope, or a contour speed.
    """
    for name in leading.FIELDS:
        if not hasattr(rate, name):
            raise ValueError("leading rate is missing " + name)
    _refuse_chi(rate, "leading rate")
    projected = leading.prolong(pair.grid, leading.decode(pair, rate))
    _refuse_chi(projected, "projected leading rate")
    diagnostic = {
        "state_flow": STATE_FLOW,
        "state_flow_version": STATE_FLOW_VERSION,
        "unprojected_rate_is_state_flow": False,
        "unprojected_rates_present": False,
        "geometry_projection_gap_max": None,
        "field_projection_gap_max": None,
        "full_projection_gap_max": None,
    }
    for name in leading.FIELDS:
        diagnostic[name] = None
    raw = None if bundle is None else bundle.get("unprojected_rates")
    if raw is None:
        return projected, diagnostic
    if len(tuple(raw)) != 6:
        raise ValueError("leading fine rate must have six fields")
    raw_state = leading.State(*raw)
    _refuse_chi(raw_state, "unprojected leading rate")
    gaps = {name: _component_gap(getattr(raw_state, name), getattr(projected, name)) for name in leading.FIELDS}
    geometry = max(gaps["Q"], gaps["r"], gaps["p_Q"], gaps["p_r"])
    field = max(gaps["phi0"], gaps["phi1"])
    diagnostic.update(gaps)
    diagnostic.update({
        "unprojected_rates_present": True,
        "geometry_projection_gap_max": geometry,
        "field_projection_gap_max": field,
        "full_projection_gap_max": max(geometry, field),
    })
    return projected, diagnostic


def _masked_peak(samples, rate, nodes):
    if nodes.size == 0:
        raise ValueError("region peak needs at least one node")
    subset_rate = None if rate is None else np.asarray(rate, dtype=float)[nodes]
    selected = regions.nodal_peak(np.asarray(samples, dtype=float)[nodes], subset_rate)
    selected["index"] = int(nodes[int(selected["peak_index"])])
    return selected


def observe(pair, state, time, control_mode="coupled", bundle=None):
    """JSON row for one leading state. Channels are scalar; cuts carry no ``-pV``."""
    _require_pair(pair)
    _refuse_chi(state, "state")
    for name in leading.FIELDS:
        if not np.isfinite(np.asarray(getattr(state, name))).all():
            raise ValueError(name + " has a nonfinite phase or sample")
    label = observer_control_mode(control_mode)
    rate_mode = geometry_rate_mode(label)
    instant = float(time)
    if not math.isfinite(instant):
        raise ValueError("observer time must be finite")
    before = _state_token(state)
    rate, bundle, jets, passed = _bundle_rate(pair, state, rate_mode, bundle)
    _refuse_chi(bundle["fine_state"], "fine state")
    if _state_token(state) != before:
        raise RuntimeError("parent observation mutated the state")
    fine = bundle["fine_state"]
    system = bundle["fine_system"]
    source = bundle["source"]
    grid = pair.grid
    weights = _weights(pair, state, fine)
    if system.occupations.shape != weights.shape or not np.allclose(system.occupations, weights):
        raise ValueError("fine occupations are not the holder's actual weights")
    if float(system.dx) != float(grid.dx_q):
        raise ValueError("fine spacing and carrier spacing differ; refusing a second dx")
    if system.derivative is not grid.derivative:
        raise ValueError("fine derivative is not the carrier Fourier derivative")
    fine_rate, projection = projected_state_rate(pair, rate, bundle)
    curvature, velocity, _acceleration = jets
    jet_gap = max(_component_gap(getattr(velocity, name), getattr(fine_rate, name)) for name in leading.FIELDS)
    if jet_gap > 1.0e-8:
        raise ValueError("passed metric jets do not match the prolonged decoded leading rate")
    for name in ("force_L", "force_Q", "force_beta"):
        if name not in source:
            raise ValueError("leading source is missing " + name)
    ledger = regional.matter_ledger(system, fine)
    for name in ("force_L", "force_Q", "force_beta"):
        gap = float(np.max(np.abs(np.asarray(ledger["source"][name]) - np.asarray(source[name]))))
        if gap > 1.0e-8:
            raise ValueError("matter ledger source does not match the passed leading source")
    terms = regional.proper_balance_terms(system, fine, fine_rate, ledger)
    slope = regions.analytic_normal_energy_slope(system, fine, fine_rate, source)
    fields = {
        "energy": np.asarray(ledger["normal_energy_nodal"], dtype=float),
        "flux": np.asarray(terms["flux_nodal"], dtype=float),
        "pressure": np.asarray(terms["proper_pressure_work"], dtype=float),
        "lapse": np.asarray(terms["momentum_lapse_work"], dtype=float),
        "slope": np.asarray(slope, dtype=float),
    }
    partition = region_partition(grid.length)
    masks = region_masks(grid)
    mass = np.sum((np.abs(fine.phi0) ** 2 + np.abs(fine.phi1) ** 2) * weights[None, :], axis=1)
    proper = np.asarray(fine.r, dtype=float) * np.asarray(fine.Q, dtype=float)
    if float(np.min(proper)) <= 0.0:
        raise ValueError("proper metric factor rQ must stay positive")
    per_proper = (mass / float(grid.dx_q)) / proper
    if float(np.min(per_proper)) < -1.0e-10:
        raise ValueError("proper-length probability went negative")
    matter_present = bool(float(np.sum(weights)) > 0.0 and float(np.sum(mass)) > 0.0)
    per_proper_rate = regions._probability_rate(fine, system, fine_rate)
    contours = regions.density_contours(grid, per_proper, per_proper_rate)
    global_region = region_of_signed_distance(float(masks["signed_distance"][int(contours["peak_index"])]))
    child_peak = _masked_peak(per_proper, per_proper_rate, np.flatnonzero(masks["child"]))
    ledgers = {name: _fixed_ledger(grid, fields, partition[name]) for name in REGION_ORDER}
    spacing = float(grid.dx_q)
    probability = {
        name: _piece_integral(grid, mass / spacing, partition[name]) for name in REGION_ORDER
    }
    probability["total"] = float(sum(probability[name] for name in REGION_ORDER))
    projections = {}
    for key, values in (
        ("shift", ledger["flux_shift"]),
        ("proper", ledger["flux_proper"]),
        ("shift_pressure", ledger["flux_shift_pressure"]),
    ):
        density = np.asarray(values, dtype=float) / spacing
        projections[key] = {
            name: _piece_integral(grid, density, partition[name]) for name in REGION_ORDER
        }
    tides = curvature["tides"]
    base = curvature["base"]
    metric = {
        "model": "leading_Einstein_EFT_diagnostic",
        "chi_rate_used": False,
        "chi_substituted_for_curvature": False,
        "R": {name: _summary(tides["R4"], masks[name]) for name in REGION_ORDER},
        "W": {name: _summary(tides["owned_W"], masks[name]) for name in REGION_ORDER},
        "tides": {
            channel: {name: _summary(tides[channel], masks[name]) for name in REGION_ORDER}
            for channel in TIDE_CHANNELS
        },
        "base": {
            channel: {name: _summary(base[channel], masks[name]) for name in REGION_ORDER}
            for channel in BASE_CHANNELS
        },
    }
    cut_locations = [float((partition["centre"] + distance) % float(grid.length)) for distance in CUT_DISTANCES]
    clock_locations = [float(value) for value in pair.clock_locations]
    proper_clocks = extent.real_periodic_values(grid, proper, clock_locations)
    cut_clocks = extent.real_periodic_values(grid, proper, cut_locations)
    channels = {
        "probability": probability,
        "normal_energy_rate": {name: ledgers[name]["rate"] for name in REGION_ORDER},
        "fixed_boundary_flux": {name: ledgers[name]["fixed_boundary_flux"] for name in REGION_ORDER},
        "pressure": {name: ledgers[name]["pressure"] for name in REGION_ORDER},
        "lapse": {name: ledgers[name]["lapse"] for name in REGION_ORDER},
        "finite_projection_defect": {
            name: ledgers[name]["finite_projection_defect"] for name in REGION_ORDER
        },
        "pointwise_projection_defect": {
            name: ledgers[name]["pointwise_projection_defect"] for name in REGION_ORDER
        },
        "flux_projection_shift": projections["shift"],
        "flux_projection_proper": projections["proper"],
        "flux_projection_shift_pressure": projections["shift_pressure"],
        "proper_clock_rates": [float(value) for value in proper_clocks],
        "cut_proper_clock_rates": [float(value) for value in cut_clocks],
        "null_coordinate_speed": [-1.0, 1.0],
        "piston_work_included": False,
        "dx_divisions": 1,
    }
    row = {
        "api_version": API_VERSION,
        "state_flow": STATE_FLOW,
        "state_flow_version": STATE_FLOW_VERSION,
        "time": instant,
        "control_mode": label,
        "geometry_rate_mode": rate_mode,
        "source_free": label == "source_free",
        "source_free_label_preserved": label == "source_free",
        "model": "leading_Einstein_EFT_diagnostic",
        "passed_bundle": bool(passed),
        "chi_rate_used": False,
        "state_token_unchanged": _state_token(state) == before,
        "rank": int(weights.size),
        "weights": [float(value) for value in weights],
        "reference_column_shape": [int(axis) for axis in np.shape(pair.reference_columns)],
        "source_column_shape": [int(axis) for axis in np.shape(pair.source_columns)],
        "source_metadata": _metadata(pair.source_metadata),
        "geometry_metadata": _metadata(pair.geometry_metadata),
        "holder_intervals": _intervals(pair, partition),
        "partition": {
            "centre": partition["centre"],
            "length": partition["length"],
            "closure": partition["closure"],
            "protected_collar_radius": PROTECTED_COLLAR_RADIUS,
            "protected_collar_accounted_separately": False,
            "accounted_parent_is_annulus": True,
            "overlapping_parent_sum_used": False,
            "pieces": {
                name: {
                    "wrapped": bool(partition[name]["wrapped"]),
                    "intervals": [
                        [float(item[0]), float(item[1])] for item in partition[name]["pieces"]
                    ],
                }
                for name in REGION_ORDER
            },
        },
        "positive_proper_probability": bool(
            float(np.min(per_proper)) >= 0.0 and float(np.sum(mass)) > 0.0
        ),
        "proper_probability_nonnegative": bool(float(np.min(per_proper)) >= -1.0e-10),
        "matter_present": matter_present,
        "matter_curves": "flat" if not matter_present else "occupied",
        "per_proper_length_minimum": float(np.min(per_proper)),
        "probability_total_nodal": float(np.sum(mass)),
        "child_peak": {
            "index": int(child_peak["index"]),
            "x": float(grid.xi_q[int(child_peak["index"])]),
            "signed_distance": float(masks["signed_distance"][int(child_peak["index"])]),
            "value": child_peak["peak_value"],
            "rate": child_peak["peak_rate"],
            "rate_status": child_peak["peak_rate_status"],
            "region": "child",
            "is_global_peak": bool(int(child_peak["index"]) == int(contours["peak_index"])),
            "replaces_global_contours": False,
        },
        "flags": {
            "parent_peak_dominates": bool(global_region == "parent"),
            "global_peak_outside_child": bool(global_region != "child"),
            "global_peak_region": global_region,
            "dominance_forced": False,
            "child_peak_tracked_separately": True,
        },
        "contours": {
            "peak_x": contours["peak_x"],
            "peak_value": contours["peak_value"],
            "peak_index": contours["peak_index"],
            "peak_rate": contours["peak_rate"],
            "peak_rate_status": contours["peak_rate_status"],
            "flags": contours["flags"],
            "global_peak_region": global_region,
            "levels": contours["levels"],
            "convention": contours["convention"],
            "continuum_optimized_peak": False,
            "contour_certified": False,
            "child_peak_replaced_global": False,
        },
        "ledgers": ledgers,
        "channels": channels,
        "metric": metric,
        "null_crossings": null_crossings(grid.length, instant, partition["centre"]),
        "proper_clocks": {
            "locations": clock_locations,
            "rates": [float(value) for value in proper_clocks],
            "cut_locations": cut_locations,
            "cut_rates": [float(value) for value in cut_clocks],
            "definition": "dτ/dT = rQ",
            "integrated": False,
        },
        "rates": {
            "Q_dot_max": float(np.max(np.abs(fine_rate.Q))),
            "r_dot_max": float(np.max(np.abs(fine_rate.r))),
            "state_flow": STATE_FLOW,
            "unprojected_rate_is_state_flow": False,
            "geometry_jets_zero": bool(
                float(np.max(np.abs(fine_rate.Q))) == 0.0 and float(np.max(np.abs(fine_rate.r))) == 0.0
            ),
        },
        "rate_projection": projection,
        "renewal_asserted": False,
        "continuum_certified": False,
        "contour_certified": False,
        "independent_energy_parcel": False,
        "column_ancestry_used_as_energy": False,
        "piston_work_included": False,
        "physical_wall": False,
        "held_population_measured": False,
        "constraint_preservation_claim": False,
    }
    if not row["state_token_unchanged"]:
        raise RuntimeError("parent observation mutated the state")
    return regions.jsonable(row)
