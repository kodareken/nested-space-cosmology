"""Spherical null expansions of the saved autonomous T=0.05 geometry.

Postprocessing only. The episode NPZ supplies the quadrature samples of
``r``, ``Q``, the lifted proper velocity, and ``K_perp``. The static clock
``L`` and ``beta`` are the coupling calibration on the same nodes. Nothing
here steps the state, rebuilds the column source, or recomputes inheritance.

Line element, signature ``+---``::

    g = r^2 [ L^2 dt^2 - Q^2 (dx + beta dt)^2 - dOmega^2 ]
    N = r L
    q = r Q
    n = (d_t - beta d_x) / N
    e1 = d_x / q
    k_plus = n + e1
    k_minus = n - e1

    theta_plus/minus = 2/r * ((r_dot - beta r_x)/N ± r_x/q)

    theta_plus * theta_minus = 4/r^2 * g^{ab} d_a r d_b r

``k_plus`` and ``k_minus`` are not rescaled by ``1/sqrt(2)``, so
``g(k_plus, k_minus) = 2``. A positive reciprocal rescaling
``k_plus -> lambda k_plus``, ``k_minus -> k_minus/lambda`` leaves the
product and each sign unchanged. Both expansions negative is trapped,
both positive is anti-trapped, opposite signs are untrapped, and a zero
expansion is the marginal indicator. Those are pointwise scalar labels on
this finite periodic chart. They are not a global event horizon.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np

from .nsc_spherical_coupling import (
    CALIBRATION,
    PERIOD,
    periodic_derivative,
    profile_coordinate,
)

SCHEMA = "NSC-SPHERICAL-NULL-EXPANSION-v1"
CPU_BUDGET_S = 20.0
JSON_BYTE_LIMIT = 1024 * 1024
MARGIN_FACTOR = 10.0
_LAB_ROOT = Path(__file__).resolve().parents[2]
_MODULE_PATH = Path(__file__).resolve()
EPISODE_JSON = _LAB_ROOT / "results" / "development" / "nsc-spherical-feedback-episode-v1.json"
EPISODE_NPZ = _LAB_ROOT / "results" / "development" / "nsc-spherical-feedback-episode-v1.npz"
RECORD_PATH = _LAB_ROOT / "results" / "development" / "nsc-spherical-null-expansion-v1.json"

RUNS = (
    {"name": "nf512_dt_0_0005", "nf": 512, "dt": 0.0005, "role": "primary"},
    {"name": "nf512_dt_0_00025", "nf": 512, "dt": 0.00025, "role": "time"},
    {"name": "nf256_dt_0_0005", "nf": 256, "dt": 0.0005, "role": "space"},
    {"name": "nf256_dt_0_00025", "nf": 256, "dt": 0.00025, "role": "supporting_time"},
)

SOURCE_FILES = {
    "coupling": _LAB_ROOT / "src" / "recursive_horizons" / "nsc_spherical_coupling.py",
    "galerkin": _LAB_ROOT / "src" / "recursive_horizons" / "nsc_spherical_galerkin_coupling.py",
    "regional_ledger": _LAB_ROOT / "src" / "recursive_horizons" / "nsc_regional_energy_exchange.py",
    "feedback_action": _LAB_ROOT / "src" / "recursive_horizons" / "nsc_spherical_feedback_action.py",
    "conformal_source": _LAB_ROOT / "src" / "recursive_horizons" / "nsc_conformal_adm_source.py",
    "episode_driver": _LAB_ROOT / "scripts" / "derive_nsc_spherical_feedback_episode.py",
    "episode_json": EPISODE_JSON,
    "episode_npz": EPISODE_NPZ,
}

EQUATIONS = {
    "signature": "+---",
    "line_element": "g = r^2 [ L^2 dt^2 - Q^2 (dx + beta dt)^2 - dOmega^2 ]",
    "lapse_and_radial_metric": "N = r L, q = r Q",
    "normal": "n = (d_t - beta d_x) / N",
    "radial_leg": "e1 = d_x / q",
    "null_legs": "k_plus = n + e1, k_minus = n - e1",
    "null_cross_norm": "g(k_plus, k_minus) = 2",
    "expansions": "theta_plus/minus = 2/r * ((r_dot - beta r_x)/N ± r_x/q)",
    "product": "theta_plus * theta_minus = 4/r^2 * g^{ab} d_a r d_b r",
    "observer_link": "K_perp = (r_dot - beta r_x) / (N r), so theta_plus + theta_minus = 4 K_perp",
    "classes": (
        "trapped: both expansions negative; anti-trapped: both positive; "
        "untrapped: opposite signs; marginal indicator: an expansion zero"
    ),
    "boost": (
        "For lambda > 0, k_plus -> lambda k_plus and k_minus -> k_minus/lambda. "
        "The product is unchanged and each sign is unchanged."
    ),
    "angular_sector": "r is constant on each symmetry sphere, so dOmega does not enter g^{ab} d_a r d_b r",
}


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def jsonable(value):
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, np.ndarray):
        return jsonable(value.tolist())
    if isinstance(value, np.floating):
        number = float(value)
        if not np.isfinite(number):
            raise ValueError("non-finite measurement")
        return number
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, float) and not np.isfinite(value):
        raise ValueError("non-finite measurement")
    return value


def tx_covariant(r, length_density, radial_factor, beta):
    """Covariant (t, x) block of the line element. Signature +---."""
    radius = np.asarray(r, dtype=float)
    lapse_factor = np.asarray(length_density, dtype=float)
    conformal = np.asarray(radial_factor, dtype=float)
    shift = np.asarray(beta, dtype=float)
    area = radius * radius
    conformal_square = conformal * conformal
    g_tt = area * (lapse_factor * lapse_factor - conformal_square * shift * shift)
    g_tx = -area * conformal_square * shift
    g_xx = -area * conformal_square
    return g_tt, g_tx, g_xx


def tx_inverse(r, length_density, radial_factor, beta):
    """Inverse (t, x) block. Raises if the lapse or radial metric vanishes."""
    radius = np.asarray(r, dtype=float)
    lapse_factor = np.asarray(length_density, dtype=float)
    conformal = np.asarray(radial_factor, dtype=float)
    shift = np.asarray(beta, dtype=float)
    if np.any(radius <= 0.0) or np.any(lapse_factor <= 0.0) or np.any(conformal <= 0.0):
        raise ValueError("positive r, L, and Q are required")
    lapse = radius * lapse_factor
    radial_metric = radius * conformal
    inverse_tt = 1.0 / lapse ** 2
    inverse_tx = -shift / lapse ** 2
    inverse_xx = -1.0 / radial_metric ** 2 + shift ** 2 / lapse ** 2
    return inverse_tt, inverse_tx, inverse_xx, lapse, radial_metric


def inverse_gradient_square(r, r_t, r_x, length_density, radial_factor, beta):
    """g^{ab} partial_a r partial_b r from the inverse metric, angular derivatives zero."""
    inverse_tt, inverse_tx, inverse_xx, _lapse, _radial = tx_inverse(
        r, length_density, radial_factor, beta
    )
    time_derivative = np.asarray(r_t, dtype=float)
    space_derivative = np.asarray(r_x, dtype=float)
    return (
        inverse_tt * time_derivative ** 2
        + 2.0 * inverse_tx * time_derivative * space_derivative
        + inverse_xx * space_derivative ** 2
    )


def _quadratic_form(g_tt, g_tx, g_xx, component_t, component_x):
    return (
        g_tt * component_t * component_t
        + 2.0 * g_tx * component_t * component_x
        + g_xx * component_x * component_x
    )


def frame_norms(r, length_density, radial_factor, beta):
    """Norms of n, e1, and k_plus/minus on the covariant metric."""
    g_tt, g_tx, g_xx = tx_covariant(r, length_density, radial_factor, beta)
    _inverse_tt, _inverse_tx, _inverse_xx, lapse, radial_metric = tx_inverse(
        r, length_density, radial_factor, beta
    )
    normal_t = 1.0 / lapse
    normal_x = -np.asarray(beta, dtype=float) / lapse
    radial_t = np.zeros_like(normal_t)
    radial_x = 1.0 / radial_metric
    plus_t = normal_t + radial_t
    plus_x = normal_x + radial_x
    minus_t = normal_t - radial_t
    minus_x = normal_x - radial_x
    cross = (
        g_tt * plus_t * minus_t
        + g_tx * (plus_t * minus_x + plus_x * minus_t)
        + g_xx * plus_x * minus_x
    )
    return {
        "g_nn": _quadratic_form(g_tt, g_tx, g_xx, normal_t, normal_x),
        "g_e1_e1": _quadratic_form(g_tt, g_tx, g_xx, radial_t, radial_x),
        "g_k_plus_k_plus": _quadratic_form(g_tt, g_tx, g_xx, plus_t, plus_x),
        "g_k_minus_k_minus": _quadratic_form(g_tt, g_tx, g_xx, minus_t, minus_x),
        "g_k_plus_k_minus": cross,
    }


def null_expansions(r, r_t, r_x, length_density, radial_factor, beta):
    """theta_plus/minus from k = n ± e1 acting on the areal radius."""
    _inverse_tt, _inverse_tx, _inverse_xx, lapse, radial_metric = tx_inverse(
        r, length_density, radial_factor, beta
    )
    radius = np.asarray(r, dtype=float)
    normal = (np.asarray(r_t, dtype=float) - np.asarray(beta, dtype=float) * np.asarray(r_x, dtype=float)) / lapse
    spatial = np.asarray(r_x, dtype=float) / radial_metric
    theta_plus = 2.0 / radius * (normal + spatial)
    theta_minus = 2.0 / radius * (normal - spatial)
    return theta_plus, theta_minus, theta_plus * theta_minus


def expansions_from_stored_normal(r, normal_velocity, r_x, radial_factor):
    """theta from the stored normal derivative (r_dot - beta r_x)/N and r_x/q."""
    radius = np.asarray(r, dtype=float)
    conformal = np.asarray(radial_factor, dtype=float)
    if np.any(radius <= 0.0) or np.any(conformal <= 0.0):
        raise ValueError("positive r and Q are required")
    radial_metric = radius * conformal
    spatial = np.asarray(r_x, dtype=float) / radial_metric
    normal = np.asarray(normal_velocity, dtype=float)
    theta_plus = 2.0 / radius * (normal + spatial)
    theta_minus = 2.0 / radius * (normal - spatial)
    return theta_plus, theta_minus, theta_plus * theta_minus, spatial


def boosted_expansions(theta_plus, theta_minus, scale):
    """Reciprocal rescaling of the two null legs by one positive factor."""
    factor = float(scale)
    if factor <= 0.0:
        raise ValueError("reciprocal null rescaling uses a positive factor")
    scaled_plus = factor * np.asarray(theta_plus, dtype=float)
    scaled_minus = np.asarray(theta_minus, dtype=float) / factor
    return scaled_plus, scaled_minus, scaled_plus * scaled_minus


def spectral_derivative(values, length):
    """Periodic spectral derivative. The Nyquist symbol is zero, as in the owner."""
    samples = np.asarray(values, dtype=float)
    count = samples.shape[-1]
    if count < 4 or count % 2:
        raise ValueError("periodic derivative expects an even point count")
    if not np.isfinite(length) or length <= 0.0:
        raise ValueError("positive period required")
    spectrum = np.fft.fft(samples, axis=-1)
    modes = np.fft.fftfreq(count) * count
    wavenumber = 2.0 * np.pi * modes / float(length)
    wavenumber = np.array(wavenumber, dtype=float, copy=True)
    wavenumber[count // 2] = 0.0
    derived = np.fft.ifft(spectrum * (1j * wavenumber), axis=-1)
    return np.real(derived)


def derivative_owner_gap(points=64, length=PERIOD):
    """FFT derivative against ``periodic_derivative``, including a Nyquist mode."""
    dense = periodic_derivative(points, length)
    coordinate = np.arange(points) * (float(length) / points)
    field = (
        np.sin(2.0 * np.pi * coordinate / length)
        + 0.2 * np.cos(6.0 * np.pi * coordinate / length)
        + 0.3 * np.cos(np.pi * np.arange(points))
    )
    fft_field = spectral_derivative(field, length)
    dense_field = dense @ field
    nyquist = np.cos(np.pi * np.arange(points))
    return {
        "points": int(points),
        "field_gap": float(np.max(np.abs(fft_field - dense_field))),
        "nyquist_symbol_max": float(np.max(np.abs(spectral_derivative(nyquist, length)))),
        "dense_nyquist_symbol_max": float(np.max(np.abs(dense @ nyquist))),
    }


def static_clock(count, length=PERIOD):
    """Calibration lapse factor and shift on the quadrature nodes. Not evolved."""
    count = int(count)
    spacing = float(length) / count
    coordinate = np.arange(count, dtype=float) * spacing
    profile = profile_coordinate(coordinate, length)
    scale = np.exp(np.log(CALIBRATION["Omega"]) * profile)
    length_density = CALIBRATION["b0"] * scale
    shift = CALIBRATION["beta0"] * scale
    return coordinate, length_density, shift


def tail_derivative_bound(values, length, mode_max):
    """L1 bound on the derivative contribution of Fourier modes above ``mode_max``."""
    samples = np.asarray(values, dtype=float)
    if samples.ndim == 1:
        samples = samples[None, :]
    count = samples.shape[-1]
    spectrum = np.fft.fft(samples, axis=-1) / count
    modes = np.fft.fftfreq(count) * count
    wavenumber = 2.0 * np.pi * np.abs(modes) / float(length)
    tail = np.abs(modes) > float(mode_max)
    if not np.any(tail):
        return 0.0
    return float(np.max(np.sum(np.abs(spectrum[:, tail]) * wavenumber[tail], axis=1)))


def geometry_mode_max(fermions):
    degree = int(fermions) - 1
    if degree < 3 or degree % 2 == 0:
        raise ValueError("geometry degree count must be the odd nf - 1")
    return degree // 2


def circular_intervals(mask):
    """Inclusive index runs on a circle. Empty if no sample is selected."""
    selected = np.asarray(mask, dtype=bool)
    count = int(selected.size)
    if count == 0 or not np.any(selected):
        return []
    indices = np.flatnonzero(selected)
    splits = np.where(np.diff(indices) > 1)[0]
    groups = np.split(indices, splits + 1)
    if len(groups) > 1 and int(groups[0][0]) == 0 and int(groups[-1][-1]) == count - 1:
        groups = [np.concatenate((groups[-1], groups[0]))] + groups[1:-1]
    intervals = []
    for group in groups:
        start = int(group[0])
        end = int(group[-1])
        width = int(group.size)
        wraps = bool(width > 1 and end < start)
        intervals.append({
            "start_index": start,
            "end_index": end,
            "count": width,
            "wraps": wraps,
        })
    return intervals


def _interval_coordinates(intervals, spacing, marker):
    located = []
    count_hint = None
    for interval in intervals:
        start = interval["start_index"]
        end = interval["end_index"]
        width = interval["count"]
        if interval["wraps"]:
            contains = bool(marker >= start or marker <= end)
        else:
            contains = bool(start <= marker <= end)
        located.append({
            **interval,
            "x_start": float(start * spacing),
            "x_end": float(end * spacing),
            "node_span": float((width - 1) * spacing),
            "contains_marker": contains,
        })
        count_hint = width
    return located, count_hint


def product_crossings(product):
    """Adjacent samples whose product changes sign, including exact zeros."""
    values = np.asarray(product, dtype=float)
    count = int(values.size)
    crossings = []
    for index in range(count):
        neighbor = (index + 1) % count
        left = float(values[index])
        right = float(values[neighbor])
        if left == 0.0 or left * right < 0.0:
            crossings.append({
                "left_index": index,
                "right_index": neighbor,
                "left_product": left,
                "right_product": right,
            })
    return crossings


def _classify_pair(theta_plus, theta_minus, margin):
    plus = np.asarray(theta_plus, dtype=float)
    minus = np.asarray(theta_minus, dtype=float)
    margin = float(margin)
    plus_negative = plus < -margin
    plus_positive = plus > margin
    minus_negative = minus < -margin
    minus_positive = minus > margin
    trapped = plus_negative & minus_negative
    anti = plus_positive & minus_positive
    untrapped = (plus_positive & minus_negative) | (plus_negative & minus_positive)
    within = ~(trapped | anti | untrapped)
    return trapped, anti, untrapped, within


def _extremum(radius, normal, theta_plus, theta_minus, product, r_x, spacing, kind):
    if kind == "maximum":
        index = int(np.argmax(radius))
    else:
        index = int(np.argmin(radius))
    return {
        "kind": kind,
        "index": index,
        "x": float(index * spacing),
        "r": float(radius[index]),
        "r_x": float(r_x[index]),
        "normal_velocity": float(normal[index]),
        "theta_plus": float(theta_plus[index]),
        "theta_minus": float(theta_minus[index]),
        "product": float(product[index]),
    }


def _critical_estimates(radius, normal, r_x, spacing, margin):
    """Linear zeros of r_x. At those points both expansions equal 2 v_n / r."""
    slope = np.asarray(r_x, dtype=float)
    count = int(slope.size)
    estimates = []
    for index in range(count):
        neighbor = (index + 1) % count
        left = float(slope[index])
        right = float(slope[neighbor])
        if left == 0.0 or left * right >= 0.0:
            continue
        weight = abs(left) / (abs(left) + abs(right))
        radius_zero = (1.0 - weight) * float(radius[index]) + weight * float(radius[neighbor])
        velocity = (1.0 - weight) * float(normal[index]) + weight * float(normal[neighbor])
        theta = 2.0 * velocity / radius_zero
        if abs(theta) <= margin:
            label = "marginal_within_margin"
        elif theta < 0.0:
            label = "both_negative"
        else:
            label = "both_positive"
        estimates.append({
            "left_index": index,
            "right_index": neighbor,
            "x": float((index + weight) * spacing),
            "r": float(radius_zero),
            "normal_velocity": float(velocity),
            "theta_both": float(theta),
            "label": label,
            "interpolation": "linear between adjacent samples",
        })
    return estimates


def summarize_frame(time, radius, conformal, normal, k_perp, r_x, length_density, shift, margin):
    theta_plus, theta_minus, product, spatial = expansions_from_stored_normal(
        radius, normal, r_x, conformal
    )
    raw_trapped, raw_anti, raw_untrapped, raw_within = _classify_pair(theta_plus, theta_minus, 0.0)
    strict_trapped, strict_anti, strict_untrapped, strict_within = _classify_pair(
        theta_plus, theta_minus, margin
    )
    spacing = float(PERIOD / radius.size)
    marker_max = int(np.argmax(radius))
    marker_min = int(np.argmin(radius))
    trapped_intervals, _width = _interval_coordinates(
        circular_intervals(raw_trapped), spacing, marker_max
    )
    anti_intervals, _anti_width = _interval_coordinates(
        circular_intervals(raw_anti), spacing, marker_min
    )
    norms = frame_norms(radius, length_density, conformal, shift)
    same_sign = raw_trapped | raw_anti
    min_same_sign = None
    if np.any(same_sign):
        min_same_sign = float(np.min(np.minimum(
            np.abs(theta_plus[same_sign]), np.abs(theta_minus[same_sign])
        )))
    observer_sum_gap = float(np.max(np.abs(theta_plus + theta_minus - 4.0 * k_perp)))
    spatial_difference_gap = float(np.max(np.abs(
        theta_plus - theta_minus - 4.0 * spatial / radius
    )))
    crossings = product_crossings(product)
    return {
        "time": float(time),
        "nodes": int(radius.size),
        "r_min": float(np.min(radius)),
        "r_max": float(np.max(radius)),
        "Q_min": float(np.min(conformal)),
        "Q_max": float(np.max(conformal)),
        "normal_min": float(np.min(normal)),
        "normal_max": float(np.max(normal)),
        "product_min": float(np.min(product)),
        "product_max": float(np.max(product)),
        "theta_plus_min": float(np.min(theta_plus)),
        "theta_plus_max": float(np.max(theta_plus)),
        "theta_minus_min": float(np.min(theta_minus)),
        "theta_minus_max": float(np.max(theta_minus)),
        "raw_counts": {
            "trapped": int(np.sum(raw_trapped)),
            "anti_trapped": int(np.sum(raw_anti)),
            "untrapped": int(np.sum(raw_untrapped)),
            "marginal_indicator": int(np.sum(raw_within)),
        },
        "strict_counts": {
            "trapped": int(np.sum(strict_trapped)),
            "anti_trapped": int(np.sum(strict_anti)),
            "untrapped": int(np.sum(strict_untrapped)),
            "within_margin": int(np.sum(strict_within)),
        },
        "same_sign_min_abs_expansion": min_same_sign,
        "product_crossing_count": len(crossings),
        "nodal_marginal_count": int(np.sum(product == 0.0)),
        "trapped_intervals": trapped_intervals,
        "anti_trapped_intervals": anti_intervals,
        "areal_maximum": _extremum(
            radius, normal, theta_plus, theta_minus, product, r_x, spacing, "maximum"
        ),
        "areal_minimum": _extremum(
            radius, normal, theta_plus, theta_minus, product, r_x, spacing, "minimum"
        ),
        "critical_point_estimates": _critical_estimates(
            radius, normal, r_x, spacing, margin
        ),
        "observer_sum_gap": observer_sum_gap,
        "spatial_difference_gap": spatial_difference_gap,
        "norm_gaps": {
            "g_nn_minus_one": float(np.max(np.abs(norms["g_nn"] - 1.0))),
            "g_e1_e1_plus_one": float(np.max(np.abs(norms["g_e1_e1"] + 1.0))),
            "g_k_plus": float(np.max(np.abs(norms["g_k_plus_k_plus"]))),
            "g_k_minus": float(np.max(np.abs(norms["g_k_minus_k_minus"]))),
            "g_k_cross_minus_two": float(np.max(np.abs(norms["g_k_plus_k_minus"] - 2.0))),
        },
        "_theta_plus": theta_plus,
        "_theta_minus": theta_minus,
        "_product": product,
        "_labels": (
            np.where(raw_trapped, -1, np.where(raw_anti, 1, np.where(raw_untrapped, 2, 3)))
        ),
    }


def _strip_private(frame):
    return {key: value for key, value in frame.items() if not key.startswith("_")}


def analyze_run(payload, spec, margin):
    name = spec["name"]
    fermions = int(spec["nf"])
    radius = np.asarray(payload[name + "_frame_r"], dtype=float)
    conformal = np.asarray(payload[name + "_frame_Q"], dtype=float)
    normal = np.asarray(payload[name + "_frame_proper"], dtype=float)
    k_perp = np.asarray(payload[name + "_frame_K_perp"], dtype=float)
    times = np.asarray(payload[name + "_frame_time"], dtype=float)
    if radius.ndim != 2 or conformal.shape != radius.shape or normal.shape != radius.shape:
        raise ValueError(f"{name} frame shapes disagree")
    if k_perp.shape != radius.shape or times.shape != (radius.shape[0],):
        raise ValueError(f"{name} clock or time shape disagrees")
    expected = 4 * fermions
    if radius.shape[1] != expected:
        raise ValueError(f"{name} quadrature {radius.shape[1]} is not 4*nf")
    expected_times = np.arange(0.0, 0.05 + 1e-12, 0.005)
    if times.shape != expected_times.shape or float(np.max(np.abs(times - expected_times))) > 1e-9:
        raise ValueError(f"{name} frames are not the 0.005 observer cadence through T=0.05")
    coordinate, length_density, shift = static_clock(radius.shape[1], PERIOD)
    r_x = spectral_derivative(radius, PERIOD)
    proper_gap = float(np.max(np.abs(k_perp * radius - normal)))
    mode_max = geometry_mode_max(fermions)
    tail_bound = tail_derivative_bound(radius, PERIOD, mode_max)
    frames = [
        summarize_frame(
            times[index], radius[index], conformal[index], normal[index], k_perp[index],
            r_x[index], length_density, shift, margin,
        )
        for index in range(radius.shape[0])
    ]
    positive = all(
        frame["r_min"] > 0.0 and frame["Q_min"] > 0.0 for frame in frames
    )
    return {
        "name": name,
        "role": spec["role"],
        "nf": fermions,
        "nq": int(radius.shape[1]),
        "ng": fermions - 1,
        "dt": float(spec["dt"]),
        "length": float(PERIOD),
        "dx": float(PERIOD / radius.shape[1]),
        "frames_stored": int(radius.shape[0]),
        "positive_r": bool(np.min(radius) > 0.0),
        "positive_Q": bool(np.min(conformal) > 0.0),
        "positive_geometry": bool(positive and np.min(length_density) > 0.0),
        "L_min": float(np.min(length_density)),
        "L_max": float(np.max(length_density)),
        "beta_min": float(np.min(shift)),
        "beta_max": float(np.max(shift)),
        "observer_proper_gap": proper_gap,
        "geometry_mode_max": int(mode_max),
        "high_mode_derivative_bound": tail_bound,
        "initial_normal_abs_max": float(np.max(np.abs(normal[0]))),
        "frames": frames,
    }


def _max_abs_pair(left, right, key):
    return float(max(
        np.max(np.abs(a[key] - b[key]))
        for a, b in zip(left, right)
    ))


def compare_runs(runs):
    by_name = {run["name"]: run for run in runs}
    fine = by_name["nf512_dt_0_0005"]
    fine_half = by_name["nf512_dt_0_00025"]
    coarse = by_name["nf256_dt_0_0005"]
    coarse_half = by_name["nf256_dt_0_00025"]
    time_fine = _max_abs_pair(fine["frames"], fine_half["frames"], "_theta_plus")
    time_fine_product = _max_abs_pair(fine["frames"], fine_half["frames"], "_product")
    time_coarse = _max_abs_pair(coarse["frames"], coarse_half["frames"], "_theta_plus")
    time_coarse_product = _max_abs_pair(coarse["frames"], coarse_half["frames"], "_product")
    count_keys = ("trapped", "anti_trapped", "untrapped")
    time_counts_match = True
    for primary, half in ((fine, fine_half), (coarse, coarse_half)):
        for frame, other in zip(primary["frames"], half["frames"]):
            if any(frame["raw_counts"][key] != other["raw_counts"][key] for key in count_keys):
                time_counts_match = False
    disagreements = []
    space_theta = 0.0
    space_product = 0.0
    for frame, other in zip(fine["frames"], coarse["frames"]):
        if abs(float(frame["time"]) - float(other["time"])) > 1e-12:
            raise ValueError("fine and coarse frame times differ")
        labels = frame["_labels"][::2]
        if labels.shape != other["_labels"].shape:
            raise ValueError("shared-node shape is not half the fine grid")
        disagreements.append(int(np.sum(labels != other["_labels"])))
        space_theta = max(space_theta, float(np.max(np.abs(frame["_theta_plus"][::2] - other["_theta_plus"]))))
        space_product = max(space_product, float(np.max(np.abs(frame["_product"][::2] - other["_product"]))))
    return {
        "dt_nf512_theta_plus_max_abs_gap": time_fine,
        "dt_nf512_product_max_abs_gap": time_fine_product,
        "dt_nf256_theta_plus_max_abs_gap": time_coarse,
        "dt_nf256_product_max_abs_gap": time_coarse_product,
        "dt_raw_counts_identical": bool(time_counts_match),
        "shared_node_sign_disagreements_by_time": disagreements,
        "shared_node_sign_disagreements_max": int(max(disagreements)),
        "space_theta_plus_max_abs_gap_on_shared_nodes": float(space_theta),
        "space_product_max_abs_gap_on_shared_nodes": float(space_product),
    }


def sign_margin(runs, comparison):
    theta_gap = max(
        comparison["dt_nf512_theta_plus_max_abs_gap"],
        comparison["dt_nf256_theta_plus_max_abs_gap"],
        comparison["space_theta_plus_max_abs_gap_on_shared_nodes"],
    )
    initial_floor = max(
        2.0 * run["initial_normal_abs_max"] / run["frames"][0]["r_min"] for run in runs
    )
    margin = MARGIN_FACTOR * max(theta_gap, initial_floor)
    return {
        "margin_factor": MARGIN_FACTOR,
        "comparison_theta_gap": float(theta_gap),
        "initial_normal_expansion_floor": float(initial_floor),
        "absolute_expansion_margin": float(margin),
        "rule": (
            "A sample expansion is strict only outside ten times the larger of "
            "the measured dt/shared-node theta gap and the initial 2|v_n|/r floor. "
            "A sample inside the margin is a marginal indicator, not a horizon."
        ),
    }


def _analytic_cases():
    return (
        {
            "name": "both_negative",
            "r": 2.0,
            "L": 1.0,
            "Q": 1.0,
            "beta": 0.0,
            "r_t": -0.5,
            "r_x": 0.0,
            "theta_plus": -0.25,
            "theta_minus": -0.25,
        },
        {
            "name": "both_positive",
            "r": 2.0,
            "L": 1.0,
            "Q": 1.0,
            "beta": 0.0,
            "r_t": 0.5,
            "r_x": 0.0,
            "theta_plus": 0.25,
            "theta_minus": 0.25,
        },
        {
            "name": "opposite_signs",
            "r": 2.0,
            "L": 1.0,
            "Q": 1.0,
            "beta": 0.0,
            "r_t": 0.0,
            "r_x": 0.5,
            "theta_plus": 0.25,
            "theta_minus": -0.25,
        },
        {
            "name": "shifted_null_frame",
            "r": 2.0,
            "L": 3.0,
            "Q": 4.0,
            "beta": 0.5,
            "r_t": -1.25,
            "r_x": 0.25,
            "theta_plus": None,
            "theta_minus": None,
        },
    )


def independent_controls():
    """Manufactured metric, norm, product, boost, and sign checks. No saved state."""
    product_gaps = []
    norm_gaps = []
    boost_product_gaps = []
    boost_sign_holds = True
    negative_scale_flips = True
    sign_gaps = []
    generator = np.random.default_rng(17)
    samples = list(_analytic_cases())
    for _index in range(8):
        samples.append({
            "name": "random",
            "r": float(generator.uniform(0.3, 3.0)),
            "L": float(generator.uniform(0.2, 2.0)),
            "Q": float(generator.uniform(0.2, 2.0)),
            "beta": float(generator.uniform(-1.0, 2.0)),
            "r_t": float(generator.uniform(-2.0, 2.0)),
            "r_x": float(generator.uniform(-2.0, 2.0)),
            "theta_plus": None,
            "theta_minus": None,
        })
    for sample in samples:
        theta_plus, theta_minus, product = null_expansions(
            sample["r"], sample["r_t"], sample["r_x"], sample["L"], sample["Q"], sample["beta"]
        )
        gradient = inverse_gradient_square(
            sample["r"], sample["r_t"], sample["r_x"], sample["L"], sample["Q"], sample["beta"]
        )
        predicted = 4.0 / sample["r"] ** 2 * gradient
        product_gaps.append(abs(float(product - predicted)))
        norms = frame_norms(sample["r"], sample["L"], sample["Q"], sample["beta"])
        norm_gaps.append(max(
            abs(float(norms["g_nn"] - 1.0)),
            abs(float(norms["g_e1_e1"] + 1.0)),
            abs(float(norms["g_k_plus_k_plus"])),
            abs(float(norms["g_k_minus_k_minus"])),
            abs(float(norms["g_k_plus_k_minus"] - 2.0)),
        ))
        for factor in (0.5, 2.0, 7.5):
            scaled_plus, scaled_minus, scaled_product = boosted_expansions(
                theta_plus, theta_minus, factor
            )
            scaled_plus = np.asarray(scaled_plus, dtype=float)
            scaled_minus = np.asarray(scaled_minus, dtype=float)
            boost_product_gaps.append(abs(float(scaled_product - product)))
            plus_values = np.asarray(theta_plus, dtype=float)
            minus_values = np.asarray(theta_minus, dtype=float)
            if np.any((plus_values != 0.0) & (np.sign(scaled_plus) != np.sign(plus_values))):
                boost_sign_holds = False
            if np.any((minus_values != 0.0) & (np.sign(scaled_minus) != np.sign(minus_values))):
                boost_sign_holds = False
        flipped_plus = -1.0 * float(theta_plus)
        flipped_minus = -1.0 * float(theta_minus)
        if theta_plus != 0.0 and flipped_plus * float(theta_plus) >= 0.0:
            negative_scale_flips = False
        if theta_minus != 0.0 and flipped_minus * float(theta_minus) >= 0.0:
            negative_scale_flips = False
        if sample["theta_plus"] is not None:
            sign_gaps.append(abs(float(theta_plus) - sample["theta_plus"]))
            sign_gaps.append(abs(float(theta_minus) - sample["theta_minus"]))
    derivative = derivative_owner_gap()
    arc = static_clock(8)
    # xi = 0 and xi = 2 lie on the packet arc where the profile coordinate is x.
    return {
        "inverse_metric_product_gap": float(max(product_gaps)),
        "null_norm_gap": float(max(norm_gaps)),
        "boost_product_gap": float(max(boost_product_gaps)),
        "positive_boost_preserves_each_sign": bool(boost_sign_holds),
        "negative_unit_rescaling_flips_both_signs": bool(negative_scale_flips),
        "analytic_sign_gap": float(max(sign_gaps)),
        "derivative_owner": derivative,
        "clock_arc_L_at_0": float(arc[1][0]),
        "clock_arc_beta_at_0": float(arc[2][0]),
        "clock_arc_L_at_2": float(arc[1][2]),
        "clock_arc_beta_at_2": float(arc[2][2]),
        "samples": len(samples),
        "saved_state_used": False,
    }


def _episode_metadata():
    record = json.loads(EPISODE_JSON.read_text())
    runs = {}
    for spec in RUNS:
        name = spec["name"]
        run = record["runs"][name]
        runs[name] = {
            "attained_T": float(run["attained_T"]),
            "positive_r": bool(run["positive_r"]),
            "positive_Q": bool(run["positive_Q"]),
            "initial_proper_min": float(run["initial"]["proper_min"]),
            "initial_proper_max": float(run["initial"]["proper_max"]),
            "final_proper_min": float(run["final"]["proper_min"]),
            "final_proper_max": float(run["final"]["proper_max"]),
        }
    return {
        "schema": record["schema"],
        "verdict": record["verdict"],
        "npz_sha256_recorded": record["npz_sha256"],
        "hashes_after": {
            key: record["hashes_after"][key]
            for key in (
                "coupling",
                "galerkin",
                "regional_module",
                "feedback_action",
                "conformal_source",
                "driver",
            )
        },
        "setup_length": {
            "nf512": float(record["setup_nf512"]["length"]),
            "nf256": float(record["setup_nf256"]["length"]),
        },
        "runs": runs,
    }


def _pattern(runs, margin_info):
    changes = True
    positive = True
    strict_same = True
    initial_same_sign = False
    later_both = True
    arcs_track_extrema = True
    for run in runs:
        frames = run["frames"]
        positive = positive and run["positive_geometry"]
        initial = frames[0]["raw_counts"]
        if initial["trapped"] or initial["anti_trapped"]:
            initial_same_sign = True
            changes = False
        if frames[0]["trapped_intervals"] or frames[0]["anti_trapped_intervals"]:
            arcs_track_extrema = False
        for frame in frames[1:]:
            counts = frame["raw_counts"]
            strict = frame["strict_counts"]
            if counts["trapped"] == 0 or counts["anti_trapped"] == 0:
                later_both = False
                changes = False
            if strict["trapped"] == 0 or strict["anti_trapped"] == 0:
                strict_same = False
            trapped_ok = (
                len(frame["trapped_intervals"]) == 1
                and frame["trapped_intervals"][0]["contains_marker"]
                and not frame["trapped_intervals"][0]["wraps"]
            )
            anti_ok = (
                len(frame["anti_trapped_intervals"]) == 1
                and frame["anti_trapped_intervals"][0]["contains_marker"]
                and not frame["anti_trapped_intervals"][0]["wraps"]
            )
            if not trapped_ok or not anti_ok:
                arcs_track_extrema = False
        if frames[0]["strict_counts"]["trapped"] or frames[0]["strict_counts"]["anti_trapped"]:
            strict_same = False
    return {
        "positive_geometry": bool(positive),
        "initial_samples_have_same_sign_expansion": bool(initial_same_sign),
        "later_frames_have_both_kinds": bool(later_both),
        "arcs_track_areal_extrema": bool(arcs_track_extrema),
        "sample_character_changes": bool(
            changes and later_both and arcs_track_extrema and not initial_same_sign and positive
        ),
        "strict_margin_agrees": bool(strict_same and changes and arcs_track_extrema),
        "absolute_expansion_margin": margin_info["absolute_expansion_margin"],
    }


def _clearance(runs, comparison):
    gap = max(
        comparison["dt_nf512_theta_plus_max_abs_gap"],
        comparison["space_theta_plus_max_abs_gap_on_shared_nodes"],
    )
    minima = [
        frame["same_sign_min_abs_expansion"]
        for run in runs
        for frame in run["frames"]
        if frame["same_sign_min_abs_expansion"] is not None
    ]
    smallest = float(min(minima))
    return {
        "smallest_same_sign_abs_expansion": smallest,
        "theta_comparison_gap": float(gap),
        "clearance_ratio": float(smallest / gap) if gap > 0.0 else None,
    }


def _finding(runs, comparison, pattern, clearance):
    fine = next(run for run in runs if run["name"] == "nf512_dt_0_0005")
    initial = fine["frames"][0]
    final = fine["frames"][-1]
    first = fine["frames"][1]
    return (
        "On every saved positive-geometry run, every quadrature sample at T=0 "
        "has opposite null expansions. The initial normal velocity reaches at most "
        f"{fine['initial_normal_abs_max']:.6e} on nf512, so a continuum areal "
        "critical point is a marginal surface within that noise: both expansions "
        "equal 2 v_n/r there. From the stored frame T=0.005 through T=0.05, each "
        "run has one both-negative arc at the areal maximum and one both-positive "
        "arc at the areal minimum, with the complementary samples untrapped. "
        f"The primary final counts are {final['raw_counts']['trapped']} both-negative, "
        f"{final['raw_counts']['anti_trapped']} both-positive, and "
        f"{final['raw_counts']['untrapped']} opposite-sign nodes out of {final['nodes']}. "
        f"At T=0.005 those counts are {first['raw_counts']['trapped']} and "
        f"{first['raw_counts']['anti_trapped']}. Halving dt moves theta_plus by at most "
        f"{comparison['dt_nf512_theta_plus_max_abs_gap']:.6e}. Shared nf512 and nf256 "
        f"nodes disagree in sign at {comparison['shared_node_sign_disagreements_max']} "
        "samples. The smallest same-sign expansion exceeds the comparison gap by "
        f"{clearance['clearance_ratio']:.6g}. Initial product max on the primary "
        f"grid is {initial['product_max']:.6e}."
    )


def _status(pattern, comparison):
    if not pattern["positive_geometry"]:
        return "POSITIVE_CHART_FAILED"
    if comparison["shared_node_sign_disagreements_max"] > 0:
        return "MEASURED_WITH_SHARED_NODE_SIGN_CONFLICT"
    if pattern["sample_character_changes"] and pattern["strict_margin_agrees"]:
        return "MEASURED_SAMPLE_TRAPPING_CHARACTER_CHANGES"
    if pattern["sample_character_changes"]:
        return "MEASURED_SAMPLE_CHANGE_INSIDE_SIGN_MARGIN"
    return "NO_SAMPLE_CHARACTER_CHANGE"


def build_record():
    """Read the saved episode and classify spherical expansions. Does not write."""
    started = time.process_time()
    wall = time.perf_counter()
    hashes_before = {name: sha256_file(path) for name, path in SOURCE_FILES.items()}
    hashes_before["module"] = sha256_file(_MODULE_PATH)
    metadata = _episode_metadata()
    controls = independent_controls()
    with np.load(EPISODE_NPZ, allow_pickle=False) as payload:
        provisional_margin = 0.0
        runs = [analyze_run(payload, spec, provisional_margin) for spec in RUNS]
    comparison = compare_runs(runs)
    margin_info = sign_margin(runs, comparison)
    with np.load(EPISODE_NPZ, allow_pickle=False) as payload:
        runs = [analyze_run(payload, spec, margin_info["absolute_expansion_margin"]) for spec in RUNS]
    comparison = compare_runs(runs)
    pattern = _pattern(runs, margin_info)
    clearance = _clearance(runs, comparison)
    metadata_gaps = []
    for run in runs:
        saved = metadata["runs"][run["name"]]
        initial = run["frames"][0]
        final = run["frames"][-1]
        metadata_gaps.append(abs(initial["normal_min"] - saved["initial_proper_min"]))
        metadata_gaps.append(abs(initial["normal_max"] - saved["initial_proper_max"]))
        metadata_gaps.append(abs(final["normal_min"] - saved["final_proper_min"]))
        metadata_gaps.append(abs(final["normal_max"] - saved["final_proper_max"]))
        if saved["positive_r"] != run["positive_r"] or saved["positive_Q"] != run["positive_Q"]:
            raise ValueError(f"episode metadata chart flag disagrees with stored {run['name']}")
        if abs(saved["attained_T"] - 0.05) > 1e-12:
            raise ValueError(f"{run['name']} did not attain T=0.05")
    hashes_after = {name: sha256_file(path) for name, path in SOURCE_FILES.items()}
    episode_unchanged = hashes_before["episode_json"] == hashes_after["episode_json"] and (
        hashes_before["episode_npz"] == hashes_after["episode_npz"]
    )
    hash_match = {
        "episode_npz_matches_metadata": hashes_after["episode_npz"] == metadata["npz_sha256_recorded"],
        "coupling_matches_episode_hashes_after": hashes_after["coupling"] == metadata["hashes_after"]["coupling"],
        "galerkin_matches_episode_hashes_after": hashes_after["galerkin"] == metadata["hashes_after"]["galerkin"],
        "ledger_matches_episode_hashes_after": (
            hashes_after["regional_ledger"] == metadata["hashes_after"]["regional_module"]
        ),
        "feedback_action_matches_episode_hashes_after": (
            hashes_after["feedback_action"] == metadata["hashes_after"]["feedback_action"]
        ),
        "conformal_source_matches_episode_hashes_after": (
            hashes_after["conformal_source"] == metadata["hashes_after"]["conformal_source"]
        ),
        "episode_driver_matches_episode_hashes_after": (
            hashes_after["episode_driver"] == metadata["hashes_after"]["driver"]
        ),
    }
    status = _status(pattern, comparison)
    cpu = time.process_time() - started
    record = {
        "schema": SCHEMA,
        "status": status,
        "role": (
            "Kinematic spherical null expansion of the stored T=0.05 geometry. "
            "Independent of the chi curvature identification."
        ),
        "evolution_performed": False,
        "inheritance_recomputed": False,
        "source_recomputed": False,
        "chi_used": False,
        "historical_1e-8_used_as_veto": False,
        "toy_B_or_lambda_comparison": False,
        "global_event_horizon_claimed": False,
        "black_hole_birth_claimed": False,
        "child_region_claimed": False,
        "regeneration_claimed": False,
        "sampling_is_continuous_horizon_proof": False,
        "equations": EQUATIONS,
        "chart": {
            "manifold": "finite periodic S^1 x S^2",
            "period": float(PERIOD),
            "coordinate": "x on [0, period), quadrature nodes of the episode observer",
            "sphere": "round symmetry sphere of areal radius r(t, x)",
            "scope": (
                "A scalar trapping label at a sample is a local property of the "
                "areal gradient. It does not identify a global causal interior, "
                "an exterior, or regeneration."
            ),
        },
        "clock": {
            "owner": "nsc_spherical_coupling static calibration, the episode gauge",
            "L": "b0 * Omega**s(x)",
            "beta": "beta0 * Omega**s(x)",
            "profile": "s=x on [0, 4], owner C^infty bridge on [4, period]",
            "normal_velocity": "stored frame_proper = (r_dot_lifted - beta r_x) / N",
            "K_perp": "stored frame_K_perp = normal_velocity / r",
            "ledger": "K_perp = (r_dot - beta r_x) / (N r) in nsc_regional_energy_exchange",
        },
        "pattern": pattern,
        "finding": _finding(runs, comparison, pattern, clearance),
        "resolved": [
            "Initial quadrature samples are opposite-sign on every saved run while r and Q stay positive.",
            "By each stored frame from T=0.005 through T=0.05, one both-negative arc lies at the areal maximum and one both-positive arc lies at the areal minimum.",
            "That split is unchanged by dt halving and by the shared nodes of nf512 and nf256.",
            "The product, null norms, and positive boost identity hold on the manufactured controls.",
        ],
        "open": [
            "The first time the product becomes positive lies somewhere in the unsampled interval (0, 0.005).",
            "Each product zero lies between adjacent nodes, so its continuum coordinate is known only to one grid spacing.",
            "A global event horizon, black-hole birth, a child region, and regeneration are outside this sample.",
            "The sign of 2 v_n/r exactly at the initial areal critical points lies inside the initial normal-velocity noise.",
            "No statement is made after T=0.05, or by identifying the expansions with chi.",
        ],
        "local_consequence": (
            "theta_plus * theta_minus = 4/r^2 g^{ab} partial_a r partial_b r is the "
            "causal character of the areal gradient on this chart. Its zero is the "
            "spherical gradient boundary. The saved samples move that boundary from "
            "the initial areal critical points, which are marginal within the normal "
            "velocity noise, onto the edges of a both-negative neighborhood of the "
            "maximum and a both-positive neighborhood of the minimum."
        ),
        "margin": margin_info,
        "clearance": clearance,
        "comparisons": comparison,
        "controls": controls,
        "metadata_normal_gap": float(max(metadata_gaps)),
        "source_bindings": {
            "hashes": hashes_after,
            "module_sha256": hashes_before["module"],
            "episode_metadata_verdict": metadata["verdict"],
            "episode_metadata_schema": metadata["schema"],
            "episode_npz_sha256_recorded_in_metadata": metadata["npz_sha256_recorded"],
            "hash_agreement": hash_match,
            "setup_length": metadata["setup_length"],
            "episode_files_unchanged": bool(episode_unchanged),
        },
        "runs": [{**run, "frames": [_strip_private(frame) for frame in run["frames"]]} for run in runs],
        "cpu_budget_seconds": CPU_BUDGET_S,
        "cpu_seconds": float(cpu),
        "wall_seconds": float(time.perf_counter() - wall),
        "cpu_budget_exceeded": bool(cpu > CPU_BUDGET_S),
    }
    if not all(hash_match.values()):
        record["status"] = "SOURCE_BINDING_MISMATCH"
    if not episode_unchanged:
        record["status"] = "EPISODE_BYTES_CHANGED_DURING_POSTPROCESS"
    return jsonable(record)


def write_record(record, path=RECORD_PATH):
    payload = json.dumps(jsonable(record), indent=2, allow_nan=False) + "\n"
    encoded = payload.encode("utf-8")
    if len(encoded) >= JSON_BYTE_LIMIT:
        raise ValueError(f"JSON record is {len(encoded)} bytes")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(encoded)
    return path
