"""Assessment of one stored coupled-transfer episode. No evolution.

The running campaign and the source belong to the driver job
``coupled-transfer-event-episode``. This module only reads a stored chart,
a realized time jet, and the owned normal-energy ledger. It does not step
the state, change the gauge, add a force, or reset a field.

The 2-metric is the owned chart

    h = L^2 dt^2 - Q^2 (dx + beta dt)^2.

Its Ricci scalar uses the Christoffel contraction already owned by
``nsc_spherical_local_history.spherical_invariants`` and
``nsc_finite_terms.direct_curvature_identities``:

    Gamma^a_bc = (1/2) g^{ad} (d_c g_bd + d_b g_cd - d_d g_bc)
    R_bd = d_a Gamma^a_db - d_d Gamma^a_ab
           + Gamma^a_ae Gamma^e_db - Gamma^a_de Gamma^e_ab
    R_h = g^{bd} R_bd.

That contraction equals

    K = (Q_dot - d_x(beta Q)) / (L Q)
    R_h = 2/(L Q) d_x(L_x/Q) - 2*(((d_t K) - beta d_x K)/L + K^2).

``d_t K`` is built from the stored realized rate of ``Q_dot``. An
Euler-Lagrange substitute for ``p_chi`` is not accepted as that rate.
``chi^2/(3 r^4)`` stays an auxiliary proxy. The metric Weyl scalar is
``(R_h - 2)^2/(3 r^4)``.

Canonical columns are sqrt(dx) samples of u = r sqrt(q) psi. The flat
ell2 weight |Phi|^2 is the probability element. Proper arc phase uses
ds = q dx and does not multiply that weight by q. F_L/r shell samples
are a separate declared structure measure, not mode probability.
"""
from __future__ import annotations

import ast
import json
import zipfile
from pathlib import Path

import numpy as np

DRIVER_JOB = "coupled-transfer-event-episode"
OWNERS = {
    "action": "nsc_spherical_feedback_action",
    "geometry": "nsc_spherical_coupling",
    "observer": "nsc_spherical_null_expansion",
    "energy_ledger": "nsc_regional_energy_exchange",
}
ACCEPTED_JET_ORIGINS = ("realized_increment", "finite_ode_derivative")
REJECTED_JET_ORIGINS = (
    "el_p_chi",
    "el_p_chi_dot",
    "chi_equals_Rh_minus_2",
    "substituted_chi",
)
# Recorded here so a caller cannot mistake them for gates. They are not read.
PROXY_THRESHOLDS_NOT_APPLIED = {
    "reversal_at_least": 0.1,
    "leader_share_at_least": 0.5,
}
INDICATOR_KEYS = (
    "arc_location",
    "bridge_location",
    "packet_resultant",
    "occupation_concentration",
    "throughput_ratio",
    "work_integral",
    "content",
)
# Galerkin columns are sqrt(dx) samples of u = r sqrt(q) psi. Their flat
# ell2 sum is the Hilbert weight. Nodal |Phi|^2 is already the probability
# element; proper length q dx is only the arc phase.
CANONICAL_HALF_DENSITY = "canonical_half_density"
POSITIVE_SHELL_ENERGY = "positive_shell_energy"
ADMISSIBLE_PACKET_MEASURES = (CANONICAL_HALF_DENSITY, POSITIVE_SHELL_ENERGY)
# Same string as nsc_conformal_adm_source.HALF_DENSITY. The Hilbert weight is
# the flat ell2 sum of the sqrt(dx) samples of this u, not q|u|^2.
CANONICAL_FRAME = "u=r sqrt(q) psi"
NEIGHBOR_AGREEMENT_RTOL = 1e-8
NEIGHBOR_AGREEMENT_ATOL = 1e-10
# Roundoff gate for a circular resultant. Not a window-share requirement.
NUMERICAL_RESULTANT_FLOOR = 1e-8
_LAB = Path(__file__).resolve().parents[2]
DEFAULT_EPISODE_NPZ = _LAB / "results" / "development" / "nsc-spherical-feedback-episode-v1.npz"
REGENERATION_JSON = _LAB / "results" / "development" / "nsc-regeneration-episode-v1.json"
REGENERATION_NPZ = _LAB / "results" / "development" / "nsc-regeneration-episode-v1.npz"
REGENERATION_RUNS = (
    "nf256_dtmax_0_0005",
    "nf256_dtmax_0_00025",
    "nf512_dtmax_0_0005",
    "nf512_dtmax_0_00025",
)
REGENERATION_WINDOWS = np.array(
    [[0.0, 2.0], [2.0, 4.0], [4.0, 6.0], [6.0, 8.0]],
    dtype=float,
)
PACKET_EDGE = 3.96875
_FRAME_FIELDS_PRESENT = (
    "r", "Q", "chi", "proper", "K_r", "K_perp", "rho", "current", "weyl_C2",
)
_SAVED_RUN = "nf512_dt_0_0005"

INPUT_SCHEMA = {
    "driver_job": DRIVER_JOB,
    "binding": {
        "action_owner": OWNERS["action"],
        "geometry_owner": OWNERS["geometry"],
        "observer_owner": OWNERS["observer"],
        "energy_ledger_owner": OWNERS["energy_ledger"],
        "source_id": "nonempty string shared with the ledger",
        "phi_bound": True,
        "driver_job": DRIVER_JOB,
    },
    "chart": {
        "period": "positive finite period",
        "x": "(n,) periodic samples, n even and at least 4",
        "times": "(n_t,) strictly increasing finite window, n_t >= 3 for an episode",
        "L": "scalar, (n,), or (n_t, n). A length-n vector is spatial, not a time series",
        "beta": "same layout as L",
        "Q": "(n_t, n)",
        "r": "scalar, (n,), or (n_t, n)",
        "chi": "optional auxiliary (n_t, n). Absence is not chi = 0. Presence is not R_h - 2",
        "gauge_hold": (
            "True only when L and beta are the owned controls and their time "
            "derivatives are declared zero. Missing rates are not filled in."
        ),
        "L_dot": "required unless gauge_hold",
        "beta_dot": "required unless gauge_hold",
    },
    "time_jet": {
        "origin": "realized_increment or finite_ode_derivative",
        "projection_kept": True,
        "projection": "name of the projection that was kept",
        "Q_dot": "(n_t, n) realized partial_t Q",
        "Q_dot_rate": "(n_t, n) realized partial_t of Q_dot",
        "neighbor": {
            "dt": "positive step of the stored neighbor increment",
            "Q_dot_next": "(n_t, n) Q_dot at the neighbor step",
            "Q_next": "required for realized_increment: Q at the neighbor step",
            "projection_kept": True,
        },
    },
    "ledger": {
        "source_id": "same string as binding.source_id",
        "flux_is_budget_term": (
            "True. flux_budget_term is the signed term in "
            "d(content)/dt = flux_budget_term + pressure_work + lapse_exchange. "
            "It is the window reduction of minus the owned flux divergence."
        ),
        "normal_energy": "(n_t, n_windows) regional shell content",
        "flux_budget_term": "(n_t, n_windows)",
        "pressure_work": "(n_t, n_windows)",
        "lapse_exchange": "(n_t, n_windows)",
        "window_edges": "(n_windows, 2) periodic coordinate edges",
        "packet_weight_kind": (
            "canonical_half_density: nodal |Phi|^2 of u=r*sqrt(q)*psi, "
            "the flat ell2 weight, not multiplied by q. "
            "positive_shell_energy: declared nonnegative FL/r-class samples, "
            "a structure measure, not mode probability."
        ),
        "packet_weight": "(n_t, n) nonnegative nodal samples of the declared kind",
        "phi_basis_occupations": "optional (n_basis,) or (n_t, n_basis), unused for the arc",
        "field_energy": "optional (n_t,) diagnostic",
    },
    "comparison": {
        "label": "echoed and not used as a location",
        "refinement_indicator": {key: "nonnegative absolute movement" for key in INDICATOR_KEYS},
        "actual_bound": "ignored. This helper does not accept a bound.",
        "paired_episode": "optional second stored episode; indicators are absolute differences",
    },
}


def spectral_dx(values, period):
    """Periodic spectral derivative. The Nyquist symbol is zero, as in the geometry owner."""
    samples = np.asarray(values, dtype=float)
    count = samples.shape[-1]
    if count < 4 or count % 2:
        raise ValueError("periodic derivative expects an even point count of at least 4")
    if not np.isfinite(period) or period <= 0.0:
        raise ValueError("positive period required")
    spectrum = np.fft.fft(samples, axis=-1)
    wavenumber = 2.0 * np.pi * np.fft.fftfreq(count) * count / float(period)
    wavenumber = np.array(wavenumber, dtype=float, copy=True)
    wavenumber[count // 2] = 0.0
    derived = np.fft.ifft(spectrum * (1j * wavenumber), axis=-1)
    return np.real(derived)


def metric_components(L, Q, beta):
    """Covariant (t, x) block of h. Signature +---."""
    lapse = np.asarray(L, dtype=float)
    radial = np.asarray(Q, dtype=float)
    shift = np.asarray(beta, dtype=float)
    g_tt = lapse * lapse - radial * radial * shift * shift
    g_tx = -radial * radial * shift
    g_xx = -radial * radial
    return g_tt, g_tx, g_xx


def ricci_scalar_from_christoffel(g, dg, ddg):
    """R_h at one point from g_ij, d_a g_ij, and d_e d_c g_ij.

    ``dg[i, j, a] = d_a g_ij`` and ``ddg[e, i, j, c] = d_e d_c g_ij``.
    The contraction is the owned 2D Ricci scalar, not a sectional formula.
    """
    metric = np.asarray(g, dtype=float)
    first = np.asarray(dg, dtype=float)
    second = np.asarray(ddg, dtype=float)
    if metric.shape != (2, 2) or first.shape != (2, 2, 2) or second.shape != (2, 2, 2, 2):
        raise ValueError("Christoffel jets must be 2x2, 2x2x2, and 2x2x2x2")
    det = metric[0, 0] * metric[1, 1] - metric[0, 1] * metric[1, 0]
    if det >= 0.0 or not np.isfinite(det):
        raise ValueError("Lorentzian 2-metric required")
    inverse = np.array(
        [[metric[1, 1], -metric[0, 1]], [-metric[1, 0], metric[0, 0]]],
        dtype=float,
    ) / det
    gamma = np.zeros((2, 2, 2))
    dgamma = np.zeros((2, 2, 2, 2))
    for a in range(2):
        for b in range(2):
            for c in range(2):
                for d in range(2):
                    # Γ^a_bc = (1/2) g^{ad}(∂_b g_cd + ∂_c g_bd - ∂_d g_bc).
                    polynomial = first[c, d, b] + first[b, d, c] - first[b, c, d]
                    gamma[a, b, c] += 0.5 * inverse[a, d] * polynomial
                    for e in range(2):
                        d_polynomial = second[e, c, d, b] + second[e, b, d, c] - second[e, b, c, d]
                        d_inverse = 0.0
                        for p in range(2):
                            for q in range(2):
                                d_inverse += -inverse[a, p] * first[p, q, e] * inverse[q, d]
                        dgamma[e, a, b, c] += 0.5 * (d_inverse * polynomial + inverse[a, d] * d_polynomial)
    ricci = np.zeros((2, 2))
    for b in range(2):
        for d in range(2):
            for a in range(2):
                ricci[b, d] += dgamma[a, a, d, b] - dgamma[d, a, a, b]
                for e in range(2):
                    ricci[b, d] += (
                        gamma[a, a, e] * gamma[e, d, b] - gamma[a, d, e] * gamma[e, a, b]
                    )
    return float(np.einsum("bd,bd->", inverse, ricci))


def direct_rh_from_K(L, Q, beta, K, K_t, K_x, L_x, L_xx, Q_x):
    """R_h from K and its realized derivatives. L_xx and Q_x build d_x(L_x/Q)."""
    lapse = np.asarray(L, dtype=float)
    radial = np.asarray(Q, dtype=float)
    shift = np.asarray(beta, dtype=float)
    extrinsic = np.asarray(K, dtype=float)
    extrinsic_t = np.asarray(K_t, dtype=float)
    extrinsic_x = np.asarray(K_x, dtype=float)
    lapse_x = np.asarray(L_x, dtype=float)
    lapse_xx = np.asarray(L_xx, dtype=float)
    radial_x = np.asarray(Q_x, dtype=float)
    if np.any(lapse <= 0.0) or np.any(radial <= 0.0):
        raise ValueError("positive L and Q are required")
    slope = lapse_x / radial
    slope_x = (lapse_xx * radial - lapse_x * radial_x) / radial**2
    convective = (extrinsic_t - shift * extrinsic_x) / lapse
    curvature = 2.0 / (lapse * radial) * slope_x - 2.0 * (convective + extrinsic**2)
    return curvature


def _positive_chart_mask(L, Q, r):
    return (
        np.isfinite(L).all()
        and np.isfinite(Q).all()
        and np.isfinite(r).all()
        and np.min(L) > 0.0
        and np.min(Q) > 0.0
        and np.min(r) > 0.0
    )


def direct_rh_grid(L, Q, beta, Q_dot, Q_dot_rate, period, *, L_dot, beta_dot):
    """Grid R_h. Time derivatives of the gauge are explicit arguments."""
    lapse = np.asarray(L, dtype=float)
    radial = np.asarray(Q, dtype=float)
    shift = np.asarray(beta, dtype=float)
    radial_t = np.asarray(Q_dot, dtype=float)
    radial_tt = np.asarray(Q_dot_rate, dtype=float)
    lapse_t = np.asarray(L_dot, dtype=float)
    shift_t = np.asarray(beta_dot, dtype=float)
    lapse_x = spectral_dx(lapse, period)
    radial_x = spectral_dx(radial, period)
    shift_x = spectral_dx(shift, period)
    radial_tx = spectral_dx(radial_t, period)
    lapse_tx = spectral_dx(lapse_t, period)
    shift_tx = spectral_dx(shift_t, period)
    lapse_xx = spectral_dx(lapse_x, period)
    beta_Q_x = shift_x * radial + shift * radial_x
    numerator = radial_t - beta_Q_x
    denominator = lapse * radial
    extrinsic = numerator / denominator
    beta_Q_tx = shift_tx * radial + shift_t * radial_x + shift_x * radial_t + shift * radial_tx
    numerator_t = radial_tt - beta_Q_tx
    denominator_t = lapse_t * radial + lapse * radial_t
    extrinsic_t = (numerator_t * denominator - numerator * denominator_t) / denominator**2
    extrinsic_x = spectral_dx(extrinsic, period)
    curvature = direct_rh_from_K(
        lapse, radial, shift, extrinsic, extrinsic_t, extrinsic_x, lapse_x, lapse_xx, radial_x,
    )
    return {
        "K": extrinsic,
        "K_t": extrinsic_t,
        "K_x": extrinsic_x,
        "R_h": curvature,
        "L_x": lapse_x,
    }


def weyl_actual(R_h, r):
    """(R_h - 2)^2 / (3 r^4). This is not the chi proxy."""
    radius = np.asarray(r, dtype=float)
    if np.any(radius <= 0.0):
        raise ValueError("positive areal radius required")
    return (np.asarray(R_h, dtype=float) - 2.0) ** 2 / (3.0 * radius**4)


def weyl_proxy_from_chi(chi, r):
    """chi^2 / (3 r^4). Equal to the metric Weyl scalar only on the shell chi = R_h - 2."""
    radius = np.asarray(r, dtype=float)
    if np.any(radius <= 0.0):
        raise ValueError("positive areal radius required")
    return np.asarray(chi, dtype=float) ** 2 / (3.0 * radius**4)


def proper_clocks(r, L, Q):
    """N = r L and q = r Q."""
    radius = np.asarray(r, dtype=float)
    return radius * np.asarray(L, dtype=float), radius * np.asarray(Q, dtype=float)


def circular_distance(left, right, period):
    delta = (float(left) - float(right) + 0.5 * period) % period - 0.5 * period
    return abs(delta)


def _phase_moment(phase, weight, harmonic):
    total = float(np.sum(weight))
    angle = float(harmonic) * np.asarray(phase, dtype=float)
    cosine = float(np.sum(weight * np.cos(angle)) / total)
    sine = float(np.sum(weight * np.sin(angle)) / total)
    resultant = float(np.hypot(cosine, sine))
    raw = float(np.atan2(sine, cosine) % (2.0 * np.pi))
    axis = raw / float(harmonic)
    return resultant, axis


def _proper_phase(x, proper_q, period):
    """Phase 2π s/S of the proper radial arc ds = q dx. q does not weight probability."""
    coordinate = np.asarray(x, dtype=float)
    radial = np.asarray(proper_q, dtype=float)
    if coordinate.shape != radial.shape or coordinate.size < 2:
        raise ValueError("proper arc needs q on the same nodes as x")
    if np.any(radial <= 0.0) or not np.isfinite(radial).all():
        raise ValueError("positive proper radial metric q required")
    gaps = np.diff(coordinate)
    if np.any(gaps <= 0.0):
        raise ValueError("x must be strictly increasing on one period")
    closing = float(coordinate[0] + float(period) - coordinate[-1])
    if closing <= 0.0:
        raise ValueError("x must lie on one positive period")
    steps = np.append(gaps, closing)
    element = radial * steps
    total = float(np.sum(element))
    cumulative = np.cumsum(element) - element
    phase = 2.0 * np.pi * cumulative / total
    return phase, total


def _phase_to_coordinate(phase, node_phase, coordinate, period):
    delta = (np.asarray(node_phase) - float(phase) + np.pi) % (2.0 * np.pi) - np.pi
    index = int(np.argmin(np.abs(delta)))
    return float(coordinate[index])


def packet_geometry(x, period, weight, proper_q):
    """Canonical nodal weight on the proper-arc phase.

    ``weight`` is the probability element, or a declared shell-energy element.
    It is not multiplied by q. The phase uses ds = q dx.
    """
    samples = np.asarray(weight, dtype=float)
    coordinate = np.asarray(x, dtype=float)
    if samples.shape != coordinate.shape:
        raise ValueError("packet weight must share the nodes of x")
    if np.any(samples < 0.0) or not np.isfinite(samples).all():
        raise ValueError("packet measure must be nonnegative")
    total = float(np.sum(samples))
    phase, proper_length = _proper_phase(coordinate, proper_q, period)
    if total <= 0.0:
        return {
            "resultant": 0.0,
            "second_resultant": 0.0,
            "arc_location": None,
            "bridge_location": None,
            "proper_phase": None,
            "location_kind": "unresolved",
            "mass": total,
            "occupation_concentration": 0.0,
            "proper_length": proper_length,
            "q_multiplied_into_measure": False,
        }
    first, first_axis = _phase_moment(phase, samples, 1)
    second, second_axis = _phase_moment(phase, samples, 2)
    if first > NUMERICAL_RESULTANT_FLOOR and first >= second:
        chosen = first_axis
        kind = "proper_arc_first_moment"
        bridge_phase = (chosen + np.pi) % (2.0 * np.pi)
    elif second > NUMERICAL_RESULTANT_FLOOR:
        chosen = second_axis
        kind = "proper_arc_second_harmonic_axis"
        bridge_phase = (chosen + 0.5 * np.pi) % (2.0 * np.pi)
    else:
        chosen = None
        kind = "unresolved"
        bridge_phase = None
    arc = None if chosen is None else _phase_to_coordinate(chosen, phase, coordinate, period)
    bridge = None if bridge_phase is None else _phase_to_coordinate(
        bridge_phase, phase, coordinate, period,
    )
    fraction = samples / total
    concentration = float(samples.size * np.sum(fraction * fraction))
    return {
        "resultant": first,
        "second_resultant": second,
        "arc_location": arc,
        "bridge_location": bridge,
        "proper_phase": None if chosen is None else float(chosen),
        "location_kind": kind,
        "mass": total,
        "occupation_concentration": concentration,
        "proper_length": proper_length,
        "q_multiplied_into_measure": False,
    }


def _broadcast_space(array, n_t, n, name):
    if array is None:
        raise ValueError(f"missing {name}")
    values = np.asarray(array, dtype=float)
    if values.ndim == 0:
        return np.broadcast_to(values, (n_t, n)).copy()
    if values.shape == (n,):
        return np.broadcast_to(values, (n_t, n)).copy()
    if values.shape == (n_t, n):
        return values.copy()
    raise ValueError(f"{name} must be scalar, length {n}, or shape {(n_t, n)}")


def _as_windows(array, n_t, name):
    values = np.asarray(array, dtype=float)
    if values.ndim != 2 or values.shape[0] != n_t:
        raise ValueError(f"{name} must have shape {(n_t, 'n_windows')}")
    return values


def _window_midpoint(start, end, period):
    length = (float(end) - float(start)) % float(period)
    return (float(start) + 0.5 * length) % float(period)


def _nodes_in_window(x, start, end, period):
    delta = (np.asarray(x, dtype=float) - float(start)) % float(period)
    length = (float(end) - float(start)) % float(period)
    if length <= 0.0:
        return np.zeros(np.shape(x), dtype=bool)
    return delta < length


def _dominant_windows(weight, x, edges, period):
    """Windows that carry the leading nodal measure. Ties stay labeled, not averaged away."""
    masses = []
    for index in range(edges.shape[0]):
        mask = _nodes_in_window(x, edges[index, 0], edges[index, 1], period)
        masses.append(float(np.sum(np.asarray(weight, dtype=float)[mask])))
    if not masses:
        return []
    peak = max(masses)
    if peak <= 0.0:
        return []
    tolerance = max(1e-8 * peak, 1e-12)
    return [index for index, mass in enumerate(masses) if peak - mass <= tolerance]


def _trapz(values, times):
    return float(np.trapezoid(np.asarray(values, dtype=float), np.asarray(times, dtype=float)))


def _rejection(reasons, **fields):
    result = {
        "completed": False,
        "metric_completed": False,
        "goal_met": None,
        "closed_regime": False,
        "actual_bound": None,
        "constraint_consistent_solution": False,
        "reasons": list(reasons),
        "programme_thresholds_not_applied": dict(PROXY_THRESHOLDS_NOT_APPLIED),
        "proxy_agreement_is_metric_verification": False,
        "el_substitution_accepted": False,
        "coarse_frame_difference_used_as_rate": False,
        "manufactured_force_added": False,
    }
    result.update(fields)
    return result


def _binding_reasons(episode):
    reasons = []
    binding = episode.get("binding")
    if not isinstance(binding, dict):
        return ["missing_source_binding"]
    declared_job = binding.get("driver_job")
    episode_job = episode.get("driver_job")
    if not isinstance(declared_job, str) or not declared_job.strip() or episode_job not in (None, declared_job):
        reasons.append("driver_job_mismatch")
    for key, owner in OWNERS.items():
        if binding.get(f"{key}_owner") != owner:
            reasons.append(f"owner_mismatch_{key}")
    source = binding.get("source_id")
    if not isinstance(source, str) or not source.strip():
        reasons.append("missing_source_binding")
    if binding.get("phi_bound") is not True:
        reasons.append("phi_not_bound")
    return reasons


def _jet_reasons(jet):
    reasons = []
    if not isinstance(jet, dict):
        return ["missing_time_jet"]
    origin = jet.get("origin")
    if origin in REJECTED_JET_ORIGINS or jet.get("substitute_chi_for_curvature") is True:
        return ["el_substitution_rejected"]
    if origin not in ACCEPTED_JET_ORIGINS:
        reasons.append("jet_origin_rejected")
    if jet.get("projection_kept") is not True:
        reasons.append("projection_not_kept")
    projection = jet.get("projection")
    if not isinstance(projection, str) or not projection.strip():
        reasons.append("projection_not_kept")
    if "Q_dot" not in jet or jet.get("Q_dot") is None:
        reasons.append("missing_Q_dot")
    if "Q_dot_rate" not in jet or jet.get("Q_dot_rate") is None:
        reasons.append("missing_Q_dot_rate")
    neighbor = jet.get("neighbor")
    if not isinstance(neighbor, dict):
        reasons.append("missing_neighbor_increment")
        return reasons
    if neighbor.get("projection_kept") is not True:
        reasons.append("projection_not_kept")
    if neighbor.get("dt") is None or neighbor.get("Q_dot_next") is None:
        reasons.append("missing_neighbor_increment")
    if origin == "realized_increment" and neighbor.get("Q_next") is None:
        reasons.append("missing_neighbor_increment")
    return reasons


def _neighbor_agrees(jet, Q, Q_dot, Q_dot_rate):
    neighbor = jet["neighbor"]
    step = np.asarray(neighbor["dt"], dtype=float)
    if np.any(step <= 0.0) or not np.isfinite(step).all():
        return False
    rate_increment = (np.asarray(neighbor["Q_dot_next"], dtype=float) - Q_dot) / step
    if not np.allclose(rate_increment, Q_dot_rate, rtol=NEIGHBOR_AGREEMENT_RTOL, atol=NEIGHBOR_AGREEMENT_ATOL):
        return False
    if jet["origin"] == "realized_increment":
        increment = (np.asarray(neighbor["Q_next"], dtype=float) - Q) / step
        if not np.allclose(increment, Q_dot, rtol=NEIGHBOR_AGREEMENT_RTOL, atol=NEIGHBOR_AGREEMENT_ATOL):
            return False
    return True


def _indicators_from_pair(primary, other):
    feature_key = {"content": "content_change"}
    indicators = {}
    for key in INDICATOR_KEYS:
        source = feature_key.get(key, key)
        if key in ("arc_location", "bridge_location"):
            # Uncertainty of the displacement. Unequal grids do not share node labels.
            source = "arc_displacement" if key == "arc_location" else "bridge_displacement"
            left = primary.get(source)
            right = other.get(source)
            if left is None and right is None:
                indicators[key] = 0.0
            elif left is None or right is None:
                return None
            else:
                indicators[key] = abs(float(left) - float(right))
        else:
            if primary[source] is None or other[source] is None:
                return None
            indicators[key] = abs(float(primary[source]) - float(other[source]))
    return indicators


def _read_indicators(episode, paired):
    supplied = None
    comparison = episode.get("comparison")
    if isinstance(comparison, dict) and isinstance(comparison.get("refinement_indicator"), dict):
        raw = comparison["refinement_indicator"]
        if all(key in raw and np.isfinite(float(raw[key])) and float(raw[key]) >= 0.0 for key in INDICATOR_KEYS):
            supplied = {key: float(raw[key]) for key in INDICATOR_KEYS}
    if paired is None:
        return supplied
    merged = dict(paired)
    if supplied is not None:
        for key in INDICATOR_KEYS:
            merged[key] = max(merged[key], supplied[key])
    return merged


def _features(episode):
    """Metric, packet, and ledger scalars. Raises ValueError on a bad jet match."""
    chart = episode["chart"]
    jet = episode["time_jet"]
    ledger = episode["ledger"]
    times = np.asarray(chart["times"], dtype=float)
    coordinate = np.asarray(chart["x"], dtype=float)
    n_t = int(times.size)
    n = int(coordinate.size)
    ledger_times = times
    if ledger.get("times") is not None:
        ledger_times = np.asarray(ledger["times"], dtype=float)
        if ledger_times.ndim != 1 or ledger_times.size < 2 or np.any(np.diff(ledger_times) <= 0.0):
            raise ValueError("ledger times must be strictly increasing")
    period = float(chart["period"])
    if times.ndim != 1 or n_t < 2 or np.any(np.diff(times) <= 0.0):
        raise ValueError("times must be strictly increasing")
    gauge_hold = chart.get("gauge_hold") is True
    if gauge_hold:
        L_dot = np.zeros((n_t, n))
        beta_dot = np.zeros((n_t, n))
        gauge_time_derivative = "declared_hold"
    else:
        if "L_dot" not in chart or "beta_dot" not in chart:
            raise ValueError("missing_gauge_time_derivative")
        L_dot = _broadcast_space(chart["L_dot"], n_t, n, "L_dot")
        beta_dot = _broadcast_space(chart["beta_dot"], n_t, n, "beta_dot")
        gauge_time_derivative = "stored"
    L = _broadcast_space(chart["L"], n_t, n, "L")
    beta = _broadcast_space(chart["beta"], n_t, n, "beta")
    Q = _broadcast_space(chart["Q"], n_t, n, "Q")
    r = _broadcast_space(chart["r"], n_t, n, "r")
    Q_dot = _broadcast_space(jet["Q_dot"], n_t, n, "Q_dot")
    Q_dot_rate = _broadcast_space(jet["Q_dot_rate"], n_t, n, "Q_dot_rate")
    if not _neighbor_agrees(jet, Q, Q_dot, Q_dot_rate):
        raise ValueError("neighbor_increment_disagrees")
    if not _positive_chart_mask(L, Q, r):
        raise ValueError("positive_chart_failed")
    metric = direct_rh_grid(
        L, Q, beta, Q_dot, Q_dot_rate, period, L_dot=L_dot, beta_dot=beta_dot,
    )
    clocks, physical_q = proper_clocks(r, L, Q)
    weight = np.asarray(ledger["packet_weight"], dtype=float)
    if weight.shape != (n_t, n) or np.any(weight < 0.0) or not np.isfinite(weight).all():
        raise ValueError("missing_spatial_packet")
    locations = [
        packet_geometry(coordinate, period, weight[index], physical_q[index])
        for index in range(n_t)
    ]
    edges = np.asarray(ledger["window_edges"], dtype=float)
    ledger_count = int(ledger_times.size)
    energy = _as_windows(ledger["normal_energy"], ledger_count, "normal_energy")
    flux = _as_windows(ledger["flux_budget_term"], ledger_count, "flux_budget_term")
    pressure = _as_windows(ledger["pressure_work"], ledger_count, "pressure_work")
    lapse_work = _as_windows(ledger["lapse_exchange"], ledger_count, "lapse_exchange")
    if edges.shape != (energy.shape[1], 2):
        raise ValueError("window_edges do not match the ledger")
    if not all(np.isfinite(array).all() for array in (energy, flux, pressure, lapse_work, edges)):
        raise ValueError("nonfinite_fields")
    content_rate = np.gradient(energy, ledger_times, axis=0)
    predicted = flux + pressure + lapse_work
    residual = content_rate - predicted
    integrated_residual = (energy[-1] - energy[0]) - np.trapezoid(predicted, ledger_times, axis=0)
    regions = _dominant_windows(weight[0], coordinate, edges, period)
    if not regions:
        transferred = None
        work = None
        content_scale = None
        region_midpoint = None
        content_change = None
        flux_end = None
        rate_end = None
    else:
        selected = np.asarray(regions, dtype=int)
        flux_series = np.sum(flux[:, selected], axis=1)
        energy_series = np.sum(energy[:, selected], axis=1)
        work_series = np.sum(pressure[:, selected] + lapse_work[:, selected], axis=1)
        rate_series = np.sum(content_rate[:, selected], axis=1)
        transferred = _trapz(flux_series, ledger_times)
        work = _trapz(work_series, ledger_times)
        content_scale = abs(_trapz(energy_series, ledger_times) / (ledger_times[-1] - ledger_times[0]))
        content_change = float(energy_series[-1] - energy_series[0])
        flux_end = float(flux_series[-1])
        rate_end = float(rate_series[-1])
        midpoints = [
            _window_midpoint(edges[index, 0], edges[index, 1], period) for index in regions
        ]
        region_midpoint = float(min(midpoints))
    weighted_lapse = np.sum(weight * clocks, axis=1) / np.maximum(np.sum(weight, axis=1), 1e-30)
    proper_time = _trapz(weighted_lapse, times)
    if ledger.get("proper_time_elapsed") is not None:
        proper_time = float(ledger["proper_time_elapsed"])
    duration = float(ledger_times[-1] - ledger_times[0])
    ratio = None
    if transferred is not None and content_scale not in (None, 0.0) and proper_time > 0.0:
        ratio = abs(transferred) / (content_scale * proper_time)
    end = locations[-1]
    displacement = None
    bridge_displacement = None
    if locations[0]["arc_location"] is not None and end["arc_location"] is not None:
        displacement = circular_distance(end["arc_location"], locations[0]["arc_location"], period)
        bridge_displacement = circular_distance(
            end["bridge_location"], locations[0]["bridge_location"], period,
        )
    resultant = float(min(item["resultant"] for item in locations))
    second_resultant = float(min(item["second_resultant"] for item in locations))
    occupation_concentration = float(min(item["occupation_concentration"] for item in locations))
    if resultant > NUMERICAL_RESULTANT_FLOOR:
        location_width = period * float(np.sqrt(-2.0 * np.log(min(resultant, 1.0 - 1e-15)))) / (2.0 * np.pi)
    else:
        location_width = period / max(float(np.sqrt(max(occupation_concentration, 1e-30))), 1e-8)
    packet_R = [
        float(np.sum(weight[index] * metric["R_h"][index]) / max(float(np.sum(weight[index])), 1e-30))
        for index in range(n_t)
    ]
    chi = chart.get("chi")
    proxy_gap = None
    proxy = None
    actual = weyl_actual(metric["R_h"], r)
    if chi is not None:
        proxy = weyl_proxy_from_chi(_broadcast_space(chi, n_t, n, "chi"), r)
        proxy_gap = actual - proxy
    basis = ledger.get("phi_basis_occupations")
    basis_peak = None
    if basis is not None:
        occupations = np.asarray(basis, dtype=float)
        if occupations.ndim == 2:
            occupations = occupations[-1]
        if occupations.size == n:
            raise ValueError("basis_occupations_have_spatial_length")
        basis_peak = int(np.argmax(np.abs(occupations)))
    return {
        "period": period,
        "times": times,
        "n_t": n_t,
        "duration": duration,
        "proper_time": proper_time,
        "gauge_time_derivative": gauge_time_derivative,
        "K": metric["K"],
        "R_h": metric["R_h"],
        "R_h_start": packet_R[0],
        "R_h_end": packet_R[-1],
        "R_h_change": abs(packet_R[-1] - packet_R[0]),
        "weyl_actual_max": float(np.max(np.abs(actual))),
        "weyl_proxy_max": None if proxy is None else float(np.max(np.abs(proxy))),
        "proxy_gap_max": None if proxy_gap is None else float(np.max(np.abs(proxy_gap))),
        "arc_location": None if end["arc_location"] is None else float(end["arc_location"]),
        "bridge_location": None if end["bridge_location"] is None else float(end["bridge_location"]),
        "arc_start": locations[0]["arc_location"],
        "bridge_start": locations[0]["bridge_location"],
        "arc_displacement": displacement,
        "bridge_displacement": bridge_displacement,
        "packet_resultant": resultant,
        "second_resultant": second_resultant,
        "occupation_concentration": occupation_concentration,
        "occupation_excess": occupation_concentration - 1.0,
        "measure_mass": float(locations[0]["mass"]),
        "measure_kind": ledger.get("packet_weight_kind"),
        "measure_is_mode_probability": ledger.get("packet_weight_kind") == CANONICAL_HALF_DENSITY,
        "shell_energy_is_mode_probability": False,
        "q_multiplied_into_measure": False,
        "location_kind": end["location_kind"],
        "proper_phase": end["proper_phase"],
        "proper_length": float(locations[0]["proper_length"]),
        "location_width": location_width,
        "region_count": len(regions),
        "region_midpoint": region_midpoint,
        "transferred_energy": transferred,
        "content_scale": content_scale,
        "content_change": content_change,
        "throughput_ratio": ratio,
        "work_integral": work,
        "budget_residual_max": float(np.max(np.abs(integrated_residual))),
        "budget_pointwise_residual_max": float(np.max(np.abs(residual))),
        "budget_residual_ratio": (
            float(np.max(np.abs(integrated_residual))) / max(content_scale * proper_time, 1e-30)
            if content_scale not in (None, 0.0) and proper_time > 0.0
            else None
        ),
        "flux_budget_end": flux_end,
        "content_rate_end": rate_end,
        "arrived_on_initial_bridge": (
            displacement is not None
            and locations[0]["bridge_location"] is not None
            and circular_distance(end["arc_location"], locations[0]["bridge_location"], period)
            <= 1e-8
        ),
        "basis_peak": basis_peak,
        "metric_from_realized_jet": True,
        "Q_max": float(np.max(Q)),
        "chi_max": None if chi is None else float(np.max(np.abs(_broadcast_space(chi, n_t, n, "chi")))),
        "input_label": None if not isinstance(episode.get("comparison"), dict) else episode["comparison"].get("label"),
    }


def _ledger_reasons(episode, n_t):
    reasons = []
    ledger = episode.get("ledger")
    if not isinstance(ledger, dict):
        return ["incomplete_ledger"]
    required = (
        "normal_energy", "flux_budget_term", "pressure_work", "lapse_exchange",
        "window_edges", "packet_weight", "source_id",
    )
    if any(name not in ledger or ledger.get(name) is None for name in required):
        reasons.append("incomplete_ledger")
    if ledger.get("flux_is_budget_term") is not True:
        reasons.append("flux_sign_ambiguous")
    binding = episode.get("binding") if isinstance(episode.get("binding"), dict) else {}
    if ledger.get("source_id") != binding.get("source_id"):
        reasons.append("source_id_mismatch")
    if "packet_weight" not in ledger or ledger.get("packet_weight") is None:
        reasons.append("missing_spatial_packet")
    if ledger.get("packet_weight_kind") not in ADMISSIBLE_PACKET_MEASURES:
        reasons.append("packet_measure_rejected")
    if "phi_basis_occupations" in ledger and "packet_weight" not in ledger:
        reasons.append("missing_spatial_packet")
    del n_t
    return reasons


def _prepare(episode):
    if not isinstance(episode, dict):
        return _rejection(["missing_episode"])
    if episode.get("manufactured_reset") or episode.get("added_force") or episode.get("pump"):
        return _rejection(["manufactured_force_rejected"])
    reasons = _binding_reasons(episode)
    reasons.extend(_jet_reasons(episode.get("time_jet")))
    chart = episode.get("chart")
    if not isinstance(chart, dict):
        reasons.append("positive_chart_failed")
    else:
        for name in ("period", "x", "times", "L", "beta", "Q", "r"):
            if name not in chart:
                reasons.append("positive_chart_failed")
                break
    n_t = 0
    if isinstance(chart, dict) and "times" in chart:
        n_t = int(np.asarray(chart["times"]).size)
    reasons.extend(_ledger_reasons(episode, n_t))
    # An EL rejection stops before any curvature is formed.
    if "el_substitution_rejected" in reasons:
        return _rejection(["el_substitution_rejected"])
    blocking = [reason for reason in reasons if reason != "two_cut_drift"]
    if blocking:
        return _rejection(blocking)
    try:
        features = _features(episode)
    except ValueError as exc:
        return _rejection([str(exc)])
    return {"ok": True, "features": features}


def _bound_decision(episode):
    """A mathematical bound is a supplied domination, not a refinement difference."""
    declared = episode.get("mathematical_bound")
    if not isinstance(declared, dict):
        return {
            "actual_bound": None,
            "closed_regime": False,
            "bound_status": "not_supplied",
        }
    value = declared.get("value")
    name = str(declared.get("name") or "mathematical_bound")
    if value is None:
        return {
            "actual_bound": None,
            "closed_regime": False,
            "bound_status": name,
        }
    number = float(value)
    proven = declared.get("proven") is True and np.isfinite(number)
    return {
        "actual_bound": number if proven else None,
        "closed_regime": bool(proven),
        "bound_status": name if proven else "unproven",
    }


def assess(episode, comparison_episode=None):
    """Decide the stored episode. A missing rate or source binding is not a zero."""
    prepared = _prepare(episode)
    if not prepared.get("ok"):
        return prepared
    features = prepared["features"]
    paired = None
    if comparison_episode is not None:
        other = _prepare(comparison_episode)
        if not other.get("ok"):
            result = _rejection(
                ["comparison_episode_rejected"],
                metric_completed=True,
                metric=_metric_public(features),
            )
            return result
        paired = _indicators_from_pair(features, other["features"])
        if paired is None:
            return _rejection(["missing_comparison_uncertainty"], metric_completed=True, metric=_metric_public(features))
    indicators = _read_indicators(episode, paired)
    metric = _metric_public(features)
    if indicators is None:
        return _rejection(
            ["missing_comparison_uncertainty"],
            metric_completed=True,
            metric=metric,
            packet=_packet_public(features),
            ledger=_ledger_public(features),
        )
    move = features["arc_displacement"]
    bridge_move = features["bridge_displacement"]
    # A relocation has to leave the packet's own proper width. A smaller drift
    # stays the same structure. The width is this sample's circular scale.
    resolved_move = (
        move is not None
        and bridge_move is not None
        and (
            (move > indicators["arc_location"] and move > features["location_width"])
            or (bridge_move > indicators["bridge_location"] and bridge_move > features["location_width"])
        )
    )
    localized = (
        features["occupation_excess"]
        > indicators["occupation_concentration"] + NUMERICAL_RESULTANT_FLOOR
    )
    residual_ratio = features["budget_residual_ratio"]
    if residual_ratio is None:
        residual_ratio = indicators["throughput_ratio"] + 1.0
    throughput_uncertainty = max(indicators["throughput_ratio"], residual_ratio)
    throughput = (
        features["throughput_ratio"] is not None
        and features["throughput_ratio"] > throughput_uncertainty
    )
    work_resolved = (
        features["work_integral"] is not None
        and abs(features["work_integral"]) > indicators["work_integral"]
    )
    meaningful = (
        features["n_t"] >= 3
        and features["duration"] > 0.0
        and features["proper_time"] > 0.0
        and episode.get("arbitrary_cut") is not True
        and (throughput or resolved_move)
    )
    bounded = (
        features["proper_time"] > 0.0
        and np.isfinite(features["duration"])
        and np.isfinite(features["R_h"]).all()
        and np.isfinite(features["K"]).all()
    )
    location_known = (
        indicators["arc_location"] <= features["location_width"]
        and indicators["bridge_location"] <= features["location_width"]
    )
    renewed = bool(
        meaningful and bounded and localized and resolved_move and throughput and work_resolved
    )
    maintained = bool(
        meaningful
        and bounded
        and localized
        and throughput
        and work_resolved
        and not resolved_move
        and location_known
    )
    goal = bool(renewed or maintained)
    bound = _bound_decision(episode)
    return {
        "completed": True,
        "metric_completed": True,
        "goal_met": goal,
        "renewed_structure": renewed,
        "maintained_structure": maintained,
        "localized": bool(localized),
        "homogeneous": not bool(localized),
        "single_harmonic_resultant_is_localization": False,
        "throughflow": bool(throughput),
        "coupled_feedback": bool(work_resolved and features["metric_from_realized_jet"]),
        "meaningful_transfer_episode": bool(meaningful),
        "arbitrary_two_cut": features["n_t"] < 3 or episode.get("arbitrary_cut") is True,
        "positive_chart": True,
        "admissible_fields": True,
        "complete_ledger": True,
        "bounded_event": bool(bounded and meaningful),
        "flux_is_content_rate": False,
        "closed_regime": bound["closed_regime"],
        "actual_bound": bound["actual_bound"],
        "bound_status": bound["bound_status"],
        "constraint_consistent_solution": False,
        "refinement_indicator": indicators,
        "actual_bound_accepted": bound["closed_regime"],
        "reasons": [],
        "programme_thresholds_not_applied": dict(PROXY_THRESHOLDS_NOT_APPLIED),
        "proxy_agreement_is_metric_verification": False,
        "el_substitution_accepted": False,
        "coarse_frame_difference_used_as_rate": False,
        "manufactured_force_added": False,
        "spatial_arc_uses_basis": False,
        "global_Q_or_chi_max_used": False,
        "input_label": features["input_label"],
        "metric": metric,
        "packet": _packet_public(features),
        "ledger": _ledger_public(features),
        "comparison_uncertainty": {
            "source": "paired_episode" if paired is not None else "stored_indicator",
            "indicator": indicators,
            "actual_bound": bound["actual_bound"],
        },
    }


def _metric_public(features):
    return {
        "formula_K": "(Q_dot - d_x(beta Q)) / (L Q)",
        "formula_R_h": "2/(L Q) d_x(L_x/Q) - 2*(((d_t K) - beta d_x K)/L + K^2)",
        "christoffel_contraction": (
            "R_bd = d_a Gamma^a_db - d_d Gamma^a_ab "
            "+ Gamma^a_ae Gamma^e_db - Gamma^a_de Gamma^e_ab"
        ),
        "K_max": float(np.max(np.abs(features["K"]))),
        "R_h_mean": float(np.mean(features["R_h"])),
        "R_h_start": features["R_h_start"],
        "R_h_end": features["R_h_end"],
        "R_h_change": features["R_h_change"],
        "weyl_actual_max": features["weyl_actual_max"],
        "weyl_proxy_max": features["weyl_proxy_max"],
        "proxy_gap_max": features["proxy_gap_max"],
        "proxy_is_not_metric": True,
        "N": "r L",
        "q": "r Q",
        "gauge_time_derivative": features["gauge_time_derivative"],
        "metric_from_realized_jet": True,
    }


def _packet_public(features):
    return {
        "arc_start": features["arc_start"],
        "arc_end": features["arc_location"],
        "bridge_start": features["bridge_start"],
        "bridge_end": features["bridge_location"],
        "arc_displacement": features["arc_displacement"],
        "bridge_displacement": features["bridge_displacement"],
        "resultant": features["packet_resultant"],
        "second_resultant": features["second_resultant"],
        "occupation_concentration": features["occupation_concentration"],
        "occupation_excess": features["occupation_excess"],
        "measure_mass": features["measure_mass"],
        "measure_kind": features["measure_kind"],
        "measure_is_mode_probability": features["measure_is_mode_probability"],
        "shell_energy_is_mode_probability": False,
        "q_multiplied_into_measure": False,
        "location_kind": features["location_kind"],
        "proper_phase": features["proper_phase"],
        "proper_length": features["proper_length"],
        "region_count": features["region_count"],
        "location_width": features["location_width"],
        "arrived_on_initial_bridge": features["arrived_on_initial_bridge"],
        "region_midpoint": features["region_midpoint"],
        "basis_peak": features["basis_peak"],
        "proper_time": features["proper_time"],
        "coordinate_duration": features["duration"],
        "Q_max_diagnostic": features["Q_max"],
        "chi_max_diagnostic": features["chi_max"],
    }


def _ledger_public(features):
    return {
        "content_rate_end": features["content_rate_end"],
        "flux_budget_end": features["flux_budget_end"],
        "flux_is_content_rate": False,
        "transferred_energy": features["transferred_energy"],
        "content_scale": features["content_scale"],
        "content_change": features["content_change"],
        "throughput_ratio": features["throughput_ratio"],
        "work_integral": features["work_integral"],
        "pressure_and_lapse_included": True,
        "budget_residual_max": features["budget_residual_max"],
        "proper_time": features["proper_time"],
        "coordinate_duration": features["duration"],
    }


def _npy_header(blob):
    if not blob.startswith(b"\x93NUMPY"):
        raise ValueError("stored array is not npy")
    major = blob[6]
    if major == 1:
        header_len = int.from_bytes(blob[8:10], "little")
        start = 10
    else:
        header_len = int.from_bytes(blob[8:12], "little")
        start = 12
    header = blob[start:start + header_len].decode("latin1").strip()
    meta = ast.literal_eval(header)
    return tuple(meta["shape"])


def inspect_saved_episode(path=None):
    """Say whether the immutable episode payload can complete this assessment.

    Frame differences are not promoted to a realized rate, and missing rates
    are not replaced by zero. Endpoint ``Q_dot`` without a neighbor increment
    of that rate is not enough.
    """
    source = DEFAULT_EPISODE_NPZ if path is None else Path(path)
    report = {
        "path": str(source),
        "present": source.is_file(),
        "completed": False,
        "metric_completed": False,
        "goal_met": None,
        "closed_regime": False,
        "actual_bound": None,
        "stored_frame_weyl_is_chi_proxy": True,
        "global_Q_or_chi_max_used": False,
        "coarse_frame_difference_used_as_rate": False,
        "missing_rates_set_to_zero": False,
    }
    if not source.is_file():
        report["missing_primitive"] = ["episode_payload"]
        return report
    with zipfile.ZipFile(source) as archive:
        names = set(archive.namelist())

        def shape_of(key):
            if key not in names:
                return None
            return _npy_header(archive.read(key)[:256])

        frame_Q = shape_of(f"{_SAVED_RUN}_frame_Q.npy")
        windows = shape_of(f"{_SAVED_RUN}_window_normal.npy")
        endpoint_rate = shape_of(f"{_SAVED_RUN}_final_rate_Q_dot.npy")
        endpoint_p_chi = shape_of(f"{_SAVED_RUN}_final_p_chi.npy")
    available_frames = [
        name for name in _FRAME_FIELDS_PRESENT
        if f"{_SAVED_RUN}_frame_{name}.npy" in names
    ]
    missing = [
        "frame Q_dot",
        "frame Q_dot_rate",
        "neighbor increment of Q_dot",
        "frame p_chi",
        "frame L",
        "frame beta",
        "projection_kept on the time jet",
        "packet arc and periodic bridge series",
    ]
    report.update({
        "run": _SAVED_RUN,
        "available": {
            "frame_fields": available_frames,
            "frame_Q_shape": frame_Q,
            "window_normal_shape": windows,
            "window_proper_flux": f"{_SAVED_RUN}_window_proper_flux.npy" in names,
            "window_proper_work": f"{_SAVED_RUN}_window_proper_work.npy" in names,
            "window_lapse_work": f"{_SAVED_RUN}_window_lapse_work.npy" in names,
            "endpoint_Q_dot_shape": endpoint_rate,
            "endpoint_p_chi_shape": endpoint_p_chi,
            "endpoint_rate_is_a_single_slice": True,
            "endpoint_rate_matches_frame_nodes": (
                frame_Q is not None
                and endpoint_rate is not None
                and frame_Q[-1] == endpoint_rate[-1]
            ),
        },
        "missing_primitive": missing,
        "reason": "missing_realized_rate",
        "endpoint_Q_dot_without_rate_increment_is_not_enough": True,
    })
    return report


def _regeneration_geometry_map(fermions):
    """Owned maps. Quadrature Q is A_g times the coarse state; Phi uses U_f."""
    from recursive_horizons.nsc_spherical_coupling import PERIOD
    from recursive_horizons.nsc_spherical_galerkin_coupling import (
        antiperiodic_interpolation,
        canonical_column_map,
        periodic_interpolation,
    )
    from recursive_horizons.nsc_spherical_null_expansion import static_clock

    nf = int(fermions)
    ng = nf - 1
    nq = 4 * nf
    coordinate, lapse, shift = static_clock(nq, PERIOD)
    return {
        "nf": nf,
        "ng": ng,
        "nq": nq,
        "period": float(PERIOD),
        "x": coordinate,
        "L": lapse,
        "beta": shift,
        "A_g": periodic_interpolation(ng, nq, PERIOD),
        "U_f": canonical_column_map(antiperiodic_interpolation(nf, nq, PERIOD), nf, nq),
    }


def _nodal_probability(phi0, phi1, occupations, columns):
    frames = []
    fermion = []
    for index in range(phi0.shape[0]):
        coarse = np.sum(
            (np.abs(phi0[index]) ** 2 + np.abs(phi1[index]) ** 2) * occupations,
            axis=-1,
        )
        fermion.append(coarse)
        lifted0 = columns @ phi0[index]
        lifted1 = columns @ phi1[index]
        frames.append(np.sum((np.abs(lifted0) ** 2 + np.abs(lifted1) ** 2) * occupations, axis=-1))
    return np.stack(fermion), np.stack(frames)


def _case_episode(arrays, prefix, geometry, occupations, protocol):
    from recursive_horizons.nsc_regeneration_controls import CAR_GAP_MAX

    def take(name):
        key = prefix + name
        if key not in arrays:
            raise KeyError(key)
        return arrays[key]

    present = np.asarray(take("frame_increment_present"), dtype=bool)
    if not bool(np.all(present)):
        raise ValueError("frame_increment_present is false; the missing neighbor is not filled with zero")
    indicator = np.asarray(take("frame_indicator_Q_dot"), dtype=float)
    increment_dt = np.asarray(take("frame_increment_dt"), dtype=float)
    if not np.isfinite(indicator).all() or not np.isfinite(increment_dt).all() or np.any(increment_dt <= 0.0):
        raise ValueError("frame_indicator_Q_dot or frame_increment_dt is not a finite positive jet")
    coarse_q = np.asarray(take("frame_coarse_Q"), dtype=float)
    quad_q = np.asarray(take("frame_quad_Q"), dtype=float)
    prolonged_q = np.stack([geometry["A_g"] @ coarse_q[index] for index in range(coarse_q.shape[0])])
    quad_gap = float(np.max(np.abs(quad_q - prolonged_q)))
    if quad_gap > 1e-8:
        raise ValueError(
            f"frame_quad_Q is not the owned prolongation of frame_coarse_Q (gap {quad_gap})"
        )
    q_dot = np.stack([geometry["A_g"] @ take("frame_Q_dot")[index] for index in range(coarse_q.shape[0])])
    increment = np.asarray(take("frame_increment_Q_dot"), dtype=float)
    q_dot_rate = np.stack([geometry["A_g"] @ indicator[index] for index in range(indicator.shape[0])])
    q_dot_next = np.stack([
        geometry["A_g"] @ (take("frame_Q_dot")[index] + increment[index])
        for index in range(increment.shape[0])
    ])
    phi0 = take("frame_phi0")
    phi1 = take("frame_phi1")
    fermion_weight, probability = _nodal_probability(phi0, phi1, occupations, geometry["U_f"])
    fermion_mass = np.sum(fermion_weight, axis=1)
    if float(np.max(np.abs(fermion_mass - float(np.sum(occupations))))) > 1e-8:
        raise ValueError("fermion-grid |Phi|^2 does not sum to the six-mode occupation")
    gram_gap = 0.0
    for index in range(phi0.shape[0]):
        gram = phi0[index].conj().T @ phi0[index] + phi1[index].conj().T @ phi1[index]
        gram_gap = max(gram_gap, float(np.max(np.abs(gram - np.eye(gram.shape[0])))))
    times = np.asarray(take("frame_time"), dtype=float)
    ledger_times = np.asarray(take("time"), dtype=float)
    episode = {
        "driver_job": protocol,
        "mathematical_bound": {
            "name": "propagated_initial_bound_null",
            "value": None,
            "proven": False,
        },
        "binding": {
            "action_owner": OWNERS["action"],
            "geometry_owner": OWNERS["geometry"],
            "observer_owner": OWNERS["observer"],
            "energy_ledger_owner": OWNERS["energy_ledger"],
            "source_id": protocol,
            "phi_bound": True,
            "driver_job": protocol,
        },
        "chart": {
            "period": geometry["period"],
            "x": geometry["x"],
            "times": times,
            "L": geometry["L"],
            "beta": geometry["beta"],
            "Q": quad_q,
            "r": np.asarray(take("frame_quad_r"), dtype=float),
            "chi": np.asarray(take("frame_quad_chi"), dtype=float),
            "gauge_hold": True,
        },
        "time_jet": {
            "origin": "finite_ode_derivative",
            "projection_kept": True,
            "projection": "galerkin_pullback",
            "Q_dot": q_dot,
            "Q_dot_rate": q_dot_rate,
            "neighbor": {
                "dt": increment_dt.reshape((-1, 1)),
                "Q_dot_next": q_dot_next,
                "projection_kept": True,
            },
        },
        "ledger": {
            "source_id": protocol,
            "flux_is_budget_term": True,
            "times": ledger_times,
            "normal_energy": np.asarray(take("window_normal"), dtype=float),
            "flux_budget_term": np.asarray(take("observer_boundary"), dtype=float),
            "pressure_work": np.asarray(take("observer_pressure"), dtype=float),
            "lapse_exchange": np.asarray(take("observer_lapse"), dtype=float),
            "window_edges": REGENERATION_WINDOWS.copy(),
            "packet_weight_kind": CANONICAL_HALF_DENSITY,
            "packet_weight": probability,
            "phi_basis_occupations": np.asarray(occupations, dtype=float),
            "proper_time_elapsed": float(take("proper_clock")[-1]),
        },
        "comparison": {"label": prefix.rstrip("_")},
    }
    return episode, {
        "fermion_occupation": [float(value) for value in fermion_mass],
        "quadrature_occupation": [float(np.sum(row)) for row in probability],
        "gram_gap": gram_gap,
        "gram_admissible": bool(gram_gap <= float(CAR_GAP_MAX)),
        "gram_margin": float(CAR_GAP_MAX),
        "quad_prolongation_gap": quad_gap,
        "rate_source": "galerkin_pullback_neighbor_increment",
        "el_p_chi_dot_used": False,
        "packet_flux_end": float(take("packet_flux")[-1]),
        "reservoir_flux_end": float(take("reservoir_flux")[-1]),
        "leader_share_start": float(take("leader_share")[0]),
        "leader_share_end": float(take("leader_share")[-1]),
        "leader_content_change": float(take("leader_content")[-1] - take("leader_content")[0]),
        "leader_clock": float(take("leader_clock")[-1]),
        "proper_clock": float(take("proper_clock")[-1]),
        "positive_chart_series": bool(
            np.all(take("positive_L") > 0.0)
            and np.all(take("positive_Q") > 0.0)
            and np.all(take("positive_r") > 0.0)
        ),
    }


def _pair_run(name):
    other = {
        "nf256_dtmax_0_0005": "nf256_dtmax_0_00025",
        "nf256_dtmax_0_00025": "nf256_dtmax_0_0005",
        "nf512_dtmax_0_0005": "nf512_dtmax_0_00025",
        "nf512_dtmax_0_00025": "nf512_dtmax_0_0005",
    }
    return other[name]


def assess_saved_episode(path=None, record_path=None):
    """Read the completed four-case regeneration episode and assess it.

    The NPZ and JSON are not rewritten. A missing jet is not replaced by zero.
    """
    from recursive_horizons.nsc_spherical_coupling import OCCUPATIONS

    npz_path = REGENERATION_NPZ if path is None else Path(path)
    json_path = REGENERATION_JSON if record_path is None else Path(record_path)
    if not npz_path.is_file() or not json_path.is_file():
        return _rejection(["episode_payload"], metric_completed=False)
    record = json.loads(json_path.read_text())
    protocol = str(record.get("protocol_id") or "")
    if protocol != "nsc-regeneration-episode-v1":
        return _rejection([f"protocol_id:{protocol}"])
    arrays = np.load(npz_path)
    try:
        built = {}
        checks = {}
        maps = {}
        for name in REGENERATION_RUNS:
            fermions = 256 if name.startswith("nf256") else 512
            if fermions not in maps:
                maps[fermions] = _regeneration_geometry_map(fermions)
            episode, check = _case_episode(arrays, name + "_", maps[fermions], OCCUPATIONS, protocol)
            built[name] = episode
            checks[name] = check
    except (KeyError, ValueError) as exc:
        return _rejection([f"mapping:{exc}"], metric_completed=False)
    cases = {}
    for name in REGENERATION_RUNS:
        cases[name] = assess(built[name], comparison_episode=built[_pair_run(name)])
        cases[name]["record"] = checks[name]
    fine = cases["nf512_dtmax_0_00025"]
    spatial = assess(
        built["nf512_dtmax_0_00025"],
        comparison_episode=built["nf256_dtmax_0_00025"],
    )
    time_indicator = fine.get("refinement_indicator") or {}
    space_indicator = spatial.get("refinement_indicator") or {}
    movement = {}
    for key in INDICATOR_KEYS:
        if key in time_indicator and key in space_indicator:
            movement[key] = max(float(time_indicator[key]), float(space_indicator[key]))
    event = None
    runs = record.get("runs") if isinstance(record.get("runs"), dict) else {}
    fine_record = runs.get("nf512_dtmax_0_00025") if isinstance(runs, dict) else None
    if isinstance(fine_record, dict):
        frames = fine_record.get("expansion_frames") or []
        if frames:
            last = frames[-1]
            intervals = last.get("anti_trapped_intervals") or []
            if intervals:
                event = {
                    "stop_reason": fine_record.get("stop_reason"),
                    "bridge_arc_x_start": intervals[0].get("x_start"),
                    "packet_edge": PACKET_EDGE,
                    "entered_packet_edge": intervals[0].get("x_start") == PACKET_EDGE,
                }
    return {
        "protocol_id": protocol,
        "completed": all(bool(case.get("completed")) for case in cases.values()),
        "cases": cases,
        "primary": "nf512_dtmax_0_00025",
        "time_refinement_indicator": time_indicator,
        "spatial_refinement_indicator": space_indicator,
        "refinement_indicator": movement,
        "geometric_event": event,
        "actual_bound": None,
        "bound_status": "propagated_initial_bound_null",
        "closed_regime": False,
        "indicator_is_not_a_curvature_bound": True,
        "q_multiplied_into_measure": False,
        "el_substitution_accepted": False,
    }
