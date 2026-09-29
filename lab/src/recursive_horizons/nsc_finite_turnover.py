"""Finite channel turnover on the inherited nested window.

The stationary covariance, mode currents and reduced filter are calculated
here. Geometry is an input. This does not close the local incoming gate.
"""
from __future__ import annotations

import numpy as np
import sympy as sp

from .nsc_nested_qualities import finite_window

SCHEMA = "NSC-FINITE-CHANNEL-TURNOVER-v1"
PASS = "PASS_FINITE_CHANNEL_TURNOVER"
FAIL = "FAIL_FINITE_CHANNEL_TURNOVER"

_REQUIRED_METRICS = (
    "G_cut_exact",
    "total_inventory_exact",
    "sign_control_value_exact",
    "probe_response_change_fro",
    "omitted_memory_error_fro",
    "omitted_drive_error_fro",
    "dropped_cross_covariance_error_fro",
    "half_identity_response_change_fro",
    "far_probe_response_change_fro",
)


def _pair(value):
    value = sp.simplify(value)
    return [str(sp.together(sp.re(value))), str(sp.together(sp.im(value)))]


def _rational(value):
    value = sp.simplify(value)
    if value.is_integer:
        return str(sp.Integer(value))
    if not value.is_rational:
        raise ValueError(f"expected a rational value, got {value}")
    return str(sp.Rational(value))


def _entry_abs(value):
    value = sp.expand(value)
    modulus = sp.simplify(sp.sqrt(sp.simplify(sp.re(value) ** 2 + sp.im(value) ** 2)))
    if not modulus.is_rational:
        raise ValueError("window entry modulus is not rational")
    return sp.Rational(modulus)


def _onsite_projectors(local):
    root = sp.sqrt(29)
    larger = sp.Rational(3, 2) + root / 10
    smaller = sp.Rational(3, 2) - root / 10
    plus = sp.simplify((local - smaller * sp.eye(2)) / (larger - smaller))
    minus = sp.simplify(sp.eye(2) - plus)
    return larger, smaller, plus, minus


def _current(operator, state, region_a, proj_a, region_b, proj_b):
    coupling = operator[2 * region_a:2 * region_a + 2, 2 * region_b:2 * region_b + 2]
    coherence = state[2 * region_b:2 * region_b + 2, 2 * region_a:2 * region_a + 2]
    amplitude = sp.trace(proj_a * coupling * proj_b * coherence)
    return sp.simplify(2 * sp.im(amplitude))


def _algebraic_abs(value):
    value = sp.simplify(value)
    if value == 0:
        return sp.Integer(0)
    if value.is_positive:
        return value
    if value.is_negative:
        return sp.simplify(-value)
    raise ValueError(f"current has no exact sign: {value}")


def _as_numpy(matrix):
    return np.array(
        [[complex(matrix[i, j]) for j in range(matrix.cols)] for i in range(matrix.rows)],
        dtype=np.complex128,
    )


def _filtered(operator, state, spectral, *, memory=True, drive=True):
    """Reduced filter [S^{-1}, S^{-1} V (zI-J_EE)^{-1}]. No dagger on V."""
    retained = operator[:2, :2]
    exterior = operator[2:, 2:]
    coupling = operator[:2, 2:]
    resolvent = np.linalg.inv(spectral * np.eye(exterior.shape[0], dtype=np.complex128) - exterior)
    self_energy = coupling @ resolvent @ coupling.conj().T
    local = spectral * np.eye(2, dtype=np.complex128) - retained
    if memory:
        local = local - self_energy
    local_inverse = np.linalg.inv(local)
    if drive:
        drive_block = local_inverse @ coupling @ resolvent
    else:
        drive_block = np.zeros((2, exterior.shape[0]), dtype=np.complex128)
    amplitude = np.concatenate((local_inverse, drive_block), axis=1)
    return amplitude @ state @ amplitude.conj().T, retained, coupling, amplitude


def _frobenius(matrix):
    return float(np.linalg.norm(matrix, ord="fro"))


def compute_result():
    """Calculate the finite turnover witness. Checks are booleans; metrics are computed."""
    local = sp.Matrix([[1, sp.I / 5], [-sp.I / 5, 2]])
    link = sp.Matrix([[sp.Rational(1, 4), sp.I / 7], [sp.Rational(1, 9), sp.Rational(1, 6)]])
    omega = sp.Rational(3, 2)
    spectral = 1 + 2 * sp.I
    depth = 3
    bound = sp.Integer(8)
    window = finite_window(local, link, omega, 0, depth)
    denominator = 4 * bound ** 3
    covariance = sp.eye(window.rows) / 2 - window ** 3 / denominator
    larger, smaller, plus, minus = _onsite_projectors(local)
    projectors = {"plus": plus, "minus": minus}
    modes = [(region, name) for region in range(depth) for name in ("plus", "minus")]
    currents = {
        (left, right): _current(window, covariance, left[0], projectors[left[1]], right[0], projectors[right[1]])
        for left in modes
        for right in modes
    }
    row_sums = [sum(_entry_abs(window[row, col]) for col in range(window.cols)) for row in range(window.rows)]
    max_row_sum = max(row_sums)
    cut_terms = [
        currents[((0, left_mode), (region, right_mode))]
        for left_mode in ("plus", "minus")
        for region in (1, 2)
        for right_mode in ("plus", "minus")
    ]
    net_cut = sp.simplify(sum(cut_terms))
    exchange = sp.simplify(sum(_algebraic_abs(term) for term in cut_terms))
    one_way = sp.simplify(sum(term for term in cut_terms if term.is_positive))
    inventory = sp.simplify(sp.trace(covariance))
    energy = sp.simplify(sp.trace(covariance * window))

    sign_state = sp.zeros(6, 1)
    sign_state[0, 0] = 1 / sp.sqrt(2)
    sign_state[2, 0] = sp.I / sp.sqrt(2)
    sign_covariance = sp.simplify(sign_state * sign_state.H)
    region = sp.diag(1, 1, 0, 0, 0, 0)
    sign_derivative = sp.simplify(-sp.I * sp.trace(region * (window * sign_covariance - sign_covariance * window)))
    exterior = sp.eye(6) - region
    regional_current = sp.simplify(sp.trace(
        sign_covariance * (-sp.I * (region * window * exterior - exterior * window * region))
    ))
    incoming = sp.simplify(sum(
        _current(window, sign_covariance, 0, projectors[left_mode], region_index, projectors[right_mode])
        for left_mode in ("plus", "minus")
        for region_index in (1, 2)
        for right_mode in ("plus", "minus")
    ))

    shifted = finite_window(local, link, omega, 1, depth)
    shifted_bound = omega * bound
    shifted_covariance = sp.eye(shifted.rows) / 2 - shifted ** 3 / (4 * shifted_bound ** 3)
    shifted_currents = {
        (left, right): _current(shifted, shifted_covariance, left[0], projectors[left[1]], right[0], projectors[right[1]])
        for left in modes
        for right in modes
    }

    operator = _as_numpy(window)
    state = _as_numpy(covariance)
    spectral_value = complex(spectral)
    baseline, retained, coupling, amplitude = _filtered(operator, state, spectral_value)
    full_rows = np.linalg.inv(spectral_value * np.eye(6) - operator)[:2, :]
    memory_free, _, _, _ = _filtered(operator, state, spectral_value, memory=False, drive=True)
    drive_free, _, _, _ = _filtered(operator, state, spectral_value, memory=True, drive=False)
    blocked = state.copy()
    blocked[:2, 2:] = 0
    blocked[2:, :2] = 0
    blocked_covariance = amplitude @ blocked @ amplitude.conj().T
    half = amplitude @ (0.5 * np.eye(6, dtype=np.complex128)) @ amplitude.conj().T
    state_snapshot = state.copy()

    def _probe(index, delta):
        probed = sp.Matrix(window)
        probed[index, index] += delta
        unchanged = probed[:2, :2] == window[:2, :2] and probed[:2, 2:] == window[:2, 2:]
        probed_covariance, _, _, _ = _filtered(_as_numpy(probed), state, spectral_value)
        return unchanged, _frobenius(probed_covariance - baseline)

    near_fixed, near_change = _probe(2, sp.Rational(1, 8))
    far_fixed, far_change = _probe(5, sp.Rational(1, 10))
    shifted_operator = _as_numpy(shifted)
    base_resolvent = np.linalg.inv(spectral_value * np.eye(6) - operator)
    shifted_resolvent = np.linalg.inv(complex(omega * spectral) * np.eye(6) - shifted_operator)
    spectrum = np.linalg.eigvalsh((state + state.conj().T) / 2)

    checks = {
        "hermitian_window": bool(window == window.H),
        "dimension_six": bool(window.rows == 6 and window.cols == 6),
        "no_corner_block": bool(window[:2, 4:] == sp.zeros(2) and window[4:, :2] == sp.zeros(2, 2)),
        "link_j02_is_one_quarter": bool(sp.simplify(window[0, 2] - sp.Rational(1, 4)) == 0),
        "row_sum_bound_M8": bool(max_row_sum <= bound),
        "cubic_denominator_is_2048": bool(denominator == 2048),
        "commutator_zero": bool(sp.expand(window * covariance - covariance * window) == sp.zeros(6)),
        "onsite_projectors_resolve_identity": bool(sp.simplify(plus + minus - sp.eye(2)) == sp.zeros(2)),
        "onsite_projectors_match_eigenvalues": bool(
            sp.simplify(plus * local - larger * plus) == sp.zeros(2)
            and sp.simplify(minus * local - smaller * minus) == sp.zeros(2)
        ),
        "currents_antisymmetric": bool(all(
            sp.simplify(currents[(left, right)] + currents[(right, left)]) == 0
            for left in modes for right in modes
        )),
        "mode_divergence_zero": bool(all(
            sp.simplify(sum(currents[(mode, other)] for other in modes)) == 0
            for mode in modes
        )),
        "corner_currents_zero": bool(all(
            currents[((0, left_mode), (2, right_mode))] == 0
            for left_mode in ("plus", "minus")
            for right_mode in ("plus", "minus")
        )),
        "intraregion_mode_currents_zero": bool(all(
            currents[((region, "plus"), (region, "minus"))] == 0
            for region in range(depth)
        )),
        "interregion_exchange_positive": bool(exchange.is_positive),
        "cut_net_current_zero": bool(net_cut == 0),
        "one_way_exchange_is_half_cut": bool(sp.simplify(2 * one_way - exchange) == 0),
        "inventory_rational": bool(inventory.is_rational and inventory.is_positive),
        "energy_rational": bool(energy.is_real),
        "spectrum_inside_quarter_and_three_quarters": bool(spectrum.min() >= 0.25 - 1e-9 and spectrum.max() <= 0.75 + 1e-9),
        "full_filter_matches_reduced": bool(np.allclose(amplitude, full_rows, atol=1e-9, rtol=0)),
        "omitted_memory_changes_response": bool(_frobenius(memory_free - baseline) > 0),
        "omitted_drive_changes_response": bool(_frobenius(drive_free - baseline) > 0),
        "dropped_cross_changes_response": bool(_frobenius(blocked_covariance - baseline) > 0),
        "half_identity_changes_response": bool(_frobenius(half - baseline) > 0),
        "near_probe_changes_response": bool(near_change > 0),
        "near_probe_preserves_jaa_v": bool(near_fixed),
        "far_probe_changes_response": bool(far_change > 0),
        "far_probe_preserves_jaa_v": bool(far_fixed),
        "probes_keep_baseline_covariance": bool(np.array_equal(state, state_snapshot)),
        "sign_derivative_matches_regional_current": bool(sp.simplify(sign_derivative - regional_current) == 0),
        "sign_derivative_matches_incoming_mode_sum": bool(sp.simplify(sign_derivative - incoming) == 0),
        "sign_control_is_one_quarter": bool(sp.simplify(sign_derivative - sp.Rational(1, 4)) == 0),
        "shifted_window_is_omega_times_base": bool(sp.expand(shifted - omega * window) == sp.zeros(6)),
        "shifted_bound_is_12": bool(shifted_bound == 12),
        "shifted_covariance_matches_base": bool(sp.expand(shifted_covariance - covariance) == sp.zeros(6)),
        "shifted_currents_scale_by_omega": bool(all(
            sp.simplify(shifted_currents[pair] - omega * currents[pair]) == 0
            for pair in currents
        )),
        "normalized_resolvent_scales": bool(np.allclose(
            complex(omega) * shifted_resolvent, base_resolvent, atol=1e-8, rtol=0
        )),
    }
    metrics = {
        "G_cut_exact": _rational(exchange),
        "G_one_way_exact": _rational(one_way),
        "total_inventory_exact": _rational(inventory),
        "total_energy_exact": _rational(sp.re(energy)) if energy.is_real else str(energy),
        "max_row_sum_exact": _rational(max_row_sum),
        "sign_control_value_exact": _rational(sign_derivative),
        "probe_response_change_fro": near_change,
        "omitted_memory_error_fro": _frobenius(memory_free - baseline),
        "omitted_drive_error_fro": _frobenius(drive_free - baseline),
        "dropped_cross_covariance_error_fro": _frobenius(blocked_covariance - baseline),
        "half_identity_response_change_fro": _frobenius(half - baseline),
        "far_probe_response_change_fro": far_change,
        "covariance_eigenvalue_min": float(spectrum.min()),
        "covariance_eigenvalue_max": float(spectrum.max()),
    }
    missing = [name for name in _REQUIRED_METRICS if name not in metrics]
    if missing:
        raise RuntimeError(f"missing metrics: {missing}")
    assumptions = [
        "Finite coherent channel currents are evaluated on one fixed six-dimensional nested window.",
        "Geometry, Omega=3/2 and the path of three regions are prescribed inputs.",
        "The stationary state is the prescribed covariance C=I/2-J^3/2048, with M=8.",
        "Mode observables are the spectral projectors of the declared onsite H, embedded in each region.",
        "Mean-current activity is not particle production, heat, entropy production, dynamic-radius stabilization, an antimatter identification, or an eternity proof.",
        "The filtered matrix F C F^dagger is a resolvent-weighted covariance/response, not an equal-time occupation or a stress tensor.",
        "Operator self-energy is not identified with C. Probes keep J_AA, V and the baseline C fixed.",
        "The local incoming gate remains OPEN.",
    ]
    result = {
        "schema": SCHEMA,
        "verdict": PASS if all(checks.values()) else FAIL,
        "physical_local_gate": "OPEN",
        "inputs": {
            "H": [_pair(local[i, j]) for i in range(2) for j in range(2)],
            "B": [_pair(link[i, j]) for i in range(2) for j in range(2)],
            "omega": str(omega),
            "depth": depth,
            "first": 0,
            "spectral_parameter": "1+2*I",
            "bound_M": str(bound),
            "covariance_formula": "I/2 - J**3/(4*M**3)",
            "cubic_denominator": str(denominator),
            "dimension": 6,
            "near_probe": {"index": [2, 2], "delta": "1/8"},
            "far_probe": {"index": [5, 5], "delta": "1/10"},
            "sign_control": "(e0+I*e2)/sqrt(2)",
            "shifted_bound_M": str(shifted_bound),
            "row_sums": [_rational(value) for value in row_sums],
        },
        "checks": checks,
        "metrics": metrics,
        "assumptions": assumptions,
    }
    return result
