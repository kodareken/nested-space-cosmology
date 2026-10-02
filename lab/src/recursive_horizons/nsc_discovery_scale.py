"""Initial areal radius selected by the same-action lapse constraint.

The sector is constant ``Q = Q0``, ``chi = 0`` and zero canonical momenta.
It is the initial constraint, not a dynamical selection and not an evolution.
The owned operator is

    (-4 d_x^2 + Q0^2) y = G(x) / y^3,    y = sqrt(r) > 0,

    G = Q0^2 r_mag^2 + Q0 rho / (8 pi A),
    rho = M (K / Q0 + kappa Re S1) / dx,    M = 4 kappa once.

A unit coefficient on ``d_x^2`` does not reproduce ``hamilton_constraint``.
The uniform source uses the measured integral of ``rho``:

    r_flat^2 = r_mag^2 + (integral rho dx) / (8 pi A Q0 L).

Shape changes the contrast of ``G``. A covariant coordinate pullback has
areal-radius weight zero, so it does not select a different physical ``r``.
No initial-radius Newton is run and no production record is written.
"""
from __future__ import annotations

import math
import os

for _thread_var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_var, "1")

import numpy as np

from . import nsc_nested_parent_child as model
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin


SCHEMA = "NSC-DISCOVERY-SCALE-v1"
SELECTION = "initial_lapse_constraint"
DYNAMICAL_SELECTION = False
INITIAL_SOLVER_CALLED = False
PRODUCTION_RECORD_WRITTEN = False
SECOND_DERIVATIVE_COEFFICIENT = 4.0
REJECTED_SECOND_DERIVATIVE_COEFFICIENT = 1.0
HELD_OUT_CARRIER_NU = 1.15
HELD_OUT_CHIRAL_PHASE = 0.37
HELD_OUT_NF = 64
SAVED_LAYOUTS = ("original", "separated")
ONE_PERCENT_PASS_DECLARED = False

RESIDUAL_DOMAINS = {
    "source_positivity": (
        "G>0 on the period, G = Q0^2 r_mag^2 + Q0 rho / (8 pi A); "
        "required before the minimum principle is stated"
    ),
    "lapse": (
        "Hamilton constraint at constant Q0, chi=0, and zero canonical momenta; "
        "not the evolved lapse"
    ),
    "shift": (
        "momentum-constraint residual equals the retained current j = -M Pmom / dx; "
        "the mean is not deleted and beta=0 does not solve it"
    ),
}

EQUATION = (
    "(-4 d_x^2 + Q0^2) y = G/y^3, "
    "G = Q0^2 r_mag^2 + Q0 rho/(8 pi A), "
    "rho = M(K/Q0 + kappa Re S1)/dx, M=4 kappa once"
)
UNIFORM_FORMULA = "r_flat^2 = r_mag^2 + (integral rho dx) / (8 pi A Q0 L)"
CONDITIONAL_BOUND = (
    "|r(x)^2 - r_flat^2| <= ||G - mean(G)||_inf / Q0^2 "
    "when a positive solution of the initial lapse constraint exists and G>0"
)
CLOCK_FACTOR = "dτ/dt = r Q0"
PROPER_RATIO = "dτ/(Q0 dt) = r"


def owned_couplings():
    """Audited couplings that enter ``G``. Nothing here is fitted to a radius."""
    coefficients = coupling.locked_coefficients()
    area = float(coefficients["A"])
    flux = float(coefficients["flux"])
    gauge = float(coefficients["C_F"])
    weyl = float(coefficients["C_W"])
    q_value = float(coupling.CALIBRATION["b0"] / coupling.CALIBRATION["a0"])
    if not math.isfinite(q_value) or q_value <= 0.0:
        raise ValueError("owned gauge ratio Q0 is not positive")
    if area == 0.0 or not math.isfinite(area):
        raise ValueError("owned area coefficient A is not a finite nonzero constant")
    multiplicity = int(coupling.sector_multiplicity(coupling.KAPPA))
    if multiplicity != 4 * int(coupling.KAPPA):
        raise ValueError("multiplicity must remain M=4 kappa once")
    rmag2 = float(coupling.magnetic_radius_square(coefficients))
    if rmag2 != gauge * flux ** 2 / (4.0 * area):
        raise ValueError("magnetic radius drifted from C_F flux^2 / (4 A)")
    return {
        "Q": q_value,
        "A": area,
        "C_F": gauge,
        "C_W": weyl,
        "flux": flux,
        "kappa": int(coupling.KAPPA),
        "multiplicity": multiplicity,
        "rmag2": rmag2,
        "eight_pi_A": 8.0 * math.pi * area,
        "length": float(coupling.PERIOD),
        "beta0": float(coupling.CALIBRATION["beta0"]),
        "Omega": float(coupling.CALIBRATION["Omega"]),
        "carrier_k_saved": float(coupling.CARRIER_K),
        "saved_minus_phase": float(coupling.CALIBRATION["phase"]),
        "hubble_parameter_used": False,
        "magnetic_field_fitted": False,
        "C_W_enters_G_at_zero_chi": False,
        "beta_enters_rho": False,
        "prescribed_lapse_profile_enters_rho": False,
    }


def bracket_G(rho, constants):
    """``G = Q0^2 r_mag^2 + Q0 rho / (8 pi A)`` with every factor kept."""
    rho = np.asarray(rho, dtype=float)
    return (
        constants["Q"] ** 2 * constants["rmag2"]
        + constants["Q"] * rho / constants["eight_pi_A"]
    )


def source_density(K, S1, q_value, spacing, multiplicity, kappa):
    """``rho = M (K/Q0 + kappa Re S1) / dx``. ``K`` and ``S1`` are nodal."""
    if spacing == 0.0:
        raise ValueError("source density requires a nonzero sample spacing")
    return multiplicity * (np.asarray(K, dtype=float) / q_value + kappa * np.real(S1)) / spacing


def uniform_radius_square(rmag2, rho_integral, area, q_value, length):
    """Uniform-source flat radius. Couplings stay in the expression."""
    length = float(length)
    area = float(area)
    q_value = float(q_value)
    if length <= 0.0 or area == 0.0 or q_value == 0.0:
        raise ValueError("uniform radius requires nonzero period, A and Q0")
    if not all(math.isfinite(value) for value in (rmag2, rho_integral, area, q_value, length)):
        raise ValueError("uniform radius received a nonfinite input")
    return float(rmag2 + float(rho_integral) / (8.0 * math.pi * area * q_value * length))


def _nodal_mean(values):
    values = np.asarray(values, dtype=float)
    return float(np.mean(values))


def _max_abs(values):
    return float(np.max(np.abs(np.asarray(values))))


def derived_algebraic_residual(system, radius, rho, constants=None):
    """Owner lapse density at constant ``Q0`` and ``chi = p = 0``.

    ``C = (8 pi A / Q0) [D^2(r^2) - 3 (Dr)^2 - Q0^2 r^2 + G]``.
    This is the discrete ``D(F)`` form. It does not assume a product rule.
    """
    if constants is None:
        constants = owned_couplings()
    radius = np.asarray(radius, dtype=float)
    rho = np.asarray(rho, dtype=float)
    if radius.shape != rho.shape or radius.ndim != 1:
        raise ValueError("radius and rho must be matching nodal vectors")
    derivative = system.derivative
    radius_derivative = derivative @ radius
    radius_square_second = derivative @ (derivative @ (radius * radius))
    source = bracket_G(rho, constants)
    q_value = constants["Q"]
    residual = (constants["eight_pi_A"] / q_value) * (
        radius_square_second
        - 3.0 * radius_derivative ** 2
        - q_value ** 2 * radius ** 2
        + source
    )
    return residual, source


def derived_y_residual(system, y, rho, constants=None, *, coefficient=SECOND_DERIVATIVE_COEFFICIENT):
    """Strong residual from ``(-(coefficient) y'' + Q0^2 y - G/y^3)``.

    The owned lapse constraint uses ``coefficient = 4``. The value ``1`` is the
    rejected reading of the same equation.
    """
    if constants is None:
        constants = owned_couplings()
    y = np.asarray(y, dtype=float)
    rho = np.asarray(rho, dtype=float)
    if y.shape != rho.shape or np.min(y) <= 0.0:
        raise ValueError("y residual requires y > 0 on the same nodes as rho")
    source = bracket_G(rho, constants)
    second = system.derivative @ (system.derivative @ y)
    strong = -float(coefficient) * second + constants["Q"] ** 2 * y - source / y ** 3
    residual = -(constants["eight_pi_A"] / constants["Q"]) * y ** 3 * strong
    return residual, source, strong


def product_rule_gap(system, y):
    """``D^2(r^2) - 3 (Dr)^2 - 4 y^3 D^2 y`` for ``r = y^2``."""
    y = np.asarray(y, dtype=float)
    derivative = system.derivative
    radius = y * y
    radius_derivative = derivative @ radius
    radius_square_second = derivative @ (derivative @ (radius * radius))
    y_second = derivative @ (derivative @ y)
    return radius_square_second - 3.0 * radius_derivative ** 2 - 4.0 * y ** 3 * y_second


def integrated_identity(system, y, source):
    """``∫ [4 (y')^2 + Q0^2 y^2] dx = ∫ G / y^2 dx`` on a periodic grid.

    The identity is the weak form of the owned equation. It holds exactly for
    a strong solution when the periodic derivative is skew and the weights
    are uniform. It is not a dynamical energy balance.
    """
    y = np.asarray(y, dtype=float)
    source = np.asarray(source, dtype=float)
    if np.min(y) <= 0.0:
        raise ValueError("integrated identity requires y > 0")
    constants = owned_couplings()
    slope = system.derivative @ y
    left = float(system.dx * np.sum(4.0 * slope ** 2 + constants["Q"] ** 2 * y ** 2))
    right = float(system.dx * np.sum(source / y ** 2))
    return {"left": left, "right": right, "gap": abs(left - right)}


def max_principle_bounds(source, q_value):
    """Conditional bounds for a positive solution of the owned equation.

    At a maximum, ``y_max^4 <= G(x_max)/Q0^2``. At a minimum, ``G>0`` gives
    ``y_min^4 >= G_min/Q0^2``. No bound is emitted when that sign fails.
    """
    source = np.asarray(source, dtype=float)
    q_value = float(q_value)
    g_min = float(np.min(source))
    g_max = float(np.max(source))
    lower_ok = bool(g_min > 0.0 and q_value > 0.0)
    upper_ok = bool(g_max > 0.0 and q_value > 0.0)
    return {
        "G_min": g_min,
        "G_max": g_max,
        "lower_applicable": lower_ok,
        "upper_applicable": upper_ok,
        "y4_lower": (g_min / q_value ** 2) if lower_ok else None,
        "y4_upper": (g_max / q_value ** 2) if upper_ok else None,
        "domain": RESIDUAL_DOMAINS["source_positivity"],
    }


def conditional_radius_gap(source, q_value):
    """Comparison of any positive solution with the uniform ``r_flat``.

    The gap is ``||G-mean(G)||_inf / Q0^2``. It is stated only for ``G>0``.
    It is not a compulsory percent test of a numerical solve.
    """
    source = np.asarray(source, dtype=float)
    q_value = float(q_value)
    g_mean = _nodal_mean(source)
    contrast = float(np.max(np.abs(source - g_mean)))
    positive = bool(np.min(source) > 0.0 and q_value > 0.0 and g_mean > 0.0)
    return {
        "claimed": positive,
        "G_mean": g_mean,
        "G_contrast_max": contrast,
        "absolute_gap": (contrast / q_value ** 2) if positive else None,
        "relative_gap": (contrast / g_mean) if positive else None,
        "one_percent_pass_declared": ONE_PERCENT_PASS_DECLARED,
        "statement": CONDITIONAL_BOUND if positive else (
            "uniform comparison bound withheld: G>0 and Q0>0 are required"
        ),
    }


def preparation_role(layout, *, carrier_frequency=None, chiral_family=False):
    """Saved original and separated layouts are never the held-out family."""
    if layout in SAVED_LAYOUTS:
        return "previously_saved_layout"
    if layout != "override":
        raise ValueError("source layout must be original, separated, or override")
    if chiral_family and carrier_frequency is not None and float(carrier_frequency) == HELD_OUT_CARRIER_NU:
        return "held_out_initial_preparation"
    return "override_not_held_out"


def areal_pullback_weight(scale=2.0):
    """Owned affine pullback assigns weight 0 to the areal radius."""
    scale = float(scale)
    if scale <= 0.0:
        raise ValueError("positive coordinate scale required")
    jets = {
        "r": np.array([1.25]),
        "Q": np.array([0.4]),
        "L": np.array([0.4]),
    }
    mapped = model.functional_pullback(jets, scale)
    radius_ratio = float(mapped["r"][0] / jets["r"][0])
    radial_weight = math.log(radius_ratio) / math.log(scale) if radius_ratio > 0.0 else None
    clock_ratio = float((mapped["r"][0] * mapped["Q"][0]) / (jets["r"][0] * jets["Q"][0]))
    return {
        "scale": scale,
        "areal_radius_weight": radial_weight,
        "areal_radius_unchanged": bool(np.array_equal(mapped["r"], jets["r"])),
        "Q_weight": 1.0,
        "Q_scaled": bool(np.allclose(mapped["Q"], scale * jets["Q"])),
        "clock_rate_weight": math.log(clock_ratio) / math.log(scale),
        "proper_ratio_weight": radial_weight,
    }


def pull_source_invariants(measured, scale):
    """Covariant chart change ``x = scale * xi``.

    ``Q`` has weight 1, the coordinate period has weight -1 in the new chart,
    and ``∫ rho dx`` is invariant. ``r_flat`` therefore has weight 0.
    """
    scale = float(scale)
    if not math.isfinite(scale) or scale <= 0.0:
        raise ValueError("positive finite coordinate scale required")
    pulled = dict(measured)
    pulled["Q"] = float(measured["Q"]) * scale
    pulled["length"] = float(measured["length"]) / scale
    pulled["rho_integral"] = float(measured["rho_integral"])
    pulled["rho_mean"] = float(measured["rho_mean"]) * scale
    return pulled


def prediction_from_invariants(measured):
    """Flat radius, clock and proper ratio from measured source invariants.

    The probe radius used to read ``rho`` is not copied into the prediction.
    Hubble and the magnetic field are not free fit parameters.
    """
    for name in ("rmag2", "rho_integral", "A", "Q", "length", "G_min", "G_max", "G_mean"):
        if name not in measured:
            raise ValueError(f"prediction is missing measured invariant {name}")
    radius_square = uniform_radius_square(
        measured["rmag2"], measured["rho_integral"], measured["A"], measured["Q"], measured["length"],
    )
    positive_radius = bool(math.isfinite(radius_square) and radius_square > 0.0)
    radius = math.sqrt(radius_square) if positive_radius else None
    q_value = float(measured["Q"])
    source_positive = bool(measured["G_min"] > 0.0)
    contrast = float(measured["G_max"] - measured["G_mean"])
    contrast = max(contrast, float(measured["G_mean"] - measured["G_min"]))
    relative = None
    if source_positive and measured["G_mean"] > 0.0:
        relative = float(measured["G_contrast_max"]) / float(measured["G_mean"])
    exact_uniform = bool(measured["G_contrast_max"] == 0.0 and source_positive and positive_radius)
    return {
        "formula": UNIFORM_FORMULA,
        "r_flat2": radius_square,
        "radius": radius,
        "clock_rate": None if radius is None else radius * q_value,
        "proper_ratio": radius,
        "clock_definition": CLOCK_FACTOR,
        "proper_ratio_definition": PROPER_RATIO,
        "source_positivity_domain": source_positive,
        "lapse_domain": RESIDUAL_DOMAINS["lapse"],
        "shift_domain": RESIDUAL_DOMAINS["shift"],
        "exact_for_uniform_positive_G": exact_uniform,
        "conditional_absolute_gap": (
            float(measured["G_contrast_max"]) / q_value ** 2 if source_positive else None
        ),
        "conditional_relative_gap": relative,
        "conditional_bound": CONDITIONAL_BOUND if source_positive else None,
        "one_percent_pass_declared": ONE_PERCENT_PASS_DECLARED,
        "hand_set_hubble": False,
        "hand_set_magnetic_field": False,
        "solver_called": False,
        "probe_radius_used_as_prediction": False,
        "dynamical_selection": False,
        "contrast_span": contrast,
    }


def _real_envelopes(xi, spacing):
    """Owned bump lobes, then an explicit real Löwdin. Spinor columns are not repaired."""
    norms = coupling._lobe_norms()
    packets = []
    for region in range(coupling.PACKET_COUNT):
        left = region * coupling.ELL
        even = coupling._sample_lobe(xi, left, "even", spacing, norms)
        odd = coupling._sample_lobe(xi, left + coupling.ELL, "odd", spacing, norms)
        even_norm = np.linalg.norm(even)
        odd_norm = np.linalg.norm(odd)
        if min(even_norm, odd_norm) <= 0.0:
            raise ValueError("held-out lobe is unresolved")
        packets.append((even / even_norm + odd / odd_norm) / np.sqrt(2.0))
    real = np.column_stack(packets)
    envelopes, envelope_gram = coupling.lowdin(real)
    imaginary = float(np.max(np.abs(np.imag(envelopes))))
    if imaginary > 1e-8:
        raise ValueError("envelope Löwdin left the real packet; refusing a complex repair")
    envelopes = np.real(envelopes)
    gap = float(np.max(np.abs(envelopes.T @ envelopes - np.eye(coupling.PACKET_COUNT))))
    if gap > 1e-8:
        raise ValueError("real envelope Gram failed; refusing to repair the frame")
    return envelopes, {
        "envelope_lowdin_explicit": True,
        "envelope_gram_before_lowdin_max": float(np.max(np.abs(envelope_gram - np.eye(3)))),
        "envelope_orthonormal_gap": gap,
        "spinor_frame_repaired": False,
        "child_frame_repaired": False,
    }


def chiral_phase_columns(xi, spacing, *, carrier_frequency, phase_sign):
    """Chiral plus/minus pair at one carrier frequency and one phase sign.

    Plus is ``(1, i) e^{i ν x}`` times the real envelope. Minus is
    ``e^{i σ α}`` times the conjugate, with ``σ = ±1`` and ``α`` fixed.
    A column phase preserves the Gaussian ``C = Φ diag(c) Φ†`` and the Gram.
    The carrier frequency is not a column phase, so it can change ``K``.
    """
    phase_sign = float(phase_sign)
    if phase_sign not in (1.0, -1.0):
        raise ValueError("chiral phase sign must be +1 or -1")
    carrier_frequency = float(carrier_frequency)
    if not math.isfinite(carrier_frequency) or carrier_frequency == 0.0:
        raise ValueError("carrier frequency must be a finite nonzero number")
    phase = phase_sign * HELD_OUT_CHIRAL_PHASE
    if abs(phase - float(coupling.CALIBRATION["phase"])) < 1e-12:
        raise ValueError("held-out phase collides with the saved minus-column phase")
    xi = np.asarray(xi, dtype=float)
    envelopes, envelope_report = _real_envelopes(xi, spacing)
    count = xi.size
    phi0 = np.empty((count, 6), dtype=complex)
    phi1 = np.empty((count, 6), dtype=complex)
    minus_phase = np.exp(1j * phase)
    for region in range(coupling.PACKET_COUNT):
        local = xi - region * coupling.ELL
        carrier = np.exp(1j * carrier_frequency * local) * envelopes[:, region]
        plus0 = carrier / np.sqrt(2.0)
        plus1 = 1j * carrier / np.sqrt(2.0)
        phi0[:, 2 * region] = plus0
        phi1[:, 2 * region] = plus1
        phi0[:, 2 * region + 1] = minus_phase * np.conjugate(plus0)
        phi1[:, 2 * region + 1] = minus_phase * np.conjugate(plus1)
    gram = phi0.conj().T @ phi0 + phi1.conj().T @ phi1
    gram_gap = float(np.max(np.abs(gram - np.eye(6))))
    if gram_gap > 1e-9:
        raise ValueError("chiral columns are not orthonormal; refusing to repair them")
    role = preparation_role(
        "override",
        carrier_frequency=carrier_frequency,
        chiral_family=True,
    )
    preparation = dict(
        envelope_report,
        layout="override",
        role=role,
        held_out=role == "held_out_initial_preparation",
        previously_saved_layout=False,
        carrier_frequency=carrier_frequency,
        chiral_phase=phase,
        phase_sign=phase_sign,
        phase_preserves_gaussian=True,
        gram_max=gram_gap,
        mean_current_deleted=False,
        saved_carrier_reused=bool(carrier_frequency == float(coupling.CARRIER_K)),
    )
    return phi0, phi1, preparation


def build_override_pair(phase_sign, *, carrier_frequency, nf=HELD_OUT_NF):
    """Pass the new columns through the existing ``columns_override`` gate."""
    nf = int(nf)
    spacing = float(coupling.PERIOD) / nf
    xi = np.arange(nf, dtype=float) * spacing
    phi0, phi1, preparation = chiral_phase_columns(
        xi, spacing, carrier_frequency=carrier_frequency, phase_sign=phase_sign,
    )
    pair = model.build_pair(
        nf,
        source_layout="override",
        columns_override=(phi0, phi1),
        occupations=coupling.OCCUPATIONS.copy(),
    )
    if pair.source_metadata.get("layout") != "override":
        raise ValueError("constructor did not keep the override layout")
    if preparation["gram_max"] > 1e-9 or pair.source_metadata["gram_max"] > 1e-9:
        raise ValueError("constructor Gram gate failed")
    preparation = dict(
        preparation,
        constructor_layout=pair.source_metadata["layout"],
        constructor_gram_max=float(pair.source_metadata["gram_max"]),
        child_source_equals_reference=bool(pair.source_metadata["child_source_equals_reference"]),
        source_equals_reference=bool(pair.source_metadata["source_equals_reference"]),
    )
    if preparation["child_source_equals_reference"] or preparation["source_equals_reference"]:
        raise ValueError("override collapsed onto a previously saved source")
    if preparation["child_frame_repaired"]:
        raise ValueError("child frame was repaired")
    return pair, preparation


def _column_signs(system, phi0, phi1):
    signs = []
    for index in range(6):
        weights = np.zeros(6, dtype=float)
        weights[index] = 1.0
        kinetic, polar, current = coupling.column_moments(phi0, phi1, weights, system.momentum)
        signs.append({
            "index": index,
            "chirality": "plus" if index % 2 == 0 else "minus",
            "K_sum": float(np.sum(kinetic)),
            "S1_real_sum": float(np.sum(np.real(polar))),
            "S1_abs_max": float(np.max(np.abs(np.real(polar)))),
            "Pmom_sum": float(np.sum(current)),
        })
    return signs


def measure_source_invariants(pair, preparation, phi0=None, phi1=None):
    """Read ``rho`` and the current at the probe radius. Do not solve ``r``.

    The probe ``r = 1`` is only a place to evaluate the column source.
    ``force_L`` does not depend on it. The current mean is the mean that
    ``source_from_columns`` returns.
    """
    constants = owned_couplings()
    grid = pair.grid
    if phi0 is None:
        phi0 = pair.source_phi0
        phi1 = pair.source_phi1
    nodal = galerkin.blank_state(grid, phi0, phi1)
    if _max_abs(nodal.chi) != 0.0 or _max_abs(nodal.p_Q) != 0.0:
        raise ValueError("probe left the initial lapse sector")
    fine = galerkin.prolong_state(grid, nodal)
    conformal = galerkin.active_fine_system(grid, fine)
    source = coupling.source_from_columns(conformal, fine)
    if source["multiplicity_applied_once"] is not True:
        raise ValueError("source multiplicity was not applied once")
    prescribed = coupling.source_from_columns(grid.fine, fine)
    if _max_abs(prescribed["force_L"] - source["force_L"]) != 0.0:
        raise ValueError("ambient lapse or shift changed rho")
    varied = fine.copy()
    varied.r = fine.r * 1.7
    other = coupling.source_from_columns(conformal, varied)
    rho_variation = _max_abs(other["force_L"] - source["force_L"])
    spacing = float(grid.dx_q)
    kinetic, polar, pmom = coupling.column_moments(
        fine.phi0, fine.phi1, conformal.occupations, conformal.momentum,
    )
    rho = source["force_L"] / spacing
    rebuilt = source_density(
        kinetic, polar, constants["Q"], spacing, constants["multiplicity"], constants["kappa"],
    )
    if _max_abs(rebuilt - rho) > 1e-12:
        raise ValueError("rho sign does not match M(K/Q0 + kappa Re S1)/dx")
    current = source["force_beta"] / spacing
    current_from_momentum = -constants["multiplicity"] * pmom / spacing
    if _max_abs(current_from_momentum - current) > 1e-12:
        raise ValueError("shift source sign does not match j = -M Pmom / dx")
    # Equal CAR weights on a chiral pair can cancel. That cancellation is kept.
    if preparation.get("mean_current_deleted"):
        raise ValueError("mean current was deleted")
    source_bracket = bracket_G(rho, constants)
    occupations = np.asarray(conformal.occupations, dtype=float)
    if occupations.shape != (6,) or np.min(occupations) <= 0.0 or np.max(occupations) > 1.0:
        raise ValueError("CAR occupation weights left (0, 1]")
    if not np.array_equal(occupations, coupling.OCCUPATIONS):
        raise ValueError("measurement changed the owned CAR occupations")
    length = float(grid.length)
    rho_integral = float(spacing * np.sum(rho))
    return {
        "Q": constants["Q"],
        "A": constants["A"],
        "C_F": constants["C_F"],
        "C_W": constants["C_W"],
        "flux": constants["flux"],
        "kappa": constants["kappa"],
        "multiplicity": constants["multiplicity"],
        "rmag2": constants["rmag2"],
        "length": length,
        "dx": spacing,
        "rho": rho,
        "rho_integral": rho_integral,
        "rho_mean": _nodal_mean(rho),
        "rho_min": float(np.min(rho)),
        "rho_max": float(np.max(rho)),
        "rho_contrast_max": float(np.max(np.abs(rho - _nodal_mean(rho)))),
        "G": source_bracket,
        "G_min": float(np.min(source_bracket)),
        "G_max": float(np.max(source_bracket)),
        "G_mean": _nodal_mean(source_bracket),
        "G_contrast_max": float(np.max(np.abs(source_bracket - _nodal_mean(source_bracket)))),
        "current": current,
        "current_mean": _nodal_mean(current),
        "current_max": _max_abs(current),
        "current_integral": float(spacing * np.sum(current)),
        "shift_residual_mean_subtracted": False,
        "K_sum": float(np.sum(kinetic)),
        "S1_real_max": float(np.max(np.abs(np.real(polar)))),
        "Pmom_sum": float(np.sum(pmom)),
        "column_signs": _column_signs(conformal, fine.phi0, fine.phi1),
        "gram_max": float(preparation["gram_max"]),
        "occupations": occupations.tolist(),
        "occupations_car_admissible": True,
        "rho_independent_of_r_max": rho_variation,
        "multiplicity_applied_once": True,
        "ambient_shift_changes_force_L": False,
        "probe_radius": 1.0,
        "probe_radius_is_prediction": False,
        "layout": preparation["layout"],
        "held_out": bool(preparation["held_out"]),
        "role": preparation["role"],
        "carrier_frequency": float(preparation["carrier_frequency"]),
        "phase_sign": float(preparation["phase_sign"]),
        "child_frame_repaired": False,
        "child_source_equals_reference": bool(preparation["child_source_equals_reference"]),
        "lapse_sector": True,
        "domains": dict(RESIDUAL_DOMAINS),
    }


def constant_radius_owner_comparison(system, rho, radius):
    """Owned residual of a constant radius against ``rho - mean(rho)``."""
    count = int(system.points)
    residual = coupling._radius_residual(system, np.full(count, float(radius)), np.asarray(rho, dtype=float))
    contrast = np.asarray(rho, dtype=float) - _nodal_mean(rho)
    gap = _max_abs(residual - contrast)
    scale = max(1.0, _max_abs(residual), _max_abs(contrast))
    return {
        "owner_residual_max": _max_abs(residual),
        "contrast_max": _max_abs(contrast),
        "identity_gap": gap,
        "relative_identity_gap": gap / scale,
    }


def _member_record(pair, preparation):
    measured = measure_source_invariants(pair, preparation)
    predicted = prediction_from_invariants(measured)
    comparison = None
    if predicted["radius"] is not None:
        comparison = constant_radius_owner_comparison(pair.grid.fine, measured["rho"], predicted["radius"])
    return {
        "role": preparation["role"],
        "held_out": bool(preparation["held_out"]),
        "carrier_frequency": float(preparation["carrier_frequency"]),
        "phase_sign": float(preparation["phase_sign"]),
        "chiral_phase": float(preparation["chiral_phase"]),
        "gram_max": float(preparation["constructor_gram_max"]),
        "column_gram_max": float(preparation["gram_max"]),
        "child_frame_repaired": False,
        "child_source_equals_reference": bool(preparation["child_source_equals_reference"]),
        "source_equals_reference": bool(preparation["source_equals_reference"]),
        "mean_current_deleted": False,
        "shift_residual_mean": measured["current_mean"],
        "shift_residual_mean_subtracted": False,
        "K_sum": measured["K_sum"],
        "S1_real_max": measured["S1_real_max"],
        "Pmom_sum": measured["Pmom_sum"],
        "column_signs": measured["column_signs"],
        "occupations": measured["occupations"],
        "rho_integral": measured["rho_integral"],
        "rho_mean": measured["rho_mean"],
        "rho_contrast_max": measured["rho_contrast_max"],
        "rho_independent_of_r_max": measured["rho_independent_of_r_max"],
        "G_min": measured["G_min"],
        "G_max": measured["G_max"],
        "G_mean": measured["G_mean"],
        "G_contrast_max": measured["G_contrast_max"],
        "source_positivity": bool(measured["G_min"] > 0.0),
        "prediction": predicted,
        "owner_comparison": comparison,
        "measured": measured,
    }


def future_heldout_experiment(nf=HELD_OUT_NF):
    """Measure the withheld chiral family and predict from those invariants.

    Both signs go through ``columns_override``. A same-constructor carrier at
    the saved wavenumber is a frequency control, not a held-out layout and
    not the original or separated source. No radius solve is started.
    """
    if float(HELD_OUT_CARRIER_NU) == float(coupling.CARRIER_K):
        raise ValueError("held-out carrier must differ from the saved carrier")
    members = []
    first_pair = None
    for sign in (1.0, -1.0):
        pair, preparation = build_override_pair(
            sign, carrier_frequency=HELD_OUT_CARRIER_NU, nf=nf,
        )
        if preparation["role"] != "held_out_initial_preparation":
            raise ValueError("chiral family was not marked held out")
        if first_pair is None:
            first_pair = pair
        members.append(_member_record(pair, preparation))
    spacing = float(coupling.PERIOD) / int(nf)
    xi = np.arange(int(nf), dtype=float) * spacing
    _phi0, _phi1, control_preparation = chiral_phase_columns(
        xi, spacing, carrier_frequency=float(coupling.CARRIER_K), phase_sign=1.0,
    )
    if control_preparation["held_out"] or control_preparation["role"] == "held_out_initial_preparation":
        raise ValueError("saved carrier frequency was labeled held out")
    control_preparation = dict(
        control_preparation,
        constructor_layout="override",
        constructor_gram_max=control_preparation["gram_max"],
        child_source_equals_reference=False,
        source_equals_reference=False,
    )
    # The control is measured on the held-out quadrature grid. It is not a
    # second saved layout, and it does not repair that grid's source columns.
    control_measured = measure_source_invariants(first_pair, control_preparation, _phi0, _phi1)
    control_predicted = prediction_from_invariants(control_measured)
    plus, minus = members
    phase_gap = abs(plus["rho_integral"] - minus["rho_integral"])
    kinetic_gap = abs(plus["K_sum"] - control_measured["K_sum"])
    return {
        "schema": SCHEMA,
        "selection": SELECTION,
        "dynamical_selection": DYNAMICAL_SELECTION,
        "initial_solver_called": INITIAL_SOLVER_CALLED,
        "production_record_written": PRODUCTION_RECORD_WRITTEN,
        "validation": "later initial source-consistent solver; not run",
        "held_out_carrier_nu": HELD_OUT_CARRIER_NU,
        "held_out_chiral_phase": HELD_OUT_CHIRAL_PHASE,
        "saved_layouts_used_as_held_out": False,
        "members": members,
        "phase_rho_integral_gap": phase_gap,
        "frequency_control": {
            "role": control_preparation["role"],
            "held_out": False,
            "carrier_frequency": float(coupling.CARRIER_K),
            "K_sum": control_measured["K_sum"],
            "rho_integral": control_measured["rho_integral"],
            "G_contrast_max": control_measured["G_contrast_max"],
            "prediction": control_predicted,
            "measured_on_heldout_grid_without_replacing_its_source": True,
        },
        "frequency_changes_kinetic_invariant": bool(kinetic_gap > 0.0),
        "kinetic_sum_gap": kinetic_gap,
        "domains": dict(RESIDUAL_DOMAINS),
    }


def _shape_controls(system, constants, rho_mean):
    """Same mean source, two shapes. ``r_flat`` stays; the contrast does not."""
    uniform = np.full(system.points, rho_mean)
    shaped = rho_mean + 0.05 * np.cos(4.0 * math.pi * system.xi / system.length)
    records = []
    for name, rho in (("uniform", uniform), ("shaped", shaped)):
        integral = float(system.dx * np.sum(rho))
        radius_square = uniform_radius_square(
            constants["rmag2"], integral, constants["A"], constants["Q"], system.length,
        )
        source = bracket_G(rho, constants)
        gap = conditional_radius_gap(source, constants["Q"])
        records.append({
            "name": name,
            "r_flat2": radius_square,
            "rho_contrast_max": float(np.max(np.abs(rho - _nodal_mean(rho)))),
            "conditional_absolute_gap": gap["absolute_gap"],
            "bound_claimed": gap["claimed"],
        })
    return {
        "same_flat_radius": bool(records[0]["r_flat2"] == records[1]["r_flat2"] or abs(
            records[0]["r_flat2"] - records[1]["r_flat2"]) <= 1e-12),
        "contrast_increased": bool(records[1]["rho_contrast_max"] > records[0]["rho_contrast_max"]),
        "cases": records,
    }


def analytic_controls(points=32, rho_value=0.15):
    """Uniform identity, contrast identity, integral identity, and maximum principle."""
    system, _state = coupling.build_system(int(points))
    constants = owned_couplings()
    q_value = constants["Q"]
    if abs(q_value - system.calibration["b0"] / system.calibration["a0"]) > 0.0:
        raise ValueError("system Q0 drifted from the owned gauge ratio")
    uniform_rho = np.full(system.points, float(rho_value))
    uniform_integral = float(system.dx * np.sum(uniform_rho))
    uniform_square = uniform_radius_square(
        constants["rmag2"], uniform_integral, constants["A"], q_value, system.length,
    )
    if uniform_square <= 0.0:
        raise ValueError("manufactured uniform source left the positive chart")
    uniform_radius = math.sqrt(uniform_square)
    uniform_residual = coupling._radius_residual(
        system, np.full(system.points, uniform_radius), uniform_rho,
    )
    variable = uniform_rho + 0.04 * np.cos(2.0 * math.pi * system.xi / system.length)
    variable_integral = float(system.dx * np.sum(variable))
    variable_square = uniform_radius_square(
        constants["rmag2"], variable_integral, constants["A"], q_value, system.length,
    )
    variable_comparison = constant_radius_owner_comparison(
        system, variable, math.sqrt(variable_square),
    )
    y = 1.5 + 0.02 * np.cos(2.0 * math.pi * system.xi / system.length)
    second = system.derivative @ (system.derivative @ y)
    manufactured = y ** 3 * (-SECOND_DERIVATIVE_COEFFICIENT * second + q_value ** 2 * y)
    bounds = max_principle_bounds(manufactured, q_value)
    identity = integrated_identity(system, y, manufactured)
    y_max = float(np.max(y))
    y_min = float(np.min(y))
    principle_holds = bool(
        bounds["upper_applicable"]
        and bounds["lower_applicable"]
        and y_max ** 4 <= bounds["y4_upper"] * (1.0 + 1e-12)
        and y_min ** 4 >= bounds["y4_lower"] * (1.0 - 1e-12)
    )
    return {
        "uniform_residual_max": _max_abs(uniform_residual),
        "uniform_r_flat2": uniform_square,
        "uniform_radius": uniform_radius,
        "uniform_clock_rate": uniform_radius * q_value,
        "uniform_proper_ratio": uniform_radius,
        "uniform_formula_gap": abs(
            uniform_square
            - (constants["rmag2"] + float(rho_value) / (8.0 * math.pi * constants["A"] * q_value))
        ),
        "variable_owner_comparison": variable_comparison,
        "variable_r_flat2": variable_square,
        "manufactured_G_min": float(np.min(manufactured)),
        "integrated_identity_gap": identity["gap"],
        "max_principle_holds": principle_holds,
        "max_principle": bounds,
        "shape": _shape_controls(system, constants, float(rho_value)),
        "one_percent_pass_declared": ONE_PERCENT_PASS_DECLARED,
    }


def coordinate_scale_control(scale=2.0):
    """Covariant scaling leaves ``r_flat`` fixed. A pure period stretch does not."""
    constants = owned_couplings()
    measured = {
        "rmag2": constants["rmag2"],
        "rho_integral": 0.2 * constants["length"],
        "rho_mean": 0.2,
        "A": constants["A"],
        "Q": constants["Q"],
        "length": constants["length"],
        "G_min": 1.0,
        "G_max": 1.0,
        "G_mean": 1.0,
        "G_contrast_max": 0.0,
    }
    base = prediction_from_invariants(measured)["r_flat2"]
    pulled = pull_source_invariants(measured, scale)
    pulled_square = uniform_radius_square(
        pulled["rmag2"], pulled["rho_integral"], pulled["A"], pulled["Q"], pulled["length"],
    )
    wrong_period = uniform_radius_square(
        measured["rmag2"], measured["rho_integral"], measured["A"], measured["Q"],
        measured["length"] * float(scale),
    )
    weights = areal_pullback_weight(scale)
    ratio = pulled_square / base
    weight = 0.0 if abs(ratio - 1.0) <= 1e-15 else math.log(ratio) / (2.0 * math.log(float(scale)))
    return {
        "areal_radius_weight": weight,
        "pulled_r_flat2": pulled_square,
        "base_r_flat2": base,
        "period_only_r_flat2": wrong_period,
        "period_only_changes_radius": bool(abs(wrong_period - base) > 1e-12),
        "pullback": weights,
    }


def owned_projection_agreement(nf=HELD_OUT_NF):
    """Projected owner residual versus the derived algebraic and ``y`` forms."""
    grid = galerkin.build_grid(int(nf), gauge="conformal")
    constants = owned_couplings()
    radius = 1.15 + 0.04 * np.cos(2.0 * math.pi * grid.xi_g / grid.length)
    rho = 0.12 + 0.03 * np.cos(4.0 * math.pi * grid.xi_q / grid.length)
    projected, full = galerkin.projected_radius_operator(grid, radius, rho)
    fine_radius = galerkin.prolong_geometry(grid, radius)
    if np.min(fine_radius) <= 0.0:
        raise ValueError("projection sample left the positive chart")
    algebraic, _source = derived_algebraic_residual(grid.fine, fine_radius, rho, constants)
    y = np.sqrt(fine_radius)
    y_residual, _source_y, _strong = derived_y_residual(grid.fine, y, rho, constants)
    reviewer, _source_r, _strong_r = derived_y_residual(
        grid.fine, y, rho, constants, coefficient=REJECTED_SECOND_DERIVATIVE_COEFFICIENT,
    )
    scaled_gap = (constants["eight_pi_A"] / constants["Q"]) * product_rule_gap(grid.fine, y)
    return {
        "algebraic_owner_gap": _max_abs(full - algebraic),
        "projected_pull_gap": _max_abs(projected - galerkin.pull_geometry(grid, algebraic)),
        "defect_explains_y_gap": _max_abs((full - y_residual) - scaled_gap),
        "y_form_gap": _max_abs(full - y_residual),
        "reviewer_unit_coefficient_gap": _max_abs(full - reviewer),
        "owned_coefficient": SECOND_DERIVATIVE_COEFFICIENT,
        "rejected_coefficient": REJECTED_SECOND_DERIVATIVE_COEFFICIENT,
    }


def _plain(value):
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items() if key not in {"rho", "G", "current", "measured"}}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    return value


def derivation_report(nf=HELD_OUT_NF):
    """Analytic controls plus the withheld preparation. No solve and no file."""
    constants = owned_couplings()
    analytic = analytic_controls()
    projection = owned_projection_agreement(nf)
    heldout = future_heldout_experiment(nf)
    coordinate = coordinate_scale_control()
    return {
        "schema": SCHEMA,
        "equation": EQUATION,
        "uniform_formula": UNIFORM_FORMULA,
        "conditional_bound": CONDITIONAL_BOUND,
        "selection": SELECTION,
        "dynamical_selection": DYNAMICAL_SELECTION,
        "initial_solver_called": INITIAL_SOLVER_CALLED,
        "production_record_written": PRODUCTION_RECORD_WRITTEN,
        "one_percent_pass_declared": ONE_PERCENT_PASS_DECLARED,
        "couplings": constants,
        "domains": dict(RESIDUAL_DOMAINS),
        "analytic": analytic,
        "projection": projection,
        "coordinate": coordinate,
        "heldout": heldout,
    }


def public_summary(report=None):
    """JSON-ready scalars. Nodal source arrays stay out of the summary."""
    if report is None:
        report = derivation_report()
    return _plain(report)
