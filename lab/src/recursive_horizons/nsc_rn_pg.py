"""Coupled flat-slice Painlevé–Gullstrand state for classical RN.

The geometry is solved from the common Einstein–Maxwell–Dirac action before
the areal gauge is imposed, then reduced to q = 1 with areal radius r.
Source-free Reissner–Nordström is the stationary point

    β = N √(2 m_Q / r − Q² / r²),    N(r_out) = 1.

The lapse is not held at 1 once the shell carries momentum. The inner mass
is a stock advanced by the q-equation flux, not a second boundary value.
No evidence writer and no production trajectory live here.
"""
from __future__ import annotations

import importlib
import math
import os

for _thread_name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_thread_name, "1")

import numpy as np
from scipy import sparse

from . import nsc_rn_observables as observables
from . import nsc_rn_reference as reference
from . import nsc_rn_source as source

API_VERSION = "nsc-rn-pg-v2"
# Integer monopole q=1. A and C_F*flux² match the audited source action.
# Flux is not solved from r_m. G_N = 1/(16 π A) is not set to 1.
ACTION_A = float(reference.AUDITED_A)
ACTION_FLUX = 1.0
_AUDITED_P2 = float(reference.magnetic_radius_square(
    reference.AUDITED_A, reference.AUDITED_C_F, reference.AUDITED_FLUX,
))
ACTION_C_F = _AUDITED_P2 * 4.0 * ACTION_A / (ACTION_FLUX * ACTION_FLUX)
A4 = ACTION_A
NEWTON_G = float(reference.newton_constant(ACTION_A))
SHELL_KAPPA = reference.KAPPA
JET_PROVENANCE = "one_sided_eighth_and_centered_fourth_nodal_derivative"
JET_SCOPE = "spatial_radial_polynomial_one_sided"
TIME_JET_PROVENANCE = "characteristic_rays_and_constraint_linearization"

# Centered fourth-order nodal jets. One-sided jets are eighth-order Fornberg
# weights: a sixth-order second jet leaves an O(10^-2) excision error on the
# 257-point chart because β branches at P²/(2M), about half a magnetic radius
# inside the inner boundary.
_CENTER_D1 = np.array((1.0 / 12.0, -2.0 / 3.0, 0.0, 2.0 / 3.0, -1.0 / 12.0))
_CENTER_D2 = np.array((-1.0 / 12.0, 4.0 / 3.0, -2.5, 4.0 / 3.0, -1.0 / 12.0))
_ONE_SIDED_D1_WIDTH = 9
_ONE_SIDED_D2_WIDTH = 10
_SIGMA1 = np.array(((0.0, 1.0), (1.0, 0.0)), dtype=complex)
_PLUS = np.array((1.0, 1.0j)) / math.sqrt(2.0)
_MINUS = np.array((1.0, -1.0j)) / math.sqrt(2.0)


class ChartAdmissionError(ValueError):
    """The flat-slice chart or its CFL restriction does not admit the stage."""

    def __init__(self, blocker):
        super().__init__(blocker)
        self.blocker = str(blocker)


def formulas():
    """Continuum relations implemented by this owner. V4 = 0, no Weyl pole."""
    return {
        "api_version": API_VERSION,
        "pre_gauge_chart": "ds^2 = N^2 dt^2 - q^2 (dx + beta dt)^2 - r^2 dOmega^2",
        "areal_chart": observables.FORMULAS["chart"],
        "action": "S = ∫ sqrt|g| [-A R_L - C_F F^2] + S_Dirac, V4 = 0",
        "newton_constant": "G_N = 1/(16 pi A), shared with the audited action, not forced to 1",
        "integer_flux": 1.0,
        "flux_fitted": False,
        "flux_convention": (
            "integer q=1 with C_F rescaled by (audited_flux/1)^2 so classical P^2 is unchanged; "
            "not the audited flux=4 quantized charge sector"
        ),
        "quantized_charge_sector_matches_audited_flux": False,
        "classical_P2_preserved": True,
        "metric_charge": "r_m^2 = C_F flux^2 / (4 A), derived from the action",
        "charge_corrected_mass": "m_Q = r v^2/2 + P^2/(2r), v = beta/N",
        "misner_sharp": observables.FORMULAS["misner_sharp"],
        "hamiltonian_constraint": "m_Q' = 4 pi G r^2 (rho - v j) = G (F_N + v F_beta)",
        "lapse_condition": "log(N)' = 4 pi G r j / v = -G F_beta / (r v), N(r_out) = 1",
        "inner_mass_rate": (
            "m_Q_t = G N [2 v F_N + (1+v^2) F_beta + v (r/N) F_r], second v F_N required"
        ),
        "intrinsic_stress_cache": (
            "F_N, F_beta and F_r/N are cached at N=1, beta=0; "
            "the Dirac operator uses the reconstructed N and beta"
        ),
        "lapse_time_constraint": "N(r_out)=1 implies N_t(r_out)=0; N is not a free zero",
        "characteristics": "lambda_pm = -beta ± N",
        "dirac_principal": "partial_t phi_c + lambda_c partial_r phi_c = angular and SBP corrections",
        "angular": (
            "characteristic frame: partial_t chi_+ gets -alpha chi_-, "
            "partial_t chi_- gets +alpha chi_+, alpha = kappa N/r"
        ),
        "shell": "kappa = 1, multiplicity 4 counted once",
        "electric_current": 0,
        "magnetic_flux_fixed": True,
        "outer_boundary": "SAT on the incoming characteristic only",
        "excision": "pure outflow, no field set to zero",
        "cfl": observables.FORMULAS["full_cfl"],
        "weyl_pole": False,
        "default_outer": "r_out = 32 r_m, r_m = sqrt(P^2), not r_minus",
        "forces": "nsc_rn_source.source_evaluate, multiplicity 4 once",
        "q_adm": 1.0,
        "physical_rho": "F_N / (4 pi q r^2), distinct from the angular-integrated rho",
        "spatial_jets": JET_SCOPE,
        "time_jets": TIME_JET_PROVENANCE,
        "dynamic_R4_from_spatial_jets": False,
        "probability_stocks_separate_from_mass_flux_and_sat": True,
    }


def _load(module_name):
    try:
        return importlib.import_module(module_name)
    except ImportError:
        return None


def sbp_operator(points, spacing):
    """Source-owned diagonal-norm SBP on a uniform nonperiodic interval."""
    step = float(spacing)
    count = int(points)
    if count < 8 or not math.isfinite(step) or step <= 0.0:
        raise ChartAdmissionError("SBP grid needs at least eight positive spacings")
    pack = source.radial_sbp(count, 0.0, step * (count - 1))
    if pack["closure_residual"] > 1e-10:
        raise ChartAdmissionError("SBP quadrature failed its boundary closure")
    derivative = np.asarray(pack["derivative"], dtype=float)
    return np.asarray(pack["weights"], dtype=float), sparse.csr_matrix(derivative), derivative


def _fornberg(width, derivative, at):
    """Polynomial weights for d^m/dx^m at node `at` on nodes 0..width-1."""
    nodes = np.arange(width, dtype=float)
    order = int(derivative)
    coeff = np.zeros((order + 1, width, width))
    coeff[0, 0, 0] = 1.0
    c1 = 1.0
    target = float(at)
    for column in range(1, width):
        c2 = 1.0
        for row in range(column):
            gap = nodes[column] - nodes[row]
            c2 *= gap
            for level in range(min(order, column) + 1):
                previous = coeff[level, column - 1, row]
                lower = coeff[level - 1, column - 1, row] if level else 0.0
                coeff[level, column, row] = ((nodes[column] - target) * previous - level * lower) / gap
        for level in range(min(order, column) + 1):
            previous = coeff[level, column - 1, column - 1]
            lower = coeff[level - 1, column - 1, column - 1] if level else 0.0
            coeff[level, column, column] = (c1 / c2) * (
                level * lower - (nodes[column - 1] - target) * previous
            )
        c1 = c2
    weights = coeff[order, width - 1].copy()
    for power in range(width):
        moment = float(np.dot(weights, (nodes - target) ** power))
        expected = float(math.factorial(order)) if power == order else 0.0
        if abs(moment - expected) > 1e-8 * (1.0 + abs(expected)):
            raise RuntimeError("nodal stencil failed its polynomial moment")
    return weights


def _boundary_table(max_width, derivative):
    table = {}
    for width in range(5, max_width + 1):
        table[width] = (
            _fornberg(width, derivative, 0),
            _fornberg(width, derivative, 1),
        )
    return table


_LEFT_D1 = _boundary_table(_ONE_SIDED_D1_WIDTH, 1)
_LEFT_D2 = _boundary_table(_ONE_SIDED_D2_WIDTH, 2)


def _differentiate(values, spacing):
    """First and second nodal derivatives. Ends are one-sided, not D@D.

    Interior nodes use the centered fourth-order stencil. The two nodes at
    each end use an eighth-order one-sided stencil on the longest window
    that fits. A shorter grid keeps the same centered rule and narrows only
    the one-sided window.
    """
    sample = np.asarray(values, dtype=float)
    step = float(spacing)
    count = sample.size
    if count < 5 or not math.isfinite(step) or step <= 0.0:
        raise ChartAdmissionError("nodal jets need at least five positive spacings")
    first = np.empty(count)
    second = np.empty(count)
    width1 = min(_ONE_SIDED_D1_WIDTH, count)
    width2 = min(_ONE_SIDED_D2_WIDTH, count)
    left1 = _LEFT_D1[width1]
    left2 = _LEFT_D2[width2]
    last = count - 1
    for index in range(count):
        if 2 <= index <= count - 3:
            window = sample[index - 2:index + 3]
            first[index] = np.dot(_CENTER_D1, window) / step
            second[index] = np.dot(_CENTER_D2, window) / step**2
        elif index < 2:
            first[index] = np.dot(left1[index], sample[:width1]) / step
            second[index] = np.dot(left2[index], sample[:width2]) / step**2
        else:
            mirror = last - index
            first[index] = -np.dot(left1[mirror], sample[-width1:][::-1]) / step
            second[index] = np.dot(left2[mirror], sample[-width2:][::-1]) / step**2
    return first, second


def action_model(mass, r_m=None):
    """Shared Einstein–Maxwell data. Flux stays the integer 1.

    r_m = sqrt(C_F flux²/(4A)) comes out of the action. A requested r_m that
    disagrees is refused. Integer q=1 keeps classical P² by rescaling the
    audited Maxwell coefficient with (audited flux/1)². That is not the
    flux=4 quantized charge sector.
    """
    if ACTION_FLUX != 1.0 or abs(ACTION_FLUX - round(ACTION_FLUX)) > 0.0:
        raise ChartAdmissionError("magnetic flux is the integer monopole 1")
    try:
        model = source.action_reference(
            mass, A=ACTION_A, C_F=ACTION_C_F, flux=ACTION_FLUX,
        )
    except ValueError as error:
        raise ChartAdmissionError(str(error)) from error
    audited = source.action_reference(mass)
    if abs(model.A - audited.A) > 0.0 or abs(model.G_N - audited.G_N) > 1e-12 * audited.G_N:
        raise ChartAdmissionError("Newton constant is not the shared action")
    if abs(model.magnetic_r2 - audited.magnetic_r2) > 1e-10 * max(1.0, audited.magnetic_r2):
        raise ChartAdmissionError("P^2 drifted from the source action")
    if abs(model.G_N - 1.0) < 1e-8:
        raise ChartAdmissionError("G_N was forced to 1")
    if model.flux != ACTION_FLUX:
        raise ChartAdmissionError("flux was replaced after the action was fixed")
    derived = math.sqrt(model.magnetic_r2)
    if r_m is not None:
        requested = float(r_m)
        if not math.isfinite(requested) or requested <= 0.0:
            raise ChartAdmissionError("r_m must be a positive magnetic radius")
        if abs(requested - derived) > 1e-8 * max(1.0, derived):
            raise ChartAdmissionError(
                "r_m is derived from C_F flux^2/(4A); refusing to fit the integer flux"
            )
    if model.horizons["r_minus"] is None or model.extremal:
        raise ChartAdmissionError("flat-slice RN needs a non-extremal horizon pair")
    return model


def build_grid(mass, *, r_m=1.0, points=33, r_out=None):
    """Excision at the horizon midpoint and outer radius 32 r_m by default.

    Child and parent are disjoint areal intervals on this one chart.
    """
    model = action_model(mass, r_m)
    derived_rm = math.sqrt(model.magnetic_r2)
    horizons = model.horizons
    r_minus = float(horizons["r_minus"])
    r_plus = float(horizons["r_plus"])
    gap = r_plus - r_minus
    outer = 32.0 * float(r_m) if r_out is None else float(r_out)
    excision = 0.5 * (r_minus + r_plus)
    child = (r_plus + gap / 16.0, r_plus + gap / 2.0)
    parent = (child[1], outer)
    if not (r_minus < excision < r_plus < child[0] < child[1] < outer):
        raise ChartAdmissionError("child and parent cuts are not ordered")
    count = int(points)
    radius = np.linspace(excision, outer, count)
    spacing = float(radius[1] - radius[0])
    weights, q_matrix, derivative = sbp_operator(count, spacing)
    return {
        "radius": radius,
        "weights": weights,
        "q_matrix": q_matrix,
        "derivative": derivative,
        "spacing": spacing,
        "points": count,
        "mass": float(model.mass),
        "charge": math.sqrt(model.magnetic_r2),
        "magnetic_r2": float(model.magnetic_r2),
        "r_m": derived_rm,
        "A": float(model.A),
        "C_F": float(model.C_F),
        "flux": float(model.flux),
        "G_N": float(model.G_N),
        "V4": 0.0,
        "flux_fitted": False,
        "integer_flux": ACTION_FLUX,
        "r_out": outer,
        "r_minus": r_minus,
        "r_plus": r_plus,
        "model": model,
        "cuts": {
            "excision": excision,
            "child": child,
            "parent": parent,
            "outer": outer,
            "same_areal_radius": False,
        },
        "child_interval": child,
        "parent_interval": parent,
        "same_areal_radius": False,
        "periodic": False,
        "kappa": SHELL_KAPPA,
        "electric_current": 0.0,
        "magnetic_flux_fixed": True,
        "weyl_pole": False,
        "areal_gauge_q": 1.0,
    }


def _shell_weight():
    """The source applies this factor once. Callers do not multiply again."""
    return float(source.shell_multiplicity(SHELL_KAPPA))


def _to_sigma(phi):
    return _PLUS[:, None, None] * phi[0] + _MINUS[:, None, None] * phi[1]


def _nodal_derivative(grid, phi):
    derivative = np.empty_like(phi)
    operator = grid["derivative"]
    for component in range(2):
        for column in range(phi.shape[2]):
            derivative[component, :, column] = operator @ phi[component, :, column]
    return derivative


def matter_densities(phi, grid, occupations, *, kappa=SHELL_KAPPA, occupation_matrix=None):
    """Intrinsic stress at N=1, β=0. The Dirac operator is not evaluated here.

    At q_ADM=1, F_N, F_β and F_r/N do not depend on N or β. Physical ρ is
    F_N/(4π r²). The angular-integrated ρ is larger by 4π and is not reused
    as ρ. Multiplicity 4 is already inside the forces.
    """
    count = grid["points"]
    ones = np.ones(count)
    zeros = np.zeros(count)
    evaluated = source.source_evaluate(
        phi,
        occupations,
        ones,
        zeros,
        ones,
        grid["radius"],
        grid["weights"],
        grid["derivative"],
        kappa=kappa,
        occupation_matrix=occupation_matrix,
    )
    if not evaluated["multiplicity_applied_once"] or int(evaluated["multiplicity"]) != 4:
        raise ChartAdmissionError("source shell factor was not applied once")
    q_adm = np.asarray(evaluated["q_adm"], dtype=float)
    if np.max(np.abs(q_adm - 1.0)) > 1e-12:
        raise ChartAdmissionError("canonical chart requires q_ADM = 1")
    force_n = np.asarray(evaluated["forces"]["N"], dtype=float)
    force_beta = np.asarray(evaluated["forces"]["beta"], dtype=float)
    force_r_unit = np.asarray(evaluated["forces"]["r"], dtype=float)
    physical_rho = np.asarray(evaluated["physical"]["rho"], dtype=float)
    integrated_rho = np.asarray(evaluated["angular_integrated"]["rho"], dtype=float)
    volume = grid["radius"] ** 2
    rho_gap = np.max(np.abs(physical_rho * float(source.FOUR_PI) * volume - force_n))
    if rho_gap > 1e-8 * (1.0 + np.max(np.abs(force_n))):
        raise ChartAdmissionError("physical rho is not F_N / (4 pi r^2)")
    return {
        "eps_N": force_n,
        "eps_beta": force_beta,
        "force_r_unit": force_r_unit,
        "angular": -force_r_unit * grid["radius"],
        "electric_current": np.asarray(evaluated["electric_current"], dtype=float),
        "source": "nsc_rn_source.source_evaluate",
        "multiplicity_included": True,
        "multiplicity": int(evaluated["multiplicity"]),
        "q_adm": 1.0,
        "intrinsic_stress_cached_at_unit_lapse": True,
        "stress_iteration": False,
        "physical_rho": physical_rho,
        "angular_integrated_rho": integrated_rho,
        "physical_minus_integrated_factor": float(source.FOUR_PI),
    }


def _positive_shift_factor(mass_value, radius, charge2):
    argument = 2.0 * mass_value / radius - charge2 / radius**2
    if argument <= 0.0 or not math.isfinite(argument):
        raise ChartAdmissionError("pg shift is not real on the flat slice")
    return math.sqrt(argument), argument


def _mass_slope(mass_value, index, grid, eps_n, eps_beta, weight):
    radius = float(grid["radius"][index])
    shift_factor, _argument = _positive_shift_factor(mass_value, radius, grid["charge"] ** 2)
    charge2 = grid["charge"] ** 2
    log_slope = -NEWTON_G * weight * float(eps_beta[index]) / (radius * shift_factor)
    enclosed = 2.0 * mass_value - charge2 / radius
    slope = NEWTON_G * weight * float(eps_n[index]) - log_slope * enclosed
    return slope, log_slope, shift_factor


def _mass_slope_mu_derivative(mass_value, index, grid, eps_beta, weight):
    radius = float(grid["radius"][index])
    shift_factor, argument = _positive_shift_factor(mass_value, radius, grid["charge"] ** 2)
    charge2 = grid["charge"] ** 2
    coefficient = NEWTON_G * weight * float(eps_beta[index])
    dlog_dmu = coefficient / (radius**2 * shift_factor**3)
    log_slope = -coefficient / (radius * shift_factor)
    enclosed = 2.0 * mass_value - charge2 / radius
    return -dlog_dmu * enclosed - 2.0 * log_slope, argument


def _stress_channels(densities, radius):
    """F_N, F_β and unit-lapse F_r. Injected angular data use F_r = -angular/r."""
    force_n = np.asarray(densities["eps_N"], dtype=float)
    force_beta = np.asarray(densities["eps_beta"], dtype=float)
    if "force_r_unit" in densities:
        force_r_unit = np.asarray(densities["force_r_unit"], dtype=float)
    else:
        angular = np.asarray(densities["angular"], dtype=float)
        force_r_unit = -angular / radius
    return force_n, force_beta, force_r_unit


def _time_flux(lapse, shift_factor, radius, force_n, force_beta, force_r_unit, weight):
    """m_Q_t = G N [2 v F_N + (1+v²) F_β + v r F_r(N=1)].

    The second v F_N is the radial-pressure copy of the energy. Dropping it
    is a different flux. F_r(N=1) stands for (1/N) F_r(actual).
    """
    velocity = np.asarray(shift_factor, dtype=float)
    bracket = (
        2.0 * velocity * force_n
        + (1.0 + velocity * velocity) * force_beta
        + velocity * radius * force_r_unit
    )
    return NEWTON_G * float(weight) * np.asarray(lapse, dtype=float) * bracket


def reconstruct_geometry(grid, inner_mass, eps_n, eps_beta, angular, *, weight, kappa=SHELL_KAPPA,
                         force_r_unit=None):
    """Integrate m_Q from the inner stock and log N from N(r_out)=1.

    The cache is not iterated. v = β/N is the shift factor of the enclosed mass.
    """
    radius = grid["radius"]
    count = grid["points"]
    channels = {
        "eps_N": eps_n,
        "eps_beta": eps_beta,
        "angular": angular,
    }
    if force_r_unit is not None:
        channels["force_r_unit"] = force_r_unit
    force_n, force_beta, force_r_unit = _stress_channels(channels, radius)
    if not (force_n.shape == force_beta.shape == force_r_unit.shape == (count,)):
        raise ChartAdmissionError("matter densities do not match the radial grid")
    mass = np.empty(count)
    slope = np.empty(count)
    log_slope = np.empty(count)
    shift_factor = np.empty(count)
    mass[0] = float(inner_mass)
    if not math.isfinite(mass[0]):
        raise ChartAdmissionError("inner mass stock is not finite")
    slope[0], log_slope[0], shift_factor[0] = _mass_slope(
        mass[0], 0, grid, force_n, force_beta, weight,
    )
    for index in range(count - 1):
        step = float(radius[index + 1] - radius[index])
        guess = mass[index] + step * slope[index]
        for _attempt in range(12):
            next_slope, _log, _factor = _mass_slope(
                guess, index + 1, grid, force_n, force_beta, weight,
            )
            residual = guess - mass[index] - 0.5 * step * (slope[index] + next_slope)
            derivative = 1.0 - 0.5 * step * _mass_slope_mu_derivative(
                guess, index + 1, grid, force_beta, weight,
            )[0]
            if derivative == 0.0 or not math.isfinite(derivative):
                raise ChartAdmissionError("mass quadrature lost its Newton derivative")
            correction = residual / derivative
            guess -= correction
            if abs(correction) <= 1e-12 * (1.0 + abs(guess)):
                break
        else:
            raise ChartAdmissionError("mass quadrature did not converge")
        mass[index + 1] = guess
        slope[index + 1], log_slope[index + 1], shift_factor[index + 1] = _mass_slope(
            mass[index + 1], index + 1, grid, force_n, force_beta, weight,
        )
    log_lapse = np.zeros(count)
    for index in range(count - 2, -1, -1):
        step = float(radius[index + 1] - radius[index])
        log_lapse[index] = log_lapse[index + 1] - 0.5 * step * (
            log_slope[index] + log_slope[index + 1]
        )
    lapse = np.exp(log_lapse)
    if abs(float(lapse[-1]) - 1.0) > 1e-12:
        raise ChartAdmissionError("outer lapse normalization is not 1")
    shift = lapse * shift_factor
    q_mass_rate = _time_flux(
        lapse, shift_factor, radius, force_n, force_beta, force_r_unit, weight,
    )
    omitted = _time_flux(
        lapse, shift_factor, radius, 0.5 * force_n, force_beta, force_r_unit, weight,
    )
    charge2 = grid["charge"] ** 2
    defined = 0.5 * radius * shift_factor ** 2 + 0.5 * charge2 / radius
    hamiltonian = NEWTON_G * float(weight) * (force_n + shift_factor * force_beta)
    lapse_gradient = -NEWTON_G * float(weight) * force_beta / (radius * shift_factor)
    return {
        "mass": mass,
        "mass_slope": slope,
        "log_lapse_slope": log_slope,
        "shift_factor": shift_factor,
        "lapse": lapse,
        "shift": shift,
        "q_equation_mass_rate": q_mass_rate,
        "time_flux_without_second_v_fn": omitted,
        "second_v_fn_term": q_mass_rate - omitted,
        "force_n": force_n,
        "force_beta": force_beta,
        "force_r_unit": force_r_unit,
        "mass_definition_residual": mass - defined,
        "mass_slope_residual": slope - hamiltonian,
        "log_slope_residual": log_slope - lapse_gradient,
        "kinetic": force_n + radius * force_r_unit,
    }


def metric_jets(radius, lapse, shift):
    """Spatial nodal jets only. These arrays are not a dynamical R4."""
    spacing = float(radius[1] - radius[0])
    if np.max(np.abs(np.diff(radius) - spacing)) > 1e-9 * spacing:
        raise ChartAdmissionError("jet differentiation expects uniform areal spacing")
    lapse_r, lapse_rr = _differentiate(lapse, spacing)
    shift_r, shift_rr = _differentiate(shift, spacing)
    return {
        "lapse_r": lapse_r,
        "lapse_rr": lapse_rr,
        "shift_r": shift_r,
        "shift_rr": shift_rr,
        "jet_provenance": JET_PROVENANCE,
        "scope": JET_SCOPE,
        "spatial_only": True,
        "boundary": "polynomial_one_sided",
        "time_derivatives_included": False,
        "dynamic_R4_not_claimed": True,
        "static_R4_used": False,
    }


def admit_chart(grid, lapse, shift):
    """Pure outflow at excision, one incoming family at the outer boundary."""
    if np.any(lapse <= 0.0) or not np.all(np.isfinite(lapse)) or not np.all(np.isfinite(shift)):
        raise ChartAdmissionError("lapse must be positive and the shift finite")
    outgoing, ingoing = observables.chart_speeds(lapse, shift)
    if not (outgoing[0] < 0.0 and ingoing[0] < 0.0):
        raise ChartAdmissionError("excision is not pure outflow")
    if not (outgoing[-1] > 0.0 and ingoing[-1] < 0.0):
        raise ChartAdmissionError("outer boundary does not have exactly one incoming characteristic")
    cfl = observables.full_cfl_dt(lapse, shift, grid["spacing"])
    if cfl <= 0.0:
        raise ChartAdmissionError("full coupled CFL step is not positive")
    return {
        "admitted": True,
        "cfl_dt": cfl,
        "excision_pure_outflow": True,
        "outer_incoming_only": True,
        "outgoing_speed": outgoing,
        "ingoing_speed": ingoing,
        "blocker": None,
        "whole_field_set_to_zero": False,
        "components_clamped": False,
        "excision_outflow_clamp": False,
        "periodic": False,
    }


def _split(grid, coeff, component):
    operator = grid["derivative"]
    return 0.5 * (coeff * (operator @ component) + operator @ (coeff * component))


def dirac_rates(grid, phi, lapse, shift, *, kappa=SHELL_KAPPA):
    """Characteristic SBP evolution plus outer incoming SAT. Excision is outflow.

    U† σ1 U = σ2, so -iH sends χ+ to -α χ- and χ- to +α χ+. The ingoing
    angular term is therefore +α χ+.
    """
    array = observables.require_phi(phi, grid["points"])
    radius = grid["radius"]
    angular_scale = lapse * float(kappa) / radius
    rates = np.empty_like(array)
    for column in range(array.shape[2]):
        plus = array[0, :, column]
        minus = array[1, :, column]
        plus_rate = -_split(grid, lapse, plus) + _split(grid, shift, plus) - angular_scale * minus
        minus_rate = _split(grid, lapse, minus) + _split(grid, shift, minus) + angular_scale * plus
        incoming = -shift[-1] - lapse[-1]
        minus_rate[-1] += (incoming / grid["weights"][-1]) * minus[-1]
        rates[0, :, column] = plus_rate
        rates[1, :, column] = minus_rate
    sat_debit_rate = float(-(-shift[-1] - lapse[-1]) * np.sum(np.abs(array[1, -1, :]) ** 2))
    return rates, sat_debit_rate


def _density_variation(grid, phi, occupations, phi_rate, *, kappa=SHELL_KAPPA,
                       occupation_matrix=None):
    """Directional derivative of the cached forces at fixed occupations."""
    base = matter_densities(
        phi, grid, occupations, kappa=kappa, occupation_matrix=occupation_matrix,
    )
    step = 1e-7
    shifted = matter_densities(
        phi + step * phi_rate, grid, occupations, kappa=kappa,
        occupation_matrix=occupation_matrix,
    )
    return {
        "eps_N": (shifted["eps_N"] - base["eps_N"]) / step,
        "eps_beta": (shifted["eps_beta"] - base["eps_beta"]) / step,
        "angular": (shifted["angular"] - base["angular"]) / step,
    }


def _probability_boundary(phi, lapse, shift):
    """Coordinate probability current. Mass flux and the SAT debit are not this."""
    outgoing = -np.asarray(shift, dtype=float) + np.asarray(lapse, dtype=float)
    ingoing = -np.asarray(shift, dtype=float) - np.asarray(lapse, dtype=float)
    current = (
        outgoing * np.sum(np.abs(phi[0]) ** 2, axis=-1)
        + ingoing * np.sum(np.abs(phi[1]) ** 2, axis=-1)
    )
    return {
        "excision_outflow": float(-current[0]),
        "outer_outflow": float(current[-1]),
    }


def _linearised_rates(grid, geometry, deps_n, deps_beta, inner_rate, weight):
    """∂t of the reconstructed mass and lapse.

    N(r_out)=1 is held, so N_t(r_out)=0 is that boundary constraint. It is
    not a free choice to freeze the lapse.
    """
    count = grid["points"]
    radius = grid["radius"]
    charge2 = grid["charge"] ** 2
    dmu = np.zeros(count)
    df = np.empty(count)
    dg = np.empty(count)
    dmu[0] = float(inner_rate)
    for index in range(count):
        shift_factor = float(geometry["shift_factor"][index])
        mass_value = float(geometry["mass"][index])
        enclosed = 2.0 * mass_value - charge2 / float(radius[index])
        coefficient = NEWTON_G * weight
        eps_beta = float(geometry["eps_beta_bare"][index])
        log_slope = -coefficient * eps_beta / (float(radius[index]) * shift_factor)
        dlog_deps = -coefficient * float(deps_beta[index]) / (float(radius[index]) * shift_factor)
        dlog_dmu = coefficient * eps_beta / (float(radius[index]) ** 2 * shift_factor ** 3)
        inhom_f = coefficient * float(deps_n[index]) - dlog_deps * enclosed
        df_dmu = -dlog_dmu * enclosed - 2.0 * log_slope
        if index == 0:
            df[0] = inhom_f + df_dmu * dmu[0]
            dg[0] = dlog_deps + dlog_dmu * dmu[0]
            continue
        step = float(radius[index] - radius[index - 1])
        denominator = 1.0 - 0.5 * step * df_dmu
        if denominator == 0.0 or not math.isfinite(denominator):
            raise ChartAdmissionError("linearised mass rate is singular")
        dmu[index] = (
            dmu[index - 1] + 0.5 * step * (df[index - 1] + inhom_f)
        ) / denominator
        df[index] = inhom_f + df_dmu * dmu[index]
        dg[index] = dlog_deps + dlog_dmu * dmu[index]
    dlog_lapse = np.zeros(count)
    for index in range(count - 2, -1, -1):
        step = float(radius[index + 1] - radius[index])
        dlog_lapse[index] = dlog_lapse[index + 1] - 0.5 * step * (dg[index] + dg[index + 1])
    lapse_rate = geometry["lapse"] * dlog_lapse
    shift_rate = lapse_rate * geometry["shift_factor"] + geometry["lapse"] * dmu / (
        radius * geometry["shift_factor"]
    )
    return dmu, lapse_rate, shift_rate


def _frequency_metadata(killing_frequency):
    """Positive scalar recorded as metadata. It does not prove a Killing mode."""
    try:
        frequency = float(killing_frequency)
    except (TypeError, ValueError) as error:
        raise ValueError("Killing frequency must be a positive finite scalar") from error
    if not math.isfinite(frequency) or frequency <= 0.0:
        raise ValueError(
            "Killing frequency must be a positive finite scalar; "
            "metadata only, not an exterior-mode proof"
        )
    return {
        "value": frequency,
        "metadata_only": True,
        "proves_exterior_killing_frequency": False,
    }


def _probability_norm(phi, weights):
    weights = np.asarray(weights, dtype=float)
    density = np.sum(np.abs(np.asarray(phi)) ** 2, axis=0)
    if density.ndim == 2:
        density = np.sum(density, axis=1)
    return float(np.real(np.sum(weights * density)))


def _spd_gram(gram):
    hermitian = 0.5 * (gram + gram.conj().T)
    eigenvalues = np.linalg.eigvalsh(hermitian)
    return bool(np.min(np.real(eigenvalues)) > 1e-14)


def _occupation_record(phi, weights, occupations, covariance):
    """Weights or a coefficient matrix. Stored covariance is the physical CAR.

    Multiplicity and rank are not divisors of that matrix.
    """
    rank = phi.shape[2]
    if occupations is None and covariance is None:
        occupations = np.zeros(rank)
    coefficients = None
    if occupations is not None:
        supplied = np.asarray(occupations)
        if supplied.ndim == 1:
            weights_c = np.asarray(supplied, dtype=float)
            if weights_c.shape != (rank,) or np.any(weights_c < 0.0) or np.any(weights_c >= 1.0):
                raise ChartAdmissionError("Gaussian occupations must lie in [0, 1)")
            if not np.all(np.isfinite(weights_c)):
                raise ChartAdmissionError("Gaussian occupations must be finite")
            coefficients = np.diag(weights_c.astype(complex))
        elif supplied.ndim == 2:
            coefficients = np.asarray(supplied, dtype=complex)
            if coefficients.shape != (rank, rank) or not np.all(np.isfinite(coefficients)):
                raise ChartAdmissionError("weighted columns must be a finite rank-by-rank matrix")
            if not np.allclose(coefficients, coefficients.conj().T, atol=1e-10):
                raise ChartAdmissionError("weighted columns must be Hermitian")
            spectrum = np.linalg.eigvalsh(0.5 * (coefficients + coefficients.conj().T))
            if np.min(spectrum) < -1e-10 or np.max(spectrum) >= 1.0:
                raise ChartAdmissionError("occupation eigenvalues must lie in [0, 1)")
            weights_c = np.real(np.diag(coefficients))
        else:
            raise ChartAdmissionError("occupations must be weights or weighted columns")
    gram = source.weighted_gram(phi, weights)
    if coefficients is None:
        car = np.asarray(covariance, dtype=complex)
        if car.shape != (rank, rank) or not np.all(np.isfinite(car)):
            raise ChartAdmissionError("physical CAR must be a finite rank-by-rank matrix")
        if not np.allclose(car, car.conj().T, atol=1e-10):
            raise ChartAdmissionError("physical CAR must be Hermitian")
        if not np.allclose(gram, np.eye(rank), atol=1e-8):
            raise ChartAdmissionError("physical CAR without weights requires orthonormal columns")
        coefficients = car
        weights_c = np.real(np.diag(coefficients))
        physical = 0.5 * (car + car.conj().T)
        car_defined = True
    elif _spd_gram(gram):
        physical = np.asarray(source.physical_car(gram, coefficients), dtype=complex)
        car_defined = True
        if covariance is not None:
            supplied_car = np.asarray(covariance, dtype=complex)
            if supplied_car.shape != physical.shape or np.max(np.abs(supplied_car - physical)) > 1e-8:
                raise ChartAdmissionError("supplied covariance is not the physical CAR")
    else:
        physical = np.array(coefficients, dtype=complex, copy=True)
        car_defined = False
    return {
        "occupations": np.array(weights_c, dtype=float, copy=True),
        "occupation_matrix": None if occupations is not None and np.ndim(occupations) == 1
        else np.array(coefficients, dtype=complex, copy=True),
        "covariance": physical,
        "covariance_kind": "physical_car" if car_defined else "coefficients_gram_singular",
        "covariance_is_physical_car": car_defined,
        "multiplicity_in_covariance": False,
        "covariance_divided_by_rank": False,
        "multiplicity_over_rank_used": False,
    }


def make_state(grid, phi=None, *, occupations=None, covariance=None, inner_mass=None,
               killing_frequency=1.0):
    """Canonical φ, fixed weights or weighted columns, inner mass and clocks.

    The stored covariance is the physical CAR when the Gram is positive
    definite. It is not divided by multiplicity or rank. The shell factor 4
    is applied once later, inside the source. A positive frequency scalar is
    metadata, not a proof that φ is an exterior Killing mode.
    """
    if phi is None:
        phi = np.zeros((2, grid["points"], 1), dtype=complex)
    array = observables.require_phi(phi, grid["points"])
    if "weights" not in grid:
        raise ChartAdmissionError("make_state needs the SBP weights")
    record = _occupation_record(array, grid["weights"], occupations, covariance)
    _shell_weight()
    frequency = _frequency_metadata(killing_frequency)
    mass0 = grid["mass"] if inner_mass is None else float(inner_mass)
    return {
        "phi": np.array(array, dtype=complex, copy=True),
        "occupations": record["occupations"],
        "occupation_matrix": record["occupation_matrix"],
        "covariance": record["covariance"],
        "covariance_kind": record["covariance_kind"],
        "covariance_is_physical_car": record["covariance_is_physical_car"],
        "multiplicity_in_covariance": False,
        "covariance_divided_by_rank": False,
        "multiplicity_over_rank_used": False,
        "occupations_fixed": True,
        "inner_mass": mass0,
        "charge": grid["charge"],
        "killing_frequency": frequency["value"],
        "killing_frequency_metadata": frequency,
        "multiplicity": int(_shell_weight()),
        "kappa": SHELL_KAPPA,
        "clocks": {"t": 0.0, "killing": 0.0, "normal": 0.0},
        "stocks": {"excision_flux": 0.0, "outer_flux": 0.0, "sat_debit": 0.0},
        "probability_stocks": {
            "remaining": _probability_norm(array, grid["weights"]),
            "excision_outflow": 0.0,
            "outer_outflow": 0.0,
        },
        "electric_current": 0.0,
        "magnetic_flux_fixed": True,
        "flux_fitted": False,
    }


def _copy_state(state, phi, inner_mass, clocks, stocks, probability):
    matrix = state["occupation_matrix"]
    return {
        "phi": np.array(phi, dtype=complex, copy=True),
        "occupations": np.array(state["occupations"], dtype=float, copy=True),
        "occupation_matrix": None if matrix is None else np.array(matrix, dtype=complex, copy=True),
        "covariance": np.array(state["covariance"], dtype=complex, copy=True),
        "covariance_kind": state["covariance_kind"],
        "covariance_is_physical_car": state["covariance_is_physical_car"],
        "multiplicity_in_covariance": False,
        "covariance_divided_by_rank": False,
        "multiplicity_over_rank_used": False,
        "occupations_fixed": True,
        "inner_mass": float(inner_mass),
        "charge": state["charge"],
        "killing_frequency": state["killing_frequency"],
        "killing_frequency_metadata": dict(state["killing_frequency_metadata"]),
        "multiplicity": state["multiplicity"],
        "kappa": state["kappa"],
        "clocks": dict(clocks),
        "stocks": dict(stocks),
        "probability_stocks": dict(probability),
        "electric_current": 0.0,
        "magnetic_flux_fixed": True,
        "flux_fitted": False,
    }


def stage_rates(state, grid, *, densities=None):
    """One right-hand side. Geometry is reconstructed before the Dirac rate."""
    array = observables.require_phi(state["phi"], grid["points"])
    # Source forces already include the shell factor once. Do not iterate them.
    weight = 1.0
    matrix = state.get("occupation_matrix")
    if densities is None:
        densities = matter_densities(
            array, grid, state["occupations"], kappa=state["kappa"],
            occupation_matrix=matrix,
        )
    else:
        densities = dict(densities)
        densities.setdefault("source", "injected_bare_densities")
        densities.setdefault("electric_current", np.zeros(grid["points"]))
        densities.setdefault("angular", np.zeros(grid["points"]))
    if np.any(np.asarray(densities["electric_current"]) != 0.0):
        raise ChartAdmissionError("electric current is not part of the first matter shell")
    geometry = reconstruct_geometry(
        grid,
        state["inner_mass"],
        densities["eps_N"],
        densities["eps_beta"],
        densities["angular"],
        weight=weight,
        kappa=state["kappa"],
        force_r_unit=densities.get("force_r_unit"),
    )
    geometry["eps_beta_bare"] = np.asarray(densities["eps_beta"], dtype=float)
    admission = admit_chart(grid, geometry["lapse"], geometry["shift"])
    phi_rate, sat_rate = dirac_rates(
        grid, array, geometry["lapse"], geometry["shift"], kappa=state["kappa"],
    )
    if densities.get("source") == "injected_bare_densities":
        deps = {
            "eps_N": np.zeros(grid["points"]),
            "eps_beta": np.zeros(grid["points"]),
            "angular": np.zeros(grid["points"]),
        }
    else:
        deps = _density_variation(
            grid, array, state["occupations"], phi_rate, kappa=state["kappa"],
            occupation_matrix=matrix,
        )
    inner_rate = float(geometry["q_equation_mass_rate"][0])
    mass_rate, lapse_rate, shift_rate = _linearised_rates(
        grid, geometry, deps["eps_N"], deps["eps_beta"], inner_rate, weight,
    )
    if abs(float(lapse_rate[-1])) > 1e-9:
        raise ChartAdmissionError("N_t at the outer boundary is fixed by N=1")
    jets = metric_jets(grid["radius"], geometry["lapse"], geometry["shift"])
    probability = _probability_boundary(array, geometry["lapse"], geometry["shift"])
    time_jets = {
        "lapse_t": lapse_rate,
        "shift_t": shift_rate,
        "mass_t": mass_rate,
        "provenance": TIME_JET_PROVENANCE,
        "static_R4_used": False,
        "spatial_jets_used": False,
    }
    outgoing = admission["outgoing_speed"]
    ingoing = admission["ingoing_speed"]
    killing_rate = 0.0
    metric_tt = geometry["lapse"][-1] ** 2 - geometry["shift"][-1] ** 2
    if metric_tt > 0.0:
        killing_rate = math.sqrt(metric_tt)
    balance = mass_rate - geometry["q_equation_mass_rate"]
    misner = observables.misner_sharp(grid["radius"], geometry["lapse"], geometry["shift"])
    return {
        "phi_rate": phi_rate,
        "inner_mass_rate": inner_rate,
        "mass_rate": mass_rate,
        "q_equation_mass_rate": geometry["q_equation_mass_rate"],
        "balance_residual": balance,
        "lapse": geometry["lapse"],
        "shift": geometry["shift"],
        "mass": geometry["mass"],
        "misner_sharp": misner,
        "lapse_rate": lapse_rate,
        "shift_rate": shift_rate,
        "mass_slope": geometry["mass_slope"],
        "log_lapse_slope": geometry["log_lapse_slope"],
        "jets": jets,
        "time_jets": time_jets,
        "admission": admission,
        "cfl_dt": admission["cfl_dt"],
        "sat_debit_rate": sat_rate,
        "outer_flux_rate": float(geometry["q_equation_mass_rate"][-1]),
        "excision_flux_rate": inner_rate,
        "probability_excision_outflow_rate": probability["excision_outflow"],
        "probability_outer_outflow_rate": probability["outer_outflow"],
        "probability_remaining": _probability_norm(array, grid["weights"]),
        "dynamic_contract": {
            "mass_definition_residual": float(np.max(np.abs(geometry["mass_definition_residual"]))),
            "mass_slope_residual": float(np.max(np.abs(geometry["mass_slope_residual"]))),
            "log_slope_residual": float(np.max(np.abs(geometry["log_slope_residual"]))),
            "second_v_fn_term_max": float(np.max(np.abs(geometry["second_v_fn_term"]))),
            "time_flux_uses_second_v_fn": True,
            "static_R4_used": False,
            "intrinsic_stress_cached_at_unit_lapse": bool(
                densities.get("intrinsic_stress_cached_at_unit_lapse", False)
            ),
            "stress_iteration": False,
            "q_adm": 1.0,
        },
        "time_derivative_constraints": {
            "N_outer": 1.0,
            "N_t_outer": float(lapse_rate[-1]),
            "N_t_outer_fixed_by_normalization": True,
            "interior_lapse_not_held_at_one": True,
        },
        "killing_rate": killing_rate,
        "normal_rate": float(geometry["lapse"][0]),
        "shell_factor": weight,
        "electric_current": 0.0,
        "magnetic_flux_fixed": True,
        "density_source": densities.get("source"),
        "outgoing_speed": outgoing,
        "ingoing_speed": ingoing,
        "jet_provenance": jets["jet_provenance"],
        "stationary_candidate": False,
        "weyl_pole": False,
        "fitted_restoring_force": False,
        "lapse_held_at_one": False,
        "flux_fitted": False,
        "G_N": float(grid.get("G_N", NEWTON_G)),
        "integer_flux": ACTION_FLUX,
    }


def _advance(state, rate, step, grid):
    phi = state["phi"] + step * rate["phi_rate"]
    clocks = {
        "t": state["clocks"]["t"] + step * 1.0,
        "killing": state["clocks"]["killing"] + step * rate["killing_rate"],
        "normal": state["clocks"]["normal"] + step * rate["normal_rate"],
    }
    stocks = {
        "excision_flux": state["stocks"]["excision_flux"] + step * rate["excision_flux_rate"],
        "outer_flux": state["stocks"]["outer_flux"] + step * rate["outer_flux_rate"],
        "sat_debit": state["stocks"]["sat_debit"] + step * rate["sat_debit_rate"],
    }
    probability = {
        "remaining": _probability_norm(phi, grid["weights"]),
        "excision_outflow": (
            state["probability_stocks"]["excision_outflow"]
            + step * rate["probability_excision_outflow_rate"]
        ),
        "outer_outflow": (
            state["probability_stocks"]["outer_outflow"]
            + step * rate["probability_outer_outflow_rate"]
        ),
    }
    return _copy_state(
        state,
        phi,
        state["inner_mass"] + step * rate["inner_mass_rate"],
        clocks,
        stocks,
        probability,
    )


def rk4_step(state, grid, dt, *, densities=None):
    """Four reconstructions. dt must respect the full characteristic CFL."""
    step = float(dt)
    if not math.isfinite(step) or step <= 0.0:
        raise ChartAdmissionError("RK4 step must be positive")
    first = stage_rates(state, grid, densities=densities)
    if step > first["cfl_dt"] * (1.0 + 1e-12):
        raise ChartAdmissionError("step exceeds the full coupled CFL")
    if densities is not None:
        # Injected sources do not travel with φ; every stage sees the same bare density.
        stage = lambda current: stage_rates(current, grid, densities=densities)
    else:
        stage = lambda current: stage_rates(current, grid)
    second = stage(_advance(state, first, 0.5 * step, grid))
    third = stage(_advance(state, second, 0.5 * step, grid))
    fourth = stage(_advance(state, third, step, grid))
    clocks = {
        "t": state["clocks"]["t"] + step,
        "killing": state["clocks"]["killing"] + (step / 6.0) * (
            first["killing_rate"] + 2.0 * second["killing_rate"]
            + 2.0 * third["killing_rate"] + fourth["killing_rate"]
        ),
        "normal": state["clocks"]["normal"] + (step / 6.0) * (
            first["normal_rate"] + 2.0 * second["normal_rate"]
            + 2.0 * third["normal_rate"] + fourth["normal_rate"]
        ),
    }
    stocks = {
        "excision_flux": state["stocks"]["excision_flux"] + (step / 6.0) * (
            first["excision_flux_rate"] + 2.0 * second["excision_flux_rate"]
            + 2.0 * third["excision_flux_rate"] + fourth["excision_flux_rate"]
        ),
        "outer_flux": state["stocks"]["outer_flux"] + (step / 6.0) * (
            first["outer_flux_rate"] + 2.0 * second["outer_flux_rate"]
            + 2.0 * third["outer_flux_rate"] + fourth["outer_flux_rate"]
        ),
        "sat_debit": state["stocks"]["sat_debit"] + (step / 6.0) * (
            first["sat_debit_rate"] + 2.0 * second["sat_debit_rate"]
            + 2.0 * third["sat_debit_rate"] + fourth["sat_debit_rate"]
        ),
    }
    phi = state["phi"] + (step / 6.0) * (
        first["phi_rate"] + 2.0 * second["phi_rate"] + 2.0 * third["phi_rate"] + fourth["phi_rate"]
    )
    inner = state["inner_mass"] + (step / 6.0) * (
        first["inner_mass_rate"] + 2.0 * second["inner_mass_rate"]
        + 2.0 * third["inner_mass_rate"] + fourth["inner_mass_rate"]
    )
    probability = {
        "remaining": _probability_norm(phi, grid["weights"]),
        "excision_outflow": state["probability_stocks"]["excision_outflow"] + (step / 6.0) * (
            first["probability_excision_outflow_rate"] + 2.0 * second["probability_excision_outflow_rate"]
            + 2.0 * third["probability_excision_outflow_rate"] + fourth["probability_excision_outflow_rate"]
        ),
        "outer_outflow": state["probability_stocks"]["outer_outflow"] + (step / 6.0) * (
            first["probability_outer_outflow_rate"] + 2.0 * second["probability_outer_outflow_rate"]
            + 2.0 * third["probability_outer_outflow_rate"] + fourth["probability_outer_outflow_rate"]
        ),
    }
    return _copy_state(state, phi, inner, clocks, stocks, probability)


def sourcefree_rates(mass, *, r_m=1.0, points=33, killing_frequency=1.0, r_out=None):
    """Vacuum RN on this chart. Rates vanish and β = N √(2M/r − P²/r²)."""
    grid = build_grid(mass, r_m=r_m, points=points, r_out=r_out)
    model = grid["model"]
    shift_target = np.asarray(model.beta_rn(grid["radius"]), dtype=float)
    lapse_target = np.asarray(model.vacuum_lapse(grid["radius"]), dtype=float)
    state = make_state(grid, killing_frequency=killing_frequency, inner_mass=grid["mass"])
    rates = stage_rates(state, grid)
    return grid, state, rates, lapse_target, shift_target


def geometry_for_observer(rates, grid, *, stationary):
    """Spatial jets for a static check. They do not carry a dynamical R4."""
    jets = rates["jets"]
    if jets.get("spatial_only") is not True or jets.get("time_derivatives_included") is not False:
        raise ChartAdmissionError("observer geometry was given jets that claim time derivatives")
    return {
        "radius": grid["radius"],
        "lapse": rates["lapse"],
        "shift": rates["shift"],
        "lapse_r": jets["lapse_r"],
        "shift_r": jets["shift_r"],
        "lapse_rr": jets["lapse_rr"],
        "shift_rr": jets["shift_rr"],
        "jet_provenance": jets["jet_provenance"],
        "scope": jets.get("scope", JET_SCOPE),
        "spatial_only": True,
        "time_derivatives_included": False,
        "dynamic_R4_from_these_jets": False,
        "static_R4_used_as_time_jet": False,
        "electric_current": 0,
        "magnetic_flux_fixed": True,
        "stationary": bool(stationary),
        "coupled_lapse_replaced": False,
    }


def directional_metric_check(state, grid, epsilon=1.0e-5):
    """Central difference of one stage against the constraint linearization.

    Static R4 is not an input. Spatial jets are not differentiated in time.
    """
    step = float(epsilon)
    if not math.isfinite(step) or step <= 0.0:
        raise ChartAdmissionError("directional check needs a positive step")
    base = stage_rates(state, grid)

    def _probe(sign):
        phi = state["phi"] + (sign * step) * base["phi_rate"]
        inner = state["inner_mass"] + (sign * step) * base["inner_mass_rate"]
        clocks = state["clocks"]
        stocks = state["stocks"]
        probability = state["probability_stocks"]
        probed = _copy_state(state, phi, inner, clocks, stocks, probability)
        return stage_rates(probed, grid)

    ahead = _probe(1.0)
    behind = _probe(-1.0)
    lapse_fd = (ahead["lapse"] - behind["lapse"]) / (2.0 * step)
    shift_fd = (ahead["shift"] - behind["shift"]) / (2.0 * step)
    mass_fd = (ahead["mass"] - behind["mass"]) / (2.0 * step)
    return {
        "lapse_residual": float(np.max(np.abs(lapse_fd - base["lapse_rate"]))),
        "shift_residual": float(np.max(np.abs(shift_fd - base["shift_rate"]))),
        "mass_residual": float(np.max(np.abs(mass_fd - base["mass_rate"]))),
        "provenance": TIME_JET_PROVENANCE,
        "static_R4_used": False,
        "spatial_jets_used_as_time_derivatives": False,
    }
