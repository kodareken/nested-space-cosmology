"""Machine-readable evidence for the version-0.9 controlled model.

The record assembled here joins six deliberately limited artifacts: the
CCT-1 closed-de-Sitter core, the IM-1 finite code isometry, classical linear
mode transfer, the DST-1 scalar/Kottler effective-limit synthesis, and a
curvature-sign discriminator conditional on an added globally closed-child
branch postulate, plus conditional literal-subclass rejection gates.  The
published Gaussian summaries in the last artifact are not rerun likelihoods.
Nothing here is evidence for a black-hole-to-child embedding or an observation.
"""

from __future__ import annotations

from math import cos, pi, sin, sqrt
from typing import Any

from .core import C, SOLAR_MASS
from .domain import (
    DomainScalarSpec,
    exchange_ledger,
    integer_damped_period_benchmark,
    kottler_balance_radius,
    kottler_radial_acceleration,
    monomial_rapid_oscillation,
    quadratic_cycle_average,
)
from .falsification import (
    DESI_DR2_FLAT_WCDM_SUMMARIES,
    BoundaryScalingSpec,
    approximate_gaussian_w_gate,
    boundary_acceleration_gate,
    disformal_photon_cone,
    gw170817_cone_gate,
    tensor_photon_delta,
)
from .information import (
    TransitionCodeSpec,
    analytic_child_entropy,
    analytic_radiation_entropy,
    apply_isometry,
    no_cloning_overlap_contrast,
    partial_trace_child,
    partial_trace_radiation,
    state_density,
    validate_transition,
    von_neumann_entropy,
)
from .perturbations import core_transfer_matrix, transfer_convergence
from .transition import (
    apparent_horizon_chi,
    classify_surface,
    closed_curvature_discriminator,
    curvature_invariants,
    friedmann_residual,
    hubble_derivative,
    hubble_parameter,
    scale_factor,
)


def _matrix_record(matrix: object) -> dict[str, float]:
    """Serialize a perturbation ``TransferMatrix`` without exposing internals."""

    return {
        "m11": float(getattr(matrix, "m11")),
        "m12": float(getattr(matrix, "m12")),
        "m21": float(getattr(matrix, "m21")),
        "m22": float(getattr(matrix, "m22")),
        "determinant": float(getattr(matrix, "determinant")),
    }


def _published_w_summary_gates(proposed_w: float) -> dict[str, dict[str, object]]:
    """Compare a fixed literal scaling to published DESI summaries only.

    These are deliberately independent one-dimensional Gaussian summaries.  They
    are not combined, and they do not reconstruct a DESI, CMB, or supernova
    likelihood.
    """

    return {
        summary.label: dict(approximate_gaussian_w_gate(proposed_w, summary))
        for summary in DESI_DR2_FLAT_WCDM_SUMMARIES
    }


def _literal_boundary_scaling_record(q: float, s: float) -> dict[str, object]:
    """Return one separately conserved literal-boundary subclass test vector."""

    spec = BoundaryScalingSpec(q=q, s=s, separately_conserved=True)
    acceleration = dict(boundary_acceleration_gate(spec))
    effective_w = float(acceleration["w_effective"])
    return {
        "q": q,
        "s": s,
        "separately_conserved": True,
        "density_scale_exponent": acceleration["density_scale_exponent"],
        "effective_w": effective_w,
        "accelerates": acceleration["accelerates"],
        "acceleration_classification": acceleration["classification"],
        "desi_dr2_flat_wcdm_summary_gates": _published_w_summary_gates(effective_w),
    }


def reproduce_controlled_model() -> dict[str, Any]:
    """Return the complete JSON-safe v0.11 controlled-model evidence record."""

    length = 1.0
    sample_times = (-2.0, -1.0, 0.0, 1.0, 2.0)
    residuals = [abs(friedmann_residual(time, length)) for time in sample_times]
    invariants = curvature_invariants(length)

    # Three fixed, predeclared examples exercise the actual null-expansion
    # classifier.  They are not a statistical sample or observational data.
    regions = {
        "contracting_equator_t_over_L_minus_2": classify_surface(
            -2.0, pi / 2.0, length
        ),
        "bounce_quarter_sphere": classify_surface(0.0, pi / 4.0, length),
        "expanding_equator_t_over_L_plus_2": classify_surface(
            2.0, pi / 2.0, length
        ),
    }
    example_time = 1.0
    example_scale = scale_factor(example_time, length)
    example_hubble = hubble_parameter(example_time, length)

    spec = TransitionCodeSpec(d=2, p=0.5)
    input_density = state_density([1.0 / sqrt(2.0), 1j / sqrt(2.0)])
    output_density = apply_isometry(spec, input_density)
    radiation = partial_trace_child(spec, output_density)
    child = partial_trace_radiation(spec, output_density)
    input_entropy = von_neumann_entropy(input_density)

    scalar_eta = 0.7
    scalar_n = 2
    scalar_mass_radius = sqrt(2.0)
    scalar = core_transfer_matrix(
        sector="scalar",
        n=scalar_n,
        eta_core=scalar_eta,
        mass_radius=scalar_mass_radius,
        steps=4_096,
    )
    scalar_frequency = scalar_n + 1.0
    scalar_elapsed = 2.0 * scalar_eta
    scalar_exact = (
        cos(scalar_frequency * scalar_elapsed),
        sin(scalar_frequency * scalar_elapsed) / scalar_frequency,
        -scalar_frequency * sin(scalar_frequency * scalar_elapsed),
        cos(scalar_frequency * scalar_elapsed),
    )
    scalar_error = max(
        abs(value - expected)
        for value, expected in zip(
            (scalar.m11, scalar.m12, scalar.m21, scalar.m22),
            scalar_exact,
            strict=True,
        )
    )

    tensor_records = transfer_convergence(
        sector="tensor",
        n=3,
        eta_core=0.7,
        initial_steps=128,
        refinements=5,
    )
    tensor_final = tensor_records[-1]

    domain_spec = DomainScalarSpec(
        vacuum_energy=2.0,
        mass=10.0,
        decay_rate=0.4,
    )
    cycle_average = quadratic_cycle_average(domain_spec, amplitude=2.0)
    quadratic_monomial = monomial_rapid_oscillation(2.0)
    spectator_endpoint = integer_damped_period_benchmark(
        domain_spec,
        amplitude=2.0,
        hubble=0.2,
        periods=9,
    )
    transfer_ledger = exchange_ledger(
        domain_spec,
        displacement=2.0,
        velocity=3.0,
        hubble=0.4,
        radiation_density=7.0,
    )
    kottler_lambda = 1.0e-52
    balance_radius = kottler_balance_radius(SOLAR_MASS, kottler_lambda)
    balance_acceleration = kottler_radial_acceleration(
        SOLAR_MASS,
        balance_radius,
        kottler_lambda,
    )

    return {
        "schema_version": 3,
        "project_version": "0.11.0",
        "scope": (
            "Controlled closed-de-Sitter core, finite code-subspace isometry, "
            "classical linear-mode transfer, scalar effective-fluid and Kottler "
            "balance benchmarks, and a sign discriminator conditional on an added "
            "globally closed-child postulate, plus fixed relative-cone and boundary-"
            "scaling rejection gates; no black-hole embedding, dark-sector "
            "identification, varying-constant detection, likelihood rerun, or empirical confirmation."
        ),
        "covariant_core": {
            "model_id": "CCT-1",
            "classification": "controlled_covariant_core_not_global_black_hole_transition",
            "action": "Einstein_gravity_plus_canonical_scalar_at_positive_stationary_point",
            "matter_condition": "ordinary_matter_absent_vacuum_energy_absorbed_into_V0",
            "units": "c=hbar=M_Pl=L=1_for_machine_benchmark",
            "a_at_bounce": scale_factor(0.0, length),
            "H_at_bounce": hubble_parameter(0.0, length),
            "Hdot_at_bounce": hubble_derivative(0.0, length),
            "max_abs_friedmann_residual": max(residuals),
            "curvature_invariants": invariants,
            "surface_examples": regions,
            "apparent_horizons_at_t_over_L_1": apparent_horizon_chi(
                example_time, length
            ),
            "logical_limit": (
                "Global de Sitter in closed slicing has constant curvature; "
                "the scale-factor minimum is not derived transition microphysics."
            ),
        },
        "information_map": {
            "model_id": "IM-1",
            "classification": "finite_code_subspace_mathematical_demonstration",
            "code_dimension": spec.d,
            "p": spec.p,
            "p_status": "chosen_symmetric_benchmark_not_derived",
            "audit": validate_transition(spec, input_density),
            "input_entropy_nats": input_entropy,
            "joint_output_entropy_nats": von_neumann_entropy(output_density),
            "radiation_entropy_nats": von_neumann_entropy(radiation),
            "child_entropy_nats": von_neumann_entropy(child),
            "analytic_radiation_entropy_nats": analytic_radiation_entropy(
                spec, input_entropy
            ),
            "analytic_child_entropy_nats": analytic_child_entropy(
                spec, input_entropy
            ),
            "no_cloning_overlap_contrast": no_cloning_overlap_contrast(),
            "logical_limit": (
                "Assumes a finite tensor-factor code approximation; does not "
                "derive gravitational factorization, dynamics, or H-SAT."
            ),
        },
        "linear_perturbations": {
            "classification": "classical_core_transfer_not_primordial_spectrum",
            "scalar_exact_regression": {
                "n": scalar_n,
                "mass_times_L": scalar_mass_radius,
                "eta_core": scalar_eta,
                "steps": 4_096,
                "matrix": _matrix_record(scalar),
                "max_abs_error_from_exact_rotation": scalar_error,
            },
            "tensor_refinement": {
                "n": 3,
                "eta_core": 0.7,
                "steps": [record.steps for record in tensor_records],
                "final_matrix": _matrix_record(tensor_final.matrix),
                "final_determinant_error": tensor_final.determinant_error,
                "final_difference_from_previous": tensor_final.difference_from_previous,
            },
            "logical_limit": (
                "No adiabatic in-state, curvature-mode conversion, reheating, "
                "n_s, A_s, Bogoliubov interpretation, or likelihood is supplied."
            ),
        },
        "domain_synthesis": {
            "model_id": "DST-1",
            "classification": (
                "controlled_effective_fluid_and_kottler_balance_benchmarks_"
                "not_dark_sector_or_recursive_transition"
            ),
            "framework_context": (
                "DSF-1 is an operational hypothesis about dimensionless domain "
                "settings and relative characteristic cones, not a computed "
                "varying-constant result."
            ),
            "scalar_potential": "V=V0+m^2*varphi^2/2",
            "quadratic_cycle_average": cycle_average,
            "quadratic_monomial_limit": quadratic_monomial,
            "fixed_h_spectator_endpoint": spectator_endpoint,
            "exchange_ledger": transfer_ledger,
            "kottler_balance": {
                "mass_kg": SOLAR_MASS,
                "lambda_m^-2": kottler_lambda,
                "balance_radius_m": balance_radius,
                "acceleration_at_balance_m_s^-2": balance_acceleration,
                "radial_slope_at_balance_s^-2": kottler_lambda * C * C,
                "equation": "a_r=-G*M/r^2+Lambda*c^2*r/3",
                "stability": "unstable_radial_separator",
                "interpretation": (
                    "Weak-field competition between compact-mass attraction and "
                    "positive-Lambda vacuum acceleration at an unstable sign-change "
                    "scale; not a stable equilibrium, not an accretion wall, and "
                    "not evidence for recursive domains."
                ),
            },
            "logical_limit": (
                "The background identities do not derive the observed dark-energy "
                "scale, dark-matter abundance or clustering, microscopic reheating, "
                "dimensionless setting variation, or a global parent-child solution."
            ),
        },
        "independent_discriminator": {
            "id": "D-KSIGN-1",
            "classification": "conditional_sign_prediction_for_strict_closed_branch",
            "quantity": "Omega_K=-K*c^2/(a^2*H^2)",
            "prediction": "Omega_K<0_when_H_is_nonzero",
            "required_branch_postulate": "observable_child_remains_globally_closed_with_K=+1",
            "illustrative_t_over_L": example_time,
            "illustrative_omega_K_not_present_day_prediction": closed_curvature_discriminator(
                example_scale, example_hubble
            ),
            "independence": (
                "Conditional on the added K=+1 global branch postulate, the sign "
                "is independent of H-SAT and of the curvature data that would test "
                "it; CCT-1 alone does not derive the post-core matching."
            ),
            "falsifier": "robust_positive_Omega_K_across_stated_models_and_systematics",
            "nonconfirmation": (
                "An observational result statistically consistent with zero does "
                "not confirm the branch; later expansion can dilute its magnitude "
                "and other closed models share the sign."
            ),
        },
        "literal_model_gates": {
            "model_id": "DSF-2_DEB-1",
            "classification": (
                "conditional_literal_subclass_rejection_gates_not_model_confirmation"
            ),
            "scope": (
                "Fixed control and rejection vectors for a late-universe relative "
                "tensor/photon cone and separately conserved literal boundary "
                "scalings. Published DESI DR2 values are one-dimensional Gaussian "
                "posterior summaries only, not a rerun or combination of likelihoods."
            ),
            "source_caveat": (
                "The GW170817 interval is a low-redshift common-path relative-cone "
                "constraint subject to its source-emission-delay assumptions. DESI "
                "summaries depend on the stated flat constant-w CDM datasets, priors, "
                "epoch, fiducial sum_mnu=0.06_eV baseline, and parameterization; "
                "they do not test an arbitrary dark sector."
            ),
            "external_benchmarks": {
                "GW170817_GRB170817A": {
                    "source": "https://arxiv.org/abs/1710.05834",
                    "observable": "delta_T=c_T/c_gamma-1",
                    "interval": [-3.0e-15, 7.0e-16],
                    "scope": "associated low-redshift multimessenger path with bounded intrinsic emission delay",
                },
                "DESI_DR2_flat_wCDM": {
                    "source": "https://arxiv.org/abs/2503.14738v3",
                    "arxiv_version": "2503.14738v3",
                    "published_doi": "10.1103/tr6y-kpc6",
                    "location": "Table_V",
                    "scope": "published marginalized means and 68_percent sigmas in flat constant-w CDM with fiducial sum_mnu=0.06_eV",
                    "use": "separate approximate_1.96_sigma_summary_gates_not_a_likelihood",
                },
            },
            "nonimplications": (
                "No raw local c variation, parent speed, external origin, global "
                "transition, generic DSF-1 confirmation, or dark-energy identity "
                "follows from these conditional gates."
            ),
            "relative_cone": {
                "quantity": "delta_T=c_T/c_gamma-1",
                "null_control": {
                    **gw170817_cone_gate(tensor_photon_delta(1.0, 1.0)),
                    "rejected": False,
                },
                "positive_1e_minus_14": {
                    **gw170817_cone_gate(1.0e-14),
                    "rejected": True,
                },
                "negative_1e_minus_14": {
                    **gw170817_cone_gate(-1.0e-14),
                    "rejected": True,
                },
                "disformal_stationary_null_control": {
                    **disformal_photon_cone(
                        conformal_factor=1.0,
                        disformal_factor=1.0,
                        kinetic_x=0.0,
                    ),
                    "relative_cone_gate": gw170817_cone_gate(0.0),
                },
            },
            "boundary_scalings": {
                "constant_density": _literal_boundary_scaling_record(q=0.0, s=1.0),
                "constant_tension_surface": _literal_boundary_scaling_record(q=1.0, s=1.0),
                "inverse_area_density": _literal_boundary_scaling_record(q=2.0, s=1.0),
            },
        },
    }
