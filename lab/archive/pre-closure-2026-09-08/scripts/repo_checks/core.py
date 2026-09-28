"""Stable repository-check orchestration and CLI."""

from __future__ import annotations

from ._shared import *  # noqa: F401,F403
from .public import *  # noqa: F401,F403
from .calibration import *  # noqa: F401,F403
from .progression import *  # noqa: F401,F403
from .temporal_diagnostics import *  # noqa: F401,F403
from .temporal_selection import *  # noqa: F401,F403
from .catalog import *  # noqa: F401,F403



def check_results() -> list[str]:
    path = REPOSITORY / "results" / "core-identities.json"
    if not path.exists():
        return ["core result is absent; run scripts/reproduce_core.py"]
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return [f"cannot parse core result: {exc}"]

    failures: list[str] = []
    if record.get("project_version") != "0.11.0":
        failures.append("core result project_version is not 0.11.0")
    standard = record.get("standard_cosmology", {})
    ratio = standard.get("hubble_sphere_compactness")
    if not isinstance(ratio, (int, float)) or abs(ratio - 1.0) > 1e-13:
        failures.append("Hubble-sphere compactness identity is outside tolerance")
    conditional = record.get("conditional_horizon_saturation", {})
    if conditional.get("classification") != "conditional_inference_not_independent_prediction":
        failures.append("H-SAT result is missing its conditional-inference label")
    status = record.get("logical_status", {})
    if status.get("lambda_closure") != "algebraic_resubstitution_not_prediction":
        failures.append("Lambda closure is not explicitly labeled as resubstitution")

    controlled_path = REPOSITORY / "results" / "controlled-model.json"
    if not controlled_path.exists():
        failures.append(
            "controlled-model result is absent; run scripts/reproduce_controlled_model.py"
        )
        return failures
    try:
        controlled = json.loads(controlled_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse controlled-model result: {exc}")
        return failures
    if controlled.get("project_version") != "0.11.0":
        failures.append("controlled-model project_version is not 0.11.0")
    if controlled.get("schema_version") != 3:
        failures.append("controlled-model schema_version is not 3")
    covariant = controlled.get("covariant_core", {})
    if covariant.get("classification") != (
        "controlled_covariant_core_not_global_black_hole_transition"
    ):
        failures.append("CCT-1 result is missing its global-transition limitation")
    if covariant.get("matter_condition") != (
        "ordinary_matter_absent_vacuum_energy_absorbed_into_V0"
    ):
        failures.append("CCT-1 exact core is missing its vacuum/matter condition")
    information = controlled.get("information_map", {})
    if information.get("p_status") != "chosen_symmetric_benchmark_not_derived":
        failures.append("IM-1 p is not labeled as a chosen benchmark")
    information_audit = information.get("audit", {})
    if information_audit.get("is_recoverable_on_code_range") is not True:
        failures.append("IM-1 full-joint code-range recovery check is absent")
    perturbations = controlled.get("linear_perturbations", {})
    if perturbations.get("classification") != (
        "classical_core_transfer_not_primordial_spectrum"
    ):
        failures.append("perturbation result is missing its non-spectrum limitation")
    domain = controlled.get("domain_synthesis", {})
    if domain.get("classification") != (
        "controlled_effective_fluid_and_kottler_balance_benchmarks_"
        "not_dark_sector_or_recursive_transition"
    ):
        failures.append("DST-1 result is missing its dark-sector/transition limitation")
    exchange = domain.get("exchange_ledger", {})
    residual = exchange.get("total_continuity_residual")
    if not isinstance(residual, (int, float)) or abs(residual) > 1e-12:
        failures.append("DST-1 exchange ledger does not close within tolerance")
    kottler = domain.get("kottler_balance", {})
    if "not an accretion wall" not in str(kottler.get("interpretation", "")):
        failures.append("Kottler balance is missing its non-wall limitation")
    if kottler.get("stability") != "unstable_radial_separator":
        failures.append("Kottler balance is missing its unstable-separator label")
    discriminator = controlled.get("independent_discriminator", {})
    if discriminator.get("prediction") != "Omega_K<0_when_H_is_nonzero":
        failures.append("strict closed-branch curvature-sign prediction is absent")
    if discriminator.get("required_branch_postulate") != (
        "observable_child_remains_globally_closed_with_K=+1"
    ):
        failures.append("curvature sign is missing its added global branch postulate")
    gates = controlled.get("literal_model_gates", {})
    if gates.get("classification") != (
        "conditional_literal_subclass_rejection_gates_not_model_confirmation"
    ):
        failures.append("literal-model gates are absent or overclassified")
    benchmarks = gates.get("external_benchmarks", {})
    if benchmarks.get("GW170817_GRB170817A", {}).get("source") != (
        "https://arxiv.org/abs/1710.05834"
    ):
        failures.append("relative-cone gate is missing primary-source provenance")
    if benchmarks.get("DESI_DR2_flat_wCDM", {}).get("source") != (
        "https://arxiv.org/abs/2503.14738v3"
    ):
        failures.append("boundary summary gates are missing DESI DR2 provenance")
    if benchmarks.get("DESI_DR2_flat_wCDM", {}).get("published_doi") != (
        "10.1103/tr6y-kpc6"
    ):
        failures.append("boundary summary gates are missing the DESI publication DOI")
    cone = gates.get("relative_cone", {})
    if cone.get("null_control", {}).get("rejected") is not False:
        failures.append("relative-cone null control is not retained")
    if cone.get("positive_1e_minus_14", {}).get("rejected") is not True:
        failures.append("relative-cone positive rejection vector is absent")
    if cone.get("negative_1e_minus_14", {}).get("rejected") is not True:
        failures.append("relative-cone negative rejection vector is absent")
    boundary = gates.get("boundary_scalings", {})
    if boundary.get("constant_density", {}).get("effective_w") != -1.0:
        failures.append("constant-density boundary control is absent")
    surface_w = boundary.get("constant_tension_surface", {}).get("effective_w")
    if not isinstance(surface_w, (int, float)) or abs(surface_w - (-2.0 / 3.0)) > 1e-15:
        failures.append("surface-energy boundary scaling is absent")
    if boundary.get("inverse_area_density", {}).get("accelerates") is not False:
        failures.append("inverse-area non-acceleration rejection is absent")

    bao_path = REPOSITORY / "results" / "desi-dr2-bao-profile.json"
    try:
        bao = json.loads(bao_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse DESI DR2 BAO profile result: {exc}")
        return failures
    if bao.get("project_version") != "0.11.0" or bao.get("schema_version") != 1:
        failures.append("DESI DR2 BAO profile has stale version or schema")
    if bao.get("artifact") != "DESI_DR2_BAO_only_flat_constant_w_profile":
        failures.append("DESI DR2 BAO profile is missing its bounded artifact label")
    best = bao.get("best_fit", {})
    if not isinstance(best.get("w"), (int, float)) or not -0.93 < best["w"] < -0.90:
        failures.append("DESI DR2 BAO best-fit w is outside the audited range")
    fixed = bao.get("fixed_w_candidates", {})
    surface_delta = fixed.get("-0.6666666666666666", {}).get("delta_chi2_from_best")
    inverse_area_delta = fixed.get("-0.3333333333333333", {}).get("delta_chi2_from_best")
    if not isinstance(surface_delta, (int, float)) or not 10.0 < surface_delta < 12.5:
        failures.append("DESI BAO surface-scaling rejection result is absent")
    if not isinstance(inverse_area_delta, (int, float)) or not 65.0 < inverse_area_delta < 75.0:
        failures.append("DESI BAO inverse-area rejection result is absent")
    if "declared_search" not in bao or "preregistered_search" in bao:
        failures.append("DESI BAO search bounds are not labeled as declared post-publication choices")

    timing_path = REPOSITORY / "results" / "gw170817-timing.json"
    try:
        timing = json.loads(timing_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse GW170817 timing result: {exc}")
        return failures
    if timing.get("project_version") != "0.11.0" or timing.get("schema_version") != 1:
        failures.append("GW170817 timing artifact has stale version or schema")
    if timing.get("classification") != "published_timing_input_reconstruction_not_raw_reanalysis":
        failures.append("GW170817 timing artifact is missing its non-reanalysis boundary")
    timing_validation = timing.get("validation", {})
    delay = timing_validation.get("reconstructed_geocentric_delay_seconds")
    if not isinstance(delay, (int, float)) or abs(delay - 1.737774) > 1.0e-9:
        failures.append("GW170817 timing reconstruction is outside tolerance")
    if timing_validation.get("published_rounded_delta_v_over_v_em_interval") != [
        -3.0e-15,
        7.0e-16,
    ]:
        failures.append("GW170817 published rounded relative-speed interval is absent")

    gmf_path = REPOSITORY / "results" / "gmf-1-spherical-matching.json"
    try:
        gmf = json.loads(gmf_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse GMF-1 matching result: {exc}")
        return failures
    if gmf.get("project_version") != "0.11.0" or gmf.get("schema_version") != 1:
        failures.append("GMF-1 matching result has stale version or schema")
    if gmf.get("artifact") != "GMF-1_spherical_matching_and_flux_closure_audit":
        failures.append("GMF-1 matching result is missing its bounded artifact label")
    if gmf.get("classification") != "spherical_matching_obstruction_and_local_shell_diagnostic":
        failures.append("GMF-1 matching result is missing its bounded obstruction/diagnostic classification")
    if gmf.get("global_solution_constructed") is not False:
        failures.append("GMF-1 must explicitly retain global_solution_constructed=false")
    if gmf.get("thin_shell_distributional") is not True:
        failures.append("GMF-1 must label a thin Israel shell as distributional")
    if gmf.get("derived_exchange_law") is not False:
        failures.append("GMF-1 must explicitly retain derived_exchange_law=false")
    if gmf.get("child_topology_proven") is not False:
        failures.append("GMF-1 must explicitly retain child_topology_proven=false")
    mass_gate = gmf.get("darmois_mass_gate", {})
    if mass_gate.get("same_lambda_nonzero_parent_mass_direct_match_allowed") is not False:
        failures.append("GMF-1 direct same-Lambda nonzero-mass no-shell branch is not rejected")
    os_branch = gmf.get("oppenheimer_snyder", {}).get("endpoint_status", {})
    if os_branch.get("smooth_match") is not True or os_branch.get("singularity_present_at_a_zero") is not True:
        failures.append("GMF-1 Oppenheimer-Snyder smooth finite-mass reference must retain its singularity")
    os_record = gmf.get("oppenheimer_snyder", {})
    density_growth = os_record.get("density_growth_a_to_a_over_100")
    ricci_initial = os_record.get("ricci_scalar_at_a")
    ricci_final = os_record.get("ricci_scalar_at_a_over_100")
    if (
        not isinstance(density_growth, (int, float))
        or not isfinite(density_growth)
        or abs(density_growth - 1.0e6) > 1.0e-6
    ):
        failures.append("GMF-1 does not calculate the OS dust singular scaling")
    if (
        not isinstance(ricci_initial, (int, float))
        or not isinstance(ricci_final, (int, float))
        or not isfinite(ricci_initial)
        or not isfinite(ricci_final)
        or ricci_initial == 0.0
        or abs(ricci_final / ricci_initial - 1.0e6) > 1.0e-6
    ):
        failures.append("GMF-1 does not calculate the OS Ricci-scalar divergence")
    shell = gmf.get("prescribed_timelike_shell", {})
    if shell.get("orientation_gate", {}).get("accepted") is not True:
        failures.append("GMF-1 local shell control does not pass its oriented equation")
    static = shell.get("static_gate", {})
    if static.get("orientation_accepted_at_rest") is not True:
        failures.append("GMF-1 static gate does not retain the unsquared orientation check")
    if static.get("actually_static") is not False:
        failures.append("GMF-1 prescribed at-rest shell must not be promoted to equilibrium")
    closure = gmf.get("closure_status", {})
    if closure.get("bulk_Q_derived") is not False:
        failures.append("GMF-1 closure must retain bulk_Q_derived=false")

    ec_path = REPOSITORY / "results" / "einstein-cartan-collapse.json"
    try:
        ec = json.loads(ec_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse EC-1 collapse result: {exc}")
        return failures
    if ec.get("project_version") != "0.11.0" or ec.get("schema_version") != 1:
        failures.append("EC-1 collapse result has stale version or schema")
    if ec.get("model_id") != "EC-1":
        failures.append("EC-1 collapse result is missing its model identifier")
    if ec.get("artifact") != "einstein_cartan_homogeneous_spin_fluid_cap":
        failures.append("EC-1 collapse result is missing its bounded artifact label")
    if ec.get("classification") != (
        "reduced_action_motivated_homogeneous_effective_fluid_benchmark"
    ):
        failures.append("EC-1 collapse result is missing its reduced-benchmark classification")
    for key in (
        "global_solution_constructed",
        "smooth_static_vacuum_boundary_proven",
        "child_topology_proven",
        "derived_external_Q",
        "late_dark_energy_derived",
        "variable_c_derived",
        "radial_stability_proven",
        "microscopic_dirac_stability_or_eft_control_proven",
    ):
        if ec.get(key) is not False:
            failures.append(f"EC-1 must explicitly retain {key}=false")
    if ec.get("homogeneous_spin_fluid_closure_assumed") is not True:
        failures.append("EC-1 must explicitly label its homogeneous spin-fluid closure as assumed")
    if ec.get("single_interior_metric_causal_sequence") is not True:
        failures.append("EC-1 must explicitly retain its one-cap causal-sequence boundary")

    parameters = ec.get("parameters", {})
    for key, expected in (("A", 10.0), ("B", 1.0)):
        value = parameters.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value - expected) > 1e-12:
            failures.append(f"EC-1 benchmark parameter {key} is not the audited value")
    chi_boundary = parameters.get("chi_boundary")
    if not isinstance(chi_boundary, (int, float)) or not isfinite(chi_boundary) or abs(chi_boundary - 1.0471975511965976) > 1e-12:
        failures.append("EC-1 benchmark chi_boundary is not pi/3")
    turning = ec.get("turning_points", {})
    for key, expected in (("a_min", 0.31783724519578294), ("a_max", 3.146264369941972)):
        value = turning.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value - expected) > 1e-10:
            failures.append(f"EC-1 {key} is outside the audited finite-bounce value")
    bounce = ec.get("bounce", {})
    if not isinstance(bounce.get("acceleration"), (int, float)) or bounce["acceleration"] <= 0.0:
        failures.append("EC-1 bounce must retain positive acceleration")
    if not isinstance(bounce.get("friedmann_F"), (int, float)) or abs(bounce["friedmann_F"]) > 1e-9:
        failures.append("EC-1 bounce Friedmann constraint does not close")
    bounce_invariants = bounce.get("curvature_invariants", {})
    if not isinstance(bounce_invariants, dict) or set(bounce_invariants) != {
        "ricci_scalar",
        "ricci_tensor_squared",
        "kretschmann_scalar",
        "weyl_tensor_squared",
    } or not all(
        isinstance(value, (int, float)) and isfinite(value)
        for value in bounce_invariants.values()
    ):
        failures.append("EC-1 bounce must retain finite reduced curvature invariants")
    bounce_fluid = bounce.get("fluid", {})
    for key, expected in (("rho_plus_p", -75.61459384415677), ("rho_plus_3p", -229.20698880751104)):
        value = bounce_fluid.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value - expected) > 1e-9:
            failures.append(f"EC-1 bounce {key} is outside the audited defocusing value")
    sequence = ec.get("null_expansion_sequence", {})
    expected_regions = [
        "normal",
        "marginal",
        "trapped",
        "marginal",
        "normal",
        "marginal",
        "anti_trapped",
        "marginal",
        "normal",
    ]
    actual_regions = [entry.get("region") for entry in sequence.get("sequence", []) if isinstance(entry, dict)]
    if actual_regions != expected_regions:
        failures.append("EC-1 does not retain its audited one-cap causal sequence")
    if sequence.get("single_interior_metric_causal_sequence") is not True:
        failures.append("EC-1 causal sequence is not labeled as one interior metric")
    proper_time = ec.get("proper_time", {})
    half_cycle = proper_time.get("finer")
    if not isinstance(half_cycle, (int, float)) or not isfinite(half_cycle) or abs(half_cycle - 3.1974562013456302) > 1e-10:
        failures.append("EC-1 half-cycle proper time is outside the audited value")
    if not isinstance(proper_time.get("fine_to_finer_difference"), (int, float)) or proper_time["fine_to_finer_difference"] > 1e-10:
        failures.append("EC-1 half-cycle quadrature does not retain convergence")
    boundary = ec.get("boundary", {})
    for key in ("work_residual_collapse", "work_residual_expansion"):
        value = boundary.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value) > 1e-12:
            failures.append(f"EC-1 boundary {key} does not close")
    drift = boundary.get("static_parent_mass_drift", {})
    ratio = drift.get("bounce_to_initial_parent_mass_ratio")
    if not isinstance(ratio, (int, float)) or not isfinite(ratio) or abs(ratio - 0.10102051443368686) > 1e-10:
        failures.append("EC-1 static-boundary mass drift is outside the audited value")
    if drift.get("fixed_vacuum_parent_can_remain_smooth_without_flux_or_layer") is not False:
        failures.append("EC-1 must reject fixed-vacuum static boundary closure")

    preflight_path = REPOSITORY / "results" / "gmf-1b-preflight.json"
    try:
        preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse GMF-1B preflight result: {exc}")
        return failures
    if preflight.get("project_version") != "0.11.0" or preflight.get("schema_version") != 1:
        failures.append("GMF-1B preflight result has stale version or schema")
    if preflight.get("model_id") != "GMF-1B-PF1":
        failures.append("GMF-1B preflight result is missing its model identifier")
    if preflight.get("artifact") != "gmf1b_focusing_and_thermal_spin_fluid_preflight":
        failures.append("GMF-1B preflight result is missing its bounded artifact label")
    if preflight.get("classification") != (
        "conditional_material_law_and_null_focusing_rejection_gates"
    ):
        failures.append("GMF-1B preflight result is missing its restricted-gate classification")
    for key in (
        "global_solution_constructed",
        "general_recursive_horizons_no_go_proven",
        "canonical_scalar_regular_bridge_supported_under_stated_assumptions",
        "full_einstein_cartan_dirac_rejected",
        "ec1_background_algebra_rejected",
        "external_Q_derived",
        "late_dark_energy_derived",
        "variable_c_derived",
    ):
        if preflight.get(key) is not False:
            failures.append(f"GMF-1B preflight must explicitly retain {key}=false")
    focusing = preflight.get("canonical_scalar_null_focusing", {})
    if focusing.get("canonical_null_convergence") is not True:
        failures.append("GMF-1B preflight must retain canonical scalar null convergence")
    if focusing.get("same_generator_defocusing_to_zero_excluded_under_stated_assumptions") is not True:
        failures.append("GMF-1B preflight must retain the restricted canonical-scalar bridge rejection")
    if focusing.get("negative_expansion_can_rise_to_zero_before_caustic_or_endpoint") is not False:
        failures.append("GMF-1B preflight must not permit the rejected canonical-scalar bridge")
    regimes = preflight.get("thermal_ec_closure_regimes", {})
    expected_regimes = {
        "causal_example": "causal_barotrope",
        "degenerate_example": "degenerate_zero_sound_speed",
        "gradient_unstable_example": "gradient_unstable",
        "singular_example": "singular_enthalpy_and_density_derivative",
        "superluminal_example": "superluminal_relative_to_metric_cone",
    }
    for key, classification in expected_regimes.items():
        if regimes.get(key, {}).get("classification") != classification:
            failures.append(f"GMF-1B preflight thermal regime {key} is not classified correctly")
    bounce_gate = preflight.get("ec1_bounce_gate", {})
    for key, expected in (("A", 10.0), ("B", 1.0)):
        value = bounce_gate.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value - expected) > 1e-12:
            failures.append(f"GMF-1B preflight benchmark {key} is not the audited value")
    if bounce_gate.get("z_bounce_strictly_between_one_half_and_one") is not True:
        failures.append("GMF-1B preflight loses the EC-1 bounce thermal interval")
    if bounce_gate.get("healthy_causal_barotrope_at_bounce") is not False:
        failures.append("GMF-1B preflight must reject the naive thermal barotrope at bounce")
    if bounce_gate.get("thermal_classification_at_bounce") != (
        "superluminal_relative_to_metric_cone"
    ):
        failures.append("GMF-1B preflight loses the audited bounce characteristic class")
    bounce_cs2 = bounce_gate.get("cs2_at_bounce")
    if (
        not isinstance(bounce_cs2, (int, float))
        or not isfinite(bounce_cs2)
        or abs(bounce_cs2 - 2.3750044297870208) > 1e-12
    ):
        failures.append("GMF-1B preflight bounce sound speed is outside the audited value")
    if bounce_gate.get("continuous_low_z_to_bounce_branch_healthy_causal") is not False:
        failures.append("GMF-1B preflight must retain the unhealthy continuous-branch gate")
    if bounce_gate.get("naive_thermal_averaged_perfect_fluid_rejected_as_gmf1b_material_law") is not True:
        failures.append("GMF-1B preflight must retain its thermal-closure rejection")
    if bounce_gate.get("full_einstein_cartan_dirac_rejected") is not False:
        failures.append("GMF-1B preflight must not reject full Einstein-Cartan-Dirac")
    if bounce_gate.get("ec1_background_algebra_rejected") is not False:
        failures.append("GMF-1B preflight must not reject EC-1 background algebra")

    symmetry_path = REPOSITORY / "results" / "gmf-1b-ecd-symmetry.json"
    try:
        symmetry = json.loads(symmetry_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse GMF-1B-ECD-SYM1 result: {exc}")
        return failures
    if symmetry.get("project_version") != "0.11.0" or symmetry.get("schema_version") != 1:
        failures.append("GMF-1B-ECD-SYM1 result has stale version or schema")
    if symmetry.get("model_id") != "GMF-1B-ECD-SYM1":
        failures.append("GMF-1B-ECD-SYM1 result is missing its model identifier")
    if symmetry.get("artifact") != (
        "minimal_ecd_spherical_spinor_pair_axial_current_and_taper_preflight"
    ):
        failures.append("GMF-1B-ECD-SYM1 result is missing its bounded artifact label")
    if symmetry.get("classification") != (
        "local_action_locked_axial_current_and_metric_null_cone_gate"
    ):
        failures.append("GMF-1B-ECD-SYM1 result is missing its local preflight classification")
    conventions = symmetry.get("conventions", {})
    if conventions.get("metric_signature") != "-+++" or conventions.get("orientation_0123") != 1:
        failures.append("GMF-1B-ECD-SYM1 does not retain its frozen metric/orientation convention")
    for key, expected in (
        ("gravitational_coupling_kappa", 1.0),
        ("contact_action_coefficient_times_kappa", 3.0 / 16.0),
        ("hehl_datta_cubic_coefficient_times_kappa", 3.0 / 8.0),
    ):
        value = conventions.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value - expected) > 1e-15:
            failures.append(f"GMF-1B-ECD-SYM1 convention {key} is outside the audited value")
    fixture = symmetry.get("unit_normalization_algebraic_fixture", {})
    for key, expected in (
        ("F", 0.75),
        ("G", 0.5),
        ("A_hat_0", 13.0 / 8.0),
        ("A_hat_r", 5.0 / 8.0),
        ("axial_squared", -9.0 / 4.0),
        ("contact_scalar", -27.0 / 64.0),
    ):
        value = fixture.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value - expected) > 1e-14:
            failures.append(f"GMF-1B-ECD-SYM1 algebraic fixture {key} is outside tolerance")
    finster = symmetry.get("finster_full_dirac_singlet_control", {})
    finster_residual = finster.get("maximum_direct_axial_residual")
    if finster.get("pauli_doublet_cancellation_certified_to_binary64_roundoff") is not True:
        failures.append("GMF-1B-ECD-SYM1 lacks the direct Finster Pauli-doublet certificate")
    if (
        not isinstance(finster_residual, (int, float))
        or not isfinite(finster_residual)
        or finster_residual > 1e-12
    ):
        failures.append("GMF-1B-ECD-SYM1 Finster axial residual is outside tolerance")
    if finster.get("gamma5_definition") != "i*gamma0*gamma1*gamma2*gamma3":
        failures.append("GMF-1B-ECD-SYM1 Finster control has the wrong gamma-five convention")
    if finster.get("single_spinor_1_axial_current_nonzero") is not True:
        failures.append("GMF-1B-ECD-SYM1 lacks a nonzero single-spinor Finster regression")
    single_spinor_residual = finster.get("single_spinor_1_max_imaginary_residual")
    if (
        not isinstance(single_spinor_residual, (int, float))
        or not isfinite(single_spinor_residual)
        or single_spinor_residual > 1e-12
    ):
        failures.append("GMF-1B-ECD-SYM1 Finster single-spinor bilinear is outside tolerance")
    if finster.get("axial_current_nonzero") is not False:
        failures.append("GMF-1B-ECD-SYM1 Finster control must retain axial cancellation")
    if finster.get("algebraic_torsion_source_nonzero") is not False:
        failures.append("GMF-1B-ECD-SYM1 Finster control must have zero algebraic torsion source")
    if finster.get("axial_contact_invariant_nonzero") is not False:
        failures.append("GMF-1B-ECD-SYM1 Finster control must have zero contact invariant")
    ventrella = symmetry.get("ventrella_left_chiral_pair", {})
    if ventrella.get("axial_current_survives_spherical_pair_sum") is not True:
        failures.append("GMF-1B-ECD-SYM1 loses the Ventrella chiral-pair axial source")
    if ventrella.get("algebraic_torsion_source_survives_spherical_pair_sum") is not True:
        failures.append("GMF-1B-ECD-SYM1 loses the Ventrella chiral-pair torsion source")
    if ventrella.get("axial_contact_invariant_nonzero_when_both_amplitudes_populated") is not True:
        failures.append("GMF-1B-ECD-SYM1 loses the populated-pair contact invariant")
    one_amplitude_controls = ventrella.get("one_amplitude_zero_controls", {})
    for label in ("F_zero", "G_zero"):
        control = one_amplitude_controls.get(label, {})
        if control.get("axial_current_nonzero") is not True:
            failures.append(f"GMF-1B-ECD-SYM1 {label} must retain its nonzero null axial current")
        if control.get("algebraic_torsion_source_nonzero") is not True:
            failures.append(f"GMF-1B-ECD-SYM1 {label} must retain its algebraic torsion source")
        if control.get("axial_contact_invariant_nonzero") is not False:
            failures.append(f"GMF-1B-ECD-SYM1 {label} must have zero axial contact invariant")
        if control.get("contact_action_density_nonzero") is not False:
            failures.append(f"GMF-1B-ECD-SYM1 {label} must have zero contact action density")
    gamma_certificate = ventrella.get("gamma_matrix_certificate", {})
    gamma_residual = gamma_certificate.get("maximum_formula_residual")
    if gamma_certificate.get("formula_certified_to_binary64_roundoff") is not True:
        failures.append("GMF-1B-ECD-SYM1 lacks its direct Ventrella gamma certificate")
    if (
        not isinstance(gamma_residual, (int, float))
        or not isfinite(gamma_residual)
        or gamma_residual > 1e-12
    ):
        failures.append("GMF-1B-ECD-SYM1 Ventrella gamma residual is outside tolerance")
    cone_gate = symmetry.get("principal_cone_gate", {})
    if cone_gate.get("outgoing_null", {}).get("metric_null_characteristic") is not True:
        failures.append("GMF-1B-ECD-SYM1 outgoing Dirac characteristic is not metric-null")
    if cone_gate.get("ingoing_null", {}).get("metric_null_characteristic") is not True:
        failures.append("GMF-1B-ECD-SYM1 ingoing Dirac characteristic is not metric-null")
    if cone_gate.get("non_null_control", {}).get("metric_null_characteristic") is not False:
        failures.append("GMF-1B-ECD-SYM1 non-null cone control is not retained")
    if cone_gate.get("outgoing_null", {}).get("algebraic_torsion_changes_metric_cone") is not False:
        failures.append("GMF-1B-ECD-SYM1 must retain the lower-order contact/cone boundary")
    annular = symmetry.get("annular_packet_gate", {})
    for key in (
        "c_infinity_compact_support",
        "centre_spinor_zero_buffer",
        "outer_spinor_zero_buffer",
        "sampled_axial_bilinears_and_contact_scalar_finite",
        "nonzero_interior_axial_contact_invariant_present",
        "spinor_boundary_axial_current_zero",
    ):
        if annular.get(key) is not True:
            failures.append(f"GMF-1B-ECD-SYM1 annular packet gate loses {key}=true")
    if annular.get("sample_count") != 257:
        failures.append("GMF-1B-ECD-SYM1 annular packet sample count is not 257")
    for key in (
        "full_ecd_contact_stress_verified",
        "full_ecd_source_so3_equivariance_verified",
        "hehl_datta_angular_closure_verified",
        "curvature_invariants_verified",
        "constraint_solution_constructed",
        "constraint_compatible_vacuum_region_verified",
        "horizon_regular_evolution_constructed",
        "global_solution_constructed",
        "child_topology_proven",
        "external_Q_derived",
        "late_dark_energy_derived",
        "variable_c_derived",
        "quantum_ensemble_equivalence_proven",
    ):
        if symmetry.get(key) is not False:
            failures.append(f"GMF-1B-ECD-SYM1 must explicitly retain {key}=false")

    identity_path = REPOSITORY / "results" / "gmf-1b-ecd-identity.json"
    try:
        identity_record = json.loads(identity_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse GMF-1B-ECD-INT1 result: {exc}")
        return failures
    if (
        identity_record.get("project_version") != "0.11.0"
        or identity_record.get("schema_version") != 2
    ):
        failures.append("GMF-1B-ECD-INT1 result has stale version or schema")
    if identity_record.get("model_id") != "GMF-1B-ECD-INT1":
        failures.append("GMF-1B-ECD-INT1 result is missing its model identifier")
    if identity_record.get("artifact") != "ecd_interaction_action_gamma_identity_gate":
        failures.append("GMF-1B-ECD-INT1 result is missing its bounded artifact label")
    if identity_record.get("classification") != (
        "bounded_ecd_contact_identity_and_tested_local_reconstruction_preflight"
    ):
        failures.append("GMF-1B-ECD-INT1 result is missing its restricted classification")

    exact_fixture = identity_record.get("exact_fixture", {})
    signature_bridge = exact_fixture.get("signature_bridge", {})
    if signature_bridge.get("canonical_signature") != "+---" or signature_bridge.get(
        "vc_signature"
    ) != "-+++":
        failures.append("GMF-1B-ECD-INT1 does not retain its explicit signature bridge")
    if signature_bridge.get("gamma_conversion") != "Gamma^a=-i*gamma_VC^a":
        failures.append("GMF-1B-ECD-INT1 does not retain its gamma conversion")
    if signature_bridge.get("canonical_action_defined_before_translation") is not True:
        failures.append("GMF-1B-ECD-INT1 lacks the canonical action provenance gate")
    if signature_bridge.get("signature_bridge_closed") is not True:
        failures.append("GMF-1B-ECD-INT1 signature bridge is not closed")
    if signature_bridge.get("line_element_physically_changed") is not False:
        failures.append("GMF-1B-ECD-INT1 must not change the physical line element")
    signature_residual = signature_bridge.get("maximum_signature_bridge_residual")
    if (
        not isinstance(signature_residual, (int, float))
        or not isfinite(signature_residual)
        or signature_residual < 0.0
        or signature_residual > 1e-12
    ):
        failures.append("GMF-1B-ECD-INT1 signature bridge residual is outside tolerance")
    contact_density = exact_fixture.get("reduced_contact_density")
    if (
        not isinstance(contact_density, (int, float))
        or not isfinite(contact_density)
        or abs(contact_density - (-27.0 / 64.0)) > 1e-12
    ):
        failures.append("GMF-1B-ECD-INT1 contact density is outside the audited fixture")
    action_sources = exact_fixture.get("action_sources", {})
    if not isinstance(action_sources.get("common"), (int, float)) or not isfinite(
        action_sources["common"]
    ):
        failures.append("GMF-1B-ECD-INT1 lacks finite reduced action sources")
    for key, expected in (("source_F", -9.0 / 16.0), ("source_G", -27.0 / 32.0)):
        component = action_sources.get(key, {})
        real = component.get("real") if isinstance(component, dict) else None
        imag = component.get("imag") if isinstance(component, dict) else None
        if (
            isinstance(real, bool)
            or not isinstance(real, (int, float))
            or not isfinite(real)
            or abs(real - expected) > 1e-12
            or isinstance(imag, bool)
            or not isinstance(imag, (int, float))
            or not isfinite(imag)
            or abs(imag) > 1e-12
        ):
            failures.append(f"GMF-1B-ECD-INT1 {key} is outside the audited fixture")

    weak_identity = exact_fixture.get("identity", {})
    if weak_identity.get("weak_form_identity_closed") is not True:
        failures.append("GMF-1B-ECD-INT1 weak action identity is not closed")
    if weak_identity.get("all_tested_angle_reconstructions_closed") is not True:
        failures.append("GMF-1B-ECD-INT1 tested-angle local reconstruction is not retained")
    for key in (
        "maximum_weak_identity_residual",
        "maximum_tested_angle_reconstruction_residual",
        "maximum_coefficient_residual",
    ):
        residual = weak_identity.get(key)
        if (
            not isinstance(residual, (int, float))
            or not isfinite(residual)
            or residual < 0.0
            or residual > 1e-12
        ):
            failures.append(f"GMF-1B-ECD-INT1 {key} is outside tolerance")

    null_contact = identity_record.get("null_contact_gate", {})
    if null_contact.get("validated_metric_null") is not True:
        failures.append("GMF-1B-ECD-INT1 null-contact vector is not metric-null")
    for key in ("metric_norm", "contact_null_contraction"):
        value = null_contact.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value) > 1e-12:
            failures.append(f"GMF-1B-ECD-INT1 {key} is outside tolerance")
    stress = identity_record.get("contact_tetrad_stress")
    expected_diagonal = (27.0 / 64.0, -27.0 / 64.0, -27.0 / 64.0, -27.0 / 64.0)
    if not isinstance(stress, list) or len(stress) != 4 or any(
        not isinstance(row, list) or len(row) != 4 for row in stress
    ):
        failures.append("GMF-1B-ECD-INT1 contact stress fixture is not a 4 by 4 matrix")
    else:
        for row in range(4):
            for column in range(4):
                value = stress[row][column]
                expected = expected_diagonal[row] if row == column else 0.0
                if (
                    isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or not isfinite(value)
                    or abs(value - expected) > 1e-12
                    or abs(value - stress[column][row]) > 1e-12
                ):
                    failures.append(
                        "GMF-1B-ECD-INT1 contact stress fixture is malformed"
                    )
                    break
            else:
                continue
            break
    if null_contact.get("contact_sector_tetrad_variation_complete") is not True:
        failures.append("GMF-1B-ECD-INT1 contact tetrad variation is not retained")
    for key in ("full_ecd_stress_derived", "contact_term_alone_derives_null_defocusing"):
        if null_contact.get(key) is not False:
            failures.append(f"GMF-1B-ECD-INT1 must explicitly retain {key}=false")

    one_amplitude = identity_record.get("one_amplitude_semantics", {}).get("axial", {})
    if one_amplitude.get("axial_current_nonzero") is not True:
        failures.append("GMF-1B-ECD-INT1 loses the one-amplitude axial-current control")
    if one_amplitude.get("algebraic_torsion_source_nonzero") is not True:
        failures.append("GMF-1B-ECD-INT1 loses the one-amplitude torsion control")
    for key in ("axial_contact_invariant_nonzero", "contact_action_density_nonzero"):
        if one_amplitude.get(key) is not False:
            failures.append(f"GMF-1B-ECD-INT1 one-amplitude control must retain {key}=false")

    for key in (
        "complete_effective_stress_tensor_derived",
        "full_ecd_effective_stress_verified",
        "full_ecd_tetrad_dependence_varied",
        "torsion_free_dirac_stress_derived",
        "full_ecd_source_so3_equivariance_verified",
        "full_ecd_stress_verified",
        "constraint_solution_constructed",
        "constraint_compatible_vacuum_region_verified",
        "null_constraint_solution_constructed",
        "global_solution_constructed",
        "external_Q_derived",
        "late_dark_energy_derived",
        "variable_c_derived",
    ):
        if identity_record.get("nonclaims", {}).get(key) is not False:
            failures.append(f"GMF-1B-ECD-INT1 must explicitly retain {key}=false")

    torsion_path = REPOSITORY / "results" / "gmf-1b-ecd-torsion.json"
    try:
        torsion_record = json.loads(torsion_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse GMF-1B-ECD-TOR1 result: {exc}")
        return failures
    if (
        torsion_record.get("project_version") != "0.11.0"
        or torsion_record.get("schema_version") != 2
    ):
        failures.append("GMF-1B-ECD-TOR1 result has stale version or schema")
    if torsion_record.get("model_id") != "GMF-1B-ECD-TOR1":
        failures.append("GMF-1B-ECD-TOR1 result is missing its model identifier")
    if torsion_record.get("artifact") != "ecd_spin_tensor_torsion_and_effective_contact_crosscheck":
        failures.append("GMF-1B-ECD-TOR1 result is missing its bounded artifact label")
    if torsion_record.get("classification") != (
        "spin_tensor_torsion_invariant_and_reduced_contact_tetrad_stress_gate"
    ):
        failures.append("GMF-1B-ECD-TOR1 result is missing its restricted classification")

    axial = torsion_record.get("axial_current", {})
    for key, expected in (
        ("A_hat_0", 13.0 / 8.0),
        ("A_hat_r", 5.0 / 8.0),
        ("axial_squared", -9.0 / 4.0),
    ):
        value = axial.get(key)
        if (
            not isinstance(value, (int, float))
            or not isfinite(value)
            or abs(value - expected) > 1e-12
        ):
            failures.append(f"GMF-1B-ECD-TOR1 axial current {key} is outside tolerance")

    spin = torsion_record.get("spin_tensor", {})
    if spin.get("hodge_dual_verified") is not True:
        failures.append("GMF-1B-ECD-TOR1 spin tensor Hodge dual is not verified")
    if spin.get("totally_antisymmetric") is not True:
        failures.append("GMF-1B-ECD-TOR1 spin tensor is not labeled totally antisymmetric")
    components = spin.get("components", {})
    for key, expected in (("S_012", 5.0 / 16.0), ("S_123", 13.0 / 16.0)):
        value = components.get(key)
        if (
            not isinstance(value, (int, float))
            or not isfinite(value)
            or abs(value - expected) > 1e-12
        ):
            failures.append(f"GMF-1B-ECD-TOR1 spin tensor component {key} is outside tolerance")

    torsion = torsion_record.get("torsion", {})
    if torsion.get("relation_verified") is not True:
        failures.append("GMF-1B-ECD-TOR1 torsion-squared relation is not verified")
    if (
        not isinstance(torsion.get("torsion_squared"), (int, float))
        or not isfinite(torsion["torsion_squared"])
        or abs(torsion["torsion_squared"] - 27.0 / 8.0) > 1e-12
    ):
        failures.append("GMF-1B-ECD-TOR1 torsion squared is outside tolerance")

    contact = torsion_record.get("contact_derivation", {})
    if contact.get("coefficient_verified") is not True:
        failures.append("GMF-1B-ECD-TOR1 contact coefficient is not verified")
    if (
        not isinstance(contact.get("full_effective_contact_lagrangian"), (int, float))
        or not isfinite(contact["full_effective_contact_lagrangian"])
        or abs(contact["full_effective_contact_lagrangian"] - (-27.0 / 64.0)) > 1e-12
    ):
        failures.append("GMF-1B-ECD-TOR1 contact Lagrangian is outside tolerance")
    if contact.get("derived_after_eliminating_connection_everywhere") is not True:
        failures.append("GMF-1B-ECD-TOR1 does not retain the full connection elimination")
    for key in (
        "einstein_hilbert_piece_alone_claimed_to_equal_full_contact",
        "separate_connection_pieces_added_to_reduced_action",
    ):
        if contact.get(key) is not False:
            failures.append(f"GMF-1B-ECD-TOR1 must explicitly retain {key}=false")

    stress = torsion_record.get("contact_tetrad_stress", {})
    if stress.get("interpretation") != "positive_cosmological_constant_w_minus_1":
        failures.append("GMF-1B-ECD-TOR1 stress is missing its w=-1 interpretation")
    if stress.get("rho_plus_p_negative") is not False:
        failures.append("GMF-1B-ECD-TOR1 must retain rho_plus_p_negative=false")
    if stress.get("contact_sector_tetrad_variation_complete") is not True:
        failures.append("GMF-1B-ECD-TOR1 contact tetrad variation is not retained")
    for key, expected in (
        ("rho", 27.0 / 64.0),
        ("p", -27.0 / 64.0),
        ("w", -1.0),
        ("rho_plus_p", 0.0),
        ("rho_plus_3p", -27.0 / 32.0),
    ):
        value = stress.get(key)
        if (
            not isinstance(value, (int, float))
            or not isfinite(value)
            or abs(value - expected) > 1e-12
        ):
            failures.append(f"GMF-1B-ECD-TOR1 stress {key} is outside tolerance")

    for key in (
        "full_ecd_tetrad_dependence_varied",
        "full_ecd_effective_stress_verified",
        "complete_effective_stress_tensor_derived",
        "homogeneous_spin_fluid_closure_derived",
        "metric_null_defocusing_derived",
        "rho_plus_p_negative_derived",
        "dark_energy_value_fixed",
        "bounce_constructed",
        "constraint_solution_constructed",
        "regular_centre_constructed",
        "global_solution_constructed",
        "child_topology_proven",
        "external_Q_derived",
        "late_dark_energy_derived",
        "variable_c_derived",
    ):
        if torsion_record.get("nonclaims", {}).get(key) is not False:
            failures.append(f"GMF-1B-ECD-TOR1 must explicitly retain {key}=false")

    kinetic_path = REPOSITORY / "results" / "gmf-1b-ecd-kinetic.json"
    try:
        kinetic_record = json.loads(kinetic_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse GMF-1B-ECD-KIN1 result: {exc}")
        return failures
    if (
        kinetic_record.get("project_version") != "0.11.0"
        or kinetic_record.get("schema_version") != 2
    ):
        failures.append("GMF-1B-ECD-KIN1 result has stale version or schema")
    if kinetic_record.get("model_id") != "GMF-1B-ECD-KIN1":
        failures.append("GMF-1B-ECD-KIN1 result is missing its model identifier")
    if kinetic_record.get("artifact") != "ecd_contact_tetrad_variation_and_free_dirac_control":
        failures.append("GMF-1B-ECD-KIN1 result is missing its bounded artifact label")
    if kinetic_record.get("classification") != (
        "corrected_contact_stress_and_bounded_torsion_free_dirac_preflight"
    ):
        failures.append("GMF-1B-ECD-KIN1 result is missing its restricted classification")

    conventions = kinetic_record.get("conventions", {})
    if conventions.get("canonical_signature") != "+---" or conventions.get("vc_signature") != "-+++":
        failures.append("GMF-1B-ECD-KIN1 does not retain its signature bridge labels")
    if conventions.get("gamma_conversion") != "Gamma^a=-i*gamma_VC^a":
        failures.append("GMF-1B-ECD-KIN1 does not retain its gamma conversion")
    if conventions.get("orientation_0123") != 1:
        failures.append("GMF-1B-ECD-KIN1 does not retain its orientation convention")

    anticommutator = kinetic_record.get("anticommutator_identity", {})
    if anticommutator.get("verified_over_all_64_index_triples") is not True:
        failures.append("GMF-1B-ECD-KIN1 anticommutator identity is not verified")
    if anticommutator.get("representation_independent") is not True:
        failures.append("GMF-1B-ECD-KIN1 identity is not labeled representation-independent")
    residual = anticommutator.get("maximum_residual")
    if not isinstance(residual, (int, float)) or not isfinite(residual) or abs(residual) > 1e-12:
        failures.append("GMF-1B-ECD-KIN1 anticommutator residual is outside tolerance")

    action = kinetic_record.get("effective_action_bookkeeping", {})
    for key, expected in (
        ("canonical_J_squared", 9.0 / 4.0),
        ("vc_A_squared", -9.0 / 4.0),
        ("canonical_full_contact_coefficient_over_kappa_J2", -3.0 / 16.0),
        ("vc_full_contact_coefficient_over_kappa_A2", 3.0 / 16.0),
        ("canonical_full_contact_lagrangian", -27.0 / 64.0),
        ("vc_full_contact_lagrangian", -27.0 / 64.0),
        ("hehl_datta_equation_coefficient_over_kappa", 3.0 / 8.0),
    ):
        value = action.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value - expected) > 1e-12:
            failures.append(f"GMF-1B-ECD-KIN1 action {key} is outside tolerance")
    for key in (
        "signature_translation_verified",
        "field_variation_factor_two_verified",
        "connection_eliminated_everywhere_before_reduced_action_used",
        "combined_onshell_coefficient_sign_resolved",
    ):
        if action.get(key) is not True:
            failures.append(f"GMF-1B-ECD-KIN1 action loses {key}=true")
    if action.get("separate_connection_pieces_added_to_reduced_action") is not False:
        failures.append("GMF-1B-ECD-KIN1 must not double count connection pieces")

    cancellation = kinetic_record.get("tetrad_current_invariant_cancellation", {})
    if cancellation.get("coframe_directions_checked") != 16:
        failures.append("GMF-1B-ECD-KIN1 does not check all 16 coframe directions")
    for key in ("baseline_invariant_residual", "maximum_total_variation_residual"):
        value = cancellation.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value) > 1e-12:
            failures.append(f"GMF-1B-ECD-KIN1 tetrad cancellation {key} is outside tolerance")
    fixed_partial = cancellation.get("maximum_fixed_coordinate_current_metric_partial")
    if not isinstance(fixed_partial, (int, float)) or not isfinite(fixed_partial) or fixed_partial <= 1e-3:
        failures.append("GMF-1B-ECD-KIN1 fixed-current partial control is not nonzero")
    for key in (
        "fixed_coordinate_current_partial_is_nonzero",
        "current_response_cancels_metric_partial",
        "internal_current_invariant_is_tetrad_independent",
    ):
        if cancellation.get(key) is not True:
            failures.append(f"GMF-1B-ECD-KIN1 tetrad cancellation loses {key}=true")
    if cancellation.get("fixed_coordinate_current_anisotropic_stress_is_physical") is not False:
        failures.append("GMF-1B-ECD-KIN1 must retract the fixed-current stress")

    contact = kinetic_record.get("contact_stress_null_energy", {})
    for key, expected in (
        ("rho", 27.0 / 64.0),
        ("p", -27.0 / 64.0),
        ("w", -1.0),
        ("rho_plus_p", 0.0),
        ("maximum_null_contraction_residual", 0.0),
    ):
        value = contact.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value - expected) > 1e-12:
            failures.append(f"GMF-1B-ECD-KIN1 contact stress {key} is outside tolerance")
    if contact.get("null_vectors_checked") != 2:
        failures.append("GMF-1B-ECD-KIN1 contact gate does not retain two null probes")
    for key in ("all_contact_null_contractions_zero", "contact_sector_tetrad_variation_complete"):
        if contact.get(key) is not True:
            failures.append(f"GMF-1B-ECD-KIN1 contact gate loses {key}=true")
    for key in ("contact_null_energy_condition_violated", "full_ecd_effective_stress_verified"):
        if contact.get(key) is not False:
            failures.append(f"GMF-1B-ECD-KIN1 contact gate must retain {key}=false")

    plane_wave = kinetic_record.get("free_dirac_plane_wave_control", {})
    for key in (
        "momentum_norm",
        "massless_dirac_residual",
        "current_equals_two_momentum_residual",
        "null_probe_norm_residual",
    ):
        value = plane_wave.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value) > 1e-12:
            failures.append(f"GMF-1B-ECD-KIN1 plane-wave {key} is outside tolerance")
    contractions = plane_wave.get("null_contractions", {})
    for key, expected in (("parallel", 0.0), ("opposite", 8.0)):
        value = contractions.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value - expected) > 1e-12:
            failures.append(f"GMF-1B-ECD-KIN1 plane-wave {key} contraction is outside tolerance")
    for key in (
        "positive_frequency_plane_wave_nec_nonnegative",
        "positive_frequency_plane_wave_control_only",
    ):
        if plane_wave.get(key) is not True:
            failures.append(f"GMF-1B-ECD-KIN1 plane-wave loses {key}=true")
    if plane_wave.get("general_classical_dirac_nec_proven") is not False:
        failures.append("GMF-1B-ECD-KIN1 must not claim a general Dirac NEC theorem")

    for key in (
        "full_ecd_tetrad_dependence_varied",
        "torsion_free_dirac_stress_fully_derived",
        "full_ecd_effective_stress_verified",
        "general_classical_dirac_nec_proven",
        "rho_plus_p_negative_derived",
        "homogeneous_bounce_rho_plus_p_negative_derived",
        "metric_null_defocusing_derived",
        "bounce_constructed",
        "global_solution_constructed",
    ):
        if kinetic_record.get("nonclaims", {}).get(key) is not False:
            failures.append(f"GMF-1B-ECD-KIN1 must explicitly retain {key}=false")

    fgc_path = REPOSITORY / "results" / "fgc-1-action-gate.json"
    try:
        fgc_record = json.loads(fgc_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse FGC-1-ACT1 result: {exc}")
        return failures
    if (
        fgc_record.get("project_version") != "0.11.0"
        or fgc_record.get("schema_version") != 1
    ):
        failures.append("FGC-1-ACT1 result has stale version or schema")
    if fgc_record.get("artifact_id") != "FGC-1-ACT1":
        failures.append("FGC-1-ACT1 result is missing its artifact identifier")
    if fgc_record.get("classification") != (
        "action_and_background_algebra_not_principal_symbol_or_collapse_result"
    ):
        failures.append("FGC-1-ACT1 result is missing its restricted classification")
    if fgc_record.get("source_config") != "configs/fgc/fgc-1-action-gate.toml":
        failures.append("FGC-1-ACT1 result does not name the frozen source config")
    fgc_config = REPOSITORY / "configs" / "fgc" / "fgc-1-action-gate.toml"
    try:
        config_digest = sha256(fgc_config.read_bytes()).hexdigest()
    except OSError as exc:
        failures.append(f"cannot hash FGC-1-ACT1 config: {exc}")
    else:
        if fgc_record.get("source_config_sha256") != config_digest:
            failures.append("FGC-1-ACT1 config hash does not match its source")
    implementation_sources = fgc_record.get("implementation_sources_sha256", {})
    for relative in (
        "src/recursive_horizons/fgc/action.py",
        "scripts/reproduce_fgc_action.py",
    ):
        try:
            source_digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-ACT1 source {relative}: {exc}")
            continue
        if implementation_sources.get(relative) != source_digest:
            failures.append(f"FGC-1-ACT1 source hash is stale for {relative}")

    fixtures = fgc_record.get("fixtures", [])
    if not isinstance(fixtures, list) or [
        fixture.get("model_id") if isinstance(fixture, dict) else None
        for fixture in fixtures
    ] != ["GR-0", "SGB-L", "FGC-QR"]:
        failures.append("FGC-1-ACT1 does not retain the frozen branch order")
    else:
        for fixture in fixtures:
            action_gate = fixture.get("action_gate", {})
            for key in (
                "candidate_covariant_action_declared",
                "scalar_equations_derived",
                "metric_equation_variationally_defined",
                "separate_matter_and_regulator_fields",
            ):
                if action_gate.get(key) is not True:
                    failures.append(
                        f"FGC-1-ACT1 {fixture['model_id']} loses {key}=true"
                    )
            weak = fixture.get("weak_field_kinetic_reference", {})
            if weak.get("full_coupled_quadratic_action_diagonalized") is not False:
                failures.append(
                    f"FGC-1-ACT1 {fixture['model_id']} must keep the full kinetic gate false"
                )
            if weak.get("nonlinear_ghost_freedom_proven") is not False:
                failures.append(
                    f"FGC-1-ACT1 {fixture['model_id']} must keep ghost freedom unproved"
                )
            low_gradient = fixture.get("low_gradient_reference_gate", {})
            if low_gradient.get("coupled_metric_chi_background_verified") is not False:
                failures.append(
                    f"FGC-1-ACT1 {fixture['model_id']} must keep the coupled background unverified"
                )
            if any(value is not False for value in fixture.get("nonclaims", {}).values()):
                failures.append(
                    f"FGC-1-ACT1 {fixture['model_id']} has a promoted nonclaim"
                )
        qr_activation = fixtures[2].get("schwarzschild_linear_activation_control", {})
        activation_mass = qr_activation.get("linear_effective_mass_squared_at_activation")
        if (
            isinstance(activation_mass, bool)
            or not isinstance(activation_mass, (int, float))
            or not isfinite(activation_mass)
            or abs(activation_mass) > 1e-12
        ):
            failures.append("FGC-1-ACT1 activation mass residual is outside tolerance")
        if qr_activation.get("interpretation") != (
            "linear_instability_scale_not_transition_or_eft_cutoff"
        ):
            failures.append("FGC-1-ACT1 promotes its linear activation scale")

    algebraic_checks = fgc_record.get("verified_algebraic_checks", {})
    for key in (
        "all_required_nonclaims_false",
        "fgc_qr_inside_sample_is_linearly_tachyonic",
        "fgc_qr_zero_regulator_solves_declared_scalar_equation",
        "frozen_branch_order",
        "gr0_exact_action_at_zero_regulator",
        "schwarzschild_activation_zero_within_relative_tolerance_1e_12",
        "sgb_l_is_curvature_sourced_at_zero_regulator",
    ):
        if algebraic_checks.get(key) is not True:
            failures.append(f"FGC-1-ACT1 loses algebraic check {key}=true")
    gate_status = fgc_record.get("gate_status", {})
    if gate_status.get("declared_action_and_scalar_algebra_certificate_passed") is not True:
        failures.append("FGC-1-ACT1 algebra certificate is not marked passed")
    for key in ("full_fgc1_health_gate_passed", "publication_local_defocusing_gate_passed"):
        if gate_status.get(key) is not False:
            failures.append(f"FGC-1-ACT1 must retain {key}=false")
    required_fgc_nonclaims = (
        "full_covariant_metric_variation_expanded_and_cross_checked",
        "coupled_metric_chi_background_verified",
        "spherical_principal_symbol_derived",
        "nonlinear_strong_hyperbolicity_proven",
        "metric_null_defocusing_derived",
        "finite_collapse_solution_constructed",
        "singularity_resolution_proven",
        "child_spacetime_constructed",
        "dark_matter_derived",
        "dark_energy_derived",
        "variable_speed_of_light_derived",
        "particle_spectrum_derived",
        "theory_of_everything_derived",
    )
    for key in required_fgc_nonclaims:
        if fgc_record.get("nonclaims", {}).get(key) is not False:
            failures.append(f"FGC-1-ACT1 must explicitly retain {key}=false")

    var1_path = REPOSITORY / "results" / "fgc-1-metric-variation.json"
    try:
        var1_record = json.loads(var1_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse FGC-1-VAR1 result: {exc}")
        return failures
    if (
        var1_record.get("project_version") != "0.11.0"
        or var1_record.get("schema_version") != 1
    ):
        failures.append("FGC-1-VAR1 result has stale version or schema")
    if var1_record.get("artifact_id") != "FGC-1-VAR1":
        failures.append("FGC-1-VAR1 result is missing its artifact identifier")
    if var1_record.get("artifact") != (
        "fgc1_covariant_metric_variation_and_exact_identity_certificate"
    ):
        failures.append("FGC-1-VAR1 result is missing its artifact label")
    if var1_record.get("classification") != (
        "covariant_metric_equation_derivation_and_crosschecks_not_principal_symbol_or_collapse"
    ):
        failures.append("FGC-1-VAR1 result is missing its restricted classification")
    if var1_record.get("source_config") != (
        "configs/fgc/fgc-1-metric-variation.toml"
    ):
        failures.append("FGC-1-VAR1 does not name its frozen source config")

    expected_var1_config_sources = (
        "configs/fgc/fgc-1-metric-variation.toml",
        "configs/fgc/fgc-1-action-gate.toml",
    )
    var1_config_hashes = var1_record.get("source_configs_sha256", {})
    if set(var1_config_hashes) != set(expected_var1_config_sources):
        failures.append("FGC-1-VAR1 source-config hash ledger has unexpected keys")
    for relative in expected_var1_config_sources:
        try:
            source_digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-VAR1 config {relative}: {exc}")
            continue
        if var1_config_hashes.get(relative) != source_digest:
            failures.append(f"FGC-1-VAR1 config hash is stale for {relative}")

    expected_var1_implementation_sources = (
        "src/recursive_horizons/fgc/action.py",
        "src/recursive_horizons/fgc/metric_variation.py",
        "scripts/reproduce_fgc_metric_variation.py",
    )
    var1_implementation_hashes = var1_record.get(
        "implementation_sources_sha256", {}
    )
    if set(var1_implementation_hashes) != set(
        expected_var1_implementation_sources
    ):
        failures.append("FGC-1-VAR1 implementation hash ledger has unexpected keys")
    for relative in expected_var1_implementation_sources:
        try:
            source_digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-VAR1 source {relative}: {exc}")
            continue
        if var1_implementation_hashes.get(relative) != source_digest:
            failures.append(f"FGC-1-VAR1 source hash is stale for {relative}")

    var1_document_relative = "docs/fgc-metric-variation.md"
    if var1_record.get("derivation_document") != var1_document_relative:
        failures.append("FGC-1-VAR1 does not name its derivation document")
    try:
        var1_document_digest = sha256(
            (REPOSITORY / var1_document_relative).read_bytes()
        ).hexdigest()
    except OSError as exc:
        failures.append(f"cannot hash FGC-1-VAR1 derivation document: {exc}")
    else:
        if var1_record.get("derivation_document_sha256") != var1_document_digest:
            failures.append("FGC-1-VAR1 derivation-document hash is stale")

    expected_var1_fixture_contract = {
        "dimension": 4,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "algorithm": "lcg64_v1",
        "seed": 20260821,
        "fixture_count": 8,
        "integer_bound": 3,
        "kulkarni_nomizu_pair_count": 2,
    }
    if var1_record.get("fixture_contract") != expected_var1_fixture_contract:
        failures.append("FGC-1-VAR1 fixture contract is not the frozen exact contract")
    expected_var1_coefficient_order = [
        "riemann_hessian",
        "metric_ricci_hessian",
        "ricci_b_hessian_a",
        "ricci_a_hessian_b",
        "ricci_ab_box",
        "scalar_metric_box",
        "scalar_hessian",
    ]
    if var1_record.get("expanded_gb_coefficient_order") != (
        expected_var1_coefficient_order
    ):
        failures.append("FGC-1-VAR1 expanded coefficient order has drifted")
    expected_primary_source_crosscheck = {
        "source": "Thaalba_et_al_arXiv_2306.01695_equations_1_and_3",
        "url": "https://arxiv.org/abs/2306.01695",
        "mapping": "f_T=alpha_T*phi^2/2_and_double_epsilon_term_equals_-8*P_acbd*nabla^c*nabla^d*f_T",
        "coefficient_and_sign_match": True,
        "unsafe_normalization_warning": "Lara_et_al_arXiv_2403.08705_equation_6_is_not_used_as_coefficient_authority",
    }
    if var1_record.get("primary_source_crosscheck") != (
        expected_primary_source_crosscheck
    ):
        failures.append("FGC-1-VAR1 primary-source mapping has drifted")

    var1_certificate = var1_record.get("variation_certificate", {})
    convention_equations = var1_certificate.get("convention_locked_equations", {})
    expected_var1_coefficients = {
        "riemann_hessian": "-8",
        "metric_ricci_hessian": "8",
        "ricci_b_hessian_a": "-8",
        "ricci_a_hessian_b": "-8",
        "ricci_ab_box": "8",
        "scalar_metric_box": "-4",
        "scalar_hessian": "4",
    }
    if convention_equations.get("compact_gb_source") != (
        "T_GB_ab=-8*P_acbd*nabla^c*nabla^d*f"
    ):
        failures.append("FGC-1-VAR1 compact Gauss-Bonnet source has drifted")
    if convention_equations.get("full_metric_equation") != (
        "F*G_ab+(g_ab*Box-nabla_a*nabla_b)F=T_phi_ab+T_chi_ab-8*P_acbd*nabla^c*nabla^d*f"
    ):
        failures.append("FGC-1-VAR1 full metric equation has drifted")
    if convention_equations.get("expanded_gb_coefficients") != (
        expected_var1_coefficients
    ):
        failures.append("FGC-1-VAR1 expanded Gauss-Bonnet coefficients have drifted")

    expected_fixture_ids = [
        "lcg64-00-draw-01",
        "lcg64-01-draw-02",
        "lcg64-02-draw-03",
        "lcg64-03-draw-04",
        "lcg64-04-draw-05",
        "lcg64-05-draw-08",
        "lcg64-06-draw-09",
        "lcg64-07-draw-10",
    ]
    expected_fixture_hashes = [
        "6551b8447a137c622cb2aa057ef3a7ae470c9a1911465e42d65d2188ffe94006",
        "d17d4611c389098558595e84c82485d9d2853904c243f1956e918a6f1ee8e925",
        "8261a44a5b2aebc2e4d75e4d3873ad5dbc4289611b6b57cc59da5a4349571b5f",
        "4687abf3ceaa7c6508e740b2b3213237cc93aa4706b461d96d6be94c43dcff0a",
        "a27edf21354514cf6931ef8be6a7c2a6a66151329be0fb08c0b3c955064e3780",
        "c0bfd8f9e4b28e2b2b1d301806ec25b76aebe15a75c89cd008c3d43fdda6ff98",
        "0d80982995f7dd66d51d7e17a8f374ead579eef73239d8680ef2af8c05e782aa",
        "3bf066766a6c5cd0eb0b1ac486d913d46a51e448df3fb2eef74b484b5690ce76",
    ]
    var1_fixtures = var1_certificate.get("fixture_records", [])
    if not isinstance(var1_fixtures, list) or [
        fixture.get("fixture_id") if isinstance(fixture, dict) else None
        for fixture in var1_fixtures
    ] != expected_fixture_ids:
        failures.append("FGC-1-VAR1 fixture identifiers have drifted")
    elif [fixture.get("fixture_source_sha256") for fixture in var1_fixtures] != (
        expected_fixture_hashes
    ):
        failures.append("FGC-1-VAR1 exact fixture hashes have drifted")
    else:
        expected_validation_keys = {
            "metric_is_minkowski_diagonal",
            "riemann_first_pair_antisymmetric",
            "riemann_second_pair_antisymmetric",
            "riemann_pair_exchange_symmetric",
            "riemann_first_bianchi",
            "ricci_symmetric",
            "nonzero_curvature",
            "nonzero_gauss_bonnet_invariant",
            "nonzero_scalar_field",
            "nonzero_scalar_gradient",
            "nonzero_scalar_hessian",
            "generic_gb_source_nonzero",
            "generic_gb_source_symmetric",
            "p_trace_max_exact_residual",
            "noether_curvature_contraction_max_exact_residual",
        }
        for fixture in var1_fixtures:
            fixture_id = fixture["fixture_id"]
            if fixture.get("generic_compact_expanded_max_exact_residual") != "0":
                failures.append(f"FGC-1-VAR1 {fixture_id} has a GB expansion residual")
            if fixture.get("trace_identity_exact_residual") != "0":
                failures.append(f"FGC-1-VAR1 {fixture_id} has a trace residual")
            validation = fixture.get("validation", {})
            if set(validation) != expected_validation_keys:
                failures.append(
                    f"FGC-1-VAR1 {fixture_id} validation ledger has unexpected keys"
                )
            for key, value in validation.items():
                expected = "0" if key.endswith("_residual") else True
                if value != expected:
                    failures.append(
                        f"FGC-1-VAR1 {fixture_id} loses validation {key}={expected!r}"
                    )
            branches = fixture.get("branches", [])
            if [
                branch.get("model_id") if isinstance(branch, dict) else None
                for branch in branches
            ] != ["GR-0", "SGB-L", "FGC-QR"]:
                failures.append(f"FGC-1-VAR1 {fixture_id} loses branch order")
                continue
            for branch in branches:
                if branch.get("gb_compact_expanded_max_exact_residual") != "0":
                    failures.append(
                        f"FGC-1-VAR1 {fixture_id} {branch['model_id']} has a GB residual"
                    )
                if branch.get("nonminimal_compact_expanded_max_exact_residual") != "0":
                    failures.append(
                        f"FGC-1-VAR1 {fixture_id} {branch['model_id']} has a nonminimal residual"
                    )

    expected_branch_coverage = {
        "GR-0": {"gb_zero": True, "nonminimal_zero": True},
        "SGB-L": {"gb_nonzero": True, "nonminimal_zero": True},
        "FGC-QR": {"gb_nonzero": True, "nonminimal_nonzero": True},
    }
    if var1_certificate.get("branch_coverage") != expected_branch_coverage:
        failures.append("FGC-1-VAR1 branch source coverage is incomplete")
    var1_checks = var1_certificate.get("verified_exact_checks", {})
    expected_var1_checks = {
        "all_algebraic_curvature_fixtures_valid",
        "gb_metric_source_symmetric",
        "compact_expanded_gb_source_identity",
        "quadratic_nonminimal_source_identity",
        "p_trace_identity",
        "trace_identity",
        "noether_curvature_contraction_identity",
        "flrw_lapse_sign_and_factor",
        "all_branch_source_coverage",
        "constant_f_topological_bulk_source_zero",
        "all_var1_nonclaims_false",
    }
    if set(var1_checks) != expected_var1_checks:
        failures.append("FGC-1-VAR1 exact-check ledger has unexpected keys")
    for key in expected_var1_checks:
        if var1_checks.get(key) is not True:
            failures.append(f"FGC-1-VAR1 loses exact check {key}=true")
    var1_aggregate = var1_certificate.get("aggregate", {})
    if var1_aggregate.get("fixture_count") != 8:
        failures.append("FGC-1-VAR1 aggregate fixture count is not eight")
    if var1_aggregate.get("maximum_exact_residual") != "0":
        failures.append("FGC-1-VAR1 aggregate exact residual is nonzero")
    expected_residual_keys = {
        "compact_expanded_gb_source",
        "quadratic_nonminimal_source",
        "p_trace",
        "trace",
        "noether_curvature_contraction",
        "flrw_lapse",
    }
    residuals_by_identity = var1_aggregate.get(
        "maximum_exact_residual_by_identity", {}
    )
    if set(residuals_by_identity) != expected_residual_keys or any(
        value != "0" for value in residuals_by_identity.values()
    ):
        failures.append("FGC-1-VAR1 per-identity residual ledger is incomplete or nonzero")
    flrw_check = var1_certificate.get("flrw_lapse_crosscheck", {})
    if (
        flrw_check.get("compact_T_GB_00") != "-405/7"
        or flrw_check.get("independent_lapse_variation_T_GB_00") != "-405/7"
        or flrw_check.get("boundary_reduced_lagrangian") != "-1323/25"
        or flrw_check.get("exact_automatic_lapse_derivative") != "3969/25"
        or flrw_check.get("scale_factor_velocity") != "21/10"
        or flrw_check.get("stress_conversion_factor") != "-N^2/a^3"
        or flrw_check.get("max_exact_residual") != "0"
    ):
        failures.append("FGC-1-VAR1 FLRW lapse sign/factor cross-check has drifted")

    var1_gate_status = var1_record.get("gate_status", {})
    for key in (
        "act1_action_and_scalar_algebra_certificate_passed",
        "var1_covariant_metric_equation_derived_and_cross_checked",
    ):
        if var1_gate_status.get(key) is not True:
            failures.append(f"FGC-1-VAR1 must retain {key}=true")
    for key in (
        "full_fgc1_health_gate_passed",
        "publication_local_defocusing_gate_passed",
    ):
        if var1_gate_status.get(key) is not False:
            failures.append(f"FGC-1-VAR1 must retain {key}=false")
    required_var1_nonclaims = (
        "tensor_cas_or_second_independent_functional_variation_performed",
        "differential_bianchi_or_full_divergence_fixture_evaluated",
        "boundary_or_junction_variation_completed",
        "coupled_metric_chi_background_verified",
        "spherical_field_equations_derived",
        "spherical_principal_symbol_derived",
        "nonlinear_strong_hyperbolicity_proven",
        "full_coupled_quadratic_action_diagonalized",
        "nonlinear_ghost_freedom_proven",
        "constraint_propagation_proven",
        "metric_null_defocusing_derived",
        "finite_collapse_solution_constructed",
        "singularity_resolution_proven",
        "child_spacetime_constructed",
        "shear_robustness_proven",
        "global_flux_entropy_closure_derived",
        "dark_matter_derived",
        "dark_energy_derived",
        "variable_speed_of_light_derived",
        "particle_spectrum_derived",
        "theory_of_everything_derived",
        "observational_confirmation_obtained",
    )
    if set(var1_record.get("nonclaims", {})) != set(required_var1_nonclaims):
        failures.append("FGC-1-VAR1 nonclaim ledger has unexpected keys")
    for key in required_var1_nonclaims:
        if var1_record.get("nonclaims", {}).get(key) is not False:
            failures.append(f"FGC-1-VAR1 must explicitly retain {key}=false")
    if var1_certificate.get("nonclaims") != var1_record.get("nonclaims"):
        failures.append("FGC-1-VAR1 certificate and top-level nonclaims differ")

    red1_path = REPOSITORY / "results" / "fgc-1-hyp1-reduction.json"
    try:
        red1_record = json.loads(red1_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse FGC-1-HYP1-RED1 result: {exc}")
        return failures
    if (
        red1_record.get("project_version") != "0.11.0"
        or red1_record.get("schema_version") != 1
    ):
        failures.append("FGC-1-HYP1-RED1 result has stale version or schema")
    if red1_record.get("artifact_id") != "FGC-1-HYP1-RED1":
        failures.append("FGC-1-HYP1-RED1 result is missing its artifact identifier")
    if red1_record.get("artifact") != (
        "fgc1_spherical_reduction_and_principal_part_preflight_certificate"
    ):
        failures.append("FGC-1-HYP1-RED1 result is missing its artifact label")
    if red1_record.get("classification") != (
        "exact_spherical_reduction_and_principal_part_preflight_not_hyperbolicity_or_evolution"
    ):
        failures.append("FGC-1-HYP1-RED1 has promoted its restricted classification")
    if red1_record.get("development_status") != "pre_publication_gate":
        failures.append("FGC-1-HYP1-RED1 development status has drifted")
    if red1_record.get("source_config") != (
        "configs/fgc/fgc-1-hyp1-reduction.toml"
    ):
        failures.append("FGC-1-HYP1-RED1 does not name its frozen source config")

    expected_red1_config_sources = {
        "configs/fgc/fgc-1-action-gate.toml",
        "configs/fgc/fgc-1-metric-variation.toml",
        "configs/fgc/fgc-1-hyp1-reduction.toml",
    }
    red1_config_hashes = red1_record.get("source_configs_sha256", {})
    if set(red1_config_hashes) != expected_red1_config_sources:
        failures.append("FGC-1-HYP1-RED1 source-config hash ledger has unexpected keys")
    for relative in expected_red1_config_sources:
        try:
            source_digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-RED1 config {relative}: {exc}")
            continue
        if red1_config_hashes.get(relative) != source_digest:
            failures.append(f"FGC-1-HYP1-RED1 config hash is stale for {relative}")

    expected_red1_implementation_sources = {
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/__init__.py",
        "scripts/reproduce_fgc_hyp1_reduction.py",
    }
    red1_implementation_hashes = red1_record.get(
        "implementation_sources_sha256", {}
    )
    if set(red1_implementation_hashes) != expected_red1_implementation_sources:
        failures.append("FGC-1-HYP1-RED1 implementation hash ledger has unexpected keys")
    for relative in expected_red1_implementation_sources:
        try:
            source_digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-RED1 source {relative}: {exc}")
            continue
        if red1_implementation_hashes.get(relative) != source_digest:
            failures.append(f"FGC-1-HYP1-RED1 source hash is stale for {relative}")

    red1_document_relative = "docs/fgc-hyp1-reduction.md"
    if red1_record.get("derivation_document") != red1_document_relative:
        failures.append("FGC-1-HYP1-RED1 does not name its derivation document")
    try:
        red1_document_digest = sha256(
            (REPOSITORY / red1_document_relative).read_bytes()
        ).hexdigest()
    except OSError as exc:
        failures.append(f"cannot hash FGC-1-HYP1-RED1 derivation document: {exc}")
    else:
        if red1_record.get("derivation_document_sha256") != red1_document_digest:
            failures.append("FGC-1-HYP1-RED1 derivation-document hash is stale")

    expected_red1_formulation = {
        "gauge_id": "generalized_radial_adm_dynamic_lambda_unfixed_areal_radius_v1",
        "primary_formulation": "generalized_radial_adm_dynamic_lambda_unfixed_areal_radius",
        "crosscheck_formulation": "general_2plus2_warped_product",
        "time_slices": "spacelike_horizon_penetrating_compatible_when_regular",
        "radial_gauge": "unfixed_in_red1",
        "physical_metric": "g_ab_matter_coupled_metric",
        "equation_modification": "none_unredefined_covariant_equations",
    }
    if red1_record.get("formulation") != expected_red1_formulation:
        failures.append("FGC-1-HYP1-RED1 formulation has drifted")
    expected_red1_fields = [
        "alpha",
        "shift",
        "lambda",
        "areal_radius",
        "phi",
        "chi",
    ]
    expected_red1_equations = [
        "metric_tt",
        "metric_tr",
        "metric_rr",
        "metric_theta_theta",
        "scalar_phi",
        "scalar_chi",
    ]
    expected_red1_ordering = {
        "field_order": expected_red1_fields,
        "independent_equation_projection_order": expected_red1_equations,
        "candidate_constraint_projection_order": ["metric_tt", "metric_tr"],
        "principal_variable_order": expected_red1_fields,
    }
    if red1_record.get("ordering") != expected_red1_ordering:
        failures.append("FGC-1-HYP1-RED1 ordering has drifted")

    expected_red1_fixture_ids = [
        "GR0_regular",
        "SGBL_constant_f",
        "FGCQR_schwarzschild_exterior",
        "FGCQR_activated_generic",
    ]
    expected_red1_models = ["GR-0", "SGB-L", "FGC-QR", "FGC-QR"]
    red1_input_fixtures = red1_record.get("fixtures", [])
    if not isinstance(red1_input_fixtures, list) or [
        fixture.get("fixture_id") if isinstance(fixture, dict) else None
        for fixture in red1_input_fixtures
    ] != expected_red1_fixture_ids:
        failures.append("FGC-1-HYP1-RED1 input fixture identifiers have drifted")
    elif [fixture.get("model_id") for fixture in red1_input_fixtures] != (
        expected_red1_models
    ):
        failures.append("FGC-1-HYP1-RED1 input branch order has drifted")

    red1_certificate = red1_record.get("reduction_certificate", {})
    if red1_certificate.get("schema_version") != 1 or red1_certificate.get(
        "artifact_id"
    ) != "FGC-1-HYP1-RED1":
        failures.append("FGC-1-HYP1-RED1 nested certificate identity has drifted")
    if red1_certificate.get("classification") != (
        "exact_local_jet_reduction_not_first_order_or_hyperbolicity"
    ):
        failures.append("FGC-1-HYP1-RED1 nested classification has drifted")

    expected_red1_checks = {
        "all_direct_warped_regular_jet_residuals_zero",
        "all_core_nonclaims_false",
        "frozen_fixture_count_retained",
        "all_fixture_principal_matrices_six_by_eighteen",
        "all_chi_principal_factors_metric_null",
        "activated_mixing_control_present_and_exercised",
        "pg_schwarzschild_vacuum_metric_residual_zero",
        "flat_flrw_analytic_curvature_control_exact",
    }
    red1_checks = red1_certificate.get("verified_exact_checks", {})
    if set(red1_checks) != expected_red1_checks:
        failures.append("FGC-1-HYP1-RED1 exact-check ledger has unexpected keys")
    for key in expected_red1_checks:
        if red1_checks.get(key) is not True:
            failures.append(f"FGC-1-HYP1-RED1 loses exact check {key}=true")

    expected_red1_aggregate = {
        "fixture_count": 4,
        "all_direct_warped_exact": True,
        "all_principal_matrices_six_by_eighteen": True,
        "all_chi_principal_factors_metric_null": True,
        "schwarzschild_control_present_and_exact": True,
        "activated_mixing_control_present_and_exercised": True,
        "flat_flrw_curvature_control_exact": True,
        "all_nonclaims_false": True,
    }
    if red1_certificate.get("aggregate") != expected_red1_aggregate:
        failures.append("FGC-1-HYP1-RED1 aggregate has drifted")

    flrw = red1_certificate.get("analytic_curvature_controls", {}).get(
        "spatially_flat_flrw", {}
    )
    expected_flrw = {
        "inputs": {
            "scale_factor": "2",
            "hubble": "3/2",
            "hubble_derivative": "-2/3",
            "radial_coordinate": "5",
        },
        "expected_R": "23",
        "actual_R": "23",
        "expected_GB": "171/2",
        "actual_GB": "171/2",
        "exact_match": True,
    }
    if flrw != expected_flrw:
        failures.append("FGC-1-HYP1-RED1 flat-FLRW control has drifted")

    expected_second_derivatives = ["dtt", "dtr", "drr"]
    expected_red1_columns = [
        f"{field}.{derivative}"
        for field in expected_red1_fields
        for derivative in expected_second_derivatives
    ]
    expected_red1_matrix_hashes = {
        "GR0_regular": "376c436c698213276244b98e4f7b2b00f399ca75fda09e23e94dad5f7061725d",
        "SGBL_constant_f": "a3cd57f7b03e7cd727fa4f5b16e212fa90ddd415513ce846a3b173f013ad4d91",
        "FGCQR_schwarzschild_exterior": "0087c6f8482bf92739156011aa6005aa28ba8395b5904a7960080fafec485d0c",
        "FGCQR_activated_generic": "03de951411c91edc408ca551dc4368867b84add96665907185544cc0bb2d4c02",
    }
    expected_red1_core_nonclaims = {
        "first_order_reduction_derived",
        "evolution_constraint_split_derived",
        "constraint_propagation_proven",
        "kinetic_matrix_invertible_on_retained_domain",
        "physical_characteristic_polynomial_derived",
        "all_characteristics_real_on_retained_domain",
        "complete_characteristic_basis_or_symmetrizer_proven",
        "strong_hyperbolicity_proven",
        "evolution_authorized",
        "metric_null_defocusing_derived",
    }
    red1_fixtures = red1_certificate.get("fixtures", [])
    if not isinstance(red1_fixtures, list) or [
        fixture.get("fixture_id") if isinstance(fixture, dict) else None
        for fixture in red1_fixtures
    ] != expected_red1_fixture_ids:
        failures.append("FGC-1-HYP1-RED1 certificate fixture identifiers have drifted")
    else:
        if [fixture.get("branch") for fixture in red1_fixtures] != expected_red1_models:
            failures.append("FGC-1-HYP1-RED1 certificate branch order has drifted")
        for fixture in red1_fixtures:
            fixture_id = fixture["fixture_id"]
            residuals = fixture.get("direct_warped_max_exact_residuals", {})
            expected_residuals = {
                "riemann",
                "ricci",
                "R",
                "GB",
                "connection",
                "hessian_phi",
                "hessian_chi",
                "metric_equation",
                "scalar_phi",
                "scalar_chi",
            }
            if set(residuals) != expected_residuals or any(
                value != "0" for value in residuals.values()
            ):
                failures.append(
                    f"FGC-1-HYP1-RED1 {fixture_id} direct/warped ledger is incomplete or nonzero"
                )
            principal = fixture.get("principal_part", {})
            matrix = principal.get("matrix", [])
            if (
                principal.get("field_order") != expected_red1_fields
                or principal.get("second_derivative_order")
                != expected_second_derivatives
                or principal.get("column_order") != expected_red1_columns
                or principal.get("equation_order") != expected_red1_equations
                or principal.get("row_count") != 6
                or principal.get("column_count") != 18
                or not isinstance(matrix, list)
                or len(matrix) != 6
                or any(not isinstance(row, list) or len(row) != 18 for row in matrix)
            ):
                failures.append(
                    f"FGC-1-HYP1-RED1 {fixture_id} principal matrix contract has drifted"
                )
            else:
                canonical = json.dumps(matrix, separators=(",", ":"), ensure_ascii=True)
                matrix_digest = sha256(canonical.encode("ascii")).hexdigest()
                if principal.get("matrix_sha256") != matrix_digest:
                    failures.append(
                        f"FGC-1-HYP1-RED1 {fixture_id} matrix digest does not match its matrix"
                    )
                if matrix_digest != expected_red1_matrix_hashes[fixture_id]:
                    failures.append(
                        f"FGC-1-HYP1-RED1 {fixture_id} frozen matrix has drifted"
                    )
                for row in matrix:
                    for value in row:
                        if not isinstance(value, str):
                            failures.append(
                                f"FGC-1-HYP1-RED1 {fixture_id} matrix contains a non-exact value"
                            )
                            break
                        try:
                            numerator, separator, denominator = value.partition("/")
                            parsed_numerator = int(numerator)
                            parsed_denominator = int(denominator) if separator else 1
                            if parsed_denominator <= 0:
                                raise ValueError
                            divisor = gcd(abs(parsed_numerator), parsed_denominator)
                            canonical_value = (
                                str(parsed_numerator // divisor)
                                if parsed_denominator // divisor == 1
                                else f"{parsed_numerator // divisor}/{parsed_denominator // divisor}"
                            )
                        except (TypeError, ValueError, ZeroDivisionError):
                            failures.append(
                                f"FGC-1-HYP1-RED1 {fixture_id} matrix has invalid rational encoding"
                            )
                            break
                        if value != canonical_value:
                            failures.append(
                                f"FGC-1-HYP1-RED1 {fixture_id} matrix has noncanonical rational encoding"
                            )
                            break
            chi_control = fixture.get("chi_principal_control", {})
            if (
                chi_control.get("chi_principal_factor_is_metric_null") is not True
                or chi_control.get("chi_second_derivatives_absent_from_other_equations")
                is not True
                or chi_control.get("actual_chi_equation_coefficients")
                != chi_control.get("expected_metric_null_coefficients")
            ):
                failures.append(
                    f"FGC-1-HYP1-RED1 {fixture_id} loses the independent-chi metric-null control"
                )
            core_nonclaims = fixture.get("nonclaims", {})
            if set(core_nonclaims) != expected_red1_core_nonclaims or any(
                value is not False for value in core_nonclaims.values()
            ):
                failures.append(
                    f"FGC-1-HYP1-RED1 {fixture_id} has promoted a core nonclaim"
                )

        schwarzschild = red1_fixtures[2]
        if schwarzschild.get("schwarzschild_vacuum_control") != {
            "metric_residual_zero": True,
            "scalar_residuals_zero": True,
            "maximum_abs_metric_residual": "0",
            "absolute_phi_residual": "0",
            "absolute_chi_residual": "0",
            "ricci_scalar_zero": True,
            "gauss_bonnet_invariant": "3/16384",
        }:
            failures.append("FGC-1-HYP1-RED1 Schwarzschild control has drifted")
        if schwarzschild.get("invariants") != {
            "F": "4",
            "R": "0",
            "GB": "3/16384",
        }:
            failures.append("FGC-1-HYP1-RED1 Schwarzschild invariants have drifted")
        activated = red1_fixtures[3]
        if activated.get("activated_mixing_control") != {
            "metric_equations_depend_on_phi_second_jets": True,
            "phi_equation_depends_on_metric_second_jets": True,
            "nonzero_phi_value_gradient_and_hessian": True,
        }:
            failures.append("FGC-1-HYP1-RED1 activated mixing control has drifted")
        expected_source_sector_controls = [
            {"nonminimal_term_zero": True, "gb_bulk_residual_term_zero": True},
            {"nonminimal_term_zero": True, "gb_bulk_residual_term_zero": True},
            {"nonminimal_term_zero": True, "gb_bulk_residual_term_zero": True},
            {"nonminimal_term_zero": False, "gb_bulk_residual_term_zero": False},
        ]
        if [
            fixture.get("source_sector_controls") for fixture in red1_fixtures
        ] != expected_source_sector_controls:
            failures.append("FGC-1-HYP1-RED1 source-sector controls have drifted")

    expected_red1_gate_status = {
        "spherical_covariant_equations_evaluated_at_local_jets": True,
        "spherical_reduction_crosschecked_at_regular_jets": True,
        "independent_and_candidate_constraint_projection_order_frozen": True,
        "radial_constraint_projections_identified": False,
        "uneliminated_second_order_principal_part_extracted": True,
        "full_fgc1_health_gate_passed": False,
        "publication_local_defocusing_gate_passed": False,
    }
    if red1_record.get("gate_status") != expected_red1_gate_status:
        failures.append("FGC-1-HYP1-RED1 gate-status ledger has drifted")

    expected_red1_nonclaims = {
        "first_order_reduction_constructed",
        "constraint_propagation_proven",
        "kinetic_matrix_invertible_on_retained_domain",
        "all_characteristics_real_on_retained_domain",
        "complete_characteristic_basis_or_symmetrizer_proven",
        "nonlinear_strong_hyperbolicity_proven",
        "full_coupled_quadratic_action_diagonalized",
        "nonlinear_ghost_freedom_proven",
        "evolution_authorized",
        "finite_collapse_solution_constructed",
        "metric_null_defocusing_derived",
        "singularity_resolution_proven",
        "child_spacetime_constructed",
        "shear_robustness_proven",
        "global_flux_entropy_closure_derived",
        "dark_matter_derived",
        "dark_energy_derived",
        "variable_speed_of_light_derived",
        "particle_spectrum_derived",
        "theory_of_everything_derived",
        "observational_confirmation_obtained",
    }
    red1_nonclaims = red1_record.get("nonclaims", {})
    if set(red1_nonclaims) != expected_red1_nonclaims or any(
        value is not False for value in red1_nonclaims.values()
    ):
        failures.append("FGC-1-HYP1-RED1 top-level nonclaim ledger has drifted")
    core_nonclaims = red1_certificate.get("nonclaims", {})
    if set(core_nonclaims) != expected_red1_core_nonclaims or any(
        value is not False for value in core_nonclaims.values()
    ):
        failures.append("FGC-1-HYP1-RED1 certificate nonclaim ledger has drifted")

    sym1_path = REPOSITORY / "results" / "fgc-1-hyp1-symbol.json"
    try:
        sym1_record = json.loads(sym1_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse FGC-1-HYP1-SYM1 result: {exc}")
        return failures
    if (
        sym1_record.get("project_version") != "0.11.0"
        or sym1_record.get("schema_version") != 1
        or sym1_record.get("artifact_id") != "FGC-1-HYP1-SYM1"
    ):
        failures.append("FGC-1-HYP1-SYM1 identity, version, or schema has drifted")
    if sym1_record.get("classification") != (
        "exact_pointwise_covariant_scalar_quotient_and_adm_formulation_preflight_not_full_hyperbolicity_or_evolution"
    ):
        failures.append("FGC-1-HYP1-SYM1 has promoted its restricted classification")
    if sym1_record.get("development_status") != (
        "pointwise_physical_scalar_quotient_complete_full_hyp1_open"
    ):
        failures.append("FGC-1-HYP1-SYM1 development status has drifted")

    expected_sym1_configs = {
        "action": "configs/fgc/fgc-1-action-gate.toml",
        "variation": "configs/fgc/fgc-1-metric-variation.toml",
        "reduction": "configs/fgc/fgc-1-hyp1-reduction.toml",
        "symbol": "configs/fgc/fgc-1-hyp1-symbol.toml",
    }
    if sym1_record.get("source_configs") != expected_sym1_configs:
        failures.append("FGC-1-HYP1-SYM1 source config map has drifted")
    sym1_config_hashes = sym1_record.get("source_config_sha256", {})
    if set(sym1_config_hashes) != set(expected_sym1_configs.values()):
        failures.append("FGC-1-HYP1-SYM1 config hash ledger has unexpected keys")
    for relative in expected_sym1_configs.values():
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-SYM1 config {relative}: {exc}")
            continue
        if sym1_config_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-SYM1 config hash is stale for {relative}")

    expected_sym1_sources = {
        "src/recursive_horizons/fgc/exact_linear_algebra.py",
        "src/recursive_horizons/fgc/spherical_symbol.py",
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/__init__.py",
        "scripts/reproduce_fgc_action.py",
        "scripts/reproduce_fgc_metric_variation.py",
        "scripts/reproduce_fgc_hyp1_reduction.py",
        "scripts/reproduce_fgc_hyp1_symbol.py",
    }
    sym1_source_hashes = sym1_record.get("implementation_sha256", {})
    if set(sym1_source_hashes) != expected_sym1_sources:
        failures.append("FGC-1-HYP1-SYM1 implementation hash ledger has unexpected keys")
    for relative in expected_sym1_sources:
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-SYM1 source {relative}: {exc}")
            continue
        if sym1_source_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-SYM1 source hash is stale for {relative}")

    sym1_document = "docs/fgc-hyp1-symbol.md"
    if sym1_record.get("derivation_document") != sym1_document:
        failures.append("FGC-1-HYP1-SYM1 owner document has drifted")
    try:
        digest = sha256((REPOSITORY / sym1_document).read_bytes()).hexdigest()
    except OSError as exc:
        failures.append(f"cannot hash FGC-1-HYP1-SYM1 owner document: {exc}")
    else:
        if sym1_record.get("derivation_document_sha256") != digest:
            failures.append("FGC-1-HYP1-SYM1 owner-document hash is stale")

    sym1_certificate = sym1_record.get("symbol_certificate", {})
    canonical_certificate = json.dumps(
        sym1_certificate, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    if sym1_record.get("symbol_certificate_sha256") != sha256(
        canonical_certificate.encode("utf-8")
    ).hexdigest():
        failures.append("FGC-1-HYP1-SYM1 certificate digest does not match its payload")
    expected_sym1_checks = {
        "right_gauge_null_identities",
        "left_bianchi_null_identities",
        "metric_schur_identities",
        "scalar_schur_offdiagonal_zero",
        "quotient_atlas_scalar_determinants_agree",
        "scalar_determinant_factorization",
        "metric_blocks_invertible_at_physical_roots",
        "positive_metric_null_discriminants",
        "positive_regulator_discriminants",
        "exact_positive_scalar_companion_symmetrizers",
        "candidate_constraints_free_of_normal_accelerations",
        "candidate_constraints_free_of_gauge_source_second_jets",
        "pointwise_metric_phi_kinetic_blocks_nonsingular",
        "naive_fixed_gauge_gr_comparator_rejected",
    }
    sym1_checks = sym1_certificate.get("verified_exact_checks", {})
    if set(sym1_checks) != expected_sym1_checks or any(
        value is not True for value in sym1_checks.values()
    ):
        failures.append("FGC-1-HYP1-SYM1 exact-check ledger has drifted")
    if sym1_certificate.get("all_declared_exact_checks_pass") is not True:
        failures.append("FGC-1-HYP1-SYM1 aggregate exact gate is not true")
    expected_sym1_fixture_ids = [
        "GR0_regular",
        "SGBL_constant_f",
        "FGCQR_schwarzschild_exterior",
        "FGCQR_activated_generic",
    ]
    expected_sym1_branches = ["GR-0", "SGB-L", "FGC-QR", "FGC-QR"]
    if sym1_record.get("input_fixture_ids") != expected_sym1_fixture_ids:
        failures.append("FGC-1-HYP1-SYM1 fixture order has drifted")
    if sym1_record.get("input_branch_order") != expected_sym1_branches:
        failures.append("FGC-1-HYP1-SYM1 branch order has drifted")

    sym1_fixtures = sym1_certificate.get("fixture_records", [])
    if not isinstance(sym1_fixtures, list) or [
        item.get("fixture_id") if isinstance(item, dict) else None
        for item in sym1_fixtures
    ] != expected_sym1_fixture_ids:
        failures.append("FGC-1-HYP1-SYM1 certificate fixtures have drifted")
    else:
        if [item.get("branch") for item in sym1_fixtures] != expected_sym1_branches:
            failures.append("FGC-1-HYP1-SYM1 certificate branches have drifted")
        for item in sym1_fixtures:
            fixture_id = item["fixture_id"]
            for key in (
                "right_gauge_null_identity_exact",
                "left_bianchi_null_identity_exact",
                "quotient_atlas_scalar_determinants_agree_exactly",
                "quotient_atlas_covers_patch_boundary",
                "scalar_determinant_factorization_exact",
                "metric_block_invertible_at_physical_roots",
            ):
                if item.get(key) is not True:
                    failures.append(f"FGC-1-HYP1-SYM1 {fixture_id} loses {key}=true")
            atlas = item.get("quotient_atlas", {})
            if set(atlas) != {"rr_patch", "tt_patch"}:
                failures.append(f"FGC-1-HYP1-SYM1 {fixture_id} quotient atlas has drifted")
            else:
                for patch in atlas.values():
                    schur = patch.get("metric_schur", {})
                    if (
                        schur.get("schur_determinant_identity") is not True
                        or schur.get("scalar_offdiagonal_numerators_zero") is not True
                    ):
                        failures.append(
                            f"FGC-1-HYP1-SYM1 {fixture_id} loses an exact Schur identity"
                        )
            adm = item.get("adm_normal_preflight", {})
            if (
                adm.get("candidate_constraints_free_of_normal_accelerations")
                is not True
                or adm.get("candidate_constraints_free_of_gauge_source_second_jets")
                is not True
                or adm.get("kinetic_determinant") in {None, "0"}
            ):
                failures.append(f"FGC-1-HYP1-SYM1 {fixture_id} ADM preflight has drifted")
            for factor_name in ("metric_null_certificate", "regulator_certificate"):
                factor = item.get(factor_name, {})
                try:
                    discriminant = Fraction(factor.get("discriminant", "0"))
                    determinant_value = Fraction(factor.get("symmetrizer_determinant", "0"))
                except (ValueError, ZeroDivisionError):
                    discriminant = determinant_value = Fraction(0)
                if (
                    discriminant <= 0
                    or determinant_value <= 0
                    or factor.get("product_is_symmetric") is not True
                    or factor.get("symmetrizer_positive_definite") is not True
                ):
                    failures.append(
                        f"FGC-1-HYP1-SYM1 {fixture_id} {factor_name} health control has drifted"
                    )

        controls = sym1_fixtures[:3]
        if any(
            item.get("regulator_cone_proportional_to_metric") is not True
            or item.get("cone_resultant") != "0"
            for item in controls
        ):
            failures.append("FGC-1-HYP1-SYM1 control cones no longer coincide")
        activated = sym1_fixtures[3]
        if (
            activated.get("regulator_cone_proportional_to_metric") is not False
            or activated.get("cones_share_no_root") is not True
            or activated.get("cone_resultant") in {None, "0"}
            or activated.get("metric_null_certificate", {}).get("roots")
            != [{"exact": "-31/42"}, {"exact": "59/42"}]
        ):
            failures.append("FGC-1-HYP1-SYM1 activated cone control has drifted")

    naive_sym1 = sym1_certificate.get("naive_fixed_gauge_gr_comparator", {})
    if (
        naive_sym1.get("zero_speed_geometric_multiplicity") != 1
        or naive_sym1.get("zero_speed_generalized_multiplicity") != 4
        or naive_sym1.get("defective_zero_speed_sector") is not True
        or naive_sym1.get("strong_hyperbolicity_rejected") is not True
    ):
        failures.append("FGC-1-HYP1-SYM1 rejected comparator has drifted")

    expected_sym1_gate_status = {
        "exact_covariant_physical_scalar_quotient_passed": True,
        "pointwise_scalar_companion_symmetrizers_passed": True,
        "adm_principal_constraint_kinetic_preflight_passed": True,
        "naive_fixed_gauge_comparator_rejected": True,
        "full_fgc1_hyp1_health_gate_passed": False,
        "evolution_authorized": False,
        "publication_local_defocusing_gate_passed": False,
    }
    if sym1_record.get("gate_status") != expected_sym1_gate_status:
        failures.append("FGC-1-HYP1-SYM1 gate-status ledger has drifted")
    expected_sym1_nonclaims = {
        "nonlinear_first_order_evolution_system_derived",
        "lower_order_radial_constraints_derived",
        "gauge_preservation_proven",
        "constraint_propagation_proven",
        "regular_center_boundary_system_derived",
        "open_retained_eft_domain_proven",
        "uniform_full_system_symmetrizer_proven",
        "full_strong_hyperbolicity_proven",
        "evolution_authorized",
        "collapse_solution_derived",
        "metric_null_affine_defocusing_derived",
    }
    sym1_nonclaims = sym1_record.get("nonclaims", {})
    if set(sym1_nonclaims) != expected_sym1_nonclaims or any(
        value is not False for value in sym1_nonclaims.values()
    ):
        failures.append("FGC-1-HYP1-SYM1 nonclaim ledger has drifted")
    if sym1_certificate.get("nonclaims") != sym1_nonclaims:
        failures.append("FGC-1-HYP1-SYM1 nested and top-level nonclaims differ")

    mhg1_path = REPOSITORY / "results" / "fgc-1-hyp1-modified-harmonic.json"
    try:
        mhg1_record = json.loads(mhg1_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse FGC-1-HYP1-MHG1 result: {exc}")
        return failures
    if (
        mhg1_record.get("project_version") != "0.11.0"
        or mhg1_record.get("schema_version") != 1
        or mhg1_record.get("artifact_id") != "FGC-1-HYP1-MHG1"
    ):
        failures.append("FGC-1-HYP1-MHG1 identity, version, or schema has drifted")
    if mhg1_record.get("classification") != (
        "exact_pointwise_modified_harmonic_spherical_principal_and_gauge_constraint_preflight_not_open_domain_or_evolution"
    ):
        failures.append("FGC-1-HYP1-MHG1 has promoted its restricted classification")
    if mhg1_record.get("development_status") != (
        "complete_pointwise_mhg_principal_preflight_full_open_domain_hyp1_open"
    ):
        failures.append("FGC-1-HYP1-MHG1 development status has drifted")

    expected_mhg1_configs = {
        "action": "configs/fgc/fgc-1-action-gate.toml",
        "variation": "configs/fgc/fgc-1-metric-variation.toml",
        "reduction": "configs/fgc/fgc-1-hyp1-reduction.toml",
        "symbol": "configs/fgc/fgc-1-hyp1-symbol.toml",
        "modified_harmonic": "configs/fgc/fgc-1-hyp1-modified-harmonic.toml",
    }
    if mhg1_record.get("source_configs") != expected_mhg1_configs:
        failures.append("FGC-1-HYP1-MHG1 source config map has drifted")
    mhg1_config_hashes = mhg1_record.get("source_config_sha256", {})
    if set(mhg1_config_hashes) != set(expected_mhg1_configs.values()):
        failures.append("FGC-1-HYP1-MHG1 config hash ledger has unexpected keys")
    for relative in expected_mhg1_configs.values():
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-MHG1 config {relative}: {exc}")
            continue
        if mhg1_config_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-MHG1 config hash is stale for {relative}")

    expected_mhg1_sources = {
        "src/recursive_horizons/fgc/exact_linear_algebra.py",
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/spherical_symbol.py",
        "src/recursive_horizons/fgc/modified_harmonic.py",
        "src/recursive_horizons/fgc/__init__.py",
        "scripts/reproduce_fgc_action.py",
        "scripts/reproduce_fgc_metric_variation.py",
        "scripts/reproduce_fgc_hyp1_reduction.py",
        "scripts/reproduce_fgc_hyp1_symbol.py",
        "scripts/reproduce_fgc_hyp1_modified_harmonic.py",
    }
    mhg1_source_hashes = mhg1_record.get("implementation_sha256", {})
    if set(mhg1_source_hashes) != expected_mhg1_sources:
        failures.append("FGC-1-HYP1-MHG1 implementation hash ledger has unexpected keys")
    for relative in expected_mhg1_sources:
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-MHG1 source {relative}: {exc}")
            continue
        if mhg1_source_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-MHG1 source hash is stale for {relative}")

    mhg1_document = "docs/fgc-hyp1-modified-harmonic.md"
    if mhg1_record.get("derivation_document") != mhg1_document:
        failures.append("FGC-1-HYP1-MHG1 owner document has drifted")
    try:
        digest = sha256((REPOSITORY / mhg1_document).read_bytes()).hexdigest()
    except OSError as exc:
        failures.append(f"cannot hash FGC-1-HYP1-MHG1 owner document: {exc}")
    else:
        if mhg1_record.get("derivation_document_sha256") != digest:
            failures.append("FGC-1-HYP1-MHG1 owner-document hash is stale")

    mhg1_certificate = mhg1_record.get("modified_harmonic_certificate", {})
    canonical_certificate = json.dumps(
        mhg1_certificate, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    if mhg1_record.get("modified_harmonic_certificate_sha256") != sha256(
        canonical_certificate.encode("utf-8")
    ).hexdigest():
        failures.append("FGC-1-HYP1-MHG1 certificate digest does not match its payload")
    expected_mhg1_checks = {
        "positive_effective_planck_coefficients",
        "gauge_fixing_symbols_nonzero",
        "complete_determinant_factorizations",
        "auxiliary_cones_nested_and_separated",
        "gauge_constraint_hat_cone_principal_propagation",
        "coordinate_time_kinetic_matrices_invertible",
        "first_order_characteristic_reductions_exact",
        "pointwise_complete_real_first_order_eigenbases",
        "direct_x_zero_var1_regularity_controls",
        "sym1_physical_factor_regression",
    }
    mhg1_checks = mhg1_certificate.get("verified_exact_checks", {})
    if set(mhg1_checks) != expected_mhg1_checks or any(
        value is not True for value in mhg1_checks.values()
    ):
        failures.append("FGC-1-HYP1-MHG1 exact-check ledger has drifted")
    if mhg1_certificate.get("all_declared_exact_checks_pass") is not True:
        failures.append("FGC-1-HYP1-MHG1 aggregate exact gate is not true")

    expected_mhg1_fixture_ids = [
        "GR0_regular",
        "SGBL_constant_f",
        "FGCQR_schwarzschild_exterior",
        "FGCQR_activated_generic",
    ]
    expected_mhg1_branches = ["GR-0", "SGB-L", "FGC-QR", "FGC-QR"]
    if mhg1_record.get("input_fixture_ids") != expected_mhg1_fixture_ids:
        failures.append("FGC-1-HYP1-MHG1 fixture order has drifted")
    if mhg1_record.get("input_branch_order") != expected_mhg1_branches:
        failures.append("FGC-1-HYP1-MHG1 branch order has drifted")
    formulation = mhg1_certificate.get("formulation", {})
    if (
        formulation.get("physical_equations") != "unredefined_ACT1_VAR1"
        or formulation.get("tilde_normal_factor") != "4"
        or formulation.get("hat_normal_factor") != "9"
        or formulation.get("first_order_reduction")
        != "standard_dt_dr_reduction_of_complete_six_field_spherical_second_order_system"
    ):
        failures.append("FGC-1-HYP1-MHG1 formulation contract has drifted")

    mhg1_fixtures = mhg1_certificate.get("fixture_records", [])
    if not isinstance(mhg1_fixtures, list) or [
        item.get("fixture_id") if isinstance(item, dict) else None
        for item in mhg1_fixtures
    ] != expected_mhg1_fixture_ids:
        failures.append("FGC-1-HYP1-MHG1 certificate fixtures have drifted")
    else:
        for item in mhg1_fixtures:
            fixture_id = item["fixture_id"]
            factorization = item.get("determinant_factorization", {})
            first_order = item.get("standard_first_order_principal_system", {})
            propagation = item.get(
                "gauge_constraint_propagation_principal_symbol", {}
            )
            if (
                item.get("effective_planck_coefficient") in {None, "0"}
                or item.get("gauge_fixing_nonzero") is not True
                or factorization.get("exact") is not True
                or factorization.get("tilde_multiplicity") != 2
                or factorization.get("hat_multiplicity") != 2
                or factorization.get("nonzero_constant_quotient") in {None, "0"}
                or item.get("auxiliary_cones_separated_from_all_physical_modes")
                is not True
                or item.get("gauge_constraint_propagation_hat_null_exact") is not True
                or first_order.get("coordinate_time_kinetic_determinant") in {None, "0"}
                or first_order.get("A_time_determinant") in {None, "0"}
                or item.get(
                    "coordinate_time_kinetic_matches_characteristic_leading_coefficient"
                )
                is not True
                or item.get("pointwise_complete_real_first_order_eigenbasis")
                is not True
            ):
                failures.append(
                    f"FGC-1-HYP1-MHG1 {fixture_id} loses a required exact principal gate"
                )
            expected_propagation_keys = {
                "input_constraint_components",
                "output_divergence_components",
                "full_projector_contraction",
                "expected_full_hat_wave_operator",
                "expected_hat_wave_factor",
                "expected_hat_wave_factor_monic",
                "spherical_projector_contraction",
                "full_projector_contraction_exact",
            }
            wave_factor = propagation.get("expected_hat_wave_factor")
            if (
                set(propagation) != expected_propagation_keys
                or propagation.get("input_constraint_components")
                != ["H^t", "H^r", "H^theta", "H^phi"]
                or propagation.get("output_divergence_components")
                != ["nu=t", "nu=r", "nu=theta", "nu=phi"]
                or propagation.get("full_projector_contraction_exact") is not True
                or propagation.get("full_projector_contraction")
                != propagation.get("expected_full_hat_wave_operator")
                or propagation.get("expected_hat_wave_factor_monic")
                != factorization.get("hat_null_monic")
                or not isinstance(wave_factor, list)
                or propagation.get("spherical_projector_contraction")
                != [[wave_factor, ["0"]], [["0"], wave_factor]]
            ):
                failures.append(
                    f"FGC-1-HYP1-MHG1 {fixture_id} projector-derived gauge propagation has drifted"
                )
            if (
                first_order.get("characteristic_pencil_determinant")
                != item.get("complete_mhg_determinant")
                or item.get(
                    "first_order_characteristic_determinant_equals_second_order"
                )
                is not True
            ):
                failures.append(
                    f"FGC-1-HYP1-MHG1 {fixture_id} first-order pencil determinant is not derived exactly"
                )
            cone_checks = item.get("auxiliary_cone_checks", {})
            if set(cone_checks) != {
                "physical_roots_are_tilde_timelike",
                "physical_roots_are_hat_timelike",
                "auxiliary_cones_do_not_intersect",
            } or any(value is not True for value in cone_checks.values()):
                failures.append(
                    f"FGC-1-HYP1-MHG1 {fixture_id} auxiliary cone checks have drifted"
                )
            modes = item.get("rational_characteristic_modes", {})
            expected_mode_dimensions = {
                "tilde_gauge": 2,
                "hat_constraint": 2,
                "physical_metric": 1 if fixture_id == "FGCQR_activated_generic" else 2,
            }
            if set(modes) != set(expected_mode_dimensions):
                failures.append(
                    f"FGC-1-HYP1-MHG1 {fixture_id} characteristic mode ledger has drifted"
                )
            else:
                for name, expected_dimension in expected_mode_dimensions.items():
                    mode = modes[name]
                    if (
                        mode.get("expected_kernel_dimension") != expected_dimension
                        or mode.get("kernel_dimensions")
                        != [expected_dimension, expected_dimension]
                        or mode.get("second_order_kernel_dimensions")
                        != [expected_dimension, expected_dimension]
                        or mode.get("first_order_kernel_dimensions")
                        != [expected_dimension, expected_dimension]
                        or mode.get("second_to_first_order_kernel_map_exact")
                        is not True
                        or mode.get("semisimple_at_each_rational_root") is not True
                    ):
                        failures.append(
                            f"FGC-1-HYP1-MHG1 {fixture_id} {name} eigenspace has drifted"
                        )
            resultants = item.get("cone_resultants", {})
            for name in (
                "tilde_vs_hat",
                "tilde_vs_physical",
                "hat_vs_physical",
                "tilde_vs_regulator",
                "hat_vs_regulator",
            ):
                if resultants.get(name) in {None, "0"}:
                    failures.append(
                        f"FGC-1-HYP1-MHG1 {fixture_id} loses {name} separation"
                    )

        x_zero_ids = [
            item["fixture_id"]
            for item in mhg1_fixtures
            if item.get("x_zero_direct_var1_coefficients_regular") is True
        ]
        if x_zero_ids != expected_mhg1_fixture_ids[:3]:
            failures.append("FGC-1-HYP1-MHG1 X=0 regularity controls have drifted")
        activated_mhg1 = mhg1_fixtures[3]
        activated_modes = activated_mhg1.get("rational_characteristic_modes", {})
        if (
            activated_modes.get("tilde_gauge", {}).get("roots")
            != ["-17/84", "73/84"]
            or activated_modes.get("hat_constraint", {}).get("roots")
            != ["-1/42", "29/42"]
            or activated_modes.get("physical_metric", {}).get("roots")
            != ["-31/42", "59/42"]
            or activated_mhg1.get("cone_resultants", {}).get(
                "physical_vs_regulator"
            )
            in {None, "0"}
        ):
            failures.append("FGC-1-HYP1-MHG1 activated exact cone control has drifted")

    expected_mhg1_gate_status = {
        "exact_pointwise_modified_harmonic_principal_gate_passed": True,
        "gauge_constraint_principal_propagation_gate_passed": True,
        "pointwise_complete_real_first_order_eigenbasis_passed": True,
        "direct_x_zero_regular_coefficient_gate_passed": True,
        "complete_lower_order_first_order_formulation_passed": False,
        "uniform_open_domain_hyperbolicity_gate_passed": False,
        "full_fgc1_hyp1_health_gate_passed": False,
        "evolution_authorized": False,
        "publication_local_defocusing_gate_passed": False,
    }
    if mhg1_record.get("gate_status") != expected_mhg1_gate_status:
        failures.append("FGC-1-HYP1-MHG1 gate-status ledger has drifted")
    expected_mhg1_nonclaims = {
        "complete_lower_order_first_order_sources_derived",
        "gauge_constraint_lower_order_coefficients_serialized",
        "reduction_constraint_propagation_proven",
        "uniform_open_domain_symmetrizer_proven",
        "open_retained_eft_domain_proven",
        "boundary_or_regular_center_system_derived",
        "evolution_authorized",
        "collapse_solution_derived",
        "metric_null_affine_defocusing_derived",
        "singularity_resolution_derived",
    }
    mhg1_nonclaims = mhg1_record.get("nonclaims", {})
    if set(mhg1_nonclaims) != expected_mhg1_nonclaims or any(
        value is not False for value in mhg1_nonclaims.values()
    ):
        failures.append("FGC-1-HYP1-MHG1 nonclaim ledger has drifted")
    if mhg1_certificate.get("nonclaims") != mhg1_nonclaims:
        failures.append("FGC-1-HYP1-MHG1 nested and top-level nonclaims differ")

    ref1_path = REPOSITORY / "results" / "fgc-1-hyp1-mhg-reference.json"
    try:
        ref1_record = json.loads(ref1_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse FGC-1-HYP1-MHG2-REF1 result: {exc}")
        return failures
    if (
        ref1_record.get("project_version") != "0.11.0"
        or ref1_record.get("schema_version") != 1
        or ref1_record.get("artifact_id") != "FGC-1-HYP1-MHG2-REF1"
    ):
        failures.append("FGC-1-HYP1-MHG2-REF1 identity, version, or schema has drifted")
    if ref1_record.get("classification") != (
        "exact_reference_connection_modified_harmonic_residual_not_implicit_solve_propagation_domain_or_evolution"
    ):
        failures.append("FGC-1-HYP1-MHG2-REF1 has promoted its restricted classification")
    if ref1_record.get("development_status") != (
        "complete_reference_gauge_residual_implicit_propagation_and_domain_gates_open"
    ):
        failures.append("FGC-1-HYP1-MHG2-REF1 development status has drifted")

    expected_ref1_configs = {
        "action": "configs/fgc/fgc-1-action-gate.toml",
        "variation": "configs/fgc/fgc-1-metric-variation.toml",
        "reduction": "configs/fgc/fgc-1-hyp1-reduction.toml",
        "symbol": "configs/fgc/fgc-1-hyp1-symbol.toml",
        "modified_harmonic": "configs/fgc/fgc-1-hyp1-modified-harmonic.toml",
        "reference": "configs/fgc/fgc-1-hyp1-mhg-reference.toml",
    }
    if ref1_record.get("source_configs") != expected_ref1_configs:
        failures.append("FGC-1-HYP1-MHG2-REF1 source config map has drifted")
    ref1_config_hashes = ref1_record.get("source_config_sha256", {})
    if set(ref1_config_hashes) != set(expected_ref1_configs.values()):
        failures.append("FGC-1-HYP1-MHG2-REF1 config hash ledger has unexpected keys")
    for relative in expected_ref1_configs.values():
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-MHG2-REF1 config {relative}: {exc}")
            continue
        if ref1_config_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-MHG2-REF1 config hash is stale for {relative}")

    expected_ref1_sources = {
        "src/recursive_horizons/fgc/exact_linear_algebra.py",
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/spherical_symbol.py",
        "src/recursive_horizons/fgc/modified_harmonic.py",
        "src/recursive_horizons/fgc/reference_connection.py",
        "src/recursive_horizons/fgc/modified_harmonic_reference.py",
        "src/recursive_horizons/fgc/__init__.py",
        "scripts/reproduce_fgc_action.py",
        "scripts/reproduce_fgc_metric_variation.py",
        "scripts/reproduce_fgc_hyp1_reduction.py",
        "scripts/reproduce_fgc_hyp1_symbol.py",
        "scripts/reproduce_fgc_hyp1_modified_harmonic.py",
        "scripts/reproduce_fgc_hyp1_mhg_reference.py",
    }
    ref1_source_hashes = ref1_record.get("implementation_sha256", {})
    if set(ref1_source_hashes) != expected_ref1_sources:
        failures.append("FGC-1-HYP1-MHG2-REF1 implementation hash ledger has unexpected keys")
    for relative in expected_ref1_sources:
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-MHG2-REF1 source {relative}: {exc}")
            continue
        if ref1_source_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-MHG2-REF1 source hash is stale for {relative}")

    ref1_document = "docs/fgc-hyp1-mhg-reference.md"
    if ref1_record.get("derivation_document") != ref1_document:
        failures.append("FGC-1-HYP1-MHG2-REF1 owner document has drifted")
    try:
        digest = sha256((REPOSITORY / ref1_document).read_bytes()).hexdigest()
    except OSError as exc:
        failures.append(f"cannot hash FGC-1-HYP1-MHG2-REF1 owner document: {exc}")
    else:
        if ref1_record.get("derivation_document_sha256") != digest:
            failures.append("FGC-1-HYP1-MHG2-REF1 owner-document hash is stale")

    ref1_certificate = ref1_record.get("reference_gauge_certificate", {})
    canonical_certificate = json.dumps(
        ref1_certificate, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    if ref1_record.get("reference_gauge_certificate_sha256") != sha256(
        canonical_certificate.encode("utf-8")
    ).hexdigest():
        failures.append("FGC-1-HYP1-MHG2-REF1 certificate digest does not match its payload")
    expected_ref1_checks = {
        "reference_annulus_strictly_positive",
        "positive_effective_planck_coefficients",
        "declared_gauge_surface_jets_exact",
        "declared_off_gauge_jets_exercise_extension",
        "full_residual_decomposition_exact",
        "scalar_equations_unmodified",
        "extension_tensors_symmetric_and_spherical",
        "reference_extension_directional_symbol_matches_mhg1_gauge_block",
        "schwarzschild_unredefined_solution_control_retained",
        "contravariant_spherical_gauge_components_only",
    }
    ref1_checks = ref1_certificate.get("verified_exact_checks", {})
    if set(ref1_checks) != expected_ref1_checks or any(
        value is not True for value in ref1_checks.values()
    ):
        failures.append("FGC-1-HYP1-MHG2-REF1 exact-check ledger has drifted")
    if ref1_certificate.get("all_declared_exact_checks_pass") is not True:
        failures.append("FGC-1-HYP1-MHG2-REF1 aggregate exact gate is not true")

    expected_ref1_fixture_ids = expected_mhg1_fixture_ids
    if ref1_record.get("input_fixture_ids") != expected_ref1_fixture_ids:
        failures.append("FGC-1-HYP1-MHG2-REF1 fixture order has drifted")
    if ref1_record.get("input_branch_order") != expected_mhg1_branches:
        failures.append("FGC-1-HYP1-MHG2-REF1 branch order has drifted")
    reference = ref1_record.get("reference", {})
    formulation = ref1_certificate.get("formulation", {})
    if (
        reference.get("reference_id") != "flat_spherical_annulus"
        or reference.get("radial_domain_minimum") != "1/2"
        or reference.get("center_included") is not False
        or formulation.get("physical_equations") != "unredefined_ACT1_VAR1"
        or formulation.get("reference_connection") != "flat_spherical_annulus"
        or formulation.get("tilde_normal_factor") != "4"
        or formulation.get("hat_normal_factor") != "9"
        or formulation.get("scalar_equations") != "unmodified"
    ):
        failures.append("FGC-1-HYP1-MHG2-REF1 formulation contract has drifted")

    ref1_fixtures = ref1_certificate.get("fixture_records", [])
    if not isinstance(ref1_fixtures, list) or [
        item.get("fixture_id") if isinstance(item, dict) else None
        for item in ref1_fixtures
    ] != expected_ref1_fixture_ids:
        failures.append("FGC-1-HYP1-MHG2-REF1 certificate fixtures have drifted")
    else:
        for index, item in enumerate(ref1_fixtures):
            fixture_id = item["fixture_id"]
            gauge_surface = index < 2
            if (
                item.get("effective_planck_coefficient") in {None, "0"}
                or item.get("full_equals_unredefined_plus_extension") is not True
                or item.get("scalar_equations_unmodified") is not True
                or item.get("extension_tensor_symmetric") is not True
                or item.get("equatorial_angular_isotropy") is not True
                or item.get("reference_extension_principal_regression_exact") is not True
                or item.get("reference_extension_principal_symbol")
                != item.get("mhg1_principal_gauge_symbol")
                or item.get("constraint_up", [None, None, None, None])[2:]
                != ["0", "0"]
            ):
                failures.append(
                    f"FGC-1-HYP1-MHG2-REF1 {fixture_id} loses a required exact residual gate"
                )
            if gauge_surface:
                if (
                    item.get("gauge_surface_at_evaluation_jet") is not True
                    or item.get("gauge_surface_extension_zero") is not True
                    or item.get("gauge_surface_full_equals_unredefined") is not True
                    or item.get("extension_nonzero") is not False
                ):
                    failures.append(
                        f"FGC-1-HYP1-MHG2-REF1 {fixture_id} gauge-surface control has drifted"
                    )
            elif (
                item.get("gauge_surface_at_evaluation_jet") is not False
                or item.get("extension_nonzero") is not True
            ):
                failures.append(
                    f"FGC-1-HYP1-MHG2-REF1 {fixture_id} off-gauge control has drifted"
                )
        schwarzschild_ref1 = ref1_fixtures[2]
        if (
            schwarzschild_ref1.get("constraint_up", [])[:2]
            != ["3/32", "1/64"]
            or any(
                value != "0"
                for value in schwarzschild_ref1.get(
                    "unredefined_residual_vector", []
                )
            )
        ):
            failures.append("FGC-1-HYP1-MHG2-REF1 Schwarzschild control has drifted")

    expected_ref1_gate_status = {
        "complete_reference_gauge_residual_passed": True,
        "exact_gauge_surface_equivalence_controls_passed": True,
        "mhg1_principal_gauge_block_regression_passed": True,
        "implicit_second_time_derivative_branch_passed": False,
        "complete_gauge_propagation_passed": False,
        "complete_first_order_reduction_passed": False,
        "uniform_open_domain_hyperbolicity_gate_passed": False,
        "full_fgc1_hyp1_health_gate_passed": False,
        "evolution_authorized": False,
        "publication_local_defocusing_gate_passed": False,
    }
    if ref1_record.get("gate_status") != expected_ref1_gate_status:
        failures.append("FGC-1-HYP1-MHG2-REF1 gate-status ledger has drifted")
    expected_ref1_nonclaims = {
        "implicit_second_time_derivative_branch_proven",
        "complete_lower_order_first_order_sources_derived",
        "gauge_constraint_lower_order_propagation_proven",
        "reduction_constraint_propagation_proven",
        "uniform_open_domain_symmetrizer_proven",
        "open_retained_eft_domain_proven",
        "boundary_or_regular_center_system_derived",
        "evolution_authorized",
        "collapse_solution_derived",
        "metric_null_affine_defocusing_derived",
        "singularity_resolution_derived",
    }
    ref1_nonclaims = ref1_record.get("nonclaims", {})
    if set(ref1_nonclaims) != expected_ref1_nonclaims or any(
        value is not False for value in ref1_nonclaims.values()
    ):
        failures.append("FGC-1-HYP1-MHG2-REF1 nonclaim ledger has drifted")
    if ref1_certificate.get("nonclaims") != ref1_nonclaims:
        failures.append("FGC-1-HYP1-MHG2-REF1 nested and top-level nonclaims differ")

    imp1_path = REPOSITORY / "results" / "fgc-1-hyp1-mhg-implicit.json"
    try:
        imp1_record = json.loads(imp1_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse FGC-1-HYP1-MHG3-IMP1 result: {exc}")
        return failures
    if (
        imp1_record.get("project_version") != "0.11.0"
        or imp1_record.get("schema_version") != 1
        or imp1_record.get("artifact_id") != "FGC-1-HYP1-MHG3-IMP1"
    ):
        failures.append("FGC-1-HYP1-MHG3-IMP1 identity, version, or schema has drifted")
    if imp1_record.get("classification") != (
        "exact_local_implicit_coordinate_time_acceleration_branch_at_flat_fgcqr_reference_gauge_root_not_open_domain_or_evolution"
    ):
        failures.append("FGC-1-HYP1-MHG3-IMP1 has promoted its restricted classification")
    if imp1_record.get("development_status") != (
        "flat_fgcqr_local_implicit_acceleration_branch_complete_propagation_domain_and_evolution_gates_open"
    ):
        failures.append("FGC-1-HYP1-MHG3-IMP1 development status has drifted")

    expected_imp1_configs = {
        "action": "configs/fgc/fgc-1-action-gate.toml",
        "variation": "configs/fgc/fgc-1-metric-variation.toml",
        "reduction": "configs/fgc/fgc-1-hyp1-reduction.toml",
        "symbol": "configs/fgc/fgc-1-hyp1-symbol.toml",
        "modified_harmonic": "configs/fgc/fgc-1-hyp1-modified-harmonic.toml",
        "reference": "configs/fgc/fgc-1-hyp1-mhg-reference.toml",
        "implicit": "configs/fgc/fgc-1-hyp1-mhg-implicit.toml",
    }
    if imp1_record.get("source_configs") != expected_imp1_configs:
        failures.append("FGC-1-HYP1-MHG3-IMP1 source config map has drifted")
    imp1_config_hashes = imp1_record.get("source_config_sha256", {})
    if set(imp1_config_hashes) != set(expected_imp1_configs.values()):
        failures.append("FGC-1-HYP1-MHG3-IMP1 config hash ledger has unexpected keys")
    for relative in expected_imp1_configs.values():
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-MHG3-IMP1 config {relative}: {exc}")
            continue
        if imp1_config_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-MHG3-IMP1 config hash is stale for {relative}")

    expected_imp1_sources = {
        "src/recursive_horizons/fgc/exact_linear_algebra.py",
        "src/recursive_horizons/fgc/exact_tangent.py",
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/spherical_symbol.py",
        "src/recursive_horizons/fgc/modified_harmonic.py",
        "src/recursive_horizons/fgc/reference_connection.py",
        "src/recursive_horizons/fgc/modified_harmonic_reference.py",
        "src/recursive_horizons/fgc/modified_harmonic_implicit.py",
        "src/recursive_horizons/fgc/__init__.py",
        "scripts/reproduce_fgc_action.py",
        "scripts/reproduce_fgc_metric_variation.py",
        "scripts/reproduce_fgc_hyp1_reduction.py",
        "scripts/reproduce_fgc_hyp1_symbol.py",
        "scripts/reproduce_fgc_hyp1_modified_harmonic.py",
        "scripts/reproduce_fgc_hyp1_mhg_reference.py",
        "scripts/reproduce_fgc_hyp1_mhg_implicit.py",
    }
    imp1_source_hashes = imp1_record.get("implementation_sha256", {})
    if set(imp1_source_hashes) != expected_imp1_sources:
        failures.append("FGC-1-HYP1-MHG3-IMP1 implementation hash ledger has unexpected keys")
    for relative in expected_imp1_sources:
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-MHG3-IMP1 source {relative}: {exc}")
            continue
        if imp1_source_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-MHG3-IMP1 source hash is stale for {relative}")

    imp1_document = "docs/fgc-hyp1-mhg-implicit.md"
    if imp1_record.get("derivation_document") != imp1_document:
        failures.append("FGC-1-HYP1-MHG3-IMP1 owner document has drifted")
    try:
        digest = sha256((REPOSITORY / imp1_document).read_bytes()).hexdigest()
    except OSError as exc:
        failures.append(f"cannot hash FGC-1-HYP1-MHG3-IMP1 owner document: {exc}")
    else:
        if imp1_record.get("derivation_document_sha256") != digest:
            failures.append("FGC-1-HYP1-MHG3-IMP1 owner-document hash is stale")

    imp1_certificate = imp1_record.get("implicit_acceleration_certificate", {})
    canonical_imp1_certificate = json.dumps(
        imp1_certificate, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    if imp1_record.get("implicit_acceleration_certificate_sha256") != sha256(
        canonical_imp1_certificate.encode("utf-8")
    ).hexdigest():
        failures.append("FGC-1-HYP1-MHG3-IMP1 certificate digest does not match its payload")

    imp1_formulation = imp1_certificate.get("formulation", {})
    if (
        imp1_formulation.get("physical_equations") != "unredefined_ACT1_VAR1"
        or imp1_formulation.get("gauge_equations")
        != "REF1_full_reference_connection_modified_harmonic"
        or imp1_formulation.get("acceleration_variables")
        != "base_metric_and_scalar_coordinate_time_second_derivatives"
        or imp1_formulation.get("differentiation_method")
        != "exact_first_tangent_forward_ad"
        or imp1_formulation.get("implicit_theorem")
        != "finite_dimensional_real_implicit_function_theorem"
        or imp1_formulation.get("coordinate_radius") != "4"
        or imp1_formulation.get("reference_radial_domain_minimum") != "1/2"
        or imp1_formulation.get("tilde_normal_factor") != "4"
        or imp1_formulation.get("hat_normal_factor") != "9"
    ):
        failures.append("FGC-1-HYP1-MHG3-IMP1 formulation contract has drifted")

    imp1_solution = imp1_certificate.get("reference_solution", {})
    expected_imp1_parameters = {
        "planck_mass": "2",
        "beta": "-1/4",
        "mu": "3",
        "g4": "1/2",
        "alpha": "0",
        "eta": "1/2",
    }
    def _all_exact_zero(value):
        if isinstance(value, list):
            return all(_all_exact_zero(item) for item in value)
        return value == "0"

    if (
        imp1_record.get("input_fixture_id") != "FGCQR_flat_reference_vacuum"
        or imp1_record.get("input_model_id") != "FGC-QR"
        or imp1_solution.get("fixture_id") != "FGCQR_flat_reference_vacuum"
        or imp1_solution.get("model_id") != "FGC-QR"
        or imp1_solution.get("action_parameters") != expected_imp1_parameters
        or imp1_solution.get("base_metric_determinant") != "-1"
        or imp1_solution.get("inverse_metric_g00") != "-1"
        or imp1_solution.get("effective_planck_coefficient") != "4"
        or not _all_exact_zero(imp1_solution.get("constraint_up"))
        or not _all_exact_zero(
            imp1_solution.get("covariant_constraint_derivative")
        )
        or not _all_exact_zero(
            imp1_solution.get("gauge_extension_residual_vector")
        )
        or not _all_exact_zero(imp1_solution.get("full_residual_vector"))
    ):
        failures.append("FGC-1-HYP1-MHG3-IMP1 exact FGC-QR root has drifted")

    imp1_smoothness = imp1_certificate.get("smoothness_domain_at_root", {})
    expected_imp1_smoothness = {
        "base_metric_determinant_negative",
        "coordinate_time_covector_physical_timelike",
        "areal_radius_positive",
        "effective_planck_coefficient_positive",
        "reference_radius_strictly_inside_annulus",
        "auxiliary_inverse_metrics_lorentzian",
    }
    if set(imp1_smoothness) != expected_imp1_smoothness or any(
        value is not True for value in imp1_smoothness.values()
    ):
        failures.append("FGC-1-HYP1-MHG3-IMP1 smoothness-domain root checks have drifted")

    imp1_acceleration = imp1_certificate.get("acceleration_map", {})
    expected_imp1_acceleration_order = [
        "h_tt.dtt",
        "h_tr.dtt",
        "h_rr.dtt",
        "areal_radius.dtt",
        "phi.dtt",
        "chi.dtt",
    ]
    expected_imp1_equation_order = [
        "metric_tt_mhg",
        "metric_tr_mhg",
        "metric_rr_mhg",
        "metric_theta_theta_mhg",
        "scalar_phi",
        "scalar_chi",
    ]
    expected_imp1_jacobian = [
        ["36", "0", "9", "9", "0", "0"],
        ["0", "72", "0", "0", "0", "0"],
        ["4", "0", "1", "-1", "0", "0"],
        ["64", "0", "-16", "0", "0", "0"],
        ["0", "0", "0", "0", "-1", "0"],
        ["0", "0", "0", "0", "0", "-1"],
    ]
    inverse = imp1_acceleration.get("jacobian_inverse")
    inverse_shape_exact = (
        isinstance(inverse, list)
        and len(inverse) == 6
        and all(isinstance(row, list) and len(row) == 6 for row in inverse)
    )
    if (
        imp1_acceleration.get("acceleration_order")
        != expected_imp1_acceleration_order
        or imp1_acceleration.get("equation_order") != expected_imp1_equation_order
        or imp1_acceleration.get("base_acceleration_vector") != ["0"] * 6
        or not _all_exact_zero(
            imp1_acceleration.get("base_full_residual_vector")
        )
        or imp1_acceleration.get("jacobian") != expected_imp1_jacobian
        or imp1_acceleration.get("mhg1_coordinate_time_kinetic_block")
        != expected_imp1_jacobian
        or imp1_acceleration.get("jacobian_determinant") != "-165888"
        or imp1_acceleration.get("jacobian_rank") != 6
        or not inverse_shape_exact
        or imp1_acceleration.get("seeded_primals_reproduce_base_residual") is not True
        or imp1_acceleration.get("left_inverse_identity_exact") is not True
        or imp1_acceleration.get("right_inverse_identity_exact") is not True
        or imp1_acceleration.get("differentiation_method")
        != "exact_first_tangent_forward_ad"
    ):
        failures.append("FGC-1-HYP1-MHG3-IMP1 acceleration map has drifted")

    expected_imp1_checks = {
        "fgcqr_branch_retained",
        "fgcqr_action_parameters_exact",
        "flat_base_metric_two_jet_exact",
        "areal_radius_matches_reference_chart_exact",
        "vacuum_scalar_two_jets_exact",
        "base_coordinate_time_accelerations_zero",
        "complete_ref1_residual_zero_at_root",
        "reference_gauge_surface_exact_at_root",
        "reference_gauge_extension_zero_at_root",
        "all_rational_smooth_domain_denominators_regular_at_root",
        "exact_forward_tangent_primals_reproduce_root",
        "full_residual_acceleration_jacobian_matches_mhg1_kinetic_block",
        "acceleration_jacobian_full_rank",
        "acceleration_jacobian_determinant_nonzero",
        "exact_two_sided_inverse_identities",
    }
    imp1_checks = imp1_certificate.get("verified_exact_checks", {})
    if set(imp1_checks) != expected_imp1_checks or any(
        value is not True for value in imp1_checks.values()
    ):
        failures.append("FGC-1-HYP1-MHG3-IMP1 exact-check ledger has drifted")
    if (
        imp1_certificate.get("all_declared_exact_checks_pass") is not True
        or imp1_certificate.get(
            "local_smooth_implicit_coordinate_time_acceleration_branch_proven"
        )
        is not True
        or imp1_certificate.get("implicit_branch_neighborhood")
        != "exists_but_not_quantified"
    ):
        failures.append("FGC-1-HYP1-MHG3-IMP1 local IFT conclusion has drifted")

    expected_imp1_gate_status = {
        "flat_fgcqr_local_implicit_acceleration_branch_passed": True,
        "quantified_implicit_branch_neighborhood_passed": False,
        "complete_gauge_propagation_passed": False,
        "complete_first_order_reduction_passed": False,
        "uniform_open_domain_hyperbolicity_gate_passed": False,
        "full_fgc1_hyp1_health_gate_passed": False,
        "evolution_authorized": False,
        "publication_local_defocusing_gate_passed": False,
    }
    if imp1_record.get("gate_status") != expected_imp1_gate_status:
        failures.append("FGC-1-HYP1-MHG3-IMP1 gate-status ledger has drifted")
    expected_imp1_nonclaims = {
        "quantified_implicit_branch_neighborhood_proven",
        "nonzero_activated_solution_on_branch_derived",
        "complete_lower_order_first_order_sources_derived",
        "gauge_constraint_lower_order_propagation_proven",
        "reduction_constraint_propagation_proven",
        "physical_initial_constraints_solved",
        "uniform_open_domain_symmetrizer_proven",
        "open_retained_eft_domain_proven",
        "boundary_or_regular_center_system_derived",
        "evolution_authorized",
        "collapse_solution_derived",
        "metric_null_affine_defocusing_derived",
        "singularity_resolution_derived",
    }
    imp1_nonclaims = imp1_record.get("nonclaims", {})
    if set(imp1_nonclaims) != expected_imp1_nonclaims or any(
        value is not False for value in imp1_nonclaims.values()
    ):
        failures.append("FGC-1-HYP1-MHG3-IMP1 nonclaim ledger has drifted")
    if imp1_certificate.get("nonclaims") != imp1_nonclaims:
        failures.append("FGC-1-HYP1-MHG3-IMP1 nested and top-level nonclaims differ")

    prop1_path = REPOSITORY / "results" / "fgc-1-hyp1-mhg-propagation.json"
    try:
        prop1_record = json.loads(prop1_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse FGC-1-HYP1-MHG4-PROP1 result: {exc}")
        return failures
    if (
        prop1_record.get("project_version") != "0.11.0"
        or prop1_record.get("schema_version") != 1
        or prop1_record.get("artifact_id") != "FGC-1-HYP1-MHG4-PROP1"
    ):
        failures.append("FGC-1-HYP1-MHG4-PROP1 identity, version, or schema has drifted")
    if prop1_record.get("classification") != (
        "exact_lower_order_reference_gauge_operator_conditional_on_noether_and_scalar_equations_not_direct_metric_residual_divergence_or_constraint_propagation"
    ):
        failures.append("FGC-1-HYP1-MHG4-PROP1 has promoted its restricted classification")
    if prop1_record.get("development_status") != (
        "conditional_lower_order_reference_gauge_operator_complete_metric_derived_constraint_propagation_first_order_domain_and_evolution_gates_open"
    ):
        failures.append("FGC-1-HYP1-MHG4-PROP1 development status has drifted")

    expected_prop1_configs = {
        "action": "configs/fgc/fgc-1-action-gate.toml",
        "variation": "configs/fgc/fgc-1-metric-variation.toml",
        "reduction": "configs/fgc/fgc-1-hyp1-reduction.toml",
        "symbol": "configs/fgc/fgc-1-hyp1-symbol.toml",
        "modified_harmonic": "configs/fgc/fgc-1-hyp1-modified-harmonic.toml",
        "reference": "configs/fgc/fgc-1-hyp1-mhg-reference.toml",
        "implicit": "configs/fgc/fgc-1-hyp1-mhg-implicit.toml",
        "propagation": "configs/fgc/fgc-1-hyp1-mhg-propagation.toml",
    }
    if prop1_record.get("source_configs") != expected_prop1_configs:
        failures.append("FGC-1-HYP1-MHG4-PROP1 source config map has drifted")
    prop1_config_hashes = prop1_record.get("source_config_sha256", {})
    if set(prop1_config_hashes) != set(expected_prop1_configs.values()):
        failures.append("FGC-1-HYP1-MHG4-PROP1 config hash ledger has unexpected keys")
    for relative in expected_prop1_configs.values():
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-MHG4-PROP1 config {relative}: {exc}")
            continue
        if prop1_config_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-MHG4-PROP1 config hash is stale for {relative}")

    expected_prop1_sources = {
        "src/recursive_horizons/fgc/exact_linear_algebra.py",
        "src/recursive_horizons/fgc/exact_tangent.py",
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/spherical_symbol.py",
        "src/recursive_horizons/fgc/modified_harmonic.py",
        "src/recursive_horizons/fgc/reference_connection.py",
        "src/recursive_horizons/fgc/modified_harmonic_reference.py",
        "src/recursive_horizons/fgc/modified_harmonic_implicit.py",
        "src/recursive_horizons/fgc/modified_harmonic_propagation.py",
        "src/recursive_horizons/fgc/__init__.py",
        "scripts/reproduce_fgc_action.py",
        "scripts/reproduce_fgc_metric_variation.py",
        "scripts/reproduce_fgc_hyp1_reduction.py",
        "scripts/reproduce_fgc_hyp1_symbol.py",
        "scripts/reproduce_fgc_hyp1_modified_harmonic.py",
        "scripts/reproduce_fgc_hyp1_mhg_reference.py",
        "scripts/reproduce_fgc_hyp1_mhg_implicit.py",
        "scripts/reproduce_fgc_hyp1_mhg_propagation.py",
    }
    prop1_source_hashes = prop1_record.get("implementation_sha256", {})
    if set(prop1_source_hashes) != expected_prop1_sources:
        failures.append("FGC-1-HYP1-MHG4-PROP1 implementation hash ledger has unexpected keys")
    for relative in expected_prop1_sources:
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-MHG4-PROP1 source {relative}: {exc}")
            continue
        if prop1_source_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-MHG4-PROP1 source hash is stale for {relative}")

    prop1_document = "docs/fgc-hyp1-mhg-propagation.md"
    if prop1_record.get("derivation_document") != prop1_document:
        failures.append("FGC-1-HYP1-MHG4-PROP1 owner document has drifted")
    try:
        digest = sha256((REPOSITORY / prop1_document).read_bytes()).hexdigest()
    except OSError as exc:
        failures.append(f"cannot hash FGC-1-HYP1-MHG4-PROP1 owner document: {exc}")
    else:
        if prop1_record.get("derivation_document_sha256") != digest:
            failures.append("FGC-1-HYP1-MHG4-PROP1 owner-document hash is stale")

    prop1_certificate = prop1_record.get("gauge_propagation_certificate", {})
    canonical_prop1_certificate = json.dumps(
        prop1_certificate, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    if prop1_record.get("gauge_propagation_certificate_sha256") != sha256(
        canonical_prop1_certificate.encode("utf-8")
    ).hexdigest():
        failures.append("FGC-1-HYP1-MHG4-PROP1 certificate digest does not match its payload")

    prop1_formulation = prop1_certificate.get("formulation", {})
    if (
        prop1_formulation.get("physical_equations") != "unredefined_ACT1_VAR1"
        or prop1_formulation.get("gauge_equations")
        != "REF1_full_reference_connection_modified_harmonic"
        or prop1_formulation.get("operator")
        != "nabla_mu[F*hat_P_alpha^(beta_mu_nu)*nabla_beta_C^alpha]"
        or prop1_formulation.get("operator_input")
        != "independent_spherically_symmetric_contravariant_gauge_vector_two_jet"
        or prop1_formulation.get("noether_identity")
        != "nabla_mu_E^mu_nu+E_phi*nabla^nu_phi+E_chi*nabla^nu_chi=0"
        or prop1_formulation.get("off_scalar_shell_source")
        != "E_phi*nabla^nu_phi+E_chi*nabla^nu_chi"
        or prop1_formulation.get("homogeneous_only_on_both_scalar_equations") is not True
        or prop1_formulation.get("direct_full_metric_residual_divergence_evaluated") is not False
        or prop1_formulation.get("reference_connection_enters_through_metric_derived_C_not_as_independent_forcing") is not True
        or prop1_formulation.get("gauge_vector_order") != ["C^t", "C^r"]
        or prop1_formulation.get("gauge_jet_order")
        != ["value", "dt", "dr", "dtt", "dtr", "drr"]
        or prop1_formulation.get("output_order")
        != ["nu=t", "nu=r", "nu=theta", "nu=phi"]
        or prop1_formulation.get("flat_coordinate_radius") != "4"
        or prop1_formulation.get("activated_coordinate_radius") != "9/2"
        or prop1_formulation.get("tilde_normal_factor") != "4"
        or prop1_formulation.get("hat_normal_factor") != "9"
    ):
        failures.append("FGC-1-HYP1-MHG4-PROP1 formulation contract has drifted")
    if (
        prop1_record.get("input_fixture_ids")
        != ["FGCQR_flat_reference_vacuum", "FGCQR_activated_generic"]
        or prop1_record.get("input_model_ids") != ["FGC-QR", "FGC-QR"]
    ):
        failures.append("FGC-1-HYP1-MHG4-PROP1 input fixture ledger has drifted")

    expected_prop1_checks = {
        "implicit_flat_fgcqr_root_predecessor_exact",
        "flat_zero_gauge_operator_exact",
        "flat_zero_two_routes_agree_exactly",
        "flat_nonzero_probe_exercises_operator",
        "flat_nonzero_two_routes_agree_exactly",
        "activated_two_routes_agree_exactly",
        "activated_gradient_F_contribution_nonzero",
        "activated_projector_derivative_contribution_nonzero",
        "activated_connection_contribution_nonzero",
        "activated_curvature_contribution_nonzero",
        "activated_first_gauge_jet_coefficients_nonzero",
        "activated_second_gauge_jet_coefficients_nonzero",
        "exact_linear_homogeneity",
        "exact_linear_superposition",
        "spherical_principal_columns_match_mhg1_exactly",
        "angular_covariant_effects_retained_without_independent_angular_C",
        "operator_input_explicitly_independent_not_metric_derived",
        "reference_connection_is_not_an_independent_operator_source",
    }
    prop1_checks = prop1_certificate.get("verified_exact_checks", {})
    if set(prop1_checks) != expected_prop1_checks or any(
        value is not True for value in prop1_checks.values()
    ):
        failures.append("FGC-1-HYP1-MHG4-PROP1 exact-check ledger has drifted")
    if (
        prop1_certificate.get("all_declared_exact_checks_pass") is not True
        or prop1_certificate.get(
            "conditional_lower_order_reference_gauge_operator_derived"
        )
        is not True
    ):
        failures.append("FGC-1-HYP1-MHG4-PROP1 conditional conclusion has drifted")

    prop1_flat = prop1_certificate.get("flat_root_control", {})
    prop1_activated = prop1_certificate.get("activated_operator_probe", {})
    operator_controls = (
        prop1_flat.get("zero_probe", {}),
        prop1_flat.get("nonzero_independent_probe", {}),
        prop1_activated.get("operator", {}),
        prop1_activated.get("zero_independent_probe_control", {}),
    )
    if any(
        control.get("coordinate_divergence")
        != control.get("expanded_divergence")
        or control.get("two_routes_agree_exactly") is not True
        or control.get("coordinate_route_uses_raw_partial_jet_not_covariant_second")
        is not True
        or control.get("expanded_route_uses_covariant_second") is not True
        for control in operator_controls
    ):
        failures.append("FGC-1-HYP1-MHG4-PROP1 direct and expanded routes have drifted")
    activated_operator = prop1_activated.get("operator", {})
    if (
        not _all_exact_zero(
            prop1_flat.get("zero_probe", {}).get("coordinate_divergence")
        )
        or not _all_exact_zero(
            prop1_activated.get("zero_independent_probe_control", {}).get(
                "coordinate_divergence"
            )
        )
        or _all_exact_zero(
            prop1_flat.get("nonzero_independent_probe", {}).get(
                "coordinate_divergence"
            )
        )
        or _all_exact_zero(activated_operator.get("coordinate_divergence"))
        or any(
            _all_exact_zero(activated_operator.get(name))
            for name in (
                "expanded_connection_wave_contribution",
                "expanded_curvature_commutator_term",
                "expanded_gradient_F_term",
                "expanded_projector_derivative_term",
            )
        )
    ):
        failures.append("FGC-1-HYP1-MHG4-PROP1 zero/nonzero operator controls have drifted")

    prop1_principal = prop1_certificate.get("principal_regression", {})
    coefficients = prop1_principal.get("C_value_first_second_partial_blocks", {})
    radial = prop1_principal.get("radial_principal_polynomial", [])
    mhg1_principal = prop1_principal.get("mhg1_full_four_index_principal", [])
    coefficient_shape_exact = (
        set(coefficients) == {"value", "dt", "dr", "dtt", "dtr", "drr"}
        and all(
            isinstance(block, list)
            and len(block) == 4
            and all(isinstance(row, list) and len(row) == 2 for row in block)
            for block in coefficients.values()
        )
    )
    spherical_principal_exact = (
        isinstance(radial, list)
        and isinstance(mhg1_principal, list)
        and len(radial) == len(mhg1_principal) == 4
        and all(
            radial[nu][alpha] == mhg1_principal[nu][alpha]
            for nu in range(4)
            for alpha in range(2)
        )
    )
    if (
        not coefficient_shape_exact
        or not spherical_principal_exact
        or prop1_principal.get("represented_C_input_order") != ["C^t", "C^r"]
        or prop1_principal.get("coefficient_block_shape") != [4, 2]
        or prop1_principal.get("second_partial_C_order")
        != ["dtt", "dtr", "drr"]
        or prop1_principal.get("spherical_t_r_columns_match_mhg1_exactly")
        is not True
        or prop1_principal.get("angular_independent_C_columns_not_represented")
        is not True
    ):
        failures.append("FGC-1-HYP1-MHG4-PROP1 principal regression has drifted")

    expected_prop1_gate_status = {
        "conditional_lower_order_reference_gauge_operator_passed": True,
        "direct_metric_derived_gauge_constraint_propagation_passed": False,
        "complete_gauge_propagation_passed": False,
        "complete_first_order_reduction_passed": False,
        "uniform_open_domain_hyperbolicity_gate_passed": False,
        "full_fgc1_hyp1_health_gate_passed": False,
        "evolution_authorized": False,
        "publication_local_defocusing_gate_passed": False,
    }
    if prop1_record.get("gate_status") != expected_prop1_gate_status:
        failures.append("FGC-1-HYP1-MHG4-PROP1 gate-status ledger has drifted")
    expected_prop1_nonclaims = {
        "direct_metric_derived_full_residual_divergence_evaluated",
        "metric_derived_gauge_constraint_propagation_proven",
        "constraint_preserving_initial_boundary_value_problem_proven",
        "complete_lower_order_first_order_sources_derived",
        "reduction_constraint_propagation_proven",
        "physical_initial_constraints_solved",
        "quantified_implicit_branch_neighborhood_proven",
        "uniform_open_domain_symmetrizer_proven",
        "open_retained_eft_domain_proven",
        "boundary_or_regular_center_system_derived",
        "evolution_authorized",
        "collapse_solution_derived",
        "metric_null_affine_defocusing_derived",
        "singularity_resolution_derived",
    }
    prop1_nonclaims = prop1_record.get("nonclaims", {})
    if set(prop1_nonclaims) != expected_prop1_nonclaims or any(
        value is not False for value in prop1_nonclaims.values()
    ):
        failures.append("FGC-1-HYP1-MHG4-PROP1 nonclaim ledger has drifted")
    if prop1_certificate.get("nonclaims") != prop1_nonclaims:
        failures.append("FGC-1-HYP1-MHG4-PROP1 nested and top-level nonclaims differ")

    fo1_path = REPOSITORY / "results" / "fgc-1-hyp1-fo1-rc1.json"
    try:
        fo1_record = json.loads(fo1_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse FGC-1-HYP1-FO1-RC1 result: {exc}")
        return failures
    expected_fo1_classification = (
        "exact_local_first_order_differential_algebraic_lift_with_flat_root_complete_linearization_and_kinematic_radial_reduction_constraint_identity_not_nonlinear_evolution"
    )
    if (
        fo1_record.get("project_version") != "0.11.0"
        or fo1_record.get("schema_version") != 1
        or fo1_record.get("artifact_id") != "FGC-1-HYP1-FO1-RC1"
    ):
        failures.append("FGC-1-HYP1-FO1-RC1 identity, version, or schema has drifted")
    if fo1_record.get("classification") != expected_fo1_classification:
        failures.append("FGC-1-HYP1-FO1-RC1 has promoted its restricted classification")
    if fo1_record.get("development_status") != (
        "local_first_order_dae_lift_flat_root_linearization_and_kinematic_radial_constraint_identity_complete_nonlinear_constraint_domain_and_evolution_gates_open"
    ):
        failures.append("FGC-1-HYP1-FO1-RC1 development status has drifted")

    expected_fo1_configs = {
        "action": "configs/fgc/fgc-1-action-gate.toml",
        "variation": "configs/fgc/fgc-1-metric-variation.toml",
        "reduction": "configs/fgc/fgc-1-hyp1-reduction.toml",
        "symbol": "configs/fgc/fgc-1-hyp1-symbol.toml",
        "modified_harmonic": "configs/fgc/fgc-1-hyp1-modified-harmonic.toml",
        "reference": "configs/fgc/fgc-1-hyp1-mhg-reference.toml",
        "implicit": "configs/fgc/fgc-1-hyp1-mhg-implicit.toml",
        "propagation": "configs/fgc/fgc-1-hyp1-mhg-propagation.toml",
        "first_order": "configs/fgc/fgc-1-hyp1-fo1-rc1.toml",
    }
    if fo1_record.get("source_configs") != expected_fo1_configs:
        failures.append("FGC-1-HYP1-FO1-RC1 source config map has drifted")
    fo1_config_hashes = fo1_record.get("source_config_sha256", {})
    if set(fo1_config_hashes) != set(expected_fo1_configs.values()):
        failures.append("FGC-1-HYP1-FO1-RC1 config hash ledger has unexpected keys")
    for relative in expected_fo1_configs.values():
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-FO1-RC1 config {relative}: {exc}")
            continue
        if fo1_config_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-FO1-RC1 config hash is stale for {relative}")

    expected_fo1_sources = {
        "src/recursive_horizons/fgc/exact_linear_algebra.py",
        "src/recursive_horizons/fgc/exact_tangent.py",
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/spherical_symbol.py",
        "src/recursive_horizons/fgc/modified_harmonic.py",
        "src/recursive_horizons/fgc/reference_connection.py",
        "src/recursive_horizons/fgc/modified_harmonic_reference.py",
        "src/recursive_horizons/fgc/modified_harmonic_implicit.py",
        "src/recursive_horizons/fgc/modified_harmonic_propagation.py",
        "src/recursive_horizons/fgc/modified_harmonic_first_order.py",
        "src/recursive_horizons/fgc/__init__.py",
        "scripts/reproduce_fgc_action.py",
        "scripts/reproduce_fgc_metric_variation.py",
        "scripts/reproduce_fgc_hyp1_reduction.py",
        "scripts/reproduce_fgc_hyp1_symbol.py",
        "scripts/reproduce_fgc_hyp1_modified_harmonic.py",
        "scripts/reproduce_fgc_hyp1_mhg_reference.py",
        "scripts/reproduce_fgc_hyp1_mhg_implicit.py",
        "scripts/reproduce_fgc_hyp1_mhg_propagation.py",
        "scripts/reproduce_fgc_hyp1_fo1_rc1.py",
    }
    fo1_source_hashes = fo1_record.get("implementation_sha256", {})
    if set(fo1_source_hashes) != expected_fo1_sources:
        failures.append("FGC-1-HYP1-FO1-RC1 implementation hash ledger has unexpected keys")
    for relative in expected_fo1_sources:
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-FO1-RC1 source {relative}: {exc}")
            continue
        if fo1_source_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-FO1-RC1 source hash is stale for {relative}")

    fo1_document = "docs/fgc-hyp1-fo1-rc1.md"
    if fo1_record.get("derivation_document") != fo1_document:
        failures.append("FGC-1-HYP1-FO1-RC1 owner document has drifted")
    try:
        digest = sha256((REPOSITORY / fo1_document).read_bytes()).hexdigest()
    except OSError as exc:
        failures.append(f"cannot hash FGC-1-HYP1-FO1-RC1 owner document: {exc}")
    else:
        if fo1_record.get("derivation_document_sha256") != digest:
            failures.append("FGC-1-HYP1-FO1-RC1 owner-document hash is stale")

    fo1_certificate = fo1_record.get("first_order_certificate", {})
    canonical_fo1_certificate = json.dumps(
        fo1_certificate, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    if fo1_record.get("first_order_certificate_sha256") != sha256(
        canonical_fo1_certificate.encode("utf-8")
    ).hexdigest():
        failures.append("FGC-1-HYP1-FO1-RC1 certificate digest does not match its payload")

    expected_fields = [
        "h_tt",
        "h_tr",
        "h_rr",
        "areal_radius",
        "phi",
        "chi",
    ]
    expected_state_order = [
        f"{group}.{field}"
        for group in ("u", "p", "q")
        for field in expected_fields
    ]
    expected_equation_order = (
        [f"dt_u_minus_p.{field}" for field in expected_fields]
        + [
            "metric_tt_mhg",
            "metric_tr_mhg",
            "metric_rr_mhg",
            "metric_theta_theta_mhg",
            "scalar_phi",
            "scalar_chi",
        ]
        + [f"dt_q_minus_dr_p.{field}" for field in expected_fields]
    )
    expected_argument_order = [
        f"{group}.{field}"
        for group in ("u", "p", "q", "p_t", "p_r", "q_r")
        for field in expected_fields
    ]
    fo1_formulation = fo1_certificate.get("formulation", {})
    if (
        fo1_formulation.get("physical_equations") != "unredefined_ACT1_VAR1"
        or fo1_formulation.get("gauge_equations")
        != "REF1_full_reference_connection_modified_harmonic"
        or fo1_formulation.get("first_order_kind")
        != "exact_local_first_order_differential_algebraic_lift"
        or fo1_formulation.get("state_definition")
        != "U=(u,p=partial_t_u,q=partial_r_u)"
        or fo1_formulation.get("state_order") != expected_state_order
        or fo1_formulation.get("equation_order") != expected_equation_order
        or fo1_formulation.get("physical_argument_order") != expected_argument_order
        or fo1_formulation.get("radial_reduction_constraint")
        != "C_r=q-partial_r_u"
        or fo1_formulation.get("reference_radial_domain_minimum") != "1/2"
        or fo1_formulation.get("tilde_normal_factor") != "4"
        or fo1_formulation.get("hat_normal_factor") != "9"
    ):
        failures.append("FGC-1-HYP1-FO1-RC1 formulation contract has drifted")
    if (
        fo1_record.get("input_fixture_ids")
        != ["FGCQR_flat_reference_vacuum", "FGCQR_activated_generic"]
        or fo1_record.get("input_model_ids") != ["FGC-QR", "FGC-QR"]
    ):
        failures.append("FGC-1-HYP1-FO1-RC1 input fixture ledger has drifted")

    expected_fo1_checks = {
        "frozen_eighteen_variable_and_equation_orders_exact",
        "flat_fgcqr_action_parameters_match_imp1_exact",
        "flat_first_order_two_jet_roundtrip_exact",
        "activated_first_order_two_jet_roundtrip_exact",
        "flat_complete_ref1_residual_reproduced_exact",
        "activated_complete_ref1_residual_reproduced_exact",
        "flat_all_kinematic_rows_and_reduction_constraints_zero_exact",
        "activated_all_kinematic_rows_and_reduction_constraints_zero_exact",
        "kinematic_radial_reduction_constraint_time_derivative_zero_exact",
        "flat_complete_ref1_residual_zero_exact",
        "flat_all_argument_tangent_primals_reproduce_residual",
        "activated_all_argument_tangent_primals_reproduce_residual",
        "flat_weighted_tangent_superposition_exact",
        "activated_weighted_tangent_superposition_exact",
        "flat_acceleration_block_matches_imp1_exact",
        "flat_acceleration_block_has_imp1_determinant_and_rank",
        "flat_implicit_branch_derivative_identity_exact",
        "flat_linearized_p_q_principal_matrix_matches_mhg1_exact",
        "activated_control_is_off_shell_and_nontrivial",
        "no_nonlinear_acceleration_solve_performed",
    }
    fo1_checks = fo1_certificate.get("verified_exact_checks", {})
    if set(fo1_checks) != expected_fo1_checks or any(
        value is not True for value in fo1_checks.values()
    ):
        failures.append("FGC-1-HYP1-FO1-RC1 exact-check ledger has drifted")
    if (
        fo1_certificate.get("all_declared_exact_checks_pass") is not True
        or fo1_certificate.get("exact_local_first_order_dae_lift_derived") is not True
        or fo1_certificate.get(
            "complete_flat_root_linearized_first_order_map_derived"
        )
        is not True
        or fo1_certificate.get(
            "kinematic_radial_reduction_constraint_identity_derived"
        )
        is not True
    ):
        failures.append("FGC-1-HYP1-FO1-RC1 exact conclusions have drifted")

    fo1_flat = fo1_certificate.get("controls", {}).get("flat", {})
    fo1_activated = fo1_certificate.get("controls", {}).get(
        "activated_off_shell", {}
    )
    flat_jacobian = fo1_flat.get("argument_jacobian", {}).get("jacobian", [])
    if (
        not _all_exact_zero(fo1_flat.get("complete_ref1_residual"))
        or fo1_activated.get("solution_status") != "off_shell_control_not_a_solution"
        or _all_exact_zero(fo1_activated.get("complete_ref1_residual"))
        or not isinstance(flat_jacobian, list)
        or len(flat_jacobian) != 6
        or any(not isinstance(row, list) or len(row) != 36 for row in flat_jacobian)
        or fo1_activated.get("acceleration_block_rank") != 6
    ):
        failures.append("FGC-1-HYP1-FO1-RC1 flat/activated controls have drifted")

    fo1_linearized = fo1_certificate.get(
        "flat_root_linearized_first_order_system", {}
    )
    branch_derivative = fo1_linearized.get("implicit_branch_derivative", [])
    principal_18 = fo1_linearized.get("linearized_radial_principal_matrix", [])
    source_18 = fo1_linearized.get("linearized_lower_order_source_matrix", [])
    principal_12 = fo1_linearized.get("p_q_radial_principal_matrix", [])
    mhg1_12 = fo1_linearized.get(
        "mhg1_standard_first_order_principal_matrix", []
    )
    if (
        fo1_linearized.get("acceleration_jacobian_determinant") != "-165888"
        or fo1_linearized.get("acceleration_jacobian_rank") != 6
        or fo1_linearized.get("implicit_branch_derivative_identity_exact") is not True
        or fo1_linearized.get("p_q_principal_matrix_matches_mhg1_exact") is not True
        or principal_12 != mhg1_12
        or len(branch_derivative) != 6
        or any(len(row) != 30 for row in branch_derivative)
        or len(principal_18) != 18
        or any(len(row) != 18 for row in principal_18)
        or len(source_18) != 18
        or any(len(row) != 18 for row in source_18)
    ):
        failures.append("FGC-1-HYP1-FO1-RC1 flat-root linearized system has drifted")

    expected_fo1_gate_status = {
        "exact_local_first_order_dae_lift_passed": True,
        "complete_flat_root_linearized_first_order_map_passed": True,
        "kinematic_radial_reduction_constraint_identity_passed": True,
        "explicit_nonlinear_acceleration_map_passed": False,
        "complete_metric_gauge_physical_constraint_propagation_passed": False,
        "uniform_open_domain_hyperbolicity_gate_passed": False,
        "full_fgc1_hyp1_health_gate_passed": False,
        "evolution_authorized": False,
        "publication_local_defocusing_gate_passed": False,
    }
    if fo1_record.get("gate_status") != expected_fo1_gate_status:
        failures.append("FGC-1-HYP1-FO1-RC1 gate-status ledger has drifted")
    expected_fo1_nonclaims = {
        "explicit_nonlinear_acceleration_map_implemented",
        "solved_acceleration_map_away_from_flat_root_derived",
        "quantified_implicit_branch_neighborhood_proven",
        "complete_nonlinear_lower_order_first_order_sources_derived",
        "complete_quasilinear_first_order_evolution_system_derived",
        "metric_derived_gauge_constraint_propagation_proven",
        "complete_reduction_constraint_system_propagation_proven",
        "physical_hamiltonian_momentum_constraints_derived",
        "physical_initial_constraints_solved",
        "constraint_preserving_initial_boundary_value_problem_proven",
        "uniform_open_domain_symmetrizer_proven",
        "open_retained_eft_domain_proven",
        "boundary_or_regular_center_system_derived",
        "evolution_authorized",
        "collapse_solution_derived",
        "metric_null_affine_defocusing_derived",
        "singularity_resolution_derived",
    }
    fo1_nonclaims = fo1_record.get("nonclaims", {})
    if set(fo1_nonclaims) != expected_fo1_nonclaims or any(
        value is not False for value in fo1_nonclaims.values()
    ):
        failures.append("FGC-1-HYP1-FO1-RC1 nonclaim ledger has drifted")
    if fo1_certificate.get("nonclaims") != fo1_nonclaims:
        failures.append("FGC-1-HYP1-FO1-RC1 nested and top-level nonclaims differ")

    qift1_path = REPOSITORY / "results" / "fgc-1-hyp1-dom1-qift1.json"
    try:
        qift1_source = qift1_path.read_text(encoding="utf-8")
        qift1_record = json.loads(
            qift1_source, object_pairs_hook=_reject_duplicate_json_pairs
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot parse FGC-1-HYP1-DOM1-QIFT1 result: {exc}")
        return failures
    if not isinstance(qift1_record, dict):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 result must be a JSON object")
        return failures
    canonical_qift1_source = (
        json.dumps(qift1_record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    if qift1_source != canonical_qift1_source:
        failures.append("FGC-1-HYP1-DOM1-QIFT1 result JSON is not canonical")
    expected_qift1_classification = (
        "exact_rational_interval_uniform_contraction_certificate_for_a_"
        "full_dimensional_local_ref1_implicit_acceleration_branch_not_a_"
        "pde_health_or_eft_theorem"
    )
    if (
        qift1_record.get("project_version") != "0.11.0"
        or qift1_record.get("schema_version") != 1
        or qift1_record.get("artifact_id") != "FGC-1-HYP1-DOM1-QIFT1"
    ):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 identity, version, or schema has drifted")
    if qift1_record.get("classification") != expected_qift1_classification:
        failures.append("FGC-1-HYP1-DOM1-QIFT1 has promoted its restricted classification")
    if qift1_record.get("development_status") != (
        "full_dimensional_local_implicit_branch_box_complete_constraint_"
        "hyperbolicity_eft_and_evolution_gates_open"
    ):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 development status has drifted")

    expected_qift1_configs = {
        "action": "configs/fgc/fgc-1-action-gate.toml",
        "variation": "configs/fgc/fgc-1-metric-variation.toml",
        "reduction": "configs/fgc/fgc-1-hyp1-reduction.toml",
        "symbol": "configs/fgc/fgc-1-hyp1-symbol.toml",
        "modified_harmonic": "configs/fgc/fgc-1-hyp1-modified-harmonic.toml",
        "reference": "configs/fgc/fgc-1-hyp1-mhg-reference.toml",
        "implicit": "configs/fgc/fgc-1-hyp1-mhg-implicit.toml",
        "propagation": "configs/fgc/fgc-1-hyp1-mhg-propagation.toml",
        "first_order": "configs/fgc/fgc-1-hyp1-fo1-rc1.toml",
        "quantified_domain": "configs/fgc/fgc-1-hyp1-dom1-qift1.toml",
    }
    if qift1_record.get("source_configs") != expected_qift1_configs:
        failures.append("FGC-1-HYP1-DOM1-QIFT1 source config map has drifted")
    qift1_config_hashes = qift1_record.get("source_config_sha256", {})
    if set(qift1_config_hashes) != set(expected_qift1_configs.values()):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 config hash ledger has unexpected keys")
    for relative in expected_qift1_configs.values():
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-DOM1-QIFT1 config {relative}: {exc}")
            continue
        if qift1_config_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-DOM1-QIFT1 config hash is stale for {relative}")

    expected_qift1_sources = {
        "src/recursive_horizons/fgc/exact_linear_algebra.py",
        "src/recursive_horizons/fgc/exact_tangent.py",
        "src/recursive_horizons/fgc/exact_interval.py",
        "src/recursive_horizons/fgc/interval_tangent.py",
        "src/recursive_horizons/fgc/spherical_reduction.py",
        "src/recursive_horizons/fgc/spherical_symbol.py",
        "src/recursive_horizons/fgc/modified_harmonic.py",
        "src/recursive_horizons/fgc/reference_connection.py",
        "src/recursive_horizons/fgc/modified_harmonic_reference.py",
        "src/recursive_horizons/fgc/modified_harmonic_implicit.py",
        "src/recursive_horizons/fgc/modified_harmonic_propagation.py",
        "src/recursive_horizons/fgc/modified_harmonic_first_order.py",
        "src/recursive_horizons/fgc/modified_harmonic_quantified_domain.py",
        "src/recursive_horizons/fgc/__init__.py",
        "scripts/reproduce_fgc_action.py",
        "scripts/reproduce_fgc_metric_variation.py",
        "scripts/reproduce_fgc_hyp1_reduction.py",
        "scripts/reproduce_fgc_hyp1_symbol.py",
        "scripts/reproduce_fgc_hyp1_modified_harmonic.py",
        "scripts/reproduce_fgc_hyp1_mhg_reference.py",
        "scripts/reproduce_fgc_hyp1_mhg_implicit.py",
        "scripts/reproduce_fgc_hyp1_mhg_propagation.py",
        "scripts/reproduce_fgc_hyp1_fo1_rc1.py",
        "scripts/reproduce_fgc_hyp1_dom1_qift1.py",
    }
    qift1_source_hashes = qift1_record.get("implementation_sha256", {})
    if set(qift1_source_hashes) != expected_qift1_sources:
        failures.append("FGC-1-HYP1-DOM1-QIFT1 implementation hash ledger has unexpected keys")
    for relative in expected_qift1_sources:
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-DOM1-QIFT1 source {relative}: {exc}")
            continue
        if qift1_source_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-DOM1-QIFT1 source hash is stale for {relative}")

    qift1_document = "docs/fgc-hyp1-dom1-qift1.md"
    if qift1_record.get("derivation_document") != qift1_document:
        failures.append("FGC-1-HYP1-DOM1-QIFT1 owner document has drifted")
    try:
        digest = sha256((REPOSITORY / qift1_document).read_bytes()).hexdigest()
    except OSError as exc:
        failures.append(f"cannot hash FGC-1-HYP1-DOM1-QIFT1 owner document: {exc}")
    else:
        if qift1_record.get("derivation_document_sha256") != digest:
            failures.append("FGC-1-HYP1-DOM1-QIFT1 owner-document hash is stale")

    qift1_certificate = qift1_record.get("quantified_branch_certificate", {})
    canonical_qift1_certificate = json.dumps(
        qift1_certificate, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    if qift1_record.get("quantified_branch_certificate_sha256") != sha256(
        canonical_qift1_certificate.encode("utf-8")
    ).hexdigest():
        failures.append("FGC-1-HYP1-DOM1-QIFT1 certificate digest does not match its payload")

    expected_qift1_parameter_order = [
        f"{group}.{field}"
        for group in ("u", "p", "q", "p_r", "q_r")
        for field in ("h_tt", "h_tr", "h_rr", "areal_radius", "phi", "chi")
    ]
    expected_qift1_acceleration_order = [
        f"p_t.{field}"
        for field in ("h_tt", "h_tr", "h_rr", "areal_radius", "phi", "chi")
    ]
    qift1_formulation = qift1_certificate.get("formulation", {})
    if (
        qift1_formulation.get("residual") != "complete_REF1_residual_R(a;z)=0"
        or qift1_formulation.get("parameter_definition")
        != "z=(u,p,q,p_r,q_r)_in_frozen_base_field_order"
        or qift1_formulation.get("acceleration_definition")
        != "a=p_t_in_frozen_base_field_order"
        or qift1_formulation.get("parameter_order") != expected_qift1_parameter_order
        or qift1_formulation.get("acceleration_order") != expected_qift1_acceleration_order
        or qift1_formulation.get("parameter_dimension") != 30
        or qift1_formulation.get("acceleration_dimension") != 6
        or qift1_formulation.get("coordinate_radius_is_fixed_not_a_box_axis") is not True
        or qift1_formulation.get("coordinate_radius") != "4"
        or qift1_formulation.get("tilde_normal_factor") != "4"
        or qift1_formulation.get("hat_normal_factor") != "9"
    ):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 formulation contract has drifted")
    if (
        qift1_record.get("input_fixture_ids") != ["FGCQR_flat_reference_vacuum"]
        or qift1_record.get("input_model_ids") != ["FGC-QR"]
    ):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 input fixture ledger has drifted")

    parameter_box = qift1_certificate.get("parameter_box", {})
    acceleration_box = qift1_certificate.get("acceleration_box", {})
    if (
        parameter_box.get("dimension") != 30
        or parameter_box.get("uniform_half_width") != "1/65536"
        or parameter_box.get("all_axes_nonzero_width") is not True
        or parameter_box.get("closed_box_has_nonempty_open_interior_in_local_jet_space")
        is not True
        or len(parameter_box.get("center", [])) != 30
        or acceleration_box.get("dimension") != 6
        or acceleration_box.get("uniform_half_width") != "1/128"
        or acceleration_box.get("center") != ["0"] * 6
    ):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 declared boxes have drifted")

    def _parse_interval(value):
        if not isinstance(value, dict) or set(value) != {"lower", "upper"}:
            return None
        try:
            lower = Fraction(value["lower"])
            upper = Fraction(value["upper"])
        except (TypeError, ValueError, ZeroDivisionError):
            return None
        return (lower, upper) if lower <= upper else None

    interval_data = qift1_certificate.get("interval_residual_and_jacobian", {})
    residual_box = interval_data.get("residual_at_center_acceleration_parameter_box", [])
    jacobian_box = interval_data.get("acceleration_jacobian_box", [])
    if (
        interval_data.get("parameter_dimension") != 30
        or interval_data.get("acceleration_dimension") != 6
        or interval_data.get("parameter_half_width") != "1/65536"
        or interval_data.get("acceleration_half_width") != "1/128"
        or interval_data.get("all_parameter_axes_have_nonzero_width") is not True
        or len(residual_box) != 6
        or any(_parse_interval(value) is None for value in residual_box)
        or len(jacobian_box) != 6
        or any(
            not isinstance(row, list)
            or len(row) != 6
            or any(_parse_interval(value) is None for value in row)
            for row in jacobian_box
        )
    ):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 interval residual/Jacobian shape has drifted")
    else:
        if any(not (parsed[0] <= 0 <= parsed[1]) for parsed in map(_parse_interval, residual_box)):
            failures.append("FGC-1-HYP1-DOM1-QIFT1 parameter residual box excludes the exact root")

    qift1_krawczyk = qift1_certificate.get("krawczyk_contraction_certificate", {})
    try:
        contraction_bound = Fraction(
            qift1_krawczyk.get("contraction_infinity_norm_upper_bound", "nan")
        )
        inverse_bound = Fraction(
            qift1_krawczyk.get(
                "uniform_acceleration_jacobian_inverse_infinity_norm_upper_bound",
                "nan",
            )
        )
        minimum_inclusion_margin = Fraction(
            qift1_krawczyk.get("minimum_strict_interior_inclusion_margin", "nan")
        )
    except (ValueError, ZeroDivisionError):
        contraction_bound = Q(2)
        inverse_bound = Q(0)
        minimum_inclusion_margin = Q(0)
    image_box = qift1_krawczyk.get("krawczyk_image_box", [])
    image_intervals = [_parse_interval(value) for value in image_box]
    if (
        qift1_krawczyk.get("norm")
        != "unscaled_componentwise_infinity_norm_in_frozen_field_order"
        or qift1_krawczyk.get("theorem_route")
        != "uniform_banach_contraction_with_krawczyk_box_inclusion"
        or not Q(0) <= contraction_bound < Q(1, 250)
        or inverse_bound <= 0
        or minimum_inclusion_margin <= Q(3, 512)
        or qift1_krawczyk.get("contraction_strictly_less_than_one") is not True
        or qift1_krawczyk.get("krawczyk_image_strictly_inside_acceleration_box")
        is not True
        or qift1_krawczyk.get("unique_acceleration_root_for_every_declared_parameter_point")
        is not True
        or qift1_krawczyk.get(
            "fixed_point_iterations_converge_from_every_point_in_declared_acceleration_box"
        )
        is not True
        or len(image_intervals) != 6
        or any(value is None for value in image_intervals)
        or any(
            not (-Q(1, 512) < value[0] <= value[1] < Q(1, 512))
            for value in image_intervals
            if value is not None
        )
    ):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 contraction/inclusion certificate has drifted")

    regular_margins = qift1_certificate.get("regular_domain_ledger", {}).get(
        "strict_positive_margins", {}
    )
    expected_margin_bounds = {
        "areal_radius_positive_lower_margin": Q(3),
        "base_metric_determinant_negative_lower_margin": Q(99, 100),
        "coordinate_time_inverse_metric_negative_lower_margin": Q(99, 100),
        "effective_planck_coefficient_positive_lower_margin": Q(3),
        "tilde_auxiliary_determinant_negative_lower_margin": Q(1, 65),
        "hat_auxiliary_determinant_negative_lower_margin": Q(1, 29),
    }
    if set(regular_margins) != set(expected_margin_bounds) | {
        "reference_annulus_coordinate_margin"
    }:
        failures.append("FGC-1-HYP1-DOM1-QIFT1 regular-margin ledger has unexpected keys")
    else:
        for name, bound in expected_margin_bounds.items():
            try:
                value = Fraction(regular_margins[name])
            except (TypeError, ValueError, ZeroDivisionError):
                value = Q(0)
            if value <= bound:
                failures.append(f"FGC-1-HYP1-DOM1-QIFT1 margin {name} lost its strict bound")
        if regular_margins.get("reference_annulus_coordinate_margin") != "7/2":
            failures.append("FGC-1-HYP1-DOM1-QIFT1 reference-annulus margin has drifted")

    expected_qift1_checks = {
        "frozen_parameter_and_acceleration_orders_exact",
        "all_thirty_parameter_axes_have_nonzero_width",
        "flat_center_complete_ref1_residual_zero_exact",
        "flat_center_acceleration_jacobian_matches_imp1_exact",
        "interval_acceleration_jacobian_contains_exact_center",
        "all_regular_chart_margins_strictly_positive",
        "contraction_infinity_norm_strictly_below_one",
        "krawczyk_image_strictly_inside_acceleration_box",
        "uniform_unique_acceleration_root_for_every_parameter_point",
    }
    qift1_checks = qift1_certificate.get("verified_exact_checks", {})
    if set(qift1_checks) != expected_qift1_checks or any(
        value is not True for value in qift1_checks.values()
    ):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 exact-check ledger has drifted")
    if (
        qift1_certificate.get("all_declared_exact_checks_pass") is not True
        or qift1_certificate.get(
            "quantified_full_dimensional_local_implicit_branch_box_derived"
        )
        is not True
    ):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 exact conclusion has drifted")

    expected_qift1_gate_status = {
        "quantified_full_dimensional_local_implicit_branch_box_passed": True,
        "complete_metric_gauge_physical_constraint_propagation_passed": False,
        "uniform_radial_strong_hyperbolicity_box_passed": False,
        "retained_eft_validity_envelope_passed": False,
        "full_fgc1_hyp1_health_gate_passed": False,
        "evolution_authorized": False,
        "publication_local_defocusing_gate_passed": False,
    }
    if qift1_record.get("gate_status") != expected_qift1_gate_status:
        failures.append("FGC-1-HYP1-DOM1-QIFT1 gate-status ledger has drifted")
    expected_qift1_nonclaims = {
        "exact_closed_form_nonlinear_acceleration_map_derived",
        "full_quasilinear_first_order_evolution_system_derived",
        "metric_derived_gauge_constraint_propagation_proven",
        "complete_reduction_constraint_system_propagation_proven",
        "physical_hamiltonian_momentum_constraints_derived",
        "physical_initial_constraints_solved",
        "constraint_preserving_initial_boundary_value_problem_proven",
        "uniform_open_domain_symmetrizer_proven",
        "uniform_radial_strong_hyperbolicity_box_proven",
        "retained_eft_cutoff_and_omitted_operator_envelope_defined",
        "open_retained_eft_domain_proven",
        "boundary_or_regular_center_system_derived",
        "evolution_authorized",
        "collapse_solution_derived",
        "metric_null_affine_defocusing_derived",
        "singularity_resolution_derived",
    }
    qift1_nonclaims = qift1_record.get("nonclaims", {})
    if set(qift1_nonclaims) != expected_qift1_nonclaims or any(
        value is not False for value in qift1_nonclaims.values()
    ):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 nonclaim ledger has drifted")
    if qift1_certificate.get("nonclaims") != qift1_nonclaims:
        failures.append("FGC-1-HYP1-DOM1-QIFT1 nested and top-level nonclaims differ")

    def _contains_float(value):
        if isinstance(value, float):
            return True
        if isinstance(value, dict):
            return any(_contains_float(item) for item in value.values())
        if isinstance(value, list):
            return any(_contains_float(item) for item in value)
        return False

    if _contains_float(qift1_certificate):
        failures.append("FGC-1-HYP1-DOM1-QIFT1 certificate contains floating-point data")

    comp1_path = REPOSITORY / "results" / "fgc-1-hyp1-con1-comp1.json"
    try:
        comp1_source = comp1_path.read_text(encoding="utf-8")
        comp1_record = json.loads(
            comp1_source, object_pairs_hook=_reject_duplicate_json_pairs
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot parse FGC-1-HYP1-CON1-COMP1 result: {exc}")
        return failures
    if not isinstance(comp1_record, dict):
        failures.append("FGC-1-HYP1-CON1-COMP1 result must be a JSON object")
        return failures
    canonical_comp1_source = (
        json.dumps(comp1_record, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    if comp1_source != canonical_comp1_source:
        failures.append("FGC-1-HYP1-CON1-COMP1 result JSON is not canonical")

    expected_comp1_classification = (
        "exact_local_activated_metric_defined_compatibility_witness_"
        "conditional_on_external_qift1_acceleration_root"
    )
    if (
        comp1_record.get("project_version") != "0.11.0"
        or comp1_record.get("schema_version") != 1
        or comp1_record.get("artifact_id") != "FGC-1-HYP1-CON1-COMP1"
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 identity, version, or schema has drifted")
    if comp1_record.get("classification") != expected_comp1_classification:
        failures.append("FGC-1-HYP1-CON1-COMP1 has promoted its restricted classification")
    if comp1_record.get("development_status") != (
        "exact_local_compatible_constraint_datum_complete_propagation_initial_"
        "data_hyperbolicity_eft_and_evolution_gates_open"
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 development status has drifted")

    expected_comp1_configs = dict(expected_qift1_configs)
    expected_comp1_configs["compatible_data"] = (
        "configs/fgc/fgc-1-hyp1-con1-comp1.toml"
    )
    if comp1_record.get("source_configs") != expected_comp1_configs:
        failures.append("FGC-1-HYP1-CON1-COMP1 source config map has drifted")
    comp1_config_hashes = comp1_record.get("source_config_sha256", {})
    if set(comp1_config_hashes) != set(expected_comp1_configs.values()):
        failures.append("FGC-1-HYP1-CON1-COMP1 config hash ledger has unexpected keys")
    for relative in expected_comp1_configs.values():
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-CON1-COMP1 config {relative}: {exc}")
            continue
        if comp1_config_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-CON1-COMP1 config hash is stale for {relative}")

    expected_comp1_predecessors = {
        "results/fgc-1-action-gate.json",
        "results/fgc-1-metric-variation.json",
        "results/fgc-1-hyp1-reduction.json",
        "results/fgc-1-hyp1-symbol.json",
        "results/fgc-1-hyp1-modified-harmonic.json",
        "results/fgc-1-hyp1-mhg-reference.json",
        "results/fgc-1-hyp1-mhg-implicit.json",
        "results/fgc-1-hyp1-mhg-propagation.json",
        "results/fgc-1-hyp1-fo1-rc1.json",
        "results/fgc-1-hyp1-dom1-qift1.json",
    }
    comp1_predecessor_hashes = comp1_record.get("source_results_sha256", {})
    if set(comp1_predecessor_hashes) != expected_comp1_predecessors:
        failures.append("FGC-1-HYP1-CON1-COMP1 predecessor-result ledger has unexpected keys")
    for relative in expected_comp1_predecessors:
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-CON1-COMP1 predecessor {relative}: {exc}")
            continue
        if comp1_predecessor_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-CON1-COMP1 predecessor hash is stale for {relative}")

    expected_comp1_sources = expected_qift1_sources | {
        "src/recursive_horizons/fgc/modified_harmonic_constraints.py",
        "scripts/reproduce_fgc_hyp1_con1_comp1.py",
    }
    comp1_source_hashes = comp1_record.get("implementation_sha256", {})
    if set(comp1_source_hashes) != expected_comp1_sources:
        failures.append("FGC-1-HYP1-CON1-COMP1 implementation hash ledger has unexpected keys")
    for relative in expected_comp1_sources:
        try:
            digest = sha256((REPOSITORY / relative).read_bytes()).hexdigest()
        except OSError as exc:
            failures.append(f"cannot hash FGC-1-HYP1-CON1-COMP1 source {relative}: {exc}")
            continue
        if comp1_source_hashes.get(relative) != digest:
            failures.append(f"FGC-1-HYP1-CON1-COMP1 source hash is stale for {relative}")

    comp1_document = "docs/fgc-hyp1-con1-comp1.md"
    if comp1_record.get("derivation_document") != comp1_document:
        failures.append("FGC-1-HYP1-CON1-COMP1 owner document has drifted")
    try:
        digest = sha256((REPOSITORY / comp1_document).read_bytes()).hexdigest()
    except OSError as exc:
        failures.append(f"cannot hash FGC-1-HYP1-CON1-COMP1 owner document: {exc}")
    else:
        if comp1_record.get("derivation_document_sha256") != digest:
            failures.append("FGC-1-HYP1-CON1-COMP1 owner-document hash is stale")

    comp1_certificate = comp1_record.get("compatible_constraint_certificate", {})
    canonical_comp1_certificate = json.dumps(
        comp1_certificate, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    if comp1_record.get("compatible_constraint_certificate_sha256") != sha256(
        canonical_comp1_certificate.encode("utf-8")
    ).hexdigest():
        failures.append("FGC-1-HYP1-CON1-COMP1 certificate digest does not match its payload")
    if (
        comp1_certificate.get("artifact_id") != "FGC-1-HYP1-CON1-COMP1"
        or comp1_certificate.get("classification") != expected_comp1_classification
        or comp1_certificate.get("all_declared_exact_checks_pass") is not True
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 certificate identity or conclusion has drifted")

    expected_comp1_acceleration_order = [
        f"p_t.{field}"
        for field in ("h_tt", "h_tr", "h_rr", "areal_radius", "phi", "chi")
    ]
    expected_comp1_projection_order = [
        "H=E_tt-2sE_tr+s^2E_rr",
        "M=E_tr-sE_rr",
    ]
    comp1_formulation = comp1_record.get("formulation", {})
    if (
        comp1_formulation.get("physical_equations") != "unredefined_ACT1_VAR1"
        or comp1_formulation.get("gauge_constraint")
        != "C^a=-tilde_g^(bc)(Gamma^a_bc-bar_Gamma^a_bc)"
        or comp1_formulation.get("acceleration_order")
        != expected_comp1_acceleration_order
        or comp1_formulation.get("projection_order")
        != expected_comp1_projection_order
        or comp1_formulation.get("normal_gauge_extension_rows")
        != ["metric_tt_mhg", "metric_tr_mhg"]
        or comp1_formulation.get("normal_gauge_derivative_order")
        != ["nabla_t_C^t", "nabla_t_C^r"]
        or comp1_formulation.get("tilde_normal_factor") != 4
        or comp1_formulation.get("hat_normal_factor") != 9
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 formulation contract has drifted")
    expected_comp1_proof_contract = {
        "require_metric_defined_C_zero_exact",
        "require_all_non_normal_nabla_C_zero_exact",
        "require_fo1_reduction_rows_zero_exact",
        "require_normal_shift_zero_exact",
        "require_unredefined_Ett_and_Etr_zero_for_every_acceleration",
        "require_ref1_scalar_equations_unmodified_exact",
        "require_exact_normal_gauge_to_extension_map_invertible",
        "require_full_mhg_root_implies_normal_nabla_C_zero",
        "require_gauge_extension_zero_at_full_mhg_root",
        "require_unredefined_full_equations_at_local_root",
        "require_no_floating_point_or_sampled_proof",
    }
    comp1_proof_contract = comp1_record.get("proof_contract", {})
    if set(comp1_proof_contract) != expected_comp1_proof_contract or any(
        value is not True for value in comp1_proof_contract.values()
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 proof contract has drifted")

    comp1_datum_contract = comp1_record.get("activated_local_parameter_datum", {})
    if comp1_datum_contract != {
        "fixture_id": "FGCQR_qift1_activated_local_parameter",
        "model_id": "FGC-QR",
        "phi_value": "1/131072",
        "phi_radial_derivative": "1/131072",
        "derived_second_derivatives": ["h_tt.drr", "areal_radius.drr"],
        "datum_status": "exact_local_parameter_datum_not_an_initial_slice_or_solution_family",
    }:
        failures.append("FGC-1-HYP1-CON1-COMP1 activated-datum contract has drifted")
    comp1_datum = comp1_certificate.get("activated_datum", {})
    comp1_state = comp1_datum.get("state", {})
    if (
        comp1_datum.get("activated_phi") != "1/131072"
        or comp1_datum.get("activated_phi_radial_derivative") != "1/131072"
        or comp1_datum.get("h_tt_drr")
        != "-584115552257/18889465931341141901312"
        or comp1_datum.get("areal_radius_drr")
        != "-584115552257/4722366482835285475328"
        or comp1_datum.get("strictly_inside_qift_parameter_box") is not True
        or comp1_datum.get("qift_acceleration_root_is_external_predecessor_fact")
        is not True
        or comp1_datum.get("qift_acceleration_root_solved_here") is not False
        or comp1_state.get("branch") != "FGC-QR"
        or comp1_state.get("phi", {}).get("value") != "1/131072"
        or comp1_state.get("phi", {}).get("dr") != "1/131072"
        or comp1_state.get("h_tt", {}).get("drr")
        != "-584115552257/18889465931341141901312"
        or comp1_state.get("areal_radius", {}).get("drr")
        != "-584115552257/4722366482835285475328"
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 exact activated datum has drifted")

    comp1_projection = comp1_certificate.get("physical_projections", {})
    comp1_polynomial = comp1_certificate.get("acceleration_polynomial", {})
    comp1_coefficients = comp1_polynomial.get("coefficients", [])
    expected_comp1_indices = {()} | {(i,) for i in range(6)} | {
        (i, j) for i in range(6) for j in range(i, 6)
    }
    coefficient_indices = {
        tuple(item.get("indices", []))
        for item in comp1_coefficients
        if isinstance(item, dict)
    }
    if (
        comp1_projection.get("H") != "0"
        or comp1_projection.get("M") != "0"
        or comp1_projection.get("source_is_unredefined_residual") is not True
        or comp1_projection.get("projection_order") != expected_comp1_projection_order
        or comp1_polynomial.get("acceleration_order")
        != expected_comp1_acceleration_order
        or comp1_polynomial.get("projection_order")
        != expected_comp1_projection_order
        or comp1_polynomial.get("maximum_total_degree") != 2
        or comp1_polynomial.get("coefficient_count_per_projection") != 28
        or comp1_polynomial.get("evaluation_count") != 28
        or len(comp1_coefficients) != 28
        or coefficient_indices != expected_comp1_indices
        or any(
            not isinstance(item, dict)
            or item.get("H") != "0"
            or item.get("M") != "0"
            or item.get("power") != len(item.get("indices", []))
            for item in comp1_coefficients
        )
        or comp1_polynomial.get("all_56_coefficients_zero") is not True
        or comp1_polynomial.get("all_H_coefficients_zero") is not True
        or comp1_polynomial.get("all_M_coefficients_zero") is not True
        or any(
            value is not True
            for value in comp1_polynomial.get("degree_bound_guard", {}).values()
        )
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 acceleration-polynomial proof has drifted")

    comp1_gauge = comp1_certificate.get("gauge", {})
    comp1_kinematic = comp1_certificate.get("kinematic", {})
    comp1_ref1_composition = comp1_certificate.get("ref1_composition", {})
    zero_six = ["0"] * 6
    if (
        comp1_gauge.get("constraint_up") != ["0"] * 4
        or comp1_gauge.get("constraint_down") != ["0"] * 4
        or comp1_kinematic.get("u_time_definition_residual") != zero_six
        or comp1_kinematic.get("radial_reduction_constraint") != zero_six
        or comp1_kinematic.get("mixed_partial_compatibility_residual") != zero_six
        or comp1_kinematic.get(
            "radial_reduction_constraint_time_derivative_on_definition_shell"
        )
        != zero_six
        or comp1_ref1_composition.get("equation_order")
        != [
            "metric_tt_mhg",
            "metric_tr_mhg",
            "metric_rr_mhg",
            "metric_theta_theta_mhg",
            "scalar_phi",
            "scalar_chi",
        ]
        or comp1_ref1_composition.get("scalar_equations_unmodified") is not True
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 gauge or kinematic exact rows have drifted")

    expected_comp1_checks = {
        "activated_parameter_point_strictly_inside_qift1_z_box",
        "activated_principal_control_uses_exact_compatible_state",
        "activated_principal_metric_and_regulator_discriminants_positive",
        "activated_principal_physical_and_regulator_cones_share_no_root",
        "all_56_unredefined_constraint_polynomial_coefficients_zero",
        "conditional_qift_root_implies_normal_gauge_compatibility",
        "metric_defined_C_zero",
        "metric_defined_all_non_normal_nabla_C_components_zero",
        "nabla_r_Ct_and_Cr_zero",
        "normal_shift_zero",
        "normal_gauge_extension_map_diagonal_nonzero",
        "normal_gauge_extension_map_nonsingular",
        "physical_projection_source_is_unredefined",
        "qift_root_is_external_predecessor_not_solved_here",
        "reduction_rows_independent_of_accelerations_by_repacking",
        "reduction_rows_zero",
        "ref1_scalar_equations_unmodified",
        "unredefined_Ett_and_Etr_zero_for_every_acceleration",
        "unredefined_Etr_zero",
        "unredefined_H_zero",
        "unredefined_M_zero",
    }
    comp1_checks = comp1_certificate.get("local_exact_checks", {})
    if set(comp1_checks) != expected_comp1_checks or any(
        value is not True for value in comp1_checks.values()
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 local exact-check ledger has drifted")

    comp1_normal_map = comp1_certificate.get("normal_gauge_extension_map", {})
    if (
        comp1_normal_map.get("matrix")
        != [
            ["-2473901162487/137438953472", "0"],
            ["0", "2473901162487/137438953472"],
        ]
        or comp1_normal_map.get("determinant")
        != "-6120186961754529976025169/18889465931478580854784"
        or comp1_normal_map.get("determinant_nonzero") is not True
        or comp1_normal_map.get("diagonal_entries_nonzero") is not True
        or comp1_normal_map.get("off_diagonal_entries_zero") is not True
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 normal gauge-extension map has drifted")

    comp1_cones = comp1_certificate.get("activated_principal_cone_control", {})
    if (
        comp1_cones.get("classification")
        != "exact_pointwise_activated_principal_cone_control_not_an_open_domain_hyperbolicity_certificate"
        or comp1_cones.get("fixture_id")
        != "CON1_COMP1_activated_compatible_principal_control"
        or comp1_cones.get("compatible_state_pullback_exact") is not True
        or comp1_cones.get("metric_null_discriminant") != "4"
        or comp1_cones.get("regulator_factor")
        != [
            "-151115727450454257303553/151115727449629623582728",
            "0",
            "1298074214621900990924908803391487/1298074214614817441201214220926976",
        ]
        or comp1_cones.get("physical_vs_regulator_resultant")
        != "302231454904756805304321/1684996666647875129859515217356021951040479227080081906456724504576"
        or comp1_cones.get("metric_null_discriminant_positive") is not True
        or comp1_cones.get("regulator_discriminant_positive") is not True
        or comp1_cones.get("physical_vs_regulator_resultant_nonzero") is not True
        or comp1_cones.get("cones_share_no_root") is not True
        or comp1_cones.get("regulator_cone_proportional_to_metric") is not False
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 pointwise cone control has drifted")

    expected_comp1_premises = {
        "datum_strictly_inside_qift1_parameter_box",
        "all_non_normal_nabla_C_components_zero_exact",
        "fo1_reduction_rows_zero_exact",
        "metric_defined_C_zero_exact",
        "normal_gauge_to_extension_map_invertible_exact",
        "normal_shift_zero_exact",
        "qift1_unique_acceleration_root_available",
        "qift_root_is_external_predecessor_not_solved_here",
        "ref1_scalar_equations_unmodified_exact",
        "unredefined_Ett_and_Etr_zero_for_every_acceleration",
        "unredefined_H_w_w_and_M_w_r_acceleration_independent_and_zero",
    }
    comp1_theorem = comp1_record.get("theorem_composition", {})
    comp1_premises = comp1_theorem.get("premises", {})
    if (
        set(comp1_premises) != expected_comp1_premises
        or any(value is not True for value in comp1_premises.values())
        or comp1_theorem.get("all_machine_checked_premises_pass") is not True
        or comp1_theorem.get("full_mhg_root_implies_normal_nabla_C_zero") is not True
        or comp1_theorem.get("gauge_extension_zero_at_full_mhg_root") is not True
        or comp1_theorem.get("unredefined_full_equations_at_local_root") is not True
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 theorem composition has drifted")

    expected_comp1_gate_status = {
        "exact_activated_local_compatible_constraint_datum_passed": True,
        "metric_derived_gauge_constraint_propagation_passed": False,
        "complete_reduction_constraint_system_propagation_passed": False,
        "physical_initial_constraints_solved": False,
        "uniform_radial_strong_hyperbolicity_box_passed": False,
        "retained_eft_validity_envelope_passed": False,
        "full_fgc1_hyp1_health_gate_passed": False,
        "evolution_authorized": False,
        "publication_local_defocusing_gate_passed": False,
    }
    if comp1_record.get("gate_status") != expected_comp1_gate_status:
        failures.append("FGC-1-HYP1-CON1-COMP1 gate-status ledger has drifted")
    expected_comp1_nonclaims = {
        "qift_acceleration_root_solved_here",
        "exact_closed_form_nonlinear_acceleration_map_derived",
        "initial_slice_or_constraint_manifold_derived",
        "metric_derived_gauge_constraint_propagation_proven",
        "complete_reduction_constraint_system_propagation_proven",
        "physical_hamiltonian_momentum_constraint_propagation_proven",
        "physical_initial_constraints_solved",
        "nontrivial_compatible_initial_data_family_derived",
        "constraint_preserving_initial_boundary_value_problem_proven",
        "complete_nonlinear_lower_order_first_order_sources_derived",
        "uniform_open_domain_symmetrizer_proven",
        "uniform_radial_strong_hyperbolicity_box_proven",
        "retained_eft_cutoff_and_omitted_operator_envelope_defined",
        "open_retained_eft_domain_proven",
        "boundary_or_regular_center_system_derived",
        "evolution_authorized",
        "collapse_solution_derived",
        "metric_null_affine_defocusing_derived",
        "singularity_resolution_derived",
    }
    comp1_nonclaims = comp1_record.get("nonclaims", {})
    if set(comp1_nonclaims) != expected_comp1_nonclaims or any(
        value is not False for value in comp1_nonclaims.values()
    ):
        failures.append("FGC-1-HYP1-CON1-COMP1 nonclaim ledger has drifted")
    if comp1_certificate.get("nonclaims") != comp1_nonclaims:
        failures.append("FGC-1-HYP1-CON1-COMP1 nested and top-level nonclaims differ")
    if _contains_float(comp1_certificate):
        failures.append("FGC-1-HYP1-CON1-COMP1 certificate contains floating-point data")

    nested_path = REPOSITORY / "results" / "nested-gradient-spectrum.json"
    try:
        nested_record = json.loads(nested_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        failures.append(f"cannot parse NGS-0 result: {exc}")
        return failures
    if (
        nested_record.get("project_version") != "0.11.0"
        or nested_record.get("schema_version") != 1
    ):
        failures.append("NGS-0 result has stale version or schema")
    if nested_record.get("model_id") != "NGS-0":
        failures.append("NGS-0 result is missing its model identifier")
    if nested_record.get("artifact") != (
        "nested_gradient_spectral_recurrence_probability_controls"
    ):
        failures.append("NGS-0 result is missing its bounded artifact label")
    if nested_record.get("classification") != (
        "conditional_mathematical_controls_not_nested_cosmology"
    ):
        failures.append("NGS-0 result is missing its restricted classification")

    band = nested_record.get("finite_band_selection", {})
    for key, expected in (
        ("k_star_squared", 1.0),
        ("k_star", 1.0),
        ("ell_star", 1.0),
        ("stationary_derivative_residual", 0.0),
        ("energy_coefficient_at_k_star", -1.0),
    ):
        value = band.get(key)
        if not isinstance(value, (int, float)) or not isfinite(value) or abs(value - expected) > 1e-12:
            failures.append(f"NGS-0 finite-band {key} is outside tolerance")
    if band.get("nonzero_finite_scale_selected") is not True:
        failures.append("NGS-0 finite-band selector does not retain its nonzero scale")

    recurrences = nested_record.get("recurrence_controls", {})
    for fixture, expected, golden_expected in (
        ("golden_fixture", 0.5 * (1.0 + sqrt(5.0)), True),
        ("non_golden_fixture", 1.0 + sqrt(2.0), False),
    ):
        recurrence = recurrences.get(fixture, {})
        ratio = recurrence.get("positive_ratio")
        residual = recurrence.get("scaled_characteristic_residual")
        if not isinstance(ratio, (int, float)) or not isfinite(ratio) or abs(ratio - expected) > 1e-12:
            failures.append(f"NGS-0 {fixture} ratio is outside tolerance")
        if not isinstance(residual, (int, float)) or not isfinite(residual) or abs(residual) > 1e-12:
            failures.append(f"NGS-0 {fixture} scaled characteristic residual is outside tolerance")
        if recurrence.get("golden_ratio_fixture") is not golden_expected:
            failures.append(f"NGS-0 {fixture} loses its golden/non-golden label")
        if recurrence.get("coefficients_derived_from_physics") is not False:
            failures.append(f"NGS-0 {fixture} must not claim derived recurrence coefficients")

    opportunity = nested_record.get("eternal_opportunity_control", {})
    probability = opportunity.get("probability_at_least_one")
    if (
        not isinstance(probability, (int, float))
        or not isfinite(probability)
        or abs(probability - 0.632120742768355) > 1e-12
    ):
        failures.append("NGS-0 finite occurrence probability is outside tolerance")
    normalization = opportunity.get("normalization_residual")
    if not isinstance(normalization, (int, float)) or not isfinite(normalization) or abs(normalization) > 1e-12:
        failures.append("NGS-0 occurrence normalization is outside tolerance")
    for key in ("fixed_probability_assumed", "independence_assumed"):
        if opportunity.get(key) is not True:
            failures.append(f"NGS-0 occurrence control loses {key}=true")
    for key in ("ergodicity_derived", "logical_necessity_proven"):
        if opportunity.get(key) is not False:
            failures.append(f"NGS-0 occurrence control must retain {key}=false")

    for key in (
        "covariant_fgc_action_derived",
        "nested_spacetime_constructed",
        "finite_cross_domain_spectrum_derived",
        "dark_matter_identified",
        "dark_energy_identified",
        "universal_hertz_increment_derived",
        "golden_ratio_cosmology_derived",
        "observers_guaranteed",
        "metaphysical_plan_proven",
    ):
        if nested_record.get("nonclaims", {}).get(key) is not False:
            failures.append(f"NGS-0 must explicitly retain {key}=false")

    post_comp1_records = (
        (
            "results/fgc-1-eft0-led1.json",
            "FGC-1-EFT0-LED1",
            ("conditional_declared_local_component_controls_passed",),
            ("retained_eft_validity_envelope_passed", "initial_data_or_evolution_authorized"),
        ),
        (
            "results/fgc-1-hyp1-dom2-mode1.json",
            "FGC-1-HYP1-DOM2-MODE1",
            ("exact_pointwise_ref1_branch_mode_controls_passed",),
            ("uniform_radial_strong_hyperbolicity_box_passed", "evolution_authorized"),
        ),
        (
            "results/fgc-1-hyp1-dom3-uhyp1.json",
            "FGC-1-HYP1-DOM3-UHYP1",
            ("compact_nonflat_radial_strong_hyperbolicity_on_implicit_branch_graph_passed",),
            ("multidirectional_strong_hyperbolicity_proven", "evolution_authorized"),
        ),
        (
            "results/fgc-1-hyp1-con2-mprop1.json",
            "FGC-1-HYP1-CON2-MPROP1",
            ("local_differential_identity", "complete_kinematic_1plus1_reduction_subsidiary"),
            (),
        ),
        (
            "results/fgc-1-hyp1-con3-cau1.json",
            "FGC-1-HYP1-CON3-CAU1",
            ("conditional_boundary_free_gauge_Cauchy_uniqueness_theorem_derived",),
            ("constraint_preserving_ACT1_IBVP_proven", "evolution_authorized"),
        ),
        (
            "results/fgc-1-hyp1-bnd1-md1.json",
            "FGC-1-HYP1-BND1-MD1",
            ("uniform_frozen_radial_main_system_boundary_dissipation_passed",),
            ("constraint_preserving_ACT1_IBVP_proven", "evolution_authorized"),
        ),
        (
            "results/fgc-1-def0-obs1.json",
            "FGC-1-DEF0-OBS1",
            ("exact_metric_null_observable_contract_and_controls_passed",),
            ("metric_null_affine_defocusing_DEF1_derived", "evolution_authorized"),
        ),
        (
            "results/fgc-1-eft1-open1.json",
            "FGC-1-EFT1-OPEN1",
            ("audit_definition_and_fail_closed_composition_passed",),
            ("retained_EFT_open_run_envelope_passed", "evolution_authorized", "collapse_authorized", "DEF1_authorized"),
        ),
        (
            "results/fgc-1-src1-nl1.json",
            "FGC-1-SRC1-NL1",
            ("nonlinear_REF1_source_and_branch_solver_verified",),
            (),
        ),
        (
            "results/fgc-1-con4-phy1.json",
            "FGC-1-CON4-PHY1",
            ("physical_gauge_reduction_constraint_system_closed",),
            (),
        ),
        (
            "results/fgc-1-ctr1-reg1.json",
            "FGC-1-CTR1-REG1",
            ("regular_center_formulation_verified",),
            (),
        ),
        (
            "results/fgc-1-run1-sym1.json",
            "FGC-1-RUN1-SYM1",
            ("authorization_contract_and_protocol_validation_passed",),
            (
                "classical_spherical_diagnostic_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "FGCQR_holdout_execution_authorized",
            ),
        ),
        (
            "results/fgc-1-id0-pref1.json",
            "FGC-1-ID0-PREF1",
            ("PROTO1_initial_data_preflight_obstruction_verified",),
            (
                "PROTO1_initial_data_calibration_authorized",
                "PROTO1_FGCQR_holdout_execution_authorized",
                "classical_spherical_diagnostic_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
            ),
        ),
        (
            "results/fgc-1-id1-fam1.json",
            "FGC-1-ID1-FAM1",
            ("nonzero_width_finite_mass_constraint_compatible_family_constructed",),
            (),
        ),
        (
            "results/fgc-1-dom4-run1.json",
            "FGC-1-DOM4-RUN1",
            ("nonzero_classical_spherical_run_envelope_passed",),
            (),
        ),
        (
            "results/fgc-1-hyp2-md1.json",
            "FGC-1-HYP2-MD1",
            ("quantitative_all_covector_weak_coupling_health_envelope_passed",),
            (),
        ),
        (
            "results/fgc-1-bnd2-cp1.json",
            "FGC-1-BND2-CP1",
            ("spherical_boundary_or_domain_of_dependence_control_passed",),
            (),
        ),
        (
            "results/fgc-1-hlt1-mon1.json",
            "FGC-1-HLT1-MON1",
            ("classical_health_and_typed_stop_monitoring_verified",),
            (),
        ),
        (
            "results/fgc-1-num1-val1.json",
            "FGC-1-NUM1-VAL1",
            ("independent_solver_validation_and_measurement_contract_passed",),
            (),
        ),
        (
            "results/fgc-1-cal0-pref2.json",
            "FGC-1-CAL0-PREF2",
            ("PROTO3_pre_holdout_numerical_contract_obstruction_verified",),
            (
                "PROTO3_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
            ),
        ),
        (
            "results/fgc-1-pro4-frz1.json",
            "FGC-1-PRO4-FRZ1",
            (
                "PROTO4_outcome_neutral_protocol_frozen",
                "PROTO4_immutable_lineage_verified",
                "PROTO4_branch_specific_stop_partition_verified",
                "PROTO4_convergent_admission_contract_frozen",
            ),
            (
                "PROTO4_resolved_holdout_manifest_authorized",
                "HLT2_successor_monitor_implemented",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
            ),
        ),
        (
            "results/fgc-1-hlt2-mon2.json",
            "FGC-1-HLT2-MON2",
            (
                "PROTO4_successor_admission_monitor_implemented",
                "PROTO4_branch_applicability_executable",
                "PROTO4_weighted_spectral_admission_executable",
                "PROTO4_common_event_constraint_admission_executable",
                "declared_compact_initial_profile_weighted_spectral_control_passed",
                "FGCQR_candidate_health_definition_bound",
            ),
            (
                "SGBL_candidate_health_definition_bound",
                "constraint_to_DEF1_observable_stability_map_supplied",
                "fresh_GR0_calibration_authorized",
                "SGBL_comparison_execution_authorized",
                "FGCQR_holdout_execution_authorized",
                "classical_spherical_diagnostic_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
            ),
        ),
        (
            "results/fgc-1-id2-all1.json",
            "FGC-1-ID2-ALL1",
            (
                "all_amplitude_all_case_static_initial_admission_completed",
                "all_case_static_ledger_completed",
                "all_seven_GR0_candidates_statically_adjudicated",
                "all_thirty_five_FGCQR_slices_statically_adjudicated",
                "at_least_one_dynamic_calibration_candidate_remains",
                "fresh_GR0_dynamic_calibration_candidate_set_frozen",
            ),
            (
                "fresh_GR0_dynamic_calibration_authorized",
                "fresh_GR0_dynamic_calibration_completed",
                "SGBL_comparison_execution_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "classical_spherical_diagnostic_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
            ),
        ),
        (
            "results/fgc-1-cal1-pref3.json",
            "FGC-1-CAL1-PREF3",
            (
                "PROTO4_original_ID2_states_fail_literal_common_event_admission",
                "PROTO4_projected_states_still_fail_literal_common_event_admission",
                "PROTO4_semidiscrete_common_event_contract_obstructed",
                "PROTO5_premise_revision_required",
                "prospective_semidiscrete_repair_passes_every_t0_composition",
            ),
            (
                "PROTO4_fresh_GR0_dynamic_calibration_authorized",
                "fresh_GR0_dynamic_calibration_completed",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "classical_spherical_diagnostic_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
            ),
        ),
        (
            "results/fgc-1-pro5-frz1.json",
            "FGC-1-PRO5-FRZ1",
            (
                "PROTO5_outcome_neutral_protocol_frozen",
                "PROTO5_immutable_lineage_verified",
                "PROTO5_native_semidiscrete_initialization_contract_frozen",
                "PROTO5_method_owned_constraint_admission_frozen",
                "PROTO5_roundoff_zero_classification_contract_frozen",
            ),
            (
                "PROTO5_successor_runtime_compositor_implemented",
                "PROTO5_fresh_GR0_dynamic_calibration_authorized",
                "PROTO5_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "FGCQR_mechanism_rejected",
                "general_gradient_route_rejected",
            ),
        ),
        (
            "results/fgc-1-hlt3-mon3.json",
            "FGC-1-HLT3-MON3",
            (
                "PROTO5_successor_runtime_compositor_implemented",
                "PROTO5_fresh_GR0_dynamic_calibration_authorized",
            ),
            (
                "PROTO5_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "FGCQR_mechanism_rejected",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-cal2-pref4.json",
            "FGC-1-CAL2-PREF4",
            (
                "PROTO5_GR0_source_stop_ownership_contract_obstructed",
                "PROTO6_premise_revision_required",
                "both_amplitudes_reached_only_t0_common_event",
                "both_amplitudes_share_threshold_adjacent_source_stop",
                "last_accepted_fine_state_passes_original_source_gate",
                "one_half_proposal_passes_unchanged_source_gate",
                "one_eighth_proposal_passes_unchanged_source_gate",
                "unaccepted_full_proposal_misses_original_source_gate",
            ),
            (
                "fresh_GR0_dynamic_calibration_completed",
                "PROTO5_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "kinetic_condition_stop_reached",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-pro6-frz1.json",
            "FGC-1-PRO6-FRZ1",
            (
                "PROTO6_outcome_neutral_protocol_frozen",
                "PROTO6_immutable_lineage_verified",
                "PROTO6_accepted_state_source_terminal_contract_frozen",
                "PROTO6_internal_trial_source_retry_contract_frozen",
                "PROTO6_non_source_gates_inherited_unchanged",
            ),
            (
                "PROTO6_successor_runtime_compositor_implemented",
                "PROTO6_fresh_GR0_dynamic_calibration_authorized",
                "PROTO6_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-hlt4-mon4.json",
            "FGC-1-HLT4-MON4",
            (
                "PROTO6_successor_runtime_compositor_implemented",
                "PROTO6_fresh_GR0_dynamic_calibration_authorized",
            ),
            (
                "PROTO6_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-cal3-pref5.json",
            "FGC-1-CAL3-PREF5",
            (
                "PROTO6_unaccepted_candidate_endpoint_stop_ownership_contract_obstructed",
                "PROTO7_general_transaction_revision_required",
            ),
            (
                "fresh_GR0_dynamic_calibration_completed",
                "PROTO6_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-pro7-frz1.json",
            "FGC-1-PRO7-FRZ1",
            (
                "PROTO7_outcome_neutral_protocol_frozen",
                "PROTO7_immutable_lineage_verified",
                "PROTO7_accepted_state_source_terminal_contract_frozen",
                "PROTO7_unaccepted_proposal_source_retry_contract_frozen",
                "PROTO7_non_source_retry_veto_contract_frozen",
                "PROTO7_complete_failure_observability_contract_frozen",
                "PROTO7_non_source_gates_inherited_unchanged",
            ),
            (
                "PROTO7_successor_runtime_compositor_implemented",
                "PROTO7_fresh_GR0_dynamic_calibration_authorized",
                "PROTO7_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-hlt5-mon5.json",
            "FGC-1-HLT5-MON5",
            (
                "PROTO7_successor_runtime_compositor_implemented",
                "PROTO7_fresh_GR0_dynamic_calibration_authorized",
            ),
            (
                "PROTO7_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-cal4-pref6.json",
            "FGC-1-CAL4-PREF6",
            (
                "PROTO7_common_event_ownership_contract_obstructed",
                "PROTO8_common_event_admission_revision_required",
            ),
            (
                "fresh_GR0_dynamic_calibration_completed",
                "FGCQR_holdout_execution_authorized",
                "classical_spherical_diagnostic_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "FGCQR_mechanism_rejected",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-pro8-frz1.json",
            "FGC-1-PRO8-FRZ1",
            (
                "PROTO8_outcome_neutral_protocol_frozen",
                "PROTO8_immutable_lineage_verified",
                "PROTO8_evolution_owned_constraint_contract_frozen",
                "PROTO8_accumulated_roundoff_order_contract_frozen",
                "PROTO8_finest_pair_nested_spectral_contract_frozen",
                "PROTO8_physical_and_numerical_thresholds_inherited_unchanged",
            ),
            (
                "PROTO8_successor_runtime_compositor_implemented",
                "PROTO8_fresh_GR0_dynamic_calibration_authorized",
                "PROTO8_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-hlt6-mon6.json",
            "FGC-1-HLT6-MON6",
            (
                "PROTO8_successor_runtime_compositor_implemented",
                "PROTO8_fresh_GR0_dynamic_calibration_authorized",
            ),
            (
                "PROTO8_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-cal5-pref7.json",
            "FGC-1-CAL5-PREF7",
            (
                "PROTO8_campaign_terminated_normally",
                "PROTO8_common_event_ownership_repair_exercised",
                "PROTO9_resolution_ladder_revision_required",
            ),
            (
                "PROTO8_GR0_case_eligible",
                "PROTO9_frozen",
                "fresh_GR0_dynamic_calibration_completed",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-pro9-frz1.json",
            "FGC-1-PRO9-FRZ1",
            (),
            (
                "PROTO9_successor_runtime_implemented",
                "PROTO9_fresh_GR0_dynamic_calibration_authorized",
                "PROTO9_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-hlt7-mon7.json",
            "FGC-1-HLT7-MON7",
            ("PROTO9_successor_runtime_implemented",),
            (
                "PROTO9_fresh_GR0_dynamic_calibration_authorized",
                "PROTO9_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-cal6-pref8.json",
            "FGC-1-CAL6-PREF8",
            (
                "PROTO9_preflight_conditioning_diagnosis_completed",
                "PROTO10_spectral_contract_revision_required",
            ),
            (
                "PROTO9_ratio_only_failure_is_candidate_or_gradient_evidence",
                "PROTO10_frozen",
                "PROTO10_successor_runtime_implemented",
                "PROTO10_fresh_GR0_dynamic_calibration_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-pro10-frz1.json",
            "FGC-1-PRO10-FRZ1",
            (),
            (
                "PROTO10_successor_runtime_implemented",
                "PROTO10_fresh_GR0_dynamic_calibration_authorized",
                "PROTO10_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "FGCQR_mechanism_rejected",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-hlt8-mon8.json",
            "FGC-1-HLT8-MON8",
            (
                "PROTO10_successor_runtime_implemented",
                "PROTO10_fresh_GR0_dynamic_calibration_authorized",
            ),
            (
                "PROTO10_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "FGCQR_mechanism_rejected",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-cal7-pref9.json",
            "FGC-1-CAL7-PREF9",
            (
                "PROTO10_campaign_terminated_normally",
                "PROTO10_guarded_common_event_contract_exercised",
                "PROTO10_center_adjacent_binary64_source_floor_localized",
                "PROTO11_well_balanced_reference_map_required",
            ),
            (
                "PROTO10_GR0_case_eligible",
                "PROTO11_frozen",
                "fresh_GR0_dynamic_calibration_completed",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-pro11-frz1.json",
            "FGC-1-PRO11-FRZ1",
            (),
            (
                "PROTO11_successor_runtime_implemented",
                "PROTO11_fresh_GR0_dynamic_calibration_authorized",
                "PROTO11_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "FGCQR_mechanism_rejected",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-hlt9-mon9.json",
            "FGC-1-HLT9-MON9",
            (
                "PROTO11_successor_runtime_implemented",
                "PROTO11_fresh_GR0_dynamic_calibration_authorized",
            ),
            (
                "PROTO11_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "FGCQR_mechanism_rejected",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-cal8-pref10.json",
            "FGC-1-CAL8-PREF10",
            (
                "PROTO11_campaign_terminated_normally",
                "PROTO11_both_amplitudes_reached_evolved_common_event",
                "PROTO11_evolved_affine_source_arithmetic_obstruction_localized",
                "PROTO11_direct_coarse_phi_derivative_veto_localized",
                "SRC2_affine_source_arithmetic_preflight_required",
            ),
            (
                "PROTO11_GR0_case_eligible",
                "PROTO12_frozen",
                "fresh_GR0_dynamic_calibration_completed",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-src2-pref11.json",
            "FGC-1-SRC2-PREF11",
            (
                "PROTO11_source_wall_replayed",
                "captured_GR0_point0_exact_affine_system_nonsingular",
                "captured_GR0_point0_exact_root_exists",
                "captured_GR0_point0_binary64_residual_evaluator_incomplete",
                "captured_GR0_point0_wall_is_a_floating_cancellation_instrument_failure",
                "cancellation_resistant_source_evaluator_required",
            ),
            (
                "CAP1_frozen_arithmetic_routes_cleared_the_wall",
                "captured_GR0_point0_physical_or_structural_breakdown_demonstrated",
                "amplitude_three_spectral_veto_cleared",
                "PROTO12_frozen",
                "GR0_calibration_completed",
                "classical_spherical_diagnostic_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-src3.json",
            "FGC-1-SRC3",
            (
                "SRC3_reference_covariant_binary64_evaluator_derived",
                "SRC3_exact_reference_bitwise_preserved",
                "SRC3_compact_and_independent_exact_oracle_controls_passed",
                "SRC3_captured_full_grid_source_gate_passed",
                "SRC3_unredefined_equations_and_raw_threshold_preserved",
                "SRC3_runtime_source_instrument_authorized_for_a_separately_frozen_GR0_protocol",
            ),
            (
                "amplitude_three_spectral_veto_cleared",
                "PROTO12_frozen",
                "fresh_GR0_calibration_completed",
                "classical_spherical_diagnostic_authorized",
                "SGBL_execution_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-rsp1-frz1.json",
            "FGC-1-RSP1-FRZ1",
            (
                "RSP1_runtime_implemented",
                "RSP1_execution_authorized",
            ),
            (
                "amplitude_three_spectral_veto_cleared",
                "PROTO12_frozen",
                "fresh_GR0_calibration_completed",
                "classical_spherical_diagnostic_authorized",
                "SGBL_execution_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-rsp1-pref12.json",
            "FGC-1-RSP1-PREF12",
            (
                "RSP1_study_terminated_normally",
                "RSP1_all_six_members_reached_endpoint",
                "RSP1_terminal_checkpoint_recomputed",
                "RSP1_source_retry_total_zero",
                "amplitude_three_spectral_veto_cleared_for_successor_design",
                "PROTO12_design_may_begin",
            ),
            (
                "RSP1_generic_all_field_raw_spectral_admission_passed",
                "PROTO12_frozen",
                "fresh_GR0_dynamic_calibration_completed",
                "GR0_case_eligible",
                "classical_spherical_diagnostic_authorized",
                "SGBL_execution_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-src4-vec1.json",
            "FGC-1-SRC4-VEC1",
            (
                "SRC4_tensor_contraction_evaluator_derived",
                "SRC4_exact_reference_and_exact_oracle_controls_passed",
                "SRC4_independent_SRC3_differential_controls_passed",
                "SRC4_CAP1_full_grid_strict_source_gate_passed",
                "SRC4_runtime_source_instrument_authorized_for_a_separately_frozen_GR0_protocol",
            ),
            (
                "wall_clock_speedup_is_a_scientific_result",
                "PROTO12_frozen",
                "fresh_GR0_calibration_completed",
                "classical_spherical_diagnostic_authorized",
                "SGBL_execution_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-pro12-frz1.json",
            "FGC-1-PRO12-FRZ1",
            (
                "PROTO12_frozen",
                "PROTO12_SRC4_source_backend_contract_frozen",
                "PROTO12_pairwise_spectral_classifier_derived",
                "RSP1_pairwise_conditioning_recomputed",
            ),
            (
                "RSP1_generic_all_field_raw_spectral_admission_passed",
                "PROTO12_successor_runtime_implemented",
                "PROTO12_fresh_GR0_dynamic_calibration_authorized",
                "PROTO12_resolved_holdout_manifest_authorized",
                "fresh_GR0_dynamic_calibration_completed",
                "classical_spherical_diagnostic_authorized",
                "SGBL_execution_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-hlt10-mon10.json",
            "FGC-1-HLT10-MON10",
            (
                "PROTO12_successor_runtime_implemented",
                "PROTO12_fresh_GR0_dynamic_calibration_authorized",
            ),
            (
                "PROTO12_resolved_holdout_manifest_authorized",
                "classical_spherical_diagnostic_authorized",
                "SGBL_execution_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
        (
            "results/fgc-1-cal9-pref13.json",
            "FGC-1-CAL9-PREF13",
            (
                "PROTO12_campaign_terminated_normally",
                "PROTO12_both_amplitudes_reached_late_evolved_common_events",
                "SRC4_arithmetic_wall_cleared_for_both_amplitudes",
                "PROTO12_dual_amplitude_SSPRK3_radial_momentum_order_veto_localized",
                "RSP2_high_ladder_constraint_preflight_required",
            ),
            (
                "PROTO12_GR0_case_eligible",
                "PROTO13_frozen",
                "fresh_GR0_dynamic_calibration_completed",
                "classical_spherical_diagnostic_authorized",
                "SGBL_execution_authorized",
                "FGCQR_holdout_execution_authorized",
                "FGCQR_mechanism_rejected",
                "retained_EFT_evolution_authorized",
                "physical_transition_claim_authorized",
                "general_gradient_route_rejected",
                "singularity_resolution_derived",
                "child_domain_or_topology_derived",
                "dark_sector_mechanism_derived",
                "varying_locally_measured_c_derived",
            ),
        ),
    )
    for relative, artifact_id, required_true, required_false in post_comp1_records:
        failures.extend(
            _check_hash_bound_post_comp1_result(
                relative,
                artifact_id,
                required_true_gates=required_true,
                required_false_gates=required_false,
            )
        )
    failures.extend(_check_cal7_pref9_result())
    failures.extend(_check_pro11_frz1_result())
    failures.extend(_check_hlt9_mon9_result())
    failures.extend(_check_cal8_pref10_result())
    failures.extend(_check_src2_cap1_contract())
    failures.extend(_check_src2_pref11_result())
    failures.extend(_check_src3_result())
    failures.extend(_check_rsp1_frz1_result())
    failures.extend(_check_rsp1_pref12_result())
    failures.extend(_check_src4_vec1_result())
    failures.extend(_check_pro12_frz1_result())
    failures.extend(_check_hlt10_mon10_result())
    failures.extend(_check_cal9_pref13_result())
    failures.extend(_check_rsp2_frz1_result())
    failures.extend(_check_rsp2_pref14_result())
    failures.extend(_check_pro13_frz1_result())
    failures.extend(_check_hlt11_mon11_result())
    failures.extend(_check_cal10_pref15_result())
    failures.extend(_check_tdg1_frz1_result())
    failures.extend(_check_tdg1_pref16_result())
    failures.extend(_check_tdg2_frz1_result())
    failures.extend(_check_tdg2_pref17_result())
    failures.extend(_check_tdg3_frz1_result())
    failures.extend(_check_tdg3_pref18_result())
    failures.extend(_check_tdg4_frz1_result())
    failures.extend(_check_tdg4_pref19_result())
    failures.extend(_check_tdg5_frz1_result())
    failures.extend(_check_tdg5_pref20_result())
    failures.extend(_check_tdg5_imp1_result())
    failures.extend(_check_tdg6_frz1_result())
    failures.extend(_check_tdg6_pref21_result())
    failures.extend(_check_tdg6_imp2_result())
    failures.extend(_check_pro14_frz1_result())
    failures.extend(_check_hlt12_mon12_result())
    failures.extend(_check_cal11_pref22_result())
    failures.extend(_check_tdg7_frz1_result())
    failures.extend(_check_tdg7_pref23_result())
    failures.extend(_check_tdg7_imp3_result())
    failures.extend(_check_pro15_frz1_result())
    failures.extend(_check_hlt13_mon13_result())
    failures.extend(_check_pro16_frz1_result())
    failures.extend(_check_pro16_pref24_result())
    failures.extend(_check_pro17_frz1_result())
    failures.extend(_check_pro17_pref25_result())
    failures.extend(_check_hlt14_mon14_result())
    failures.extend(_check_pro18_frz1_result())
    failures.extend(_check_pro18_pref26_result(run_raw_reproduction=False))
    failures.extend(_check_pro18_auth1_result())
    failures.extend(_check_pro18_pref27_result())
    failures.extend(_check_pro19_frz1_result())
    failures.extend(_check_hlt16_mon16_result())
    failures.extend(_check_pro19_prelaunch_result())
    failures.extend(_check_pro19_cfl1_compact_result())
    failures.extend(_check_pro19_sid1_compact_result())
    failures.extend(_check_pro19_sid2_compact_result())
    failures.extend(_check_pro19_sid3_compact_result())
    failures.extend(_check_pro19_pref28_compact_result())
    failures.extend(_check_tdg8_successor_prelaunch_result())
    failures.extend(_check_tdg8_rcv1_prelaunch_result())
    failures.extend(_check_tdg8_rcv2_prelaunch_result())
    failures.extend(_check_tdg8_rcv3_prelaunch_result())
    failures.extend(_check_tdg8_rcv3_pref1_result())
    failures.extend(_check_tdg8_rcv3_auth1_result())
    failures.extend(_check_tdg8_rcv3_rec1_auth1_result())
    failures.extend(_check_tdg8_rcv3_pref2_result())
    failures.extend(_check_tdg9_ar1_auth1_result())
    failures.extend(_check_tdg9_ar1_pref1_result())
    failures.extend(_check_tdg9_loc1_frz1_result())
    failures.extend(_check_tdg9_loc2_frz1_result())
    failures.extend(_check_tdg9_loc2_pref2_result())
    failures.extend(_check_tdg9_ti1_frz1_result())
    failures.extend(_check_tdg9_ti2_frz1_result())
    failures.extend(_check_tdg9_ti2_pref1_result())
    failures.extend(_check_tdg9_ac1_frz1_result())
    failures.extend(_check_tdg9_ac1_pref1_result())
    failures.extend(_check_tdg9_ur1_frz1_result())
    failures.extend(_check_tdg9_ur1_pref1_result())
    failures.extend(_check_tdg10_qa1_frz1_result())
    failures.extend(_check_tdg10_qa1_pref1_result())
    failures.extend(_check_tdg10_qa2_frz1_result())
    failures.extend(_check_tdg10_qa2_pref1_result())
    failures.extend(_check_tdg11_msel1_pref1_result())
    failures.extend(_check_tdg11_imp1_result())
    failures.extend(_check_hlt17_srcq1_rec1_pref1_result())
    failures.extend(_check_pro20_ev1_pref1_result())
    failures.extend(_check_def1_stab1_frz1_result())
    failures.extend(_check_def1_stab1_pref1_result())
    failures.extend(_check_sgb1_ctl1_sol1_frz1_result())
    failures.extend(_check_sgb1_ctl1_sol1_pref1_result())

    dom4_path = REPOSITORY / "results/fgc-1-dom4-run1.json"
    hyp2_path = REPOSITORY / "results/fgc-1-hyp2-md1.json"
    try:
        dom4 = json.loads(
            dom4_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
        hyp2 = json.loads(
            hyp2_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect DOM4/HYP2 certificates: {exc}")
        return failures
    dom4_payload = dom4.get("artifact_payload", {})
    dom4_quantitative = dom4_payload.get("quantitative_evidence", {})
    dom4_continuation = dom4_quantitative.get("seed_branch_continuation", {})
    dom4_stress = dom4_quantitative.get("stress_initial_slice", {})
    dom4_boundary = dom4_quantitative.get("epistemic_boundary", {})
    domain_hash = dom4_quantitative.get("domain_definition_sha256")
    if (
        dom4.get("classification")
        != "classical_spherical_run_domain_certificate"
        or dom4.get("generated_by") != "scripts/reproduce_fgc_dom4_run1.py"
        or domain_hash
        != "fc27f0f6617e6739807af17fca6219da5d9064ea1a3a5a782778532cb6f5164e"
        or dom4_quantitative.get("parameter_container", {}).get(
            "dimensionless_coordinate_volume"
        )
        != "9/4"
        or dom4_quantitative.get("parameter_container", {}).get(
            "all_PROTO3_candidate_and_holdout_coordinates_contained"
        )
        is not True
        or [
            item.get("phi_seed_factor")
            for item in dom4_continuation.get("records", [])
            if isinstance(item, dict)
        ]
        != ["0", "1", "2"]
        or dom4_continuation.get("strict_nonsingular_branch_margins") is not True
        or dom4_stress.get("regular_center") is not True
        or dom4_stress.get("finite_mass") is not True
        or dom4_stress.get("no_initial_trapped_sphere") is not True
        or dom4_stress.get("inside_protocol_compactness_window") is not True
        or dom4_boundary.get("local_open_state_neighborhoods_are_not_a_global_run_tube")
        is not True
        or dom4_boundary.get("time_evolution_performed") is not False
        or dom4_boundary.get("FGCQR_holdout_outcome_inspected") is not False
    ):
        failures.append("FGC-1-DOM4-RUN1 domain, witness, or nonclaim contract has drifted")

    hyp2_payload = hyp2.get("artifact_payload", {})
    hyp2_quantitative = hyp2_payload.get("quantitative_evidence", {})
    hyp2_crosscheck = hyp2_quantitative.get("independent_exact_radial_crosscheck", {})
    hyp2_extrema = hyp2_quantitative.get("witness_extrema", {})
    hyp2_boundary = hyp2_quantitative.get("epistemic_boundary", {})
    hyp2_records = hyp2_quantitative.get("witness_records", [])
    if (
        hyp2.get("classification")
        != "classical_all_covector_early_kill_certificate"
        or hyp2.get("generated_by") != "scripts/reproduce_fgc_hyp2_md1.py"
        or hyp2_quantitative.get("domain_definition_sha256") != domain_hash
        or hyp2_crosscheck.get("independent_exact_spherical_symbol_agreement")
        is not True
        or not isinstance(hyp2_crosscheck.get("maximum_absolute_root_difference"), (int, float))
        or hyp2_crosscheck.get("maximum_absolute_root_difference") > 2e-13
        or not isinstance(hyp2_extrema.get("maximum_companion_deformation_operator_2"), (int, float))
        or hyp2_extrema.get("maximum_companion_deformation_operator_2") >= 1 / 4096
        or not isinstance(hyp2_extrema.get("maximum_action_energy_deformation_operator_2"), (int, float))
        or hyp2_extrema.get("maximum_action_energy_deformation_operator_2") >= 1 / 8192
        or not isinstance(hyp2_extrema.get("minimum_cluster_resolvent_singular_margin"), (int, float))
        or hyp2_extrema.get("minimum_cluster_resolvent_singular_margin") <= 0
        or not isinstance(hyp2_extrema.get("minimum_physical_energy_coercivity_lower"), (int, float))
        or hyp2_extrema.get("minimum_physical_energy_coercivity_lower") <= 0
        or not isinstance(hyp2_records, list)
        or len(hyp2_records) != 4
        or any(
            not isinstance(record, dict)
            or record.get("passed") is not True
            or not all(record.get("predicates", {}).values())
            for record in hyp2_records
        )
        or hyp2_boundary.get("trajectory_containment_is_deferred_to_HLT1")
        is not True
        or hyp2_boundary.get("retained_EFT_remainder_or_frequency_control_supplied")
        is not False
        or hyp2_boundary.get("time_evolution_performed") is not False
        or hyp2_boundary.get("FGCQR_holdout_outcome_inspected") is not False
    ):
        failures.append("FGC-1-HYP2-MD1 all-covector margin or nonclaim contract has drifted")

    bnd2_path = REPOSITORY / "results/fgc-1-bnd2-cp1.json"
    hlt1_path = REPOSITORY / "results/fgc-1-hlt1-mon1.json"
    try:
        bnd2 = json.loads(
            bnd2_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
        hlt1 = json.loads(
            hlt1_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect BND2/HLT1 certificates: {exc}")
        return failures
    bnd2_payload = bnd2.get("artifact_payload", {})
    bnd2_quantitative = bnd2_payload.get("quantitative_evidence", {})
    bnd2_move = bnd2_quantitative.get("outer_boundary_move_control", {})
    bnd2_budget = bnd2_quantitative.get("causal_budget_contract", {})
    bnd2_runtime = bnd2_budget.get("runtime_controls", {})
    bnd2_boundary = bnd2_quantitative.get("epistemic_boundary", {})
    if (
        bnd2.get("classification")
        != "spherical_domain_of_dependence_and_boundary_certificate"
        or bnd2.get("generated_by") != "scripts/reproduce_fgc_bnd2_cp1.py"
        or any(
            bnd2_payload.get(key) is not True
            for key in (
                "evolving_physical_and_auxiliary_cones_included",
                "SBP_stencil_reach_included",
                "measured_affine_interval_outside_boundary_domain_of_dependence",
                "outer_boundary_move_check_passed",
                "incoming_constraint_characteristics_controlled",
            )
        )
        or bnd2_move.get("smallest_reference_margin_over_required_buffer")
        != "689/40"
        or len(bnd2_move.get("grid_records", [])) != 9
        or bnd2_budget.get("candidate_endpoint_speed_is_evaluated_explicitly")
        is not True
        or bnd2_runtime.get("accepted_control", {}).get(
            "accepted_endpoint_speed_upper"
        )
        != 1.3
        or bnd2_runtime.get("injected_failure", {}).get("typed_reason")
        != "boundary_causal_buffer"
        or bnd2_runtime.get("injected_failure", {}).get("trial_rejected")
        is not True
        or bnd2_runtime.get(
            "last_accepted_state_unchanged_by_failed_pure_transaction"
        )
        is not True
        or bnd2_boundary.get(
            "complete_nonlinear_constraint_preserving_IBVP_proved"
        )
        is not False
        or bnd2_boundary.get("time_evolution_performed") is not False
        or bnd2_boundary.get("FGCQR_holdout_outcome_inspected") is not False
    ):
        failures.append("FGC-1-BND2-CP1 causal-budget or nonclaim contract has drifted")

    hlt1_payload = hlt1.get("artifact_payload", {})
    hlt1_quantitative = hlt1_payload.get("quantitative_evidence", {})
    hlt1_injections = hlt1_quantitative.get("typed_stop_injections", {})
    hlt1_transaction = hlt1_injections.get("transaction_control", {})
    hlt1_boundary = hlt1_quantitative.get("epistemic_boundary", {})
    if (
        hlt1.get("classification")
        != "classical_health_and_typed_stop_monitor_certificate"
        or hlt1.get("generated_by") != "scripts/reproduce_fgc_hlt1_mon1.py"
        or any(
            hlt1_payload.get(key) is not True
            for key in (
                "canonical_dimensionless_norms_executable",
                "all_protocol_thresholds_implemented_strictly",
                "every_accepted_stage_monitored",
                "nonmonotone_internal_stage_abscissae_supported",
                "strict_global_stage_transaction_order_enforced",
                "stop_occurs_before_threshold_crossing",
                "first_failed_premise_is_preserved",
                "scale_stops_not_promoted_to_EFT_validity",
            )
        )
        or len(hlt1_injections.get("records", {})) != 26
        or hlt1_injections.get("every_reason_independently_stopped") is not True
        or hlt1_transaction.get("simultaneous_first_reason")
        != "nonpositive_lapse"
        or hlt1_transaction.get("first_failure_immutable_on_later_attempt")
        is not True
        or hlt1_transaction.get("SSPRK3_nonmonotone_stage_abscissae_supported")
        is not True
        or hlt1_transaction.get("strict_global_stage_transaction_order_enforced")
        is not True
        or hlt1_transaction.get(
            "ordering_error_does_not_mutate_last_accepted_state"
        )
        is not True
        or hlt1_boundary.get(
            "scale_thresholds_are_classical_stop_proxies_not_EFT_remainder_control"
        )
        is not True
        or hlt1_boundary.get("time_evolution_performed") is not False
        or hlt1_boundary.get("FGCQR_holdout_outcome_inspected") is not False
    ):
        failures.append("FGC-1-HLT1-MON1 typed-stop or nonclaim contract has drifted")

    num1_path = REPOSITORY / "results/fgc-1-num1-val1.json"
    try:
        num1 = json.loads(
            num1_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-NUM1-VAL1 certificate: {exc}")
        return failures
    num1_payload = num1.get("artifact_payload", {})
    num1_quantitative = num1_payload.get("quantitative_evidence", {})
    num1_methods = num1_quantitative.get("method_controls", [])
    num1_source = num1_quantitative.get("exact_scalar_and_batched_source", {})
    num1_center = num1_quantitative.get("regular_center_acceleration_limit", {})
    num1_runtime = num1_quantitative.get("runtime_transaction", {})
    num1_affine = num1_quantitative.get("affine_measurement", {})
    num1_affine_records = num1_affine.get("records", [])
    num1_boundary = num1_quantitative.get("epistemic_boundary", {})
    required_num1_payload = {
        "Minkowski_preservation_passed",
        "regular_center_manufactured_solution_passed",
        "GR0_scalar_control_passed",
        "exact_vs_floating_backend_passed",
        "constraint_convergence_passed",
        "boundary_isolation_passed",
        "all_typed_stops_injected",
        "checkpoint_restart_equivalence_passed",
        "direct_and_Raychaudhuri_routes_agree",
        "primary_and_comparator_methods_passed",
    }
    expected_num1_methods = {
        "fourth_order_diagonal_norm_SBP_with_second_order_boundary_closure_plus_RK4",
        "second_order_diagonal_norm_SBP_plus_SSPRK3",
    }
    method_contract_passed = (
        isinstance(num1_methods, list)
        and len(num1_methods) == 2
        and {item.get("method") for item in num1_methods if isinstance(item, dict)}
        == expected_num1_methods
        and all(
            isinstance(item, dict)
            and item.get("all_method_controls_passed") is True
            and item.get("Minkowski", {}).get("bitwise_preserved") is True
            and item.get("Minkowski", {}).get("maximum_binary64_drift") == 0.0
            and item.get("regular_center_manufactured_solution", {}).get(
                "minimum_constraint_order", 0.0
            )
            > 1.5
            and item.get("GR0_scalar", {}).get("minimum_constraint_order", 0.0)
            > 1.5
            and item.get("checkpoint_restart", {}).get(
                "canonical_roundtrip_bitwise"
            )
            is True
            and item.get("checkpoint_restart", {}).get(
                "continuous_and_restart_endpoint_bitwise_equal"
            )
            is True
            and item.get("boundary_isolation", {}).get(
                "maximum_measurement_difference"
            )
            == 0.0
            for item in num1_methods
        )
    )
    affine_contract_passed = (
        isinstance(num1_affine_records, list)
        and len(num1_affine_records) == 2
        and {item.get("method") for item in num1_affine_records if isinstance(item, dict)}
        == expected_num1_methods
        and all(
            isinstance(item, dict)
            and abs(item.get("initial_normalization_residual", 1.0)) <= 1.0e-10
            and item.get("maximum_null_residual", 1.0) <= 1.0e-12
            and item.get("maximum_affine_geodesic_residual", 1.0) <= 1.0e-10
            and item.get("maximum_direct_Raychaudhuri_disagreement", 1.0)
            <= 1.0e-8
            and item.get("radial_screen_shear_squared") == 0.0
            and item.get("hypersurface_twist_squared") == 0.0
            for item in num1_affine_records
        )
        and num1_affine.get("cross_method_theta_difference_infinity", 1.0)
        <= 1.0e-8
    )
    if (
        num1.get("classification")
        != "independent_solver_and_affine_measurement_validation_certificate"
        or num1.get("generated_by") != "scripts/reproduce_fgc_num1_val1.py"
        or any(num1_payload.get(key) is not True for key in required_num1_payload)
        or not method_contract_passed
        or num1_quantitative.get("typed_stop_count") != 26
        or num1_runtime.get("all_transactions_passed") is not True
        or num1_center.get("minimum_observed_order", 0.0) < 5.0
        or num1_source.get("scalar_vs_batch_ADM_root_absolute_error", 1.0)
        > 1.0e-12
        or num1_source.get("exact_vs_batch_ADM_root_absolute_error", 1.0)
        > 1.0e-12
        or num1_source.get("batch_REF1_residual_infinity", 1.0) > 1.0e-12
        or num1_source.get("batch_affine_in_accelerations_verified") is not True
        or not affine_contract_passed
        or num1_boundary.get("controls_are_independent_of_FGCQR_holdout_outcomes")
        is not True
        or num1_boundary.get("FGCQR_trajectory_evaluated") is not False
        or num1_boundary.get("collapse_evolution_performed") is not False
        or num1_boundary.get("trapped_interval_measured") is not False
        or num1_boundary.get("regulator_activation_measured") is not False
        or num1_boundary.get("defocusing_margin_measured") is not False
        or num1_boundary.get("retained_EFT_validity_inferred") is not False
    ):
        failures.append(
            "FGC-1-NUM1-VAL1 numerical-control or scientific-boundary contract has drifted"
        )

    cal0_path = REPOSITORY / "results/fgc-1-cal0-pref2.json"
    try:
        cal0 = json.loads(
            cal0_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-CAL0-PREF2 certificate: {exc}")
        return failures
    cal0_payload = cal0.get("artifact_payload", {})
    cal0_spectral = cal0_payload.get("spectral_input_obstruction", {})
    cal0_records = cal0_spectral.get("records", [])
    cal0_applicability = cal0_payload.get("branch_applicability_obstruction", {})
    cal0_boundary = cal0_payload.get("epistemic_boundary", {})
    amplitude_groups: dict[str, list[dict[str, object]]] = {}
    if isinstance(cal0_records, list):
        for item in cal0_records:
            if isinstance(item, dict) and isinstance(item.get("amplitude"), str):
                amplitude_groups.setdefault(item["amplitude"], []).append(item)
    tails_converge = all(
        len(items) == 3
        and [item.get("point_count") for item in items] == [1025, 2049, 4097]
        and all(item.get("alias_band_occupied") is True for item in items)
        and all(item.get("support_bin") == item.get("nyquist_bin") for item in items)
        and all(
            items[index + 1].get("top_eighth_field_energy_fraction", 1.0)
            < items[index].get("top_eighth_field_energy_fraction", 0.0)
            and items[index + 1].get(
                "top_eighth_gradient_energy_fraction", 1.0
            )
            < items[index].get("top_eighth_gradient_energy_fraction", 0.0)
            for index in range(2)
        )
        for items in amplitude_groups.values()
    )
    if (
        cal0.get("classification")
        != "pre_holdout_protocol_numerical_contract_obstruction"
        or cal0.get("generated_by") != "scripts/reproduce_fgc_cal0_pref2.py"
        or cal0_spectral.get("estimator") != "HLT1_windowed_spectral_support"
        or cal0_spectral.get("all_frozen_resolutions_fail_at_initial_time")
        is not True
        or cal0_spectral.get("nonzero_amplitude_scaling_invariance_verified")
        is not True
        or len(cal0_records) != 6
        or set(amplitude_groups) != {"2", "5"}
        or not tails_converge
        or cal0_applicability.get("calibration_branch") != "GR-0"
        or cal0_applicability.get("PROTO3_all_RUN1_stops_apply") is not True
        or cal0_applicability.get("target_branch_owners")
        != {"HLT1": "FGC-QR", "DOM4": "FGC-QR", "HYP2": "FGC-QR"}
        or cal0_applicability.get("branch_applicability_map_present") is not False
        or cal0_applicability.get("scope_conflation_verified") is not True
        or cal0_boundary.get("protocol_design_route_closed")
        != "PROTO3_as_executable_resolved_holdout_contract"
        or cal0_boundary.get("FGCQR_evolution_outcome_inspected") is not False
        or cal0_boundary.get("SGBL_evolution_outcome_inspected") is not False
        or cal0_boundary.get("FGCQR_action_or_mechanism_tested") is not False
        or cal0_boundary.get("FGCQR_action_or_mechanism_rejected") is not False
        or cal0_boundary.get("general_gradient_mechanism_rejected") is not False
    ):
        failures.append(
            "FGC-1-CAL0-PREF2 protocol-obstruction or scientific-boundary contract has drifted"
        )

    pro4_path = REPOSITORY / "results/fgc-1-pro4-frz1.json"
    try:
        pro4 = json.loads(
            pro4_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-PRO4-FRZ1 certificate: {exc}")
        return failures
    pro4_payload = pro4.get("artifact_payload", {})
    pro4_protocol = pro4_payload.get("protocol_certificate", {})
    pro4_lineage = pro4_payload.get("immutable_lineage", {})
    pro4_predecessor = pro4_lineage.get("predecessor_protocol", {})
    pro4_diagnosis = pro4_lineage.get("diagnosis", {})
    pro4_inherited = pro4_payload.get("inherited_physical_contract", {})
    pro4_cases = pro4_inherited.get("immutable_case_coordinates", {})
    pro4_revision = pro4_payload.get("numerical_contract_revision", {})
    pro4_namespace = pro4_payload.get("namespace_contract", {})
    pro4_boundary = pro4_payload.get("epistemic_boundary", {})
    if (
        pro4.get("classification")
        != "outcome_neutral_successor_protocol_freeze_certificate"
        or pro4.get("generated_by") != "scripts/reproduce_fgc_pro4_frz1.py"
        or pro4_protocol.get("artifact_id") != "FGC-2-SF1-PROTO4"
        or pro4_protocol.get("protocol_version") != 4
        or pro4_protocol.get("outcome_neutral_contract_validated") is not True
        or pro4_protocol.get("premise_revision_only") is not True
        or pro4_protocol.get("resolved_holdout_manifest_present") is not False
        or pro4_lineage.get("checkpoint_commit")
        != "b5aa4934b58afd505b22e008202285c475942892"
        or pro4_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or pro4_predecessor.get("artifact_id") != "FGC-2-SF1-PROTO3"
        or pro4_predecessor.get("blob_sha256")
        != "57b28eb9abe80c81f9e47ae6bea6ce2d69e103542df500e97c6c93805a8ac13b"
        or pro4_diagnosis.get("artifact_id") != "FGC-1-CAL0-PREF2"
        or pro4_diagnosis.get("blob_sha256")
        != "d952e737d22d63b33e8d87e366538925bd28acbf82a3062563de49867a229138"
        or pro4_diagnosis.get("PROTO3_obstruction_verified") is not True
        or pro4_diagnosis.get("FGCQR_outcome_inspected") is not False
        or pro4_diagnosis.get("SGBL_outcome_inspected") is not False
        or pro4_inherited.get("inheritance_is_by_exact_PROTO3_blob_hash") is not True
        or pro4_inherited.get("action_couplings_profiles_and_equations_changed")
        is not False
        or pro4_inherited.get(
            "null_observable_outcome_rules_and_robustness_changed"
        )
        is not False
        or pro4_cases.get("amplitude_candidates")
        != ["2", "5/2", "3", "7/2", "4", "9/2", "5"]
        or pro4_cases.get("held_out_case_ids")
        != [
            "FGCQR-CENTRAL",
            "FGCQR-SEED-HALF",
            "FGCQR-SEED-DOUBLE",
            "FGCQR-WIDTH-SEVEN-EIGHTHS",
            "FGCQR-WIDTH-NINE-EIGHTHS",
        ]
        or pro4_revision.get("legacy_stop_count") != 26
        or pro4_revision.get("universal_runtime_stop_count") != 9
        or pro4_revision.get("candidate_branch_only_stop_count") != 11
        or pro4_revision.get("replaced_stop_count") != 6
        or pro4_revision.get("partition_is_exact_and_disjoint") is not True
        or pro4_revision.get("last_bin_boolean_is_diagnostic_only") is not True
        or pro4_revision.get(
            "three_nested_common_event_admission_required"
        )
        is not True
        or pro4_revision.get(
            "Richardson_error_enters_conservative_observable_error_sum"
        )
        is not True
        or pro4_namespace.get("resolved_holdout_manifest")
        != "configs/fgc/fgc-1-pro4-hld1.toml"
        or pro4_namespace.get("output_namespace_contents_accessed_by_this_certificate")
        is not False
        or pro4_boundary.get("GR0_exploratory_outcomes_were_seen_before_freeze")
        is not True
        or pro4_boundary.get(
            "GR0_exploratory_runs_admissible_as_calibration_evidence"
        )
        is not False
        or pro4_boundary.get("FGCQR_evolution_outcome_inspected") is not False
        or pro4_boundary.get("SGBL_evolution_outcome_inspected") is not False
        or pro4_boundary.get("HLT2_monitor_implemented") is not False
        or pro4_boundary.get("FGCQR_action_or_mechanism_tested") is not False
        or pro4_boundary.get("general_gradient_route_rejected") is not False
        or any(value is not False for value in pro4.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-PRO4-FRZ1 lineage, numerical-contract, or nonclaim boundary has drifted"
        )

    hlt2_path = REPOSITORY / "results/fgc-1-hlt2-mon2.json"
    try:
        hlt2 = json.loads(
            hlt2_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-HLT2-MON2 certificate: {exc}")
        return failures
    hlt2_payload = hlt2.get("artifact_payload", {})
    hlt2_branches = hlt2_payload.get("branch_applicability", {})
    hlt2_GR0 = hlt2_branches.get("GR0", {})
    hlt2_FGCQR = hlt2_branches.get("FGCQR_with_exact_owner", {})
    hlt2_SGBL = hlt2_branches.get("SGBL_without_owner", {})
    hlt2_evidence = hlt2_payload.get("quantitative_evidence", {})
    hlt2_compact = hlt2_evidence.get(
        "declared_compact_initial_profile_control", {}
    )
    hlt2_compact_records = hlt2_compact.get("records", [])
    hlt2_nested = hlt2_compact.get("nested_admission", {})
    hlt2_constraints = hlt2_evidence.get("constraint_and_error_controls", {})
    hlt2_boundary = hlt2_payload.get("epistemic_boundary", {})
    hlt2_tail_ratios = []
    for key in (
        "field_power_tail_ratios",
        "derivative_power_tail_ratios",
    ):
        table = hlt2_nested.get(key, {})
        if isinstance(table, dict):
            for values in table.values():
                if isinstance(values, list):
                    hlt2_tail_ratios.extend(values)
    hlt2_field_budgets = [
        budget
        for item in hlt2_compact_records
        if isinstance(item, dict)
        for budget in (
            item.get("field_budgets", {}).values()
            if isinstance(item.get("field_budgets"), dict)
            else ()
        )
        if isinstance(budget, dict)
    ]
    if (
        hlt2.get("classification")
        != "outcome_neutral_PROTO4_successor_admission_monitor_certificate"
        or hlt2.get("generated_by") != "scripts/reproduce_fgc_hlt2_mon2.py"
        or hlt2.get("scope_bindings", {}).get("protocol_artifact_id")
        != "FGC-2-SF1-PROTO4"
        or hlt2.get("scope_bindings", {}).get("calibration_branch") != "GR-0"
        or hlt2.get("scope_bindings", {}).get("comparison_branch") != "SGB-L"
        or hlt2.get("scope_bindings", {}).get("holdout_branch") != "FGC-QR"
        or hlt2_GR0.get("runtime_monitor_complete") is not True
        or len(hlt2_GR0.get("universal_stop_ids", [])) != 9
        or hlt2_GR0.get("candidate_branch_stop_ids") != []
        or len(hlt2_GR0.get("replaced_stop_ids", [])) != 6
        or hlt2_FGCQR.get("runtime_monitor_complete") is not True
        or len(hlt2_FGCQR.get("candidate_branch_stop_ids", [])) != 11
        or hlt2_FGCQR.get("candidate_definition_owner")
        != "FGC-1-SRC1-NL1+FGC-1-DOM4-RUN1+FGC-1-HYP2-MD1"
        or hlt2_SGBL.get("runtime_monitor_complete") is not False
        or hlt2_SGBL.get("candidate_branch_stop_ids") != []
        or hlt2_branches.get("SGBL_rejected_FGCQR_cross_owner") is not True
        or [item.get("point_count") for item in hlt2_compact_records]
        != [1025, 2049, 4097]
        or any(
            item.get("every_individual_budget_passed") is not True
            for item in hlt2_compact_records
            if isinstance(item, dict)
        )
        or len(hlt2_field_budgets) != 18
        or any(
            budget.get("individual_admission_passed") is not True
            for budget in hlt2_field_budgets
        )
        or not any(
            budget.get("last_bin_alias_diagnostic") is True
            for budget in hlt2_field_budgets
        )
        or hlt2_nested.get("admission_passed") is not True
        or len(hlt2_tail_ratios) != 24
        or any(
            not isinstance(value, (int, float))
            or isinstance(value, bool)
            or value < 0
            or value >= 0.25
            for value in hlt2_tail_ratios
        )
        or hlt2_constraints.get(
            "manufactured_second_order_common_event", {}
        ).get("admission_passed")
        is not True
        or hlt2_constraints.get("nonconvergent_injection", {}).get(
            "admission_passed"
        )
        is not False
        or hlt2_constraints.get(
            "common_event_misalignment_injection", {}
        ).get("admission_passed")
        is not False
        or hlt2_constraints.get(
            "DEF1_missing_affine_alignment_injection", {}
        ).get("common_event_alignment_passed")
        is not False
        or hlt2_constraints.get("DEF1_affine_alignment_control", {}).get(
            "common_event_alignment_passed"
        )
        is not True
        or hlt2_constraints.get("Richardson_additive_control", {}).get(
            "subtraction_available"
        )
        is not False
        or hlt2_constraints.get("Richardson_additive_control", {}).get(
            "passed_common_event_required"
        )
        is not True
        or hlt2_constraints.get("Richardson_additive_control", {}).get(
            "misaligned_common_event_rejected"
        )
        is not True
        or hlt2_boundary.get("dynamic_trajectory_used") is not False
        or hlt2_boundary.get("SGBL_health_definition_missing") is not True
        or hlt2_boundary.get("constraint_to_Raychaudhuri_stability_map_missing")
        is not True
        or hlt2_boundary.get("mechanism_question_answered") is not False
        or any(value is not False for value in hlt2.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-HLT2-MON2 branch, spectral, constraint, or nonclaim contract has drifted"
        )

    id2_path = REPOSITORY / "results/fgc-1-id2-all1.json"
    try:
        id2 = json.loads(
            id2_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-ID2-ALL1 certificate: {exc}")
        return failures
    id2_payload = id2.get("artifact_payload", {})
    id2_partition = id2_payload.get("protocol_partition", {})
    id2_summary = id2_payload.get("static_decision_summary", {})
    id2_eligibility = id2_payload.get("dynamic_calibration_eligibility", [])
    id2_GR0 = id2_payload.get("GR0_candidate_ledger", [])
    id2_FGCQR = id2_payload.get("FGCQR_all_amplitude_all_case_ledger", [])
    id2_boundary = id2_payload.get("epistemic_boundary", {})
    id2_next = id2_payload.get("next_stage_contract", {})
    expected_id2_eligibility = [
        (
            "2",
            False,
            ["FGCQR-WIDTH-NINE-EIGHTHS:peak_compactness_below_protocol_minimum"],
        ),
        ("5/2", True, []),
        ("3", True, []),
        ("7/2", False, ["PROTO4_nested_weighted_spectral_admission"]),
        (
            "4",
            False,
            ["FGCQR-WIDTH-SEVEN-EIGHTHS:individual_weighted_spectral_budget"],
        ),
        (
            "9/2",
            False,
            [
                "PROTO4_nested_weighted_spectral_admission",
                "FGCQR-WIDTH-SEVEN-EIGHTHS:individual_weighted_spectral_budget",
            ],
        ),
        (
            "5",
            False,
            [
                "PROTO4_nested_weighted_spectral_admission",
                "FGCQR-WIDTH-SEVEN-EIGHTHS:peak_compactness_above_protocol_maximum",
            ],
        ),
    ]
    actual_id2_eligibility = [
        (
            item.get("amplitude"),
            item.get("eligible_for_fresh_dynamic_GR0_calibration"),
            item.get("ordered_ineligibility_reasons"),
        )
        for item in id2_eligibility
        if isinstance(item, dict)
    ]
    id2_state_hashes = [
        resolution.get("state_sha256")
        for branch_ledger in (id2_GR0, id2_FGCQR)
        for candidate in branch_ledger
        if isinstance(candidate, dict)
        for method in candidate.get("methods", [])
        if isinstance(method, dict)
        for resolution in method.get("resolutions", [])
        if isinstance(resolution, dict)
    ]
    id2_min_compactness = id2_summary.get(
        "minimum_observed_FGCQR_peak_compactness"
    )
    id2_max_compactness = id2_summary.get(
        "maximum_observed_FGCQR_peak_compactness"
    )
    if (
        id2.get("classification")
        != "pre_calibration_all_amplitude_all_case_static_initial_admission_certificate"
        or id2.get("generated_by") != "scripts/reproduce_fgc_id2_all1.py"
        or id2.get("scope_bindings", {}).get("protocol_artifact_id")
        != "FGC-2-SF1-PROTO4"
        or id2.get("scope_bindings", {}).get("physical_predecessor_protocol_artifact_id")
        != "FGC-2-SF1-PROTO3"
        or id2.get("scope_bindings", {}).get("calibration_branch") != "GR-0"
        or id2.get("scope_bindings", {}).get("holdout_branch") != "FGC-QR"
        or id2_partition.get("GR0_candidate_count") != 7
        or id2_partition.get("FGCQR_static_slice_count") != 35
        or id2_partition.get("method_count") != 2
        or id2_partition.get("resolution_count") != 3
        or id2_partition.get("complete_grid_state_hash_count") != 252
        or id2_partition.get("cartesian_product_completed_before_dynamic_calibration")
        is not True
        or id2_partition.get("every_complete_grid_state_hash_unique") is not True
        or len(id2_GR0) != 7
        or len(id2_FGCQR) != 35
        or len(id2_state_hashes) != 252
        or len(set(id2_state_hashes)) != 252
        or any(
            not isinstance(value, str)
            or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)
            for value in id2_state_hashes
        )
        or actual_id2_eligibility != expected_id2_eligibility
        or any(
            item.get("dynamic_outcome_used") is not False
            for item in id2_eligibility
            if isinstance(item, dict)
        )
        or any(
            item.get("dynamic_trajectory_read") is not False
            for item in id2_GR0
            if isinstance(item, dict)
        )
        or any(
            item.get("FGCQR_evolution_outcome_inspected") is not False
            for item in id2_FGCQR
            if isinstance(item, dict)
        )
        or id2_summary.get("FGCQR_static_slice_pass_count") != 31
        or id2_summary.get("FGCQR_static_slice_fail_count") != 4
        or id2_summary.get("FGCQR_nested_spectral_diagnostic_pass_count") != 17
        or id2_summary.get("FGCQR_nested_spectral_diagnostic_fail_count") != 18
        or id2_summary.get("eligible_amplitudes_in_frozen_order") != ["5/2", "3"]
        or id2_summary.get("ineligible_amplitudes_in_frozen_order")
        != ["2", "7/2", "4", "9/2", "5"]
        or not isinstance(id2_min_compactness, (int, float))
        or isinstance(id2_min_compactness, bool)
        or not 0.0953 < id2_min_compactness < 0.0955
        or not isinstance(id2_max_compactness, (int, float))
        or isinstance(id2_max_compactness, bool)
        or not 0.7938 < id2_max_compactness < 0.7941
        or id2_summary.get("amplitude_2_widest_case_below_compactness_floor")
        is not True
        or id2_summary.get("amplitude_5_narrowest_case_above_compactness_ceiling")
        is not True
        or id2_summary.get(
            "amplitudes_4_and_9_over_2_narrow_case_miss_an_individual_spectral_budget"
        )
        is not True
        or id2_summary.get(
            "nested_FGCQR_spectral_result_is_disclosed_but_not_a_PROTO3_static_five_boolean_decision"
        )
        is not True
        or id2_next.get("fresh_dynamic_GR0_candidates_must_run_in_this_order")
        != ["5/2", "3"]
        or id2_next.get("dynamic_calibration_not_performed_here") is not True
        or id2_next.get("development_runs_cannot_enter_evidence") is not True
        or id2_boundary.get("dynamic_trajectory_used") is not False
        or id2_boundary.get("fresh_GR0_calibration_completed") is not False
        or id2_boundary.get("PRO4_HLD1_resolved") is not False
        or id2_boundary.get("mechanism_question_answered") is not False
        or id2_boundary.get("thresholds_changed_after_inspecting_static_results")
        is not False
        or id2_boundary.get("static_input_exclusion_is_not_FGCQR_mechanism_rejection")
        is not True
        or id2_boundary.get("SGBL_branch_health_definition_missing") is not True
        or id2_boundary.get("constraint_to_Raychaudhuri_stability_map_missing")
        is not True
    ):
        failures.append(
            "FGC-1-ID2-ALL1 Cartesian ledger, eligibility, hash, or nonclaim contract has drifted"
        )

    cal1_path = REPOSITORY / "results/fgc-1-cal1-pref3.json"
    try:
        cal1 = json.loads(
            cal1_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-CAL1-PREF3 certificate: {exc}")
        return failures
    cal1_payload = cal1.get("artifact_payload", {})
    cal1_diagnosis = cal1_payload.get("diagnosis", {})
    cal1_boundary = cal1_payload.get("epistemic_boundary", {})
    cal1_candidates = cal1_payload.get("candidate_compositions", [])
    cal1_methods = [
        (candidate.get("amplitude"), method)
        for candidate in cal1_candidates
        if isinstance(candidate, dict)
        for method in candidate.get("methods", [])
        if isinstance(method, dict)
    ]
    cal1_resolutions = [
        resolution
        for _amplitude, method in cal1_methods
        for resolution in method.get("resolutions", [])
        if isinstance(resolution, dict)
    ]
    cal1_state_hashes = [
        resolution.get(key)
        for resolution in cal1_resolutions
        for key in (
            "original_ID2_state_sha256",
            "projected_semidiscrete_state_sha256",
        )
    ]
    expected_cal1_methods = {
        ("5/2", "RK4", 4),
        ("5/2", "SSPRK3", 2),
        ("3", "RK4", 4),
        ("3", "SSPRK3", 2),
    }
    if (
        cal1.get("classification")
        != "pre_trajectory_semidiscrete_test_contract_obstruction_not_mechanism_result"
        or cal1.get("generated_by") != "scripts/reproduce_fgc_cal1_pref3.py"
        or cal1.get("scope_bindings", {}).get("target_protocol")
        != "FGC-2-SF1-PROTO4"
        or cal1.get("scope_bindings", {}).get("calibration_branch") != "GR-0"
        or cal1.get("scope_bindings", {}).get("eligible_amplitudes")
        != ["5/2", "3"]
        or cal1.get("scope_bindings", {}).get("methods") != ["RK4", "SSPRK3"]
        or cal1.get("scope_bindings", {}).get("resolutions")
        != [1025, 2049, 4097]
        or cal1.get("scope_bindings", {}).get("dynamic_output_namespace_accessed")
        is not False
        or cal1_diagnosis
        != {
            "exact_zero_only_order_rule_turns_binary64_roundoff_into_fake_nonconvergence": True,
            "missing_continuum_to_semidiscrete_q_map": True,
            "new_protocol_version_required": True,
            "physical_equations_or_candidate_outcome_implicated": False,
            "single_method_independent_constraint_guards_reject_declared_second_order_comparator_at_t0": True,
        }
        or cal1_boundary
        != {
            "FGCQR_outcome_read": False,
            "collapse_or_trapping_classified": False,
            "mechanism_question_answered": False,
            "trajectory_read": False,
        }
        or [candidate.get("amplitude") for candidate in cal1_candidates]
        != ["5/2", "3"]
        or len(cal1_methods) != 4
        or {
            (amplitude, method.get("method"), method.get("spatial_order"))
            for amplitude, method in cal1_methods
        }
        != expected_cal1_methods
        or any(
            method.get("original_ID2_state_under_literal_PROTO4", {}).get(
                "admission_passed"
            )
            is not False
            or method.get("projected_state_under_literal_PROTO4", {}).get(
                "admission_passed"
            )
            is not False
            or method.get("projected_state_under_prospective_repair", {}).get(
                "admission_passed"
            )
            is not True
            or method.get("projected_state_under_prospective_repair", {}).get(
                "common_event_alignment_passed"
            )
            is not True
            or method.get("projected_state_under_prospective_repair", {}).get(
                "monotone_refinement_passed"
            )
            is not True
            or method.get("projected_state_under_prospective_repair", {}).get(
                "finest_pair_order_passed"
            )
            is not True
            or method.get("projected_state_under_prospective_repair", {}).get(
                "point_counts"
            )
            != [1025, 2049, 4097]
            for _amplitude, method in cal1_methods
        )
        or len(cal1_resolutions) != 12
        or any(
            resolution.get("u_bitwise_preserved") is not True
            or resolution.get("p_bitwise_preserved") is not True
            or not isinstance(resolution.get("q_projection_change_infinity"), (int, float))
            or resolution.get("q_projection_change_infinity") <= 0.0
            or not isinstance(resolution.get("original_raw_global_constraint"), (int, float))
            or resolution.get("original_raw_global_constraint") < 0.0
            or not isinstance(resolution.get("projected_raw_global_constraint"), (int, float))
            or resolution.get("projected_raw_global_constraint") < 0.0
            or len(resolution.get("normalized_roundoff_zero_enclosure", {})) != 10
            or any(
                not isinstance(value, (int, float)) or value <= 0.0
                for value in resolution.get(
                    "normalized_roundoff_zero_enclosure", {}
                ).values()
            )
            or len(resolution.get("original_raw_component_constraints", {})) != 10
            or len(resolution.get("projected_raw_component_constraints", {})) != 10
            for resolution in cal1_resolutions
        )
        or len(cal1_state_hashes) != 24
        or len(set(cal1_state_hashes)) != 24
        or any(
            not isinstance(value, str)
            or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)
            for value in cal1_state_hashes
        )
        or any(value is not False for value in cal1.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-CAL1-PREF3 semidiscrete composition, prospective repair, hash, or nonclaim contract has drifted"
        )

    pro5_path = REPOSITORY / "results/fgc-1-pro5-frz1.json"
    try:
        pro5 = json.loads(
            pro5_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-PRO5-FRZ1 certificate: {exc}")
        return failures
    pro5_payload = pro5.get("artifact_payload", {})
    pro5_protocol = pro5_payload.get("protocol_certificate", {})
    pro5_lineage = pro5_payload.get("immutable_lineage", {})
    pro5_predecessor = pro5_lineage.get("predecessor_protocol", {})
    pro5_diagnosis = pro5_lineage.get("diagnosis", {})
    pro5_repair = pro5_payload.get("frozen_repair", {})
    pro5_boundary = pro5_payload.get("epistemic_boundary", {})
    if (
        pro5.get("classification")
        != "immutable_outcome_neutral_PROTO5_semidiscrete_protocol_freeze_not_run_authorization"
        or pro5.get("generated_by") != "scripts/reproduce_fgc_pro5_frz1.py"
        or pro5.get("scope_bindings", {}).get("target_protocol")
        != "FGC-2-SF1-PROTO5"
        or pro5.get("scope_bindings", {}).get("predecessor_protocol")
        != "FGC-2-SF1-PROTO4"
        or pro5.get("scope_bindings", {}).get("diagnosis_artifact")
        != "FGC-1-CAL1-PREF3"
        or pro5.get("scope_bindings", {}).get(
            "fresh_GR0_calibration_trajectory_inspected"
        )
        is not False
        or pro5_protocol.get("artifact_id") != "FGC-2-SF1-PROTO5"
        or pro5_protocol.get("protocol_version") != 5
        or pro5_protocol.get("frozen") is not True
        or pro5_protocol.get("outcome_neutral_contract_validated") is not True
        or pro5_protocol.get("eligible_calibration_amplitudes") != ["5/2", "3"]
        or pro5_protocol.get("native_semidiscrete_q_initialization_validated")
        is not True
        or pro5_protocol.get("physical_u_and_p_unchanged") is not True
        or pro5_protocol.get("method_owned_raw_constraint_guards_validated")
        is not True
        or pro5_protocol.get("roundoff_zero_classification_only_validated")
        is not True
        or pro5_protocol.get("monotone_three_grid_finest_pair_order_validated")
        is not True
        or pro5_protocol.get("resolved_holdout_manifest")
        != "configs/fgc/fgc-1-pro5-hld1.toml"
        or pro5_protocol.get("calibration_output_root")
        != "runs/fgc-2-sf1/proto5/calibration"
        or pro5_protocol.get("holdout_output_root")
        != "runs/fgc-2-sf1/proto5/holdout"
        or any(value is not False for value in pro5_protocol.get("claims", {}).values())
        or pro5_lineage.get("checkpoint_commit")
        != "0781d7d218e9c129858b134c3181fad07ac692b0"
        or pro5_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or pro5_predecessor.get("artifact_id") != "FGC-2-SF1-PROTO4"
        or pro5_predecessor.get("protocol_version") != 4
        or pro5_predecessor.get("claims_all_false") is not True
        or pro5_diagnosis.get("artifact_id") != "FGC-1-CAL1-PREF3"
        or pro5_diagnosis.get("PROTO4_obstruction_verified") is not True
        or pro5_diagnosis.get("premise_revision_required") is not True
        or pro5_diagnosis.get("eligible_amplitudes") != ["5/2", "3"]
        or pro5_diagnosis.get("prospective_t0_compositions_passed") != 4
        or pro5_diagnosis.get("trajectory_read") is not False
        or pro5_diagnosis.get("FGCQR_outcome_read") is not False
        or pro5_repair
        != {
            "all_three_resolutions_must_refine_monotonically": True,
            "comparator_guards": {"coarsest": "1/10", "finest": "1/200"},
            "eligible_amplitudes_in_order": ["5/2", "3"],
            "minimum_finest_pair_order": "3/2",
            "only_q_is_projected": True,
            "primary_guards": {"coarsest": "1/50", "finest": "1/1000"},
            "raw_norms_decide_magnitude_guards": True,
            "roundoff_enclosure_is_zero_classification_only": True,
            "roundoff_operation_budget": 4096,
            "u_and_p_must_be_bitwise_preserved": True,
        }
        or pro5_boundary
        != {
            "FGCQR_outcome_read": False,
            "SGBL_outcome_read": False,
            "fresh_calibration_trajectory_read": False,
            "initial_common_event_diagnostics_were_seen": True,
            "mechanism_question_answered": False,
            "physical_experiment_changed": False,
        }
        or any(value is not False for value in pro5.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-PRO5-FRZ1 lineage, semidiscrete repair, namespace, or nonclaim contract has drifted"
        )

    hlt3_path = REPOSITORY / "results/fgc-1-hlt3-mon3.json"
    try:
        hlt3 = json.loads(
            hlt3_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-HLT3-MON3 certificate: {exc}")
        return failures
    hlt3_payload = hlt3.get("artifact_payload", {})
    hlt3_lineage = hlt3_payload.get("immutable_lineage", {})
    hlt3_runtime = hlt3_payload.get("runtime_composition", {})
    hlt3_temporal = hlt3_payload.get("temporal_spectral_controls", {})
    hlt3_plan = hlt3_payload.get("frozen_run_plan", {})
    hlt3_inputs = hlt3_payload.get("frozen_run_inputs", [])
    hlt3_events = hlt3_payload.get("initial_common_events", [])
    hlt3_namespace = hlt3_payload.get("namespace_precondition", {})
    hlt3_boundary = hlt3_payload.get("epistemic_boundary", {})
    expected_hlt3_stop_ids = [
        "nonpositive_lapse",
        "nonpositive_radial_metric",
        "nonpositive_areal_radius_away_from_center",
        "hat_cone_not_Lorentzian",
        "newton_residual_limit",
        "newton_iteration_limit",
        "newton_residual_not_monotonic",
        "kinetic_condition_limit",
        "boundary_causal_buffer",
    ]
    expected_hlt3_inputs = [
        (amplitude, method, point_count)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
        for point_count in (1025, 2049, 4097)
    ]
    observed_hlt3_inputs = [
        (item.get("amplitude"), item.get("method"), item.get("point_count"))
        for item in hlt3_inputs
        if isinstance(item, dict)
    ]
    hlt3_hash_values = [
        item.get(name)
        for item in hlt3_inputs
        if isinstance(item, dict)
        for name in (
            "continuum_ID2_state_sha256",
            "projected_state_sha256",
            "expanded_run_config_sha256",
        )
    ]
    expected_hlt3_events = [
        (amplitude, method, 0.0)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
    ]
    observed_hlt3_events = [
        (item.get("amplitude"), item.get("method"), item.get("coordinate_time"))
        for item in hlt3_events
        if isinstance(item, dict)
    ]
    hlt3_namespace_records = hlt3_namespace.get("records", [])
    hlt3_cfl = hlt3_runtime.get("CFL_retry_control", {})
    hlt3_recovery = hlt3_runtime.get("event_log_recovery_controls", {})
    expected_hlt3_implementation = {
        "scripts/reproduce_fgc_hlt3_mon3.py",
        "scripts/run_fgc_gr0_calibration.py",
        "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto5_runtime.py",
    }
    if (
        hlt3.get("classification")
        != "pre_trajectory_PROTO5_GR0_runtime_compositor_and_fresh_calibration_authorization"
        or hlt3.get("generated_by") != "scripts/reproduce_fgc_hlt3_mon3.py"
        or hlt3.get("scope_bindings", {}).get("target_protocol")
        != "FGC-2-SF1-PROTO5"
        or hlt3.get("scope_bindings", {}).get("calibration_branch") != "GR-0"
        or hlt3.get("scope_bindings", {}).get("fresh_GR0_trajectory_read")
        is not False
        or hlt3.get("scope_bindings", {}).get("FGCQR_trajectory_read") is not False
        or hlt3_lineage.get("checkpoint_commit")
        != "cf48cf729ff1c3a6effce6ecd3fdbeaf8d9645bf"
        or hlt3_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or hlt3_lineage.get("freeze_gate_passed") is not True
        or hlt3_lineage.get("diagnosis_is_pre_trajectory") is not True
        or hlt3_lineage.get("eligible_amplitudes") != ["5/2", "3"]
        or hlt3_lineage.get("HLT2_monitor_gate_passed") is not True
        or hlt3_runtime.get("universal_stop_ids") != expected_hlt3_stop_ids
        or hlt3_runtime.get("candidate_action_stop_count") != 0
        or hlt3_runtime.get("all_controls_passed") is not True
        or hlt3_runtime.get("late_failure_rolled_back") is not True
        or hlt3_runtime.get("first_scientific_failure_immutable") is not True
        or hlt3_runtime.get("boundary_failure_rolled_back") is not True
        or hlt3_recovery.get("actions")
        != {
            "absent": "restored_missing_log_from_checkpoint",
            "ahead": "archived_uncommitted_tail_and_restored_checkpoint",
            "behind": "completed_checkpoint_owned_log_write",
            "matching": "event_log_already_matches_checkpoint",
        }
        or hlt3_recovery.get("divergent_history_rejected") is not True
        or hlt3_recovery.get("all_recovery_controls_passed") is not True
        or hlt3_cfl.get("retry_required") is not True
        or hlt3_cfl.get("scientific_failure_latched") is not False
        or hlt3_cfl.get("causal_or_stage_acceptance_advanced") is not False
        or hlt3_temporal.get("minimum_sample_count") != 64
        or hlt3_temporal.get("low_frequency_control_passed") is not True
        or hlt3_temporal.get("alias_injection_rejected") is not True
        or hlt3_plan.get("artifact_id") != "FGC-1-CAL2-RUN1-PLAN"
        or hlt3_plan.get("campaign_id")
        != "9df7046b47a044df920bfe192e6b63fbda7953ab97d9174850927afec10177be"
        or hlt3_plan.get("ordered_amplitudes") != ["5/2", "3"]
        or hlt3_plan.get("expanded_input_count") != 12
        or hlt3_plan.get("first_eligible_candidate_stops_later_candidates")
        is not True
        or observed_hlt3_inputs != expected_hlt3_inputs
        or any(
            item.get("trajectory_advanced") is not False
            or item.get("u_bitwise_preserved") is not True
            or item.get("p_bitwise_preserved") is not True
            or item.get("q_initialized_by_native_SBP_operator") is not True
            for item in hlt3_inputs
            if isinstance(item, dict)
        )
        or len(hlt3_hash_values) != 36
        or any(
            not isinstance(value, str)
            or len(value) != 64
            or any(character not in "0123456789abcdef" for character in value)
            for value in hlt3_hash_values
        )
        or observed_hlt3_events != expected_hlt3_events
        or any(
            item.get("admission_passed") is not True
            or item.get("trajectory_advanced") is not False
            for item in hlt3_events
            if isinstance(item, dict)
        )
        or hlt3_namespace.get("both_fresh_namespaces_empty") is not True
        or hlt3_namespace.get("authorization_created_no_namespace") is not True
        or len(hlt3_namespace_records) != 2
        or any(
            item.get("empty") is not True
            or item.get("entry_count") != 0
            or item.get("created_by_authorization") is not False
            for item in hlt3_namespace_records
            if isinstance(item, dict)
        )
        or hlt3_payload.get("decision")
        != "HLT3_implements_the_PROTO5_GR0_runtime_contract_and_authorizes_only_the_fresh_calibration_campaign"
        or hlt3_boundary.get("fresh_GR0_trajectory_read") is not False
        or hlt3_boundary.get("fresh_GR0_trajectory_advanced") is not False
        or hlt3_boundary.get("fresh_GR0_calibration_is_now_authorized_but_not_completed")
        is not True
        or hlt3_boundary.get("trapped_or_collapse_outcome_classified") is not False
        or hlt3_boundary.get("SGBL_or_FGCQR_outcome_read") is not False
        or hlt3_boundary.get("resolved_holdout_manifest_exists") is not False
        or hlt3_boundary.get("mechanism_question_answered") is not False
        or set(hlt3.get("implementation_sha256", {}))
        != expected_hlt3_implementation
        or any(value is not False for value in hlt3.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-HLT3-MON3 runtime, input-freeze, namespace, recovery, or nonclaim contract has drifted"
        )

    cal2_path = REPOSITORY / "results/fgc-1-cal2-pref4.json"
    try:
        cal2 = json.loads(
            cal2_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-CAL2-PREF4 certificate: {exc}")
        return failures
    cal2_payload = cal2.get("artifact_payload", {})
    cal2_stops = cal2_payload.get("amplitude_stop_evidence", [])
    cal2_replay = cal2_payload.get("diagnostic_replay", {})
    cal2_accepted = cal2_replay.get("accepted_state_source_evidence", {})
    cal2_proposals = cal2_replay.get("proposals", [])
    cal2_decision = cal2_payload.get("decision", {})
    cal2_repair = cal2_decision.get("prospective_repair", {})
    cal2_boundary = cal2_payload.get("epistemic_boundary", {})
    cal2_campaign = cal2_payload.get("immutable_campaign", {})
    expected_cal2_raw_hashes = {
        "runs/fgc-2-sf1/proto5/calibration/campaign-result.json":
            "9cbd2cdd44a85bb7c3d65b2c13c7b6f192e4d38375f8e4efc642208717be2296",
        "runs/fgc-2-sf1/proto5/calibration/events.jsonl":
            "022ca3629134eb710df931aeef03ccba91c288856473696bf7a6ddeaf9e1b2d1",
        "runs/fgc-2-sf1/proto5/calibration/latest-checkpoint.npz":
            "8607800ada752c7571b00e975399b8c038c145fa8635b3ffaedf4a010c1e7f12",
        "runs/fgc-2-sf1/proto5/calibration/manifest.json":
            "f47c9c54c9ddccb9c1b97807512c3ed1ae93f6692030e0d583ac8b71c6ba55d0",
    }
    expected_cal2_repair = {
        "accepted_state_source_failure_remains_terminal": True,
        "amplitude_order_unchanged": True,
        "candidate_action_health_stops_unchanged": True,
        "constraint_spectral_trapped_and_boundary_rules_unchanged": True,
        "fresh_output_namespace_and_new_protocol_version_required": True,
        "maximum_retries_inherited": 32,
        "minimum_step_size_inherited": "1/1073741824",
        "only_internal_unaccepted_stage_source_failure_becomes_retryable": True,
        "physical_inputs_unchanged": True,
        "raw_source_residual_tolerance_unchanged": True,
        "retry_factor_inherited": "1/2",
    }
    expected_cal2_implementation = {
        "scripts/reproduce_fgc_cal2_pref4.py",
        "src/recursive_horizons/fgc/evolution/cal2_source_diagnosis.py",
    }
    observed_cal2_stops = [
        (
            item.get("amplitude"),
            item.get("reason"),
            item.get("last_fully_completed_common_event_index"),
            item.get("observed_residual_infinity"),
            item.get("raw_residual_tolerance"),
        )
        for item in cal2_stops
        if isinstance(item, dict)
    ]
    observed_cal2_proposals = [
        (
            item.get("step_factor"),
            item.get("all_stage_raw_source_gates_passed"),
        )
        for item in cal2_proposals
        if isinstance(item, dict)
    ]
    cal2_proposal_conditions = [
        item.get("maximum_kinetic_condition_infinity")
        for item in cal2_proposals
        if isinstance(item, dict)
    ]
    cal2_proposal_corrections = [
        item.get("maximum_relative_acceleration_correction")
        for item in cal2_proposals
        if isinstance(item, dict)
    ]
    if (
        cal2.get("classification")
        != "post_calibration_numerical_trial_stop_ownership_obstruction_not_mechanism_result"
        or cal2.get("generated_by") != "scripts/reproduce_fgc_cal2_pref4.py"
        or cal2.get("scope_bindings")
        != {
            "FGCQR_trajectory_read": False,
            "GR0_calibration_trajectory_read": True,
            "SGBL_trajectory_read": False,
            "calibration_branch": "GR-0",
            "first_nonzero_common_event_completed": False,
            "target_protocol": "FGC-2-SF1-PROTO5",
        }
        or cal2.get("raw_campaign_sha256") != expected_cal2_raw_hashes
        or observed_cal2_stops
        != [
            ("5/2", "newton_residual_limit", 0, 1.0511642392237476e-12, 1e-12),
            ("3", "newton_residual_limit", 0, 1.0511640449698101e-12, 1e-12),
        ]
        or cal2_accepted.get("raw_source_gate_passed") is not True
        or not isinstance(cal2_accepted.get("raw_source_residual_infinity"), (int, float))
        or cal2_accepted.get("raw_source_residual_infinity") >= 1e-12
        or cal2_replay.get("failed_unit_proposal_stage_names") != ["rk4_k4"]
        or cal2_replay.get("first_retry_factor_that_passes_every_stage") != "1/2"
        or cal2_replay.get("failure_localizes_to_first_annular_node") is not True
        or observed_cal2_proposals
        != [("1", False), ("1/2", True), ("1/4", True), ("1/8", True)]
        or len(cal2_proposal_conditions) != 4
        or any(
            not isinstance(value, (int, float)) or value >= 1e10
            for value in cal2_proposal_conditions
        )
        or len(cal2_proposal_corrections) != 4
        or any(
            not isinstance(value, (int, float)) or value >= 1e-12
            for value in cal2_proposal_corrections
        )
        or cal2_decision.get("PROTO5_GR0_source_stop_ownership_contract_obstructed")
        is not True
        or cal2_decision.get("new_protocol_required") != "FGC-2-SF1-PROTO6"
        or cal2_repair != expected_cal2_repair
        or cal2_campaign.get("terminal_classification")
        != "calibration_failed_no_eligible_GR0_case"
        or cal2_campaign.get("event_log_line_count") != 2
        or cal2_campaign.get("manifest", {}).get("git_commit")
        != "fb7d1dc2499438ac639e6a1686cd3f5561fa93bc"
        or cal2_campaign.get("terminal_checkpoint", {}).get(
            "completed_common_event_index"
        )
        != 0
        or cal2_boundary
        != {
            "GR0_amplitude_selected": False,
            "GR0_calibration_completed": False,
            "SGBL_or_FGCQR_trajectory_read": False,
            "candidate_or_general_gradient_route_rejected": False,
            "dynamic_trapped_sphere_classified": False,
            "mechanism_question_answered": False,
        }
        or set(cal2.get("implementation_sha256", {}))
        != expected_cal2_implementation
        or any(value is not False for value in cal2.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-CAL2-PREF4 campaign, replay, repair, or nonclaim contract has drifted"
        )

    pro6_path = REPOSITORY / "results/fgc-1-pro6-frz1.json"
    try:
        pro6 = json.loads(
            pro6_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-PRO6-FRZ1 certificate: {exc}")
        return failures
    pro6_payload = pro6.get("artifact_payload", {})
    pro6_protocol = pro6_payload.get("protocol_certificate", {})
    pro6_lineage = pro6_payload.get("immutable_lineage", {})
    pro6_predecessor = pro6_lineage.get("predecessor_protocol", {})
    pro6_diagnosis = pro6_lineage.get("diagnosis", {})
    pro6_repair = pro6_payload.get("frozen_repair", {})
    pro6_boundary = pro6_payload.get("epistemic_boundary", {})
    expected_pro6_repair = {
        "accepted_state_source_failure_is_terminal": True,
        "all_non_source_stops_are_unchanged": True,
        "eligible_amplitudes_in_order": ["5/2", "3"],
        "maximum_retry_count": 32,
        "minimum_step_size": "1/1073741824",
        "only_internal_unaccepted_stage_source_failure_is_retryable": True,
        "raw_source_residual_tolerance_is_unchanged": True,
        "retry_factor": "1/2",
        "retry_restarts_from_bitwise_last_accepted_state": True,
        "source_refinement_algorithm_is_unchanged": True,
    }
    expected_pro6_implementation = {
        "scripts/reproduce_fgc_pro6_frz1.py",
        "src/recursive_horizons/fgc/evolution/protocol_v5.py",
        "src/recursive_horizons/fgc/evolution/protocol_v6.py",
    }
    if (
        pro6.get("classification")
        != "immutable_outcome_neutral_PROTO6_source_stop_ownership_protocol_freeze_not_run_authorization"
        or pro6.get("generated_by") != "scripts/reproduce_fgc_pro6_frz1.py"
        or pro6.get("scope_bindings", {}).get("target_protocol")
        != "FGC-2-SF1-PROTO6"
        or pro6.get("scope_bindings", {}).get("predecessor_protocol")
        != "FGC-2-SF1-PROTO5"
        or pro6.get("scope_bindings", {}).get("diagnosis_artifact")
        != "FGC-1-CAL2-PREF4"
        or pro6.get("scope_bindings", {}).get(
            "PROTO5_GR0_calibration_trajectory_disclosed"
        )
        is not True
        or pro6.get("scope_bindings", {}).get("first_nonzero_common_event_completed")
        is not False
        or pro6.get("scope_bindings", {}).get("FGCQR_evolution_outcomes_inspected")
        is not False
        or pro6.get("scope_bindings", {}).get("mechanism_question_answered")
        is not False
        or pro6_protocol.get("artifact_id") != "FGC-2-SF1-PROTO6"
        or pro6_protocol.get("protocol_version") != 6
        or pro6_protocol.get("frozen") is not True
        or pro6_protocol.get("outcome_neutral_contract_validated") is not True
        or pro6_protocol.get("accepted_state_source_failure_terminal") is not True
        or pro6_protocol.get("internal_unaccepted_stage_source_failure_retryable")
        is not True
        or pro6_protocol.get("raw_source_residual_tolerance_unchanged") is not True
        or pro6_protocol.get("all_non_source_stops_unchanged") is not True
        or pro6_protocol.get("retry_factor") != "1/2"
        or pro6_protocol.get("maximum_retry_count") != 32
        or pro6_protocol.get("minimum_step_size") != "1/1073741824"
        or pro6_protocol.get("resolved_holdout_manifest")
        != "configs/fgc/fgc-1-pro6-hld1.toml"
        or pro6_protocol.get("calibration_output_root")
        != "runs/fgc-2-sf1/proto6/calibration"
        or pro6_protocol.get("holdout_output_root")
        != "runs/fgc-2-sf1/proto6/holdout"
        or any(value is not False for value in pro6_protocol.get("claims", {}).values())
        or pro6_lineage.get("checkpoint_commit")
        != "47f4c83bae2df26c51c3544eed8c06507a0cd9eb"
        or pro6_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or pro6_predecessor.get("artifact_id") != "FGC-2-SF1-PROTO5"
        or pro6_predecessor.get("protocol_version") != 5
        or pro6_predecessor.get("claims_all_false") is not True
        or pro6_diagnosis.get("artifact_id") != "FGC-1-CAL2-PREF4"
        or pro6_diagnosis.get("PROTO5_source_stop_ownership_obstructed") is not True
        or pro6_diagnosis.get("premise_revision_required") is not True
        or pro6_diagnosis.get("accepted_state_source_gate_passed") is not True
        or pro6_diagnosis.get("unit_proposal_source_gate_missed") is not True
        or pro6_diagnosis.get("first_retry_factor_passed") != "1/2"
        or pro6_diagnosis.get("first_nonzero_common_event_completed") is not False
        or pro6_diagnosis.get("FGCQR_outcome_read") is not False
        or pro6_diagnosis.get("mechanism_question_answered") is not False
        or pro6_repair != expected_pro6_repair
        or pro6_boundary
        != {
            "FGCQR_outcome_read": False,
            "PROTO5_GR0_calibration_trajectory_was_seen": True,
            "PROTO6_trajectory_read": False,
            "SGBL_outcome_read": False,
            "first_nonzero_common_event_completed": False,
            "mechanism_question_answered": False,
            "physical_experiment_changed": False,
        }
        or pro6_payload.get("decision")
        != "PROTO6_is_frozen_but_HLT4_fresh_calibration_and_every_candidate_outcome_remain_closed"
        or set(pro6.get("implementation_sha256", {}))
        != expected_pro6_implementation
        or any(value is not False for value in pro6.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-PRO6-FRZ1 lineage, source-stop ownership, namespace, or nonclaim contract has drifted"
        )

    hlt4_path = REPOSITORY / "results/fgc-1-hlt4-mon4.json"
    try:
        hlt4 = json.loads(
            hlt4_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-HLT4-MON4 certificate: {exc}")
        return failures
    hlt4_payload = hlt4.get("artifact_payload", {})
    hlt4_lineage = hlt4_payload.get("immutable_lineage", {})
    hlt4_runtime = hlt4_payload.get("runtime_controls", {})
    hlt4_internal = hlt4_runtime.get("internal_trial_control", {})
    hlt4_accepted = hlt4_runtime.get("accepted_state_control", {})
    hlt4_endpoint = hlt4_runtime.get("candidate_endpoint_control", {})
    hlt4_simultaneous = hlt4_runtime.get("simultaneous_failure_control", {})
    hlt4_healthy = hlt4_runtime.get("healthy_control", {})
    hlt4_rollback = hlt4_runtime.get("cross_member_rollback_control", {})
    hlt4_equivalence = hlt4_payload.get("accepted_source_operator_equivalence", {})
    hlt4_inputs = hlt4_payload.get("frozen_run_inputs", [])
    hlt4_prechecks = hlt4_payload.get("accepted_state_source_prechecks", [])
    hlt4_events = hlt4_payload.get("initial_common_events", [])
    hlt4_plan = hlt4_payload.get("frozen_run_plan", {})
    hlt4_namespace = hlt4_payload.get("namespace_precondition", {})
    hlt4_boundary = hlt4_payload.get("epistemic_boundary", {})
    expected_hlt4_order = [
        (amplitude, method, point_count)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
        for point_count in (1025, 2049, 4097)
    ]
    observed_hlt4_order = [
        (item.get("amplitude"), item.get("method"), item.get("point_count"))
        for item in hlt4_inputs
        if isinstance(item, dict)
    ]
    observed_hlt4_prechecks = [
        (item.get("amplitude"), item.get("method"), item.get("point_count"))
        for item in hlt4_prechecks
        if isinstance(item, dict)
    ]
    expected_hlt4_events = [
        (amplitude, method, 0.0)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
    ]
    observed_hlt4_events = [
        (item.get("amplitude"), item.get("method"), item.get("coordinate_time"))
        for item in hlt4_events
        if isinstance(item, dict)
    ]
    expected_hlt4_implementation = {
        "scripts/reproduce_fgc_hlt4_mon4.py",
        "scripts/run_fgc_gr0_calibration.py",
        "scripts/run_fgc_gr0_calibration_v6.py",
        "src/recursive_horizons/fgc/evolution/cal2_source_diagnosis.py",
        "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto5_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto6_runtime.py",
    }
    if (
        hlt4.get("classification")
        != "pre_trajectory_PROTO6_GR0_source_retry_runtime_and_fresh_calibration_authorization"
        or hlt4.get("generated_by") != "scripts/reproduce_fgc_hlt4_mon4.py"
        or hlt4.get("scope_bindings", {}).get("target_protocol")
        != "FGC-2-SF1-PROTO6"
        or hlt4.get("scope_bindings", {}).get("calibration_branch") != "GR-0"
        or hlt4.get("scope_bindings", {}).get("PROTO6_trajectory_read") is not False
        or hlt4.get("scope_bindings", {}).get("FGCQR_trajectory_read") is not False
        or hlt4_lineage.get("checkpoint_commit")
        != "f3a32ff9faddfa4dba0e0f5e192d58ef94d5e094"
        or hlt4_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or hlt4_lineage.get("protocol_freeze_gate_passed") is not True
        or hlt4_lineage.get("predecessor_runtime_gate_passed") is not True
        or hlt4_lineage.get("static_input_gate_passed") is not True
        or hlt4_runtime.get("all_controls_passed") is not True
        or hlt4_accepted.get("terminal_reason") != "newton_residual_limit"
        or hlt4_accepted.get("first_failure_latched") is not True
        or hlt4_accepted.get("accepted_stage_count") != 0
        or hlt4_internal.get("failed_stage_name") != "rk4_k4"
        or hlt4_internal.get("state_unchanged") is not True
        or hlt4_internal.get("causal_debit_unchanged") is not True
        or hlt4_internal.get("preaccept_not_called") is not True
        or hlt4_endpoint.get("terminal_reason") != "newton_residual_limit"
        or hlt4_endpoint.get("first_failure_latched") is not True
        or hlt4_simultaneous.get("reason") != "nonpositive_lapse"
        or "newton_residual_limit" not in hlt4_simultaneous.get("failures", [])
        or hlt4_healthy.get("accepted") is not True
        or hlt4_healthy.get("accepted_stage_count") != 5
        or hlt4_healthy.get("preaccept_observed_before_commit") is not True
        or hlt4_rollback.get("restored_value") != 3
        or hlt4_rollback.get("incomplete_key_set_rejected") is not True
        or hlt4_equivalence.get("all_controls_passed") is not True
        or len(hlt4_equivalence.get("controls", [])) != 2
        or any(
            item.get("binary64_RHS_arrays_identical") is not True
            for item in hlt4_equivalence.get("controls", [])
            if isinstance(item, dict)
        )
        or observed_hlt4_order != expected_hlt4_order
        or any(
            item.get("trajectory_advanced") is not False
            or item.get("source_stop_ownership")
            != "accepted_terminal_internal_source_only_trial_retry"
            for item in hlt4_inputs
            if isinstance(item, dict)
        )
        or observed_hlt4_prechecks != expected_hlt4_order
        or any(
            item.get("accepted_state_source_gate_passed") is not True
            or item.get("trajectory_advanced") is not False
            or not isinstance(item.get("source_residual_infinity"), (int, float))
            or item.get("source_residual_infinity")
            >= item.get("source_residual_maximum", 0)
            for item in hlt4_prechecks
            if isinstance(item, dict)
        )
        or observed_hlt4_events != expected_hlt4_events
        or any(
            item.get("admission_passed") is not True
            or item.get("trajectory_advanced") is not False
            for item in hlt4_events
            if isinstance(item, dict)
        )
        or hlt4_plan.get("artifact_id") != "FGC-1-CAL3-RUN1-PLAN"
        or hlt4_plan.get("campaign_id")
        != "62a0f63e0c381da444890bdb5b2fe7bfeeb2c9b8ba84b07a3f9f77669af6f1b8"
        or hlt4_plan.get("ordered_amplitudes") != ["5/2", "3"]
        or hlt4_plan.get("expanded_input_count") != 12
        or hlt4_namespace.get("both_fresh_namespaces_empty") is not True
        or hlt4_namespace.get("authorization_created_no_namespace") is not True
        or hlt4_payload.get("decision")
        != "HLT4_implements_only_PROTO6_source_retry_ownership_and_authorizes_one_fresh_GR0_recalibration"
        or hlt4_boundary.get("PROTO6_trajectory_read") is not False
        or hlt4_boundary.get("PROTO6_trajectory_advanced") is not False
        or hlt4_boundary.get("fresh_GR0_recalibration_is_now_authorized_but_not_completed")
        is not True
        or hlt4_boundary.get("trapped_or_collapse_outcome_classified") is not False
        or hlt4_boundary.get("SGBL_or_FGCQR_outcome_read") is not False
        or hlt4_boundary.get("mechanism_question_answered") is not False
        or set(hlt4.get("implementation_sha256", {}))
        != expected_hlt4_implementation
        or any(value is not False for value in hlt4.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-HLT4-MON4 runtime ownership, input freeze, namespace, or nonclaim contract has drifted"
        )

    cal3_path = REPOSITORY / "results/fgc-1-cal3-pref5.json"
    try:
        cal3 = json.loads(
            cal3_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-CAL3-PREF5 certificate: {exc}")
        return failures
    cal3_payload = cal3.get("artifact_payload", {})
    cal3_campaign = cal3_payload.get("immutable_campaign", {})
    cal3_stops = cal3_campaign.get("amplitude_stop_records", [])
    cal3_replay = cal3_payload.get("deterministic_terminal_replay", {})
    cal3_terminal = cal3_replay.get("terminal_evidence", {})
    cal3_terminal_calls = cal3_terminal.get("calls", [])
    cal3_smaller = cal3_replay.get("one_further_halving_control", {})
    cal3_smaller_calls = cal3_smaller.get("calls", [])
    cal3_decision = cal3_payload.get("decision", {})
    cal3_repair = cal3_decision.get("prospective_repair", {})
    cal3_boundary = cal3_payload.get("epistemic_boundary", {})
    expected_cal3_repair = {
        "accepted_state_source_failure_remains_terminal": True,
        "amplitude_order_unchanged": True,
        "any_non_source_failure_prevents_retry": True,
        "candidate_endpoint_is_unaccepted_until_transaction_commit": True,
        "constraint_spectral_trapped_boundary_and_physical_rules_unchanged": True,
        "every_source_only_failure_in_an_unaccepted_proposal_is_retryable": True,
        "fresh_output_namespace_and_new_protocol_version_required": True,
        "maximum_retries_inherited": 32,
        "minimum_step_size_inherited": "1/1073741824",
        "physical_inputs_unchanged": True,
        "raw_source_residual_tolerance_unchanged": True,
        "retry_factor_inherited": "1/2",
        "terminal_and_retry_exhaustion_evidence_must_name_member_stage_and_residual": True,
    }
    expected_cal3_campaign_hashes = {
        "campaign_result_sha256":
            "fb50beb09821763964445f952a080215f7a8ba02a7e6dc70b3a0d4b76cfd3397",
        "checkpoint_sha256":
            "b5c0220ba5373cc182a852b071b724dc425af1c9914741e9ecd67eb4df34c9d2",
        "event_log_sha256":
            "a4885ebe857f0b395012d0a09552cad775c670686c1d146c80b93afa97ec2bef",
        "manifest_sha256":
            "8a9f92748e677c6a768ca61f38fd1390d6704d25ceb9421e6744e1a8b8490104",
    }
    observed_cal3_stops = [
        (
            item.get("amplitude"),
            item.get("eligible"),
            item.get("last_completed_common_event_index"),
            item.get("stop", {}).get("reason"),
            item.get("stop", {}).get("cross_member_state_rolled_back"),
        )
        for item in cal3_stops
        if isinstance(item, dict)
    ]
    observed_cal3_terminal_calls = [
        (
            item.get("stage_name"),
            item.get("source_raw_gate_passed"),
            item.get("source_residual_infinity"),
        )
        for item in cal3_terminal_calls
        if isinstance(item, dict)
    ]
    observed_cal3_smaller_calls = [
        (item.get("stage_name"), item.get("source_raw_gate_passed"))
        for item in cal3_smaller_calls
        if isinstance(item, dict)
    ]
    if (
        cal3.get("classification")
        != "post_PROTO6_GR0_unaccepted_candidate_endpoint_source_stop_ownership_diagnosis"
        or cal3.get("generated_by") != "scripts/reproduce_fgc_cal3_pref5.py"
        or cal3.get("scope_bindings")
        != {
            "FGCQR_trajectory_read": False,
            "GR0_calibration_trajectory_read": True,
            "SGBL_trajectory_read": False,
            "calibration_branch": "GR-0",
            "collapse_or_trapped_outcome_classified": False,
            "first_nonzero_common_event_completed": False,
            "mechanism_question_answered": False,
            "role": "post_calibration_unaccepted_candidate_endpoint_source_stop_ownership_diagnosis",
            "target_protocol": "FGC-2-SF1-PROTO6",
        }
        or cal3_campaign.get("authorization_commit")
        != "f8412dd202886dc60650acb887721be3fc62a737"
        or cal3_campaign.get("campaign_id")
        != "62a0f63e0c381da444890bdb5b2fe7bfeeb2c9b8ba84b07a3f9f77669af6f1b8"
        or any(
            cal3_campaign.get(key) != value
            for key, value in expected_cal3_campaign_hashes.items()
        )
        or cal3_campaign.get("terminal_classification")
        != "calibration_failed_no_eligible_GR0_case"
        or cal3_campaign.get("event_log_line_count") != 4
        or cal3_campaign.get("terminal_checkpoint_common_event_index") != 0
        or cal3_campaign.get("cross_member_rollback_proved") is not True
        or observed_cal3_stops
        != [
            ("5/2", False, 0, "newton_residual_limit", True),
            ("3", False, 0, "newton_residual_limit", True),
        ]
        or cal3_replay.get("accepted_attempt_count") != 36
        or cal3_replay.get("CFL_retry_count") != 37
        or cal3_replay.get("source_retry_count") != 1
        or cal3_replay.get("public_retry_record_reproduced_exactly") is not True
        or cal3_replay.get("last_accepted_state_passes_raw_gate") is not True
        or cal3_replay.get("all_terminal_failures_are_source_only") is not True
        or cal3_replay.get("failed_terminal_stages")
        != ["rk4_k4", "candidate_endpoint"]
        or cal3_replay.get("maximum_terminal_kinetic_condition")
        != 2835841.6197018796
        or cal3_replay.get("maximum_terminal_relative_acceleration_correction")
        != 4.1208681855221474e-14
        or cal3_terminal.get("reason") != "newton_residual_limit"
        or cal3_terminal.get("last_accepted_time") != 0.058593708693194996
        or cal3_terminal.get("last_accepted_transaction_serial") != 180
        or cal3_terminal.get("terminal_trial_step_size")
        != 0.0008138005612409849
        or observed_cal3_terminal_calls
        != [
            ("accepted_state_source_precheck", True, 5.370340537498194e-13),
            ("rk4_k1", True, 5.370340537498194e-13),
            ("rk4_k2", True, 4.59448208096713e-13),
            ("rk4_k3", True, 4.593694160895005e-13),
            ("rk4_k4", False, 1.3639071551646946e-12),
            ("candidate_endpoint", False, 1.3639179325880144e-12),
        ]
        or cal3_smaller.get("step_size") != 0.00040690028062049247
        or cal3_smaller.get("accepted") is not True
        or cal3_smaller.get("retry_returned") is not False
        or observed_cal3_smaller_calls
        != [
            ("accepted_state_source_precheck", True),
            ("rk4_k1", True),
            ("rk4_k2", True),
            ("rk4_k3", True),
            ("rk4_k4", True),
            ("candidate_endpoint", True),
        ]
        or cal3_decision.get(
            "PROTO6_unaccepted_candidate_endpoint_stop_ownership_contract_obstructed"
        )
        is not True
        or cal3_decision.get("new_protocol_required") != "FGC-2-SF1-PROTO7"
        or cal3_repair != expected_cal3_repair
        or cal3_boundary
        != {
            "GR0_amplitude_selected": False,
            "GR0_calibration_completed": False,
            "SGBL_or_FGCQR_trajectory_read": False,
            "candidate_or_general_gradient_route_rejected": False,
            "dynamic_trapped_sphere_classified": False,
            "first_nonzero_common_event_completed": False,
            "mechanism_question_answered": False,
        }
        or set(cal3.get("implementation_sha256", {}))
        != {"scripts/reproduce_fgc_cal3_pref5.py"}
        or any(value is not False for value in cal3.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-CAL3-PREF5 campaign, terminal replay, proposed repair, or nonclaim contract has drifted"
        )

    pro7_path = REPOSITORY / "results/fgc-1-pro7-frz1.json"
    try:
        pro7 = json.loads(
            pro7_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-PRO7-FRZ1 certificate: {exc}")
        return failures
    pro7_payload = pro7.get("artifact_payload", {})
    pro7_protocol = pro7_payload.get("protocol_certificate", {})
    pro7_lineage = pro7_payload.get("immutable_lineage", {})
    pro7_predecessor = pro7_lineage.get("predecessor_protocol", {})
    pro7_diagnosis = pro7_lineage.get("diagnosis", {})
    pro7_repair = pro7_payload.get("frozen_repair", {})
    pro7_boundary = pro7_payload.get("epistemic_boundary", {})
    expected_pro7_repair = {
        "accepted_state_source_failure_is_terminal": True,
        "all_non_source_stops_are_unchanged": True,
        "all_source_only_unaccepted_proposal_failures_are_retryable": True,
        "any_non_source_failure_vetoes_retry": True,
        "candidate_endpoint_is_unaccepted_until_transaction_commit": True,
        "complete_proposal_is_previewed_before_classification": True,
        "eligible_amplitudes_in_order": ["5/2", "3"],
        "maximum_retry_count": 32,
        "minimum_step_size": "1/1073741824",
        "raw_source_residual_tolerance_is_unchanged": True,
        "retry_factor": "1/2",
        "retry_restarts_from_bitwise_last_accepted_state": True,
        "source_refinement_algorithm_is_unchanged": True,
        "terminal_and_retry_exhaustion_evidence_is_complete": True,
    }
    expected_pro7_implementation = {
        "scripts/reproduce_fgc_pro7_frz1.py",
        "src/recursive_horizons/fgc/evolution/protocol_v6.py",
        "src/recursive_horizons/fgc/evolution/protocol_v7.py",
    }
    if (
        pro7.get("classification")
        != "immutable_outcome_neutral_PROTO7_general_transaction_stop_ownership_protocol_freeze_not_run_authorization"
        or pro7.get("generated_by") != "scripts/reproduce_fgc_pro7_frz1.py"
        or pro7.get("scope_bindings")
        != {
            "FGCQR_evolution_outcomes_inspected": False,
            "PROTO6_GR0_calibration_trajectory_disclosed": True,
            "SGBL_evolution_outcomes_inspected": False,
            "calibration_branch": "GR-0",
            "diagnosis_artifact": "FGC-1-CAL3-PREF5",
            "first_nonzero_common_event_completed": False,
            "freeze_role": "premise_only_post_calibration_general_unaccepted_proposal_ownership_protocol_and_lineage_certificate",
            "holdout_branch": "FGC-QR",
            "mechanism_question_answered": False,
            "predecessor_protocol": "FGC-2-SF1-PROTO6",
            "target_protocol": "FGC-2-SF1-PROTO7",
        }
        or pro7_protocol.get("artifact_id") != "FGC-2-SF1-PROTO7"
        or pro7_protocol.get("protocol_version") != 7
        or pro7_protocol.get("frozen") is not True
        or pro7_protocol.get("outcome_neutral_contract_validated") is not True
        or pro7_protocol.get("accepted_state_source_failure_terminal") is not True
        or pro7_protocol.get("all_source_only_unaccepted_proposal_failures_retryable")
        is not True
        or pro7_protocol.get("candidate_endpoint_in_unaccepted_proposal") is not True
        or pro7_protocol.get("non_source_failure_vetoes_retry") is not True
        or pro7_protocol.get("complete_terminal_observability_required") is not True
        or pro7_protocol.get("raw_source_residual_tolerance_unchanged") is not True
        or pro7_protocol.get("all_non_source_stops_unchanged") is not True
        or pro7_protocol.get("retry_factor") != "1/2"
        or pro7_protocol.get("maximum_retry_count") != 32
        or pro7_protocol.get("minimum_step_size") != "1/1073741824"
        or pro7_protocol.get("resolved_holdout_manifest")
        != "configs/fgc/fgc-1-pro7-hld1.toml"
        or pro7_protocol.get("calibration_output_root")
        != "runs/fgc-2-sf1/proto7/calibration"
        or pro7_protocol.get("holdout_output_root")
        != "runs/fgc-2-sf1/proto7/holdout"
        or any(value is not False for value in pro7_protocol.get("claims", {}).values())
        or pro7_lineage.get("checkpoint_commit")
        != "27c63e7aca1439e75e18b576c6da414e458280ff"
        or pro7_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or pro7_predecessor.get("artifact_id") != "FGC-2-SF1-PROTO6"
        or pro7_predecessor.get("protocol_version") != 6
        or pro7_predecessor.get("claims_all_false") is not True
        or pro7_diagnosis.get("artifact_id") != "FGC-1-CAL3-PREF5"
        or pro7_diagnosis.get("PROTO6_transaction_ownership_obstructed") is not True
        or pro7_diagnosis.get("PROTO7_general_transaction_revision_required") is not True
        or pro7_diagnosis.get("accepted_state_source_gate_passed") is not True
        or pro7_diagnosis.get("terminal_failures_source_only") is not True
        or pro7_diagnosis.get("failed_terminal_stages")
        != ["rk4_k4", "candidate_endpoint"]
        or pro7_diagnosis.get("one_further_halving_complete_transaction_passed")
        is not True
        or pro7_diagnosis.get("first_nonzero_common_event_completed") is not False
        or pro7_diagnosis.get("FGCQR_outcome_read") is not False
        or pro7_diagnosis.get("mechanism_question_answered") is not False
        or pro7_repair != expected_pro7_repair
        or pro7_boundary
        != {
            "FGCQR_outcome_read": False,
            "PROTO6_GR0_calibration_trajectory_was_seen": True,
            "PROTO7_trajectory_read": False,
            "SGBL_outcome_read": False,
            "first_nonzero_common_event_completed": False,
            "mechanism_question_answered": False,
            "physical_experiment_changed": False,
        }
        or pro7_payload.get("decision")
        != "PROTO7_is_frozen_but_HLT5_fresh_calibration_and_every_candidate_outcome_remain_closed"
        or set(pro7.get("implementation_sha256", {}))
        != expected_pro7_implementation
        or any(value is not False for value in pro7.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-PRO7-FRZ1 lineage, transaction ownership, observability, or nonclaim contract has drifted"
        )

    hlt5_path = REPOSITORY / "results/fgc-1-hlt5-mon5.json"
    try:
        hlt5 = json.loads(
            hlt5_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-HLT5-MON5 certificate: {exc}")
        return failures
    hlt5_payload = hlt5.get("artifact_payload", {})
    hlt5_lineage = hlt5_payload.get("immutable_lineage", {})
    hlt5_runtime = hlt5_payload.get("runtime_controls", {})
    hlt5_runner = hlt5_payload.get("runner_controls", {})
    hlt5_equivalence = hlt5_payload.get("accepted_source_operator_equivalence", {})
    hlt5_plan = hlt5_payload.get("frozen_run_plan", {})
    hlt5_inputs = hlt5_payload.get("frozen_run_inputs", [])
    hlt5_prechecks = hlt5_payload.get("accepted_state_source_prechecks", [])
    hlt5_events = hlt5_payload.get("initial_common_events", [])
    hlt5_namespace = hlt5_payload.get("namespace_precondition", {})
    retry_positions = hlt5_runtime.get(
        "source_only_retry_controls_at_every_proposal_position", []
    )
    veto_positions = hlt5_runtime.get(
        "non_source_veto_controls_at_every_proposal_position", []
    )
    expected_hlt5_implementation = {
        "scripts/reproduce_fgc_hlt5_mon5.py",
        "scripts/run_fgc_gr0_calibration.py",
        "scripts/run_fgc_gr0_calibration_v6.py",
        "scripts/run_fgc_gr0_calibration_v7.py",
        "src/recursive_horizons/fgc/evolution/cal2_source_diagnosis.py",
        "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto5_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto6_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto7_runtime.py",
    }
    if (
        hlt5.get("classification")
        != "pre_trajectory_PROTO7_GR0_general_transaction_runtime_and_fresh_calibration_authorization"
        or hlt5.get("generated_by") != "scripts/reproduce_fgc_hlt5_mon5.py"
        or hlt5.get("scope_bindings")
        != {
            "FGCQR_trajectory_read": False,
            "PROTO7_trajectory_read": False,
            "SGBL_trajectory_read": False,
            "calibration_branch": "GR-0",
            "candidate_action_health_values_invented_for_GR0": False,
            "collapse_or_trapped_outcome_classified": False,
            "role": "pre_trajectory_general_unaccepted_proposal_runtime_compositor_and_fresh_input_freeze",
            "runtime_owner": "FGC-1-HLT5-MON5",
            "target_protocol": "FGC-2-SF1-PROTO7",
        }
        or hlt5_lineage.get("checkpoint_commit")
        != "ec8ac628528b822c600de35eb2ca6b2b49506051"
        or hlt5_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or hlt5_lineage.get("protocol", {}).get("artifact_id")
        != "FGC-2-SF1-PROTO7"
        or hlt5_lineage.get("protocol_freeze_artifact_id")
        != "FGC-1-PRO7-FRZ1"
        or hlt5_lineage.get("predecessor_runtime_artifact_id")
        != "FGC-1-HLT4-MON4"
        or hlt5_lineage.get("static_input_artifact_id") != "FGC-1-ID2-ALL1"
        or hlt5_runtime.get("all_controls_passed") is not True
        or [item.get("observed_stage") for item in retry_positions]
        != ["rk4_k1", "rk4_k2", "rk4_k3", "rk4_k4", "candidate_endpoint"]
        or any(
            item.get("retryable") is not True
            or item.get("state_unchanged") is not True
            or item.get("causal_debit_unchanged") is not True
            for item in retry_positions
        )
        or len(veto_positions) != 5
        or any(
            item.get("reason") != "nonpositive_lapse"
            or item.get("non_source_failure_vetoed_retry") is not True
            for item in veto_positions
        )
        or hlt5_runtime.get("mixed_chronological_priority_control", {}).get(
            "selected_reason"
        )
        != "newton_residual_limit"
        or hlt5_runner.get("all_controls_passed") is not True
        or hlt5_runner.get("candidate_endpoint_retry_then_accept", {}).get(
            "serialized_retry_count"
        )
        != 1
        or hlt5_runner.get("retry_exhaustion", {}).get(
            "serialized_retry_count"
        )
        != 2
        or hlt5_equivalence.get("all_controls_passed") is not True
        or hlt5_plan.get("artifact_id") != "FGC-1-CAL4-RUN1-PLAN"
        or hlt5_plan.get("campaign_id")
        != "f1e5d99d061f73417e570e7bb927a0f82458e1002138c68b02360b725d0b3f03"
        or hlt5_plan.get("ordered_amplitudes") != ["5/2", "3"]
        or hlt5_plan.get("expanded_input_count") != 12
        or len(hlt5_inputs) != 12
        or len(hlt5_prechecks) != 12
        or any(
            item.get("accepted_state_source_gate_passed") is not True
            or item.get("trajectory_advanced") is not False
            for item in hlt5_prechecks
        )
        or len(hlt5_events) != 4
        or any(
            item.get("admission_passed") is not True
            or item.get("trajectory_advanced") is not False
            for item in hlt5_events
        )
        or hlt5_namespace.get("both_fresh_namespaces_empty") is not True
        or hlt5_namespace.get("authorization_created_no_namespace") is not True
        or [item.get("path") for item in hlt5_namespace.get("records", [])]
        != [
            "runs/fgc-2-sf1/proto7/calibration",
            "runs/fgc-2-sf1/proto7/holdout",
        ]
        or hlt5_payload.get("decision")
        != "HLT5_implements_only_PROTO7_transaction_ownership_and_authorizes_one_fresh_GR0_recalibration"
        or hlt5_payload.get("epistemic_boundary", {}).get("PROTO7_trajectory_read")
        is not False
        or hlt5_payload.get("epistemic_boundary", {}).get(
            "fresh_GR0_recalibration_is_now_authorized_but_not_completed"
        )
        is not True
        or set(hlt5.get("implementation_sha256", {}))
        != expected_hlt5_implementation
        or any(value is not False for value in hlt5.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-HLT5-MON5 lineage, runtime controls, input freeze, authorization, or nonclaim contract has drifted"
        )

    cal4_path = REPOSITORY / "results/fgc-1-cal4-pref6.json"
    try:
        cal4 = json.loads(
            cal4_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-CAL4-PREF6 certificate: {exc}")
        return failures
    cal4_payload = cal4.get("artifact_payload", {})
    cal4_campaign = cal4_payload.get("immutable_campaign", {})
    cal4_diagnosis = cal4_payload.get("common_event_ownership_diagnosis", {})
    cal4_amp5 = cal4_diagnosis.get("amplitude_5over2", {})
    cal4_amp3 = cal4_diagnosis.get("amplitude_3", {})
    cal4_decision = cal4_payload.get("decision", {})
    cal4_boundary = cal4_payload.get("epistemic_boundary", {})
    expected_cal4_implementation = {
        "scripts/reproduce_fgc_cal4_pref6.py",
        "scripts/run_fgc_gr0_calibration_v7.py",
        "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py",
    }
    if (
        cal4_campaign.get("authorization_commit")
        != "b8f1921985191b9ad488ecfd2f183c19936996b1"
        or cal4_campaign.get("campaign_id")
        != "f1e5d99d061f73417e570e7bb927a0f82458e1002138c68b02360b725d0b3f03"
        or cal4_campaign.get("terminal_classification")
        != "calibration_failed_no_eligible_GR0_case"
        or cal4_campaign.get("first_nonzero_common_event_completed_by_both_amplitudes")
        is not True
        or cal4_campaign.get("event_log_line_count") != 8
        or cal4_amp5.get("both_methods_prospectively_admitted") is not True
        or cal4_amp3.get("both_methods_prospectively_admitted") is not False
        or set(cal4_amp5.get("methods", {})) != {"RK4", "SSPRK3"}
        or set(cal4_amp3.get("methods", {})) != {"RK4", "SSPRK3"}
        or any(
            method.get("prospective_common_event_admission_passed") is not True
            or method.get("old_full_domain_reduction_failures_boundary_localized")
            is not True
            or method.get(
                "owned_reduction_residuals_within_accumulated_roundoff_enclosure"
            )
            is not True
            for method in cal4_amp5.get("methods", {}).values()
        )
        or cal4_decision.get("PROTO7_common_event_ownership_contract_obstructed")
        is not True
        or cal4_decision.get("new_protocol_required") != "FGC-2-SF1-PROTO8"
        or cal4_boundary
        != {
            "GR0_amplitude_selected": False,
            "GR0_calibration_completed": False,
            "PROTO7_transaction_repair_succeeded_through_t1": True,
            "SGBL_or_FGCQR_trajectory_read": False,
            "candidate_or_general_gradient_route_rejected": False,
            "dynamic_trapped_sphere_classified": False,
            "mechanism_question_answered": False,
        }
        or set(cal4.get("implementation_sha256", {}))
        != expected_cal4_implementation
    ):
        failures.append(
            "FGC-1-CAL4-PREF6 campaign, ownership diagnosis, successor decision, or epistemic boundary has drifted"
        )

    pro8_path = REPOSITORY / "results/fgc-1-pro8-frz1.json"
    try:
        pro8 = json.loads(
            pro8_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-PRO8-FRZ1 certificate: {exc}")
        return failures
    pro8_payload = pro8.get("artifact_payload", {})
    pro8_protocol = pro8_payload.get("protocol_certificate", {})
    pro8_lineage = pro8_payload.get("immutable_lineage", {})
    pro8_repair = pro8_payload.get("frozen_repair", {})
    pro8_boundary = pro8_payload.get("epistemic_boundary", {})
    expected_pro8_implementation = {
        "scripts/reproduce_fgc_pro8_frz1.py",
        "src/recursive_horizons/fgc/evolution/protocol_v7.py",
        "src/recursive_horizons/fgc/evolution/protocol_v8.py",
    }
    if (
        pro8.get("classification")
        != "immutable_outcome_neutral_PROTO8_common_event_diagnostic_ownership_protocol_freeze_not_run_authorization"
        or pro8.get("generated_by") != "scripts/reproduce_fgc_pro8_frz1.py"
        or pro8_lineage.get("checkpoint_commit")
        != "219ba9493899029bc51bd313e3c5f689f3809c2c"
        or pro8_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or pro8_lineage.get("predecessor_protocol", {}).get("artifact_id")
        != "FGC-2-SF1-PROTO7"
        or pro8_lineage.get("diagnosis", {}).get("artifact_id")
        != "FGC-1-CAL4-PREF6"
        or pro8_lineage.get("diagnosis", {}).get(
            "first_nonzero_common_event_completed_by_both_amplitudes"
        )
        is not True
        or pro8_lineage.get("diagnosis", {}).get(
            "amplitude_5over2_prospectively_admitted"
        )
        is not True
        or pro8_lineage.get("diagnosis", {}).get(
            "amplitude_3_prospectively_admitted"
        )
        is not False
        or pro8_protocol.get("artifact_id") != "FGC-2-SF1-PROTO8"
        or pro8_protocol.get("protocol_version") != 8
        or pro8_protocol.get("outcome_neutral_contract_validated") is not True
        or pro8_protocol.get("RK4_excluded_outer_rows") != 7
        or pro8_protocol.get("SSPRK3_excluded_outer_rows") != 5
        or pro8_protocol.get("all_physical_and_numerical_thresholds_unchanged")
        is not True
        or pro8_repair.get("accumulated_roundoff_changes_order_classification_only")
        is not True
        or pro8_repair.get("medium_and_finest_absolute_spectral_budgets_required")
        is not True
        or pro8_repair.get("all_adjacent_nested_tail_ratios_required") is not True
        or pro8_boundary
        != {
            "PROTO7_GR0_calibration_trajectory_was_seen": True,
            "first_nonzero_common_event_completed": True,
            "PROTO8_trajectory_read": False,
            "FGCQR_outcome_read": False,
            "SGBL_outcome_read": False,
            "mechanism_question_answered": False,
            "physical_experiment_changed": False,
            "numerical_threshold_changed": False,
        }
        or pro8_payload.get("decision")
        != "PROTO8_is_frozen_but_HLT6_fresh_calibration_and_every_candidate_outcome_remain_closed"
        or set(pro8.get("implementation_sha256", {}))
        != expected_pro8_implementation
        or any(value is not False for value in pro8.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-PRO8-FRZ1 lineage, ownership freeze, unchanged-threshold boundary, or nonclaim contract has drifted"
        )

    hlt6_path = REPOSITORY / "results/fgc-1-hlt6-mon6.json"
    try:
        hlt6 = json.loads(
            hlt6_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-HLT6-MON6 certificate: {exc}")
        return failures
    hlt6_payload = hlt6.get("artifact_payload", {})
    hlt6_lineage = hlt6_payload.get("immutable_lineage", {})
    hlt6_transaction = hlt6_payload.get("inherited_PROTO7_transaction", {})
    hlt6_controls = hlt6_payload.get("common_event_adapter_controls", {})
    hlt6_spectrum = hlt6_controls.get("spatial_spectrum", {})
    hlt6_plan = hlt6_payload.get("frozen_run_plan", {})
    hlt6_inputs = hlt6_payload.get("frozen_run_inputs", [])
    hlt6_prechecks = hlt6_payload.get("accepted_state_source_prechecks", [])
    hlt6_events = hlt6_payload.get("initial_common_events", [])
    hlt6_namespace = hlt6_payload.get("namespace_precondition", {})
    hlt6_boundary = hlt6_payload.get("epistemic_boundary", {})
    expected_hlt6_implementation = {
        "scripts/reproduce_fgc_hlt6_mon6.py",
        "scripts/run_fgc_gr0_calibration.py",
        "scripts/run_fgc_gr0_calibration_v6.py",
        "scripts/run_fgc_gr0_calibration_v7.py",
        "scripts/run_fgc_gr0_calibration_v8.py",
        "src/recursive_horizons/fgc/evolution/cal2_source_diagnosis.py",
        "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py",
        "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto5_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto6_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto7_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto8_runtime.py",
    }
    expected_hlt6_scope = {
        "FGCQR_trajectory_read": False,
        "PROTO7_calibration_history_disclosed": True,
        "PROTO8_trajectory_read": False,
        "SGBL_trajectory_read": False,
        "calibration_branch": "GR-0",
        "candidate_action_health_values_invented_for_GR0": False,
        "collapse_or_trapped_outcome_classified": False,
        "role": "pre_trajectory_common_event_diagnostic_ownership_runtime_compositor_and_fresh_input_freeze",
        "runtime_owner": "FGC-1-HLT6-MON6",
        "target_protocol": "FGC-2-SF1-PROTO8",
    }
    if (
        hlt6.get("classification")
        != "pre_trajectory_PROTO8_GR0_common_event_ownership_runtime_and_fresh_calibration_authorization"
        or hlt6.get("generated_by") != "scripts/reproduce_fgc_hlt6_mon6.py"
        or hlt6.get("scope_bindings") != expected_hlt6_scope
        or hlt6_lineage.get("checkpoint_commit")
        != "b4c3cdcd216e305ef55e2f2700fbee4111afc95c"
        or hlt6_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or hlt6_lineage.get("protocol", {}).get("artifact_id")
        != "FGC-2-SF1-PROTO8"
        or hlt6_lineage.get("protocol", {}).get("protocol_version") != 8
        or hlt6_lineage.get("freeze", {}).get("artifact_id")
        != "FGC-1-PRO8-FRZ1"
        or hlt6_lineage.get("freeze", {}).get("claims_all_false") is not True
        or hlt6_lineage.get("predecessor_runtime", {}).get("artifact_id")
        != "FGC-1-HLT5-MON5"
        or hlt6_lineage.get("predecessor_runtime", {}).get(
            "implementation_files_hash_matched"
        )
        != 9
        or hlt6_lineage.get("static_input", {}).get("artifact_id")
        != "FGC-1-ID2-ALL1"
        or hlt6_transaction.get("inherited_implementation_files_hash_matched")
        != 9
        or hlt6_transaction.get(
            "source_solver_transaction_retry_and_non_source_stops_unchanged"
        )
        is not True
        or hlt6_controls.get("all_controls_passed") is not True
        or hlt6_controls.get("outer_projector_defect", {}).get(
            "public_but_not_owned_control_passed"
        )
        is not True
        or hlt6_controls.get("interior_defect", {}).get(
            "interior_veto_control_passed"
        )
        is not True
        or any(
            hlt6_spectrum.get(key) is not True
            for key in (
                "coarse_witness_finest_pair_control_passed",
                "medium_absolute_miss_vetoed",
                "noncontracting_tail_vetoed",
            )
        )
        or hlt6_plan.get("artifact_id") != "FGC-1-CAL5-RUN1-PLAN"
        or hlt6_plan.get("campaign_id")
        != "721dd0f590c2546d2fd9c5b8b2e7d0185938da9bf1ea38e45a8ddde2de4e0ffb"
        or hlt6_plan.get("ordered_amplitudes") != ["5/2", "3"]
        or hlt6_plan.get("expanded_input_count") != 12
        or len(hlt6_inputs) != 12
        or any(
            item.get("trajectory_advanced") is not False
            or item.get("u_bitwise_preserved") is not True
            or item.get("p_bitwise_preserved") is not True
            or item.get("q_initialized_by_native_SBP_operator") is not True
            or item.get("source_stop_ownership")
            != "accepted_terminal_all_unaccepted_proposal_source_only_retry"
            for item in hlt6_inputs
        )
        or len(hlt6_prechecks) != 12
        or any(
            item.get("accepted_state_source_gate_passed") is not True
            or item.get("trajectory_advanced") is not False
            for item in hlt6_prechecks
        )
        or len(hlt6_events) != 4
        or any(
            item.get("admission_passed") is not True
            or item.get("trajectory_advanced") is not False
            or item.get("accepted_stage_counts") != [0, 0, 0]
            or item.get("finest_pair_spatial_individual_budgets_passed")
            is not True
            or item.get("all_adjacent_spatial_tail_ratios_passed") is not True
            or item.get("excluded_outer_rows")
            != ([7, 7, 7] if item.get("method") == "RK4" else [5, 5, 5])
            for item in hlt6_events
        )
        or hlt6_namespace.get("both_fresh_namespaces_empty") is not True
        or hlt6_namespace.get("authorization_created_no_namespace") is not True
        or [item.get("path") for item in hlt6_namespace.get("records", [])]
        != [
            "runs/fgc-2-sf1/proto8/calibration",
            "runs/fgc-2-sf1/proto8/holdout",
        ]
        or hlt6_payload.get("decision")
        != "HLT6_implements_only_PROTO8_common_event_ownership_and_authorizes_one_fresh_GR0_recalibration"
        or hlt6_boundary.get("PROTO8_trajectory_read") is not False
        or hlt6_boundary.get("PROTO8_trajectory_advanced") is not False
        or hlt6_boundary.get(
            "fresh_GR0_recalibration_is_now_authorized_but_not_completed"
        )
        is not True
        or hlt6_boundary.get("mechanism_question_answered") is not False
        or set(hlt6.get("implementation_sha256", {}))
        != expected_hlt6_implementation
        or any(value is not False for value in hlt6.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-HLT6-MON6 lineage, adapter controls, input freeze, authorization, or nonclaim contract has drifted"
        )

    cal5_path = REPOSITORY / "results/fgc-1-cal5-pref7.json"
    try:
        cal5 = json.loads(
            cal5_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-CAL5-PREF7 certificate: {exc}")
        return failures
    cal5_payload = cal5.get("artifact_payload", {})
    cal5_campaign = cal5_payload.get("immutable_campaign", {})
    cal5_diagnosis = cal5_payload.get("resolution_ladder_diagnosis", {})
    cal5_5 = cal5_diagnosis.get("amplitude_5over2_stop", {})
    cal5_3 = cal5_diagnosis.get("amplitude_3_stop", {})
    cal5_decision = cal5_payload.get("decision", {})
    cal5_revision = cal5_decision.get("prospective_revision", {})
    cal5_boundary = cal5_payload.get("epistemic_boundary", {})
    expected_cal5_implementation = {
        "scripts/reproduce_fgc_cal5_pref7.py",
        "src/recursive_horizons/fgc/evolution/cal5_resolution_diagnosis.py",
    }
    cal5_methods_valid = True
    for amplitude, record, event_index, event_time in (
        ("5/2", cal5_5, 2, 0.125),
        ("3", cal5_3, 1, 0.0625),
    ):
        if (
            record.get("amplitude") != amplitude
            or record.get("event_index") != event_index
            or record.get("coordinate_time") != event_time
            or record.get("one_level_resolution_shift_is_a_discriminating_test")
            is not True
            or record.get("one_level_resolution_shift_is_already_admitted")
            is not False
        ):
            cal5_methods_valid = False
            continue
        methods = record.get("methods", {})
        if set(methods) != {"RK4", "SSPRK3"}:
            cal5_methods_valid = False
            continue
        for method, item in methods.items():
            if (
                item.get("method") != method
                or item.get("PROTO8_admission_passed") is not False
                or item.get("constraint_admission_passed") is not False
                or item.get("prospective_2049_coarsest_guard_passed") is not True
                or item.get("conditional_8193_finest_guard_passed") is not True
                or item.get("conditional_projection_is_admission_evidence")
                is not False
                or item.get("medium_and_fine_absolute_spectral_budgets_passed")
                is not True
                or item.get("new_4097_to_8193_tail_is_unmeasured") is not True
            ):
                cal5_methods_valid = False
        if amplitude == "5/2" and any(
            methods[method].get("current_nested_tail_failures") != []
            for method in methods
        ):
            cal5_methods_valid = False
        if amplitude == "3":
            expected_ratios = {
                "RK4": 0.3267172814251121,
                "SSPRK3": 0.26675704013240803,
            }
            for method, ratio in expected_ratios.items():
                if methods[method].get("current_nested_tail_failures") != [
                    {
                        "family": "derivative",
                        "field": "phi_over_Lambda",
                        "maximum": 0.25,
                        "pair": "2049_to_4097",
                        "ratio": ratio,
                    }
                ]:
                    cal5_methods_valid = False
    if (
        cal5.get("classification")
        != "post_PROTO8_GR0_resolution_ladder_adequacy_diagnosis"
        or cal5.get("generated_by") != "scripts/reproduce_fgc_cal5_pref7.py"
        or cal5_campaign.get("authorization_commit")
        != "d667c7465786883b1e5d6afad168b64e7b8921d0"
        or cal5_campaign.get("campaign_id")
        != "721dd0f590c2546d2fd9c5b8b2e7d0185938da9bf1ea38e45a8ddde2de4e0ffb"
        or cal5_campaign.get("event_log_line_count") != 9
        or cal5_campaign.get("terminal_classification")
        != "calibration_failed_no_eligible_GR0_case"
        or cal5_campaign.get("selected_amplitude") is not None
        or cal5_campaign.get("terminal_checkpoint_amplitude") != "3"
        or cal5_campaign.get("terminal_checkpoint_event_index") != 1
        or cal5_diagnosis.get("public_event_count") != 9
        or cal5_diagnosis.get("common_event_count") != 5
        or cal5_diagnosis.get("source_only_retry_count") != 4
        or cal5_diagnosis.get("current_point_counts") != [1025, 2049, 4097]
        or cal5_diagnosis.get("smallest_prospective_point_counts")
        != [2049, 4097, 8193]
        or cal5_diagnosis.get("threshold_changes_required") is not False
        or cal5_diagnosis.get("new_8193_data_required_before_admission")
        is not True
        or cal5_diagnosis.get("candidate_or_physical_trajectory_read") is not False
        or cal5_methods_valid is not True
        or cal5_decision.get(
            "PROTO8_ownership_repair_exercised_by_evolved_data"
        )
        is not True
        or cal5_decision.get("PROTO8_GR0_case_eligible") is not False
        or cal5_decision.get(
            "current_failure_is_a_candidate_or_gradient_obstruction"
        )
        is not False
        or cal5_decision.get("new_protocol_required") != "FGC-2-SF1-PROTO9"
        or cal5_revision.get("only_resolution_ladder_may_change") is not True
        or cal5_revision.get("point_counts") != [2049, 4097, 8193]
        or any(
            cal5_revision.get(key) is not True
            for key in (
                "physical_inputs_unchanged",
                "amplitude_order_unchanged",
                "evolution_equations_unchanged",
                "source_solver_and_raw_residual_tolerance_unchanged",
                "methods_CFL_and_retry_rules_unchanged",
                "constraint_ownership_magnitude_and_order_rules_unchanged",
                "spectral_ownership_budgets_and_tail_rules_unchanged",
                "trapped_sign_boundary_and_physical_rules_unchanged",
                "fresh_output_namespace_required",
            )
        )
        or any(value is not False for value in cal5_boundary.values())
        or set(cal5.get("implementation_sha256", {}))
        != expected_cal5_implementation
        or any(value is not False for value in cal5.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-CAL5-PREF7 campaign lineage, fixed-threshold resolution diagnosis, or nonclaim contract has drifted"
        )

    pro9_path = REPOSITORY / "results/fgc-1-pro9-frz1.json"
    try:
        pro9 = json.loads(
            pro9_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-PRO9-FRZ1 certificate: {exc}")
        return failures
    pro9_payload = pro9.get("artifact_payload", {})
    pro9_lineage = pro9_payload.get("immutable_lineage", {})
    pro9_predecessor = pro9_lineage.get("predecessor_protocol", {})
    pro9_diagnosis = pro9_lineage.get("diagnosis", {})
    pro9_protocol = pro9_payload.get("protocol_validation", {})
    pro9_revision = pro9_payload.get("resolution_revision", {})
    pro9_namespace = pro9_payload.get("namespace_precondition", {})
    pro9_boundary = pro9_payload.get("epistemic_boundary", {})
    expected_pro9_implementation = {
        "scripts/reproduce_fgc_pro9_frz1.py",
        "src/recursive_horizons/fgc/evolution/protocol_v9.py",
    }
    expected_pro9_claims = {
        "FGCQR_holdout_execution_authorized",
        "FGCQR_mechanism_rejected",
        "PROTO9_fresh_GR0_dynamic_calibration_authorized",
        "PROTO9_resolved_holdout_manifest_authorized",
        "PROTO9_successor_runtime_implemented",
        "child_domain_or_topology_derived",
        "classical_spherical_diagnostic_authorized",
        "dark_sector_mechanism_derived",
        "general_gradient_route_rejected",
        "physical_transition_claim_authorized",
        "retained_EFT_evolution_authorized",
        "singularity_resolution_derived",
        "varying_locally_measured_c_derived",
    }
    if (
        pro9.get("classification")
        != "outcome_neutral_PROTO9_resolution_ladder_protocol_freeze"
        or pro9.get("generated_by") != "scripts/reproduce_fgc_pro9_frz1.py"
        or pro9_lineage.get("checkpoint_commit")
        != "a2c929e33ec8f2ee5464b086a57db8ed7f9277f0"
        or pro9_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or pro9_lineage.get("predecessor_protocol_sha256")
        != "5de88456ed041f8a4293150f49b4caea057560dd796c8ea2b9a8110ecf7770ab"
        or pro9_predecessor.get("artifact_id") != "FGC-2-SF1-PROTO8"
        or pro9_predecessor.get("protocol_version") != 8
        or pro9_predecessor.get("frozen") is not True
        or pro9_diagnosis.get("artifact_id") != "FGC-1-CAL5-PREF7"
        or pro9_diagnosis.get("classification")
        != "post_PROTO8_GR0_resolution_ladder_adequacy_diagnosis"
        or pro9_diagnosis.get("campaign_id")
        != "721dd0f590c2546d2fd9c5b8b2e7d0185938da9bf1ea38e45a8ddde2de4e0ffb"
        or pro9_diagnosis.get("terminal_classification")
        != "calibration_failed_no_eligible_GR0_case"
        or pro9_diagnosis.get("selected_amplitude") is not None
        or pro9_diagnosis.get("current_point_counts") != [1025, 2049, 4097]
        or pro9_diagnosis.get("prospective_point_counts")
        != [2049, 4097, 8193]
        or pro9_diagnosis.get("new_8193_data_required") is not True
        or pro9_diagnosis.get("claims_all_candidate_and_physical_false")
        is not True
        or pro9_diagnosis.get("config_sha256")
        != "02ab7bbe58d42661cbf40017acf4e7d935a82d47f3a4f5cf2ea52a03a56a1732"
        or pro9_diagnosis.get("result_sha256")
        != "66404f44dae76021eb98268e6c48e72f17058423837177c410445d20c178b6fa"
        or pro9_protocol.get("artifact_id") != "FGC-2-SF1-PROTO9"
        or pro9_protocol.get("protocol_version") != 9
        or pro9_protocol.get("frozen") is not True
        or pro9_protocol.get("outcome_neutral_contract_validated") is not True
        or pro9_protocol.get("premise_revision_only") is not True
        or pro9_protocol.get("predecessor_protocol_artifact_id")
        != "FGC-2-SF1-PROTO8"
        or pro9_protocol.get("diagnosis_artifact_id") != "FGC-1-CAL5-PREF7"
        or pro9_protocol.get("diagnosis_checkpoint_commit")
        != "a2c929e33ec8f2ee5464b086a57db8ed7f9277f0"
        or pro9_protocol.get("eligible_calibration_amplitudes") != ["5/2", "3"]
        or pro9_protocol.get("old_point_counts") != [1025, 2049, 4097]
        or pro9_protocol.get("point_counts") != [2049, 4097, 8193]
        or pro9_protocol.get("only_resolution_ladder_changes") is not True
        or pro9_protocol.get("constraint_rules_unchanged") is not True
        or pro9_protocol.get("spectral_rules_unchanged") is not True
        or pro9_protocol.get("physical_and_numerical_thresholds_unchanged")
        is not True
        or pro9_protocol.get("new_8193_data_required") is not True
        or pro9_protocol.get("conditional_projection_is_admission") is not False
        or pro9_protocol.get(
            "PROTO8_calibration_trajectory_inspected_before_revision"
        )
        is not True
        or pro9_protocol.get("SGBL_outcomes_inspected_before_revision")
        is not False
        or pro9_protocol.get("FGCQR_outcomes_inspected_before_revision")
        is not False
        or pro9_protocol.get("mechanism_question_answered_before_revision")
        is not False
        or pro9_protocol.get("calibration_output_root")
        != "runs/fgc-2-sf1/proto9/calibration"
        or pro9_protocol.get("holdout_output_root")
        != "runs/fgc-2-sf1/proto9/holdout"
        or pro9_protocol.get("resolved_holdout_manifest")
        != "configs/fgc/fgc-1-pro9-hld1.toml"
        or pro9_protocol.get("resolved_holdout_manifest_present") is not False
        or set(pro9_protocol.get("claims", {})) != expected_pro9_claims
        or any(value is not False for value in pro9_protocol.get("claims", {}).values())
        or pro9_revision.get("old_point_counts") != [1025, 2049, 4097]
        or pro9_revision.get("point_counts") != [2049, 4097, 8193]
        or pro9_revision.get("only_resolution_ladder_changes") is not True
        or pro9_revision.get("constraint_rules_unchanged") is not True
        or pro9_revision.get("spectral_rules_unchanged") is not True
        or pro9_revision.get(
            "all_other_physical_and_numerical_rules_unchanged"
        )
        is not True
        or pro9_revision.get("conditional_projection_is_admission") is not False
        or pro9_revision.get("new_8193_data_exist") is not False
        or pro9_namespace.get("both_fresh_namespaces_absent") is not True
        or pro9_namespace.get("freeze_created_no_namespace") is not True
        or pro9_namespace.get("records")
        != [
            {
                "exists": False,
                "path": "runs/fgc-2-sf1/proto9/calibration",
            },
            {
                "exists": False,
                "path": "runs/fgc-2-sf1/proto9/holdout",
            },
        ]
        or pro9_payload.get("decision")
        != "PROTO9_is_frozen_but_HLT7_fresh_calibration_and_every_candidate_outcome_remain_closed"
        or set(pro9_boundary)
        != {
            "GR0_calibration_completed",
            "PROTO9_trajectory_advanced",
            "PROTO9_trajectory_read",
            "SGBL_or_FGCQR_trajectory_read",
            "mechanism_question_answered",
            "new_8193_state_constructed",
        }
        or any(value is not False for value in pro9_boundary.values())
        or set(pro9.get("implementation_sha256", {}))
        != expected_pro9_implementation
        or set(pro9.get("gate_status", {})) != expected_pro9_claims
        or any(value is not False for value in pro9.get("gate_status", {}).values())
        or any(value is not False for value in pro9.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-PRO9-FRZ1 immutable lineage, resolution-only revision, namespace, or nonclaim contract has drifted"
        )

    hlt7_path = REPOSITORY / "results/fgc-1-hlt7-mon7.json"
    try:
        hlt7 = json.loads(
            hlt7_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-HLT7-MON7 certificate: {exc}")
        return failures
    hlt7_payload = hlt7.get("artifact_payload", {})
    hlt7_lineage = hlt7_payload.get("immutable_lineage", {})
    hlt7_protocol = hlt7_lineage.get("protocol", {})
    hlt7_freeze = hlt7_lineage.get("freeze", {})
    hlt7_predecessor = hlt7_lineage.get("predecessor_runtime", {})
    hlt7_controls = hlt7_payload.get("resolution_adapter_controls", {})
    hlt7_inherited = hlt7_payload.get("inherited_runtime", {})
    hlt7_plan = hlt7_payload.get("frozen_run_plan", {})
    hlt7_inputs = hlt7_payload.get("frozen_run_inputs", [])
    hlt7_prechecks = hlt7_payload.get("accepted_state_source_prechecks", [])
    hlt7_events = hlt7_payload.get("initial_common_events", [])
    hlt7_namespace = hlt7_payload.get("namespace_precondition", {})
    hlt7_boundary = hlt7_payload.get("epistemic_boundary", {})
    expected_hlt7_implementation = {
        "scripts/reproduce_fgc_hlt7_mon7.py",
        "scripts/run_fgc_gr0_calibration.py",
        "scripts/run_fgc_gr0_calibration_v6.py",
        "scripts/run_fgc_gr0_calibration_v7.py",
        "scripts/run_fgc_gr0_calibration_v8.py",
        "scripts/run_fgc_gr0_calibration_v9.py",
        "src/recursive_horizons/fgc/evolution/cal2_source_diagnosis.py",
        "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py",
        "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto5_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto6_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto7_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto8_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto9_runtime.py",
        "src/recursive_horizons/fgc/evolution/protocol_v9.py",
    }
    expected_hlt7_order = [
        (amplitude, method, count)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
        for count in (2049, 4097, 8193)
    ]
    hlt7_inputs_valid = (
        len(hlt7_inputs) == 12
        and [
            (item.get("amplitude"), item.get("method"), item.get("point_count"))
            for item in hlt7_inputs
        ]
        == expected_hlt7_order
        and all(
            item.get("branch") == "GR-0"
            and item.get("u_bitwise_preserved") is True
            and item.get("p_bitwise_preserved") is True
            and item.get("q_initialized_by_native_SBP_operator") is True
            and item.get("conditional_projection_used") is False
            and item.get("trajectory_advanced") is False
            and item.get("source_stop_ownership")
            == "accepted_terminal_all_unaccepted_proposal_source_only_retry"
            and item.get("new_8193_state_constructed")
            is (item.get("point_count") == 8193)
            and item.get("overlap_projected_state_matches_HLT6")
            is (item.get("point_count") in {2049, 4097})
            for item in hlt7_inputs
        )
    )
    hlt7_events_valid = (
        len(hlt7_events) == 4
        and [
            (item.get("amplitude"), item.get("method")) for item in hlt7_events
        ]
        == [
            ("5/2", "RK4"),
            ("5/2", "SSPRK3"),
            ("3", "RK4"),
            ("3", "SSPRK3"),
        ]
        and all(
            item.get("point_counts") == [2049, 4097, 8193]
            and item.get("accepted_stage_counts") == [0, 0, 0]
            and item.get("excluded_outer_rows")
            == ([7, 7, 7] if item.get("method") == "RK4" else [5, 5, 5])
            and item.get("constraint_admission_passed") is True
            and item.get("constraint_coarsest_guard_passed") is True
            and item.get("constraint_finest_guard_passed") is True
            and item.get("constraint_monotone_refinement_passed") is True
            and item.get("constraint_finest_pair_order_passed") is True
            and item.get("constraint_minimum_finite_finest_pair_order", 0.0) > 1.5
            and item.get("finest_pair_spatial_individual_budgets_passed") is True
            and item.get("all_adjacent_spatial_tail_ratios_passed") is False
            and item.get("new_8193_state_present") is True
            and item.get("conditional_projection_used") is False
            and item.get("admission_passed") is False
            and item.get("trajectory_advanced") is False
            and item.get("field_power_tail_ratios", {})
            .get("lambda_minus_1", [0.0, 0.0])[1]
            > 0.25
            and item.get("derivative_power_tail_ratios", {})
            .get("lambda_minus_1", [0.0, 0.0])[1]
            > 0.25
            for item in hlt7_events
        )
    )
    if (
        hlt7.get("classification")
        != "pre_trajectory_PROTO9_resolution_ladder_runtime_with_t0_nested_tail_obstruction"
        or hlt7.get("generated_by") != "scripts/reproduce_fgc_hlt7_mon7.py"
        or hlt7_lineage.get("checkpoint_commit")
        != "833910d1b31ed820745d61e0ce7976965e37b487"
        or hlt7_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or hlt7_protocol.get("artifact_id") != "FGC-2-SF1-PROTO9"
        or hlt7_protocol.get("protocol_version") != 9
        or hlt7_protocol.get("point_counts") != [2049, 4097, 8193]
        or hlt7_freeze.get("artifact_id") != "FGC-1-PRO9-FRZ1"
        or hlt7_freeze.get("claims_all_false") is not True
        or hlt7_freeze.get("result_sha256")
        != "cccda72803bd4b6bbdbc091a90e9746e1cc4b7952ef9f4dd9b684aee962b892d"
        or hlt7_predecessor.get("artifact_id") != "FGC-1-HLT6-MON6"
        or hlt7_predecessor.get("implementation_files_hash_matched") != 12
        or hlt7_predecessor.get("result_sha256")
        != "cda19a43fc08d6100511be046fb12a178adb29848429aa656a07160cf6787cbc"
        or hlt7_lineage.get("predecessor_run_plan_sha256")
        != "7967c1a785e9b6afcac3d516a426db4ed3652d58ab87fddb8b8874d824129e75"
        or hlt7_controls.get("old_PROTO8_ladder_rejected") is not True
        or hlt7_controls.get("reordered_ladder_rejected") is not True
        or hlt7_controls.get("conditional_projection_path_exists") is not False
        or hlt7_controls.get("exact_new_ladder_tested_on_physical_t0_inputs")
        is not True
        or hlt7_controls.get("all_controls_passed") is not True
        or hlt7_controls.get("inherited_PROTO8_adapter_controls", {}).get(
            "all_controls_passed"
        )
        is not True
        or hlt7_inherited.get(
            "PROTO8_common_event_and_PROTO7_evolution_unchanged"
        )
        is not True
        or hlt7_inherited.get("inherited_implementation_files_hash_matched")
        != 12
        or hlt7_inherited.get("process_local_runner_bindings_are_exception_safe")
        is not True
        or hlt7_plan.get("artifact_id") != "FGC-1-CAL6-RUN1-PLAN"
        or hlt7_plan.get("run_plan_sha256")
        != "ed0daf5be8353274131ad20ea13534fed491aa238a4234f07648e2a3bd6cafb3"
        or hlt7_plan.get("campaign_id")
        != "744258580cee7e7fb3c932470e354ad33d481bffb742ae3d3ff9c317ba68cb4c"
        or hlt7_plan.get("point_counts") != [2049, 4097, 8193]
        or hlt7_plan.get("expanded_input_count") != 12
        or hlt7_plan.get("new_8193_input_count") != 4
        or hlt7_plan.get("t0_common_event_evaluated_count") != 4
        or hlt7_plan.get("t0_common_event_passed_count") != 0
        or hlt7_inputs_valid is not True
        or len(hlt7_prechecks) != 12
        or any(
            item.get("accepted_state_source_gate_passed") is not True
            or item.get("trajectory_advanced") is not False
            or item.get("source_residual_infinity", 1.0)
            >= item.get("source_residual_maximum", 0.0)
            for item in hlt7_prechecks
        )
        or hlt7_events_valid is not True
        or hlt7_namespace.get("both_fresh_namespaces_empty") is not True
        or hlt7_namespace.get("authorization_created_no_namespace") is not True
        or [item.get("path") for item in hlt7_namespace.get("records", [])]
        != [
            "runs/fgc-2-sf1/proto9/calibration",
            "runs/fgc-2-sf1/proto9/holdout",
        ]
        or any(
            item.get("empty") is not True
            or item.get("created_by_authorization") is not False
            for item in hlt7_namespace.get("records", [])
        )
        or hlt7_payload.get("decision")
        != "HLT7_implements_the_PROTO9_resolution_binding_but_keeps_calibration_closed_by_the_t0_nested_tail_veto"
        or hlt7_boundary.get("PROTO9_trajectory_read") is not False
        or hlt7_boundary.get("PROTO9_trajectory_advanced") is not False
        or hlt7_boundary.get("conditional_8193_projection_used_as_evidence")
        is not False
        or hlt7_boundary.get("four_new_8193_t0_states_constructed") is not True
        or hlt7_boundary.get("fresh_GR0_recalibration_authorized") is not False
        or hlt7_boundary.get("t0_constraint_admission_passed_for_all_four_cases")
        is not True
        or hlt7_boundary.get(
            "t0_absolute_spatial_budgets_passed_for_all_four_cases"
        )
        is not True
        or hlt7_boundary.get("t0_nested_tail_admission_passed_for_any_case")
        is not False
        or hlt7_boundary.get("mechanism_question_answered") is not False
        or set(hlt7.get("implementation_sha256", {}))
        != expected_hlt7_implementation
        or set(hlt7.get("gate_status", {})) != expected_pro9_claims
        or hlt7.get("gate_status", {}).get("PROTO9_successor_runtime_implemented")
        is not True
        or any(
            value is not False
            for key, value in hlt7.get("gate_status", {}).items()
            if key != "PROTO9_successor_runtime_implemented"
        )
        or any(value is not False for value in hlt7.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-HLT7-MON7 lineage, direct-input, t0 obstruction, or nonclaim contract has drifted"
        )

    cal6_pref8_path = REPOSITORY / "results/fgc-1-cal6-pref8.json"
    try:
        cal6_pref8 = json.loads(
            cal6_pref8_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-CAL6-PREF8 certificate: {exc}")
        return failures
    pref8_payload = cal6_pref8.get("artifact_payload", {})
    pref8_immutable = pref8_payload.get("immutable_preflight", {})
    pref8_contract = pref8_payload.get("diagnostic_contract", {})
    pref8_cases = pref8_payload.get("case_diagnostics", [])
    pref8_aggregate = pref8_payload.get("aggregate", {})
    pref8_decision = pref8_payload.get("decision", {})
    pref8_revision = pref8_payload.get("prospective_revision", {})
    pref8_boundary = pref8_payload.get("epistemic_boundary", {})
    expected_pref8_fields = {
        "alpha_minus_1",
        "shift_over_r",
        "lambda_minus_1",
        "R_over_r_minus_1",
        "phi_over_Lambda",
        "chi_over_Lambda",
    }
    expected_pref8_case_order = [
        ("5/2", "RK4"),
        ("5/2", "SSPRK3"),
        ("3", "RK4"),
        ("3", "SSPRK3"),
    ]
    hlt7_input_hashes = {
        (item.get("amplitude"), item.get("method"), item.get("point_count")):
        item.get("projected_state_sha256")
        for item in hlt7_inputs
    }
    hlt7_events_by_case = {
        (item.get("amplitude"), item.get("method")): item
        for item in hlt7_events
    }
    pref8_failed_metrics = 0
    pref8_saturated_metrics = 0
    pref8_direct_metrics = 0
    pref8_cases_valid = (
        len(pref8_cases) == 4
        and [
            (item.get("amplitude"), item.get("method")) for item in pref8_cases
        ]
        == expected_pref8_case_order
    )
    if pref8_cases_valid:
        for case in pref8_cases:
            key = (case.get("amplitude"), case.get("method"))
            grids = case.get("grid_diagnostics", [])
            raw = case.get("raw_PROTO9_admission", {})
            convergence = case.get("whole_profile_convergence", {})
            successor = case.get("prospective_resolved_or_saturated_admission", {})
            raw_field = raw.get("field_power_tail_ratios", {})
            raw_derivative = raw.get("derivative_power_tail_ratios", {})
            field_classes = successor.get("finest_pair_field_classification", {})
            derivative_classes = successor.get(
                "finest_pair_derivative_classification", {}
            )
            if (
                case.get("point_counts") != [2049, 4097, 8193]
                or case.get("coordinate_time") != 0.0
                or case.get("trajectory_read") is not False
                or len(grids) != 3
                or [item.get("point_count") for item in grids]
                != [2049, 4097, 8193]
                or [item.get("spectral_sample_count") for item in grids]
                != [385, 769, 1537]
                or any(
                    item.get("projected_state_sha256")
                    != hlt7_input_hashes.get(
                        (key[0], key[1], item.get("point_count"))
                    )
                    or set(item.get("round_trip_interpolation_infinity_by_field", {}))
                    != expected_pref8_fields
                    or set(item.get("spectral_budgets_by_field", {}))
                    != expected_pref8_fields
                    or set(item.get("tail_sensitivity_by_field", {}))
                    != expected_pref8_fields
                    for item in grids
                )
                or raw.get("admission_passed") is not False
                or raw.get("every_nested_tail_ratio_passed") is not False
                or raw.get("finest_pair_individual_budgets_passed") is not True
                or raw_field
                != hlt7_events_by_case.get(key, {}).get("field_power_tail_ratios")
                or raw_derivative
                != hlt7_events_by_case.get(key, {}).get(
                    "derivative_power_tail_ratios"
                )
                or set(raw_field) != expected_pref8_fields
                or set(raw_derivative) != expected_pref8_fields
                or convergence.get("every_field_contracted") is not True
                or any(
                    value is not True
                    for value in convergence.get(
                        "complete_profile_contraction_by_field", {}
                    ).values()
                )
                or len(convergence.get("adjacent_differences", [])) != 2
                or successor.get("admission_passed") is not True
                or successor.get("saturation_used") is not True
                or successor.get("every_profile_contracted") is not True
                or successor.get("every_tail_resolved_or_saturated") is not True
                or successor.get("finest_pair_individual_budgets_passed") is not True
                or any(
                    value is not True
                    for value in successor.get(
                        "coarse_to_medium_direct_passed_by_field", {}
                    ).values()
                )
                or any(
                    value is not True
                    for value in successor.get(
                        "coarse_to_medium_direct_passed_by_derivative", {}
                    ).values()
                )
                or set(field_classes) != expected_pref8_fields
                or set(derivative_classes) != expected_pref8_fields
                or any(
                    classification
                    not in {"directly_resolved", "diagnostically_saturated"}
                    for classification in (*field_classes.values(), *derivative_classes.values())
                )
            ):
                pref8_cases_valid = False
                break
            for grid in grids[1:]:
                for name in expected_pref8_fields:
                    budget = grid["spectral_budgets_by_field"][name]
                    witness = grid["tail_sensitivity_by_field"][name]
                    round_trip = grid["round_trip_interpolation_infinity_by_field"][name]
                    if budget.get("individual_admission_passed") is not True:
                        pref8_cases_valid = False
                    if round_trip > 0.0 and (
                        witness.get("diagnostically_saturated") is not True
                        or witness.get("erasure_within_round_trip_scale") is not True
                        or witness.get("field_tail_within_round_trip_scale") is not True
                        or witness.get(
                            "derivative_tail_within_nyquist_round_trip_scale"
                        )
                        is not True
                        or witness.get("top_band_erasure_perturbation_infinity", 1.0)
                        > round_trip
                    ):
                        pref8_cases_valid = False
            for name in expected_pref8_fields:
                for ratios, classes in (
                    (raw_field, field_classes),
                    (raw_derivative, derivative_classes),
                ):
                    finest = ratios[name][1]
                    if finest is None or finest >= 0.25:
                        pref8_failed_metrics += 1
                    if classes[name] == "diagnostically_saturated":
                        pref8_saturated_metrics += 1
                    elif classes[name] == "directly_resolved":
                        pref8_direct_metrics += 1
    expected_pref8_implementation = {
        "scripts/reproduce_fgc_cal6_pref8.py",
        "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py",
        "src/recursive_horizons/fgc/evolution/proto4_admission.py",
        "src/recursive_horizons/fgc/evolution/proto8_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto9_runtime.py",
        "src/recursive_horizons/fgc/evolution/spectral_sensitivity.py",
    }
    expected_pref8_true_claims = {
        "PROTO9_preflight_conditioning_diagnosis_completed",
        "PROTO10_spectral_contract_revision_required",
    }
    expected_pref8_false_claims = {
        "PROTO9_ratio_only_failure_is_candidate_or_gradient_evidence",
        "PROTO10_frozen",
        "PROTO10_successor_runtime_implemented",
        "PROTO10_fresh_GR0_dynamic_calibration_authorized",
        "classical_spherical_diagnostic_authorized",
        "FGCQR_holdout_execution_authorized",
        "FGCQR_mechanism_rejected",
        "retained_EFT_evolution_authorized",
        "physical_transition_claim_authorized",
        "general_gradient_route_rejected",
        "singularity_resolution_derived",
        "child_domain_or_topology_derived",
        "dark_sector_mechanism_derived",
        "varying_locally_measured_c_derived",
    }
    pref8_gates = cal6_pref8.get("gate_status", {})
    if (
        cal6_pref8.get("classification")
        != "pre_trajectory_proper_spectral_ratio_conditioning_diagnosis"
        or cal6_pref8.get("generated_by")
        != "scripts/reproduce_fgc_cal6_pref8.py"
        or pref8_immutable.get("checkpoint_commit")
        != "234c5a5cc2c0b6bdc16d22146cdb9a9c6ae5b98d"
        or pref8_immutable.get("checkpoint_is_ancestor_of_HEAD") is not True
        or pref8_immutable.get("artifact_id") != "FGC-1-HLT7-MON7"
        or pref8_immutable.get("PROTO9_successor_runtime_implemented") is not True
        or pref8_immutable.get("PROTO9_fresh_GR0_dynamic_calibration_authorized")
        is not False
        or pref8_immutable.get("trajectory_read") is not False
        or pref8_contract.get("point_counts") != [2049, 4097, 8193]
        or set(pref8_contract.get("field_order", [])) != expected_pref8_fields
        or pref8_contract.get(
            "coarse_to_medium_field_and_derivative_ratios_must_pass_directly"
        )
        is not True
        or pref8_contract.get(
            "round_trip_interpolation_diagnostic_is_not_a_continuum_error_bound"
        )
        is not True
        or pref8_contract.get("no_epsilon_floor_or_new_numeric_tolerance_is_permitted")
        is not True
        or pref8_cases_valid is not True
        or pref8_failed_metrics != 16
        or pref8_saturated_metrics != 16
        or pref8_direct_metrics != 32
        or pref8_aggregate.get("case_count") != 4
        or pref8_aggregate.get("raw_failed_finest_pair_metric_count") != 16
        or pref8_aggregate.get(
            "prospective_diagnostically_saturated_metric_count"
        )
        != 16
        or pref8_aggregate.get(
            "prospective_directly_resolved_finest_pair_metric_count"
        )
        != 32
        or pref8_aggregate.get("all_raw_PROTO9_admissions_failed") is not True
        or pref8_aggregate.get(
            "all_prospective_resolved_or_saturated_admissions_passed"
        )
        is not True
        or not 0.0
        < pref8_aggregate.get("maximum_top_band_erasure_over_round_trip", 1.0)
        < 0.006
        or not 0.0
        < pref8_aggregate.get(
            "maximum_complete_profile_finest_pair_difference_ratio", 1.0
        )
        < 1.0
        or not 0.0
        < pref8_aggregate.get("maximum_round_trip_adjacent_ratio", 1.0)
        < 1.0
        or pref8_aggregate.get(
            "minimum_medium_or_fine_field_absolute_budget_headroom", 0.0
        )
        < 8.0e5
        or pref8_aggregate.get(
            "minimum_medium_or_fine_derivative_absolute_budget_headroom", 0.0
        )
        < 7.0e4
        or pref8_decision.get(
            "PROTO9_ratio_only_failure_is_candidate_or_gradient_evidence"
        )
        is not False
        or pref8_decision.get(
            "HLT7_finest_pair_ratio_is_conditioned_below_the_diagnostic_map_scale"
        )
        is not True
        or pref8_decision.get("new_protocol_required") != "FGC-2-SF1-PROTO10"
        or pref8_decision.get("current_campaign_authorized") is not False
        or pref8_revision.get("new_protocol_required") != "FGC-2-SF1-PROTO10"
        or pref8_revision.get("only_finest_pair_nested_tail_interpretation_may_change")
        is not True
        or pref8_revision.get("resolution_ladder_unchanged") is not True
        or pref8_revision.get("absolute_spectral_budgets_and_direct_ratio_ceiling_unchanged")
        is not True
        or pref8_revision.get("coarse_to_medium_direct_ratio_veto_unchanged")
        is not True
        or pref8_boundary
        != {
            "diagnostic_saturation_proves_physical_resolution": False,
            "mechanism_question_answered": False,
            "round_trip_is_a_continuum_error_bound": False,
            "tail_erasure_is_a_physical_field_perturbation": False,
            "trajectory_or_candidate_outcome_read": False,
        }
        or set(cal6_pref8.get("implementation_sha256", {}))
        != expected_pref8_implementation
        or set(pref8_gates) != expected_pref8_true_claims | expected_pref8_false_claims
        or any(pref8_gates.get(name) is not True for name in expected_pref8_true_claims)
        or any(pref8_gates.get(name) is not False for name in expected_pref8_false_claims)
        or any(value is not False for value in cal6_pref8.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-CAL6-PREF8 immutable lineage, raw-veto, conditioning, successor, or nonclaim contract has drifted"
        )

    pro10_path = REPOSITORY / "results/fgc-1-pro10-frz1.json"
    try:
        pro10 = json.loads(
            pro10_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-PRO10-FRZ1 certificate: {exc}")
        return failures
    pro10_payload = pro10.get("artifact_payload", {})
    pro10_lineage = pro10_payload.get("immutable_lineage", {})
    pro10_predecessor = pro10_lineage.get("predecessor_protocol", {})
    pro10_diagnosis = pro10_lineage.get("diagnosis", {})
    pro10_protocol = pro10_payload.get("protocol_validation", {})
    pro10_revision = pro10_payload.get("spectral_revision", {})
    pro10_witness = pro10_revision.get("immutable_design_witness", {})
    pro10_namespace = pro10_payload.get("namespace_precondition", {})
    pro10_boundary = pro10_payload.get("epistemic_boundary", {})
    expected_pro10_implementation = {
        "scripts/reproduce_fgc_pro10_frz1.py",
        "src/recursive_horizons/fgc/evolution/protocol_v10.py",
    }
    expected_pro10_claims = {
        "PROTO10_successor_runtime_implemented",
        "PROTO10_fresh_GR0_dynamic_calibration_authorized",
        "PROTO10_resolved_holdout_manifest_authorized",
        "classical_spherical_diagnostic_authorized",
        "FGCQR_holdout_execution_authorized",
        "retained_EFT_evolution_authorized",
        "physical_transition_claim_authorized",
        "FGCQR_mechanism_rejected",
        "general_gradient_route_rejected",
        "singularity_resolution_derived",
        "child_domain_or_topology_derived",
        "dark_sector_mechanism_derived",
        "varying_locally_measured_c_derived",
    }
    expected_pro10_namespaces = [
        {"path": "runs/fgc-2-sf1/proto10/calibration", "exists": False},
        {"path": "runs/fgc-2-sf1/proto10/holdout", "exists": False},
    ]
    if (
        pro10.get("classification")
        != "outcome_neutral_PROTO10_finest_pair_spectral_interpretation_freeze"
        or pro10.get("generated_by") != "scripts/reproduce_fgc_pro10_frz1.py"
        or pro10_lineage.get("checkpoint_commit")
        != "dd3f5adb22f08aa7950b74cef2595aaca2e8cd6e"
        or pro10_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or pro10_lineage.get("predecessor_protocol_sha256")
        != "81ce8426b006518d58accae50e43fcbfe942d013dfd010a716ff265b57afb48e"
        or pro10_predecessor.get("artifact_id") != "FGC-2-SF1-PROTO9"
        or pro10_predecessor.get("protocol_version") != 9
        or pro10_predecessor.get("point_counts") != [2049, 4097, 8193]
        or pro10_predecessor.get("frozen") is not True
        or pro10_diagnosis.get("artifact_id") != "FGC-1-CAL6-PREF8"
        or pro10_diagnosis.get("classification")
        != "pre_trajectory_proper_spectral_ratio_conditioning_diagnosis"
        or pro10_diagnosis.get("config_sha256")
        != "940238188ac8c76527165be6d48fcd1ef0e9b84565a9fe7a7e53b12a8e8559ea"
        or pro10_diagnosis.get("result_sha256")
        != "3d4624cca60363c9daf6f65f9393d959d7a7fafbcf6eac7f8dc9c2f237c24f5e"
        or pro10_diagnosis.get("case_count") != 4
        or pro10_diagnosis.get("raw_failed_finest_pair_metric_count") != 16
        or pro10_diagnosis.get("diagnostically_saturated_metric_count") != 16
        or pro10_diagnosis.get("directly_resolved_finest_pair_metric_count") != 32
        or pro10_diagnosis.get("all_prospective_design_admissions_passed")
        is not True
        or pro10_diagnosis.get("current_campaign_authorized") is not False
        or pro10_diagnosis.get("claims_candidate_and_physical_remain_false")
        is not True
        or not 0.0
        < pro10_diagnosis.get("maximum_top_band_erasure_over_round_trip", 1.0)
        < 0.006
        or not 0.0
        < pro10_diagnosis.get(
            "maximum_complete_profile_finest_pair_difference_ratio", 1.0
        )
        < 1.0
        or not 0.0
        < pro10_diagnosis.get("maximum_round_trip_adjacent_ratio", 1.0)
        < 1.0
        or pro10_protocol.get("artifact_id") != "FGC-2-SF1-PROTO10"
        or pro10_protocol.get("protocol_version") != 10
        or pro10_protocol.get("frozen") is not True
        or pro10_protocol.get("outcome_neutral_contract_validated") is not True
        or pro10_protocol.get("premise_revision_only") is not True
        or pro10_protocol.get("predecessor_protocol_artifact_id")
        != "FGC-2-SF1-PROTO9"
        or pro10_protocol.get("diagnosis_artifact_id") != "FGC-1-CAL6-PREF8"
        or pro10_protocol.get("diagnosis_checkpoint_commit")
        != "dd3f5adb22f08aa7950b74cef2595aaca2e8cd6e"
        or pro10_protocol.get("eligible_calibration_amplitudes") != ["5/2", "3"]
        or pro10_protocol.get("point_counts") != [2049, 4097, 8193]
        or pro10_protocol.get("resolution_ladder_unchanged") is not True
        or pro10_protocol.get("only_finest_pair_spectral_interpretation_changes")
        is not True
        or pro10_protocol.get("maximum_nested_tail_ratio") != "1/4"
        or pro10_protocol.get("coarse_to_medium_direct_ratio_veto") is not True
        or pro10_protocol.get("medium_and_fine_absolute_budgets_unchanged")
        is not True
        or pro10_protocol.get("finest_pair_direct_ratio_remains_primary")
        is not True
        or pro10_protocol.get("guarded_saturation_requires_tail_erasure_map_witness")
        is not True
        or pro10_protocol.get("guarded_saturation_requires_round_trip_contraction")
        is not True
        or pro10_protocol.get(
            "guarded_saturation_requires_complete_profile_contraction"
        )
        is not True
        or pro10_protocol.get("round_trip_is_a_continuum_error_bound") is not False
        or pro10_protocol.get("diagnostic_saturation_is_physical_resolution_evidence")
        is not False
        or pro10_protocol.get("new_numeric_tolerance_added") is not False
        or pro10_protocol.get("physical_and_numerical_thresholds_unchanged")
        is not True
        or pro10_protocol.get("resolved_holdout_manifest")
        != "configs/fgc/fgc-1-pro10-hld1.toml"
        or pro10_protocol.get("resolved_holdout_manifest_present") is not False
        or pro10_protocol.get("calibration_output_root")
        != "runs/fgc-2-sf1/proto10/calibration"
        or pro10_protocol.get("holdout_output_root")
        != "runs/fgc-2-sf1/proto10/holdout"
        or set(pro10_protocol.get("claims", {})) != expected_pro10_claims
        or any(value is not False for value in pro10_protocol.get("claims", {}).values())
        or pro10_revision.get("point_counts") != [2049, 4097, 8193]
        or pro10_revision.get("resolution_ladder_unchanged") is not True
        or pro10_revision.get("maximum_nested_tail_ratio") != "1/4"
        or pro10_revision.get("raw_ratio_ceiling_unchanged") is not True
        or pro10_revision.get("coarse_to_medium_direct_ratio_veto_unchanged")
        is not True
        or pro10_revision.get("medium_and_fine_absolute_budgets_unchanged")
        is not True
        or pro10_revision.get("finest_pair_direct_route_remains_primary")
        is not True
        or pro10_revision.get("finest_pair_guarded_saturation_route_added")
        is not True
        or pro10_revision.get("tail_erasure_map_witness_required") is not True
        or pro10_revision.get("round_trip_contraction_required") is not True
        or pro10_revision.get(
            "complete_profile_infinity_and_RMS_contraction_required"
        )
        is not True
        or pro10_revision.get("round_trip_is_a_continuum_error_bound") is not False
        or pro10_revision.get("diagnostic_saturation_is_physical_resolution_evidence")
        is not False
        or pro10_revision.get("new_epsilon_floor_or_numeric_tolerance_added")
        is not False
        or pro10_witness.get("raw_failed_finest_pair_metric_count") != 16
        or pro10_witness.get("diagnostically_saturated_metric_count") != 16
        or pro10_witness.get("directly_resolved_finest_pair_metric_count") != 32
        or pro10_witness.get("is_calibration_or_mechanism_evidence") is not False
        or pro10_witness.get("maximum_top_band_erasure_over_round_trip")
        != pro10_diagnosis.get("maximum_top_band_erasure_over_round_trip")
        or pro10_witness.get("maximum_complete_profile_finest_pair_difference_ratio")
        != pro10_diagnosis.get(
            "maximum_complete_profile_finest_pair_difference_ratio"
        )
        or pro10_witness.get("maximum_round_trip_adjacent_ratio")
        != pro10_diagnosis.get("maximum_round_trip_adjacent_ratio")
        or pro10_namespace.get("records") != expected_pro10_namespaces
        or pro10_namespace.get("both_fresh_namespaces_absent") is not True
        or pro10_namespace.get("freeze_created_no_namespace") is not True
        or pro10_payload.get("decision")
        != "PROTO10_is_frozen_but_HLT8_fresh_calibration_and_every_candidate_outcome_remain_closed"
        or pro10_boundary
        != {
            "GR0_calibration_completed": False,
            "PROTO10_trajectory_advanced": False,
            "PROTO10_trajectory_read": False,
            "SGBL_or_FGCQR_trajectory_read": False,
            "fresh_state_constructed": False,
            "mechanism_question_answered": False,
        }
        or set(pro10.get("implementation_sha256", {}))
        != expected_pro10_implementation
        or set(pro10.get("gate_status", {})) != expected_pro10_claims
        or any(value is not False for value in pro10.get("gate_status", {}).values())
        or any(value is not False for value in pro10.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-PRO10-FRZ1 immutable lineage, guarded spectral contract, namespace, or nonclaim boundary has drifted"
        )

    hlt8_path = REPOSITORY / "results/fgc-1-hlt8-mon8.json"
    try:
        hlt8 = json.loads(
            hlt8_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-HLT8-MON8 certificate: {exc}")
        return failures
    hlt8_payload = hlt8.get("artifact_payload", {})
    hlt8_lineage = hlt8_payload.get("immutable_lineage", {})
    hlt8_protocol = hlt8_lineage.get("protocol", {})
    hlt8_freeze = hlt8_lineage.get("freeze", {})
    hlt8_predecessor = hlt8_lineage.get("predecessor_runtime", {})
    hlt8_diagnosis = hlt8_lineage.get("diagnosis", {})
    hlt8_adapter = hlt8_payload.get("runtime_adapter_controls", {})
    hlt8_attacks = hlt8_payload.get("adversarial_guard_controls", {})
    hlt8_inherited = hlt8_payload.get("inherited_runtime", {})
    hlt8_plan = hlt8_payload.get("frozen_run_plan", {})
    hlt8_inputs = hlt8_payload.get("frozen_run_inputs", [])
    hlt8_prechecks = hlt8_payload.get("accepted_state_source_prechecks", [])
    hlt8_events = hlt8_payload.get("initial_common_events", [])
    hlt8_namespace = hlt8_payload.get("namespace_precondition", {})
    hlt8_boundary = hlt8_payload.get("epistemic_boundary", {})
    expected_hlt8_implementation = {
        "scripts/reproduce_fgc_hlt8_mon8.py",
        "scripts/run_fgc_gr0_calibration.py",
        "scripts/run_fgc_gr0_calibration_v6.py",
        "scripts/run_fgc_gr0_calibration_v7.py",
        "scripts/run_fgc_gr0_calibration_v8.py",
        "scripts/run_fgc_gr0_calibration_v9.py",
        "scripts/run_fgc_gr0_calibration_v10.py",
        "src/recursive_horizons/fgc/evolution/cal2_source_diagnosis.py",
        "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py",
        "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto5_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto6_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto7_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto8_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto9_runtime.py",
        "src/recursive_horizons/fgc/evolution/proto10_runtime.py",
        "src/recursive_horizons/fgc/evolution/protocol_v10.py",
        "src/recursive_horizons/fgc/evolution/spectral_sensitivity.py",
    }
    expected_hlt8_order = [
        (amplitude, method, count)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
        for count in (2049, 4097, 8193)
    ]
    predecessor_hlt7_hashes = {
        (item.get("amplitude"), item.get("method"), item.get("point_count")): item.get(
            "projected_state_sha256"
        )
        for item in hlt7_inputs
    }
    hlt8_inputs_valid = (
        len(hlt8_inputs) == 12
        and [
            (item.get("amplitude"), item.get("method"), item.get("point_count"))
            for item in hlt8_inputs
        ]
        == expected_hlt8_order
        and all(
            item.get("projected_state_sha256")
            == predecessor_hlt7_hashes.get(
                (item.get("amplitude"), item.get("method"), item.get("point_count"))
            )
            and item.get("projected_state_matches_HLT7") is True
            and item.get("diagnostic_saturation_is_state_mutation") is False
            and item.get("trajectory_advanced") is False
            for item in hlt8_inputs
        )
    )
    hlt7_event_map = {
        (item.get("amplitude"), item.get("method")): item for item in hlt7_events
    }
    hlt8_events_valid = (
        len(hlt8_events) == 4
        and {
            (item.get("amplitude"), item.get("method")) for item in hlt8_events
        }
        == {
            ("5/2", "RK4"),
            ("5/2", "SSPRK3"),
            ("3", "RK4"),
            ("3", "SSPRK3"),
        }
        and all(
            item.get("raw_PROTO9_event_projection")
            == hlt7_event_map.get((item.get("amplitude"), item.get("method")))
            and item.get("raw_PROTO9_admission_passed") is False
            and item.get("guarded_PROTO10_admission_passed") is True
            and item.get("diagnostic_saturation_used") is True
            and item.get("diagnostic_saturation_is_physical_resolution") is False
            and item.get("round_trip_is_a_continuum_error_bound") is False
            and item.get("trajectory_advanced") is False
            and item.get("guarded_PROTO10_common_event", {}).get(
                "raw_PROTO9_admission_passed"
            )
            is False
            and item.get("guarded_PROTO10_common_event", {}).get(
                "admission_passed"
            )
            is True
            and item.get("guarded_PROTO10_common_event", {})
            .get("spatial_spectral_admission", {})
            .get("saturation_used")
            is True
            for item in hlt8_events
        )
    )
    expected_hlt8_claims = expected_pro10_claims
    hlt8_gates = hlt8.get("gate_status", {})
    hlt8_namespace_records = hlt8_namespace.get("records", [])
    if (
        hlt8.get("classification")
        != "pre_trajectory_PROTO10_guarded_spectral_runtime_and_fresh_GR0_calibration_authorization"
        or hlt8.get("generated_by") != "scripts/reproduce_fgc_hlt8_mon8.py"
        or hlt8_lineage.get("checkpoint_commit")
        != "3f7ac298a027209e5bdb2d59c8c82e04c25211ab"
        or hlt8_lineage.get("checkpoint_is_ancestor_of_HEAD") is not True
        or hlt8_protocol.get("artifact_id") != "FGC-2-SF1-PROTO10"
        or hlt8_protocol.get("protocol_version") != 10
        or hlt8_protocol.get("point_counts") != [2049, 4097, 8193]
        or hlt8_protocol.get("maximum_nested_tail_ratio") != "1/4"
        or hlt8_freeze.get("artifact_id") != "FGC-1-PRO10-FRZ1"
        or hlt8_freeze.get("claims_all_false") is not True
        or hlt8_freeze.get("result_sha256")
        != "f807eabf24e1bfec3dd9a30cacf489ac812356ddd60a08881d42563cd3e17ae3"
        or hlt8_predecessor.get("artifact_id") != "FGC-1-HLT7-MON7"
        or hlt8_predecessor.get("raw_t0_event_count") != 4
        or hlt8_predecessor.get("result_sha256")
        != "e6eb5832302eea39f7adf7fee2996cf659979ebbc58317b8f6d3b7e092f789c5"
        or hlt8_diagnosis.get("artifact_id") != "FGC-1-CAL6-PREF8"
        or hlt8_diagnosis.get("raw_failed_finest_pair_metric_count") != 16
        or hlt8_diagnosis.get("prospective_saturated_metric_count") != 16
        or hlt8_lineage.get("predecessor_run_plan_sha256")
        != "ed0daf5be8353274131ad20ea13534fed491aa238a4234f07648e2a3bd6cafb3"
        or hlt8_lineage.get("inherited_implementation_files_hash_matched") != 20
        or hlt8_adapter.get("old_PROTO8_ladder_rejected") is not True
        or hlt8_adapter.get("reordered_ladder_rejected") is not True
        or hlt8_adapter.get("correct_ladder_direct_spectral_control_passed")
        is not True
        or hlt8_adapter.get("direct_finest_pair_route_used_saturation") is not False
        or hlt8_adapter.get("all_controls_passed") is not True
        or hlt8_attacks.get("coarse_direct_failure_admitted") is not False
        or hlt8_attacks.get("missing_tail_erasure_map_witness_admitted") is not False
        or hlt8_attacks.get("missing_complete_profile_contraction_admitted") is not False
        or hlt8_attacks.get("missing_round_trip_contraction_admitted") is not False
        or hlt8_attacks.get("maximum_nested_tail_ratio") != 0.25
        or hlt8_attacks.get("new_epsilon_floor_added") is not False
        or hlt8_attacks.get("all_adversarial_controls_passed") is not True
        or hlt8_inherited.get(
            "PROTO9_raw_common_event_and_PROTO7_evolution_unchanged"
        )
        is not True
        or hlt8_inherited.get("process_local_runner_bindings_are_exception_safe")
        is not True
        or hlt8_plan.get("artifact_id") != "FGC-1-CAL7-RUN1-PLAN"
        or hlt8_plan.get("run_plan_sha256")
        != sha256(
            (REPOSITORY / "configs/fgc/fgc-1-cal7-run1.toml").read_bytes()
        ).hexdigest()
        or not isinstance(hlt8_plan.get("campaign_id"), str)
        or len(hlt8_plan.get("campaign_id", "")) != 64
        or hlt8_plan.get("point_counts") != [2049, 4097, 8193]
        or hlt8_plan.get("expanded_input_count") != 12
        or hlt8_plan.get("projected_state_hash_match_count") != 12
        or hlt8_plan.get("raw_PROTO9_t0_passed_count") != 0
        or hlt8_plan.get("guarded_PROTO10_t0_passed_count") != 4
        or hlt8_plan.get("t0_common_event_evaluated_count") != 4
        or hlt8_inputs_valid is not True
        or len(hlt8_prechecks) != 12
        or any(
            item.get("accepted_state_source_gate_passed") is not True
            or not 0.0 <= item.get("source_residual_infinity", 1.0) < 1.0e-12
            or item.get("source_residual_maximum") != 1.0e-12
            or item.get("trajectory_advanced") is not False
            for item in hlt8_prechecks
        )
        or hlt8_events_valid is not True
        or hlt8_namespace.get("both_fresh_namespaces_empty") is not True
        or hlt8_namespace.get("authorization_created_no_namespace") is not True
        or [item.get("path") for item in hlt8_namespace_records]
        != [
            "runs/fgc-2-sf1/proto10/calibration",
            "runs/fgc-2-sf1/proto10/holdout",
        ]
        or any(
            item.get("empty") is not True
            or item.get("entry_count") != 0
            or item.get("created_by_authorization") is not False
            for item in hlt8_namespace_records
        )
        or hlt8_payload.get("decision")
        != "HLT8_implements_the_PROTO10_compositor_and_authorizes_one_fresh_GR0_calibration_but_keeps_every_candidate_and_physical_claim_closed"
        or hlt8_boundary.get("PROTO10_trajectory_read") is not False
        or hlt8_boundary.get("PROTO10_trajectory_advanced") is not False
        or hlt8_boundary.get("all_twelve_source_prechecks_passed") is not True
        or hlt8_boundary.get(
            "all_four_raw_PROTO9_admissions_failed_and_remain_public"
        )
        is not True
        or hlt8_boundary.get("all_four_guarded_PROTO10_admissions_passed")
        is not True
        or hlt8_boundary.get("fresh_GR0_recalibration_authorized") is not True
        or hlt8_boundary.get("fresh_GR0_recalibration_completed") is not False
        or hlt8_boundary.get("diagnostic_saturation_is_physical_resolution")
        is not False
        or hlt8_boundary.get("round_trip_is_a_continuum_error_bound") is not False
        or hlt8_boundary.get("mechanism_question_answered") is not False
        or set(hlt8.get("implementation_sha256", {}))
        != expected_hlt8_implementation
        or set(hlt8_gates) != expected_hlt8_claims
        or hlt8_gates.get("PROTO10_successor_runtime_implemented") is not True
        or hlt8_gates.get("PROTO10_fresh_GR0_dynamic_calibration_authorized")
        is not True
        or any(
            value is not False
            for key, value in hlt8_gates.items()
            if key
            not in {
                "PROTO10_successor_runtime_implemented",
                "PROTO10_fresh_GR0_dynamic_calibration_authorized",
            }
        )
        or any(value is not False for value in hlt8.get("nonclaims", {}).values())
    ):
        failures.append(
            "FGC-1-HLT8-MON8 lineage, raw/guarded compositor, authorization, or nonclaim contract has drifted"
        )

    ctr1_path = REPOSITORY / "results/fgc-1-ctr1-reg1.json"
    try:
        ctr1 = json.loads(
            ctr1_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-CTR1-REG1 certificate: {exc}")
        return failures
    ctr1_payload = ctr1.get("artifact_payload", {})
    ctr1_quantitative = ctr1_payload.get("quantitative_evidence", {})
    ctr1_minkowski = ctr1_quantitative.get("Minkowski_control", {})
    ctr1_smooth = ctr1_quantitative.get("smooth_nontrivial_control", {})
    ctr1_formal = ctr1_smooth.get("formal_certificate", {})
    ctr1_curvature = ctr1_formal.get("curvature_series", {})
    ctr1_defects = ctr1_quantitative.get("negative_controls", {})
    ctr1_conical = ctr1_defects.get("conical_defect", {})
    ctr1_odd = ctr1_defects.get("odd_scalar_parity", {})
    ctr1_first_grid = ctr1_smooth.get("first_grid_point_convergence", {})
    required_ctr1_payload = {
        "R_equals_r_times_A_and_v_equals_r_times_V",
        "all_removable_singular_limits_derived_analytically",
        "negative_power_guard_edge_fails_closed",
        "finite_curvature_series_controls_passed",
        "first_grid_point_convergence_passed",
        "parity_and_elementary_flatness_enforced",
    }
    if (
        ctr1.get("classification") != "regular_center_formulation_certificate"
        or ctr1.get("generated_by") != "scripts/reproduce_fgc_ctr1_reg1.py"
        or any(ctr1_payload.get(key) is not True for key in required_ctr1_payload)
        or ctr1_minkowski.get("exact_zero_residual_gauge_and_curvature") is not True
        or ctr1_minkowski.get("annular_regular_equations") != ["0"] * 6
        or ctr1_minkowski.get("annular_regular_gauge") != ["0"] * 2
        or ctr1_formal.get("center_limits")
        != [
            "-9293007317/2822400000",
            "61885569858289/70785792000000",
            "-14836176781/62092800000",
            "-14836176781/62092800000",
            "-1178969917/2956800000",
            "247/4500",
        ]
        or ctr1_smooth.get("center_radial_angular_metric_isotropy_exact") is not True
        or ctr1_smooth.get("nonzero_finite_Kretschmann_control") is not True
        or ctr1_curvature.get("Ricci_scalar", {}).get("0") != "2701/2400"
        or ctr1_curvature.get("Kretschmann", {}).get("0") != "18910513/17280000"
        or ctr1_conical.get("rejected_before_certificate") is not True
        or ctr1_conical.get("analytic_ledger", {}).get(
            "Ricci_scalar_r_minus_2_from_conical_defect"
        )
        != "-42/121"
        or ctr1_odd.get("rejected_before_certificate") is not True
        or ctr1_odd.get("analytic_ledger", {}).get(
            "box_phi_r_minus_1_from_odd_scalar_term"
        )
        != "2/17"
        or ctr1_first_grid.get("all_regular_equation_and_gauge_controls_passed")
        is not True
        or ctr1_first_grid.get("not_a_PDE_discretization_convergence_claim")
        is not True
    ):
        failures.append("FGC-1-CTR1-REG1 exact centre/defect/convergence contract has drifted")
    if _contains_float(ctr1):
        failures.append("FGC-1-CTR1-REG1 certificate contains floating-point data")

    run1_path = REPOSITORY / "results/fgc-1-run1-sym1.json"
    try:
        run1 = json.loads(
            run1_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-RUN1-SYM1 decision: {exc}")
        return failures
    run1_audit = run1.get("scoped_run_authorization_audit", {})
    scoped = run1_audit.get("authorization_audits", {}).get(
        "classical_spherical_diagnostic", {}
    )
    if (
        scoped.get("passed_predicate_count") != 7
        or scoped.get("required_predicate_count") != 8
        or scoped.get("authorized") is not False
        or scoped.get("missing_predicate_ids")
        != ["frozen_outcome_neutral_protocol_and_resolved_holdout_hash"]
    ):
        failures.append("FGC-1-RUN1-SYM1 must retain its current 7/8 scoped stop")
    protocol = run1_audit.get("protocol_validation", {})
    if (
        protocol.get("artifact_id") != "FGC-2-SF1-PROTO3"
        or protocol.get("protocol_version") != 3
        or protocol.get("predecessor_protocol_artifact_id")
        != "FGC-2-SF1-PROTO2"
        or protocol.get("preflight_artifact_id") != "FGC-1-ID0-PREF1"
        or protocol.get("diagnosis_artifact_id") != "FGC-1-ID1-FAM1"
        or protocol.get("amendment_validated") is not True
        or protocol.get("PROTO2_amendment_validated") is not True
        or protocol.get("PROTO1_and_PROTO2_preserved") is not True
        or protocol.get("premise_revision_only") is not True
        or protocol.get("FGCQR_outcomes_inspected_before_revision") is not False
        or protocol.get("SGBL_outcomes_inspected_before_revision") is not False
        or protocol.get("explicit_phi_unit_normal_momentum") != "0"
        or protocol.get("amplitude_candidates")
        != ["2", "5/2", "3", "7/2", "4", "9/2", "5"]
        or protocol.get("exact_widest_center_buffer_over_L0") != "39/16"
        or protocol.get("widest_case_coarsest_grid_intervals") != 78
        or protocol.get(
            "static_all_case_initial_premises_required_before_GR0_dynamic_eligibility"
        )
        is not True
        or protocol.get("semantic_holdout_contract_sha256")
        != "d18950ffea629e64c3b67d5f6fec7f989255d7db4e9b8aecdd681988d64ca9ee"
        or protocol.get("outcome_neutral_contract_validated") is not True
        or protocol.get("held_out_outcomes_inspected_before_freeze") is not False
        or protocol.get("resolved_holdout_manifest_present") is not False
    ):
        failures.append("FGC-1-RUN1-SYM1 protocol boundary has drifted")
    diagnosis = run1.get("protocol_amendment_diagnosis", {})
    if (
        diagnosis.get("artifact_id") != "FGC-1-ID1-FAM1"
        or diagnosis.get("historical_checkpoint_commit")
        != "fae2c5b76ba53680f9e70188f79784b160826734"
        or diagnosis.get("immutable_PROTO2_config_sha256")
        != "091d6d7dc8109d21d9297ea6eb37a7296e0ea45fc5e05269bba1b488159f09ba"
        or diagnosis.get("immutable_ID1_result_sha256")
        != "aeda6543854ca54456dc58648ffc6824be297eb3930f0fbe7c8b8beae5e8075f"
        or diagnosis.get(
            "PROTO2_width_nine_eighths_case_unqualified_at_checkpoint"
        )
        is not True
        or diagnosis.get("PROTO3_revision_is_premise_only") is not True
        or diagnosis.get("FGCQR_evolution_outcomes_inspected_before_revision")
        is not False
        or diagnosis.get("SGBL_evolution_outcomes_inspected_before_revision")
        is not False
    ):
        failures.append("FGC-1-RUN1-SYM1 immutable PROTO3 diagnosis has drifted")
    preflight = run1.get("protocol_amendment_preflight", {})
    if (
        preflight.get("artifact_id") != "FGC-1-ID0-PREF1"
        or preflight.get("historical_checkpoint_commit")
        != "d4f0cc8f58408ef4ee12fb32fe3619916e231795"
        or preflight.get("PROTO1_initial_data_preflight_obstruction_verified")
        is not True
        or preflight.get("PROTO2_revision_is_premise_only") is not True
        or preflight.get("PROTO3_preserves_PROTO2_ID0_lineage") is not True
        or preflight.get("FGCQR_outcomes_inspected_before_revision") is not False
        or preflight.get("SGBL_outcomes_inspected_before_revision") is not False
    ):
        failures.append("FGC-1-RUN1-SYM1 protocol-amendment preflight has drifted")
    future = run1_audit.get("future_evidence_status", {})
    expected_future = {
        "protocol_holdout",
        "run_domain",
        "multidirectional_health",
        "nonlinear_source",
        "physical_constraints",
        "regular_center",
        "initial_data",
        "boundary_control",
        "health_monitor",
        "numerical_validation",
    }
    nonlinear = future.get("nonlinear_source", {})
    constraints = future.get("physical_constraints", {})
    regular_center = future.get("regular_center", {})
    initial_data = future.get("initial_data", {})
    run_domain = future.get("run_domain", {})
    multidirectional = future.get("multidirectional_health", {})
    boundary_control = future.get("boundary_control", {})
    health_monitor = future.get("health_monitor", {})
    numerical_validation = future.get("numerical_validation", {})
    protocol_holdout = future.get("protocol_holdout", {})
    absent = {
        key: value
        for key, value in future.items()
        if key
        not in {
            "nonlinear_source",
            "physical_constraints",
            "regular_center",
            "initial_data",
            "run_domain",
            "multidirectional_health",
            "boundary_control",
            "health_monitor",
            "numerical_validation",
        }
    }
    if (
        set(future) != expected_future
        or not isinstance(nonlinear, dict)
        or nonlinear.get("status") != "scope_bound_hash_bound_gate_passed"
        or nonlinear.get("passed") is not True
        or nonlinear.get("declared_run_envelope_sha256")
        != protocol.get("semantic_holdout_contract_sha256")
        or not isinstance(constraints, dict)
        or constraints.get("status") != "scope_bound_hash_bound_gate_passed"
        or constraints.get("passed") is not True
        or constraints.get("declared_run_envelope_sha256")
        != protocol.get("semantic_holdout_contract_sha256")
        or not isinstance(regular_center, dict)
        or regular_center.get("status") != "scope_bound_hash_bound_gate_passed"
        or regular_center.get("passed") is not True
        or regular_center.get("declared_run_envelope_sha256")
        != protocol.get("semantic_holdout_contract_sha256")
        or not isinstance(initial_data, dict)
        or initial_data.get("status") != "scope_bound_hash_bound_gate_passed"
        or initial_data.get("passed") is not True
        or initial_data.get("declared_run_envelope_sha256")
        != protocol.get("semantic_holdout_contract_sha256")
        or not isinstance(run_domain, dict)
        or run_domain.get("status") != "scope_bound_hash_bound_gate_passed"
        or run_domain.get("passed") is not True
        or run_domain.get("declared_run_envelope_sha256")
        != protocol.get("semantic_holdout_contract_sha256")
        or not isinstance(multidirectional, dict)
        or multidirectional.get("status")
        != "scope_bound_hash_bound_gate_passed"
        or multidirectional.get("passed") is not True
        or multidirectional.get("declared_run_envelope_sha256")
        != protocol.get("semantic_holdout_contract_sha256")
        or not isinstance(boundary_control, dict)
        or boundary_control.get("status") != "scope_bound_hash_bound_gate_passed"
        or boundary_control.get("passed") is not True
        or boundary_control.get("declared_run_envelope_sha256")
        != protocol.get("semantic_holdout_contract_sha256")
        or not isinstance(health_monitor, dict)
        or health_monitor.get("status") != "scope_bound_hash_bound_gate_passed"
        or health_monitor.get("passed") is not True
        or health_monitor.get("declared_run_envelope_sha256")
        != protocol.get("semantic_holdout_contract_sha256")
        or not isinstance(numerical_validation, dict)
        or numerical_validation.get("status")
        != "scope_bound_hash_bound_gate_passed"
        or numerical_validation.get("passed") is not True
        or numerical_validation.get("declared_run_envelope_sha256")
        != protocol.get("semantic_holdout_contract_sha256")
        or not isinstance(protocol_holdout, dict)
        or protocol_holdout.get("artifact_id") != "FGC-1-PRO3-HLD1"
        or protocol_holdout.get("status") != "absent"
        or protocol_holdout.get("passed") is not False
        or any(
            not isinstance(value, dict)
            or value.get("status") != "absent"
            or value.get("passed") is not False
            for value in absent.values()
        )
    ):
        failures.append(
            "FGC-1-RUN1-SYM1 must compose DOM4, HYP2, SRC1, CON4, CTR1, ID1, BND2, HLT1, and NUM1 while retaining only PRO3-HLD1 absent"
        )

    id0_path = REPOSITORY / "results/fgc-1-id0-pref1.json"
    try:
        id0 = json.loads(
            id0_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-ID0-PREF1 certificate: {exc}")
        return failures
    id0_payload = id0.get("artifact_payload", {})
    id0_diagnosis = id0_payload.get("protocol_diagnosis", {})
    id0_exact = id0_payload.get("exact_constraint_specialization", {})
    id0_analytic = id0_payload.get("analytic_entire_box_obstruction", {})
    id0_bound = id0_analytic.get("bound", {})
    id0_comparisons = id0_analytic.get("comparisons", {})
    id0_numerical = id0_payload.get("independent_numerical_preflight", {})
    id0_original = id0_numerical.get("original_box", {})
    id0_repair = id0_numerical.get("premise_repair_scan", {})
    id0_boundary = id0_payload.get("epistemic_boundary", {})
    id0_history = id0.get("historical_preflight_checkpoint", {})
    id0_contract = id0.get("certificate_contract", {})
    exact_fixtures = id0_exact.get("fixtures", {})
    repair_records = id0_repair.get("records", [])
    maximum_original = id0_original.get("maximum_peak_compactness")
    if (
        id0.get("classification")
        != "pre_holdout_GR0_protocol_obstruction_not_an_FGCQR_mechanism_result"
        or id0.get("generated_by") != "scripts/reproduce_fgc_id0_pref1.py"
        or id0_diagnosis.get("phi_unit_normal_momentum_frozen_by_PROTO1") is not False
        or id0_diagnosis.get("minimal_preflight_completion") != "Pi_phi=0"
        or id0_diagnosis.get(
            "missing_phi_momentum_alone_prevents_executable_initial_data_contract"
        )
        is not True
        or id0_exact.get("fixture_count") != 2
        or id0_exact.get("all_exact_pair_equalities_passed") is not True
        or set(exact_fixtures) != {"fixture_A", "fixture_B"}
        or any(
            not isinstance(fixture, dict)
            or fixture.get("exact_pair_equality") is not True
            or fixture.get("source_is_unredefined_ACT1_VAR1") is not True
            for fixture in exact_fixtures.values()
        )
        or id0_bound.get("compactness_absolute_bound")
        != "3927454576397784264080597/276701161105643274240000000"
        or id0_bound.get("inverse_metric_bootstrap_closed") is not True
        or id0_comparisons.get("compactness_strictly_below_separator") is not True
        or id0_comparisons.get("separator_strictly_below_protocol_floor") is not True
        or not isinstance(maximum_original, (int, float))
        or not isfinite(maximum_original)
        or not 0.00042 < maximum_original < 0.00044
        or id0_original.get("all_candidates_below_declared_upper_bound") is not True
        or len(repair_records) != 7
        or id0_repair.get("all_candidates_inside_original_initial_compactness_window")
        is not True
        or id0_repair.get("FGCQR_outcomes_used") is not False
        or id0_boundary.get("protocol_design_route_closed") != "PROTO1_as_written"
        or id0_boundary.get("FGCQR_action_or_mechanism_tested") is not False
        or id0_boundary.get("FGCQR_action_or_mechanism_rejected") is not False
        or id0_boundary.get("general_gradient_mechanism_rejected") is not False
        or id0_history.get("commit")
        != "d4f0cc8f58408ef4ee12fb32fe3619916e231795"
        or id0_history.get("commit_is_ancestor_of_HEAD") is not True
        or id0_history.get("original_ID0_canonical_ledgers_verified_at_commit")
        is not True
        or id0_history.get("original_RUN1_two_of_eight_stop_verified_at_commit")
        is not True
        or id0_history.get("FGCQR_holdout_outcomes_inspected_before_ID0")
        is not False
        or id0_contract.get("historical_pre_PROTO2_checkpoint_verified") is not True
    ):
        failures.append("FGC-1-ID0-PREF1 exact protocol-obstruction boundary has drifted")

    id1_path = REPOSITORY / "results/fgc-1-id1-fam1.json"
    try:
        id1 = json.loads(
            id1_path.read_text(encoding="utf-8"),
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (OSError, json.JSONDecodeError, ValueError) as exc:
        failures.append(f"cannot inspect FGC-1-ID1-FAM1 certificate: {exc}")
        return failures
    id1_payload = id1.get("artifact_payload", {})
    id1_quantitative = id1_payload.get("quantitative_evidence", {})
    id1_exact = id1_quantitative.get("exact_full_evaluator_controls", [])
    id1_gr0 = id1_quantitative.get("exact_zero_seed_GR0_reduction", {})
    id1_convergence = id1_quantitative.get("central_common_grid_convergence", {})
    id1_family = id1_quantitative.get("family_box", {})
    id1_extrema = id1_family.get("extrema", {})
    id1_inequalities = id1_family.get("required_inequalities", {})
    id1_continuation = id1_quantitative.get("phi_seed_continuation_from_GR0", {})
    id1_boundary = id1_quantitative.get("epistemic_boundary", {})
    id1_constraint_gauge = id1_quantitative.get("constraint_and_gauge_contract", {})
    id1_true_payload = {
        "nonzero_width_family",
        "finite_Misner_Sharp_mass",
        "independent_chi_pulse",
        "explicit_nonzero_phi_seed",
        "physical_and_gauge_constraints_solved",
        "exact_outer_vacuum_buffer",
        "continuation_failure_reasons_typed",
    }
    rk_order = id1_convergence.get("RK4", {}).get(
        "observed_common_grid_profile_order"
    )
    ssp_order = id1_convergence.get("SSPRK3", {}).get(
        "observed_common_grid_profile_order"
    )
    maximum_residual = id1_extrema.get("maximum_constraint_residual_infinity")
    minimum_determinant = id1_extrema.get(
        "minimum_abs_constraint_jacobian_determinant"
    )
    maximum_condition = id1_extrema.get(
        "maximum_constraint_jacobian_condition_infinity"
    )
    minimum_compactness = id1_extrema.get("minimum_peak_compactness")
    maximum_compactness = id1_extrema.get("maximum_peak_compactness")
    if (
        id1.get("classification")
        != "finite_mass_constraint_compatible_data_family_certificate"
        or id1.get("generated_by") != "scripts/reproduce_fgc_id1_fam1.py"
        or any(id1_payload.get(key) is not True for key in id1_true_payload)
        or len(id1_exact) != 2
        or any(
            not isinstance(control, dict)
            or control.get("exact_constraint_pair_equality") is not True
            or control.get("exact_metric_defined_gauge_zero") is not True
            or control.get("source_is_unredefined_ACT1_VAR1") is not True
            for control in id1_exact
        )
        or id1_gr0.get("exact_pair_equality") is not True
        or not isinstance(rk_order, (int, float))
        or not isfinite(rk_order)
        or rk_order <= 3.0
        or not isinstance(ssp_order, (int, float))
        or not isfinite(ssp_order)
        or ssp_order <= 2.5
        or id1_family.get("exact_positive_parameter_volume") != "3/8388608"
        or id1_family.get("tensor_grid_member_count") != 27
        or len(id1_family.get("records", [])) != 27
        or not isinstance(id1_inequalities, dict)
        or not id1_inequalities
        or any(value is not True for value in id1_inequalities.values())
        or id1_family.get(
            "smooth_local_open_family_follows_from_regular_ODE_continuous_dependence"
        )
        is not True
        or id1_family.get("whole_declared_box_interval_enclosure_proven") is not False
        or id1_family.get("PROTO3_width_nine_eighths_case_qualified_here") is not False
        or not isinstance(maximum_residual, (int, float))
        or not isfinite(maximum_residual)
        or maximum_residual >= 1e-12
        or not isinstance(minimum_determinant, (int, float))
        or minimum_determinant <= 3.0
        or not isinstance(maximum_condition, (int, float))
        or maximum_condition >= 20.0
        or not isinstance(minimum_compactness, (int, float))
        or minimum_compactness < 0.1
        or not isinstance(maximum_compactness, (int, float))
        or maximum_compactness > 0.75
        or len(id1_continuation.get("records", [])) != 5
        or id1_continuation.get("zero_seed_constraints_reduce_exactly_to_GR0")
        is not True
        or id1_continuation.get(
            "all_nonzero_stages_retain_strict_constraint_branch_margins"
        )
        is not True
        or id1_constraint_gauge.get(
            "Hamiltonian_and_momentum_solved_at_every_reported_point"
        )
        is not True
        or id1_constraint_gauge.get(
            "metric_defined_Ct_and_Cr_solved_by_exact_algebraic_rates"
        )
        is not True
        or id1_boundary.get("initial_hypersurfaces_constructed") is not True
        or any(
            id1_boundary.get(key) is not False
            for key in (
                "time_evolution_performed",
                "collapse_or_trapped_surface_derived",
                "Raychaudhuri_defocusing_tested",
                "FGCQR_holdout_outcomes_inspected",
                "retained_EFT_validity_established",
            )
        )
        or id1.get("scope_bindings", {}).get("declared_run_envelope_sha256")
        != protocol.get("semantic_holdout_contract_sha256")
        or any(value is not False for value in id1.get("nonclaims", {}).values())
    ):
        failures.append("FGC-1-ID1-FAM1 family and epistemic boundary has drifted")
    return failures




def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--only-hlt14-mon14",
        action="store_true",
        help=(
            "verify the compact public repository surface plus the "
            "FGC-1-HLT14-MON14 synthetic-only PROTO17 adapter/runtime "
            "qualification; do not open a production namespace"
        ),
    )
    parser.add_argument(
        "--only-pro18-frz1",
        action="store_true",
        help=(
            "verify the compact public repository surface plus the "
            "FGC-1-PRO18-FRZ1 static production-prelaunch overlay; do not "
            "open historical campaign data or re-observe a future namespace"
        ),
    )
    parser.add_argument(
        "--only-pro18-pref26",
        action="store_true",
        help=(
            "verify the FGC-1-PRO18-PREF26 read-only raw-source binder; it "
            "opens only declared historical source files and never a future "
            "PROTO17 output root"
        ),
    )
    parser.add_argument(
        "--only-pro18-auth1",
        action="store_true",
        help=(
            "verify the FGC-1-PRO18-AUTH1 compact calibration authority "
            "content; do not authenticate a launch, open raw predecessors, "
            "or observe a future PROTO17 output root"
        ),
    )
    parser.add_argument(
        "--only-pro18-pref27",
        action="store_true",
        help=(
            "verify compact FGC-1-PRO18-PREF27 evidence and its historical "
            "optional generation-zero live-store shape; this explicit temporal "
            "audit is expected to reject a legitimately progressed successor "
            "store and does not replay GEN1, raw reimport, or refusal controls"
        ),
    )
    parser.add_argument("--only-pro19-frz1", action="store_true", help="verify compact PRO19 first-edge freeze only")
    parser.add_argument(
        "--only-hlt16-mon16",
        action="store_true",
        help=(
            "verify the sealed compact HLT16 durable-runtime qualification; "
            "do not open or advance the authenticated run store"
        ),
    )
    parser.add_argument(
        "--only-pro19-prelaunch",
        action="store_true",
        help=(
            "verify the sealed PRO19 committed-image manifest and compact "
            "prelaunch evidence; do not run the first event"
        ),
    )
    parser.add_argument(
        "--only-pro19-cfl1-prelaunch",
        action="store_true",
        help=(
            "verify the corrected committed image and exact generation-six "
            "continuation certificate without taking a writer lease"
        ),
    )
    parser.add_argument(
        "--only-pro19-sid1-prelaunch",
        action="store_true",
        help=(
            "verify the one-time FGC-1-PRO19-SID1-FRZ1 generation-seven "
            "suffix anchor before recovery; this live temporal audit must not "
            "be rerun after its authorized successor checkpoint exists"
        ),
    )
    parser.add_argument(
        "--only-pro19-sid2-postrecovery",
        action="store_true",
        help=(
            "verify the one-time FGC-1-PRO19-SID2-PREF1 generation-eight "
            "post-recovery binder and exact live-store identity; this "
            "explicit temporal audit does not authorize trajectory resume"
        ),
    )
    parser.add_argument(
        "--only-pro19-sid3-prelaunch",
        action="store_true",
        help=(
            "require and compactly verify the complete "
            "FGC-1-PRO19-SID3-AUTH1 committed-image authority bundle; do "
            "not open or mutate the authenticated generation-eight store"
        ),
    )
    parser.add_argument(
        "--only-pro19-pref28-postattempt",
        action="store_true",
        help=(
            "verify the FGC-1-PRO19-PREF28 compact certificate and exact "
            "read-only generation-nine terminal store; do not diagnose, "
            "repair, resume, or open a candidate branch"
        ),
    )
    parser.add_argument(
        "--only-tdg8-successor-prelaunch",
        action="store_true",
        help=(
            "require and compactly verify the FGC-1-TDG8-RUN1-AUTH1 "
            "successor authority; do not open either campaign, run the event, "
            "or re-observe destination absence"
        ),
    )
    parser.add_argument(
        "--only-tdg8-rcv1-prelaunch",
        action="store_true",
        help=(
            "require and compactly verify the FGC-1-TDG8-RCV1-AUTH1 "
            "recovery/resume authority; do not open or mutate the live run "
            "store"
        ),
    )
    parser.add_argument(
        "--only-tdg8-rcv2-prelaunch",
        action="store_true",
        help=(
            "require and compactly verify the FGC-1-TDG8-RCV2-AUTH1 "
            "bounded retry-chain authority; do not open or mutate the live "
            "run store"
        ),
    )
    parser.add_argument(
        "--only-tdg8-rcv3-prelaunch",
        action="store_true",
        help=(
            "require and compactly verify the FGC-1-TDG8-RCV3-FRZ1 "
            "generation-nine prefix projection freeze; do not open either "
            "run store, re-observe destination absence, or bootstrap"
        ),
    )
    parser.add_argument(
        "--only-tdg8-rcv3-pref1",
        action="store_true",
        help=(
            "require and compactly verify the FGC-1-TDG8-RCV3-PREF1 "
            "post-install binder; do not open the source or installed run "
            "store"
        ),
    )
    parser.add_argument(
        "--only-tdg8-rcv3-auth1",
        action="store_true",
        help=(
            "require and compactly verify the FGC-1-TDG8-RCV3-AUTH1 "
            "same-event execution authority; do not open or mutate the "
            "installed run store"
        ),
    )
    parser.add_argument(
        "--only-tdg8-rcv3-rec1-auth1",
        action="store_true",
        help=(
            "require and compactly verify the FGC-1-TDG8-RCV3-REC1-AUTH1 "
            "exact-recovery/same-event authority; do not open, reconcile, or "
            "mutate the installed run store"
        ),
    )
    parser.add_argument(
        "--only-tdg8-rcv3-pref2",
        action="store_true",
        help=(
            "require and compactly verify the FGC-1-TDG8-RCV3-PREF2 "
            "outcome-neutral terminal binder; do not open or mutate the "
            "installed run store"
        ),
    )
    parser.add_argument(
        "--only-tdg9-ar1-auth1",
        action="store_true",
        help=(
            "require and compactly verify the prospective "
            "FGC-1-TDG9-AR1-AUTH1 environment/selection/committed-image-bound "
            "exact-arithmetic authority; do not open "
            "the campaign or diagnostic output namespace, replay a proposal, "
            "or authorize continuation"
        ),
    )
    parser.add_argument(
        "--only-tdg9-ar1-pref1",
        action="store_true",
        help=(
            "require and compactly verify the independently bound "
            "FGC-1-TDG9-AR1-PREF1 terminal; do not open the campaign or raw "
            "diagnostic namespace or rerun an exact evaluator"
        ),
    )
    parser.add_argument(
        "--only-tdg9-loc1-frz1",
        action="store_true",
        help=(
            "require and compactly verify prospective FGC-1-TDG9-LOC1-FRZ1; "
            "do not open the campaign or diagnostic namespace or localize a cubic"
        ),
    )
    parser.add_argument(
        "--only-tdg9-loc2-frz1",
        action="store_true",
        help=(
            "require and compactly verify FGC-1-TDG9-LOC1-ERR1 and the "
            "prospective FGC-1-TDG9-LOC2-FRZ1 authority; do not open the "
            "campaign or diagnostic namespace or localize a cubic"
        ),
    )
    parser.add_argument(
        "--only-tdg9-loc2-pref2",
        action="store_true",
        help=(
            "require and compactly verify the independently bound "
            "FGC-1-TDG9-LOC2-PREF2 terminal; do not open the campaign or raw "
            "LOC2 namespace or rerun a localizer"
        ),
    )
    parser.add_argument(
        "--only-tdg9-ti1-frz1",
        action="store_true",
        help=(
            "require and compactly verify the historical invalid "
            "FGC-1-TDG9-TI1-FRZ1 authority; do not open "
            "the campaign store, raw namespace, or construct a shadow proposal"
        ),
    )
    parser.add_argument(
        "--only-tdg9-ti2-frz1",
        action="store_true",
        help=(
            "require and compactly verify the prospective fingerprint-only "
            "FGC-1-TDG9-TI2-FRZ1 recovery; do not open the campaign store, "
            "raw namespace, restore a live predecessor, or construct a shadow"
        ),
    )
    parser.add_argument(
        "--only-tdg9-ti2-pref1",
        action="store_true",
        help=(
            "require and compactly verify the independently bound "
            "FGC-1-TDG9-TI2-PREF1 terminal; do not open the campaign store, "
            "raw namespace, restore a predecessor, or construct a shadow"
        ),
    )
    parser.add_argument(
        "--only-tdg9-ac1-frz1",
        action="store_true",
        help=(
            "require and compactly verify the prospective retry-3 "
            "FGC-1-TDG9-AC1-FRZ1 all-channel audit; do not open the "
            "campaign store, raw namespace, restore a predecessor, or "
            "construct a shadow"
        ),
    )
    parser.add_argument(
        "--only-tdg9-ac1-pref1",
        action="store_true",
        help=(
            "require and compactly verify the independently bound "
            "FGC-1-TDG9-AC1-PREF1 terminal; do not open the campaign store, "
            "raw namespace, restore a predecessor, or construct a shadow"
        ),
    )
    parser.add_argument(
        "--only-tdg9-ur1-frz1",
        action="store_true",
        help=(
            "require and compactly verify the prospective retry-3 "
            "FGC-1-TDG9-UR1-FRZ1 u:R envelope-owner diagnostic; do not open "
            "the campaign store, raw namespace, restore a live predecessor, "
            "or construct a shadow"
        ),
    )
    parser.add_argument(
        "--only-tdg9-ur1-pref1",
        action="store_true",
        help=(
            "require and compactly verify the independently bound "
            "FGC-1-TDG9-UR1-PREF1 terminal; do not open the campaign store, "
            "raw namespace, restore a predecessor, or construct a shadow"
        ),
    )
    parser.add_argument(
        "--only-tdg10-qa1-frz1",
        action="store_true",
        help=(
            "require and compactly verify the prospective retry-3 "
            "FGC-1-TDG10-QA1-FRZ1 exact complete-C qualification; do not "
            "open the campaign store, raw namespace, restore a live "
            "predecessor, or construct a shadow"
        ),
    )
    parser.add_argument(
        "--only-tdg10-qa1-pref1",
        action="store_true",
        help=(
            "require and compactly verify the independently bound retry-3 "
            "FGC-1-TDG10-QA1-PREF1 terminal; do not open the campaign store, "
            "raw namespace, restore a predecessor, construct a shadow, or "
            "run the live binder"
        ),
    )
    parser.add_argument(
        "--only-tdg10-qa2-frz1",
        action="store_true",
        help=(
            "require and compactly verify the prospective retry-4/5 "
            "FGC-1-TDG10-QA2-FRZ1 two-width exact complete-C robustness "
            "qualification; do not open the campaign store, raw namespace, "
            "restore a live predecessor, or construct a shadow"
        ),
    )
    parser.add_argument(
        "--only-tdg10-qa2-pref1",
        action="store_true",
        help=(
            "require and compactly verify the independently bound retry-4/5 "
            "FGC-1-TDG10-QA2-PREF1 terminal; do not open the campaign store, "
            "raw namespace, restore a predecessor, construct a shadow, run "
            "the live binder, or launch QA2"
        ),
    )
    parser.add_argument(
        "--only-tdg11-msel1-frz1", action="store_true",
        help="Verify only the compact TDG11 selection freeze; no Git, raw store, shadow, or live environment access.",
    )
    parser.add_argument(
        "--only-tdg11-msel1-pref1",
        action="store_true",
        help=(
            "require and compactly verify the independently bound "
            "FGC-1-TDG11-MSEL1-PREF1 result through binder.verify_compact; "
            "do not open the campaign store, raw namespace, restore a "
            "predecessor, construct a shadow, or run live --bind --write-result"
        ),
    )
    parser.add_argument(
        "--only-tdg11-imp1",
        action="store_true",
        help="Verify only the compact IMP1 implementation certificate; no synthetic qualification or live observations.",
    )
    parser.add_argument(
        "--only-def1-stab1-frz1",
        action="store_true",
        help=(
            "require and compactly verify the candidate-blind "
            "FGC-1-DEF1-STAB1-FRZ1 instrument freeze; do not reconstruct "
            "from live owners, inspect Git, open runs/, or read a trajectory"
        ),
    )
    parser.add_argument(
        "--only-def1-stab1-pref1",
        action="store_true",
        help=(
            "require and compactly verify the independent "
            "FGC-1-DEF1-STAB1-PREF1 error-map readiness binder; do not "
            "reconstruct from live owners, inspect Git, open runs/, or "
            "read a trajectory"
        ),
    )
    parser.add_argument(
        "--only-sgb1-ctl1-sol1-frz1",
        action="store_true",
        help=(
            "require and compactly verify the prospective SOL1-FRZ1 scoped "
            "continuation freeze; do not reconstruct, inspect Git, or open runs/"
        ),
    )
    parser.add_argument(
        "--only-sgb1-ctl1-sol1-pref1",
        action="store_true",
        help=(
            "require and compactly verify the independently bound SOL1-PREF1 "
            "scoped outcome; do not reconstruct, inspect Git, or open runs/"
        ),
    )
    parser.add_argument(
        "--only-hlt13-mon13",
        action="store_true",
        help=(
            "verify the compact public repository surface plus the "
            "FGC-1-HLT13-MON13 certificate and synthetic runtime binding; "
            "do not traverse historical campaign data"
        ),
    )
    parser.add_argument(
        "--only-pro16-frz1",
        action="store_true",
        help=(
            "verify the compact public repository surface plus the "
            "FGC-1-PRO16-FRZ1 transition/genesis freeze; do not traverse "
            "historical campaign data"
        ),
    )
    parser.add_argument(
        "--only-pro16-pref24",
        action="store_true",
        help=(
            "verify the compact public repository surface plus the "
            "FGC-1-PRO16-PREF24 acyclicity/under-specification binder; "
            "do not traverse historical campaign data"
        ),
    )
    parser.add_argument(
        "--only-pro17-frz1",
        action="store_true",
        help=(
            "verify the compact public repository surface plus the "
            "FGC-1-PRO17-FRZ1 pure construction freeze; do not traverse "
            "historical campaign data or a production namespace"
        ),
    )
    parser.add_argument(
        "--only-pro17-pref25",
        action="store_true",
        help=(
            "verify the compact public repository surface plus the "
            "FGC-1-PRO17-PREF25 independent construction theorem; do not "
            "open raw campaign data or a production namespace"
        ),
    )
    args = parser.parse_args(argv)
    if args.only_def1_stab1_frz1 and any(
        value for name, value in vars(args).items()
        if name.startswith("only_") and name != "only_def1_stab1_frz1"
    ):
        parser.error("only one focused compact certificate may be selected")
    if args.only_def1_stab1_pref1 and any(
        value for name, value in vars(args).items()
        if name.startswith("only_") and name != "only_def1_stab1_pref1"
    ):
        parser.error("only one focused compact certificate may be selected")
    if args.only_sgb1_ctl1_sol1_frz1 and any(
        value for name, value in vars(args).items()
        if name.startswith("only_") and name != "only_sgb1_ctl1_sol1_frz1"
    ):
        parser.error("only one focused compact certificate may be selected")
    if args.only_sgb1_ctl1_sol1_pref1 and any(
        value for name, value in vars(args).items()
        if name.startswith("only_") and name != "only_sgb1_ctl1_sol1_pref1"
    ):
        parser.error("only one focused compact certificate may be selected")
    if args.only_tdg11_imp1 and any(
        value for name, value in vars(args).items()
        if name.startswith("only_") and name != "only_tdg11_imp1"
    ):
        parser.error("only one focused compact certificate may be selected")
    if args.only_tdg11_msel1_frz1 and any(value for name, value in vars(args).items() if name.startswith("only_") and name != "only_tdg11_msel1_frz1"):
        parser.error("only one focused compact certificate may be selected")
    if args.only_tdg11_msel1_pref1 and any(value for name, value in vars(args).items() if name.startswith("only_") and name != "only_tdg11_msel1_pref1"):
        parser.error("only one focused compact certificate may be selected")
    if sum((args.only_hlt14_mon14, args.only_hlt13_mon13, args.only_pro16_frz1, args.only_pro16_pref24, args.only_pro17_frz1, args.only_pro17_pref25, args.only_pro18_frz1, args.only_pro18_pref26, args.only_pro18_auth1, args.only_pro18_pref27, args.only_pro19_frz1, args.only_hlt16_mon16, args.only_pro19_prelaunch, args.only_pro19_cfl1_prelaunch, args.only_pro19_sid1_prelaunch, args.only_pro19_sid2_postrecovery, args.only_pro19_sid3_prelaunch, args.only_pro19_pref28_postattempt, args.only_tdg8_successor_prelaunch, args.only_tdg8_rcv1_prelaunch, args.only_tdg8_rcv2_prelaunch, args.only_tdg8_rcv3_prelaunch, args.only_tdg8_rcv3_pref1, args.only_tdg8_rcv3_auth1, args.only_tdg8_rcv3_rec1_auth1, args.only_tdg8_rcv3_pref2, args.only_tdg9_ar1_auth1, args.only_tdg9_ar1_pref1, args.only_tdg9_loc1_frz1, args.only_tdg9_loc2_frz1, args.only_tdg9_loc2_pref2, args.only_tdg9_ti1_frz1, args.only_tdg9_ti2_frz1, args.only_tdg9_ti2_pref1, args.only_tdg9_ac1_frz1, args.only_tdg9_ac1_pref1, args.only_tdg9_ur1_frz1, args.only_tdg9_ur1_pref1, args.only_tdg10_qa1_frz1, args.only_tdg10_qa1_pref1, args.only_tdg10_qa2_frz1, args.only_tdg10_qa2_pref1, args.only_tdg11_msel1_pref1)) > 1:
        parser.error("only one focused compact certificate may be selected")
    if args.only_sgb1_ctl1_sol1_pref1:
        checks = ((
            "FGC-1-SGB1-CTL1-SOL1-PREF1 compact binder",
            _check_sgb1_ctl1_sol1_pref1_result,
        ),)
    elif args.only_sgb1_ctl1_sol1_frz1:
        checks = ((
            "FGC-1-SGB1-CTL1-SOL1-FRZ1 compact freeze",
            _check_sgb1_ctl1_sol1_frz1_result,
        ),)
    elif args.only_def1_stab1_frz1:
        checks = (
            (
                "FGC-1-DEF1-STAB1-FRZ1 compact freeze",
                _check_def1_stab1_frz1_result,
            ),
        )
    elif args.only_def1_stab1_pref1:
        checks = (
            (
                "FGC-1-DEF1-STAB1-PREF1 compact binder",
                _check_def1_stab1_pref1_result,
            ),
        )
    elif args.only_tdg11_imp1:
        checks = (
            ("required files", check_required),
            ("TDG11 compact qualified implementation", _check_tdg11_imp1_result),
        )
    elif args.only_tdg11_msel1_frz1:
        checks = (("TDG11 compact three-route selection freeze", _check_tdg11_msel1_frz1_result),)
    elif args.only_tdg11_msel1_pref1:
        checks = (
            ("required files", check_required),
            (
                "TDG11 compact independent selection binder",
                _check_tdg11_msel1_pref1_result,
            ),
        )
    elif args.only_tdg9_ar1_auth1:
        checks = (
            ("required files", check_required),
            (
                "TDG9 compact exact-arithmetic authority",
                _check_tdg9_ar1_auth1_result,
            ),
        )
    elif args.only_tdg9_ar1_pref1:
        checks = (
            ("required files", check_required),
            (
                "TDG9 compact exact-arithmetic terminal binder",
                _check_tdg9_ar1_pref1_result,
            ),
        )
    elif args.only_tdg9_loc1_frz1:
        checks = (
            ("required files", check_required),
            ("TDG9 compact exact localization authority", _check_tdg9_loc1_frz1_result),
        )
    elif args.only_tdg9_loc2_frz1:
        checks = (
            ("required files", check_required),
            ("TDG9 compact LOC1 diagnosis and LOC2 authority", _check_tdg9_loc2_frz1_result),
        )
    elif args.only_tdg9_loc2_pref2:
        checks = (
            ("required files", check_required),
            (
                "TDG9 compact LOC2 terminal binder",
                _check_tdg9_loc2_pref2_result,
            ),
        )
    elif args.only_tdg9_ti1_frz1:
        checks = (
            ("required files", check_required),
            (
                "TDG9 compact same-state tableau authority",
                _check_tdg9_ti1_frz1_result,
            ),
        )
    elif args.only_tdg9_ti2_frz1:
        checks = (
            ("required files", check_required),
            (
                "TDG9 compact fingerprint recovery authority",
                _check_tdg9_ti2_frz1_result,
            ),
        )
    elif args.only_tdg9_ti2_pref1:
        checks = (
            ("required files", check_required),
            (
                "TDG9 compact tableau-shadow terminal binder",
                _check_tdg9_ti2_pref1_result,
            ),
        )
    elif args.only_tdg9_ac1_frz1:
        checks = (
            ("required files", check_required),
            (
                "TDG9 compact retry-3 all-channel production-classifier audit",
                _check_tdg9_ac1_frz1_result,
            ),
        )
    elif args.only_tdg9_ac1_pref1:
        checks = (
            ("required files", check_required),
            (
                "TDG9 compact retry-3 all-channel production-classifier binder",
                _check_tdg9_ac1_pref1_result,
            ),
        )
    elif args.only_tdg9_ur1_frz1:
        checks = (
            ("required files", check_required),
            (
                "TDG9 compact retry-3 u:R envelope-owner diagnostic",
                _check_tdg9_ur1_frz1_result,
            ),
        )
    elif args.only_tdg9_ur1_pref1:
        checks = (
            ("required files", check_required),
            (
                "TDG9 compact retry-3 u:R envelope-owner binder",
                _check_tdg9_ur1_pref1_result,
            ),
        )
    elif args.only_tdg10_qa1_frz1:
        checks = (
            ("required files", check_required),
            (
                "TDG10 compact retry-3 exact complete-C qualification",
                _check_tdg10_qa1_frz1_result,
            ),
        )
    elif args.only_tdg10_qa1_pref1:
        checks = (
            ("required files", check_required),
            (
                "TDG10 compact retry-3 exact complete-C terminal binder",
                _check_tdg10_qa1_pref1_result,
            ),
        )
    elif args.only_tdg10_qa2_frz1:
        checks = (
            ("required files", check_required),
            (
                "TDG10 compact retry-4/5 exact complete-C two-width robustness",
                _check_tdg10_qa2_frz1_result,
            ),
        )
    elif args.only_tdg10_qa2_pref1:
        checks = (
            ("required files", check_required),
            (
                "TDG10 compact retry-4/5 exact complete-C terminal binder",
                _check_tdg10_qa2_pref1_result,
            ),
        )
    elif args.only_tdg8_rcv3_pref2:
        checks = (
            ("required files", check_required),
            (
                "TDG8 compact RCV3 terminal binder",
                _check_tdg8_rcv3_pref2_result,
            ),
        )
    elif args.only_tdg8_rcv3_rec1_auth1:
        checks = (
            ("required files", check_required),
            (
                "TDG8 compact RCV3 exact-recovery authority",
                _check_tdg8_rcv3_rec1_auth1_result,
            ),
        )
    elif args.only_tdg8_rcv3_auth1:
        checks = (
            ("required files", check_required),
            (
                "TDG8 compact RCV3 same-event authority",
                _check_tdg8_rcv3_auth1_result,
            ),
        )
    elif args.only_tdg8_rcv3_pref1:
        checks = (
            ("required files", check_required),
            (
                "TDG8 compact installed projection binder",
                _check_tdg8_rcv3_pref1_result,
            ),
        )
    elif args.only_tdg8_rcv3_prelaunch:
        checks = (
            ("required files", check_required),
            (
                "TDG8 compact generation-nine projection freeze",
                _check_tdg8_rcv3_prelaunch_result,
            ),
        )
    elif args.only_tdg8_rcv2_prelaunch:
        checks = (
            ("required files", check_required),
            (
                "TDG8 compact bounded retry-chain authority",
                _check_tdg8_rcv2_prelaunch_result,
            ),
        )
    elif args.only_tdg8_rcv1_prelaunch:
        checks = (
            ("required files", check_required),
            (
                "TDG8 compact retry recovery/resume authority",
                _check_tdg8_rcv1_prelaunch_result,
            ),
        )
    elif args.only_tdg8_successor_prelaunch:
        checks = (
            ("required files", check_required),
            (
                "TDG8 compact fresh-successor authority",
                _check_tdg8_successor_prelaunch_result,
            ),
        )
    elif args.only_pro19_pref28_postattempt:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            (
                "PRO19 exact PREF28 post-attempt binder",
                _check_pro19_pref28_postattempt_result,
            ),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro19_sid3_prelaunch:
        checks = (
            ("required files", check_required),
            (
                "PRO19 compact SID3 committed-image authority",
                lambda: _check_pro19_sid3_compact_result(required=True),
            ),
            ("Git/data boundary", check_git_data_boundary),
        )
    elif args.only_pro19_sid2_postrecovery:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            (
                "PRO19 exact SID2 post-recovery binder",
                _check_pro19_sid2_postrecovery_result,
            ),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro19_sid1_prelaunch:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("PRO19 exact SID1 pre-recovery authority", _check_pro19_sid1_prelaunch_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro19_cfl1_prelaunch:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("PRO19 exact CFL1 continuation", _check_pro19_cfl1_prelaunch_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro19_prelaunch:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("PRO19 sealed prelaunch evidence", _check_pro19_prelaunch_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_hlt16_mon16:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("FGC-1-HLT16-MON16 compact qualification", _check_hlt16_mon16_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro19_frz1:
        checks = (("required files", check_required), ("FGC-1-PRO19-FRZ1 compact freeze", _check_pro19_frz1_result))
    elif args.only_pro18_pref27:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            (
                "FGC-1-PRO18-PREF27 compact and historical live evidence",
                lambda: _check_pro18_pref27_result(inspect_live_store=True),
            ),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro18_auth1:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("FGC-1-PRO18-AUTH1 compact content authority", _check_pro18_auth1_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro18_pref26:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("FGC-1-PRO18-PREF26 raw-dependent certificate", _check_pro18_pref26_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro18_frz1:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("FGC-1-PRO18-FRZ1 compact certificate", _check_pro18_frz1_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_hlt14_mon14:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("FGC-1-HLT14-MON14 compact certificate", _check_hlt14_mon14_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_hlt13_mon13:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("FGC-1-HLT13-MON13 compact certificate", _check_hlt13_mon13_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro16_frz1:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("FGC-1-PRO16-FRZ1 compact certificate", _check_pro16_frz1_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro16_pref24:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("FGC-1-PRO16-PREF24 compact certificate", _check_pro16_pref24_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro17_frz1:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("FGC-1-PRO17-FRZ1 compact certificate", _check_pro17_frz1_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    elif args.only_pro17_pref25:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("FGC-1-PRO17-PREF25 compact certificate", _check_pro17_pref25_result),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    else:
        checks = (
            ("required files", check_required),
            ("local Markdown links", check_links),
            ("public text hygiene", check_public_text),
            ("FGC, NGS, DSF, and epistemic alignment", check_fgc_alignment),
            ("machine-readable result", check_results),
            ("Git/data boundary", check_git_data_boundary),
            ("reference sequence", check_reference_sequence),
        )
    # Focused compact artifact modes preserve their historical I/O boundary.
    # Catalog regeneration inspects Git and therefore belongs to the full
    # repository audit (and the explicit current-development route), not to
    # a raw/store/shadow/Git-blind compact binder check.
    if not any(
        value for name, value in vars(args).items() if name.startswith("only_")
    ):
        checks = (*checks, ("artifact catalog and result index", check_artifact_catalog))
        if (REPOSITORY / "configs/fgc/fgc-1-tdg11-msel1-frz1.toml").exists():
            checks = (*checks, ("TDG11 compact three-route selection freeze", _check_tdg11_msel1_frz1_result))
    failures: list[str] = []
    for label, check in checks:
        findings = check()
        if findings:
            print(f"FAIL {label}")
            failures.extend(findings)
            for finding in findings:
                print(f"  - {finding}")
        else:
            print(f"PASS {label}")
    if failures:
        print(f"repository check failed with {len(failures)} finding(s)", file=sys.stderr)
        return 1
    print("repository check complete: canonical public surface is internally consistent")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

__all__ = tuple(name for name in globals() if not name.startswith("__"))
