"""Fail-closed retained-EFT run-authorization audit for FGC-1.

The existing EFT0 ledger proves a few exact component-amplitude inequalities
under an externally supplied cutoff.  Radial hyperbolicity and frozen modal
boundary certificates add genuine PDE evidence, but none of those facts bounds
proper frequency, completes the Wilsonian operator basis, controls the total
omitted remainder, or proves a constraint-complete quasilinear problem.

This module composes those records into one explicit conjunction.  It does not
infer missing evidence and it never treats a local component box as an
evolution-invariant EFT domain.
"""

from __future__ import annotations

from typing import Any, Mapping


EFT1_OPEN1_ARTIFACT_ID = "FGC-1-EFT1-OPEN1"


def _mapping(name: str, value: object) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return value


def _artifact(name: str, value: Mapping[str, Any], expected_id: str) -> None:
    if value.get("artifact_id") != expected_id:
        raise ValueError(f"{name} artifact identity differs")
    nonclaims = _mapping(f"{name}.nonclaims", value.get("nonclaims"))
    if any(item is not False for item in nonclaims.values()):
        raise ValueError(f"{name} contains a promoted nonclaim")


def _gate(name: str, value: Mapping[str, Any], key: str, expected: bool) -> None:
    status = _mapping(f"{name}.gate_status", value.get("gate_status"))
    if status.get(key) is not expected:
        raise ValueError(f"{name} gate {key} must be {str(expected).lower()}")


def _require_boolean(name: str, value: object) -> bool:
    if not isinstance(value, bool):
        raise ValueError(f"{name} must be boolean")
    return value


def retained_eft_open_run_audit(
    eft0: Mapping[str, Any],
    uhyp1: Mapping[str, Any],
    bnd1: Mapping[str, Any],
    con3: Mapping[str, Any],
) -> dict[str, Any]:
    """Compose the current exact records into a necessary-condition audit.

    The authorization predicate is an all-of conjunction.  Passing evidence is
    preserved as partial progress, while every unavailable premise is named
    and mapped to a deterministic stop decision.
    """

    eft0 = _mapping("EFT0-LED1", eft0)
    uhyp1 = _mapping("DOM3-UHYP1", uhyp1)
    bnd1 = _mapping("BND1-MD1", bnd1)
    con3 = _mapping("CON3-CAU1", con3)
    _artifact("EFT0-LED1", eft0, "FGC-1-EFT0-LED1")
    _artifact("DOM3-UHYP1", uhyp1, "FGC-1-HYP1-DOM3-UHYP1")
    _artifact("BND1-MD1", bnd1, "FGC-1-HYP1-BND1-MD1")
    _artifact("CON3-CAU1", con3, "FGC-1-HYP1-CON3-CAU1")

    _gate(
        "EFT0-LED1",
        eft0,
        "conditional_declared_local_component_controls_passed",
        True,
    )
    _gate(
        "EFT0-LED1",
        eft0,
        "retained_eft_validity_envelope_passed",
        False,
    )
    _gate(
        "DOM3-UHYP1",
        uhyp1,
        "compact_nonflat_radial_strong_hyperbolicity_on_implicit_branch_graph_passed",
        True,
    )
    _gate("DOM3-UHYP1", uhyp1, "multidirectional_strong_hyperbolicity_proven", False)
    _gate(
        "BND1-MD1",
        bnd1,
        "uniform_frozen_radial_main_system_boundary_dissipation_passed",
        True,
    )
    _gate("BND1-MD1", bnd1, "constraint_preserving_ACT1_IBVP_proven", False)
    _gate(
        "CON3-CAU1",
        con3,
        "conditional_boundary_free_gauge_Cauchy_uniqueness_theorem_derived",
        True,
    )
    _gate("CON3-CAU1", con3, "constraint_preserving_ACT1_IBVP_proven", False)

    ledger = _mapping(
        "EFT0-LED1.conditional_local_ledger_certificate",
        eft0.get("conditional_local_ledger_certificate"),
    )
    cutoff = _mapping("EFT0-LED1.cutoff", ledger.get("cutoff"))
    if cutoff.get("external_assumption") is not True:
        raise ValueError("EFT0-LED1 cutoff must remain an explicit external assumption")
    if cutoff.get("inferred_from_fixture") is not False:
        raise ValueError("EFT0-LED1 cutoff must not be inferred from fixture values")
    controls = _mapping("EFT0-LED1.local_controls", ledger.get("local_controls"))
    if not controls or any(
        _mapping(f"EFT0-LED1.local_controls.{name}", control).get("passes") is not True
        for name, control in controls.items()
    ):
        raise ValueError("EFT0-LED1 local component controls are not all passed")
    basis = _mapping("EFT0-LED1.operator_basis", ledger.get("operator_basis"))
    frequency = _mapping("EFT0-LED1.frequency_control", ledger.get("frequency_control"))
    component_box = _mapping("EFT0-LED1.component_box", ledger.get("component_box"))

    complete_basis = _require_boolean(
        "complete_wilsonian_basis", basis.get("complete_wilsonian_basis")
    )
    omitted_remainders = _require_boolean(
        "representative_omitted_remainders_locally_evaluated",
        basis.get("representative_omitted_remainders_locally_evaluated"),
    )
    frequency_available = _require_boolean("frequency.available", frequency.get("available"))
    full_jet_box = _require_boolean("component_box.full_jet_box", component_box.get("full_jet_box"))

    predicates = (
        {
            "id": "declared_cutoff_and_matching_scale",
            "passed": True,
            "evidence": "EFT0-LED1 records Lambda=16 as an external assumption",
            "failure_consequence": "no dimensionless EFT hierarchy can be evaluated",
        },
        {
            "id": "cutoff_not_inferred_from_fixture",
            "passed": True,
            "evidence": "EFT0-LED1 explicitly separates Lambda from code-unit fixture couplings",
            "failure_consequence": "the EFT scale would be circularly assigned",
        },
        {
            "id": "declared_local_component_controls",
            "passed": True,
            "evidence": "all exact EFT0-LED1 F and scalar-amplitude inequalities pass",
            "failure_consequence": "the retained local action leaves its declared component bounds",
        },
        {
            "id": "complete_symmetry_reduced_basis_through_target_order",
            "passed": complete_basis,
            "evidence": "EFT0-LED1 labels its omitted set representative_only",
            "failure_consequence": "unlisted allowed operators can dominate the claimed truncation error",
        },
        {
            "id": "coefficient_bounds_or_matching_for_complete_omitted_basis",
            "passed": False,
            "evidence": "only representative order-one Wilson bounds are assumed",
            "failure_consequence": "no bound exists on the total omitted contribution",
        },
        {
            "id": "invariant_proper_frequency_and_wavenumber_bound",
            "passed": frequency_available,
            "evidence": frequency.get("reason"),
            "failure_consequence": "short-wavelength data can violate the cutoff inside a small component box",
        },
        {
            "id": "covariant_field_gradient_and_curvature_bounds",
            "passed": full_jet_box,
            "evidence": "EFT0-LED1 explicitly records full_jet_box=false",
            "failure_consequence": "derivative and curvature expansion parameters are uncontrolled",
        },
        {
            "id": "uniform_total_omitted_remainder_bound",
            "passed": omitted_remainders,
            "evidence": basis.get("representative_omitted_remainders_reason"),
            "failure_consequence": "the retained equations have no quantified truncation error",
        },
        {
            "id": "multidirectional_full_system_strong_hyperbolicity",
            "passed": False,
            "evidence": "DOM3-UHYP1 proves only the compact radial branch graph",
            "failure_consequence": "a radial mode theorem cannot authorize a general Cauchy evolution",
        },
        {
            "id": "constraint_complete_local_existence_or_IBVP",
            "passed": False,
            "evidence": "BND1-MD1 is frozen main-system dissipation and CON3-CAU1 is boundary-free conditional gauge uniqueness",
            "failure_consequence": "there is no proved physical/gauge constraint-preserving evolution problem",
        },
        {
            "id": "nonzero_spacetime_run_domain_with_regular_center_or_declared_boundaries",
            "passed": False,
            "evidence": "the current compatible solution witness is one local jet at r=4 on a reference annulus",
            "failure_consequence": "no collapse initial slice or finite run region exists",
        },
        {
            "id": "evolution_bootstrap_monitor_and_fail_closed_stop_rule",
            "passed": False,
            "evidence": "no nonlinear evolution has been constructed on which to propagate the bounds",
            "failure_consequence": "even initially valid inequalities are not known to remain valid",
        },
    )
    missing = tuple(item["id"] for item in predicates if item["passed"] is not True)
    authorization = all(item["passed"] is True for item in predicates)
    if authorization or not missing:
        raise ValueError("current evidence was unexpectedly promoted to an open-run envelope")

    return {
        "artifact_id": EFT1_OPEN1_ARTIFACT_ID,
        "classification": "fail_closed_retained_EFT_open_run_authorization_audit_not_an_EFT_validity_or_evolution_certificate",
        "authorization_logic": "all_required_predicates_must_pass; one_false_predicate_stops_the_run",
        "partial_positive_evidence": {
            "externally_declared_cutoff": True,
            "exact_local_component_controls": True,
            "compact_nonflat_radial_strong_hyperbolicity": True,
            "uniform_frozen_radial_main_system_maximal_dissipation": True,
            "conditional_boundary_free_gauge_Cauchy_uniqueness": True,
        },
        "required_predicates": list(predicates),
        "passed_predicate_count": len(predicates) - len(missing),
        "required_predicate_count": len(predicates),
        "missing_predicate_ids": list(missing),
        "retained_EFT_open_run_envelope_passed": False,
        "evolution_authorized": False,
        "stop_mask": {
            "construct_or_tune_collapse_data_as_an_FGC_result": False,
            "run_nonlinear_FGC_collapse_evolution": False,
            "promote_point_Raychaudhuri_diagnostic_to_DEF1": False,
            "claim_finite_transition_surface": False,
            "claim_singularity_resolution": False,
        },
        "minimum_next_evidence": [
            "complete symmetry-reduced operator basis through a declared target order with matched or bounded coefficients",
            "covariant field-gradient, curvature, and proper-frequency bounds on a nonzero spacetime run domain",
            "uniform total omitted-remainder estimate below a predeclared tolerance",
            "multidirectional full-system hyperbolicity plus physical/gauge constraint-complete local existence or IBVP",
            "regular-centre or declared-boundary data and a bootstrap monitor that stops before any bound is crossed",
        ],
        "nonclaims": {
            "complete_Wilsonian_EFT_defined": False,
            "UV_completion_derived": False,
            "radiative_stability_proven": False,
            "retained_EFT_open_run_envelope_proven": False,
            "constraint_complete_IBVP_proven": False,
            "initial_data_constructed": False,
            "evolution_authorized": False,
            "collapse_solution_derived": False,
            "affine_null_defocusing_derived": False,
            "finite_invariant_transition_surface_derived": False,
            "singularity_resolution_derived": False,
        },
    }
