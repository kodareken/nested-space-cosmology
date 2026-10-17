#!/usr/bin/env python3
"""Authenticate prediction declarations and saved-value projections; no replay.

An optional NSC-P001 finite-amplitude extension is checked from the saved
comparison and its three endpoint receipts. Products use a relative tolerance.
Subtractions of radii also allow a few ulps of those operands, because a
one-ulp hypot change can dominate a tiny difference and a zero difference has
no relative scale. Identifiers, booleans, paths, hashes, saved residuals,
alpha, accuracy constants, and copied fields stay exact. Convergence booleans
are not relaxed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re

from analyze_nsc_atlas_fixed_endpoint import compare
from recursive_horizons.nsc_atlas_endpoint_comparison import (
    ENDPOINT_REFERENCE_RADIUS, ENDPOINT_SOURCE_X, FIDUCIAL_LENGTH_METRES,
)

ROOT = Path(__file__).resolve().parents[1]
REGISTER = "results/prediction-register.json"
DOCUMENT = "docs/prediction-register.md"
SCHEMA = "nsc-prediction-register-v1"
FINITE_SCHEMA = "nsc-p001-finite-amplitude-v1"
# Products and quotients may differ by this relative band across libm builds.
# A difference of two radii instead takes an absolute floor of a few ulps of
# its own operands: one ulp of hypot near radius 0.55 is enormous next to a
# tiny or zero paired gap, while remaining negligible beside a real gap.
RECOMPUTE_REL_TOL = 1e-12
RADIUS_DIFFERENCE_ULPS = 4
ENVELOPE_STATUS = "finite_amplitude_numerical_estimate_with_resolution_envelope"
UNCONVERGED_STATUS = "unconverged_numerical_estimate"
COMPARISON_META = ("analysis_source_commit", "analysis_source_sha256", "inputs")
EXACT_COMPARISON_KEYS = (
    "all_endpoints_completed", "endpoint_radius_pair_difference_contracts",
    "independent_observational_test", "numerical_method_note", "operator_revision_note",
    "physical_evolutions", "prediction_status", "radial_accuracy",
    "requested_mass_and_flux_contraction", "sample_cuts", "schema", "shared_alpha",
    "source_commits", "status",
)
OBSERVABLE_DERIVED = (
    "circumference_metres", "difference_from_quiet_reference_metres", "grid_values_metres",
    "quiet_reference_metres", "resolution_half_width_metres", "resolution_range_metres",
    "value_metres",
)
REFINEMENT_DERIVED = ("fine_over_coarse",)
FINITE_FIELDS = (
    "comparison_status", "endpoint_radius_pair_difference_contracts", "event", "evidence",
    "limits", "physical_observable", "prediction_status", "radial_accuracy",
    "requested_mass_and_flux_contraction", "schema", "shared_alpha",
)
LIMITS = dict(
    empirical_confirmation=False, calibration="none", physical_source_selected=False,
    physical_emitter_selected=False, observational_validation=False,
    rigorous_error_enclosure=False, historical_coefficients_predict_this_event=False,
    SI_scale_role="declared_fiducial_input", measured_astrophysical_scale=False,
)
EVIDENCE = {
    "cap_wall_owner": "docs/nsc-atlas-cap-wall-response.md",
    "wall_owner": "docs/nsc-atlas-wall-evolution.md",
    "model": "results/development/nsc-atlas-cap-wall-response-v1/model-specification.json",
    "cap_wall": "results/development/nsc-atlas-cap-wall-response-v1/cap-wall.json",
    "cap_wall_manifest": "results/development/nsc-atlas-cap-wall-response-v1/manifest.json",
    "wall_field": "results/development/nsc-atlas-wall-evolution-v1/field/response.json",
    "wall_field_manifest": "results/development/nsc-atlas-wall-evolution-v1/field/manifest.json",
}
COEFFICIENTS = ("mu_plus", "nu_minus", "rho", "eta")
IDENTITY_POINTERS = {
    "delta_mean_K": "/exact_algebra/zero_residuals/linear_mean_curvature_variation_vanishes",
    "delta_a_plus_minus_delta_a_minus": "/exact_algebra/zero_residuals/temporal_jump_from_angular_and_mixed_flux",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def local_file(root, relative):
    require(isinstance(relative, str) and "\\" not in relative and ":" not in relative,
            "evidence path must be repository-relative POSIX text")
    parts = PurePosixPath(relative)
    require(bool(parts.parts) and not parts.is_absolute() and ".." not in parts.parts,
            "evidence path escapes repository")
    path = root.joinpath(*parts.parts)
    require(not any(parent.is_symlink() for parent in (path, *path.parents) if parent != root.parent),
            "symlink evidence is not permitted")
    require(path.resolve().is_relative_to(root.resolve()) and path.is_file(),
            "missing or escaping evidence: " + relative)
    return path


def read_json(raw):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, "duplicate JSON key: " + key)
            result[key] = value
        return result
    def nonfinite(value):
        raise ValueError("nonfinite JSON value: " + value)
    return json.loads(raw, object_pairs_hook=unique, parse_constant=nonfinite)


def exact(actual, expected, label):
    # JSON equality distinguishes booleans from numeric coefficients.
    require(json.dumps(actual, sort_keys=True, allow_nan=False) ==
            json.dumps(expected, sort_keys=True, allow_nan=False), label + " differs")


def number(value, label):
    require(type(value) in (int, float) and math.isfinite(value), label + " must be finite real")
    return value


def pointer(value, path):
    require(isinstance(path, str) and path.startswith("/"), "invalid JSON pointer")
    for part in path[1:].split("/"):
        key = part.replace("~1", "/").replace("~0", "~")
        value = value[int(key)] if isinstance(value, list) else value[key]
    return value


def bound_bytes(root, row, label):
    exact(sorted(row), ["bytes", "path", "sha256"], label + " binding fields")
    require(type(row["bytes"]) is int and row["bytes"] > 0
            and isinstance(row["sha256"], str)
            and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is not None,
            "invalid byte identity: " + label)
    raw = local_file(root, row["path"]).read_bytes()
    require(len(raw) == row["bytes"] and hashlib.sha256(raw).hexdigest() == row["sha256"],
            "evidence bytes changed: " + label)
    return raw


def authenticate(root, bindings):
    exact(set_as_list(bindings), sorted(EVIDENCE), "evidence inventory")
    records = {}
    for name, expected_path in EVIDENCE.items():
        row = bindings[name]
        exact(sorted(row), ["bytes", "path", "sha256"], name + " binding fields")
        exact(row["path"], expected_path, name + " path")
        raw = bound_bytes(root, row, name)
        if row["path"].endswith(".json"):
            records[name] = read_json(raw)
    for manifest, output in (("cap_wall_manifest", "cap_wall"),
                             ("cap_wall_manifest", "model"),
                             ("wall_field_manifest", "wall_field")):
        row = bindings[output]
        exact(records[manifest]["outputs"][Path(row["path"]).name],
              {key:row[key] for key in ("bytes", "sha256")}, output + " original manifest")
    for manifest, owner in (("cap_wall_manifest", "cap_wall_owner"),
                            ("wall_field_manifest", "wall_owner")):
        exact(records[manifest]["source_bindings"][bindings[owner]["path"]],
              bindings[owner]["sha256"], owner + " original source binding")
    return records


def set_as_list(value):
    require(isinstance(value, dict), "expected named fields")
    return sorted(value)


def subtraction_abs_tol(operands):
    """Absolute slack for one subtraction, measured from that subtraction's operands."""
    require(len(operands) >= 2, "difference tolerance needs the operands")
    total = 0.0
    for value in operands:
        number(value, "difference operand")
        total += math.ulp(float(value))
    return RADIUS_DIFFERENCE_ULPS * total


def agree_float(actual, expected, label, abs_tol):
    require(type(actual) in (int, float) and type(expected) in (int, float)
            and not isinstance(actual, bool) and not isinstance(expected, bool)
            and math.isfinite(float(actual)) and math.isfinite(float(expected)),
            label + " must be a finite float")
    require(abs_tol >= 0 and math.isfinite(abs_tol), label + " tolerance")
    require(math.isclose(float(actual), float(expected), rel_tol=RECOMPUTE_REL_TOL, abs_tol=abs_tol),
            label + " exceeds recompute tolerance")


def agree_numbers(actual, expected, label, abs_tol=0.0):
    """Match structure exactly. Floats use relative slack plus an optional absolute floor."""
    if isinstance(expected, bool) or isinstance(actual, bool):
        exact(actual, expected, label)
    elif isinstance(expected, str) or isinstance(expected, int) or expected is None:
        exact(actual, expected, label)
    elif isinstance(expected, float):
        agree_float(actual, expected, label, abs_tol)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), label + " length")
        for index, item in enumerate(expected):
            agree_numbers(actual[index], item, label + "[" + str(index) + "]", abs_tol)
    elif isinstance(expected, dict):
        require(isinstance(actual, dict), label + " type")
        exact(sorted(actual), sorted(expected), label + " keys")
        for key, item in expected.items():
            agree_numbers(actual[key], item, label + "." + key, abs_tol)
    else:
        raise ValueError("unsupported recomputed value: " + label)


def agree_pair_differences(actual, expected, radii, label):
    require(isinstance(actual, list) and isinstance(expected, list) and len(radii) >= 2,
            label + " type")
    require(len(actual) == len(expected) == len(radii) - 1, label + " length")
    for index, want in enumerate(expected):
        agree_float(actual[index], want, label + "[" + str(index) + "]",
                    subtraction_abs_tol((radii[index], radii[index + 1])))


def agree_refinement_map(actual, expected, label):
    require(isinstance(actual, dict), label + " type")
    exact(sorted(actual), sorted(expected), label + " names")
    for name, row in expected.items():
        require(isinstance(actual[name], dict), label + " row")
        exact(sorted(actual[name]), sorted(row), label + " " + name + " fields")
        for key, value in row.items():
            target = label + "." + name + "." + key
            if key in REFINEMENT_DERIVED:
                agree_numbers(actual[name][key], value, target)
            else:
                exact(actual[name][key], value, target)


def agree_observable(actual, expected):
    require(isinstance(actual, dict), "physical observable type")
    exact(sorted(actual), sorted(expected), "physical observable fields")
    grid = expected["grid_values_metres"]
    envelope_abs = {
        "resolution_half_width_metres": subtraction_abs_tol(grid),
        "resolution_range_metres": subtraction_abs_tol(grid),
        "difference_from_quiet_reference_metres": subtraction_abs_tol(
            (expected["value_metres"], expected["quiet_reference_metres"])),
    }
    for key, value in expected.items():
        if key in envelope_abs:
            agree_numbers(actual[key], value, "physical_observable." + key, envelope_abs[key])
        elif key in OBSERVABLE_DERIVED:
            agree_numbers(actual[key], value, "physical_observable." + key)
        else:
            exact(actual[key], value, "physical_observable." + key)


def assert_comparison_matches(saved, recomputed):
    require(isinstance(saved, dict), "comparison type")
    require(set(COMPARISON_META) <= set(saved), "comparison provenance missing")
    body = {key: value for key, value in saved.items() if key not in COMPARISON_META}
    exact(sorted(body), sorted(recomputed), "comparison fields")
    for key in EXACT_COMPARISON_KEYS:
        exact(body[key], recomputed[key], "comparison " + key)
    agree_refinement_map(body["residual_maxima"], recomputed["residual_maxima"], "residual_maxima")
    agree_refinement_map(body["accumulated_transport_gaps"], recomputed["accumulated_transport_gaps"],
                         "accumulated_transport_gaps")
    observable = recomputed["physical_observable"]
    ell = observable["scale"]["ell_metres"]
    number(ell, "fiducial length")
    require(float(ell) > 0, "fiducial length must be positive")
    radii = [value / float(ell) for value in observable["grid_values_metres"]]
    agree_pair_differences(body["endpoint_radius_pair_differences"],
                           recomputed["endpoint_radius_pair_differences"], radii,
                           "endpoint_radius_pair_differences")
    agree_observable(body["physical_observable"], observable)
    require(re.fullmatch(r"[0-9a-f]{40,64}", saved["analysis_source_commit"]) is not None,
            "analysis commit must be hexadecimal")
    require(re.fullmatch(r"[0-9a-f]{64}", saved["analysis_source_sha256"]) is not None,
            "analysis source hash must be SHA256")
    require(saved["prediction_status"] in (ENVELOPE_STATUS, UNCONVERGED_STATUS),
            "unknown finite-amplitude status")
    passes = (saved["requested_mass_and_flux_contraction"]
              and saved["endpoint_radius_pair_difference_contracts"])
    expected_status = ENVELOPE_STATUS if passes else UNCONVERGED_STATUS
    exact(saved["prediction_status"], expected_status, "finite-amplitude status")


def load_comparison(root, binding):
    comparison = read_json(bound_bytes(root, binding, "finite-amplitude comparison"))
    require(isinstance(comparison.get("inputs"), list) and len(comparison["inputs"]) == 3,
            "three analysis inputs required")
    found = {}
    for row in comparison["inputs"]:
        require(isinstance(row, dict) and row.get("path") not in found, "duplicate analysis input")
        raw = bound_bytes(root, row, "analysis input")
        found[row["path"]] = (row, read_json(raw))
    return comparison, found


def endpoint_receipts(found, endpoints):
    require(isinstance(endpoints, list) and len(endpoints) == 3, "three endpoint bindings required")
    receipts = []
    for index, row in enumerate(endpoints):
        exact(sorted(row), ["bytes", "level", "path", "sha256"], "endpoint binding fields")
        require(type(row["level"]) is int and row["level"] == index, "endpoint levels must be 0, 1, 2")
        require(row["path"] in found, "endpoint path is not an analysis input: " + str(row["path"]))
        saved, record = found[row["path"]]
        exact({key: row[key] for key in ("bytes", "path", "sha256")}, saved, "endpoint analysis binding")
        exact(record["level"], row["level"], "endpoint level")
        receipts.append(record)
    exact(sorted(row["path"] for row in endpoints), sorted(found), "analysis input paths")
    return receipts


def finite_event():
    return dict(kind="fixed_incoming_label", affine_source_x=ENDPOINT_SOURCE_X, normalized_V=1.0,
                reference_radius_over_ell=ENDPOINT_REFERENCE_RADIUS,
                distinct_from_perturbative_rho=True)


def finite_amplitude_block(comparison, comparison_binding, endpoint_rows):
    observable = comparison["physical_observable"]
    exact(observable["event"]["affine_source_x"], ENDPOINT_SOURCE_X, "incoming-label event")
    exact(observable["event"]["normalized_V"], 1.0, "incoming-label normalized V")
    exact(observable["event"]["reference_radius_over_ell"], ENDPOINT_REFERENCE_RADIUS,
          "incoming-label reference radius")
    exact(observable["scale"]["ell_metres"], FIDUCIAL_LENGTH_METRES, "fiducial length")
    exact(observable["scale"]["role"], "declared_fiducial_input", "fiducial scale role")
    exact(observable["scale"]["measured_astrophysical_scale"], False, "measured scale claim")
    exact(observable["scale"]["fitted_to_output"], False, "fitted scale claim")
    exact(observable["observational_validation"], False, "observational validation claim")
    exact(observable["uncertainty"]["kind"], "observed_three_grid_resolution_envelope",
          "resolution-envelope kind")
    exact(observable["uncertainty"]["confidence_level"], None, "statistical confidence claim")
    exact(observable["uncertainty"]["rigorous_error_enclosure"], False, "total-error promotion")
    return dict(
        schema=FINITE_SCHEMA, event=finite_event(), comparison_status=comparison["status"],
        prediction_status=comparison["prediction_status"],
        requested_mass_and_flux_contraction=comparison["requested_mass_and_flux_contraction"],
        endpoint_radius_pair_difference_contracts=comparison["endpoint_radius_pair_difference_contracts"],
        shared_alpha=comparison["shared_alpha"], radial_accuracy=comparison["radial_accuracy"],
        physical_observable=observable, limits=dict(LIMITS),
        evidence=dict(comparison=comparison_binding, endpoints=endpoint_rows),
    )


def check_finite_amplitude(root, entry):
    block = entry["finite_amplitude"]
    require(isinstance(block, dict), "finite-amplitude extension type")
    exact(sorted(block), sorted(FINITE_FIELDS), "finite-amplitude fields")
    exact(block["schema"], FINITE_SCHEMA, "finite-amplitude schema")
    exact(block["event"], finite_event(), "finite-amplitude event")
    exact(block["limits"], LIMITS, "finite-amplitude limits")
    evidence = block["evidence"]
    exact(sorted(evidence), ["comparison", "endpoints"], "finite-amplitude evidence fields")
    comparison, found = load_comparison(root, evidence["comparison"])
    receipts = endpoint_receipts(found, evidence["endpoints"])
    recomputed = compare(receipts)
    assert_comparison_matches(comparison, recomputed)
    expected = finite_amplitude_block(comparison, evidence["comparison"], evidence["endpoints"])
    exact(block["comparison_status"], expected["comparison_status"], "comparison status copy")
    exact(block["prediction_status"], expected["prediction_status"], "finite-amplitude status copy")
    exact(block["requested_mass_and_flux_contraction"], expected["requested_mass_and_flux_contraction"],
          "contraction status copy")
    exact(block["endpoint_radius_pair_difference_contracts"],
          expected["endpoint_radius_pair_difference_contracts"], "radius-contraction status copy")
    exact(block["shared_alpha"], expected["shared_alpha"], "alpha copy")
    exact(block["radial_accuracy"], expected["radial_accuracy"], "radial accuracy copy")
    exact(block["physical_observable"], expected["physical_observable"], "physical observable copy")
    require(block["prediction_status"] != ENVELOPE_STATUS
            or (block["requested_mass_and_flux_contraction"]
                and block["endpoint_radius_pair_difference_contracts"]),
            "unconverged comparison labeled as a resolution-envelope estimate")


def check_p001(entry, records, root):
    exact(entry["registered_on"], "2026-10-16", "retrospective registration date")
    exact(entry["registration"], dict(kind="retrospective", evidence_previously_calculated=True,
          preregistered=False, withheld_observations=False), "registration status")
    exact(entry["status"], "closed_conditional_perturbative_result", "claim scope")
    model, cap, wall = (records[key] for key in ("model", "cap_wall", "wall_field"))
    exact(entry["model_id"], model["id"], "model identity")
    expected_parameters = dict(h="H_ell", b="k_ell", m="M_over_ell", xi="xi")
    exact(set_as_list(entry["parameters"]), sorted(expected_parameters), "parameter inventory")
    for name, source in expected_parameters.items():
        saved = model["parameters"][source]
        exact(entry["parameters"][name], dict(value=saved["value"], model_parameter=source,
                                             role=saved["role"]), "model parameter " + name)
        if name in cap["domain"]:
            require(entry["parameters"][name]["value"] == cap["domain"][name], "cap parameter differs")
    exact(entry["domain"], dict(reference_radius_over_ell=[.45,.55], branch="expanding_timelike",
        source="original_C3_compact_characteristic_preparation", quiet_cap_and_bottom=True,
        no_return=True, normal_orientation="cap_to_exterior", frame="Jordan"), "finite domain")
    exact(entry["domain"]["reference_radius_over_ell"], cap["domain"]["x"], "cap interval")
    exact(entry["domain"]["reference_radius_over_ell"], model["domain"]["wall_radius_over_ell"], "model interval")
    exact(model["preparation"]["source_regularity"], "C3_compact_polynomial", "original source")
    require(model["preparation"]["bottom_scalar"] == model["preparation"]["earlier_cap_scalar"] == 0,
            "quiet reference data required")
    exact(entry["expansion"]["alpha"], "kappa_4*epsilon^2/ell^2", "amplitude convention")
    exact(entry["expansion"]["units"], "dimensionless_coefficients", "coefficient units")
    exact(entry["prediction_point"], dict(reference_radius_over_ell=.55, normalized_V=1.,
        evidence="wall_field", pointer="/solutions/1/wall_samples/8"), "prediction point")
    sample = pointer(wall, entry["prediction_point"]["pointer"])
    exact(sample["x"], .55, "saved reference endpoint")
    exact(sample["t"], 1., "saved normalized V endpoint")
    exact(set_as_list(entry["coefficients"]), sorted(COEFFICIENTS), "coefficient inventory")
    for key in COEFFICIENTS:
        row = entry["coefficients"][key]
        expected = dict(value=sample[key], field_refinement_maximum_absolute_difference=
            wall["wall_differences"][key]["maximum_absolute_difference"])
        for field, value in row.items(): number(value, key + "." + field)
        exact(row, expected, key + " saved-value projection")
    exact(set_as_list(entry["identities"]), sorted(IDENTITY_POINTERS), "identity inventory")
    for name, path in IDENTITY_POINTERS.items():
        row = entry["identities"][name]
        exact(row["pointer"], path, name + " source pointer")
        exact(row["evidence"], "cap_wall", name + " evidence")
        exact(row["order"], "coefficient_of_alpha", name + " order")
        exact(row["value"], pointer(cap, path), name + " exact residual")
        exact(row["value"], "0", name + " zero identity")
    exact(entry["identities"]["delta_a_plus_minus_delta_a_minus"]["common_value_formula"],
          "6*m*rho/(b*x^4)-2*(mu_plus-nu_minus)/(b*x^3)", "acceleration formula")
    exact(entry["uncertainty"], dict(kind="sampled_field_grid_sensitivity", field_nodes=[129,257],
        sample_count=9, source_pointer="/wall_differences/{coefficient}/maximum_absolute_difference",
        scope="maximum over nine common cuts; not an endpoint-specific or total error bound",
        other_settings_held_fixed=True, rigorous_error_enclosure=False,
        finite_alpha_selected=False, finite_amplitude_remainder_bound=None), "uncertainty scope")
    exact(wall["refinement_is_rigorous_error_enclosure"], False, "saved uncertainty scope")
    exact([settings["field_nodes"] for settings in wall["settings"]], [129,257], "field refinement")
    exact([{k:v for k,v in settings.items() if k != "field_nodes"} for settings in wall["settings"]][0],
          {k:v for k,v in wall["settings"][1].items() if k != "field_nodes"}, "fixed settings")
    require(all(len(solution["wall_samples"]) == 9 for solution in wall["solutions"]), "nine cuts required")
    # The selected ell is a declared fiducial input. It is not a measured scale,
    # and it is recorded only together with the finite-amplitude extension.
    exact(entry["observational_test"], dict(status="pending", calibration="none",
        parameters_chosen_as_model_inputs=True, data_refs=[], physical_source_selected=False,
        SI_scale_selected="finite_amplitude" in entry, future_comparison_uncertainty=None),
          "observational status")
    exact(model["calibration"], dict(observation_refs=[], parameter_ids=[], status="not_performed"),
          "model calibration status")
    exact(entry["open_limits"], ["nonlinear_finite_amplitude_closure", "uniform_remainder_and_total_numerical_error",
        "physical_preparation_and_observational_identification", "formation_stability_and_multi_region_composition"],
          "remaining limits")
    require(set(entry["rejection"]) == {"fixed_assumption_identity_failure", "independent_prediction_test"}
            and all(isinstance(value, str) and value.strip() for value in entry["rejection"].values()),
            "explicit rejection conditions required")
    if "finite_amplitude" in entry:
        check_finite_amplitude(root, entry)


def check_register(root, register):
    root = Path(root).resolve()
    exact(register["schema"], SCHEMA, "register schema")
    exact(register["document"], DOCUMENT, "register document")
    document = local_file(root, DOCUMENT).read_text(encoding="utf-8")
    entries = register["entries"]
    require(isinstance(entries, list) and len(entries) == 1, "one currently registered prediction required")
    entry = entries[0]
    exact(entry["id"], "NSC-P001", "prediction identity")
    require("## NSC-P001" in document, "prediction missing from reader-facing register")
    records = authenticate(root, entry["evidence"])
    check_p001(entry, records, root)
    result = dict(predictions=1, evidence_files=len(EVIDENCE), coefficient_projections=4,
                  exact_identity_projections=2, scientific_replays=0)
    if "finite_amplitude" in entry:
        result["finite_amplitude_files"] = 4
        result["finite_amplitude_endpoints"] = 3
    return result


def check(root=ROOT):
    root = Path(root).resolve()
    register = read_json(local_file(root, REGISTER).read_bytes())
    return check_register(root, register)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args(argv)
    try:
        result = check(args.root)
    except (ValueError, OSError, KeyError, TypeError, IndexError) as error:
        print("prediction register check failed: " + str(error))
        return 1
    print("prediction register check passed: " + json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
