#!/usr/bin/env python3
"""Authenticate prediction declarations and saved-value projections; no replay."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re

ROOT = Path(__file__).resolve().parents[1]
REGISTER = "results/prediction-register.json"
DOCUMENT = "docs/prediction-register.md"
SCHEMA = "nsc-prediction-register-v1"
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


def authenticate(root, bindings):
    exact(set_as_list(bindings), sorted(EVIDENCE), "evidence inventory")
    records = {}
    for name, expected_path in EVIDENCE.items():
        row = bindings[name]
        exact(sorted(row), ["bytes", "path", "sha256"], name + " binding fields")
        exact(row["path"], expected_path, name + " path")
        require(type(row["bytes"]) is int and row["bytes"] > 0
                and isinstance(row["sha256"], str)
                and re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is not None,
                "invalid byte identity: " + name)
        raw = local_file(root, row["path"]).read_bytes()
        require(len(raw) == row["bytes"] and hashlib.sha256(raw).hexdigest() == row["sha256"],
                "evidence bytes changed: " + name)
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


def check_p001(entry, records):
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
    exact(entry["observational_test"], dict(status="pending", calibration="none",
        parameters_chosen_as_model_inputs=True, data_refs=[], physical_source_selected=False,
        SI_scale_selected=False, future_comparison_uncertainty=None), "observational status")
    exact(model["calibration"], dict(observation_refs=[], parameter_ids=[], status="not_performed"),
          "model calibration status")
    exact(entry["open_limits"], ["nonlinear_finite_amplitude_closure", "uniform_remainder_and_total_numerical_error",
        "physical_preparation_and_observational_identification", "formation_stability_and_multi_region_composition"],
          "remaining limits")
    require(set(entry["rejection"]) == {"fixed_assumption_identity_failure", "independent_prediction_test"}
            and all(isinstance(value, str) and value.strip() for value in entry["rejection"].values()),
            "explicit rejection conditions required")


def check(root=ROOT):
    root = Path(root).resolve()
    register = read_json(local_file(root, REGISTER).read_bytes())
    exact(register["schema"], SCHEMA, "register schema")
    exact(register["document"], DOCUMENT, "register document")
    document = local_file(root, DOCUMENT).read_text(encoding="utf-8")
    entries = register["entries"]
    require(isinstance(entries, list) and len(entries) == 1, "one currently registered prediction required")
    entry = entries[0]
    exact(entry["id"], "NSC-P001", "prediction identity")
    require("## NSC-P001" in document, "prediction missing from reader-facing register")
    records = authenticate(root, entry["evidence"])
    check_p001(entry, records)
    return dict(predictions=1, evidence_files=len(EVIDENCE), coefficient_projections=4,
                exact_identity_projections=2, scientific_replays=0)


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
