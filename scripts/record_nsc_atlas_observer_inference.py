#!/usr/bin/env python3
"""Record the analytic preparation inverse without replaying a field solver.

Predecessor manifests and every stored output are authenticated as immutable
inputs. Their source bindings are historical provenance, not a request to
compare historical owners with the current checkout or import old producers.
This checker authenticates the declaration and its proof; it is not a proof
assistant and does not numerically reconstruct a source.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-observer-inference-v1"
MANIFEST_SCHEMA = "nsc-atlas-observer-inference-manifest-v1"
OWNER = "docs/nsc-atlas-observer-inference.md"
PRODUCER = "scripts/record_nsc_atlas_observer_inference.py"
TESTS = "tests/test_nsc_atlas_observer_inference.py"
RELATIVE_OUTPUT = "results/development/nsc-atlas-observer-inference-v1"
DEFAULT = ROOT / RELATIVE_OUTPUT
PREREQUISITES = {
    "v7": {
        "directory": "results/development/nsc-cosmological-atlas-v7",
        "manifest_sha256": "08e3f07a4c0a910852d1fed81949cf334b97bb3497d9ccb3c1aa611259d0bef2",
        "source_commit": "7fee894ad5a4da2ba378f724e45e926401f82670",
        "outputs": 8,
    },
    "v8": {
        "directory": "results/development/nsc-cosmological-atlas-v8",
        "manifest_sha256": "c4d4d0d1b1ac667886eb5f20cbcb11205eead9975c7915bbce3bdb8c87aaadb5",
        "source_commit": "1154baf290d901e73058392e143d0acb37d22951",
        "outputs": 7,
    },
}


def canonical_bytes(value: dict) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True,
                       allow_nan=False) + "\n").encode("utf-8")


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def read_relative(relative: str) -> bytes:
    path = ROOT / relative
    if (Path(relative).is_absolute() or ".." in Path(relative).parts
            or not path.resolve().is_relative_to(ROOT.resolve()) or not path.is_file()):
        raise ValueError("missing or unsafe prerequisite/source: " + relative)
    return path.read_bytes()


def source_bindings() -> dict[str, str]:
    return {path: sha(read_relative(path)) for path in sorted((OWNER, PRODUCER, TESTS))}


def authenticate_prerequisites() -> tuple[dict, dict]:
    records, bindings = {}, {}
    for key, specification in PREREQUISITES.items():
        directory = specification["directory"]
        manifest_path = directory + "/manifest.json"
        raw = read_relative(manifest_path)
        if sha(raw) != specification["manifest_sha256"]:
            raise ValueError("pinned prerequisite manifest mismatch: " + key)
        manifest = json.loads(raw)
        if (manifest.get("schema") != "nsc-cosmological-atlas-" + key
                or len(manifest.get("outputs", {})) != specification["outputs"]):
            raise ValueError("prerequisite schema or inventory mismatch: " + key)
        payloads = {}
        for name, binding in manifest["outputs"].items():
            if Path(name).name != name or name in {".", ".."}:
                raise ValueError("unsafe prerequisite output: " + name)
            raw_output = read_relative(directory + "/" + name)
            if sha(raw_output) != binding["sha256"] or len(raw_output) != binding["bytes"]:
                raise ValueError("prerequisite artifact binding mismatch: " + key + "/" + name)
            if name.endswith(".json"):
                payloads[name] = json.loads(raw_output)
        if manifest["input_sha256"] != manifest["outputs"]["atlas.json"]["sha256"]:
            raise ValueError("prerequisite atlas binding mismatch: " + key)
        records[key] = payloads
        bindings[key] = {
            "manifest": manifest_path,
            "manifest_sha256": specification["manifest_sha256"],
            "source_commit": specification["source_commit"],
            "authenticated_outputs": deepcopy(manifest["outputs"]),
            "verification": "stored manifest and all output bytes; no source replay or numerical solve",
        }
    return records, bindings


def build_record() -> tuple[dict, dict]:
    records, prerequisites = authenticate_prerequisites()
    preparation = records["v7"]["preparation-access.json"]
    response = records["v8"]["prepared-response.json"]
    if records["v8"]["preparation-access.json"] != preparation:
        raise ValueError("v8 no longer retains the v7 preparation certificate")
    if (preparation["status"] != "FINITE_RETAINED_IV_PREPARATION_AND_OPERATOR_CERTIFIED"
            or response["status"] != "FINITE_PREPARED_CLASSICAL_RESPONSE_WITH_POSITIVE_TAIL"):
        raise ValueError("prerequisite result status differs")
    if any(response[key] != preparation[key] for key in ("inputs", "preparation")):
        raise ValueError("prerequisite field inputs or preparation differ")
    geometry, observer, operator = (preparation[key] for key in ("geometry", "observer", "operator"))
    if any(response["domain"][key] != value for key, value in geometry.items()):
        raise ValueError("prerequisite retained domains differ")
    if response["domain"]["arrival_window_over_ell"] != observer["arrival_window_over_ell"]:
        raise ValueError("prerequisite observer windows differ")
    if (not 0 <= operator["analytic_norm_bound"] < 1
            or not preparation["checks"]["analytic_contraction_independent_of_quadrature"]
            or not preparation["checks"]["analytic_no_return_independent_of_quadrature"]
            or observer["analytic_no_return_log_margin_lower_bound"] <= 0
            or not preparation["preparation"]["corner_compatible"]):
        raise ValueError("the stored inverse prerequisites are not certified")
    scenario = records["v8"]["atlas.json"]["scenarios"][-1]
    if scenario["id"] != "finite_iv_prepared_response":
        raise ValueError("the predecessor scenario is not the fixed prepared response")
    connection = deepcopy(scenario["connection_records"][0])
    connection.update(
        id="finite-IV-preparation-inverse", status="conditional", waveform_evaluated=False,
        assumptions=["fixed v7/v8 geometry and linear conformal scalar law",
                     "C2 source and wall traces vanish on the fixed initial segment",
                     "zero initial cap, no centre return, and the known monotone observer clock"],
        owner_refs=[OWNER], evidence_refs=[RELATIVE_OUTPUT + "/inference.json"],
        independent_check_refs=[OWNER], prerequisite_ids=["observer-response"],
        unlock_ids=["observer-inference"],
        missing=["an independently constrained physical preparation law",
                 "a second compatible geometry and a source-linked discriminator"],
        prescribed_boundary="U=U0 in retained IV; the inverse treats the compatible preparation class",
    )
    record = {
        "schema": SCHEMA,
        "status": "CONDITIONAL_FIXED_OPERATOR_PREPARATION_INVERSE",
        "classification": "derived_conditional",
        "owner": OWNER,
        "predecessor_scenario": scenario["id"],
        "prerequisite_records": {
            "preparation": PREREQUISITES["v7"]["directory"] + "/preparation-access.json",
            "response": PREREQUISITES["v8"]["directory"] + "/prepared-response.json",
            "programme": PREREQUISITES["v8"]["directory"] + "/programme.json",
        },
        "inputs": deepcopy(preparation["inputs"]),
        "domain": deepcopy(response["domain"]),
        "preparation_class": {
            "regularity": "g,h in C2([Va,Vb]); Y in C1([theta_a,theta_b])",
            "regularity_scope": "C2 reduced field on the retained exterior wedge and C1 central observable; full four-dimensional scalar C2 at the centre is not asserted for every C2 input",
            "initial_zero_segment": [geometry["Va"], preparation["preparation"]["support_V"][0]],
            "initial_zero_segment_origin": "Vs equals the stored v7 input support's left endpoint",
            "source_and_wall_condition": "g=h=0 on [Va,Vs]",
            "observer_condition": "Y=0 on theta([Va,Vs]); theta is the fixed monotone clock",
            "late_support_restriction": None,
            "sign_restriction": None,
            "original_pulse_evidence_unchanged": True,
            "inverse_preserves_original_pulse_subclass": False,
        },
        "borrowed_certificate": {
            "analytic_norm_bound": operator["analytic_norm_bound"],
            "bound_proof": operator["bound_proof"],
            "kernel_equation": operator["equation"],
            "wall_condition": operator["wall_condition"],
            "bottom_data": preparation["preparation"]["bottom_data"],
            "cap_data": preparation["preparation"]["cap_data"],
            "observer": deepcopy(observer),
        },
        "analytic_result": {
            "field_normalization": "f=chi/epsilon; g and h are normalized source and wall traces",
            "observer_normalization": "theta=T/ell; Y=ell phi_C/epsilon",
            "forward": "f=Eg-Af; (Eg)(U,V)=g(V)",
            "inverse_integral": "(Bf)(U,V)=integral[U,Uw(V)]dp integral[L(p),V]dq K(p,q)f(p,q)",
            "inverse": "f=Eh+Bf=(I-B)^(-1)Eh; g(V)=f(U0,V)",
            "lower_limit": operator["lower_limit"],
            "inverse_norm_bound": "||B||<=q_*<1; ||g||_infinity<=||h||_infinity/(1-q_*)",
            "observer_to_wall": "h(V)=(1/2) integral[theta_a,theta(V)] Y(s) ds",
            "observer_to_source_stability": "||g1-g2||_infinity<=(theta_b-theta_a)||Y1-Y2||_infinity/[2(1-q_*)]",
            "fixed_operator_bijection": "the declared g, h and Y classes are in bijection",
            "verification_method": "differentiate the inverse equation to recover PDE, wall, source and bottom; use the independent Q algebraic factorization",
            "two_model_implication": "a shared observed trace gives a conditional range-overlap question after clocks, normalization and allowed preparations are specified",
        },
        "preparation_class_test": {
            "source_cutoff": preparation["preparation"]["support_V"][1],
            "range_membership": "S(I-B)^(-1)E[(1/2) integral[theta_a,theta(V)]Y(s)ds] belongs to the physically permitted source class P",
            "restricted_source_class": "nonzero nonnegative C2 source supported in [Vs,Ve]",
            "post_input_conditions": ["h(V)<0", "h'(V)>0", "Y(theta(V))>0",
                                      "integral[theta_a,theta(V)]Y(s)ds=2h(V)<0"],
            "post_input_domain": "Ve<V<=Vb only",
            "exclusion_witness": "a smooth nonzero nonnegative observer bump supported strictly between theta(Vs) and theta(Ve)",
            "exclusion_argument": "its integrated h is nonnegative and becomes a positive constant; positive B implies f>=Eh and g>=h>0 after Ve, violating the source cutoff",
            "scope": "exact within-model class exclusion, not a simulated waveform or a second physical geometry",
        },
        "independent_analytic_review": {
            "status": "independently_checked_in_declared_domain",
            "reference": OWNER + "#borrowing-independent-review-and-evidence",
            "checked": ["B inverse recovers PDE, wall, source and bottom",
                        "Q algebraic cross-check and central normalization",
                        "C2 source/wall and C1 observer classes with the initially quiet corner",
                        "positive compact-source post-input signs",
                        "positive observer-bump exclusion of the late source cutoff"],
            "limit": "analytic review of this finite application; provenance tests do not prove the theorem or audit all cited external results",
        },
        "connection_record": connection,
        "knowledge_boundary": {
            "fixed_operator_source_inverse": True,
            "numerically_reconstructed_source": False,
            "second_physical_geometry_constructed": False,
            "observer_equivalence_established": False,
            "geometry_identified_from_one_unknown_preparation": False,
            "right_parent_causal_exclusion_changed": False,
            "physical_preparation_law_established": False,
        },
        "next_unlock": "constrain the preparation independently before comparing source-linked observables or a second compatible geometry",
        "execution": {
            "scientific_numerical_solves": 0,
            "new_waveforms": 0,
            "numerically_evaluated_inverse": False,
            "new_scientific_numbers": False,
            "check_scope": "provenance and exact analytic declaration, not machine verification of theorem truth",
        },
        "limits": [
            "The original compact positive pulse is one retained preparation, not the full inverse class.",
            "No inverse-preserved late compact support or sign condition is asserted.",
            "A fixed-model source inverse does not identify the geometry or establish equivalence of two physical models.",
            "Zero cap and the known observer clock fix the otherwise free integration constant and labels.",
            "C2 regularity of the full four-dimensional scalar at the centre is not claimed for every C2 reduced input; the original C3 pulse and smooth witness have stronger regularity.",
            "No detector dynamics, physical source generation, late returns or global evolution are supplied.",
        ],
    }
    return record, prerequisites


def create(directory: Path) -> dict:
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only evidence: use --check or a fresh directory")
    before = source_bindings()
    record, prerequisites = build_record()
    payload = canonical_bytes(record)
    if source_bindings() != before:
        raise RuntimeError("source changed during evidence recording")
    manifest = {
        "schema": MANIFEST_SCHEMA,
        "source_bindings": before,
        "prerequisite_bindings": prerequisites,
        "outputs": {"inference.json": {"sha256": sha(payload), "bytes": len(payload)}},
        "scope": "analytic declaration and proof binding; no numerical field replay",
    }
    directory.mkdir(parents=True, exist_ok=False)
    with (directory / "inference.json").open("xb") as handle:
        handle.write(payload)
    with (directory / "manifest.json").open("xb") as handle:
        handle.write(canonical_bytes(manifest))
    return {"created": directory.name, "outputs": 1, "scientific_numerical_solves": 0}


def check(directory: Path) -> dict:
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_bytes())
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("inference manifest schema mismatch")
    if manifest.get("source_bindings") != source_bindings():
        raise ValueError("inference source binding mismatch")
    record, prerequisites = build_record()
    if manifest.get("prerequisite_bindings") != prerequisites:
        raise ValueError("inference prerequisite binding mismatch")
    if (set(manifest.get("outputs", {})) != {"inference.json"}
            or {path.name for path in directory.iterdir()} != {"inference.json", "manifest.json"}):
        raise ValueError("inference output inventory mismatch")
    payload = (directory / "inference.json").read_bytes()
    binding = manifest["outputs"]["inference.json"]
    if sha(payload) != binding["sha256"] or len(payload) != binding["bytes"]:
        raise ValueError("inference artifact binding mismatch")
    if payload != canonical_bytes(record):
        raise ValueError("inference declaration mismatch")
    return {"checked": directory.name, "outputs": 1, "sources": len(manifest["source_bindings"]),
            "authenticated_predecessor_outputs": sum(len(row["authenticated_outputs"]) for row in prerequisites.values()),
            "scientific_numerical_solves": 0}


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    print(json.dumps(check(args.output) if args.check else create(args.output), sort_keys=True))


if __name__ == "__main__":
    main()
