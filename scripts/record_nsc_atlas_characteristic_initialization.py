#!/usr/bin/env python3
"""Creation-only characteristic constraint record with explicit predecessor bytes.

--check authenticates stored inputs, bound sources and output without rerunning
quadratures. --check-numerics repeats the admitted bounded calculation.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

from validation_source_successors import verify_saved_bindings, verify_saved_source

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-characteristic-initialization-v1"
MANIFEST_SCHEMA = SCHEMA+"-manifest"
OWNER = "docs/nsc-atlas-characteristic-initialization.md"
PRODUCER = "scripts/record_nsc_atlas_characteristic_initialization.py"
MODULE = "src/recursive_horizons/nsc_atlas_characteristic_initialization.py"
TESTS = "tests/test_nsc_atlas_characteristic_initialization.py"
SOURCE_PATHS = (OWNER, PRODUCER, MODULE, TESTS)
RELATIVE_OUTPUT = "results/development/nsc-atlas-characteristic-initialization-v1"
DEFAULT = ROOT / RELATIVE_OUTPUT
ADAPTERS = tuple("src/recursive_horizons/"+name+".py" for name in (
    "nsc_atlas_causal_reference", "nsc_atlas_causal_observers", "nsc_atlas_connection_slice",
    "nsc_atlas_moving_boundaries", "nsc_atlas_field_preparation", "nsc_atlas_field_transfer"))
PREREQUISITES = {
    "v7": {"directory": "results/development/nsc-cosmological-atlas-v7",
           "manifest_sha256": "08e3f07a4c0a910852d1fed81949cf334b97bb3497d9ccb3c1aa611259d0bef2",
           "schema": "nsc-cosmological-atlas-v7", "outputs": 8,
           "source_commit": "7fee894ad5a4da2ba378f724e45e926401f82670"},
    "v8": {"directory": "results/development/nsc-cosmological-atlas-v8",
           "manifest_sha256": "c4d4d0d1b1ac667886eb5f20cbcb11205eead9975c7915bbce3bdb8c87aaadb5",
           "schema": "nsc-cosmological-atlas-v8", "outputs": 7,
           "source_commit": "1154baf290d901e73058392e143d0acb37d22951"},
    "scalar_wall": {"directory": "results/development/nsc-atlas-scalar-wall-consistency-v1",
                    "manifest_sha256": "61f4f2287b2f70196bf1a5fbc773bf529864514d38b46c80a44d9726a8d8b801",
                    "schema": "nsc-atlas-scalar-wall-consistency-manifest-v1", "outputs": 1},
}
IDENTITIES = ("improved_null_stress_constraint", "first_moment_derivative",
              "second_moment_derivative", "einstein_affine_focusing", "einstein_frame_moment")
SAMPLE_KEYS = {"V", "x", "p", "p_prime", "p_second", "h", "h_prime", "h_second",
               "h_direct", "h_prime_direct", "h_einstein", "integral_p_squared",
               "integral_x_dp_squared", "weighted_integral_x_dp_squared"}


def canonical_bytes(value):
    return (json.dumps(value, sort_keys=True, ensure_ascii=False, indent=2,
                       allow_nan=False)+"\n").encode("utf-8")


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def read_relative(relative):
    path = ROOT / relative
    if (Path(relative).is_absolute() or ".." in Path(relative).parts
            or not path.resolve().is_relative_to(ROOT.resolve()) or not path.is_file()):
        raise ValueError("missing or unsafe prerequisite/source: "+relative)
    return path.read_bytes()


def source_bindings():
    return {path: sha(read_relative(path)) for path in sorted(SOURCE_PATHS)}


def _require_saved_source(source, expected, label):
    try:
        verify_saved_source(source, expected, ROOT)
    except ValueError as error:
        raise ValueError(label) from error


def authenticate_prerequisites():
    records, bindings = {}, {}
    for key, specification in PREREQUISITES.items():
        relative = specification["directory"]+"/manifest.json"
        raw = read_relative(relative)
        if sha(raw) != specification["manifest_sha256"]:
            raise ValueError("pinned prerequisite manifest mismatch: "+key)
        manifest = json.loads(raw)
        if (manifest.get("schema") != specification["schema"]
                or len(manifest.get("outputs", {})) != specification["outputs"]):
            raise ValueError("prerequisite schema or inventory mismatch: "+key)
        payloads = {}
        for name, binding in manifest["outputs"].items():
            if Path(name).name != name or name in {".", ".."}:
                raise ValueError("unsafe prerequisite output: "+name)
            content = read_relative(specification["directory"]+"/"+name)
            if sha(content) != binding["sha256"] or len(content) != binding["bytes"]:
                raise ValueError("prerequisite artifact binding mismatch: "+key+"/"+name)
            if name.endswith(".json"):
                payloads[name] = json.loads(content)
        if key in ("v7", "v8") and manifest["input_sha256"] != manifest["outputs"]["atlas.json"]["sha256"]:
            raise ValueError("prerequisite atlas binding mismatch")
        if key == "v8":
            for adapter in ADAPTERS:
                _require_saved_source(adapter, manifest["source_bindings"][adapter],
                    "current adapter differs from explicitly borrowed v8 source: "+adapter)
        records[key] = payloads
        bindings[key] = {"manifest": relative, "manifest_sha256": sha(raw),
                         "authenticated_outputs": deepcopy(manifest["outputs"]),
                         "verification": "pinned manifests and all stored output bytes; no field replay"}
        if "source_commit" in specification:
            bindings[key]["source_commit"] = specification["source_commit"]
    return records, bindings


def declaration():
    records, bindings = authenticate_prerequisites()
    preparation = records["v7"]["preparation-access.json"]
    response = records["v8"]["prepared-response.json"]
    shell = records["scalar_wall"]["consistency.json"]
    if records["v8"]["preparation-access.json"] != preparation:
        raise ValueError("v8 preparation differs from original v7")
    if any(response[key] != preparation[key] for key in ("inputs", "preparation")):
        raise ValueError("prerequisite inputs or preparation differ")
    if shell["inputs"] != preparation["inputs"] or shell["domain"] != response["domain"]:
        raise ValueError("scalar-wall and field domains differ")
    if preparation["preparation"]["regularity"] != "C3 compact polynomial; derivative C2; not C-infinity":
        raise ValueError("the original C3 source is required")
    return {
        "schema": SCHEMA, "status": "QUADRATIC_SOURCE_LEG_NULL_CONSTRAINT_INITIALIZED",
        "owner": OWNER, "claim_class": "derived_and_numerically_evaluated_conditional",
        "inputs": deepcopy(preparation["inputs"]),
        "reference_domain": deepcopy(preparation["geometry"]),
        "borrowed": {"null_constraints": "2503.24162v3 Eqs. 16, 29a, 86-92",
                     "double_null_geometric_identity": "2005.04090v1 Eqs. 2.1, 2.6-2.7; its minimal matter source is not substituted for the conformal source",
                     "conformal_stress_and_frame": "1910.07534v2 Eqs. 2.11-2.14, 2.21-2.24",
                     "preparation": "original v7 C3 pulse, zero bottom and earlier cap; v8 leading field",
                     "methods": "one-dimensional Gauss quadrature, integration by parts and smooth-exterior frame transformation"},
        "initial_data": {
            "source": "U=U0, Va<=V<=Vb, retained Schwarzschild sector IV",
            "bottom": "V=Va, U0<=U<=Ua; scalar zero and unchanged reference null data",
            "corner": "same background radius and both reference null expansions; h(.4)=h'(.4)=0",
            "gauge": "A=A0 on source; x=R0/ell is the fixed corner-anchored affine coordinate there",
            "radius": "R/ell=x+alpha h(x)+O(alpha^2); radius is not independently fixed along the source",
            "field": "phi=(epsilon/ell)p(x)+O(epsilon^3); p=B(4(V(x)-Vmid)/DeltaV)/x",
            "parameter": "alpha=kappa_4 epsilon^2/ell^2; kappa_4=8 pi G in c=hbar=1 units",
            "amplitude": "no epsilon is selected; no finite-amplitude remainder bound is supplied",
            "theory": "canonical scalar kinetic coefficient one, xi=1/6, zero potential and scalar-independent exterior vacuum; no new wall scalar coupling",
        },
        "equations": {
            "exact_null_constraint": "(1-kappa_4 xi phi^2) R''=-(kappa_4 R/2)[(1-2xi)phi'^2-2xi phi phi'']; affine derivatives",
            "quadratic_constraint": "h''=-x p'^2/3+x p p''/6",
            "first_integral": "h'=-integral(s p'^2)/2+x p p'/6-p^2/12",
            "second_integral": "h=-integral((x-s)s p'^2)/2+x p^2/12-integral(p^2)/6",
            "einstein_affine_map": "xbar=x-alpha integral(p^2)/6+O(alpha^2)",
            "einstein_radius": "Rbar/ell=x+alpha(h-x p^2/12)+O(alpha^2)",
            "einstein_coefficient": "hE=h-x p^2/12+integral(p^2)/6=-integral((x-s)s p'^2)/2; hE''=-x p'^2/2",
            "integration_limits": "corner .4 to x, with integrands zero outside the original compact support",
        },
        "analytic_consequences": {
            "post_support": "h'=-integral(s p'^2)/2<0 and h=-integral((x-s)s p'^2)/2-integral(p^2)/6<0 for the nonzero compact source",
            "meaning": "pulse forces nonzero radius and affine radial-derivative corrections on the initial source leg",
            "source_tangent_stress": "after support phi=phi_V=0 on U0, so this leg's T_VV=0 and h' stays constant",
            "transverse_caveat": "chi_U may remain nonzero after source support; T_UU need not vanish; no vacuum or Schwarzschild post-pulse region is inferred",
            "frame_domain": "smooth exterior only with positive conformal factor; no transformation across the material wall is used",
            "unit_scaling": "h and h' are dimensionless at fixed geometry. At fixed alpha, delta R scales as ell. At fixed epsilon and G, alpha scales as ell^-2; this changes physical source strength, not merely units.",
        },
        "knowledge_boundary": {
            "full_characteristic_evolution_solved": False, "mass_shift_determined": False,
            "corrected_wall_solved": False, "interior_observer_correction_computed": False,
            "finite_amplitude_or_error_bound_established": False, "stability_established": False,
            "post_source_region_claimed_vacuum": False, "physical_source_generator_established": False,
            "h_prime_identified_with_null_expansion": False,
        },
        "remaining_equations": ["mixed characteristic Einstein equations and transverse metric data in the wedge",
                                "scalar response at the next perturbative order",
                                "perturbed material wall, matching and normal-force equations",
                                "observer clock and interior response after a coupled solution"],
        "next_unlock": "compatible source-leg null data for a separately admitted coupled exterior/wall problem",
        "verification_scope": "numerical refinement diagnoses integration; it is not a rigorous enclosure or a metric error bound",
    }, bindings


def validate_calculation(calculation, reference_domain):
    samples = calculation["samples"]
    if set(samples) != SAMPLE_KEYS or set(calculation["last"]) != SAMPLE_KEYS:
        raise ValueError("characteristic sample inventory mismatch")
    for key, values in samples.items():
        if (len(values) != 129 or not all(isinstance(v, (int, float)) and math.isfinite(v) for v in values)
                or calculation["last"][key] != values[-1]):
            raise ValueError("characteristic samples/last mismatch")
    if samples["h"][0] != 0 or samples["h_prime"][0] != 0 or not samples["h"][-1] < 0 or not samples["h_prime"][-1] < 0:
        raise ValueError("characteristic corner or analytic sign mismatch")
    geometry = calculation["geometry"]
    for key in ("U0", "Ua", "Ub", "Va", "Vb"):
        if not math.isclose(geometry[key], reference_domain[key], rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError("calculated characteristic differs from borrowed reference domain")
    if geometry["x_corner"] != .4 or not .4 < geometry["x_last"] < .55 or geometry["order"] != 64:
        raise ValueError("calculated source radius or geometry order differs")
    for i, x in enumerate(samples["x"]):
        p, dp = samples["p"][i], samples["p_prime"][i]
        I0, I1, I2 = (samples[key][i] for key in (
            "integral_p_squared", "integral_x_dp_squared", "weighted_integral_x_dp_squared"))
        for name, expected in (("h", -I2/2+x*p*p/12-I0/6),
                               ("h_prime", -I1/2+x*p*dp/6-p*p/12), ("h_einstein", -I2/2)):
            if abs(samples[name][i]-expected) > 1e-8+1e-8*abs(expected):
                raise ValueError("stored characteristic moment identity mismatch")
        for name in ("h", "h_prime"):
            if abs(samples[name][i]-samples[name+"_direct"][i]) > 1e-8+1e-8*abs(samples[name][i]):
                raise ValueError("stored direct characteristic identity mismatch")
    execution = calculation["execution"]
    if (execution["displayed_source_samples"] != 129
            or not 0 < execution["scalar_integrand_evaluations"] <= 20000
            or execution["geometric_integrand_evaluations"] != 1440
            or execution["total_integrand_evaluations"] != execution["scalar_integrand_evaluations"]+1440
            or execution["total_integrand_evaluations"] > 20000
            or execution["budget_guard_counts_geometry_too"] is not True
            or execution["maximum_scalar_integrand_evaluations"] != 20000
            or execution["cpu_seconds_max"] != 60 or execution["quadrature_orders"] != [32, 64]
            or execution["field_solves"] != 0 or execution["bulk_einstein_evolutions"] != 0):
        raise ValueError("characteristic execution budget/scope mismatch")
    checks = calculation["checks"]
    if (checks["refinement_is_rigorous_enclosure"] is not False
            or checks["tolerance_absolute"] != 1e-8 or checks["tolerance_relative"] != 1e-8):
        raise ValueError("characteristic refinement interpretation mismatch")


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only evidence: use --check or a fresh directory")
    from recursive_horizons.nsc_atlas_characteristic_initialization import (
        Budget, exact_algebra_checks, initialization_calculation,
    )
    budget = Budget()
    before = source_bindings()
    record, prerequisites = declaration()
    record["exact_symbolic_residuals"] = exact_algebra_checks()
    if record["exact_symbolic_residuals"] != dict.fromkeys(IDENTITIES, "0"):
        raise ValueError("exact characteristic certificate mismatch")
    record["calculation"] = initialization_calculation(budget)
    validate_calculation(record["calculation"], record["reference_domain"])
    if source_bindings() != before:
        raise RuntimeError("source changed during evidence recording")
    payload = canonical_bytes(record)
    manifest = {"schema": MANIFEST_SCHEMA, "source_bindings": before,
                "prerequisite_bindings": prerequisites,
                "outputs": {"initialization.json": {"sha256": sha(payload), "bytes": len(payload)}},
                "scope": "quadratic initial null constraint only; no field or bulk Einstein evolution"}
    budget.check()
    directory.mkdir(parents=True, exist_ok=False)
    for name, content in (("initialization.json", payload), ("manifest.json", canonical_bytes(manifest))):
        with (directory / name).open("xb") as handle:
            handle.write(content)
    return {"created": directory.name, "outputs": 1,
            "scalar_integrand_evaluations": record["calculation"]["execution"]["scalar_integrand_evaluations"],
            "total_integrand_evaluations": budget.evaluations, "cpu_guard_passed": True}


def check(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_bytes())
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("characteristic manifest schema mismatch")
    try:
        verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCE_PATHS)
    except ValueError as error:
        raise ValueError("characteristic source binding mismatch") from error
    header, prerequisites = declaration()
    if manifest.get("prerequisite_bindings") != prerequisites:
        raise ValueError("characteristic prerequisite binding mismatch")
    if (set(manifest.get("outputs", {})) != {"initialization.json"}
            or {p.name for p in directory.iterdir()} != {"manifest.json", "initialization.json"}):
        raise ValueError("characteristic output inventory mismatch")
    payload = (directory / "initialization.json").read_bytes()
    binding = manifest["outputs"]["initialization.json"]
    if sha(payload) != binding["sha256"] or len(payload) != binding["bytes"]:
        raise ValueError("characteristic artifact binding mismatch")
    record = json.loads(payload)
    if (set(record) != set(header) | {"calculation", "exact_symbolic_residuals"}
            or any(record[key] != value for key, value in header.items())
            or record["exact_symbolic_residuals"] != dict.fromkeys(IDENTITIES, "0")):
        raise ValueError("characteristic declaration mismatch")
    validate_calculation(record["calculation"], record["reference_domain"])
    return {"checked": directory.name, "outputs": 1, "sources": len(manifest["source_bindings"]),
            "authenticated_predecessor_outputs": sum(len(v["authenticated_outputs"]) for v in prerequisites.values()),
            "scope": "stored evidence authentication; no quadrature replay"}


def check_numerics(directory):
    check(directory)
    from recursive_horizons.nsc_atlas_characteristic_initialization import (
        Budget, initialization_calculation, exact_algebra_checks,
    )
    budget = Budget()
    if exact_algebra_checks() != dict.fromkeys(IDENTITIES, "0"):
        raise ValueError("exact characteristic certificate mismatch")
    fresh = initialization_calculation(budget)
    stored = json.loads((Path(directory) / "initialization.json").read_bytes())["calculation"]
    maximum = 0.
    for key in SAMPLE_KEYS:
        for a, b in zip(fresh["samples"][key], stored["samples"][key]):
            maximum = max(maximum, abs(a-b))
            if abs(a-b) > 1e-8+1e-8*abs(b):
                raise ValueError("characteristic numerical replay disagreement: "+key)
    budget.check()
    return {"replayed": Path(directory).name, "maximum_sample_difference": maximum,
            "scalar_integrand_evaluations": fresh["execution"]["scalar_integrand_evaluations"],
            "total_integrand_evaluations": budget.evaluations, "cpu_guard_passed": True}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--check-numerics", action="store_true")
    args = parser.parse_args(argv)
    result = check_numerics(args.output) if args.check_numerics else check(args.output) if args.check else create(args.output)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
