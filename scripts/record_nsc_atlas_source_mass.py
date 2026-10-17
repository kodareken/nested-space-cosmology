#!/usr/bin/env python3
"""Seal source-leg invariant mass evidence; --check does not rerun quadratures."""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path

from validation_source_successors import verify_saved_bindings, verify_saved_source

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-source-mass-v1"
MANIFEST_SCHEMA = SCHEMA+"-manifest"
OWNER = "docs/nsc-atlas-source-mass.md"
PRODUCER = "scripts/record_nsc_atlas_source_mass.py"
MODULE = "src/recursive_horizons/nsc_atlas_source_mass.py"
TESTS = "tests/test_nsc_atlas_source_mass.py"
SOURCE_PATHS = (OWNER, PRODUCER, MODULE, TESTS)
RELATIVE_OUTPUT = "results/development/nsc-atlas-source-mass-v1"
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
    "initialization": {"directory": "results/development/nsc-atlas-characteristic-initialization-v1",
                       "manifest_sha256": "bb1f0bdf9a2d651f5a2bb7de7531d2c3006521dadbd00fbc49778b17bbe70c65",
                       "schema": "nsc-atlas-characteristic-initialization-v1-manifest", "outputs": 1},
}
IDENTITIES = ("background_transverse_constraint", "transverse_scalar_transport",
              "original_kernel_change_of_variable", "mixed_radius_integral", "mixed_radius_equation",
              "jordan_mixed_improved_stress", "mixed_and_frame_mass", "jordan_improved_stress_mass_transport",
              "post_support_mass_from_radius", "null_boost_kernel_normalization")
CONTROLS = ("omitting_transverse_frame_term_changes_mass", "dropping_improved_mixed_stress_changes_transport")
SAMPLE_KEYS = {"x", "V", "p", "p_prime", "p_second", "h", "h_prime", "d", "s", "s_kernel",
               "s_prime", "k", "mu_einstein", "mu_jordan", "mu_geometric", "frame_mass_difference",
               "mu_jordan_prime_stress", "mu_jordan_prime_frame", "T_xx_normalized", "T_Ux_over_C_normalized",
               "integral_p_over_x_squared", "integral_p_squared", "integral_x_a_minus_x_dp_squared",
               "integral_V_kernel_b"}


def canonical_bytes(value):
    return (json.dumps(value, sort_keys=True, ensure_ascii=False, indent=2, allow_nan=False)+"\n").encode("utf-8")


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
        borrowed = ADAPTERS if key == "v8" else tuple(manifest["source_bindings"]) if key == "initialization" else ()
        for source in borrowed:
            _require_saved_source(source, manifest["source_bindings"][source],
                "current source differs from explicitly borrowed context: "+source)
        records[key] = payloads
        bindings[key] = {"manifest": relative, "manifest_sha256": sha(raw),
                         "authenticated_outputs": deepcopy(manifest["outputs"]),
                         "verification": "pinned manifests and all output bytes; borrowed adapters and initialization sources authenticated; no field replay"}
        if "source_commit" in specification:
            bindings[key]["source_commit"] = specification["source_commit"]
    return records, bindings


def expected_certificate():
    return {"zero_residuals": dict.fromkeys(IDENTITIES, "0"), "nonzero_controls": dict.fromkeys(CONTROLS, True)}


def declaration():
    records, bindings = authenticate_prerequisites()
    preparation = records["v7"]["preparation-access.json"]
    response = records["v8"]["prepared-response.json"]
    initialization = records["initialization"]["initialization.json"]
    shell = records["scalar_wall"]["consistency.json"]
    if records["v8"]["preparation-access.json"] != preparation:
        raise ValueError("v8 preparation differs from original v7")
    if any(response[key] != preparation[key] for key in ("inputs", "preparation")):
        raise ValueError("prerequisite inputs or preparation differ")
    if (initialization["inputs"] != preparation["inputs"] or shell["inputs"] != preparation["inputs"]
            or initialization["reference_domain"] != preparation["geometry"] or shell["domain"] != response["domain"]):
        raise ValueError("initialization, scalar-wall and field domains differ")
    if preparation["preparation"]["regularity"] != "C3 compact polynomial; derivative C2; not C-infinity":
        raise ValueError("the original C3 source is required")
    return {
        "schema": SCHEMA, "status": "QUADRATIC_SOURCE_LEG_INVARIANT_MASS_RESPONSE",
        "owner": OWNER, "claim_class": "derived_and_numerically_evaluated_conditional",
        "inputs": deepcopy(preparation["inputs"]), "reference_domain": deepcopy(initialization["calculation"]["geometry"]),
        "borrowed": {"mass_and_mixed_equations": "2005.04090v1 Eqs. 1.3, 2.4, 2.14-2.15; canonical normalization restored",
                     "geometric_mass_and_transverse_transport": "2503.24162v3 Eqs. 21, 84c, 91",
                     "general_matter_mass_law": "Hayward, gr-qc/9710089v2 Eqs. 2.9d, 2.10b; spherical geometric mass with full conformal stress",
                     "frame_transformation": "Faraoni-Vitagliano, 1401.1189v1 Eq. 7; both null derivatives retained",
                     "scalar_frame": "1910.07534v2; smooth-exterior conformal field and metric transformation",
                     "numerics": "Gauss 32/64 quadrature, analytic source derivatives and unchanged initialization adapters"},
        "conventions": {
            "region": "same finite source U=U0 in retained white-hole sector IV; no added future domain",
            "metric": "ds^2/ell^2=-C dU dx+r^2 dOmega^2 on the source, C=-4a/U0, a=.8, m=.4",
            "initial_data": "same C3 source, zero bottom, original corner radius and both expansions, fixed source lapse; zero corner transverse correction",
            "field": "phi=(epsilon/ell)p+O(epsilon^3); s=ell phi_U/(C epsilon) at leading order",
            "parameter": "alpha=kappa_4 epsilon^2/ell^2; kappa_4=8 pi G in c=hbar=1 units",
            "mass": "geometric mass length mathfrak_m=R(1-g^ab partial_a R partial_b R)/2; mathfrak_m_J/ell=.4+alpha mu_jordan+O(alpha^2)",
            "interpretation": "scalar under coordinate changes; not generally conformally invariant, asymptotic energy or whole-exterior mass",
            "frame_domain": "smooth exterior, positive conformal factor; no transformed material wall or cap sources used",
            "amplitude": "formal quadratic coefficient; no nonzero amplitude or perturbative remainder bound selected",
        },
        "equations": {
            "background_transverse": "d=r0_U/C=(a/x-1)/4>0; (x d)'=-1/4",
            "transverse_field": "s' = -(s+d p')/x; s=-d p/x-a/(4x) integral(p/t^2); s(.4)=0",
            "original_kernel_check": "s=(-integral_V(K b)/C-d p)/x; original chi_U/epsilon=-integral_V(K b)",
            "minimal_frame_mass": "mu_einstein=-integral[t(a-t)p'^2]/4; mu_einstein'=-x^2 d p'^2",
            "jordan_frame_mass": "mu_jordan=mu_einstein+m p^2/12+x^2 p(s+d p')/3",
            "transverse_geometry": "k=-d h/x+I0/(24x)+d p^2/6+x p s/6; r_U/C=d+alpha k+O(alpha^2)",
            "geometric_check": "mu_geometric=m h/x+2x(k+d h')",
            "jordan_mixed_check": "(x k+d h)'=x^2(p' s-p s')/6",
            "improved_stress_transport": "mu_jordan'=x^2[(p' s-p s')/3-d(2p'^2-p p'')/3]",
            "post_source": "mu_jordan=mu_einstein=-integral_support[t(a-t)p'^2]/4<0; constant on remaining finite source leg",
        },
        "meaning": {
            "sign": "0<t<a and d>0 in this sector IV; positive minimal null stress gives the declared negative coefficient, not an untrapped-exterior monotonicity theorem",
            "coordinate_boost": "Unew=bU,Vnew=V/b,Cnew=C/b and phi_Unew=phi_U/b; s and geometric mass unchanged",
            "unit_scaling": "dimensionless coefficient fixed; formal mass-length correction ell alpha mu. Changing ell at fixed epsilon and G changes alpha rather than only units",
            "post_source": "vanishing p and p' makes the two frame masses coincide here; nonzero transverse s can persist, so no vacuum future is inferred",
        },
        "knowledge_boundary": {
            "whole_exterior_mass_determined": False, "asymptotic_energy_computed": False,
            "bulk_einstein_evolution_solved": False, "corrected_wall_solved": False,
            "interior_observer_correction_computed": False, "finite_amplitude_validity_bound_established": False,
            "post_source_region_claimed_vacuum": False, "stability_established": False,
            "physical_source_generator_established": False,
        },
        "next_unlock": "transport initialized sourced geometry toward the moving wall, with a separately admitted domain and dynamical junction check",
        "verification_scope": "exact local identities and numerical refinement diagnostics; no rigorous quadrature enclosure or finite-amplitude error bound",
    }, bindings


def validate_calculation(calculation, reference):
    samples = calculation["samples"]
    if set(samples) != SAMPLE_KEYS or set(calculation["last"]) != SAMPLE_KEYS:
        raise ValueError("source mass sample inventory mismatch")
    for key, values in samples.items():
        if (len(values) != 129 or not all(isinstance(v, (int, float)) and math.isfinite(v) for v in values)
                or calculation["last"][key] != values[-1]):
            raise ValueError("source mass samples/last mismatch")
    for key in ("U0", "Ua", "Ub", "Va", "Vb", "x_corner", "x_last"):
        if not math.isclose(calculation["geometry"][key], reference[key], rel_tol=1e-10, abs_tol=1e-10):
            raise ValueError("source mass geometry differs from initialized reference")
    if not calculation["last"]["mu_jordan"] < 0 or samples["mu_jordan"][0] != 0 or samples["k"][0] != 0:
        raise ValueError("source mass corner or analytic sign mismatch")
    for left, right in (("s", "s_kernel"), ("mu_jordan", "mu_geometric"),
                        ("mu_jordan_prime_frame", "mu_jordan_prime_stress")):
        if any(abs(a-b) > 1e-8+1e-8*abs(a) for a, b in zip(samples[left], samples[right])):
            raise ValueError("stored source mass independent identity mismatch")
    for i, x in enumerate(samples["x"]):
        p, dp, s = (samples[key][i] for key in ("p", "p_prime", "s"))
        expected = samples["mu_einstein"][i]+.4*p*p/12+x*x*p*(s+samples["d"][i]*dp)/3
        if abs(samples["mu_jordan"][i]-expected) > 1e-8+1e-8*abs(expected):
            raise ValueError("stored source mass frame relation mismatch")
    e = calculation["execution"]
    count_keys = ("geometry_integrand_evaluations", "reused_constraint_integrand_evaluations",
                  "mass_and_transverse_x_integrand_evaluations", "original_V_kernel_integrand_evaluations")
    if (e["displayed_source_samples"] != 129 or e["quadrature_orders"] != [32, 64]
            or e["total_integrand_evaluations"] != sum(e[k] for k in count_keys)
            or not 0 < e["total_integrand_evaluations"] <= 100000
            or e["maximum_integrand_evaluations"] != 100000 or e["cpu_seconds_max"] != 60
            or e["field_solves"] != 0 or e["bulk_einstein_evolutions"] != 0):
        raise ValueError("source mass budget/scope mismatch")
    checks = calculation["checks"]
    if (checks["refinement_is_rigorous_enclosure"] is not False
            or checks["tolerance_absolute"] != 1e-8 or checks["tolerance_relative"] != 1e-8):
        raise ValueError("source mass refinement interpretation mismatch")


def create(directory):
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only evidence: use --check or a fresh directory")
    from recursive_horizons.nsc_atlas_source_mass import Budget, exact_algebra_checks, mass_calculation, MAX_EVALUATIONS
    budget = Budget(maximum_evaluations=MAX_EVALUATIONS)
    before = source_bindings()
    record, prerequisites = declaration()
    record["exact_symbolic_certificate"] = exact_algebra_checks()
    if record["exact_symbolic_certificate"] != expected_certificate():
        raise ValueError("exact source mass certificate mismatch")
    record["calculation"] = mass_calculation(budget)
    validate_calculation(record["calculation"], record["reference_domain"])
    if source_bindings() != before:
        raise RuntimeError("source changed during evidence recording")
    payload = canonical_bytes(record)
    manifest = {"schema": MANIFEST_SCHEMA, "source_bindings": before, "prerequisite_bindings": prerequisites,
                "outputs": {"mass.json": {"sha256": sha(payload), "bytes": len(payload)}},
                "scope": "quadratic source-leg geometric mass only; no coupled bulk or wall evolution"}
    budget.check()
    directory.mkdir(parents=True, exist_ok=False)
    for name, content in (("mass.json", payload), ("manifest.json", canonical_bytes(manifest))):
        with (directory / name).open("xb") as handle:
            handle.write(content)
    return {"created": directory.name, "outputs": 1, "total_integrand_evaluations": budget.evaluations,
            "cpu_guard_passed": True}


def check(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_bytes())
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("source mass manifest schema mismatch")
    try:
        verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCE_PATHS)
    except ValueError as error:
        raise ValueError("source mass source binding mismatch") from error
    header, prerequisites = declaration()
    if manifest.get("prerequisite_bindings") != prerequisites:
        raise ValueError("source mass prerequisite binding mismatch")
    if (set(manifest.get("outputs", {})) != {"mass.json"}
            or {p.name for p in directory.iterdir()} != {"manifest.json", "mass.json"}):
        raise ValueError("source mass output inventory mismatch")
    payload = (directory / "mass.json").read_bytes()
    binding = manifest["outputs"]["mass.json"]
    if sha(payload) != binding["sha256"] or len(payload) != binding["bytes"]:
        raise ValueError("source mass artifact binding mismatch")
    record = json.loads(payload)
    if (set(record) != set(header) | {"calculation", "exact_symbolic_certificate"}
            or any(record[key] != value for key, value in header.items())
            or record["exact_symbolic_certificate"] != expected_certificate()):
        raise ValueError("source mass declaration mismatch")
    validate_calculation(record["calculation"], record["reference_domain"])
    return {"checked": directory.name, "outputs": 1, "sources": len(manifest["source_bindings"]),
            "authenticated_predecessor_outputs": sum(len(v["authenticated_outputs"]) for v in prerequisites.values()),
            "scope": "stored evidence authentication; no quadrature replay"}


def check_numerics(directory):
    check(directory)
    from recursive_horizons.nsc_atlas_source_mass import Budget, mass_calculation, exact_algebra_checks, MAX_EVALUATIONS
    budget = Budget(maximum_evaluations=MAX_EVALUATIONS)
    if exact_algebra_checks() != expected_certificate():
        raise ValueError("exact source mass certificate mismatch")
    fresh = mass_calculation(budget)
    stored = json.loads((Path(directory) / "mass.json").read_bytes())["calculation"]
    maximum = 0.
    for key in SAMPLE_KEYS:
        for a, b in zip(fresh["samples"][key], stored["samples"][key]):
            maximum = max(maximum, abs(a-b))
            if abs(a-b) > 1e-8+1e-8*abs(b):
                raise ValueError("source mass numerical replay disagreement: "+key)
    budget.check()
    return {"replayed": Path(directory).name, "maximum_sample_difference": maximum,
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
