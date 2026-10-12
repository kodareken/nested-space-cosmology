#!/usr/bin/env python3
"""Record finite scalar/wall compatibility without replaying the field solver.

Creation runs exact symbolic checks; ordinary --check authenticates their bound
declaration and sources using only the standard library. --check-algebra repeats
the symbolic identities. Neither route solves the bulk Einstein equations.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-scalar-wall-consistency-v1"
MANIFEST_SCHEMA = "nsc-atlas-scalar-wall-consistency-manifest-v1"
OWNER = "docs/nsc-atlas-scalar-wall-consistency.md"
PRODUCER = "scripts/record_nsc_atlas_scalar_wall_consistency.py"
TESTS = "tests/test_nsc_atlas_scalar_wall_consistency.py"
RELATIVE_OUTPUT = "results/development/nsc-atlas-scalar-wall-consistency-v1"
DEFAULT = ROOT / RELATIVE_OUTPUT
SYMPY_VERSION = "1.14.0"
PREREQUISITES = {
    "v7": {
        "directory": "results/development/nsc-cosmological-atlas-v7",
        "manifest_sha256": "08e3f07a4c0a910852d1fed81949cf334b97bb3497d9ccb3c1aa611259d0bef2",
        "source_commit": "7fee894ad5a4da2ba378f724e45e926401f82670", "outputs": 8,
    },
    "v8": {
        "directory": "results/development/nsc-cosmological-atlas-v8",
        "manifest_sha256": "c4d4d0d1b1ac667886eb5f20cbcb11205eead9975c7915bbce3bdb8c87aaadb5",
        "source_commit": "1154baf290d901e73058392e143d0acb37d22951", "outputs": 7,
    },
}
FORMULAS = {
    "surface_coefficient": "xi*(2*phi*d - 2*k*phi**2)",
    "surface_after_scalar_jump": "2*xi*(6*xi - 1)*k*phi**2",
    "amm_material_coefficient": "-4*k*f + 2*d*fp",
    "flux_jump": "(1 - 2*xi)*d*v - 2*xi*phi*dd - 2*xi*k*phi*v",
    "flux_after_scalar_jump": "4*xi*(1 - 6*xi)*k*phi*v",
    "surface_tangential_derivative": "xi*(2*v*d + 2*phi*dd - 4*k*phi*v)",
}
IDENTITY_NAMES = (
    "improved_surface_matches_component_distribution",
    "surface_is_negative_amm_residual",
    "surface_after_scalar_jump",
    "flux_matches_product_rule_components",
    "flux_after_scalar_jump",
    "general_coupling_surface_flux_balance",
    "conformal_surface_vanishes",
    "conformal_amm_residual_vanishes",
    "conformal_flux_vanishes",
    "common_scalar_action_normalization_surface",
    "common_scalar_action_normalization_junction",
    "zero_field_surface",
    "zero_field_flux",
    "normal_and_side_reversal_curvature_jump",
    "normal_and_side_reversal_surface",
    "normal_and_side_reversal_flux",
    "amplitude_homogeneity_surface",
    "amplitude_homogeneity_flux",
    "length_rescaling_surface",
    "length_rescaling_flux",
)
CONTROL_NAMES = (
    "changed_scalar_jump_leaves_surface_residual",
    "nonconformal_nonminimal_surface_is_not_identically_zero",
    "nonconformal_nonminimal_flux_is_not_identically_zero",
)


def canonical_bytes(value: dict) -> bytes:
    return (json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True,
                       allow_nan=False) + "\n").encode("utf-8")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def read_relative(relative: str) -> bytes:
    path = ROOT / relative
    if (Path(relative).is_absolute() or ".." in Path(relative).parts
            or not path.resolve().is_relative_to(ROOT.resolve()) or not path.is_file()):
        raise ValueError("missing or unsafe prerequisite/source: " + relative)
    return path.read_bytes()


def source_bindings() -> dict:
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
            content = read_relative(directory + "/" + name)
            if sha(content) != binding["sha256"] or len(content) != binding["bytes"]:
                raise ValueError("prerequisite artifact binding mismatch: " + key + "/" + name)
            if name.endswith(".json"):
                payloads[name] = json.loads(content)
        if manifest["input_sha256"] != manifest["outputs"]["atlas.json"]["sha256"]:
            raise ValueError("prerequisite atlas binding mismatch: " + key)
        records[key] = payloads
        bindings[key] = {
            "manifest": manifest_path, "manifest_sha256": specification["manifest_sha256"],
            "source_commit": specification["source_commit"],
            "authenticated_outputs": deepcopy(manifest["outputs"]),
            "verification": "stored manifest and all output bytes; no source replay or numerical solve",
        }
    return records, bindings


def expected_certificate() -> dict:
    return {"sympy_version": SYMPY_VERSION,
            "zero_residuals": {name: "0" for name in IDENTITY_NAMES},
            "nonzero_controls": {name: True for name in CONTROL_NAMES},
            "arithmetic": "exact symbolic expressions; no floating-point scientific evaluation"}


def exact_algebra_checks() -> dict:
    """Check component and variational routes, not the imported GR theorem."""
    started = time.process_time()
    import sympy as sp
    if sp.__version__ != SYMPY_VERSION:
        raise RuntimeError("exact checks require pinned SymPy " + SYMPY_VERSION)
    phi, xi, k, d, dd, v, kappa, amplitude, scale = sp.symbols(
        "phi xi k d dd v kappa amplitude scale", nonzero=True)
    f_symbol, fp_symbol = sp.symbols("f fp")
    symbols = dict(zip(("phi", "xi", "k", "d", "dd", "v", "f", "fp"),
                       (phi, xi, k, d, dd, v, f_symbol, fp_symbol)))
    expressions = {key: sp.sympify(value, locals=symbols) for key, value in FORMULAS.items()}
    surface, flux = expressions["surface_coefficient"], expressions["flux_jump"]
    f = 1 / (2*kappa) - xi*phi**2/2
    normalization = sp.symbols("normalization", positive=True)
    normalized_f = 1/(2*kappa)-normalization*xi*phi**2/2
    # Independent variational route: derive f' from the fixed action.
    amm = expressions["amm_material_coefficient"].subs({f_symbol: f, fp_symbol: sp.diff(f, phi)})
    amm_residual = amm + 2*k/kappa  # sigma=2k/kappa, unchanged material tension.
    junction = {d: 6*xi*k*phi, dd: 6*xi*k*v}
    # Independent component route: Gsing_ab=-C_ab=-2k gamma_ab,
    # [normal(phi^2)]=2 phi d; scalar Hessian uses the product rule.
    component_surface = xi*(-2*k*phi**2 + 2*phi*d)
    qp, qm, dp, dm, kp, km = sp.symbols("qp qm dp dm kp km")
    component_flux_raw = (qp*v - xi*(2*qp*v + 2*phi*(dp-kp*v))) - (
        qm*v - xi*(2*qm*v + 2*phi*(dm-km*v)))
    component_flux = component_flux_raw.subs({qp: qm+d, dp: dm+dd, kp: km-k})
    surface_j = surface.subs(junction)
    flux_j = flux.subs(junction)
    conformal = {xi: sp.Rational(1, 6)}
    zero_field = {phi: 0, v: 0, d: 0, dd: 0}
    reversal = {qp: -qm, qm: -qp, dp: -dm, dm: -dp, kp: -km, km: -kp}
    surface_traces = surface.subs(d, qp-qm)
    scaling = {phi: amplitude*phi, v: amplitude*v, d: amplitude*d, dd: amplitude*dd}
    units = {k: k/scale, d: d/scale, v: v/scale, dd: dd/scale**2}
    residuals = {
        "improved_surface_matches_component_distribution": surface-component_surface,
        "surface_is_negative_amm_residual": surface+amm_residual,
        "surface_after_scalar_jump": surface_j-expressions["surface_after_scalar_jump"],
        "flux_matches_product_rule_components": flux-component_flux,
        "flux_after_scalar_jump": flux_j-expressions["flux_after_scalar_jump"],
        "general_coupling_surface_flux_balance": (
            expressions["surface_tangential_derivative"]+flux).subs(junction),
        "conformal_surface_vanishes": surface_j.subs(conformal),
        "conformal_amm_residual_vanishes": amm_residual.subs(junction).subs(conformal),
        "conformal_flux_vanishes": flux_j.subs(conformal),
        "common_scalar_action_normalization_surface": (
            -4*k*normalized_f+2*d*sp.diff(normalized_f, phi)+2*k/kappa+normalization*surface),
        "common_scalar_action_normalization_junction": (
            normalization*d-6*normalization*xi*k*phi).subs(junction),
        "zero_field_surface": surface.subs(zero_field),
        "zero_field_flux": flux.subs(zero_field),
        "normal_and_side_reversal_curvature_jump": (kp-km).xreplace(reversal)-(kp-km),
        "normal_and_side_reversal_surface": surface_traces.xreplace(reversal)-surface_traces,
        "normal_and_side_reversal_flux": component_flux_raw.xreplace(reversal)-component_flux_raw,
        "amplitude_homogeneity_surface": surface.xreplace(scaling)-amplitude**2*surface,
        "amplitude_homogeneity_flux": flux.xreplace(scaling)-amplitude**2*flux,
        "length_rescaling_surface": surface.xreplace(units)-surface/scale,
        "length_rescaling_flux": flux.xreplace(units)-flux/scale**2,
    }
    residuals = {name: sp.expand(value) for name, value in residuals.items()}
    if set(residuals) != set(IDENTITY_NAMES) or any(value != 0 for value in residuals.values()):
        raise ValueError("scalar/wall exact identity failed: " + str(residuals))
    controls = {
        "changed_scalar_jump_leaves_surface_residual": surface.subs(d, 0).subs(conformal),
        "nonconformal_nonminimal_surface_is_not_identically_zero": surface_j.subs(xi, sp.Rational(1, 4)),
        "nonconformal_nonminimal_flux_is_not_identically_zero": flux_j.subs(xi, sp.Rational(1, 4)),
    }
    if any(sp.expand(value) == 0 for value in controls.values()):
        raise ValueError("scalar/wall negative control unexpectedly vanished")
    if time.process_time() - started > 60:
        raise RuntimeError("exact algebra exceeded the admitted 60 CPU seconds")
    return expected_certificate()


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
    if any(response["domain"][key] != value for key, value in preparation["geometry"].items()):
        raise ValueError("prerequisite retained domains differ")
    if response["domain"]["arrival_window_over_ell"] != preparation["observer"]["arrival_window_over_ell"]:
        raise ValueError("prerequisite observer windows differ")
    if preparation["preparation"]["regularity"] != "C3 compact polynomial; derivative C2; not C-infinity":
        raise ValueError("the original C3 preparation is not retained")
    return {
        "schema": SCHEMA, "status": "CONDITIONAL_CONFORMAL_SCALAR_SURFACE_AND_TANGENTIAL_FLUX_COMPATIBILITY",
        "classification": "derived_conditional", "owner": OWNER,
        "inputs": deepcopy(preparation["inputs"]), "domain": deepcopy(response["domain"]),
        "preparation_regularity": preparation["preparation"]["regularity"],
        "borrowed": {"source": "1910.07534v2", "equations": ["2.10-2.14", "3.42-3.46"],
                     "use": "non-null variational scalar-tensor junction and field stress; fixed-background application",
                     "stored_evidence_use": "existing wall, conformal probe and regular finite preparation; no waveform derivatives used"},
        "conventions": {
            "signature": "(-+++)", "normal": "common spacelike normal from cap (-) to exterior (+)",
            "bracket": "plus trace minus minus trace using the same normal",
            "extrinsic_curvature": "K_ab=(1/2) L_n gamma_ab",
            "jump": "[K_ab]=-k gamma_ab; [K]=-3k; C_ab=[K_ab]-gamma_ab[K]=2k gamma_ab",
            "material_wall": "S_ab=-sigma gamma_ab; k=kappa_4 sigma/2; kappa_4=8 pi G",
            "action": "f(phi) R -(1/2)(grad phi)^2; f=1/(2 kappa_4)-xi phi^2/2; scalar potential zero",
            "split_action_boundary_terms": "the borrowed variational construction includes its 2 f K boundary terms",
            "stress_accounting": "xi G_mu_nu phi^2 included once in the improved scalar Hilbert stress on the right",
            "field_traces": "[phi]=0; q_pm=n phi_pm; d=[q]; v_a=D_a phi; dd=D_a d",
            "scalar_jump": "d=6 xi k phi; D_a d=6 xi k v_a for constant k",
            "orientation_check": "reverse the common normal and exchange plus/minus side labels together; one-sided normal derivatives and extrinsic curvatures negate and swap, their jumps and the flux jump stay unchanged, and positive k=kappa_4 sigma/2 stays fixed",
        },
        "assumptions": ["four dimensions; xi=1/6; no direct scalar-wall action or new surface source",
                        "zero scalar background; phi=epsilon phi_1; quadratic stress evaluated on the unchanged background",
                        "continuous scalar and induced metric, with sufficient piecewise regular traces for the stated products",
                        "original C3 finite preparation; full stress regularity is not inferred for the inverse's broader C2 class",
                        "both smooth bulk G_na=0; scalar-independent background vacuum stress",
                        "f>0 for weak-field interpretation; no quantitative amplitude bound is supplied"],
        "symbolic_identities": deepcopy(FORMULAS),
        "identity_interpretation": {
            "surface_coefficient": "Sphi_ab = surface_coefficient gamma_ab; Tsing_ab=Sphi_ab delta_Sigma, with proper-distance delta",
            "amm_material_coefficient": "coefficient required of material S_ab; compare with -sigma under k=kappa_4 sigma/2",
            "flux_jump": "[T_na] from one-sided regular traces, with [K_a^b]=-k delta_a^b and G_na=0",
            "general_coupling_balance": "D^a Sphi_ab + [T_nb]=0 after the borrowed scalar jump; algebraic diagnostic on the fixed background, not a new nonconformal solution",
            "normalization_control": "multiply the scalar kinetic and nonminimal terms by the same positive constant; surface stress scales by that constant and the scalar jump is unchanged; canonical coefficient one remains the chosen model",
            "length_rescaling_control": "hold phi fixed while rescaling lengths: k,d,v scale as 1/s and dd as 1/s^2; this checks derivative homogeneity, not a complete natural-unit rescaling or comparison of different physical couplings",
        },
        "exact_symbolic_certificate": expected_certificate(),
        "independent_review": {
            "scope": "independent distributional-stress and variational-junction derivations; common-normal tangential-flux calculation",
            "reference": OWNER, "limit": "exact algebra checks identities once their geometric and distributional assumptions are supplied"},
        "knowledge_boundary": {
            "additional_singular_scalar_surface_source_required_at_conformal_coupling": False,
            "nonzero_tangential_scalar_flux_jump_at_conformal_coupling": False,
            "normal_force_condition_computed": False, "bulk_einstein_response_solved": False,
            "fixed_background_with_nonzero_probe_claimed_exact": False,
            "quantitative_amplitude_bound_established": False, "stability_established": False,
            "physical_preparation_law_established": False,
        },
        "next_unlock": "a separately admitted preparation-energy or sourced-geometry question; keep normal-force and bulk response explicit",
        "execution": {"scientific_numerical_solves": 0, "new_waveforms": 0, "new_floating_point_scientific_values": False,
                      "exact_algebra_cpu_seconds_max": 60, "creation_runs_exact_algebra": True,
                      "ordinary_check_scope": "stored prerequisite, source and declaration authentication; use --check-algebra to repeat exact identities"},
    }, prerequisites


def create(directory: Path) -> dict:
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only evidence: use --check or a fresh directory")
    before = source_bindings()
    record, prerequisites = build_record()
    if exact_algebra_checks() != record["exact_symbolic_certificate"]:
        raise ValueError("exact symbolic certificate mismatch")
    payload = canonical_bytes(record)
    if source_bindings() != before:
        raise RuntimeError("source changed during evidence recording")
    manifest = {"schema": MANIFEST_SCHEMA, "source_bindings": before,
                "prerequisite_bindings": prerequisites,
                "outputs": {"consistency.json": {"sha256": sha(payload), "bytes": len(payload)}},
                "scope": "analytic shell compatibility; no numerical field replay or bulk Einstein solve"}
    directory.mkdir(parents=True, exist_ok=False)
    for name, content in (("consistency.json", payload), ("manifest.json", canonical_bytes(manifest))):
        with (directory / name).open("xb") as handle:
            handle.write(content)
    return {"created": directory.name, "outputs": 1, "scientific_numerical_solves": 0}


def check(directory: Path) -> dict:
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_bytes())
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("scalar/wall manifest schema mismatch")
    if manifest.get("source_bindings") != source_bindings():
        raise ValueError("scalar/wall source binding mismatch")
    record, prerequisites = build_record()
    if manifest.get("prerequisite_bindings") != prerequisites:
        raise ValueError("scalar/wall prerequisite binding mismatch")
    if (set(manifest.get("outputs", {})) != {"consistency.json"}
            or {path.name for path in directory.iterdir()} != {"consistency.json", "manifest.json"}):
        raise ValueError("scalar/wall output inventory mismatch")
    payload = (directory / "consistency.json").read_bytes()
    binding = manifest["outputs"]["consistency.json"]
    if sha(payload) != binding["sha256"] or len(payload) != binding["bytes"]:
        raise ValueError("scalar/wall artifact binding mismatch")
    if payload != canonical_bytes(record):
        raise ValueError("scalar/wall declaration mismatch")
    return {"checked": directory.name, "outputs": 1, "sources": len(manifest["source_bindings"]),
            "authenticated_predecessor_outputs": sum(len(row["authenticated_outputs"]) for row in prerequisites.values()),
            "scientific_numerical_solves": 0}


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--check-algebra", action="store_true")
    args = parser.parse_args(argv)
    result = exact_algebra_checks() if args.check_algebra else check(args.output) if args.check else create(args.output)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
