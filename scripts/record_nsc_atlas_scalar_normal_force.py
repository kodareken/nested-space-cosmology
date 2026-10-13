#!/usr/bin/env python3
"""Exact scalar normal-stress compatibility on the retained reference wall.

Ordinary --check authenticates declarations and evidence using only stdlib.
Explicit --check-algebra authenticates first and then repeats the new identities.
The original surface/tangential record and its dated knowledge flags stay intact.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-scalar-normal-force-v1"
MANIFEST_SCHEMA = SCHEMA+"-manifest"
OWNER = "docs/nsc-atlas-scalar-normal-force.md"
PRODUCER = "scripts/record_nsc_atlas_scalar_normal_force.py"
TESTS = "tests/test_nsc_atlas_scalar_normal_force.py"
RELATIVE_OUTPUT = "results/development/nsc-atlas-scalar-normal-force-v1"
DEFAULT = ROOT / RELATIVE_OUTPUT
PREDECESSOR_DIRECTORY = "results/development/nsc-atlas-scalar-wall-consistency-v1"
PREDECESSOR_MANIFEST_SHA = "61f4f2287b2f70196bf1a5fbc773bf529864514d38b46c80a44d9726a8d8b801"
PREDECESSOR_PRODUCER = "scripts/record_nsc_atlas_scalar_wall_consistency.py"
SYMPY_VERSION = "1.14.0"
CPU_SECONDS = 60
FORMULAS = {
    "normal_projection": "p**2/2-v2/2+xi*(phi**2*G+2*K*phi*p+lap)",
    "normal_jump": "pbar*d+xi*phi**2*deltaG+2*xi*phi*(Kbar*d+deltaK*pbar)",
    "after_scalar_and_geometric_jumps": "2*xi*(6*xi-1)*k*Kbar*phi**2",
    "gauss_jump": "-2*k*Kbar",
    "vacuum_gauss_jump": "3*H**2",
    "vacuum_mean_curvature": "-3*H**2/(2*k)",
}
IDENTITIES = (
    "gaussian_normal_hilbert_projection", "one_sided_jump_matches_product_rule",
    "gauss_constraint_from_tensor_and_trace_jumps", "general_coupling_reduction",
    "reference_beta_trace_jump", "reference_beta_mean_curvature",
    "general_jump_equals_original_surface_contraction", "conformal_normal_jump_vanishes",
    "israel_gauss_material_force_sign", "distributional_normal_divergence_sign",
    "vacuum_force_balance", "omitted_cap_einstein_term_false_residual",
    "normal_and_side_reversal_jump", "normal_and_side_reversal_scalar_jump",
    "normal_and_side_reversal_curvature_jump", "zero_field_jump",
    "quadratic_amplitude_homogeneity", "fixed_field_derivative_homogeneity",
    "fixed_vacuum_tension_geometric_remainder",
)
CONTROLS = (
    "wrong_material_force_sign_breaks_vacuum_balance",
    "omitting_cap_einstein_tensor_produces_false_force",
    "omitting_extrinsic_improvement_changes_jump",
    "altering_scalar_jump_changes_force",
    "nonconformal_diagnostic_is_not_identically_zero",
)


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
    return {path: sha(read_relative(path)) for path in sorted((OWNER, PRODUCER, TESTS))}


def authenticate_predecessor():
    """Authenticate original code bytes before reusing its existing checker."""
    relative = PREDECESSOR_DIRECTORY+"/manifest.json"
    raw = read_relative(relative)
    if sha(raw) != PREDECESSOR_MANIFEST_SHA:
        raise ValueError("pinned scalar-wall predecessor manifest mismatch")
    manifest = json.loads(raw)
    if (manifest.get("schema") != "nsc-atlas-scalar-wall-consistency-manifest-v1"
            or set(manifest.get("outputs", {})) != {"consistency.json"}
            or len(manifest.get("source_bindings", {})) != 3
            or PREDECESSOR_PRODUCER not in manifest["source_bindings"]):
        raise ValueError("scalar-wall predecessor schema/inventory mismatch")
    for source, expected in manifest["source_bindings"].items():
        if sha(read_relative(source)) != expected:
            raise ValueError("current scalar-wall predecessor source differs: "+source)
    payload = read_relative(PREDECESSOR_DIRECTORY+"/consistency.json")
    binding = manifest["outputs"]["consistency.json"]
    if sha(payload) != binding["sha256"] or len(payload) != binding["bytes"]:
        raise ValueError("scalar-wall predecessor artifact binding mismatch")
    spec = importlib.util.spec_from_file_location("authenticated_scalar_wall_predecessor", ROOT / PREDECESSOR_PRODUCER)
    predecessor = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(predecessor)
    verified = predecessor.check(ROOT / PREDECESSOR_DIRECTORY)
    if verified["authenticated_predecessor_outputs"] != 15 or verified["sources"] != 3:
        raise ValueError("unexpected scalar-wall prerequisite chain")
    record = json.loads(payload)
    return record, {
        "manifest": relative, "manifest_sha256": PREDECESSOR_MANIFEST_SHA,
        "authenticated_outputs": deepcopy(manifest["outputs"]),
        "authenticated_sources": deepcopy(manifest["source_bindings"]),
        "original_prerequisite_bindings": deepcopy(manifest["prerequisite_bindings"]),
        "authenticated_predecessor_outputs": 16,
        "verification": "original manifest/output/source bytes checked before invoking its stdlib checker; all fifteen v7/v8 outputs authenticated; no algebra or field replay",
    }


def expected_certificate():
    return {"sympy_version": SYMPY_VERSION, "zero_residuals": dict.fromkeys(IDENTITIES, "0"),
            "nonzero_controls": dict.fromkeys(CONTROLS, True),
            "arithmetic": "exact symbolic expressions; no floating-point scientific evaluation"}


def exact_algebra_checks(original_surface_formula, started=None):
    started = time.process_time() if started is None else started
    import sympy as sp
    if sp.__version__ != SYMPY_VERSION:
        raise RuntimeError("exact checks require pinned SymPy "+SYMPY_VERSION)
    names = "phi xi k d pbar Kbar deltaK deltaG p K G v2 lap H kappa scale amplitude"
    symbols = dict(zip(names.split(), sp.symbols(names)))
    phi, xi, k, d, pbar, Kbar, deltaK, deltaG, p, K, G, v2, lap, H, kappa, scale, amplitude = symbols.values()
    expressions = {key: sp.sympify(value, locals=symbols) for key, value in FORMULAS.items()}
    projection, jump = expressions["normal_projection"], expressions["normal_jump"]
    nn = sp.symbols("normal_second_phi_squared")
    box = nn+2*K*phi*p+lap
    hilbert_projection = p*p-(p*p+v2)/2+xi*(phi*phi*G+box-nn)
    pp, pm, Kp, Km, Gp, Gm = sp.symbols("pplus pminus Kplus Kminus Gplus Gminus")
    one_sided = projection.subs({p: pp, K: Kp, G: Gp})-projection.subs({p: pm, K: Km, G: Gm})
    means = {pp: pbar+d/2, pm: pbar-d/2, Kp: Kbar+deltaK/2, Km: Kbar-deltaK/2,
             Gp: (Gp+Gm)/2+deltaG/2, Gm: (Gp+Gm)/2-deltaG/2}
    component_jump = one_sided.subs(means, simultaneous=True)
    jumps = {d: 6*xi*k*phi, deltaK: -3*k, deltaG: -2*k*Kbar}
    reduced = jump.subs(jumps)
    original_surface = sp.sympify(original_surface_formula, locals=symbols)
    conformal = {xi: sp.Rational(1, 6)}
    sigma = 2*k/kappa
    contraction = -sigma*Kbar
    gauss = (2*Kbar*(-3*k)-2*(-k)*Kbar)/2
    vacuum_mean = expressions["vacuum_mean_curvature"]
    vacuum = {Kbar: vacuum_mean, deltaG: 3*H*H, deltaK: -3*k, d: k*phi}
    no_einstein = (jump-xi*phi*phi*deltaG).subs(vacuum).subs(conformal)
    reversal = {pp: -pm, pm: -pp, Kp: -Km, Km: -Kp, Gp: Gm, Gm: Gp}
    scaled_field = {phi: amplitude*phi, pbar: amplitude*pbar, d: amplitude*d}
    scaled_derivatives = {pbar: pbar/scale, d: d/scale, Kbar: Kbar/scale,
                          deltaK: deltaK/scale, deltaG: deltaG/scale**2}
    variation = sp.symbols("delta_mean_curvature")
    radius, mass = sp.symbols("radius mass", nonzero=True)
    beta_minus = mass/(k*radius**2)+(k*k-H*H)*radius/(2*k)
    beta_plus = mass/(k*radius**2)-(k*k+H*H)*radius/(2*k)
    trace_minus = sp.diff(beta_minus, radius)+2*beta_minus/radius
    trace_plus = sp.diff(beta_plus, radius)+2*beta_plus/radius
    residuals = {
        "gaussian_normal_hilbert_projection": projection-hilbert_projection,
        "one_sided_jump_matches_product_rule": jump-component_jump,
        "gauss_constraint_from_tensor_and_trace_jumps": gauss-expressions["gauss_jump"],
        "reference_beta_trace_jump": trace_plus-trace_minus+3*k,
        "reference_beta_mean_curvature": (trace_plus+trace_minus)/2-vacuum_mean,
        "general_coupling_reduction": reduced-expressions["after_scalar_and_geometric_jumps"],
        "general_jump_equals_original_surface_contraction": reduced-original_surface*Kbar,
        "conformal_normal_jump_vanishes": reduced.subs(conformal),
        "israel_gauss_material_force_sign": contraction-gauss/kappa,
        "distributional_normal_divergence_sign": gauss/kappa-contraction,
        "vacuum_force_balance": contraction.subs(Kbar, vacuum_mean)-expressions["vacuum_gauss_jump"]/kappa,
        "omitted_cap_einstein_term_false_residual": no_einstein+H*H*phi*phi/2,
        "normal_and_side_reversal_jump": one_sided.xreplace(reversal)+one_sided,
        "normal_and_side_reversal_scalar_jump": (pp-pm).xreplace(reversal)-(pp-pm),
        "normal_and_side_reversal_curvature_jump": (Kp-Km).xreplace(reversal)-(Kp-Km),
        "zero_field_jump": jump.subs({phi: 0, pbar: 0, d: 0}),
        "quadratic_amplitude_homogeneity": jump.xreplace(scaled_field)-amplitude**2*jump,
        "fixed_field_derivative_homogeneity": jump.xreplace(scaled_derivatives)-jump/scale**2,
        "fixed_vacuum_tension_geometric_remainder":
            -sigma*(Kbar+variation)+sigma*Kbar+sigma*variation-reduced.subs(conformal),
    }
    simplified = {name: sp.simplify(value) for name, value in residuals.items()}
    if set(simplified) != set(IDENTITIES) or any(value != 0 for value in simplified.values()):
        raise ValueError("scalar normal-force exact identity failed: "+str(simplified))
    controls = {
        "wrong_material_force_sign_breaks_vacuum_balance": contraction.subs(Kbar, vacuum_mean)+3*H*H/kappa,
        "omitting_cap_einstein_tensor_produces_false_force": no_einstein,
        "omitting_extrinsic_improvement_changes_jump":
            (jump-2*xi*phi*(Kbar*d+deltaK*pbar)).subs(vacuum).subs(conformal),
        "altering_scalar_jump_changes_force": jump.subs({d: 0, deltaK: -3*k, deltaG: -2*k*Kbar}).subs(conformal),
        "nonconformal_diagnostic_is_not_identically_zero": reduced.subs(xi, sp.Rational(1, 4)),
    }
    if set(controls) != set(CONTROLS) or any(sp.simplify(value) == 0 for value in controls.values()):
        raise ValueError("scalar normal-force negative control unexpectedly vanished")
    if time.process_time()-started > CPU_SECONDS:
        raise RuntimeError("normal-force creation exceeded admitted process CPU budget")
    return expected_certificate()


def build_record():
    predecessor, bindings = authenticate_predecessor()
    if (predecessor["status"] != "CONDITIONAL_CONFORMAL_SCALAR_SURFACE_AND_TANGENTIAL_FLUX_COMPATIBILITY"
            or predecessor["inputs"]["xi"] != 1/6
            or predecessor["preparation_regularity"] != "C3 compact polynomial; derivative C2; not C-infinity"
            or predecessor["knowledge_boundary"]["normal_force_condition_computed"] is not False):
        raise ValueError("original scalar-wall model or dated knowledge boundary differs")
    return {
        "schema": SCHEMA, "status": "CONDITIONAL_CONFORMAL_SCALAR_NORMAL_STRESS_COMPATIBILITY",
        "classification": "derived_conditional", "owner": OWNER,
        "inputs": deepcopy(predecessor["inputs"]), "domain": deepcopy(predecessor["domain"]),
        "preparation_regularity": predecessor["preparation_regularity"],
        "borrowed": {"hilbert_stress_and_scalar_material_jumps": "1910.07534v2 Eqs. 2.11-2.14, 3.42-3.46",
                     "normal_gauss_and_shell_context": "gr-qc/0507063v2 Eqs. 5-6; printed evolution sign is checked against our conventions instead of copied",
                     "previous_surface_formula": predecessor["symbolic_identities"]["surface_after_scalar_jump"]},
        "conventions": {"signature": "(-+++)", "normal": "common spacelike normal from cap (-) to exterior (+)",
                        "bracket": "plus trace minus minus trace; bars are arithmetic means",
                        "curvature": "K_ab=(1/2)L_n gamma_ab; [K_ab]=-k gamma_ab; [K]=-3k",
                        "scalar": "p_pm=n phi_pm; d=[p]=6 xi k phi, hence k phi at xi=1/6",
                        "stress": "full canonical positive-kinetic improved scalar stress; xi G_nn phi^2 retained exactly once",
                        "force_sign": "S^ab mean(K_ab)=[T_nn]; bulk divergence contributes +[T_nn] delta, surface contributes -S^ab mean(K_ab) delta",
                        "vacuum": "Gplus_nn=0; Gminus_nn=-3H^2; [G_nn]=3H^2; mean(K)=-3H^2/(2k)",
                        "material": "S_ab=-sigma gamma_ab; sigma=2k/kappa_4; kappa_4=8 pi G",
                        "order": "quadratic scalar forcing on unchanged reference geometry, not an exact sourced spacetime"},
        "formulas": deepcopy(FORMULAS),
        "derived_result": {
            "general_diagnostic": "[Tphi_nn]=2 xi(6xi-1) k mean(K) phi^2 = mean(K) Sphi_coefficient; symbolic coupling diagnostic, not an evolved nonconformal model",
            "physical_model": "[Tphi_nn]=0 at xi=1/6; individual one-sided normal stresses are not thereby zero",
            "full_balance": "-sigma mean(K)=[Tvac_nn]+[Tphi_nn]",
            "remaining_geometric_condition": "fixed vacuum and tension imply delta mean(K)=0 at quadratic order; this combines metric and wall-embedding corrections and does not imply delta R_wall=0",
            "cap_omission_control": "omitting the cap Einstein term falsely gives -H^2 phi^2/2",
            "intrinsic_terms": "the shared -v_a v^a/2 and xi D^2(phi^2) cancel in the jump under the declared trace regularity",
        },
        "orientation_control": "reverse normal and swap sides together: pbar, Kbar, deltaG and normal-stress jump negate; d, deltaK and positive k remain unchanged",
        "scaling_control": "amplitude scales phi,pbar,d together; derivative homogeneity holds at fixed phi when first normal derivatives/curvatures scale as inverse length and G_nn as inverse length squared",
        "exact_symbolic_certificate": expected_certificate(),
        "knowledge_boundary": {"scalar_normal_stress_jump_computed": True,
            "coupled_geometric_normal_force_equation_solved": False, "bulk_einstein_response_solved": False,
            "corrected_wall_trajectory_solved": False, "fixed_background_with_nonzero_probe_claimed_exact": False,
            "one_sided_pressure_curves_evaluated": False, "finite_amplitude_validity_bound_established": False,
            "stability_established": False},
        "execution": {"scientific_numerical_solves": 0, "new_waveforms": 0,
                      "new_floating_point_scientific_values": False, "cpu_seconds_max": CPU_SECONDS,
                      "ordinary_check_scope": "stored output, source and declaration authentication only",
                      "explicit_algebra_replay_authenticates_first": True},
        "next_unlock": "transport sourced geometry and test its coupled wall equations under a separate admission; no pressure simulation is needed to repeat this identity",
    }, bindings


def create(directory):
    started = time.process_time()
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only evidence: use --check or a fresh directory")
    before = source_bindings()
    record, bindings = build_record()
    if exact_algebra_checks(record["borrowed"]["previous_surface_formula"], started) != record["exact_symbolic_certificate"]:
        raise ValueError("exact normal-force certificate mismatch")
    if source_bindings() != before or authenticate_predecessor()[1] != bindings:
        raise RuntimeError("source or prerequisite changed during evidence recording")
    payload = canonical_bytes(record)
    manifest = {"schema": MANIFEST_SCHEMA, "source_bindings": before,
                "predecessor_binding": bindings,
                "outputs": {"normal-force.json": {"sha256": sha(payload), "bytes": len(payload)}},
                "scope": "exact local scalar normal-stress jump; coupled metric and wall corrections remain unsolved"}
    if time.process_time()-started > CPU_SECONDS:
        raise RuntimeError("normal-force creation exceeded admitted process CPU budget")
    directory.mkdir(parents=True, exist_ok=False)
    for name, content in (("normal-force.json", payload), ("manifest.json", canonical_bytes(manifest))):
        with (directory / name).open("xb") as handle:
            handle.write(content)
    return {"created": directory.name, "outputs": 1, "scientific_numerical_solves": 0, "cpu_guard_passed": True}


def check(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_bytes())
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("normal-force manifest schema mismatch")
    if manifest.get("source_bindings") != source_bindings():
        raise ValueError("normal-force source binding mismatch")
    record, bindings = build_record()
    if manifest.get("predecessor_binding") != bindings:
        raise ValueError("normal-force predecessor binding mismatch")
    if (set(manifest.get("outputs", {})) != {"normal-force.json"}
            or {p.name for p in directory.iterdir()} != {"normal-force.json", "manifest.json"}):
        raise ValueError("normal-force output inventory mismatch")
    payload = (directory / "normal-force.json").read_bytes()
    binding = manifest["outputs"]["normal-force.json"]
    if sha(payload) != binding["sha256"] or len(payload) != binding["bytes"]:
        raise ValueError("normal-force artifact binding mismatch")
    if payload != canonical_bytes(record):
        raise ValueError("normal-force declaration mismatch")
    return {"checked": directory.name, "outputs": 1, "sources": 3,
            "authenticated_predecessor_outputs": 16, "scientific_numerical_solves": 0}


def check_algebra(directory):
    started = time.process_time()
    check(directory)
    record = json.loads((Path(directory) / "normal-force.json").read_bytes())
    result = exact_algebra_checks(record["borrowed"]["previous_surface_formula"], started)
    if result != record["exact_symbolic_certificate"]:
        raise ValueError("normal-force algebra replay differs from authenticated record")
    return {"replayed": Path(directory).name, "exact_zero_residuals": len(IDENTITIES),
            "nonzero_controls": len(CONTROLS), "cpu_guard_passed": True}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--check-algebra", action="store_true")
    args = parser.parse_args(argv)
    result = check_algebra(args.output) if args.check_algebra else check(args.output) if args.check else create(args.output)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
