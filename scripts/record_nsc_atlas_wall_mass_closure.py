#!/usr/bin/env python3
"""Exact exterior-side wall mass closure; no numerical mass evolution.

Ordinary --check authenticates bytes and declarations using only stdlib.
Explicit --check-algebra authenticates before repeating the exact identities.
"""
from __future__ import annotations

import argparse
from copy import deepcopy
import hashlib
import importlib.util
import json
from pathlib import Path
import time

from validation_source_successors import verify_saved_bindings, verify_saved_source

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-wall-mass-closure-v1"
MANIFEST_SCHEMA = SCHEMA+"-manifest"
OWNER = "docs/nsc-atlas-wall-mass-closure.md"
PRODUCER = "scripts/record_nsc_atlas_wall_mass_closure.py"
TESTS = "tests/test_nsc_atlas_wall_mass_closure.py"
SOURCE_PATHS = (OWNER, PRODUCER, TESTS)
RELATIVE_OUTPUT = "results/development/nsc-atlas-wall-mass-closure-v1"
DEFAULT = ROOT / RELATIVE_OUTPUT
SYMPY_VERSION = "1.14.0"
CPU_SECONDS = 60
PREDECESSORS = {
    "source_mass": {
        "directory": "results/development/nsc-atlas-source-mass-v1",
        "manifest_sha256": "2ab5a389669aebb1414633691f4416ecb5bdbabae47aa35823ecd5720e7230b5",
        "schema": "nsc-atlas-source-mass-v1-manifest", "output": "mass.json", "sources": 4,
        "producer": "scripts/record_nsc_atlas_source_mass.py", "chain_outputs": 17,
        "role": "mass normalization, quiet characteristic anchor and authenticated leading field context",
    },
    "normal_force": {
        "directory": "results/development/nsc-atlas-scalar-normal-force-v1",
        "manifest_sha256": "c15a3744aaf7d1b1eeba3636efa5a46bba9b7e582db9c1f1186f392ec225445a",
        "schema": "nsc-atlas-scalar-normal-force-v1-manifest", "output": "normal-force.json", "sources": 3,
        "producer": "scripts/record_nsc_atlas_scalar_normal_force.py", "chain_outputs": 16,
        "role": "boundary-context status only; normal-force identity is not an algebraic premise of exterior mass closure",
    },
}
FORMULAS = {
    "tangent_field": "(j-v*p)/x",
    "normal_field": "(j-beta*p)/x",
    "minimal_mass_density": "x**2*(beta*a*b-v*(a**2+b**2)/2)/2",
    "frame_endpoint": "m*p**2/12+x**2*p*(v*a-beta*b)/6",
    "frame_endpoint_from_trace": "(v-beta)*w*j/6+(1/(6*x)-m/(4*x**2))*w**2",
    "improved_t_un": "2*a*b/3-p*Hun/3",
    "improved_t_nn": "a**2/6+b**2/2-p*Hnn/3",
    "wave_hessian": "Huu+2*(v*a-beta*b)/x",
    "null_mass_density": "-x**2*(xV*pU**2*Us+xU*pV**2*Vs)/C",
}
IDENTITIES = (
    "tangent_trace_from_reduced_field", "normal_trace_from_no_return_condition",
    "no_return_retains_transverse_scalar_derivative", "endpoint_from_first_order_wall_trace",
    "both_null_fluxes_equal_minimal_wall_density", "full_conformal_mass_endpoint",
    "trapped_sector_minimal_density_decomposition",
    "jordan_improved_stress_equals_integral_plus_endpoint_derivative",
    "wall_acceleration_cancels", "zero_field_density", "quiet_anchor_endpoint",
    "quadratic_density_homogeneity", "quadratic_endpoint_homogeneity",
    "independent_null_coordinate_reparametrizations", "consistent_normal_reversal_density",
    "consistent_normal_reversal_endpoint", "beta_zero_density_regular",
    "beta_zero_endpoint_regular", "physical_mass_unit_rescaling",
    "fixed_coupling_length_change_is_not_unit_rescaling",
    "constant_background_mass_makes_cut_displacement_higher_order",
)
CONTROLS = (
    "omitting_first_null_flux_changes_density", "omitting_second_null_flux_changes_density",
    "omitting_transverse_frame_term_changes_mass", "omitting_frame_endpoint_breaks_jordan_transport",
    "reversing_field_normal_without_geometry_changes_density",
    "zeroing_phi_U_after_chi_U_zero_omits_nonzero_flux",
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
    return {path: sha(read_relative(path)) for path in sorted(SOURCE_PATHS)}


def _require_saved_source(source, expected, label):
    try:
        verify_saved_source(source, expected, ROOT)
    except FileNotFoundError as error:
        raise ValueError("missing or unsafe prerequisite/source: "+source) from error
    except ValueError as error:
        raise ValueError(label) from error


def authenticate_predecessors():
    """Check both immediate source sets before loading either existing checker."""
    records, bindings, pending = {}, {}, []
    for key, specification in PREDECESSORS.items():
        relative = specification["directory"]+"/manifest.json"
        raw = read_relative(relative)
        if sha(raw) != specification["manifest_sha256"]:
            raise ValueError("pinned predecessor manifest mismatch: "+key)
        manifest = json.loads(raw)
        if (manifest.get("schema") != specification["schema"]
                or set(manifest.get("outputs", {})) != {specification["output"]}
                or len(manifest.get("source_bindings", {})) != specification["sources"]
                or specification["producer"] not in manifest["source_bindings"]):
            raise ValueError("predecessor schema/inventory mismatch: "+key)
        for source, expected in manifest["source_bindings"].items():
            _require_saved_source(source, expected, "current predecessor source differs: "+source)
        payload = read_relative(specification["directory"]+"/"+specification["output"])
        binding = manifest["outputs"][specification["output"]]
        if sha(payload) != binding["sha256"] or len(payload) != binding["bytes"]:
            raise ValueError("predecessor artifact binding mismatch: "+key)
        records[key] = json.loads(payload)
        bindings[key] = {
            "manifest": relative, "manifest_sha256": sha(raw), "role": specification["role"],
            "authenticated_sources": deepcopy(manifest["source_bindings"]),
            "authenticated_outputs": deepcopy(manifest["outputs"]),
            "original_chain_bindings": deepcopy(manifest.get("prerequisite_bindings", manifest.get("predecessor_binding"))),
            "authenticated_chain_outputs_including_immediate": specification["chain_outputs"]+1,
            "verification": "pinned manifest, output and current bound sources before stdlib checker; no algebra, waveform or quadrature replay",
        }
        pending.append((key, specification))
    for key, specification in pending:
        spec = importlib.util.spec_from_file_location("authenticated_wall_mass_"+key, ROOT / specification["producer"])
        predecessor = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(predecessor)
        verified = predecessor.check(ROOT / specification["directory"])
        if (verified["sources"] != specification["sources"]
                or verified["authenticated_predecessor_outputs"] != specification["chain_outputs"]):
            raise ValueError("unexpected predecessor chain: "+key)
    return records, bindings


def expected_certificate():
    return {"sympy_version": SYMPY_VERSION, "zero_residuals": dict.fromkeys(IDENTITIES, "0"),
            "nonzero_controls": dict.fromkeys(CONTROLS, True),
            "arithmetic": "exact symbolic expressions; no numerical mass curve or floating-point scientific evaluation"}


def exact_algebra_checks(started=None):
    started = time.process_time() if started is None else started
    import sympy as sp
    if sp.__version__ != SYMPY_VERSION:
        raise RuntimeError("exact checks require pinned SymPy "+SYMPY_VERSION)
    names = "x m p a b v beta j w acc Huu Hun Hnn C xU xV pU pV Us Vs amplitude ru rv ell kappa epsilon mu scale alpha"
    symbols = dict(zip(names.split(), sp.symbols(names)))
    (x, m, p, a, b, v, beta, j, w, acc, Huu, Hun, Hnn, C, xU, xV,
     pU, pV, Us, Vs, amplitude, ru, rv, ell, kappa, epsilon, mu, scale, alpha) = symbols.values()
    expr = {key: sp.sympify(value, locals=symbols) for key, value in FORMULAS.items()}
    L, Q = expr["minimal_mass_density"], expr["frame_endpoint"]
    # Differentiate along the timelike reference before imposing its norm relation.
    derivatives = {x: v, p: a, a: Huu+acc*b, b: Hun+acc*a,
                   v: acc*beta-m/x**2, beta: acc*v}
    Qs = sum(sp.diff(Q, variable)*rate for variable, rate in derivatives.items())
    jordan = x*x*(beta*expr["improved_t_un"]-v*expr["improved_t_nn"])/2
    jordan = jordan.subs(Hnn, expr["wave_hessian"])
    null = expr["null_mass_density"]
    null_legs = {xU: (v-beta)/(2*Us), xV: (v+beta)/(2*Vs),
                 pU: (a-b)/(2*Us), pV: (a+b)/(2*Vs)}
    null_on_wall = null.subs(null_legs).subs(C, 1/(Us*Vs))
    omega = 1-alpha*p*p/12
    omega_u, omega_n = -alpha*p*a/6, -alpha*p*b/6
    # Full smooth-exterior mass transform, including both contractions. Its
    # gradient-squared term starts at alpha^2 but is kept before expansion.
    transformed_mass = (omega*(m+alpha*mu)-x*x*(-omega_u*v+omega_n*beta)
                        -x**3*(-omega_u**2+omega_n**2)/(2*omega))
    field_scale = {p: amplitude*p, a: amplitude*a, b: amplitude*b}
    coordinate_change = {C: C/(ru*rv), xU: xU/ru, xV: xV/rv,
                         pU: pU/ru, pV: pV/rv, Us: ru*Us, Vs: rv*Vs}
    physical_correction = kappa*epsilon**2*mu/ell
    wall_traces = {a: expr["tangent_field"], b: expr["normal_field"]}
    cut, displacement = sp.symbols("cut displacement")
    coefficient = sp.Function("mass_coefficient")
    shifted_mass = m+alpha*coefficient(cut+alpha*displacement)
    residuals = {
        "tangent_trace_from_reduced_field": expr["tangent_field"]-(j/x-w*v/x**2).subs(w, x*p),
        "normal_trace_from_no_return_condition": expr["normal_field"]-(j/x-w*beta/x**2).subs(w, x*p),
        "no_return_retains_transverse_scalar_derivative": (a-b+p*(v-beta)/x).subs(wall_traces),
        "endpoint_from_first_order_wall_trace": Q.subs(wall_traces).subs(p, w/x)-expr["frame_endpoint_from_trace"],
        "both_null_fluxes_equal_minimal_wall_density": null_on_wall-L,
        "trapped_sector_minimal_density_decomposition": L+x*x*((v+beta)*(a-b)**2+(v-beta)*(a+b)**2)/8,
        "full_conformal_mass_endpoint": sp.diff(transformed_mass, alpha).subs(alpha, 0)-(mu-Q),
        "jordan_improved_stress_equals_integral_plus_endpoint_derivative": L+Qs-jordan,
        "wall_acceleration_cancels": sp.diff(Qs, acc),
        "zero_field_density": L.subs({p: 0, a: 0, b: 0}),
        "quiet_anchor_endpoint": Q.subs({p: 0, a: 0, b: 0}),
        "quadratic_density_homogeneity": L.xreplace(field_scale)-amplitude**2*L,
        "quadratic_endpoint_homogeneity": Q.xreplace(field_scale)-amplitude**2*Q,
        "independent_null_coordinate_reparametrizations": null.xreplace(coordinate_change)-null,
        "consistent_normal_reversal_density": L.xreplace({beta: -beta, b: -b})-L,
        "consistent_normal_reversal_endpoint": Q.xreplace({beta: -beta, b: -b})-Q,
        "beta_zero_density_regular": L.subs(beta, 0)+x*x*v*(a*a+b*b)/4,
        "beta_zero_endpoint_regular": Q.subs(beta, 0)-m*p*p/12-x*x*p*v*a/6,
        "physical_mass_unit_rescaling": physical_correction.xreplace({ell: scale*ell, kappa: scale**2*kappa})-scale*physical_correction,
        "fixed_coupling_length_change_is_not_unit_rescaling": physical_correction.subs(ell, scale*ell)-physical_correction/scale,
        "constant_background_mass_makes_cut_displacement_higher_order": sp.diff(shifted_mass, alpha).subs(alpha, 0)-coefficient(cut),
    }

    def reduce(expression):
        # Equivalent to beta^2-v^2=1-2m/x; applied after differentiation.
        return sp.factor(expression.subs(m, x*(1-beta*beta+v*v)/2))

    simplified = {name: reduce(value) for name, value in residuals.items()}
    if set(simplified) != set(IDENTITIES) or any(value != 0 for value in simplified.values()):
        raise ValueError("wall mass exact identity failed: "+str(simplified))
    first = -x*x*xV*pU*pU*Us/C
    second = -x*x*xU*pV*pV*Vs/C
    controls = {
        "omitting_first_null_flux_changes_density": second.subs(null_legs).subs(C, 1/(Us*Vs))-L,
        "omitting_second_null_flux_changes_density": first.subs(null_legs).subs(C, 1/(Us*Vs))-L,
        "omitting_transverse_frame_term_changes_mass": x*x*p*beta*b/6,
        "omitting_frame_endpoint_breaks_jordan_transport": L-jordan,
        "reversing_field_normal_without_geometry_changes_density": L.subs(b, -b)-L,
        "zeroing_phi_U_after_chi_U_zero_omits_nonzero_flux": first.subs(null_legs).subs(C, 1/(Us*Vs)).subs(wall_traces),
    }
    if set(controls) != set(CONTROLS) or any(reduce(value) == 0 for value in controls.values()):
        raise ValueError("wall mass negative control unexpectedly vanished")
    if time.process_time()-started > CPU_SECONDS:
        raise RuntimeError("wall mass creation exceeded admitted process CPU budget")
    return expected_certificate()


def build_record():
    records, bindings = authenticate_predecessors()
    mass, force = records["source_mass"], records["normal_force"]
    if (mass["status"] != "QUADRATIC_SOURCE_LEG_INVARIANT_MASS_RESPONSE"
            or force["status"] != "CONDITIONAL_CONFORMAL_SCALAR_NORMAL_STRESS_COMPATIBILITY"
            or mass["inputs"] != force["inputs"] or mass["inputs"]["xi"] != 1/6
            or force["preparation_regularity"] != "C3 compact polynomial; derivative C2; not C-infinity"
            or mass["knowledge_boundary"]["corrected_wall_solved"] is not False
            or force["knowledge_boundary"]["corrected_wall_trajectory_solved"] is not False):
        raise ValueError("predecessor model or dated knowledge boundary differs")
    return {
        "schema": SCHEMA, "status": "CONDITIONAL_EXTERIOR_WALL_MASS_FUNCTIONAL_CLOSED",
        "classification": "derived_conditional", "owner": OWNER,
        "inputs": deepcopy(mass["inputs"]), "domain": deepcopy(force["domain"]),
        "preparation_regularity": force["preparation_regularity"],
        "borrowed": {
            "both_null_mass_laws": "An--Lim, 2005.04090v1 Eqs. 2.14-2.15; canonical source-mass normalization retained",
            "smooth_exterior_frame": "Aviles--Maeda--Martinez, 1910.07534v2 Eqs. 2.21-2.24; not continued through the material wall",
            "full_geometric_mass_transform": "Faraoni--Vitagliano, 1401.1189v1 Eq. 7; both derivative contractions retained",
            "improved_jordan_transport": "Hayward, gr-qc/9710089v2 Eq. 2.9d; canonical improved conformal stress",
            "normal_force_record_role": PREDECESSORS["normal_force"]["role"],
        },
        "conventions": {
            "region": "one-sided smooth exterior limit of the original wall, x=R/ell in [0.45,0.55], retained sector IV",
            "clock": "s=tau_wall/ell; v=dx/ds>0; beta=beta_plus; beta^2-v^2=1-2m/x; m=2/5",
            "normal": "outward cap-to-exterior normal (-Udot,Vdot); beta=n R on the reference",
            "traces": "w=chi_wall/epsilon; j=dw/ds=ell Vdot J/epsilon; p=w/x; a=dp/ds; b=n_dimensionless p",
            "no_return": "chi_U=0 on the retained wall, hence n chi=dot chi; no second independently prescribed characteristic",
            "mass": "mathfrak_m_J/ell=m+alpha mu_J+O(alpha^2); alpha=kappa_4 epsilon^2/ell^2; formal quadratic order only",
            "anchor": "quiet lower wall event with zero scalar data and zero geometric mass correction; no imported source-end value",
        },
        "formulas": deepcopy(FORMULAS),
        "derived_result": {
            "closure": "mu_J(s)=integral from s_a to s of L(q) dq+Q(s), with Q(s_a)=0",
            "inputs_required": "existing leading wall traces w,j and background x,v,beta, clock; no unknown Hessian or perturbed metric enters this coefficient",
            "independent_route": "L+dQ/ds=x^2(beta t_un-v t_nn)/2 with full Jordan improved stress; acceleration cancels",
            "null_route": "both Einstein-frame null fluxes along the same timelike path give L; add the full smooth-exterior frame endpoint Q",
            "regularity": "L and Q have no division by beta; their beta=0 values remain finite at x>0",
            "minimal_sign": "v>|beta| in the retained trapped sector makes L nonpositive; no late Jordan-mass sign is inferred",
            "cut_identification": "constant background exterior mass makes an order-alpha cut displacement invisible to this coefficient; no vanishing wall displacement is inferred",
        },
        "scaling": {
            "amplitude": "rescaling all leading field traces together multiplies L and Q quadratically; varying epsilon alone is represented by alpha",
            "units": "physical mass-length correction ell alpha mu_J; unit relabelling scales ell and kappa_4 as length and length squared, respectively",
            "coordinate": "independent regular null-coordinate changes leave both-flux density invariant; they do not alter wall proper time",
        },
        "exact_symbolic_certificate": expected_certificate(),
        "knowledge_boundary": {
            "exterior_wall_mass_functional_closed": True, "numerical_wall_mass_curve_evaluated": False,
            "cap_mass_determined": False, "material_shell_mass_determined": False,
            "bulk_einstein_evolution_solved": False, "corrected_wall_trajectory_solved": False,
            "interior_observer_correction_computed": False, "finite_amplitude_validity_bound_established": False,
            "stability_established": False, "full_coupled_solution_established": False,
        },
        "execution": {"scientific_numerical_solves": 0, "new_waveforms": 0, "numerical_mass_integrals": 0,
                      "new_floating_point_scientific_values": False, "cpu_seconds_max": CPU_SECONDS,
                      "ordinary_check_scope": "stored output, source and declaration authentication only",
                      "explicit_algebra_replay_authenticates_first": True},
        "next_unlock": "separately admitted finite mass quadrature and exterior input to coupled wall/cap matching; v8 waveform convergence does not certify this integral",
    }, bindings


def create(directory):
    started = time.process_time()
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only evidence: use --check or a fresh directory")
    before = source_bindings()
    record, bindings = build_record()
    if exact_algebra_checks(started) != record["exact_symbolic_certificate"]:
        raise ValueError("exact wall mass certificate mismatch")
    if source_bindings() != before or authenticate_predecessors()[1] != bindings:
        raise RuntimeError("source or prerequisite changed during evidence recording")
    payload = canonical_bytes(record)
    manifest = {"schema": MANIFEST_SCHEMA, "source_bindings": before, "predecessor_bindings": bindings,
                "outputs": {"wall-mass.json": {"sha256": sha(payload), "bytes": len(payload)}},
                "scope": "exact exterior-side quadratic mass functional; no numerical evaluation or coupled wall solution"}
    if time.process_time()-started > CPU_SECONDS:
        raise RuntimeError("wall mass creation exceeded admitted process CPU budget")
    directory.mkdir(parents=True, exist_ok=False)
    for name, content in (("wall-mass.json", payload), ("manifest.json", canonical_bytes(manifest))):
        with (directory / name).open("xb") as handle:
            handle.write(content)
    return {"created": directory.name, "outputs": 1, "scientific_numerical_solves": 0, "cpu_guard_passed": True}


def check(directory):
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_bytes())
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("wall mass manifest schema mismatch")
    try:
        verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCE_PATHS)
    except ValueError as error:
        raise ValueError("wall mass source binding mismatch") from error
    record, bindings = build_record()
    if manifest.get("predecessor_bindings") != bindings:
        raise ValueError("wall mass predecessor binding mismatch")
    if (set(manifest.get("outputs", {})) != {"wall-mass.json"}
            or {p.name for p in directory.iterdir()} != {"wall-mass.json", "manifest.json"}):
        raise ValueError("wall mass output inventory mismatch")
    payload = (directory / "wall-mass.json").read_bytes()
    binding = manifest["outputs"]["wall-mass.json"]
    if sha(payload) != binding["sha256"] or len(payload) != binding["bytes"]:
        raise ValueError("wall mass artifact binding mismatch")
    if payload != canonical_bytes(record):
        raise ValueError("wall mass declaration mismatch")
    return {"checked": directory.name, "outputs": 1, "sources": 3,
            "authenticated_immediate_predecessor_outputs": 2, "scientific_numerical_solves": 0}


def check_algebra(directory):
    started = time.process_time()
    check(directory)
    record = json.loads((Path(directory) / "wall-mass.json").read_bytes())
    if exact_algebra_checks(started) != record["exact_symbolic_certificate"]:
        raise ValueError("wall mass algebra replay differs from authenticated record")
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
