#!/usr/bin/env python3
"""Regenerate the exact compact radial FGC-1-HYP1-DOM3-UHYP1 record."""
from __future__ import annotations
import argparse, json, sys, tomllib
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
# Exact interval frame enclosures can contain integers with more than the
# interpreter's default decimal-conversion digit cap.  They are deliberate
# machine-proof data, not untrusted user-supplied integers.
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hyp1-dom3-uhyp1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hyp1-dom3-uhyp1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hyp1-dom3-uhyp1.md"

from scripts.reproduce_fgc_hyp1_dom2_mode1 import load_config as load_mode1, load_canonical_result as load_mode1_result
from scripts.reproduce_fgc_hyp1_dom1_qift1 import load_config as load_qift1, load_canonical_result as load_qift1_result
from scripts.reproduce_fgc_hyp1_con1_comp1 import load_config as load_comp1, load_canonical_result as load_comp1_result
from recursive_horizons.fgc.modified_harmonic_constraints import activated_compatible_state
from recursive_horizons.fgc.modified_harmonic_uniform_domain import uniform_radial_hyperbolicity_certificate

Q = Fraction
NONCLAIMS = {name: False for name in (
    "multidirectional_strong_hyperbolicity_proven", "constraint_propagation_proven", "IBVP_proven",
    "retained_EFT_domain_proven", "evolution_authorized", "collapse_solution_derived",
    "metric_null_affine_defocusing_derived", "singularity_resolution_derived")}
IMPLEMENTATION = tuple(REPOSITORY / p for p in (
    "src/recursive_horizons/fgc/exact_interval.py", "src/recursive_horizons/fgc/interval_tangent.py",
    "src/recursive_horizons/fgc/exact_interval_krawczyk.py", "src/recursive_horizons/fgc/exact_interval_linear_algebra.py",
    "src/recursive_horizons/fgc/modified_harmonic_interval_principal.py",
    "src/recursive_horizons/fgc/modified_harmonic_auxiliary_identity.py",
    "src/recursive_horizons/fgc/modified_harmonic_modes.py", "src/recursive_horizons/fgc/modified_harmonic_physical_modes.py",
    "src/recursive_horizons/fgc/modified_harmonic_uniform_domain.py", "scripts/reproduce_fgc_hyp1_dom3_uhyp1.py"))

def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected: raise ValueError(f"{name} keys differ: missing={sorted(expected-set(value))}, extra={sorted(set(value)-expected)}")
def _text(q: Fraction) -> str: return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"
def _frac(name: str, value: Any) -> Fraction:
    if not isinstance(value, str) or not value: raise ValueError(f"{name} must be a canonical rational string")
    try: q=Q(value)
    except (ValueError,ZeroDivisionError) as exc: raise ValueError(f"{name} must be a canonical rational string") from exc
    if value != _text(q): raise ValueError(f"{name} must be a canonical rational string")
    if q <= 0: raise ValueError(f"{name} must be positive")
    return q
def _path(name: str, value: Any, expected: str) -> Path:
    if not isinstance(value,str) or value != expected: raise ValueError(f"{name} must name {expected}")
    p=(REPOSITORY/Path(value)).resolve()
    if p.relative_to(REPOSITORY).as_posix()!=value or not p.is_file(): raise ValueError(f"{name} must be canonical, traversal free, and existing")
    return p
def _sha(p: Path) -> str: return sha256(p.read_bytes()).hexdigest()
def _rel(p: Path) -> str: return p.resolve().relative_to(REPOSITORY).as_posix()
def _serial(x: Any) -> Any:
    if isinstance(x,Fraction): return _text(x)
    if is_dataclass(x): return _serial(asdict(x))
    if isinstance(x,Mapping): return {str(k):_serial(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [_serial(v) for v in x]
    if x is None or isinstance(x,(str,bool,int)): return x
    raise TypeError(f"unsupported exact certificate type {type(x).__name__}")
def _pairs(pairs):
    out={}
    for k,v in pairs:
        if k in out: raise ValueError(f"duplicate JSON object key: {k}")
        out[k]=v
    return out
def _canonical_json_text(value: Any) -> str: return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+"\n"

def _serialized_sha256(value: Any) -> str:
    return sha256(_canonical_json_text(value).encode("utf-8")).hexdigest()

def _compact_krawczyk(value: Mapping[str, Any]) -> dict[str, Any]:
    keep = (
        "conditional_theorem", "dimension", "input_enclosure_proven_here",
        "krawczyk_image_strictly_inside_displacement_box",
        "minimum_strict_componentwise_inclusion_margin", "norm",
        "rho_infinity_upper_bound", "rho_strictly_below_one",
        "strict_componentwise_inclusion_margins",
        "uniform_inverse_infinity_norm_upper_bound",
    )
    answer = {key: value[key] for key in keep}
    answer["full_exact_krawczyk_payload_sha256"] = _serialized_sha256(value)
    return answer

def _compact_mode(value: Mapping[str, Any]) -> dict[str, Any]:
    keep = (
        "speed_box", "source", "full_six_row_reason", "sixth_row_reason",
        "free_column", "scale", "center_speed", "center_seed_speed",
    )
    answer = {key: value[key] for key in keep if key in value}
    if "linear_pivot_krawczyk" in value:
        answer["linear_pivot_krawczyk"] = _compact_krawczyk(value["linear_pivot_krawczyk"])
    if "nonlinear_full_five_row_krawczyk" in value:
        answer["nonlinear_full_five_row_krawczyk"] = _compact_krawczyk(value["nonlinear_full_five_row_krawczyk"])
    answer["full_exact_mode_payload_sha256"] = _serialized_sha256(value)
    return answer

def _compact_interval_principal(value: Mapping[str, Any]) -> dict[str, Any]:
    box = value["box"]
    kraw = value["complete_REF1_acceleration_Krawczyk"]
    kinetic = value["kinetic_neumann_inverse"]
    return {
        "classification": value["classification"],
        "box": {
            key: box[key]
            for key in (
                "acceleration_half_width", "all_30_parameter_axes_have_nonzero_width",
                "all_6_acceleration_axes_have_nonzero_width", "argument_order",
                "differentiation_method", "equation_order", "parameter_half_width",
                "seeded_primals_reproduce_same_complete_residual_box", "state_order",
            )
        },
        "complete_REF1_acceleration_Krawczyk": {
            key: kraw[key]
            for key in (
                "contraction_infinity_norm_upper_bound", "contraction_strictly_less_than_one",
                "fixed_point_iterations_converge_from_every_point_in_declared_acceleration_box",
                "inverse_center_jacobian_infinity_norm",
                "krawczyk_image_strictly_inside_acceleration_box",
                "minimum_strict_interior_inclusion_margin", "norm",
                "strict_interior_inclusion_margins", "theorem_route",
                "uniform_acceleration_jacobian_inverse_infinity_norm_upper_bound",
                "unique_acceleration_root_for_every_declared_parameter_point",
            )
        },
        "kinetic_neumann_inverse": {
            key: kinetic[key]
            for key in (
                "exact_center_inverse_infinity_norm", "inverse_difference_infinity_upper_bound",
                "kinetic_inverse_regular_over_entire_box", "neumann_rho_infinity_upper_bound",
            )
        },
        "orders": value["orders"],
        "regularity": value["regularity"],
        "nonclaims": value["nonclaims"],
        "omitted_exact_box_payload_sha256": _serialized_sha256(box),
        "omitted_exact_acceleration_Krawczyk_payload_sha256": _serialized_sha256(kraw),
        "omitted_exact_kinetic_Neumann_payload_sha256": _serialized_sha256(kinetic),
        "omitted_exact_solved_branch_principal_box_sha256": _serialized_sha256(value["solved_branch_principal_box"]),
        "full_exact_interval_principal_payload_sha256": _serialized_sha256(value),
    }

def _compact_uniform_certificate(value: Mapping[str, Any]) -> dict[str, Any]:
    frame = value["eigenframe"]
    neumann = frame["neumann_inverse"]
    modes = value["modes"]
    return {
        "classification": value["classification"],
        "domain": value["domain"],
        "center_root": value["center_root"],
        "interval_principal": _compact_interval_principal(value["interval_principal"]),
        "cone_quadratics": value["cone_quadratics"],
        "universal_auxiliary_identities": value["universal_auxiliary_identities"],
        "point_atlases": value["point_atlases"],
        "modes": {
            sector: [_compact_mode(mode) for mode in modes[sector]]
            for sector in ("tilde", "hat", "physical_chi", "regulator")
        } | {"correlated_eigenmode_routes": modes["correlated_eigenmode_routes"]},
        "eigenframe": {
            "all_enclosed_frames_invertible": frame["all_enclosed_frames_invertible"],
            "real_smooth_radial_eigenframe": frame["real_smooth_radial_eigenframe"],
            "neumann_inverse": {
                key: neumann[key]
                for key in (
                    "center_inverse_infinity_norm", "dimension",
                    "every_enclosed_matrix_invertible", "inverse_tail_entrywise_bound",
                    "method", "rho_infinity", "rho_strictly_below_one",
                )
            },
            "exact_center_frame_sha256": _serialized_sha256(frame["exact_center_frame"]),
            "interval_columns_sha256": _serialized_sha256(frame["interval_columns"]),
            "full_exact_neumann_payload_sha256": _serialized_sha256(neumann),
            "full_exact_eigenframe_payload_sha256": _serialized_sha256(frame),
        },
        "radial_symmetrizer": value["radial_symmetrizer"],
        "nonclaims": value["nonclaims"],
    }

def _compact_existing_result(path: Path) -> dict[str, Any]:
    value = load_canonical_result(path)
    full = value.get("uniform_radial_certificate")
    if not isinstance(full, Mapping) or "full_exact_eigenframe_payload_sha256" in full.get("eigenframe", {}):
        raise ValueError("compact-existing requires the full unprojected UHYP1 result")
    full_text = _canonical_json_text(full)
    value["full_uniform_radial_certificate_sha256"] = sha256(full_text.encode("utf-8")).hexdigest()
    value["full_uniform_radial_certificate_canonical_json_bytes"] = len(full_text.encode("utf-8"))
    value["uniform_radial_certificate"] = _compact_uniform_certificate(full)
    value["implementation_sha256"][_rel(Path(__file__))] = _sha(Path(__file__))
    value["derivation_document_sha256"] = _sha(OWNER_DOCUMENT)
    return value

def load_config(path: Path=DEFAULT_CONFIG) -> dict[str, Any]:
    raw=tomllib.loads(path.resolve().read_text(encoding="utf-8"))
    _keys("root",raw,{"schema_version","artifact_id","project_version","metric_signature","riemann_convention","mode1_config","mode1_result","qift1_config","qift1_result","comp1_config","comp1_result","reference","box","proof_contract","open_gates"})
    if {k:raw[k] for k in ("schema_version","artifact_id","project_version","metric_signature","riemann_convention")} != {"schema_version":1,"artifact_id":"FGC-1-HYP1-DOM3-UHYP1","project_version":"0.11.0","metric_signature":"-+++","riemann_convention":"plus_partial_mu_gamma_nu"}: raise ValueError("UHYP1 identity/version/conventions are frozen")
    paths={k:_path(k,raw[k],v) for k,v in {"mode1_config":"configs/fgc/fgc-1-hyp1-dom2-mode1.toml","mode1_result":"results/fgc-1-hyp1-dom2-mode1.json","qift1_config":"configs/fgc/fgc-1-hyp1-dom1-qift1.toml","qift1_result":"results/fgc-1-hyp1-dom1-qift1.json","comp1_config":"configs/fgc/fgc-1-hyp1-con1-comp1.toml","comp1_result":"results/fgc-1-hyp1-con1-comp1.json"}.items()}
    mode,qift,comp=load_mode1(paths["mode1_config"]),load_qift1(paths["qift1_config"]),load_comp1(paths["comp1_config"])
    def convention(x: Any, key: str) -> Any:
        return x.raw[key] if hasattr(x, "raw") else getattr(x, key)
    if any(convention(x,"project_version") != "0.11.0" or convention(x,"metric_signature") != "-+++" or convention(x,"riemann_convention") != "plus_partial_mu_gamma_nu" for x in (mode,qift,comp)): raise ValueError("predecessor conventions differ")
    if not load_mode1_result(paths["mode1_result"])["gate_status"].get("exact_pointwise_ref1_branch_mode_controls_passed") or not load_qift1_result(paths["qift1_result"])["gate_status"].get("quantified_full_dimensional_local_implicit_branch_box_passed") or not load_comp1_result(paths["comp1_result"])["gate_status"].get("exact_activated_local_compatible_constraint_datum_passed"): raise ValueError("required predecessor gate is not passed")
    _keys("reference",raw["reference"],{"reference_id","radial_domain_minimum","coordinate_radius","tilde_normal_factor","hat_normal_factor"})
    if raw["reference"]["reference_id"]!="flat_spherical_annulus" or raw["reference"]["radial_domain_minimum"]!="1/2": raise ValueError("reference is frozen")
    ref={k:_frac(f"reference.{k}",raw["reference"][k]) for k in ("coordinate_radius","tilde_normal_factor","hat_normal_factor")}
    if ref != {"coordinate_radius":Q(4),"tilde_normal_factor":Q(4),"hat_normal_factor":Q(9)}: raise ValueError("reference factors are frozen")
    _keys("box",raw["box"],{"parameter_half_width","acceleration_half_width","root_bracket_half_width","hat_mode_displacement_width","regulator_mode_displacement_width"})
    box={k:_frac(f"box.{k}",v) for k,v in raw["box"].items()}
    if box != {"parameter_half_width":Q(1,2**300),"acceleration_half_width":Q(1,2**290),"root_bracket_half_width":Q(1,2**200),"hat_mode_displacement_width":Q(1,2**140),"regulator_mode_displacement_width":Q(1,2**68)}: raise ValueError("UHYP1 widths are frozen")
    _keys("proof_contract",raw["proof_contract"],{"certificate","radial_only","all_axes_strictly_nonzero","exact_rational_intervals_only","no_determinant_factor_or_interval_residual_zero_shortcut"})
    if raw["proof_contract"] != {"certificate":"complete_REF1_interval_branch_plus_explicit_sector_modes_plus_frame_Neumann","radial_only":True,"all_axes_strictly_nonzero":True,"exact_rational_intervals_only":True,"no_determinant_factor_or_interval_residual_zero_shortcut":True}: raise ValueError("proof contract is frozen")
    if raw["open_gates"] != NONCLAIMS: raise ValueError("open gate contract is frozen")
    return {"raw":raw,"paths":paths,"box":box,"reference":ref,"source_sha256":_sha(path.resolve())}

def record(config_path: Path=DEFAULT_CONFIG) -> dict[str,Any]:
    cfg=load_config(config_path); qift=load_qift1(cfg["paths"]["qift1_config"])
    datum=activated_compatible_state(qift.flat_fixture)
    cert=uniform_radial_hyperbolicity_certificate(datum["state"],reference=datum["reference"],coordinate_radius=cfg["reference"]["coordinate_radius"],parameter_half_width=cfg["box"]["parameter_half_width"],acceleration_half_width=cfg["box"]["acceleration_half_width"],root_bracket_half_width=cfg["box"]["root_bracket_half_width"],mode_displacement_width=cfg["box"]["hat_mode_displacement_width"],regulator_displacement_width=cfg["box"]["regulator_mode_displacement_width"])
    if not cert["eigenframe"]["all_enclosed_frames_invertible"] or cert["eigenframe"]["neumann_inverse"]["rho_infinity"] >= 1: raise ValueError("UHYP1 frame gate failed")
    full=_serial(cert); full_text=_canonical_json_text(full)
    return {"schema_version":1,"artifact_id":"FGC-1-HYP1-DOM3-UHYP1","classification":"exact_rational_compact_nonflat_radial_REF1_strong_hyperbolicity_on_implicit_branch_graph","project_version":"0.11.0","generated_by":_rel(Path(__file__)),"derivation_document":_rel(OWNER_DOCUMENT),"derivation_document_sha256":_sha(OWNER_DOCUMENT),"source_config_sha256":{_rel(config_path.resolve()):cfg["source_sha256"]},"predecessor_sha256":{_rel(p):_sha(p) for p in cfg["paths"].values()},"implementation_sha256":{_rel(p):_sha(p) for p in IMPLEMENTATION},"frozen_configuration":_serial(cfg["raw"]),"full_uniform_radial_certificate_sha256":sha256(full_text.encode("utf-8")).hexdigest(),"full_uniform_radial_certificate_canonical_json_bytes":len(full_text.encode("utf-8")),"uniform_radial_certificate":_compact_uniform_certificate(full),"gate_status":{"compact_nonflat_radial_strong_hyperbolicity_on_implicit_branch_graph_passed":True,"multidirectional_strong_hyperbolicity_proven":False,"constraint_propagation_proven":False,"IBVP_proven":False,"retained_EFT_domain_proven":False,"evolution_authorized":False},"nonclaims":NONCLAIMS}
def load_canonical_result(path: Path=DEFAULT_OUTPUT)->dict[str,Any]:
    try: text=path.read_text(encoding="utf-8"); value=json.loads(text,object_pairs_hook=_pairs)
    except (OSError,json.JSONDecodeError,ValueError) as exc: raise ValueError("UHYP1 result must be valid unique-key JSON") from exc
    if not isinstance(value,dict) or text != _canonical_json_text(value): raise ValueError("UHYP1 result must use canonical sorted JSON")
    return value
def main()->int:
    ap=argparse.ArgumentParser(); ap.add_argument("--config",type=Path,default=DEFAULT_CONFIG); ap.add_argument("--output",type=Path,default=DEFAULT_OUTPUT); ap.add_argument("--compact-existing-full-result",action="store_true"); ns=ap.parse_args(); payload=_compact_existing_result(ns.output) if ns.compact_existing_full_result else record(ns.config); ns.output.write_text(_canonical_json_text(payload),encoding="utf-8"); return 0
if __name__=="__main__": raise SystemExit(main())
