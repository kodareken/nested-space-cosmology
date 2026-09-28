#!/usr/bin/env python3
"""Regenerate the exact pointwise FGC-1-HYP1-DOM2-MODE1 certificate."""
from __future__ import annotations
import argparse, json, sys, tomllib
from dataclasses import asdict, dataclass, is_dataclass
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from typing import Any, Mapping

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hyp1-dom2-mode1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hyp1-dom2-mode1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hyp1-dom2-mode1.md"
from scripts.reproduce_fgc_hyp1_fo1_rc1 import IMPLEMENTATION_FILES as FO1_FILES, load_config as load_fo1
from scripts.reproduce_fgc_hyp1_dom1_qift1 import load_config as load_qift, load_canonical_result as load_qift_result
from scripts.reproduce_fgc_hyp1_con1_comp1 import load_config as load_comp, load_canonical_result as load_comp_result
from recursive_horizons.fgc.modified_harmonic_constraints import activated_compatible_state
from recursive_horizons.fgc.reference_connection import flat_spherical_annulus_reference
from recursive_horizons.fgc.spherical_reduction import state_from_generalized_adm_pg_fixture
from recursive_horizons.fgc.modified_harmonic_modes import MODE1_DERIVATIVE_STATE_ORDER, auxiliary_polynomial_mode_bases, exact_comp1_acceleration_root, mode1_point_certificate
from recursive_horizons.fgc.modified_harmonic_physical_modes import physical_mode_certificate
from recursive_horizons.fgc.reference_principal_identity import reference_principal_identity_certificate

Q = Fraction
IMPLEMENTATION_FILES = tuple(dict.fromkeys(FO1_FILES + (REPOSITORY / "src/recursive_horizons/fgc/modified_harmonic_constraints.py", REPOSITORY / "src/recursive_horizons/fgc/modified_harmonic_modes.py", REPOSITORY / "src/recursive_horizons/fgc/modified_harmonic_physical_modes.py", REPOSITORY / "src/recursive_horizons/fgc/reference_principal_identity.py", Path(__file__).resolve())))
NONCLAIMS = {name: False for name in ("uniform_open_domain_symmetrizer_proven", "uniform_radial_strong_hyperbolicity_box_proven", "metric_derived_gauge_constraint_propagation_proven", "complete_reduction_constraint_system_propagation_proven", "physical_initial_constraints_solved", "constraint_preserving_initial_boundary_value_problem_proven", "retained_eft_cutoff_and_omitted_operator_envelope_defined", "open_retained_eft_domain_proven", "evolution_authorized", "collapse_solution_derived", "metric_null_affine_defocusing_derived", "singularity_resolution_derived")}

def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected: raise ValueError(f"{name} keys differ")
def _path(name: str, value: Any) -> tuple[Path,str]:
    if not isinstance(value,str) or not value: raise ValueError(f"{name} must be a repository-relative path")
    rel=Path(value)
    if rel.is_absolute(): raise ValueError(f"{name} must be repository relative")
    p=(REPOSITORY/rel).resolve()
    try: canon=p.relative_to(REPOSITORY)
    except ValueError as exc: raise ValueError(f"{name} must stay inside repository") from exc
    if canon != rel or not p.is_file(): raise ValueError(f"{name} must be canonical, traversal free, and existing")
    return p,canon.as_posix()
def _frac(value: Any) -> Fraction:
    if isinstance(value,bool) or not isinstance(value,(int,str)): raise ValueError("fraction must be canonical")
    x=Q(value)
    if isinstance(value,str) and value != (str(x.numerator) if x.denominator==1 else f"{x.numerator}/{x.denominator}"): raise ValueError("fraction must be canonical")
    return x
def _serial(value: Any) -> Any:
    if isinstance(value,Fraction): return str(value.numerator) if value.denominator==1 else f"{value.numerator}/{value.denominator}"
    if is_dataclass(value): return _serial(asdict(value))
    if isinstance(value,Mapping): return {str(k):_serial(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [_serial(v) for v in value]
    if value is None or isinstance(value,(str,bool,int)): return value
    raise ValueError(f"unsupported certificate value {type(value).__name__}")
def _sha(path: Path) -> str: return sha256(path.read_bytes()).hexdigest()
def _rel(path: Path) -> str: return path.resolve().relative_to(REPOSITORY).as_posix()
def _pairs(pairs):
    out={}
    for k,v in pairs:
        if k in out: raise ValueError(f"duplicate JSON object key: {k}")
        out[k]=v
    return out
def _canonical(value: Any) -> str: return json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+"\n"

@dataclass(frozen=True,slots=True)
class Config:
    raw: dict[str,Any]; paths: dict[str,Path]; relatives: dict[str,str]; source_sha256: str

def load_config(path: Path=DEFAULT_CONFIG) -> Config:
    source=path.resolve().read_bytes(); raw=tomllib.loads(source.decode());
    _keys("root",raw,{"schema_version","artifact_id","project_version","metric_signature","riemann_convention","first_order_config","first_order_result","quantified_domain_config","quantified_domain_result","compatible_data_config","compatible_data_result","reference","controls","formulation","proof_contract","open_gates"})
    if {k:raw[k] for k in ("schema_version","artifact_id","project_version","metric_signature","riemann_convention")} != {"schema_version":1,"artifact_id":"FGC-1-HYP1-DOM2-MODE1","project_version":"0.11.0","metric_signature":"-+++","riemann_convention":"plus_partial_mu_gamma_nu"}: raise ValueError("MODE1 identity/version/conventions are frozen")
    expected={"first_order_config":"configs/fgc/fgc-1-hyp1-fo1-rc1.toml","first_order_result":"results/fgc-1-hyp1-fo1-rc1.json","quantified_domain_config":"configs/fgc/fgc-1-hyp1-dom1-qift1.toml","quantified_domain_result":"results/fgc-1-hyp1-dom1-qift1.json","compatible_data_config":"configs/fgc/fgc-1-hyp1-con1-comp1.toml","compatible_data_result":"results/fgc-1-hyp1-con1-comp1.json"}
    paths={}; rels={}
    for key,want in expected.items():
        p,r=_path(key,raw[key]);
        if r!=want: raise ValueError(f"{key} must name the frozen predecessor")
        paths[key],rels[key]=p,r
    fo1=load_fo1(paths["first_order_config"]); qift=load_qift(paths["quantified_domain_config"]); comp=load_comp(paths["compatible_data_config"])
    for item in (fo1,qift,comp):
        if item.project_version != raw["project_version"] or item.metric_signature != raw["metric_signature"] or item.riemann_convention != raw["riemann_convention"]: raise ValueError("predecessor conventions differ")
    if load_qift_result(paths["quantified_domain_result"]).get("gate_status",{}).get("quantified_full_dimensional_local_implicit_branch_box_passed") is not True or load_comp_result(paths["compatible_data_result"]).get("gate_status",{}).get("exact_activated_local_compatible_constraint_datum_passed") is not True: raise ValueError("required predecessor gate is not passed")
    ref=raw["reference"]; _keys("reference",ref,{"reference_id","radial_domain_minimum","center_included"})
    if ref != {"reference_id":"flat_spherical_annulus","radial_domain_minimum":"1/2","center_included":False}: raise ValueError("MODE1 reference is frozen")
    controls=raw["controls"]; _keys("controls",controls,{"flat_fixture_id","activated_fixture_id","activated_solution_status","comp1_fixture_id"})
    if controls != {"flat_fixture_id":"FGCQR_flat_reference_vacuum","activated_fixture_id":"FGCQR_activated_generic","activated_solution_status":"off_shell_control_not_a_solution","comp1_fixture_id":"FGCQR_qift1_activated_local_parameter"}: raise ValueError("MODE1 controls are frozen")
    form=raw["formulation"]; _keys("formulation",form,{"principal_relation","derivative_state_order","auxiliary_mode_basis_method","reference_principal_identity","physical_mode_basis_method","tilde_normal_factor","hat_normal_factor"})
    if form["principal_relation"]!="direct_REF1_branch_principal_equals_MHG1_standard_first_order_principal" or form["auxiliary_mode_basis_method"]!="exact_rref_nullspace_at_exact_auxiliary_polynomial_roots" or form["reference_principal_identity"]!="universal_formal_REF1_connection_variation_equals_MHG_projector_principal_identity" or form["physical_mode_basis_method"]!="exact_rr_schur_physical_and_regulator_quadratic_unit_pivot_charts_at_COMP1_root" or form["derivative_state_order"]!=list(MODE1_DERIVATIVE_STATE_ORDER) or _frac(form["tilde_normal_factor"])!=4 or _frac(form["hat_normal_factor"])!=9: raise ValueError("MODE1 formulation is frozen")
    proof=raw["proof_contract"]; _keys("proof_contract",proof,{"require_exact_ref1_mhg1_branch_equality_on_all_controls","require_exact_rref_auxiliary_bases_on_all_controls","require_exact_auxiliary_basis_residuals_zero","require_polynomial_tilde_identity_when_core_exposes_it","require_comp1_hat_quotient_atlas_identity_when_core_exposes_it","require_zero_speed_lift_when_core_exposes_it","require_exact_comp1_acceleration_root_and_zero_ref1_residual","require_universal_formal_ref1_principal_identity","require_comp1_physical_metric_and_regulator_mode_charts","require_no_floating_point_or_sampled_proof"})
    if any(v is not True for v in proof.values()) or raw["open_gates"] != NONCLAIMS: raise ValueError("MODE1 proof/open-gate contract is frozen")
    return Config(dict(raw),paths,rels,sha256(source).hexdigest())

def record(config_path: Path=DEFAULT_CONFIG) -> dict[str,Any]:
    cfg=load_config(config_path); fo1=load_fo1(cfg.paths["first_order_config"]); qift=load_qift(cfg.paths["quantified_domain_config"])
    ref=flat_spherical_annulus_reference(radial_domain_minimum=Q(1,2)); args={"reference":ref,"tilde_normal_factor":Q(4),"hat_normal_factor":Q(9)}
    comp1_base_state=activated_compatible_state(qift.flat_fixture)["state"]
    comp1_root=exact_comp1_acceleration_root(comp1_base_state,coordinate_radius=Q(4),qift_acceleration_half_width=qift.acceleration_half_width,**args)
    if not comp1_root["solved_full_residual_zero"] or not comp1_root["root_strictly_inside_qift_acceleration_box"]: raise ValueError("MODE1 exact COMP1 root certificate failed")
    comp1_solved_state=comp1_root.pop("solved_state")
    reference_identity=reference_principal_identity_certificate()
    physical_modes=physical_mode_certificate(comp1_solved_state)
    if not reference_identity["connection_variation_route_equals_projector_route"] or not physical_modes["all_exact_checks_pass"]: raise ValueError("MODE1 formal or physical mode control failed")
    controls=[("flat",state_from_generalized_adm_pg_fixture(fo1.flat_fixture),qift.coordinate_radius),("activated_off_shell",state_from_generalized_adm_pg_fixture(fo1.activated_fixture),fo1.coordinate_radii["FGCQR_activated_generic"]),("comp1_exact_root",comp1_solved_state,Q(4))]
    records={}
    for name,state,radius in controls:
        point=mode1_point_certificate(state,coordinate_radius=radius,**args)
        if name == "comp1_exact_root":
            point["auxiliary_polynomial_mode_bases"]=auxiliary_polynomial_mode_bases(state,tilde_normal_factor=Q(4),hat_normal_factor=Q(9))
        records[name]=point
    if not all(v["ref1_branch_equals_mhg1"] and v["all_auxiliary_modes_pointwise_semisimple"] for v in records.values()): raise ValueError("MODE1 point control failed")
    polynomial=records["comp1_exact_root"]["auxiliary_polynomial_mode_bases"]
    if not polynomial["tilde"]["all_residues_zero"] or not polynomial["hat"]["all_residues_zero"] or not polynomial["hat"]["free_coordinate_identity_exact"]: raise ValueError("MODE1 COMP1 polynomial mode control failed")
    optional={"comp1_exact_root":{"tilde_polynomial_identity":polynomial["tilde"]["all_residues_zero"],"hat_quotient_atlas_identity":polynomial["hat"]["free_coordinate_identity_exact"]}}
    source_paths=tuple(fo1.source_paths.values())+(cfg.paths["first_order_config"],cfg.paths["quantified_domain_config"],cfg.paths["compatible_data_config"],config_path.resolve())
    result_paths=(REPOSITORY/"results/fgc-1-action-gate.json",REPOSITORY/"results/fgc-1-metric-variation.json",REPOSITORY/"results/fgc-1-hyp1-reduction.json",REPOSITORY/"results/fgc-1-hyp1-symbol.json",REPOSITORY/"results/fgc-1-hyp1-modified-harmonic.json",REPOSITORY/"results/fgc-1-hyp1-mhg-reference.json",REPOSITORY/"results/fgc-1-hyp1-mhg-implicit.json",REPOSITORY/"results/fgc-1-hyp1-mhg-propagation.json",cfg.paths["first_order_result"],cfg.paths["quantified_domain_result"],cfg.paths["compatible_data_result"])
    return {"schema_version":1,"artifact_id":"FGC-1-HYP1-DOM2-MODE1","classification":"exact_pointwise_REF1_branch_principal_and_auxiliary_mode_basis_certificate_not_uniform_domain_or_evolution","project_version":"0.11.0","generated_by":_rel(Path(__file__)),"source_configs":cfg.relatives,"source_config_sha256":{_rel(p):_sha(p) for p in source_paths},"source_results_sha256":{_rel(p):_sha(p) for p in result_paths},"derivation_document":_rel(OWNER_DOCUMENT),"derivation_document_sha256":_sha(OWNER_DOCUMENT),"implementation_sha256":{_rel(p):_sha(p) for p in IMPLEMENTATION_FILES},"formulation":_serial(cfg.raw["formulation"]),"universal_reference_principal_identity":_serial(reference_identity),"controls":_serial(records),"exact_comp1_acceleration_root":_serial(comp1_root),"comp1_physical_mode_certificate":_serial(physical_modes),"optional_core_identities":_serial(optional),"gate_status":{"exact_pointwise_ref1_branch_mode_controls_passed":True,"uniform_radial_strong_hyperbolicity_box_passed":False,"open_retained_eft_domain_proven":False,"evolution_authorized":False},"nonclaims":NONCLAIMS}
def load_canonical_result(path: Path=DEFAULT_OUTPUT)->dict[str,Any]:
    try: source=path.read_text(); value=json.loads(source,object_pairs_hook=_pairs)
    except (OSError,json.JSONDecodeError,ValueError) as exc: raise ValueError("MODE1 result must be valid unique-key JSON") from exc
    if not isinstance(value,dict) or source != _canonical(value): raise ValueError("MODE1 result must use canonical sorted JSON")
    return value
def main()->int:
    parser=argparse.ArgumentParser(); parser.add_argument("--config",type=Path,default=DEFAULT_CONFIG); parser.add_argument("--output",type=Path,default=DEFAULT_OUTPUT); ns=parser.parse_args(); ns.output.write_text(_canonical(record(ns.config))); return 0
if __name__=="__main__": raise SystemExit(main())
