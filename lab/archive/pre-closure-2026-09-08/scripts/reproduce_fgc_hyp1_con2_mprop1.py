#!/usr/bin/env python3
"""Regenerate the exact local CON2-MPROP1 subsidiary-identity record."""
from __future__ import annotations
import argparse
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import re
import sys
import tomllib
from typing import Any, Mapping

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hyp1-con2-mprop1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hyp1-con2-mprop1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hyp1-con2-mprop1.md"
Q = Fraction

from scripts.reproduce_fgc_hyp1_dom1_qift1 import load_config as load_qift1_config
from recursive_horizons.fgc.modified_harmonic_constraints import activated_compatible_state
from recursive_horizons.fgc.modified_harmonic_metric_propagation import SphericalThirdJetState, metric_derived_gauge_propagation_identity
from recursive_horizons.fgc.modified_harmonic_reduction_subsidiary import ReductionDifferentialFieldJet, ReductionDifferentialState, reduction_subsidiary_identity
from recursive_horizons.fgc.reference_connection import SPHERICAL_FLAT_REFERENCE_ID, flat_spherical_annulus_reference
from recursive_horizons.fgc.spherical_reduction import BASE_FIELD_ORDER

IMPLEMENTATION_FILES = tuple(REPOSITORY / path for path in (
    "src/recursive_horizons/fgc/third_jet.py", "src/recursive_horizons/fgc/reference_connection.py", "src/recursive_horizons/fgc/reference_connection_second.py",
    "src/recursive_horizons/fgc/modified_harmonic_metric_propagation.py", "src/recursive_horizons/fgc/modified_harmonic_reduction_subsidiary.py", Path(__file__).relative_to(REPOSITORY).as_posix(),
))

def _keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if not isinstance(value, Mapping) or set(value) != expected: raise ValueError(f"{name} keys differ")
def _path(name: str, value: Any) -> tuple[Path,str]:
    if not isinstance(value,str) or not value: raise ValueError(f"{name} must be a repository-relative path")
    p=Path(value)
    if p.is_absolute(): raise ValueError(f"{name} must be repository relative")
    resolved=(REPOSITORY/p).resolve()
    try: rel=resolved.relative_to(REPOSITORY).as_posix()
    except ValueError as exc: raise ValueError(f"{name} must stay inside repository") from exc
    if rel!=value or not resolved.is_file(): raise ValueError(f"{name} must be canonical and existing")
    return resolved,rel
def _fraction(name: str, value: Any) -> Fraction:
    if not isinstance(value,str) or not value: raise ValueError(f"{name} must be canonical rational string")
    try: out=Q(value)
    except (ValueError,ZeroDivisionError) as exc: raise ValueError(f"{name} must be canonical rational string") from exc
    text=str(out.numerator) if out.denominator==1 else f"{out.numerator}/{out.denominator}"
    if text!=value: raise ValueError(f"{name} must be canonical rational string")
    return out
def _sha(path: Path) -> str: return sha256(path.read_bytes()).hexdigest()
def _json_object_pairs(pairs):
    out={}
    for key,value in pairs:
        if key in out: raise ValueError("JSON object must have unique-key members")
        out[key]=value
    return out
def _encode(value: Any) -> Any:
    if isinstance(value,Fraction): return str(value.numerator) if value.denominator==1 else f"{value.numerator}/{value.denominator}"
    if is_dataclass(value): return _encode(asdict(value))
    if isinstance(value,Mapping): return {str(k):_encode(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)): return [_encode(v) for v in value]
    if isinstance(value,bool) or value is None or isinstance(value,str): return value
    if isinstance(value,int): return value
    raise TypeError(f"cannot canonically encode {type(value).__name__}")
def _canonical_json_text(value: Any) -> str: return json.dumps(_encode(value),sort_keys=True,indent=2,ensure_ascii=False)+"\n"
def load_canonical_result(path: Path=DEFAULT_OUTPUT)->dict[str,Any]:
    raw=path.read_text(encoding="utf-8"); value=json.loads(raw,object_pairs_hook=_json_object_pairs)
    if _canonical_json_text(value)!=raw: raise ValueError("result is not canonical sorted JSON")
    def restore(item: Any) -> Any:
        if isinstance(item,dict): return {key:restore(value) for key,value in item.items()}
        if isinstance(item,list): return [restore(value) for value in item]
        if isinstance(item,str) and re.fullmatch(r"-?(?:0|[1-9][0-9]*)(?:/(?:[1-9][0-9]*))?",item): return Q(item)
        return item
    return restore(value)
def load_config(path: Path=DEFAULT_CONFIG)->dict[str,Any]:
    raw=tomllib.loads(path.read_text(encoding="utf-8")); _keys("root",raw,{"schema_version","artifact_id","project_version","metric_signature","riemann_convention","predecessors","reference","controls","proof_contract","open_gates"})
    if (raw["schema_version"],raw["artifact_id"],raw["project_version"],raw["metric_signature"],raw["riemann_convention"]) != (1,"FGC-1-HYP1-CON2-MPROP1","0.11.0","-+++","plus_partial_mu_gamma_nu"): raise ValueError("frozen scalar contract differs")
    _keys("predecessors",raw["predecessors"],{"fo1_config","fo1_result","qift1_config","qift1_result","comp1_config","comp1_result"})
    paths={key:_path(f"predecessors.{key}",value) for key,value in raw["predecessors"].items()}
    expected_predecessors={"fo1":"FGC-1-HYP1-FO1-RC1","qift1":"FGC-1-HYP1-DOM1-QIFT1","comp1":"FGC-1-HYP1-CON1-COMP1"}
    for prefix,artifact_id in expected_predecessors.items():
        result_path,_=paths[f"{prefix}_result"]; config_path,_=paths[f"{prefix}_config"]
        predecessor=load_canonical_result(result_path)
        if predecessor.get("artifact_id")!=artifact_id or predecessor.get("project_version")!=raw["project_version"]: raise ValueError("predecessor result provenance differs")
        source_ledger=predecessor.get("source_config_sha256")
        if not isinstance(source_ledger,Mapping) or source_ledger.get(config_path.relative_to(REPOSITORY).as_posix())!=_sha(config_path): raise ValueError("predecessor result is stale relative to its frozen config")
    _keys("reference",raw["reference"],{"reference_id","radial_domain_minimum","coordinate_radius","center_included","tilde_normal_factor","hat_normal_factor"})
    r=raw["reference"]
    if r["reference_id"]!=SPHERICAL_FLAT_REFERENCE_ID or r["center_included"] is not False: raise ValueError("reference contract differs")
    for key,expected in (("radial_domain_minimum",Q(1,2)),("coordinate_radius",Q(4)),("tilde_normal_factor",Q(4)),("hat_normal_factor",Q(9))):
        if _fraction(f"reference.{key}",r[key])!=expected: raise ValueError("reference rational contract differs")
    _keys("controls",raw["controls"],{"flat_id","activated_id","third_derivative_formula"})
    if raw["controls"]!={"flat_id":"flat_zero_third_jet_reference","activated_id":"activated_comp1_deterministic_third_jet","third_derivative_formula":"dttt=i/17;dttr=-i/19;dtrr=i/23;drrr=-i/29_in_BASE_FIELD_ORDER_i=1_to_6"}: raise ValueError("controls contract differs")
    _keys("proof_contract",raw["proof_contract"],{"metric_derived_C_through_second_partials","act1_noether_sign_checked","prop1_two_route_equality_checked","complete_kinematic_1plus1_reduction_subsidiary_checked"})
    if not all(raw["proof_contract"].values()): raise ValueError("proof contract is not frozen true")
    expected_open={"Cauchy_uniqueness_and_zero_constraint_propagation_proven","physical_initial_constraint_hypersurface_solved","physical_Hamiltonian_momentum_constraint_propagation_proven","constraint_preserving_boundary_conditions_derived","boundary_stability_or_Kreiss_estimate_proven","initial_boundary_value_problem_proven","evolution_authorized","collapse_solution_derived","finite_invariant_transition_surface_derived","singularity_resolution_derived"}
    _keys("open_gates",raw["open_gates"],expected_open)
    if any(raw["open_gates"].values()): raise ValueError("open gate is promoted")
    return {
        "raw": raw,
        "paths": paths,
        "config_path": path.resolve(),
        "source_sha256": _sha(path),
    }
def _third_table(*, activated: bool)->dict[str,dict[str,Fraction]]:
    return {field:{"dttt":Q(0) if not activated else Q(i,17),"dttr":Q(0) if not activated else -Q(i,19),"dtrr":Q(0) if not activated else Q(i,23),"drrr":Q(0) if not activated else -Q(i,29)} for i,field in enumerate(BASE_FIELD_ORDER,1)}
def _reduction_state(state)->ReductionDifferentialState:
    fields={}
    for field in BASE_FIELD_ORDER:
        j=getattr(state,field); fields[field]=ReductionDifferentialFieldJet(u_t=j.dt,p=j.dt,q=j.dr,u_r=j.dr,q_t=j.dtr,p_r=j.dtr,u_tr=j.dtr)
    return ReductionDifferentialState(**fields)
def record(config_path: Path=DEFAULT_CONFIG)->dict[str,Any]:
    config=load_config(config_path); raw=config["raw"]; qift=load_qift1_config(config["paths"]["qift1_config"][0])
    datum=activated_compatible_state(qift.flat_fixture)
    reference=flat_spherical_annulus_reference(radial_domain_minimum=Q(1,2))
    controls={}
    for name,state,active in ((raw["controls"]["flat_id"],qift.flat_fixture,False),(raw["controls"]["activated_id"],datum["state"],True)):
        base = datum["state"] if active else __import__('recursive_horizons.fgc.spherical_reduction',fromlist=['state_from_generalized_adm_pg_fixture']).state_from_generalized_adm_pg_fixture(state)
        third=SphericalThirdJetState.from_spherical_state(base,_third_table(activated=active))
        metric=metric_derived_gauge_propagation_identity(third,reference=reference,coordinate_radius=Q(4),tilde_normal_factor=Q(4),hat_normal_factor=Q(9))
        reduction=reduction_subsidiary_identity(_reduction_state(base))
        if not metric["local_metric_derived_gauge_subsidiary_identity_derived"] or not metric["metric_derived_extension_operator"]["two_routes_agree_exactly"] or not reduction["off_shell_identity_exact"]: raise ValueError("core local identity gate failed")
        controls[name]={"third_derivatives":_third_table(activated=active),"metric_derived_gauge_subsidiary":metric,"reduction_subsidiary":reduction}
    source_files={p.relative_to(REPOSITORY).as_posix():_sha(p) for p in IMPLEMENTATION_FILES}
    predecessor_hashes={rel:_sha(path) for path,rel in config["paths"].values()}
    result={"schema_version":1,"artifact_id":raw["artifact_id"],"project_version":raw["project_version"],"classification":"exact_local_metric_derived_gauge_and_kinematic_reduction_subsidiary_identities_not_Cauchy_or_IBVP","generated_by":"scripts/reproduce_fgc_hyp1_con2_mprop1.py","source_config_sha256":{config["config_path"].relative_to(REPOSITORY).as_posix():config["source_sha256"]},"derivation_document":"docs/fgc-hyp1-con2-mprop1.md","derivation_document_sha256":_sha(OWNER_DOCUMENT),"predecessor_sha256":predecessor_hashes,"implementation_sha256":source_files,"proof_contract":raw["proof_contract"],"controls":controls,"gate_status":{"local_differential_identity":True,"metric_derived_C_through_second_partials":True,"ACT1_noether_sign":True,"PROP1_two_route_equality":True,"complete_kinematic_1plus1_reduction_subsidiary":True},"nonclaims":raw["open_gates"]}
    return result
def main()->None:
    parser=argparse.ArgumentParser(); parser.add_argument("--config",type=Path,default=DEFAULT_CONFIG); parser.add_argument("--output",type=Path,default=DEFAULT_OUTPUT); args=parser.parse_args(); args.output.write_text(_canonical_json_text(record(args.config)),encoding="utf-8")
if __name__=="__main__": main()
