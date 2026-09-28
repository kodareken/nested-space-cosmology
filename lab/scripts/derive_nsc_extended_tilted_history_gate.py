#!/usr/bin/env python3
"""Compose the three tilted/frequency-mixing owners without inventing a selector."""
from __future__ import annotations

import argparse,hashlib,json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/"src"))
from recursive_horizons.nsc_extended_tilted_history_gate import compose_extended_gate

OUTPUT=ROOT/"results/development/nsc-extended-tilted-history-gate.json"
INPUTS=(
 "results/development/nsc-general-ks-reference.json",
 "results/development/nsc-general-ks-local-history.json",
 "results/development/nsc-tilted-landau-interface.json",
 "results/development/nsc-general-ks-same-action-history.json",
 "results/development/nsc-mode-resolved-cauchy-state.json",
 "results/development/scale-binding.json",
 "results/development/compact-boundary-action.json",
)
SOURCES=("src/recursive_horizons/nsc_extended_tilted_history_gate.py","scripts/derive_nsc_extended_tilted_history_gate.py","docs/nsc-extended-tilted-history-gate.md")

def hashes(paths):return {p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
def compare(a,b,path="$",atol=3e-13,rtol=3e-13):
 if isinstance(a,dict):
  if not isinstance(b,dict) or a.keys()!=b.keys():raise AssertionError(f"keys differ {path}")
  for k in a:compare(a[k],b[k],f"{path}/{k}",atol,rtol)
 elif isinstance(a,list):
  if not isinstance(b,list) or len(a)!=len(b):raise AssertionError(f"list differs {path}")
  for i,(x,y) in enumerate(zip(a,b)):compare(x,y,f"{path}/{i}",atol,rtol)
 elif isinstance(a,float):
  if isinstance(b,bool) or not isinstance(b,(int,float)) or abs(a-b)>atol+rtol*abs(a):raise AssertionError(f"number differs {path}")
 elif type(a) is not type(b) or a!=b:raise AssertionError(f"value differs {path}")

def calculate():
 reference,local,interface,old,state,scale,boundary=(json.loads((ROOT/p).read_text()) for p in INPUTS)
 gate=compose_extended_gate(reference,local,interface,old)
 if not gate["components_expose_residuals"]:raise ArithmeticError("a component owner is not executable")
 if gate["extended_functional_complete"]:raise ArithmeticError("record expects the unresolved extended branch")
 branch=scale["development_branch"];locked=state["locked_inputs"]
 if not (locked["magnetic_flux"]==branch["magnetic_flux"]==4 and locked["Omega"]==branch["omega"] and locked["zeta"]==branch["zeta"] and locked["V_full"]==0):raise ValueError("locked branch changed")
 lr=local["residuals"];ir=interface["residuals"];rr=reference["residuals"]
 return {
  "schema":"NSC-EXTENDED-TILTED-HISTORY-GATE-v1",
  "status":"COMPONENT CONSTRUCTION PASS / EXTENDED STATIONARITY OPEN: all three owners expose residuals, but the physical mixing kernel and Weyl endpoint completion are not selected",
  "source_hashes":hashes(SOURCES),"input_hashes":hashes(INPUTS),
  "locked_inputs":{"A":locked["A"],"magnetic_flux":locked["magnetic_flux"],"Omega":locked["Omega"],"zeta":locked["zeta"],"V_full":locked["V_full"],"refit":False},
  "component_gates":{
   "GeneralKSFourthOrderReferenceHistory":{"status":reference["status"],"pass":gate["reference_pass"]},
   "GeneralKSLocalInducedHistory":{"status":local["status"],"bulk_pass":gate["local_bulk_pass"],"endpoint_pass":gate["weyl_endpoint_completion_pass"]},
   "TransmittingTiltedLandauInterface":{"status":interface["status"],"construction_pass":gate["tilted_interface_construction_pass"],"selector_pass":gate["physical_kernel_selector_pass"]},
  },
  "residual_vector":{
   "reference_bloch":rr["maximum_bloch_recursion_residual"],
   "reference_order_normalization":rr["maximum_order_normalization_residual"],
   "reference_parity_beta":rr["maximum_parity_completed_beta_projection"],
   "local_node_gradient":lr["directional_gradient_check"]["maximum_absolute"],
   "local_weyl_free_endpoint":lr["weyl_endpoint_completion"]["maximum_absolute"],
   "interface_variation":ir["maximum_interface_variation_max_abs"],
   "interface_unitarity":ir["maximum_canonical_unitarity_max_abs"],
   "interface_CAR":ir["maximum_CAR_spectrum_transport_max_abs"],
  },
  "tolerances":{
   "reference":3e-8,"local_node_gradient":lr["directional_gradient_check"]["tolerance"],
   "stationarity":lr["weyl_endpoint_completion"]["tolerance"],"interface":ir["declared_tolerance"],
  },
  "frequency_mixing_family":{
   "full_rank":interface["rank_nonuniqueness_certificate"]["full_block_map_rank"],"rank_defect":interface["rank_nonuniqueness_certificate"]["full_block_map_rank_defect"],
   "constraint_rank":interface["rank_nonuniqueness_certificate"]["weighted_isometry_constraint_jacobian_rank"],"admissible_nullity":interface["rank_nonuniqueness_certificate"]["admissible_family_real_dimension_and_linearized_nullity"],
   "physical_Vc_selected":False,
  },
  "composition":gate,
  "blocking_residuals":[
   {"owner":"GeneralKSLocalInducedHistory","channel":"free-endpoint Weyl momentum","value":lr["weyl_endpoint_completion"]["maximum_absolute"],"tolerance":lr["weyl_endpoint_completion"]["tolerance"]},
   {"owner":"TransmittingTiltedLandauInterface","channel":"physical frequency-mixing kernel V_c","value":None,"admissible_family_dimension":interface["rank_nonuniqueness_certificate"]["admissible_family_real_dimension_and_linearized_nullity"]},
  ],
  "exact_next_owner":{
   "name":"ModeResolvedTransmittingBoundaryHistoryAction",
   "already_declared_role":"Gamma_rest/interface term in pi_parent+pi_child+delta Gamma_rest/delta h=0",
   "must_determine":["one V_c per retained channel from the existing link/boundary dynamics","tilted embedding and frequency mixing","metric/interface variation that completes the Weyl endpoint momentum"],
   "new_physical_term_required":False,
  },
  "extended_stationarity":{
   "delta_Gamma_assembled":False,"optimizer_started":False,"existence_claimed":False,"nonexistence_claimed":False,
   "finite_stress":None,"null_signs":None,"updated_constraints":None,
   "reason":"the action derivative is not single-valued until V_c and the matching endpoint variation are selected by the same boundary action",
  },
  "regression":{"homogeneous_no_interface_nonexistence_still_passes":gate["killed_homogeneous_regression_pass"]},
  "gate":{"three_records_real":True,"extended_composition_ready":False,"coupled_evolution_reopened":False,"finite_stress_fabricated":False},
  "scope":{"metric_timestep_started":False,"diagnostic_kernel_selected_physical":False,"A_q_Omega_zeta_or_V_full_changed":False,"new_counterflow_dark_fluid_or_compensator":False,"old_generators_rerun":False},
  "comparison":{"fields":"all","exact":"schema, strings, booleans and source/input hashes","float_atol":3e-13,"float_rtol":3e-13,"exceptions":[]},
 }

def main():
 p=argparse.ArgumentParser(description=__doc__);g=p.add_mutually_exclusive_group();g.add_argument("--check",action="store_true");g.add_argument("--output",type=Path);a=p.parse_args();r=calculate()
 if a.check:compare(json.loads(OUTPUT.read_text()),r);print("extended tilted history component gate reproduced")
 elif a.output:
  if a.output.exists():raise FileExistsError("refusing overwrite")
  a.output.write_text(json.dumps(r,indent=2,sort_keys=True)+"\n");print(a.output)
 else:print(json.dumps(r,indent=2,sort_keys=True))
if __name__=="__main__":main()
