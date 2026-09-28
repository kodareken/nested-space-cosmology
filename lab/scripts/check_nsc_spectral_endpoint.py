#!/usr/bin/env python3
"""Test the missing local/full spectral source interface at the child endpoint."""
import argparse
import json
from math import pi, sqrt, exp
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.special import exp1
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from recursive_horizons.nsc_spectral_endpoint import (
    endpoint,geometry_scales,local_endpoint,sphere_spectrum,
)
from recursive_horizons.nsc_warped_source import WarpedCompactWeight
from recursive_horizons.nsc_finite_terms import _symbolic_basis
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare,native
from check_nsc_vacuum_charge_matching import hashes

OUTPUT=ROOT/"results/development/spectral-endpoint.json"
SOURCES=("scripts/check_nsc_spectral_endpoint.py","src/recursive_horizons/nsc_spectral_endpoint.py",
         "docs/nsc-spectral-endpoint.md","src/recursive_horizons/nsc_lorentzian.py",
         "src/recursive_horizons/nsc_finite_terms.py","scripts/check_nsc_compact_boundary_action.py",
         "scripts/check_nsc_compact_casimir.py","scripts/check_nsc_vacuum_charge_matching.py")
INPUTS=("results/development/compact-matching.json","results/development/horizon-source.json",
        "results/development/curvature-eft.json","results/nsc-5-clock-horizon.json")


def calculate():
    match,horizon,_,clock=[authenticated_record(p) for p in INPUTS]
    parameters=match["conventions"]["parameters"]
    coefficients=next(c for c in match["matched_coefficients"] if c["matching_cutoff"]==1.)
    cutoff=parameters["Lambda"]
    def weight(modes=34,points=112,path=None):
        return WarpedCompactWeight(cutoff=cutoff,interval=parameters["ell"],modes=modes,points=points,
                                   amplitude=parameters["amplitude"],path=parameters["warp_path"] if path is None else path)
    exact={}
    def zero(name,value):
        assert sp.simplify(value)==0,name
        exact[name]="0"
    # Bind the already-checked sectional contractions to the quasi-global chart.
    jets,_,sections,*_= _symbolic_basis()
    F,Fp,Fpp,r,rp,rpp=sp.symbols("F Fp Fpp r rp rpp",positive=True)
    values=((sp.sqrt(F),Fp/(2*sp.sqrt(F)),Fpp/(2*sp.sqrt(F))-Fp**2/(4*F**sp.Rational(3,2))),
            (1/sp.sqrt(F),-Fp/(2*F**sp.Rational(3,2)),0),(r,rp,rpp))
    substitution={a:b for row,v in zip(jets,values) for a,b in zip(row,v)}
    targets=(Fpp/2,Fp*rp/(2*r),F*rpp/r+Fp*rp/(2*r),(1-F*rp**2)/r**2)
    for i,(old,new) in enumerate(zip(sections,targets)):
        zero(f"stored_section_to_quasiglobal_{i}",old.subs(substitution)-new)
    # The trace factor uses BOTH signs of the published self-adjoint spectrum.
    m=sp.Symbol("m",integer=True,positive=True)
    zero("signed_to_squared_spectrum_multiplicity",2*4*(m-1)*m*(m+1)/6-sp.Rational(4,3)*(m**3-m))

    throat=geometry_scales(0.)
    h=geometry_scales(horizon["benchmark_Unruh_per_abs_q"]["reused_horizon_rho"])
    assert clock["child_geometry"]["asymptotic_Lambda"]=="9pi"
    assert clock["child_geometry"]["asymptotic_directional_H"]=="sqrt(3pi)"
    child_H2=3*pi
    child={"R_L":-36*pi,"Riemann_squared":216*pi*pi,
           "curvature_scale_squared":child_H2,"largest_section_magnitude":child_H2}
    kk_squared=(pi/match["proper_coordinate_zero_limit"]["proper_length"])**2
    scales=[]
    for name,g in (("horizon",h),("neck",throat),("child_limit",child)):
        scales.append({"location":name,"geometry":g,
                       "largest_section_over_cutoff_squared":g["largest_section_magnitude"]/cutoff**2,
                       "largest_section_over_matching_squared":g["largest_section_magnitude"]/coefficients["matching_cutoff"]**2,
                       "largest_section_over_first_compact_squared":g["largest_section_magnitude"]/kk_squared})
    child_radius=1/sqrt(child_H2)
    local_root=brentq(lambda a:local_endpoint(a,coefficients)["d_log_radius"],.2,.8,xtol=1e-13)
    radii=(child_radius,local_root,1.,2.,4.)
    w=weight()
    comparisons=[]
    for a in radii:
        full=endpoint(w,a);local=local_endpoint(a,coefficients)
        comparisons.append({"full":full,"local_through_a4_plus_exact_light":local,
                            "local_log_derivative_error":local["d_log_radius"]-full["d_log_radius"]})
    assert comparisons[0]["full"]["d_log_radius"]>0
    assert comparisons[0]["local_through_a4_plus_exact_light"]["d_log_radius"]<0
    assert abs(comparisons[1]["local_through_a4_plus_exact_light"]["d_log_radius"])<1e-10
    assert comparisons[1]["full"]["d_log_radius"]>.08

    controls=[]
    for a in (child_radius,local_root,1.,4.):
        compact=[endpoint(weight(n),a) for n in (18,26,34)]
        refined=endpoint(weight(34,144),a)
        sphere_refined=endpoint(w,a,last_level=compact[-1]["last_level"]+8)
        # Independent derivative of the total finite spectral sum; two steps
        # remove the leading central-difference error without changing inputs.
        step=1e-4
        derivatives=[]
        for e in (step,step/2):
            derivatives.append((endpoint(w,a*exp(e))["action"]-endpoint(w,a*exp(-e))["action"])/(2*e))
        independent=(4*derivatives[1]-derivatives[0])/3
        base=compact[-1]
        fd_error=independent-base["d_log_radius"]
        compact_delta={k:compact[-1][k]-compact[-2][k] for k in ("action","d_log_radius")}
        assert abs(fd_error)<3e-8*max(1,abs(base["d_log_radius"]))
        for k in compact_delta:
            assert abs(compact_delta[k])<3e-8*max(1,abs(base[k]))
            assert abs(refined[k]-base[k])<2e-11*max(1,abs(base[k]))
            assert abs(sphere_refined[k]-base[k])<2e-12*max(1,abs(base[k]))
        controls.append({"radius":a,"compact_18_26_34":[{k:r[k] for k in ("action","d_log_radius")} for r in compact],
                         "compact_last_difference":compact_delta,
                         "quadrature_112_to_144":{k:refined[k]-base[k] for k in compact_delta},
                         "sphere_extra_eight_levels":{k:sphere_refined[k]-base[k] for k in compact_delta},
                         "central_difference_steps":derivatives,"extrapolated_log_derivative":independent,
                         "independent_derivative_error":fd_error})

    # A direct sum of the known unwarped tower checks the NEW S4 composition.
    flat=weight(path=0.)
    flat_result=endpoint(flat,1.)
    y,deg=sphere_spectrum(1.,flat_result["last_level"])
    p=np.arange(flat.modes+1);multiplicity=np.array([1]+[2]*flat.modes)
    energies=y[:,None]+(p*pi/flat.interval)**2
    flat_action=.5*np.sum(deg[:,None]*multiplicity[None,:]*exp1(energies/cutoff**2))
    flat_slope=np.sum(deg[:,None]*multiplicity[None,:]*y[:,None]*np.exp(-energies/cutoff**2)/energies)
    assert abs(flat_action-flat_result["action"])<2e-11
    assert abs(flat_slope-flat_result["d_log_radius"])<2e-11
    return native({"schema":"NSC-SPECTRAL-ENDPOINT-v1","status":"full compact-weight endpoint changes the local source decision",
       "source_hashes":hashes(SOURCES),"input_hashes":hashes(INPUTS),"exact_residuals":exact,
       "conventions":{"sphere":"S4 radius a; eigenvalues +/-m/a, m>=2; squared multiplicity 4*(m³-m)/3",
                      "functional":"Gamma5=1/2 Tr4 h(D4²); no additional 5D copy factor",
                      "variation":"partial_log_a at fixed Lambda,ell,warp,state prescription; simultaneous scale stationarity not inferred",
                      "local_comparator":"exact canonical light determinant + matched local complement through a4, each included once",
                      "units":"hbar=c=1; sphere/profile lengths in the stored development units; physical scales not fitted",
                      "isotropic_action_response":"(partial_log_a Gamma)/(4 Vol(S4)); not separately added matter on top of the same induced gravity"},
       "reused_parameters":parameters,"reused_coefficients":coefficients,"first_compact_squared":kk_squared,
       "geometry_scale_comparison":scales,"endpoint_comparisons":comparisons,
       "local_stationary_control":{"radius":local_root,"bracket":[.2,.8],"a0_a2_radius_estimate":sqrt(6*coefficients["A_Dirac"]/coefficients["V_Dirac"]),
                                    "full_log_derivative":comparisons[1]["full"]["d_log_radius"],"is_full_stationary_solution":False},
       "numerical_controls":controls,"direct_flat_S4_composition":{"action_error":float(flat_action-flat_result["action"]),"log_derivative_error":float(flat_slope-flat_result["d_log_radius"])},
       "uncertainty":{"sphere_tails":"analytic continuum comparison bounds, evaluated in floating point",
                      "compact":"18/26/34 Galerkin changes measured; no rigorous retained-eigenvalue error bound claimed",
                      "derivative":"generalized-eigenvalue variation versus differentiated total action, with two step sizes",
                      "local_expansion":"the horizon ratios alone do not prove all derivative, gauge, state or nonlocal remainder bounds"},
       "state_and_scope":{"global_black_universe_continued_to_S4":False,"magnetic_flux_sector_extended_to_S4":False,
                          "de_Sitter_invariant_state_selected_by_recursion":False,"full_quantum_compensator_or_measure_derived":False,
                          "all_curved_radii_no_stationarity_theorem":False,"self_sourced_geometry_or_cosmological_flux_derived":False,
                          "partial_metric_variation_may_be_omitted_in_joint_solve":False},
       "decision":"use the full spectral source at child-scale curvature; do not accept the local-truncation stationary radius as a solution of the same determinant",
       "primary_source":{"url":"https://arxiv.org/html/gr-qc/9505009v1","equations":["3.51","3.52","3.54"],"use":"published sphere spectrum and both-sign multiplicity, not a new spectrum derivation"},
       "comparison":{"fields":"all","float_atol":3e-9,"float_rtol":3e-8,"exact":"strings, keys, structures, source/input hashes","exceptions":[]}})


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true");group.add_argument("--output",type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError("refusing to overwrite recorded evidence")
    expected=json.loads(OUTPUT.read_text()) if args.check else None
    if args.check:
        compare(expected["source_hashes"],hashes(SOURCES));compare(expected["input_hashes"],hashes(INPUTS))
    result=calculate()
    if args.check:compare(expected,result);print("spectral endpoint: all fields reproduced; prior generators were not run")
    elif args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open("x") as stream:json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write("\n")
        print(args.output)
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()
