#!/usr/bin/env python3
"""Reproduce the finite-radius canonical massless angular source only."""
import argparse
import json
from math import pi
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from recursive_horizons.nsc_angular_stress import (
    angular_reference,physical_source,smooth_integration_window,profile_jets,
    adiabatic_coefficients,curvature_squared_tensor,static_cylinder_density,
)
from recursive_horizons.nsc_covariant_operator import cylinder_zero_winding_energy
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare,native
from check_nsc_vacuum_charge_matching import hashes

OUTPUT=ROOT/"results/development/angular-stress.json"
SOURCES=("scripts/check_nsc_angular_stress.py","src/recursive_horizons/nsc_angular_stress.py",
         "src/recursive_horizons/_angular_transport.cpp","docs/nsc-angular-stress.md","pyproject.toml",
         "src/recursive_horizons/nsc_covariant_operator.py","scripts/check_nsc_compact_boundary_action.py",
         "scripts/check_nsc_compact_casimir.py","scripts/check_nsc_vacuum_charge_matching.py")
INPUTS=("results/development/massless-reference.json","results/nsc-9-covariant-source.json")
WINDOWS=((4096,6144),(4096,8192),(5120,8192),(6144,8192))
SPECS={
    "base":dict(angular_max=1024,momentum_max=4096,momentum_points=2048,steps=20000,backend="native"),
    "angular_tail":dict(angular_start=1025,angular_max=2048,momentum_max=4096,momentum_points=2048,steps=40000,backend="native"),
    "momentum_tail":dict(angular_max=2048,momentum_min=4096,momentum_max=8192,momentum_points=1024,steps=40000,backend="native",momentum_windows=WINDOWS),
    "momentum_coarse":dict(angular_max=2048,momentum_min=4096,momentum_max=8192,momentum_points=512,steps=40000,backend="native",momentum_windows=WINDOWS),
    "time_coarse":dict(angular_max=512,momentum_max=2048,momentum_points=1536,steps=8000,backend="native"),
    "time_fine":dict(angular_max=512,momentum_max=2048,momentum_points=1536,steps=16000,backend="native"),
}
COMPONENTS=("rho","p_parallel","p_sphere","radial_null")


def source_from_terms(base,energy,pressure,angular_window):
    weights=smooth_integration_window(np.arange(1,len(energy)+1),*angular_window)
    row=dict(base)
    row["bar_rho"]=float(base["static_density"]+base["restored_curvature_rho"]+weights@energy)
    row["bar_p_parallel"]=float(-base["static_density"]+base["restored_curvature_p_parallel"]+weights@pressure)
    return {"bar_rho":row["bar_rho"],"bar_p_parallel":row["bar_p_parallel"],"physical":physical_source(row)}


def assemble(base,angular_tail,momentum_tail):
    rows=[]
    for b,t,h in zip(base["rows"],angular_tail["rows"],momentum_tail["rows"]):
        assert b["q"]==t["q"]==h["q"]
        energy=np.r_[b["angular_energy_remainders"],t["angular_energy_remainders"]]
        pressure=np.r_[b["angular_pressure_remainders"],t["angular_pressure_remainders"]]
        assert len(energy)==len(pressure)==2048
        variations=[]
        for window in h["windowed_momentum_remainders"]:
            e=energy+window["angular_energy_remainders"]
            p=pressure+window["angular_pressure_remainders"]
            for angular in ((512,1024),(1024,1536),(1024,2048),(1536,2048)):
                variations.append({"angular_window":list(angular),
                    "momentum_window":[window["begin"],window["end"]],
                    **source_from_terms(b,e,p,angular)})
        refined=[v for v in variations if v["angular_window"][0]>=1024]
        spread={key:max(v["physical"][key] for v in refined)-min(v["physical"][key] for v in refined) for key in COMPONENTS}
        rows.append({"selected":variations[-1],"resolved_window_spread":spread,"window_variations":variations})
    return rows


def validate_cached_run(run,spec):
    """Development output may prepare a candidate; --check always reruns it."""
    defaults={"angular_start":1,"momentum_min":0.,"momentum_windows":[]}
    for key,value in spec.items():
        actual=run["parameters"].get(key,defaults.get(key))
        assert native(actual)==native(value),(key,actual,value)
    assert [r["q"] for r in run["rows"]]==[.1,.2,.5,1.,pi/2]


def calculate(completed_development_runs=None):
    reference,_=[authenticated_record(p) for p in INPUTS]
    runs={}
    for name,spec in SPECS.items():
        if completed_development_runs is not None and name in completed_development_runs:
            run=completed_development_runs[name];validate_cached_run(run,spec)
        else:
            print("Evaluating angular-source block: "+name,file=sys.stderr,flush=True)
            run=angular_reference(**spec)
        runs[name]=run
    rows=assemble(runs["base"],runs["angular_tail"],runs["momentum_tail"])
    coarse=assemble(runs["base"],runs["angular_tail"],runs["momentum_coarse"])
    momentum_control=[]
    for a,b in zip(coarse,rows):
        momentum_control.append({"q":b["selected"]["physical"]["q"],
            "difference":{key:b["selected"]["physical"][key]-a["selected"]["physical"][key] for key in COMPONENTS}})
    time_control=[]
    for a,b in zip(runs["time_coarse"]["rows"],runs["time_fine"]["rows"]):
        x=source_from_terms(a,np.array(a["angular_energy_remainders"]),np.array(a["angular_pressure_remainders"]),(256,512))
        y=source_from_terms(b,np.array(b["angular_energy_remainders"]),np.array(b["angular_pressure_remainders"]),(256,512))
        time_control.append({"q":b["q"],"difference":{key:y["physical"][key]-x["physical"][key] for key in COMPONENTS}})
    assert max(abs(v) for row in momentum_control for v in row["difference"].values())<2e-6
    assert max(abs(v) for row in time_control for v in row["difference"].values())<2e-4
    assert max(v for row in rows for v in row["resolved_window_spread"].values())<.002

    exact={}
    def zero(name,value):
        assert sp.simplify(value)==0,name
        exact[name]="0"
    r,mu,zp,gamma=sp.symbols("r mu zeta_prime_minus3 EulerGamma",positive=True)
    V=-(2*zp+(sp.log(mu**2*r**2)+1-gamma)/120)/(sp.pi*r*r)
    rho=V/(4*sp.pi*r*r);p_sphere=-sp.diff(V,r)/(8*sp.pi*r)
    zero("joint_cylinder_normalization_trace",2*rho-2*p_sphere+1/(240*sp.pi**2*r**4))
    A4=sp.Symbol("a4_2");s=sp.Symbol("s")
    zeta_series=A4*r*r/sp.pi*(sp.Rational(1,2)+s*(gamma+sp.log(r)))
    finite=(sp.diff(zeta_series,s).subs(s,0)+sp.log(mu**2)*zeta_series.subs(s,0)-gamma*zeta_series.subs(s,0))/2
    zero("joint_angular_harmonic_finite_part",finite-A4*r*r/(2*sp.pi)*(sp.log(mu*r)+gamma/2))
    coefficient=4*(1-3*sp.pi/2)
    assert coefficient.is_negative
    exact["positive_Einstein_coefficient_requires_negative_neck_null"]="4*(1-3*pi/2)<0"

    static=[]
    for cutoff in (4.,8.,16.):
        raw=cylinder_zero_winding_energy(1.,1.,cutoff,int(10*cutoff))
        sub=cutoff**4/(4*pi)-cutoff**2/(12*pi)-np.log(cutoff**2)/(120*pi)
        measured=(raw-sub)/(4*pi)
        static.append({"cutoff":cutoff,"subtracted_density":measured,"exact_density":static_cylinder_density(),"difference":measured-static_cylinder_density()})
    assert abs(static[-1]["difference"])<2.1e-7
    adiabatic=[]
    for q in (.2,1.,pi/2):
        jets=profile_jets(q);he,hz=curvature_squared_tensor(jets)
        value=[quad(lambda p:adiabatic_coefficients(p,1.,jets)[i]/pi,0,np.inf,epsabs=1e-10)[0] for i in (2,3)]
        target=[he/(480*pi),-hz/(480*pi)]
        assert max(abs(x-y) for x,y in zip(value,target))<2e-9
        adiabatic.append({"q":q,"integrated_E4_P4":value,"heat_target":target,"difference":[x-y for x,y in zip(value,target)]})

    backend_args=dict(points=(.2,pi/2),angular_max=12,momentum_points=64,momentum_max=128.,steps=1000)
    python_run=angular_reference(**backend_args);native_run=angular_reference(**backend_args,backend="native")
    backend_error=max(abs(a[key]-b[key]) for a,b in zip(python_run["rows"],native_run["rows"]) for key in ("bar_rho","bar_p_parallel"))
    assert backend_error<2e-10
    centers=(.2,.5,1.,pi/2);eps=.00025
    points=sorted({q+d for q in centers for d in (-2*eps,-eps,0.,eps,2*eps)})
    ward_run=angular_reference(points=points,angular_max=24,momentum_points=256,momentum_max=256,steps=12000,backend="native")
    physical=[physical_source(row) for row in ward_run["rows"]]
    ward=[]
    for index,q in enumerate(centers):
        b=ward_run["rows"][5*index:5*index+5];g=physical[5*index:5*index+5]
        db=(b[0]["bar_rho"]-8*b[1]["bar_rho"]+8*b[3]["bar_rho"]-b[4]["bar_rho"])/(12*eps)
        dg=(g[0]["rho"]-8*g[1]["rho"]+8*g[3]["rho"]-g[4]["rho"])/(12*eps)
        w,w1=b[2]["W_jets"][:2];m=g[2]
        br=db+w1/(2*w)*(b[2]["bar_rho"]+b[2]["bar_p_parallel"])
        gr=-np.sqrt(w)*np.sin(q)*dg+m["H_parallel"]*(m["rho"]+m["p_parallel"])+2*m["H_sphere"]*(m["rho"]+m["p_sphere"])
        assert max(abs(br),abs(gr))<2e-6
        ward.append({"q":q,"barred_residual":br,"physical_residual":gr})
    max_norm=max(row["unitarity_defect"] for run in runs.values() for row in run["rows"])
    assert max_norm<1e-10
    neck=rows[-1]["selected"]["physical"]
    return native({"schema":"NSC-ANGULAR-STRESS-v1","status":"computed finite-radius canonical massless reference tensor with observed numerical convergence",
      "source_hashes":hashes(SOURCES),"input_hashes":hashes(INPUTS),"exact_checks":exact,
      "reference":{"width":.1,"state":"the unchanged smooth auxiliary ultrastatic reference from NSC-MASSLESS-REFERENCE-v1",
                   "field":"one untwisted massless 4D Dirac, all S2 angular channels in the convergent sum",
                   "parent_state_matched":False,"massive_compact_tower_included":False,
                   "units":"hbar=c=L_throat=1; normalization mu=1; no observational parameter fitted"},
      "method":{"angular_degeneracy":"4*kappa, kappa>=1; 2D spin counted inside each Dirac channel",
                "subtraction":"2D orders0,2 plus fourth-order convergence subtraction; joint4D proper-time finite part restored",
                "harmonic_finite_part":"ln(mu*R)+EulerGamma/2",
                "sphere_pressure":"exact4D trace after the two independent base-metric responses",
                "physical_map":"integrated anomaly with independent lapse and longitudinal-metric variations",
                "windows":"C-infinity windows on the already convergent remainder, fixed in comoving labels; no physical cutoff substituted",
                "selected_angular_window":[1536,2048],"selected_momentum_window":[6144,8192],
                "native":"C++17 fourth-order unitary transport; Python control retained"},
      "run_specifications":SPECS,"source_rows":rows,
      "numerical_checks":{"momentum_512_to_1024_nodes":momentum_control,"time_8000_to_16000_steps_probe":time_control,
        "static_proper_time_normalization":static,"fourth_order_adiabatic_integrals":adiabatic,
        "native_python_difference":backend_error,"maximum_unitarity_defect":max_norm,"independent_density_derivative_Ward":ward,
        "uncertainty_scope":"window spreads and resolution differences are measured sensitivities, not rigorous infinite-tail error bounds"},
      "neck_source_budget":{"reference_tensor":{key:neck[key] for key in COMPONENTS},
        "reference_window_spread":rows[-1]["resolved_window_spread"],
        "geometric_null_requirement":"T_total,kk * L_throat^4 = 4*a_EH*(1-3*pi/2), a_EH=A_EH*L_throat²>0",
        "coefficient_of_a_EH":float(coefficient),"reference_alone_has_required_null_sign":bool(neck["radial_null"]<0),
        "required_remaining_null":"4*a_EH*(1-3*pi/2)-reference_null; this is a residual requirement, not an inserted source",
        "Einstein_coefficient_fitted":False,"complete_source_budget_closed":False},
      "remaining_owners":["parent/reference covariance and boundary flux","actual finite-Lambda CTP conversion/complement",
                          "massive compact and gauge contributions from the same functional","independent complete metric/scale equations"],
      "scope":{"full_finite_cutoff_Dirac_stress_claimed":False,"auxiliary_reference_tuned_for_a_neck":False,
               "asymptotic_result_substituted_at_neck":False,"self_sourcing_or_observational_prediction_derived":False,
               "new_general_adiabatic_or_conformal_theorem_claimed":False},
      "comparison":{"fields":"all","float_atol":3e-9,"float_rtol":3e-8,"exact":"keys, strings, structures, source/input hashes and scope","exceptions":[]}})


def main():
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true");group.add_argument("--output",type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():raise FileExistsError("refusing to overwrite recorded evidence")
    expected=json.loads(OUTPUT.read_text()) if args.check else None
    if args.check:compare(expected["source_hashes"],hashes(SOURCES));compare(expected["input_hashes"],hashes(INPUTS))
    result=calculate()
    if args.check:compare(expected,result);print("angular stress: every source/control field reproduced; previous scientific generators were not run")
    elif args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open("x") as stream:json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False);stream.write("\n")
        print(args.output)
    else:print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":main()
