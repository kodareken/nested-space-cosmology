#!/usr/bin/env python3
"""Read-only leading static virial; explicit immutable record creation only."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import sympy as sp
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_leading_einstein as leading
from recursive_horizons import nsc_discovery_dynamic_preparation as preparation
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons import nsc_spherical_coupling as coupling

ROOT=leading.ROOT
SCHEMA="NSC-DISCOVERY-LEADING-VIRIAL-v1"
LEADING=leading.OUTPUT/"manifest.json"
LEADING_PRODUCER="e6e784b32213ea2368b7fc863c37c72fba6aa4e3"
OUTPUT=leading.LAB/"results/development/nsc-discovery-leading-virial-v1.json"
OWNED=("lab/scripts/assess_nsc_discovery_leading_virial.py",
       "lab/tests/test_nsc_discovery_leading_virial.py",
       "lab/docs/nsc-discovery-leading-virial.md")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def authenticate_leading():
    """Bind the actual e6 leading evaluator/data; no trajectory replay or Git."""
    record=json.loads(LEADING.read_text())
    if record.get("schema")!=leading.SCHEMA or record.get("producing_commit")!=LEADING_PRODUCER:
        raise ValueError("leading virial requires the declared e6 leading/data producer")
    runtime={p:digest for p,digest in record["source_hashes"].items() if p.startswith("lab/src/")}
    for path,digest in {**runtime,**record["input_hashes"]}.items():
        if sha256(ROOT/path)!=digest:
            raise ValueError("actual leading source/input authentication failed: "+path)
    return {"producing_commit":record["producing_commit"],"manifest_sha256":sha256(LEADING),
        "runtime_source_hashes":runtime,"input_hashes":record["input_hashes"],
        "trajectory_replayed":False,"scope":"current e6 evaluator and immutable preparation inputs"}


def verify_producer(commit,hashes):
    """Root's explicit production freeze; reads immutable blobs only."""
    resolved=subprocess.check_output(["git","rev-parse","--verify",commit+"^{commit}"],cwd=ROOT,text=True).strip()
    for path,digest in hashes.items():
        value=subprocess.check_output(["git","show",resolved+":"+path],cwd=ROOT)
        if hashlib.sha256(value).hexdigest()!=digest:
            raise ValueError("assessment producer blob mismatch: "+path)
    return resolved


def symbolic_checks():
    x=sp.symbols("x",real=True);g,mag=sp.symbols("g mag",positive=True)
    r,Q=(sp.Function(name)(x) for name in ("r","Q"))
    F=-g*r*r/2;V=g*r*r-mag;Z=-3*g;u=sp.log(Q)
    C=Z*sp.diff(r,x)**2/Q-Q*V-2*sp.diff(sp.diff(F,x)/Q,x)
    pr=2*Z*sp.diff(r,x,2)+2*g*Q**2*r-2*g*r*sp.diff(u,x,2)
    HE=-g*(3*sp.diff(r,x)**2+Q**2*r**2+2*r*sp.diff(r,x)*sp.diff(u,x))
    radial_gradient=-(sp.diff(HE,r)-sp.diff(sp.diff(HE,sp.diff(r,x)),x))
    boundary=r*sp.diff(r,x)+r*r*sp.diff(u,x)
    # Matter is added afterwards: Q rho. It has no fixed-Q r derivative.
    residual=sp.simplify(Q*C+pr*r/2-mag*Q**2+g*sp.diff(boundary,x))
    BR={r:1,Q:1/sp.cos(x),mag:g}
    return {"leading_radial_Euler_derivative":sp.simplify(radial_gradient-pr)==0,
        "continuum_boundary_identity":residual==0,
        "BR_electrovac_constraint":sp.simplify(C.subs(BR).doit())==0,
        "BR_electrovac_radial_force":sp.simplify(pr.subs(BR).doit())==0,
        "BR_boundary_charge":sp.simplify(boundary.subs(BR).doit()-sp.tan(x))==0,
        "BR_axial_curvature":sp.simplify(2*sp.diff(-sp.log(sp.cos(x)),x,2)*sp.cos(x)**2-2)==0}


def finite_control():
    """Actual NF32 positive-eigenstate SBP probe, deliberately not solved IC."""
    with threadpool_limits(limits=1):
        pair=preparation.make_pair(32);grid=pair.grid
        Q=.6+.03*np.cos(2*np.pi*grid.xi_g/grid.length)
        state,spectral=preparation.spectral_state(pair,Q)
        state.r=1.3+.08*np.sin(2*np.pi*grid.xi_g/grid.length)
        fine=galerkin.prolong_state(grid,state)
        system=galerkin.active_fine_system(grid,fine)
        source=coupling.source_from_columns(system,fine)
        C,D=leading.constraint_arrays(grid,fine,system,source)
        rate=leading.fine_rates(system,fine,source)
        _,_,mag=leading.coefficients(system)
        dx=grid.dx_q
        E=float(np.sum(fine.Q*source["force_L"]))
        Emag=float(mag*dx*np.sum(fine.Q**2))
        radial=float(-.5*dx*np.dot(fine.r,rate[3]))
        constraint=float(-dx*np.dot(fine.Q,C))
        coarse_radial=float(-.5*grid.dx_g*np.dot(state.r,galerkin.pull_geometry(grid,rate[3])))
        coarse_constraint=float(-grid.dx_g*np.dot(state.Q,galerkin.pull_geometry(grid,C)))
        frame=np.vstack((state.phi0,state.phi1))
        gram=frame.conj().T@frame
        weighted=np.sqrt(pair.weights)[:,None]*gram*np.sqrt(pair.weights)[None,:]
        return {"label":"constructor-prepared NF32 positive-spectrum eigenstate SBP probe; not source-solved initial data",
            "nf":grid.nf,"ng":grid.ng,"nq":grid.nq,"period":grid.length,"gauge":"conformal",
            "Q_declaration":"0.6+0.03*cos(2pi*x/8)","r_declaration":"1.3+0.08*sin(2pi*x/8)",
            "zero_momenta":True,"chi_or_p_chi_used_by_leading_evaluator":False,
            "initial_constraints_solved":False,"static_criticality_imposed":False,"trajectory_performed":False,
            "occupations":pair.weights.tolist(),"source_trace":float(np.sum(pair.weights)),
            "multiplicity":system.multiplicity,"eigenvalues":spectral["eigenvalues"].tolist(),
            "gram_max":spectral["gram_max"],"commutator_max":spectral["commutator_max"],
            "CAR_eigenvalues":np.linalg.eigvalsh(weighted).tolist(),
            "source_energy":E,"spectral_energy":float(system.multiplicity*np.dot(pair.weights,spectral["eigenvalues"])),
            "magnetic_energy":Emag,"lhs":-E-Emag,"rhs_radial":radial,"rhs_constraint":constraint,
            "rhs":radial+constraint,"SBP_identity_indicator":-E-Emag-radial-constraint,
            "coarse_rhs":coarse_radial+coarse_constraint,
            "adjoint_pairing_indicator":radial+constraint-coarse_radial-coarse_constraint,
            "full_lapse_residual_max":float(np.max(np.abs(C))),"full_shift_residual_max":float(np.max(np.abs(D))),
            "continuum_or_roundoff_error_certificate":False}


def assess(*,producer_commit=None):
    start=time.process_time();binding=authenticate_leading()
    hashes={path:sha256(ROOT/path) for path in sorted(set(OWNED)|set(binding["runtime_source_hashes"]))}
    producer=verify_producer(producer_commit,hashes) if producer_commit else None
    inputs={str(LEADING.relative_to(ROOT)):sha256(LEADING),**binding["input_hashes"],
        "lab/results/development/nsc-subgap-history-response.json":sha256(leading.LAB/"results/development/nsc-subgap-history-response.json")}
    checks=symbolic_checks();probe=finite_control()
    if hashes!={path:sha256(ROOT/path) for path in hashes} or any(sha256(ROOT/p)!=digest for p,digest in inputs.items()):
        raise RuntimeError("assessment source or producing input changed")
    return {"schema":SCHEMA,"status":"STATIC_POSITIVE_ENERGY_PERIODIC_CLASS_EXCLUDED" if all(checks.values()) else "UNRESOLVED_SYMBOLIC_IDENTITY",
        "symbolic_checks":checks,"finite_control":probe,"leading_binding":binding,
        "equations":{"periodic":"-<Q rho>-mag<Q^2>=-0.5<r p_r_dot>-<Q C>",
            "positive_spectrum":"<Q rho>=M Tr(Cov H_Q)>0; M=4kappa once",
            "boundary":"E_D+mag integral Q^2=integral QC+0.5 integral r p_r_dot+8pi A [r r_x+r^2(logQ)_x]_a^b",
            "static_boundary":"E_D+mag integral Q^2=8pi A [r r_x+r^2(logQ)_x]_a^b",
            "BR":"r=1,Q=L=sec x,R_h=2,mag=8pi A; source=0; boundary=8pi A[tan x]"},
        "assumptions":{"positive_r_Q":True,"static_metric_beta_zero_L_equals_Q":True,
            "static_velocities_force_zero_leading_momenta":True,"leading_Einstein_action":True,
            "massless_four_dimensional_Dirac_trace":True,"Maxwell_tracefree":True,
            "angular_kappa_over_r_is_not_fundamental_4D_mass":True,
            "source_energy_positive_or_nonzero_magnetic_energy":True,"periodic_no_boundary_term":True},
        "scope":{"STATIC_class_only":True,"finite_SBP_identity":True,"continuum_boundary_formula_separate":True,
            "dynamic_or_oscillatory_regimes_excluded":False,"all_sources_excluded":False,"NSC_verdict":False,
            "BR_matter_solution_constructed":False,"exterior_pressure_or_wall_added":False,
            "source_solved_trajectory_generated":False,"numerical_error_certificate":False},
        "missing_periodic_term":"the signed boundary traction/gradient charge 8pi A[rr_x+r^2(logQ)_x]; static Dirichlet data alone do not supply physical exterior matching",
        "source_hashes":hashes,"input_hashes":inputs,"producing_commit":producer,
        "cpu_seconds":time.process_time()-start,"evidence_written":False}


def write_record(record,path=OUTPUT):
    target=Path(path).expanduser().resolve()
    if target.exists():
        raise FileExistsError("leading virial record already exists")
    if not record.get("producing_commit"):
        raise ValueError("creation requires a frozen producing commit")
    if not all(record["symbolic_checks"].values()):
        raise ValueError("unresolved symbolic identity")
    verify_producer(record["producing_commit"],record["source_hashes"])
    for key,digest in record["input_hashes"].items():
        if sha256(ROOT/key)!=digest:
            raise ValueError("producing input changed")
    target.parent.mkdir(parents=True,exist_ok=True)
    result=dict(record,evidence_written=True)
    with target.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write("\n")
    return result


def check_record(path=OUTPUT):
    target=Path(path).expanduser().resolve();record=json.loads(target.read_text())
    if record.get("schema")!=SCHEMA:
        raise ValueError("leading virial schema mismatch")
    verify_producer(record["producing_commit"],record["source_hashes"])
    for key,digest in {**record["source_hashes"],**record["input_hashes"]}.items():
        if sha256(ROOT/key)!=digest:
            raise ValueError("leading virial source/input mismatch")
    fresh=assess()
    for key in ("status","symbolic_checks","equations","assumptions","scope","missing_periodic_term","leading_binding"):
        if record[key]!=fresh[key]:
            raise ValueError("leading virial symbolic replay mismatch: "+key)
    for key,value in fresh["finite_control"].items():
        previous=record["finite_control"][key]
        if isinstance(value,float):
            if not np.isclose(value,previous,rtol=1e-10,atol=1e-10):
                raise ValueError("leading virial finite replay mismatch: "+key)
        elif isinstance(value,list) and value and isinstance(value[0],float):
            if not np.allclose(value,previous,rtol=1e-10,atol=1e-10):
                raise ValueError("leading virial array replay mismatch: "+key)
        elif value!=previous:
            raise ValueError("leading virial finite metadata mismatch: "+key)
    return {"ok":True,"bytes_written":0,"trajectory_recomputed":False,"record_sha256":sha256(target)}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--producer-commit")
    parser.add_argument("--write",type=Path,nargs="?",const=OUTPUT)
    parser.add_argument("--check",type=Path,nargs="?",const=OUTPUT)
    args=parser.parse_args(argv)
    if args.check is not None:
        if args.write is not None or args.producer_commit:
            parser.error("check cannot write or change its producer")
        result=check_record(args.check)
    else:
        if args.write is not None and not args.producer_commit:
            parser.error("creation requires producer-commit")
        if args.write is not None and args.write.exists():
            raise FileExistsError("record already exists")
        result=assess(producer_commit=args.producer_commit)
        if args.write is not None:
            result=write_record(result,args.write)
    print(json.dumps(result,indent=2,allow_nan=False))
    return 0


if __name__=="__main__":
    raise SystemExit(main())
