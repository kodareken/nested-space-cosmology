"""Read-only constrained flat-patch assessment of the current local action.

This is an analytic nonperiodic patch, not new initial data or evolution.
The extra pole is assessed when the local truncation is treated as exact.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import subprocess
import time

import sympy as sp

from . import nsc_spherical_feedback_action as action

LAB=Path(__file__).resolve().parents[2]
ROOT=LAB.parent
SCHEMA="NSC-DISCOVERY-CURVATURE-SECTOR-v1"
OUTPUT=LAB/"results/development/nsc-discovery-curvature-sector-v1.json"
COEFFICIENTS=LAB/"results/development/nsc-subgap-history-response.json"
SOURCE_PATHS=(
    "lab/src/recursive_horizons/nsc_discovery_curvature_sector.py",
    "lab/scripts/assess_nsc_discovery_curvature_sector.py",
    "lab/tests/test_nsc_discovery_curvature_sector.py",
    "lab/docs/nsc-discovery-curvature-sector.md",
    "lab/src/recursive_horizons/nsc_spherical_feedback_action.py",
    "lab/src/recursive_horizons/nsc_discovery_tidal.py",
    "lab/docs/nsc-spherical-feedback-action.md",
    "lab/docs/nsc-spherical-conformal-gauge.md",
    "lab/docs/nsc-spherical-action.md",
    "lab/docs/nsc-curvature-eft.md",
)
EQUATIONS={
    "patch":"r=x>0; Q=L=1/x; beta=chi=p_Q=p_r=p_chi=Phi=flux=0",
    "metric":"ds^2=dt^2-dx^2-x^2 dOmega^2; N=rQ=1",
    "perturbations":"a=delta r; q=delta logQ; c=delta chi; g=8pi A; f=2alpha",
    "momenta":["delta p_chi=2f q_t", "delta p_r=-6g a_t-2gx q_t", "delta p_Q=2fx c_t-2gx^2 a_t"],
    "radial":"Box a=c/(6x)",
    "auxiliary":"Box q=-(c+4q)/(2x^2)",
    "metric_equation":"Box c=(gx/f)Box a+(2g/f)(q+a/x-a_x)-2c/x^2",
    "b":"b=a_x-a/x-q; T_c=c_x/x+c/x^2",
    "lapse":"delta C/(2gx^2)=b_x+b/x+(f/g)(c/x^3-c_x/x^2-c_xx/x)",
    "shift":"delta D/(2gx)=partial_t[b-(f/g)T_c]",
    "constraint_solution":"b=(f/g)T_c+K/x; partial_t K=partial_x K=0",
    "constrained_chi":"Box c+2c_x/x+4c/x^2-m_star^2 c=-2gK/(fx)",
    "mass":"m_star^2=g/(6f)=2pi A/(3alpha)=-A/(2C_W)",
    "Einstein_Coulomb":"c_E=12K/x; K=-delta Schwarzschild_mass",
    "physical_variable":"psi=(c-c_E)/x",
    "physical_equation":"psi_tt-psi_xx+(6/x^2-m_star^2)psi=0",
    "large_R_symbol_variables":["q","a/R","c/R^2","p_Q/R^3","p_r/R","p_chi"],
    "large_R_characteristic":"(lambda^2+k^2)^2 (lambda^2+k^2-m_star^2)",
    "physical_eigenvector":"(q,a/R,c/R^2,p_Q/R^3,p_r/R,p_chi)=(-3f cbar/g,f cbar/g,cbar,0,0,-6f^2 lambda cbar/g)",
    "proper_dispersion":"omega_proper^2=k_proper^2-m_star^2; d_tau=dt; d_ell=dx",
    "four_dimensional_action":"S_bulk=integral sqrt|g| [-A R-C_W C^2] modulo the recorded boundary terms",
    "TT_factor":"Box (Box+A/(2C_W)) = Box (Box-m_star^2)",
    "Ricci_Weyl_identity":"C^2=E4+2(Ricci^2-R^2/3)",
}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    return {path:sha256(ROOT/path) for path in SOURCE_PATHS}


def verify_producer(commit,hashes):
    """Read immutable blobs when root explicitly binds a production record."""
    resolved=subprocess.check_output(["git","rev-parse","--verify",commit+"^{commit}"],cwd=ROOT,text=True).strip()
    for path,digest in hashes.items():
        raw=subprocess.check_output(["git","show",resolved+":"+path],cwd=ROOT)
        if hashlib.sha256(raw).hexdigest()!=digest:
            raise ValueError("producer blob mismatch: "+path)
    return resolved


def large_radius_symbol(g,f,k):
    """Proper-flat limit of the six owned canonical geometric rates."""
    return sp.Matrix([
        [0,0,0,0,0,1/(2*f)],
        [0,0,0,0,-1/(6*g),-1/(6*f)],
        [0,0,0,1/(2*f),-1/(6*f),-g/(6*f**2)],
        [0,2*g*k**2,-2*f*k**2,0,0,0],
        [2*g*k**2,6*g*k**2,0,0,0,0],
        [-2*f*k**2,0,-f,0,0,0],
    ])


def symbolic_checks():
    """Derive from the actual first-order density, momenta and constraints."""
    t,x,eps=sp.symbols("t x eps",real=True)
    g,f=sp.symbols("g f",nonzero=True)
    r,u,chi=(sp.Function(name)(t,x) for name in ("r","u","chi"))
    Q=sp.exp(u);A=g/(8*sp.pi);CW=-3*f/(8*sp.pi)
    density=action.first_order_density(Q,sp.diff(Q,x),Q,sp.diff(Q,t),sp.diff(Q,x),
        r,sp.diff(r,t),sp.diff(r,x),chi,sp.diff(chi,t),sp.diff(chi,x),0,0,A,CW,0,0)
    density=sp.simplify(density)
    def euler(field):
        return sp.simplify(sp.diff(density,field)
            -sp.diff(sp.diff(density,sp.diff(field,t)),t)
            -sp.diff(sp.diff(density,sp.diff(field,x)),x))
    Er,Eu,Ec=(euler(field) for field in (r,u,chi))
    a,q,c=(sp.Function(name)(t,x) for name in ("a","q","c"))
    perturb={r:x+eps*a,u:-sp.log(x)+eps*q,chi:eps*c}
    box=lambda v:sp.diff(v,t,2)-sp.diff(v,x,2)
    linear=lambda v:sp.simplify(sp.diff(v.subs(perturb).doit(),eps).subs(eps,0))
    checks={
        "flat_Euler_background":all(sp.simplify(v.subs({r:x,u:-sp.log(x),chi:0}).doit())==0 for v in (Er,Eu,Ec)),
        "radial_Euler_linearization":sp.simplify(linear(Er)-2*g*(3*box(a)+x*box(q)+2*q/x))==0,
        "auxiliary_Euler_linearization":sp.simplify(linear(Ec)+2*f*(box(q)+(c+4*q)/(2*x**2)))==0,
        "metric_Euler_linearization":sp.simplify(linear(Eu)
            -(2*g*x*box(a)+4*g*(q+a/x-sp.diff(a,x))-4*f*c/x**2-2*f*box(c)))==0,
    }
    Qp=(1+eps*q)/x;rp=x+eps*a;cp=eps*c
    pQ=2*f*x*sp.diff(c,t)-2*g*x**2*sp.diff(a,t)
    pr=-6*g*sp.diff(a,t)-2*g*x*sp.diff(q,t)
    pc=2*f*sp.diff(q,t)
    momenta=action.momenta(Qp,Qp,sp.diff(Qp,t),sp.diff(Qp,x),rp,sp.diff(rp,t),sp.diff(rp,x),
        cp,sp.diff(cp,t),sp.diff(cp,x),0,0,A,CW)
    for name,owned,target in zip(("p_Q","p_r","p_chi"),momenta,(pQ,pr,pc)):
        checks["canonical_"+name]=sp.simplify(sp.diff(owned,eps).subs(eps,0)-target)==0
    C=action.lapse_constraint(eps*pQ,eps*pr,eps*pc,Qp,rp,cp,sp.diff(rp,x),sp.diff(rp,x,2),
        sp.diff(cp,x),sp.diff(cp,x,2),sp.diff(Qp,x),A,CW,0,0)
    D=action.shift_constraint(eps*sp.diff(pQ,x),eps*pr,eps*pc,Qp,sp.diff(rp,x),sp.diff(cp,x))
    checks["flat_constraint_background"]=sp.simplify(C.subs(eps,0))==0 and sp.simplify(D.subs(eps,0))==0
    b=sp.diff(a,x)-a/x-q;T=sp.diff(c,x)/x+c/x**2
    expectedC=2*g*x**2*(sp.diff(b,x)+b/x+(f/g)*(c/x**3-sp.diff(c,x)/x**2-sp.diff(c,x,2)/x))
    expectedD=2*g*x*sp.diff(b-f*T/g,t)
    checks["lapse_constraint_linearization"]=sp.simplify(sp.diff(C,eps).subs(eps,0)-expectedC)==0
    checks["shift_constraint_linearization"]=sp.simplify(sp.diff(D,eps).subs(eps,0)-expectedD)==0
    K=sp.symbols("K");m2=g/(6*f);psi=sp.Function("psi")(t,x)
    constrained=box(c)+2*sp.diff(c,x)/x+4*c/x**2-m2*c+2*g*K/(f*x)
    constrained_q=sp.diff(a,x)-a/x-f*T/g-K/x
    metric_wave=box(c)-(g*x/f)*box(a)-(2*g/f)*(q+a/x-sp.diff(a,x))+2*c/x**2
    checks["constraint_reduced_chi_equation"]=sp.simplify(metric_wave.subs(box(a),c/(6*x)).subs(q,constrained_q)-constrained)==0
    checks["Coulomb_constraint_solution"]=sp.simplify((sp.diff(b,x)+b/x
        +(f/g)*(c/x**3-sp.diff(c,x)/x**2-sp.diff(c,x,2)/x)).subs(q,sp.diff(a,x)-a/x-f*T/g-K/x).doit())==0
    checks["physical_radial_equation"]=sp.simplify(constrained.subs(c,x*psi).doit()/x
        -(box(psi)+(6/x**2-m2)*psi+2*g*K/(f*x**2)))==0
    checks["Einstein_Coulomb_separation"]=sp.simplify((box(psi)+(6/x**2-m2)*psi+2*g*K/(f*x**2)).subs(psi,12*K/x**2).doit())==0
    compatibility=(box(constrained_q)+(c+4*constrained_q)/(2*x**2)).doit()
    replacement={sp.diff(a,t,2,x):sp.diff(a,x,3)+sp.diff(c,x)/(6*x)-c/(6*x**2),
        sp.diff(a,t,2):sp.diff(a,x,2)+c/(6*x),
        sp.diff(c,t,2,x):sp.diff(c,x,3)+m2*sp.diff(c,x)-2*sp.diff(c,x,2)/x-2*sp.diff(c,x)/x**2+8*c/x**3+2*g*K/(f*x**2),
        sp.diff(c,t,2):sp.diff(c,x,2)+m2*c-2*sp.diff(c,x)/x-4*c/x**2-2*g*K/(f*x)}
    checks["auxiliary_equation_constraint_compatibility"]=sp.simplify(compatibility.subs(replacement,simultaneous=True))==0
    k,lam,cv=sp.symbols("k lambda cbar",nonzero=True)
    J=large_radius_symbol(g,f,k);polynomial=J.charpoly();z=polynomial.gen
    checks["large_R_characteristic_factorization"]=sp.simplify(polynomial.as_expr()-(z*z+k*k)**2*(z*z+k*k-m2))==0
    vector=sp.Matrix([-3*f*cv/g,f*cv/g,cv,0,0,-6*f**2*lam*cv/g])
    checks["constrained_massive_symbol_eigenvector"]=all(sp.simplify(v.subs(lam**2,m2-k**2))==0 for v in J*vector-lam*vector)
    checks["nongauge_chi_background"]=sp.diff(sp.Integer(0),x)==0 and sp.diff(sp.Integer(0),t)==0
    # This check compares the covariant spherical coefficient, not a guessed
    # signature continuation of a Euclidean heat-kernel coefficient.
    checks["Lorentz_Weyl_coefficient_map"]=sp.simplify(4*sp.pi*(-CW)/3-f/2)==0
    return checks


def proper_pole_scales(A,C_W):
    if not math.isfinite(A) or not math.isfinite(C_W) or A<=0 or C_W==0:
        raise ValueError("positive Einstein coefficient and nonzero Weyl chart required")
    alpha=-4*math.pi*C_W/3;mstar=-A/(2*C_W)
    return {"A":float(A),"C_W":float(C_W),"alpha":alpha,"f":2*alpha,
        "actual_Lorentz_Weyl_coefficient":-float(C_W),
        "m_star_square":mstar,"physical_spin2_pole_mass_square":-mstar,
        "long_wavelength_growth_rate":math.sqrt(mstar) if mstar>0 else None,
        "long_wavelength_growth_time":1/math.sqrt(mstar) if mstar>0 else None,
        "tachyonic_extra_pole":mstar>0,"opposite_Einstein_pole_residue":True,
        "proper_clock":"d_tau=dt because N=rQ=1 in this exact patch",
        "proper_radial_distance":"d_ell=dx in this exact patch",
        "units":"locked local-action proper units; no external clock calibration asserted"}


def assess(*,producer_commit=None):
    """Pure symbolic calculation; no trajectory, solver, or evidence write."""
    start=time.process_time();before=source_hashes();input_hash=sha256(COEFFICIENTS)
    producer=verify_producer(producer_commit,before) if producer_commit else None
    inputs=json.loads(COEFFICIENTS.read_text())["locked_inputs"]
    scales=proper_pole_scales(float(inputs["A"]),float(inputs["C_Weyl"]))
    checks=symbolic_checks()
    if before!=source_hashes() or input_hash!=sha256(COEFFICIENTS):
        raise RuntimeError("assessment owners or immutable coefficient input changed during calculation")
    return {"schema":SCHEMA,"status":"CONSTRAINED_FLAT_EXTRA_POLE_CONFIRMED" if all(checks.values()) else "UNSATISFIED_SYMBOLIC_IDENTITY",
        "equations":EQUATIONS,"checks":checks,"proper_scales":scales,
        "convention":{"signature":"(+---)",
            "Riemann":"R^a_bcd=d_c Gamma^a_db-d_d Gamma^a_cb+Gamma^a_ce Gamma^e_db-Gamma^a_de Gamma^e_cb",
            "Ricci":"Ricci_bd=R^a_bad", "physical_extra_mass_square":"mu2=A/(2C_W)",
            "ghost_scope":"opposite residue relative to the healthy Einstein pole in this local truncation"},
        "scope":{"local_truncated_action_treated_as_exact":True,"analytic_nonperiodic_flat_patch":True,
            "periodic_whole_domain_initial_data":False,"source_and_flux_zero_in_patch":True,
            "code_declared_action_sign_mismatch_found":False,"all_source_instability_proved":False,
            "full_nonlocal_spectral_action_mode_health_assessed":False,"EFT_pole_validity_established":False,
            "order_reduction_implemented":False,"NSC_verdict":False,"trajectory_performed":False},
        "EFT_connection":{"owner":"lab/docs/nsc-curvature-eft.md",
            "meaning":"that first-order perturbative branch retains leading initial data rather than the extra modes of a resummed fourth-order truncation",
            "not_a_proof_of_full_spectral_mode_health":True},
        "source_hashes":before,"input_hashes":{str(COEFFICIENTS.relative_to(ROOT)):input_hash},
        "producing_commit":producer,"sympy_version":sp.__version__,
        "cpu_seconds":time.process_time()-start,"evidence_written":False}


def write_record(record,path=OUTPUT):
    """Root's explicit immutable create; refuse unfrozen or failed assessments."""
    target=Path(path).expanduser().resolve()
    if target.exists():
        raise FileExistsError("curvature-sector record already exists; use a successor")
    if record.get("producing_commit") is None:
        raise ValueError("production record needs its frozen producing commit")
    if record.get("schema")!=SCHEMA or not all(record["checks"].values()):
        raise ValueError("unresolved symbolic assessment cannot be published as confirmed")
    verify_producer(record["producing_commit"],record["source_hashes"])
    for path,digest in record["input_hashes"].items():
        if sha256(ROOT/path)!=digest:
            raise ValueError("assessment input changed before creation")
    result=dict(record,evidence_written=True)
    target.parent.mkdir(parents=True,exist_ok=True)
    with target.open("x") as stream:
        json.dump(result,stream,indent=2,allow_nan=False);stream.write("\n")
    return result


def check_record(path=OUTPUT):
    """Authenticate the immutable producer and recompute symbolic statements."""
    target=Path(path).expanduser().resolve();old=json.loads(target.read_text())
    if old.get("schema")!=SCHEMA:
        raise ValueError("curvature-sector schema mismatch")
    verify_producer(old["producing_commit"],old["source_hashes"])
    if source_hashes()!=old["source_hashes"]:
        raise ValueError("live assessment differs from frozen producer; pinned replay required")
    for path,digest in old["input_hashes"].items():
        if sha256(ROOT/path)!=digest:
            raise ValueError("curvature-sector immutable input mismatch")
    fresh=assess()
    for key in ("status","equations","checks","proper_scales","convention","scope","EFT_connection"):
        if fresh[key]!=old[key]:
            raise ValueError("curvature-sector symbolic replay mismatch: "+key)
    return {"ok":True,"bytes_written":0,"trajectory_performed":False,"symbolic_recomputed":True,
            "record_sha256":sha256(target),"producing_commit":old["producing_commit"]}
