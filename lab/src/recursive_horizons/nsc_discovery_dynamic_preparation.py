"""Positive general-Q, chi=0 Cauchy preparation using the owned finite action.

The source is selected once per initial geometry. No evolution, static force
condition, current deletion, or evidence write occurs in the constructors.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import io
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigh, solve
from threadpoolctl import threadpool_limits

from . import nsc_nested_parent_child as nested
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin

LAB = Path(__file__).resolve().parents[2]
ROOT = LAB.parent
SCHEMA = "NSC-DISCOVERY-DYNAMIC-PREPARATION-v1"
BALANCE_INPUT = LAB / "results/development/nsc-discovery-stationary-balance-seed-v2.json"
DEFAULT_OUTPUT = LAB / "results/development/nsc-discovery-dynamic-preparation-v1"
WEIGHTS = np.array([.75, .75, .5, .5, .25, .25])
OWNERS = ("lab/src/recursive_horizons/nsc_discovery_dynamic_preparation.py",
          "lab/scripts/derive_nsc_discovery_dynamic_preparation.py",
          "lab/tests/test_nsc_discovery_dynamic_preparation.py",
          "lab/docs/nsc-discovery-dynamic-preparation.md",
          "lab/src/recursive_horizons/nsc_nested_parent_child.py",
          "lab/src/recursive_horizons/nsc_spherical_coupling.py",
          "lab/src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
          "lab/src/recursive_horizons/nsc_spherical_feedback_action.py")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def maximum(value):
    return float(np.max(np.abs(value)))


def make_pair(nf=128):
    columns = (np.eye(nf, 6, dtype=complex), np.zeros((nf, 6), complex))
    return nested.build_pair(nf, coarse_modes=1, child_details=2,
                             source_layout="override", columns_override=columns,
                             occupations=WEIGHTS)


def declared_q(grid, *, harmonic=8, lift=1.1, pattern_fraction=1.):
    """Use saved balance geometry only as a declaration of Q mean and shape.

    Finite cosines avoid exponentiation aliases. Tests may explicitly override
    the harmonic; that is a declared test preparation, not saved n=8 evidence.
    """
    if not np.isfinite(lift) or lift <= 0 or not 0 <= pattern_fraction <= 1:
        raise ValueError("positive lift and pattern fraction in [0,1] required")
    if int(harmonic) != harmonic or harmonic < 1 or harmonic > (grid.ng-1)//2:
        raise ValueError("pattern harmonic must be resolved")
    saved = json.loads(BALANCE_INPUT.read_text())
    payload = BALANCE_INPUT.with_name(saved["payload"])
    if sha256(payload) != saved["payload_sha256"]:
        raise ValueError("saved balance geometry payload hash mismatch")
    row = saved["trials"][saved["selected_trial_index"]]
    prefix = row["array_prefix"]
    with np.load(payload, allow_pickle=False) as arrays:
        original = arrays[prefix+"Q"]
        fourier = np.fft.fft(original/np.mean(original))/len(original)
    count = min(7, (grid.ng-1)//(2*harmonic))
    coefficients = 2*fourier[np.arange(1,count+1)*row["harmonic"]].real
    mean = lift*saved["measurements"]["Q_mean"]
    cosine = np.cos(2*np.pi*grid.xi_g[:,None]*harmonic*np.arange(1,count+1)/grid.length)
    Q = mean*(1+pattern_fraction*cosine@coefficients)
    if np.min(galerkin.prolong_geometry(grid, Q)) <= 0:
        raise ValueError("declared Q leaves the fine positive chart")
    metadata = {"mean_Q":float(mean), "harmonic":int(harmonic),
                "lift":float(lift), "pattern_fraction":float(pattern_fraction),
                "cosine_coefficients":(pattern_fraction*coefficients).tolist(),
                "Q_is_declared_initial_input":True,
                "saved_geometry_is_holding_evidence":False}
    return Q, metadata


def spectral_state(pair, Q):
    """Real positive standing eigenmodes of actual U_f† H(Q) U_f, eta=1."""
    grid = pair.grid
    Q = np.asarray(Q, float)
    if Q.shape != (grid.ng,) or not np.isfinite(Q).all():
        raise ValueError("Q must be a finite geometry vector")
    zeros = np.zeros(grid.ng)
    columns = np.zeros((grid.nf,6), complex)
    state = coupling.CauchyState(Q.copy(), np.ones(grid.ng), zeros.copy(),
                                zeros.copy(), zeros.copy(), zeros.copy(),
                                columns.copy(), columns.copy())
    H = nested.hamiltonian(pair, nested.encode_state(pair,state))
    if maximum(H.imag) > 1e-10 or maximum(H-H.conj().T) > 1e-10:
        raise ValueError("actual standing Hamiltonian is not real Hermitian")
    levels, vectors = eigh((H.real+H.real.T)/2)
    positive = np.flatnonzero(levels > 1e-10*max(1.,maximum(H)))
    if len(positive) < 7:
        raise ValueError("fewer than seven resolved positive eigenmodes")
    ids = positive[:6]
    phi = vectors[:,ids].astype(complex)
    state.phi0,state.phi1 = phi[:grid.nf],phi[grid.nf:]
    covariance = (phi*WEIGHTS)@phi.conj().T
    spectral = {"eigenvalues":levels[ids], "gram_max":maximum(phi.conj().T@phi-np.eye(6)),
                "commutator_max":maximum(H@covariance-covariance@H),
                "eigen_residual_max":maximum(H@phi-phi*levels[ids]),
                "unequal_weight_gap":float(np.min(np.diff(levels[positive[:7]])[[1,3,5]])),
                "CAR_max":float(np.max(WEIGHTS)), "source_trace":float(np.sum(WEIGHTS))}
    return state,spectral


def radius_residual_jacobian(grid, state, rho):
    """Owned sampled C and its exact general-Q radial derivative/pullback."""
    fine = galerkin.prolong_state(grid,state)
    system = galerkin.active_fine_system(grid,fine)
    D,U = system.derivative,grid.A_g
    full = coupling.hamilton_constraint(system,fine)+rho
    Fr = -8*np.pi*system.A*fine.r
    Z = float(coupling.feedback_Z(system.A))
    image = (2*Z*(D@fine.r/fine.Q)[:,None]*(D@U)
             -16*np.pi*system.A*(fine.Q*fine.r)[:,None]*U
             -2*D@((D@(Fr[:,None]*U))/fine.Q[:,None]))
    return galerkin.pull_geometry(grid,full),full,galerkin.pull_geometry(grid,image)


def _positive_newton(values, oracle, lift, deadline, *, maximum_iterations=32):
    """Bounded preparation Newton with positivity and radius/seed correction stop."""
    values = np.array(values,copy=True)
    history=[]
    for iteration in range(maximum_iterations):
        if time.process_time() > deadline:
            return values,{"converged":False,"blocker":"preparation CPU cap", "history":history}
        residual,J = oracle(values)
        try:
            correction=solve(J,-residual)
        except np.linalg.LinAlgError:
            return values,{"converged":False,"blocker":"singular preparation Jacobian", "history":history}
        floor=max(1e-10,32*float(np.max(np.spacing(np.maximum(1.,np.abs(values))))))
        history.append({"iteration":iteration,"residual_max":maximum(residual),
                        "correction_max":maximum(correction),"correction_floor":floor})
        if maximum(correction) <= floor:
            return values,{"converged":True,"blocker":None,"history":history}
        accepted=False
        for halving in range(24):
            trial=values+2.**(-halving)*correction
            if not np.isfinite(trial).all() or np.min(lift(trial)) <= 0:
                continue
            other,_=oracle(trial)
            if np.linalg.norm(other) < np.linalg.norm(residual):
                values=trial;accepted=True;break
        if not accepted:
            return values,{"converged":False,"blocker":"positive preparation Newton stalled", "history":history}
    return values,{"converged":False,"blocker":"finite preparation iteration cap", "history":history}


def prepare_state(pair,Q,*,cpu_limit=15.):
    """Convex positive y seed, then exact finite C correction. No time stepping."""
    if not 0 < cpu_limit <= 30:
        raise ValueError("preparation CPU limit must be in (0,30]")
    start=time.process_time();deadline=start+cpu_limit;grid=pair.grid
    state,spectral=spectral_state(pair,Q)
    fine=galerkin.prolong_state(grid,state);system=galerkin.active_fine_system(grid,fine)
    source=coupling.source_from_columns(system,fine)
    rho=source["force_L"]/grid.dx_q
    mag=2*np.pi*system.C_F*system.flux**2
    B=(rho+mag*fine.Q)/(8*np.pi*system.A)
    if np.min(B) <= 0:
        raise ValueError("positive convex lapse forcing fails on the actual carrier")
    D,U=system.derivative,grid.A_g;DU=D@U
    def seed_oracle(y):
        fy=U@y
        residual=-4*D@((D@fy)/fine.Q)+fine.Q*fy-B/fy**3
        J=-4*D@(DU/fine.Q[:,None])+(fine.Q+3*B/fy**4)[:,None]*U
        return galerkin.pull_geometry(grid,residual),galerkin.pull_geometry(grid,J)
    y=np.full(grid.ng,float(np.mean(B/fine.Q))**.25)
    y,seed_report=_positive_newton(y,seed_oracle,lambda v:U@v,deadline)
    state.r=y*y
    if seed_report["converged"]:
        def exact_oracle(r):
            probe=state.copy();probe.r=r
            res,_,J=radius_residual_jacobian(grid,probe,rho)
            return res,J
        state.r,exact_report=_positive_newton(state.r,exact_oracle,lambda v:U@v,deadline)
    else:
        exact_report={"converged":False,"blocker":"convex seed unresolved", "history":[]}
    rate,bundle=galerkin.compose_fine_hamiltonian(grid,state)
    fine=bundle["fine_state"];source=bundle["source"]
    constraints=coupling.constraint_residuals(bundle["fine_system"],fine,source)
    retained=galerkin.pull_geometry(grid,constraints["hamilton"])
    tail=np.vstack((source["image0"],source["image1"]))-np.vstack((fine.phi0,fine.phi1))*spectral["eigenvalues"]
    positive=system.multiplicity*np.sum((np.abs(fine.phi0)**2+np.abs(fine.phi1)**2)*(WEIGHTS*spectral["eigenvalues"]),axis=1)/(grid.dx_q*fine.Q)
    report={"converged":bool(seed_report["converged"] and exact_report["converged"]),
            "blocker":exact_report["blocker"],"convex_seed":seed_report,"exact_correction":exact_report,
            "cpu_seconds":time.process_time()-start,"rho_min":float(np.min(constraints["rho"])),
            "convex_forcing_min":float(np.min(B)),"Q_min":float(np.min(fine.Q)),
            "r_min":float(np.min(fine.r)),"r_mean":float(np.mean(fine.r)),
            "retained_lapse_max":maximum(retained),"full_lapse_max":constraints["hamilton_max"],
            "unresolved_lapse_max":maximum(constraints["hamilton"]-U@retained),
            "full_shift_max":constraints["momentum_max"],
            "source_current_mean":float(np.mean(constraints["current"])),
            "fine_eigen_tail_max":maximum(tail),"pointwise_source_tail_max":maximum(constraints["rho"]-positive),
            "source_current_deleted":False,"spectral":{k:v for k,v in spectral.items() if k!="eigenvalues"},
            "actual_dynamic_force_max":{name:maximum(bundle["unprojected_rates"][i]) for name,i in (("p_Q",3),("p_r",4),("p_chi",5))},
            "continuum_initial_data_certified":False,"holding_claim":False,"evolution_performed":False}
    return state,report,bundle,constraints,spectral


def build_record(*,execute=False,nf=128,harmonic=8,lift=1.1,cpu_limit=30.):
    """Pure constructor/record builder. Writes require a separate explicit call."""
    if not 0 < cpu_limit <= 30:
        raise ValueError("whole preparation CPU limit must be in (0,30]")
    balance_metadata=json.loads(BALANCE_INPUT.read_text())
    balance_payload=BALANCE_INPUT.with_name(balance_metadata["payload"])
    report={"schema":SCHEMA,"status":"PREPARATION_PREFLIGHT","nf":int(nf),"harmonic":int(harmonic),
            "lift":float(lift),"cpu_limit_seconds":float(cpu_limit),"evidence_written":False,
            "cases":{},"source_hashes":{p:sha256(ROOT/p) for p in OWNERS},
            "input_hashes":{str(BALANCE_INPUT.relative_to(ROOT)):sha256(BALANCE_INPUT),
                            str(balance_payload.relative_to(ROOT)):sha256(balance_payload)},
            "scope":{"initial_constraint_preparation_only":True,"static_criticality_required":False,
                     "evolution_performed":False,"continuum_initial_data_certified":False}}
    if not execute:
        return report,{}
    start=time.process_time();arrays={}
    with threadpool_limits(limits=1):
        pair=make_pair(nf);grid=pair.grid
        report.update(ng=grid.ng,nq=grid.nq,period=grid.length,gauge="conformal",
                      occupations=WEIGHTS.tolist(),source_eta=1.,source_trace=3.,
                      multiplicity=grid.fine.multiplicity,locked_coefficients=grid.fine.coefficients)
        for name,fraction in (("uniform",0.),("pattern",1.)):
            remaining=cpu_limit-(time.process_time()-start)
            if remaining<=0:
                report["cases"][name]={"converged":False,"blocker":"whole preparation CPU cap",
                                       "checkpoint_available":False};break
            Q,declaration=declared_q(grid,harmonic=harmonic,lift=lift,pattern_fraction=fraction)
            state,measured,bundle,constraints,spectral=prepare_state(pair,Q,cpu_limit=remaining)
            report["cases"][name]=dict(measured,declaration=declaration,checkpoint_available=True)
            for key in nested.STATE_NAMES:
                arrays[name+"_nodal_"+key]=np.array(getattr(state,key),copy=True)
                arrays[name+"_"+key]=np.array(getattr(nested.encode_state(pair,state),key),copy=True)
            for key in ("rho","current","hamilton","momentum"):
                arrays[name+"_fine_"+key]=constraints[key].copy()
            arrays[name+"_eigenvalues"]=spectral["eigenvalues"].copy()
        arrays["geometry_map"]=pair.geometry_map.copy()
    report["cpu_seconds"]=time.process_time()-start
    report["status"]="FINITE_INITIAL_DATA" if len(report["cases"])==2 and all(v["converged"] for v in report["cases"].values()) else "UNRESOLVED_INITIAL_PREPARATION"
    return report,arrays


def output_paths(path):
    stem=Path(path).expanduser().resolve()
    if stem.suffix in (".json",".npz"):
        stem=stem.with_suffix("")
    return stem.with_suffix(".json"),stem.with_suffix(".npz")


def write_record(report,arrays,path=DEFAULT_OUTPUT):
    """Creation-only checkpoint, including actual nodal and nested states."""
    jp,payload_path=output_paths(path)
    if jp.exists() or payload_path.exists():
        raise FileExistsError("checkpoint already exists; use a successor")
    if not arrays:
        raise ValueError("preflight has no state checkpoint")
    stream=io.BytesIO();np.savez_compressed(stream,**arrays)
    payload=stream.getvalue();bound=dict(report,evidence_written=True,payload=payload_path.name,
                                       payload_sha256=hashlib.sha256(payload).hexdigest())
    jp.parent.mkdir(parents=True,exist_ok=True)
    with payload_path.open("xb") as f:
        f.write(payload)
    with jp.open("x") as f:
        json.dump(bound,f,indent=2,allow_nan=False);f.write("\n")
    return bound


def load_case(path,case="pattern"):
    jp,payload=output_paths(path);record=json.loads(jp.read_text())
    if sha256(payload)!=record["payload_sha256"]:
        raise ValueError("dynamic preparation payload hash mismatch")
    with np.load(payload,allow_pickle=False) as arrays:
        if case+"_nodal_phi0" not in arrays:
            raise ValueError("requested case has no prepared checkpoint")
        phi0=arrays[case+"_nodal_phi0"].copy();phi1=arrays[case+"_nodal_phi1"].copy()
        pair=nested.build_pair(record["nf"],coarse_modes=1,child_details=2,source_layout="override",
                               columns_override=(phi0,phi1),occupations=WEIGHTS)
        state=coupling.CauchyState(**{key:arrays[case+"_nodal_"+key].copy() for key in nested.STATE_NAMES})
        if maximum(pair.geometry_map-arrays["geometry_map"])>1e-12:
            raise ValueError("checkpoint geometry frame mismatch")
        encoded=nested.encode_state(pair,state)
        for key in nested.STATE_NAMES:
            if maximum(getattr(encoded,key)-arrays[case+"_"+key])>1e-10:
                raise ValueError("checkpoint nested state mismatch")
    return pair,encoded,record


def check_record(path):
    jp,payload=output_paths(path);record=json.loads(jp.read_text())
    if record.get("schema")!=SCHEMA:
        raise ValueError("dynamic preparation schema mismatch")
    for p,digest in {**record["source_hashes"],**record["input_hashes"]}.items():
        if sha256(ROOT/p)!=digest:
            raise ValueError("checkpoint source/input hash mismatch: "+p)
    with np.load(payload,allow_pickle=False) as arrays:
        for case,measured in record["cases"].items():
            if not measured.get("checkpoint_available",True):
                if measured.get("converged"):
                    raise ValueError("converged case lacks its checkpoint")
                continue
            pair,encoded,_=load_case(path,case)
            state=nested.reconstruct_state(pair,encoded)
            _,bundle=galerkin.compose_fine_hamiltonian(pair.grid,state)
            con=coupling.constraint_residuals(bundle["fine_system"],bundle["fine_state"],bundle["source"])
            for key in ("rho","current","hamilton","momentum"):
                if maximum(con[key]-arrays[case+"_fine_"+key])>1e-7:
                    raise ValueError("checkpoint owned constraint/source mismatch")
            if abs(con["hamilton_max"]-measured["full_lapse_max"])>1e-7:
                raise ValueError("checkpoint constraint summary mismatch")
    return {"ok":True,"bytes_written":0,"solver_executed":False,"cases":list(record["cases"])}
