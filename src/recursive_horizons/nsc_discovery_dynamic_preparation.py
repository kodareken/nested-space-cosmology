"""Positive general-Q, chi=0 Cauchy preparation using the owned finite action.

The source is selected once per initial geometry. No evolution, static force
condition, current deletion, or evidence write occurs in the constructors.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import io
import json
import subprocess
import sys
from contextlib import contextmanager
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


def prepare_state(pair,Q,*,cpu_limit=15.,rotation_angle=0.,columns_override=None):
    """Convex positive y seed, then exact finite C correction. No time stepping."""
    if not 0 < cpu_limit <= 30:
        raise ValueError("preparation CPU limit must be in (0,30]")
    start=time.process_time();deadline=start+cpu_limit;grid=pair.grid
    state,spectral=spectral_state(pair,Q)
    if columns_override is not None:
        state.phi0,state.phi1=(np.array(v,dtype=complex,copy=True) for v in columns_override)
        frame=np.vstack((state.phi0,state.phi1))
        if frame.shape!=(2*grid.nf,6) or maximum(frame.conj().T@frame-np.eye(6))>1e-9 or maximum(frame.imag)>1e-10:
            raise ValueError("coherence reference must be a real orthonormal frozen frame")
    if not np.isfinite(rotation_angle):
        raise ValueError("rotation angle must be finite")
    phi=np.vstack((state.phi0,state.phi1))
    spectral["rotation_reference_frame_sha256"]=hashlib.sha256(np.ascontiguousarray(phi).tobytes()).hexdigest()
    covariance_before=(phi*WEIGHTS)@phi.conj().T
    if rotation_angle:
        left,right=phi[:,0].copy(),phi[:,4].copy()
        phi[:,0]=np.cos(rotation_angle)*left+np.sin(rotation_angle)*right
        phi[:,4]=-np.sin(rotation_angle)*left+np.cos(rotation_angle)*right
        state.phi0,state.phi1=phi[:grid.nf],phi[grid.nf:]
    covariance=(phi*WEIGHTS)@phi.conj().T
    H=nested.hamiltonian(pair,nested.encode_state(pair,state))
    spectral["commutator_max"]=maximum(H@covariance-covariance@H)
    spectral["source_covariance_change_frobenius"]=float(np.linalg.norm(covariance-covariance_before))
    spectral["frozen_reference_order_and_sign_preserved"]=True
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
            "fine_eigen_tail_max":maximum(tail) if rotation_angle==0. else None,
            "physical_column_stationarity_defect_max":maximum(tail),
            "pointwise_source_tail_max":maximum(constraints["rho"]-positive) if rotation_angle==0. else None,
            "source_current_deleted":False,"spectral":{k:v for k,v in spectral.items() if k!="eigenvalues"},
            "rotation_angle":float(rotation_angle),"stationary_spectral_control":rotation_angle==0.,
            "rotation_is_charge_conjugation":False,
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


def check_record(path,*,observed_binding=None):
    jp,payload=output_paths(path);record=json.loads(jp.read_text())
    if record.get("schema")!=SCHEMA:
        raise ValueError("dynamic preparation schema mismatch")
    authenticate_preparation(path,observed_binding=observed_binding)
    for p in OWNERS[4:]:
        if p in record["source_hashes"] and sha256(ROOT/p)!=record["source_hashes"][p]:
            raise ValueError("historical physical owner changed; live numerical replay unavailable")
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


def _git_hashes(commit,hashes):
    """Read immutable Git blobs; never update a historical producer pin."""
    resolved=subprocess.check_output(["git","rev-parse","--verify",commit+"^{commit}"],cwd=ROOT,text=True).strip()
    for path,digest in hashes.items():
        blob=subprocess.check_output(["git","show",resolved+":"+path],cwd=ROOT)
        if hashlib.sha256(blob).hexdigest()!=digest:
            raise ValueError("producer Git blob does not match: "+path)
    return resolved


def authenticate_preparation(path,*,observed_binding=None):
    jp,payload=output_paths(path);record=json.loads(jp.read_text())
    if sha256(payload)!=record["payload_sha256"]:
        raise ValueError("frozen preparation payload mismatch")
    for p,digest in record["input_hashes"].items():
        if sha256(ROOT/p)!=digest:
            raise ValueError("frozen preparation input mismatch: "+p)
    commit=record.get("producing_commit")
    binding=None
    if observed_binding is not None:
        binding=json.loads(Path(observed_binding).read_text()) if isinstance(observed_binding,(str,Path)) else dict(observed_binding)
        if binding.get("record_sha256")!=sha256(jp) or binding.get("payload_sha256")!=sha256(payload):
            raise ValueError("external observed binding does not name this exact preparation")
        commit=binding["producing_commit"]
    if commit:
        _git_hashes(commit,record["source_hashes"])
    else:
        for p,digest in record["source_hashes"].items():
            if sha256(ROOT/p)!=digest:
                raise ValueError("historical preparation needs an external observed producer binding")
    return {"record_sha256":sha256(jp),"payload_sha256":sha256(payload),
            "producing_commit":commit,"external_post_run_observation":binding is not None,
            "old_record_modified":False}


def dependency_hashes():
    """Observed local import closure, including the frozen diagnostic consumer."""
    from . import nsc_discovery_episode, nsc_discovery_extent, nsc_discovery_tidal
    paths=set(OWNERS)
    for module in list(sys.modules.values()):
        file=getattr(module,"__file__",None)
        if file is None:
            continue
        path=Path(file).resolve()
        if path.is_relative_to(LAB/"src") and path.suffix==".py":
            paths.add(str(path.relative_to(ROOT)))
    return {p:sha256(ROOT/p) for p in sorted(paths)}


def build_coherent_successor(initial,*,observed_binding=None,execute=False,cpu_limit=30.,producer_commit=None):
    """Preserve both frozen controls; solve only the declared rotated source."""
    binding=authenticate_preparation(initial,observed_binding=observed_binding)
    jp,payload=output_paths(initial);old=json.loads(jp.read_text())
    report=dict(old,status="COHERENT_PREPARATION_PREFLIGHT",evidence_written=False,
                predecessor_binding=binding,source_hashes=dependency_hashes(),
                input_hashes={str(jp):sha256(jp),str(payload):sha256(payload)},cases={})
    report.pop("payload",None);report.pop("payload_sha256",None)
    if producer_commit:
        report["producing_commit"]=_git_hashes(producer_commit,report["source_hashes"])
    elif execute:
        raise ValueError("successor execution requires a frozen producing commit")
    if not execute:
        return report,{}
    if not 0<cpu_limit<=30:
        raise ValueError("coherent preparation budget must be in (0,30]")
    start=time.process_time()
    with np.load(payload,allow_pickle=False) as original:
        arrays={key:original[key].copy() for key in original.files}
    report["cases"]={case:dict(old["cases"][case],frozen_control_copied=True) for case in ("uniform","pattern")}
    pair,encoded,_=load_case(initial,"pattern")
    frozen=nested.reconstruct_state(pair,encoded)
    with threadpool_limits(limits=1):
        state,measured,bundle,constraints,spectral=prepare_state(pair,frozen.Q,cpu_limit=cpu_limit,
            rotation_angle=np.pi/6,columns_override=(frozen.phi0,frozen.phi1))
    report["cases"]["coherent"]=dict(measured,checkpoint_available=True,
        declaration=dict(old["cases"]["pattern"]["declaration"],coherence_rotation_columns=[0,4],coherence_angle=float(np.pi/6)))
    coherent_pair=replace(pair,source_phi0=state.phi0.copy(),source_phi1=state.phi1.copy(),
                          source_columns=np.vstack((state.phi0,state.phi1)))
    encoded=nested.encode_state(coherent_pair,state)
    for key in nested.STATE_NAMES:
        arrays["coherent_nodal_"+key]=np.array(getattr(state,key),copy=True)
        arrays["coherent_"+key]=np.array(getattr(encoded,key),copy=True)
    for key in ("rho","current","hamilton","momentum"):
        arrays["coherent_fine_"+key]=constraints[key].copy()
    arrays["coherent_eigenvalues"]=spectral["eigenvalues"].copy()
    base=arrays["pattern_fine_rho"];changed=constraints["rho"]
    report["cases"]["coherent"]["source_density_difference"]={
        "maximum":maximum(changed-base),"relative_contrast":float((np.max(changed)-np.min(changed))/np.mean(changed)),
        "source_covariance_changed":True,"charge_conjugation":False}
    report["cpu_seconds"]=time.process_time()-start
    report["predecessor_preparation_cpu_seconds"]=float(old.get("cpu_seconds",0.))
    report["status"]="FINITE_INITIAL_DATA" if measured["converged"] else "UNRESOLVED_INITIAL_PREPARATION"
    if report["source_hashes"]!=dependency_hashes():
        raise ValueError("producer closure changed during coherent preparation")
    return report,arrays


EPISODE_STATIONS=(.3,1.,3.)
EPISODE_SCHEMA="NSC-DISCOVERY-DYNAMIC-EPISODE-v1"
EXTENT_OBSERVER_SHA256="0b9adea69e33f4e1683d50f44d6961f76ab71cea4533d64c025d5c8aa7080170"
_observer_normal_clocks=None


def _dynamic_stations(stations=None):
    values=EPISODE_STATIONS if stations is None else tuple(float(v) for v in stations)
    if not values or list(values)!=sorted(set(values)) or any(v not in EPISODE_STATIONS for v in values):
        raise ValueError("dynamic diagnostic stations are .3,1,3")
    return values


def _pinned_scalar_row(pair,state,time_value,control_mode="coupled"):
    """Use frozen extent real-integral/tidal observer, never the dynamic atlas."""
    from . import nsc_discovery_extent as extent
    from . import nsc_discovery_tidal as tidal
    from . import nsc_regional_energy_exchange as regional
    started=time.process_time()
    if sha256(extent.__file__)!=EXTENT_OBSERVER_SHA256:
        raise ValueError("period-aware observer differs from its explicit freeze")
    if _observer_normal_clocks is None:
        raise ValueError("owned episode normal-clock context is absent")
    row,_=extent.observe(pair,state,time_value,_observer_normal_clocks,control_mode)
    jets=tidal.analytic_accelerations(pair,state,control_mode)
    fine,system,bundle=jets["fine"],jets["bundle"]["fine_system"],jets["bundle"]
    source=bundle["source"];rate=jets["rate"];grid=pair.grid
    lifted=coupling.CauchyRate(*[galerkin.prolong_geometry(grid,getattr(rate,k)) for k in
        ("Q","r","chi","p_Q","p_r","p_chi")],
        galerkin.prolong_columns(grid,rate.phi0),galerkin.prolong_columns(grid,rate.phi1),
        rate.force_L,rate.force_Q,rate.force_beta,rate.fieldwork_power)
    ledger=regional.matter_ledger(system,fine)
    terms=regional.proper_balance_terms(system,fine,lifted,ledger)
    integral=lambda value,interval:extent.real_interval_integral(grid,value/grid.dx_q,interval)
    endpoint=lambda interval:extent.real_periodic_values(grid,terms["flux_nodal"]/grid.dx_q,interval)
    intervals={"left":(0.,1.),"child":(1.,3.),"right":(3.,4.),"ambient":(4.,grid.length)}
    accounting={}
    for name,interval in intervals.items():
        flux=endpoint(interval)
        accounting[name]={"interval":list(interval),"normal_energy":integral(ledger["normal_energy_nodal"],interval),
            "pressure_work":integral(terms["proper_pressure_work"],interval),
            "lapse_work":integral(terms["momentum_lapse_work"],interval),
            "boundary_flux":float(flux[0]-flux[1]),"physical_wall":False}
    statistics=lambda value:{"min":float(np.min(value)),"max":float(np.max(value)),"max_abs":maximum(value)}
    mass=np.sum((np.abs(fine.phi0)**2+np.abs(fine.phi1)**2)*pair.weights,axis=1)
    proper_q=fine.r*fine.Q
    probability_density=mass/(grid.dx_q*proper_q)
    contrast=float((np.max(probability_density)-np.min(probability_density))/np.mean(probability_density))
    flat=contrast<=1e-8
    proper_length=extent.real_interval_integral(grid,proper_q,(0.,grid.length))
    inverse_width=extent.real_interval_integral(grid,probability_density**2*proper_q,(0.,grid.length))
    participation_width=float(np.sum(mass)**2/inverse_width)
    peak=int(np.argmax(probability_density));x=grid.xi_q
    distance=(x-x[peak]+grid.length/2)%grid.length-grid.length/2
    childflux=endpoint(pair.child_interval);parentflux=endpoint(pair.parent_interval)
    row.update(control_mode=control_mode,coordinate_fieldwork=float(rate.fieldwork_power),
        pressure_work=sum(v["pressure_work"] for v in accounting.values()),
        lapse_work=sum(v["lapse_work"] for v in accounting.values()),
        boundary_child=float(childflux[0]-childflux[1]),boundary_parent=float(parentflux[0]-parentflux[1]),
        localization_peak_x=None if flat else float(x[peak]),localization_coordinate_width=float(np.sqrt(np.sum(mass*distance**2)/np.sum(mass))),
        localization_width_domain="coordinate second moment; proper lengths separately owned",partition_accounting=accounting,
        clock_values_carried_in_checkpoint=True,observation_cpu_seconds=time.process_time()-started,
        quadrature="trapezoid_coordinate_time",dense_propagator_stored=False,
        energy_total=row["total_energy"]["total"],energy_field=row["total_energy"]["field"],
        energy_gravity=row["total_energy"]["gravity"],regime="finite observed state; holding unresolved")
    row["tagged_source_interpretation"]="fixed occupied column slots 2,3; no regional ancestry claim"
    row["all_column_localization"]={"probability_per_proper_length":statistics(probability_density),
        "contrast":contrast,"flat":flat,"flat_roundoff_indicator_threshold":1e-8,
        "peak_x":None if flat else float(x[peak]),"participation_proper_width":participation_width,
        "carrier_proper_length":proper_length,"participation_ratio":participation_width/proper_length,
        "primary_source_measurement":True,"continuum_error_bound":None}
    row["actual_rates"]={key:statistics(getattr(lifted,key)) for key in ("Q","r","chi","p_Q","p_r","p_chi")}
    row["actual_source_forces"]={key:statistics(source[key]) for key in ("force_Q","force_L","force_beta")}
    row["accounting_scope"]="disjoint instantaneous terms; integrated trapezoid work is not a propagated-error certificate"
    return row


def _pinned_observe(pair,state,time_value,*,control_mode="coupled"):
    row=_pinned_scalar_row(pair,state,time_value,control_mode)
    return {"time":float(time_value),"external_observation":row,"stability":row,
            "observer_module":"nsc_discovery_extent frozen consumer", "notes":[]}


def _real_clock_rates(pair,state):
    from . import nsc_discovery_extent as extent
    if sha256(extent.__file__)!=EXTENT_OBSERVER_SHA256:
        raise ValueError("period-aware clock sampler differs from its freeze")
    return np.asarray(extent.real_metrics(pair,state)["clock_rates"],dtype=float)


def _event_aware_advance(pair,state,**kwargs):
    """Owned advance_case in cadence segments, retaining a measured sign event.

    CPU admission includes the entire segment, including observations. There
    is no second time stepper and no state correction between segments.
    """
    from . import nsc_discovery_episode as episode
    global _observer_normal_clocks
    original=_advance_owner
    current=state;options=dict(kwargs);all_rows=[];snapshots=[];reached=[]
    used=float(options.get("spent",0.));step_cpu=0.;observation_cpu=0.
    initial=None;event_info=(options.get("channel_sample") or {}).get("first_source_density_sign_change")
    while True:
        mark=float(options["coordinate_time"])
        end=min(episode.next_observation_time(mark,options.get("observation_cadence",.05)),max(options["stations"]))
        previous=options.get("until_time")
        if previous is not None:
            end=min(end,float(previous))
        start=time.process_time()
        _observer_normal_clocks=np.array(options.get("normal_clocks",np.zeros(3)),copy=True)
        result=original(pair,current,**dict(options,until_time=end,spent=used))
        _observer_normal_clocks=result["normal_clocks"].copy()
        for row in result["observations"]:
            if abs(row["time"]-result["time"])<1e-12:
                row["normal_clocks"]=result["normal_clocks"].tolist()
                row["clock_values_source"]="owned advance_case returned endpoint clocks"
        used+=time.process_time()-start
        step_cpu+=result["step_cpu_seconds"];observation_cpu+=result["observation_cpu_seconds"]
        all_rows.extend(result["observations"])
        if initial is None and all_rows:
            initial=all_rows[0]
        reached=sorted(set(reached+result["stations_reached"]))
        for snapshot in result["snapshots"]:
            snapshot["stations_reached"]=[v for v in options["stations"] if v<=snapshot["time"]+1e-12]
            if snapshot["kind"]=="station":
                if not snapshots or snapshots[-1] is not snapshot:
                    snapshots.append(snapshot)
        last=all_rows[-1] if all_rows else {}
        changed_sign=(initial is not None and initial.get("source_rho",{}).get("min",0)>0
                      and last.get("source_rho",{}).get("min",1)<=0 and event_info is None)
        continuing=result["status"]=="until_time" and result["time"]<max(options["stations"])-1e-12
        if changed_sign:
            event_info={"event":"initially positive sampled lapse source changed sign","time":result["time"],
                        "physical_instability_claimed":False,"finite_sampled_domain":True,"automatic_stop":False}
            last["source_density_sign_change_first_observation"]=True
            if result["snapshots"]:
                snapshot=result["snapshots"][-1]
                snapshot["kind"]="event";snapshot["channel_sample"]=dict(last,first_source_density_sign_change=event_info)
                if not snapshots or snapshots[-1] is not snapshot:
                    snapshots.append(snapshot)
        if event_info is not None:
            last["first_source_density_sign_change"]=event_info
            if result["channel_sample"] is not None:
                result["channel_sample"]["first_source_density_sign_change"]=event_info
        if last.get("observation_error") or last.get("CAR_min",0)<-1e-8 or last.get("CAR_max",0)>1+1e-8:
            result["status"]="diagnostic_admissibility_stop"
            result["stop"]={"stop_class":"diagnostic_admissibility_stop","reason":last.get("observation_error") or "CAR numerical domain",
                            "physical_instability_claimed":False}
            continuing=False
        if used>=float(options["cpu_allowance"]):
            result["status"]="budget_stop";result["stop"]={"stop_class":"budget_stop","includes_observation_cpu":True}
            continuing=False
        if not continuing:
            if result["snapshots"] and (not snapshots or snapshots[-1]["time"]!=result["snapshots"][-1]["time"]):
                snapshots.append(result["snapshots"][-1])
            result.update(observations=all_rows,snapshots=snapshots,stations_reached=reached,
                          cpu_seconds=used,step_cpu_seconds=step_cpu,observation_cpu_seconds=observation_cpu)
            return result
        current=result["state"]
        options.update(coordinate_time=result["time"],steps=result["steps"],normal_clocks=result["normal_clocks"],
                       work_ledger=result["work_ledger"],channel_sample=result["channel_sample"],channel_time=result["channel_time"])


@contextmanager
def episode_adapter():
    """Scoped runner hooks; retain owned rates/RK4/CFL and clock quadrature."""
    from . import nsc_discovery_episode as episode
    global _advance_owner
    if episode.advance_case is _event_aware_advance:
        # A fork worker may inherit the parent's already-installed adapter.
        yield episode
        return
    saved=(episode.normalize_stations,episode.scalar_observation_row,episode.observe,episode.advance_case,episode._clock_rates)
    _advance_owner=episode.advance_case
    episode.normalize_stations=_dynamic_stations
    episode.scalar_observation_row=_pinned_scalar_row
    episode.observe=_pinned_observe
    episode.advance_case=_event_aware_advance
    episode._clock_rates=_real_clock_rates
    try:
        yield episode
    finally:
        episode.normalize_stations,episode.scalar_observation_row,episode.observe,episode.advance_case,episode._clock_rates=saved


def _worker_initializer():
    global _worker_adapter
    _worker_adapter=episode_adapter();_worker_adapter.__enter__()


def _episode_pool(max_workers):
    return _AdapterPool(max_workers)


def _dynamic_case_worker(payload):
    from . import nsc_discovery_episode as episode
    global _observer_normal_clocks
    _,arrays=episode.load_checkpoint(payload["directory"],payload["case_id"])
    _observer_normal_clocks=np.array(arrays["normal_clocks"],copy=True)
    return episode._case_worker(payload)


class _AdapterPool:
    """Delegate the existing pool worker after installing the thin hooks."""
    def __init__(self,max_workers):
        from concurrent.futures import ProcessPoolExecutor
        self.pool=ProcessPoolExecutor(max_workers=max_workers,initializer=_worker_initializer)

    def __enter__(self):
        self.pool.__enter__();return self

    def __exit__(self,*args):
        return self.pool.__exit__(*args)

    def submit(self,function,payload):
        return self.pool.submit(_dynamic_case_worker,payload)


def prepare_episode(initial,output,*,observed_binding=None,producer_commit=None,execute=False,
                    resolved_error=1e-6,stations=EPISODE_STATIONS):
    """Freeze six exact native-basis initial checkpoints. No stepping/initial solve."""
    started=time.process_time()
    if not np.isfinite(resolved_error) or resolved_error<=0:
        raise ValueError("declared initial resolved error must be positive finite")
    binding=authenticate_preparation(initial,observed_binding=observed_binding)
    hashes=dependency_hashes()
    report={"schema":EPISODE_SCHEMA,"status":"EPISODE_PREFLIGHT","stations":list(_dynamic_stations(stations)),
            "source_hashes":hashes,"initial_binding":binding,"workers_maximum":4,"cpu_budget_seconds":300.,
            "cases":[],"evolved":False,"initial_constraints_resolved_error":float(resolved_error),
            "initial_record":str(output_paths(initial)[0]),
            "physical_event_policy":"record first finite source sign change; continue without positivity clipping or a sector veto"}
    if not execute:
        return report
    if producer_commit is None:
        raise ValueError("episode preparation requires its frozen producing commit")
    report["producing_commit"]=_git_hashes(producer_commit,hashes)
    from . import nsc_discovery_episode as episode
    directory=episode.assert_campaign_output(output)
    if episode._manifest_path(directory).exists():
        raise FileExistsError("dynamic episode successor already exists")
    jp,_=output_paths(initial);saved=json.loads(jp.read_text())
    for case in ("uniform","pattern","coherent"):
        measured=saved["cases"].get(case)
        if measured is None or not measured["converged"] or measured["full_lapse_max"]>resolved_error or measured["full_shift_max"]>resolved_error:
            raise ValueError("owned initial constraints unresolved for "+case)
    loaded={}
    for case in ("uniform","pattern","coherent"):
        pair,state,_=load_case(initial,case)
        _,bundle=nested.rates(pair,state,return_bundle=True)
        actual=coupling.constraint_residuals(bundle["fine_system"],bundle["fine_state"],bundle["source"])
        if actual["hamilton_max"]>resolved_error or actual["momentum_max"]>resolved_error:
            raise ValueError("actual owned initial constraints exceed declared error for "+case)
        if coupling.chart_failure(bundle["fine_system"],bundle["fine_state"]) is not None:
            raise ValueError("initial state is outside the actual positive chart")
        loaded[case]=(pair,state)
    directory.mkdir(parents=True,exist_ok=True)
    for case in ("uniform","pattern","coherent"):
        pair,state=loaded[case]
        basis={"W":pair.geometry_map,"source_phi0":pair.source_phi0,"source_phi1":pair.source_phi1,
               "observer_columns":pair.reference_columns,"source_weights":pair.weights}
        arrays=episode.arrays_from_state(state,basis_arrays=basis,clocks={"rates":_real_clock_rates(pair,state),"normal_clocks":np.zeros(3)})
        pins={"W_sha256":episode.array_sha256(pair.geometry_map),"source_columns_sha256":episode.array_sha256(pair.source_columns),
              "observer_columns_sha256":episode.array_sha256(pair.reference_columns),"weights_sha256":episode.array_sha256(pair.weights),
              "initial_state_sha256":episode.state_sha256(state),"initial_record_sha256":binding["record_sha256"]}
        # Phi equals its preparation at T=0 by definition. Whole-state hash pins
        # that value; episode's later Phi/source inequality is not an initial test.
        for cap in (.001,.0005):
            identifier=f"nf{pair.grid.nf}_{case}_dt{cap}"
            record={"case_id":identifier,"nf":pair.grid.nf,"coordinate_time":0.,"steps":0,
                "stations":report["stations"],"stations_reached":[],"step_cap":cap,"control_mode":"coupled","geometry":"evolving",
                "momentum_representation":episode.CANONICAL_PI,"source_pins":pins,"source_pins_before":pins,"source_pins_after":pins,
                "coarse_indices":pair.geometry_coarse_indices.tolist(),"child_indices":pair.geometry_child_indices.tolist(),
                "parent_indices":pair.geometry_parent_indices.tolist(),"source_metadata":pair.source_metadata,"geometry_metadata":pair.geometry_metadata,
                "clock_locations":list(pair.clock_locations),"clock_protocol":"endpoint_trapezoid_per_owned_RK4_step",
                "initial_state_called":False,"verify_external_pins":False,"status":"prepared","snapshot_kind":"handoff",
                "work_ledger":episode.empty_work_ledger(),"source_hashes":hashes,"producing_commit":report["producing_commit"],
                "initial_binding":binding,"preparation_case":case,"initial_phi_is_prepared_source":True}
            committed=episode.commit_checkpoint(directory,record,arrays)
            report["cases"].append({key:committed[key] for key in ("case_id","nf","npz","json","ordinal","coordinate_time","steps")})
    report.update(status="PREPARED",stage=0,workers=4,chunk_limit_bytes=episode.CHUNK_LIMIT_BYTES,child_cpu_seconds=0.,
                  memory_limit_bytes=episode.DEFAULT_MEMORY_BYTES,forecast_factor=episode.DEFAULT_FORECAST_FACTOR,
                  initial_preparation_cpu_seconds=float(saved.get("cpu_seconds",0.))+float(saved.get("predecessor_preparation_cpu_seconds",0.)),
                  adapter_preparation_cpu_seconds=time.process_time()-started,
                  observer_sha256=EXTENT_OBSERVER_SHA256,clock_sampler="frozen extent.real_metrics; same rQ normal clock")
    charged=report["initial_preparation_cpu_seconds"]+report["adapter_preparation_cpu_seconds"]
    _,ledger=episode._ledger_update(directory,pilot=charged,budget=300.)
    report["cumulative_cpu_seconds"]=ledger["spent"]
    episode._write_json(episode._manifest_path(directory),report)
    return report


def run_episode(output,*,workers=4,cpu_budget=300.,accepted_six_hour=False,max_steps=None):
    """Reuse the existing pool, admission ledger, RK4, checkpoints and clocks."""
    if not 1<=workers<=4 or not 0<cpu_budget<=21600 or cpu_budget>300 and not accepted_six_hour:
        raise ValueError("at most four workers; over300CPU requires accepted six-hour budget")
    from . import nsc_discovery_episode as episode
    manifest=episode.read_manifest(output)
    if manifest.get("schema")!=EPISODE_SCHEMA:
        raise ValueError("not a dynamic preparation episode")
    if dependency_hashes()!=manifest["source_hashes"]:
        raise ValueError("episode dependency closure changed before evolution")
    _git_hashes(manifest["producing_commit"],manifest["source_hashes"])
    live=authenticate_preparation(manifest["initial_record"])
    if live["record_sha256"]!=manifest["initial_binding"]["record_sha256"] or live["payload_sha256"]!=manifest["initial_binding"]["payload_sha256"]:
        raise ValueError("initial preparation reference changed before evolution")
    with episode_adapter():
        result=episode.run(output,workers=workers,cpu_budget_seconds=cpu_budget,max_steps=max_steps,executor=_episode_pool)
    ledger=json.loads(episode._ledger_path(output).read_text())
    result.update(schema=EPISODE_SCHEMA,cumulative_cpu_seconds=ledger["spent"],
                  physical_outcome="finite checkpointed measurements or recorded blocker; startup force is not holding")
    result["assessment"]=assess_episode(output,result)
    episode._write_json(episode._manifest_path(output),result)
    return result


def assess_episode(output,manifest=None):
    """Measured effects and matched-step indicators from stored observations."""
    from . import nsc_discovery_episode as episode
    manifest=episode.read_manifest(output) if manifest is None else manifest
    channels={"child_proper_length":("child_proper_length",),
        "child_areal_radius":("child_r_proper_mean",),"child_Q":("child_Q_proper_mean",),
        "child_probability":("all_source_child_probability",),
        "child_normal_energy":("child_matter_normal_energy",),
        "proper_probability_width":("all_column_localization","participation_proper_width"),
        "probability_contrast":("all_column_localization","contrast"),
        "radial_tide_rms":("actual_tides","R_0101","rms"),
        "angular_tide_rms":("actual_tides","R_0202","rms")}
    def values(row):
        result={}
        for name,path in channels.items():
            value=row
            for key in path:
                value=value.get(key) if isinstance(value,dict) else None
            if value is not None and np.isfinite(value):
                result[name]=float(value)
        return result
    cases={};series={};origins={}
    for entry in manifest["cases"]:
        identifier=entry["case_id"]
        rows=episode.read_observations(output,identifier)
        record,_=episode.load_checkpoint(output,identifier)
        origins[identifier]=record.get("source_pins",{}).get("initial_state_sha256")
        final=record.get("stability")
        if isinstance(final,dict) and "all_column_localization" in final:
            rows.append(final)
        ordered=sorted({float(row["time"]):row for row in rows if "time" in row}.values(),key=lambda row:row["time"])
        series[identifier]=ordered
        if len(ordered)<2 or ordered[-1]["time"]<=ordered[0]["time"]:
            cases[identifier]={"blocker":record.get("stop") or "no evolved physical observation",
                               "holding_claim":False,"evolved_effect_measured":False}
            continue
        initial,last=ordered[0],ordered[-1];a,b=values(initial),values(last)
        changes={name:b[name]-a[name] for name in a.keys()&b.keys()}
        cases[identifier]={"initial_time":initial["time"],"last_time":last["time"],"initial":a,"last":b,"changes":changes,
            "geometry_signs":{name:int(np.sign(changes[name])) for name in ("child_proper_length","child_areal_radius") if name in changes},
            "source_flat_initial":initial["all_column_localization"]["flat"],
            "source_flat_last":last["all_column_localization"]["flat"],
            "first_source_density_sign_change":last.get("first_source_density_sign_change"),
            "last_constraint_pair":last.get("constraint_pair"),"stop":record.get("stop"),
            "evolved_effect_measured":True,"holding_claim":False,"continuum_error_bound":None}
    comparisons=[]
    for preparation in ("uniform","pattern","coherent"):
        identifiers=[entry["case_id"] for entry in manifest["cases"] if "_"+preparation+"_" in entry["case_id"]]
        if len(identifiers)!=2:
            continue
        left,right=(series[key] for key in identifiers)
        right_map={float(row["time"]):row for row in right}
        for row in left:
            target=float(row["time"])
            other=next((v for t,v in right_map.items() if abs(t-target)<1e-10),None)
            if target<=0 or other is None:
                continue
            a,b=values(row),values(other)
            comparisons.append({"preparation":preparation,"time":target,
                "case_ids":identifiers,"matched_initial_state":origins[identifiers[0]] is not None and origins[identifiers[0]]==origins[identifiers[1]],
                "second_minus_first":{name:b[name]-a[name] for name in a.keys()&b.keys()},
                "finite_timestep_indicator_only":True,"spatial_convergence_claim":False})
    return {"cases":cases,"matched_step_comparisons":comparisons,
            "source_frame_matched_between_caps":True,"startup_force_is_holding":False,
            "finite_observed_domain_only":True}
