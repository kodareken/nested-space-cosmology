"""One fixed-state physical-resolution comparison; no evolution or new field solve."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from threadpoolctl import threadpool_info, threadpool_limits

from recursive_horizons.nsc_spherical_coupling import (
    CauchyState, field_energy, gravity_energy, hamilton_constraint,
    locked_coefficients, shift_constraint,
)
from recursive_horizons.nsc_spherical_galerkin_coupling import (
    build_grid, compose_fine_hamiltonian,
)
from recursive_horizons.nsc_conformal_resolution import (
    log_rates, log_state, project_geometry, project_rates, reference_projector,
)
from recursive_horizons.nsc_conformal_vacuum_generator import lambda_generator
from measure_nsc_conformal_geometric_average import (
    EPISODE, CURVATURE, EXPECTED_PAYLOAD, PRIMARY, sha256,
)

ROOT = Path(__file__).resolve().parents[1]
REFERENCE_LENGTH = 4.80028044552352
CPU_LIMIT = 200.
BYTE_LIMIT = 24_000
GEOMETRY_FIELDS = ("Q", "r", "chi", "p_Q", "p_r", "p_chi")
SOURCE_PATHS = (
    "scripts/measure_nsc_conformal_resolution.py",
    "scripts/measure_nsc_conformal_geometric_average.py",
    "src/recursive_horizons/nsc_conformal_resolution.py",
    "tests/test_nsc_conformal_resolution.py",
    "tests/test_nsc_conformal_resolution_correlation.py",
    "src/recursive_horizons/nsc_conformal_vacuum_generator.py",
    "tests/test_nsc_conformal_vacuum_generator.py",
    "src/recursive_horizons/nsc_spherical_coupling.py",
    "src/recursive_horizons/nsc_spherical_galerkin_coupling.py",
    "src/recursive_horizons/nsc_spherical_feedback_action.py",
    "src/recursive_horizons/nsc_conformal_adm_source.py",
    "src/recursive_horizons/nsc_covariant_operator.py",
)


def frozen_sources():
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    hashes = {}
    for name in SOURCE_PATHS:
        raw = (ROOT/name).read_bytes()
        committed = subprocess.check_output(["git", "show", f"{commit}:{name}"], cwd=ROOT)
        if raw != committed:
            raise ValueError("resolution method must be frozen before application: "+name)
        hashes[name] = hashlib.sha256(raw).hexdigest()
    return commit, hashes


def norm_report(actual, candidate, dx):
    actual, candidate = np.asarray(actual), np.asarray(candidate)
    residual = actual-candidate
    rms = lambda value: float(np.sqrt(np.mean(abs(value)**2)))
    a, c, e = rms(actual), rms(candidate), rms(residual)
    return {"correlation_rms": a, "candidate_rms": c, "residual_rms": e,
            "correlation_max_abs": float(np.max(abs(actual))),
            "residual_max_abs": float(np.max(abs(residual))),
            "relative_residual_within_block": e/max(a,c) if max(a,c)>0 else 0.,
            "residual_squared_canonical_integral": float(dx*np.sum(abs(residual)**2))}


def saved_rate_replay_errors(grid, rate, saved_rates):
    """The archive writer stores prolonged projected rates, not ng coefficients."""
    errors = {}
    for name in GEOMETRY_FIELDS:
        saved = np.asarray(saved_rates[name])
        if saved.shape != (grid.nq,) or not np.isfinite(saved).all():
            raise ValueError("saved rate must use the archive's quadrature grid: "+name)
        errors[name] = float(np.max(abs(grid.A_g@getattr(rate,name)-saved)))
    return errors


def evaluate(grid, state):
    rate, bundle = compose_fine_hamiltonian(grid, state)
    system, fine, source = (bundle[key] for key in ("fine_system", "fine_state", "source"))
    # Reuse the two already computed Dirac images: no third field evaluation.
    energy = gravity_energy(system, fine)+field_energy(
        system, fine, images=(source["image0"], source["image1"]))
    ham = hamilton_constraint(system, fine)+source["force_L"]/grid.dx_q
    mom = shift_constraint(system, fine)+source["force_beta"]/grid.dx_q
    constraints = {
        "hamilton_max_abs": float(np.max(abs(ham))),
        "hamilton_rms": float(np.sqrt(np.mean(ham**2))),
        "momentum_max_abs": float(np.max(abs(mom))),
        "momentum_rms": float(np.sqrt(np.mean(mom**2))),
        "r_min_quadrature": float(np.min(fine.r)),
        "Q_min_quadrature": float(np.min(fine.Q)),
    }
    return energy, rate, constraints


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path,
                        default=ROOT/".build/data/nsc-spherical-conformal-episode-v2.npz")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("use a fresh summary path")
    if not args.output.parent.is_dir():
        raise ValueError("output parent must exist")
    started = time.process_time()

    def budget():
        if time.process_time()-started > CPU_LIMIT:
            raise RuntimeError("resolution probe exceeded 200 process-CPU seconds; no automatic retry")

    commit, bindings = frozen_sources()
    old = json.loads(CURVATURE.read_text())["source_bindings"]
    episode = json.loads(EPISODE.read_text())
    if (sha256(args.geometry) != EXPECTED_PAYLOAD or sha256(EPISODE) != old["episode_v2"]
            or episode["payload_sha256"] != EXPECTED_PAYLOAD or episode["gauge"] != "conformal"
            or episode["all_target_reached"] is not True or episode["source_changed"] is not False
            or episode["midtrajectory_reset"] is not False):
        raise ValueError("original conformal episode is not authenticated")
    for key, path in (("galerkin", "src/recursive_horizons/nsc_spherical_galerkin_coupling.py"),
                      ("action", "src/recursive_horizons/nsc_spherical_feedback_action.py")):
        if bindings[path] != old[key]:
            raise ValueError("original action or projection changed")
    # One useful 8-thread linear-algebra team, leaving other live processes alone.
    with threadpool_limits(limits=8):
        grid = build_grid(512, quadrature=2048, gauge="conformal")
        budget()
        with np.load(args.geometry, allow_pickle=False) as data:
            times = data[PRIMARY+"_frame_times"]
            i0, i1 = (int(np.argmin(abs(times-t))) for t in (.05,.3))
            if max(abs(float(times[i0])-.05),abs(float(times[i1])-.3)) > 1e-12:
                raise ValueError("required reference/final frame is absent")
            r0 = data[PRIMARY+"_frames_r"][i0]
            Q0 = data[PRIMARY+"_frames_Q"][i0]
            ruler = float((grid.A_g@r0)[grid.nq//8])
            if abs(ruler-REFERENCE_LENGTH)>1e-12:
                raise ValueError("independently fixed observer ruler changed")
            state = CauchyState(**{name: np.array(data[PRIMARY+"_frames_"+name][i1], copy=True)
                                  for name in (*GEOMETRY_FIELDS,"phi0","phi1")})
            saved_rates = {name: np.array(data[PRIMARY+"_frames_rate_"+name][i1], copy=True)
                           for name in GEOMETRY_FIELDS}
            if any(value.shape != (grid.nq,) or not np.isfinite(value).all()
                   for value in saved_rates.values()):
                raise ValueError("saved projected-rate representation differs from the archive writer")
        split = reference_projector(grid, r0, Q0, REFERENCE_LENGTH)
        budget()
        coarse = project_geometry(state, split.matrix)
        H_full, full_rate, full_constraints = evaluate(grid, state)
        budget()
        H_coarse, coarse_rate, coarse_constraints = evaluate(grid, coarse)
        budget()
        replay = saved_rate_replay_errors(grid,full_rate,saved_rates)
        for name, error in replay.items():
            if error > 1e-8*max(1.,float(np.max(abs(saved_rates[name])))):
                raise ValueError("same-action rate differs from its saved counterpart: "+name)
        full_log, coarse_log = log_rates(state, full_rate), log_rates(coarse, coarse_rate)
        delta_rate = project_rates({name: full_log[name]-coarse_log[name] for name in full_log},
                                   split.matrix)
        unit = lambda_generator(grid, coarse)
        coefficient = (H_full-H_coarse)/unit["energy"]
        candidate = project_rates(unit["canonical_rates"],split.matrix)
        blocks = {name: norm_report(delta_rate[name],coefficient*candidate[name],grid.dx_g)
                  for name in delta_rate}
        for name in ("phi0","phi1"):
            difference = getattr(full_rate,name)-getattr(coarse_rate,name)
            blocks[name] = norm_report(difference,np.zeros_like(difference),1.)
            blocks[name]["pairing_note"] = "Euclidean column norm; not a geometric canonical norm or energy"
            blocks[name].pop("residual_squared_canonical_integral")
        original = log_state(state)
        unresolved = {name: float(np.sqrt(np.mean((value-split.matrix@value)**2)))
                      for name,value in original.items()}
        budget()
        pools = [{key: item.get(key) for key in ("internal_api","num_threads","prefix")}
                 for item in threadpool_info()]

    result = {
        "schema": "nsc-conformal-resolution-summary-v1", "status": "measured_instantaneous_correlation",
        "source_commit": commit, "source_sha256": bindings,
        "input": {"episode_sha256":sha256(EPISODE), "payload_sha256":EXPECTED_PAYLOAD,
                  "payload_bytes":args.geometry.stat().st_size, "case":PRIMARY},
        "domain": {"reference_T":.05,"comparison_T":.3,"period":grid.length,
                   "topology":"S1 x S2","ng":grid.ng,"nf":grid.nf,"nq":grid.nq,
                   "gauge":"conformal L=Q, beta=0", "rate_clock":"same coordinate T generator",
                   "reference_length":REFERENCE_LENGTH,"cut_mu_times_Lref2_max":1.,
                   "projector_fixed_at_reference":True,"target_observation_used":False},
        "parameters":locked_coefficients(), "occupations":grid.fine.occupations.tolist(),
        "angular_kappa":grid.fine.kappa,"multiplicity_once":grid.fine.multiplicity,
        "projector": {**split.diagnostics,"retained_modes":int(split.retained),
                      "unresolved_modes":grid.ng-int(split.retained)},
        "unresolved_canonical_rms":unresolved,
        "full": {"hamiltonian":H_full,"constraints":full_constraints},
        "coarse": {"hamiltonian":H_coarse,"constraints":coarse_constraints},
        "saved_rate_replay_max_abs":replay,
        "energy_projection": {"delta_H":H_full-H_coarse,"unit_Lambda_H_at_coarse":unit["energy"],
                              "Lambda_in_action_length_units":coefficient,
                              "Lambda_times_Lref2":coefficient*REFERENCE_LENGTH**2,
                              "formula":"H_unitLambda=8*pi*A*integral(r^4*L*Q dx), L=Q",
                              "interpretation":"energy projection only; whole-force residual is required",
                              "rho_equals":"2*A*Lambda", "A_equals":"1/(16*pi*G)"},
        "canonical_rate_comparison":blocks,
        "full_fermion_columns_retained":bool(np.array_equal(state.phi0,coarse.phi0) and
                                             np.array_equal(state.phi1,coarse.phi1)),
        "limits": {"cpu_seconds":CPU_LIMIT,"eigensolves":1,"Hamiltonian_evaluations":2,
                   "field_evolution_steps":0,"new_field_solves":0,"parameter_scans":0},
        "runtime_cpu_seconds":time.process_time()-started,"blas_pools_during_application":pools,
        "scope":"Separate Einstein/Weyl/Dirac state; instantaneous frozen environment, not retarded closure. "
                "Original finite constraint defects remain. Bare vacuum and quantum B branch absent. "
                "No matched full/coarse proper clock, P001 parent, measured fractal dimension "
                "or cosmological Lambda identification.",
        "uncertainty":"Eigen/projector and same-state replay controls; no propagated continuum, model "
                      "or observational uncertainty bound for the correlation coefficient.",
    }
    raw = json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n"
    if len(raw.encode())>BYTE_LIMIT:
        raise ValueError("resolution summary exceeds 24KB")
    args.output.write_text(raw,encoding="utf-8",newline="\n")
    print(f"Wrote {args.output}; retained {split.retained}/{grid.ng}; "
          f"energy-projected Lambda*Lref^2={coefficient*REFERENCE_LENGTH**2:.12g}; "
          f"{result['runtime_cpu_seconds']:.3f} process-CPU seconds; no evolution.")


if __name__ == "__main__":
    main()
