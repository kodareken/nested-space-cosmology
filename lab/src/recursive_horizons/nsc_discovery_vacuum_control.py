"""Bounded source-free Bertotti–Robinson control of the locked conformal action.

The canonical matter covariance is zero. The induced and magnetic terms remain
the existing locked action. This prepares a different initial radius from the
occupied discovery family and is a conditional comparison, not its replay.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import io
import json
import math
import os
from pathlib import Path
import time

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import numpy as np
from threadpoolctl import threadpool_limits

from . import nsc_discovery_tidal as tidal
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_feedback_action as action
from . import nsc_spherical_galerkin_coupling as galerkin

SCHEMA = "NSC-DISCOVERY-VACUUM-CONTROL-v1"
LAB = Path(__file__).resolve().parents[2]
OUTPUT = LAB / "results/development/nsc-discovery-vacuum-control-v1"
FIELDS = ("Q", "r", "chi", "p_Q", "p_r", "p_chi", "phi0", "phi1")
EXACT_TIMES = (0.0, 1.0, 3.0, 8.0)
RUN_TIMES = (0.0, 1.0, 3.0)
FERMION_COUNTS = (32, 64)
STEP_DT = 0.01
CPU_BUDGET = 30.0
MAX_PAYLOAD_BYTES = 1 << 20


def sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    paths = [Path(__file__).resolve(),
             LAB / "scripts/derive_nsc_discovery_vacuum_control.py",
             LAB / "tests/test_nsc_discovery_vacuum_control.py",
             LAB / "docs/nsc-discovery-vacuum-control.md"]
    paths.extend(Path(path) for path in tidal.producer_hashes())
    return {str(path.resolve()): sha256_file(path) for path in paths}


def control_grid(nf):
    """Original conformal grid with six zero occupations, no coefficient change."""
    if int(nf) not in FERMION_COUNTS:
        raise ValueError("the bounded control uses nf=32 or nf=64")
    grid = galerkin.build_grid(int(nf), gauge="conformal")
    grid.fine = replace(grid.fine, occupations=np.zeros(6))
    return grid


def parameters(grid):
    system = grid.fine
    radius2 = coupling.magnetic_radius_square(system.coefficients)
    if radius2 <= 0:
        raise ValueError("the locked magnetic radius must be positive")
    return {
        "radius": math.sqrt(radius2),
        "radius_square": radius2,
        "Q0": float(system.calibration["b0"] / system.calibration["a0"]),
        "alpha": float(action.alpha_of(system.C_W)),
        "A": float(system.A),
        "F_chi": float(action.partial_F(math.sqrt(radius2), system.A, system.C_W)[1]),
    }


def exact_values(grid, coordinate_time):
    """Closed solution and its derivatives, derived from the action momenta."""
    p = parameters(grid)
    t = float(coordinate_time)
    Q = p["Q0"] / math.cosh(p["Q0"] * t)
    D = -p["Q0"] * math.tanh(p["Q0"] * t)
    Ddot = -Q * Q
    Fr = float(action.partial_F(p["radius"], p["A"], grid.fine.C_W)[0])
    return {
        **p, "coordinate_time": t, "Q": Q, "D": D,
        "Q_t": Q * D, "Q_tt": Q * (D * D + Ddot),
        "p_Q": 0.0, "p_chi": 2 * p["F_chi"] * D, "p_r": 2 * Fr * D,
        "p_chi_t": 2 * p["F_chi"] * Ddot, "p_r_t": 2 * Fr * Ddot,
        "normal_proper_time": p["radius"] * math.atan(math.sinh(p["Q0"] * t)),
        "period_proper_length": grid.length * p["radius"] * Q,
        "child_proper_length": 2 * p["radius"] * Q,
    }


def exact_state(grid, coordinate_time):
    """Native CauchyState uses nodal momentum densities, not canonical pi."""
    v = exact_values(grid, coordinate_time)
    fill = lambda value: np.full(grid.ng, value, dtype=float)
    zero = np.zeros(grid.ng)
    columns = np.zeros((grid.nf, 6), dtype=complex)
    return coupling.CauchyState(fill(v["Q"]), fill(v["radius"]), zero.copy(),
                               fill(v["p_Q"]), fill(v["p_r"]), fill(v["p_chi"]),
                               columns.copy(), columns.copy())


def _max(values):
    return float(np.max(np.abs(values)))


def actual_metric(grid, state, rate, bundle):
    """Metric jets from the actual projected discrete rate; chi is not inserted."""
    fine = bundle["fine_state"]
    lift = lambda values: galerkin.prolong_geometry(grid, values)
    project = lambda values: lift(galerkin.pull_geometry(grid, values))
    qt, rt = lift(rate.Q), lift(rate.r)
    prt, pct = lift(rate.p_r), lift(rate.p_chi)
    Fr, fchi = action.partial_F(fine.r, grid.fine.A, grid.fine.C_W)
    qtt = project((qt * fine.p_chi + fine.Q * pct) / (2 * fchi))
    Frt = -8 * np.pi * grid.fine.A * rt
    rtt = project((prt - (Frt * fine.p_chi + Fr * pct) / fchi)
                  / (2 * action.feedback_Z(grid.fine.A)))
    dx = lambda values: tidal.metric.spectral_dx(values, grid.length)
    qx, rx = dx(fine.Q), dx(fine.r)
    package = tidal.curvature_from_local_jets(
        fine.Q, fine.r, qt, qx, qtt, dx(qt), dx(qx),
        rt, rx, rtt, dx(rt), dx(rx),
    )
    return package, {"Q_t": qt, "Q_tt": qtt, "r_t": rt, "r_tt": rtt}


def diagnostics(grid, state, coordinate_time):
    rate, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    fine = bundle["fine_state"]
    v = exact_values(grid, coordinate_time)
    expected = exact_state(grid, coordinate_time)
    rhs = {"Q": v["Q_t"], "r": 0.0, "chi": 0.0, "p_Q": 0.0,
           "p_r": v["p_r_t"], "p_chi": v["p_chi_t"], "phi0": 0.0, "phi1": 0.0}
    package, jets = actual_metric(grid, state, rate, bundle)
    tides, base = package["tides"], package["base"]
    radius = v["radius"]
    invariant_targets = {"R_h": 2.0, "R4": 0.0, "owned_W": 0.0,
                         "Ricci2": 4 / radius**4, "K": 8 / radius**4,
                         "R_0101": -1 / radius**2, "R_0202": 0.0}
    invariant_arrays = {name: tides[name] if name in tides else base[name]
                        for name in invariant_targets}
    Fr, fchi = action.partial_F(fine.r, grid.fine.A, grid.fine.C_W)
    Pi = fine.p_r - Fr * fine.p_chi / fchi
    D = fine.p_chi / (2 * fchi)
    constraints = galerkin.constraint_diagnostics(
        grid, state, source=bundle["source"], fine_state=fine,
    )
    covariance = (np.vstack((state.phi0, state.phi1)) * grid.fine.occupations[None, :])
    covariance = covariance @ np.vstack((state.phi0, state.phi1)).conj().T
    return {
        "coordinate_time": float(coordinate_time), "nf": grid.nf,
        "state_errors": {name: _max(getattr(state, name) - getattr(expected, name)) for name in FIELDS},
        "rate_residuals": {name: _max(getattr(rate, name) - rhs[name]) for name in FIELDS},
        "Q": float(np.mean(fine.Q)), "r": float(np.mean(fine.r)),
        "Q_min": float(np.min(fine.Q)), "r_min": float(np.min(fine.r)),
        "D_square_plus_Q_square_error": _max(D**2 + fine.Q**2 - v["Q0"]**2),
        "Pi_max": _max(Pi), "energy": float(galerkin.energy(grid, state)),
        "constraints": constraints,
        "source_force_max": max(_max(bundle["source"][name]) for name in ("force_L", "force_Q", "force_beta")),
        "fieldwork_power": float(rate.fieldwork_power),
        "covariance_max": _max(covariance),
        "covariance_eigenvalues": [0.0] * (2 * grid.nf) if _max(covariance) == 0 else np.linalg.eigvalsh(covariance).tolist(),
        "metric_invariants": {
            name: {"mean": float(np.mean(values)), "max_error_from_exact": _max(values - invariant_targets[name])}
            for name, values in invariant_arrays.items()
        },
        "metric_time_jet_errors": {name: _max(values - v.get(name, 0.0)) for name, values in jets.items()},
        "proper_lengths": {"period": float(grid.dx_q * np.sum(fine.Q * fine.r)),
                           "child": float(2 * np.mean(fine.Q * fine.r))},
        "exact_normal_proper_time": v["normal_proper_time"],
        "chi_substituted_for_curvature": False,
    }


def rk4_step(grid, state, dt):
    """Independent RK4 composition of the original full discrete rates."""
    add = lambda base, rate, scale: coupling.CauchyState(
        *(getattr(base, name) + scale * getattr(rate, name) for name in FIELDS)
    )
    k1 = galerkin.rates(grid, state)
    k2 = galerkin.rates(grid, add(state, k1, dt / 2))
    k3 = galerkin.rates(grid, add(state, k2, dt / 2))
    k4 = galerkin.rates(grid, add(state, k3, dt))
    return coupling.CauchyState(*(
        getattr(state, name) + dt / 6 * (
            getattr(k1, name) + 2 * getattr(k2, name) + 2 * getattr(k3, name) + getattr(k4, name)
        ) for name in FIELDS
    ))


def _snapshot(payload, grid, state, ordinal):
    prefix = f"nf{grid.nf}_station{ordinal}_"
    for name in FIELDS:
        payload[prefix + name] = np.array(getattr(state, name), copy=True)
    for name in ("p_Q", "p_r", "p_chi"):
        payload[prefix + "canonical_pi_" + name[2:]] = grid.dx_g * getattr(state, name)
    payload[f"nf{grid.nf}_dx_g"] = np.array(grid.dx_g)
    payload[f"nf{grid.nf}_times"] = np.array(RUN_TIMES)


def run_control():
    """Bounded direct residual checks and two short original-rate integrations."""
    started = time.process_time()
    before = source_hashes()
    input_before = {str(coupling.LOCKED_RECORD): sha256_file(coupling.LOCKED_RECORD)}
    payload, exact, numerical = {}, {}, {}
    with threadpool_limits(limits=1):
        for nf in FERMION_COUNTS:
            grid = control_grid(nf)
            exact[str(nf)] = [diagnostics(grid, exact_state(grid, t), t) for t in EXACT_TIMES]
            state = exact_state(grid, 0.0)
            rows = [diagnostics(grid, state, 0.0)]
            _snapshot(payload, grid, state, 0)
            steps = 0
            for ordinal, target in enumerate(RUN_TIMES[1:], 1):
                left = RUN_TIMES[ordinal - 1]
                count = int(round((target - left) / STEP_DT))
                for _ in range(count):
                    if time.process_time() - started >= CPU_BUDGET:
                        raise RuntimeError("source-free control reached its 30 CPU-second budget")
                    state = rk4_step(grid, state, STEP_DT)
                    steps += 1
                rows.append(diagnostics(grid, state, target))
                _snapshot(payload, grid, state, ordinal)
            numerical[str(nf)] = {"step_dt": STEP_DT, "steps": steps, "stations": rows}
        params = parameters(grid)
    after = source_hashes()
    input_after = {path: sha256_file(path) for path in input_before}
    if before != after or input_before != input_after:
        raise RuntimeError("a control source or locked input changed during measurement")
    cpu = time.process_time() - started
    if cpu > CPU_BUDGET:
        raise RuntimeError("source-free control exceeded its 30 CPU-second budget")
    return {
        "schema": SCHEMA, "gauge": "conformal", "metric": "g=(r Q)^2(dt^2-dx^2)-r^2 dOmega^2; signature +---",
        "parameters": params, "locked_coefficients": dict(grid.fine.coefficients),
        "momentum_representation": "native coarse nodal momentum densities p; stored canonical_pi=dx_g*p in the identity frame",
        "source": {"columns": "six zero spinor columns", "occupations": [0.0]*6,
                   "covariance": "C=0; all eigenvalues zero; 0<=C<=I",
                   "occupied_rank6_Gram_identity_required": False,
                   "absolute_continuum_vacuum_claim": False, "sea_subtracted": False,
                   "magnetic_and_induced_action_retained": True, "multiplicity": grid.fine.multiplicity},
        "scope": {"same_initial_geometry_as_discovery_family": False,
                  "comparison": "conditional source-free control with independently prepared locked magnetic radius",
                  "initial_radius_solver_called": False, "coefficients_changed": False,
                  "constraint_projection": False, "singularity_declared": False,
                  "family_instability_decided": False,
                  "interpretation": "Q and proper lengths contract while r and actual 4D tides remain constant. Proper-length contraction alone does not establish physical instability."},
        "exact_formula": "r=r_m; chi=p_Q=0; Q=Q0 sech(Q0 t); D=-Q0 tanh(Q0 t); p_chi=4 alpha D; p_r=-16 pi A r_m D",
        "exact": exact, "numerical": numerical, "cpu_seconds": cpu, "cpu_budget_seconds": CPU_BUDGET,
        "numerical_threads": 1, "source_hashes": after, "source_hashes_unchanged": before == after,
        "input_hashes": input_after, "input_hashes_unchanged": input_before == input_after,
        "python_array_versions": {"numpy": np.__version__},
    }, payload


def output_directory(path=OUTPUT):
    directory = Path(path).expanduser().resolve()
    allowed = OUTPUT.resolve()
    if directory != allowed and allowed not in directory.parents:
        raise ValueError("output must be the vacuum-control-v1 directory or its descendant")
    return directory


def write_record(report, arrays, directory=OUTPUT):
    """Creation only, bounded JSON and NPZ, after the completed bounded run."""
    directory = output_directory(directory)
    if directory.exists():
        raise FileExistsError("refusing to overwrite " + str(directory))
    stream = io.BytesIO()
    np.savez_compressed(stream, **arrays)
    payload = stream.getvalue()
    if len(payload) > MAX_PAYLOAD_BYTES:
        raise ValueError("control payload exceeded its one-MiB bound")
    record = dict(report, payload_sha256=hashlib.sha256(payload).hexdigest(),
                  payload_bytes=len(payload), payload="payload.npz")
    encoded = (json.dumps(record, indent=2, allow_nan=False) + "\n").encode()
    if len(encoded) > MAX_PAYLOAD_BYTES:
        raise ValueError("control JSON exceeded its one-MiB bound")
    directory.mkdir(parents=True, exist_ok=False)
    with (directory / "payload.npz").open("xb") as dest:
        dest.write(payload)
    with (directory / "record.json").open("xb") as dest:
        dest.write(encoded)
    return record


def check_record(directory=OUTPUT):
    """Read-only byte and covariance authentication; no rate evaluation."""
    directory = output_directory(directory)
    record = json.loads((directory / "record.json").read_text())
    if record.get("schema") != SCHEMA:
        raise ValueError("unexpected source-free control schema")
    for group in ("source_hashes", "input_hashes"):
        if not record.get(group):
            raise ValueError("record lacks " + group)
        for path, digest in record[group].items():
            if sha256_file(path) != digest:
                raise ValueError("bound hash changed: " + path)
    payload = directory / "payload.npz"
    if sha256_file(payload) != record["payload_sha256"]:
        raise ValueError("control payload hash changed")
    if payload.stat().st_size > MAX_PAYLOAD_BYTES:
        raise ValueError("control payload exceeded its size bound")
    with np.load(payload, allow_pickle=False) as arrays:
        for name in arrays.files:
            values = arrays[name]
            if not np.isfinite(values).all():
                raise ValueError("nonfinite control array: " + name)
            if name.endswith(("_phi0", "_phi1")) and np.any(values != 0):
                raise ValueError("source-free control has nonzero spinor columns")
    return {"ok": True, "schema": SCHEMA, "evolved": False, "bytes_written": 0,
            "checked_source_hashes": len(record["source_hashes"]),
            "checked_input_hashes": len(record["input_hashes"])}
