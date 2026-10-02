"""Uniform initial-rate/auxiliary controls of the SAME conformal action.

Only initial Cauchy data change. The actual PREPARED-v2 field/eigenframe,
occupations, observer, period, angular multiplicity and action stay fixed.
This is not an inhomogeneous h solver, a pump, or a static-force fit.
Production uses the existing dynamic-preparation episode adapter/pool.
"""
from __future__ import annotations

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import subprocess
import time
import types

import numpy as np
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

from . import nsc_discovery_dynamic_preparation as prep
from . import nsc_discovery_episode as episode
from . import nsc_discovery_extent as extent
from . import nsc_discovery_tidal as tidal
from . import nsc_discovery_response as response
from . import nsc_nested_parent_child as nested
from . import nsc_spherical_coupling as coupling
from . import nsc_spherical_galerkin_coupling as galerkin

LAB = Path(__file__).resolve().parents[2]
ROOT = LAB.parent
SCHEMA = "NSC-DISCOVERY-INITIAL-SECTOR-PREPARATION-v1"
EPISODE_SCHEMA = "NSC-DISCOVERY-INITIAL-SECTOR-EPISODE-v1"
BALANCED_SCHEMA = "NSC-DISCOVERY-INITIAL-SECTOR-PREPARATION-balanced-v2"
BALANCED_EPISODE_SCHEMA = "NSC-DISCOVERY-INITIAL-SECTOR-EPISODE-balanced-v2"
PREPARED = LAB / "results/development/nsc-discovery-dynamic-preparation-v2.json"
COMPARISON = LAB / "results/development/nsc-discovery-dynamic-episode-v1"
DEFAULT_PREPARATION = LAB / "results/development/nsc-discovery-initial-sector-preparation-v1"
DEFAULT_EPISODE = LAB / "results/development/nsc-discovery-initial-sector-v1"
BALANCED_PREPARATION = LAB / "results/development/nsc-discovery-initial-sector-preparation-balanced-v2"
BALANCED_EPISODE = LAB / "results/development/nsc-discovery-initial-sector-balanced-v2"
SECTORS = {"h_minus": (-.5, 0.), "h_plus": (.5, 0.), "chi_minus2": (0., -2.)}
BALANCED_NAMES = ("chi_balanced", "chi_lower", "chi_upper")
CAPS = (.001, .0005)
STATIONS = (.3, 1., 3.)
CPU_BUDGET = 300.
FORECAST_FACTOR = 1.5
OWNERS = tuple("lab/" + name for name in (
    "src/recursive_horizons/nsc_discovery_initial_sector.py",
    "scripts/derive_nsc_discovery_initial_sector.py",
    "tests/test_nsc_discovery_initial_sector.py",
    "docs/nsc-discovery-initial-sector.md"))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dependency_hashes():
    paths = set(prep.dependency_hashes()) | set(OWNERS)
    return {path: sha256(ROOT / path) for path in sorted(paths)}


def uniform_radius_square(system, Q, rho, chi):
    """C=0 with zero pQ and Pi, including the actual auxiliary potential."""
    alpha = float(coupling.alpha_of(system.C_W))
    magnetic = 2 * np.pi * system.C_F * system.flux ** 2
    return (rho / Q + magnetic + alpha * (4 * chi + chi ** 2)) / (8 * np.pi * system.A)


def source_balanced_curvature(pair, base):
    """ONE conditional initial chi_ddot=0 criterion, no dynamic selection."""
    nodal = nested.reconstruct_state(pair, base)
    fine = galerkin.prolong_state(pair.grid, nodal)
    system = galerkin.active_fine_system(pair.grid, fine)
    source = coupling.source_from_columns(system, fine)
    Q = float(np.mean(fine.Q))
    rho_values = source["force_L"] / pair.grid.dx_q
    spin_values = source["S1"].real / pair.grid.dx_q
    rho, spin = float(np.mean(rho_values)), float(np.mean(spin_values))
    if np.max(np.abs(fine.Q - Q)) > 1e-10 or np.max(np.abs(rho_values - rho)) > 1e-8 * max(1., abs(rho)):
        raise ValueError("source-balanced cubic requires the frozen uniform Q/source domain")
    alpha = float(coupling.alpha_of(system.C_W))
    magnetic = float(2 * np.pi * system.C_F * system.flux ** 2)
    M, kappa = int(system.multiplicity), int(system.kappa)
    lower = rho + Q * (magnetic - 16 * alpha / 3)
    if alpha <= 0 or Q <= 0 or rho <= 0 or magnetic <= 16 * alpha / 3 or lower <= 0:
        raise ValueError("strict global derivative/positive-source hypothesis is absent")
    coefficients = np.array([alpha * Q, 4 * alpha * Q, rho + magnetic * Q,
                             6 * rho - 3 * M * kappa * spin])
    polynomial = lambda value: float(np.polyval(coefficients, value))
    bound = 1.
    for _iteration in range(64):
        if polynomial(-bound) <= 0 <= polynomial(bound):
            break
        bound *= 2
    else:
        raise ValueError("finite cubic root bracket failed")
    selected = float(brentq(polynomial, -bound, bound, xtol=np.finfo(float).tiny,
                            rtol=4 * np.finfo(float).eps))
    radius_square = float(uniform_radius_square(system, Q, rho, selected))
    if not np.isfinite(radius_square) or radius_square <= 0:
        raise ValueError("conditional cubic has no positive constraint-selected radius")
    return {"criterion": "uniform initial chi_ddot=0 with zero momenta; fixed Q and Gaussian source",
        "cubic_law": "alpha*Q*chi^3+4*alpha*Q*chi^2+(rho+mag*Q)*chi+6*rho-3*M*kappa*S=0",
        "cubic_coefficients": coefficients.tolist(), "chi_star": selected,
        "cubic_residual": polynomial(selected), "global_derivative_lower_bound": float(lower),
        "derivative_identity": "3*alpha*Q*(chi+4/3)^2+rho+Q*(mag-16*alpha/3)",
        "unique_real_root": True, "Q": Q, "rho": rho, "S": spin, "alpha": alpha,
        "f": 2 * alpha, "mag": magnetic, "M": M, "kappa": kappa,
        "rho_nonuniform_max": float(np.max(np.abs(rho_values - rho))),
        "spin_nonuniform_max": float(np.max(np.abs(spin_values - spin))),
        "constraint_radius_square": radius_square, "positive_radius": True,
        "conditional_preparation": True, "fully_dynamic_source_selection": False,
        "higher_adiabatic_jets_matched": False, "holding_claim": False,
        "perturbation_fraction": .01, "source_eta_changed": False}


def balanced_controls(criterion):
    selected = float(criterion["chi_star"])
    delta = .01 * abs(selected)
    return {"chi_balanced": (0., selected), "chi_lower": (0., selected - delta),
            "chi_upper": (0., selected + delta)}


def _record_controls(record):
    return {name: (float(value["h"]), float(value["chi"])) for name, value in record["controls"].items()}


def initial_measurement(pair, state, *, h, chi):
    """Actual projected rates/metric jets; chi is NEVER a curvature proxy."""
    nodal = nested.reconstruct_state(pair, state)
    rate, bundle = galerkin.compose_fine_hamiltonian(pair.grid, nodal)
    fine, system, source = bundle["fine_state"], bundle["fine_system"], bundle["source"]
    constraints = coupling.constraint_residuals(system, fine, source)
    retained = galerkin.pull_geometry(pair.grid, constraints["hamilton"])
    qdot = galerkin.prolong_geometry(pair.grid, rate.Q)
    rdot = galerkin.prolong_geometry(pair.grid, rate.r)
    chidot = galerkin.prolong_geometry(pair.grid, rate.chi)
    lapse = fine.r * fine.Q
    Kr = (rdot / fine.r + qdot / fine.Q) / lapse
    Kperp = rdot / (lapse * fine.r)
    Fr, f = coupling.partial_F(fine.r, system.A, system.C_W)
    f = float(f)
    Z = float(coupling.feedback_Z(system.A))
    Pi = fine.p_r - Fr * fine.p_chi / f
    first = pair.grid.dx_q * np.sum(fine.Q * fine.p_Q * fine.p_chi / (2 * f))
    second = pair.grid.dx_q * np.sum(Pi ** 2 / (4 * Z))
    columns = np.vstack((state.phi0, state.phi1))
    gram = columns.conj().T @ columns
    roots = np.sqrt(pair.weights)
    eigen = np.linalg.eigvalsh(roots[:, None] * gram * roots[None, :])
    row, _profiles = extent.observe(pair, state, 0., np.zeros(3), "coupled")
    jets = tidal.analytic_accelerations(pair, state, "coupled")
    tangent = response.StateTangent(*(getattr(rate, key) for key in nested.STATE_NAMES), np.zeros(6))
    directional = response.rate_jacobian_vector(pair.grid, nodal, tangent)
    chi_ddot = galerkin.prolong_geometry(pair.grid, directional.chi)
    Qmean = float(np.mean(fine.Q)); rmean = float(np.mean(fine.r))
    rhobar = float(np.mean(constraints["rho"]))
    spinbar = float(np.mean(source["S1"].real / pair.grid.dx_q))
    alpha = f / 2
    chi_ddot_formula = (2 * Qmean * rhobar - system.multiplicity * system.kappa * Qmean * spinbar
        + (8 / 3) * np.pi * system.A * rmean ** 2 * Qmean ** 2 * chi) / (4 * alpha)
    summary = lambda value: {"min": float(np.min(value)), "max": float(np.max(value)),
                             "max_abs": float(np.max(np.abs(value)))}
    report = {"h": float(h), "chi_initial": float(chi),
        "Q": summary(fine.Q), "r": summary(fine.r), "rho": summary(constraints["rho"]),
        "source_current_mean": float(np.mean(constraints["current"])),
        "shift_residual_mean": float(np.mean(constraints["momentum"])),
        "current_mean_deleted": False, "retained_lapse_max": float(np.max(np.abs(retained))),
        "full_lapse_max": constraints["hamilton_max"], "full_shift_max": constraints["momentum_max"],
        "held_out_lapse_max": float(np.max(np.abs(constraints["hamilton"] - pair.grid.A_g @ retained))),
        "Qdot": summary(qdot), "rdot": summary(rdot), "chidot": summary(chidot),
        "momentum_rates": {name: summary(galerkin.prolong_geometry(pair.grid, getattr(rate, name)))
                           for name in nested.MOMENTUM_NAMES},
        "actual_Qtt": summary(jets["Q_tt"]), "actual_rtt": summary(jets["r_tt"]),
        "actual_chi_ddot": dict(summary(chi_ddot), mean=float(np.mean(chi_ddot)),
                                  owner="analytic full-state rate_jacobian_vector"),
        "uniform_chi_ddot_formula": float(chi_ddot_formula),
        "chi_ddot_formula_native_gap_max": float(np.max(np.abs(chi_ddot - chi_ddot_formula))),
        "K_radial": summary(Kr), "K_angular": summary(Kperp),
        "Pi_max": float(np.max(np.abs(Pi))), "gram_max": float(np.max(np.abs(gram - np.eye(6)))),
        "CAR_nonzero_spectrum": eigen.tolist(), "CAR_min": float(min(0., np.min(eigen))),
        "CAR_max": float(np.max(eigen)), "source_trace": float(np.dot(pair.weights, np.diag(gram).real)),
        "multiplicity": int(system.multiplicity), "kappa": int(system.kappa),
        "source_coordinate_energy": row["total_energy"]["field"],
        "gravity_energy": row["total_energy"]["gravity"], "total_energy": row["total_energy"]["total"],
        "gravity_kinetic_budget": {"pQ_pchi": float(first), "Pi_squared": float(second), "total": float(first + second)},
        "actual_metric_curvature": row["actual_tides"], "chi_substituted_for_curvature": False,
        "initial_observation": row, "future_rates_fixed": False,
        "continuum_initial_data_certified": False, "holding_claim": False}
    arrays = {"fine_hamilton": constraints["hamilton"], "fine_momentum": constraints["momentum"],
              "fine_rho": constraints["rho"], "fine_current": constraints["current"],
              "fine_chi_ddot": chi_ddot}
    return report, arrays


def prepare_uniform_sector(pair, base, name, *, cpu_limit=10., balanced_criterion=None):
    """Use the frozen source. Only the declared initial geometry sector changes."""
    controls = SECTORS if balanced_criterion is None else balanced_controls(balanced_criterion)
    if name not in controls or not 0 < cpu_limit <= 30:
        raise ValueError("one declared sector and preparation CPU in (0,30] required")
    started = time.process_time()
    nodal = nested.reconstruct_state(pair, base)
    Q0 = float(np.mean(nodal.Q))
    if np.max(np.abs(nodal.Q - Q0)) > 1e-10 or np.max(np.abs(nodal.r - np.mean(nodal.r))) > 1e-9:
        raise ValueError("this bounded adapter requires the uniform prepared case")
    if np.any(nodal.chi != 0) or any(np.any(getattr(base, name) != 0) for name in nested.MOMENTUM_NAMES):
        raise ValueError("baseline must be the frozen chi0 zero-momentum preparation")
    h, chi = controls[name]
    nodal.chi = np.full(pair.grid.ng, chi)
    correction = {"performed": False, "reason": "chi0 intrinsic source-selected geometry retained"}
    if chi != 0:
        fine = galerkin.prolong_state(pair.grid, nodal)
        system = galerkin.active_fine_system(pair.grid, fine)
        rho = coupling.source_from_columns(system, fine)["force_L"] / pair.grid.dx_q
        selected = float(uniform_radius_square(system, Q0, float(np.mean(rho)), chi))
        if not np.isfinite(selected) or selected <= 0:
            raise ValueError("actual auxiliary/source constraint has no positive uniform radius")
        nodal.r = np.full(pair.grid.ng, np.sqrt(selected))
        def oracle(radius):
            trial = nodal.copy(); trial.r = radius
            residual, _full, jacobian = prep.radius_residual_jacobian(pair.grid, trial, rho)
            return residual, jacobian
        nodal.r, solved = prep._positive_newton(nodal.r, oracle,
            lambda values: pair.grid.A_g @ values, started + cpu_limit)
        if not solved["converged"]:
            raise ValueError("auxiliary/source finite radius correction failed: " + str(solved["blocker"]))
        correction = dict(solved, performed=True, analytic_radius_square=selected,
                          old_radius_copied=False, target="owned finite initial lapse constraint")
    a = -8 * np.pi * pair.grid.fine.A
    f = -8 * np.pi * pair.grid.fine.C_W / 3
    nodal.p_Q = np.zeros(pair.grid.ng)
    nodal.p_chi = np.full(pair.grid.ng, 2 * f * h)
    nodal.p_r = 2 * a * nodal.r * h
    encoded = nested.encode_state(pair, nodal)
    if not np.array_equal(encoded.phi0, base.phi0) or not np.array_equal(encoded.phi1, base.phi1):
        raise ValueError("initial-sector construction changed the frozen source")
    measured, arrays = initial_measurement(pair, encoded, h=h, chi=chi)
    measured.update(radius_preparation=correction, cpu_seconds=time.process_time() - started,
                    field_frame_changed=False, observer_changed=False, coefficients_changed=False)
    if balanced_criterion is not None:
        measured["curvature_preparation"] = {
            "criterion_satisfied_by_declaration": name == "chi_balanced",
            "chi_star": balanced_criterion["chi_star"], "chi_offset": chi - balanced_criterion["chi_star"],
            "perturbation_fraction": .01 if name != "chi_balanced" else 0.,
            "nearby_control_not_claimed_balanced": name != "chi_balanced",
            "conditional_initial_criterion_only": True, "future_curvature_fixed": False}
    for key in nested.STATE_NAMES:
        arrays["nodal_" + key] = np.array(getattr(nodal, key), copy=True)
        arrays[key] = np.array(getattr(encoded, key), copy=True)
    return encoded, measured, arrays


def comparison_references(directory=COMPARISON):
    """Bind the completed chi0/h0 uniform comparison, never copy or rerun it."""
    directory = Path(directory).resolve()
    paths = {str(directory / "manifest.json"): sha256(directory / "manifest.json")}
    cases = []
    for cap in CAPS:
        identifier = f"nf128_uniform_dt{cap:g}"
        initial, _ = episode.load_checkpoint(directory, identifier, ordinal=0)
        final, _ = episode.load_checkpoint(directory, identifier)
        if abs(float(final["coordinate_time"]) - 3.) > 1e-12:
            raise ValueError("old uniform comparison did not reach the recorded station")
        for path in sorted(directory.glob(identifier + "-*")):
            if path.suffix in (".json", ".npz", ".jsonl"):
                paths[str(path)] = sha256(path)
        cases.append({"case_id": identifier, "initial_state_sha256": initial["source_pins"]["initial_state_sha256"],
                      "last_time": float(final["coordinate_time"]), "producing_commit": final["producing_commit"]})
    return {"directory": str(directory), "cases": cases, "hashes": paths, "rerun": False,
            "sector": {"h": 0., "chi": 0.}, "producer_preserved": True}


def build_record(initial=PREPARED, *, execute=False, cpu_limit=30., producer_commit=None,
                 comparison=COMPARISON, balanced=False):
    whole_started = time.process_time()
    if not 0 < cpu_limit <= 30:
        raise ValueError("whole development/preparation CPU limit must be in (0,30]")
    binding = prep.authenticate_preparation(initial)
    references = comparison_references(comparison)
    criterion = None; controls = SECTORS
    if balanced:
        with threadpool_limits(limits=1):
            pre_pair, pre_base, _old = prep.load_case(initial, "uniform")
            criterion = source_balanced_curvature(pre_pair, pre_base)
        controls = balanced_controls(criterion)
    report = {"schema": BALANCED_SCHEMA if balanced else SCHEMA,
        "status": "SOURCE_BALANCED_CURVATURE_PREFLIGHT" if balanced else "INITIAL_SECTOR_PREFLIGHT", "evidence_written": False,
        "initial_binding": binding, "initial_record": str(prep.output_paths(initial)[0]),
        "source_case": "uniform", "source_hashes": dependency_hashes(), "cases": {},
        "controls": {name: {"h": values[0], "chi": values[1]} for name, values in controls.items()},
        "old_comparison": references, "stations": list(STATIONS), "step_caps": list(CAPS),
        "initial_h_is_declared_input": True, "future_rates_fixed": False,
        "source_reselected": False, "observer_reset": False, "coefficients_changed": False,
        "mass_changed": False, "evolution_performed": False, "continuum_initial_data_certified": False,
        "cpu_limit_seconds": float(cpu_limit)}
    if balanced:
        report.update(curvature_criterion=criterion, family="balanced-v2",
            initial_zero_momenta_declared=True,
            route_policy="if all three show curvature growth or chart loss, stop this initial-sector route; no higher-jet or alpha tuning",
            physical_curvature_selection_claim=False)
    if producer_commit is not None:
        report["producing_commit"] = prep._git_hashes(producer_commit, report["source_hashes"])
    if not execute:
        return report, {}
    if producer_commit is None:
        raise ValueError("new production preparation requires a frozen producing commit")
    start = whole_started; arrays = {}
    with threadpool_limits(limits=1):
        pair, base, old = prep.load_case(initial, "uniform")
        report.update(nf=pair.grid.nf, ng=pair.grid.ng, nq=pair.grid.nq, period=pair.grid.length,
                      gauge="conformal", occupations=pair.weights.tolist(), source_eta=old["source_eta"],
                      multiplicity=pair.grid.fine.multiplicity, locked_coefficients=pair.grid.fine.coefficients)
        report["baseline_prepared_state_sha256"] = episode.state_sha256(base)
        if any(case["initial_state_sha256"] != report["baseline_prepared_state_sha256"]
               for case in references["cases"]):
            raise ValueError("saved h0/chi0 comparison is not the same prepared uniform Cauchy state")
        report["baseline_initial_measurement"] = initial_measurement(pair, base, h=0., chi=0.)[0]
        for name in controls:
            remaining = cpu_limit - (time.process_time() - start)
            if remaining <= 0:
                raise TimeoutError("initial-sector preparation CPU cap")
            _state, measured, saved = prepare_uniform_sector(pair, base, name, cpu_limit=remaining,
                                                            balanced_criterion=criterion)
            report["cases"][name] = measured
            arrays.update({name + "_" + key: value for key, value in saved.items()})
        arrays.update(W=pair.geometry_map, source_phi0=pair.source_phi0, source_phi1=pair.source_phi1,
                      observer_columns=pair.reference_columns, source_weights=pair.weights)
        report["basis_indices"] = {"coarse": pair.geometry_coarse_indices.tolist(),
                                   "child": pair.geometry_child_indices.tolist(),
                                   "parent": pair.geometry_parent_indices.tolist()}
        report["clock_locations"] = list(pair.clock_locations)
    report.update(status="FINITE_INITIAL_SECTOR_DATA", cpu_seconds=time.process_time() - start)
    if report["source_hashes"] != dependency_hashes():
        raise ValueError("producer source closure changed during preparation")
    return report, arrays


def write_record(report, arrays, path=None):
    if path is None:
        path = BALANCED_PREPARATION if report["schema"] == BALANCED_SCHEMA else DEFAULT_PREPARATION
    return prep.write_record(report, arrays, path)


def load_sector(path, name):
    jp, payload = prep.output_paths(path); record = json.loads(jp.read_text())
    if record.get("schema") not in (SCHEMA, BALANCED_SCHEMA) or sha256(payload) != record["payload_sha256"] or name not in record["controls"]:
        raise ValueError("initial-sector record/payload mismatch")
    with np.load(payload, allow_pickle=False) as saved:
        arrays = {key: np.array(saved[key], copy=True) for key in saved.files}
    nf = int(record["nf"])
    chunk_arrays = {"W": arrays["W"], "source_weights": arrays["source_weights"],
        "observer_columns": arrays["observer_columns"], "source_phi0": arrays["source_phi0"], "source_phi1": arrays["source_phi1"]}
    indices = record["basis_indices"]
    pair = episode.pair_from_arrays(chunk_arrays, dict(nf=nf, coarse_indices=indices["coarse"],
        child_indices=indices["child"], parent_indices=indices["parent"], clock_locations=record["clock_locations"],
        source_metadata={"source": "frozen PREPARED-v2 uniform actual eigenframe", "eta": 1.},
        geometry_metadata={"one_metric": True, "all_geometry_degrees_retained": True}))
    state = nested.NestedState(*(arrays[name + "_" + key] for key in nested.STATE_NAMES))
    return pair, state, record


def authenticate_record(path):
    jp, payload = prep.output_paths(path); record = json.loads(jp.read_text())
    if record.get("schema") not in (SCHEMA, BALANCED_SCHEMA) or sha256(payload) != record["payload_sha256"]:
        raise ValueError("initial-sector binding failed")
    prep._git_hashes(record["producing_commit"], record["source_hashes"])
    old = prep.authenticate_preparation(record["initial_record"])
    if old != record["initial_binding"]:
        raise ValueError("frozen PREPARED-v2 input changed")
    for filename, expected in record["old_comparison"]["hashes"].items():
        if sha256(filename) != expected:
            raise ValueError("old chi0/h0 comparison changed")
    return record


def _historical_v1_check(path, record):
    """Execute the recorded caf producer's MODULE bytes, without rewriting data."""
    source_path = OWNERS[0]
    commit = record["producing_commit"]
    source = subprocess.check_output(["git", "show", commit + ":" + source_path], cwd=ROOT)
    expected = record["source_hashes"][source_path]
    if hashlib.sha256(source).hexdigest() != expected:
        raise ValueError("historical initial-sector producer bytes differ")
    for filename, wanted in record["source_hashes"].items():
        if filename not in OWNERS and sha256(ROOT / filename) != wanted:
            raise ValueError("historical physical dependency changed: " + filename)
    module = types.ModuleType("recursive_horizons._initial_sector_historical")
    module.__package__ = "recursive_horizons"
    module.__file__ = str(ROOT / source_path)
    exec(compile(source, module.__file__, "exec"), module.__dict__)
    result = module.check_record(path)
    return dict(result, historical_producer_source_loaded=True, historical_source_commit=commit,
                historical_source_sha256=expected, live_checker_source_used=False)


def check_record(path):
    record = authenticate_record(path)
    if record["schema"] == SCHEMA and record["source_hashes"][OWNERS[0]] != sha256(ROOT / OWNERS[0]):
        return _historical_v1_check(path, record)
    for name, declared in _record_controls(record).items():
        pair, state, _ = load_sector(path, name)
        measured, _ = initial_measurement(pair, state, h=declared[0], chi=declared[1])
        for key in ("full_lapse_max", "full_shift_max", "source_coordinate_energy", "source_trace"):
            if abs(measured[key] - record["cases"][name][key]) > 1e-8:
                raise ValueError("saved initial-sector measurement changed: " + name + "/" + key)
        if record["schema"] == BALANCED_SCHEMA:
            old = record["cases"][name]
            if abs(measured["actual_chi_ddot"]["mean"] - old["actual_chi_ddot"]["mean"]) > 1e-7:
                raise ValueError("saved balanced native acceleration changed")
            for channel in ("R_h", "owned_W"):
                if abs(measured["actual_metric_curvature"][channel]["max_abs"] -
                       old["actual_metric_curvature"][channel]["max_abs"]) > 1e-7:
                    raise ValueError("saved actual initial metric curvature changed")
    if record["schema"] == BALANCED_SCHEMA:
        with threadpool_limits(limits=1):
            pair, base, _ = prep.load_case(record["initial_record"], "uniform")
            if source_balanced_curvature(pair, base) != record["curvature_criterion"]:
                raise ValueError("saved source-conditioned cubic does not replay")
    return {"ok": True, "bytes_written": 0, "evolution_performed": False,
            "initial_source_reselected": False, "cases": list(record["controls"])}


def prepare_episode(initial, output=None, *, execute=False, producer_commit=None,
                    resolved_error=1e-6):
    started = time.process_time(); saved = authenticate_record(initial)
    controls = _record_controls(saved)
    balanced = saved["schema"] == BALANCED_SCHEMA
    if output is None:
        output = BALANCED_EPISODE if balanced else DEFAULT_EPISODE
    if not np.isfinite(resolved_error) or resolved_error <= 0:
        raise ValueError("existing initial resolved-error domain must be positive finite")
    hashes = dependency_hashes()
    report = {"schema": BALANCED_EPISODE_SCHEMA if balanced else EPISODE_SCHEMA, "status": "EPISODE_PREFLIGHT", "cases": [],
        "source_hashes": hashes, "initial_record": str(prep.output_paths(initial)[0]),
        "initial_binding": {"record_sha256": sha256(prep.output_paths(initial)[0]),
                            "payload_sha256": sha256(prep.output_paths(initial)[1])},
        "old_comparison": saved["old_comparison"], "sectors": list(controls), "stations": list(STATIONS),
        "workers_maximum": 4, "cpu_budget_seconds": CPU_BUDGET, "forecast_factor": FORECAST_FACTOR,
        "initial_constraints_resolved_error": float(resolved_error), "evolved": False,
        "future_rates_fixed": False, "source_reset": False, "observer_reset": False}
    if balanced:
        report.update(curvature_criterion=saved["curvature_criterion"], route_policy=saved["route_policy"],
                      conditional_initial_curvature_criterion_only=True, higher_adiabatic_jets_matched=False)
    if not execute:
        return report
    if producer_commit is None:
        raise ValueError("episode creation requires a frozen producing commit")
    report["producing_commit"] = prep._git_hashes(producer_commit, hashes)
    directory = episode.assert_campaign_output(output)
    if episode._manifest_path(directory).exists():
        raise FileExistsError("initial-sector episode already exists")
    loaded = {}
    for name in controls:
        pair, state, _ = load_sector(initial, name)
        measured, _ = initial_measurement(pair, state, h=controls[name][0], chi=controls[name][1])
        if measured["full_lapse_max"] > resolved_error or measured["full_shift_max"] > resolved_error:
            raise ValueError("owned initial constraints exceed the existing declared domain: " + name)
        if measured["CAR_min"] < -1e-9 or measured["CAR_max"] > 1 + 1e-9:
            raise ValueError("prepared Gaussian source left CAR domain")
        loaded[name] = pair, state
    directory.mkdir(parents=True, exist_ok=True)
    for name, (pair, state) in loaded.items():
        basis = {"W": pair.geometry_map, "source_phi0": pair.source_phi0, "source_phi1": pair.source_phi1,
                 "observer_columns": pair.reference_columns, "source_weights": pair.weights}
        arrays = episode.arrays_from_state(state, basis_arrays=basis,
            clocks={"rates": extent.real_metrics(pair, state)["clock_rates"], "normal_clocks": np.zeros(3)})
        pins = {"W_sha256": episode.array_sha256(pair.geometry_map),
            "source_columns_sha256": episode.array_sha256(pair.source_columns),
            "observer_columns_sha256": episode.array_sha256(pair.reference_columns),
            "weights_sha256": episode.array_sha256(pair.weights), "initial_state_sha256": episode.state_sha256(state)}
        for cap in CAPS:
            identifier = f"nf{pair.grid.nf}_{name}_dt{cap:g}"
            record = {"case_id": identifier, "nf": pair.grid.nf, "coordinate_time": 0., "steps": 0,
                "stations": list(STATIONS), "stations_reached": [], "step_cap": cap,
                "control_mode": "coupled", "geometry": "evolving", "momentum_representation": episode.CANONICAL_PI,
                "source_pins": pins, "source_pins_before": pins, "source_pins_after": pins,
                "coarse_indices": pair.geometry_coarse_indices.tolist(), "child_indices": pair.geometry_child_indices.tolist(),
                "parent_indices": pair.geometry_parent_indices.tolist(), "clock_locations": list(pair.clock_locations),
                "source_metadata": pair.source_metadata, "geometry_metadata": pair.geometry_metadata,
                "source_hashes": hashes, "producing_commit": report["producing_commit"],
                "initial_sector": {"name": name, "h": controls[name][0], "chi": controls[name][1]},
                "initial_measurement": saved["cases"][name], "initial_binding": report["initial_binding"],
                "initial_state_called": False, "sector_initial_preparation_performed": True,
                "verify_external_pins": False, "status": "prepared", "snapshot_kind": "handoff",
                "clock_protocol": "endpoint_trapezoid_per_owned_RK4_step",
                "future_rates_fixed": False, "work_ledger": episode.empty_work_ledger()}
            committed = episode.commit_checkpoint(directory, record, arrays)
            report["cases"].append({key: committed[key] for key in ("case_id", "nf", "npz", "json", "ordinal", "coordinate_time", "steps")})
    charged = float(saved["cpu_seconds"]) + time.process_time() - started
    _, ledger = episode._ledger_update(directory, pilot=charged, budget=CPU_BUDGET)
    report.update(status="PREPARED", stage=0, workers=4, chunk_limit_bytes=episode.CHUNK_LIMIT_BYTES,
        child_cpu_seconds=0., cumulative_cpu_seconds=ledger["spent"],
        adapter_preparation_cpu_seconds=time.process_time() - started,
        observer_sha256=prep.EXTENT_OBSERVER_SHA256, memory_limit_bytes=episode.DEFAULT_MEMORY_BYTES)
    episode._write_json(episode._manifest_path(directory), report)
    return report


def run_episode(output, *, workers=4, cpu_budget=CPU_BUDGET, max_steps=None):
    if int(workers) != workers or not 1 <= workers <= 4 or not 0 < cpu_budget <= CPU_BUDGET:
        raise ValueError("initial-sector episode allows at most four workers and 300 CPU seconds")
    manifest = episode.read_manifest(output)
    if manifest.get("schema") not in (EPISODE_SCHEMA, BALANCED_EPISODE_SCHEMA) or dependency_hashes() != manifest["source_hashes"]:
        raise ValueError("initial-sector episode source closure changed")
    prep._git_hashes(manifest["producing_commit"], manifest["source_hashes"])
    saved = authenticate_record(manifest["initial_record"])
    if sha256(manifest["initial_record"]) != manifest["initial_binding"]["record_sha256"]:
        raise ValueError("initial-sector prepared input changed")
    with prep.episode_adapter():
        result = episode.run(output, workers=int(workers), cpu_budget_seconds=cpu_budget,
            forecast_factor=FORECAST_FACTOR, max_steps=max_steps, executor=prep._episode_pool)
    ledger = json.loads(episode._ledger_path(output).read_text())
    result.update(schema=manifest["schema"], cumulative_cpu_seconds=ledger["spent"],
                  future_rates_fixed=False, source_reset=False, observer_reset=False)
    result["assessment"] = assess_episode(output, result)
    episode._write_json(episode._manifest_path(output), result)
    return result


def _channel(row, name):
    if name == "proper_probability_width":
        return row.get("all_column_localization", {}).get("participation_proper_width")
    return row.get(name)


def assess_episode(output, manifest=None):
    manifest = episode.read_manifest(output) if manifest is None else manifest
    base = prep.assess_episode(output, manifest)
    channels = ("child_proper_length", "child_r_proper_mean", "all_source_child_probability",
                "child_matter_normal_energy", "proper_probability_width")
    comparisons = []; controls = []
    nf = manifest["cases"][0]["nf"]
    old_directory = manifest["old_comparison"]["directory"]
    for path, expected in manifest["old_comparison"]["hashes"].items():
        if sha256(path) != expected:
            raise ValueError("read-only old comparison binding changed")
    for name in manifest.get("sectors", SECTORS):
        identifiers = [f"nf{nf}_{name}_dt{cap:g}" for cap in CAPS]
        series = [episode.read_observations(output, identifier) for identifier in identifiers]
        records = [episode.load_checkpoint(output, identifier)[0] for identifier in identifiers]
        origin = [record["source_pins"]["initial_state_sha256"] for record in records]
        fine = {float(row["time"]): row for row in series[1]}
        for row in series[0]:
            other = fine.get(float(row["time"]))
            if other is None or row["time"] <= 0:
                continue
            gap = {key: _channel(other, key) - _channel(row, key) for key in channels
                   if _channel(other, key) is not None and _channel(row, key) is not None}
            comparisons.append({"sector": name, "time": row["time"], "matched_initial_state": origin[0] == origin[1],
                                "second_minus_first": gap, "finite_timestep_indicator_only": True})
        old = episode.read_observations(old_directory, f"nf{nf}_uniform_dt0.0005")
        old_map = {float(row["time"]): row for row in old}
        if series[1] and old:
            start_new, start_old = series[1][0], old[0]
            for row in series[1]:
                baseline = old_map.get(float(row["time"]))
                if baseline is None or row["time"] <= 0:
                    continue
                effect = {key: (_channel(row, key) - _channel(start_new, key)) -
                               (_channel(baseline, key) - _channel(start_old, key)) for key in channels
                          if all(_channel(item, key) is not None for item in (row, baseline, start_new, start_old))}
                controls.append({"sector": name, "time": row["time"], "initial_offset_removed": True,
                                 "effects_vs_saved_chi0_h0": effect, "proper_clocks_matched": False})
    base.update(matched_step_comparisons=comparisons, saved_baseline_control_effects=controls,
                old_baseline_rerun=False, physical_h_selection_claim=False, holding_claim=False)
    if manifest.get("schema") == BALANCED_EPISODE_SCHEMA:
        evidence = {}
        for name in manifest["sectors"]:
            identifier = f"nf{nf}_{name}_dt0.0005"
            stations = []
            for filename in sorted(Path(output).glob(identifier + "-*.json")):
                if "observations" in filename.name:
                    continue
                record = json.loads(filename.read_text())
                with np.load(filename.with_suffix(".npz"), allow_pickle=False) as stored:
                    chi_nodes = stored["W"] @ stored["chi"]
                excursion = float(np.max(np.abs(chi_nodes - record["initial_sector"]["chi"])))
                stations.append({"time": float(record["coordinate_time"]), "chi_excursion_max": excursion,
                                 "status": record["status"]})
            last, _ = episode.load_checkpoint(output, identifier)
            positive = [row for row in stations if row["time"] > 0]
            # This labels only an observed finite departure trend; not a
            # spectral instability, infinite runaway, or a stability gate.
            growing = len(positive) >= 2 and all(b["chi_excursion_max"] > a["chi_excursion_max"]
                                               for a, b in zip(positive[:-1], positive[1:]))
            chart_loss = last["status"] == "chart_exit"
            evidence[name] = {"retained_stations": stations, "growing_retained_chi_departure": growing,
                              "chart_loss": chart_loss, "physical_instability_proved": False}
        limited = len(evidence) == 3 and all(row["chart_loss"] or row["growing_retained_chi_departure"]
                                           for row in evidence.values())
        base["conditional_curvature_route"] = {
            "evidence": evidence, "initial_criterion_only": True, "sector_limitation_observed": limited,
            "stop_this_route_without_higher_jet_or_alpha_tuning": limited,
            "growth_definition": "strictly increasing |chi-chi_initial| maxima at retained positive-time stations",
            "holding_or_infinite_future_claim": False,
            "scope": "finite conditional initial-sector family; route policy is not a broad theory verdict"}
    return base
