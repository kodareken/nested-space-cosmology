#!/usr/bin/env python3
"""Read saved four-dimensional curvature and normal-clock Dirac frequency scales.

Default prints a report and writes nothing. --write PATH --producer-commit SHA
creates one immutable report. --check PATH authenticates and remeasures it.
No initial data, trajectory, bound campaign, or threshold gate is generated.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import sys
import time

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_name, "1")

import numpy as np
import scipy
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_backend as backend
from recursive_horizons import nsc_discovery_episode as episode
from recursive_horizons import nsc_discovery_tidal as tidal
from recursive_horizons import nsc_nested_parent_child as nested
from recursive_horizons import nsc_spherical_coupling as coupling
from recursive_horizons import provenance

SCHEMA = "NSC-DISCOVERY-EFT-DOMAIN-v1"
LAB = episode.LAB
REPO = episode.REPO
DEFAULT_DYNAMIC = LAB/"results/development/nsc-discovery-dynamic-episode-v1"
DEFAULT_BALANCED = LAB/"results/development/nsc-discovery-initial-sector-balanced-v2"
ORDINAL_STATIONS = (0., .3, 1., 3.)
CASES = {"dynamic": ("uniform", "pattern", "coherent"),
         "balanced": ("chi_balanced", "chi_lower", "chi_upper")}


def digest_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def source_hashes():
    """Transitive local Python import closure plus the mathematical and test owners."""
    package = LAB/"src/recursive_horizons"
    pending = [Path(__file__).resolve(), LAB/"tests/test_nsc_discovery_eft_domain.py",
               LAB/"docs/nsc-discovery-eft-domain.md", LAB/"docs/nsc-curvature-eft.md", package/"__init__.py"]
    pending += [Path(module.__file__).resolve() for module in
                (backend, episode, tidal, nested, coupling, provenance)]
    closure = set()
    while pending:
        path = pending.pop().resolve()
        if path in closure:
            continue
        closure.add(path)
        if path.suffix != ".py":
            continue
        for node in ast.walk(ast.parse(path.read_text())):
            modules = []
            if isinstance(node, ast.ImportFrom):
                module = node.module or ""
                if node.level and path.parent == package:
                    modules = [module.split(".")[0]] if module else [item.name for item in node.names]
                elif module == "recursive_horizons":
                    modules = [item.name for item in node.names]
                elif module.startswith("recursive_horizons."):
                    modules = [module.split(".")[1]]
            elif isinstance(node, ast.Import):
                modules = [item.name.split(".")[1] for item in node.names
                           if item.name.startswith("recursive_horizons.")]
            for module in modules:
                candidate = package/(module+".py")
                if candidate.is_file():
                    pending.append(candidate)
    return {path.relative_to(REPO).as_posix(): digest_file(path) for path in sorted(closure)}


def authenticate_sources(commit, hashes):
    """Read-only immutable Git-blob authentication through the existing provenance owner."""
    for path, digest in hashes.items():
        provenance.resolve_pinned_source_bytes(REPO, path, digest, commit=commit)
    return {"commit": commit, "source_files": len(hashes), "authenticated": True}


def summary(values):
    array = np.asarray(values, dtype=float)
    if not np.isfinite(array).all():
        raise ValueError("nonfinite physical diagnostic")
    peak = float(np.max(abs(array)))
    return {"min": float(np.min(array)), "max": float(np.max(array)), "max_abs": peak,
            "rms": 0. if peak == 0 else peak*float(np.sqrt(np.mean((array/peak)**2)))}


def physical_curvature(profiles, cw_over_a, index):
    """Keep signed Lorentz contractions and observer components separate, all in physical units."""
    base, tides = profiles["base"], profiles["tides"]
    riemann = np.maximum.reduce([abs(base[key]) for key in ("A", "B00", "B11", "B01", "S")])
    normal_tide = np.maximum(abs(tides["R_0101"]), abs(tides["R_0202"]))
    ricci = np.maximum.reduce([abs(base["A"]-2*base["B00"]), abs(-base["A"]-2*base["B11"]),
                              abs(-2*base["B01"]), abs(base["S"]+base["B00"]-base["B11"])])
    actual = {"R4": tides["R4"], "Ricci2": base["Ricci2"], "K": base["K"], "Weyl2": tides["owned_W"],
              "R_0101": tides["R_0101"], "R_0202": tides["R_0202"],
              "Riemann_components_max": riemann, "normal_tidal_components_max": normal_tide,
              "Ricci_components_max": ricci}
    scales = {"abs_R4": abs(tides["R4"]), "Riemann_components": riemann, "normal_tidal_components": normal_tide,
              "Ricci_components": ricci, "sqrt_abs_Ricci2": np.sqrt(abs(base["Ricci2"])),
              "sqrt_abs_K": np.sqrt(abs(base["K"])), "sqrt_Weyl2": np.sqrt(abs(tides["owned_W"]))}
    components = {"R0101": -base["A"], "R0202_R0303": base["B00"], "R1212_R1313": base["B11"],
                  "R0212_R0313": base["B01"], "R2323": -base["S"]}
    numerator = float(np.max(base["R4_condition_numerator"]))
    r4_peak = float(np.max(abs(tides["R4"])))
    return {"physical_values": {key: summary(value) for key, value in actual.items()},
            "worldline_values": {key: float(value[index]) for key, value in actual.items()},
            "tetrad_components": {key: summary(value) for key, value in components.items()},
            "cw_over_a_times_physical_scale": {key: cw_over_a*float(np.max(value)) for key, value in scales.items()},
            "worldline_cw_over_a_times_scale": {key: cw_over_a*float(value[index]) for key, value in scales.items()},
            "R4_cancellation_indicator": None if r4_peak == 0 else numerator/r4_peak,
            "scalar_route_max_gap": float(np.max(abs(base["R4"]-tides["R4"]))),
            "Weyl_large_invariant_subtraction_max_gap": float(np.max(abs(base["assembled_W"]-tides["owned_W"]))),
            "chi_substituted": False, "base_A_is_R2_over_2_not_action_A": True,
            "invariant_squares_are_positive_norms": False, "indicators_are_error_bounds": False}


def frequency_scales(pair, state, jets, pole_mass, index):
    """Instantaneous frozen-generator shells, angular mass, and actual source-state variation."""
    fine, grid = jets["fine"], pair.grid
    normal_clock = fine.r*fine.Q
    hamiltonian = nested.hamiltonian(pair, state)
    hermitian_gap = float(np.max(abs(hamiltonian-hamiltonian.conj().T)))
    if hermitian_gap > 1e-9*max(1., float(np.max(abs(hamiltonian)))):
        raise ValueError("instantaneous representative Dirac generator is not Hermitian")
    # sigma3 anticommutes with the conformal Dirac block. The first nf levels are negative.
    levels = eigh(hamiltonian, eigvals_only=True, subset_by_index=(grid.nf, grid.nf+5))
    if np.min(levels) <= 0:
        raise ValueError("six lowest strictly positive instantaneous Dirac levels are unavailable")
    proper_shells = levels[:, None]/normal_clock[None, :]
    rate = jets["rate"]  # actual nodal projected Cauchy rate, not a reselected state
    phi0_t, phi1_t = grid.U_f@rate.phi0, grid.U_f@rate.phi1
    mass = np.sum((abs(fine.phi0)**2+abs(fine.phi1)**2)*pair.weights, axis=1)
    numerator = np.sum((abs(phi0_t)**2+abs(phi1_t)**2)*pair.weights, axis=1)
    supported = mass > 0
    if not np.any(supported):
        raise ValueError("source temporal RMS requires a nonzero saved covariance")
    omega = np.sqrt(numerator[supported]/mass[supported])/normal_clock[supported]
    angular = jets["bundle"]["fine_system"].kappa/fine.r
    return {"coordinate_positive_levels": levels.tolist(), "hermitian_gap": hermitian_gap,
            "normal_frequency_positive_shell_range": [float(np.min(proper_shells)), float(np.max(proper_shells))],
            "positive_shell_over_pole_range": [float(np.min(proper_shells)/pole_mass), float(np.max(proper_shells)/pole_mass)],
            "worldline_positive_shell_over_pole": (levels/(normal_clock[index]*pole_mass)).tolist(),
            "angular_kappa_over_r": summary(angular),
            "angular_over_pole_range": [float(np.min(angular)/pole_mass), float(np.max(angular)/pole_mass)],
            "source_temporal_RMS_normal_clock": summary(omega),
            "source_temporal_RMS_over_pole_range": [float(np.min(omega)/pole_mass), float(np.max(omega)/pole_mass)],
            "source_temporal_RMS_worldline_over_pole": None if mass[index] == 0 else
                float(np.sqrt(numerator[index]/mass[index])/(normal_clock[index]*pole_mass)),
            "zero_density_nodes_excluded": int(np.count_nonzero(~supported)),
            "evolved_source_occupies_six_instantaneous_modes_claimed": False,
            "temporal_RMS_is_particle_frequency_or_adiabatic_certificate": False,
            "multiplicity_in_frequency": False}


def measure_station(directory, identifier, ordinal, cw_over_a, pole_mass):
    record, arrays = episode.load_checkpoint(directory, identifier, ordinal)
    state = episode.state_from_arrays(arrays, record["momentum_representation"])
    before = episode.state_sha256(state)
    pair = backend.make_fft_pair(episode.pair_from_arrays(arrays, record))
    jets = tidal.analytic_accelerations(pair, state, "coupled")
    profiles = tidal.profiles_on_grid(pair, jets)
    fine = jets["fine"]
    location = pair.clock_locations[1]
    index = int(np.argmin(abs(pair.grid.xi_q-location)))
    columns = np.vstack((state.phi0, state.phi1))
    gram = columns.conj().T@columns
    roots = np.sqrt(pair.weights)
    car = np.linalg.eigvalsh(roots[:, None]*gram*roots[None, :])
    time_value = float(record["coordinate_time"])
    target = ORDINAL_STATIONS[ordinal]
    result = {"case_id": identifier, "ordinal": ordinal, "requested_station": target,
              "actual_time": time_value, "requested_station_reached": abs(time_value-target) < 1e-9,
              "status": record.get("status"), "stop": record.get("stop"), "nf": pair.grid.nf,
              "control_mode": record.get("control_mode"), "worldline_x": float(pair.grid.xi_q[index]),
              "normal_clocks": np.asarray(arrays["normal_clocks"]).tolist(), "r": summary(fine.r), "Q": summary(fine.Q),
              "chi_auxiliary_dimensionless": summary(fine.chi),
              "normal_lapse_rQ": summary(fine.r*fine.Q), "chart_positive": bool(np.min(fine.r)>0 and np.min(fine.Q)>0),
              "gram_max": float(np.max(abs(gram-np.eye(6)))), "CAR_nonzero_eigenvalues": car.tolist(),
              "occupation_weights": pair.weights.tolist(), "source_trace": float(np.sum(car)),
              "curvature": physical_curvature(profiles, cw_over_a, index),
              "frequencies": frequency_scales(pair, state, jets, pole_mass, index),
              "state_sha256": before, "state_unchanged": episode.state_sha256(state) == before}
    if not result["state_unchanged"]:
        raise ValueError("postprocessing changed a saved state")
    return result


def units_and_scope():
    return {"units": {"t_x_and_carrier_period": "dimensionless chart coordinates in the locked reference units",
                      "r": "physical length", "Q_L_chi_R_h": "dimensionless", "N_q_rQ": "physical length",
                      "normal_time_d_tau": "r*Q*dt; physical length (c=hbar=1)",
                      "action_A": "length^-2", "C_W": "dimensionless", "abs_C_W_over_A": "length^2",
                      "R4_and_tetrad_Riemann_Ricci": "length^-2", "Ricci2_K_Weyl2": "length^-4",
                      "normal_frequency_and_pole_mass": "length^-1", "coordinate_H_levels": "dimensionless",
                      "normal_stress_T": "length^-4"},
            "formula_R4": "(R_h-2)/r^2 - 6*(r_tt-r_xx)/(r^3*Q^2)",
            "formula_Weyl2": "(R_h-2)^2/(3*r^4); never chi^2 substitution",
            "formula_domain_indicator": "abs(C_W)/A times a physical length^-2 curvature scale",
            "formula_pole_mass": "sqrt(A/(2*abs(C_W))); a flat extra-mode comparison scale, not a matched cutoff certificate",
            "formula_normal_shell_frequency": "instantaneous coordinate H eigenvalue/(rQ) at the stated normal observer",
            "formula_source_temporal_RMS": "sqrt(sum_j c_j |U_f dot(phi_j)|^2 / sum_j c_j |U_f phi_j|^2)/(rQ)",
            "formula_local_uniform_dispersion": "sqrt((k/(rQ))^2+(kappa/r)^2), k=2*pi*(m+1/2)/period",
            "state_scope": "saved nf128 half-step coupled Gaussian columns, fixed action/flux/kappa/weights; no reselected source",
            "shell_scope": "six lowest positive instantaneous representative-generator modes; not asserted occupied after evolution",
            "validity_requirements": ["all relevant |c_i curvature|/A small", "|c_i T|/A^2 small",
                                      "external scales below matched local expansion scales", "finite-regulator/state/boundary Jacobian treatment"],
            "unmeasured_requirements": ["full stress and higher-derivative scale budget", "matched nonlocal remainder",
                                         "quantum Jacobian", "propagated/continuum error enclosure"],
            "interpretation": "indicators of loss of separation for the local curvature EFT assumptions; no invalidity conclusion for the full nonlocal theory",
            "raw_R_h_used_as_physical_curvature": False, "new_trajectory": False, "geometry_evolved": False,
            "threshold_gate": None, "EFT_validity_certified": False, "full_nonlocal_theory_invalid": False}


def build_report(dynamic_directory=DEFAULT_DYNAMIC, balanced_directory=DEFAULT_BALANCED):
    start = time.process_time()
    directories = {"dynamic": Path(dynamic_directory).expanduser().resolve(),
                   "balanced": Path(balanced_directory).expanduser().resolve()}
    inputs, saved_producers, files, preparations = {}, [], [], {}
    for family, directory in directories.items():
        manifest = directory/"manifest.json"
        if not manifest.is_file():
            raise FileNotFoundError("saved family manifest unavailable: "+str(manifest))
        files.append(manifest)
        for case in CASES[family]:
            identifier = "nf128_"+case+"_dt0.0005"
            for ordinal in range(4):
                record_path = directory/(identifier+f"-{ordinal:06d}.json")
                record = json.loads(record_path.read_text())
                files += [record_path, record_path.with_suffix(".npz")]
                binding = {"commit": record["producing_commit"], "hashes": record["source_hashes"]}
                if binding not in saved_producers:
                    saved_producers.append(binding)
        metadata = json.loads(manifest.read_text())
        initial = Path(metadata["initial_record"])
        if not initial.is_absolute():
            initial = LAB/initial
        original = json.loads(initial.read_text())
        files += [initial, initial.parent/original["payload"]]
        preparations[family] = {"record": str(initial), "source_case": original.get("source_case"),
            "harmonic": original.get("harmonic"), "lift": original.get("lift"),
            "curvature_criterion": original.get("curvature_criterion"),
            "source_eta": original.get("source_eta"), "multiplicity": original.get("multiplicity"),
            "case_declarations": {case: {key: details[key] for key in
                ("rotation_angle", "stationary_spectral_control", "source_current_deleted", "holding_claim")
                if key in details} for case, details in original["cases"].items()}}
    inputs = {str(path): digest_file(path) for path in sorted(set(files))}
    sources = source_hashes()
    coefficients = coupling.locked_coefficients()
    factor = abs(coefficients["C_W"])/coefficients["A"]
    mass = float(np.sqrt(1/(2*factor)))
    rows = []
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        for family, directory in directories.items():
            for case in CASES[family]:
                for ordinal in range(4):
                    row = measure_station(directory, "nf128_"+case+"_dt0.0005", ordinal, factor, mass)
                    row.update(family=family, preparation_case=case)
                    rows.append(row)
    if any(digest_file(path) != digest for path, digest in inputs.items()) or source_hashes() != sources:
        raise ValueError("input or measurement source changed during finite postprocessing")
    return {"schema": SCHEMA, "scope": units_and_scope(), "locked_coefficients": coefficients,
            "abs_CW_over_A": factor, "extra_flat_pole_mass": mass,
            "directories": {key: str(path) for key, path in directories.items()},
            "initial_preparations": preparations,
            "rows": rows, "input_hashes": inputs, "measurement_source_hashes": sources,
            "saved_producer_bindings": saved_producers, "measurement_producer_commit": None,
            "input_and_source_hashes_unchanged": True, "record_written": False,
            "runtime": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
                        "numerical_threads": 1}, "cpu_seconds": time.process_time()-start}


def write_report(report, output, producer_commit):
    destination = Path(output).expanduser().resolve()
    if destination.exists():
        raise FileExistsError("refusing to overwrite exact EFT diagnostic output")
    authenticate_sources(producer_commit, report["measurement_source_hashes"])
    for binding in report["saved_producer_bindings"]:
        authenticate_sources(binding["commit"], binding["hashes"])
    if source_hashes() != report["measurement_source_hashes"] or any(
            digest_file(path) != expected for path, expected in report["input_hashes"].items()):
        raise ValueError("bound input or source changed before immutable record creation")
    stored = dict(report, measurement_producer_commit=producer_commit, record_written=True)
    stored["record_sha256"] = hashlib.sha256(canonical(stored)).hexdigest()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("x") as stream:
        stream.write(json.dumps(stored, indent=2, allow_nan=False)+"\n")
    destination.chmod(0o444)
    return stored


def check_report(path):
    path = Path(path).expanduser().resolve()
    if path.stat().st_mode & 0o222:
        raise ValueError("EFT diagnostic record is not immutable")
    record = json.loads(path.read_text())
    digest = record.pop("record_sha256")
    if record.get("schema") != SCHEMA or hashlib.sha256(canonical(record)).hexdigest() != digest:
        raise ValueError("EFT diagnostic record digest/schema changed")
    if record["scope"] != units_and_scope() or record["input_and_source_hashes_unchanged"] is not True:
        raise ValueError("EFT diagnostic scope or input-preservation contract changed")
    authenticate_sources(record["measurement_producer_commit"], record["measurement_source_hashes"])
    for binding in record["saved_producer_bindings"]:
        authenticate_sources(binding["commit"], binding["hashes"])
    if source_hashes() != record["measurement_source_hashes"]:
        raise ValueError("current measurement source differs from authenticated producer; refusing substitute replay")
    for name, expected in record["input_hashes"].items():
        if digest_file(name) != expected:
            raise ValueError("EFT diagnostic bound input changed: "+name)
    repeated = build_report(record["directories"]["dynamic"], record["directories"]["balanced"])
    if repeated["input_hashes"] != record["input_hashes"] or len(repeated["rows"]) != len(record["rows"]):
        raise ValueError("EFT diagnostic replay input inventory changed")
    def compare(first, second, name="rows"):
        if isinstance(first, dict):
            if first.keys() != second.keys(): raise ValueError("replay keys changed: "+name)
            for key in first: compare(first[key], second[key], name+"."+key)
        elif isinstance(first, list):
            if len(first) != len(second): raise ValueError("replay length changed: "+name)
            for i, (a, b) in enumerate(zip(first, second)): compare(a, b, name+"."+str(i))
        elif isinstance(first, float):
            if not np.isclose(first, second, rtol=2e-9, atol=2e-12):
                raise ValueError("finite numerical replay changed: "+name)
        elif first != second:
            raise ValueError("replay changed: "+name)
    for key in ("locked_coefficients", "abs_CW_over_A", "extra_flat_pole_mass", "initial_preparations", "rows"):
        compare(repeated[key], record[key], key)
    return {"schema": SCHEMA, "ok": True, "rows_checked": len(record["rows"]),
            "immutable_closure_authenticated": True, "numeric_replay": True,
            "evolved": False, "bytes_written": 0, "cpu_seconds": repeated["cpu_seconds"]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", metavar="PATH")
    mode.add_argument("--check", metavar="PATH")
    parser.add_argument("--producer-commit")
    parser.add_argument("--dynamic-dir", type=Path, default=DEFAULT_DYNAMIC)
    parser.add_argument("--balanced-dir", type=Path, default=DEFAULT_BALANCED)
    args = parser.parse_args(argv)
    try:
        if args.check:
            if args.producer_commit:
                parser.error("--producer-commit belongs only to explicit --write")
            result = check_report(args.check)
        else:
            if bool(args.write) != bool(args.producer_commit):
                parser.error("explicit --write requires --producer-commit; default is read-only")
            report = build_report(args.dynamic_dir, args.balanced_dir)
            result = write_report(report, args.write, args.producer_commit) if args.write else report
        print(json.dumps(result, indent=2, allow_nan=False))
        return 0
    except (ValueError, RuntimeError, FileNotFoundError, FileExistsError, KeyError) as error:
        print(type(error).__name__+": "+str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
