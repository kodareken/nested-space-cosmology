"""Finite connected geometric-force vertices of an authenticated live slice.

No state or geometry is evolved. The representative Q vertex has kappa but
no angular multiplicity. Gaussian fermion fields obey Wick's theorem;
their quadratic force observables can have a nonzero third cumulant.
The total angular stress-noise mapping is deliberately not supplied.
"""
from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
import scipy
from threadpoolctl import threadpool_limits

from . import nsc_discovery_backend as backend
from . import nsc_discovery_coupled_memory as memory
from . import nsc_discovery_episode as episode
from . import nsc_influence as influence
from . import nsc_nested_parent_child as model
from . import nsc_spherical_coupling as coupling


LAB = Path(__file__).resolve().parents[2]
REPO = LAB.parent
CLI = LAB / "scripts/derive_nsc_discovery_geometric_vertices.py"
DEFAULT_OUTPUT = LAB / "results/development/nsc-discovery-geometric-vertices-v1.json"
SCHEMA = "NSC-DISCOVERY-GEOMETRIC-VERTICES-v1"
CPU_LIMIT = 10.
COUNTING_H = (.04, .02)
HELD_OUT_COUNTING_FIELD = .03


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_digest(values):
    values = np.ascontiguousarray(values)
    hasher = hashlib.sha256()
    hasher.update(str(values.dtype).encode())
    hasher.update(str(values.shape).encode())
    hasher.update(values.tobytes())
    return hasher.hexdigest()


def computational_paths():
    """Local Python import closure, plus all data read by this calculation.

    This is the calculation's byte closure, not a replay of the historical
    campaigns that originally produced its authenticated inputs.
    """
    package = LAB / "src/recursive_horizons"
    pending = [Path(__file__).resolve(), CLI.resolve()]
    visited = set()
    while pending:
        path = pending.pop()
        if path in visited:
            continue
        if not path.is_file():
            raise FileNotFoundError("computational source is unavailable: " + str(path))
        visited.add(path)
        for node in ast.walk(ast.parse(path.read_text())):
            names = []
            if isinstance(node, ast.ImportFrom):
                if node.level == 1:
                    names = [node.module] if node.module else [alias.name for alias in node.names]
                elif node.module and node.module.startswith("recursive_horizons"):
                    stem = node.module.removeprefix("recursive_horizons").lstrip(".")
                    names = [stem] if stem else [alias.name for alias in node.names]
            elif isinstance(node, ast.Import):
                names = [alias.name.removeprefix("recursive_horizons.") for alias in node.names
                         if alias.name.startswith("recursive_horizons.")]
            for name in names:
                target = package / (name.replace(".", "/") + ".py")
                if target.is_file() and target.resolve() not in visited:
                    pending.append(target.resolve())
    visited.update((episode.V1_JSON.resolve(), episode.V1_NPZ.resolve(),
                    episode.BASIS_JSON.resolve(), episode.BASIS_NPZ.resolve(),
                    coupling.LOCKED_RECORD.resolve(), (LAB / "pyproject.toml").resolve()))
    return tuple(sorted(visited))


def source_bindings():
    return {path.relative_to(REPO).as_posix(): digest(path) for path in computational_paths()}


def authenticate_git(commit, bindings):
    """Authenticate every declared byte at one explicit frozen Git commit."""
    if not isinstance(commit, str) or len(commit) != 40 or any(c not in "0123456789abcdef" for c in commit):
        raise ValueError("an explicit full lowercase frozen Git commit is required")
    requests = [commit + ":" + path for path in bindings]
    completed = subprocess.run(["git", "cat-file", "--batch"], cwd=REPO,
                               input=("\n".join(requests) + "\n").encode(),
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    data, offset = completed.stdout, 0
    objects = {}
    for path, request in zip(bindings, requests):
        end = data.index(b"\n", offset)
        header = data[offset:end].decode()
        offset = end + 1
        if header.endswith(" missing"):
            raise ValueError("calculation dependency is absent at the frozen commit: " + path)
        object_id, kind, size = header.split()
        if kind != "blob":
            raise ValueError("calculation dependency is not a Git blob: " + path)
        blob = data[offset:offset + int(size)]
        offset += int(size) + 1
        if hashlib.sha256(blob).hexdigest() != bindings[path]:
            raise ValueError("current calculation bytes differ from the frozen commit: " + path)
        objects[path] = object_id
    return {"commit": commit, "all_declared_bytes_authenticated": True, "blob_ids": objects}


def geometric_vertex(pair, canonical_index):
    """G_j=∂H_G/∂a_Qj; full ambient Galerkin operator, not a six-mode truncation."""
    j = int(canonical_index)
    if not 0 <= j < pair.grid.ng:
        raise ValueError("canonical Q mode is outside the retained geometry band")
    profile = pair.grid.A_g @ pair.geometry_map[:, j]
    interpolation = pair.grid.U_f @ np.eye(pair.grid.nf, dtype=complex)
    mass = pair.grid.fine.kappa * (interpolation.conj().T @ (profile[:, None] * interpolation))
    zero = np.zeros_like(mass)
    vertex = np.block([[zero, mass], [mass, zero]])
    inside = (pair.grid.xi_q >= pair.child_interval[0]) & (pair.grid.xi_q <= pair.child_interval[1])
    power = float(np.sum(abs(profile) ** 2))
    group = "child_detail" if j in pair.geometry_child_indices else "parent_coarse" if j in pair.geometry_parent_indices else "remaining_geometry"
    return vertex, {"canonical_Q_index": j, "geometry_group": group,
        "direction_convention": "unit increment of a_Qj; nodal direction W[:,j]",
        "mixed_coupling": "delta a_Qj * phi_dagger G_j phi",
        "profile_sha256": array_digest(profile), "profile_min": float(np.min(profile)),
        "profile_max": float(np.max(profile)), "profile_squared_power": power,
        "sampled_profile_power_outside_child_fraction": float(np.sum(abs(profile[~inside]) ** 2) / power),
        "parent_interval": list(pair.parent_interval), "child_interval": list(pair.child_interval),
        "exact_compact_support_claimed": False, "representative_vertex_has_M": False,
        "kappa": int(pair.grid.fine.kappa), "operator_dimension": vertex.shape[0],
        "vertex_sha256": array_digest(vertex)}


def bilinear_cumulants(covariance, vertex):
    """First three cumulants of f=c†G c, in one finite Gaussian channel."""
    covariance = np.asarray(covariance, dtype=complex)
    vertex = np.asarray(vertex, dtype=complex)
    if covariance.shape != vertex.shape or covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
        raise ValueError("covariance and vertex must be matching finite matrices")
    if not np.allclose(vertex, vertex.conj().T, rtol=0, atol=1e-11):
        raise ValueError("physical geometric vertex must be Hermitian")
    g2, cg = vertex @ vertex, covariance @ vertex
    terms = (np.trace(covariance @ (g2 @ vertex)),
             -3 * np.trace(cg @ covariance @ g2), 2 * np.trace(cg @ cg @ cg))
    return {"kappa1": float(np.trace(cg).real),
        "kappa2": float(np.trace(covariance @ g2 - cg @ cg).real),
        "kappa3": float(sum(terms).real), "kappa3_imaginary_residual": float(abs(sum(terms).imag)),
        "kappa3_trace_terms": [{"real": float(t.real), "imag": float(t.imag)} for t in terms],
        "elementary_fermion_state_is_Gaussian": True, "bilinear_force_is_not_assumed_Gaussian": True}


def counting_logdet(columns, weights, vertex, amplitude, *, eigensystem=None):
    """log det(I-C+C exp(-i lambda G)), using the exact thin Sylvester determinant.

    exp(-i lambda G) still acts in the FULL ambient band. The six-dimensional
    determinant comes from the source rank and does not discard G leakage.
    """
    columns, weights = np.asarray(columns, dtype=complex), np.asarray(weights, dtype=float)
    if columns.ndim != 2 or weights.shape != (columns.shape[1],) or vertex.shape != (columns.shape[0], columns.shape[0]):
        raise ValueError("counting field needs the full source columns and matching weights")
    values, basis = np.linalg.eigh(vertex) if eigensystem is None else eigensystem
    image = (basis * np.exp(-1j * float(amplitude) * values)[None, :]) @ (basis.conj().T @ columns)
    kernel = np.eye(len(weights)) + weights[:, None] * (columns.conj().T @ (image - columns))
    sign, logarithm = np.linalg.slogdet(kernel)
    if sign == 0:
        raise ValueError("counting characteristic is singular on the declared branch")
    return complex(logarithm, np.angle(sign))


def independent_counting_check(columns, weights, vertex, cumulants):
    eigensystem = np.linalg.eigh(vertex)
    phase = lambda value: counting_logdet(columns, weights, vertex, value, eigensystem=eigensystem).imag
    estimates = []
    for h in COUNTING_H:
        third = (phase(2 * h) - 2 * phase(h) + 2 * phase(-h) - phase(-2 * h)) / (2 * h ** 3)
        estimates.append({"h": h, "third_derivative": float(third),
            "gap_from_trace_kappa3": abs(float(third) - cumulants["kappa3"])})
    lam = HELD_OUT_COUNTING_FIELD
    observed_odd = (phase(lam) - phase(-lam)) / 2 + lam * cumulants["kappa1"]
    predicted_cubic = lam ** 3 * cumulants["kappa3"] / 6
    return {"characteristic": "det(I-C+C exp(-i lambda G))",
        "independent_method": "finite counting-field logdet derivative; full-band exponential on thin source",
        "derivative_steps": estimates, "held_out_lambda": lam,
        "held_out_not_a_derivative_step": lam not in COUNTING_H,
        "held_out_odd_phase_after_linear_mean": float(observed_odd),
        "predicted_cubic_odd_phase": float(predicted_cubic),
        "held_out_gap": float(abs(observed_odd - predicted_cubic)),
        "held_out_gap_over_cubic_effect": None if predicted_cubic == 0 else float(abs(observed_odd - predicted_cubic) / abs(predicted_cubic)),
        "odd_remainder_order": "O(lambda^5); no dynamical force term is installed",
        "principal_branch_operator_norm_times_largest_lambda": float(np.max(abs(eigensystem[0])) * (2 * max(COUNTING_H))),
        "elementary_nongaussian_interaction_added": False}


def occupied_empty_cut(covariance, hamiltonian, child_vertex, parent_vertex):
    """Nonzero instantaneous retarded cross-slope, with an occupied/empty cut.

    A=i[H,G_child] is the local Heisenberg derivative, not a new interaction.
    No time propagation, stationary spectral assumption, or pole fit occurs.
    """
    a = 1j * (hamiltonian @ child_vertex - child_vertex @ hamiltonian)
    greater = influence.connected(covariance, a, parent_vertex)
    lesser = influence.connected(covariance, parent_vertex, a)
    commutator = np.trace(covariance @ (a @ parent_vertex - parent_vertex @ a))
    encode = lambda z: {"real": float(z.real), "imag": float(z.imag)}
    return {"A": "i[H_G,G_child5]", "B": "G_parent0",
        "greater": encode(greater), "lesser": encode(lesser), "commutator": encode(commutator),
        "factorization_identity": "Tr(C A (I-C) B)-Tr(C B (I-C) A)=Tr(C[A,B])",
        "cut_identity_gap": float(abs(greater - lesser - commutator)),
        "instantaneous_retarded_cross_slope": float((-1j * commutator).real),
        "retarded_convention": "H_int=+J B; chi_AB=-i theta <[A(t),B(s)]>",
        "stationarity_assumed": False, "physical_pole_claimed": False,
        "new_geometric_force_added": False}


def calculate():
    """Default: bounded saved-slice calculation only; no evidence writes or evolution."""
    started = time.process_time()
    before = source_bindings()
    with threadpool_limits(limits=1), backend.fft_thread_limit(1):
        loaded = memory.load_matched_handoff(128, backend_name="fft")
        pair, state = loaded["pair"], loaded["state"]
        columns = np.vstack((state.phi0, state.phi1))
        weights = np.asarray(pair.weights)
        covariance = (columns * weights) @ columns.conj().T
        gram = columns.conj().T @ columns
        roots = np.sqrt(weights)
        eigenvalues = np.linalg.eigvalsh(roots[:, None] * gram * roots[None, :])
        if eigenvalues[0] < -1e-10 or eigenvalues[-1] > 1 + 1e-9:
            raise ValueError("actual saved covariance is not CAR admissible")
        selected = (int(pair.geometry_parent_indices[0]), int(pair.geometry_child_indices[0]), int(pair.geometry_child_indices[1]))
        if selected != (0, 5, 6):
            raise ValueError("authenticated geometry basis no longer names the declared modes 0,5,6")
        force = model.source_geometry_forces(pair, state)["Q_energy_gradient"]
        vertices, reports = {}, {}
        multiplicity = int(pair.grid.fine.multiplicity)
        for j in selected:
            vertex, profile = geometric_vertex(pair, j)
            vertices[j] = vertex
            cumulants = bilinear_cumulants(covariance, vertex)
            independent = independent_counting_check(columns, weights, vertex, cumulants)
            reports[str(j)] = {"physical_operator": profile, "representative_channel": cumulants,
                "counting_field_check": independent, "mean_force_gradient": float(force[j]),
                "mean_M_trace_CG": multiplicity * cumulants["kappa1"],
                "mean_force_mapping_gap": float(abs(force[j] - multiplicity * cumulants["kappa1"]))}
            if time.process_time() - started > CPU_LIMIT:
                raise RuntimeError("saved-slice control exceeds its admitted 10 CPU seconds")
        hamiltonian = model.hamiltonian(pair, state)
        cut = occupied_empty_cut(covariance, hamiltonian, vertices[5], vertices[0])
        stationarity_gap = float(np.linalg.norm(hamiltonian @ covariance - covariance @ hamiltonian))
    after = source_bindings()
    if after != before:
        raise RuntimeError("computational source or authenticated input changed during the calculation")
    elapsed = time.process_time() - started
    if elapsed > CPU_LIMIT:
        raise RuntimeError("saved-slice control exceeds its admitted 10 CPU seconds")
    return {"schema": SCHEMA, "status": "MEASURED_FINITE_GEOMETRIC_INFLUENCE_VERTEX",
        "case": loaded["case"], "nf": 128, "coordinate_time": .3, "source_pins": loaded["pins"],
        "state_sha256": episode.state_sha256(state), "hamiltonian_sha256": array_digest(hamiltonian),
        "source_covariance_sha256": array_digest(covariance), "basis_W_sha256": array_digest(pair.geometry_map),
        "weights": weights.tolist(), "covariance_nonzero_eigenvalues": eigenvalues.tolist(),
        "CAR_admissible": True, "column_gram_gap": float(np.max(abs(gram - np.eye(6)))),
        "source_rank": 6, "operator_dimension": 2 * pair.grid.nf, "full_ambient_band_retained": True,
        "covariance_stationarity_gap": stationarity_gap, "vertices": reports, "occupied_empty_cut": cut,
        "multiplicity_scope": {"M": multiplicity, "mean": "M Tr(C G), verified against the existing source gradient",
            "derived_cumulants": "one representative channel of the existing Gaussian determinant",
            "scaled_observable_assumption": "if f_scaled=M f on that same channel, kappa3(f_scaled)=M^3 kappa3(f)",
            "independent_angular_copy_assumption": "if rho_total is a product of M identical channel states, kappa3(sum f_s)=M kappa3(f)",
            "angular_product_state_derived_here": False, "total_stress_noise_mapping": None},
        "evolution_ran": False, "initial_state_called": False, "new_dynamics": False,
        "physical_frequency_pole_claimed": False, "cosmohedron_equivalence_claimed": False,
        "scope": "finite geometric bilinear influence cumulants and an instantaneous causal cut on the actual saved slice",
        "source_bindings_before": before, "source_bindings_after": after,
        "runtime": {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__},
        "cpu_seconds": elapsed, "cpu_limit_seconds": CPU_LIMIT}


def write_record(output, commit):
    """Exclusive calculation record, after all bytes match the explicit frozen commit."""
    path = Path(output).resolve()
    if path.exists():
        raise FileExistsError("geometric-vertex record exists; use --check or a named successor")
    record = calculate()
    record["frozen_git"] = authenticate_git(commit, record["source_bindings_before"])
    if source_bindings() != record["source_bindings_before"]:
        raise RuntimeError("source changed before the exclusive record write")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return record


def check_record(output):
    """Read-only byte authentication and bounded calculation replay; no evolution."""
    record = json.loads(Path(output).read_text())
    if record.get("schema") != SCHEMA:
        raise ValueError("geometric-vertex record has the wrong schema")
    bindings = source_bindings()
    if bindings != record["source_bindings_before"] or bindings != record["source_bindings_after"]:
        raise ValueError("geometric-vertex input or producer bytes changed")
    frozen = authenticate_git(record["frozen_git"]["commit"], bindings)
    if frozen != record["frozen_git"]:
        raise ValueError("frozen Git object closure changed")
    replay = calculate()

    def compare(saved, fresh, key="record"):
        if isinstance(fresh, dict):
            if set(saved) != set(fresh):
                raise ValueError("saved calculation fields changed: " + key)
            for name in fresh:
                if name != "cpu_seconds":
                    compare(saved[name], fresh[name], key + "." + name)
        elif isinstance(fresh, list):
            if len(saved) != len(fresh):
                raise ValueError("saved calculation shape changed: " + key)
            for i, item in enumerate(fresh):
                compare(saved[i], item, key + "." + str(i))
        elif isinstance(fresh, float):
            if not np.isfinite(saved) or not np.isclose(saved, fresh, rtol=1e-8, atol=1e-12):
                raise ValueError("saved geometric-vertex number changed: " + key)
        elif saved != fresh:
            raise ValueError("saved geometric-vertex scope changed: " + key)

    compare({key: value for key, value in record.items() if key != "frozen_git"}, replay)
    return record
