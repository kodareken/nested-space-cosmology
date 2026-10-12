"""Fixed nested observers and conditional response of one live finite action.

The child is the middle T0 pair; its parent contains all six T0 modes.
Child, parent detail, and ambient complement form one orthonormal frame.
The Hamiltonian is the representative Dirac generator, WITHOUT the source
copy multiplicity. Multiplicity enters energy/variational traces only.

This module does not produce geometry. ``retained_response`` consumes the
same evolved H(g(t)) and actual segment-start source as its full reference.
It preserves the entire ambient source, including every covariance cross
term. No time-indexed dense exterior propagator is stored.
"""
from __future__ import annotations

from dataclasses import dataclass
import time

import numpy as np

from .nsc_coupled_local_response import evolve_full_columns
from .nsc_evolving_reduction import evolve_retained_region, fixed_observer_frame


SECTORS = ("child", "detail", "ambient")


def _matrix(values, label):
    matrix = np.asarray(values, dtype=np.complex128)
    if matrix.ndim != 2 or not np.isfinite(matrix).all():
        raise ValueError(label + " must be a finite matrix")
    return matrix


def _hermitian(values, dimension):
    matrix = _matrix(values, "Hamiltonian")
    if matrix.shape != (dimension, dimension):
        raise ValueError("Hamiltonian dimension differs from the fixed nested frame")
    scale = max(1., float(np.linalg.norm(matrix)))
    if np.linalg.norm(matrix - matrix.conj().T) > 1e-10 * scale:
        raise ValueError("representative Hamiltonian must be Hermitian")
    return (matrix + matrix.conj().T) / 2


def _source(columns, weights, dimension):
    columns = _matrix(columns, "source columns")
    weights = np.asarray(weights, dtype=float)
    if columns.shape[0] != dimension or columns.shape[1] == 0 or weights.shape != (columns.shape[1],):
        raise ValueError("source columns and weights have incompatible dimensions")
    if not np.isfinite(weights).all() or np.any(weights < 0) or np.any(weights > 1):
        raise ValueError("finite Gaussian weights must lie in [0,1]")
    gram = columns.conj().T @ columns
    spectrum = np.linalg.eigvalsh(np.sqrt(weights)[:, None] * gram * np.sqrt(weights)[None, :])
    if spectrum[0] < -1e-10 or spectrum[-1] > 1 + 1e-9:
        raise ValueError("weighted source covariance is not admissible")
    return columns, weights


@dataclass(frozen=True)
class NestedFrame:
    """Child+detail+ambient completion; ``parent`` keeps original T0 order."""
    child: np.ndarray
    detail: np.ndarray
    ambient: np.ndarray
    parent: np.ndarray
    frame: np.ndarray
    child_indices: tuple[int, ...]
    detail_indices: tuple[int, ...]

    @property
    def dimension(self):
        return self.frame.shape[0]

    @property
    def slices(self):
        child = self.child.shape[1]
        parent = self.parent.shape[1]
        return {"child": slice(0, child), "detail": slice(child, parent),
                "ambient": slice(parent, self.dimension)}


def nested_frame(reference_columns, child_indices=(2, 3)):
    """Preserve both original pair phases and all six parent columns exactly.

    Spatial inclusion is a separately declared physical-domain statement;
    finite Fourier observer tails are not replaced by an exact support claim.
    """
    parent = np.array(_matrix(reference_columns, "original T0 observer columns"), copy=True)
    if parent.shape[1] != 6 or parent.shape[0] <= 6:
        raise ValueError("the nested preparation requires six modes and a nonempty ambient complement")
    indices = tuple(child_indices)
    if len(indices) != 2 or any(isinstance(i, bool) or not isinstance(i, (int, np.integer)) for i in indices):
        raise ValueError("child_indices must contain two distinct integer mode indices")
    if len(set(indices)) != 2 or min(indices) < 0 or max(indices) >= 6:
        raise ValueError("child indices must name two distinct original parent modes")
    detail_indices = tuple(i for i in range(6) if i not in indices)
    child = parent[:, indices].copy()
    detail = parent[:, detail_indices].copy()
    reordered = np.concatenate((child, detail), axis=1)
    completed = fixed_observer_frame(reordered)
    ambient = completed["exterior_basis"]
    return NestedFrame(child, detail, ambient, parent, completed["frame"], indices, detail_indices)


def nested_blocks(hamiltonian, hierarchy):
    """Geometry-derived blocks; no guessed link or inserted scale/copy factor."""
    matrix = _hermitian(hamiltonian, hierarchy.dimension)
    transformed = hierarchy.frame.conj().T @ matrix @ hierarchy.frame
    slices = hierarchy.slices
    return {"matrix": transformed, **{
        left + "_" + right: transformed[slices[left], slices[right]]
        for left in SECTORS for right in SECTORS
    }}


def _trace_parts(operator, columns, weights, hierarchy, multiplicity):
    if not np.isfinite(multiplicity) or multiplicity <= 0:
        raise ValueError("source multiplicity must be positive")
    matrix = _hermitian(operator, hierarchy.dimension)
    columns, weights = _source(columns, weights, hierarchy.dimension)
    transformed = hierarchy.frame.conj().T @ matrix @ hierarchy.frame
    coefficients = hierarchy.frame.conj().T @ columns
    diagonal = np.zeros(3)
    interaction = np.zeros((3, 3))
    currents = np.zeros((3, 3))
    population = np.zeros(3)
    slices = tuple(hierarchy.slices[name] for name in SECTORS)
    for i, si in enumerate(slices):
        zi = coefficients[si]
        population[i] = float(np.sum(weights * np.sum(abs(zi) ** 2, axis=0)))
        value = np.sum(weights * np.sum(zi.conj() * (transformed[si, si] @ zi), axis=0))
        diagonal[i] = multiplicity * value.real
        for j in range(i + 1, 3):
            sj = slices[j]
            zj = coefficients[sj]
            value = np.sum(weights * np.sum(zi.conj() * (transformed[si, sj] @ zj), axis=0))
            interaction[i, j] = 2 * multiplicity * value.real
            currents[i, j], currents[j, i] = 2 * value.imag, -2 * value.imag
    direct = float(multiplicity * np.sum(weights * np.sum(columns.conj() * (matrix @ columns), axis=0)).real)
    assembled = float(np.sum(diagonal) + np.sum(interaction))
    return coefficients, transformed, {"sectors": SECTORS, "population": population,
        "pair_current_into_row": currents, "population_derivative": currents.sum(axis=1),
        "diagonal_energy": diagonal, "interaction_energy_upper": interaction,
        "source_energy": direct, "assembled_energy": assembled, "energy_closure_gap": abs(direct - assembled),
        "parent_internal_energy": float(diagonal[0] + diagonal[1] + interaction[0, 1]),
        "ambient_internal_energy": float(diagonal[2]),
        "parent_ambient_interaction_energy": float(interaction[0, 2] + interaction[1, 2]),
        "multiplicity": float(multiplicity)}


def block_accounting(hamiltonian, columns, weights, hierarchy, *, multiplicity=4.):
    """Full source trace and reciprocal probability currents at one live state.

    Pair currents are representative canonical probabilities, without M.
    Interaction energies are counted once, in the upper triangle. Parent
    energy already includes child: adding their complete energies is invalid.
    Only small retained covariance blocks are returned, never a dense history.
    """
    coefficients, transformed, report = _trace_parts(hamiltonian, columns, weights, hierarchy, multiplicity)
    weights = np.asarray(weights, dtype=float)
    rate = -1j * transformed @ coefficients
    direct_rates = np.array([2 * np.sum(weights * np.sum(coefficients[sl].conj() * rate[sl], axis=0)).real
                             for sl in hierarchy.slices.values()])
    report["current_closure_gap"] = float(np.max(abs(direct_rates - report["population_derivative"])))
    c, d, e = (coefficients[hierarchy.slices[name]] for name in SECTORS)
    report["child_covariance"] = (c * weights) @ c.conj().T
    report["detail_covariance"] = (d * weights) @ d.conj().T
    report["child_detail_covariance"] = (c * weights) @ d.conj().T
    report["child_ambient_covariance_norm"] = float(np.linalg.norm((c * weights) @ e.conj().T))
    report["detail_ambient_covariance_norm"] = float(np.linalg.norm((d * weights) @ e.conj().T))
    report["full_covariance_included"] = True
    return report


def force_accounting(hamiltonian_derivatives, columns, weights, hierarchy, *, multiplicity=4.):
    """M Tr(C ∂H/∂q), retaining link derivatives and all source cross terms.

    Inputs are derivatives with respect to the caller's canonical geometry
    coordinates. Returned traces are energy derivatives, not forces divided
    by a grid spacing or a substitute for the geometric Hamiltonian terms.
    """
    reports = {}
    for name, derivative in hamiltonian_derivatives.items():
        parts = _trace_parts(derivative, columns, weights, hierarchy, multiplicity)[2]
        # Probability currents of ∂H are not the live generator's currents.
        reports[name] = {key: parts[key] for key in ("sectors", "diagonal_energy", "interaction_energy_upper",
            "source_energy", "assembled_energy", "energy_closure_gap", "multiplicity")}
        reports[name].update(kind="variational_energy_derivative", value=parts["source_energy"], full_covariance_included=True)
    return reports


def nested_schur(hamiltonian, hierarchy, z=1 + 1j):
    """Eliminate ambient then parent detail, with the full resolvent as control."""
    z = complex(z)
    if not np.isfinite(z) or z.imag <= 0:
        raise ValueError("Schur control requires z in the upper half-plane")
    h = nested_blocks(hamiltonian, hierarchy)["matrix"]
    dimension, child, parent = hierarchy.dimension, hierarchy.child.shape[1], hierarchy.parent.shape[1]
    inverse = z * np.eye(dimension) - h
    pp, pe, ep, ee = inverse[:parent, :parent], inverse[:parent, parent:], inverse[parent:, :parent], inverse[parent:, parent:]
    parent_inverse = pp - pe @ np.linalg.solve(ee, ep)
    cc, cd, dc, dd = parent_inverse[:child, :child], parent_inverse[:child, child:], parent_inverse[child:, :child], parent_inverse[child:, child:]
    child_inverse = cc - cd @ np.linalg.solve(dd, dc)
    joined_inverse = inverse[:child, :child] - inverse[:child, child:] @ np.linalg.solve(inverse[child:, child:], inverse[child:, :child])
    rhs = np.eye(dimension)[:, :child]
    direct = np.linalg.solve(inverse, rhs)[:child]
    nested = np.linalg.solve(child_inverse, np.eye(child))
    joined = np.linalg.solve(joined_inverse, np.eye(child))
    ambient_omitted = np.linalg.solve(pp, np.eye(parent)[:, :child])[:child]
    return {"z": z, "parent_inverse_response": parent_inverse, "child_inverse_response": child_inverse,
            "direct_child_response": direct, "nested_child_response": nested, "joined_child_response": joined,
            "nested_direct_gap": float(np.max(abs(nested - direct))),
            "joined_direct_gap": float(np.max(abs(joined - direct))),
            "nested_joined_inverse_gap": float(np.max(abs(child_inverse - joined_inverse))),
            "ambient_omission_response_gap": float(np.max(abs(ambient_omitted - direct))),
            "all_ambient_links_retained": True}


def _covariance(amplitudes, weights):
    return np.einsum("tak,k,tbk->tab", amplitudes, weights, amplitudes.conj())


def _diagonal(covariance):
    return np.diagonal(covariance, axis1=1, axis2=2).real


def _movement(values):
    return float(np.max(abs(values)))


def _comparison(approximate, reference, effect):
    error = _movement(np.asarray(approximate) - np.asarray(reference))
    return {"movement": error, "effect": float(effect),
            "fraction": None if effect == 0 else error / effect,
            "within_one_percent": bool(effect > 0 and error <= .01 * effect)}


def response_allocation_plan(hierarchy, output_times, source_rank=6):
    """Array sizes only; CPU must be measured before a production campaign."""
    if output_times < 2 or source_rank < 1:
        raise ValueError("need at least two output times and a nonempty source")
    d, nt, r = hierarchy.dimension, int(output_times), int(source_rank)
    return {"dimension": d, "source_rank": r, "child_split_columns": 3 * r,
            "child_history_bytes": nt * (d - 2) * (3 * r) * 16,
            "parent_history_bytes": nt * (d - 6) * r * 16,
            "full_split_column_frames_bytes": nt * d * (3 * r) * 16,
            "one_dense_frame_bytes": d * d * 16,
            "dense_propagator_history_bytes": 0,
            "cpu_forecast_available": False}


def _allocation(result):
    allocation = result["allocation"]
    if allocation["time_indexed_exterior_propagator_bytes"] != 0 or result["exterior_propagator"] is not None:
        raise RuntimeError("nested response stored a dense exterior propagator history")
    return {**allocation, "history_shape": list(result["history"].shape),
            "history_bytes": int(result["history"].nbytes)}


def retained_response(hamiltonian, times, initial_columns, weights, hierarchy, *, substeps=2,
                      reference_substeps=4, autonomous_columns=None, cpu_limit_s=60.):
    """Same-evolved-source child/parent reduction and independent controls.

    ``hamiltonian`` must be the actual evolved metric schedule. The initial
    columns are the actual transported source at this segment boundary, while
    ``hierarchy`` stays the T0 frame. No observer reset or geometry evolution
    occurs here. Three-way source splitting distinguishes parent detail from
    the ambient complement. Every omission is a conditional negative control,
    not a separately coupled geometry or an additive share of the effect.
    """
    columns, weights = _source(initial_columns, weights, hierarchy.dimension)
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or len(times) < 2 or not np.isfinite(times).all() or np.any(np.diff(times) <= 0):
        raise ValueError("response needs a finite increasing time grid")
    if not np.allclose(np.diff(times), np.diff(times)[0], rtol=0, atol=1e-12):
        raise ValueError("streamed reduction requires a uniform time grid")
    if any(isinstance(x, bool) or not isinstance(x, (int, np.integer)) or x < 1 for x in (substeps, reference_substeps)):
        raise ValueError("midpoint substeps must be positive integers")
    if reference_substeps <= substeps:
        raise ValueError("independent reference must use more midpoint substeps")
    if not callable(hamiltonian) or not np.isfinite(cpu_limit_s) or cpu_limit_s <= 0:
        raise ValueError("need a callable actual H(g(t)) and a positive CPU ceiling")
    if autonomous_columns is not None:
        autonomous_columns = np.asarray(autonomous_columns, dtype=np.complex128)
        if autonomous_columns.shape != (len(times),) + columns.shape or not np.isfinite(autonomous_columns).all():
            raise ValueError("autonomous source frames differ from the response domain")
        if not np.array_equal(autonomous_columns[0], columns):
            raise ValueError("full/reduced initial source differs from the actual autonomous segment")
    started = time.process_time()

    def checked(mark):
        if time.process_time() - started > cpu_limit_s:
            raise RuntimeError("nested response exceeded its admitted CPU ceiling")
        return _hermitian(hamiltonian(mark), hierarchy.dimension)

    components = [basis @ (basis.conj().T @ columns) for basis in (hierarchy.child, hierarchy.detail, hierarchy.ambient)]
    split = np.concatenate(components, axis=1)
    if np.max(abs(sum(components) - columns)) > 1e-9:
        raise ValueError("nested source split did not reconstruct the entire preparation")
    rank = columns.shape[1]
    streamed = evolve_retained_region(checked, hierarchy.child, times, split, weights=np.tile(weights, 3),
                                      backend="streamed", propagator_substeps=int(substeps))
    parent = evolve_retained_region(checked, hierarchy.parent, times, columns, weights=weights,
                                   backend="streamed", propagator_substeps=int(substeps))
    # A separate full-column equation supplies each control's own reference.
    full = evolve_full_columns(checked, times, split, substeps=int(substeps))
    refined = evolve_full_columns(checked, times, split, substeps=int(reference_substeps))
    projected = np.einsum("ij,tjk->tik", hierarchy.child.conj().T, full)
    reduced = streamed["amplitudes"]

    def assemble(projected_pieces):
        c, d, e = (projected_pieces[:, :, i * rank:(i + 1) * rank] for i in range(3))
        total = _covariance(c + d + e, weights)
        marginal = _covariance(c, weights) + _covariance(d + e, weights)
        child_detail_cross = _covariance(c + d, weights) - _covariance(c, weights) - _covariance(d, weights)
        control_covariances = {"outside_drive_off": _covariance(c, weights),
            "outside_cross_removed": marginal, "parent_detail_drive_off": _covariance(c + e, weights),
            "child_detail_cross_removed": total - child_detail_cross}
        return {"covariance": total, "occupation": _diagonal(total),
                **{name: _diagonal(value) for name, value in control_covariances.items()},
                **{name + "_covariance": value for name, value in control_covariances.items()}}

    child = assemble(reduced)
    child_full = assemble(projected)
    projected_refined = np.einsum("ij,tjk->tik", hierarchy.child.conj().T, refined)
    child_refined = assemble(projected_refined)
    child_effect = _movement(child_full["occupation"] - child_full["occupation"][0])
    child["full_covariance"] = child_full["covariance"]
    child["full_covariance_refined"] = child_refined["covariance"]
    child["occupation_full"] = child_full["occupation"]
    child["occupation_full_refined"] = child_refined["occupation"]
    child["coherence"] = child["covariance"][:, 0, 1]
    child["coherence_full"] = child_full["covariance"][:, 0, 1]
    coherence_effect = _movement(child["coherence_full"] - child["coherence_full"][0])
    child["coherence_error"] = _comparison(child["coherence"], child["coherence_full"], coherence_effect)
    child["reduction_error"] = _comparison(child["occupation"], child_full["occupation"], child_effect)
    child["full_midpoint_indicator"] = _comparison(child_full["occupation"], child_refined["occupation"], child_effect)
    child["covariance_error_max"] = float(np.max(np.linalg.norm(child["covariance"] - child_full["covariance"], axis=(1, 2))))
    child["controls"] = {}
    for name in ("outside_drive_off", "outside_cross_removed", "parent_detail_drive_off", "child_detail_cross_removed"):
        effect = _movement(child_full[name] - child_full["occupation"])
        reduced_effect = child[name] - child["occupation"]
        full_effect = child_full[name] - child_full["occupation"]
        refined_effect = child_refined[name] - child_refined["occupation"]
        child["controls"][name] = {"occupation_movement": effect,
            "own_reference_error": _comparison(reduced_effect, full_effect, effect),
            "control_curve_reference_error": _comparison(child[name], child_full[name], effect),
            "reference_midpoint_indicator": _comparison(full_effect, refined_effect, effect)}
        coherence_effect = _movement(child_full[name + "_covariance"][:, 0, 1] - child["coherence_full"])
        child["controls"][name]["coherence_movement"] = coherence_effect
        child["controls"][name]["coherence_reference_error"] = _comparison(
            child[name + "_covariance"][:, 0, 1] - child["coherence"],
            child_full[name + "_covariance"][:, 0, 1] - child["coherence_full"], coherence_effect)
        child["controls"][name]["admissible_full_covariance_claimed"] = name != "child_detail_cross_removed"
        child[name + "_full"] = child_full[name]
        child[name + "_full_covariance"] = child_full[name + "_covariance"]
        child[name + "_full_refined"] = child_refined[name]
    reconstructed = sum(full[:, :, i * rank:(i + 1) * rank] for i in range(3))
    reconstructed_refined = sum(refined[:, :, i * rank:(i + 1) * rank] for i in range(3))
    parent_projection = np.einsum("ij,tjk->tik", hierarchy.parent.conj().T, reconstructed)
    parent_projection_refined = np.einsum("ij,tjk->tik", hierarchy.parent.conj().T, reconstructed_refined)
    parent_covariance = _covariance(parent["amplitudes"], weights)
    parent_full_covariance = _covariance(parent_projection, weights)
    parent_refined_covariance = _covariance(parent_projection_refined, weights)
    parent_effect = _movement(_diagonal(parent_full_covariance) - _diagonal(parent_full_covariance)[0])
    parent_report = {"covariance": parent_covariance, "full_covariance": parent_full_covariance,
        "full_covariance_refined": parent_refined_covariance,
        "occupation": _diagonal(parent_covariance), "occupation_full": _diagonal(parent_full_covariance),
        "reduction_error": _comparison(_diagonal(parent_covariance), _diagonal(parent_full_covariance), parent_effect),
        "full_midpoint_indicator": _comparison(_diagonal(parent_full_covariance), _diagonal(parent_refined_covariance), parent_effect)}
    parent_indices = np.asarray(hierarchy.child_indices)
    nested_gap = _movement(parent_full_covariance[:, parent_indices[:, None], parent_indices] - child_full["covariance"])
    initial = block_accounting(checked(float(times[0])), columns, weights, hierarchy)
    if not np.array_equal(streamed["local_basis"], hierarchy.child) or not np.array_equal(parent["local_basis"], hierarchy.parent):
        raise RuntimeError("a nested T0 observer was reset or rephased")
    report = {"times": times.copy(), "child": child, "parent": parent_report,
        "nested_full_covariance_gap": nested_gap,
        "initial_child_detail_cross_norm": float(np.linalg.norm(initial["child_detail_covariance"])),
        "initial_child_ambient_cross_norm": initial["child_ambient_covariance_norm"],
        "initial_detail_column_norm": float(np.linalg.norm(components[1])),
        "initial_ambient_column_norm": float(np.linalg.norm(components[2])),
        "allocation": {"child": _allocation(streamed), "parent": _allocation(parent)},
        "allocation_plan": response_allocation_plan(hierarchy, len(times), rank),
        "observer_preserved": True, "all_source_components_retained": True,
        "geometry_rerun": False, "reduction_feeds_geometry": False,
        "scope": "conditional reduction of the same live coupled trajectory; parent detail distinct from ambient",
        "controls_additive": False, "state_error_certified": False}
    if autonomous_columns is not None:
        actual = np.einsum("ij,tjk->tik", hierarchy.child.conj().T, autonomous_columns)
        actual_covariance = _covariance(actual, weights)
        child["autonomous_covariance"] = actual_covariance
        child["conditional_versus_autonomous"] = _comparison(child_full["occupation"], _diagonal(actual_covariance), child_effect)
    report["cpu_seconds"] = time.process_time() - started
    if report["cpu_seconds"] > cpu_limit_s:
        raise RuntimeError("nested response exceeded its admitted CPU ceiling")
    return report
